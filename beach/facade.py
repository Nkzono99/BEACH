"""Session facade for one Fortran output directory."""

from __future__ import annotations

import functools
import inspect
import textwrap
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Callable, Iterator, Mapping

import numpy as np

from .analysis._selection import _build_mesh_selection, _mesh_element_indices
from .analysis.context import RunContext
from .analysis.coulomb import calc_coulomb
from .analysis.field_kernel import (
    _require_total_field_config,
    calc_object_forces_kernel,
    field_kernel_from_result,
    field_kernel_options_from_result,
)
from .analysis.field_lines import compute_electric_field_points, trace_field_lines
from .analysis.mobility import analyze_coulomb_mobility
from .analysis.object_interaction import ObjectInteractionSnapshot
from .analysis.potential import (
    _uses_saved_mesh_potential,
    compute_potential_mesh,
    compute_potential_points,
    compute_potential_slices,
)
from .analysis.scene import BeachScene
from .io.results import load_run
from .io.types import RunResult
from .native.kernel import FieldKernel
from .plot.animation import animate_history_mesh
from .plot.boxplot import plot_mesh_source_boxplot
from .plot.field_lines import plot_field_lines_3d
from .plot.matrix import plot_coulomb_force_matrix
from .plot.mesh import plot_charge_mesh, plot_charges, plot_potential_mesh
from .plot.slices import plot_potential_slices

_NO_KERNEL = object()
_PLAN_ARGUMENTS = ("periodic2", "theta", "leaf_max", "order")


def _forward(
    func: Callable[..., Any],
    *,
    kernel_step: object = _NO_KERNEL,
    needs_kernel: Callable[["Beach", Mapping[str, Any]], bool] | None = None,
) -> Callable[..., Any]:
    """Expose ``func`` as a ``Beach`` method that keeps its signature and defaults.

    The method passes the session as ``result`` and the session's
    ``library_path`` when none is given. With ``kernel_step`` it also passes a
    session kernel for the call's plan arguments (see :meth:`Beach._lend_kernel`),
    unless the caller
    gives ``kernel``, ``config_path``, or another ``library_path`` (those calls
    build their own kernel), or ``needs_kernel`` says this call evaluates none.
    """

    signature = inspect.signature(func)
    parameters = list(signature.parameters.values())
    method_signature = signature.replace(
        parameters=[
            inspect.Parameter("self", inspect.Parameter.POSITIONAL_OR_KEYWORD),
            *parameters[1:],
        ]
    )
    accepted = method_signature.parameters

    @functools.wraps(func)
    def method(self: "Beach", *args: Any, **kwargs: Any) -> Any:
        bound = method_signature.bind(self, *args, **kwargs)
        bound.apply_defaults()
        arguments = dict(bound.arguments)
        del arguments["self"]
        library_path = arguments.get("library_path")
        own_plan = (
            arguments.get("kernel") is not None
            or arguments.get("config_path") is not None
            or (library_path is not None and Path(library_path) != self.library_path)
        )
        if "library_path" in accepted and library_path is None:
            arguments["library_path"] = self.library_path
        if (
            kernel_step is not _NO_KERNEL
            and not own_plan
            and (needs_kernel is None or needs_kernel(self, arguments))
        ):
            _require_total_field_config(self.context, operation=func.__name__)
            with self._lend_kernel(
                step=arguments.get("step", kernel_step),
                **{name: arguments[name] for name in _PLAN_ARGUMENTS if name in accepted},
            ) as kernel:
                arguments["kernel"] = kernel
                return func(self, **arguments)
        return func(self, **arguments)

    method.__signature__ = method_signature  # type: ignore[attr-defined]
    target = f"{func.__module__}.{func.__qualname__}"
    note = f"Session method for :func:`{target}` with this run as ``result``."
    if "library_path" in accepted:
        note += " ``library_path`` defaults to the session's."
    if kernel_step is not _NO_KERNEL:
        note += (
            " A kernel-evaluating call receives the session's cached kernel for its "
            "plan arguments as ``kernel``, holding the call's step charges; passing "
            "``kernel``, ``config_path``, or another ``library_path`` builds a kernel "
            "for the call instead."
        )
    method.__doc__ = f"{inspect.cleandoc(func.__doc__ or '')}\n\nNotes\n-----\n{textwrap.fill(note, 88)}"
    return method


def _potential_needs_kernel(run: "Beach", arguments: Mapping[str, Any]) -> bool:
    return not _uses_saved_mesh_potential(
        run,
        source=arguments["source"],
        periodic2=arguments["periodic2"],
        reference_point=arguments["reference_point"],
        kernel=None,
    )


class Beach:
    """Session over one Fortran output directory.

    A ``Beach`` owns the loaded result, its :class:`RunContext` (the one place that
    resolves ``beach.toml``), and the native field kernels built for it. Methods that
    need a kernel reuse :meth:`kernel`, so one plan serves several figures. Close the
    session, or use it as a context manager, to release those kernels.

    Parameters
    ----------
    output_dir : str or pathlib.Path, default "outputs/latest"
        Directory containing Fortran output files.
    config_path : str, pathlib.Path, or None, default None
        ``beach.toml`` used by config-aware analyses. ``None`` searches
        ``output_dir``, its parent, and its grandparent.
    library_path : str, pathlib.Path, or None, default None
        Native field-kernel shared library. ``None`` uses ``BEACH_FIELD_KERNEL_LIB``
        or the default search.
    """

    def __init__(
        self,
        output_dir: str | Path = "outputs/latest",
        *,
        config_path: str | Path | None = None,
        library_path: str | Path | None = None,
    ) -> None:
        self.output_dir = Path(output_dir)
        self.library_path = None if library_path is None else Path(library_path)
        self._config_path = None if config_path is None else Path(config_path)
        self._context: RunContext | None = None
        self._kernels: dict[tuple[object, ...], FieldKernel] = {}

    @property
    def context(self) -> RunContext:
        """Return the run context (loads the result on first access).

        Returns
        -------
        RunContext
            Loaded result and its lazily resolved ``beach.toml``.
        """

        if self._context is None:
            self._context = RunContext(
                result=load_run(self.output_dir),
                requested_config_path=self._config_path,
            )
        return self._context

    @property
    def result(self) -> RunResult:
        """Return the loaded run result.

        Returns
        -------
        RunResult
            Parsed result for ``output_dir``.
        """

        return self.context.result

    def reload(self) -> RunResult:
        """Close the cached kernels and reload output files from disk.

        Returns
        -------
        RunResult
            Refreshed parsed result.
        """

        self.close()
        self._context = None
        return self.result

    def kernel(
        self,
        *,
        step: int | None = -1,
        periodic2: Mapping[str, object] | None = None,
        theta: float | None = None,
        leaf_max: int | None = None,
        order: int | None = None,
    ) -> FieldKernel:
        """Return the session's native field kernel for one plan (heavy on first use).

        The kernel is built with
        :func:`beach.analysis.field_kernel.field_kernel_from_result` and cached by
        its charge values, resolved field-kernel options, and the session's
        ``library_path``. Equivalent defaults and explicit options, and steps
        with the same charges, share one kernel. The session owns the kernel: do not close
        it, and restore its charges after :meth:`FieldKernel.update_charges`.
        :meth:`close` and :meth:`reload` release every cached kernel.

        Parameters
        ----------
        step : int or None, default -1
            History step whose charges load the kernel. ``-1`` selects the latest
            history snapshot, ``None`` the final charges.
        periodic2 : mapping or None, default None
            Two-axis periodic setting. ``None`` resolves it from the run.
        theta, leaf_max, order : float, int, int, or None, default None
            FMM plan overrides. ``None`` uses the run's receipt or config.

        Returns
        -------
        FieldKernel
            Open kernel owned by this session.
        """

        key = self._kernel_key(
            step, periodic2=periodic2, theta=theta, leaf_max=leaf_max, order=order
        )
        kernel = self._kernels.get(key)
        if kernel is None or kernel.closed:
            kernel = field_kernel_from_result(
                self,
                step=step,
                periodic2=periodic2,
                theta=theta,
                leaf_max=leaf_max,
                order=order,
                library_path=self.library_path,
            )
            self._kernels[key] = kernel
        return kernel

    def _kernel_key(self, step: int | None, **plan: Any) -> tuple[object, ...]:
        charges = np.asarray(self.context.charges_at(step), dtype=np.float64)
        options = field_kernel_options_from_result(self, **plan)
        return (charges.tobytes(), options, self.library_path)

    @contextmanager
    def _lend_kernel(self, *, step: int | None, **plan: Any) -> Iterator[FieldKernel]:
        """Lend a session kernel holding the ``step`` charges for one method call.

        A session method needs the plan, not a handle of its own: when no kernel
        holds these charges but one with the same plan is open, the call borrows
        it with the ``step`` charges and gets its own charges back afterwards. So
        methods called for many steps share one plan instead of caching one per
        step. :meth:`kernel` handles keep their charges outside the call.
        """

        key = self._kernel_key(step, **plan)
        kernel = self._kernels.get(key)
        if kernel is None or kernel.closed:
            kernel = next(
                (
                    candidate
                    for candidate_key, candidate in self._kernels.items()
                    if candidate_key[1:] == key[1:] and not candidate.closed
                ),
                None,
            )
            if kernel is not None:
                charges = kernel.source_charges
                kernel.update_charges(self.context.charges_at(step))
                try:
                    yield kernel
                finally:
                    kernel.update_charges(charges)
                return
            kernel = self.kernel(step=step, **plan)
        yield kernel

    def close(self) -> None:
        """Close every cached native kernel. The loaded result stays available."""

        kernels = list(self._kernels.values())
        self._kernels.clear()
        for kernel in kernels:
            kernel.close()

    def __enter__(self) -> "Beach":
        return self

    def __exit__(self, *_exc_info: object) -> None:
        self.close()

    @property
    def mesh_ids(self) -> tuple[int, ...]:
        """Return available mesh ids.

        Returns
        -------
        tuple of int
            Sorted unique mesh identifiers.
        """

        ids = np.unique(self.context.mesh_ids_or_default())
        return tuple(int(v) for v in ids)

    def get_mesh(self, *mesh_ids: int, step: int | None = -1):
        """Return mesh selections for requested ids.

        Parameters
        ----------
        *mesh_ids : int
            One or more mesh identifiers.
        step : int or None, default -1
            History step used to select charges. ``-1`` selects latest history,
            ``None`` uses final charges.

        Returns
        -------
        MeshSelection or tuple[MeshSelection, ...]
            Single selection when one id is given, tuple otherwise.

        Raises
        ------
        ValueError
            If no mesh ids are provided or a mesh id is invalid.
        """

        if len(mesh_ids) == 0:
            raise ValueError("at least one mesh id must be provided.")
        selections = tuple(
            _build_mesh_selection(self.result, (int(mesh_id),), step=step)
            for mesh_id in mesh_ids
        )
        if len(selections) == 1:
            return selections[0]
        return selections

    def get_mesh_charge(self, *mesh_ids: int, step: int | None = -1):
        """Return per-element charges for requested mesh ids.

        Parameters
        ----------
        *mesh_ids : int
            One or more mesh identifiers.
        step : int or None, default -1
            History step used to select charges.

        Returns
        -------
        numpy.ndarray or tuple[numpy.ndarray, ...]
            Charge array for each requested mesh.
        """

        if len(mesh_ids) == 0:
            raise ValueError("at least one mesh id must be provided.")
        indices = [
            _mesh_element_indices(self.context, (int(mesh_id),))[1]
            for mesh_id in mesh_ids
        ]
        charges = np.asarray(self.context.charges_at(step))
        # Integer indexing copies, so each returned array is independent.
        selected = tuple(charges[index] for index in indices)
        if len(selected) == 1:
            return selected[0]
        return selected

    calc_coulomb = _forward(calc_coulomb, kernel_step=-1)
    analyze_coulomb_mobility = _forward(analyze_coulomb_mobility)
    object_interaction_snapshot = _forward(ObjectInteractionSnapshot.from_result)
    calc_object_forces_kernel = _forward(calc_object_forces_kernel, kernel_step=-1)
    scene = _forward(BeachScene.from_result)
    plot_mesh = _forward(plot_charge_mesh)
    plot_bar = _forward(plot_charges)
    compute_potential = _forward(
        compute_potential_mesh,
        kernel_step=None,
        needs_kernel=_potential_needs_kernel,
    )
    compute_potential_points = _forward(compute_potential_points, kernel_step=None)
    compute_potential_slices = _forward(compute_potential_slices, kernel_step=None)
    plot_potential = _forward(
        plot_potential_mesh,
        kernel_step=None,
        needs_kernel=_potential_needs_kernel,
    )
    plot_potential_slices = _forward(plot_potential_slices, kernel_step=None)
    plot_mesh_source_boxplot = _forward(plot_mesh_source_boxplot)
    plot_coulomb_force_matrix = _forward(
        plot_coulomb_force_matrix,
        kernel_step=-1,
        needs_kernel=lambda _run, arguments: arguments["workers"] == 1,
    )
    animate_mesh = _forward(
        animate_history_mesh,
        kernel_step=-1,
        needs_kernel=lambda _run, arguments: str(arguments["quantity"]).lower()
        == "potential",
    )
    compute_electric_field = _forward(compute_electric_field_points, kernel_step=-1)
    trace_field_lines = _forward(trace_field_lines, kernel_step=-1)
    plot_field_lines = _forward(plot_field_lines_3d, kernel_step=-1)


for _name, _member in vars(Beach).items():
    if hasattr(_member, "__wrapped__"):
        _member.__name__ = _name
        _member.__qualname__ = f"Beach.{_name}"
del _name, _member

"""Triangle-panel potential evaluation through the native BEACH field kernel."""

from __future__ import annotations

from contextlib import AbstractContextManager, nullcontext
from pathlib import Path
from typing import Iterable, Literal, Mapping

import numpy as np

from beach.io._numeric import _points
from beach.io.types import FortranRunResult, PotentialSlice2D
from beach.native.kernel import FieldKernel, FieldKernelOptions
from beach.native.periodic2 import Periodic2Config, coerce_periodic2

from .context import RunContext
from .field_kernel import (
    _options_from_result,
    _require_complete_total_kernel,
    _require_run_kernel,
    _require_total_field_config,
    _require_total_field_reconstruction,
    _restoring_charges,
)
from .mesh import _triangle_centers
from .periodic import auto_periodic2_from_result


def compute_potential_mesh(
    result: FortranRunResult | object,
    *,
    periodic2: Mapping[str, object] | None = None,
    reference_point: Iterable[float] | str | None = None,
    config_path: str | Path | None = None,
    library_path: str | Path | None = None,
    source: Literal["auto", "saved", "kernel"] = "auto",
    kernel: FieldKernel | None = None,
) -> np.ndarray:
    """Evaluate triangle-P0 potential at every panel centroid (heavy unless saved values are used).

    Parameters
    ----------
    result : FortranRunResult, RunContext, or Beach
        Run to evaluate.
    periodic2 : mapping or None, default None
        Two-axis periodic setting. ``None`` resolves it from the run.
    reference_point : iterable of float, {"species1_injection_center"}, or None, default None
        Point whose potential is subtracted. ``None`` returns absolute potential.
    config_path, library_path : str, pathlib.Path, or None, default None
        Explicit ``beach.toml`` and native kernel library.
    source : {"auto", "saved", "kernel"}, default "auto"
        ``"saved"`` returns the simulator's ``mesh_potential.csv`` (light) and
        rejects ``periodic2``, ``reference_point``, and ``kernel``. ``"kernel"``
        always evaluates the native kernel (heavy). ``"auto"`` uses the saved values
        when they exist and none of ``periodic2``, ``reference_point``, or
        ``kernel`` is given, and the kernel otherwise.
    kernel : FieldKernel or None, default None
        Open kernel built from this run, for example ``Beach.kernel(step=None)``.
        It is evaluated with its own plan and charges and is not closed;
        ``periodic2`` and ``library_path`` then do not apply. ``None`` builds a
        kernel for this call and closes it.

    Returns
    -------
    numpy.ndarray
        Potential in volts with shape ``(mesh_nelem,)``.
    """

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    context.require_triangle_source_model()
    if _uses_saved_mesh_potential(
        context,
        source=source,
        periodic2=periodic2,
        reference_point=reference_point,
        kernel=kernel,
    ):
        return np.asarray(resolved.mesh_potential_v, dtype=float).copy()

    reference_xyz = _resolve_reference_point(
        resolved,
        reference_point,
        context=context,
    )
    points = _triangle_centers(context.require_triangles())
    with _potential_kernel(
        context,
        periodic2=periodic2,
        library_path=library_path,
        kernel=kernel,
    ) as evaluator:
        potential_v = evaluator.eval_phi(points)
        if reference_xyz is not None:
            potential_v = potential_v - float(
                evaluator.eval_phi(reference_xyz.reshape(1, 3))[0]
            )
    return potential_v


def _uses_saved_mesh_potential(
    result: FortranRunResult | object,
    *,
    source: str,
    periodic2: Mapping[str, object] | None,
    reference_point: Iterable[float] | str | None,
    kernel: FieldKernel | None,
) -> bool:
    """Decide between saved ``mesh_potential.csv`` and kernel evaluation."""

    if source == "kernel":
        return False
    if source not in {"auto", "saved"}:
        raise ValueError('source must be "auto", "saved", or "kernel".')
    context = RunContext.from_value(result)
    saved = context.result.mesh_potential_v is not None
    reference_xyz = _resolve_reference_point(
        context.result,
        reference_point,
        context=context,
    )
    unmodified = periodic2 is None and reference_xyz is None and kernel is None
    if source == "auto":
        return saved and unmodified
    if not saved:
        raise ValueError('source="saved" requires mesh_potential.csv.')
    if not unmodified:
        raise ValueError(
            'source="saved" cannot be combined with periodic2, reference_point, '
            'or kernel; use source="kernel".'
        )
    return True


def compute_potential_points(
    result: FortranRunResult | object,
    points: np.ndarray,
    *,
    chunk_size: int = 2048,
    periodic2: Mapping[str, object] | None = None,
    reference_point: Iterable[float] | str | None = None,
    config_path: str | Path | None = None,
    library_path: str | Path | None = None,
    kernel: FieldKernel | None = None,
) -> np.ndarray:
    """Evaluate triangle-P0 potential at arbitrary points (heavy: native kernel).

    Parameters
    ----------
    result : FortranRunResult, RunContext, or Beach
        Run to evaluate.
    points : numpy.ndarray
        Sampling points with shape ``(n_points, 3)`` in meters.
    chunk_size : int, default 2048
        Number of points evaluated per kernel call.
    periodic2 : mapping or None, default None
        Two-axis periodic setting. ``None`` resolves it from the run.
    reference_point : iterable of float, {"species1_injection_center"}, or None, default None
        Point whose potential is subtracted. ``None`` returns absolute potential.
    config_path, library_path : str, pathlib.Path, or None, default None
        Explicit ``beach.toml`` and native kernel library.
    kernel : FieldKernel or None, default None
        Open kernel built from this run, for example ``Beach.kernel(step=None)``.
        It is evaluated with its own plan and charges and is not closed;
        ``periodic2`` and ``library_path`` then do not apply. ``None`` builds a
        kernel for this call and closes it.

    Returns
    -------
    numpy.ndarray
        Potential in volts with shape ``(n_points,)``.
    """

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    context.require_triangle_source_model()
    if chunk_size <= 0:
        raise ValueError("chunk_size must be > 0.")

    sample_points = _points(points)
    if sample_points.shape[0] == 0:
        return np.empty(0, dtype=float)

    reference_xyz = _resolve_reference_point(
        resolved,
        reference_point,
        context=context,
    )
    with _potential_kernel(
        context,
        periodic2=periodic2,
        library_path=library_path,
        kernel=kernel,
    ) as evaluator:
        potential_v = _eval_phi_chunks(evaluator, sample_points, chunk_size)
        if reference_xyz is not None:
            potential_v = potential_v - float(
                evaluator.eval_phi(reference_xyz.reshape(1, 3))[0]
            )
    return potential_v


def compute_potential_slices(
    result: FortranRunResult | object,
    *,
    box_min: Iterable[float],
    box_max: Iterable[float],
    grid_n: int = 200,
    xy_z: float | None = None,
    yz_x: float | None = None,
    xz_y: float | None = None,
    chunk_size: int = 2048,
    periodic2: Mapping[str, object] | None = None,
    reference_point: Iterable[float] | str | None = None,
    config_path: str | Path | None = None,
    library_path: str | Path | None = None,
    kernel: FieldKernel | None = None,
) -> dict[str, PotentialSlice2D]:
    """Sample triangle-P0 potential on XY, YZ, and XZ slices (heavy: native kernel).

    The planes cross ``box_min``-``box_max`` at ``xy_z``, ``yz_x``, and ``xz_y``
    (box center when ``None``) on a ``grid_n`` x ``grid_n`` grid. The other
    arguments, including ``kernel``, are those of :func:`compute_potential_points`.
    """

    if grid_n < 2:
        raise ValueError("grid_n must be >= 2.")
    min_corner, max_corner = _coerce_box_bounds(box_min, box_max)
    x = np.linspace(min_corner[0], max_corner[0], grid_n, dtype=float)
    y = np.linspace(min_corner[1], max_corner[1], grid_n, dtype=float)
    z = np.linspace(min_corner[2], max_corner[2], grid_n, dtype=float)
    z_xy = _coerce_slice_coordinate(
        xy_z,
        lower=min_corner[2],
        upper=max_corner[2],
        name="xy_z",
    )
    x_yz = _coerce_slice_coordinate(
        yz_x,
        lower=min_corner[0],
        upper=max_corner[0],
        name="yz_x",
    )
    y_xz = _coerce_slice_coordinate(
        xz_y,
        lower=min_corner[1],
        upper=max_corner[1],
        name="xz_y",
    )

    xx, yy = np.meshgrid(x, y, indexing="xy")
    yy2, zz = np.meshgrid(y, z, indexing="xy")
    xx2, zz2 = np.meshgrid(x, z, indexing="xy")
    points_xy = np.column_stack(
        (xx.reshape(-1), yy.reshape(-1), np.full(xx.size, z_xy))
    )
    points_yz = np.column_stack(
        (np.full(yy2.size, x_yz), yy2.reshape(-1), zz.reshape(-1))
    )
    points_xz = np.column_stack(
        (xx2.reshape(-1), np.full(xx2.size, y_xz), zz2.reshape(-1))
    )
    all_points = np.concatenate((points_xy, points_yz, points_xz), axis=0)
    all_potential = compute_potential_points(
        result,
        all_points,
        chunk_size=chunk_size,
        periodic2=periodic2,
        reference_point=reference_point,
        config_path=config_path,
        library_path=library_path,
        kernel=kernel,
    )
    plane_size = grid_n * grid_n
    potential_xy = all_potential[:plane_size].reshape(grid_n, grid_n)
    potential_yz = all_potential[plane_size : 2 * plane_size].reshape(
        grid_n, grid_n
    )
    potential_xz = all_potential[2 * plane_size :].reshape(grid_n, grid_n)

    return {
        "xy": PotentialSlice2D(
            plane="xy",
            axis_u="x",
            axis_v="y",
            fixed_axis="z",
            fixed_value_m=float(z_xy),
            u_values_m=x,
            v_values_m=y,
            potential_V=potential_xy,
        ),
        "yz": PotentialSlice2D(
            plane="yz",
            axis_u="y",
            axis_v="z",
            fixed_axis="x",
            fixed_value_m=float(x_yz),
            u_values_m=y,
            v_values_m=z,
            potential_V=potential_yz,
        ),
        "xz": PotentialSlice2D(
            plane="xz",
            axis_u="x",
            axis_v="z",
            fixed_axis="y",
            fixed_value_m=float(y_xz),
            u_values_m=x,
            v_values_m=z,
            potential_V=potential_xz,
        ),
    }


def _potential_history(
    charges_history: np.ndarray,
    triangles: np.ndarray,
    *,
    periodic2: Periodic2Config | None = None,
    reference_point: np.ndarray | None = None,
    field_options: FieldKernelOptions | None = None,
    library_path: str | Path | None = None,
) -> np.ndarray:
    """Evaluate panel potential history with one reusable native plan."""

    charge_frames = np.asarray(charges_history, dtype=np.float64)
    panels = np.asarray(triangles, dtype=np.float64)
    if charge_frames.ndim != 2:
        raise ValueError("charges_history must have shape (nelem, nframe).")
    if panels.ndim != 3 or panels.shape[1:] != (3, 3):
        raise ValueError("triangles must have shape (nelem, 3, 3).")
    if charge_frames.shape[0] != panels.shape[0]:
        raise ValueError("charges_history and triangles must have the same nelem.")
    if charge_frames.shape[1] == 0:
        return np.empty_like(charge_frames)

    centers = _triangle_centers(panels)
    if field_options is not None and periodic2 is not None:
        raise ValueError("field_options and periodic2 cannot be used together.")
    options = (
        FieldKernelOptions(periodic2=periodic2)
        if field_options is None
        else field_options
    )
    _require_complete_total_kernel(
        options,
        operation="potential history recomputation",
    )
    with FieldKernel(
        panels,
        charge_frames[:, 0],
        options=options,
        library_path=library_path,
    ) as kernel:
        return _potential_frames(kernel, centers, charge_frames, reference_point)


def _potential_frames(
    kernel: FieldKernel,
    centers: np.ndarray,
    charge_frames: np.ndarray,
    reference_point: np.ndarray | None,
) -> np.ndarray:
    """Evaluate centroid potential for each charge frame on one kernel plan."""

    values = np.empty_like(charge_frames)
    reference = None if reference_point is None else _reference_array(reference_point)
    for frame in range(charge_frames.shape[1]):
        kernel.update_charges(charge_frames[:, frame])
        values[:, frame] = kernel.eval_phi(centers)
        if reference is not None:
            values[:, frame] -= float(kernel.eval_phi(reference.reshape(1, 3))[0])
    return values


def _run_potential_history(
    context: RunContext,
    charges_history: np.ndarray,
    *,
    periodic2: Mapping[str, object] | None,
    reference_point: Iterable[float] | str | None,
    library_path: str | Path | None,
    kernel: FieldKernel | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Recompute one run's panel potential history and its resolved reference.

    A caller's ``kernel`` evaluates the frames with its own plan; its source
    charges are restored afterwards.
    """

    resolved = context.result
    context.require_triangle_source_model()
    if kernel is not None:
        _require_run_kernel(
            context,
            kernel,
            operation="potential history recomputation",
        )
        resolved_reference = _resolve_reference_point(
            resolved,
            reference_point,
            context=context,
        )
        centers = _triangle_centers(context.require_triangles())
        frames = np.asarray(charges_history, dtype=np.float64)
        with _restoring_charges(kernel):
            values = _potential_frames(kernel, centers, frames, resolved_reference)
        return values, resolved_reference
    periodic_cfg = coerce_periodic2(
        periodic2,
        allow_cached_kneq0=True,
    )
    if periodic_cfg is None:
        periodic_cfg = auto_periodic2_from_result(
            resolved,
            allow_cached_kneq0=True,
            context=context,
        )
    resolved_reference = _resolve_reference_point(
        resolved,
        reference_point,
        context=context,
    )
    _require_total_field_config(
        context,
        operation="potential history recomputation",
    )
    field_options = _options_from_result(
        resolved,
        periodic2=periodic_cfg,
        theta=None,
        leaf_max=None,
        order=None,
        config_path=context.requested_config_path,
        context=context,
    )
    _require_total_field_reconstruction(
        context,
        field_options,
        operation="potential history recomputation",
    )
    values = _potential_history(
        charges_history,
        context.require_triangles(),
        reference_point=resolved_reference,
        field_options=field_options,
        library_path=library_path,
    )
    return values, resolved_reference


def _potential_kernel(
    context: RunContext,
    *,
    periodic2: Mapping[str, object] | None,
    library_path: str | Path | None,
    kernel: FieldKernel | None,
) -> AbstractContextManager[FieldKernel]:
    """Borrow a caller's kernel or build one that the ``with`` block closes."""

    if kernel is None:
        return _panel_kernel(context, periodic2=periodic2, library_path=library_path)
    _require_run_kernel(context, kernel, operation="potential recomputation")
    return nullcontext(kernel)


def _panel_kernel(
    context: RunContext,
    *,
    periodic2: Mapping[str, object] | None,
    library_path: str | Path | None,
) -> FieldKernel:
    resolved = context.result
    context.require_triangle_source_model()
    triangles = context.require_triangles()
    if periodic2 is not None:
        coerce_periodic2(periodic2, allow_cached_kneq0=True)
    _require_total_field_config(
        context,
        operation="potential recomputation",
    )
    options = _options_from_result(
        resolved,
        periodic2=periodic2,
        theta=None,
        leaf_max=None,
        order=None,
        config_path=context.requested_config_path,
        context=context,
    )
    _require_total_field_reconstruction(
        context,
        options,
        operation="potential recomputation",
    )
    return FieldKernel(
        triangles,
        resolved.charges,
        options=options,
        library_path=library_path,
    )


def _eval_phi_chunks(
    kernel: FieldKernel,
    points: np.ndarray,
    chunk_size: int,
) -> np.ndarray:
    values = np.empty(points.shape[0], dtype=float)
    for start in range(0, points.shape[0], chunk_size):
        stop = min(start + chunk_size, points.shape[0])
        values[start:stop] = kernel.eval_phi(points[start:stop])
    return values


def _coerce_box_bounds(
    box_min: Iterable[float],
    box_max: Iterable[float],
) -> tuple[np.ndarray, np.ndarray]:
    min_corner = np.asarray(list(box_min), dtype=float)
    max_corner = np.asarray(list(box_max), dtype=float)
    if min_corner.shape != (3,) or max_corner.shape != (3,):
        raise ValueError("box_min and box_max must contain exactly 3 values.")
    if not np.all(np.isfinite(min_corner)) or not np.all(np.isfinite(max_corner)):
        raise ValueError("box_min and box_max must contain finite values.")
    if np.any(max_corner <= min_corner):
        raise ValueError("box_max must be greater than box_min on all axes.")
    return min_corner, max_corner


def _coerce_slice_coordinate(
    value: float | None,
    *,
    lower: float,
    upper: float,
    name: str,
) -> float:
    if value is None:
        return 0.5 * (lower + upper)
    coordinate = float(value)
    if not np.isfinite(coordinate) or coordinate < lower or coordinate > upper:
        raise ValueError(
            f"{name} must be within [{lower:.6e}, {upper:.6e}] "
            f"but got {coordinate:.6e}."
        )
    return coordinate


def _resolve_reference_point(
    resolved: FortranRunResult,
    reference_point: Iterable[float] | str | None,
    *,
    context: RunContext | None = None,
) -> np.ndarray | None:
    if reference_point is None:
        return None
    if isinstance(reference_point, str):
        key = reference_point.strip().lower()
        if key in {"species1_injection_center", "species1", "default"}:
            return _species1_injection_center_from_result(
                resolved,
                context=context,
            )
        raise ValueError(
            'reference_point string must be "species1_injection_center" '
            "or a 3D coordinate."
        )
    return _reference_array(reference_point)


def _reference_array(value: Iterable[float] | np.ndarray) -> np.ndarray:
    point = np.asarray(list(value), dtype=float)
    if point.shape != (3,):
        raise ValueError("reference_point must contain exactly 3 values.")
    if not np.all(np.isfinite(point)):
        raise ValueError("reference_point must contain finite values.")
    return point


def _species1_injection_center_from_result(
    resolved: FortranRunResult,
    *,
    context: RunContext | None = None,
) -> np.ndarray | None:
    config = (context or RunContext.from_value(resolved)).config
    if config is None:
        return None
    particles = config.get("particles")
    if not isinstance(particles, Mapping):
        return None
    species = particles.get("species")
    if not isinstance(species, list) or not species:
        return None
    first_species = species[0]
    if not isinstance(first_species, Mapping):
        return None
    if (
        "inject_face" not in first_species
        or "pos_low" not in first_species
        or "pos_high" not in first_species
    ):
        return None
    pos_low = np.asarray(
        _coerce_vec3(
            first_species["pos_low"],
            name="particles.species[0].pos_low",
        )
    )
    pos_high = np.asarray(
        _coerce_vec3(
            first_species["pos_high"],
            name="particles.species[0].pos_high",
        )
    )
    return 0.5 * (pos_low + pos_high)


def _wrap_periodic2_points(
    points: np.ndarray,
    *,
    axes: tuple[int, int],
    lengths: tuple[float, float],
    origins: tuple[float, float],
) -> np.ndarray:
    """Wrap points into a periodic cell without changing non-periodic axes."""

    wrapped = np.asarray(points, dtype=float).copy()
    wrapped[:, axes[0]] = origins[0] + np.mod(
        wrapped[:, axes[0]] - origins[0],
        lengths[0],
    )
    wrapped[:, axes[1]] = origins[1] + np.mod(
        wrapped[:, axes[1]] - origins[1],
        lengths[1],
    )
    return wrapped


def _coerce_vec3(value: object, *, name: str) -> tuple[float, float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError(f"{name} must contain exactly 3 values.")
    try:
        result = (float(value[0]), float(value[1]), float(value[2]))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain numeric values.") from exc
    if not all(np.isfinite(component) for component in result):
        raise ValueError(f"{name} must contain finite values.")
    return result

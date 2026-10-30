"""Native field kernels and their options resolved from one BEACH run."""

from __future__ import annotations

import warnings
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np

from beach.io._numeric import _vec3
from beach.io.types import RunResult
from beach.native.kernel import (
    FieldKernel,
    FieldKernelOptions,
    _coerce_resolved_field_solver,
)
from beach.native.periodic2 import coerce_periodic2
from beach.tooling._cache import default_periodic_cache_dir

from .context import RunContext
from .mesh import _triangle_centers
from .panel_quadrature import panel_target_quadrature
from .periodic import _sim_with_separated_boundaries, periodic2_from_sim


@dataclass(frozen=True)
class KernelObjectForceRecord:
    """Object-level net force/torque computed by the Fortran field kernel."""

    mesh_id: int
    step: int | None
    total_charge_C: float
    center_m: np.ndarray
    force_N: np.ndarray
    torque_Nm: np.ndarray


def field_kernel_from_result(
    result: RunResult | object,
    *,
    step: int | None = -1,
    periodic2: Mapping[str, object] | None = None,
    theta: float | None = None,
    leaf_max: int | None = None,
    order: int | None = None,
    config_path: str | Path | None = None,
    library_path: str | Path | None = None,
) -> FieldKernel:
    """Build a native kernel from one BEACH output directory (heavy).

    The kernel holds the run's panels with the ``step`` charges (``-1`` latest
    history snapshot, ``None`` final ``charges.csv``). The caller owns it and closes
    it; :meth:`beach.Beach.kernel` returns a session-owned, cached one instead.
    """

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    if resolved.receipts.field_reconstruction is not None and resolved.receipts.field_reconstruction.upper_fourier_switch_active:
        warnings.warn(
            "The simulator uses upper Fourier evaluation near H; FieldKernel reproduces only the FMM component.",
            RuntimeWarning, stacklevel=2,
        )
    context.require_triangle_source_model()
    triangles = context.require_triangles()
    charges = context.charges_at(step)
    options = _options_from_result(
        resolved,
        periodic2=periodic2,
        theta=theta,
        leaf_max=leaf_max,
        order=order,
        config_path=config_path,
        context=context,
    )
    return FieldKernel(
        triangles,
        charges,
        options=options,
        library_path=library_path,
    )


def calc_object_forces_kernel(
    result: RunResult | object,
    *,
    step: int | None = -1,
    target_mesh_ids: int | Iterable[int] | None = None,
    periodic2: Mapping[str, object] | None = None,
    theta: float | None = None,
    leaf_max: int | None = None,
    order: int | None = None,
    config_path: str | Path | None = None,
    library_path: str | Path | None = None,
    kernel: FieldKernel | None = None,
) -> tuple[KernelObjectForceRecord, ...]:
    """Compute object-wise net force with the native field kernel (heavy).

    For each target object, its own source charges are zeroed before evaluating
    the target panel lattice is integrated with seventh-order Gauss-Duffy
    quadrature. A configured uniform ``sim.e0`` is added to the target field,
    matching the simulator pusher semantics while avoiding central-cell
    primary self-force contamination.

    ``kernel`` is an open kernel built from this run, for example
    ``Beach.kernel()``. Its plan is used as given, so ``periodic2``, ``theta``,
    ``leaf_max``, ``order``, and ``library_path`` do not apply; its source charges
    are replaced by the ``step`` charges during the call and restored afterwards,
    and it is not closed. ``None`` builds a kernel for this call and closes it.
    """

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    context.require_triangle_source_model()
    triangles = context.require_triangles()
    centers = _triangle_centers(triangles)
    charges = context.charges_at(step)
    mesh_ids = context.mesh_ids_or_default()
    available_ids = tuple(int(v) for v in np.unique(mesh_ids))
    if target_mesh_ids is None:
        target_ids = available_ids
    elif isinstance(target_mesh_ids, (int, np.integer)):
        target_ids = (int(target_mesh_ids),)
    else:
        target_ids = tuple(dict.fromkeys(int(v) for v in target_mesh_ids))
    missing = [mid for mid in target_ids if mid not in available_ids]
    if missing:
        raise ValueError(
            f"unknown mesh id(s): {missing}. available={list(available_ids)}"
        )

    operation = "calc_object_forces_kernel"
    if kernel is None:
        _require_total_field_config(context, operation=operation)
        options = _options_from_result(
            resolved,
            periodic2=periodic2,
            theta=theta,
            leaf_max=leaf_max,
            order=order,
            config_path=config_path,
            context=context,
        )
        _require_total_field_reconstruction(context, options, operation=operation)
        evaluator = FieldKernel(
            triangles, charges, options=options, library_path=library_path
        )
    else:
        _require_run_kernel(context, kernel, operation=operation)
        evaluator = _restoring_charges(kernel)
    with evaluator as kernel:
        return _object_force_records(
            kernel,
            triangles=triangles,
            centers=centers,
            charges=charges,
            mesh_ids=mesh_ids,
            target_ids=target_ids,
            step=step,
        )


def _object_force_records(
    kernel: FieldKernel,
    *,
    triangles: np.ndarray,
    centers: np.ndarray,
    charges: np.ndarray,
    mesh_ids: np.ndarray,
    target_ids: Iterable[int],
    step: int | None,
) -> tuple[KernelObjectForceRecord, ...]:
    """Net force/torque on each target object from all other source charges.

    The target's own charges are zeroed in ``kernel`` and its panels integrated
    with seventh-order Gauss-Duffy quadrature about the target centroid. The
    caller owns ``kernel`` and restores its charges.
    """

    records: list[KernelObjectForceRecord] = []
    for mesh_id in target_ids:
        mask = mesh_ids == mesh_id
        if not np.any(mask):
            continue
        source_q = np.asarray(charges, dtype=np.float64).copy()
        source_q[mask] = 0.0
        kernel.update_charges(source_q)
        target_q = np.asarray(charges[mask], dtype=np.float64)
        target_points, target_weights, _ = panel_target_quadrature(
            triangles[mask],
            target_q,
            7,
        )
        center = np.asarray(centers[mask].mean(axis=0), dtype=np.float64)
        force, torque = kernel.force_on_charges(
            target_points,
            target_weights,
            origin=center,
        )
        records.append(
            KernelObjectForceRecord(
                mesh_id=mesh_id,
                step=step,
                total_charge_C=float(np.sum(target_q)),
                center_m=center,
                force_N=force,
                torque_Nm=torque,
            )
        )
    return tuple(records)


def field_kernel_options_from_result(
    result: RunResult | object,
    *,
    periodic2: Mapping[str, object] | None = None,
    theta: float | None = None,
    leaf_max: int | None = None,
    order: int | None = None,
    config_path: str | Path | None = None,
) -> FieldKernelOptions:
    """Resolve field-kernel options from a BEACH result and optional config (light)."""

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    context.require_triangle_source_model()
    return _options_from_result(
        resolved,
        periodic2=periodic2,
        theta=theta,
        leaf_max=leaf_max,
        order=order,
        config_path=config_path,
        context=context,
    )


def _options_from_result(
    resolved: RunResult,
    *,
    periodic2: Mapping[str, object] | None,
    theta: float | None,
    leaf_max: int | None,
    order: int | None,
    config_path: str | Path | None,
    context: RunContext | None = None,
) -> FieldKernelOptions:
    run_context = context or RunContext.from_value(resolved, config_path=config_path)
    if (
        periodic2 is None
        and _periodic_nonzero_mode_backend(run_context)
        == "panel_spectral_reference"
    ):
        raise ValueError(
            "field-kernel options cannot automatically translate a "
            "saved/configured periodic2.nonzero_mode_backend="
            "panel_spectral_reference model to a "
            "finite-image plan."
        )
    sim = (
        None
        if run_context.sim is None
        else _sim_with_separated_boundaries(
            run_context.sim,
            domain=run_context.domain,
            field_boundary=run_context.field_boundary,
        )
    )
    periodic_cfg = coerce_periodic2(periodic2, allow_cached_kneq0=True)
    if periodic_cfg is None and sim is not None:
        allow_historical_root_oracle = run_context.requested_config_path is None
        periodic_cfg = coerce_periodic2(
            periodic2_from_sim(
                sim,
                allow_cached_kneq0=True,
                allow_historical_root_oracle=allow_historical_root_oracle,
            ),
            allow_cached_kneq0=True,
            allow_historical_root_oracle=allow_historical_root_oracle,
        )
    resolved_field_solver = _resolved_field_solver(run_context)
    if resolved_field_solver == "treecode":
        raise ValueError(
            "field reconstruction cannot automatically reproduce the simulator's "
            "resolved treecode backend; use simulator-written field/potential output "
            "or rerun with direct or fmm."
        )
    if periodic_cfg is not None:
        # An explicit or configured finite-periodic model is a native FMM model.
        # panel_spectral_reference was rejected above instead of being translated.
        resolved_field_solver = "fmm"
    receipt = run_context.result.receipts.field_reconstruction
    receipt_order = 4 if receipt is None else receipt.fmm_expansion_order
    resolved_order = receipt_order if order is None else int(order)
    resolved_theta = float(
        theta if theta is not None else (sim or {}).get("tree_theta", 0.5)
    )
    resolved_leaf_max = int(
        leaf_max if leaf_max is not None else (sim or {}).get("tree_leaf_max", 16)
    )
    box_min: tuple[float, float, float] | None = None
    box_max: tuple[float, float, float] | None = None
    if sim is not None and "box_min" in sim and "box_max" in sim:
        box_min = tuple(float(v) for v in sim["box_min"])  # type: ignore[index]
        box_max = tuple(float(v) for v in sim["box_max"])  # type: ignore[index]
    return FieldKernelOptions(
        resolved_field_solver=resolved_field_solver,
        external_e0=_external_e0_from_sim(sim),
        theta=resolved_theta,
        leaf_max=resolved_leaf_max,
        order=resolved_order,
        periodic2=periodic_cfg,
        box_min=box_min,
        box_max=box_max,
        periodic_cache_dir=str(
            (sim or {}).get("field_periodic_cache_dir") or default_periodic_cache_dir()
        ),
        periodic_generation_tolerance=float(
            (sim or {}).get("field_periodic_generation_tolerance", 1.0e-8)
        ),
    )


def _require_run_kernel(
    context: RunContext,
    kernel: FieldKernel,
    *,
    operation: str,
) -> None:
    """Check that a caller's kernel holds this run's panels and the total field."""

    nsource = kernel.source_count
    if nsource != context.result.summary.mesh_nelem:
        raise ValueError(
            f"{operation}: the kernel holds {nsource} source panels but the run has "
            f"{context.result.summary.mesh_nelem}; build it from this run with Beach.kernel() "
            "or field_kernel_from_result()."
        )
    _require_total_field_reconstruction(context, kernel.options, operation=operation)


@contextmanager
def _restoring_charges(kernel: FieldKernel) -> Iterator[FieldKernel]:
    """Lend a caller's kernel and put its source charges back afterwards."""

    charges = kernel.source_charges
    try:
        yield kernel
    finally:
        kernel.update_charges(charges)


def _require_total_field_reconstruction(
    context: RunContext,
    options: FieldKernelOptions,
    *,
    operation: str,
) -> None:
    """Reject component-only kernels where a public API promises total fields."""

    _require_total_field_config(context, operation=operation)
    _require_complete_total_kernel(options, operation=operation)


def _require_total_field_config(
    context: RunContext,
    *,
    operation: str,
) -> None:
    """Validate configuration before resolving a total-field kernel."""

    receipt = context.result.receipts.field_reconstruction
    if receipt is not None:
        if receipt.resolved_field_solver == "treecode":
            raise ValueError(
                f"{operation} cannot automatically reproduce the simulator's "
                "resolved treecode backend; use simulator-written field/potential "
                "output or rerun with direct or fmm."
            )
        if receipt.periodic_nonzero_mode_backend == "panel_spectral_reference":
            raise ValueError(
                f"{operation} cannot reconstruct the simulator's "
                "panel_spectral_reference split-periodic field with the native "
                "finite-image kernel; use simulator-written potential output or "
                "the dedicated periodic reference validation workflow."
            )
        return
    config = context.config
    if config is None:
        raise ValueError(
            f"{operation} requires the run's beach.toml to reconstruct the "
            "simulator total field; pass config_path or keep beach.toml near "
            "the output directory."
        )
    if _periodic_nonzero_mode_backend(context) == "panel_spectral_reference":
        raise ValueError(
            f"{operation} cannot reconstruct the simulator's "
            "panel_spectral_reference split-periodic field with the native "
            "finite-image kernel; use simulator-written potential output or "
            "the dedicated periodic reference validation workflow."
        )
    if _resolved_field_solver(context) == "treecode":
        raise ValueError(
            f"{operation} cannot automatically reproduce the simulator's resolved "
            "treecode backend; use simulator-written field/potential output or "
            "rerun with direct or fmm."
        )
    sim = context.sim
    if (
        context.requested_config_path is not None
        and sim is not None
        and str(sim.get("field_periodic_far_correction", "")).strip().lower()
        == "m2l_root_oracle"
    ):
        raise ValueError(
            'periodic2.far_correction "m2l_root_oracle" was removed; use '
            '"none" for total-field post-processing.'
        )


def _periodic_nonzero_mode_backend(context: RunContext) -> str | None:
    receipt = context.result.receipts.field_reconstruction
    if receipt is not None:
        return receipt.periodic_nonzero_mode_backend
    config = context.config
    if config is None:
        return None
    periodic2 = config.get("periodic2")
    if not isinstance(periodic2, Mapping):
        return None
    value = periodic2.get("nonzero_mode_backend")
    return None if value is None else str(value).strip().lower()


def _resolved_field_solver(context: RunContext) -> str:
    """Resolve the simulator backend for automatic result reconstruction."""

    receipt = context.result.receipts.field_reconstruction
    if receipt is not None:
        return _coerce_resolved_field_solver(receipt.resolved_field_solver)
    sim = context.sim
    if sim is None:
        # Low-level field_kernel_from_result historically permits a result without
        # its config. High-level total-field APIs reject that case separately.
        return "fmm"
    requested = str(sim.get("field_solver", "auto")).strip().lower()
    if requested in {"direct", "treecode", "fmm"}:
        return requested
    if requested != "auto":
        raise ValueError(
            "sim.field_solver must be direct, treecode, fmm, or auto for field "
            "reconstruction."
        )
    tree_min_nelem = int(sim.get("tree_min_nelem", 256))
    if tree_min_nelem < 1:
        raise ValueError("sim.tree_min_nelem must be >= 1 for field reconstruction.")
    return "fmm" if context.result.summary.mesh_nelem >= tree_min_nelem else "direct"


def _require_complete_total_kernel(
    options: FieldKernelOptions,
    *,
    operation: str,
) -> None:
    """Reject a native kernel option set that represents only one component."""

    resolved_field_solver = _coerce_resolved_field_solver(
        options.resolved_field_solver
    )
    if resolved_field_solver == "treecode":
        raise ValueError(
            f"{operation} cannot automatically reproduce the simulator's resolved "
            "treecode backend."
        )
    periodic_cfg = coerce_periodic2(
        options.periodic2,
        allow_cached_kneq0=True,
    )
    if resolved_field_solver == "direct" and periodic_cfg is not None:
        raise ValueError(
            f"{operation} cannot combine exact-direct reconstruction with a "
            "periodic field plan."
        )
    if periodic_cfg is not None and periodic_cfg.far_correction == "cached_kneq0":
        raise ValueError(
            f"{operation} cannot treat cached_kneq0 as the simulator total "
            "field: the native kernel returns only the k!=0 component and does "
            "not compose the physical periodic zero mode. Use saved "
            "mesh_potential.csv where applicable or ObjectInteractionSnapshot "
            "for supported force workflows."
        )


def _external_e0_from_sim(
    sim: Mapping[str, object] | None,
) -> tuple[float, float, float]:
    if sim is None:
        return (0.0, 0.0, 0.0)
    has_vector = "e0" in sim
    has_abs = "e0_abs" in sim
    has_phi_xy = "e0_phi_xy_deg" in sim
    has_phi_z = "e0_phi_z_deg" in sim
    if has_vector and (has_abs or has_phi_xy or has_phi_z):
        raise ValueError(
            "sim.e0 cannot be combined with sim.e0_abs/e0_phi_xy_deg/e0_phi_z_deg."
        )
    if has_vector:
        e0 = _vec3(sim["e0"], name="e0")
        return tuple(float(v) for v in e0)
    if (has_phi_xy or has_phi_z) and not has_abs:
        raise ValueError("sim.e0_phi_xy_deg/e0_phi_z_deg require sim.e0_abs.")
    if not has_abs:
        return (0.0, 0.0, 0.0)

    e0_abs = float(sim.get("e0_abs", 0.0))
    phi_xy = np.deg2rad(float(sim.get("e0_phi_xy_deg", 0.0)))
    phi_z = np.deg2rad(float(sim.get("e0_phi_z_deg", 0.0)))
    if not np.isfinite(e0_abs) or e0_abs < 0.0:
        raise ValueError("sim.e0_abs must be finite and >= 0.")
    if not np.isfinite(phi_xy) or not np.isfinite(phi_z):
        raise ValueError("sim.e0_phi_xy_deg/e0_phi_z_deg must be finite.")
    return (
        float(e0_abs * np.cos(phi_z) * np.cos(phi_xy)),
        float(e0_abs * np.cos(phi_z) * np.sin(phi_xy)),
        float(e0_abs * np.sin(phi_z)),
    )

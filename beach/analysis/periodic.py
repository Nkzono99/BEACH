"""Periodic2 settings resolved from a run's receipt or configuration."""

from __future__ import annotations

from typing import Mapping

from beach.io.types import FortranRunResult
from beach.native.periodic2 import (
    Periodic2Config,
    _coerce_vec3,
    coerce_periodic2,
    normalize_periodic2_far_correction,
)

from .context import RunContext


def auto_periodic2_from_result(
    resolved: FortranRunResult,
    *,
    allow_cached_kneq0: bool = False,
    context: RunContext | None = None,
) -> Periodic2Config | None:
    """Load and normalize periodic2 settings associated with one result."""

    run_context = context or RunContext.from_value(resolved)
    if run_context.sim is None:
        return None
    allow_historical_root_oracle = run_context.requested_config_path is None
    sim = _sim_with_separated_boundaries(
        run_context.sim,
        domain=run_context.domain,
        field_boundary=run_context.field_boundary,
    )
    periodic2 = periodic2_from_sim(
        sim,
        allow_cached_kneq0=allow_cached_kneq0,
        allow_historical_root_oracle=allow_historical_root_oracle,
    )
    return coerce_periodic2(
        periodic2,
        allow_cached_kneq0=allow_cached_kneq0,
        allow_historical_root_oracle=allow_historical_root_oracle,
    )



def _sim_with_separated_boundaries(
    sim: Mapping[str, object],
    *,
    domain: Mapping[str, object] | None,
    field_boundary: Mapping[str, object] | None,
) -> dict[str, object]:
    """Compose the legacy post-processing view without changing public TOML."""

    resolved = dict(sim)
    if field_boundary is not None:
        resolved["field_bc_mode"] = field_boundary.get("mode", "free")
    if domain is not None:
        if "box_min" in domain:
            resolved["box_min"] = domain["box_min"]
        if "box_max" in domain:
            resolved["box_max"] = domain["box_max"]
        periodic_axes = domain.get("periodic_axes", [])
        for axis in ("x", "y", "z"):
            mode = "periodic" if axis in periodic_axes else "open"
            resolved[f"bc_{axis}_low"] = mode
            resolved[f"bc_{axis}_high"] = mode
    return resolved



def periodic2_from_sim(
    sim: Mapping[str, object],
    *,
    allow_cached_kneq0: bool = False,
    allow_historical_root_oracle: bool = False,
) -> dict[str, object] | None:
    """Translate a simulator table into the external periodic2 mapping shape."""

    field_bc_mode = str(sim.get("field_bc_mode", "free")).strip().lower()
    if field_bc_mode != "periodic2":
        return None

    if "box_min" not in sim or "box_max" not in sim:
        raise ValueError(
            'periodic2 potential requires "sim.box_min" and "sim.box_max".'
        )
    box_min = _coerce_vec3(sim["box_min"], name="sim.box_min")
    box_max = _coerce_vec3(sim["box_max"], name="sim.box_max")

    periodic_axes: list[int] = []
    for axis_idx, axis_name in enumerate(("x", "y", "z")):
        low = _canonical_boundary_mode(sim.get(f"bc_{axis_name}_low", "open"))
        high = _canonical_boundary_mode(sim.get(f"bc_{axis_name}_high", "open"))
        if (low == "periodic") != (high == "periodic"):
            raise ValueError("periodic2 requires bc_low(axis)=bc_high(axis)=periodic.")
        if low == "periodic":
            periodic_axes.append(axis_idx)

    if len(periodic_axes) != 2:
        raise ValueError(
            'sim.field_bc_mode="periodic2" requires exactly two periodic axes.'
        )

    lengths = [box_max[axis] - box_min[axis] for axis in periodic_axes]
    if lengths[0] <= 0.0 or lengths[1] <= 0.0:
        raise ValueError("periodic2 requires positive box length on periodic axes.")

    ewald_layers = int(sim.get("field_periodic_ewald_layers", 4))
    far_correction, ewald_layers = normalize_periodic2_far_correction(
        sim.get("field_periodic_far_correction", "none"),
        ewald_layers=ewald_layers,
        allow_cached_kneq0=allow_cached_kneq0,
        allow_historical_root_oracle=allow_historical_root_oracle,
    )

    return {
        "axes": tuple(periodic_axes),
        "lengths": tuple(lengths),
        "origins": tuple(box_min[axis] for axis in periodic_axes),
        "image_layers": int(sim.get("field_periodic_image_layers", 1)),
        "far_correction": far_correction,
        "ewald_alpha": float(sim.get("field_periodic_ewald_alpha", 0.0)),
        "ewald_layers": ewald_layers,
    }



def _canonical_boundary_mode(value: object) -> str:
    mode = str(value).strip().lower()
    if mode in {"open", "outflow", "escape"}:
        return "open"
    return mode

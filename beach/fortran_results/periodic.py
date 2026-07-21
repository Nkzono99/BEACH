"""Typed periodic2 configuration parsing shared by result utilities."""

from __future__ import annotations

import math
from typing import Mapping, NamedTuple, TypeAlias

from .context import RunContext
from .types import FortranRunResult


class Periodic2Config(NamedTuple):
    """Validated periodic2 settings with legacy tuple compatibility."""

    axes: tuple[int, int]
    lengths: tuple[float, float]
    origins: tuple[float, float]
    image_layers: int
    far_correction: str
    ewald_alpha: float
    ewald_layers: int


Periodic2Tuple: TypeAlias = tuple[
    tuple[int, int],
    tuple[float, float],
    tuple[float, float],
    int,
    str,
    float,
    int,
]
Periodic2Input: TypeAlias = Mapping[str, object] | Periodic2Tuple | Periodic2Config


def coerce_periodic2(
    periodic2: Periodic2Input | None,
    *,
    allow_cached_kneq0: bool = False,
) -> Periodic2Config | None:
    """Validate external periodic2 input and return its canonical typed form."""

    if periodic2 is None:
        return None
    if isinstance(periodic2, tuple) and len(periodic2) == 7:
        (
            axes,
            lengths,
            origins,
            image_layers,
            far_correction,
            ewald_alpha,
            ewald_layers,
        ) = periodic2
        periodic2 = {
            "axes": axes,
            "lengths": lengths,
            "origins": origins,
            "image_layers": image_layers,
            "far_correction": far_correction,
            "ewald_alpha": ewald_alpha,
            "ewald_layers": ewald_layers,
        }
    if not isinstance(periodic2, Mapping):
        raise ValueError("periodic2 must be a mapping or None.")

    if "axes" not in periodic2 or "lengths" not in periodic2:
        raise ValueError('periodic2 requires "axes" and "lengths".')
    axes_obj = periodic2["axes"]
    lengths_obj = periodic2["lengths"]
    if not isinstance(axes_obj, (list, tuple)) or len(axes_obj) != 2:
        raise ValueError("periodic2.axes must be a length-2 sequence.")
    if not isinstance(lengths_obj, (list, tuple)) or len(lengths_obj) != 2:
        raise ValueError("periodic2.lengths must be a length-2 sequence.")

    axes = (int(axes_obj[0]), int(axes_obj[1]))
    if axes[0] == axes[1] or any(axis < 0 or axis > 2 for axis in axes):
        raise ValueError(
            "periodic2.axes must contain two distinct axis indices in {0,1,2}."
        )

    lengths = (float(lengths_obj[0]), float(lengths_obj[1]))
    if any((not math.isfinite(length)) or length <= 0.0 for length in lengths):
        raise ValueError("periodic2.lengths must be finite and positive.")

    if "origins" in periodic2:
        origins_obj = periodic2["origins"]
        if not isinstance(origins_obj, (list, tuple)) or len(origins_obj) != 2:
            raise ValueError("periodic2.origins must be a length-2 sequence.")
        origins = (float(origins_obj[0]), float(origins_obj[1]))
    elif "box_min" in periodic2:
        box_min = _coerce_vec3(periodic2["box_min"], name="periodic2.box_min")
        origins = (box_min[axes[0]], box_min[axes[1]])
    else:
        origins = (0.0, 0.0)
    if any(not math.isfinite(origin) for origin in origins):
        raise ValueError("periodic2.origins must be finite.")

    image_layers = int(periodic2.get("image_layers", 1))
    if image_layers < 0:
        raise ValueError("periodic2.image_layers must be >= 0.")

    ewald_layers = int(periodic2.get("ewald_layers", 4))
    if ewald_layers < 0:
        raise ValueError("periodic2.ewald_layers must be >= 0.")

    far_correction, ewald_layers = normalize_periodic2_far_correction(
        periodic2.get("far_correction", "none"),
        ewald_layers=ewald_layers,
        allow_cached_kneq0=allow_cached_kneq0,
    )

    ewald_alpha = float(periodic2.get("ewald_alpha", 0.0))
    if (not math.isfinite(ewald_alpha)) or ewald_alpha < 0.0:
        raise ValueError("periodic2.ewald_alpha must be finite and >= 0.")
    if far_correction == "m2l_root_oracle" and ewald_layers < 1:
        raise ValueError("periodic2.ewald_layers must be >= 1 for m2l_root_oracle.")

    return Periodic2Config(
        axes=axes,
        lengths=lengths,
        origins=origins,
        image_layers=image_layers,
        far_correction=far_correction,
        ewald_alpha=ewald_alpha,
        ewald_layers=ewald_layers,
    )


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
    periodic2 = periodic2_from_sim(
        run_context.sim,
        allow_cached_kneq0=allow_cached_kneq0,
    )
    return coerce_periodic2(periodic2, allow_cached_kneq0=allow_cached_kneq0)


def periodic2_from_sim(
    sim: Mapping[str, object],
    *,
    allow_cached_kneq0: bool = False,
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


def normalize_periodic2_far_correction(
    raw_far_correction: object,
    *,
    ewald_layers: int,
    allow_cached_kneq0: bool = False,
) -> tuple[str, int]:
    """Normalize the configured far-correction policy."""

    far_correction = str(raw_far_correction).strip().lower()
    allowed = {"auto", "none", "m2l_root_oracle"}
    if allow_cached_kneq0:
        allowed.add("cached_kneq0")
    if far_correction not in allowed:
        expected = '"auto", "none", or "m2l_root_oracle"'
        if allow_cached_kneq0:
            expected = '"auto", "none", "m2l_root_oracle", or "cached_kneq0"'
        raise ValueError(f"periodic2.far_correction must be {expected}.")
    if far_correction == "auto":
        return "none", ewald_layers
    return far_correction, ewald_layers


def _coerce_vec3(value: object, *, name: str) -> tuple[float, float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError(f"{name} must contain exactly 3 values.")
    try:
        return float(value[0]), float(value[1]), float(value[2])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain numeric values.") from exc


def _canonical_boundary_mode(value: object) -> str:
    mode = str(value).strip().lower()
    if mode in {"open", "outflow", "escape"}:
        return "open"
    return mode

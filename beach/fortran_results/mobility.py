"""Coulomb mobility analysis helpers."""

from __future__ import annotations

from math import pi
from pathlib import Path
from typing import Iterable

import numpy as np

from .coulomb import calc_coulomb
from .mesh import _triangle_centers
from .objects import normalize_kind_filter, resolve_object_specs
from .selection import _build_mesh_selection, _resolve_result
from .types import (
    CoulombMobilityAnalysis,
    CoulombMobilityRecord,
    FortranRunResult,
)


def analyze_coulomb_mobility(
    result: FortranRunResult | object,
    *,
    step: int | None = -1,
    softening: float = 0.0,
    config_path: str | Path | None = None,
    gravity: Iterable[float] = (0.0, 0.0, -9.81),
    support_normal: Iterable[float] | None = None,
    support_kinds: Iterable[str] | None = ("plane",),
    target_kinds: Iterable[str] | None = None,
    density_kg_m3: float | None = None,
    mu_static: float | None = None,
    mu_roll: float | None = None,
    adhesion_force_N: float = 0.0,
) -> CoulombMobilityAnalysis:
    """Analyze per-object Coulomb mobility indicators.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.
    step : int or None, default -1
        History batch step used for charge snapshot.
        ``None`` uses final charges from ``charges.csv``.
    softening : float, default 0.0
        Softening length used in Coulomb-force evaluation.
    config_path : str, pathlib.Path, or None, default None
        Optional ``beach.toml`` path used for object labels and geometry metadata.
    gravity : iterable of float, default ``(0, 0, -9.81)``
        Gravity vector in meters per second squared.
    support_normal : iterable of float or None, default None
        Support-normal direction. ``None`` uses ``-gravity`` when possible,
        otherwise ``+z``.
    support_kinds : iterable of str or None, default ``("plane",)``
        Object kinds treated as fixed supports and excluded from default targets.
    target_kinds : iterable of str or None, default None
        Explicit target kinds to analyze. ``None`` means all non-support objects,
        falling back to all objects when no non-support object exists.
    density_kg_m3 : float or None, default None
        Bulk density used to estimate weight for volumetric templates.
    mu_static : float or None, default None
        Static-friction coefficient for slide-ratio estimation.
    mu_roll : float or None, default None
        Rolling-resistance coefficient for roll-ratio estimation.
    adhesion_force_N : float, default 0.0
        Additional normal adhesion resisting lift.

    Returns
    -------
    CoulombMobilityAnalysis
        Per-object force/torque summary and mobility indicators.

    Raises
    ------
    ValueError
        If arguments are invalid or no target objects are available.
    """

    resolved = _resolve_result(result)
    if softening < 0.0:
        raise ValueError("softening must be >= 0.")
    if density_kg_m3 is not None and density_kg_m3 <= 0.0:
        raise ValueError("density_kg_m3 must be > 0 when specified.")
    if mu_static is not None and mu_static < 0.0:
        raise ValueError("mu_static must be >= 0 when specified.")
    if mu_roll is not None and mu_roll < 0.0:
        raise ValueError("mu_roll must be >= 0 when specified.")
    if adhesion_force_N < 0.0:
        raise ValueError("adhesion_force_N must be >= 0.")

    gravity_vec = _coerce_vec3(gravity, name="gravity")
    support_normal_vec = _resolve_support_normal(gravity_vec, support_normal)
    support_kind_filter = normalize_kind_filter(support_kinds) or set()
    target_kind_filter = normalize_kind_filter(target_kinds)
    object_specs = resolve_object_specs(resolved, config_path=config_path)

    if target_kind_filter is None:
        target_specs = [
            spec for spec in object_specs if spec.kind not in support_kind_filter
        ]
        if len(target_specs) == 0:
            target_specs = list(object_specs)
    else:
        target_specs = [
            spec for spec in object_specs if spec.kind in target_kind_filter
        ]

    if len(target_specs) == 0:
        available_kinds = sorted({spec.kind for spec in object_specs})
        raise ValueError(
            "target_kinds did not match any objects. "
            f"available={available_kinds}"
        )

    records: list[CoulombMobilityRecord] = []
    gravity_support_mag = max(0.0, float(-np.dot(gravity_vec, support_normal_vec)))
    gravity_mag = float(np.linalg.norm(gravity_vec))

    for spec in target_specs:
        source_mesh_ids = [
            other.mesh_id for other in object_specs if other.mesh_id != spec.mesh_id
        ]
        if source_mesh_ids:
            interaction = calc_coulomb(
                resolved,
                target=spec.mesh_id,
                source=source_mesh_ids,
                step=step,
                softening=softening,
                torque_origin="target_center",
            )
            force = interaction.force_on_a_N
            torque = interaction.torque_on_a_Nm
            record_step = interaction.step
        else:
            force = np.zeros(3, dtype=float)
            torque = np.zeros(3, dtype=float)
            record_step = step

        selection = _build_mesh_selection(resolved, (spec.mesh_id,), step=step)
        center = _triangle_centers(selection.triangles).mean(axis=0)

        force_normal = float(np.dot(force, support_normal_vec))
        force_tangent_vec = force - force_normal * support_normal_vec
        force_tangent = float(np.linalg.norm(force_tangent_vec))
        torque_normal = float(np.dot(torque, support_normal_vec))
        torque_tangent_vec = torque - torque_normal * support_normal_vec
        torque_tangent = float(np.linalg.norm(torque_tangent_vec))

        notes: list[str] = []
        mass_kg, radius_m, geom_notes = _estimate_mass_and_radius(
            spec.kind,
            spec.template,
            density_kg_m3=density_kg_m3,
        )
        notes.extend(geom_notes)

        weight_support: float | None = None
        resisting_normal: float | None = None
        effective_normal_load: float | None = None
        lift_ratio: float | None = None
        slide_ratio: float | None = None
        roll_ratio: float | None = None

        can_estimate_weight = mass_kg is not None or gravity_mag <= np.finfo(float).tiny
        if can_estimate_weight:
            if mass_kg is None:
                weight_support = 0.0
            else:
                weight_support = mass_kg * gravity_support_mag
            resisting_normal = weight_support + adhesion_force_N
            effective_normal_load = max(resisting_normal - force_normal, 0.0)
            lift_ratio = _safe_ratio(max(force_normal, 0.0), resisting_normal)

            if mu_static is not None:
                slide_limit = mu_static * effective_normal_load
                slide_ratio = _safe_ratio(force_tangent, slide_limit)

            if mu_roll is not None:
                if radius_m is None:
                    notes.append("roll_ratio requires characteristic radius")
                else:
                    roll_limit = mu_roll * effective_normal_load * radius_m
                    roll_ratio = _safe_ratio(torque_tangent, roll_limit)
        else:
            notes.append("mass estimate unavailable; lift/slide/roll ratios skipped")

        records.append(
            CoulombMobilityRecord(
                mesh_id=spec.mesh_id,
                label=spec.label,
                kind=spec.kind,
                step=record_step,
                center_m=center,
                force_N=force,
                torque_Nm=torque,
                force_normal_N=force_normal,
                force_tangent_N=force_tangent,
                torque_normal_Nm=torque_normal,
                torque_tangent_Nm=torque_tangent,
                mass_kg=mass_kg,
                characteristic_radius_m=radius_m,
                weight_support_N=weight_support,
                resisting_normal_N=resisting_normal,
                effective_normal_load_N=effective_normal_load,
                lift_ratio=lift_ratio,
                slide_ratio=slide_ratio,
                roll_ratio=roll_ratio,
                notes=tuple(dict.fromkeys(notes)),
            )
        )

    return CoulombMobilityAnalysis(
        step=records[0].step if records else step,
        softening=float(softening),
        gravity_m_s2=gravity_vec,
        support_normal_m=support_normal_vec,
        support_kinds=tuple(sorted(support_kind_filter)),
        density_kg_m3=density_kg_m3,
        mu_static=mu_static,
        mu_roll=mu_roll,
        adhesion_force_N=float(adhesion_force_N),
        records=tuple(records),
    )


def _coerce_vec3(value: Iterable[float], *, name: str) -> np.ndarray:
    vector = np.asarray(list(value), dtype=float)
    if vector.shape != (3,):
        raise ValueError(f"{name} must contain exactly 3 values.")
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain finite values.")
    return vector


def _resolve_support_normal(
    gravity: np.ndarray,
    support_normal: Iterable[float] | None,
) -> np.ndarray:
    if support_normal is not None:
        normal = _coerce_vec3(support_normal, name="support_normal")
    else:
        gravity_mag = float(np.linalg.norm(gravity))
        if gravity_mag > np.finfo(float).tiny:
            normal = -gravity
        else:
            normal = np.array([0.0, 0.0, 1.0], dtype=float)

    normal_mag = float(np.linalg.norm(normal))
    if normal_mag <= np.finfo(float).tiny:
        raise ValueError("support_normal must not be the zero vector.")
    return normal / normal_mag


def _estimate_mass_and_radius(
    kind: str,
    template: object,
    *,
    density_kg_m3: float | None,
) -> tuple[float | None, float | None, list[str]]:
    notes: list[str] = []
    if density_kg_m3 is None:
        notes.append("density_kg_m3 is not specified")
        return None, None, notes
    if not isinstance(template, dict):
        notes.append("template geometry metadata is unavailable")
        return None, None, notes

    kind_key = str(kind).strip().lower()
    if kind_key == "sphere":
        radius = _read_positive_float(template.get("radius"))
        if radius is None:
            notes.append("sphere radius is unavailable")
            return None, None, notes
        volume = 4.0 * pi * radius * radius * radius / 3.0
        return density_kg_m3 * volume, radius, notes

    if kind_key == "box":
        size = _read_vec(template.get("size"), expected=3)
        if size is None:
            notes.append("box size is unavailable")
            return None, None, notes
        volume = float(size[0] * size[1] * size[2])
        return density_kg_m3 * volume, None, notes

    if kind_key == "cylinder":
        radius = _read_positive_float(template.get("radius"))
        height = _read_positive_float(template.get("height"))
        if radius is None or height is None:
            notes.append("cylinder radius/height is unavailable")
            return None, None, notes
        volume = pi * radius * radius * height
        return density_kg_m3 * volume, radius, notes

    notes.append(f"mass estimate is not implemented for kind={kind_key}")
    return None, None, notes


def _read_positive_float(value: object) -> float | None:
    try:
        scalar = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(scalar) or scalar <= 0.0:
        return None
    return scalar


def _read_vec(value: object, *, expected: int) -> np.ndarray | None:
    if not isinstance(value, (list, tuple)) or len(value) != expected:
        return None
    try:
        vector = np.asarray([float(item) for item in value], dtype=float)
    except (TypeError, ValueError):
        return None
    if not np.all(np.isfinite(vector)) or np.any(vector <= 0.0):
        return None
    return vector


def _safe_ratio(numerator: float, denominator: float | None) -> float | None:
    if denominator is None:
        return None
    if denominator <= np.finfo(float).tiny:
        if numerator <= np.finfo(float).tiny:
            return 0.0
        return float("inf")
    return float(numerator / denominator)

"""Validate particle sources, injection, and particle boundaries."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

from ._shared import (
    _FACE_KEYS,
    _RESERVOIR_SOURCE_MODES,
    ConfigValidationError,
    _maybe_vec3,
)


def _validate_grid_flux_keys(species_table: Mapping[str, Any], *, index: int) -> None:
    has_flux = "particle_flux_m2_s" in species_table
    has_current = "current_density_a_m2" in species_table
    if has_flux and has_current:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] cannot define both "
            "particle_flux_m2_s and current_density_a_m2."
        )
    if not has_flux and not has_current:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            'velocity_distribution="grid" and requires particle_flux_m2_s '
            "or current_density_a_m2."
        )
    if has_flux:
        value = species_table["particle_flux_m2_s"]
        if (
            not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            or float(value) <= 0.0
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].particle_flux_m2_s "
                "must be finite and > 0."
            )
    if has_current:
        value = species_table["current_density_a_m2"]
        if (
            not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            or float(value) == 0.0
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].current_density_a_m2 "
                "must be finite and non-zero."
            )
        q_particle = species_table.get("q_particle", -1.602176634e-19)
        if (
            not isinstance(q_particle, (int, float))
            or not math.isfinite(float(q_particle))
            or float(q_particle) == 0.0
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].q_particle "
                "must be finite and non-zero for current_density_a_m2."
            )


def _validate_velocity_grid_forbidden(
    species_table: Mapping[str, Any],
    *,
    index: int,
    source_mode: str,
) -> None:
    if (
        str(species_table.get("velocity_distribution", "maxwellian")).strip().lower()
        != "maxwellian"
    ):
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            f'source_mode="{source_mode}" and cannot define velocity_distribution="grid".'
        )
    if (
        str(species_table.get("velocity_grid_sampling", "auto")).strip().lower()
        != "auto"
    ):
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            f'source_mode="{source_mode}" and cannot define velocity_grid_sampling.'
        )
    for key in ("velocity_grid_path", "particle_flux_m2_s", "current_density_a_m2"):
        if key in species_table:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] uses "
                f'source_mode="{source_mode}" and cannot define {key}.'
            )


def _validate_face_source_common(
    species_table: Mapping[str, Any],
    *,
    index: int,
    source_mode: str,
    use_box: bool,
    batch_duration: float,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> None:
    if not use_box:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            f'source_mode="{source_mode}" and requires [domain].'
        )
    if batch_duration <= 0.0:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            f'source_mode="{source_mode}" and requires batch_duration > 0.'
        )
    inject_face = species_table.get("inject_face")
    if not isinstance(inject_face, str) or inject_face == "":
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            f'source_mode="{source_mode}" and requires inject_face.'
        )
    if box_min is not None and box_max is not None:
        _validate_face_bounds(
            species_table,
            index=index,
            inject_face=inject_face,
            box_min=box_min,
            box_max=box_max,
        )


def _validate_face_bounds(
    species_table: Mapping[str, Any],
    *,
    index: int,
    inject_face: str,
    box_min: Sequence[float],
    box_max: Sequence[float],
) -> None:
    pos_low = _maybe_vec3(
        species_table.get("pos_low"),
        name=f"particles.species[{index}].pos_low",
    )
    pos_high = _maybe_vec3(
        species_table.get("pos_high"),
        name=f"particles.species[{index}].pos_high",
    )
    if pos_low is None or pos_high is None:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] must define pos_low "
            "and pos_high on the inject_face."
        )

    axis_by_face = {
        "x_low": (0, float(box_min[0])),
        "x_high": (0, float(box_max[0])),
        "y_low": (1, float(box_min[1])),
        "y_high": (1, float(box_max[1])),
        "z_low": (2, float(box_min[2])),
        "z_high": (2, float(box_max[2])),
    }
    if inject_face not in axis_by_face:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] has invalid inject_face={inject_face!r}."
        )
    axis, boundary = axis_by_face[inject_face]
    if pos_low[axis] != boundary or pos_high[axis] != boundary:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] pos_low/pos_high must "
            f"lie on inject_face={inject_face!r}."
        )
    for other_axis in range(3):
        if other_axis == axis:
            continue
        if pos_low[other_axis] > pos_high[other_axis]:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] pos_low must be <= "
                "pos_high along the inject-face coordinates."
            )
        low_bound = float(box_min[other_axis])
        high_bound = float(box_max[other_axis])
        if not (low_bound <= pos_low[other_axis] <= high_bound):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] pos_low is outside the box."
            )
        if not (low_bound <= pos_high[other_axis] <= high_bound):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] pos_high is outside the box."
            )


def _validate_species_particle_boundary(
    species: Mapping[str, Any],
    *,
    index: int,
    periodic_axes: set[str],
    global_boundary: Mapping[str, Any],
) -> dict[str, str]:
    face_keys = {
        "x_low",
        "x_high",
        "y_low",
        "y_high",
        "z_low",
        "z_high",
    }
    raw_boundary = species.get("boundary", {})
    if not isinstance(raw_boundary, Mapping):
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}].boundary must be a table."
        )
    unknown = set(raw_boundary) - face_keys
    if unknown:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}].boundary has "
            "unsupported key(s): " + ", ".join(sorted(unknown)) + "."
        )

    effective: dict[str, str] = {}
    for face in face_keys:
        axis = face[0]
        action = raw_boundary.get(face, "inherit")
        if action not in {
            "inherit",
            "open",
            "reflect",
            "redistributed_reflect",
        }:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].boundary."
                f'{face} must be "inherit", "open", "reflect", or '
                '"redistributed_reflect".'
            )
        if axis in periodic_axes:
            if action != "inherit":
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}].boundary."
                    f"{face} cannot override a periodic domain face."
                )
            effective[face] = "periodic"
        elif action == "inherit":
            effective[face] = str(global_boundary.get(face, "open"))
        else:
            effective[face] = str(action)
    return effective


def _validate_species_boundary_inflow(
    species: Mapping[str, Any],
    *,
    index: int,
    periodic_axes: set[str],
    effective_boundary: Mapping[str, str],
) -> tuple[str, ...]:
    raw_inflow = species.get("boundary_inflow", {})
    if not isinstance(raw_inflow, Mapping):
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}].boundary_inflow "
            "must be a table."
        )
    unknown = set(raw_inflow) - _FACE_KEYS
    if unknown:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}].boundary_inflow has "
            "unsupported key(s): " + ", ".join(sorted(unknown)) + "."
        )

    enabled_faces: list[str] = []
    for face in sorted(raw_inflow):
        if raw_inflow[face] != "reservoir":
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].boundary_inflow."
                f'{face} must be "reservoir".'
            )
        if face[0] in periodic_axes:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].boundary_inflow."
                f"{face} cannot inject through a periodic domain face."
            )
        if effective_boundary.get(face, "open") != "open":
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].boundary_inflow."
                f"{face} requires an open particle boundary action."
            )
        enabled_faces.append(face)
    return tuple(enabled_faces)


def _validate_reservoir_injection_common(
    species_table: Mapping[str, Any],
    *,
    index: int,
    use_box: bool,
    batch_duration: float,
) -> None:
    if not use_box:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] reservoir injection "
            "requires [domain]."
        )
    if batch_duration <= 0.0:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] reservoir injection "
            "requires batch_duration > 0."
        )
    source_mode = species_table.get("source_mode", "volume_seed")
    if source_mode != "volume_seed" and species_table.get("boundary_inflow"):
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] cannot combine "
            f'source_mode="{source_mode}" with boundary_inflow.'
        )
    if species_table.get("boundary_inflow"):
        if "source_normal" in species_table:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].source_normal "
                'is only valid for source_mode="plane_source".'
            )
        if "inject_face" in species_table:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].inject_face "
                "is not used by boundary_inflow."
            )


def _validate_reservoir_physics(
    species_table: Mapping[str, Any],
    *,
    index: int,
    source_mode: str,
) -> None:
    has_weight = "w_particle" in species_table
    has_target = "target_macro_particles_per_batch" in species_table
    if has_weight == has_target:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] reservoir injection "
            "requires exactly one of w_particle and target_macro_particles_per_batch."
        )
    if has_weight:
        weight = species_table["w_particle"]
        if (
            not isinstance(weight, (int, float))
            or isinstance(weight, bool)
            or not math.isfinite(float(weight))
            or float(weight) <= 0.0
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].w_particle "
                "must be finite and > 0."
            )
    else:
        target = species_table["target_macro_particles_per_batch"]
        if (
            not isinstance(target, int)
            or isinstance(target, bool)
            or (target < 1 and target != -1)
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}]."
                "target_macro_particles_per_batch must be > 0 or -1."
            )

    if source_mode in _RESERVOIR_SOURCE_MODES:
        if "npcls_per_step" in species_table:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] uses "
                f'source_mode="{source_mode}" and cannot define npcls_per_step.'
            )
        for key in (
            "emit_current_density_a_m2",
            "rays_per_batch",
            "ray_direction",
            "deposit_opposite_charge_on_emit",
        ):
            if key in species_table:
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}] uses "
                    f'source_mode="{source_mode}" and cannot define {key}.'
                )

    velocity_distribution = (
        str(species_table.get("velocity_distribution", "maxwellian")).strip().lower()
    )
    if velocity_distribution == "grid":
        npcls_per_step = species_table.get("npcls_per_step", 0)
        if (
            species_table.get("boundary_inflow")
            and isinstance(npcls_per_step, int)
            and not isinstance(npcls_per_step, bool)
            and npcls_per_step > 0
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] grid "
                "boundary_inflow cannot be combined with positive npcls_per_step."
            )
        if "velocity_grid_path" not in species_table:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] uses "
                'velocity_distribution="grid" and requires velocity_grid_path.'
            )
        _validate_grid_flux_keys(species_table, index=index)
        for key in (
            "number_density_cm3",
            "number_density_m3",
            "temperature_k",
            "temperature_ev",
        ):
            if key in species_table:
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}] uses "
                    f'velocity_distribution="grid" and cannot define {key}.'
                )
        return

    _validate_velocity_grid_forbidden(
        species_table,
        index=index,
        source_mode=source_mode,
    )
    density_keys = [
        key
        for key in ("number_density_cm3", "number_density_m3")
        if key in species_table
    ]
    if len(density_keys) != 1:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] reservoir injection "
            "requires exactly one of number_density_cm3 and number_density_m3."
        )
    density = species_table[density_keys[0]]
    if (
        not isinstance(density, (int, float))
        or isinstance(density, bool)
        or not math.isfinite(float(density))
        or float(density) <= 0.0
    ):
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}].{density_keys[0]} "
            "must be finite and > 0."
        )
    for key in ("temperature_k", "temperature_ev"):
        if key not in species_table:
            continue
        temperature = species_table[key]
        if (
            not isinstance(temperature, (int, float))
            or isinstance(temperature, bool)
            or not math.isfinite(float(temperature))
            or float(temperature) < 0.0
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].{key} "
                "must be finite and >= 0."
            )


def _validate_plane_source_geometry(
    species_table: Mapping[str, Any],
    *,
    index: int,
    use_box: bool,
    box_min: Sequence[float] | None,
    box_max: Sequence[float] | None,
) -> None:
    if not use_box or box_min is None or box_max is None:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            'source_mode="plane_source" and requires [domain].'
        )
    pos_low = _maybe_vec3(
        species_table.get("pos_low"),
        name=f"particles.species[{index}].pos_low",
    )
    pos_high = _maybe_vec3(
        species_table.get("pos_high"),
        name=f"particles.species[{index}].pos_high",
    )
    source_normal = _maybe_vec3(
        species_table.get("source_normal"),
        name=f"particles.species[{index}].source_normal",
    )
    if pos_low is None or pos_high is None or source_normal is None:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] plane_source "
            "requires pos_low, pos_high, and source_normal."
        )

    zero_axes = [
        axis
        for axis in range(3)
        if math.isclose(pos_low[axis], pos_high[axis], rel_tol=0.0, abs_tol=1.0e-12)
    ]
    if len(zero_axes) != 1:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] plane_source "
            "pos_low/pos_high must define an axis-aligned zero-thickness rectangle."
        )
    normal_axis = zero_axes[0]
    for axis in range(3):
        low_bound = float(box_min[axis])
        high_bound = float(box_max[axis])
        if axis == normal_axis:
            if not (low_bound < pos_low[axis] < high_bound):
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}] plane_source "
                    "must lie strictly inside the box along its normal axis."
                )
        elif not (low_bound <= pos_low[axis] < pos_high[axis] <= high_bound):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] plane_source "
                "must have positive in-box extent along both tangential axes."
            )

    for axis, component in enumerate(source_normal):
        invalid_component = (
            abs(component) <= 1.0e-12
            if axis == normal_axis
            else abs(component) > 1.0e-12
        )
        if invalid_component:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].source_normal "
                "must be a non-zero axis-aligned vector along the plane normal axis."
            )


def _validate_particle_species(
    species: list[Mapping[str, Any]], *, use_box: bool, resolved_batch_duration: float,
    box_min: list[float] | None, box_max: list[float] | None,
    periodic_axes: set[str], global_particle_boundary: Mapping[str, Any],
    automatic_current_species: set[str],
) -> tuple[bool, int]:
    """Validate species and report which source families contribute particles."""
    uses_face_sources = False
    total_npcls_per_step = 0
    species_keys: set[str] = set()
    for index, species_table in enumerate(species, start=1):
        species_key = species_table.get("species_key", f"species_{index}")
        if species_key in species_keys:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].species_key "
                f"duplicates {species_key!r}; species_key values must be unique."
            )
        species_keys.add(species_key)
        if not species_table.get("enabled", True):
            continue
        if species_table.get("q_particle", -1.602176634e-19) == 0.0:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].q_particle "
                "must be non-zero for an enabled species."
            )
        if not use_box and any(
            action != "inherit"
            for action in species_table.get("boundary", {}).values()
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].boundary "
                "requires a finite [domain]."
            )
        source_mode = species_table.get("source_mode", "volume_seed")
        effective_particle_boundary = _validate_species_particle_boundary(
            species_table,
            index=index,
            periodic_axes=periodic_axes,
            global_boundary=global_particle_boundary,
        )
        boundary_inflow_faces = _validate_species_boundary_inflow(
            species_table,
            index=index,
            periodic_axes=periodic_axes,
            effective_boundary=effective_particle_boundary,
        )
        has_reservoir_injection = source_mode in _RESERVOIR_SOURCE_MODES or bool(
            boundary_inflow_faces
        )
        if has_reservoir_injection:
            uses_face_sources = True
            _validate_reservoir_injection_common(
                species_table,
                index=index,
                use_box=use_box,
                batch_duration=resolved_batch_duration,
            )
            _validate_reservoir_physics(
                species_table,
                index=index,
                source_mode=source_mode,
            )
        surface_charge_closure = species_table.get("surface_charge_closure", "explicit")
        fixed_absorbed_present = "target_absorbed_current_a" in species_table
        fixed_emission_present = "target_emission_current_a" in species_table
        if surface_charge_closure == "explicit" and (
            fixed_absorbed_present or fixed_emission_present
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] target "
                'surface currents require surface_charge_closure="fixed_current".'
            )
        if surface_charge_closure == "fixed_current":
            if (
                species_table.get("enabled", True) is not False
                and (
                    not math.isfinite(resolved_batch_duration)
                    or resolved_batch_duration <= 0.0
                )
            ):
                raise ConfigValidationError(
                    "BEACH constraint error: sim.batch_duration must be > 0 "
                    "for fixed_current."
                )
            species_key = str(species_table.get("species_key", f"species_{index}"))
            if not (
                fixed_absorbed_present
                or fixed_emission_present
                or species_key in automatic_current_species
            ):
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}] "
                    "fixed_current requires at least one target current."
                )
            q_particle = species_table.get("q_particle", -1.602176634e-19)
            if not isinstance(q_particle, (int, float)) or isinstance(q_particle, bool):
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}]."
                    "q_particle must be numeric for fixed_current."
                )
            if fixed_absorbed_present:
                target = species_table["target_absorbed_current_a"]
                if (
                    not isinstance(target, (int, float))
                    or isinstance(target, bool)
                    or not math.isfinite(float(target))
                    or float(target) * float(q_particle) < 0.0
                ):
                    raise ConfigValidationError(
                        f"BEACH constraint error: particles.species[{index}]."
                        "target_absorbed_current_a must be finite and have the "
                        "same sign as q_particle."
                    )
            if fixed_emission_present:
                target = species_table["target_emission_current_a"]
                if (
                    source_mode != "photo_raycast"
                    or species_table.get("deposit_opposite_charge_on_emit") is not True
                    or not isinstance(target, (int, float))
                    or isinstance(target, bool)
                    or not math.isfinite(float(target))
                    or float(target) * float(q_particle) > 0.0
                ):
                    raise ConfigValidationError(
                        f"BEACH constraint error: particles.species[{index}]."
                        "target_emission_current_a requires photo_raycast, an "
                        "opposite-charge emission deposit, and current sign "
                        "opposite to q_particle."
                    )
        if surface_charge_closure == "neutral_return":
            inject_face = species_table.get("inject_face")
            q_particle = species_table.get("q_particle", -1.602176634e-19)
            if (
                source_mode != "photo_raycast"
                or not isinstance(q_particle, (int, float))
                or isinstance(q_particle, bool)
                or float(q_particle) >= 0.0
                or species_table.get("deposit_opposite_charge_on_emit") is not True
                or not isinstance(inject_face, str)
                or effective_particle_boundary.get(inject_face)
                not in {"reflect", "redistributed_reflect"}
            ):
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}] "
                    "neutral_return requires a negative photo_raycast species, "
                    "deposit_opposite_charge_on_emit=true, and a reflecting "
                    "action on inject_face."
                )
        if "photo_escape_model" in species_table:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].photo_escape_model "
                "was removed; track emitted photoelectrons and use the neutral-return "
                "surface charge closure instead."
            )
        if "temperature_k" in species_table and "temperature_ev" in species_table:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] cannot define both "
                "temperature_k and temperature_ev."
            )
        if source_mode == "volume_seed":
            if not has_reservoir_injection:
                _validate_velocity_grid_forbidden(
                    species_table,
                    index=index,
                    source_mode=source_mode,
                )
            npcls_per_step = species_table.get("npcls_per_step", 0)
            total_npcls_per_step += npcls_per_step
            if (
                not has_reservoir_injection
                and "target_macro_particles_per_batch" in species_table
            ):
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}] uses "
                    'source_mode="volume_seed" and cannot define '
                    "target_macro_particles_per_batch."
                )
            continue

        if source_mode == "reservoir_face":
            _validate_face_source_common(
                species_table,
                index=index,
                source_mode=source_mode,
                use_box=use_box,
                batch_duration=resolved_batch_duration,
                box_min=box_min,
                box_max=box_max,
            )
            continue

        if source_mode == "plane_source":
            _validate_plane_source_geometry(
                species_table,
                index=index,
                use_box=use_box,
                box_min=box_min,
                box_max=box_max,
            )
            continue

        if source_mode == "photo_raycast":
            uses_face_sources = True
            _validate_velocity_grid_forbidden(
                species_table,
                index=index,
                source_mode=source_mode,
            )
            _validate_face_source_common(
                species_table,
                index=index,
                source_mode=source_mode,
                use_box=use_box,
                batch_duration=resolved_batch_duration,
                box_min=box_min,
                box_max=box_max,
            )
            current_density = species_table.get("emit_current_density_a_m2", 0.0)
            rays_per_batch = species_table.get("rays_per_batch", 0)
            if (
                not isinstance(current_density, (int, float))
                or float(current_density) <= 0.0
            ):
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}] uses "
                    'source_mode="photo_raycast" and requires emit_current_density_a_m2 > 0.'
                )
            if not isinstance(rays_per_batch, int) or rays_per_batch <= 0:
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}] uses "
                    'source_mode="photo_raycast" and requires rays_per_batch > 0.'
                )
            forbidden = ["npcls_per_step"]
            if not has_reservoir_injection:
                forbidden.extend(
                    (
                        "number_density_cm3",
                        "number_density_m3",
                        "w_particle",
                        "target_macro_particles_per_batch",
                    )
                )
            for key in forbidden:
                if key in species_table:
                    raise ConfigValidationError(
                        f"BEACH constraint error: particles.species[{index}] uses "
                        f'source_mode="photo_raycast" and cannot define {key}.'
                    )
            continue

        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] has unsupported "
            f"source_mode={source_mode!r}."
        )
    return uses_face_sources, total_npcls_per_step

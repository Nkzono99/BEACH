"""Coordinate validation of runtime tables and cross-model constraints."""

from __future__ import annotations

import copy
import math
from collections.abc import Mapping
from typing import Any

from ._field_validation import (
    _validate_field_config,
    _validate_runtime_boundary_tables,
    _validate_runtime_external_e_field,
)
from ._mesh_validation import (
    _validate_runtime_mesh,
)
from ._particle_validation import (
    _validate_particle_species,
)
from ._shared import (
    _REMOVED_SIM_KEYS,
    _REQUIRED_RUNTIME_TABLES,
    _RESERVOIR_SOURCE_MODES,
    ConfigValidationError,
    _maybe_vec3,
    _optional_runtime_table,
    _require_table,
    _validate_fragment_structure,
)
from ._surface_validation import (
    _validate_surface_current_model,
)


def validate_runtime_config(config: Mapping[str, Any]) -> None:
    """Validate the merged final config against known BEACH constraints."""

    final_config = copy.deepcopy(dict(config))
    _validate_fragment_structure(
        final_config,
        context="runtime config",
        allow_meta_keys=False,
    )
    for key in _REQUIRED_RUNTIME_TABLES:
        if key not in final_config:
            raise ConfigValidationError(
                f"BEACH constraint error: runtime config is missing top-level [{key}] table."
            )

    sim = _require_table(final_config, "sim", context="runtime config")
    particles = _require_table(final_config, "particles", context="runtime config")
    mesh = _require_table(final_config, "mesh", context="runtime config")
    _require_table(final_config, "output", context="runtime config")
    domain = _optional_runtime_table(final_config, "domain")
    field_boundary = _optional_runtime_table(final_config, "field_boundary")
    particle_boundary = _optional_runtime_table(final_config, "particle_boundary")
    reservoir = _optional_runtime_table(final_config, "reservoir")
    surface_current_model = _optional_runtime_table(
        final_config, "surface_current_model"
    )
    periodic2_config = final_config.get("periodic2", {})
    removed_sim_keys = sorted(set(sim) & _REMOVED_SIM_KEYS)
    if removed_sim_keys:
        raise ConfigValidationError(
            "BEACH constraint error: removed sim key(s): "
            + ", ".join(removed_sim_keys)
            + "."
        )
    _validate_runtime_external_e_field(sim)
    _validate_runtime_mesh(mesh)
    _validate_runtime_boundary_tables(
        domain=domain,
        field_boundary=field_boundary,
        particle_boundary=particle_boundary,
        reservoir=reservoir,
    )

    species = particles.get("species")
    if (
        not isinstance(species, list)
        or len(species) == 0
        or not all(isinstance(item, Mapping) for item in species)
    ):
        raise ConfigValidationError(
            "BEACH constraint error: particles.species must be a non-empty array of tables."
        )
    automatic_current_species = set()
    if (
        surface_current_model is not None
        and surface_current_model.get("model") == "zhao_stationary"
    ):
        automatic_current_species = {
            str(surface_current_model.get(key, ""))
            for key in (
                "electron_species",
                "ion_species",
                "photoelectron_species",
            )
        }

    resolved_batch_duration = _resolve_batch_duration(sim)
    use_box = domain is not None
    adaptive_nonzero_mode_limit, field_bc_mode = _validate_field_config(
        sim=sim, domain=domain, field_boundary=field_boundary, reservoir=reservoir,
        periodic2_config=periodic2_config, mesh=mesh,
        resolved_batch_duration=resolved_batch_duration, use_box=use_box,
    )

    box_min = (
        _maybe_vec3(domain.get("box_min"), name="domain.box_min")
        if domain is not None
        else None
    )
    box_max = (
        _maybe_vec3(domain.get("box_max"), name="domain.box_max")
        if domain is not None
        else None
    )
    periodic_axes = (
        set(domain.get("periodic_axes", [])) if domain is not None else set()
    )
    global_particle_boundary = particle_boundary or {}

    multiple_event_retry_backend = sim.get(
        "multiple_box_events_retry_backend", "none"
    )
    if multiple_event_retry_backend not in {"none", "upper_panel_fourier"}:
        raise ConfigValidationError(
            "BEACH constraint error: sim.multiple_box_events_retry_backend "
            'must be "none" or "upper_panel_fourier".'
        )
    retry_nonzero_backend = (
        periodic2_config.get("nonzero_mode_backend")
        if isinstance(periodic2_config, Mapping)
        else None
    )
    if (
        retry_nonzero_backend is None
        and sim.get("field_periodic_far_correction") == "cached_kneq0"
    ):
        retry_nonzero_backend = "cached_kneq0"
    if multiple_event_retry_backend == "upper_panel_fourier" and (
        field_bc_mode != "periodic2"
        or retry_nonzero_backend != "cached_kneq0"
    ):
        raise ConfigValidationError(
            "BEACH constraint error: upper_panel_fourier retry requires "
            'field_boundary.mode="periodic2" and '
            'periodic2.nonzero_mode_backend="cached_kneq0".'
        )

    multiple_event_policy = sim.get("multiple_box_events_policy", "abort")
    if multiple_event_policy == "soft_discard":
        count_grace = sim.get("multiple_box_events_soft_discard_count_grace", 1000)
        fraction_limit = sim.get(
            "multiple_box_events_soft_discard_fraction_limit", 1.0e-6
        )
        charge_limit = sim.get(
            "multiple_box_events_soft_discard_abs_charge_limit", 1.0e-12
        )
        if (
            not isinstance(count_grace, int)
            or isinstance(count_grace, bool)
            or count_grace < 0
        ):
            raise ConfigValidationError(
                "BEACH constraint error: soft_discard count grace must be an integer >= 0."
            )
        if (
            not isinstance(fraction_limit, (int, float))
            or isinstance(fraction_limit, bool)
            or not math.isfinite(float(fraction_limit))
            or float(fraction_limit) <= 0.0
            or float(fraction_limit) > 1.0
        ):
            raise ConfigValidationError(
                "BEACH constraint error: soft_discard fraction limit must be finite and in (0, 1]."
            )
        if (
            not isinstance(charge_limit, (int, float))
            or isinstance(charge_limit, bool)
            or not math.isfinite(float(charge_limit))
            or float(charge_limit) <= 0.0
        ):
            raise ConfigValidationError(
                "BEACH constraint error: soft_discard absolute charge limit must be finite and > 0."
            )

    if any(
        bool(item.get("enabled", True))
        and item.get("source_mode", "volume_seed") == "photo_raycast"
        for item in species
    ):
        max_bounce = sim.get("raycast_max_bounce", 16)
        if (
            not isinstance(max_bounce, int)
            or isinstance(max_bounce, bool)
            or max_bounce < 1
        ):
            raise ConfigValidationError(
                "BEACH constraint error: sim.raycast_max_bounce must be an integer >= 1 "
                "when photo_raycast is enabled."
            )

    uses_face_sources, has_volume_seed, total_npcls_per_step = _validate_particle_species(
        species, use_box=use_box, resolved_batch_duration=resolved_batch_duration,
        box_min=box_min, box_max=box_max, periodic_axes=periodic_axes,
        global_particle_boundary=global_particle_boundary,
        automatic_current_species=automatic_current_species,
    )

    _validate_surface_current_model(
        surface_current_model,
        species=species,
        sim=sim,
        domain=domain,
        field_boundary=field_boundary,
        particle_boundary=particle_boundary,
        reservoir=reservoir,
        periodic2_config=periodic2_config,
    )

    if has_volume_seed and not uses_face_sources and total_npcls_per_step < 1:
        raise ConfigValidationError(
            "BEACH constraint error: volume_seed species require total npcls_per_step >= 1."
        )
    if adaptive_nonzero_mode_limit > 0.0 and any(
        item.get("enabled", True) is True
        and item.get("source_mode", "volume_seed") == "volume_seed"
        and isinstance(item.get("npcls_per_step", 0), int)
        and item.get("npcls_per_step", 0) > 0
        for item in species
    ):
        raise ConfigValidationError(
            "BEACH constraint error: "
            "periodic2.max_nonzero_mode_potential_step requires time-scaled "
            "reservoir_face/photo_raycast sources."
        )
    if adaptive_nonzero_mode_limit > 0.0 and any(
        item.get("enabled", True) is True
        and (
            item.get("source_mode", "volume_seed") in _RESERVOIR_SOURCE_MODES
            or bool(item.get("boundary_inflow"))
        )
        and "target_macro_particles_per_batch" not in item
        for item in species
    ):
        raise ConfigValidationError(
            "BEACH constraint error: adaptive reservoir injection requires "
            "target_macro_particles_per_batch instead of fixed w_particle."
        )


def _resolve_batch_duration(sim: Mapping[str, Any]) -> float:
    dt = float(sim.get("dt", 1.0e-9))
    has_batch_duration = "batch_duration" in sim
    has_batch_duration_step = "batch_duration_step" in sim
    if has_batch_duration and has_batch_duration_step:
        raise ConfigValidationError(
            "BEACH constraint error: sim.batch_duration and sim.batch_duration_step "
            "cannot be specified together."
        )
    if has_batch_duration_step:
        return dt * float(sim["batch_duration_step"])
    return float(sim.get("batch_duration", 0.0))

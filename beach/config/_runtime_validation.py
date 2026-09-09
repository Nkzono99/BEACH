"""Coordinate validation of runtime tables and cross-model constraints."""

from __future__ import annotations

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
    _RESERVOIR_SOURCE_MODES,
    ConfigValidationError,
    _maybe_vec3,
    _optional_runtime_table,
)
from ._surface_validation import (
    _validate_surface_current_model,
)


def _validate_runtime_semantics(config: Mapping[str, Any]) -> None:
    """Check physical combinations after the shared schema/value validation."""

    sim = config.get("sim", {})
    particles = config["particles"]
    mesh = config.get("mesh", {})
    domain = _optional_runtime_table(config, "domain")
    field_boundary = _optional_runtime_table(config, "field_boundary")
    particle_boundary = _optional_runtime_table(config, "particle_boundary")
    reservoir = _optional_runtime_table(config, "reservoir")
    surface_current_model = _optional_runtime_table(
        config, "surface_current_model"
    )
    periodic2_config = config.get("periodic2")
    _validate_runtime_external_e_field(sim)
    _validate_runtime_mesh(mesh)
    _validate_runtime_boundary_tables(
        domain=domain,
        field_boundary=field_boundary,
        particle_boundary=particle_boundary,
        reservoir=reservoir,
    )

    species = particles.get("species")
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

    uses_face_sources, total_npcls_per_step = _validate_particle_species(
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

    if not uses_face_sources and total_npcls_per_step < 1:
        raise ConfigValidationError(
            "BEACH constraint error: at least one enabled particle source is required; "
            "volume_seed species require total npcls_per_step >= 1."
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
        duration = dt * float(sim["batch_duration_step"])
        if not math.isfinite(duration) or duration <= 0.0:
            raise ConfigValidationError(
                "BEACH constraint error: sim.dt * sim.batch_duration_step must "
                "produce a finite positive sim.batch_duration."
            )
        return duration
    return float(sim.get("batch_duration", 0.0))

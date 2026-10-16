"""Map grouped public input to internal reader paths and render runtime dictionaries.

The Fortran routing table and editor schema are generated from these paths.
"""
from __future__ import annotations

import copy
from collections.abc import Mapping
from typing import Any

from ._shared import ConfigValidationError

PATHS = {'sim.dt': 'particles.tracking.dt_s',
 'sim.max_step': 'particles.tracking.max_steps_per_particle',
 'sim.rng_seed': 'run.rng_seed',
 'sim.batch_count': 'run.batch_count',
 'sim.batch_duration': 'run.batch.duration_s',
 'sim.batch_duration_step': 'run.batch.duration_steps',
 'sim.tol_rel': 'output.diagnostics.charge_rel_change_threshold',
 'sim.q_floor': 'output.diagnostics.charge_floor_c',
 'sim.field_solver': 'fields.solver.method',
 'sim.field_normalization': 'fields.solver.normalization',
 'sim.field_length_scale': 'fields.solver.length_scale_m',
 'sim.tree_theta': 'fields.solver.tree.theta',
 'sim.tree_leaf_max': 'fields.solver.tree.leaf_max',
 'sim.tree_min_nelem': 'fields.solver.auto_fmm_min_elements',
 'sim.field_periodic_image_layers': 'fields.periodic.image_layers',
 'sim.field_periodic_far_correction': 'fields.periodic.backend',
 'sim.field_periodic_ewald_alpha': 'fields.periodic.ewald_alpha',
 'sim.field_periodic_ewald_layers': 'fields.periodic.ewald_layers',
 'sim.field_periodic_cache_dir': 'fields.periodic.cache_dir',
 'sim.field_periodic_generation_tolerance': 'fields.periodic.generation_tolerance',
 'sim.e0': 'fields.external.electric_v_m',
 'sim.e0_abs': 'fields.external.electric_magnitude_v_m',
 'sim.e0_phi_xy_deg': 'fields.external.electric_azimuth_deg',
 'sim.e0_phi_z_deg': 'fields.external.electric_elevation_deg',
 'sim.b0': 'fields.external.magnetic_t',
 'sim.multiple_box_events_policy': 'particles.tracking.events.policy',
 'sim.multiple_box_events_retry_backend': 'particles.tracking.events.retry_backend',
 'sim.multiple_box_events_soft_discard_count_grace': 'particles.tracking.events.discard_count_grace',
 'sim.multiple_box_events_soft_discard_fraction_limit': 'particles.tracking.events.discard_fraction_limit',
 'sim.multiple_box_events_soft_discard_abs_charge_limit': 'particles.tracking.events.discard_charge_warning_c',
 'sim.raycast_max_bounce': 'particles.raycast.max_bounces',
 'particles.secondary_emission_max_generations': 'particles.secondary_emission.max_generations',
 'periodic2.upper_fourier_switch_height': 'fields.periodic.upper_fourier.switch_height',
 'periodic2.upper_fourier_mode_layers': 'fields.periodic.upper_fourier.mode_layers',
 'periodic2.lower_boundary_model': 'fields.periodic.lower_boundary_model',
 'periodic2.reference_mode_layers': 'fields.periodic.reference.mode_layers',
 'periodic2.panel_quadrature_order': 'fields.periodic.reference.panel_quadrature_order',
 'periodic2.max_nonzero_mode_potential_step': 'run.batch.adaptive.max_nonzero_mode_potential_step_v',
 'surface_current_model.model': 'sheath.closure',
 'surface_current_model.zhao_branch': 'sheath.zhao.branch',
 'surface_current_model.electron_species': 'sheath.species.electron',
 'surface_current_model.ion_species': 'sheath.species.ion',
 'surface_current_model.photoelectron_species': 'sheath.species.photoelectron',
 'surface_current_model.secondary_electron_species': 'sheath.species.secondary_electron',
 'surface_current_model.solar_elevation_deg': 'sheath.photoelectrons.solar_elevation_deg',
 'surface_current_model.photoelectron_ref_density_m3': 'sheath.photoelectrons.ref_density_m3',
 'surface_current_model.photoelectron_source_scale': 'sheath.photoelectrons.source_scale',
 'surface_current_model.reference_area_m2': 'sheath.reference_area_m2',
 'surface_current_model.outflow_refresh_batches': 'sheath.coupling.outflow_refresh_batches',
 'surface_current_model.zhao_upstream_band_tolerance': 'sheath.zhao.upstream_band_tolerance',
 'surface_current_model.zhao_upstream_field_floor_v_m': 'sheath.zhao.upstream_field_floor_v_m',
 'surface_current_model.zhao_ion_pressure_factor': 'sheath.zhao.ion_pressure_factor',
 'output.write_files': 'output.enabled',
 'output.dir': 'output.dir',
 'output.write_mesh_potential': 'output.final.mesh_potential',
 'output.write_potential_history': 'output.history.potential',
 'output.history_stride': 'output.history.stride_batches',
 'output.checkpoint_stride': 'output.checkpoint.stride_batches',
 'output.restart_from': 'run.restart.from',
 'field_boundary.mode': 'fields.boundary',
 'particle_boundary.x_low': 'particles.boundary.x_low',
 'particle_boundary.x_high': 'particles.boundary.x_high',
 'particle_boundary.y_low': 'particles.boundary.y_low',
 'particle_boundary.y_high': 'particles.boundary.y_high',
 'particle_boundary.z_low': 'particles.boundary.z_low',
 'particle_boundary.z_high': 'particles.boundary.z_high',
 'reservoir.inflow_model': 'particles.reservoir.inflow_model',
 'reservoir.phi_infty': 'particles.reservoir.phi_infty_v',
 'reservoir.face_potential_grid_n': 'particles.reservoir.face_potential_grid_n'}

SPECIES_PATHS = {'species_key': 'species_key',
 'enabled': 'enabled',
 'q_particle': 'charge_c',
 'm_particle': 'mass_kg',
 'number_density_cm3': 'distribution.number_density_cm3',
 'number_density_m3': 'distribution.number_density_m3',
 'temperature_ev': 'distribution.temperature_ev',
 'temperature_k': 'distribution.temperature_k',
 'drift_velocity': 'distribution.drift_velocity_m_s',
 'source_mode': 'source.mode',
 'pos_low': 'source.pos_low',
 'pos_high': 'source.pos_high',
 'inject_face': 'source.inject_face',
 'ray_direction': 'source.ray_direction',
 'inject_region_mode': 'source.region_mode',
 'uv_low': 'source.uv_low',
 'uv_high': 'source.uv_high',
 'emit_current_density_a_m2': 'source.emit_current_density_a_m2',
 'deposit_opposite_charge_on_emit': 'source.deposit_opposite_charge_on_emit',
 'normal_drift_speed': 'source.normal_drift_speed_m_s',
 'npcls_per_step': 'sampling.volume_macro_particles_per_batch',
 'w_particle': 'sampling.weight',
 'target_macro_particles_per_batch': 'sampling.target_macro_particles_per_batch',
 'rays_per_batch': 'sampling.rays_per_batch',
 'surface_charge_closure': 'charging.closure',
 'target_absorbed_current_a': 'charging.target_absorbed_current_a',
 'target_emission_current_a': 'charging.target_emission_current_a',
 'boundary': 'boundary',
 'boundary_inflow': 'inflow'}

SECONDARY_EMISSION_RULES = {
    'enabled': {'type': 'boolean', 'default': False},
    'model': {'type': 'string', 'enum': ['whipple', 'constant'], 'default': 'whipple'},
    'yield_max': {'type': 'number', 'default': 1.5},
    'energy_max_ev': {'type': 'number', 'default': 500.0},
    'const_yield': {'type': 'number', 'default': 0.0},
    'target_species': {'type': 'string', 'maxLength': 64},
}
SPECIES_PATHS.update({f'secondary_emission_{key}': f'surface_interaction.secondary_emission.{key}'
                      for key in SECONDARY_EMISSION_RULES})

TOP_LEVEL_ORDER = ("run", "domain", "mesh", "particles", "fields", "sheath", "output")
PRESERVED = {"domain": "domain", "mesh.groups": "mesh.groups", "mesh.parts": "mesh.parts"}
CHOICES = {
    "sheath.closure": {"none": "none", "zero_current": "zhao_stationary"},
}


READER_DEFINITIONS = {
    "sim": "sim", "domain": "domain", "mesh": "mesh", "particles": "particles",
    "field_boundary": "fieldBoundary", "particle_boundary": "particleBoundary",
    "reservoir": "reservoir", "surface_current_model": "surfaceCurrentModel",
    "periodic2": "periodic2", "output": "output",
}


def reader_schema(schema: Mapping[str, Any]) -> dict[str, Any]:
    """Build internal reader rules for schema generation and key-contract checks."""
    return {"type": "object", "additionalProperties": False,
            "properties": {key: {"$ref": f"#/$defs/{name}"}
                           for key, name in READER_DEFINITIONS.items()},
            "$defs": schema["$defs"]}


def is_grouped(config: Mapping[str, Any]) -> bool:
    """Recognize a grouped document, including partial boundary-only documents."""
    def has(table, keys):
        return isinstance(table, Mapping) and bool({str(k).lower() for k in table} & keys)
    if has(config, {"sim", "field_boundary", "particle_boundary", "surface_current_model", "periodic2"}):
        return False
    if has(config, {"run", "fields", "sheath"}):
        return True
    lowered = {str(k).lower(): v for k, v in config.items()}
    particles = lowered.get("particles", {})
    if has(particles, {"boundary", "reservoir", "tracking", "raycast", "secondary_emission"}):
        return True
    if isinstance(particles, Mapping):
        species = next((v for k,v in particles.items() if str(k).lower()=="species"), [])
        if isinstance(species, list) and any(has(item, {"charge_c", "mass_kg", "distribution",
                "source", "sampling", "inflow", "charging", "surface_interaction"}) for item in species):
            return True
    return "mesh" in lowered or has(lowered.get("output"),
            {"enabled", "history", "final", "checkpoint", "diagnostics"})


def _put(document, path, value):
    parts = path.split(".")
    for part in parts[:-1]:
        document = document.setdefault(part, {})
    if parts[-1] in document:
        raise ConfigValidationError(f"duplicate setting: {path}")
    document[parts[-1]] = copy.deepcopy(value)


def _leaves(document, prefix="", preserved=()):
    for key, value in document.items():
        path = f"{prefix}.{key}" if prefix else key
        if isinstance(value, Mapping) and path not in preserved:
            yield from _leaves(value, path, preserved)
        else:
            yield path, value


def to_runtime_layout(config: Mapping[str, Any]) -> dict[str, Any]:
    """Normalize grouped input without changing physical values or units."""
    if not is_grouped(config):
        raise ConfigValidationError(
            "flat beach.toml layout is not supported; use the grouped layout "
            "(run/domain/mesh/particles/fields/sheath/output)"
        )
    reverse = {new: old for old, new in PATHS.items()}
    reverse.update(PRESERVED)
    result = {}
    preserved = {*PRESERVED, "particles.species", "run.restart"}
    for path, value in _leaves(config, preserved=preserved):
        if path == "particles.species":
            species = []
            for item in value:
                flat = {}
                inverse = {new: old for old, new in SPECIES_PATHS.items()}
                for key, val in _leaves(item, preserved={"boundary", "inflow"}):
                    if key not in inverse:
                        raise ConfigValidationError(f"unknown or mixed-layout key: particles.species.{key}")
                    flat[inverse[key]] = copy.deepcopy(val)
                if "source" in item and "mode" not in item["source"]:
                    raise ConfigValidationError("particles.species.source requires mode")
                if "source" not in item and item.get("sampling", {}).get("volume_macro_particles_per_batch", 0) != 0:
                    raise ConfigValidationError("volume samples require particles.species.source")
                species.append(flat)
            _put(result, "particles.species", species)
        elif path == "run.restart":
            _put(result, "output.resume", True)
            for key, val in value.items():
                if key != "from":
                    raise ConfigValidationError(f"unknown run.restart key: {key}")
                _put(result, "output.restart_from", val)
        elif path == "fields.periodic.backend":
            if value not in {"cached_kneq0", "panel_spectral_reference", "finite_images"}:
                raise ConfigValidationError("invalid fields.periodic.backend")
            _put(result, "sim.field_periodic_far_correction",
                 "cached_kneq0" if value == "cached_kneq0" else "none")
            if value != "finite_images":
                _put(result, "periodic2.nonzero_mode_backend", value)
                _put(result, "periodic2.zero_mode_policy", "exclude_k0")
        elif path in reverse:
            if path in CHOICES:
                try:
                    value = CHOICES[path][value]
                except (KeyError, TypeError) as exc:
                    raise ConfigValidationError(f"invalid {path}") from exc
            _put(result, reverse[path], value)
        else:
            raise ConfigValidationError(f"unknown or mixed-layout key: {path}")
    periodic = config.get("fields", {}).get("periodic", {})
    if periodic.get("backend", "finite_images") == "finite_images" and (
            "lower_boundary_model" in periodic or "reference" in periodic or "upper_fourier" in periodic):
        raise ConfigValidationError("split periodic settings require an explicit split backend")
    # An empty source table has semantics, even though it has no leaves.
    for item in config.get("particles", {}).get("species", []):
        if "source" in item and "mode" not in item["source"]:
            raise ConfigValidationError("particles.species.source requires mode")
    return result


def to_grouped_layout(config: Mapping[str, Any]) -> dict[str, Any]:
    """Render an internal runtime dictionary using the grouped public input layout."""
    if is_grouped(config):
        return copy.deepcopy(dict(config))
    result = {}
    paths = {**PATHS, **{old: new for new, old in PRESERVED.items()}}
    for path, value in _leaves(config, preserved={*PRESERVED.values(), "particles.species", "periodic2"}):
        if path == "particles.species":
            species = []
            for item in value:
                grouped = {}
                for key, val in item.items():
                    if key not in SPECIES_PATHS:
                        raise ConfigValidationError(f"unknown species key: {key}")
                    # Default volume_seed/count=0 means no standalone source.
                    if item.get("source_mode", "volume_seed") == "volume_seed" and item.get("npcls_per_step", 0) == 0:
                        if SPECIES_PATHS[key].startswith("source.") or key == "npcls_per_step":
                            continue
                    _put(grouped, SPECIES_PATHS[key], val)
                species.append(grouped)
            _put(result, "particles.species", species)
        elif path == "periodic2":
            for key, val in value.items():
                if key in {"nonzero_mode_backend", "zero_mode_policy"}:
                    continue
                old = "periodic2." + key
                if old not in paths:
                    raise ConfigValidationError(f"unknown key: {old}")
                _put(result, paths[old], val)
        elif path == "output.resume":
            if value:
                result.setdefault("run", {}).setdefault("restart", {})
        elif path == "output.restart_from" and not config.get("output", {}).get("resume", False):
            # An inactive internal restart path must not turn a fresh run into a resume.
            if value:
                raise ConfigValidationError("set output.resume=true before rendering restart_from")
        elif path == "sim.field_periodic_far_correction":
            continue  # Select the effective typed backend below.
        elif path in paths:
            new = paths[path]
            if new in CHOICES:
                value = {old: new for new, old in CHOICES[new].items()}[value]
            _put(result, new, value)
        else:
            raise ConfigValidationError(f"unknown runtime key: {path}")
    periodic = config.get("periodic2")
    far = config.get("sim", {}).get("field_periodic_far_correction", "none")
    if isinstance(periodic, Mapping) or far == "cached_kneq0":
        backend = periodic.get("nonzero_mode_backend", "panel_spectral_reference") if periodic is not None else far
        _put(result, "fields.periodic.backend", backend)
        # Materialize the effective lower boundary of internal cached settings.
        lower = periodic.get("lower_boundary_model", "e_bottom_zero") if periodic is not None else "e_bottom_zero"
        if "lower_boundary_model" not in result["fields"]["periodic"]:
            _put(result, "fields.periodic.lower_boundary_model", lower)
    elif "field_periodic_far_correction" in config.get("sim", {}):
        _put(result, "fields.periodic.backend", "finite_images")
    # Internal source mode defaults also apply to volume samples with no explicit mode.
    for item in result.get("particles", {}).get("species", []):
        if "source" in item or item.get("sampling", {}).get("volume_macro_particles_per_batch", 0) > 0:
            item.setdefault("source", {}).setdefault("mode", "volume_seed")
    return {key: result[key] for key in TOP_LEVEL_ORDER if key in result}

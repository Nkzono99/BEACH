"""Validate stationary-sheath and matching-plane surface-current models."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

from ._particle_validation import (
    _validate_species_particle_boundary,
)
from ._shared import (
    ConfigValidationError,
    _maybe_vec3,
)


def _validate_surface_current_model(
    model_config: Mapping[str, Any] | None,
    *,
    species: list[Mapping[str, Any]],
    sim: Mapping[str, Any],
    domain: Mapping[str, Any] | None,
    field_boundary: Mapping[str, Any] | None,
    particle_boundary: Mapping[str, Any] | None,
    reservoir: Mapping[str, Any] | None,
    periodic2_config: object,
) -> None:
    if model_config is None:
        return
    model = model_config.get("model")
    if model is None or model == "none":
        if set(model_config) - {"model"}:
            raise ConfigValidationError(
                'BEACH constraint error: surface_current_model.model="none" cannot '
                "use Zhao or matching-plane settings."
            )
        return
    matching_keys = {
        "response_backend",
        "response_table_path",
        "zhao_root_selection",
        "implicit_zero_mode",
        "coupling_rtol",
        "coupling_atol",
        "coupling_max_iterations",
        "coupling_relaxation",
    }
    if model == "matching_plane_quasistatic":
        _validate_matching_plane_model(
            model_config,
            species=species,
            sim=sim,
            domain=domain,
            field_boundary=field_boundary,
            particle_boundary=particle_boundary,
            reservoir=reservoir,
            periodic2_config=periodic2_config,
        )
        return
    if model != "zhao_stationary":
        raise ConfigValidationError(
            'BEACH constraint error: surface_current_model.model must be "none", '
            '"zhao_stationary", or "matching_plane_quasistatic".'
        )
    if matching_keys.intersection(model_config):
        raise ConfigValidationError(
            "BEACH constraint error: Zhao surface-current model cannot use "
            "matching-plane-specific settings."
        )
    zhao_branch = model_config.get("zhao_branch", "auto")
    if zhao_branch not in {"auto", "a", "b", "c"}:
        raise ConfigValidationError(
            "BEACH constraint error: surface_current_model.zhao_branch must be "
            '"auto", "a", "b", or "c".'
        )
    source_scale = model_config.get("photoelectron_source_scale", 1.0)
    if (
        not isinstance(source_scale, (int, float))
        or isinstance(source_scale, bool)
        or not math.isfinite(float(source_scale))
        or float(source_scale) < 0.0
    ):
        raise ConfigValidationError(
            "BEACH constraint error: surface_current_model.photoelectron_source_scale "
            "must be finite and >= 0."
    )
    photoelectron_active = float(source_scale) > 0.0
    if not photoelectron_active and zhao_branch not in {"auto", "c"}:
        raise ConfigValidationError(
            "BEACH constraint error: photoelectron_source_scale=0 requires "
            'surface_current_model.zhao_branch to be "auto" or "c".'
        )
    photoelectron_keys = (
        "photoelectron_species",
        "solar_elevation_deg",
        "photoelectron_ref_density_m3",
    )
    if photoelectron_active:
        for key in ("solar_elevation_deg", "photoelectron_ref_density_m3"):
            value = model_config.get(key)
            if (
                not isinstance(value, (int, float))
                or isinstance(value, bool)
                or not math.isfinite(float(value))
                or float(value) <= 0.0
            ):
                raise ConfigValidationError(
                    f"BEACH constraint error: surface_current_model.{key} "
                    "must be finite and > 0."
                )
        if float(model_config["solar_elevation_deg"]) > 90.0:
            raise ConfigValidationError(
                "BEACH constraint error: surface_current_model.solar_elevation_deg "
                "must not exceed 90."
            )
    elif any(key in model_config for key in photoelectron_keys):
        raise ConfigValidationError(
            "BEACH constraint error: photoelectron_source_scale=0 requires omitting "
            "all photoelectron-specific Zhao settings."
        )
    if "reference_area_m2" in model_config:
        area = model_config["reference_area_m2"]
        if (
            not isinstance(area, (int, float))
            or isinstance(area, bool)
            or not math.isfinite(float(area))
            or float(area) <= 0.0
        ):
            raise ConfigValidationError(
                "BEACH constraint error: surface_current_model.reference_area_m2 "
                "must be finite and > 0."
            )
    elif domain is None:
        raise ConfigValidationError(
            "BEACH constraint error: surface_current_model requires reference_area_m2 "
            "or a finite domain."
        )

    by_key = {
        str(item.get("species_key", f"species_{index}")): item
        for index, item in enumerate(species, start=1)
        if item.get("enabled", True) is True
    }
    roles = ["electron", "ion"]
    if photoelectron_active:
        roles.append("photoelectron")
    role_keys = {role: model_config.get(f"{role}_species") for role in roles}
    if any(not isinstance(value, str) or not value for value in role_keys.values()):
        raise ConfigValidationError(
            "BEACH constraint error: surface_current_model requires its active species roles."
        )
    if len(set(role_keys.values())) != len(roles):
        raise ConfigValidationError(
            "BEACH constraint error: surface_current_model species references must be distinct."
        )
    if any(value not in by_key for value in role_keys.values()):
        raise ConfigValidationError(
            "BEACH constraint error: surface_current_model references an unknown or disabled species."
        )
    selected = {role: by_key[key] for role, key in role_keys.items()}
    for role, item in selected.items():
        if item.get("surface_charge_closure", "explicit") != "fixed_current":
            raise ConfigValidationError(
                f"BEACH constraint error: Zhao {role} species requires "
                'surface_charge_closure="fixed_current".'
            )
        if "target_absorbed_current_a" in item or "target_emission_current_a" in item:
            raise ConfigValidationError(
                "BEACH constraint error: automatic surface-current species cannot "
                "also specify manual target currents."
            )
    q_e = float(selected["electron"].get("q_particle", -1.602176634e-19))
    q_i = float(selected["ion"].get("q_particle", -1.602176634e-19))
    if q_e >= 0.0 or q_i <= 0.0:
        raise ConfigValidationError(
            "BEACH constraint error: Zhao roles require negative electron and positive ion species."
        )
    q_pe = None
    if photoelectron_active:
        q_pe = float(selected["photoelectron"].get("q_particle", -1.602176634e-19))
        if q_pe >= 0.0:
            raise ConfigValidationError(
                "BEACH constraint error: Zhao photoelectron species requires negative charge."
            )
    elementary_charge = 1.602176634e-19
    role_charges = [q_e, q_i]
    if q_pe is not None:
        role_charges.append(q_pe)
    if any(
        abs(abs(charge) - elementary_charge) > 1.0e-6 * elementary_charge
        for charge in role_charges
    ):
        raise ConfigValidationError(
            "BEACH constraint error: Zhao stationary current requires singly charged "
            "active species."
        )
    if photoelectron_active:
        electron_mass = float(selected["electron"].get("m_particle", 9.10938356e-31))
        photoelectron_mass = float(
            selected["photoelectron"].get("m_particle", 9.10938356e-31)
        )
        if abs(photoelectron_mass - electron_mass) > 1.0e-6 * electron_mass:
            raise ConfigValidationError(
                "BEACH constraint error: Zhao stationary current requires matching "
                "ambient-electron and photoelectron masses."
            )
        photo = selected["photoelectron"]
        if (
            photo.get("source_mode", "volume_seed") != "photo_raycast"
            or photo.get("deposit_opposite_charge_on_emit") is not True
            or photo.get("inject_face") != "z_high"
        ):
            raise ConfigValidationError(
                "BEACH constraint error: Zhao photoelectron species requires photo_raycast, "
                "opposite-charge emission deposit, and inject_face=z_high."
            )
        photo_boundary = photo.get("boundary", {})
        if not isinstance(photo_boundary, Mapping):
            photo_boundary = {}
        photo_z_high = photo_boundary.get("z_high", "inherit")
        if photo_z_high == "inherit":
            photo_z_high = (particle_boundary or {}).get("z_high", "open")
        if photo_z_high != "open":
            raise ConfigValidationError(
                "BEACH constraint error: Zhao photoelectron kinetic closure requires "
                "an open z-high particle boundary."
            )
    for role in ("electron", "ion"):
        item = selected[role]
        boundary_inflow = item.get("boundary_inflow", {})
        is_reservoir = (
            isinstance(boundary_inflow, Mapping)
            and boundary_inflow.get("z_high") == "reservoir"
        ) or (
            item.get("source_mode") == "reservoir_face"
            and item.get("inject_face") == "z_high"
        )
        drift = item.get("drift_velocity", [0.0, 0.0, -8.0e5])
        try:
            drift_z = (
                float(drift[2])
                if isinstance(drift, Sequence)
                and not isinstance(drift, (str, bytes))
                and len(drift) == 3
                else math.nan
            )
        except (TypeError, ValueError):
            drift_z = math.nan
        if not is_reservoir or not math.isfinite(drift_z) or drift_z >= 0.0:
            raise ConfigValidationError(
                f"BEACH constraint error: Zhao {role} species requires inward z-high "
                "reservoir drift."
            )

    ion = selected["ion"]
    if "number_density_cm3" in ion:
        ion_density_m3 = float(ion["number_density_cm3"]) * 1.0e6
    else:
        ion_density_m3 = float(ion.get("number_density_m3", 0.0))
    if not math.isfinite(ion_density_m3) or ion_density_m3 <= 0.0:
        raise ConfigValidationError(
            "BEACH constraint error: Zhao ion species requires a positive number density."
        )

    def temperature_k(item: Mapping[str, Any]) -> float:
        if "temperature_ev" in item:
            return float(item["temperature_ev"]) * 1.160451812e4
        return float(item.get("temperature_k", 2.0e4))

    electron_temperature_k = temperature_k(selected["electron"])
    ion_temperature_k = temperature_k(ion)
    if not math.isfinite(electron_temperature_k) or electron_temperature_k <= 0.0:
        raise ConfigValidationError(
            "BEACH constraint error: Zhao electron temperature must be positive."
        )
    if photoelectron_active:
        photoelectron_temperature_k = temperature_k(selected["photoelectron"])
        if (
            not math.isfinite(photoelectron_temperature_k)
            or photoelectron_temperature_k <= 0.0
        ):
            raise ConfigValidationError(
                "BEACH constraint error: Zhao photoelectron temperature must be positive."
            )
    if (
        not math.isfinite(ion_temperature_k)
        or ion_temperature_k > 0.1 * electron_temperature_k
    ):
        raise ConfigValidationError(
            "BEACH constraint error: Zhao stationary current requires cold ions "
            "with T_i <= 0.1 T_e."
        )


def _validate_matching_plane_model(
    model_config: Mapping[str, Any],
    *,
    species: list[Mapping[str, Any]],
    sim: Mapping[str, Any],
    domain: Mapping[str, Any] | None,
    field_boundary: Mapping[str, Any] | None,
    particle_boundary: Mapping[str, Any] | None,
    reservoir: Mapping[str, Any] | None,
    periodic2_config: object,
) -> None:
    stationary_zhao_keys = {
        "solar_elevation_deg",
        "photoelectron_ref_density_m3",
        "photoelectron_source_scale",
        "reference_area_m2",
    }
    if stationary_zhao_keys.intersection(model_config):
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic cannot use "
            "stationary-Zhao source settings or reference_area_m2."
        )

    response_backend = model_config.get("response_backend", "table")
    implicit_zero_mode = model_config.get("implicit_zero_mode", False)
    if (
        not isinstance(response_backend, str)
        or response_backend not in {"table", "zhao_online"}
    ):
        raise ConfigValidationError(
            "BEACH constraint error: surface_current_model.response_backend must "
            'be "table" or "zhao_online".'
        )
    if response_backend == "table":
        if {"zhao_branch", "zhao_root_selection"}.intersection(model_config):
            raise ConfigValidationError(
                "BEACH constraint error: matching_plane_quasistatic "
                'response_backend="table" cannot use Zhao-specific settings such '
                "as zhao_branch or zhao_root_selection."
            )
        response_path = model_config.get("response_table_path")
        if (
            not isinstance(response_path, str)
            or not response_path.strip()
            or len(response_path) > 256
        ):
            raise ConfigValidationError(
                "BEACH constraint error: matching_plane_quasistatic "
                'response_backend="table" requires a non-empty '
                "response_table_path of at most 256 characters."
            )
    else:
        if "response_table_path" in model_config:
            raise ConfigValidationError(
                "BEACH constraint error: matching_plane_quasistatic "
                'response_backend="zhao_online" cannot use response_table_path.'
            )
        zhao_branch = model_config.get("zhao_branch", "auto")
        if (
            not isinstance(zhao_branch, str)
            or zhao_branch not in {"auto", "a", "b", "c"}
        ):
            raise ConfigValidationError(
                "BEACH constraint error: surface_current_model.zhao_branch must be "
                '"auto", "a", "b", or "c".'
            )
        root_selection = model_config.get("zhao_root_selection", "require_unique")
        if (
            not isinstance(root_selection, str)
            or root_selection
            not in {"require_unique", "minimum_energy", "continuation"}
        ):
            raise ConfigValidationError(
                "BEACH constraint error: surface_current_model.zhao_root_selection "
                'must be "require_unique", "minimum_energy", or "continuation".'
            )
        if root_selection == "continuation" and (
            zhao_branch != "a" or implicit_zero_mode is not True
        ):
            raise ConfigValidationError(
                "BEACH constraint error: surface_current_model."
                'zhao_root_selection="continuation" requires '
                'response_backend="zhao_online", zhao_branch="a", and '
                "implicit_zero_mode=true."
            )
    for key, default in (("coupling_rtol", 1.0e-4), ("coupling_relaxation", 0.5)):
        value = model_config.get(key, default)
        if (
            not isinstance(value, (int, float))
            or isinstance(value, bool)
            or not math.isfinite(float(value))
            or not 0.0 < float(value) <= 1.0
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: surface_current_model.{key} must be "
                "finite and in (0, 1]."
            )
    coupling_atol = model_config.get("coupling_atol", [0.0, 0.0, 0.0, 0.0])
    if (
        not isinstance(coupling_atol, Sequence)
        or isinstance(coupling_atol, (str, bytes))
        or len(coupling_atol) != 4
        or any(
            not isinstance(value, (int, float))
            or isinstance(value, bool)
            or not math.isfinite(float(value))
            or float(value) < 0.0
            for value in coupling_atol
        )
    ):
        raise ConfigValidationError(
            "BEACH constraint error: surface_current_model.coupling_atol must "
            "contain four finite values >= 0."
        )
    if response_backend == "zhao_online" and any(
        float(value) > 0.0 for value in coupling_atol[2:4]
    ):
        raise ConfigValidationError(
            "BEACH constraint error: zhao_online coupling_atol must be zero on "
            "inactive ambient-outward feedback axes."
        )
    max_iterations = model_config.get("coupling_max_iterations", 20)
    if (
        not isinstance(max_iterations, int)
        or isinstance(max_iterations, bool)
        or max_iterations < 1
    ):
        raise ConfigValidationError(
            "BEACH constraint error: surface_current_model.coupling_max_iterations "
            "must be an integer >= 1."
        )

    if domain is None or field_boundary is None:
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic requires a "
            "periodic2 [domain] box."
        )
    if field_boundary.get("mode", "free") != "periodic2":
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic requires "
            'field_boundary.mode="periodic2".'
        )
    if set(domain.get("periodic_axes", [])) != {"x", "y"}:
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic requires x/y "
            "periodic and z-open box topology."
        )
    if not isinstance(periodic2_config, Mapping):
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic requires an "
            "explicit split-zero-mode [periodic2] table."
        )
    if periodic2_config.get("nonzero_mode_backend") not in {
        "cached_kneq0",
        "panel_spectral_reference",
    }:
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic requires a split "
            "periodic2 nonzero-mode backend."
        )
    if periodic2_config.get("zero_mode_policy") != "exclude_k0":
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic requires "
            'periodic2.zero_mode_policy="exclude_k0".'
        )
    if periodic2_config.get("lower_boundary_model") not in {
        "e_bottom_zero",
        "symmetric_vacuum",
    }:
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic requires a "
            "supported periodic2 lower boundary model."
        )
    if (
        implicit_zero_mode
        and periodic2_config.get("lower_boundary_model") != "e_bottom_zero"
    ):
        raise ConfigValidationError(
            "BEACH constraint error: surface_current_model.implicit_zero_mode "
            'requires periodic2.lower_boundary_model="e_bottom_zero".'
        )

    if "e0" in sim:
        e0 = _maybe_vec3(sim.get("e0"), name="sim.e0")
        e0_is_zero = e0 is not None and all(value == 0.0 for value in e0)
    else:
        e0_abs = sim.get("e0_abs", 0.0)
        e0_is_zero = (
            isinstance(e0_abs, (int, float))
            and not isinstance(e0_abs, bool)
            and math.isfinite(float(e0_abs))
            and float(e0_abs) == 0.0
        )
    b0 = _maybe_vec3(sim.get("b0", [0.0, 0.0, 0.0]), name="sim.b0")
    if not e0_is_zero or b0 is None or any(value != 0.0 for value in b0):
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic requires "
            "sim.e0=sim.b0=[0,0,0]."
        )
    if (
        reservoir is not None
        and reservoir.get("inflow_model", "source_vdf") != "source_vdf"
    ):
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic cannot use the "
            "generic reservoir potential model."
        )
    if (particle_boundary or {}).get("ordinary_open_model", "escape") != "escape":
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic requires "
            'particle_boundary.ordinary_open_model="escape".'
        )
    by_key = {
        str(item.get("species_key", f"species_{index}")): item
        for index, item in enumerate(species, start=1)
        if item.get("enabled", True) is True
    }
    role_keys = {
        role: model_config.get(f"{role}_species") for role in ("electron", "ion")
    }
    if any(not isinstance(value, str) or not value for value in role_keys.values()):
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic requires electron "
            "and ion species roles."
        )
    photoelectron_key = model_config.get("photoelectron_species")
    if photoelectron_key is not None:
        if not isinstance(photoelectron_key, str) or not photoelectron_key:
            raise ConfigValidationError(
                "BEACH constraint error: matching_plane_quasistatic "
                "photoelectron_species must be a non-empty string when provided."
            )
        role_keys["photoelectron"] = photoelectron_key
    photoelectron_active = "photoelectron" in role_keys
    if len(set(role_keys.values())) != len(role_keys):
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic species references "
            "must be distinct."
        )
    if any(value not in by_key for value in role_keys.values()):
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic references an "
            "unknown or disabled species."
        )
    selected = {role: by_key[key] for role, key in role_keys.items()}
    for role, item in selected.items():
        if item.get("surface_charge_closure", "explicit") != "explicit":
            raise ConfigValidationError(
                f"BEACH constraint error: matching_plane_quasistatic {role} species "
                'requires surface_charge_closure="explicit".'
            )
    if any(
        item.get("surface_charge_closure", "explicit") == "fixed_current"
        or "target_absorbed_current_a" in item
        or "target_emission_current_a" in item
        for item in species
        if item.get("enabled", True) is True
    ):
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic cannot use manual "
            "fixed_current targets on any enabled species."
        )
    if len(by_key) != len(role_keys):
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic requires exactly "
            "its enabled electron, ion, and optional photoelectron roles."
        )

    charges: dict[str, float] = {}
    for role, item in selected.items():
        raw_charge = item.get("q_particle", -1.602176634e-19)
        if (
            not isinstance(raw_charge, (int, float))
            or isinstance(raw_charge, bool)
            or not math.isfinite(float(raw_charge))
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: matching_plane_quasistatic {role} charge "
                "must be finite and numeric."
            )
        charges[role] = float(raw_charge)
    if charges["electron"] >= 0.0 or charges["ion"] <= 0.0:
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic requires a negative "
            "electron and positive ion species."
        )
    if photoelectron_active and charges["photoelectron"] >= 0.0:
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic requires a negative "
            "photoelectron species."
        )

    for role in ("electron", "ion"):
        ambient = selected[role]
        boundary_inflow = ambient.get("boundary_inflow", {})
        npcls_per_step = ambient.get("npcls_per_step", 0)
        if (
            ambient.get("source_mode", "volume_seed") != "volume_seed"
            or not isinstance(npcls_per_step, int)
            or isinstance(npcls_per_step, bool)
            or npcls_per_step != 0
            or not isinstance(boundary_inflow, Mapping)
            or dict(boundary_inflow) != {"z_high": "reservoir"}
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: matching_plane_quasistatic {role} species "
                "requires volume_seed, npcls_per_step=0, and only z-high "
                'boundary_inflow="reservoir".'
            )
    if photoelectron_active:
        photo = selected["photoelectron"]
        if (
            photo.get("source_mode", "volume_seed") != "photo_raycast"
            or photo.get("deposit_opposite_charge_on_emit") is not True
            or photo.get("inject_face") != "z_high"
        ):
            raise ConfigValidationError(
                "BEACH constraint error: matching_plane_quasistatic photoelectrons "
                "require opposite-deposit photo_raycast from z_high."
            )

    expected_boundaries = {
        "x_low": "periodic",
        "x_high": "periodic",
        "y_low": "periodic",
        "y_high": "periodic",
        "z_low": "open",
        "z_high": "open",
    }
    species_indices = {
        str(item.get("species_key", f"species_{index}")): index
        for index, item in enumerate(species, start=1)
    }
    for role, item in selected.items():
        effective_boundary = _validate_species_particle_boundary(
            item,
            index=species_indices[str(role_keys[role])],
            periodic_axes={"x", "y"},
            global_boundary=particle_boundary or {},
        )
        if effective_boundary != expected_boundaries:
            raise ConfigValidationError(
                f"BEACH constraint error: matching_plane_quasistatic {role} species "
                "requires x/y periodic and z-low/z-high open particle boundaries."
            )

    if response_backend == "zhao_online":
        _validate_matching_plane_zhao_online(selected, charges)


def _validate_matching_plane_zhao_online(
    selected: Mapping[str, Mapping[str, Any]],
    charges: Mapping[str, float],
) -> None:
    elementary_charge = 1.602176634e-19
    if any(
        not math.isfinite(charge)
        or abs(abs(charge) - elementary_charge) > 1.0e-6 * elementary_charge
        for charge in charges.values()
    ):
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic zhao_online "
            "requires singly charged role species."
        )

    def finite_float(value: object) -> float | None:
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            return None
        numeric = float(value)
        return numeric if math.isfinite(numeric) else None

    electron_mass = finite_float(
        selected["electron"].get("m_particle", 9.10938356e-31)
    )
    photoelectron = selected.get("photoelectron")
    photoelectron_mass = (
        finite_float(photoelectron.get("m_particle", 9.10938356e-31))
        if photoelectron is not None
        else None
    )
    ion_mass = finite_float(selected["ion"].get("m_particle", 9.10938356e-31))
    if (
        electron_mass is None
        or electron_mass <= 0.0
        or ion_mass is None
        or ion_mass <= 0.0
        or (
            photoelectron is not None
            and (photoelectron_mass is None or photoelectron_mass <= 0.0)
        )
    ):
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic zhao_online "
            "requires positive role-species masses."
        )
    if photoelectron is not None and (
        photoelectron_mass is None
        or abs(photoelectron_mass - electron_mass) > 1.0e-6 * electron_mass
    ):
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic zhao_online "
            "requires matching ambient-electron and photoelectron masses."
        )

    def temperature_k(item: Mapping[str, Any]) -> float | None:
        if "temperature_ev" in item:
            value = finite_float(item["temperature_ev"])
            return None if value is None else value * 1.160451812e4
        return finite_float(item.get("temperature_k", 2.0e4))

    electron_temperature = temperature_k(selected["electron"])
    photoelectron_temperature = (
        temperature_k(photoelectron) if photoelectron is not None else None
    )
    ion_temperature = temperature_k(selected["ion"])
    if electron_temperature is None or electron_temperature <= 0.0:
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic zhao_online "
            "requires a positive electron temperature."
        )
    if photoelectron is not None and (
        photoelectron_temperature is None or photoelectron_temperature <= 0.0
    ):
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic zhao_online "
            "requires a positive photoelectron temperature."
        )
    if (
        ion_temperature is None
        or ion_temperature < 0.0
        or ion_temperature > 0.1 * electron_temperature
    ):
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic zhao_online "
            "requires cold ions with T_i <= 0.1 T_e."
        )

    for role in ("electron", "ion"):
        drift = selected[role].get("drift_velocity", [0.0, 0.0, -8.0e5])
        drift_components: list[float] | None = None
        if (
            isinstance(drift, Sequence)
            and not isinstance(drift, (str, bytes))
            and len(drift) == 3
        ):
            parsed = [finite_float(component) for component in drift]
            if all(component is not None for component in parsed):
                drift_components = [float(component) for component in parsed]
        if drift_components is None or drift_components[2] >= 0.0:
            raise ConfigValidationError(
                "BEACH constraint error: matching_plane_quasistatic zhao_online "
                "requires finite drift vectors and positive ambient inward drift at z-high."
            )

    ion = selected["ion"]
    if "number_density_cm3" in ion:
        ion_density = finite_float(ion["number_density_cm3"])
        if ion_density is not None:
            ion_density *= 1.0e6
    else:
        ion_density = finite_float(ion.get("number_density_m3", 0.0))
    if ion_density is None or ion_density <= 0.0:
        raise ConfigValidationError(
            "BEACH constraint error: matching_plane_quasistatic zhao_online "
            "requires a positive ion number density."
        )

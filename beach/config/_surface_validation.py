"""Validate the Zhao zero-current surface-current model."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

from ._shared import ConfigValidationError


def _validate_surface_current_model(
    model_config: Mapping[str, Any] | None,
    *,
    species: list[Mapping[str, Any]],
    domain: Mapping[str, Any] | None,
    field_boundary: Mapping[str, Any] | None,
    particle_boundary: Mapping[str, Any] | None,
    periodic2_config: object,
) -> None:
    if model_config is None:
        return
    model = model_config.get("model")
    if model is None or model == "none":
        if set(model_config) - {"model"}:
            raise ConfigValidationError(
                'BEACH constraint error: surface_current_model.model="none" cannot '
                "use Zhao settings."
            )
        return
    if model != "zhao_stationary":
        raise ConfigValidationError(
            'BEACH constraint error: surface_current_model.model must be "none" or '
            '"zhao_stationary".'
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
    refresh_batches = model_config.get("outflow_refresh_batches", 0)
    if (
        not isinstance(refresh_batches, int)
        or isinstance(refresh_batches, bool)
        or refresh_batches < 0
    ):
        raise ConfigValidationError(
            "BEACH constraint error: surface_current_model.outflow_refresh_batches "
            "must be an integer >= 0."
        )
    if refresh_batches > 0:
        _validate_zhao_outflow_refresh(
            photoelectron_active=photoelectron_active,
            domain=domain,
            field_boundary=field_boundary,
            periodic2_config=periodic2_config,
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
        # A nondrifting electron reservoir is a valid Zhao state; cold ions need a beam.
        outward = drift_z > 0.0 if role == "electron" else drift_z >= 0.0
        if not is_reservoir or not math.isfinite(drift_z) or outward:
            raise ConfigValidationError(
                f"BEACH constraint error: Zhao {role} species requires inward z-high "
                "reservoir drift."
            )

    def number_density_m3(item: Mapping[str, Any]) -> float:
        if "number_density_cm3" in item:
            return float(item["number_density_cm3"]) * 1.0e6
        return float(item.get("number_density_m3", 0.0))

    ion = selected["ion"]
    ion_density_m3 = number_density_m3(ion)
    if not math.isfinite(ion_density_m3) or ion_density_m3 <= 0.0:
        raise ConfigValidationError(
            "BEACH constraint error: Zhao ion species requires a positive number density."
        )
    # Both densities name the solar wind at infinity; the root sets the injected electron density.
    electron_density_m3 = number_density_m3(selected["electron"])
    if not abs(electron_density_m3 - ion_density_m3) <= 1.0e-9 * ion_density_m3:
        raise ConfigValidationError(
            "BEACH constraint error: Zhao ambient electron and ion species must share "
            "the solar-wind number density; the injected electron reservoir density "
            "is derived from the Zhao root."
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


def _validate_zhao_outflow_refresh(
    *,
    photoelectron_active: bool,
    domain: Mapping[str, Any] | None,
    field_boundary: Mapping[str, Any] | None,
    periodic2_config: object,
) -> None:
    """Require a periodic cell whose z-high mean potential can carry the outer wall potential."""
    prefix = "BEACH constraint error: surface_current_model.outflow_refresh_batches"
    if not photoelectron_active:
        raise ConfigValidationError(
            f"{prefix} requires photoelectron_source_scale > 0."
        )
    if (
        domain is None
        or field_boundary is None
        or field_boundary.get("mode", "free") != "periodic2"
    ):
        raise ConfigValidationError(f"{prefix} requires a periodic2 [domain] box.")
    if set(domain.get("periodic_axes", [])) != {"x", "y"}:
        raise ConfigValidationError(f"{prefix} requires x/y periodic axes.")
    # The schema restricts an explicit [periodic2] table to split backends with exclude_k0.
    if not isinstance(periodic2_config, Mapping):
        raise ConfigValidationError(
            f"{prefix} requires an explicit split-zero-mode [periodic2] table."
        )

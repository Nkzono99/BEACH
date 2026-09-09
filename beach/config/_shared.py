"""Shared configuration errors, key sets, and structural/value validation."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

CONFIG_FILENAME = "beach.toml"

SCHEMA_BASE_URL = "https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas"

BEACH_SCHEMA_URL = f"{SCHEMA_BASE_URL}/beach.schema.json"

TOP_LEVEL_CONFIG_ORDER = (
    "sim",
    "domain",
    "field_boundary",
    "particle_boundary",
    "reservoir",
    "surface_current_model",
    "periodic2",
    "particles",
    "mesh",
    "output",
)

SCHEMA_DIRECTIVE = f"#:schema {BEACH_SCHEMA_URL}"

_RESERVED_TOP_LEVEL_KEYS = frozenset(
    {"schema_version", "title", "use_presets", "override", "base_case"}
)

_FACE_SOURCE_MODES = frozenset({"reservoir_face", "photo_raycast"})

_RESERVOIR_SOURCE_MODES = frozenset({"reservoir_face", "plane_source"})

_FACE_KEYS = frozenset({"x_low", "x_high", "y_low", "y_high", "z_low", "z_high"})

_REMOVED_SIM_KEYS = frozenset(
    {
        "reservoir_potential_model",
        "open_boundary_model",
        "field_bc_mode",
        "phi_infty",
        "injection_face_phi_grid_n",
        "use_box",
        "box_min",
        "box_max",
        "box_origin",
        "box_size",
        "bc_x_low",
        "bc_x_high",
        "bc_y_low",
        "bc_y_high",
        "bc_z_low",
        "bc_z_high",
        "sheath_alpha_deg",
        "sheath_photoelectron_ref_density_cm3",
        "sheath_electron_drift_mode",
        "sheath_ion_drift_mode",
        "sheath_injection_model",
        "sheath_reference_coordinate",
        "multiple_box_events_soft_discard_count_limit",
        "softening",
    }
)

class ConfigError(ValueError):
    """Base error for BEACH config handling."""


class ConfigValidationError(ConfigError):
    """Raised when ``beach.toml`` violates known BEACH constraints."""


def _validate_legacy_keys(document: Mapping[str, Any]) -> None:
    """Keep migration diagnostics; the schema owns the set of accepted keys."""

    forbidden = sorted(set(document) & _RESERVED_TOP_LEVEL_KEYS)
    if forbidden:
        raise ConfigError(
            "config error: reserved top-level key(s) are not allowed: "
            + ", ".join(forbidden)
        )
    sim = document.get("sim")
    if isinstance(sim, Mapping):
        removed = sorted(set(sim) & _REMOVED_SIM_KEYS)
        if removed:
            raise ConfigValidationError(
                "BEACH constraint error: removed sim key(s): " + ", ".join(removed) + "."
            )


def _validate_high_level_fragment(
    document: Mapping[str, Any],
    *,
    context: str,
) -> None:
    domain = document.get("domain")
    if isinstance(domain, Mapping):
        if "box_origin" in domain and "box_min" in domain:
            raise ConfigError(
                f"{context} error: domain.box_origin and domain.box_min cannot be specified "
                "in the same fragment."
            )
        if "box_size" in domain and "box_max" in domain:
            raise ConfigError(
                f"{context} error: domain.box_size and domain.box_max cannot be specified "
                "in the same fragment."
            )

    particles = document.get("particles")
    species = particles.get("species") if isinstance(particles, Mapping) else None
    if isinstance(species, list):
        for index, item in enumerate(species, start=1):
            if not isinstance(item, Mapping):
                continue
            if not any(
                key in item for key in ("inject_region_mode", "uv_low", "uv_high")
            ):
                continue
            source_mode = item.get("source_mode", "volume_seed")
            if (
                not isinstance(source_mode, str)
                or source_mode not in _FACE_SOURCE_MODES
            ):
                raise ConfigError(
                    f"{context} error: particles.species[{index}] uses inject_region_mode/uv_* "
                    'but source_mode must be "reservoir_face" or "photo_raycast".'
                )


def _optional_runtime_table(
    document: Mapping[str, Any], key: str
) -> Mapping[str, Any] | None:
    value = document.get(key)
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise ConfigValidationError(f"BEACH constraint error: [{key}] must be a table.")
    return value


def _maybe_vec3(value: object, *, name: str) -> list[float] | None:
    if value is None:
        return None
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != 3
    ):
        raise ConfigValidationError(
            f"BEACH constraint error: {name} must be a 3-element array."
        )
    out: list[float] = []
    for item in value:
        if not isinstance(item, (int, float)) or isinstance(item, bool):
            raise ConfigValidationError(
                f"BEACH constraint error: {name} must contain only numeric values."
            )
        numeric = float(item)
        if not math.isfinite(numeric):
            raise ConfigValidationError(
                f"BEACH constraint error: {name} must contain only finite values."
            )
        out.append(numeric)
    return out


def _is_array_of_tables(value: object) -> bool:
    return (
        isinstance(value, list)
        and len(value) > 0
        and all(isinstance(item, Mapping) for item in value)
    )

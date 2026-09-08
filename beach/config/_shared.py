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

_REQUIRED_RUNTIME_TABLES = ("sim", "particles", "mesh", "output")

SCHEMA_DIRECTIVE = f"#:schema {BEACH_SCHEMA_URL}"

_FRAGMENT_TOP_LEVEL_KEYS = frozenset(TOP_LEVEL_CONFIG_ORDER)

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


def _validate_id_fields(
    items: list[Any],
    *,
    context: str,
    table_name: str,
) -> None:
    """Validate that ``id`` fields, when present, are non-empty strings."""
    for index, item in enumerate(items, start=1):
        if not isinstance(item, Mapping):
            continue
        item_id = item.get("id")
        if item_id is None:
            continue
        if not isinstance(item_id, str) or not item_id:
            raise ConfigError(
                f"{context} error: {table_name}[{index}].id must be a non-empty string."
            )


def _validate_fragment_structure(
    document: Mapping[str, Any],
    *,
    context: str,
    allow_meta_keys: bool,
) -> None:
    unknown_keys = [key for key in document if key not in _FRAGMENT_TOP_LEVEL_KEYS]
    if not allow_meta_keys:
        forbidden = [key for key in document if key in _RESERVED_TOP_LEVEL_KEYS]
        if forbidden:
            raise ConfigError(
                f"{context} error: reserved top-level key(s) are not allowed: "
                + ", ".join(sorted(forbidden))
            )
    if unknown_keys:
        raise ConfigError(
            f"{context} error: unsupported top-level key(s): "
            + ", ".join(sorted(unknown_keys))
        )

    for key in _FRAGMENT_TOP_LEVEL_KEYS.intersection(document):
        value = document[key]
        if not isinstance(value, Mapping):
            raise ConfigError(
                f"{context} error: top-level key {key!r} must be a table."
            )

    particles = document.get("particles")
    if isinstance(particles, Mapping) and "species" in particles:
        species = particles["species"]
        if not isinstance(species, list) or not all(
            isinstance(item, Mapping) for item in species
        ):
            raise ConfigError(
                f"{context} error: particles.species must be an array of tables."
            )
        _validate_id_fields(species, context=context, table_name="particles.species")

    mesh = document.get("mesh")
    if isinstance(mesh, Mapping) and "templates" in mesh:
        templates = mesh["templates"]
        if not isinstance(templates, list) or not all(
            isinstance(item, Mapping) for item in templates
        ):
            raise ConfigError(
                f"{context} error: mesh.templates must be an array of tables."
            )
        _validate_id_fields(templates, context=context, table_name="mesh.templates")


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


def _require_table(
    document: Mapping[str, Any],
    key: str,
    *,
    context: str,
) -> dict[str, Any]:
    value = document.get(key)
    if not isinstance(value, Mapping):
        raise ConfigValidationError(
            f"BEACH constraint error: {context} requires [{key}] to be a table."
        )
    return dict(value)


def _optional_runtime_table(
    document: Mapping[str, Any], key: str
) -> dict[str, Any] | None:
    value = document.get(key)
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise ConfigValidationError(f"BEACH constraint error: [{key}] must be a table.")
    return dict(value)


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

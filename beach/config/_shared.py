"""Shared configuration errors, key sets, and structural/value validation."""

from __future__ import annotations

from collections.abc import Mapping
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

_FACE_SOURCE_MODES = frozenset({"photo_raycast"})


class ConfigError(ValueError):
    """Base error for BEACH config handling."""


class ConfigValidationError(ConfigError):
    """Raised when ``beach.toml`` violates known BEACH constraints."""


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
                    'but source_mode must be "photo_raycast".'
                )


def _is_array_of_tables(value: object) -> bool:
    return (
        isinstance(value, list)
        and len(value) > 0
        and all(isinstance(item, Mapping) for item in value)
    )

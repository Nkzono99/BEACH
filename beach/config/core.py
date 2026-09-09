"""Public BEACH configuration loading, rendering, and semantic comparison."""

from __future__ import annotations

import copy
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ._authoring import (
    _resolve_domain_high_level,
    _resolve_mesh_high_level,
    _resolve_sim_high_level,
    _resolve_species_high_level,
    _resolved_box_bounds,
)
from ._runtime_validation import (
    _validate_runtime_semantics,
)
from ._shared import (
    BEACH_SCHEMA_URL as BEACH_SCHEMA_URL,
)
from ._shared import (
    CONFIG_FILENAME as CONFIG_FILENAME,
)
from ._shared import (
    SCHEMA_BASE_URL as SCHEMA_BASE_URL,
)
from ._shared import (
    SCHEMA_DIRECTIVE as SCHEMA_DIRECTIVE,
)
from ._shared import (
    TOP_LEVEL_CONFIG_ORDER as TOP_LEVEL_CONFIG_ORDER,
)
from ._shared import (
    ConfigError as ConfigError,
)
from ._shared import (
    ConfigValidationError as ConfigValidationError,
)
from ._shared import (
    _is_array_of_tables,
    _validate_legacy_keys,
    _validate_high_level_fragment,
)
from ._toml import load_toml_file, render_toml_document
from .schema import (
    ConfigSchemaError,
    load_schema,
    prepare_schema_document,
    reject_basic_schema_errors,
    validation_errors,
)


def default_config() -> dict[str, Any]:
    """Return the small multi-batch surface-charging tutorial config."""

    return {
        "sim": {
            "dt": 5.0e-8,
            "batch_count": 20,
            "max_step": 80,
            "rng_seed": 12345,
            "field_solver": "direct",
        },
        "domain": {
            "box_min": [0.0, 0.0, 0.0],
            "box_max": [1.0, 1.0, 1.0],
            "periodic_axes": [],
        },
        "field_boundary": {"mode": "free"},
        "particle_boundary": {"z_low": "open", "z_high": "open"},
        "particles": {
            "species": [
                {
                    "source_mode": "volume_seed",
                    "q_particle": -1.602176634e-19,
                    "m_particle": 9.10938356e-31,
                    "w_particle": 2.0e5,
                    "npcls_per_step": 200,
                    "pos_low": [0.15, 0.15, 0.8],
                    "pos_high": [0.85, 0.85, 0.8],
                    "drift_velocity": [0.0, 0.0, -1.0e6],
                    "temperature_k": 0.0,
                },
            ]
        },
        "mesh": {
            "mode": "template",
            "templates": [
                {
                    "kind": "plane",
                    "enabled": True,
                    "surface_model": "insulator",
                    "surface_side": "normal_plus",
                    "size_x": 0.8,
                    "size_y": 0.8,
                    "nx": 12,
                    "ny": 12,
                    "center": [0.5, 0.5, 0.2],
                }
            ],
        },
        "output": {
            "write_files": True,
            "dir": "outputs/tutorial",
            "history_stride": 1,
            "write_mesh_potential": True,
            "write_potential_history": True,
        },
    }


def load_config_file(path: str | Path) -> dict[str, Any]:
    """Load, normalize, and validate one direct ``beach.toml`` file."""

    return _load_config_file(path)


def _load_config_file(
    path: str | Path, *, schema: Mapping[str, Any] | None = None,
    additional_schema: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    return _validate_config_document(
        load_toml_file(path), schema=schema, additional_schema=additional_schema,
        config_path=Path(path), normalize=True,
    )


def _resolve_surface_response_table_path(
    config: dict[str, Any], *, config_path: Path
) -> None:
    """Resolve a relative outer-response path from the config directory."""

    model = config.get("surface_current_model")
    if not isinstance(model, dict):
        return
    raw_path = model.get("response_table_path")
    if not isinstance(raw_path, str) or not raw_path.strip():
        return
    if Path(raw_path).is_absolute():
        return
    model["response_table_path"] = str(config_path.parent / raw_path)


def normalize_config_document(config: Mapping[str, Any]) -> dict[str, Any]:
    """Resolve high-level authoring notation and validate the runtime config."""

    return _validate_config_document(config, normalize=True)


def validate_runtime_config(config: Mapping[str, Any]) -> None:
    """Validate a runtime document without changing the caller's mapping."""

    _validate_config_document(config, normalize=False)


def _validate_config_document(
    config: Mapping[str, Any], *, normalize: bool,
    schema: Mapping[str, Any] | None = None, config_path: Path | None = None,
    additional_schema: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    if schema is None:
        schema, _ = load_schema()
    document = prepare_schema_document(config, schema)
    _validate_legacy_keys(document)
    phase = "authoring" if normalize else "normalized"
    input_errors = validation_errors(document, schema)
    reject_basic_schema_errors(input_errors, phase=phase)
    if normalize:
        _validate_high_level_fragment(document, context="config")
        resolved = normalize_high_level_config(document)
    else:
        resolved = document
    if config_path is not None:
        _resolve_surface_response_table_path(resolved, config_path=config_path)
    # Authoring arithmetic and path resolution can create new invalid values.
    resolved = prepare_schema_document(resolved, schema)
    resolved_errors = validation_errors(resolved, schema)
    reject_basic_schema_errors(resolved_errors, phase="normalized")
    _validate_runtime_semantics(resolved)
    if input_errors:
        raise ConfigSchemaError(input_errors, phase=phase)
    if resolved_errors:
        raise ConfigSchemaError(resolved_errors, phase="normalized")
    if additional_schema is not None:
        for candidate, candidate_phase in ((document, phase), (resolved, "normalized")):
            errors = validation_errors(candidate, additional_schema)
            if errors:
                raise ConfigSchemaError(errors, phase=candidate_phase)
    return resolved


def _strip_id_fields(config: dict[str, Any]) -> None:
    """Remove direct-config-only ``id`` keys before writing the final config."""
    _ARRAY_OF_TABLE_PATHS = (
        ("particles", "species"),
        ("mesh", "templates"),
    )
    for path in _ARRAY_OF_TABLE_PATHS:
        table = config
        for segment in path[:-1]:
            table = table.get(segment, {})
            if not isinstance(table, Mapping):
                break
        else:
            items = table.get(path[-1])
            if isinstance(items, list):
                for item in items:
                    if isinstance(item, dict):
                        item.pop("id", None)


def normalize_high_level_config(config: Mapping[str, Any]) -> dict[str, Any]:
    """Resolve high-level spatial notation into runtime beach.toml values."""

    resolved = copy.deepcopy(dict(config))
    sim = resolved.get("sim")
    if isinstance(sim, Mapping):
        resolved["sim"] = _resolve_sim_high_level(dict(sim))
    domain = resolved.get("domain")
    if isinstance(domain, Mapping):
        resolved["domain"] = _resolve_domain_high_level(dict(domain))
    box_min, box_max = _resolved_box_bounds(resolved)

    particles = resolved.get("particles")
    if isinstance(particles, Mapping):
        species = particles.get("species")
        if isinstance(species, list):
            resolved_species = [
                _resolve_species_high_level(
                    dict(item), box_min=box_min, box_max=box_max
                )
                for item in species
            ]
            particles_dict = dict(particles)
            particles_dict["species"] = resolved_species
            resolved["particles"] = particles_dict

    mesh = resolved.get("mesh")
    if isinstance(mesh, Mapping):
        resolved["mesh"] = _resolve_mesh_high_level(
            dict(mesh), box_min=box_min, box_max=box_max
        )

    _strip_id_fields(resolved)
    return resolved


def dump_beach_toml(
    config: Mapping[str, Any],
    *,
    source_config: str | Path | None = None,
) -> str:
    """Dump one validated config to ``beach.toml`` text."""

    header_comments = [SCHEMA_DIRECTIVE, "# Generated by beachx config init"]
    if source_config is not None:
        header_comments.append(f"# source_config={source_config}")
    return render_toml_document(
        config,
        header_comments=header_comments,
        top_level_order=TOP_LEVEL_CONFIG_ORDER,
    )


def semantic_diff(left: Any, right: Any) -> list[str]:
    """Return human-readable semantic differences between two TOML payloads."""

    lines: list[str] = []
    _append_semantic_diff(lines, (), left, right)
    return lines


def _append_semantic_diff(
    lines: list[str],
    path: tuple[str | int, ...],
    left: Any,
    right: Any,
) -> None:
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        ordered_keys: list[str] = list(left.keys())
        for key in right.keys():
            if key not in left:
                ordered_keys.append(key)
        for key in ordered_keys:
            in_left = key in left
            in_right = key in right
            next_path = (*path, key)
            if in_left and in_right:
                _append_semantic_diff(lines, next_path, left[key], right[key])
            elif in_left:
                lines.append(
                    f"- {_format_path(next_path)} = {_summarize_value(left[key])}"
                )
            else:
                lines.append(
                    f"+ {_format_path(next_path)} = {_summarize_value(right[key])}"
                )
        return

    if _is_array_of_tables(left) and _is_array_of_tables(right):
        left_items = list(left)
        right_items = list(right)
        limit = max(len(left_items), len(right_items))
        for index in range(limit):
            next_path = (*path, index)
            if index < len(left_items) and index < len(right_items):
                _append_semantic_diff(
                    lines, next_path, left_items[index], right_items[index]
                )
            elif index < len(left_items):
                lines.append(
                    f"- {_format_path(next_path)} = {_summarize_value(left_items[index])}"
                )
            else:
                lines.append(
                    f"+ {_format_path(next_path)} = {_summarize_value(right_items[index])}"
                )
        return

    if left != right:
        lines.append(
            f"{_format_path(path)}: {_summarize_value(left)} -> {_summarize_value(right)}"
        )


def _format_path(path: tuple[str | int, ...]) -> str:
    if not path:
        return "<root>"
    fragments: list[str] = []
    for part in path:
        if isinstance(part, int):
            fragments[-1] = f"{fragments[-1]}[{part}]"
        else:
            fragments.append(part)
    return ".".join(fragments)


def _summarize_value(value: Any) -> str:
    if isinstance(value, str):
        return repr(value)
    return str(value)

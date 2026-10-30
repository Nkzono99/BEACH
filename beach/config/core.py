"""Public BEACH configuration loading, rendering, and semantic comparison."""

from __future__ import annotations

import copy
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ._layout import TOP_LEVEL_ORDER, is_grouped, to_grouped_layout, to_runtime_layout

from ._authoring import (
    _resolve_box_shorthand,
    _resolve_mesh_high_level,
    _resolve_species_high_level,
    _resolved_box_bounds,
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
    _validate_high_level_fragment,
)
from ._toml import load_toml_file, render_toml_document
from .schema import (
    ConfigSchemaError,
    load_schema,
    prepare_schema_document,
    reject_basic_schema_errors,
    schema_validator,
)


def default_config() -> dict[str, Any]:
    """Return the official tutorial config: a floating plate charged by an electron-proton plasma."""

    def plasma_species(key: str, charge: float, mass: float, target: int) -> dict[str, Any]:
        return {
            "species_key": key,
            "q_particle": charge,
            "m_particle": mass,
            "number_density_cm3": 1.0,
            "temperature_ev": 1.0,
            "drift_velocity": [0.0, 0.0, 0.0],
            "target_macro_particles_per_batch": target,
            "boundary_inflow": {"z_high": "reservoir"},
        }

    return {
        "sim": {
            "batch_count": 10000,
            "rng_seed": 12345,
            "batch_duration": 5.0e-4,
            "dt": 1.0e-8,
            "max_step": 3000,
            "field_solver": "fmm",
            "field_periodic_far_correction": "cached_kneq0",
            "field_periodic_image_layers": 1,
        },
        "domain": {
            "box_min": [0.0, 0.0, 0.0],
            "box_max": [0.1, 0.1, 0.11],
            "periodic_axes": ["x", "y"],
        },
        "mesh": {
            "parts": [
                {
                    "kind": "plane",
                    "enabled": True,
                    "surface_model": "conductor",
                    "surface_side": "normal_plus",
                    "size_x": 0.1,
                    "size_y": 0.1,
                    "nx": 2,
                    "ny": 2,
                    "center": [0.05, 0.05, 0.01],
                }
            ],
        },
        "particle_boundary": {"z_low": "open", "z_high": "open"},
        "reservoir": {"inflow_model": "source_vdf"},
        "particles": {
            "species": [
                plasma_species("electron", -1.602176634e-19, 9.1093837139e-31, 1500),
                # -1 reuses the electron weight.
                plasma_species("proton", 1.602176634e-19, 1.67262192369e-27, -1),
            ]
        },
        "field_boundary": {"mode": "periodic2"},
        "periodic2": {
            "nonzero_mode_backend": "cached_kneq0",
            "zero_mode_policy": "exclude_k0",
            "lower_boundary_model": "e_bottom_zero",
        },
        "output": {
            "write_files": True,
            "dir": "outputs/tutorial",
            "history_stride": 10,
            "write_potential_history": True,
            "checkpoint_stride": 5000,
            "write_mesh_potential": True,
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


def normalize_config_document(config: Mapping[str, Any]) -> dict[str, Any]:
    """Resolve high-level authoring notation and validate the runtime config."""

    return _validate_config_document(config, normalize=True)


def _validate_config_document(
    config: Mapping[str, Any], *, normalize: bool,
    schema: Mapping[str, Any] | None = None, config_path: Path | None = None,
    additional_schema: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    if schema is None:
        schema, _ = load_schema()
    # Reject flat input before schema diagnostics or internal conversion.
    if not is_grouped(config):
        to_runtime_layout(config)
    document = prepare_schema_document(config, schema)
    phase = "authoring" if normalize else "normalized"
    schema_errors_of = schema_validator(schema)
    input_errors = schema_errors_of(document)
    reject_basic_schema_errors(input_errors, phase=phase)
    runtime_document = to_runtime_layout(document)
    if normalize:
        _validate_high_level_fragment(runtime_document, context="config")
        resolved = normalize_high_level_config(runtime_document)
    else:
        resolved = runtime_document
    # Authoring arithmetic can create new invalid values.
    resolved_grouped = prepare_schema_document(to_grouped_layout(resolved), schema)
    resolved_errors = schema_errors_of(resolved_grouped)
    reject_basic_schema_errors(resolved_errors, phase="normalized")
    if input_errors:
        raise ConfigSchemaError(input_errors, phase=phase)
    if resolved_errors:
        raise ConfigSchemaError(resolved_errors, phase="normalized")
    if additional_schema is not None:
        additional_errors_of = schema_validator(additional_schema)
        for candidate, candidate_phase in ((document, phase), (resolved_grouped, "normalized")):
            errors = additional_errors_of(candidate)
            if errors:
                raise ConfigSchemaError(errors, phase=candidate_phase)
    return resolved


def _strip_id_fields(config: dict[str, Any]) -> None:
    """Remove direct-config-only ``id`` keys before writing the final config."""
    _ARRAY_OF_TABLE_PATHS = (
        ("particles", "species"),
        ("mesh", "parts"),
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
    """Resolve high-level spatial notation in an internal runtime dictionary."""

    resolved = copy.deepcopy(dict(config))
    sim = resolved.get("sim")
    if isinstance(sim, Mapping):
        resolved["sim"] = _resolve_box_shorthand(dict(sim))
    domain = resolved.get("domain")
    if isinstance(domain, Mapping):
        resolved["domain"] = _resolve_box_shorthand(dict(domain))
    box_min, box_max = _resolved_box_bounds(resolved)

    particles = resolved.get("particles")
    if isinstance(particles, Mapping):
        species = particles.get("species")
        if isinstance(species, list):
            resolved_species = [
                _resolve_species_high_level(
                    dict(item), species_idx=index, box_min=box_min, box_max=box_max
                )
                for index, item in enumerate(species, start=1)
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
        to_grouped_layout(config),
        header_comments=header_comments,
        top_level_order=TOP_LEVEL_ORDER,
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

"""BEACH config loading and validation."""

from __future__ import annotations

import copy
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from ._toml import load_toml_file, render_toml_document

CONFIG_FILENAME = "beach.toml"
SCHEMA_BASE_URL = "https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas"
BEACH_SCHEMA_URL = f"{SCHEMA_BASE_URL}/beach.schema.json"
TOP_LEVEL_CONFIG_ORDER = (
    "sim",
    "domain",
    "field_boundary",
    "particle_boundary",
    "reservoir",
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
_FACE_KEYS = frozenset(
    {"x_low", "x_high", "y_low", "y_high", "z_low", "z_high"}
)
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
        "softening",
    }
)


class ConfigError(ValueError):
    """Base error for BEACH config handling."""


class ConfigValidationError(ConfigError):
    """Raised when ``beach.toml`` violates known BEACH constraints."""


def default_config() -> dict[str, Any]:
    """Return a small runnable ``beach.toml`` document."""

    return {
        "sim": {
            "dt": 1.0e-7,
            "batch_count": 1,
            "max_step": 10,
            "rng_seed": 12345,
            "field_solver": "fmm",
            "field_periodic_image_layers": 1,
            "field_periodic_far_correction": "none",
        },
        "domain": {
            "box_min": [0.0, 0.0, 0.0],
            "box_max": [1.0, 1.0, 1.0],
            "periodic_axes": ["x", "y"],
        },
        "field_boundary": {"mode": "periodic2"},
        "particle_boundary": {"z_low": "open", "z_high": "open"},
        "particles": {
            "species": [
                {
                    "source_mode": "volume_seed",
                    "q_particle": -1.602176634e-19,
                    "m_particle": 9.10938356e-31,
                    "w_particle": 1.0,
                    "npcls_per_step": 1,
                    "pos_low": [0.5, 0.5, 0.8],
                    "pos_high": [0.5, 0.5, 0.8],
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
                    "size_x": 1.0,
                    "size_y": 1.0,
                    "nx": 4,
                    "ny": 4,
                    "center": [0.5, 0.5, 0.2],
                }
            ],
        },
        "output": {
            "write_files": True,
            "dir": "outputs/latest",
            "history_stride": 1,
        },
    }


def load_config_file(path: str | Path) -> dict[str, Any]:
    """Load, normalize, and validate one direct ``beach.toml`` file."""

    raw = load_toml_file(path)
    return normalize_config_document(raw)


def normalize_config_document(config: Mapping[str, Any]) -> dict[str, Any]:
    """Resolve high-level authoring notation and validate the runtime config."""

    _validate_fragment_structure(
        config,
        context="config",
        allow_meta_keys=False,
    )
    _validate_high_level_fragment(config, context="config")
    normalized = normalize_high_level_config(config)
    validate_runtime_config(normalized)
    return normalized


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


def _resolve_sim_high_level(sim: dict[str, Any]) -> dict[str, Any]:
    box_origin = sim.pop("box_origin", None)
    box_size = sim.pop("box_size", None)
    if (box_origin is None) != (box_size is None):
        raise ConfigError(
            "high-level config error: sim.box_origin and sim.box_size must be specified together."
        )
    if box_origin is not None and box_size is not None:
        origin = _coerce_numeric_sequence(box_origin, length=3, name="sim.box_origin")
        size = _coerce_numeric_sequence(box_size, length=3, name="sim.box_size")
        if any(component <= 0.0 for component in size):
            raise ConfigError(
                "high-level config error: sim.box_size components must be > 0."
            )
        sim["box_min"] = origin
        sim["box_max"] = [origin[i] + size[i] for i in range(3)]
    return sim


def _resolve_domain_high_level(domain: dict[str, Any]) -> dict[str, Any]:
    box_origin = domain.pop("box_origin", None)
    box_size = domain.pop("box_size", None)
    if (box_origin is None) != (box_size is None):
        raise ConfigError(
            "high-level config error: domain.box_origin and domain.box_size "
            "must be specified together."
        )
    if box_origin is not None and box_size is not None:
        origin = _coerce_numeric_sequence(
            box_origin, length=3, name="domain.box_origin"
        )
        size = _coerce_numeric_sequence(box_size, length=3, name="domain.box_size")
        if any(component <= 0.0 for component in size):
            raise ConfigError(
                "high-level config error: domain.box_size components must be > 0."
            )
        domain["box_min"] = origin
        domain["box_max"] = [origin[i] + size[i] for i in range(3)]
    return domain


def _resolve_species_high_level(
    species: dict[str, Any],
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> dict[str, Any]:
    mode = species.pop("inject_region_mode", None)
    uv_low = species.pop("uv_low", None)
    uv_high = species.pop("uv_high", None)
    source_mode = species.get("source_mode", "volume_seed")
    if mode is None:
        if uv_low is not None or uv_high is not None:
            raise ConfigError(
                'high-level config error: uv_low/uv_high require inject_region_mode="face_fraction".'
            )
        return species
    if not isinstance(source_mode, str) or source_mode not in _FACE_SOURCE_MODES:
        raise ConfigError(
            "high-level config error: inject_region_mode is only supported for "
            'source_mode="reservoir_face" or "photo_raycast".'
        )
    if not isinstance(mode, str):
        raise ConfigError(
            "high-level config error: inject_region_mode must be a string."
        )
    if mode == "absolute":
        if uv_low is not None or uv_high is not None:
            raise ConfigError(
                'high-level config error: inject_region_mode="absolute" cannot use uv_low/uv_high.'
            )
        return species
    if mode != "face_fraction":
        raise ConfigError(
            f"high-level config error: unsupported inject_region_mode={mode!r}."
        )
    if "pos_low" in species or "pos_high" in species:
        raise ConfigError(
            "high-level config error: face_fraction injection cannot be combined with pos_low/pos_high."
        )
    if box_min is None or box_max is None:
        raise ConfigError(
            'high-level config error: inject_region_mode="face_fraction" requires domain.box_min/box_max.'
        )
    if uv_low is None or uv_high is None:
        raise ConfigError(
            "high-level config error: face_fraction injection requires uv_low and uv_high."
        )
    uv_low_vec = _coerce_numeric_sequence(uv_low, length=2, name="uv_low")
    uv_high_vec = _coerce_numeric_sequence(uv_high, length=2, name="uv_high")
    if any(value < 0.0 or value > 1.0 for value in [*uv_low_vec, *uv_high_vec]):
        raise ConfigError(
            "high-level config error: uv_low/uv_high must be inside [0, 1]."
        )
    if any(uv_low_vec[i] > uv_high_vec[i] for i in range(2)):
        raise ConfigError(
            "high-level config error: uv_low must be <= uv_high component-wise."
        )
    inject_face = species.get("inject_face")
    if not isinstance(inject_face, str) or inject_face == "":
        raise ConfigError(
            "high-level config error: face_fraction injection requires inject_face."
        )
    pos_low, pos_high = _resolve_face_fraction_region(
        inject_face=inject_face,
        uv_low=uv_low_vec,
        uv_high=uv_high_vec,
        box_min=box_min,
        box_max=box_max,
    )
    species["pos_low"] = pos_low
    species["pos_high"] = pos_high
    return species


def _resolve_mesh_high_level(
    mesh: dict[str, Any],
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> dict[str, Any]:
    groups_raw = mesh.pop("groups", {})
    group_definitions = _normalize_group_definitions(groups_raw)
    templates = mesh.get("templates")
    if not isinstance(templates, list):
        return mesh

    resolved_templates: list[dict[str, Any]] = []
    for template in templates:
        if not isinstance(template, Mapping):
            raise ConfigError(
                "high-level config error: mesh.templates entries must be tables."
            )
        template_dict = dict(template)
        if "group" in template_dict:
            resolved_templates.append(
                _resolve_grouped_template(
                    template_dict,
                    group_definitions=group_definitions,
                    box_min=box_min,
                    box_max=box_max,
                )
            )
        else:
            resolved_templates.append(
                _resolve_direct_template(
                    template_dict,
                    box_min=box_min,
                    box_max=box_max,
                )
            )
    mesh["templates"] = resolved_templates
    return mesh


def _resolve_direct_template(
    template: dict[str, Any],
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> dict[str, Any]:
    placement_mode = template.pop("placement_mode", "absolute")
    anchor = template.pop("anchor", None)
    offset = template.pop("offset", None)
    offset_frac = template.pop("offset_frac", None)

    if placement_mode == "absolute":
        if anchor is not None or offset is not None or offset_frac is not None:
            raise ConfigError(
                'high-level config error: placement_mode="absolute" cannot use anchor/offset/offset_frac.'
            )
        if "center" not in template:
            raise ConfigError(
                "high-level config error: direct mesh template requires center."
            )
    elif placement_mode == "box_anchor":
        if "center" in template:
            raise ConfigError(
                'high-level config error: placement_mode="box_anchor" cannot be combined with center.'
            )
        template["center"] = _resolve_anchor_position(
            anchor=anchor,
            offset=offset,
            offset_frac=offset_frac,
            box_min=box_min,
            box_max=box_max,
        )
    else:
        raise ConfigError(
            f"high-level config error: unsupported placement_mode={placement_mode!r}."
        )

    _resolve_template_size_high_level(template, box_min=box_min, box_max=box_max)
    return template


def _resolve_grouped_template(
    template: dict[str, Any],
    *,
    group_definitions: Mapping[str, Mapping[str, Any]],
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> dict[str, Any]:
    group_name = template.pop("group")
    if not isinstance(group_name, str) or group_name == "":
        raise ConfigError(
            "high-level config error: mesh template group must be a non-empty string."
        )
    if group_name not in group_definitions:
        raise ConfigError(
            f'high-level config error: mesh template references undefined group "{group_name}".'
        )
    for forbidden_key in (
        "center",
        "placement_mode",
        "anchor",
        "offset",
        "offset_frac",
        "size_mode",
        "size_frac",
    ):
        if forbidden_key in template:
            raise ConfigError(
                f'high-level config error: grouped mesh template "{group_name}" cannot define {forbidden_key}.'
            )
    center_local = template.pop("center_local", None)
    if center_local is None:
        raise ConfigError(
            f'high-level config error: grouped mesh template "{group_name}" requires center_local.'
        )
    center_local_vec = _coerce_numeric_sequence(
        center_local, length=3, name="center_local"
    )
    group_origin = _resolve_group_origin(
        dict(group_definitions[group_name]),
        box_min=box_min,
        box_max=box_max,
    )
    group_scale = _resolve_group_scale(
        dict(group_definitions[group_name]),
        box_min=box_min,
        box_max=box_max,
    )
    template["center"] = [
        group_origin[i] + group_scale * center_local_vec[i] for i in range(3)
    ]
    _scale_template_lengths(template, group_scale)
    return template


def _resolve_template_size_high_level(
    template: dict[str, Any],
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> None:
    size_mode = template.pop("size_mode", "absolute")
    size_frac = template.pop("size_frac", None)
    if size_mode == "absolute":
        if size_frac is not None:
            raise ConfigError(
                'high-level config error: size_frac requires size_mode="box_fraction".'
            )
        return
    if size_mode != "box_fraction":
        raise ConfigError(
            f"high-level config error: unsupported size_mode={size_mode!r}."
        )
    if size_frac is None:
        raise ConfigError(
            'high-level config error: size_mode="box_fraction" requires size_frac.'
        )
    box_size = _require_box_size(box_min=box_min, box_max=box_max, context="size_frac")
    kind = template.get("kind", "plane")
    if not isinstance(kind, str):
        raise ConfigError(
            "high-level config error: mesh template kind must be a string."
        )
    if kind in {"plane", "plane_hole", "plate_hole"}:
        frac = _coerce_numeric_sequence(size_frac, length=2, name="size_frac")
        template["size_x"] = frac[0] * box_size[0]
        template["size_y"] = frac[1] * box_size[1]
        return
    if kind == "box":
        frac = _coerce_numeric_sequence(size_frac, length=3, name="size_frac")
        template["size"] = [frac[i] * box_size[i] for i in range(3)]
        return
    if kind == "sphere":
        frac = _coerce_numeric_scalar(size_frac, name="size_frac")
        template["radius"] = frac * min(box_size)
        return
    if kind == "cylinder":
        frac = _coerce_numeric_sequence(size_frac, length=2, name="size_frac")
        template["radius"] = frac[0] * min(box_size[0], box_size[1])
        template["height"] = frac[1] * box_size[2]
        return
    raise ConfigError(
        f'high-level config error: size_mode="box_fraction" is not supported for kind={kind!r}.'
    )


def _resolve_face_fraction_region(
    *,
    inject_face: str,
    uv_low: Sequence[float],
    uv_high: Sequence[float],
    box_min: Sequence[float],
    box_max: Sequence[float],
) -> tuple[list[float], list[float]]:
    pos_low = [float(component) for component in box_min]
    pos_high = [float(component) for component in box_max]
    if inject_face == "x_low":
        boundary = float(box_min[0])
        ranges = ((1, box_min[1], box_max[1]), (2, box_min[2], box_max[2]))
        axis = 0
    elif inject_face == "x_high":
        boundary = float(box_max[0])
        ranges = ((1, box_min[1], box_max[1]), (2, box_min[2], box_max[2]))
        axis = 0
    elif inject_face == "y_low":
        boundary = float(box_min[1])
        ranges = ((0, box_min[0], box_max[0]), (2, box_min[2], box_max[2]))
        axis = 1
    elif inject_face == "y_high":
        boundary = float(box_max[1])
        ranges = ((0, box_min[0], box_max[0]), (2, box_min[2], box_max[2]))
        axis = 1
    elif inject_face == "z_low":
        boundary = float(box_min[2])
        ranges = ((0, box_min[0], box_max[0]), (1, box_min[1], box_max[1]))
        axis = 2
    elif inject_face == "z_high":
        boundary = float(box_max[2])
        ranges = ((0, box_min[0], box_max[0]), (1, box_min[1], box_max[1]))
        axis = 2
    else:
        raise ConfigError(
            f"high-level config error: invalid inject_face={inject_face!r}."
        )

    pos_low[axis] = boundary
    pos_high[axis] = boundary
    for uv_index, (coordinate_index, low_bound, high_bound) in enumerate(ranges):
        span = float(high_bound) - float(low_bound)
        pos_low[coordinate_index] = float(low_bound) + uv_low[uv_index] * span
        pos_high[coordinate_index] = float(low_bound) + uv_high[uv_index] * span
    return pos_low, pos_high


def _normalize_group_definitions(groups_raw: object) -> dict[str, dict[str, Any]]:
    if groups_raw in ({}, None):
        return {}
    if not isinstance(groups_raw, Mapping):
        raise ConfigError(
            "high-level config error: mesh.groups must be a table of group definitions."
        )
    normalized: dict[str, dict[str, Any]] = {}
    for name, value in groups_raw.items():
        if not isinstance(value, Mapping):
            raise ConfigError(
                f"high-level config error: mesh.groups.{name} must be a table."
            )
        normalized[name] = dict(value)
    return normalized


def _resolve_group_origin(
    group: dict[str, Any],
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> list[float]:
    placement_mode = group.pop("placement_mode", "absolute")
    anchor = group.pop("anchor", None)
    offset = group.pop("offset", None)
    offset_frac = group.pop("offset_frac", None)
    group.pop("scale", None)
    group.pop("scale_from", None)
    group.pop("scale_factor", None)
    if placement_mode == "absolute":
        if anchor is not None:
            raise ConfigError(
                'high-level config error: group placement_mode="absolute" cannot use anchor.'
            )
        base = [0.0, 0.0, 0.0]
        offset_vec = _resolve_offset_vector(
            offset=offset,
            offset_frac=offset_frac,
            box_min=box_min,
            box_max=box_max,
            context="group absolute placement",
            allow_zero=True,
        )
        return [base[i] + offset_vec[i] for i in range(3)]
    if placement_mode != "box_anchor":
        raise ConfigError(
            f"high-level config error: unsupported group placement_mode={placement_mode!r}."
        )
    return _resolve_anchor_position(
        anchor=anchor,
        offset=offset,
        offset_frac=offset_frac,
        box_min=box_min,
        box_max=box_max,
    )


def _resolve_group_scale(
    group: dict[str, Any],
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> float:
    scale = group.get("scale")
    scale_from = group.get("scale_from")
    scale_factor = group.get("scale_factor")
    if scale is not None and (scale_from is not None or scale_factor is not None):
        raise ConfigError(
            "high-level config error: group scale cannot combine scale with scale_from/scale_factor."
        )
    if scale is not None:
        scale_value = _coerce_numeric_scalar(scale, name="scale")
        if scale_value <= 0.0:
            raise ConfigError("high-level config error: group scale must be > 0.")
        return scale_value
    if scale_from is None and scale_factor is None:
        return 1.0
    if scale_from is None or scale_factor is None:
        raise ConfigError(
            "high-level config error: group scale_from and scale_factor must be specified together."
        )
    if not isinstance(scale_from, str):
        raise ConfigError("high-level config error: scale_from must be a string.")
    factor = _coerce_numeric_scalar(scale_factor, name="scale_factor")
    if factor <= 0.0:
        raise ConfigError("high-level config error: scale_factor must be > 0.")
    reference = _resolve_scale_reference(scale_from, box_min=box_min, box_max=box_max)
    return factor * reference


def _resolve_scale_reference(
    scale_from: str,
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> float:
    box_size = _require_box_size(box_min=box_min, box_max=box_max, context="scale_from")
    refs = {
        "box_x": box_size[0],
        "box_y": box_size[1],
        "box_z": box_size[2],
        "box_min_xy": min(box_size[0], box_size[1]),
        "box_max_xy": max(box_size[0], box_size[1]),
        "box_min_xyz": min(box_size),
        "box_max_xyz": max(box_size),
    }
    if scale_from not in refs:
        raise ConfigError(
            f"high-level config error: unsupported scale_from={scale_from!r}."
        )
    return refs[scale_from]


def _resolve_anchor_position(
    *,
    anchor: object,
    offset: object,
    offset_frac: object,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> list[float]:
    if not isinstance(anchor, str) or anchor == "":
        raise ConfigError(
            "high-level config error: box_anchor placement requires anchor."
        )
    base = _resolve_anchor(anchor, box_min=box_min, box_max=box_max)
    offset_vec = _resolve_offset_vector(
        offset=offset,
        offset_frac=offset_frac,
        box_min=box_min,
        box_max=box_max,
        context=f"anchor {anchor}",
        allow_zero=True,
    )
    return [base[i] + offset_vec[i] for i in range(3)]


def _resolve_anchor(
    anchor: str,
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> list[float]:
    if box_min is None or box_max is None:
        raise ConfigError(
            "high-level config error: box_anchor placement requires domain.box_min/box_max."
        )
    center = [(box_min[i] + box_max[i]) * 0.5 for i in range(3)]
    anchors = {
        "box_center": center,
        "x_low_face_center": [box_min[0], center[1], center[2]],
        "x_high_face_center": [box_max[0], center[1], center[2]],
        "y_low_face_center": [center[0], box_min[1], center[2]],
        "y_high_face_center": [center[0], box_max[1], center[2]],
        "z_low_face_center": [center[0], center[1], box_min[2]],
        "z_high_face_center": [center[0], center[1], box_max[2]],
    }
    if anchor not in anchors:
        raise ConfigError(f"high-level config error: unsupported anchor={anchor!r}.")
    return list(anchors[anchor])


def _resolve_offset_vector(
    *,
    offset: object,
    offset_frac: object,
    box_min: list[float] | None,
    box_max: list[float] | None,
    context: str,
    allow_zero: bool,
) -> list[float]:
    if offset is not None and offset_frac is not None:
        raise ConfigError(
            f"high-level config error: {context} cannot combine offset and offset_frac."
        )
    if offset is None and offset_frac is None:
        return [0.0, 0.0, 0.0] if allow_zero else []
    if offset is not None:
        return _coerce_numeric_sequence(offset, length=3, name="offset")
    frac = _coerce_numeric_sequence(offset_frac, length=3, name="offset_frac")
    box_size = _require_box_size(
        box_min=box_min, box_max=box_max, context="offset_frac"
    )
    return [frac[i] * box_size[i] for i in range(3)]


def _scale_template_lengths(template: dict[str, Any], scale: float) -> None:
    for key in ("size_x", "size_y", "radius", "inner_radius", "height"):
        if key in template:
            template[key] = scale * float(template[key])
    if "size" in template:
        template["size"] = [scale * float(value) for value in template["size"]]


def _resolved_box_bounds(
    config: Mapping[str, Any],
) -> tuple[list[float] | None, list[float] | None]:
    domain = config.get("domain")
    if isinstance(domain, Mapping):
        box_min = domain.get("box_min")
        box_max = domain.get("box_max")
        name_prefix = "domain"
    else:
        sim = config.get("sim")
        if not isinstance(sim, Mapping):
            return None, None
        box_min = sim.get("box_min")
        box_max = sim.get("box_max")
        name_prefix = "sim"
    if box_min is None or box_max is None:
        return None, None
    return (
        _coerce_numeric_sequence(box_min, length=3, name=f"{name_prefix}.box_min"),
        _coerce_numeric_sequence(box_max, length=3, name=f"{name_prefix}.box_max"),
    )


def _require_box_size(
    *,
    box_min: list[float] | None,
    box_max: list[float] | None,
    context: str,
) -> list[float]:
    if box_min is None or box_max is None:
        raise ConfigError(
            f"high-level config error: {context} requires domain.box_min/box_max."
        )
    box_size = [box_max[i] - box_min[i] for i in range(3)]
    if any(component <= 0.0 for component in box_size):
        raise ConfigError(
            f"high-level config error: {context} requires positive box dimensions."
        )
    return box_size


def _coerce_numeric_sequence(value: object, *, length: int, name: str) -> list[float]:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != length
    ):
        raise ConfigError(
            f"high-level config error: {name} must be a {length}-element numeric array."
        )
    out: list[float] = []
    for index in range(length):
        item = value[index]
        if not isinstance(item, (int, float)) or isinstance(item, bool):
            raise ConfigError(
                f"high-level config error: {name} must contain only numeric values."
            )
        numeric = float(item)
        if not math.isfinite(numeric):
            raise ConfigError(
                f"high-level config error: {name} must contain only finite values."
            )
        out.append(numeric)
    return out


def _coerce_numeric_scalar(value: object, *, name: str) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ConfigError(f"high-level config error: {name} must be numeric.")
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ConfigError(f"high-level config error: {name} must be finite.")
    return numeric


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

    resolved_batch_duration = _resolve_batch_duration(sim)
    use_box = domain is not None
    adaptive_nonzero_mode_limit = 0.0
    if isinstance(periodic2_config, Mapping):
        raw_adaptive_limit = periodic2_config.get(
            "max_nonzero_mode_potential_step", 0.0
        )
        if (
            not isinstance(raw_adaptive_limit, (int, float))
            or isinstance(raw_adaptive_limit, bool)
            or not math.isfinite(float(raw_adaptive_limit))
            or float(raw_adaptive_limit) < 0.0
        ):
            raise ConfigValidationError(
                "BEACH constraint error: "
                "periodic2.max_nonzero_mode_potential_step must be finite and >= 0."
            )
        adaptive_nonzero_mode_limit = float(raw_adaptive_limit)

    field_bc_mode = (
        field_boundary.get("mode", "free")
        if field_boundary is not None
        else "free"
    )
    field_solver = sim.get("field_solver", "auto")
    if field_solver not in {"direct", "treecode", "fmm", "auto"}:
        raise ConfigValidationError(
            'BEACH constraint error: sim.field_solver must be "direct", "treecode", '
            '"fmm", or "auto".'
        )
    phi_infty = reservoir.get("phi_infty", 0.0) if reservoir is not None else 0.0
    if (
        not isinstance(phi_infty, (int, float))
        or isinstance(phi_infty, bool)
        or not math.isfinite(phi_infty)
    ):
        raise ConfigValidationError(
            "BEACH constraint error: reservoir.phi_infty must be finite."
        )
    if field_bc_mode not in {"free", "periodic2"}:
        raise ConfigValidationError(
            'BEACH constraint error: field_boundary.mode must be "free" or "periodic2".'
        )
    if field_bc_mode == "periodic2":
        supported_lower_boundaries = {"e_bottom_zero", "symmetric_vacuum"}
        if isinstance(periodic2_config, Mapping):
            lower_boundary_model = periodic2_config.get("lower_boundary_model")
            if (
                lower_boundary_model is not None
                and lower_boundary_model not in supported_lower_boundaries
            ):
                raise ConfigValidationError(
                    "BEACH constraint error: periodic2.lower_boundary_model must be "
                    '"e_bottom_zero" or "symmetric_vacuum".'
                )
        split_reference = (
            field_solver == "direct"
            and isinstance(periodic2_config, Mapping)
            and periodic2_config.get("nonzero_mode_backend")
            == "panel_spectral_reference"
            and periodic2_config.get("zero_mode_policy") == "exclude_k0"
            and periodic2_config.get("lower_boundary_model")
            in supported_lower_boundaries
        )
        if field_solver != "fmm" and not split_reference:
            raise ConfigValidationError(
                'BEACH constraint error: field_boundary.mode="periodic2" requires field_solver="fmm" '
                "or the direct panel_spectral_reference split model."
            )
        if not use_box:
            raise ConfigValidationError(
                'BEACH constraint error: field_boundary.mode="periodic2" requires [domain].'
            )
        periodic_axes = set(domain.get("periodic_axes", [])) if domain else set()
        if periodic_axes != {"x", "y"}:
            raise ConfigValidationError(
                'BEACH constraint error: field_boundary.mode="periodic2" requires '
                'domain.periodic_axes=["x", "y"].'
            )
    if adaptive_nonzero_mode_limit > 0.0:
        if (
            not isinstance(periodic2_config, Mapping)
            or periodic2_config.get("nonzero_mode_backend") != "cached_kneq0"
        ):
            raise ConfigValidationError(
                "BEACH constraint error: "
                "periodic2.max_nonzero_mode_potential_step requires "
                'nonzero_mode_backend="cached_kneq0".'
            )
        if not math.isfinite(resolved_batch_duration) or resolved_batch_duration <= 0.0:
            raise ConfigValidationError(
                "BEACH constraint error: "
                "periodic2.max_nonzero_mode_potential_step requires a finite positive "
                "sim.batch_duration or resolved batch_duration_step."
            )
    if field_bc_mode != "free" and _mesh_has_surface_model(mesh, "conductor"):
        raise ConfigValidationError(
            'BEACH constraint error: surface_model="conductor" currently requires '
            'field_boundary.mode="free".'
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
    periodic_axes = set(domain.get("periodic_axes", [])) if domain is not None else set()
    global_particle_boundary = particle_boundary or {}

    uses_face_sources = False
    has_volume_seed = False
    total_npcls_per_step = 0
    for index, item in enumerate(species, start=1):
        species_table = dict(item)
        source_mode = species_table.get("source_mode", "volume_seed")
        if not isinstance(source_mode, str):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] source_mode must be a string."
            )
        effective_particle_boundary = _validate_species_particle_boundary(
            species_table,
            index=index,
            periodic_axes=periodic_axes,
            global_boundary=global_particle_boundary,
        )
        boundary_inflow_faces = _validate_species_boundary_inflow(
            species_table,
            index=index,
            periodic_axes=periodic_axes,
            effective_boundary=effective_particle_boundary,
        )
        has_reservoir_injection = (
            source_mode in _RESERVOIR_SOURCE_MODES
            or bool(boundary_inflow_faces)
        )
        if has_reservoir_injection:
            uses_face_sources = True
            _validate_reservoir_injection_common(
                species_table,
                index=index,
                use_box=use_box,
                batch_duration=resolved_batch_duration,
            )
            _validate_reservoir_physics(
                species_table,
                index=index,
                source_mode=source_mode,
            )
        surface_charge_closure = species_table.get(
            "surface_charge_closure", "explicit"
        )
        if surface_charge_closure not in {"explicit", "neutral_return"}:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}]."
                'surface_charge_closure must be "explicit" or "neutral_return".'
            )
        if surface_charge_closure == "neutral_return":
            inject_face = species_table.get("inject_face")
            q_particle = species_table.get("q_particle", -1.602176634e-19)
            if (
                source_mode != "photo_raycast"
                or not isinstance(q_particle, (int, float))
                or isinstance(q_particle, bool)
                or float(q_particle) >= 0.0
                or species_table.get("deposit_opposite_charge_on_emit") is not True
                or not isinstance(inject_face, str)
                or effective_particle_boundary.get(inject_face)
                not in {"reflect", "redistributed_reflect"}
            ):
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}] "
                    "neutral_return requires a negative photo_raycast species, "
                    "deposit_opposite_charge_on_emit=true, and a reflecting "
                    "action on inject_face."
                )
        if "photo_escape_model" in species_table:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].photo_escape_model "
                "was removed; track emitted photoelectrons and use the neutral-return "
                "surface charge closure instead."
            )
        if "temperature_k" in species_table and "temperature_ev" in species_table:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] cannot define both "
                "temperature_k and temperature_ev."
            )
        velocity_distribution = (
            str(species_table.get("velocity_distribution", "maxwellian"))
            .strip()
            .lower()
        )
        velocity_grid_pdf_kind = (
            str(species_table.get("velocity_grid_pdf_kind", "phase_space"))
            .strip()
            .lower()
        )
        velocity_grid_sampling = (
            str(species_table.get("velocity_grid_sampling", "auto")).strip().lower()
        )
        if velocity_distribution not in {"maxwellian", "grid"}:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] has unsupported "
                f"velocity_distribution={velocity_distribution!r}."
            )
        if velocity_grid_pdf_kind not in {"phase_space", "flux_weighted"}:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] has unsupported "
                f"velocity_grid_pdf_kind={velocity_grid_pdf_kind!r}."
            )
        if velocity_grid_sampling not in {"auto", "rectilinear", "discrete"}:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] has unsupported "
                f"velocity_grid_sampling={velocity_grid_sampling!r}."
            )

        if source_mode == "volume_seed":
            has_volume_seed = True
            if not has_reservoir_injection:
                _validate_velocity_grid_forbidden(
                    species_table,
                    index=index,
                    source_mode=source_mode,
                )
            npcls_per_step = species_table.get("npcls_per_step", 0)
            if not isinstance(npcls_per_step, int):
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}].npcls_per_step "
                    "must be an integer."
                )
            total_npcls_per_step += npcls_per_step
            if (
                not has_reservoir_injection
                and "target_macro_particles_per_batch" in species_table
            ):
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}] uses "
                    'source_mode="volume_seed" and cannot define '
                    "target_macro_particles_per_batch."
                )
            continue

        if source_mode == "reservoir_face":
            _validate_face_source_common(
                species_table,
                index=index,
                source_mode=source_mode,
                use_box=use_box,
                batch_duration=resolved_batch_duration,
                box_min=box_min,
                box_max=box_max,
            )
            continue

        if source_mode == "plane_source":
            _validate_plane_source_geometry(
                species_table,
                index=index,
                use_box=use_box,
                box_min=box_min,
                box_max=box_max,
            )
            continue

        if source_mode == "photo_raycast":
            uses_face_sources = True
            _validate_velocity_grid_forbidden(
                species_table,
                index=index,
                source_mode=source_mode,
            )
            _validate_face_source_common(
                species_table,
                index=index,
                source_mode=source_mode,
                use_box=use_box,
                batch_duration=resolved_batch_duration,
                box_min=box_min,
                box_max=box_max,
            )
            current_density = species_table.get("emit_current_density_a_m2", 0.0)
            rays_per_batch = species_table.get("rays_per_batch", 0)
            if (
                not isinstance(current_density, (int, float))
                or float(current_density) <= 0.0
            ):
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}] uses "
                    'source_mode="photo_raycast" and requires emit_current_density_a_m2 > 0.'
                )
            if not isinstance(rays_per_batch, int) or rays_per_batch <= 0:
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}] uses "
                    'source_mode="photo_raycast" and requires rays_per_batch > 0.'
                )
            forbidden = ["npcls_per_step"]
            if not has_reservoir_injection:
                forbidden.extend(
                    (
                        "number_density_cm3",
                        "number_density_m3",
                        "w_particle",
                        "target_macro_particles_per_batch",
                    )
                )
            for key in forbidden:
                if key in species_table:
                    raise ConfigValidationError(
                        f"BEACH constraint error: particles.species[{index}] uses "
                        f'source_mode="photo_raycast" and cannot define {key}.'
                    )
            continue

        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] has unsupported "
            f"source_mode={source_mode!r}."
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


def _validate_runtime_external_e_field(sim: Mapping[str, Any]) -> None:
    has_vector = "e0" in sim
    has_abs = "e0_abs" in sim
    has_phi_xy = "e0_phi_xy_deg" in sim
    has_phi_z = "e0_phi_z_deg" in sim
    if has_vector and (has_abs or has_phi_xy or has_phi_z):
        raise ConfigValidationError(
            "BEACH constraint error: sim.e0 cannot be combined with "
            "sim.e0_abs/e0_phi_xy_deg/e0_phi_z_deg."
        )
    if has_vector:
        e0 = sim.get("e0")
        if (
            not isinstance(e0, Sequence)
            or isinstance(e0, (str, bytes))
            or len(e0) != 3
            or not all(
                isinstance(v, (int, float)) and math.isfinite(float(v)) for v in e0
            )
        ):
            raise ConfigValidationError(
                "BEACH constraint error: sim.e0 must contain 3 finite values."
            )
        return
    if (has_phi_xy or has_phi_z) and not has_abs:
        raise ConfigValidationError(
            "BEACH constraint error: sim.e0_phi_xy_deg/e0_phi_z_deg require sim.e0_abs."
        )
    if has_abs:
        e0_abs = sim.get("e0_abs")
        if (
            not isinstance(e0_abs, (int, float))
            or not math.isfinite(float(e0_abs))
            or float(e0_abs) < 0.0
        ):
            raise ConfigValidationError(
                "BEACH constraint error: sim.e0_abs must be finite and >= 0."
            )
        for key in ("e0_phi_xy_deg", "e0_phi_z_deg"):
            value = sim.get(key, 0.0)
            if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                raise ConfigValidationError(
                    f"BEACH constraint error: sim.{key} must be finite."
                )


def _validate_grid_flux_keys(species_table: Mapping[str, Any], *, index: int) -> None:
    has_flux = "particle_flux_m2_s" in species_table
    has_current = "current_density_a_m2" in species_table
    if has_flux and has_current:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] cannot define both "
            "particle_flux_m2_s and current_density_a_m2."
        )
    if not has_flux and not has_current:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            'velocity_distribution="grid" and requires particle_flux_m2_s '
            "or current_density_a_m2."
        )
    if has_flux:
        value = species_table["particle_flux_m2_s"]
        if (
            not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            or float(value) <= 0.0
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].particle_flux_m2_s "
                "must be finite and > 0."
            )
    if has_current:
        value = species_table["current_density_a_m2"]
        if (
            not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            or float(value) == 0.0
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].current_density_a_m2 "
                "must be finite and non-zero."
            )
        q_particle = species_table.get("q_particle", -1.602176634e-19)
        if (
            not isinstance(q_particle, (int, float))
            or not math.isfinite(float(q_particle))
            or float(q_particle) == 0.0
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].q_particle "
                "must be finite and non-zero for current_density_a_m2."
            )


def _validate_velocity_grid_forbidden(
    species_table: Mapping[str, Any],
    *,
    index: int,
    source_mode: str,
) -> None:
    if (
        str(species_table.get("velocity_distribution", "maxwellian")).strip().lower()
        != "maxwellian"
    ):
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            f'source_mode="{source_mode}" and cannot define velocity_distribution="grid".'
        )
    if (
        str(species_table.get("velocity_grid_sampling", "auto")).strip().lower()
        != "auto"
    ):
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            f'source_mode="{source_mode}" and cannot define velocity_grid_sampling.'
        )
    for key in ("velocity_grid_path", "particle_flux_m2_s", "current_density_a_m2"):
        if key in species_table:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] uses "
                f'source_mode="{source_mode}" and cannot define {key}.'
            )


def _validate_face_source_common(
    species_table: Mapping[str, Any],
    *,
    index: int,
    source_mode: str,
    use_box: bool,
    batch_duration: float,
    box_min: list[float] | None,
    box_max: list[float] | None,
) -> None:
    if not use_box:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            f'source_mode="{source_mode}" and requires [domain].'
        )
    if batch_duration <= 0.0:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            f'source_mode="{source_mode}" and requires batch_duration > 0.'
        )
    inject_face = species_table.get("inject_face")
    if not isinstance(inject_face, str) or inject_face == "":
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            f'source_mode="{source_mode}" and requires inject_face.'
        )
    if box_min is not None and box_max is not None:
        _validate_face_bounds(
            species_table,
            index=index,
            inject_face=inject_face,
            box_min=box_min,
            box_max=box_max,
        )


def _validate_face_bounds(
    species_table: Mapping[str, Any],
    *,
    index: int,
    inject_face: str,
    box_min: Sequence[float],
    box_max: Sequence[float],
) -> None:
    pos_low = _maybe_vec3(
        species_table.get("pos_low"),
        name=f"particles.species[{index}].pos_low",
    )
    pos_high = _maybe_vec3(
        species_table.get("pos_high"),
        name=f"particles.species[{index}].pos_high",
    )
    if pos_low is None or pos_high is None:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] must define pos_low "
            "and pos_high on the inject_face."
        )

    axis_by_face = {
        "x_low": (0, float(box_min[0])),
        "x_high": (0, float(box_max[0])),
        "y_low": (1, float(box_min[1])),
        "y_high": (1, float(box_max[1])),
        "z_low": (2, float(box_min[2])),
        "z_high": (2, float(box_max[2])),
    }
    if inject_face not in axis_by_face:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] has invalid inject_face={inject_face!r}."
        )
    axis, boundary = axis_by_face[inject_face]
    if pos_low[axis] != boundary or pos_high[axis] != boundary:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] pos_low/pos_high must "
            f"lie on inject_face={inject_face!r}."
        )
    for other_axis in range(3):
        if other_axis == axis:
            continue
        if pos_low[other_axis] > pos_high[other_axis]:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] pos_low must be <= "
                "pos_high along the inject-face coordinates."
            )
        low_bound = float(box_min[other_axis])
        high_bound = float(box_max[other_axis])
        if not (low_bound <= pos_low[other_axis] <= high_bound):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] pos_low is outside the box."
            )
        if not (low_bound <= pos_high[other_axis] <= high_bound):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] pos_high is outside the box."
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
        raise ConfigValidationError(
            f"BEACH constraint error: [{key}] must be a table."
        )
    return dict(value)


def _validate_runtime_boundary_tables(
    *,
    domain: Mapping[str, Any] | None,
    field_boundary: Mapping[str, Any] | None,
    particle_boundary: Mapping[str, Any] | None,
    reservoir: Mapping[str, Any] | None,
) -> None:
    periodic_axes: set[str] = set()
    if domain is not None:
        unknown = set(domain) - {"box_min", "box_max", "periodic_axes"}
        if unknown:
            raise ConfigValidationError(
                "BEACH constraint error: unsupported domain key(s): "
                + ", ".join(sorted(unknown))
                + "."
            )
        if "box_min" not in domain or "box_max" not in domain:
            raise ConfigValidationError(
                "BEACH constraint error: [domain] requires box_min and box_max."
            )
        box_min = _maybe_vec3(domain.get("box_min"), name="domain.box_min")
        box_max = _maybe_vec3(domain.get("box_max"), name="domain.box_max")
        assert box_min is not None and box_max is not None
        if any(box_max[i] <= box_min[i] for i in range(3)):
            raise ConfigValidationError(
                "BEACH constraint error: domain.box_max must be greater than "
                "domain.box_min on every axis."
            )
        raw_axes = domain.get("periodic_axes", [])
        if (
            not isinstance(raw_axes, list)
            or not all(isinstance(axis, str) for axis in raw_axes)
            or len(raw_axes) != len(set(raw_axes))
            or not set(raw_axes) <= {"x", "y", "z"}
        ):
            raise ConfigValidationError(
                "BEACH constraint error: domain.periodic_axes must contain unique "
                'axis names from "x", "y", and "z".'
            )
        periodic_axes = set(raw_axes)

    if field_boundary is not None:
        unknown = set(field_boundary) - {"mode"}
        if unknown:
            raise ConfigValidationError(
                "BEACH constraint error: unsupported field_boundary key(s): "
                + ", ".join(sorted(unknown))
                + "."
            )
        if field_boundary.get("mode", "free") not in {"free", "periodic2"}:
            raise ConfigValidationError(
                'BEACH constraint error: field_boundary.mode must be "free" or '
                '"periodic2".'
            )

    if particle_boundary is not None:
        face_keys = {
            "x_low",
            "x_high",
            "y_low",
            "y_high",
            "z_low",
            "z_high",
        }
        unknown = set(particle_boundary) - face_keys - {"ordinary_open_model"}
        if unknown:
            raise ConfigValidationError(
                "BEACH constraint error: unsupported particle_boundary key(s): "
                + ", ".join(sorted(unknown))
                + "."
            )
        for face in face_keys & set(particle_boundary):
            if particle_boundary[face] not in {
                "open",
                "reflect",
                "redistributed_reflect",
            }:
                raise ConfigValidationError(
                    f"BEACH constraint error: particle_boundary.{face} must be "
                    '"open", "reflect", or "redistributed_reflect".'
                )
            if face[0] in periodic_axes:
                raise ConfigValidationError(
                    f"BEACH constraint error: particle_boundary.{face} cannot "
                    "override a periodic domain face."
                )
        if particle_boundary.get("ordinary_open_model", "escape") not in {
            "escape",
            "potential_barrier",
        }:
            raise ConfigValidationError(
                "BEACH constraint error: particle_boundary.ordinary_open_model must "
                'be "escape" or "potential_barrier".'
            )

    if reservoir is not None:
        unknown = set(reservoir) - {
            "inflow_model",
            "phi_infty",
            "face_potential_grid_n",
        }
        if unknown:
            raise ConfigValidationError(
                "BEACH constraint error: unsupported reservoir key(s): "
                + ", ".join(sorted(unknown))
                + "."
            )
        if reservoir.get("inflow_model", "source_vdf") not in {
            "source_vdf",
            "infinity_barrier",
        }:
            raise ConfigValidationError(
                "BEACH constraint error: reservoir.inflow_model must be "
                '"source_vdf" or "infinity_barrier".'
            )
        grid_n = reservoir.get("face_potential_grid_n", 3)
        if not isinstance(grid_n, int) or isinstance(grid_n, bool) or grid_n < 1:
            raise ConfigValidationError(
                "BEACH constraint error: reservoir.face_potential_grid_n must be >= 1."
            )


def _validate_species_particle_boundary(
    species: Mapping[str, Any],
    *,
    index: int,
    periodic_axes: set[str],
    global_boundary: Mapping[str, Any],
) -> dict[str, str]:
    face_keys = {
        "x_low",
        "x_high",
        "y_low",
        "y_high",
        "z_low",
        "z_high",
    }
    raw_boundary = species.get("boundary", {})
    if not isinstance(raw_boundary, Mapping):
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}].boundary must be a table."
        )
    unknown = set(raw_boundary) - face_keys
    if unknown:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}].boundary has "
            "unsupported key(s): "
            + ", ".join(sorted(unknown))
            + "."
        )

    effective: dict[str, str] = {}
    for face in face_keys:
        axis = face[0]
        action = raw_boundary.get(face, "inherit")
        if action not in {
            "inherit",
            "open",
            "reflect",
            "redistributed_reflect",
        }:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].boundary."
                f'{face} must be "inherit", "open", "reflect", or '
                '"redistributed_reflect".'
            )
        if axis in periodic_axes:
            if action != "inherit":
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}].boundary."
                    f"{face} cannot override a periodic domain face."
                )
            effective[face] = "periodic"
        elif action == "inherit":
            effective[face] = str(global_boundary.get(face, "open"))
        else:
            effective[face] = str(action)
    return effective


def _validate_species_boundary_inflow(
    species: Mapping[str, Any],
    *,
    index: int,
    periodic_axes: set[str],
    effective_boundary: Mapping[str, str],
) -> tuple[str, ...]:
    raw_inflow = species.get("boundary_inflow", {})
    if not isinstance(raw_inflow, Mapping):
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}].boundary_inflow "
            "must be a table."
        )
    unknown = set(raw_inflow) - _FACE_KEYS
    if unknown:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}].boundary_inflow has "
            "unsupported key(s): "
            + ", ".join(sorted(unknown))
            + "."
        )

    enabled_faces: list[str] = []
    for face in sorted(raw_inflow):
        if raw_inflow[face] != "reservoir":
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].boundary_inflow."
                f'{face} must be "reservoir".'
            )
        if face[0] in periodic_axes:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].boundary_inflow."
                f"{face} cannot inject through a periodic domain face."
            )
        if effective_boundary.get(face, "open") != "open":
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].boundary_inflow."
                f"{face} requires an open particle boundary action."
            )
        enabled_faces.append(face)
    return tuple(enabled_faces)


def _validate_reservoir_injection_common(
    species_table: Mapping[str, Any],
    *,
    index: int,
    use_box: bool,
    batch_duration: float,
) -> None:
    if not use_box:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] reservoir injection "
            "requires [domain]."
        )
    if batch_duration <= 0.0:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] reservoir injection "
            "requires batch_duration > 0."
        )
    source_mode = species_table.get("source_mode", "volume_seed")
    if source_mode != "volume_seed" and species_table.get("boundary_inflow"):
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] cannot combine "
            f'source_mode="{source_mode}" with boundary_inflow.'
        )
    if species_table.get("boundary_inflow"):
        if "source_normal" in species_table:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].source_normal "
                'is only valid for source_mode="plane_source".'
            )
        if "inject_face" in species_table:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].inject_face "
                "is not used by boundary_inflow."
            )


def _validate_reservoir_physics(
    species_table: Mapping[str, Any],
    *,
    index: int,
    source_mode: str,
) -> None:
    has_weight = "w_particle" in species_table
    has_target = "target_macro_particles_per_batch" in species_table
    if has_weight == has_target:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] reservoir injection "
            "requires exactly one of w_particle and target_macro_particles_per_batch."
        )
    if has_weight:
        weight = species_table["w_particle"]
        if (
            not isinstance(weight, (int, float))
            or isinstance(weight, bool)
            or not math.isfinite(float(weight))
            or float(weight) <= 0.0
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].w_particle "
                "must be finite and > 0."
            )
    else:
        target = species_table["target_macro_particles_per_batch"]
        if (
            not isinstance(target, int)
            or isinstance(target, bool)
            or (target < 1 and target != -1)
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}]."
                "target_macro_particles_per_batch must be > 0 or -1."
            )

    if source_mode in _RESERVOIR_SOURCE_MODES:
        if "npcls_per_step" in species_table:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] uses "
                f'source_mode="{source_mode}" and cannot define npcls_per_step.'
            )
        for key in (
            "emit_current_density_a_m2",
            "rays_per_batch",
            "ray_direction",
            "deposit_opposite_charge_on_emit",
        ):
            if key in species_table:
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}] uses "
                    f'source_mode="{source_mode}" and cannot define {key}.'
                )

    velocity_distribution = (
        str(species_table.get("velocity_distribution", "maxwellian"))
        .strip()
        .lower()
    )
    if velocity_distribution == "grid":
        npcls_per_step = species_table.get("npcls_per_step", 0)
        if (
            species_table.get("boundary_inflow")
            and isinstance(npcls_per_step, int)
            and not isinstance(npcls_per_step, bool)
            and npcls_per_step > 0
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] grid "
                "boundary_inflow cannot be combined with positive npcls_per_step."
            )
        if "velocity_grid_path" not in species_table:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] uses "
                'velocity_distribution="grid" and requires velocity_grid_path.'
            )
        _validate_grid_flux_keys(species_table, index=index)
        for key in (
            "number_density_cm3",
            "number_density_m3",
            "temperature_k",
            "temperature_ev",
        ):
            if key in species_table:
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}] uses "
                    f'velocity_distribution="grid" and cannot define {key}.'
                )
        return

    _validate_velocity_grid_forbidden(
        species_table,
        index=index,
        source_mode=source_mode,
    )
    density_keys = [
        key
        for key in ("number_density_cm3", "number_density_m3")
        if key in species_table
    ]
    if len(density_keys) != 1:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] reservoir injection "
            "requires exactly one of number_density_cm3 and number_density_m3."
        )
    density = species_table[density_keys[0]]
    if (
        not isinstance(density, (int, float))
        or isinstance(density, bool)
        or not math.isfinite(float(density))
        or float(density) <= 0.0
    ):
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}].{density_keys[0]} "
            "must be finite and > 0."
        )
    for key in ("temperature_k", "temperature_ev"):
        if key not in species_table:
            continue
        temperature = species_table[key]
        if (
            not isinstance(temperature, (int, float))
            or isinstance(temperature, bool)
            or not math.isfinite(float(temperature))
            or float(temperature) < 0.0
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].{key} "
                "must be finite and >= 0."
            )


def _validate_plane_source_geometry(
    species_table: Mapping[str, Any],
    *,
    index: int,
    use_box: bool,
    box_min: Sequence[float] | None,
    box_max: Sequence[float] | None,
) -> None:
    if not use_box or box_min is None or box_max is None:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] uses "
            'source_mode="plane_source" and requires [domain].'
        )
    pos_low = _maybe_vec3(
        species_table.get("pos_low"),
        name=f"particles.species[{index}].pos_low",
    )
    pos_high = _maybe_vec3(
        species_table.get("pos_high"),
        name=f"particles.species[{index}].pos_high",
    )
    source_normal = _maybe_vec3(
        species_table.get("source_normal"),
        name=f"particles.species[{index}].source_normal",
    )
    if pos_low is None or pos_high is None or source_normal is None:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] plane_source "
            "requires pos_low, pos_high, and source_normal."
        )

    zero_axes = [
        axis
        for axis in range(3)
        if math.isclose(pos_low[axis], pos_high[axis], rel_tol=0.0, abs_tol=1.0e-12)
    ]
    if len(zero_axes) != 1:
        raise ConfigValidationError(
            f"BEACH constraint error: particles.species[{index}] plane_source "
            "pos_low/pos_high must define an axis-aligned zero-thickness rectangle."
        )
    normal_axis = zero_axes[0]
    for axis in range(3):
        low_bound = float(box_min[axis])
        high_bound = float(box_max[axis])
        if axis == normal_axis:
            if not (low_bound < pos_low[axis] < high_bound):
                raise ConfigValidationError(
                    f"BEACH constraint error: particles.species[{index}] plane_source "
                    "must lie strictly inside the box along its normal axis."
                )
        elif not (
            low_bound <= pos_low[axis] < pos_high[axis] <= high_bound
        ):
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}] plane_source "
                "must have positive in-box extent along both tangential axes."
            )

    for axis, component in enumerate(source_normal):
        invalid_component = (
            abs(component) <= 1.0e-12
            if axis == normal_axis
            else abs(component) > 1.0e-12
        )
        if invalid_component:
            raise ConfigValidationError(
                f"BEACH constraint error: particles.species[{index}].source_normal "
                "must be a non-zero axis-aligned vector along the plane normal axis."
            )


def _validate_runtime_mesh(mesh: Mapping[str, Any]) -> None:
    mode = mesh.get("mode", "template")
    if not isinstance(mode, str) or mode not in {"auto", "obj", "template"}:
        raise ConfigValidationError(
            'BEACH constraint error: mesh.mode must be "auto", "obj", or "template".'
        )
    if mode in {"auto", "obj"}:
        _validate_surface_side(
            mesh.get("surface_side"),
            name="mesh.surface_side",
        )
    _validate_surface_model(
        mesh.get("surface_model", "insulator"),
        name="mesh.surface_model",
    )
    _validate_epsilon_r(mesh.get("epsilon_r", 1.0), name="mesh.epsilon_r")
    templates = mesh.get("templates")
    if templates is None:
        return
    if not isinstance(templates, list) or not all(
        isinstance(item, Mapping) for item in templates
    ):
        raise ConfigValidationError(
            "BEACH constraint error: mesh.templates must be an array of tables."
        )
    for index, item in enumerate(templates, start=1):
        _validate_runtime_template(dict(item), index=index)


def _validate_runtime_template(template: Mapping[str, Any], *, index: int) -> None:
    enabled = template.get("enabled", True)
    if not isinstance(enabled, bool):
        raise ConfigValidationError(
            f"BEACH constraint error: mesh.templates[{index}].enabled must be boolean."
        )

    kind_value = template.get("kind", "plane")
    if not isinstance(kind_value, str):
        raise ConfigValidationError(
            f"BEACH constraint error: mesh.templates[{index}].kind must be a string."
        )
    kind = kind_value.strip().lower() or "plane"
    _validate_surface_model(
        template.get("surface_model", "insulator"),
        name=f"mesh.templates[{index}].surface_model",
    )
    _validate_epsilon_r(
        template.get("epsilon_r", 1.0),
        name=f"mesh.templates[{index}].epsilon_r",
    )
    if enabled:
        _validate_surface_side(
            template.get("surface_side"),
            name=f"mesh.templates[{index}].surface_side",
        )

    if "center" in template:
        _maybe_vec3(template.get("center"), name=f"mesh.templates[{index}].center")

    if kind == "plane":
        _positive_template_scalar(template, index=index, key="size_x", default=1.0)
        _positive_template_scalar(template, index=index, key="size_y", default=1.0)
        return

    if kind in {"plate_hole", "plane_hole"}:
        size_x = _positive_template_scalar(
            template, index=index, key="size_x", default=1.0
        )
        size_y = _positive_template_scalar(
            template, index=index, key="size_y", default=1.0
        )
        radius = _positive_template_scalar(
            template, index=index, key="radius", default=0.2
        )
        if radius >= 0.5 * min(size_x, size_y):
            raise ConfigValidationError(
                f"BEACH constraint error: mesh.templates[{index}] radius must be smaller "
                "than half of min(size_x, size_y)."
            )
        return

    if kind == "disk":
        _positive_template_scalar(template, index=index, key="radius", default=0.5)
        return

    if kind == "annulus":
        radius = _positive_template_scalar(
            template, index=index, key="radius", default=0.5
        )
        inner_radius = _nonnegative_template_scalar(
            template,
            index=index,
            key="inner_radius",
            default=0.25,
        )
        if inner_radius >= radius:
            raise ConfigValidationError(
                f"BEACH constraint error: mesh.templates[{index}].inner_radius must be "
                "smaller than radius."
            )
        return

    if kind == "box":
        size = _maybe_vec3(
            template.get("size", [1.0, 1.0, 1.0]),
            name=f"mesh.templates[{index}].size",
        )
        if size is None:
            raise ConfigValidationError(
                f"BEACH constraint error: mesh.templates[{index}].size must be a 3-element array."
            )
        if any(component <= 0.0 for component in size):
            raise ConfigValidationError(
                f"BEACH constraint error: mesh.templates[{index}].size must be positive on all axes."
            )
        return

    if kind == "cylinder":
        _positive_template_scalar(template, index=index, key="radius", default=0.5)
        _positive_template_scalar(template, index=index, key="height", default=1.0)
        return

    if kind == "sphere":
        _positive_template_scalar(template, index=index, key="radius", default=0.5)
        return

    raise ConfigValidationError(
        f"BEACH constraint error: mesh.templates[{index}] has unsupported kind={kind_value!r}."
    )


def _validate_surface_model(value: object, *, name: str) -> None:
    if not isinstance(value, str):
        raise ConfigValidationError(f"BEACH constraint error: {name} must be a string.")
    if value not in {"insulator", "conductor", "dielectric"}:
        raise ConfigValidationError(
            f'BEACH constraint error: {name} must be "insulator", "conductor", or "dielectric".'
        )


def _validate_surface_side(value: object, *, name: str) -> None:
    if not isinstance(value, str):
        raise ConfigValidationError(
            f"BEACH constraint error: {name} must be specified as a string."
        )
    if value not in {"normal_plus", "normal_minus", "outward_closed"}:
        raise ConfigValidationError(
            f'BEACH constraint error: {name} must be "normal_plus", '
            '"normal_minus", or "outward_closed".'
        )


def _mesh_has_surface_model(mesh: Mapping[str, Any], target: str) -> bool:
    mode = str(mesh.get("mode", "template")).strip().lower()
    if mode != "template" and mesh.get("surface_model", "insulator") == target:
        return True
    if mode == "obj":
        return False
    templates = mesh.get("templates")
    if not isinstance(templates, list):
        return False
    for template in templates:
        if not isinstance(template, Mapping):
            continue
        if not bool(template.get("enabled", True)):
            continue
        if template.get("surface_model", "insulator") == target:
            return True
    return False


def _validate_nonnegative_finite_number(value: object, *, name: str) -> None:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ConfigValidationError(f"BEACH constraint error: {name} must be numeric.")
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ConfigValidationError(f"BEACH constraint error: {name} must be finite.")
    if numeric < 0.0:
        raise ConfigValidationError(f"BEACH constraint error: {name} must be >= 0.")


def _validate_epsilon_r(value: object, *, name: str) -> None:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ConfigValidationError(f"BEACH constraint error: {name} must be numeric.")
    if not math.isfinite(float(value)):
        raise ConfigValidationError(f"BEACH constraint error: {name} must be finite.")
    if float(value) < 1.0:
        raise ConfigValidationError(f"BEACH constraint error: {name} must be >= 1.")


def _positive_template_scalar(
    template: Mapping[str, Any],
    *,
    index: int,
    key: str,
    default: float,
) -> float:
    value = _template_scalar(template, index=index, key=key, default=default)
    if value <= 0.0:
        raise ConfigValidationError(
            f"BEACH constraint error: mesh.templates[{index}].{key} must be > 0."
        )
    return value


def _nonnegative_template_scalar(
    template: Mapping[str, Any],
    *,
    index: int,
    key: str,
    default: float,
) -> float:
    value = _template_scalar(template, index=index, key=key, default=default)
    if value < 0.0:
        raise ConfigValidationError(
            f"BEACH constraint error: mesh.templates[{index}].{key} must be >= 0."
        )
    return value


def _template_scalar(
    template: Mapping[str, Any],
    *,
    index: int,
    key: str,
    default: float,
) -> float:
    raw = template.get(key, default)
    if not isinstance(raw, (int, float)) or isinstance(raw, bool):
        raise ConfigValidationError(
            f"BEACH constraint error: mesh.templates[{index}].{key} must be numeric."
        )
    value = float(raw)
    if not math.isfinite(value):
        raise ConfigValidationError(
            f"BEACH constraint error: mesh.templates[{index}].{key} must be finite."
        )
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


def _is_array_of_tables(value: object) -> bool:
    return (
        isinstance(value, list)
        and len(value) > 0
        and all(isinstance(item, Mapping) for item in value)
    )

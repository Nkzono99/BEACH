"""Lower spatial authoring notation to runtime configuration values."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

from ._shared import (
    _FACE_SOURCE_MODES,
    ConfigError,
)


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

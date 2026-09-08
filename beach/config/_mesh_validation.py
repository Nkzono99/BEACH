"""Validate mesh templates, surface materials, and surface sides."""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

from ._shared import (
    ConfigValidationError,
    _maybe_vec3,
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
    if "epsilon_r" in mesh:
        raise ConfigValidationError(
            "BEACH constraint error: mesh.epsilon_r was removed because dielectric polarization is not implemented."
        )
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
    if "epsilon_r" in template:
        raise ConfigValidationError(
            f"BEACH constraint error: mesh.templates[{index}].epsilon_r was removed because "
            "dielectric polarization is not implemented."
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
    if value == "dielectric":
        raise ConfigValidationError(
            f'BEACH constraint error: {name}="dielectric" is not implemented; '
            'use "insulator" for charge accumulation.'
        )
    if value not in {"insulator", "conductor"}:
        raise ConfigValidationError(
            f'BEACH constraint error: {name} must be "insulator" or "conductor".'
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

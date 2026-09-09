"""Validate field boundaries, periodic solvers, and external fields."""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

from ._mesh_validation import (
    _mesh_has_surface_model,
)
from ._shared import (
    ConfigValidationError,
)
from .schema import load_schema


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
    if (has_phi_xy or has_phi_z) and not has_abs:
        raise ConfigValidationError(
            "BEACH constraint error: sim.e0_phi_xy_deg/e0_phi_z_deg require sim.e0_abs."
        )


def _validate_runtime_boundary_tables(
    *,
    domain: Mapping[str, Any] | None,
    field_boundary: Mapping[str, Any] | None,
    particle_boundary: Mapping[str, Any] | None,
    reservoir: Mapping[str, Any] | None,
) -> None:
    """Validate relationships between already schema-checked boundary tables."""

    periodic_axes: set[str] = set()
    if domain is not None:
        if "box_min" not in domain or "box_max" not in domain:
            raise ConfigValidationError(
                "BEACH constraint error: [domain] requires box_min and box_max."
            )
        box_min = domain["box_min"]
        box_max = domain["box_max"]
        if any(box_max[i] <= box_min[i] for i in range(3)):
            raise ConfigValidationError(
                "BEACH constraint error: domain.box_max must be greater than "
                "domain.box_min on every axis."
            )
        periodic_axes = set(domain.get("periodic_axes", []))
    if particle_boundary is not None:
        faces = {"x_low", "x_high", "y_low", "y_high", "z_low", "z_high"}
        if domain is None and faces.intersection(particle_boundary):
            raise ConfigValidationError(
                "BEACH constraint error: [particle_boundary] requires a finite [domain]."
            )
        for face in faces:
            if face in particle_boundary and face[0] in periodic_axes:
                raise ConfigValidationError(
                    f"BEACH constraint error: particle_boundary.{face} cannot "
                    "override a periodic domain face."
                )


def _validate_field_config(
    *, sim: Mapping[str, Any], domain: Mapping[str, Any] | None,
    field_boundary: Mapping[str, Any] | None, reservoir: Mapping[str, Any] | None,
    periodic2_config: object, mesh: Mapping[str, Any], resolved_batch_duration: float,
    use_box: bool,
) -> tuple[float, str]:
    """Validate field and periodic solver choices; return the adaptive limit and boundary mode."""
    if isinstance(periodic2_config, Mapping):
        schema, _ = load_schema()
        periodic_defaults = {
            key: rule["default"]
            for key, rule in schema["properties"]["periodic2"]["properties"].items()
            if "default" in rule
        }
        periodic2_config = {**periodic_defaults, **periodic2_config}
    adaptive_nonzero_mode_limit = float(periodic2_config.get(
        "max_nonzero_mode_potential_step", 0.0
    )) if isinstance(periodic2_config, Mapping) else 0.0

    field_bc_mode = (
        field_boundary.get("mode", "free") if field_boundary is not None else "free"
    )
    field_solver = sim.get("field_solver", "auto")
    if isinstance(periodic2_config, Mapping):
        if field_bc_mode != "periodic2" or not use_box:
            raise ConfigValidationError(
                'BEACH constraint error: [periodic2] requires field_boundary.mode="periodic2" and [domain].'
            )
        nonzero_backend = periodic2_config["nonzero_mode_backend"]
        if nonzero_backend == "panel_spectral_reference" and field_solver != "direct":
            raise ConfigValidationError(
                'BEACH constraint error: panel_spectral_reference requires sim.field_solver="direct".'
            )
        if nonzero_backend == "cached_kneq0" and (
            field_solver != "fmm"
            or sim.get("field_periodic_far_correction", "none") != "cached_kneq0"
        ):
            raise ConfigValidationError(
                'BEACH constraint error: periodic2 cached_kneq0 requires sim.field_solver="fmm" '
                'and sim.field_periodic_far_correction="cached_kneq0".'
            )
    if field_solver != "direct":
        tree_theta = sim.get("tree_theta", 0.5)
        if not 0.0 < tree_theta <= 1.0:
            raise ConfigValidationError(
                "BEACH constraint error: sim.tree_theta must satisfy 0 < theta <= 1 "
                "for tree-capable solvers."
            )
        for key, default in (("tree_leaf_max", 16), ("tree_min_nelem", 256)):
            value = sim.get(key, default)
            if value < 1:
                raise ConfigValidationError(
                    f"BEACH constraint error: sim.{key} must be an integer >= 1 "
                    "for tree-capable solvers."
                )
    if (
        field_bc_mode == "periodic2"
        and sim.get("field_periodic_far_correction", "none") == "cached_kneq0"
    ):
        image_layers = sim.get("field_periodic_image_layers", 1)
        ewald_layers = sim.get("field_periodic_ewald_layers", 4)
        ewald_alpha = sim.get("field_periodic_ewald_alpha", 0.0)
        if image_layers < 1:
            raise ConfigValidationError(
                "BEACH constraint error: sim.field_periodic_image_layers must be an integer >= 1 "
                "for cached_kneq0."
            )
        if ewald_layers < 1:
            raise ConfigValidationError(
                "BEACH constraint error: sim.field_periodic_ewald_layers must be an integer >= 1 "
                "for cached_kneq0."
            )
        if ewald_alpha < 0.0:
            raise ConfigValidationError(
                "BEACH constraint error: sim.field_periodic_ewald_alpha must be finite and >= 0 "
                "for cached_kneq0."
            )
        tolerance = sim.get("field_periodic_generation_tolerance", 1.0e-8)
        if tolerance <= 0.0:
            raise ConfigValidationError(
                "BEACH constraint error: sim.field_periodic_generation_tolerance must be "
                "finite and > 0 for cached_kneq0."
            )
        cache_dir = sim.get("field_periodic_cache_dir", ".beach_cache/periodic2")
        if not cache_dir:
            raise ConfigValidationError(
                "BEACH constraint error: sim.field_periodic_cache_dir must be a non-empty "
                "string for cached_kneq0."
            )
    if field_bc_mode == "periodic2":
        if sim.get("field_periodic_image_layers", 1) < 0:
            raise ConfigValidationError(
                "BEACH constraint error: sim.field_periodic_image_layers must be >= 0 for periodic2."
            )
        supported_lower_boundaries = {"e_bottom_zero", "symmetric_vacuum"}
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
    return adaptive_nonzero_mode_limit, field_bc_mode

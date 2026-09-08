"""Validate field boundaries, periodic solvers, and external fields."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

from ._mesh_validation import (
    _mesh_has_surface_model,
)
from ._shared import (
    ConfigValidationError,
    _maybe_vec3,
)


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


def _validate_field_config(
    *, sim: Mapping[str, Any], domain: Mapping[str, Any] | None,
    field_boundary: Mapping[str, Any] | None, reservoir: Mapping[str, Any] | None,
    periodic2_config: object, mesh: Mapping[str, Any], resolved_batch_duration: float,
    use_box: bool,
) -> tuple[float, str]:
    """Validate field and periodic solver choices; return the adaptive limit and boundary mode."""
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
        field_boundary.get("mode", "free") if field_boundary is not None else "free"
    )
    field_solver = sim.get("field_solver", "auto")
    if field_solver not in {"direct", "treecode", "fmm", "auto"}:
        raise ConfigValidationError(
            'BEACH constraint error: sim.field_solver must be "direct", "treecode", '
            '"fmm", or "auto".'
        )
    if field_solver != "direct":
        tree_theta = sim.get("tree_theta", 0.5)
        if (
            not isinstance(tree_theta, (int, float))
            or isinstance(tree_theta, bool)
            or not math.isfinite(float(tree_theta))
            or not 0.0 < float(tree_theta) <= 1.0
        ):
            raise ConfigValidationError(
                "BEACH constraint error: sim.tree_theta must satisfy 0 < theta <= 1 "
                "for tree-capable solvers."
            )
        for key, default in (("tree_leaf_max", 16), ("tree_min_nelem", 256)):
            value = sim.get(key, default)
            if not isinstance(value, int) or isinstance(value, bool) or value < 1:
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
        if (
            not isinstance(image_layers, int)
            or isinstance(image_layers, bool)
            or image_layers < 1
        ):
            raise ConfigValidationError(
                "BEACH constraint error: sim.field_periodic_image_layers must be an integer >= 1 "
                "for cached_kneq0."
            )
        if (
            not isinstance(ewald_layers, int)
            or isinstance(ewald_layers, bool)
            or ewald_layers < 1
        ):
            raise ConfigValidationError(
                "BEACH constraint error: sim.field_periodic_ewald_layers must be an integer >= 1 "
                "for cached_kneq0."
            )
        if (
            not isinstance(ewald_alpha, (int, float))
            or isinstance(ewald_alpha, bool)
            or not math.isfinite(float(ewald_alpha))
            or float(ewald_alpha) < 0.0
        ):
            raise ConfigValidationError(
                "BEACH constraint error: sim.field_periodic_ewald_alpha must be finite and >= 0 "
                "for cached_kneq0."
            )
        tolerance = sim.get("field_periodic_generation_tolerance", 1.0e-8)
        if (
            not isinstance(tolerance, (int, float))
            or isinstance(tolerance, bool)
            or not math.isfinite(float(tolerance))
            or float(tolerance) <= 0.0
        ):
            raise ConfigValidationError(
                "BEACH constraint error: sim.field_periodic_generation_tolerance must be "
                "finite and > 0 for cached_kneq0."
            )
        cache_dir = sim.get("field_periodic_cache_dir", ".beach_cache/periodic2")
        if not isinstance(cache_dir, str) or not cache_dir:
            raise ConfigValidationError(
                "BEACH constraint error: sim.field_periodic_cache_dir must be a non-empty "
                "string for cached_kneq0."
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
    return adaptive_nonzero_mode_limit, field_bc_mode

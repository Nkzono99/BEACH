"""Read and validate the field-reconstruction receipt from summary metadata."""

from __future__ import annotations

import numpy as np

from ._summary_values import _parse_nonnegative_finite_float, _parse_nonnegative_int
from .types import FieldReconstructionReceipt


def _parse_field_reconstruction_receipt(
    data: dict[str, str],
) -> FieldReconstructionReceipt | None:
    schema_key = "field_reconstruction_schema_version"
    if schema_key not in data:
        return None
    schema_version = _parse_nonnegative_int(data[schema_key], key=schema_key)
    if schema_version != 2:
        raise ValueError(
            "summary.txt field_reconstruction_schema_version is not supported."
        )

    required = (
        "field_reconstruction_resolved_field_solver",
        "field_reconstruction_fmm_expansion_order",
        "field_reconstruction_field_bc_mode",
        "field_reconstruction_tree_theta",
        "field_reconstruction_tree_leaf_max",
        "field_reconstruction_e0_V_m",
        "field_reconstruction_use_box",
        "field_reconstruction_box_min_m",
        "field_reconstruction_box_max_m",
        "field_reconstruction_boundary_low",
        "field_reconstruction_boundary_high",
        "field_reconstruction_periodic_image_layers",
        "field_reconstruction_periodic_far_correction",
        "field_reconstruction_periodic_nonzero_mode_backend",
        "field_reconstruction_periodic_zero_mode_policy",
        "field_reconstruction_periodic_lower_boundary_model",
        "field_reconstruction_periodic_reference_mode_layers",
        "field_reconstruction_periodic_panel_quadrature_order",
        "field_reconstruction_periodic_ewald_alpha",
        "field_reconstruction_periodic_ewald_layers",
        "field_reconstruction_periodic_cache_dir",
        "field_reconstruction_periodic_generation_tolerance",
    )
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(
            "summary.txt field reconstruction receipt requires "
            + ", ".join(missing)
            + "."
        )

    resolved_field_solver = data[
        "field_reconstruction_resolved_field_solver"
    ].strip().lower()
    if resolved_field_solver not in {"direct", "treecode", "fmm"}:
        raise ValueError(
            "summary.txt field_reconstruction_resolved_field_solver must be "
            "direct, treecode, or fmm."
        )
    fmm_expansion_order = _parse_positive_int(
        data["field_reconstruction_fmm_expansion_order"],
        key="field_reconstruction_fmm_expansion_order",
    )
    field_bc_mode = data["field_reconstruction_field_bc_mode"].strip().lower()
    if field_bc_mode not in {"free", "periodic2"}:
        raise ValueError(
            "summary.txt field_reconstruction_field_bc_mode must be free or periodic2."
        )
    tree_theta = _parse_positive_finite_float(
        data["field_reconstruction_tree_theta"],
        key="field_reconstruction_tree_theta",
    )
    if tree_theta > 1.0:
        raise ValueError(
            "summary.txt field_reconstruction_tree_theta must be <= 1."
        )
    tree_leaf_max = _parse_positive_int(
        data["field_reconstruction_tree_leaf_max"],
        key="field_reconstruction_tree_leaf_max",
    )
    use_box = _parse_required_bool(
        data["field_reconstruction_use_box"],
        key="field_reconstruction_use_box",
    )
    boundary_low = _parse_int3(
        data["field_reconstruction_boundary_low"],
        key="field_reconstruction_boundary_low",
    )
    boundary_high = _parse_int3(
        data["field_reconstruction_boundary_high"],
        key="field_reconstruction_boundary_high",
    )
    if any(value not in {0, 1, 2, 3} for value in (*boundary_low, *boundary_high)):
        raise ValueError(
            "summary.txt field reconstruction boundaries contain an unknown code."
        )
    box_min_m = _parse_float3(
        data["field_reconstruction_box_min_m"],
        key="field_reconstruction_box_min_m",
    )
    box_max_m = _parse_float3(
        data["field_reconstruction_box_max_m"],
        key="field_reconstruction_box_max_m",
    )
    if use_box and any(high <= low for low, high in zip(box_min_m, box_max_m)):
        raise ValueError(
            "summary.txt field reconstruction box_max must exceed box_min."
        )
    if field_bc_mode == "periodic2":
        periodic_axes = [
            axis
            for axis, (low, high) in enumerate(zip(boundary_low, boundary_high))
            if low == 2 and high == 2
        ]
        has_one_sided_periodic = any(
            (low == 2) != (high == 2)
            for low, high in zip(boundary_low, boundary_high)
        )
        if not use_box or has_one_sided_periodic or len(periodic_axes) != 2:
            raise ValueError(
                "summary.txt periodic2 field reconstruction requires a box and "
                "exactly two paired periodic axes."
            )
    periodic_far_correction = data[
        "field_reconstruction_periodic_far_correction"
    ].strip().lower()
    if periodic_far_correction not in {"auto", "none", "cached_kneq0"}:
        raise ValueError(
            "summary.txt field reconstruction far correction is not supported."
        )
    periodic_nonzero_mode_backend = data[
        "field_reconstruction_periodic_nonzero_mode_backend"
    ].strip().lower()
    periodic_zero_mode_policy = data[
        "field_reconstruction_periodic_zero_mode_policy"
    ].strip().lower()
    periodic_lower_boundary_model = data[
        "field_reconstruction_periodic_lower_boundary_model"
    ].strip().lower()
    periodic_reference_mode_layers = _parse_positive_int(
        data["field_reconstruction_periodic_reference_mode_layers"],
        key="field_reconstruction_periodic_reference_mode_layers",
    )
    periodic_panel_quadrature_order = _parse_positive_int(
        data["field_reconstruction_periodic_panel_quadrature_order"],
        key="field_reconstruction_periodic_panel_quadrature_order",
    )
    if periodic_panel_quadrature_order < 2:
        raise ValueError(
            "summary.txt periodic panel quadrature order must be >= 2."
        )
    _validate_periodic_reconstruction_model(
        resolved_field_solver=resolved_field_solver,
        field_bc_mode=field_bc_mode,
        far_correction=periodic_far_correction,
        nonzero_mode_backend=periodic_nonzero_mode_backend,
        zero_mode_policy=periodic_zero_mode_policy,
        lower_boundary_model=periodic_lower_boundary_model,
    )
    periodic_image_layers = _parse_nonnegative_int(
        data["field_reconstruction_periodic_image_layers"],
        key="field_reconstruction_periodic_image_layers",
    )
    periodic_ewald_layers = _parse_nonnegative_int(
        data["field_reconstruction_periodic_ewald_layers"],
        key="field_reconstruction_periodic_ewald_layers",
    )
    if periodic_far_correction == "cached_kneq0" and (
        periodic_image_layers < 1 or periodic_ewald_layers < 1
    ):
        raise ValueError(
            "summary.txt cached_kneq0 field reconstruction requires positive "
            "image and Ewald layer counts."
        )
    periodic_cache_dir = data[
        "field_reconstruction_periodic_cache_dir"
    ].strip()
    if not periodic_cache_dir:
        raise ValueError(
            "summary.txt field reconstruction cache directory must not be empty."
        )
    periodic_ewald_alpha = _parse_nonnegative_finite_float(
        data["field_reconstruction_periodic_ewald_alpha"],
        key="field_reconstruction_periodic_ewald_alpha",
    )
    periodic_generation_tolerance = _parse_positive_finite_float(
        data["field_reconstruction_periodic_generation_tolerance"],
        key="field_reconstruction_periodic_generation_tolerance",
    )
    return FieldReconstructionReceipt(
        schema_version=schema_version,
        resolved_field_solver=resolved_field_solver,
        fmm_expansion_order=fmm_expansion_order,
        field_bc_mode=field_bc_mode,
        tree_theta=tree_theta,
        tree_leaf_max=tree_leaf_max,
        external_e0_v_m=_parse_float3(
            data["field_reconstruction_e0_V_m"],
            key="field_reconstruction_e0_V_m",
        ),
        use_box=use_box,
        box_min_m=box_min_m,
        box_max_m=box_max_m,
        boundary_low=boundary_low,
        boundary_high=boundary_high,
        periodic_image_layers=periodic_image_layers,
        periodic_far_correction=periodic_far_correction,
        periodic_nonzero_mode_backend=periodic_nonzero_mode_backend,
        periodic_zero_mode_policy=periodic_zero_mode_policy,
        periodic_lower_boundary_model=periodic_lower_boundary_model,
        periodic_reference_mode_layers=periodic_reference_mode_layers,
        periodic_panel_quadrature_order=periodic_panel_quadrature_order,
        periodic_ewald_alpha=periodic_ewald_alpha,
        periodic_ewald_layers=periodic_ewald_layers,
        periodic_cache_dir=periodic_cache_dir,
        periodic_generation_tolerance=periodic_generation_tolerance,
    )


def _validate_periodic_reconstruction_model(
    *,
    resolved_field_solver: str,
    field_bc_mode: str,
    far_correction: str,
    nonzero_mode_backend: str,
    zero_mode_policy: str,
    lower_boundary_model: str,
) -> None:
    if field_bc_mode == "free":
        if far_correction not in {"auto", "none"}:
            raise ValueError(
                "summary.txt free-space reconstruction has an incompatible "
                "periodic far correction."
            )
        expected = ("not_applicable", "not_applicable", "not_applicable")
    elif nonzero_mode_backend == "legacy_finite_images":
        if resolved_field_solver != "fmm":
            raise ValueError(
                "summary.txt legacy finite-image reconstruction requires the fmm solver."
            )
        expected = (
            "legacy_finite_images",
            "legacy_not_decomposed",
            "legacy_implicit",
        )
        if far_correction not in {"auto", "none"}:
            raise ValueError(
                "summary.txt legacy finite-image reconstruction has an "
                "incompatible far correction."
            )
    elif nonzero_mode_backend in {"cached_kneq0", "panel_spectral_reference"}:
        if zero_mode_policy != "exclude_k0" or lower_boundary_model not in {
            "e_bottom_zero",
            "symmetric_vacuum",
        }:
            raise ValueError(
                "summary.txt split-periodic reconstruction has an incompatible "
                "zero-mode contract."
            )
        if nonzero_mode_backend == "cached_kneq0" and far_correction != "cached_kneq0":
            raise ValueError(
                "summary.txt cached periodic backend and far correction disagree."
            )
        if nonzero_mode_backend == "cached_kneq0" and resolved_field_solver != "fmm":
            raise ValueError(
                "summary.txt cached periodic reconstruction requires the fmm solver."
            )
        if nonzero_mode_backend == "panel_spectral_reference" and far_correction not in {
            "auto",
            "none",
        }:
            raise ValueError(
                "summary.txt panel spectral reference has an incompatible far correction."
            )
        if (
            nonzero_mode_backend == "panel_spectral_reference"
            and resolved_field_solver != "direct"
        ):
            raise ValueError(
                "summary.txt panel spectral reference reconstruction requires the direct solver."
            )
        return
    else:
        raise ValueError(
            "summary.txt periodic nonzero-mode reconstruction backend is not supported."
        )

    actual = (
        nonzero_mode_backend,
        zero_mode_policy,
        lower_boundary_model,
    )
    if actual != expected:
        raise ValueError(
            "summary.txt field reconstruction periodic model is inconsistent."
        )


def _parse_positive_int(value: str, *, key: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise ValueError(f"summary.txt {key} must be > 0.")
    return parsed


def _parse_positive_finite_float(value: str, *, key: str) -> float:
    parsed = float(value)
    if not np.isfinite(parsed) or parsed <= 0.0:
        raise ValueError(f"summary.txt {key} must be finite and > 0.")
    return parsed


def _parse_required_bool(value: str, *, key: str) -> bool:
    parsed = value.strip().lower()
    if parsed in {"t", "true"}:
        return True
    if parsed in {"f", "false"}:
        return False
    raise ValueError(f"summary.txt {key} must be true or false.")


def _parse_float3(value: str, *, key: str) -> tuple[float, float, float]:
    fields = value.replace(",", " ").split()
    if len(fields) != 3:
        raise ValueError(f"summary.txt {key} must contain exactly three values.")
    parsed = tuple(float(field) for field in fields)
    if not all(np.isfinite(component) for component in parsed):
        raise ValueError(f"summary.txt {key} values must be finite.")
    return parsed


def _parse_int3(value: str, *, key: str) -> tuple[int, int, int]:
    fields = value.replace(",", " ").split()
    if len(fields) != 3:
        raise ValueError(f"summary.txt {key} must contain exactly three values.")
    parsed = tuple(int(field) for field in fields)
    return parsed

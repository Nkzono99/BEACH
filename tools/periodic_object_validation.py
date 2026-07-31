#!/usr/bin/env python3
"""Stage and verify archived/new periodic-object validation runs."""

from __future__ import annotations

import argparse
from collections import Counter
import copy
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Iterable, Mapping, Sequence

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10
    import tomli as tomllib

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from beach import (  # noqa: E402
    AdhesionProfile,
    Beach,
    FieldKernel,
    FortranRunResult,
    FortranChargeHistory,
    ObjectInteractionSnapshot,
    finite_shell_convergence,
    field_kernel_build_info,
)
from beach.fortran_results.constants import K_COULOMB  # noqa: E402


JOB_TEMPLATE = REPO_ROOT / "examples/job_scripts/periodic_object_validation_sysa.sh"
ANALYSIS_JOB_TEMPLATE = (
    REPO_ROOT / "examples/job_scripts/periodic_object_analysis_sysa.sh"
)
DEFAULT_PARTITION = "gr20001a"
DEFAULT_EWALD_LAYERS = 4
DEFAULT_GENERATION_TOLERANCE = 1.0e-8
TARGET_GEOMETRY_RADIUS_RELATIVE_TOLERANCE = 5.0e-3
TARGET_GEOMETRY_REPRESENTATION = "periodic2_mesh_connected"
PHYSICAL_OBJECT_COMPONENTS = (
    "other_objects_all_images",
    "target_periodic_images",
    "external_uniform",
    "total_external",
)
ADDITIVE_OBJECT_COMPONENTS = PHYSICAL_OBJECT_COMPONENTS[:-1]
CACHED_NUMERICAL_COMPONENTS = (
    "periodic_kneq0",
    "physical_k0",
    "primary_free_subtraction",
    "cached_kneq0_trace_correction",
)
ORACLE_UNIFORM_RELATIVE_TOLERANCE = 1.2e-1
ORACLE_COSINE_FINE_RELATIVE_TOLERANCE = 8.0e-2
ORACLE_COSINE_SAMPLE_ABS_Z_M = (0.25, 0.50)
ORACLE_COSINE_EXPECTED_DECAY_RATIO = math.exp(-math.pi * 0.25)
ORACLE_COSINE_DECAY_RATIO_RELATIVE_TOLERANCE = 1.8e-1
ORACLE_QUADRATURE_RELATIVE_TOLERANCE = 1.0e-2
ORACLE_EFFECTIVE_FIELD_CONTRACT = {
    "requested_periodic_model": "infinite_physical",
    "configured_far_correction": "none",
    "effective_far_correction": "cached_kneq0",
    "nonzero_mode_backend": "cached_kneq0",
    "zero_mode_policy": "exclude_k0",
    "lower_boundary_model": "e_bottom_zero",
}
ORACLE_CACHE_EVALUATION_LABELS = ("uniform_4", "cosine_4", "cosine_8")
ORACLE_CACHE_EVALUATION_GROUPS = {
    "4": ("uniform_4", "cosine_4"),
    "8": ("cosine_8",),
}
BUILD_INFO_KEYS = (
    "build_info_schema_version",
    "build_version",
    "build_version_mode",
    "build_source_commit",
    "build_id",
)
SYS_MODULE_PATTERN = re.compile(r"(?<![A-Za-z0-9_])Sys(?:CL|A|B|C|G)/[^\s()]+")
INTEL_MODULE_PATTERN = re.compile(r"(?<![A-Za-z0-9_])intel(?:mpi)?/[^\s()]+")
SAFE_STAGE_PATH_PATTERN = re.compile(r"[A-Za-z0-9_./:+-]+\Z")
PRODUCTION_FIELD_EXECUTION_CONTRACT = {
    "field_backend": "fmm",
    "field_normalization": "si",
    "field_source_model": "triangle_p0",
    "field_kernel_id": "triangle_p0_exact_p2m_near",
}
EXPECTED_RESOURCES = {
    "partition": DEFAULT_PARTITION,
    "mpi_processes": 6,
    "openmp_threads": 112,
    "cores_per_process": 112,
}
RECEIPT_POLICY = {
    "schema_version": 1,
    "write_mode": "create_once_compare",
    "hash_algorithm": "sha256",
    "output_inventory": "all_regular_files",
    "strict_requires_existing": True,
}
STRICT_CASES = (
    "cache_prime",
    "smoke_finite_configured",
    "smoke_infinite_physical",
    "full_finite_configured_140000",
    "full_finite_configured_280000",
    "full_infinite_physical_140000",
    "full_infinite_physical_280000",
)
EXPECTED_SUBMITTED_JOBS = {
    "smoke": ("beach-periodic-smoke", "p=6:t=112:c=112"),
    "finite_140000": ("beach-finite-140k", "p=6:t=112:c=112"),
    "finite_280000": ("beach-finite-280k", "p=6:t=112:c=112"),
    "infinite_140000": ("beach-infinite-140k", "p=6:t=112:c=112"),
    "infinite_280000": ("beach-infinite-280k", "p=6:t=112:c=112"),
    "analysis": ("beach-periodic-analysis", "p=1:t=28:c=28"),
}
CASE_PRODUCER_JOB_ROLES = {
    "cache_prime": "smoke",
    "smoke_finite_configured": "smoke",
    "smoke_infinite_physical": "smoke",
    "full_finite_configured_140000": "finite_140000",
    "full_finite_configured_280000": "finite_280000",
    "full_infinite_physical_140000": "infinite_140000",
    "full_infinite_physical_280000": "infinite_280000",
}
PRODUCER_ROLE_CASES = {
    role: tuple(
        case_name
        for case_name, producer_role in CASE_PRODUCER_JOB_ROLES.items()
        if producer_role == role
    )
    for role in EXPECTED_SUBMITTED_JOBS
    if role != "analysis"
}
EXPECTED_CASE_GRAPH: dict[str, dict[str, Any]] = {
    "cache_prime": {
        "role": "cache_prime",
        "periodic_model": "infinite_physical",
        "batch_count": 1,
        "history_stride": 1,
        "resume": False,
        "restart_case": None,
        "cache_expectation": "miss",
        "config_relative": "input/cache_prime/infinite_physical.toml",
        "output_relative": "run/cache_prime/infinite_physical",
    },
    "smoke_finite_configured": {
        "role": "smoke",
        "periodic_model": "finite_configured",
        "batch_count": 100,
        "history_stride": 10,
        "resume": False,
        "restart_case": None,
        "cache_expectation": "none",
        "config_relative": "input/smoke/finite_configured.toml",
        "output_relative": "run/smoke/finite_configured",
    },
    "smoke_infinite_physical": {
        "role": "smoke",
        "periodic_model": "infinite_physical",
        "batch_count": 100,
        "history_stride": 10,
        "resume": False,
        "restart_case": None,
        "cache_expectation": "hit",
        "config_relative": "input/smoke/infinite_physical.toml",
        "output_relative": "run/smoke/infinite_physical",
    },
    "full_finite_configured": {
        "role": "canonical",
        "periodic_model": "finite_configured",
        "batch_count": 280000,
        "history_stride": 1000,
        "resume": False,
        "restart_case": None,
        "cache_expectation": "none",
        "config_relative": "input/full/finite_configured.toml",
        "output_relative": "run/full/finite_configured/canonical",
    },
    "full_infinite_physical": {
        "role": "canonical",
        "periodic_model": "infinite_physical",
        "batch_count": 280000,
        "history_stride": 1000,
        "resume": False,
        "restart_case": None,
        "cache_expectation": "hit",
        "config_relative": "input/full/infinite_physical.toml",
        "output_relative": "run/full/infinite_physical/canonical",
    },
    "full_finite_configured_140000": {
        "role": "segment",
        "periodic_model": "finite_configured",
        "batch_count": 140000,
        "history_stride": 1000,
        "resume": False,
        "restart_case": None,
        "cache_expectation": "none",
        "config_relative": "input/full/segments/finite_configured_140000.toml",
        "output_relative": "run/full/finite_configured/140000",
    },
    "full_finite_configured_280000": {
        "role": "segment",
        "periodic_model": "finite_configured",
        "batch_count": 280000,
        "history_stride": 1000,
        "resume": True,
        "restart_case": "full_finite_configured_140000",
        "cache_expectation": "none",
        "config_relative": "input/full/segments/finite_configured_280000.toml",
        "output_relative": "run/full/finite_configured/280000",
    },
    "full_infinite_physical_140000": {
        "role": "segment",
        "periodic_model": "infinite_physical",
        "batch_count": 140000,
        "history_stride": 1000,
        "resume": False,
        "restart_case": None,
        "cache_expectation": "hit",
        "config_relative": "input/full/segments/infinite_physical_140000.toml",
        "output_relative": "run/full/infinite_physical/140000",
    },
    "full_infinite_physical_280000": {
        "role": "segment",
        "periodic_model": "infinite_physical",
        "batch_count": 280000,
        "history_stride": 1000,
        "resume": True,
        "restart_case": "full_infinite_physical_140000",
        "cache_expectation": "hit",
        "config_relative": "input/full/segments/infinite_physical_280000.toml",
        "output_relative": "run/full/infinite_physical/280000",
    },
}
ARCHIVE_REQUIRED_ANALYSIS_OUTPUTS = (
    "work/latest/summary.txt",
    "work/latest/charges.csv",
    "work/latest/mesh_triangles.csv",
    "work/latest/mesh_sources.csv",
    "work/latest/mesh_potential.csv",
    "work/latest/charge_history.csv",
)
PRODUCTION_MECHANICS_INPUTS = (
    "input/release_kernel_base.toml",
    "analysis/local_release/release_model_summary.json",
)
LEGACY_ESTIMATOR_INPUTS = (
    "analysis/local_release/moving_sphere_model_summary.json",
    "analysis/local_release/moving_sphere_force_curves.csv",
    "analysis/local_release/moving_sphere_release_summary.csv",
    "analysis/local_release/force_timeseries.csv",
)
LEGACY_NATIVE_KEYS = ((149001, 7), (180001, 6), (279001, 6), (279001, 7))
LEGACY_CURVE_FIELDS = (
    "target_mesh_id",
    "batch",
    "time_s",
    "target_charge_multiplier",
    "displacement_z_m",
    "f_coulomb_z_n",
    "f_resist_n",
    "net_force_z_n",
)
LEGACY_RELEASE_FIELDS = (
    "target_mesh_id",
    "batch",
    "time_s",
    "processed_particles",
    "rel_change",
    "target_charge_multiplier",
    "radius_m",
    "z_max_m",
    "z_samples",
    "q_target_base_c",
    "q_target_effective_c",
    "q_source_c",
    "f_coulomb_initial_z_n",
    "f_coulomb_final_z_n",
    "f_coulomb_max_z_n",
    "f_coulomb_min_z_n",
    "f_adh_n",
    "f_gravity_n",
    "f_resist_n",
    "w_coulomb_signed_j",
    "w_resist_j",
    "w_net_signed_j",
    "w_net_positive_part_j",
    "v_release_signed_net_m_per_s",
    "v_release_positive_part_m_per_s",
    "crossed_force_threshold",
    "first_crossing_z_m",
)
LEGACY_TIMESERIES_FIELDS = (
    "top_mesh_id",
    "radius_m",
    "local_charge_multiplier",
    "batch",
    "time_s",
    "processed_particles",
    "rel_change",
    "q_base_signed_c",
    "q_base_abs_c",
    "q_effective_abs_c",
    "f_direct_object_z_n",
    "f_local_pair_proxy_n",
    "f_adh_n",
    "f_gravity_n",
    "f_resist_n",
    "net_direct_z_n",
    "net_local_pair_n",
    "w_elec_proxy_j",
    "w_barrier_j",
    "w_release_proxy_j",
    "v_release_proxy_m_per_s",
)
ALLOWED_DIFFERENCES = [
    "output.dir",
    "output.resume",
    "output.restart_from",
    "output.history_stride",
    "sim.batch_count",
    "sim.field_periodic_far_correction",
    "sim.field_periodic_image_layers(default=1 explicit normalization)",
    "sim.field_periodic_ewald_layers(default=4 normalization)",
    "sim.field_periodic_cache_dir",
    "sim.field_periodic_generation_tolerance",
]
CSV_SCHEMAS: dict[str, tuple[str, ...]] = {
    "run_summary.csv": (
        "case",
        "status",
        "source_version",
        "batches",
        "mesh_nelem",
        "mesh_count",
        "processed_particles",
        "absorbed",
        "escaped",
        "escaped_boundary",
        "survived_max_step",
        "last_rel_change",
        "total_charge_C",
        "model_fingerprint",
        "mesh_fingerprint",
        "species_fingerprint",
        "message",
    ),
    "charge_history_pair.csv": (
        "case",
        "batch",
        "processed_particles",
        "rel_change",
        "total_charge_C",
        "source_segment",
    ),
    "particle_ledger_pair.csv": (
        "case",
        "status",
        "species_idx",
        "batch",
        "injected_count",
        "absorbed_count",
        "escaped_count",
        "residual_C",
    ),
    "mesh_potential_pair.csv": ("case", "elem_idx", "potential_V", "status"),
    "snapshot_manifest.csv": (
        "snapshot_id",
        "case",
        "sample_kind",
        "resolved_batch",
        "run_final_batch",
        "is_latest_history",
        "mesh_fingerprint",
        "field_source_model",
        "field_kernel_id",
        "scope",
        "mesh_id",
        "processed_particles",
        "rel_change",
        "charge_vector_sha256",
        "source_file",
        "source_file_sha256",
        "charge_sum_C",
        "element_charge_l1_C",
        "element_charge_l2_C",
        "element_charge_linf_C",
        "status",
        "message",
    ),
    "object_wrench.csv": (
        "case",
        "periodic_model",
        "effective_far_correction",
        "zero_mode_policy",
        "lower_boundary_model",
        "periodic_cache_hit",
        "periodic_cache_build_count",
        "periodic_cache_fingerprint",
        "periodic_cache_path",
        "periodic_cache_path_sha256",
        "periodic_cache_prime_receipt_sha256",
        "step_selector",
        "resolved_batch",
        "mesh_id",
        "model_radius_m",
        "radius_source",
        "geometry_radius_m",
        "radius_relative_mismatch",
        "radius_relative_tolerance",
        "target_geometry_representation",
        "component",
        "component_kind",
        "force_x_N",
        "force_y_N",
        "force_z_N",
        "torque_x_Nm",
        "torque_y_Nm",
        "torque_z_Nm",
        "torque_origin_policy",
        "torque_origin_x_m",
        "torque_origin_y_m",
        "torque_origin_z_m",
        "potential_energy_J",
        "total_charge_C",
        "status",
        "message",
    ),
    "object_path_curves.csv": (
        "case",
        "periodic_model",
        "effective_far_correction",
        "zero_mode_policy",
        "lower_boundary_model",
        "step_selector",
        "resolved_batch",
        "mesh_id",
        "model_radius_m",
        "radius_source",
        "geometry_radius_m",
        "radius_relative_mismatch",
        "radius_relative_tolerance",
        "target_geometry_representation",
        "component",
        "displacement_m",
        "force_x_N",
        "force_y_N",
        "force_z_N",
        "torque_x_Nm",
        "torque_y_Nm",
        "torque_z_Nm",
        "torque_origin_policy",
        "torque_origin_x_m",
        "torque_origin_y_m",
        "torque_origin_z_m",
        "electrostatic_work_J",
        "potential_difference_work_J",
        "gravity_work_J",
        "adhesion_work_J",
        "available_energy_J",
        "electrostatic_only_speed_m_s",
        "gravity_corrected_speed_m_s",
        "speed_m_s",
        "eta_translation",
        "status",
    ),
    "object_path_summary.csv": (
        "case",
        "periodic_model",
        "effective_far_correction",
        "zero_mode_policy",
        "lower_boundary_model",
        "step_selector",
        "resolved_batch",
        "mesh_id",
        "radius_m",
        "radius_source",
        "geometry_radius_m",
        "radius_relative_mismatch",
        "radius_relative_tolerance",
        "target_geometry_representation",
        "torque_origin_policy",
        "initial_torque_origin_x_m",
        "initial_torque_origin_y_m",
        "initial_torque_origin_z_m",
        "mass_kg",
        "gravity_m_s2",
        "adhesion_model",
        "adhesion_force_N",
        "adhesion_work_J",
        "adhesion_equivalent_range_m",
        "eta_translation",
        "eta_translation_sensitivity",
        "energy_tolerance_J",
        "path_start_m",
        "path_end_m",
        "path_status",
        "potential_work_available",
        "numerically_qualified",
        "endpoint_work_J",
        "max_force_z_N",
        "work_relative_mismatch",
        "path_relative_tolerance",
        "path_force_absolute_tolerance_N",
        "path_work_absolute_tolerance_J",
        "path_max_refinement",
        "peak_force_N",
        "force_floor_to_peak_ratio",
        "force_margin_over_absolute_tolerance",
        "force_resolution_qualified",
        "speed_status",
        "barrier_free_from_rest",
        "minimum_available_energy_J",
        "minimum_energy_margin_over_tolerance",
        "barrier_decision_margin_J",
        "barrier_decision_margin_over_tolerance",
        "barrier_resolution_qualified",
        "first_inaccessible_displacement_m",
        "endpoint_available_energy_J",
        "endpoint_positive",
        "endpoint_reachable_from_rest",
        "endpoint_speed_m_s",
        "endpoint_speed_sensitivity_m_s",
        "maximum_reachable_speed_m_s",
        "instantaneous_force_margin_N",
        "barrier_interpretation",
        "status",
        "message",
    ),
    "finite_shell_convergence.csv": (
        "case",
        "periodic_model",
        "effective_far_correction",
        "zero_mode_policy",
        "lower_boundary_model",
        "step_selector",
        "resolved_batch",
        "mesh_id",
        "model_radius_m",
        "radius_source",
        "geometry_radius_m",
        "radius_relative_mismatch",
        "radius_relative_tolerance",
        "target_geometry_representation",
        "path_start_m",
        "path_end_m",
        "image_layers",
        "force_increment_error_N",
        "work_increment_error_J",
        "force_tail_proxy_N",
        "work_tail_proxy_J",
        "converged",
        "reference_model",
        "reference_force_error_N",
        "reference_work_error_J",
        "reference_converged",
        "selected_image_layers",
        "status",
    ),
    "comparison_matrix.csv": (
        "comparison",
        "comparison_kind",
        "sample_kind",
        "resolved_batch",
        "scope",
        "mesh_id",
        "left_snapshot_id",
        "right_snapshot_id",
        "left_effective_far_correction",
        "right_effective_far_correction",
        "metric",
        "left_value",
        "right_value",
        "signed_difference",
        "absolute_difference",
        "relative_difference",
        "status",
        "interpretation",
    ),
    "legacy_estimator_comparison.csv": (
        "comparison_kind",
        "estimator",
        "legacy_self_policy",
        "current_self_policy",
        "batch",
        "mesh_id",
        "displacement_m",
        "component",
        "quantity",
        "legacy_value",
        "current_value",
        "signed_difference",
        "absolute_difference",
        "relative_difference",
        "status",
        "interpretation",
    ),
}


class ValidationError(RuntimeError):
    """Raised when a staged input or output violates the validation contract."""


def _lexical_absolute(path: str | Path) -> Path:
    raw = os.fspath(path)
    value = Path(raw)
    if value.is_absolute() and raw != str(value):
        raise ValidationError(f"path is not lexically canonical: {raw!r}")
    return value if value.is_absolute() else Path.cwd() / value


def _safe_stage_path(path: str | Path, *, label: str) -> Path:
    raw = os.fspath(path)
    if not raw or SAFE_STAGE_PATH_PATTERN.fullmatch(raw) is None:
        raise ValidationError(
            f"unsafe {label} path; allowed characters are [A-Za-z0-9_./:+-]: {raw!r}"
        )
    return Path(raw)


def _canonical_relative_path(relative: str | Path, *, label: str) -> Path:
    raw = os.fspath(relative)
    value = Path(raw)
    if (
        not raw
        or value.is_absolute()
        or "\\" in raw
        or value.as_posix() != raw
        or any(part in {".", ".."} for part in value.parts)
    ):
        raise ValidationError(f"{label} is not a canonical relative path: {raw!r}")
    return value


def _reject_existing_symlink_chain(root: Path, leaf: Path, *, label: str) -> None:
    try:
        relative = leaf.relative_to(root)
    except ValueError as exc:
        raise ValidationError(f"{label} is outside its declared root: {leaf}") from exc
    for depth in range(len(relative.parts) + 1):
        candidate = root.joinpath(*relative.parts[:depth])
        if candidate.is_symlink():
            raise ValidationError(
                f"{label} has a symlink in its root-to-leaf path: {candidate}"
            )


def _reject_existing_symlink_ancestors(path: Path, *, label: str) -> None:
    absolute = _lexical_absolute(path)
    anchor = Path(absolute.anchor)
    relative = absolute.relative_to(anchor)
    for depth in range(len(relative.parts) + 1):
        candidate = anchor.joinpath(*relative.parts[:depth])
        if candidate.is_symlink():
            raise ValidationError(
                f"{label} has an existing symlink ancestor: {candidate}"
            )


def _require_expected_path(
    root: Path,
    declared: str | Path,
    relative: str | Path | None,
    *,
    label: str,
) -> Path:
    root = Path(root)
    if not root.is_absolute() or ".." in root.parts:
        raise ValidationError(f"{label} root is not a canonical absolute path: {root}")
    if relative is None:
        expected = root
    else:
        expected = root / _canonical_relative_path(relative, label=f"{label} key")
    declared_raw = os.fspath(declared)
    actual = Path(declared_raw)
    if declared_raw != str(expected) or actual != expected:
        raise ValidationError(
            f"{label} path is not the canonical expected path: {actual} != {expected}"
        )
    _reject_existing_symlink_chain(root, expected, label=label)
    return expected


def _require_descendant_path(
    root: Path,
    declared: str | Path,
    *,
    prefix: str | Path | None = None,
    label: str,
) -> Path:
    declared_raw = os.fspath(declared)
    actual = Path(declared_raw)
    try:
        relative = actual.relative_to(root)
    except ValueError as exc:
        raise ValidationError(
            f"{label} is outside the validation root: {actual}"
        ) from exc
    if not relative.parts or any(part in {".", ".."} for part in relative.parts):
        raise ValidationError(f"{label} is not a canonical descendant path: {actual}")
    if prefix is not None:
        expected_prefix = _canonical_relative_path(prefix, label=f"{label} prefix")
        if relative.parts[: len(expected_prefix.parts)] != expected_prefix.parts:
            raise ValidationError(
                f"{label} is outside the expected {expected_prefix} subtree: {actual}"
            )
    return _require_expected_path(root, declared_raw, relative, label=label)


def _is_same_or_descendant(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _load_toml(path: Path) -> dict[str, Any]:
    try:
        with path.open("rb") as stream:
            value = tomllib.load(stream)
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise ValidationError(f"failed to read TOML {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValidationError(f"TOML root must be a table: {path}")
    return value


def _write_toml(path: Path, data: Mapping[str, Any]) -> None:
    import tomli_w

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(tomli_w.dumps(dict(data)), encoding="utf-8")


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _write_json_atomic(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(
            json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as exc:
        raise ValidationError(f"failed to hash {path}: {exc}") from exc
    return digest.hexdigest()


def _snapshot_python_source(root: Path, commit: str) -> dict[str, Any]:
    snapshot = root / "source" / commit
    beach_destination = snapshot / "beach"
    tool_destination = snapshot / "tools/periodic_object_validation.py"
    snapshot.mkdir(parents=True, exist_ok=False)
    shutil.copytree(
        REPO_ROOT / "beach",
        beach_destination,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo"),
    )
    for directory in ("src", "app"):
        shutil.copytree(
            REPO_ROOT / directory,
            snapshot / directory,
            ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.i90"),
        )
    for filename in ("fpm.toml", "build.sh", "Makefile", "pyproject.toml"):
        shutil.copy2(REPO_ROOT / filename, snapshot / filename)
    tool_destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(Path(__file__).resolve(), tool_destination)
    files = sorted(path for path in snapshot.rglob("*") if path.is_file())
    hashes = {path.relative_to(snapshot).as_posix(): _sha256(path) for path in files}
    hash_file = snapshot / "source_files.sha256"
    hash_file.write_text(
        "".join(f"{digest}  {relative}\n" for relative, digest in hashes.items()),
        encoding="utf-8",
    )
    for path in sorted(snapshot.rglob("*"), reverse=True):
        path.chmod(0o555 if path.is_dir() else 0o444)
    snapshot.chmod(0o555)
    return {
        "root": str(snapshot),
        "tool": str(tool_destination),
        "hash_file": str(hash_file),
        "hash_file_sha256": _sha256(hash_file),
        "files": hashes,
    }


def _git_provenance() -> dict[str, Any]:
    def run(*args: str) -> str:
        result = subprocess.run(
            ["git", *args],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        return result.stdout.strip() if result.returncode == 0 else ""

    status = run("status", "--porcelain=v1")
    return {
        "repository": str(REPO_ROOT),
        "commit": run("rev-parse", "HEAD") or "unknown",
        "describe": run("describe", "--tags", "--always", "--dirty") or "unknown",
        "dirty": bool(status),
        "status_porcelain": status.splitlines(),
    }


def _parse_build_info(text: str, *, label: str) -> dict[str, Any]:
    fields: dict[str, str] = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        key, separator, value = line.partition("=")
        key = key.strip()
        value = value.strip()
        if not separator or not key or not value or key in fields:
            raise ValidationError(f"{label} build origin payload is malformed")
        fields[key] = value
    if tuple(fields) != BUILD_INFO_KEYS:
        raise ValidationError(
            f"{label} build origin keys differ from the canonical contract"
        )
    try:
        schema_version = int(fields["build_info_schema_version"])
    except ValueError as exc:
        raise ValidationError(f"{label} build origin schema is invalid") from exc
    if schema_version != 1:
        raise ValidationError(f"{label} build origin schema is unsupported")
    return {
        "build_info_schema_version": schema_version,
        **{key: fields[key] for key in BUILD_INFO_KEYS[1:]},
    }


def _binary_build_info(binary: Path) -> dict[str, Any]:
    try:
        result = subprocess.run(
            [str(binary), "--build-info"],
            text=True,
            capture_output=True,
            check=False,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ValidationError(
            f"failed to read executable build origin from {binary}: {exc}"
        ) from exc
    if result.returncode != 0 or result.stderr.strip():
        raise ValidationError(
            "executable build origin query failed: "
            f"exit={result.returncode}, stderr={result.stderr.strip()!r}"
        )
    return _parse_build_info(result.stdout, label="executable")


def _library_build_info(library: Path) -> dict[str, Any]:
    try:
        info = field_kernel_build_info(library)
    except Exception as exc:
        raise ValidationError(
            f"failed to read field-kernel library build origin from {library}: {exc}"
        ) from exc
    return {
        "build_info_schema_version": int(info.schema_version),
        "build_version": str(info.version),
        "build_version_mode": str(info.version_mode),
        "build_source_commit": str(info.source_commit),
        "build_id": str(info.build_id),
    }


def _validate_production_build_origin(
    binary: Path,
    library: Path,
    source_commit: str,
) -> dict[str, Any]:
    binary_info = _binary_build_info(binary)
    library_info = _library_build_info(library)
    if binary_info != library_info:
        raise ValidationError(
            "artifact build origin mismatch between executable and field-kernel library"
        )
    try:
        valid_commit = len(source_commit) == 40 and int(source_commit, 16) >= 0
    except ValueError:
        valid_commit = False
    if (
        not valid_commit
        or binary_info["build_source_commit"] != source_commit
        or binary_info["build_id"] != f"{source_commit}:clean"
        or binary_info["build_version_mode"] not in {"git", "override"}
    ):
        raise ValidationError(
            "artifact build origin does not identify the staged clean source commit"
        )
    return binary_info


def _verify_staged_build_origin(
    manifest: Mapping[str, Any],
    *,
    binary: Path,
    library: Path,
) -> dict[str, Any]:
    expected = manifest.get("build_origin")
    if not isinstance(expected, Mapping) or set(expected) != set(BUILD_INFO_KEYS):
        raise ValidationError("staged build origin contract is missing")
    source_commit = str(manifest.get("source", {}).get("commit", ""))
    actual = _validate_production_build_origin(binary, library, source_commit)
    if actual != dict(expected):
        raise ValidationError("staged artifact build origin no longer matches manifest")
    return actual


def _require_strict_build_origin_contract(manifest: Mapping[str, Any]) -> None:
    execution = manifest.get("execution")
    build_origin = manifest.get("build_origin")
    if (
        not isinstance(execution, Mapping)
        or execution.get("require_clean_source") is not True
        or not isinstance(build_origin, Mapping)
        or set(build_origin) != set(BUILD_INFO_KEYS)
    ):
        raise ValidationError(
            "complete analysis requires a clean build origin contract"
        )


def _archive_version(archive_run: Path) -> str:
    path = archive_run / "work/SIMULATOR_VERSION.txt"
    if not path.exists():
        return "unknown"
    for line in reversed(
        path.read_text(encoding="utf-8", errors="replace").splitlines()
    ):
        value = line.strip()
        if value and not value.startswith(("executable:", "20")):
            return value
    return "unknown"


def _archive_executable_state(archive_run: Path) -> dict[str, Any]:
    manifest_path = archive_run / "manifest.toml"
    result: dict[str, Any] = {
        "declared_path": None,
        "declared_sha256": None,
        "current_sha256": None,
        "exact_executable_available": False,
    }
    if not manifest_path.exists():
        return result
    manifest = _load_toml(manifest_path)
    source = manifest.get("simulator_source", {})
    if not isinstance(source, dict):
        return result
    declared_path = source.get("executable")
    declared_hash = str(source.get("exe_hash", ""))
    if declared_hash.startswith("sha256:"):
        declared_hash = declared_hash.split(":", 1)[1]
    result["declared_path"] = declared_path
    result["declared_sha256"] = declared_hash or None
    if isinstance(declared_path, str) and Path(declared_path).is_file():
        current_hash = _sha256(Path(declared_path))
        result["current_sha256"] = current_hash
        result["exact_executable_available"] = bool(
            declared_hash and current_hash == declared_hash
        )
    return result


def _archive_job_resources(archive_run: Path) -> dict[str, int]:
    manifest_path = archive_run / "manifest.toml"
    if not manifest_path.is_file():
        raise ValidationError(f"archive manifest does not exist: {manifest_path}")
    manifest = _load_toml(manifest_path)
    job = manifest.get("job")
    if not isinstance(job, dict):
        raise ValidationError("archive manifest is missing [job] metadata")
    expected = {"processes": 6, "threads": 112, "cores": 112}
    actual: dict[str, int] = {}
    for key, value in expected.items():
        try:
            actual[key] = int(job[key])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValidationError(
                f"archived MPI resources are missing/invalid: job.{key}"
            ) from exc
        if actual[key] != value:
            raise ValidationError(
                f"archived MPI resources differ: job.{key}={actual[key]} != {value}"
            )
    return actual


def _archive_output_expectations(archive_run: Path) -> dict[str, int | None]:
    summary_path = archive_run / "work/latest/summary.txt"
    if not summary_path.is_file():
        return {"mesh_nelem": None, "mesh_count": None}
    summary = _summary(summary_path)
    try:
        mesh_nelem = int(summary["mesh_nelem"])
        mesh_count = int(summary.get("mesh_count", "0")) or None
    except (KeyError, ValueError) as exc:
        raise ValidationError(
            f"archive summary has invalid mesh metadata: {exc}"
        ) from exc
    if mesh_nelem < 1 or (mesh_count is not None and mesh_count < 1):
        raise ValidationError("archive summary mesh metadata must be positive")
    return {"mesh_nelem": mesh_nelem, "mesh_count": mesh_count}


def _production_mechanics_number(
    values: Mapping[str, Any],
    key: str,
    *,
    positive: bool = False,
    nonnegative: bool = False,
) -> float:
    raw = values.get(key)
    if isinstance(raw, bool) or not isinstance(raw, (int, float)):
        raise ValidationError(f"production mechanics {key} must be numeric")
    number = float(raw)
    if not math.isfinite(number):
        raise ValidationError(f"production mechanics {key} must be finite")
    if positive and number <= 0.0:
        raise ValidationError(f"production mechanics {key} must be positive")
    if nonnegative and number < 0.0:
        raise ValidationError(f"production mechanics {key} must be nonnegative")
    return number


def _validate_production_release_mechanics(archive_run: Path) -> None:
    paths = {
        relative: _require_expected_path(
            archive_run,
            archive_run / relative,
            relative,
            label=f"production mechanics {relative}",
        )
        for relative in PRODUCTION_MECHANICS_INPUTS
    }
    missing = [relative for relative, path in paths.items() if not path.is_file()]
    if missing:
        raise ValidationError(
            "production mechanics file is missing: " + ", ".join(missing)
        )

    summary_path = paths["analysis/local_release/release_model_summary.json"]
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValidationError(
            f"production mechanics release model summary is invalid: {exc}"
        ) from exc
    if not isinstance(summary, dict) or not isinstance(
        summary.get("assumptions"), dict
    ):
        raise ValidationError(
            "production mechanics release model summary requires assumptions"
        )
    assumptions = summary["assumptions"]
    _production_mechanics_number(assumptions, "radius_m", positive=True)
    _production_mechanics_number(assumptions, "dust_density_kg_m3", positive=True)
    _production_mechanics_number(assumptions, "moon_gravity_m_s2", nonnegative=True)
    energy_partition = _production_mechanics_number(assumptions, "energy_partition")
    if not 0.0 <= energy_partition <= 1.0:
        raise ValidationError(
            "production mechanics energy_partition must lie in [0, 1]"
        )

    adhesion_path = paths["input/release_kernel_base.toml"]
    try:
        parsed = _load_toml(adhesion_path)
    except ValidationError as exc:
        raise ValidationError(
            f"production mechanics adhesion TOML is invalid: {exc}"
        ) from exc
    adhesion = parsed.get("adhesion")
    if not isinstance(adhesion, dict) or not isinstance(adhesion.get("model"), str):
        raise ValidationError(
            "production mechanics adhesion TOML requires [adhesion].model"
        )
    model = str(adhesion["model"]).strip().lower()
    if model not in {"none", "vdw_work"}:
        raise ValidationError(
            f"production mechanics adhesion model is unsupported: {model!r}"
        )
    if model == "none":
        return
    values = {
        key: _production_mechanics_number(adhesion, key, positive=True)
        for key in (
            "hamaker_constant",
            "contact_distance",
            "cutoff_distance",
            "roughness_factor",
            "contact_count",
            "peel_factor",
        )
    }
    geometry = adhesion.get("contact_geometry")
    if geometry not in {"sphere_sphere", "sphere_plane"}:
        raise ValidationError(
            "production mechanics contact_geometry must be sphere_sphere or sphere_plane"
        )
    if values["cutoff_distance"] <= values["contact_distance"]:
        raise ValidationError(
            "production mechanics cutoff_distance must exceed contact_distance"
        )


def _archive_analysis_inputs(archive_run: Path) -> dict[str, dict[str, str]]:
    inputs: dict[str, dict[str, str]] = {}
    for relative in (
        *PRODUCTION_MECHANICS_INPUTS,
        *LEGACY_ESTIMATOR_INPUTS,
    ):
        path = _require_expected_path(
            archive_run,
            archive_run / relative,
            relative,
            label=f"archive analysis input {relative}",
        )
        if path.is_file():
            inputs[relative] = {"path": str(path), "sha256": _sha256(path)}
    return inputs


def _archive_analysis_outputs(archive_run: Path) -> dict[str, dict[str, str]]:
    output = archive_run / "work/latest"
    candidates = [
        output / Path(relative).name for relative in ARCHIVE_REQUIRED_ANALYSIS_OUTPUTS
    ]
    candidates.extend(
        path
        for pattern in (
            "charge_ledger.csv",
            "rng_state*.txt",
            "macro_residuals*.csv",
        )
        for path in sorted(output.glob(pattern))
    )
    files: dict[str, dict[str, str]] = {}
    for path in candidates:
        if path.is_file():
            relative = path.relative_to(archive_run).as_posix()
            path = _require_expected_path(
                archive_run,
                path,
                relative,
                label=f"archive analysis output {relative}",
            )
            files[relative] = {"path": str(path), "sha256": _sha256(path)}
    return files


def _verify_declared_files(
    files: Mapping[str, Any],
    *,
    root: Path,
    label: str,
) -> None:
    for name, value in files.items():
        if not isinstance(value, Mapping):
            raise ValidationError(f"{label} metadata is invalid: {name}")
        path = _require_expected_path(
            root,
            str(value.get("path", "")),
            str(name),
            label=f"{label} {name}",
        )
        if not path.is_file() or _sha256(path) != value.get("sha256"):
            raise ValidationError(f"{label} hash mismatch: {name}")


def _require_archive_contract(config: Mapping[str, Any], archive_run: Path) -> None:
    sim = config.get("sim")
    output = config.get("output")
    if not isinstance(sim, dict) or not isinstance(output, dict):
        raise ValidationError("archived input requires [sim] and [output] tables")
    if "field" in config:
        raise ValidationError(
            "archive uses removed point-source configuration table: [field]"
        )
    checks = {
        "sim.batch_count": (sim.get("batch_count"), 280000),
        "sim.rng_seed": (sim.get("rng_seed"), 12345),
        "sim.field_bc_mode": (sim.get("field_bc_mode"), "periodic2"),
        "sim.field_solver": (sim.get("field_solver"), "fmm"),
        "sim.field_normalization": (sim.get("field_normalization", "si"), "si"),
        "sim.field_periodic_image_layers": (
            sim.get("field_periodic_image_layers", 1),
            1,
        ),
        "sim.field_periodic_far_correction": (
            sim.get("field_periodic_far_correction", "none"),
            "none",
        ),
        "sim.bc_x_low": (sim.get("bc_x_low"), "periodic"),
        "sim.bc_x_high": (sim.get("bc_x_high"), "periodic"),
        "sim.bc_y_low": (sim.get("bc_y_low"), "periodic"),
        "sim.bc_y_high": (sim.get("bc_y_high"), "periodic"),
        "sim.bc_z_low": (sim.get("bc_z_low"), "open"),
        "sim.bc_z_high": (sim.get("bc_z_high"), "open"),
    }
    for name, (actual, expected) in checks.items():
        if actual != expected:
            raise ValidationError(
                f"archive {archive_run} has {name}={actual!r}; expected {expected!r}"
            )
    removed_keys = sorted(
        key
        for key in ("softening", "field_source_model", "field_kernel_id")
        if key in sim
    )
    if removed_keys:
        raise ValidationError(
            "archive uses removed point-source configuration keys: "
            + ", ".join(f"sim.{key}" for key in removed_keys)
        )


def _case_config(
    archive: Mapping[str, Any],
    *,
    periodic_model: str,
    batch_count: int,
    history_stride: int,
    output_dir: Path,
    cache_dir: Path,
    resume: bool,
    restart_from: Path | None,
) -> dict[str, Any]:
    config = copy.deepcopy(dict(archive))
    sim = config.setdefault("sim", {})
    output = config.setdefault("output", {})
    if not isinstance(sim, dict) or not isinstance(output, dict):
        raise ValidationError("[sim] and [output] must be TOML tables")
    sim["batch_count"] = int(batch_count)
    sim["field_periodic_image_layers"] = 1
    if periodic_model == "infinite_physical":
        sim["field_periodic_far_correction"] = "cached_kneq0"
        sim["field_periodic_ewald_layers"] = DEFAULT_EWALD_LAYERS
        sim["field_periodic_cache_dir"] = str(cache_dir)
        sim["field_periodic_generation_tolerance"] = DEFAULT_GENERATION_TOLERANCE
    elif periodic_model == "finite_configured":
        sim["field_periodic_far_correction"] = "none"
    else:
        raise ValidationError(f"unknown periodic model: {periodic_model}")
    output["dir"] = str(output_dir)
    output["history_stride"] = int(history_stride)
    output["resume"] = bool(resume)
    if restart_from is None:
        output.pop("restart_from", None)
    else:
        output["restart_from"] = str(restart_from)
    return config


def _case_spec(
    root: Path,
    *,
    name: str,
    config_path: Path,
    output_dir: Path,
    periodic_model: str,
    batch_count: int,
    history_stride: int,
    resume: bool,
    restart_from: Path | None,
    cache_expectation: str,
    role: str,
) -> dict[str, Any]:
    return {
        "name": name,
        "role": role,
        "config_path": str(config_path),
        "output_dir": str(output_dir),
        "periodic_model": periodic_model,
        "batch_count": batch_count,
        "history_stride": history_stride,
        "resume": resume,
        "restart_from": None if restart_from is None else str(restart_from),
        "cache_expectation": cache_expectation,
        "validation_root": str(root),
    }


def _render_job(
    *,
    root: Path,
    binary: Path,
    job_name: str,
    source_commit: str,
    source_snapshot: Mapping[str, Any],
    walltime: str,
    commands: str,
) -> str:
    template = JOB_TEMPLATE.read_text(encoding="utf-8")
    replacements = {
        "@PARTITION@": DEFAULT_PARTITION,
        "@WALLTIME@": walltime,
        "@JOB_NAME@": job_name,
        "@LOG_DIR@": str(root / "submit/logs"),
        "@VALIDATION_ROOT@": str(root),
        "@BINARY@": str(binary),
        "@BINARY_SHA256@": _sha256(binary),
        "@TOOL@": str(source_snapshot["tool"]),
        "@REPO_ROOT@": str(source_snapshot["root"]),
        "@SOURCE_HASH_FILE@": str(source_snapshot["hash_file"]),
        "@SOURCE_HASH_SHA256@": str(source_snapshot["hash_file_sha256"]),
        "@SOURCE_COMMIT@": source_commit,
        "@COMMANDS@": commands,
    }
    for marker, value in replacements.items():
        template = template.replace(marker, value)
    if "@" in "".join(line for line in template.splitlines() if "@" in line):
        unresolved = [line for line in template.splitlines() if "@" in line]
        raise ValidationError(f"unresolved job template markers: {unresolved}")
    return template


def _run_case_commands(
    case_names: Sequence[str],
    cases: Mapping[str, Any],
    *,
    producer_job_role: str,
) -> str:
    if any(
        CASE_PRODUCER_JOB_ROLES.get(case_name) != producer_job_role
        for case_name in case_names
    ):
        raise ValidationError(
            f"producer job role {producer_job_role!r} does not cover {tuple(case_names)!r}"
        )
    lines = [
        "run_case() {",
        '  local name="$1" config="$2" output="$3" batches="$4" config_sha="$5" restart_dir="$6" restart_batches="$7"',
        '  if [ -n "${restart_dir}" ]; then',
        '    python3.11 "${TOOL}" verify-run --case-dir "${restart_dir}" --expected-batches "${restart_batches}" --require-existing-receipt',
        "  fi",
        '  if [ -e "${output}" ] && [ -n "$(find "${output}" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]; then',
        "    printf 'refusing non-empty output: %s\\n' \"${output}\" >&2",
        "    return 2",
        "  fi",
        '  mkdir -p "$(dirname "${output}")"',
        '  if ! printf \'%s  %s\\n\' "${config_sha}" "${config}" | sha256sum --check - | tee -a "${HASH_LOG}"; then',
        "    printf 'input hash mismatch before execution: %s\\n' \"${config}\" >&2",
        "    return 2",
        "  fi",
        '  /usr/bin/time -v -o "${VALIDATION_ROOT}/provenance/time/${name}.txt" srun "${BINARY}" "${config}"',
        '  python3.11 "${TOOL}" verify-run --case-dir "${output}" '
        ' --expected-batches "${batches}"'
        f' --producer-job-role "{producer_job_role}"',
        "}",
    ]
    for case_name in case_names:
        case = cases[case_name]
        restart_dir = ""
        restart_batches = "0"
        if case.get("restart_from") is not None:
            restart_dir = str(case["restart_from"])
            previous = next(
                value for value in cases.values() if value["output_dir"] == restart_dir
            )
            restart_batches = str(previous["batch_count"])
        lines.append(
            "run_case "
            + " ".join(
                [
                    f'"{case_name}"',
                    f'"{case["config_path"]}"',
                    f'"{case["output_dir"]}"',
                    f'"{case["batch_count"]}"',
                    f'"{case["config_sha256"]}"',
                    f'"{restart_dir}"',
                    f'"{restart_batches}"',
                ]
            )
        )
    return "\n".join(lines)


def _write_jobs(
    root: Path,
    archive_path: Path,
    binary: Path,
    cases: Mapping[str, Any],
    *,
    source_commit: str,
    source_snapshot: Mapping[str, Any],
    analysis_library: Mapping[str, Any] | None,
) -> list[Path]:
    submit = root / "submit"
    jobs: list[tuple[str, str, str, str, Sequence[str]]] = [
        (
            "smoke_sysa.sh",
            "beach-periodic-smoke",
            "1:00:00",
            "smoke",
            (
                "cache_prime",
                "smoke_finite_configured",
                "smoke_infinite_physical",
            ),
        ),
        (
            "full_finite_140000_sysa.sh",
            "beach-finite-140k",
            "12:00:00",
            "finite_140000",
            ("full_finite_configured_140000",),
        ),
        (
            "full_finite_280000_sysa.sh",
            "beach-finite-280k",
            "12:00:00",
            "finite_280000",
            ("full_finite_configured_280000",),
        ),
        (
            "full_infinite_140000_sysa.sh",
            "beach-infinite-140k",
            "14:00:00",
            "infinite_140000",
            ("full_infinite_physical_140000",),
        ),
        (
            "full_infinite_280000_sysa.sh",
            "beach-infinite-280k",
            "14:00:00",
            "infinite_280000",
            ("full_infinite_physical_280000",),
        ),
    ]
    paths: list[Path] = []
    for filename, job_name, walltime, producer_job_role, case_names in jobs:
        commands = _run_case_commands(
            case_names,
            cases,
            producer_job_role=producer_job_role,
        )
        if filename == "smoke_sysa.sh" and analysis_library is not None:
            commands = (
                f'python3.11 "{source_snapshot["tool"]}" probe-library '
                f'--library "{analysis_library["staged_path"]}"\n' + commands
            )
        path = submit / filename
        path.write_text(
            _render_job(
                root=root,
                binary=binary,
                job_name=job_name,
                source_commit=source_commit,
                source_snapshot=source_snapshot,
                walltime=walltime,
                commands=commands,
            ),
            encoding="utf-8",
        )
        path.chmod(0o755)
        paths.append(path)
    if analysis_library is not None:
        analysis_template = ANALYSIS_JOB_TEMPLATE.read_text(encoding="utf-8")
        replacements = {
            "@PARTITION@": DEFAULT_PARTITION,
            "@LOG_DIR@": str(root / "submit/logs"),
            "@VALIDATION_ROOT@": str(root),
            "@ARCHIVE_RUN@": str(archive_path),
            "@TOOL@": str(source_snapshot["tool"]),
            "@REPO_ROOT@": str(source_snapshot["root"]),
            "@SOURCE_COMMIT@": source_commit,
            "@SOURCE_HASH_FILE@": str(source_snapshot["hash_file"]),
            "@SOURCE_HASH_SHA256@": str(source_snapshot["hash_file_sha256"]),
            "@LIBRARY@": str(analysis_library["staged_path"]),
            "@LIBRARY_SHA256@": str(analysis_library["sha256"]),
        }
        for marker, value in replacements.items():
            analysis_template = analysis_template.replace(marker, value)
        unresolved = [line for line in analysis_template.splitlines() if "@" in line]
        if unresolved:
            raise ValidationError(
                f"unresolved analysis job template markers: {unresolved}"
            )
        analysis_path = submit / "analysis_sysa.sh"
        analysis_path.write_text(analysis_template, encoding="utf-8")
        analysis_path.chmod(0o755)
        paths.append(analysis_path)
    chain = submit / "submit_chain.sh"
    chain_text = """#!/bin/bash
set -euo pipefail
current_sys="$(module -t list 2>&1 | grep -E '^Sys(A|B|C|CL|G)/' | head -n 1 || true)"
if [ -n "${current_sys}" ] && [ "${current_sys}" != "SysA/2022" ]; then
  module switch "${current_sys}" SysA/2022
elif [ -z "${current_sys}" ]; then
  module load SysA/2022
fi
module list 2>&1 | grep -q 'SysA/2022' || {
  printf '%s\n' 'SysA/2022 module is required for gr20001a submission.' >&2
  exit 2
}
spartition | awk '$1 == "gr20001a" && $2 == "UP" { found=1 } END { exit !found }' || {
  printf '%s\n' 'gr20001a is not visible and UP.' >&2
  exit 2
}
if command -v qgroup >/dev/null 2>&1; then
  qgroup
fi
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VALIDATION_ROOT="@VALIDATION_ROOT@"
TOOL="@TOOL@"
SOURCE_ROOT="@REPO_ROOT@"
unset PYTHONHOME
export PYTHONNOUSERSITE=1
export PYTHONPATH="${SOURCE_ROOT}"
export PYTHONDONTWRITEBYTECODE=1
JOB_IDS="${HERE}/job_ids.json"
SUBMIT_LOCK="${HERE}/.submit_chain.lock"
if [ -e "${JOB_IDS}" ] || [ -L "${JOB_IDS}" ]; then
  printf 'refusing to resubmit an existing validation chain: %s\n' "${JOB_IDS}" >&2
  exit 2
fi
if ! mkdir "${SUBMIT_LOCK}" 2>/dev/null; then
  printf 'another validation-chain submission is in progress: %s\n' "${SUBMIT_LOCK}" >&2
  exit 2
fi
cleanup_lock() {
  rmdir "${SUBMIT_LOCK}" 2>/dev/null || true
}
trap cleanup_lock EXIT
python3.11 "${TOOL}" verify-inputs --validation-root "${VALIDATION_ROOT}"
smoke=""
finite_140=""
infinite_140=""
finite_280=""
infinite_280=""
analysis=""
submission_complete=false
write_job_ids() {
  local temporary="${HERE}/job_ids.json.tmp.$$"
  printf '{\n  "submission_complete": %s,\n  "smoke": "%s",\n  "finite_140000": "%s",\n  "finite_280000": "%s",\n  "infinite_140000": "%s",\n  "infinite_280000": "%s",\n  "analysis": "%s"\n}\n' \
    "${submission_complete}" "${smoke}" "${finite_140}" "${finite_280}" \
    "${infinite_140}" "${infinite_280}" "${analysis}" > "${temporary}"
  mv "${temporary}" "${JOB_IDS}"
}
finish_submission() {
  local exit_code=$?
  write_job_ids || true
  cleanup_lock
  return "${exit_code}"
}
trap finish_submission EXIT
write_job_ids
smoke=$(sbatch --parsable "${HERE}/smoke_sysa.sh"); smoke=${smoke%%;*}
write_job_ids
finite_140=$(sbatch --parsable --dependency=afterok:${smoke} "${HERE}/full_finite_140000_sysa.sh"); finite_140=${finite_140%%;*}
write_job_ids
infinite_140=$(sbatch --parsable --dependency=afterok:${smoke} "${HERE}/full_infinite_140000_sysa.sh"); infinite_140=${infinite_140%%;*}
write_job_ids
finite_280=$(sbatch --parsable --dependency=afterok:${finite_140} "${HERE}/full_finite_280000_sysa.sh"); finite_280=${finite_280%%;*}
write_job_ids
infinite_280=$(sbatch --parsable --dependency=afterok:${infinite_140} "${HERE}/full_infinite_280000_sysa.sh"); infinite_280=${infinite_280%%;*}
write_job_ids
if [ -f "${HERE}/analysis_sysa.sh" ]; then
  analysis=$(sbatch --parsable --dependency=afterok:${finite_280}:${infinite_280} "${HERE}/analysis_sysa.sh"); analysis=${analysis%%;*}
  write_job_ids
fi
submission_complete=true
write_job_ids
trap - EXIT
cleanup_lock
printf 'smoke=%s finite_140=%s finite_280=%s infinite_140=%s infinite_280=%s analysis=%s\n' \
  "${smoke}" "${finite_140}" "${finite_280}" "${infinite_140}" "${infinite_280}" "${analysis}"
"""
    chain.write_text(
        chain_text.replace("@VALIDATION_ROOT@", str(root))
        .replace("@TOOL@", str(source_snapshot["tool"]))
        .replace("@REPO_ROOT@", str(source_snapshot["root"])),
        encoding="utf-8",
    )
    chain.chmod(0o755)
    paths.append(chain)
    return paths


def stage_validation(
    archive_run: str | Path,
    validation_root: str | Path,
    binary: str | Path,
    library: str | Path | None = None,
    require_clean_source: bool = False,
) -> dict[str, Any]:
    """Create deterministic configs, provenance, and SysA job scripts."""

    archive_path = _safe_stage_path(archive_run, label="archive run").resolve()
    raw_root = _safe_stage_path(validation_root, label="validation root")
    _reject_existing_symlink_ancestors(raw_root, label="validation root")
    root = raw_root.resolve()
    binary_path = _safe_stage_path(binary, label="binary").resolve()
    library_path = (
        None
        if library is None
        else _safe_stage_path(library, label="analysis library").resolve()
    )
    if _is_same_or_descendant(root, REPO_ROOT.resolve()) or _is_same_or_descendant(
        root, archive_path
    ):
        raise ValidationError(
            f"validation root must be outside the repository and archive: {root}"
        )
    archive_input = _require_expected_path(
        archive_path,
        archive_path / "input/beach.toml",
        "input/beach.toml",
        label="archived input",
    )
    if not archive_input.is_file():
        raise ValidationError(f"archived input does not exist: {archive_input}")
    if not binary_path.is_file():
        raise ValidationError(f"binary does not exist: {binary_path}")
    if library_path is not None and not library_path.is_file():
        raise ValidationError(f"field-kernel library does not exist: {library_path}")
    if root.exists() and any(root.iterdir()):
        raise ValidationError(f"validation root must be empty before staging: {root}")
    git = _git_provenance()
    commit = str(git["commit"])
    if require_clean_source and bool(git.get("dirty")):
        raise ValidationError(
            "production staging requires a clean source worktree: "
            + "; ".join(str(value) for value in git.get("status_porcelain", []))
        )
    if require_clean_source and library_path is None:
        raise ValidationError(
            "production staging requires an analysis library built with the executable"
        )
    build_origin: dict[str, Any] | None = None
    if require_clean_source:
        assert library_path is not None
        build_origin = _validate_production_build_origin(
            binary_path,
            library_path,
            commit,
        )
    archive_config = _load_toml(archive_input)
    _require_archive_contract(archive_config, archive_path)
    if require_clean_source:
        _validate_production_release_mechanics(archive_path)
    archive_resources = _archive_job_resources(archive_path)
    archive_output = _archive_output_expectations(archive_path)
    archive_analysis_inputs = _archive_analysis_inputs(archive_path)
    archive_analysis_outputs = _archive_analysis_outputs(archive_path)
    archive_mesh_source_contract = _mesh_source_contract(
        archive_path / "work/latest/mesh_sources.csv",
        expected_meshes=archive_output.get("mesh_count"),
        expected_elements=archive_output.get("mesh_nelem"),
    )
    if (root / "manifest.json").exists():
        raise ValidationError(f"validation root is already staged: {root}")

    for relative in (
        "archive",
        "bin",
        "cache/periodic2",
        "cache/oracles",
        "input/source",
        "input/cache_prime",
        "input/smoke",
        "input/full/segments",
        "run/cache_prime",
        "run/smoke",
        "run/full/finite_configured",
        "run/full/infinite_physical",
        "submit/logs",
        "provenance/time",
        "provenance/oracles",
        "analysis",
    ):
        (root / relative).mkdir(parents=True, exist_ok=True)

    source_copy = root / "input/source/archive_beach.toml"
    shutil.copy2(archive_input, source_copy)
    staged_binary = root / "bin" / commit / "beach"
    staged_binary.parent.mkdir(parents=True, exist_ok=True)
    if binary_path != staged_binary:
        shutil.copy2(binary_path, staged_binary)
    staged_binary.chmod(0o555)
    analysis_library: dict[str, Any] | None = None
    if library_path is not None:
        staged_library = root / "lib" / commit / library_path.name
        staged_library.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(library_path, staged_library)
        staged_library.chmod(0o555)
        analysis_library = {
            "provided_path": str(library_path),
            "staged_path": str(staged_library),
            "sha256": _sha256(staged_library),
        }

    cache_dir = root / "cache/periodic2"
    case_definitions = [
        (
            "cache_prime",
            root / "input/cache_prime/infinite_physical.toml",
            root / "run/cache_prime/infinite_physical",
            "infinite_physical",
            1,
            1,
            False,
            None,
            "miss",
            "cache_prime",
        ),
        (
            "smoke_finite_configured",
            root / "input/smoke/finite_configured.toml",
            root / "run/smoke/finite_configured",
            "finite_configured",
            100,
            10,
            False,
            None,
            "none",
            "smoke",
        ),
        (
            "smoke_infinite_physical",
            root / "input/smoke/infinite_physical.toml",
            root / "run/smoke/infinite_physical",
            "infinite_physical",
            100,
            10,
            False,
            None,
            "hit",
            "smoke",
        ),
        (
            "full_finite_configured",
            root / "input/full/finite_configured.toml",
            root / "run/full/finite_configured/canonical",
            "finite_configured",
            280000,
            1000,
            False,
            None,
            "none",
            "canonical",
        ),
        (
            "full_infinite_physical",
            root / "input/full/infinite_physical.toml",
            root / "run/full/infinite_physical/canonical",
            "infinite_physical",
            280000,
            1000,
            False,
            None,
            "hit",
            "canonical",
        ),
        (
            "full_finite_configured_140000",
            root / "input/full/segments/finite_configured_140000.toml",
            root / "run/full/finite_configured/140000",
            "finite_configured",
            140000,
            1000,
            False,
            None,
            "none",
            "segment",
        ),
        (
            "full_finite_configured_280000",
            root / "input/full/segments/finite_configured_280000.toml",
            root / "run/full/finite_configured/280000",
            "finite_configured",
            280000,
            1000,
            True,
            root / "run/full/finite_configured/140000",
            "none",
            "segment",
        ),
        (
            "full_infinite_physical_140000",
            root / "input/full/segments/infinite_physical_140000.toml",
            root / "run/full/infinite_physical/140000",
            "infinite_physical",
            140000,
            1000,
            False,
            None,
            "hit",
            "segment",
        ),
        (
            "full_infinite_physical_280000",
            root / "input/full/segments/infinite_physical_280000.toml",
            root / "run/full/infinite_physical/280000",
            "infinite_physical",
            280000,
            1000,
            True,
            root / "run/full/infinite_physical/140000",
            "hit",
            "segment",
        ),
    ]
    cases: dict[str, dict[str, Any]] = {}
    for definition in case_definitions:
        (
            name,
            config_path,
            output_dir,
            periodic_model,
            batch_count,
            history_stride,
            resume,
            restart_from,
            cache_expectation,
            role,
        ) = definition
        config = _case_config(
            archive_config,
            periodic_model=periodic_model,
            batch_count=batch_count,
            history_stride=history_stride,
            output_dir=output_dir,
            cache_dir=cache_dir,
            resume=resume,
            restart_from=restart_from,
        )
        _write_toml(config_path, config)
        spec = _case_spec(
            root,
            name=name,
            config_path=config_path,
            output_dir=output_dir,
            periodic_model=periodic_model,
            batch_count=batch_count,
            history_stride=history_stride,
            resume=resume,
            restart_from=restart_from,
            cache_expectation=cache_expectation,
            role=role,
        )
        spec["config_sha256"] = _sha256(config_path)
        cases[name] = spec

    source_snapshot = _snapshot_python_source(root, commit)
    jobs = _write_jobs(
        root,
        archive_path,
        staged_binary,
        cases,
        source_commit=commit,
        source_snapshot=source_snapshot,
        analysis_library=analysis_library,
    )
    executable_state = _archive_executable_state(archive_path)
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "created_at": _utc_now(),
        "validation_root": str(root),
        "archive": {
            "run": str(archive_path),
            "input": str(archive_input),
            "input_sha256": _sha256(archive_input),
            "source_copy": str(source_copy),
            "version": _archive_version(archive_path),
            "manifest_sha256": _sha256(archive_path / "manifest.toml"),
            "job_resources": archive_resources,
            "expected_output": archive_output,
            "analysis_inputs": archive_analysis_inputs,
            "analysis_outputs": archive_analysis_outputs,
            "mesh_source_contract": archive_mesh_source_contract,
            **executable_state,
            "provenance_limitation": (
                "archived output spans parent/child executables; exact archived "
                "binary is not assumed available"
            ),
        },
        "source": git,
        "build_origin": build_origin,
        "source_snapshot": source_snapshot,
        "binary": {
            "provided_path": str(binary_path),
            "staged_path": str(staged_binary),
            "sha256": _sha256(staged_binary),
        },
        "analysis_library": analysis_library,
        "resources": dict(EXPECTED_RESOURCES),
        "modules": [
            part for part in os.environ.get("LOADEDMODULES", "").split(":") if part
        ],
        "allowed_differences": list(ALLOWED_DIFFERENCES),
        "execution": {
            "require_clean_source": require_clean_source,
            "fresh_resume": False,
            "segment_batches": [140000, 280000],
            "immutable_segment_outputs": True,
            "canonical_finite_representation": (
                "archived physics unchanged; output path and strict fresh-start "
                "resume=false are operational differences"
            ),
            "restart_policy": (
                "280000 segments read verified 140000 checkpoints through "
                "output.restart_from and write separate output directories"
            ),
        },
        "receipt_policy": dict(RECEIPT_POLICY),
        "cases": cases,
        "scripts": {
            path.name: {"path": str(path), "sha256": _sha256(path)} for path in jobs
        },
    }
    _write_json(root / "manifest.json", manifest)
    _write_json(root / "provenance/source.json", git)
    _write_json(
        root / "archive/source.json",
        {"archive": manifest["archive"], "interpretation": manifest["execution"]},
    )
    return manifest


def _normalize_config(config: Mapping[str, Any]) -> dict[str, Any]:
    normalized = copy.deepcopy(dict(config))
    sim = normalized.setdefault("sim", {})
    if isinstance(sim, dict):
        sim.setdefault("field_periodic_ewald_layers", DEFAULT_EWALD_LAYERS)
    return normalized


def _first_difference(left: Any, right: Any, path: str = "") -> str | None:
    if isinstance(left, dict) and isinstance(right, dict):
        keys = sorted(set(left) | set(right))
        for key in keys:
            child = f"{path}.{key}" if path else str(key)
            if key not in left:
                return f"{child}: unexpected key"
            if key not in right:
                return f"{child}: missing key"
            difference = _first_difference(left[key], right[key], child)
            if difference is not None:
                return difference
        return None
    if isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            return f"{path}: length {len(left)} != {len(right)}"
        for index, (left_item, right_item) in enumerate(zip(left, right)):
            difference = _first_difference(left_item, right_item, f"{path}[{index}]")
            if difference is not None:
                return difference
        return None
    if type(left) is not type(right) or left != right:
        return f"{path}: {left!r} != {right!r}"
    return None


def _load_manifest(root: Path) -> dict[str, Any]:
    path = _require_expected_path(
        root,
        root / "manifest.json",
        "manifest.json",
        label="validation manifest",
    )
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValidationError(
            f"failed to read validation manifest {path}: {exc}"
        ) from exc
    if not isinstance(value, dict):
        raise ValidationError("validation manifest root must be an object")
    return value


def _validate_case_graph(root: Path, cases: Mapping[str, Any]) -> None:
    if set(cases) != set(EXPECTED_CASE_GRAPH):
        raise ValidationError(
            "staged case graph differs from the fixed validation contract"
        )
    config_paths: set[Path] = set()
    output_paths: set[Path] = set()
    for name, expected in EXPECTED_CASE_GRAPH.items():
        value = cases.get(name)
        if not isinstance(value, Mapping):
            raise ValidationError(f"staged case graph mismatch for {name}")
        expected_config = root / str(expected["config_relative"])
        expected_output = root / str(expected["output_relative"])
        restart_case = expected["restart_case"]
        expected_restart = (
            None
            if restart_case is None
            else root / str(EXPECTED_CASE_GRAPH[str(restart_case)]["output_relative"])
        )
        try:
            actual_config = _require_expected_path(
                root,
                str(value["config_path"]),
                str(expected["config_relative"]),
                label=f"case {name} config",
            )
            actual_output = _require_expected_path(
                root,
                str(value["output_dir"]),
                str(expected["output_relative"]),
                label=f"case {name} output",
            )
            actual_restart = (
                None
                if value.get("restart_from") is None
                else _require_expected_path(
                    root,
                    str(value["restart_from"]),
                    str(EXPECTED_CASE_GRAPH[str(restart_case)]["output_relative"]),
                    label=f"case {name} restart",
                )
            )
            _require_expected_path(
                root,
                str(value.get("validation_root", "")),
                None,
                label=f"case {name} validation root",
            )
            matches = (
                value.get("name") == name
                and value.get("role") == expected["role"]
                and value.get("periodic_model") == expected["periodic_model"]
                and int(value.get("batch_count", -1)) == int(expected["batch_count"])
                and int(value.get("history_stride", -1))
                == int(expected["history_stride"])
                and value.get("resume") is expected["resume"]
                and value.get("cache_expectation") == expected["cache_expectation"]
                and actual_config == expected_config
                and actual_output == expected_output
                and actual_restart == expected_restart
            )
        except (KeyError, TypeError, ValueError, OSError, ValidationError) as exc:
            raise ValidationError(
                f"staged case graph mismatch for {name}: {exc}"
            ) from exc
        if not matches:
            raise ValidationError(f"staged case graph mismatch for {name}")
        if actual_config in config_paths or actual_output in output_paths:
            raise ValidationError(f"staged case graph path alias for {name}")
        config_paths.add(actual_config)
        output_paths.add(actual_output)


def verify_inputs(
    validation_root: str | Path,
    *,
    require_empty_outputs: bool = True,
) -> dict[str, Any]:
    """Fail closed unless every staged input matches its declared transformation."""

    root = _lexical_absolute(validation_root)
    manifest = _load_manifest(root)
    _require_expected_path(
        root,
        str(manifest.get("validation_root", "")),
        None,
        label="validation root",
    )
    for relative in (
        "cache/periodic2",
        "cache/oracles",
        "run",
        "submit",
        "provenance",
        "provenance/time",
        "provenance/oracles",
        "analysis",
    ):
        _require_expected_path(
            root,
            root / relative,
            relative,
            label=f"validation {relative}",
        )
    if manifest.get("resources") != EXPECTED_RESOURCES:
        raise ValidationError(
            f"MPI resource metadata mismatch: {manifest.get('resources')!r}"
        )
    if manifest.get("receipt_policy") != RECEIPT_POLICY:
        raise ValidationError(
            "receipt policy differs from the fixed validation contract"
        )
    source_snapshot = manifest.get("source_snapshot")
    if not isinstance(source_snapshot, dict):
        raise ValidationError("source snapshot metadata is missing")
    snapshot_root = Path(str(source_snapshot.get("root", "")))
    source_commit = str(manifest.get("source", {}).get("commit", ""))
    snapshot_root = _require_expected_path(
        root,
        snapshot_root,
        f"source/{source_commit}",
        label="source snapshot root",
    )
    _require_expected_path(
        root,
        str(source_snapshot.get("tool", "")),
        f"source/{source_commit}/tools/periodic_object_validation.py",
        label="source snapshot tool",
    )
    source_files = source_snapshot.get("files")
    if not isinstance(source_files, dict) or not source_files:
        raise ValidationError("source snapshot file hashes are missing")
    actual_source_files: set[str] = set()
    for path in snapshot_root.rglob("*"):
        if path.is_symlink():
            raise ValidationError(f"source snapshot contains a symlink: {path}")
        if path.is_file() and path.name != "source_files.sha256":
            actual_source_files.add(path.relative_to(snapshot_root).as_posix())
    if actual_source_files != set(source_files):
        missing = sorted(set(source_files) - actual_source_files)
        extra = sorted(actual_source_files - set(source_files))
        raise ValidationError(
            f"source snapshot inventory mismatch: missing={missing}, extra={extra}"
        )
    for relative, expected_hash in source_files.items():
        path = snapshot_root / str(relative)
        if not path.is_file() or _sha256(path) != expected_hash:
            raise ValidationError(f"source snapshot hash mismatch: {relative}")
    hash_file = _require_expected_path(
        root,
        str(source_snapshot.get("hash_file", "")),
        f"source/{source_commit}/source_files.sha256",
        label="source snapshot hash manifest",
    )
    if not hash_file.is_file() or _sha256(hash_file) != source_snapshot.get(
        "hash_file_sha256"
    ):
        raise ValidationError("source snapshot hash manifest is missing")
    archive = manifest.get("archive", {})
    if not isinstance(archive, Mapping):
        raise ValidationError("archive metadata is invalid")
    archive_root = Path(str(archive.get("run", "")))
    _require_expected_path(
        archive_root,
        archive_root,
        None,
        label="archive run",
    )
    archive_input = _require_expected_path(
        archive_root,
        str(archive.get("input", "")),
        "input/beach.toml",
        label="archived input",
    )
    if not archive_input.is_file() or _sha256(archive_input) != archive.get(
        "input_sha256"
    ):
        raise ValidationError("archived input hash mismatch")
    _require_expected_path(
        root,
        str(archive.get("source_copy", "")),
        "input/source/archive_beach.toml",
        label="archived input source copy",
    )
    archive_manifest = _require_expected_path(
        archive_root,
        archive_root / "manifest.toml",
        "manifest.toml",
        label="archive manifest",
    )
    if not archive_manifest.is_file() or _sha256(archive_manifest) != archive.get(
        "manifest_sha256"
    ):
        raise ValidationError("archive manifest hash mismatch")
    if archive.get("job_resources") != {"processes": 6, "threads": 112, "cores": 112}:
        raise ValidationError("archived MPI resources differ from the staged contract")
    analysis_inputs = archive.get("analysis_inputs", {})
    if not isinstance(analysis_inputs, dict):
        raise ValidationError("archive analysis input metadata is invalid")
    _verify_declared_files(
        analysis_inputs,
        root=archive_root,
        label="archive analysis input",
    )
    analysis_outputs = archive.get("analysis_outputs", {})
    if not isinstance(analysis_outputs, dict):
        raise ValidationError("archive analysis output metadata is invalid")
    _verify_declared_files(
        analysis_outputs,
        root=archive_root,
        label="archive analysis output",
    )
    mesh_source_contract = archive.get("mesh_source_contract")
    if not isinstance(mesh_source_contract, Mapping):
        raise ValidationError("archive mesh source contract is missing")
    declared_mesh_source_path = _require_expected_path(
        archive_root,
        str(mesh_source_contract.get("path", "")),
        "work/latest/mesh_sources.csv",
        label="archive mesh source contract",
    )
    expected_output = archive.get("expected_output", {})
    if not isinstance(expected_output, Mapping):
        raise ValidationError("archive output contract is invalid")
    current_mesh_source_contract = _mesh_source_contract(
        declared_mesh_source_path,
        expected_meshes=expected_output.get("mesh_count"),
        expected_elements=expected_output.get("mesh_nelem"),
    )
    if current_mesh_source_contract != dict(mesh_source_contract):
        raise ValidationError("archive mesh source contract no longer matches staging")
    staged_binary = _require_expected_path(
        root,
        str(manifest.get("binary", {}).get("staged_path", "")),
        f"bin/{source_commit}/beach",
        label="staged binary",
    )
    if not staged_binary.is_file() or _sha256(staged_binary) != manifest.get(
        "binary", {}
    ).get("sha256"):
        raise ValidationError("staged binary hash mismatch")
    analysis_library = manifest.get("analysis_library")
    staged_library: Path | None = None
    if analysis_library is not None:
        if not isinstance(analysis_library, dict):
            raise ValidationError("analysis library metadata is invalid")
        staged_library = _require_descendant_path(
            root,
            str(analysis_library.get("staged_path", "")),
            prefix=f"lib/{source_commit}",
            label="staged analysis library",
        )
        if not staged_library.is_file() or _sha256(
            staged_library
        ) != analysis_library.get("sha256"):
            raise ValidationError("staged analysis library hash mismatch")
    requires_clean_source = bool(
        manifest.get("execution", {}).get("require_clean_source", False)
    )
    if requires_clean_source and manifest.get("build_origin") is None:
        raise ValidationError("production manifest has no build origin contract")
    if manifest.get("build_origin") is not None:
        if staged_library is None:
            raise ValidationError("build origin verification requires a staged library")
        _verify_staged_build_origin(
            manifest,
            binary=staged_binary,
            library=staged_library,
        )
    if requires_clean_source:
        missing_mechanics = sorted(
            set(PRODUCTION_MECHANICS_INPUTS) - set(analysis_inputs)
        )
        if missing_mechanics:
            raise ValidationError(
                "production mechanics metadata is missing: "
                + ", ".join(missing_mechanics)
            )
        _validate_production_release_mechanics(archive_root)
    archive_config = _load_toml(archive_input)
    cache_dir = root / "cache/periodic2"
    cases = manifest.get("cases")
    if not isinstance(cases, dict):
        raise ValidationError("manifest cases must be an object")
    _validate_case_graph(root, cases)
    for name, value in sorted(cases.items()):
        if not isinstance(value, dict):
            raise ValidationError(f"case {name} metadata must be an object")
        path = Path(str(value["config_path"]))
        actual = _load_toml(path)
        restart_value = value.get("restart_from")
        expected = _case_config(
            archive_config,
            periodic_model=str(value["periodic_model"]),
            batch_count=int(value["batch_count"]),
            history_stride=int(value["history_stride"]),
            output_dir=Path(str(value["output_dir"])),
            cache_dir=cache_dir,
            resume=bool(value["resume"]),
            restart_from=None if restart_value is None else Path(str(restart_value)),
        )
        difference = _first_difference(
            _normalize_config(actual), _normalize_config(expected)
        )
        if difference is not None:
            raise ValidationError(f"input mismatch for {name}: {difference}")
        if _sha256(path) != value.get("config_sha256"):
            raise ValidationError(f"input hash mismatch for {name}")
        output_dir = Path(str(value["output_dir"]))
        if require_empty_outputs and output_dir.exists() and any(output_dir.iterdir()):
            raise ValidationError(f"fresh output is not empty for {name}: {output_dir}")
    forbidden = ("--ntasks", "--cpus-per-task", "mpiexec", "mpirun", "srun -n")
    scripts = manifest.get("scripts")
    if not isinstance(scripts, Mapping):
        raise ValidationError("staged script metadata is invalid")
    expected_scripts = {
        "smoke_sysa.sh",
        "full_finite_140000_sysa.sh",
        "full_finite_280000_sysa.sh",
        "full_infinite_140000_sysa.sh",
        "full_infinite_280000_sysa.sh",
        "submit_chain.sh",
    }
    if analysis_library is not None:
        expected_scripts.add("analysis_sysa.sh")
    if set(scripts) != expected_scripts:
        raise ValidationError("staged script set differs from the fixed job graph")
    case_coverage = {
        "smoke_sysa.sh": (
            "cache_prime",
            "smoke_finite_configured",
            "smoke_infinite_physical",
        ),
        "full_finite_140000_sysa.sh": ("full_finite_configured_140000",),
        "full_finite_280000_sysa.sh": ("full_finite_configured_280000",),
        "full_infinite_140000_sysa.sh": ("full_infinite_physical_140000",),
        "full_infinite_280000_sysa.sh": ("full_infinite_physical_280000",),
    }
    script_producer_roles = {
        "smoke_sysa.sh": "smoke",
        "full_finite_140000_sysa.sh": "finite_140000",
        "full_finite_280000_sysa.sh": "finite_280000",
        "full_infinite_140000_sysa.sh": "infinite_140000",
        "full_infinite_280000_sysa.sh": "infinite_280000",
    }
    for name, script in scripts.items():
        if not isinstance(script, Mapping):
            raise ValidationError(f"job script metadata is invalid: {name}")
        path = _require_expected_path(
            root,
            str(script["path"]),
            f"submit/{name}",
            label=f"job script {name}",
        )
        if not path.is_file() or _sha256(path) != script.get("sha256"):
            raise ValidationError(f"job script hash mismatch: {name}")
        text = path.read_text(encoding="utf-8")
        required_python_environment = (
            "unset PYTHONHOME",
            "export PYTHONNOUSERSITE=1",
            'export PYTHONPATH="${SOURCE_ROOT}"',
        )
        if any(token not in text for token in required_python_environment) or (
            "${PYTHONPATH:+" in text
        ):
            raise ValidationError(
                f"job script {name} lacks the isolated Python environment contract"
            )
        if name == "analysis_sysa.sh":
            bad = [token for token in forbidden if token in text]
            if bad:
                raise ValidationError(
                    f"job script {name} uses forbidden launch tokens: {bad}"
                )
            required = (
                "#SBATCH --rsc p=1:t=28:c=28",
                "srun python3.11",
                "probe-periodic-oracles",
                "--require-complete",
            )
            if any(token not in text for token in required):
                raise ValidationError(
                    f"job script {name} lacks the strict SysA analysis contract"
                )
            if text.index("probe-periodic-oracles") > text.index(" analyze "):
                raise ValidationError(
                    f"job script {name} runs strict analysis before its plane oracles"
                )
        elif name.endswith("_sysa.sh"):
            bad = [token for token in forbidden if token in text]
            if bad:
                raise ValidationError(
                    f"job script {name} uses forbidden launch tokens: {bad}"
                )
            if (
                "#SBATCH --rsc p=6:t=112:c=112" not in text
                or 'srun "${BINARY}"' not in text
            ):
                raise ValidationError(
                    f"job script {name} lacks the SysA direct-srun contract"
                )
            for case_name in case_coverage[name]:
                if f'run_case "{case_name}"' not in text:
                    raise ValidationError(
                        f"job script {name} omits staged case {case_name}"
                    )
            producer_role = script_producer_roles[name]
            self_verify = [
                line
                for line in text.splitlines()
                if "verify-run" in line and '--expected-batches "${batches}"' in line
            ]
            if len(self_verify) != 1 or (
                f'--producer-job-role "{producer_role}"' not in self_verify[0]
            ):
                raise ValidationError(
                    f"job script {name} has an invalid producer job role"
                )
            parent_verify = [
                line
                for line in text.splitlines()
                if "verify-run" in line and "--require-existing-receipt" in line
            ]
            if any("--producer-job-role" in line for line in parent_verify):
                raise ValidationError(
                    f"job script {name} binds a parent receipt to the child producer job role"
                )
        elif name == "submit_chain.sh":
            required_edges = (
                "finite_140=$(sbatch --parsable --dependency=afterok:${smoke}",
                "infinite_140=$(sbatch --parsable --dependency=afterok:${smoke}",
                "finite_280=$(sbatch --parsable --dependency=afterok:${finite_140}",
                "infinite_280=$(sbatch --parsable --dependency=afterok:${infinite_140}",
                "--dependency=afterok:${finite_280}:${infinite_280}",
                "refusing to resubmit an existing validation chain",
            )
            if any(token not in text for token in required_edges):
                raise ValidationError("submit_chain.sh differs from the fixed job DAG")
    return {"status": "ok", "case_count": len(cases), "validation_root": str(root)}


def _find_validation_root(case_dir: Path) -> Path:
    for candidate in (case_dir, *case_dir.parents):
        if (candidate / "manifest.json").is_file():
            return candidate
    raise ValidationError(
        f"cannot locate manifest.json above case directory: {case_dir}"
    )


def _summary(path: Path) -> dict[str, str]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise ValidationError(f"failed to read {path}: {exc}") from exc
    result: dict[str, str] = {}
    for line in lines:
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        normalized_key = key.strip()
        if normalized_key in result:
            raise ValidationError(
                f"summary contains duplicate key {normalized_key!r}: {path}"
            )
        result[normalized_key] = value.strip()
    return result


def _summary_int(summary: Mapping[str, str], key: str) -> int:
    try:
        value = int(summary[key])
    except (KeyError, ValueError) as exc:
        raise ValidationError(f"summary missing/invalid integer {key}") from exc
    if value < 0:
        raise ValidationError(f"summary {key} must be nonnegative")
    return value


def _summary_float(summary: Mapping[str, str], key: str) -> float:
    try:
        value = float(summary[key])
    except (KeyError, ValueError) as exc:
        raise ValidationError(f"summary missing/invalid float {key}") from exc
    if not math.isfinite(value):
        raise ValidationError(f"summary {key} must be finite")
    return value


def _summary_bool(summary: Mapping[str, str], key: str) -> bool:
    value = summary.get(key, "").strip().lower()
    if value in {"t", "true"}:
        return True
    if value in {"f", "false"}:
        return False
    raise ValidationError(f"summary missing/invalid boolean {key}")


def _summary_build_origin(summary: Mapping[str, str]) -> dict[str, Any]:
    return {
        "build_info_schema_version": _summary_int(summary, "build_info_schema_version"),
        **{key: summary.get(key, "") for key in BUILD_INFO_KEYS[1:]},
    }


def _validate_charges(path: Path, expected_rows: int) -> None:
    try:
        with path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None or not {"elem_idx", "charge_C"}.issubset(
                reader.fieldnames
            ):
                raise ValidationError("charges.csv is missing required columns")
            seen: set[int] = set()
            for row in reader:
                try:
                    index = int(row["elem_idx"])
                    charge = float(row["charge_C"])
                except (TypeError, ValueError) as exc:
                    raise ValidationError(
                        "charges.csv contains malformed values"
                    ) from exc
                if index in seen or not math.isfinite(charge):
                    raise ValidationError(
                        "charges.csv contains duplicate/nonfinite values"
                    )
                seen.add(index)
    except OSError as exc:
        raise ValidationError(f"failed to read {path}: {exc}") from exc
    if seen != set(range(1, expected_rows + 1)):
        raise ValidationError("charges.csv element index coverage mismatch")


def _validate_mesh_potential(path: Path, expected_rows: int) -> None:
    try:
        with path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None or not {
                "elem_idx",
                "potential_V",
            }.issubset(reader.fieldnames):
                raise ValidationError("mesh_potential.csv is missing required columns")
            seen: set[int] = set()
            for row in reader:
                try:
                    index = int(row["elem_idx"])
                    potential = float(row["potential_V"])
                except (TypeError, ValueError) as exc:
                    raise ValidationError(
                        "mesh_potential.csv contains malformed values"
                    ) from exc
                if index in seen or not math.isfinite(potential):
                    raise ValidationError(
                        "mesh_potential.csv contains duplicate/nonfinite values"
                    )
                seen.add(index)
    except OSError as exc:
        raise ValidationError(f"failed to read mesh_potential.csv: {exc}") from exc
    if seen != set(range(1, expected_rows + 1)):
        raise ValidationError("mesh_potential.csv element index coverage mismatch")


def _validate_mesh_triangles(path: Path, expected_rows: int) -> Counter[int]:
    try:
        with path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None or "elem_idx" not in reader.fieldnames:
                raise ValidationError("mesh_triangles.csv is missing elem_idx")
            numeric_columns = [
                name
                for name in reader.fieldnames
                if name not in {"elem_idx", "mesh_id"}
            ]
            if len(numeric_columns) < 9:
                raise ValidationError(
                    "mesh_triangles.csv is missing triangle coordinates"
                )
            seen: set[int] = set()
            mesh_counts: Counter[int] = Counter()
            for row in reader:
                try:
                    index = int(row["elem_idx"])
                    values = [float(row[name]) for name in numeric_columns]
                    mesh_id = int(row["mesh_id"]) if "mesh_id" in row else 1
                except (TypeError, ValueError) as exc:
                    raise ValidationError(
                        "mesh_triangles.csv contains malformed values"
                    ) from exc
                if (
                    index in seen
                    or mesh_id < 1
                    or any(not math.isfinite(value) for value in values)
                ):
                    raise ValidationError(
                        "mesh_triangles.csv contains duplicate/nonfinite values"
                    )
                seen.add(index)
                mesh_counts[mesh_id] += 1
    except OSError as exc:
        raise ValidationError(f"failed to read {path}: {exc}") from exc
    if seen != set(range(1, expected_rows + 1)):
        raise ValidationError("mesh_triangles.csv element index coverage mismatch")
    return mesh_counts


def _mesh_source_contract(
    path: Path,
    *,
    expected_meshes: int | None = None,
    expected_elements: int | None = None,
) -> dict[str, Any]:
    required = {
        "mesh_id",
        "source_kind",
        "template_kind",
        "surface_model",
        "epsilon_r",
        "elem_count",
    }
    try:
        with path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValidationError("mesh_sources.csv is missing required columns")
            ordered_mesh_ids: list[int] = []
            by_mesh_id: dict[str, dict[str, Any]] = {}
            element_total = 0
            for row in reader:
                try:
                    mesh_id = int(row["mesh_id"])
                    element_count = int(row["elem_count"])
                    epsilon_r = float(row["epsilon_r"])
                except (TypeError, ValueError) as exc:
                    raise ValidationError(
                        "mesh_sources.csv contains malformed values"
                    ) from exc
                if (
                    mesh_id < 1
                    or str(mesh_id) in by_mesh_id
                    or element_count < 1
                    or not math.isfinite(epsilon_r)
                    or epsilon_r < 1.0
                    or not row["source_kind"].strip()
                    or not row["template_kind"].strip()
                    or not row["surface_model"].strip()
                ):
                    raise ValidationError("mesh_sources.csv violates the mesh contract")
                ordered_mesh_ids.append(mesh_id)
                by_mesh_id[str(mesh_id)] = {
                    "source_kind": row["source_kind"].strip().lower(),
                    "template_kind": row["template_kind"].strip().lower(),
                    "surface_model": row["surface_model"].strip().lower(),
                    "epsilon_r": epsilon_r,
                    "elem_count": element_count,
                }
                element_total += element_count
    except OSError as exc:
        raise ValidationError(f"failed to read mesh_sources.csv: {exc}") from exc
    if expected_meshes is not None and len(ordered_mesh_ids) != expected_meshes:
        raise ValidationError("mesh_sources.csv mesh count mismatch")
    if expected_elements is not None and element_total != expected_elements:
        raise ValidationError("mesh_sources.csv element count mismatch")
    return {
        "path": str(path.resolve()),
        "sha256": _sha256(path),
        "ordered_mesh_ids": ordered_mesh_ids,
        "by_mesh_id": by_mesh_id,
    }


def _validate_mesh_sources(
    path: Path,
    *,
    expected_meshes: int,
    expected_elements: int,
    archived_contract: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    contract = _mesh_source_contract(
        path,
        expected_meshes=expected_meshes,
        expected_elements=expected_elements,
    )
    if archived_contract is not None:
        actual_semantics = {
            "ordered_mesh_ids": contract["ordered_mesh_ids"],
            "by_mesh_id": contract["by_mesh_id"],
        }
        expected_semantics = {
            "ordered_mesh_ids": archived_contract.get("ordered_mesh_ids"),
            "by_mesh_id": archived_contract.get("by_mesh_id"),
        }
        if actual_semantics != expected_semantics:
            raise ValidationError(
                "mesh_sources.csv differs from the archived mesh source contract"
            )
    return contract


def _validate_charge_ledger(
    path: Path,
    *,
    expected_batches: int,
    expected_species: int,
    processed: int,
    absorbed: int,
    escaped_boundary: int,
    survived_max_step: int,
    soft_discarded: int,
) -> None:
    count_fields = (
        "injected_count",
        "emitted_count",
        "absorbed_count",
        "escaped_count",
        "discarded_unresolved_count",
    )
    required = {
        "batch",
        "species_idx",
        "injected_from_remote_C",
        "emitted_from_surface_C",
        "absorbed_on_surface_C",
        "escaped_to_infinity_C",
        "discarded_unresolved_C",
        *count_fields,
    }
    seen: set[int] = set()
    totals = {field: 0 for field in count_fields}
    try:
        with path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValidationError("charge_ledger.csv is missing required columns")
            for row in reader:
                try:
                    batch = int(row["batch"])
                    species = int(row["species_idx"])
                    charges = [
                        float(row[name]) for name in required if name.endswith("_C")
                    ]
                    counts = {field: int(row[field]) for field in count_fields}
                except (TypeError, ValueError) as exc:
                    raise ValidationError(
                        "charge_ledger.csv contains malformed values"
                    ) from exc
                if (
                    batch != expected_batches
                    or species in seen
                    or species < 1
                    or any(not math.isfinite(value) for value in charges)
                    or any(value < 0 for value in counts.values())
                ):
                    raise ValidationError(
                        "charge_ledger.csv violates the ledger contract"
                    )
                seen.add(species)
                for field, value in counts.items():
                    totals[field] += value
    except OSError as exc:
        raise ValidationError(f"failed to read {path}: {exc}") from exc
    if seen != set(range(1, expected_species + 1)):
        raise ValidationError("charge_ledger.csv species coverage mismatch")
    if totals["injected_count"] + totals["emitted_count"] != processed:
        raise ValidationError("charge_ledger processed particle count mismatch")
    if totals["absorbed_count"] != absorbed:
        raise ValidationError("charge_ledger absorbed particle count mismatch")
    if totals["escaped_count"] != escaped_boundary:
        raise ValidationError("charge_ledger outcome counts mismatch")
    if totals["discarded_unresolved_count"] != survived_max_step + soft_discarded:
        raise ValidationError("charge_ledger unresolved particle count mismatch")


def _validate_checkpoint(output: Path, world_size: int) -> None:
    required = [output / "macro_residuals.csv"]
    required.extend(
        output / f"rng_state_rank{rank:05d}.txt" for rank in range(world_size)
    )
    missing = [
        str(path) for path in required if not path.is_file() or path.stat().st_size == 0
    ]
    if missing:
        raise ValidationError(f"checkpoint files are missing/incomplete: {missing}")


def _validate_cached_electrostatic_state(summary: Mapping[str, str]) -> None:
    if not _summary_bool(summary, "electrostatic_split_periodic_active"):
        raise ValidationError("cached electrostatic state is not active")
    if summary.get("electrostatic_status") != "cached_kneq0":
        raise ValidationError("cached electrostatic status mismatch")
    try:
        _summary_float(summary, "gauss_residual_C")
    except ValidationError as exc:
        raise ValidationError(f"incomplete cached electrostatic state: {exc}") from exc


def _validate_history(
    path: Path,
    *,
    mesh_nelem: int,
    expected_batches: int,
    stride: int,
    restart_batch: int,
) -> None:
    first = restart_batch + 1
    expected_sequence = [
        batch
        for batch in range(first, expected_batches + 1)
        if (batch - 1) % stride == 0
    ]
    actual_sequence: list[int] = []
    current_batch: int | None = None
    seen_elements: set[int] = set()

    def finish_batch() -> None:
        if current_batch is None:
            return
        if seen_elements != set(range(1, mesh_nelem + 1)):
            raise ValidationError(
                f"charge_history.csv element coverage mismatch at batch {current_batch}"
            )
        actual_sequence.append(current_batch)

    try:
        with path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            required = {
                "batch",
                "processed_particles",
                "rel_change",
                "elem_idx",
                "charge_C",
            }
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValidationError("charge_history.csv is missing required columns")
            for row in reader:
                try:
                    batch = int(row["batch"])
                    processed = int(row["processed_particles"])
                    rel_change = float(row["rel_change"])
                    element = int(row["elem_idx"])
                    charge = float(row["charge_C"])
                except (TypeError, ValueError) as exc:
                    raise ValidationError(
                        "charge_history.csv contains malformed values"
                    ) from exc
                if (
                    processed < 0
                    or rel_change < 0.0
                    or not all(math.isfinite(value) for value in (rel_change, charge))
                ):
                    raise ValidationError("charge_history.csv contains invalid values")
                if current_batch is None:
                    current_batch = batch
                elif batch != current_batch:
                    finish_batch()
                    current_batch = batch
                    seen_elements = set()
                if element in seen_elements:
                    raise ValidationError(
                        f"charge_history.csv duplicate element at batch {batch}"
                    )
                seen_elements.add(element)
        finish_batch()
    except OSError as exc:
        raise ValidationError(f"failed to read {path}: {exc}") from exc
    if actual_sequence != expected_sequence:
        raise ValidationError(
            f"history batch sequence mismatch: {actual_sequence} != {expected_sequence}"
        )


def _matching_case(
    manifest: Mapping[str, Any], case_dir: Path
) -> tuple[str, dict[str, Any]]:
    for name, value in manifest.get("cases", {}).items():
        if Path(str(value["output_dir"])) == case_dir:
            return str(name), dict(value)
    raise ValidationError(f"case directory is not declared in manifest: {case_dir}")


def _canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")


def _json_sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_json_bytes(value)).hexdigest()


def _receipt_payload_sha256(receipt: Mapping[str, Any]) -> str:
    payload = dict(receipt)
    payload.pop("receipt_payload_sha256", None)
    return _json_sha256(payload)


def _output_inventory(output: Path) -> dict[str, str]:
    if output.is_symlink():
        raise ValidationError(f"output inventory rejects symlink root: {output}")
    files: dict[str, str] = {}
    try:
        paths = sorted(output.rglob("*"))
    except OSError as exc:
        raise ValidationError(f"failed to inventory output {output}: {exc}") from exc
    for path in paths:
        if path.is_symlink():
            raise ValidationError(f"output inventory rejects symlink: {path}")
        if path.is_file():
            files[path.relative_to(output).as_posix()] = _sha256(path)
    if not files:
        raise ValidationError(f"output inventory is empty: {output}")
    return files


def _load_execution_receipt(path: Path) -> dict[str, Any]:
    if path.is_symlink():
        raise ValidationError(f"execution receipt must not be a symlink: {path}")
    try:
        raw = path.read_bytes()
        value = json.loads(raw.decode("utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValidationError(f"invalid execution receipt {path}: {exc}") from exc
    except UnicodeDecodeError as exc:
        raise ValidationError(
            f"invalid execution receipt encoding {path}: {exc}"
        ) from exc
    if not isinstance(value, dict) or value.get("receipt_schema_version") != 1:
        raise ValidationError(f"invalid execution receipt schema: {path}")
    digest = value.get("receipt_payload_sha256")
    if not isinstance(digest, str) or digest != _receipt_payload_sha256(value):
        raise ValidationError(f"execution receipt self-hash mismatch: {path}")
    canonical = (
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    ).encode("utf-8")
    if raw != canonical:
        raise ValidationError(f"execution receipt is not canonical JSON: {path}")
    return value


def _assert_dependency_receipt_current(
    path: Path,
    *,
    validation_root: Path,
    expected_case: str,
    expected_case_dir: Path,
    current_manifest_sha256: str,
    label: str,
) -> dict[str, Any]:
    path = _require_expected_path(
        validation_root,
        path,
        f"provenance/verified/{expected_case}.json",
        label=f"{label} receipt",
    )
    if not path.is_file():
        raise ValidationError(f"{label} receipt is missing: {path}")
    receipt = _load_execution_receipt(path)
    if receipt.get("case") != expected_case:
        raise ValidationError(f"{label} receipt case mismatch")
    if receipt.get("manifest_sha256") != current_manifest_sha256:
        raise ValidationError(f"{label} receipt manifest mismatch")
    case_dir = _require_descendant_path(
        validation_root,
        str(receipt.get("case_dir", "")),
        prefix="run",
        label=f"{label} receipt case directory",
    )
    if case_dir != expected_case_dir:
        raise ValidationError(f"{label} receipt case directory mismatch")
    if _output_inventory(case_dir) != receipt.get("outputs"):
        raise ValidationError(f"{label} outputs differ from immutable receipt")
    return receipt


def _receipt_comparable(receipt: Mapping[str, Any]) -> dict[str, Any]:
    value = dict(receipt)
    value.pop("verified_at", None)
    value.pop("receipt_payload_sha256", None)
    return value


def _load_submission_journal(
    root: Path,
    *,
    require_complete: bool,
) -> tuple[Path, bool, dict[str, str]]:
    path = _require_expected_path(
        root,
        root / "submit/job_ids.json",
        "submit/job_ids.json",
        label="job submission journal",
    )
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValidationError(
            f"job submission journal is missing or invalid: {exc}"
        ) from exc
    expected_keys = {"submission_complete", *EXPECTED_SUBMITTED_JOBS}
    if not isinstance(value, Mapping) or set(value) != expected_keys:
        raise ValidationError("job submission journal has an unexpected key set")
    complete = value.get("submission_complete")
    if type(complete) is not bool:
        raise ValidationError("job submission journal completion flag must be boolean")
    identifiers: dict[str, str] = {}
    for name in EXPECTED_SUBMITTED_JOBS:
        identifier = value.get(name)
        if not isinstance(identifier, str) or (
            identifier != "" and not identifier.isdigit()
        ):
            raise ValidationError("job submission journal has invalid job IDs")
        identifiers[name] = identifier
    assigned = [value for value in identifiers.values() if value]
    if len(set(assigned)) != len(assigned):
        raise ValidationError("job submission journal has duplicate job IDs")
    if complete and len(assigned) != len(EXPECTED_SUBMITTED_JOBS):
        raise ValidationError("completed job submission journal has missing job IDs")
    if require_complete and not complete:
        raise ValidationError("job submission journal is incomplete")
    if require_complete and len(assigned) != len(EXPECTED_SUBMITTED_JOBS):
        raise ValidationError("job submission journal is incomplete")
    return path, complete, identifiers


def _execution_producer_binding(
    *,
    root: Path,
    manifest: Mapping[str, Any],
    case_name: str,
    case: Mapping[str, Any],
    config_path: Path,
    producer_job_role: str | None,
    existing_receipt: Mapping[str, Any] | None,
) -> dict[str, str] | None:
    if manifest.get("build_origin") is None:
        return None
    expected_role = CASE_PRODUCER_JOB_ROLES.get(case_name)
    if expected_role is None:
        raise ValidationError(
            f"case {case_name} is not executable in the production job graph"
        )
    if producer_job_role is not None and producer_job_role != expected_role:
        raise ValidationError(
            f"producer job role mismatch for {case_name}: "
            f"{producer_job_role!r} != {expected_role!r}"
        )
    if existing_receipt is None and producer_job_role is None:
        raise ValidationError(
            f"clean production receipt requires a producer job role for {case_name}"
        )

    _journal_path, _complete, identifiers = _load_submission_journal(
        root,
        require_complete=existing_receipt is not None,
    )
    job_id = identifiers[expected_role]
    if not job_id:
        raise ValidationError(f"producer job ID is missing for role {expected_role}")
    expected = {
        "job_role": expected_role,
        "job_id": job_id,
        "config_sha256": str(case.get("config_sha256", "")),
    }
    if existing_receipt is not None:
        recorded = existing_receipt.get("execution_producer")
        if not isinstance(recorded, Mapping) or dict(recorded) != expected:
            raise ValidationError(
                f"execution producer binding differs from immutable receipt for {case_name}"
            )
    if producer_job_role is not None:
        current_job = os.environ.get("SLURM_JOB_ID")
        if current_job is None or current_job != job_id:
            raise ValidationError(
                f"SLURM_JOB_ID differs from producer job {job_id} for {case_name}"
            )

    job_name = EXPECTED_SUBMITTED_JOBS[expected_role][0]
    token = f"{job_id}.{job_name}"
    hash_path = _require_expected_path(
        root,
        root / "provenance/hashes" / f"{token}.sha256",
        f"provenance/hashes/{token}.sha256",
        label=f"producer {expected_role} config hash log",
    )
    if not hash_path.is_file():
        raise ValidationError(
            f"producer config hash log is missing for {case_name}: {hash_path}"
        )
    try:
        hash_lines = set(hash_path.read_text(encoding="utf-8").splitlines())
    except OSError as exc:
        raise ValidationError(
            f"failed to read producer config hash log for {case_name}: {exc}"
        ) from exc
    if f"{config_path}: OK" not in hash_lines:
        raise ValidationError(
            f"producer config hash log does not contain the exact case path for {case_name}"
        )
    return expected


def _publish_execution_receipt(path: Path, report: Mapping[str, Any]) -> dict[str, Any]:
    value = dict(report)
    value["receipt_payload_sha256"] = _receipt_payload_sha256(value)
    serialized = (
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    ).encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(serialized)
            stream.flush()
            os.fsync(stream.fileno())
        temporary.chmod(0o444)
        try:
            os.link(temporary, path)
        except FileExistsError:
            existing = _load_execution_receipt(path)
            if _receipt_comparable(existing) != _receipt_comparable(value):
                raise ValidationError(
                    f"current output differs from immutable execution receipt: {path}"
                )
            return existing
        directory_descriptor = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    finally:
        temporary.unlink(missing_ok=True)
    return value


def verify_run(
    case_dir: str | Path,
    expected_batches: int,
    *,
    require_existing_receipt: bool = False,
    producer_job_role: str | None = None,
) -> dict[str, Any]:
    """Verify one completed output against its staged case and diagnostics."""

    requested_output = _lexical_absolute(case_dir)
    root = _find_validation_root(requested_output)
    manifest = _load_manifest(root)
    _require_expected_path(
        root,
        str(manifest.get("validation_root", "")),
        None,
        label="validation root",
    )
    cases = manifest.get("cases")
    if not isinstance(cases, Mapping):
        raise ValidationError("manifest cases must be an object")
    _validate_case_graph(root, cases)
    name, case = _matching_case(manifest, requested_output)
    output = _require_expected_path(
        root,
        requested_output,
        str(EXPECTED_CASE_GRAPH[name]["output_relative"]),
        label=f"case {name} output",
    )
    receipt_path = _require_expected_path(
        root,
        root / "provenance/verified" / f"{name}.json",
        f"provenance/verified/{name}.json",
        label=f"case {name} execution receipt",
    )
    existing_receipt = (
        _load_execution_receipt(receipt_path) if receipt_path.is_file() else None
    )
    if require_existing_receipt and existing_receipt is None:
        raise ValidationError(
            f"strict verification requires existing execution receipt: {name}"
        )
    manifest_sha256 = _sha256(root / "manifest.json")
    expected = int(expected_batches)
    if expected != int(case["batch_count"]):
        raise ValidationError(
            f"expected batch argument {expected} differs from staged case {case['batch_count']}"
        )
    summary = _summary(output / "summary.txt")
    expected_build_origin = manifest.get("build_origin")
    if expected_build_origin is not None:
        if _summary_build_origin(summary) != expected_build_origin:
            raise ValidationError(
                "summary build origin differs from the staged artifact build origin"
            )
    if _summary_int(summary, "checkpoint_schema_version") != 3:
        raise ValidationError("checkpoint_schema_version must be 3")
    mesh_nelem = _summary_int(summary, "mesh_nelem")
    mesh_count = _summary_int(summary, "mesh_count")
    world_size = _summary_int(summary, "mpi_world_size")
    processed = _summary_int(summary, "processed_particles")
    absorbed = _summary_int(summary, "absorbed")
    escaped = _summary_int(summary, "escaped")
    batches = _summary_int(summary, "batches")
    escaped_boundary = _summary_int(summary, "escaped_boundary")
    survived = _summary_int(summary, "survived_max_step")
    soft_discarded = (
        _summary_int(summary, "multiple_box_events_soft_discarded")
        if "multiple_box_events_soft_discarded" in summary
        else 0
    )
    _summary_float(summary, "last_rel_change")
    if mesh_nelem < 1 or mesh_count < 1:
        raise ValidationError("mesh counts must be positive")
    expected_output = manifest.get("archive", {}).get("expected_output", {})
    if expected_output.get("mesh_nelem") is not None and mesh_nelem != int(
        expected_output["mesh_nelem"]
    ):
        raise ValidationError("mesh_nelem differs from the archived geometry")
    if expected_output.get("mesh_count") is not None and mesh_count != int(
        expected_output["mesh_count"]
    ):
        raise ValidationError("mesh_count differs from the archived geometry")
    if world_size != EXPECTED_RESOURCES["mpi_processes"]:
        raise ValidationError(f"mpi_world_size mismatch: {world_size}")
    if batches != expected:
        raise ValidationError(f"batches mismatch: {batches} != {expected}")
    if processed != absorbed + escaped:
        raise ValidationError("processed_particles != absorbed + escaped")
    if escaped != escaped_boundary + survived + soft_discarded:
        raise ValidationError(
            "escaped != escaped_boundary + survived_max_step + multiple_box_events_soft_discarded"
        )
    for key in ("model_fingerprint", "mesh_fingerprint", "species_fingerprint"):
        value = summary.get(key, "")
        if len(value) != 16:
            raise ValidationError(f"summary {key} must be a 16-character fingerprint")
    for key, expected_value in PRODUCTION_FIELD_EXECUTION_CONTRACT.items():
        actual_value = summary.get(key)
        if actual_value != expected_value:
            raise ValidationError(
                f"summary {key}={actual_value!r}; expected {expected_value!r}"
            )
    residual = _summary_float(summary, "charge_ledger_residual_C")
    if abs(residual) > 1.0e-12:
        raise ValidationError(f"charge ledger residual is too large: {residual}")
    if any("fallback" in value.lower() for value in summary.values()):
        raise ValidationError("summary reports a forbidden fallback")
    config_path = Path(str(case["config_path"]))
    case_config = _load_toml(config_path)
    output_config = case_config.get("output", {})
    if not isinstance(output_config, dict):
        raise ValidationError(f"output configuration is invalid for {name}")
    _validate_charges(output / "charges.csv", mesh_nelem)
    write_mesh_potential = bool(output_config.get("write_mesh_potential", False))
    if write_mesh_potential:
        _validate_mesh_potential(output / "mesh_potential.csv", mesh_nelem)
    triangle_mesh_counts = _validate_mesh_triangles(
        output / "mesh_triangles.csv", mesh_nelem
    )
    source_contract = _validate_mesh_sources(
        output / "mesh_sources.csv",
        expected_meshes=mesh_count,
        expected_elements=mesh_nelem,
        archived_contract=manifest.get("archive", {}).get("mesh_source_contract"),
    )
    source_mesh_counts = Counter(
        {
            int(mesh_id): int(value["elem_count"])
            for mesh_id, value in source_contract["by_mesh_id"].items()
        }
    )
    if triangle_mesh_counts != source_mesh_counts:
        raise ValidationError(
            "mesh_sources.csv and mesh_triangles.csv per-mesh element counts differ"
        )
    ledger_species = _summary_int(summary, "charge_ledger_nspecies")
    if _summary_int(summary, "charge_ledger_batch_count") != batches:
        raise ValidationError("charge ledger batch count differs from summary batches")
    _validate_charge_ledger(
        output / "charge_ledger.csv",
        expected_batches=batches,
        expected_species=ledger_species,
        processed=processed,
        absorbed=absorbed,
        escaped_boundary=escaped_boundary,
        survived_max_step=survived,
        soft_discarded=soft_discarded,
    )
    _validate_checkpoint(output, world_size)
    restart_batch = 0
    if case.get("restart_from") is not None:
        restart_batch = _summary_int(
            _summary(Path(str(case["restart_from"])) / "summary.txt"), "batches"
        )
    _validate_history(
        output / "charge_history.csv",
        mesh_nelem=mesh_nelem,
        expected_batches=batches,
        stride=int(case["history_stride"]),
        restart_batch=restart_batch,
    )

    if _sha256(config_path) != case.get("config_sha256"):
        raise ValidationError(f"input hash mismatch for {name}")
    execution_producer = _execution_producer_binding(
        root=root,
        manifest=manifest,
        case_name=name,
        case=case,
        config_path=config_path,
        producer_job_role=producer_job_role,
        existing_receipt=existing_receipt,
    )
    binary = Path(str(manifest["binary"]["staged_path"]))
    if _sha256(binary) != manifest["binary"]["sha256"]:
        raise ValidationError("binary hash mismatch")
    if expected_build_origin is not None:
        analysis_library = manifest.get("analysis_library")
        if not isinstance(analysis_library, Mapping):
            raise ValidationError(
                "build origin verification requires an analysis library"
            )
        _verify_staged_build_origin(
            manifest,
            binary=binary,
            library=Path(str(analysis_library.get("staged_path", ""))),
        )

    expectation = str(case["cache_expectation"])
    dependencies: dict[str, Any] = {}
    cache: dict[str, Any] = {
        "expectation": expectation,
        "hit": None,
        "build_count": None,
        "fingerprint": None,
        "path": None,
        "path_sha256": None,
    }
    backend = summary.get("periodic2_nonzero_mode_backend", "")
    if case["periodic_model"] == "finite_configured":
        if backend != "legacy_finite_images":
            raise ValidationError(f"finite case backend mismatch: {backend!r}")
        if summary.get("periodic2_zero_mode_policy") != "legacy_not_decomposed":
            raise ValidationError("finite case zero-mode policy mismatch")
        if summary.get("periodic2_lower_boundary_model") != "legacy_implicit":
            raise ValidationError("finite case lower-boundary model mismatch")
        if "periodic2_cache_hit" in summary and _summary_bool(
            summary, "periodic2_cache_hit"
        ):
            raise ValidationError("finite case unexpectedly reports a cache hit")
    else:
        if backend != "cached_kneq0":
            raise ValidationError(f"infinite case backend mismatch: {backend!r}")
        if summary.get("periodic2_zero_mode_policy") != "exclude_k0":
            raise ValidationError("infinite case zero-mode policy mismatch")
        if summary.get("periodic2_lower_boundary_model") != "e_bottom_zero":
            raise ValidationError("infinite case lower-boundary model mismatch")
        _validate_cached_electrostatic_state(summary)
        tolerance = _summary_float(summary, "periodic2_generation_tolerance")
        if tolerance != DEFAULT_GENERATION_TOLERANCE:
            raise ValidationError("periodic generation tolerance mismatch")
        hit = _summary_bool(summary, "periodic2_cache_hit")
        build_count = _summary_int(summary, "periodic2_operator_build_count")
        fingerprint = summary.get("periodic2_cache_fingerprint", "")
        if not fingerprint:
            raise ValidationError("cached case has an empty cache fingerprint")
        cache_path = _require_descendant_path(
            root,
            summary.get("periodic2_cache_path", ""),
            prefix="cache/periodic2",
            label="cached operator path",
        )
        if not cache_path.is_file() or cache_path.stat().st_size == 0:
            raise ValidationError(f"cached operator is missing or empty: {cache_path}")
        cache_hash = _sha256(cache_path)
        cache.update(
            hit=hit,
            build_count=build_count,
            fingerprint=fingerprint,
            path=str(cache_path),
            path_sha256=cache_hash,
        )
        if expectation == "miss" and (hit or build_count != 1):
            raise ValidationError(
                "cache prime must report a cache miss and build_count=1"
            )
        if expectation == "hit" and (not hit or build_count != 0):
            raise ValidationError(
                "warm cached case must report a cache hit and build_count=0"
            )
        prime_dir = Path(str(manifest["cases"]["cache_prime"]["output_dir"]))
        if expectation == "hit":
            prime_summary_path = prime_dir / "summary.txt"
            if not prime_summary_path.is_file():
                raise ValidationError(
                    "warm cached case requires the cache prime summary"
                )
            prime_summary = _summary(prime_summary_path)
            prime_fingerprint = prime_summary.get("periodic2_cache_fingerprint", "")
            if fingerprint != prime_fingerprint:
                raise ValidationError("warm cache fingerprint differs from cache prime")
            prime_receipt_path = root / "provenance/verified/cache_prime.json"
            prime_receipt = _assert_dependency_receipt_current(
                prime_receipt_path,
                validation_root=root,
                expected_case="cache_prime",
                expected_case_dir=prime_dir,
                current_manifest_sha256=manifest_sha256,
                label="cache prime",
            )
            prime_cache = prime_receipt.get("cache", {})
            if not isinstance(prime_cache, dict):
                raise ValidationError("cache prime receipt has invalid cache metadata")
            if prime_cache.get("path_sha256") != cache_hash:
                raise ValidationError(
                    "warm cached operator hash differs from the verified cache prime"
                )
            dependencies["cache_prime"] = {
                "case": "cache_prime",
                "receipt_sha256": _sha256(prime_receipt_path),
            }

    restart_from = case.get("restart_from")
    if restart_from is not None:
        restart_output = Path(str(restart_from))
        parent_name, _parent_case = _matching_case(manifest, restart_output)
        parent_receipt_path = root / "provenance/verified" / f"{parent_name}.json"
        _assert_dependency_receipt_current(
            parent_receipt_path,
            validation_root=root,
            expected_case=parent_name,
            expected_case_dir=restart_output,
            current_manifest_sha256=manifest_sha256,
            label="restart parent",
        )
        dependencies["restart"] = {
            "case": parent_name,
            "receipt_sha256": _sha256(parent_receipt_path),
        }
        previous = _summary(restart_output / "summary.txt")
        for key in (
            "model_fingerprint",
            "mesh_fingerprint",
            "species_fingerprint",
            *PRODUCTION_FIELD_EXECUTION_CONTRACT,
            *BUILD_INFO_KEYS,
        ):
            if previous.get(key) != summary.get(key):
                raise ValidationError(f"restart {key} differs from prior segment")
    required_outputs = [
        output / "summary.txt",
        output / "charges.csv",
        output / "mesh_triangles.csv",
        output / "charge_ledger.csv",
        output / "macro_residuals.csv",
        output / "charge_history.csv",
        output / "mesh_sources.csv",
    ]
    if write_mesh_potential:
        required_outputs.append(output / "mesh_potential.csv")
    required_outputs.extend(
        output / f"rng_state_rank{rank:05d}.txt" for rank in range(world_size)
    )
    for path in required_outputs:
        if not path.is_file():
            raise ValidationError(f"required output is missing: {path}")
    output_hashes = _output_inventory(output)
    report = {
        "receipt_schema_version": 1,
        "status": "ok",
        "case": name,
        "case_dir": str(output),
        "batches": batches,
        "mesh_nelem": mesh_nelem,
        "mesh_count": mesh_count,
        "manifest_sha256": manifest_sha256,
        "case_spec_sha256": _json_sha256(case),
        "fingerprints": {
            key: summary[key]
            for key in (
                "model_fingerprint",
                "mesh_fingerprint",
                "species_fingerprint",
            )
        },
        "field_execution_contract": {
            key: summary[key] for key in PRODUCTION_FIELD_EXECUTION_CONTRACT
        },
        "build_origin": expected_build_origin,
        "config": {
            "path": str(config_path),
            "sha256": case.get("config_sha256"),
        },
        "binary": {
            "path": str(binary),
            "sha256": manifest["binary"]["sha256"],
        },
        "cache": cache,
        "dependencies": dependencies,
        "outputs": output_hashes,
        "verified_at": _utc_now(),
    }
    if execution_producer is not None:
        report["execution_producer"] = execution_producer
    if existing_receipt is not None:
        if _receipt_comparable(existing_receipt) != _receipt_comparable(report):
            raise ValidationError(
                f"current output differs from immutable execution receipt: {receipt_path}"
            )
        return existing_receipt
    return _publish_execution_receipt(receipt_path, report)


def _write_csv(
    path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def _case_result(
    case: str,
    output: Path,
    config_path: Path,
    source_version: str,
) -> tuple[dict[str, Any], Beach | None]:
    if not (output / "summary.txt").is_file() or not (output / "charges.csv").is_file():
        return {
            "case": case,
            "status": "missing",
            "source_version": source_version,
            "message": f"output not available: {output}",
        }, None
    try:
        run = Beach(output, config_path=config_path)
        result = run.result
        charges = [float(value) for value in result.charges]
        row = {
            "case": case,
            "status": "available",
            "source_version": source_version,
            "batches": result.batches,
            "mesh_nelem": result.mesh_nelem,
            "mesh_count": len(set(int(v) for v in result.mesh_ids))
            if result.mesh_ids is not None
            else "",
            "processed_particles": result.processed_particles,
            "absorbed": result.absorbed,
            "escaped": result.escaped,
            "escaped_boundary": result.escaped_boundary,
            "survived_max_step": result.survived_max_step,
            "last_rel_change": result.last_rel_change,
            "total_charge_C": math.fsum(charges),
            "model_fingerprint": result.model_fingerprint or "",
            "mesh_fingerprint": result.mesh_fingerprint or "",
            "species_fingerprint": result.species_fingerprint or "",
            "message": "",
        }
        return row, run
    except Exception as exc:
        return {
            "case": case,
            "status": "invalid",
            "source_version": source_version,
            "message": str(exc),
        }, None


def _history_rows(case: str, paths: Sequence[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    seen_batches: set[int] = set()
    for path in paths:
        if not path.is_file():
            continue
        current_batch: int | None = None
        processed = ""
        rel_change = ""
        charge_sum = 0.0
        with path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            required = {"batch", "processed_particles", "rel_change", "charge_C"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValidationError(f"history is missing required columns: {path}")
            for raw in reader:
                batch = int(raw["batch"])
                if current_batch is not None and batch != current_batch:
                    if current_batch in seen_batches:
                        raise ValidationError(
                            f"duplicate history batch {current_batch} for {case}"
                        )
                    rows.append(
                        {
                            "case": case,
                            "batch": current_batch,
                            "processed_particles": processed,
                            "rel_change": rel_change,
                            "total_charge_C": charge_sum,
                            "source_segment": str(path.parent),
                        }
                    )
                    seen_batches.add(current_batch)
                    charge_sum = 0.0
                current_batch = batch
                processed = raw["processed_particles"]
                rel_change = raw["rel_change"]
                charge = float(raw["charge_C"])
                if not math.isfinite(charge):
                    raise ValidationError(f"nonfinite history charge in {path}")
                charge_sum += charge
        if current_batch is not None:
            if current_batch in seen_batches:
                raise ValidationError(
                    f"duplicate history batch {current_batch} for {case}"
                )
            rows.append(
                {
                    "case": case,
                    "batch": current_batch,
                    "processed_particles": processed,
                    "rel_change": rel_change,
                    "total_charge_C": charge_sum,
                    "source_segment": str(path.parent),
                }
            )
            seen_batches.add(current_batch)
    rows.sort(key=lambda row: int(row["batch"]))
    return rows


def _comparison_rows(run_rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    by_case = {str(row["case"]): row for row in run_rows}
    definitions = [
        (
            "archived_v1_3_to_new_finite",
            "archive_version_drift",
            "archived_v1_3",
            "new_finite_configured",
            "version_and_reproducibility_drift",
        ),
        (
            "new_finite_to_new_infinite",
            "actual_end_to_end",
            "new_finite_configured",
            "new_infinite_physical",
            "boundary_model_effect",
        ),
    ]
    metrics = ("processed_particles", "absorbed", "escaped", "total_charge_C")
    rows: list[dict[str, Any]] = []
    for label, kind, left_name, right_name, interpretation in definitions:
        left = by_case[left_name]
        right = by_case[right_name]
        if left.get("status") != "available" or right.get("status") != "available":
            rows.append(
                {
                    "comparison": label,
                    "comparison_kind": kind,
                    "metric": "availability",
                    "left_value": left.get("status"),
                    "right_value": right.get("status"),
                    "status": "missing_comparison",
                    "interpretation": interpretation,
                }
            )
            continue
        for metric in metrics:
            left_value = float(left[metric])
            right_value = float(right[metric])
            absolute = right_value - left_value
            scale = max(abs(left_value), abs(right_value))
            rows.append(
                {
                    "comparison": label,
                    "comparison_kind": kind,
                    "sample_kind": "final",
                    "metric": metric,
                    "left_value": left_value,
                    "right_value": right_value,
                    "signed_difference": absolute,
                    "absolute_difference": abs(absolute),
                    "relative_difference": (
                        0.0 if scale == 0.0 else abs(absolute) / scale
                    ),
                    "status": "computed",
                    "interpretation": interpretation,
                }
            )
    return rows


def _charge_vector_comparison_rows(
    *,
    comparison: str,
    comparison_kind: str,
    sample_kind: str,
    resolved_batch: int,
    scope: str,
    mesh_id: int | str,
    left_snapshot_id: str,
    right_snapshot_id: str,
    left: np.ndarray,
    right: np.ndarray,
    interpretation: str,
) -> list[dict[str, Any]]:
    left_value = np.asarray(left, dtype=np.float64)
    right_value = np.asarray(right, dtype=np.float64)
    if (
        left_value.ndim != 1
        or right_value.shape != left_value.shape
        or left_value.size == 0
        or not np.all(np.isfinite(left_value))
        or not np.all(np.isfinite(right_value))
    ):
        raise ValidationError("charge comparison requires matching finite vectors")
    base = {
        "comparison": comparison,
        "comparison_kind": comparison_kind,
        "sample_kind": sample_kind,
        "resolved_batch": resolved_batch,
        "scope": scope,
        "mesh_id": mesh_id,
        "left_snapshot_id": left_snapshot_id,
        "right_snapshot_id": right_snapshot_id,
        "status": "computed",
        "interpretation": interpretation,
    }
    rows: list[dict[str, Any]] = []
    left_sum = float(math.fsum(float(value) for value in left_value))
    right_sum = float(math.fsum(float(value) for value in right_value))
    signed = right_sum - left_sum
    sum_scale = max(abs(left_sum), abs(right_sum))
    rows.append(
        {
            **base,
            "metric": "charge_sum_C",
            "left_value": left_sum,
            "right_value": right_sum,
            "signed_difference": signed,
            "absolute_difference": abs(signed),
            "relative_difference": 0.0 if sum_scale == 0.0 else abs(signed) / sum_scale,
        }
    )
    difference = right_value - left_value
    for metric, order in (
        ("element_charge_l1_C", 1),
        ("element_charge_l2_C", 2),
        ("element_charge_linf_C", np.inf),
    ):
        left_norm = float(np.linalg.norm(left_value, ord=order))
        right_norm = float(np.linalg.norm(right_value, ord=order))
        difference_norm = float(np.linalg.norm(difference, ord=order))
        scale = max(left_norm, right_norm)
        rows.append(
            {
                **base,
                "metric": metric,
                "left_value": left_norm,
                "right_value": right_norm,
                "signed_difference": "",
                "absolute_difference": difference_norm,
                "relative_difference": (
                    0.0 if scale == 0.0 else difference_norm / scale
                ),
            }
        )
    return rows


def _charge_vector_sha256(values: np.ndarray) -> str:
    canonical = np.ascontiguousarray(np.asarray(values, dtype="<f8"))
    return hashlib.sha256(canonical.tobytes()).hexdigest()


def _collect_charge_snapshots(
    *,
    case: str,
    run: Beach | None,
    history_paths: Sequence[Path],
) -> tuple[list[dict[str, Any]], dict[tuple[str, int], dict[str, Any]]]:
    if run is None:
        return (
            [
                {
                    "snapshot_id": "",
                    "case": case,
                    "status": "missing",
                    "message": "run output is unavailable",
                }
            ],
            {},
        )
    result = run.result
    mesh_ids = getattr(result, "mesh_ids", None)
    if mesh_ids is None:
        raise ValidationError(f"snapshot comparison requires mesh_ids for {case}")
    mesh_ids_array = np.asarray(mesh_ids, dtype=np.int64)
    if mesh_ids_array.shape != (int(result.mesh_nelem),):
        raise ValidationError(f"snapshot mesh_ids shape mismatch for {case}")
    history_snapshots: dict[int, dict[str, Any]] = {}
    source_hashes: dict[Path, str] = {}
    for path in history_paths:
        if not path.is_file():
            continue
        history = FortranChargeHistory(path, mesh_nelem=int(result.mesh_nelem))
        batches = np.asarray(history.batch_indices, dtype=np.int64)
        processed = np.asarray(history.processed_particles_by_batch, dtype=np.int64)
        rel_change = np.asarray(history.rel_change_by_batch, dtype=np.float64)
        if not (batches.shape == processed.shape == rel_change.shape):
            raise ValidationError(f"history metadata shape mismatch: {path}")
        source_hashes[path] = _sha256(path)
        for index, batch_value in enumerate(batches):
            batch = int(batch_value)
            if batch in history_snapshots:
                raise ValidationError(
                    f"duplicate history batch {batch} across segments for {case}"
                )
            history_snapshots[batch] = {
                "charges": np.asarray(history.get_step(batch), dtype=np.float64),
                "processed_particles": int(processed[index]),
                "rel_change": float(rel_change[index]),
                "source_file": path,
            }
    latest_history = max(history_snapshots, default=None)
    snapshots: dict[tuple[str, int], dict[str, Any]] = {}
    for batch, value in sorted(history_snapshots.items()):
        snapshots[("history", batch)] = {
            **value,
            "sample_kind": "history",
            "resolved_batch": batch,
            "is_latest_history": batch == latest_history,
        }
    final_batch = int(result.batches)
    final_charges = np.asarray(result.charges, dtype=np.float64)
    snapshots[("final", final_batch)] = {
        "charges": final_charges,
        "processed_particles": int(result.processed_particles),
        "rel_change": float(result.last_rel_change),
        "source_file": Path(getattr(result, "directory", ".")) / "charges.csv",
        "sample_kind": "final",
        "resolved_batch": final_batch,
        "is_latest_history": False,
    }
    final_source = Path(getattr(result, "directory", ".")) / "charges.csv"
    if final_source.is_file():
        source_hashes[final_source] = _sha256(final_source)

    rows: list[dict[str, Any]] = []
    unique_meshes = sorted(set(int(value) for value in mesh_ids_array))
    for (_kind, _batch), snapshot in snapshots.items():
        charges = np.asarray(snapshot["charges"], dtype=np.float64)
        if charges.shape != mesh_ids_array.shape or not np.all(np.isfinite(charges)):
            raise ValidationError(f"snapshot charges are invalid for {case}")
        snapshot_id = f"{case}:{snapshot['sample_kind']}:{snapshot['resolved_batch']}"
        snapshot["snapshot_id"] = snapshot_id
        source_file = Path(snapshot["source_file"])
        for scope, mesh_id, mask in (
            ("run", "", np.ones(charges.size, dtype=bool)),
            *(
                (
                    "object",
                    mesh_id,
                    mesh_ids_array == mesh_id,
                )
                for mesh_id in unique_meshes
            ),
        ):
            selected = charges[mask]
            rows.append(
                {
                    "snapshot_id": snapshot_id,
                    "case": case,
                    "sample_kind": snapshot["sample_kind"],
                    "resolved_batch": snapshot["resolved_batch"],
                    "run_final_batch": final_batch,
                    "is_latest_history": snapshot["is_latest_history"],
                    "mesh_fingerprint": getattr(result, "mesh_fingerprint", None) or "",
                    "field_source_model": getattr(
                        result, "field_source_model", "unknown"
                    ),
                    "field_kernel_id": getattr(result, "field_kernel_id", None) or "",
                    "scope": scope,
                    "mesh_id": mesh_id,
                    "processed_particles": snapshot["processed_particles"],
                    "rel_change": snapshot["rel_change"],
                    "charge_vector_sha256": _charge_vector_sha256(selected),
                    "source_file": str(source_file),
                    "source_file_sha256": source_hashes.get(source_file, ""),
                    "charge_sum_C": math.fsum(float(value) for value in selected),
                    "element_charge_l1_C": float(np.linalg.norm(selected, ord=1)),
                    "element_charge_l2_C": float(np.linalg.norm(selected, ord=2)),
                    "element_charge_linf_C": float(
                        np.linalg.norm(selected, ord=np.inf)
                    ),
                    "status": "available",
                    "message": "",
                }
            )
    return rows, snapshots


def _geometry_comparison(left_result: Any, right_result: Any) -> dict[str, Any]:
    left_mesh_ids_value = getattr(left_result, "mesh_ids", None)
    right_mesh_ids_value = getattr(right_result, "mesh_ids", None)
    left_triangles_value = getattr(left_result, "triangles", None)
    right_triangles_value = getattr(right_result, "triangles", None)
    if any(
        value is None
        for value in (
            left_mesh_ids_value,
            right_mesh_ids_value,
            left_triangles_value,
            right_triangles_value,
        )
    ):
        return {
            "status": "missing_geometry",
            "mesh_id_order_matches": False,
            "max_coordinate_difference_m": None,
            "coordinate_tolerance_m": None,
        }
    left_mesh_ids = np.asarray(left_mesh_ids_value, dtype=np.int64)
    right_mesh_ids = np.asarray(right_mesh_ids_value, dtype=np.int64)
    left_triangles = np.asarray(left_triangles_value, dtype=np.float64)
    right_triangles = np.asarray(right_triangles_value, dtype=np.float64)
    if (
        left_mesh_ids.ndim != 1
        or right_mesh_ids.ndim != 1
        or left_triangles.shape != (left_mesh_ids.size, 3, 3)
        or right_triangles.shape != (right_mesh_ids.size, 3, 3)
        or not np.all(np.isfinite(left_triangles))
        or not np.all(np.isfinite(right_triangles))
    ):
        return {
            "status": "invalid_geometry",
            "mesh_id_order_matches": False,
            "max_coordinate_difference_m": None,
            "coordinate_tolerance_m": None,
        }
    if left_mesh_ids.shape != right_mesh_ids.shape or not np.array_equal(
        left_mesh_ids, right_mesh_ids
    ):
        return {
            "status": "mesh_id_order_mismatch",
            "mesh_id_order_matches": False,
            "max_coordinate_difference_m": None,
            "coordinate_tolerance_m": None,
        }
    if left_triangles.shape != right_triangles.shape:
        return {
            "status": "triangle_shape_mismatch",
            "mesh_id_order_matches": True,
            "max_coordinate_difference_m": None,
            "coordinate_tolerance_m": None,
        }
    combined = np.concatenate(
        (left_triangles.reshape(-1, 3), right_triangles.reshape(-1, 3)),
        axis=0,
    )
    length_scale = float(np.max(np.ptp(combined, axis=0)))
    tolerance = max(
        1.0e-18,
        64.0 * np.finfo(np.float64).eps * length_scale,
    )
    max_difference = float(np.max(np.abs(left_triangles - right_triangles)))
    return {
        "status": "match" if max_difference <= tolerance else "coordinate_mismatch",
        "mesh_id_order_matches": True,
        "max_coordinate_difference_m": max_difference,
        "coordinate_tolerance_m": tolerance,
    }


def _analysis_geometry_validation(
    runs: Mapping[str, Beach | None],
) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for label, left_case, right_case in (
        (
            "archived_v1_3_to_new_finite",
            "archived_v1_3",
            "new_finite_configured",
        ),
        (
            "new_finite_to_new_infinite",
            "new_finite_configured",
            "new_infinite_physical",
        ),
    ):
        left = runs.get(left_case)
        right = runs.get(right_case)
        if left is None or right is None:
            comparison = {
                "status": "not_evaluated",
                "mesh_id_order_matches": False,
                "max_coordinate_difference_m": None,
                "coordinate_tolerance_m": None,
            }
        else:
            comparison = _geometry_comparison(left.result, right.result)
        result[label] = {
            "left_case": left_case,
            "right_case": right_case,
            **comparison,
        }
    return result


def _charge_snapshot_comparisons(
    snapshots: Mapping[str, Mapping[tuple[str, int], Mapping[str, Any]]],
    runs: Mapping[str, Beach | None],
) -> list[dict[str, Any]]:
    definitions = (
        (
            "archived_v1_3_to_new_finite",
            "archive_version_drift",
            "archived_v1_3",
            "new_finite_configured",
            "version_and_reproducibility_drift",
        ),
        (
            "new_finite_to_new_infinite",
            "boundary_history_response_common_evaluator",
            "new_finite_configured",
            "new_infinite_physical",
            "boundary_model_charging_response",
        ),
    )
    rows: list[dict[str, Any]] = []
    for label, kind, left_case, right_case, interpretation in definitions:
        left_run = runs.get(left_case)
        right_run = runs.get(right_case)
        if left_run is None or right_run is None:
            continue
        left_result = left_run.result
        right_result = right_run.result
        left_mesh_ids = np.asarray(left_result.mesh_ids, dtype=np.int64)
        geometry = _geometry_comparison(left_result, right_result)
        geometry_matches = geometry["status"] == "match"
        common = sorted(
            set(snapshots.get(left_case, {})) & set(snapshots.get(right_case, {})),
            key=lambda value: (value[0] == "final", value[1]),
        )
        for sample_key in common:
            left_snapshot = snapshots[left_case][sample_key]
            right_snapshot = snapshots[right_case][sample_key]
            for metric in ("processed_particles", "rel_change"):
                left_value = float(left_snapshot[metric])
                right_value = float(right_snapshot[metric])
                signed = right_value - left_value
                scale = max(abs(left_value), abs(right_value))
                rows.append(
                    {
                        "comparison": label,
                        "comparison_kind": kind,
                        "sample_kind": sample_key[0],
                        "resolved_batch": sample_key[1],
                        "scope": "run",
                        "mesh_id": "",
                        "left_snapshot_id": left_snapshot["snapshot_id"],
                        "right_snapshot_id": right_snapshot["snapshot_id"],
                        "metric": metric,
                        "left_value": left_value,
                        "right_value": right_value,
                        "signed_difference": signed,
                        "absolute_difference": abs(signed),
                        "relative_difference": (
                            0.0 if scale == 0.0 else abs(signed) / scale
                        ),
                        "status": "computed",
                        "interpretation": interpretation,
                    }
                )
            if not geometry_matches:
                rows.append(
                    {
                        "comparison": label,
                        "comparison_kind": kind,
                        "sample_kind": sample_key[0],
                        "resolved_batch": sample_key[1],
                        "metric": "element_charge_norms",
                        "left_value": geometry.get("max_coordinate_difference_m"),
                        "right_value": geometry.get("coordinate_tolerance_m"),
                        "status": "incomparable_geometry",
                        "interpretation": (
                            f"{interpretation}; geometry={geometry['status']}"
                        ),
                    }
                )
                continue
            left_charges = np.asarray(left_snapshot["charges"], dtype=np.float64)
            right_charges = np.asarray(right_snapshot["charges"], dtype=np.float64)
            scopes = [("run", "", np.ones(left_charges.size, dtype=bool))]
            scopes.extend(
                (
                    "object",
                    mesh_id,
                    left_mesh_ids == mesh_id,
                )
                for mesh_id in sorted(set(int(value) for value in left_mesh_ids))
            )
            for scope, mesh_id, mask in scopes:
                rows.extend(
                    _charge_vector_comparison_rows(
                        comparison=label,
                        comparison_kind=kind,
                        sample_kind=sample_key[0],
                        resolved_batch=sample_key[1],
                        scope=scope,
                        mesh_id=mesh_id,
                        left_snapshot_id=str(left_snapshot["snapshot_id"]),
                        right_snapshot_id=str(right_snapshot["snapshot_id"]),
                        left=left_charges[mask],
                        right=right_charges[mask],
                        interpretation=interpretation,
                    )
                )
    return rows


def _physics_comparison_rows(
    wrench_rows: Sequence[Mapping[str, Any]],
    path_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    wrench_index = {
        (
            str(row.get("case")),
            str(row.get("periodic_model")),
            int(row.get("resolved_batch", 0)),
            int(row.get("mesh_id", 0)),
        ): row
        for row in wrench_rows
        if row.get("status") == "available"
        and row.get("component") == "total_external"
        and str(row.get("mesh_id", "")).isdigit()
    }
    path_index = {
        (
            str(row.get("case")),
            str(row.get("periodic_model")),
            int(row.get("resolved_batch", 0)),
            int(row.get("mesh_id", 0)),
        ): row
        for row in path_rows
        if row.get("status") == "available" and str(row.get("mesh_id", "")).isdigit()
    }
    definitions = (
        (
            "archive_configured_to_infinite_same_charge",
            "frozen_field_override",
            ("archived_v1_3", "configured"),
            ("archived_v1_3", "infinite_physical"),
            "same archived charge snapshot; field closure only",
        ),
        (
            "new_finite_configured_to_infinite_same_charge",
            "frozen_field_override",
            ("new_finite_configured", "configured"),
            ("new_finite_configured", "infinite_physical"),
            "same new-finite charge snapshot; field closure only",
        ),
        (
            "archive_to_new_finite_configured",
            "archive_version_drift",
            ("archived_v1_3", "configured"),
            ("new_finite_configured", "configured"),
            "version and runtime drift under configured finite evaluation",
        ),
        (
            "new_finite_to_new_infinite_common_evaluator",
            "boundary_history_response_common_evaluator",
            ("new_finite_configured", "infinite_physical"),
            ("new_infinite_physical", "configured"),
            "charging-history response under the same cached infinite evaluator",
        ),
        (
            "new_finite_to_new_infinite_end_to_end",
            "actual_end_to_end",
            ("new_finite_configured", "configured"),
            ("new_infinite_physical", "configured"),
            "charging-history and field-closure effects combined",
        ),
    )
    rows: list[dict[str, Any]] = []
    wrench_metrics = (
        "force_x_N",
        "force_y_N",
        "force_z_N",
        "torque_x_Nm",
        "torque_y_Nm",
        "torque_z_Nm",
    )
    path_metrics = (
        "endpoint_work_J",
        "max_force_z_N",
        "minimum_available_energy_J",
        "minimum_energy_margin_over_tolerance",
        "barrier_decision_margin_J",
        "barrier_decision_margin_over_tolerance",
        "endpoint_available_energy_J",
        "endpoint_speed_m_s",
        "endpoint_speed_sensitivity_m_s",
        "maximum_reachable_speed_m_s",
        "instantaneous_force_margin_N",
        "force_margin_over_absolute_tolerance",
    )

    def append_scalar(
        *,
        label: str,
        kind: str,
        left: Mapping[str, Any],
        right: Mapping[str, Any],
        metric: str,
        batch: int,
        mesh_id: int,
        interpretation: str,
        status: str,
    ) -> None:
        left_value = float(left[metric])
        right_value = float(right[metric])
        signed = right_value - left_value
        scale = max(abs(left_value), abs(right_value))
        sample_kind = "final" if left.get("step_selector") == "final" else "history"
        rows.append(
            {
                "comparison": label,
                "comparison_kind": kind,
                "sample_kind": sample_kind,
                "resolved_batch": batch,
                "scope": "object",
                "mesh_id": mesh_id,
                "left_snapshot_id": (f"{left['case']}:{sample_kind}:{batch}"),
                "right_snapshot_id": (f"{right['case']}:{sample_kind}:{batch}"),
                "left_effective_far_correction": left.get(
                    "effective_far_correction", ""
                ),
                "right_effective_far_correction": right.get(
                    "effective_far_correction", ""
                ),
                "metric": metric,
                "left_value": left_value,
                "right_value": right_value,
                "signed_difference": signed,
                "absolute_difference": abs(signed),
                "relative_difference": (0.0 if scale == 0.0 else abs(signed) / scale),
                "status": status,
                "interpretation": interpretation,
            }
        )

    def append_categorical(
        *,
        label: str,
        kind: str,
        left: Mapping[str, Any],
        right: Mapping[str, Any],
        metric: str,
        batch: int,
        mesh_id: int,
        interpretation: str,
        status: str,
    ) -> None:
        sample_kind = "final" if left.get("step_selector") == "final" else "history"
        rows.append(
            {
                "comparison": label,
                "comparison_kind": kind,
                "sample_kind": sample_kind,
                "resolved_batch": batch,
                "scope": "object",
                "mesh_id": mesh_id,
                "left_snapshot_id": f"{left['case']}:{sample_kind}:{batch}",
                "right_snapshot_id": f"{right['case']}:{sample_kind}:{batch}",
                "left_effective_far_correction": left.get(
                    "effective_far_correction", ""
                ),
                "right_effective_far_correction": right.get(
                    "effective_far_correction", ""
                ),
                "metric": metric,
                "left_value": left.get(metric, ""),
                "right_value": right.get(metric, ""),
                "signed_difference": "",
                "absolute_difference": "",
                "relative_difference": "",
                "status": status,
                "interpretation": interpretation,
            }
        )

    for label, kind, left_model, right_model, interpretation in definitions:
        left_keys = {
            (batch, mesh_id): row
            for (case, model, batch, mesh_id), row in wrench_index.items()
            if (case, model) == left_model
        }
        right_keys = {
            (batch, mesh_id): row
            for (case, model, batch, mesh_id), row in wrench_index.items()
            if (case, model) == right_model
        }
        for batch, mesh_id in sorted(set(left_keys) & set(right_keys)):
            left = left_keys[(batch, mesh_id)]
            right = right_keys[(batch, mesh_id)]
            for metric in wrench_metrics:
                append_scalar(
                    label=label,
                    kind=kind,
                    left=left,
                    right=right,
                    metric=metric,
                    batch=batch,
                    mesh_id=mesh_id,
                    interpretation=interpretation,
                    status="computed",
                )
            left_path = path_index.get((*left_model, batch, mesh_id))
            right_path = path_index.get((*right_model, batch, mesh_id))
            if left_path is None or right_path is None:
                continue
            qualified = bool(left_path.get("numerically_qualified")) and bool(
                right_path.get("numerically_qualified")
            )
            endpoints_reachable = bool(
                left_path.get("endpoint_reachable_from_rest")
            ) and bool(right_path.get("endpoint_reachable_from_rest"))
            for metric in path_metrics:
                if not qualified:
                    status = "numerically_unqualified"
                elif metric.startswith("endpoint_speed"):
                    status = (
                        "conditional_local_release_proxy"
                        if endpoints_reachable
                        else "mechanically_inaccessible"
                    )
                else:
                    status = "computed"
                append_scalar(
                    label=label,
                    kind=kind,
                    left=left_path,
                    right=right_path,
                    metric=metric,
                    batch=batch,
                    mesh_id=mesh_id,
                    interpretation=interpretation,
                    status=status,
                )
            categorical_status = (
                "numerically_unqualified"
                if not qualified
                else "conditional_local_release_proxy"
                if endpoints_reachable
                else "mechanically_inaccessible"
            )
            for metric in (
                "barrier_free_from_rest",
                "endpoint_reachable_from_rest",
                "barrier_interpretation",
            ):
                append_categorical(
                    label=label,
                    kind=kind,
                    left=left_path,
                    right=right_path,
                    metric=metric,
                    batch=batch,
                    mesh_id=mesh_id,
                    interpretation=interpretation,
                    status=categorical_status,
                )
    return rows


def _comparison_artifact_contract(
    comparison_rows: Sequence[Mapping[str, Any]],
    snapshot_rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    required_kinds = {
        "archive_version_drift",
        "frozen_field_override",
        "boundary_history_response_common_evaluator",
        "actual_end_to_end",
    }
    actual_kinds = {
        str(row.get("comparison_kind"))
        for row in comparison_rows
        if row.get("comparison_kind")
    }
    missing_kinds = sorted(required_kinds - actual_kinds)
    unexpected_kinds = sorted(actual_kinds - required_kinds)
    required_metrics = {
        "force_z_N",
        "endpoint_work_J",
        "minimum_available_energy_J",
        "barrier_free_from_rest",
        "endpoint_reachable_from_rest",
    }
    missing_metrics: dict[str, list[str]] = {}
    for kind in sorted(required_kinds):
        metrics = {
            str(row.get("metric"))
            for row in comparison_rows
            if row.get("comparison_kind") == kind
        }
        missing = sorted(required_metrics - metrics)
        if missing:
            missing_metrics[kind] = missing
    known_snapshots = {
        str(row.get("snapshot_id"))
        for row in snapshot_rows
        if row.get("status") == "available" and row.get("snapshot_id")
    }
    unresolved_snapshots: list[dict[str, str]] = []
    semantic_errors: list[str] = []
    bad_status_rows: list[dict[str, Any]] = []
    numerical_status_rows: list[dict[str, Any]] = []
    expected_effective = {
        "frozen_field_override": ("none", "cached_kneq0"),
        "boundary_history_response_common_evaluator": (
            "cached_kneq0",
            "cached_kneq0",
        ),
        "actual_end_to_end": ("none", "cached_kneq0"),
    }
    for row in comparison_rows:
        kind = str(row.get("comparison_kind", ""))
        left_snapshot = str(row.get("left_snapshot_id", ""))
        right_snapshot = str(row.get("right_snapshot_id", ""))
        for side, snapshot_id in (
            ("left", left_snapshot),
            ("right", right_snapshot),
        ):
            if snapshot_id and snapshot_id not in known_snapshots:
                unresolved_snapshots.append(
                    {
                        "comparison": str(row.get("comparison", "")),
                        "side": side,
                        "snapshot_id": snapshot_id,
                    }
                )
        status = str(row.get("status", ""))
        if status in {
            "missing_comparison",
            "incomparable_geometry",
            "invalid",
            "not_evaluated",
        }:
            bad_status_rows.append(
                {
                    "comparison": row.get("comparison"),
                    "metric": row.get("metric"),
                    "status": status,
                }
            )
        if status == "numerically_unqualified":
            numerical_status_rows.append(
                {
                    "comparison": row.get("comparison"),
                    "metric": row.get("metric"),
                }
            )
        left_effective = str(row.get("left_effective_far_correction", ""))
        right_effective = str(row.get("right_effective_far_correction", ""))
        if kind in expected_effective and (left_effective or right_effective):
            if (left_effective, right_effective) != expected_effective[kind]:
                semantic_errors.append(
                    f"{kind} effective tuple is ({left_effective}, {right_effective})"
                )
        if (
            kind == "frozen_field_override"
            and left_snapshot
            and right_snapshot
            and left_snapshot != right_snapshot
        ):
            semantic_errors.append(
                "frozen_field_override compares different charge snapshots"
            )
    complete = not any(
        (
            missing_kinds,
            unexpected_kinds,
            missing_metrics,
            unresolved_snapshots,
            semantic_errors,
            bad_status_rows,
        )
    )
    return {
        "status": "complete" if complete else "incomplete",
        "required_comparison_kinds": sorted(required_kinds),
        "missing_comparison_kinds": missing_kinds,
        "unexpected_comparison_kinds": unexpected_kinds,
        "missing_metrics": missing_metrics,
        "unresolved_snapshots": unresolved_snapshots,
        "semantic_errors": semantic_errors,
        "bad_status_rows": bad_status_rows,
        "numerically_unqualified_rows": numerical_status_rows,
    }


def _placeholder_object_rows(
    case_rows: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    wrenches: list[dict[str, Any]] = []
    paths: list[dict[str, Any]] = []
    for row in case_rows:
        if row.get("status") == "available":
            status = "not_evaluated"
            message = "physics evaluation is deferred when target mesh/native kernel is unavailable"
        else:
            status = str(row.get("status"))
            message = str(row.get("message", ""))
        wrenches.append(
            {
                "case": row["case"],
                "periodic_model": "configured",
                "mesh_id": "",
                "component": "",
                "status": status,
                "message": message,
            }
        )
        paths.append(
            {
                "case": row["case"],
                "periodic_model": "configured",
                "mesh_id": "",
                "path_status": "not_evaluated",
                "speed_status": "contact_model_required",
                "status": status,
                "message": message,
            }
        )
    return wrenches, paths


def _load_json_if_exists(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _load_json_object_if_exists(path: Path, *, label: str) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValidationError(f"failed to read {label} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValidationError(f"{label} root must be an object: {path}")
    return value


def _toml_section(path: Path, section: str) -> dict[str, Any]:
    if not path.is_file():
        return {}
    parsed = _load_toml(path)
    value = parsed.get(section, {})
    if not isinstance(value, dict):
        raise ValidationError(f"TOML section [{section}] must be a table: {path}")
    return value


def _legacy_finite_number(value: Any, *, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"legacy estimator {label} is invalid") from exc
    if not math.isfinite(number):
        raise ValidationError(f"legacy estimator {label} is not finite")
    return number


def _legacy_integer(value: float, *, label: str) -> int:
    number = int(value)
    if float(number) != value:
        raise ValidationError(f"legacy estimator {label} is not an integer")
    return number


def _legacy_selected_csv_rows(
    path: Path,
    *,
    fields: Sequence[str],
    key_fields: Sequence[str],
    selector: Any,
    optional_fields: Sequence[str] = (),
) -> list[dict[str, float | None]]:
    selected: list[dict[str, float | None]] = []
    seen: set[tuple[float | None, ...]] = set()
    optional = set(optional_fields)
    try:
        with path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            if tuple(reader.fieldnames or ()) != tuple(fields):
                raise ValidationError(f"legacy estimator schema mismatch: {path.name}")
            for line_number, row in enumerate(reader, start=2):
                if None in row or any(row.get(field) is None for field in fields):
                    raise ValidationError(
                        f"legacy estimator row shape mismatch: "
                        f"{path.name}:{line_number}"
                    )
                parsed: dict[str, float | None] = {}
                for field in fields:
                    raw = str(row[field]).strip()
                    if field in optional and not raw:
                        parsed[field] = None
                    else:
                        parsed[field] = _legacy_finite_number(
                            raw,
                            label=f"{path.name}:{line_number}:{field}",
                        )
                key = tuple(parsed[field] for field in key_fields)
                if key in seen:
                    raise ValidationError(
                        f"legacy estimator duplicate key in {path.name}: {key}"
                    )
                seen.add(key)
                if selector(parsed):
                    selected.append(parsed)
    except OSError as exc:
        raise ValidationError(f"legacy estimator input cannot be read: {path}") from exc
    return selected


def _legacy_model_contract(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ValidationError(f"legacy estimator input is missing: {path}")
    value = _load_json_object_if_exists(path, label="legacy estimator model summary")
    if (
        value.get("status") != "ok"
        or value.get("model") != "moving_top_mesh_pairwise_coulomb"
        or Path(str(value.get("force_curves_csv", ""))).name
        != "moving_sphere_force_curves.csv"
        or Path(str(value.get("release_summary_csv", ""))).name
        != "moving_sphere_release_summary.csv"
    ):
        raise ValidationError("legacy estimator model summary contract is invalid")
    radius = _legacy_finite_number(value.get("radius_m"), label="model radius")
    z_max = _legacy_finite_number(value.get("z_max_m"), label="model z_max")
    softening = _legacy_finite_number(
        value.get("pair_softening_m"), label="model pair softening"
    )
    samples = _legacy_integer(
        _legacy_finite_number(value.get("z_samples"), label="model z_samples"),
        label="model z_samples",
    )
    if (
        radius <= 0.0
        or samples < 2
        or softening != 0.0
        or not math.isclose(z_max, 2.0 * radius, rel_tol=1.0e-12, abs_tol=1.0e-18)
    ):
        raise ValidationError("legacy estimator model geometry contract is invalid")
    return {"radius_m": radius, "z_max_m": z_max, "z_samples": samples}


def _legacy_comparison_row(
    *,
    comparison_kind: str,
    estimator: str,
    batch: int,
    mesh_id: int,
    displacement_m: float,
    component: str,
    quantity: str,
    legacy_value: float,
    current_value: float,
    interpretation: str,
) -> dict[str, Any]:
    signed = current_value - legacy_value
    scale = max(abs(legacy_value), abs(current_value))
    return {
        "comparison_kind": comparison_kind,
        "estimator": estimator,
        "legacy_self_policy": "exclude_target_entirely",
        "current_self_policy": "exclude_primary_keep_images",
        "batch": batch,
        "mesh_id": mesh_id,
        "displacement_m": displacement_m,
        "component": component,
        "quantity": quantity,
        "legacy_value": legacy_value,
        "current_value": current_value,
        "signed_difference": signed,
        "absolute_difference": abs(signed),
        "relative_difference": 0.0 if scale == 0.0 else abs(signed) / scale,
        "status": "computed",
        "interpretation": interpretation,
    }


def _legacy_estimator_audit_impl(
    archive_run: Path,
    wrench_rows: Sequence[Mapping[str, Any]],
    curve_rows: Sequence[Mapping[str, Any]],
    path_rows: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    paths = {
        Path(relative).name: archive_run / relative
        for relative in LEGACY_ESTIMATOR_INPUTS
    }
    missing = sorted(name for name, path in paths.items() if not path.is_file())
    if missing:
        raise ValidationError(
            "legacy estimator input set is incomplete: " + ", ".join(missing)
        )
    model = _legacy_model_contract(paths["moving_sphere_model_summary.json"])
    expected_keys = set(LEGACY_NATIVE_KEYS)

    def selected_curve(row: Mapping[str, float | None]) -> bool:
        key = (
            _legacy_integer(float(row["batch"]), label="curve batch"),
            _legacy_integer(float(row["target_mesh_id"]), label="curve mesh"),
        )
        return row["target_charge_multiplier"] == 1.0 and key in expected_keys

    def selected_release(row: Mapping[str, float | None]) -> bool:
        key = (
            _legacy_integer(float(row["batch"]), label="release batch"),
            _legacy_integer(float(row["target_mesh_id"]), label="release mesh"),
        )
        return row["target_charge_multiplier"] == 1.0 and key in expected_keys

    def selected_timeseries(row: Mapping[str, float | None]) -> bool:
        key = (
            _legacy_integer(float(row["batch"]), label="timeseries batch"),
            _legacy_integer(float(row["top_mesh_id"]), label="timeseries mesh"),
        )
        return row["local_charge_multiplier"] == 1.0 and key in expected_keys

    legacy_curves = _legacy_selected_csv_rows(
        paths["moving_sphere_force_curves.csv"],
        fields=LEGACY_CURVE_FIELDS,
        key_fields=(
            "target_mesh_id",
            "batch",
            "target_charge_multiplier",
            "displacement_z_m",
        ),
        selector=selected_curve,
    )
    legacy_releases = _legacy_selected_csv_rows(
        paths["moving_sphere_release_summary.csv"],
        fields=LEGACY_RELEASE_FIELDS,
        key_fields=("target_mesh_id", "batch", "target_charge_multiplier"),
        selector=selected_release,
        optional_fields=("first_crossing_z_m",),
    )
    legacy_timeseries = _legacy_selected_csv_rows(
        paths["force_timeseries.csv"],
        fields=LEGACY_TIMESERIES_FIELDS,
        key_fields=("top_mesh_id", "batch", "local_charge_multiplier"),
        selector=selected_timeseries,
    )

    release_index: dict[tuple[int, int], Mapping[str, float | None]] = {}
    for row in legacy_releases:
        key = (int(float(row["batch"])), int(float(row["target_mesh_id"])))
        release_index[key] = row
    timeseries_index: dict[tuple[int, int], Mapping[str, float | None]] = {}
    for row in legacy_timeseries:
        key = (int(float(row["batch"])), int(float(row["top_mesh_id"])))
        timeseries_index[key] = row
    curve_groups: dict[tuple[int, int], list[Mapping[str, float | None]]] = {
        key: [] for key in expected_keys
    }
    for row in legacy_curves:
        key = (int(float(row["batch"])), int(float(row["target_mesh_id"])))
        curve_groups[key].append(row)
    if set(release_index) != expected_keys or set(timeseries_index) != expected_keys:
        raise ValidationError(
            "legacy estimator selected batch/mesh coverage is invalid"
        )

    current_paths: dict[tuple[int, int], Mapping[str, Any]] = {}
    for row in path_rows:
        if (
            row.get("case") == "archived_v1_3"
            and row.get("periodic_model") == "configured"
            and row.get("status") == "available"
        ):
            key = (int(row.get("resolved_batch", 0)), int(row.get("mesh_id", 0)))
            if key in expected_keys:
                if key in current_paths:
                    raise ValidationError(
                        f"legacy estimator duplicate current path: {key}"
                    )
                current_paths[key] = row
    current_wrenches: dict[tuple[int, int], Mapping[str, Any]] = {}
    for row in wrench_rows:
        if (
            row.get("case") == "archived_v1_3"
            and row.get("periodic_model") == "configured"
            and row.get("component") == "total_external"
            and row.get("status") == "available"
        ):
            key = (int(row.get("resolved_batch", 0)), int(row.get("mesh_id", 0)))
            if key in expected_keys:
                if key in current_wrenches:
                    raise ValidationError(
                        f"legacy estimator duplicate current wrench: {key}"
                    )
                current_wrenches[key] = row
    if set(current_paths) != expected_keys or set(current_wrenches) != expected_keys:
        raise ValidationError("legacy estimator current native coverage is invalid")

    components = (
        "other_objects_all_images",
        "target_periodic_images",
        "total_external",
    )
    current_curves: dict[tuple[int, int, str, float], Mapping[str, Any]] = {}
    for row in curve_rows:
        if (
            row.get("case") != "archived_v1_3"
            or row.get("periodic_model") != "configured"
            or row.get("component") not in components
            or row.get("status") not in {"available", "converged"}
        ):
            continue
        key_base = (int(row.get("resolved_batch", 0)), int(row.get("mesh_id", 0)))
        if key_base not in expected_keys:
            continue
        radius = _legacy_finite_number(
            current_paths[key_base].get("radius_m"), label="current radius"
        )
        displacement = _legacy_finite_number(
            row.get("displacement_m"), label="current displacement"
        )
        endpoint: float | None = None
        for candidate in (0.0, 2.0 * radius):
            if math.isclose(displacement, candidate, rel_tol=1.0e-12, abs_tol=1.0e-18):
                endpoint = candidate
                break
        if endpoint is None:
            continue
        key = (*key_base, str(row["component"]), endpoint)
        if key in current_curves:
            raise ValidationError(f"legacy estimator duplicate current curve: {key}")
        current_curves[key] = row
    expected_current_curves = {
        (*key, component, displacement)
        for key in LEGACY_NATIVE_KEYS
        for component in components
        for displacement in (0.0, 2.0 * model["radius_m"])
    }
    if set(current_curves) != expected_current_curves:
        raise ValidationError("legacy estimator current curve coverage is invalid")

    comparison_rows: list[dict[str, Any]] = []
    for key in LEGACY_NATIVE_KEYS:
        batch, mesh_id = key
        path = current_paths[key]
        wrench = current_wrenches[key]
        release = release_index[key]
        timeseries = timeseries_index[key]
        radius = _legacy_finite_number(path.get("radius_m"), label="current radius")
        if not math.isclose(
            radius, model["radius_m"], rel_tol=1.0e-12, abs_tol=1.0e-18
        ):
            raise ValidationError(f"legacy estimator radius mapping mismatch: {key}")
        current_charge = _legacy_finite_number(
            wrench.get("total_charge_C"), label="current target charge"
        )
        legacy_charge = float(release["q_target_base_c"])
        timeseries_charge = float(timeseries["q_base_signed_c"])
        if not (
            math.isclose(
                current_charge, legacy_charge, rel_tol=1.0e-12, abs_tol=1.0e-30
            )
            and math.isclose(
                current_charge,
                timeseries_charge,
                rel_tol=1.0e-12,
                abs_tol=1.0e-30,
            )
        ):
            raise ValidationError(f"legacy estimator charge mapping mismatch: {key}")
        if (
            _legacy_integer(float(release["z_samples"]), label="release z_samples")
            != model["z_samples"]
            or not math.isclose(
                float(release["radius_m"]),
                model["radius_m"],
                rel_tol=1.0e-12,
                abs_tol=1.0e-18,
            )
            or not math.isclose(
                float(release["z_max_m"]),
                model["z_max_m"],
                rel_tol=1.0e-12,
                abs_tol=1.0e-18,
            )
            or not math.isclose(
                float(timeseries["radius_m"]),
                model["radius_m"],
                rel_tol=1.0e-12,
                abs_tol=1.0e-18,
            )
        ):
            raise ValidationError(f"legacy estimator summary geometry mismatch: {key}")
        group = curve_groups[key]
        if len(group) != model["z_samples"]:
            raise ValidationError(
                f"legacy estimator curve sample coverage mismatch: {key}"
            )
        legacy_endpoints: dict[float, Mapping[str, float | None]] = {}
        for row in group:
            displacement = float(row["displacement_z_m"])
            for endpoint in (0.0, model["z_max_m"]):
                if math.isclose(
                    displacement, endpoint, rel_tol=1.0e-12, abs_tol=1.0e-18
                ):
                    if endpoint in legacy_endpoints:
                        raise ValidationError(
                            f"legacy estimator duplicate curve endpoint: {key}"
                        )
                    legacy_endpoints[endpoint] = row
        if set(legacy_endpoints) != {0.0, model["z_max_m"]}:
            raise ValidationError(
                f"legacy estimator curve endpoints are missing: {key}"
            )
        if not (
            math.isclose(
                float(legacy_endpoints[0.0]["f_coulomb_z_n"]),
                float(release["f_coulomb_initial_z_n"]),
                rel_tol=1.0e-12,
                abs_tol=1.0e-30,
            )
            and math.isclose(
                float(legacy_endpoints[model["z_max_m"]]["f_coulomb_z_n"]),
                float(release["f_coulomb_final_z_n"]),
                rel_tol=1.0e-12,
                abs_tol=1.0e-30,
            )
        ):
            raise ValidationError(
                f"legacy estimator curve-summary endpoint mismatch: {key}"
            )
        ordered_curve = sorted(group, key=lambda row: float(row["displacement_z_m"]))
        displacements = [float(row["displacement_z_m"]) for row in ordered_curve]
        if any(right <= left for left, right in zip(displacements, displacements[1:])):
            raise ValidationError(
                f"legacy estimator curve displacement order is invalid: {key}"
            )
        integrated_work = sum(
            0.5
            * (float(left["f_coulomb_z_n"]) + float(right["f_coulomb_z_n"]))
            * (float(right["displacement_z_m"]) - float(left["displacement_z_m"]))
            for left, right in zip(ordered_curve, ordered_curve[1:])
        )
        if not math.isclose(
            integrated_work,
            float(release["w_coulomb_signed_j"]),
            rel_tol=1.0e-12,
            abs_tol=1.0e-30,
        ):
            raise ValidationError(
                f"legacy estimator curve-summary work mismatch: {key}"
            )

        current_total_z0 = _legacy_finite_number(
            current_curves[(*key, "total_external", 0.0)].get("force_z_N"),
            label="current total z0 force",
        )
        for estimator, field in (
            ("direct_object_z", "f_direct_object_z_n"),
            ("local_pair_proxy", "f_local_pair_proxy_n"),
        ):
            comparison_rows.append(
                _legacy_comparison_row(
                    comparison_kind="legacy_force_timeseries",
                    estimator=estimator,
                    batch=batch,
                    mesh_id=mesh_id,
                    displacement_m=0.0,
                    component="total_external",
                    quantity="force_z_N",
                    legacy_value=float(timeseries[field]),
                    current_value=current_total_z0,
                    interpretation=(
                        "legacy heuristic proxy; diagnostic difference only"
                        if estimator == "local_pair_proxy"
                        else "legacy direct object estimator versus current native finite"
                    ),
                )
            )
        for displacement in (0.0, model["z_max_m"]):
            legacy_force = float(legacy_endpoints[displacement]["f_coulomb_z_n"])
            for component in components:
                current_force = _legacy_finite_number(
                    current_curves[(*key, component, displacement)].get("force_z_N"),
                    label="current component force",
                )
                comparison_rows.append(
                    _legacy_comparison_row(
                        comparison_kind="legacy_moving_sphere_force",
                        estimator="moving_top_mesh_pairwise_coulomb",
                        batch=batch,
                        mesh_id=mesh_id,
                        displacement_m=displacement,
                        component=component,
                        quantity="force_z_N",
                        legacy_value=legacy_force,
                        current_value=current_force,
                        interpretation=(
                            "self policies differ; diagnostic decomposition, not a closeness gate"
                        ),
                    )
                )
        comparison_rows.append(
            _legacy_comparison_row(
                comparison_kind="legacy_moving_sphere_work",
                estimator="moving_top_mesh_pairwise_coulomb",
                batch=batch,
                mesh_id=mesh_id,
                displacement_m=model["z_max_m"],
                component="total_external",
                quantity="electrostatic_work_J",
                legacy_value=float(release["w_coulomb_signed_j"]),
                current_value=_legacy_finite_number(
                    path.get("endpoint_work_J"), label="current endpoint work"
                ),
                interpretation=(
                    "electrostatic work only; resistance and speed models are not equated"
                ),
            )
        )

    expected_row_count = len(LEGACY_NATIVE_KEYS) * 9
    if len(comparison_rows) != expected_row_count:
        raise ValidationError("legacy estimator comparison coverage is incomplete")
    return comparison_rows, {
        "status": "complete",
        "expected_native_keys": [list(key) for key in LEGACY_NATIVE_KEYS],
        "covered_native_keys": [list(key) for key in LEGACY_NATIVE_KEYS],
        "comparison_row_count": len(comparison_rows),
        "closeness_is_a_gate": False,
        "input_sha256": {name: _sha256(path) for name, path in sorted(paths.items())},
    }


def _legacy_estimator_audit(
    archive_run: Path,
    wrench_rows: Sequence[Mapping[str, Any]],
    curve_rows: Sequence[Mapping[str, Any]],
    path_rows: Sequence[Mapping[str, Any]],
    *,
    strict: bool,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    try:
        return _legacy_estimator_audit_impl(
            archive_run,
            wrench_rows,
            curve_rows,
            path_rows,
        )
    except ValidationError as exc:
        if strict:
            raise
        message = str(exc)
        return [
            {
                "status": "unavailable",
                "interpretation": message,
            }
        ], {
            "status": "unavailable",
            "expected_native_keys": [list(key) for key in LEGACY_NATIVE_KEYS],
            "covered_native_keys": [],
            "comparison_row_count": 0,
            "closeness_is_a_gate": False,
            "message": message,
        }


def _release_parameters(
    archive_run: Path,
    geometry_radius_m: float,
    *,
    allow_geometry_radius_fallback: bool = True,
    radius_relative_tolerance: float = TARGET_GEOMETRY_RADIUS_RELATIVE_TOLERANCE,
) -> dict[str, Any]:
    release_summary = _load_json_object_if_exists(
        archive_run / "analysis/local_release/release_model_summary.json",
        label="release model summary",
    )
    assumptions_value = release_summary.get("assumptions", {})
    if not isinstance(assumptions_value, dict):
        raise ValidationError("release model summary assumptions must be an object")
    assumptions = assumptions_value
    if not math.isfinite(geometry_radius_m) or geometry_radius_m <= 0.0:
        raise ValidationError("target geometry radius must be finite and positive")
    if not math.isfinite(radius_relative_tolerance) or radius_relative_tolerance < 0.0:
        raise ValidationError(
            "radius relative tolerance must be finite and nonnegative"
        )
    if "radius_m" not in assumptions:
        if not allow_geometry_radius_fallback:
            raise ValidationError(
                "declared dust radius is missing from "
                "release_model_summary.json assumptions.radius_m"
            )
        radius_m = float(geometry_radius_m)
        radius_source = "probe.vertex_bounding_radius_m_fallback"
    else:
        declared_radius = assumptions["radius_m"]
        if isinstance(declared_radius, bool):
            raise ValidationError("declared dust radius must be finite and positive")
        try:
            radius_m = float(declared_radius)
        except (TypeError, ValueError) as exc:
            raise ValidationError(
                "declared dust radius must be finite and positive"
            ) from exc
        if not math.isfinite(radius_m) or radius_m <= 0.0:
            raise ValidationError("declared dust radius must be finite and positive")
        radius_source = "release_model_summary.assumptions.radius_m"
    radius_relative_mismatch = abs(geometry_radius_m - radius_m) / radius_m
    if radius_relative_mismatch > radius_relative_tolerance:
        raise ValidationError(
            "declared/geometry radius mismatch exceeds tolerance: "
            f"declared={radius_m:.17g} m, geometry={geometry_radius_m:.17g} m, "
            f"relative_mismatch={radius_relative_mismatch:.17g}, "
            f"tolerance={radius_relative_tolerance:.17g}"
        )
    density = float(assumptions.get("dust_density_kg_m3", 3000.0))
    gravity = float(assumptions.get("moon_gravity_m_s2", 1.62))
    eta_translation_sensitivity = float(assumptions.get("energy_partition", 0.5))
    eta_translation = 1.0
    energy_tolerance = 1.0e-18
    if not math.isfinite(density) or density <= 0.0:
        raise ValidationError("dust density must be finite and positive")
    if not math.isfinite(gravity) or gravity < 0.0:
        raise ValidationError("gravity must be finite and nonnegative")
    if (
        not math.isfinite(eta_translation_sensitivity)
        or not 0.0 <= eta_translation_sensitivity <= 1.0
    ):
        raise ValidationError("energy_partition must lie in [0, 1]")
    mass = 4.0 * math.pi * density * radius_m**3 / 3.0
    adhesion_config = _toml_section(
        archive_run / "input/release_kernel_base.toml", "adhesion"
    )
    adhesion: AdhesionProfile = AdhesionProfile.none()
    adhesion_force = 0.0
    adhesion_work = 0.0
    adhesion_equivalent_range = 0.0
    adhesion_model = "none"
    model = str(adhesion_config.get("model", "none")).strip().lower()
    if model not in {"none", "vdw_work"}:
        raise ValidationError(f"unsupported adhesion model: {model!r}")
    try:
        if model == "vdw_work":
            hamaker = float(adhesion_config["hamaker_constant"])
            contact = float(adhesion_config["contact_distance"])
            cutoff = float(adhesion_config["cutoff_distance"])
            roughness = float(adhesion_config["roughness_factor"])
            contacts = float(adhesion_config["contact_count"])
            peel = float(adhesion_config["peel_factor"])
            geometry = str(adhesion_config.get("contact_geometry", "sphere_sphere"))
            values = {
                "hamaker_constant": hamaker,
                "contact_distance": contact,
                "cutoff_distance": cutoff,
                "roughness_factor": roughness,
                "contact_count": contacts,
                "peel_factor": peel,
            }
            invalid = [
                name
                for name, value in values.items()
                if not math.isfinite(value) or value <= 0.0
            ]
            if invalid:
                raise ValidationError(
                    "invalid vdw adhesion positive parameters: " + ", ".join(invalid)
                )
            if cutoff <= contact:
                raise ValidationError(
                    "invalid vdw adhesion: cutoff_distance must exceed contact_distance"
                )
            if geometry not in {"sphere_sphere", "sphere_plane"}:
                raise ValidationError(
                    f"invalid vdw adhesion contact_geometry: {geometry!r}"
                )
            effective_radius = (
                radius_m / 2.0 if geometry == "sphere_sphere" else radius_m
            )
            coefficient = contacts * roughness * hamaker * effective_radius / 6.0
            adhesion_force = coefficient / contact**2
            adhesion_work = peel * coefficient * (1.0 / contact - 1.0 / cutoff)
            equivalent_range = adhesion_work / adhesion_force
            if not all(
                math.isfinite(value) and value > 0.0
                for value in (adhesion_force, adhesion_work, equivalent_range)
            ):
                raise ValidationError("invalid vdw adhesion derived force/work/range")
            adhesion = AdhesionProfile.finite_range_constant(
                adhesion_force, equivalent_range
            )
            adhesion_equivalent_range = equivalent_range
            adhesion_model = (
                "vdw_work_equivalent_constant_"
                "initial_force_and_total_work_preserving_not_1_over_s2"
            )
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValidationError(f"invalid vdw adhesion configuration: {exc}") from exc
    return {
        "radius_m": radius_m,
        "radius_source": radius_source,
        "geometry_radius_m": float(geometry_radius_m),
        "radius_relative_mismatch": radius_relative_mismatch,
        "radius_relative_tolerance": float(radius_relative_tolerance),
        "density_kg_m3": density,
        "mass_kg": mass,
        "gravity_m_s2": gravity,
        "eta_translation": eta_translation,
        "eta_translation_sensitivity": eta_translation_sensitivity,
        "energy_tolerance_J": energy_tolerance,
        "adhesion": adhesion,
        "adhesion_model": adhesion_model,
        "adhesion_force_N": adhesion_force,
        "adhesion_work_J": adhesion_work,
        "adhesion_equivalent_range_m": adhesion_equivalent_range,
    }


def _probe_geometry_contract(probe: Any) -> tuple[float, str]:
    representation = str(getattr(probe, "target_geometry_representation", ""))
    if representation != TARGET_GEOMETRY_REPRESENTATION:
        raise ValidationError(
            "target geometry representation mismatch: expected "
            f"{TARGET_GEOMETRY_REPRESENTATION!r}, got {representation!r}"
        )
    raw_radius = getattr(probe, "vertex_bounding_radius_m", None)
    if isinstance(raw_radius, (bool, np.bool_)) or np.asarray(raw_radius).ndim != 0:
        raise ValidationError("probe vertex bounding radius must be a scalar")
    try:
        radius = float(raw_radius)
    except (AttributeError, TypeError, ValueError) as exc:
        raise ValidationError("probe vertex bounding radius is unavailable") from exc
    if not math.isfinite(radius) or radius <= 0.0:
        raise ValidationError(
            "probe vertex bounding radius must be finite and positive"
        )
    return radius, representation


def _validate_geometry_metadata(
    metadata: Mapping[str, Any],
    *,
    probe_radius_m: float,
) -> None:
    representation = str(metadata.get("target_geometry_representation", ""))
    if representation != TARGET_GEOMETRY_REPRESENTATION:
        raise ValidationError(
            "wrench/path target geometry representation mismatch: expected "
            f"{TARGET_GEOMETRY_REPRESENTATION!r}, got {representation!r}"
        )
    raw_radius = metadata.get("vertex_bounding_radius_m")
    if isinstance(raw_radius, bool):
        raise ValidationError("wrench/path metadata radius is invalid")
    try:
        metadata_radius = float(raw_radius)
    except (TypeError, ValueError) as exc:
        raise ValidationError("wrench/path metadata radius is invalid") from exc
    if not math.isfinite(metadata_radius) or not math.isclose(
        metadata_radius,
        probe_radius_m,
        rel_tol=1.0e-12,
        abs_tol=0.0,
    ):
        raise ValidationError(
            "wrench/path metadata radius does not match probe vertex bounding radius"
        )


def _torque_origin_provenance(
    metadata: Mapping[str, Any],
    *,
    point_count: int | None = None,
) -> tuple[str, np.ndarray]:
    policy = str(metadata.get("torque_origin_policy", ""))
    if policy not in {
        "moving_geometric_area_centroid",
        "fixed_origin",
        "fixed_explicit",
    }:
        raise ValidationError(f"invalid torque origin policy: {policy!r}")
    try:
        origin = np.asarray(metadata.get("torque_origin_m"), dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValidationError("torque origin metadata is not numeric") from exc
    expected_shape = (3,) if point_count is None else (point_count, 3)
    if origin.shape != expected_shape or not np.all(np.isfinite(origin)):
        raise ValidationError(
            f"torque origin metadata must have shape {expected_shape} and be finite"
        )
    return policy, origin


def _validate_wrench_torque_origin(wrench: Any, metadata_origin_m: np.ndarray) -> None:
    try:
        actual_origin = np.asarray(wrench.torque_origin_m, dtype=np.float64)
    except (AttributeError, TypeError, ValueError) as exc:
        raise ValidationError("wrench torque origin is unavailable") from exc
    if actual_origin.shape != (3,) or not np.all(np.isfinite(actual_origin)):
        raise ValidationError("wrench torque origin must have shape (3,) and be finite")
    if not np.allclose(actual_origin, metadata_origin_m, rtol=1.0e-12, atol=0.0):
        raise ValidationError(
            "wrench metadata torque origin does not match wrench.torque_origin_m"
        )


def _validate_production_torque_origin_contract(
    *,
    probe: Any,
    wrench_policy: str,
    wrench_origin_m: Any,
    path_policy: str,
    path_origin_m: Any,
    displacement_m: Any,
) -> None:
    if (
        wrench_policy != "moving_geometric_area_centroid"
        or path_policy != "moving_geometric_area_centroid"
    ):
        raise ValidationError(
            "production torque origin policy must be moving_geometric_area_centroid"
        )
    try:
        centroid = np.asarray(probe.geometric_area_centroid_m, dtype=np.float64)
        wrench_origin = np.asarray(wrench_origin_m, dtype=np.float64)
        path_origins = np.asarray(path_origin_m, dtype=np.float64)
        displacement = np.asarray(displacement_m, dtype=np.float64)
    except (AttributeError, TypeError, ValueError) as exc:
        raise ValidationError(
            "production torque origin contract is unavailable"
        ) from exc
    if centroid.shape != (3,) or wrench_origin.shape != (3,):
        raise ValidationError("production wrench/centroid origins must have shape (3,)")
    if displacement.ndim != 1 or path_origins.shape != (displacement.size, 3):
        raise ValidationError("production path torque origins have the wrong shape")
    if not all(
        np.all(np.isfinite(value))
        for value in (centroid, wrench_origin, path_origins, displacement)
    ):
        raise ValidationError("production torque origin contract must be finite")
    expected = np.broadcast_to(centroid, path_origins.shape).copy()
    expected[:, 2] += displacement
    if not np.allclose(wrench_origin, centroid, rtol=1.0e-12, atol=1.0e-15):
        raise ValidationError(
            "production wrench torque origin is not the target centroid"
        )
    if not np.allclose(path_origins, expected, rtol=1.0e-12, atol=1.0e-15):
        raise ValidationError(
            "production moving torque origins do not track the displacement grid"
        )
    if not np.allclose(path_origins[0], wrench_origin, rtol=1.0e-12, atol=1.0e-15):
        raise ValidationError("wrench/path initial torque origins do not match")


def _component_array(value: Any, *, shape: tuple[int, ...], label: str) -> np.ndarray:
    try:
        result = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"{label} is not numeric") from exc
    if result.shape != shape or not np.all(np.isfinite(result)):
        raise ValidationError(f"{label} must have shape {shape} and be finite")
    return result


def _require_additive_identity(
    total: np.ndarray,
    additive: Sequence[np.ndarray],
    *,
    label: str,
) -> None:
    try:
        total_values = np.asarray(total, dtype=np.float64)
        stacked = np.stack(additive, axis=0).astype(np.float64, copy=False)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"{label} additive components are invalid") from exc
    if (
        stacked.shape[1:] != total_values.shape
        or not np.all(np.isfinite(total_values))
        or not np.all(np.isfinite(stacked))
    ):
        raise ValidationError(f"{label} additive components are invalid")

    expected = np.empty_like(total_values)
    for index in np.ndindex(total_values.shape):
        values = [
            float(stacked[(component, *index)]) for component in range(stacked.shape[0])
        ]
        magnitude = max((abs(value) for value in values), default=0.0)
        if magnitude == 0.0:
            summed = 0.0
        else:
            summed = magnitude * math.fsum(value / magnitude for value in values)
        if not math.isfinite(summed):
            raise ValidationError(f"{label} additive components have a non-finite sum")
        expected[index] = summed

    scale = np.maximum.reduce(
        (
            np.abs(total_values),
            np.abs(expected),
            np.full_like(total_values, np.finfo(float).tiny),
        )
    )
    normalized_delta = np.abs(total_values / scale - expected / scale)
    if np.any(normalized_delta > 1.0e-12):
        raise ValidationError(f"{label} total is not the sum of additive components")


def _validate_strict_wrench_component_contract(
    components: Mapping[str, Any],
    metadata: Mapping[str, Any],
    *,
    effective_far_correction: str,
    wrench_force_N: Any,
    wrench_torque_Nm: Any,
) -> None:
    if set(components) != set(PHYSICAL_OBJECT_COMPONENTS):
        raise ValidationError(
            "strict physical wrench component set must be exactly "
            + ", ".join(PHYSICAL_OBJECT_COMPONENTS)
        )
    force: dict[str, np.ndarray] = {}
    torque: dict[str, np.ndarray] = {}
    energy: dict[str, float] = {}
    for name in PHYSICAL_OBJECT_COMPONENTS:
        component = components[name]
        force[name] = _component_array(
            component.force_N, shape=(3,), label=f"{name} force"
        )
        torque[name] = _component_array(
            component.torque_Nm, shape=(3,), label=f"{name} torque"
        )
        try:
            energy[name] = float(component.potential_energy_J)
        except (AttributeError, TypeError, ValueError) as exc:
            raise ValidationError(f"{name} potential energy is invalid") from exc
        if not math.isfinite(energy[name]):
            raise ValidationError(f"{name} potential energy is invalid")
    _require_additive_identity(
        force["total_external"],
        [force[name] for name in ADDITIVE_OBJECT_COMPONENTS],
        label="additive physical wrench force",
    )
    _require_additive_identity(
        torque["total_external"],
        [torque[name] for name in ADDITIVE_OBJECT_COMPONENTS],
        label="additive physical wrench torque",
    )
    _require_additive_identity(
        np.asarray([energy["total_external"]]),
        [np.asarray([energy[name]]) for name in ADDITIVE_OBJECT_COMPONENTS],
        label="additive physical wrench potential energy",
    )
    if not np.array_equal(
        _component_array(wrench_force_N, shape=(3,), label="wrench total force"),
        force["total_external"],
    ) or not np.array_equal(
        _component_array(wrench_torque_Nm, shape=(3,), label="wrench total torque"),
        torque["total_external"],
    ):
        raise ValidationError(
            "wrench total arrays do not match total_external components"
        )
    available_numerical = {
        name
        for name in CACHED_NUMERICAL_COMPONENTS
        if isinstance(metadata.get(name), Mapping)
    }
    if effective_far_correction == "cached_kneq0":
        if available_numerical != set(CACHED_NUMERICAL_COMPONENTS):
            raise ValidationError(
                "strict cached numerical wrench component set must be exactly "
                + ", ".join(CACHED_NUMERICAL_COMPONENTS)
            )
        required_numerical = CACHED_NUMERICAL_COMPONENTS
    else:
        if available_numerical != {"primary_free_subtraction"}:
            raise ValidationError(
                "strict finite numerical wrench component set must be exactly "
                "primary_free_subtraction"
            )
        required_numerical = ("primary_free_subtraction",)
    for name in required_numerical:
        component = metadata[name]
        _component_array(
            component.get("force_N"),
            shape=(3,),
            label=f"numerical {name} force",
        )
        _component_array(
            component.get("torque_Nm"),
            shape=(3,),
            label=f"numerical {name} torque",
        )
        try:
            potential = float(component.get("potential_energy_J"))
        except (TypeError, ValueError) as exc:
            raise ValidationError(
                f"numerical {name} potential energy is invalid"
            ) from exc
        if not math.isfinite(potential):
            raise ValidationError(f"numerical {name} potential energy is invalid")


def _validate_strict_path_component_contract(
    path: Any,
    *,
    expected_start_m: float,
    expected_end_m: float,
) -> None:
    try:
        displacement = np.asarray(path.displacement_m, dtype=np.float64)
        expected_start = float(expected_start_m)
        expected_end = float(expected_end_m)
    except (AttributeError, TypeError, ValueError) as exc:
        raise ValidationError("path displacement contract is unavailable") from exc
    if (
        displacement.ndim != 1
        or displacement.size < 2
        or not np.all(np.isfinite(displacement))
        or not math.isfinite(expected_start)
        or not math.isfinite(expected_end)
        or expected_end <= expected_start
        or np.any(np.diff(displacement) <= 0.0)
        or not math.isclose(
            float(displacement[0]),
            expected_start,
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        )
        or not math.isclose(
            float(displacement[-1]),
            expected_end,
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        )
    ):
        raise ValidationError(
            "path displacement must be finite, strictly increasing, and preserve the required endpoints"
        )
    force_mapping = path.component_force_N
    torque_mapping = path.component_torque_Nm
    if set(force_mapping) != set(PHYSICAL_OBJECT_COMPONENTS):
        raise ValidationError(
            "strict path force component set must be exactly "
            + ", ".join(PHYSICAL_OBJECT_COMPONENTS)
        )
    if set(torque_mapping) != set(PHYSICAL_OBJECT_COMPONENTS):
        raise ValidationError(
            "strict path torque component set must be exactly "
            + ", ".join(PHYSICAL_OBJECT_COMPONENTS)
        )
    npoint = displacement.size
    shape = (npoint, 3)
    force = {
        name: _component_array(
            force_mapping[name], shape=shape, label=f"{name} path force"
        )
        for name in PHYSICAL_OBJECT_COMPONENTS
    }
    torque = {
        name: _component_array(
            torque_mapping[name], shape=shape, label=f"{name} path torque"
        )
        for name in PHYSICAL_OBJECT_COMPONENTS
    }
    _require_additive_identity(
        force["total_external"],
        [force[name] for name in ADDITIVE_OBJECT_COMPONENTS],
        label="additive physical path force",
    )
    _require_additive_identity(
        torque["total_external"],
        [torque[name] for name in ADDITIVE_OBJECT_COMPONENTS],
        label="additive physical path torque",
    )
    if not np.array_equal(
        _component_array(path.force_N, shape=shape, label="path total force"),
        force["total_external"],
    ) or not np.array_equal(
        _component_array(path.torque_Nm, shape=shape, label="path total torque"),
        torque["total_external"],
    ):
        raise ValidationError(
            "path total arrays do not match total_external components"
        )


def _validate_shell_reference_contract(shell: Any) -> None:
    try:
        layers = np.asarray(shell.image_layers)
        increment = np.asarray(shell.increment_converged)
        reference = np.asarray(shell.reference_converged)
        reference_force = np.asarray(shell.reference_force_error_N, dtype=np.float64)
        reference_work = np.asarray(shell.reference_work_error_J, dtype=np.float64)
        force_increment = np.asarray(shell.force_increment_error_N, dtype=np.float64)
        work_increment = np.asarray(shell.work_increment_error_J, dtype=np.float64)
        force_tail = np.asarray(shell.force_tail_proxy_N, dtype=np.float64)
        work_tail = np.asarray(shell.work_tail_proxy_J, dtype=np.float64)
    except (AttributeError, TypeError, ValueError) as exc:
        raise ValidationError("finite shell reference contract is unavailable") from exc
    if (
        layers.ndim != 1
        or layers.size == 0
        or not np.issubdtype(layers.dtype, np.integer)
    ):
        raise ValidationError(
            "finite shell image_layers must be consecutive and unique"
        )
    if int(layers[0]) != 0:
        raise ValidationError("finite shell image_layers must start at zero")
    if np.any(np.diff(layers) != 1):
        raise ValidationError(
            "finite shell image_layers must be consecutive and unique"
        )
    if increment.dtype != np.dtype(bool) or reference.dtype != np.dtype(bool):
        raise ValidationError(
            "finite shell convergence gate arrays must have boolean dtype"
        )
    if increment.shape != (layers.size - 1,):
        raise ValidationError("finite shell increment_converged has the wrong shape")
    increment_shape = (layers.size - 1,)
    increment_errors = (force_increment, work_increment, force_tail, work_tail)
    if any(value.shape != increment_shape for value in increment_errors) or any(
        not np.all(np.isfinite(value)) or np.any(value < 0.0)
        for value in increment_errors
    ):
        raise ValidationError(
            "finite shell increment/tail errors must have the right shape and be nonnegative finite"
        )
    if str(getattr(shell, "reference_model", "")) != "infinite_physical":
        raise ValidationError("finite shell reference_model must be infinite_physical")
    if (
        reference.shape != layers.shape
        or reference_force.shape != layers.shape
        or reference_work.shape != layers.shape
        or not np.all(np.isfinite(reference_force))
        or not np.all(np.isfinite(reference_work))
    ):
        raise ValidationError(
            "finite shell physical reference arrays have the wrong shape"
        )
    if np.any(reference_force < 0.0) or np.any(reference_work < 0.0):
        raise ValidationError(
            "finite shell physical reference errors must be nonnegative"
        )
    selected = getattr(shell, "selected_image_layers", None)
    status = str(getattr(shell, "status", ""))
    if status not in {"converged", "not_converged"}:
        raise ValidationError("finite shell status is invalid")
    if status == "converged":
        if selected is None or int(selected) != int(layers[-1]) or layers.size < 3:
            raise ValidationError("finite shell selected_image_layers is inconsistent")
        if not np.all(increment[-2:]):
            raise ValidationError(
                "finite shell selected layer requires two successive increment gates"
            )
        if not np.all(reference[-2:]):
            raise ValidationError(
                "finite shell selected layer requires two successive physical reference gates"
            )
    elif selected is not None:
        raise ValidationError("non-converged finite shell cannot select an image layer")


def _radius_provenance_fields(
    mechanics: Mapping[str, Any],
    target_geometry_representation: str,
) -> dict[str, Any]:
    return {
        "model_radius_m": mechanics["radius_m"],
        "radius_source": mechanics["radius_source"],
        "geometry_radius_m": mechanics["geometry_radius_m"],
        "radius_relative_mismatch": mechanics["radius_relative_mismatch"],
        "radius_relative_tolerance": mechanics["radius_relative_tolerance"],
        "target_geometry_representation": target_geometry_representation,
    }


def _wrench_row(
    *,
    case: str,
    periodic_model: str,
    effective_far_correction: str,
    zero_mode_policy: str,
    lower_boundary_model: str,
    periodic_cache: Mapping[str, Any] | None,
    step_selector: int | str,
    resolved_batch: int,
    mesh_id: int,
    mechanics: Mapping[str, Any],
    target_geometry_representation: str,
    component: str,
    component_kind: str,
    force: Any,
    torque: Any,
    torque_origin_policy: str,
    torque_origin_m: Any,
    total_charge: float,
    potential_energy: Any,
) -> dict[str, Any]:
    force_value = np.asarray(force, dtype=np.float64)
    torque_value = np.asarray(torque, dtype=np.float64)
    origin_value = np.asarray(torque_origin_m, dtype=np.float64)
    if force_value.shape != (3,) or torque_value.shape != (3,):
        raise ValidationError("wrench component vectors must have shape (3,)")
    if not np.all(np.isfinite(force_value)) or not np.all(np.isfinite(torque_value)):
        raise ValidationError("wrench component vectors must be finite")
    if origin_value.shape != (3,) or not np.all(np.isfinite(origin_value)):
        raise ValidationError("wrench torque origin must have shape (3,) and be finite")
    try:
        charge_value = float(total_charge)
        potential_value = None if potential_energy is None else float(potential_energy)
    except (TypeError, ValueError) as exc:
        raise ValidationError("wrench scalar fields must be finite") from exc
    if not math.isfinite(charge_value) or (
        potential_value is not None and not math.isfinite(potential_value)
    ):
        raise ValidationError("wrench scalar fields must be finite")
    return {
        "case": case,
        "periodic_model": periodic_model,
        "effective_far_correction": effective_far_correction,
        "zero_mode_policy": zero_mode_policy,
        "lower_boundary_model": lower_boundary_model,
        "periodic_cache_hit": ("" if periodic_cache is None else periodic_cache["hit"]),
        "periodic_cache_build_count": (
            "" if periodic_cache is None else periodic_cache["build_count"]
        ),
        "periodic_cache_fingerprint": (
            "" if periodic_cache is None else periodic_cache["fingerprint"]
        ),
        "periodic_cache_path": (
            "" if periodic_cache is None else periodic_cache["path"]
        ),
        "periodic_cache_path_sha256": (
            "" if periodic_cache is None else periodic_cache["path_sha256"]
        ),
        "periodic_cache_prime_receipt_sha256": (
            ""
            if periodic_cache is None
            else periodic_cache.get("prime_receipt_sha256", "")
        ),
        "step_selector": step_selector,
        "resolved_batch": resolved_batch,
        "mesh_id": mesh_id,
        **_radius_provenance_fields(mechanics, target_geometry_representation),
        "component": component,
        "component_kind": component_kind,
        "force_x_N": force_value[0],
        "force_y_N": force_value[1],
        "force_z_N": force_value[2],
        "torque_x_Nm": torque_value[0],
        "torque_y_Nm": torque_value[1],
        "torque_z_Nm": torque_value[2],
        "torque_origin_policy": torque_origin_policy,
        "torque_origin_x_m": origin_value[0],
        "torque_origin_y_m": origin_value[1],
        "torque_origin_z_m": origin_value[2],
        "potential_energy_J": "" if potential_value is None else potential_value,
        "total_charge_C": charge_value,
        "status": "available",
        "message": "",
    }


def _finite_shell_rows(
    *,
    case: str,
    periodic_model: str,
    effective_far_correction: str,
    zero_mode_policy: str,
    lower_boundary_model: str,
    step_selector: int | str,
    resolved_batch: int,
    mesh_id: int,
    shell: Any,
    mechanics: Mapping[str, Any],
    target_geometry_representation: str,
    displacement_m: Any,
) -> list[dict[str, Any]]:
    _validate_shell_reference_contract(shell)
    path_grid = np.asarray(displacement_m, dtype=np.float64)
    if path_grid.ndim != 1 or path_grid.size < 2 or not np.all(np.isfinite(path_grid)):
        raise ValidationError("finite shell displacement grid is invalid")
    if not math.isclose(
        float(path_grid[0]), 0.0, rel_tol=0.0, abs_tol=0.0
    ) or not math.isclose(
        float(path_grid[-1]),
        2.0 * float(mechanics["radius_m"]),
        rel_tol=1.0e-12,
        abs_tol=0.0,
    ):
        raise ValidationError("finite shell displacement grid must span 0..2R")
    rows: list[dict[str, Any]] = []
    for index, layer in enumerate(shell.image_layers):
        rows.append(
            {
                "case": case,
                "periodic_model": periodic_model,
                "effective_far_correction": effective_far_correction,
                "zero_mode_policy": zero_mode_policy,
                "lower_boundary_model": lower_boundary_model,
                "step_selector": step_selector,
                "resolved_batch": resolved_batch,
                "mesh_id": mesh_id,
                **_radius_provenance_fields(mechanics, target_geometry_representation),
                "path_start_m": path_grid[0],
                "path_end_m": path_grid[-1],
                "image_layers": int(layer),
                "force_increment_error_N": (
                    "" if index == 0 else shell.force_increment_error_N[index - 1]
                ),
                "work_increment_error_J": (
                    "" if index == 0 else shell.work_increment_error_J[index - 1]
                ),
                "force_tail_proxy_N": (
                    "" if index == 0 else shell.force_tail_proxy_N[index - 1]
                ),
                "work_tail_proxy_J": (
                    "" if index == 0 else shell.work_tail_proxy_J[index - 1]
                ),
                "converged": (
                    "" if index == 0 else bool(shell.increment_converged[index - 1])
                ),
                "reference_model": (
                    "" if shell.reference_model is None else shell.reference_model
                ),
                "reference_force_error_N": (
                    ""
                    if shell.reference_force_error_N is None
                    else shell.reference_force_error_N[index]
                ),
                "reference_work_error_J": (
                    ""
                    if shell.reference_work_error_J is None
                    else shell.reference_work_error_J[index]
                ),
                "reference_converged": (
                    ""
                    if shell.reference_converged is None
                    else bool(shell.reference_converged[index])
                ),
                "selected_image_layers": (
                    ""
                    if shell.selected_image_layers is None
                    else shell.selected_image_layers
                ),
                "status": shell.status,
            }
        )
    return rows


def _expected_effective_far_correction(case: str, periodic_model: str) -> str:
    if periodic_model == "infinite_physical" or case == "new_infinite_physical":
        return "cached_kneq0"
    if periodic_model == "configured":
        return "none"
    raise ValidationError(f"unsupported requested periodic model: {periodic_model}")


def _periodic_cache_provenance(
    metadata: Mapping[str, Any],
    *,
    effective_far_correction: str,
    validation_root: Path,
    verified_prime: Mapping[str, Any] | None,
    file_hashes: dict[Path, str],
) -> dict[str, Any] | None:
    if effective_far_correction != "cached_kneq0":
        return None
    raw = metadata.get("periodic_cache")
    if not isinstance(raw, Mapping):
        raise ValidationError("cached evaluator did not report periodic cache metadata")
    hit = raw.get("hit")
    try:
        build_count = int(raw.get("build_count"))
    except (TypeError, ValueError) as exc:
        raise ValidationError("periodic cache build count is invalid") from exc
    fingerprint = str(raw.get("fingerprint") or "")
    cache_path = _require_descendant_path(
        validation_root,
        str(raw.get("path") or ""),
        prefix="cache/periodic2",
        label="analysis periodic cache path",
    )
    if not bool(hit) or build_count != 0:
        raise ValidationError(
            "cached analysis evaluator must reuse the staged operator"
        )
    if not fingerprint:
        raise ValidationError("cached analysis evaluator has no cache fingerprint")
    if not cache_path.is_file() or cache_path.stat().st_size == 0:
        raise ValidationError(f"analysis periodic cache is missing: {cache_path}")
    path_sha256 = file_hashes.get(cache_path)
    if path_sha256 is None:
        path_sha256 = _sha256(cache_path)
        file_hashes[cache_path] = path_sha256
    prime_receipt_sha256 = ""
    if verified_prime is not None:
        if (
            fingerprint != verified_prime.get("fingerprint")
            or str(cache_path) != verified_prime.get("path")
            or path_sha256 != verified_prime.get("path_sha256")
        ):
            raise ValidationError(
                "analysis periodic cache differs from the verified cache prime"
            )
        prime_receipt_sha256 = str(verified_prime["receipt_sha256"])
    return {
        "hit": True,
        "build_count": build_count,
        "fingerprint": fingerprint,
        "path": str(cache_path),
        "path_sha256": path_sha256,
        "prime_receipt_sha256": prime_receipt_sha256,
    }


def _verified_cache_prime_contract(validation_root: Path) -> dict[str, Any]:
    receipt_path = _require_expected_path(
        validation_root,
        validation_root / "provenance/verified/cache_prime.json",
        "provenance/verified/cache_prime.json",
        label="verified cache-prime receipt",
    )
    receipt = _load_execution_receipt(receipt_path)
    cache = receipt.get("cache")
    if not isinstance(cache, Mapping):
        raise ValidationError("verified cache-prime receipt has invalid cache metadata")
    fingerprint = str(cache.get("fingerprint") or "")
    path = str(cache.get("path") or "")
    path_sha256 = str(cache.get("path_sha256") or "")
    if not fingerprint or not path or len(path_sha256) != 64:
        raise ValidationError("verified cache-prime identity is incomplete")
    cache_path = _require_descendant_path(
        validation_root,
        path,
        prefix="cache/periodic2",
        label="verified cache-prime path",
    )
    return {
        "fingerprint": fingerprint,
        "path": str(cache_path),
        "path_sha256": path_sha256,
        "receipt_sha256": _sha256(receipt_path),
    }


def _evaluate_object_physics_at_step(
    *,
    archive_run: Path,
    validation_root: Path,
    library: Path,
    run_rows: Sequence[Mapping[str, Any]],
    runs: Mapping[str, Beach | None],
    step: int | None,
    evaluation_target_ids: Sequence[int],
    run_shell: bool,
    allow_geometry_radius_fallback: bool = True,
    require_complete_contract: bool = False,
    verified_cache_prime: Mapping[str, Any] | None = None,
    cache_file_hashes: dict[Path, str] | None = None,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, Any],
]:
    wrench_rows: list[dict[str, Any]] = []
    curve_rows: list[dict[str, Any]] = []
    path_rows: list[dict[str, Any]] = []
    shell_rows: list[dict[str, Any]] = []
    if cache_file_hashes is None:
        cache_file_hashes = {}
    successes = 0
    failures: list[dict[str, str]] = []
    effective_models: dict[tuple[str, str], str] = {}
    target_ids = (6, 7)
    models = ("configured", "infinite_physical")

    for row in run_rows:
        case = str(row["case"])
        run = runs[case]
        if run is None:
            for model in models:
                wrench_rows.append(
                    {
                        "case": case,
                        "periodic_model": model,
                        "status": row["status"],
                        "message": row.get("message", ""),
                    }
                )
                path_rows.append(
                    {
                        "case": case,
                        "periodic_model": model,
                        "path_status": "not_evaluated",
                        "speed_status": "not_evaluated",
                        "status": row["status"],
                        "message": row.get("message", ""),
                    }
                )
            continue
        result = run.result
        available_ids = (
            set(int(value) for value in np.asarray(result.mesh_ids).tolist())
            if result.mesh_ids is not None
            else set()
        )
        missing_ids = sorted(set(target_ids) - available_ids)
        if missing_ids:
            message = "required target mesh ids 6/7 are unavailable: " + ", ".join(
                str(value) for value in missing_ids
            )
            failures.append({"case": case, "message": message})
            wrench_rows.append(
                {
                    "case": case,
                    "periodic_model": "configured",
                    "status": "not_evaluated",
                    "message": message,
                }
            )
            path_rows.append(
                {
                    "case": case,
                    "periodic_model": "configured",
                    "path_status": "not_evaluated",
                    "speed_status": "not_evaluated",
                    "status": "not_evaluated",
                    "message": message,
                }
            )
            continue
        selected_ids = tuple(int(value) for value in evaluation_target_ids)
        case_models = ("configured",) if case == "new_infinite_physical" else models
        for periodic_model in case_models:
            try:
                snapshot_context = run.object_interaction_snapshot(
                    step=step,
                    periodic_model=periodic_model,
                    cache_dir=validation_root / "cache/periodic2",
                    generation_tolerance=DEFAULT_GENERATION_TOLERANCE,
                    library_path=library,
                )
                with snapshot_context as snapshot:
                    for mesh_id in selected_ids:
                        try:
                            probe = snapshot.object_probe(mesh_id)
                            geometry_radius, geometry_representation = (
                                _probe_geometry_contract(probe)
                            )
                            mechanics = _release_parameters(
                                archive_run,
                                geometry_radius,
                                allow_geometry_radius_fallback=(
                                    allow_geometry_radius_fallback
                                ),
                            )
                            radius = mechanics["radius_m"]
                            initial_grid = np.linspace(0.0, 2.0 * radius, 17)
                            wrench = probe.wrench(components=True)
                            metadata = wrench.numerical_metadata
                            if not isinstance(metadata, Mapping):
                                raise ValidationError(
                                    "wrench numerical metadata is unavailable"
                                )
                            _validate_geometry_metadata(
                                metadata,
                                probe_radius_m=geometry_radius,
                            )
                            torque_origin_policy, torque_origin_m = (
                                _torque_origin_provenance(metadata)
                            )
                            _validate_wrench_torque_origin(wrench, torque_origin_m)
                            effective = str(
                                metadata.get("effective_far_correction", "unknown")
                            )
                            expected_effective = _expected_effective_far_correction(
                                case, periodic_model
                            )
                            if effective != expected_effective:
                                raise ValidationError(
                                    "effective far correction mismatch for "
                                    f"{case}/{periodic_model}: expected "
                                    f"{expected_effective}, got {effective}"
                                )
                            previous_effective = effective_models.get(
                                (case, periodic_model)
                            )
                            if (
                                previous_effective is not None
                                and previous_effective != effective
                            ):
                                raise ValidationError(
                                    "effective far correction changed across target meshes"
                                )
                            effective_models[(case, periodic_model)] = effective
                            cache_provenance = _periodic_cache_provenance(
                                metadata,
                                effective_far_correction=effective,
                                validation_root=validation_root,
                                verified_prime=verified_cache_prime,
                                file_hashes=cache_file_hashes,
                            )
                            physical_components = dict(wrench.components)
                            if require_complete_contract:
                                _validate_strict_wrench_component_contract(
                                    physical_components,
                                    metadata,
                                    effective_far_correction=effective,
                                    wrench_force_N=wrench.force_N,
                                    wrench_torque_Nm=wrench.torque_Nm,
                                )
                            elif "total_external" not in physical_components:
                                physical_components["total_external"] = SimpleNamespace(
                                    force_N=wrench.force_N,
                                    torque_Nm=wrench.torque_Nm,
                                    potential_energy_J=None,
                                )
                            for (
                                component_name,
                                component,
                            ) in physical_components.items():
                                wrench_rows.append(
                                    _wrench_row(
                                        case=case,
                                        periodic_model=periodic_model,
                                        effective_far_correction="unknown",
                                        zero_mode_policy="unknown",
                                        lower_boundary_model="unknown",
                                        periodic_cache=cache_provenance,
                                        step_selector=(
                                            "final" if step is None else step
                                        ),
                                        resolved_batch=(
                                            int(result.batches)
                                            if step is None
                                            else step
                                        ),
                                        mesh_id=mesh_id,
                                        mechanics=mechanics,
                                        target_geometry_representation=(
                                            geometry_representation
                                        ),
                                        component=component_name,
                                        component_kind=(
                                            "physical_total"
                                            if component_name == "total_external"
                                            else "physical_additive"
                                        ),
                                        force=component.force_N,
                                        torque=component.torque_Nm,
                                        torque_origin_policy=torque_origin_policy,
                                        torque_origin_m=torque_origin_m,
                                        total_charge=float(wrench.total_charge_C),
                                        potential_energy=getattr(
                                            component, "potential_energy_J", None
                                        ),
                                    )
                                )
                            for component_name in (
                                "periodic_kneq0",
                                "physical_k0",
                                "primary_free_subtraction",
                                "cached_kneq0_trace_correction",
                            ):
                                component = metadata.get(component_name)
                                if not isinstance(component, Mapping):
                                    continue
                                wrench_rows.append(
                                    _wrench_row(
                                        case=case,
                                        periodic_model=periodic_model,
                                        effective_far_correction="unknown",
                                        zero_mode_policy="unknown",
                                        lower_boundary_model="unknown",
                                        periodic_cache=cache_provenance,
                                        step_selector=(
                                            "final" if step is None else step
                                        ),
                                        resolved_batch=(
                                            int(result.batches)
                                            if step is None
                                            else step
                                        ),
                                        mesh_id=mesh_id,
                                        mechanics=mechanics,
                                        target_geometry_representation=(
                                            geometry_representation
                                        ),
                                        component=f"numerical:{component_name}",
                                        component_kind=(
                                            "numerical_diagnostic_included"
                                            if component_name
                                            == "cached_kneq0_trace_correction"
                                            else "numerical_decomposition"
                                        ),
                                        force=component["force_N"],
                                        torque=component["torque_Nm"],
                                        torque_origin_policy=torque_origin_policy,
                                        torque_origin_m=torque_origin_m,
                                        total_charge=float(wrench.total_charge_C),
                                        potential_energy=component.get(
                                            "potential_energy_J"
                                        ),
                                    )
                                )

                            path = probe.vertical_path(
                                initial_grid,
                                adaptive=True,
                                relative_tolerance=5.0e-3,
                                force_absolute_tolerance_N=1.0e-12,
                                work_absolute_tolerance_J=1.0e-18,
                                max_refinement=8,
                                components=True,
                            )
                            path_metadata = path.numerical_metadata
                            if not isinstance(path_metadata, Mapping):
                                raise ValidationError(
                                    "path numerical metadata is unavailable"
                                )
                            _validate_geometry_metadata(
                                path_metadata,
                                probe_radius_m=geometry_radius,
                            )
                            path_torque_policy, path_torque_origins = (
                                _torque_origin_provenance(
                                    path_metadata,
                                    point_count=len(path.displacement_m),
                                )
                            )
                            _validate_production_torque_origin_contract(
                                probe=probe,
                                wrench_policy=torque_origin_policy,
                                wrench_origin_m=torque_origin_m,
                                path_policy=path_torque_policy,
                                path_origin_m=path_torque_origins,
                                displacement_m=path.displacement_m,
                            )
                            if require_complete_contract:
                                _validate_strict_path_component_contract(
                                    path,
                                    expected_start_m=0.0,
                                    expected_end_m=2.0 * radius,
                                )
                            release = path.evaluate_release(
                                mass_kg=mechanics["mass_kg"],
                                gravity_m_s2=mechanics["gravity_m_s2"],
                                adhesion=mechanics["adhesion"],
                                eta_translation=mechanics["eta_translation"],
                                energy_tolerance_J=mechanics["energy_tolerance_J"],
                            )
                            release_sensitivity = path.evaluate_release(
                                mass_kg=mechanics["mass_kg"],
                                gravity_m_s2=mechanics["gravity_m_s2"],
                                adhesion=mechanics["adhesion"],
                                eta_translation=mechanics[
                                    "eta_translation_sensitivity"
                                ],
                                energy_tolerance_J=mechanics["energy_tolerance_J"],
                            )
                            potential_work = path.potential_difference_work_J
                            potential_work_available = potential_work is not None
                            peak_force = float(
                                np.max(np.linalg.norm(path.force_N, axis=1))
                            )
                            force_floor_to_peak = (
                                1.0e-12 / peak_force if peak_force > 0.0 else math.inf
                            )
                            force_margin_over_tolerance = (
                                abs(release.instantaneous_force_margin_N) / 1.0e-12
                            )
                            force_resolution_qualified = (
                                force_floor_to_peak <= 5.0e-3
                                and force_margin_over_tolerance > 1.0
                            )
                            barrier_decision_margin = (
                                release.minimum_available_energy_J
                                + mechanics["energy_tolerance_J"]
                            )
                            barrier_decision_ratio = (
                                abs(barrier_decision_margin)
                                / mechanics["energy_tolerance_J"]
                                if mechanics["energy_tolerance_J"] > 0.0
                                else math.inf
                            )
                            barrier_resolution_qualified = (
                                barrier_decision_ratio > 1.0e-1
                            )
                            numerically_qualified = (
                                path.status == "converged"
                                and potential_work_available
                                and force_resolution_qualified
                                and barrier_resolution_qualified
                            )
                            for index, displacement in enumerate(path.displacement_m):
                                curve_rows.append(
                                    {
                                        "case": case,
                                        "periodic_model": periodic_model,
                                        "mesh_id": mesh_id,
                                        **_radius_provenance_fields(
                                            mechanics, geometry_representation
                                        ),
                                        "component": "total_external",
                                        "displacement_m": displacement,
                                        "force_x_N": path.force_N[index, 0],
                                        "force_y_N": path.force_N[index, 1],
                                        "force_z_N": path.force_N[index, 2],
                                        "torque_x_Nm": path.torque_Nm[index, 0],
                                        "torque_y_Nm": path.torque_Nm[index, 1],
                                        "torque_z_Nm": path.torque_Nm[index, 2],
                                        "torque_origin_policy": path_torque_policy,
                                        "torque_origin_x_m": path_torque_origins[
                                            index, 0
                                        ],
                                        "torque_origin_y_m": path_torque_origins[
                                            index, 1
                                        ],
                                        "torque_origin_z_m": path_torque_origins[
                                            index, 2
                                        ],
                                        "electrostatic_work_J": path.electrostatic_work_J[
                                            index
                                        ],
                                        "potential_difference_work_J": (
                                            ""
                                            if potential_work is None
                                            else potential_work[index]
                                        ),
                                        "gravity_work_J": release.gravity_work_J[index],
                                        "adhesion_work_J": release.adhesion_work_J[
                                            index
                                        ],
                                        "available_energy_J": release.available_energy_J[
                                            index
                                        ],
                                        "electrostatic_only_speed_m_s": release.electrostatic_only_speed_m_s[
                                            index
                                        ],
                                        "gravity_corrected_speed_m_s": release.gravity_corrected_speed_m_s[
                                            index
                                        ],
                                        "speed_m_s": release.speed_m_s[index],
                                        "eta_translation": mechanics["eta_translation"],
                                        "status": path.status,
                                    }
                                )
                                for (
                                    component_name,
                                    values,
                                ) in path.component_force_N.items():
                                    if component_name == "total_external":
                                        continue
                                    torques = path.component_torque_Nm[component_name]
                                    curve_rows.append(
                                        {
                                            "case": case,
                                            "periodic_model": periodic_model,
                                            "mesh_id": mesh_id,
                                            **_radius_provenance_fields(
                                                mechanics, geometry_representation
                                            ),
                                            "component": component_name,
                                            "displacement_m": displacement,
                                            "force_x_N": values[index, 0],
                                            "force_y_N": values[index, 1],
                                            "force_z_N": values[index, 2],
                                            "torque_x_Nm": torques[index, 0],
                                            "torque_y_Nm": torques[index, 1],
                                            "torque_z_Nm": torques[index, 2],
                                            "torque_origin_policy": (
                                                path_torque_policy
                                            ),
                                            "torque_origin_x_m": (
                                                path_torque_origins[index, 0]
                                            ),
                                            "torque_origin_y_m": (
                                                path_torque_origins[index, 1]
                                            ),
                                            "torque_origin_z_m": (
                                                path_torque_origins[index, 2]
                                            ),
                                            "status": path.status,
                                        }
                                    )
                            path_rows.append(
                                {
                                    "case": case,
                                    "periodic_model": periodic_model,
                                    "mesh_id": mesh_id,
                                    "radius_m": radius,
                                    "radius_source": mechanics["radius_source"],
                                    "geometry_radius_m": mechanics["geometry_radius_m"],
                                    "radius_relative_mismatch": mechanics[
                                        "radius_relative_mismatch"
                                    ],
                                    "radius_relative_tolerance": mechanics[
                                        "radius_relative_tolerance"
                                    ],
                                    "target_geometry_representation": (
                                        geometry_representation
                                    ),
                                    "torque_origin_policy": path_torque_policy,
                                    "initial_torque_origin_x_m": (
                                        path_torque_origins[0, 0]
                                    ),
                                    "initial_torque_origin_y_m": (
                                        path_torque_origins[0, 1]
                                    ),
                                    "initial_torque_origin_z_m": (
                                        path_torque_origins[0, 2]
                                    ),
                                    "mass_kg": mechanics["mass_kg"],
                                    "gravity_m_s2": mechanics["gravity_m_s2"],
                                    "adhesion_model": mechanics["adhesion_model"],
                                    "adhesion_force_N": mechanics["adhesion_force_N"],
                                    "adhesion_work_J": mechanics["adhesion_work_J"],
                                    "adhesion_equivalent_range_m": mechanics[
                                        "adhesion_equivalent_range_m"
                                    ],
                                    "eta_translation": mechanics["eta_translation"],
                                    "eta_translation_sensitivity": mechanics[
                                        "eta_translation_sensitivity"
                                    ],
                                    "energy_tolerance_J": mechanics[
                                        "energy_tolerance_J"
                                    ],
                                    "path_start_m": path.displacement_m[0],
                                    "path_end_m": path.displacement_m[-1],
                                    "path_status": path.status,
                                    "potential_work_available": potential_work_available,
                                    "numerically_qualified": numerically_qualified,
                                    "endpoint_work_J": path.electrostatic_work_J[-1],
                                    "max_force_z_N": float(np.max(path.force_N[:, 2])),
                                    "work_relative_mismatch": path.work_relative_mismatch,
                                    "path_relative_tolerance": 5.0e-3,
                                    "path_force_absolute_tolerance_N": 1.0e-12,
                                    "path_work_absolute_tolerance_J": 1.0e-18,
                                    "path_max_refinement": 8,
                                    "peak_force_N": peak_force,
                                    "force_floor_to_peak_ratio": force_floor_to_peak,
                                    "force_margin_over_absolute_tolerance": (
                                        force_margin_over_tolerance
                                    ),
                                    "force_resolution_qualified": (
                                        force_resolution_qualified
                                    ),
                                    "speed_status": (
                                        "numerically_unqualified"
                                        if not numerically_qualified
                                        else "conditional_local_release_proxy"
                                        if release.endpoint_reachable_from_rest
                                        else "mechanically_inaccessible"
                                    ),
                                    "barrier_free_from_rest": release.barrier_free_from_rest,
                                    "minimum_available_energy_J": (
                                        release.minimum_available_energy_J
                                    ),
                                    "minimum_energy_margin_over_tolerance": (
                                        release.minimum_available_energy_J
                                        / mechanics["energy_tolerance_J"]
                                        if mechanics["energy_tolerance_J"] > 0.0
                                        else ""
                                    ),
                                    "barrier_decision_margin_J": (
                                        barrier_decision_margin
                                    ),
                                    "barrier_decision_margin_over_tolerance": (
                                        barrier_decision_ratio
                                    ),
                                    "barrier_resolution_qualified": (
                                        barrier_resolution_qualified
                                    ),
                                    "first_inaccessible_displacement_m": (
                                        ""
                                        if release.first_inaccessible_displacement_m
                                        is None
                                        else release.first_inaccessible_displacement_m
                                    ),
                                    "endpoint_available_energy_J": (
                                        release.endpoint_available_energy_J
                                    ),
                                    "endpoint_positive": release.endpoint_positive,
                                    "endpoint_reachable_from_rest": release.endpoint_reachable_from_rest,
                                    "endpoint_speed_m_s": release.endpoint_speed_m_s,
                                    "endpoint_speed_sensitivity_m_s": (
                                        release_sensitivity.endpoint_speed_m_s
                                    ),
                                    "maximum_reachable_speed_m_s": release.maximum_reachable_speed_m_s,
                                    "instantaneous_force_margin_N": release.instantaneous_force_margin_N,
                                    "barrier_interpretation": (
                                        "profile_dependent_conservative_sufficient_if_true"
                                        if mechanics["adhesion_model"] != "none"
                                        else "no_adhesion_profile"
                                    ),
                                    "status": "available",
                                    "message": "",
                                }
                            )
                            if periodic_model == "infinite_physical" and run_shell:
                                shell = finite_shell_convergence(
                                    snapshot,
                                    probe,
                                    initial_grid,
                                    max_layers=12,
                                    relative_tolerance=1.0e-2,
                                    force_floor_N=1.0e-12,
                                    work_floor_J=1.0e-18,
                                )
                                shell_rows.extend(
                                    _finite_shell_rows(
                                        case=case,
                                        periodic_model=periodic_model,
                                        effective_far_correction="unknown",
                                        zero_mode_policy="unknown",
                                        lower_boundary_model="unknown",
                                        step_selector=(
                                            "final" if step is None else step
                                        ),
                                        resolved_batch=(
                                            int(result.batches)
                                            if step is None
                                            else step
                                        ),
                                        mesh_id=mesh_id,
                                        shell=shell,
                                        mechanics=mechanics,
                                        target_geometry_representation=(
                                            geometry_representation
                                        ),
                                        displacement_m=initial_grid,
                                    )
                                )
                            successes += 1
                        except Exception as exc:
                            message = str(exc)
                            failures.append(
                                {
                                    "case": case,
                                    "periodic_model": periodic_model,
                                    "mesh_id": str(mesh_id),
                                    "message": message,
                                }
                            )
                            wrench_rows.append(
                                {
                                    "case": case,
                                    "periodic_model": periodic_model,
                                    "mesh_id": mesh_id,
                                    "status": "invalid",
                                    "message": message,
                                }
                            )
                            path_rows.append(
                                {
                                    "case": case,
                                    "periodic_model": periodic_model,
                                    "mesh_id": mesh_id,
                                    "path_status": "not_evaluated",
                                    "speed_status": "not_evaluated",
                                    "status": "invalid",
                                    "message": message,
                                }
                            )
            except Exception as exc:
                message = str(exc)
                failures.append(
                    {
                        "case": case,
                        "periodic_model": periodic_model,
                        "message": message,
                    }
                )
                for mesh_id in selected_ids:
                    wrench_rows.append(
                        {
                            "case": case,
                            "periodic_model": periodic_model,
                            "mesh_id": mesh_id,
                            "status": "invalid",
                            "message": message,
                        }
                    )
                    path_rows.append(
                        {
                            "case": case,
                            "periodic_model": periodic_model,
                            "mesh_id": mesh_id,
                            "path_status": "not_evaluated",
                            "speed_status": "not_evaluated",
                            "status": "invalid",
                            "message": message,
                        }
                    )
    step_selector: int | str = "final" if step is None else step
    batches_by_case = {
        case: int(getattr(run.result, "batches", 0))
        for case, run in runs.items()
        if run is not None
    }
    for rows in (wrench_rows, curve_rows, path_rows, shell_rows):
        for item in rows:
            case = str(item.get("case", ""))
            model = str(item.get("periodic_model", ""))
            effective = effective_models.get((case, model), "unknown")
            item["step_selector"] = step_selector
            item["resolved_batch"] = (
                batches_by_case.get(case, 0) if step is None else step
            )
            item["effective_far_correction"] = effective
            item["zero_mode_policy"] = (
                "exclude_k0"
                if effective == "cached_kneq0"
                else "legacy_not_decomposed"
                if effective == "none"
                else "unknown"
            )
            item["lower_boundary_model"] = (
                "e_bottom_zero"
                if effective == "cached_kneq0"
                else "legacy_implicit"
                if effective == "none"
                else "not_applicable"
                if effective == "free"
                else "unknown"
            )
    status = "available" if successes > 0 and not failures else "partial"
    if successes == 0:
        status = "not_evaluated"
    return (
        wrench_rows,
        curve_rows,
        path_rows,
        shell_rows,
        {"status": status, "successful_target_models": successes, "failures": failures},
    )


def _evaluate_new_infinite_shell_reference(
    *,
    archive_run: Path,
    validation_root: Path,
    library: Path,
    run: Beach,
    allow_geometry_radius_fallback: bool = True,
    verified_cache_prime: Mapping[str, Any] | None = None,
    cache_file_hashes: dict[Path, str] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """Evaluate finite-shell convergence for the new infinite run's final charge."""

    result = run.result
    if cache_file_hashes is None:
        cache_file_hashes = {}
    target_ids = (6, 7)
    available_ids = (
        set(int(value) for value in np.asarray(result.mesh_ids).tolist())
        if result.mesh_ids is not None
        else set()
    )
    missing_ids = sorted(set(target_ids) - available_ids)
    if missing_ids:
        return [], [
            {
                "case": "new_infinite_physical",
                "periodic_model": "infinite_physical",
                "message": "required target mesh ids 6/7 are unavailable: "
                + ", ".join(str(value) for value in missing_ids),
            }
        ]

    rows: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    try:
        snapshot_context = run.object_interaction_snapshot(
            step=None,
            periodic_model="infinite_physical",
            cache_dir=validation_root / "cache/periodic2",
            generation_tolerance=DEFAULT_GENERATION_TOLERANCE,
            library_path=library,
        )
        with snapshot_context as snapshot:
            for mesh_id in target_ids:
                try:
                    probe = snapshot.object_probe(mesh_id)
                    geometry_radius, _geometry_representation = (
                        _probe_geometry_contract(probe)
                    )
                    mechanics = _release_parameters(
                        archive_run,
                        geometry_radius,
                        allow_geometry_radius_fallback=(allow_geometry_radius_fallback),
                    )
                    if verified_cache_prime is not None:
                        wrench = probe.wrench(components=False)
                        metadata = wrench.numerical_metadata
                        if not isinstance(metadata, Mapping):
                            raise ValidationError(
                                "shell reference cache metadata is unavailable"
                            )
                        _validate_geometry_metadata(
                            metadata,
                            probe_radius_m=geometry_radius,
                        )
                        effective = str(
                            metadata.get("effective_far_correction", "unknown")
                        )
                        if effective != "cached_kneq0":
                            raise ValidationError(
                                "shell reference did not use cached_kneq0"
                            )
                        _periodic_cache_provenance(
                            metadata,
                            effective_far_correction=effective,
                            validation_root=validation_root,
                            verified_prime=verified_cache_prime,
                            file_hashes=cache_file_hashes,
                        )
                    radius = mechanics["radius_m"]
                    initial_grid = np.linspace(0.0, 2.0 * radius, 17)
                    shell = finite_shell_convergence(
                        snapshot,
                        probe,
                        initial_grid,
                        max_layers=12,
                        relative_tolerance=1.0e-2,
                        force_floor_N=1.0e-12,
                        work_floor_J=1.0e-18,
                    )
                    rows.extend(
                        _finite_shell_rows(
                            case="new_infinite_physical",
                            periodic_model="infinite_physical",
                            effective_far_correction="cached_kneq0",
                            zero_mode_policy="exclude_k0",
                            lower_boundary_model="e_bottom_zero",
                            step_selector="final",
                            resolved_batch=int(result.batches),
                            mesh_id=mesh_id,
                            shell=shell,
                            mechanics=mechanics,
                            target_geometry_representation=(_geometry_representation),
                            displacement_m=initial_grid,
                        )
                    )
                except Exception as exc:
                    failures.append(
                        {
                            "case": "new_infinite_physical",
                            "periodic_model": "infinite_physical",
                            "mesh_id": str(mesh_id),
                            "message": str(exc),
                        }
                    )
    except Exception as exc:
        failures.append(
            {
                "case": "new_infinite_physical",
                "periodic_model": "infinite_physical",
                "message": str(exc),
            }
        )
    return rows, failures


def _evaluate_object_physics(
    *,
    archive_run: Path,
    validation_root: Path,
    library: Path,
    run_rows: Sequence[Mapping[str, Any]],
    runs: Mapping[str, Beach | None],
    require_history: bool = False,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, Any],
]:
    plans: tuple[tuple[int | None, tuple[int, ...]], ...] = (
        (149001, (7,)),
        (180001, (6,)),
        (279001, (6, 7)),
        (None, (6, 7)),
    )
    wrench_rows: list[dict[str, Any]] = []
    curve_rows: list[dict[str, Any]] = []
    path_rows: list[dict[str, Any]] = []
    shell_rows: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    successes = 0
    evaluated_snapshots: list[dict[str, Any]] = []
    verified_cache_prime = (
        _verified_cache_prime_contract(validation_root) if require_history else None
    )
    cache_file_hashes: dict[Path, str] = {}

    for step, target_ids in plans:
        eligible_rows: list[Mapping[str, Any]] = []
        eligible_runs: dict[str, Beach | None] = {}
        for row in run_rows:
            case = str(row["case"])
            run = runs[case]
            if step is None or run is None:
                eligible_rows.append(row)
                eligible_runs[case] = run
                continue
            history = getattr(run.result, "history", None)
            if history is None or not bool(history.has_data):
                if require_history:
                    failures.append(
                        {
                            "case": case,
                            "step": str(step),
                            "message": "required charge history is unavailable",
                        }
                    )
                continue
            available = set(
                int(value)
                for value in np.asarray(history.batch_indices, dtype=np.int64)
            )
            if step not in available:
                failures.append(
                    {
                        "case": case,
                        "step": str(step),
                        "message": f"required charge-history batch {step} is unavailable",
                    }
                )
                continue
            eligible_rows.append(row)
            eligible_runs[case] = run
        if not eligible_rows:
            continue
        evaluated = _evaluate_object_physics_at_step(
            archive_run=archive_run,
            validation_root=validation_root,
            library=library,
            run_rows=eligible_rows,
            runs=eligible_runs,
            step=step,
            evaluation_target_ids=target_ids,
            run_shell=step is None,
            allow_geometry_radius_fallback=not require_history,
            require_complete_contract=require_history,
            verified_cache_prime=verified_cache_prime,
            cache_file_hashes=cache_file_hashes,
        )
        step_wrench, step_curves, step_paths, step_shells, step_status = evaluated
        wrench_rows.extend(step_wrench)
        curve_rows.extend(step_curves)
        path_rows.extend(step_paths)
        shell_rows.extend(step_shells)
        step_failures = step_status.get("failures", [])
        if isinstance(step_failures, list):
            failures.extend(
                dict(value) for value in step_failures if isinstance(value, Mapping)
            )
        step_successes = int(step_status.get("successful_target_models", 0))
        successes += step_successes
        evaluated_snapshots.append(
            {
                "step_selector": "final" if step is None else step,
                "target_mesh_ids": list(target_ids),
                "successful_target_models": step_successes,
            }
        )

    new_infinite_run = runs.get("new_infinite_physical")
    if new_infinite_run is not None:
        reference_rows, reference_failures = _evaluate_new_infinite_shell_reference(
            archive_run=archive_run,
            validation_root=validation_root,
            library=library,
            run=new_infinite_run,
            allow_geometry_radius_fallback=not require_history,
            verified_cache_prime=verified_cache_prime,
            cache_file_hashes=cache_file_hashes,
        )
        shell_rows.extend(reference_rows)
        failures.extend(reference_failures)
        evaluated_snapshots.append(
            {
                "step_selector": "final",
                "target_mesh_ids": [6, 7],
                "periodic_model": "infinite_physical",
                "diagnostic": "finite_shell_reference_for_new_infinite_charge",
                "successful_target_models": len(
                    {int(row["mesh_id"]) for row in reference_rows}
                ),
            }
        )

    status = "available" if successes > 0 and not failures else "partial"
    if successes == 0:
        status = "not_evaluated"
    return (
        wrench_rows,
        curve_rows,
        path_rows,
        shell_rows,
        {
            "status": status,
            "successful_target_models": successes,
            "failures": failures,
            "required_history_batches": [149001, 180001, 279001],
            "evaluated_snapshots": evaluated_snapshots,
        },
    )


def _physics_review_lines(
    wrench_rows: Sequence[Mapping[str, Any]],
    path_rows: Sequence[Mapping[str, Any]],
    shell_rows: Sequence[Mapping[str, Any]],
) -> str:
    available_paths = [row for row in path_rows if row.get("status") == "available"]
    if not available_paths:
        return (
            "- fixed-discretization path/work/shell gate: "
            "object force/path/shell の実評価結果はまだありません。"
        )
    wrench_index = {
        (
            str(row.get("case")),
            str(row.get("periodic_model")),
            int(row.get("resolved_batch", 0)),
            str(row.get("mesh_id")),
            str(row.get("component")),
        ): row
        for row in wrench_rows
        if row.get("status") == "available"
    }
    shell_status = {
        (
            str(row.get("case")),
            str(row.get("periodic_model")),
            int(row.get("resolved_batch", 0)),
            str(row.get("mesh_id")),
        ): str(row.get("status"))
        for row in shell_rows
    }
    lines: list[str] = []
    for path in available_paths:
        case = str(path["case"])
        model = str(path["periodic_model"])
        batch = int(path["resolved_batch"])
        mesh = str(path["mesh_id"])
        total = wrench_index.get((case, model, batch, mesh, "total_external"), {})
        replicas = wrench_index.get(
            (case, model, batch, mesh, "target_periodic_images"), {}
        )
        zero_mode = wrench_index.get(
            (case, model, batch, mesh, "numerical:physical_k0"), {}
        )
        lines.append(
            "- "
            f"{case} / {model} / batch {batch} / mesh {mesh}: "
            f"Fz={total.get('force_z_N', 'unavailable')} N; "
            f"target replicas Fz={replicas.get('force_z_N', 'unavailable')} N; "
            f"physical_k0 Fz={zero_mode.get('force_z_N', 'unavailable')} N; "
            f"endpoint work={path.get('endpoint_work_J')} J; "
            f"potential work available={path.get('potential_work_available')}; "
            "fixed-discretization path/work/shell gate: path_work="
            f"{path.get('numerically_qualified')}; "
            "barrier on fixed saved discretization="
            f"{path.get('barrier_free_from_rest') if path.get('numerically_qualified') else 'not_claimed_unqualified'}; "
            f"instantaneous margin={path.get('instantaneous_force_margin_N')} N; "
            "shell="
            f"{shell_status.get((case, model, batch, mesh), 'not_evaluated')}"
        )
    return "\n".join(lines)


def _expected_object_evaluation_keys() -> set[tuple[str, str, int, int]]:
    states = (
        (149001, 7),
        (180001, 6),
        (279001, 6),
        (279001, 7),
        (280000, 6),
        (280000, 7),
    )
    keys: set[tuple[str, str, int, int]] = set()
    for case in (
        "archived_v1_3",
        "new_finite_configured",
        "new_infinite_physical",
    ):
        models = (
            ("configured",)
            if case == "new_infinite_physical"
            else ("configured", "infinite_physical")
        )
        for model in models:
            keys.update((case, model, batch, mesh_id) for batch, mesh_id in states)
    return keys


def _local_model_numerical_qualification(
    physics_evaluation: Mapping[str, Any],
    wrench_rows: Sequence[Mapping[str, Any]],
    path_rows: Sequence[Mapping[str, Any]],
    shell_rows: Sequence[Mapping[str, Any]],
    case_status: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    evaluated_paths = [row for row in path_rows if row.get("status") == "available"]
    expected_path_keys = _expected_object_evaluation_keys()
    path_keys = [
        (
            str(row.get("case")),
            str(row.get("periodic_model")),
            int(row.get("resolved_batch", 0)),
            int(row.get("mesh_id", 0)),
        )
        for row in evaluated_paths
    ]
    path_counts = Counter(path_keys)
    actual_path_keys = set(path_keys)
    missing_path_keys = sorted(expected_path_keys - actual_path_keys)
    unexpected_path_keys = sorted(actual_path_keys - expected_path_keys)
    duplicate_path_keys = sorted(
        key for key, count in path_counts.items() if count != 1
    )
    total_wrenches = [
        row
        for row in wrench_rows
        if row.get("status") == "available" and row.get("component") == "total_external"
    ]
    wrench_keys = [
        (
            str(row.get("case")),
            str(row.get("periodic_model")),
            int(row.get("resolved_batch", 0)),
            int(row.get("mesh_id", 0)),
        )
        for row in total_wrenches
    ]
    wrench_counts = Counter(wrench_keys)
    actual_wrench_keys = set(wrench_keys)
    missing_wrench_keys = sorted(expected_path_keys - actual_wrench_keys)
    unexpected_wrench_keys = sorted(actual_wrench_keys - expected_path_keys)
    duplicate_wrench_keys = sorted(
        key for key, count in wrench_counts.items() if count != 1
    )
    unqualified_paths = [
        {
            "case": row.get("case"),
            "periodic_model": row.get("periodic_model"),
            "resolved_batch": row.get("resolved_batch"),
            "mesh_id": row.get("mesh_id"),
            "path_status": row.get("path_status"),
        }
        for row in evaluated_paths
        if not bool(row.get("numerically_qualified"))
    ]
    expected_cases = {
        "archived_v1_3",
        "new_finite_configured",
        "new_infinite_physical",
    }
    missing_cases = sorted(
        case
        for case in expected_cases
        if case_status.get(case, {}).get("status") != "available"
    )
    expected_path_count = len(expected_path_keys)
    missing_path_count = max(0, expected_path_count - len(evaluated_paths))
    shell_groups = {
        (
            str(row.get("case")),
            str(row.get("periodic_model")),
            int(row.get("resolved_batch", 0)),
            int(row.get("mesh_id", 0)),
        ): str(row.get("status"))
        for row in shell_rows
        if str(row.get("mesh_id", "")).isdigit()
    }
    expected_shell_groups = {
        (case, "infinite_physical", 280000, mesh_id)
        for case in expected_cases
        for mesh_id in (6, 7)
    }
    missing_shell_groups = sorted(expected_shell_groups - set(shell_groups))
    shell_failures = [
        {
            "case": key[0],
            "periodic_model": key[1],
            "resolved_batch": key[2],
            "mesh_id": key[3],
            "status": value,
        }
        for key, value in shell_groups.items()
        if value != "converged"
    ]
    qualified = (
        physics_evaluation.get("status") == "available"
        and not missing_cases
        and len(evaluated_paths) == expected_path_count
        and not missing_path_keys
        and not unexpected_path_keys
        and not duplicate_path_keys
        and not missing_wrench_keys
        and not unexpected_wrench_keys
        and not duplicate_wrench_keys
        and not unqualified_paths
        and not missing_shell_groups
        and not shell_failures
    )
    return {
        "status": "qualified" if qualified else "not_qualified",
        "status_semantics": "path_work_shell_on_fixed_saved_discretization",
        "verified_numerical_axes": [
            "fixed_saved_discretization_coverage",
            "path_integration",
            "work_potential_consistency",
            "force_and_barrier_decision_resolution",
            "finite_shell_to_infinite_reference",
        ],
        "unverified_numerical_axes": [
            "saved_sphere_mesh_refinement",
            "source_discretization_refinement",
            "sphere_absolute_force_error",
            "sphere_absolute_torque_error",
        ],
        "saved_sphere_mesh_refinement_status": "not_evaluated",
        "source_discretization_refinement_status": "not_evaluated",
        "sphere_absolute_force_error_status": "not_evaluated",
        "sphere_absolute_torque_error_status": "not_evaluated",
        "plane_oracle_used_as_sphere_error_bar": False,
        "structurally_complete": physics_evaluation.get("status") == "available",
        "expected_path_count": expected_path_count,
        "evaluated_path_count": len(evaluated_paths),
        "missing_path_count": missing_path_count,
        "missing_path_keys": [list(value) for value in missing_path_keys],
        "unexpected_path_keys": [list(value) for value in unexpected_path_keys],
        "duplicate_path_keys": [list(value) for value in duplicate_path_keys],
        "missing_wrench_keys": [list(value) for value in missing_wrench_keys],
        "unexpected_wrench_keys": [list(value) for value in unexpected_wrench_keys],
        "duplicate_wrench_keys": [list(value) for value in duplicate_wrench_keys],
        "missing_cases": missing_cases,
        "unqualified_paths": unqualified_paths,
        "shell_group_count": len(shell_groups),
        "missing_shell_groups": [list(value) for value in missing_shell_groups],
        "shell_failures": shell_failures,
        "claim_scope": (
            "fixed_saved_mesh_and_source_discretization; "
            "local_frozen_field_0_to_2R_path_work_shell_only; "
            "not full_discretization_or_escape_to_infinity"
        ),
        "remaining_physical_conditions": [
            "frozen charge during motion",
            "E_bottom=0 closure selection",
            "single paired random seed",
            "equivalent finite-range adhesion profile",
            "nonlocal plasma closure not included in escape energy",
        ],
    }


def _verify_complete_runs(
    root: Path,
    manifest: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    reports: dict[str, dict[str, Any]] = {}
    cases = manifest.get("cases", {})
    if not isinstance(cases, dict):
        raise ValidationError("complete analysis requires staged case metadata")
    for name in STRICT_CASES:
        case = cases.get(name)
        if not isinstance(case, dict):
            raise ValidationError(f"complete analysis requires staged case {name}")
        try:
            report = verify_run(
                case["output_dir"],
                int(case["batch_count"]),
                require_existing_receipt=True,
            )
        except (KeyError, TypeError, ValueError, ValidationError) as exc:
            raise ValidationError(
                f"complete analysis requires verified case {name}: {exc}"
            ) from exc
        reports[name] = report
    _verify_pair_fingerprint_contract(reports)
    return reports


def _verify_pair_fingerprint_contract(
    reports: Mapping[str, Mapping[str, Any]],
) -> None:
    def values(case_names: Sequence[str], key: str) -> set[str]:
        result: set[str] = set()
        for case_name in case_names:
            fingerprints = reports.get(case_name, {}).get("fingerprints")
            if not isinstance(fingerprints, Mapping):
                raise ValidationError(f"pair {key} metadata is missing for {case_name}")
            value = str(fingerprints.get(key, ""))
            if len(value) != 16:
                raise ValidationError(f"pair {key} is invalid for {case_name}")
            result.add(value)
        return result

    mesh_values = values(STRICT_CASES, "mesh_fingerprint")
    if len(mesh_values) != 1:
        raise ValidationError("pair mesh_fingerprint differs across validation runs")
    species_values = values(STRICT_CASES, "species_fingerprint")
    if len(species_values) != 1:
        raise ValidationError("pair species_fingerprint differs across validation runs")
    finite_cases = (
        "smoke_finite_configured",
        "full_finite_configured_140000",
        "full_finite_configured_280000",
    )
    infinite_cases = (
        "cache_prime",
        "smoke_infinite_physical",
        "full_infinite_physical_140000",
        "full_infinite_physical_280000",
    )
    finite_models = values(finite_cases, "model_fingerprint")
    if len(finite_models) != 1:
        raise ValidationError("pair finite model_fingerprint differs within model")
    infinite_models = values(infinite_cases, "model_fingerprint")
    if len(infinite_models) != 1:
        raise ValidationError("pair infinite model_fingerprint differs within model")
    if finite_models == infinite_models:
        raise ValidationError(
            "pair finite/infinite model_fingerprint must differ with the boundary model"
        )


def _verify_submission_provenance(
    root: Path,
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    job_ids_path, _complete, identifiers = _load_submission_journal(
        root,
        require_complete=True,
    )
    current_job = os.environ.get("SLURM_JOB_ID")
    if current_job and current_job != identifiers["analysis"]:
        raise ValidationError(
            "strict analysis SLURM_JOB_ID differs from the submitted analysis job"
        )

    source = manifest.get("source")
    source_snapshot = manifest.get("source_snapshot")
    binary = manifest.get("binary")
    library = manifest.get("analysis_library")
    if not all(
        isinstance(value, Mapping)
        for value in (source, source_snapshot, binary, library)
    ):
        raise ValidationError("submission provenance metadata is incomplete")
    commit = str(source["commit"])
    log_hashes: dict[str, dict[str, str]] = {}
    for name, (job_name, resources) in EXPECTED_SUBMITTED_JOBS.items():
        job_id = identifiers[name]
        token = f"{job_id}.{job_name}"
        module_path = _require_expected_path(
            root,
            root / "provenance/modules" / f"{token}.txt",
            f"provenance/modules/{token}.txt",
            label=f"submitted {name} module log",
        )
        hash_path = _require_expected_path(
            root,
            root / "provenance/hashes" / f"{token}.sha256",
            f"provenance/hashes/{token}.sha256",
            label=f"submitted {name} hash log",
        )
        for label, path in (("module", module_path), ("hash", hash_path)):
            if not path.is_file() or path.stat().st_size == 0:
                raise ValidationError(
                    f"submitted {name} job is missing its {label} log: {path}"
                )
        module_text = module_path.read_text(encoding="utf-8", errors="replace")
        sys_modules = set(SYS_MODULE_PATTERN.findall(module_text))
        if sys_modules != {"SysA/2022"}:
            raise ValidationError(
                f"submitted {name} job module log must contain SysA/2022 only"
            )
        intel_modules = set(INTEL_MODULE_PATTERN.findall(module_text))
        if intel_modules != {"intel/2023.2", "intelmpi/2023.2"}:
            raise ValidationError(
                f"submitted {name} job module log must contain exact "
                "intel/2023.2 and intelmpi/2023.2 modules"
            )
        hash_lines = set(
            hash_path.read_text(encoding="utf-8", errors="replace").splitlines()
        )
        required_hash_paths = [str(source_snapshot["hash_file"])]
        if name == "analysis":
            required_hash_paths.append(str(library["staged_path"]))
        else:
            required_hash_paths.append(str(binary["staged_path"]))
            required_hash_paths.extend(
                str(manifest["cases"][case_name]["config_path"])
                for case_name in PRODUCER_ROLE_CASES[name]
            )
        required_hash_lines = {
            f"{required_path}: OK" for required_path in required_hash_paths
        }
        if not required_hash_lines.issubset(hash_lines):
            raise ValidationError(f"submitted {name} job hash log is incomplete")
        log_hashes[name] = {
            "module_log_sha256": _sha256(module_path),
            "hash_log_sha256": _sha256(hash_path),
        }
        if name == "analysis":
            continue
        status_path = _require_expected_path(
            root,
            root / "provenance/jobs" / f"{token}.status",
            f"provenance/jobs/{token}.status",
            label=f"submitted {name} status",
        )
        if not status_path.is_file():
            raise ValidationError(f"submitted {name} job status is missing")
        status = _summary(status_path)
        if (
            status.get("job_id") != job_id
            or status.get("job_name") != job_name
            or status.get("source_commit") != commit
            or status.get("resources") != resources
            or status.get("exit_code") != "0"
        ):
            raise ValidationError(f"submitted {name} job status is invalid")
        log_hashes[name]["status_sha256"] = _sha256(status_path)
    return {
        "status": "verified",
        "job_ids_path": str(job_ids_path),
        "job_ids_sha256": _sha256(job_ids_path),
        "job_ids": identifiers,
        "logs": log_hashes,
    }


def analyze_validation(
    archive_run: str | Path,
    validation_root: str | Path,
    *,
    library: str | Path,
    require_complete: bool = False,
) -> dict[str, Any]:
    """Write comparison artifacts and transactionally publish strict results."""

    if not require_complete:
        return _analyze_validation_impl(
            archive_run,
            validation_root,
            library=library,
            require_complete=False,
        )
    root = _lexical_absolute(validation_root)
    final_analysis = _require_expected_path(
        root,
        root / "analysis",
        "analysis",
        label="strict analysis output",
    )
    if final_analysis.exists() and any(final_analysis.iterdir()):
        raise ValidationError(
            f"strict analysis output is already published: {final_analysis}"
        )
    temporary = Path(tempfile.mkdtemp(prefix=".analysis.", dir=root))
    expected_artifacts = {
        "run_summary.csv",
        "charge_history_pair.csv",
        "particle_ledger_pair.csv",
        "mesh_potential_pair.csv",
        "snapshot_manifest.csv",
        "object_wrench.csv",
        "object_path_curves.csv",
        "object_path_summary.csv",
        "finite_shell_convergence.csv",
        "comparison_matrix.csv",
        "legacy_estimator_comparison.csv",
        "analysis_summary.json",
        "review_ja.md",
        "artifacts.json",
    }
    try:
        report = _analyze_validation_impl(
            archive_run,
            validation_root,
            library=library,
            require_complete=True,
            analysis_directory=temporary,
        )
        actual_artifacts = {path.name for path in temporary.iterdir() if path.is_file()}
        if actual_artifacts != expected_artifacts:
            raise ValidationError(
                "strict analysis artifact set mismatch: "
                f"missing={sorted(expected_artifacts - actual_artifacts)}, "
                f"extra={sorted(actual_artifacts - expected_artifacts)}"
            )
        artifact_manifest = _load_json_object_if_exists(
            temporary / "artifacts.json",
            label="strict artifact manifest",
        )
        declared_artifacts = artifact_manifest.get("artifacts")
        if not isinstance(declared_artifacts, Mapping) or set(
            declared_artifacts
        ) != expected_artifacts - {"artifacts.json"}:
            raise ValidationError("strict artifact manifest has an invalid file set")
        for name, metadata in declared_artifacts.items():
            path = temporary / name
            if (
                not isinstance(metadata, Mapping)
                or metadata.get("sha256") != _sha256(path)
                or int(metadata.get("size_bytes", -1)) != path.stat().st_size
            ):
                raise ValidationError(f"strict artifact manifest mismatch: {name}")
        for path in temporary.iterdir():
            if path.is_symlink() or not path.is_file() or path.stat().st_size == 0:
                raise ValidationError(f"strict analysis artifact is invalid: {path}")
            with path.open("rb") as stream:
                os.fsync(stream.fileno())
        directory_descriptor = os.open(temporary, os.O_RDONLY)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
        if final_analysis.exists():
            final_analysis.rmdir()
        os.replace(temporary, final_analysis)
        parent_descriptor = os.open(root, os.O_RDONLY)
        try:
            os.fsync(parent_descriptor)
        finally:
            os.close(parent_descriptor)
        qualification_status = report.get(
            "numerical_qualification_for_local_frozen_model", {}
        ).get("status")
        if qualification_status != "qualified":
            raise ValidationError(
                "strict local frozen-model numerical qualification failed; "
                f"diagnostics published at {final_analysis}"
            )
        comparison_status = report.get("comparison_artifact_contract", {}).get("status")
        if comparison_status != "complete":
            raise ValidationError(
                "strict comparison artifact contract failed; "
                f"diagnostics published at {final_analysis}"
            )
        legacy_status = report.get("legacy_estimator_audit", {}).get("status")
        if legacy_status != "complete":
            raise ValidationError(
                "strict legacy estimator audit failed; "
                f"diagnostics published at {final_analysis}"
            )
        return report
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def _analyze_validation_impl(
    archive_run: str | Path,
    validation_root: str | Path,
    *,
    library: str | Path,
    require_complete: bool = False,
    analysis_directory: Path | None = None,
) -> dict[str, Any]:
    """Write stable comparison artifacts, including explicit missing-run states."""

    archive_path = _lexical_absolute(archive_run)
    root = _lexical_absolute(validation_root)
    library_path = _lexical_absolute(library)
    manifest = _load_manifest(root)
    _require_expected_path(
        root,
        str(manifest.get("validation_root", "")),
        None,
        label="validation root",
    )
    _require_expected_path(
        archive_path,
        archive_path,
        None,
        label="archive run",
    )
    declared_archive = Path(str(manifest.get("archive", {}).get("run", "")))
    if archive_path != declared_archive:
        raise ValidationError(
            f"archive run differs from staged manifest: {archive_path} != {declared_archive}"
        )
    if not library_path.is_file():
        raise ValidationError(f"field-kernel library does not exist: {library_path}")
    if require_complete:
        _require_strict_build_origin_contract(manifest)
    strict_reports: dict[str, dict[str, Any]] = {}
    submission_provenance: dict[str, Any] | None = None
    periodic_oracle_validation: dict[str, Any] | None = None
    if require_complete:
        verify_inputs(root, require_empty_outputs=False)
        staged_library = manifest.get("analysis_library")
        if not isinstance(staged_library, dict):
            raise ValidationError(
                "complete analysis requires a staged analysis library"
            )
        staged_path = _require_descendant_path(
            root,
            str(staged_library.get("staged_path", "")),
            prefix=f"lib/{manifest.get('source', {}).get('commit', '')}",
            label="staged analysis library",
        )
        if library_path != staged_path:
            raise ValidationError(
                "complete analysis must use the staged analysis library"
            )
        if not staged_path.is_file() or _sha256(staged_path) != staged_library.get(
            "sha256"
        ):
            raise ValidationError("staged analysis library hash mismatch")
        try:
            submission_provenance = _verify_submission_provenance(root, manifest)
        except ValidationError as exc:
            raise ValidationError(
                f"complete analysis requires submission provenance: {exc}"
            ) from exc
        periodic_oracle_validation = _verify_periodic_oracle_receipt(
            root,
            library_path,
            expected_job_id=str(submission_provenance["job_ids"]["analysis"]),
        )
        strict_reports = _verify_complete_runs(root, manifest)
    analysis_inputs = manifest.get("archive", {}).get("analysis_inputs", {})
    if not isinstance(analysis_inputs, dict):
        raise ValidationError("archive analysis input metadata is invalid")
    _verify_declared_files(
        analysis_inputs,
        root=archive_path,
        label="archive analysis input",
    )
    declared_legacy_inputs = {
        relative for relative in analysis_inputs if relative in LEGACY_ESTIMATOR_INPUTS
    }
    if require_complete and declared_legacy_inputs != set(LEGACY_ESTIMATOR_INPUTS):
        raise ValidationError(
            "complete analysis requires the exact legacy estimator input set"
        )
    if require_complete:
        missing_mechanics = sorted(
            set(PRODUCTION_MECHANICS_INPUTS) - set(analysis_inputs)
        )
        if missing_mechanics:
            raise ValidationError(
                "production mechanics metadata is missing: "
                + ", ".join(missing_mechanics)
            )
        _validate_production_release_mechanics(archive_path)
    analysis_outputs = manifest.get("archive", {}).get("analysis_outputs", {})
    if not isinstance(analysis_outputs, dict):
        raise ValidationError("archive analysis output metadata is invalid")
    _verify_declared_files(
        analysis_outputs,
        root=archive_path,
        label="archive analysis output",
    )
    if require_complete:
        missing_archive_outputs = [
            relative
            for relative in ARCHIVE_REQUIRED_ANALYSIS_OUTPUTS
            if relative not in analysis_outputs
        ]
        if missing_archive_outputs:
            raise ValidationError(
                "complete analysis requires fixed archive outputs: "
                + ", ".join(missing_archive_outputs)
            )
    analysis = (
        _require_expected_path(
            root,
            root / "analysis",
            "analysis",
            label="analysis output",
        )
        if analysis_directory is None
        else _require_descendant_path(
            root,
            analysis_directory,
            label="temporary analysis output",
        )
    )
    analysis.mkdir(parents=True, exist_ok=True)
    cases = [
        (
            "archived_v1_3",
            archive_path / "work/latest",
            archive_path / "input/beach.toml",
            str(manifest.get("archive", {}).get("version", "unknown")),
        ),
        (
            "new_finite_configured",
            Path(manifest["cases"]["full_finite_configured_280000"]["output_dir"]),
            Path(manifest["cases"]["full_finite_configured_280000"]["config_path"]),
            str(manifest.get("source", {}).get("describe", "unknown")),
        ),
        (
            "new_infinite_physical",
            Path(manifest["cases"]["full_infinite_physical_280000"]["output_dir"]),
            Path(manifest["cases"]["full_infinite_physical_280000"]["config_path"]),
            str(manifest.get("source", {}).get("describe", "unknown")),
        ),
    ]
    run_rows: list[dict[str, Any]] = []
    runs: dict[str, Beach | None] = {}
    for name, output, config, version in cases:
        row, run = _case_result(name, output, config, version)
        run_rows.append(row)
        runs[name] = run
    if require_complete:
        incomplete = [
            str(row["case"]) for row in run_rows if row.get("status") != "available"
        ]
        if incomplete:
            raise ValidationError(
                "complete analysis requires all analysis cases to be available: "
                + ", ".join(incomplete)
            )
    geometry_validation = _analysis_geometry_validation(runs)
    if require_complete:
        invalid_geometry = [
            label
            for label, value in geometry_validation.items()
            if value.get("status") != "match"
        ]
        if invalid_geometry:
            raise ValidationError(
                "complete analysis requires comparable geometry: "
                + ", ".join(invalid_geometry)
            )
    _write_csv(analysis / "run_summary.csv", CSV_SCHEMAS["run_summary.csv"], run_rows)

    history_sources = {
        "archived_v1_3": [archive_path / "work/latest/charge_history.csv"],
        "new_finite_configured": [
            root / "run/full/finite_configured/140000/charge_history.csv",
            root / "run/full/finite_configured/280000/charge_history.csv",
        ],
        "new_infinite_physical": [
            root / "run/full/infinite_physical/140000/charge_history.csv",
            root / "run/full/infinite_physical/280000/charge_history.csv",
        ],
    }
    snapshot_rows: list[dict[str, Any]] = []
    snapshots: dict[str, dict[tuple[str, int], dict[str, Any]]] = {}
    for case_name, _output, _config, _version in cases:
        try:
            rows_for_case, snapshots_for_case = _collect_charge_snapshots(
                case=case_name,
                run=runs[case_name],
                history_paths=history_sources[case_name],
            )
        except ValidationError as exc:
            if require_complete:
                raise
            rows_for_case = [
                {
                    "case": case_name,
                    "status": "invalid",
                    "message": str(exc),
                }
            ]
            snapshots_for_case = {}
        snapshot_rows.extend(rows_for_case)
        snapshots[case_name] = snapshots_for_case
    _write_csv(
        analysis / "snapshot_manifest.csv",
        CSV_SCHEMAS["snapshot_manifest.csv"],
        snapshot_rows,
    )

    history_rows: list[dict[str, Any]] = []
    history_rows.extend(
        _history_rows(
            "archived_v1_3", [archive_path / "work/latest/charge_history.csv"]
        )
    )
    for model, case_name in (
        ("finite_configured", "new_finite_configured"),
        ("infinite_physical", "new_infinite_physical"),
    ):
        history_rows.extend(
            _history_rows(
                case_name,
                [
                    root / f"run/full/{model}/140000/charge_history.csv",
                    root / f"run/full/{model}/280000/charge_history.csv",
                ],
            )
        )
    _write_csv(
        analysis / "charge_history_pair.csv",
        CSV_SCHEMAS["charge_history_pair.csv"],
        history_rows,
    )

    ledger_rows: list[dict[str, Any]] = []
    potential_rows: list[dict[str, Any]] = []
    for row in run_rows:
        run = runs[str(row["case"])]
        if run is None:
            ledger_rows.append({"case": row["case"], "status": row["status"]})
            potential_rows.append({"case": row["case"], "status": row["status"]})
            continue
        result = run.result
        ledger = result.charge_ledger
        if ledger:
            for entry in ledger:
                ledger_rows.append(
                    {
                        "case": row["case"],
                        "status": "available",
                        "species_idx": entry.species_idx,
                        "batch": entry.batch,
                        "injected_count": entry.injected_count,
                        "absorbed_count": entry.absorbed_count,
                        "escaped_count": entry.escaped_count,
                        "residual_C": result.charge_ledger_residual_c,
                    }
                )
        else:
            ledger_rows.append(
                {
                    "case": row["case"],
                    "status": "not_recorded",
                    "residual_C": result.charge_ledger_residual_c or "",
                }
            )
        if result.mesh_potential_v is None:
            potential_rows.append({"case": row["case"], "status": "not_recorded"})
        else:
            for index, value in enumerate(result.mesh_potential_v, start=1):
                potential_rows.append(
                    {
                        "case": row["case"],
                        "elem_idx": index,
                        "potential_V": float(value),
                        "status": "available",
                    }
                )
    _write_csv(
        analysis / "particle_ledger_pair.csv",
        CSV_SCHEMAS["particle_ledger_pair.csv"],
        ledger_rows,
    )
    _write_csv(
        analysis / "mesh_potential_pair.csv",
        CSV_SCHEMAS["mesh_potential_pair.csv"],
        potential_rows,
    )

    (
        wrench_rows,
        object_curve_rows,
        path_summary_rows,
        shell_rows,
        physics_evaluation,
    ) = _evaluate_object_physics(
        archive_run=archive_path,
        validation_root=root,
        library=library_path,
        run_rows=run_rows,
        runs=runs,
        require_history=require_complete,
    )
    case_status = {
        str(row["case"]): {
            "status": row["status"],
            "message": row.get("message", ""),
        }
        for row in run_rows
    }
    if require_complete and physics_evaluation.get("status") != "available":
        raise ValidationError(
            "complete analysis physics status must be available with no failures"
        )
    local_qualification = _local_model_numerical_qualification(
        physics_evaluation,
        wrench_rows,
        path_summary_rows,
        shell_rows,
        case_status,
    )
    legacy_rows, legacy_audit = _legacy_estimator_audit(
        archive_path,
        wrench_rows,
        object_curve_rows,
        path_summary_rows,
        strict=require_complete,
    )
    _write_csv(
        analysis / "object_wrench.csv", CSV_SCHEMAS["object_wrench.csv"], wrench_rows
    )
    _write_csv(
        analysis / "object_path_curves.csv",
        CSV_SCHEMAS["object_path_curves.csv"],
        object_curve_rows,
    )
    _write_csv(
        analysis / "object_path_summary.csv",
        CSV_SCHEMAS["object_path_summary.csv"],
        path_summary_rows,
    )
    _write_csv(
        analysis / "finite_shell_convergence.csv",
        CSV_SCHEMAS["finite_shell_convergence.csv"],
        shell_rows,
    )
    comparisons = _comparison_rows(run_rows)
    comparisons.extend(_charge_snapshot_comparisons(snapshots, runs))
    comparisons.extend(_physics_comparison_rows(wrench_rows, path_summary_rows))
    comparison_contract = _comparison_artifact_contract(
        comparisons,
        snapshot_rows,
    )
    _write_csv(
        analysis / "comparison_matrix.csv",
        CSV_SCHEMAS["comparison_matrix.csv"],
        comparisons,
    )
    _write_csv(
        analysis / "legacy_estimator_comparison.csv",
        CSV_SCHEMAS["legacy_estimator_comparison.csv"],
        legacy_rows,
    )

    report: dict[str, Any] = {
        "schema_version": 2,
        "created_at": _utc_now(),
        "archive_run": str(archive_path),
        "validation_root": str(root),
        "library": str(library_path),
        "library_sha256": _sha256(library_path) if library_path.is_file() else None,
        "cases": case_status,
        "physics_evaluation": physics_evaluation,
        "geometry_validation": geometry_validation,
        "comparison_artifact_contract": comparison_contract,
        "legacy_estimator_audit": legacy_audit,
        "numerical_qualification_for_local_frozen_model": local_qualification,
        "snapshot_count": sum(len(values) for values in snapshots.values()),
        "strict_validation": {
            "required": require_complete,
            "verified_cases": strict_reports,
            "submission_provenance": submission_provenance,
            "periodic_plane_oracles": periodic_oracle_validation,
        },
        "comparison_semantics": {
            "archive_version_drift": "archived output to current finite model",
            "frozen_field_override": "same charge snapshot with only field closure changed",
            "boundary_history_response_common_evaluator": "different charging histories under the same infinite evaluator",
            "actual_end_to_end": "charging-history and field-closure effects combined",
        },
        "limitations": [
            "archived output spans parent/child executables and the exact binaries are unavailable",
            "one paired seed is not a release probability estimate",
            "failed object force/path/shell evaluations remain explicit invalid/not_evaluated rows",
            "frozen charge and contact/adhesion assumptions must be stated separately",
            "vdw_work release uses an equivalent constant-force profile that preserves initial force and total work, not the original 1/s^2 barrier shape",
            "E_bottom=0 plane and infinite_physical results are closure-dependent and are not universal free-space self forces",
            "0..2R frozen-field work and speed are local detachment diagnostics, not escape-to-infinity energy or speed",
            "saved sphere mesh and source discretization refinement remain not_evaluated",
            "plane-oracle errors are not used as sphere force or torque error bars",
        ],
    }
    _write_json(analysis / "analysis_summary.json", report)
    finite_status = case_status["new_finite_configured"]["status"]
    infinite_status = case_status["new_infinite_physical"]["status"]
    physics_review = _physics_review_lines(wrench_rows, path_summary_rows, shell_rows)
    review = f"""# 周期境界 object-force 比較レビュー

## 実行状態

- 旧 archive: {case_status["archived_v1_3"]["status"]}
- 新 finite: {"未実行" if finite_status == "missing" else finite_status}
- 新 infinite: {"未実行" if infinite_status == "missing" else infinite_status}

## Object force・work・release

{physics_review}

- structural physics evaluation: {physics_evaluation.get("status")}
- fixed-discretization path/work/shell gate: {local_qualification.get("status")}
- 保存済み sphere mesh/source refinement: {local_qualification.get("saved_sphere_mesh_refinement_status")}
- sphere absolute force error: {local_qualification.get("sphere_absolute_force_error_status")}
- sphere absolute torque error: {local_qualification.get("sphere_absolute_torque_error_status")}
- 平面 oracle の誤差を sphere error bar に流用しません。
- claim scope: {local_qualification.get("claim_scope")}

## Legacy estimator audit

- status: {legacy_audit.get("status")}
- compared rows: {legacy_audit.get("comparison_row_count")}
- matching batches/meshes: {legacy_audit.get("covered_native_keys")}
- numerical closeness is a gate: {legacy_audit.get("closeness_is_a_gate")}

旧 moving-sphere/direct/proxy estimator は target 全体を self source から除外し、current
native evaluator は primary だけを除いて target periodic images を保持します。したがって
差分は機械可読に保存しますが、値の近さ自体を合否条件にはしません。

## 解釈

旧結果は二つの実行バイナリにまたがり、当時の実体も保存されていません。したがって
`archive -> new finite` は境界条件だけでなく version/restart/runtime drift を含みます。
`new finite -> new infinite` だけを境界モデル差として扱います。

`infinite_physical` は有限像の単純な延長ではなく、条件付き周期和と
`E_bottom=0` zero-mode closure を同時に選びます。一様平面 oracle はこの closure の下での
Maxwell traction であり、自由空間の普遍的 self-force ではありません。上方一定場が残る
非中性 case では `0..2R` の work/speed は局所 frozen-field 指標であり、無限遠 escape
energy/speed ではありません。

`not_evaluated` / `invalid` の項目は結論に使いません。実評価済みでも frozen charge、
接触・接着 bracket、一つの乱数 seed という制約を残します。
保存済み sphere mesh/source discretization の refinement と sphere force/torque の絶対誤差は
`not_evaluated` です。平面 oracle の誤差を sphere error bar に流用しません。
`vdw_work` は初期力と全仕事を保存する等価な定数力 profile であり、元の `1/s^2`
障壁形状そのものではありません。
"""
    (analysis / "review_ja.md").write_text(review, encoding="utf-8")
    artifacts = {
        path.name: {"sha256": _sha256(path), "size_bytes": path.stat().st_size}
        for path in sorted(analysis.iterdir())
        if path.is_file() and path.name != "artifacts.json"
    }
    _write_json(
        analysis / "artifacts.json", {"schema_version": 1, "artifacts": artifacts}
    )
    return report


def probe_library(library: str | Path) -> dict[str, Any]:
    """Load the native ABI and execute one finite field/potential evaluation."""

    library_path = Path(library).resolve()
    if not library_path.is_file():
        raise ValidationError(f"field-kernel library does not exist: {library_path}")
    source_triangles = np.array(
        [[[0.0, -5.0e-5, -5.0e-5], [0.0, 5.0e-5, -5.0e-5], [0.0, 0.0, 5.0e-5]]],
        dtype=np.float64,
    )
    source_charges = np.array([1.0e-15], dtype=np.float64)
    targets = np.array([[1.0e-3, 0.0, 0.0]], dtype=np.float64)
    try:
        with FieldKernel(
            source_triangles,
            source_charges,
            library_path=library_path,
        ) as kernel:
            electric_field = np.asarray(kernel.eval_e(targets), dtype=np.float64)
            potential = np.asarray(kernel.eval_phi(targets), dtype=np.float64)
    except Exception as exc:
        raise ValidationError(f"field-kernel ABI smoke failed: {exc}") from exc
    if (
        electric_field.shape != (1, 3)
        or potential.shape != (1,)
        or not np.all(np.isfinite(electric_field))
        or not np.all(np.isfinite(potential))
        or float(np.linalg.norm(electric_field[0])) <= 0.0
        or float(abs(potential[0])) <= 0.0
    ):
        raise ValidationError("field-kernel ABI smoke returned invalid values")
    return {
        "status": "ok",
        "library": str(library_path),
        "library_sha256": _sha256(library_path),
        "electric_field_V_m": electric_field[0].tolist(),
        "potential_V": float(potential[0]),
    }


def _oracle_plane_triangles(
    cells_per_axis: int,
    *,
    length: float = 2.0,
) -> np.ndarray:
    triangles: list[np.ndarray] = []
    spacing = length / cells_per_axis
    for iy in range(cells_per_axis):
        for ix in range(cells_per_axis):
            x0 = ix * spacing
            x1 = (ix + 1) * spacing
            y0 = iy * spacing
            y1 = (iy + 1) * spacing
            p00 = np.array([x0, y0, 0.0])
            p10 = np.array([x1, y0, 0.0])
            p11 = np.array([x1, y1, 0.0])
            p01 = np.array([x0, y1, 0.0])
            triangles.extend((np.array([p00, p10, p11]), np.array([p00, p11, p01])))
    return np.asarray(triangles)


def _oracle_panel_result(
    directory: Path,
    triangles: np.ndarray,
    charges: np.ndarray,
    *,
    mesh_ids: np.ndarray | None = None,
) -> FortranRunResult:
    return FortranRunResult(
        directory=directory,
        mesh_nelem=len(charges),
        processed_particles=0,
        absorbed=0,
        escaped=0,
        batches=0,
        escaped_boundary=0,
        survived_max_step=0,
        last_rel_change=0.0,
        charges=np.asarray(charges, dtype=np.float64),
        triangles=np.asarray(triangles, dtype=np.float64),
        mesh_ids=(
            np.ones(len(charges), dtype=np.int64)
            if mesh_ids is None
            else np.asarray(mesh_ids, dtype=np.int64)
        ),
        field_source_model="triangle_p0",
        field_kernel_id="triangle_p0_exact_p2m_near",
    )


def _oracle_panel_config_data() -> dict[str, Any]:
    return {
        "sim": {
            "field_solver": "fmm",
            "field_bc_mode": "periodic2",
            "bc_x_low": "periodic",
            "bc_x_high": "periodic",
            "bc_y_low": "periodic",
            "bc_y_high": "periodic",
            "bc_z_low": "open",
            "bc_z_high": "open",
            "box_min": [0.0, 0.0, -1.0],
            "box_max": [2.0, 2.0, 1.0],
            "tree_theta": 0.2,
            "tree_leaf_max": 16,
            "tree_order": 4,
            "field_periodic_image_layers": 1,
            "field_periodic_far_correction": "none",
            "field_periodic_ewald_layers": 4,
            "field_periodic_generation_tolerance": 1.0e-8,
            "e0": [0.0, 0.0, 0.0],
        },
        "periodic2": {
            "nonzero_mode_backend": "cached_kneq0",
            "zero_mode_policy": "exclude_k0",
            "lower_boundary_model": "e_bottom_zero",
        },
    }


def _oracle_panel_config(path: Path) -> None:
    _write_toml(path, _oracle_panel_config_data())


def _oracle_nonnegative_number(value: Any, *, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"periodic plane-oracle {label} is invalid") from exc
    if not math.isfinite(number) or number < 0.0:
        raise ValidationError(f"periodic plane-oracle {label} is invalid")
    return number


def _oracle_integer(value: Any, *, label: str) -> int:
    try:
        number = int(value)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"periodic plane-oracle {label} is invalid") from exc
    if isinstance(value, float) and not value.is_integer():
        raise ValidationError(f"periodic plane-oracle {label} is invalid")
    return number


def _verify_periodic_oracle_metrics(
    kernel_oracle: Mapping[str, Any],
    *,
    label: str,
) -> None:
    if kernel_oracle.get("effective_field_contract") != ORACLE_EFFECTIVE_FIELD_CONTRACT:
        raise ValidationError(
            f"periodic plane-oracle {label} effective field contract is invalid"
        )
    cache_diagnostics = kernel_oracle.get("cache_diagnostics")
    if (
        not isinstance(cache_diagnostics, Mapping)
        or cache_diagnostics.get("hit") is not True
        or type(cache_diagnostics.get("build_count")) is not int
        or cache_diagnostics.get("build_count") != 0
        or not str(cache_diagnostics.get("fingerprint") or "")
        or not str(cache_diagnostics.get("path") or "")
        or len(str(cache_diagnostics.get("sha256") or "")) != 64
    ):
        raise ValidationError(
            f"periodic plane-oracle {label} cache diagnostics are invalid"
        )
    uniform = kernel_oracle.get("uniform_plane")
    cosine = kernel_oracle.get("neutral_cosine_plane")
    if not isinstance(uniform, Mapping) or not isinstance(cosine, Mapping):
        raise ValidationError(f"periodic plane-oracle {label} metrics are missing")
    uniform_tolerance = _oracle_nonnegative_number(
        uniform.get("relative_tolerance"),
        label=f"{label} uniform tolerance",
    )
    quadrature_tolerance = _oracle_nonnegative_number(
        uniform.get("quadrature_relative_tolerance"),
        label=f"{label} quadrature tolerance",
    )
    if (
        uniform.get("closure") != "e_bottom_zero"
        or uniform.get("interpretation")
        != "Maxwell traction under E_bottom=0; not universal free-space self force"
        or _oracle_integer(
            uniform.get("cells_per_axis"),
            label=f"{label} uniform cells_per_axis",
        )
        != 4
        or uniform_tolerance != ORACLE_UNIFORM_RELATIVE_TOLERANCE
        or quadrature_tolerance != ORACLE_QUADRATURE_RELATIVE_TOLERANCE
        or _oracle_integer(
            uniform.get("field_cells_per_axis"),
            label=f"{label} uniform field_cells_per_axis",
        )
        != 4
    ):
        raise ValidationError(
            f"periodic uniform-plane oracle {label} contract is invalid"
        )
    uniform_errors = [
        _oracle_nonnegative_number(
            uniform.get("below_absolute_error_V_m"),
            label="uniform below-plane error",
        ),
        _oracle_nonnegative_number(
            uniform.get("nonzero_relative_error"),
            label="uniform nonzero-field error",
        ),
        _oracle_nonnegative_number(
            uniform.get("transverse_absolute_error_V_m"),
            label="uniform transverse-field error",
        ),
        _oracle_nonnegative_number(
            uniform.get("potential_gauge_absolute_error_V"),
            label="uniform potential gauge error",
        ),
        _oracle_nonnegative_number(
            uniform.get("potential_nonzero_relative_error"),
            label="uniform nonzero potential error",
        ),
    ]
    force_errors = uniform.get("force_relative_errors")
    transverse_errors = uniform.get("force_transverse_relative_errors")
    torque_errors = uniform.get("torque_relative_errors")
    if not all(
        isinstance(values, list) and len(values) >= 1
        for values in (force_errors, transverse_errors, torque_errors)
    ):
        raise ValidationError(
            f"periodic uniform-plane oracle {label} wrench errors are invalid"
        )
    wrench_errors = [
        _oracle_nonnegative_number(value, label="uniform wrench error")
        for values in (force_errors, transverse_errors, torque_errors)
        for value in values
    ]
    quadrature_error = _oracle_nonnegative_number(
        uniform.get("quadrature_relative_difference"),
        label="uniform quadrature difference",
    )
    if (
        max((*uniform_errors, *wrench_errors)) > ORACLE_UNIFORM_RELATIVE_TOLERANCE
        or quadrature_error > ORACLE_QUADRATURE_RELATIVE_TOLERANCE
    ):
        raise ValidationError(
            f"periodic uniform-plane oracle {label} exceeds its thresholds: "
            f"field_errors={uniform_errors}, wrench_errors={wrench_errors}, "
            f"quadrature_error={quadrature_error}"
        )

    cosine_tolerance = _oracle_nonnegative_number(
        cosine.get("fine_relative_tolerance"),
        label=f"{label} cosine tolerance",
    )
    sample_abs_z = cosine.get("sample_abs_z_m")
    if not isinstance(sample_abs_z, list) or len(sample_abs_z) != 2:
        raise ValidationError("periodic cosine-plane sample heights are invalid")
    parsed_sample_abs_z = [
        _oracle_nonnegative_number(value, label="cosine sample height")
        for value in sample_abs_z
    ]
    expected_decay_ratio = _oracle_nonnegative_number(
        cosine.get("expected_decay_ratio"),
        label="cosine expected decay ratio",
    )
    decay_ratio_tolerance = _oracle_nonnegative_number(
        cosine.get("decay_ratio_relative_tolerance"),
        label="cosine decay ratio tolerance",
    )
    if (
        cosine.get("expected_decay") != "exp(-k*abs(z-z0))"
        or cosine_tolerance != ORACLE_COSINE_FINE_RELATIVE_TOLERANCE
        or parsed_sample_abs_z != list(ORACLE_COSINE_SAMPLE_ABS_Z_M)
        or not math.isclose(
            expected_decay_ratio,
            ORACLE_COSINE_EXPECTED_DECAY_RATIO,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        )
        or decay_ratio_tolerance != ORACLE_COSINE_DECAY_RATIO_RELATIVE_TOLERANCE
    ):
        raise ValidationError("periodic cosine-plane oracle contract is invalid")
    errors = cosine.get("errors")
    if not isinstance(errors, list) or len(errors) != 2:
        raise ValidationError("periodic cosine-plane errors are invalid")
    parsed: list[tuple[int, float, float, float, float, float, float, float]] = []
    for value in errors:
        if not isinstance(value, Mapping):
            raise ValidationError("periodic cosine-plane error row is invalid")
        parsed.append(
            (
                _oracle_integer(
                    value.get("cells_per_axis"),
                    label=f"{label} cosine cells_per_axis",
                ),
                _oracle_nonnegative_number(
                    value.get("field_relative_error"),
                    label="cosine field error",
                ),
                _oracle_nonnegative_number(
                    value.get("potential_relative_error"),
                    label="cosine potential error",
                ),
                _oracle_nonnegative_number(
                    value.get("field_decay_ratio"),
                    label="cosine field decay ratio",
                ),
                _oracle_nonnegative_number(
                    value.get("potential_decay_ratio"),
                    label="cosine potential decay ratio",
                ),
                _oracle_nonnegative_number(
                    value.get("field_decay_ratio_relative_error"),
                    label="cosine field decay ratio error",
                ),
                _oracle_nonnegative_number(
                    value.get("potential_decay_ratio_relative_error"),
                    label="cosine potential decay ratio error",
                ),
                _oracle_nonnegative_number(
                    value.get("charge_neutrality_ratio"),
                    label="cosine charge neutrality",
                ),
            )
        )
    ratio_records_are_consistent = all(
        math.isclose(
            field_ratio_error,
            abs(field_ratio - expected_decay_ratio) / expected_decay_ratio,
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        )
        and math.isclose(
            potential_ratio_error,
            abs(potential_ratio - expected_decay_ratio) / expected_decay_ratio,
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        )
        for (
            _cells,
            _field_error,
            _potential_error,
            field_ratio,
            potential_ratio,
            field_ratio_error,
            potential_ratio_error,
            _charge_neutrality,
        ) in parsed
    )
    if (
        [value[0] for value in parsed] != [4, 8]
        or parsed[1][1] >= parsed[0][1]
        or parsed[1][2] >= parsed[0][2]
        or parsed[1][1] > ORACLE_COSINE_FINE_RELATIVE_TOLERANCE
        or parsed[1][2] > ORACLE_COSINE_FINE_RELATIVE_TOLERANCE
        or not ratio_records_are_consistent
        or max(parsed[1][5], parsed[1][6])
        > ORACLE_COSINE_DECAY_RATIO_RELATIVE_TOLERANCE
        or max(value[7] for value in parsed) > 1.0e-12
    ):
        raise ValidationError(
            f"periodic cosine-plane oracle {label} exceeds its thresholds"
        )


def _verify_periodic_oracle_receipt(
    validation_root: Path,
    library: Path,
    *,
    expected_job_id: str | None = None,
    _candidate_receipt: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    receipt_path = _require_expected_path(
        validation_root,
        validation_root / "provenance/oracles/periodic_plane.json",
        "provenance/oracles/periodic_plane.json",
        label="periodic plane-oracle receipt",
    )
    if _candidate_receipt is None:
        if not receipt_path.is_file():
            raise ValidationError("periodic plane-oracle receipt is missing")
        receipt = _load_execution_receipt(receipt_path)
    else:
        receipt = dict(_candidate_receipt)
    manifest_path = _require_expected_path(
        validation_root,
        validation_root / "manifest.json",
        "manifest.json",
        label="validation manifest",
    )
    manifest = _load_manifest(validation_root)
    _require_expected_path(
        validation_root,
        str(manifest.get("validation_root", "")),
        None,
        label="validation root",
    )
    library = _require_descendant_path(
        validation_root,
        library,
        prefix="lib",
        label="periodic plane-oracle library",
    )
    cache_dir = _require_expected_path(
        validation_root,
        validation_root / "cache/oracles",
        "cache/oracles",
        label="periodic plane-oracle cache directory",
    )
    config_path = _require_expected_path(
        validation_root,
        validation_root / "provenance/oracles/periodic_plane.toml",
        "provenance/oracles/periodic_plane.toml",
        label="periodic plane-oracle config",
    )
    if (
        not config_path.is_file()
        or Path(str(receipt.get("config", ""))) != config_path
        or receipt.get("config_sha256") != _sha256(config_path)
    ):
        raise ValidationError("periodic plane-oracle config no longer matches staging")
    if (
        receipt.get("oracle_schema_version") != 3
        or receipt.get("status") != "qualified"
        or receipt.get("manifest_sha256") != _sha256(manifest_path)
        or receipt.get("library_sha256") != _sha256(library)
        or Path(str(receipt.get("library", ""))) != library
        or Path(str(receipt.get("cache_dir", ""))) != cache_dir
        or receipt.get("cache_files") != _output_inventory(cache_dir)
        or (
            expected_job_id is not None
            and str(receipt.get("execution_job_id", "")) != expected_job_id
        )
    ):
        raise ValidationError("periodic plane-oracle receipt no longer matches staging")
    expected_build_origin = manifest.get("build_origin")
    if expected_build_origin is not None:
        if (
            receipt.get("library_build_origin") != expected_build_origin
            or _library_build_info(library) != expected_build_origin
        ):
            raise ValidationError(
                "periodic plane-oracle library build origin no longer matches staging"
            )
    difference = _first_difference(
        _load_toml(config_path),
        _oracle_panel_config_data(),
    )
    if difference is not None:
        raise ValidationError(f"periodic plane-oracle config mismatch: {difference}")
    kernel_configs = receipt.get("kernel_configs")
    expected_configs = {
        "triangle_p0": (
            config_path,
            _oracle_panel_config_data(),
        ),
    }
    if not isinstance(kernel_configs, Mapping) or set(kernel_configs) != set(
        expected_configs
    ):
        raise ValidationError("periodic plane-oracle kernel config set is invalid")
    for label, (expected_path, expected_data) in expected_configs.items():
        metadata = kernel_configs[label]
        if not isinstance(metadata, Mapping):
            raise ValidationError(
                f"periodic plane-oracle {label} config metadata is invalid"
            )
        expected_path = _require_expected_path(
            validation_root,
            str(metadata.get("path", "")),
            expected_path.relative_to(validation_root).as_posix(),
            label=f"periodic plane-oracle {label} config",
        )
        if not expected_path.is_file() or metadata.get("sha256") != _sha256(
            expected_path
        ):
            raise ValidationError(
                f"periodic plane-oracle {label} config no longer matches staging"
            )
        config_difference = _first_difference(
            _load_toml(expected_path),
            expected_data,
        )
        if config_difference is not None:
            raise ValidationError(
                f"periodic plane-oracle {label} config mismatch: {config_difference}"
            )
    kernel_oracles = receipt.get("kernel_oracles")
    if not isinstance(kernel_oracles, Mapping):
        raise ValidationError("periodic plane-oracle kernel_oracles are missing")
    expected_models = {
        "triangle_p0": ("triangle_p0", "triangle_p0_exact_p2m_near"),
    }
    missing_models = sorted(set(expected_models) - set(kernel_oracles))
    extra_models = sorted(set(kernel_oracles) - set(expected_models))
    if missing_models or extra_models:
        raise ValidationError(
            "periodic plane-oracle kernel model set mismatch: "
            f"missing={missing_models}, extra={extra_models}"
        )
    verified_kernel_oracles: dict[str, Any] = {}
    for label, (source_model, kernel_id) in expected_models.items():
        value = kernel_oracles[label]
        if not isinstance(value, Mapping):
            raise ValidationError(f"periodic plane-oracle {label} record is invalid")
        if value.get("field_source_model") != source_model:
            raise ValidationError(
                f"periodic plane-oracle {label} field_source_model is invalid"
            )
        if value.get("field_kernel_id") != kernel_id:
            raise ValidationError(
                f"periodic plane-oracle {label} field_kernel_id is invalid"
            )
        _verify_periodic_oracle_metrics(value, label=label)
        diagnostics = value["cache_diagnostics"]
        assert isinstance(diagnostics, Mapping)
        cache_identities = value.get("cache_identities")
        if (
            not isinstance(cache_identities, Mapping)
            or set(cache_identities) != {"4", "8"}
            or not all(
                isinstance(identity, Mapping) for identity in cache_identities.values()
            )
        ):
            raise ValidationError(
                f"periodic plane-oracle {label} cache identity groups are invalid"
            )
        if _first_difference(diagnostics, cache_identities["4"]) is not None:
            raise ValidationError(
                f"periodic plane-oracle {label} cache diagnostics alias is invalid"
            )
        canonical_by_group: dict[str, tuple[str, Path, str]] = {}
        for group, identity in cache_identities.items():
            if set(identity) != {
                "hit",
                "build_count",
                "fingerprint",
                "path",
                "sha256",
            }:
                raise ValidationError(
                    f"periodic plane-oracle {label} cache group {group} schema is invalid"
                )
            cache_path = _require_descendant_path(
                validation_root,
                str(identity.get("path", "")),
                prefix="cache/oracles",
                label=(f"periodic plane-oracle {label} cache group {group} path"),
            )
            if (
                identity.get("hit") is not True
                or type(identity.get("build_count")) is not int
                or identity.get("build_count") != 0
                or not str(identity.get("fingerprint", ""))
                or not cache_path.is_file()
                or identity.get("sha256") != _sha256(cache_path)
            ):
                raise ValidationError(
                    f"periodic plane-oracle {label} cache group {group} identity is invalid"
                )
            canonical_by_group[str(group)] = (
                str(identity["fingerprint"]),
                cache_path,
                str(identity["sha256"]),
            )
        if (
            canonical_by_group["4"][0] == canonical_by_group["8"][0]
            or canonical_by_group["4"][1] == canonical_by_group["8"][1]
        ):
            raise ValidationError(
                f"periodic plane-oracle {label} 4/8 cache identities collide"
            )
        cache_evaluations = value.get("cache_evaluations")
        expected_cache_labels = ORACLE_CACHE_EVALUATION_LABELS
        if (
            not isinstance(cache_evaluations, list)
            or not all(isinstance(row, Mapping) for row in cache_evaluations)
            or tuple(row.get("label") for row in cache_evaluations)
            != expected_cache_labels
        ):
            raise ValidationError(
                f"periodic plane-oracle {label} cache evaluation labels are invalid"
            )
        cache_group_for_label = {
            evaluation_label: group
            for group, evaluation_labels in ORACLE_CACHE_EVALUATION_GROUPS.items()
            for evaluation_label in evaluation_labels
        }
        for evaluation in cache_evaluations:
            if set(evaluation) != {
                "label",
                "hit",
                "build_count",
                "fingerprint",
                "path",
                "sha256",
            }:
                raise ValidationError(
                    f"periodic plane-oracle {label} cache evaluation schema is invalid"
                )
            evaluation_path = _require_descendant_path(
                validation_root,
                str(evaluation.get("path", "")),
                prefix="cache/oracles",
                label=f"periodic plane-oracle {label} cache evaluation path",
            )
            evaluation_identity = (
                str(evaluation.get("fingerprint", "")),
                evaluation_path,
                str(evaluation.get("sha256", "")),
            )
            group = cache_group_for_label.get(str(evaluation.get("label", "")))
            if (
                group is None
                or evaluation.get("hit") is not True
                or type(evaluation.get("build_count")) is not int
                or evaluation.get("build_count") != 0
                or evaluation_identity != canonical_by_group[group]
                or not evaluation_path.is_file()
                or evaluation.get("sha256") != _sha256(evaluation_path)
            ):
                raise ValidationError(
                    f"periodic plane-oracle {label} cache evaluation identity is invalid"
                )
        if label == "triangle_p0":
            uniform = value["uniform_plane"]
            assert isinstance(uniform, Mapping)
            wrench_rows = uniform.get("wrench_refinement")
            force_errors = uniform.get("force_relative_errors")
            transverse_errors = uniform.get("force_transverse_relative_errors")
            torque_errors = uniform.get("torque_relative_errors")
            if (
                uniform.get("target_integration") != "gauss_duffy_order_3_and_7"
                or uniform.get("quadrature_order") != [3, 7]
                or not isinstance(wrench_rows, list)
                or len(wrench_rows) != 2
                or not all(isinstance(row, Mapping) for row in wrench_rows)
                or not all(
                    isinstance(errors, list) and len(errors) == 2
                    for errors in (force_errors, transverse_errors, torque_errors)
                )
            ):
                raise ValidationError(
                    "periodic triangle_p0 integration structure is invalid"
                )
            parsed_triangle_wrenches = [
                (
                    _oracle_integer(
                        row.get("cells_per_axis"),
                        label="triangle_p0 wrench cells",
                    ),
                    _oracle_nonnegative_number(
                        row.get("force_relative_error"),
                        label="triangle_p0 wrench force error",
                    ),
                    _oracle_nonnegative_number(
                        row.get("transverse_relative_error"),
                        label="triangle_p0 wrench transverse error",
                    ),
                    _oracle_nonnegative_number(
                        row.get("torque_relative_error"),
                        label="triangle_p0 wrench torque error",
                    ),
                    _oracle_nonnegative_number(
                        row.get("component_consistency_relative_error"),
                        label="triangle_p0 component consistency",
                    ),
                    _oracle_nonnegative_number(
                        row.get("other_objects_normalized_absolute"),
                        label="triangle_p0 other-object component",
                    ),
                    _oracle_nonnegative_number(
                        row.get("external_uniform_normalized_absolute"),
                        label="triangle_p0 external component",
                    ),
                    _oracle_nonnegative_number(
                        row.get("total_minus_images_normalized_absolute"),
                        label="triangle_p0 total/image identity",
                    ),
                    _oracle_nonnegative_number(
                        row.get("primary_free_subtraction_normalized_absolute"),
                        label="triangle_p0 primary subtraction",
                    ),
                )
                for row in wrench_rows
            ]
            if (
                [row[0] for row in parsed_triangle_wrenches] != [4, 4]
                or [row[1] for row in parsed_triangle_wrenches]
                != [float(item) for item in force_errors]
                or [row[2] for row in parsed_triangle_wrenches]
                != [float(item) for item in transverse_errors]
                or [row[3] for row in parsed_triangle_wrenches]
                != [float(item) for item in torque_errors]
                or any(
                    max(other, external, total_minus_images) > 1.0e-12
                    or component > ORACLE_UNIFORM_RELATIVE_TOLERANCE
                    or not math.isclose(
                        component,
                        max(other, external, total_minus_images, primary),
                        rel_tol=0.0,
                        abs_tol=1.0e-15,
                    )
                    for (
                        _cells,
                        _force,
                        _transverse,
                        _torque,
                        component,
                        other,
                        external,
                        total_minus_images,
                        primary,
                    ) in parsed_triangle_wrenches
                )
            ):
                raise ValidationError(
                    "periodic triangle_p0 integration structure is invalid"
                )
        verified_kernel_oracles[label] = dict(value)
    return {
        "status": "qualified",
        "receipt_path": (str(receipt_path) if _candidate_receipt is None else None),
        "receipt_sha256": (
            _sha256(receipt_path) if _candidate_receipt is None else None
        ),
        "execution_job_id": receipt.get("execution_job_id"),
        "kernel_oracles": verified_kernel_oracles,
    }


def _run_periodic_plane_kernel_oracle(
    *,
    root: Path,
    config_path: Path,
    cache_dir: Path,
    library_path: Path,
) -> dict[str, Any]:
    eps0 = 1.0 / (4.0 * math.pi * K_COULOMB)

    area_xy = 4.0
    plane_triangles_by_cells = {
        cells: _oracle_plane_triangles(cells) for cells in (4, 8)
    }

    expected_field_z = np.array([0.0, 0.5, 1.0])
    canonical_cache_identities: dict[str, dict[str, Any]] = {}
    cache_evaluations: list[dict[str, Any]] = []
    observed_effective_field_contract: dict[str, str] | None = None
    oracle_sim = _load_toml(config_path).get("sim", {})
    if not isinstance(oracle_sim, Mapping):
        raise ValidationError("periodic plane-oracle sim config is invalid")
    cache_groups = ORACLE_CACHE_EVALUATION_GROUPS
    cache_group_for_label = {
        label: group for group, labels in cache_groups.items() for label in labels
    }

    def record_effective_contract(snapshot: ObjectInteractionSnapshot) -> None:
        nonlocal observed_effective_field_contract
        periodic2 = snapshot._options.periodic2
        effective_far_correction = "free" if periodic2 is None else str(periodic2[4])
        has_zero_mode = snapshot._zero_mode is not None
        actual = {
            "requested_periodic_model": str(snapshot.periodic_model),
            "configured_far_correction": str(
                oracle_sim.get("field_periodic_far_correction", "none")
            ),
            "effective_far_correction": effective_far_correction,
            "nonzero_mode_backend": effective_far_correction,
            "zero_mode_policy": "exclude_k0" if has_zero_mode else "unavailable",
            "lower_boundary_model": (
                "e_bottom_zero" if has_zero_mode else "unavailable"
            ),
        }
        if actual != ORACLE_EFFECTIVE_FIELD_CONTRACT:
            raise ValidationError(
                "periodic triangle_p0 plane-oracle resolved an invalid "
                f"effective field contract: {actual}"
            )
        if (
            observed_effective_field_contract is not None
            and actual != observed_effective_field_contract
        ):
            raise ValidationError(
                "periodic plane-oracle effective field contract changed during execution"
            )
        observed_effective_field_contract = actual

    def record_cache_evaluation(label: str, kernel: FieldKernel) -> None:
        diagnostics = kernel.diagnostics()
        cache_path = diagnostics.periodic_cache_path
        if (
            diagnostics.periodic_cache_hit is not True
            or diagnostics.periodic_operator_build_count != 0
            or not diagnostics.periodic_cache_fingerprint
            or cache_path is None
            or not cache_path.is_file()
        ):
            raise ValidationError(
                "periodic triangle_p0 plane-oracle cache evaluation "
                f"{label} was not a warm reuse"
            )
        evaluation = {
            "label": label,
            "hit": True,
            "build_count": 0,
            "fingerprint": diagnostics.periodic_cache_fingerprint,
            "path": str(cache_path.resolve()),
            "sha256": _sha256(cache_path),
        }
        identity = {key: evaluation[key] for key in ("fingerprint", "path", "sha256")}
        group = cache_group_for_label.get(label)
        if group is None:
            raise ValidationError(
                "periodic triangle_p0 plane-oracle cache evaluation "
                f"{label} has no resolution group"
            )
        canonical = canonical_cache_identities.get(group)
        if canonical is None:
            canonical_cache_identities[group] = {
                "hit": True,
                "build_count": 0,
                **identity,
            }
        elif identity != {
            key: canonical[key] for key in ("fingerprint", "path", "sha256")
        }:
            raise ValidationError(
                "periodic triangle_p0 plane-oracle cache evaluation "
                f"{label} changed its resolution-group operator identity"
            )
        cache_evaluations.append(evaluation)

    def prime_cache(cells_per_axis: int) -> None:
        triangles = plane_triangles_by_cells[cells_per_axis]
        charges = np.full(
            triangles.shape[0],
            eps0 * area_xy / triangles.shape[0],
        )
        with ObjectInteractionSnapshot.from_result(
            _oracle_panel_result(
                root,
                triangles,
                charges,
            ),
            step=None,
            config_path=config_path,
            periodic_model="infinite_physical",
            cache_dir=cache_dir,
            library_path=library_path,
        ) as snapshot:
            diagnostics = snapshot._periodic.diagnostics()
            if (
                diagnostics.periodic_cache_hit not in (True, False)
                or diagnostics.periodic_operator_build_count
                != (0 if diagnostics.periodic_cache_hit else 1)
                or not diagnostics.periodic_cache_fingerprint
                or diagnostics.periodic_cache_path is None
                or not diagnostics.periodic_cache_path.is_file()
            ):
                raise ValidationError(
                    "periodic triangle_p0 plane-oracle cache prime failed"
                )

    def evaluate_uniform(
        cells_per_axis: int,
        points: np.ndarray,
        *,
        cache_label: str,
        quadrature_orders: Sequence[int] | None,
    ) -> tuple[np.ndarray, np.ndarray, list[Any], float]:
        triangles = plane_triangles_by_cells[cells_per_axis]
        triangle_area = area_xy / triangles.shape[0]
        charges = np.full(triangles.shape[0], eps0 * triangle_area)
        with ObjectInteractionSnapshot.from_result(
            _oracle_panel_result(
                root,
                triangles,
                charges,
            ),
            step=None,
            config_path=config_path,
            periodic_model="infinite_physical",
            cache_dir=cache_dir,
            library_path=library_path,
        ) as snapshot:
            record_effective_contract(snapshot)
            field = snapshot._periodic.eval_e(points)
            potential = snapshot._periodic.eval_phi(points)
            if snapshot._zero_mode is None:
                raise ValidationError("uniform plane oracle has no declared zero mode")
            zero_phi, zero_ez = snapshot._zero_mode.eval(
                points[:, 2],
                trace="plus",
            )
            field[:, 2] += zero_ez
            potential += zero_phi
            if quadrature_orders is None:
                evaluated_wrenches = [snapshot.object_probe(1).wrench()]
            else:
                evaluated_wrenches = [
                    snapshot.object_probe(1, quadrature_order=order).wrench()
                    for order in quadrature_orders
                ]
            record_cache_evaluation(cache_label, snapshot._periodic)
        return field, potential, evaluated_wrenches, float(np.sum(charges))

    prime_cache(4)
    prime_cache(8)

    uniform_points = np.array(
        [
            [0.37, 0.61, -0.25],
            [0.37, 0.61, 0.0],
            [0.37, 0.61, 0.25],
        ]
    )
    uniform_field, uniform_potential, wrenches, total_charge = evaluate_uniform(
        4,
        uniform_points,
        cache_label="uniform_4",
        quadrature_orders=(3, 7),
    )
    wrench_cells = [4, 4]
    below_absolute_error = float(abs(uniform_field[0, 2]))
    nonzero_relative_error = float(
        np.max(
            np.abs(uniform_field[1:, 2] - expected_field_z[1:])
            / np.abs(expected_field_z[1:])
        )
    )
    transverse_absolute_error = float(np.max(np.abs(uniform_field[:, :2])))
    expected_potential = np.array([0.0, 0.0, -0.25])
    potential_gauge_absolute_error = float(
        np.max(np.abs(uniform_potential[:2] - expected_potential[:2]))
    )
    potential_nonzero_relative_error = float(
        abs(uniform_potential[2] - expected_potential[2]) / abs(expected_potential[2])
    )
    expected_force = total_charge / 2.0
    force_relative_errors = [
        float(abs(wrench.force_N[2] - expected_force) / abs(expected_force))
        for wrench in wrenches
    ]
    force_transverse_relative_errors = [
        float(np.linalg.norm(wrench.force_N[:2]) / abs(expected_force))
        for wrench in wrenches
    ]
    torque_relative_errors = [
        float(np.linalg.norm(wrench.torque_Nm) / (abs(expected_force) * 2.0))
        for wrench in wrenches
    ]
    component_metrics: list[dict[str, float]] = []
    for wrench in wrenches:
        other = wrench.components["other_objects_all_images"].force_N
        external = wrench.components["external_uniform"].force_N
        total = wrench.components["total_external"].force_N
        images = wrench.components["target_periodic_images"].force_N
        primary = np.asarray(
            wrench.numerical_metadata["primary_free_subtraction"]["force_N"],
            dtype=float,
        )
        normalized_components = {
            "other_objects_normalized_absolute": float(
                np.linalg.norm(other) / abs(expected_force)
            ),
            "external_uniform_normalized_absolute": float(
                np.linalg.norm(external) / abs(expected_force)
            ),
            "total_minus_images_normalized_absolute": float(
                np.linalg.norm(total - images) / abs(expected_force)
            ),
            "primary_free_subtraction_normalized_absolute": float(
                np.linalg.norm(primary) / abs(expected_force)
            ),
        }
        normalized_components["component_consistency_relative_error"] = max(
            normalized_components.values()
        )
        component_metrics.append(normalized_components)
    wrench_refinement = [
        {
            "cells_per_axis": cells,
            "force_relative_error": force_error,
            "transverse_relative_error": transverse_error,
            "torque_relative_error": torque_error,
            **components,
        }
        for cells, force_error, transverse_error, torque_error, components in zip(
            wrench_cells,
            force_relative_errors,
            force_transverse_relative_errors,
            torque_relative_errors,
            component_metrics,
        )
    ]
    quadrature_relative_difference = (
        0.0
        if len(wrenches) == 1
        else float(
            np.linalg.norm(wrenches[0].force_N - wrenches[1].force_N)
            / abs(expected_force)
        )
    )
    length = 2.0
    wave_number = 2.0 * math.pi / length
    sigma_amplitude = 2.0 * eps0
    cosine_xy = np.array(
        [
            [0.20, 0.73],
            [0.55, 0.73],
            [1.10, 0.73],
            [1.65, 0.73],
        ]
    )
    sample_abs_z = np.asarray(ORACLE_COSINE_SAMPLE_ABS_Z_M, dtype=float)
    positive_cosine_points = np.concatenate(
        [
            np.column_stack(
                (
                    cosine_xy,
                    np.full(cosine_xy.shape[0], height),
                )
            )
            for height in sample_abs_z
        ]
    )
    cosine_points = positive_cosine_points
    decay = np.exp(-wave_number * np.abs(cosine_points[:, 2]))
    phase = wave_number * cosine_points[:, 0]
    expected_cosine_field = np.zeros_like(cosine_points)
    expected_cosine_field[:, 0] = sigma_amplitude / (2.0 * eps0) * np.sin(phase) * decay
    expected_cosine_field[:, 2] = (
        np.sign(cosine_points[:, 2])
        * sigma_amplitude
        / (2.0 * eps0)
        * np.cos(phase)
        * decay
    )
    expected_cosine_potential = (
        sigma_amplitude / (2.0 * eps0 * wave_number) * np.cos(phase) * decay
    )
    cosine_errors: list[dict[str, float | int]] = []
    for cells_per_axis in (4, 8):
        triangles = plane_triangles_by_cells[cells_per_axis]
        centers = np.mean(triangles, axis=1)
        triangle_area = length**2 / triangles.shape[0]
        charges = sigma_amplitude * np.cos(wave_number * centers[:, 0]) * triangle_area
        with ObjectInteractionSnapshot.from_result(
            _oracle_panel_result(
                root,
                triangles,
                charges,
            ),
            step=None,
            config_path=config_path,
            periodic_model="infinite_physical",
            cache_dir=cache_dir,
            library_path=library_path,
        ) as snapshot:
            record_effective_contract(snapshot)
            field = snapshot._periodic.eval_e(cosine_points)
            potential = snapshot._periodic.eval_phi(cosine_points)
            record_cache_evaluation(f"cosine_{cells_per_axis}", snapshot._periodic)
        positive_field = field[: positive_cosine_points.shape[0]]
        positive_potential = potential[: positive_cosine_points.shape[0]]
        field_amplitudes = np.linalg.norm(
            positive_field.reshape(sample_abs_z.size, -1),
            axis=1,
        )
        potential_amplitudes = np.linalg.norm(
            positive_potential.reshape(sample_abs_z.size, -1),
            axis=1,
        )
        field_decay_ratio = float(field_amplitudes[1] / field_amplitudes[0])
        potential_decay_ratio = float(potential_amplitudes[1] / potential_amplitudes[0])
        cosine_errors.append(
            {
                "cells_per_axis": cells_per_axis,
                "field_relative_error": float(
                    np.linalg.norm(field - expected_cosine_field)
                    / np.linalg.norm(expected_cosine_field)
                ),
                "potential_relative_error": float(
                    np.linalg.norm(potential - expected_cosine_potential)
                    / np.linalg.norm(expected_cosine_potential)
                ),
                "field_decay_ratio": field_decay_ratio,
                "potential_decay_ratio": potential_decay_ratio,
                "field_decay_ratio_relative_error": float(
                    abs(field_decay_ratio - ORACLE_COSINE_EXPECTED_DECAY_RATIO)
                    / ORACLE_COSINE_EXPECTED_DECAY_RATIO
                ),
                "potential_decay_ratio_relative_error": float(
                    abs(potential_decay_ratio - ORACLE_COSINE_EXPECTED_DECAY_RATIO)
                    / ORACLE_COSINE_EXPECTED_DECAY_RATIO
                ),
                "charge_neutrality_ratio": float(
                    abs(np.sum(charges)) / np.sum(np.abs(charges))
                ),
            }
        )

    if set(canonical_cache_identities) != {"4", "8"}:
        raise ValidationError(
            "periodic triangle_p0 plane-oracle did not reuse both "
            "resolution-group operator caches"
        )
    if tuple(row["label"] for row in cache_evaluations) != (
        ORACLE_CACHE_EVALUATION_LABELS
    ):
        raise ValidationError(
            "periodic triangle_p0 plane-oracle cache evaluation coverage is invalid"
        )
    for group, labels in cache_groups.items():
        actual_labels = {
            row["label"] for row in cache_evaluations if row["label"] in labels
        }
        if actual_labels != set(labels):
            raise ValidationError(
                f"periodic triangle_p0 plane-oracle cache group {group} "
                "coverage is invalid"
            )
    coarse_identity = canonical_cache_identities["4"]
    fine_identity = canonical_cache_identities["8"]
    if (
        coarse_identity["fingerprint"] == fine_identity["fingerprint"]
        or coarse_identity["path"] == fine_identity["path"]
    ):
        raise ValidationError(
            "periodic triangle_p0 plane-oracle 4/8 resolution cache "
            "identities must be distinct"
        )
    if observed_effective_field_contract is None:
        raise ValidationError("periodic triangle_p0 plane-oracle has no field contract")
    result: dict[str, Any] = {
        "field_source_model": "triangle_p0",
        "field_kernel_id": "triangle_p0_exact_p2m_near",
        "effective_field_contract": observed_effective_field_contract,
        "cache_diagnostics": dict(coarse_identity),
        "cache_identities": {
            group: dict(identity)
            for group, identity in canonical_cache_identities.items()
        },
        "cache_evaluations": cache_evaluations,
        "uniform_plane": {
            "closure": "e_bottom_zero",
            "cells_per_axis": 4,
            "below_absolute_error_V_m": below_absolute_error,
            "nonzero_relative_error": nonzero_relative_error,
            "transverse_absolute_error_V_m": transverse_absolute_error,
            "potential_gauge_absolute_error_V": potential_gauge_absolute_error,
            "potential_nonzero_relative_error": potential_nonzero_relative_error,
            "force_relative_errors": force_relative_errors,
            "force_transverse_relative_errors": force_transverse_relative_errors,
            "torque_relative_errors": torque_relative_errors,
            "quadrature_relative_difference": quadrature_relative_difference,
            "quadrature_relative_tolerance": ORACLE_QUADRATURE_RELATIVE_TOLERANCE,
            "relative_tolerance": ORACLE_UNIFORM_RELATIVE_TOLERANCE,
            "field_cells_per_axis": 4,
            "target_integration": "gauss_duffy_order_3_and_7",
            "quadrature_order": [3, 7],
            "interpretation": "Maxwell traction under E_bottom=0; not universal free-space self force",
            "wrench_refinement": wrench_refinement,
        },
        "neutral_cosine_plane": {
            "errors": cosine_errors,
            "sample_abs_z_m": list(ORACLE_COSINE_SAMPLE_ABS_Z_M),
            "expected_decay_ratio": ORACLE_COSINE_EXPECTED_DECAY_RATIO,
            "decay_ratio_relative_tolerance": (
                ORACLE_COSINE_DECAY_RATIO_RELATIVE_TOLERANCE
            ),
            "fine_relative_tolerance": ORACLE_COSINE_FINE_RELATIVE_TOLERANCE,
            "expected_decay": "exp(-k*abs(z-z0))",
        },
    }
    _verify_periodic_oracle_metrics(result, label="triangle_p0")
    return result


def probe_periodic_oracles(
    validation_root: str | Path,
    library: str | Path,
) -> dict[str, Any]:
    """Run analytic infinite-periodic plane oracles with the staged library."""

    root = _lexical_absolute(validation_root)
    library_path = _lexical_absolute(library)
    manifest = _load_manifest(root)
    _require_expected_path(
        root,
        str(manifest.get("validation_root", "")),
        None,
        label="validation root",
    )
    staged_library = manifest.get("analysis_library")
    if not isinstance(staged_library, Mapping):
        raise ValidationError("periodic oracles require a staged analysis library")
    staged_library_path = _require_descendant_path(
        root,
        str(staged_library.get("staged_path", "")),
        prefix="lib",
        label="staged analysis library",
    )
    if (
        library_path != staged_library_path
        or not library_path.is_file()
        or _sha256(library_path) != staged_library.get("sha256")
    ):
        raise ValidationError("periodic oracles must use the staged analysis library")
    if (
        manifest.get("build_origin") is not None
        and _library_build_info(library_path) != manifest["build_origin"]
    ):
        raise ValidationError(
            "periodic oracle library build origin differs from the staged source"
        )
    receipt_path = _require_expected_path(
        root,
        root / "provenance/oracles/periodic_plane.json",
        "provenance/oracles/periodic_plane.json",
        label="periodic plane-oracle receipt",
    )
    if receipt_path.exists() or receipt_path.is_symlink():
        return _verify_periodic_oracle_receipt(root, library_path)
    cache_dir = _require_expected_path(
        root,
        root / "cache/oracles",
        "cache/oracles",
        label="periodic plane-oracle cache directory",
    )
    if cache_dir.exists() and any(cache_dir.iterdir()):
        raise ValidationError(
            "periodic oracle cache must be empty before first execution"
        )
    cache_dir.mkdir(parents=True, exist_ok=True)
    config_path = _require_expected_path(
        root,
        root / "provenance/oracles/periodic_plane.toml",
        "provenance/oracles/periodic_plane.toml",
        label="periodic plane-oracle config",
    )
    if any(config_path.parent.iterdir()):
        raise ValidationError(
            "periodic oracle provenance directory must be empty without a receipt"
        )
    _oracle_panel_config(config_path)
    kernel_oracles = {
        "triangle_p0": _run_periodic_plane_kernel_oracle(
            root=root,
            config_path=config_path,
            cache_dir=cache_dir,
            library_path=library_path,
        ),
    }
    report = {
        "receipt_schema_version": 1,
        "oracle_schema_version": 3,
        "status": "qualified",
        "manifest_sha256": _sha256(root / "manifest.json"),
        "library": str(library_path),
        "library_sha256": _sha256(library_path),
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "kernel_configs": {
            "triangle_p0": {
                "path": str(config_path),
                "sha256": _sha256(config_path),
            },
        },
        "cache_dir": str(cache_dir),
        "cache_files": _output_inventory(cache_dir),
        "execution_job_id": os.environ.get("SLURM_JOB_ID", "manual"),
        "kernel_oracles": kernel_oracles,
        "verified_at": _utc_now(),
    }
    if manifest.get("build_origin") is not None:
        report["library_build_origin"] = manifest["build_origin"]
    _verify_periodic_oracle_receipt(
        root,
        library_path,
        _candidate_receipt=report,
    )
    return _publish_execution_receipt(receipt_path, report)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest="command", required=True)
    stage = subcommands.add_parser(
        "stage", help="stage deterministic validation inputs"
    )
    stage.add_argument("--archive-run", required=True, type=Path)
    stage.add_argument("--validation-root", required=True, type=Path)
    stage.add_argument("--binary", required=True, type=Path)
    stage.add_argument("--library", type=Path)
    stage.add_argument("--require-clean-source", action="store_true")
    verify_inputs_parser = subcommands.add_parser(
        "verify-inputs", help="verify staged inputs"
    )
    verify_inputs_parser.add_argument("--validation-root", required=True, type=Path)
    verify_run_parser = subcommands.add_parser(
        "verify-run", help="verify one completed case"
    )
    verify_run_parser.add_argument("--case-dir", required=True, type=Path)
    verify_run_parser.add_argument("--expected-batches", required=True, type=int)
    verify_run_parser.add_argument("--require-existing-receipt", action="store_true")
    verify_run_parser.add_argument(
        "--producer-job-role",
        choices=tuple(PRODUCER_ROLE_CASES),
    )
    analyze = subcommands.add_parser("analyze", help="write comparison artifacts")
    analyze.add_argument("--archive-run", required=True, type=Path)
    analyze.add_argument("--validation-root", required=True, type=Path)
    analyze.add_argument("--library", required=True, type=Path)
    analyze.add_argument("--require-complete", action="store_true")
    probe = subcommands.add_parser(
        "probe-library", help="load and execute a native field-kernel smoke"
    )
    probe.add_argument("--library", required=True, type=Path)
    periodic_oracles = subcommands.add_parser(
        "probe-periodic-oracles",
        help="qualify the staged infinite-periodic kernel against plane oracles",
    )
    periodic_oracles.add_argument("--validation-root", required=True, type=Path)
    periodic_oracles.add_argument("--library", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "stage":
            result = stage_validation(
                args.archive_run,
                args.validation_root,
                args.binary,
                library=args.library,
                require_clean_source=args.require_clean_source,
            )
        elif args.command == "verify-inputs":
            result = verify_inputs(args.validation_root)
        elif args.command == "verify-run":
            result = verify_run(
                args.case_dir,
                args.expected_batches,
                require_existing_receipt=args.require_existing_receipt,
                producer_job_role=args.producer_job_role,
            )
        elif args.command == "analyze":
            result = analyze_validation(
                args.archive_run,
                args.validation_root,
                library=args.library,
                require_complete=args.require_complete,
            )
            physics_status = result.get("physics_evaluation", {}).get("status")
            if physics_status == "partial" or (
                args.require_complete and physics_status != "available"
            ):
                raise ValidationError(
                    "one or more available object physics evaluations failed; "
                    "see analysis_summary.json"
                )
            qualification_status = result.get(
                "numerical_qualification_for_local_frozen_model", {}
            ).get("status")
            if args.require_complete and qualification_status != "qualified":
                raise ValidationError(
                    "local frozen-model numerical qualification failed; "
                    "see analysis_summary.json"
                )
        elif args.command == "probe-library":
            result = probe_library(args.library)
        elif args.command == "probe-periodic-oracles":
            result = probe_periodic_oracles(
                args.validation_root,
                args.library,
            )
        else:  # pragma: no cover
            raise AssertionError(args.command)
    except ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

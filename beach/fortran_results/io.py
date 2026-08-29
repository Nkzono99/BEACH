"""I/O helpers for Fortran output directories."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from beach.summary import CORE_SUMMARY_REQUIRED_KEYS, load_summary_file

from .history import FortranChargeHistory
from .types import (
    ChargeLedgerEntry,
    FieldReconstructionReceipt,
    FortranRunResult,
    MatchingPlaneHistoryEntry,
    MatchingPlaneState,
    MeshSource,
)


_ADAPTIVE_NONZERO_MODE_SUMMARY_KEYS = (
    "simulated_time_s",
    "adaptive_nonzero_mode_rejected_trials",
    "adaptive_nonzero_mode_last_batch_duration_s",
    "adaptive_nonzero_mode_last_potential_step_V",
    "adaptive_nonzero_mode_omp_threads",
)

_MATCHING_PLANE_HISTORY_HEADER = (
    "batch",
    "simulated_time_s",
    "D_H_C_m2",
    "phi_H_V",
    "electron_inward_flux_m2_s",
    "ion_inward_flux_m2_s",
    "electron_access_potential_V",
    "ion_access_potential_V",
    "photoelectron_barrier_potential_V",
    "photoelectron_outward_flux_m2_s",
    "photoelectron_mean_normal_energy_eV",
    "electron_outward_flux_m2_s",
    "ion_outward_flux_m2_s",
    "photoelectron_return_flux_m2_s",
    "photoelectron_escape_flux_m2_s",
    "iterations",
    "residual",
)


def load_fortran_result(directory: str | Path) -> FortranRunResult:
    """Load one Fortran output directory into a structured result object.

    Parameters
    ----------
    directory : str or pathlib.Path
        Output directory that contains ``summary.txt`` and ``charges.csv``.

    Returns
    -------
    FortranRunResult
        Parsed output container with optional mesh/history metadata.

    Raises
    ------
    FileNotFoundError
        If required files are missing.
    ValueError
        If required summary keys are missing or CSV values are malformed.
    """

    out_dir = Path(directory)
    summary = load_summary_file(
        out_dir / "summary.txt",
        required_keys=CORE_SUMMARY_REQUIRED_KEYS,
    )
    mesh_nelem = _parse_nonnegative_int(summary["mesh_nelem"], key="mesh_nelem")
    q_values = _load_charges(out_dir / "charges.csv", mesh_nelem=mesh_nelem)

    triangles, mesh_ids = _load_triangles_if_exists(
        out_dir / "mesh_triangles.csv",
        mesh_nelem=mesh_nelem,
    )
    mesh_sources = _load_mesh_sources_if_exists(out_dir / "mesh_sources.csv")
    mesh_potential_v = _load_mesh_potential_if_exists(
        out_dir / "mesh_potential.csv", mesh_nelem=mesh_nelem
    )
    history_path = out_dir / "charge_history.csv"
    history: FortranChargeHistory | None = None

    if history_path.exists():
        history = FortranChargeHistory(history_path, mesh_nelem=mesh_nelem)
    charge_ledger = _load_charge_ledger_if_exists(out_dir / "charge_ledger.csv")
    _validate_adaptive_nonzero_mode_summary_contract(summary)
    field_reconstruction = _parse_field_reconstruction_receipt(summary)
    matching_plane_state = _parse_matching_plane_state(summary)
    matching_plane_history = None
    if matching_plane_state is not None:
        matching_plane_history = _load_matching_plane_history_if_exists(
            out_dir / "matching_plane_history.csv"
        )
    retry_attempted = _parse_nonnegative_int(
        summary.get("multiple_box_events_retry_attempted", "0"),
        key="multiple_box_events_retry_attempted",
    )
    retry_resolved = _parse_nonnegative_int(
        summary.get("multiple_box_events_retry_resolved", "0"),
        key="multiple_box_events_retry_resolved",
    )
    if retry_resolved > retry_attempted:
        raise ValueError(
            "multiple_box_events_retry_resolved must not exceed "
            "multiple_box_events_retry_attempted"
        )

    return FortranRunResult(
        directory=out_dir,
        mesh_nelem=mesh_nelem,
        processed_particles=_parse_nonnegative_int(
            summary["processed_particles"],
            key="processed_particles",
        ),
        absorbed=_parse_nonnegative_int(summary["absorbed"], key="absorbed"),
        escaped=_parse_nonnegative_int(summary["escaped"], key="escaped"),
        batches=_parse_nonnegative_int(summary["batches"], key="batches"),
        escaped_boundary=_parse_nonnegative_int(
            summary.get("escaped_boundary", "0"),
            key="escaped_boundary",
        ),
        survived_max_step=_parse_nonnegative_int(
            summary.get("survived_max_step", "0"),
            key="survived_max_step",
        ),
        multiple_box_events_retry_attempted=retry_attempted,
        multiple_box_events_retry_resolved=retry_resolved,
        multiple_box_events_soft_discarded=_parse_nonnegative_int(
            summary.get("multiple_box_events_soft_discarded", "0"),
            key="multiple_box_events_soft_discarded",
        ),
        multiple_box_events_soft_discarded_abs_charge_c=_parse_nonnegative_finite_float(
            summary.get("multiple_box_events_soft_discarded_abs_charge_C", "0"),
            key="multiple_box_events_soft_discarded_abs_charge_C",
        ),
        last_rel_change=_parse_nonnegative_finite_float(
            summary["last_rel_change"],
            key="last_rel_change",
        ),
        simulated_time_s=_parse_optional_nonnegative_finite_float(
            summary, "simulated_time_s"
        ),
        adaptive_nonzero_mode_rejected_trials=(
            _parse_optional_nonnegative_int(
                summary, "adaptive_nonzero_mode_rejected_trials"
            )
            or 0
        ),
        adaptive_nonzero_mode_last_batch_duration_s=(
            _parse_optional_nonnegative_finite_float(
                summary, "adaptive_nonzero_mode_last_batch_duration_s"
            )
        ),
        adaptive_nonzero_mode_last_potential_step_v=(
            _parse_optional_nonnegative_finite_float(
                summary, "adaptive_nonzero_mode_last_potential_step_V"
            )
        ),
        adaptive_nonzero_mode_omp_threads=_parse_optional_nonnegative_int(
            summary, "adaptive_nonzero_mode_omp_threads"
        ),
        periodic2_max_nonzero_mode_potential_step_v=(
            _parse_optional_nonnegative_finite_float(
                summary, "periodic2_max_nonzero_mode_potential_step_V"
            )
        ),
        charges=q_values,
        triangles=triangles,
        mesh_ids=mesh_ids,
        mesh_sources=mesh_sources,
        mesh_potential_v=mesh_potential_v,
        history=history,
        checkpoint_schema_version=_parse_optional_nonnegative_int(
            summary, "checkpoint_schema_version"
        ),
        model_fingerprint=summary.get("model_fingerprint"),
        mesh_fingerprint=summary.get("mesh_fingerprint"),
        species_fingerprint=summary.get("species_fingerprint"),
        charge_ledger_residual_c=_parse_optional_finite_float(
            summary, "charge_ledger_residual_C"
        ),
        charge_ledger=charge_ledger,
        field_source_model=summary.get("field_source_model", "unknown").strip().lower(),
        field_kernel_id=summary.get("field_kernel_id"),
        periodic2_cache_hit=_parse_optional_bool(summary, "periodic2_cache_hit"),
        periodic2_operator_build_count=_parse_optional_nonnegative_int(
            summary, "periodic2_operator_build_count"
        ),
        periodic2_cache_fingerprint=summary.get("periodic2_cache_fingerprint"),
        periodic2_cache_path=summary.get("periodic2_cache_path"),
        periodic2_generation_tolerance=_parse_optional_finite_float(
            summary, "periodic2_generation_tolerance"
        ),
        field_reconstruction=field_reconstruction,
        matching_plane_state=matching_plane_state,
        matching_plane_history=matching_plane_history,
    )


def list_fortran_runs(root: str | Path) -> list[Path]:
    """List directories that look like valid Fortran output runs.

    Parameters
    ----------
    root : str or pathlib.Path
        Root directory to scan.

    Returns
    -------
    list[pathlib.Path]
        Directories containing both ``summary.txt`` and ``charges.csv``.
    """

    root_path = Path(root)
    runs: list[Path] = []
    for path in sorted(root_path.glob("*")):
        if not path.is_dir():
            continue
        if (path / "summary.txt").exists() and (path / "charges.csv").exists():
            runs.append(path)
    return runs


def _parse_nonnegative_int(value: str, *, key: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise ValueError(f"summary.txt {key} must be >= 0.")
    return parsed


def _parse_nonnegative_finite_float(value: str, *, key: str) -> float:
    parsed = float(value)
    if not np.isfinite(parsed) or parsed < 0.0:
        raise ValueError(f"summary.txt {key} must be finite and >= 0.")
    return parsed


def _parse_optional_nonnegative_int(data: dict[str, str], key: str) -> int | None:
    if key not in data:
        return None
    return _parse_nonnegative_int(data[key], key=key)


def _parse_optional_finite_float(data: dict[str, str], key: str) -> float | None:
    if key not in data:
        return None
    parsed = float(data[key])
    if not np.isfinite(parsed):
        raise ValueError(f"summary.txt {key} must be finite.")
    return parsed


def _parse_optional_nonnegative_finite_float(
    data: dict[str, str], key: str
) -> float | None:
    if key not in data:
        return None
    return _parse_nonnegative_finite_float(data[key], key=key)


def _parse_optional_bool(data: dict[str, str], key: str) -> bool | None:
    if key not in data:
        return None
    value = data[key].strip().lower()
    if value in {"t", "true"}:
        return True
    if value in {"f", "false"}:
        return False
    raise ValueError(f"summary.txt {key} must be true or false.")


def _parse_matching_plane_state(data: dict[str, str]) -> MatchingPlaneState | None:
    valid = _parse_optional_bool(data, "matching_plane_state_valid")
    if valid is None or not valid:
        return None

    required = (
        "matching_plane_displacement_C_m2",
        "matching_plane_phi_V",
        "matching_plane_electron_inward_flux_m2_s",
        "matching_plane_ion_inward_flux_m2_s",
        "matching_plane_electron_access_potential_V",
        "matching_plane_ion_access_potential_V",
        "matching_plane_photoelectron_barrier_potential_V",
        "matching_plane_photoelectron_outward_flux_m2_s",
        "matching_plane_photoelectron_mean_normal_energy_eV",
        "matching_plane_electron_outward_flux_m2_s",
        "matching_plane_ion_outward_flux_m2_s",
        "matching_plane_photoelectron_return_flux_m2_s",
        "matching_plane_photoelectron_escape_flux_m2_s",
        "matching_plane_iterations",
        "matching_plane_residual",
    )
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(
            "summary.txt valid matching-plane state requires "
            + ", ".join(missing)
            + "."
        )

    state = MatchingPlaneState(
        displacement_c_m2=_parse_optional_finite_float(
            data, "matching_plane_displacement_C_m2"
        ),
        phi_v=_parse_optional_finite_float(data, "matching_plane_phi_V"),
        electron_inward_flux_m2_s=_parse_nonnegative_finite_float(
            data["matching_plane_electron_inward_flux_m2_s"],
            key="matching_plane_electron_inward_flux_m2_s",
        ),
        ion_inward_flux_m2_s=_parse_nonnegative_finite_float(
            data["matching_plane_ion_inward_flux_m2_s"],
            key="matching_plane_ion_inward_flux_m2_s",
        ),
        electron_access_potential_v=_parse_optional_finite_float(
            data, "matching_plane_electron_access_potential_V"
        ),
        ion_access_potential_v=_parse_optional_finite_float(
            data, "matching_plane_ion_access_potential_V"
        ),
        photoelectron_barrier_potential_v=_parse_optional_finite_float(
            data, "matching_plane_photoelectron_barrier_potential_V"
        ),
        photoelectron_outward_flux_m2_s=_parse_nonnegative_finite_float(
            data["matching_plane_photoelectron_outward_flux_m2_s"],
            key="matching_plane_photoelectron_outward_flux_m2_s",
        ),
        photoelectron_mean_normal_energy_ev=_parse_nonnegative_finite_float(
            data["matching_plane_photoelectron_mean_normal_energy_eV"],
            key="matching_plane_photoelectron_mean_normal_energy_eV",
        ),
        electron_outward_flux_m2_s=_parse_nonnegative_finite_float(
            data["matching_plane_electron_outward_flux_m2_s"],
            key="matching_plane_electron_outward_flux_m2_s",
        ),
        ion_outward_flux_m2_s=_parse_nonnegative_finite_float(
            data["matching_plane_ion_outward_flux_m2_s"],
            key="matching_plane_ion_outward_flux_m2_s",
        ),
        photoelectron_return_flux_m2_s=_parse_nonnegative_finite_float(
            data["matching_plane_photoelectron_return_flux_m2_s"],
            key="matching_plane_photoelectron_return_flux_m2_s",
        ),
        photoelectron_escape_flux_m2_s=_parse_nonnegative_finite_float(
            data["matching_plane_photoelectron_escape_flux_m2_s"],
            key="matching_plane_photoelectron_escape_flux_m2_s",
        ),
        iterations=_parse_nonnegative_int(
            data["matching_plane_iterations"], key="matching_plane_iterations"
        ),
        residual=_parse_nonnegative_finite_float(
            data["matching_plane_residual"], key="matching_plane_residual"
        ),
    )
    _validate_matching_plane_state(state, context="summary.txt")
    return state


def _load_matching_plane_history_if_exists(
    path: Path,
) -> tuple[MatchingPlaneHistoryEntry, ...] | None:
    if not path.exists():
        return None

    entries: list[MatchingPlaneHistoryEntry] = []
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != _MATCHING_PLANE_HISTORY_HEADER:
            raise ValueError(
                "matching_plane_history.csv header does not match the 17-column contract."
            )
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(row.get(key) is None for key in _MATCHING_PLANE_HISTORY_HEADER):
                raise ValueError(
                    f"matching_plane_history.csv line {line_number} must have exactly 17 columns."
                )
            state = MatchingPlaneState(
                displacement_c_m2=_parse_history_finite(row, "D_H_C_m2", line_number),
                phi_v=_parse_history_finite(row, "phi_H_V", line_number),
                electron_inward_flux_m2_s=_parse_history_nonnegative(
                    row, "electron_inward_flux_m2_s", line_number
                ),
                ion_inward_flux_m2_s=_parse_history_nonnegative(
                    row, "ion_inward_flux_m2_s", line_number
                ),
                electron_access_potential_v=_parse_history_finite(
                    row, "electron_access_potential_V", line_number
                ),
                ion_access_potential_v=_parse_history_finite(
                    row, "ion_access_potential_V", line_number
                ),
                photoelectron_barrier_potential_v=_parse_history_finite(
                    row, "photoelectron_barrier_potential_V", line_number
                ),
                photoelectron_outward_flux_m2_s=_parse_history_nonnegative(
                    row, "photoelectron_outward_flux_m2_s", line_number
                ),
                photoelectron_mean_normal_energy_ev=_parse_history_nonnegative(
                    row, "photoelectron_mean_normal_energy_eV", line_number
                ),
                electron_outward_flux_m2_s=_parse_history_nonnegative(
                    row, "electron_outward_flux_m2_s", line_number
                ),
                ion_outward_flux_m2_s=_parse_history_nonnegative(
                    row, "ion_outward_flux_m2_s", line_number
                ),
                photoelectron_return_flux_m2_s=_parse_history_nonnegative(
                    row, "photoelectron_return_flux_m2_s", line_number
                ),
                photoelectron_escape_flux_m2_s=_parse_history_nonnegative(
                    row, "photoelectron_escape_flux_m2_s", line_number
                ),
                iterations=_parse_history_nonnegative_int(row, "iterations", line_number),
                residual=_parse_history_nonnegative(row, "residual", line_number),
            )
            _validate_matching_plane_state(
                state, context=f"matching_plane_history.csv line {line_number}"
            )
            batch = _parse_history_nonnegative_int(row, "batch", line_number)
            if batch <= 0 or (entries and batch <= entries[-1].batch):
                raise ValueError(
                    "matching_plane_history.csv batch values must be positive and strictly increasing."
                )
            entries.append(
                MatchingPlaneHistoryEntry(
                    batch=batch,
                    simulated_time_s=_parse_history_nonnegative(
                        row, "simulated_time_s", line_number
                    ),
                    state=state,
                )
            )
    return tuple(entries)


def _parse_history_finite(row: dict[str, str | None], key: str, line: int) -> float:
    raw = row.get(key)
    assert raw is not None
    parsed = float(raw)
    if not np.isfinite(parsed):
        raise ValueError(f"matching_plane_history.csv line {line} {key} must be finite.")
    return parsed


def _parse_history_nonnegative(
    row: dict[str, str | None], key: str, line: int
) -> float:
    parsed = _parse_history_finite(row, key, line)
    if parsed < 0.0:
        raise ValueError(
            f"matching_plane_history.csv line {line} {key} must be >= 0."
        )
    return parsed


def _parse_history_nonnegative_int(
    row: dict[str, str | None], key: str, line: int
) -> int:
    raw = row.get(key)
    assert raw is not None
    try:
        parsed = int(raw)
    except ValueError as exc:
        raise ValueError(
            f"matching_plane_history.csv line {line} {key} must be an integer."
        ) from exc
    if parsed < 0:
        raise ValueError(
            f"matching_plane_history.csv line {line} {key} must be >= 0."
        )
    return parsed


def _validate_matching_plane_state(
    state: MatchingPlaneState, *, context: str
) -> None:
    if state.iterations <= 0:
        raise ValueError(f"{context} matching-plane iterations must be > 0.")
    budget_scale = max(
        1.0,
        state.photoelectron_outward_flux_m2_s,
        state.photoelectron_return_flux_m2_s,
        state.photoelectron_escape_flux_m2_s,
    )
    budget_error = abs(
        state.photoelectron_outward_flux_m2_s
        - state.photoelectron_return_flux_m2_s
        - state.photoelectron_escape_flux_m2_s
    )
    if budget_error > np.sqrt(np.finfo(float).eps) * budget_scale:
        raise ValueError(f"{context} matching-plane photoelectron budget is inconsistent.")


def _parse_optional_fingerprint(data: dict[str, str], key: str) -> str | None:
    if key not in data:
        return None
    value = data[key].strip()
    if len(value) != 16 or any(char not in "0123456789ABCDEF" for char in value):
        raise ValueError(
            f"summary.txt {key} must be exactly 16 uppercase hexadecimal characters."
        )
    return value


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
    return parsed  # type: ignore[return-value]


def _parse_int3(value: str, *, key: str) -> tuple[int, int, int]:
    fields = value.replace(",", " ").split()
    if len(fields) != 3:
        raise ValueError(f"summary.txt {key} must contain exactly three values.")
    parsed = tuple(int(field) for field in fields)
    return parsed  # type: ignore[return-value]


def _validate_adaptive_nonzero_mode_summary_contract(data: dict[str, str]) -> None:
    limit = _parse_optional_nonnegative_finite_float(
        data, "periodic2_max_nonzero_mode_potential_step_V"
    )
    if limit is None or limit == 0.0:
        return
    missing = [key for key in _ADAPTIVE_NONZERO_MODE_SUMMARY_KEYS if key not in data]
    if missing:
        raise ValueError(
            "summary.txt adaptive nonzero-mode output requires "
            + ", ".join(missing)
            + "."
        )
    omp_threads = _parse_nonnegative_int(
        data["adaptive_nonzero_mode_omp_threads"],
        key="adaptive_nonzero_mode_omp_threads",
    )
    if omp_threads == 0:
        raise ValueError(
            "summary.txt adaptive_nonzero_mode_omp_threads must be > 0 "
            "when adaptive nonzero-mode progression is enabled."
        )


def _load_charge_ledger_if_exists(path: Path) -> tuple[ChargeLedgerEntry, ...] | None:
    if not path.exists():
        return None
    charge_fields = (
        "injected_from_remote_C",
        "emitted_from_surface_C",
        "absorbed_on_surface_C",
        "escaped_to_infinity_C",
        "discarded_unresolved_C",
    )
    count_fields = (
        "injected_count",
        "emitted_count",
        "absorbed_count",
        "escaped_count",
        "discarded_unresolved_count",
    )
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"batch", "species_idx", *charge_fields, *count_fields}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("charge_ledger.csv is missing required columns.")
        rows: list[ChargeLedgerEntry] = []
        for row in reader:
            batch = int(row["batch"])
            species_idx = int(row["species_idx"])
            charges = [float(row[name]) for name in charge_fields]
            counts = [int(row[name]) for name in count_fields]
            neutral_return_correction_c = float(
                row.get("neutral_return_correction_C") or 0.0
            )
            neutral_return_weight_scale = float(
                row.get("neutral_return_weight_scale") or 1.0
            )
            neutral_return_unresolved_fraction = float(
                row.get("neutral_return_unresolved_fraction") or 0.0
            )
            fixed_absorbed_target_charge_c = float(
                row.get("fixed_absorbed_target_charge_C") or 0.0
            )
            fixed_absorbed_applied_charge_c = float(
                row.get("fixed_absorbed_applied_charge_C")
                or fixed_absorbed_target_charge_c
            )
            fixed_absorbed_weight_scale = float(
                row.get("fixed_absorbed_weight_scale") or 1.0
            )
            fixed_emission_target_charge_c = float(
                row.get("fixed_emission_target_charge_C") or 0.0
            )
            fixed_emission_applied_charge_c = float(
                row.get("fixed_emission_applied_charge_C")
                or fixed_emission_target_charge_c
            )
            fixed_emission_weight_scale = float(
                row.get("fixed_emission_weight_scale") or 1.0
            )
            fixed_escape_target_charge_c = float(
                row.get("fixed_escape_target_charge_C") or 0.0
            )
            fixed_escape_applied_charge_c = float(
                row.get("fixed_escape_applied_charge_C") or 0.0
            )
            fixed_escape_correction_c = float(
                row.get("fixed_escape_correction_C") or 0.0
            )
            fixed_current_correction_c = float(
                row.get("fixed_current_correction_C") or 0.0
            )
            if batch < 0 or species_idx < 1 or any(count < 0 for count in counts):
                raise ValueError("charge_ledger.csv indices and counts are invalid.")
            neutral_values = (
                neutral_return_correction_c,
                neutral_return_weight_scale,
                neutral_return_unresolved_fraction,
            )
            fixed_values = (
                fixed_absorbed_target_charge_c,
                fixed_absorbed_applied_charge_c,
                fixed_absorbed_weight_scale,
                fixed_emission_target_charge_c,
                fixed_emission_applied_charge_c,
                fixed_emission_weight_scale,
                fixed_escape_target_charge_c,
                fixed_escape_applied_charge_c,
                fixed_escape_correction_c,
                fixed_current_correction_c,
            )
            if not np.all(np.isfinite((*charges, *neutral_values, *fixed_values))):
                raise ValueError("charge_ledger.csv charge values must be finite.")
            rows.append(
                ChargeLedgerEntry(
                    batch,
                    species_idx,
                    *charges,
                    *counts,
                    neutral_return_correction_c=neutral_return_correction_c,
                    neutral_return_weight_scale=neutral_return_weight_scale,
                    neutral_return_unresolved_fraction=neutral_return_unresolved_fraction,
                    fixed_absorbed_target_charge_c=fixed_absorbed_target_charge_c,
                    fixed_absorbed_applied_charge_c=fixed_absorbed_applied_charge_c,
                    fixed_absorbed_weight_scale=fixed_absorbed_weight_scale,
                    fixed_emission_target_charge_c=fixed_emission_target_charge_c,
                    fixed_emission_applied_charge_c=fixed_emission_applied_charge_c,
                    fixed_emission_weight_scale=fixed_emission_weight_scale,
                    fixed_escape_target_charge_c=fixed_escape_target_charge_c,
                    fixed_escape_applied_charge_c=fixed_escape_applied_charge_c,
                    fixed_escape_correction_c=fixed_escape_correction_c,
                    fixed_current_correction_c=fixed_current_correction_c,
                )
            )
    return tuple(rows)


def _load_charges(path: Path, *, mesh_nelem: int) -> np.ndarray:
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.size == 0:
        if mesh_nelem != 0:
            raise ValueError("charges.csv is empty but mesh_nelem > 0.")
        return np.empty(0, dtype=float)
    if data.ndim == 1:
        data = data[None, :]
    if data.shape[1] < 2:
        raise ValueError("charges.csv must contain elem_idx and charge_C columns.")
    if data.shape[0] != mesh_nelem:
        raise ValueError(
            f"charges.csv row count ({data.shape[0]}) does not match mesh_nelem ({mesh_nelem})."
        )
    expected = np.arange(1, mesh_nelem + 1, dtype=np.int64)
    if not np.all(np.isfinite(data[:, 0])) or not np.array_equal(data[:, 0], expected):
        raise ValueError("charges.csv elem_idx column must be 1..mesh_nelem in order.")
    q_values = np.asarray(data[:, 1], dtype=float)
    if not np.all(np.isfinite(q_values)):
        raise ValueError("charges.csv charge_C values must be finite.")
    return q_values


def _load_triangles_if_exists(
    path: Path,
    *,
    mesh_nelem: int,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    if not path.exists():
        return None, None

    with path.open("r", encoding="utf-8") as stream:
        header_line = stream.readline().strip()
    header = [name.strip() for name in header_line.split(",")] if header_line else []
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.size == 0:
        if mesh_nelem != 0:
            raise ValueError("mesh_triangles.csv is empty but mesh_nelem > 0.")
        return None, None
    if data.ndim == 1:
        data = data[None, :]
    if data.shape[1] < 10:
        raise ValueError(
            "mesh_triangles.csv must contain elem_idx and triangle vertices."
        )
    if data.shape[0] != mesh_nelem:
        raise ValueError(
            "mesh_triangles.csv row count "
            f"({data.shape[0]}) does not match mesh_nelem ({mesh_nelem})."
        )
    expected = np.arange(1, mesh_nelem + 1, dtype=np.int64)
    if not np.all(np.isfinite(data[:, 0])) or not np.array_equal(data[:, 0], expected):
        raise ValueError(
            "mesh_triangles.csv elem_idx column must be 1..mesh_nelem in order."
        )
    verts = data[:, 1:10].reshape(-1, 3, 3)
    if not np.all(np.isfinite(verts)):
        raise ValueError("mesh_triangles.csv vertex values must be finite.")
    mesh_ids: np.ndarray | None = None
    if "mesh_id" in header:
        idx = header.index("mesh_id")
        if 0 <= idx < data.shape[1]:
            if not np.all(np.isfinite(data[:, idx])) or not np.array_equal(
                data[:, idx],
                data[:, idx].astype(np.int64),
            ):
                raise ValueError("mesh_triangles.csv mesh_id values must be integer.")
            mesh_ids = data[:, idx].astype(np.int64)
            if np.any(mesh_ids <= 0):
                raise ValueError("mesh_triangles.csv mesh_id values must be positive.")
    return verts, mesh_ids


def _load_mesh_sources_if_exists(path: Path) -> dict[int, MeshSource] | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
    if not rows:
        return None
    out: dict[int, MeshSource] = {}
    for row in rows:
        mesh_id = int(row["mesh_id"])
        if mesh_id <= 0:
            raise ValueError("mesh_sources.csv mesh_id values must be positive.")
        elem_count = int(row.get("elem_count", "0"))
        if elem_count < 0:
            raise ValueError("mesh_sources.csv elem_count values must be >= 0.")
        surface_model = (row.get("surface_model") or "insulator").strip().lower()
        if surface_model not in {"insulator", "conductor", "dielectric"}:
            raise ValueError(
                'mesh_sources.csv surface_model must be "insulator", "conductor", or "dielectric".'
            )
        epsilon_r = float(row.get("epsilon_r") or "1.0")
        if not np.isfinite(epsilon_r) or epsilon_r < 1.0:
            raise ValueError(
                "mesh_sources.csv epsilon_r values must be finite and >= 1."
            )
        out[mesh_id] = MeshSource(
            mesh_id=mesh_id,
            source_kind=row.get("source_kind", ""),
            template_kind=row.get("template_kind", ""),
            elem_count=elem_count,
            surface_model=surface_model,
            epsilon_r=epsilon_r,
        )
    return out


def _load_mesh_potential_if_exists(path: Path, *, mesh_nelem: int) -> np.ndarray | None:
    if not path.exists():
        return None

    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.size == 0:
        if mesh_nelem != 0:
            raise ValueError("mesh_potential.csv is empty but mesh_nelem > 0.")
        return np.empty(0, dtype=float)
    if data.ndim == 1:
        data = data[None, :]
    if data.shape[1] < 2:
        raise ValueError(
            "mesh_potential.csv must contain elem_idx and potential_V columns."
        )
    if data.shape[0] != mesh_nelem:
        raise ValueError(
            f"mesh_potential.csv row count ({data.shape[0]}) does not match mesh_nelem ({mesh_nelem})."
        )
    expected = np.arange(1, mesh_nelem + 1, dtype=np.int64)
    if not np.all(np.isfinite(data[:, 0])) or not np.array_equal(data[:, 0], expected):
        raise ValueError(
            "mesh_potential.csv elem_idx column must be 1..mesh_nelem in order."
        )
    potential = np.asarray(data[:, 1], dtype=float)
    if not np.all(np.isfinite(potential)):
        raise ValueError("mesh_potential.csv potential_V values must be finite.")
    return potential

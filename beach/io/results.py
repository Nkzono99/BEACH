"""Load mesh, charges, histories, and metadata from a Fortran output directory."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .summary import (
    CORE_SUMMARY_REQUIRED_KEYS,
    _parse_summary_fields,
    load_summary_file,
)

from ._field_reconstruction_io import _parse_field_reconstruction_receipt
from ._matching_plane_io import (
    _load_matching_plane_history_if_exists,
    _parse_matching_plane_state,
)
from ._summary_values import (
    _parse_nonnegative_int,
    _parse_optional_nonnegative_finite_float,
)
from ._secondary_emission_io import load_secondary_emission_history
from .csv_history import (
    load_fixed_current_history,
    load_matching_plane_spectrum_history,
    load_top_reference_history,
)
from .performance import load_performance_profile
from .history import BatchHistory
from .types import (
    ChargeLedgerEntry,
    MeshSource,
    RunResult,
    RunSummary,
    RunMesh,
    RunHistory,
    RunReceipts,
    RunDiagnostics,
    RunField,
)

_ADAPTIVE_NONZERO_MODE_SUMMARY_KEYS = (
    "simulated_time_s",
    "adaptive_nonzero_mode_rejected_trials",
    "adaptive_nonzero_mode_last_batch_duration_s",
    "adaptive_nonzero_mode_last_potential_step_V",
    "adaptive_nonzero_mode_omp_threads",
)


def load_run(directory: str | Path) -> RunResult:
    """Read one output directory without constructing a native kernel (light).

    Parameters
    ----------
    directory : str or pathlib.Path
        Output directory containing ``summary.txt`` and ``charges.csv``.

    Returns
    -------
    RunResult
        Counts, geometry, receipts, diagnostics and lazy charge/potential histories.
        Missing optional CSV products are ``None``; header-only diagnostic CSVs
        are empty tuples. Persisted outputs from released versions remain readable.

    Raises
    ------
    FileNotFoundError
        Required files are missing.
    ValueError
        Required summary keys are missing or persisted values are malformed.
    """

    out_dir = Path(directory)
    summary = load_summary_file(
        out_dir / "summary.txt", required_keys=CORE_SUMMARY_REQUIRED_KEYS
    )
    groups = _parse_summary_fields(summary)
    core = RunSummary(**groups["summary"])
    diagnostics = RunDiagnostics(**groups["diagnostics"])
    mesh_nelem = core.mesh_nelem
    charges = _load_charges(out_dir / "charges.csv", mesh_nelem=mesh_nelem)
    triangles, mesh_ids = _load_triangles_if_exists(
        out_dir / "mesh_triangles.csv", mesh_nelem=mesh_nelem
    )
    mesh = RunMesh(
        triangles=triangles,
        mesh_ids=mesh_ids,
        mesh_sources=_load_mesh_sources_if_exists(out_dir / "mesh_sources.csv"),
        mesh_potential_v=_load_mesh_potential_if_exists(
            out_dir / "mesh_potential.csv", mesh_nelem=mesh_nelem
        ),
    )
    histories = {}
    for quantity, column in (("charge", "charge_C"), ("potential", "potential_V")):
        path = out_dir / f"{quantity}_history.csv"
        histories[quantity] = (
            BatchHistory(path, mesh_nelem=mesh_nelem, value_column=column)
            if path.exists()
            else None
        )
    _validate_adaptive_nonzero_mode_summary_contract(summary)
    receipts = RunReceipts(
        **groups["receipts"],
        field_reconstruction=_parse_field_reconstruction_receipt(summary),
        matching_plane_state=_parse_matching_plane_state(summary),
    )
    soft_discarded = diagnostics.multiple_box_events_soft_discarded
    if soft_discarded > core.processed_particles:
        raise ValueError(
            "multiple_box_events_soft_discarded must not exceed processed_particles"
        )
    soft_fraction = (
        soft_discarded / core.processed_particles
        if core.processed_particles > 0
        else 0.0
    )
    if "multiple_box_events_soft_discard_fraction" in summary:
        reported = diagnostics.multiple_box_events_soft_discard_fraction
        if reported > 1.0:
            raise ValueError(
                "summary.txt multiple_box_events_soft_discard_fraction must be <= 1."
            )
        if not np.isclose(
            reported,
            soft_fraction,
            rtol=16.0 * np.finfo(float).eps,
            atol=np.finfo(float).tiny,
        ):
            raise ValueError(
                "summary.txt multiple_box_events_soft_discard_fraction must match "
                "multiple_box_events_soft_discarded / processed_particles."
            )
    if (
        diagnostics.multiple_box_events_retry_resolved
        > diagnostics.multiple_box_events_retry_attempted
    ):
        raise ValueError(
            "multiple_box_events_retry_resolved must not exceed multiple_box_events_retry_attempted"
        )
    from dataclasses import replace

    diagnostics = replace(
        diagnostics, multiple_box_events_soft_discard_fraction=soft_fraction
    )
    profile_path = out_dir / "performance_profile.csv"
    return RunResult(
        directory=out_dir,
        summary=core,
        charges=charges,
        mesh=mesh,
        history=RunHistory(**histories),
        receipts=receipts,
        diagnostics=diagnostics,
        field=RunField(**groups["field"]),
        ledger=_load_charge_ledger_if_exists(out_dir / "charge_ledger.csv"),
        matching_plane_history=(
            _load_matching_plane_history_if_exists(
                out_dir / "matching_plane_history.csv"
            )
            if receipts.matching_plane_state is not None
            else None
        ),
        secondary_emission_history=load_secondary_emission_history(
            out_dir / "secondary_emission_history.csv"
        ),
        top_reference_history=load_top_reference_history(
            out_dir / "top_reference_history.csv"
        ),
        fixed_current_history=load_fixed_current_history(
            out_dir / "fixed_current_history.csv"
        ),
        matching_plane_spectrum_history=load_matching_plane_spectrum_history(
            out_dir / "matching_plane_spectrum_history.csv"
        ),
        performance_profile=load_performance_profile(profile_path)
        if profile_path.exists()
        else None,
    )


def list_runs(root: str | Path) -> list[Path]:
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

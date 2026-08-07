"""I/O helpers for Fortran output directories."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from beach.summary import CORE_SUMMARY_REQUIRED_KEYS, load_summary_file

from .history import FortranChargeHistory
from .types import ChargeLedgerEntry, FortranRunResult, MeshSource


_ADAPTIVE_NONZERO_MODE_SUMMARY_KEYS = (
    "simulated_time_s",
    "adaptive_nonzero_mode_rejected_trials",
    "adaptive_nonzero_mode_last_batch_duration_s",
    "adaptive_nonzero_mode_last_potential_step_V",
    "adaptive_nonzero_mode_omp_threads",
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


def _parse_optional_fingerprint(data: dict[str, str], key: str) -> str | None:
    if key not in data:
        return None
    value = data[key].strip()
    if len(value) != 16 or any(char not in "0123456789ABCDEF" for char in value):
        raise ValueError(
            f"summary.txt {key} must be exactly 16 uppercase hexadecimal characters."
        )
    return value


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

"""I/O helpers for Fortran output directories."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

import numpy as np

from .history import FortranChargeHistory
from .types import FortranRunResult, MeshSource


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
    summary = _load_summary(out_dir / "summary.txt")
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
        last_rel_change=_parse_nonnegative_finite_float(
            summary["last_rel_change"],
            key="last_rel_change",
        ),
        charges=q_values,
        triangles=triangles,
        mesh_ids=mesh_ids,
        mesh_sources=mesh_sources,
        mesh_potential_v=mesh_potential_v,
        history=history,
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


def _load_summary(path: Path) -> dict[str, str]:
    data: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        data[key.strip()] = value.strip()
    _ensure_keys(
        data,
        [
            "mesh_nelem",
            "processed_particles",
            "absorbed",
            "escaped",
            "batches",
            "last_rel_change",
        ],
    )
    return data


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
        raise ValueError("mesh_triangles.csv must contain elem_idx and triangle vertices.")
    if data.shape[0] != mesh_nelem:
        raise ValueError(
            "mesh_triangles.csv row count "
            f"({data.shape[0]}) does not match mesh_nelem ({mesh_nelem})."
        )
    expected = np.arange(1, mesh_nelem + 1, dtype=np.int64)
    if not np.all(np.isfinite(data[:, 0])) or not np.array_equal(data[:, 0], expected):
        raise ValueError("mesh_triangles.csv elem_idx column must be 1..mesh_nelem in order.")
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
            raise ValueError("mesh_sources.csv epsilon_r values must be finite and >= 1.")
        out[mesh_id] = MeshSource(
            mesh_id=mesh_id,
            source_kind=row.get("source_kind", ""),
            template_kind=row.get("template_kind", ""),
            elem_count=elem_count,
            surface_model=surface_model,
            epsilon_r=epsilon_r,
        )
    return out


def _load_mesh_potential_if_exists(
    path: Path, *, mesh_nelem: int
) -> np.ndarray | None:
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
        raise ValueError("mesh_potential.csv must contain elem_idx and potential_V columns.")
    if data.shape[0] != mesh_nelem:
        raise ValueError(
            f"mesh_potential.csv row count ({data.shape[0]}) does not match mesh_nelem ({mesh_nelem})."
        )
    expected = np.arange(1, mesh_nelem + 1, dtype=np.int64)
    if not np.all(np.isfinite(data[:, 0])) or not np.array_equal(data[:, 0], expected):
        raise ValueError("mesh_potential.csv elem_idx column must be 1..mesh_nelem in order.")
    potential = np.asarray(data[:, 1], dtype=float)
    if not np.all(np.isfinite(potential)):
        raise ValueError("mesh_potential.csv potential_V values must be finite.")
    return potential


def _ensure_keys(data: dict[str, str], required: Iterable[str]) -> None:
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(f"Missing keys in summary: {missing}")

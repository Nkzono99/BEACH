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
    mesh_nelem = int(summary["mesh_nelem"])
    charges = np.loadtxt(out_dir / "charges.csv", delimiter=",", skiprows=1)
    if charges.ndim == 1:
        charges = charges[None, :]
    q_values = charges[:, 1]

    triangles, mesh_ids = _load_triangles_if_exists(out_dir / "mesh_triangles.csv")
    mesh_sources = _load_mesh_sources_if_exists(out_dir / "mesh_sources.csv")
    history_path = out_dir / "charge_history.csv"
    history: FortranChargeHistory | None = None

    if history_path.exists():
        history = FortranChargeHistory(history_path, mesh_nelem=mesh_nelem)

    return FortranRunResult(
        directory=out_dir,
        mesh_nelem=mesh_nelem,
        processed_particles=int(summary["processed_particles"]),
        absorbed=int(summary["absorbed"]),
        escaped=int(summary["escaped"]),
        batches=int(summary["batches"]),
        escaped_boundary=int(summary.get("escaped_boundary", 0)),
        survived_max_step=int(summary.get("survived_max_step", 0)),
        last_rel_change=float(summary["last_rel_change"]),
        charges=q_values,
        triangles=triangles,
        mesh_ids=mesh_ids,
        mesh_sources=mesh_sources,
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


def _load_triangles_if_exists(path: Path) -> tuple[np.ndarray | None, np.ndarray | None]:
    if not path.exists():
        return None, None

    with path.open("r", encoding="utf-8") as stream:
        header_line = stream.readline().strip()
    header = [name.strip() for name in header_line.split(",")] if header_line else []
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.ndim == 1:
        data = data[None, :]
    verts = data[:, 1:10].reshape(-1, 3, 3)
    mesh_ids: np.ndarray | None = None
    if "mesh_id" in header:
        idx = header.index("mesh_id")
        if 0 <= idx < data.shape[1]:
            mesh_ids = data[:, idx].astype(np.int64)
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
        out[mesh_id] = MeshSource(
            mesh_id=mesh_id,
            source_kind=row.get("source_kind", ""),
            template_kind=row.get("template_kind", ""),
            elem_count=int(row.get("elem_count", "0")),
        )
    return out


def _ensure_keys(data: dict[str, str], required: Iterable[str]) -> None:
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(f"Missing keys in summary: {missing}")

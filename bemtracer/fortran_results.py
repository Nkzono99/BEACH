"""Utilities for reading and visualizing Fortran simulation outputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class FortranRunResult:
    """Container for one Fortran simulation output directory."""

    directory: Path
    mesh_nelem: int
    processed_particles: int
    absorbed: int
    escaped: int
    batches: int
    last_rel_change: float
    charges: np.ndarray


def load_fortran_result(directory: str | Path) -> FortranRunResult:
    """Load summary and per-element charges written by the Fortran executable."""

    out_dir = Path(directory)
    summary = _load_summary(out_dir / "summary.txt")
    charges = np.loadtxt(out_dir / "charges.csv", delimiter=",", skiprows=1)
    if charges.ndim == 1:
        charges = charges[None, :]
    q_values = charges[:, 1]

    return FortranRunResult(
        directory=out_dir,
        mesh_nelem=int(summary["mesh_nelem"]),
        processed_particles=int(summary["processed_particles"]),
        absorbed=int(summary["absorbed"]),
        escaped=int(summary["escaped"]),
        batches=int(summary["batches"]),
        last_rel_change=float(summary["last_rel_change"]),
        charges=q_values,
    )


def list_fortran_runs(root: str | Path) -> list[Path]:
    """List directories that look like Fortran output runs under ``root``."""

    root_path = Path(root)
    runs: list[Path] = []
    for path in sorted(root_path.glob("*")):
        if not path.is_dir():
            continue
        if (path / "summary.txt").exists() and (path / "charges.csv").exists():
            runs.append(path)
    return runs


def plot_charges(result: FortranRunResult):
    """Quick bar plot of per-element charge values."""

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 3))
    x = np.arange(result.charges.size)
    ax.bar(x, result.charges)
    ax.set_xlabel("element index")
    ax.set_ylabel("charge [C]")
    ax.set_title(f"Fortran charge distribution: {result.directory}")
    fig.tight_layout()
    return fig, ax


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


def _ensure_keys(data: dict[str, str], required: Iterable[str]) -> None:
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(f"Missing keys in summary: {missing}")

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
    escaped_boundary: int
    survived_max_step: int
    last_rel_change: float
    charges: np.ndarray
    triangles: np.ndarray | None = None
    charge_history: np.ndarray | None = None
    processed_particles_by_batch: np.ndarray | None = None
    rel_change_by_batch: np.ndarray | None = None
    batch_indices: np.ndarray | None = None


def load_fortran_result(directory: str | Path) -> FortranRunResult:
    """Load summary and per-element outputs written by the Fortran executable."""

    out_dir = Path(directory)
    summary = _load_summary(out_dir / "summary.txt")
    charges = np.loadtxt(out_dir / "charges.csv", delimiter=",", skiprows=1)
    if charges.ndim == 1:
        charges = charges[None, :]
    q_values = charges[:, 1]

    triangles = _load_triangles_if_exists(out_dir / "mesh_triangles.csv")
    charge_history, processed_particles_by_batch, rel_change_by_batch, batch_indices = _load_charge_history_if_exists(
        out_dir / "charge_history.csv", mesh_nelem=int(summary["mesh_nelem"])
    )

    return FortranRunResult(
        directory=out_dir,
        mesh_nelem=int(summary["mesh_nelem"]),
        processed_particles=int(summary["processed_particles"]),
        absorbed=int(summary["absorbed"]),
        escaped=int(summary["escaped"]),
        batches=int(summary["batches"]),
        escaped_boundary=int(summary.get("escaped_boundary", 0)),
        survived_max_step=int(summary.get("survived_max_step", 0)),
        last_rel_change=float(summary["last_rel_change"]),
        charges=q_values,
        triangles=triangles,
        charge_history=charge_history,
        processed_particles_by_batch=processed_particles_by_batch,
        rel_change_by_batch=rel_change_by_batch,
        batch_indices=batch_indices,
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


def plot_charge_mesh(result: FortranRunResult, *, cmap: str = "coolwarm"):
    """Plot mesh triangles in 3D, colored by per-element surface charge density."""

    if result.triangles is None:
        raise ValueError("mesh_triangles.csv is not found. Re-run Fortran with latest output format.")

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")

    sigma = _surface_charge_density(result.charges, result.triangles)
    sigma_max_abs = float(np.max(np.abs(sigma)))
    if sigma_max_abs <= np.finfo(float).tiny:
        sigma_max_abs = 1.0

    norm = plt.Normalize(vmin=-sigma_max_abs, vmax=sigma_max_abs)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    facecolors = sm.to_rgba(sigma)

    mesh = Poly3DCollection(
        result.triangles,
        facecolors=facecolors,
        edgecolor=(0.0, 0.0, 0.0, 0.45),
        linewidth=0.35,
        alpha=0.88,
    )
    ax.add_collection3d(mesh)

    pts = result.triangles.reshape(-1, 3)
    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)
    center = (mins + maxs) * 0.5
    radius = float(np.max(maxs - mins)) * 0.5
    radius = max(radius, 1.0e-12)

    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)
    ax.set_box_aspect((1.0, 1.0, 1.0))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.view_init(elev=24.0, azim=-58.0)
    ax.set_title(f"Surface charge density mesh: {result.directory}")

    sm.set_array(sigma)
    fig._beach_color_mappable = sm
    fig.colorbar(sm, ax=ax, shrink=0.75, label="surface charge density [C/m^2]")
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



def _load_charge_history_if_exists(
    path: Path, *, mesh_nelem: int
) -> tuple[np.ndarray | None, np.ndarray | None, np.ndarray | None, np.ndarray | None]:
    if not path.exists():
        return None, None, None, None
    if path.stat().st_size == 0:
        return None, None, None, None
    with path.open("r", encoding="utf-8") as stream:
        stream.readline()
        for line in stream:
            if line.strip():
                break
        else:
            return None, None, None, None

    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.ndim == 1:
        data = data[None, :]

    batches = data[:, 0].astype(np.int64)
    processed = data[:, 1].astype(np.int64)
    rel_change = data[:, 2]
    elem_idx = data[:, 3].astype(np.int64)
    charges = data[:, 4]

    batch_indices = np.unique(batches)
    n_snapshots = batch_indices.size
    history = np.zeros((mesh_nelem, n_snapshots), dtype=float)
    processed_by_batch = np.zeros(n_snapshots, dtype=np.int64)
    rel_by_batch = np.zeros(n_snapshots, dtype=float)

    for col, batch in enumerate(batch_indices):
        mask = batches == batch
        batch_elem_idx = elem_idx[mask] - 1
        history[batch_elem_idx, col] = charges[mask]
        processed_by_batch[col] = processed[mask][0]
        rel_by_batch[col] = rel_change[mask][0]

    return history, processed_by_batch, rel_by_batch, batch_indices


def _load_triangles_if_exists(path: Path) -> np.ndarray | None:
    if not path.exists():
        return None

    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.ndim == 1:
        data = data[None, :]
    verts = data[:, 1:10].reshape(-1, 3, 3)
    return verts


def _surface_charge_density(charges: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    areas = _triangle_areas(triangles)
    min_area = np.finfo(float).tiny
    return charges / np.maximum(areas, min_area)


def _triangle_areas(triangles: np.ndarray) -> np.ndarray:
    edge_a = triangles[:, 1, :] - triangles[:, 0, :]
    edge_b = triangles[:, 2, :] - triangles[:, 0, :]
    cross = np.cross(edge_a, edge_b)
    return 0.5 * np.linalg.norm(cross, axis=1)


def _ensure_keys(data: dict[str, str], required: Iterable[str]) -> None:
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(f"Missing keys in summary: {missing}")

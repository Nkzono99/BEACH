"""Utilities for reading and visualizing Fortran simulation outputs."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from matplotlib.animation import FuncAnimation

K_COULOMB = 8.9875517923e9


@dataclass(frozen=True)
class MeshSource:
    """Metadata for one source mesh written by Fortran output."""

    mesh_id: int
    source_kind: str
    template_kind: str
    elem_count: int


@dataclass(frozen=True)
class MeshSelection:
    """One mesh (or merged meshes) extracted from a run result."""

    directory: Path
    mesh_ids: tuple[int, ...]
    elem_indices: np.ndarray
    triangles: np.ndarray
    charges: np.ndarray
    step: int | None = None


@dataclass(frozen=True)
class CoulombInteraction:
    """Coulomb interaction summary between two mesh groups."""

    group_a_mesh_ids: tuple[int, ...]
    group_b_mesh_ids: tuple[int, ...]
    step: int | None
    softening: float
    torque_origin_m: np.ndarray
    force_on_a_N: np.ndarray
    force_on_b_N: np.ndarray
    torque_on_a_Nm: np.ndarray
    torque_on_b_Nm: np.ndarray
    mean_force_on_a_per_element_N: np.ndarray
    mean_torque_on_a_per_element_Nm: np.ndarray


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
    mesh_ids: np.ndarray | None = None
    mesh_sources: dict[int, MeshSource] | None = None
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

    triangles, mesh_ids = _load_triangles_if_exists(out_dir / "mesh_triangles.csv")
    mesh_sources = _load_mesh_sources_if_exists(out_dir / "mesh_sources.csv")
    charge_history, processed_particles_by_batch, rel_change_by_batch, batch_indices = (
        _load_charge_history_if_exists(
            out_dir / "charge_history.csv", mesh_nelem=int(summary["mesh_nelem"])
        )
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
        mesh_ids=mesh_ids,
        mesh_sources=mesh_sources,
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


class Beach:
    """Facade for one Fortran output directory with lazy result loading."""

    def __init__(self, output_dir: str | Path = "outputs/latest") -> None:
        self.output_dir = Path(output_dir)
        self._result: FortranRunResult | None = None

    @property
    def result(self) -> FortranRunResult:
        if self._result is None:
            self._result = load_fortran_result(self.output_dir)
        return self._result

    def reload(self) -> FortranRunResult:
        self._result = load_fortran_result(self.output_dir)
        return self._result

    @property
    def mesh_ids(self) -> tuple[int, ...]:
        ids = np.unique(_mesh_ids_or_default(self.result))
        return tuple(int(v) for v in ids)

    def get_mesh(self, *mesh_ids: int, step: int | None = -1):
        if len(mesh_ids) == 0:
            raise ValueError("at least one mesh id must be provided.")
        selections = tuple(
            _build_mesh_selection(self.result, (int(mesh_id),), step=step)
            for mesh_id in mesh_ids
        )
        if len(selections) == 1:
            return selections[0]
        return selections

    def get_mesh_charge(self, *mesh_ids: int, step: int | None = -1):
        selection = self.get_mesh(*mesh_ids, step=step)
        if isinstance(selection, tuple):
            return tuple(mesh.charges.copy() for mesh in selection)
        return selection.charges.copy()

    def calc_coulomb(
        self,
        group_a: int | MeshSelection | Iterable[int | MeshSelection],
        group_b: int | MeshSelection | Iterable[int | MeshSelection],
        *,
        step: int | None = -1,
        softening: float = 0.0,
        torque_origin: Literal["group_a_center", "group_b_center", "origin"] = (
            "group_a_center"
        ),
    ) -> CoulombInteraction:
        return calc_coulomb(
            self.result,
            group_a,
            group_b,
            step=step,
            softening=softening,
            torque_origin=torque_origin,
        )

    def plot_mesh(self, *, cmap: str = "coolwarm"):
        return plot_charge_mesh(self.result, cmap=cmap)

    def plot_bar(self):
        return plot_charges(self.result)

    def compute_potential(
        self,
        *,
        softening: float = 0.0,
        self_term: str = "area_equivalent",
    ) -> np.ndarray:
        return compute_potential_mesh(
            self.result,
            softening=softening,
            self_term=self_term,
        )

    def plot_potential(
        self,
        *,
        softening: float = 0.0,
        self_term: str = "area_equivalent",
        cmap: str = "viridis",
    ):
        return plot_potential_mesh(
            self.result,
            softening=softening,
            self_term=self_term,
            cmap=cmap,
        )

    def animate_mesh(
        self,
        output_path: str | Path | None = None,
        *,
        quantity: Literal["charge", "potential"] = "charge",
        fps: int = 10,
        frame_stride: int = 1,
        total_frames: int | None = None,
        cmap: str | None = None,
        softening: float = 0.0,
        self_term: str = "area_equivalent",
    ) -> Path | FuncAnimation:
        return animate_history_mesh(
            self.result,
            output_path=output_path,
            quantity=quantity,
            fps=fps,
            frame_stride=frame_stride,
            total_frames=total_frames,
            cmap=cmap,
            softening=softening,
            self_term=self_term,
        )


def _resolve_result(result: FortranRunResult | Beach) -> FortranRunResult:
    if isinstance(result, Beach):
        return result.result
    return result


def calc_coulomb(
    result: FortranRunResult | Beach,
    group_a: int | MeshSelection | Iterable[int | MeshSelection],
    group_b: int | MeshSelection | Iterable[int | MeshSelection],
    *,
    step: int | None = -1,
    softening: float = 0.0,
    torque_origin: Literal["group_a_center", "group_b_center", "origin"] = (
        "group_a_center"
    ),
) -> CoulombInteraction:
    """Compute Coulomb force/torque between two mesh groups."""

    resolved = _resolve_result(result)
    if softening < 0.0:
        raise ValueError("softening must be >= 0.")

    sel_a = _coerce_group_selection(resolved, group_a, step=step)
    sel_b = _coerce_group_selection(resolved, group_b, step=step)
    if sel_a.elem_indices.size == 0:
        raise ValueError("group_a does not contain any mesh elements.")
    if sel_b.elem_indices.size == 0:
        raise ValueError("group_b does not contain any mesh elements.")

    if torque_origin == "group_a_center":
        origin = _triangle_centers(sel_a.triangles).mean(axis=0)
    elif torque_origin == "group_b_center":
        origin = _triangle_centers(sel_b.triangles).mean(axis=0)
    elif torque_origin == "origin":
        origin = np.zeros(3, dtype=float)
    else:
        raise ValueError(
            "torque_origin must be one of {'group_a_center', 'group_b_center', 'origin'}."
        )

    centers_a = _triangle_centers(sel_a.triangles)
    centers_b = _triangle_centers(sel_b.triangles)
    force_a, torque_a = _pairwise_force_torque(
        centers_a,
        sel_a.charges,
        centers_b,
        sel_b.charges,
        origin=origin,
        softening=softening,
    )
    force_b = -force_a
    torque_b = -torque_a

    return CoulombInteraction(
        group_a_mesh_ids=sel_a.mesh_ids,
        group_b_mesh_ids=sel_b.mesh_ids,
        step=sel_a.step,
        softening=softening,
        torque_origin_m=origin,
        force_on_a_N=force_a,
        force_on_b_N=force_b,
        torque_on_a_Nm=torque_a,
        torque_on_b_Nm=torque_b,
        mean_force_on_a_per_element_N=force_a / float(sel_a.elem_indices.size),
        mean_torque_on_a_per_element_Nm=torque_a / float(sel_a.elem_indices.size),
    )


def plot_charges(result: FortranRunResult | Beach):
    """Quick bar plot of per-element charge values."""

    import matplotlib.pyplot as plt

    resolved = _resolve_result(result)
    fig, ax = plt.subplots(figsize=(8, 3))
    x = np.arange(resolved.charges.size)
    ax.bar(x, resolved.charges)
    ax.set_xlabel("element index")
    ax.set_ylabel("charge [C]")
    ax.set_title(f"Fortran charge distribution: {resolved.directory}")
    fig.tight_layout()
    return fig, ax


def plot_charge_mesh(result: FortranRunResult | Beach, *, cmap: str = "coolwarm"):
    """Plot mesh triangles in 3D, colored by per-element surface charge density."""

    resolved = _resolve_result(result)
    triangles = _require_triangles(resolved)
    sigma = _surface_charge_density(resolved.charges, triangles)
    return _plot_scalar_mesh(
        triangles,
        sigma,
        title=f"Surface charge density mesh: {resolved.directory}",
        colorbar_label="surface charge density [C/m^2]",
        cmap=cmap,
    )


def compute_potential_mesh(
    result: FortranRunResult | Beach,
    *,
    softening: float = 0.0,
    self_term: str = "area_equivalent",
) -> np.ndarray:
    """Compute one electric potential value per triangle centroid."""

    resolved = _resolve_result(result)
    if softening < 0.0:
        raise ValueError("softening must be >= 0.")

    triangles = _require_triangles(resolved)
    centers = _triangle_centers(triangles)
    offdiag_kernel = _potential_offdiag_kernel(centers, softening=softening)
    self_coeff = _self_potential_coefficients(
        triangles, self_term=self_term, softening=softening
    )
    potential = offdiag_kernel @ resolved.charges + self_coeff * resolved.charges
    return K_COULOMB * potential


def plot_potential_mesh(
    result: FortranRunResult | Beach,
    *,
    softening: float = 0.0,
    self_term: str = "area_equivalent",
    cmap: str = "viridis",
):
    """Plot mesh triangles in 3D, colored by per-element electric potential."""

    resolved = _resolve_result(result)
    triangles = _require_triangles(resolved)
    phi = compute_potential_mesh(resolved, softening=softening, self_term=self_term)
    return _plot_scalar_mesh(
        triangles,
        phi,
        title=f"Electric potential mesh: {resolved.directory}",
        colorbar_label="potential [V]",
        cmap=cmap,
    )


def animate_history_mesh(
    result: FortranRunResult | Beach,
    output_path: str | Path | None = None,
    *,
    quantity: Literal["charge", "potential"] = "charge",
    fps: int = 10,
    frame_stride: int = 1,
    total_frames: int | None = None,
    cmap: str | None = None,
    softening: float = 0.0,
    self_term: str = "area_equivalent",
) -> Path | FuncAnimation:
    """Render charge-history snapshots as a GIF animation on the 3D mesh."""

    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    resolved = _resolve_result(result)
    if fps <= 0:
        raise ValueError("fps must be > 0.")
    if frame_stride <= 0:
        raise ValueError("frame_stride must be > 0.")
    if total_frames is not None and total_frames <= 0:
        raise ValueError("total_frames must be > 0.")
    if total_frames is not None and frame_stride != 1:
        raise ValueError("frame_stride and total_frames cannot be used together.")

    triangles = _require_triangles(resolved)
    charge_history = _require_charge_history(resolved)
    quantity_key = quantity.lower()
    self_term_key = self_term.replace("-", "_")
    frame_cols = _select_frame_columns(
        charge_history.shape[1],
        frame_stride=frame_stride,
        total_frames=total_frames,
    )
    charges_sampled = charge_history[:, frame_cols]

    if quantity_key == "charge":
        values_history = _surface_charge_density_history(charges_sampled, triangles)
        colorbar_label = "surface charge density [C/m^2]"
        title_prefix = "Surface charge density history"
        use_cmap = "coolwarm" if cmap is None else cmap
    elif quantity_key == "potential":
        values_history = _potential_history(
            charges_sampled,
            triangles,
            softening=softening,
            self_term=self_term_key,
        )
        colorbar_label = "potential [V]"
        title_prefix = "Electric potential history"
        use_cmap = "viridis" if cmap is None else cmap
    else:
        raise ValueError("quantity must be one of {'charge', 'potential'}.")

    max_abs = float(np.max(np.abs(values_history)))
    if max_abs <= np.finfo(float).tiny:
        max_abs = 1.0
    norm = plt.Normalize(vmin=-max_abs, vmax=max_abs)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=use_cmap)

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    mesh = Poly3DCollection(
        triangles,
        facecolors=sm.to_rgba(values_history[:, 0]),
        edgecolor=(0.0, 0.0, 0.0, 0.45),
        linewidth=0.35,
        alpha=0.88,
    )
    ax.add_collection3d(mesh)
    _configure_mesh_axes(ax, triangles)

    def _title_for_frame(frame_idx: int) -> str:
        col = int(frame_cols[frame_idx])
        suffix: list[str] = [f"frame={frame_idx + 1}/{frame_cols.size}"]
        if resolved.batch_indices is not None and col < resolved.batch_indices.size:
            suffix.append(f"batch={int(resolved.batch_indices[col])}")
        if (
            resolved.processed_particles_by_batch is not None
            and col < resolved.processed_particles_by_batch.size
        ):
            suffix.append(
                f"processed={int(resolved.processed_particles_by_batch[col])}"
            )
        if (
            resolved.rel_change_by_batch is not None
            and col < resolved.rel_change_by_batch.size
        ):
            suffix.append(f"rel={resolved.rel_change_by_batch[col]:.3e}")
        return f"{title_prefix}: {resolved.directory}\n" + " ".join(suffix)

    ax.set_title(_title_for_frame(0))
    sm.set_array(values_history[:, 0])
    fig.colorbar(sm, ax=ax, shrink=0.75, label=colorbar_label)
    fig.tight_layout()

    def _update(frame_idx: int):
        mesh.set_facecolor(sm.to_rgba(values_history[:, frame_idx]))
        ax.set_title(_title_for_frame(frame_idx))
        return (mesh,)

    animation = FuncAnimation(
        fig,
        _update,
        frames=frame_cols.size,
        interval=1000.0 / float(fps),
        blit=False,
    )

    if output_path is None:
        return animation

    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    animation.save(out_path, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return out_path


def _select_frame_columns(
    n_snapshots: int,
    *,
    frame_stride: int,
    total_frames: int | None,
) -> np.ndarray:
    if n_snapshots <= 0:
        raise ValueError("n_snapshots must be > 0.")

    if total_frames is None:
        return np.arange(0, n_snapshots, frame_stride, dtype=np.int64)

    if total_frames >= n_snapshots:
        return np.arange(n_snapshots, dtype=np.int64)

    if total_frames == 1:
        return np.array([0], dtype=np.int64)

    numerators = np.arange(total_frames, dtype=np.int64) * (n_snapshots - 1)
    return numerators // (total_frames - 1)


def _coerce_group_selection(
    result: FortranRunResult,
    group: int | MeshSelection | Iterable[int | MeshSelection],
    *,
    step: int | None,
) -> MeshSelection:
    mesh_ids: list[int] = []
    mesh_steps: list[int] = []

    def _add_item(item: int | MeshSelection) -> None:
        if isinstance(item, MeshSelection):
            if item.directory != result.directory:
                raise ValueError("mesh selections must come from the same output directory.")
            mesh_ids.extend(item.mesh_ids)
            if item.step is not None:
                mesh_steps.append(item.step)
            return
        mesh_ids.append(int(item))

    if isinstance(group, MeshSelection) or isinstance(group, (int, np.integer)):
        _add_item(group)
    else:
        for item in group:
            _add_item(item)

    if not mesh_ids:
        raise ValueError("mesh group is empty.")

    requested_step = step
    if requested_step is None and mesh_steps:
        unique_steps = sorted(set(mesh_steps))
        if len(unique_steps) > 1:
            raise ValueError("mesh selections in one group must share the same step.")
        requested_step = unique_steps[0]

    return _build_mesh_selection(result, tuple(mesh_ids), step=requested_step)


def _build_mesh_selection(
    result: FortranRunResult, mesh_ids: tuple[int, ...], *, step: int | None
) -> MeshSelection:
    unique_mesh_ids = tuple(dict.fromkeys(int(mid) for mid in mesh_ids))
    all_mesh_ids = _mesh_ids_or_default(result)
    available = set(int(mid) for mid in np.unique(all_mesh_ids))
    missing = [mid for mid in unique_mesh_ids if mid not in available]
    if missing:
        raise ValueError(f"unknown mesh id(s): {missing}. available={sorted(available)}")

    triangles = _require_triangles(result)
    charges = _charges_for_step(result, step=step)
    mask = np.isin(all_mesh_ids, np.asarray(unique_mesh_ids, dtype=np.int64))
    elem_indices = np.flatnonzero(mask)

    return MeshSelection(
        directory=result.directory,
        mesh_ids=unique_mesh_ids,
        elem_indices=elem_indices,
        triangles=triangles[elem_indices],
        charges=charges[elem_indices],
        step=step,
    )


def _charges_for_step(result: FortranRunResult, *, step: int | None) -> np.ndarray:
    if step is None:
        return result.charges

    history = result.charge_history
    batch_indices = result.batch_indices
    if step == -1:
        if history is None or history.size == 0:
            return result.charges
        return history[:, -1]

    if history is None or history.size == 0:
        raise ValueError(
            "charge_history.csv is required when step is specified and must not be empty."
        )
    if batch_indices is None:
        raise ValueError("batch indices are unavailable although charge history exists.")
    cols = np.flatnonzero(batch_indices == step)
    if cols.size == 0:
        available = [int(v) for v in batch_indices]
        raise ValueError(f"step={step} is not found in history. available={available}")
    return history[:, int(cols[0])]


def _mesh_ids_or_default(result: FortranRunResult) -> np.ndarray:
    if result.mesh_ids is None or result.mesh_ids.size != result.mesh_nelem:
        return np.ones(result.mesh_nelem, dtype=np.int64)
    return result.mesh_ids.astype(np.int64, copy=False)


def _pairwise_force_torque(
    centers_a: np.ndarray,
    charges_a: np.ndarray,
    centers_b: np.ndarray,
    charges_b: np.ndarray,
    *,
    origin: np.ndarray,
    softening: float,
) -> tuple[np.ndarray, np.ndarray]:
    force = np.zeros(3, dtype=float)
    torque = np.zeros(3, dtype=float)
    eps2 = softening * softening
    min_dist2 = np.finfo(float).tiny

    for i in range(centers_a.shape[0]):
        delta = centers_a[i] - centers_b
        dist2 = np.sum(delta * delta, axis=1) + eps2
        inv_r3 = 1.0 / (np.maximum(dist2, min_dist2) * np.sqrt(np.maximum(dist2, min_dist2)))
        coeff = K_COULOMB * charges_a[i] * charges_b * inv_r3
        f_i = np.sum(coeff[:, None] * delta, axis=0)
        force += f_i
        torque += np.cross(centers_a[i] - origin, f_i)

    return force, torque


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


def _surface_charge_density(charges: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    areas = _triangle_areas(triangles)
    min_area = np.finfo(float).tiny
    return charges / np.maximum(areas, min_area)


def _surface_charge_density_history(
    charges_history: np.ndarray, triangles: np.ndarray
) -> np.ndarray:
    areas = _triangle_areas(triangles)
    min_area = np.finfo(float).tiny
    return charges_history / np.maximum(areas[:, None], min_area)


def _triangle_centers(triangles: np.ndarray) -> np.ndarray:
    return triangles.mean(axis=1)


def _triangle_areas(triangles: np.ndarray) -> np.ndarray:
    edge_a = triangles[:, 1, :] - triangles[:, 0, :]
    edge_b = triangles[:, 2, :] - triangles[:, 0, :]
    cross = np.cross(edge_a, edge_b)
    return 0.5 * np.linalg.norm(cross, axis=1)


def _potential_offdiag_kernel(centers: np.ndarray, *, softening: float) -> np.ndarray:
    delta = centers[:, None, :] - centers[None, :, :]
    dist2 = np.sum(delta * delta, axis=2) + softening * softening
    min_dist2 = np.finfo(float).tiny
    inv_r = 1.0 / np.sqrt(np.maximum(dist2, min_dist2))
    np.fill_diagonal(inv_r, 0.0)
    return inv_r


def _potential_history(
    charges_history: np.ndarray,
    triangles: np.ndarray,
    *,
    softening: float,
    self_term: str,
) -> np.ndarray:
    if softening < 0.0:
        raise ValueError("softening must be >= 0.")

    centers = _triangle_centers(triangles)
    offdiag_kernel = _potential_offdiag_kernel(centers, softening=softening)
    self_coeff = _self_potential_coefficients(
        triangles, self_term=self_term, softening=softening
    )
    potential = offdiag_kernel @ charges_history + self_coeff[:, None] * charges_history
    return K_COULOMB * potential


def _self_potential_coefficients(
    triangles: np.ndarray, *, self_term: str, softening: float
) -> np.ndarray:
    if self_term == "exclude":
        return np.zeros(triangles.shape[0], dtype=float)

    if self_term == "softened_point":
        if softening <= 0.0:
            raise ValueError("softening must be > 0 when self_term='softened_point'.")
        return np.full(triangles.shape[0], 1.0 / softening, dtype=float)

    if self_term == "area_equivalent":
        areas = _triangle_areas(triangles)
        min_area = np.finfo(float).tiny
        safe_areas = np.maximum(areas, min_area)
        return 2.0 * np.sqrt(np.pi) / np.sqrt(safe_areas)

    raise ValueError(
        "self_term must be one of {'area_equivalent', 'exclude', 'softened_point'}."
    )


def _require_triangles(result: FortranRunResult) -> np.ndarray:
    if result.triangles is None:
        raise ValueError(
            "mesh_triangles.csv is not found. Re-run Fortran with latest output format."
        )
    return result.triangles


def _require_charge_history(result: FortranRunResult) -> np.ndarray:
    if result.charge_history is None or result.charge_history.size == 0:
        raise ValueError(
            "charge_history.csv is not found or empty. Enable history output and rerun."
        )
    return result.charge_history


def _configure_mesh_axes(ax, triangles: np.ndarray) -> None:
    pts = triangles.reshape(-1, 3)
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


def _plot_scalar_mesh(
    triangles: np.ndarray,
    values: np.ndarray,
    *,
    title: str,
    colorbar_label: str,
    cmap: str,
):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")

    max_abs = float(np.max(np.abs(values)))
    if max_abs <= np.finfo(float).tiny:
        max_abs = 1.0

    norm = plt.Normalize(vmin=-max_abs, vmax=max_abs)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    facecolors = sm.to_rgba(values)

    mesh = Poly3DCollection(
        triangles,
        facecolors=facecolors,
        edgecolor=(0.0, 0.0, 0.0, 0.45),
        linewidth=0.35,
        alpha=0.88,
    )
    ax.add_collection3d(mesh)
    _configure_mesh_axes(ax, triangles)
    ax.set_title(title)

    sm.set_array(values)
    fig._beach_color_mappable = sm
    fig.colorbar(sm, ax=ax, shrink=0.75, label=colorbar_label)
    fig.tight_layout()
    return fig, ax


def _ensure_keys(data: dict[str, str], required: Iterable[str]) -> None:
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(f"Missing keys in summary: {missing}")

"""Mesh animation helpers for history outputs."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, TYPE_CHECKING

import numpy as np

from .mesh import _configure_mesh_axes, _surface_charge_density_history
from .potential import _potential_history
from .selection import _require_history, _require_triangles, _resolve_result
from .types import FortranRunResult

if TYPE_CHECKING:
    from matplotlib.animation import FuncAnimation


def animate_history_mesh(
    result: FortranRunResult | object,
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
    """Render mesh history as an animation.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.
    output_path : str or pathlib.Path or None, default None
        Output GIF path. If ``None``, return a live ``FuncAnimation`` object.
    quantity : {"charge", "potential"}, default "charge"
        Quantity used for per-face coloring.
    fps : int, default 10
        Frames per second.
    frame_stride : int, default 1
        Use every ``frame_stride``-th snapshot.
    total_frames : int or None, default None
        Evenly sampled frame count. Cannot be combined with ``frame_stride != 1``.
    cmap : str or None, default None
        Matplotlib colormap name. ``None`` uses quantity-specific defaults.
    softening : float, default 0.0
        Softening length in meters (used for potential mode).
    self_term : {"area_equivalent", "exclude", "softened_point"}, default "area_equivalent"
        Potential self-term model (used for potential mode).

    Returns
    -------
    pathlib.Path or matplotlib.animation.FuncAnimation
        Written GIF path when ``output_path`` is provided, otherwise animation object.

    Raises
    ------
    ValueError
        If arguments are invalid or history/geometry data is unavailable.
    """

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
    history = _require_history(resolved)
    charge_history = history.as_array()
    batch_indices = history.batch_indices
    processed_by_batch = history.processed_particles_by_batch
    rel_change_by_batch = history.rel_change_by_batch
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
        if batch_indices is not None and col < batch_indices.size:
            suffix.append(f"batch={int(batch_indices[col])}")
        if processed_by_batch is not None and col < processed_by_batch.size:
            suffix.append(f"processed={int(processed_by_batch[col])}")
        if rel_change_by_batch is not None and col < rel_change_by_batch.size:
            suffix.append(f"rel={rel_change_by_batch[col]:.3e}")
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

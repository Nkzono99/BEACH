"""Mesh animation helpers for history outputs."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Literal, Mapping, TYPE_CHECKING

import numpy as np

from beach.analysis.context import RunContext
from beach.analysis.mesh import (
    _maybe_apply_periodic2_mesh,
    _maybe_replicate_periodic2,
    _surface_charge_density_history,
)
from beach.analysis.potential import _run_potential_history
from beach.io.types import RunResult

from ._defaults import (
    DEFAULT_CHARGE_CMAP,
    DEFAULT_MESH_VIEW_AZIM,
    DEFAULT_MESH_VIEW_ELEV,
    DEFAULT_POTENTIAL_CMAP,
)
from ._mesh3d import _configure_mesh_axes

if TYPE_CHECKING:
    from matplotlib.animation import FuncAnimation

    from beach.analysis.field_kernel import FieldKernel


def animate_history_mesh(
    result: RunResult | object,
    output_path: str | Path | None = None,
    *,
    quantity: Literal["charge", "potential"] = "charge",
    fps: int = 10,
    frame_stride: int = 1,
    total_frames: int | None = None,
    cmap: str | None = None,
    periodic2: Mapping[str, object] | None = None,
    apply_periodic2_mesh: bool = False,
    periodic2_repeat: int = 0,
    reference_point: Iterable[float] | str | None = None,
    view_elev: float | None = None,
    view_azim: float | None = None,
    library_path: str | Path | None = None,
    kernel: FieldKernel | None = None,
) -> Path | FuncAnimation:
    """Render mesh history as an animation (heavy for ``quantity="potential"``).

    Parameters
    ----------
    result : RunResult or Beach-like object
        Run result or object exposing ``result`` as ``RunResult``.
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
    periodic2 : mapping or None, default None
        Two-axis periodic setting for potential mode. ``None`` の場合は
        出力ディレクトリ近傍の ``beach.toml`` から自動判定する。
    apply_periodic2_mesh : bool, default False
        ``True`` の場合、triangle 重心を周期セルへ wrap する平行移動を各 face に
        適用して描画する。
    periodic2_repeat : int, default 0
        周期イメージの複製数。``n`` を指定すると各周期軸方向に ``-n`` から ``+n``
        までの ``(2n+1)^2`` コピーを生成して描画する。``0`` は複製なし。
    reference_point : iterable of float, {"species1_injection_center"}, or None, default None
        基準電位を差し引く参照点。``None`` は絶対電位、``"species1_injection_center"``
        は species1 の注入面中心を基準にした電位差を描く。
    view_elev : float or None, default None
        3D view elevation angle in degrees. ``None`` uses the default (24.0).
    view_azim : float or None, default None
        3D view azimuth angle in degrees. ``None`` uses the default (-58.0).
    library_path : str, pathlib.Path, or None, default None
        Native kernel library for potential mode.
    kernel : FieldKernel or None, default None
        Open kernel built from this run that evaluates potential mode with its own
        plan; its source charges are stepped through the history frames and then
        restored, and it is not closed. ``periodic2`` then only wraps and
        replicates the drawn mesh. ``None`` builds a kernel for potential mode.

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

    context = RunContext.from_value(result)
    resolved = context.result
    if fps <= 0:
        raise ValueError("fps must be > 0.")
    if frame_stride <= 0:
        raise ValueError("frame_stride must be > 0.")
    if total_frames is not None and total_frames <= 0:
        raise ValueError("total_frames must be > 0.")
    if total_frames is not None and frame_stride != 1:
        raise ValueError("frame_stride and total_frames cannot be used together.")

    triangles = context.require_triangles()
    plot_triangles = _maybe_apply_periodic2_mesh(
        resolved,
        triangles,
        periodic2=periodic2,
        apply_periodic2_mesh=apply_periodic2_mesh,
        context=context,
    )
    history = context.require_history()
    batch_indices = history.batch_indices
    processed_by_batch = history.processed_particles_by_batch
    rel_change_by_batch = history.rel_change_by_batch
    quantity_key = quantity.lower()
    frame_cols = _select_frame_columns(
        batch_indices.size,
        frame_stride=frame_stride,
        total_frames=total_frames,
    )
    # Read only the selected snapshots instead of the whole (nelem, nsnapshot) history.
    charges_sampled = history._select_columns(frame_cols)

    if quantity_key == "charge":
        values_history = _surface_charge_density_history(charges_sampled, triangles)
        colorbar_label = "surface charge density [C/m^2]"
        title_prefix = "Surface charge density history"
        use_cmap = DEFAULT_CHARGE_CMAP if cmap is None else cmap
    elif quantity_key == "potential":
        values_history, resolved_reference = _run_potential_history(
            context,
            charges_sampled,
            periodic2=periodic2,
            reference_point=reference_point,
            library_path=library_path,
            kernel=kernel,
        )
        colorbar_label = "potential [V]" if resolved_reference is None else "potential difference [V]"
        title_prefix = "Electric potential history" if resolved_reference is None else "Electric potential difference history"
        use_cmap = DEFAULT_POTENTIAL_CMAP if cmap is None else cmap
    else:
        raise ValueError("quantity must be one of {'charge', 'potential'}.")

    plot_triangles, values_history = _maybe_replicate_periodic2(
        resolved,
        plot_triangles,
        values_history,
        periodic2=periodic2,
        periodic2_repeat=periodic2_repeat,
        context=context,
    )

    max_abs = float(np.max(np.abs(values_history)))
    if max_abs <= np.finfo(float).tiny:
        max_abs = 1.0
    norm = plt.Normalize(vmin=-max_abs, vmax=max_abs)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=use_cmap)

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    mesh = Poly3DCollection(
        plot_triangles,
        facecolors=sm.to_rgba(values_history[:, 0]),
        edgecolor=(0.0, 0.0, 0.0, 0.45),
        linewidth=0.35,
        alpha=0.88,
    )
    ax.add_collection3d(mesh)
    _configure_mesh_axes(
        ax,
        plot_triangles,
        view_elev=DEFAULT_MESH_VIEW_ELEV if view_elev is None else view_elev,
        view_azim=DEFAULT_MESH_VIEW_AZIM if view_azim is None else view_azim,
    )

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

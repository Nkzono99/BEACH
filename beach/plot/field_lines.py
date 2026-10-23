"""3D field-line figures over the source mesh."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Mapping

import numpy as np

from beach.analysis.context import RunContext
from beach.analysis.field_lines import trace_field_lines
from beach.analysis.mesh import _surface_charge_density
from beach.io.types import FortranRunResult

from ._defaults import (
    DEFAULT_CHARGE_CMAP,
    DEFAULT_FIELD_LINE_CMAP,
    DEFAULT_MESH_VIEW_AZIM,
    DEFAULT_MESH_VIEW_ELEV,
)
from ._mesh3d import _configure_mesh_axes

if TYPE_CHECKING:
    from beach.analysis.field_kernel import FieldKernel


def plot_field_lines_3d(
    result: FortranRunResult | object,
    seed_points: np.ndarray,
    *,
    ds: float | None = None,
    max_steps: int = 500,
    periodic2: Mapping[str, object] | None = None,
    direction: str = "both",
    box_min: Iterable[float] | None = None,
    box_max: Iterable[float] | None = None,
    show_mesh: bool = True,
    mesh_alpha: float = 0.25,
    mesh_cmap: str = DEFAULT_CHARGE_CMAP,
    line_color: str | None = None,
    line_cmap: str = DEFAULT_FIELD_LINE_CMAP,
    line_width: float = 1.2,
    view_elev: float = DEFAULT_MESH_VIEW_ELEV,
    view_azim: float = DEFAULT_MESH_VIEW_AZIM,
    title: str = "Electric field lines",
    figsize: tuple[float, float] = (9, 7),
    config_path: str | Path | None = None,
    library_path: str | Path | None = None,
    kernel: FieldKernel | None = None,
):
    """Plot triangle-P0 field lines and an optional surface mesh (heavy: native kernel).

    The tracing arguments, including ``kernel``, are those of
    :func:`beach.trace_field_lines`. ``show_mesh`` overlays the mesh colored by
    surface charge density with ``mesh_cmap`` and ``mesh_alpha``; lines use
    ``line_color`` or, when it is ``None``, ``line_cmap``. Returns
    ``(figure, axes)``.
    """

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    context = RunContext.from_value(result, config_path=config_path)
    resolved = context.result
    lines = trace_field_lines(
        context,
        seed_points,
        ds=ds,
        max_steps=max_steps,
        periodic2=periodic2,
        direction=direction,
        box_min=box_min,
        box_max=box_max,
        config_path=config_path,
        library_path=library_path,
        kernel=kernel,
    )
    triangles = context.require_triangles()
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")
    if show_mesh:
        density = _surface_charge_density(resolved.charges, triangles)
        max_abs = max(float(np.max(np.abs(density))), np.finfo(float).tiny)
        norm = plt.Normalize(vmin=-max_abs, vmax=max_abs)
        mapper = plt.cm.ScalarMappable(norm=norm, cmap=mesh_cmap)
        facecolors = mapper.to_rgba(density)
        facecolors[:, 3] = mesh_alpha
        ax.add_collection3d(
            Poly3DCollection(
                triangles,
                facecolors=facecolors,
                edgecolor=(0.0, 0.0, 0.0, 0.15),
                linewidth=0.2,
            )
        )

    if line_color is None:
        cmap = plt.get_cmap(line_cmap)
        colors = [
            cmap(index / max(len(lines) - 1, 1))
            for index in range(len(lines))
        ]
    else:
        colors = [line_color] * len(lines)
    for line, color in zip(lines, colors):
        if line.shape[0] < 2:
            continue
        ax.plot(
            line[:, 0],
            line[:, 1],
            line[:, 2],
            color=color,
            linewidth=line_width,
        )
        midpoint = line.shape[0] // 2
        if 0 < midpoint < line.shape[0] - 1:
            delta = line[midpoint + 1] - line[midpoint]
            ax.quiver(
                *line[midpoint],
                *delta,
                color=color,
                arrow_length_ratio=0.4,
                linewidth=line_width * 1.5,
            )

    seeds = np.asarray(seed_points, dtype=float)
    if seeds.ndim == 1:
        seeds = seeds.reshape(1, 3)
    ax.scatter(
        seeds[:, 0],
        seeds[:, 1],
        seeds[:, 2],
        c="red",
        s=20,
        zorder=5,
        depthshade=False,
    )
    _configure_mesh_axes(
        ax,
        triangles,
        view_elev=view_elev,
        view_azim=view_azim,
    )
    ax.set_title(title)
    fig.tight_layout()
    return fig, ax

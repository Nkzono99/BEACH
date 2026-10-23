"""3D mesh axes and scalar-mesh rendering shared by the plotting modules."""

from __future__ import annotations

import numpy as np

from ._defaults import DEFAULT_MESH_VIEW_AZIM, DEFAULT_MESH_VIEW_ELEV


def _configure_mesh_axes(
    ax,
    triangles: np.ndarray,
    *,
    view_elev: float = DEFAULT_MESH_VIEW_ELEV,
    view_azim: float = DEFAULT_MESH_VIEW_AZIM,
    axis_unit: str = "m",
) -> None:
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
    unit, _scale = _resolve_axis_unit(axis_unit)
    ax.set_xlabel(f"x [{unit}]")
    ax.set_ylabel(f"y [{unit}]")
    ax.set_zlabel(f"z [{unit}]")
    ax.view_init(elev=view_elev, azim=view_azim)


def _resolve_axis_unit(axis_unit: str) -> tuple[str, float]:
    unit = str(axis_unit).strip().lower()
    if unit == "m":
        return "m", 1.0
    if unit in {"um", "micron", "microns", "micrometer", "micrometers"}:
        return "um", 1.0e6
    if unit in {"nm", "nanometer", "nanometers"}:
        return "nm", 1.0e9
    raise ValueError("axis_unit must be one of {'m', 'um', 'nm'}.")


def _plot_scalar_mesh(
    triangles: np.ndarray,
    values: np.ndarray,
    *,
    title: str,
    colorbar_label: str,
    cmap: str,
    view_elev: float = DEFAULT_MESH_VIEW_ELEV,
    view_azim: float = DEFAULT_MESH_VIEW_AZIM,
    axis_unit: str = "m",
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

    _unit, axis_scale = _resolve_axis_unit(axis_unit)
    plot_triangles = triangles * axis_scale

    mesh = Poly3DCollection(
        plot_triangles,
        facecolors=facecolors,
        edgecolor=(0.0, 0.0, 0.0, 0.45),
        linewidth=0.35,
        alpha=0.88,
    )
    ax.add_collection3d(mesh)
    _configure_mesh_axes(
        ax,
        plot_triangles,
        view_elev=view_elev,
        view_azim=view_azim,
        axis_unit=axis_unit,
    )
    ax.set_title(title)

    sm.set_array(values)
    fig._beach_color_mappable = sm
    fig.colorbar(sm, ax=ax, shrink=0.75, label=colorbar_label)
    fig.tight_layout()
    return fig, ax

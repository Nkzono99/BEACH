"""Geometry and mesh scalar-field helpers."""

from __future__ import annotations

import numpy as np


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

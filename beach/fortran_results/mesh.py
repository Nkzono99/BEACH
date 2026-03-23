"""Geometry and mesh scalar-field helpers."""

from __future__ import annotations

from typing import TYPE_CHECKING, Mapping

import numpy as np

if TYPE_CHECKING:
    from .types import FortranRunResult

DEFAULT_MESH_VIEW_ELEV = 24.0
DEFAULT_MESH_VIEW_AZIM = -58.0


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


def _wrap_periodic2_triangles_by_centroid(
    triangles: np.ndarray,
    *,
    axes: tuple[int, int],
    lengths: tuple[float, float],
    origins: tuple[float, float],
) -> np.ndarray:
    shifted = np.asarray(triangles, dtype=float).copy()
    centers = _triangle_centers(shifted)
    shift = np.zeros((shifted.shape[0], 3), dtype=float)

    for axis, length, origin in zip(axes, lengths, origins):
        wrapped_centers = origin + np.mod(centers[:, axis] - origin, length)
        shift[:, axis] = wrapped_centers - centers[:, axis]

    shifted += shift[:, None, :]
    return shifted


def _maybe_apply_periodic2_mesh(
    resolved: "FortranRunResult",
    triangles: np.ndarray,
    *,
    periodic2: Mapping[str, object] | None = None,
    apply_periodic2_mesh: bool = False,
) -> np.ndarray:
    if not apply_periodic2_mesh:
        return triangles

    from .potential import _auto_periodic2_from_result, _coerce_periodic2

    periodic_cfg = _coerce_periodic2(periodic2)
    if periodic_cfg is None:
        periodic_cfg = _auto_periodic2_from_result(resolved)
    if periodic_cfg is None:
        raise ValueError(
            "apply_periodic2_mesh requires periodic2 settings or nearby beach.toml "
            'with sim.field_bc_mode="periodic2".'
        )

    axes, lengths, origins, _nimg, _far_correction, _alpha, _ewald_layers = periodic_cfg
    return _wrap_periodic2_triangles_by_centroid(
        triangles,
        axes=axes,
        lengths=lengths,
        origins=origins,
    )


def _configure_mesh_axes(
    ax,
    triangles: np.ndarray,
    *,
    view_elev: float = DEFAULT_MESH_VIEW_ELEV,
    view_azim: float = DEFAULT_MESH_VIEW_AZIM,
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
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.view_init(elev=view_elev, azim=view_azim)


def _plot_scalar_mesh(
    triangles: np.ndarray,
    values: np.ndarray,
    *,
    title: str,
    colorbar_label: str,
    cmap: str,
    view_elev: float = DEFAULT_MESH_VIEW_ELEV,
    view_azim: float = DEFAULT_MESH_VIEW_AZIM,
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
    _configure_mesh_axes(
        ax,
        triangles,
        view_elev=view_elev,
        view_azim=view_azim,
    )
    ax.set_title(title)

    sm.set_array(values)
    fig._beach_color_mappable = sm
    fig.colorbar(sm, ax=ax, shrink=0.75, label=colorbar_label)
    fig.tight_layout()
    return fig, ax

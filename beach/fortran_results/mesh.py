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


def _wrap_periodic2_triangles_by_mesh_centroid(
    triangles: np.ndarray,
    mesh_ids: np.ndarray,
    *,
    axes: tuple[int, int],
    lengths: tuple[float, float],
    origins: tuple[float, float],
) -> np.ndarray:
    shifted = np.asarray(triangles, dtype=float).copy()
    mesh_ids_array = np.asarray(mesh_ids, dtype=np.int64)
    if mesh_ids_array.shape != (shifted.shape[0],):
        raise ValueError("mesh_ids must have one value per triangle.")

    for mesh_id in np.unique(mesh_ids_array):
        mask = mesh_ids_array == mesh_id
        for axis, length, origin in zip(axes, lengths, origins):
            coords = shifted[mask, :, axis]
            center, concentration = _periodic_coordinate_mean(
                coords,
                length=length,
                origin=origin,
            )
            if concentration < 1.0e-3:
                continue
            shifted[mask, :, axis] = center + np.mod(
                coords - center + 0.5 * length,
                length,
            ) - 0.5 * length

        mesh_center = _triangle_centers(shifted[mask]).mean(axis=0)
        shift = np.zeros(3, dtype=float)
        for axis, length, origin in zip(axes, lengths, origins):
            wrapped_center = origin + np.mod(mesh_center[axis] - origin, length)
            shift[axis] = wrapped_center - mesh_center[axis]
        shifted[mask] += shift
    return shifted


def _periodic_coordinate_mean(
    coords: np.ndarray,
    *,
    length: float,
    origin: float,
) -> tuple[float, float]:
    angles = 2.0 * np.pi * (np.asarray(coords, dtype=float) - origin) / length
    mean_sin = float(np.mean(np.sin(angles)))
    mean_cos = float(np.mean(np.cos(angles)))
    concentration = float(np.hypot(mean_cos, mean_sin))
    mean_angle = float(np.arctan2(mean_sin, mean_cos))
    if mean_angle < 0.0:
        mean_angle += 2.0 * np.pi
    return origin + mean_angle * length / (2.0 * np.pi), concentration


def _maybe_apply_periodic2_mesh(
    resolved: "FortranRunResult",
    triangles: np.ndarray,
    *,
    periodic2: Mapping[str, object] | None = None,
    apply_periodic2_mesh: bool = False,
    periodic2_mesh_mode: str = "centroid",
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
    mode = str(periodic2_mesh_mode).strip().lower()
    if mode in {"centroid", "face", "triangle"}:
        return _wrap_periodic2_triangles_by_centroid(
            triangles,
            axes=axes,
            lengths=lengths,
            origins=origins,
        )
    if mode in {"mesh", "object"}:
        if resolved.mesh_ids is None:
            raise ValueError("periodic2_mesh_mode='mesh' requires mesh_ids.")
        return _wrap_periodic2_triangles_by_mesh_centroid(
            triangles,
            resolved.mesh_ids,
            axes=axes,
            lengths=lengths,
            origins=origins,
        )
    raise ValueError("periodic2_mesh_mode must be one of {'centroid', 'mesh'}.")


def _replicate_periodic2(
    triangles: np.ndarray,
    values: np.ndarray,
    *,
    axes: tuple[int, int],
    lengths: tuple[float, float],
    n_repeat: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Tile triangles and per-element values over periodic images.

    Parameters
    ----------
    triangles : numpy.ndarray
        Triangle vertices with shape ``(nelem, 3, 3)``.
    values : numpy.ndarray
        Per-element scalar values with shape ``(nelem,)`` or
        ``(nelem, n_snapshots)`` for history data.
    axes : tuple of int
        Two periodic axis indices.
    lengths : tuple of float
        Box lengths along each periodic axis.
    n_repeat : int
        Number of periodic images in each direction.
        ``n_repeat=1`` produces offsets ``-1, 0, +1`` → 3×3 = 9 copies.

    Returns
    -------
    tuple of numpy.ndarray
        ``(replicated_triangles, replicated_values)``.
    """
    if n_repeat <= 0:
        return triangles, values

    all_triangles: list[np.ndarray] = []
    all_values: list[np.ndarray] = []
    for ix in range(-n_repeat, n_repeat + 1):
        for iy in range(-n_repeat, n_repeat + 1):
            shifted = triangles.copy()
            shifted[:, :, axes[0]] += ix * lengths[0]
            shifted[:, :, axes[1]] += iy * lengths[1]
            all_triangles.append(shifted)
            all_values.append(values)

    return np.concatenate(all_triangles, axis=0), np.concatenate(all_values, axis=0)


def _maybe_replicate_periodic2(
    resolved: "FortranRunResult",
    triangles: np.ndarray,
    values: np.ndarray,
    *,
    periodic2: Mapping[str, object] | None = None,
    periodic2_repeat: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Replicate triangles and values for periodic display when requested.

    Returns the original arrays unchanged when ``periodic2_repeat <= 0``.
    """
    if periodic2_repeat <= 0:
        return triangles, values

    from .potential import _auto_periodic2_from_result, _coerce_periodic2

    periodic_cfg = _coerce_periodic2(periodic2)
    if periodic_cfg is None:
        periodic_cfg = _auto_periodic2_from_result(resolved)
    if periodic_cfg is None:
        raise ValueError(
            "periodic2_repeat requires periodic2 settings or nearby beach.toml "
            'with sim.field_bc_mode="periodic2".'
        )

    axes, lengths, _origins, _nimg, _far_correction, _alpha, _ewald_layers = periodic_cfg
    return _replicate_periodic2(
        triangles,
        values,
        axes=axes,
        lengths=lengths,
        n_repeat=periodic2_repeat,
    )


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

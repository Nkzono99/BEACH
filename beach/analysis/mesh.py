"""Triangle geometry, surface densities, and periodic2 mesh wrapping."""

from __future__ import annotations

from typing import Mapping

import numpy as np

from beach.io.types import FortranRunResult
from beach.native.periodic2 import coerce_periodic2

from .context import RunContext
from .periodic import auto_periodic2_from_result


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
    resolved: FortranRunResult,
    triangles: np.ndarray,
    *,
    periodic2: Mapping[str, object] | None = None,
    apply_periodic2_mesh: bool = False,
    periodic2_mesh_mode: str = "centroid",
    context: RunContext | None = None,
) -> np.ndarray:
    if not apply_periodic2_mesh:
        return triangles

    periodic_cfg = coerce_periodic2(periodic2)
    if periodic_cfg is None:
        periodic_cfg = auto_periodic2_from_result(resolved, context=context)
    if periodic_cfg is None:
        raise ValueError(
            "apply_periodic2_mesh requires periodic2 settings or nearby beach.toml "
            'with sim.field_bc_mode="periodic2".'
        )

    mode = str(periodic2_mesh_mode).strip().lower()
    if mode in {"centroid", "face", "triangle"}:
        return _wrap_periodic2_triangles_by_centroid(
            triangles,
            axes=periodic_cfg.axes,
            lengths=periodic_cfg.lengths,
            origins=periodic_cfg.origins,
        )
    if mode in {"mesh", "object"}:
        if resolved.mesh_ids is None:
            raise ValueError("periodic2_mesh_mode='mesh' requires mesh_ids.")
        return _wrap_periodic2_triangles_by_mesh_centroid(
            triangles,
            resolved.mesh_ids,
            axes=periodic_cfg.axes,
            lengths=periodic_cfg.lengths,
            origins=periodic_cfg.origins,
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
    resolved: FortranRunResult,
    triangles: np.ndarray,
    values: np.ndarray,
    *,
    periodic2: Mapping[str, object] | None = None,
    periodic2_repeat: int = 0,
    context: RunContext | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Replicate triangles and values for periodic display when requested.

    Returns the original arrays unchanged when ``periodic2_repeat <= 0``.
    """
    if periodic2_repeat <= 0:
        return triangles, values

    periodic_cfg = coerce_periodic2(periodic2)
    if periodic_cfg is None:
        periodic_cfg = auto_periodic2_from_result(resolved, context=context)
    if periodic_cfg is None:
        raise ValueError(
            "periodic2_repeat requires periodic2 settings or nearby beach.toml "
            'with sim.field_bc_mode="periodic2".'
        )

    return _replicate_periodic2(
        triangles,
        values,
        axes=periodic_cfg.axes,
        lengths=periodic_cfg.lengths,
        n_repeat=periodic2_repeat,
    )

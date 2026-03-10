"""Potential computation helpers for mesh and sampled points."""

from __future__ import annotations

from typing import Iterable

import numpy as np

from .constants import K_COULOMB
from .mesh import _triangle_areas, _triangle_centers
from .selection import _require_triangles, _resolve_result
from .types import FortranRunResult, PotentialSlice2D


def compute_potential_mesh(
    result: FortranRunResult | object,
    *,
    softening: float = 0.0,
    self_term: str = "area_equivalent",
) -> np.ndarray:
    """Compute electric potential at each triangle centroid.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.
    softening : float, default 0.0
        Softening length in meters for off-diagonal interactions.
    self_term : {"area_equivalent", "exclude", "softened_point"}, default "area_equivalent"
        Self-interaction model at each centroid.

    Returns
    -------
    numpy.ndarray
        Potential values in volts with shape ``(mesh_nelem,)``.

    Raises
    ------
    ValueError
        If softening is negative, triangles are unavailable, or ``self_term`` is invalid.
    """

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


def compute_potential_points(
    result: FortranRunResult | object,
    points: np.ndarray,
    *,
    softening: float = 0.0,
    chunk_size: int = 2048,
) -> np.ndarray:
    """Compute electric potential at arbitrary 3D sampling points.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.
    points : numpy.ndarray
        Sampling points with shape ``(n_points, 3)`` in meters.
    softening : float, default 0.0
        Softening length in meters.
    chunk_size : int, default 2048
        Number of points processed per chunk.

    Returns
    -------
    numpy.ndarray
        Potential values in volts with shape ``(n_points,)``.

    Raises
    ------
    ValueError
        If argument shapes/ranges are invalid.
    """

    resolved = _resolve_result(result)
    if softening < 0.0:
        raise ValueError("softening must be >= 0.")
    if chunk_size <= 0:
        raise ValueError("chunk_size must be > 0.")

    sample_points = np.asarray(points, dtype=float)
    if sample_points.ndim != 2 or sample_points.shape[1] != 3:
        raise ValueError("points must have shape (n_points, 3).")
    if sample_points.shape[0] == 0:
        return np.empty(0, dtype=float)

    triangles = _require_triangles(resolved)
    centers = _triangle_centers(triangles)
    eps2 = softening * softening
    min_dist2 = np.finfo(float).tiny

    potential = np.empty(sample_points.shape[0], dtype=float)
    for start in range(0, sample_points.shape[0], chunk_size):
        stop = min(start + chunk_size, sample_points.shape[0])
        delta = sample_points[start:stop, None, :] - centers[None, :, :]
        dist2 = np.sum(delta * delta, axis=2) + eps2
        inv_r = 1.0 / np.sqrt(np.maximum(dist2, min_dist2))
        potential[start:stop] = inv_r @ resolved.charges

    return K_COULOMB * potential


def compute_potential_slices(
    result: FortranRunResult | object,
    *,
    box_min: Iterable[float],
    box_max: Iterable[float],
    grid_n: int = 200,
    xy_z: float | None = None,
    yz_x: float | None = None,
    xz_y: float | None = None,
    softening: float = 0.0,
    chunk_size: int = 2048,
) -> dict[str, PotentialSlice2D]:
    """Sample potential on XY/YZ/XZ slices inside a simulation box.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.
    box_min : iterable of float
        Lower simulation-box corner ``[x, y, z]`` in meters.
    box_max : iterable of float
        Upper simulation-box corner ``[x, y, z]`` in meters.
    grid_n : int, default 200
        Grid size per slice axis (``grid_n x grid_n``).
    xy_z : float or None, default None
        Z coordinate for the XY slice. ``None`` uses box center.
    yz_x : float or None, default None
        X coordinate for the YZ slice. ``None`` uses box center.
    xz_y : float or None, default None
        Y coordinate for the XZ slice. ``None`` uses box center.
    softening : float, default 0.0
        Softening length in meters.
    chunk_size : int, default 2048
        Number of sampling points processed per chunk.

    Returns
    -------
    dict[str, PotentialSlice2D]
        Slice data keyed by ``"xy"``, ``"yz"``, and ``"xz"``.

    Raises
    ------
    ValueError
        If box bounds, grid settings, or slice coordinates are invalid.
    """

    if grid_n < 2:
        raise ValueError("grid_n must be >= 2.")

    resolved = _resolve_result(result)
    min_corner, max_corner = _coerce_box_bounds(box_min, box_max)
    x = np.linspace(min_corner[0], max_corner[0], grid_n, dtype=float)
    y = np.linspace(min_corner[1], max_corner[1], grid_n, dtype=float)
    z = np.linspace(min_corner[2], max_corner[2], grid_n, dtype=float)

    z_xy = _coerce_slice_coordinate(
        xy_z,
        lower=min_corner[2],
        upper=max_corner[2],
        name="xy_z",
    )
    x_yz = _coerce_slice_coordinate(
        yz_x,
        lower=min_corner[0],
        upper=max_corner[0],
        name="yz_x",
    )
    y_xz = _coerce_slice_coordinate(
        xz_y,
        lower=min_corner[1],
        upper=max_corner[1],
        name="xz_y",
    )

    xx, yy = np.meshgrid(x, y, indexing="xy")
    points_xy = np.column_stack((xx.reshape(-1), yy.reshape(-1), np.full(xx.size, z_xy)))
    potential_xy = compute_potential_points(
        resolved,
        points_xy,
        softening=softening,
        chunk_size=chunk_size,
    ).reshape(grid_n, grid_n)

    yy2, zz = np.meshgrid(y, z, indexing="xy")
    points_yz = np.column_stack((np.full(yy2.size, x_yz), yy2.reshape(-1), zz.reshape(-1)))
    potential_yz = compute_potential_points(
        resolved,
        points_yz,
        softening=softening,
        chunk_size=chunk_size,
    ).reshape(grid_n, grid_n)

    xx2, zz2 = np.meshgrid(x, z, indexing="xy")
    points_xz = np.column_stack((xx2.reshape(-1), np.full(xx2.size, y_xz), zz2.reshape(-1)))
    potential_xz = compute_potential_points(
        resolved,
        points_xz,
        softening=softening,
        chunk_size=chunk_size,
    ).reshape(grid_n, grid_n)

    return {
        "xy": PotentialSlice2D(
            plane="xy",
            axis_u="x",
            axis_v="y",
            fixed_axis="z",
            fixed_value_m=float(z_xy),
            u_values_m=x,
            v_values_m=y,
            potential_V=potential_xy,
        ),
        "yz": PotentialSlice2D(
            plane="yz",
            axis_u="y",
            axis_v="z",
            fixed_axis="x",
            fixed_value_m=float(x_yz),
            u_values_m=y,
            v_values_m=z,
            potential_V=potential_yz,
        ),
        "xz": PotentialSlice2D(
            plane="xz",
            axis_u="x",
            axis_v="z",
            fixed_axis="y",
            fixed_value_m=float(y_xz),
            u_values_m=x,
            v_values_m=z,
            potential_V=potential_xz,
        ),
    }


def _coerce_box_bounds(
    box_min: Iterable[float], box_max: Iterable[float]
) -> tuple[np.ndarray, np.ndarray]:
    min_corner = np.asarray(list(box_min), dtype=float)
    max_corner = np.asarray(list(box_max), dtype=float)
    if min_corner.shape != (3,) or max_corner.shape != (3,):
        raise ValueError("box_min and box_max must contain exactly 3 values.")
    if np.any(max_corner <= min_corner):
        raise ValueError("box_max must be greater than box_min on all axes.")
    return min_corner, max_corner


def _coerce_slice_coordinate(
    value: float | None, *, lower: float, upper: float, name: str
) -> float:
    if value is None:
        return 0.5 * (lower + upper)
    coordinate = float(value)
    if coordinate < lower or coordinate > upper:
        raise ValueError(
            f"{name} must be within [{lower:.6e}, {upper:.6e}] but got {coordinate:.6e}."
        )
    return coordinate


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
        triangles,
        self_term=self_term,
        softening=softening,
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

"""Target-side quadrature for uniformly charged triangle panels."""

from __future__ import annotations

import operator

import numpy as np


def panel_target_quadrature(
    triangles_m: np.ndarray,
    element_charges_C: np.ndarray,
    order: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return Gauss-Duffy points and charge weights for triangle P0 targets."""

    quadrature_order = _quadrature_order(order)
    triangles = np.asarray(triangles_m, dtype=np.float64)
    charges = np.asarray(element_charges_C, dtype=np.float64)
    if triangles.ndim != 3 or triangles.shape[1:] != (3, 3):
        raise ValueError("triangles_m must have shape (n_triangles, 3, 3).")
    if triangles.shape[0] == 0:
        raise ValueError("triangles_m must contain at least one triangle.")
    if charges.shape != (triangles.shape[0],):
        raise ValueError(
            f"element_charges_C must have shape ({triangles.shape[0]},)."
        )
    if not np.all(np.isfinite(triangles)) or not np.all(np.isfinite(charges)):
        raise ValueError("triangle coordinates and element charges must be finite.")
    edge1 = triangles[:, 1] - triangles[:, 0]
    edge2 = triangles[:, 2] - triangles[:, 0]
    double_area = np.linalg.norm(np.cross(edge1, edge2), axis=1)
    if np.any(double_area <= 0.0) or not np.all(np.isfinite(double_area)):
        raise ValueError("triangles_m must contain non-degenerate triangles.")

    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    unit_nodes = 0.5 * (nodes + 1.0)
    unit_weights = 0.5 * weights
    u, v = np.meshgrid(unit_nodes, unit_nodes, indexing="ij")
    wu, wv = np.meshgrid(unit_weights, unit_weights, indexing="ij")
    u = u.reshape(-1)
    v = v.reshape(-1)
    reference_weight = (2.0 * wu * wv * (1.0 - u.reshape(wu.shape))).reshape(-1)
    reference_weight /= np.sum(reference_weight)

    bary1 = u
    bary2 = (1.0 - u) * v
    points = (
        triangles[:, None, 0, :]
        + bary1[None, :, None] * edge1[:, None, :]
        + bary2[None, :, None] * edge2[:, None, :]
    )
    charge_weights = charges[:, None] * reference_weight[None, :]
    charge_weights[:, -1] = charges - np.sum(charge_weights[:, :-1], axis=1)
    element_index = np.repeat(
        np.arange(triangles.shape[0], dtype=np.int64),
        quadrature_order * quadrature_order,
    )

    points_out = np.ascontiguousarray(points.reshape(-1, 3))
    charges_out = np.ascontiguousarray(charge_weights.reshape(-1))
    points_out.setflags(write=False)
    charges_out.setflags(write=False)
    element_index.setflags(write=False)
    return points_out, charges_out, element_index


def _quadrature_order(order: int) -> int:
    try:
        value = operator.index(order)
    except TypeError as exc:
        raise ValueError("panel target quadrature order must be integer 3 or 7.") from exc
    if value not in {3, 7}:
        raise ValueError("panel target quadrature order must be 3 or 7.")
    return value

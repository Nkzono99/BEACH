"""Coulomb force/torque computation utilities."""

from __future__ import annotations

from typing import Iterable, Literal, Mapping

import numpy as np

from .constants import K_COULOMB
from .mesh import _triangle_centers
from .selection import _coerce_group_selection, _resolve_result
from .types import CoulombInteraction, FortranRunResult, MeshSelection


def calc_coulomb(
    result: FortranRunResult | object,
    target: int | MeshSelection | Iterable[int | MeshSelection],
    source: int | MeshSelection | Iterable[int | MeshSelection],
    *,
    step: int | None = -1,
    softening: float = 0.0,
    torque_origin: Literal[
        "target_center",
        "source_center",
        "origin",
        "group_a_center",
        "group_b_center",
    ] = "target_center",
    periodic2: Mapping[str, object] | None = None,
) -> CoulombInteraction:
    """Compute Coulomb force/torque where target receives interaction from source.

    Parameters
    ----------
    result : FortranRunResult or Beach-like object
        Run result or object exposing ``result`` as ``FortranRunResult``.
    target : int, MeshSelection, or iterable of those
        Target mesh group (group A).
    source : int, MeshSelection, or iterable of those
        Source mesh group (group B).
    step : int or None, default -1
        History step used to read charges. ``-1`` selects latest history,
        ``None`` uses final charges from ``charges.csv``.
    softening : float, default 0.0
        Softening length in meters.
    torque_origin : {"target_center", "source_center", "origin", "group_a_center", "group_b_center"}, default "target_center"
        Reference point for torque computation.
    periodic2 : mapping or None, default None
        Two-axis periodic boundary setting. Uses a mapping with keys
        ``axes`` (length-2, 0-based axis indices), ``lengths`` (length-2,
        positive box lengths), and optional ``origins``, ``box_min``,
        ``image_layers`` (int, default 1).
        If ``None``, ``result.directory`` 近傍の ``beach.toml`` を探索し、
        ``sim.field_bc_mode="periodic2"`` なら自動適用する。

    Returns
    -------
    CoulombInteraction
        Aggregated force/torque summary for both groups.

    Raises
    ------
    ValueError
        If selection is empty, softening is negative, or arguments are invalid.
    """

    from .potential import _auto_periodic2_from_result, _coerce_periodic2

    resolved = _resolve_result(result)
    if softening < 0.0:
        raise ValueError("softening must be >= 0.")

    sel_target = _coerce_group_selection(resolved, target, step=step)
    sel_source = _coerce_group_selection(resolved, source, step=step)
    if sel_target.elem_indices.size == 0:
        raise ValueError("target does not contain any mesh elements.")
    if sel_source.elem_indices.size == 0:
        raise ValueError("source does not contain any mesh elements.")

    if torque_origin == "group_a_center":
        torque_origin = "target_center"
    elif torque_origin == "group_b_center":
        torque_origin = "source_center"

    if torque_origin == "target_center":
        origin = _triangle_centers(sel_target.triangles).mean(axis=0)
    elif torque_origin == "source_center":
        origin = _triangle_centers(sel_source.triangles).mean(axis=0)
    elif torque_origin == "origin":
        origin = np.zeros(3, dtype=float)
    else:
        raise ValueError(
            "torque_origin must be one of {'target_center', 'source_center', 'origin'}."
        )

    periodic_cfg = _coerce_periodic2(periodic2)
    if periodic_cfg is None:
        periodic_cfg = _auto_periodic2_from_result(resolved)

    centers_target = _triangle_centers(sel_target.triangles)
    centers_source = _triangle_centers(sel_source.triangles)
    force_target, torque_target = _pairwise_force_torque(
        centers_target,
        sel_target.charges,
        centers_source,
        sel_source.charges,
        origin=origin,
        softening=softening,
        periodic2=periodic_cfg,
    )
    force_source = -force_target
    torque_source = -torque_target

    return CoulombInteraction(
        group_a_mesh_ids=sel_target.mesh_ids,
        group_b_mesh_ids=sel_source.mesh_ids,
        step=sel_target.step,
        softening=softening,
        torque_origin_m=origin,
        force_on_a_N=force_target,
        force_on_b_N=force_source,
        torque_on_a_Nm=torque_target,
        torque_on_b_Nm=torque_source,
        mean_force_on_a_per_element_N=force_target / float(sel_target.elem_indices.size),
        mean_torque_on_a_per_element_Nm=torque_target
        / float(sel_target.elem_indices.size),
    )


def _pairwise_force_torque(
    centers_a: np.ndarray,
    charges_a: np.ndarray,
    centers_b: np.ndarray,
    charges_b: np.ndarray,
    *,
    origin: np.ndarray,
    softening: float,
    periodic2: tuple | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute pairwise Coulomb force and torque.

    When *periodic2* is provided, image shells of the source charges are
    included so that the nearest-cell summation accounts for periodic
    boundary conditions along two axes.
    """

    if periodic2 is not None:
        return _pairwise_force_torque_periodic2(
            centers_a,
            charges_a,
            centers_b,
            charges_b,
            origin=origin,
            softening=softening,
            periodic2=periodic2,
        )

    force = np.zeros(3, dtype=float)
    torque = np.zeros(3, dtype=float)
    eps2 = softening * softening
    min_dist2 = np.finfo(float).tiny

    for i in range(centers_a.shape[0]):
        delta = centers_a[i] - centers_b
        dist2 = np.sum(delta * delta, axis=1) + eps2
        inv_r3 = 1.0 / (
            np.maximum(dist2, min_dist2) * np.sqrt(np.maximum(dist2, min_dist2))
        )
        coeff = K_COULOMB * charges_a[i] * charges_b * inv_r3
        f_i = np.sum(coeff[:, None] * delta, axis=0)
        force += f_i
        torque += np.cross(centers_a[i] - origin, f_i)

    return force, torque


def _pairwise_force_torque_periodic2(
    centers_a: np.ndarray,
    charges_a: np.ndarray,
    centers_b: np.ndarray,
    charges_b: np.ndarray,
    *,
    origin: np.ndarray,
    softening: float,
    periodic2: tuple,
) -> tuple[np.ndarray, np.ndarray]:
    """Pairwise Coulomb force/torque with periodic2 image-shell summation."""

    axes, lengths, origins, nimg, _far_correction, _alpha, _ewald_layers = periodic2
    axis1, axis2 = axes
    l1, l2 = lengths

    force = np.zeros(3, dtype=float)
    torque = np.zeros(3, dtype=float)
    eps2 = softening * softening
    min_dist2 = np.finfo(float).tiny

    for ix in range(-nimg, nimg + 1):
        for iy in range(-nimg, nimg + 1):
            shifted_b = centers_b.copy()
            shifted_b[:, axis1] += float(ix) * l1
            shifted_b[:, axis2] += float(iy) * l2

            for i in range(centers_a.shape[0]):
                delta = centers_a[i] - shifted_b
                dist2 = np.sum(delta * delta, axis=1) + eps2
                inv_r3 = 1.0 / (
                    np.maximum(dist2, min_dist2)
                    * np.sqrt(np.maximum(dist2, min_dist2))
                )
                coeff = K_COULOMB * charges_a[i] * charges_b * inv_r3
                f_i = np.sum(coeff[:, None] * delta, axis=0)
                force += f_i
                torque += np.cross(centers_a[i] - origin, f_i)

    return force, torque

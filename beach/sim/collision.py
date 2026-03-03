from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .mesh import BEMMesh


@dataclass
class HitInfo:
    mesh_id: int
    elem_idx: int
    t: float
    pos: np.ndarray
    normal: np.ndarray


def _segment_triangle_intersect(p0, p1, v0, v1, v2, eps=1e-12):
    """Möller–Trumbore with segment constraint t in [0,1]."""
    d = p1 - p0
    e1 = v1 - v0
    e2 = v2 - v0
    h = np.cross(d, e2)
    a = np.dot(e1, h)
    if abs(a) < eps:
        return None
    f = 1.0 / a
    s = p0 - v0
    u = f * np.dot(s, h)
    if u < 0.0 or u > 1.0:
        return None
    q = np.cross(s, e1)
    v = f * np.dot(d, q)
    if v < 0.0 or (u + v) > 1.0:
        return None
    t = f * np.dot(e2, q)
    if t < 0.0 or t > 1.0:
        return None
    hit = p0 + t * d
    return float(t), hit


def find_first_hit_in_mesh(
    mesh: BEMMesh, p0: np.ndarray, p1: np.ndarray
) -> tuple[Optional[int], Optional[float], Optional[np.ndarray]]:
    """Returns (elem_idx, t, hitpos) for earliest intersection, or (None, None, None)."""
    if mesh.nelem == 0:
        return None, None, None

    seg_min = np.minimum(p0, p1)
    seg_max = np.maximum(p0, p1)
    cand = np.where(
        np.all(mesh.bb_max >= seg_min, axis=1) & np.all(mesh.bb_min <= seg_max, axis=1)
    )[0]

    best_t = None
    best_hit = None
    best_idx = None
    for i in cand:
        out = _segment_triangle_intersect(p0, p1, mesh.v0[i], mesh.v1[i], mesh.v2[i])
        if out is None:
            continue
        t, hit = out
        if (best_t is None) or (t < best_t):
            best_t, best_hit, best_idx = t, hit, int(i)

    return best_idx, best_t, best_hit


def find_first_hit(
    bem_list: list[BEMMesh], p0: np.ndarray, p1: np.ndarray
) -> Optional[HitInfo]:
    """Checks all meshes and returns earliest hit among them."""
    best: Optional[HitInfo] = None
    for mid, mesh in enumerate(bem_list):
        idx, t, hit = find_first_hit_in_mesh(mesh, p0, p1)
        if idx is None:
            continue
        if (best is None) or (t < best.t):
            best = HitInfo(
                mesh_id=mid,
                elem_idx=idx,
                t=t,
                pos=hit,
                normal=mesh.normals[idx].copy() if mesh.nelem else np.zeros(3),
            )
    return best

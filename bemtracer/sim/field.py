from __future__ import annotations

from typing import Optional, Protocol

import numpy as np

from .constants import K_COULOMB
from .mesh import BEMMesh


def _triangle_subcentroids(v0, v1, v2, n_sub: int) -> np.ndarray:
    if n_sub <= 1:
        return ((v0 + v1 + v2) / 3.0)[None, :]

    e1 = v1 - v0
    e2 = v2 - v0
    pts = []
    n = n_sub
    for i in range(n):
        for j in range(n - i):
            p00 = v0 + (i / n) * e1 + (j / n) * e2
            p10 = v0 + ((i + 1) / n) * e1 + (j / n) * e2
            p01 = v0 + (i / n) * e1 + ((j + 1) / n) * e2
            pts.append((p00 + p10 + p01) / 3.0)

            if j < (n - 1 - i):
                p11 = v0 + ((i + 1) / n) * e1 + ((j + 1) / n) * e2
                pts.append((p10 + p11 + p01) / 3.0)

    return np.asarray(pts, dtype=float)


def _E_point_charges(
    r: np.ndarray, centers: np.ndarray, charges: np.ndarray, softening: float
) -> np.ndarray:
    r = np.asarray(r, dtype=float)
    if centers.shape[0] == 0:
        return np.zeros_like(r, dtype=float)

    if r.ndim == 1:
        dr = r[None, :] - centers
        r2 = np.einsum("ij,ij->i", dr, dr) + softening * softening
        inv_r3 = 1.0 / (np.sqrt(r2) * r2)
        E = K_COULOMB * np.einsum("i,ij->j", charges * inv_r3, dr)
        return E

    dr = r[..., None, :] - centers[None, ...]
    r2 = np.einsum("...ij,...ij->...i", dr, dr) + softening * softening
    inv_r3 = 1.0 / (np.sqrt(r2) * r2)
    E = K_COULOMB * np.einsum("...i,...ij->...j", charges[None, ...] * inv_r3, dr)
    return E


class MagneticFieldModel(Protocol):
    def B(self, r: np.ndarray, t: float) -> np.ndarray: ...


class ZeroB:
    def B(self, r: np.ndarray, t: float) -> np.ndarray:
        return np.zeros(3, dtype=float)


class UniformB:
    def __init__(self, B0: np.ndarray):
        self.B0 = np.asarray(B0, dtype=float).reshape(3)

    def B(self, r: np.ndarray, t: float) -> np.ndarray:
        return self.B0


class BEMField:
    """Electric field from bem_list charges."""

    def __init__(
        self,
        bem_list: list[BEMMesh],
        use_hybrid: bool = True,
        r_switch_factor: float = 3.0,
        n_sub: int = 2,
        softening_factor: float = 0.1,
    ):
        self.bem_list = bem_list
        self.use_hybrid = use_hybrid
        self.r_switch_factor = float(r_switch_factor)
        self.n_sub = int(n_sub)
        self.softening_factor = float(softening_factor)
        self._centers_all: Optional[np.ndarray] = None
        self._charges_all: Optional[np.ndarray] = None
        self._softening: Optional[float] = None

    def rebuild_cache(self) -> None:
        centers_all = []
        charges_all = []
        h_refs = []
        for mesh in self.bem_list:
            if mesh.nelem == 0:
                continue
            centers_all.append(mesh.centers)
            charges_all.append(mesh.charges())
            h_refs.append(np.median(mesh.h_elem) if mesh.nelem > 0 else 1.0)

        self._centers_all = (
            np.vstack(centers_all) if centers_all else np.zeros((0, 3), dtype=float)
        )
        self._charges_all = (
            np.concatenate(charges_all) if charges_all else np.zeros((0,), dtype=float)
        )
        h_ref = float(np.median(np.array(h_refs))) if len(h_refs) > 0 else 1.0
        self._softening = self.softening_factor * h_ref

    def E(self, r: np.ndarray) -> np.ndarray:
        if (
            self._centers_all is None
            or self._charges_all is None
            or self._softening is None
        ):
            self.rebuild_cache()

        centers_all = self._centers_all
        charges_all = self._charges_all
        soft = float(self._softening)
        E_far = _E_point_charges(r, centers_all, charges_all, softening=soft)
        if not self.use_hybrid:
            return E_far

        E_corr = np.zeros(3, dtype=float)
        for mesh in self.bem_list:
            if mesh.nelem == 0:
                continue

            h_ref = float(np.median(mesh.h_elem))
            r_switch = self.r_switch_factor * h_ref
            idxs = mesh.near_elements_by_center(r, r_switch=r_switch)
            if idxs.size == 0:
                continue

            centers = mesh.centers[idxs]
            q = mesh.charges()[idxs]
            E_approx_near = _E_point_charges(r, centers, q, softening=soft)

            E_exact_near = np.zeros(3, dtype=float)
            n_sub = max(1, self.n_sub)
            for ei in idxs:
                qj = float(mesh.elements[int(ei)].q)
                if qj == 0.0:
                    continue
                pts = _triangle_subcentroids(
                    mesh.v0[ei], mesh.v1[ei], mesh.v2[ei], n_sub=n_sub
                )
                qk = (qj / pts.shape[0]) * np.ones((pts.shape[0],), dtype=float)
                E_exact_near += _E_point_charges(r, pts, qk, softening=soft)

            E_corr += E_exact_near - E_approx_near

        return E_far + E_corr


def calc_electric_field(x, y, z, bem_list: list[BEMMesh], **kwargs):
    field = BEMField(bem_list, **kwargs)
    r = np.stack([np.asarray(x), np.asarray(y), np.asarray(z)], axis=-1)

    if r.ndim == 1:
        E = field.E(r)
        return E[0], E[1], E[2]

    Eout = np.zeros_like(r, dtype=float)
    it = np.nditer(r[..., 0], flags=["multi_index"])
    while not it.finished:
        idx = it.multi_index
        Eout[idx] = field.E(r[idx])
        it.iternext()
    return Eout[..., 0], Eout[..., 1], Eout[..., 2]

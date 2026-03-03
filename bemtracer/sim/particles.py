from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np

from .collision import HitInfo
from .mesh import BEMMesh


@dataclass
class Particle:
    x: np.ndarray
    v: np.ndarray
    q: float
    m: float
    w: float = 1.0
    alive: bool = True


class InteractionModel(Protocol):
    def on_collision(self, p: Particle, hit: HitInfo) -> tuple[bool, float]: ...


class AbsorptionInteraction:
    """Absorb particle at collision. Deposit particle charge to surface."""

    def on_collision(self, p: Particle, hit: HitInfo) -> tuple[bool, float]:
        dq = p.q * p.w
        return False, dq


class SurfaceChargeModel(Protocol):
    def begin_batch(self, nelem_total: int) -> None: ...

    def deposit(self, hit: HitInfo, dq: float) -> None: ...

    def commit_batch(self, bem_list: list[BEMMesh]) -> tuple[float, float]: ...


class InsulatorSurfaceCharge:
    """Insulator: simply accumulate dq on hit elements; no redistribution."""

    def __init__(self):
        self._dq_per_mesh: list[np.ndarray] = []

    def begin_batch(self, nelem_total: int) -> None:
        self._dq_per_mesh = []

    def _ensure_buffers(self, bem_list: list[BEMMesh]) -> None:
        if self._dq_per_mesh:
            return
        self._dq_per_mesh = [np.zeros(mesh.nelem, dtype=float) for mesh in bem_list]

    def deposit(self, hit: HitInfo, dq: float) -> None:
        self._dq_per_mesh[hit.mesh_id][hit.elem_idx] += dq

    def commit_batch(self, bem_list: list[BEMMesh]) -> tuple[float, float]:
        self._ensure_buffers(bem_list)

        dq_all = (
            np.concatenate(self._dq_per_mesh)
            if self._dq_per_mesh
            else np.zeros((0,), dtype=float)
        )
        norm_dq = float(np.linalg.norm(dq_all))

        for mesh, dq in zip(bem_list, self._dq_per_mesh):
            mesh.add_charges(dq)

        q_all = (
            np.concatenate([mesh.charges() for mesh in bem_list])
            if bem_list
            else np.zeros((0,), dtype=float)
        )
        norm_q = float(np.linalg.norm(q_all))

        self._dq_per_mesh = [np.zeros(mesh.nelem, dtype=float) for mesh in bem_list]
        return norm_dq, norm_q

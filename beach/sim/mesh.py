from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class BEMElement:
    """Triangular boundary element with total charge q [C]."""

    v0: np.ndarray
    v1: np.ndarray
    v2: np.ndarray
    q: float = 0.0

    def centroid(self) -> np.ndarray:
        return (self.v0 + self.v1 + self.v2) / 3.0

    def area(self) -> float:
        return 0.5 * np.linalg.norm(np.cross(self.v1 - self.v0, self.v2 - self.v0))

    def normal(self) -> np.ndarray:
        n = np.cross(self.v1 - self.v0, self.v2 - self.v0)
        nn = np.linalg.norm(n)
        return n / nn if nn > 0 else n


class BEMMesh:
    """Triangular surface mesh with per-element charge."""

    def __init__(self, elements: list[BEMElement]):
        self.elements = elements
        self.nelem = len(elements)

        v0 = (
            np.stack([e.v0 for e in elements], axis=0)
            if self.nelem
            else np.zeros((0, 3), dtype=float)
        )
        v1 = (
            np.stack([e.v1 for e in elements], axis=0)
            if self.nelem
            else np.zeros((0, 3), dtype=float)
        )
        v2 = (
            np.stack([e.v2 for e in elements], axis=0)
            if self.nelem
            else np.zeros((0, 3), dtype=float)
        )
        self.v0 = v0
        self.v1 = v1
        self.v2 = v2

        self.centers = (
            np.stack([e.centroid() for e in elements], axis=0)
            if self.nelem
            else np.zeros((0, 3), dtype=float)
        )
        self.areas = (
            np.array([e.area() for e in elements], dtype=float)
            if self.nelem
            else np.zeros((0,), dtype=float)
        )
        self.normals = (
            np.stack([e.normal() for e in elements], axis=0)
            if self.nelem
            else np.zeros((0, 3), dtype=float)
        )

        self.bb_min = (
            np.minimum(np.minimum(v0, v1), v2)
            if self.nelem
            else np.zeros((0, 3), dtype=float)
        )
        self.bb_max = (
            np.maximum(np.maximum(v0, v1), v2)
            if self.nelem
            else np.zeros((0, 3), dtype=float)
        )

        self.h_elem = (
            np.sqrt(np.maximum(self.areas, 0.0))
            if self.nelem
            else np.zeros((0,), dtype=float)
        )

    def charges(self) -> np.ndarray:
        return np.array([e.q for e in self.elements], dtype=float)

    def add_charges(self, dq: np.ndarray) -> None:
        """Add dq [C] to each element. dq shape must be (nelem,)."""
        if dq.shape != (self.nelem,):
            raise ValueError(f"dq must have shape ({self.nelem},), got {dq.shape}")
        for i, d in enumerate(dq):
            if d != 0.0:
                self.elements[i].q += float(d)

    def near_elements_by_center(self, x: np.ndarray, r_switch: float) -> np.ndarray:
        """Cheap near query: centroid distance < r_switch (+ element size margin)."""
        if self.nelem == 0:
            return np.zeros((0,), dtype=int)
        dr = self.centers - x[None, :]
        d2 = np.einsum("ij,ij->i", dr, dr)
        margin = 0.5 * self.h_elem + 1e-15
        return np.where(d2 <= (r_switch + margin) ** 2)[0]

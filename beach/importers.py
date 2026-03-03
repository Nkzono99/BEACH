from __future__ import annotations

from pathlib import Path
from typing import Protocol

import numpy as np

from .sim import BEMMesh
from .geometry import mesh_from_vertices_faces


class MeshImporter(Protocol):
    """Interface for geometry import backends."""

    def load(self, path: str | Path) -> BEMMesh: ...


class OBJImporter:
    """Minimal Wavefront OBJ importer (v/f only)."""

    def load(self, path: str | Path) -> BEMMesh:
        p = Path(path)
        vertices: list[list[float]] = []
        faces: list[list[int]] = []

        with p.open("r", encoding="utf-8") as f:
            for raw in f:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                if line.startswith("v "):
                    _, x, y, z, *_ = line.split()
                    vertices.append([float(x), float(y), float(z)])
                elif line.startswith("f "):
                    tokens = line.split()[1:]
                    idx = []
                    for tok in tokens:
                        head = tok.split("/")[0]
                        if not head:
                            raise ValueError(f"unsupported face token in OBJ: {tok}")
                        vi = int(head)
                        if vi < 0:
                            vi = len(vertices) + vi + 1
                        idx.append(vi - 1)
                    if len(idx) < 3:
                        continue
                    for i in range(1, len(idx) - 1):
                        faces.append([idx[0], idx[i], idx[i + 1]])

        if not vertices or not faces:
            raise ValueError(f"no vertices/faces found in OBJ: {p}")

        return mesh_from_vertices_faces(np.asarray(vertices), np.asarray(faces))


class FreeCADImporter:
    """Planned importer hook for FreeCAD mesh export formats."""

    def load(self, path: str | Path) -> BEMMesh:
        raise NotImplementedError(
            "FreeCAD direct import is not implemented yet. "
            "Use FreeCAD to export OBJ, then load with OBJImporter for now."
        )


def load_mesh(path: str | Path, importer: MeshImporter | None = None) -> BEMMesh:
    """Load a mesh from file using either explicit importer or extension-based default."""
    if importer is not None:
        return importer.load(path)

    suffix = Path(path).suffix.lower()
    if suffix == ".obj":
        return OBJImporter().load(path)

    raise ValueError(
        f"unsupported mesh format: {suffix}. "
        "Pass a custom importer, or convert to OBJ."
    )

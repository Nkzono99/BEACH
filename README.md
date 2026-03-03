# bemtracer

Boundary Element Method (BEM) surface-charging + test-particle prototype in Python.

This repository implements:
- Triangular boundary mesh with per-element charge (insulator accumulation mode)
- Test particle pusher (Boris, E+B ready)
- Segment–triangle collision (first-hit)
- Batched charge deposition (`npcls_per_step`) with electric-field cache rebuild
- Hybrid electric field: far-field centroid point-charge + near-field triangle subdivision correction (difference form)
- Geometry templates: plane / box / cylinder / sphere
- Import extension points for external CAD mesh workflow (OBJ available, FreeCAD hook prepared)

## Quickstart

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
source .venv/bin/activate
pip install -e .
python examples/minimal_run.py
# テンプレート形状の生成
python examples/template_meshes.py
# charge 分布を 3D 可視化する例（matplotlib が必要）
python examples/visualize_bem_list_3d.py
```

## Template geometry API

```python
from bemtracer import make_plane, make_box, make_cylinder, make_sphere

plane = make_plane(size_x=1.0, size_y=0.8, nx=20, ny=16)
box = make_box(size=(1.0, 0.8, 0.6), nx=6, ny=6, nz=4)
cyl = make_cylinder(radius=0.5, height=2.0, n_theta=64, n_z=24, cap=True)
sph = make_sphere(radius=0.5, n_lon=64, n_lat=32)
```

## Mesh import API (extensible)

```python
from bemtracer import load_mesh, OBJImporter

mesh = load_mesh("model.obj")
mesh2 = OBJImporter().load("another.obj")
```

`FreeCADImporter` is intentionally left as a design hook for future direct integration.
For now, export from FreeCAD to OBJ and load via `OBJImporter`.

## Files
- `SPEC.md` : specification (v0.1)
- `bemtracer/core.py` : simulator core
- `bemtracer/geometry.py` : reusable template geometry builders
- `bemtracer/importers.py` : file import interface and OBJ importer
- `examples/minimal_run.py` : minimal runnable example

## Notes
- v0.1 focuses on correctness and extensibility. Performance optimizations (BVH/Octree/FMM) are out of scope for now.
- Offset evaluation (δ-shift) is intentionally *not* used in v0.1.

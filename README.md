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
# ランダム注入（位置一様 + シフテッドマクスウェル速度）
python examples/random_injection_shifted_maxwell.py
```

## Template geometry API

```python
from bemtracer import make_plane, make_box, make_cylinder, make_sphere

plane = make_plane(size_x=1.0, size_y=0.8, nx=20, ny=16)
box = make_box(size=(1.0, 0.8, 0.6), nx=6, ny=6, nz=4)
cyl = make_cylinder(radius=0.5, height=2.0, n_theta=64, n_z=24, cap=True)
sph = make_sphere(radius=0.5, n_lon=64, n_lat=32)
```


## Custom particle injection API

```python
import numpy as np
from bemtracer import (
    RandomBeamInjector,
    ShiftedMaxwellVelocitySampler,
    UniformPositionSampler,
)

rng = np.random.default_rng(123)
pos_sampler = UniformPositionSampler(
    low=np.array([0.1, 0.1, 0.45]),
    high=np.array([0.4, 0.4, 0.55]),
    rng=rng,
)
vel_sampler = ShiftedMaxwellVelocitySampler(
    drift_velocity=np.array([0.0, 0.0, -2.0e5]),
    temperature_K=2.0e4,
    rng=rng,
)
injector = RandomBeamInjector(
    q=-1.602176634e-19,
    m=9.10938356e-31,
    w=1.0,
    position_sampler=pos_sampler,
    velocity_sampler=vel_sampler,
)
particles = injector.sample(500)
```

- `UniformPositionSampler`: 直方体領域で一様サンプル
- `ShiftedMaxwellVelocitySampler`: 各速度成分を `N(0, sigma^2)` + ドリフトで生成
  - `temperature_K` 指定時は `sigma = sqrt(k_B T / m)`
  - あるいは `thermal_speed` を直接指定可能
- 既存の固定ビーム相当には `FixedBeamInjector` を利用可能

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

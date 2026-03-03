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
- `bemtracer/sim/` : simulator core (責務別に分割した実装)
- `bemtracer/geometry.py` : reusable template geometry builders
- `bemtracer/importers.py` : file import interface and OBJ importer
- `examples/minimal_run.py` : minimal runnable example

## Notes
- v0.1 focuses on correctness and extensibility. Performance optimizations (BVH/Octree/FMM) are out of scope for now.
- Offset evaluation (δ-shift) is intentionally *not* used in v0.1.

## Fortran + OpenMP implementation (fpm)

Python版とは別に、コア計算をFortranで実装した最小構成を追加しています。

- `fpm.toml`: fpmビルド設定
- `src/*.f90`: Fortranモジュール（mesh / field / boris pusher / collision / simulator）
- `app/main.f90`: 吸収+絶縁体蓄積モードのデモ実行

### Build / Run

```bash
fpm run --profile release --flag "-fopenmp"
```


Fortranデモは以下を実行します（引数なし時は `examples/simple_plate.obj` がある場合はOBJ読込、なければ平面テンプレート生成）。

```bash
fpm run --profile release --flag "-fopenmp"
```

設定をパラメータファイルから与える場合は、実行ファイルへ `TOML` を渡せます。

```bash
fpm run --profile release --flag "-fopenmp" -- examples/fortran_config.toml
```

- `examples/fortran_config.toml` は `[[mesh.templates]]` を複数並べることで、template形状を合成した境界を作成できます（同じ `kind` の重複指定も可能）。

Fortran側でもPython版に揃えて、`bem_templates` モジュールで plane/box/cylinder/sphere の境界テンプレート生成を提供し、`bem_importers` モジュールでOBJ（三角形/多角形面の三角形分割対応）を読み込めます。
さらに `bem_injection` モジュールで、`sample_uniform_positions`（位置一様サンプリング）と `sample_shifted_maxwell_velocities`（シフテッド・マクスウェル速度サンプリング）および `init_random_beam_particles` を提供し、Python版と同様のランダム注入をFortranでも利用できます。

OpenMPスレッド数は環境変数で制御できます。

```bash
OMP_NUM_THREADS=8 fpm run --profile release --flag "-fopenmp"
```

## Fortran結果のファイル出力とPython連携

Fortran実行時に、最終結果を `summary.txt` と `charges.csv` として出力できます（デフォルト有効）。

```toml
[output]
write_files = true
dir = "outputs/latest"
```

- `summary.txt`: 粒子処理統計（absorbed / escaped など）
- `charges.csv`: 要素ごとの最終電荷

Python側では `bemtracer.fortran_results` で読み込み・管理・可視化できます。

```python
from bemtracer import load_fortran_result, list_fortran_runs, plot_charges

result = load_fortran_result("outputs/latest")
print(result.absorbed, result.escaped, result.charges.sum())

fig, ax = plot_charges(result)
```

CLIサンプル: 

```bash
python examples/inspect_fortran_output.py outputs/latest --save outputs/latest/charges.png
```

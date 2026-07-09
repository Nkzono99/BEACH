title: Python Post-processing API Reference

Lang: [English](PythonPostprocessAPI.en.md) | [日本語](PythonPostprocessAPI.md)

# Python Post-processing API Reference

The BEACH Python package (`beach`) is the post-processing layer for loading, analyzing, and visualizing Fortran simulation results.
It reads files emitted by the Fortran runtime, such as `summary.txt`, `charges.csv`, and `mesh_triangles.csv`, and performs potential reconstruction, Coulomb force calculations, electric-field calculations, field-line tracing, and 3D visualization on the Python side.

## 1. Package Structure

| Module | Role |
|---|---|
| `beach.fortran_results.io` | Reads output directories (`load_fortran_result`, `list_fortran_runs`) |
| `beach.fortran_results.facade` | High-level facade class `Beach` |
| `beach.fortran_results.potential` | Potential reconstruction (`compute_potential_mesh`, `compute_potential_points`, `compute_potential_slices`) |
| `beach.fortran_results.coulomb` | Coulomb force/torque calculation (`calc_coulomb`) |
| `beach.fortran_results.kernel` | Shared-library calls to the Fortran FMM field kernel (`FieldKernel`, `calc_object_forces_kernel`) |
| `beach.fortran_results.scene` | Temporary object translation/rotation and field-kernel evaluation of edited scenes |
| `beach.fortran_results.field_lines` | Electric-field calculation, field-line tracing, and 3D plotting (`compute_electric_field_points`, `trace_field_lines`, `plot_field_lines_3d`) |
| `beach.fortran_results.mobility` | Coulomb mobility analysis (`analyze_coulomb_mobility`) |
| `beach.fortran_results.plotting` | Plotting utilities (`plot_charge_mesh`, `plot_charges`, `plot_potential_mesh`, etc.) |
| `beach.fortran_results.animation` | History animation (`animate_history_mesh`) |
| `beach.fortran_results.history` | Batch-step access to `charge_history.csv` (`FortranChargeHistory`) |
| `beach.fortran_results.types` | Public data types (`FortranRunResult`, `CoulombInteraction`, etc.) |
| `beach.fortran_results.constants` | Physical constants (`K_COULOMB`) |

All public symbols can be imported from the `beach` top level and from `beach.fortran_results`.

```python
from beach import Beach, calc_coulomb, compute_electric_field_points, trace_field_lines
```

## 2. `Beach` Facade Class

This high-level interface binds one output directory and provides the main analysis and visualization methods.

```python
b = Beach("outputs/latest")
```

You can explicitly specify the configuration file as in `Beach("outputs/latest", config_path="path/to/beach.toml")`. When `config_path=None`, BEACH searches automatically in order through `output_dir/beach.toml`, the parent directory, and the grandparent directory. Config-aware object/kernel analyses resolve object kind/order, `sim.softening`, periodic2, and tree parameters from this `beach.toml`.

### 2.1 Constructor

| Parameter | Type | Default | Description |
|---|---|---|---|
| `output_dir` | `str \| Path` | `"outputs/latest"` | Fortran output directory |

### 2.2 Properties

| Name | Return type | Description |
|---|---|---|
| `result` | `FortranRunResult` | Loaded result, loaded lazily |
| `mesh_ids` | `tuple[int, ...]` | Available mesh IDs |

### 2.3 Methods

| Method | Delegates to | Summary |
|---|---|---|
| `reload()` | `load_fortran_result` | Reload from disk |
| `get_mesh(*mesh_ids, step)` | Internal | Get a `MeshSelection` by mesh ID |
| `get_mesh_charge(*mesh_ids, step)` | Internal | Get an element charge array by mesh ID |
| `calc_coulomb(target, source, ...)` | `calc_coulomb` | Coulomb force/torque calculation |
| `calc_object_forces_kernel(...)` | `calc_object_forces_kernel` | Per-object resultant-force calculation using the Fortran field kernel |
| `scene(step, ...)` | `BeachScene.from_result` | What-if scene with temporary object translation/rotation |
| `analyze_coulomb_mobility(...)` | `analyze_coulomb_mobility` | Per-object mobility analysis |
| `compute_potential(...)` | `compute_potential_mesh` | Potential reconstruction at centroids |
| `compute_potential_points(points, ...)` | `compute_potential_points` | Potential at arbitrary points |
| `compute_potential_slices(...)` | `compute_potential_slices` | Potential on XY/YZ/XZ sections |
| `compute_electric_field(points, ...)` | `compute_electric_field_points` | Electric-field vectors at arbitrary points |
| `trace_field_lines(seed_points, ...)` | `trace_field_lines` | RK4 tracing of electric field lines |
| `plot_mesh(...)` | `plot_charge_mesh` | 3D mesh plot of charge density |
| `plot_potential(...)` | `plot_potential_mesh` | 3D mesh plot of potential |
| `plot_potential_slices(...)` | `plot_potential_slices` | Potential-section plots |
| `plot_field_lines(seed_points, ...)` | `plot_field_lines_3d` | 3D plot of electric field lines |
| `plot_bar()` | `plot_charges` | Bar chart of element charge |
| `plot_mesh_source_boxplot(...)` | `plot_mesh_source_boxplot` | Box plot by mesh source |
| `plot_coulomb_force_matrix(...)` | `plot_coulomb_force_matrix` | Coulomb force matrix plot |
| `animate_mesh(...)` | `animate_history_mesh` | Charge/potential history animation |

## 3. Loading Results

### 3.1 `load_fortran_result(directory)`

Loads a Fortran output directory and returns a `FortranRunResult`.

```python
from beach import load_fortran_result

result = load_fortran_result("outputs/latest")
print(f"要素数: {result.mesh_nelem}, バッチ数: {result.batches}")
print(f"吸収: {result.absorbed}, 脱出: {result.escaped}")
```

Required files: `summary.txt`, `charges.csv`
Optional files: `mesh_triangles.csv`, `mesh_sources.csv`, `charge_history.csv`, `mesh_potential.csv`

### 3.2 `FortranRunResult` Type

| Field | Type | Description |
|---|---|---|
| `directory` | `Path` | Output directory path |
| `mesh_nelem` | `int` | Number of mesh elements |
| `processed_particles` | `int` | Number of processed particles |
| `absorbed` | `int` | Number of absorbed particles |
| `escaped` | `int` | Number of escaped particles |
| `batches` | `int` | Number of processed batches |
| `escaped_boundary` | `int` | Number of particles that escaped through boundaries |
| `survived_max_step` | `int` | Number of particles that reached max_step |
| `last_rel_change` | `float` | Final relative charge change |
| `charges` | `ndarray (mesh_nelem,)` | Element charge array [C] |
| `triangles` | `ndarray (mesh_nelem, 3, 3) \| None` | Triangle vertex coordinates [m] |
| `mesh_ids` | `ndarray (mesh_nelem,) \| None` | Element mesh IDs |
| `mesh_sources` | `dict[int, MeshSource] \| None` | Mesh kind, surface model, and epsilon_r metadata |
| `mesh_potential_v` | `ndarray (mesh_nelem,) \| None` | Centroid potentials output by Fortran [V] |
| `history` | `FortranChargeHistory \| None` | Charge history accessor |

## 4. Potential Reconstruction

### 4.1 `compute_potential_mesh(result, *, softening, self_term, periodic2, reference_point)`

Reconstructs the potential at triangle centroids. If Fortran has already output `mesh_potential.csv` and the conditions match, that output is preferred.

| Parameter | Type | Default | Unit | Description |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (required) | - | Result object |
| `softening` | `float \| None` | `None` | m | When `None`, automatically references `sim.softening` |
| `self_term` | `str` | `"auto"` | - | Self interaction: `"auto"` / `"area_equivalent"` / `"exclude"` / `"softened_point"` |
| `periodic2` | `Mapping \| None` | `None` | - | Two-axis periodic settings, described below. `None` enables automatic detection |
| `reference_point` | `Iterable[float] \| str \| None` | `None` | m | Reference potential point. `"species1_injection_center"` means the species 1 injection-face center |

Return value: `ndarray (mesh_nelem,)` [V]

### 4.2 `compute_potential_points(result, points, *, softening, chunk_size, periodic2, reference_point)`

Computes the potential at arbitrary 3D points.

| Parameter | Type | Default | Unit | Description |
|---|---|---|---|---|
| `points` | `ndarray (n_points, 3)` | (required) | m | Sampling point coordinates |
| `softening` | `float \| None` | `None` | m | Automatic when `None` |
| `chunk_size` | `int` | `2048` | - | Chunk size |
| `periodic2` | `Mapping \| None` | `None` | - | Automatic detection when `None` |
| `reference_point` | `Iterable[float] \| str \| None` | `None` | m | Reference potential point |

Return value: `ndarray (n_points,)` [V]

### 4.3 `compute_potential_slices(result, *, box_min, box_max, grid_n, xy_z, yz_x, xz_y, ...)`

Computes potential sections on the XY/YZ/XZ planes.

Return value: `dict[str, PotentialSlice2D]` (keys: `"xy"`, `"yz"`, `"xz"`)

## 5. Coulomb Force/Torque Calculation

### 5.1 `calc_coulomb(result, target, source, *, step, softening, torque_origin, periodic2)`

Computes the Coulomb force and torque that the target mesh group receives from the source mesh group.

| Parameter | Type | Default | Unit | Description |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (required) | - | Result object |
| `target` | `int \| MeshSelection \| Iterable` | (required) | - | Target mesh group (group A) |
| `source` | `int \| MeshSelection \| Iterable` | (required) | - | Source mesh group (group B) |
| `step` | `int \| None` | `-1` | - | History step. `-1` means latest; `None` means final charge |
| `softening` | `float` | `0.0` | m | Softening length |
| `torque_origin` | `str` | `"target_center"` | - | Torque reference point: `"target_center"` / `"source_center"` / `"origin"` |
| `periodic2` | `Mapping \| None` | `None` | - | Two-axis periodic boundary settings. `None` enables automatic detection |

Return value: `CoulombInteraction`

### Periodic Coulomb Sum Using the `periodic2` Parameter

When `periodic2` is specified, source-charge image shells `ix in [-nimg, nimg], iy in [-nimg, nimg]` are generated in the two periodic directions, and the Coulomb force/torque is computed as a nearest-cell sum.

When `periodic2=None` (default), BEACH searches for `beach.toml` near the output directory and automatically applies periodic settings if `sim.field_bc_mode="periodic2"` is set. This uses the same automatic-detection logic as functions such as `compute_potential_mesh`.

### 5.2 `CoulombInteraction` Type

| Field | Type | Unit | Description |
|---|---|---|---|
| `group_a_mesh_ids` | `tuple[int, ...]` | - | Target mesh IDs |
| `group_b_mesh_ids` | `tuple[int, ...]` | - | Source mesh IDs |
| `step` | `int \| None` | - | History step used |
| `softening` | `float` | m | Softening length used |
| `torque_origin_m` | `ndarray (3,)` | m | Torque reference point |
| `force_on_a_N` | `ndarray (3,)` | N | Net force acting on group A |
| `force_on_b_N` | `ndarray (3,)` | N | Net force acting on group B |
| `torque_on_a_Nm` | `ndarray (3,)` | N m | Net torque acting on group A |
| `torque_on_b_Nm` | `ndarray (3,)` | N m | Net torque acting on group B |
| `mean_force_on_a_per_element_N` | `ndarray (3,)` | N | Mean force per target element |
| `mean_torque_on_a_per_element_Nm` | `ndarray (3,)` | N m | Mean torque per target element |

### Example

```python
from beach import Beach

b = Beach("outputs/latest")

# mesh_id=0 が mesh_id=1 から受ける Coulomb 力
interaction = b.calc_coulomb(target=0, source=1)
print(f"Force on target: {interaction.force_on_a_N} [N]")
print(f"Torque on target: {interaction.torque_on_a_Nm} [N m]")

# periodic2 を明示指定する場合
interaction_p = b.calc_coulomb(
    target=0, source=1,
    periodic2={"axes": [0, 1], "lengths": [0.01, 0.01], "image_layers": 2},
)
```

## 6. Electric Field Calculation

### 6.1 `compute_electric_field_points(result, points, *, softening, chunk_size, periodic2)`

Computes electric-field vectors at arbitrary 3D points directly from surface charges using Coulomb's law.

Formula: `E(r) = K * sum_j q_j * (r - r_j) / |r - r_j|^3`

Here, `r_j` is the centroid of triangle element `j`, and `q_j` is the element charge.

| Parameter | Type | Default | Unit | Description |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (required) | - | Result object |
| `points` | `ndarray (n_points, 3)` | (required) | m | Sampling point coordinates |
| `softening` | `float \| None` | `None` | m | Softening length. When `None`, automatically references `sim.softening` |
| `chunk_size` | `int` | `2048` | - | Chunk size |
| `periodic2` | `Mapping \| None` | `None` | - | Two-axis periodic settings. `None` enables automatic detection |

Return value: `ndarray (n_points, 3)` [V/m]

In periodic2 mode, sampling points are wrapped into the periodic cell, and source-charge image shells are superposed with `ix in [-nimg, nimg], iy in [-nimg, nimg]`.

### Example

```python
import numpy as np
from beach import Beach

b = Beach("outputs/latest")

# グリッド点での電場計算
x = np.linspace(0.0, 0.01, 50)
y = np.linspace(0.0, 0.01, 50)
xx, yy = np.meshgrid(x, y)
zz = np.full_like(xx, 0.005)
points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])

efield = b.compute_electric_field(points)
print(f"Shape: {efield.shape}")   # (2500, 3)
print(f"E-field [V/m]: {efield[0]}")
```

## 7. Electric Field-Line Tracing

### 7.1 `trace_field_lines(result, seed_points, *, ds, max_steps, softening, periodic2, direction, box_min, box_max)`

Traces electric field lines from seed points in the electric-field direction, or in the opposite direction, using RK4 integration.

| Parameter | Type | Default | Unit | Description |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (required) | - | Result object |
| `seed_points` | `ndarray (n_seeds, 3)` | (required) | m | Starting point coordinates for field lines |
| `ds` | `float \| None` | `None` | m | Integration step size. When `None`, automatically set from average mesh edge length x 0.5 |
| `max_steps` | `int` | `500` | - | Maximum number of integration steps in each direction |
| `softening` | `float \| None` | `None` | m | Softening length. Automatic when `None` |
| `periodic2` | `Mapping \| None` | `None` | - | Two-axis periodic settings. `None` enables automatic detection |
| `direction` | `str` | `"both"` | - | Trace direction: `"forward"` (electric-field direction) / `"backward"` (opposite direction) / `"both"` (both directions) |
| `box_min` | `Iterable[float] \| None` | `None` | m | Lower boundary-box limit. Stop when a field line exits this box |
| `box_max` | `Iterable[float] \| None` | `None` | m | Upper boundary-box limit |

Return value: `list[ndarray]` -- each element is a field-line coordinate array with shape `(n_points_i, 3)` [m]

#### Tracing Algorithm

- Advances by `ds` along the unit vector of the electric field using fourth-order Runge-Kutta (RK4)
- Stops if the electric-field norm becomes smaller than `1e-30` at any RK4 stage
- When `direction="both"`, joins the forward and backward traces, removing the duplicated seed point
- Stops once a trace exceeds `box_min` / `box_max`

#### Limitations

- This is a direct-sum calculation on the Python side, so large meshes with tens of thousands of elements or more can take a long time
- Fortran treecode/fmm is not used. Python-side electric-field calculation is only a direct sum over point charges at element centroids
- In periodic2, only explicit image shells are reconstructed; oracle residuals (Ewald far correction) are not reproduced on the Python side

### Example

```python
import numpy as np
from beach import Beach

b = Beach("outputs/latest")

# シード点を手動指定
seeds = np.array([
    [0.005, 0.005, 0.008],
    [0.005, 0.005, 0.002],
])

lines = b.trace_field_lines(seeds, max_steps=1000)
print(f"力線数: {len(lines)}")
for i, line in enumerate(lines):
    print(f"  力線 {i}: {line.shape[0]} 点")
```

### 7.2 `plot_field_lines_3d(result, seed_points, *, ...)`

Draws electric field lines in 3D, optionally overlaying the mesh surface colored by charge density.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `result` | `FortranRunResult \| object` | (required) | Result object |
| `seed_points` | `ndarray (n_seeds, 3)` | (required) | Seed points [m] |
| `ds` | `float \| None` | `None` | Integration step size [m] |
| `max_steps` | `int` | `500` | Maximum number of steps |
| `softening` | `float \| None` | `None` | Softening length [m] |
| `periodic2` | `Mapping \| None` | `None` | Periodic settings |
| `direction` | `str` | `"both"` | Trace direction |
| `box_min` | `Iterable[float] \| None` | `None` | Lower boundary-box limit [m] |
| `box_max` | `Iterable[float] \| None` | `None` | Upper boundary-box limit [m] |
| `show_mesh` | `bool` | `True` | Show mesh overlay |
| `mesh_alpha` | `float` | `0.25` | Mesh-face transparency |
| `mesh_cmap` | `str` | `"coolwarm"` | Mesh surface charge density colormap |
| `line_color` | `str \| None` | `None` | Fixed field-line color. `None` enables coloring by `line_cmap` |
| `line_cmap` | `str` | `"plasma"` | Field-line colormap, used when `line_color=None` |
| `line_width` | `float` | `1.2` | Field-line width |
| `view_elev` | `float` | `24.0` | Elevation angle [deg] |
| `view_azim` | `float` | `-58.0` | Azimuth angle [deg] |
| `title` | `str` | `"Electric field lines"` | Plot title |
| `figsize` | `tuple[float, float]` | `(9, 7)` | Figure size [inch] |

Return value: `(figure, axes)` -- matplotlib Figure / Axes3D

#### Drawn Elements

- Each field line is drawn as a line. When `line_color=None`, each line is colored using `line_cmap`
- Direction arrows (quiver) are drawn near the midpoint of each field line
- Seed points are drawn as red scatter points
- When `show_mesh=True`, the triangle mesh is colored by surface charge density `q / A` and overlaid semi-transparently

### Example

```python
import numpy as np
from beach import Beach

b = Beach("outputs/latest")

seeds = np.array([
    [0.005, 0.005, 0.008],
    [0.003, 0.007, 0.006],
    [0.007, 0.003, 0.006],
])

fig, ax = b.plot_field_lines(
    seeds,
    max_steps=800,
    direction="both",
    show_mesh=True,
    mesh_alpha=0.3,
    line_cmap="viridis",
    view_elev=30,
    view_azim=-45,
)
fig.savefig("field_lines.png", dpi=150)
```

## 8. `periodic2` Parameter Specification

The common `periodic2` parameter is available for potential reconstruction, Coulomb force calculation, electric-field calculation, and field-line tracing.

### 8.1 Automatic Detection (`periodic2=None`)

With the default `None`, BEACH searches for `beach.toml` near the output directory and automatically applies periodic boundary settings when `sim.field_bc_mode="periodic2"` is set.
If no configuration file is found, or if `field_bc_mode` is not `periodic2`, calculations are performed in free-space mode.

### 8.2 Explicit Specification

Specify `periodic2` as a `Mapping` with the following keys.

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `axes` | `list[int]` (length 2) | required | - | 0-based indices of the periodic axes, for example `[0, 1]` for the x and y axes |
| `lengths` | `list[float]` (length 2) | required | - | Box length on each periodic axis [m]. Must be positive |
| `origins` | `list[float]` (length 2) | - | `[0.0, 0.0]` | Box origin on each periodic axis [m] |
| `box_min` | `list[float]` (length 3) | - | - | Alternative to `origins`. Extracts periodic-axis origins from the lower limit of the 3D box |
| `image_layers` | `int` | - | `1` | Number of image-shell layers. Evaluates `[-N, N]` on each periodic axis |
| `far_correction` | `str` | - | `"none"` | `"auto"` / `"none"` / `"m2l_root_oracle"`. `auto` is treated as `"none"` for compatibility. Kept on the Python side for configuration compatibility, but the oracle residual itself is not reproduced |
| `ewald_alpha` | `float` | - | `0.0` | Ewald decomposition parameter, reserved |
| `ewald_layers` | `int` | - | `4` | Ewald truncation depth, reserved |

If both `origins` and `box_min` are specified, `origins` takes precedence.

### Example

```python
p2 = {
    "axes": [0, 1],
    "lengths": [0.01, 0.01],
    "image_layers": 2,
}
potential = b.compute_potential(periodic2=p2)
efield = b.compute_electric_field(points, periodic2=p2)
interaction = b.calc_coulomb(target=0, source=1, periodic2=p2)
lines = b.trace_field_lines(seeds, periodic2=p2)
```

### 8.3 Limitations of the Python-side periodic2 Implementation

- On the Python side, periodic sums are reconstructed only by direct sums over explicit image shells
- The Ewald far correction using explicit `m2l_root_oracle` in the Fortran FMM is not reproduced on the Python side. `far_correction` is kept for configuration compatibility, but it does not affect Python direct-sum calculations
- Larger `image_layers` improves accuracy but increases cost by a factor of `(2*N+1)^2`

## 9. Coulomb Mobility Analysis

### 9.1 `analyze_coulomb_mobility(result, *, step, softening, config_path, gravity, support_normal, ...)`

Analyzes the tendency of each object to slide, roll, or lift under Coulomb forces.

Return value: `CoulombMobilityAnalysis` (stores a tuple of `CoulombMobilityRecord` in `.records`)

## 10. Fortran Field Kernel Integration

### 10.1 `FieldKernel`

Loads `build/libbeach_field_kernel.so`, generated by `make build-kernel`, through `ctypes` and evaluates electric fields and potentials using the same Fortran FMM core as the simulation. It reads softening, periodic2, and tree settings from `config_path` or an automatically discovered `beach.toml`.

```python
from beach import Beach, FieldKernel

run = Beach("outputs/latest")
with FieldKernel.from_result(run) as kernel:
    e = kernel.eval_e([[0.0, 0.0, 0.01]])
```

If the shared library is located elsewhere, specify `library_path=` or the environment variable `BEACH_FIELD_KERNEL_LIB`.

### 10.2 `calc_object_forces_kernel(result, ...)`

For each object, zeros its own source charges and calculates `sum(q_i E_not_self(r_i))` and torque. This path uses the Fortran kernel and passes the periodic2 settings explicitly defined in `beach.toml`, including `m2l_root_oracle`, so it is a diagnostic closer to the simulation-side field definition than the Python direct-sum `calc_coulomb`.

```python
from beach import Beach

run = Beach("outputs/latest")
records = run.calc_object_forces_kernel()
for record in records:
    print(record.mesh_id, record.total_charge_C, record.force_N)
```

### 10.3 `BeachScene`

`Beach.scene()` is a what-if view that temporarily edits an already charged mesh on the Python side.
`move` / `rotate` returns a new scene; each element keeps its charge attached to the same element while only centroids and vertices are rigidly transformed. Calling `calc_object_forces_kernel` afterward passes the edited geometry to the Fortran field kernel as source/target geometry.

```python
from beach import Beach

run = Beach("outputs/latest", config_path="examples/beach.toml")
scene = run.scene()
moved = scene.move(2, by=[1.0e-3, 0.0, 0.0]).rotate(
    2,
    axis=[0.0, 0.0, 1.0],
    angle_deg=15.0,
)
records = moved.calc_object_forces_kernel(target_mesh_ids=[2])
print(records[0].force_N, records[0].torque_Nm)
```

By default, rigid transformations on the Python side are processed with NumPy. To use Numba, install the optional dependency with `pip install ".[accel]"` and specify `run.scene(transform_backend="numba")`. The main calculation that defines the meaning of FMM, periodic2, and far correction is still performed by the Fortran kernel.

## 11. Visualization Functions

### 11.1 Charge/Potential Mesh Plots

```python
fig, ax = b.plot_mesh(cmap="coolwarm")                        # 電荷密度
fig, ax = b.plot_potential(reference_point="species1_injection_center")  # 電位
```

### 11.2 Potential Sections

```python
fig, axes = b.plot_potential_slices(
    box_min=[0, 0, 0], box_max=[0.01, 0.01, 0.01],
    xy_z=0.005,
)
```

### 11.3 History Animation

```python
gif_path = b.animate_mesh("charge_animation.gif", quantity="charge", fps=10)
```

### 11.4 Coulomb Force Matrix

```python
fig, ax = b.plot_coulomb_force_matrix(component="z")
```

## 12. CLI Commands

### 12.1 Unified CLI (`beachx`)

Since v1.0.0, the unified `beachx` CLI is recommended.

| Command | Description |
|---|---|
| `beachx inspect <output_dir>` | Show an output-directory summary |
| `beachx animate <output_dir>` | Generate an animated GIF of charge/potential history |
| `beachx workload <config.toml>` | Estimate workload |
| `beachx slices <output_dir>` | Plot potential sections |
| `beachx profile <output_dir>` | Plot a performance profile |
| `beachx coulomb <output_dir>` | Plot a Coulomb force matrix |
| `beachx mobility <output_dir>` | Run Coulomb mobility analysis |
| `beachx kernel-forces <output_dir>` | Output per-object resultant-force CSV using the Fortran field kernel |
| `beachx lint <config.toml>` | Check TOML / JSON Schema / BEACH constraints |
| `beachx config render <config.toml>` | Expand high-level notation into final `beach.toml` keys |
| `beachx config validate <config.toml>` | Validate a configuration file |
| `beachx model close-pack` | Generate a close-packed model |

### 12.2 Legacy CLI (Deprecated)

The legacy entry points below are retained for backward compatibility, but they may be removed in a future version.

| Command | Description |
|---|---|
| `beach-inspect <output_dir>` | Show an output-directory summary |
| `beach-animate-history <output_dir>` | Generate an animated GIF of charge/potential history |
| `beach-estimate-workload <config.toml>` | Estimate workload |
| `beach-plot-potential-slices <output_dir>` | Plot potential sections |
| `beach-plot-performance-profile <output_dir>` | Plot a performance profile |
| `beach-plot-coulomb-force-matrix <output_dir>` | Plot a Coulomb force matrix |

## 12. Physical Constants

| Symbol | Value | Unit | Description |
|---|---|---|---|
| `K_COULOMB` | `8.9875517923e9` | N m^2 / C^2 | Coulomb constant |

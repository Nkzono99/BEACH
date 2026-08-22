title: Python Post-processing API Reference

Lang: [English](PythonPostprocessAPI.en.md) | [日本語](PythonPostprocessAPI.md)

# Python Post-processing API Reference

The BEACH Python package (`beach`) is the post-processing layer for loading, analyzing, and visualizing Fortran simulation results.
It reads files emitted by the Fortran runtime, such as `summary.txt`, `charges.csv`, and `mesh_triangles.csv`, and performs potential reconstruction, Coulomb force calculations, electric-field calculations, field-line tracing, and 3D visualization on the Python side.

If you only need the first plots, start with
[Post-processing Tutorial](PostprocessTutorial.en.html).

## 1. Package Structure

| Module | Role |
|---|---|
| `beach.fortran_results.io` | Reads output directories (`load_fortran_result`, `list_fortran_runs`) |
| `beach.fortran_results.facade` | High-level facade class `Beach` |
| `beach.fortran_results.potential` | Potential reconstruction (`compute_potential_mesh`, `compute_potential_points`, `compute_potential_slices`) |
| `beach.fortran_results.coulomb` | Coulomb force/torque calculation (`calc_coulomb`) |
| `beach.fortran_results.kernel` | Shared-library calls to the native Fortran field kernel (`FieldKernel`, `calc_object_forces_kernel`) |
| `beach.fortran_results.object_interaction` | Object force, torque, and vertical paths against frozen source charges (`ObjectInteractionSnapshot`, `ObjectProbe`) |
| `beach.fortran_results.detachment` | Immutable results for path work, adhesion, gravity, speed, and the from-rest barrier |
| `beach.fortran_results.periodic_force_oracle` | Finite periodic-image shells and `E_bottom=0` closure convergence diagnostics |
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

You can explicitly specify the configuration file as in `Beach("outputs/latest", config_path="path/to/beach.toml")`. When
`config_path=None`, BEACH searches automatically in order through `output_dir/beach.toml`, the parent directory, and the
grandparent directory. Config-aware object analyses resolve object kind/order and related metadata from this `beach.toml`.
For new outputs, total-field reconstruction reads the resolved boundary, box,
uniform-field, periodic2, actual solver, and tree/FMM parameters from the
`field_reconstruction_*` receipt in `summary.txt`. When this receipt is present, a nearby or explicitly
selected `beach.toml` does not override its field settings. Only legacy outputs
without the receipt require the original `beach.toml` rather than guessing
boundary or external-field policies. A configuration is also not required
when reading an unchanged simulator-written `mesh_potential.csv`.

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
| `object_interaction_snapshot(...)` | `ObjectInteractionSnapshot.from_result` | Primary-only self exclusion and detachment paths against frozen sources |
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
Optional files: `mesh_triangles.csv`, `mesh_sources.csv`, `charge_history.csv`, `mesh_potential.csv`,
`matching_plane_history.csv`

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
| `matching_plane_state` | `MatchingPlaneState \| None` | Last accepted matching-plane state; `None` when the summary marks it invalid |
| `matching_plane_history` | `tuple[MatchingPlaneHistoryEntry, ...] \| None` | Matching-plane history; `None` when the state is invalid or the CSV is absent |

`MatchingPlaneState` exposes named fields for $D_H$, $\Phi_H$, inward responses, outward feedback, PE return / escape,
iteration count, and convergence residual. Each `MatchingPlaneHistoryEntry` contains `batch`, `simulated_time_s`, and the
corresponding `state: MatchingPlaneState`. A header-only `matching_plane_history.csv` produces the empty tuple `()`.
Thus `None` means that no history file or no valid matching-plane state exists, whereas an empty tuple means that the
state is valid but the history contains no rows.

### 3.3 `FortranChargeHistory`

Fortran `charge_history.csv` files contain dense snapshots. Every recorded batch
must contain every `elem_idx=1..mesh_nelem` exactly once. Element rows may appear
in any order, but batch groups must be strictly increasing, `processed_particles`
and `rel_change` must be constant within a batch, and both `charge_C` and
`rel_change` must be finite.

`load_fortran_result(...)` and history-accessor construction remain lazy. The
first access to a history property such as `batch_indices`, to `get_step(...)`,
or to `as_array()` builds the byte-offset index and validates every batch in the
entire CSV. Missing or duplicate elements, out-of-range indices, non-finite
values, inconsistent metadata, or decreasing batches raise `ValueError` with
the affected batch and defect. A missing element is corruption, not a physical
charge of `0 C`.

After validation, individual charge vectors are still loaded on demand, and the
full history matrix is not materialized until `as_array()` is called.
`FortranChargeHistory.from_arrays(...)` retains its existing behavior as a
trusted in-memory path whose caller supplies the dense matrix.

## 4. Potential Reconstruction

### 4.1 `compute_potential_mesh(result, *, periodic2, reference_point, config_path, library_path)`

Evaluates potential at triangle centroids with the native P0 triangle panel
kernel. If Fortran has already output `mesh_potential.csv` and the conditions
match, that output is preferred.

| Parameter | Type | Default | Unit | Description |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (required) | - | Result object |
| `periodic2` | `Mapping \| None` | `None` | - | Two-axis periodic settings, described below. `None` enables automatic detection |
| `reference_point` | `Iterable[float] \| str \| None` | `None` | m | Reference potential point. `"species1_injection_center"` means the species 1 injection-face center |
| `config_path` | `str \| Path \| None` | `None` | - | Explicit `beach.toml` path |
| `library_path` | `str \| Path \| None` | `None` | - | Explicit native field-kernel library path |

Return value: `ndarray (mesh_nelem,)` [V]

### 4.2 `compute_potential_points(result, points, *, chunk_size, periodic2, reference_point, config_path, library_path)`

Computes the potential at arbitrary 3D points.

| Parameter | Type | Default | Unit | Description |
|---|---|---|---|---|
| `points` | `ndarray (n_points, 3)` | (required) | m | Sampling point coordinates |
| `chunk_size` | `int` | `2048` | - | Chunk size |
| `periodic2` | `Mapping \| None` | `None` | - | Automatic detection when `None` |
| `reference_point` | `Iterable[float] \| str \| None` | `None` | m | Reference potential point |
| `config_path` | `str \| Path \| None` | `None` | - | Explicit `beach.toml` path |
| `library_path` | `str \| Path \| None` | `None` | - | Explicit native field-kernel library path |

Return value: `ndarray (n_points,)` [V]

### 4.3 `compute_potential_slices(result, *, box_min, box_max, grid_n, xy_z, yz_x, xz_y, ...)`

Computes potential sections on the XY/YZ/XZ planes.

Return value: `dict[str, PotentialSlice2D]` (keys: `"xy"`, `"yz"`, `"xz"`)

## 5. Coulomb Force/Torque Calculation

### 5.1 `calc_coulomb(result, target, source, *, step, torque_origin, periodic2, quadrature_order, config_path, library_path)`

Computes the Coulomb force and torque that the target mesh group receives from the source mesh group.

| Parameter | Type | Default | Unit | Description |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (required) | - | Result object |
| `target` | `int \| MeshSelection \| Iterable` | (required) | - | Target mesh group (group A) |
| `source` | `int \| MeshSelection \| Iterable` | (required) | - | Source mesh group (group B) |
| `step` | `int \| None` | `-1` | - | History step. `-1` means latest; `None` means final charge |
| `torque_origin` | `str` | `"target_center"` | - | Torque reference point: `"target_center"` / `"source_center"` / `"origin"` |
| `periodic2` | `Mapping \| None` | `None` | - | Two-axis periodic boundary settings. `None` enables automatic detection |
| `quadrature_order` | `int` | `7` | - | Gauss-Duffy integration order for target panels |
| `config_path` | `str \| Path \| None` | `None` | - | Explicit `beach.toml` path |
| `library_path` | `str \| Path \| None` | `None` | - | Explicit native field-kernel library path |

Return value: `CoulombInteraction`

### Periodic Coulomb Sum Using the `periodic2` Parameter

When `periodic2` is specified, the native field kernel is constructed with
that two-axis periodic configuration for the force and torque evaluation.

When `periodic2=None` (default), BEACH searches for `beach.toml` near the output directory. If
`field_boundary.mode="periodic2"`, it builds periodic settings from `domain.periodic_axes` and `domain.box_min` /
`domain.box_max`. Functions such as `compute_potential_mesh` use the same automatic detection.

### 5.2 `CoulombInteraction` Type

| Field | Type | Unit | Description |
|---|---|---|---|
| `group_a_mesh_ids` | `tuple[int, ...]` | - | Target mesh IDs |
| `group_b_mesh_ids` | `tuple[int, ...]` | - | Source mesh IDs |
| `step` | `int \| None` | - | History step used |
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

### 6.1 `compute_electric_field_points(result, points, *, chunk_size, periodic2, config_path, library_path)`

Computes electric-field vectors at arbitrary 3D points with the native P0
triangle panel kernel.

| Parameter | Type | Default | Unit | Description |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (required) | - | Result object |
| `points` | `ndarray (n_points, 3)` | (required) | m | Sampling point coordinates |
| `chunk_size` | `int` | `2048` | - | Chunk size |
| `periodic2` | `Mapping \| None` | `None` | - | Two-axis periodic settings. `None` enables automatic detection |
| `config_path` | `str \| Path \| None` | `None` | - | Explicit `beach.toml` path |
| `library_path` | `str \| Path \| None` | `None` | - | Explicit native field-kernel library path |

Return value: `ndarray (n_points, 3)` [V/m]

In periodic2 mode, the native kernel evaluates finite images or the cached
nonzero mode according to the configuration.

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

### 7.1 `trace_field_lines(result, seed_points, *, ds, max_steps, periodic2, direction, box_min, box_max, config_path, library_path)`

Traces electric field lines from seed points in the electric-field direction, or in the opposite direction, using RK4 integration.

| Parameter | Type | Default | Unit | Description |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (required) | - | Result object |
| `seed_points` | `ndarray (n_seeds, 3)` | (required) | m | Starting point coordinates for field lines |
| `ds` | `float \| None` | `None` | m | Integration step size. When `None`, automatically set from average mesh edge length x 0.5 |
| `max_steps` | `int` | `500` | - | Maximum number of integration steps in each direction |
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

- The native field-kernel shared library is required.
- Every RK4 stage evaluates the field kernel, so cost grows with the seed count and `max_steps`.
- Tracing does not stop on mesh collision; it stops at the field threshold, `max_steps`, or the configured box.

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

With the default `None`, BEACH first restores the resolved field boundary,
periodic axes, lengths, and origins from the output's
`field_reconstruction_*` receipt. For a legacy output without the receipt, it
searches for `beach.toml` near the output directory. When
`field_boundary.mode="periodic2"`, it resolves the periodic geometry from
`domain.periodic_axes` and `domain.box_min` / `domain.box_max`. If no
configuration file is found for a legacy output, high-level total-field
reconstruction stops instead of silently falling back to free space.
Free-space reconstruction is used when the field boundary resolves to `free`.
A split-periodic `panel_spectral_reference` receipt is not substituted with the
finite-image native kernel: high-level total-field reconstruction stops and
directs the caller to simulator-written potential output or the dedicated
periodic-reference validation workflow.
When the receipt's actual solver is `direct`, automatic reconstruction uses the
non-periodic exact-direct kernel. For `fmm`, it uses the recorded FMM expansion
order. A resolved `treecode` stops all high-level field, potential, and force
APIs instead of being substituted with FMM, even when `periodic2` is explicit.

### 8.2 Explicit Specification

Specify `periodic2` as a `Mapping` with the following keys.

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `axes` | `list[int]` (length 2) | required | - | 0-based indices of the periodic axes, for example `[0, 1]` for the x and y axes |
| `lengths` | `list[float]` (length 2) | required | - | Box length on each periodic axis [m]. Must be positive |
| `origins` | `list[float]` (length 2) | - | `[0.0, 0.0]` | Box origin on each periodic axis [m] |
| `box_min` | `list[float]` (length 3) | - | - | Alternative to `origins`. Extracts periodic-axis origins from the lower limit of the 3D box |
| `image_layers` | `int` | - | `1` | Number of image-shell layers. Evaluates `[-N, N]` on each periodic axis |
| `far_correction` | `str` | - | `"none"` | Explicit inputs accept `"auto"` / `"none"`. `auto` is treated as `"none"` for compatibility |
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
- The native `cached_kneq0` operator returns only the $k\ne0$ component. General potential, field, and force APIs reject it rather than presenting it as the total field. Use a saved `mesh_potential.csv` or `ObjectInteractionSnapshot`, which explicitly composes the physical zero mode
- Removed `m2l_root_oracle` metadata is accepted only when `periodic2=None` auto-discovers historical metadata near the output, and is normalized to `none` for the finite image shell. Explicit `periodic2` values and config paths are rejected. Python does not reproduce the Ewald correction
- Larger `image_layers` improves accuracy but increases cost by a factor of `(2*N+1)^2`

## 9. Coulomb Mobility Analysis

### 9.1 `analyze_coulomb_mobility(result, *, step, config_path, library_path, gravity, support_normal, ...)`

Analyzes the tendency of each object to slide, roll, or lift under Coulomb forces.

Return value: `CoulombMobilityAnalysis` (stores a tuple of `CoulombMobilityRecord` in `.records`)

## 10. Fortran Field Kernel Integration

### 10.1 `FieldKernel`

Loads `build/libbeach_field_kernel.so`, generated by `make build-kernel`,
through `ctypes` and evaluates P0 triangle panel fields and potentials using the
same Fortran field core as the simulation. It reads periodic2 and tree settings
and the actual solver from the `field_reconstruction_*` receipt for new outputs,
or from `config_path` or an automatically discovered `beach.toml` for legacy
outputs.

```python
from beach import Beach, FieldKernel

run = Beach("outputs/latest")
with FieldKernel.from_result(run) as kernel:
    e = kernel.eval_e([[0.0, 0.0, 0.01]])
    phi = kernel.eval_phi([[0.0, 0.0, 0.01]])
```

If the shared library is located elsewhere, specify `library_path=` or the environment variable `BEACH_FIELD_KERNEL_LIB`.

`eval_e()` and `eval_phi()` use the free, finite-periodic, or cached-periodic
configuration with which the plan was built.
Automatic construction maps resolved `direct` to exact direct and resolved
`fmm` to FMM; resolved `treecode` fails closed. The receipt's FMM expansion
order is the default, and only an explicit `order=` overrides it. With resolved
`direct`, `eval_e()`, `eval_phi()`, and `force_on_charges()` all use exact direct
and compose uniform `sim.e0` exactly once.
For a cached plan, `eval_e()` and `eval_phi()` return the `cached_kneq0`
$k\ne0$ component plus `sim.e0`, not the simulator total field
with its physical zero mode. `FieldKernel` is the low-level API for handling
that component explicitly. `eval_e_direct()` and
`eval_phi_direct()` are **non-periodic exact-direct** diagnostics over the same
source geometry and charges. The direct methods reject periodic plans and do
not add uniform `sim.e0`. They are intended for small FMM accuracy or primary
free-space subtraction oracles, not as an infinite-periodic replacement.

### 10.2 `calc_object_forces_kernel(result, ...)`

For each object, this API zeros its own source charges and calculates
`sum(q_i E_not_self(r_i))` and torque. Its existing self policy is
`exclude_target_lattice`: in a periodic calculation it removes both the target's
primary source and that object's periodic images.
This high-level API rejects an uncomposed `cached_kneq0` field.

To retain an object's own periodic images when evaluating detachment, use
`ObjectInteractionSnapshot.object_probe(...)` as described below. Its self
policy is fixed to `exclude_primary_keep_images`.
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

### 10.4 `ObjectInteractionSnapshot` and Frozen-Source Paths

This API freezes all saved source geometry and charge once, then moves only the
selected central-cell object as an independent target probe. The fixed
`exclude_primary_keep_images` policy subtracts only the target's
central-cell free-space primary contribution. It retains that target's periodic
images, all images of other objects, and the uniform field.

```python
import numpy as np
from beach import AdhesionProfile, Beach

run = Beach("outputs/latest", config_path="beach.toml")
with run.object_interaction_snapshot(
    periodic_model="infinite_physical",
    cache_dir=".beach_cache/periodic2",
) as snapshot:
    probe = snapshot.object_probe(6)
    wrench = probe.wrench()
    path = probe.vertical_path(np.linspace(0.0, 2.0e-4, 65))

release = path.evaluate_release(
    mass_kg=2.0e-12,
    gravity_m_s2=9.80665,
    adhesion=AdhesionProfile.finite_range_constant(
        force_N=1.0e-10,
        range_m=2.0e-6,
    ),
)
print(wrench.force_N, wrench.torque_Nm)
print(path.status, path.work_relative_mismatch)
print(release.barrier_free_from_rest, release.endpoint_speed_m_s)
```

`periodic_model` has the following meanings.

| Value | Field definition |
|---|---|
| `"configured"` | Uses the run's `beach.toml` unchanged. The result can therefore be free space, a finite shell with `far_correction="none"`, or cached periodic |
| `"infinite_physical"` | For an x/y-periodic run, combines cached `k != 0` with the physical `k = 0` mode selected by `[periodic2].lower_boundary_model`. Cache generation or reuse must succeed |

A complete `beach.toml` with `domain.box_min` and `domain.box_max` is required.

When an x/y-periodic mesh crosses a cell seam, the snapshot keeps all-source
geometry in the saved representation so it remains identical to the simulation
and cache identity. `object_probe()` unwraps only the selected mesh into one
periodically connected branch. It uses that same target geometry for quadrature,
central-primary removal, rigid transforms, the area centroid, and the bounding
radius. Inspect it through `probe.target_geometry_representation`,
`target_triangles_m`, `geometric_area_centroid_m`, `vertex_bounding_center_m`,
and `vertex_bounding_radius_m`.

`ObjectWrench.components` preserves these physical contributions:

| Key | Meaning |
|---|---|
| `other_objects_all_images` | Other objects and all of their periodic images |
| `target_periodic_images` | The target's own periodic images, with its central primary removed |
| `external_uniform` | Configured uniform external field |
| `total_external` | Sum of the three rows above; equals `ObjectWrench.force_N` / `torque_Nm` |

Kernel and zero-mode entries in `numerical_metadata` are a numerical
decomposition of the same total, not additional physical forces. Torque depends
on its reference point, so retain `ObjectWrench.torque_origin_m` and
`numerical_metadata["torque_origin_policy"]` with every force record.
`vertical_path()` records the reference point at every height in
`numerical_metadata["torque_origin_m"]`.

Target integration uses order-7 Gauss-Duffy area integration by default.
`target_integration="centroid_compatibility"` is an explicitly labelled
comparison with historical centroid integration and is never selected by
`auto`.

Mechanical force on a surface uses the principal-value (PV) zero-mode trace.
The simulator particle pusher uses the one-sided `zero_mode_trace_plus` value at
the surface, so the two traces serve different purposes. For cached results,
`numerical_metadata["cached_kneq0_trace_correction"]` is a diagnostic that has
**already been folded into `periodic_kneq0`** to form the PV decomposition. Do
not add it again to `force_N` or `periodic_kneq0`.

`vertical_path()` holds source geometry and charges at their initial positions
and translates only target quadrature in `+z`. Every target triangle vertex
must remain inside the non-periodic box/interface. Because the environment is
frozen, its potential energy is `U_env = sum_i(q_i phi_env(r_i))`, with no factor
of `1/2`. Always inspect the agreement between
`electrostatic_work_J = integral(F_z dh)` and
`potential_difference_work_J = U_env(0)-U_env(h)`, adaptive-refinement
diagnostics, and `status`.

`AdhesionProfile.finite_range_constant(F, d)` supplies a resisting force for
`0 <= h < d` and adhesion work that saturates at `F*d`. `evaluate_release()`
checks the full continuous piecewise path after gravity, adhesion, and optional
dissipation. `endpoint_positive=True` does not prove that an intermediate
negative barrier is accessible from rest; use `barrier_free_from_rest` and
`first_inaccessible_displacement_m`. Speed arrays clamp negative available
energy to zero.

For a non-neutral infinite x/y-periodic cell, the `E_bottom=0` zero mode can
leave a constant far field and a linear potential. A finite-box endpoint work
or speed remains well defined, but must not be presented as an escape energy or
terminal speed at infinity.

### 10.5 Finite-Image Shell Oracle

`finite_shell_wrench()` uses the native Fortran finite-image kernel with
`far_correction="none"`. It preserves both the raw symmetric shell and the
`E_bottom=0` result corrected by `Q_cell/(2 epsilon0 A_xy) e_z`; `selected`
points to the requested closure. The source representation is native,
canonical, and unwrapped. Do not replicate images again in Python or wrap the
moved target into the primary cell.

```python
from beach import finite_shell_convergence, finite_shell_wrench

row = finite_shell_wrench(
    snapshot, probe, transform=None, image_layers=1, closure="e_bottom_zero"
)
shells = finite_shell_convergence(
    snapshot, probe, path.displacement_m, max_layers=12
)
```

`finite_shell_convergence()` does not select from adjacent-shell increments
alone: two consecutive small increments previously produced false convergence.
The current API uses `force_tail_proxy_N` / `work_tail_proxy_J`; for an
`infinite_physical` snapshot it also sets
`reference_model="infinite_physical"` and requires
`reference_force_error_N`, `reference_work_error_J`, and
`reference_converged`. `increment_converged` is the combined tail-proxy and
physical-reference gate. Only two consecutive true values select the corrected
path. With `status="not_converged"`, both `selected_image_layers` and
`selected_path` are `None`.

Opt-in cached infinite-periodic oracles include (1) the `E_bottom=0` analytic
uniform non-neutral triangle plane: zero field below, `sigma/epsilon0` above,
surface PV `sigma/(2 epsilon0)`, and total object force
`Q^2/(2 epsilon0 A)`; and (2) the neutral `sigma_0 cos(kx)` `k != 0` sheet,
whose field and potential decay as `exp(-k |z-z0|)`, together with triangle-mesh
refinement. These cache-generating checks are not part of the lightweight test
path.

`ObjectForcePath.status="converged"` and
`DetachmentResult.numerically_qualified=True` only establish the requested path
integration tolerance. They do not establish physical qualification until
mesh, quadrature, FMM/cache, shell, path endpoint, charge snapshot, and seed
sensitivity have been checked.

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
| `beachx object-detachment <output_dir>` | Analyze a frozen-source wrench, detachment path, work, and speed while retaining periodic images |
| `beachx lint <config.toml>` | Check TOML / JSON Schema / BEACH constraints |
| `beachx config validate <config.toml>` | Validate a configuration file |
| `beachx model close-pack` | Generate a close-packed model |

The ordinary `beachx inspect` summary prints `potential_min` / `potential_max`
only from a precomputed `mesh_potential.csv`. When that file is absent, inspect
does not implicitly reconstruct potential in Python and omits those two lines.
A corrupt, non-finite, or wrong-length `mesh_potential.csv` is an input error,
not a condition for falling back to reconstruction.

With `--recompute-potential`, the summary uses the existing
`Beach.compute_potential` path regardless of whether a precomputed array is
available. This explicit path can be $O(N^2)$ with mesh size.
`--save-potential-mesh` and `--show` are also explicit potential-plot requests
and retain their existing ability to calculate potential without the flag. If
`--recompute-potential` and a potential plot are requested together, the current
plot API cannot consume the summary's computed array, so the two paths evaluate
independently and may duplicate the potential calculation.

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

`beach-inspect` follows the same precomputed-only ordinary-summary and explicit
`--recompute-potential` rules as `beachx inspect`.

## 13. Physical Constants

| Symbol | Value | Unit | Description |
|---|---|---|---|
| `K_COULOMB` | `8.9875517923e9` | N m^2 / C^2 | Coulomb constant |

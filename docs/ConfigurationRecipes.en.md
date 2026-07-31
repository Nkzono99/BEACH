title: Design a Simulation Case

Lang: [English](ConfigurationRecipes.en.md) | [日本語](ConfigurationRecipes.md)

# Design a Simulation Case

This task guide starts from the official tutorial case and replaces one mesh, particle source, or boundary setting at a time.
See the [input parameter reference](Parameters.en.html) for every key and
[Create and Validate `beach.toml`](Configuration.en.html) for file generation and validation.

## Common procedure

**Prerequisite:** Install BEACH and create a working directory.

**Action:**

```bash
beachx config init beach.toml
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

**Expected output:** `lint` accepts the configuration, and the run creates
`outputs/latest/summary.txt` and `charges.csv`.

**Interpretation:** Successful completion proves that the configuration and execution path work. It does not establish
numerical convergence or physical validity; use [Validate Results](ValidationGuide.en.html) for those checks.

**Next choices:**

| Goal | Main section to change |
| --- | --- |
| Use built-in geometry | `[mesh]`, `[[mesh.templates]]` |
| Use OBJ geometry | `[mesh]` |
| Select initial particles, reservoir inflow, or photoelectrons | `[[particles.species]]` |
| Select box, field, and particle boundaries | `[domain]`, `[field_boundary]`, `[particle_boundary]`, `[reservoir]` |
| Save time histories | `[output]` |
| Continue from a checkpoint | `[output]` |

Make one type of change at a time and rerun `beachx lint beach.toml` after each change.

## Add built-in meshes

**Prerequisite:** Use `[mesh].mode="template"`. Each `[[mesh.templates]]` entry receives a separate `mesh_id`.

**Action:** Add a shape after the existing template.

| `kind` | Shape | Dimensions | Resolution |
| --- | --- | --- | --- |
| `plane` | XY rectangle | `size_x`, `size_y` | `nx`, `ny` |
| `plate_hole`, `plane_hole` | Rectangle with a circular hole | `size_x`, `size_y`, `radius` | `n_theta`, `n_r` |
| `disk` | Disk | `radius` | `n_theta`, `n_r` |
| `annulus` | Annulus | `radius`, `inner_radius` | `n_theta`, `n_r` |
| `box` | Closed rectangular surface | `size` | `nx`, `ny`, `nz` |
| `cylinder` | z-axis cylinder | `radius`, `height`, `cap` | `n_theta`, `n_z` |
| `sphere` | Sphere | `radius` | `n_lon`, `n_lat` |

```toml
[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
enabled = true
surface_model = "insulator"
size_x = 1.0
size_y = 1.0
nx = 20
ny = 20
center = [0.5, 0.5, 0.02]

[[mesh.templates]]
kind = "sphere"
enabled = true
surface_model = "insulator"
center = [0.5, 0.5, 0.55]
radius = 0.12
n_lon = 24
n_lat = 12
```

**Expected output:** Both `mesh_id` values appear in `mesh_triangles.csv`, and `mesh_sources.csv` identifies their templates.

**Interpretation:** Increasing the resolution also increases field-evaluation and collision costs. Check placement and collisions
with a coarse mesh first. Use `surface_model="insulator"` for the standard charging model.
`conductor` is supported only with `field_boundary.mode="free"`, and `dielectric` `epsilon_r` is currently metadata; polarization is not solved.

**Next choice:** See the
[built-in geometry reference](Parameters.en.html#meshtemplates-built-in-geometries) for element counts and constraints.

## Replace the mesh with an OBJ file

**Prerequisite:** Put an OBJ file containing triangular faces where the run can read it.

**Action:**

```toml
[mesh]
mode = "obj"
obj_path = "mesh/object.obj"
surface_model = "insulator"
obj_scale = 1.0
obj_rotation = [0.0, 0.0, 0.0]
obj_offset = [0.0, 0.0, 0.0]
```

**Expected output:** The complete OBJ file appears as one `mesh_id` in `mesh_triangles.csv`.

**Interpretation:** Transformations are applied in scale → rotation → offset order.

**Next choice:** Prefer multiple built-in templates when independent objects must have separate `mesh_id` values.

## Select a particle source

**Prerequisite:** Select one `source_mode` per `[[particles.species]]` entry. When changing modes, replace the complete entry
so that keys specific to another mode do not remain.

| `source_mode` | Use | Main required keys |
| --- | --- | --- |
| `volume_seed` | Small tests with initial particles in the box | `npcls_per_step`, `pos_low`, `pos_high` |
| `plane_source` | One-way flux from an internal rectangle | `pos_low`, `pos_high`, `source_normal`, and density or grid flux |
| `reservoir_face` (deprecated) | Legacy local reservoir plane | `inject_face`, `pos_low`, `pos_high` |
| `photo_raycast` | Raycast photoelectron emission from surfaces | `rays_per_batch`, `emit_current_density_a_m2`, `ray_direction` |

Add inflow from an external reservoir with `[particles.species.boundary_inflow]`, not another `source_mode`.

### `volume_seed`

```toml
[[particles.species]]
source_mode = "volume_seed"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
npcls_per_step = 100
w_particle = 1.0
pos_low = [0.1, 0.1, 0.6]
pos_high = [0.9, 0.9, 0.9]
drift_velocity = [0.0, 0.0, -1.0e6]
temperature_k = 0.0
```

**Expected output:** Each batch creates `npcls_per_step` macro-particles.

**Interpretation:** This source does not derive its count from a physical face flux.

### Boundary-reservoir inflow

**Prerequisite:** Configure `[domain]` and a positive `sim.batch_duration`. Select the inflow barrier in `[reservoir]`.

```toml
[reservoir]
inflow_model = "source_vdf"
phi_infty = 0.0
face_potential_grid_n = 3
```

```toml
[[particles.species]]
source_mode = "volume_seed"
npcls_per_step = 0
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
target_macro_particles_per_batch = 300
drift_velocity = [0.0, 0.0, -4.0e5]

[particles.species.boundary_inflow]
z_high = "reservoir"
```

**Expected output:** External-reservoir flux determines particle weight, and particles enter through the complete z-high face.

**Interpretation:** `target_macro_particles_per_batch` fixes the computational particle count per batch. Use `w_particle`
instead to specify the weight directly. Do not set both `temperature_k` and `temperature_ev`.

**Next choice:** See
[reservoir inflow through simulation boundaries](ReservoirInjection.en.html) for flux, weight, velocity distributions, and
`infinity_barrier`.

### `plane_source`

**Prerequisite:** Configure `[domain]`, a positive `sim.batch_duration`, and an axis-aligned rectangle inside the box.

```toml
[[particles.species]]
source_mode = "plane_source"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
target_macro_particles_per_batch = 300
pos_low = [0.2, 0.2, 2.0]
pos_high = [0.8, 0.8, 2.0]
source_normal = [0.0, 0.0, -1.0]
```

**Expected output:** One-way flux is generated in the -z direction from the rectangle at z=2.0.

**Interpretation:** Because this is not an external boundary,
`reservoir.inflow_model="infinity_barrier"` and `phi_infty` do not apply. Legacy
`source_mode = "reservoir_face"` remains compatible, but new cases should select `boundary_inflow` or `plane_source` by intent.

### `photo_raycast`

**Prerequisite:** Set a positive `sim.batch_duration`.

```toml
[[particles.species]]
source_mode = "photo_raycast"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
temperature_ev = 2.2
emit_current_density_a_m2 = 4.5e-6
rays_per_batch = 1000
deposit_opposite_charge_on_emit = true
inject_face = "z_high"
inject_region_mode = "face_fraction"
uv_low = [0.0, 0.0]
uv_high = [1.0, 1.0]
ray_direction = [0.0, 0.0, -1.0]
```

**Expected output:** Rays that first hit a mesh generate photoelectrons.

**Interpretation:** `rays_per_batch` is the number of illumination rays; the generated particle count depends on the hit rate.
`deposit_opposite_charge_on_emit=true` leaves opposite charge on the emitting surface.

**Next choice:** See [Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for emission, reabsorption,
and the closed-PE `neutral_return` closure.

For closed PE with a per-species reflection on the injection face:

```toml
[particles.species.boundary]
z_high = "reflect"
```

This table belongs to the preceding `[[particles.species]]`. With
`surface_charge_closure="neutral_return"`, the effective `inject_face` action must be `reflect` or
`redistributed_reflect`. The example uses `reflect`, which preserves tangential position. See
[Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for uniform in-plane redistribution of return positions.

## Use a two-axis periodic boundary

**Prerequisite:** Configure box geometry and x/y periodicity in `[domain]`, and use `field_solver="fmm"`.
Select periodicity only through `domain.periodic_axes`.

**Action:**

```toml
[sim]
field_solver = "fmm"

[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 4.0]
periodic_axes = ["x", "y"]

[field_boundary]
mode = "periodic2"

[particle_boundary]
z_low = "open"
z_high = "open"
ordinary_open_model = "escape"
```

**Expected output:** The field, collision, and raycast paths share x/y periodic topology, while the z faces use the selected particle actions.

**Interpretation:** `[particle_boundary]` cannot select `periodic`; `[domain]` determines periodic-face behavior.

**Next choice:** Select the operator in `[periodic2]`. Use
[Finite periodic2 configuration](FinitePeriodicConfiguration.en.html) and
[`examples/periodic2_closed_photoelectron.toml`](../examples/periodic2_closed_photoelectron.toml)
for the reference boundary-reservoir + closed-PE case.

## Save histories

**Action:**

```toml
[output]
dir = "outputs/latest"
history_stride = 1
write_mesh_potential = true
write_potential_history = true
```

**Expected output:** The run creates `charge_history.csv` and `potential_history.csv`.

**Interpretation:** Potential history reevaluates the potential at every saved point and can be expensive for large meshes.

**Next choice:** Start with charge history and enable `write_potential_history` only when relative potential is required.

## Resume from a checkpoint

**Prerequisite:** The source directory must contain a valid checkpoint set.

**Action:**

```toml
[output]
dir = "outputs/continuation"
resume = true
restart_from = "../previous/outputs/latest"
```

**Expected output:** BEACH loads the checkpoint and writes new outputs to `outputs/continuation`.

**Interpretation:** `sim.batch_count` is the cumulative target. If the checkpoint has `batches=100` and the new
`batch_count=150`, the resumed run executes 50 additional batches.

**Next choice:** To continue in the same directory, omit `restart_from` and set `dir` to the checkpoint directory.

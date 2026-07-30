title: Design a Simulation Case

Lang: [English](ConfigurationRecipes.en.md) | [日本語](ConfigurationRecipes.md)

# Design a Simulation Case

Build common cases by replacing the mesh and particle sources in the official beginner `beach.toml`.
This page focuses on choosing the physical setup. For the complete key reference, see
[Input Parameters Reference](Parameters.en.html). See [Create and Validate `beach.toml`](Configuration.en.html) for generating
and checking the file.

Each snippet below is a replacement or diff relative to the [official beginner case](Tutorial.en.html).
Some snippets are not standalone configurations.

## Official Run Procedure

```bash
beachx config init
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

Pass `beach.toml` directly to the Fortran executable. Box-relative coordinate and placement parameters are converted to physical
coordinates while loading.

## Recipes

| Recipe | Main section | Use case |
| --- | --- | --- |
| Add built-in meshes | `[mesh]`, `[[mesh.templates]]` | Combine planes, spheres, boxes, cylinders, and other shapes |
| OBJ mesh | `[mesh]` | Use external geometry |
| Choose particle injection | `[[particles.species]]` | Switch inflow, initial particles, or photoemission |
| Two-periodic-axis boundary (finite image sum) | `[sim]` | Include a selected range of periodic images |
| Advanced outer-sheath coupling | `[periodic2]`, `[external_boundary]` | Couple infinite periodicity, `kinetic_1d`, and UV photoelectrons |
| History output | `[output]` | Visualize time evolution |
| Resume run | `[output]` | Continue from checkpoint |

## Add Built-In Meshes

With `[mesh].mode="template"`, each `[[mesh.templates]]` entry adds one shape. Every enabled template receives a distinct
`mesh_id`. To keep the official beginner plane and add a sphere, append a new entry after the existing `[[mesh.templates]]` entry.

The available shapes and their main size and resolution keys are:

| `kind` | Shape | Size | Resolution |
| --- | --- | --- | --- |
| `plane` | XY rectangle | `size_x`, `size_y` | `nx`, `ny` |
| `plate_hole`, `plane_hole` | Rectangle with a circular hole | `size_x`, `size_y`, `radius` | `n_theta`, `n_r` |
| `disk` | Disk | `radius` | `n_theta`, `n_r` |
| `annulus` | Concentric ring | `radius`, `inner_radius` | `n_theta`, `n_r` |
| `box` | Closed rectangular-box surface | `size = [sx, sy, sz]` | `nx`, `ny`, `nz` |
| `cylinder` | z-axis cylinder | `radius`, `height`, `cap` | `n_theta`, `n_z` |
| `sphere` | Sphere surface | `radius` | `n_lon`, `n_lat` |

Plane example:

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
```

Append a sphere to this plane:

```toml
[[mesh.templates]]
kind = "sphere"
enabled = true
surface_model = "insulator"
center = [0.5, 0.5, 0.55]
radius = 0.12
n_lon = 24
n_lat = 12
```

`center` is the shape center. Increasing subdivision counts raises the element count and the cost of field evaluation and
collision detection. Start with a coarse mesh, confirm placement and collisions, then refine it. See the
[built-in shape reference](Parameters.en.html#meshtemplates-built-in-shapes) for element counts and constraints.

Use `surface_model="insulator"` for the ordinary charging workflow. `conductor` supports only `field_bc_mode="free"` and cannot
be combined with `periodic2`. For `dielectric`, `epsilon_r` is currently metadata; BEACH does not yet solve dielectric polarization.

### Use an OBJ File

Replace `[mesh]` with OBJ mode for geometry that the built-in shapes cannot represent.

```toml
[mesh]
mode = "obj"
obj_path = "mesh/object.obj"
surface_model = "insulator"
obj_scale = 1.0
obj_rotation = [0.0, 0.0, 0.0]
obj_offset = [0.0, 0.0, 0.0]
```

Transforms apply in scale → rotation → offset order. The whole OBJ file receives one `mesh_id`. Prefer multiple built-in
templates when separate objects need distinct identities.

## Choose Particle Injection

Each `[[particles.species]]` entry adds one species. Edit the existing entry to replace the official beginner electron source,
or append entries when, for example, electrons and ions must enter together.
When changing `source_mode`, replace the complete existing entry with the corresponding example. Keeping mode-specific keys such
as `npcls_per_step`, `w_particle`, or `number_density_*` from another mode produces a validation error.

| `source_mode` | Use case | Main keys |
| --- | --- | --- |
| `volume_seed` | Small tests with particles initially placed in the box | `npcls_per_step`, `pos_low`, `pos_high` |
| `reservoir_face` | Normal inflow from a face, using Maxwellian-style inputs | `number_density_cm3`, `temperature_ev`, `inject_face`, `target_macro_particles_per_batch` |
| `photo_raycast` | Raycast photoelectron emission from surfaces | `rays_per_batch`, `emit_current_density_a_m2`, `ray_direction` |

`volume_seed` example:

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

`volume_seed` creates `npcls_per_step` particles in every batch. It does not derive the count from a physical surface flux.

Before using `reservoir_face` or `photo_raycast`, add a positive `batch_duration` for the intended physical timescale to the
existing `[sim]` table. The following is a minimal `reservoir_face` example when `batch_duration = 1.0e-5`.

```toml
[[particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
target_macro_particles_per_batch = 300
inject_face = "z_high"
inject_region_mode = "face_fraction"
uv_low = [0.0, 0.0]
uv_high = [1.0, 1.0]
drift_velocity = [0.0, 0.0, -4.0e5]
```

`reservoir_face` requires `sim.use_box=true`, positive `sim.batch_duration`, and `inject_face`.
With `inject_region_mode="face_fraction"`, `uv_low` / `uv_high` specify the aperture as fractions of the injection face.
`target_macro_particles_per_batch` fixes the computational particle count per batch and solves the particle weight from the physical inflow.
Use `w_particle` when you want to specify the weight directly.
`temperature_k` and `temperature_ev` cannot be specified together.

`photo_raycast` example:

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

`photo_raycast` also requires positive `sim.batch_duration`. `rays_per_batch` is the number of illumination rays; the first-hit
rate determines the number of emitted particles. `deposit_opposite_charge_on_emit=true` leaves the opposite charge at the source.
See [Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for emission, reabsorption, and open-face handling.

See the [particle-source overview](ParticleSourcesBoundaries.en.html) for source selection and common post-creation processing.
See [`reservoir_face` inflow and velocity sampling](ReservoirInjection.en.html) for reservoir flux, weight, and velocity distributions.

## Use a Finite Image Sum on Two Periodic Axes

`periodic2` treats exactly two axes as periodic field boundaries. This recipe sums only the primary cell and the periodic images
selected by `field_periodic_image_layers`.

```toml
[sim]
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 4.0]

bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"

field_solver = "fmm"
field_bc_mode = "periodic2"
field_periodic_image_layers = 1
field_periodic_far_correction = "none"
```

Requirements:

- `sim.use_box = true`
- exactly two periodic axes
- each periodic axis has `box_max - box_min > 0`
- `sim.field_solver = "fmm"`
- `field_periodic_image_layers >= 0`

With `field_periodic_image_layers = 1`, the field source contains the primary cell and its surrounding images, for a total of
$3\times3$ cells. `field_periodic_far_correction = "none"` does not replace more distant periodic images with an Ewald sum or a
cached operator. Increase the image layer until target quantities such as field, particle flux, and charging distribution stop
changing.

This is not the infinite-periodic solution obtained by summing the complete lattice. See
[Finite-image periodic2 configuration](FinitePeriodicConfiguration.en.html) for image-layer meaning, a local solar-wind
reservoir at z-high, photoelectron-only `neutral_return`, top-plane potential reference, and its distinction from scalar
potential barriers. The integrated runnable example is
[`examples/periodic2_closed_photoelectron.toml`](../examples/periodic2_closed_photoelectron.toml).
Use the next recipe when the case needs an infinite-periodic operator and an external sheath.

## Use the Coupled-Calculation Guide for Advanced Outer Sheaths

Infinite-periodic `periodic2`, an external `kinetic_1d` sheath, and reservoir inflow and return share potential and particle-transfer
state across several sections. Do not assemble that advanced setup from fragments on this page. Use
[Infinite-periodic periodic2 with outer plasma](InfinitePeriodicOuterConfiguration.en.html) as the canonical coupled-calculation guide.
The standard small contract fixture is
[`examples/periodic2_kinetic_outer.toml`](../examples/periodic2_kinetic_outer.toml).

To include UV photoelectrons in mean outer-sheath density, follow “Select mean outer density separately from tracked photoelectrons”
in the same coupled-calculation guide. See [Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for local
`photo_raycast` settings, source charge, and reabsorption checks. `kinetic_mean`, tracked return, and surface deposition have
different responsibilities, so this page does not present them as independent TOML fragments.

## History Output

```toml
[output]
dir = "outputs/latest"
history_stride = 1
write_mesh_potential = true
write_potential_history = true
```

`write_potential_history` evaluates potential at every history output, so it can be expensive for large meshes.
With `sim.use_box=true`, the same batches also record the z-high plane mean in
`top_reference_history.csv`. Inspect charge history first, then enable potential history when relative potential is needed.

## Resume Runs

Continue from the same output directory:

```toml
[output]
dir = "outputs/latest"
resume = true
```

Read checkpoints from a previous directory and write new outputs elsewhere:

```toml
[output]
dir = "outputs/continuation"
resume = true
restart_from = "../previous/outputs/latest"
```

`sim.batch_count` is the cumulative target batch count. If a checkpoint has `batches=100` and the new config has `batch_count=150`, only 50 additional batches run.

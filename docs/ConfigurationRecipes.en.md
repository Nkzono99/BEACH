title: Configuration Recipes

Lang: [English](ConfigurationRecipes.en.md) | [日本語](ConfigurationRecipes.md)

# Configuration Recipes

Build common cases by replacing the required tables or sections in the official beginner `beach.toml`.
For the complete key reference, see [Input Parameters Reference](Parameters.en.html). See [Edit configuration](Configuration.en.html)
for creating and checking the file.

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
| Minimal plane-mesh run | `[mesh]`, `[[mesh.templates]]` | First smoke test |
| Choose particle injection | `[[particles.species]]` | Switch inflow, initial particles, or photoemission |
| Two-periodic-axis boundary (finite image sum) | `[sim]` | Include a selected range of periodic images |
| Infinite periodic + external kinetic sheath | `[sim]`, `[periodic2]`, `[outer_plasma]`, `[coupling]` | Production lunar-regolith setup |
| Add UV photoelectrons to the outer sheath | `[outer_plasma]`, `[[particles.species]]` | Consistent emission, return, and outer space charge |
| OBJ mesh | `[mesh]` | External geometry |
| History output | `[output]` | Visualize time evolution |
| Resume run | `[output]` | Continue from checkpoint |

## Change Plane-Mesh Resolution

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

Increasing `nx` / `ny` increases the element count. Start small, inspect outputs and particle counts, then increase resolution.

## Choose Particle Injection

| `source_mode` | Use case | Main keys |
| --- | --- | --- |
| `volume_seed` | Small tests with particles initially placed in the box | `npcls_per_step`, `pos_low`, `pos_high` |
| `reservoir_face` | Normal inflow from a face, using Maxwellian-style inputs | `number_density_cm3`, `temperature_ev`, `inject_face`, `target_macro_particles_per_batch` |
| `photo_raycast` | Raycast photoelectron emission from surfaces | `rays_per_batch`, `emit_current_density_a_m2`, `ray_direction` |

Minimal `reservoir_face` example:

```toml
[[particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
target_macro_particles_per_batch = 300
inject_face = "z_high"
pos_low = [0.0, 0.0, 4.0]
pos_high = [1.0, 1.0, 4.0]
drift_velocity = [0.0, 0.0, -4.0e5]
```

`target_macro_particles_per_batch` fixes the computational particle count per batch and solves the particle weight from the physical inflow.
Use `w_particle` when you want to specify the weight directly.

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
[Finite-image periodic2 configuration](FinitePeriodicConfiguration.en.html) for image-layer meaning, scalar potential barriers,
and convergence checks. Use the next recipe when the case needs an infinite-periodic operator and an external sheath.

## Infinite Periodic + External Kinetic Sheath

Use this `kinetic_1d` composition as the standard starting point for an outer sheath in BEACH. Consider
[advanced rough-surface linear screening](UnifiedLinearResponse.en.html) only for specialized cases that must screen lateral
fields near a rough surface linearly.

The production lunar-regolith setup solves the `k != 0` surface-charge component with the infinite-periodic cache and
connects the area-mean `k = 0` component to an external 1D kinetic sheath. Add these main sections to an existing box,
mesh, and ambient electron/ion configuration:

```toml
[sim]
b0 = [0.0, 0.0, 0.0]
field_solver = "fmm"
field_bc_mode = "periodic2"
field_periodic_far_correction = "cached_kneq0"
field_periodic_cache_dir = ".beach_cache/periodic2"
reservoir_potential_model = "none"
sheath_injection_model = "none"

[periodic2]
nonzero_mode_backend = "cached_kneq0"
zero_mode_policy = "exclude_k0"
lower_boundary_model = "symmetric_vacuum"

[outer_plasma]
model = "kinetic_1d"
photoelectron_density_model = "none"
return_model = "kinetic_1d_profile_return"
interface_z = 9.899494936611664e-4
infinity_potential = 0.0
debye_length = 10.5132
thermal_voltage = 10.0

[coupling]
update_mode = "explicit"
particle_transfer_mode = "electrostatic_1d_instant_return"
outer_update_stride = 1
field_evolution_timescale = 6.060915267313266e-8
max_frozen_field_ratio = 0.1
outer_queue_enabled = false
```

Set `interface_z` to the z-high box face. Negative and positive `reservoir_face` species provide the infinity ambient
electron and ion VDFs, so define both at z-high and satisfy quasineutrality and the ion Bohm-entry condition.
`infinity_potential = 0` fixes the infinity-potential gauge. Inflow acceleration and outgoing-particle escape/return
use the same kinetic-profile difference `phi_interface - phi_infinity`.

Do not copy the example values for `debye_length`, `thermal_voltage`, or `field_evolution_timescale`; derive them for
the target plasma and timescale. Validate first with `outer_update_stride = 1`. Before increasing it to a production
value such as `100`, confirm that surface potential, absorbed/escaped flux, charge balance, and detachment force are
unchanged. Unsupported non-monotonic branches, sub-Bohm ions, and frozen-field-limit violations stop without falling
back to another model. See `examples/periodic2_kinetic_outer.toml` for the complete small fixture and
[Input Parameters Reference](Parameters.en.html) for the full contract.

## Add UV Photoelectrons to the Outer Sheath

For UV emission, solve the mean photoelectron density in the external region with `kinetic_mean`, while tracked
`photo_raycast` particles update surface charge. Change and add the following entries:

```toml
[outer_plasma]
model = "kinetic_1d"
photoelectron_density_model = "kinetic_mean"
return_model = "kinetic_1d_profile_return"
interface_z = 9.899494936611664e-4
infinity_potential = 0.0
debye_length = 10.5132
thermal_voltage = 10.0

[[particles.species]]
enabled = true
source_mode = "photo_raycast"
emit_current_density_a_m2 = 4.5e-6
rays_per_batch = 5000
deposit_opposite_charge_on_emit = true
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
temperature_ev = 2.2
normal_drift_speed = 0.0
inject_face = "z_high"
ray_direction = [0.0, 0.0, -1.0]
```

The first negative `photo_raycast` species supplies the emission flux and temperature for the mean closure.
`deposit_opposite_charge_on_emit = true` is required.
`kinetic_mean` supplies only the outer profile and does not add the return current to the surface a second time.
Compare UV-off and UV-on runs with the same mesh, batch duration, and ambient inflow. Inspect
`outer_plasma_profile.csv`, solver residual and species currents in `summary.txt`, and the charge ledger.

`sim.sheath_injection_model = "zhao_*"` is the older inflow-distribution correction, not the external
`kinetic_1d` Poisson profile. BEACH rejects combining them or combining the kinetic model with
`reservoir_potential_model`, so both remain `"none"` in this recipe.

## OBJ Mesh

```toml
[mesh]
mode = "obj"
obj_path = "mesh/object.obj"
surface_model = "insulator"
```

OBJ mode assigns the same `surface_model` to all elements.
Use template meshes when you need separate objects with different surface models, or check the current limitations before relying on mixed OBJ surfaces.

## History Output

```toml
[output]
dir = "outputs/latest"
history_stride = 1
write_mesh_potential = true
write_potential_history = true
```

`write_potential_history` evaluates potential at every history output, so it can be expensive for large meshes.
Inspect charge history first, then enable potential history when needed.

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

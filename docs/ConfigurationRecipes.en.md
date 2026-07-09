title: Configuration Recipes

Lang: [English](ConfigurationRecipes.en.md) | [日本語](ConfigurationRecipes.md)

# Configuration Recipes

This page shows how to modify `beach.toml` for common cases.
For the complete key reference, see [Input Parameters Reference](Parameters.en.html). For high-level notation, see [`beachx config` / High-Level Notation Guide](Configuration.en.html).

## Official Run Procedure

```bash
beachx config init
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

Pass `beach.toml` directly to the Fortran executable. The Fortran parser normalizes high-level notation while loading the file.

## Recipes

| Recipe | Main section | Use case |
| --- | --- | --- |
| Minimal plane-mesh run | `[mesh]`, `[[mesh.templates]]` | First smoke test |
| Choose particle injection | `[[particles.species]]` | Switch inflow, initial particles, or photoemission |
| Two-periodic-axis boundary | `[sim]` | Periodic cell simulations |
| OBJ mesh | `[mesh]` | External geometry |
| History output | `[output]` | Visualize time evolution |
| Resume run | `[output]` | Continue from checkpoint |

## Minimal Plane-Mesh Run

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

## Two-Periodic-Axis Boundary

`periodic2` treats exactly two axes as periodic field boundaries. In the current implementation it requires `field_solver = "fmm"`.

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

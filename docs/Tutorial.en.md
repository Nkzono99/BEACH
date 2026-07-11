title: Ten-Minute Tutorial

Lang: [English](Tutorial.en.md) | [日本語](Tutorial.md)

# Ten-Minute Tutorial

The official beginner case tracks one electron into an insulating plane and verifies absorption and charge deposition.
It intentionally excludes FMM, periodic boundaries, photoelectrons, conductors, and dielectric models.

## 1. Create the configuration

```bash
mkdir beach-tutorial
cd beach-tutorial
beachx config init beach.toml
beachx lint beach.toml
```

The generated file is identical to
[`examples/tutorial_insulator.toml`](https://github.com/Nkzono99/BEACH/blob/main/examples/tutorial_insulator.toml).
The complete configuration is:

```toml
[sim]
dt = 1.0e-7
batch_count = 1
max_step = 10
softening = 1.0e-6
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 1.0]
bc_x_low = "reflect"
bc_x_high = "reflect"
bc_y_low = "reflect"
bc_y_high = "reflect"
bc_z_low = "open"
bc_z_high = "open"
rng_seed = 12345
open_boundary_model = "escape"
field_solver = "direct"
field_bc_mode = "free"
field_periodic_far_correction = "none"

[particles]
[[particles.species]]
source_mode = "volume_seed"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
w_particle = 1.0
npcls_per_step = 1
pos_low = [0.5, 0.5, 0.8]
pos_high = [0.5, 0.5, 0.8]
drift_velocity = [0.0, 0.0, -1.0e6]
temperature_k = 0.0

[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
enabled = true
surface_model = "insulator"
size_x = 1.0
size_y = 1.0
nx = 4
ny = 4
center = [0.5, 0.5, 0.2]

[output]
write_files = true
dir = "outputs/latest"
history_stride = 1
```

## 2. Run

```bash
beach beach.toml
beachx inspect outputs/latest
```

A successful run creates `outputs/latest/summary.txt` and `charges.csv`. This deterministic case should report
`batches=1` and `processed_particles=1`, with the particle absorbed by the plane.

```text
processed_particles=1
absorbed=1
batches=1
survived_max_step=0
```

## 3. Inspect the result

```bash
beachx inspect outputs/latest --save-mesh outputs/latest/charge.png
```

Only the hit element receives negative charge. This verifies the execution path; it is not a steady-state charging result.

![Surface charge density from the official beginner case](media/tutorial_insulator_charge.png)

## 4. First parameters to change

| Goal | Keys |
| --- | --- |
| More particles | `npcls_per_step` |
| Broader launch region | `pos_low`, `pos_high` |
| Different velocity | `drift_velocity`, `temperature_k` |
| Finer surface | `nx`, `ny` |
| Longer run | `batch_count`, `max_step` |
| Thinner history | `history_stride` |

Continue with [Configuration Recipes](ConfigurationRecipes.en.html), then complete
[Validating Simulation Results](ValidationGuide.en.html) before quantitative use.

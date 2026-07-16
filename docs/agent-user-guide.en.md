title: BEACH Agent User Guide

Lang: [English](agent-user-guide.en.md) | [日本語](agent-user-guide.md)

# BEACH Agent User Guide

> Reference guide for AI agents operating BEACH simulations.
> Intended to be loaded from CLAUDE.md with `@import docs/agent-user-guide.en.md`.
> [BEACH Documentation](index.en.html) is the entry point for normal use and execution instructions.

---

## Overview

BEACH (BEM + Accumulated CHarge) is a hybrid boundary-element-method plus particle-tracking simulator for charge accumulation on insulator surfaces.

- **Fortran core**: particle dynamics, electric-field solver, collision detection, charge deposition
- **Python layer**: configuration management, post-processing, visualization
- **Version**: 1.6.2

---

## Quick Start

### Minimal Run Procedure

```bash
# 1. Install
pip install beach-bem

# 2. Create a configuration file
beachx config init

# 3. Validate the configuration
beachx lint beach.toml

# 4. Run the simulation
beach beach.toml

# 5. Inspect results
beachx inspect outputs/latest
```

## Build and Run

### Requirements

| Item | Requirement |
|------|------|
| Fortran compiler | gfortran, or compatible |
| fpm | v0.10+ |
| Python | 3.10+ |
| Main Python dependencies | matplotlib >= 3.8, numpy >= 1.24 |

### Build Commands

```bash
make check                     # development build check with fixed dev version
make build                     # build with git describe version
make run CONFIG=beach.toml     # run with fixed dev version
make install-generic           # portable gfortran build
make install-camphor           # Intel compiler optimized build
make test                      # L1: Python + quick Fortran tests
```

`make check` / `make test` / `make run` use `BEACH_VERSION_MODE=dev` and keep the version macro passed to Fortran stable. Because the fpm compile-flag hash does not change when the git hash changes, incremental compilation during development is easier to reuse. When an executable with the git hash is required, use `make build VERSION_MODE=git` or `make install`.

### Low-level Direct fpm Run

```bash
fpm run --profile release --flag "-fopenmp" -- beach.toml
```

### MPI Parallel Run

```bash
FPM_FC=mpiifx fpm run --profile release \
  --flag "-fpp -DUSE_MPI -qopenmp" \
  --runner "mpirun -n 4" -- beach.toml
```

### Tests

```bash
make test-l0      # L0: static/schema/build check
make test         # L1: normal development loop
make test-l2      # L2: contract/integration
make test-l3      # L3: cumulative L0-L3 verification
make test-heavy   # heavy Fortran targets only
make test-fortran-far-correction  # oracle far-correction correctness
make test-fortran-far-correction-diagnostics  # assertion-free diagnostics
make test-fortran-benchmark  # release-profile runtime benchmark
make test-field-kernel-cache  # opt-in native cache/plane-oracle receipt gate
make test-full    # unfiltered fpm test
make test-mpi     # MPI テスト
pytest -q         # Python テストのみ
```

`make test` is the L1 alias and is the normal default for the inner AI/development loop.
Long-running FMM targets are run explicitly with `make test-l3` / `make test-heavy` / `make test-fortran-heavy` / `make test-full`. Use `make test-fortran-far-correction` for `m2l_root_oracle` correctness, `make test-fortran-far-correction-diagnostics` for assertion-free output, and `make test-fortran-benchmark` for release-profile timing.
Run the shared-kernel cache contract and native periodic plane-oracle receipt
with the opt-in `make test-field-kernel-cache` target. It passes the built
library's absolute path to the test and is not included in L1/L2/L3 or
`make test-physics-release`.
Individual targets can be checked with `FPM_ACTION=test ./build.sh --target <name>`.

---

## Configuration Parameters

### Configuration Workflow

1. **beach.toml**: the normal file to edit, read directly by the Fortran executable
2. **beachx lint**: validates TOML parsing, JSON Schema, coordinate and placement parameters, and known constraints
3. **Fortran parser**: converts box-relative settings into physical coordinates such as `box_min` / `box_max` / `center`

### [sim] Section - Simulation Basics

| Parameter | Type | Default | Description |
|------------|------|-----------|------|
| `dt` | float | 1.0e-9 | Time step [s] |
| `rng_seed` | int | 12345 | Random seed |
| `batch_count` | int | 1 | Number of batches in normal runs. During resume, cumulative target batch count |
| `batch_duration` | float | - | Batch duration [s], mutually exclusive with `batch_duration_step` |
| `batch_duration_step` | float | - | Resolved as `batch_duration = dt * batch_duration_step` |
| `max_step` | int | 400 | Maximum integration steps per particle |
| `tol_rel` | float | 1.0e-8 | Monitoring metric, not an early-stop condition |
| `q_floor` | float | 1.0e-30 | Denominator floor for relative-change calculation |
| `softening` | float | 1.0e-6 | Electric-field softening length [m] |

### Electric-field Solver Parameters

| Parameter | Type | Default | Choices | Description |
|------------|------|-----------|--------|------|
| `field_solver` | string | "auto" | direct, treecode, fmm, auto | Electric-field evaluation method |
| `field_bc_mode` | string | "free" | free, periodic2 | Normal periodic2 uses FMM; the Direct split reference is the exception |
| `field_periodic_image_layers` | int | 1 | >= 0 | Number of periodic2 image-shell layers |
| `field_periodic_far_correction` | string | "none" | auto, none, m2l_root_oracle, cached_kneq0 | `cached_kneq0` is the production infinite-periodic nonzero mode |
| `field_periodic_ewald_alpha` | float | 0.0 | >= 0 | Ewald splitting parameter, 0 means automatic |
| `field_periodic_ewald_layers` | int | 4 | >= 0 | Ewald shell depth |
| `field_periodic_cache_dir` | string | ".beach_cache/periodic2" | non-empty | Versioned operator cache |
| `field_periodic_generation_tolerance` | float | 1e-8 | > 0 | Generation tolerance in the cache identity |
| `tree_theta` | float | 0.5 | (0, 1] | Tree-method MAC parameter |
| `tree_leaf_max` | int | 16 | >= 1 | Maximum number of elements per leaf node |
| `tree_min_nelem` | int | 256 | >= 1 | Threshold for auto -> treecode switch |

**auto estimation table** (when tree_theta / tree_leaf_max are unspecified):

| Element count | theta | leaf_max |
|--------|-------|----------|
| < 1500 | 0.40 | 12 |
| 1500-9999 | 0.50 | 16 |
| 10000-49999 | 0.58 | 20 |
| 50000+ | 0.65 | 24 |

### Magnetic Field and Potential

| Parameter | Type | Default | Description |
|------------|------|-----------|------|
| `b0` | float[3] | [0, 0, 0] | Uniform magnetic field [T] |
| `reservoir_potential_model` | string | "none" | none, infinity_barrier |
| `phi_infty` | float | 0.0 | Infinity reference potential [V] |
| `injection_face_phi_grid_n` | int | 3 | Injection-face potential grid resolution NxN |
| `raycast_max_bounce` | int | 16 | Maximum number of ray reflections for photo_raycast |

### Sheath Model

| Parameter | Type | Default | Description |
|------------|------|-----------|------|
| `sheath_injection_model` | string | "none" | none, zhao_auto, zhao_a, zhao_b, zhao_c, floating_no_photo |
| `sheath_alpha_deg` | float | 60.0 | Solar elevation angle [deg] |
| `sheath_photoelectron_ref_density_cm3` | float | 64.0 | Reference photoelectron density [cm^-3] |
| `sheath_reference_coordinate` | float | - | Sheath reference-plane coordinate [m] |
| `sheath_electron_drift_mode` | string | "normal" | normal, full |
| `sheath_ion_drift_mode` | string | "normal" | normal, full |

### Box Boundary Conditions

| Parameter | Type | Default | Description |
|------------|------|-----------|------|
| `use_box` | bool | false | Enable box boundary |
| `box_min` | float[3] | [-1, -1, -1] | Lower limit [m] |
| `box_max` | float[3] | [1, 1, 1] | Upper limit [m] |
| `bc_{x,y,z}_{low,high}` | string | "open" | open, reflect, periodic |

**periodic2 constraint**: exactly two axes are periodic, and the remaining one axis is open/reflect.

### [[particles.species]] Section - Particle Species (at least one required)

#### Common Parameters

| Parameter | Type | Default | Description |
|------------|------|-----------|------|
| `enabled` | bool | true | Enable this species |
| `source_mode` | string | "volume_seed" | volume_seed, reservoir_face, photo_raycast |
| `q_particle` | float | -1.602e-19 | Particle charge [C] |
| `m_particle` | float | 9.109e-31 | Particle mass [kg] |
| `pos_low` | float[3] | [-0.4, -0.4, 0.2] | Lower generation-position limit [m] |
| `pos_high` | float[3] | [0.4, 0.4, 0.5] | Upper generation-position limit [m] |
| `drift_velocity` | float[3] | [0, 0, -8e5] | Drift velocity [m/s] |
| `temperature_k` | float | 20000 | Temperature [K], mutually exclusive with `temperature_ev` |
| `temperature_ev` | float | - | Temperature [eV] |

#### volume_seed Mode

| Parameter | Type | Default | Description |
|------------|------|-----------|------|
| `npcls_per_step` | int | 0 | Number of particles generated per batch |
| `w_particle` | float | 1.0 | Macro-particle weight |

**Constraint**: sum of `npcls_per_step` across all species >= 1.

#### reservoir_face Mode (Physical Flux Injection)

| Parameter | Type | Default | Description |
|------------|------|-----------|------|
| `number_density_cm3` | float | - | Upstream density [cm^-3], mutually exclusive with `number_density_m3` |
| `number_density_m3` | float | - | Upstream density [m^-3] |
| `w_particle` | float | - | Macro-particle weight, mutually exclusive with `target_macro_particles_per_batch` |
| `target_macro_particles_per_batch` | int | - | Target number of macro particles per batch; -1 reuses the species[1] weight |
| `inject_face` | string | required | x_low, x_high, y_low, y_high, z_low, z_high |

**Constraint**: `use_box = true` and `batch_duration > 0` are required. `pos_low/pos_high` are placed on the specified face.

#### photo_raycast Mode (Photoelectron Emission)

| Parameter | Type | Default | Description |
|------------|------|-----------|------|
| `emit_current_density_a_m2` | float | required | Emission current density [A/m^2] |
| `rays_per_batch` | int | required | Number of rays per batch |
| `deposit_opposite_charge_on_emit` | bool | false | Deposit opposite-sign charge on the emitting element |
| `normal_drift_speed` | float | 0.0 | Normal-direction drift speed [m/s] |
| `ray_direction` | float[3] | inward normal | Ray direction vector |

A `photo_raycast` species using tracked outer transfer must set `deposit_opposite_charge_on_emit=true`, regardless of whether
the histogram is enabled.

### [mesh] Section - Geometry

| Parameter | Type | Default | Description |
|------------|------|-----------|------|
| `mode` | string | "auto" | auto, obj, template |
| `obj_path` | string | examples/simple_plate.obj | OBJ file path |
| `surface_model` | string | "insulator" | Surface model for the whole OBJ mesh (`insulator`, `conductor`, `dielectric`) |
| `epsilon_r` | float | 1.0 | Relative permittivity for the whole OBJ mesh (`>= 1`) |
| `obj_scale` | float | 1.0 | Scale factor |
| `obj_rotation` | float[3] | [0, 0, 0] | Rotation angles [deg], extrinsic x->y->z |
| `obj_offset` | float[3] | [0, 0, 0] | Translation [m] |

**Transform order**: scale -> rotate -> offset

#### [[mesh.templates]] - Procedural Mesh Generation

Common: `enabled` (bool), `kind` (enum), `surface_model` (enum), `epsilon_r` (float), `center` (float[3])

`conductor` is redistributed as an equipotential floating conductor per mesh_id. The current implementation supports only `sim.field_bc_mode = "free"`. OBJ input reads the entire file as `mesh_id = 1`, so separated conductor parts in a single OBJ are also treated as the same floating conductor. To treat them as independent conductors, split their mesh_id values, for example by using template input. `dielectric` is currently metadata that retains per-object `epsilon_r`; physical branches for dielectric polarization are an extension point.

| kind | Main parameters |
|------|---------------|
| `plane` | `size_x`, `size_y`, `nx`, `ny` |
| `plate_hole` | `size_x`, `size_y`, `radius`, `n_theta`, `n_r` |
| `disk` | `radius`, `n_theta`, `n_r` |
| `annulus` | `radius`, `inner_radius`, `n_theta`, `n_r` |
| `box` | `size` (float[3]), `nx`, `ny`, `nz` |
| `cylinder` | `radius`, `height`, `n_theta`, `n_z`, `cap`, `cap_top`, `cap_bottom` |
| `sphere` | `radius`, `n_lon`, `n_lat` |

### [output] Section - File Output

| Parameter | Type | Default | Description |
|------------|------|-----------|------|
| `write_files` | bool | true | Enable file output |
| `write_mesh_potential` | bool | false | Output element potentials to mesh_potential.csv |
| `write_potential_history` | bool | false | Output potential history |
| `dir` | string | "outputs/latest" | Output directory |
| `history_stride` | int | 1 | History output interval [batches], disabled at 0 |
| `resume` | bool | false | Resume from a checkpoint |
| `restart_from` | string | none | Checkpoint read source when `resume=true`. New output is saved to `dir` |

---

## Output File Formats

Output destination: the directory specified by `output.dir`, default `outputs/latest/`

### Required Output Files

| File | Format | Contents |
|----------|------|------|
| `summary.txt` | Text (key-value) | Run metadata and statistics |
| `charges.csv` | CSV: `elem_idx, charge_C` | Final element charge |
| `mesh_triangles.csv` | CSV: `elem_idx, v0x, v0y, v0z, v1x, v1y, v1z, v2x, v2y, v2z, charge_C, mesh_id` | Triangle vertices, charge, and mesh_id |
| `mesh_sources.csv` | CSV: `mesh_id, source_kind, template_kind, surface_model, epsilon_r, elem_count` | Mesh source metadata |
| `rng_state.txt` | Text | Random state, for resume |
| `charge_ledger.csv` | CSV | Per-species charge balance and restartable cumulative state |

### Optional Output Files

| File | Condition | Format |
|----------|------|------|
| `charge_history.csv` | `history_stride > 0` | CSV: `batch, elem_idx, charge_C` |
| `potential_history.csv` | `write_potential_history = true` and `history_stride > 0` | CSV: `batch, elem_idx, potential_V` |
| `mesh_potential.csv` | `write_mesh_potential = true` | CSV: `elem_idx, potential_V` |
| `macro_residuals.csv` | When reservoir_face is used | CSV: injection residual state |
| `outer_plasma_profile.csv` | A ready `kinetic_1d` / `unified_linear_response` outer state | CSV: outer profile and conditional checkpoint |
| `photoelectron_histogram.csv` | `photoelectron_histogram_enabled=true` | CSV: previous-batch and cumulative histogram; a conditional checkpoint |
| `performance_profile.csv` | When the `BEACH_PROFILE=1` environment variable is set | CSV: measured times for each region |

### Additional Files During MPI Runs

- `rng_state_rank00000.txt`, `rng_state_rank00001.txt`, ...
- `macro_residuals.csv` remains global; the root writes one file rather than one per rank

The canonical conditions and resume requirements are in the [Output Guide](OutputGuide.en.html#resume-outputs).

---

## Python CLI (beachx)

### config - Configuration Management

```bash
beachx lint [beach.toml]                           # check schema and semantic constraints
beachx config init [beach.toml]                    # create a new beach.toml
beachx config validate [beach.toml]                # validate coordinate/placement parameters and constraints
beachx config diff left.toml right.toml            # compare configurations
```

### inspect - Result Inspection

```bash
beachx inspect [output_dir]                        # print summary
beachx inspect outputs/latest --show               # display plots
beachx inspect outputs/latest --save-bar charges.png
beachx inspect outputs/latest --save-mesh charges_mesh.png
beachx inspect outputs/latest --save-potential-mesh potential_mesh.png
beachx inspect outputs/latest --apply-periodic2-mesh
beachx inspect outputs/latest --periodic2-repeat 1
```

### animate - Create Animations

```bash
beachx animate [output_dir]                        # 履歴アニメーション
beachx animate outputs/latest --quantity charge     # 電荷 (デフォルト)
beachx animate outputs/latest --quantity potential   # ポテンシャル
beachx animate outputs/latest --save-gif charge.gif
beachx animate outputs/latest --total-frames 200
beachx animate outputs/latest --frame-stride 2 --fps 15
```

### slices - Potential Sections

```bash
beachx slices [output_dir]
beachx slices outputs/latest --grid-n 200
beachx slices outputs/latest --xy-z 0.5 --yz-x 0.5 --xz-y 0.5
beachx slices outputs/latest --vmin -20 --vmax 20
beachx slices outputs/latest --save slices.png
```

### coulomb - Coulomb Force Matrix

```bash
beachx coulomb [output_dir]
beachx coulomb outputs/latest --component z
beachx coulomb outputs/latest --target-kinds sphere
beachx coulomb outputs/latest --save forces.png
```

### mobility - Coulomb Mobility Analysis

```bash
beachx mobility [output_dir]
beachx mobility outputs/latest --density-kg-m3 2500
beachx mobility outputs/latest --mu-static 0.4
beachx mobility outputs/latest --save-csv mobility.csv
```

### workload - Workload Estimation

```bash
beachx workload beach.toml
beachx workload beach.toml --threads 8
beachx workload beach.toml --mpi-ranks 4 --mpi-rank 0 \
  --macro-residuals outputs/latest/macro_residuals.csv
```

With MPI options, `total_particles` is the estimate for the selected rank.
`global_total_particles` is the all-rank total, while `local_reservoir_particles` and
`global_reservoir_particles` show reservoir injection after and before rank distribution.
The reservoir residual is updated once from the global expectation, so the global sequence is independent of MPI size.

---

## Python API

```python
from beach import Beach

run = Beach("outputs/latest")

# 電荷分布の可視化
fig, ax = run.plot_charges(step=-1)
fig, ax = run.plot_mesh()
fig, ax = run.plot_mesh_source_boxplot(quantity="charge", step=-1)

# ポテンシャル解析
fig, ax = run.plot_potential()
slices = run.compute_potential_slices(grid_n=200)
fig, axes = run.plot_potential_slices(slices)

# クーロン力解析
fig, ax = run.plot_coulomb_force_matrix(component="z")
result = run.analyze_coulomb_mobility(density_kg_m3=2500, mu_static=0.4)

# 電場・電界線
field = run.compute_electric_field_points(points)
lines = run.trace_field_lines(seed_points, direction="forward")
fig = run.plot_field_lines_3d(lines)

# アニメーション
run.animate_mesh(quantity="charge", save_path="charge.gif")
```

### Main Function List

| Function | Purpose |
|------|------|
| `plot_mesh()` | 3D mesh visualization |
| `plot_potential()` | 3D potential mesh |
| `plot_charges()` | Bar chart of element charge |
| `plot_mesh_source_boxplot()` | Box plot by mesh source |
| `animate_mesh()` | Mesh animation from history |
| `compute_potential_mesh()` | Element-potential reconstruction |
| `compute_potential_points()` | Potential evaluation at arbitrary points |
| `compute_potential_slices()` | Section-data calculation |
| `compute_electric_field_points()` | Electric-field vector evaluation |
| `plot_potential_slices()` | 2D section plots |
| `plot_coulomb_force_matrix()` | Force-matrix heatmap |
| `analyze_coulomb_mobility()` | Mobility analysis |
| `trace_field_lines()` | Field-line integration |
| `plot_field_lines_3d()` | Field-line plus mesh drawing |

---

## Coordinate and placement helper parameters

`box_origin` / `box_size`, face fractions, box-relative template placement, and group scaling are ordinary configuration
parameters. See [Input Parameters Reference](Parameters.en.html#coordinate-and-placement-helper-parameters) for calculated
targets, invalid combinations, and the conditions that replace explicit dimensions.

---

## Configuration Examples

### Minimal Configuration (volume_seed)

```toml
[sim]
dt = 1.0e-8
batch_count = 10
max_step = 500
use_box = true
box_min = [0, 0, 0]
box_max = [1, 1, 1]

[[particles.species]]
source_mode = "volume_seed"
npcls_per_step = 100

[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
size_x = 1.0
size_y = 1.0
nx = 10
ny = 10

[output]
dir = "outputs/test"
```

### Periodic Boundary + volume seed

```toml
[sim]
batch_count = 200
field_solver = "fmm"
field_bc_mode = "periodic2"
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 10.0]
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"

[[particles.species]]
source_mode = "volume_seed"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
npcls_per_step = 10

[[particles.species]]
source_mode = "volume_seed"
q_particle = 1.602176634e-19
m_particle = 1.672482821616e-27
npcls_per_step = 10
```

### Complete Post-processing Pipeline

```bash
beach beach.toml
beachx inspect outputs/latest --save-mesh charges.png
beachx animate outputs/latest --quantity charge --save-gif charge.gif
beachx slices outputs/latest --save slices.png
beachx coulomb outputs/latest --save forces.png
beachx mobility outputs/latest --density-kg-m3 2500 --mu-static 0.4
```

---

## Parameter Constraint Checklist

| Feature | Required conditions |
|------|---------|
| `reservoir_face` | `use_box=true`, `batch_duration>0`, `inject_face` specified |
| `photo_raycast` | `use_box=true`, `batch_duration>0`, `emit_current_density_a_m2>0`, `rays_per_batch>=1` |
| `periodic2` | Normally `field_solver=fmm`; the Direct split reference is the exception. Exactly two periodic axes and `use_box=true` |
| Sheath model | Compatible with `reservoir_potential_model = "none"` |
| Resume | `write_files=true`, checkpoint file exists, in the `restart_from` directory when specified, MPI size matches |
| Performance profiling | Environment variable `BEACH_PROFILE=1` |
| MPI run | Compiled with `-DUSE_MPI`, MPI compiler wrapper used |

The canonical solver, kernel, and boundary table is [Field evaluation](FieldSolvers.en.html#solver-and-field-boundary-compatibility).

---

## Project Structure

```
BEACH/
├── app/main.f90              # Fortran entry point
├── src/                      # Fortran library modules
│   ├── config/               #   configuration parser
│   ├── core/                 #   simulator main loop
│   ├── mesh/                 #   triangle mesh management
│   ├── particles/            #   particle dynamics (Boris pusher)
│   ├── physics/              #   electric-field solvers (direct, treecode, FMM)
│   └── runtime/              #   output and restart
├── beach/                    # Python package
│   ├── cli/                  #   CLI subcommands
│   ├── config/               #   configuration validation
│   └── fortran_results/      #   post-processing and visualization
├── schemas/                  # JSON Schema for validation
│   └── beach.schema.json     #   beach.toml schema
├── examples/                 # sample configurations and scripts
├── tests/                    # test suite
├── docs/                     # documentation
├── fpm.toml                  # Fortran package manifest
├── pyproject.toml            # Python package metadata
├── Makefile                  # build automation
├── SPEC.md                   # Fortran implementation specification
└── AGENTS.md                 # agent configuration
```

---

## Detailed Documentation References

| Document | Contents |
|-------------|------|
| `SPEC.md` | Fortran implementation specification, authoritative |
| `docs/OutputGuide.en.md` | How to read output files |
| `docs/ConfigurationRecipes.en.md` | Common configuration recipes |
| `docs/Parameters.en.md` | Detailed parameter specification |
| `docs/Execution.en.md` | Normal runs, workload estimation, and resume |
| `docs/Workflow.en.md` | Development, testing, and HPC operation |
| `docs/Algorithms.en.md` | Computational model and batch-loop overview |
| `docs/SurfaceModels.en.md` | Batch charge commit and surface models |
| `docs/ParticleSourcesBoundaries.en.md` | Particle-source overview |
| `docs/ReservoirInjection.en.md` | Reservoir flux, velocity distributions, and potential mapping |
| `docs/PhotoelectronEmission.en.md` | Photoelectron emission and lifecycle |
| `docs/SheathInjectionClosures.en.md` | Zhao and floating source-VDF corrections |
| `docs/ParticleTrackingCollision.en.md` | Particle-update flow |
| `docs/BorisPusher.en.md` | Boris velocity and position update |
| `docs/ParticleEvents.en.md` | Triangle collisions and box-boundary events |
| `docs/FieldSolvers.en.md` | Field-evaluation selection and common settings |
| `docs/DirectSolver.en.md` | Direct field and potential evaluation |
| `docs/Treecode.en.md` | Treecode construction, accuracy, and constraints |
| `docs/PeriodicElectrostatics.en.md` | periodic2 field and zero mode |
| `docs/FinitePeriodicConfiguration.en.md` | Integrated finite-image and scalar-boundary configuration |
| `docs/InfinitePeriodicOuterConfiguration.en.md` | Integrated infinite-periodic and outer-plasma configuration |
| `docs/OuterPlasmaModels.en.md` | Choose the standard outer sheath or advanced rough-surface screening |
| `docs/KineticOuterPlasma.en.md` | Standard and recommended self-consistent kinetic 1-D outer sheath |
| `docs/UnifiedLinearResponse.en.md` | Advanced rough-surface linear screening |
| `docs/ParticleEscapeReturn.en.md` | Open boundaries, 1-D return, and 3-D outer orbits |
| `docs/FMM.en.md` | FMM selection and verification |
| `docs/FMMCore.en.md` | FMM internals and Ewald |
| `docs/BatchDurationStability.en.md` | `batch_duration` stability |
| `docs/Configuration.en.md` | Configuration creation, lint, and schema |
| `docs/PostprocessTutorial.en.md` | Post-processing tutorial |
| `docs/PythonPostprocessAPI.en.md` | Python API reference |
| `schemas/beach.schema.json` | JSON Schema for IDE validation |

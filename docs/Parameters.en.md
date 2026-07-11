title: BEACH Input Parameters Reference

Lang: [English](Parameters.en.md) | [日本語](Parameters.md)

# Input Parameters Reference

This document is the parameter reference for the final `beach.toml` read by the
Fortran runtime. Unless otherwise noted, units are SI units.

For first-time configuration work, start with
[Configuration Recipes](ConfigurationRecipes.en.html).

High-level notation normalized by the Fortran parser is summarized in
[Configuration](Configuration.en.html). This document focuses on runtime keys
used after normalization.

| Related document | Contents |
|---|---|
| [Configuration Recipes](ConfigurationRecipes.en.html) | Task-oriented configuration steps and tuning points |
| [Configuration](Configuration.en.html) | `beachx config`, high-level notation, schema/lint |
| [Algorithms](Algorithms.en.html) | Entry point to BEM field calculation, particle push, collision, and accumulated-charge procedure |
| [Workflow](Workflow.en.html) | Execution, development, testing, and KUDPC notes |
| [FMMCore](FMMCore.en.html) | Numerical algorithm for `field_solver="fmm"` |

---

## Table of Contents

- [Input Parameters Reference](#input-parameters-reference)
  - [Table of Contents](#table-of-contents)
  - [Loading Rules](#loading-rules)
  - [Units and Coordinates](#units-and-coordinates)
  - [Minimal Example](#minimal-example)
  - [Section List](#section-list)
  - [Detailed Parameter Reference](#detailed-parameter-reference)
    - [`[sim]`: Run Control and Field Calculation](#sim-run-control-and-field-calculation)
    - [`[[particles.species]]`: Particle Species](#particlesspecies-particle-species)
    - [`sim.sheath_injection_model`: Sheath Injection Correction](#simsheath_injection_model-sheath-injection-correction)
    - [`[mesh]`: Mesh Input](#mesh-mesh-input)
    - [`[[mesh.templates]]`: Built-In Shapes](#meshtemplates-built-in-shapes)
    - [`[output]`: Output and Resume](#output-output-and-resume)
  - [Relationship to High-Level Notation](#relationship-to-high-level-notation)
  - [Validation Rules](#validation-rules)

---

## Loading Rules

| Item | Specification |
|---|---|
| Explicit input | `beach path/to/beach.toml` |
| Default input | With no argument, `beach.toml` in the current directory |
| Developer run | `fpm run -- path/to/beach.toml` behaves the same |
| Format | TOML. Multi-line arrays are supported |
| Unknown keys | Unknown section names and key names are errors |
| schema | `schemas/beach.schema.json` |
| lint | `beachx lint beach.toml` |

To use an editor schema, put a comment directive at the beginning of
`beach.toml`. The Fortran parser does not accept regular keys before the first
section, so do not use `"$schema" = "..."`.

```toml
#:schema ../schemas/beach.schema.json
```

You can also specify the GitHub Raw URL.

```toml
#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.schema.json
```

Relative paths are interpreted relative to that `beach.toml` file. If the file
is placed in a deeper directory such as `outputs/.../beach.toml`, adjust the path
for example to `../../schemas/beach.schema.json`.

---

## Units and Coordinates

| Kind | Representative keys | Unit / direction |
|---|---|---|
| Time | `dt`, `batch_duration` | seconds |
| Length | `box_min`, `box_max`, `pos_low`, `pos_high` | m |
| Charge | `q_particle`, element charge output | C |
| Mass | `m_particle` | kg |
| Velocity | `drift_velocity`, `ray_direction` | m/s. `ray_direction` is a direction vector |
| Electric field | `e0`, `e0_abs` | V/m |
| Magnetic field | `b0` | T |
| Density | `number_density_cm3`, `number_density_m3` | cm^-3 or m^-3 |
| Temperature | `temperature_k`, `temperature_ev` | K or eV. They cannot both be specified |
| Angle | `e0_phi_xy_deg`, `e0_phi_z_deg`, `sheath_alpha_deg` | degree |

`*_low` / `*_high` are lower and upper bounds on each axis. `inject_face` is one
of `x_low`, `x_high`, `y_low`, `y_high`, `z_low`, or `z_high`.

---

## Minimal Example

For first-time use, `source_mode = "reservoir_face"` is recommended because it
lets you specify physical inflow directly.

```toml
[sim]
dt = 2.0e-8
batch_duration_step = 60000.0
batch_count = 100
max_step = 10000
softening = 1.0e-6
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 10.0]
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
rng_seed = 12345
field_solver = "fmm"
field_bc_mode = "periodic2"
field_periodic_far_correction = "none"

[[particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
target_macro_particles_per_batch = 5000
inject_face = "z_high"
pos_low = [0.0, 0.0, 10.0]
pos_high = [1.0, 1.0, 10.0]
drift_velocity = [0.0, 0.0, -4.0e5]

[[particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = 1.602176634e-19
m_particle = 1.672482821616e-27
target_macro_particles_per_batch = -1
inject_face = "z_high"
pos_low = [0.0, 0.0, 10.0]
pos_high = [1.0, 1.0, 10.0]
drift_velocity = [0.0, 0.0, -4.0e5]

[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
enabled = true
size_x = 1.0
size_y = 1.0
nx = 20
ny = 20
center = [0.5, 0.5, 0.02]

[output]
write_files = true
write_mesh_potential = false
dir = "outputs/latest"
history_stride = 1
```

---

## Section List

| Section | Required | Contents |
|---|---:|---|
| `[sim]` | conditional | Time step, batch count, field solver, boundaries, external fields, sheath correction |
| `[particles]` | yes | Container for `[[particles.species]]` |
| `[[particles.species]]` | yes | Species, injection mode, velocity distribution, macro-particle weight |
| `[mesh]` | no | Selection of OBJ or built-in template input |
| `[[mesh.templates]]` | no | Built-in shapes used with `mode="template"` |
| `[output]` | no | Output directory, history, checkpoint resume |

`[sim]` is required when using `reservoir_face` or `photo_raycast`.
At least one `[[particles.species]]` entry is required.

---

## Detailed Parameter Reference

### `[sim]`: Run Control and Field Calculation

#### Run Control

| Key | Type | Default | Description |
|---|---|---:|---|
| `dt` | float | `1.0e-9` | Time step [s] |
| `rng_seed` | int | `12345` | Random seed |
| `batch_count` | int | `1` | Number of batches processed in a normal run. With `output.resume=true`, this is the cumulative target batch count |
| `batch_duration` | float | `0.0` | Physical time per batch [s] |
| `batch_duration_step` | float | `0.0` | Resolved as `batch_duration = dt * batch_duration_step` |
| `max_step` | int | `400` | Maximum number of pushes per particle |
| `tol_rel` | float | `1.0e-8` | Monitored relative-change value. Not used as a stop condition |
| `q_floor` | float | `1.0e-30` | Lower bound for the denominator in `rel_change` calculations |

Specifying both `batch_duration` and `batch_duration_step` is an error. For
`reservoir_face` / `photo_raycast`, the resolved `batch_duration > 0` is required.

#### Field Solver

`field_solver` selects the method for computing the Coulomb electric field at
evaluation points from boundary-element charges. The table below lists the
corresponding parameters for each option.

| `field_solver` | Use case | Supported field boundary |
|---|---|---|
| `direct` | Exact all-to-all evaluation for small element counts | `field_bc_mode="free"` |
| `treecode` | Approximate evaluation for medium and larger cases | `field_bc_mode="free"` |
| `fmm` | Large-scale evaluation, `periodic2`, FMM core validation | `field_bc_mode="free"` / `"periodic2"` |
| `auto` | Automatically choose direct / treecode based on element count | `field_bc_mode="free"` |

Common keys:

| Key | Type | Default | Description |
|---|---|---:|---|
| `softening` | float | `1.0e-6` | Softening length for Coulomb field calculation [m] |
| `field_solver` | string | `"auto"` | `direct` / `treecode` / `fmm` / `auto` |
| `field_normalization` | string | `"si"` | `si` / `box` / `mesh` / `length` |
| `field_length_scale` | float | `1.0` | Length scale used with `field_normalization="length"` [m] |

`field_normalization` only changes normalization of coordinates, softening, and
periodic cells inside field calculations. Output electric fields and potentials
are converted back to SI.

| `field_normalization` | Length scale |
|---|---|
| `si` | Use input SI coordinates as-is |
| `box` | Maximum width of `sim.box_max - sim.box_min`. Requires `sim.use_box=true` |
| `mesh` | Maximum width of the mesh bounding box. Falls back to `field_length_scale` if the mesh is empty |
| `length` | `field_length_scale` |

##### `field_solver = "direct"`

Sums all source elements directly. There is no approximation error, and the
computational complexity is `O(MN)`, where `M` is the number of evaluation points
and `N` is the number of elements.

| Key | Type | Default | Description |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | Specify `"direct"` |
| `softening` | float | `1.0e-6` | Relaxes singularity when a source and evaluation point are close |
| `field_normalization` | string | `"si"` | Normalize coordinates before direct evaluation |
| `field_length_scale` | float | `1.0` | Used with `field_normalization="length"` or mesh fallback |
| `field_bc_mode` | string | `"free"` | Only `"free"` is supported for `direct` |

`tree_theta`, `tree_leaf_max`, and `tree_min_nelem` are not used for `direct`.

##### `field_solver = "treecode"`

Builds a source octree. Distant nodes are evaluated with a multipole
approximation, and near nodes are evaluated with a direct sum. Unlike FMM, it
does not use local expansion and traverses the tree for each evaluation point.

| Key | Type | Default | Description |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | Specify `"treecode"` |
| `softening` | float | `1.0e-6` | Softening for near direct sums and multipole evaluation |
| `field_normalization` | string | `"si"` | Normalize coordinates before tree construction |
| `field_length_scale` | float | `1.0` | Used with `field_normalization="length"` or mesh fallback |
| `tree_theta` | float | `0.5` | MAC parameter. `0 < theta <= 1`. Larger values are faster and coarser |
| `tree_leaf_max` | int | `16` | Maximum number of sources per leaf node. `>= 1` |
| `field_bc_mode` | string | `"free"` | Only `"free"` is supported for `treecode` |

`tree_min_nelem` is a threshold for `field_solver="auto"`, so it does not switch
anything when `treecode` is specified explicitly.

##### `field_solver = "fmm"`

Uses the simulator-independent Coulomb FMM core. It separates the source
geometry plan from the state updated for each charge update, and evaluates with
P2M/M2M/M2L/L2L/L2P plus near direct sums. See
[Coulomb FMM Core Details](FMMCore.en.html) for details.

| Key | Type | Default | Description |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | Specify `"fmm"` |
| `softening` | float | `1.0e-6` | Softening for near direct sums and FMM evaluation |
| `field_normalization` | string | `"si"` | Normalize coordinates before FMM plan construction |
| `field_length_scale` | float | `1.0` | Used with `field_normalization="length"` or mesh fallback |
| `tree_theta` | float | `0.5` | MAC parameter for near/far classification. `0 < theta <= 1` |
| `tree_leaf_max` | int | `16` | Maximum number of sources per source-tree leaf node. `>= 1` |
| `field_bc_mode` | string | `"free"` | `free` / `periodic2` |
| `field_periodic_image_layers` | int | `1` | Number of near image layers for `periodic2` |
| `field_periodic_far_correction` | string | `"none"` | Far correction for `periodic2`. `none` / `auto` / `m2l_root_oracle` |
| `field_periodic_ewald_alpha` | float | `0.0` | Ewald decomposition parameter for `m2l_root_oracle` |
| `field_periodic_ewald_layers` | int | `4` | Real-space outer shell / reciprocal cutoff depth for the Ewald oracle |

`field_periodic_*` only has meaning when `field_bc_mode="periodic2"`.
`tree_min_nelem` is not used when `fmm` is specified explicitly.

##### `field_solver = "auto"`

Uses direct evaluation when the element count is less than `tree_min_nelem`, and
treecode otherwise. `auto` does not switch to FMM.

| Key | Type | Default | Description |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | Specify `"auto"` |
| `softening` | float | `1.0e-6` | Softening for direct / treecode |
| `field_normalization` | string | `"si"` | Common normalization used before automatic selection |
| `field_length_scale` | float | `1.0` | Used with `field_normalization="length"` or mesh fallback |
| `tree_min_nelem` | int | `256` | Element-count threshold for switching to treecode. `>= 1` |
| `tree_theta` | float | `0.5` | MAC parameter when treecode is selected |
| `tree_leaf_max` | int | `16` | Maximum number of sources per leaf node when treecode is selected |
| `field_bc_mode` | string | `"free"` | Only `"free"` is supported for `auto` |

If `tree_theta` and `tree_leaf_max` are not specified explicitly, the following
values are used based on the element count.

| Element count `nelem` | `tree_theta` | `tree_leaf_max` |
|---:|---:|---:|
| `< 1500` | `0.40` | `12` |
| `1500 <= nelem < 10000` | `0.50` | `16` |
| `10000 <= nelem < 50000` | `0.58` | `20` |
| `50000 <= nelem` | `0.65` | `24` |

#### Field Boundary and periodic2

| Key | Type | Default | Description |
|---|---|---:|---|
| `field_bc_mode` | string | `"free"` | `free` / `periodic2` |
| `field_periodic_image_layers` | int | `1` | Number of near image layers for `periodic2` |
| `field_periodic_far_correction` | string | `"none"` | `none` / `auto` / `m2l_root_oracle` |
| `field_periodic_ewald_alpha` | float | `0.0` | Ewald decomposition parameter for `m2l_root_oracle` |
| `field_periodic_ewald_layers` | int | `4` | Outer shell / reciprocal cutoff depth for the Ewald oracle |

`periodic2` is available only with `field_solver="fmm"`. It is valid when
`sim.use_box=true`, exactly two axes are `periodic`, and the remaining axis is
open. Periodic images are considered not only for field evaluation, but also for
collision and `photo_raycast` raycasting.

`field_periodic_far_correction="auto"` is accepted for compatibility and
currently behaves the same as `none`. `m2l_root_oracle` is a diagnostic build-time
mode that fits the exact periodic Ewald residual to the root operator, and is
high cost.

#### External Fields

| Key | Type | Default | Description |
|---|---|---:|---|
| `e0` | float[3] | `[0,0,0]` | Uniform external electric field [V/m] |
| `e0_abs` | float | unspecified | Magnitude of the uniform external electric field [V/m] |
| `e0_phi_xy_deg` | float | `0.0` | Azimuthal angle in the xy plane when `e0_abs` is specified [deg] |
| `e0_phi_z_deg` | float | `0.0` | Elevation angle from the xy plane when `e0_abs` is specified [deg] |
| `b0` | float[3] | `[0,0,0]` | Uniform magnetic field [T] |

The uniform external electric field can be specified directly as
`e0 = [Ex, Ey, Ez]`, or by `e0_abs` and angles. Mixing the two forms is an error.

#### Inflow Helpers and Sheath Correction

| Key | Type | Default | Description |
|---|---|---:|---|
| `reservoir_potential_model` | string | `"none"` | `none` / `infinity_barrier` |
| `phi_infty` | float | `0.0` | Reference potential at infinity [V] |
| `open_boundary_model` | string | `"escape"` | `escape` / `potential_barrier` |
| `injection_face_phi_grid_n` | int | `3` | `N x N` evaluation grid for injection-face average potential |
| `raycast_max_bounce` | int | `16` | Maximum number of reflections for `photo_raycast` |
| `sheath_injection_model` | string | `"none"` | `none` / `zhao_auto` / `zhao_a` / `zhao_b` / `zhao_c` / `floating_no_photo` |
| `sheath_alpha_deg` | float | `60.0` | Solar elevation angle for the Zhao sheath [deg] |
| `sheath_photoelectron_ref_density_cm3` | float | `64.0` | Reference photoelectron density for the Zhao sheath [cm^-3] |
| `sheath_reference_coordinate` | float | unspecified | Reference plane position for the sheath 1D coordinate [m] |
| `sheath_electron_drift_mode` | string | `"normal"` | `normal` / `full` |
| `sheath_ion_drift_mode` | string | `"normal"` | `normal` / `full` |

Currently, `sheath_injection_model != "none"` is used together with
`reservoir_potential_model="none"`. See
[`sim.sheath_injection_model`](#simsheath_injection_model-sheath-injection-correction)
for details.

#### Computational Domain and Particle Boundaries

| Key | Type | Default | Description |
|---|---|---:|---|
| `use_box` | bool | `false` | Enable box boundaries |
| `box_min` | float[3] | `[-1,-1,-1]` | Lower box bounds [m] |
| `box_max` | float[3] | `[1,1,1]` | Upper box bounds [m] |
| `bc_x_low`, `bc_x_high` | string | `"open"` | Particle boundaries at the lower and upper x limits |
| `bc_y_low`, `bc_y_high` | string | `"open"` | Particle boundaries at the lower and upper y limits |
| `bc_z_low`, `bc_z_high` | string | `"open"` | Particle boundaries at the lower and upper z limits |

Particle boundaries are `open`, `reflect`, or `periodic`. `open` also accepts
`outflow` and `escape` as synonyms.

With `open_boundary_model="potential_barrier"`, particles crossing an `open`
face are tested against a local BEM potential barrier. The barrier is
`q_particle * (phi_infty - phi_boundary)`, where `phi_boundary` is evaluated at
the boundary crossing point. If the barrier is positive and exceeds the normal
kinetic energy `0.5 * m_particle * v_normal^2`, the normal velocity is reflected;
otherwise the particle escapes. The uniform external field `e0` is not included
in this barrier because it does not define a finite potential relative to
infinity.

For `periodic2`, the mesh is translated at runtime to a canonical unwrapped
representation for collision before ray-triangle tests. Raw vertices may lie
outside the box along periodic axes, but triangles are not folded modulo the box
vertex by vertex.

---

### `[[particles.species]]`: Particle Species

At least one `[[particles.species]]` entry is required. The keys and constraints
used depend on `source_mode`.

#### Common Keys

| Key | Type | Default | Description |
|---|---|---:|---|
| `enabled` | bool | `true` | Enable the species |
| `source_mode` | string | `"volume_seed"` | `volume_seed` / `reservoir_face` / `photo_raycast` |
| `q_particle` | float | `-1.602176634e-19` | Particle charge [C] |
| `m_particle` | float | `9.10938356e-31` | Particle mass [kg] |
| `pos_low` | float[3] | `[-0.4,-0.4,0.2]` | Lower position bounds [m] |
| `pos_high` | float[3] | `[0.4,0.4,0.5]` | Upper position bounds [m] |
| `drift_velocity` | float[3] | `[0,0,-8e5]` | Drift velocity [m/s] |
| `temperature_k` | float | `2.0e4` | Temperature [K] |
| `temperature_ev` | float | unspecified | Temperature [eV]. Mutually exclusive with `temperature_k` |
| `velocity_distribution` | string | `"maxwellian"` | `maxwellian` / `grid` |
| `inject_face` | string | unspecified | Injection face. Required for `reservoir_face` / `photo_raycast` |

#### `source_mode = "volume_seed"`

| Key | Type | Default | Description |
|---|---|---:|---|
| `npcls_per_step` | int | `0` | Number of macro particles generated per batch |
| `w_particle` | float | `1.0` | Macro-particle weight |

Constraints:

| Condition | Details |
|---|---|
| Particle count | The sum of `npcls_per_step` over all enabled species must be at least 1 |
| Automatic weight resolution | `target_macro_particles_per_batch` cannot be used |

#### `source_mode = "reservoir_face"`

| Key | Type | Description |
|---|---|---|
| `number_density_cm3`, `number_density_m3` | float | Upstream density. Specify exactly one of them |
| `w_particle` | float | Macro-particle weight. Positive value |
| `target_macro_particles_per_batch` | int | For automatic `w_particle` resolution. `>0` or `-1` |
| `velocity_grid_path` | string | CSV path for `velocity_distribution="grid"` |
| `velocity_grid_pdf_kind` | string | `phase_space` / `flux_weighted` |
| `velocity_grid_sampling` | string | `auto` / `rectilinear` / `discrete` |
| `particle_flux_m2_s`, `current_density_a_m2` | float | Inflow amount for the `grid` distribution. Specify exactly one of them |

Basic constraints:

| Condition | Details |
|---|---|
| Domain | `sim.use_box=true` is required |
| Time | `sim.batch_duration > 0` is required |
| Injection face | `inject_face` is required |
| Injection range | `pos_low` / `pos_high` must lie on the specified face |
| Weight | `w_particle` and `target_macro_particles_per_batch` cannot both be specified |
| Weight sharing | `target_macro_particles_per_batch=-1` is allowed only for species 2 and later. It shares species 1 `w_particle` |

For the Maxwellian distribution, the one-sided flux of a drifting Maxwellian is
computed from `number_density_*` and temperature.

For the grid distribution, the CSV at `velocity_grid_path` is read. Required
columns are `vx_m_s, vy_m_s, vz_m_s, f`. `f` is normalized internally so that
`sum f = 1`. In this case, `number_density_*` / `temperature_*` are not used;
the inflow amount is set by `particle_flux_m2_s` or `current_density_a_m2`.
`current_density_a_m2` is converted to particle flux as `abs(J / q_particle)`.

| `velocity_grid_sampling` | Behavior |
|---|---|
| `auto` | Trilinear interpolation for a complete rectilinear grid. Discrete sampling for incomplete grids or scattered points |
| `rectilinear` | Force rectilinear-grid interpolation. Error if the grid is not rectilinear |
| `discrete` | Sample CSV rows directly |

| `velocity_grid_pdf_kind` | Behavior |
|---|---|
| `phase_space` | Sample with `v_n f(v)`, multiplying by inward normal velocity `v_n` through the inflow face |
| `flux_weighted` | Treat CSV `f` as an already flux-weighted distribution |

For either PDF, only velocities with `v_n > 0` are used. The relative path in
`velocity_grid_path` is based on the runtime current directory. Currently,
`velocity_distribution="grid"` is valid only when
`sim.sheath_injection_model="none"`.

The particle count is determined as follows.

```text
n_macro_expected = gamma_in * A * batch_duration / w_particle
n_injected = floor(residual + n_macro_expected)
```

The residual is carried over to the next batch. When
`target_macro_particles_per_batch > 0`, `w_particle` is calculated automatically
to approach that value.

#### `source_mode = "photo_raycast"`

| Key | Type | Default | Description |
|---|---|---:|---|
| `emit_current_density_a_m2` | float | `0.0` | Emission current density referenced to the ray-normal plane [A/m^2] |
| `rays_per_batch` | int | `0` | Number of launched rays per batch |
| `deposit_opposite_charge_on_emit` | bool | `false` | Deposit opposite-sign charge on the emitting element |
| `photo_escape_model` | string | `"none"` | `none` / `boltzmann_cutoff` |
| `normal_drift_speed` | float | `0.0` | Drift speed along the emission normal [m/s] |
| `ray_direction` | float[3] | inward normal of injection face | Ray direction |

Constraints:

| Condition | Details |
|---|---|
| Domain | `sim.use_box=true` is required |
| Time | `sim.batch_duration > 0` is required |
| Emission amount | `emit_current_density_a_m2 > 0`, `rays_per_batch > 0` are required |
| Injection face | `inject_face` is required |
| Particle properties | `q_particle` is nonzero, and `m_particle > 0` |
| Ray direction | Must be normalizable, and its dot product with the inward normal of the injection face must be positive |
| Unavailable keys | `npcls_per_step`, `number_density_*`, `w_particle`, `target_macro_particles_per_batch` |

The weight when one ray hits is:

```text
w_hit = J_perp * A_perp * batch_duration / (|q_particle| * rays_per_batch)
```

The actual number of generated particles depends on the hit rate, so the number
generated per batch is at most `rays_per_batch`. With `field_bc_mode="periodic2"`,
emission starts from the hit coordinate wrapped to the primary cell even when a
periodic image is hit.

With `photo_escape_model="boltzmann_cutoff"`, the barrier is evaluated using the
center potential of the emitting element excluding its self contribution.

```text
barrier = max(phi_emit - phi_infty, 0)
escape_factor = exp(-|q_particle| * barrier / (k_B * T_PE))
```

When `deposit_opposite_charge_on_emit=true`, the same effective weight is used
for the opposite-sign charge left on the emitting element. Suppressed
photoelectron current is treated as an immediate return.

---

### `sim.sheath_injection_model`: Sheath Injection Correction

`sim.sheath_injection_model` is a shared setting that groups existing
`reservoir_face` / `photo_raycast` species and overrides sheath-aware fluxes and
normal-velocity cutoffs.

| Value | Details |
|---|---|
| `none` | No correction |
| `zhao_auto` | Automatically choose Zhao Type A/B/C branches based on solar elevation angle |
| `zhao_a`, `zhao_b`, `zhao_c` | Use Zhao 1D photoelectron sheath conditions with the specified branch |
| `floating_no_photo` | Simple floating sheath without photoelectrons |

For Zhao models, the following species are detected automatically.

| Target | Detection condition |
|---|---|
| solar-wind electron | First negative-charge `reservoir_face` species |
| ion | First positive-charge `reservoir_face` species |
| photoelectron | First negative-charge `photo_raycast` species |

Effects of Zhao models:

| Target | Override |
|---|---|
| electron reservoir | Replace effective density with Zhao solution `n_swe_inf` and apply `vmin_normal` according to the barrier |
| ion reservoir | If `sheath_reference_coordinate` is specified, update to local density, local normal velocity, and cold-beam approximation |
| photoelectron | Replace `emit_current_density_a_m2` with the Zhao free photoelectron current, and treat `normal_drift_speed=0` |

`floating_no_photo` solves a negative floating potential from the current balance
of the first negative-charge / positive-charge `reservoir_face` species. It
applies a cutoff to the electron reservoir species and treats emitted current as
0 even if a `photo_raycast` species exists.

Notes:

- Zhao models reuse `temperature_*`, `number_density_*`, `drift_velocity`,
  `m_particle`, and `q_particle` as background plasma conditions.
- `sheath_reference_coordinate` is the reference plane position along the normal
  axis of the shared `inject_face`.
- For example, if `inject_face="z_high"` and `sheath_reference_coordinate=0.02`,
  the plane `z=0.02` is treated as `z_sheath=0`.
- If `sheath_reference_coordinate` is not specified, only the shared cutoff-based
  correction is applied.
- In the Fortran implementation, Type A reconstructs the local profile with
  first-order integration, and Type B/C do so with first-order integration on a
  monotone branch.
- `zhao_auto` tries branch solutions in the order `C -> A -> B` for
  `alpha < 20 deg`, and `A -> B -> C` otherwise.

---

### `[mesh]`: Mesh Input

| Key | Type | Default | Description |
|---|---|---:|---|
| `mode` | string | `"auto"` | `auto` / `obj` / `template` |
| `obj_path` | string | `"examples/simple_plate.obj"` | OBJ file path |
| `surface_model` | string | `"insulator"` | Surface model for the whole OBJ |
| `epsilon_r` | float | `1.0` | Relative permittivity for the whole OBJ. `>= 1` |
| `obj_scale` | float | `1.0` | Uniform scale after loading the OBJ |
| `obj_rotation` | float[3] | `[0,0,0]` | Rotation angle after loading the OBJ [deg] |
| `obj_offset` | float[3] | `[0,0,0]` | Translation after loading the OBJ [m] |

With `mode="auto"`, an OBJ is used if `obj_path` exists; otherwise a template is
used. The OBJ transformation order is `scale -> rotate -> offset`.

```text
v_new = R(rotation) * (v_old * obj_scale) + obj_offset
```

For OBJ input, the whole file is read as `mesh_id=1`. Even if one OBJ contains
separate `conductor` parts, they are treated as the same floating conductor. To
treat them as independent conductors, split `mesh_id` values with template input
or another method.

Supported OBJ input:

| Item | Supported |
|---|---|
| Newlines | LF / CRLF |
| Face lines | `f v`, `f v/vt`, `f v/vt/vn`, `f v//vn` |
| Polygons | Quadrilaterals and larger polygons are fan-triangulated |

---

### `[[mesh.templates]]`: Built-In Shapes

Common keys:

| Key | Type | Default | Description |
|---|---|---:|---|
| `enabled` | bool | `true` | Enable the template |
| `kind` | string | `"plane"` | `plane` / `plate_hole` / `plane_hole` / `disk` / `annulus` / `box` / `cylinder` / `sphere` |
| `surface_model` | string | `"insulator"` | `insulator` / `conductor` / `dielectric` |
| `surface_side` | string | unset | Vacuum side for `triangle_p0`: `normal_plus` / `normal_minus` / `outward_closed` |
| `epsilon_r` | float | `1.0` | Relative permittivity. `>= 1` |
| `center` | float[3] | `[0,0,0]` | Shape center [m] |

When `[[mesh.templates]]` entries are written, the number of templates actually
used is determined by the number of definitions. Disabled templates are not added
to the mesh and do not consume a `mesh_id`.

Overview of `kind`:

| `kind` | Generated shape | Reference plane / axis |
|---|---|---|
| `plane` | Rectangular plane | XY plane, `z=center[3]` |
| `plate_hole`, `plane_hole` | Rectangular plane with a circular center hole | XY plane, hole center is `center` |
| `disk` | Disk | XY plane, center is `center` |
| `annulus` | Concentric ring | XY plane, center is `center` |
| `box` | Closed rectangular-box surface | Six faces parallel to the axes |
| `cylinder` | Cylinder side surface and optional top/bottom caps | z-axis direction |
| `sphere` | Sphere surface | Center is `center` |

#### `kind = "plane"`

Splits a rectangle on the XY plane into `nx * ny` rectangular cells, then splits
each cell into two triangles.

| Key | Type | Default | Description |
|---|---|---:|---|
| `center` | float[3] | `[0,0,0]` | Plane center `[x, y, z]` [m] |
| `size_x` | float | `1.0` | Size in the x direction [m]. `> 0` |
| `size_y` | float | `1.0` | Size in the y direction [m]. `> 0` |
| `nx` | int | `1` | Number of divisions in the x direction. `>= 1` |
| `ny` | int | `1` | Number of divisions in the y direction. `>= 1` |

The number of elements is `2 * nx * ny`.

#### `kind = "plate_hole"` / `"plane_hole"`

A rectangular plate on the XY plane with a circular center hole removed.
`plane_hole` is an alias for `plate_hole`. The hole boundary is approximated by a
polygon with `n_theta` divisions, and the region from the hole edge to the outer
boundary is divided into `n_r` layers.

| Key | Type | Default | Description |
|---|---|---:|---|
| `center` | float[3] | `[0,0,0]` | Plate center and hole center `[x, y, z]` [m] |
| `size_x` | float | `1.0` | Size in the x direction [m]. `> 0` |
| `size_y` | float | `1.0` | Size in the y direction [m]. `> 0` |
| `radius` | float | `0.5` | Hole radius [m]. At runtime, `0 < radius < min(size_x, size_y) / 2` |
| `n_theta` | int | `24` | Circumferential divisions of the hole boundary. `>= 3` |
| `n_r` | int | `4` | Radial divisions from the hole edge to the outer boundary. `>= 1` |

The outer boundary matches the rectangular boundary. It is an error if the
circular hole radius is at least the half-width or half-height. The common
default `radius=0.5` hits this constraint with the default
`size_x=size_y=1.0`, so specify `radius` explicitly for `plate_hole`.

#### `kind = "disk"`

A disk on the XY plane. The interior is divided in polar coordinates and
triangulated from the center toward the outer boundary.

| Key | Type | Default | Description |
|---|---|---:|---|
| `center` | float[3] | `[0,0,0]` | Disk center `[x, y, z]` [m] |
| `radius` | float | `0.5` | Disk radius [m]. `> 0` |
| `n_theta` | int | `24` | Circumferential divisions. `>= 3` |
| `n_r` | int | `4` | Radial divisions. `>= 1` |

Internally this uses the same generation path as `annulus` with
`inner_radius=0`.

#### `kind = "annulus"`

A concentric ring on the XY plane. The region from inner radius to outer radius
is divided into `n_r` layers.

| Key | Type | Default | Description |
|---|---|---:|---|
| `center` | float[3] | `[0,0,0]` | Ring center `[x, y, z]` [m] |
| `radius` | float | `0.5` | Outer radius [m]. `> 0` |
| `inner_radius` | float | `0.25` | Inner radius [m]. `0 <= inner_radius < radius` |
| `n_theta` | int | `24` | Circumferential divisions. `>= 3` |
| `n_r` | int | `4` | Radial divisions. `>= 1` |

`inner_radius=0` is accepted, but `kind="disk"` is clearer when creating a disk.

#### `kind = "box"`

A closed rectangular-box surface. All six faces are triangulated, and vertex
ordering is set so normals point outward.

| Key | Type | Default | Description |
|---|---|---:|---|
| `center` | float[3] | `[0,0,0]` | Box center `[x, y, z]` [m] |
| `size` | float[3] | `[1,1,1]` | Size in the x, y, and z directions [m]. Each component `> 0` |
| `nx` | int | `1` | Number of divisions in the x direction. `>= 1` |
| `ny` | int | `1` | Number of divisions in the y direction. `>= 1` |
| `nz` | int | `1` | Number of divisions in the z direction. `>= 1` |

The number of elements is `4 * (nx * ny + ny * nz + nx * nz)`. This counts two
opposite faces for each axis pair after splitting each rectangular face cell
into two triangles.

#### `kind = "cylinder"`

A cylinder along the z axis. The side surface is split into `n_theta * n_z`
rectangular cells, and top/bottom caps are added as needed.

| Key | Type | Default | Description |
|---|---|---:|---|
| `center` | float[3] | `[0,0,0]` | Cylinder center `[x, y, z]` [m] |
| `radius` | float | `0.5` | Cylinder radius [m]. `> 0` |
| `height` | float | `1.0` | Height in the z direction [m]. `> 0` |
| `n_theta` | int | `24` | Circumferential divisions. `>= 3` |
| `n_z` | int | `1` | Axial divisions. `>= 1` |
| `cap` | bool | `true` | Enable both top and bottom caps together |
| `cap_top` | bool | value of `cap` | Top cap. Overrides `cap` when specified |
| `cap_bottom` | bool | value of `cap` | Bottom cap. Overrides `cap` when specified |

The cylinder extends from `z = center[3] - height/2` to
`z = center[3] + height/2`. The side surface has `2 * n_theta * n_z` elements.
Each enabled cap adds `n_theta` triangles.

#### `kind = "sphere"`

A sphere surface based on longitude and latitude divisions.

| Key | Type | Default | Description |
|---|---|---:|---|
| `center` | float[3] | `[0,0,0]` | Sphere center `[x, y, z]` [m] |
| `radius` | float | `0.5` | Sphere radius [m]. `> 0` |
| `n_lon` | int | `24` | Longitude divisions. `>= 3` |
| `n_lat` | int | `12` | Latitude divisions. `>= 2` |

The number of elements is `2 * n_lon * (n_lat - 1)`. Near the poles, each cell
uses one triangle; other latitude bands use two triangles.

Surface model:

| `surface_model` | Behavior |
|---|---|
| `insulator` | Accumulate the charge of colliding particles on elements |
| `conductor` | Redistribute element charge for each `mesh_id` floating conductor so it becomes equipotential while conserving total charge |
| `dielectric` | Store `epsilon_r` as metadata. Current field calculation and charge accumulation do not yet branch for dielectric polarization |

At present, `conductor` redistributes charge using direct Coulomb coefficients
with `field_bc_mode="free"`. It cannot be combined with
`field_bc_mode="periodic2"`. For cases with many conductor elements, the dense
solve increases additional per-batch cost.

---

### `[output]`: Output and Resume

| Key | Type | Default | Description |
|---|---|---:|---|
| `write_files` | bool | `true` | Enable/disable file output |
| `write_mesh_potential` | bool | `false` | Output `mesh_potential.csv` |
| `write_potential_history` | bool | `false` | Output `potential_history.csv` |
| `dir` | string | `"outputs/latest"` | Output directory |
| `history_stride` | int | `1` | Output interval for history CSV [batch] |
| `resume` | bool | `false` | Resume from an existing checkpoint |
| `restart_from` | string | none | Checkpoint source when `resume=true` |

Output files:

| File | Condition / contents |
|---|---|
| `summary.txt` | Run statistics and configuration summary |
| `charges.csv` | Final element charges |
| `mesh_triangles.csv` | Element geometry. Includes the `mesh_id` column |
| `mesh_sources.csv` | Original mesh kind, surface model, `epsilon_r`, and element count for each `mesh_id` |
| `mesh_potential.csv` | When `write_mesh_potential=true` |
| `charge_history.csv` | When `history_stride > 0` |
| `potential_history.csv` | When `write_potential_history=true` and `history_stride > 0` |
| `rng_state.txt` | Random-number state |
| `macro_residuals.csv` | Residual carry-over for macro-particle counts |
| `charge_ledger.csv` | Per-species signed-charge flux, counts, and restartable cumulative values |

`mesh_potential.csv` records the potential [V] at each element centroid. The
self term uses `1/softening` when `softening > 0`; otherwise it uses an
area-equivalent radius approximation. For `periodic2`, the explicit image shell
is added. The exact Ewald residual is also added only when
`field_periodic_far_correction="m2l_root_oracle"`.

`potential_history.csv` records the potential for each element with the same
`history_stride` as `charge_history.csv`. The format is
`batch, elem_idx, potential_V`. Enabling it increases computational cost because
`field_solver%refresh` and `compute_mesh_potential` run for each history output.

Requirements for `resume=true`:

| Condition | Details |
|---|---|
| Output | `write_files=true` is required |
| Source | If `restart_from` is unspecified, use `output.dir`; otherwise use `restart_from` |
| Required files | `summary.txt`, `charges.csv`, `rng_state.txt` |
| Optional files | `macro_residuals.csv`; `charge_ledger.csv` is required when schema-v2 ledger metadata is present |
| Behavior | If a required checkpoint is missing, stop instead of falling back to a new run |

`restart_from` changes only the checkpoint read source. New output is always
written to `output.dir`.

During MPI execution:

| File | Contents |
|---|---|
| `rng_state_rankNNNNN.txt` | Random-number state per rank |
| `macro_residuals_rankNNNNN.csv` | Residuals per rank |

`mpi_world_size` in `summary.txt` must match the current number of ranks. Schema v2 also requires matching model, ordered-mesh, and ordered-species fingerprints.

`[[particles.species]].species_key` is a stable identifier for restart fingerprints. When omitted it becomes `species_<1-based index>`; explicit keys must be unique across species.

---

## Relationship to High-Level Notation

The following keys are included in the schema, but they are high-level notation
normalized to runtime keys by the Fortran parser. The parser treats them as the
keys in the right column while loading the file.

| High-level key | Resolution / use |
|---|---|
| `sim.box_origin`, `sim.box_size` | `sim.box_min`, `sim.box_max` |
| `inject_region_mode`, `uv_low`, `uv_high` | `pos_low`, `pos_high` on `inject_face` |
| `mesh.groups`, template `group`, `center_local` | Real coordinates for each template |
| template `placement_mode`, `anchor`, `offset`, `offset_frac` | `center` |
| template `size_mode`, `size_frac` | `size_x`, `size_y`, `size`, `radius`, etc. |

For details, examples, and lint behavior for high-level notation, see
[Configuration](Configuration.en.html).

---

## Validation Rules

| Item | Rule |
|---|---|
| Unknown keys | All are errors |
| `[particles]` | Used only as the container for `[[particles.species]]`. Do not write `key = value` directly under it |
| Old keys | Old names are treated as unknown keys |
| Type | Validated by both the schema and the Fortran parser |
| Value range | `beachx lint` and the runtime parser validate known constraints |

Before running, this is recommended.

```bash
beachx lint beach.toml
```
### `[field]`: element kernel

`element_kernel="point"` is the compatibility default. `element_kernel="triangle_p0"` treats each `q_elem` as total charge distributed with constant density over its triangle. Phase 1 requires `sim.field_solver="direct"`, `sim.field_bc_mode="free"`, `sim.softening=0`, and insulator-only surfaces. Set `[mesh].surface_side` for OBJ input or `surface_side` on every enabled template. `outward_closed` is valid only for consistently oriented, closed two-manifold components.

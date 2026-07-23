title: BEACH Input Parameters Reference

Lang: [English](Parameters.en.md) | [日本語](Parameters.md)

# Input Parameters Reference

This document is the parameter reference for `beach.toml` read by the Fortran runtime.
Unless otherwise noted, units are SI units.

For first-time configuration work, start with
[Design a Simulation Case](ConfigurationRecipes.en.html).

Coordinate and placement helpers are listed as ordinary input keys. See
[Coordinate and placement helper parameters](#coordinate-and-placement-helper-parameters) for the values each helper calculates
or replaces.

| Related document | Contents |
|---|---|
| [Design a Simulation Case](ConfigurationRecipes.en.html) | Task-oriented configuration steps and tuning points |
| [Create and Validate `beach.toml`](Configuration.en.html) | `beachx config`, schema, and lint |
| [Algorithms](Algorithms.en.html) | Entry point to BEM field calculation, particle push, collision, and accumulated-charge procedure |
| [Workflow](Workflow.en.html) | Execution, development, testing, and KUDPC notes |
| [FMM](FMM.en.html) | Selection and accuracy checks for `field_solver="fmm"` |

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

To use an editor schema, put a comment directive with the GitHub Raw URL at the beginning of `beach.toml`.
The Fortran parser does not accept regular keys before the first section, so do not use `"$schema" = "..."`.

```toml
#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.schema.json
```

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

## Official Beginner Case

For the first run, use the [Ten-Minute Tutorial](Tutorial.en.html) and
`examples/tutorial_insulator.toml`. `beachx config init` generates the same
configuration. This beginner case uses a finite image sum for x/y `periodic2` with
`field_periodic_image_layers=1` and `field_periodic_far_correction="none"`.

Add `reservoir_face`, an infinite-periodic correction, and an outer sheath later from
[Design a Simulation Case](ConfigurationRecipes.en.html), after the beginner case
runs successfully.

---

## TOML Hierarchy and Section List

`[sim]`, `[field]`, `[particles]`, `[mesh]`, `[periodic2]`, `[external_boundary]`,
and `[output]` form the normal public configuration. Legacy `[outer_plasma]` and
`[coupling]` remain compatibility inputs, but cannot be mixed with
`[external_boundary]`. Their nesting is:

```text
beach.toml
├── [sim]
├── [field]
├── [particles]
│   └── [[particles.species]]       # one or more array-of-table entries
├── [mesh]
│   ├── [mesh.groups.<name>]        # named child table
│   └── [[mesh.templates]]          # zero or more array-of-table entries
├── [periodic2]
├── [external_boundary]
│   ├── [external_boundary.field]
│   ├── [external_boundary.particles]
│   └── [external_boundary.ordinary_open]
├── [outer_plasma]                  # legacy/runtime representation
├── [coupling]                      # legacy/runtime representation
└── [output]
```

Paths such as `sim.dt` and `external_boundary.field.model` mean “table name.key” in the prose. In TOML, write
`dt = ...` under `[sim]` and `model = ...` under `[external_boundary.field]`.

| TOML table | Parent | Cardinality / requirement | Contents |
|---|---|---|---|
| `[sim]` | root | conditional | Time step, batch count, field solver, boundaries, external fields, sheath correction |
| `[field]` | root | optional | Element-charge discretization kernel |
| `[particles]` | root | required | Container for `[[particles.species]]`; do not put ordinary keys directly under it |
| `[[particles.species]]` | `[particles]` | one or more | Species, injection mode, velocity distribution, macro-particle weight |
| `[mesh]` | root | optional | Selection of OBJ or built-in template input |
| `[mesh.groups.<name>]` | `[mesh]` | zero or more | Placement and scale shared by multiple templates |
| `[[mesh.templates]]` | `[mesh]` | zero or more | Built-in shapes used with `mode="template"` |
| `[periodic2]` | root | conditional | Nonzero mode, zero mode, and lower-boundary model for split periodic2 |
| `[external_boundary.field]` | `[external_boundary]` | required with facade | External-field model and physical/diagnostic parameters |
| `[external_boundary.particles]` | `[external_boundary]` | required with facade | z-high particle handling, inflow, timing, and orbit guards |
| `[external_boundary.ordinary_open]` | `[external_boundary]` | optional | Open faces not owned by the outer model; defaults to `escape` |
| `[outer_plasma]` | root | legacy only | Normalized external-field and return settings |
| `[coupling]` | root | legacy only | Normalized outer refresh and particle-transfer settings |
| `[output]` | root | optional | Output directory, history, checkpoint resume |

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
| `direct` | Exact all-to-all evaluation for small element counts and split references | `free`, or a constrained `periodic2` split reference |
| `treecode` | Approximate evaluation for medium and larger cases | `field_bc_mode="free"` |
| `fmm` | Large-scale evaluation, `periodic2`, FMM core validation | `field_bc_mode="free"` / `"periodic2"` |
| `auto` | Select direct / treecode for point sources or direct / FMM for triangle P0 sources based on element count | `field_bc_mode="free"` |

Use the canonical [solver and field-boundary compatibility table](FieldSolvers.en.html#solver-and-field-boundary-compatibility) to choose the combination.

Common keys:

| Key | Type | Default | Description |
|---|---|---:|---|
| `softening` | float | `1.0e-6` | Softening length for the point kernel [m]. Must be `0` for `triangle_p0` |
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
| `softening` | float | `1.0e-6` | Relaxes point-source singularities; `triangle_p0` requires `0` |
| `field_normalization` | string | `"si"` | Normalize coordinates before direct evaluation |
| `field_length_scale` | float | `1.0` | Used with `field_normalization="length"` or mesh fallback |
| `field_bc_mode` | string | `"free"` | Normally `free`; only the `triangle_p0` split reference may use `periodic2` |

`tree_theta`, `tree_leaf_max`, and `tree_min_nelem` are not used for `direct`.

##### `field_solver = "treecode"`

Builds a source octree. Distant nodes are evaluated with a monopole
approximation, and near nodes are evaluated directly with the selected source kernel. Unlike FMM, it
does not use local expansion and traverses the tree for each evaluation point.

| Key | Type | Default | Description |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | Specify `"treecode"` |
| `softening` | float | `1.0e-6` | Softening for point near sums and monopoles; `triangle_p0` requires `0` |
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
See [FMM](FMM.en.html) for selection and verification, and
[Coulomb FMM core details](FMMCore.en.html) for implementation internals.

| Key | Type | Default | Description |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | Specify `"fmm"` |
| `softening` | float | `1.0e-6` | Used by point-source near sums and FMM evaluation; `triangle_p0` requires `0` |
| `field_normalization` | string | `"si"` | Normalize coordinates before FMM plan construction |
| `field_length_scale` | float | `1.0` | Used with `field_normalization="length"` or mesh fallback |
| `tree_theta` | float | `0.5` | MAC parameter for near/far classification. `0 < theta <= 1` |
| `tree_leaf_max` | int | `16` | Maximum number of sources per source-tree leaf node. `>= 1` |
| `field_bc_mode` | string | `"free"` | `free` / `periodic2` |
| `field_periodic_image_layers` | int | `1` | Number of near image layers for `periodic2` |
| `field_periodic_far_correction` | string | `"none"` | `none` / `auto` / production `cached_kneq0` |
| `field_periodic_ewald_alpha` | float | `0.0` | Ewald parameter for operator generation |
| `field_periodic_ewald_layers` | int | `4` | Real-space / reciprocal cutoff depth for the Ewald generator |
| `field_periodic_cache_dir` | string | `".beach_cache/periodic2"` | Versioned periodic operator cache directory |
| `field_periodic_generation_tolerance` | float | `1e-8` | Generation tolerance included in the cache fingerprint |

`field_periodic_*` only has meaning when `field_bc_mode="periodic2"`.
`tree_min_nelem` is not used when `fmm` is specified explicitly.

##### `field_solver = "auto"`

For `element_kernel="point"`, auto uses direct evaluation below `tree_min_nelem`
and treecode otherwise. For `element_kernel="triangle_p0"`, the same threshold
selects direct or FMM.

| Key | Type | Default | Description |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | Specify `"auto"` |
| `softening` | float | `1.0e-6` | Used by point direct / treecode; `triangle_p0` requires `0` |
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

#### Field Boundary

| Key | Type | Default | Description |
|---|---|---:|---|
| `field_bc_mode` | string | `"free"` | `free` / `periodic2` |
| `field_periodic_image_layers` | int | `1` | Number of near image layers for `periodic2` |
| `field_periodic_far_correction` | string | `"none"` | `none` / `auto` / `cached_kneq0` |
| `field_periodic_ewald_alpha` | float | `0.0` | Ewald parameter for operator generation |
| `field_periodic_ewald_layers` | int | `4` | Outer shell / reciprocal cutoff depth for the Ewald generator |
| `field_periodic_cache_dir` | string | `".beach_cache/periodic2"` | Versioned periodic operator cache directory |
| `field_periodic_generation_tolerance` | float | `1e-8` | Generation tolerance included in the cache fingerprint |

### `[external_boundary]`: Public External-Condition Configuration

The canonical authoring surface for external conditions in new configuration
files is this `[external_boundary]` facade. The later `[outer_plasma]`,
`[coupling]`, and legacy `[sim]` selectors remain only to interpret existing
files and the normalized runtime representation.

Configure external conditions by three responsibilities instead of combining
internal implementation selectors:

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "absorbing_maxwellian"
debye_length = 2.0e-3
thermal_voltage = 3.0

[external_boundary.particles]
mode = "same_batch"
field_evolution_timescale = 1.0
```

- `field` selects the external plasma-response potential/field model. Linear
  and kinetic models cover beyond z-high; unified covers the rough surface to the far region.
- `particles` selects z-high crossing behavior and ownership of the inflow VDF.
- `ordinary_open` applies only to open faces not owned by the outer model.

`field` and `particles` are required when the facade is present. `ordinary_open`
is optional and defaults to `model="escape"`. Omit `[external_boundary]`
entirely for an ordinary case that does not need these controls.

#### `[external_boundary.field]`

| Key | Type | Default | Description |
|---|---|---:|---|
| `model` | string | required | `none` / `linear_debye` / `kinetic_1d` / `unified_linear_response` |
| `kinetic_closure` | string | `absorbing_maxwellian` | For `kinetic_1d` only: `absorbing_maxwellian` / `zhao_charge_driven` |
| `zhao_branch` | string | `auto` | For `zhao_charge_driven` only: `auto` / `a` / `b` / `c` |
| `photoelectron_source_scale` | float | `1` | Analytic source multiplier for `zhao_charge_driven`; use `0` without UV |
| `photoelectron_density_model` | string | `none` | Optional mean density for `kinetic_1d + absorbing_maxwellian`: `none` / `kinetic_mean` |
| `photoelectron_histogram_enabled` | bool | `false` | Set to `true` only to enable the histogram for `linear_debye + same_batch`; otherwise omit the key |
| `infinity_potential` | float | `0` | Reference potential at infinity [V]; accepted only for `linear_debye` |
| `debye_length` | float | required for active models | Length scale for `linear_debye`, `kinetic_1d`, or `unified_linear_response` [m] |
| `thermal_voltage` | float | required for active models | Potential scale for `linear_debye`, `kinetic_1d`, or `unified_linear_response` [V] |
| `unified_grid_points` | int | `129` | Zero-mode Poisson grid points for `unified_linear_response` |
| `accessible_fraction_tolerance` | float | `0.1` | Rough-surface accessible-fraction convergence tolerance for `unified_linear_response` |
| `max_linearity_ratio` | float | `0.25` | Linearity limit for `linear_debye` / `unified_linear_response` |
| `max_gap_ratio` | float | `5` | Interface-to-mesh gap limit for `linear_debye` / `kinetic_1d` |
| `max_local_charge_ratio` | float | `50` | Local plasma-charge estimate limit for `linear_debye` / `kinetic_1d` |
| `photoelectron_histogram_bins` | int | `32` | Number of bins when the `linear_debye + same_batch` histogram is enabled |
| `photoelectron_histogram_energy_max` | float | required with histogram | Upper edge of the enabled histogram [J] |
| `photoelectron_ambient_charge_scale` | float | required with histogram | Ambient signed-charge scale for the enabled histogram [C] |
| `max_photoelectron_charge_ratio` | float | `0.1` | Photoelectron charge-ratio limit for the enabled histogram |

Do not specify `interface_z` or a particle return model. `interface_z` is
derived from `sim.box_max[2]`; the return model is derived from `field.model`
and `particles.mode`. With `model="none"`, no field key other than `model` is allowed.

Start with only the row matching the model and closure. Add other conditional keys only when enabling the corresponding
diagnostic or feature. A key outside the selected row is rejected as a no-op even when its value equals the default.

| Selection | Required or normally configured keys | Conditionally added keys |
|---|---|---|
| `none` | `model` only | none allowed |
| `linear_debye` | `debye_length`, `thermal_voltage` | `infinity_potential`, linearity/gap/local-charge limits, and `same_batch` histogram |
| `kinetic_1d + absorbing_maxwellian` | `debye_length`, `thermal_voltage` | `kinetic_closure`, `photoelectron_density_model`, and gap/local-charge limits |
| `kinetic_1d + zhao_charge_driven` | `debye_length`, `thermal_voltage`, `kinetic_closure`, and required `sim.sheath_*` physics values | source scale, branch, and gap/local-charge limits |
| `unified_linear_response` | `debye_length`, `thermal_voltage` | grid, accessible fraction, and linearity limit |

For `kinetic_closure="zhao_charge_driven"`, `debye_length` and `thermal_voltage` remain schema-required but are not the
physical scales of the Zhao root or profile. Photoemitting Zhao derives those scales from $T_{pe}$ and the reference density;
the no-photo case derives them from ambient $T_e$ and $n_\infty$.

The two configured values are reference inputs for split-interface gap, lateral-field, and local-charge diagnostics.
Normally omit defaults or model-fixed values such as
`zhao_branch="auto"`, zero gauge, and no separate density model.

#### `[external_boundary.particles]`

| Key | Type | Default | Description |
|---|---|---:|---|
| `mode` | string | required | `local_source` / `same_batch` / `zhao_queue` |
| `inflow_model` | string | `auto` | `auto` / `source_vdf` / `infinity_barrier` / `legacy_sheath` |
| `legacy_sheath_model` | string | conditional | With `legacy_sheath`: `floating_no_photo` / `zhao_auto` / `zhao_a` / `zhao_b` / `zhao_c` |
| `steady_start_mode` | string | `none` | Specify only to enable `zhao_floating` with `kinetic_1d + zhao_charge_driven + same_batch` |
| `steady_start_mesh_id` | int | `1` | Mesh ID used with `steady_start_mode="zhao_floating"` |
| `outer_update_stride` | int | `1` | Refresh interval accepted only for `local_source` / `same_batch` with `linear_debye` / `kinetic_1d` [batches] |
| `field_evolution_timescale` | float | `0` | Frozen-field diagnostic timescale for `same_batch` / `zhao_queue` [s] |
| `max_frozen_field_ratio` | float | `0.1` | Frozen-field applicability limit for `same_batch` / `zhao_queue` |
| `outer_orbit_dt` | float | `0` | Outer-orbit step for unified 3-D `same_batch`; specify a positive value |
| `outer_orbit_max_steps` | int | `100000` | Maximum number of outer-orbit steps for unified 3-D `same_batch` |
| `outer_orbit_energy_tolerance` | float | `1e-4` | Relative outer-orbit energy-error limit for unified 3-D `same_batch` |

`mode` selects only whether z-high particles enter an external trajectory and when its result is applied. Select the
reservoir-inflow distribution or correction independently with `inflow_model`.

| `mode` | Meaning | Compatible field |
|---|---|---|
| `local_source` | Retain no external trajectory; handle z-high crossings with `ordinary_open` | all |
| `same_batch` | Classify return or escape for z-high crossings within the same batch | `linear_debye` / `kinetic_1d` / `unified_linear_response` |
| `zhao_queue` | Hold particles in the Zhao reservoir queue and reinject in a later batch | `kinetic_1d` + `zhao_charge_driven` + `zhao_branch="auto"` |

The additional requirements by mode are:

| `particles.mode` and field | Additional requirements |
|---|---|
| `local_source` + `none` / unified | add no transport, time, or orbit keys |
| `local_source` + linear / kinetic | add only `outer_update_stride` when needed |
| `same_batch` + linear / kinetic | `field_evolution_timescale > 0` and `inflow_model="auto"`; optionally add update-stride and time guards |
| `same_batch` + unified | `field_evolution_timescale > 0` and `outer_orbit_dt > 0`; optionally add time and orbit guards |
| `zhao_queue` | `sim.batch_duration > 0`, `field_evolution_timescale > 0`, and a positive photoelectron source; update stride is fixed internally to 1 |

Write `steady_start_*` only for Zhao `same_batch`, `outer_orbit_*` only for unified `same_batch`, and time guards only for
`same_batch` / `zhao_queue`. A key with no effect in the selected mode is an error.

`inflow_model="auto"` delegates inflow to the same 1-D profile for tracked
`linear_debye` / `kinetic_1d`; otherwise it resolves to `source_vdf`.
Those tracked 1-D configurations cannot stack another inflow correction.
`unified_linear_response` owns the outer orbit but not the inflow VDF, so it can
use `source_vdf`, `infinity_barrier`, or `legacy_sheath`.

#### `[external_boundary.ordinary_open]`

| Key | Type | Default | Description |
|---|---|---:|---|
| `model` | string | `escape` | `escape` / `potential_barrier` |

With `local_source`, `ordinary_open` also handles z-high. With `same_batch` or
`zhao_queue`, the outer model owns z-high and `ordinary_open` handles only the
remaining open faces. The inflow choice
`particles.inflow_model="infinity_barrier"` and the outflow choice
`ordinary_open.model="potential_barrier"` may be selected independently.

The facade is normalized at load time to the legacy `[sim]` / `[outer_plasma]` /
`[coupling]` runtime representation. Mixing both syntaxes in one file is an error.
See [Configure the External Boundary](OuterPlasmaModels.en.html)
for the selection workflow and
[Boundary Configuration Migration](BoundaryConfigurationMigration.en.html) for
the legacy mapping.

### `[periodic2]`: Nonzero Mode, Zero Mode, and Lower Boundary

`[periodic2]` is a top-level table, not a child of `[sim]`.

| Key | Default | Meaning |
|---|---:|---|
| `nonzero_mode_backend` | required | `panel_spectral_reference` / `cached_kneq0` |
| `zero_mode_policy` | required | `exclude_k0` |
| `lower_boundary_model` | required | `symmetric_vacuum` / `e_bottom_zero` |
| `reference_mode_layers` | `4` | Fourier-mode cutoff |
| `panel_quadrature_order` | `12` | panel-area quadrature order |
| `interface_sample_n` | `5` | diagnostic samples along each interface axis |
| `interface_phi_tolerance` | `1e-3` | upper bound on the nonzero-mode potential ratio |
| `interface_field_tolerance` | `1e-3` | upper bound on the nonzero-mode field ratio |

Legacy `periodic2` uses `field_solver="fmm"`. The small-system split reference instead explicitly selects
`field_solver="direct"`, `nonzero_mode_backend="panel_spectral_reference"`, `zero_mode_policy="exclude_k0"`, and a lower
boundary model. `symmetric_vacuum` is the parameter-free homogeneous-vacuum closure; `e_bottom_zero` remains available for
old-run reproduction.

### `[outer_plasma]`: Compatibility and Normalized-Runtime Reference

> **This is not an authoring API for new configurations.** Use `[external_boundary.field]` in new files.
> The raw names below are retained for reading legacy files and the normalized runtime contract;
> they cannot be added alongside `[external_boundary]`.

The facade derives `interface_z` and `return_model`. Legacy `[outer_plasma]` is a top-level table, not a child of
`[periodic2]`.

| Key | Default | Meaning |
|---|---:|---|
| `model` | required | `linear_debye` / `kinetic_1d` / `unified_linear_response` |
| `interface_z` | required | upper-z interface; initially the top of the box |
| `infinity_potential` | `0` | reference potential at infinity [V] |
| `debye_length` | required | length scale for linear/`absorbing_maxwellian` tails and split diagnostics; not the Zhao profile's physical scale |
| `thermal_voltage` | required | potential scale for linearity and split diagnostics; not the Zhao profile's physical scale |
| `unified_grid_points` | `129` | unified zero-mode Poisson grid points; at least 17 |
| `accessible_fraction_tolerance` | `0.1` | maximum accessible-fraction change after doubling rough-surface samples along both axes |
| `max_linearity_ratio` | `0.25` | upper bound on `abs(phi_t-phi_inf)/thermal_voltage` |
| `max_gap_ratio` | `5` | upper bound on `(z_t-z_mesh,max)/lambda` |
| `max_local_charge_ratio` | `50` | upper bound on the local mean-plasma charge estimate ratio |
| `kinetic_closure` | `absorbing_maxwellian` | density/VDF closure for `kinetic_1d`: `absorbing_maxwellian` / `zhao_charge_driven` |
| `zhao_branch` | `auto` | Zhao closure branch: `auto` / `a` / `b` / `c` |
| `photoelectron_source_scale` | `1` | independent multiplier for the analytic Zhao photoelectron source; use `0` without UV; distinct from queue occupancy $\eta$ |
| `photoelectron_density_model` | `none` | `none` / `kinetic_mean`; the latter adds mean photoelectron density to `kinetic_1d` |
| `photoelectron_histogram_enabled` | `false` | enable the outward z-high photoelectron histogram and applicability check |
| `photoelectron_histogram_bins` | `32` | normal-kinetic-energy histogram bins |
| `photoelectron_histogram_energy_max` | required with histogram | histogram upper edge [J]; positive |
| `photoelectron_ambient_charge_scale` | required with histogram | ambient signed-charge scale for linear-model applicability [C] |
| `max_photoelectron_charge_ratio` | `0.1` | upper bound on `abs(Q_pe,batch)/Q_ambient,scale` |
| `return_model` | `none` | ID of the analytic 1-D return or unified explicit 3-D orbit |

Applicability failures do not fall back to a legacy model. See `examples/periodic2_linear_outer_reference.toml`.
The MPI-global photoelectron histogram stores signed charge, kinetic energy, tangential momentum, and count in
`photoelectron_histogram.csv`; it does not control particle return. See `examples/periodic2_photoelectron_return.toml`.

#### Normalized `kinetic_1d` Contract (Standard and Recommended)

Selecting `external_boundary.field.model="kinetic_1d"` through the public facade normalizes the runtime
`outer_plasma.model` to `kinetic_1d`. Production runs use it with `cached_kneq0`. Negative and positive z-high
`reservoir_face` species define the infinity electron and ion VDFs.

| Item | Contract |
| --- | --- |
| Gauge | `phi(infinity)=0`; reject nonzero `infinity_potential` |
| Far boundary | `absorbing_maxwellian` uses a Robin tail of length `debye_length`; photoemitting Zhao instant return derives $\lambda_{D,pe}$ from $T_{pe},n_{ref}$; no-photo Zhao derives $\lambda_{D,e}$ from ambient $T_e,n_\infty$; queue mode ends at $L=10\lambda_{D,pe}$ |
| Closure | default `absorbing_maxwellian`, or `zhao_charge_driven`, which retains the accumulated-charge interface-field condition |
| Supported branch | monotone for `absorbing_maxwellian`; Zhao A/B/C, including nonmonotone Type A, for `zhao_charge_driven` |
| Unsupported | virtual cathodes, trapped populations, and sub-Bohm inflow under `absorbing_maxwellian` |
| Nonlinear solve | analytic bordered-tridiagonal Jacobian plus branch-preserving Newton |
| Recovery | pseudo-transient and interface-field continuation; accept only after the original Poisson residual passes |
| Fallback | never return another sheath model or a held previous profile as a converged solution |

The public `external_boundary.particles.mode="same_batch"` normalizes to
`return_model="kinetic_1d_profile_return"`,
`particle_transfer_mode="electrostatic_1d_instant_return"`, and `outer_queue_enabled=false`.
Do not author these raw IDs; after resolution they denote the following instant map.

1. Map the infinity VDF to the interface using the refreshed `phi_interface-phi_infinity`.
2. Use the same discrete profile and Robin tail to classify escape or a turning point.
3. Construct the state corresponding to the analytically integrated round trip
   and return the particle at the same simulation time and in the same batch.

Scope of the instant map:

- The model targets stationary and quasistationary sheaths and supports mean current and detachment force after equilibration.
- It does not represent delayed return current during UV turn-on or other transients.
- `tau_outer/field_evolution_timescale` bounds the quasistationary approximation.
- With `tau_outer/batch_duration >= 1`, do not treat batch history as a physical return-current time history.

Set `photoelectron_source_scale=0` when UV is absent. This path requires no enabled `photo_raycast` species,
`sheath_photoelectron_ref_density_cm3`, or $T_{pe}$. It uses the quasineutral-region $n_\infty,T_e,u_e,u_i$ from the
z-high ambient electron and ion species to derive the incoming-electron reservoir normalization, cutoff, and velocity map
from the same Zhao VDF.

$E_I=0$ is the flat Type-B/C junction and $E_I<0$ selects Type C. Zero current remains a diagnostic,
not a per-batch root condition, so a stationary charge state can recover the legacy no-photo Type-C floating root.

The public `external_boundary.particles.mode="zhao_queue"` normalizes to runtime `outer_queue_enabled=true`, so
outer-flight delay enters the batch history for cases such as strong-UV turn-on. This mode requires
`external_boundary.field.kinetic_closure="zhao_charge_driven"`, a positive `batch_duration` resolved either directly or from
`dt * batch_duration_step`.

Runtime fixes the update stride to 1 and disables the histogram, so write neither key in the public input.

1. Pop due events from each rank-local queue at the batch start.
2. Divide the remaining global photoelectron inventory by horizontal area, then refresh the profile and Zhao population scale
   $\eta$ that match the finite column over $0\le z\le10\lambda_{D,pe}$.
3. Advance fresh sources and due returns, resolve each outward particle as return before $L$ or reservoir escape at $L$, then
   enqueue it at $t_{due}=t_{mid}+\tau_{outer}$ with the batch midpoint as its crossing time.

Despite its output name `outer_photoelectron_population_fraction`, $\eta$ is an occupancy scale relative to the stationary
reference population, not a probability. The solver follows the path connected to $\eta=0$ over $0\le\eta\le16$, permits
$\eta>1$.

It does not clamp, fall back to a target-independent full-population solution, or jump to a disconnected branch. Queue mode
requires `zhao_branch="auto"`; only continuous branch transitions satisfying the degeneracy condition and a column that
increases monotonically with $\eta$ are accepted. Paths containing a fold and unreachable targets stop the run.

Events are released only at batch starts, and a terminal state is not reintegrated after the outer field changes. Queue mode
does not use a Robin tail outside $L$; reaching $L$ is absorption/escape into the reservoir.

For each event, `tau_outer` plus the quantization delay to the next batch-start poll and the half-batch midpoint uncertainty must not exceed
`max_frozen_field_ratio * field_evolution_timescale`. Configuration validation applies the same limit to `batch_duration`.

Check convergence in `batch_duration`, tracked-particle count, horizontal area, effective-interface location, and profile grid.

Combination constraints:

- `reservoir_potential_model`, Zhao injection correction, and nonzero `b0` are rejected.
- `kinetic_closure="zhao_charge_driven"` requires `model="kinetic_1d"`, `infinity_potential=0`, and
  `photoelectron_density_model="none"`. It cannot be combined with legacy `sheath_injection_model` or
  `reservoir_potential_model` corrections.
- Explicit `zhao_branch="a"`, `"b"`, or `"c"` requires `zhao_charge_driven`; `auto` searches the available branches.
- `zhao_charge_driven` requires quasineutral ambient electron/ion densities and positive inward electron and ion drifts,
  and rejects `sheath_reference_coordinate`.
- With `photoelectron_source_scale>0`, exactly one negative `photo_raycast` species and a positive
  `sheath_photoelectron_ref_density_cm3` are required. With `photoelectron_source_scale=0`, an enabled `photo_raycast`
  species is rejected and queue mode is unavailable.
- Only `sheath_electron_drift_mode="normal"` and `sheath_ion_drift_mode="normal"` are accepted.
- The negative `photo_raycast` species used by Zhao requires `normal_drift_speed=0`, and ion temperature must satisfy the
  cold-ion condition $T_i\le0.1T_e$.
- Photoemitting Zhao uses $T_{pe}$ and the $\lambda_{D,pe}$ derived from photoelectron temperature and reference density.
  No-photo Zhao switches the same normalization to ambient $T_e,n_\infty$ and reports the derived $\lambda_{D,e}$.
- `debye_length` and `thermal_voltage` do not change the Zhao root or profile. They are still required as reference inputs for
  the split-interface `interface_eta_gap`, lateral potential/field, and local-charge diagnostics.
- Tracked `photo_raycast.emit_current_density_a_m2` must agree within 1% with
  `photoelectron_source_scale * |q_{pe}|n_{ref}\sin(\alpha)v_{th,pe}/(2\sqrt{\pi})`, where
  $v_{th,pe}=\sqrt{2T_{pe}/m_{pe}}$ and $T_{pe}$ is in joules.
- The analytic raw current enters the tracked-source consistency check and current-density diagnostics, but not the root,
  surface charge, or ledger. Only tracked emission and reabsorption update the latter two, and $\eta$ does not scale the raw
  photoelectron emission-current term in the current diagnostic.
- This first version is an effective-plane approximation. It does not self-consistently connect `ray_direction` or the VDF
  arriving from a rough surface to the Zhao outer population. `ray_direction` and `sheath_alpha_deg` independently specify
  illumination-ray sampling of emitting surfaces and the analytic source, respectively.
- Test Zhao convergence with the profile grid, effective interface location, tracked-particle count, and `dt` or batch
  resolution. Do not interpret changes in `debye_length` or `thermal_voltage` as a profile-convergence test.
- With the default `absorbing_maxwellian`, `photoelectron_density_model="kinetic_mean"` uses the first negative
  `photo_raycast` species to supply a half-Maxwellian flux.
- The mean closure supplies only the outer profile. Tracked particles update surface charge, so statistical return current is not added again.
- Every tracked-return species requires `deposit_opposite_charge_on_emit=true`.

See [Particle Escape and Return](ParticleEscapeReturn.en.html).

The default-closure example is
[`periodic2_kinetic_outer.toml`](../examples/periodic2_kinetic_outer.toml). The photoemitting charge-driven example is
[`periodic2_zhao_charge_driven_outer.toml`](../examples/periodic2_zhao_charge_driven_outer.toml), and the no-photo example is
[`periodic2_zhao_no_photo_outer.toml`](../examples/periodic2_zhao_no_photo_outer.toml).

The transient-queue example is
[`periodic2_zhao_transient_outer.toml`](../examples/periodic2_zhao_transient_outer.toml).

The transient-queue example is an expected-fail guard fixture that rejects a long flight at its stated physical timescale,
not a successful physical-validation example.

The model assumptions are
documented in [ADR 0001](adr/0001-kinetic-outer-plasma.md) and
[ADR 0003](adr/0003-zhao-charge-driven-outer-closure.md).

#### Normalized `unified_linear_response` Contract (Advanced and Specialized)

`unified_linear_response` is not a higher-accuracy replacement for `kinetic_1d`. Select it as rough-surface linear screening only
when roughness and plasma response occupy the same region, no split window is available, and the linearity gate passes.

| Item | Contract |
| --- | --- |
| Zero mode | one Poisson grid from surface projection to the far boundary |
| Rough surface | include plasma-accessible area and linear mean-plasma charge |
| Nonzero mode | join a $\sqrt{k^2+\kappa^2}$ tail above the highest surface point |
| `interface_z` | particle ownership plane, not a field truncation plane |
| Geometry/kernel | require single-valued topography and `triangle_p0` |
| Photoelectron mean | require `photoelectron_density_model="none"` |
| Particle transfer | `none` or `electrostatic_3d_explicit_orbit` |
| 3D orbit | require `b0=0`, fixed step, and energy/frozen-field error contracts |
| Applicability | fail closed beyond the configured linearity bound |

Numerical and applicability rules:

- `unified_grid_points >= 17` is required; the default is `129`.
- `accessible_fraction_tolerance` bounds the maximum accessible-fraction change after doubling height samples along both periodic axes.
- The solve uses the refined samples and stops during initialization when the tolerance is exceeded.
- Production studies should demonstrate grid refinement of the reported observables.

See `examples/periodic2_unified_linear_response.toml`, `examples/periodic2_unified_explicit_orbit.toml`, and
`docs/adr/0002-unified-periodic-outer-domain.md`.

### `[coupling]`: Compatibility and Normalized-Runtime Reference

> **This is not an authoring API for new configurations.** Use `[external_boundary.particles]` in new files.
> The raw names below are retained for reading legacy files and the normalized runtime contract.

The facade derives `update_mode`, `particle_transfer_mode`, and `outer_queue_enabled` from its model and mode.

| Key | Type | Default | Description |
| --- | --- | ---: | --- |
| `update_mode` | string | `"explicit"` | Only `explicit` is supported; refresh the outer profile at explicit update points |
| `particle_transfer_mode` | string | `"none"` | Facade-derived `none` / `electrostatic_1d_instant_return` / `electrostatic_3d_explicit_orbit` |
| `steady_start_mode` | string | `"none"` | `none` / `zhao_floating`; initialize a fresh run from a Zhao zero-current stationary root and its matching plane charge |
| `steady_start_mesh_id` | int | `1` | `mesh_id` of the horizontal plane that receives the `zhao_floating` charge in proportion to triangle area |
| `outer_update_stride` | int | `1` | Batch interval between outer-profile refreshes |
| `field_evolution_timescale` | float | `0` | Frozen-field comparison time [s]; positive for 1-D return |
| `max_frozen_field_ratio` | float | `0.1` | Limit on instant `tau_outer`, or queue `tau_outer` plus next-poll delay and half-batch midpoint uncertainty, divided by `field_evolution_timescale`; queue mode also applies it to `batch_duration` |
| `outer_orbit_dt` | float | `0` | Fixed 3-D outer-orbit step [s]; positive in 3-D mode |
| `outer_orbit_max_steps` | int | `100000` | 3-D outer-orbit step limit; reaching it stops instead of discarding |
| `outer_orbit_energy_tolerance` | float | `1e-4` | Relative total-energy error limit for a 3-D outer orbit |
| `outer_queue_enabled` | bool | `false` | In the supported Zhao configuration, retain outer flight across batches and close the transient queued photoelectron column |

Runtime return/transfer pairs resolved by the facade:

| `outer_plasma.model` | `outer_plasma.return_model` | `coupling.particle_transfer_mode` |
| --- | --- | --- |
| `linear_debye` | `electrostatic_1d_instant_return` | `electrostatic_1d_instant_return` |
| `kinetic_1d` | `kinetic_1d_profile_return` | `electrostatic_1d_instant_return` |
| `unified_linear_response` | `electrostatic_3d_explicit_orbit` | `electrostatic_3d_explicit_orbit` |

The return and transfer strings are intentionally different for `kinetic_1d_profile_return`. With active 1-D transfer,
`linear_debye` and `kinetic_1d` also own inflow through the same profile, so normalized runtime state has both
`reservoir_potential_model` and `sheath_injection_model` set to `none`. `infinity_potential` is fixed to zero for `kinetic_1d` and
`unified_linear_response`.

Stationary warm start:

In the public facade, select `external_boundary.particles.mode="same_batch"`,
`steady_start_mode="zhao_floating"`, and `steady_start_mesh_id` when needed.
The facade derives the raw transfer ID and `outer_queue_enabled`.

Before the first batch of a fresh run, `zhao_floating` solves the Zhao zero-current stationary root from the configured
infinity reservoir and UV source. It constructs the kinetic profile with `phi(infinity)=0`. For horizontal area $A$ and
stationary-root interface field $E_I$, the total charge placed on the selected plane is

$$
Q_{seed}=
\begin{cases}
2\epsilon_0 A E_I, & \texttt{symmetric_vacuum},\\
\epsilon_0 A E_I, & \texttt{e_bottom_zero}.
\end{cases}
$$

The charge is distributed over triangles of `steady_start_mesh_id` in proportion to area; all other meshes start with zero
charge. Thus, in a plane-plus-sphere input where mesh 1 is the plane, `steady_start_mesh_id=1` seeds only the plane and leaves
the sphere neutral.

The first outer refresh, corrected inflow from the infinity reservoir, and instant return or escape at the interface all use
this same profile. Analytic current is not deposited separately; subsequent surface-charge updates still come only from
tracked particles.

This mode explicitly selects an initial condition on a stationary branch instead of time-integrating the physical relaxation
from an uncharged state. Use it for stationary and quasistationary observables, not for claims about UV turn-on or delayed
return current. The queue transient closure remains a separate mode for that purpose.

Public-facade `zhao_floating` requires:

- `external_boundary.field.model="kinetic_1d"` +
  `kinetic_closure="zhao_charge_driven"` + `external_boundary.particles.mode="same_batch"`;
- `zero_mode_policy="exclude_k0"` and a supported lower-boundary model;
- zero initial charge on every mesh for a fresh run; with `output.resume=true`, restore checkpoint mesh charge and outer
  state without reseeding;
- `mesh.mode="template"`, with the selected mesh forming one coplanar horizontal plane that covers periodic-cell area $A$
  and lies below the outer interface.

A warm start alone does not establish physical uniqueness or stability. Publication use should confirm that an independently
relaxed or perturbed seed returns to the same stationary observables.

Transfer rules:

- Normalized `outer_plasma.return_model` and `coupling.particle_transfer_mode` form the matching pair shown above.
- The 1-D path supports only the open z-high interface, x/y periodic wrapping, and `b0=0`.
- `kinetic_1d` requires exactly one enabled negative and one enabled positive z-high `reservoir_face` species.
- Instant mode requires a positive `field_evolution_timescale` and uses `max_frozen_field_ratio` as an applicability limit.
- Public-facade `zhao_queue` requires `kinetic_1d` + `zhao_charge_driven` + `zhao_branch="auto"`, a positive
  `batch_duration` resolved directly or from `dt * batch_duration_step`; runtime fixes `outer_update_stride` to 1.
  Each event's `tau_outer`, delay to the next batch-start poll, and half-batch midpoint uncertainty are bounded by
  `max_frozen_field_ratio * field_evolution_timescale`; the same bound applies to `batch_duration`.
- Queue mode exposes no photoelectron-histogram setting. Persistent queuing remains unavailable for a 3-D explicit orbit.
- See `examples/periodic2_outer_particle_transfer.toml` and `examples/periodic2_zhao_transient_outer.toml`.

Photoelectron-histogram rules:

- Enabling the histogram requires `photoelectron_histogram_energy_max` [J] and `photoelectron_ambient_charge_scale` [C].
- Every `photo_raycast` species using tracked outer transfer requires `deposit_opposite_charge_on_emit=true`.
- The histogram is diagnostic only. `field.model` and `particles.mode` resolve return and escape.
- The public facade supports it only for `linear_debye + same_batch`.
- `photoelectron_density_model="kinetic_mean"` cannot be enabled simultaneously because it requires another return-model branch.

### Combined periodic2 and External-Boundary Constraints

Periodic2 requires `sim.use_box=true` and exactly two periodic axes. A configuration with outer transfer uses x/y as those axes
and an open z-high interface. The same periodicity applies to field evaluation, collision, and `photo_raycast`.

| Far correction | Meaning |
| --- | --- |
| `none` | finite-image approximation |
| `auto` | compatibility input; currently `none` |
| `cached_kneq0` | production nonzero mode using a reusable versioned operator |

The removed `m2l_root_oracle` value is rejected at startup with guidance to use `cached_kneq0`.

With `cached_kneq0`, the `exclude_k0` provider adds the physical zero mode exactly once.
`symmetric_vacuum` uses $\pm Q/(2\epsilon_0A)$ above and below; `e_bottom_zero` uses zero below and $Q/(\epsilon_0A)$ above.

### `[sim]`: External Fields and Physical Values Used by Boundaries

`phi_infty` and the `sheath_*` values are physical inputs used with the new
`[external_boundary]` facade. In contrast, `reservoir_potential_model`,
`open_boundary_model`, and `sheath_injection_model` are legacy selectors. Do not
write those three with the facade; let `external_boundary.particles` and
`external_boundary.ordinary_open` derive them.

Their raw descriptions below remain for reading legacy inputs and normalized
runtime state; they are not an alternative authoring API for new files.

| Key | Type | Default | Description |
|---|---|---:|---|
| `e0` | float[3] | `[0,0,0]` | Uniform external electric field [V/m] |
| `e0_abs` | float | unspecified | Magnitude of the uniform external electric field [V/m] |
| `e0_phi_xy_deg` | float | `0.0` | Azimuthal angle in the xy plane when `e0_abs` is specified [deg] |
| `e0_phi_z_deg` | float | `0.0` | Elevation angle from the xy plane when `e0_abs` is specified [deg] |
| `b0` | float[3] | `[0,0,0]` | Uniform magnetic field [T] |

The uniform external electric field can be specified directly as
`e0 = [Ex, Ey, Ez]`, or by `e0_abs` and angles. Mixing the two forms is an error.

#### Physical Values and Legacy-Selector Compatibility Reference

The three rows marked as legacy selectors are compatibility inputs. All other
rows are physical values or common boundary-processing values that the facade
may still use.

| Key | Type | Default | Description |
|---|---|---:|---|
| `reservoir_potential_model` | string | `"none"` | Legacy selector: `none` / `infinity_barrier` |
| `phi_infty` | float | `0.0` | Reference potential at infinity [V] |
| `open_boundary_model` | string | `"escape"` | Legacy selector: `escape` / `potential_barrier` |
| `multiple_box_events_policy` | string | `"abort"` | `abort` / `soft_discard`; action when a step exceeds the box-event limit |
| `multiple_box_events_soft_discard_count_limit` | int | `1000` | abort after the cumulative soft-discard count exceeds this value |
| `multiple_box_events_soft_discard_abs_charge_limit` | float | `1e-12` | abort after cumulative absolute soft-discarded macro charge [C] exceeds this value |
| `injection_face_phi_grid_n` | int | `3` | `N x N` evaluation grid for injection-face average potential |
| `raycast_max_bounce` | int | `16` | Maximum number of reflections for `photo_raycast` |
| `sheath_injection_model` | string | `"none"` | Legacy selector: `none` / `zhao_auto` / `zhao_a` / `zhao_b` / `zhao_c` / `floating_no_photo` |
| `sheath_alpha_deg` | float | `60.0` | Solar elevation angle for the Zhao sheath [deg] |
| `sheath_photoelectron_ref_density_cm3` | float | `64.0` | Reference photoelectron density for the Zhao sheath [cm^-3] |
| `sheath_reference_coordinate` | float | unspecified | Reference plane position for the sheath 1D coordinate [m] |
| `sheath_electron_drift_mode` | string | `"normal"` | `normal` / `full` |
| `sheath_ion_drift_mode` | string | `"normal"` | `normal` / `full` |

Compatibility-runtime combinations and details:

- In the public facade, `particles.inflow_model="legacy_sheath"` and
  `"infinity_barrier"` are mutually exclusive. In normalized runtime state,
  this appears as `sheath_injection_model != "none"` with
  `reservoir_potential_model="none"`.
- See [Sheath Injection Closures](SheathInjectionClosures.en.html) for the
  legacy-sheath values and behavior.
- See [`reservoir_face` Inflow and Velocity Sampling](ReservoirInjection.en.html) for flux and velocity.
- See [Sheath Injection Closures](SheathInjectionClosures.en.html) and
  [Particle Escape and Return](ParticleEscapeReturn.en.html) for sheath correction, reflection, and return.

Evaluation with
`external_boundary.particles.inflow_model="infinity_barrier"` (represented as
`reservoir_potential_model="infinity_barrier"` at normalized runtime):

- Evaluate injection-face average potential from the field and potential refreshed at the beginning of each batch.
- Include the point / `triangle_p0` kernel, periodic2 terms, zero mode, outer profile, and `e0` under the particle-motion convention.
- The same `N x N` grid gives the population standard deviation, minimum, and maximum without extra potential evaluations.
- For a Maxwellian reservoir, the MPI root warns on the first and final batch when
  `abs(q_particle) * phi_std > 0.1 * (k_B*T + 0.5*m_particle*u_normal^2)`.

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

Decision with `external_boundary.ordinary_open.model="potential_barrier"`
(represented as `open_boundary_model="potential_barrier"` at normalized
runtime):

1. Evaluate BEM potential `phi_boundary` at the crossing point with the particle-motion snapshot convention, including local `e0` potential.
2. Compute the barrier `q_particle * (phi_infty - phi_boundary)`.
3. If it is positive and exceeds `0.5 * m_particle * v_normal^2`, reverse the normal velocity. Otherwise, the particle escapes.

A uniform field has no finite potential at infinity. When it is combined with this model, supply `phi_infty` as a consistent
effective reservoir reference.

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
`external_boundary.particles.inflow_model!="legacy_sheath"`. This corresponds
to `sim.sheath_injection_model="none"` in normalized runtime state.

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

Each emitted photoelectron uses `w_hit` as its weight and is tracked as an
ordinary particle. Surface return is absorbed as an ordinary collision.
`external_boundary.ordinary_open.model` controls ordinary open faces, while
transport resolved from `external_boundary.particles.mode` controls return or
escape at the z-high interface.

---

### Legacy Runtime Selector `sim.sheath_injection_model`

For new files, use
`external_boundary.particles.inflow_model="legacy_sheath"` and
`legacy_sheath_model`. The details below remain for interpreting legacy inputs
and normalized runtime values. `sim.sheath_injection_model` groups existing
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

#### `[[mesh.templates]]`: Built-In Shapes

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

##### `kind = "plane"`

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

##### `kind = "plate_hole"` / `"plane_hole"`

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

##### `kind = "disk"`

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

##### `kind = "annulus"`

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

##### `kind = "box"`

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

##### `kind = "cylinder"`

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

##### `kind = "sphere"`

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

`conductor` constraints:

- Redistribute charge with direct Coulomb coefficients under `field_bc_mode="free"`.
- Do not combine it with `field_bc_mode="periodic2"`.
- Large conductor element counts increase per-batch dense-solve cost.

---

### `[field]`: Element Kernel

`element_kernel="point"` is the compatibility default, and `sim.softening` applies to this point kernel.

Rules for `element_kernel="triangle_p0"`:

| Item | Rule |
| --- | --- |
| Source | Treat each `q_elem` as a constant surface-charge density over its triangle |
| Solver | `direct` / `treecode` / `fmm` / `auto`; `auto` selects direct / FMM with `tree_min_nelem` |
| Requirements | `sim.softening=0` and insulator-only surfaces |
| Treecode | exact panel near + monopole far |
| FMM | exact panel near + exact triangle P2M |
| Surface side | `[mesh].surface_side` for OBJ; `surface_side` on every enabled template |

Use `outward_closed` only for consistently oriented, closed two-manifold components.

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
| `outer_plasma_profile.csv` | Profile for a ready `kinetic_1d` / `unified_linear_response` outer state; a conditional checkpoint |
| `photoelectron_histogram.csv` | Previous-batch and cumulative histogram when `photoelectron_histogram_enabled=true`; a conditional checkpoint |
| `outer_event_queue.csv` / `outer_event_queue_rankNNNNN.csv` | Active events when `outer_queue_enabled=true`; the former is serial and the latter is one conditional checkpoint per MPI rank |
| `mesh_potential.csv` | When `write_mesh_potential=true` |
| `charge_history.csv` | When `history_stride > 0` |
| `potential_history.csv` | When `write_potential_history=true` and `history_stride > 0` |
| `performance_profile.csv` | When `BEACH_PROFILE=1` |
| `rng_state.txt` / `rng_state_rankNNNNN.txt` | Serial or MPI rank-local random-number state |
| `macro_residuals.csv` | One MPI-global macro-particle residual file |
| `charge_ledger.csv` | Per-species signed-charge flux, counts, and restartable cumulative values |

When the histogram state is ready, `summary.txt` adds:

| Group | Keys |
| --- | --- |
| Histogram definition | `photoelectron_histogram_bins`, `photoelectron_histogram_energy_max_J` |
| Progress | `photoelectron_last_completed_batch` |
| Cumulative | `photoelectron_cumulative_signed_charge_C`, `photoelectron_cumulative_kinetic_energy_J`, `photoelectron_cumulative_count` |
| Previous batch | `photoelectron_previous_signed_current_A`, `photoelectron_previous_charge_ratio` |
| Applicability | `photoelectron_max_charge_ratio`, `photoelectron_linear_applicability_status` |

`coupling_steady_start_mode`, `coupling_steady_start_mesh_id`, and `coupling_outer_queue_enabled` are always written to
`summary.txt`. At a fresh `zhao_floating` startup, one standard-output record beginning with
`zhao_steady_start_branch=...` reports the resolved branch, $E_I$, $Q_{seed}$, and mesh ID. Resume does not reseed and instead
reports the restored batch count as `zhao_steady_start_restored_after_batches=...`.

Only when `coupling_outer_queue_enabled=T` is the following queue state added:

| Group | Keys |
| --- | --- |
| Closure | `outer_photoelectron_population_fraction` |
| Column | `outer_photoelectron_column_per_area_m2`, `outer_photoelectron_column_target_per_area_m2`, `outer_photoelectron_column_residual_per_area_m2` |
| Queue stock | `outer_queue_event_count`, `outer_queue_signed_charge_C`, `outer_queue_fingerprint` |

See [Configuration-specific output](OutputGuide.en.html#locate-configuration-specific-values) to locate these values.

[Files Used for Resume](OutputGuide.en.html#files-used-for-resume) consolidates column definitions and conditional checkpoint requirements.

Evaluation rules for `mesh_potential.csv`:

- Record potential [V] at each element centroid.
- Use `1/softening` for the self term when `softening > 0`; otherwise use an area-equivalent radius approximation.
- Add the explicit image shell for `periodic2`.
- Add the cached nonzero mode and boundary-specific `k=0` for `cached_kneq0`.

`potential_history.csv`:

- Use the same `history_stride` as `charge_history.csv`.
- Write `batch, elem_idx, potential_V`.
- Run `field_solver%refresh` and `compute_mesh_potential` for each history output, which increases computational cost.

Requirements for `resume=true`:

| Condition | Details |
|---|---|
| Output | `write_files=true` is required |
| Source | If `restart_from` is unspecified, use `output.dir`; otherwise use `restart_from` |
| Required files | `summary.txt`, `charges.csv`, and either serial `rng_state.txt` or every MPI `rng_state_rankNNNNN.txt` |
| Conditional files | `charge_ledger.csv` with ledger metadata, `outer_plasma_profile.csv` for a ready outer state, `photoelectron_histogram.csv` when the histogram is enabled, and serial `outer_event_queue.csv` or every MPI `outer_event_queue_rankNNNNN.csv` when the queue is enabled |
| Optional state | Restore the global residual when `macro_residuals.csv` exists |
| Behavior | If a required checkpoint is missing, stop instead of falling back to a new run |

`restart_from` changes only the checkpoint read source. New output is always written to `output.dir`.

During MPI execution:

| File | Contents |
|---|---|
| `rng_state_rankNNNNN.txt` | Random-number state per rank |
| `outer_event_queue_rankNNNNN.csv` | Rank-local active events for the transient Zhao queue; every rank writes one |
| `macro_residuals.csv` | One global residual shared by all ranks and written by the root |

Resume consistency rules:

- Reject checkpoints with legacy `macro_residuals_rankNNNNN.csv` instead of converting them implicitly.
- Match `mpi_world_size` in `summary.txt` to the current rank count.
- Schema v2/v3/v4 requires matching model, ordered-mesh, and ordered-species fingerprints.
- Schema v3 outer profiles require `field_V_m` and `charge_density_C_m3`.
- A schema-v4 queue resume stops unless the transient Zhao state, queue-file schema 2, rank, world size, completed batch,
  global count, signed charge, and all-rank queue fingerprint match.
- `[[particles.species]].species_key` is stable. Omission yields `species_<1-based index>`; explicit values must be unique.

---

## Coordinate and Placement Helper Parameters

These are ordinary TOML parameters, but the loader uses them to calculate the physical coordinates or dimensions in the third
column. Distinguish combinations that fail validation from those that intentionally replace an explicit value.

| Key and value | Type / condition | Calculated value | Relationship to a directly specified key |
|---|---|---|---|
| `sim.box_origin`, `sim.box_size` | float[3]. Specify both, with `box_size > 0` | `box_min = box_origin`; `box_max = box_origin + box_size` | Combining either with `box_min` or `box_max` is an error |
| `inject_region_mode="face_fraction"`, `uv_low`, `uv_high` | `uv_*` are float[2] in `[0,1]`; only for `reservoir_face` / `photo_raycast` | `pos_low`, `pos_high` on `inject_face` | Combining with either `pos_low` or `pos_high` is an error |
| Template `placement_mode="box_anchor"`, `anchor`, and either `offset` or `offset_frac` | Anchor is the box center or a face center; `offset` is in meters and `offset_frac` is relative to box size | Template `center` | Combining with `center` is an error; combining `offset` and `offset_frac` is also an error |
| Template `size_mode="box_fraction"`, `size_frac` | Supported for `plane` / `plane_hole` / `plate_hole` / `box` / `sphere` / `cylinder` | Shape-specific dimensions listed below | Calculated dimensions replace corresponding explicitly specified dimension keys without an error |
| Group `placement_mode`, `anchor`, and either `offset` or `offset_frac` | `absolute` or `box_anchor` | Shared `group_origin` for grouped templates | Combining `offset` with `offset_frac` is an error; specifying `anchor` in `absolute` mode is also an error |
| Template `group`, `center_local` | Requires `[mesh.groups.<name>]` | `center = group_origin + group_scale * center_local` | Combining with `center`, `placement_mode`, `anchor`, `offset`, `offset_frac`, `size_mode`, or `size_frac` is an error |
| Group `scale` or `scale_from` + `scale_factor` | `scale > 0`; `scale_from` names a box-size reference | Scale explicitly specified template length keys | Replaces the explicit dimensions with `scale * input`; dimensions left at defaults are not scaled |

The dimensions replaced by `size_mode="box_fraction"` depend on `kind`.

| `kind` | `size_frac` | Replaced keys |
|---|---|---|
| `plane`, `plane_hole`, `plate_hole` | float[2] | `size_x`, `size_y` |
| `box` | float[3] | `size` |
| `sphere` | float | `radius`, referenced to the minimum of the three box dimensions |
| `cylinder` | float[2] | `radius`, `height`; radius uses the minimum x/y box size and height uses z size |

Helper-parameter choices:

- Group scaling applies only to explicitly specified `size_x`, `size_y`, `size`, `radius`, `inner_radius`, and `height`.
- `anchor` accepts `box_center` or each axis's `*_low_face_center` / `*_high_face_center`.
- `scale_from` accepts `box_x`, `box_y`, `box_z`, `box_min_xy`, `box_max_xy`, `box_min_xyz`, or `box_max_xyz`.
- `placement_mode="absolute"`, `size_mode="absolute"`, and `inject_region_mode="absolute"` preserve direct values.

Run `beachx lint` from [Create and Validate `beach.toml`](Configuration.en.html) before execution.

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

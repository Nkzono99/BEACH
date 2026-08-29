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
| Developer run | `fpm run --target beach -- path/to/beach.toml` behaves the same |
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
| Length | `domain.box_min`, `domain.box_max`, `pos_low`, `pos_high` | m |
| Charge | `q_particle`, element charge output | C |
| Mass | `m_particle` | kg |
| Velocity | `drift_velocity`, `ray_direction` | m/s. `ray_direction` is a direction vector |
| Electric field | `e0`, `e0_abs` | V/m |
| Magnetic field | `b0` | T |
| Density | `number_density_cm3`, `number_density_m3` | cm^-3 or m^-3 |
| Temperature | `temperature_k`, `temperature_ev` | K or eV. They cannot both be specified |
| Angle | `e0_phi_xy_deg`, `e0_phi_z_deg` | degree |

`*_low` / `*_high` are lower and upper bounds on each axis. `inject_face` is one
of `x_low`, `x_high`, `y_low`, `y_high`, `z_low`, or `z_high`.

---

## Official Beginner Case

For the first run, use the [Ten-Minute Tutorial](Tutorial.en.html) and
`examples/tutorial_insulator.toml`. `beachx config init` generates the same
configuration. This beginner case uses `domain.periodic_axes=["x","y"]` and
`field_boundary.mode="periodic2"`.

Add boundary-reservoir inflow, `plane_source`, closed photoelectrons, and an infinite-periodic correction later from
[Design a Simulation Case](ConfigurationRecipes.en.html), after the beginner case
runs successfully.

---

## TOML Hierarchy and Section List

`[sim]`, `[domain]`, `[field_boundary]`, `[particle_boundary]`, `[reservoir]`,
`[surface_current_model]`, `[particles]`, `[mesh]`, `[periodic2]`, and `[output]` form the public configuration.

```text
beach.toml
├── [sim]
├── [domain]
├── [field_boundary]
├── [particle_boundary]
├── [reservoir]
├── [surface_current_model]
├── [particles]
│   └── [[particles.species]]       # one or more array-of-table entries
│       ├── [particles.species.boundary_inflow]
│       └── [particles.species.boundary]
├── [mesh]
│   ├── [mesh.groups.<name>]        # named child table
│   └── [[mesh.templates]]          # zero or more array-of-table entries
├── [periodic2]
└── [output]
```

Paths such as `sim.dt` and `domain.periodic_axes` mean “table name.key” in the prose.

| TOML table | Parent | Cardinality / requirement | Contents |
|---|---|---|---|
| `[sim]` | root | conditional | Time step, batch count, field solver, and external fields |
| `[domain]` | root | with a box | Box geometry and periodic axes |
| `[field_boundary]` | root | optional | `free` / `periodic2` field closure |
| `[particle_boundary]` | root | optional | Global particle actions on nonperiodic faces |
| `[reservoir]` | root | optional | External-reservoir inflow barrier and reference potential |
| `[surface_current_model]` | root | optional | External sheath closure resolving per-species `fixed_current` targets or a matching-plane response |
| `[particles]` | root | required | Container for `[[particles.species]]`; do not put ordinary keys directly under it |
| `[[particles.species]]` | `[particles]` | one or more | Species, injection mode, velocity distribution, macro-particle weight |
| `[particles.species.boundary]` | latest `[[particles.species]]` | optional | Nonperiodic-face overrides for that species |
| `[particles.species.boundary_inflow]` | latest `[[particles.species]]` | optional | Nonperiodic faces receiving external-reservoir inflow |
| `[mesh]` | root | optional | Selection of OBJ or built-in template input |
| `[mesh.groups.<name>]` | `[mesh]` | zero or more | Placement and scale shared by multiple templates |
| `[[mesh.templates]]` | `[mesh]` | zero or more | Built-in shapes used with `mode="template"` |
| `[periodic2]` | root | conditional | Nonzero mode, zero mode, and lower-boundary model for split periodic2 |
| `[output]` | root | optional | Output directory, history, checkpoint resume |

`[sim]` and `[domain]` are required when using `boundary_inflow`, `plane_source`, `reservoir_face`, or `photo_raycast`.
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
| `batch_duration` | float | `0.0` | Physical time per batch [s]. Under adaptive nonzero-mode progression, the maximum trial width of each accepted batch |
| `batch_duration_step` | float | `0.0` | Resolved as `batch_duration = dt * batch_duration_step` |
| `max_step` | int | `400` | Maximum number of pushes per particle |
| `tol_rel` | float | `1.0e-8` | Monitored relative-change value. Not used as a stop condition |
| `q_floor` | float | `1.0e-30` | Lower bound for the denominator in `rel_change` calculations |
| `multiple_box_events_policy` | string | `"abort"` | `abort` / `soft_discard` after the per-step boundary-event limit |
| `multiple_box_events_soft_discard_count_limit` | int | `1000` | Stop limit for cumulative soft discards |
| `multiple_box_events_soft_discard_abs_charge_limit` | float | `1.0e-12` | Stop limit for cumulative soft-discard absolute charge [C] |
| `raycast_max_bounce` | int | `16` | Maximum bounce count for `photo_raycast` |

Specifying both `batch_duration` and `batch_duration_step` is an error. For
`boundary_inflow` / `plane_source` / `reservoir_face` / `photo_raycast`, the resolved `batch_duration > 0` is required.

#### Field Solver

`field_solver` selects the method for computing the Coulomb electric field at
evaluation points from boundary-element charges. The table below lists the
corresponding parameters for each option.

| `field_solver` | Use case | Supported field boundary |
|---|---|---|
| `direct` | Exact all-to-all evaluation for small element counts and split references | `free`, or a constrained `periodic2` split reference |
| `treecode` | Approximate evaluation for medium and larger cases | `field_boundary.mode="free"` |
| `fmm` | Large-scale evaluation, `periodic2`, FMM core validation | `field_boundary.mode="free"` / `"periodic2"` |
| `auto` | Select direct / FMM based on element count | `field_boundary.mode="free"` |

Use the canonical [solver and field-boundary compatibility table](FieldSolvers.en.html#solver-and-field-boundary-compatibility) to choose the combination.
Element charge is always evaluated as a constant density on its triangle (a P0 panel); there is no configurable kernel.

Common keys:

| Key | Type | Default | Description |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | `direct` / `treecode` / `fmm` / `auto` |
| `field_normalization` | string | `"si"` | `si` / `box` / `mesh` / `length` |
| `field_length_scale` | float | `1.0` | Length scale used with `field_normalization="length"` [m] |
| `field_periodic_image_layers` | int | `1` | Near-image shell for `periodic2` |
| `field_periodic_far_correction` | string | `"none"` | `none` / `auto` / `cached_kneq0` |
| `field_periodic_ewald_alpha` | float | `0.0` | Ewald split used for cache generation; `0` selects automatically |
| `field_periodic_ewald_layers` | int | `4` | Real/reciprocal shell depth for cache generation |
| `field_periodic_cache_dir` | string | `".beach_cache/periodic2"` | Operator-cache directory |
| `field_periodic_generation_tolerance` | float | `1.0e-8` | Generation tolerance included in cache identity |

`field_periodic_far_correction="auto"` is accepted for compatibility and is
currently treated the same as `none`.

`field_normalization` only changes normalization of coordinates and
periodic cells inside field calculations. Output electric fields and potentials
are converted back to SI.

| `field_normalization` | Length scale |
|---|---|
| `si` | Use input SI coordinates as-is |
| `box` | Maximum width of `domain.box_max - domain.box_min`. Requires `[domain]` |
| `mesh` | Maximum width of the mesh bounding box. Falls back to `field_length_scale` if the mesh is empty |
| `length` | `field_length_scale` |

##### `field_solver = "direct"`

Sums all source elements directly. There is no approximation error, and the
computational complexity is `O(MN)`, where `M` is the number of evaluation points
and `N` is the number of elements.

| Key | Type | Default | Description |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | Specify `"direct"` |
| `field_normalization` | string | `"si"` | Normalize coordinates before direct evaluation |
| `field_length_scale` | float | `1.0` | Used with `field_normalization="length"` or mesh fallback |

`tree_theta`, `tree_leaf_max`, and `tree_min_nelem` are not used for `direct`.

##### `field_solver = "treecode"`

Builds a source octree. Distant nodes are evaluated with a monopole
approximation, and near nodes are evaluated with the analytic panel kernel. Unlike FMM, it
does not use local expansion and traverses the tree for each evaluation point.

| Key | Type | Default | Description |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | Specify `"treecode"` |
| `field_normalization` | string | `"si"` | Normalize coordinates before tree construction |
| `field_length_scale` | float | `1.0` | Used with `field_normalization="length"` or mesh fallback |
| `tree_theta` | float | `0.5` | MAC parameter. `0 < theta <= 1`. Larger values are faster and coarser |
| `tree_leaf_max` | int | `16` | Maximum number of sources per leaf node. `>= 1` |

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
| `field_normalization` | string | `"si"` | Normalize coordinates before FMM plan construction |
| `field_length_scale` | float | `1.0` | Used with `field_normalization="length"` or mesh fallback |
| `tree_theta` | float | `0.5` | MAC parameter for near/far classification. `0 < theta <= 1` |
| `tree_leaf_max` | int | `16` | Maximum number of sources per source-tree leaf node. `>= 1` |
`tree_min_nelem` is not used when `fmm` is specified explicitly.

##### `field_solver = "auto"`

Auto uses direct evaluation below `tree_min_nelem` and FMM at and above the threshold.

| Key | Type | Default | Description |
|---|---|---:|---|
| `field_solver` | string | `"auto"` | Specify `"auto"` |
| `field_normalization` | string | `"si"` | Common normalization used before automatic selection |
| `field_length_scale` | float | `1.0` | Used with `field_normalization="length"` or mesh fallback |
| `tree_min_nelem` | int | `256` | Element-count threshold for switching to FMM. `>= 1` |
| `tree_theta` | float | `0.5` | Near/far MAC parameter when FMM is selected |
| `tree_leaf_max` | int | `16` | Maximum number of sources per leaf node when FMM is selected |

If `tree_theta` and `tree_leaf_max` are not specified explicitly, the following
values are used based on the element count.

| Element count `nelem` | `tree_theta` | `tree_leaf_max` |
|---:|---:|---:|
| `< 1500` | `0.40` | `12` |
| `1500 <= nelem < 10000` | `0.50` | `16` |
| `10000 <= nelem < 50000` | `0.58` | `20` |
| `50000 <= nelem` | `0.65` | `24` |

### `[domain]`: Box Geometry and Periodic Topology

```toml
[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 4.0]
periodic_axes = ["x", "y"]
```

| Key | Type | Default | Description |
|---|---|---:|---|
| `box_min`, `box_max` | float[3] | none | Lower and upper box bounds [m]. Specify both |
| `box_origin`, `box_size` | float[3] | none | Origin and positive size [m]. Specify both |
| `periodic_axes` | string[] | `[]` | Unique entries from `"x"` / `"y"` / `"z"` |

Do not combine the two geometry forms. When `[domain]` is present, one complete pair is required.
`domain.periodic_axes` is the only public setting that selects periodicity.
With `field_boundary.mode="periodic2"`, `periodic_axes=["x","y"]` is required.

### `[field_boundary]`: Field Closure

```toml
[field_boundary]
mode = "periodic2"
```

| Key | Type | Default | Description |
|---|---|---:|---|
| `mode` | string | `"free"` | `free` / `periodic2` |

For `periodic2`, `[domain]` is required. Add `[periodic2]` when selecting the nonzero- and zero-mode operators explicitly.
Follow the
[solver and field-boundary compatibility table](FieldSolvers.en.html#solver-and-field-boundary-compatibility).

### `[particle_boundary]`: Global Particle Boundaries

```toml
[particle_boundary]
z_low = "open"
z_high = "open"
ordinary_open_model = "escape"
```

| Key | Type | Default | Description |
|---|---|---:|---|
| `x_low`, `x_high` | string | follow domain | `open` / `reflect` / `redistributed_reflect` on nonperiodic x faces |
| `y_low`, `y_high` | string | follow domain | `open` / `reflect` / `redistributed_reflect` on nonperiodic y faces |
| `z_low`, `z_high` | string | follow domain | `open` / `reflect` / `redistributed_reflect` on nonperiodic z faces |
| `ordinary_open_model` | string | `"escape"` | `escape` / `potential_barrier` on effective open faces |

An omitted face inherits domain topology and is open when nonperiodic.
This table cannot override a periodic face.
`reflect` reverses only normal velocity and preserves tangential velocity and the event position's tangential components.

`redistributed_reflect` applies the same velocity action and, for a single-face event, uniformly resamples both in-plane
coordinates over the box span excluding its end guards. At a simultaneous edge or corner event, only axes outside the event
mask are resampled.

`potential_barrier` evaluates `q_particle * (phi_infty - phi_boundary)` from the crossing-point potential
and `reservoir.phi_infty`; it reflects only a positive barrier that exceeds the normal kinetic energy.

### `[reservoir]`: External Reservoir Conditions

```toml
[reservoir]
inflow_model = "source_vdf"
phi_infty = 0.0
face_potential_grid_n = 3
```

| Key | Type | Default | Description |
|---|---|---:|---|
| `inflow_model` | string | `"source_vdf"` | `source_vdf` / `infinity_barrier` |
| `phi_infty` | float | `0.0` | Infinity-reference potential used by barriers [V] |
| `face_potential_grid_n` | int | `3` | `N x N` grid for injection-face average potential. `>=1` |

`infinity_barrier` adjusts the `boundary_inflow` normal VDF using the face-average potential and `phi_infty`.
It does not apply to an internal `plane_source` or deprecated `reservoir_face`.
A uniform field has no finite potential at infinity, so use an effective reservoir reference when combining it with this model.

### `[surface_current_model]`: External Sheath Closure

When this top-level table is omitted, configure each species manually with `target_*_current_a`.
`model="zhao_stationary"` solves the Zhao A/B/C zero-current stationary root for a planar, collisionless, unmagnetized
external sheath and resolves ambient-electron, ion, PE-emission, PE-escape, and PE-return currents plus a z-high kinetic
barrier once before batching.

`model="matching_plane_quasistatic"` couples the box top quasistatically to an outer sheath supplied by the selected
response backend. `response_backend="table"` is the default; `"zhao_online"` solves a finite-$H$, charge-driven Zhao
A/B/C response inside BEACH.

For `zhao_stationary`, `photoelectron_source_scale=0.0` creates no PE channels and uses the ambient-electron/ion
zero-current root.

```toml
[surface_current_model]
model = "zhao_stationary"
zhao_branch = "auto"
electron_species = "solar_wind_electron"
ion_species = "solar_wind_ion"
photoelectron_species = "photoelectron"
solar_elevation_deg = 60.0
photoelectron_ref_density_m3 = 6.4e7
photoelectron_source_scale = 1.0
# reference_area_m2 = 1.0e-8
```

For a no-photoelectron closure, omit every PE-specific key and explicitly set the source scale to zero.

```toml
[surface_current_model]
model = "zhao_stationary"
zhao_branch = "auto"
electron_species = "solar_wind_electron"
ion_species = "solar_wind_ion"
photoelectron_source_scale = 0.0
```

| Key | Type | Default | Description |
|---|---|---:|---|
| `model` | string | `"none"` | `none` / `zhao_stationary` / `matching_plane_quasistatic` |
| `response_backend` | string | `"table"` | Matching response source: `table` / `zhao_online` |
| `zhao_branch` | string | `"auto"` | `auto` / `a` / `b` / `c`; branch for stationary or online Zhao; no-PE stationary accepts only `auto` / `c` |
| `electron_species` | string | required | `species_key` for ambient electrons |
| `ion_species` | string | required | `species_key` for cold ions |
| `photoelectron_species` | string | required when PE is enabled | `species_key` for the `photo_raycast` tracking PE emission and return; omitting it disables PE in matching |
| `solar_elevation_deg` | float | required with stationary PE | Solar elevation $\alpha$ used by the Zhao source; $0<\alpha\le90$ degrees |
| `photoelectron_ref_density_m3` | float | required with stationary PE | Reference PE density $n_{pe,ref}$ [m^-3] |
| `photoelectron_source_scale` | float | `1.0` | Stationary-Zhao $s_{UV}$ in $n_{pe,0}=s_{UV}n_{pe,ref}\sin\alpha$; `0.0` disables PE |
| `reference_area_m2` | float | domain x-y area | Area converting Zhao current densities to total currents [m^2]; forbidden for matching |
| `response_table_path` | string | required for table matching | Outer-sheath response CSV v1; relative paths use the `beach.toml` directory and absolute paths are retained; forbidden online |
| `implicit_zero_mode` | bool | `false` | Apply backward Euler only to the matching-table mean $D_H$; requires `e_bottom_zero`, at least two $D_H$ nodes, and singleton feedback axes |
| `coupling_rtol` | float | `1.0e-4` | Relative matching fixed-point tolerance; finite $0<r\le1$ |
| `coupling_atol` | float[4] | `[0.0, 0.0, 0.0, 0.0]` | Per-feedback-component absolute tolerances, ordered as outward PE flux [m^-2 s^-1], PE mean normal energy [eV], outward electron flux [m^-2 s^-1], and outward ion flux [m^-2 s^-1]; values must be finite and nonnegative, with zero on inactive components |
| `coupling_max_iterations` | int | `20` | Maximum matching fixed-point iterations; `>=1` |
| `coupling_relaxation` | float | `0.5` | Matching update relaxation; finite $0<\omega\le1$ |

#### Zhao stationary closure

The referenced species must be enabled, distinct, and set `surface_charge_closure="fixed_current"`.
They cannot also specify manual `target_absorbed_current_a` or `target_emission_current_a`. Ambient electrons and ions must
enter inward from the z-high reservoir.

With `photoelectron_source_scale=0.0`, omit `photoelectron_species`, `solar_elevation_deg`, and
`photoelectron_ref_density_m3`, and set `zhao_branch` to `"auto"` or `"c"`. BEACH solves the Zhao Type C
$J_e+J_i=0$ closure and creates only the electron/ion absorption targets and z-high kinetic map; it creates no PE
emission, return, or escape targets.

When enabled, PE requires a negative `photo_raycast`, `inject_face="z_high"`, `deposit_opposite_charge_on_emit=true`, and an effective
open z-high particle boundary. All three roles must be singly charged, the ambient-electron and PE masses must match,
and $T_e>0$, $T_{pe}>0$, and $T_i\le0.1T_e$ must hold.

The external closure is unmagnetized, so `sim.b0=[0,0,0]` is also required. Because Zhao uses its own 0 V reservoir,
it cannot be combined with `reservoir.inflow_model="infinity_barrier"`.

The `ion_species` `number_density_*` supplies the ion density at infinity; the stationary root solves the ambient-electron
density. The configured electron density and PE `emit_current_density_a_m2` are sampling inputs for the raw Monte Carlo maps,
not fixed-current targets. The resolved electron density is recorded in `summary.txt`.

For signed current densities $J_e<0$, $J_i>0$, $J_{emit}>0$, and $J_{escape}>0$, PE return is
$J_{return}=J_{escape}-J_{emit}\le0$. BEACH multiplies them by area $A$, sends electron, ion, and return currents to their
absorption channels, and sends emission to the PE source-reaction channel.

The outward PE particle current $-AJ_{escape}$ is recorded as an external escape target. It is not deposited onto
surface elements. Its difference from the raw boundary-escape charge is recorded as an external-boundary closure
correction, so PE continuity and the stationary surface-current budget can be checked independently.

The element distribution and raw escape statistics measured by BEACH raycasting and trajectory tracking are not
overwritten. Outputs distinguish raw / target / applied for absorption and emission, and raw / target / applied /
correction for escape.

When Zhao is active, the selected ambient species z-high inflow VDF is energy-mapped from a 0 V infinity reservoir to the
current mean face potential. Type A electrons must also cross $\phi_m$ as an access bottleneck; Type B/C electrons and ions
use 0 V as the outside bottleneck.

The PE emission VDF remains the configured surface half-Maxwellian from `temperature_ev` and `normal_drift_speed`. At
outward z-high crossings, PE uses the remaining barrier to $\phi_m$ for Type A and to 0 V for
Type B/C: sub-barrier particles return and the energetic tail escapes.

Fixed-current scaling then preserves these raw maps while matching the Zhao channel totals. This model-specific map applies
only to the referenced species on z-high and is independent of the generic `reservoir.inflow_model`.

This model does not solve the field, space charge, Debye shielding, return orbit, or return delay outside the box, and it does
not update currents from the surface potential during batches. A return at z-high omits the distance and time to the outer
turning point. It is a fixed stationary external-current closure, not a self-consistent transient outer sheath. See
`examples/periodic2_zhao_fixed_current.toml` for a complete case.

#### Matching-plane quasistatic closure

Configure the table backend as follows. Omitting `response_backend` also selects `"table"`.

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "table"
response_table_path = "responses/outer_sheath.csv"
implicit_zero_mode = false
electron_species = "solar_wind_electron"
ion_species = "solar_wind_ion"
photoelectron_species = "photoelectron"
coupling_rtol = 1.0e-4
coupling_atol = [0.0, 0.0, 0.0, 0.0] # default: relative tolerance only
coupling_max_iterations = 20
coupling_relaxation = 0.5
```

Do not specify a response path for the online Zhao backend.

```toml
[surface_current_model]
model = "matching_plane_quasistatic"
response_backend = "zhao_online"
zhao_branch = "auto"
electron_species = "solar_wind_electron"
ion_species = "solar_wind_ion"
photoelectron_species = "photoelectron"
coupling_rtol = 1.0e-4
coupling_atol = [0.0, 0.05, 0.0, 0.0] # PE-energy tolerance [eV] for finite-ray sampling
coupling_max_iterations = 20
coupling_relaxation = 0.5
```

For matching without PE, omit `photoelectron_species` and the PE species. Online Zhao connects to the outer reservoir
with a zero PE population. For the table backend, set both PE flux and PE energy input axes to zero.

For each active feedback component $j$, with backend scale $s_j$, `coupling_rtol` $r$, and
`coupling_atol[j]` $a_j$, BEACH tests $|X_{raw,j}-X_j|\le\max(r s_j,a_j)$. The `0.05` eV above is an
example for finite-ray sampling; select it from a convergence study that varies the macro-particle count.

The third and fourth components are inactive for online Zhao because they represent outward ambient-electron and ion
fluxes. A singleton table-feedback axis is also inactive. BEACH rejects a nonzero absolute tolerance on any inactive
component.

It rescales an absolute-tolerance-dominated residual component so that an accepted trial still reports
`matching_plane_residual <= coupling_rtol`.

The matching plane $H$ is the z component of `domain.box_max`, and its coupling area is the domain x-y area. Every mesh
vertex must lie strictly below $H$.
`reference_area_m2` cannot replace that area.

For the table backend, `response_table_path` must point to a nonlinear response CSV v1 produced by an external Zhao
calculation or a 1D PIC model.
A synthetic table is suitable only for wiring and interpolation tests; it is not
a production-valid physics input. Resolve a relative path from the config-file directory and keep the resolved path at no
more than 256 characters. `beachx lint` applies the same resolved-path check.

For every query, the online backend uses $E_H=D_H/\epsilon_0$ as a boundary condition and solves a finite-$H$ Sagdeev
A/B/C root that connects $H$ to an upstream reservoir at 0 V with zero field. This is not the wall zero-current root of
`zhao_stationary`, and it does not impose $J=0$ while the surface charge is evolving.

`zhao_branch="auto"` searches the applicable branches; `"a"`, `"b"`, and `"c"` select one branch explicitly.
`auto` stops if it detects multiple physical roots or if a numerical failure in a compatible branch prevents a unique
selection. Its v1 multiple-root check clusters results from a finite multistart set rather than proving root isolation.

$H$ fixes the origin of the outer half-space interface, the zero-mode gauge, and the PE-moment measurement plane. Planar
translational symmetry makes its absolute coordinate inactive in the Sagdeev equation; this backend does not solve a
wall-to-$H$ distance constraint.

The outward PE number flux and mean normal energy are reduced to a half-Maxwellian that reproduces those two moments.
When the PE flux is zero, the PE population remains zero. The configured PE temperature supplies the numerical scale
when that role exists; otherwise the ambient-electron temperature is used.

The outward ambient-electron and ion axes are transparent and inactive in the fixed-point residual. The online solver is
stateless: it retains no outer inventory, previous-root continuation, outer flight time, or delayed-return queue.

If the configured branch policy finds no physical root or the solve does not converge, BEACH fails closed. It does not
silently change an explicitly selected branch or switch backend.

This model requires all of the following configuration invariants:

- `field_boundary.mode="periodic2"`, `domain.periodic_axes=["x","y"]`, and a z-nonperiodic open box;
- an explicit split `[periodic2]` configuration, with `nonzero_mode_backend` set to `cached_kneq0` or
  `panel_spectral_reference`, `zero_mode_policy="exclude_k0"`, and `lower_boundary_model` set to
  `e_bottom_zero` or `symmetric_vacuum`;
- `sim.e0=[0,0,0]`, `sim.b0=[0,0,0]`, no generic reservoir potential model, and
  `particle_boundary.ordinary_open_model="escape"`; `sim.multiple_box_events_policy` is `"abort"` or bounded
  `"soft_discard"` with count and absolute-charge limits;
- only distinct electron, ion, and optional photoelectron roles enabled, each with
  `surface_charge_closure="explicit"`;
- negative electron charge and positive ion charge; electron and ion roles use
  `source_mode="volume_seed"`, `npcls_per_step=0`, and only z-high `boundary_inflow="reservoir"`;
- when specified, a negative photoelectron role using `source_mode="photo_raycast"`, `inject_face="z_high"`, and
  `deposit_opposite_charge_on_emit=true`; and
- effective `periodic` x/y and `open` z-low/z-high particle boundaries for every role.
- no manual `fixed_current` target.

The backend-specific constraints are:

- the table backend requires `response_table_path` and forbids `zhao_branch`;
- the online backend forbids `response_table_path` and accepts
  `zhao_branch="auto"`, `"a"`, `"b"`, or `"c"`;
- the online backend requires singly charged species for every role, $T_e>0$, $0\le T_i\le0.1T_e$, positive ion
  number density, and positive inward ambient-electron and ion drift speeds at z-high (`drift_velocity` has a negative
  z component). When PE is specified, equal ambient-electron and PE masses and $T_{pe}>0$ are also required; and
- do not set the stationary-Zhao source keys `solar_elevation_deg`, `photoelectron_ref_density_m3`, or
  `photoelectron_source_scale` for matching.

Online Zhao is a reduced planar, collisionless, unmagnetized quasistatic closure. It does not solve a full VDF, 1D PIC,
a time-dependent outer sheath, or outer return of ambient outward populations.

With `model="none"`, do not specify any other key in the table.
The removed top-level `[outer_plasma]` and `[coupling]` tables
remain invalid; all coupling controls belong in `[surface_current_model]`. See
[Quasistatic Matching-Plane Coupling](MatchingPlaneCoupling.en.html) for the table CSV columns, online Zhao reduction,
fixed-point iteration, and validation.

### `[periodic2]`: Nonzero Mode, Zero Mode, and Lower Boundary

`[periodic2]` is a top-level table. Set `domain.periodic_axes=["x","y"]`
and `field_boundary.mode="periodic2"`.
Use `field_solver="fmm"` for production and reserve `field_solver="direct"` for the small split reference.

| Key | Default | Meaning |
|---|---:|---|
| `nonzero_mode_backend` | required | `panel_spectral_reference` / `cached_kneq0` |
| `zero_mode_policy` | required | `exclude_k0` |
| `lower_boundary_model` | required | `symmetric_vacuum` / `e_bottom_zero` |
| `max_nonzero_mode_potential_step` | `0` | Maximum accepted change of the $k\ne0$ potential [V]. Omission or `0` disables the feature |
| `reference_mode_layers` | `4` | Fourier-mode cutoff |
| `panel_quadrature_order` | `12` | panel-area quadrature order |

`max_nonzero_mode_potential_step > 0` is supported only with
`nonzero_mode_backend="cached_kneq0"`. Let $h_0$ be the resolved
`sim.batch_duration`. Each accepted batch tests the fixed ladder
$h_0,h_0/2,h_0/4,\ldots$ and accepts the first trial whose candidate charge
$\mathbf q_{\mathrm{candidate}}$ satisfies, at every panel centroid,

$$
\max_j\left|
\left[P_{k\ne0}
  \left(\mathbf q_{\mathrm{candidate}}-\mathbf q_{\mathrm{current}}\right)
\right]_j
\right|
\le
\texttt{max\_nonzero\_mode\_potential\_step}.
$$

Here, `P_{k\ne0}` is the cached nonzero-mode potential operator. A rejected
trial restores the RNG and macro-particle residuals to their pre-trial values. It does not contribute
to statistics, history, or the charge ledger.

This feature supports time-scaled `boundary_inflow`, `plane_source`, `reservoir_face`, and `photo_raycast`.
Flux-driven species must specify `target_macro_particles_per_batch`; fixed `w_particle` cannot be used.
It cannot be combined with an ordinary `volume_seed`.

`sim.batch_count` counts accepted batches, while `simulated_time_s` is the sum
of accepted widths. Dynamic teams are disabled during adaptive progression, and
MPI ranks use equal team sizes within a run. A restart may use a different team
size; the resumed run records its actual size as the new diagnostic value.

The test is a
local-voltage trust bound for the frozen-field approximation, not a
local-truncation-error guarantee. Verify time-width convergence by repeating
the calculation with this limit halved.

### Combined periodic2 Constraints

Periodic2 requires `[domain]`, `periodic_axes=["x","y"]`, and `field_boundary.mode="periodic2"`.
`examples/periodic2_closed_photoelectron.toml` is the reference combination of x/y periodicity, a boundary reservoir, and
closed photoelectrons. The same periodicity applies to field evaluation, collision, and `photo_raycast`.

| `nonzero_mode_backend` | Meaning |
| --- | --- |
| `panel_spectral_reference` | small-system split reference |
| `cached_kneq0` | production nonzero mode using a reusable versioned operator |

With `cached_kneq0`, the `exclude_k0` provider adds the physical zero mode exactly once.
`symmetric_vacuum` uses $\pm Q/(2\epsilon_0A)$ above and below; `e_bottom_zero` uses zero below and $Q/(\epsilon_0A)$ above.

### `[sim]`: External Fields

| Key | Type | Default | Description |
|---|---|---:|---|
| `e0` | float[3] | `[0,0,0]` | Uniform external electric field [V/m] |
| `e0_abs` | float | unspecified | Magnitude of the uniform external electric field [V/m] |
| `e0_phi_xy_deg` | float | `0.0` | Azimuthal angle in the xy plane when `e0_abs` is specified [deg] |
| `e0_phi_z_deg` | float | `0.0` | Elevation angle from the xy plane when `e0_abs` is specified [deg] |
| `b0` | float[3] | `[0,0,0]` | Uniform magnetic field [T] |

The uniform external electric field can be specified directly as
`e0 = [Ex, Ey, Ez]`, or by `e0_abs` and angles. Mixing the two forms is an error.

---

### `[[particles.species]]`: Particle Species

At least one `[[particles.species]]` entry is required. The keys and constraints
used depend on `source_mode`.

#### Common Keys

| Key | Type | Default | Description |
|---|---|---:|---|
| `enabled` | bool | `true` | Enable the species |
| `source_mode` | string | `"volume_seed"` | `volume_seed` / `plane_source` / `photo_raycast` / deprecated `reservoir_face` |
| `q_particle` | float | `-1.602176634e-19` | Particle charge [C] |
| `m_particle` | float | `9.10938356e-31` | Particle mass [kg] |
| `pos_low` | float[3] | `[-0.4,-0.4,0.2]` | Lower position bounds [m] |
| `pos_high` | float[3] | `[0.4,0.4,0.5]` | Upper position bounds [m] |
| `drift_velocity` | float[3] | `[0,0,-8e5]` | Drift velocity [m/s] |
| `temperature_k` | float | `2.0e4` | Temperature [K] |
| `temperature_ev` | float | unspecified | Temperature [eV]. Mutually exclusive with `temperature_k` |
| `velocity_distribution` | string | `"maxwellian"` | `maxwellian` / `grid` |
| `inject_face` | string | unspecified | Illumination-aperture face for `photo_raycast`. Also required by deprecated `reservoir_face` |
| `source_normal` | float[3] | unspecified | One-way `plane_source` normal. A nonzero axis-aligned vector |
| `boundary` | table | unspecified | Per-species six-face overrides in `[particles.species.boundary]` |
| `boundary_inflow` | table | unspecified | Per-species reservoir inflow faces in `[particles.species.boundary_inflow]` |
| `surface_charge_closure` | string | `"explicit"` | Surface-source charge closure. `explicit` / `fixed_current` / `neutral_return` |
| `target_absorbed_current_a` | float | unspecified | Signed absorbed target current [A] for `fixed_current`; its sign matches `q_particle` |
| `target_emission_current_a` | float | unspecified | Signed emission-reaction target current [A] for `fixed_current`; its sign opposes `q_particle` |

#### `[particles.species.boundary]`: Per-Species Overrides

This table belongs to the preceding `[[particles.species]]` entry.

```toml
[particles.species.boundary]
z_high = "reflect"
```

| Key | Type | Default | Description |
|---|---|---:|---|
| `x_low`, `x_high` | string | `"inherit"` | `inherit` / `open` / `reflect` / `redistributed_reflect` |
| `y_low`, `y_high` | string | `"inherit"` | `inherit` / `open` / `reflect` / `redistributed_reflect` |
| `z_low`, `z_high` | string | `"inherit"` | `inherit` / `open` / `reflect` / `redistributed_reflect` |

`inherit` uses the global action from `[particle_boundary]`. A periodic face cannot be overridden.
`inject_face` selects generation for `photo_raycast` and legacy `reservoir_face`; the species boundary controls the
trajectory after generation. `boundary_inflow` creates particles from outside and does not override outward actions.

`surface_charge_closure="neutral_return"` is accepted only for a negative
`photo_raycast` species with `deposit_opposite_charge_on_emit=true` and
an effective `reflect` or `redistributed_reflect` action on its `inject_face`. It scales resolved return deposits by one global
species factor so the photoelectron contribution to total surface charge is
zero.

It cannot be used when an actual escape or `soft_discard` occurs. The default
`"explicit"` commits ordinary tracked
charge without this closure. BEACH stops without correction when the unresolved
fraction exceeds the fixed 5% limit.

`surface_charge_closure="fixed_current"` preserves the element distribution measured by trajectory tracking while matching
the species total to an external current closure. Supply targets manually on each species or calculate them with the top-level
`[surface_current_model]`. Configure a manual absorbed channel as

```toml
surface_charge_closure = "fixed_current"
target_absorbed_current_a = -2.0e-6
```

BEACH multiplies raw absorbed charge $R_s$ by
$I_{s,\mathrm{abs}}^{\mathrm{target}}\Delta t/R_s$. For `photo_raycast` with
`deposit_opposite_charge_on_emit=true`, `target_emission_current_a` independently scales the source reaction channel. BEACH
does not scale the net PE current formed by subtracting these two channels. A nonzero target with an empty raw channel fails
closed.

`fixed_current` rescales the elementwise empirical distribution of the raw samples. If only one raw hit exists, that
element receives the entire target. The required hit count and acceptable scale depend on the mesh and the error tolerance
of the measured quantity, so BEACH does not impose a fixed minimum count or maximum scale.

Inspect `absorbed_count` / `emitted_count`, raw charge, and `fixed_*_weight_scale` in `charge_ledger.csv`, then test
convergence of the elementwise
charge distribution across macro-particle or `rays_per_batch` counts, batch widths, and RNG seeds. A small ledger residual
checks charge balance; it does not establish statistical accuracy of the spatial map.

`fixed_current` and `neutral_return` are mutually exclusive for one species. When an external model supplies a PE-return VDF
through a separate species at the top face, set `target_absorbed_current_a` on that return species and do not also use full
reflection or `neutral_return`; those combinations count the same return current twice.

#### `source_mode = "volume_seed"`

| Key | Type | Default | Description |
|---|---|---:|---|
| `npcls_per_step` | int | `0` | Number of macro particles generated per batch |
| `w_particle` | float | `1.0` | Macro-particle weight |

Constraints:

| Condition | Details |
|---|---|
| Particle count | Without boundary inflow, the sum of `npcls_per_step` over enabled species must be at least 1 |
| Automatic weight resolution | A `volume_seed` without `boundary_inflow` cannot use `target_macro_particles_per_batch` |

When the species has `boundary_inflow`, `npcls_per_step=0` is allowed. For a Maxwell distribution, a positive value combines
a volume seed with boundary inflow for the same species. Velocity-grid boundary inflow cannot use a positive value.

#### `[particles.species.boundary_inflow]`

This table belongs to the preceding `[[particles.species]]` entry.

```toml
[particles.species.boundary_inflow]
z_high = "reservoir"
```

| Key | Type | Default | Description |
|---|---|---:|---|
| `x_low`, `x_high` | string | omitted | `reservoir`; omission disables inflow |
| `y_low`, `y_high` | string | omitted | `reservoir`; omission disables inflow |
| `z_low`, `z_high` | string | omitted | `reservoir`; omission disables inflow |

`reservoir` injects the external VDF across the complete selected box face. It cannot be assigned to periodic faces.
Outward actions remain controlled independently by `[particle_boundary]` and `[particles.species.boundary]`, and each
inflow face must have an effective `open` action.

A Maxwell distribution uses `number_density_*` and temperature. A velocity grid uses `particle_flux_m2_s` or
`current_density_a_m2`. `w_particle` and `target_macro_particles_per_batch` are mutually exclusive.
Only `reservoir.inflow_model="infinity_barrier"` adjusts inflow using mean face potential and `phi_infty`.

With multiple faces, `target_macro_particles_per_batch` is the total target across all inflow faces, and one
macro-particle remainder is retained for each species-face pair.

The initial implementation can combine this table only with `source_mode="volume_seed"`. It cannot share a species with
`plane_source`, `photo_raycast`, or deprecated `reservoir_face`. The responsibilities remain separate to support a future
extension to multiple sources. See [reservoir inflow](ReservoirInjection.en.html) for flux and barrier details.

#### `source_mode = "plane_source"`

`plane_source` generates one-way flux along `source_normal` from an axis-aligned rectangle inside the box.
It uses the same particle-count and velocity-distribution keys as deprecated `reservoir_face`.

| Condition | Details |
|---|---|
| Domain | `[domain]` is required |
| Time | `sim.batch_duration > 0` is required |
| Rectangle | `pos_low` / `pos_high` are equal on exactly one axis and have positive extent on the other two |
| Placement | The normal coordinate lies strictly between box faces; tangential ranges stay inside and may reach box bounds |
| Direction | `source_normal` is nonzero along the zero-thickness axis. It is normalized internally; positive or negative unit input is recommended |
| External correction | `[reservoir]` settings `infinity_barrier`, `phi_infty`, and `face_potential_grid_n` do not apply |

For a Maxwell distribution, one-way flux is calculated from density and temperature. For a velocity grid, specify
`particle_flux_m2_s` or `current_density_a_m2`. Rectangle area, `batch_duration`, and `w_particle` determine the
macro-particle count, and the fractional remainder is carried to the next batch.

#### `source_mode = "reservoir_face"` (deprecated)

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
| Domain | `[domain]` is required |
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
`velocity_grid_path` is based on the runtime current directory.

The CSV contents are included in the ordered-species fingerprint. Changing a CSV in place therefore makes an existing
checkpoint incompatible instead of silently resuming with a different distribution.

The CSV is read once at its first use in a run. Later sampling and checkpoint fingerprints share that in-memory snapshot, so
replacing the file on disk cannot switch the distribution partway through a run. A new process reads the replacement as a new
physical model.

`velocity_distribution="grid"` follows the same face, time, and flux constraints as other `reservoir_face` sources.

The particle count is determined as follows.

```text
n_macro_expected = gamma_in * A * batch_duration / w_particle
n_injected = floor(residual + n_macro_expected)
```

The residual is carried over to the next batch. When
`target_macro_particles_per_batch > 0`, `w_particle` is calculated automatically
to approach that value.

This mode preserves existing-case behavior. BEACH does not silently convert it to `boundary_inflow` or `plane_source`.
Use `boundary_inflow` for a new external-plasma condition and `plane_source` for an internal rectangle.

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
| Domain | `[domain]` is required |
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
generated per batch is at most `rays_per_batch`. With `field_boundary.mode="periodic2"`,
emission starts from the hit coordinate wrapped to the primary cell even when a
periodic image is hit.

Each emitted photoelectron uses `w_hit` as its weight and is tracked as an
ordinary particle. Surface return is absorbed as an ordinary collision.

When a species boundary is `inherit`, it uses `[particle_boundary]`.
For closed PE, set its `inject_face` to `reflect` or `redistributed_reflect` in `[particles.species.boundary]`.

---

### `[mesh]`: Mesh Input

| Key | Type | Default | Description |
|---|---|---:|---|
| `mode` | string | `"template"` | `auto` / `obj` / `template` |
| `obj_path` | string | `"examples/simple_plate.obj"` | OBJ file path |
| `surface_model` | string | `"insulator"` | Surface model for the whole OBJ |
| `surface_side` | string | required with `mode="obj"` or `"auto"` | Vacuum side of OBJ panels: `normal_plus` / `normal_minus` / `outward_closed` |
| `obj_scale` | float | `1.0` | Uniform scale after loading the OBJ |
| `obj_rotation` | float[3] | `[0,0,0]` | Rotation angle after loading the OBJ [deg] |
| `obj_offset` | float[3] | `[0,0,0]` | Translation after loading the OBJ [m] |

With `mode="auto"`, an OBJ is used if `obj_path` exists; otherwise a template is
used. `surface_side` is still required so either branch has a complete OBJ
contract. The OBJ transformation order is `scale -> rotate -> offset`.

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
| `surface_model` | string | `"insulator"` | `insulator` / `conductor` |
| `surface_side` | string | required when `enabled=true` | Vacuum side of the panel: `normal_plus` / `normal_minus` / `outward_closed` |
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

`dielectric` and `epsilon_r` were removed from input because they only acted as metadata aliases without solving
polarization. Dielectric polarization, permittivity interface conditions, and internal fields remain unimplemented.

`conductor` constraints:

- Redistribute charge with direct Coulomb coefficients under `field_boundary.mode="free"`.
- Do not combine it with `field_boundary.mode="periodic2"`.
- Large conductor element counts increase per-batch dense-solve cost.

---

### Fixed element-source rules

Element sources use an implicit P0 triangle panel. The former `[field]` table
and `sim.softening` key have been removed; leaving either in an input fails as
an unknown table or key.

| Item | Rule |
| --- | --- |
| Source | Treat each `q_elem` as a constant surface-charge density over its triangle |
| Solver | `direct` / `treecode` / `fmm` / `auto`; `auto` selects direct / FMM with `tree_min_nelem` |
| Surface models | Common source discretization for `insulator` and `conductor` |
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
| `write_potential_history` | bool | `false` | Output `potential_history.csv`; with `[domain]`, also output same-batch `top_reference_history.csv` |
| `dir` | string | `"outputs/latest"` | Output directory |
| `history_stride` | int | `1` | Output interval for history CSV [batch] |
| `checkpoint_stride` | int | `0` | Restart-checkpoint interval [accepted batches]. `0` disables periodic output |
| `resume` | bool | `false` | Resume from an existing checkpoint |
| `restart_from` | string | none | Checkpoint source when `resume=true` |

Output files:

| File | Condition / contents |
|---|---|
| `summary.txt` | Run statistics and configuration summary |
| `charges.csv` | Final element charges |
| `mesh_triangles.csv` | Element geometry. Includes the `mesh_id` column |
| `mesh_sources.csv` | Original mesh kind, surface model, compatibility `epsilon_r` column, and element count; current inputs write `epsilon_r=1` |
| `mesh_potential.csv` | When `write_mesh_potential=true` |
| `charge_history.csv` | When `history_stride > 0` |
| `potential_history.csv` | When `write_potential_history=true` and `history_stride > 0` |
| `top_reference_history.csv` | Above plus `[domain]`; mean, standard deviation, minimum, and maximum potential across the full z-high face |
| `matching_plane_history.csv` | With `surface_current_model.model="matching_plane_quasistatic"` and `history_stride > 0`; selected accepted-batch coupling state and residual |
| `performance_profile.csv` | When `BEACH_PROFILE=1` |
| `rng_state.txt` / `rng_state_rankNNNNN.txt` | Serial or MPI rank-local random-number state |
| `macro_residuals.csv` | One MPI-global macro-particle residual file, distinguished by species and face |
| `charge_ledger.csv` | Per-species signed-charge flux, counts, and restartable cumulative values |
| `checkpoint_complete.txt` | Completion manifest for schema v8+ final output and each periodic slot |
| `checkpoint_latest.txt` | With `checkpoint_stride > 0`; advisory index of the latest normally published periodic slot |

See [Configuration-specific output](OutputGuide.en.html#locate-configuration-specific-values) to locate these values.

[Files Used for Resume](OutputGuide.en.html#files-used-for-resume) consolidates column definitions and conditional checkpoint requirements.

Evaluation rules for `mesh_potential.csv`:

- Record potential [V] at each element centroid.
- Evaluate the self term with the analytic P0 panel kernel.
- Add the explicit image shell for `periodic2`.
- Add the cached nonzero mode and boundary-specific `k=0` for `cached_kneq0`.

`potential_history.csv`:

- Use the same `history_stride` as `charge_history.csv`.
- Write `batch, elem_idx, potential_V`.
- Run `field_solver%refresh` and `compute_mesh_potential` for each history output, which increases computational cost.

`top_reference_history.csv`:

- Uses the same post-commit snapshot and batch as `potential_history.csv`.
- Writes `batch,simulated_time_s,z_high_m,sample_n,potential_mean_V,potential_std_V,potential_min_V,potential_max_V`.
- Uses a cell-centered grid over the full box z-high face with `sample_n=reservoir.face_potential_grid_n`.
- Records a relative-potential diagnostic, not infinity or plasma potential.

`matching_plane_history.csv`:

- Uses the same phase as `charge_history.csv`, writing accepted batch 1 and then every `history_stride` batches.
- Columns are `batch,simulated_time_s,D_H_C_m2,phi_H_V`, electron/ion inward fluxes and access potentials, the PE
  barrier potential, four outward feedback values, PE return/escape fluxes, and `iterations,residual`.
- $D_H$ and $\Phi_H$ describe the pre-commit state used for that batch's trajectories; `simulated_time_s` is after
  accepting the trial.

Requirements for `resume=true`:

| Condition | Details |
|---|---|
| Output | `write_files=true` is required |
| Source | If `restart_from` is unspecified, use `output.dir`; otherwise use `restart_from` |
| Required files | `summary.txt`, `charges.csv`, and either serial `rng_state.txt` or every MPI `rng_state_rankNNNNN.txt`; schema v8+ also requires `checkpoint_complete.txt` |
| Conditional files | `charge_ledger.csv` when ledger metadata is present |
| Conditional state | Schema v8+ requires `macro_residuals.csv` when declared by `checkpoint_complete.txt`; legacy schemas restore it when present |
| Behavior | If a required checkpoint is missing, stop instead of falling back to a new run |

`restart_from` changes only the checkpoint read source. New output is always written to `output.dir`.

With `checkpoint_stride > 0`, BEACH alternates between `checkpoints/slot0` and `slot1` after accepted-batch commit.
It atomically marks each directory complete with `checkpoint_complete.txt` only after all files are closed, then switches
`checkpoint_latest.txt` for periodic output. This preserves the previous slot if writing is interrupted.

On resume, BEACH compares the final output and both periodic slots below `output.dir` or `restart_from`, then selects the
loadable complete checkpoint with the largest `batches` value. It scans complete slot manifests directly when
`checkpoint_latest.txt` is missing, malformed, or stale. Final output remains restartable regardless of `checkpoint_stride`.

During MPI execution:

| File | Contents |
|---|---|
| `rng_state_rankNNNNN.txt` | Random-number state per rank |
| `macro_residuals.csv` | One global residual shared by all ranks and written by the root |

Resume consistency rules:

- Reject checkpoints with legacy `macro_residuals_rankNNNNN.csv` instead of converting them implicitly.
- Match `mpi_world_size` in `summary.txt` to the current rank count.
- Exclude checkpoint schemas unsupported by the current loader when recovering either periodic slot.
- Schema v2 and later require matching model, ordered-mesh, and ordered-species fingerprints.
- The model fingerprint also includes trajectory-contract versions for the boundary-event velocity aligned with the chord and
  discrete electric-field work, and for the surface-injection position. Checkpoints from
  the former contract are intentionally rejected instead of switching motion
  rules during a resumed run.
- `tree_theta` and `tree_leaf_max` include both their values and whether they were explicitly specified. Equal raw values still
  represent different solver contracts when one uses automatic estimates and the other uses explicit overrides.
- Schema v8 and later require `checkpoint_complete.txt` to match the summary batch, MPI world size, and conditional-file declarations.
- Schema v5 restores neutral-return correction, scale, and unresolved fraction from `charge_ledger.csv`.
- Schema v6 writes `macro_residuals.csv` as `species_idx,face,residual`. `face=0` denotes the legacy source and
  `1..6` denote boundary faces. The older two-column `species_idx,residual` form remains readable.
- `[[particles.species]].species_key` is stable. Omission yields `species_<1-based index>`; explicit values must be unique.

---

## Coordinate and Placement Helper Parameters

These are ordinary TOML parameters, but the loader uses them to calculate the physical coordinates or dimensions in the third
column. Distinguish combinations that fail validation from those that intentionally replace an explicit value.

| Key and value | Type / condition | Calculated value | Relationship to a directly specified key |
|---|---|---|---|
| `domain.box_origin`, `domain.box_size` | float[3]. Specify both, with `box_size > 0` | `box_min = box_origin`; `box_max = box_origin + box_size` | Combining either with `domain.box_min` or `domain.box_max` is an error |
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

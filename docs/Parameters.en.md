title: BEACH Input Parameters Reference

Lang: [English](Parameters.en.md) | [日本語](Parameters.md)

# Input Parameters Reference

This document is the parameter reference for `beach.toml` read by the Fortran runtime.
Unless otherwise noted, units are SI units.

This page lists every input key and retains its type, default, unit, range, required or mutually exclusive conditions,
and a one-sentence effect. Model derivations, algorithms, design rationale, output interpretation, and operational
examples are kept in dedicated explanations and guides linked from this reference.

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

Numbers and array components must be finite unless a key explicitly permits otherwise.
`*_low` / `*_high` are lower and upper bounds on each axis. `inject_face` is one
of `x_low`, `x_high`, `y_low`, `y_high`, `z_low`, or `z_high`.

---

## Official Beginner Case

For the first run, use the [Ten-Minute Tutorial](Tutorial.en.html) and
`examples/tutorial_insulator.toml`. `beachx config init` generates the same
configuration. This beginner case tracks 200 macro-electrons per batch for 20 surface-charge updates and uses
`field_solver="direct"` with `field_boundary.mode="free"`.

Add boundary-reservoir inflow, `plane_source`, periodic boundaries, closed photoelectrons, and an infinite-periodic correction later from
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
| `dt` | float | `1.0e-9` | Time step [s]. Finite and `>0` |
| `rng_seed` | int | `12345` | Random seed |
| `batch_count` | int | `1` | Batch count, or cumulative target count when resuming. `>=1` |
| `batch_duration` | float | `0.0` | Physical time per batch [s]. Finite; must be `>0` for flux-driven sources |
| `batch_duration_step` | float | `0.0` | Resolves `batch_duration=dt*batch_duration_step`. `>0`; mutually exclusive with `batch_duration` |
| `max_step` | int | `400` | Maximum pushes per particle. `>=1` |
| `tol_rel` | float | `1.0e-8` | Monitored relative-change value. Finite and `>=0`; not a stop condition |
| `q_floor` | float | `1.0e-30` | Denominator floor for `rel_change`. Finite and `>0` |
| `multiple_box_events_policy` | string | `"abort"` | `abort` / `soft_discard` after the per-step boundary-event limit |
| `multiple_box_events_retry_backend` | string | `"none"` | Retry after `multiple_box_events`: `none` / `upper_panel_fourier` |
| `multiple_box_events_soft_discard_count_grace` | int | `1000` | Count grace before enforcing the cumulative soft-discard fraction. Must be `>= 0` |
| `multiple_box_events_soft_discard_fraction_limit` | float | `1.0e-6` | Stop limit for the cumulative soft-discard fraction. Must satisfy `0 < value <= 1` |
| `multiple_box_events_soft_discard_abs_charge_limit` | float | `1.0e-12` | Cumulative soft-discard absolute-charge limit [C]. Finite and `>0` |
| `raycast_max_bounce` | int | `16` | Maximum `photo_raycast` bounce count. `>=1` when enabled |

Specifying both `batch_duration` and `batch_duration_step` is an error. For
`boundary_inflow` / `plane_source` / `reservoir_face` / `photo_raycast`, the resolved `batch_duration > 0` is required.

`upper_panel_fourier` is valid only for a `cached_kneq0` `periodic2` configuration. `soft_discard` stops after the count
grace when its cumulative fraction exceeds the fraction limit, or whenever cumulative absolute charge exceeds its
independent limit. See [Particle Events](ParticleEvents.en.html#stop-when-the-query-cannot-be-completed) for the retry
validity domain and diagnostics.

#### External Fields

| Key | Type | Default | Description |
|---|---|---:|---|
| `e0` | float[3] | `[0,0,0]` | Uniform external electric field [V/m] |
| `e0_abs` | float | unspecified | Uniform external electric-field magnitude [V/m]. `>=0` |
| `e0_phi_xy_deg` | float | `0.0` | Azimuth in the xy plane when `e0_abs` is specified [deg]. Explicit use requires `e0_abs` |
| `e0_phi_z_deg` | float | `0.0` | Elevation from the xy plane when `e0_abs` is specified [deg]. Explicit use requires `e0_abs` |
| `b0` | float[3] | `[0,0,0]` | Uniform magnetic field [T] |

Direct `e0` input is mutually exclusive with the `e0_abs` and angle form.

#### Field Solver

`field_solver` selects how Coulomb fields are evaluated from element charges.

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
| `field_length_scale` | float | `1.0` | Length used by `field_normalization="length"` [m]. `>0` |
| `field_periodic_image_layers` | int | `1` | Near-image shell for `periodic2`. `>=0`; `>=1` with `cached_kneq0` |
| `field_periodic_far_correction` | string | `"none"` | `none` / `auto` / `cached_kneq0` |
| `field_periodic_ewald_alpha` | float | `0.0` | Ewald split for cache generation. `>=0`; `0` selects automatically |
| `field_periodic_ewald_layers` | int | `4` | Real/reciprocal shell depth. `>=1` when far correction is enabled |
| `field_periodic_cache_dir` | string | `".beach_cache/periodic2"` | Operator-cache directory. Must be nonempty with `cached_kneq0` |
| `field_periodic_generation_tolerance` | float | `1.0e-8` | Cache-identity tolerance. Finite and `>0` with `cached_kneq0` |
| `tree_theta` | float | element-count heuristic below when omitted | Treecode/FMM MAC parameter. `0<theta<=1`; larger is faster and coarser |
| `tree_leaf_max` | int | element-count heuristic below when omitted | Maximum sources per treecode/FMM leaf. `>=1` |
| `tree_min_nelem` | int | `256` | Element count where `auto` switches from direct to FMM. `>=1` |

`field_periodic_far_correction="auto"` is accepted for compatibility and is
currently treated the same as `none`.

`field_normalization` changes only internal coordinates; output fields and potentials remain SI.

| `field_normalization` | Length scale |
|---|---|
| `si` | Use input SI coordinates as-is |
| `box` | Maximum width of `domain.box_max - domain.box_min`. Requires `[domain]` |
| `mesh` | Maximum width of the mesh bounding box. Falls back to `field_length_scale` if the mesh is empty |
| `length` | `field_length_scale` |

Mode-specific keys:

| `field_solver` | Evaluation | Additional keys used |
|---|---|---|
| `direct` | `O(MN)` reference sum over every source | none |
| `treecode` | Distant monopoles plus nearby panel kernels | `tree_theta`, `tree_leaf_max` |
| `fmm` | Approximate large-scale Coulomb FMM | `tree_theta`, `tree_leaf_max`, periodic2 keys when applicable |
| `auto` | Direct below `tree_min_nelem`; FMM at and above it | `tree_min_nelem`; `tree_theta` and `tree_leaf_max` when FMM is selected |

See [Use FMM](FMM.en.html) for selection and verification and
[Coulomb FMM core details](FMMCore.en.html) for implementation internals.

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
| `box_min`, `box_max` | float[3] | none | Lower and upper box bounds [m]. Specify both, with `box_max > box_min` on every axis |
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

An omitted face inherits domain topology and is `open` when nonperiodic. Periodic faces cannot be overridden.
See [Particle Events](ParticleEvents.en.html#open-reflect-redistributed_reflect-and-periodic-faces) for reflection,
redistributed reflection, and potential-barrier definitions.

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

`infinity_barrier` applies only to `boundary_inflow`. With a uniform field, set `phi_infty` to a consistent effective
reservoir reference. See [Inject Particles Through a Boundary](ReservoirInjection.en.html#3-choose-source_vdf-or-infinity_barrier)
for the correction and validity domain.

### `[surface_current_model]`: External Sheath Closure

When omitted, species-level settings provide manual current targets. Select a model only when using an external sheath
closure.

| `model` | Effect | Details |
|---|---|---|
| `none` | Do not use an external-sheath closure | Configure species `surface_charge_closure` |
| `zhao_stationary` | Resolve fixed currents and the z-high energy barrier from a Zhao stationary root | [Zhao Stationary Closure](ZhaoStationaryClosure.en.html) |
| `matching_plane_quasistatic` | Couple the box top quasistatically to an outer-sheath response | [Use Quasistatic Matching-Plane Coupling](MatchingPlaneCoupling.en.html) |

| Key | Type | Default | Description |
|---|---|---:|---|
| `model` | string | `"none"` | `none` / `zhao_stationary` / `matching_plane_quasistatic` |
| `response_backend` | string | `"table"` | Matching response source: `table` / `zhao_online` |
| `zhao_branch` | string | `"auto"` | `auto` / `a` / `b` / `c`; branch for stationary or online Zhao; no-PE stationary accepts only `auto` / `c` |
| `electron_species` | string | unspecified | Ambient-electron `species_key`; required for Zhao / matching; 1–64 characters |
| `ion_species` | string | unspecified | Cold-ion `species_key`; required for Zhao / matching; 1–64 characters |
| `photoelectron_species` | string | required when PE is enabled | PE `species_key`; 1–64 characters; omission disables PE in matching |
| `solar_elevation_deg` | float | required with stationary PE | Solar elevation $\alpha$ used by the Zhao source; $0<\alpha\le90$ degrees |
| `photoelectron_ref_density_m3` | float | required with stationary PE | Reference PE density $n_{pe,ref}$ [m^-3]. `>0` |
| `photoelectron_source_scale` | float | `1.0` | Stationary-Zhao $s_{UV}$. `>=0`; `0` disables PE |
| `reference_area_m2` | float | domain x-y area | Area converting Zhao current densities to total currents [m^2]. `>0`; forbidden for matching |
| `response_table_path` | string | required for table matching | Outer-sheath response CSV v1. Resolved length 1–256 characters; forbidden online |
| `implicit_zero_mode` | bool | `false` | Apply backward Euler only to the matching-table mean $D_H$; requires `e_bottom_zero`, at least two $D_H$ nodes, and singleton feedback axes |
| `coupling_rtol` | float | `1.0e-4` | Relative matching fixed-point tolerance; finite $0<r\le1$ |
| `coupling_atol` | float[4] | `[0.0, 0.0, 0.0, 0.0]` | Per-feedback-component absolute tolerances, ordered as outward PE flux [m^-2 s^-1], PE mean normal energy [eV], outward electron flux [m^-2 s^-1], and outward ion flux [m^-2 s^-1]; values must be finite and nonnegative, with zero on inactive components |
| `coupling_max_iterations` | int | `20` | Maximum matching fixed-point iterations; `>=1` |
| `coupling_relaxation` | float | `0.5` | Matching update relaxation; finite $0<\omega\le1$ |

#### Zhao stationary closure

Input constraints:

| Item | Required condition |
|---|---|
| Role species | Enabled, mutually distinct, and `surface_charge_closure="fixed_current"` |
| Ambient electron / ion | Enter inward from the z-high reservoir; no manual `target_*_current_a` |
| No PE | `photoelectron_source_scale=0.0`; omit PE-specific keys; use `zhao_branch="auto"` or `"c"` |
| With PE | Negative `photo_raycast`, `inject_face="z_high"`, `deposit_opposite_charge_on_emit=true`, effective z-high boundary `open` |
| Species properties | Singly charged; equal ambient-electron and PE masses; $T_e>0$, $T_{pe}>0$, $T_i\le0.1T_e$ |
| External field | `sim.b0=[0,0,0]`; `reservoir.inflow_model="infinity_barrier"` is forbidden |

Without PE, Type C produces only electron/ion absorption targets satisfying $J_e+J_i=0$ and the z-high kinetic-barrier
map; it produces no PE emission, return, or escape target.
`ion_species.number_density_*` is the ion density at infinity. Electron density and PE emission-current density are
sampling inputs; the closure determines current targets.

This is a stationary-current closure, not a transient outer-sheath solve. See
[Zhao Stationary Closure](ZhaoStationaryClosure.en.html) for current, barrier, PE-return, and output definitions.

The complete case is `examples/periodic2_zhao_fixed_current.toml`.

#### Matching-plane quasistatic closure

`response_backend="table"` uses an external response CSV; `"zhao_online"` uses the finite-$H$ Zhao response implemented
in BEACH. Without PE, omit `photoelectron_species`; table PE-flux and PE-energy input axes must also be zero. The matching
plane coordinate $H$ is the z component of `domain.box_max`, its area is the domain x-y area, and every mesh vertex must lie below $H$.

All rows below are required:

| Item | Required condition |
|---|---|
| Box / field | x/y periodic, z nonperiodic and open; `field_boundary.mode="periodic2"`; explicit `[periodic2]` split |
| Split | `nonzero_mode_backend` is `cached_kneq0` or `panel_spectral_reference`; `zero_mode_policy="exclude_k0"`; lower boundary is `e_bottom_zero` or `symmetric_vacuum` |
| External field / open face | `sim.e0=sim.b0=[0,0,0]`; `ordinary_open_model="escape"`; no generic reservoir-potential model |
| Event policy | `abort`, or `soft_discard` with fraction, count-grace, and absolute-charge limits |
| Roles | Only distinct, enabled electron / ion / optional PE roles; each uses `surface_charge_closure="explicit"` |
| Electron / ion source | Negative / positive charge; `source_mode="volume_seed"`; `npcls_per_step=0`; only z-high `boundary_inflow="reservoir"` |
| PE source | Negative `photo_raycast`; `inject_face="z_high"`; `deposit_opposite_charge_on_emit=true` |
| Particle boundaries | `periodic` x/y and `open` z-low/z-high for every role |
| Current targets | No manual `fixed_current` target |

| Backend | Required, forbidden, and species constraints |
|---|---|
| `table` | Requires `response_table_path`; forbids `zhao_branch` |
| `zhao_online` | Forbids `response_table_path`; accepts `zhao_branch` `auto` / `a` / `b` / `c` |
| `zhao_online` species | Every role singly charged; $T_e>0$; $0\le T_i\le0.1T_e$; positive ion density; negative z component of electron / ion `drift_velocity`; with PE, equal electron/PE masses and $T_{pe}>0$ |
| All matching backends | Forbid stationary-only `solar_elevation_deg`, `photoelectron_ref_density_m3`, and `photoelectron_source_scale` |

With `model="none"`, do not specify another key. Removed `[outer_plasma]` and `[coupling]` tables remain invalid.
See [Quasistatic Matching-Plane Coupling](MatchingPlaneCoupling.en.html) for model selection, physical meaning, and
applicability limits. The [matching-plane numerical and response-table reference](MatchingPlaneReference.en.html)
defines the CSV, implicit update, and fixed-point contracts.

### `[periodic2]`: Nonzero Mode, Zero Mode, and Lower Boundary

`[periodic2]` is a top-level table. Set `domain.periodic_axes=["x","y"]`
and `field_boundary.mode="periodic2"`.
Use `field_solver="fmm"` for production and reserve `field_solver="direct"` for the small split reference.

| Key | Type | Default | Meaning and constraints |
|---|---|---:|---|
| `nonzero_mode_backend` | string | required | `panel_spectral_reference` / `cached_kneq0` |
| `zero_mode_policy` | string | required | `exclude_k0` |
| `lower_boundary_model` | string | required | `symmetric_vacuum` / `e_bottom_zero` |
| `max_nonzero_mode_potential_step` | float | `0` | Accepted $k\ne0$ potential-change limit [V]. `>=0`; `0` disables; `cached_kneq0` only |
| `reference_mode_layers` | int | `4` | Fourier-mode cutoff. `>=1` |
| `panel_quadrature_order` | int | `12` | Panel-area quadrature order. `>=2` |

The adaptive potential limit supports `boundary_inflow`, `plane_source`, `reservoir_face`, and `photo_raycast`. The first
three require `target_macro_particles_per_batch` and forbid fixed `w_particle`; a positive `volume_seed` is also forbidden.
`sim.batch_count` counts accepted batches; `simulated_time_s` sums accepted widths. See
[How to choose `batch_duration`](BatchDurationStability.en.html)
for acceptance, rollback, and convergence checks.

### Combined periodic2 Constraints

Periodic2 requires `[domain]`, `periodic_axes=["x","y"]`, and `field_boundary.mode="periodic2"`.
`examples/periodic2_closed_photoelectron.toml` is the reference combination of x/y periodicity, a boundary reservoir, and
closed photoelectrons. The same periodicity applies to field evaluation, collision, and `photo_raycast`.

| `nonzero_mode_backend` | Meaning |
| --- | --- |
| `panel_spectral_reference` | small-system split reference |
| `cached_kneq0` | production nonzero mode using a reusable versioned operator |

See [periodic2 electrostatics](PeriodicElectrostatics.en.html) for the nonzero/zero-mode split, the role of
`exclude_k0`, and the mean field selected by the lower-boundary model.

### `[[particles.species]]`: Particle Species

At least one `[[particles.species]]` entry is required. The keys and constraints
used depend on `source_mode`.

#### Common Keys

| Key | Type | Default | Description |
|---|---|---:|---|
| `species_key` | string | `"species_<1-based index>"` | Stable ID. 1–64 characters and unique across species |
| `enabled` | bool | `true` | Enable the species |
| `source_mode` | string | `"volume_seed"` | `volume_seed` / `plane_source` / `photo_raycast` / deprecated `reservoir_face` |
| `q_particle` | float | `-1.602176634e-19` | Particle charge [C]. Nonzero for an enabled species |
| `m_particle` | float | `9.10938356e-31` | Particle mass [kg]. `>0` |
| `pos_low` | float[3] | `[-0.4,-0.4,0.2]` | Lower position bounds [m] |
| `pos_high` | float[3] | `[0.4,0.4,0.5]` | Upper position bounds [m] |
| `drift_velocity` | float[3] | `[0,0,-8e5]` | Drift velocity [m/s] |
| `temperature_k` | float | `2.0e4` | Temperature [K]. `>=0` |
| `temperature_ev` | float | unspecified | Temperature [eV]. `>=0`; mutually exclusive with `temperature_k` |
| `velocity_distribution` | string | `"maxwellian"` | `maxwellian` / `grid` |
| `inject_face` | string | unspecified | Illumination-aperture face for `photo_raycast`. Also required by deprecated `reservoir_face` |
| `source_normal` | float[3] | unspecified | One-way `plane_source` normal. A nonzero axis-aligned vector |
| `boundary` | table | unspecified | Per-species six-face overrides in `[particles.species.boundary]` |
| `boundary_inflow` | table | unspecified | Per-species reservoir inflow faces in `[particles.species.boundary_inflow]` |
| `surface_charge_closure` | string | `"explicit"` | Surface-source charge closure. `explicit` / `fixed_current` / `neutral_return` |
| `target_absorbed_current_a` | float | unspecified | Signed absorbed target current [A] for `fixed_current`; zero or the same sign as `q_particle` |
| `target_emission_current_a` | float | unspecified | Signed emission-reaction target current [A] for `fixed_current`; zero or the sign opposite `q_particle` |

Common keys for flux-driven sources (`boundary_inflow`, `plane_source`, and legacy `reservoir_face`):

| Key | Type | Default | Description and constraints |
|---|---|---:|---|
| `number_density_cm3` | float | unspecified | Upstream Maxwellian density [cm^-3]. `>0`; mutually exclusive with `number_density_m3` |
| `number_density_m3` | float | unspecified | Upstream Maxwellian density [m^-3]. `>0`; mutually exclusive with `number_density_cm3` |
| `w_particle` | float | exactly one of this and `target_macro_particles_per_batch` required | Macro-particle weight. `>0` |
| `target_macro_particles_per_batch` | int | exactly one of this and `w_particle` required | Macro-particle target. `>0`, or `-1` on species 2+ to share species 1 weight |
| `velocity_grid_path` | string | unspecified | Nonempty CSV path for `velocity_distribution="grid"`; columns `vx_m_s,vy_m_s,vz_m_s,f` |
| `velocity_grid_pdf_kind` | string | `"phase_space"` | `phase_space` / `flux_weighted` |
| `velocity_grid_sampling` | string | `"auto"` | `auto` / `rectilinear` / `discrete` |
| `particle_flux_m2_s` | float | unspecified | Incoming grid-distribution particle flux [m^-2 s^-1]. `>0`; mutually exclusive with `current_density_a_m2` |
| `current_density_a_m2` | float | unspecified | Incoming grid-distribution current density [A/m^2]. Nonzero and mutually exclusive with `particle_flux_m2_s` |

Specify exactly one of `w_particle` and `target_macro_particles_per_batch`.
For a Maxwell distribution, specify exactly one density form and a temperature. For a grid distribution, specify the CSV
and exactly one flux form. See [Inject Particles Through a Boundary](ReservoirInjection.en.html#2-choose-maxwell-or-velocity-grid-input)
for distribution and CSV-sampling semantics.

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

| `surface_charge_closure` | Effect | Required conditions |
|---|---|---|
| `explicit` | Commit tracked charge directly | Default |
| `neutral_return` | Correct the PE contribution to zero total surface-charge change | Negative `photo_raycast`; emission reaction enabled; injection face reflects; no escape / `soft_discard`; unresolved fraction $\le5\%$ |
| `fixed_current` | Preserve the sampled element distribution while scaling its total to a target current | Positive `batch_duration`; manual `target_*_current_a` or `[surface_current_model]`; a raw channel for every nonzero target |

`fixed_current` and `neutral_return` are mutually exclusive for one species. See
[Surface-Charge Update Numerical Specification](SurfaceChargeNumerics.en.html),
[Photoelectron Emission and Lifecycle](PhotoelectronEmission.en.html), and
[Output Format Reference](OutputReference.en.html#charge-ledger) for scaling, PE-return double-counting constraints, and
sample-count convergence checks.

#### `source_mode = "volume_seed"`

| Key | Type | Default | Description |
|---|---|---:|---|
| `npcls_per_step` | int | `0` | Macro particles generated per batch. `>=0` |
| `w_particle` | float | `1.0` | Macro-particle weight. `>0` |

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

`reservoir` injects across a complete selected nonperiodic box face, whose effective particle boundary must be `open`.
It combines only with `source_mode="volume_seed"`; on multiple faces, `target_macro_particles_per_batch` is the total
across all inflow faces. See [Inject Particles Through a Boundary](ReservoirInjection.en.html) for flux and barrier details.

#### `source_mode = "plane_source"`

`plane_source` generates flux along `source_normal` from an axis-aligned rectangle inside the box.

| Condition | Details |
|---|---|
| Domain | `[domain]` is required |
| Time | `sim.batch_duration > 0` is required |
| Rectangle | `pos_low` / `pos_high` are equal on exactly one axis and have positive extent on the other two |
| Placement | The normal coordinate lies strictly between box faces; tangential ranges stay inside and may reach box bounds |
| Direction | `source_normal` is nonzero along the zero-thickness axis. It is normalized internally; positive or negative unit input is recommended |
| External correction | `[reservoir]` settings `infinity_barrier`, `phi_infty`, and `face_potential_grid_n` do not apply |

Use the common flux-driven-source keys above for its flux and velocity distribution.

#### `source_mode = "reservoir_face"` (deprecated)

Use the common flux-driven-source keys above. Additional compatibility constraints:

| Condition | Details |
|---|---|
| Domain | `[domain]` is required |
| Time | `sim.batch_duration > 0` is required |
| Injection face | `inject_face` is required |
| Injection range | `pos_low` / `pos_high` must lie on the specified face |
| Weight | `w_particle` and `target_macro_particles_per_batch` cannot both be specified |
| Weight sharing | `target_macro_particles_per_batch=-1` is allowed only for species 2 and later. It shares species 1 `w_particle` |

For new cases, use `boundary_inflow` for an external plasma and `plane_source` for an internal rectangle. BEACH does not
silently convert the legacy mode. See
[Inject Particles Through a Boundary](ReservoirInjection.en.html#6-migrate-from-legacy-reservoir_face) for CSV and
weighting contracts.

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

See [Photoelectron Emission and Lifecycle](PhotoelectronEmission.en.html) for ray weighting, periodic images,
reabsorption, and closed-PE configuration.

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
| `history_stride` | int | `1` | History CSV interval [batch]. `>=0`; `0` disables |
| `checkpoint_stride` | int | `0` | Restart-checkpoint interval [accepted batches]. `>=0`; `0` disables periodic output; positive values require `write_files=true` |
| `resume` | bool | `false` | Resume from an existing checkpoint |
| `restart_from` | string | none | Checkpoint source when `resume=true`; requires `write_files=true` |

See [Output Format Reference](OutputReference.en.html) for file-generation conditions, columns, potential conventions,
matching-plane state, and ledger interpretation.

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

See [Files Used for Resume](OutputReference.en.html#files-used-for-resume) for periodic-slot selection, MPI-required files,
fingerprints, and schema compatibility.

---

## Coordinate and Placement Helper Parameters

These keys convert box-relative values to physical coordinates and dimensions.

| Key | Type | Default | Effect and constraints |
|---|---|---:|---|
| `domain.box_origin` | float[3] | unspecified | Sets `box_min`; pair with `box_size`; mutually exclusive with `box_min` / `box_max` |
| `domain.box_size` | float[3] | unspecified | Sets `box_max=box_origin+box_size`; every component `>0` |
| `inject_region_mode` | string | `"absolute"` | `absolute` / `face_fraction`; `face_fraction` only for `reservoir_face` / `photo_raycast` |
| `uv_low`, `uv_high` | float[2] | unspecified | Both required for `face_fraction`; components in `[0,1]`; mutually exclusive with `pos_low` / `pos_high` |
| Template `placement_mode` | string | `"absolute"` | `absolute` / `box_anchor` |
| Template `anchor` | string | unspecified | Box center or one of six face centers for `box_anchor` |
| Template `offset` | float[3] | unspecified | Offset from anchor [m]; mutually exclusive with `offset_frac` |
| Template `offset_frac` | float[3] | unspecified | Box-size-relative offset; mutually exclusive with `offset` |
| Template `size_mode` | string | `"absolute"` | `absolute` / `box_fraction` |
| Template `size_frac` | float / float[2] / float[3] | unspecified | Required for `box_fraction`; every component `>0`; replaces dimensions by kind as listed below |
| Template `group` | string | unspecified | Name in `[mesh.groups.<name>]`; nonempty |
| Template `center_local` | float[3] | unspecified | Required with `group`; `center=group_origin+group_scale*center_local` |
| Group `placement_mode` | string | `"absolute"` | `absolute` / `box_anchor` |
| Group `anchor` | string | unspecified | Box center or one of six face centers for `box_anchor` |
| Group `offset` | float[3] | unspecified | Offset added to group origin [m]; mutually exclusive with `offset_frac` |
| Group `offset_frac` | float[3] | unspecified | Box-size-relative offset; mutually exclusive with `offset` |
| Group `scale` | float | `1.0` | Multiplies group coordinates and explicit dimensions. `>0`; mutually exclusive with `scale_from` / `scale_factor` |
| Group `scale_from` | string | unspecified | Box-size reference; specify together with `scale_factor` |
| Group `scale_factor` | float | unspecified | Positive multiplier for `scale_from`. `>0` |

With `group`, a template cannot also use `center`, direct-placement keys, `size_mode`, or `size_frac`.

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

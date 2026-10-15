title: BEACH Input Parameters Reference

Lang: [English](Parameters.en.md) | [日本語](Parameters.md)

# Input Parameters Reference

Configuration files accept only the seven-group layout (`run / domain / mesh / particles / fields / sheath / output`). See `examples/` for complete input files.

Look up `beach.toml` keys, types, defaults, units, and required or mutually exclusive conditions here.
Units are SI unless stated otherwise. Start a new configuration with [Case design](ConfigurationRecipes.en.html),
then follow [Configuration validation](Configuration.en.html) after editing it.

| What to look up | Section |
|---|---|
| Run control and restart | [`run`](#run) |
| Domain geometry | [`domain`](#domain) |
| Surface mesh | [`mesh`](#mesh) |
| Particle tracking, sources and boundaries | [`particles`](#particles) |
| Field solver, periodic and external fields | [`fields`](#fields) |
| Outer sheath | [`sheath`](#sheath) |
| Output and diagnostics | [`output`](#output) |

You do not need to copy every entry into your configuration. Select the settings required by your sources, boundaries,
and solver; omitted settings use the listed defaults. For entries marked required or unspecified, also check the
applicability conditions below each table.

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
| Validation | `beachx lint beach.toml` checks keys, types, and ranges with the packaged schema, then combinations and derived values with `beach --check-config`. The Fortran preflight owns the combination constraints |
| Fortran validation | For development and diagnostics, `beach --check-config beach.toml` checks loading, normalization, and minimum model prerequisites; the path is required. It does not replace lint before running |
| Basic values | Reals and array components must be finite; integer fields do not accept reals or booleans |
| Integer storage range | Signed 32-bit `-2147483648..2147483647`; out-of-range values are rejected before conversion. Lint checks field-specific ranges |
| Strings | Values exceeding their Fortran destination length are errors; no implicit truncation |

To use an editor schema, put a comment directive with the GitHub Raw URL at the beginning of `beach.toml`.
The Fortran parser does not accept regular keys before the first section, so do not use `"$schema" = "..."`.

```toml
#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.schema.json
```

---

## Units and Coordinates

| Kind | Representative keys | Unit / direction |
|---|---|---|
| Time | `particles.tracking.dt_s`, `run.batch.duration_s` | seconds |
| Length | `domain.box_min`, `domain.box_max`, `source.pos_low`, `source.pos_high` | m |
| Charge | `charge_c`, element charge output | C |
| Mass | `mass_kg` | kg |
| Velocity | `distribution.drift_velocity_m_s`, `source.ray_direction` | m/s. `source.ray_direction` is a direction vector |
| Electric field | `fields.external.electric_v_m`, `fields.external.electric_magnitude_v_m` | V/m |
| Magnetic field | `fields.external.magnetic_t` | T |
| Density | `distribution.number_density_cm3`, `distribution.number_density_m3` | cm^-3 or m^-3 |
| Temperature | `distribution.temperature_k`, `distribution.temperature_ev` | K or eV. They cannot both be specified |
| Angle | `fields.external.electric_azimuth_deg`, `fields.external.electric_elevation_deg` | degree |

Numbers and array components must be finite. Disabled species still undergo checks on the types, finiteness, and string
lengths of supplied values, including the same integer storage range.
Schema / Python own detailed input diagnostics for enums, declared ranges, and settings supplied to inactive features.
Fortran does not repeat all these diagnostics; it resolves derived values and references and checks minimum prerequisites
for execution. Identical acceptance of invalid input is not guaranteed.
Physical combination checks follow each setting's applicability conditions.
`*_low` / `*_high` are lower and upper bounds on each axis. `source.inject_face` is one
of `x_low`, `x_high`, `y_low`, `y_high`, `z_low`, or `z_high`.

---

## Official Beginner Case

For the first run, use the [Ten-Minute Tutorial](Tutorial.en.html) and
`examples/tutorial_insulator.toml`. `beachx config init` generates the same
configuration. This beginner case tracks 200 macro-electrons per batch for 20 surface-charge updates and uses
`fields.solver.method="direct"` with `fields.boundary="free"`.

Add boundary-reservoir inflow, periodic boundaries, closed photoelectrons, and an infinite-periodic correction later from
[Design a Simulation Case](ConfigurationRecipes.en.html), after the beginner case
runs successfully.

---

## TOML Hierarchy and Section List

The public configuration consists of seven groups. Paths such as `run.batch.duration_s` mean “table.key”.

```text
beach.toml
├── [run]
│   ├── [run.batch]
│   │   └── [run.batch.adaptive]
│   └── [run.restart]
├── [domain]
├── [mesh]
│   ├── [mesh.obj]
│   ├── [mesh.groups.<name>]
│   └── [[mesh.templates]]
├── [particles]
│   ├── [particles.tracking]
│   │   └── [particles.tracking.events]
│   ├── [particles.raycast]
│   ├── [particles.boundary]
│   ├── [particles.reservoir]
│   ├── [particles.secondary_emission]
│   └── [[particles.species]]
│       ├── [particles.species.distribution]
│       ├── [particles.species.source]
│       ├── [particles.species.sampling]
│       ├── [particles.species.charging]
│       ├── [particles.species.boundary]
│       ├── [particles.species.inflow]
│       └── [particles.species.surface_interaction]
│           └── [particles.species.surface_interaction.secondary_emission]
├── [fields]
│   ├── [fields.solver]
│   │   └── [fields.solver.tree]
│   ├── [fields.periodic]
│   │   └── [fields.periodic.reference]
│   └── [fields.external]
├── [sheath]
│   ├── [sheath.zhao]
│   ├── [sheath.species]
│   ├── [sheath.photoelectrons]
│   └── [sheath.coupling]
└── [output]
    ├── [output.diagnostics]
    ├── [output.final]
    ├── [output.history]
    └── [output.checkpoint]
```

| TOML table | Parent | Cardinality / requirement | Contents |
|---|---|---|---|
| `[run]` | root | conditional | Random seed and batch count |
| `[run.batch]` | `[run]` | positive duration for flux sources | Physical time per batch |
| `[run.batch.adaptive]` | `[run.batch]` | optional; `cached_kneq0` only | Nonzero-mode potential-change limit |
| `[run.restart]` | `[run]` | present for restart | Checkpoint input source |
| `[domain]` | root | with a box | Box geometry and periodic axes |
| `[mesh]` | root | optional | Selection of OBJ or built-in templates |
| `[mesh.obj]` | `[mesh]` | with OBJ input | OBJ file, surface, and placement |
| `[mesh.groups.<name>]` | `[mesh]` | zero or more | Placement and scale shared by templates |
| `[[mesh.templates]]` | `[mesh]` | zero or more | Built-in shapes for `mesh.mode="template"` |
| `[particles]` | root | required | Species, tracking, boundaries, SEE generation control |
| `[particles.tracking]` | `[particles]` | optional | Particle time step and tracking limit |
| `[particles.tracking.events]` | `[particles.tracking]` | optional | Boundary-event retries and discards |
| `[particles.raycast]` | `[particles]` | optional | Photoelectron-ray bounce limit |
| `[particles.boundary]` | `[particles]` | optional | Global particle actions on nonperiodic faces |
| `[particles.reservoir]` | `[particles]` | optional | Reservoir inflow barrier and reference potential |
| `[particles.secondary_emission]` | `[particles]` | optional | Secondary-emission generation limit |
| `[[particles.species]]` | `[particles]` | one or more | Species identity, charge, and mass |
| `[particles.species.distribution]` | `[[particles.species]]` | when defining a distribution | Density, temperature, and drift |
| `[particles.species.source]` | `[[particles.species]]` | with an independent particle source | Generation or emission mode and position |
| `[particles.species.sampling]` | `[[particles.species]]` | with independent sources or boundary inflow | Macro count, weight, and ray count |
| `[particles.species.charging]` | `[[particles.species]]` | optional | Surface-current correction and targets |
| `[particles.species.boundary]` | `[[particles.species]]` | optional | Nonperiodic-face overrides for this species |
| `[particles.species.inflow]` | `[[particles.species]]` | optional | Nonperiodic reservoir inflow faces |
| `[particles.species.surface_interaction]` | `[[particles.species]]` | optional | Surface interaction on absorption |
| `[particles.species.surface_interaction.secondary_emission]` | `[particles.species.surface_interaction]` | optional | SEE settings on the incident species |
| `[fields]` | root | optional | `free` / `periodic2` field boundary |
| `[fields.solver]` | `[fields]` | optional | Field solver and length normalization |
| `[fields.solver.tree]` | `[fields.solver]` | optional | Treecode / FMM accuracy and subdivision |
| `[fields.periodic]` | `[fields]` | with periodic2 | Periodic backend and lower boundary |
| `[fields.periodic.reference]` | `[fields.periodic]` | with the split reference | Fourier cutoff and panel quadrature |
| `[fields.external]` | `[fields]` | optional | Uniform external electric and magnetic fields |
| `[sheath]` | root | optional | Closure from the outer-sheath zero-current root |
| `[sheath.zhao]` | `[sheath]` | with zero-current closure | Zhao branch and acceptance criteria |
| `[sheath.species]` | `[sheath]` | with zero-current closure | Electron, ion, PE, and SEE roles |
| `[sheath.photoelectrons]` | `[sheath]` | with zero-current closure | Outer-sheath photoelectron source |
| `[sheath.coupling]` | `[sheath]` | optional | Outer-root refresh from observed outflow |
| `[output]` | root | optional | Output enablement and directory |
| `[output.diagnostics]` | `[output]` | optional | Charge-change monitoring |
| `[output.final]` | `[output]` | optional | Final potential output |
| `[output.history]` | `[output]` | optional | History quantities and interval |
| `[output.checkpoint]` | `[output]` | optional | Periodic checkpoint interval |

---

## Detailed Parameter Reference

### `[run]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `rng_seed` | int | `12345` | Random seed |
| `batch_count` | int | `1` | Batch count, or cumulative target count when resuming. `>=1` |

#### `[run.batch]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `duration_s` | float | `0.0` | Physical time per batch [s]. Finite; must be `>0` for flux-driven sources |
| `duration_steps` | float | `0.0` | Resolves `run.batch.duration_s=particles.tracking.dt_s*run.batch.duration_steps`. `>0`; mutually exclusive with `run.batch.duration_s` |

Specifying both `run.batch.duration_s` and `run.batch.duration_steps` is an error. For
`inflow` / `photo_raycast`, the resolved `run.batch.duration_s > 0` is required.
See [Case design](ConfigurationRecipes.en.html#7-set-time-scales-and-particle-sampling) to choose steps, tracking limits,
and particle counts together.

##### `[run.batch.adaptive]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `max_nonzero_mode_potential_step_v` | float | `0` | Accepted $k\ne0$ potential-change limit [V]. `>=0`; `0` disables; `cached_kneq0` only |

The adaptive potential limit supports `inflow` and `photo_raycast`. Boundary inflow
requires `sampling.target_macro_particles_per_batch` and forbids fixed `sampling.weight`; a positive `volume_seed` is also forbidden.
`run.batch_count` counts accepted batches; `simulated_time_s` sums accepted widths. See
[How to choose `run.batch.duration_s`](BatchDurationStability.en.html)
for acceptance, rollback, and convergence checks.

#### `[run.restart]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `from` | string | `output.dir` | Checkpoint input directory; requires `output.enabled=true`. New output is written to `output.dir` |

Presence of `[run.restart]` enables restart; an empty table reads from `output.dir`.

`[run.restart]`:

| Condition | Details |
|---|---|
| Output | `output.enabled=true` is required |
| Source | If `run.restart.from` is unspecified, use `output.dir`; otherwise use `run.restart.from` |
| Required files | `summary.txt`, `charges.csv`, and either serial `rng_state.txt` or every MPI `rng_state_rankNNNNN.txt`; schema v8+ also requires `checkpoint_complete.txt` |
| Conditional files | `charge_ledger.csv` when ledger metadata is present |
| Conditional state | Schema v8+ requires `macro_residuals.csv` when declared by `checkpoint_complete.txt`; legacy schemas restore it when present |
| Behavior | If a required checkpoint is missing, stop instead of falling back to a new run |

`run.restart.from` changes only the checkpoint read source. New output is always written to `output.dir`.

See [Files Used for Resume](OutputReference.en.html#files-used-for-resume) for periodic-slot selection, MPI-required files,
mesh identity, and schema compatibility.

### `[domain]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `box_min`, `box_max` | float[3] | none | Lower and upper box bounds [m]. Specify both, with `box_max > box_min` on every axis |
| `box_origin`, `box_size` | float[3] | none | Origin and positive size [m]. Specify both; `box_max=box_origin+box_size`, all `box_size>0` |
| `periodic_axes` | string[] | `[]` | Unique entries from `"x"` / `"y"` / `"z"` |

Do not combine the two geometry forms. When `[domain]` is present, one complete pair is required.
`domain.periodic_axes` is the only public setting that selects periodicity.
With `fields.boundary="periodic2"`, `periodic_axes=["x","y"]` is required.

### `[mesh]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `mode` | string | `"template"` | `auto` / `obj` / `template` |

#### `[mesh.obj]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `path` | string | `"examples/simple_plate.obj"` | OBJ file path |
| `surface_model` | string | `"insulator"` | Surface model for the whole OBJ |
| `surface_side` | string | required with `mode="obj"` or `"auto"` | Vacuum side of OBJ panels: `normal_plus` / `normal_minus` / `outward_closed` |
| `scale` | float | `1.0` | Uniform scale after loading the OBJ |
| `rotation` | float[3] | `[0,0,0]` | Rotation angle after loading the OBJ [deg] |
| `offset` | float[3] | `[0,0,0]` | Translation after loading the OBJ [m] |

With `mode="auto"`, an OBJ is used if `mesh.obj.path` exists; otherwise a template is
used. `surface_side` is still required so either branch has a complete OBJ
contract. The OBJ transformation order is `scale -> rotate -> offset`.

```text
v_new = R(mesh.obj.rotation) * (v_old * mesh.obj.scale) + mesh.obj.offset
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

#### `[mesh.groups.<name>]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `placement_mode` | string | `"absolute"` | `absolute` / `box_anchor` |
| `anchor` | string | unspecified | Box center or one of six face centers for `box_anchor` |
| `offset` | float[3] | unspecified | Offset added to group origin [m]; mutually exclusive with `offset_frac` |
| `offset_frac` | float[3] | unspecified | Box-size-relative offset; mutually exclusive with `offset` |
| `scale` | float | `1.0` | Multiplies group coordinates and explicit dimensions. `>0`; mutually exclusive with `scale_from` / `scale_factor` |
| `scale_from` | string | unspecified | Box-size reference; specify together with `scale_factor` |
| `scale_factor` | float | unspecified | Positive multiplier for `scale_from`. `>0` |

#### `[[mesh.templates]]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `enabled` | bool | `true` | Enable the template |
| `kind` | string | `"plane"` | `plane` / `plate_hole` / `plane_hole` / `disk` / `annulus` / `box` / `cylinder` / `sphere` |
| `surface_model` | string | `"insulator"` | `insulator` / `conductor` |
| `surface_side` | string | required when `enabled=true` | Vacuum side of the panel: `normal_plus` / `normal_minus` / `outward_closed` |
| `center` | float[3] | `[0,0,0]` | Shape center [m] |
| `placement_mode` | string | `"absolute"` | `absolute` / `box_anchor` |
| `anchor` | string | unspecified | Box center or one of six face centers for `box_anchor` |
| `offset` | float[3] | unspecified | Offset from anchor [m]; mutually exclusive with `offset_frac` |
| `offset_frac` | float[3] | unspecified | Box-size-relative offset; mutually exclusive with `offset` |
| `size_mode` | string | `"absolute"` | `absolute` / `box_fraction` |
| `size_frac` | float / float[2] / float[3] | unspecified | Required for `box_fraction`; every component `>0`; replaces dimensions by kind as listed below |
| `group` | string | unspecified | Name in `[mesh.groups.<name>]`; nonempty |
| `center_local` | float[3] | unspecified | Required with `group`; `center=group_origin+group_scale*center_local` |

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

- Equalizes the potential returned by the field snapshot of the batch; available in free space and in `periodic2`.
- The influence matrix is assembled at initialization with one field refresh per conductor element. Large conductor element counts increase the initialization and dense-solve cost.

### `[particles]`

At least one `[[particles.species]]` entry is required.

#### `[particles.tracking]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `dt_s` | float | `1.0e-9` | Time step [s]. Finite and `>0` |
| `max_steps_per_particle` | int | `400` | Maximum pushes per particle. `>=1` |

##### `[particles.tracking.events]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `policy` | string | `"abort"` | `abort` / `soft_discard` after the per-step boundary-event limit |
| `retry_backend` | string | `"none"` | Retry after `multiple_box_events`: `none` / `upper_panel_fourier` |
| `discard_count_grace` | int | `1000` | Count grace before enforcing the cumulative soft-discard fraction. Must be `>= 0` |
| `discard_fraction_limit` | float | `1.0e-6` | Stop limit for the cumulative soft-discard fraction. Must satisfy `0 < value <= 1` |
| `discard_charge_warning_c` | float | `1.0e-12` | Warning threshold for cumulative soft-discard absolute charge [C]. Finite and `>0` |

Use these controls when adjusting boundary retries, particle discard, or photoelectron ray tracing.

`upper_panel_fourier` is valid only for a `cached_kneq0` `periodic2` configuration. `soft_discard` stops only after the
count grace when its cumulative fraction exceeds the fraction limit. Because cumulative absolute charge grows with
run duration, crossing its threshold emits a warning instead. See
[Particle Events](ParticleEvents.en.html#stop-when-the-query-cannot-be-completed) for the retry validity domain and
diagnostics.

#### `[particles.raycast]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `max_bounces` | int | `16` | Maximum `photo_raycast` bounce count. `>=1` when enabled |

#### `[particles.boundary]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `x_low`, `x_high` | string | follow domain | `open` / `reflect` / `redistributed_reflect` on nonperiodic x faces |
| `y_low`, `y_high` | string | follow domain | `open` / `reflect` / `redistributed_reflect` on nonperiodic y faces |
| `z_low`, `z_high` | string | follow domain | `open` / `reflect` / `redistributed_reflect` on nonperiodic z faces |

An omitted face inherits domain topology and is `open` when nonperiodic. Periodic faces cannot be overridden.
See [Particle Events](ParticleEvents.en.html#open-reflect-redistributed_reflect-and-periodic-faces) for reflection
and redistributed reflection definitions.

#### `[particles.reservoir]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `inflow_model` | string | `"source_vdf"` | `source_vdf` / `infinity_barrier` |
| `phi_infty_v` | float | `0.0` | Infinity-reference potential used by barriers [V] |
| `face_potential_grid_n` | int | `3` | `N x N` grid for injection-face average potential. `>=1` |

`infinity_barrier` applies only to `inflow`. With a uniform field, set `particles.reservoir.phi_infty_v` to a consistent effective
reservoir reference. See [Inject Particles Through a Boundary](ReservoirInjection.en.html#3-choose-the-inflow-mapping)
for the correction and validity domain.

#### `[particles.secondary_emission]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `max_generations` | int | `3` | `>=0`. Primaries are generation 0. Hits at the cap do not emit and report the would-be count |

### `[[particles.species]]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `species_key` | string | `"species_<1-based index>"` | Stable ID. 1–64 characters and unique across species |
| `enabled` | bool | `true` | Enable the species |
| `charge_c` | float | `-1.602176634e-19` | Particle charge [C]. Nonzero for an enabled species |
| `mass_kg` | float | `9.10938356e-31` | Particle mass [kg]. `>0` |

#### `[particles.species.distribution]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `drift_velocity_m_s` | float[3] | `[0,0,-8e5]` | Drift velocity [m/s] |
| `temperature_k` | float | `2.0e4` | Temperature [K]. `>=0` |
| `temperature_ev` | float | unspecified | Temperature [eV]. `>=0`; mutually exclusive with `distribution.temperature_k` |
| `number_density_cm3` | float | unspecified | Upstream Maxwellian density [cm^-3]. `>0`; mutually exclusive with `distribution.number_density_m3` |
| `number_density_m3` | float | unspecified | Upstream Maxwellian density [m^-3]. `>0`; mutually exclusive with `distribution.number_density_cm3` |

For boundary inflow, specify exactly one density form and a temperature; the inflow velocities follow a drifting Maxwellian.
See [Inject Particles Through a Boundary](ReservoirInjection.en.html#2-choose-the-distribution) for the distribution.

#### `[particles.species.source]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `mode` | string | required with `[particles.species.source]` | `volume_seed` / `photo_raycast` / `secondary_emission` |
| `pos_low` | float[3] | `[-0.4,-0.4,0.2]` | Lower position bounds [m] |
| `pos_high` | float[3] | `[0.4,0.4,0.5]` | Upper position bounds [m] |
| `inject_face` | string | unspecified | Illumination-aperture face for `photo_raycast` |
| `emit_current_density_a_m2` | float | `0.0` | Emission current density referenced to the ray-normal plane [A/m^2] |
| `deposit_opposite_charge_on_emit` | bool | `false` | Deposit opposite-sign charge on the emitting element |
| `normal_drift_speed_m_s` | float | `0.0` | Drift speed along the emission normal [m/s] |
| `ray_direction` | float[3] | inward normal of injection face | Ray direction |
| `region_mode` | string | `"absolute"` | `absolute` / `face_fraction`; `face_fraction` only for `photo_raycast` |
| `uv_low`, `uv_high` | float[2] | unspecified | Both required for `face_fraction`; components in `[0,1]`; mutually exclusive with `source.pos_low` / `source.pos_high` |

##### `mode = "photo_raycast"`

Constraints:

| Condition | Details |
|---|---|
| Domain | `[domain]` is required |
| Time | `run.batch.duration_s > 0` is required |
| Emission amount | `source.emit_current_density_a_m2 > 0`, `sampling.rays_per_batch > 0` are required |
| Injection face | `source.inject_face` is required |
| Particle properties | `charge_c` is nonzero, and `mass_kg > 0` |
| Ray direction | Must be normalizable, and its dot product with the inward normal of the injection face must be positive |
| Unavailable keys | `sampling.volume_macro_particles_per_batch`, `distribution.number_density_cm3` / `distribution.number_density_m3`, `sampling.weight`, `sampling.target_macro_particles_per_batch` |

See [Photoelectron Emission and Lifecycle](PhotoelectronEmission.en.html) for ray weighting, periodic images,
reabsorption, and closed-PE configuration.

#### `[particles.species.sampling]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `weight` | float | `1.0` | Macro-particle weight. `>0`; boundary inflow requires exactly one of this and `target_macro_particles_per_batch` |
| `target_macro_particles_per_batch` | int | exactly one of this and `sampling.weight` required | Macro-particle target. `>0`, or `-1` on species 2+ to share species 1 weight |
| `volume_macro_particles_per_batch` | int | `0` | Macro-particles generated by this species per batch. `>=0` |
| `rays_per_batch` | int | `0` | Number of launched rays per batch |

For boundary inflow, specify exactly one of `sampling.weight` and `sampling.target_macro_particles_per_batch`.

Constraints:

| Condition | Details |
|---|---|
| Particle count | Without boundary inflow, the sum of `sampling.volume_macro_particles_per_batch` over enabled species must be at least 1 |
| Automatic weight resolution | A `volume_seed` without `inflow` cannot use `sampling.target_macro_particles_per_batch` |

For boundary inflow alone, omit `source` and `sampling.volume_macro_particles_per_batch`; no volume population is generated.
To add volume generation with a Maxwell distribution, specify `source.mode="volume_seed"` and a positive
`sampling.volume_macro_particles_per_batch`. Volume generation cannot be combined with velocity-grid boundary inflow.

#### `[particles.species.charging]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `closure` | string | `"explicit"` | Surface-source charge closure. `explicit` / `fixed_current` / `neutral_return` |
| `target_absorbed_current_a` | float | unspecified | Signed absorbed target current [A] for `fixed_current`; zero or the same sign as `charge_c` |
| `target_emission_current_a` | float | unspecified | Signed emission-reaction target current [A] for `fixed_current`; zero or the sign opposite `charge_c` |

| `charging.closure` | Effect | Required conditions |
|---|---|---|
| `explicit` | Commit tracked charge directly | Default |
| `neutral_return` | Correct the PE contribution to zero total surface-charge change | Negative `photo_raycast`; emission reaction enabled; injection face reflects; no escape / `soft_discard`; unresolved fraction $\le5\%$ |
| `fixed_current` | Preserve the sampled element distribution while scaling its total to a target current | Positive `run.batch.duration_s`; manual `target_*_current_a` or `[sheath]`; a raw channel for every nonzero target |

`fixed_current` and `neutral_return` are mutually exclusive for one species. See
[Surface-Charge Update Numerical Specification](SurfaceChargeNumerics.en.html),
[Photoelectron Emission and Lifecycle](PhotoelectronEmission.en.html), and
[Output Format Reference](OutputReference.en.html#charge-ledger) for scaling, PE-return double-counting constraints, and
sample-count convergence checks.

#### `[particles.species.boundary]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `x_low`, `x_high` | string | `"inherit"` | `inherit` / `open` / `reflect` / `redistributed_reflect` |
| `y_low`, `y_high` | string | `"inherit"` | `inherit` / `open` / `reflect` / `redistributed_reflect` |
| `z_low`, `z_high` | string | `"inherit"` | `inherit` / `open` / `reflect` / `redistributed_reflect` |

This table belongs to the preceding `[[particles.species]]` entry.

`inherit` uses the global action from `[particles.boundary]`. A periodic face cannot be overridden.
`source.inject_face` selects generation for `photo_raycast`; the species boundary controls the
trajectory after generation. `inflow` creates particles from outside and does not override outward actions.

#### `[particles.species.inflow]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `x_low`, `x_high` | string | omitted | `reservoir`; omission disables inflow |
| `y_low`, `y_high` | string | omitted | `reservoir`; omission disables inflow |
| `z_low`, `z_high` | string | omitted | `reservoir`; omission disables inflow |

`reservoir` injects across a complete selected nonperiodic box face, whose effective particle boundary must be `open`.
For boundary inflow alone, omit `source` and `sampling.volume_macro_particles_per_batch`.
To add volume generation with a Maxwell distribution, specify `source.mode="volume_seed"` and a positive
`sampling.volume_macro_particles_per_batch`. Only a `volume_seed` source can be combined with inflow; on multiple faces,
`sampling.target_macro_particles_per_batch` is the total
across all inflow faces. See [Inject Particles Through a Boundary](ReservoirInjection.en.html) for flux and barrier details.

#### `[particles.species.surface_interaction]`

##### `[particles.species.surface_interaction.secondary_emission]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `enabled` | bool | `false` | emit after absorption |
| `model` | string | `"whipple"` | `whipple` / `constant` |
| `yield_max` | float | `1.5` | dimensionless; `>0` for enabled Whipple |
| `energy_max_ev` | float | `500.0` | [eV]; `>0` for enabled Whipple |
| `const_yield` | float | `0.0` | dimensionless; `>=0` for enabled constant |
| `target_species` | string | unset | required when enabled, resolves to an enabled `secondary_emission` species |

For `source.mode="secondary_emission"`, omit all `sampling.*`, `inflow.*`, `source.emit_current_density_a_m2`, ray count, and `charging.closure`.
Charge must be negative, mass positive, and temperature nonnegative; weight is inherited from the incident macro.
The sheath role additionally requires charge -e, the ambient-electron mass, and an open z-high face.
See [Secondary electron emission](SecondaryElectronEmission.en.html) for the model and provisional Phase 1 current allocation.

### `[fields]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `boundary` | string | `"free"` | `free` / `periodic2` |

For `periodic2`, `[domain]` is required. Add `[fields.periodic]` when selecting the nonzero- and zero-mode operators explicitly.
Follow the
[solver and field-boundary compatibility table](FieldSolvers.en.html#solver-and-field-boundary-compatibility).

#### `[fields.solver]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `method` | string | `"auto"` | `direct` / `treecode` / `fmm` / `auto` |
| `normalization` | string | `"si"` | `si` / `box` / `mesh` / `length` |
| `length_scale_m` | float | `1.0` | Length used by `fields.solver.normalization="length"` [m]. `>0` |
| `auto_fmm_min_elements` | int | `256` | Element count where `auto` switches from direct to FMM. `>=1` |

`fields.solver.method` selects how Coulomb fields are evaluated from element charges.

| `fields.solver.method` | Use case | Supported field boundary |
|---|---|---|
| `direct` | Exact all-to-all evaluation for small element counts and split references | `free`, or a constrained `periodic2` split reference |
| `treecode` | Approximate evaluation for medium and larger cases | `fields.boundary="free"` |
| `fmm` | Large-scale evaluation, `periodic2`, FMM core validation | `fields.boundary="free"` / `"periodic2"` |
| `auto` | Select direct / FMM based on element count | `fields.boundary="free"` |

Use the canonical [solver and field-boundary compatibility table](FieldSolvers.en.html#solver-and-field-boundary-compatibility) to choose the combination.
Element charge is always evaluated as a constant density on its triangle (a P0 panel); there is no configurable kernel.

currently treated the same as `none`.

`fields.solver.normalization` changes only internal coordinates; output fields and potentials remain SI.

| `fields.solver.normalization` | Length scale |
|---|---|
| `si` | Use input SI coordinates as-is |
| `box` | Maximum width of `domain.box_max - domain.box_min`. Requires `[domain]` |
| `mesh` | Maximum width of the mesh bounding box. Falls back to `fields.solver.length_scale_m` if the mesh is empty |
| `length` | `fields.solver.length_scale_m` |

Mode-specific keys:

| `fields.solver.method` | Evaluation | Additional keys used |
|---|---|---|
| `direct` | `O(MN)` reference sum over every source | none |
| `treecode` | Distant monopoles plus nearby panel kernels | `fields.solver.tree.theta`, `fields.solver.tree.leaf_max` |
| `fmm` | Approximate large-scale Coulomb FMM | `fields.solver.tree.theta`, `fields.solver.tree.leaf_max`, periodic2 keys when applicable |
| `auto` | Direct below `fields.solver.auto_fmm_min_elements`; FMM at and above it | `fields.solver.auto_fmm_min_elements`; `fields.solver.tree.theta` and `fields.solver.tree.leaf_max` when FMM is selected |

See [Use FMM](FMM.en.html) for selection and verification and
[Coulomb FMM core details](FMMCore.en.html) for implementation internals.

If `fields.solver.tree.theta` and `fields.solver.tree.leaf_max` are not specified explicitly, the following
values are used based on the element count.

| Element count `nelem` | `fields.solver.tree.theta` | `fields.solver.tree.leaf_max` |
|---:|---:|---:|
| `< 1500` | `0.40` | `12` |
| `1500 <= nelem < 10000` | `0.50` | `16` |
| `10000 <= nelem < 50000` | `0.58` | `20` |
| `50000 <= nelem` | `0.65` | `24` |

##### `[fields.solver.tree]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `theta` | float | element-count heuristic above when omitted | Treecode/FMM MAC parameter. `0<theta<=1`; larger is faster and coarser |
| `leaf_max` | int | element-count heuristic above when omitted | Maximum sources per treecode/FMM leaf. `>=1` |

#### `[fields.periodic]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `image_layers` | int | `1` | Near-image shell for `periodic2`. `>=0`; `>=1` with `cached_kneq0` |
| `backend` | string | `"finite_images"` | `finite_images` / `cached_kneq0` / `panel_spectral_reference` |
| `ewald_alpha` | float | `0.0` | Ewald split for cache generation. `>=0`; `0` selects automatically |
| `ewald_layers` | int | `4` | Real/reciprocal shell depth. `>=1` when far correction is enabled |
| `cache_dir` | string | `".beach_cache/periodic2"` | Operator-cache directory. Must be nonempty with `cached_kneq0` |
| `generation_tolerance` | float | `1.0e-8` | Cache-identity tolerance. Finite and `>0` with `cached_kneq0` |
| `lower_boundary_model` | string | required | `symmetric_vacuum` / `e_bottom_zero` |

`fields.boundary="periodic2"` requires `domain.periodic_axes=["x","y"]`. Use `fields.solver.method="fmm"` for production and `"direct"` only for the small split reference.

##### `[fields.periodic.reference]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `mode_layers` | int | `4` | Fourier-mode cutoff. `>=1` |
| `panel_quadrature_order` | int | `12` | Panel-area quadrature order. `>=2` |

#### Combined periodic2 Constraints

Periodic2 requires `[domain]`, `periodic_axes=["x","y"]`, and `fields.boundary="periodic2"`.
`examples/closed_photoelectron.toml` is the reference combination of x/y periodicity, a boundary reservoir, and
closed photoelectrons. The same periodicity applies to field evaluation, collision, and `photo_raycast`.

| `fields.periodic.backend` | Meaning |
| --- | --- |
| `panel_spectral_reference` | small-system split reference |
| `cached_kneq0` | production nonzero mode using a reusable versioned operator |

See [periodic2 electrostatics](PeriodicElectrostatics.en.html) for the nonzero/zero-mode split, the role of
`exclude_k0`, and the mean field selected by the lower-boundary model.

#### `[fields.external]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `electric_v_m` | float[3] | `[0,0,0]` | Uniform external electric field [V/m] |
| `electric_magnitude_v_m` | float | unspecified | Uniform external electric-field magnitude [V/m]. `>=0` |
| `electric_azimuth_deg` | float | `0.0` | Azimuth in the xy plane when `fields.external.electric_magnitude_v_m` is specified [deg]. Explicit use requires `fields.external.electric_magnitude_v_m` |
| `electric_elevation_deg` | float | `0.0` | Elevation from the xy plane when `fields.external.electric_magnitude_v_m` is specified [deg]. Explicit use requires `fields.external.electric_magnitude_v_m` |
| `magnetic_t` | float[3] | `[0,0,0]` | Uniform magnetic field [T] |

Direct `fields.external.electric_v_m` input is mutually exclusive with the `fields.external.electric_magnitude_v_m` and angle form.

#### Fixed element-source rules

| Item | Rule |
| --- | --- |
| Source | Treat each `q_elem` as a constant surface-charge density over its triangle |
| Solver | `direct` / `treecode` / `fmm` / `auto`; `auto` selects direct / FMM with `fields.solver.auto_fmm_min_elements` |
| Surface models | Common source discretization for `insulator` and `conductor` |
| Treecode | exact panel near + monopole far |
| FMM | exact panel near + exact triangle P2M |
| Surface side | `mesh.obj.surface_side` for OBJ; `surface_side` on every enabled template |

Use `outward_closed` only for consistently oriented, closed two-manifold components.

### `[sheath]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `closure` | string | `"none"` | `none` / `zero_current` |
| `reference_area_m2` | float | domain x-y area | Area converting current densities to total currents [m^2]; `>0` |

When omitted, species-level settings provide manual current targets. Select a model only when using an external sheath
closure.

| `sheath.closure` | Effect | Details |
|---|---|---|
| `none` | Do not use an external-sheath closure | Configure species `charging.closure` |
| `zero_current` | Resolve species currents and the z-high energy barrier from the Zhao zero-current root of an outer 1-D sheath | [Zhao Closure](ZhaoStationaryClosure.en.html) |

Input constraints:

| Item | Required condition |
|---|---|
| Role species | Enabled and mutually distinct; ambient electron, ion, and PE use `fixed_current`; SEE omits `charging.closure` |
| Ambient electron / ion | Enter from the z-high reservoir; the electron `distribution.drift_velocity_m_s` z component is at most 0 and the ion component is negative. Give electrons the same drift as the ions with `sheath.zhao.upstream_field_floor_v_m=2.0e-3` (not needed with zero electron drift; [why](ZhaoUpstreamBand.en.html)); no manual `target_*_current_a` |
| No PE | `sheath.photoelectrons.source_scale=0.0`; omit PE-specific keys; use `sheath.zhao.branch="auto"` or `"c"` |
| With PE | Negative `photo_raycast`, `source.inject_face="z_high"`, `source.deposit_opposite_charge_on_emit=true`, effective z-high boundary `open`. Startup warns (without stopping) when `source.emit_current_density_a_m2` (per ray-normal area) times the sine of the ray elevation differs by more than 5% from the outer source $e\,s_{UV}\,n_{pe,ref}\sin\alpha\,v_{th}/(2\sqrt\pi)$, or when the ray elevation differs from $\alpha$ by more than 1 degree |
| Species properties | Singly charged; equal ambient-electron and PE masses; $T_e>0$, $T_{pe}>0$, $T_i\ge0$ with a supersonic ion entry $m_i u_z^2/e >$ `sheath.zhao.ion_pressure_factor` $\cdot T_i$ |
| External field | `fields.external.magnetic_t=[0,0,0]`; `particles.reservoir.inflow_model="infinity_barrier"` is forbidden |
| Field | `fields.boundary="periodic2"`, x/y periodic axes, and an explicit split `[fields.periodic]` table (`fields.periodic.backend` set to `cached_kneq0` or `panel_spectral_reference`), so that the top-face mean potential can be fixed to the wall potential $\phi_0$ |
| Outer-root refresh | `outflow_refresh_batches>0` requires PE or SEE |

Without PE, Type C produces only electron/ion absorption targets satisfying $J_e+J_i=0$ and the z-high kinetic-barrier
map; it produces no PE emission, return, or escape target. SEE permits refresh without PE; the initial root alone uses a zero-source Type C. `distribution.number_density_m3` of the species selected by `sheath.species.ion` is the solar-wind density at infinity; write the same value in
`sheath.species.electron` (a different value is a configuration error). The electron inflow is injected at the upstream electron
Maxwellian density that the Zhao root solves from quasi-neutrality at infinity
([reservoir density](ZhaoStationaryClosure.en.html#inject-electrons-at-the-upstream-density)). PE emission-current density is a sampling input; the
closure determines current targets.

With a split `[fields.periodic]` table, the z-high plane-mean potential is fixed to the outer wall potential $\phi_0$; with x/y
periodic axes, particles returned by the outer barrier re-enter at a cell-uniform position. `outflow_refresh_batches>0`
is a weak coupling that quasi-statically re-solves the outer root from the observed PE + SEE outflow. BEACH does not solve a
transient outer sheath. See the [Zhao Closure](ZhaoStationaryClosure.en.html) for current, barrier, PE-return, refresh,

The complete case is `examples/zero_current.toml`; the outer-root refresh case is
`examples/zero_current_refresh.toml`.

#### `[sheath.zhao]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `branch` | string | `"auto"` | `auto` / `a` / `b` / `c`; without PE only `auto` / `c` |
| `upstream_band_tolerance` | float | `0.0` | For A/C roots with electron drift, the accepted width of negative $E^2$ next to upstream as a fraction of $\lvert\phi_m\rvert$ (A) or $\lvert\phi_0\rvert$ (C). `0<=x<1`; 0 disables this criterion, and with `sheath.zhao.upstream_field_floor_v_m` also 0 drifting A/C are rejected. The examples accept by the field floor instead of the width. See [Electron drift and the upstream band](ZhaoUpstreamBand.en.html) |
| `upstream_field_floor_v_m` | float | `0.0` | Places the upstream boundary where the field has fallen to this external level [V/m]. For A/C roots with electron drift, accepts the root when the field deficit $\sqrt{-\min E^2}$ inside the upstream band does not exceed it. The solar-wind motional field $\lvert\mathbf{u}\times\mathbf{B}\rvert$ (about 2 mV/m for 400 km/s and 5 nT) is the natural scale. `>=0`; 0 disables this criterion, and with `sheath.zhao.upstream_band_tolerance` also 0 drifting A/C are rejected. See [Electron drift and the upstream band](ZhaoUpstreamBand.en.html) |
| `ion_pressure_factor` | float | `1.0` | Pressure coefficient of the ion fluid in the outer root (Sana & Mishra 2026): 1 isothermal, 3 one-dimensional adiabatic. The ion temperature is that of `sheath.species.ion`; `>0`. For the solar wind (entry energy about 1 keV) it changes the root by less than 0.1 mV. See [Connecting to the outer sheath](ZhaoStationaryClosure.en.html#assumptions) |

#### `[sheath.species]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `electron` | string | unset | Ambient-electron `species_key`; required, 1–64 characters |
| `ion` | string | unset | Solar-wind ion `species_key`; required, 1–64 characters. Its temperature also enters the ion fluid of the outer root |
| `photoelectron` | string | required with PE | PE `species_key`; 1–64 characters |
| `secondary_electron` | string | unset | enabled SEE species; optional, 1–64 characters; fully tracked |

#### `[sheath.photoelectrons]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `solar_elevation_deg` | float | required with PE | Solar elevation $\alpha$ for the surface emission; $0<\alpha\le90$ degree |
| `ref_density_m3` | float | required with PE | PE reference density $n_{pe,ref}$ [m^-3]; `>0` |
| `source_scale` | float | `1.0` | $s_{UV}$; `>=0`, and 0 disables PE |

#### `[sheath.coupling]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `outflow_refresh_batches` | int | `0` | Re-solve the outer root every N accepted batches from the PE + SEE outflow spectrum observed at z-high; `>=0`, and 0 keeps the initial root |

### `[output]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `enabled` | bool | `true` | Enable/disable file output |
| `dir` | string | `"outputs/latest"` | Output directory |

See [Output Format Reference](OutputReference.en.html) for file-generation conditions, columns, potential conventions,
outer-sheath state, and ledger interpretation.

Requirements for

#### `[output.diagnostics]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `charge_rel_change_threshold` | float | `1.0e-8` | Monitored relative-change value. Finite and `>=0`; not a stop condition |
| `charge_floor_c` | float | `1.0e-30` | Denominator floor for `rel_change`. Finite and `>0` |

These settings diagnose charge changes. The run ends at `run.batch_count`; changing these values does not enable early stopping.

#### `[output.final]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `mesh_potential` | bool | `false` | Output `mesh_potential.csv` |

#### `[output.history]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `potential` | bool | `false` | Output `potential_history.csv`; with `[domain]`, also output same-batch `top_reference_history.csv` |
| `stride_batches` | int | `1` | History CSV interval [batch]. `>=0`; `0` disables |

#### `[output.checkpoint]`

| Key | Type | Default | Description / conditions |
|---|---|---:|---|
| `stride_batches` | int | `0` | Restart-checkpoint interval [accepted batches]. `>=0`; `0` disables periodic output; positive values require `output.enabled=true` |

## Coordinate and Placement Helper Parameters

The placement keys in `domain`, `particles.species.source`, `[mesh.groups.<name>]`, and `[[mesh.templates]]` convert box-relative values to physical coordinates and dimensions. Types and defaults are listed with each table above.

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
- `placement_mode="absolute"`, `size_mode="absolute"`, and `source.region_mode="absolute"` preserve direct values.

Run `beachx lint` from [Create and Validate `beach.toml`](Configuration.en.html) before execution.

---

## Validation Rules

| Item | Rule |
|---|---|
| Unknown keys | All are errors |
| `[particles]` | Contains species and `particles.secondary_emission.max_generations`; standalone injection of SEE is forbidden |
| Type | Validated by both the schema and the Fortran parser |
| Enums, ranges, and settings for inactive features | `beachx lint` provides detailed diagnostics; Fortran checks minimum prerequisites for execution |

Before running, this is recommended.

```bash
beachx lint beach.toml
```

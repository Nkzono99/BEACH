title: Migrate External-Boundary Configuration

Lang: [English](BoundaryConfigurationMigration.en.md) | [日本語](BoundaryConfigurationMigration.md)

# Migrate External-Boundary Configuration

This page maps legacy external-boundary keys under `[sim]`, `[outer_plasma]`, and `[coupling]` to the
`[external_boundary]` authoring facade. Legacy syntax remains accepted as compatibility input, while new examples and guides
use the facade as the canonical form.

## Migration Rules

- Do not put `[external_boundary]` and legacy `[outer_plasma]` / `[coupling]` in the same file.
- With the facade, do not set `sim.reservoir_potential_model`, `sim.sheath_injection_model`, or
  `sim.open_boundary_model`.
- `sim.phi_infty`, physical legacy-sheath parameters, box, species, periodic2, and solver settings remain usable.
- Legacy and new syntax lower to the same runtime contract; the authoring form is not included in the physical fingerprint.

Mixing remains an error even when both sides contain the same value. Last-wins behavior would hide configuration mistakes.

## Key Mapping

| Legacy setting | New setting |
| --- | --- |
| `outer_plasma.model` | `external_boundary.field.model` |
| `outer_plasma.kinetic_closure` | `external_boundary.field.kinetic_closure` |
| Other outer physical and diagnostic keys | Same key under `external_boundary.field` |
| Coupling time-scale, orbit, and steady-start keys | Same key under `external_boundary.particles` |
| `sim.open_boundary_model` | `external_boundary.ordinary_open.model` |
| `sim.reservoir_potential_model="infinity_barrier"` | `particles.inflow_model="infinity_barrier"` |
| `sim.sheath_injection_model=...` | `inflow_model="legacy_sheath"` plus `legacy_sheath_model=...` |
| `outer_plasma.return_model` | Removed; derived from field and particle mode |
| `coupling.particle_transfer_mode` | Removed; derived from `particles.mode` |
| `coupling.outer_queue_enabled` | Removed; derived from `particles.mode="zhao_queue"` |
| `coupling.update_mode` | Removed; currently fixed to `explicit` |
| `outer_plasma.interface_z` | Removed; derived from `sim.box_max[2]` |

## Scalar Inflow and Ordinary Open Faces

Legacy:

```toml
[sim]
reservoir_potential_model = "infinity_barrier"
open_boundary_model = "potential_barrier"
phi_infty = 0.0
```

New:

```toml
[sim]
phi_infty = 0.0

[external_boundary.field]
model = "none"

[external_boundary.particles]
mode = "local_source"
inflow_model = "infinity_barrier"

[external_boundary.ordinary_open]
model = "potential_barrier"
```

`infinity_barrier` owns inflow, while `potential_barrier` owns ordinary open faces. Either may still be enabled alone.

## Legacy Zhao Inflow Correction

Legacy:

```toml
[sim]
sheath_injection_model = "zhao_auto"
open_boundary_model = "escape"
sheath_alpha_deg = 60.0
sheath_photoelectron_ref_density_cm3 = 64.0
```

New:

```toml
[sim]
sheath_alpha_deg = 60.0
sheath_photoelectron_ref_density_cm3 = 64.0

[external_boundary.field]
model = "none"

[external_boundary.particles]
mode = "local_source"
inflow_model = "legacy_sheath"
legacy_sheath_model = "zhao_auto"

[external_boundary.ordinary_open]
model = "escape"
```

This remains a static source-VDF correction. It is not automatically converted to
`kinetic_closure="zhao_charge_driven"`.

## Tracked Kinetic 1-D Return

Legacy:

```toml
[outer_plasma]
model = "kinetic_1d"
kinetic_closure = "absorbing_maxwellian"
return_model = "kinetic_1d_profile_return"
interface_z = 1.0
debye_length = 0.2
thermal_voltage = 2.0

[coupling]
update_mode = "explicit"
particle_transfer_mode = "electrostatic_1d_instant_return"
field_evolution_timescale = 1.0
max_frozen_field_ratio = 0.1
outer_queue_enabled = false
```

New:

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "absorbing_maxwellian"
debye_length = 0.2
thermal_voltage = 2.0

[external_boundary.particles]
mode = "same_batch"
inflow_model = "auto"
field_evolution_timescale = 1.0
max_frozen_field_ratio = 0.1

[external_boundary.ordinary_open]
model = "escape"
```

Use `field.model="linear_debye"` with `particles.mode="local_source"` for linear field-only operation, or
`particles.mode="same_batch"` for linear 1-D return. The same distinction preserves kinetic field-only operation.

## Zhao Queue

Replace legacy `outer_queue_enabled=true` with:

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "zhao_charge_driven"
debye_length = 10.5
thermal_voltage = 10.0

[external_boundary.particles]
mode = "zhao_queue"
inflow_model = "auto"
field_evolution_timescale = 2.0e-5
max_frozen_field_ratio = 0.2
```

The queue requires automatic Zhao branch selection, update stride 1, and no legacy histogram. Fixed values are not repeated;
contradictory explicit values are errors.

## Unified 3-D Orbit

Replace the legacy matching explicit 3-D return/transfer pair with `same_batch`:

```toml
[external_boundary.field]
model = "unified_linear_response"
debye_length = 0.2
thermal_voltage = 10.0
unified_grid_points = 129

[external_boundary.particles]
mode = "same_batch"
inflow_model = "source_vdf"
field_evolution_timescale = 1.0e-4
max_frozen_field_ratio = 0.1
outer_orbit_dt = 1.0e-9
outer_orbit_max_steps = 10000
outer_orbit_energy_tolerance = 1.0e-3
```

The unified 3-D orbit does not own inflow, so `source_vdf`, `infinity_barrier`, or a legacy sheath may be selected separately.
Use `particles.mode="local_source"` for field-only operation.

After migration, run `beachx lint path/to/beach.toml`, then inspect `summary.txt` for the resolved inflow, ordinary-open rule,
interface transport, and particle mode.

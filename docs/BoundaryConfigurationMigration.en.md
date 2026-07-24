title: Migrate External-Boundary Configuration

Lang: [English](BoundaryConfigurationMigration.en.md) | [日本語](BoundaryConfigurationMigration.md)

# Migrate External-Boundary Configuration

This page explains how to move supported settings from legacy `[outer_plasma]` / `[coupling]` tables into the
`[external_boundary]` authoring facade and how to handle inputs that used removed external-boundary models.
Use `[external_boundary]` as the canonical syntax for new input.

## Current Choices

The public configuration vocabulary is limited to:

| Responsibility | Choices |
| --- | --- |
| `external_boundary.field.model` | `none` / `kinetic_1d` / `unified_linear_response` |
| `external_boundary.particles.mode` | `local_source` / `same_batch` / `zhao_queue` |
| `external_boundary.particles.inflow_model` | `auto` / `source_vdf` / `infinity_barrier` |
| `external_boundary.ordinary_open.model` | `escape` / `potential_barrier` |

`potential_barrier` handles outflow through ordinary open faces; `infinity_barrier` handles reservoir inflow.
They are independent choices. In a tracked `kinetic_1d` mode, one kinetic profile owns both inflow and return,
so use `inflow_model="auto"`.

## Removed Settings

The analytic linear Debye outer field, static sheath source correction, and the linear-model-only photoelectron histogram
have been removed. Their facade values and legacy runtime selector are errors and are not converted to another model.

| Removed input or artifact | Current behavior |
| --- | --- |
| `external_boundary.field.model="linear_debye"` / `outer_plasma.model="linear_debye"` | Configuration error; there is no alias to another field model |
| `external_boundary.particles.inflow_model="legacy_sheath"`, `legacy_sheath_model`, and `sim.sheath_injection_model` | Configuration error; the static source correction is not run |
| `external_boundary.field.infinity_potential` / `outer_plasma.infinity_potential` | Configuration error; the kinetic and unified infinity gauge is fixed internally at zero |
| `photoelectron_histogram_*` and `photoelectron_histogram.csv` | Configuration, output, and checkpoint state have been removed |

The following are choices by physical intent, not equivalent replacements:

| Previous intent | Current configuration to consider | Important difference |
| --- | --- | --- |
| Reduced 1-D field or return | `kinetic_1d` + `local_source` or `same_batch` | Solves species VDFs and a nonlinear Poisson profile; it does not reproduce the analytic linear solution |
| 3-D screening over a rough surface | `unified_linear_response` | A linear 3-D response, not a replacement for a 1-D sheath |
| Inflow without source-VDF correction | `field.model="none"` + `inflow_model="source_vdf"` | Uses the configured VDF directly as the face distribution |
| Inflow barrier from an explicit infinity potential | `field.model="none"` + `inflow_model="infinity_barrier"` | A scalar energy barrier, not a static current-balance closure |
| Zhao sheath closed by accumulated charge | `kinetic_1d` + `kinetic_closure="zhao_charge_driven"` | Refreshes the profile from surface charge each batch and closes inflow and return with the same profile |
| Linear-model-only photoelectron histogram | No direct replacement | Derive a needed distribution separately from particle or event output |

BEACH 1.5 or 1.6 configurations that use a removed value are errors in the current parser.
Checkpoints created with a removed model or a non-default removed feature cannot be resumed because their runtime contract and model fingerprint differ. Checkpoint-v4 runs using `kinetic_1d` or `unified_linear_response` with the retired features left at their former defaults remain fingerprint-compatible after the retired keys are removed from the configuration.
Start the new configuration from its initial state.

## Move Current Raw Settings into the Facade

- Do not put `[external_boundary]` and legacy `[outer_plasma]` / `[coupling]` in the same file.
- Move `sim.open_boundary_model` and `sim.reservoir_potential_model` into the facade.
- Box, species, periodic2, solver, `sim.phi_infty`, and other physical or numerical settings remain usable.
- Values from the two syntaxes are not merged. Mixed syntax is a configuration error.

| Legacy setting | Facade |
| --- | --- |
| `outer_plasma.model` | `external_boundary.field.model` |
| `outer_plasma.kinetic_closure` | `external_boundary.field.kinetic_closure` |
| Other current outer physics and diagnostic keys | Same key under `external_boundary.field` |
| Time-scale, orbit, and steady-start keys under `coupling` | Same key under `external_boundary.particles` |
| `sim.open_boundary_model` | `external_boundary.ordinary_open.model` |
| `sim.reservoir_potential_model="infinity_barrier"` | `external_boundary.particles.inflow_model="infinity_barrier"` |
| `outer_plasma.return_model` | Omit; derived from the field and particle mode |
| `coupling.particle_transfer_mode` | Omit; derived from `particles.mode` |
| `coupling.outer_queue_enabled` | Omit; derived from `particles.mode="zhao_queue"` |
| `coupling.update_mode` | Omit; current runtime is fixed to `explicit` |
| `outer_plasma.interface_z` | Omit; derived from `sim.box_max[2]` |

## Scalar Inflow and Ordinary Open Faces

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

`infinity_barrier` owns inflow; `potential_barrier` owns ordinary open faces. Either can also be used alone.

## Kinetic 1-D Tracked Return

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

Use `particles.mode="local_source"` for field-only operation, or `particles.mode="same_batch"` to close return and escape
with the kinetic profile.

## Zhao Queue

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

The queue requires `zhao_branch="auto"` and update stride 1. Do not repeat fixed values; contradictory explicit values are errors.

## Unified 3-D Orbit

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

A unified 3-D orbit does not own inflow, so select `source_vdf` or `infinity_barrier` separately.
Use `particles.mode="local_source"` for field-only operation.

After migration, run `beachx lint path/to/beach.toml`. After a run, inspect `summary.txt` for the resolved inflow,
ordinary-open rule, interface transport, and particle mode.

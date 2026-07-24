title: Configure the External Boundary

Lang: [English](OuterPlasmaModels.en.md) | [日本語](OuterPlasmaModels.md)

# Configure the External Boundary

This page configures the external plasma-response field, particles crossing z-high, and all other open faces.
Ordinary input edits only the tables under `[external_boundary]`. BEACH converts these choices into implementation-level
return, transfer, and queue settings, so users do not have to assemble those internal identifiers.

## Choose Three Responsibilities First

```text
[external_boundary]
├── [external_boundary.field]          external plasma-response field
├── [external_boundary.particles]      z-high particles and reservoir inflow
└── [external_boundary.ordinary_open]  open faces not owned by the outer model
```

`field` and `particles` are required. `ordinary_open` is optional and defaults to `escape`.

The standard self-consistent 1-D outer sheath is:

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "absorbing_maxwellian"
debye_length = 0.2
thermal_voltage = 2.0

[external_boundary.particles]
mode = "same_batch"
field_evolution_timescale = 1.0
```

BEACH derives kinetic-profile inflow mapping, 1-D return, and same-batch return from this input.
Do not specify `interface_z` or the explicit update mode: they resolve to `sim.box_max[2]` and `explicit`.
Omit `ordinary_open`, as above, when ordinary open faces simply escape.

Omit `[external_boundary]` entirely when all open faces simply escape and no inflow correction or external field is needed.
Otherwise, start with the row closest to the intended calculation.

| Goal | `field` | `particles.mode` | `inflow_model` |
| --- | --- | --- | --- |
| No external field or a scalar barrier | `none` | `local_source` | `source_vdf` / `infinity_barrier` |
| Standard self-consistent 1-D sheath | `kinetic_1d` + `absorbing_maxwellian` | `same_batch` | `auto` |
| Stationary or quasistationary Zhao sheath closed by accumulated charge | `kinetic_1d` + `zhao_charge_driven` | `same_batch` | `auto` |
| Zhao transient with outer-flight delay | `kinetic_1d` + `zhao_charge_driven` | `zhao_queue` | `auto` |

Add the following only when ordinary outflow needs a scalar barrier:

```toml
[external_boundary.ordinary_open]
model = "potential_barrier"
```

## Choose the External Field

`external_boundary.field.model` selects only the external plasma-response field.
The kinetic model solves beyond the interface. Whether particles enter the external region is
a separate choice under `external_boundary.particles.mode`.

| `field.model` | Position | Main use |
| --- | --- | --- |
| `none` | No external field | Ordinary open faces, the configured source VDF, and scalar barriers |
| `kinetic_1d` | Supported self-consistent 1-D sheath | Close the reservoir VDF, mean sheath, inflow, and return with one profile |

Field-only configurations remain supported:

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "absorbing_maxwellian"
debye_length = 0.2
thermal_voltage = 10.0

[external_boundary.particles]
mode = "local_source"
inflow_model = "source_vdf"
```

## Choose Interface-Particle Handling

| `particles.mode` | Particle leaving through z-high | Use |
| --- | --- | --- |
| `local_source` | Apply `ordinary_open`; retain no external trajectory | No external field, field-only, or scalar inflow |
| `same_batch` | Evaluate return with the kinetic 1-D profile | Tracked return in steady and quasi-steady studies |
| `zhao_queue` | Hold Zhao photoelectron return/escape events until their due batch | Delayed current during strong-UV turn-on |

`particles.mode` selects only external transport of z-high particles. Select the reservoir-inflow VDF independently with
`inflow_model` below.

Combine `same_batch` with `field.model="kinetic_1d"`. Return and escape are evaluated on the discrete kinetic 1-D profile.

`zhao_queue` is not a generic delayed transport option. It is available only with `kinetic_1d` and
`kinetic_closure="zhao_charge_driven"`.

## Choose Reservoir Inflow

`external_boundary.particles.inflow_model` controls the upstream VDF supplied to `reservoir_face`.

| `inflow_model` | Behavior |
| --- | --- |
| `auto` | Use the field profile for tracked 1-D models; otherwise use the configured source VDF |
| `source_vdf` | Treat the configured VDF as the face distribution |
| `infinity_barrier` | Correct accessibility and arrival speed from the mean face potential and `sim.phi_infty` |

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "zhao_charge_driven"
# Reference scales for split-interface diagnostics, not Zhao profile scales
debye_length = 10.5
thermal_voltage = 10.0

[external_boundary.particles]
mode = "same_batch" # use "zhao_queue" for the transient queue
field_evolution_timescale = 1.0
```

For photoemitting Zhao, the physical length and potential scales come from $T_{pe}$ and the reference density; the no-photo
case derives them from ambient $T_e$ and $n_\infty$. `debye_length` and `thermal_voltage` do not change the Zhao root or profile,
but remain required reference inputs for split-interface gap, lateral-field, and local-charge diagnostics.
Switching to `mode="zhao_queue"` also requires a positive `sim.batch_duration` and photoelectron source.
Do not add an update-stride key; queue mode fixes it internally to 1.

For tracked `kinetic_1d`, including `zhao_queue`, the same 1-D profile owns inflow.
Only `inflow_model="auto"` is accepted, preventing a second `infinity_barrier` correction.

## Choose Ordinary Open Faces

`external_boundary.ordinary_open.model` applies to open faces not owned by the outer model.

| `ordinary_open.model` | Behavior |
| --- | --- |
| `escape` | Remove the particle as permanent outflow |
| `potential_barrier` | Decide reflection or escape from local potential and normal kinetic energy |

`potential_barrier` has a different responsibility from the reservoir-inflow `infinity_barrier`. Either or both may be used.
When a tracked outer model owns z-high, `ordinary_open` still applies to the remaining open faces.

## Configurations Rejected as Errors

BEACH does not repair contradictory input or silently fall back to another model.

- `field.model="none"` with `same_batch` or `zhao_queue`
- A tracked kinetic 1-D profile with a local inflow correction
- `zhao_queue` with a non-Zhao closure
- `zhao_queue` with a manual branch or steady start
- A key with no effect for the selected field or particle mode, even when explicitly set to its default
- Unsupported magnetic field, species, periodic2, zero-mode, or time-scale choices
- Mixing `[external_boundary]` with legacy `[outer_plasma]` or `[coupling]`

Physical inputs and numerical guards are not guessed. Specify species, Debye length, temperature, field time scale,
and periodic2 backend as required by the selected model; contradictory values remain errors.
For `zhao_queue`, update stride is fixed internally to 1; do not write it in the public input.

## Inspect the Resolved Contract

`summary.txt` records the resolved inflow, ordinary-open rule, interface transport, and particle mode actually used at runtime,
not just the authoring syntax.

Use [External-Boundary Configuration Migration](BoundaryConfigurationMigration.en.html) to convert legacy
`[sim]` / `[outer_plasma]` / `[coupling]` input. Model physics and validation are documented separately in
[Outer Sheath: Kinetic 1-D](KineticOuterPlasma.en.html),
[Open Boundaries, Escape, and Return](ParticleEscapeReturn.en.html).

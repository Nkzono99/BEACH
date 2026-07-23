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
| Scalar barrier or static inflow correction | `none` | `local_source` | `source_vdf` / `infinity_barrier` / `legacy_sheath` |
| Reduced 1-D field or return reference | `linear_debye` | `local_source` / `same_batch` | select locally / `auto` when tracked |
| Standard self-consistent 1-D sheath | `kinetic_1d` + `absorbing_maxwellian` | `same_batch` | `auto` |
| Stationary or quasistationary Zhao sheath closed by accumulated charge | `kinetic_1d` + `zhao_charge_driven` | `same_batch` | `auto` |
| Zhao transient with outer-flight delay | `kinetic_1d` + `zhao_charge_driven` | `zhao_queue` | `auto` |
| Advanced linear 3-D response over a rough surface | `unified_linear_response` | `local_source` / `same_batch` | select locally |

Add the following only when ordinary outflow needs a scalar barrier:

```toml
[external_boundary.ordinary_open]
model = "potential_barrier"
```

## Choose the External Field

`external_boundary.field.model` selects only the external plasma-response field.
The linear and kinetic models solve beyond the interface; unified solves from
the rough surface to the far region. Whether particles enter the external region is
a separate choice under `external_boundary.particles.mode`.

| `field.model` | Position | Main use |
| --- | --- | --- |
| `none` | No external field | Ordinary open faces, scalar barriers, and legacy inflow corrections |
| `linear_debye` | Simplified 1-D reference | Field-only calculation or reduced 1-D return |
| `kinetic_1d` | **Standard and recommended** self-consistent 1-D sheath | Close the reservoir VDF, mean sheath, inflow, and return with one profile |
| `unified_linear_response` | Advanced, specialized 3-D linear response | Cases without a split window between roughness and plasma response |

`unified_linear_response` is not a higher-accuracy replacement for `kinetic_1d`. It does not solve species VDFs or floating
current balance; it constructs linear screening from a rough surface to the far region.

Field-only configurations remain supported:

```toml
[external_boundary.field]
model = "linear_debye"
infinity_potential = 0.0
debye_length = 0.2
thermal_voltage = 10.0

[external_boundary.particles]
mode = "local_source"
inflow_model = "source_vdf"
```

## Choose Interface-Particle Handling

| `particles.mode` | Particle leaving through z-high | Use |
| --- | --- | --- |
| `local_source` | Apply `ordinary_open`; retain no external trajectory | No external field, field-only, or scalar/legacy inflow |
| `same_batch` | Evaluate the field model's 1-D return or 3-D orbit | Tracked return in steady and quasi-steady studies |
| `zhao_queue` | Hold Zhao photoelectron return/escape events until their due batch | Delayed current during strong-UV turn-on |

`particles.mode` selects only external transport of z-high particles. Select the reservoir-inflow VDF independently with
`inflow_model` below.

The field model determines the concrete same-batch calculation.

| `field.model` | Resolved result for `particles.mode="same_batch"` |
| --- | --- |
| `linear_debye` | Analytic linear 1-D profile |
| `kinetic_1d` | Discrete kinetic 1-D profile |
| `unified_linear_response` | Explicit 3-D outer orbit |

`zhao_queue` is not a generic delayed transport option. It is available only with `kinetic_1d` and
`kinetic_closure="zhao_charge_driven"`.

## Choose Reservoir Inflow

`external_boundary.particles.inflow_model` controls the upstream VDF supplied to `reservoir_face`.

| `inflow_model` | Behavior |
| --- | --- |
| `auto` | Use the field profile for tracked 1-D models; otherwise use the configured source VDF |
| `source_vdf` | Treat the configured VDF as the face distribution |
| `infinity_barrier` | Correct accessibility and arrival speed from the mean face potential and `sim.phi_infty` |
| `legacy_sheath` | Select a static floating or Zhao VDF correction with `legacy_sheath_model` |

Configure legacy Zhao explicitly:

```toml
[external_boundary.field]
model = "none"

[external_boundary.particles]
mode = "local_source"
inflow_model = "legacy_sheath"
legacy_sheath_model = "zhao_auto"
```

This is a static source-VDF correction. A self-consistent Zhao outer sheath is a different configuration:

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
Do not add update-stride or histogram keys; queue mode fixes them internally.

For tracked `linear_debye` and `kinetic_1d`, including `zhao_queue`, the same 1-D profile owns inflow.
Only `inflow_model="auto"` is accepted, preventing a second `infinity_barrier` or legacy sheath correction.
A 3-D orbit does not own inflow, so `unified_linear_response + same_batch` may select local inflow separately.

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
- A tracked linear or kinetic 1-D profile with a local inflow correction
- `zhao_queue` with a non-Zhao closure
- `zhao_queue` with a manual branch, steady start, or unsupported histogram
- A key with no effect for the selected field or particle mode, even when explicitly set to its default
- Unsupported magnetic field, species, periodic2, zero-mode, or time-scale choices
- Mixing `[external_boundary]` with legacy `[outer_plasma]` or `[coupling]`

Physical inputs and numerical guards are not guessed. Specify species, Debye length, temperature, field time scale, orbit step,
and periodic2 backend as required by the selected model; contradictory values remain errors.
For `zhao_queue`, update stride is fixed internally to 1 and the histogram is disabled; write neither in the public input.

## Inspect the Resolved Contract

`summary.txt` records the resolved inflow, ordinary-open rule, interface transport, and particle mode actually used at runtime,
not just the authoring syntax. Legacy and new inputs that resolve to the same contract retain the same model fingerprint.

Use [External-Boundary Configuration Migration](BoundaryConfigurationMigration.en.html) to convert legacy
`[sim]` / `[outer_plasma]` / `[coupling]` input. Model physics and validation are documented separately in
[Outer Sheath: Kinetic 1-D](KineticOuterPlasma.en.html),
[Advanced Rough-Surface Linear Screening](UnifiedLinearResponse.en.html), and
[Open Boundaries, Escape, and Return](ParticleEscapeReturn.en.html).

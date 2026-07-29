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
Do not specify `interface_z` or an update mode. The interface is derived from `sim.box_max[2]`; update mode is normally
`explicit` and becomes internal `implicit_mean` for the ambient-linear or strong-photoemission Zhao compositions described
below.
Omit `ordinary_open`, as above, when ordinary open faces simply escape.

Omit `[external_boundary]` entirely when all open faces simply escape and no inflow correction or external field is needed.
Otherwise, start with the row closest to the intended calculation.

| Goal | `field` | `particles.mode` | `inflow_model` |
| --- | --- | --- | --- |
| No external field or a scalar barrier | `none` | `local_source` | `source_vdf` / `infinity_barrier` |
| Standard self-consistent 1-D sheath | `kinetic_1d` + `absorbing_maxwellian` | `same_batch` | `auto` |
| Linear response separating local weak-photoemission structure from mean charging | `kinetic_1d` + `ambient_linear_debye` | `same_batch` | `auto` |
| Empirical-CDF strong-photoemission quasisteady sheath with a virtual cathode | `kinetic_1d` + `zhao_charge_driven` | `same_batch` | `auto` |
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

For weak photoemission under a linear-response approximation, select:

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "ambient_linear_debye"
debye_length = 0.2
thermal_voltage = 2.0

[external_boundary.particles]
mode = "same_batch"
field_evolution_timescale = 1.0
```

This closure responds with the ambient plasma only,
$\phi(z)=\lambda_D E_I\exp(-z/\lambda_D)$.
Do not write `photoelectron_density_model` in public TOML. The facade resolves its internal value to `none` and rejects the
explicit key even when its value is `"none"`.
With an enabled negative `photo_raycast` species, BEACH automatically selects internal `implicit_mean` without adding a
public key. The initial 3-D trace holds batch-start $k\ne0$ fixed and records emission, reabsorption before the interface, and
ambient absorption on their elements. Photoelectron z-high crossings are deferred until the updated mean field is available.
Their measured outward current and tracked ambient-absorption current enter an analytic-Maxwellian backward-Euler solver that
determines total charge, $\phi_I$, and continuous escape/return fractions. Configured `emit_current_density_a_m2` determines
tracked 3-D emission weights but is not reused as an independent PE mean source.

After the mean solve, BEACH traces each full-weight crossing ray once with energy resolution on the kinetic 1-D profile.
A return includes outer flight time and lateral displacement and is followed to a terminal state under mean $k=0$ and
batch-start $k\ne0$. A ray that recrosses z-high after return is sent through the same profile mapper again.

Analytic return charge is normalized over the source legs (temporary neutralization at emission elements) and actual-hit
destination legs of rays that are ultimately reabsorbed, using one `return_weight_scale`. The transaction remains zero-sum;
pending emission countercharge is retained and
the destination deposits are committed once.
Mean total charge and retarding potential update within the current batch; sampled local redistribution changes $k\ne0$ in
the next batch. Discrete ray classification is therefore not iterated back into the mean solve.

A mean-solver failure, positive analytic return with no reabsorbed sample, charge mismatch, or a ray that does not terminate
as absorbed or escaped within the allowed trace stops the run. Deferred PE rays are quasistationary shadows: flight time and
frozen-field ratio remain diagnostic, but an over-limit ratio does not stop the run. Ordinary `same_batch` particles and
ambient species remain fail-closed on an over-limit ratio. Photoelectrons do not enter the 1-D mean density or outer space
charge; nonlinear photoelectron density and outer-cloud inventory remain unsupported. Escape and
reabsorption charge in the ledger follow analytic closure totals, while integer counts describe terminal ray samples.
Standard output records `transaction_residual_C`, `mean_solver_iterations`, `sample_escape_fraction`, and
`return_weight_scale`.

This composition requires `outer_update_stride=1`, positive `batch_duration`,
`deposit_opposite_charge_on_emit=true`, exactly one negative `photo_raycast` species, and
`photo_raycast.normal_drift_speed=0`; the analytic Maxwellian escape fraction does not include normal emission drift.
This path adds no
mesh-mode-specific requirement.

This quasistationary shadow sampling does not time-evolve delayed return current during UV turn-on. That use requires a
separately designed delayed inventory or queue compatible with `ambient_linear_debye`.
Use this closure only while photoelectron space charge remains small compared with the ambient linear response: it cannot
represent a virtual cathode, a space-charge-limited or inverse sheath, or trapped populations. See
[`periodic2_ambient_linear_photo_outer.toml`](../examples/periodic2_ambient_linear_photo_outer.toml)
for a complete example.

### Connect strong photoemission to a Zhao profile with an empirical CDF

Photoemitting `zhao_charge_driven + same_batch` makes BEACH select a separate `implicit_mean` path automatically. It does
not use the ambient path's analytic-Maxwellian backward Euler update or common `return_weight_scale`. BEACH forms exact
charge-weighted groups of measured normal energy ${\cal E}_n=m_{pe}v_z^2/2$ at z-high and solves the generalized charge root

$$
Q=Q_{\rm base}+C_{\rm esc}[B(Q)]
$$

with the full-profile barrier $B(Q)$, including a Type-A virtual cathode. A fractional escape weight is used only for an
equal-energy group lying on the barrier.

Surface `emit_current_density_a_m2` determines emitted rays and macro weights. The source amplitude passed to Zhao for each
batch is instead resolved from measured charge that actually reaches the interface divided by $A\Delta t$. Reabsorption over
the rough surface and transmission to the interface therefore affect source normalization.

Starting from the committed A/B/C root, the solver follows
$\mathbf G(\mathbf y;E_I(\lambda),n_{pe,0}(\lambda))=\mathbf 0$ with adaptive pseudo-arclength continuation. It fails closed
at a parameter fold; Type A/B requires $E_I>0$, Type C requires $E_I<0$, and a $\lambda=1$ event corrector lands at the
target. For candidate
paths at fixed final source, it checks the tangent barrier slope with respect to prescribed charge at every adaptively
accepted point and stops on a negative slope or oppositely directed secant. Two paths from the common root to $Q_0$ and
$Q_M$ check the complete charge interval before a binary search locates the first true order-predicate index. Pure-root
selection therefore costs $O(\log M)$ connected candidate solves. The measured-CDF charge residual is recomputed at final
acceptance. The local continuation guard is not an interval-arithmetic mathematical proof over the continuum between
accepted points. A branch change, loss of rank, fold, barrier-slope reversal, inconsistent order-predicate bracket, failure
to bracket a marginal energy, or violation of the frozen-cohort 0.25 trust region stops without falling back to another
branch or the old profile. This internal solver adds no TOML key.

Only interface-source amplitude and the normal-energy CDF are measured. Zhao density populations remain analytic; this path
does not solve an arbitrary trapped VDF, collisions, magnetization, or PIC. MPI gathers every interface sample to the root,
so root memory and gather traffic scale with the number of crossing rays. See
[Outer Sheath: Kinetic 1-D](KineticOuterPlasma.en.html#close-strong-photoemission-with-the-measured-interface-energy-cdf)
for the equations, trust region, and ray-terminal contract.

## Choose Interface-Particle Handling

| `particles.mode` | Particle leaving through z-high | Use |
| --- | --- | --- |
| `local_source` | Apply `ordinary_open`; retain no external trajectory | No external field, field-only, or scalar inflow |
| `same_batch` | Evaluate return with the kinetic 1-D profile | Tracked return in steady and quasi-steady studies |
| `zhao_queue` | Hold Zhao photoelectron return/escape events until their due batch | Delayed current during strong-UV turn-on |

`particles.mode` selects only external transport of z-high particles. Select the reservoir-inflow VDF independently with
`inflow_model` below.

Combine `same_batch` with `field.model="kinetic_1d"`. Return and escape are evaluated on the discrete kinetic 1-D profile.
The `implicit_mean` photoelectron path above defers an interface crossing until after the mean-field update, then classifies
that particle through the same profile mapper. Ambient and strong-photoemission Zhao use different mean solvers and ray
weights. Same-batch return for other species is unchanged.

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
mode = "same_batch"
steady_start_mode = "zhao_floating"
steady_start_mesh_id = 1
field_evolution_timescale = 1.0
```

For photoemitting Zhao, the physical length and potential scales come from $T_{pe}$ and the reference density; the no-photo
case derives them from ambient $T_e$ and $n_\infty$. `debye_length` and `thermal_voltage` do not change the Zhao root or profile,
but remain required reference inputs for split-interface gap, lateral-field, and local-charge diagnostics.
Strong-photoemission `same_batch` requires `steady_start_mesh_id` to select a horizontal `mesh.mode="template"` plane that
covers the complete periodic-cell area and can receive the uniform seed charge.
When switching to `mode="zhao_queue"`, remove `steady_start_mode` and `steady_start_mesh_id`, and provide a positive
`sim.batch_duration` and photoelectron source.
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
- Photoemitting Zhao `same_batch` without a `zhao_floating` steady start
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

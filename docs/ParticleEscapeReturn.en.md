title: Particle escape and return

Lang: [日本語](ParticleEscapeReturn.md) | [English](ParticleEscapeReturn.en.md)

# Particle escape and return

Treatment of a particle reaching an open boundary can be organized into five models.

1. `escape`: remove the particle immediately.
2. `potential_barrier`: decide reflection or escape from a scalar potential at the crossing.
3. `linear_debye`: map escape or return through an analytic 1-D Debye profile.
4. `kinetic_1d`: map escape or return through a solved discrete 1-D profile.
5. `unified_linear_response`: integrate an orbit in an external 3-D field to decide escape or return.

This list compares boundary-processing algorithms; it is not an operational ranking. The standard path for a self-consistent
outer sheath is `kinetic_1d` with its matching 1-D transfer. `unified_linear_response` with an explicit 3-D orbit is an advanced
path for cases that require both rough-surface linear screening and 3-D motion outside the box.

These are not five alternatives under one configuration key. `escape` and `potential_barrier` are values of
`sim.open_boundary_model` for ordinary open faces. The other three are values of `outer_plasma.model` when z-high is an
ownership interface to an outer region, combined with corresponding `return_model` and `particle_transfer_mode` settings.

```text
particle crosses an open face
+-- face is not owned by outer transfer
|   +-- escape                    remove unconditionally
|   +-- potential_barrier         reflect or escape at a scalar barrier
+-- z-high is owned by outer transfer
    +-- linear_debye              analytic 1-D return
    +-- kinetic_1d                discrete 1-D profile return
    +-- unified_linear_response   explicit 3-D outer orbit
```

## Compare the five models

| Model | External state | Particle decision | Return time | Typical use |
| --- | --- | --- | --- | --- |
| `escape` | None | Always escape | None | Simple finite box |
| `potential_barrier` | Scalar potential at the crossing | Compare normal energy with a barrier | None | Low-cost local reflection |
| `linear_debye` | Analytic exponential 1-D profile | Conserved energy | Analytic expression | Reduced quasisteady 1-D outer plasma |
| `kinetic_1d` | Converged discrete 1-D profile | Conserved energy and turning-point search | Profile integration | **Standard:** self-consistent mean sheath |
| `unified_linear_response` | 3-D field with zero and nonzero modes | Time-integrate an outer orbit | Measured from the orbit | **Advanced:** linear 3-D response over a rough surface |

All five treatments are independent of particle source. Reservoir particles, photoelectrons, and `volume_seed` particles receive
the same boundary treatment when they cross the same face in the same state. See
[Choosing boundary and outer-domain models](OuterPlasmaModels.en.html) for selecting the external field itself.

## 1. `escape`: remove a particle at an open face

This is the simplest model. Select it for an open face that is not owned by outer transfer.

```toml
[sim]
open_boundary_model = "escape"
```

The particle is removed at the boundary crossing. Its macro charge $qw$ is recorded in species-resolved
`escaped_to_infinity`, while surface charge `q_elem` remains unchanged. There is no external field, turning point, flight time,
or return state.

At a corner crossing multiple faces simultaneously, ordinary event rules still treat the combination deterministically. See
[Particle collision and boundary events](ParticleEvents.en.html) for combinations with reflecting or periodic faces and
reintegration of the step remainder.

## 2. `potential_barrier`: decide reflection at a scalar barrier

This reduced model uses only the scalar potential at the crossing. It selects reflection or escape from a local energy condition
without constructing a spatial profile outside the box.

```toml
[sim]
open_boundary_model = "potential_barrier"
phi_infty = 0.0
```

### Compare only normal energy at the crossing

For crossing-point potential $\phi_b$ and outward normal speed $v_n>0$, the potential barrier to infinity is

$$
U_b=q(\phi_\infty-\phi_b).
$$

If

$$
\frac12m v_n^2<U_b\quad\text{and}\quad U_b>0,
$$

the normal velocity is reversed and the step remainder is tracked. Otherwise the particle escapes. Tangential velocity is
unchanged.

### What this model does not represent

Its only external state is the scalar potential at the crossing. It has no external $E(\mathbf x)$, turning position, flight
time, or space charge. It is not generalized to a corner crossing multiple open faces; that case fails with
`unsupported_barrier_corner`.

The crossing potential follows the same snapshot convention as particle motion and therefore includes the local potential of
`sim.e0`. Because a uniform field has no finite potential at infinity, a configuration combining `sim.e0` with this model must
supply `phi_infty` as a consistent effective reservoir reference.

## 3. `linear_debye`: map return through an analytic 1-D profile

This model represents the potential difference outside the interface by an exponential profile with Debye length $\lambda_D$.
It makes z-high an ownership interface and obtains escape or return, round-trip time, and tangential displacement analytically.

```toml
[outer_plasma]
model = "linear_debye"
return_model = "electrostatic_1d_instant_return"

[coupling]
particle_transfer_mode = "electrostatic_1d_instant_return"
```

### Separate escape and return with conserved energy

For interface potential $\phi_I$, infinity potential $\phi_\infty$, and outward interface normal speed $v_{n,I}$,

$$
v_{n,\infty}^2
=v_{n,I}^2+\frac{2q(\phi_I-\phi_\infty)}{m}.
$$

If $v_{n,\infty}^2\ge0$, the particle can reach infinity and escapes. A negative value implies a turning point in the
exponential profile and therefore return. This is the reverse of the infinity-to-interface map in
[`reservoir_face` inflow and velocity sampling](ReservoirInjection.en.html).

### Derive return time from a linear-Debye profile

For a particle that cannot escape, define

$$
D=-v_{n,\infty}^2>0.
$$

The round-trip time is

$$
\tau_\mathrm{outer}
=\frac{4\lambda_D}{\sqrt D}
\tan^{-1}\left(\frac{v_{n,I}}{\sqrt D}\right).
$$

Return reverses only normal velocity to $-v_{n,I}$ and preserves tangential velocity. Tangential position advances by
$\mathbf v_t\tau_\mathrm{outer}$, then x/y are wrapped into the primary periodic cell.

This model is inexpensive but assumes an exponential external potential. Use `kinetic_1d` when a self-consistent VDF or
nonlinear sheath is required.

## 4. `kinetic_1d`: obtain return from a discrete sheath profile

This model solves a nonlinear 1-D Poisson problem from ambient electron and ion VDFs. It uses the converged $\phi(z)$ at the
start of each batch for both inflow and outflow.

```toml
[outer_plasma]
model = "kinetic_1d"
return_model = "kinetic_1d_profile_return"

[coupling]
particle_transfer_mode = "electrostatic_1d_instant_return"
```

### Find the turning point on a discrete kinetic profile

At every grid point, the model evaluates

$$
v_n^2(z)=v_{n,I}^2+\frac{2q[\phi_I-\phi(z)]}{m}.
$$

The mapper scans outward from the interface. The first interval where $v_n^2$ changes sign brackets a turning point, which is
linearly interpolated for return. A particle escapes only when it can traverse the entire discrete profile and far Robin tail.

The one-way time over a positive-speed interval $a\to b$ is accumulated as

$$
\Delta t=\frac{2\Delta z}{\sqrt{v_{n,a}^2}+\sqrt{v_{n,b}^2}}.
$$

The turning interval is integrated only to the crossing fraction. If the turning point lies above the grid, the far Robin
exponential tail is integrated analytically. Velocity reversal, tangential displacement, and periodic wrapping after the round
trip follow the `linear_debye` treatment.

The current Zhao closure also uses the photoelectron reference Debye length as this off-grid tail scale. It generally differs
from the true asymptotic `abs(phi_end/E_end)`, so near-separatrix flight times are provisional. Quantify this difference before
production use and add an independent tail scale to the outer state when necessary.

### Scan the full barrier of a nonmonotone Type A profile

A Type A profile from `kinetic_closure="zhao_charge_driven"` contains an intermediate potential minimum, so the endpoint
potential difference alone cannot determine accessibility. For an outgoing particle, the mapper scans every profile interval
in z order and uses the first inaccessible interval as the turning point.

For inflow from infinity with normal speed $v_{n,\infty}$, every grid point is checked using

$$
v_n^2(z)=v_{n,\infty}^2+\frac{2q[\phi_\infty-\phi(z)]}{m}.
$$

If this value is negative at any point, the particle cannot reach the interface. This is equivalent to using the most
restrictive potential-energy barrier along the path instead of only the endpoint difference.

### Validate the physical profile branch

The profile is checked for finite values, strictly increasing z coordinates, consistency with the interface point, and the
potential shape required by its closure and resolved branch.

| Closure / branch | Accepted potential shape |
| --- | --- |
| `absorbing_maxwellian` | Monotone increasing or monotone decreasing over the full grid |
| Zhao `A` | Nonincreasing to an interior minimum, then nondecreasing |
| Zhao `B` | Nonincreasing over the full grid |
| Zhao `C` | Nondecreasing over the full grid |
| Zhao `0` | Flat bootstrap |

Other nonmonotone profiles are rejected. An interval that fails to bracket a physical turning point, or a nonpositive Robin tail when one
is required, stops as an invalid model. The instant path uses that tail; the Zhao queue path ends at the finite $L$ boundary
described below. The model can represent a self-consistent mean sheath, but assumes a plane-averaged,
electrostatic, collisionless, unmagnetized 1-D profile.
See [Outer field: kinetic 1D](KineticOuterPlasma.en.html) for details.

## 5. `unified_linear_response`: integrate an external 3-D orbit

This model constructs one linear response field from a rough surface to the far plane, including the zero mode and screened
nonzero modes. `unified_linear_response` alone does not enable particle return. The following 3-D orbit settings explicitly
track particles outside the ownership interface.

```toml
[outer_plasma]
model = "unified_linear_response"
return_model = "electrostatic_3d_explicit_orbit"

[coupling]
particle_transfer_mode = "electrostatic_3d_explicit_orbit"
```

### Advance an outer orbit in the unified 3-D field

The batch-fixed field is integrated with fixed-step velocity Verlet using `outer_orbit_dt`:

$$
\mathbf v^{n+1/2}=\mathbf v^n+\frac{q\mathbf E(\mathbf x^n)}{2m}\Delta t_o,
$$

$$
\mathbf x^{n+1}=\mathbf x^n+\mathbf v^{n+1/2}\Delta t_o,
\qquad
\mathbf v^{n+1}=\mathbf v^{n+1/2}+\frac{q\mathbf E(\mathbf x^{n+1})}{2m}\Delta t_o.
$$

An inward recrossing of the ownership interface returns the particle; an outward crossing of the far plane at the top of the
unified grid escapes. Event position and velocity are linearly interpolated within the final outer step.

### Check orbit convergence and energy

Failure to reach either plane within `outer_orbit_max_steps` requires a persistent queue and stops. The method also checks
relative change in total energy

$$
\mathcal E=\frac12m|\mathbf v|^2+q\phi(\mathbf x)
$$

between the initial and event states against `outer_orbit_energy_tolerance`. Check `outer_orbit_dt` convergence of return and
escape fractions, flight time, and energy error. See
[Outer field: unified linear response](UnifiedLinearResponse.en.html) for field construction and applicability.

## Processing common to outer transfer

### Make z-high the particle-ownership interface

Only an open z-high face with active `coupling.particle_transfer_mode` passes crossing data to an outer model.

| `particle_transfer_mode` | Corresponding model | Particle treatment |
| --- | --- | --- |
| `none` | No outer transfer | Ordinary open boundary |
| `electrostatic_1d_instant_return` | `linear_debye` or `kinetic_1d` | Energy-based escape/return map; the matching Zhao configuration can delay it through the queue |
| `electrostatic_3d_explicit_orbit` | `unified_linear_response` | Time-integrated orbit in a batch-fixed 3-D field |

Both 1-D and 3-D transfer require open z-high and `sim.b0=0`; current models do not include an external magnetic orbit.

### Pass the boundary-crossing state to the outer model

[Particle collision and boundary events](ParticleEvents.en.html) select the first trajectory event between a mesh hit and a box
face. The outer-model crossing record contains:

- position and outward velocity at the interface;
- crossing time within the local Boris step;
- `dt_remaining` after the crossing.

In instant mode, a returned particle is placed just inside the interface, and ordinary Boris/event handling reintegrates only
`dt_remaining`. In queue mode, the resolved return position and velocity are saved in an event record and appended after fresh
sources as a local-domain particle in the batch where it becomes due. Outer flight time is a separate diagnostic from the
local step remainder.

### Keep outer flight outside global simulation time in instant mode

“Instant” in the 1-D model means outer flight affects the state map but does not advance global simulation time. The current 3-D
explicit orbit uses the same convention. A particle returns at the simulation time of the local step in which it left, and
outward and returned charge are recorded in the same batch.

This approximation targets a steady or quasisteady outer plasma. Use the delayed-return queue in the next section for UV
turn-on, abrupt plasma changes, or short-pulse transients.

### Start the instant path from a stationary Zhao profile

`steady_start_mode="zhao_floating"` initializes a fresh run with the Zhao zero-current stationary root and the plane charge
that produces its $E_I$. It does not change the return algorithm. The same kinetic profile connected to
`phi(infinity)=0` supplies both the inflow barrier from the infinity reservoir and instant escape or return outside the
interface. On resume, BEACH restores the checkpoint profile and mesh charge without reseeding the stationary root. Because
this mode bypasses the transient from an uncharged state, do not use it to evaluate delayed return current. See
[Outer field: kinetic 1D](KineticOuterPlasma.en.html#start-stationary-studies-from-the-zhao-zero-current-root) for the
configuration and charge relation.

### Queue outer flight for the transient Zhao closure

For a case that must put outer-flight delay into the batch history, such as strong-UV turn-on, use this Zhao composition.

```toml
[sim]
batch_duration = 2.5e-7

[outer_plasma]
model = "kinetic_1d"
kinetic_closure = "zhao_charge_driven"
photoelectron_histogram_enabled = false
return_model = "kinetic_1d_profile_return"

[coupling]
particle_transfer_mode = "electrostatic_1d_instant_return"
outer_update_stride = 1
field_evolution_timescale = 2.0e-5
max_frozen_field_ratio = 0.2
outer_queue_enabled = true
```

Combinations with `linear_debye`, `absorbing_maxwellian`, the explicit 3-D orbit, or the legacy photoelectron histogram are
rejected. The configuration must also satisfy
`batch_duration <= max_frozen_field_ratio * field_evolution_timescale`.
`outer_queue_enabled=true` connects tracked-particle outer flight and the Zhao photoelectron population through one conserved
inventory. Each rank retains its local events, while the photoelectron macro-particle number is summed over MPI as the Zhao
closure input. After due events are removed at a batch start, the target column is

$$
N_{pe,q}=\frac{1}{A_{xy}}
\sum_{j\in\text{queued photoelectron}}w_j.
$$

Over the finite control volume $L=10\lambda_{D,pe}$, the Zhao solver finds a population scale $\eta$ satisfying

$$
N_{pe,Zhao}(\eta)=
\int_0^L\left[n_{pe,f}(z;\eta)+n_{pe,c}(z;\eta)\right]dz
=N_{pe,q}.
$$

Despite the output name `outer_photoelectron_population_fraction`, $\eta$ is an occupancy scale relative to the stationary
reference population, not a probability. The solver follows a physical path connected to $\eta=0$ over
$0\le\eta\le16$ and permits transient overshoot above one. It does not clamp to `[0,1]`, fall back to a target-independent
full-population solution, or jump to a disconnected branch. Queue mode requires `zhao_branch="auto"`; only continuous A/B or
other branch transitions satisfying the degeneracy condition are allowed. The current bisection additionally requires the
column to increase monotonically with $\eta$ and does not support folds. A target without a connected, monotone path stops with
`no_physical_solution`; a zero target uses exactly $\eta=0$.
$\eta$ scales photoelectron density, infinity quasineutrality, and Sagdeev terms, but not the raw photoelectron emission-current
term in the current diagnostic. That analytic raw current enters the tracked-source consistency check and current-density
diagnostics, but not the root, surface charge, or ledger.

The same $0\le z\le L$ interval is the queue's particle-ownership domain. A turning point before $L$ creates a return event;
reaching $L$ creates an escape event absorbed by the exterior reservoir. Queue mode does not extend a Robin tail beyond $L$ to
classify return.

One batch advances this state in the following order.

1. At the start time $t_b=(b-1)\Delta t_b$ of batch $b$, pop rank-local events that are due.
2. Form $N_{pe,q}$ from the remaining global photoelectron inventory, then refresh the Zhao profile and $\eta$.
3. Generate fresh sources, append due return particles, and run the local particle loop. Count a due escape in this batch's
   `escaped_to_infinity`.
4. For an outward z-high crossing, use the current profile to resolve interface return or reservoir escape at $L$ and
   $\tau_{outer}$. Enqueue it
   at $t_{due}=(b-\tfrac12)\Delta t_b+\tau_{outer}$ using the batch midpoint as its crossing time.
5. Commit surface charge and correct the Zhao profile and $\eta$ with the post-enqueue inventory. This state is the next
   batch's continuation seed and the checkpoint state, so straight and split-resume runs execute the same per-batch sequence.

BEACH particles within one batch do not share a synchronized physical time, so the crossing time is represented by the batch
midpoint. Events are released only at batch starts, quantizing return and escape to `batch_duration`. The terminal state and
due time resolved at enqueue are not reintegrated after the outer field changes. This closure represents flight delay and an
outer photoelectron column; it is not a time-dependent Vlasov--Poisson solve, an outer collision model, or an energy-resolved
cloud evolution.

Halve `batch_duration` and double `batch_count` to compare the same final physical time. At minimum, verify convergence of
$\eta$, the column residual, return/escape current, surface charge, and detachment force. The profile is resampled to a fixed
128-point grid over $0\le z\le10\lambda_{D,pe}$, so a production column-grid study also requires exposing the point count as
an input. Independently refine tracked-particle count, horizontal area, and effective-interface location.

Queue state is written to `outer_event_queue.csv` in serial and one `outer_event_queue_rankNNNNN.csv` per MPI rank. It stores
active phase-space records, terminal outcomes, due times, and `next_event_id`. Queue-file schema 2 stores a
`local_fingerprint` of each rank-local payload, while `outer_queue_fingerprint` in the summary binds the queue contents and
ordering across all ranks. A queue-enabled resume requires every rank file and rejects schema, rank, world-size,
completed-batch, global-event-count, signed-charge, or fingerprint mismatches fail closed.

### Check whether the field can remain frozen during outer flight

In instant mode, snapshot validity over a flight is measured by

$$
\epsilon_\mathrm{ad}
=\frac{\tau_\mathrm{outer}}{\texttt{field\_evolution\_timescale}}
$$

and must not exceed `max_frozen_field_ratio`. If $\tau_\mathrm{outer}/\mathrm{batch\_duration}\gtrsim1$, a converged long-time
mean can still be useful under a strong steady-state assumption, but per-batch return current is not temporally correct.

Zhao queue mode still requires positive `field_evolution_timescale` and `max_frozen_field_ratio`. Because events are released
only at batch starts, it includes the quantization delay $\delta_{poll}$ from $t_{due}$ to the first batch-start poll and the
half-batch bound on uncertainty from approximating the crossing time by the batch midpoint:

$$
\frac{\tau_{outer}+\delta_{poll}+\Delta t_b/2}{\texttt{field\_evolution\_timescale}}
\le\texttt{max\_frozen\_field\_ratio}.
$$

An over-limit event stops before enqueue; it is not discarded and does not fall back to instant return. Configuration
validation additionally bounds one batch interval:

$$
\texttt{batch\_duration}
\le
\texttt{max\_frozen\_field\_ratio}\,
\texttt{field\_evolution\_timescale}.
$$

A violation stops the run. Persistent queuing remains unavailable for the explicit 3-D orbit; an orbit that does not finish within
its in-batch limit stops.

## Verify model results with diagnostics

Species-resolved output separates `interface_outward_gross`, `interface_returned_gross`, and `escaped_to_infinity`. It also
reports maximum `outer_flight_time`, frozen-field ratio, and 3-D orbit energy error.

For Zhao queue mode, inspect `outer_photoelectron_population_fraction`, `outer_photoelectron_column_per_area_m2`,
`outer_photoelectron_column_target_per_area_m2`, `outer_photoelectron_column_residual_per_area_m2`,
`outer_queue_event_count`, and `outer_queue_signed_charge_C` in `summary.txt`.
`charge_ledger_outer_flight_charge_before_C` and `charge_ledger_outer_flight_charge_after_C` include queue stock in the charge
conservation residual.

Gross outward minus returned equals net escape only when transfer coverage and the charge-balance interval match for that
species. See [Inspect Output Files](OutputGuide.en.html) for these fields and the `charge_ledger.csv` format.

## Code reference

- `escape` and `potential_barrier`: [`bem_particle_stepper.f90`](../src/runtime/simulator/bem_particle_stepper.f90)
- `linear_debye` and `kinetic_1d`: [`bem_outer_plasma_interface.f90`](../src/physics/outer_plasma/bem_outer_plasma_interface.f90)
- Delayed event queue: [`bem_outer_event_queue.f90`](../src/runtime/coupling/bem_outer_event_queue.f90)
- Queue checkpoint: [`bem_outer_event_queue_io.f90`](../src/runtime/coupling/bem_outer_event_queue_io.f90)
- `unified_linear_response` 3-D orbit: [`bem_outer_plasma_orbit.f90`](../src/physics/outer_plasma/bem_outer_plasma_orbit.f90)
- Interface transfer and diagnostic aggregation: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- Model-combination validation: [`bem_physics_config_types.f90`](../src/config/bem_physics_config_types.f90)

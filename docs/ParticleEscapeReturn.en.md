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

A particle that traverses the profile escapes. Otherwise, the first interval where $v_n^2$ changes sign brackets a turning
point, which is linearly interpolated for return.

The one-way time over a positive-speed interval $a\to b$ is accumulated as

$$
\Delta t=\frac{2\Delta z}{\sqrt{v_{n,a}^2}+\sqrt{v_{n,b}^2}}.
$$

The turning interval is integrated only to the crossing fraction. If the turning point lies above the grid, the far Robin
exponential tail is integrated analytically. Velocity reversal, tangential displacement, and periodic wrapping after the round
trip follow the `linear_debye` treatment.

### Validate the physical profile branch

The profile is checked for finite values, monotonicity, and consistency with the interface point. A nonmonotone profile, an
interval that fails to bracket a physical turning point, or a nonpositive Robin tail stops as an invalid model. The model can
represent a self-consistent mean sheath, but assumes a plane-averaged, electrostatic, collisionless, unmagnetized 1-D profile.
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
| `electrostatic_1d_instant_return` | `linear_debye` or `kinetic_1d` | Direct energy-based escape/return map |
| `electrostatic_3d_explicit_orbit` | `unified_linear_response` | Time-integrated orbit in a batch-fixed 3-D field |

Both 1-D and 3-D transfer require open z-high and `sim.b0=0`; current models do not include an external magnetic orbit.

### Pass the boundary-crossing state to the outer model

[Particle collision and boundary events](ParticleEvents.en.html) select the first trajectory event between a mesh hit and a box
face. The outer-model crossing record contains:

- position and outward velocity at the interface;
- crossing time within the local Boris step;
- `dt_remaining` after the crossing.

A returned particle is placed just inside the interface, and ordinary Boris/event handling reintegrates only `dt_remaining`.
Outer flight time is a separate diagnostic from this local step remainder.

### Keep outer flight outside global simulation time

“Instant” in the 1-D model means outer flight affects the state map but does not advance global simulation time. The current 3-D
explicit orbit uses the same convention. A particle returns at the simulation time of the local step in which it left, and
outward and returned charge are recorded in the same batch.

This approximation targets a steady or quasisteady outer plasma. UV turn-on, abrupt plasma changes, and short-pulse transients
require a delayed-return queue that retains past outward current.

### Check whether the field can remain frozen during outer flight

Snapshot validity over a flight is measured by

$$
\epsilon_\mathrm{ad}
=\frac{\tau_\mathrm{outer}}{\texttt{field\_evolution\_timescale}}
$$

and must not exceed `max_frozen_field_ratio`. If $\tau_\mathrm{outer}/\mathrm{batch\_duration}\gtrsim1$, a converged long-time
mean can still be useful under a strong steady-state assumption, but per-batch return current is not temporally correct.

A persistent delayed-return queue is not implemented. Configuration validation rejects `outer_queue_enabled=true`, and
trajectories requiring a queue stop.

## Verify model results with diagnostics

Species-resolved output separates `interface_outward_gross`, `interface_returned_gross`, and `escaped_to_infinity`. It also
reports maximum `outer_flight_time`, frozen-field ratio, and 3-D orbit energy error.

Gross outward minus returned equals net escape only when transfer coverage and the charge-balance interval match for that
species. See [Inspect Output Files](OutputGuide.en.html) for these fields and the `charge_ledger.csv` format.

## Code reference

- `escape` and `potential_barrier`: [`bem_particle_stepper.f90`](../src/runtime/simulator/bem_particle_stepper.f90)
- `linear_debye` and `kinetic_1d`: [`bem_outer_plasma_interface.f90`](../src/physics/outer_plasma/bem_outer_plasma_interface.f90)
- `unified_linear_response` 3-D orbit: [`bem_outer_plasma_orbit.f90`](../src/physics/outer_plasma/bem_outer_plasma_orbit.f90)
- Interface transfer and diagnostic aggregation: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- Model-combination validation: [`bem_physics_config_types.f90`](../src/config/bem_physics_config_types.f90)

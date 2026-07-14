title: Particle escape and return

Lang: [日本語](ParticleEscapeReturn.md) | [English](ParticleEscapeReturn.en.md)

# Particle escape and return

A particle reaching an open boundary can be removed, reflected by a scalar barrier, or transferred to an outer-plasma region for
a return-or-escape decision. The choice is independent of particle source. Reservoir particles, photoelectrons, and
`volume_seed` particles receive the same boundary treatment when they cross the same face in the same state.

## Pass the boundary-crossing state to the outer model

[Particle collision and boundary events](ParticleEvents.en.html) select the first trajectory event between a mesh hit and box faces. An ordinary open
face handles escape at the crossing position and velocity. Only the z-high face transfers crossing data to an outer model when
`particle_transfer_mode` is active.

The crossing record contains:

- position and outward velocity at the interface;
- crossing time within the local Boris step;
- `dt_remaining` after the crossing.

When the outer model returns a particle, it is placed just inside the interface and ordinary Boris/event handling reintegrates
only `dt_remaining`. Outer flight time is a separate diagnostic from this local step remainder.

## Remove particles at open faces with no outer state

With `open_boundary_model="escape"` on an open face not owned by outer transfer, the particle is removed immediately. Its macro
charge $q w$ is recorded in species-resolved `escaped_to_infinity`, while surface charge `q_elem` remains unchanged.

At a corner crossing multiple faces simultaneously, ordinary event rules still treat the combination deterministically. See
[Particle collision and boundary events](ParticleEvents.en.html) for combinations with reflect or periodic faces and reintegration of the step remainder.

## Compare only normal energy in a scalar-barrier model

`open_boundary_model="potential_barrier"` is a reduced model using only crossing-point potential $\phi_b$ and
`sim.phi_infty`. For outward normal speed $v_n>0$, the barrier to infinity is

$$
U_b=q(\phi_\infty-\phi_b).
$$

If

$$
\frac12m v_n^2<U_b\quad\text{and}\quad U_b>0,
$$

the normal velocity is reversed and the step remainder is tracked. Otherwise the particle escapes. Tangential velocity is
unchanged.

The external state of this model consists only of the scalar potential at the crossing. Its result is reflect or escape; it has
no external $E(\mathbf x)$, turning position, flight time, or space charge. It is not generalized to a corner crossing multiple open faces; that case fails with
`unsupported_barrier_corner`. The crossing potential follows the same snapshot convention as particle motion and therefore
includes the local potential of `sim.e0`. Because a uniform field has no finite potential at infinity, a configuration combining
`sim.e0` with this model must supply `phi_infty` as a consistent effective reservoir reference.

## Make z-high the particle-ownership interface

`coupling.particle_transfer_mode` connects the z-high open face, used as a particle-ownership interface, to an outer region.

| `particle_transfer_mode` | Corresponding field/return configuration | Particle treatment |
| --- | --- | --- |
| `none` | No outer transfer | Ordinary open boundary |
| `electrostatic_1d_instant_return` | `linear_debye` or `kinetic_1d` | Direct energy-based escape/return map |
| `electrostatic_3d_explicit_orbit` | `unified_linear_response` | Time-integrated orbit in a frozen 3-D field |

Both 1-D and 3-D transfer require open z-high and `sim.b0=0`; current outer-particle models do not include an external magnetic
orbit. Returned x/y positions are wrapped into the primary periodic cell.

## Link inflow and outflow with one energy equation

For interface potential $\phi_I$, infinity potential $\phi_\infty$, and outward interface normal speed $v_{n,I}$,

$$
v_{n,\infty}^2
=v_{n,I}^2+\frac{2q(\phi_I-\phi_\infty)}{m}.
$$

If $v_{n,\infty}^2\ge0$, the particle can reach infinity and escapes. A negative value implies a turning point in the outer
profile and therefore return. This is the reverse of the infinity-to-interface map in
[Reservoir injection](ReservoirInjection.en.html). Using the same $\phi_I-\phi_\infty$ in both directions ties inflow and outflow
to one outer model.

## Derive return time from a linear-Debye profile

With `outer_plasma.model="linear_debye"`, the potential outside the interface is an exponential profile of Debye length
$\lambda_D$. For a particle that cannot escape, define

$$
D=-v_{n,\infty}^2>0.
$$

The implemented round-trip time is

$$
\tau_\mathrm{outer}
=\frac{4\lambda_D}{\sqrt D}
\tan^{-1}\left(\frac{v_{n,I}}{\sqrt D}\right).
$$

Return reverses only normal velocity to $-v_{n,I}$ and preserves tangential velocity. Tangential position advances by
$\mathbf v_t\tau_\mathrm{outer}$, then x/y are wrapped periodically.

This analytic state map derives outer residence time and tangential displacement from conserved energy and constructs the return
state.

## Find the turning point on a discrete kinetic profile

`outer_plasma.model="kinetic_1d"` with `return_model="kinetic_1d_profile_return"` uses the discrete $\phi(z)$ converged at the
start of the batch. At every grid point it evaluates

$$
v_n^2(z)=v_{n,I}^2+\frac{2q[\phi_I-\phi(z)]}{m}
$$

and linearly interpolates a turning point in the first interval where the sign changes. The one-way time over a positive-speed
interval $a\to b$ is accumulated as

$$
\Delta t=\frac{2\Delta z}{\sqrt{v_{n,a}^2}+\sqrt{v_{n,b}^2}}.
$$

The turning interval is integrated only to the crossing fraction. If the turning point lies above the grid, the far Robin
exponential tail is integrated analytically. Return position and velocity after the round trip follow the linear-Debye treatment.

The profile is checked for finite values, monotonicity, and consistency with the interface point. A nonmonotone profile, an
interval that fails to bracket a physical turning point, or a nonpositive Robin tail stops as an invalid model.

## Advance an outer orbit in the unified 3-D field

`unified_linear_response + electrostatic_3d_explicit_orbit` tracks the outer particle in the same electrostatic snapshot that
combines zero and screened nonzero modes. Fixed-step velocity Verlet with `outer_orbit_dt` is

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

Failure to reach either plane within `outer_orbit_max_steps` requires a persistent queue and stops. The method also checks relative
change in total energy

$$
\mathcal E=\frac12m|\mathbf v|^2+q\phi(\mathbf x)
$$

between the initial and event states against `outer_orbit_energy_tolerance`. Check `outer_orbit_dt` convergence of return and
escape fractions, flight time, and energy error.

## Keep outer flight outside global simulation time

“Instant” in the 1-D model means outer flight affects the state map but does not advance global simulation time. The current 3-D
explicit orbit has the same global-time convention. A particle returns at the simulation time of the local step in which it left,
and outward and returned charge are recorded in the same batch.

This approximation targets a steady or quasi-steady outer plasma. It preserves energy-based escape/return and long-time mean
current while eliminating the outer particle stock from the state. UV turn-on, abrupt plasma changes, and short-pulse transients
require a delayed-return queue that retains past outward current.

## Check whether the field can remain frozen during outer flight

Snapshot validity over a flight is measured by

$$
\epsilon_\mathrm{ad}
=\frac{\tau_\mathrm{outer}}{\texttt{field\_evolution\_timescale}}
$$

and must not exceed `max_frozen_field_ratio`. Set the denominator `field_evolution_timescale` to the physical evolution timescale
of surface potential or the outer profile. If $\tau_\mathrm{outer}/\mathrm{batch\_duration}\gtrsim1$, a converged long-time mean can
still be useful under a strong steady-state assumption, but per-batch return current is not temporally correct.

A persistent delayed-return queue for particles violating the gate is not implemented. Configuration validation rejects
`outer_queue_enabled=true`, and trajectories requiring a queue stop.

## Separate interface crossings from final escape

Species-resolved output separates `interface_outward_gross`, `interface_returned_gross`, and `escaped_to_infinity`. It also reports
maximum `outer_flight_time`, frozen-field ratio, and 3-D orbit energy error. Gross outward minus returned equals net escape only
when transfer coverage and the charge-balance interval match for that species. See
[Reading output files](OutputGuide.en.html) for these fields and the `charge_ledger.csv` file format.

## Code reference

- Box events and scalar barrier: [`bem_particle_stepper.f90`](../src/runtime/simulator/bem_particle_stepper.f90)
- Linear and kinetic 1-D return: [`bem_outer_plasma_interface.f90`](../src/physics/outer_plasma/bem_outer_plasma_interface.f90)
- Unified 3-D orbit: [`bem_outer_plasma_orbit.f90`](../src/physics/outer_plasma/bem_outer_plasma_orbit.f90)
- Interface transfer and diagnostic aggregation: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- Return-model combination validation: [`bem_physics_config_types.f90`](../src/config/bem_physics_config_types.f90)

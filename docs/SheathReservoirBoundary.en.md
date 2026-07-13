title: Outer Sheath and Reservoir Particle Boundaries

Lang: [English](SheathReservoirBoundary.en.md) | [日本語](SheathReservoirBoundary.md)

# Outer Sheath and Reservoir Particle Boundaries

This document covers `reservoir_face` inflow, electrostatic acceleration and
deceleration, inaccessible particles, outgoing reflection/return/escape, and
Zhao-family injection corrections. See
[periodic2 Zero Mode and Outer Plasma](PeriodicZeroModeOuterPlasma.en.md) for the
zero mode and outer Poisson solvers.

## 1. Distinguish the models called "sheath"

BEACH contains several reduced models with different responsibilities.

| Feature | Input | What it changes | Supplies spatial $E(z)$ to the particle pusher |
| --- | --- | --- | --- |
| `reservoir_potential_model="infinity_barrier"` | face-average potential and `phi_infty` | inflow flux, normal cutoff, and face speed | no |
| `sheath_injection_model="zhao_*"` | background species, solar elevation, photoelectron reference density | reservoir density/VDF cutoff, ion drift, photoemission current | no |
| `floating_no_photo` | electron/ion inflow fluxes | electron cutoff | no |
| `outer_plasma.model="kinetic_1d"` | infinity VDFs and surface zero-mode field | self-consistent 1D potential, field, and density profile | used in the outer domain |
| `unified_linear_response` | surface field, accessible area, linear plasma response | zero/nonzero field from local to far domain | yes |
| `open_boundary_model="potential_barrier"` | crossing potential and `phi_infty` | outgoing escape or reflection | no |

Zhao is not another solver for the `kinetic_1d` problem. It pre-corrects the
injection distribution and does not solve a batch-dependent outer Poisson problem
from the evolving surface charge.

## 2. Common energy convention

Let $\phi_\infty$ be the infinity potential and $\phi_f$ the reservoir-face or
interface potential. One-dimensional mapping preserves tangential velocity and
normal electrostatic energy:

$$
\frac12m v_{n,f}^2+q\phi_f
=\frac12m v_{n,\infty}^2+q\phi_\infty,
$$

so

$$
v_{n,f}^2=v_{n,\infty}^2-B,\qquad
B=\frac{2q(\phi_f-\phi_\infty)}{m}.
$$

| $B$ | Result |
| ---: | --- |
| $B>0$ | deceleration toward the face; particles with $v_{n,\infty}<\sqrt B$ cannot reach it |
| $B=0$ | unchanged normal speed |
| $B<0$ | acceleration toward the face, with $v_{n,f}>v_{n,\infty}$ |

An inaccessible infinity particle is never instantiated as a simulation
particle. This upstream-VDF accessibility cutoff is different from specularly
reflecting an already tracked particle at a box face.

## 3. Basic `reservoir_face` sampling

For inward face normal $\mathbf n$, `reservoir_face` samples a flux-weighted
normal component from a drifting Maxwellian. With upstream cutoff $v_{\min}$,

$$
\Gamma_{\mathrm{in}}=
\int_{v_n\ge v_{\min}}v_n f(\mathbf v)\,d^3v
$$

sets the batch macro-particle count. With a potential barrier,
$v_{\min}=\sqrt{\max(B,0)}$. The sampler first selects accessible upstream
velocities and then maps their normal speed to the face using the energy equation
above. Tangential components are retained.

Positions are uniform over the face rectangle. After a short velocity jitter,
periodic coordinates are wrapped and nonperiodic coordinates are clamped inside
the box. `infinity_barrier` and outer-profile inflow therefore perform both
upstream flux selection and normal acceleration/deceleration.

## 4. Legacy face-average `infinity_barrier`

### 4.1 Face potential

At the beginning of each batch, the refreshed electrostatic snapshot is sampled
on an `injection_face_phi_grid_n x injection_face_phi_grid_n` cell-centered grid:

$$
\bar\phi_f=\frac{1}{N^2}\sum_{a,b}\phi(\mathbf x_{ab}).
$$

The snapshot uses the configured point/`triangle_p0` kernel, periodic2, physical
`k=0`, outer state, and `e0` conventions. The scalar barrier is
$B=2q(\bar\phi_f-\phi_\infty)/m$.

### 4.2 What it does not model

This model uses only the end-point potential difference. It does not solve an
intermediate $E(z)$, turning location, flight time, or space charge. It also
reduces any lateral face variation to one average. It is a legacy/reduced model
for a small or nearly uniform aperture, not the production default for a rough
surface sheath.

## 5. Inflow through an outer profile

With `particle_transfer_mode="electrostatic_1d_instant_return"`, z-high
`reservoir_face` species represent infinity VDFs.

- `kinetic_1d` obtains $B$ from the refreshed
  `phi_interface-phi_infinity` each batch.
- The linear-Debye reference obtains its interface difference from surface
  zero-mode charge.
- Accessible upstream particles are sampled and energy-mapped to interface speed.

The `kinetic_1d` Poisson solve uses the surface-zero-mode `interface_field` as a
Neumann condition. Inflow speed mapping uses the resulting potential difference,
not $E_I$ directly: field determines the profile, while potential difference
determines particle energy change.

## 6. Outgoing escape and return

### 6.1 Legacy `open_boundary_model="potential_barrier"`

For crossing potential $\phi_b$ and outward normal speed $v_n$, the barrier to
infinity is

$$
U_b=q(\phi_\infty-\phi_b).
$$

If $U_b>\tfrac12mv_n^2$, only the normal velocity is reversed; otherwise the
particle escapes. General simultaneous crossings of multiple open faces are not
implemented and fail closed. Uniform `e0` is excluded because it does not define
a finite infinity reference potential.

### 6.2 Linear-Debye instant return

For an outward interface crossing,

$$
v_{n,\infty}^2=v_{n,I}^2+\frac{2q(\phi_I-\phi_\infty)}{m}.
$$

A nonnegative value means escape. A negative value implies a turning point in the
exponential Debye profile. The analytic round-trip time advances tangential
position by $\mathbf v_t\tau_\mathrm{outer}$, wraps x/y, and reverses normal speed.

### 6.3 Kinetic-profile return

`kinetic_1d_profile_return` uses the converged discrete profile:

$$
v_n^2(z)=v_{n,I}^2+\frac{2q[\phi_I-\phi(z)]}{m}.
$$

The first segment where this changes sign brackets the turning point. The far
Robin exponential tail is integrated analytically after the grid end. A round
trip that is too long relative to `field_evolution_timescale` violates the
frozen-field contract and stops instead of being accepted.

This is not a trajectory integrator that advances through the outer domain with
small time steps. It performs the following reduced mapping:

1. Receive the same-time position and velocity at the interface crossing.
2. If $v_{n,\infty}^2$ is nonnegative, classify the particle as escaping to infinity.
3. Otherwise scan $v_n^2(z)$ over the discrete profile using conserved energy.
4. Interpolate the turning point in the first segment with $v_n^2\le0$.
5. Accumulate one-way time over each piecewise-linear potential segment using
   $$
   \Delta t=\frac{2\Delta z}{\sqrt{v_{n,a}^2}+\sqrt{v_{n,b}^2}},
   $$
   then add the analytic Robin-tail contribution when needed.
6. Construct the interface return state corresponding to round-trip time
   $\tau_\mathrm{outer}$.

On return, only normal velocity is reversed. Tangential velocity is retained,
the tangential position advances by $\mathbf v_t\tau_\mathrm{outer}$ and is
wrapped in x/y, and only the remaining `dt` of the intercepted local step is
reintegrated by the ordinary stepper. Outer flight is not added to global simulation time.

The model is therefore more detailed than an immediate specular reflection
based only on interface speed and a scalar barrier, but it is not explicit
time-stepped outer-orbit tracking. Use
`unified_linear_response + electrostatic_3d_explicit_orbit` for the latter.

#### Why return is immediate and where the approximation applies

In `instant_return`, the outer flight time affects the mapped particle state but
does not advance global simulation time. A turning particle returns at the same
simulation time at which it crossed outward, and outward and returned charge are
recorded in the same batch. The implementation does not wait for either `dt` or
`batch_duration` before returning it.

This is a reduced closure that eliminates a stationary or quasistationary sheath
from the particle domain. In a stationary collisionless electrostatic profile,
total energy determines whether a particle escapes or returns. Once the system
is stationary, mean return current at the surface does not depend on individual
round-trip times. The `kinetic_1d` outgoing/returning density closure includes
the residence-time contribution to outer space charge in the Poisson solve.
Immediate return is therefore the standard choice for stationary potential,
long-time mean current balance, and detachment force after equilibration.

It is not a transient sheath model. After UV turn-on, an abrupt plasma change,
or a short pulse, the physical return current depends on earlier outward current.
The current model does not retain that delay, the net charge temporarily stored
in the outer domain, or delay-driven overshoot and oscillation. Do not use it to
infer transient current or rise time; evaluate quasistationary results after the
initial transient instead.

The quasistatic criterion is

$$
\epsilon_\mathrm{ad}=\frac{\tau_\mathrm{outer}}{\tau_\mathrm{field}},
$$

where `field_evolution_timescale` supplies $\tau_\mathrm{field}$ and
`max_frozen_field_ratio` bounds $\epsilon_\mathrm{ad}$. This compares flight time
with a physical field-evolution time, not with the numerical `dt`. If
$\tau_\mathrm{outer}/\mathrm{batch\_duration}\gtrsim1$, batch-resolved return
current is not temporally faithful. Long-time means after equilibration may
still be useful when $\epsilon_\mathrm{ad}\ll1$, but the batch history must not
be interpreted as a physical transient. Persistent delayed-return queues are
not implemented, and `outer_queue_enabled=true` is rejected.

### 6.4 Unified 3D orbit

`electrostatic_3d_explicit_orbit` advances particles with fixed-step
velocity-Verlet in the combined zero mode and screened nonzero tail. A crossing
back through the ownership plane returns to the local domain; an outward far-plane
crossing escapes. Energy error, flight time, frozen-field ratio, and step count
are checked.

## 7. Zhao-family injection correction

### 7.1 What it solves

Zhao closures combine analytic solar-wind electron, ion, and photoelectron
densities with current constraints. A nonlinear root solve determines surface
potential $\phi_0$, potential minimum $\phi_m$ when needed, and effective
solar-wind electron density $n_{\mathrm{swe},\infty}$. Solar elevation enters as

$$
n_{\mathrm{phe},0}=n_{\mathrm{phe,ref}}\sin\alpha.
$$

| Branch | Implemented potential/population structure |
| --- | --- |
| A | nonmonotone $\phi_0>0$, $\phi_m<0$; captured photoelectrons below the minimum and reflected solar-wind electrons above it |
| B | monotone $\phi_0>0$; photoelectron capture without reflected solar-wind electrons |
| C | monotone $\phi_0<0$; reflected solar-wind electrons and zero photoelectron cutoff |

`zhao_auto` tries C, A, B for `alpha < 20 deg`, otherwise A, B, C. Failure to
find the requested branch stops rather than falling back to another physical model.

### 7.2 Values applied to species

| Species | Override |
| --- | --- |
| first negative `reservoir_face` | effective density and branch-dependent normal cutoff |
| first positive `reservoir_face` | local ion density and cold-beam normal speed when `sheath_reference_coordinate` is set |
| first negative `photo_raycast` | free-photoelectron current density, branch-dependent emission cutoff, and zero normal drift |

When `sheath_reference_coordinate` is set, that plane is sheath coordinate
$z_s=0$. The Zhao 1D profile is sampled at the distance to the shared reservoir
face, reconstructing free/reflected solar-wind electrons, free/captured
photoelectrons, and ion speed. Since this is already a local VDF construction,
the generic barrier energy shift is not applied again.

### 7.3 Important non-guarantees

The reconstructed Zhao $E(z)$ is not added to the particle-pusher field snapshot.
It is not an outer field self-consistent with arbitrary geometry or evolving
`q_elem`. Zhao is an injection/photoemission closure for literature comparison.
Use `kinetic_1d` or `unified_linear_response` when particle orbits must share the
same outer potential profile.

## 8. `floating_no_photo` and reduced photoelectron closure

`floating_no_photo` uses bisection to find a negative potential where cutoff
electron flux equals ion inflow flux, then applies the electron cutoff. It does
not construct a spatial sheath profile.

`photo_escape_model="boltzmann_cutoff"` computes only

$$
f_\mathrm{escape}=\exp\left[-\frac{|q|\max(\phi_\mathrm{emit}-\phi_\infty,0)}{k_BT_\mathrm{PE}}\right].
$$

It is an instantaneous reduced closure with no returned-particle trajectory,
reabsorption position, or delay.

## 9. Compatibility and recommended use

| Goal | Recommended configuration |
| --- | --- |
| infinite-periodic regolith plus self-consistent 1D sheath | `cached_kneq0` + `kinetic_1d` + `kinetic_1d_profile_return` |
| rough-surface linear check without a split window | `unified_linear_response` plus explicit 3D orbit |
| reproduce an older scalar face barrier | `infinity_barrier` |
| compare with the Zhao literature closure | use one `zhao_*` model alone |
| simple no-photoelectron current balance | `floating_no_photo` |

`kinetic_1d_profile_return` rejects `reservoir_potential_model` and Zhao
corrections to avoid applying the same potential difference or cutoff twice.
Current 1D instant-return and explicit-3D-orbit modes require `b0=0`.

## 10. Implementation map

| Operation | Fortran implementation |
| --- | --- |
| reservoir flux and velocity sampling | `src/particles/bem_injection.f90` |
| face-average potential and energy-shift setup | `src/config/bem_app_config_runtime.f90` |
| Zhao runtime coupling | `src/physics/sheath/bem_sheath_runtime.f90` |
| Zhao core/root/profile | `src/physics/sheath/bem_sheath_model_core.f90` |
| 1D outer return | `src/physics/outer_plasma/bem_outer_plasma_interface.f90` |
| unified 3D orbit | `src/physics/outer_plasma/bem_outer_plasma_orbit.f90` |
| legacy open-face reflection | `src/runtime/simulator/bem_particle_stepper.f90` |

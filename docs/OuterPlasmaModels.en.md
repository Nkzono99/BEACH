title: Outer plasma models

Lang: [English](OuterPlasmaModels.en.md) | [日本語](OuterPlasmaModels.md)

# Outer plasma models

An outer-plasma model defines which state closes the region outside the finite particle domain. A closure that modifies only
the inflow VDF, a one-dimensional profile connected at an interface, and a three-dimensional outer field provide different
information to reservoir inflow and outward particles. A simple open boundary carries no outer state and removes particles at
the crossing.

## Classify models by the state carried outside the box

| Model | What it solves | Spatial field for the particle pusher? |
| --- | --- | --- |
| `infinity_barrier` | Inflow cutoff from face-average potential | No |
| `floating_no_photo` | Reduced electron/ion inflow balance | No |
| `zhao_*` | Injection-VDF correction from an analytic sheath closure | No |
| `linear_debye` | Exponential zero-mode response | Used for outer return |
| `kinetic_1d` | A 1D profile from infinity VDFs and Poisson's equation | In the outer region |
| `unified_linear_response` | Linear 1D response including rough surfaces | Combined from local to far field |

`floating_no_photo` and `zhao_*` return source-VDF corrections. See
[Sheath injection closures](SheathInjectionClosures.en.html) for branches, cutoffs, and applicability.

## Split models connect the local field to a 1D profile

`linear_debye` and `kinetic_1d` connect the local mesh domain to a one-dimensional profile beyond an ownership
interface. The surface zero-mode field supplies the interface condition, and the outer potential difference
controls inflow acceleration and outward-particle escape or return.

`kinetic_1d` accepts only solutions satisfying the original Poisson residual, monotonic branch, ion
accessibility, Bohm entry, and infinity quasineutrality. Failure does not silently select another model.

See [Kinetic 1-D outer plasma](KineticOuterPlasma.en.html) for the Poisson problem, VDF density closures, Newton
solve, and physical acceptance conditions.

## The unified model solves one field from surface to far region

`unified_linear_response` places the rough-surface height range, plane-averaged surface source, accessible area,
and linear Debye response on one 1D grid. It is useful when no vacuum split window exists between the surface and
interface. The configuration that also solves species VDFs, a Bohm condition, and current balance is `kinetic_1d`.

See [Unified linear response](UnifiedLinearResponse.en.html) for accessible fraction, the zero-mode Poisson solve,
nonzero-mode reflection and transmission, and the linearity gate.

## Use the same outer state to determine escape and return

| Transfer | Behavior |
| --- | --- |
| Instant return | Derive turning point and flight time from energy and map back at the same simulation time |
| Kinetic-profile return | Integrate turning point and flight time over the converged discrete $\phi(z)$ |
| Explicit 3D orbit | Advance through the combined zero/nonzero outer field |

Instant return is a reduced steady or quasisteady closure. Because outer-flight delay does not advance global
simulation time, it is not suitable for startup, pulses, or transient currents varying on the flight-time scale.

See [Particle escape and return](ParticleEscapeReturn.en.html) for the interface energy criterion, flight time in
linear and kinetic profiles, the 3-D orbit, and frozen-field gate.

## Select a configuration from the required outer physics

- Simple finite box: no outer model and `open_boundary_model="escape"`
- Finite-image scalar reservoir barrier: `infinity_barrier`
- Compare an analytic injection closure: `zhao_*`
- Infinite-periodic surface with self-consistent 1D sheath: `kinetic_1d`
- Linear response including roughness: `unified_linear_response`
- Important lateral variation of outer trajectories: explicit 3D orbit

## Check the acceptance conditions of the selected model

Inspect Poisson residuals, boundary conditions, integrated charge, current balance, grid refinement, and
frozen-field ratio in addition to solver status. See ADR 0001/0002 and
[Validate simulation results](ValidationGuide.en.html).

See [Reservoir injection](ReservoirInjection.en.html) for inflow and
[Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for its source coupling.

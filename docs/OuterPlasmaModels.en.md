title: Outer plasma models

Lang: [English](OuterPlasmaModels.en.md) | [日本語](OuterPlasmaModels.md)

# Outer plasma models

An outer-plasma model defines how the region outside the finite particle domain is closed. Distinguish a simple
open boundary, injection-VDF correction, a one-dimensional Poisson solve, and a three-dimensional outer orbit.

## Model classes

| Model | What it solves | Spatial field for the particle pusher? |
| --- | --- | --- |
| `infinity_barrier` | Inflow cutoff from face-average potential | No |
| `floating_no_photo` | Reduced electron/ion inflow balance | No |
| `zhao_*` | Injection-VDF correction from an analytic sheath closure | No |
| `linear_debye` | Exponential zero-mode response | Used for outer return |
| `kinetic_1d` | A 1D profile from infinity VDFs and Poisson's equation | In the outer region |
| `unified_linear_response` | Linear 1D response including rough surfaces | Combined from local to far field |

## Split outer models

`linear_debye` and `kinetic_1d` connect the local mesh domain to a one-dimensional profile beyond an ownership
interface. The surface zero-mode field supplies the interface condition, and the outer potential difference
controls inflow acceleration and outward-particle escape or return.

`kinetic_1d` accepts only solutions satisfying the original Poisson residual, monotonic branch, ion
accessibility, Bohm entry, and infinity quasineutrality. Failure does not silently select another model.

## Unified linear response

`unified_linear_response` places the rough-surface height range, plane-averaged surface source, accessible area,
and linear Debye response on one 1D grid. It is useful when no vacuum split window exists between the surface and
interface. It is not a nonlinear VDF sheath solver and does not solve a Bohm condition or floating-current balance.

## Particle return

| Transfer | Behavior |
| --- | --- |
| Instant return | Derive turning point and flight time from energy and map back at the same simulation time |
| Kinetic-profile return | Integrate turning point and flight time over the converged discrete $\phi(z)$ |
| Explicit 3D orbit | Advance through the combined zero/nonzero outer field |

Instant return is a reduced steady or quasisteady closure. Because outer-flight delay does not advance global
simulation time, it is not suitable for startup, pulses, or transient currents varying on the flight-time scale.

## Selection guide

- Simple finite box: no outer model and `open_boundary_model="escape"`
- Reproduce a legacy scalar barrier: `infinity_barrier`
- Compare an analytic injection closure: `zhao_*`
- Infinite-periodic surface with self-consistent 1D sheath: `kinetic_1d`
- Linear response including roughness: `unified_linear_response`
- Important lateral variation of outer trajectories: explicit 3D orbit

## Required diagnostics

Inspect Poisson residuals, boundary conditions, integrated charge, current balance, grid refinement, and
frozen-field ratio in addition to solver status. See ADR 0001/0002 and
[Validate simulation results](ValidationGuide.en.html).

Detailed solver contracts remain in the legacy combined
[periodic2 zero mode and outer plasma](PeriodicZeroModeOuterPlasma.en.html#4-coupling-to-outer-plasma-models),
with particle mapping in [Outer sheath and reservoir particle boundaries](SheathReservoirBoundary.en.html).

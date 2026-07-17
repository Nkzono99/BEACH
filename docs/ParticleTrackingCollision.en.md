title: Particle update

Lang: [日本語](ParticleTrackingCollision.md) | [English](ParticleTrackingCollision.en.md)

# Particle update

BEACH stores particle position and velocity at the same time level. Each step first constructs a next-time candidate, then
commits the trajectory only as far as the earliest mesh or box event.

```text
(xⁿ, vⁿ)
    │ sample the field once at a predicted midpoint
    ▼
Boris velocity update + trapezoidal position update
    │
    ▼
candidate (xⁿ⁺¹, vⁿ⁺¹)
    │ compare mesh hit and box-face crossing
    ├─ mesh first ─────── absorb
    ├─ open face first ── escape or pass to the outer interface
    └─ reflect/periodic ─ re-integrate the remaining time
```

| Stage | Details |
| --- | --- |
| Field sample, Boris rotation, and position update | [Boris particle update](BorisPusher.en.html) |
| Triangle collision, box faces, periodic images, and event ordering | [Particle collision and boundary events](ParticleEvents.en.html) |
| Absorption and element charge delta after an event | [Surface charge update](SurfaceModels.en.html) |
| Escape, return, and outer transfer at an open face | [Particle sources and boundaries](ParticleSourcesBoundaries.en.html), [outer-plasma models](OuterPlasmaModels.en.html) |

## State committed by one step

A step ends in an ordinary next-time state, surface absorption, box escape, an outer-interface crossing, or an incomplete status.
Reflect and periodic events keep the particle alive and advance the remaining time with the same update method.

For a mesh hit, the hit position and element index are committed and the candidate endpoint is discarded. The absorbed charge
does not immediately change the field. It is accumulated in a thread-local element delta and committed at the end of the batch.

## Step limit

A particle advances at most `sim.max_step` times while waiting for absorption or escape. A particle still unclassified at the
limit is counted as `survived_max_step`, not silently reassigned to absorption or escape.

`sim.dt` is one particle-step interval, and `sim.max_step * sim.dt` is the maximum tracked time for one particle.
`batch_duration`, in contrast, connects particle supply to surface-charge updates and is not the particle-step interval.

## First checks

- halve `sim.dt` and check stability of trajectories, hit elements, and absorbed counts
- verify that `survived_max_step` is too small to affect the conclusion
- exercise trajectories crossing a periodic seam, corner, or post-reflection remainder
- reconcile absorbed, escaped, and unresolved particle counts with the charge ledger

See [Configuration parameters](Parameters.en.html) for settings and [Inspect Output Files](OutputGuide.en.html) for result categories.

## Code reference

- Step candidates and event ordering: [`bem_particle_stepper.f90`](../src/runtime/simulator/bem_particle_stepper.f90)
- Batch loop tracking particles: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- Step regression tests: [`test_particle_stepper.f90`](../tests/fortran/test_particle_stepper.f90)

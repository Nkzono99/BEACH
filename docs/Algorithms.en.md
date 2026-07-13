title: Computational model overview

Lang: [English](Algorithms.en.md) | [日本語](Algorithms.md)

# Computational model overview

BEACH couples the electric field produced by charge on triangular surfaces with charged-particle motion in
batches. It is not a grid PIC code: boundary-element charges act as the sources for electric-field and
potential evaluation.

This page describes only the overall computation. Equations, discretizations, and implementation details are
kept in the linked pages.

## Persistent state

The main state carried between batches is:

| State | Contents |
| --- | --- |
| Surface mesh | Triangle vertices, centroids, normals, areas, and surface models |
| Element charge | Total accumulated charge `q_elem` on each triangle |
| Injection state | Values carried to the next batch, such as fractional reservoir macro-particle counts |
| Statistics | Absorption, escape, unresolved particles, and completed batches |
| Restart data | RNG state and model/mesh/species fingerprints |

Particles are normally created and tracked within one batch. Charge from particles absorbed by the surface is
added to the element charge and contributes to the field in the next batch.

## One-batch flow

```text
Current surface charge
        ↓
Evaluate field and potential
        ↓
Generate particles
        ↓
Advance particles ── select the first box or triangle event
        ↓
Accumulate charge deltas from absorbed particles
        ↓
Commit element charge
        ↓
Update statistics, history, and restart state
```

Charge deltas are committed at the end of the batch rather than immediately for each particle. Particles in
one batch therefore share a field produced by the surface charge at the start of that batch. See the
[batch coupling algorithm](BatchAlgorithm.en.html) for the detailed ordering.

## Physical models

### Surface interaction

The primary model absorbs particles at the surface and accumulates their charge on the hit insulating element.
Floating-conductor equipotential relaxation is available, while `epsilon_r` on a `dielectric` is currently
metadata rather than an independent polarization boundary condition.

See [Surface charging models](SurfaceModels.en.html) for post-collision deposition and surface models.

### Particle sources and outer boundaries

Particles can be generated from volume seeds, reservoir faces, and photo ray casts. A particle moving through
an open boundary can escape, reflect, or return according to the selected boundary model.

See [Particle sources and boundaries](ParticleSourcesBoundaries.en.html) for generation and box boundaries, and
[Outer plasma models](OuterPlasmaModels.en.html) for sheath and return closures.

### Periodic domains and outer plasma

`periodic2` is a two-axis-periodic field boundary. It distinguishes near images, the infinite-periodic nonzero
mode, the surface-average zero mode, and an optional outer-plasma response.

See [periodic2 field evaluation](PeriodicElectrostatics.en.html) and
[Outer plasma models](OuterPlasmaModels.en.html).

## Numerical methods

### Particle tracking

Velocity is updated with the Boris method and position with a trapezoidal update using same-time states. Each
step checks the motion segment against box boundaries and mesh triangles and accepts only the earliest event.

### Field and potential

Triangle charge is evaluated either as a point charge at the centroid or as a constant density over the
triangle. Direct, treecode, and FMM evaluation are available according to problem size and boundary conditions.
See [Field solvers and boundary conditions](FieldSolvers.en.html).

### Time scales

`sim.dt` is the particle integration step. `batch_duration` instead connects particle supply with surface-charge
updates. See [`batch_duration` stability and steady value](BatchDurationStability.en.html) for stability and
convergence checks.

## Documentation map

| Question | Page |
| --- | --- |
| What is the update order within one batch? | [Batch coupling algorithm](BatchAlgorithm.en.html) |
| How is absorbed charge stored on surfaces? | [Surface charging models](SurfaceModels.en.html) |
| How are particles generated and handled at box boundaries? | [Particle sources and boundaries](ParticleSourcesBoundaries.en.html) |
| How are particles advanced and collided with triangles? | [Particle tracking and collision](ParticleTrackingCollision.en.html) |
| How should direct, treecode, or FMM be selected? | [Field solvers and boundary conditions](FieldSolvers.en.html) |
| How are periodic images and the zero mode combined? | [periodic2 field evaluation](PeriodicElectrostatics.en.html) |
| How do sheath, outer plasma, escape, and return work? | [Outer plasma models](OuterPlasmaModels.en.html) |
| How should FMM be selected and verified? | [Field evaluation with FMM](FMM.en.html) |
| How is the FMM core implemented? | [Coulomb FMM core details](FMMCore.en.html) |
| How should discretization convergence be checked? | [Validate simulation results](ValidationGuide.en.html) |

Use the [Configuration parameters](Parameters.en.html) when looking for input keys.

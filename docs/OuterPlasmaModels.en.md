title: Choosing Boundary and Outer-Domain Models

Lang: [English](OuterPlasmaModels.en.md) | [日本語](OuterPlasmaModels.md)

# Choosing Boundary and Outer-Domain Models

Connecting a finite particle box to an external reservoir involves five stages: particle creation, inflow correction, outflow
boundary handling, outer-field construction, and particle motion outside the box. These are not alternatives under one `model`
key; they are separate stages that must be combined consistently.

## Distinguish five configuration stages

| Stage | Main configuration | Role |
| --- | --- | --- |
| Particle creation | `particles.species[].source_mode` | Create macro particles from `reservoir_face`, `photo_raycast`, or another source |
| Inflow correction | `sim.reservoir_potential_model` or `sim.sheath_injection_model` | Set accessibility, density, drift, and cutoff of the upstream VDF |
| Outflow boundary | `sim.open_boundary_model` | Select unconditional escape or reflection by a scalar barrier |
| Outer field | `outer_plasma.model` | Construct a 1-D exterior profile or a response field from the surface to the far region |
| Outside particle | `coupling.particle_transfer_mode` | Map an interface crossing through 1-D return or a 3-D outer orbit |

`reservoir_face` is a particle source. `infinity_barrier` and sheath-injection closures modify the VDF supplied to that source.
`potential_barrier` acts in the opposite direction on particles leaving through an open face. `kinetic_1d` and
`unified_linear_response` retain more spatial information than a scalar correction.

## Separate particle creation from inflow-VDF correction

```text
upstream VDF
   |
   +-- no correction
   +-- infinity_barrier         energy map from mean face potential
   +-- sheath_injection_model  analytic density, drift, and cutoff closure
   +-- kinetic_1d profile      energy map through solved outer potential
             |
       reservoir_face
       integrate flux and create particles
```

`source_mode="reservoir_face"` converts an upstream VDF into a flux-weighted distribution and creates particles on a box face.
The sampler does not solve a sheath by itself. The diagram shows the reservoir-inflow path. Zhao-family closures additionally
correct the density, cutoff, and emitted current of the matching `photo_raycast` source.

| Inflow model | Resolved source information | Spatial profile |
| --- | --- | --- |
| No correction | Treat configured VDF as the face distribution | None |
| `infinity_barrier` | Accessibility cutoff and face speed from mean face potential and `phi_infty` | None |
| `floating_no_photo` | Electron cutoff balancing electron and ion inflow | None |
| `zhao_*` | Electron, ion, and photoelectron VDF corrections from an analytic branch | Analytic closure only |
| `kinetic_1d` | Map inflow to the interface through a converged outer Poisson profile | Discrete 1-D profile |

`floating_no_photo` and `zhao_*` are selected with `sim.sheath_injection_model`. They preprocess source VDFs; they do not provide
the spatial field used to push generated particles. See [Reservoir Injection](ReservoirInjection.en.html) for sampling and
[Inflow-VDF Sheath Closures](SheathInjectionClosures.en.html) for analytic closures.

## Separate open-boundary handling from outer transfer

A particle reaching an open face follows either an ordinary open-boundary rule or transfer into an outer region.

| Outflow model | Action after crossing | State outside the box |
| --- | --- | --- |
| `open_boundary_model="escape"` | Remove the particle immediately | None |
| `open_boundary_model="potential_barrier"` | Reflect or escape from crossing potential and normal kinetic energy | Scalar potential only |
| `electrostatic_1d_instant_return` | Compute escape, turning point, and round-trip time from a 1-D profile | Analytic or discrete 1-D profile |
| `electrostatic_3d_explicit_orbit` | Integrate an orbit in the batch-fixed external 3-D field | Zero/nonzero 3-D field |

`potential_barrier` is independent of particle source. Reservoir particles, photoelectrons, and volume-seeded particles receive
the same decision when they cross the same face in the same state. When z-high is an outer ownership interface, escape and return
on that face are owned by the corresponding outer model rather than `open_boundary_model`.

See [Open Boundaries, Escape, and Return](ParticleEscapeReturn.en.html) for decision equations, return time, 3-D orbits, and the
quasi-steady approximation.

## Choose an outer model by required field complexity

| `outer_plasma.model` | Solved state | Appropriate use |
| --- | --- | --- |
| `none` | No outer field | Simple open boundaries, scalar barriers, and analytic injection closures |
| `linear_debye` | Exponential zero-mode response from an interface | Reduced 1-D instant return |
| `kinetic_1d` | Nonlinear 1-D Poisson profile with VDF closures | Self-consistent inflow, current, and escape/return |
| `unified_linear_response` | Linear zero/nonzero response from a rough surface to the far region | Rough surfaces without a clean vacuum split window |

`kinetic_1d` solves a monotone branch satisfying ambient electron/ion far VDFs, the Bohm condition, Poisson residual, and far
quasineutrality. `unified_linear_response` does not solve species VDFs or current balance; it adds a linear Debye-Hueckel response
and plasma-accessible area to the field.

Use [Outer Field: kinetic 1D](KineticOuterPlasma.en.html) when a self-consistent mean sheath is required. Use
[Outer Field: unified linear response](UnifiedLinearResponse.en.html) when roughness and plasma response occupy the same region
but a linear approximation is acceptable.

## Select a typical composition

| Goal | Inflow | Outflow | Outer field and transfer |
| --- | --- | --- | --- |
| Simple finite box | `reservoir_face`, no correction | `escape` | None |
| Finite images with scalar barriers | `infinity_barrier` | `potential_barrier` | None |
| Zhao literature closure | `sheath_injection_model="zhao_*"` | Ordinary open boundary | None |
| Self-consistent 1-D sheath | Inflow mapped through kinetic profile | Kinetic profile return | `kinetic_1d` with 1-D transfer |
| Linear response over a rough surface | Source prescribed at the face | Open or 3-D outer orbit | `unified_linear_response` |

Complete examples are provided in [Finite Periodic Configuration](FinitePeriodicConfiguration.en.html) and
[Infinite-Periodic Configuration with Outer Plasma](InfinitePeriodicOuterConfiguration.en.html).

## Do not stack corrections with the same responsibility

- Do not combine `sheath_injection_model` with `reservoir_potential_model`.
- Kinetic profile return uses the same profile for inflow; do not add Zhao or `infinity_barrier`.
- `unified_linear_response` does not determine source VDFs or floating-current balance. Define any reservoir distribution at the face.
- Do not interpret z-high outer transfer and scalar reflection by `potential_barrier` as the same operation.
- A failed model does not silently fall back to a simpler closure.

Inspect the Poisson residual, boundary conditions, integrated charge, current balance, grid refinement, and frozen-field ratio as
applicable to the selected model. See [Validating Simulation Results](ValidationGuide.en.html) for the common workflow.

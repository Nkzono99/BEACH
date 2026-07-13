title: Finite-image periodic2 configuration

Lang: [日本語](FinitePeriodicConfiguration.md) | [English](FinitePeriodicConfiguration.en.md)

# Finite-image periodic2 configuration

This configuration defines the surface field of an x/y-periodic box as the sum of the primary cell and a selected finite image
shell. It carries no spatial outer-plasma profile; inflow, outflow, and photoelectrons can use reduced closures based on scalar
potential.

## Keep field, inflow, and outflow within one finite-image model

| Process | Choice |
| --- | --- |
| Surface field | `field_bc_mode="periodic2"`, finite image sum, no far correction |
| Source | `volume_seed`, `reservoir_face`, or `photo_raycast` |
| Reservoir correction | None or legacy `infinity_barrier` |
| Open outflow | Unconditional `escape` or legacy `potential_barrier` |
| Photoelectron | Ordinary tracking or reduced `boltzmann_cutoff` escape |
| Outer Poisson/profile | None |

See [periodic2 electrostatics](PeriodicElectrostatics.en.html) for the finite-image kernel, difference from Ewald far correction,
and collision-image search.

## Let the image shell define the physical field range

For image layer $N$, sources from $(2N+1)^2$ cells are added. This can approximate an infinite-periodic sum, but the result at a
configured $N$ is a finite-image model. For a nonneutral cell, increasing image layer does not automatically choose potential gauge
or the far boundary along z.

`field_periodic_far_correction="none"` and compatibility `auto` add no Ewald or cached far operator. Use this path only when field
observables, trajectories, absorption locations, and charging distributions converge under image-layer refinement.

## Correct reservoir inflow with face-average potential

Without correction, the configured upstream VDF is sampled as a flux-weighted distribution on the face with no potential cutoff or
acceleration.

With `reservoir_potential_model="infinity_barrier"`, the batch-start finite-image snapshot gives mean aperture potential
$\bar\phi_f`. Relative to `phi_infty`, it:

- selects accessible normal speeds from the upstream VDF;
- maps accepted normal velocity to the face by the same potential difference.

See [Reservoir injection](ReservoirInjection.en.html) for equations and `injection_face_phi_grid_n`.

This uses one face-average scalar and has no intermediate $E(z)$, turning position, flight time, or space charge. Because
$\bar\phi_f$ can change with image layer, converge face potential as well as particle flux.

## Select escape or scalar reflection at open faces

`open_boundary_model="escape"` removes a particle crossing an open face. `potential_barrier` compares crossing-point potential to
`phi_infty` and reflects only particles whose normal energy is insufficient.

Reservoir `infinity_barrier` and outflow `potential_barrier` can use the same energy convention, but the former filters an upstream
VDF with a face average while the latter reflects an already created particle at its individual crossing. They are not forward and
reverse maps through one spatial sheath profile. See [Particle escape and return](ParticleEscapeReturn.en.html), including corner
limitations.

## Use in-box tracking or a reduced cutoff for photoelectrons

With `photo_escape_model="none"`, a photoelectron created at a ray hit is tracked as an ordinary particle until reabsorption or
open escape. Return after leaving the box is not represented.

`boltzmann_cutoff` derives escape fraction from local source potential and reduces tracked macro weight. It is a reduced closure
that never creates the nonescaping trajectory or reabsorption location. Source reaction charge uses the same reduced weight.

See [Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for ray visibility, surface charge, and tracked versus
reduced escape.

## Configure the finite-image model

```toml
[sim]
field_bc_mode = "periodic2"
field_solver = "fmm"
field_periodic_image_layers = 1
field_periodic_far_correction = "none"
use_box = true
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
open_boundary_model = "escape"
reservoir_potential_model = "none"
```

For `infinity_barrier`, specify `phi_infty` and `injection_face_phi_grid_n`. A reduced photoelectron cutoff needs the same
`phi_infty` gauge. See [Configuration parameters](Parameters.en.html).

The runnable [`examples/periodic2_basic/beach.toml`](../examples/periodic2_basic/beach.toml) is the minimal numerical-path check.
Choose image layers and sampling counts for a production case from the convergence checks below.

## Choose the finite-image model for its valid range

Suitable for:

- small comparisons with explicit image layer;
- regression tests of free/periodic boundaries;
- reproduction of legacy scalar-barrier results;
- a near-image reference for Ewald/cached configurations.

Not suitable for:

- claiming an infinite-periodic surface without layer refinement;
- cases requiring self-consistent outer space charge or a sheath potential profile;
- transients where outer flight time, delayed return, or a 3-D external orbit matters.

## Converge from image shells through charging distributions

1. Increase `field_periodic_image_layers` through $N,N+1,N+2$.
2. Compare reservoir flux, absorption/escape rates, and charging distribution in addition to point $\phi,\mathbf E$.
3. Increase `injection_face_phi_grid_n` when using `infinity_barrier`.
4. Independently converge `rays_per_batch` for photo-raycast and `dt` for particle tracking.
5. Inspect charge ledger and `discarded_unresolved`.

If convergence is slow or a nonneutral cell requires a physical far closure, move to
[Infinite-periodic periodic2 with outer plasma](InfinitePeriodicOuterConfiguration.en.html).

title: Finite-image periodic2 configuration

Lang: [日本語](FinitePeriodicConfiguration.md) | [English](FinitePeriodicConfiguration.en.md)

# Finite-image periodic2 configuration

This configuration defines the surface field of an x/y-periodic box as the sum of the primary cell and a selected finite image
shell. It carries no spatial outer-plasma profile; inflow, outflow, and photoelectrons can use reduced closures based on scalar
potential.

## Typical configuration used on this page

| Process | Choice |
| --- | --- |
| Surface field | `field_bc_mode="periodic2"`, finite image sum, no far correction |
| Source | `reservoir_face` and `photo_raycast`, plus `volume_seed` when an initial population is needed |
| Reservoir correction | Map the upstream VDF to the face with `infinity_barrier` |
| Open outflow | Classify reflection or escape at each crossing with `potential_barrier` |
| Photoelectron | Track as an ordinary particle inside the box |
| Outer Poisson/profile | None |

See [periodic2 electrostatics](PeriodicElectrostatics.en.html) for the finite-image kernel, difference from Ewald far correction,
and collision-image search.

## Parameters that enable this configuration

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
open_boundary_model = "potential_barrier"
reservoir_potential_model = "infinity_barrier"
phi_infty = 0.0
injection_face_phi_grid_n = 5
```

| Parameter | Typical value | Meaning in this configuration |
| --- | --- | --- |
| `field_bc_mode` | `"periodic2"` | Add x/y periodic images to the field source |
| `field_periodic_image_layers` | Start at `1`, then converge | Set the finite image-shell range |
| `field_periodic_far_correction` | `"none"` | Keep a finite-image model without an Ewald or cached far operator |
| `reservoir_potential_model` | `"infinity_barrier"` | Use face-average potential for upstream accessibility and face velocity |
| `open_boundary_model` | `"potential_barrier"` | Use energy at each crossing to classify reflection or escape |
| `phi_infty` | Match the physical gauge | Shared reference potential for reservoir and open outflow |
| `injection_face_phi_grid_n` | Start at `5`, then converge | Samples per periodic axis for face mean and variation |

The reservoir inflow and open outflow therefore share `phi_infty` as their reference potential.

Use the following combination for a photoelectron source:

```toml
[[particles.species]]
source_mode = "photo_raycast"
emit_current_density_a_m2 = 2.0e-4
rays_per_batch = 500
deposit_opposite_charge_on_emit = true
q_particle = -1.602176634e-19
m_particle = 9.1093837139e-31
temperature_ev = 1.5
normal_drift_speed = 1.0e5
```

The emitted photoelectron is tracked as an ordinary particle inside the box. Surface return is absorbed, while an open-face
crossing is classified by `potential_barrier`.
`deposit_opposite_charge_on_emit=true` adds the opposite reaction charge to the emitting element. See
[Configuration parameters](Parameters.en.html) for individual constraints.

The runnable [`examples/periodic2_basic/beach.toml`](../examples/periodic2_basic/beach.toml) checks only the minimal field path.
Add the reservoir and photoelectron fragments above for this typical configuration, then choose image layers and sampling counts
from the convergence checks below.

## Choose how many neighboring periodic cells contribute to the field

An x/y-periodic system repeats the primary cell to its left, right, front, and back. This configuration does not sum every copy in
the infinite lattice. It includes only copies through the selected image layer.

| `field_periodic_image_layers` | Cells included in the field |
| --- | --- |
| `0` | Primary cell only ($1\times1$) |
| `1` | One surrounding layer ($3\times3=9$ cells) |
| `2` | Two surrounding layers ($5\times5=25$ cells) |

Image layer $N$ therefore does not refine the mesh. It sets how many cells away charge remains part of the interaction. Periodic
copies outside that range do not contribute to this finite-image field.

Increase the layer until the quantities of interest—field, particle absorption and escape rates, and final charging
distribution—change negligibly. Continued changes mean that the current layer cuts off too much of the distant-cell influence.

`field_periodic_far_correction="none"` and compatibility alias `"auto"` do not replace omitted cells with an Ewald sum or cached
operator. For a cell with nonzero net charge, increasing the layer alone also does not define the infinite-periodic potential gauge
or the far boundary along z. Use the infinite-periodic plus outer-plasma configuration when that physical far closure is required.

## Correct reservoir inflow with face-average potential

With `reservoir_potential_model="infinity_barrier"`, the batch-start finite-image snapshot gives mean aperture potential
$\bar\phi_f$. Relative to `phi_infty`, it:

- selects accessible normal speeds from the upstream VDF;
- maps accepted normal velocity to the face by the same potential difference.

See [`reservoir_face` inflow and velocity sampling](ReservoirInjection.en.html) for equations and `injection_face_phi_grid_n`.

This uses one face-average scalar and has no intermediate $E(z)$, turning position, flight time, or space charge. Because
$\bar\phi_f$ can change with image layer, converge face potential as well as particle flux.

The same `N x N` samples used by the mean also accumulate the potential population standard deviation, minimum, and maximum, so
this diagnostic adds no potential evaluations. For a Maxwellian reservoir, the MPI root warns on the first and final batch when
the energy represented by the local-potential standard deviation exceeds 10% of the characteristic thermal plus normal-drift
energy.

## Classify crossing-point energy at open faces

`open_boundary_model="potential_barrier"` compares crossing-point potential to `phi_infty`, reflects particles whose normal
energy is insufficient, and removes particles that can cross the barrier as escape.

Reservoir `infinity_barrier` and outflow `potential_barrier` can use the same energy convention, but the former filters an upstream
VDF with a face average while the latter reflects an already created particle at its individual crossing. They are not forward and
reverse maps through one spatial sheath profile. See [Particle escape and return](ParticleEscapeReturn.en.html), including corner
limitations.

## Track photoelectrons normally inside the box

A photoelectron created at a ray hit is tracked as an ordinary particle until reabsorption or an open-face crossing. The same
`potential_barrier` used for other particles classifies a photoelectron reaching an open face. Return after escape from the box
is not represented.

See [Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for ray visibility, surface charge, and tracked escape.

## Choose the finite-image model for its valid range

Suitable for:

- small comparisons with explicit image layer;
- regression tests of free/periodic boundaries;
- finite-image reservoir comparisons using scalar barriers;
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

title: Photoelectron emission and lifecycle

Lang: [日本語](PhotoelectronEmission.md) | [English](PhotoelectronEmission.en.md)

# Photoelectron emission and lifecycle

`source_mode="photo_raycast"` emits particles from the first surface hit by an illumination ray. Mesh ray casting determines
the emission location and surface state stores source charge, but the photoelectron becomes an ordinary particle after creation.
It uses the same field snapshot, Boris update, collision query, and box boundaries.

## Track emission through reabsorption in the same batch

1. Launch rays from an illumination aperture on a box face.
2. Find the first mesh hit while applying box boundary conditions.
3. Create a photoelectron from the plasma-facing normal of the hit element.
4. Optionally record $-qw$ on the emission element.
5. Track it as an ordinary particle and hand it to the common escape/return treatment after it reaches a box boundary.
6. Commit emission and absorption charge to the surface at batch end.

Emission and reabsorption can occur in the same batch, but the field snapshot is not refreshed between them. Net surface charge
created by emission therefore changes the field starting in the next batch.

## Locate the emitting surface with illumination rays

`inject_face` and `pos_low`/`pos_high` define a rectangular aperture on a box face. `ray_direction` must point inward from the
aperture and defaults to the inward face normal. For aperture area $A$, inward normal $\mathbf n_\mathrm{in}$, and normalized ray
direction $\hat{\mathbf d}$, the area projected perpendicular to the rays is

$$
A_\mathrm{proj}=A\left|\hat{\mathbf d}\cdot\mathbf n_\mathrm{in}\right|.
$$

Ray origins are uniform in the aperture. Each ray searches the first triangle hit on successive segments from its current
position to the next box face.

| Reached box face | Ray treatment |
| --- | --- |
| `open` | Leave the box and terminate without emission |
| `reflect` | Reverse the normal direction component and continue |
| `periodic` | Wrap to the opposite face and continue |

With `field_bc_mode="periodic2"`, first-hit search includes periodic images and wraps the emission position into the primary
cell. The hit element itself must lie inside the simulation box. A ray that still has no hit after `raycast_max_bounce` creates
no particle. An incomplete collision-query status such as an image limit, index range, or stalled DDA fails the batch instead of
silently treating the ray as a miss.

## Assign emitted current to each ray

For emission current density $J_\mathrm{emit}>0$, physical-particle charge $q\ne0$, and global ray count $N_\mathrm{ray}$,
the macro-particle weight created by a hit ray is

$$
w_\mathrm{hit}
=\frac{J_\mathrm{emit}A_\mathrm{proj}\,\Delta t_\mathrm{batch}}
{|q|N_\mathrm{ray}}.
$$

A missed ray creates no particle, so shadowing and apparent surface area enter total emission through hit probability.
`w_particle` and `target_macro_particles_per_batch` are not accepted for `photo_raycast`.

`rays_per_batch` is the sampling count for the illumination integral, not the physical emission amount. Increasing it reduces
$w_\mathrm{hit}$ and Monte Carlo noise in visibility and emission location. Results must be checked for convergence with ray
count.

## Build the emitted state from the surface normal

If the stored triangle normal points along the incident ray, it is reversed to form emission normal $\mathbf n_s$ facing back
toward the illumination side. Position is offset $10^{-12}$ m along $\mathbf n_s$ to avoid immediate self-intersection with the
source element.

With $\sigma=\sqrt{k_\mathrm{B}T/m}$, velocity is sampled in local basis
$(\mathbf n_s,\mathbf t_1,\mathbf t_2)$:

- normal speed follows a flux-weighted half-range Maxwell distribution with `normal_drift_speed`;
- both tangential components are zero-mean Gaussians with standard deviation $\sigma$;
- when Zhao or another closure supplies $v_{\min}$, normal speed satisfies $v_n\ge v_{\min}$;
- Gaussian sampling is truncated at $6\sigma$.

Normal velocity is positive, so a new particle leaves the surface toward the illumination side. Its tracked orbit and common
box-boundary or outer-transfer treatment then decide whether it returns or escapes.

## Check charge balance across emission, reabsorption, and escape

With `deposit_opposite_charge_on_emit=true`, source element $i$ receives

$$
\Delta q_{i,\mathrm{emit}}=-q w.
$$

For an electron $q<0$, this leaves positive charge. The difference is accumulated separately as `photo_emission_dq`, combined
with collision deposition after MPI all-reduce, and included in the same batch commit.

If the emitted particle is reabsorbed on element $j$, ordinary absorption deposits $+qw$ there. Return to the source element
cancels emission charge; return to another element transfers net charge across the surface. The current insulator model applies
no subsequent lateral surface conduction.

## Use the common escape and return treatment after emission

A ray hit always creates a photoelectron with weight $w_\mathrm{hit}$. There is no emission-time setting that multiplies this
weight by an escape fraction. Surface return is handled as an ordinary collision, while a particle reaching an open face uses
the same `external_boundary.ordinary_open` or outer particle mode as reservoir and `volume_seed` particles.

For a finite box without a solved external region,
`external_boundary.ordinary_open.model="potential_barrier"` classifies reflection or escape
from crossing-point potential and normal kinetic energy. With a self-consistent external sheath,
`external_boundary.field.model` and `external_boundary.particles.mode` are canonical. See
[Particle escape and return](ParticleEscapeReturn.en.html) for the scalar barrier, 1-D outer-profile return, and unified 3-D
explicit orbit.

Every `photo_raycast` species using tracked outer transfer requires `deposit_opposite_charge_on_emit=true`, regardless of whether
the histogram is enabled, to close emission and return charge balance.

## Include photoelectrons in the mean outer-plasma density

`external_boundary.field.photoelectron_density_model="kinetic_mean"` uses the temperature and emission current density of the first negative
`photo_raycast` species as a plane-averaged source in the 1-D Poisson density closure. Its outgoing and turning-return populations
contribute to outer space charge, but it neither replaces surface absorption of tracked particles nor deposits an extra
statistical return current on the surface.

The mean density model and histogram have separate responsibilities, but their
currently supported model/mode combinations differ, so they cannot be enabled together.

Tracked photoelectrons transferred through z-high use the same source-independent escape/return treatment as other particles.
See [Particle escape and return](ParticleEscapeReturn.en.html) for the quasi-steady approximation that omits outer flight from
global time and the 3-D explicit orbit, and [Outer plasma models](OuterPlasmaModels.en.html) for field construction.

Legacy Zhao, selected with
`external_boundary.particles.inflow_model="legacy_sheath"`, applies only
branch-dependent emission current density, normal cutoff, and drift to source
sampling. It does not construct a spatial $E(z)$, so generated particles
advance in the separately composed batch-fixed field.

In contrast,
`external_boundary.field.kinetic_closure="zhao_charge_driven"` constructs a
self-consistent 1-D outer profile that preserves the interface field set by
accumulated charge. `external_boundary.particles.mode="zhao_queue"`
additionally closes the photoelectron population from tracked outer inventory.
With `mode="same_batch"` or `"zhao_queue"`, that profile controls inflow and
return or escape at the z-high interface. This path is separate from the legacy
Zhao source correction.

## Save a photoelectron histogram at the outer interface

`external_boundary.field.photoelectron_histogram_enabled=true` bins outward `photo_raycast` crossings at z-high by normal kinetic energy.
It writes previous-batch and cumulative signed charge, total kinetic energy, tangential momentum, and count to
`photoelectron_histogram.csv`. The run stops if signed charge crossing the z-high interface outward, relative to the configured
ambient charge scale, exceeds the applicability limit.

This switch enables diagnostics and the applicability check only. It does not
change particle return or escape, which remain controlled by
`external_boundary.field.model` and `external_boundary.particles.mode`. Every
`photo_raycast` species using tracked outer transfer must set
`deposit_opposite_charge_on_emit=true`.

## Check convergence of photoelectron emission

Increase `rays_per_batch` and verify convergence of hit fraction, emitted current, and charging distribution. When reabsorption
location matters, reduce `dt` and verify that the result is unchanged. See [Inspect Output Files](OutputGuide.en.html) for the
species-resolved charge balance across emission, absorption, and escape and for diagnostics specific to outer return.

## Code reference

- Ray propagation, hit, emission velocity, and weight: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- Source creation and emission charge difference: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- Outer-transfer and histogram compatibility validation: [`bem_app_config_parser.f90`](../src/config/app_config_parser/bem_app_config_parser.f90)
- Kinetic mean photoelectron density: [`bem_outer_plasma_kinetic.f90`](../src/physics/outer_plasma/bem_outer_plasma_kinetic.f90)
- Kinetic mean runtime assembly: [`bem_outer_plasma_kinetic_runtime.f90`](../src/runtime/bem_outer_plasma_kinetic_runtime.f90)
- Photoelectron histogram and applicability check: [`bem_outer_plasma_photoelectron.f90`](../src/physics/outer_plasma/bem_outer_plasma_photoelectron.f90)

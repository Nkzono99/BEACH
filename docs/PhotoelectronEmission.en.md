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
5. Track it as an ordinary particle through reabsorption, infinity escape, or outer return.
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

Normal velocity is positive, so a new particle leaves the surface toward the illumination side. A tracked orbit or selected
reduced closure then decides whether it returns or escapes.

## Record emission and reabsorption in surface charge

With `deposit_opposite_charge_on_emit=true`, source element $i$ receives

$$
\Delta q_{i,\mathrm{emit}}=-q w.
$$

For an electron $q<0$, this leaves positive charge. The difference is accumulated separately as `photo_emission_dq`, combined
with collision deposition after MPI all-reduce, and included in the same batch commit.

If the emitted particle is reabsorbed on element $j$, ordinary absorption deposits $+qw$ there. Return to the source element
cancels emission charge; return to another element transfers net charge across the surface. The current insulator model applies
no subsequent lateral surface conduction.

## Represent motion outside the box with one model

Exactly one of the following approximations applies to a photoelectron moving out of the box. The selected configuration
determines the tracked weight and how much of the return position and flight time is resolved.

| Configuration | Created tracked weight | Representation of return |
| --- | --- | --- |
| `photo_escape_model="none"`, no outer transfer | $w_\mathrm{hit}$ | Track inside the box; ordinary escape at an open face |
| `photo_escape_model="boltzmann_cutoff"` | $f_\mathrm{esc}w_\mathrm{hit}$ | Immediate reduced closure that never creates the nonescaping fraction |
| 1-D outer-profile return | $w_\mathrm{hit}$ | Classify individual interface crossings by conserved energy |
| Unified 3-D explicit orbit | $w_\mathrm{hit}$ | Follow individual trajectories in the 3-D outer zero/nonzero-mode field |

Legacy Boltzmann cutoff evaluates the local emission potential $\phi_\mathrm{emit}$ without the primary self term:

$$
f_\mathrm{esc}=
\exp\left[-\frac{|q|\max(\phi_\mathrm{emit}-\phi_\infty,0)}{k_\mathrm{B}T_\mathrm{PE}}\right].
$$

It multiplies particle weight by this factor; at $T_\mathrm{PE}=0$ a positive barrier gives zero. It has no returning-particle
trajectory, reabsorption location, or flight time and does not solve which element receives the nonescaping fraction. Opposite
source charge uses the same reduced weight.

Configurations with tracked individual return require `deposit_opposite_charge_on_emit=true` and reject a legacy non-`none`
`photo_escape_model`, preventing duplicate application of an escape cutoff and tracked return.

## Connect outer-plasma mean density to tracked particles

`outer_plasma.photoelectron_closure="kinetic_mean"` uses the temperature and emission current density of the first negative
`photo_raycast` species as a plane-averaged source in the 1-D Poisson density closure. Its outgoing and turning-return populations
contribute to outer space charge, but it neither replaces surface absorption of tracked particles nor deposits an extra
statistical return current on the surface.

With `individual_return` or `kinetic_mean + kinetic_1d_profile_return`, tracked photoelectrons crossing z-high are classified as
escape or return by conserved energy. Interface outward and returned charge and an emission-velocity histogram are diagnostics.
See [Particle escape and return](ParticleEscapeReturn.en.html) for the quasi-steady approximation that omits outer flight from
global time and the 3-D explicit orbit, and [Outer plasma models](OuterPlasmaModels.en.html) for field construction.

Zhao models supply emission current density, normal cutoff, and drift according to the selected branch. Tracked particles
advance in the ordinary field snapshot rather than the Zhao profile $E(z)$.

## Check charge balance across emission, reabsorption, and escape

The species-resolved ledger separates emission from the surface, surface absorption, infinity escape, unresolved discard, and
gross outward and returned charge at the outer interface. At minimum, compare:

- hit fraction, emission current, and charging distribution as `rays_per_batch` increases;
- reabsorption location and escape/return fraction as `dt` decreases;
- charge balance among surface emission $-qw$, tracked absorption or escape, and surface stock between batches;
- maximum flight time, frozen-field ratio, and energy error when outer return is active.

## Code reference

- Ray propagation, hit, emission velocity, and weight: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- Reduced escape factor and source charge difference: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- Tracked-return compatibility validation: [`bem_app_config_parser.f90`](../src/config/app_config_parser/bem_app_config_parser.f90)
- Ordinary tracking, outer transfer, and charge-balance accounting: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- Kinetic mean photoelectron closure: [`bem_outer_plasma_photoelectron.f90`](../src/physics/outer_plasma/bem_outer_plasma_photoelectron.f90)

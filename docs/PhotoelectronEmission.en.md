title: Photoelectron emission and lifecycle

Lang: [日本語](PhotoelectronEmission.md) | [English](PhotoelectronEmission.en.md)

# Photoelectron emission and lifecycle

`source_mode="photo_raycast"` emits particles from the first surface hit by an illumination ray. This page covers ray casting,
emission amount and velocity, source charge, and the closed-PE `neutral_return` closure. After creation, a photoelectron uses
the ordinary fixed field, Boris update, collision query, and box boundaries.

## Track emission through reabsorption in the same batch

1. Launch rays from an illumination aperture on a box face.
2. Find the first mesh hit while applying box boundary conditions.
3. Create a photoelectron from the plasma-facing normal of the hit element.
4. Optionally record $-qw$ on the emission element.
5. Track it as an ordinary particle and hand it to the common escape/return treatment after it reaches a box boundary.
6. Commit emission and absorption charge to the surface at batch end.

Emission and reabsorption can occur in one batch, but the field is not refreshed between them. Net surface charge changes the
field starting in the next batch.

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
| Nonperiodic face | Leave the box and terminate without emission |
| Face on `domain.periodic_axes` | Wrap to the opposite face and continue |

With `field_boundary.mode="periodic2"`, first-hit search includes periodic images and wraps the emission position into the primary
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

`rays_per_batch` is the sample count for the illumination integral, not the physical emission amount. Increasing it reduces
$w_\mathrm{hit}$ and Monte Carlo noise in visibility and emission location. Check convergence with ray count.

## Build the emitted state from the surface normal

If the stored triangle normal points along the incident ray, it is reversed to form emission normal $\mathbf n_s$ facing back
toward the illumination side. Position is offset $10^{-12}$ m along $\mathbf n_s$ to avoid immediate self-intersection with the
source element.

With $\sigma=\sqrt{k_\mathrm{B}T/m}$, velocity is sampled in local basis
$(\mathbf n_s,\mathbf t_1,\mathbf t_2)$:

- normal speed follows a flux-weighted half-range Maxwell distribution with `normal_drift_speed`;
- both tangential components are zero-mean Gaussians with standard deviation $\sigma$;
- Gaussian sampling is truncated at $6\sigma$.

Normal velocity is positive, so a new particle leaves toward the illumination side. Its tracked orbit and common box boundaries
then determine reabsorption, escape, or local reflection.

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

A ray hit always creates a photoelectron with weight $w_\mathrm{hit}$. No emission-time setting multiplies the weight by an
escape fraction. Surface return is an ordinary collision; an open-face event uses the same
`particle_boundary.ordinary_open_model` treatment as other sources.

### Close only photoelectrons at the injection face

Use local reflection to evaluate photoelectron surface redistribution on closed orbits:

```toml
[particle_boundary]
z_high = "open"
ordinary_open_model = "escape"

[[particles.species]]
species_key = "photoelectron"
source_mode = "photo_raycast"
inject_face = "z_high"
deposit_opposite_charge_on_emit = true
surface_charge_closure = "neutral_return"

[particles.species.boundary]
z_high = "reflect"
```

In this example, a z-high crossing reverses only normal velocity, preserves tangential velocity, and reintegrates the step
remainder. Ambient species using the default `inherit` in `[particles.species.boundary]` follow the global open contract.
All six species faces accept `inherit`, `open`, `reflect`, or `redistributed_reflect`. Closed PE requires the effective action
on the face named by `inject_face` to be `reflect` or `redistributed_reflect`; a face on `domain.periodic_axes` cannot be
overridden.

Ordinary `reflect` also preserves tangential position. Only when the in-plane position of returning photoelectrons should be
uniformly redistributed, replace the baseline value above with:

```toml
[particles.species.boundary]
z_high = "redistributed_reflect"
```

`redistributed_reflect` applies ordinary reflection to velocity and, for a single face, uniformly resamples only the two
in-plane coordinates over the box span excluding its end guards. This option generalizes the horizontal-position randomization
used for top-boundary photoelectron return by
[Zimmerman et al. (2016)](https://doi.org/10.1002/2016JE005049) to any nonperiodic face and simultaneous face event. It does not
add a self-consistent outer sheath. See
[Particle collision and boundary events](ParticleEvents.en.html) for simultaneous-event rules.

Species-level reflection alone closes only the orbit. A particle that does
not return by `max_step` remains unresolved. With
`surface_charge_closure="neutral_return"`, BEACH measures emitted charge $S<0$
and resolved absorbed charge $R<0$ globally for the batch, then multiplies
deposits at measured return destinations by $S/R$. Combined with source
reaction charge, the photoelectron contribution to total surface charge is
exactly zero. The approximation assigns unresolved particles the same
destination distribution as resolved returns.

Raw `absorbed_on_surface_C` and `discarded_unresolved_C` are not replaced.
`charge_ledger.csv` records the correction,
`neutral_return_weight_scale`, and `neutral_return_unresolved_fraction`
separately. A run stops on nonzero emission without resolved return, actual
escape, `soft_discard`, or a charge-sign mismatch.
An unresolved fraction above 5% also stops before correction because it is
outside this closure's fixed applicability range.

Full reflection is an artificial mirror at the injection face of a finite box, not a self-consistent sheath or quasineutrality solution.
`neutral_return` is a statistical zero-net-photoelectron-current closure, not a resolved trajectory for every long-lived
particle. Converge `max_step`, `dt`, ray count, and batch width until `abs(weight_scale-1)` and the unresolved fraction are small. See
[Finite-image periodic2 configuration](FinitePeriodicConfiguration.en.html)
for the integrated field, inflow, and potential-reference setup.

For coupling to an external current model, use `surface_charge_closure="fixed_current"` and specify
`target_emission_current_a` and `target_absorbed_current_a` as separate signed surface currents. BEACH uniformly rescales the
raycast emission-origin distribution and the tracked return-destination distribution independently; it never uses the small
net PE current, the difference of these two large channels, as a scale denominator. If the external model injects a return
VDF, keep the top face open and do not also enable `neutral_return` on the same or another return channel.

With `[surface_current_model] model="zhao_stationary"`, BEACH independently calculates PE emission, escape, and return
from the Zhao zero-current stationary root. Emission and return are applied to surface channels; escape is recorded as
an external-boundary target that is not deposited onto the surface. The closure never uses the PE net current as a
scaling denominator, so it remains stable when large currents nearly cancel. Raw escape statistics are preserved for
comparison with target / applied / correction values.

This stability concerns total-current normalization, not the statistical accuracy of the raw spatial map. If there is one
return hit, that element receives the entire return target. BEACH imposes no fixed minimum hit count, so inspect the raw
ledger counts and `fixed_*_weight_scale`, then test the elementwise distribution across `rays_per_batch`, batch widths, and
RNG seeds. Closure of the Zhao zero-current budget does not replace this convergence check.

The emission VDF remains the configured surface half-Maxwellian. Keep z-high open; the outward crossing test uses $\phi_m$
for Type A and 0 V for Type B/C as the outside barrier potential, together with the local crossing potential and normal
kinetic energy, to split return from escape. This reflection contracts the outer turning point onto z-high: BEACH does not
solve the outside sheath field, space charge, return distance, or return delay. See
`examples/periodic2_zhao_fixed_current.toml` for the complete setup.

See [Particle escape and local return](ParticleEscapeReturn.en.html) to choose between
`particle_boundary.ordinary_open_model` and closed PE.

## Check convergence of photoelectron emission

Increase `rays_per_batch` and verify convergence of hit fraction, emitted current, and charging distribution. When reabsorption
location matters, reduce `dt` and verify that the result is unchanged. See [Inspect Output Files](OutputGuide.en.html) for the
species-resolved charge balance and closed-PE correction diagnostics.

## Code reference

- Ray propagation, hit, emission velocity, and weight: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- Source creation and emission charge difference: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)

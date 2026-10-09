title: Photoelectron emission and charge closures

Lang: [English](PhotoelectronEmission.en.md) | [日本語](PhotoelectronEmission.md)

# Photoelectron emission and charge closures

Photoelectrons are emitted from illuminated surfaces and either return to a surface in the field of the cell or leave the box.
BEACH traces illumination rays to find the emission sites (`source.mode="photo_raycast"`) and tracks the emitted photoelectrons
in the same field and with the same method as any other particle. This page explains how the emitted amount and the emission
velocity are set, the reaction charge left on the emitter, and the photoelectron charge closures (closed photoelectrons and
fixed current).

## From emission to absorption

1. Rays start from an illumination aperture placed on a box face.
2. Find the first triangle each ray hits.
3. Emit a photoelectron from the hit element toward the illuminated side.
4. Record the reaction charge on the emitting element.
5. Track the photoelectron like any other particle. When it reaches a box face, the treatment of that face applies ([Particles at box boundaries](ParticleEscapeReturn.en.html)).
6. Commit the emitted and absorbed charge to the surface at the end of the batch.

Emission and re-absorption can occur within one batch, but the field is not changed during the batch. The net surface charge
affects the field from the next batch.

## Find the emitting surface by illumination

`source.inject_face` and `source.pos_low` / `source.pos_high` define a rectangular aperture on a box face. `source.ray_direction`
must point from the aperture into the box; if omitted it is the inward normal of the face. With aperture area $A$, inward normal
$\mathbf n_\mathrm{in}$, and unit ray direction $\hat{\mathbf d}$, the area projected perpendicular to the rays is

$$
A_\mathrm{proj}=A\left|\hat{\mathbf d}\cdot\mathbf n_\mathrm{in}\right|
$$

Ray origins are sampled uniformly in the aperture. Each ray searches for the first triangle it hits in each segment up to the
next box face. A ray that reaches a non-periodic face ends without emission; a ray that reaches a periodic face wraps to the
opposite side and continues. In a periodic field the first hit is searched including periodic images, and the emission position
is mapped back into the cell. A ray that does not hit within `particles.raycast.max_bounces` wraps creates no particle.

## Emission current per ray

With emission current density $J_\mathrm{emit}>0$, physical particle charge $q$, total ray count over all MPI ranks
$N_\mathrm{ray}$, and batch duration $\Delta t_\mathrm{batch}$, the macro-particle weight created by a hitting ray is

$$
w_\mathrm{hit}
=\frac{J_\mathrm{emit}A_\mathrm{proj}\,\Delta t_\mathrm{batch}}
{|q|N_\mathrm{ray}}
$$

Rays that miss create no particle, so shadowing and apparent area enter the emitted amount through the hit fraction.
`sampling.rays_per_batch` is a sample count, not an emitted amount; increasing it lowers $w_\mathrm{hit}$ and reduces the
statistical error of the emission positions. `sampling.weight` and `sampling.target_macro_particles_per_batch` are not used.

## Emission velocity

Of the two normals of the hit triangle, the one facing away from the incident ray is the emission normal $\mathbf n_s$. The
emission position is shifted by $10^{-12}$ m along $\mathbf n_s$ from the hit point to avoid an immediate re-collision with the
same element. From the temperature, $\sigma=\sqrt{k_\mathrm{B}T/m}$ is computed and the velocity is sampled in the local basis
$(\mathbf n_s,\mathbf t_1,\mathbf t_2)$.

- Normal velocity: a flux-weighted half-Maxwellian with drift `source.normal_drift_speed_m_s`.
- Two tangential components: Gaussian with zero mean and standard deviation $\sigma$ (truncated at $6\sigma$).

## Reaction charge

If `source.deposit_opposite_charge_on_emit=true`, the emitting element $i$ receives

$$
\Delta q_{i,\mathrm{emit}}=-q w
$$

For electrons $q<0$, so positive charge remains on the surface. When the photoelectron is absorbed by element $j$, $+qw$ is
added to $j$ as an ordinary absorption. Returning to the same element cancels the two; returning to another element moves
charge within the surface.

## Photoelectron charge closures

Because the cell is shorter than the Debye length, many photoelectrons that reach the top face would in reality be pushed back by
the outer sheath. With an open top face, photoelectron escape is overcounted. Choose how to compensate with the species
`charging.closure`.

| `charging.closure` | Treatment of photoelectrons | When to use |
|---|---|---|
| Omitted | As tracked. Only the treatment of the top face decides | Baseline for comparison |
| `neutral_return` | Reflect at the top face and match the net photoelectron current to zero | Study redistribution within the surface |
| `fixed_current` | Match emission and return currents to targets given from outside | Combine with an external current model |

For periodic-surface cases, the choice is compared in [Set up a periodic surface in plasma](PeriodicPlasmaSurface.en.html).

### Closed photoelectrons (`neutral_return`)

In the photoelectron species, make the illuminated face (`source.inject_face`) reflecting and set `neutral_return`.

```toml
[particles.species.charging]
closure = "neutral_return"

[particles.species.boundary]
z_high = "reflect"
```

Reflection reverses only the normal velocity and keeps the tangential velocity and the position. With `"redistributed_reflect"`
the velocity is reversed in the same way and only the return position is chosen again uniformly on the face (the same idea as
the top-boundary photoelectron return of [Zimmerman et al. (2016)](https://doi.org/10.1002/2016JE005049)).

Even with reflection, some photoelectrons do not return within the per-particle step limit. `neutral_return` sums the emitted
photoelectron charge $S<0$ and the charge absorbed by the surface (returned) $R<0$ of one batch over all MPI ranks, and scales
the charge deposited at the return sites by $S/R$.

$$
(-S)+\frac{S}{R}R=0
$$

Together with the reaction charge on the emitters, the change of total surface charge by photoelectrons is exactly zero.
Unreturned photoelectrons are approximated as returning with the same distribution as those that returned in the same batch.
Although the total charge is zero, charge moved to surfaces at different heights still leaves a plane-averaged vertical dipole.

The batch is not accepted and the run stops in the following cases.

- There is emission but no return.
- A photoelectron escaped through an open face or was discarded as a particle crossing several box faces at once (soft discard).
- A value is not finite or the signs are inconsistent.
- The unreturned fraction exceeds 5%.

### Fixed current (`fixed_current`)

Match the emission and return currents to targets set by an external model. The targets are signed as contributions to surface
charging and are given separately.

```toml
[particles.species.charging]
closure = "fixed_current"
target_emission_current_a = 4.5e-15
target_absorbed_current_a = -3.7e-15
```

BEACH scales the emitter distribution and the return distribution uniformly and separately to the targets. The net current, a
difference of two large currents, is never used as the denominator of a scale, so the closure stays stable when emission and
return nearly cancel. Keep the top face open and do not combine it with `neutral_return`. With the outer-sheath connection
([Zhao stationary sheath](ZhaoStationaryClosure.en.html)), the targets are set automatically from the zero-current root of the
outer sheath.

A stable scale does not make the per-element distribution statistically accurate. If only one return is tracked, the whole
return target is deposited on that element.

## Outputs to check

| Output | What to check |
|---|---|
| `charge_ledger.csv` | Emitted, absorbed, and escaped charge and counts. For `neutral_return`, `neutral_return_weight_scale` and `neutral_return_unresolved_fraction`; for `fixed_current`, `fixed_*_weight_scale` |
| `fixed_current_history.csv` | Time variation of the `fixed_current` scale (`target_over_tracked`) |

In a run whose scale is far from 1, the closure assumption rather than tracking sets the charge distribution.

## Check convergence

- Increase `sampling.rays_per_batch` and confirm that the hit fraction, emission current, and charging distribution do not change.
- If you evaluate return positions, reduce `particles.tracking.dt_s`, increase `max_steps_per_particle`, and confirm that they do not change.
- For `neutral_return`, revise the step limit and the box height until the scale is near 1 and the unreturned fraction is small.

## Scope

- The emission velocity is the configured surface half-Maxwellian. Surface material and work-function distributions are not modeled.
- Reflection at the top face (`neutral_return`) places a mirror above a finite box; it does not solve the outer sheath or quasi-neutrality.
- On an insulating surface, returned charge stays on its element and does not conduct along the surface.

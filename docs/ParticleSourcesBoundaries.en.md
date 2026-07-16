title: Particle-source overview

Lang: [日本語](ParticleSourcesBoundaries.md) | [English](ParticleSourcesBoundaries.en.md)

# Particle-source overview

Use this page to choose `source_mode` from where particles originate and to understand the processing shared after creation.
The numerical details and physical models for each source are documented on their dedicated pages.

## Choose a source for the intended task

Set one `source_mode` for every `[[particles.species]]` entry.

| `source_mode` | What determines particle count | Creation location | Suitable use |
| --- | --- | --- | --- |
| `volume_seed` | `npcls_per_step` | Between `pos_low` and `pos_high` | Initial populations, orbit tests, and specified counts |
| `reservoir_face` | Inflow flux, aperture area, and `batch_duration` | A selected box face | Continuous inflow from an external reservoir |
| `photo_raycast` | Current density, projected area, `batch_duration`, and ray count | First surface hit by a ray | Illumination-driven surface emission |

Choose `volume_seed` to specify a particle count directly, `reservoir_face` for physical inflow across a face, or
`photo_raycast` for emission from an illuminated surface.

## Create a specified population with `volume_seed`

`volume_seed` creates `npcls_per_step` particles in every batch. Position is uniform in the rectangular volume
`[pos_low, pos_high]`, and velocity follows the drifting Maxwell distribution

$$
\mathbf v=\mathbf u+\sigma\mathbf Z,
\qquad
\sigma=\sqrt{\frac{k_\mathrm{B}T}{m}}.
$$

If `thermal_speed` is supplied, it overrides the $\sigma$ derived from temperature. Each standard-normal component is truncated
at $6\sigma$.

This mode specifies particle count directly; it does not derive the count from a physical surface flux.

## Create continuous external inflow with `reservoir_face`

`reservoir_face` converts an upstream velocity distribution function (VDF) outside the box into particles entering through a
selected face. Physical flux determines the count in each batch, and normal velocity is drawn from a flux-weighted distribution
appropriate for particles crossing a surface.

After selecting this mode, use the following pages for the distinct calculations:

- Convert the upstream VDF into inflow, macro-particle count, initial position, and face-arrival velocity:
  [`reservoir_face` inflow and velocity sampling](ReservoirInjection.en.html)
- Select a consistent combination of `infinity_barrier`, outer plasma, and outflow boundary:
  [Selecting boundary and outer-domain models](OuterPlasmaModels.en.html)
- Correct the upstream VDF with a Zhao-family closure or `floating_no_photo`:
  [Inflow-VDF sheath closures](SheathInjectionClosures.en.html)

`reservoir_face` itself does not solve an outer sheath or particle trajectories outside the box.

## Emit from illuminated surfaces with `photo_raycast`

`photo_raycast` launches rays from an illumination aperture on a box face and propagates them under the box boundary conditions.
It emits particles from the first in-box element hit and samples velocity from a flux-weighted Maxwell distribution relative to
the element normal.

See [Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for the opposite charge left at the emission source and
the photoelectron-specific reduced escape model.

## Enter the same particle tracker after creation

All source modes enter the same particle state and tracking path after creation. The main state comprises position $\mathbf x$,
velocity $\mathbf v$, physical-particle charge $q$ and mass $m$, macro-particle weight $w$, and species ID. The charge of one
tracked macro particle is $qw$.

| Outcome within a batch | Treatment |
| --- | --- |
| Absorbed by the mesh | Deposit $qw$ on the hit element |
| Escapes to infinity through an open face | Remove and include in species-resolved escape |
| Returns from an outer region | Map the same particle to the interface and reintegrate the remaining step |
| Alive after `max_step` | Discard and report as unresolved at batch end |

See [Boris particle update](BorisPusher.en.html) for particle advancement,
[particle collision and boundary events](ParticleEvents.en.html) for mesh-versus-box ordering, and
[particle escape and return](ParticleEscapeReturn.en.html) for behavior outside the box.

## Evaluate sources from the batch-start field

Sources are evaluated after the field, potential, and outer-plasma state are refreshed at the start of a batch. Reservoir velocity
corrections and photoelectron escape fractions therefore see surface charge committed through the preceding batch. Created
particles move through the snapshot fixed for the current batch.

The only source operation that changes surface charge at creation is the optional opposite charge left by `photo_raycast`.
That difference is committed at batch end with absorbed charge and does not alter the current batch field. See
[Model overview](Algorithms.en.html) for the complete order.

## Preserve generated amounts across MPI and restart

The global `reservoir_face` count and macro-particle remainder are determined once on the root rank and then split across ranks.
`photo_raycast` `rays_per_batch` is also a global total. Expected inflow and emission therefore do not depend on MPI world size,
although random sequences and individual trajectories may change with world size.

The per-species reservoir `macro_residual` is saved in `macro_residuals.csv` and restored on resume. See
[`reservoir_face` inflow and velocity sampling](ReservoirInjection.en.html) for the residual calculation and diagnostics.

## Code reference

- Particle distributions and ray casting: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- Per-source batch creation and macro-particle residuals: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- Source input validation: [`bem_app_config_parser_validate.f90`](../src/config/app_config_parser/bem_app_config_parser_validate.f90)
- Charge-balance accounting and batch tracking: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- Macro-particle residual checkpoints: [`bem_restart.f90`](../src/runtime/bem_restart.f90)

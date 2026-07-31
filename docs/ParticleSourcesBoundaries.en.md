title: Particle-source overview

Lang: [日本語](ParticleSourcesBoundaries.md) | [English](ParticleSourcesBoundaries.en.md)

# Particle-source overview

Choose `source_mode` for each `[[particles.species]]` from its origin and governing physical quantity. Tracking is shared after creation.

## Choose a source for the intended task

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

This mode does not derive count from a physical surface flux.

## Create continuous external inflow with `reservoir_face`

`reservoir_face` converts a supplied reservoir VDF into particles entering through a selected face. Physical flux determines the
count in each batch, and normal velocity is drawn from a flux-weighted distribution
appropriate for particles crossing a surface.

See [`reservoir_face` inflow and velocity sampling](ReservoirInjection.en.html) for flux, macro-particle count, initial
position, and face-arrival velocity. See [particle escape and local return](ParticleEscapeReturn.en.html) to combine
`infinity_barrier` with open-face treatment. `reservoir_face` does not solve outside-box trajectories or a self-consistent sheath.

## Emit from illuminated surfaces with `photo_raycast`

`photo_raycast` launches rays from an illumination aperture on a box face and propagates them under the box boundary conditions.
It emits particles from the first in-box element hit and samples velocity from a flux-weighted Maxwell distribution relative to
the element normal.

See [Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for source reaction charge, reabsorption, and closed PE.

## Enter the same particle tracker after creation

All source modes enter the same particle state and tracking path after creation. The main state comprises position $\mathbf x$,
velocity $\mathbf v$, physical-particle charge $q$ and mass $m$, macro-particle weight $w$, and species ID. The charge of one
tracked macro particle is $qw$.

| Outcome within a batch | Treatment |
| --- | --- |
| Absorbed by the mesh | Deposit $qw$ on the hit element |
| Escapes through an open face | Remove and include in species-resolved escape |
| `reflect`, `redistributed_reflect`, or periodic box action | Reintegrate the remaining step for the same particle |
| Alive after `max_step` | Discard and report as unresolved at batch end |

Configure global nonperiodic actions in `[particle_boundary]`, species overrides in `[particles.species.boundary]`, and periodic
axes in `domain.periodic_axes`. See [Boris particle update](BorisPusher.en.html) for advancement,
[particle collision and boundary events](ParticleEvents.en.html) for mesh-versus-box ordering, and
[particle escape and local return](ParticleEscapeReturn.en.html) for open-face treatment.

## Evaluate sources from the batch-start field

Sources are evaluated after the field and potential snapshot is built at batch start. Reservoir velocity correction and particle
motion therefore see surface charge committed through the preceding batch and use one fixed snapshot during the current batch.

The only source operation that creates a surface-charge delta is the optional opposite charge left by `photo_raycast`.
It is committed at batch end with absorbed charge and does not alter the current field. See
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

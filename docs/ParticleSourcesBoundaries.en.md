title: Particle sources

Lang: [日本語](ParticleSourcesBoundaries.md) | [English](ParticleSourcesBoundaries.en.md)

# Particle sources

Particle tracking in each batch starts by creating a new population for every `[[particles.species]]` entry. Once created,
all sources enter the same particle state and use the [Boris particle update](BorisPusher.en.html) and
[particle-event](ParticleEvents.en.html) processing.

| `source_mode` | What determines particle count | Creation location | Main use |
| --- | --- | --- | --- |
| `volume_seed` | `npcls_per_step` | Between `pos_low` and `pos_high` | Initial particles and orbit tests |
| `reservoir_face` | Inflow flux, area, and `batch_duration` | A selected box face | Continuous inflow from an external reservoir |
| `photo_raycast` | Current density, projected area, `batch_duration`, and ray count | First surface hit by a ray | Illumination-driven surface emission |

## Create the batch particles after refreshing the snapshot

Sources are evaluated after the field and outer-plasma snapshots have been refreshed at the start of a batch. Reservoir velocity
corrections and reduced photoelectron escape fractions therefore see surface charge committed through the preceding batch.
Created particles move through that same snapshot until absorption, escape, or `max_step`.

The only source operation that changes surface charge at creation is the optional opposite charge left at a `photo_raycast`
emission element. That difference is committed at the end of the batch together with absorption charge and does not alter the
field seen during the current batch. See [Model overview](Algorithms.en.html) for the complete order.

## Create a specified initial population with `volume_seed`

`volume_seed` creates `npcls_per_step` particles in every batch. Position is uniform in the rectangular volume
`[pos_low, pos_high]`, and velocity is a drifting Maxwell distribution,

$$
\mathbf v=\mathbf u+\sigma\mathbf Z,
\qquad
\sigma=\sqrt{\frac{k_\mathrm{B}T}{m}}.
$$

If `thermal_speed` is supplied it overrides the $\sigma$ derived from temperature. Each standard-normal component is truncated at
$6\sigma$.

`npcls_per_step` directly sets the injected population for this source. Use `reservoir_face` for continuous inflow derived from
a physical flux.

## Build a continuous inflow flux with `reservoir_face`

`reservoir_face` converts a density, temperature, and drift or a velocity grid outside the box into flux crossing a selected face
inward. Particle count is derived from physical flux, and normal velocity is sampled from a flux-weighted distribution. In a
configuration with a potential difference between the upstream reservoir and the face, the same difference selects accessible
upstream particles and maps their velocity to the face by energy conservation.

See [Reservoir injection](ReservoirInjection.en.html) for equations, velocity grids, macro-particle residuals,
`infinity_barrier`, and coupling to an outer profile.

## Emit particles from illuminated surfaces with `photo_raycast`

`photo_raycast` launches rays from an illumination aperture on a box face, propagates them under box boundary conditions, and
emits from the first in-box element hit. Emission velocity follows a flux-weighted Maxwell distribution relative to the element
normal.

After emission, a photoelectron enters the common particle state and uses the same field, Boris update, mesh collisions, and box
boundaries as every other particle. Source charge, reduced escape closure, and outer-sheath return are documented together in
[Photoelectron emission and lifecycle](PhotoelectronEmission.en.html).

## Enter the common particle state and charge ledger after creation

After creation, the main stored quantities are position $\mathbf x$, velocity $\mathbf v$, physical-particle charge $q$ and mass
$m$, macro-particle weight $w$, and species ID. The charge of a tracked macro particle is $q w$; surface absorption and the charge
ledger use that value.

| Batch outcome | Treatment |
| --- | --- |
| Absorbed by the mesh | Deposit $q w$ on the hit element |
| Escapes to infinity through an open face | Remove and count in species-resolved escape |
| Returns from an outer region | Map the same particle to the interface and reintegrate the step remainder |
| Alive after `max_step` | Discard and report as unresolved at batch end |

See [Particle collision and boundary events](ParticleEvents.en.html) for mesh-versus-box ordering,
[Outer plasma models](OuterPlasmaModels.en.html) for the external field model,
[Particle escape and return](ParticleEscapeReturn.en.html) for the particle map, and
[Surface charge update](SurfaceModels.en.html) for charge commit.

## Preserve global generation amounts across MPI and restart

The global `reservoir_face` count and macro-particle remainder are determined once on the root and then split across ranks, so the
expected inflow does not change with MPI world size. A per-species `macro_residual` is saved in `macro_residuals.csv` and restored
on resume.

`photo_raycast` `rays_per_batch` is also a global total. Each ray weight is divided by the global ray count, and emitted,
absorbed, and escaped charge-ledger values are global after MPI all-reduce. Expected inflow is world-size independent, while
random sequences and individual trajectories can change with world size.

## Code reference

- Particle distributions and ray casting: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- Per-source batch creation and macro residuals: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- Source input validation: [`bem_app_config_parser_validate.f90`](../src/config/app_config_parser/bem_app_config_parser_validate.f90)
- Charge ledger and batch tracking: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- Macro-residual checkpoint: [`bem_restart.f90`](../src/runtime/bem_restart.f90)

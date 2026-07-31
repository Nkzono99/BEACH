title: Particle sources and boundary inflow

Lang: [日本語](ParticleSourcesBoundaries.md) | [English](ParticleSourcesBoundaries.en.md)

# Particle sources and boundary inflow

`source_mode` selects how particles are created inside the simulation domain.
`[particles.species.boundary_inflow]` adds particles crossing the boundary from outside. These settings have separate roles.

## Choose a setting from the creation location

| Setting | What determines particle count | Creation location | Suitable use |
| --- | --- | --- | --- |
| `source_mode="volume_seed"` | `npcls_per_step` | Volume from `pos_low` to `pos_high` | Initial populations, orbit tests, and specified counts |
| `source_mode="plane_source"` | Inflow flux, rectangle area, and `batch_duration` | Axis-aligned rectangle inside the box | Explicit one-way internal surface source |
| `source_mode="photo_raycast"` | Current density, projected area, `batch_duration`, and ray count | First surface hit by a ray | Illumination-driven surface emission |
| `[particles.species.boundary_inflow]` | Reservoir flux, box-face area, and `batch_duration` | Nonperiodic simulation boundary | Continuous inflow from an external plasma reservoir |

`boundary_inflow` is not a `source_mode`. For a boundary-only species, use `source_mode="volume_seed"` with
`npcls_per_step=0`; for a Maxwell distribution, set `npcls_per_step>0` if the same species also needs an initial population.
Velocity-grid boundary inflow cannot use a positive `npcls_per_step`. In the initial
implementation, a flux-driven `plane_source` or legacy `reservoir_face` cannot be combined with `boundary_inflow` for the
same species.

## Create a specified population with `volume_seed`

`volume_seed` creates `npcls_per_step` particles in every batch. Position is uniform in the rectangular volume
`[pos_low, pos_high]`, and velocity follows a drifting Maxwell distribution. If `thermal_speed` is supplied, it overrides
the value derived from temperature.

This source does not derive count from physical surface flux. `npcls_per_step=0` is valid when the same species uses only
`boundary_inflow`. A positive value can accompany Maxwell boundary inflow, but not velocity-grid inflow.

## Emit from an internal rectangle with `plane_source`

`plane_source` creates particles from the axis-aligned rectangle selected by `pos_low` and `pos_high`, directed along
`source_normal`.

```toml
[[particles.species]]
source_mode = "plane_source"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.1093837139e-31
target_macro_particles_per_batch = 300
pos_low = [0.2, 0.2, 2.0]
pos_high = [0.8, 0.8, 2.0]
source_normal = [0.0, 0.0, -1.0]
```

`pos_low` and `pos_high` must be equal on exactly one axis and have positive extent on the other two. The normal coordinate
must be strictly between the box faces, while tangential ranges may reach the box bounds. `source_normal` is a nonzero vector
along the zero-thickness axis. Positive or negative
unit vectors are recommended in configurations; the implementation normalizes magnitudes such as `[2,0,0]`.

Particle count comes from Maxwell reservoir density and temperature, or the specified velocity-grid flux, multiplied by
rectangle area and `batch_duration`. Normal velocity follows the one-way flux distribution along `source_normal`.
`[reservoir]` settings `infinity_barrier`, `phi_infty`, and `face_potential_grid_n` do not apply to an internal plane.

## Inflow from an external reservoir with `boundary_inflow`

Set one or more of the six faces in `[particles.species.boundary_inflow]` to `"reservoir"`.

```toml
[[particles.species]]
species_key = "solar_wind_electron"
source_mode = "volume_seed"
npcls_per_step = 0
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.1093837139e-31
target_macro_particles_per_batch = 300
drift_velocity = [0.0, 0.0, -4.0e5]

[particles.species.boundary_inflow]
z_high = "reservoir"
```

Particles enter across the complete selected box face. Multiple nonperiodic faces can be enabled; periodic faces cannot.
Configure outward actions independently in `[particle_boundary]` or `[particles.species.boundary]`; each inflow face must
have an effective `open` action.

See [Reservoir inflow through simulation boundaries](ReservoirInjection.en.html) for flux, macro-particle count,
per-face remainders, and `infinity_barrier`. This model does not solve outside-box trajectories or a self-consistent sheath.

## Emit from illuminated surfaces with `photo_raycast`

`photo_raycast` launches rays from an illumination aperture on a box face and propagates them under the box boundary
conditions. It emits particles from the first in-box element hit and samples velocity from a flux-weighted Maxwell
distribution relative to the element normal.

See [Photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for source reaction charge, reabsorption, and closed PE.

## Treat legacy `reservoir_face` as compatibility input

`source_mode="reservoir_face"` is a deprecated input whose behavior is retained for existing cases. `inject_face` and
`pos_low` / `pos_high` select an aperture on a box face. BEACH does not silently convert it to `boundary_inflow` or
`plane_source`.

Use `boundary_inflow` for new external-plasma cases and `plane_source` for explicit internal rectangles. When moving a
legacy aperture to complete-face boundary inflow, account for the change in physical inflow caused by the area change.

## Enter the same particle tracker after creation

All creation paths enter the same particle tracker.

| Outcome within a batch | Treatment |
| --- | --- |
| Absorbed by the mesh | Deposit macro-particle charge on the hit element |
| Escapes through an open face | Remove and include in species-resolved escape |
| `reflect`, `redistributed_reflect`, or periodic box action | Reintegrate the remaining step for the same particle |
| Alive after `max_step` | Discard and report as unresolved at batch end |

Configure global nonperiodic actions in `[particle_boundary]`, species overrides in `[particles.species.boundary]`, and
periodic axes in `domain.periodic_axes`. See [Boris particle update](BorisPusher.en.html) for advancement,
[particle collision and boundary events](ParticleEvents.en.html) for mesh-versus-box ordering, and
[particle escape and local return](ParticleEscapeReturn.en.html) for open-face treatment.

## Preserve inflow across MPI and restart

The root rank determines global reservoir-inflow and `plane_source` counts and remainders before splitting particles across
ranks. `photo_raycast` `rays_per_batch` is also a global total. Expected inflow and emission therefore do not depend on MPI
world size, although random sequences and individual trajectories may change with world size.

Boundary-inflow remainders are stored per species-face pair in checkpoints and restored on resume. See
[reservoir inflow](ReservoirInjection.en.html) for residual calculation and diagnostics.

## Code reference

- Particle distributions and ray casting: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- Per-source and boundary-inflow batch creation: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- Input validation: [`bem_app_config_parser_validate.f90`](../src/config/app_config_parser/bem_app_config_parser_validate.f90)
- Charge-balance accounting and batch tracking: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- Macro-particle remainder checkpoints: [`bem_restart.f90`](../src/runtime/bem_restart.f90)

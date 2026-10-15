title: Choose where particles enter

Lang: [日本語](ParticleSourcesBoundaries.md) | [English](ParticleSourcesBoundaries.en.md)

# Choose where particles enter

This page helps you choose a particle source from its creation location and physical supply mechanism. `source.mode`
creates particles inside the domain. `[particles.species.inflow]` brings particles across a box boundary from
outside; it has a separate role.

> **Decision rule:** use `volume_seed` for a specified count, `inflow` for an external plasma reservoir,
> and `photo_raycast` for emission from illuminated surfaces.

After reading this page, you should be able to choose one path and continue to the relevant detailed guide or reference.

## Choose from the creation location

| Goal | Setting | What determines particle count | Creation location |
| --- | --- | --- | --- |
| Create a specified initial population or orbit test | `source.mode="volume_seed"` | `sampling.volume_macro_particles_per_batch` | Between `source.pos_low` and `source.pos_high` |
| Bring particles in from plasma outside the box | `inflow` | Reservoir flux × box area × `run.batch.duration_s` | Nonperiodic box face |
| Emit from a surface reached by light | `source.mode="photo_raycast"` | Current density, projected area, `run.batch.duration_s` | First surface hit by a ray |

For boundary inflow alone, specify the physical parameters and `inflow`. Omit `source.mode` and
`sampling.volume_macro_particles_per_batch`; specify them when using volume generation or another source. Omitting them creates no volume
population. See [Input Parameters](Parameters.en.html) for defaults and combination constraints.

## Place a specified count: `volume_seed`

`volume_seed` creates `sampling.volume_macro_particles_per_batch` particles in every batch. Positions are sampled from the box enclosed by
`source.pos_low` and `source.pos_high`; velocities follow the configured drifting Maxwell distribution.

It does not derive count from a physical surface flux. A run can therefore use this source with `run.batch.duration_s=0`,
but no physical duration in seconds is assigned to that batch. The official
[10-minute tutorial](Tutorial.en.html) uses this path to demonstrate inter-batch surface-charging feedback.

## Enter from an external reservoir: `inflow`

`inflow` represents plasma entering through an entire nonperiodic box face. Separately configure the outward
action on every selected face as `open`. Outgoing particles escape at an ordinary open face;
[outer-sheath return](ParticleEscapeReturn.en.html#outer-sheath-barrier) is determined by the Zhao kinetic barrier.

```toml
[[particles.species]]
# density, temperature, charge, mass, and macro-particle weight

[particles.species.inflow]
z_high = "reservoir"
```

See [inflow through a simulation boundary](ReservoirInjection.en.html) for the Maxwell distribution,
per-face flux and remainder, and correction from an upstream plasma potential to the boundary. This model does not solve
trajectories or a self-consistent sheath outside the box. If the case needs an outer sheath response itself, consider
the [Zhao closure](ZhaoStationaryClosure.en.html).

## Emit from illuminated surfaces: `photo_raycast`

`photo_raycast` launches rays through an illumination aperture on a box face and emits particles from the first triangle
hit. See [photoelectron emission and lifecycle](PhotoelectronEmission.en.html) for reaction charge left at the source,
reabsorption, and closed-photoelectron return.

## Tracking is common after creation

Every creation path enters the same particle tracker.

| Outcome | BEACH treatment |
| --- | --- |
| Absorbed by a triangle | Add charge to the hit element's pending change |
| Escapes through an open face | Remove it and count a species-resolved escape |
| Reflected or crosses a periodic boundary | Continue tracking the same particle for the remaining time |
| Still alive at `particles.tracking.max_steps_per_particle` | Count an unresolved survivor and discard it at batch end |

See [particle escape and return](ParticleEscapeReturn.en.html) for outward box boundaries and
[particle collision and boundary events](ParticleEvents.en.html) for event ordering.

## Where to go next

- Combine geometry, source, boundaries, and solver into a case: [Design a case](ConfigurationRecipes.en.html)
- Configure boundary flux and `particles.reservoir.phi_infty_v`: [Inflow through a simulation boundary](ReservoirInjection.en.html)
- Treat charge after a surface collision: [How surfaces charge](SurfaceModels.en.html)
- Search the complete key and constraint set: [Input parameters](Parameters.en.html)

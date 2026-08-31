title: Choose where particles enter

Lang: [日本語](ParticleSourcesBoundaries.md) | [English](ParticleSourcesBoundaries.en.md)

# Choose where particles enter

This page helps you choose a particle source from its creation location and physical supply mechanism. `source_mode`
creates particles inside the domain. `[particles.species.boundary_inflow]` brings particles across a box boundary from
outside; it has a separate role.

> **Decision rule:** use `volume_seed` for a specified count, `plane_source` for continuous supply from an internal
> rectangle, `boundary_inflow` for an external plasma reservoir, and `photo_raycast` for emission from illuminated
> surfaces.

After reading this page, you should be able to choose one path and continue to the relevant detailed guide or reference.

## Choose from the creation location

| Goal | Setting | What determines particle count | Creation location |
| --- | --- | --- | --- |
| Create a specified initial population or orbit test | `source_mode="volume_seed"` | `npcls_per_step` | Between `pos_low` and `pos_high` |
| Supply particles one way from an internal surface | `source_mode="plane_source"` | Flux × area × `batch_duration` | Axis-aligned rectangle |
| Bring particles in from plasma outside the box | `boundary_inflow` | Reservoir flux × box area × `batch_duration` | Nonperiodic box face |
| Emit from a surface reached by light | `source_mode="photo_raycast"` | Current density, projected area, `batch_duration` | First surface hit by a ray |

`source_mode` and `boundary_inflow` are separate settings. In the current schema, a species that uses boundary inflow
alone must also specify `source_mode="volume_seed"` and `npcls_per_step=0`. This does not create particles in a volume.

## Place a specified count: `volume_seed`

`volume_seed` creates `npcls_per_step` particles in every batch. Positions are sampled from the box enclosed by
`pos_low` and `pos_high`; velocities follow the configured drifting Maxwell distribution.

It does not derive count from a physical surface flux. A run can therefore use this source with `batch_duration=0`,
but no physical duration in seconds is assigned to that batch. The official
[10-minute tutorial](Tutorial.en.html) uses this path to demonstrate inter-batch surface-charging feedback.

## Supply particles from an internal surface: `plane_source`

`plane_source` supplies particles along `source_normal` from a rectangle inside the box. `pos_low` and `pos_high`
must share one coordinate and span a positive area in the other two directions.

```toml
source_mode = "plane_source"
pos_low = [0.2, 0.2, 0.7]
pos_high = [0.8, 0.8, 0.7]
source_normal = [0.0, 0.0, -1.0]
```

Flux, rectangle area, and `batch_duration` determine the particle count. Because the rectangle is not an outer plasma
boundary, it does not use the potential correction from infinity defined by `reservoir.phi_infty`.

## Enter from an external reservoir: `boundary_inflow`

`boundary_inflow` represents plasma entering through an entire nonperiodic box face. Separately configure the outward
action on every selected face as `open`.

```toml
[[particles.species]]
source_mode = "volume_seed"
npcls_per_step = 0
# density, temperature, charge, mass, and macro-particle weight

[particles.species.boundary_inflow]
z_high = "reservoir"
```

See [inflow through a simulation boundary](ReservoirInjection.en.html) for Maxwell or velocity-grid distributions,
per-face flux and remainder, and correction from an upstream plasma potential to the boundary. This model does not solve
trajectories or a self-consistent sheath outside the box. If the case needs an outer sheath response itself, consider
[quasistatic matching-plane coupling](MatchingPlaneCoupling.en.html).

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
| Still alive at `sim.max_step` | Count an unresolved survivor and discard it at batch end |

See [particle escape and return](ParticleEscapeReturn.en.html) for outward box boundaries and
[particle collision and boundary events](ParticleEvents.en.html) for event ordering.

## When reading legacy `reservoir_face`

`source_mode="reservoir_face"` is a deprecated compatibility input. Use `boundary_inflow` for a new external-plasma
case and `plane_source` for an internal rectangle. Moving a legacy aperture to a complete box face changes area and total
inflow, so BEACH does not convert it silently.

The [input parameter reference](Parameters.en.html) is the canonical list of keys, combinations, MPI behavior, and
checkpoint behavior.

## Where to go next

- Combine geometry, source, boundaries, and solver into a case: [Design a case](ConfigurationRecipes.en.html)
- Configure boundary flux and `phi_infty`: [Inflow through a simulation boundary](ReservoirInjection.en.html)
- Treat charge after a surface collision: [How surfaces charge](SurfaceModels.en.html)
- Search the complete key and constraint set: [Input parameters](Parameters.en.html)

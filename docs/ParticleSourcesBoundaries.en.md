title: Particle sources and boundaries

Lang: [English](ParticleSourcesBoundaries.en.md) | [日本語](ParticleSourcesBoundaries.md)

# Particle sources and boundaries

This page organizes how particles enter the simulation and what happens at box boundaries. See
[Outer plasma models](OuterPlasmaModels.en.html) for the outer-domain Poisson models themselves.

## Particle sources

| `source_mode` | Purpose | Position |
| --- | --- | --- |
| `volume_seed` | Initial particles and simple orbit tests | Between `pos_low` and `pos_high` |
| `reservoir_face` | Inflow from a reservoir outside the box | Selected box face |
| `photo_raycast` | Emission from illuminated surfaces | First surface hit by a ray |

### Volume seed

`npcls_per_step` particles are generated in the selected region with a shifted Maxwellian velocity. This source
does not derive particle count from a physical flux; use `reservoir_face` for quantitative steady inflow.

### Reservoir face

Expected macro-particle count follows from area $A$, inflow flux $\Gamma_\mathrm{in}$, `batch_duration`, and
macro-particle weight $w$.

$$
N_\mathrm{macro,expected}=\frac{\Gamma_\mathrm{in}A\,\mathrm{batch\_duration}}{w}
$$

The fractional count is carried to the next batch as a species residual. Under MPI, the global count and residual
are determined once before rank distribution, so expected global inflow does not depend on world size.

Velocity is sampled from a flux-weighted Maxwellian or a velocity grid. If a potential barrier exists between
the upstream reservoir and the face, accessible normal velocities are selected and mapped by energy conservation.

### Photo raycast

Rays start at an injection face and emit a particle from the first hit element. Ray count, projected area,
current density, and `batch_duration` determine macro-particle weight. See
[Surface charging models](SurfaceModels.en.html#photoemission-ledger) for source-element charge bookkeeping.

## Box boundaries

| Boundary | Particle behavior |
| --- | --- |
| `open` | Transfer to an outer model or remove as escape |
| `reflect` | Reverse normal velocity and integrate the remaining time |
| `periodic` | Wrap to the opposite face and integrate the remaining time |

When a mesh hit and box crossing occur in one step, the earlier trajectory event wins. See
[Particle tracking and collision](ParticleTrackingCollision.en.html).

## Open-boundary handling

- `open_boundary_model="escape"`: count escape immediately.
- Legacy potential barrier: decide escape or reflection from boundary and infinity potentials.
- 1D instant return: decide escape or return from conserved energy in an outer profile.
- 3D explicit orbit: track motion in the outer three-dimensional field.

Unsupported combinations stop during validation to prevent applying the same potential shift or cutoff twice.

## Selection guide

- Small deterministic test: `volume_seed`
- Steady plasma inflow: `reservoir_face`
- Surface photoemission: `photo_raycast`
- Simple finite box: `open_boundary_model="escape"`
- Self-consistent outer sheath: pair the matching outer and return models

See [Outer plasma models](OuterPlasmaModels.en.html) and
[Configuration parameters](Parameters.en.html). Detailed formulas remain in the legacy combined
[Outer sheath and reservoir particle boundaries](SheathReservoirBoundary.en.html).

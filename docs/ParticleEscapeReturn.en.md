title: Particle escape and local return

Lang: [日本語](ParticleEscapeReturn.md) | [English](ParticleEscapeReturn.en.md)

# Particle escape and local return

BEACH does not evolve plasma outside the box. An open face locally applies
`external_boundary.ordinary_open.model="escape" | "potential_barrier"`. The treatment is shared by `reservoir_face`,
`photo_raycast`, and `volume_seed`. Closed PE alone combines a species-specific z-high reflection with a charge closure.

## 1. `escape`: remove a particle at an open face

```toml
[external_boundary.field]
model = "none"

[external_boundary.particles]
mode = "local_source"

[external_boundary.ordinary_open]
model = "escape"
```

The particle is removed at the crossing. Its macro charge $qw$ is recorded in the species-resolved
`escaped_to_infinity` ledger term, while surface charge `q_elem` is unchanged.

See [Particle collision and boundary events](ParticleEvents.en.html) for simultaneous face crossings and
reintegration of the step remainder after reflecting or periodic faces.

## 2. `potential_barrier`: decide reflection at a scalar barrier

`potential_barrier` is a reduced model that uses only scalar potential at the crossing.

```toml
[sim]
phi_infty = 0.0

[external_boundary.field]
model = "none"

[external_boundary.particles]
mode = "local_source"

[external_boundary.ordinary_open]
model = "potential_barrier"
```

For crossing-point potential $\phi_b$ and outward normal speed $v_n>0$, the barrier to infinity is

$$
U_b=q(\phi_\infty-\phi_b).
$$

If

$$
\frac12m v_n^2<U_b\quad\text{and}\quad U_b>0,
$$

the normal velocity is reversed and the step remainder is tracked. Otherwise the particle escapes. Tangential velocity
is unchanged.

It does not retain an electric field outside the open face, a turning position, flight time, or space charge. A corner that
crosses multiple open faces stops as `ambiguous_open_corner`.

The crossing potential follows the batch-start fixed field used for particle motion and includes the local potential of
`sim.e0`. Because a uniform field has no finite potential at infinity, a configuration that combines `sim.e0` with this
model must use a consistent effective reservoir reference for `phi_infty`.

## Return of closed photoelectrons

To keep photoelectrons inside the box, set `z_high_boundary="reflect"` and
`surface_charge_closure="neutral_return"` on the corresponding negative `photo_raycast` species. A reflected
photoelectron remains the same tracked particle, and the emission countercharge is corrected at the end of the batch.
See [Photoelectron emission](PhotoelectronEmission.en.html) for the configuration and charge-balance requirements.

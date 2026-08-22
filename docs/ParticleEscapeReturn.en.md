title: Particle escape and local return

Lang: [日本語](ParticleEscapeReturn.md) | [English](ParticleEscapeReturn.en.md)

# Particle escape and local return

BEACH does not evolve plasma outside the box. `[particle_boundary]` makes each nonperiodic face `open`, `reflect`, or
`redistributed_reflect`, and open
faces apply `ordinary_open_model="escape" | "potential_barrier"`. The treatment is shared by boundary-reservoir inflow,
`plane_source`, `photo_raycast`, `volume_seed`, and deprecated `reservoir_face`. Periodicity belongs to
`domain.periodic_axes`; particle-boundary tables cannot specify
`periodic`.
The fully qualified key is `particle_boundary.ordinary_open_model`.

## 1. `escape`: remove a particle at an open face

```toml
[particle_boundary]
x_low = "open"
x_high = "open"
y_low = "open"
y_high = "open"
z_low = "open"
z_high = "open"
ordinary_open_model = "escape"
```

The particle is removed at the crossing. Its macro charge $qw$ is recorded in the species-resolved
`escaped_to_infinity` ledger term, while surface charge `q_elem` is unchanged.

See [Particle collision and boundary events](ParticleEvents.en.html) for simultaneous face crossings and
reintegration of the step remainder after reflecting or periodic faces.

## 2. `potential_barrier`: decide reflection at a scalar barrier

`potential_barrier` is a reduced model that uses only scalar potential at the crossing.

```toml
[particle_boundary]
ordinary_open_model = "potential_barrier"

[reservoir]
phi_infty = 0.0
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
model must use a consistent effective reservoir reference for `reservoir.phi_infty`.

## Outer barrier at a matching plane

`matching_plane_quasistatic` applies electron and ion access potentials plus the maximum outer PE barrier returned by
the response table at z-high, independently of the ordinary-open model. The outward decision uses the same scalar-energy
test, but its reference values update in every fixed-point iteration. A reflected PE contributes to return, a transmitted
PE contributes to escape, and BEACH compares their sum with the measured outward crossings.

Return reverses normal velocity immediately at the same z-high position. This quasistatic approximation retains no outer
turning position, flight time, or lateral displacement. Do not stack a manual `potential_barrier` or closed-PE reflection
on a matching-plane case. See [Quasistatic Matching-Plane Coupling](MatchingPlaneCoupling.en.html) for response columns and scope.

## Return of closed photoelectrons

For closed PE, set the face matching `inject_face` to `reflect` or `redistributed_reflect` in
`[particles.species.boundary]` for the negative `photo_raycast` species, and set
`surface_charge_closure="neutral_return"`. Both actions reverse normal velocity and preserve tangential velocity, but they
treat the return position differently.

| Action | Return position |
| --- | --- |
| `reflect` | Preserve the tangential position of the boundary event |
| `redistributed_reflect` | For one face, uniformly redistribute both in-plane coordinates over the box span excluding its end guards |

At a simultaneous edge or corner event, `redistributed_reflect` relocates only axes outside the event mask. See
[Particle collision and boundary events](ParticleEvents.en.html) for details. Keeping the global face open lets ambient species
with `inherit` follow the ordinary-open contract. See [Photoelectron emission](PhotoelectronEmission.en.html) for the
configuration and charge-balance requirements.

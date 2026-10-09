title: Particles at box boundaries

Lang: [English](ParticleEscapeReturn.en.md) | [日本語](ParticleEscapeReturn.md)

# Particles at box boundaries

When a particle reaches a box face, it wraps to the opposite side on a periodic face and follows the treatment chosen for the
face on a non-periodic face. The treatment is the same regardless of the source that created the particle. BEACH does not solve
the plasma outside the box, so the treatment chosen here decides in its place whether a particle that leaves the face comes back.

## Choose the treatment of each face

| Treatment | Setting | What happens to the particle |
|---|---|---|
| Periodic | `domain.periodic_axes` | Moves to the opposite face with unchanged velocity |
| Open | `"open"` | Escapes. With `open_model="potential_barrier"`, a potential barrier decides between reflection and escape |
| Reflect | `"reflect"` | Reverses the normal velocity. Position and tangential velocity are unchanged |
| Reflect with redistribution | `"redistributed_reflect"` | Reverses the normal velocity and chooses the in-plane position again uniformly |

```toml
[particles.boundary]
z_low = "open"
z_high = "open"
open_model = "escape"      # treatment of open faces; the default
```

Periodic axes are shared by the field and the particles, so neither `[particles.boundary]` nor per-species settings can change
them. Omitted non-periodic faces are open. To change the treatment for one species, choose `inherit` (default), `open`,
`reflect`, or `redistributed_reflect` in `[particles.species.boundary]`.

## Escape

A particle that crosses an open face is removed at the crossing point. Its charge is counted as escaped charge of its species
(`escaped_to_infinity_C` in `charge_ledger.csv`), and the surface charge does not change.

## Potential barrier (`potential_barrier`)

Set the upstream potential $\phi_\infty$ of the external plasma (`particles.reservoir.phi_infty_v`), and decide whether a
particle leaving an open face can reach upstream from the potential $\phi_b$ at its crossing point.

$$
\Delta U=q(\phi_\infty-\phi_b),\qquad
K_n=\frac12 m v_n^2
$$

$v_n>0$ is the outward normal velocity. If $\Delta U>0$ and $K_n<\Delta U$, the normal velocity is reversed and tracking
continues; otherwise the particle escapes. The tangential velocity is unchanged.

```toml
[particles.boundary]
open_model = "potential_barrier"

[particles.reservoir]
phi_infty_v = 0.0
```

$\phi_b$ is evaluated with the field frozen at the start of the batch, including any external field. A uniform external field
has no potential at infinity, so when you use one, set $\phi_\infty$ consistently with that reference. At a corner where several
open faces are crossed at once the test is undefined and the run stops.

To filter incoming particles with the same upstream potential, set the inflow mapping to `infinity_barrier`
([Inject through a boundary](ReservoirInjection.en.html#3-choose-the-inflow-mapping)).

## Outer-sheath barrier

With the outer-sheath connection, electrons and photoelectrons leaving through the top face are tested against the barrier given
by the zero-current root of the outer sheath. The test has the same form as the potential barrier; the barrier potential and the
choice of the return position differ ([Connecting to the outer sheath](ZhaoStationaryClosure.en.html#particles-entering-and-leaving-through-the-top-face)).

## Reflect and reflect with redistribution

Both reverse only the normal velocity and keep the tangential velocity. `reflect` keeps the position. For reflection at a single
face, `redistributed_reflect` chooses both in-plane coordinates again uniformly within the face span excluding small margins at
the edges. When a particle reaches several faces at once at a corner or edge, only the axes not belonging to the faces reached
are chosen again.

Reflection of photoelectrons at the top face is used by closed photoelectrons
([Photoelectron emission and charge closures](PhotoelectronEmission.en.html#closed-photoelectrons-neutral_return)).
The order of decisions when several faces are reached at once is in [Collision and boundary events](ParticleEvents.en.html).

## Outputs to check

| Output | What to check |
|---|---|
| `summary.txt` | `escaped_boundary` (particles that left through a boundary) and `particle_ordinary_open_model` (treatment used for open faces) |
| `charge_ledger.csv` | `escaped_count` and `escaped_to_infinity_C` per species |

`escaped` in `summary.txt` also includes particles whose fate was not decided within the step limit (`survived_max_step`).
Read the number that actually left through a boundary from `escaped_boundary`.

## Scope

- The field outside an open face, turning positions, flight times, and space charge are not solved. The potential barrier is an
  approximation that decides only from the potential difference to upstream.
- Reflection places a mirror on a box face. It does not represent an external plasma or sheath.

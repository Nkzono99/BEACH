title: Surface charging models

Lang: [English](SurfaceModels.en.md) | [日本語](SurfaceModels.md)

# Surface charging models

This page describes what happens after a particle hits a triangular surface and how element charge is carried
into the next batch. The primary model accumulates absorbed-particle charge on an insulating surface.

## Common interaction

The default particle-surface interaction is absorption. When particle `p` first hits element `i`, tracking ends
and its weighted charge is added to a batch-local delta.

$$
\Delta q_i \mathrel{+}=q_p w_p
$$

`delta_q_elem` is committed to `q_elem` after all particles in the batch are processed. Charge deposited in one
batch therefore contributes to the field starting in the next batch.

## Surface models

| `surface_model` | Current behavior | Note |
| --- | --- | --- |
| `insulator` | Retain charge on the hit element | Primary scope |
| `conductor` | Preserve object charge while redistributing it toward equal element potential | Only with `field_bc_mode="free"` |
| `dielectric` | Retain charge and output `epsilon_r` as metadata | Polarization is not implemented |

### Insulator

`q_elem` stores total element charge in coulombs. With `field.element_kernel="triangle_p0"`, it remains a total
charge and is divided by area only during field evaluation. Lateral conduction and leakage are not included.

### Floating conductor

Conductor elements with one `mesh_id` form a floating object. After charge commit, BEACH solves for an element
charge redistribution that preserves object charge while approximating equal potential. This is not a grounded
fixed-potential boundary. Check the support matrix in [Configuration parameters](Parameters.en.html).

### Dielectric metadata

`dielectric` and `epsilon_r` currently preserve material metadata only. BEACH does not impose dielectric
interface conditions or solve self-consistent polarization charge.

## Photoemission ledger

With `photo_raycast` and `deposit_opposite_charge_on_emit=true`, charge opposite to the emitted particle is added
to the source element as `photo_emission_dq`. It is tracked separately from later collision deposition and merged
at batch commit.

## Outputs to inspect

- `charges.csv`: committed element charge
- `charge_history.csv`: charge by batch
- `summary.txt`: absorption and charge ledger
- `mesh_sources.csv`: `mesh_id` and surface model

Check charge balance and mesh refinement with [Validate simulation results](ValidationGuide.en.html). The legacy
combined implementation detail remains in
[Particle tracking and surface charge accumulation](ParticleChargeLoop.en.html#11-charge-deposition-and-surface-models).

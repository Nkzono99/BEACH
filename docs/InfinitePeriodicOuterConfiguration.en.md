title: Infinite-periodic periodic2 with outer plasma

Lang: [日本語](InfinitePeriodicOuterConfiguration.md) | [English](InfinitePeriodicOuterConfiguration.en.md)

# Infinite-periodic periodic2 with outer plasma

This configuration assembles an x/y infinite-periodic surface field and a z-direction outer-plasma closure in one electrostatic
snapshot. Near images, an Ewald-generated far `k\ne0` operator, physical `k=0`, and outer response are composed without duplication,
and the same outer potential controls reservoir inflow and particle return.

The standard composition is split `kinetic_1d`, which solves a nonlinear mean sheath from species VDFs. The unified composition
is an advanced option limited to rough surfaces that have no split window and require linear screening.

## Assign one owner to each field component

| Component | Owner |
| --- | --- |
| Primary and near images | FMM or spectral base field |
| Infinite-periodic far `k\ne0` | `cached_kneq0`, or `panel_spectral_reference` for small validation |
| Surface `k=0` | Periodic zero-mode plan/state |
| Outer mean response | `kinetic_1d` or unified zero-mode solve |
| Outer nonzero response | Reflection/transmission correction of the unified model |
| Reservoir velocity map | Outer interface potential difference |
| Outward escape/return | The same outer profile or unified 3-D field |

See [periodic2 electrostatics](PeriodicElectrostatics.en.html) for periodic operators and `k=0` equations.

## Select standard split kinetic or advanced linear unified response

| Configuration | Mean plasma | Nonzero mode | Particle transfer |
| --- | --- | --- | --- |
| **Standard: split kinetic** | Nonlinear `kinetic_1d` VDF closure | Assumed decayed by the interface | `particles.mode="local_source"`, or the 1-D profile with `"same_batch"` |
| Advanced: unified linear | Linear Poisson with accessible fraction | Joined to screened modes at response start | `particles.mode="local_source"`, or a 3-D orbit with `"same_batch"` |

Split kinetic represents species VDFs, Bohm entry, and mean photoelectron density, but assumes a split window containing local
surface field below the interface. See [Kinetic 1-D outer plasma](KineticOuterPlasma.en.html).

Unified linear starts plasma response across the roughness range but uses a linear Debye closure and does not solve species VDFs or
Bohm entry. See [Unified linear response](UnifiedLinearResponse.en.html).

## Share one batch-start snapshot between inflow and return

1. Refresh FMM/source multipoles and surface zero mode from committed `q_elem`.
2. Apply the cached far operator to the current root multipole.
3. In split mode, solve the outer profile on selected strides from interface field; in unified mode, solve every snapshot refresh.
4. Update potential gauge, Gauss residual, and interface or linearity diagnostics.
5. Use outer state to determine global z-high reservoir counts and interface velocities.
6. Perform photo ray casting and record source reaction charge in the batch difference.
7. Track all particles through the immutable snapshot.
8. Map z-high outward crossings to escape or return with the same outer state.
9. MPI all-reduce surface absorption and emission differences and commit at batch end.

The operator and outer Poisson problem are not rebuilt after each hit. Sources and return share one batch-start outer state, while
surface charge changes the field starting in the next batch.

## Use one energy equation for kinetic inflow and return

The first negative and positive z-high `reservoir_face` species define electron and ion VDFs at infinity. The converged
$\phi_I-\phi_\infty$:

- selects accessible $v_{n,\infty}$ and accelerates or decelerates them to $v_{n,I}$ on inflow;
- derives $v_{n,\infty}^2$ from outward $v_{n,I}$ and classifies escape or turning return.

For tracked split kinetic, `external_boundary.particles.inflow_model="auto"`
delegates inflow to the same profile. Do not combine it with
`inflow_model="infinity_barrier"`. See
[`reservoir_face` inflow and velocity sampling](ReservoirInjection.en.html) and [Particle escape and return](ParticleEscapeReturn.en.html).

## Select mean outer density separately from tracked photoelectrons

Photoelectron outer density and tracked return are separate choices.

| Choice | Outer density | Tracked particle |
| --- | --- | --- |
| `photoelectron_density_model="none"` | No mean photoelectron density | Ordinary source and orbit only |
| `photoelectron_density_model="kinetic_mean"` with `particles.mode="local_source"` | Stationary mean outgoing/returning density | Ordinary open handling at z-high |
| `photoelectron_density_model="kinetic_mean"` with `particles.mode="same_batch"` | Same mean density | Individual interface crossings also return or escape through the profile |
| Unified with `particles.mode="same_batch"` | No mean photoelectron closure | Individual 3-D outer trajectories |

Tracked return requires `deposit_opposite_charge_on_emit=true`. The mean density model neither
replaces tracked surface deposition nor adds a statistical return deposit. See
[Photoelectron emission and lifecycle](PhotoelectronEmission.en.html).

## Use the cached nonzero operator in production

FMM production makes ownership explicit:

```toml
[sim]
field_solver = "fmm"
field_bc_mode = "periodic2"
field_periodic_image_layers = 1
field_periodic_far_correction = "cached_kneq0"
field_periodic_ewald_layers = 4

[field]
element_kernel = "triangle_p0"

[periodic2]
nonzero_mode_backend = "cached_kneq0"
zero_mode_policy = "exclude_k0"
lower_boundary_model = "symmetric_vacuum"
```

This block alone does not enable outer plasma. Add kinetic or unified
`[external_boundary.field]` and `[external_boundary.particles]` configuration. See the
minimal cache example [`examples/periodic2_cached_panel.toml`](../examples/periodic2_cached_panel.toml).

Small references use Direct plus `panel_spectral_reference`.
[`examples/periodic2_kinetic_outer.toml`](../examples/periodic2_kinetic_outer.toml) is a small contract fixture for the standard
kinetic composition. [`examples/periodic2_zhao_transient_outer.toml`](../examples/periodic2_zhao_transient_outer.toml) is an
expected-fail fixture in which the transient Zhao queue's physical-timescale guard rejects a long flight.
[`examples/periodic2_unified_linear_response.toml`](../examples/periodic2_unified_linear_response.toml)
checks the analytic limit and applicability of advanced linear screening. Neither selects a large production backend by itself.

## Apply each closure exactly once

- Symmetric `k=0` inside `cached_kneq0` versus physical `k=0` in the snapshot.
- Kinetic interface-potential map versus finite-image `infinity_barrier`.
- Profile return versus finite-image open `potential_barrier` at z-high.
- `photoelectron_density_model="kinetic_mean"` outer density versus a fictitious statistical return deposit on the surface.
- Unified base nonzero field versus the incident mode after reflection/transmission replacement.

Validation fails closed on major unsupported combinations, but numerical convergence and physical applicability still require user
verification.

## Hold the batch-start field fixed during outer flight

Both 1-D instant return and 3-D explicit orbit calculate outer flight time but do not add it to global simulation time. Flight time
relative to `field_evolution_timescale` must remain below `max_frozen_field_ratio`. The matching `kinetic_1d` +
`zhao_charge_driven` configuration can instead retain 1-D events across batches and count return/escape when each event becomes
due. Its outer domain ends at $L=10\lambda_{D,pe}$: a particle reaching $L$ escapes into the reservoir, and no Robin tail beyond
$L$ is used to classify return. The terminal state is resolved with the field at enqueue and is not reintegrated after that
field changes. Queue mode applies the frozen-field bound to each event's `tau_outer`, its delay to the next batch-start poll,
and the half-batch bound on midpoint crossing-time uncertainty. Configuration validation applies the same bound to
`batch_duration`. Persistent queuing remains unavailable for a
3-D explicit orbit.
See
[Particle escape and return](ParticleEscapeReturn.en.html#queue-outer-flight-for-the-transient-zhao-closure).

## Check convergence and balance component by component

| Target | Parameters to vary | Quantities to compare |
| --- | --- | --- |
| Near/far periodic field | Image layer, Ewald layer, cache cold/warm | $\phi,\mathbf E$, force, operator residual |
| Physical `k=0` | Lower closure and height/grid refinement | Gauss residual and interface field |
| Kinetic profile | Debye length, source sampling, outer stride | $\phi_I$, currents, nonlinear residual |
| Unified profile | `unified_grid_points`, height sampling, mode layer | Linearity, accessible fraction, Gauss residual |
| Reservoir | Macro target and batch duration | Inflow current and macro residual |
| Photoelectron | Ray count, `dt`, outer return | Emission, reabsorption, and escape/return charge |
| Transient Zhao queue | `batch_duration`, time-scale guard, ray count, area, interface location | $\eta$, column residual, return/escape current, force |
| Batch coupling | `batch_duration` | Steady surface charge and current balance |

Inspect `summary.txt`, `outer_plasma_profile.csv`, `charge_ledger.csv`, and `charges.csv` together. Checking only field, particle
count, or surface charge can miss component duplication or charge loss.

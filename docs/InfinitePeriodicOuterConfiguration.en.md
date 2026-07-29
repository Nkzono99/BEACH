title: Infinite-periodic periodic2 with outer plasma

Lang: [日本語](InfinitePeriodicOuterConfiguration.md) | [English](InfinitePeriodicOuterConfiguration.en.md)

# Infinite-periodic periodic2 with outer plasma

This configuration assembles an x/y infinite-periodic surface field and a z-direction outer-plasma closure in one electrostatic
snapshot. Near images, an Ewald-generated far `k\ne0` operator, physical `k=0`, and outer response are composed without duplication,
and the same outer potential controls reservoir inflow and particle return.

The supported outer-sheath composition is split `kinetic_1d`, which constructs a mean sheath from species VDFs or the
ambient linear response.

## Assign one owner to each field component

| Component | Owner |
| --- | --- |
| Primary and near images | FMM or spectral base field |
| Infinite-periodic far `k\ne0` | `cached_kneq0`, or `panel_spectral_reference` for small validation |
| Surface `k=0` | Periodic zero-mode plan/state |
| Outer mean response | `kinetic_1d` |
| Reservoir velocity map | Outer interface potential difference |
| Outward escape/return | The same outer profile |

See [periodic2 electrostatics](PeriodicElectrostatics.en.html) for periodic operators and `k=0` equations.

## Inspect the split-kinetic responsibilities

| Configuration | Mean plasma | Nonzero mode | Particle transfer |
| --- | --- | --- | --- |
| Split kinetic | Nonlinear `kinetic_1d` VDF closure | Assumed decayed by the interface | `particles.mode="local_source"`, or the 1-D profile with `"same_batch"` |

Split kinetic represents species VDFs, Bohm entry, and mean photoelectron density, but assumes a split window containing local
surface field below the interface. See [Kinetic 1-D outer plasma](KineticOuterPlasma.en.html).

## Share one batch-start snapshot between inflow and return

1. Refresh FMM/source multipoles and surface zero mode from committed `q_elem`.
2. Apply the cached far operator to the current root multipole.
3. Solve the outer profile on selected strides from interface field.
4. Update potential gauge, Gauss residual, and interface diagnostics.
5. Use outer state to determine global z-high reservoir counts and interface velocities.
6. Perform photo ray casting and record source reaction charge in the batch difference.
7. Track all particles through the immutable snapshot.
8. Map z-high outward crossings to escape or return with the same outer state.
9. MPI all-reduce surface absorption and emission differences and commit at batch end.

The operator and outer Poisson problem are not rebuilt after each hit. On the ordinary explicit path, sources and return share
one batch-start outer state, while surface charge changes the field starting in the next batch. The implicit mean path for
`ambient_linear_debye`, described below, instead retains batch-start $k\ne0$ after the initial particle trace while updating
$k=0$ and the outer state with a continuous-Maxwellian mean solve, then uses energy-resolved rays once for local redistribution.

## Use one energy equation for kinetic inflow and return

The first negative and positive z-high `reservoir_face` species define electron and ion VDFs at infinity. The converged
$\phi_I-\phi_\infty$:

- selects accessible $v_{n,\infty}$ and accelerates or decelerates them to $v_{n,I}$ on inflow;
- derives $v_{n,\infty}^2$ from outward $v_{n,I}$ and classifies escape or turning return.

For tracked split kinetic, `external_boundary.particles.inflow_model="auto"`
delegates inflow to the same profile. Do not combine it with
`inflow_model="infinity_barrier"`. Under `implicit_mean`, photoelectron interface crossings are held until the mean-field
update and then passed to the same individual outflow map. See
[`reservoir_face` inflow and velocity sampling](ReservoirInjection.en.html) and [Particle escape and return](ParticleEscapeReturn.en.html).

## Select mean outer density separately from tracked photoelectrons

Photoelectron outer density and tracked return are separate choices.

| Choice | Outer density | Tracked particle |
| --- | --- | --- |
| `none` (`ambient_linear_debye`: internally resolved, no public key) | No mean photoelectron density | Ordinary source and orbit; trace to the interface with `ambient_linear_debye + same_batch` |
| `photoelectron_density_model="kinetic_mean"` with `particles.mode="local_source"` | Stationary mean outgoing/returning density | Ordinary open handling at z-high |
| `photoelectron_density_model="kinetic_mean"` with `particles.mode="same_batch"` | Same mean density | Individual interface crossings also return or escape through the profile |

Tracked return requires `deposit_opposite_charge_on_emit=true`. The mean density model neither
replaces tracked surface deposition nor adds a statistical return deposit. See
[Photoelectron emission and lifecycle](PhotoelectronEmission.en.html).

## Advance local $k\ne0$ and mean $k=0$ separately for the ambient linear response

The combination of `kinetic_closure="ambient_linear_debye"`, `particles.mode="same_batch"`, and an enabled negative
`photo_raycast` species automatically resolves to internal `implicit_mean`. No update mode is added to public configuration.

1. Trace existing 3-D particles with the batch-start field and deposit local emission, reabsorption before the interface,
   and ambient absorption on elements, while deferring each photoelectron z-high crossing with its full macro weight and
   source provenance.
2. Measure every surface-charge delta staged by the first trace, separate the ambient electron/ion and PE outward
   components, and retain the remainder as `J_other`. Use an analytic-Maxwellian backward-Euler solver to determine mean
   total charge, $\phi_I$, and continuous escape/return fractions.
3. Apply the analytic return fraction to each crossing's source countercharge, temporarily neutralize it at the source,
   and refresh only $k=0$ and the outer profile from the resulting element distribution with the mean-solved total charge.
4. Pass each full-weight crossing through the common kinetic 1-D profile mapper. Classify return or escape from normal
   energy; include outer flight and lateral displacement for return, and remap a later z-high recrossing through the same driver.
5. Normalize analytic return charge with one scale over the source legs (temporary neutralization at emission elements) and
   actual-hit destination legs of rays that are ultimately reabsorbed, forming a zero-sum transaction.
6. Commit the normalized actual deposits exactly once and retain the updated $k=0$ and outer profile for the next batch.

The $k\ne0$ operator remains at its batch-start value. Mean total charge and retarding potential update within the current
batch, while sampled local redistribution changes $k\ne0$ only after commit in the next batch. Not feeding discrete ray
classification back into the mean solve avoids a return/escape two-cycle near the separatrix. A nonzero transaction residual
is an error. `emit_current_density_a_m2` determines tracked 3-D emission weights, while PE mean-source
amplitude is measured from actual interface crossings. `J_other_A_m2` reports the remaining measured surface current,
including extra species and other-face or unresolved outcomes. The analytic-Maxwellian closure determines mean return/escape
charge totals; energy-resolved rays sample orbit and landing distributions. The `BEACH implicit-mean` record reports
`transaction_residual_C`, `mean_solver_iterations`, `sample_escape_fraction`, and `return_weight_scale`.

This path requires `outer_update_stride=1`, positive `batch_duration`, and
`deposit_opposite_charge_on_emit=true`, together with exactly one negative `photo_raycast` species and
`photo_raycast.normal_drift_speed=0`. The analytic Maxwellian escape fraction does not include normal emission drift.
It adds no
mesh-mode-specific requirement. Omit
`photoelectron_density_model` from public TOML; it resolves internally to `none`. No public update-mode or return-kernel
setting is added. The equations and diagnostics are documented in
[Kinetic 1-D outer plasma](KineticOuterPlasma.en.html#separate-the-ambient-linear-debye-response-from-tracked-photoelectrons).
The run fails closed if the mean solver does not converge, analytic return is positive but no reabsorbed sample is available,
transaction charge does not balance, or a deferred ray does not terminate as absorbed or escaped within its allowed trace.
Deferred PE rays are quasistationary shadows: outer flight time and frozen-field ratio remain diagnostic, but an over-limit
ratio is not a stopping condition. This ambient linear path does not support nonlinear photoelectron density or an
outer-cloud inventory.

A UV-only configuration without ambient electron and ion reservoirs does not satisfy this infinity closure. With
`field.model="none"` and z-high escape, it is a finite-box transient control, not a quasineutral stationary outer sheath.

## Couple a measured energy CDF to a nonmonotone Zhao profile under strong PE

The combination of `kinetic_closure="zhao_charge_driven"`, `particles.mode="same_batch"`, a negative
`photo_raycast` species, and `steady_start_mode="zhao_floating"` also uses internal `implicit_mean`; it adds no public
update-mode key. Instead of the ambient linear closure's Maxwellian return fraction, each rank contributes normal kinetic
energy and positive macro-charge magnitude for its z-high crossings to rank 0. BEACH combines the charge-weighted empirical
CDF with the full-path Zhao barrier $B(Q)$ and solves

$$
Q=Q_{\mathrm{base}}+Q_{\mathrm{escape}}\!\left[B(Q)\right].
$$

Macro particles with identical energy form one group. Only the group intersected by the separatrix can receive a fractional
escape weight. The solved weights are broadcast to all ranks and applied as one zero-sum transaction over each crossing's
source leg and the actual post-return hit destination. `emit_current_density_a_m2` controls tracked 3-D surface emission,
whereas the Zhao mean-source scale is resolved separately from measured interface-crossing current.

Every candidate profile is restricted to one connected path from the previously committed root,

$$
G_b\!\left(y;E_I(\lambda),n_{\mathrm{pe0}}(\lambda)\right)=0,
$$

followed by pseudo-arclength continuation. After constructing the common source-and-field anchor, the fixed-source slice
adaptively checks $dB/dQ\ge0$ from its tangent. A vanishing lambda tangent at a fold, fixed-parameter Jacobian rank loss,
barrier decrease, or branch/domain endpoint stops without fallback. Type A/B requires $E_I>0$, Type C requires $E_I<0$,
and a lambda-equals-one event corrector lands the target instead of an ordinary fixed-field Newton solve. Paths from the
common root to $Q_0$ and $Q_M$ check the complete charge interval before a binary search locates the first true value of a
monotone order predicate. Pure-root selection costs $O(\log M)$ connected candidate solves. This local continuation guard
is finite precision, not a global interval-arithmetic proof. Final acceptance also recomputes the measured-CDF charge
residual and stops on an inconsistent order-predicate or marginal bracket.

The $k\ne0$ operator remains at its batch-start value during the first 3-D trace; mean charge and the Zhao outer state update
within the same batch, while sampled local redistribution enters $k\ne0$ only after commit. A deferred ray may cross outward
once and return at most once. Recrossing or disagreement between the solved weight and terminal outcome fails closed. This
quasistationary path does not resolve outer-flight delay or retain an outer inventory between batches; use `zhao_queue` when
those effects are required.

## Use the cached nonzero operator in production

FMM production makes ownership explicit:

```toml
[sim]
field_solver = "fmm"
field_bc_mode = "periodic2"
field_periodic_image_layers = 1
field_periodic_far_correction = "cached_kneq0"
field_periodic_ewald_layers = 4

[periodic2]
nonzero_mode_backend = "cached_kneq0"
zero_mode_policy = "exclude_k0"
lower_boundary_model = "symmetric_vacuum"
```

This block alone does not enable outer plasma. Add kinetic
`[external_boundary.field]` and `[external_boundary.particles]` configuration. See the
minimal cache example [`examples/periodic2_cached_panel.toml`](../examples/periodic2_cached_panel.toml).

Small references use Direct plus `panel_spectral_reference`.
[`examples/periodic2_kinetic_outer.toml`](../examples/periodic2_kinetic_outer.toml) is a small contract fixture for the standard
kinetic composition. [`examples/periodic2_zhao_transient_outer.toml`](../examples/periodic2_zhao_transient_outer.toml) is an
expected-fail fixture in which the transient Zhao queue's physical-timescale guard rejects a long flight.
Neither selects a large production backend by itself.

## Apply each closure exactly once

- Symmetric `k=0` inside `cached_kneq0` versus physical `k=0` in the snapshot.
- Kinetic interface-potential map versus finite-image `infinity_barrier`.
- Profile return versus finite-image open `potential_barrier` at z-high.
- `photoelectron_density_model="kinetic_mean"` outer density versus a fictitious statistical return deposit on the surface.

Validation fails closed on major unsupported combinations, but numerical convergence and physical applicability still require user
verification.

## Hold the selected field fixed during outer flight

The 1-D instant return calculates outer flight time but does not add it to global simulation time. The ordinary explicit path
holds the batch-start field fixed and stops fail-closed when flight time relative to `field_evolution_timescale` exceeds
`max_frozen_field_ratio`. An `implicit_mean` photoelectron instead uses quasistationary shadow sampling with mean-solved
$k=0$ and batch-start $k\ne0$ fixed. Its flight time and ratio remain in the same diagnostics, but an over-limit shadow does
not stop the run. Limits for ordinary `same_batch` particles and ambient species are unchanged. The matching `kinetic_1d` +
`zhao_charge_driven` configuration can instead retain 1-D events across batches and count return/escape when each event becomes
due. Its outer domain ends at $L=10\lambda_{D,pe}$: a particle reaching $L$ escapes into the reservoir, and no Robin tail beyond
$L$ is used to classify return. The terminal state is resolved with the field at enqueue and is not reintegrated after that
field changes. Queue mode applies the frozen-field bound to each event's `tau_outer`, its delay to the next batch-start poll,
and the half-batch bound on midpoint crossing-time uncertainty. Configuration validation applies the same bound to
`batch_duration`.
See
[Particle escape and return](ParticleEscapeReturn.en.html#queue-outer-flight-for-the-transient-zhao-closure).

`implicit_mean` does not solve delayed return current during UV turn-on or retain an outer inventory between batches. Those
effects require a separately designed delayed inventory or queue compatible with the time-dependent closure.

## Check convergence and balance component by component

| Target | Parameters to vary | Quantities to compare |
| --- | --- | --- |
| Near/far periodic field | Image layer, Ewald layer, cache cold/warm | $\phi,\mathbf E$, force, operator residual |
| Physical `k=0` | Lower closure and height/grid refinement | Gauss residual and interface field |
| Kinetic profile | Debye length, source sampling, outer stride | $\phi_I$, currents, nonlinear residual |
| Reservoir | Macro target and batch duration | Inflow current and macro residual |
| Photoelectron | Ray count, `dt`, outer return | Emission, reabsorption, and escape/return charge |
| Implicit mean $k=0$ | `batch_duration`, ray count, macro-particle count, lower boundary | $\phi_I$, species currents, `mean_solver_iterations`, `transaction_residual_C`, `sample_escape_fraction`, `return_weight_scale` |
| Transient Zhao queue | `batch_duration`, time-scale guard, ray count, area, interface location | $\eta$, column residual, return/escape current, force |
| Batch coupling | `batch_duration` | Steady surface charge and current balance |

Inspect `summary.txt`, `outer_plasma_profile.csv`, `charge_ledger.csv`, and `charges.csv` together. Checking only field, particle
count, or surface charge can miss component duplication or charge loss.

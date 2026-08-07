title: `batch_duration` Stability and Steady Value

Lang: [English](BatchDurationStability.en.md) | [日本語](BatchDurationStability.md)

# `batch_duration` Stability and Steady Value

This task guide explains how to select `sim.batch_duration` and check the time-width sensitivity of charging results.
In a fixed-width run, it is the physical time and surface-charge update width of one batch.
When `sim.batch_duration_step` is set, `sim.batch_duration = sim.dt * sim.batch_duration_step`.

> `sim.tol_rel` is a monitoring and output value. The current implementation does not use it for early stopping;
> the run continues to the accepted batch count set by `sim.batch_count`.

## Select a fixed width

### Prerequisites

- `boundary_inflow`, `plane_source`, `reservoir_face`, and `photo_raycast` require a positive `sim.batch_duration`.
- Set `history_stride > 0` to save `charge_history.csv` for comparison.
- Use the same mesh, particle distributions, RNG seed, and OpenMP/MPI layout.

### Action

1. Start with a conservative, small `batch_duration` or `batch_duration_step`.
2. Run three cases at 1/2, 1, and 2 times the reference value.
3. Compare the runs near the same physical time.
4. Check `last_rel_change`, total charge, and absorbed/escaped counts in `summary.txt`,
   together with element charge in `charge_history.csv`.

### Expected output

Each run creates `summary.txt` and `charge_history.csv`, allowing comparison of:

- final surface-charge distribution
- total charge and local-potential range
- oscillation, divergence, and Monte Carlo jitter in charge history
- `simulated_time_s` and accepted batch count

### Interpretation

| Observation | Decision |
| --- | --- |
| Final charge and history nearly match at 1/2 and 1 times | The reference width is a practical candidate |
| Increasing the width causes oscillation or divergence | Lower `batch_duration` |
| Final charge changes strongly with width | Recompute with a smaller width |
| Noise obscures the history | Adjust `w_particle` or `target_macro_particles_per_batch` |
| Change continues at the end of the run | Increase `batch_count` |

Successful completion alone does not establish stability or a steady state. This comparison is a step-size sensitivity check,
not Richardson extrapolation, because it does not assume a power law for the error.

### Next choices

- The reference and half-width runs agree: retain the smaller width as the validation baseline, or use the reference width when cost matters.
- Deterministic oscillation remains: lower `batch_duration`.
- Noise dominates: adjust macro-particle count or weight before changing the time width.
- A `cached_kneq0` run must bound the local-potential change of one batch: use adaptive progression below.

## Use adaptive $k\ne0$ progression

### Prerequisites

This path requires:

- `[periodic2].nonzero_mode_backend = "cached_kneq0"`
- time-scaled `boundary_inflow`, `plane_source`, `reservoir_face`, or `photo_raycast`
- `target_macro_particles_per_batch`, rather than fixed `w_particle`, for reservoir inflow and surface sources
- a positive `sim.batch_duration`

This path does not support `volume_seed`.

### Action

```toml
[periodic2]
nonzero_mode_backend = "cached_kneq0"
max_nonzero_mode_potential_step = 1.0e-2 # V
```

Let $h_0$ be the resolved `sim.batch_duration`. Each accepted batch tests
$h_0,h_0/2,h_0/4,\ldots$ in order. BEACH evaluates the $k\ne0$ potential produced by the difference
between candidate and batch-start charge at every panel centroid, and accepts the first trial whose maximum absolute value
does not exceed the limit.

### Expected output

- A rejected trial rolls back the RNG and macro-particle residuals and does not appear in statistics, history, or the charge ledger.
- `simulated_time_s` in `summary.txt` is the actual physical end time.
- `batch_count` counts accepted batches, so physical time is not generally
  `batch_count * batch_duration`.

Rejected trials within one run use a fixed OpenMP team size. A restart may use a different team size, so compare charge
distributions and accepted widths numerically instead of requiring bitwise identity across the restart.

### Interpretation

`max_nonzero_mode_potential_step` is a local-potential trust bound for freezing the $k\ne0$ field within a batch.
It does not guarantee local truncation error or an order of global accuracy, and it does not control stability of the $k=0$ update.

For target-count reservoirs and fixed-`rays_per_batch` photo sources, halving the trial width also halves macro-particle charge.
A limit-halving comparison therefore mixes time-discretization changes with Monte Carlo variance changes. Use the same RNG seed
and report both charge-distribution norms and particle statistics.

### Next choices

1. Halve `max_nonzero_mode_potential_step`.
2. Compare surface charge, local-potential range, total charge, and particle statistics near the same `simulated_time_s`.
3. For a fixed-width control, omit the key or set it to `0`.

## Theory needed for interpretation

Let $q_j$ be the accumulated charge of insulator element $j$, $A_j$ its area, and
$J_j(\mathbf q)$ its incident charge flux. The mean model is

$$
\frac{dq_j}{dt}=J_j(\mathbf q)A_j.
$$

One batch can be viewed as an explicit update that freezes the field at batch start:

$$
\mathbf q^{n+1}
=\mathbf q^n+\Delta t_b\,\mathbf J(\mathbf q^n)\mathbf A+\boldsymbol\eta^n,
$$

where $\boldsymbol\eta^n$ is Monte Carlo error.

The fixed point of the mean update satisfies $\mathbf J(\mathbf q^\ast)=0$, so the fixed point itself does not depend on
$\Delta t_b$ when the mean model converges stably. Actual runs still contain finite-sample and finite-time stopping errors,
and their observed results can retain time-width dependence.

The general linear stability condition near a fixed point is

$$
\rho(I+\Delta t_b M)<1.
$$

Only when the dominant eigenvalues are real negative and the fastest response can be represented by $\tau_{\min}$,
$\Delta t_b<2\tau_{\min}$ is a non-divergence guide and $\Delta t_b<\tau_{\min}$ is a monotone-convergence guide.
These are not general BEACH CFL conditions.

The inverse plasma frequency $\omega_{pe}^{-1}$ and charging time
$\tau_\text{charge}\sim C_\text{eff}/G_\text{eff}$ provide separate physical scales, but the actual limit also depends on
geometry, potential, and the inflow distribution. Select the final value with the step-size sensitivity check.

## Related documents

- [Input parameter reference](Parameters.en.html) — `sim.batch_duration` and `sim.batch_duration_step`
- [boundary-reservoir inflow and velocity sampling](ReservoirInjection.en.html) — particle count and weight
- [Validate Results](ValidationGuide.en.html) — numerical convergence and physical validity
- [Computational model overview](Algorithms.en.html) — batch loop

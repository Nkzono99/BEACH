title: How to choose `batch_duration`

Lang: [English](BatchDurationStability.en.md) | [日本語](BatchDurationStability.md)

# How to choose `batch_duration`

This page answers how large the physical time represented by one batch may be.
A single run cannot establish a safe value. The basic method is to compare fixed widths at 1/2, 1, and 2 times a
reference value and find a range that does not change the result at the same physical time.

After reading this page, you can run a fixed-width comparison and decide whether to retain that width or use adaptive
`cached_kneq0` progression.

> `sim.tol_rel` is a monitoring and output value. The current implementation does not use it for early stopping;
> the run continues to the accepted batch count set by `sim.batch_count`.

## Choose the path first

| Goal | Path |
| --- | --- |
| Check time-width sensitivity for a general case | Compare fixed widths at 1/2, 1, and 2 times |
| Bound the local-potential change of one `cached_kneq0` batch | Use `max_nonzero_mode_potential_step` |
| Check overall stability or global accuracy including $k=0$ | Run a fixed-width comparison; adaptive progression alone is insufficient |

`sim.batch_duration` is the physical time of one batch and the time width over which surface charge is updated.
When `sim.batch_duration_step` is used, the resolved width is

$$
\texttt{batch\_duration}=\texttt{dt}\times\texttt{batch\_duration\_step}.
$$

The two keys are mutually exclusive. `sim.dt` is the particle-push time step, while `batch_duration` controls how long
surface-charge changes are accumulated before they are applied.

`boundary_inflow`, `plane_source`, `reservoir_face`, `photo_raycast`, and
`surface_charge_closure="fixed_current"` require a positive resolved `batch_duration`.
See the [Input parameter reference](Parameters.en.html) for the complete input contract.

## Compare fixed widths

### Hold the comparison conditions fixed

- Set `write_files=true` and `history_stride > 0` to save `summary.txt` and `charge_history.csv`.
- Keep the mesh, particle distributions, RNG seed, OpenMP thread count, and MPI rank count fixed.
- Give each run a separate `output.dir`.
- Change `batch_count` together with the width so that the runs are compared near the same `simulated_time_s`.

For a reference width $h$ and reference batch count $N$, use:

| Run | `batch_duration` | `batch_count` | `output.dir` |
| --- | ---: | ---: | --- |
| half | $h/2$ | $2N$ | `outputs/batch-half` |
| reference | $h$ | $N$ | `outputs/batch-reference` |
| double | $2h$ | $N/2$ | `outputs/batch-double` |

Choose an even $N$. For example, with $h=1.0\times10^{-7}$ s and $N=100$, the half, reference, and double settings are
`(5.0e-8, 200)`, `(1.0e-7, 100)`, and `(2.0e-7, 50)`, respectively.
These numbers demonstrate the edit and are not recommended physical defaults.
When using `batch_duration_step`, keep `dt` fixed and scale the step value by the same factors of 1/2, 1, and 2.

### Create and run the configurations

Copy the original configuration three times, then change only the corresponding `[sim]` and `[output]` values from the
table above.

```bash
cp beach.toml batch-half.toml
cp beach.toml batch-reference.toml
cp beach.toml batch-double.toml
```

Validate the inputs, then run them locally or on a compute node.

```bash
beachx lint batch-half.toml
beachx lint batch-reference.toml
beachx lint batch-double.toml

beach batch-half.toml
beach batch-reference.toml
beach batch-double.toml
```

See [Run a simulation](Execution.en.html) for choosing an execution environment.

### Check the expected outputs

Each run needs `summary.txt`, `charges.csv`, and `charge_history.csv`.
Check the completed batch count and physical end time with:

```bash
for output_dir in outputs/batch-half outputs/batch-reference outputs/batch-double; do
  beachx inspect "$output_dir"
  grep -E '^(batches|last_rel_change|simulated_time_s)=' "$output_dir/summary.txt"
done
```

The expected state is:

- All three runs complete with exit code `0`.
- `batches` equals `sim.batch_count` in each configuration.
- `simulated_time_s` reaches the same physical end time in all three runs.
- The `batch` column in `charge_history.csv` can be converted with each run's width to compare element charge at
  corresponding physical times.

Successful completion does not establish numerical stability or a steady state. `last_rel_change` is also a diagnostic
for comparing history and final state, not a stopping condition.

### Choose the width

Compare the final surface-charge distribution, total charge, local-potential range, absorbed/escaped counts, and history
oscillation. Define “agreement” from the accuracy required for the study's quantities of interest before comparing runs.

| Observation | Decision |
| --- | --- |
| Half and reference agree within the chosen tolerance | The reference is a practical candidate; retain half as the validation baseline |
| The reference or double run oscillates or diverges | Lower `batch_duration` |
| Final charge changes systematically with width | Repeat the comparison at smaller widths |
| Monte Carlo noise obscures the history | Adjust `w_particle` or `target_macro_particles_per_batch` first |
| Change continues at the end | Increase `batch_count` while keeping the physical end time aligned |

This is a step-size sensitivity check, not Richardson extrapolation. It assumes neither a power law for the error nor a
particular convergence order. [`batch_duration` theory](BatchDurationTheory.en.html) explains why.

## Use adaptive $k\ne0$ progression

### Cases where it applies

This path is an advanced option for an existing periodic2 case and requires all of the following:

- `[periodic2].nonzero_mode_backend = "cached_kneq0"`
- a positive `sim.batch_duration`, or a `sim.batch_duration_step` that resolves to a positive value
- time-scaled `boundary_inflow`, `plane_source`, `reservoir_face`, or `photo_raycast`
- `target_macro_particles_per_batch`, rather than fixed `w_particle`, for reservoir inflow and `plane_source`

A `volume_seed` with positive `npcls_per_step` cannot be used on this path.

### Set the limit and run

Copy the original periodic2 configuration, then change the corresponding existing values under `[periodic2]` and
`[output]` in `adaptive.toml`.

```bash
cp beach.toml adaptive.toml
```

```toml
[periodic2]
nonzero_mode_backend = "cached_kneq0"
max_nonzero_mode_potential_step = 1.0e-2 # V

[output]
dir = "outputs/adaptive"
```

`1.0e-2` V is an input example. Set it from the acceptable local-potential change, then compare it with a run at half the
limit.

```bash
beachx lint adaptive.toml
beach adaptive.toml
beachx inspect outputs/adaptive
```

Let $h_0$ be the resolved `sim.batch_duration`. For every accepted batch, BEACH tests
$h_0,h_0/2,h_0/4,\ldots$. It evaluates the $k\ne0$ potential change produced by the candidate charge at every panel
centroid and accepts the first width whose maximum absolute value does not exceed the limit.

### Check the expected outputs

```bash
grep -E '^(batches|simulated_time_s|periodic2_max_nonzero_mode_potential_step_V|adaptive_nonzero_mode_rejected_trials|adaptive_nonzero_mode_last_batch_duration_s|adaptive_nonzero_mode_last_potential_step_V|adaptive_nonzero_mode_omp_threads)=' \
  outputs/adaptive/summary.txt
```

- `adaptive_nonzero_mode_last_batch_duration_s` is the last accepted width.
- `adaptive_nonzero_mode_last_potential_step_V` is the last accepted $k\ne0$ potential change and does not exceed the
  configured limit.
- `adaptive_nonzero_mode_rejected_trials=0` is valid. The field is the cumulative number of halvings needed to satisfy
  the limit.
- `simulated_time_s` accumulates accepted widths and is not generally equal to
  `batch_count * batch_duration`.

A rejected retry restores the RNG and macro-particle residuals and does not update statistics, history, or the charge
ledger. The run stops if one batch still exceeds the limit after 24 halvings. In that case, examine the field model,
particle statistics, and charge change as well as reducing the maximum width.

### Decide whether to retain the limit

1. Create a run with half the `max_nonzero_mode_potential_step`.
2. Compare surface charge, local-potential range, total charge, and particle statistics near the same `simulated_time_s`.
3. If they agree within the chosen tolerance, retain the larger limit as a practical candidate.
4. For a fixed-width control, omit the key or set it to `0`.

`max_nonzero_mode_potential_step` is a local-potential trust bound for freezing the $k\ne0$ field within a batch.
It does not guarantee local truncation error, an order of global accuracy, or stability of the $k=0$ update. Adaptive
progression therefore still requires a fixed-width or limit-halving sensitivity check.

For target-count reservoirs and fixed-`rays_per_batch` photo sources, halving the trial width also halves macro-particle
charge. A limit-halving comparison therefore includes both time-discretization and Monte Carlo variance changes. Use the
same RNG seed and report charge-distribution norms together with particle statistics.

Rejected retries within one run use a fixed OpenMP team size. A restart may use a different team size, so check numerical
agreement of charge distributions and accepted widths rather than bitwise identity across the restart.

## Related documents

- [`batch_duration` theory](BatchDurationTheory.en.html) — fixed points, linear stability, and applicability limits
- [Input parameter reference](Parameters.en.html) — key types and complete constraints
- [Boundary-reservoir inflow and velocity sampling](ReservoirInjection.en.html) — particle count and weight
- [Output Format Reference](OutputReference.en.html#adaptive-batch-diagnostics) — adaptive receipts
- [Validate Results](ValidationGuide.en.html) — numerical convergence and physical validity
- [BEACH computational cycle](Algorithms.en.html) — update order with a frozen field within each batch

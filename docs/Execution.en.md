title: Run a simulation

Lang: [English](Execution.en.md) | [日本語](Execution.md)

# Run a simulation

Handle an existing `beach.toml` in the order of input validation, execution, and output inspection. Estimate workload before
large runs, and resume from an existing checkpoint with the same configuration.
Start with [Installation](Installation.en.html) if BEACH is not installed, or the
[10-minute tutorial](Tutorial.en.html) if you do not have a case yet.

## Basic flow

```bash
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

1. `beachx lint` checks TOML, JSON Schema, coordinate and placement combinations, and known constraints.
2. `beach` runs the simulation from the configuration file.
3. `beachx inspect` reads the results under `output.dir`.

With no argument, `beach` reads `beach.toml` from the current directory.

## Parallel execution

Set the OpenMP thread count through the environment.

```bash
OMP_NUM_THREADS=8 beach beach.toml
```

Launch an MPI build through an MPI runner.

```bash
mpirun -n 4 beach beach.toml
```

See [Development and testing](Workflow.en.html) for combined MPI/OpenMP execution, compiler settings, and
KUDPC-specific operation.

Internally, OpenMP distributes particle indices with `dynamic, 1` to reduce imbalance from different particle
lifetimes. Collision charge is accumulated thread-locally in `dq_thread(nelem, nth)` before being combined.
MPI splits particle generation and tracking across ranks while each rank retains the same mesh and `q_elem`.
The batch-end allreduces of `dq(nelem)` and outcome counts keep committed surface charge and statistics identical
on every rank. The root writes normal results, history, and the global macro residual, while RNG state is saved per rank.

## Estimate the workload

For `reservoir_face` and `photo_raycast`, particle counts per batch are dynamic, so estimate the workload first.

```bash
beachx workload beach.toml --threads 8
```

Include MPI rank distribution when needed.

```bash
beachx workload beach.toml \
  --threads 8 \
  --mpi-ranks 4 \
  --mpi-rank 0
```

## Profiling

Set `BEACH_PROFILE=1` to record coarse phase timings.

```bash
BEACH_PROFILE=1 OMP_NUM_THREADS=8 beach beach.toml
beachx profile outputs/latest/performance_profile.csv \
  --save outputs/latest/performance_profile.png
```

For scaling comparisons, use `rank_max_s` on the `simulation_total` row of `performance_profile.csv`.
The measured regions cover initialization, batch preparation, field refresh, particle tracking, charge commit,
MPI reduction, statistics and history updates, and result and checkpoint output. They appear as `load_or_init`,
`field_solver_init`, `prepare_batch`, `field_refresh`, `particle_batch`, `commit_charge`, `mpi_reduce`,
`stats_update`, `history_write`, `write_results`, and `write_checkpoint`.

## Resume a run

Resume from the same output directory with:

```toml
[output]
dir = "outputs/latest"
resume = true
```

Use `restart_from` to separate the checkpoint from the new output directory.

```toml
[output]
dir = "outputs/continuation"
resume = true
restart_from = "../parent_run/outputs/latest"
```

`sim.batch_count` is the cumulative target, not an additional batch count. If a checkpoint has `batches=100`
and the new target is `batch_count=150`, BEACH runs 50 more batches. An MPI resume requires the current rank
count to match the saved run.
An adaptive $k\ne0$ restart also requires the actual OpenMP team size saved by
the checkpoint. Unequal team sizes across MPI ranks or a checkpoint mismatch
fail fast.

## Checks after execution

- The process exit code is zero.
- `batches` in `summary.txt` equals `sim.batch_count`.
- Review `absorbed`, `escaped_boundary`, and `survived_max_step`.
- Do not treat `tol_rel` as an automatic stopping condition.

Continue with [Inspect Output Files](OutputGuide.en.html), the
[Post-processing tutorial](PostprocessTutorial.en.html), and
[Validate simulation results](ValidationGuide.en.html).

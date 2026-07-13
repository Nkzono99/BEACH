title: Run a simulation

Lang: [English](Execution.en.md) | [日本語](Execution.md)

# Run a simulation

This page covers checking and running an existing `beach.toml`, estimating its workload, and resuming a run.
Start with [Installation](Installation.en.html) if BEACH is not installed, or the
[10-minute tutorial](Tutorial.en.html) if you do not have a case yet.

## Basic flow

```bash
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

1. `beachx lint` checks TOML, JSON Schema, high-level notation, and known constraints.
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

## Checks after execution

- The process exit code is zero.
- `batches` in `summary.txt` equals `sim.batch_count`.
- Review `absorbed`, `escaped_boundary`, and `survived_max_step`.
- Do not treat `tol_rel` as an automatic stopping condition.

Continue with [Output files](OutputGuide.en.html), the
[Post-processing tutorial](PostprocessTutorial.en.html), and
[Validate simulation results](ValidationGuide.en.html).

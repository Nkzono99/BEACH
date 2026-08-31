title: Run and resume a simulation

Lang: [English](Execution.en.md) | [日本語](Execution.md)

# Run and resume a simulation

This page explains how to run an existing `beach.toml` safely and continue from a checkpoint when needed. Proceed through
input validation, workload estimation, execution, and output inspection. For a restart, preserve the source output and
increase the cumulative `batch_count` in a separate configuration.

Start with [Installation](Installation.en.html) if BEACH is not installed, or the
[10-minute tutorial](Tutorial.en.html) if you do not have a case yet.

Do not run `beach` or `mpirun` directly on a KUDPC login node. Run the commands below either locally or inside a
compute-node allocation. On KUDPC, use `tssrun` for a short check and `srun` inside `sbatch` for a long run; see the
[`examples/job_scripts/`](https://github.com/Nkzono99/BEACH/tree/main/examples/job_scripts) examples.
On other HPC systems, follow the site's login-node execution policy.

## Basic flow

This example uses the official tutorial's `output.dir="outputs/tutorial"`.
For another case, replace the final argument with that case's `output.dir`.

```bash
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/tutorial
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

Match the MPI build, launcher, and compiler modules to the execution environment. Developers changing MPI / OpenMP state
ownership or reduction should see [runtime architecture](Architecture.en.html).

## Estimate the workload

For `boundary_inflow`, `plane_source`, deprecated `reservoir_face`, and `photo_raycast`, particle counts per batch are
dynamic, so estimate the workload first.

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
beachx profile outputs/tutorial/performance_profile.csv \
  --save outputs/tutorial/performance_profile.png
```

For scaling comparisons, use `rank_max_s` on the `simulation_total` row of `performance_profile.csv`. The profile separates
initialization, field updates, particle tracking, charge updates, MPI, and output phases.

## Resume once from a checkpoint

A normally completed final output remains restartable even with `checkpoint_stride=0`.
For the first resume check, preserve the original `outputs/tutorial` directory and write the continuation to a separate
directory.

### Extend the official tutorial by one batch

**Prerequisite:** Complete the [10-minute tutorial](Tutorial.en.html) and confirm that the working directory contains
`beach.toml` and `outputs/tutorial`. Read the completed batch count:

```bash
grep '^batches=' outputs/tutorial/summary.txt
```

The official tutorial reports `batches=20`. For a general case, call this completed count `B`.

Keep the original configuration and create a resume configuration:

```bash
cp beach.toml resume.toml
```

Change the corresponding existing values under `[sim]` and `[output]` in `resume.toml`.
This complete example extends the official tutorial from 20 batches to 21:

```toml
[sim]
batch_count = 21

[output]
write_files = true
dir = "outputs/resumed"
resume = true
restart_from = "outputs/tutorial"
```

Leave the other settings in `resume.toml` unchanged. `restart_from` is the checkpoint input, while `dir` is the new
output directory. From the same tutorial working directory, validate the input, run it, and inspect the new output:

```bash
beachx lint resume.toml
beach resume.toml
beachx inspect outputs/resumed
```

On success, the `beach` and `beachx inspect` output includes:

```text
resuming_from_batches=20
...
batches=21 ...
```

For another case with completed count `B`, set `batch_count` to the integer `B+1` and check for
`resuming_from_batches=B` and `batches=B+1`.
`sim.batch_count` is the cumulative target, not the number of additional batches.

### Periodic checkpoints for a long run

Enable periodic checkpoints by accepted-batch count:

```toml
[output]
dir = "outputs/tutorial"
checkpoint_stride = 1000
```

This saves restart state to two alternating slots every 1000 accepted batches. `0` disables periodic saves while
preserving the final checkpoint written after normal completion.

### Continue in the same directory and resume MPI runs

To continue writing into the checkpoint directory, omit `restart_from` and set `dir` to that directory.
Use the separate-directory procedure above when you want to preserve the original final output for comparison.

```toml
[output]
write_files = true
dir = "outputs/tutorial"
resume = true
```

In general, a checkpoint at `batches=100` with a new `batch_count=150` runs 50 additional batches.
For an MPI resume, match the checkpoint `mpi_world_size` to the current rank count.
Adaptive $k\ne0$ retries require equal actual OpenMP team sizes across MPI ranks within the current run.
A restart may use a different team size; BEACH records the resumed run's actual team size as a new diagnostic value.

[Files used for resume](OutputReference.en.html#files-used-for-resume) is the source of truth for required files and checkpoint
selection. If a resume is rejected, use
[Troubleshooting](Troubleshooting.en.html#cannot-resume-from-a-checkpoint) to check required files, fingerprints,
MPI world size, and the cumulative `batch_count`.

## Checks after execution

- The process exit code is zero.
- `batches` in `summary.txt` equals `sim.batch_count`.
- Review `absorbed`, `escaped_boundary`, and `survived_max_step`.
- Do not treat `tol_rel` as an automatic stopping condition.

Continue with [Inspect Output Files](OutputGuide.en.html), the
[Post-processing tutorial](PostprocessTutorial.en.html), and
[Validate simulation results](ValidationGuide.en.html).

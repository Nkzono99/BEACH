title: Development and Operations Workflow

Lang: [English](Workflow.en.md) | [日本語](Workflow.md)

# Development and Operations Workflow

This task guide covers source development from local execution through HPC verification.
Fortran runs the simulation; Python provides configuration utilities, post-processing, and visualization.
For normal use without source changes, see [Installation](Installation.en.html) and [Run a Simulation](Execution.en.html).

## Run a case in an installed environment

**Prerequisite:** Install `beach-bem` and make `beach` and `beachx` available on `PATH`.

```bash
python -m pip install -U pip setuptools wheel
python -m pip install beach-bem
```

To try the development version directly:

```bash
python -m pip install "git+https://github.com/Nkzono99/BEACH.git"
```

**Action:**

```bash
mkdir beach-tutorial
cd beach-tutorial
beachx config init beach.toml
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

**Expected output:** The run creates `outputs/latest/summary.txt`, `charges.csv`, and mesh information.

**Interpretation:** Successful completion proves that the configuration, build, and execution path work.
It does not establish steady-state convergence or physical validity.

**Next choices:** Use [Design a Simulation Case](ConfigurationRecipes.en.html) to change the case,
[Inspect Output Files](OutputGuide.en.html) to interpret files, and [Validate Results](ValidationGuide.en.html) to assess validity.

## Create a development environment

**Prerequisites:** Install Python, `make`, a Fortran compiler, and `fpm`.

```bash
make --version
gfortran --version
fpm --version
python --version
```

**Action:**

```bash
python -m pip install -U pip setuptools wheel
python -m pip install -e . --no-build-isolation
make check
```

**Expected output:** The Python package is installed in editable mode and the Fortran source compiles.

**Interpretation:** `make check` is the lightweight development build check. It uses `BEACH_VERSION_MODE=dev`,
so fpm can reuse incremental compilation when only the Git hash changes.

**Next choices:**

```bash
make run CONFIG=examples/beach.toml
make build VERSION_MODE=plain
make build VERSION_MODE=git
make install-generic
```

Normally use `make run` and `make check` through `build.sh`. Reserve direct fpm execution for low-level checks:

```bash
fpm run --profile release --flag "-fopenmp" -- examples/beach.toml
```

## Test a change

**Prerequisite:** Identify the changed area and the required test tier. Do not run multiple `fpm test` commands concurrently.

**Action:**

```bash
make test-l0      # static/schema/build
make test         # L1: normal development loop
make test-l2      # L2: contract/integration
make test-l3      # L3: cumulative verification including heavy FMM
```

To run one Fortran target:

```bash
FPM_ACTION=test ./build.sh --target test_version
```

**Expected output:** Each command completes without a nonzero status and its targets pass.

**Interpretation:**

- L0: `git diff --check`, JSON schema, and `make check`
- L1: L0 + Python tests + lightweight Fortran targets
- L2: L1 + integration targets such as the C/kernel contract
- L3: L2 + heavy FMM targets

The heavy `test_dynamics_fmm` and `test_coulomb_fmm_core_basic` targets are not part of L1.

**Next choices:**

```bash
make test-heavy
make test-fortran-far-correction
make test-fortran-benchmark
make test-field-kernel-cache
make test-full
```

`test-field-kernel-cache` is an opt-in native cache/plane-oracle receipt gate and is not part of L1/L2/L3.
Use [Physics Release Verification](PhysicsReleaseVerification.en.html) for release decisions.

## Run with OpenMP or MPI

**Prerequisite:** Prepare an OpenMP build, or an MPI build made with an MPI compiler.

**Action:**

```bash
OMP_NUM_THREADS=8 beach beach.toml
mpirun -n 4 beach examples/beach.toml
```

To add coarse phase profiling:

```bash
BEACH_PROFILE=1 OMP_NUM_THREADS=8 beach examples/beach.toml
beachx profile outputs/latest/performance_profile.csv \
  --save outputs/latest/performance_profile.png
```

**Expected output:** The run creates the normal simulation outputs and, when profiling is enabled,
`performance_profile.csv`.

**Interpretation:** Use `rank_max_s` on the `simulation_total` row for scaling comparisons.
For an MPI restart, the checkpoint `mpi_world_size` must match the current rank count.

**Next choice:** To test only the hybrid path:

```bash
FPM_FC=mpiifx \
fpm test --target test_mpi_hybrid \
  --flag "-fpp -DUSE_MPI -qopenmp" \
  --runner "mpirun -n 2"
```

For an Intel `ifx` OpenMP build, use `OPENMP_FLAG=-qopenmp`.

## Run and test on KUDPC

**Prerequisites:** Inspect `hostname`, `module list`, and, when available, `spartition` and `qgroup`
to identify the host and allocation.

**Action:** On a login node, limit work to editing, short log inspection, `make check`, job submission, and monitoring.
Run `make test*`, `fpm test`, long simulations, and benchmarks on a compute node.

- Short interactive check: `tssrun`
- Batch execution: `srun` inside an `sbatch` job

Example for 112 OpenMP threads per rank on SysA:

```bash
export OMP_NUM_THREADS=112
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
srun beach beach.toml
```

**Expected output:** The simulation or test completes inside a Slurm allocation and leaves a job log and normal outputs.

**Interpretation:** Fix thread placement for performance comparison and reproducibility. Login-node behavior is not a measure
of compute-node performance.

**Next choice:** Use the repository KUDPC plugin and
`examples/job_scripts/camphor_mpi_hybrid_job.sh` for environment and job configuration.

## Estimate workload before a run

**Prerequisite:** In a case with `boundary_inflow`, `plane_source`, deprecated `reservoir_face`, or `photo_raycast`, the
configuration determines particles per batch.

**Action:**

```bash
beachx workload examples/beach.toml --threads 8
```

To include one MPI rank and checkpoint residuals:

```bash
beachx workload examples/beach.toml \
  --threads 8 \
  --mpi-ranks 4 \
  --mpi-rank 0 \
  --macro-residuals outputs/latest/macro_residuals.csv
```

**Expected output:** The command reports rank-local `total_particles` and all-rank `global_total_particles`.

**Interpretation:** This is a workload estimate, not a wall-time forecast or a physical-validity check.

**Next choice:** If the estimate is large, begin with a smoke case using a coarser mesh, fewer macro-particles, and a smaller `batch_count`.

## Resume from a checkpoint

**Prerequisite:** The source directory contains a checkpoint set including `summary.txt`, `charges.csv`, and RNG state.

**Action:**

```toml
[output]
dir = "outputs/continuation"
resume = true
restart_from = "../parent_run/outputs/latest"
```

**Expected output:** BEACH loads the checkpoint and writes new outputs to `outputs/continuation`.

**Interpretation:** `sim.batch_count` is the cumulative target. If the checkpoint has `batches=100` and the new
`batch_count=150`, the resumed run adds 50 batches. A target below the completed batch count is rejected.

**Next choice:** To continue in the same directory, omit `restart_from` and set `dir` to the checkpoint directory.

## Contracts to check after a change

- A normal run proceeds to the batch count set by `sim.batch_count`.
- `sim.tol_rel` is a monitoring and output value, not an early-stop condition.
- Configure box geometry and periodic axes in `[domain]`, and the field closure in `[field_boundary]`.
- Configure global particle faces in `[particle_boundary]` and species overrides in `[particles.species.boundary]`.
- Put the external-reservoir inflow model and `phi_infty` in `[reservoir]`, and select inflow faces in
  `[particles.species.boundary_inflow]`.
- The standard v1.0 surface model is insulator accumulation.
- Boundary reservoir + closed PE is a reduced closure inside the finite box; it does not solve an external region self-consistently.
- Execution success, numerical convergence, and physical validity require separate checks.

When changing a public API, configuration, or output, update the Japanese and English documentation, examples, schema,
and corresponding tests in the same change.

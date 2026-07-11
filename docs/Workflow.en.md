title: Execution and Development Workflow

Lang: [English](Workflow.en.md) | [日本語](Workflow.md)

# Execution and Development Workflow

The project is centered on the **Fortran runtime**. Python handles post-processing, visualization, and utility
workflows. For normal users, the recommended path is to install `beach-bem` and run the `beach` command.

## 1. User Setup (Recommended)

### 1.1 Check Tools

```bash
make --version
gfortran --version
fpm --version
python --version
```

### 1.2 Install from PyPI

```bash
python -m pip install -U pip setuptools wheel
python -m pip install beach-bem
```

During `pip install`, `make install` installs both the Python CLI and the Fortran runtime binary.
Pip builds use `INSTALL_PROFILE=auto` by default and fall back to `generic` when needed.
Set `BEACH_PIP_FALLBACK_GENERIC=0` to disable that fallback.

```bash
export PATH="$HOME/.local/bin:$PATH"
```

To try the development version directly from Git:

```bash
python -m pip install "git+https://github.com/Nkzono99/BEACH.git"
```

## 2. Development Workflow

### 2.1 Python Editable Install

```bash
python -m pip install -U pip setuptools wheel
python -m pip install -e . --no-build-isolation
```

### 2.2 Fortran Runtime (`make`)

```bash
make check
make run CONFIG=examples/beach.toml
```

`make check` is the standard development build check. It uses `BEACH_VERSION_MODE=dev` to pass a stable
version string such as `1.5.0-dev` to the Fortran side, so changes to the git hash do not invalidate fpm's
compile-flag hash and incremental builds stay reusable.

`make build` and `make install` embed a git-hash version by default. Override the version mode when needed.

```bash
make build VERSION_MODE=dev
make build VERSION_MODE=plain
make build VERSION_MODE=git
```

Install profiles can also be selected explicitly.

```bash
make install-generic
make install-camphor
```

### 2.3 Direct `fpm` Execution

```bash
fpm run --profile release --flag "-fopenmp" -- examples/beach.toml
```

For normal development, prefer `make run` / `make check` through `build.sh`. The wrapper passes
`__BEACH_VERSION__` and `__BEACH_VERSION_MODE__` consistently.

### 2.4 Tests

```bash
make test-l0      # L0: static/schema/build check
make test         # L1: normal development loop
make test-l2      # L2: contract/integration
make test-l3      # L3: cumulative L0-L3 verification
make test-physics-release  # HPC: minimal release correctness + MPI manifest
make test-heavy   # heavy Fortran targets only
make test-fortran-far-correction  # oracle far-correction correctness
make test-fortran-far-correction-diagnostics  # assertion-free diagnostics
make test-fortran-benchmark  # release-profile runtime benchmark
make test-full    # unfiltered fpm test
```

The test suite is tiered for the development loop.

- L0: `git diff --check`, JSON schema parse check, `make check`
- L1: L0 + Python tests + lightweight Fortran test targets (`make test` / `make test-l1`)
- L2: L1 + contract/integration targets such as the C field-kernel contract
- L3: L2 + heavy FMM targets for release gates, nightly checks, and integration before `main`

`make test-fortran` aliases the lightweight Fortran targets. Heavy FMM targets such as
`test_dynamics_fmm` and `test_coulomb_fmm_core_basic` are excluded from normal `make test` and must be run with
`make test-l3`, `make test-heavy`, `make test-fortran-heavy`, or `make test-full`.
The `m2l_root_oracle` correctness tests are opt-in through `make test-fortran-far-correction`.
The assertion-free `test_periodic2_flat_oracle_diag` is separated under
`make test-fortran-far-correction-diagnostics`. Runtime comparison uses the release profile through
`make test-fortran-benchmark` instead of running inside a debug correctness test.
Tiered tests under Intel `ifx` / `mpiifx` suppress only the
`arg_temp_created` check by default because each expected array temporary can
otherwise emit a full stack trace. Other debug checks, including bounds
checks, remain enabled. To inspect array temporaries themselves, override the
flags explicitly, for example with
`FORTRAN_TEST_FLAGS="-qopenmp" make test-fortran-heavy FPM_FC=mpiifx`.

Run a single target with:

```bash
FPM_ACTION=test ./build.sh --target test_version
```

On KUDPC login nodes, do not run `make test*`, `fpm test`, or equivalent build/test payloads directly.
Submit them to compute nodes with `tssrun` or `sbatch`.

## 3. Run Flow

Usually, edit `beach.toml` and pass that file directly to `beach`. See
[`beachx config` / High-Level Notation Guide](Configuration.en.html) for the high-level notation layer.

1. Prepare `beach.toml`.
2. Check it with `beachx lint beach.toml`.
3. Run the simulation with `beach beach.toml`.
4. Inspect files under `output.dir`.
5. Visualize with Python CLI commands or the `Beach` API.

See [Input Parameters Reference](Parameters.en.html) for the `beach.toml` specification.

### 3.1 Shortest Example

```bash
mkdir run_periodic2
cd run_periodic2
beachx config init
beachx lint beach.toml
beach beach.toml
```

### 3.2 Direct `beach.toml` Use

1. Prepare `beach.toml` (see [Input Parameters Reference](Parameters.en.html)).
2. Use high-level notation if useful; the Fortran parser resolves it while loading.
3. Run the simulation with `beach beach.toml`.
4. Inspect `output.dir`.
5. Visualize with a Python CLI command or the `Beach` API.

## 4. Run Commands

### 4.1 Recommended: `beach`

```bash
beach beach.toml
```

Without arguments, `beach` reads `beach.toml` from the current directory.

### 4.2 Thread Count

```bash
OMP_NUM_THREADS=8 beach beach.toml
```

### 4.3 MPI + OpenMP

After installing an MPI build, launch `beach` through the MPI runner.

```bash
mpirun -n 4 beach examples/beach.toml
```

### 4.4 Profiling

Enable coarse phase profiling with:

```bash
BEACH_PROFILE=1 OMP_NUM_THREADS=8 beach examples/beach.toml
```

For scaling comparisons, use the `rank_max_s` value on the `simulation_total` row in `performance_profile.csv`.

Visualization example:

```bash
beachx profile outputs/latest/performance_profile.csv \
  --save outputs/latest/performance_profile.png
```

## 5. Output Files

Primary outputs:

- `summary.txt`
- `charges.csv`
- `mesh_potential.csv` when `write_mesh_potential = true`
- `mesh_triangles.csv`
- `mesh_sources.csv`
- `charge_history.csv` when `history_stride > 0`
- `potential_history.csv` when `write_potential_history = true` and `history_stride > 0`
- `performance_profile.csv` when `BEACH_PROFILE=1`
- `rng_state.txt`
- `macro_residuals.csv`

`mesh_triangles.csv` includes `mesh_id` for each element. `mesh_sources.csv` maps each `mesh_id` to the source
template kind, surface model, `epsilon_r`, and element count. `conductor` is relaxed as a floating conductor
only with `field_bc_mode = "free"`. `dielectric` is metadata-only in the current implementation; `summary.txt`
also prints a note when it appears. Enabling `mesh_potential.csv` stores centroid potential [V] with the same
element ordering.

MPI runs (`world_size > 1`) write only the RNG state per rank. The reservoir residual is shared by all ranks,
so the root rank writes one `macro_residuals.csv`.

- `rng_state_rank00000.txt`, `rng_state_rank00001.txt`, ...
- `macro_residuals.csv`

## 6. Workload Estimation

For `reservoir_face` / `photo_raycast`, particle counts per batch are dynamic, so estimating the workload first
is recommended.

```bash
beachx workload examples/beach.toml --threads 8
```

Rank-local estimate with residuals:

```bash
beachx workload examples/beach.toml \
  --threads 8 \
  --mpi-ranks 4 \
  --mpi-rank 0 \
  --macro-residuals outputs/latest/macro_residuals.csv
```

`total_particles` is for the selected rank and `global_total_particles` is the all-rank total.
For reservoir injection, compare `local_reservoir_particles` with `global_reservoir_particles`.

## 7. Resume Runs

```toml
[output]
dir = "outputs/latest"
resume = true
```

Rerunning `beach` with the same `output.dir` reads `summary.txt`, `charges.csv`, and RNG state and continues
from the checkpoint. Use `restart_from` when the checkpoint directory and the new output directory differ.

```toml
[output]
dir = "outputs/continuation"
resume = true
restart_from = "../parent_run/outputs/latest"
```

In that case, checkpoint files are read from `restart_from`, while the new `summary.txt`, `charges.csv`,
history files, and RNG state are written under `dir`.

`sim.batch_count` is the cumulative target batch count. If an existing checkpoint has `batches=100` and the new
`batch_count=150`, only 50 additional batches run. If `batch_count` is smaller than the completed checkpoint
batch count, execution stops.

For MPI resume, the `mpi_world_size` in `summary.txt` must match the current number of ranks.

## 8. Python Post-processing

### 8.1 CLI

```bash
beachx inspect outputs/latest \
  --save-bar outputs/latest/charges_bar.png \
  --save-mesh outputs/latest/charges_mesh.png \
  --save-potential-mesh outputs/latest/potential_mesh.png \
  --potential-self-term area-equivalent

# Draw a sim.field_bc_mode = "periodic2" mesh wrapped into the periodic cell.
beachx inspect outputs/latest \
  --save-mesh outputs/latest/charges_mesh_periodic.png \
  --save-potential-mesh outputs/latest/potential_mesh_periodic.png \
  --apply-periodic2-mesh

# Tile the periodic mesh by n layers. For 1, this draws 3x3 = 9 copies.
beachx inspect outputs/latest \
  --save-mesh outputs/latest/charges_mesh_tiled.png \
  --periodic2-repeat 1 \
  --apply-periodic2-mesh

beachx animate outputs/latest \
  --quantity charge \
  --save-gif outputs/latest/charge_history.gif \
  --total-frames 200

beachx slices outputs/latest \
  --grid-n 200 \
  --vmin -20 --vmax 20 \
  --save outputs/latest/potential_slices.png

beachx coulomb outputs/latest \
  --component z \
  --save outputs/latest/coulomb_force_z.png

beachx mobility outputs/latest \
  --density-kg-m3 2500 \
  --mu-static 0.4 \
  --save-csv outputs/latest/mobility_summary.csv

# Use the same field kernel as the Fortran FMM core to output net charge, force, and torque per object.
make build-kernel
beachx kernel-forces outputs/latest \
  --save-csv outputs/latest/object_forces_kernel.csv

# Retain periodic images while excluding only the central-cell primary self field.
beachx object-detachment outputs/latest \
  --config beach.toml \
  --target-mesh-id 6 \
  --periodic-model infinite-physical \
  --z-max-m 2.0e-4 \
  --z-points 65 \
  --mass-kg 2.0e-12 \
  --gravity-m-s2 9.80665 \
  --output-dir outputs/latest/object_detachment
```

`beachx coulomb` reads object kind and order from nearby `beach.toml` `mesh.templates` when available, and by
default places all objects along the target axis for visualization. Use `--target-kinds sphere` to restrict the
target set. `beachx mobility` treats `plane` as the support by default and writes object force/torque plus
`lift_ratio`, `slide_ratio`, and `roll_ratio` to CSV. Mass-derived indicators need `--density-kg-m3` and geometry
from `beach.toml`. `beachx kernel-forces` calls the Fortran FMM core through `libbeach_field_kernel`, using
`sim.softening`, `sim.field_bc_mode`, periodic2, and tree settings from `beach.toml`. Build the library with
`make build-kernel`; if it is elsewhere, pass `--library` or set `BEACH_FIELD_KERNEL_LIB`. Use
`--config path/to/beach.toml` when no config exists near the output directory.

`kernel-forces` is the legacy `exclude_target_lattice` diagnostic, which also
removes the target object's periodic images. `object-detachment` uses
`exclude_primary_keep_images`. It writes the instantaneous wrench, frozen-source
vertical path and work, and a from-rest barrier including gravity and optional
finite-range adhesion to four artifacts. `configured` preserves the run's
policy, while `infinite-physical` uses cached `k != 0` plus the x/y-periodic
`E_bottom=0` zero mode. The CLI defaults to lunar gravity, `1.62 m/s^2`; the
example above explicitly uses Earth gravity, `9.80665 m/s^2`. A cold cached
operator or a long path can be expensive;
on KUDPC submit this analysis to a compute node instead of running it on a login
node.

A successful CLI invocation establishes artifact generation only. It is not a
physical qualification until path status, work/potential agreement,
quadrature, shell/cache, and endpoint sensitivity have been checked. A
finite-height speed in a non-neutral periodic cell is not escape speed at
infinity.

Legacy aliases such as `beach-inspect`, `beach-animate-history`, `beach-plot-coulomb-force-matrix`,
`beach-plot-potential-slices`, `beach-estimate-workload`, and `beach-plot-performance-profile` remain available
for now, but are deprecated.

### 8.2 Python API

```python
from beach import Beach

beach = Beach("outputs/latest")
print(beach.result.absorbed, beach.result.escaped)
# History is always lazy-loaded.
history_step10 = beach.result.history_at(10)
if beach.result.history is not None:
    print(beach.result.history.batch_indices)

beach.plot_bar()
beach.plot_mesh()
beach.plot_potential()
beach.plot_mesh(apply_periodic2_mesh=True)
beach.plot_potential(apply_periodic2_mesh=True)
beach.plot_potential_slices(
    box_min=[0.0, 0.0, 0.0],
    box_max=[1.0, 1.0, 10.0],
    grid_n=200,
    vmin=-20.0,
    vmax=20.0,
)
beach.animate_mesh("outputs/latest/charge_history.gif", quantity="charge", total_frames=200)

mesh1 = beach.get_mesh(1)
mesh2, mesh3 = beach.get_mesh(2, 3)
mesh1_step10 = beach.get_mesh(1, step=10)
charge_step10 = beach.get_mesh_charge(1, step=10)

interaction = beach.calc_coulomb(target=[mesh1, mesh2], source=[mesh3], step=10)
print(interaction.force_on_a_N, interaction.torque_on_a_Nm)

fig_force, ax_force = beach.plot_coulomb_force_matrix(
    component="z",
)
fig_force.savefig("outputs/latest/coulomb_force_z.png", dpi=150)

mobility = beach.analyze_coulomb_mobility(
    density_kg_m3=2500.0,
    mu_static=0.4,
)
for record in mobility.records:
    print(record.label, record.lift_ratio, record.slide_ratio)
```

## 9. MPI Path Test (Developers)

```bash
FPM_FC=mpiifort \
fpm test --target test_mpi_hybrid \
  --flag "-fpp -DUSE_MPI -qopenmp" \
  --runner "mpirun -n 2"
```

## 10. Easy-to-Misread Implementation Behaviors

- A normal run advances exactly `sim.batch_count` batches. A resume run advances from the checkpoint batch count
  until `sim.batch_count` is reached.
- `sim.tol_rel` is a monitoring value. It is not used as an early-stop condition in the current implementation.
- The Fortran electric field uses element-centroid point-charge approximation plus `sim.softening`.

For a camphor MPI job example, see `examples/job_scripts/camphor_mpi_hybrid_job.sh`.
`test-physics-release` sequentially runs the L1 convergence subset, L3-heavy, far-correction correctness,
MPI ledger, and MPI periodic-cache gates without repeating the full portable L2 suite.
It records the commit, dirty state, host, compilers, status, elapsed time, and peak RSS for
each gate in `build/physics-release/manifest.txt` by default. It refuses KUDPC
login nodes and selects `srun` for the MPI payload inside a Slurm allocation.
Override the output with `PHYSICS_RELEASE_MANIFEST=/path/to/manifest.txt`.
The same directory receives `convergence.csv`; see
[Physics release verification](PhysicsReleaseVerification.en.md).
It also receives `test_l3-target-timings.csv` and `far_correction-target-timings.csv`, containing profile,
status, and elapsed seconds for each selected Fortran target.

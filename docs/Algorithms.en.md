title: BEACH Algorithm Overview

Lang: [English](Algorithms.en.md) | [日本語](Algorithms.md)

# BEACH Algorithm Overview

This document describes the numerical algorithms and execution order used by the current Fortran implementation of BEACH.
BEACH is not a grid PIC code. It is a surface charging simulator that couples a BEM-like Coulomb field evaluation, with charges on triangular boundary elements as sources, to test-particle tracking in batches.

Implementation links point to the current file and the main symbol. If line numbers drift in later edits, prefer the linked file and symbol name.

---

## Algorithm Document Structure

| Document | Contents |
|---|---|
| [Algorithm Overview](Algorithms.en.html) | BEACH computation model, initialization, and batch loop |
| [Field Solvers and Boundary Conditions](FieldSolvers.en.html) | Coulomb field, direct/treecode/FMM switching, and periodic2 field boundary |
| [Particle Tracking and Charge Accumulation](ParticleChargeLoop.en.html) | Particle generation, Boris pusher, collision, charge deposition, statistics, and restart |
| [Coulomb FMM Core Details](FMMCore.en.html) | FMM core API, tree construction, M2L, and periodic2 Ewald/oracle |
| [`batch_duration` Stability](BatchDurationStability.en.html) | Batch time step, steady value, linear stability, and Monte Carlo noise |


## 1. BEACH computation model

**Source**:
[`bem_simulator`](../src/runtime/simulator/bem_simulator.f90),
[`bem_simulator_loop`](../src/runtime/simulator/bem_simulator_loop.f90),
[`bem_particle_stepper`](../src/runtime/simulator/bem_particle_stepper.f90),
[`bem_field_solver`](../src/physics/field_solver/bem_field_solver.f90),
[`bem_injection`](../src/particles/bem_injection.f90)

The main state in BEACH is the charge `q_elem(i)` on triangular mesh elements plus the particles generated for each batch.
The mesh geometry is fixed. Particles are absorbed when they hit a surface. The absorbed particle charge is accumulated on the hit element and affects the field calculation in the next batch.

### 1.1 State variables

The main data types are defined in [`bem_types`](../src/core/bem_types.f90).

| Type | Main contents | Role |
| --- | --- | --- |
| `sim_config` | `dt`, `batch_count`, `max_step`, `field_solver`, `field_bc_mode`, `box_min/max`, `bc_low/high` | Time advance, field calculation, boundary conditions |
| `mesh_type` | `v0/v1/v2`, `centers`, `normals`, `bb_min/max`, `q_elem`, `elem_surface_model` | Triangular boundary elements and charge |
| `particles_soa` | `x`, `v`, `q`, `m`, `w`, `species_id`, `alive` | SoA representation of particles in a batch |
| `injection_state` | `macro_residual(:)` | Carries fractional `reservoir_face` macro-particle counts across batches and restarts |
| `sim_stats` | `processed_particles`, `absorbed`, `escaped`, `survived_max_step`, `batches`, `last_rel_change` | Run statistics |

### 1.2 Physical approximation

Each triangular element contributes to the field as a point charge `q_i` at its centroid `c_i`.
The current implementation uses a centroid point-charge approximation, not an exact BEM integration over a continuous surface charge distribution.

Softened Coulomb kernel:

$$
G_\epsilon(\mathbf{r}) = \frac{1}{\sqrt{\lVert\mathbf{r}\rVert^2 + \epsilon^2}}
$$

Potential:

$$
\phi(\mathbf{x}) =
k_\mathrm{c}
\sum_j q_j G_\epsilon(\mathbf{x} - \mathbf{c}_j)
$$

Electric field:

$$
\mathbf{E}(\mathbf{x}) =
k_\mathrm{c}
\sum_j q_j
\frac{\mathbf{x} - \mathbf{c}_j}
{\left(\lVert\mathbf{x} - \mathbf{c}_j\rVert^2 + \epsilon^2\right)^{3/2}}
$$

Here `k_coulomb` is the Coulomb constant from [`bem_constants`](../src/core/bem_constants.f90).
The field applied to particles is the boundary-element field plus the uniform external field `sim.e0`.

### 1.3 Execution unit

BEACH advances batches up to `sim.batch_count`. Each batch means:

1. Generate particles for each species from the current surface charge state.
2. Refresh the field solver with the current `q_elem`.
3. Track each particle for at most `sim.max_step` steps.
4. Accumulate charge from collided particles into element-wise deltas `dq`.
5. Add `dq` to `q_elem` and apply surface-model relaxation if needed.
6. Update statistics and history.

`sim.tol_rel` is an output and monitoring value. In the current Fortran implementation it is not an early-stop condition.

---

## 2. Execution entry point and initialization

**Source**:
[`app/main.f90`](../app/main.f90),
[`bem_app_config_runtime`](../src/config/bem_app_config_runtime.f90),
[`bem_mesh`](../src/mesh/bem_mesh.f90),
[`bem_restart`](../src/runtime/bem_restart.f90)

### 2.1 Main program order

`app/main.f90` is the CLI entry point. Its high-level order is:

| Order | Processing | Main implementation |
| --- | --- | --- |
| 1 | Handle CLI options that finish before config loading, such as `--version` | `handle_early_cli` |
| 2 | Initialize MPI and the performance profile | `mpi_initialize`, `perf_configure_from_env` |
| 3 | Load config, build mesh, read restart or initialize RNG seed | `load_or_init_run_state` |
| 4 | Open history CSV writers | `open_history_writer`, `open_potential_history_writer` |
| 5 | Run the batch simulation | `run_absorption_insulator` |
| 6 | Write the summary and final CSV files | `print_run_summary`, `write_result_files` |
| 7 | Write RNG state and macro residuals as checkpoints | `write_rng_state_file`, `write_macro_residuals_file` |
| 8 | Write performance profile and shut down MPI | `perf_write_outputs`, `mpi_shutdown` |

If a config path is passed explicitly, BEACH uses it. Otherwise it uses `beach.toml` in the current directory.
If no config file exists, BEACH runs with `default_app_config` defaults.

### 2.2 Mesh construction

The mesh is built from templates or an OBJ file according to `mesh.mode`. For template meshes,
[`build_template_mesh`](../src/config/bem_app_config_runtime.f90) does the following:

1. Iterate over `mesh.templates`.
2. Skip templates with `enabled=false`.
3. Dispatch by `kind` to `make_plane`, `make_box`, `make_cylinder`, `make_sphere`, and related builders.
4. Assign a `mesh_id` per template.
5. Expand `surface_model` and `epsilon_r` to element arrays.
6. Concatenate all template triangle arrays and pass them to `init_mesh`.

`init_mesh` precomputes:

- element centroids `centers(:, i)`
- element normals `normals(:, i)`
- element AABBs `bb_min/max(:, i)`
- representative length `h_elem(i) = sqrt(area_i)`
- initial charge `q_elem(i)`
- collision grid

### 2.3 periodic2 collision mesh

When `sim.field_bc_mode="periodic2"`, `prepare_periodic2_collision_mesh` translates each triangle to a canonical position along the two periodic axes.
This is separate from the periodic image sum used in the field calculation. It stabilizes the primitive-cell mesh for collision tests.
Element indices remain those of the base element, so a hit on a periodic image is still deposited on the base element.

### 2.4 Restart

When `output.resume=true`, `load_restart_checkpoint` reads:

- `summary.txt`: completed batch count and statistics
- `charges.csv`: element charge
- `rng_state.txt` or `rng_state_rankNNNNN.txt`: RNG state
- `macro_residuals.csv` or rank-local residual files: reservoir fractional particle counts
- `charge_ledger.csv`: cumulative signed-charge ledger for schema v2

`sim.batch_count` is the cumulative target batch count. If the checkpoint has `batches=100` and `batch_count=150`, the resumed run executes only 50 batches.
Schema v2 stores model, ordered-mesh, and ordered-species fingerprints and checks them before mutating runtime state. Legacy schemas are accepted only for implemented legacy point-source models.

The particle-integration contract is `particle_time_centering="same_time_midpoint_boris"`. Pure-E, pure-B, time-reversal, smooth-field second-order convergence, and batch-restart-continuity regressions detect changes to this contract.

---

## 3. Batch loop

**Source**:
[`run_absorption_insulator`](../src/runtime/simulator/bem_simulator_loop.f90),
[`prepare_batch_state`](../src/runtime/simulator/bem_simulator_loop.f90),
[`process_particle_batch`](../src/runtime/simulator/bem_simulator_loop.f90),
[`commit_batch_charge`](../src/runtime/simulator/bem_simulator_loop.f90)

### 3.1 Loop skeleton

When `run_absorption_insulator` receives `initial_stats`, it resumes from `initial_stats%batches`.

```text
final_batch_idx = sim.batch_count
batch_count_this_run = final_batch_idx - stats.batches

field_solver.init(mesh, sim)

for local_batch_idx = 1..batch_count_this_run:
    prepare_batch_state(...)
    field_solver.refresh(mesh)
    process_particle_batch(...)
    commit_batch_charge(...)
    count_batch_outcomes(...)
    MPI allreduce statistics
    accumulate_batch_stats(...)
    write charge/potential history when stride matches
```

### 3.2 Batch state

`prepare_batch_state` prepares:

- current batch number `batch_idx = stats%batches + 1`
- particles from `init_particle_batch_from_config`
- `dq_thread(nelem, nth)`: charge delta per OpenMP thread
- `photo_emission_dq(nelem)`: source-side opposite charge from `photo_raycast`
- `escaped_boundary_flag(:)` and `absorbed_flag(:)`

Splitting `dq_thread` by thread lets collision deposition be accumulated without atomics.
When a `photo_raycast` collision query is incomplete, the sampler does not stop inside OpenMP. It returns the lowest
ray / bounce and status through `init_particle_batch_from_config` to `prepare_batch_state`. Before field refresh or charge
processing, the main loop selects the lowest-rank species / ray / bounce / status and every rank stops with the same diagnostic.
Incomplete particle arrays and `photo_emission_dq` from the failing rank are not used by later processing.

### 3.3 Particle processing

`process_particle_batch` advances each particle for at most `sim.max_step` time steps.

| Order | Processing |
| --- | --- |
| 1 | Read the same-time current position `x0` and velocity `v0` |
| 2 | Form the predicted midpoint `x_mid = x0 + 0.5*v0*dt` in `build_particle_step_candidate` |
| 3 | Evaluate the boundary-element field with `field_solver%eval_e(mesh, x_mid, e_mid)` and add the uniform external field `sim.e0` exactly once |
| 4 | Use `boris_push` with the midpoint field to compute candidate velocity `v1` and trapezoidal candidate position `x1` |
| 5 | Issue one collision query on `x0 -> x1`; if `x1` is inside the box, commit that result |
| 6 | For a box crossing, compare the mesh-hit and first-face parameters, giving mesh hits priority on ties |
| 7 | Open ends at the event; reflect/periodic rebuilds the remainder once and checks that chord for mesh hits |
| 8 | Deposit `q * w` on a hit; fail closed without committing state at a second box event |
| 9 | If the particle survives, update same-time `x` and `v` and continue to the next step |

`build_particle_step_candidate` evaluates the spatial field exactly once at the predicted midpoint without modifying
the field solver. `boris_update_velocity(v, q, m, dt, e, b, v_new)` is a public pure procedure that performs the
electric half kick, magnetic rotation, and electric half kick for the velocity update. The existing public call
`boris_push(x, v, q, m, dt, e, b, x_new, v_new)` keeps its signature, delegates its velocity calculation to this
procedure, and updates position with `x_new = x + 0.5*(v + v_new)*dt`. Input and output positions and velocities are
same-time states. Predicted-midpoint spatial-field sampling and the trapezoidal position update make the candidate
kinematics second-order accurate.

When the candidate endpoint is strictly inside the box, `advance_particle_step` completes with one field evaluation and
one collision query and does no additional event geometry. Crossing steps use `find_first_boundary_event` and
`apply_escape_reflect_periodic_event` to apply simultaneous corner/edge faces together. If a reflected or periodic
remainder reaches another face without an earlier mesh hit, it returns `particle_step_multiple_box_events` instead of
entering an unbounded loop. The existing `apply_box_boundary` remains for photo rays.
If a periodic2 full-chord query reaches a range limit beyond the box, the production loop falls back to a query truncated
at the first box event.

For a single open face, `potential_barrier` retains the legacy scalar-energy formula evaluated at the event position and
interpolated velocity. Multiple open faces fail closed; a shared-potential/gauge physical model is deferred.

If `BEACH_WARN_LONG_PARTICLE_STEPS` is set to a positive integer, BEACH prints diagnostics for long-lived particles at that step interval.

The collision statuses are `collision_query_ok=0`, `collision_query_image_limit=1`, and
`collision_query_index_range=2`. A public call that omits `status` does not treat an incomplete query as a miss; it stops
fail closed. The main loop aggregates collision and boundary-event failures in one envelope, leaves the OpenMP region,
then shares the selected rank's particle/step and `dt/x/v` so every rank stops with the same diagnostic.

### 3.4 Charge commit

`commit_batch_charge` merges thread-local deltas and photo-emission deltas.

$$
\Delta q_i =
\sum_{\mathrm{thread}} \Delta q_{i,\mathrm{thread}}
{}+ \Delta q_{i,\mathrm{photo}}
$$

In MPI runs, `mpi_allreduce_sum_real_dp_array` sums `dq` across ranks.
Then BEACH applies:

$$
q_i^{new} = q_i^{old} + \Delta q_i
$$

and runs surface-model relaxation. The monitoring value is computed from the actually changed charge:

$$
\operatorname{rel\_change} =
\frac{\lVert\mathbf{q}^{\mathrm{new}} - \mathbf{q}^{\mathrm{old}}\rVert_2}
{\max\left(\lVert\mathbf{q}^{\mathrm{new}}\rVert_2, q_\mathrm{floor}\right)}
$$

This value is stored in `stats%last_rel_change` and written to history output.

---

# Stage 2 Particle Events And MPI Injection Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Preserve the synchronized second-order Boris step, order mesh and box events correctly with negligible no-crossing overhead, and make reservoir count/residual behavior MPI-global across runtime, restart, and Python estimates.

**Architecture:** Normal particle steps retain one midpoint field evaluation and one collision query. Only a candidate endpoint outside the box enters boundary geometry; open crossings terminate, while reflect/periodic crossings rebuild the remaining interval once. A second box crossing in the same outer `dt` fails closed with context instead of adding adaptive or unbounded event loops. Reservoir expected counts are computed once globally, split across ranks, and checkpointed as one root-owned residual vector.

**Tech Stack:** Fortran 2008, fpm, OpenMP, optional MPI, Python 3.10+, pytest, KUDPC Slurm.

## Global Constraints

- Work on branch `codex/review-remediation`; never rewrite the user's pre-existing commit history.
- Run Fortran/Python payload tests on KUDPC compute nodes through `tssrun -p eb`; do not run them directly on `laurel31`.
- Do not run multiple fpm controllers concurrently because they share `build/`.
- Use test-first RED/GREEN cycles and commit each independently reviewable task.
- Preserve the public `boris_push`, `apply_box_boundary`, `run_absorption_insulator`, and `load_restart_checkpoint` call surfaces.
- Keep `open_boundary_model="potential_barrier"` as a legacy experimental closure. Do not add a new barrier model or public barrier API in Stage 2.
- Do not add adaptive-substep configuration or a general event loop. A second box event in one outer step is an explicit numerical failure advising a smaller `sim.dt`.
- Default no-crossing work remains one field evaluation and one collision query; representative runtime slowdown must be no more than 3% outside measurement noise.
- Global reservoir count/residual sequences must not depend on MPI world size. Rank ownership may differ.
- MPI-size-changing restart remains rejected because rank-owned RNG streams are not portable.

---

## Current Baseline

- [x] **P1:** `e1470d3 refactor: separate the Boris velocity update`
- [x] **P2:** `f8b36fd fix: use same-time second-order Boris kinematics`

The production candidate is the symmetric drift/kick/rotate/kick/drift form:

```fortran
x_mid = x0 + 0.5_dp*v0*dt
call field_solver%eval_e(mesh, x_mid, e_mid)
e_mid = e_mid + sim%e0
call boris_push(x0, v0, q, m, dt, e_mid, bfield, x1, v1)
```

No later Stage 2 task may restagger stored velocity or add another field evaluation to the no-crossing path.

## Dependency Order

```text
P1 -> P2 -> B1 -> B2 -> PERF
                  |
M1 -> M2 -> M3 ---+-> final L1/MPI verification
```

## Task B1: Lightweight First Box Event Geometry

**Files:**
- Modify: `src/physics/bem_boundary.f90`
- Modify: `tests/fortran/test_boundary.f90`
- Modify: `SPEC.md`
- Modify: `docs/Algorithms.md`
- Modify: `docs/Algorithms.en.md`

**Interfaces:**

```fortran
integer(i32), parameter, public :: boundary_event_ok = 0_i32
integer(i32), parameter, public :: boundary_event_invalid_geometry = 1_i32

type, public :: boundary_event_type
  logical :: has_event = .false.
  real(dp) :: fraction = 0.0_dp
  integer(i32) :: face_mask = 0_i32
  integer(i32) :: face_bc(6) = -1_i32
end type boundary_event_type

subroutine find_first_boundary_event(cfg, x0, x1, event, status)
  type(sim_config), intent(in) :: cfg
  real(dp), intent(in) :: x0(3), x1(3)
  type(boundary_event_type), intent(out) :: event
  integer(i32), intent(out) :: status
end subroutine find_first_boundary_event

subroutine apply_escape_reflect_periodic_event(cfg, event, x, v, alive, escaped, status)
  type(sim_config), intent(in) :: cfg
  type(boundary_event_type), intent(in) :: event
  real(dp), intent(inout) :: x(3), v(3)
  logical, intent(inout) :: alive
  logical, intent(out) :: escaped
  integer(i32), intent(out) :: status
end subroutine apply_escape_reflect_periodic_event
```

Face bits are `x-low=1`, `x-high=2`, `y-low=4`, `y-high=8`, `z-low=16`, `z-high=32`. All faces within `64*epsilon(1.0_dp)` of the minimum fraction form one simultaneous event. Detection requires a finite start inside the closed box and finite endpoints/bounds with positive spans. Event application is transactional; invalid input does not mutate state. Reflect and periodic destinations use `nearest(face, inward_sign)`.

`apply_escape_reflect_periodic_event` treats any open face as terminal escape. It rejects `open_boundary_model="potential_barrier"`; B2 owns the existing single-face legacy energy decision privately. `apply_box_boundary` remains unchanged for photo-ray and compatibility callers.

- [x] Add compile-time and behavioral tests for single-face, corner, mixed reflect/periodic, endpoint-on-face, inward/outward boundary start, disabled box, invalid geometry, and state non-mutation.
- [x] Run RED (job `22256702`, missing `boundary_event_type`):

```bash
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 bash -lc \
  'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_boundary'
```

- [x] Implement the additive event geometry/action API without changing legacy callers.
- [x] Run GREEN with the same command (job `22256723`, 20 tests / 74 assertions), format the modified Fortran file, and run `git diff --check`.
- [x] Commit as `feat: add lightweight box boundary events`.

## Task B2: Mesh/Box Earliest Ordering With One Remainder

**Files:**
- Modify: `src/runtime/simulator/bem_particle_stepper.f90`
- Modify: `src/runtime/simulator/bem_simulator.f90`
- Modify: `src/runtime/simulator/bem_simulator_loop.f90`
- Modify: `tests/fortran/test_particle_stepper.f90`
- Modify: `tests/fortran/test_simulator.f90`
- Modify: `SPEC.md`
- Modify: `docs/ParticleChargeLoop.md`
- Modify: `docs/ParticleChargeLoop.en.md`
- Modify: `CHANGELOG.md`

**Interfaces:**

```fortran
integer(i32), parameter, public :: particle_step_ok = 0_i32
integer(i32), parameter, public :: particle_step_invalid_boundary = 1001_i32
integer(i32), parameter, public :: particle_step_multiple_box_events = 1002_i32
integer(i32), parameter, public :: particle_step_unsupported_barrier_corner = 1003_i32

type, public :: particle_step_result
  real(dp) :: x(3) = 0.0_dp
  real(dp) :: v(3) = 0.0_dp
  logical :: absorbed = .false.
  logical :: escaped_boundary = .false.
  integer(i32) :: elem_idx = -1_i32
  integer(i32) :: status = particle_step_ok
  integer(i32) :: field_eval_count = 0_i32
  integer(i32) :: collision_query_count = 0_i32
end type particle_step_result

subroutine advance_particle_step(mesh, sim, field_solver, bfield, x0, v0, q, m, dt, result)
  type(mesh_type), intent(in) :: mesh
  type(sim_config), intent(in) :: sim
  type(field_solver_type), intent(inout) :: field_solver
  real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt
  type(particle_step_result), intent(out) :: result
end subroutine advance_particle_step
```

**Algorithm:**

1. Build the full candidate once.
2. Fast path: if `use_box=false` or all candidate coordinates are inside the closed box, query `x0 -> x1` once and accept/absorb.
3. Crossing path: find the first box event, query only `x0 -> x_event`, and give mesh hits at the event-time tie priority over the box.
4. Default open terminates at the event. A single-face legacy potential-barrier event evaluates the existing scalar closure at the event point/velocity; simultaneous multi-open barrier events return `particle_step_unsupported_barrier_corner` without inventing new physics.
5. Reflect/periodic transforms the event state and rebuilds the remaining interval once using the same synchronized Boris candidate.
6. If the remainder endpoint is inside, query its chord and accept/absorb.
7. If it is outside, find the second box event and query only up to it. Absorb an earlier/tied mesh hit; otherwise return `particle_step_multiple_box_events` with state uncommitted and advise reducing `sim.dt`.
8. Propagate collision statuses unchanged and aggregate all particle-step failures outside the OpenMP region with batch/rank/particle/step context.

- [x] Add RED fixtures for mesh-before-box absorption, box-before-mesh open escape, reflected-remainder absorption, periodic remainder, simultaneous corner action, second-crossing failure, and the no-crossing `field_eval_count=1`/`collision_query_count=1` contract.
- [x] Run RED sequentially for `test_particle_stepper` (job `22256730`) and `test_simulator` (job `22256739`) through `tssrun -p eb`.
- [x] Implement `advance_particle_step` and migrate only the simulator particle loop; leave photo-ray boundary handling on the legacy API.
- [x] Run GREEN for `test_boundary`, `test_particle_stepper`, `test_simulator`, and `test_dynamics_basic` sequentially (job `22256746`).
- [x] Update docs, format, run `git diff --check`, and commit as `fix: order mesh and box boundary events`.

## Task M1: One Global Reservoir Count Per Species

**Files:**
- Modify: `src/core/bem_mpi.F90`
- Modify: `src/config/bem_app_config_runtime.f90`
- Modify: `tests/fortran/test_reservoir_injection.f90`
- Modify: `tests/fortran/test_mpi_hybrid.f90`
- Modify: `SPEC.md`
- Modify: `docs/ParticleChargeLoop.md`
- Modify: `docs/ParticleChargeLoop.en.md`

**Interfaces:**

```fortran
subroutine mpi_bcast_i32_array(ctx, values, root)
subroutine mpi_bcast_real_dp_array(ctx, values, root)
```

Root computes each enabled `reservoir_face` count using the full `batch_duration` and one global residual. It broadcasts `global_counts(:)` and the updated `macro_residual(:)`, after which every rank sets `counts_max(s)=mpi_split_count(global_counts(s), rank, size)`. Volume and photo sources retain their current global-count splitting. Synthetic `mpi_rank`/`mpi_size` calls without a real `mpi_context` deterministically recompute the same global count locally before splitting.

- [x] Add a fractional sequence whose expected global count is below one per batch and assert the same global sequence for 1, 2, and 4 ranks.
- [x] Extend `test_mpi_hybrid` to sum local reservoir counts and compare the global sequence.
- [x] Run RED on `test_reservoir_injection` (job `22256752`), then implement the root/broadcast path and run GREEN (job `22256757`).
- [x] Build the MPI target (job `22256760`) and run the two-rank executable directly through `tssrun --rsc p=2` without a nested MPI launcher (job `22256761`).
- [ ] Update docs and commit as `fix: distribute global reservoir counts across ranks`.

## Task M2: One Root-Owned Global Residual Checkpoint

**Files:**
- Modify: `src/runtime/bem_restart.f90`
- Modify: `app/main.f90`
- Modify: `tests/fortran/test_restart.f90`
- Modify: `tests/fortran/test_mpi_hybrid.f90`
- Modify: `SPEC.md`
- Modify: `docs/Workflow.md`
- Modify: `docs/Workflow.en.md`

`restart_macro_residual_path(...)` always returns `<out_dir>/macro_residuals.csv`. Only root writes it; each rank continues to write its own RNG checkpoint. On restart root reads the global residual and broadcasts it after validation. Rank-suffixed residual files are rejected as legacy ambiguous state rather than merged. Existing MPI-size equality remains mandatory.

- [x] Change restart tests first to require the global path for serial and MPI arguments, root-only writes, and broadcast-equivalent restoration.
- [x] Run RED on `test_restart` (job `22256768`).
- [x] Implement root ownership and load broadcast while preserving existing public signatures.
- [x] Run GREEN on `test_restart` (job `22256772`), then the real two-rank restart fixture (build `22256773`, run `22256775`).
- [ ] Update docs and commit as `fix: checkpoint one global reservoir residual`.

## Task M3: Python Workload Estimator Parity

**Files:**
- Modify: `beach/cli/estimate_fortran_workload.py`
- Modify: `tests/python/test_estimate_fortran_workload.py`
- Modify: `docs/Workflow.md`
- Modify: `docs/Workflow.en.md`
- Modify: `docs/agent-user-guide.md`
- Modify: `docs/agent-user-guide.en.md`

For each reservoir species, compute the uninterrupted global integer/residual sequence once, then split each batch count with the same quotient/remainder rule as `mpi_split_count`. Report both global and selected-rank reservoir totals. Resume reads one `macro_residuals.csv`; rank-suffixed residuals are not accepted.

- [x] Replace the old MPI-scaled reservoir expectations with hard-coded global sequences for 1, 2, and 4 ranks, plus uninterrupted/resumed equality.
- [x] Run RED with Python 3.11 on a compute node (job `22256780`; the earlier job `22256777` exposed the unsupported system Python 3.6 environment).
- [x] Implement the global algorithm without new dependencies and run GREEN plus ruff (job `22256781`; 35 tests passed).
- [ ] Update docs and commit as `fix: estimate global reservoir counts before rank splitting`.

## Task PERF: No-Crossing Performance Gate

**Files:**
- Modify only if the gate fails: the B1/B2 files above.
- Record results: `.superpowers/sdd/stage2-performance-report.md`.

- [ ] Run a fixed direct-field, no-box-crossing scenario at `f8b36fd` and at Stage 2 HEAD with identical compiler, OpenMP threads, seed, particles, and steps.
- [ ] Use at least five measured repetitions after one warmup and compare medians.
- [ ] Require `field_eval_count=1`, `collision_query_count=1`, identical physical counters, and median slowdown `<=3%` unless the observed run-to-run noise interval overlaps the difference.
- [ ] If the gate fails, profile before changing code; retain the endpoint-inside branch as the first boundary operation.

## Final Verification

- [ ] Run `make fmt-check-fortran` and `git diff --check` locally only if they do not execute payload tests.
- [ ] Run `make test-l1` through one `tssrun -p eb` controller allocation.
- [ ] Run the focused two-rank MPI integration executable through direct `tssrun --rsc p=2`.
- [ ] Confirm no new adaptive keys, no new public barrier API, and no general event loop exist.
- [ ] Review the complete branch diff from `2e47428` and record remaining physical limitations.

## Explicitly Deferred

- Adaptive particle substeps and their configuration/validation.
- More than one surviving box event per outer `sim.dt`.
- General multi-face potential-barrier physics and a shared potential snapshot.
- Atomic generation checkpoints, mesh/config fingerprints, and MPI-size-changing restart (Stage 3).
- FMM/far-correction redesign and shared potential gauge (later stages).

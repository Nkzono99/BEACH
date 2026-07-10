# Stage 2 Particle And MPI Correctness Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Track every step with the checkboxes below and obtain an independent review at each commit boundary.

**Goal:** Correct same-time particle integration, earliest mesh/box event handling, physically scaled adaptive particle substeps, and MPI-global reservoir residual semantics without changing the public simulator entry point.

**Architecture:** Separate the Boris velocity operation from the runtime particle-step engine, then make one internal stepper own candidate construction, event splitting, adaptive substeps, and fail-closed status propagation. Separately, compute reservoir counts from one global budget, broadcast the resulting state, and persist one root-owned residual while keeping rank-owned RNG state.

**Tech Stack:** Fortran 2008, fpm, OpenMP, conditional MPI (`mpif.h` through `bem_mpi`), Python 3, pytest, Ruff, KUDPC `tssrun`.

## Baseline

This is the implementation-ready Stage 2 plan for repository HEAD `70c5073` (`fix: propagate photo ray collision failures`). Stage 1 safety and data-integrity work is complete at this baseline. Implementers must branch from this exact behavior or first re-audit every frozen interface below if HEAD has moved.

## Source Of Truth And Current State

This plan implements all Stage 2 particle and MPI work in `docs/superpowers/specs/2026-07-10-review-remediation-design.md` against the Fortran and Python layout at `70c5073`, including gyro-angle, acceleration-change, and element-length adaptive substeps. Stage 1 already supplied shell-free output creation, fail-closed periodic collision queries, mixed-sign tree descent, strict history validation, precomputed-only inspect behavior, and the documented far-correction default; none of those completed changes is reopened here.

Current behavior that the tasks below must replace:

- `bem_pusher::boris_push` documents same-time `(x,v)` but advances position with `x + v_new*dt`.
- `process_particle_batch` samples `E` at `x0`, computes a complete candidate step, checks the complete chord for a mesh hit, and only then folds or kills an overshot particle at a box boundary.
- Reflection therefore omits the crossing-time velocity, the reflected remainder of the step, and collisions on that remainder.
- `reservoir_face` divides `batch_duration` by MPI size and updates one residual per rank, so low expected counts can burst simultaneously on every rank.
- MPI restart writes and reads rank-specific macro residual files even though the corrected state must be global.
- The Python workload estimator mirrors the current per-rank residual algorithm.

The Stage 1 collision-failure path is already complete at `70c5073`. Main particle processing returns `(status, particle, step)`; photo injection returns `(status, species, ray, bounce)`. `run_absorption_insulator` performs a failure-count Allreduce, selects lowest-rank metadata with `mpi_select_lowest_rank_i32_values`, and stops before field refresh or charge commit as appropriate. Stage 2 extends these paths without replacing, reordering, or bypassing them.

## Frozen Stage 1 Interfaces

The following declarations and metadata layouts are the Stage 2 entry contract. Re-audit callers before changing any of them.

```fortran
module subroutine run_absorption_insulator( &
  mesh, app, stats, history_unit, history_stride, initial_stats, inject_state, mpi, &
  mesh_potential_v, potential_history_unit)
  type(mesh_type), intent(inout) :: mesh
  type(app_config), intent(in) :: app
  type(sim_stats), intent(out) :: stats
  integer, intent(in), optional :: history_unit
  integer(i32), intent(in), optional :: history_stride
  type(sim_stats), intent(in), optional :: initial_stats
  type(injection_state), intent(inout), optional :: inject_state
  type(mpi_context), intent(in), optional :: mpi
  real(dp), allocatable, intent(out), optional :: mesh_potential_v(:)
  integer, intent(in), optional :: potential_history_unit
end subroutine run_absorption_insulator

module subroutine process_particle_batch( &
  mesh, app, field_solver, pcls_batch, dq_thread, escaped_boundary_flag, absorbed_flag, &
  bfield, batch_idx, mpi_rank, collision_failure_status, &
  collision_failure_particle, collision_failure_step)
  type(mesh_type), intent(in) :: mesh
  type(app_config), intent(in) :: app
  type(field_solver_type), intent(inout) :: field_solver
  type(particles_soa), intent(inout) :: pcls_batch
  real(dp), intent(inout) :: dq_thread(:, :)
  logical, intent(inout) :: escaped_boundary_flag(:), absorbed_flag(:)
  real(dp), intent(in) :: bfield(3)
  integer(i32), intent(in) :: batch_idx, mpi_rank
  integer(i32), intent(out) :: collision_failure_status
  integer(i32), intent(out) :: collision_failure_particle, collision_failure_step
end subroutine process_particle_batch

module subroutine prepare_batch_state( &
  mesh, app, stats, batch_idx, dq_thread, pcls_batch, escaped_boundary_flag, absorbed_flag, &
  photo_emission_dq, mpi, inject_state, collision_failure_status, collision_failure_species, &
  collision_failure_ray, collision_failure_bounce)
  type(mesh_type), intent(in) :: mesh
  type(app_config), intent(in) :: app
  type(sim_stats), intent(in) :: stats
  integer(i32), intent(out) :: batch_idx
  real(dp), intent(inout) :: dq_thread(:, :)
  type(particles_soa), intent(out) :: pcls_batch
  logical, allocatable, intent(inout) :: escaped_boundary_flag(:), absorbed_flag(:)
  real(dp), intent(out) :: photo_emission_dq(:)
  type(mpi_context), intent(in) :: mpi
  type(injection_state), intent(inout), optional :: inject_state
  integer(i32), intent(out) :: collision_failure_status, collision_failure_species
  integer(i32), intent(out) :: collision_failure_ray, collision_failure_bounce
end subroutine prepare_batch_state

subroutine init_particle_batch_from_config( &
  cfg, batch_idx, pcls, state, mesh, photo_emission_dq, mpi_rank, mpi_size, mpi, &
  collision_failure_status, collision_failure_species, collision_failure_ray, &
  collision_failure_bounce)
  type(app_config), intent(in) :: cfg
  integer(i32), intent(in) :: batch_idx
  type(particles_soa), intent(out) :: pcls
  type(injection_state), intent(inout), optional :: state
  type(mesh_type), intent(in), optional :: mesh
  real(dp), intent(out), optional :: photo_emission_dq(:)
  integer(i32), intent(in), optional :: mpi_rank, mpi_size
  type(mpi_context), intent(in), optional :: mpi
  integer(i32), intent(out), optional :: collision_failure_status
  integer(i32), intent(out), optional :: collision_failure_species
  integer(i32), intent(out), optional :: collision_failure_ray, collision_failure_bounce
end subroutine init_particle_batch_from_config

subroutine mpi_select_lowest_rank_i32_values( &
  ctx, local_present, local_values, selected_rank, selected_values)
  type(mpi_context), intent(in) :: ctx
  logical, intent(in) :: local_present
  integer(i32), intent(in) :: local_values(:)
  integer(i32), intent(out) :: selected_rank
  integer(i32), intent(out) :: selected_values(:)
end subroutine mpi_select_lowest_rank_i32_values
```

The existing collision constants remain `collision_query_ok=0_i32`, `collision_query_image_limit=1_i32`, and `collision_query_index_range=2_i32`. Main failure selection packs `[status, particle, step]`; photo failure selection packs `[species, ray, bounce, status]`. Both use a failure-presence `mpi_allreduce_sum_i32_scalar` before lowest-rank selection. Do not change these metadata orders while implementing Stage 2.

## Global Constraints

- Review findings and implementer reports must be written in Japanese.
- Particle state is same-time `(x^n,v^n)` at every public step boundary.
- Preserve the existing positional signature and public availability of `bem_pusher::boris_push`.
- Preserve `bem_simulator::run_absorption_insulator`, `process_particle_batch`, `prepare_batch_state`, and `init_particle_batch_from_config` exactly as frozen above, including all main-particle and photo-ray failure outputs.
- Preserve `bem_boundary::apply_box_boundary` for existing callers such as photo ray casting. New event APIs are additive.
- `max_step` continues to count configured physical steps, not internal boundary events.
- A mesh hit still means absorption by default: kill the particle and deposit exactly `q_particle*w_particle` once on the base element.
- Do not turn incomplete collision work, an event-loop limit, or invalid event geometry into no-hit or escape.
- Do not execute `error stop` from the OpenMP particle region. Return status, aggregate context, and stop after the region.
- Preserve both Stage 1 collective protocols: failure-presence Allreduce, lowest-rank metadata selection, then a common stop before any invalid state is committed. Reuse `mpi_select_lowest_rank_i32_values(ctx, local_present, local_values, selected_rank, selected_values)` without changing its signature or either metadata layout.
- Global reservoir count and residual must be independent of MPI world size. Only rank ownership of those particles may change.
- Rank-specific RNG state and same-MPI-size resume remain unchanged. MPI-size-changing resume remains rejected.
- Do not introduce Stage 6 potential-kernel or Zhao changes here. The existing open-potential formula is evaluated at the corrected crossing state, but its field/gauge unification remains a later stage.
- Ignore generated `*.i90` files.
- Before any payload allocation, run `hostname -f 2>/dev/null || hostname`, `module list`, `spartition`, and `command -v qgroup >/dev/null 2>&1 && qgroup`; use the exact `eb` commands below only when that queue is visible and authorized.

## Dependency Order

```text
Stage 1 main/photo fail-closed collision status at 70c5073
  -> P1 Boris velocity refactor
  -> P2 same-time midpoint candidate
  -> B1 boundary event primitives
  -> B2 mesh/box event stepper and simulator integration
  -> B3 physically scaled adaptive substeps

M1 global reservoir count/distribution
  -> M2 global residual checkpoint semantics
  -> M3 Python workload estimator parity

P/B and M work are logically independent, but should be committed in the order above to avoid overlapping simulator and test baselines.
```

Do not run two `fpm test` commands concurrently; they share `build/`.

---

## Implementation Tasks

### Task 1: P1 - Separate Boris Velocity Update Without Changing The Public Call

#### Prerequisite

Start from `70c5073` with the frozen Stage 1 interfaces above unchanged.

#### Files

- Modify `src/physics/bem_pusher.f90`
- Modify `tests/fortran/test_dynamics_basic.f90`
- Modify `docs/Algorithms.md`, `docs/Algorithms.en.md`, and `SPEC.md`

#### Interfaces

Add a reusable procedure:

```fortran
pure subroutine boris_update_velocity(v, q, m, dt, e, b, v_new)
  real(dp), intent(in) :: v(3), q, m, dt, e(3), b(3)
  real(dp), intent(out) :: v_new(3)
end subroutine boris_update_velocity
```

Keep the current public procedure exactly callable as:

```fortran
call boris_push(x, v, q, m, dt, e, b, x_new, v_new)
```

In this first refactor commit, `boris_push` delegates its velocity calculation to `boris_update_velocity`; do not yet claim that the production simulator is second order.

#### Invariants

- The Boris rotation formula and velocity result are bitwise-equivalent to the current implementation for the same inputs.
- `dt=0` is identity for both velocity and position.
- Pure magnetic motion preserves `sum(v*v)` to roundoff.
- No particle state representation becomes half-step staggered.

#### RED Test And Command

Add direct tests for `boris_update_velocity` before the procedure exists:

- electric-only velocity update;
- pure-B speed preservation;
- forward `dt`, then backward `-dt`, recovers the original velocity.

The focused target must fail to compile or link before implementation.

```bash
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_dynamics_basic'
```

Expected RED: compile or link failure naming the missing `boris_update_velocity` symbol.

#### GREEN Test And Command

```bash
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_dynamics_basic'
```

Expected GREEN: `test_dynamics_basic` passes, including exact velocity equivalence, electric-only update, pure-B speed preservation, and forward/backward recovery.

#### Execution Steps

- [ ] Add the three RED cases and an exact `boris_push` versus `boris_update_velocity` velocity comparison.
- [ ] Run the RED command and record the missing-symbol failure.
- [ ] Extract the existing velocity statements unchanged into the pure helper; make `boris_push` call it and retain `x_new = x + v_new*dt` in this commit.
- [ ] Update the listed algorithm/SPEC wording without claiming second-order production motion, then run the GREEN command.
- [ ] Format the modified Fortran files, run `git diff --check`, and commit only the Task 1 (P1) files.

#### Commit Boundary

`refactor: separate the Boris velocity update`

This commit is behavior-preserving and green on its own.

### Task 2: P2 - Adopt Same-Time Second-Order Candidate Kinematics

#### Prerequisite

Task 1 (P1) commit `refactor: separate the Boris velocity update` is present and green.

#### Files

- Modify `src/physics/bem_pusher.f90`
- Create `src/runtime/simulator/bem_particle_stepper.f90`
- Modify `src/runtime/simulator/bem_simulator.f90`
- Modify `src/runtime/simulator/bem_simulator_loop.f90`
- Create `tests/fortran/test_particle_stepper.f90`
- Modify `tests/fortran/test_dynamics_basic.f90`
- Modify `fpm.toml` and the `FORTRAN_L1_TARGETS` list in `Makefile`
- Update `docs/Algorithms.md`, `docs/Algorithms.en.md`, `docs/ParticleChargeLoop.md`, `docs/ParticleChargeLoop.en.md`, and `SPEC.md`
- Update `CHANGELOG.md`

#### Interfaces

Change the implementation contract of the existing `boris_push` without changing its signature:

```text
v_new = Boris(v, E_mid, B_mid, dt)
x_new = x + 0.5*(v + v_new)*dt
```

Add an internal runtime procedure in `bem_particle_stepper`:

```fortran
subroutine build_particle_step_candidate( &
  mesh, sim, field_solver, bfield, x0, v0, q, m, dt, x1, v1)
  type(mesh_type), intent(in) :: mesh
  type(sim_config), intent(in) :: sim
  type(field_solver_type), intent(inout) :: field_solver
  real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt
  real(dp), intent(out) :: x1(3), v1(3)
end subroutine build_particle_step_candidate
```

It computes `x_mid = x0 + 0.5*v0*dt`, evaluates the BEM field at `x_mid`, adds `sim%e0`, and passes that midpoint field to `boris_push`.

`process_particle_batch` calls this procedure instead of evaluating `E(x0)` directly. The event behavior remains unchanged until Task 4 (B2).

#### Invariants

- Public particle inputs and outputs are same-time states.
- Constant electric acceleration gives the analytic displacement exactly up to roundoff.
- Uniform-B velocity norm is preserved and orbit position is second-order convergent.
- The complete production step evaluates the spatial electric field once, at the predicted midpoint.
- Uniform `e0` is included exactly once.
- With `q=0`, `E=0`, or `B=0`, no new special time convention is introduced.

#### RED Tests And Commands

1. Change the existing electric-only displacement expectation from `0.5` to the analytic `0.25` for `q=2`, `m=1`, `E=1`, `dt=0.5`, `x0=v0=0`.
2. Add a pure-B fixed-final-time convergence test. Compare against
   `v=(cos(omega*t),-sin(omega*t),0)` and
   `x=(sin(omega*t)/omega,(cos(omega*t)-1)/omega,0)` for the current cross-product sign. For both refinements require an error ratio in `[3.2,4.8]`.
3. In `test_particle_stepper`, use a fixed softened one-element charged mesh and compare integrations with `dt`, `dt/2`, and `dt/4` against a `dt/64` direct reference. Require both position-error ratios to be at least `3.2`. This test fails if the runtime samples only `E(x0)`.

```bash
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_dynamics_basic'
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_particle_stepper'
```

Expected RED: the corrected constant-E displacement assertion fails under the old position update; the new stepper target also fails until `bem_particle_stepper` exists.

#### GREEN Tests And Commands

```bash
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_dynamics_basic'
tssrun -p eb -t 0:15:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_particle_stepper'
tssrun -p eb -t 0:15:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_simulator'
```

Run the three commands sequentially; never overlap fpm controllers. Expected GREEN: all three targets pass and the two convergence ratios meet their bounds.

#### Execution Steps

- [ ] Change the analytic displacement expectation and add the pure-B convergence RED in `test_dynamics_basic`.
- [ ] Register `test_particle_stepper`, add its midpoint-field RED, and run both RED commands.
- [ ] Change only `boris_push` position output to the trapezoidal formula and implement the typed candidate interface above.
- [ ] Replace the production `E(x0)`/`boris_push` block with `build_particle_step_candidate`; retain the frozen Stage 1 failure outputs and leave legacy boundary handling unchanged.
- [ ] Update all listed docs and `CHANGELOG.md`, run the three GREEN commands, format, and run `git diff --check`.
- [ ] Commit only Task 2 (P2) files with the message below.

#### Compatibility Note

The signature is preserved, but trajectories change. Treat this as a numerical bug fix, update the scientific contract docs, and record it in release notes when Stage 2 is released.

#### Commit Boundary

`fix: use same-time second-order Boris kinematics`

The commit includes both trapezoidal position and midpoint field sampling so production is not left with a partially second-order claim.

### Task 3: B1 - Add Boundary Event Geometry And Simultaneous-Face Actions

#### Prerequisite

Task 2 (P2) is the branch predecessor; the new boundary API remains additive and does not wire production behavior in this task.

#### Files

- Modify `src/physics/bem_boundary.f90`
- Modify `tests/fortran/test_boundary.f90`
- Update boundary sections in `docs/Algorithms.md`, `docs/Algorithms.en.md`, `docs/ParticleChargeLoop.md`, `docs/ParticleChargeLoop.en.md`, and `SPEC.md`

#### Interfaces

Add public constants, one public type, and additive procedures with these exact declarations:

```fortran
integer(i32), parameter, public :: boundary_event_ok = 0_i32
integer(i32), parameter, public :: boundary_event_invalid_geometry = 1_i32
integer(i32), parameter, public :: boundary_face_x_low  = 1_i32
integer(i32), parameter, public :: boundary_face_x_high = 2_i32
integer(i32), parameter, public :: boundary_face_y_low  = 4_i32
integer(i32), parameter, public :: boundary_face_y_high = 8_i32
integer(i32), parameter, public :: boundary_face_z_low  = 16_i32
integer(i32), parameter, public :: boundary_face_z_high = 32_i32

type, public :: boundary_event_type
  logical :: has_event = .false.
  real(dp) :: fraction = 0.0_dp
  integer(i32) :: face_mask = 0_i32
  integer(i32) :: face_bc(6) = -1_i32
end type

subroutine find_first_boundary_event(cfg, x0, x1, event, status)
  type(sim_config), intent(in) :: cfg
  real(dp), intent(in) :: x0(3), x1(3)
  type(boundary_event_type), intent(out) :: event
  integer(i32), intent(out) :: status
end subroutine find_first_boundary_event

subroutine apply_boundary_event( &
  cfg, event, x, v, alive, escaped, status, &
  q_particle, m_particle, phi_boundary)
  type(sim_config), intent(in) :: cfg
  type(boundary_event_type), intent(in) :: event
  real(dp), intent(inout) :: x(3), v(3)
  logical, intent(inout) :: alive
  logical, intent(out) :: escaped
  integer(i32), intent(out) :: status
  real(dp), intent(in), optional :: q_particle, m_particle, phi_boundary
end subroutine apply_boundary_event
```

`face_bc` stores the configured `bc_open`/`bc_reflect`/`bc_periodic` value for each bit position in x-low, x-high, y-low, y-high, z-low, z-high order. A corner or edge crossing sets every face whose fraction matches the minimum within a `64*epsilon` scale-aware tolerance.

Keep `apply_box_boundary` public with its existing positional and optional arguments. Existing callers continue to compile and retain the legacy overshoot-folding behavior until migrated explicitly.

#### Deterministic Face Rules

- Detect all simultaneous faces before applying any action; never depend on axis loop order.
- Reflect every crossed reflective normal component exactly once.
- Wrap every crossed periodic axis to the opposite face.
- For the default open model, any crossed open face escapes.
- For `potential_barrier`, evaluate each crossed open normal using the event velocity. If any open component is allowed to escape, the particle escapes; otherwise reflect the blocked open components.
- Use `nearest(boundary, inward_sign)` for the post-event representable point instead of a fixed SI epsilon.
- A start exactly on a boundary and moving inward is not an event. A start on a boundary and moving outward is an event at fraction zero and must be resolved without an infinite loop.

#### Invariants

- `event%fraction` is finite and in `[0,1]`.
- Face masks are invariant under axis permutation of equivalent configurations.
- Reflection preserves speed and tangential components.
- Periodic wrapping preserves velocity.
- Detection and action do not inspect or modify mesh charge.
- `apply_box_boundary` remains source compatible.
- Non-finite endpoints, non-positive box spans, invalid face masks, or a non-finite/out-of-range fraction return `boundary_event_invalid_geometry`; none is converted to `has_event=.false.`.
- A `potential_barrier` event containing an open face requires finite `q_particle`, positive finite `m_particle`, and finite `phi_boundary`; missing/invalid values return `boundary_event_invalid_geometry` rather than being interpreted as escape.
- `apply_boundary_event` verifies that every set face bit has the same `face_bc` value captured from `cfg`; mismatched event/config pairs are invalid.

#### RED Tests And Command

Add tests before implementation for:

- a single earliest face crossing;
- simultaneous x-high/y-high reflection at a corner;
- mixed reflect/periodic simultaneous crossing;
- no crossing;
- boundary-start inward and outward cases;
- potential-barrier decision using a supplied crossing velocity;
- face-order invariance.

Also assert the exact face-bit values, `face_bc` contents, and explicit invalid-geometry status for a NaN endpoint, a non-positive box span, a mismatched event/config pair, and a `potential_barrier` open event missing `phi_boundary`.

```bash
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_boundary'
```

Expected RED: compilation fails because the new constants/type/procedures do not exist.

#### GREEN Test And Command

```bash
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_boundary'
```

Expected GREEN: all legacy `apply_box_boundary` cases and all new event/status cases pass.

#### Execution Steps

- [ ] Add the exact constant/type/API compile-time RED plus the deterministic geometry/action cases above.
- [ ] Run the RED command and record the missing-symbol failure.
- [ ] Implement finite geometry validation, simultaneous-face detection, `face_bc` capture, and deterministic action order.
- [ ] Use `nearest` to place reflected/wrapped states one representable value inside the destination face; do not use the legacy fixed `1e-12` offset in the new API.
- [ ] Update the listed docs, run GREEN, format `bem_boundary.f90`, and run `git diff --check`.
- [ ] Commit only Task 3 (B1) files with the message below.

#### Commit Boundary

`feat: add simultaneous box-boundary events`

The new API is not wired into production yet, so this commit is additive and independently green.

### Task 4: B2 - Integrate Earliest Mesh/Box Events And Reflected Remainders

#### Prerequisites

Task 2 (P2) supplies `bem_particle_stepper` and midpoint candidates; Task 3 (B1) supplies the additive boundary event type and procedures. The frozen Stage 1 main/photo failure protocols remain unchanged.

#### Files

- Expand `src/runtime/simulator/bem_particle_stepper.f90`
- Modify `src/runtime/simulator/bem_simulator.f90`
- Modify `src/runtime/simulator/bem_simulator_loop.f90`
- Modify `src/physics/bem_collision.f90`
- Modify `src/particles/bem_injection.f90`
- Modify `src/config/bem_app_config_runtime.f90`
- Modify `tests/fortran/test_particle_stepper.f90`
- Modify `tests/fortran/test_simulator.f90`
- Modify `tests/fortran/test_dynamics_basic.f90`
- Modify `tests/fortran/test_injection_sampling.f90`
- Modify `tests/fortran/test_mpi_hybrid.f90`
- Update `docs/Algorithms.md`, `docs/Algorithms.en.md`, `docs/ParticleChargeLoop.md`, `docs/ParticleChargeLoop.en.md`, and `SPEC.md`
- Update `CHANGELOG.md`

#### Interfaces

Add result/status types internal to `bem_particle_stepper`:

```fortran
integer(i32), parameter :: particle_step_ok = collision_query_ok
integer(i32), parameter :: particle_step_event_limit = 1001_i32
integer(i32), parameter :: particle_step_invalid_event = 1002_i32
integer(i32), parameter :: max_particle_boundary_events = 64_i32

type :: particle_step_result
  real(dp) :: x(3) = 0.0_dp
  real(dp) :: v(3) = 0.0_dp
  logical :: absorbed = .false.
  logical :: escaped_boundary = .false.
  integer(i32) :: elem_idx = -1_i32
  integer(i32) :: status = particle_step_ok
  integer(i32) :: event_count = 0_i32
end type

subroutine advance_particle_step( &
  mesh, sim, field_solver, bfield, x0, v0, q, m, dt, result)
  type(mesh_type), intent(in) :: mesh
  type(sim_config), intent(in) :: sim
  type(field_solver_type), intent(inout) :: field_solver
  real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt
  type(particle_step_result), intent(out) :: result
end subroutine advance_particle_step
```

`bem_particle_stepper` is an internal library module with `private` default; export `particle_step_result`, the three `particle_step_*` status constants, `build_particle_step_candidate`, and `advance_particle_step` for simulator use and focused tests. Keep `max_particle_boundary_events` private.

Add `collision_query_invalid_geometry = 3_i32` to `bem_collision`. Non-finite segment endpoints return that status through the existing optional trailing `status`; omission remains fail-closed. `particle_step_result%status` is one namespace: collision codes `1`, `2`, and `3` propagate unchanged, while event failures use `1001` and `1002`. Do not add a redundant `collision_status` field.

`process_particle_batch` retains this current positional signature and dummy types:

```fortran
module subroutine process_particle_batch( &
  mesh, app, field_solver, pcls_batch, dq_thread, escaped_boundary_flag, absorbed_flag, &
  bfield, batch_idx, mpi_rank, collision_failure_status, &
  collision_failure_particle, collision_failure_step)
  type(mesh_type), intent(in) :: mesh
  type(app_config), intent(in) :: app
  type(field_solver_type), intent(inout) :: field_solver
  type(particles_soa), intent(inout) :: pcls_batch
  real(dp), intent(inout) :: dq_thread(:, :)
  logical, intent(inout) :: escaped_boundary_flag(:), absorbed_flag(:)
  real(dp), intent(in) :: bfield(3)
  integer(i32), intent(in) :: batch_idx, mpi_rank
  integer(i32), intent(out) :: collision_failure_status
  integer(i32), intent(out) :: collision_failure_particle, collision_failure_step
end subroutine process_particle_batch
```

The three existing trailing outputs now carry any nonzero particle-step status, particle index, and outer physical-step index. The internal names may be generalized, but the dummy list, order, and kinds remain unchanged. `run_absorption_insulator` must continue to use the existing failure-presence Allreduce and `mpi_select_lowest_rank_i32_values` before charge commit. Extend the main status-name formatter for `invalid_geometry`, `event_limit`, and `invalid_event`. Because `find_first_hit` is shared, also teach the existing photo finalizers and stop formatter the `collision_query_invalid_geometry` name; keep the photo metadata tuple `[species, ray, bounce, status]`, its pre-field-refresh collective/stop, and every public photo signature unchanged.

#### Event Algorithm

For each remaining portion of the configured step:

1. Build a same-time midpoint candidate for all `dt_remaining`.
2. Find the earliest box crossing on the candidate chord.
3. If a box event exists, restrict the collision query to `x0 -> x_boundary`; otherwise query `x0 -> x_candidate`. This prevents a long periodic overshoot from hitting the image-enumeration guard before it can be split at the periodic face.
4. A mesh hit on the restricted segment, including its endpoint within tolerance, wins and returns absorption.
5. If there is no mesh hit and no box event, accept the candidate state.
6. At a box event, use `v_event = v0 + fraction*(v_candidate-v0)` as the second-order dense event velocity, apply all simultaneous faces, reduce `dt_remaining` by the consumed fraction, and rebuild the midpoint candidate from the transformed event state.
7. For an open potential barrier, evaluate potential at the exact event point and use `v_event`, not the full-step endpoint velocity.
8. Exceeding an internal event cap returns `particle_step_event_limit`; it never silently accepts, escapes, or survives.

Every call is transactional with respect to the caller's particle state: `advance_particle_step` works in local state and `process_particle_batch` copies `result%x/result%v`, sets flags, or deposits only when `result%status == particle_step_ok`.

Tie policy: a mesh hit at the same time as a box face is absorption. This preserves the absorption-only v1 default for surfaces on a box face.

#### Invariants

- Physical time consumed by all pieces sums to the requested `dt` within roundoff.
- A boundary action is applied at the crossing point, never at an overshoot endpoint.
- Reflection and wrap preserve the unconsumed time.
- A reflected-path mesh hit deposits on the base element exactly once.
- `max_step` counts outer physical steps only; `event_count` is a separate numerical guard.
- Collision `image_limit` and `index_range` statuses propagate unchanged into `particle_step_result`.
- Collision `invalid_geometry` and boundary invalid geometry are distinct statuses and both fail closed.
- Non-finite particle state/field candidates, negative `dt`, or non-positive `m` return `particle_step_invalid_event`; `dt == 0` is an identity success.
- All `error stop` for event/collision failure occurs after `!$omp end parallel` and includes batch, rank, particle, outer step, and status.
- Existing no-boundary/no-hit trajectories match Task 2 (P2).

#### RED Tests And Commands

Add deterministic fixtures:

1. **Reflected-path absorption:** box x in `[0,1]`, reflective x-high, plane at x=0.5, particle starts x=0.8 with vx=1 and `dt=1`. It crosses x-high at 0.2, reflects, and must hit the plane during the remaining 0.8. Current code reports no hit.
2. **Mesh before boundary:** the same geometry with a plane before the first face must absorb before any boundary action.
3. **Corner reflection:** simultaneous x/y reflection consumes one event and returns the analytic zero-field final state.
4. **Periodic remainder:** wrap at one face, continue the remaining time, and hit a primary-cell element after wrapping.
5. **Open escape:** default `open_boundary_model="escape"` still kills at the first open event and sets only `escaped_boundary`.
6. **Potential crossing velocity:** an accelerating particle is reflected/escaped according to `v_event`, not `v_candidate`.
7. **Event limit:** a pathological multi-crossing step returns explicit failure and the simulator child process reports context outside OpenMP.
8. **Invalid geometry:** a non-finite collision endpoint returns `collision_query_invalid_geometry`; status omission stops in a child process, while a non-finite boundary fraction, negative `dt`, and non-positive mass map to `particle_step_invalid_event`. A zero `dt` returns the input state with `particle_step_ok`.

Add a focused photo caller case in `test_injection_sampling` that propagates `collision_query_invalid_geometry` through the existing `(status, ray, bounce)` outputs. Keep the existing simulator photo failure case and require its global message to name `invalid_geometry` while preserving `[species, ray, bounce, status]` ordering.

Update `test_simulator` with an end-to-end reflected-path fixture and assert:

- `absorbed=1`, `escaped_boundary=0`;
- deposited charge is exactly `q*w`;
- default insulator absorption behavior is unchanged.

Extend `test_mpi_hybrid`'s existing lowest-rank metadata case with `particle_step_event_limit` as the selected status. This verifies reuse of the frozen Stage 1 helper without introducing a second collective protocol.

```bash
tssrun -p eb -t 0:15:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_particle_stepper'
tssrun -p eb -t 0:15:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_simulator'
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_dynamics_basic'
```

Expected RED: the reflected-remainder and periodic-remainder assertions fail under overshoot folding, and the invalid-geometry status symbol is initially missing.

#### GREEN Tests And Commands

```bash
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_boundary'
tssrun -p eb -t 0:15:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_particle_stepper'
tssrun -p eb -t 0:15:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_simulator'
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_dynamics_basic'
tssrun -p eb -t 0:15:00 --rsc p=1:t=4:c=4 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_FC=mpiifort FPM_ACTION=test FPM_PROFILE=debug FPM_FFLAGS="-fpp -DUSE_MPI -qopenmp" ./build.sh --target test_mpi_hybrid --runner true'
tssrun -p eb -t 0:10:00 --rsc p=2:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && exe=$(find build -type f -path "*/mpiifort_*/test/test_mpi_hybrid" -printf "%T@ %p\n" | sort -nr | head -1 | cut -d" " -f2-); test -n "$exe"; OMP_NUM_THREADS=2 exec "$exe"'
```

Run every controller/allocation sequentially. The MPI build uses one controller process; the runtime uses direct `tssrun` with `p=2`, not nested `mpirun`.

#### Execution Steps

- [ ] Add collision invalid-geometry RED cases to `test_dynamics_basic` and the nine event/result RED fixtures to `test_particle_stepper`/`test_simulator`.
- [ ] Add the generic status selection assertion to `test_mpi_hybrid`, then run the three serial RED commands.
- [ ] Implement the unified status constants/result type and the earliest-event loop without mutating caller-visible state on failure.
- [ ] Replace the legacy production candidate/hit/boundary block with one `advance_particle_step` call per outer step; preserve deterministic lowest particle/step selection and suppress survivor warnings for every failure status.
- [ ] Rename the main internal stop formatter to particle-step terminology, map all five nonzero main codes, and keep the stop after both MPI collectives and before `commit_batch_charge`.
- [ ] Map collision `invalid_geometry` in the existing photo injection/config/simulator formatters without changing their signatures, metadata order, or collective placement.
- [ ] Remove `find_open_boundary_probe` only after the exact crossing point drives potential evaluation in the new stepper.
- [ ] Update the listed docs and `CHANGELOG.md`, run all GREEN commands sequentially, format, and run `git diff --check`.
- [ ] Commit only Task 4 (B2) files with the message below.

#### Commit Boundary

`fix: advance particles through earliest boundary events`

Production wiring, reflected-path collision, status propagation, tests, and docs land atomically. Remove `find_open_boundary_probe` and its helper only in this commit, after the new crossing API is in use.

### Task 5: B3 - Add Physically Scaled Adaptive Particle Substeps

#### Prerequisite

Task 4 (B2) supplies the transactional earliest-event interval engine and the existing main failure envelope `[status, particle, outer_step]`.

#### Files

- Modify `src/core/bem_types.f90`
- Modify `src/config/bem_app_config_types.f90`
- Modify `src/config/app_config_parser/bem_app_config_parser.f90`
- Expand `src/runtime/simulator/bem_particle_stepper.f90`
- Modify `src/runtime/simulator/bem_simulator.f90`
- Modify `src/runtime/simulator/bem_simulator_loop.f90`
- Modify `tests/fortran/test_app_config_parser.f90`
- Modify `tests/fortran/test_particle_stepper.f90`
- Modify `tests/fortran/test_simulator.f90`
- Modify `tests/fortran/test_mpi_hybrid.f90`
- Modify `schemas/beach.schema.json`
- Modify `examples/beach.toml`
- Update `docs/Algorithms.md`, `docs/Algorithms.en.md`, `docs/ParticleChargeLoop.md`, `docs/ParticleChargeLoop.en.md`, `docs/Parameters.md`, `docs/Parameters.en.md`, and `SPEC.md`
- Update `CHANGELOG.md`

#### Interfaces

Add these fields to `sim_config` and parse the same `[sim]` keys:

```fortran
real(dp) :: particle_max_gyro_angle_rad = 0.2_dp
real(dp) :: particle_max_relative_acceleration_change = 0.1_dp
real(dp) :: particle_max_displacement_fraction = 0.25_dp
integer(i32) :: particle_max_substeps = 4096_i32
```

All three real thresholds must be finite and strictly positive; `particle_max_gyro_angle_rad <= 0.5*pi`, `particle_max_relative_acceleration_change <= 1`, and `particle_max_displacement_fraction <= 1`. `particle_max_substeps` must be positive. Use the defaults above in `sim_config`, `init_default_app_config`, the JSON schema, and the example.

Extend the internal result/status contract without changing any frozen simulator entry point:

```fortran
integer(i32), parameter :: particle_step_substep_limit = 1003_i32

type :: particle_step_result
  real(dp) :: x(3) = 0.0_dp
  real(dp) :: v(3) = 0.0_dp
  real(dp) :: hit_time = -1.0_dp
  logical :: absorbed = .false.
  logical :: escaped_boundary = .false.
  integer(i32) :: elem_idx = -1_i32
  integer(i32) :: status = particle_step_ok
  integer(i32) :: event_count = 0_i32
  integer(i32) :: substep_count = 0_i32
end type

subroutine advance_particle_step( &
  mesh, sim, field_solver, bfield, x0, v0, q, m, dt, result, element_length)
  type(mesh_type), intent(in) :: mesh
  type(sim_config), intent(in) :: sim
  type(field_solver_type), intent(inout) :: field_solver
  real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt
  type(particle_step_result), intent(out) :: result
  real(dp), intent(in), optional :: element_length
end subroutine advance_particle_step
```

The trailing optional argument preserves every Task 4 call. Production computes `element_length = minval(mesh%h_elem)` once before entering the OpenMP particle region, validates that it is finite and positive, and passes it by keyword. Focused tests may pass a fixture value directly. Do not scan all elements inside a particle or substep loop.

#### Adaptive Algorithm

For a trial interval `dt_trial`, build the Task 2 midpoint candidate and evaluate Lorentz acceleration at the same-time start and candidate end:

```text
a0 = (q/m) * (E(x0) + cross(v0, B))
a1 = (q/m) * (E(x1) + cross(v1, B))
gyro_angle = abs(q) * norm(B) * dt_trial / m
relative_acceleration_change =
  norm(a1-a0) / max(norm(a0), norm(a1), sqrt(tiny(real_dp)))
displacement_ratio = norm(x1-x0) / element_length
```

Accept a trial only when all three values are within their configured limits. For `q=0`, define gyro angle and relative acceleration change as zero. Cache the candidate endpoint acceleration as the next start acceleration only when the accepted interval has no boundary event and its actual final `(x,v)` equals that candidate endpoint. After reflection or periodic wrap, evaluate acceleration again at the transformed actual state before constructing the next trial; absorption or open escape has no next trial.

Refactor the Task 4 implementation so a private primitive advances an already accepted trial only to its first mesh/box event; it must return the time consumed, transformed/terminal state, and event status without advancing the post-reflection/wrap remainder. Keep `advance_particle_step` as the driver and do not expose this primitive as a new public API.

Start with `dt_trial=dt_remaining`. On rejection, halve `dt_trial` and rebuild the candidate; never accept a criterion violation. When a trial passes all three criteria:

1. If there is no event, consume all `dt_trial`, increment `substep_count`, and keep that accepted size for the next trial, clipped to `dt_remaining`.
2. If absorption or open escape is first, consume only the event fraction of `dt_trial`, set terminal state and `hit_time` as applicable, and return.
3. If reflection or periodic wrap is first, consume only the crossing fraction of `dt_trial`, apply the transform, increment `event_count`, discard the old candidate and its unconsumed remainder, and build a new trial from the transformed actual `(x,v)` for the unchanged remaining physical time. The new trial must pass gyro-angle, acceleration-change, and element-length criteria independently.
4. A positive-time accepted segment, including one truncated by an event, increments `substep_count` once. A fraction-zero boundary transform increments only `event_count`; the existing event cap prevents a zero-time loop.

Accumulate consumed time explicitly and stop immediately on any existing event/collision failure. Candidate endpoint acceleration may be cached only after case 1. Cases 2-4 never reuse the discarded endpoint acceleration.

If another accepted interval would exceed `sim%particle_max_substeps`, return `particle_step_substep_limit` transactionally. A non-finite criterion, a trial time that cannot be halved to a smaller positive representable value, or an invalid element length returns `particle_step_invalid_event`. `hit_time` is elapsed physical time from the start of the outer step to absorption and remains `-1` when there is no hit.

The outer `process_particle_batch` loop is unchanged: `max_step` still counts configured physical steps, not accepted substeps or boundary events. `particle_step_substep_limit` uses the existing OpenMP-local selection, MPI failure-presence Allreduce, lowest-rank metadata selection, and common stop before charge commit. Add no new collective protocol.

#### Invariants

- Gyro rotation per accepted interval is at most the configured angle.
- Relative Lorentz-acceleration change per accepted interval is at most the configured limit.
- Candidate displacement per accepted interval is at most the configured fraction of the mesh-wide minimum `h_elem`.
- Accepted interval durations sum to the requested outer `dt` within roundoff.
- No reflected or wrapped remainder advances until a fresh trial from the transformed actual state passes all three adaptive criteria.
- `max_step` and failure metadata continue to report the outer physical-step index.
- Existing collision, boundary-event, and adaptive failures share one transactional `particle_step_result` envelope.
- A failure never deposits charge, mutates caller-visible particle state, or becomes `survived_max_step`.
- For smooth fields, position/velocity retain second-order convergence; near-surface hit element and `hit_time` converge as the outer `dt -> 0`.

#### RED Tests And Commands

Add:

1. parser/schema cases for all four exact keys, defaults, and invalid zero/non-finite/out-of-range values;
2. a pure-B case whose outer step violates only the gyro-angle limit and therefore requires multiple accepted substeps;
3. a spatially varying electric-field case that violates only the relative-acceleration-change limit;
4. a zero-field high-speed case that violates only the element-displacement limit;
5. a `particle_max_substeps=1` case that returns `particle_step_substep_limit` without changing the input state;
6. a simulator child-process case that reports batch/rank/particle/outer-step context for that status;
7. a curved pure-B trajectory near two competing triangles. Compare outer `dt`, `dt/2`, and `dt/4` against a `dt/64` reference; require the last two refinements to select the reference element, monotonically reduce hit-time error, and give a finest error no greater than `0.1*(dt/4)`;
8. a `test_mpi_hybrid` selection case proving status `1003` reuses the frozen lowest-rank tuple.
9. a reflection and a periodic-wrap fixture where the transformed remainder violates one criterion at the pre-event trial size; require a fresh smaller trial, correct `substep_count/event_count`, and `hit_time` accumulated from only consumed segments.

```bash
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_app_config_parser'
tssrun -p eb -t 0:20:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_particle_stepper'
tssrun -p eb -t 0:15:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_simulator'
```

Expected RED: the parser rejects the new keys, the stepper does not expose adaptive counts/status/hit time, and the old one-chord outer step violates at least one hard fixture bound.

#### GREEN Tests And Commands

```bash
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_app_config_parser'
tssrun -p eb -t 0:20:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_particle_stepper'
tssrun -p eb -t 0:15:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_simulator'
tssrun -p eb -t 0:15:00 --rsc p=1:t=4:c=4 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_FC=mpiifort FPM_ACTION=test FPM_PROFILE=debug FPM_FFLAGS="-fpp -DUSE_MPI -qopenmp" ./build.sh --target test_mpi_hybrid --runner true'
tssrun -p eb -t 0:10:00 --rsc p=2:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && exe=$(find build -type f -path "*/mpiifort_*/test/test_mpi_hybrid" -printf "%T@ %p\n" | sort -nr | head -1 | cut -d" " -f2-); test -n "$exe"; OMP_NUM_THREADS=2 exec "$exe"'
```

Run sequentially. Expected GREEN: all three independent criteria force refinement in their isolated fixtures, the cap fails closed through the existing envelope, and hit element/time converge to the direct reference.

#### Execution Steps

- [ ] Add parser/schema/default RED cases and the nine focused adaptive/status fixtures.
- [ ] Run the three serial RED commands and record the exact old-behavior failures.
- [ ] Add the four validated `sim_config` controls and update schema/example/docs.
- [ ] Refactor the Task 4 engine into a private first-event primitive behind the adaptive outer-step driver; never let it advance a reflected/wrapped remainder internally.
- [ ] Cache the mesh-wide minimum element length outside OpenMP, re-evaluate all three criteria after every boundary transform, and reuse endpoint acceleration only across event-free full accepted trials.
- [ ] Implement deterministic halving, accumulated event/substep counts, `hit_time`, and transactional `particle_step_substep_limit`.
- [ ] Reuse the existing main failure tuple and collective order; extend only its status-name mapping.
- [ ] Run all GREEN commands sequentially, format, and run `git diff --check`.
- [ ] Commit only Task 5 (B3) files with the message below.

#### Commit Boundary

`fix: adapt particle steps to physical scales`

### Task 6: M1 - Compute One Global Reservoir Count And Split It Across Ranks

#### Prerequisite

Use the frozen 13-argument `init_particle_batch_from_config` contract and existing `bem_mpi` helpers from `70c5073`; when following the prescribed commit order, Task 4 (B2) is already green.

#### Files

- Modify `src/core/bem_mpi.F90`
- Modify `src/config/bem_app_config_runtime.f90`
- Modify `tests/fortran/test_reservoir_injection.f90`
- Modify `tests/fortran/test_external_field_velocity_grid.f90`
- Modify `tests/fortran/test_templates_importers_runtime.f90`
- Modify `tests/fortran/test_mpi_hybrid.f90`
- Update reservoir/MPI semantics in `docs/ParticleChargeLoop.md`, `docs/ParticleChargeLoop.en.md`, `docs/Workflow.md`, `docs/Workflow.en.md`, and `SPEC.md`

#### Interfaces

Add MPI wrappers with root defaulting to rank zero and these exact signatures:

```fortran
integer(i32), parameter, public :: mpi_collective_ok = 0_i32
integer(i32), parameter, public :: mpi_collective_invalid_context = 1_i32
integer(i32), parameter, public :: mpi_collective_invalid_root = 2_i32
integer(i32), parameter, public :: mpi_collective_count_range = 3_i32
integer(i32), parameter, public :: mpi_collective_backend_error = 4_i32

subroutine mpi_bcast_i32_array(ctx, values, status, root)
  type(mpi_context), intent(in) :: ctx
  integer(i32), intent(inout) :: values(:)
  integer(i32), intent(out) :: status
  integer(i32), intent(in), optional :: root
end subroutine mpi_bcast_i32_array

subroutine mpi_bcast_real_dp_array(ctx, values, status, root)
  type(mpi_context), intent(in) :: ctx
  real(dp), intent(inout) :: values(:)
  integer(i32), intent(out) :: status
  integer(i32), intent(in), optional :: root
end subroutine mpi_bcast_real_dp_array

subroutine mpi_allreduce_min_i32_scalar(ctx, value, status)
  type(mpi_context), intent(in) :: ctx
  integer(i32), intent(inout) :: value
  integer(i32), intent(out) :: status
end subroutine mpi_allreduce_min_i32_scalar

subroutine mpi_allreduce_max_i32_scalar(ctx, value, status)
  type(mpi_context), intent(in) :: ctx
  integer(i32), intent(inout) :: value
  integer(i32), intent(out) :: status
end subroutine mpi_allreduce_max_i32_scalar
```

All four wrappers initialize `status=mpi_collective_ok`, never `error stop`, range-check array counts before converting to the default MPI integer kind, and map validation/backend failures to the exact status namespace above. The broadcast wrappers validate `0 <= root < ctx%size`, use root `0_i32` when omitted, and are no-ops after validation when MPI is disabled or `ctx%size == 1`. The i32 min/max wrappers are in-place and no-op when MPI is disabled. Keep every pre-Stage-2 `bem_mpi` public signature, including `mpi_select_lowest_rank_i32_values`, unchanged.

Add a private non-stopping count helper and status namespace in `bem_app_config_runtime`:

```fortran
integer(i32), parameter :: reservoir_count_ok = 0_i32
integer(i32), parameter :: reservoir_count_invalid_input = 1_i32
integer(i32), parameter :: reservoir_count_invalid_residual = 2_i32
integer(i32), parameter :: reservoir_count_overflow = 3_i32

subroutine compute_macro_particles_for_species_checked( &
  sim, spec, residual, count, status, vmin_normal, number_density_override, &
  w_particle_override, temperature_k_override, drift_velocity_override, &
  particle_flux_override)
  type(sim_config), intent(in) :: sim
  type(particle_species_spec), intent(in) :: spec
  real(dp), intent(inout) :: residual
  integer(i32), intent(out) :: count
  integer(i32), intent(out) :: status
  real(dp), intent(in), optional :: vmin_normal
  real(dp), intent(in), optional :: number_density_override
  real(dp), intent(in), optional :: w_particle_override
  real(dp), intent(in), optional :: temperature_k_override
  real(dp), intent(in), optional :: drift_velocity_override(3)
  real(dp), intent(in), optional :: particle_flux_override
end subroutine compute_macro_particles_for_species_checked
```

The checked helper always initializes `count=0_i32`, returns one explicit status, and never calls `error stop`. It uses the full `sim%batch_duration`; the old private helper may remain as a compatibility wrapper for internal focused tests, but production M1 counting must call the checked helper.

Preserve the complete frozen signature of `init_particle_batch_from_config`, including the four trailing optional photo collision outputs. Redefine `injection_state%macro_residual` as a global per-species residual on every rank. Reservoir validation must reach one common decision before photo sampling and must not repurpose or alter photo failure metadata.

```fortran
subroutine init_particle_batch_from_config( &
  cfg, batch_idx, pcls, state, mesh, photo_emission_dq, mpi_rank, mpi_size, mpi, &
  collision_failure_status, collision_failure_species, collision_failure_ray, &
  collision_failure_bounce)
  type(app_config), intent(in) :: cfg
  integer(i32), intent(in) :: batch_idx
  type(particles_soa), intent(out) :: pcls
  type(injection_state), intent(inout), optional :: state
  type(mesh_type), intent(in), optional :: mesh
  real(dp), intent(out), optional :: photo_emission_dq(:)
  integer(i32), intent(in), optional :: mpi_rank, mpi_size
  type(mpi_context), intent(in), optional :: mpi
  integer(i32), intent(out), optional :: collision_failure_status
  integer(i32), intent(out), optional :: collision_failure_species
  integer(i32), intent(out), optional :: collision_failure_ray, collision_failure_bounce
end subroutine init_particle_batch_from_config
```

For real MPI execution:

1. Derive `has_enabled_reservoir` from the identical parsed config on every rank. Allocate `global_counts(cfg%n_particle_species)` and set it to zero. Steps 2-5 apply only when that flag is true; a volume-seed/photo-only config skips the residual protocol and continues to allow absent `state`.
2. For an enabled reservoir, set `local_length=-1_i32` when state is absent/unallocated and otherwise use `size(state%macro_residual, kind=i32)`. Copy it to `min_length` and `max_length`, call the new i32 min/max Allreduces, and check each explicit collective status on every rank before entering the next collective. Require `min_length == max_length == cfg%n_particle_species`. Because these are fixed-size scalar collectives, rank-dependent allocation cannot create a count mismatch. All ranks stop through the same collective/size-status branch before any array broadcast.
3. Validate every enabled reservoir's inputs and starting residual locally. Pack fixed `[status, species]` metadata and reuse `mpi_select_lowest_rank_i32_values`; all ranks stop on the selected validation status before root counting or photo sampling.
4. Root calls `compute_macro_particles_for_species_checked` for each enabled reservoir, stopping the loop on its first nonzero `(status,species)`. Broadcast that fixed two-element result with `mpi_bcast_i32_array(..., collective_status)`. Every rank first checks `collective_status`, then returns from the root-only counting phase into the same caller branch and issues the same mapped `error stop` before photo sampling when count status is nonzero.
5. Only after common success, broadcast exactly `cfg%n_particle_species` entries for `global_counts`, check the explicit wrapper status on all ranks, then broadcast the updated global residual vector and check again. No caller proceeds to photo sampling after a nonzero count or collective status.
6. Each rank receives `mpi_split_count(global_count, rank, size)` particles.

For the existing `mpi_rank`/`mpi_size` test interface without a communicator, every synthetic rank deterministically recomputes the same global count from a separately initialized, identical input state, then splits it. Do not sequentially reuse one mutated state across synthetic ranks. Document that these calls require identical starting global residuals.

Remove `batch_duration_scale=1/n_ranks` from the production reservoir path. Volume-seed and photo-ray counts retain their existing global-count splitting.

#### Invariants

- For every species and batch, `sum(local_count(rank)) == global_count`.
- The global count sequence and residual sequence are identical for MPI sizes 1, 2, and 4.
- Every rank holds the same finite residual in `[0,1)` after broadcast.
- No rank enters a residual broadcast until all ranks have confirmed the same exact species-vector length.
- Expected counts below one produce one global particle when the accumulated global budget crosses one; they never produce an MPI-size burst.
- Rank RNG streams and sampled particle states may differ, but counts obey the global contract.
- Target-weight and velocity-grid reservoir branches use the same global counting path.
- Volume-seed/photo-only configurations retain the existing optional-state contract and enter no residual collective.

#### RED Tests And Commands

1. Unit sequence with global expected count `0.25`: four batches must produce global counts `[0,0,0,1]` and residuals `[0.25,0.5,0.75,0]`.
2. Synthetic rank calls for sizes 1, 2, and 4 must sum to that same sequence.
3. Extend `test_mpi_hybrid` with the same reservoir. Use `mpi_allreduce_sum_i32_scalar` to assert the global sequence. The current per-rank residual code produces no particle at batch four for multiple ranks and therefore fails RED.
4. Assert min/max residual across ranks are equal after every batch.
5. Give one rank a residual vector of the wrong length and assert equal min/max-derived size status on every rank without entering an array broadcast.
6. Force each checked-helper failure code directly, and force a root count failure under MPI; require the common mapped stop before photo sampling.
7. Cover both Maxwellian-density and velocity-grid flux count paths.
8. Run volume-seed-only and photo-only MPI fixtures with `state` absent; require unchanged counts and no residual collective failure.

```bash
tssrun -p eb -t 0:15:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_reservoir_injection'
tssrun -p eb -t 0:15:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_external_field_velocity_grid'
tssrun -p eb -t 0:15:00 --rsc p=1:t=4:c=4 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_FC=mpiifort FPM_ACTION=test FPM_PROFILE=debug FPM_FFLAGS="-fpp -DUSE_MPI -qopenmp" ./build.sh --target test_mpi_hybrid --runner true'
```

Expected RED: the serial/synthetic global sequence assertions expose per-rank duration scaling, and the MPI target initially fails to compile because the broadcast wrappers are absent.

#### GREEN Tests And Commands

```bash
tssrun -p eb -t 0:15:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_reservoir_injection'
tssrun -p eb -t 0:15:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_templates_importers_runtime'
tssrun -p eb -t 0:15:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_external_field_velocity_grid'
tssrun -p eb -t 0:15:00 --rsc p=1:t=4:c=4 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_FC=mpiifort FPM_ACTION=test FPM_PROFILE=debug FPM_FFLAGS="-fpp -DUSE_MPI -qopenmp" ./build.sh --target test_mpi_hybrid --runner true'
tssrun -p eb -t 0:10:00 --rsc p=2:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && exe=$(find build -type f -path "*/mpiifort_*/test/test_mpi_hybrid" -printf "%T@ %p\n" | sort -nr | head -1 | cut -d" " -f2-); test -n "$exe"; OMP_NUM_THREADS=2 exec "$exe"'
tssrun -p eb -t 0:10:00 --rsc p=4:t=1:c=1 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && exe=$(find build -type f -path "*/mpiifort_*/test/test_mpi_hybrid" -printf "%T@ %p\n" | sort -nr | head -1 | cut -d" " -f2-); test -n "$exe"; OMP_NUM_THREADS=1 exec "$exe"'
```

Run all allocations sequentially. Compile through one `p=1` controller allocation, then run the MPI executable directly with `tssrun p=2` and `p=4`; do not nest `mpirun` under `tssrun`.

#### Execution Steps

- [ ] Add serial and synthetic rank RED sequences with a fresh identical `injection_state` per synthetic rank.
- [ ] Extend `test_mpi_hybrid` with the global count/residual sequence and cross-rank min/max checks, then run RED.
- [ ] Add the two typed broadcast wrappers plus i32 scalar min/max wrappers, the five-value collective status namespace, and public declarations without changing existing MPI helpers.
- [ ] Add the explicit-status checked count helper, cover all four statuses, and keep production count failures out of `error stop` until their fixed metadata has been broadcast.
- [ ] In production communicator mode, validate rank lengths with min/max, compute reservoir budgets only on root, broadcast fixed failure metadata, then broadcast counts/residuals and split the global count; in synthetic mode, recompute once per independently initialized call.
- [ ] Remove only the production `batch_duration_scale=1/n_ranks` argument; retain the optional helper argument for source compatibility.
- [ ] Update docs, run all GREEN commands sequentially, format, and run `git diff --check`.
- [ ] Commit only Task 6 (M1) files with the message below.

#### Commit Boundary

`fix: distribute global reservoir counts across MPI ranks`

This commit changes runtime count semantics but not checkpoint file semantics. A run must not be advertised as exactly resumable across this intermediate commit; Task 7 (M2) follows immediately before Stage 2 is released.

### Task 7: M2 - Store And Restore One Global Reservoir Residual

#### Prerequisite

Task 6 (M1) has made `injection_state%macro_residual` global on every rank and added typed broadcast plus i32 scalar min/max wrappers with the exact signatures used below.

#### Files

- Modify `src/runtime/bem_restart.f90`
- Modify `tests/fortran/test_restart.f90`
- Modify `tests/fortran/test_mpi_hybrid.f90`
- Update `docs/Workflow.md`, `docs/Workflow.en.md`, `docs/Parameters.md`, `docs/Parameters.en.md`, and `SPEC.md`

#### Interfaces

Keep the public signatures of `load_restart_checkpoint`, `write_macro_residuals_file`, and `restart_macro_residual_path`, including all existing optional `mpi_rank`, `mpi_size`, and `mpi` arguments:

```fortran
subroutine write_macro_residuals_file(out_dir, state, mpi_rank, mpi_size, mpi)
  character(len=*), intent(in) :: out_dir
  type(injection_state), intent(in) :: state
  integer(i32), intent(in), optional :: mpi_rank, mpi_size
  type(mpi_context), intent(in), optional :: mpi
end subroutine write_macro_residuals_file

function restart_macro_residual_path(out_dir, mpi_rank, mpi_size, mpi) result(path)
  character(len=*), intent(in) :: out_dir
  integer(i32), intent(in), optional :: mpi_rank, mpi_size
  type(mpi_context), intent(in), optional :: mpi
  character(len=1024) :: path
end function restart_macro_residual_path

subroutine load_restart_checkpoint( &
  out_dir, mesh, stats, has_restart, state, mpi_rank, mpi_size, mpi, require_checkpoint)
  character(len=*), intent(in) :: out_dir
  type(mesh_type), intent(inout) :: mesh
  type(sim_stats), intent(out) :: stats
  logical, intent(out) :: has_restart
  type(injection_state), intent(inout), optional :: state
  integer(i32), intent(in), optional :: mpi_rank, mpi_size
  type(mpi_context), intent(in), optional :: mpi
  logical, intent(in), optional :: require_checkpoint
end subroutine load_restart_checkpoint
```

Their corrected contract is one root-owned `<out_dir>/macro_residuals.csv` independent of MPI rank. `restart_macro_residual_path` returns that path for serial and MPI calls. Continue writing rank-specific RNG files through the unchanged `restart_rng_state_path` contract.

Use these private I/O status values; private readers/writers return a status and never `error stop` before the production ranks have synchronized it:

```fortran
integer(i32), parameter :: macro_residual_io_ok = 0_i32
integer(i32), parameter :: macro_residual_io_missing = 1_i32
integer(i32), parameter :: macro_residual_io_legacy = 2_i32
integer(i32), parameter :: macro_residual_io_open = 3_i32
integer(i32), parameter :: macro_residual_io_format = 4_i32
integer(i32), parameter :: macro_residual_io_value = 5_i32
integer(i32), parameter :: macro_residual_io_size = 6_i32
```

Implementation rules:

- When the optional `state` is absent, retain the existing no-residual-I/O path. When it is present, derive `local_length=-1_i32` for unallocated state, obtain `min_length` and `max_length` with the Task 6 checked i32 scalar wrappers, and branch commonly on each wrapper status. Require equal positive lengths; map any mismatch to `macro_residual_io_size`. Then validate finiteness and `[0,1)` without stopping and use `mpi_select_lowest_rank_i32_values` so every rank observes the same lowest-rank validation failure before file I/O.
- Root alone attempts the write through a private non-stopping writer. Broadcast its one-element I/O status with `mpi_bcast_i32_array(..., collective_status)`; all ranks check collective status first, then either return success or issue the same mapped I/O error. Non-root performs no file write.
- During load, all ranks continue loading shared summary/charges and their rank-owned RNG file as today. Root alone parses the global residual through a private non-stopping reader, broadcasts its one-element I/O status with an explicit collective status, and broadcasts the residual vector only when both are successful. A nonzero collective or I/O status follows one common error path on all ranks.
- `summary.txt` continues to reject an MPI world-size mismatch, preserving exact same-size resume only.
- Add a private `legacy_restart_macro_residual_path(out_dir, mpi_rank)` helper solely for detection. After summary/charges/RNG have established that a checkpoint exists, if the global residual file is absent, use an i32 Allreduce over each rank's legacy-file presence, set the common status to `macro_residual_io_legacy` when any exists, and otherwise use `macro_residual_io_missing` whenever `state` is present. A directory with no checkpoint still returns `has_restart=.false.` under the existing `require_checkpoint=.false.` contract.
- Serial checkpoints keep the same file name and numeric content.
- If both the global file and stale legacy files exist, load the global file and ignore legacy files; never write or update a legacy path.

The path behavior change is a correctness fix despite a preserved procedure signature. Document it explicitly.

#### Invariants

- Exactly one macro residual vector exists per committed logical state.
- All ranks resume with identical residual values.
- Residual validation remains finite and `[0,1)` for every species.
- Writer and loader reject unequal residual-vector lengths collectively before broadcasting residual values.
- Rank RNG files and `mpi_world_size` enforcement are unchanged.
- Missing or legacy-only MPI residual state cannot silently reset to zero.
- No root-only residual parse/write error strands another rank inside a broadcast.

#### RED Tests And Commands

- Under MPI, set intentionally different in-memory rank residuals and call the writer. The resulting checkpoint must contain the root global value only, and every rank must load that same value.
- Assert only `macro_residuals.csv` is created for a new MPI checkpoint.
- Assert rank-specific RNG files still exist.
- Assert a legacy rank-only residual fixture is rejected with a migration message.
- Allocate a different residual-vector length on one rank and require the common `macro_residual_io_size` path without entering a real-array broadcast.
- Use the existing serial child-process pattern to verify a corrupt root file maps to `macro residual format`; in `test_mpi_hybrid`, verify the root status broadcast and common branch condition before the normal successful load.
- Run an uninterrupted global-count sequence and a save/resume sequence; global counts and final residual must match exactly for the same MPI size.

```bash
tssrun -p eb -t 0:15:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_restart'
tssrun -p eb -t 0:15:00 --rsc p=1:t=4:c=4 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_FC=mpiifort FPM_ACTION=test FPM_PROFILE=debug FPM_FFLAGS="-fpp -DUSE_MPI -qopenmp" ./build.sh --target test_mpi_hybrid --runner true'
tssrun -p eb -t 0:10:00 --rsc p=2:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && exe=$(find build -type f -path "*/mpiifort_*/test/test_mpi_hybrid" -printf "%T@ %p\n" | sort -nr | head -1 | cut -d" " -f2-); test -n "$exe"; OMP_NUM_THREADS=2 exec "$exe"'
```

Expected RED: `test_restart` exposes the old rank-suffixed path, and the MPI test sees distinct per-rank restored residuals instead of the root vector.

#### GREEN Tests And Commands

```bash
tssrun -p eb -t 0:15:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_ACTION=test ./build.sh --target test_restart'
tssrun -p eb -t 0:15:00 --rsc p=1:t=4:c=4 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_FC=mpiifort FPM_ACTION=test FPM_PROFILE=debug FPM_FFLAGS="-fpp -DUSE_MPI -qopenmp" ./build.sh --target test_mpi_hybrid --runner true'
tssrun -p eb -t 0:10:00 --rsc p=2:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && exe=$(find build -type f -path "*/mpiifort_*/test/test_mpi_hybrid" -printf "%T@ %p\n" | sort -nr | head -1 | cut -d" " -f2-); test -n "$exe"; OMP_NUM_THREADS=2 exec "$exe"'
tssrun -p eb -t 0:10:00 --rsc p=4:t=1:c=1 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && exe=$(find build -type f -path "*/mpiifort_*/test/test_mpi_hybrid" -printf "%T@ %p\n" | sort -nr | head -1 | cut -d" " -f2-); test -n "$exe"; OMP_NUM_THREADS=1 exec "$exe"'
```

Run sequentially with the same compile/direct-runtime separation as Task 6 (M1).

#### Execution Steps

- [ ] Replace rank-suffixed expectations in `test_restart` with the global path and add root/non-root synthetic writer, missing-file, and legacy-only RED cases.
- [ ] Change `test_mpi_hybrid` to write distinct in-memory rank values but require root's value on every rank after load; retain rank-specific RNG assertions.
- [ ] Run all RED commands and record both the path-contract and cross-rank residual failures.
- [ ] Make the writer root-owned, make the public path global, implement private legacy detection, and replace root-only I/O stops with the seven-value status protocol above.
- [ ] Allreduce residual-vector length min/max with explicit wrapper status, select any rank-local value failure, then broadcast root I/O status before branching; broadcast residual values only on common success, and assert each common branch in the MPI target.
- [ ] Verify summary MPI-size rejection occurs before residual state is accepted, update docs, and run all GREEN commands.
- [ ] Format, run `git diff --check`, and commit only Task 7 (M2) files with the message below.

#### Commit Boundary

`fix: checkpoint one global reservoir residual`

Atomic generation manifests remain Stage 3 scope. This commit corrects state ownership only.

### Task 8: M3 - Mirror Global Residual Semantics In The Workload Estimator

#### Prerequisite

Tasks 6 (M1) and 7 (M2) define the Fortran global count/residual and `macro_residuals.csv` contracts that this Python parity task must reproduce.

#### Files

- Modify `beach/cli/estimate_fortran_workload.py`
- Modify `tests/python/test_estimate_fortran_workload.py`
- Update `docs/Workflow.md`, `docs/Workflow.en.md`, `docs/PythonPostprocessAPI.md`, and `docs/PythonPostprocessAPI.en.md`

#### Interfaces

Keep the exact callable signature and every existing result key:

```python
def estimate_workload(
    config: dict[str, Any],
    threads: int,
    initial_residuals: list[float] | None = None,
    mpi_ranks: int = 1,
    mpi_rank: int = 0,
    completed_batches: int = 0,
) -> dict[str, Any]:
```

For each reservoir species and batch:

1. Compute the full global expected count without dividing by `mpi_ranks`.
2. Update the single global residual.
3. Floor one global macro count.
4. Split it with `_split_count_for_rank` for `species_per_batch`.

Add result keys `global_species_per_batch: list[list[int]]` and `global_batch_totals: list[int]`. Existing `species_per_batch`, `batch_totals`, `batch_thread_min`, `batch_thread_max`, and `total_particles` remain local to the selected rank. Existing `final_residuals` becomes explicitly global and must be identical for every queried rank.

#### Invariants

- Summing `species_per_batch` over all rank queries reproduces `global_species_per_batch`.
- `global_species_per_batch`, `global_batch_totals`, and `final_residuals` do not depend on `mpi_rank` or `mpi_ranks`.
- `total_particles` remains the selected rank's local workload estimate.
- Resume from `macro_residuals.csv` follows the same sequence as Fortran.

#### RED Tests And Command

- Replace the current per-rank reservoir scaling expectation with a fractional global sequence that distinguishes sizes 1, 2, and 4.
- Query every rank, sum local counts, and compare with hard-coded global counts/residuals.
- Compare uninterrupted and resumed sequences from the same initial global residual.
- Keep volume-seed and photo-ray split tests unchanged.

```bash
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && python -m pytest -q tests/python/test_estimate_fortran_workload.py'
```

Expected RED: the old reservoir test produces per-rank residual/count sequences and the new global result keys are absent.

#### GREEN Test And Command

```bash
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && python -m pytest -q tests/python/test_estimate_fortran_workload.py'
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && ruff check beach/cli/estimate_fortran_workload.py tests/python/test_estimate_fortran_workload.py'
```

Expected GREEN: all estimator tests pass for rank counts 1, 2, and 4, including uninterrupted/resumed residual equality.

#### Execution Steps

- [ ] Replace the per-rank reservoir expectation with hard-coded global sequences and add missing-key/global-invariance RED assertions.
- [ ] Run RED and record the old divided-duration sequence.
- [ ] Build global species counts first, update one residual vector, derive global batch totals, then split each species count for the selected rank.
- [ ] Keep every old local result key and add only the two typed global keys; update CLI output labels and the four exact docs listed above.
- [ ] Run GREEN and `ruff check beach/cli/estimate_fortran_workload.py tests/python/test_estimate_fortran_workload.py` inside the same `p=1` compute allocation.
- [ ] Run `git diff --check` and commit only Task 8 (M3) files with the message below.

#### Commit Boundary

`fix: estimate global reservoir counts before MPI splitting`

## Stage 2 Final Verification

Run only after all focused targets are green. Run the static commands on the control host, then run every payload allocation sequentially:

```bash
git diff --check
make fmt-check-fortran
tssrun -p eb -t 0:30:00 --rsc p=1:t=4:c=4 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 VERSION_MODE=dev make test-l1'
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && python -m pytest -q tests/python/test_estimate_fortran_workload.py'
tssrun -p eb -t 0:10:00 --rsc p=1:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && ruff check .'
tssrun -p eb -t 0:15:00 --rsc p=1:t=4:c=4 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && OMP_NUM_THREADS=2 BEACH_VERSION_MODE=dev FPM_FC=mpiifort FPM_ACTION=test FPM_PROFILE=debug FPM_FFLAGS="-fpp -DUSE_MPI -qopenmp" ./build.sh --target test_mpi_hybrid --runner true'
tssrun -p eb -t 0:10:00 --rsc p=2:t=2:c=2 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && exe=$(find build -type f -path "*/mpiifort_*/test/test_mpi_hybrid" -printf "%T@ %p\n" | sort -nr | head -1 | cut -d" " -f2-); test -n "$exe"; OMP_NUM_THREADS=2 exec "$exe"'
tssrun -p eb -t 0:10:00 --rsc p=4:t=1:c=1 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH && exe=$(find build -type f -path "*/mpiifort_*/test/test_mpi_hybrid" -printf "%T@ %p\n" | sort -nr | head -1 | cut -d" " -f2-); test -n "$exe"; OMP_NUM_THREADS=1 exec "$exe"'
```

The MPI build is a single controller process. The two MPI runs launch the executable directly with Slurm process counts and do not nest another launcher.

Do not run L3 FMM or far-correction targets for these changes unless a failure shows an actual dependency. Record exact Slurm job IDs, compiler/MPI modules, commands, exit codes, and focused assertion summaries in each implementer report.

## Acceptance Checklist

- [ ] `boris_push` remains source compatible and uses same-time trapezoidal position.
- [ ] Production field sampling is at the predicted midpoint.
- [ ] Constant-E and pure-B analytic/convergence contracts pass.
- [ ] Boundary detection returns a simultaneous face mask independent of axis order.
- [ ] Reflection and periodic wrap consume only time up to the event and advance the remainder.
- [ ] Mesh collision is queried only up to the first boundary and includes reflected/wrapped remainders.
- [ ] Default absorption deposits exactly once and existing public simulator APIs remain unchanged.
- [ ] Numerical event failures leave the OpenMP region before stopping.
- [ ] Existing photo collision failures still select `[species, ray, bounce, status]` globally and stop before field refresh or charge commit.
- [ ] `image_limit`, `index_range`, collision `invalid_geometry`, event-limit, and invalid-event statuses all retain deterministic batch/rank/particle/step context.
- [ ] Every accepted adaptive interval satisfies the configured gyro-angle, relative-acceleration-change, and element-displacement bounds.
- [ ] Adaptive cap failure uses the existing transactional failure envelope, while `max_step` remains an outer physical-step count.
- [ ] Curved near-surface trajectories converge to the same hit element and hit time as `dt -> 0`.
- [ ] Reservoir count/residual sequences are independent of MPI world size.
- [ ] New checkpoints contain one root-owned global residual and rank-owned RNG state.
- [ ] Same-size uninterrupted and resumed MPI sequences agree.
- [ ] Python workload estimates reproduce the Fortran global and local sequences.
- [ ] Documentation states the same-time state, earliest-event order, and global residual contracts.

## Explicitly Deferred

- Unified potential snapshot, periodic potential gauge, photo conditional VDF, and Zhao prescribed fields (Stage 6).
- Atomic checkpoint generations, mesh/config fingerprints, and history transaction repair (Stage 3).
- MPI-size-changing non-bitwise resume.

## Self-review

This review is against `70c5073` and must be repeated if implementation starts from a different baseline.

- **仕様 coverage:** P2 が同時刻 Boris と midpoint 場評価、B1 が同時 face の幾何と作用、B2 が mesh/box の最早順序、反射・周期 remainder、collision `invalid_geometry`、Stage 1 の fail-closed 統合を担当する。B3 が gyro angle・加速度変化・要素代表長の全基準、outer-step不変、hit element/time収束を担当する。M1/M2/M3 は global count、明示status、residual broadcast、root checkpoint、Python parity を担当し、Stage 2 scopeを全て実装する。
- **Scope completion:** production adaptive policyはTask 5 (B3)へ入り、Stage 2から延期していない。Explicitly DeferredはStage 3/6とMPI-size-changing resumeだけである。
- **記述具体性:** 後から埋める記述、対象不明の文書、値未定の status、launcher のメタ変数は残っていない。全8 task に exact path、型付き interface、RED failure、GREEN command、checkbox step、単一 commit boundary がある。
- **型/signature consistency:** P1/P2 は9引数 `boris_push` を維持する。B2/B3 は frozen `process_particle_batch` の末尾3個の i32 output と、photo側4 metadataを維持し、全step statusを単一result fieldで伝播する。M1 は `init_particle_batch_from_config` の全13引数を維持し、new collective/count helperだけにmandatory status outputを持たせる。M2 はM1のchecked collective signaturesを再利用し、M3はPython callableの完全なsignatureと既存local result keyを維持する。
- **MPI/KUDPC safety:** 全 fpm/build controller は `p=1` で逐次実行する。MPI compile と direct `tssrun p=2/p=4` runtime を分離し、controller の重複起動と Slurm allocation 内の nested `mpirun` を避ける。
- **Commit isolation:** P1, P2, B1, B2, B3, M1, M2, M3 はそれぞれ独立レビュー可能である。M1 直後だけ exact resume を告知できないことを明記し、Stage 2 release 前に M2 を必須としている。

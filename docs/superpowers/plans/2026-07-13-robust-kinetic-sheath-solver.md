# Robust Kinetic Sheath Solver Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the `kinetic_1d` sheath solve pass the previously failing lunar interface field and complete the 100,000-batch validation without changing its physical equations or acceptance criteria.

**Architecture:** Replace the finite-difference dense Newton system with analytic closure derivatives and an `O(N)` bordered-tridiagonal solve. Globalize the same nonlinear residual using branch-preserving Newton, pseudo-transient recovery, and adaptive interface-field continuation.

**Tech Stack:** Fortran 2008, fpm, existing BEACH test support, KUDPC SysA Slurm.

## Global Constraints

- Keep the kinetic density closures, Bohm condition, monotonic branch, Neumann/Robin boundaries, and `residual_tolerance` unchanged.
- Never return a different physical model as a fallback kinetic solution.
- Classify physical inaccessibility separately from exhausted numerical recovery.
- Run Fortran tests and simulations on a KUDPC compute node, never directly on the login node.
- Do not run multiple `fpm test` controllers concurrently because they share `build/`.

---

### Task 1: Closure Derivative Contract

**Files:**
- Modify: `src/physics/outer_plasma/bem_outer_plasma_kinetic.f90`
- Modify: `tests/fortran/test_outer_plasma_kinetic_core.f90`

**Interfaces:**
- Produces: closure evaluation routines that expose `d_density_d_phi` and, where applicable, `d_density_d_phi_interface`.
- Preserves: existing public density values and status classifications.

- [x] **Step 1: Add finite-difference derivative assertions**

  Add central-difference checks for absorbing electrons, cold drifting ions,
  and both emitted-photoelectron branches at regular interior points.  Use a
  perturbation `h = sqrt(epsilon(1.0_dp))*max(1.0_dp, abs(phi))` and compare
  against the analytic outputs with a scale-aware tolerance.

- [x] **Step 2: Run the core target on SysA and verify RED**

  Submit one controller process that runs
  `FPM_ACTION=test ./build.sh --target test_outer_plasma_kinetic_core`.
  Expected result: compilation fails because the new derivative arguments do
  not yet exist.

- [x] **Step 3: Implement analytic local and interface derivatives**

  Differentiate the exact piecewise formulas in
  `eval_absorbing_maxwellian_density`, `eval_cold_drift_density`, and
  `eval_emitted_maxwellian_density`.  Keep density/status results identical and
  expose interface derivatives only for closures that depend on `phi(1)`.

- [x] **Step 4: Re-run the core target and verify GREEN**

  Expected result: all closure density and derivative assertions pass.

- [x] **Step 5: Commit**

  Commit the test and closure contract together as
  `feat: add analytic kinetic closure derivatives`.

### Task 2: Bordered-Tridiagonal Analytic Newton Step

**Files:**
- Modify: `src/physics/outer_plasma/bem_outer_plasma_kinetic.f90`
- Modify: `tests/fortran/test_outer_plasma_kinetic.f90`

**Interfaces:**
- Consumes: analytic closure derivatives from Task 1.
- Produces: `kinetic_residual_jacobian` and a private bordered-tridiagonal linear solver.
- Removes: production use of `numerical_kinetic_jacobian` and `solve_dense_system`.

- [x] **Step 1: Add Jacobian and difficult-field regression tests**

  Expose a test-only/public Jacobian action routine, compare `J*v` with a
  centered residual difference for a feasible nonlinear profile, and add a
  lunar-options solve at `interface_field=-0.7289832458_dp` seeded from a
  nearby converged field.  Assert `outer_plasma_ok`, monotonicity, and residual
  at or below `1e-8`.

- [x] **Step 2: Run the kinetic target on SysA and verify RED**

  Run `FPM_ACTION=test ./build.sh --target test_outer_plasma_kinetic`.
  Expected result: the Jacobian action API is missing or the difficult-field
  solve reproduces the line-search failure.

- [x] **Step 3: Assemble the analytic structured Jacobian**

  Fill lower, diagonal, upper, and first-column arrays in the same loop that
  evaluates the residual.  Include derivative contributions from every charge
  density and exact derivatives of both boundary rows.

- [x] **Step 4: Implement the bordered solve**

  Solve the tridiagonal base for both `rhs` and the first-column border, then
  apply the scalar rank-one correction.  Reject non-finite pivots and a
  scale-aware near-zero correction denominator.

- [x] **Step 5: Replace the dense Newton path and verify GREEN**

  Remove dense Jacobian allocation and call the structured assembly/solve from
  `solve_outer_plasma_kinetic`.  Re-run the core and kinetic targets; expected
  result is PASS, including the Jacobian action check.

- [x] **Step 6: Commit**

  Commit as `perf: use analytic bordered kinetic Jacobian`.

### Task 3: Pseudo-Transient Recovery

**Files:**
- Modify: `src/physics/outer_plasma/bem_outer_plasma_kinetic.f90`
- Modify: `tests/fortran/test_outer_plasma_kinetic.f90`

**Interfaces:**
- Consumes: structured residual/Jacobian and bordered solver from Task 2.
- Produces: branch-preserving pseudo-transient recovery inside the nonlinear solve.

- [x] **Step 1: Add a recovery-path regression**

  Use the observed lunar failure path: the analytic-Newton-only implementation
  fails while constructing the `-0.70 V/m` source profile, then the recovery
  implementation converges through `-0.7289832458 V/m`.  Assert the original
  residual tolerance and monotonic branch without adding test-only controls.

- [x] **Step 2: Run the kinetic target on SysA and verify RED**

  Expected result: compilation fails because recovery diagnostics and controls
  do not yet exist, or the poor seed ends in numerical failure.

- [x] **Step 3: Add adaptive diagonal regularization**

  When branch-preserving Newton backtracking is exhausted, solve the same
  Jacobian with the pseudo-time diagonal shift.  Increase pseudo-time after an
  accepted residual reduction, reduce it after rejection, and return to
  unshifted Newton when possible.  Check convergence only with the original
  residual.

- [x] **Step 4: Preserve failure classification**

  Propagate closure accessibility failures as physical failures.  Report
  `outer_plasma_numerical_failure` only after the recovery iteration and
  pseudo-time limits are exhausted, with a message naming the exhausted path.

- [x] **Step 5: Re-run kinetic tests and verify GREEN**

  Expected result: existing behavior remains green and the forced recovery
  case converges within tolerance.

- [x] **Step 6: Commit**

  Commit as `fix: recover stalled kinetic sheath solves`.

### Task 4: Adaptive Interface-Field Continuation

**Files:**
- Modify: `src/physics/outer_plasma/bem_outer_plasma_kinetic.f90`
- Modify: `tests/fortran/test_outer_plasma_kinetic.f90`

**Interfaces:**
- Consumes: robust single-field nonlinear solve from Task 3.
- Produces: adaptive continuation from a supplied profile's reconstructed field to the requested field.

- [x] **Step 1: Add large field-change and branch-crossing regressions**

  Seed the lunar solve from a known lower-magnitude field, request a change
  that requires at least one continuation subdivision, and assert final-field
  equality, monotonicity, residual tolerance, and nonzero continuation count.
  Retain the existing sign-crossing regression.

- [x] **Step 2: Run the kinetic target on SysA and verify RED**

  Expected result: the continuation diagnostic is missing or the direct jump
  fails before reaching the target field.

- [x] **Step 3: Implement adaptive field stepping**

  Factor the fixed-field nonlinear solve into a private routine.  Starting from
  the reconstructed previous field, attempt an intermediate target; double a
  successful increment up to the remaining interval and halve a numerically
  failed increment.  Retry from the last converged profile and bound the number
  of halvings.

- [x] **Step 4: Re-run kinetic tests and verify GREEN**

  Expected result: the requested target, not an intermediate state, is returned
  and all old/new kinetic tests pass.

- [x] **Step 5: Commit**

  Commit as `fix: adapt kinetic sheath field continuation`.

### Task 5: Documentation and Tiered Verification

**Files:**
- Modify: `SPEC.md`
- Modify: `docs/USER_GUIDE.md`
- Modify: `docs/superpowers/plans/2026-07-13-robust-kinetic-sheath-solver.md`
- Create: `_handoff/kinetic_sheath_robust_solver_validation.md`

**Interfaces:**
- Consumes: completed solver and tests.
- Produces: user-facing solver semantics and reproducible validation evidence.

- [ ] **Step 1: Document solver and failure semantics**

  State that `kinetic_1d` uses analytic structured Newton with
  pseudo-transient and field continuation, always validates the original
  residual, and never falls back to another physical model.

- [ ] **Step 2: Run static/build checks on the login node**

  Run `git diff --check` and `make check`.  Expected result: both succeed; no
  simulation or Fortran test executable is launched locally.

- [ ] **Step 3: Submit L1 then L2 on SysA**

  Use one `sbatch` controller job at a time with `p=1`, running `make test-l1`
  and then `make test-l2`.  Expected result: both exit zero.

- [ ] **Step 4: Run the 3,000-batch smoke validation**

  Reuse the established lunar-regolith validation input and process count.
  Record job ID, commit, elapsed time, maximum sheath residual, final interface
  field, and any recovery/continuation counts in the handoff report.

- [ ] **Step 5: Run the 100,000-batch validation**

  Submit on an available authorized SysA partition using the previously
  established core count or a modest multiple justified by queue state.  Verify
  that the run passes batch 14,400, completes all 100,000 batches, and never
  accepts a sheath residual above `1e-8`.

- [ ] **Step 6: Review and commit documentation**

  Mark completed checklist entries with observed evidence and commit as
  `docs: record robust kinetic sheath validation`.

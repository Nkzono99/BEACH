# All Phases Completion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the remaining Phase 7 restart/MPI and Phase 8 ownership-interface invariance evidence gaps, then publish a requirement-by-requirement completion audit for physical-redesign Phases 0-9.

**Architecture:** Preserve the existing snapshot as the sole owner of electrostatic and outer-plasma state. Checkpoint schema v3 persists the complete restartable outer profile and scalar solver state while retaining a read-only schema-v2 migration path. Phase 8 invariance is tested through the production simulator so the same physical trajectory gives the same flux ledger and surface charge with two different ownership-interface heights.

**Tech Stack:** Fortran 2008, fpm, MPI, Python 3.11/pytest, Starlight Markdown, KUDPC Slurm.

## Global Constraints

- Run test controllers only on KUDPC compute nodes through `tssrun` or `sbatch`; never run them on `laurel31`.
- Keep `b0=0` and `electrostatic_3d_explicit_orbit` fail-closed contracts unchanged.
- Do not add adaptive particle substeps or a general boundary-event loop.
- Preserve schema-v2 checkpoint reading; write new checkpoints as schema v3.
- Use TDD: each behavior change starts with a focused failing test.
- Do not run multiple `fpm test` commands concurrently in the shared worktree.

---

### Task 1: Complete Outer-Plasma Restart State

**Files:**
- Modify: `tests/fortran/test_restart.f90`
- Modify: `tests/fortran/test_electrostatic_unified.f90`
- Modify: `src/physics/bem_electrostatic_snapshot.f90`
- Modify: `src/runtime/bem_output_writer.f90`
- Modify: `src/runtime/bem_restart.f90`
- Modify: `docs/OutputGuide.md`
- Modify: `docs/OutputGuide.en.md`

**Interfaces:**
- Extends `electrostatic_restart_state_type` with profile field/charge-density arrays and solver scalar state.
- `outer_plasma_profile.csv` schema v3 columns are `point,z_m,potential_V,field_V_m,charge_density_C_m3`.
- `load_kinetic_outer_profile` accepts both the schema-v2 three-column profile and schema-v3 five-column profile.

- [ ] **Step 1: Write failing schema-v3 restart tests**

  Assert that `write_result_files` emits checkpoint schema 3 and five outer-profile columns, then that `load_restart_checkpoint` restores field, charge density, nonlinear iterations/residual, interface field, integrated charge, and individual current components.

- [ ] **Step 2: Run RED**

  Run:

  ```bash
  FPM_ACTION=test ./build.sh --target test_restart
  FPM_ACTION=test ./build.sh --target test_electrostatic_unified
  ```

  Expected: failures because the profile only contains `z,potential` and restart scalars are absent.

- [ ] **Step 3: Implement schema-v3 persistence and v2 migration**

  Extend `electrostatic_diagnostics_type` and `electrostatic_restart_state_type`; export/restore every field from `outer_plasma_state_type` needed by a held outer state. Write/read the new CSV columns and summary keys. Accept schema 2 as a migration input, reconstruct derivative arrays conservatively, and require exact v3 fields for newly written checkpoints.

- [ ] **Step 4: Run GREEN and regression targets**

  Run the two focused targets plus `test_output_writer_io` on one compute-node controller. Expected: all pass.

- [ ] **Step 5: Commit**

  ```bash
  git commit -m "fix: restore complete outer plasma state"
  ```

### Task 2: Prove Kinetic MPI State Determinism

**Files:**
- Modify: `tests/fortran/test_mpi_hybrid.f90`

**Interfaces:**
- Adds `run_mpi_kinetic_outer_test(mpi)` using the production `electrostatic_snapshot_type` kinetic path.
- Compares profile, residual, iteration count, and current diagnostics across ranks after root solve/broadcast and after a held-state refresh.

- [ ] **Step 1: Write the failing MPI contract test**

  Configure physical electron/ion `z_high` reservoirs, resolve `kinetic_outer_plasma_options_type`, initialize the split periodic snapshot with `mpi=mpi`, and assert all-rank equality for profile and scalar diagnostics.

- [ ] **Step 2: Run RED with two ranks**

  Build with `MPI_FC=mpiifx`, then run the executable using direct `tssrun --rsc p=2`. Expected: the new held/restart scalar equality assertions fail until Task 1 state is complete.

- [ ] **Step 3: Correct any missing broadcast fields**

  Keep root-only nonlinear solve ownership. Add only missing scalar/array broadcasts in `solve_kinetic_collective`; do not replicate the solver on ranks.

- [ ] **Step 4: Run GREEN**

  Expected: both ranks report the kinetic test as passing with identical values.

- [ ] **Step 5: Commit**

  ```bash
  git commit -m "test: prove kinetic outer mpi determinism"
  ```

### Task 3: Prove Ownership-Interface Invariance End To End

**Files:**
- Modify: `tests/fortran/test_simulator.f90`
- Modify if the test exposes a defect: `src/runtime/simulator/bem_simulator_loop.f90`
- Modify if the test exposes a defect: `src/physics/outer_plasma/bem_outer_plasma_orbit.f90`

**Interfaces:**
- Adds a simulator test that runs one physical charged-particle case with two `box_max(3)`/`outer_plasma.interface_z` ownership heights while keeping the geometry-derived response plane and far-field solution fixed.
- Compares absorbed/escaped flux, interface gross ledger after accounting for ownership location, final surface charge, energy error, and production-valid diagnostics.

- [ ] **Step 1: Write the failing invariance test**

  Use a negatively charged rough panel and an outward positive test particle whose electrostatic turning point lies between the low and high ownership planes. The low case must hand off and return through the explicit 3D orbit; the high case must turn locally. Both must hit the same absorbing surface and deposit the same charge.

- [ ] **Step 2: Run RED**

  Run `FPM_ACTION=test ./build.sh --target test_simulator` on a compute node. Expected: failure identifies any ownership-dependent ledger, remainder, or surface-charge behavior.

- [ ] **Step 3: Implement the minimum production correction if needed**

  Preserve the shared snapshot and fixed physical response plane. Correct only ownership accounting or return-state composition exposed by the test; do not introduce a second field model.

- [ ] **Step 4: Run GREEN and no-crossing performance regression**

  Run `test_simulator`, `test_particle_stepper`, and the existing release performance target. Expected: invariance passes and the ordinary no-crossing path remains one field evaluation/one collision query.

- [ ] **Step 5: Commit**

  ```bash
  git commit -m "test: close ownership interface invariance"
  ```

### Task 4: Publish The Completion Audit And Release Evidence

**Files:**
- Create: `docs/PhysicsRedesignCompletionAudit.md`
- Create: `docs/PhysicsRedesignCompletionAudit.en.md`
- Modify: `tools/sync_starlight_docs.py`
- Modify: `docs-site/astro.config.mjs`
- Modify: `docs/PhysicsReleaseVerification.md`
- Modify: `docs/PhysicsReleaseVerification.en.md`
- Test: `tests/python/test_docs_sync.py`

**Interfaces:**
- One row per Phase 0-9 completion condition, with source/test/output evidence and explicit scope exclusions.
- Records which earlier review-stage proposals were superseded by the later authoritative physical-redesign plan rather than silently claiming implementation.

- [ ] **Step 1: Add failing docs registration assertions**

  Require the Japanese/English completion-audit pages in the sync map and Starlight Numerics/Developer navigation.

- [ ] **Step 2: Run RED**

  Run `python3.11 -m pytest -q tests/python/test_docs_sync.py` on a compute node. Expected: missing registration failure.

- [ ] **Step 3: Write the evidence matrix and register pages**

  Cite exact test targets, convergence rows, MPI rank count, cache cold/warm constraints, and checkpoint schema. Mark dielectric polarization, PIC space charge, magnetic outer orbits, and statistical return as scope exclusions, not unfinished phases.

- [ ] **Step 4: Run docs checks and Starlight build**

  Run docs sync, `--check`, Ruff, and `npm --prefix docs-site run build` on a compute node. Expected: 43 pages or the current registered count plus two.

- [ ] **Step 5: Commit**

  ```bash
  git commit -m "docs: audit physical redesign phase completion"
  ```

### Task 5: Full Verification And Integration

**Files:**
- Verify: `build/physics-release/manifest.txt`
- Verify: `build/physics-release/convergence.csv`

- [ ] **Step 1: Run static, docs, and focused checks**

  Run `git diff --check`, Ruff, schema checks, docs sync, and Starlight build on a compute node.

- [ ] **Step 2: Run the full physics release gate**

  Submit one two-rank `sbatch` allocation that runs `make test-physics-release` with `MPI_FC=mpiifx` and `BEACH_RELEASE_MPI_RUNNER=srun`.

- [ ] **Step 3: Verify the manifest**

  Require final `status=passed`, passed L3/far/MPI/MPI-cache states, RSS below 8 GiB, and all six convergence categories.

- [ ] **Step 4: Review the complete branch diff**

  Confirm no unresolved Critical/Important correctness issue, no unrelated change, and clean worktree.

- [ ] **Step 5: Merge and push**

  Fast-forward `main`, push `origin/main`, and verify local `main`, remote-tracking `origin/main`, and the pushed hash are identical.

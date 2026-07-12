# Kinetic Sheath Particle Coupling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `kinetic_1d` field, reservoir inflow, escape, and return use one MPI-global potential profile referenced to infinity.

**Architecture:** Refresh the electrostatic snapshot before each batch is injected and pass its read-only outer state into the injection runtime. Add a kinetic-profile mapper beside the existing linear-Debye mapper; select it by an explicit return-model identifier and keep unsupported mixed models fail closed.

**Tech Stack:** Fortran 2008, fpm, MPI, OpenMP, TOML schema/docs, KUDPC SysA Slurm tests.

## Global Constraints

- Preserve existing `linear_debye`, Zhao, and legacy barrier behavior.
- Support only monotonic, collisionless, electrostatic, unmagnetized 1D kinetic profiles.
- Do not silently fall back when an outer state or profile is invalid.
- Keep reservoir counts MPI-global and rank-count invariant.
- Run Fortran payload tests only on a KUDPC compute node.

---

### Task 1: Snapshot-Owned Reservoir Inflow

**Files:**
- Modify: `src/runtime/simulator/bem_simulator.f90`
- Modify: `src/runtime/simulator/bem_simulator_loop.f90`
- Modify: `src/config/bem_app_config_runtime.f90`
- Test: `tests/fortran/test_reservoir_injection.f90`
- Test: `tests/fortran/test_simulator.f90`

**Interfaces:**
- Extend `init_particle_batch_from_config(..., outer_state)` with optional `type(outer_plasma_state_type), intent(in)`.
- Extend `prepare_batch_state(..., outer_state)` and pass `snapshot%outer` after refresh.
- Add `kinetic_inflow_velocity_correction(outer_state, species, vmin_normal, barrier_normal)`.

- [ ] Add a reservoir test whose mesh total charge gives a deliberately different linear-Debye estimate from `outer_state%interface_potential`; assert particle count and interface speed follow the outer state.
- [ ] Run `FPM_ACTION=test ./build.sh --target test_reservoir_injection` on SysA and verify the new assertion fails for the old total-charge path.
- [ ] Move `outer_coupler%refresh(snapshot, mesh, batch_idx)` before `prepare_batch_state` in the batch loop.
- [ ] Pass the ready outer state through `prepare_batch_state` into injection and compute `delta_phi=interface_potential-infinity_potential`.
- [ ] Reuse the existing flux cutoff and barrier energy shift so the accessible infinity flux determines macro count and sampled velocity obeys `v_i^2=v_inf^2-2*q*delta_phi/m`.
- [ ] Reject missing, non-ready, non-finite, or non-kinetic outer state when kinetic transfer is active.
- [ ] Re-run the focused reservoir and simulator tests and verify they pass.

### Task 2: Piecewise-Linear Kinetic Return Mapper

**Files:**
- Modify: `src/physics/outer_plasma/bem_outer_plasma_interface.f90`
- Modify: `src/runtime/simulator/bem_simulator.f90`
- Modify: `src/runtime/simulator/bem_simulator_loop.f90`
- Test: `tests/fortran/test_outer_plasma_interface.f90`

**Interfaces:**
- Add `map_outer_particle_kinetic_profile(state, box_min, box_max, charge, mass, crossing, field_timescale, max_frozen_field_ratio, queue_enabled, outcome)`.
- Add private helpers that validate monotonic profile order, locate the first turning segment, and integrate `dz/sqrt(a+b*z)` exactly on each piecewise-linear segment.

- [ ] Add failing tests for constant-field return time, smooth-profile convergence, escape, periodic tangential shift, and invalid non-monotonic state.
- [ ] Run the focused interface target on SysA and verify the new kinetic cases fail before implementation.
- [ ] Implement endpoint escape from normal-energy conservation and the first monotonic turning point.
- [ ] Implement analytic piecewise-linear flight-time integration with a finite square-root turning endpoint.
- [ ] Preserve existing frozen-field, queue, x/y wrapping, and normal-energy residual contracts.
- [ ] Dispatch to the kinetic mapper only for `outer_plasma.return_model="kinetic_1d_profile_return"`; retain the existing linear mapper unchanged.
- [ ] Re-run the focused interface and simulator targets and verify they pass.

### Task 3: Configuration And Photoelectron Ownership

**Files:**
- Modify: `src/config/bem_physics_config_types.f90`
- Modify: `src/runtime/simulator/bem_simulator_loop.f90`
- Modify: `src/runtime/bem_model_fingerprint.f90` only if the existing return-model feed is insufficient
- Test: `tests/fortran/test_physics_config_types.f90`
- Test: `tests/fortran/test_simulator.f90`
- Test: `tests/fortran/test_mpi_hybrid.f90`

**Interfaces:**
- Accept `kinetic_1d_profile_return` only with `outer.model="kinetic_1d"` and `coupling.particle_transfer_mode="electrostatic_1d_instant_return"`.
- Permit `photoelectron_closure="kinetic_mean"` with this tracked return mode without creating a statistical return particle or ledger deposit.

- [ ] Add failing validation tests for the valid kinetic return combination and rejected linear/kinetic identifier mismatches.
- [ ] Add a ledger test proving kinetic mean diagnostics do not add surface return charge beyond tracked emission and absorption.
- [ ] Update validation to reject simultaneous kinetic transfer with `reservoir_potential_model` or `sheath_injection_model`.
- [ ] Enable outgoing histogram auditing for kinetic-mean tracked photoelectrons without using the `individual_return` closure as a second physical model.
- [ ] Verify model fingerprint already includes `return_model`; add no redundant fingerprint field.
- [ ] Re-run focused config, simulator, and MPI targets on SysA.

### Task 4: Restart, Documentation, And Examples

**Files:**
- Modify: `tests/fortran/test_restart.f90`
- Modify: `tests/fortran/test_mpi_hybrid.f90`
- Modify: `examples/periodic2_kinetic_outer.toml`
- Modify: `docs/Parameters.md`
- Modify: `docs/Parameters.en.md`
- Modify: `docs/Algorithms.md`
- Modify: `SPEC.md`
- Modify: `schemas/beach.schema.json`

**Interfaces:**
- The example selects `kinetic_1d_profile_return`, positive frozen-field controls, and no legacy/Zhao injection correction.
- Existing checkpointed outer profile and model fingerprint remain the restart contract.

- [ ] Add a restart regression comparing the first resumed injection batch with an uninterrupted run using the restored kinetic profile.
- [ ] Add an MPI rank-invariance regression for accessible reservoir counts under a nonzero kinetic barrier.
- [ ] Update schema enum/conditional validation for the new return identifier.
- [ ] Update Japanese and English parameter docs, algorithm flow, SPEC, and the kinetic example.
- [ ] Run schema parsing and focused restart/MPI tests on SysA.

### Task 5: Verification And Integration

**Files:**
- Review all modified files.

- [ ] Run `git diff --check` and `make fmt-check-fortran`.
- [ ] Run `make test-l1` on SysA.
- [ ] Run `make test-l2` on SysA.
- [ ] Run the relevant heavy periodic/FMM and MPI targets without running multiple `fpm test` commands concurrently.
- [ ] Inspect logs for skips, warnings, non-finite diagnostics, and rank-dependent results.
- [ ] Review the final diff against the design specification and ADR 0001.
- [ ] Commit the verified implementation on `feat/periodic-object-detachment`.
- [ ] Merge the feature branch into `main`, rerun the merge gate, and push `main`.

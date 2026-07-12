# Periodic cache cold-build optimization implementation plan

**Goal:** Remove redundant QR work and use assigned OpenMP cores while
preserving the cached periodic operator contract.

**Architecture:** Introduce an internal reusable least-squares factorization,
move one target's construction into a target-local routine, and parallelize
independent targets. Cache locking, publication, fingerprinting, and MPI
broadcast remain unchanged.

**Tech Stack:** Fortran 2008, OpenMP, fpm, Intel MPI on KUDPC SysA.

## Global constraints

- Run Fortran tests and benchmarks only in a SysA Slurm allocation.
- Do not run multiple `fpm test` commands concurrently.
- Preserve cache format and fingerprint.
- Write a failing test before each production change.

### Task 1: Reusable least-squares factorization

**Files:**
- Modify: `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90`
- Modify: `tests/fortran/test_periodic2_operator_cache.f90` or add a focused
  internal test target following existing test exposure patterns.

1. Add a test that prepares one matrix, solves several RHS values, and compares
   every solution with the existing one-shot solver. Require an exposed
   diagnostic counter to show one factorization preparation.
2. Run the focused target on SysA and confirm the assertion fails because the
   reusable preparation interface does not exist.
3. Add a private factorization type containing `q`, `r`, and `col_scale`, plus
   prepare and solve routines. Retain the current regularization and arithmetic.
4. Change target construction to prepare once and solve all proxy RHS values.
5. Run focused periodic operator, infinite operator, and cache tests on SysA.

### Task 2: OpenMP target parallelism

**Files:**
- Modify: `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90`
- Modify: `tests/fortran/test_periodic2_operator_cache.f90`

1. Add a test that builds the same cold operator with one and multiple OpenMP
   threads in separate cache directories and compares operator behavior,
   checksums, build counts, and warm reloads.
2. Run it on SysA and confirm the multi-thread generation diagnostic fails.
3. Extract a target-local builder whose scratch state is not shared, then add an
   OpenMP parallel-do over target indices with static scheduling.
4. Keep cache I/O and MPI operations outside the parallel region.
5. Run the focused test with thread checking enabled, then run L1.

### Task 3: Cold/warm scaling benchmark

**Files:**
- Modify: `tools/periodic_object_validation.py` only if a reproducible
  cold-build-only benchmark entry point is missing.
- Create benchmark artifacts below the validation root, not in the repository.

1. Stage unique cache directories for thread counts 1, 8, 28, 56, and 112.
2. Submit short SysA benchmark jobs to an authorized, available partition.
3. Record cold wall time, cache checksum/fingerprint, warm wall time, speedup,
   core efficiency, and maximum RSS.
4. Select the production thread count at the point where marginal speedup no
   longer justifies additional cores.
5. If cold time remains above five minutes at 112 threads, profile Ewald versus
   least-squares time and write a separate design for fused Ewald evaluation.

### Task 4: Release verification

1. Run L1, L3/heavy, and far-correction opt-in gates sequentially on SysA.
2. Run the two-rank MPI cache concurrency test and a six-rank production-style
   cold/warm smoke.
3. Update performance documentation with measured results and operational
   guidance.
4. Commit the implementation and benchmark documentation before returning to
   the full finite/infinite comparison DAG.

# Periodic cache cold-build optimization implementation plan

**Goal:** Remove redundant QR work and use assigned OpenMP cores while
preserving the cached periodic operator contract.

**Architecture:** Introduce an internal reusable least-squares factorization,
distribute target slices across MPI ranks, parallelize proxy columns within
each target using OpenMP, and assemble the operator with an allreduce. Cache
locking, publication, and fingerprinting remain single-writer operations.

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

### Task 2: Hybrid MPI/OpenMP generation

**Files:**
- Modify: `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90`
- Modify: `tests/fortran/test_periodic2_operator_cache.f90`

1. Extend the MPI cache test to require that cold local target counts sum to the
   global target count, differ by at most one, and are nonzero when ranks do not
   outnumber targets. Require every warm local target count to be zero.
2. Run it with two ranks on SysA and confirm it fails because only rank zero
   currently generates the complete operator.
3. Make cache-hit/miss selection collective while retaining the rank-zero
   filesystem lock and recheck. Distribute targets cyclically, allreduce the
   zero-filled partial operator, and let rank zero publish it.
4. Extract a per-target builder and add an OpenMP parallel-do over proxy columns
   with static scheduling. Make all RHS and coefficient buffers thread-private.
5. Add a serial-versus-multithread cache comparison for operator behavior,
   checksums, build counts, and warm reloads.
6. Run focused serial and MPI tests with thread checking enabled, then run L1.

### Task 3: Cold/warm scaling benchmark

**Files:**
- Modify: `tools/periodic_object_validation.py` only if a reproducible
  cold-build-only benchmark entry point is missing.
- Create benchmark artifacts below the validation root, not in the repository.

1. Stage unique cache directories for `1x1`, `1x112`, `2x112`, `4x112`, and
   `6x112` rank/thread layouts, plus a fixed-core comparison where practical.
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

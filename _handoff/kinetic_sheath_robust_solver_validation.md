# Kinetic Sheath Robust Solver Validation

## Scope

This report records validation of the analytic bordered-tridiagonal Newton,
pseudo-transient recovery, and adaptive interface-field continuation added to
the `kinetic_1d` outer-plasma solver.

The physical equations, monotonic branch, Bohm condition, Neumann/Robin
boundaries, and final residual tolerance are unchanged. No alternate physical
model is used as a fallback.

## Implementation

- Branch: `feature/robust-kinetic-sheath-solver`
- Design commit: `1cf584b`
- Closure derivative commit: `a8221c6`
- Structured Jacobian and recovery commit: `e673aac`
- Adaptive continuation commit: `f56122f`

## Focused SysA Evidence

| Job | Target | Result |
|---|---|---|
| `8113020` | pre-change kinetic core + solver baseline | 7/7 core and 9/9 solver tests passed |
| `8113022` | closure derivative RED | expected compile failure on missing derivative API |
| `8113028` | closure derivative GREEN | 10 tests, 24 assertions passed |
| `8113031` | analytic Jacobian RED | expected compile failure on missing Jacobian action API |
| `8113034` | analytic Jacobian without recovery | Jacobian action passed; lunar `-0.70 V/m` source reproduced Newton line-search failure |
| `8113040` | pseudo-transient GREEN | 10 core and 11 solver tests passed; `-0.70` and `-0.7289832458 V/m` converged |
| `8113041` | continuation RED | expected compile failure on missing continuation diagnostic |
| `8113047` | continuation GREEN | 10 core and 11 solver tests passed, 41 solver assertions |
| `8113048` | `make check` | debug OpenMP build passed on SysA |
| `8113052` | documentation mirror regression | passed after canonical/plugin synchronization |
| `8113053` | L1 | passed: Python 555 passed / 36 skipped and all light Fortran targets passed |
| `8113065` | L2 | passed: L1 plus C field-kernel and periodic zero-mode contracts |
| `8113081` | release MPI build | passed with Intel 2023.2 MPI/OpenMP |

The difficult-field regression first solves the lunar ambient state at zero
field, reaches `-0.70 V/m` through multiple continuation steps, then reaches
the former runtime failure field `-0.72898324579369622 V/m`. The final state is
monotonic and its original residual is at most `1e-8`.

## Runtime Chain

- Validation root: `/LARGE1/gr20001/b36291/codex-tmp/beach-robust-kinetic-validation-20260713`
- Binary SHA-256: `054736111fb133b6d84175d8ff7c5fac405a332887a2b3f3953e0671b85e1439`
- Profiling diagnostic: jobs `8113082` / `8113083` were cancelled because
  `BEACH_PROFILE=1` changed the runtime comparison condition; the partial smoke
  output is retained as `output/smoke.profiled-cancelled-8113082`.
- 3,000-batch smoke: job `8113092`, `gr20001a`, `p=6:t=112:c=112`, one-hour limit
- 100,000-batch full: job `8113093`, `gr20001a`, `p=12:t=112:c=112`, ten-hour limit
- Dependency: `afterok:8113092`; a failed smoke prevents the full run from starting

## Pending Runtime Evidence

- 3,000-batch smoke result: pending
- 100,000-batch lunar-regolith result: pending

## Runtime Performance Incident

The first release binary for this branch was built through the Makefile default
`MPI_FC=mpiifort`, and its launch specified `OMP_NUM_THREADS=112` without the
thread-placement variables used by the prior production job.  It completed the
3,000-batch smoke in `00:50:26`, compared with `00:17:37` for the previous
`mpiifx` binary with `OMP_PROC_BIND=spread` and `OMP_PLACES=cores`, despite
producing identical particle counts, charge history, and sheath state to
numerical precision.  The dependent 100,000-batch job `8113093` reached only
about 3,700 batches in 4.6 hours and was cancelled with its partial output
retained as `output/full.mpiifort-cancelled-8113093`.

The old production binary is byte-for-byte identical to the retained `mpiifx`
build, while the slow binary was built with classic `ifort 2021.10`.  Because
compiler and affinity differed together, follow-up jobs isolate both effects:
job `8113792` repeats the 3,000-batch current-source `mpiifx` smoke with the old
affinity, and job `8113793` compares current-source `mpiifort` and `mpiifx`
binaries under identical affinity on a 300-batch fixture.

Job `8113793` completed successfully.  With six MPI ranks, 112 OpenMP threads
per rank, `OMP_PROC_BIND=spread`, and `OMP_PLACES=cores`, elapsed times were
`364.28 s` for `mpiifort` and `131.84 s` for `mpiifx`, a `2.76x` speedup.  The
profile attributes the largest compiler difference to `particle_batch`
(`301.60 s` versus `77.57 s` rank maximum), not to the kinetic solve in
`field_refresh`.  The two summaries differ only in build/cache provenance.
Job `8113816` is the matched `mpiifx` no-affinity measurement.

The build default and release-gate fallback now select `mpiifx 2023.2.4`;
`mpiifort` remains an explicit override.  Production documentation now requires
`OMP_PROC_BIND=spread` and `OMP_PLACES=cores` for the 112-thread SysA layout.
Interface-field continuation was also corrected to try the requested field
first and halve only after a numerical failure, avoiding unconditional
intermediate solves.

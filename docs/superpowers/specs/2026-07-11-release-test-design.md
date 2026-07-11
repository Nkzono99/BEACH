# Release Test Design

## Goal

Keep the existing cumulative L0-L3 developer interface while making the HPC
release gate shorter, measurable per Fortran target, and explicit about the
difference between correctness tests, diagnostics, and performance benchmarks.

## Test tiers

- `test-l0` through `test-l3` remain cumulative and backward compatible.
- Portable CI continues to own `test-l2`.
- The HPC release gate runs only the L1 tests that emit required convergence
  rows, the L3-heavy tests, far-correction correctness tests, and MPI tests.
- A program without numerical assertions is a diagnostic and is not a release
  correctness gate.
- Runtime comparisons execute as non-installable fpm examples under the release
  profile; debug-profile correctness tests do not claim a performance contract.

## Timing artifacts

Every explicitly selected Fortran target is executed sequentially because fpm
targets share `build/`. A reusable shell runner records target name, profile,
status, and elapsed seconds in a CSV when a timing path is supplied. The release
manifest points to the per-gate CSV files; target timing does not introduce a
machine-dependent hard failure threshold.

## Far-correction fixtures

`test_periodic2_flat_oracle_diag` remains available through a diagnostic target
but is removed from the release gate until it has a justified total-field error
contract. The cache test uses a unique temporary cache directory so it can
perform one cold generation, one warm load, and one corrupt-cache rebuild
without a probe generation solely to discover the cache path.

## Compatibility

Existing public simulation APIs and the cumulative `make test-l*` commands do
not change. Documentation and Python contract tests are updated together with
new Makefile entry points and manifest keys.

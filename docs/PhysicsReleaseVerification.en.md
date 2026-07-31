# Physics release verification

Quantitative use of the periodic2 model requires both portable CI and the HPC release gate:

```bash
make test-l2
make test-physics-release
```

The HPC gate sequentially runs the three L1 tests needed for convergence rows, L3-heavy tests,
far-correction correctness, MPI ledger, and MPI periodic-cache concurrency. It does not repeat the full L2 suite
owned by portable CI.
Its manifest records elapsed time and GNU time peak RSS for every gate and enforces an 8 GiB default budget,
overridable with `BEACH_RELEASE_MAX_RSS_KB`. It also regenerates `convergence.csv` from standardized test output
and fails when a required convergence category is missing.
Per-target Fortran timings are written to `test_l3-target-timings.csv` and
`far_correction-target-timings.csv`.

Run periodic far-correction correctness checks with `make test-fortran-far-correction`.
Runtime comparison is separated from debug correctness and runs under the release profile with
`make test-fortran-benchmark`.

The reference values and interpretation are maintained in `docs/PhysicsReleaseVerification.md`.
They are small correctness fixtures and do not replace mesh, time-step, FMM, outer-grid, and sampling studies
for a production simulation.

# Physics release verification

Quantitative use of the new periodic2 and outer-plasma models requires both portable CI and the HPC release gate:

```bash
make test-l2
make test-physics-release
```

The HPC gate runs L3, far-correction, MPI ledger, and MPI periodic-cache concurrency sequentially.
Its manifest records elapsed time and GNU time peak RSS for every gate and enforces an 8 GiB default budget,
overridable with `BEACH_RELEASE_MAX_RSS_KB`. It also regenerates `convergence.csv` from standardized test output
and fails when a required convergence category is missing.

The reference values and interpretation are maintained in `docs/PhysicsReleaseVerification.md`.
They are small correctness fixtures and do not replace mesh, time-step, FMM, outer-grid, and sampling studies
for a production simulation.

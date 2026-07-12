# Periodic cache cold-build optimization design

## Goal

Reduce the `cached_kneq0` cold-build time without changing the generated
operator's mathematical definition, cache identity, or warm-load behavior.
The reference case is the archived regolith geometry with FMM order 4,
64 target anchor nodes, 280 proxy points, and 840 check points. Its observed
cold-build wall time is 31 minutes 24 seconds on SysA while six MPI ranks and
112 OpenMP threads per rank are allocated.

## Current bottlenecks

For each target node, `field_matrix` is constant across all proxy sources.
The current implementation nevertheless allocates and factors the same
regularized `2520 x 34` least-squares matrix once per proxy point. This causes
280 identical QR factorizations per target, or 17,920 factorizations for the
reference case.

Only MPI rank zero builds the operator. The other ranks wait for the final
broadcast, and the build loop has no OpenMP region. Increasing MPI ranks or
`OMP_NUM_THREADS` therefore does not accelerate a cache miss.

## Design

### Reusable least-squares factorization

Split the current solver into two operations:

1. Prepare a regularized, column-scaled QR factorization from `field_matrix`.
2. Solve any number of right-hand sides using that immutable factorization.

The factorization is prepared once per target and reused for all 280 proxy
right-hand sides. The same modified Gram-Schmidt implementation, ridge term,
column scaling, and triangular solve are retained. This deliberately avoids a
LAPACK migration in the first optimization stage so numerical drift is limited
to roundoff from unchanged operations.

### OpenMP target parallelism

Each target produces a disjoint `periodic_root_operator(:,:,target_idx)`.
Parallelize the outer target loop with OpenMP. All scratch arrays used inside a
target iteration are private or allocated inside a target-local helper. The
plan and Ewald tables are read-only during generation. Cache publication and
MPI broadcast remain serial after the parallel region.

MPI distribution is not included in this stage. Rank zero remains the sole
cache producer under the existing filesystem lock, which preserves the
single-writer and MPI-concurrency contracts. OpenMP consumes the cores already
assigned to rank zero; production jobs need not increase MPI ranks for cold
generation.

## Numerical contract

- FMM order, proxy/check points, target nodes, Ewald layers, ridge, and QR
  tolerance do not change.
- Cache format, generator tag, fingerprint, and checksum do not change.
- Serial and multi-thread cold builds must produce operators equal within a
  tight floating-point tolerance and must agree at field/potential probes.
- Cold build count remains one on the root rank; warm load remains a cache hit.
- A failed allocation or invalid factorization continues to fail closed.

Byte-for-byte cache identity is desirable but is not required for OpenMP
execution because independent target slices may finish in a different order;
each slice's arithmetic order remains deterministic. The checksum must validate
the published result.

## Verification and benchmark

Add a focused factorization-reuse unit test that counts one preparation and
multiple solves and compares them with the legacy one-shot solver. Extend the
periodic cache test to compare cold and warm fields and cache metadata.

Run focused Fortran targets and L1 on a SysA compute node. Then benchmark an
isolated cold build of the archived reference configuration with 1, 8, 28,
56, and 112 OpenMP threads on one MPI rank where the queue permits. Report wall
time, speedup, and core efficiency. Finally repeat the normal six-rank launch
to verify MPI broadcast and production compatibility.

## Deferred optimization

Fusing field and potential Ewald evaluation, replacing the QR implementation
with LAPACK, and distributing targets across MPI ranks are deferred until this
stage is measured. They change more arithmetic or communication behavior and
should only be added if QR reuse plus OpenMP does not reduce the cold build to
an operationally acceptable duration.

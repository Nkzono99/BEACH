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

### Hybrid MPI and OpenMP parallelism

Each target produces a disjoint `periodic_root_operator(:,:,target_idx)`.
Distribute target indices cyclically across MPI ranks. Each rank initializes
only its assigned operator slices and leaves the other slices zero. An
`MPI_Allreduce(SUM)` assembles the complete operator on every rank before rank
zero publishes it.

Within each assigned target, prepare the immutable QR factorization once and
parallelize the 280 independent proxy-source columns with OpenMP static
scheduling. Every thread owns its field/potential RHS, coefficients, and local
operator column. The plan, Ewald tables, check points, field matrix, and QR
factorization are read-only in the parallel region. This decomposition has
enough work for 112 threads even when six MPI ranks divide 64 target nodes;
target-only OpenMP parallelism would expose only about 10 or 11 tasks per rank.

Rank zero alone acquires and holds the existing filesystem lock. It first
rechecks the cache after acquiring the lock and broadcasts whether generation
is required. On a miss all ranks generate their assigned slices collectively;
rank zero then publishes and releases the lock. On a hit rank zero broadcasts
the loaded operator without entering generation. This preserves the
single-writer contract across concurrent jobs.

## Numerical contract

- FMM order, proxy/check points, target nodes, Ewald layers, ridge, and QR
  tolerance do not change.
- Cache format, generator tag, fingerprint, and checksum do not change.
- Serial and multi-thread cold builds must produce operators equal within a
  tight floating-point tolerance and must agree at field/potential probes.
- Cold build count remains one on the root rank; warm load remains a cache hit.
- On a cold build, the sum of per-rank local target counts equals the global
  target count and work is balanced to within one target. A warm load reports
  zero locally generated targets on every rank.
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
isolated cold build of the archived reference configuration across a bounded
hybrid matrix such as `1x1`, `1x112`, `2x112`, `4x112`, and `6x112`, subject to
partition availability. Report wall time, speedup, and core efficiency. Include
a fixed-core comparison where practical to distinguish MPI scaling from total
core-count scaling.

## Deferred optimization

Fusing field and potential Ewald evaluation and replacing the QR implementation
with LAPACK are deferred until this stage is measured. They change more
arithmetic behavior and should only be added if QR reuse plus hybrid MPI/OpenMP
does not reduce the cold build to an operationally acceptable duration.

---
name: BEACH performance architecture overview
description: Key hotspots, OpenMP patterns, and field solver modes (direct/treecode/FMM) in the Fortran core
type: project
---

BEACH Fortran core has three field solver modes: direct O(N), treecode O(N log N), and FMM O(N).

**Why:** Field solver choice is the primary determinant of per-particle cost; direct mode is the default fallback.
**How to apply:** Performance review of field evaluation must consider which solver mode is active. The treecode/FMM paths have their own OpenMP-parallel upward/downward passes.

Key performance-critical code paths (2026-03-25 audit):
- **Particle batch loop**: `bem_simulator_loop.f90` lines 135-169, OpenMP parallel do over particles
- **Field evaluation**: Called per-particle per-step inside the parallel region. In FMM mode, `eval_point` does tree traversal (serial per call). In direct mode, `electric_field_at` uses `!$omp simd`.
- **Collision detection**: `find_first_hit` called per-particle per-step. Uses DDA grid or linear scan. No OpenMP inside.
- **FMM state update**: `bem_coulomb_fmm_state_ops.f90` has P2M/M2M/M2L/L2L passes, all with `!$omp do schedule(static)` inside a single `!$omp parallel` region.
- **Tree refresh (treecode)**: `bem_field_solver_tree.f90` has level-by-level `!$omp parallel do`.

Data layout: SoA for particles (`x(3,n)`, `v(3,n)`, `q(n)`, `m(n)`). Mesh uses separate `center_x(:)`, `center_y(:)`, `center_z(:)` for SIMD-friendly access. Vertex data is `(3,nelem)` column-major.

title: Physics redesign completion audit

Lang: [日本語](PhysicsRedesignCompletionAudit.md) | [English](PhysicsRedesignCompletionAudit.en.md)

# Physics redesign completion audit

This page traces implementation and verification evidence for the periodic2, P0-panel, outer-plasma, and particle
handoff redesign. The design source of truth is
`_handoff/BEACH_periodic2_sheath_panel_redesign_implementation_plan.md`. Because `_handoff` is not published, the
operational conclusions are preserved here.

## Status convention

- **Complete** means that the production path, fail-closed applicability rules, and relevant test tier exist.
- **Out of scope** denotes physics excluded by design, not an unfinished phase.
- Quantitative use still requires case-specific convergence checks from the [validation guide](ValidationGuide.html).
- Evidence uses analytic solutions, conservation, convergence, and MPI/restart determinism, not solver agreement alone.

## Phase 0–9 evidence matrix

| Phase | Status | Production owner | Main contract test / gate |
| --- | --- | --- | --- |
| 0 contracts, ledger, migration | Complete | typed physics config, `bem_charge_ledger`, model fingerprints, versioned restart | `test_charge_ledger`, `test_model_fingerprint`, `test_restart`, `test_simulator`, Python result tests |
| 1 direct P0 panel and surface side | Complete | `bem_panel_geometry`, `bem_panel_kernel`, `bem_panel_self_terms`, per-element vacuum side | `test_panel_kernel` analytic/oracle, continuity, jump, and self tests; `test_panel_moments` transform, permutation, and scaling tests |
| 2 separated `k0` and linear outer state | Complete | `bem_periodic_zero_mode_*`, panel-spectral `k!=0` reference, electrostatic snapshot, linear Debye state | `test_periodic_zero_mode`, `test_outer_plasma_linear`, `test_electrostatic_snapshot`, far-correction oracle |
| 3 ambient/interface transaction | Complete | earliest mesh/box event, typed crossings and outcomes, reservoir flux map, outer coupler | `test_outer_plasma_interface`, `test_interface_particle_buffer`, `test_particle_stepper`, `test_simulator`, MPI global reservoir count |
| 4 photoelectron transfer | Complete | individual-return histogram, previous-batch ownership, signed charge ledger | `test_outer_plasma_photoelectron`, simulator emission/return transaction, MPI global histogram |
| 5 full panel FMM | Complete | panel-aware topology, near subtract/add, exact panel P2M, public C kernel contract | `test_panel_geometry_near`, `test_panel_near_correction`, `test_coulomb_fmm_core_panel`, `test_dynamics_panel_fmm`, L2 C/kernel contract |
| 6 production infinite-periodic operator | Complete | versioned `K_periodic,k!=0 - K_shell` cache, checksums/fingerprints, MPI target distribution, rank-local OpenMP, collective assembly | `test_periodic2_operator_cache`, `test_periodic2_infinite_operator`, `test_periodic2_cached_snapshot`, two-rank cache MPI test |
| 7 nonlinear kinetic outer sheath | Complete | stretched-grid Poisson/Robin solver, branch/applicability status, root-only collective solve, restartable profile | kinetic core/coupled/runtime tests, `test_restart`, two-rank kinetic snapshot broadcast |
| 8 unified outer domain | Complete | local mean/accessibility, screened nonzero tail, unified zero mode, explicit electrostatic 3D outer orbit | local-mean, nonzero-tail, unified-snapshot, outer-orbit, and simulator interface-height surface-charge tests |
| 9 production promotion | Complete | portable L2 plus HPC L3/far/MPI gates, RSS budget, convergence artifact, bilingual docs/schema/example sync | `make test-physics-release` manifest and `convergence.csv`, Starlight build, docs sync check |

The original implementation plan defines Phases 0–8 and defines production release criteria in section 9. This audit
calls that final promotion Phase 9; it does not introduce another physics phase.

## Restart contract

New checkpoints use schema 3 and restore the complete kinetic/unified outer state:

- full `z`, potential, electric-field, and charge-density profiles;
- solver status, nonlinear iteration count, and residual;
- interface/infinity potentials, interface field, Debye length, and linearity metric;
- integrated charge per area and electron, ion, photoelectron, and total currents.

A schema-2 three-column profile remains a read-only migration input. It is an initial guess, not a complete held
state, so the next kinetic refresh is mandatory. Missing schema-3 state fails closed.

## Interface and conservation contract

`outer_plasma.interface_z == sim.box_max[3]` is a particle ownership plane, not a field-domain truncation. The
production simulator runs one bound orbit at two ownership heights. The low plane hands the particle to the explicit
outer orbit and receives it back; the high plane turns it inside the local domain. Both runs must have identical
absorbed/escaped counts and final total surface charge, while each ownership ledger and the total charge residual close.

The pusher uses the same-time `(x^n,v^n)` midpoint Boris contract. Pure-magnetic speed conservation, time reversal,
and second-order `dt` convergence are verified. The ordinary in-box hot path does not enter an extra event loop;
multiple box crossings caused by an excessive `dt` fail closed.

## Release evidence

The supported entry points are:

```bash
make test-l2
make test-physics-release
```

The HPC gate sequentially runs the L1 convergence subset, L3-heavy, far-correction correctness, two-rank MPI hybrid,
and two-rank MPI cache-concurrency tests without repeating the portable L2 suite. A
release requires final and per-gate `status=passed`, peak RSS below the default 8 GiB budget, and all six convergence
categories: `boris_dt`, `panel_fmm_order`, `rough_panel_mesh`, `rough_outer_grid`, `rough_accessibility`, and
`outer_orbit_dt`.

Run metadata is written to `build/physics-release/manifest.txt`; numerical rows are written to
`build/physics-release/convergence.csv`, and per-target timings to the sibling `*-target-timings.csv` files.
These artifacts are regenerated for each release rather than committed as
fixed repository data.

## Superseded review stages

`docs/superpowers/specs/2026-07-10-review-remediation-design.md` remains useful prior review context, but conflicting
parts were superseded by the later physical-redesign plan:

- the point-source final baseline was replaced by the continuous triangle-P0 contract;
- general charged-wall closure was replaced by an explicit zero mode and outer provider;
- analytic periodic M2L as the sole runtime backend was replaced by the versioned cached operator;
- Zhao as the final outer closure was separated as a legacy injection model.

Non-conflicting requirements such as literal output paths, strict histories, MPI-global counts, the common
field/potential snapshot, and root-only cache I/O were retained. This page does not claim completion of every
Stage 3–6 proposal, including whole-checkpoint temporary-write/atomic publication,
generation-directory/current-manifest checkpoints, or a Zhao prescribed profile. Those stages are a separate roadmap
from physical Phases 0–9.

## Explicitly out of scope

- a self-consistent dielectric-polarization boundary condition;
- 3D volume PIC/Poisson space charge and nonlinear lateral outer coupling;
- magnetic outer orbits; the explicit outer path is electrostatic and requires `b0=0`;
- photoelectron statistical-return landing distributions, delays, and persistent queues;
- panel-aware treecode and a general softened continuous-panel kernel;
- substrate electrical closures beyond `E_bottom=0`.

Cases needing these features fail validation or applicability checks instead of silently selecting an approximation.

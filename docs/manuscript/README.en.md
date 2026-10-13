# BEACH physical models and numerical algorithms

Lang: [English](README.en.md) | [日本語](README.md)

This is a Japanese LaTeX manuscript in paper/technical-book form covering BEACH's models, algorithms,
and related research. Start with [beach_models.tex](beach_models.tex). Chapter sources are under
`sections/`, and the bibliography is [references.bib](references.bib). It includes Japanese and English
abstracts, a table of contents, equations, a computation-cycle diagram, and references.

## Contents

1. Surface charging, airless bodies, and related research
2. Model partition, governing equations, and batches
3. P0 triangular sources, analytic integration, self terms, and one-sided fields
4. Direct, treecode, Cartesian FMM, and two-periodic Ewald/zero-mode fields
5. Same-time Boris updates, collisions, box boundaries, and periodic images
6. Boundary inflow, velocity distributions, and ray-based photoelectron emission
7. Charge budgets, return closure, floating conductors, and batch stability
8. Stationary outer-sheath closure (zero-current root, ion fluid, upstream band and field floor, top-face mapping, outer-root refresh)
9. Verification, convergence, Monte Carlo errors, and reproducibility
10. Force and detachment analysis, methods and limits to report, and manuscript statements

## Writing policy

The text assumes general knowledge of plasma physics and numerical methods and keeps the models and algorithms
self-contained. It does not name input keys, source paths, output files, tests, commits, or computing environments.
Origins of the methods are given by citations. The correspondence to the implementation is kept in this table.

| Chapter | Main implementation and documents |
| --- | --- |
| 2 Computation cycle | `src/runtime/simulator/bem_simulator_loop.f90`, `docs/Algorithms.md` |
| 3 Triangle fields | `src/physics/field_solver/panel/bem_panel_kernel.f90`, `src/physics/field_solver/panel/bem_panel_self_terms.f90`, `docs/DirectSolver.md` |
| 4 Treecode and FMM | `docs/Treecode.md`, `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90`, `docs/FMMCore.md` |
| 4 Periodic far field | `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90`, `docs/PeriodicFarCorrection.md` |
| 4 Zero mode | `src/physics/field_solver/periodic/bem_periodic_zero_mode_plan.f90`, `src/physics/field_solver/periodic/bem_periodic_zero_mode_eval.f90`, `docs/PeriodicElectrostatics.md` |
| 5 Particle integration and collisions | `src/physics/bem_pusher.f90`, `src/physics/bem_collision.f90`, `src/runtime/simulator/bem_particle_stepper.f90`, `docs/BorisPusher.md`, `docs/ParticleEvents.md` |
| 6 Inflow and photoelectrons | `src/particles/injection/bem_injection_flux.f90`, `src/particles/injection/bem_photoelectron_injection.f90`, `docs/ParticleSourcesBoundaries.md`, `docs/PhotoelectronEmission.md` |
| 7 Charge update and conductors | `src/physics/bem_surface_models_conductor.f90`, `docs/SurfaceChargeNumerics.md`, `docs/BatchDurationTheory.md` |
| 8 Outer sheath | `src/runtime/sheath/bem_surface_current_model.f90`, `src/runtime/sheath/bem_zhao_outflow_refresh.f90`, `docs/ZhaoStationaryClosure.md`, `docs/ZhaoUpstreamBand.md`, sheath-model (commit pinned in `fpm.toml`) |
| 9 Verification | `docs/ValidationGuide.md`, `docs/PhysicsReleaseVerification.md`, `docs/BatchDurationStability.md` |
| 10 Forces and detachment | `docs/ObjectForcesDetachment.md`, `docs/PythonPostprocessAPI.md` |

## Build the PDF

Required: XeLaTeX, BibTeX, `xeCJK`, standard LaTeX packages, and TeX Gyre/Noto CJK fonts. The manuscript uses
Noto Serif CJK JP, Noto Sans CJK JP, and Noto Sans Mono CJK JP. To use other fonts, edit the font settings
in `beach_models.tex`. Shell escape is not required.

On a local workstation or within a compute-node allocation, run:

```bash
cd /path/to/BEACH/docs/manuscript
make
```

The output is `build/beach_models.pdf`. Auxiliary files also stay under the ignored `build/` directory.
The build runs LaTeX, BibTeX, and two further LaTeX passes to resolve citations and cross-references.

On KUDPC login nodes, submit the typesetting payload to a compute node. Check the Sys module,
`spartition`, and `qgroup`, and replace the queue below with a verified available queue:

```bash
tssrun -p <queue> -t 0:10:0 --rsc p=1:t=1:c=1 \
  bash -lc 'cd /path/to/BEACH/docs/manuscript && make'
```

### Build verification in this environment

On October 13, 2026, the manuscript produced a 30-page A4 PDF on a System C compute node. There are no unresolved
citations or references, overfull lines, or BibTeX warnings. Missing `xeCJK` and LaTeX3 bundles were unpacked from
the matching TeX Live generation into `BEACH/build/manuscript-tools/texmf2018/` for the existing TeX Live 2018
installation. The installed TeX environment was not modified. This temporary directory is ignored by Git.
In a normal TeX environment, install the required packages and use `make` as above.

To submit to a system other than the one selected by the login node's Sys module, switch it with `SLURM_CONF`:

```bash
SLURM_CONF=/etc/slurm/sysC/slurm.conf tssrun -p gr20001c -t 0:10:0 --rsc p=1:t=1:c=1 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH/docs/manuscript && TEXINPUTS=/LARGE0/gr20001/b36291/Github/BEACH/build/manuscript-tools/texmf2018/tex//: make'
```

## Sources and updates

The manuscript was first written against the working tree of October 5, 2026 (base
`57561eb2fe795d84e06e5b6abaf29ff15c1d8070`). On October 13, 2026, Chapter 8 was revised against the implementation
(commit `7e56bb5`, sheath-model `b4a2521`). On the same day all chapters were rewritten as a description of the models
and algorithms without implementation-specific statements. This revision also brought the periodic far-field operator in
Chapter 4 (node-to-node operators at depth 2 and the zero-mode corrections $s$ and $S$) and the removed coupling path in
Chapter 5 in line with the implementation. When the implementation changes, check the matching chapter in the table above.

DOIs, authors, titles, and publication details were checked using Crossref, publisher/institutional sources, and arXiv.
Content checks used available full text or abstracts; this manuscript does not claim a complete full-text reading of
every reference. Codex and Claude assisted with this draft. Authors, affiliations, contributions, funding, and competing
interests remain to be completed by the authors.

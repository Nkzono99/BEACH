# BEACH physical models and numerical algorithms

Lang: [English](README.en.md) | [日本語](README.md)

This is a Japanese LaTeX manuscript in paper/technical-book form covering BEACH's models, algorithms,
and related research. Start with [beach_models.tex](beach_models.tex). Chapter sources are under
`sections/`, and the bibliography is [references.bib](references.bib). It includes Japanese and English
abstracts, a table of contents, equations, 20 figures, a list of figures, and references.

## Contents

1. Surface charging, airless bodies, and related research
2. Model partition, governing equations, and batches
3. P0 triangular sources, analytic integration, self terms, and one-sided fields
4. Direct, treecode, Cartesian FMM, and two-periodic Ewald/zero-mode fields
5. Same-time Boris updates, collisions, box boundaries, and periodic images
6. Boundary inflow, velocity distributions, ray-based photoelectron emission, and secondary electron emission on absorption
7. Charge budgets, return closure, floating conductors, and batch stability
8. Outer-sheath model and its connection to the top face (assumptions, orbits and populations, derived densities and fluxes, first integral of Poisson's equation and Types A, B, and C, closures, admissibility, numerical algorithms, numerical examples, closure conditions and type choice, electron drift and field floor, top-face mapping, outer-root refresh)
9. Verification, convergence, Monte Carlo errors, and reproducibility
10. Force and detachment analysis, methods and limits to report, and manuscript statements

## Explanatory figures

The manuscript includes the computation-cycle diagram and the following 19 figures at the relevant points in the text.
The figures use TikZ vectors so that text and lines remain sharp when enlarged. Drawing sources are under
`figures/`, with shared colors and line styles in [figures/styles.tex](figures/styles.tex).
Edit a drawing and run `make` to update the PDF. No additional drawing software or shell escape is required.

| Chapter | Figure content | Drawing source |
| --- | --- | --- |
| 2 | Computational domain and model partition | [model_domain.tex](figures/model_domain.tex) |
| 3 | Triangle evaluation geometry and one-sided fields | [panel_geometry.tex](figures/panel_geometry.tex) |
| 4 | FMM aggregation, translation, and sharing | [fmm_translations.tex](figures/fmm_translations.tex) |
| 4 | Periodic images and zero-mode partition | [periodic_split.tex](figures/periodic_split.tex) |
| 5 | Boris update and first geometric event | [particle_step.tex](figures/particle_step.tex) |
| 6 | Volume and boundary-flux distributions | [flux_distribution.tex](figures/flux_distribution.tex) |
| 6 | Illumination shadowing and photoelectron charge budgets | [photoemission.tex](figures/photoemission.tex) |
| 7 | Batch-width stability of linear relaxation | [batch_stability.tex](figures/batch_stability.tex) |
| 8 | Electron energy diagram and orbit classes (numerical) | [sheath_model_populations.tex](figures/sheath_model_populations.tex) |
| 8 | Charge density and first integral (numerical) | [sheath_model_first_integral.tex](figures/sheath_model_first_integral.tex) |
| 8 | Outer-sheath Types A, B, and C | [sheath_types.tex](figures/sheath_types.tex) |
| 8 | Upstream band and field deficit (numerical) | [sheath_model_upstream_band.tex](figures/sheath_model_upstream_band.tex) |
| 8 | Zero-current root search | [sheath_model_algorithm.tex](figures/sheath_model_algorithm.tex) |
| 8 | Potential and field of Types A, B, and C (numerical) | [sheath_model_profiles.tex](figures/sheath_model_profiles.tex) |
| 8 | Density of each population (numerical) | [sheath_model_densities.tex](figures/sheath_model_densities.tex) |
| 8 | Electron velocity distributions (numerical) | [sheath_model_vdf.tex](figures/sheath_model_vdf.tex) |
| 8 | Zero-current roots versus solar elevation (numerical) | [sheath_model_alpha_sweep.tex](figures/sheath_model_alpha_sweep.tex) |
| 8 | Types with an admissible root (numerical) | [sheath_model_type_map.tex](figures/sheath_model_type_map.tex) |
| 8 | Outer-sheath connection to the top face | [sheath_mapping.tex](figures/sheath_mapping.tex) |

Schematic shapes, trajectories, and potential curves are not numerical solutions. The flux-distribution and
batch-width plots are analytic examples derived from the equations in the text. The Chapter 8 figures marked
"numerical" are solutions of the outer-sheath model. Captions state this distinction;
[figures/figure_trace.yaml](figures/figure_trace.yaml) records the source sections/equations, source data, and the limitations of each figure.

### Regenerating the Chapter 8 numerical figures

The numerical figures are made in two steps, computation and drawing. The generated TikZ sources and their data
(`figures/data/sheath_model/`) are committed, so `make` alone needs no Python. Only when the conditions or the
outer-sheath model change, run the following on a compute node. The computation needs the sheath-model Python
package (NumPy and SciPy).

```bash
cd /path/to/BEACH/docs/manuscript
# Numerical solutions to CSV/JSON (about 9 minutes with 32 workers, mostly the type-map grid)
PYTHONPATH=/path/to/sheath-model python tools/sheath_model_data.py figures/data/sheath_model --workers 32
# CSV to figures/sheath_model_*.tex (NumPy only)
python tools/sheath_model_figures.py figures/data/sheath_model figures
make
```

`provenance.json` records the sheath-model commit and the common physical inputs.
`sheath_model_algorithm.tex` is drawn by hand and is not generated.

## Writing policy

The text assumes general knowledge of plasma physics and numerical methods and keeps the models and algorithms
self-contained. It does not name input keys, source paths, output files, tests, commits, or computing environments.
Origins of the methods are given by citations. The correspondence to the implementation is kept in this table.

| Chapter | Main implementation and documents |
| --- | --- |
| 2 Computation cycle | `src/runtime/simulator/bem_simulator_loop.f90`, `docs/models/Algorithms.md` |
| 3 Triangle fields | `src/physics/field_solver/panel/bem_panel_kernel.f90`, `src/physics/field_solver/panel/bem_panel_self_terms.f90`, `docs/developer/DirectSolver.md` |
| 4 Treecode and FMM | `docs/developer/Treecode.md`, `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90`, `docs/developer/FMMCore.md` |
| 4 Periodic far field | `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90`, `docs/developer/PeriodicFarCorrection.md` |
| 4 Zero mode | `src/physics/field_solver/periodic/bem_periodic_zero_mode_plan.f90`, `src/physics/field_solver/periodic/bem_periodic_zero_mode_eval.f90`, `docs/models/PeriodicElectrostatics.md` |
| 5 Particle integration and collisions | `src/physics/bem_pusher.f90`, `src/physics/bem_collision.f90`, `src/runtime/simulator/bem_particle_stepper.f90`, `docs/models/BorisPusher.md`, `docs/developer/ParticleEvents.md` |
| 6 Inflow and electron emission | `src/particles/injection/bem_injection_flux.f90`, `src/particles/injection/bem_photoelectron_injection.f90`, `src/particles/bem_secondary_emission.f90`, `src/runtime/simulator/bem_simulator_particles.f90`, `docs/guide/ParticleSourcesBoundaries.md`, `docs/models/PhotoelectronEmission.md` |
| 7 Charge update and conductors | `src/physics/bem_surface_models_conductor.f90`, `docs/developer/SurfaceChargeNumerics.md`, `docs/models/BatchDurationTheory.md` |
| 8 Outer-sheath model | sheath-model (commit pinned in `fpm.toml`): `docs/model.md`, `docs/numerics.md`, `docs/fortran-api.md`, `docs/coupling.md`, `src/internal/sheath_model_admissibility.f90`, `src/sheath_model_emission.f90`, `src/sheath_model_matching.f90`, `sheath_model/_physics.py`, `sheath_model/_profile.py`, `sheath_model/search.py` |
| 8 Top-face connection | `src/runtime/sheath/bem_surface_current_model.f90`, `src/runtime/sheath/bem_zhao_outflow_refresh.f90`, `docs/models/ZhaoStationaryClosure.md`, `docs/models/ZhaoUpstreamBand.md` |
| 9 Verification | `docs/guide/ValidationGuide.md`, `docs/developer/PhysicsReleaseVerification.md`, `docs/guide/BatchDurationStability.md` |
| 10 Forces and detachment | `docs/guide/ObjectForcesDetachment.md`, `docs/reference/PythonPostprocessAPI.md` |

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

On October 15, 2026, the manuscript produced a 50-page A4 PDF with 20 figures on a System C compute node.
There are no unresolved
citations or references, overfull lines, or BibTeX warnings. Missing `xeCJK` and LaTeX3 bundles were unpacked from
the matching TeX Live generation into `BEACH/build/manuscript-tools/texmf2018/` for the existing TeX Live 2018
installation. The installed TeX environment was not modified. This temporary directory is ignored by Git.
In a normal TeX environment, install the required packages and use `make` as above.

To submit to a system other than the one selected by the login node's Sys module, switch that module and check the queues.
The following example assumes SysG is currently loaded:

```bash
module switch SysG SysC
module list
spartition
qgroup
tssrun -p gr20001c -t 0:10:0 --rsc p=1:t=1:c=1 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH/docs/manuscript && TEXINPUTS=/LARGE0/gr20001/b36291/Github/BEACH/build/manuscript-tools/texmf2018/tex//: make'
```

## Sources and updates

The manuscript was first written against the working tree of October 5, 2026 (base
`57561eb2fe795d84e06e5b6abaf29ff15c1d8070`). On October 13, 2026, Chapter 8 was revised against the implementation
(commit `7e56bb5`, sheath-model `b4a2521`). On the same day all chapters were rewritten as a description of the models
and algorithms without implementation-specific statements. This revision also brought the periodic far-field operator in
Chapter 4 (node-to-node operators at depth 2 and the zero-mode corrections $s$ and $S$) and the removed coupling path in
Chapter 5 in line with the implementation. On October 14, 2026, Chapter 8 was added to describe the outer-sheath model itself (sheath-model `26569c1`, whose
physics is identical to `b4a2521`; only `.gitignore` differs). The assumptions, orbits and populations, derivations of
densities and fluxes, first integral, closures, admissibility, numerical algorithms, and numerical examples moved there.
On October 15, 2026, the former Chapter 8 on the stationary outer-sheath closure was merged into the last section of Chapter 8,
"Connection to the BEACH top face". The duplicated upstream-band schematic was removed in favor of the numerical Chapter 8 figure. The upstream-band density expansion
(equation `eq:hlogh`) now states the factor of 2 in front of the integral that comes from slow electrons moving in both directions.
When the implementation changes, check the matching chapter in the table above.

DOIs, authors, titles, and publication details were checked using Crossref, publisher/institutional sources, and arXiv.
Content checks used available full text or abstracts; this manuscript does not claim a complete full-text reading of
every reference. Codex and Claude assisted with this draft. Authors, affiliations, contributions, funding, and competing
interests remain to be completed by the authors.

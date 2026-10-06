# BEACH physical models and numerical algorithms

Lang: [English](README.en.md) | [日本語](README.md)

This is a Japanese LaTeX manuscript in paper/technical-book form covering BEACH's models, algorithms,
and related research. Start with [beach_models.tex](beach_models.tex). Chapter sources are under
`sections/`, and the bibliography is [references.bib](references.bib). It includes Japanese and English
abstracts, a table of contents, equations, a computation-cycle diagram, and references.

## Contents

1. Surface charging, airless bodies, and related research
2. Model partition, governing equations, and batches
3. P0 triangular sources, analytic integration, self terms, and field traces
4. Direct, treecode, Cartesian FMM, and two-periodic Ewald/zero-mode fields
5. Same-time Boris updates, collisions, box boundaries, and periodic images
6. Reservoir inflow, velocity distributions, and ray-based photoelectron emission
7. Charge budgets, return closure, floating conductors, and batch stability
8. Stationary Zhao and matching-plane fixed points/implicit zero mode
9. Verification, convergence, Monte Carlo errors, and reproducibility
10. Interpreting `examples/beach.toml`, force/detachment analysis, and research use

The appendix maps statements to implementation files and records the scope of literature verification.

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

On October 6, 2026, the revised manuscript produced a 31-page A4 PDF on a System B compute node (Job 24860895).
All 12 citations, chapter/equation references, and implementation paths were checked. There are no unresolved
references or overfull lines. Missing `xeCJK` and LaTeX3 bundles were unpacked from the matching TeX Live
generation into `BEACH/build/manuscript-tools/texmf2018/` for the existing TeX Live 2018 installation.
The installed TeX environment was not modified. This temporary directory is ignored by Git.
In a normal TeX environment, install the required packages and use `make` as above.
Job 24860896 on the same day validated all 7 updated examples and their 14 boundary-inflow species, both with
omitted keys and with explicit default values. All 3 existing boundary-inflow/schema tests selected for this change passed.

The typesetting command used was the following. Before submission, the target was checked with
`module switch SysG/2022 SysB`, `module list`, `spartition`, and `qgroup`.

```bash
tssrun -p gr20001b -t 0:10:0 --rsc p=1:t=1:c=1 \
  bash -lc 'cd /LARGE0/gr20001/b36291/Github/BEACH/docs/manuscript && TEXINPUTS=/LARGE0/gr20001/b36291/Github/BEACH/build/manuscript-tools/texmf2018/tex//: make'
```

## Sources and updates

The manuscript describes the local working tree inspected on October 5, 2026, based on commit
`57561eb2fe795d84e06e5b6abaf29ff15c1d8070`, including uncommitted periodic-field/zero-mode changes.
Each chapter lists its specification and implementation sources. No new simulation or benchmark was run.
The configuration examples and particle-source explanation were revised on October 6, 2026, to omit `source_mode`
and `npcls_per_step` for boundary inflow alone.
Update the text against the implementation before citing it in research results.

DOIs, authors, titles, and publication details were checked using Crossref and publisher/institutional
sources. Content checks used available full text or abstracts; this manuscript does not claim a complete
full-text reading of every reference. Codex assisted with this draft. Authors, affiliations, contributions,
funding, and competing interests remain to be completed by the authors.

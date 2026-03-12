# AGENTS.md

Always review in Japanese.

## Project overview
This repository provides **BEACH (BEM + Accumulated CHarge)**.

- Main simulation engine: **Fortran** (`src/`, `app/`, managed by fpm)
- Python layer: post-processing / visualization / workload utilities (`beach/`, `examples/`)

Core batch loop (current Fortran implementation):
1) compute E (and optional uniform B) from boundary-element charges
2) advance particles with Boris pusher
3) detect first collision against triangle elements
4) deposit charge to hit element (insulator accumulation)
5) commit charge deltas per batch and update stats/history

v0.x focuses on insulator accumulation. Keep extension points for conductor/resistive models.

## Repo layout
- `src/`, `app/`: Fortran runtime and physics core
- `tests/fortran`: Fortran unit/integration tests (fpm targets)
- `beach/`: Python package for reading/plotting Fortran outputs
- `tests/python`: Python tests for CLI/result utilities
- `examples/`: runnable configs and helper scripts
- `docs/`, `SPEC.md`: behavior/spec documentation

## Setup
Python utilities:
- `python -m pip install "git+https://github.com/Nkzono99/BEACH.git"` (builds Python + Fortran `beach` binary via `make`)
- `python -m pip install -U pip setuptools wheel`
- `python -m pip install -e . --no-build-isolation`

Fortran execution/testing requires a Fortran compiler (`gfortran` etc.) and `fpm`.

## How to run
Fortran main (recommended):
- `python -m pip install "git+https://github.com/Nkzono99/BEACH.git"`
- `beach examples/beach.toml`

Fortran main (developer direct run):
- `make`
- `fpm run --profile release --flag "-fopenmp" -- examples/beach.toml`

Python CLI examples:
- `beach-inspect outputs/latest`
- `beach-animate-history outputs/latest --quantity charge --save-gif outputs/latest/charge_history.gif`
- `beach-estimate-workload examples/beach.toml --threads 8`

## Tests / quality
- Python tests: `pytest -q`
- Python lint: `ruff check .`
- Fortran tests: `fpm test` (or per-target via `fpm test --target <name>`)

## Coding rules
- Do not change public APIs without updating docs + examples.
- For core simulation behavior, treat Fortran implementation and docs (`SPEC.md`) as source of truth.
- Keep algorithms correctness-first; gate performance features behind flags.
- Add/extend tests when modifying field, collision, boundary, injection, or resume logic.
- Keep Python side lightweight; avoid heavy dependencies unless justified.
- Ignore `*.i90` files, as they are automatically generated backup files.

## Simulator invariants (v0.x)
- Interaction: absorption only (default).
- Surface model: insulator accumulation only.
- Batch execution: run for configured `sim.batch_count` batches.
- Per-particle advance limit: `sim.max_step`.
- `tol_rel` is currently a monitoring/output metric and **not** an early-stop condition.

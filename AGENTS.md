# AGENTS.md

Always review in Japanese.

## Project overview
This repository implements a Python prototype of a BEM-based test-particle simulator.
Core loop:
1) compute E (and optional B) from boundary element charges
2) advance particles with Boris pusher
3) detect collisions against triangle elements
4) deposit charge to hit element (insulator mode)
5) commit every `npcls_per_step` and rebuild field cache

Insulator mode only in v0.x. Keep design extensible for conductor/resistive surface models.

## Repo layout
- `src/` : library code
- `examples/` : minimal runnable scripts
- `docs/` : specs and notes

## Setup commands
- Create venv: `python -m venv .venv`
- Activate: `source .venv/bin/activate` (or Windows equivalent)
- Install: `pip install -e .[dev]`

## Setup
python -m pip install -U pip setuptools wheel
python -m pip install -e . --no-build-isolation

## How to run
- Example: `python examples/run_absorption_insulator.py`

## Tests / quality
- Run tests: `pytest -q`
- Lint: `ruff check .`
- Format: `black .`

## Coding rules
- Do not change public APIs without updating docs + examples.
- Prefer pure NumPy; avoid heavy dependencies unless justified.
- Keep algorithms correct-first; add performance features behind flags.
- Add type hints for public functions/classes.
- When modifying collision or field logic, add/extend tests.

## Simulator invariants (v0.x)
- Interaction: absorption only (default).
- Surface model: insulator accumulation only.
- Batch update: commit per `npcls_per_step`.
- Stop condition: rel_change < tol OR processed >= max_step (configurable).

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

v1.0 focuses on insulator accumulation. Keep extension points for conductor/resistive models.

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
- Fortran build check: `make check` (uses stable dev version metadata for incremental fpm builds)
- Tiered tests:
  - L0 static/build: `make test-l0` (`git diff --check`, JSON schema parse checks, `make check`)
  - L1 development loop: `make test` or `make test-l1` (Python tests + quick Fortran targets, excluding heavy FMM)
  - L2 contract/integration: `make test-l2` (L1 + C/kernel contract target)
  - L3 release/heavy: `make test-l3`, `make test-heavy`, or `make test-full` (heavy FMM/full fpm suite)
- Per-target Fortran test: `FPM_ACTION=test ./build.sh --target <name>`
- Fortran format: `fprettify -i 2` for `*.f90` / `*.F90` (`*.i90` is excluded). Prefer `make fmt-fortran`, `make fmt-check-fortran`, and `pre-commit install`.
- KUDPC のログインノード（`camphor*`, `laurel*`, `cinnamon*`, `gardenia*`）、アプリケーションサーバ、ファイル転送サーバでは、`make test*` / `fpm test` / 長時間の Fortran テストターゲットを直接実行しない。KUDPC plugin のホスト判定または `hostname` でノード種別を確認し、Slurm 割当内でない場合は計算ノードに投入する。
- KUDPC で全体テストや長時間テストが必要な場合は、計算ノード上で実行する。短い対話確認は `tssrun`、バッチ実行は `sbatch` ジョブ内の `srun` を使う。ログインノードでは編集、軽いログ確認、`make check` 程度のビルド確認、ジョブ投入・監視までに留める。
- `make test` は開発用 L1 とし、重い FMM 系（`test_dynamics_fmm`, `test_coulomb_fmm_core_basic`, `test_coulomb_fmm_core_periodic`, `test_periodic2_flat_oracle_diag`）は `make test-l3` / `make test-heavy` / `make test-fortran-heavy` / `make test-full` で明示実行する。新しい長時間 FMM 回帰テストは、`make test` の既定経路に足す前に opt-in の L3/heavy へ分離する。
- Do not run multiple `fpm test` commands in parallel. They share the same `build/` directory and can conflict with each other.
- `fpm test` は時間がかかるため、Bash ツールの `run_in_background: true` でバックグラウンド実行し、完了通知を待つ間に他の作業（コードレビュー・編集など）を並行して進めること。ただし同時に複数の `fpm test` をバックグラウンドで走らせてはいけない。

## Coding rules
- Do not change public APIs without updating docs + examples.
- When bumping versions or preparing a release, update both `pyproject.toml` and `fpm.toml` together to keep Python and Fortran package versions in sync.
- For core simulation behavior, treat Fortran implementation and docs (`SPEC.md`) as source of truth.
- Keep algorithms correctness-first; gate performance features behind flags.
- Add/extend tests when modifying field, collision, boundary, injection, or resume logic.
- Keep Python side lightweight; avoid heavy dependencies unless justified.
- Ignore `*.i90` files, as they are automatically generated backup files.

## Simulator invariants (v1.0)
- Interaction: absorption only (default).
- Surface model: insulator accumulation only.
- Batch execution: run for configured `sim.batch_count` batches.
- Per-particle advance limit: `sim.max_step`.
- `tol_rel` is currently a monitoring/output metric and **not** an early-stop condition.

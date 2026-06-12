---
name: beach-run-diagnose
description: Diagnose failed BEACH installs, builds, runs, abnormal termination, missing outputs, restart failures, config parser errors, and suspicious simulation statistics.
---

Use this skill when a user provides a failed command, build log, runtime log, output directory listing, Slurm result, stack trace, or abnormal BEACH result.

## Response Language

- Respond in the user's language. If the request mixes languages and includes Japanese, use Japanese.
- Keep code identifiers, TOML keys, file paths, commands, log excerpts, and module names unchanged.
- If language is unclear, default to Japanese.

## Context Sources

Prefer:

- User command, log, config file, output listing, recent parameter changes, and installed version.
- Bundled references: `../../references/SPEC.md`, `../../references/README.md`, `../../references/fortran_parameter_file.md`, `../../references/config_workflow.md`, `../../references/fortran_workflow.md`, `../../references/python_postprocess_api.md`.
- Bundled docs: `../../docs/known-failure-modes.md`, `../../docs/usage-workflows.md`, `../../docs/simulator-context.md`.
- Repo root docs/source only when the full checkout is available and may be newer.

## Diagnosis Flow

1. Classify the failure phase:
   - install/package build
   - Fortran/fpm compile or link
   - config rendering or validation
   - config parsing in Fortran runtime
   - initialization
   - main-loop physics/numerics
   - output writing or Python post-processing
   - resume/continuation
2. Extract the first meaningful error. Later MPI/fpm/Python errors are often downstream.
3. Map symptoms to BEACH facts:
   - `beachx config validate` checks `case.toml` and preset composition.
   - The Fortran runtime reads rendered `beach.toml`.
   - Unknown TOML sections/keys are runtime parser errors.
   - `tol_rel` is not an early-stop condition.
   - `history_stride <= 0` means no `charge_history.csv`.
   - `output.write_mesh_potential = false` means no `mesh_potential.csv`.
   - Resume needs compatible `summary.txt`, `charges.csv`, and RNG state files.
4. Propose the smallest confirmation command.
5. If code changes are needed, scope them to the failing phase and avoid unrelated interfaces.

## Common Symptom Map

- Missing compiler or `fpm`: install/build prerequisite problem.
- `pip install beach-bem` build failure: Fortran toolchain, `make install`, or `INSTALL_PROFILE` issue.
- `Unknown key` / `Unknown section`: rendered `beach.toml` key mismatch with Fortran parser.
- `batch_duration` errors: mutually exclusive or missing resolved positive duration for `reservoir_face` / `photo_raycast`.
- Empty or tiny absorption: injection direction, injection region, mesh placement, boundary condition, or `max_step`.
- Large `survived_max_step`: particles remain alive until `sim.max_step`.
- Missing output file: output flag/cadence, early termination, wrong `output.dir`, or looking at stale directory.
- Resume failure: missing/stale files, changed mesh/config, or MPI world-size mismatch.
- Python analysis mismatch: Python direct reconstruction does not fully reproduce Fortran FMM/periodic far correction.

## Output

Use the response language and translate headings when appropriate:

```text
## 実行診断

### フェーズ
...

### 根本原因候補
...

### 確認した根拠
...

### 最小の対処
1. ...

### 追加で必要な情報
- ...
```

When a command was requested, include the important output lines in the answer because the user may not see tool output directly.

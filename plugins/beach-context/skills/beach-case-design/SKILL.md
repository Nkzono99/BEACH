---
name: beach-case-design
description: Design BEACH simulation cases, presets, parameter edits, and sweeps from physical intent while mapping choices to case.toml, beach.toml, mesh, species, solver, output, and workload controls.
---

Use this skill when a user asks what parameter controls a phenomenon, how to build a new BEACH setup, how to choose presets, or how to plan a parameter sweep.

## Response Language

- Respond in the user's language. If the request mixes languages and includes Japanese, use Japanese.
- Keep code identifiers, TOML keys, preset names, commands, and file paths unchanged.
- If language is unclear, default to Japanese.

## Source Order

1. User's scientific objective, existing config, constraints, and target output.
2. Bundled references: `../../references/config_workflow.md`, `../../references/fortran_parameter_file.md`, `../../references/SPEC.md`, `../../references/examples/periodic2_basic/case.toml`, `../../references/examples/beach.toml`, `../../references/presets/`.
3. Bundled docs: `../../docs/simulator-context.md`, `../../docs/usage-workflows.md`, `../../docs/known-failure-modes.md`.
4. Repo docs/source only when the full checkout is available and may be newer.

## Design Workflow

- Decide whether the user should edit `case.toml` + presets or a rendered `beach.toml`. Prefer `case.toml` + presets for new work.
- Identify the physical objective: reservoir plasma, photoelectron emission, object charging, periodic slab, mobility, force/torque, or workload test.
- Separate physical controls from numerical controls:
  - physical: species charge/mass, density, temperature, drift velocity, injection face, external `e0`/`b0`, mesh/object geometry, photo current.
  - numerical: `dt`, `batch_duration`, `batch_count`, `max_step`, macro-particle target, solver mode, softening, tree/FMM controls, history cadence.
- Preserve comparable diagnostics across sweep cases: output flags, history cadence, seed policy, mesh resolution, and solver settings.
- For sweeps, vary one primary physical control at a time unless the user explicitly asks for coupled scaling.
- Estimate cost risk before proposing large batch/mesh/history settings.
- Keep v1 scope clear: absorption and insulator accumulation are standard; conductor/resistive/secondary-emission models are extension points unless the repo now documents them.

## Output

For explanations:

```text
## パラメータ説明
- 対象:
- 関連セクション:
- 意味:
- 関連する制約:
- 参照:
```

For setup or sweep design:

```text
## 設計案

### 固定する条件
- ...

### 変更するパラメータ
- ...

### sweep 表
| case | parameter | value | reason |
| --- | --- | --- | --- |

### 実行前チェック
- ...
```

Use TOML snippets only for the fields that matter. Do not generate a huge full config unless the user asks for one.

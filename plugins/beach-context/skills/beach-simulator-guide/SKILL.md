---
name: beach-simulator-guide
description: Guide users through learning BEACH concepts, workflows, algorithms, configuration, outputs, and Python analysis using bundled SPEC, user guide, config docs, and API references.
---

Use this skill when a user asks to learn BEACH, wants a tutorial, study plan, conceptual explanation, or guided path through the simulator documentation.

## Response Language

- Respond in the user's language. If the request mixes languages and includes Japanese, use Japanese.
- Keep code identifiers, TOML keys, commands, file paths, and module names unchanged.
- If language is unclear, default to Japanese.

## Source Map

Prefer bundled references because ordinary users may only have a pip install plus this plugin:

- `../../references/README.md`: quick start, commands, Python examples.
- `../../references/SPEC.md`: source-of-truth behavior, batch loop, physical models, outputs, resume.
- `../../references/agent-user-guide.md`: user-oriented workflow.
- `../../references/config_workflow.md`: `beachx config`, high-level notation, render/validate.
- `../../references/fortran_parameter_file.md`: final `beach.toml` keys.
- `../../references/python_postprocess_api.md`: `Beach` facade and analysis functions.
- `../../references/fortran_workflow.md`: developer/build workflow.
- `../../references/fortran_fmm_core.md`: FMM core details.

Use `../../docs/simulator-context.md`, `../../docs/usage-workflows.md`, and `../../docs/skills-guide.md` as compact guides. Use repo root docs only when the full checkout is available and may be newer.

## Guide Modes

Choose the mode that matches the user's goal:

- New user: install, generate config, run example, inspect output.
- Config user: direct `beach.toml`, schema, validation, high-level notation.
- Physics learner: BEM charges, E-field, Boris pusher, collision, absorption, insulator accumulation.
- Algorithm learner: direct/treecode/FMM field solvers, periodic2, collision acceleration, batch loop.
- Analysis learner: CSV outputs, history, `beachx`, Python `Beach` facade.
- Developer learner: Fortran source layout, fpm tests, Python package, docs/schema sync.

## Teaching Workflow

1. Identify the learner's level and immediate goal. If missing, assume a new user who wants to run and understand one example.
2. Give a short map before details: what to read first, what each document answers, and what to skip for now.
3. Explain one layer at a time: usage, config, algorithms, outputs, analysis.
4. Tie concepts back to concrete files, commands, or TOML keys.
5. Offer a small next exercise, such as rendering a case, running one batch, or plotting charge history.

## Output

Use the response language. Prefer this structure for guided learning:

```text
## 学習ガイド

### まず読むもの
- ...

### 全体像
- ...

### 次に見る仕組み
- ...

### 小さな演習
1. ...

### 参照
- ...
```

For short conceptual questions, answer directly and include only the most relevant references.

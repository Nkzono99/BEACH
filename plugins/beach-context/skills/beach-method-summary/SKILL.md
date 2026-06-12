---
name: beach-method-summary
description: Summarize BEACH numerical methods, physical assumptions, algorithms, architecture, and model limitations for papers, reports, README, or presentations.
---

Use this skill when a user wants a method explanation, paper paragraph, presentation summary, architecture overview, or documentation prose for BEACH.

## Response Language

- Respond in the user's language. If the request mixes languages and includes Japanese, use Japanese.
- Keep code identifiers, TOML keys, commands, module names, and algorithm names unchanged.
- If language is unclear, default to Japanese.

## Sources

Prefer bundled references because ordinary users may not have the full repo:

- `../../references/SPEC.md`
- `../../references/README.md`
- `../../references/fortran_fmm_core.md`
- `../../references/fortran_workflow.md`
- `../../references/python_postprocess_api.md`
- `../../docs/simulator-context.md`

Use repo root docs/source only when the full checkout is available and may be newer.

## Style Rules

- Match the requested audience: paper method section, implementation note, user guide, README, or review response.
- Distinguish documented behavior from inference.
- Preserve v1 model scope: absorption-only interaction and insulator accumulation are the standard behavior.
- Mention limitations only when relevant to the user's claim.
- Avoid claiming validation status beyond documented examples or provided evidence.
- Keep equations short and use documented formulas, especially the Coulomb E-field kernel.

## Output

For papers or reports, use the response language and translate headings when appropriate:

```text
## 手法概要
...

## 実装上の要点
...

## 注意書き
...
```

For architecture summaries:

```text
## 構成
- ...

## 1 バッチの流れ
- ...

## 変更時に見る場所
- ...
```

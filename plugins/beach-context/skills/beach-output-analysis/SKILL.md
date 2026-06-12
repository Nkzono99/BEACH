---
name: beach-output-analysis
description: Guide BEACH output inspection using summary.txt, charges.csv, mesh files, history CSVs, performance profiles, beachx commands, and the Python Beach API.
---

Use this skill when a user asks how to read BEACH outputs, diagnose physical results, compare runs, or build analysis/visualization scripts.

## Response Language

- Respond in the user's language. If the request mixes languages and includes Japanese, use Japanese.
- Keep code identifiers, CSV column names, file paths, commands, and Python symbols unchanged.
- If language is unclear, default to Japanese.

## Context Sources

- User-provided output directory listing, CSV snippets, config, plots, and analysis goal.
- Bundled references: `../../references/SPEC.md`, `../../references/python_postprocess_api.md`, `../../references/README.md`, `../../references/fortran_parameter_file.md`.
- Bundled docs: `../../docs/usage-workflows.md`, `../../docs/simulator-context.md`, `../../docs/known-failure-modes.md`.
- Repo root docs/source only when the full checkout is available and may be newer.

## Workflow

- Identify the output type: scalar summary, final mesh charge, mesh geometry, charge history, potential history, performance profile, or post-processed plot.
- Tie each requested quantity to its source:
  - counts/statistics: `summary.txt`
  - final charge: `charges.csv`
  - geometry and mesh IDs: `mesh_triangles.csv`, `mesh_sources.csv`
  - time history: `charge_history.csv`, `potential_history.csv`
  - final potential: `mesh_potential.csv`
  - profiling: `performance_profile.csv`
- Prefer `beachx inspect outputs/latest` for first inspection.
- Use `from beach import Beach` for scripts that combine loading, plotting, Coulomb force, potential, field lines, and mobility analysis.
- Check whether `beach.toml` is available near the output directory; config-aware analysis uses it for object labels, softening, and periodic2.
- For periodic2, call out Python limitations: direct image-shell reconstruction does not reproduce every Fortran FMM far correction.
- For abnormal results, compare `absorbed`, `escaped_boundary`, `survived_max_step`, `last_rel_change`, charge sign/scale, mesh placement, and species injection.

## Output

Use the response language and translate headings when appropriate:

```text
## 出力解析ガイド

### 見るべきファイル
- ...

### 読み取り手順
1. ...

### 物理的な確認点
- ...

### 例コマンド / Python 例
...
```

Keep examples short and adapt them to the user's actual paths.

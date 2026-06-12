---
name: beach-config-review
description: Review BEACH beach.toml, case.toml, or preset files for schema, parameter, physics, injection, boundary, mesh, output, and resume consistency.
---

Use this skill when a user asks to review, validate, sanity-check, or explain a BEACH config file.

## Response Language

- Respond in the user's language. If the request mixes languages and includes Japanese, use Japanese.
- Keep code identifiers, TOML keys, file paths, commands, and schema names unchanged.
- If language is unclear, default to Japanese.

## Context Sources

Prefer sources in this order:

1. User-provided `beach.toml`, `case.toml`, preset, or generated snippet.
2. Bundled references: `../../references/fortran_parameter_file.md`, `../../references/config_workflow.md`, `../../references/SPEC.md`, `../../references/schemas/beach.schema.json`, `../../references/schemas/beach.case.schema.json`, `../../references/schemas/beach.preset.schema.json`, `../../references/examples/beach.toml`, `../../references/examples/periodic2_basic/case.toml`.
3. Bundled docs: `../../docs/simulator-context.md`, `../../docs/known-failure-modes.md`, `../../docs/usage-workflows.md`, `../../docs/skills-guide.md`.
4. Repo primary docs/source only when the full checkout is available and may be newer.

Do not invent defaults. If a value is not present and not documented in accessible sources, mark it unknown.
For `case.toml` or presets, prefer calling or suggesting `beachx config validate` when the CLI is available.
For rendered `beach.toml`, prefer schema/TOML inspection and `beachx estimate-workload <config> --threads N` when workload is relevant.

## Review Checklist

- Determine input type: `case.toml`, preset, rendered `beach.toml`, or legacy/direct config.
- Check schema directive style. Use `#:schema ...`, not `"$schema"` before the first section.
- For `case.toml`, check `schema_version`, `use_presets`, `override`, preset resolution, and `id`-based merge expectations.
- Check required sections: `[sim]`, `[[particles.species]]`, `[mesh]`, and `[output]` for rendered configs.
- Check batch controls: `dt`, `batch_count`, `batch_duration` vs `batch_duration_step`, `max_step`, and `tol_rel` expectations.
- Check source modes:
  - `volume_seed`: total `npcls_per_step > 0`, no `target_macro_particles_per_batch`.
  - `reservoir_face`: density, temperature, `inject_face`, region, drift direction, and `batch_duration > 0`.
  - `photo_raycast`: rays/current settings, source face, hit expectations, and opposite-charge deposition setting.
- Check mutually exclusive forms: `temperature_k` vs `temperature_ev`, `e0` vs `e0_abs` angle form.
- Check box and boundary consistency, especially `field_bc_mode = "periodic2"` with exactly two periodic axes.
- Check mesh placement relative to box, injection faces, and periodic primitive cell assumptions.
- Check solver controls: `field_solver`, `softening`, tree/FMM parameters, periodic far correction, and field normalization.
- Check output size risk: `history_stride`, `write_potential_history`, `write_mesh_potential`, mesh size, batch count.
- Check resume compatibility when `output.resume = true`.

## Output

Return a concise report in the response language. Translate headings when appropriate:

```text
## 設定レビュー

### 対象
- ファイル:
- 入力形式:
- 実行意図:

### 修正必須
- ...

### 確認推奨
- ...

### 問題なし / 情報
- ...

### 次の確認コマンド
- ...
```

Suggest concrete commands only when they fit the environment, for example `beachx config validate`, `beachx config render`, `beachx estimate-workload beach.toml --threads 8`, or `beach examples/beach.toml`.

---
name: beach-config-review
description: Review BEACH beach.toml files for schema, parameter, physics, injection, boundary, mesh, output, and resume consistency.
---

Use this skill when a user asks to review, validate, sanity-check, or explain a BEACH config file.

## Response Language

- Respond in the user's language. If the request mixes languages and includes Japanese, use Japanese.
- Keep code identifiers, TOML keys, file paths, commands, and schema names unchanged.
- If language is unclear, default to Japanese.

## Context Sources

Prefer sources in this order:

1. User-provided `beach.toml` or generated snippet.
2. Bundled references: `../../references/fortran_parameter_file.md`, `../../references/config_workflow.md`, `../../references/SPEC.md`, `../../references/schemas/beach.schema.json`, `../../references/examples/beach.toml`, `../../references/examples/periodic2_basic.toml`.
3. Bundled docs: `../../docs/simulator-context.md`, `../../docs/known-failure-modes.md`, `../../docs/usage-workflows.md`, `../../docs/skills-guide.md`.
4. Repo primary docs/source only when the full checkout is available and may be newer.

Do not invent defaults. If a value is not present and not documented in accessible sources, mark it unknown.
Prefer calling or suggesting `beachx lint beach.toml` when the CLI is available.
For workload questions, use schema/TOML inspection and `beachx estimate-workload <config> --threads N` when relevant.

## Review Checklist

- Determine whether the input is a direct `beach.toml` or a snippet.
- Check schema directive style. Use `#:schema ...`, not `"$schema"` before the first section.
- Reject legacy composition keys such as `schema_version`, `use_presets`, `override`, and `base_case`.
- Check the seven public groups and required `[[particles.species]]` entries; apply each table’s feature conditions.
- Check batch controls: `particles.tracking.dt_s`, `batch_count`, `run.batch.duration_s` vs `run.batch.duration_steps`, `particles.tracking.max_steps_per_particle`, and `output.diagnostics.charge_rel_change_threshold` expectations.
- Check source modes:
  - `volume_seed`: total `sampling.volume_macro_particles_per_batch > 0`, no `sampling.target_macro_particles_per_batch`.
  - `inflow`: density, temperature, nonperiodic open box faces, drift direction, and `run.batch.duration_s > 0`.
  - `photo_raycast`: rays/current settings, source face, hit expectations, and opposite-charge deposition setting.
- Check mutually exclusive forms: `distribution.temperature_k` vs `distribution.temperature_ev`, `fields.external.electric_v_m` vs `fields.external.electric_magnitude_v_m` angle form.
- Check box and boundary consistency, especially `fields.boundary = "periodic2"` with exactly two periodic axes.
- Check mesh placement relative to box, injection faces, and periodic primitive cell assumptions.
- Check solver controls: `fields.solver.method`, tree/FMM parameters, periodic far correction, and field normalization.
- Treat the P0 triangle panel as the fixed element source.
- Check output size risk: `output.history.stride_batches`, `output.history.potential`, `output.final.mesh_potential`, mesh size, batch count.
- Check resume compatibility when `[run.restart]` is present.
- Use `fields.boundary` for field closure and `sheath.closure` for the outer sheath.

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

Suggest concrete commands only when they fit the environment, for example `beachx lint beach.toml`, `Fortran parser`, `beachx estimate-workload beach.toml --threads 8`, or `beach examples/beach.toml`.

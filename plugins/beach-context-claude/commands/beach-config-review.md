---
description: "BEACH の beach.toml をレビューする"
argument-hint: "<beach.toml path or snippet>"
---

`beach-config-review` agent を使い、指定された `beach.toml` または snippet をレビューしてください。

確認する観点:
- schema と必須 section
- `run / domain / mesh / particles / fields / sheath / output`
- `mesh.kind` ごとの必須/任意 parameter
- `fields.solver.method` ごとの有効 parameter
- output size と workload risk

参照は `docs/reference/Parameters.md`, `schemas/beach.schema.json`, `SPEC.md`, `docs/models/Algorithms.md` を優先してください。
KUDPC login node では simulation や重い test は直接実行しないでください。

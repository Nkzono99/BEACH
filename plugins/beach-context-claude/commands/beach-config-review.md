---
description: "BEACH の beach.toml をレビューする"
argument-hint: "<beach.toml path or snippet>"
---

`beach-config-review` agent を使い、指定された `beach.toml` または snippet をレビューしてください。

確認する観点:
- schema と必須 section
- `sim`, `particles`, `mesh`, `field_solver`, `boundary`, `output`, `resume`
- `mesh.kind` ごとの必須/任意 parameter
- `field_solver` ごとの有効 parameter
- output size と workload risk

参照は`docs/Parameters.md`, `schemas/beach.schema.json`, `SPEC.md`, `docs/Configuration.md`,
`docs/FieldSolvers.md`, `docs/PeriodicElectrostatics.md`を優先してください。
`config init`の現行既定はx/y periodic、`field_periodic_image_layers=1`の$3\times3$有限画像和です。
KUDPC login node では simulation や重い test は直接実行しないでください。

---
description: "BEACHのcaseとparameter sweepを設計する"
argument-hint: "<physical goal/config/constraints>"
---

`beach-case-design` agentを使い、物理目的、既存config、計算資源、必要なdiagnosticから
`beach.toml`の変更またはparameter sweepを設計してください。

`docs/ConfigurationRecipes.md`, `docs/Parameters.md`, `docs/ValidationGuide.md`, `SPEC.md`,
`schemas/beach.schema.json`を優先してください。`config init`を基準にする場合、現行既定は
x/y periodic、`field_solver="fmm"`、`field_periodic_image_layers=1`、far correctionなしの
$3\times3$有限画像和です。

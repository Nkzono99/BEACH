---
description: "BEACH outputs の読み方と解析手順を組む"
argument-hint: "<output dir or analysis goal>"
---

`beach-output-analysis` agent を使い、指定 output directory または解析目的に対して見るべき file、
読み取り順、Python/CLI 例を整理してください。

出力契約は`docs/OutputGuide.md`と`SPEC.md`、Python APIは`docs/PythonPostprocessAPI.md`を優先します。

対応 file:
- `summary.txt`
- `charges.csv`
- `mesh_triangles.csv`
- `mesh_sources.csv`
- `charge_history.csv`
- `potential_history.csv`
- `mesh_potential.csv`
- `performance_profile.csv`

大きな CSV 読み込みや可視化 payload は KUDPC login node で直接実行しないでください。

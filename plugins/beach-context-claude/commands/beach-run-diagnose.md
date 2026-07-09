---
description: "BEACH の実行失敗や異常出力を診断する"
argument-hint: "<command/log/output dir/config>"
---

`beach-run-diagnose` agent を使い、失敗 command、log、output directory、config から原因を切り分けてください。

必ず行うこと:
- 失敗フェーズを分類する。
- 最初の意味ある error を抽出する。
- `SPEC.md`, `docs/Parameters.md`, `docs/Algorithms.md` の仕様と照合する。
- 最小の確認手順を提案する。

KUDPC login node では build/test/simulation/大規模解析を直接実行しないでください。

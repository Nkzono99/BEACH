title: 計算結果の妥当性確認

Lang: [日本語](ValidationGuide.md) | [English](ValidationGuide.en.md)

# 計算結果の妥当性確認

正常終了は、物理的・数値的に妥当であることを意味しません。確認を2段階に分けます。

## レベル1: 実行が完了した

- processの終了codeが0
- `summary.txt`、`charges.csv`、必要な履歴が存在する
- `batches == sim.batch_count`
- `beachx inspect`が読める
- restart時はmodel/mesh/species fingerprintが一致する

## レベル2: 結果を定量利用できる

- `absorbed`、`escaped_boundary`、`survived_max_step`の内訳を物理条件から説明できる
- `charge_ledger_residual_C`が丸め誤差範囲で、`discarded_unresolved`を別に確認した
- `tol_rel`を収束停止条件と誤解せず、履歴が十分な長さを持つ
- `sim.dt`を半分にして主要量が変わらない
- mesh解像度、FMM order/tolerance、outer gridを上げて主要量が収束する
- `batch_duration`を0.5倍/2倍にして結論が安定する
- stochastic caseはseedまたはensemble依存を確認する

## model固有の確認

| model | 必須診断 |
| --- | --- |
| periodic2 cached | cache fingerprint、cold/warm一致、zero-mode/Gauss residual |
| unified outer | accessibility refinement、linearity、outer energy/frozen-field error |
| kinetic outer | solver status、Poisson residual、Bohm/branch applicability |
| photoelectron | emission/return ledger、ambient charge ratio、histogram範囲 |

release用の小規模基準は[Physics release verification](PhysicsReleaseVerification.html)にあります。
本番caseでは同じ収束軸を自分の観測量に対して評価してください。

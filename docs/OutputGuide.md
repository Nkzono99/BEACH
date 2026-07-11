title: 出力の読み方

Lang: [日本語](OutputGuide.md) | [English](OutputGuide.en.md)

# 出力の読み方

このページは、初回実行後に `outputs/latest/` で何を見ればよいかをまとめます。
詳細な後処理コマンドは [後処理チュートリアル](PostprocessTutorial.html)、全パラメータは [入力パラメータリファレンス](Parameters.html) を参照してください。

## まず確認すること

1. `outputs/latest/summary.txt` がある。
2. `summary.txt` の `batches` が `sim.batch_count` に到達している。
3. `charges.csv` があり、要素ごとの最終電荷が出ている。
4. `beachx inspect outputs/latest` がエラーなく概要を表示する。

```bash
beachx inspect outputs/latest
```

`output.dir` を変更している場合は、`outputs/latest` を設定した出力ディレクトリに置き換えてください。

## 主要ファイル

| ファイル | いつ出るか | まず見る内容 |
| --- | --- | --- |
| `summary.txt` | 常時 | batch 数、吸収数、脱出数、最後の相対変化、MPI rank 数 |
| `charges.csv` | 常時 | 要素ごとの最終電荷 |
| `mesh_triangles.csv` | 常時 | 三角形頂点、要素 ID、`mesh_id` |
| `mesh_sources.csv` | template mesh 利用時 | `mesh_id` と template kind / surface model / 要素数の対応 |
| `mesh_potential.csv` | `output.write_mesh_potential=true` | 最終時刻の要素重心電位 |
| `charge_history.csv` | `output.history_stride > 0` | batch ごとの要素電荷履歴 |
| `potential_history.csv` | `write_potential_history=true` かつ `history_stride > 0` | batch ごとの要素重心電位履歴 |
| `performance_profile.csv` | `BEACH_PROFILE=1` | phase 別実行時間 |
| `rng_state*.txt` | 常時 | 再開用乱数状態 |
| `macro_residuals*.csv` | reservoir 系注入時 | 端数マクロ粒子数の再開用状態 |
| `charge_ledger.csv` | 常時 | species 別の signed charge flux と粒子数 |

## 成功と注意の読み分け

`summary.txt` で最初に確認する量は次のとおりです。

`field_source_model` と `field_kernel_id` は出力を生成した要素核を示します。`triangle_p0` 出力に対し、現行 Python potential/field/force/field-line estimator と point-only C kernel は誤った点電荷再構成を避けるため停止します。panel 対応診断は simulator が書いた値を利用してください。

split periodic2では`summary.txt`にinterface potential/normal field、`eta_phi_kneq0`、`eta_field_kneq0`、`eta_gap`、`eta_local_charge`、Gauss residual、outer積分電荷、最終outer更新batchを保存します。これらは物理適用性とrestart状態の一部であり、欠損したsplit checkpointは再開できません。

| 項目 | 見方 |
| --- | --- |
| `batches` | 通常実行では `sim.batch_count` と一致していれば完了 |
| `absorbed` | 表面に吸収された粒子数。帯電が進んでいるかを見る主指標 |
| `escaped` | open boundary から出た粒子数。注入・境界条件の確認に使う |
| `survived_max_step` | `sim.max_step` まで生存した粒子数。多い場合は `dt`、箱、注入条件を見直す |
| `last_rel_change` | 最終 batch の電荷変化監視値。現行実装では早期停止条件ではない |
| `charge_ledger_residual_C` | transactional 電荷保存残差。0 でも unresolved discard は別途確認する |
| `charge_ledger_discarded_unresolved_abs_C` | species 間で相殺しない max-step discard 電荷の絶対値和 |

`conductor` や `dielectric` を含む mesh では、`summary.txt` に注意書きが出る場合があります。
`dielectric` は現行実装ではメタデータであり、誘電体境界条件を解くモデルではありません。

## 履歴を使う

時間発展を見るには、`output.history_stride` を正の値にします。

```toml
[output]
history_stride = 1
write_potential_history = true
```

電位履歴は追加の場評価を行うため、要素数や履歴頻度が大きい場合は計算コストが増えます。
最初は `history_stride = 1` で小ケースを確認し、大きい計算では間引いてください。

## 次に実行するコマンド

```bash
beachx inspect outputs/latest \
  --save-mesh outputs/latest/charges_mesh.png \
  --save-potential-mesh outputs/latest/potential_mesh.png

beachx animate outputs/latest \
  --quantity charge \
  --save-gif outputs/latest/charge_history.gif
```

Python から読む場合:

```python
from beach import Beach

b = Beach("outputs/latest")
print(b.result.absorbed, b.result.escaped)
b.plot_mesh()
b.plot_potential()
```

## 再開実行の出力

`output.resume=true` の場合、`summary.txt`、`charges.csv`、`rng_state*.txt`、`macro_residuals*.csv`、schema v2 の `charge_ledger.csv` が checkpoint として使われます。schema v2 は model / mesh / species fingerprint の不一致を拒否します。
`output.restart_from` を指定すると、checkpoint は `restart_from` から読み、新しい出力は `output.dir` に書きます。

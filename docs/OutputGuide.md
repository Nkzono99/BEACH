title: 出力ファイルを調べる

# 出力ファイルを調べる

Lang: [日本語](OutputGuide.md) | [English](OutputGuide.en.md)

BEACH の出力は、最終状態、履歴、再開状態に分かれます。ファイル生成条件の機械可読な正本は
`schemas/beach.output-manifest.json` です。

## 最初に見るファイル

| 目的 | ファイル |
| --- | --- |
| 実行が完了したか、何 batch 処理したか | `summary.txt` |
| 最終的な要素電荷 | `charges.csv` |
| 三角形 geometry と mesh ID | `mesh_triangles.csv` |
| mesh ID と入力 mesh の対応 | `mesh_sources.csv` |
| species 別の電荷収支 | `charge_ledger.csv` |
| 性能内訳 | `performance_profile.csv` |

`summary.txt` には統計、設定の解決結果、build 情報、checkpoint schema、model / mesh / species fingerprint が
記録されます。境界 reservoir と通常 open 面の解決結果は `reservoir_inflow_map` と
`particle_ordinary_open_model` で確認できます。

## 履歴

| ファイル | 条件 | 主な列 |
| --- | --- | --- |
| `charge_history.csv` | `history_stride > 0` | batch、要素、電荷 |
| `potential_history.csv` | `write_potential_history=true` | batch、要素、電位 |
| `top_reference_history.csv` | 上記かつ `[domain]` の box あり | batch、時間、z-high 面の電位統計 |

`top_reference_history.csv` の基準は box の z-high 面平均です。無限遠電位やプラズマ電位ではありません。
要素相対電位は、同じ batch の `potential_history.csv` と結合し、
`potential_V - potential_mean_V` として求めます。

## mesh 電位

`write_mesh_potential=true` では `mesh_potential.csv` を出力します。

- 各要素重心の電位 [V]
- 自己項は解析的 P0 triangle panel kernel
- finite periodic2 は指定 image shell を加算
- `cached_kneq0` は cached 非零 mode と設定された物理的 zero mode を合成

## charge ledger

`charge_ledger.csv` は species ごとに次を記録します。

- 注入、表面放出、表面吸収、escape、未解決 discard の signed charge
- 各 terminal outcome の count
- closed PE の `neutral_return_correction_C`
- `neutral_return_weight_scale`
- `neutral_return_unresolved_fraction`

`summary.txt` の `charge_ledger_residual_C` は surface / local-flight / unresolved stock の変化と
外部 flux、neutral-return 補正から作る保存残差です。

closed PE では raw の吸収・未解決量を上書きしません。補正量と係数を別に記録するため、
表面総電荷が閉じていても未解決率を独立に確認できます。

## 適応 batch の診断

`periodic2.max_nonzero_mode_potential_step > 0` の場合、`summary.txt` に次を記録します。

- `simulated_time_s`
- `periodic2_max_nonzero_mode_potential_step_V`
- `adaptive_nonzero_mode_rejected_trials`
- `adaptive_nonzero_mode_last_batch_duration_s`
- `adaptive_nonzero_mode_last_potential_step_V`
- `adaptive_nonzero_mode_omp_threads`

棄却 trial は履歴や ledger へ出力されません。

## 再開に使うファイル

`output.checkpoint_stride > 0` では、accepted batch の commit 後に次の構造を更新します。

```text
outputs/latest/
├── checkpoint_latest.txt
└── checkpoints/
    ├── slot0/
    └── slot1/
```

各 slot には下表の再開用ファイル一式が入ります。非 active slot を書き終えてから
`checkpoint_latest.txt` を原子的に切り替えるため、同時に保持するのは最大 2 世代です。

| ファイル | 役割 |
| --- | --- |
| `summary.txt` | 統計、schema、fingerprint、ledger stock |
| `charges.csv` | 要素電荷 |
| `rng_state.txt` | serial の RNG |
| `rng_state_rankNNNNN.txt` | MPI rank ごとの RNG |
| `macro_residuals.csv` | global な macro 粒子数残差。species×faceを区別し、存在時に復元 |
| `charge_ledger.csv` | summary に ledger metadata がある場合に復元 |

`output.restart_from` は checkpoint の読み込み元だけを変更します。新しい出力は `output.dir` に書きます。
読み込み元に `checkpoint_latest.txt` があれば、直下の最終出力と定期 slot のうち、必須ファイルが揃い
`batches` が最大の checkpoint を自動選択します。
必須ファイルの欠落、fingerprint、mesh 要素数、species 数、MPI world size の不一致では
新規実行へ fallback せず停止します。

checkpoint schema v6の`macro_residuals.csv`は`species_idx,face,residual`です。`face=0`は従来source、
`1..6`はboundary faceを表します。旧`species_idx,residual`の2列形式も読み込めます。

## Python から読む

```python
from beach import FortranResults

result = FortranResults("outputs/latest")
print(result.summary)
print(result.charges)
```

詳細は [Python 後処理 API](PythonPostprocessAPI.md) にまとめています。
結果の物理・数値的な受理判断は
[計算結果の妥当性確認](ValidationGuide.html)を参照してください。

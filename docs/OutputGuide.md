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
`field_reconstruction_*`は、実行時に解決された`E0`、box、境界、periodic nonzero/zero-mode model、
tree設定、実際のfield solver、およびFMM展開次数のreceiptです。schema v2では
`field_reconstruction_resolved_field_solver`と`field_reconstruction_fmm_expansion_order`を必須とします。
新しい出力をPythonで電場再構成するときは近傍の`beach.toml`ではなく、このreceiptを使います。

`[surface_current_model]`のreceiptは`surface_current_model`に記録されます。`zhao_stationary`では、選択branch、
参照面積、$\phi_0$、$\phi_m$、解いたambient electron密度、およびelectron / ion / PE emission / PE escape /
PE return / netのsigned電流密度を`surface_current_model_*`で出力します。特に
`surface_current_model_photoelectron_active`はPE channelの有無を示し、falseではPE関連receiptはゼロです。
`surface_current_model_pe_return_current_density_A_m2`は負、emissionとescapeは正です。
`surface_current_model_pe_escape_particle_current_A`は外向きPE粒子が運ぶsigned電流なので負です。
二つのbudget residualは、PEのemission-return-escape連続条件と表面定常電流を独立に検証します。
`surface_current_model_kinetic_contract`とinflow access / outflow barrierの電位・face receiptは、Zhao電流に対応する
速度空間境界写像を記録します。face番号は`1..6 = x_low, x_high, y_low, y_high, z_low, z_high`です。

`matching_plane_quasistatic`では、静的な`surface_current_model_*`電流targetではなく、accepted batchの
`matching_plane_*`を読みます。`matching_plane_displacement_C_m2`と`matching_plane_phi_V`が電磁気的な整合値、
electron / ion の inward flux、access potential、PE barrier potential が応答表の出力、4つのoutward flux / energyが
固定点feedback、`matching_plane_iterations`と`matching_plane_residual`が数値的な収束receiptです。PEの
`matching_plane_photoelectron_return_flux_m2_s`と
`matching_plane_photoelectron_escape_flux_m2_s`は、同じbatchのoutward fluxとの収支確認に使います。
`surface_current_model_response_content_fingerprint`は読込済み応答表のcanonical内容を識別します。run provenanceは、
`surface_current_model_response_table_path`、`surface_current_model_matching_plane_z_m`、
`surface_current_model_electron_species`、`surface_current_model_ion_species`、
`surface_current_model_photoelectron_species`と合わせて確認します。
反復条件は`surface_current_model_coupling_rtol`、`surface_current_model_coupling_max_iterations`、
`surface_current_model_coupling_relaxation`、状態の生成元は
`surface_current_model_dynamic_state_source=accepted_batch_fixed_point`に記録されます。
`matching_plane_state_valid`がfalseなら、同じsummary内の`matching_plane_*`値をaccepted stateとして使ってはいけません。
各行の$D_H$と$\Phi_H$は、そのbatchの粒子追跡に使ったcommit前の表面電荷stateに対応します。
`simulated_time_s`はtrialを受理して進めた後の時刻なので、次batch開始時のpost-commit場と取り違えないでください。

## 履歴

| ファイル | 条件 | 主な列 |
| --- | --- | --- |
| `charge_history.csv` | `history_stride > 0` | batch、要素、電荷 |
| `potential_history.csv` | `write_potential_history=true` かつ `history_stride > 0` | batch、要素、電位 |
| `top_reference_history.csv` | 上記かつ `[domain]` の box あり | batch、時間、z-high 面の電位統計 |
| `matching_plane_history.csv` | `matching_plane_quasistatic` かつ `history_stride > 0` | batch、時間、$D_H$、$\Phi_H$、inward応答、outward / return / escape flux、反復receipt |

Python では `load_fortran_result(...)` の `matching_plane_state` と `matching_plane_history` から、列番号を
手作業で管理せず typed receipt として参照できます。

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
- `fixed_absorbed_target_charge_C` と `fixed_absorbed_weight_scale`
- `fixed_emission_target_charge_C` と `fixed_emission_weight_scale`
- 外部 closure が表面電荷へ加えた `fixed_current_correction_C`
- 互換列 `fixed_absorbed_applied_charge_C` / `fixed_emission_applied_charge_C`。各target列と同値
- 外部escapeの `fixed_escape_target_charge_C` と、同値の互換列 `fixed_escape_applied_charge_C`
- raw escapeとの差 `fixed_escape_correction_C`

`escaped_to_infinity_C`は軌道追跡で得たraw値です。Zhao固定電流ではこれを上書きせず、外部escape targetと
並べて保存します。escape補正は表面要素へdepositされないため、表面電荷更新と外部境界収支を分けて解析できます。
`*_applied_charge_C`は既存reader向けに残した出力aliasで、内部ledgerは`*_target_charge_C`だけを保持します。

`summary.txt` の `charge_ledger_residual_C` は surface / local-flight / unresolved stock の変化と
外部 flux、neutral-return 補正、fixed-current 補正から作る保存残差です。

`fixed_current`の統計診断では、対象channelの`absorbed_count`または`emitted_count`、対応するraw charge、
`fixed_*_weight_scale`を組にして読みます。countが1ならtarget全量が1標本の要素へ局所化されます。
保存残差が小さくてもこの標本分散は検出できないため、粒子数・ray数・batch幅・乱数seedを変えた
要素別電荷分布の収束を別に確認します。

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
├── checkpoint_complete.txt
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
| `macro_residuals.csv` | global な macro 粒子数残差。species×faceを区別し、schema v8 以降は manifest 宣言時に必須 |
| `charge_ledger.csv` | summary に ledger metadata がある場合に復元 |
| `checkpoint_complete.txt` | schema v8 以降の checkpoint 一式を最後に公開する完了 manifest |

`output.restart_from` は checkpoint の読み込み元だけを変更します。新しい出力は `output.dir` に書きます。
直下の最終出力と両定期 slot のうち、必須ファイルが揃う load 可能な checkpoint を比較し、
`batches` が最大のものを自動選択します。`checkpoint_latest.txt` が欠落、破損、または古い場合も、
完了 manifest を持つ slot は回収対象です。
schema v8 以降では、再開状態を書き始める前に `checkpoint_complete.txt` を `in_progress` にし、
summary、charges、全 rank の RNG、manifest で宣言した residual と ledger を閉じてから `complete` を原子的に公開します。
この schema では `checkpoint_complete.txt` 自体が再開の必須ファイルです。
直下の最終出力が中断されて新旧世代のファイルが残った場合は、その出力を選ばず完全な定期 slot へ戻ります。
必須ファイルの欠落、fingerprint、mesh 要素数、species 数、MPI world size の不一致では
新規実行へ fallback せず停止します。

checkpoint schema v6の`macro_residuals.csv`は`species_idx,face,residual`です。`face=0`は従来source、
`1..6`はboundary faceを表します。旧`species_idx,residual`の2列形式も読み込めます。
schema v9はmatching-planeのaccepted feedbackと反復receiptを`summary.txt`へ保存します。response tableの
canonical内容はmodel fingerprintに含まれるため、pathが同じでも内容を変えたcheckpointは再開できません。

## Python から読む

```python
from beach import load_fortran_result

result = load_fortran_result("outputs/latest")
print(result.charges)
print(result.matching_plane_state)
```

詳細は [Python 後処理 API](PythonPostprocessAPI.md) にまとめています。
結果の物理・数値的な受理判断は
[計算結果の妥当性確認](ValidationGuide.html)を参照してください。

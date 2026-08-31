title: 出力形式リファレンス

# 出力形式リファレンス

Lang: [日本語](OutputReference.md) | [English](OutputReference.en.md)

このページは、BEACH の出力ファイル、`summary.txt` receipt、CSV 列、checkpoint 契約を検索するための
リファレンスです。公式ケースを初めて確認する場合は、先に
[出力ファイルを調べる](OutputGuide.html)を実行してください。

ファイル生成条件の機械可読な正本は `schemas/beach.output-manifest.json` です。ここでは、その条件と
出力を読むときの意味を人間向けにまとめます。

## ファイル生成条件

| ファイル | 生成条件 | 再開での役割 |
| --- | --- | --- |
| `summary.txt` | `output.write_files=true` | 必須 |
| `charges.csv` | `output.write_files=true` | 必須 |
| `mesh_triangles.csv` | `output.write_files=true` | なし |
| `mesh_sources.csv` | OBJ または template mesh で `output.write_files=true` | なし |
| `mesh_potential.csv` | 上記かつ `output.write_mesh_potential=true` | なし |
| `charge_history.csv` | 上記かつ `output.history_stride>0` | なし |
| `potential_history.csv` | 上記かつ `output.write_potential_history=true` | なし |
| `top_reference_history.csv` | potential history の条件に加えて `[domain]` の box あり | なし |
| `matching_plane_history.csv` | `matching_plane_quasistatic` かつ `output.history_stride>0` | なし |
| `charge_ledger.csv` | charge ledger state あり | summary に ledger metadata があれば必須 |
| `rng_state.txt` / `rng_state_rankNNNNN.txt` | serial / MPI の再開状態 | 必須 |
| `macro_residuals.csv` | macro 粒子数の残差 state あり | schema v8 以降は manifest 宣言時に必須 |
| `checkpoint_complete.txt` | final output と periodic checkpoint slot | schema v8 以降の完了 manifest |
| `checkpoint_latest.txt` | `checkpoint_stride>0` で該当 accepted batch の commit 後 | periodic slot の案内 |
| `performance_profile.csv` | `BEACH_PROFILE=1` かつ `output.write_files=true` | なし |

`summary.txt` と `charges.csv` は Python reader の必須入力です。他のファイルは設定または runtime state に
応じて省略されます。

## 構成固有の値を探す

### `summary.txt` の共通 receipt

`summary.txt` には統計、設定の解決結果、build 情報、checkpoint schema、model / mesh / species fingerprint が
記録されます。境界 reservoir と通常 open 面の解決結果は `reservoir_inflow_map` と
`particle_ordinary_open_model` です。

### field と境界の receipt

`field_reconstruction_*` は、実行時に解決された `E0`、box、境界、periodic nonzero / zero-mode model、
tree 設定、実際の field solver、FMM 展開次数を記録します。schema v2 では
`field_reconstruction_resolved_field_solver` と `field_reconstruction_fmm_expansion_order` が必須です。
新しい出力を Python で電場再構成するときは、近くの `beach.toml` ではなく、この receipt を使います。

### `zhao_stationary`

モデルの物理的意味と適用限界は[Zhao stationary closure](ZhaoStationaryClosure.html)を参照してください。
`[surface_current_model]` の選択結果は `surface_current_model` に記録されます。`zhao_stationary` では、
選択 branch、参照面積、$\phi_0$、$\phi_m$、解いた ambient electron 密度、および electron / ion /
PE emission / PE escape / PE return / net の signed 電流密度を `surface_current_model_*` で出力します。

`surface_current_model_photoelectron_active` は PE channel の有無を示し、false では PE 関連 receipt は 0 です。
`surface_current_model_pe_return_current_density_A_m2` は負、emission と escape は正です。
`surface_current_model_pe_escape_particle_current_A` は、外向き PE 粒子が運ぶ signed 電流なので負です。
二つの budget residual は、PE の emission-return-escape 連続条件と表面定常電流を独立に検証します。

`surface_current_model_kinetic_contract` と inflow access / outflow barrier の電位・face receipt は、
Zhao 電流に対応する速度空間境界写像を記録します。face 番号は
`1..6 = x_low, x_high, y_low, y_high, z_low, z_high` です。

### `matching_plane_quasistatic`

`matching_plane_quasistatic` では、静的な `surface_current_model_*` 電流 target ではなく、accepted batch の
`matching_plane_*` を読みます。`matching_plane_displacement_C_m2` と `matching_plane_phi_V` が電磁気的な
整合値、electron / ion の inward flux、access potential、PE barrier potential が選択した response backend の
出力です。4 つの outward flux / energy は固定点 feedback、`matching_plane_iterations` と
`matching_plane_residual` は数値的な収束 receipt です。

`matching_plane_photoelectron_return_flux_m2_s` と
`matching_plane_photoelectron_escape_flux_m2_s` は、同じ batch の outward flux との収支確認に使います。
共通の run provenance は、`surface_current_model_response_backend`、
`surface_current_model_matching_plane_z_m`、`surface_current_model_electron_species`、
`surface_current_model_ion_species`、`surface_current_model_photoelectron_species` です。

反復条件は `surface_current_model_coupling_rtol`、`surface_current_model_coupling_atol`、
`surface_current_model_coupling_max_iterations`、`surface_current_model_coupling_relaxation`、状態の生成元は
`surface_current_model_dynamic_state_source=accepted_batch_fixed_point` に記録されます。
`surface_current_model_coupling_atol` の 4 値は、順に PE 外向き flux [m^-2 s^-1]、PE 平均法線 energy [eV]、
electron 外向き flux [m^-2 s^-1]、ion 外向き flux [m^-2 s^-1] です。既定値はすべて 0 で、inactive 成分も
0 でなければなりません。active 成分の判定閾値は `max(coupling_rtol * backend_scale, coupling_atol)` です。
絶対許容値が支配する成分は有効残差へ換算されるため、accepted state の `matching_plane_residual` は
引き続き `surface_current_model_coupling_rtol` 以下になります。

table backend だけが `surface_current_model_response_table_path` と
`surface_current_model_response_content_fingerprint` を持ち、後者は読込済み応答表の canonical 内容を
識別します。online backend は次を記録します。

- `surface_current_model_response_contract=matching_plane_zhao_online_v1`
- `surface_current_model_zhao_branch`
- `surface_current_model_outer_solver=charge_driven_finite_h_sagdeev`
- `surface_current_model_photoelectron_closure=moment_matched_half_maxwellian`
- `surface_current_model_ambient_outward_feedback=transparent`
- `surface_current_model_outer_solver_state=stateless`

`matching_plane_state_valid` が false なら、同じ summary 内の `matching_plane_*` 値を accepted state として
使ってはいけません。各行の $D_H$ と $\Phi_H$ は、その batch の粒子追跡に使った commit 前の表面電荷 state に
対応します。`simulated_time_s` は trial を受理して進めた後の時刻なので、次 batch 開始時の post-commit 場と
取り違えないでください。

## 履歴

| ファイル | 条件 | 主な列 |
| --- | --- | --- |
| `charge_history.csv` | `history_stride > 0` | batch、処理粒子数、相対変化、要素、電荷 |
| `potential_history.csv` | `write_potential_history=true` かつ `history_stride > 0` | batch、要素、電位 |
| `top_reference_history.csv` | 上記かつ `[domain]` の box あり | batch、時間、z-high 面の電位統計 |
| `matching_plane_history.csv` | matching-plane、stride 有効 | batch、時間、$D_H$、$\Phi_H$、flux、反復 |

Python では `load_fortran_result(...)` の `matching_plane_state` と `matching_plane_history` から、列番号を
手作業で管理せず typed receipt として参照できます。

`top_reference_history.csv` の基準は box の z-high 面平均です。無限遠電位や plasma 電位ではありません。
要素相対電位は、同じ batch の `potential_history.csv` と結合し、
`potential_V - potential_mean_V` として求めます。

## mesh 電位

`write_mesh_potential=true` では `mesh_potential.csv` を出力します。

- 各要素重心の電位 [V]
- 自己項は解析的 P0 triangle panel kernel
- finite periodic2 は指定 image shell を加算
- `cached_kneq0` は cached 非零 mode と設定された物理的 zero mode を合成

## charge ledger

`charge_ledger.csv` は species ごとに次の 25 列を記録します。

| 列 | 意味 |
| --- | --- |
| `batch` | ledger state が対応する batch |
| `species_idx` | 1 始まりの species index |
| `injected_from_remote_C` | 外部から注入された signed charge [C] |
| `emitted_from_surface_C` | 表面から放出された signed charge [C] |
| `absorbed_on_surface_C` | 表面へ吸収された signed charge [C] |
| `escaped_to_infinity_C` | 軌道追跡で得た raw escape charge [C] |
| `discarded_unresolved_C` | 未解決 discard の signed charge [C] |
| `neutral_return_correction_C` | closed PE の neutral-return 補正 [C] |
| `neutral_return_weight_scale` | neutral-return に使った weight scale |
| `neutral_return_unresolved_fraction` | neutral-return の未解決率 |
| `fixed_absorbed_target_charge_C` | fixed-current の吸収 target [C] |
| `fixed_absorbed_weight_scale` | 吸収標本へ適用した weight scale |
| `fixed_emission_target_charge_C` | fixed-current の放出 target [C] |
| `fixed_emission_weight_scale` | 放出標本へ適用した weight scale |
| `fixed_current_correction_C` | 外部 closure が表面電荷へ加えた補正 [C] |
| `fixed_absorbed_applied_charge_C` | 吸収 target と同値の互換 alias |
| `fixed_emission_applied_charge_C` | 放出 target と同値の互換 alias |
| `fixed_escape_target_charge_C` | 外部 escape target [C] |
| `fixed_escape_applied_charge_C` | escape target と同値の互換 alias |
| `fixed_escape_correction_C` | target と raw escape の差 [C] |
| `injected_count` | 注入したマクロ粒子数 |
| `emitted_count` | 表面から放出したマクロ粒子数 |
| `absorbed_count` | 表面へ吸収したマクロ粒子数 |
| `escaped_count` | escape したマクロ粒子数 |
| `discarded_unresolved_count` | 未解決 discard のマクロ粒子数 |

`escaped_to_infinity_C` は軌道追跡で得た raw 値です。Zhao fixed-current ではこれを上書きせず、外部 escape
target と並べて保存します。escape 補正は表面要素へ deposit されないため、表面電荷更新と外部境界収支を
分けて解析できます。`*_applied_charge_C` は既存 reader 向けの出力 alias で、内部 ledger は
`*_target_charge_C` だけを保持します。

`summary.txt` の `charge_ledger_residual_C` は surface / local-flight / unresolved stock の変化と外部 flux、
neutral-return 補正、fixed-current 補正から作る保存残差です。対応する before / after stock は
`charge_ledger_surface_charge_*_C`、`charge_ledger_local_flight_charge_*_C`、
`charge_ledger_unresolved_stock_*_C` に記録されます。

`fixed_current` の統計診断では、対象 channel の `absorbed_count` または `emitted_count`、対応する raw charge、
`fixed_*_weight_scale` を組にして読みます。count が 1 なら target 全量が 1 標本の要素へ局所化されます。
保存残差が小さくてもこの標本分散は検出できないため、粒子数、ray 数、batch 幅、乱数 seed を変えた
要素別電荷分布の収束を別に確認します。

closed PE では raw の吸収・未解決量を上書きしません。補正量と係数を別に記録するため、表面総電荷が
閉じていても未解決率を独立に確認できます。公式ケースの基本的な照合手順は
[charge ledger を確認する](OutputGuide.html#charge-ledger)を参照してください。

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
<output.dir>/
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
| `macro_residuals.csv` | global な macro 粒子数残差。species × face を区別する |
| `charge_ledger.csv` | summary に ledger metadata がある場合に復元 |
| `checkpoint_complete.txt` | schema v8 以降の checkpoint 一式を最後に公開する完了 manifest |

`output.restart_from` は checkpoint の読み込み元だけを変更します。新しい出力は `output.dir` に書きます。
直下の最終出力と両 periodic slot のうち、必須ファイルが揃う load 可能な checkpoint を比較し、`batches` が
最大のものを自動選択します。`checkpoint_latest.txt` が欠落、破損、または古い場合も、完了 manifest を持つ
slot は回収対象です。必須ファイルの欠落、fingerprint、mesh 要素数、species 数、MPI world size の不一致では、
新規実行へ fallback せず停止します。

schema v8 以降では、再開状態を書き始める前に `checkpoint_complete.txt` を `in_progress` にし、summary、
charges、全 rank の RNG、manifest で宣言した residual と ledger を閉じてから `complete` を原子的に
公開します。この schema では `checkpoint_complete.txt` 自体が再開の必須ファイルです。直下の最終出力が
中断されて新旧世代のファイルが残った場合は、その出力を選ばず完全な periodic slot へ戻ります。

checkpoint schema v6 の `macro_residuals.csv` は `species_idx,face,residual` です。`face=0` は従来 source、
`1..6` は boundary face を表します。旧 `species_idx,residual` の 2 列形式も読み込めます。
schema v9 は matching-plane の accepted feedback と反復 receipt を `summary.txt` へ保存します。
model fingerprint には、table backend では response table の canonical 内容、online backend では Zhao contract と
branch policy が含まれます。これらを変えた checkpoint は再開できません。

実際の再開操作は[実行・再開する](Execution.html)を参照してください。

## Python から読む

```python
from beach import load_fortran_result

result = load_fortran_result("outputs/tutorial")
print(result.charges)
print(result.charge_ledger)
print(result.field_reconstruction)
print(result.matching_plane_state)
```

`summary.txt` と `charges.csv` は必須です。mesh、履歴、ledger、field reconstruction、matching-plane state は
対応する出力がある場合だけ属性へ読み込まれます。`matching_plane_state_valid` が false の場合、reader は
`matching_plane_state` と `matching_plane_history` を accepted state として公開しません。

別のケースでは `outputs/tutorial` を実際の `output.dir` に置き換えます。class と属性の詳細は
[Python 後処理 API](PythonPostprocessAPI.html)、物理・数値的な受理判断は
[計算結果の妥当性確認](ValidationGuide.html)を参照してください。

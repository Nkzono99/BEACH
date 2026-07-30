title: 出力ファイルを調べる

Lang: [日本語](OutputGuide.md) | [English](OutputGuide.en.md)

# 出力ファイルを調べる

このページは、BEACH の出力ディレクトリから目的の値を見つけるためのガイドです。
計算を受理できるかの判定は[計算結果の妥当性確認](ValidationGuide.html)、図やアニメーションの作成は
[後処理チュートリアル](PostprocessTutorial.html)に分けています。

## 最初に見る 4 つの出力

```bash
beachx inspect outputs/latest
```

`output.dir` を変更した場合は、`outputs/latest` を実際の出力先に置き換えます。

| 知りたいこと | 読む出力 | 主な内容 |
| --- | --- | --- |
| 何バッチ・何秒進み、粒子がどう終了したか | `summary.txt` | `batches`, `simulated_time_s`, `absorbed`, `escaped`, `survived_max_step` |
| 粒子種ごとの電荷がどこから入り、どこへ出たか | `charge_ledger.csv` | 注入、放出、吸収、脱出、未解決破棄 |
| 最終的な表面電荷がどこに分布したか | `charges.csv` + `mesh_triangles.csv` | 要素番号、電荷、三角形座標、`mesh_id` |
| バッチごとの変化を見られるか | `charge_history.csv`, `potential_history.csv` | 要素電荷・要素重心電位の履歴 |

`beachx inspect` が表示するのは実行概要です。完了、電荷保存、時間収束、モデルの適用範囲を判定する手順は
[計算結果の妥当性確認](ValidationGuide.html)を使います。

## 出力が更新される条件

このページは `output.write_files = true` を前提とします。`false` の実行では、BEACH は `output.dir` の
`summary.txt` や `charges.csv` を作成・更新しません。ディレクトリに既存ファイルがあっても、それは過去の実行結果です。

生成条件と再開時の役割は、機械可読な
[`beach.output-manifest.json`](../schemas/beach.output-manifest.json)を正本とします。

## `summary.txt`: 実行全体の概要

`summary.txt` は `key=value` 形式です。最初は次のグループだけを探せば十分です。

| グループ | 主なキー | 意味 |
| --- | --- | --- |
| 実行識別 | `build_version`, `build_source_commit`, `model_fingerprint`, `mesh_fingerprint`, `species_fingerprint` | どのビルド・設定・メッシュ・粒子種で作った出力か |
| 規模 | `mesh_nelem`, `mesh_count`, `mpi_world_size` | 要素数、メッシュ数、MPI ランク数 |
| 粒子処理 | `processed_particles`, `absorbed`, `escaped`, `escaped_boundary`, `survived_max_step`, `multiple_box_events_soft_discarded`, `multiple_box_events_soft_discarded_abs_charge_C` | 粒子イベントの集計数と記録付きsoft discard |
| 進行 | `batches`, `simulated_time_s`, `last_rel_change` | accepted batch数、受理した物理時間の総和、最終バッチの電荷変化監視値 |
| 適応的な非零モード進行 | `periodic2_max_nonzero_mode_potential_step_V`, `adaptive_nonzero_mode_rejected_trials`, `adaptive_nonzero_mode_last_batch_duration_s`, `adaptive_nonzero_mode_last_potential_step_V`, `adaptive_nonzero_mode_omp_threads` | 設定上限、累積棄却trial数、最後に受理した幅・実測$k\ne0$電位変化、再開に必要なthread数 |
| 場計算 | `field_backend`, `field_source_model`, `field_kernel_id` | 出力を作った場ソルバーと source kernel |
| z-high電位基準 | `top_reference_available`, `top_reference_definition`, `top_reference_last_batch`, `top_reference_simulated_time_s`, `top_reference_z_high_m`, `top_reference_sample_n`, `top_reference_potential_mean_V`, `top_reference_potential_std_V`, `top_reference_potential_min_V`, `top_reference_potential_max_V` | 最終mesh電位と同じ状態で評価した全z-high面の定義・時刻・標本数・統計 |
| 光電子帰還closure | `charge_ledger_neutral_return_correction_C` | 全speciesで解決済み帰還先へ追加した有効電荷の累積和 |
| 外部境界の解決結果 | `coupling_update_mode`, `external_inflow_map`, `external_ordinary_open_model`, `external_interface_transport`, `outer_particle_mode_resolved` | 公開設定から実行時に解決された更新、流入、通常open面、z-high輸送、処理時機 |

`absorbed` は吸収イベント数であり、電荷の符号やマクロ粒子重みを含みません。電荷量は
`charge_ledger.csv`、最終分布は `charges.csv` から読みます。

`last_rel_change` は監視値です。現行実装では早期停止条件ではありません。

`periodic2_max_nonzero_mode_potential_step_V=0`は適応進行が無効であることを表します。正の場合、
`sim.batch_duration`は最大trial幅であり、`batches`はaccepted batchだけを数えます。
`simulated_time_s`は受理したtrial幅の累積値なので、一般には
`batches * sim.batch_duration`と一致しません。棄却trialはRNG、macro粒子数残差、outer/mean transactionを
rollbackされ、粒子集計、`charge_history.csv`、`potential_history.csv`、`charge_ledger.csv`へ現れません。
`adaptive_nonzero_mode_last_potential_step_V`は全panel重心で測った最後のaccepted trialの最大
$k\ne0$電位変化であり、局所打切り誤差ではありません。
`adaptive_nonzero_mode_rejected_trials`は共通ladder上の総棄却数なので、$k\ne0$上限超過に加えて
`implicit_mean` Zhaoの有限なambientまたは界面電場・障壁のtrust-region棄却も含みます。
時間幅で縮まらない計測source規格化の変化と、その他の$k=0$ closure failureは
再試行せず停止するため、この値には含まれません。内訳は実行ログで確認します。
適応runの再開は、加算順と受理ladderを再現するため
`adaptive_nonzero_mode_omp_threads`と現在の実OpenMP team sizeが全MPI rankで一致する場合だけ受理します。
適応進行が無効なrunでは、この値は`0`です。

実行中の標準出力では、棄却ごとに`BEACH adaptive-kneq0 reject`行がbatch、halving回数、trial幅を
記録します。$k\ne0$上限による棄却では実測電位変化と設定上限、`implicit_mean` Zhaoの
回復可能なtrust-region棄却では`implicit_status`と`reason`を記録するため、両者を区別できます。受理時の
`BEACH adaptive-kneq0 accept`行は、受理後の
`time_end_s`、trial幅、実測電位変化、halving回数を記録します。

外部境界の5キーは設定項目ではなく、`[external_boundary]` facadeが実行時にどう解決されたかを確認する
receiptです。

| キー | 出力値 |
| --- | --- |
| `coupling_update_mode` | `explicit` / `implicit_mean` |
| `external_inflow_map` | `source_vdf` / `infinity_barrier` / `kinetic_profile` |
| `external_ordinary_open_model` | `escape` / `potential_barrier` |
| `external_interface_transport` | `none` / `kinetic_1d` |
| `outer_particle_mode_resolved` | `local_source` / `same_batch` / `zhao_queue` |

たとえば`particles.inflow_model="auto"`は、fieldとparticle modeに応じて
`external_inflow_map=source_vdf`または`kinetic_profile`へ解決されます。
`ambient_linear_debye + same_batch` と負電荷 `photo_raycast` は
`coupling_update_mode=implicit_mean`へ自動解決され、公開TOMLへこの内部値は書きません。
再現性を確認するときは、入力facadeとこのreceiptを併せて保存します。

陰的平均更新の実行中は、標準出力の `BEACH implicit-mean` 行に batch ごとの表面総電荷、
interface 電位・電場、species 電流、追加の実測表面電流 `J_other_A_m2`、
零和transactionの `transaction_residual_C`、scalar根の `mean_solver_iterations`、
ray標本の `sample_escape_fraction`、`return_weight_scale` を記録します。
ambient linear-Debye closureでは、前者は正規化前のmacro電荷標本比、後者は
$R_{\mathrm{analytic}}/R_{\mathrm{sample}}$ です。後者は確率ではなく1を超える場合があります。
`mean_solver_iterations` はこのclosureではscalar root反復数です。

Zhao closureでは、続く `BEACH implicit-zhao` 行にbranch、virtual cathode電位 `phi_min_V`、
全profileの `barrier_J`、実測界面電流から決めた `source_scale`、marginal energyとそのescape fraction、
非線形電荷残差、再越境charge fraction、終端不一致charge fractionを記録します。この経路はenergy群ごとの
実測CDF重みを直接使うため、`return_weight_scale=1`です。Zhaoで共有表示される
`mean_solver_iterations` は、endpoint certificate、order-statistic二分探索、marginal bisectionが呼んだ
connected candidate solve数であり、
pseudo-arclength内部のNewton反復数ではありません。どちらの進行行も収束判定そのものではありません。
`summary.txt`、`charge_ledger.csv`、`charges.csv`と合わせて確認します。

## `charge_ledger.csv`: 粒子種別の電荷移動

`charge_ledger.csv` は 1 粒子種につき 1 行を持ち、符号付き電荷量 [C] とイベント数を分けて記録します。

| 列 | 記録する電荷 |
| --- | --- |
| `injected_from_remote_C` | `volume_seed` または `reservoir_face` から入った電荷 |
| `emitted_from_surface_C` | `photo_raycast` で表面から放出された追跡対象の電荷 |
| `absorbed_on_surface_C` | メッシュへ吸収された電荷 |
| `escaped_to_infinity_C` | 開放境界または外部モデルから無限遠へ出た電荷 |
| `discarded_unresolved_C` | `max_step` 到達時に生存したまま破棄された電荷 |
| `interface_outward_gross_C` | 局所領域から外部領域へ渡した総電荷 |
| `interface_returned_gross_C` | 外部領域から局所領域へ戻した総電荷 |
| `neutral_return_correction_C` | `neutral_return`が解決済み帰還先へ追加した有効電荷。raw吸収量には含めない |

`neutral_return_weight_scale`は解決済み帰還depositへ掛けた係数、`neutral_return_unresolved_fraction`は
放出電荷のうち`max_step`まで未帰還だった絶対量の割合です。closureを使わないspeciesではそれぞれ1と0です。
`absorbed_on_surface_C`と`discarded_unresolved_C`はtracked raw値のままなので、補正の大きさを独立に監査できます。
最終CSVの係数と割合は、そのrunの累積raw放出・吸収・未解決電荷から再計算した値です。

`injected_count`、`emitted_count`、`absorbed_count`、`escaped_count`、
`discarded_unresolved_count` は対応するterminalイベント数です。interface grossに独立したcount列はありません。
`summary.txt` の `charge_ledger_residual_C` と
`charge_ledger_discarded_unresolved_abs_C` は全粒子種を集約した値です。許容値の決め方は
[粒子数と電荷の収支を確認する](ValidationGuide.html#2-粒子数と電荷の収支を確認する)で扱います。

`coupling_update_mode=implicit_mean` の光電子では、電荷量とcountを別々に解釈します。
置換するのは、z-highを外向きに横切ってdeferredされた成分だけです。ambient linear-Debye closureでは、この成分の
`escaped_to_infinity_C` とreturn後の `absorbed_on_surface_C` 寄与をscalar closureの
continuous Maxwellian escape / return総量に合わせ、再吸収電荷を1回追跡したrayの
実hit分布へ `return_weight_scale` で正規化配分します。Zhao closureでは、実測界面energy CDFと非線形
$Q(\Phi_I)$解から各rayのescape / return重みを直接決めるため、共通のanalytic scaleを掛けません。interface到達前の局所再吸収、
z-lowなど別のopen面からのescape、`discarded_unresolved_C` はtracked値を保持します。
したがって `escaped_to_infinity_C` はtrackedな別open面escapeとclosureが決めたz-high escapeの和、
`absorbed_on_surface_C` はtrackedな局所再吸収とclosureが決めたz-high return後再吸収の和です。
`escaped_count` と `absorbed_count` のinterface-return由来分はray標本のterminal分類であり、
電荷比のsource of truthではありません。解析置換の対象外にあるterminal countは、通常どおり
tracked粒子が更新します。

`interface_outward_gross_C` はlocal領域からz-highを外向きに横切った実crossing、
`interface_returned_gross_C` は1D領域からlocal領域へ戻ったcrossingの符号付きmacro電荷を累積します。
各rayのcrossingはterminal return / escape総量に合うよう正規化したweightで集計し、
z-high deferred成分のclosure escape符号付き電荷を
$Q_{\mathrm{esc,z-high}}^{\mathrm{closure}}$ とすると、
`interface_returned_gross_C` = `interface_outward_gross_C` - $Q_{\mathrm{esc,z-high}}^{\mathrm{closure}}$
の電荷恒等式を保ちます。別open面からのtracked escapeも含むterminal総量
`escaped_to_infinity_C` を、この式の右辺へそのまま代入してはいけません。
ambient linear-Debye経路ではreturn後に再びz-highを横切るrayもgross列へ追加され得ます。Zhao実測CDF経路では、
この再越境を適用性違反として有意chargeがあれば停止します。どちらの場合もgross量はterminal量ではありません。

Zhao過渡queueでは、`charge_ledger_outer_flight_charge_before_C`と
`charge_ledger_outer_flight_charge_after_C`がbatch前後のouter-flight stockを記録します。

## 表面分布とメッシュを対応させる

| ファイル | 1 行の単位 | 主な列 |
| --- | --- | --- |
| `charges.csv` | 1 三角形要素 | `elem_idx`, `charge_C` |
| `mesh_triangles.csv` | 1 三角形要素 | `elem_idx`, 三頂点座標, `charge_C`, `mesh_id` |
| `mesh_sources.csv` | 1 メッシュ | `mesh_id`, `source_kind`, `surface_model`, 要素数 |
| `mesh_potential.csv` | 1 三角形要素 | `elem_idx`, `potential_V` |

`elem_idx` で `charges.csv` と `mesh_triangles.csv` を対応付けます。複数メッシュを含む場合は、
`mesh_triangles.csv` の `mesh_id` と `mesh_sources.csv` を使って対象を分けます。

`mesh_potential.csv` は `output.write_mesh_potential = true` の場合だけ生成され、最終時刻の要素重心電位を持ちます。
出力済みの電位をそのまま使う場合はこのファイルを読み、任意点で再評価する場合はnative P0 panel kernelを使う
`compute_potential_mesh` / `compute_potential_points`を使います。

## 履歴ファイル

履歴を出すには、実行前に `output.history_stride` を正の値にします。

```toml
[output]
history_stride = 1
write_potential_history = true
```

| ファイル | 生成条件 | 内容 |
| --- | --- | --- |
| `charge_history.csv` | `output.history_stride > 0` | 記録バッチごとの要素電荷 |
| `potential_history.csv` | 上記に加えて `output.write_potential_history = true` | 記録バッチごとの要素重心電位 |
| `top_reference_history.csv` | 上記に加えて `sim.use_box = true` | 同じbatchのz-high面電位統計 |

受理済み batch 1 は常に含まれ、その後は受理済み batch の番号が `history_stride` ごとに記録されます。
`top_reference_history.csv`は
`batch,simulated_time_s,z_high_m,sample_n,potential_mean_V,potential_std_V,potential_min_V,potential_max_V`
を持ちます。`potential_history.csv`と`batch`でjoinし、
`potential_V - potential_mean_V`をz-high面基準の相対電位として読みます。この平均は無限遠電位や
プラズマ電位ではありません。

旧checkpointを、`top_reference_history.csv`がまだない同じ出力先へresumeした場合、このファイルは
resume後のbatchから始まります。二つの履歴に共通する`batch`だけをinner joinします。

電位履歴は記録時に追加の場評価を行うため、
大きなメッシュでは `history_stride` を増やして出力頻度を下げます。

履歴の図とアニメーションは[後処理チュートリアル](PostprocessTutorial.html#履歴アニメーション)にあります。

## その他のファイルを目的から探す

| 目的 | ファイル | 生成条件または役割 |
| --- | --- | --- |
| 実行時間を分解する | `performance_profile.csv` | `BEACH_PROFILE=1` |
| 外部シースの格子 profile を読む | `outer_plasma_profile.csv` | `kinetic_1d` の外部状態が準備済み |
| Zhao過渡eventを再開する | `outer_event_queue.csv` / `outer_event_queue_rankNNNNN.csv` | `external_boundary.particles.mode="zhao_queue"`。serialでは前者、MPIではrankごとに後者 |
| 乱数状態から再開する | `rng_state.txt` / `rng_state_rankNNNNN.txt` | serial では前者、MPI ではランクごとに後者 |
| 端数マクロ粒子を復元する | `macro_residuals.csv` | 残差状態が割り当てられた場合。MPI でも 1 ファイル |

すべての生成には `output.write_files = true` が必要です。

## 構成固有の値を探す

高度な構成の判定基準は各モデルのページに置き、このページでは出力先だけを示します。

| 構成 | 主な出力 | 詳細 |
| --- | --- | --- |
| 有限画像 `periodic2` | periodic2構成、`charges.csv`、`top_reference_history.csv`、neutral-return台帳列 | [有限画像構成](FinitePeriodicConfiguration.html) |
| `cached_kneq0` | `periodic2_cache_hit`, `periodic2_operator_build_count`, `periodic2_cache_fingerprint`, `periodic2_cache_path` | [周期遠方補正](PeriodicFarCorrection.html) |
| 適応的な$k\ne0$進行 | `simulated_time_s`, `adaptive_nonzero_mode_rejected_trials`, `adaptive_nonzero_mode_last_batch_duration_s`, `adaptive_nonzero_mode_last_potential_step_V` | [`batch_duration`の安定性と定常値](BatchDurationStability.html#periodic2非零モードの適応的な進行) |
| `kinetic_1d` | `outer_plasma_profile.csv`, `interface_potential_V`, `gauss_residual_C`, `last_outer_update_batch` | [標準 1D kinetic 外部シース](KineticOuterPlasma.html) |
| Zhao過渡queue | `outer_event_queue.csv` / `outer_event_queue_rankNNNNN.csv`, `outer_photoelectron_population_fraction`, `outer_photoelectron_column_per_area_m2`, `outer_photoelectron_column_target_per_area_m2`, `outer_photoelectron_column_residual_per_area_m2`, `outer_queue_event_count`, `outer_queue_signed_charge_C`, `outer_queue_fingerprint` | [粒子の escape と return](ParticleEscapeReturn.html#zhao-過渡closureでouter-flightをqueueする) |
| 外部粒子移送 | `interface_outward_gross_C`, `interface_returned_gross_C`, `max_outer_flight_time_s`, `max_outer_frozen_field_ratio`, `max_outer_energy_relative_error` | [粒子の escape と return](ParticleEscapeReturn.html) |
| `implicit_mean` return shadow | `implicit_mean_last_returned_outer_flight_time_mean_s`, `implicit_mean_last_estimated_returning_photoelectron_column_charge_per_area_C_m2` | [標準 1D kinetic 外部シース](KineticOuterPlasma.html) |

`outer_infinity_potential_V` は内部の無限遠 gauge を記録する診断値であり、入力 key ではありません。
現行の kinetic 状態では 0 に固定されます。
`max_outer_energy_relative_error` は、kinetic 1D の return / escape 写像で法線運動エネルギーと静電エネルギーの
保存残差を規格化した値の最大です。
`implicit_mean` 光電子も個別の1D profile写像を使うため、outer flight time、frozen-field比、
energy保存残差に寄与します。return後のrecrossも同じ診断へ累積されます。この光電子は準定常shadowなので、
`max_outer_frozen_field_ratio` が設定上限を超えても、それ自体は run の失敗ではありません。比は shadow 軌道の
時間 scale を示しますが、UV turn-on の遅延 return current を解いたことは意味しません。通常の `same_batch` 粒子と
ambient species は、上限超過で従来どおり fail closed にします。

2 つの `implicit_mean_last_*` 値は最大値や累積値ではなく、今回の実行で最後に完了した batch だけを表します。
batch を進めなかった no-op resume では両 key を出力しません。
前者は analytic weight 適用後の return excursion を正の電荷 magnitude で重み付けした平均 outer 往復時間、
後者は同じ return excursion の $\sum W_j\tau_j/(A\Delta t)$ から得る Little の法則による正の
photoelectron column shadow 推定値です。後者は実 queue / ledger stock ではありません。
`escaped_to_infinity` outcome 自体には滞在時間を加えませんが、最終 escape 前に完了した return excursion は集計します。
`outer_integrated_charge_per_area_C_m2` が非零なら、shadow 推定値をその絶対値で割った
$\chi_{\mathrm{PE,shadow}}$ は、省略した returning-photoelectron column と 1D outer 積分電荷の相対 scale を示します。
これは解釈用の比であり、組み込みの受理閾値ではありません。

`dielectric` 要素がある場合、`summary.txt` は `surface_model_dielectric_elem_count` と
`surface_model_note=metadata_only_dielectric_present` を記録します。現行実装の `dielectric` はメタデータであり、
誘電体境界条件を解いたことを意味しません。

## 再開に使うファイル

再開用ファイルは解析出力と同じディレクトリにありますが、個別に編集しません。

| チェックポイント | 再開時の役割 |
| --- | --- |
| `summary.txt`, `charges.csv` | 常に必須 |
| `rng_state.txt` / `rng_state_rankNNNNN.txt` | serial では前者、MPI では全ランク分が必須 |
| `macro_residuals.csv` | 存在する場合に全体の端数状態を復元 |
| `charge_ledger.csv` | `summary.txt` に台帳メタデータがある場合に必須 |
| `outer_plasma_profile.csv` | 準備済みの `kinetic_1d` 状態を再開する場合に必須 |
| `outer_event_queue.csv` / `outer_event_queue_rankNNNNN.csv` | Zhao過渡queue有効時に必須。serialでは前者、MPIでは同じworld sizeの全rank分が必要 |

`output.restart_from` を指定すると、チェックポイントを `restart_from` から読み、新しい出力を `output.dir` に書きます。
再開手順と fingerprint の照合は[実行する](Execution.html#再開実行)にまとめています。
`summary.txt`から`simulated_time_s`、累積棄却trial数、最後のaccepted trial幅・電位変化、
適応進行時の実OpenMP team sizeも復元するため、
再開後の`batch_count`は引き続きaccepted batchの累積到達数です。
queue checkpointはactive phase-space record、terminal outcome、due時刻、`next_event_id`を保持し、schema、rank、
world size、完了batchの不一致をfail closedで拒否します。

## 次に読むページ

- 数値・物理的に結果を受理できるか判断する: [計算結果の妥当性確認](ValidationGuide.html)
- 図、アニメーション、Python 解析を作る: [後処理チュートリアル](PostprocessTutorial.html)
- `[output]` の設定値を調べる: [入力パラメータリファレンス](Parameters.html#output-出力と再開)

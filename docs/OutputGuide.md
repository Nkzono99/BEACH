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
| 何バッチ実行し、粒子がどう終了したか | `summary.txt` | `batches`, `absorbed`, `escaped`, `survived_max_step` |
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
| 進行 | `batches`, `last_rel_change` | 完了バッチ数と最終バッチの電荷変化監視値 |
| 場計算 | `field_backend`, `field_source_model`, `field_kernel_id` | 出力を作った場ソルバーと source kernel |

`absorbed` は吸収イベント数であり、電荷の符号やマクロ粒子重みを含みません。電荷量は
`charge_ledger.csv`、最終分布は `charges.csv` から読みます。

`last_rel_change` は監視値です。現行実装では早期停止条件ではありません。

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

対応する `*_count` 列はイベント数です。`summary.txt` の `charge_ledger_residual_C` と
`charge_ledger_discarded_unresolved_abs_C` は全粒子種を集約した値です。許容値の決め方は
[粒子数と電荷の収支を確認する](ValidationGuide.html#2-粒子数と電荷の収支を確認する)で扱います。
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
`field_source_model="triangle_p0"` の電位を調べる場合は、Python 側で点電荷として再構成せず、このファイルを使います。

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

バッチ 1 は常に含まれ、その後は `history_stride` ごとに記録されます。電位履歴は記録時に追加の場評価を行うため、
大きなメッシュでは `history_stride` を増やして出力頻度を下げます。

履歴の図とアニメーションは[後処理チュートリアル](PostprocessTutorial.html#履歴アニメーション)にあります。

## その他のファイルを目的から探す

| 目的 | ファイル | 生成条件または役割 |
| --- | --- | --- |
| 実行時間を分解する | `performance_profile.csv` | `BEACH_PROFILE=1` |
| 外部シースの格子 profile を読む | `outer_plasma_profile.csv` | `kinetic_1d` / `unified_linear_response` の外部状態が準備済み |
| 光電子のエネルギー分布を読む | `photoelectron_histogram.csv` | `outer_plasma.photoelectron_histogram_enabled=true` で状態が準備済み |
| Zhao過渡eventを再開する | `outer_event_queue.csv` / `outer_event_queue_rankNNNNN.csv` | `coupling.outer_queue_enabled=true`。serialでは前者、MPIではrankごとに後者 |
| 乱数状態から再開する | `rng_state.txt` / `rng_state_rankNNNNN.txt` | serial では前者、MPI ではランクごとに後者 |
| 端数マクロ粒子を復元する | `macro_residuals.csv` | 残差状態が割り当てられた場合。MPI でも 1 ファイル |

すべての生成には `output.write_files = true` が必要です。

## 構成固有の値を探す

高度な構成の判定基準は各モデルのページに置き、このページでは出力先だけを示します。

| 構成 | 主な出力 | 詳細 |
| --- | --- | --- |
| 有限画像 `periodic2` | `summary.txt` の periodic2 構成、`charges.csv` | [有限画像構成](FinitePeriodicConfiguration.html) |
| `cached_kneq0` | `periodic2_cache_hit`, `periodic2_operator_build_count`, `periodic2_cache_fingerprint`, `periodic2_cache_path` | [周期遠方補正](PeriodicFarCorrection.html) |
| `kinetic_1d` | `outer_plasma_profile.csv`, `interface_potential_V`, `gauss_residual_C`, `last_outer_update_batch` | [標準 1D kinetic 外部シース](KineticOuterPlasma.html) |
| Zhao過渡queue | `outer_event_queue.csv` / `outer_event_queue_rankNNNNN.csv`, `outer_photoelectron_population_fraction`, `outer_photoelectron_column_per_area_m2`, `outer_photoelectron_column_target_per_area_m2`, `outer_photoelectron_column_residual_per_area_m2`, `outer_queue_event_count`, `outer_queue_signed_charge_C`, `outer_queue_fingerprint` | [粒子の escape と return](ParticleEscapeReturn.html#zhao-過渡closureでouter-flightをqueueする) |
| `unified_linear_response` | 上記に加えて `outer_accessible_fraction_min`, `outer_accessible_fraction_max`, `outer_accessible_fraction_refinement_error` | [高度な粗面線形 screening](UnifiedLinearResponse.html) |
| 光電子 histogram | `photoelectron_histogram.csv`, `photoelectron_previous_charge_ratio`, `photoelectron_linear_applicability_status` | [光電子の放出とライフサイクル](PhotoelectronEmission.html) |
| 外部粒子移送 | `interface_outward_gross_C`, `interface_returned_gross_C`, `max_outer_flight_time_s`, `max_outer_frozen_field_ratio`, `max_outer_energy_relative_error` | [粒子の escape と return](ParticleEscapeReturn.html) |

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
| `outer_plasma_profile.csv` | 準備済みの `kinetic_1d` / `unified_linear_response` 状態を再開する場合に必須 |
| `photoelectron_histogram.csv` | 光電子 histogram 状態を再開する場合に必須 |
| `outer_event_queue.csv` / `outer_event_queue_rankNNNNN.csv` | Zhao過渡queue有効時に必須。serialでは前者、MPIでは同じworld sizeの全rank分が必要 |

`output.restart_from` を指定すると、チェックポイントを `restart_from` から読み、新しい出力を `output.dir` に書きます。
再開手順と fingerprint の照合は[実行する](Execution.html#再開実行)にまとめています。
queue checkpointはactive phase-space record、terminal outcome、due時刻、`next_event_id`を保持し、schema、rank、
world size、完了batchの不一致をfail closedで拒否します。

## 次に読むページ

- 数値・物理的に結果を受理できるか判断する: [計算結果の妥当性確認](ValidationGuide.html)
- 図、アニメーション、Python 解析を作る: [後処理チュートリアル](PostprocessTutorial.html)
- `[output]` の設定値を調べる: [入力パラメータリファレンス](Parameters.html#output-出力と再開)

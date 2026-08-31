title: 出力形式リファレンス

# 出力形式リファレンス

Lang: [日本語](OutputReference.md) | [English](OutputReference.en.md)

BEACH の出力について、**どのファイルが生成されるか、どの列・receipt を読めばよいか、再開に何が必要か**を
調べるためのリファレンスです。receipt は、入力値そのものではなく、実行時に解決された設定や状態を
`summary.txt` に `key=value` 形式で残した記録です。

初めて結果を確認する場合は、先に[出力ファイルを調べる](OutputGuide.html)の手順を実行してください。
ファイル生成条件の機械可読な正本は `schemas/beach.output-manifest.json` です。

## ファイル生成条件

各行に、そのファイルの生成条件、形式、用途、reader または再開での役割をまとめています。複雑なファイルは
リンク先で同じ情報を再掲してから、列や判定規則を展開します。

| ファイル | 生成条件 | 形式と用途 | Python reader / 再開 |
| --- | --- | --- | --- |
| `summary.txt` | `output.write_files=true` | `key=value` receipt。[完了判定と主要 key](#構成固有の値を探す) | `FortranRunResult` の各属性。再開に必須 |
| `charges.csv` | `output.write_files=true` | `elem_idx`, `charge_C`。要素ごとの最終 signed charge | `result.charges`。再開に必須 |
| `mesh_triangles.csv` | `output.write_files=true` | `elem_idx`, `v0x`, `v0y`, `v0z`, `v1x`, `v1y`, `v1z`, `v2x`, `v2y`, `v2z`, `charge_C`, `mesh_id` | `result.triangles`, `result.mesh_ids`。再開には不使用 |
| `mesh_sources.csv` | OBJ または template mesh で `output.write_files=true` | `mesh_id`, `source_kind`, `template_kind`, `surface_model`, `epsilon_r`, `elem_count` | `result.mesh_sources`。再開には不使用 |
| `mesh_potential.csv` | `output.write_files=true` かつ `output.write_mesh_potential=true` | `elem_idx`, `potential_V`。[電位の構成](#mesh-電位) | `result.mesh_potential_v`。再開には不使用 |
| `charge_history.csv` | `output.write_files=true` かつ `output.history_stride>0` | `batch`, `processed_particles`, `rel_change`, `elem_idx`, `charge_C` | `result.history`。再開には不使用 |
| `potential_history.csv` | `output.write_files=true`、`output.write_potential_history=true`、`output.history_stride>0` | `batch`, `elem_idx`, `potential_V`。[基準電位との結合](#履歴) | `FortranRunResult` 専用属性なし。CSV を直接読む |
| `top_reference_history.csv` | `output.write_files=true`、`output.write_potential_history=true`、`output.history_stride>0`、`[domain]` の box あり | `batch`, `simulated_time_s`, `z_high_m`, `sample_n`, `potential_mean_V`, `potential_std_V`, `potential_min_V`, `potential_max_V` | `FortranRunResult` 専用属性なし。CSV を直接読む |
| `matching_plane_history.csv` | `output.write_files=true`、`output.history_stride>0`、`surface_current_model.model=matching_plane_quasistatic` | [17 列の accepted state](#matching_plane_quasistatic) | `result.matching_plane_history`。再開には不使用 |
| `charge_ledger.csv` | `output.write_files=true` かつ charge ledger state あり | [species ごとの 25 列](#charge-ledger) | `result.charge_ledger`。summary に ledger metadata があれば再開に必須 |
| `rng_state.txt` / `rng_state_rankNNNNN.txt` | `output.write_files=true`。serial は前者、MPI は rank ごとに後者 | RNG の内部状態 | `FortranRunResult` 属性なし。再開に必須 |
| `macro_residuals.csv` | `output.write_files=true` かつ macro 粒子数の残差 state あり | `species_idx`, `face`, `residual` | `FortranRunResult` 属性なし。schema v8 以降は manifest 宣言時に再開に必須 |
| `checkpoint_complete.txt` | `output.write_files=true`。final output と各 periodic slot | `key=value` 完了 manifest。[意味と選択規則](#再開に使うファイル) | `FortranRunResult` 属性なし。schema v8 以降は再開に必須 |
| `checkpoint_latest.txt` | `output.write_files=true` かつ `output.checkpoint_stride>0` で該当 accepted batch の commit 後 | periodic slot の案内 | reader の候補選択 hint であり、完全性の証拠ではない |
| `performance_profile.csv` | `BEACH_PROFILE=1` かつ `output.write_files=true` | コメント行後に `region`, `calls_sum`, `calls_mean`, `rank_min_s`, `rank_mean_s`, `rank_max_s`, `imbalance_ratio` | `beach-plot-performance-profile <output_dir>`。再開には不使用 |

Python reader の必須入力は `summary.txt` と `charges.csv` です。他のファイルは設定または runtime state に
応じて省略されます。

## 構成固有の値を探す

`summary.txt` は `output.write_files=true` で生成される `key=value` ファイルです。

`load_fortran_result(...)` は主要値を `FortranRunResult` の属性へ変換し、型付きでない receipt は
`beach.summary.load_summary_file(...)` から読めます。

### 完了判定

実行完了と checkpoint の完全性は別の条件です。

| 確認対象 | 条件 | この条件が保証しないこと |
| --- | --- | --- |
| 要求した batch の完了 | `summary.txt` の `batches == sim.batch_count` | 数値収束や物理的妥当性 |
| schema v8 以降の checkpoint 完全性 | `checkpoint_complete.txt` が `state=complete` で、manifest と必須ファイルが整合する | `sim.batch_count` までの実行完了 |
| Python reader での読込成功 | `summary.txt` と `charges.csv` が reader 契約を満たす | checkpoint 完全性や物理的妥当性 |

`checkpoint_complete.txt` は書き込み transaction の完了を示す manifest です。途中の periodic checkpoint も
`state=complete` になり得るため、実行目標の到達判定には使いません。数値・物理的な受理条件は
[計算結果の妥当性確認](ValidationGuide.html)を参照してください。現行 schema の正常終了出力として扱うには、
表の最初の 2 条件を両方確認します。

### `summary.txt` の共通 receipt

以下は結果の読解と再開判断に使う主要 receipt です。`summary.txt` のすべての内部診断 key を固定 API として
列挙する節ではありません。全 producer は [`bem_output_writer.f90`](../src/runtime/bem_output_writer.f90)、
Python が型付き属性へ変換する key は [`io.py`](../beach/fortran_results/io.py)を正本とします。

| receipt | 内容 |
| --- | --- |
| `mesh_nelem`, `mesh_count`, `mpi_world_size` | mesh と MPI 構成 |
| `batches`, `processed_particles`, `absorbed`, `escaped`, `escaped_boundary`, `survived_max_step` | 実行量と主要な粒子 outcome |
| `simulated_time_s`, `last_rel_change` | 受理済みの模擬時間と最終相対変化。`last_rel_change` は早期停止条件ではない |
| `checkpoint_schema_version`, `checkpoint_stride` | checkpoint の形式世代と解決済み出力間隔 |
| `model_fingerprint`, `mesh_fingerprint`, `species_fingerprint` | 再開時に照合する構成識別子 |
| `reservoir_inflow_map`, `particle_ordinary_open_model` | reservoir と通常 open 面の解決済み境界モデル |

build receipt は次の 5 key です。

| receipt | 内容 |
| --- | --- |
| `build_info_schema_version` | build receipt の schema |
| `build_version` | BEACH version |
| `build_version_mode` | version の解決方法 |
| `build_source_commit` | source commit |
| `build_id` | build 識別子 |

複数回の box 境界 event は次の 5 key で診断します。

| receipt | 内容 |
| --- | --- |
| `multiple_box_events_retry_attempted` | retry を試した粒子数 |
| `multiple_box_events_retry_resolved` | retry で解決した粒子数 |
| `multiple_box_events_soft_discarded` | soft-discard した粒子数 |
| `multiple_box_events_soft_discard_fraction` | `processed_particles` に対する soft-discard の割合 |
| `multiple_box_events_soft_discarded_abs_charge_C` | soft-discard した絶対電荷の総量 [C] |

### 場と境界の receipt

field-reconstruction schema v2 は、Python で実行時と同じ電場を再構成するための receipt です。新しい出力では、
近くに置かれた `beach.toml` から推測せず、`result.field_reconstruction` を使います。

| receipt | 内容 |
| --- | --- |
| `field_reconstruction_schema_version` | field-reconstruction receipt の schema。現行値は 2 |
| `field_reconstruction_resolved_field_solver` | 実際に選ばれた `direct` / `treecode` / `fmm` |
| `field_reconstruction_fmm_expansion_order` | FMM 展開次数。schema v2 では上の key とともに必須 |
| `field_reconstruction_field_bc_mode` | 解決済みの `free` / `periodic2` |
| `field_reconstruction_e0_V_m` | 解決済みの一様外部電場 `E0` |
| `field_reconstruction_use_box`, `field_reconstruction_box_min_m`, `field_reconstruction_box_max_m` | box の有無と範囲 |
| `field_reconstruction_boundary_low`, `field_reconstruction_boundary_high` | 6 面の解決済み境界 |
| `field_reconstruction_tree_theta`, `field_reconstruction_tree_leaf_max` | 実行時に使った tree 設定 |

periodic field の key は次のとおりです。物理的な意味と solver の選択は
[periodic2 の電場モデル](PeriodicElectrostatics.html)を参照してください。

| receipt | 内容 |
| --- | --- |
| `field_reconstruction_periodic_image_layers` | 近傍 image shell 数 |
| `field_reconstruction_periodic_far_correction` | far-correction model |
| `field_reconstruction_periodic_nonzero_mode_backend` | 非零 mode backend |
| `field_reconstruction_periodic_zero_mode_policy` | zero-mode policy |
| `field_reconstruction_periodic_lower_boundary_model` | lower-boundary model |
| `field_reconstruction_periodic_reference_mode_layers` | reference 非零 mode の image layers |
| `field_reconstruction_periodic_panel_quadrature_order` | panel quadrature order |
| `field_reconstruction_periodic_ewald_alpha` | Ewald splitting parameter |
| `field_reconstruction_periodic_ewald_layers` | Ewald layer 数 |
| `field_reconstruction_periodic_cache_dir` | operator cache directory |
| `field_reconstruction_periodic_generation_tolerance` | operator 生成 tolerance |

### `zhao_stationary`

`[surface_current_model]` で選ぶ物理モデル、Type A/B/C、適用限界は
[Zhao stationary closure](ZhaoStationaryClosure.html)を参照してください。
この節は receipt の読み方だけを定義します。

| receipt | 読み方 |
| --- | --- |
| `surface_current_model`, `surface_current_model_zhao_branch` | 選択された model と branch |
| `surface_current_model_kinetic_contract` | Zhao 電流と対になる速度空間境界写像の契約 |
| `surface_current_model_photoelectron_active` | PE channel の有無。false なら PE 関連 receipt は 0 |
| `surface_current_model_reference_area_m2` | 電流 target の参照面積 |
| `surface_current_model_phi0_V`, `surface_current_model_phi_m_V` | 解かれた表面電位と極小電位 |
| `surface_current_model_ambient_electron_density_m3` | 解かれた ambient electron 密度 |
| `surface_current_model_electron_current_density_A_m2`, `surface_current_model_ion_current_density_A_m2` | signed electron / ion 電流密度 |
| `surface_current_model_pe_emission_current_density_A_m2`, `surface_current_model_pe_escape_current_density_A_m2`, `surface_current_model_pe_return_current_density_A_m2` | signed PE 電流密度。emission と escape は正、return は負 |
| `surface_current_model_pe_escape_particle_current_A` | 外向き PE 粒子が運ぶ signed 電流なので負 |
| `surface_current_model_net_current_density_A_m2` | 全 channel の表面正味電流密度 |
| `surface_current_model_pe_budget_residual_current_density_A_m2` | PE の emission-return-escape 連続条件の残差 |
| `surface_current_model_surface_budget_residual_current_density_A_m2` | 表面定常電流の残差 |
| `surface_current_model_electron_inflow_reservoir_potential_V`, `surface_current_model_electron_inflow_access_potential_V`, `surface_current_model_electron_inflow_face` | electron inflow の reservoir、access potential、face |
| `surface_current_model_pe_outflow_barrier_potential_V`, `surface_current_model_pe_outflow_barrier_face` | PE outflow barrier の電位と face。PE inactive では 0 |

face 番号は `1..6 = x_low, x_high, y_low, y_high, z_low, z_high` です。

### `matching_plane_quasistatic`

この model では、静的な surface-current 電流 target ではなく、下表に列挙する accepted batch ごとの
matching-plane receipt を状態値として読みます。
固定点式、応答 CSV、`implicit_zero_mode` の正本は
[matching-plane 数値・応答表リファレンス](MatchingPlaneReference.html)です。

| 項目 | `matching_plane_history.csv` の契約 |
| --- | --- |
| 生成条件 | `output.write_files=true`、`output.history_stride>0`、`surface_current_model.model=matching_plane_quasistatic` |
| 1 行 | 先頭の `batch`, `simulated_time_s` と下表の 15 状態列。accepted state だけを書く |
| 解釈 | 同じ行の flux、barrier、反復 receipt を 1 batch の固定点解として読む |
| Python | `result.matching_plane_history`。最後の accepted state は `result.matching_plane_state` |
| 再開 | CSV 自体は使わない。schema v9 は最後の accepted state を `summary.txt` に保存する |

#### 受理済み状態

| `summary.txt` の receipt | `matching_plane_history.csv` の列 | 意味 |
| --- | --- | --- |
| `matching_plane_displacement_C_m2` | `D_H_C_m2` | matching plane の変位電荷密度 $D_H$ |
| `matching_plane_phi_V` | `phi_H_V` | matching plane の電位 $\Phi_H$ |
| `matching_plane_electron_inward_flux_m2_s` | `electron_inward_flux_m2_s` | electron inward flux |
| `matching_plane_ion_inward_flux_m2_s` | `ion_inward_flux_m2_s` | ion inward flux |
| `matching_plane_electron_access_potential_V` | `electron_access_potential_V` | electron access potential |
| `matching_plane_ion_access_potential_V` | `ion_access_potential_V` | ion access potential |
| `matching_plane_photoelectron_barrier_potential_V` | `photoelectron_barrier_potential_V` | PE barrier potential |
| `matching_plane_photoelectron_outward_flux_m2_s` | `photoelectron_outward_flux_m2_s` | 固定点へ返す PE outward flux |
| `matching_plane_photoelectron_mean_normal_energy_eV` | `photoelectron_mean_normal_energy_eV` | 固定点へ返す PE 平均法線 energy |
| `matching_plane_electron_outward_flux_m2_s` | `electron_outward_flux_m2_s` | 固定点へ返す electron outward flux |
| `matching_plane_ion_outward_flux_m2_s` | `ion_outward_flux_m2_s` | 固定点へ返す ion outward flux |
| `matching_plane_photoelectron_return_flux_m2_s` | `photoelectron_return_flux_m2_s` | backend が返す PE return flux |
| `matching_plane_photoelectron_escape_flux_m2_s` | `photoelectron_escape_flux_m2_s` | backend が返す PE escape flux |
| `matching_plane_iterations` | `iterations` | 固定点反復回数 |
| `matching_plane_residual` | `residual` | 受理時の有効相対残差 |

history の先頭 2 列は `batch`, `simulated_time_s` です。したがって全体は 17 列です。
`matching_plane_state_valid` が false の summary にある上表の状態 receipt は、accepted state として
使えません。

PE 収支は、同じ batch について次式で確認します。

$$
\mathtt{photoelectron\_outward\_flux\_m2\_s}
=\mathtt{photoelectron\_return\_flux\_m2\_s}
+\mathtt{photoelectron\_escape\_flux\_m2\_s}.
$$

各行の $D_H$ と $\Phi_H$ は、その batch の粒子追跡に使った commit 前の表面電荷 state に対応します。
一方、`simulated_time_s` は trial を受理して進めた後の時刻です。次 batch 開始時の post-commit 場とは区別します。

#### 実行条件と収束条件

| receipt | 内容 |
| --- | --- |
| `surface_current_model_response_backend` | `table` または `zhao_online` |
| `surface_current_model_matching_plane_z_m` | matching plane の z 座標 |
| `surface_current_model_electron_species`, `surface_current_model_ion_species`, `surface_current_model_photoelectron_species` | 各 channel に割り当てた species |
| `surface_current_model_coupling_rtol` | active 成分の相対許容値 |
| `surface_current_model_coupling_atol` | 4 成分の絶対許容値 |
| `surface_current_model_coupling_max_iterations` | 最大反復回数 |
| `surface_current_model_coupling_relaxation` | 固定点緩和係数 |
| `surface_current_model_dynamic_state_source` | accepted state では `surface_current_model_dynamic_state_source=accepted_batch_fixed_point` |

`surface_current_model_coupling_atol` の 4 値の順序は次のとおりです。

1. PE 外向き flux [m^-2 s^-1]
2. PE 平均法線 energy [eV]
3. electron 外向き flux [m^-2 s^-1]
4. ion 外向き flux [m^-2 s^-1]

既定値はすべて 0 で、inactive 成分も 0 でなければなりません。active 成分の判定閾値は
`max(coupling_rtol * backend_scale, coupling_atol)` です。絶対許容値が支配する成分は有効残差へ換算されるため、
accepted state の `matching_plane_residual` は `surface_current_model_coupling_rtol` 以下になります。

#### Backend 固有の receipt

| backend | receipt |
| --- | --- |
| `table` | `surface_current_model_response_table_path`, `surface_current_model_response_content_fingerprint`。fingerprint は読込済み応答表の canonical 内容を識別する |
| `zhao_online` | `surface_current_model_response_contract=matching_plane_zhao_online_v1` |
| `zhao_online` | `surface_current_model_zhao_branch` |
| `zhao_online` | `surface_current_model_outer_solver=charge_driven_finite_h_sagdeev` |
| `zhao_online` | `surface_current_model_photoelectron_closure=moment_matched_half_maxwellian` |
| `zhao_online` | `surface_current_model_ambient_outward_feedback=transparent` |
| `zhao_online` | `surface_current_model_outer_solver_state=stateless` |

## 履歴

| ファイル | 生成条件 | 列と解釈 | Python |
| --- | --- | --- | --- |
| `charge_history.csv` | `output.write_files=true` かつ `output.history_stride>0` | `batch`, `processed_particles`, `rel_change`, `elem_idx`, `charge_C`。1 snapshot は全要素の電荷 | `result.history` |
| `potential_history.csv` | `output.write_files=true`、`output.write_potential_history=true`、`output.history_stride>0` | `batch`, `elem_idx`, `potential_V`。1 snapshot は全要素重心の電位 | 専用属性なし。CSV を直接読む |
| `top_reference_history.csv` | `output.write_files=true`、`output.write_potential_history=true`、`output.history_stride>0`、`[domain]` の box あり | `batch`, `simulated_time_s`, `z_high_m`, `sample_n`, `potential_mean_V`, `potential_std_V`, `potential_min_V`, `potential_max_V` | 専用属性なし。CSV を直接読む |
| `matching_plane_history.csv` | matching-plane かつ stride 有効 | [17 列の accepted state](#matching_plane_quasistatic) | `result.matching_plane_history` |

`top_reference_history.csv` の基準は box の z-high 面平均であり、無限遠電位や plasma 電位ではありません。
要素相対電位は同じ batch の `potential_history.csv` と結合し、
`potential_V - potential_mean_V` として求めます。

matching-plane では `result.matching_plane_state` と `result.matching_plane_history` を使うと、列番号を
手作業で管理せず typed receipt として参照できます。

## mesh 電位

| 項目 | `mesh_potential.csv` の契約 |
| --- | --- |
| 生成条件 | `output.write_files=true` かつ `output.write_mesh_potential=true` |
| 列 | `elem_idx`, `potential_V` |
| 解釈 | 各要素重心の電位 [V] |
| Python | `result.mesh_potential_v` |
| 再開 | 使用しない |

| field model | 電位へ含める項 |
| --- | --- |
| 共通 | 自己項を解析的 P0 triangle panel kernel で評価 |
| finite periodic2 | 指定された image shell を加算 |
| `cached_kneq0` | cached 非零 mode と設定された物理的 zero mode を合成 |

## charge ledger

| 項目 | `charge_ledger.csv` の契約 |
| --- | --- |
| 生成条件 | `output.write_files=true` かつ charge ledger state あり |
| 1 行 | 1 species の signed charge、補正、weight、粒子数。全 25 列を下表に示す |
| Python | `result.charge_ledger` の `ChargeLedgerEntry` |
| 再開 | `summary.txt` に ledger metadata がある場合は必須 |

列順は次のとおりです。

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

### `summary.txt` の companion receipt

| receipt | 内容 |
| --- | --- |
| `charge_ledger_nspecies`, `charge_ledger_batch_count` | ledger の species 数と batch |
| `charge_ledger_surface_charge_before_C`, `charge_ledger_surface_charge_after_C` | 表面電荷 stock の更新前後 |
| `charge_ledger_local_flight_charge_before_C`, `charge_ledger_local_flight_charge_after_C` | local-flight stock の更新前後 |
| `charge_ledger_unresolved_stock_before_C`, `charge_ledger_unresolved_stock_after_C` | unresolved stock の更新前後 |
| `charge_ledger_residual_C`, `charge_ledger_discarded_unresolved_abs_C` | 保存残差と unresolved 絶対電荷 |
| `charge_ledger_neutral_return_correction_C`, `charge_ledger_fixed_current_correction_C` | closure 補正の合計 |
| `charge_ledger_fixed_absorbed_applied_charge_C`, `charge_ledger_fixed_emission_applied_charge_C` | 表面 fixed-current target の合計 |
| `charge_ledger_raw_escape_charge_C`, `charge_ledger_fixed_escape_target_charge_C`, `charge_ledger_fixed_escape_applied_charge_C`, `charge_ledger_fixed_escape_correction_C` | raw escape、target、互換 alias、補正 |
| `charge_ledger_fixed_applied_surface_net_charge_C` | fixed-current が表面へ適用した正味電荷 |
| `charge_ledger_fixed_pe_continuity_residual_C` | PE active 時の fixed PE 連続条件残差 |

### 判定規則

| 確認対象 | 一緒に読む値 | 判定上の注意 |
| --- | --- | --- |
| 外部 escape | `escaped_to_infinity_C`, `fixed_escape_target_charge_C`, `fixed_escape_correction_C` | raw trajectory 値は上書きされない。escape 補正は表面要素へ deposit されない |
| 保存残差 | `charge_ledger_residual_C` と上表の surface / local-flight / unresolved の before / after | stock、外部 flux、neutral-return、fixed-current の補正を含む |
| `fixed_current` の標本分散 | `absorbed_count`, `emitted_count`, raw charge, `fixed_absorbed_weight_scale`, `fixed_emission_weight_scale` | count が 1 なら target 全量が 1 標本の要素へ局所化される |
| closed PE | raw absorption、unresolved 値、`neutral_return_correction_C`, `neutral_return_unresolved_fraction` | 総電荷が閉じても未解決率は別に確認する |

既存 reader 向けの出力 alias と、内部 ledger が保持する同値の target は次の対応です。

| 出力 alias | 内部 target |
| --- | --- |
| `fixed_absorbed_applied_charge_C` | `fixed_absorbed_target_charge_C` |
| `fixed_emission_applied_charge_C` | `fixed_emission_target_charge_C` |
| `fixed_escape_applied_charge_C` | `fixed_escape_target_charge_C` |

小さい保存残差は要素別の標本収束を保証しません。粒子数、ray 数、batch 幅、乱数 seed を変えた電荷分布の
収束を別に確認してください。基本的な照合手順は
[charge ledger を確認する](OutputGuide.html#charge-ledger)にあります。

## 適応 batch の診断

`periodic2.max_nonzero_mode_potential_step > 0` の場合、次の receipt を `summary.txt` に記録します。

| receipt | 内容 |
| --- | --- |
| `simulated_time_s` | accepted trial までの模擬時間 |
| `periodic2_max_nonzero_mode_potential_step_V` | 設定した非零 mode 電位変化上限 |
| `adaptive_nonzero_mode_rejected_trials` | 棄却した trial の累積数 |
| `adaptive_nonzero_mode_last_batch_duration_s` | 最終 accepted batch の duration |
| `adaptive_nonzero_mode_last_potential_step_V` | 最終 accepted batch の非零 mode 電位変化 |
| `adaptive_nonzero_mode_omp_threads` | adaptive restart で固定する OpenMP thread 数 |

棄却 trial は履歴や ledger へ出力されません。

## 再開に使うファイル

`output.checkpoint_stride > 0` では、accepted batch の commit 後に次の構造を更新します。
`output.checkpoint_stride=0` でも、正常終了時の最終 checkpoint は出力されます。

この節の「完全」は**再開用ファイル一式の transaction が完全**という意味です。要求した実行の完了は
`summary.txt` の `batches == sim.batch_count` で別に判定します。

```text
<output.dir>/
├── checkpoint_complete.txt
├── checkpoint_latest.txt
└── checkpoints/
    ├── slot0/
    └── slot1/
```

### 各 checkpoint の内容

| ファイル | 必須条件 | 形式と役割 |
| --- | --- | --- |
| `summary.txt` | 常に必須 | `key=value` の統計、schema、fingerprint、ledger stock |
| `charges.csv` | 常に必須 | `elem_idx`, `charge_C` の要素電荷 |
| `rng_state.txt` | serial で必須 | RNG 配列長と state |
| `rng_state_rankNNNNN.txt` | MPI で全 rank 分が必須 | rank ごとの RNG 配列長と state |
| `macro_residuals.csv` | schema v8 以降で manifest が宣言した場合 | `species_idx`, `face`, `residual`。global な macro 粒子数残差 |
| `charge_ledger.csv` | summary / manifest に ledger metadata がある場合 | [25 列の species ledger](#charge-ledger) |
| `checkpoint_complete.txt` | schema v8 以降で必須 | `schema_version`, `state`, `batches`, `mpi_world_size`, `macro_residuals_present`, `charge_ledger_present` |

### 更新と自動選択

| 段階 | 契約 |
| --- | --- |
| 更新 | 非 active slot を書き終えてから `checkpoint_latest.txt` を原子的に切り替える。保持は最大 2 世代 |
| 探索 | `output.restart_from` 直下の最終出力と `slot0` / `slot1` を検査する |
| 選択 | 必須ファイルが揃う load 可能な候補のうち、`batches` が最大のものを選ぶ |
| index 障害 | `checkpoint_latest.txt` が欠落、破損、古い場合も、完了 manifest を持つ slot を回収する |
| 不整合 | 必須ファイル、fingerprint、mesh 要素数、species 数、MPI world size が不一致なら、新規実行へ fallback せず停止する |

`output.restart_from` は読み込み元だけを変更し、新しい出力は `output.dir` に書きます。

### Schema 世代

| schema | 契約 |
| --- | --- |
| v6 | `macro_residuals.csv` は `species_idx,face,residual`。`face=0` は従来 source、`1..6` は boundary face。旧 2 列形式 `species_idx,residual` も読込可能 |
| v8 以降 | 書き始めに `state=in_progress` を公開する。summary、charges、全 rank の RNG、manifest 宣言済み residual / ledger を閉じてから `state=complete` を原子的に公開する |
| v9 | matching-plane の accepted feedback、potential、return / escape flux、反復 receipt を `summary.txt` に保存する |

schema v8 以降では `checkpoint_complete.txt` 自体が必須です。直下の最終出力に新旧世代のファイルが混在した
場合、そのディレクトリを選ばず、完全な periodic slot へ戻ります。

model fingerprint は table backend なら response table の canonical 内容、online backend なら Zhao contract と
branch policy を含みます。これらを変更した checkpoint は再開できません。

実際の操作は[実行・再開する](Execution.html)を参照してください。

## Python から読む

各ファイルの属性または専用属性の有無は[ファイル生成条件](#ファイル生成条件)の同じ行に示しています。
代表的な typed 属性は次のように読みます。

```python
from beach import load_fortran_result

result = load_fortran_result("outputs/tutorial")
print(result.charges)
print(result.charge_ledger)
print(result.field_reconstruction)
print(result.matching_plane_state)
```

| 出力 | Python 属性 |
| --- | --- |
| `charges.csv` | `result.charges` |
| `charge_ledger.csv` | `result.charge_ledger` |
| field reconstruction receipt | `result.field_reconstruction` |
| matching-plane summary | `result.matching_plane_state` |
| `matching_plane_history.csv` | `result.matching_plane_history` |

`summary.txt` と `charges.csv` は必須です。mesh、履歴、ledger、field reconstruction、matching-plane state は
対応する出力がある場合だけ読み込まれます。

`matching_plane_state_valid` が false なら、reader は
`matching_plane_state` と `matching_plane_history` を accepted state として公開しません。

別のケースでは `outputs/tutorial` を実際の `output.dir` に置き換えます。class と属性の詳細は
[Python 後処理 API](PythonPostprocessAPI.html)、物理・数値的な受理判断は
[計算結果の妥当性確認](ValidationGuide.html)を参照してください。

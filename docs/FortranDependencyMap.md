title: Fortran 依存関係マップ

Lang: [日本語](FortranDependencyMap.md) | [English](FortranDependencyMap.en.md)

# Fortran 依存関係マップ

> このページは `tools/generate_fortran_dependency_report.py` から自動生成しています。

## 概要

- ソースファイル数: 71
- モジュール数: 57
- submodule 数: 13
- program 数: 1
- 内部依存エッジ数: 266

## 全体グラフ

![Fortran モジュール依存グラフ](../media/fortran_module_dependencies.svg)

実線は `use` 依存、破線は `submodule(parent)` の親参照を表します。

## ディレクトリ別サマリ

| ディレクトリ | エンティティ数 | 内部依存数 |
| --- | ---: | ---: |
| `app` | 1 | 12 |
| `src/config` | 6 | 33 |
| `src/config/app_config_parser` | 3 | 10 |
| `src/core` | 7 | 5 |
| `src/mesh` | 3 | 11 |
| `src/particles` | 2 | 9 |
| `src/physics` | 7 | 24 |
| `src/physics/field_solver` | 5 | 20 |
| `src/physics/field_solver/fmm/api` | 4 | 8 |
| `src/physics/field_solver/fmm/internal/common` | 2 | 5 |
| `src/physics/field_solver/fmm/internal/periodic` | 6 | 24 |
| `src/physics/field_solver/fmm/internal/runtime` | 2 | 14 |
| `src/physics/field_solver/fmm/internal/tree` | 2 | 12 |
| `src/physics/panel` | 5 | 13 |
| `src/physics/periodic_zero_mode` | 3 | 9 |
| `src/runtime` | 6 | 26 |
| `src/runtime/coupling` | 1 | 1 |
| `src/runtime/simulator` | 6 | 30 |

## 被依存の多いモジュール

| エンティティ | kind | 被依存数 |
| --- | --- | ---: |
| `bem_kinds` | `module` | 53 |
| `bem_types` | `module` | 23 |
| `bem_string_utils` | `module` | 17 |
| `bem_constants` | `module` | 13 |
| `bem_panel_geometry` | `module` | 12 |
| `bem_coulomb_fmm_types` | `module` | 9 |
| `bem_app_config_types` | `module` | 8 |
| `bem_coulomb_fmm_core` | `module` | 8 |
| `bem_mpi` | `module` | 7 |
| `bem_periodic_zero_mode_plan` | `module` | 7 |

## エンティティ一覧

| エンティティ | kind | パス | 内部依存 | 概要 |
| --- | --- | --- | --- | --- |
| `main` | `program` | `app/main.f90` | `bem_kinds`, `bem_version`, `bem_types`, `bem_mpi`, `bem_performance_profile`, `bem_simulator`, `bem_restart`, `bem_output_writer`, `bem_app_config`, `bem_mesh`, `bem_charge_ledger`, `bem_electrostatic_snapshot` | 設定読込・メッシュ生成・粒子初期化・シミュレーション実行・結果出力を順に行うCLIエントリーポイント。 |
| `bem_app_config` | `module` | `src/config/bem_app_config.f90` | `bem_app_config_types`, `bem_physics_config_types`, `bem_app_config_parser`, `bem_string_utils`, `bem_app_config_runtime` | 設定型・TOMLパーサ・実行時変換ロジックを束ねる後方互換ファサード。 |
| `bem_app_config_authoring` | `module` | `src/config/bem_app_config_authoring.f90` | `bem_kinds`, `bem_app_config_types`, `bem_string_utils` | BEACH TOML の高水準 authoring キーを実行時設定へ正規化する補助モジュール。 |
| `bem_app_config_runtime` | `module` | `src/config/bem_app_config_runtime.f90` | `bem_kinds`, `bem_constants`, `bem_types`, `bem_mpi`, `bem_electrostatic_snapshot`, `bem_templates`, `bem_mesh`, `bem_panel_surface_sides`, `bem_collision`, `bem_importers`, `bem_injection`, `bem_particles`, `bem_external_boundary_contract`, `bem_app_config_types`, `bem_string_utils`, `bem_config_helpers` | `app_config` からメッシュ・粒子群を構築する実行時変換モジュール。 |
| `bem_app_config_types` | `module` | `src/config/bem_app_config_types.f90` | `bem_kinds`, `bem_types`, `bem_physics_config_types` | アプリ設定の型定義と、設定由来の粒子数計算をまとめるモジュール。 |
| `bem_config_helpers` | `module` | `src/config/bem_config_helpers.f90` | `bem_kinds`, `bem_app_config_types`, `bem_string_utils` | 設定型のヘルパー関数。パーサに依存せず下位層から利用可能。 |
| `bem_physics_config_types` | `module` | `src/config/bem_physics_config_types.f90` | `bem_kinds`, `bem_types`, `bem_string_utils` | 場・periodic2・panel の型付き設定と旧 `[sim]` 設定の正規化を定義する。 |
| `bem_app_config_parser` | `module` | `src/config/app_config_parser/bem_app_config_parser.f90` | `bem_kinds`, `bem_constants`, `bem_types`, `bem_app_config_types`, `bem_physics_config_types`, `bem_app_config_authoring`, `bem_string_utils` | TOML設定ファイルを `toml-f` で読み込み、`app_config` へ反映する。 |
| `bem_app_config_parser_finalize` | `submodule` | `src/config/app_config_parser/bem_app_config_parser_finalize.f90` | `bem_app_config_parser` | 読み込み済み TOML 設定の正規化・派生値確定・検証を実装する submodule。 |
| `bem_app_config_parser_validate` | `submodule` | `src/config/app_config_parser/bem_app_config_parser_validate.f90` | `bem_app_config_parser`, `bem_config_helpers` | `bem_app_config_parser` の入力検証・物理量導出手続きを実装する submodule。 |
| `bem_constants` | `module` | `src/core/bem_constants.f90` | `bem_kinds` | シミュレーションで使用する物理定数を定義する。 |
| `bem_external_boundary_contract` | `module` | `src/core/bem_external_boundary_contract.f90` | `bem_kinds`, `bem_string_utils` | 局所 reservoir 流入と通常 open 面の処理を実行時契約へ正規化する。 |
| `bem_kinds` | `module` | `src/core/bem_kinds.f90` | - | 倍精度実数と32bit整数のkind定義を集約する基盤モジュール。 |
| `bem_mpi` | `module` | `src/core/bem_mpi.F90` | `bem_kinds` | MPIの初期化・集約を抽象化し、非MPIビルドでは単一ランク動作へフォールバックする。 |
| `bem_string_utils` | `module` | `src/core/bem_string_utils.f90` | - | ASCII 文字列操作ユーティリティ。 |
| `bem_types` | `module` | `src/core/bem_types.f90` | `bem_kinds` | シミュレーション設定・統計・メッシュ・粒子・衝突情報の主要データ型を定義する。 |
| `bem_version` | `module` | `src/core/bem_version.F90` | - | BEACH build metadata stamped by build.sh. |
| `bem_importers` | `module` | `src/mesh/bem_importers.f90` | `bem_kinds`, `bem_types`, `bem_mesh` | OBJメッシュを走査・解析し、内部 `mesh_type` へ変換するインポートモジュール。 |
| `bem_mesh` | `module` | `src/mesh/bem_mesh.f90` | `bem_kinds`, `bem_types`, `bem_string_utils`, `bem_panel_geometry`, `bem_panel_quadrature` | 三角形メッシュ幾何量(重心・法線・AABB・代表長)を前計算して保持するモジュール。 |
| `bem_templates` | `module` | `src/mesh/bem_templates.f90` | `bem_kinds`, `bem_types`, `bem_mesh` | 平面/穴あき平面/円板/リング/箱/円柱/球テンプレートから三角形メッシュを生成するユーティリティ。 |
| `bem_injection` | `module` | `src/particles/bem_injection.f90` | `bem_kinds`, `bem_constants`, `bem_particles`, `bem_types`, `bem_boundary`, `bem_collision`, `bem_string_utils` | 乱数シード設定と粒子位置/速度サンプリングを担う粒子注入モジュール。 |
| `bem_particles` | `module` | `src/particles/bem_particles.f90` | `bem_kinds`, `bem_types` | 粒子SoAデータ構造の初期化を提供するモジュール。 |
| `bem_boundary` | `module` | `src/physics/bem_boundary.f90` | `bem_kinds`, `bem_types` | シミュレーションボックス境界（流出/反射/周期）を適用するモジュール。 |
| `bem_collision` | `module` | `src/physics/bem_collision.f90` | `bem_kinds`, `bem_types`, `bem_string_utils` | 粒子軌道セグメントと三角形要素の交差判定を提供する衝突検出モジュール。 |
| `bem_electrostatic_snapshot` | `module` | `src/physics/bem_electrostatic_snapshot.f90` | `bem_kinds`, `bem_types`, `bem_field_solver`, `bem_physics_config_types`, `bem_string_utils`, `bem_periodic_zero_mode_plan`, `bem_periodic_zero_mode_eval`, `bem_coulomb_fmm_periodic_nonzero_reference`, `bem_panel_geometry`, `bem_panel_kernel` | 1バッチ内で固定する静電場 snapshot。 |
| `bem_pusher` | `module` | `src/physics/bem_pusher.f90` | `bem_kinds` | 荷電粒子の時間発展にBoris法を適用する運動方程式ソルバ。 |
| `bem_surface_models` | `module` | `src/physics/bem_surface_models.f90` | `bem_kinds`, `bem_types`, `bem_string_utils` | 表面モデルごとの電荷更新後処理を扱うモジュール。 |
| `bem_electrostatic_snapshot_eval` | `submodule` | `src/physics/bem_electrostatic_snapshot_eval.f90` | `bem_electrostatic_snapshot` | `bem_electrostatic_snapshot` の局所電場・電位評価を実装する submodule。 |
| `bem_surface_models_conductor` | `submodule` | `src/physics/bem_surface_models_conductor.f90` | `bem_surface_models`, `bem_constants`, `bem_panel_geometry`, `bem_panel_kernel` | `bem_surface_models` の浮遊導体電荷再配分を実装する submodule。 |
| `bem_field_kernel_c` | `module` | `src/physics/field_solver/bem_field_kernel_c.f90` | `bem_constants`, `bem_coulomb_fmm_core`, `bem_kinds`, `bem_panel_geometry`, `bem_version` | C ABI wrapper for the simulator-independent Coulomb FMM field kernel. |
| `bem_field_solver` | `module` | `src/physics/field_solver/bem_field_solver.f90` | `bem_kinds`, `bem_constants`, `bem_types`, `bem_coulomb_fmm_core`, `bem_string_utils`, `bem_physics_config_types` | 粒子位置での電場評価を direct / treecode / fmm で切り替える場ソルバ。 |
| `bem_field_solver_config` | `submodule` | `src/physics/field_solver/bem_field_solver_config.f90` | `bem_field_solver`, `bem_coulomb_fmm_core`, `bem_physics_config_types` | `bem_field_solver` の初期化・設定補助手続きを実装する submodule。 |
| `bem_field_solver_eval` | `submodule` | `src/physics/field_solver/bem_field_solver_eval.f90` | `bem_field_solver`, `bem_coulomb_fmm_core`, `bem_panel_geometry`, `bem_panel_kernel` | `bem_field_solver` の電場評価と木走査ロジックを実装する submodule。 |
| `bem_field_solver_tree` | `submodule` | `src/physics/field_solver/bem_field_solver_tree.f90` | `bem_field_solver`, `bem_coulomb_fmm_core` | `bem_field_solver` の octree 構築・更新とメモリ管理を実装する submodule。 |
| `bem_coulomb_fmm_core` | `module` | `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90` | `bem_kinds`, `bem_coulomb_fmm_types` | `mesh_type` や `sim_config` に依存しない Coulomb FMM コア API。 |
| `bem_coulomb_fmm_core_build` | `submodule` | `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_build.f90` | `bem_coulomb_fmm_core`, `bem_coulomb_fmm_plan_ops` | `bem_coulomb_fmm_core` の plan 構築 API ラッパ。 |
| `bem_coulomb_fmm_core_eval` | `submodule` | `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_eval.f90` | `bem_coulomb_fmm_core`, `bem_coulomb_fmm_eval_ops` | `bem_coulomb_fmm_core` の評価 API ラッパ。 |
| `bem_coulomb_fmm_core_state` | `submodule` | `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_state.f90` | `bem_coulomb_fmm_core`, `bem_coulomb_fmm_state_ops` | `bem_coulomb_fmm_core` の state 更新 API ラッパ。 |
| `bem_coulomb_fmm_basis` | `module` | `src/physics/field_solver/fmm/internal/common/bem_coulomb_fmm_basis.f90` | `bem_kinds`, `bem_coulomb_fmm_types` | Coulomb FMM の multi-index と微分テーブル計算。 |
| `bem_coulomb_fmm_types` | `module` | `src/physics/field_solver/fmm/internal/common/bem_coulomb_fmm_types.f90` | `bem_kinds`, `bem_panel_geometry`, `bem_periodic_zero_mode_plan` | Coulomb FMM コアで共有する型定義。 |
| `bem_coulomb_fmm_periodic` | `module` | `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic.f90` | `bem_kinds`, `bem_coulomb_fmm_types` | Coulomb FMM の periodic2 境界処理。 |
| `bem_coulomb_fmm_periodic_cache` | `module` | `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_cache.f90` | `bem_filesystem`, `bem_kinds` | Versioned stream codec for cached periodic root operators. |
| `bem_coulomb_fmm_periodic_ewald` | `module` | `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_ewald.f90` | `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_periodic` | periodic2 cached operator生成用のEwald teacher。 |
| `bem_coulomb_fmm_periodic_nonzero_reference` | `module` | `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_nonzero_reference.f90` | `bem_kinds`, `bem_constants`, `bem_types`, `bem_panel_geometry`, `bem_panel_quadrature` | - |
| `bem_coulomb_fmm_periodic_root_ops` | `module` | `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90` | `bem_version`, `bem_kinds`, `bem_mpi`, `bem_filesystem`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_basis`, `bem_coulomb_fmm_periodic`, `bem_coulomb_fmm_periodic_ewald`, `bem_coulomb_fmm_periodic_cache`, `bem_regularized_qr`, `bem_coulomb_fmm_tree_utils` | periodic2 root operator の前計算。 |
| `bem_regularized_qr` | `module` | `src/physics/field_solver/fmm/internal/periodic/bem_regularized_qr.f90` | `bem_kinds` | Reusable column-scaled QR factorization for ridge-regularized least squares. |
| `bem_coulomb_fmm_eval_ops` | `module` | `src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90` | `bem_kinds`, `bem_constants`, `bem_panel_geometry`, `bem_panel_kernel`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_basis`, `bem_coulomb_fmm_periodic`, `bem_coulomb_fmm_tree_utils`, `bem_periodic_zero_mode_eval` | Coulomb FMM 電場評価。 |
| `bem_coulomb_fmm_state_ops` | `module` | `src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_state_ops.f90` | `bem_kinds`, `bem_constants`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_tree_utils`, `bem_periodic_zero_mode_plan` | Coulomb FMM state 更新と upward/downward pass。 |
| `bem_coulomb_fmm_plan_ops` | `module` | `src/physics/field_solver/fmm/internal/tree/bem_coulomb_fmm_plan_ops.f90` | `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_basis`, `bem_coulomb_fmm_periodic`, `bem_coulomb_fmm_periodic_ewald`, `bem_coulomb_fmm_periodic_root_ops`, `bem_panel_geometry`, `bem_periodic_zero_mode_plan`, `bem_coulomb_fmm_tree_utils` | Coulomb FMM plan 構築と tree トポロジ前計算。 |
| `bem_coulomb_fmm_tree_utils` | `module` | `src/physics/field_solver/fmm/internal/tree/bem_coulomb_fmm_tree_utils.f90` | `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_periodic` | Coulomb FMM tree 構造の共通ユーティリティ。 |
| `bem_panel_geometry` | `module` | `src/physics/panel/bem_panel_geometry.f90` | `bem_kinds` | Ordered triangle geometry and exact Cartesian surface moments. |
| `bem_panel_kernel` | `module` | `src/physics/panel/bem_panel_kernel.f90` | `bem_kinds`, `bem_constants`, `bem_panel_geometry`, `bem_panel_self_terms` | Analytic free-space P0 triangle potential, field, principal value, and jump. |
| `bem_panel_quadrature` | `module` | `src/physics/panel/bem_panel_quadrature.f90` | `bem_kinds`, `bem_constants`, `bem_panel_geometry` | Independent triangle cubature and Gauss-Duffy correctness oracles. |
| `bem_panel_self_terms` | `module` | `src/physics/panel/bem_panel_self_terms.f90` | `bem_kinds`, `bem_panel_geometry` | On-surface P0 triangle potential and principal-value field integrals. |
| `bem_panel_surface_sides` | `module` | `src/physics/panel/bem_panel_surface_sides.f90` | `bem_kinds`, `bem_types`, `bem_string_utils` | Resolve physical vacuum sides without changing ordered triangle winding. |
| `bem_periodic_zero_mode_c` | `module` | `src/physics/periodic_zero_mode/bem_periodic_zero_mode_c.f90` | `bem_kinds`, `bem_periodic_zero_mode_eval`, `bem_periodic_zero_mode_plan` | C ABI wrapper for the physical periodic zero mode. |
| `bem_periodic_zero_mode_eval` | `module` | `src/physics/periodic_zero_mode/bem_periodic_zero_mode_eval.f90` | `bem_kinds`, `bem_constants`, `bem_periodic_zero_mode_plan` | - |
| `bem_periodic_zero_mode_plan` | `module` | `src/physics/periodic_zero_mode/bem_periodic_zero_mode_plan.f90` | `bem_kinds`, `bem_constants`, `bem_types` | - |
| `bem_checkpoint_contract` | `module` | `src/runtime/bem_checkpoint_contract.f90` | `bem_kinds` | BEACH checkpoint metadata shared by output and restart code. |
| `bem_filesystem` | `module` | `src/runtime/bem_filesystem.f90` | - | Minimal POSIX filesystem operations used by the runtime. |
| `bem_model_fingerprint` | `module` | `src/runtime/bem_model_fingerprint.f90` | `bem_kinds`, `bem_types`, `bem_app_config_types` | Restart compatibility fingerprints for the ordered physical model contract. |
| `bem_output_writer` | `module` | `src/runtime/bem_output_writer.f90` | `bem_kinds`, `bem_types`, `bem_app_config_types`, `bem_charge_ledger`, `bem_checkpoint_contract`, `bem_electrostatic_snapshot`, `bem_external_boundary_contract`, `bem_model_fingerprint`, `bem_version`, `bem_filesystem`, `bem_string_utils` | 実行サマリ・最終CSV・履歴CSVの出力を担当するモジュール。 |
| `bem_performance_profile` | `module` | `src/runtime/bem_performance_profile.f90` | `bem_kinds`, `bem_mpi`, `bem_string_utils` | 実行フェーズごとの壁時計計測と MPI 集約出力を担う軽量プロファイラ。 |
| `bem_restart` | `module` | `src/runtime/bem_restart.f90` | `bem_kinds`, `bem_types`, `bem_app_config_types`, `bem_charge_ledger`, `bem_model_fingerprint`, `bem_physics_config_types`, `bem_checkpoint_contract`, `bem_mpi` | チェックポイントファイルの保存/復元を扱う補助モジュール。 |
| `bem_charge_ledger` | `module` | `src/runtime/coupling/bem_charge_ledger.f90` | `bem_kinds` | batch 間の signed charge stock と移送 flux から電荷収支を集計する。 |
| `bem_particle_stepper` | `module` | `src/runtime/simulator/bem_particle_stepper.f90` | `bem_kinds`, `bem_types`, `bem_electrostatic_snapshot`, `bem_pusher`, `bem_collision`, `bem_boundary`, `bem_external_boundary_contract` | 同一時刻の粒子状態から、空間電場を中点評価した1ステップ候補を構築する。 |
| `bem_simulator` | `module` | `src/runtime/simulator/bem_simulator.f90` | `bem_kinds`, `bem_types`, `bem_app_config`, `bem_app_config_runtime`, `bem_electrostatic_snapshot`, `bem_particle_stepper`, `bem_collision`, `bem_surface_models`, `bem_charge_ledger`, `bem_string_utils`, `bem_external_boundary_contract`, `bem_simulator_workspace`, `bem_mpi` | 吸着(insulator)モデルのメインループを実行し、電荷堆積と統計更新を行う。 |
| `bem_simulator_workspace` | `module` | `src/runtime/simulator/bem_simulator_workspace.f90` | `bem_kinds` | シミュレーション実行中に再利用するバッチ作業配列を管理する。 |
| `bem_simulator_io` | `submodule` | `src/runtime/simulator/bem_simulator_io.f90` | `bem_simulator`, `bem_app_config_runtime`, `bem_output_writer` | `bem_simulator` の進捗表示と履歴出力を実装する submodule。 |
| `bem_simulator_loop` | `submodule` | `src/runtime/simulator/bem_simulator_loop.f90` | `bem_simulator`, `bem_app_config_runtime`, `bem_periodic_zero_mode_plan`, `bem_performance_profile`, `bem_mpi` | `bem_simulator` の主ループと粒子処理計算を実装する submodule。 |
| `bem_simulator_stats` | `submodule` | `src/runtime/simulator/bem_simulator_stats.f90` | `bem_simulator` | `bem_simulator` のバッチ集計・統計更新処理を実装する submodule。 |

## 詳細

### `main`

- kind: `program`
- path: `app/main.f90`
- group: `app`
- 内部依存: `bem_kinds`, `bem_version`, `bem_types`, `bem_mpi`, `bem_performance_profile`, `bem_simulator`, `bem_restart`, `bem_output_writer`, `bem_app_config`, `bem_mesh`, `bem_charge_ledger`, `bem_electrostatic_snapshot`
- external dependencies: `iso_fortran_env`
- 概要: 設定読込・メッシュ生成・粒子初期化・シミュレーション実行・結果出力を順に行うCLIエントリーポイント。

### `bem_app_config`

- kind: `module`
- path: `src/config/bem_app_config.f90`
- group: `src/config`
- 内部依存: `bem_app_config_types`, `bem_physics_config_types`, `bem_app_config_parser`, `bem_string_utils`, `bem_app_config_runtime`
- external dependencies: なし
- 概要: 設定型・TOMLパーサ・実行時変換ロジックを束ねる後方互換ファサード。

### `bem_app_config_authoring`

- kind: `module`
- path: `src/config/bem_app_config_authoring.f90`
- group: `src/config`
- 内部依存: `bem_kinds`, `bem_app_config_types`, `bem_string_utils`
- external dependencies: `ieee_arithmetic`
- 概要: BEACH TOML の高水準 authoring キーを実行時設定へ正規化する補助モジュール。

### `bem_app_config_runtime`

- kind: `module`
- path: `src/config/bem_app_config_runtime.f90`
- group: `src/config`
- 内部依存: `bem_kinds`, `bem_constants`, `bem_types`, `bem_mpi`, `bem_electrostatic_snapshot`, `bem_templates`, `bem_mesh`, `bem_panel_surface_sides`, `bem_collision`, `bem_importers`, `bem_injection`, `bem_particles`, `bem_external_boundary_contract`, `bem_app_config_types`, `bem_string_utils`, `bem_config_helpers`
- external dependencies: `iso_fortran_env`, `ieee_arithmetic`
- 概要: `app_config` からメッシュ・粒子群を構築する実行時変換モジュール。

### `bem_app_config_types`

- kind: `module`
- path: `src/config/bem_app_config_types.f90`
- group: `src/config`
- 内部依存: `bem_kinds`, `bem_types`, `bem_physics_config_types`
- external dependencies: なし
- 概要: アプリ設定の型定義と、設定由来の粒子数計算をまとめるモジュール。

### `bem_config_helpers`

- kind: `module`
- path: `src/config/bem_config_helpers.f90`
- group: `src/config`
- 内部依存: `bem_kinds`, `bem_app_config_types`, `bem_string_utils`
- external dependencies: なし
- 概要: 設定型のヘルパー関数。パーサに依存せず下位層から利用可能。

### `bem_physics_config_types`

- kind: `module`
- path: `src/config/bem_physics_config_types.f90`
- group: `src/config`
- 内部依存: `bem_kinds`, `bem_types`, `bem_string_utils`
- external dependencies: `ieee_arithmetic`
- 概要: 場・periodic2・panel の型付き設定と旧 `[sim]` 設定の正規化を定義する。

### `bem_app_config_parser`

- kind: `module`
- path: `src/config/app_config_parser/bem_app_config_parser.f90`
- group: `src/config/app_config_parser`
- 内部依存: `bem_kinds`, `bem_constants`, `bem_types`, `bem_app_config_types`, `bem_physics_config_types`, `bem_app_config_authoring`, `bem_string_utils`
- external dependencies: `tomlf`, `ieee_arithmetic`
- 概要: TOML設定ファイルを `toml-f` で読み込み、`app_config` へ反映する。

### `bem_app_config_parser_finalize`

- kind: `submodule`
- path: `src/config/app_config_parser/bem_app_config_parser_finalize.f90`
- group: `src/config/app_config_parser`
- parent: `bem_app_config_parser`
- 内部依存: `bem_app_config_parser`
- external dependencies: なし
- 概要: 読み込み済み TOML 設定の正規化・派生値確定・検証を実装する submodule。

### `bem_app_config_parser_validate`

- kind: `submodule`
- path: `src/config/app_config_parser/bem_app_config_parser_validate.f90`
- group: `src/config/app_config_parser`
- parent: `bem_app_config_parser`
- 内部依存: `bem_app_config_parser`, `bem_config_helpers`
- external dependencies: なし
- 概要: `bem_app_config_parser` の入力検証・物理量導出手続きを実装する submodule。

### `bem_constants`

- kind: `module`
- path: `src/core/bem_constants.f90`
- group: `src/core`
- 内部依存: `bem_kinds`
- external dependencies: なし
- 概要: シミュレーションで使用する物理定数を定義する。

### `bem_external_boundary_contract`

- kind: `module`
- path: `src/core/bem_external_boundary_contract.f90`
- group: `src/core`
- 内部依存: `bem_kinds`, `bem_string_utils`
- external dependencies: なし
- 概要: 局所 reservoir 流入と通常 open 面の処理を実行時契約へ正規化する。

### `bem_kinds`

- kind: `module`
- path: `src/core/bem_kinds.f90`
- group: `src/core`
- 内部依存: なし
- external dependencies: `iso_fortran_env`
- 概要: 倍精度実数と32bit整数のkind定義を集約する基盤モジュール。

### `bem_mpi`

- kind: `module`
- path: `src/core/bem_mpi.F90`
- group: `src/core`
- 内部依存: `bem_kinds`
- external dependencies: なし
- 概要: MPIの初期化・集約を抽象化し、非MPIビルドでは単一ランク動作へフォールバックする。

### `bem_string_utils`

- kind: `module`
- path: `src/core/bem_string_utils.f90`
- group: `src/core`
- 内部依存: なし
- external dependencies: なし
- 概要: ASCII 文字列操作ユーティリティ。

### `bem_types`

- kind: `module`
- path: `src/core/bem_types.f90`
- group: `src/core`
- 内部依存: `bem_kinds`
- external dependencies: なし
- 概要: シミュレーション設定・統計・メッシュ・粒子・衝突情報の主要データ型を定義する。

### `bem_version`

- kind: `module`
- path: `src/core/bem_version.F90`
- group: `src/core`
- 内部依存: なし
- external dependencies: なし
- 概要: BEACH build metadata stamped by build.sh.

### `bem_importers`

- kind: `module`
- path: `src/mesh/bem_importers.f90`
- group: `src/mesh`
- 内部依存: `bem_kinds`, `bem_types`, `bem_mesh`
- external dependencies: なし
- 概要: OBJメッシュを走査・解析し、内部 `mesh_type` へ変換するインポートモジュール。

### `bem_mesh`

- kind: `module`
- path: `src/mesh/bem_mesh.f90`
- group: `src/mesh`
- 内部依存: `bem_kinds`, `bem_types`, `bem_string_utils`, `bem_panel_geometry`, `bem_panel_quadrature`
- external dependencies: `ieee_arithmetic`
- 概要: 三角形メッシュ幾何量(重心・法線・AABB・代表長)を前計算して保持するモジュール。

### `bem_templates`

- kind: `module`
- path: `src/mesh/bem_templates.f90`
- group: `src/mesh`
- 内部依存: `bem_kinds`, `bem_types`, `bem_mesh`
- external dependencies: なし
- 概要: 平面/穴あき平面/円板/リング/箱/円柱/球テンプレートから三角形メッシュを生成するユーティリティ。

### `bem_injection`

- kind: `module`
- path: `src/particles/bem_injection.f90`
- group: `src/particles`
- 内部依存: `bem_kinds`, `bem_constants`, `bem_particles`, `bem_types`, `bem_boundary`, `bem_collision`, `bem_string_utils`
- external dependencies: `iso_fortran_env`, `ieee_arithmetic`
- 概要: 乱数シード設定と粒子位置/速度サンプリングを担う粒子注入モジュール。

### `bem_particles`

- kind: `module`
- path: `src/particles/bem_particles.f90`
- group: `src/particles`
- 内部依存: `bem_kinds`, `bem_types`
- external dependencies: なし
- 概要: 粒子SoAデータ構造の初期化を提供するモジュール。

### `bem_boundary`

- kind: `module`
- path: `src/physics/bem_boundary.f90`
- group: `src/physics`
- 内部依存: `bem_kinds`, `bem_types`
- external dependencies: `ieee_arithmetic`
- 概要: シミュレーションボックス境界（流出/反射/周期）を適用するモジュール。

### `bem_collision`

- kind: `module`
- path: `src/physics/bem_collision.f90`
- group: `src/physics`
- 内部依存: `bem_kinds`, `bem_types`, `bem_string_utils`
- external dependencies: `ieee_arithmetic`, `iso_fortran_env`
- 概要: 粒子軌道セグメントと三角形要素の交差判定を提供する衝突検出モジュール。

### `bem_electrostatic_snapshot`

- kind: `module`
- path: `src/physics/bem_electrostatic_snapshot.f90`
- group: `src/physics`
- 内部依存: `bem_kinds`, `bem_types`, `bem_field_solver`, `bem_physics_config_types`, `bem_string_utils`, `bem_periodic_zero_mode_plan`, `bem_periodic_zero_mode_eval`, `bem_coulomb_fmm_periodic_nonzero_reference`, `bem_panel_geometry`, `bem_panel_kernel`
- external dependencies: `ieee_arithmetic`
- 概要: 1バッチ内で固定する静電場 snapshot。

### `bem_pusher`

- kind: `module`
- path: `src/physics/bem_pusher.f90`
- group: `src/physics`
- 内部依存: `bem_kinds`
- external dependencies: なし
- 概要: 荷電粒子の時間発展にBoris法を適用する運動方程式ソルバ。

### `bem_surface_models`

- kind: `module`
- path: `src/physics/bem_surface_models.f90`
- group: `src/physics`
- 内部依存: `bem_kinds`, `bem_types`, `bem_string_utils`
- external dependencies: なし
- 概要: 表面モデルごとの電荷更新後処理を扱うモジュール。

### `bem_electrostatic_snapshot_eval`

- kind: `submodule`
- path: `src/physics/bem_electrostatic_snapshot_eval.f90`
- group: `src/physics`
- parent: `bem_electrostatic_snapshot`
- 内部依存: `bem_electrostatic_snapshot`
- external dependencies: なし
- 概要: `bem_electrostatic_snapshot` の局所電場・電位評価を実装する submodule。

### `bem_surface_models_conductor`

- kind: `submodule`
- path: `src/physics/bem_surface_models_conductor.f90`
- group: `src/physics`
- parent: `bem_surface_models`
- 内部依存: `bem_surface_models`, `bem_constants`, `bem_panel_geometry`, `bem_panel_kernel`
- external dependencies: なし
- 概要: `bem_surface_models` の浮遊導体電荷再配分を実装する submodule。

### `bem_field_kernel_c`

- kind: `module`
- path: `src/physics/field_solver/bem_field_kernel_c.f90`
- group: `src/physics/field_solver`
- 内部依存: `bem_constants`, `bem_coulomb_fmm_core`, `bem_kinds`, `bem_panel_geometry`, `bem_version`
- external dependencies: `ieee_arithmetic`, `iso_c_binding`
- 概要: C ABI wrapper for the simulator-independent Coulomb FMM field kernel.

### `bem_field_solver`

- kind: `module`
- path: `src/physics/field_solver/bem_field_solver.f90`
- group: `src/physics/field_solver`
- 内部依存: `bem_kinds`, `bem_constants`, `bem_types`, `bem_coulomb_fmm_core`, `bem_string_utils`, `bem_physics_config_types`
- external dependencies: `ieee_arithmetic`
- 概要: 粒子位置での電場評価を direct / treecode / fmm で切り替える場ソルバ。

### `bem_field_solver_config`

- kind: `submodule`
- path: `src/physics/field_solver/bem_field_solver_config.f90`
- group: `src/physics/field_solver`
- parent: `bem_field_solver`
- 内部依存: `bem_field_solver`, `bem_coulomb_fmm_core`, `bem_physics_config_types`
- external dependencies: なし
- 概要: `bem_field_solver` の初期化・設定補助手続きを実装する submodule。

### `bem_field_solver_eval`

- kind: `submodule`
- path: `src/physics/field_solver/bem_field_solver_eval.f90`
- group: `src/physics/field_solver`
- parent: `bem_field_solver`
- 内部依存: `bem_field_solver`, `bem_coulomb_fmm_core`, `bem_panel_geometry`, `bem_panel_kernel`
- external dependencies: なし
- 概要: `bem_field_solver` の電場評価と木走査ロジックを実装する submodule。

### `bem_field_solver_tree`

- kind: `submodule`
- path: `src/physics/field_solver/bem_field_solver_tree.f90`
- group: `src/physics/field_solver`
- parent: `bem_field_solver`
- 内部依存: `bem_field_solver`, `bem_coulomb_fmm_core`
- external dependencies: なし
- 概要: `bem_field_solver` の octree 構築・更新とメモリ管理を実装する submodule。

### `bem_coulomb_fmm_core`

- kind: `module`
- path: `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90`
- group: `src/physics/field_solver/fmm/api`
- 内部依存: `bem_kinds`, `bem_coulomb_fmm_types`
- external dependencies: なし
- 概要: `mesh_type` や `sim_config` に依存しない Coulomb FMM コア API。

### `bem_coulomb_fmm_core_build`

- kind: `submodule`
- path: `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_build.f90`
- group: `src/physics/field_solver/fmm/api`
- parent: `bem_coulomb_fmm_core`
- 内部依存: `bem_coulomb_fmm_core`, `bem_coulomb_fmm_plan_ops`
- external dependencies: なし
- 概要: `bem_coulomb_fmm_core` の plan 構築 API ラッパ。

### `bem_coulomb_fmm_core_eval`

- kind: `submodule`
- path: `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_eval.f90`
- group: `src/physics/field_solver/fmm/api`
- parent: `bem_coulomb_fmm_core`
- 内部依存: `bem_coulomb_fmm_core`, `bem_coulomb_fmm_eval_ops`
- external dependencies: なし
- 概要: `bem_coulomb_fmm_core` の評価 API ラッパ。

### `bem_coulomb_fmm_core_state`

- kind: `submodule`
- path: `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_state.f90`
- group: `src/physics/field_solver/fmm/api`
- parent: `bem_coulomb_fmm_core`
- 内部依存: `bem_coulomb_fmm_core`, `bem_coulomb_fmm_state_ops`
- external dependencies: なし
- 概要: `bem_coulomb_fmm_core` の state 更新 API ラッパ。

### `bem_coulomb_fmm_basis`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/common/bem_coulomb_fmm_basis.f90`
- group: `src/physics/field_solver/fmm/internal/common`
- 内部依存: `bem_kinds`, `bem_coulomb_fmm_types`
- external dependencies: なし
- 概要: Coulomb FMM の multi-index と微分テーブル計算。

### `bem_coulomb_fmm_types`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/common/bem_coulomb_fmm_types.f90`
- group: `src/physics/field_solver/fmm/internal/common`
- 内部依存: `bem_kinds`, `bem_panel_geometry`, `bem_periodic_zero_mode_plan`
- external dependencies: なし
- 概要: Coulomb FMM コアで共有する型定義。

### `bem_coulomb_fmm_periodic`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic.f90`
- group: `src/physics/field_solver/fmm/internal/periodic`
- 内部依存: `bem_kinds`, `bem_coulomb_fmm_types`
- external dependencies: なし
- 概要: Coulomb FMM の periodic2 境界処理。

### `bem_coulomb_fmm_periodic_cache`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_cache.f90`
- group: `src/physics/field_solver/fmm/internal/periodic`
- 内部依存: `bem_filesystem`, `bem_kinds`
- external dependencies: `iso_fortran_env`
- 概要: Versioned stream codec for cached periodic root operators.

### `bem_coulomb_fmm_periodic_ewald`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_ewald.f90`
- group: `src/physics/field_solver/fmm/internal/periodic`
- 内部依存: `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_periodic`
- external dependencies: なし
- 概要: periodic2 cached operator生成用のEwald teacher。

### `bem_coulomb_fmm_periodic_nonzero_reference`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_nonzero_reference.f90`
- group: `src/physics/field_solver/fmm/internal/periodic`
- 内部依存: `bem_kinds`, `bem_constants`, `bem_types`, `bem_panel_geometry`, `bem_panel_quadrature`
- external dependencies: `ieee_arithmetic`

### `bem_coulomb_fmm_periodic_root_ops`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90`
- group: `src/physics/field_solver/fmm/internal/periodic`
- 内部依存: `bem_version`, `bem_kinds`, `bem_mpi`, `bem_filesystem`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_basis`, `bem_coulomb_fmm_periodic`, `bem_coulomb_fmm_periodic_ewald`, `bem_coulomb_fmm_periodic_cache`, `bem_regularized_qr`, `bem_coulomb_fmm_tree_utils`
- external dependencies: なし
- 概要: periodic2 root operator の前計算。

### `bem_regularized_qr`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/periodic/bem_regularized_qr.f90`
- group: `src/physics/field_solver/fmm/internal/periodic`
- 内部依存: `bem_kinds`
- external dependencies: なし
- 概要: Reusable column-scaled QR factorization for ridge-regularized least squares.

### `bem_coulomb_fmm_eval_ops`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90`
- group: `src/physics/field_solver/fmm/internal/runtime`
- 内部依存: `bem_kinds`, `bem_constants`, `bem_panel_geometry`, `bem_panel_kernel`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_basis`, `bem_coulomb_fmm_periodic`, `bem_coulomb_fmm_tree_utils`, `bem_periodic_zero_mode_eval`
- external dependencies: なし
- 概要: Coulomb FMM 電場評価。

### `bem_coulomb_fmm_state_ops`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_state_ops.f90`
- group: `src/physics/field_solver/fmm/internal/runtime`
- 内部依存: `bem_kinds`, `bem_constants`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_tree_utils`, `bem_periodic_zero_mode_plan`
- external dependencies: なし
- 概要: Coulomb FMM state 更新と upward/downward pass。

### `bem_coulomb_fmm_plan_ops`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/tree/bem_coulomb_fmm_plan_ops.f90`
- group: `src/physics/field_solver/fmm/internal/tree`
- 内部依存: `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_basis`, `bem_coulomb_fmm_periodic`, `bem_coulomb_fmm_periodic_ewald`, `bem_coulomb_fmm_periodic_root_ops`, `bem_panel_geometry`, `bem_periodic_zero_mode_plan`, `bem_coulomb_fmm_tree_utils`
- external dependencies: なし
- 概要: Coulomb FMM plan 構築と tree トポロジ前計算。

### `bem_coulomb_fmm_tree_utils`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/tree/bem_coulomb_fmm_tree_utils.f90`
- group: `src/physics/field_solver/fmm/internal/tree`
- 内部依存: `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_periodic`
- external dependencies: なし
- 概要: Coulomb FMM tree 構造の共通ユーティリティ。

### `bem_panel_geometry`

- kind: `module`
- path: `src/physics/panel/bem_panel_geometry.f90`
- group: `src/physics/panel`
- 内部依存: `bem_kinds`
- external dependencies: `ieee_arithmetic`
- 概要: Ordered triangle geometry and exact Cartesian surface moments.

### `bem_panel_kernel`

- kind: `module`
- path: `src/physics/panel/bem_panel_kernel.f90`
- group: `src/physics/panel`
- 内部依存: `bem_kinds`, `bem_constants`, `bem_panel_geometry`, `bem_panel_self_terms`
- external dependencies: なし
- 概要: Analytic free-space P0 triangle potential, field, principal value, and jump.

### `bem_panel_quadrature`

- kind: `module`
- path: `src/physics/panel/bem_panel_quadrature.f90`
- group: `src/physics/panel`
- 内部依存: `bem_kinds`, `bem_constants`, `bem_panel_geometry`
- external dependencies: なし
- 概要: Independent triangle cubature and Gauss-Duffy correctness oracles.

### `bem_panel_self_terms`

- kind: `module`
- path: `src/physics/panel/bem_panel_self_terms.f90`
- group: `src/physics/panel`
- 内部依存: `bem_kinds`, `bem_panel_geometry`
- external dependencies: なし
- 概要: On-surface P0 triangle potential and principal-value field integrals.

### `bem_panel_surface_sides`

- kind: `module`
- path: `src/physics/panel/bem_panel_surface_sides.f90`
- group: `src/physics/panel`
- 内部依存: `bem_kinds`, `bem_types`, `bem_string_utils`
- external dependencies: なし
- 概要: Resolve physical vacuum sides without changing ordered triangle winding.

### `bem_periodic_zero_mode_c`

- kind: `module`
- path: `src/physics/periodic_zero_mode/bem_periodic_zero_mode_c.f90`
- group: `src/physics/periodic_zero_mode`
- 内部依存: `bem_kinds`, `bem_periodic_zero_mode_eval`, `bem_periodic_zero_mode_plan`
- external dependencies: `ieee_arithmetic`, `iso_c_binding`
- 概要: C ABI wrapper for the physical periodic zero mode.

### `bem_periodic_zero_mode_eval`

- kind: `module`
- path: `src/physics/periodic_zero_mode/bem_periodic_zero_mode_eval.f90`
- group: `src/physics/periodic_zero_mode`
- 内部依存: `bem_kinds`, `bem_constants`, `bem_periodic_zero_mode_plan`
- external dependencies: なし

### `bem_periodic_zero_mode_plan`

- kind: `module`
- path: `src/physics/periodic_zero_mode/bem_periodic_zero_mode_plan.f90`
- group: `src/physics/periodic_zero_mode`
- 内部依存: `bem_kinds`, `bem_constants`, `bem_types`
- external dependencies: `ieee_arithmetic`

### `bem_checkpoint_contract`

- kind: `module`
- path: `src/runtime/bem_checkpoint_contract.f90`
- group: `src/runtime`
- 内部依存: `bem_kinds`
- external dependencies: なし
- 概要: BEACH checkpoint metadata shared by output and restart code.

### `bem_filesystem`

- kind: `module`
- path: `src/runtime/bem_filesystem.f90`
- group: `src/runtime`
- 内部依存: なし
- external dependencies: `iso_c_binding`
- 概要: Minimal POSIX filesystem operations used by the runtime.

### `bem_model_fingerprint`

- kind: `module`
- path: `src/runtime/bem_model_fingerprint.f90`
- group: `src/runtime`
- 内部依存: `bem_kinds`, `bem_types`, `bem_app_config_types`
- external dependencies: なし
- 概要: Restart compatibility fingerprints for the ordered physical model contract.

### `bem_output_writer`

- kind: `module`
- path: `src/runtime/bem_output_writer.f90`
- group: `src/runtime`
- 内部依存: `bem_kinds`, `bem_types`, `bem_app_config_types`, `bem_charge_ledger`, `bem_checkpoint_contract`, `bem_electrostatic_snapshot`, `bem_external_boundary_contract`, `bem_model_fingerprint`, `bem_version`, `bem_filesystem`, `bem_string_utils`
- external dependencies: なし
- 概要: 実行サマリ・最終CSV・履歴CSVの出力を担当するモジュール。

### `bem_performance_profile`

- kind: `module`
- path: `src/runtime/bem_performance_profile.f90`
- group: `src/runtime`
- 内部依存: `bem_kinds`, `bem_mpi`, `bem_string_utils`
- external dependencies: `iso_fortran_env`
- 概要: 実行フェーズごとの壁時計計測と MPI 集約出力を担う軽量プロファイラ。

### `bem_restart`

- kind: `module`
- path: `src/runtime/bem_restart.f90`
- group: `src/runtime`
- 内部依存: `bem_kinds`, `bem_types`, `bem_app_config_types`, `bem_charge_ledger`, `bem_model_fingerprint`, `bem_physics_config_types`, `bem_checkpoint_contract`, `bem_mpi`
- external dependencies: `ieee_arithmetic`
- 概要: チェックポイントファイルの保存/復元を扱う補助モジュール。

### `bem_charge_ledger`

- kind: `module`
- path: `src/runtime/coupling/bem_charge_ledger.f90`
- group: `src/runtime/coupling`
- 内部依存: `bem_kinds`
- external dependencies: なし
- 概要: batch 間の signed charge stock と移送 flux から電荷収支を集計する。

### `bem_particle_stepper`

- kind: `module`
- path: `src/runtime/simulator/bem_particle_stepper.f90`
- group: `src/runtime/simulator`
- 内部依存: `bem_kinds`, `bem_types`, `bem_electrostatic_snapshot`, `bem_pusher`, `bem_collision`, `bem_boundary`, `bem_external_boundary_contract`
- external dependencies: `ieee_arithmetic`
- 概要: 同一時刻の粒子状態から、空間電場を中点評価した1ステップ候補を構築する。

### `bem_simulator`

- kind: `module`
- path: `src/runtime/simulator/bem_simulator.f90`
- group: `src/runtime/simulator`
- 内部依存: `bem_kinds`, `bem_types`, `bem_app_config`, `bem_app_config_runtime`, `bem_electrostatic_snapshot`, `bem_particle_stepper`, `bem_collision`, `bem_surface_models`, `bem_charge_ledger`, `bem_string_utils`, `bem_external_boundary_contract`, `bem_simulator_workspace`, `bem_mpi`
- external dependencies: `iso_fortran_env`, `ieee_arithmetic`
- 概要: 吸着(insulator)モデルのメインループを実行し、電荷堆積と統計更新を行う。

### `bem_simulator_workspace`

- kind: `module`
- path: `src/runtime/simulator/bem_simulator_workspace.f90`
- group: `src/runtime/simulator`
- 内部依存: `bem_kinds`
- external dependencies: なし
- 概要: シミュレーション実行中に再利用するバッチ作業配列を管理する。

### `bem_simulator_io`

- kind: `submodule`
- path: `src/runtime/simulator/bem_simulator_io.f90`
- group: `src/runtime/simulator`
- parent: `bem_simulator`
- 内部依存: `bem_simulator`, `bem_app_config_runtime`, `bem_output_writer`
- external dependencies: なし
- 概要: `bem_simulator` の進捗表示と履歴出力を実装する submodule。

### `bem_simulator_loop`

- kind: `submodule`
- path: `src/runtime/simulator/bem_simulator_loop.f90`
- group: `src/runtime/simulator`
- parent: `bem_simulator`
- 内部依存: `bem_simulator`, `bem_app_config_runtime`, `bem_periodic_zero_mode_plan`, `bem_performance_profile`, `bem_mpi`
- external dependencies: `iso_fortran_env`
- 概要: `bem_simulator` の主ループと粒子処理計算を実装する submodule。

### `bem_simulator_stats`

- kind: `submodule`
- path: `src/runtime/simulator/bem_simulator_stats.f90`
- group: `src/runtime/simulator`
- parent: `bem_simulator`
- 内部依存: `bem_simulator`
- external dependencies: なし
- 概要: `bem_simulator` のバッチ集計・統計更新処理を実装する submodule。

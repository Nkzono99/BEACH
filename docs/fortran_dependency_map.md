title: Fortran 依存関係マップ

# Fortran 依存関係マップ

> このページは `tools/generate_fortran_dependency_report.py` から自動生成しています。

## 概要

- ソースファイル数: 47
- モジュール数: 35
- submodule 数: 11
- program 数: 1
- 内部依存エッジ数: 147

## 全体グラフ

![Fortran モジュール依存グラフ](../media/fortran_module_dependencies.svg)

実線は `use` 依存、破線は `submodule(parent)` の親参照を表します。

## ディレクトリ別サマリ

| ディレクトリ | エンティティ数 | 内部依存数 |
| --- | ---: | ---: |
| `app` | 1 | 9 |
| `src/config` | 3 | 17 |
| `src/config/app_config_parser` | 3 | 6 |
| `src/core` | 4 | 3 |
| `src/mesh` | 3 | 8 |
| `src/particles` | 3 | 9 |
| `src/physics` | 4 | 8 |
| `src/physics/field_solver` | 4 | 11 |
| `src/physics/field_solver/fmm/api` | 4 | 8 |
| `src/physics/field_solver/fmm/internal/common` | 2 | 3 |
| `src/physics/field_solver/fmm/internal/periodic` | 3 | 11 |
| `src/physics/field_solver/fmm/internal/runtime` | 2 | 9 |
| `src/physics/field_solver/fmm/internal/tree` | 2 | 10 |
| `src/physics/sheath` | 2 | 9 |
| `src/runtime` | 3 | 13 |
| `src/runtime/simulator` | 4 | 13 |

## 被依存の多いモジュール

| エンティティ | kind | 被依存数 |
| --- | --- | ---: |
| `bem_kinds` | `module` | 33 |
| `bem_types` | `module` | 17 |
| `bem_coulomb_fmm_types` | `module` | 10 |
| `bem_coulomb_fmm_core` | `module` | 8 |
| `bem_constants` | `module` | 7 |
| `bem_coulomb_fmm_periodic` | `module` | 6 |
| `bem_app_config_parser` | `module` | 5 |
| `bem_app_config_types` | `module` | 5 |
| `bem_mpi` | `module` | 5 |
| `bem_coulomb_fmm_periodic_ewald` | `module` | 4 |

## エンティティ一覧

| エンティティ | kind | パス | 内部依存 | 概要 |
| --- | --- | --- | --- | --- |
| `main` | `program` | `app/main.f90` | `bem_kinds`, `bem_types`, `bem_mpi`, `bem_performance_profile`, `bem_simulator`, `bem_restart`, `bem_output_writer`, `bem_app_config`, `bem_mesh` | 設定読込・メッシュ生成・粒子初期化・シミュレーション実行・結果出力を順に行うCLIエントリーポイント。 |
| `bem_app_config` | `module` | `src/config/bem_app_config.f90` | `bem_app_config_types`, `bem_app_config_parser`, `bem_app_config_runtime` | 設定型・TOMLパーサ・実行時変換ロジックを束ねる後方互換ファサード。 |
| `bem_app_config_runtime` | `module` | `src/config/bem_app_config_runtime.f90` | `bem_kinds`, `bem_types`, `bem_mpi`, `bem_field`, `bem_templates`, `bem_mesh`, `bem_importers`, `bem_injection`, `bem_particles`, `bem_sheath_injection_model`, `bem_app_config_types`, `bem_app_config_parser` | `app_config` からメッシュ・粒子群を構築する実行時変換モジュール。 |
| `bem_app_config_types` | `module` | `src/config/bem_app_config_types.f90` | `bem_kinds`, `bem_types` | アプリ設定の型定義と、設定由来の粒子数計算をまとめるモジュール。 |
| `bem_app_config_parser` | `module` | `src/config/app_config_parser/bem_app_config_parser.f90` | `bem_kinds`, `bem_constants`, `bem_types`, `bem_app_config_types` | TOML風設定ファイルを `app_config` へ読み込む軽量パーサ。 |
| `bem_app_config_parser_parse_utils` | `submodule` | `src/config/app_config_parser/bem_app_config_parser_parse_utils.f90` | `bem_app_config_parser` | `bem_app_config_parser` の文字列パース補助手続きを実装する submodule。 |
| `bem_app_config_parser_validate` | `submodule` | `src/config/app_config_parser/bem_app_config_parser_validate.f90` | `bem_app_config_parser` | `bem_app_config_parser` の入力検証・物理量導出手続きを実装する submodule。 |
| `bem_constants` | `module` | `src/core/bem_constants.f90` | `bem_kinds` | シミュレーションで使用する物理定数を定義する。 |
| `bem_kinds` | `module` | `src/core/bem_kinds.f90` | - | 倍精度実数と32bit整数のkind定義を集約する基盤モジュール。 |
| `bem_mpi` | `module` | `src/core/bem_mpi.F90` | `bem_kinds` | MPIの初期化・集約を抽象化し、非MPIビルドでは単一ランク動作へフォールバックする。 |
| `bem_types` | `module` | `src/core/bem_types.f90` | `bem_kinds` | シミュレーション設定・統計・メッシュ・粒子・衝突情報の主要データ型を定義する。 |
| `bem_importers` | `module` | `src/mesh/bem_importers.f90` | `bem_kinds`, `bem_types`, `bem_mesh` | OBJメッシュを走査・解析し、内部 `mesh_type` へ変換するインポートモジュール。 |
| `bem_mesh` | `module` | `src/mesh/bem_mesh.f90` | `bem_kinds`, `bem_types` | 三角形メッシュ幾何量(重心・法線・AABB・代表長)を前計算して保持するモジュール。 |
| `bem_templates` | `module` | `src/mesh/bem_templates.f90` | `bem_kinds`, `bem_types`, `bem_mesh` | 平面/穴あき平面/円板/リング/箱/円柱/球テンプレートから三角形メッシュを生成するユーティリティ。 |
| `bem_injection` | `module` | `src/particles/bem_injection.f90` | `bem_kinds`, `bem_constants`, `bem_particles`, `bem_types`, `bem_boundary`, `bem_collision` | 乱数シード設定と粒子位置/速度サンプリングを担う粒子注入モジュール。 |
| `bem_particles` | `module` | `src/particles/bem_particles.f90` | `bem_kinds`, `bem_types` | 粒子SoAデータ構造の初期化を提供するモジュール。 |
| `bem_sheath_injection_model` | `module` | `src/particles/bem_sheath_injection_model.f90` | `bem_sheath_runtime` | 互換性維持のためのシース注入ラッパモジュール。 |
| `bem_boundary` | `module` | `src/physics/bem_boundary.f90` | `bem_kinds`, `bem_types` | シミュレーションボックス境界（流出/反射/周期）を適用するモジュール。 |
| `bem_collision` | `module` | `src/physics/bem_collision.f90` | `bem_kinds`, `bem_types` | 粒子軌道セグメントと三角形要素の交差判定を提供する衝突検出モジュール。 |
| `bem_field` | `module` | `src/physics/bem_field.f90` | `bem_kinds`, `bem_constants`, `bem_types` | 境界要素に蓄積した電荷から観測点の電場を評価する場計算モジュール。 |
| `bem_pusher` | `module` | `src/physics/bem_pusher.f90` | `bem_kinds` | 荷電粒子の時間発展にBoris法を適用する運動方程式ソルバ。 |
| `bem_field_solver` | `module` | `src/physics/field_solver/bem_field_solver.f90` | `bem_kinds`, `bem_constants`, `bem_types`, `bem_field`, `bem_coulomb_fmm_core` | 粒子位置での電場評価を direct / treecode / fmm で切り替える場ソルバ。 |
| `bem_field_solver_config` | `submodule` | `src/physics/field_solver/bem_field_solver_config.f90` | `bem_field_solver`, `bem_coulomb_fmm_core` | `bem_field_solver` の初期化・設定補助手続きを実装する submodule。 |
| `bem_field_solver_eval` | `submodule` | `src/physics/field_solver/bem_field_solver_eval.f90` | `bem_field_solver`, `bem_coulomb_fmm_core` | `bem_field_solver` の電場評価と木走査ロジックを実装する submodule。 |
| `bem_field_solver_tree` | `submodule` | `src/physics/field_solver/bem_field_solver_tree.f90` | `bem_field_solver`, `bem_coulomb_fmm_core` | `bem_field_solver` の octree 構築・更新とメモリ管理を実装する submodule。 |
| `bem_coulomb_fmm_core` | `module` | `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90` | `bem_kinds`, `bem_coulomb_fmm_types` | `mesh_type` や `sim_config` に依存しない Coulomb FMM コア API。 |
| `bem_coulomb_fmm_core_build` | `submodule` | `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_build.f90` | `bem_coulomb_fmm_core`, `bem_coulomb_fmm_plan_ops` | `bem_coulomb_fmm_core` の plan 構築 API ラッパ。 |
| `bem_coulomb_fmm_core_eval` | `submodule` | `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_eval.f90` | `bem_coulomb_fmm_core`, `bem_coulomb_fmm_eval_ops` | `bem_coulomb_fmm_core` の評価 API ラッパ。 |
| `bem_coulomb_fmm_core_state` | `submodule` | `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_state.f90` | `bem_coulomb_fmm_core`, `bem_coulomb_fmm_state_ops` | `bem_coulomb_fmm_core` の state 更新 API ラッパ。 |
| `bem_coulomb_fmm_basis` | `module` | `src/physics/field_solver/fmm/internal/common/bem_coulomb_fmm_basis.f90` | `bem_kinds`, `bem_coulomb_fmm_types` | Coulomb FMM の multi-index と微分テーブル計算。 |
| `bem_coulomb_fmm_types` | `module` | `src/physics/field_solver/fmm/internal/common/bem_coulomb_fmm_types.f90` | `bem_kinds` | Coulomb FMM コアで共有する型定義。 |
| `bem_coulomb_fmm_periodic` | `module` | `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic.f90` | `bem_kinds`, `bem_coulomb_fmm_types` | Coulomb FMM の periodic2 境界処理。 |
| `bem_coulomb_fmm_periodic_ewald` | `module` | `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_ewald.f90` | `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_periodic` | periodic2 build-only Ewald oracle と fallback exact correction。 |
| `bem_coulomb_fmm_periodic_root_ops` | `module` | `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90` | `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_basis`, `bem_coulomb_fmm_periodic`, `bem_coulomb_fmm_periodic_ewald`, `bem_coulomb_fmm_tree_utils` | periodic2 root operator の前計算。 |
| `bem_coulomb_fmm_eval_ops` | `module` | `src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90` | `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_basis`, `bem_coulomb_fmm_periodic`, `bem_coulomb_fmm_periodic_ewald`, `bem_coulomb_fmm_tree_utils` | Coulomb FMM 電場評価。 |
| `bem_coulomb_fmm_state_ops` | `module` | `src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_state_ops.f90` | `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_tree_utils` | Coulomb FMM state 更新と upward/downward pass。 |
| `bem_coulomb_fmm_plan_ops` | `module` | `src/physics/field_solver/fmm/internal/tree/bem_coulomb_fmm_plan_ops.f90` | `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_basis`, `bem_coulomb_fmm_periodic`, `bem_coulomb_fmm_periodic_ewald`, `bem_coulomb_fmm_periodic_root_ops`, `bem_coulomb_fmm_tree_utils` | Coulomb FMM plan 構築と tree トポロジ前計算。 |
| `bem_coulomb_fmm_tree_utils` | `module` | `src/physics/field_solver/fmm/internal/tree/bem_coulomb_fmm_tree_utils.f90` | `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_periodic` | Coulomb FMM tree 構造の共通ユーティリティ。 |
| `bem_sheath_model_core` | `module` | `src/physics/sheath/bem_sheath_model_core.f90` | `bem_kinds`, `bem_constants`, `bem_injection` | Zhao 系シース数値モデルの core 実装。 |
| `bem_sheath_runtime` | `module` | `src/physics/sheath/bem_sheath_runtime.f90` | `bem_kinds`, `bem_constants`, `bem_types`, `bem_app_config_types`, `bem_app_config_parser`, `bem_sheath_model_core` | シース数値モデルと app_config / 注入ランタイムの橋渡しを行うモジュール。 |
| `bem_output_writer` | `module` | `src/runtime/bem_output_writer.f90` | `bem_kinds`, `bem_constants`, `bem_types`, `bem_app_config_types`, `bem_coulomb_fmm_core`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_periodic`, `bem_coulomb_fmm_periodic_ewald` | 実行サマリ・最終CSV・履歴CSVの出力を担当するモジュール。 |
| `bem_performance_profile` | `module` | `src/runtime/bem_performance_profile.f90` | `bem_kinds`, `bem_mpi` | 実行フェーズごとの壁時計計測と MPI 集約出力を担う軽量プロファイラ。 |
| `bem_restart` | `module` | `src/runtime/bem_restart.f90` | `bem_kinds`, `bem_types`, `bem_mpi` | 出力ディレクトリに保存したチェックポイントの保存/復元を扱う補助モジュール。 |
| `bem_simulator` | `module` | `src/runtime/simulator/bem_simulator.f90` | `bem_kinds`, `bem_types`, `bem_app_config`, `bem_field_solver`, `bem_pusher`, `bem_collision`, `bem_boundary`, `bem_mpi` | 吸着(insulator)モデルのメインループを実行し、電荷堆積と統計更新を行うモジュール。 |
| `bem_simulator_io` | `submodule` | `src/runtime/simulator/bem_simulator_io.f90` | `bem_simulator` | `bem_simulator` の進捗表示と履歴出力を実装する submodule。 |
| `bem_simulator_loop` | `submodule` | `src/runtime/simulator/bem_simulator_loop.f90` | `bem_simulator`, `bem_output_writer`, `bem_performance_profile` | `bem_simulator` の主ループと粒子処理計算を実装する submodule。 |
| `bem_simulator_stats` | `submodule` | `src/runtime/simulator/bem_simulator_stats.f90` | `bem_simulator` | `bem_simulator` のバッチ集計・統計更新処理を実装する submodule。 |

## 詳細

### `main`

- kind: `program`
- path: `app/main.f90`
- group: `app`
- internal dependencies: `bem_kinds`, `bem_types`, `bem_mpi`, `bem_performance_profile`, `bem_simulator`, `bem_restart`, `bem_output_writer`, `bem_app_config`, `bem_mesh`
- external dependencies: なし
- summary: 設定読込・メッシュ生成・粒子初期化・シミュレーション実行・結果出力を順に行うCLIエントリーポイント。

### `bem_app_config`

- kind: `module`
- path: `src/config/bem_app_config.f90`
- group: `src/config`
- internal dependencies: `bem_app_config_types`, `bem_app_config_parser`, `bem_app_config_runtime`
- external dependencies: なし
- summary: 設定型・TOMLパーサ・実行時変換ロジックを束ねる後方互換ファサード。

### `bem_app_config_runtime`

- kind: `module`
- path: `src/config/bem_app_config_runtime.f90`
- group: `src/config`
- internal dependencies: `bem_kinds`, `bem_types`, `bem_mpi`, `bem_field`, `bem_templates`, `bem_mesh`, `bem_importers`, `bem_injection`, `bem_particles`, `bem_sheath_injection_model`, `bem_app_config_types`, `bem_app_config_parser`
- external dependencies: `ieee_arithmetic`
- summary: `app_config` からメッシュ・粒子群を構築する実行時変換モジュール。

### `bem_app_config_types`

- kind: `module`
- path: `src/config/bem_app_config_types.f90`
- group: `src/config`
- internal dependencies: `bem_kinds`, `bem_types`
- external dependencies: なし
- summary: アプリ設定の型定義と、設定由来の粒子数計算をまとめるモジュール。

### `bem_app_config_parser`

- kind: `module`
- path: `src/config/app_config_parser/bem_app_config_parser.f90`
- group: `src/config/app_config_parser`
- internal dependencies: `bem_kinds`, `bem_constants`, `bem_types`, `bem_app_config_types`
- external dependencies: `ieee_arithmetic`
- summary: TOML風設定ファイルを `app_config` へ読み込む軽量パーサ。

### `bem_app_config_parser_parse_utils`

- kind: `submodule`
- path: `src/config/app_config_parser/bem_app_config_parser_parse_utils.f90`
- group: `src/config/app_config_parser`
- parent: `bem_app_config_parser`
- internal dependencies: `bem_app_config_parser`
- external dependencies: なし
- summary: `bem_app_config_parser` の文字列パース補助手続きを実装する submodule。

### `bem_app_config_parser_validate`

- kind: `submodule`
- path: `src/config/app_config_parser/bem_app_config_parser_validate.f90`
- group: `src/config/app_config_parser`
- parent: `bem_app_config_parser`
- internal dependencies: `bem_app_config_parser`
- external dependencies: なし
- summary: `bem_app_config_parser` の入力検証・物理量導出手続きを実装する submodule。

### `bem_constants`

- kind: `module`
- path: `src/core/bem_constants.f90`
- group: `src/core`
- internal dependencies: `bem_kinds`
- external dependencies: なし
- summary: シミュレーションで使用する物理定数を定義する。

### `bem_kinds`

- kind: `module`
- path: `src/core/bem_kinds.f90`
- group: `src/core`
- internal dependencies: なし
- external dependencies: `iso_fortran_env`
- summary: 倍精度実数と32bit整数のkind定義を集約する基盤モジュール。

### `bem_mpi`

- kind: `module`
- path: `src/core/bem_mpi.F90`
- group: `src/core`
- internal dependencies: `bem_kinds`
- external dependencies: なし
- summary: MPIの初期化・集約を抽象化し、非MPIビルドでは単一ランク動作へフォールバックする。

### `bem_types`

- kind: `module`
- path: `src/core/bem_types.f90`
- group: `src/core`
- internal dependencies: `bem_kinds`
- external dependencies: なし
- summary: シミュレーション設定・統計・メッシュ・粒子・衝突情報の主要データ型を定義する。

### `bem_importers`

- kind: `module`
- path: `src/mesh/bem_importers.f90`
- group: `src/mesh`
- internal dependencies: `bem_kinds`, `bem_types`, `bem_mesh`
- external dependencies: なし
- summary: OBJメッシュを走査・解析し、内部 `mesh_type` へ変換するインポートモジュール。

### `bem_mesh`

- kind: `module`
- path: `src/mesh/bem_mesh.f90`
- group: `src/mesh`
- internal dependencies: `bem_kinds`, `bem_types`
- external dependencies: なし
- summary: 三角形メッシュ幾何量(重心・法線・AABB・代表長)を前計算して保持するモジュール。

### `bem_templates`

- kind: `module`
- path: `src/mesh/bem_templates.f90`
- group: `src/mesh`
- internal dependencies: `bem_kinds`, `bem_types`, `bem_mesh`
- external dependencies: なし
- summary: 平面/穴あき平面/円板/リング/箱/円柱/球テンプレートから三角形メッシュを生成するユーティリティ。

### `bem_injection`

- kind: `module`
- path: `src/particles/bem_injection.f90`
- group: `src/particles`
- internal dependencies: `bem_kinds`, `bem_constants`, `bem_particles`, `bem_types`, `bem_boundary`, `bem_collision`
- external dependencies: なし
- summary: 乱数シード設定と粒子位置/速度サンプリングを担う粒子注入モジュール。

### `bem_particles`

- kind: `module`
- path: `src/particles/bem_particles.f90`
- group: `src/particles`
- internal dependencies: `bem_kinds`, `bem_types`
- external dependencies: なし
- summary: 粒子SoAデータ構造の初期化を提供するモジュール。

### `bem_sheath_injection_model`

- kind: `module`
- path: `src/particles/bem_sheath_injection_model.f90`
- group: `src/particles`
- internal dependencies: `bem_sheath_runtime`
- external dependencies: なし
- summary: 互換性維持のためのシース注入ラッパモジュール。

### `bem_boundary`

- kind: `module`
- path: `src/physics/bem_boundary.f90`
- group: `src/physics`
- internal dependencies: `bem_kinds`, `bem_types`
- external dependencies: なし
- summary: シミュレーションボックス境界（流出/反射/周期）を適用するモジュール。

### `bem_collision`

- kind: `module`
- path: `src/physics/bem_collision.f90`
- group: `src/physics`
- internal dependencies: `bem_kinds`, `bem_types`
- external dependencies: なし
- summary: 粒子軌道セグメントと三角形要素の交差判定を提供する衝突検出モジュール。

### `bem_field`

- kind: `module`
- path: `src/physics/bem_field.f90`
- group: `src/physics`
- internal dependencies: `bem_kinds`, `bem_constants`, `bem_types`
- external dependencies: なし
- summary: 境界要素に蓄積した電荷から観測点の電場を評価する場計算モジュール。

### `bem_pusher`

- kind: `module`
- path: `src/physics/bem_pusher.f90`
- group: `src/physics`
- internal dependencies: `bem_kinds`
- external dependencies: なし
- summary: 荷電粒子の時間発展にBoris法を適用する運動方程式ソルバ。

### `bem_field_solver`

- kind: `module`
- path: `src/physics/field_solver/bem_field_solver.f90`
- group: `src/physics/field_solver`
- internal dependencies: `bem_kinds`, `bem_constants`, `bem_types`, `bem_field`, `bem_coulomb_fmm_core`
- external dependencies: なし
- summary: 粒子位置での電場評価を direct / treecode / fmm で切り替える場ソルバ。

### `bem_field_solver_config`

- kind: `submodule`
- path: `src/physics/field_solver/bem_field_solver_config.f90`
- group: `src/physics/field_solver`
- parent: `bem_field_solver`
- internal dependencies: `bem_field_solver`, `bem_coulomb_fmm_core`
- external dependencies: なし
- summary: `bem_field_solver` の初期化・設定補助手続きを実装する submodule。

### `bem_field_solver_eval`

- kind: `submodule`
- path: `src/physics/field_solver/bem_field_solver_eval.f90`
- group: `src/physics/field_solver`
- parent: `bem_field_solver`
- internal dependencies: `bem_field_solver`, `bem_coulomb_fmm_core`
- external dependencies: なし
- summary: `bem_field_solver` の電場評価と木走査ロジックを実装する submodule。

### `bem_field_solver_tree`

- kind: `submodule`
- path: `src/physics/field_solver/bem_field_solver_tree.f90`
- group: `src/physics/field_solver`
- parent: `bem_field_solver`
- internal dependencies: `bem_field_solver`, `bem_coulomb_fmm_core`
- external dependencies: なし
- summary: `bem_field_solver` の octree 構築・更新とメモリ管理を実装する submodule。

### `bem_coulomb_fmm_core`

- kind: `module`
- path: `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90`
- group: `src/physics/field_solver/fmm/api`
- internal dependencies: `bem_kinds`, `bem_coulomb_fmm_types`
- external dependencies: なし
- summary: `mesh_type` や `sim_config` に依存しない Coulomb FMM コア API。

### `bem_coulomb_fmm_core_build`

- kind: `submodule`
- path: `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_build.f90`
- group: `src/physics/field_solver/fmm/api`
- parent: `bem_coulomb_fmm_core`
- internal dependencies: `bem_coulomb_fmm_core`, `bem_coulomb_fmm_plan_ops`
- external dependencies: なし
- summary: `bem_coulomb_fmm_core` の plan 構築 API ラッパ。

### `bem_coulomb_fmm_core_eval`

- kind: `submodule`
- path: `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_eval.f90`
- group: `src/physics/field_solver/fmm/api`
- parent: `bem_coulomb_fmm_core`
- internal dependencies: `bem_coulomb_fmm_core`, `bem_coulomb_fmm_eval_ops`
- external dependencies: なし
- summary: `bem_coulomb_fmm_core` の評価 API ラッパ。

### `bem_coulomb_fmm_core_state`

- kind: `submodule`
- path: `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_state.f90`
- group: `src/physics/field_solver/fmm/api`
- parent: `bem_coulomb_fmm_core`
- internal dependencies: `bem_coulomb_fmm_core`, `bem_coulomb_fmm_state_ops`
- external dependencies: なし
- summary: `bem_coulomb_fmm_core` の state 更新 API ラッパ。

### `bem_coulomb_fmm_basis`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/common/bem_coulomb_fmm_basis.f90`
- group: `src/physics/field_solver/fmm/internal/common`
- internal dependencies: `bem_kinds`, `bem_coulomb_fmm_types`
- external dependencies: なし
- summary: Coulomb FMM の multi-index と微分テーブル計算。

### `bem_coulomb_fmm_types`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/common/bem_coulomb_fmm_types.f90`
- group: `src/physics/field_solver/fmm/internal/common`
- internal dependencies: `bem_kinds`
- external dependencies: なし
- summary: Coulomb FMM コアで共有する型定義。

### `bem_coulomb_fmm_periodic`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic.f90`
- group: `src/physics/field_solver/fmm/internal/periodic`
- internal dependencies: `bem_kinds`, `bem_coulomb_fmm_types`
- external dependencies: なし
- summary: Coulomb FMM の periodic2 境界処理。

### `bem_coulomb_fmm_periodic_ewald`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_ewald.f90`
- group: `src/physics/field_solver/fmm/internal/periodic`
- internal dependencies: `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_periodic`
- external dependencies: なし
- summary: periodic2 build-only Ewald oracle と fallback exact correction。

### `bem_coulomb_fmm_periodic_root_ops`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90`
- group: `src/physics/field_solver/fmm/internal/periodic`
- internal dependencies: `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_basis`, `bem_coulomb_fmm_periodic`, `bem_coulomb_fmm_periodic_ewald`, `bem_coulomb_fmm_tree_utils`
- external dependencies: なし
- summary: periodic2 root operator の前計算。

### `bem_coulomb_fmm_eval_ops`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90`
- group: `src/physics/field_solver/fmm/internal/runtime`
- internal dependencies: `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_basis`, `bem_coulomb_fmm_periodic`, `bem_coulomb_fmm_periodic_ewald`, `bem_coulomb_fmm_tree_utils`
- external dependencies: なし
- summary: Coulomb FMM 電場評価。

### `bem_coulomb_fmm_state_ops`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_state_ops.f90`
- group: `src/physics/field_solver/fmm/internal/runtime`
- internal dependencies: `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_tree_utils`
- external dependencies: なし
- summary: Coulomb FMM state 更新と upward/downward pass。

### `bem_coulomb_fmm_plan_ops`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/tree/bem_coulomb_fmm_plan_ops.f90`
- group: `src/physics/field_solver/fmm/internal/tree`
- internal dependencies: `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_basis`, `bem_coulomb_fmm_periodic`, `bem_coulomb_fmm_periodic_ewald`, `bem_coulomb_fmm_periodic_root_ops`, `bem_coulomb_fmm_tree_utils`
- external dependencies: なし
- summary: Coulomb FMM plan 構築と tree トポロジ前計算。

### `bem_coulomb_fmm_tree_utils`

- kind: `module`
- path: `src/physics/field_solver/fmm/internal/tree/bem_coulomb_fmm_tree_utils.f90`
- group: `src/physics/field_solver/fmm/internal/tree`
- internal dependencies: `bem_kinds`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_periodic`
- external dependencies: なし
- summary: Coulomb FMM tree 構造の共通ユーティリティ。

### `bem_sheath_model_core`

- kind: `module`
- path: `src/physics/sheath/bem_sheath_model_core.f90`
- group: `src/physics/sheath`
- internal dependencies: `bem_kinds`, `bem_constants`, `bem_injection`
- external dependencies: `ieee_arithmetic`
- summary: Zhao 系シース数値モデルの core 実装。

### `bem_sheath_runtime`

- kind: `module`
- path: `src/physics/sheath/bem_sheath_runtime.f90`
- group: `src/physics/sheath`
- internal dependencies: `bem_kinds`, `bem_constants`, `bem_types`, `bem_app_config_types`, `bem_app_config_parser`, `bem_sheath_model_core`
- external dependencies: なし
- summary: シース数値モデルと app_config / 注入ランタイムの橋渡しを行うモジュール。

### `bem_output_writer`

- kind: `module`
- path: `src/runtime/bem_output_writer.f90`
- group: `src/runtime`
- internal dependencies: `bem_kinds`, `bem_constants`, `bem_types`, `bem_app_config_types`, `bem_coulomb_fmm_core`, `bem_coulomb_fmm_types`, `bem_coulomb_fmm_periodic`, `bem_coulomb_fmm_periodic_ewald`
- external dependencies: なし
- summary: 実行サマリ・最終CSV・履歴CSVの出力を担当するモジュール。

### `bem_performance_profile`

- kind: `module`
- path: `src/runtime/bem_performance_profile.f90`
- group: `src/runtime`
- internal dependencies: `bem_kinds`, `bem_mpi`
- external dependencies: `iso_fortran_env`
- summary: 実行フェーズごとの壁時計計測と MPI 集約出力を担う軽量プロファイラ。

### `bem_restart`

- kind: `module`
- path: `src/runtime/bem_restart.f90`
- group: `src/runtime`
- internal dependencies: `bem_kinds`, `bem_types`, `bem_mpi`
- external dependencies: なし
- summary: 出力ディレクトリに保存したチェックポイントの保存/復元を扱う補助モジュール。

### `bem_simulator`

- kind: `module`
- path: `src/runtime/simulator/bem_simulator.f90`
- group: `src/runtime/simulator`
- internal dependencies: `bem_kinds`, `bem_types`, `bem_app_config`, `bem_field_solver`, `bem_pusher`, `bem_collision`, `bem_boundary`, `bem_mpi`
- external dependencies: `iso_fortran_env`
- summary: 吸着(insulator)モデルのメインループを実行し、電荷堆積と統計更新を行うモジュール。

### `bem_simulator_io`

- kind: `submodule`
- path: `src/runtime/simulator/bem_simulator_io.f90`
- group: `src/runtime/simulator`
- parent: `bem_simulator`
- internal dependencies: `bem_simulator`
- external dependencies: なし
- summary: `bem_simulator` の進捗表示と履歴出力を実装する submodule。

### `bem_simulator_loop`

- kind: `submodule`
- path: `src/runtime/simulator/bem_simulator_loop.f90`
- group: `src/runtime/simulator`
- parent: `bem_simulator`
- internal dependencies: `bem_simulator`, `bem_output_writer`, `bem_performance_profile`
- external dependencies: なし
- summary: `bem_simulator` の主ループと粒子処理計算を実装する submodule。

### `bem_simulator_stats`

- kind: `submodule`
- path: `src/runtime/simulator/bem_simulator_stats.f90`
- group: `src/runtime/simulator`
- parent: `bem_simulator`
- internal dependencies: `bem_simulator`
- external dependencies: なし
- summary: `bem_simulator` のバッチ集計・統計更新処理を実装する submodule。


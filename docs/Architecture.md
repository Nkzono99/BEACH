title: 開発者向けアーキテクチャ

Lang: [日本語](Architecture.md) | [English](Architecture.en.md)

# 開発者向けアーキテクチャ

BEACH の Fortran 実装を初めて変更する開発者が、実行入口から変更対象と直接 test まで移動するための概要です。
通常の物理・数値サイクルは[BEACH の計算サイクル](Algorithms.html)、build と test の選び方は
[開発ワークフロー](Workflow.html)を参照してください。このページは全 module の一覧を再掲せず、runtime の
制御フロー、主要 state の所有者、subsystem の境界だけを扱います。

## 実行フローを追う

```mermaid
flowchart TD
    cli["app/main.f90\nCLI / MPI 初期化"]
    config["config parser / runtime\napp_config と mesh を構築"]
    restart["restart または初期 state\nq_elem / stats / residual / ledger"]
    loop["run_absorption_insulator\naccepted batch / trial loop"]
    field["electrostatic snapshot\ncommit 済み q_elem から refresh"]
    inject["source plan / injection\ntrial の particles_soa を生成"]
    step["particle step / events\nBoris → 最初の mesh / box event"]
    commit["closure / MPI reduce / commit\ndq を q_elem へ一度反映"]
    record["stats / history / checkpoint\naccepted state だけを記録"]
    final["main の最終出力\nsummary / CSV / checkpoint"]

    cli --> config --> restart --> loop
    loop --> field --> inject --> step --> commit --> record
    record -->|次の accepted batch| field
    record -->|batch_count 到達| final
```

1. [`app/main.f90`](../app/main.f90) は CLI、MPI、performance profile を初期化し、設定 path を解決します。
   `load_or_init_run_state` が設定を読み、mesh と初期 state または restart state を用意します。
2. [`bem_app_config_parser.f90`](../src/config/app_config_parser/bem_app_config_parser.f90) は TOML を
   `app_config` へ読みます。派生値と組合せ制約は parser の finalize / validate submodule で確定します。
   [`bem_app_config_mesh_runtime.f90`](../src/config/bem_app_config_mesh_runtime.f90) が template または OBJ から
   `mesh_type` を構築します。
3. `main` は [`run_absorption_insulator`](../src/runtime/simulator/bem_simulator.f90) を呼びます。interface は
   `bem_simulator.f90`、主 loop は [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)、
   粒子生成・追跡は [`bem_simulator_particles.f90`](../src/runtime/simulator/bem_simulator_particles.f90)、
   電荷反映・電流補正・台帳は [`bem_simulator_charge.f90`](../src/runtime/simulator/bem_simulator_charge.f90) が担当します。
   統計と履歴は `bem_simulator_stats.f90` と `bem_simulator_io.f90` の submodule に分かれます。
4. simulator は commit 済み `mesh%q_elem` から
   [`electrostatic_snapshot_type`](../src/physics/bem_electrostatic_snapshot.f90) を refresh します。
   同じ trial の粒子追跡中はこの snapshot を固定し、accepted commit の電荷は次の batch の refresh で初めて場へ入ります。
5. `build_particle_source_plan` と `prepare_batch_state` が、source 設定と `batch_duration` から trial 用の
   `particles_soa` を作ります。設定からの粒子構築は
   [`bem_app_config_particle_runtime.f90`](../src/config/bem_app_config_particle_runtime.f90)、分布 sampling は
   `src/particles/` が担当します。
6. `process_particle_batch` は [`bem_particle_stepper.f90`](../src/runtime/simulator/bem_particle_stepper.f90) を通して
   予測中点場、Boris 更新、候補軌道を作ります。`bem_collision.f90` と `bem_boundary.f90` が最初の mesh hit または
   box event を確定し、吸収、escape、reflect、periodic wrap 後の再積分へ分岐します。
7. hit 電荷と放出反作用電荷は [`simulator_batch_workspace_type`](../src/runtime/simulator/bem_simulator_workspace.f90) の
   thread-local `dq` に蓄積します。surface / current closure と MPI reduce が成功した accepted trial だけを
   `commit_batch_charge` が `mesh%q_elem` へ一度加え、必要なら conductor 電荷を再配分します。
8. commit 後に `sim_stats` と [`charge_ledger_type`](../src/runtime/coupling/bem_charge_ledger.f90) を更新し、履歴と定期
   checkpoint を書きます。`main` は最終的に [`bem_output_writer.f90`](../src/runtime/bem_output_writer.f90) から
   summary、CSV、最終 checkpoint を公開します。

matching-plane の初期化、試行状態、固定点判定、継続解の確定は
[`bem_matching_plane_coupling.f90`](../src/runtime/simulator/bem_matching_plane_coupling.f90) が管理します。
陰的な面平均電荷の更新と根探索は
[`bem_matching_plane_implicit.f90`](../src/physics/sheath/bem_matching_plane_implicit.f90) に分かれています。
主ループはこれらの内部状態を直接変更せず、粒子の再試行とバッチ全体の受理・電荷反映を進めます。

adaptive batch-duration は手順 5--7 を同じ batch 開始 state から再生します。matching-plane 固定点反復では、
応答と snapshot の gauge も更新して手順 4--7 を再生します。棄却 trial の候補電荷、粒子 outcome、RNG、
macro 粒子端数、outer state は accepted state にしません。受理・rollback の
詳細は[`batch_duration` の理論](BatchDurationTheory.html)と
[matching-plane 準定常連成](MatchingPlaneCoupling.html)を参照してください。

## 主要 state の所有者を確認する

| State | 所有者と lifetime | 更新規則 |
| --- | --- | --- |
| `app_config` | `main` が構築し、run 全体で保持 | parser / runtime resolution 後は simulator へ read-only で渡す |
| `mesh_type` geometry | `main` が構築し、run 全体で保持 | 頂点、panel geometry、collision index は原則不変 |
| `mesh%q_elem` | `mesh_type` が持つ canonical な表面電荷 | run 前に初期化または復元し、その後は accepted trial の `commit_batch_charge` だけが更新する |
| `electrostatic_snapshot_type` | simulator が run 中に保持する派生 cache | commit 済み `q_elem` から refresh し、粒子追跡中は固定する。正本の電荷ではない |
| `particles_soa` | 1 trial の粒子 batch | source から生成し、吸収・escape・上限到達まで追跡した後に破棄する |
| `simulator_batch_workspace_type` | simulator が再利用する作業領域 | thread-local `dq`、候補電荷、outcome flag を保持する。commit 前は canonical state ではない |
| `injection_state` と RNG | accepted batch 間で継続し、restart で復元 | macro 粒子端数と乱数列を継続する。trial 棄却時は batch 開始 state へ戻す |
| `sim_stats` | `main` と simulator が保持する累積統計 | accepted trial だけを加算し、summary / checkpoint へ保存する |
| `charge_ledger_type` | run 全体の signed charge stock / flux | accepted batch の移送だけを累積し、保存残差と未解決量を別々に保持する |
| output / checkpoint files | writer が committed state から作る serialized copy | file 自体を runtime state の所有者にせず、loader が契約検査後に復元する |

この区別により、変更時には「候補値を作る場所」と「accepted state を確定する場所」を分けて確認できます。
trial-local 配列を更新しただけで、統計、ledger、履歴、checkpoint まで更新したとみなしてはいけません。

## 公開入口から責務の担当へ進む

公開入口は呼び出し順とデータの受け渡しを管理し、個別の形式や物理分野の処理を次の実装へ委譲します。

| 処理 | 公開入口 | 実装の担当 |
| --- | --- | --- |
| 場の評価 | `bem_field_solver.f90` | `_config` は設定解決、`_tree` は treecode の木とモーメント、`_fmm` は FMM core のパネル幾何・電荷状態、`_eval` は評価方式の切り替え |
| 外部シース応答 | `bem_matching_plane_response_provider.f90` | 親 module はモデル評価とフィードバックの契約、`_mpi` は設定からの初期化・rank 間の合意・root の評価結果の配信 |
| 応答テーブル | `bem_matching_plane_response.f90` | 親 module は不変 snapshot の共有と補間、`_io` は CSV 読み込み・格子検証 |
| Fortran の結果出力 | `bem_output_writer.f90` | `_history` submodule は履歴の生成・追記、`_summary` はサマリ、`_files` はメッシュ・電荷・台帳 CSV |
| チェックポイントの再開 | `bem_restart.f90` | `_contract` は再開条件の検証、`_records` は統計・電荷・台帳の読み込み、`_injection` は乱数・マクロ粒子端数の保存と復元 |
| 設定から粒子を生成 | `bem_app_config_particle_runtime.f90` | 親 module は粒子源計画、`_batch` は MPI 配分とバッチ構築、`_sampling` は種別ごとのサンプリングと注入速度補正 |
| Python の設定処理 | [`beach/config/core.py`](../beach/config/core.py) | [`_authoring.py`](../beach/config/_authoring.py) は空間指定の展開、[`_runtime_validation.py`](../beach/config/_runtime_validation.py) は場・粒子・表面電流・メッシュの検証を順に呼ぶ |
| Python の結果読み込み | [`beach/fortran_results/io.py`](../beach/fortran_results/io.py) | 基本のメッシュ・電荷を読み、[`_matching_plane_io.py`](../beach/fortran_results/_matching_plane_io.py) と [`_field_reconstruction_io.py`](../beach/fortran_results/_field_reconstruction_io.py) に連成状態・場の再構築メタデータを委譲する |

FMM の木構造・相互作用リストは `field_solver_type%fmm_core_plan`、電荷から計算する作業状態は
`%fmm_core_state` が保持します。旧 FMM の複製 view と未使用の局所展開配列は削除しました。
診断時は core plan / state を参照し、treecode 用配列から FMM の状態を読まないでください。
`init` は既存の FMM 作業状態を解放してから再構築します。たとえば、同じ `solver` に対して
`call solver%init(mesh, sim)` を再度呼び、`sim%field_solver` を切り替えられます。
`refresh(mesh)` は電荷更新に用い、FMM では空メッシュまたは要素数の変更にも対応します。
同じ要素数で頂点座標を変更した場合は `init` で幾何を再構築してください。

### シースと外部応答の担当

`src/physics/sheath/` には、固定電流を決める定常問題、帯電状態に応答する matching-plane 問題、
その応答表を作るオフラインツールが入っています。定常問題は零電流条件を解き、matching-plane 問題は
与えられた電束密度と粒子流束から外部応答を返すため、両者の根探索は同じ問題ではありません。

| ファイル（`src/physics/sheath/` 以下） | 担当 |
| --- | --- |
| `bem_sheath_model_core.f90` | Zhao モデルの密度・電荷密度・残差式と、定常解の非線形方程式 |
| `bem_surface_current_model.f90` | 設定と定常シース解を、粒子種別の吸収・放出・流入電流へ変換 |
| `bem_surface_closure_contract.f90` | simulator が受け取る電流・境界条件のデータ型。モデル固有の解法は持たない |
| `bem_matching_plane_zhao.f90` | 電束密度と流束を入力に、A/B/C 分岐の根探索、物理解の選択、応答評価、継続解の追跡 |
| `bem_matching_plane_implicit.f90` | 硬い面平均帯電を後退 Euler で解く。応答モデルを反復評価し、根の挟み込みと必要な時間分割を行う |
| `bem_matching_plane_response_provider.f90` | table / online Zhao の共通入口と、フィードバックの範囲・尺度・収束判定 |
| `bem_matching_plane_response_provider_mpi.f90` | provider の設定解決、rank 間の設定・query 合意、root で求めた応答の配信 |
| `bem_matching_plane_response.f90` | 応答表の保持、path ごとの snapshot cache、5 次元補間、補間軸の取得 |
| `bem_matching_plane_response_io.f90` | CSV の構文・単位付き列名・格子の欠損や重複を検証して応答表を構築 |
| `bem_matching_plane_response_mpi.f90` | root が読んだ補間軸・値・高度・出典 path を全 rank に配信 |
| `bem_matching_plane_response_generator.f90` | online Zhao を格子上で評価して、実行用の応答 CSV を作るオフラインツール |
| `bem_matching_plane_zhao_atlas.f90` | A/B/C 分岐の成立範囲や失敗理由を調べるオフライン診断ツール |

通常の連成では runtime の `bem_matching_plane_coupling` が provider を使い、陰解法を選んだ場合だけ
`bem_matching_plane_implicit` を介します。generator と atlas は通常の batch loop には入りません。
`src/physics/bem_surface_models*.f90` は物体側の電荷再配分・導体条件などを担当し、外部シース応答とは別です。

読みやすさの課題は、`bem_matching_plane_zhao` に数式、Newton 法、候補解の選別、continuation が集まっていることと、
generator / atlas に query CSV の読み取り処理が重複していることです。これらは次の分離候補ですが、
分岐の成立条件や負の電場二乗を拒否する検査は物理モデルに必要です。

応答表のハッシュ照合は廃止し、root の読込結果を配信します。補間軸が必要なコードは
`call table%get_axis_data(axis_sizes, axis_values, matching_plane_z_m, status, message)` を使います。
再開時はメッシュ識別子だけを照合し、モデル・粒子種・応答内容の fingerprint は生成しません。
FMM の演算子キャッシュには、別条件の演算子を再利用しないための識別子を残しています。

### 粒子配列を直接生成する

`bem_particles` の `allocate_particles(pcls, n)` が非負の粒子数 `n` に一致する配列の確保を担当します。
位置・速度・電荷・質量は呼び出し元が埋めます。重みは 1、粒子種 ID は 0、放出元要素は -1、
`alive` は true で初期化されます。`n=0` でも長さゼロの配列を確保し、既存の粒子群を置き換えます。
たとえば、静止粒子を直接生成する場合は次のように書けます。

```fortran
use bem_kinds, only: dp, i32
use bem_types, only: particles_soa
use bem_particles, only: allocate_particles

type(particles_soa) :: pcls

call allocate_particles(pcls, 100_i32)
pcls%x = 0.0_dp
pcls%v = 0.0_dp
pcls%q = -1.602176634e-19_dp
pcls%m = 9.1093837015e-31_dp
```

既存の配列を検証してコピーする用途は引き続き `init_particles` が担当します。
注入時は分布サンプラーが生成先へ直接書き込み、バッチ構築では種別ごとの必要数だけを一時保持し、
従来の種別交互順で完成した SoA へ詰めます。乱数の消費順と放出元要素の対応を維持します。

mesh の幾何更新は `bem_panel_geometry` が計算した重心・法線・面積を利用します。
`fill_panel_quadrature(panel, position, weight)` は確保済みの `position(3,7)` と `weight(7)` に
積分点と重みを書き込みます。mesh はこの入口を使い、三角形ごとの一時配列の確保・コピーを省きます。
積分計画を新しく構築する既存の `build_panel_quadrature` も、この計算を共有します。
mesh 更新では次のように要素 `i` の配列へ書き込みます。

```fortran
call fill_panel_quadrature(panel, mesh%panel_quad_position(:, :, i), mesh%panel_quad_weight(:, i))
```

積分対象から除外する細長い三角形も、衝突判定用の幾何情報は保持します。

## Subsystem から実装と test へ移動する

表の test は変更直後に使う直接 test です。必要な累積 gate は[開発ワークフロー](Workflow.html#変更からテストを選ぶ)で
選びます。

| Subsystem | 主な source | 直接 test | 正本・解説 |
| --- | --- | --- | --- |
| CLI、config、runtime resolution | `app/main.f90`、`src/config/` | [`test_app_config_parser.f90`](../tests/fortran/test_app_config_parser.f90)、[`test_physics_config_types.f90`](../tests/fortran/test_physics_config_types.f90)、`tests/python/test_config_schema.py`、`test_config_cli.py` | [設定を編集する](Configuration.html)、[設定パラメータ](Parameters.html) |
| mesh、template、OBJ、panel geometry | `src/mesh/`、`src/physics/panel/` | [`test_templates_importers_runtime.f90`](../tests/fortran/test_templates_importers_runtime.f90)、[`test_panel_geometry_near.f90`](../tests/fortran/test_panel_geometry_near.f90)、[`test_panel_kernel.f90`](../tests/fortran/test_panel_kernel.f90) | [設定レシピ](ConfigurationRecipes.html)、[Direct](DirectSolver.html) |
| batch orchestration | `src/runtime/simulator/bem_simulator*.f90` | [`test_simulator.f90`](../tests/fortran/test_simulator.f90)、[`test_dynamics_basic.f90`](../tests/fortran/test_dynamics_basic.f90) | [`SPEC.md`](../SPEC.md)、[BEACH の計算サイクル](Algorithms.html) |
| field snapshot、Direct / Treecode / FMM、periodic2 | `bem_electrostatic_snapshot*.f90`、`src/physics/field_solver/`、`src/physics/periodic_zero_mode/` | [`test_electrostatic_snapshot.f90`](../tests/fortran/test_electrostatic_snapshot.f90)、[`test_dynamics_field_solver.f90`](../tests/fortran/test_dynamics_field_solver.f90)、`test_dynamics_fmm`、`test_periodic_zero_mode`、`test_periodic2_cached_snapshot` | [場の評価](FieldSolvers.html)、[FMM](FMM.html)、[periodic2 静電場](PeriodicElectrostatics.html) |
| particle source と injection | `bem_app_config_particle_runtime.f90`、`src/particles/` | [`test_injection_sampling.f90`](../tests/fortran/test_injection_sampling.f90)、[`test_reservoir_injection.f90`](../tests/fortran/test_reservoir_injection.f90)、[`test_external_field_velocity_grid.f90`](../tests/fortran/test_external_field_velocity_grid.f90) | [粒子をどこから入れるか](ParticleSourcesBoundaries.html)、[境界から粒子を流入させる](ReservoirInjection.html)、[光電子放出](PhotoelectronEmission.html) |
| Boris、collision、box event | `bem_particle_stepper.f90`、`bem_pusher.f90`、`bem_collision.f90`、`bem_boundary.f90` | [`test_particle_stepper.f90`](../tests/fortran/test_particle_stepper.f90)、[`test_boundary.f90`](../tests/fortran/test_boundary.f90)、`test_dynamics_basic` | [粒子更新](ParticleTrackingCollision.html)、[Boris](BorisPusher.html)、[粒子 event](ParticleEvents.html) |
| surface charge、closure、ledger | `bem_surface_models*.f90`、`src/physics/sheath/`、`bem_matching_plane_coupling.f90`、`bem_simulator_charge.f90`、`bem_charge_ledger.f90` | [`test_surface_models.f90`](../tests/fortran/test_surface_models.f90)、[`test_surface_current_model.f90`](../tests/fortran/test_surface_current_model.f90)、[`test_charge_ledger.f90`](../tests/fortran/test_charge_ledger.f90)、`test_matching_plane_simulator` | [表面はどう帯電するか](SurfaceModels.html)、[表面電荷更新の数値仕様](SurfaceChargeNumerics.html)、[matching-plane 連成](MatchingPlaneCoupling.html) |
| stats、output、checkpoint、restart | `bem_simulator_stats.f90`、`bem_simulator_io.f90`、`bem_output_writer.f90`、`bem_periodic_checkpoint.f90`、`bem_restart.f90` | [`test_output_writer_io.f90`](../tests/fortran/test_output_writer_io.f90)、[`test_output_writer_potential.f90`](../tests/fortran/test_output_writer_potential.f90)、[`test_restart.f90`](../tests/fortran/test_restart.f90) | [出力ガイド](OutputGuide.html)、[実行と再開](Execution.html)、`SPEC.md` の出力・再開契約 |
| Python reader、解析、可視化 | `beach/` | `tests/python/test_fortran_results.py`、対応する CLI / analysis test | [後処理チュートリアル](PostprocessTutorial.html)、[Python API](PythonPostprocessAPI.html) |

module 名や `use` 依存を検索するときは、自動生成した
[Fortran 依存関係マップ](FortranDependencyMap.html)と[Fortran API](https://nkzono99.github.io/BEACH/fortran/)を使います。
依存関係マップは source inventory であり、runtime の呼出順、state ownership、behavioral contract の正本ではありません。

## 正本の責務を区別する

| 情報 | 正本 | guide / reference の責務 |
| --- | --- | --- |
| 現行 simulation behavior と model scope | Fortran 実装と [`SPEC.md`](../SPEC.md) | model / numerical-method page は理由、式、適用範囲、検証方法を説明する |
| 公開 TOML の table、key、型、構造制約 | `schemas/beach.schema.json`。派生値と意味的な組合せは Fortran parser / validator | `Parameters.md` / `.en.md` は検索可能な人間向け reference、Configuration は編集手順を示す |
| output file の生成条件 | `schemas/beach.output-manifest.json` と Fortran writer | OutputGuide は column の意味、確認順、restart での役割を説明する |
| checkpoint compatibility | checkpoint contract、mesh identity、writer / loader、`SPEC.md` | Execution は安全な再開手順を示す |
| test target と tier | `fpm.toml` と `Makefile` | Workflow は変更範囲から実行すべき target へ案内する |
| site の page inventory と sidebar | `docs-site/navigation.json` | `docs/*.md` と `.en.md` が編集する source で、`docs-site/src/content/docs/` は生成物 |
| module / procedure API と依存 | Fortran source、生成した FORD API、FortranDependencyMap | Architecture は人が読む実行フローと subsystem 境界だけを維持する |

tutorial、task guide、example は正本の契約を短く適用する入口です。そこへ全 parameter や全分岐を複製せず、
該当する reference または specification へリンクします。behavior、config、output を変更するときの同期対象は
[公開契約を変更するとき](Workflow.html#公開契約を変更するとき)で確認してください。

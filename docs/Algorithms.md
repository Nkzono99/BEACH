title: BEACH アルゴリズム概要

Lang: [日本語](Algorithms.md) | [English](Algorithms.en.md)

# BEACH アルゴリズム概要

この文書は、BEACH の現行 Fortran 実装に基づいて、数値アルゴリズムと実行順序を説明します。
BEACH は格子 PIC ではなく、三角形境界要素上の電荷を source とする BEM 風の Coulomb 場評価と、
テスト粒子追跡を batch 単位で結合した表面帯電シミュレータです。

実装へのリンクは、現在のファイルと主な symbol を示します。行番号が後続変更でずれても、リンク先ファイルと symbol 名を優先してください。

---

## アルゴリズム文書の構成

| 文書 | 内容 |
|---|---|
| [アルゴリズム概要](Algorithms.html) | BEACH の計算モデル、初期化、batch ループ |
| [場ソルバーと境界条件](FieldSolvers.html) | Coulomb 場、direct/treecode/FMM 切替、periodic2 場境界 |
| [粒子追跡と表面電荷蓄積](ParticleChargeLoop.html) | 粒子生成、Boris pusher、衝突、電荷堆積、統計・再開 |
| [Coulomb FMM コア詳細](FMMCore.html) | FMM core API、tree 構築、M2L、periodic2 Ewald/oracle |
| [`batch_duration` の安定性](BatchDurationStability.html) | batch 時間幅、定常値、線形安定性、Monte Carlo ノイズ |


## 1. BEACH の計算モデル

**Source**:
[`bem_simulator`](../src/runtime/simulator/bem_simulator.f90),
[`bem_simulator_loop`](../src/runtime/simulator/bem_simulator_loop.f90),
[`bem_particle_stepper`](../src/runtime/simulator/bem_particle_stepper.f90),
[`bem_field_solver`](../src/physics/field_solver/bem_field_solver.f90),
[`bem_injection`](../src/particles/bem_injection.f90)

BEACH の主状態は、三角形メッシュ要素の電荷 `q_elem(i)` と、batch ごとに生成される粒子群です。
メッシュ要素は固定形状で、粒子は表面に衝突した時点で吸収されます。吸収された粒子の電荷は命中要素に蓄積し、次 batch の場計算へ反映されます。

### 1.1 状態変数

主要データ型は [`bem_types`](../src/core/bem_types.f90) に定義されています。

| 型 | 主な内容 | 役割 |
| --- | --- | --- |
| `sim_config` | `dt`, `batch_count`, `max_step`, `field_solver`, `field_bc_mode`, `box_min/max`, `bc_low/high` | 時間発展・場計算・境界条件 |
| `mesh_type` | `v0/v1/v2`, `centers`, `normals`, `bb_min/max`, `q_elem`, `elem_surface_model` | 三角形境界要素と電荷 |
| `particles_soa` | `x`, `v`, `q`, `m`, `w`, `species_id`, `alive` | batch 内粒子の SoA 表現 |
| `injection_state` | `macro_residual(:)` | `reservoir_face` の端数粒子数を batch 間・再開間で保持 |
| `sim_stats` | `processed_particles`, `absorbed`, `escaped`, `survived_max_step`, `batches`, `last_rel_change` | 実行統計 |

### 1.2 物理近似

各三角形要素は、重心 `c_i` にある点電荷 `q_i` として場に寄与します。要素面積上の連続電荷分布を積分する厳密 BEM ではなく、現行実装は重心点電荷近似です。

softening 付き Coulomb kernel:

$$
G_\epsilon(\mathbf{r}) = \frac{1}{\sqrt{\lVert\mathbf{r}\rVert^2 + \epsilon^2}}
$$

電位:

$$
\phi(\mathbf{x}) =
k_\mathrm{c}
\sum_j q_j G_\epsilon(\mathbf{x} - \mathbf{c}_j)
$$

電場:

$$
\mathbf{E}(\mathbf{x}) =
k_\mathrm{c}
\sum_j q_j
\frac{\mathbf{x} - \mathbf{c}_j}
{\left(\lVert\mathbf{x} - \mathbf{c}_j\rVert^2 + \epsilon^2\right)^{3/2}}
$$

ここで `k_coulomb` は [`bem_constants`](../src/core/bem_constants.f90) の Coulomb 定数です。
粒子が実際に受ける電場は、境界要素電荷による場に一様外部電場 `sim.e0` を加えたものです。

### 1.3 実行単位

BEACH は `sim.batch_count` まで batch を進めます。各 batch は次の意味を持ちます。

1. 現在の表面電荷に基づき、粒子種ごとの粒子群を生成する。
2. 現在の `q_elem` で field solver を更新する。
3. 各粒子を最大 `sim.max_step` ステップまで追跡する。
4. 衝突した粒子の電荷を要素別差分 `dq` に集計する。
5. `dq` を `q_elem` に加算し、必要なら表面モデルの緩和をかける。
6. 統計・履歴を更新する。

`sim.tol_rel` は出力・監視用の値であり、現行 Fortran 実装では早期停止条件ではありません。

---

## 2. 実行エントリーポイントと初期化

**Source**:
[`app/main.f90`](../app/main.f90),
[`bem_app_config_runtime`](../src/config/bem_app_config_runtime.f90),
[`bem_mesh`](../src/mesh/bem_mesh.f90),
[`bem_restart`](../src/runtime/bem_restart.f90)

### 2.1 main program の順序

`app/main.f90` は CLI エントリーポイントです。概略順序は次の通りです。

| 順序 | 処理 | 主な実装 |
| --- | --- | --- |
| 1 | `--version` など設定読込前に完結する CLI option を処理 | `handle_early_cli` |
| 2 | MPI と performance profile を初期化 | `mpi_initialize`, `perf_configure_from_env` |
| 3 | 設定読込、メッシュ構築、restart 読込または乱数 seed 初期化 | `load_or_init_run_state` |
| 4 | 履歴 CSV を open | `open_history_writer`, `open_potential_history_writer` |
| 5 | batch simulation を実行 | `run_absorption_insulator` |
| 6 | summary と最終 CSV を出力 | `print_run_summary`, `write_result_files` |
| 7 | RNG state と macro residual を checkpoint として出力 | `write_rng_state_file`, `write_macro_residuals_file` |
| 8 | performance profile を出力し MPI を終了 | `perf_write_outputs`, `mpi_shutdown` |

設定ファイルは、明示引数があればその path、なければカレントディレクトリの `beach.toml` を使います。
設定ファイルがない場合は `default_app_config` の既定値で走ります。

### 2.2 メッシュ構築

メッシュは `mesh.mode` に応じて template または OBJ から作られます。template の場合は
[`build_template_mesh`](../src/config/bem_app_config_runtime.f90) が次を行います。

1. `mesh.templates` を順に走査する。
2. `enabled=false` の template は飛ばす。
3. `kind` に応じて `make_plane`, `make_box`, `make_cylinder`, `make_sphere` などへ dispatch する。
4. template ごとに `mesh_id` を割り当てる。
5. `surface_model` と `epsilon_r` を要素配列へ展開する。
6. 全 template の三角形配列を連結し、`init_mesh` へ渡す。

`init_mesh` は頂点配列から次を前計算します。

- 要素重心 `centers(:, i)`
- 要素法線 `normals(:, i)`
- 要素 AABB `bb_min/max(:, i)`
- 代表長 `h_elem(i) = sqrt(area_i)`
- 初期電荷 `q_elem(i)`
- collision grid

### 2.3 periodic2 collision mesh

`sim.field_bc_mode="periodic2"` のとき、`prepare_periodic2_collision_mesh` は周期 2 軸に沿って各三角形を canonical な位置へ平行移動します。
これは場計算の周期画像和とは別で、衝突判定用の primitive cell メッシュを安定化するための処理です。
要素 index は base element のまま維持されるため、periodic image に命中しても電荷は base element に集約されます。

### 2.4 restart

`output.resume=true` のとき、`load_restart_checkpoint` は以下を読みます。

- `summary.txt`: 完了済み batch と統計
- `charges.csv`: 要素電荷
- `rng_state.txt` または `rng_state_rankNNNNN.txt`: 乱数状態
- `macro_residuals.csv` または rank 別 residual: reservoir の端数粒子数
- `charge_ledger.csv`: schema v2 の累積 signed charge ledger

`sim.batch_count` は累積到達 batch 数です。checkpoint が `batches=100` で `batch_count=150` なら、実行するのは 50 batch だけです。
schema v2 は model、ordered mesh、ordered species の fingerprint を保存し、状態を変更する前に照合します。旧 schema は実装済み legacy point-source model に限って受理します。

粒子積分の contract は `particle_time_centering="same_time_midpoint_boris"` です。pure E / pure B、time reversal、smooth field の二次収束、batch restart continuity の回帰試験を、この contract の変更検出に使います。

---

## 3. batch ループ

**Source**:
[`run_absorption_insulator`](../src/runtime/simulator/bem_simulator_loop.f90),
[`prepare_batch_state`](../src/runtime/simulator/bem_simulator_loop.f90),
[`process_particle_batch`](../src/runtime/simulator/bem_simulator_loop.f90),
[`commit_batch_charge`](../src/runtime/simulator/bem_simulator_loop.f90)

### 3.1 loop skeleton

`run_absorption_insulator` は、初期統計 `initial_stats` を受け取った場合はその `batches` から再開します。

```text
final_batch_idx = sim.batch_count
batch_count_this_run = final_batch_idx - stats.batches

field_solver.init(mesh, sim)

for local_batch_idx = 1..batch_count_this_run:
    prepare_batch_state(...)
    field_solver.refresh(mesh)
    process_particle_batch(...)
    commit_batch_charge(...)
    count_batch_outcomes(...)
    MPI allreduce statistics
    accumulate_batch_stats(...)
    write charge/potential history when stride matches
```

### 3.2 batch state

`prepare_batch_state` は次を準備します。

- 今 batch の番号 `batch_idx = stats%batches + 1`
- `init_particle_batch_from_config` による粒子群
- `dq_thread(nelem, nth)`: OpenMP thread ごとの電荷差分
- `photo_emission_dq(nelem)`: `photo_raycast` の放出元逆符号電荷
- `escaped_boundary_flag(:)` と `absorbed_flag(:)`

`dq_thread` を thread 別に分けることで、衝突時の要素電荷加算を atomic なしで集計できます。
`photo_raycast` の衝突照会が不完全な場合、sampler は OpenMP 内で停止せず、最小の ray / bounce と status を
`init_particle_batch_from_config` から `prepare_batch_state` へ返します。main loop は field refresh や電荷処理の前に
全 rank から最小 rank の species / ray / bounce / status を選び、全 rank が同じ診断で停止します。
失敗 rank の未完成 particle 配列と `photo_emission_dq` は後続処理で使用しません。

### 3.3 particle processing

`process_particle_batch` は粒子ごとに最大 `sim.max_step` 回の時間発展を行います。

| 順序 | 処理 |
| --- | --- |
| 1 | 同一時刻の現在位置 `x0` と速度 `v0` を読む |
| 2 | `build_particle_step_candidate` で予測中点 `x_mid = x0 + 0.5*v0*dt` を作る |
| 3 | `field_solver%eval_e(mesh, x_mid, e_mid)` で境界要素電荷による電場を評価し、一様外部電場 `sim.e0` を1回加える |
| 4 | 中点場を使う `boris_push` で候補速度 `v1` と台形則による候補位置 `x1` を計算 |
| 5 | `x0 -> x1` を1回collision queryし、`x1` がbox内部ならその結果を確定する |
| 6 | box crossingがあればmesh hit parameterと最初のface parameterを比較し、tieではmeshを優先する |
| 7 | openはevent点で終了し、reflect/periodicは残り時間を一度だけ再積分してmesh hitを調べる |
| 8 | hitなら`q * w`を堆積して吸収し、2回目のbox eventならstateをcommitせずfail closedにする |
| 9 | 生存していれば同時刻の`x`と`v`を更新して次stepへ進む |

`build_particle_step_candidate` は場ソルバを変更せず、予測中点で空間電場を1回だけ評価します。
`boris_update_velocity(v, q, m, dt, e, b, v_new)` は、電場の half kick、磁場回転、電場の half kick からなる
速度更新を提供する public pure procedure です。既存の public call
`boris_push(x, v, q, m, dt, e, b, x_new, v_new)` は署名を変えず、速度計算をこの procedure に委譲して
`x_new = x + 0.5*(v + v_new)*dt` で位置を更新します。入出力の位置と速度は同一時刻の状態であり、
予測中点の空間電場評価と台形位置更新により candidate kinematics は二次精度です。

production loopはcandidateとmesh queryを先に作り、候補終点がstrictなbox内部なら、追加のevent geometryなしに場評価1回・collision query 1回で
完了します。crossing時だけ `find_first_boundary_event` と `apply_escape_reflect_periodic_event` を使い、corner/edgeの
同時faceを一括処理します。reflect/periodic remainderがさらにfaceへ達し、それ以前にmesh hitがなければ、
任意回数のloopへ入らず `particle_step_multiple_box_events` を返します。既存 `apply_box_boundary` はphoto ray用に維持します。
periodic2のfull-chord queryがbox外区間でrange limitに達した場合は、最初のbox eventまでに制限したqueryへfallbackします。

`potential_barrier` は単一open faceに限り旧scalar energy式をevent位置と補間速度で評価します。複数open faceは
一般化せずfail closedであり、shared potential/gaugeに基づく物理モデルは後段の対象です。

`BEACH_WARN_LONG_PARTICLE_STEPS` を正整数で設定すると、長く生き残る粒子の診断出力を一定 step ごとに出します。

collision status は `collision_query_ok=0`、`collision_query_image_limit=1`、
`collision_query_index_range=2` です。`status` を省略した public call は不完全な照会を miss として扱わず、
fail closed で停止します。main particle loop はcollisionとboundary eventのfailureを同じenvelopeで集約し、OpenMP regionを抜けてから
MPI全rankで最小rankのparticle/stepとその`dt/x/v`を共有し、同一messageで停止します。

### 3.4 charge commit

`commit_batch_charge` は thread 別差分と photo emission 差分を統合します。

$$
\Delta q_i =
\sum_{\mathrm{thread}} \Delta q_{i,\mathrm{thread}}
{}+ \Delta q_{i,\mathrm{photo}}
$$

MPI 実行時は `mpi_allreduce_sum_real_dp_array` により rank 間で `dq` を和にします。
その後:

$$
q_i^{new} = q_i^{old} + \Delta q_i
$$

を適用し、表面モデル緩和を実行します。最後に実際に変化した電荷量で監視値を計算します。

$$
\operatorname{rel\_change} =
\frac{\lVert\mathbf{q}^{\mathrm{new}} - \mathbf{q}^{\mathrm{old}}\rVert_2}
{\max\left(\lVert\mathbf{q}^{\mathrm{new}}\rVert_2, q_\mathrm{floor}\right)}
$$

この値は `stats%last_rel_change` と履歴出力に使われます。

---

### Direct P0 triangle panel

`field.element_kernel="triangle_p0"` では `q_elem` を要素総電荷、`sigma=q_elem/area` を一定面密度として、辺対数項と signed solid angle による解析式で電位・電場を評価します。面上評価は `bem_panel_self_terms` が所有し、電位は連続、principal-value 法線場は両側極限の平均、真空側極限は `elem_vacuum_sign` で選びます。幾何、辺量、厳密一次・二次 moment、7点求積 plan は mesh 初期化時に固定し、batch 更新では電荷だけを変更します。

この経路は Phase 1 では free-space direct の correctness oracle です。treecode/FMM/periodic2 には黙って点電荷近似せず、初期化時に停止します。

### 2.8 periodic2 split reference

`panel_spectral_reference`は電位を`k!=0`のP0 panel Fourier和、厳密なtriangle-height zero mode、線形Debye outer profileへ分解します。zero modeは傾斜三角形の累積面積を区分二次式として前計算し、電場とその積分電位を`O(log N)`で評価します。interfaceでは面内格子上の非零モード電位・全電場を測り、1D outerへ切り替えられる減衰量かを検査します。これは小規模参照実装であり、production規模の周期演算子ではありません。

### 2.9 outer particle interface

z-highのbox eventは、face、event fraction、同時刻の位置・速度、remaining `dt`をtyped payloadとしてsimulatorへ返します。linear-Debye instant-return mapは法線エネルギーからescape/turningを判定し、turning軌道の解析flight timeで接線位置を進めます。return後だけ通常stepperへ戻してremaining `dt`を再積分します。interfaceを使わない通常候補では追加探索も追加field評価もありません。

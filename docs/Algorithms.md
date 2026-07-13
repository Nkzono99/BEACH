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

境界要素の source model は `field.element_kernel` で選びます。互換既定の `"point"` と、
三角形上の一定面密度を積分する `"triangle_p0"` の2系統があります。

#### Point source

`element_kernel="point"` では、各三角形要素は重心 `c_i` にある点電荷 `q_i` として場に寄与します。

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

#### P0 triangle source

`element_kernel="triangle_p0"` では `q_j` を要素総電荷、`A_j` を要素面積、
`sigma_j=q_j/A_j` を三角形 `T_j` 上の一定面密度として扱います。

$$
\phi(\mathbf{x}) =
k_\mathrm{c}\sum_j \frac{q_j}{A_j}
\int_{T_j}\frac{1}{\lVert\mathbf{x}-\mathbf{y}\rVert}\,dA_{\mathbf{y}}
$$

$$
\mathbf{E}(\mathbf{x}) =
k_\mathrm{c}\sum_j \frac{q_j}{A_j}
\int_{T_j}\frac{\mathbf{x}-\mathbf{y}}
{\lVert\mathbf{x}-\mathbf{y}\rVert^3}\,dA_{\mathbf{y}}
$$

direct はこの panel 積分を解析式で評価します。FMM は近傍を同じ解析 panel kernel、遠方 P2M を
三角形上の monomial の厳密面積平均で評価します。`triangle_p0` は `sim.softening=0` を必須とし、
softened point source へ fallback しません。

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

### 2.3 P0 triangle panel field kernel

`field.element_kernel="triangle_p0"`では、`q_elem`は要素総電荷、`sigma=q_elem/area`は一定面密度です。

| 経路 | kernel | 実装contract |
| --- | --- | --- |
| direct | `triangle_p0_exact_direct` | edge-logとsigned solid angleによるfree-space correctness oracle |
| FMM | `triangle_p0_exact_p2m_near` | 解析panel near + 三角形monomialの厳密面積平均P2M |
| auto | `triangle_p0_exact_auto` | 対応solverから上記を選択 |

面上電位は連続です。principal-value法線場は両側traceの平均で、`elem_vacuum_sign`が真空側traceを選びます。
geometry、edge量、一次・二次moment、7点求積planはmesh構築時に固定し、batchでは電荷だけを更新します。
treecodeとpoint-source `m2l_root_oracle`へはfallbackせず、非対応設定として初期化時に停止します。

### 2.4 periodic2 split reference

`panel_spectral_reference`は小規模correctness用です。

| 成分 | 評価 |
| --- | --- |
| $k\ne0$ | P0 panel Fourier和 |
| $k=0$ | triangle-heightの区分二次累積面積を前計算し、$O(\log N)$で評価 |
| outer profile | 線形Debye 1D profile |
| interface gate | 面内格子で非零モード電位と全電場の減衰を測定 |

production規模の周期演算子ではありません。

### 2.5 outer particle interface

z-high eventはface、event fraction、同時刻位置・速度、remaining `dt`をtyped payloadで返します。

| return model | 外部領域の扱い |
| --- | --- |
| linear-Debye instant return | 法線エネルギーでescape/turningを分類し、解析flight timeで接線位置を進める |
| `electrostatic_3d_explicit_orbit` | unified 3D電位・電場を使い、固定刻みvelocity-Verletで追跡 |

return後だけremaining `dt`を通常stepperで再積分します。3D orbitはreturn/far-plane crossingをstep内補間し、
energy error、flight time、frozen-field ratioを検査します。上限違反はdiscardせずfail closedです。
interface無効時は通常fast pathに探索やfield評価を追加しません。

### 2.6 periodic2 collision mesh

`sim.field_bc_mode="periodic2"` のとき、`prepare_periodic2_collision_mesh` は周期 2 軸に沿って各三角形を canonical な位置へ平行移動します。
これは場計算の周期画像和とは別で、衝突判定用の primitive cell メッシュを安定化するための処理です。
要素 index は base element のまま維持されるため、periodic image に命中しても電荷は base element に集約されます。

### 2.7 restart

`output.resume=true` のとき、`load_restart_checkpoint` は以下を読みます。

- `summary.txt`: 完了済み batch と統計
- `charges.csv`: 要素電荷
- `rng_state.txt` または `rng_state_rankNNNNN.txt`: 乱数状態
- `macro_residuals.csv` または rank 別 residual: reservoir の端数粒子数
- `charge_ledger.csv`: schema v2/v3 の累積 signed charge ledger

`sim.batch_count` は累積到達 batch 数です。checkpoint が `batches=100` で `batch_count=150` なら、実行するのは 50 batch だけです。
schema v3 は model、ordered mesh、ordered species の fingerprint に加え、outer solver の `z, phi, E, rho`、status、反復数、residual、電流を保存し、状態を変更する前に照合します。schema v2 は migration input として受理し、不完全な outer profile は再解法します。それより古い schema は実装済み legacy point-source model に限って受理します。

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
    outer_coupler.refresh(snapshot, mesh, batch_idx)
    prepare_batch_state(...)
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

`kinetic_1d_profile_return`では、直前にMPI-global更新された`snapshot.outer`を
`init_particle_batch_from_config`へ渡します。これにより無限遠reservoirのaccessible fluxと
interface速度が、そのbatchのfieldと同一の`phi_interface-phi_infinity`から決まります。

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
| 2 | `build_particle_step_candidate` で予測中点 `x_mid = x0 + 0.5*v0*dt` を作り、場評価点だけをprimitive target boxへ写像する |
| 3 | `field_solver%eval_e(mesh, x_mid, e_mid)` で境界要素電荷による電場を評価し、一様外部電場 `sim.e0` を1回加える |
| 4 | 中点場を使う `boris_push` で候補速度 `v1` と台形則による候補位置 `x1` を計算 |
| 5 | `x0 -> x1` を1回collision queryし、`x1` がbox内部ならその結果を確定する |
| 6 | box crossingがあればmesh hit parameterと最初のface parameterを比較し、tieではmeshを優先する |
| 7 | openはevent点で終了し、reflect/periodicは残り時間を最大8回再積分して各chordのmesh hitを調べる |
| 8 | hitなら`q * w`を堆積して吸収し、9回目のbox eventならstateをcommitせずfail closedにする |
| 9 | 生存していれば同時刻の`x`と`v`を更新して次stepへ進む |

粒子積分の要点は次の通りです。

| 層 | contract |
| --- | --- |
| state | 入出力は $\mathbf{x}^n,\mathbf{v}^n$ と $\mathbf{x}^{n+1},\mathbf{v}^{n+1}$ の同時刻対 |
| field sample | 予測中点で空間電場を1回評価し、`sim.e0`を1回加える |
| velocity | `boris_update_velocity` が electric half kick、magnetic rotation、electric half kickを行う |
| position | `boris_push` が $\mathbf{x}^{n+1}=\mathbf{x}^n+(\mathbf{v}^n+\mathbf{v}^{n+1})\Delta t/2$ を適用 |
| accuracy | smoothな規定場で位置・速度とも二次精度 |

これはhalf-step速度を保持するleapfrog stateではありません。pure BのBoris回転は速度ノルムを保存しますが、
衝突、open境界、batchごとの自己無撞着場更新を含む全系の厳密な長時間エネルギー保存は保証しません。
式、event処理、保存性の範囲は[粒子時間積分](ParticleChargeLoop.html#8-粒子時間積分-boris速度更新と同時刻状態)を参照してください。

event処理は通常経路とcrossing経路を分けます。

| 経路 | 追加処理 |
| --- | --- |
| strictなbox内部 | なし。場評価1回、collision query 1回 |
| box crossing | mesh/faceの最早eventを比較し、corner/edgeの同時faceを一括処理 |
| reflect/periodic remainder | 最大8 eventまで再積分 |
| 9回目が必要 | `particle_step_multiple_box_events` でfail closed |

`apply_box_boundary` はphoto ray互換用です。periodic2 full-chordがbox外区間のrange limitへ達した場合だけ、
最初のbox eventまでに制限したqueryへ切り替えます。

`potential_barrier` は単一open faceに限り旧scalar energy式をevent位置と補間速度で評価します。複数open faceは
一般化せずfail closedであり、shared potential/gaugeに基づく物理モデルは後段の対象です。

`BEACH_WARN_LONG_PARTICLE_STEPS` を正整数で設定すると、長く生き残る粒子の診断出力を一定 step ごとに出します。

collision status は `collision_query_ok=0`、`collision_query_image_limit=1`、
`collision_query_index_range=2`、`collision_query_invalid_segment=3`、`collision_query_grid_stalled=4` です。
非有限 candidate はcollision前に `particle_step_invalid_boundary` として停止します。`status` を省略した public call は不完全な照会を miss として扱わず、
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

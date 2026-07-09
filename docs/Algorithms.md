title: BEACH アルゴリズム解説

# BEACH アルゴリズム解説

この文書は、BEACH の現行 Fortran 実装に基づいて、数値アルゴリズムと実行順序を説明します。
BEACH は格子 PIC ではなく、三角形境界要素上の電荷を source とする BEM 風の Coulomb 場評価と、
テスト粒子追跡を batch 単位で結合した表面帯電シミュレータです。

実装へのリンクは、現在のファイルと主な symbol を示します。行番号が後続変更でずれても、リンク先ファイルと symbol 名を優先してください。

---

## 目次

### Part I: 全体像
1. [BEACH の計算モデル](#1-beach-の計算モデル)
2. [実行エントリーポイントと初期化](#2-実行エントリーポイントと初期化)
3. [batch ループ](#3-batch-ループ)

### Part II: 場計算
4. [境界要素電荷による Coulomb 場](#4-境界要素電荷による-coulomb-場)
5. [direct / treecode / FMM の切替](#5-direct--treecode--fmm-の切替)
6. [periodic2 場境界](#6-periodic2-場境界)

### Part III: 粒子計算
7. [粒子生成と injection state](#7-粒子生成と-injection-state)
8. [Boris pusher](#8-boris-pusher)
9. [衝突検出](#9-衝突検出)
10. [box 境界条件](#10-box-境界条件)

### Part IV: 表面電荷と出力
11. [電荷堆積と表面モデル](#11-電荷堆積と表面モデル)
12. [統計・履歴・再開](#12-統計履歴再開)
13. [並列化と性能計測](#13-並列化と性能計測)

### Part V: FMM と batch 安定性
14. [Coulomb FMM コア詳細](#14-coulomb-fmm-コア詳細)
15. [`batch_duration` の安定性と定常値](#15-batch_duration-の安定性と定常値)

---

## 1. BEACH の計算モデル

**Source**:
[`bem_simulator`](../src/runtime/simulator/bem_simulator.f90),
[`bem_simulator_loop`](../src/runtime/simulator/bem_simulator_loop.f90),
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
G_\epsilon(\mathbf{r}) = \frac{1}{\sqrt{\|\mathbf{r}\|^2 + \epsilon^2}}
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
{\left(\|\mathbf{x} - \mathbf{c}_j\|^2 + \epsilon^2\right)^{3/2}}
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

`sim.batch_count` は累積到達 batch 数です。checkpoint が `batches=100` で `batch_count=150` なら、実行するのは 50 batch だけです。

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

### 3.3 particle processing

`process_particle_batch` は粒子ごとに最大 `sim.max_step` 回の時間発展を行います。

| 順序 | 処理 |
| --- | --- |
| 1 | 現在位置 `x0` と速度 `v0` を読む |
| 2 | `field_solver%eval_e(mesh, x0, e)` で境界要素電荷による電場を評価 |
| 3 | 一様外部電場 `sim.e0` を加える |
| 4 | `boris_push` で候補位置 `x1` と候補速度 `v1` を計算 |
| 5 | `find_first_hit(mesh, x0, x1, hit, sim=sim)` で線分衝突を調べる |
| 6 | hit があれば `q * w` を命中要素の `dq_thread(elem, tid)` へ加算し、粒子を吸収終了 |
| 7 | hit がなければ `apply_box_boundary` で open / reflect / periodic を適用 |
| 8 | 粒子が生存していれば `x` と `v` を更新して次 step へ進む |

`BEACH_WARN_LONG_PARTICLE_STEPS` を正整数で設定すると、長く生き残る粒子の診断出力を一定 step ごとに出します。

### 3.4 charge commit

`commit_batch_charge` は thread 別差分と photo emission 差分を統合します。

$$
\Delta q_i =
\sum_{\mathrm{thread}} \Delta q_{i,\mathrm{thread}}
+ \Delta q_{i,\mathrm{photo}}
$$

MPI 実行時は `mpi_allreduce_sum_real_dp_array` により rank 間で `dq` を和にします。
その後:

$$
q_i^{new} = q_i^{old} + \Delta q_i
$$

を適用し、表面モデル緩和を実行します。最後に実際に変化した電荷量で監視値を計算します。

$$
\mathrm{rel\_change}
=
\frac{\|\mathbf{q}^{new} - \mathbf{q}^{old}\|_2}
{\max(\|\mathbf{q}^{new}\|_2, q_\mathrm{floor})}
$$

この値は `stats%last_rel_change` と履歴出力に使われます。

---

## 4. 境界要素電荷による Coulomb 場

**Source**:
[`bem_field_solver`](../src/physics/field_solver/bem_field_solver.f90),
[`bem_field_solver_config`](../src/physics/field_solver/bem_field_solver_config.f90),
[`bem_field_solver_eval`](../src/physics/field_solver/bem_field_solver_eval.f90)

### 4.1 direct evaluation

direct mode では、評価点 `r` に対して全要素を直接足し込みます。

$$
\mathbf{E}(\mathbf{r}) =
k_c \sum_{i=1}^{N}
q_i
\frac{\mathbf{r} - \mathbf{c}_i}
{\left(\|\mathbf{r} - \mathbf{c}_i\|^2 + \epsilon^2\right)^{3/2}}
$$

計算量は 1 評価点あたり `O(nelem)` です。
粒子 step 数と粒子数が増えると支配的になるため、要素数が大きいケースでは treecode または FMM を使います。

### 4.2 length normalization

`sim.field_normalization` により、内部計算の長さスケール `L0` を選べます。

| 値 | `L0` |
| --- | --- |
| `si` | `1 m` |
| `length` | `sim.field_length_scale` |
| `box` | `max(box_max - box_min)` |
| `mesh` | mesh bounding box の最大幅 |

内部では

$$
\mathbf{x}' = \frac{\mathbf{x} - \mathbf{x}_0}{L_0},
\quad
\epsilon' = \frac{\epsilon}{L_0}
$$

として評価し、電場は `k_c / L0^2`、電位は `k_c / L0` を掛けて SI 単位へ戻します。
入力設定と出力 CSV は SI 単位のままです。

---

## 5. direct / treecode / FMM の切替

**Source**:
[`init_field_solver`](../src/physics/field_solver/bem_field_solver_config.f90),
[`refresh_field_solver`](../src/physics/field_solver/bem_field_solver_tree.f90),
[`eval_e_field_solver`](../src/physics/field_solver/bem_field_solver_eval.f90),
[14. Coulomb FMM コア詳細](#14-coulomb-fmm-コア詳細)

### 5.1 mode selection

`sim.field_solver` は次を受けます。

| 値 | 挙動 |
| --- | --- |
| `direct` | 常に direct 和 |
| `treecode` | octree + monopole 近似 |
| `fmm` | Coulomb FMM core |
| `auto` | `nelem >= tree_min_nelem` なら treecode、それ以外は direct |

`periodic2` は現行実装では `field_solver="fmm"` が必要です。

### 5.2 treecode

treecode は要素重心を octree に分割します。

1. `elem_order` に全要素 index を並べる。
2. node 内の要素重心 AABB を求める。
3. 要素数が `leaf_max` 以下、または分割不能なら leaf とする。
4. そうでなければ中心で 8 octant に分類し、子 node を再帰構築する。

`refresh_field_solver` は bottom-up に各 node の monopole を再集計します。

$$
Q_n = \sum_{i \in n} q_i
$$

$$
\mathbf{c}_{Q,n}
=
\begin{cases}
\frac{1}{Q_n}\sum_{i \in n} q_i \mathbf{c}_i & |Q_n| > 0 \\
\mathbf{c}_{node} & Q_n \approx 0
\end{cases}
$$

評価時には node 半径 `R` と評価点までの距離 `d` から遠方採用を判定し、採用できる node は monopole として評価します。
採用できない node は子へ降り、leaf では direct 和を行います。
電荷打ち消しが大きい node は monopole 誤差が出やすいため、`abs(Q) < charge_cancellation_tol * sum(abs(q_i))` の場合は遠方採用を抑制します。

### 5.3 FMM

FMM mode は simulator 非依存の Coulomb FMM core を呼びます。
field solver adapter は次を担当します。

1. メッシュ重心を source 座標 `src_pos(3, nelem)` に変換する。
2. `build_plan(plan, src_pos, options)` で幾何依存 plan を作る。
3. `update_state(plan, state, q_elem)` で電荷依存 state を更新する。
4. 粒子位置ごとに `eval_point(plan, state, r, e)` を呼ぶ。

FMM core の P2M / M2M / M2L / L2L / L2P の詳細は
[14. Coulomb FMM コア詳細](#14-coulomb-fmm-コア詳細) を参照してください。

---

## 6. periodic2 場境界

**Source**:
[`bem_field_solver_config`](../src/physics/field_solver/bem_field_solver_config.f90),
[`bem_coulomb_fmm_periodic`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic.f90),
[`bem_collision`](../src/physics/bem_collision.f90)

`sim.field_bc_mode="periodic2"` は、3 軸のうちちょうど 2 軸を周期軸として扱います。
周期軸は `bc_low(axis) == bc_high(axis) == periodic` で判定されます。
第三軸は開放方向です。

### 6.1 validation

`periodic2` では次が必須です。

- `sim.use_box = true`
- ちょうど 2 軸が periodic
- 各周期軸の `box_max - box_min > 0`
- `sim.field_solver = "fmm"`

`field_periodic_far_correction` は `auto`, `none`, `m2l_root_oracle` を受けます。
現行実装では `auto` は互換用に `none` へ正規化されます。
`m2l_root_oracle` は明示指定時だけ有効になる診断的な遠方補正です。

### 6.2 near images and far correction

`field_periodic_image_layers = N` は、周期 2 軸の近傍画像を

$$
(i, j) \in [-N, N]^2
$$

で列挙します。FMM core は primary cell source と画像 source を組み合わせて近傍寄与を扱います。
`m2l_root_oracle` を選ぶと、build 時に Ewald residual を使って root local 補正を fit し、runtime では root local へ注入します。

### 6.3 collision side

衝突判定側の periodic2 は、場計算の FMM とは別に処理されます。
`find_first_hit_periodic2` は粒子線分と canonical mesh AABB から必要な image shift 範囲を計算し、shift した線分で base mesh との交差を調べます。
同じ `t` に複数候補がある場合は、要素 index と image shift で deterministic に tie-break します。

---

## 7. 粒子生成と injection state

**Source**:
[`init_particle_batch_from_config`](../src/config/bem_app_config_runtime.f90),
[`bem_injection`](../src/particles/bem_injection.f90),
[`bem_sheath_runtime`](../src/physics/sheath/bem_sheath_runtime.f90)

### 7.1 species processing

`init_particle_batch_from_config` は全 species を走査し、rank ごとの生成数を決めてから SoA 粒子配列へ interleave します。
MPI 実行時は `mpi_split_count` により global count を rank に分配します。

対応する `source_mode` は次の通りです。

| source_mode | 生成数 | 位置 | 速度 |
| --- | --- | --- | --- |
| `volume_seed` | `npcls_per_step` | `pos_low` から `pos_high` の一様乱数 | shifted Maxwell |
| `reservoir_face` | flux と batch_duration から動的決定 | 注入面矩形 | 流入 flux で重み付けした Maxwell または velocity grid |
| `photo_raycast` | `rays_per_batch` ray の命中数 | ray が最初に命中した表面 | 表面法線方向の flux-weighted Maxwell |

### 7.2 reservoir_face macro count

drifting Maxwellian の流入 flux は、注入面 inward normal `n` に対する法線速度成分

$$
u_n = \mathbf{u} \cdot \mathbf{n}
$$

と熱速度

$$
\sigma = \sqrt{\frac{k_B T}{m}}
$$

から計算します。実装では `flux_weighted_normal_tail(vmin, u_n, sigma)` を使い、
指定された法線速度下限 `vmin_normal` 以上の粒子だけを流入粒子として数えます。

面積 `A`、batch duration `T_b`、マクロ粒子重み `w` に対して:

$$
N_\mathrm{phys} = \Gamma_\mathrm{in} A T_b
$$

$$
N_\mathrm{macro,expected} = \frac{N_\mathrm{phys}}{w}
$$

端数は `injection_state%macro_residual(species)` に繰り越します。

$$
B = r_\mathrm{old} + N_\mathrm{macro,expected}
$$

$$
N_\mathrm{macro} = \lfloor B \rfloor,\quad
r_\mathrm{new} = B - \lfloor B \rfloor
$$

MPI 実行時は `batch_duration_scale = 1 / nrank` を使い、各 rank が global flux の担当分だけ生成します。

### 7.3 reservoir sampling

`reservoir_face` の位置は注入面上の矩形から一様乱数で選び、必要なら `position_jitter_dt=sim.dt` により面から少し進んだ位置へずらします。
速度は次のいずれかです。

- shifted Maxwellian を生成し、法線成分だけ flux-weighted distribution から再サンプルする。
- `velocity_distribution="grid"` の場合、CSV velocity grid を読み、`phase_space` なら `max(v_n,0) f(v)`、`flux_weighted` なら入力値を流入分布として扱う。

`reservoir_potential_model="infinity_barrier"` では、注入面平均電位と `phi_infty` の差から法線速度下限を補正します。

### 7.4 photo_raycast

`photo_raycast` は注入面から ray を発射し、最初に命中した mesh element から粒子を放出します。

1. 注入面矩形から ray 起点を一様サンプルする。
2. `ray_direction` を正規化し、inward 方向であることを確認する。
3. ray を box 境界まで伸ばし、その線分で `find_first_hit` を実行する。
4. mesh に命中したら、その要素法線を使って放出位置と放出速度を作る。
5. box 境界に当たり、境界条件で反射または周期 wrap した場合は `raycast_max_bounce` まで継続する。

1 hit あたりの重みは:

$$
w_\mathrm{hit}
=
\frac{J_\perp A_\perp T_b}
{|q| N_\mathrm{ray,global}}
$$

ここで `A_perp = A * abs(dot(ray_direction, inward_normal))` です。

`deposit_opposite_charge_on_emit=true` の場合、放出元要素には

$$
\Delta q_\mathrm{emit} = -q_\mathrm{particle} w_\mathrm{eff}
$$

を `photo_emission_dq` として加算します。

### 7.5 photo escape closure

`photo_escape_model="boltzmann_cutoff"` では、放出元要素自身の寄与を除いた中心電位を使います。

$$
\mathrm{barrier} =
\max(\phi_\mathrm{emit} - \phi_\infty, 0)
$$

$$
f_\mathrm{escape} =
\exp\left(
-\frac{|q|\,\mathrm{barrier}}{k_B T_\mathrm{PE}}
\right)
$$

実効重みは

$$
w_\mathrm{eff} = w_\mathrm{hit} f_\mathrm{escape}
$$

です。これは戻り photoelectron を個別追跡せず、即時中和として扱う reduced closure です。

---

## 8. Boris pusher

**Source**:
[`bem_pusher`](../src/physics/bem_pusher.f90)

粒子運動は一様磁場 `sim.b0` と、粒子位置で評価した電場 `E` による Boris 法です。

入力:

- 位置 `x`
- 速度 `v`
- 電荷 `q`
- 質量 `m`
- 時間刻み `dt`
- 電場 `E`
- 磁束密度 `B`

更新式:

$$
\mathbf{v}^- =
\mathbf{v}^n
+ \frac{q}{m}\mathbf{E}\frac{\Delta t}{2}
$$

$$
\mathbf{t} =
\frac{q}{m}\mathbf{B}\frac{\Delta t}{2}
,\quad
\mathbf{s} =
\frac{2\mathbf{t}}{1 + \|\mathbf{t}\|^2}
$$

$$
\mathbf{v}' =
\mathbf{v}^- + \mathbf{v}^- \times \mathbf{t}
$$

$$
\mathbf{v}^+ =
\mathbf{v}^- + \mathbf{v}' \times \mathbf{s}
$$

$$
\mathbf{v}^{n+1} =
\mathbf{v}^+ + \frac{q}{m}\mathbf{E}\frac{\Delta t}{2}
$$

$$
\mathbf{x}^{n+1} =
\mathbf{x}^{n} + \mathbf{v}^{n+1}\Delta t
$$

BEACH では、この `x^n -> x^{n+1}` の線分に対して衝突判定を行います。
衝突があれば粒子は吸収され、`x^{n+1}` は粒子状態に保存されません。

---

## 9. 衝突検出

**Source**:
[`bem_collision`](../src/physics/bem_collision.f90),
[`bem_mesh`](../src/mesh/bem_mesh.f90)

### 9.1 broad phase

`init_mesh` は各要素の AABB と collision grid を構築します。
要素数が小さい場合は線形探索、大きい場合は一様グリッド + 3D-DDA を使います。

collision grid:

1. 全要素 AABB の bounding box を作る。
2. 目標 `target_elems_per_cell` からセル幅を見積もる。
3. 各要素 AABB が重なるセルへ要素 index を CSR 形式で登録する。

粒子線分 `p0 -> p1` はまず grid AABB と交差判定され、通過セルだけを 3D-DDA で列挙します。
各セルに登録された要素だけを narrow phase へ渡します。

### 9.2 narrow phase: Möller-Trumbore

線分

$$
\mathbf{p}(t) = \mathbf{p}_0 + t(\mathbf{p}_1 - \mathbf{p}_0), \quad 0 \le t \le 1
$$

と三角形

$$
\mathbf{v}(u,v) =
\mathbf{v}_0 + u(\mathbf{v}_1-\mathbf{v}_0) + v(\mathbf{v}_2-\mathbf{v}_0)
$$

の交差を Möller-Trumbore 法で判定します。
条件は次の通りです。

- 三角形が退化していない。
- 線分方向と三角形面がほぼ平行でない。
- `0 <= u <= 1`
- `0 <= v`
- `u + v <= 1`
- `0 <= t <= 1`

複数要素に命中した場合は、最小 `t` の要素を採用します。

### 9.3 periodic2 collision

`periodic2` では、mesh 本体は base element だけを持ちます。
`find_first_hit_periodic2` は線分と mesh AABB から必要な image shift 範囲を計算します。

$$
n_\mathrm{min}
=
\left\lceil
\frac{\min(p_0, p_1) - \max(mesh) - tol}{L}
\right\rceil
$$

$$
n_\mathrm{max}
=
\left\lfloor
\frac{\max(p_0, p_1) - \min(mesh) + tol}{L}
\right\rfloor
$$

各 image について、線分側を `-shift` して base mesh と交差判定します。
命中位置は物理 image 座標 `hit%pos` と primary cell へ折り返した `hit%pos_wrapped` の両方を持ちます。

---

## 10. box 境界条件

**Source**:
[`bem_boundary`](../src/physics/bem_boundary.f90)

粒子が mesh に衝突せず、更新候補位置が simulation box の外へ出た場合、軸ごとの境界条件を適用します。

| 境界条件 | 処理 |
| --- | --- |
| `open` | 粒子を消滅し、`escaped_boundary` として集計 |
| `reflect` | 境界面で位置を鏡映し、法線速度成分を反転 |
| `periodic` | 反対側へ wrap |

`apply_box_boundary` は 3 軸を順に見ます。
1 step で複数周期分を跨いだ場合でも、`periodic` では `modulo` により box 内へ戻します。
`reflect` と `periodic` では、境界上にぴったり残る数値不安定を避けるため、`1e-12` 程度の微小量で box 内側へ clamp します。

---

## 11. 電荷堆積と表面モデル

**Source**:
[`commit_batch_charge`](../src/runtime/simulator/bem_simulator_loop.f90),
[`bem_surface_models`](../src/physics/bem_surface_models.f90)

### 11.1 insulator accumulation

既定の `surface_model="insulator"` では、吸収された粒子の電荷がそのまま要素に蓄積されます。

粒子 `p` が要素 `i` に衝突した場合:

$$
\Delta q_i \mathrel{+}= q_p w_p
$$

この蓄積電荷は次 batch の field solver refresh で source 電荷として使われます。

### 11.2 photo emission bookkeeping

`photo_raycast` で `deposit_opposite_charge_on_emit=true` の場合、放出元要素へ逆符号電荷を加えます。
これは粒子の後続衝突による堆積とは別の `photo_emission_dq` として batch commit 時に統合されます。

### 11.3 floating conductor relaxation

`surface_model="conductor"` は、`field_bc_mode="free"` のときだけ使えます。
conductor 要素は `mesh_id` ごとに floating conductor group として扱われます。
目的は、各 conductor object の総電荷を保存しながら、同じ object 内の要素電位を等しくすることです。

未知量:

- conductor 要素電荷 `q_i`
- conductor group ごとの等電位値 `V_g`

要素 `i` が group `g(i)` に属するとき:

$$
\sum_j A_{ij} q_j - V_{g(i)}
=
-\phi_i^\mathrm{fixed}
$$

ここで

$$
A_{ij}
=
\begin{cases}
1/\epsilon & i=j,\ \epsilon>0 \\
2\sqrt{\pi}/h_i & i=j,\ \epsilon=0 \\
1/\sqrt{\|\mathbf{c}_i-\mathbf{c}_j\|^2+\epsilon^2} & i\ne j
\end{cases}
$$

です。`phi_fixed` は conductor 以外の電荷と一様外部電場が作る、`k_coulomb` で割った電位です。

各 group には総電荷保存制約を加えます。

$$
\sum_{i \in g} q_i = Q_g^\mathrm{before}
$$

この正方線形系を部分ピボット付き Gauss 消去で解き、conductor 要素の `q_elem` を置き換えます。

### 11.4 dielectric metadata

`surface_model="dielectric"` と `epsilon_r` は現行バージョンでは metadata です。
誘電分極や誘電境界条件の場への反映は実装されていません。

---

## 12. 統計・履歴・再開

**Source**:
[`bem_simulator_stats`](../src/runtime/simulator/bem_simulator_stats.f90),
[`bem_simulator_io`](../src/runtime/simulator/bem_simulator_io.f90),
[`bem_output_writer`](../src/runtime/bem_output_writer.f90),
[`bem_restart`](../src/runtime/bem_restart.f90)

### 12.1 batch outcomes

`count_batch_outcomes` は local rank の batch 粒子を次の 5 項目に集計します。

| index | 内容 |
| --- | --- |
| 1 | batch 内粒子総数 |
| 2 | mesh に吸収された粒子数 |
| 3 | escaped として数える粒子数 |
| 4 | box open 境界から出た粒子数 |
| 5 | `max_step` まで生存した粒子数 |

`absorbed_flag` が立っていない粒子のうち、open boundary で消えたものは `escaped_boundary`、最後まで `alive` のものは `survived_max_step` です。

MPI 実行時は batch count 配列を allreduce し、root 以外の rank の粒子も統計に含めます。

### 12.2 history output

`history_stride > 0` のとき `charge_history.csv` を書きます。
出力条件は:

$$
(stats.batches - 1) \bmod history\_stride = 0
$$

なので、batch 1 は常に出力対象です。

`output.write_potential_history=true` の場合、同じ stride で `potential_history.csv` も出します。
電位履歴は、その時点の `q_elem` で field solver を refresh し、要素重心電位を計算して書きます。

### 12.3 final output

`output.write_files=true` のとき、root rank が主な最終出力を書きます。

- `summary.txt`
- `charges.csv`
- `mesh_potential.csv`（有効時）
- `mesh_triangles.csv`
- `mesh_sources.csv`

全 rank は checkpoint 用に RNG state と macro residual を保存します。
MPI 実行時は rank 別ファイル名になります。

### 12.4 restart consistency

再開時は次を検証します。

- checkpoint の `mesh_nelem` が現在の mesh 要素数と一致する。
- MPI world size が前回と一致する。
- `summary.txt` の統計値が有限・非負である。
- `charges.csv` の要素数と電荷値が妥当である。
- RNG state と macro residual が読み込める。

必須 checkpoint がない場合、`output.resume=true` では新規実行へフォールバックせず停止します。

---

## 13. 並列化と性能計測

**Source**:
[`bem_mpi`](../src/core/bem_mpi.F90),
[`bem_performance_profile`](../src/runtime/bem_performance_profile.f90),
[`bem_simulator_loop`](../src/runtime/simulator/bem_simulator_loop.f90)

### 13.1 OpenMP

粒子追跡は OpenMP で粒子 index を並列化します。

- `dq_thread(nelem, nth)` により、衝突電荷は thread local に集計する。
- schedule は `dynamic, 1` で、粒子ごとの寿命差による load imbalance を抑える。
- field solver の refresh や treecode node 集計にも一部 OpenMP loop が使われる。

### 13.2 MPI

MPI 並列は粒子生成と粒子追跡を rank 分割します。
mesh と `q_elem` は各 rank が保持し、batch commit 時に `dq` を allreduce して全 rank の電荷状態を一致させます。

主な allreduce:

- `dq(nelem)` の和
- batch outcome counts の和

root rank だけが human-readable な最終 CSV と history を書きます。
RNG state と macro residual は rank ごとに保存されます。

### 13.3 performance profile

`BEACH_PROFILE=1` を設定すると、主要 region の時間を `performance_profile.csv` に出します。
主な region:

- `load_or_init`
- `field_solver_init`
- `prepare_batch`
- `field_refresh`
- `particle_batch`
- `commit_charge`
- `mpi_reduce`
- `stats_update`
- `history_write`
- `write_results`
- `write_checkpoint`

MPI 実行時の scaling 評価では、rank 間の最大時間 `rank_max_s` を見るのが推奨です。

---

## 14. Coulomb FMM コア詳細

この節は、現行 Fortran 実装の Coulomb FMM コア
[`bem_coulomb_fmm_core` module page](../module/bem_coulomb_fmm_core.html)
と、その実装を分割した関連ファイル群の仕様とアルゴリズムをまとめたものです。

- 公開 API / 境界: `src/physics/field_solver/fmm/api/`
- 内部共通実装: `src/physics/field_solver/fmm/internal/common/`
- tree / plan 実装: `src/physics/field_solver/fmm/internal/tree/`
- state / eval 実装: `src/physics/field_solver/fmm/internal/runtime/`
- periodic2 実装: `src/physics/field_solver/fmm/internal/periodic/`

対象は simulator 非依存の内部 API で、`mesh_type` や `sim_config` を直接 `use` しません。
BEACH 側では field solver adapter がこのコアを呼び出します。

### 1. 目的

この FMM コアの目的は、固定された source 点群 `src_pos(3,n)` と可変電荷 `src_q(n)` に対して、
多数の評価点で Coulomb 電場を高速に返すことです。

現行実装の設計目標は次の通りです。

- kernel は 3D Coulomb のみ
- source 幾何と電荷更新を分離する
- `free` と `periodic2` のみを対象にする
- 近傍 direct 和もコア内部に含める
- simulator からは配列 API だけが見えるようにする

### 2. 公開 API

コアが提供する主な手続きは次の 4 つです。

```fortran
call build_plan(plan, src_pos, options)
call update_state(plan, state, src_q)
call eval_points(plan, state, target_pos, e)
call eval_point(plan, state, r, e)
```

入力・出力の意味は次の通りです。

- `src_pos(3,n)`:
  source 点の座標。`build_plan` 後は固定とみなす。
- `src_q(n)`:
  source 点の電荷。`update_state` ごとに更新できる。
- `target_pos(3,m)` または `r(3)`:
  評価点。
- `e(3,m)` または `e(3)`:
  電場ベクトル。

注意点:

- コアが返す電場には `k_coulomb` を掛けていません。BEACH の adapter 側で最後に掛けます。
- `build_plan` は幾何依存処理、`update_state` は電荷依存処理です。
- `eval_point(s)` は `plan` と `state` が ready な前提です。

#### 2.2 C ABI / Python 連携

`src/physics/field_solver/bem_field_kernel_c.f90` は、この Fortran API を
`iso_c_binding` の opaque handle API として公開します。共有ライブラリは
`make build-kernel` で `build/libbeach_field_kernel.so` に生成します。

主な C ABI:

```text
beach_kernel_create(handle)
beach_kernel_destroy(handle)
beach_kernel_build(handle, src_pos, options...)
beach_kernel_update_charges(handle, src_q)
beach_kernel_eval_e(handle, target_pos, e)
beach_kernel_eval_phi(handle, target_pos, phi)
beach_kernel_force_on_charges(handle, target_pos, target_q, origin, force, torque)
```

Python 側は `beach.fortran_results.kernel.FieldKernel` がこの ABI を `ctypes`
で呼びます。`calc_object_forces_kernel` は object 自身の source 電荷をゼロにして
`sum(q_i E_not_self(r_i))` を評価するため、自己力の混入を避けながら
`periodic2 + m2l_root_oracle` を含む同じ field kernel を使えます。
`Beach.scene()` / `BeachScene` は Python 側で object の剛体移動・回転を一時的に
適用し、編集後の重心配列を同じ ABI に渡します。剛体変換の補助処理は NumPy
既定で、任意依存の Numba backend も選べますが、場評価そのものは Fortran
kernel が担当します。

#### 2.3 BEACH adapter での使い方

BEACH の field solver adapter は、メッシュ要素重心を `src_pos` としてこのコアへ渡します。

- 初期化時は `build_plan` の直後に `update_state` を実行します。
- その後の refresh では、メッシュ幾何が変わらない通常運用を前提に既存 `plan` を再利用し、`src_q` 更新として `update_state` だけを呼びます。
- `build_plan` と legacy tree metadata の同期をやり直すのは、plan 未構築時・source 数変更時・要素数 0 件で plan/state を破棄したときだけです。

### 3. データ構造

#### 3.1 `fmm_options_type`

主な内部オプション:

- `theta`: well-separated 判定用パラメータ
- `leaf_max`: source octree の葉に許す最大 source 数
- `order`: Cartesian 展開次数
- `softening`: softened Coulomb kernel の `epsilon`
- `use_periodic2`: 2 周期軸モードの有効化
- `periodic_axes(2)`, `periodic_len(2)`: 周期軸と周期長
- `periodic_image_layers`: 近傍画像和の層数 `N`
- `periodic_far_correction`: core が受ける値は `auto`, `none`, `m2l_root_oracle`。`periodic2` 有効時の `auto` は互換用に `none` へ正規化され、`m2l_root_oracle` は明示指定時だけ有効になる
- `periodic_ewald_alpha`, `periodic_ewald_layers`: `m2l_root_oracle` の build-time Ewald fit で使う分解パラメータと打切り深さ
- `target_box_min/max`: dual-target tree を作るときの box

BEACH の adapter は現状 `order = 4` を使いますが、コア自体は可変次数を受けられます。
`periodic2` の `auto` は `none` に正規化されます。`m2l_root_oracle` は遠方補正を明示的に有効化します。

#### 3.2 `fmm_plan_type`

幾何にだけ依存する不変データです。

- 多重指数テーブル `alpha`, `deriv_alpha`
- source octree
- optional target tree
- source 葉一覧 `source_leaf_nodes`
- target 葉一覧 `leaf_nodes`
- 近傍 list `near_start/near_nodes`
- 遠方 node list `far_start/far_nodes`
- M2L pair cache `m2l_target_nodes/m2l_source_nodes`
- periodic image shift 配列
- M2L 微分表 `m2l_deriv`
- P2M 基底表 `source_p2m_basis`
- M2M/L2L の平行移動用圧縮テーブル

#### 3.3 `fmm_state_type`

電荷に依存して毎回更新されるデータです。

- `src_q(n)`
- `multipole(ncoef, nnode)`
- `local(ncoef, n_target_nodes)`
- `multipole_active(nnode)`
- `local_active(n_target_nodes)`

`multipole` は source tree ノードごとの多重極係数、`local` は target tree ノードごとの局所展開係数です。
`*_active` は zero-node を早く飛ばすための 0/1 フラグです。

### 4. 数学的定義

#### 4.1 kernel

現行コアは softening 付き Coulomb kernel を使います。

$$
G_\epsilon(\mathbf{r}) = \frac{1}{\sqrt{\|\mathbf{r}\|^2 + \epsilon^2}}
$$

$$
\phi(\mathbf{x}) = \sum_j q_j \, G_\epsilon(\mathbf{x} - \mathbf{x}_j)
$$

$$
\mathbf{E}(\mathbf{x}) = - \nabla \phi(\mathbf{x})
$$

近傍 direct 和でも遠方展開でも、同じ $G_\epsilon$ を使います。

#### 4.2 多重指数

多重指数 $\alpha = (\alpha_x, \alpha_y, \alpha_z)$ を使います。

$$
|\alpha| = \alpha_x + \alpha_y + \alpha_z
$$

$$
\alpha! = \alpha_x! \, \alpha_y! \, \alpha_z!
$$

$$
\mathbf{r}^\alpha = r_x^{\alpha_x} r_y^{\alpha_y} r_z^{\alpha_z}
$$

#### 4.3 P2M

node center を $c$ とすると、葉ノードの multipole 係数は

$$
M_\alpha(c) = \sum_{j \in \text{leaf}} q_j
\frac{(\mathbf{x}_j - \mathbf{c})^\alpha}{\alpha!}
$$

で定義します。

#### 4.4 M2M

子ノード中心 $c_{\mathrm{child}}$ の係数を親中心 $c_{\mathrm{parent}}$ に平行移動して集約します。
$\mathbf{d} = c_{\mathrm{child}} - c_{\mathrm{parent}}$ とすると

$$
M_\beta(c_{\mathrm{parent}}) =
\sum_{\alpha \le \beta}
M_\alpha(c_{\mathrm{child}})
\frac{\mathbf{d}^{\beta-\alpha}}{(\beta-\alpha)!}
$$

現行実装では $\beta - \alpha$ に対応する index と
$\mathbf{d}^{\beta-\alpha} / (\beta-\alpha)!$ を `build_plan` 時に前計算します。

#### 4.5 M2L

source node 中心 $c_s$、target node 中心 $c_t$ に対して
$R = c_t - c_s$ とします。

局所展開係数は

$$
L_\alpha(c_t) \mathrel{+}=
\sum_\beta (-1)^{|\beta|}
M_\beta(c_s)
D^{\alpha+\beta} G_\epsilon(R)
$$

で更新します。

ここで $D^\gamma$ は multi-index 微分です。
現行実装は $D^{\alpha+\beta} G_\epsilon(R)$ を `m2l_deriv(:, pair)` として pair ごとに前計算します。

#### 4.6 L2L

親中心 $c_{\mathrm{parent}}$ の局所展開を子中心 $c_{\mathrm{child}}$ へ平行移動します。
$\mathbf{d} = c_{\mathrm{child}} - c_{\mathrm{parent}}$ とすると

$$
L_\alpha(c_{\mathrm{child}}) \mathrel{+}=
\sum_{\gamma \ge \alpha}
L_\gamma(c_{\mathrm{parent}})
\frac{\mathbf{d}^{\gamma-\alpha}}{(\gamma-\alpha)!}
$$

これも `build_plan` 時に shift monomial を前計算します。

#### 4.7 L2P

評価点 $x$ が属する target leaf の中心を $c_{\mathrm{leaf}}$ とし、
$\mathbf{dr} = x - c_{\mathrm{leaf}}$ とすると

$$
E_k(x) = - \sum_{|\alpha| < p} L_{\alpha + e_k}(c_{\mathrm{leaf}}) \frac{\mathbf{dr}^\alpha}{\alpha!}
$$

で電場を評価します。
ここで $e_k$ は軸 $k$ の単位 multi-index です。

### 5. `build_plan` のアルゴリズム

`build_plan` は幾何依存処理だけを行います。

#### 5.1 source tree

source 座標の bounding box を再帰的に 8 分割して octree を作ります。
停止条件は次のどちらかです。

- source 数 `<= leaf_max`
- bbox が十分に小さく、これ以上分割しても意味がない

#### 5.2 target topology

target 側は 2 通りあります。

- `target_box` が無効:
  source tree の葉をそのまま target leaf として使う
- `target_box` が有効:
  box 全体を覆う別 target tree を作る

`periodic2` では target point を box 内に wrap してから target leaf を探します。

#### 5.3 near/far と M2L pair cache

各 target leaf について source tree を再帰走査し、
near node と far node を作ります。

well-separated 判定は

$$
(r_s + r_t)^2 < \theta_{\mathrm{eff}}^2 \, \|\mathbf{d}\|^2
$$

です。

- $r_s$: source node 半径
- $r_t$: target node 半径
- $\mathbf{d}$: node center 間ベクトル
- $\theta_{\mathrm{eff}} = \theta$ for `free` と `periodic2`

`periodic2` では $\mathbf{d}$ に minimum-image 補正を入れます。

その後、dual-tree 再帰で M2L pair cache を作り、
target node ごとの索引配列も準備します。

#### 5.4 build 時の前計算

`build_plan` の最後で、refresh ごとに変わらない量を前計算します。

- `source_parent_of`
- `parent_of`
- `source_p2m_basis`
- `m2m_term_count`, `m2m_alpha_list`, `m2m_delta_list`
- `l2l_term_count`, `l2l_gamma_list`, `l2l_delta_list`
- `source_shift_monomial`
- `target_shift_monomial`
- `shift_axis1`, `shift_axis2`
- `periodic_ewald`
- `periodic_root_operator`
- `m2l_deriv`

この前計算により `update_state` は charge-dependent な加算だけに近づきます。

#### 5.5 擬似コード

```text
build_plan(src_pos, options):
  initialize_basis_tables(order)
  build_source_tree(src_pos)
  precompute_source_p2m_basis()
  build_target_topology(target_box)
  build_interactions()
  precompute_translation_operators()
  precompute_periodic2_ewald_data()
  precompute_periodic_root_operator()
  precompute_m2l_derivatives()
```

### 6. `update_state` のアルゴリズム

`update_state` は legacy 実装の refresh に相当する処理です。
source 座標は変わらず、`src_q` だけが変わる前提です。

#### 6.1 処理順

```text
update_state(plan, state, src_q):
  ensure_state_capacity()
  copy src_q
  clear active flags
  clear multipole/local only when the tree has no source leaves or no M2L pairs
  P2M on source leaves
  M2M bottom-up
  M2L on cached pairs
  L2L top-down
  mark state ready
```

#### 6.2 OpenMP 並列化

現行実装では、次の箇所に OpenMP を入れています。

- `update_state` 全体を 1 つの parallel region で囲み、その内側で `src_q` コピーと active flag 初期化
- `P2M`: source leaf ごとのループ
- `M2M`: 同一 depth の node ループ
- `M2L`: target node ごとのループ
- `L2L`: 同一 depth の node ループ
- `build_plan` 時の translation / M2L 微分前計算

各ループは「1 node 1 thread」になりやすいように書いてあり、
共有配列への更新は node 単位で独立させています。

#### 6.3 実装上の最適化

`update_state` では次の無駄を避けています。

- $\beta - \alpha$ の multi-index 差分を毎回計算しない
- 親子 center 差分のべき乗を毎回作り直さない
- `P2M` の monomial 基底を source ごとに build 時に前計算する
- `M2M/L2L` の有効な `(alpha, delta)` 項だけを圧縮して持ち、無効項を走査しない
- `M2L` では source node の active flag を見て zero-node を pair 単位で早期 skip する
- `M2L` で target node 列へ細かく何度も書かず、thread-local な `local_acc` にためてから戻す
- `P2M` で target leaf ではなく source leaf 専用 index を使う

### 7. `eval_point(s)` のアルゴリズム

評価時の処理は次の通りです。

```text
eval_point(r):
  if plan is not built or state is not ready:
    return zero vector

  if periodic2:
    wrap r into target box

  leaf = locate_target_leaf(r)
  if leaf not found or leaf is not mapped to a leaf slot:
    use direct sum over all sources
    if periodic2 and far correction is m2l_root_oracle:
      add exact periodic Ewald correction
    return

  evaluate local expansion at leaf center
  add near direct interactions
  root local already carries periodic root correction when enabled
```

#### 7.1 葉の特定

- `periodic2` では評価点を target box 内へ wrap してから探索する
- target tree があるときは target tree の葉を使う
- target tree が無いときは source tree の葉を使う
- leaf lookup に失敗した場合、あるいは leaf が tree の葉 slot に写像できない場合は direct fallback に入る

#### 7.2 近傍 direct

near list に入った source index については direct 和を取ります。
`periodic2` では `[-N, N] x [-N, N]` の画像シフトを陽に回します。
fallback でも同じ direct kernel を使いますが、`periodic2` で明示 `m2l_root_oracle` が有効なときは oracle 補正を別途加算します。

#### 7.3 box 外 fallback

dual-target tree を使う場合、評価点が target box の外に出ることがあります。
そのときは target leaf を持たないので、全 source に対する direct 和へ fallback します。
明示 `m2l_root_oracle` では build-time Ewald fit の teacher と同じ exact periodic correction を direct fallback へ足します。

#### 7.4 root 補正の位置

`m2l_root_oracle` の root 補正は `update_state` 側で `state%local(:, root)` に注入されます。
したがって通常の leaf 評価では、`eval_point(s)` は root 補正を再計算せず、`state` に載っている local 展開をそのまま使います。

### 8. `periodic2` と遠方補正

#### 8.1 `periodic2`

`periodic2` は「ちょうど 2 軸だけ周期、残り 1 軸は開放」です。

近傍画像和は

$$
i, j \in [-N, N]
$$

の有限画像を陽に足します。

M2L でも同じ画像シフト集合を使い、各 pair の derivative を画像和で前計算します。

#### 8.2 `periodic2` Ewald（Ewald2P）補正

`bem_coulomb_fmm_periodic_ewald.f90` は、2 周期・1 開放の Coulomb field に対する Ewald 形の補正を実装します。
ここでいう `exact` は「コードが実際に評価する有限和」を指します。理論上の無限和そのものではなく、`field_periodic_image_layers = N` と `field_periodic_ewald_layers = L` で real-space / reciprocal-space の打切り深さを決める build-time oracle です。

##### 8.2.1 記法

周期軸を `a_1, a_2`、開放軸を `f` とします。
周期長、セル面積、画像集合、逆格子集合を次のように置きます。

$$
L_1 = \texttt{periodic\_len(1)},\qquad
L_2 = \texttt{periodic\_len(2)},\qquad
A = L_1 L_2
$$

$$
\mathcal I_N = \{(i,j)\in\mathbb Z^2 \mid |i|,|j|\le N\},\qquad
\mathcal K_L = \{(m,n)\in\mathbb Z^2 \mid |m|,|n|\le L,\ (m,n)\neq(0,0)\}
$$

画像シフトと逆格子ベクトルは

$$
\mathbf L_{ij} = iL_1\,\mathbf e_{a_1} + jL_2\,\mathbf e_{a_2},\qquad
\mathbf k_{mn} = 2\pi\left(\frac{m}{L_1}\mathbf e_{a_1} + \frac{n}{L_2}\mathbf e_{a_2}\right)
$$

と書けます。ソース位置を \(\mathbf s\)、評価点を \(\mathbf r\) とし、

$$
\mathbf R_{ij} = \mathbf r - \mathbf s - \mathbf L_{ij},\qquad
R_{ij} = \|\mathbf R_{ij}\|,\qquad
z = (\mathbf r - \mathbf s)\cdot \mathbf e_f
$$

を導入します。以下では \(\alpha =\) `field_periodic_ewald_alpha`、\(\epsilon =\) `softening` とします。

##### 8.2.2 実空間項

`add_screened_point_charge` が実装している screened Coulomb field は

$$
\mathbf E_\alpha(\mathbf R)
=
q\left(
\frac{\operatorname{erfc}(\alpha R)}{R^3}
+\frac{2\alpha}{\sqrt{\pi}}\frac{e^{-\alpha^2 R^2}}{R^2}
\right)\mathbf R
$$

です。これはポテンシャル

$$
\Phi_\alpha(\mathbf R) = q\,\frac{\operatorname{erfc}(\alpha R)}{R}
$$

の勾配に一致します。

`add_softened_point_charge` が使う direct kernel は

$$
\mathbf E_\epsilon(\mathbf R)
=
q\,\frac{\mathbf R}{(R^2+\epsilon^2)^{3/2}}
$$

で、通常の runtime direct path と同じ softening を使います。

実装上の real-space 補正は

$$
\mathbf E_{\mathrm{real}}
=
\sum_{(i,j)\in\mathcal I_{N+L}} \mathbf E_\alpha(\mathbf R_{ij})
- \sum_{(i,j)\in\mathcal I_N} \mathbf E_\epsilon(\mathbf R_{ij})
$$

です。実装では `r2 <= tiny(1.0d0)` の項はスキップするため、self interaction は入らない扱いです。`add_periodic2_exact_ewald_correction_single_source` に direct fallback の \(\sum_{(i,j)\in\mathcal I_N}\mathbf E_\epsilon\) を足すと、inner image の softened 部分が打ち消され、outer shell 側は screened 形に置き換わります。

##### 8.2.3 逆空間項

`add_exact_periodic2_reciprocal_space_correction` が使う逆空間項は、\((m,n)\neq(0,0)\) に対して

$$
\theta_{mn} = \mathbf k_{mn}\cdot(\mathbf r-\mathbf s),\qquad
k_{mn} = \|\mathbf k_{mn}\|
$$

$$
G^\pm_{mn}(z)
=
e^{\pm k_{mn} z}\operatorname{erfc}\!\left(\frac{k_{mn}}{2\alpha}\pm \alpha z\right)
$$

を定義すると

$$
\mathbf E_{\mathrm{rec}}
=
q \sum_{(m,n)\in\mathcal K_L}
\frac{\pi}{A}
\begin{pmatrix}
\frac{(\mathbf k_{mn})_{a_1}}{k_{mn}}\sin\theta_{mn}\,\bigl(G^+_{mn}(z)+G^-_{mn}(z)\bigr) \\
\frac{(\mathbf k_{mn})_{a_2}}{k_{mn}}\sin\theta_{mn}\,\bigl(G^+_{mn}(z)+G^-_{mn}(z)\bigr) \\
\cos\theta_{mn}\,\bigl(G^-_{mn}(z)-G^+_{mn}(z)\bigr)
\end{pmatrix}
$$

です。コードでは `term_p`, `term_m`, `pair_sum` に対応します。
この式は、逆格子の `k=0` を除いた高周波成分に対応します。

##### 8.2.4 `k=0` 項

`add_exact_periodic2_k0_correction` が実装しているゼロモード補正は

$$
\mathbf E_0
=
q\,\frac{2\pi}{A}\operatorname{erf}(\alpha z)\,\mathbf e_f
$$

です。single-source の oracle では `k=0` の電場寄与としてこの形を保持します。

##### 8.2.5 実装される補正

以上をまとめると、`add_periodic2_exact_ewald_correction_single_source` が 1 粒子分に加える補正は

$$
\mathbf E_{\mathrm{corr}}
=
\mathbf E_{\mathrm{real}}
+ \mathbf E_{\mathrm{rec}}
+ \mathbf E_0
$$

です。`add_periodic2_exact_ewald_correction_all_sources` はまずこれを全ソースに対して総和します。

##### 8.2.6 `charged_walls` total-charge 補正

非中性 slab の `charged_walls` closure に合わせて、`add_periodic2_exact_ewald_correction_all_sources` は全ソース和のあとに total-charge 補正

$$
\mathbf E_{\mathrm{walls}}(z)
=
\begin{cases}
\ \ \frac{2\pi Q_{\mathrm{tot}}}{A}\,\mathbf e_f & (z < z_{\mathrm{low}}), \\\\
\ \ 0 & (z_{\mathrm{low}} \le z \le z_{\mathrm{high}}), \\\\
-\frac{2\pi Q_{\mathrm{tot}}}{A}\,\mathbf e_f & (z > z_{\mathrm{high}})
\end{cases}
$$

を加えます。ここで `A = L_1 L_2` は周期セル面積、`Q_tot = \sum_j q_j`、`z_low/high` は `target_box_min/max` の非周期軸境界です。
この項は 2 枚の補償壁の場に対応し、slab 内では厳密に打ち消し合うため、target box 内で build する root oracle や通常の粒子前進には影響しません。影響するのは target box 外へ落ちた direct fallback 評価だけです。

`field_periodic_ewald_alpha` が `<= 0` の場合、`resolve_periodic2_ewald_alpha` は

$$
\alpha = \frac{1.2}{(N+1)\min(L_1,L_2)}
$$

を自動選択します。`min(L_1,L_2)\le 0` なら `alpha = 0` として oracle を無効化します。
また内部では `kmax = max(1, field_periodic_ewald_layers)` として逆空間の有限和を作ります。

実際の runtime direct fallback は

$$
\mathbf E_{\mathrm{fallback}}
=
\sum_{(i,j)\in\mathcal I_N} \mathbf E_\epsilon(\mathbf R_{ij})
+
\mathbf E_{\mathrm{corr}}
+
\mathbf E_{\mathrm{walls}}
$$

です。`m2l_root_oracle` の build-time fit では check points が target box 内にあるため `\mathbf E_{\mathrm{walls}} = 0` となり、teacher には single-source の `\mathbf E_{\mathrm{corr}}` だけが使われます。`periodic_root_operator` 側では定数 potential mode を使わないため、monopole column は 0 に固定されます。

##### 8.2.7 `m2l_root_oracle`

`m2l_root_oracle` は、この Ewald2P 補正を teacher にして root multipole から root local への演算子を proxy/check 点で fit する明示 opt-in の高コスト診断モードです。通常運用では `none` を使います。

- `periodic_image_layers = N`: runtime で explicit に残す近傍画像殻
- `periodic_ewald_layers = L`: build-time oracle の real-space outer shell `N < max(|i|,|j|) <= N+L` と reciprocal cutoff `|m|, |n| <= L`
- `periodic_ewald_alpha = alpha`: Ewald 分解パラメータ。`<= 0` なら自動決定
- build では exact periodic Ewald correction を check points で評価し、field residual を least-squares fit して root local 演算子を作る
- runtime では `local(:, root) += T_root_oracle * multipole(:, root)` を加えるだけなので eval path に Ewald 本体は入らない
- tree 外 fallback では direct sum に exact periodic correction を直接足して、target box 外でも periodic residual を落とさない
- fit は potential ではなく field を使い、local の定数 potential mode は 0 に固定する

### 9. 計算量の見方

固定次数 $p$、bounded interaction list を仮定すると、実用上は次のように見てよいです。

- `build_plan`: $O(n \log n)$ に近い
- `update_state`: $O(n)$ に近い
- `eval_point`: $O(\log n + n_{\mathrm{near}} \, n_{\mathrm{img}}^2)$ に近い
- `eval_points`: 上記の点評価を target ごとに並列実行

厳密な定数因子は次に強く依存します。

- `order`
- `theta`
- `leaf_max`
- `periodic_image_layers`
- `periodic_ewald_layers`
- target tree の有無

### 10. 現行実装の制約

この FMM コアは汎用 kernel FMM ではありません。

- kernel は Coulomb 固定
- simulator adapter の既定次数は `order = 4`
- source 座標は `build_plan` 後に不変とみなす
- 対応境界は `free` と `periodic2`
- `periodic2` は正確に 2 周期軸が必要
- far correction は `none`（既定）, `auto`, `m2l_root_oracle`（`periodic2` の `auto` は互換用に `none` へ正規化、`m2l_root_oracle` は明示 opt-in）
- `eval_point(s)` の返り値には `k_coulomb` を含めない

### 11. 実装との対応

主な対応箇所:

- 公開 API / ラッパ:
  `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90`,
  `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_build.f90`,
  `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_state.f90`,
  `src/physics/field_solver/fmm/api/bem_coulomb_fmm_core_eval.f90`
- shared 型定義:
  `fmm_options_type`, `fmm_plan_type`, `fmm_state_type`
  （`src/physics/field_solver/fmm/internal/common/bem_coulomb_fmm_types.f90`）
- plan 構築:
  `build_plan`
  （`src/physics/field_solver/fmm/internal/tree/bem_coulomb_fmm_plan_ops.f90`）
- charge refresh:
  `update_state`, `p2m_leaf_moments`, `m2m_upward_pass`, `m2l_accumulate`, `l2l_downward_pass`
  （`src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_state_ops.f90`）
- 評価:
  `eval_point`, `eval_points`
  （`src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90`）
- periodic2 補助:
  `has_valid_target_box`, `use_periodic2_m2l_root_oracle`,
  `use_periodic2_root_operator`, `build_periodic_shift_values`, `add_point_charge_images_field`,
  `wrap_periodic2_point`, `apply_periodic2_minimum_image`, `distance_to_source_bbox`,
  `distance_to_source_bbox_periodic`
  （`src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic.f90`）
- periodic2 Ewald/oracle:
  `resolve_periodic2_ewald_alpha`, `precompute_periodic2_ewald_data`,
  `add_periodic2_exact_ewald_correction_single_source`, `add_periodic2_exact_ewald_correction_all_sources`
  （`src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_ewald.f90`）
- periodic2 root operator:
  `precompute_periodic_root_operator`
  （`src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90`）
- BEACH adapter:
  `src/physics/field_solver/bem_field_solver_config.f90`,
  `src/physics/field_solver/bem_field_solver_tree.f90`,
  `src/physics/field_solver/bem_field_solver_eval.f90`

設計上の責務分担は次の通りです。

- コア:
  幾何前処理、展開係数更新、近傍 direct、点評価
- BEACH adapter:
  `mesh_type` から `src_pos` を作る、`q_elem` を `src_q` へ流す、`k_coulomb` を最後に掛ける

---

## 15. `batch_duration` の安定性と定常値

この節は、BEACH のバッチループにおける `sim.batch_duration`（または `sim.batch_duration_step` から決まる 1 バッチの物理時間）と、収束した壁面電荷分布の正当性・安定性の関係を理論的に整理したものです。
現行実装では、`batch_duration_step` を使う場合は `sim.batch_duration = sim.dt * sim.batch_duration_step` と解釈され、`reservoir_face` 注入では 1 バッチの物理流入量からマクロ粒子数または重みが決まります。

実装上の入口は次のとおりです。

- バッチ計算手順: [`SPEC.md` §4](https://github.com/Nkzono99/BEACH/blob/main/SPEC.md)
- パラメータ定義: [Parameters](Parameters.html) の `sim.batch_duration` / `sim.batch_duration_step`
- 注入での使用: `src/particles/bem_injection.f90`（`reservoir_face` / `photo_raycast`）
- バッチ生成・重み解決: `src/config/bem_app_config_runtime.f90`

### 1. 連続時間モデルへの還元

絶縁体壁面の各要素 `j` の蓄積電荷を $q_j(t)$、そのときの入射電荷フラックス（壁単位面積あたり）を $J_j(\mathbf q)$ とすると、吸収のみを考える基本モデルは次の ODE になります。

$$
\frac{dq_j}{dt} \;=\; J_j(\mathbf q)\, A_j
$$

ここで $A_j$ は要素面積です。$J$ は壁面電荷が作る電場に依存するため、一般には **非線形** です（$J = J(\mathbf q)$）。

BEACH の 1 バッチ更新は、この連続時間モデルを **場をバッチ先頭で凍結した explicit な更新** として評価するものとみなせます。平均的には

$$
\mathbf q^{n+1} \;=\; \mathbf q^n \;+\; \Delta t_b \cdot \mathbf J(\mathbf q^n)\,\mathbf A \;+\; \boldsymbol\eta^n
$$

と表せます。ここで

- $\Delta t_b = $ `sim.batch_duration`
- $\mathbf A$ は要素面積ベクトル
- $\boldsymbol\eta^n$ はバッチ内 Monte Carlo サンプリング誤差

です。実装上も、バッチ内で粒子は同じ場 $E(\mathbf q^n)$ を見て進み、バッチ末尾で電荷差分がまとめて壁面へ反映されます。したがって、`batch_duration` はこの explicit 更新の時間刻みに相当します。

### 2. 定常値の正当性

平均更新写像を

$$
F_{\Delta t_b}(\mathbf q) \;=\; \mathbf q \;+\; \Delta t_b\, \mathbf J(\mathbf q)\,\mathbf A
$$

と書くと、その不動点 $\mathbf q^\*$ は

$$
F_{\Delta t_b}(\mathbf q^\*) = \mathbf q^\*
\quad\Longleftrightarrow\quad
\mathbf J(\mathbf q^\*) = 0
$$

で与えられます。したがって、**平均モデルの不動点自体は $\Delta t_b$ に依存しません**。

この意味で、反復が安定に収束しており、かつ Monte Carlo 誤差が十分に平均化されているなら、`batch_duration` を変えても目指している連続時間の定常解は同じです。

ただし、ここで言えるのはあくまで **平均モデルの固定点** についてです。実際の計算では

- バッチごとの有限サンプル誤差
- 収束判定に用いる監視量の揺らぎ
- 有限バッチ数で打ち切ることによる残差

が入るため、観測される収束値には弱い `batch_duration` 依存が残りえます。したがって、厳密には

> 反復の平均的な固定点は `batch_duration` に依存しないが、有限サンプル・有限時間の実計算では小さな step-size dependence が観測されうる

と書くのが安全です。

### 3. 線形安定性

不動点 $\mathbf q^\*$ 近傍で摂動 $\delta\mathbf q^n = \mathbf q^n - \mathbf q^\*$ を考えると、平均更新の線形化は

$$
\delta \mathbf q^{n+1} \;=\; \bigl(I + \Delta t_b\, M\bigr)\,\delta \mathbf q^n,
\qquad
M_{ij} \;\equiv\; \frac{\partial (J_i A_i)}{\partial q_j}\bigg|_{\mathbf q^\*}
$$

となります。一般の多自由度系での安定条件は、この更新行列のスペクトル半径に対する条件

$$
\rho\!\left(I + \Delta t_b\, M\right) < 1
$$

であり、各固有値 $\lambda_k$ に対して

$$
|1 + \Delta t_b\, \lambda_k| < 1
$$

が必要です。これが本質的な「BEACH 版の安定条件」です。

絶縁体壁が電荷を貯めるとそれ以上同種粒子を引き寄せにくくなる（負のフィードバック）ため、$M$ の主要な固有値 $\lambda_k$ は実負（$\mathrm{Re}(\lambda_k) < 0$）と期待されます。**この実負優勢モードの仮定** が成り立つときに限り、応答時間スケール $\tau_k \equiv 1/|\lambda_k|$ を用いて、最速モードに対して

$$
0 \;<\; \Delta t_b \;<\; \frac{2}{|\lambda_{\max}|} \;=\; 2\,\tau_{\min}
$$

で発散を避けられ、

$$
0 \;<\; \Delta t_b \;<\; \frac{1}{|\lambda_{\max}|} \;=\; \tau_{\min}
$$

で単調収束（過減衰）になります。

実務上は次のように読むのが適切です。

- $\Delta t_b < 2\,\tau_{\min}$ : 実負優勢モード仮定のもとでの非発散条件
- $\Delta t_b < \tau_{\min}$ : 同仮定のもとでの単調収束条件
- 一般の結合系では、厳密には $\rho(I + \Delta t_b\, M) < 1$ が本体

この意味で、$2\tau$ / $\tau$ 条件は「BEACH 版の CFL 条件」と断言するより、**実負優勢モードを仮定した explicit Euler 型の安定目安** と呼ぶ方が正確です。

### 4. Monte Carlo ノイズとの関係

ノイズ込みの 1 モード近似として

$$
\delta q^{n+1} \;=\; \left(1 - \frac{\Delta t_b}{\tau}\right)\,\delta q^n \;+\; \xi^n
$$

を考えると、定常分散は $\xi^n$ の分散に依存します。ここで重要なのは、**$\mathrm{Var}(\xi^n)$ の $\Delta t_b$ 依存性が、注入の正規化方式に依存する** ということです。BEACH の `reservoir_face` には 2 つの方式があります。

#### 4.1 `w_particle` 固定の場合

`w_particle` を直接指定し、物理流入数だけが $\Delta t_b$ に比例して増減する方式では、1 バッチあたりの期待マクロ粒子数は

$$
N_\text{macro} \;\propto\; \Delta t_b
$$

となります。バッチ電荷増分のショットノイズ分散も概ね $\propto \Delta t_b$ とみなせ、

$$
\mathrm{Var}(\xi^n) \;\approx\; \alpha\, \Delta t_b
$$

と置けます。$\Delta t_b \ll \tau$ の極限では定常分散は `batch_duration` に強くは依存しません。

#### 4.2 `target_macro_particles_per_batch` 固定の場合

一方、`target_macro_particles_per_batch` から `w_particle` を解く方式では、`src/config/bem_app_config_runtime.f90:644` のとおり

$$
w_\text{particle} \;\propto\; \frac{\Gamma\, A\, \Delta t_b}{N_\text{target}}
$$

の形で重みが決まるため、ノイズの $\Delta t_b$ 依存は §4.1 の単純な $\mathrm{Var}(\xi^n) \propto \Delta t_b$ とは変わります。マクロ粒子数は固定され、1 粒子あたりの寄与が $\Delta t_b$ に比例するためです。

#### 4.3 実務上の解釈

理論整理としては

- `batch_duration` は主に **deterministic な安定性** のつまみ
- 統計ノイズの主な制御つまみは **`w_particle` または `target_macro_particles_per_batch`**

と分離して考えるのがよいです。特に

> 「`batch_duration` を小さくすれば必ずノイズが下がる」
> 「`batch_duration` を大きくしてもノイズはほとんど変わらない」

のいずれも一般論としては言い切れません。そこは注入の正規化方式に依存します。

### 5. $\tau_{\min}$ の物理的見積もり

$\tau_{\min}$ は、数値安定性を支配する最速の有効応答時間です。ただし、これを 1 つの物理式で一般に与えるのは難しく、幾何・電位分布・上流分布関数・注入モデルに依存します。実務上は、次の 2 つを別物として見積もるのが安全です。

#### 5.1 充電・シース緩和時間

代表的には、ある有効容量 $C_\text{eff}$ と有効コンダクタンス $G_\text{eff}$ を用いて

$$
\tau_\text{charge} \;\sim\; \frac{C_\text{eff}}{G_\text{eff}}
$$

あるいは典型電位変化 $\Delta\phi$ と有効電流 $I_\text{eff}$ を用いて

$$
\tau_\text{charge} \;\sim\; \frac{C_\text{eff}\,\Delta\phi}{I_\text{eff}}
$$

と見積もるのが自然です。これは幾何や遮蔽の影響を受ける、比較的遅い charging timescale です。

#### 5.2 プラズマ周波数の逆数

もうひとつの速い基準は

$$
\tau_{pe} \;=\; \omega_{pe}^{-1} \;=\; \sqrt{\frac{\varepsilon_0\, m_e}{n_e\, e^2}}
$$

です。電子プラズマの微視的な速い時間スケールであり、系がどこまで急峻に応答しうるかを見る基準にはなります。

ただし、$\omega_{pe}^{-1}$ をそのまま $\tau_{\min}$ の上界とみなすのは強すぎます。むしろこれは **速い側の物理基準** であり、実際に `batch_duration` を制限する有効時定数は、幾何や入射律速を含んだ $\tau_\text{charge}$ 側で決まることが少なくありません。

#### 5.3 実用上の選び方

したがって $\tau_{\min}$ については

- $\omega_{pe}^{-1}$ : 微視的な速い基準
- $\tau_\text{charge}$ : 系固有の charging / sheath relaxation timescale

を別物として見積もり、最終的には数値実験で詰めるのが安全です。

> $\omega_{pe}^{-1}$ は速い基準にすぎず、実際の安定制約は $\tau_\text{charge}$ を含む有効応答時間で決まる

とまとめておくのが理論的に無理がありません。

### 6. 実用的な使い方

1. まず物理スケールとして $\omega_{pe}^{-1}$ と $\tau_\text{charge}$ の両方を概算する。
2. 初期の `batch_duration` は、振動を避けたいなら保守的に小さめから始める。
3. `batch_duration` を 1/2 倍と 2 倍に振って、電荷履歴や監視量の振る舞いを比較する（**step-size sensitivity check**）。
4. 収束値がほぼ一致し、かつ振動や発散傾向が見えなければ、その `batch_duration` は実用上十分と判断できる。
5. ノイズが大きい場合は、まず `w_particle` または `target_macro_particles_per_batch` を調整する。`batch_duration` の変更だけでノイズを解決しようとしない。
6. `charge_history.csv` の `last_rel_change` の振動や要素電荷時系列のジッタは有用な診断量だが、これは厳密な Richardson 外挿（誤差の冪乗則を仮定する手法）というより、**step-size sensitivity check** と呼ぶのが適切である。

### 7. まとめ

| 項目 | 結論 |
|---|---|
| 定常値の正当性 | 平均更新の固定点は `batch_duration` に依存しない |
| 厳密な安定条件 | $\rho(I + \Delta t_b\, M) < 1$ |
| $2\tau$, $\tau$ 条件 | 実負優勢モードを仮定した explicit Euler の近似目安 |
| $\omega_{pe}^{-1}$ の位置づけ | 微視的な速い基準であり、一般にそのまま安定上界ではない |
| ノイズと `batch_duration` | 依存性は注入の正規化方式に依存する |
| ノイズ低減の主手段 | `w_particle` または `target_macro_particles_per_batch` の調整 |
| 実務上の確認方法 | `batch_duration` を振った step-size sensitivity check |

「平均モデルの定常値が `batch_duration` の取り方に依存しない」ことは理論的にきれいに言えます。
一般的な安定条件 $\rho(I + \Delta t_b M) < 1$ も古典的な安定性解析からそのまま従います。
残る不確かさは $\tau_{\min}$ の値そのもので、ここだけはケース依存の物理見積もりと数値実験の両輪で詰める必要があります。

### 関連文書

- [Fortran パラメータファイル仕様](Parameters.html) — `sim.batch_duration` / `sim.batch_duration_step` の指定方法
- [Fortran 中心ワークフロー](Workflow.html) — バッチループの実行制御
- `SPEC.md` — 1 バッチの計算手順と停止条件

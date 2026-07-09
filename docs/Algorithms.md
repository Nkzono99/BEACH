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
[`docs/FMMCore.md`](FMMCore.md)

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
[`docs/FMMCore.md`](FMMCore.md) を参照してください。

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

title: 粒子追跡と表面電荷蓄積

Lang: [日本語](ParticleChargeLoop.md) | [English](ParticleChargeLoop.en.md)

# 粒子追跡と表面電荷蓄積

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

MPI 実行時は root rank が全 `batch_duration` に対する `N_macro` と残差を一度だけ更新し、その値を全 rank へ
broadcast します。その後、確定した整数個数を `mpi_split_count` で分配します。このため、global count と残差履歴は
MPI world size に依存しません。

### 7.3 reservoir sampling

`reservoir_face` の位置は注入面上の矩形から一様乱数で選び、必要なら `position_jitter_dt=sim.dt` により面から少し進んだ位置へずらします。
ジッタ後の位置は、周期軸ではprimitive cellへwrapし、非周期軸ではbox面へclampして有効box内に保ちます。
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
w_\mathrm{hit} =
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
[`bem_pusher`](../src/physics/bem_pusher.f90),
[`bem_particle_stepper`](../src/runtime/simulator/bem_particle_stepper.f90)

粒子運動は一様磁場 `sim.b0` と、予測中点で評価した境界要素電場に一様外部電場 `sim.e0` を1回加えた
電場による Boris 法です。box crossing候補では場評価点だけを周期軸でwrap・非周期軸でbox面へclampし、
軌道候補座標は写像しません。位置と速度の入力・出力は同一時刻の状態です。

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
\mathbf{x}_\mathrm{mid} =
\mathbf{x}^{n} + \frac{1}{2}\mathbf{v}^{n}\Delta t
,\quad
\mathbf{E}_\mathrm{mid} =
\mathbf{E}_\mathrm{BEM}(\mathbf{x}_\mathrm{mid}) + \mathbf{E}_0
$$

以下の速度更新では $\mathbf{E}=\mathbf{E}_\mathrm{mid}$ を使います。

$$
\mathbf{v}^- =
\mathbf{v}^n
{}+ \frac{q}{m}\mathbf{E}\frac{\Delta t}{2}
$$

$$
\mathbf{t} =
\frac{q}{m}\mathbf{B}\frac{\Delta t}{2}
,\quad
\mathbf{s} =
\frac{2\mathbf{t}}{1 + \lVert\mathbf{t}\rVert^2}
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
\mathbf{x}^{n} + \frac{1}{2}
\left(\mathbf{v}^{n} + \mathbf{v}^{n+1}\right)\Delta t
$$

予測中点の空間電場評価と台形位置更新により、smooth な場で candidate kinematics は二次精度になります。
一様電場による一定加速度の変位は丸め誤差まで解析解と一致します。
BEACH は `x^n -> x^{n+1}` を1回だけ衝突判定し、候補終点がbox内部ならその結果を確定します。box crossing時は
mesh hitと最初のfaceのparameterを比較して最早順序を決め、reflect/periodic後の残り時間を最大2 eventまで再積分します。
periodic2のfull-chord照会がbox外区間でrange limitに達した場合だけ、最初のfaceまでに制限して再照会します。
衝突があれば粒子は吸収され、候補状態は保存されません。残り時間中の3回目box eventは、そこまでにhitがなければ
stateを変更せず明示failureとし、`sim.dt`の縮小を要求します。

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
n_\mathrm{min} =
\left\lceil
\frac{\operatorname{min}(p_0, p_1) - \operatorname{max}(\mathrm{mesh}) - \mathrm{tol}}{L}
\right\rceil
$$

$$
n_\mathrm{max} =
\left\lfloor
\frac{\operatorname{max}(p_0, p_1) - \operatorname{min}(\mathrm{mesh}) + \mathrm{tol}}{L}
\right\rfloor
$$

各 image について、線分側を `-shift` して base mesh と交差判定します。
命中位置は物理 image 座標 `hit%pos` と primary cell へ折り返した `hit%pos_wrapped` の両方を持ちます。

`find_first_hit` の照会状態は次の public 定数で返されます。

| 定数 | 値 | 意味 |
| --- | ---: | --- |
| `collision_query_ok` | 0 | 必要な候補をすべて調べた完了照会 |
| `collision_query_image_limit` | 1 | 1軸または直積の image 数が安全上限 4096 を超えた |
| `collision_query_index_range` | 2 | shift bound が非有限、または i64 / i32 の表現範囲外 |

`status` を要求した caller には未完了状態を返します。`status` を省略した public call は未完了を「命中なし」として
継続せず、fail closed で停止します。通常の粒子追跡は OpenMP 内で各 query の status を受け、最小 particle / step を
named critical で集約して parallel region を抜けます。その後、全 rank から最小 rank の metadata を選択し、
全 rank が同じ batch / rank / particle / step / status message で停止します。

`photo_raycast` も各 ray query の status を必ず受け、最小 ray / bounce を専用 named critical で集約します。
未完了 ray はそこで処理を止め、OpenMP 終了後に species / ray / bounce / status を batch preparation へ返します。
main loop は field refresh と photo charge 処理の前に全 rank から最小 rank の metadata を選び、同一 message で停止するため、
失敗 rank の未完成 particle 配列や放出電荷差分は使用されません。

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

`field.element_kernel="triangle_p0"` でも蓄積量 `q_elem` は要素総電荷 [C] のままです。場評価時だけ `sigma=q_elem/area` へ変換します。衝突・吸収で使う ordered triangle と、電場の one-sided trace で使う `elem_vacuum_sign` は同じ mesh geometry から生成され、三角形 winding 自体は変更しません。

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
\sum_j A_{ij} q_j - V_{g(i)} =
-\phi_i^\mathrm{fixed}
$$

ここで

$$
A_{ij} =
\begin{cases}
1/\epsilon, & i=j,\ \epsilon>0, \\
2\sqrt{\pi}/h_i, & i=j,\ \epsilon=0, \\
1/\sqrt{\lVert\mathbf{c}_i-\mathbf{c}_j\rVert^2+\epsilon^2}, & i\ne j
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
- `charge_ledger.csv`

全 rank は checkpoint 用に rank 別 RNG state を保存し、root rank は共有 macro residual を単一ファイルへ保存します。

`charge_ledger.csv` は species ごとの注入、表面放出、吸収、無限遠 escape、未解決 discard、interface 往復の signed charge と粒子数を保持します。`summary.txt` の transactional residual と、species 間で相殺しない unresolved discard 絶対値は別の妥当性指標です。

### 12.4 restart consistency

再開時は次を検証します。

- checkpoint の `mesh_nelem` が現在の mesh 要素数と一致する。
- MPI world size が前回と一致する。
- `summary.txt` の統計値が有限・非負である。
- `charges.csv` の要素数と電荷値が妥当である。
- RNG state と macro residual が読み込める。
- schema v2/v3 の model / ordered mesh / ordered species fingerprint が一致する。
- 台帳を含む schema v2/v3 では stock と species 行が欠損していない。
- schema v3 の outer profile では `z, phi, E, rho` と solver scalar state が完全である。

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

root rank だけが human-readable な最終 CSV、history、global macro residual を書きます。
RNG state は rank ごとに保存されます。

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

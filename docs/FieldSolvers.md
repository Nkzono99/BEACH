title: 場ソルバーと境界条件

Lang: [日本語](FieldSolvers.md) | [English](FieldSolvers.en.md)

# 場ソルバーと境界条件

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
{\left(\lVert\mathbf{r} - \mathbf{c}_i\rVert^2 + \epsilon^2\right)^{3/2}}
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
[Coulomb FMM コア詳細](FMMCore.html)

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
\mathbf{c}_{Q,n} =
\begin{cases}
Q_n^{-1}\sum_{i \in n} q_i \mathbf{c}_i, & |Q_n| > 0, \\
\mathbf{c}_{\mathrm{node}}, & Q_n \approx 0
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
[Coulomb FMM コア詳細](FMMCore.html) を参照してください。

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

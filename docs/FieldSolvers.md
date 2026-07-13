title: 場ソルバーと境界条件

Lang: [日本語](FieldSolvers.md) | [English](FieldSolvers.en.md)

# 場ソルバーと境界条件

## まず選び方

| やりたいこと | 推奨設定 | 注意 |
| --- | --- | --- |
| 小さいメッシュで動作確認 | `field_solver = "auto"` または `"direct"` | `direct` は厳密だが `O(nelem)` |
| 要素数が多い通常計算 | `field_solver = "fmm"` | FMM の詳細は [Coulomb FMM コア詳細](FMMCore.html) |
| 2軸周期境界を使う | `field_bc_mode = "periodic2"` と `field_solver = "fmm"` | ちょうど 2 軸を periodic にする |
| 精度確認・デバッグ | 小ケースで `direct` と `fmm` を比較 | 同じ mesh / particle 条件で比較する |

`periodic2` は現行実装では FMM 専用です。`auto` は `field_bc_mode="free"` の小・中規模ケース向けと考えてください。

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
Stage 1 では、$A_n = \sum_{i \in n}|q_i|$ としたとき、幾何条件を満たす内部 node でも
`abs(A_n - abs(Q_n)) <= 64 * epsilon(1.0d0) * max(A_n, abs(Q_n))` を満たす場合だけ monopole として採用します。
この machine-epsilon tolerance は電荷集計の丸め差だけを許容します。
mixed-sign の内部 node は leaf まで子へ降りて direct interaction を行い、same-sign node は既存の monopole 経路と性能を維持します。

将来は monopole+dipole の誤差判定を導入し、不安定な signed charge centroid 近似を再導入せずに mixed-sign node の性能を回復する予定です。

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

`field_periodic_far_correction` の選択肢は次の通りです。

| 値 | 用途 | 備考 |
| --- | --- | --- |
| `none` | 有限画像近似 | 既定値 |
| `auto` | 互換入力 | 現在は`none`へ正規化 |
| `m2l_root_oracle` | correctness診断 | 明示opt-inの高コストfit |
| `cached_kneq0` | production無限周期非零モード | x/y periodic、z open、`exclude_k0`、lower boundary modelが必須 |

`cached_kneq0` と組み合わせる平均場closureは次の2種類です。

| lower boundary model | 下側平均場 | 上側平均場 | 用途 |
| --- | ---: | ---: | --- |
| `symmetric_vacuum` | $-Q/(2\epsilon_0A)$ | $+Q/(2\epsilon_0A)$ | 対称な真空半空間 |
| `e_bottom_zero` | $0$ | $Q/(\epsilon_0A)$ | 旧計算の再現 |

### 6.2 near images and far correction

`field_periodic_image_layers = N` は、周期 2 軸の近傍画像を

$$
(i, j) \in [-N, N]^2
$$

で列挙します。runtime場は選択したbackendに応じて次のように構成されます。

| 成分 | `none` | `m2l_root_oracle` | `cached_kneq0` |
| --- | --- | --- | --- |
| primary + near images | FMM | FMM | FMM |
| 遠方非零モード | なし | build時にEwald residualをroot localへfit | versioned cached operator |
| 対称$k=0$除去 | なし | oracle内 | state更新時に構築しevalで減算 |
| 物理的$k=0$ | legacy field contract | oracle contract | snapshotのboundary providerが1回合成 |

cache missではrank 0がfilesystem lockとatomic renameを担当し、operator計算はMPI/OpenMPへ分配します。
warm field evaluationとcharge refreshにはall-source Ewald和がありません。生成手順と測定値は
[cached periodic nonzero operator](FMMCore.html#101-cached-periodic-nonzero-operator)を参照してください。

### 6.3 collision side

衝突判定側の periodic2 は、場計算の FMM とは別に処理されます。
`find_first_hit_periodic2` は粒子線分と canonical mesh AABB から必要な image shift 範囲を計算し、shift した線分で base mesh との交差を調べます。
同じ `t` に複数候補がある場合は、要素 index と image shift で deterministic に tie-break します。

---

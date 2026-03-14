# Fortran Coulomb FMM コア仕様

この文書は、現行 Fortran 実装の Coulomb FMM コア
[`src/physics/field_solver/fmm/bem_coulomb_fmm_core.f90`](../src/physics/field_solver/fmm/bem_coulomb_fmm_core.f90)
と、その実装を分割した関連ファイル群の仕様とアルゴリズムをまとめたものです。

- 公開 API / 境界: `src/physics/field_solver/fmm/`
- 内部実装: `src/physics/field_solver/fmm/internal/`

対象は simulator 非依存の内部 API で、`mesh_type` や `sim_config` を直接 `use` しません。
BEACH 側では field solver adapter がこのコアを呼び出します。

## 1. 目的

この FMM コアの目的は、固定された source 点群 `src_pos(3,n)` と可変電荷 `src_q(n)` に対して、
多数の評価点で Coulomb 電場を高速に返すことです。

現行実装の設計目標は次の通りです。

- kernel は 3D Coulomb のみ
- source 幾何と電荷更新を分離する
- `free` と `periodic2` のみを対象にする
- 近傍 direct 和もコア内部に含める
- simulator からは配列 API だけが見えるようにする

## 2. 公開 API

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

## 3. データ構造

### 3.1 `fmm_options_type`

主な内部オプション:

- `theta`: well-separated 判定用パラメータ
- `leaf_max`: source octree の葉に許す最大 source 数
- `order`: Cartesian 展開次数
- `softening`: softened Coulomb kernel の `epsilon`
- `use_periodic2`: 2 周期軸モードの有効化
- `periodic_axes(2)`, `periodic_len(2)`: 周期軸と周期長
- `periodic_image_layers`: 近傍画像和の層数 `N`
- `periodic_far_correction`: `none` または `ewald_like`
- `periodic_ewald_alpha`, `periodic_ewald_layers`: `ewald_like` 補正の係数
- `target_box_min/max`: dual-target tree を作るときの box

BEACH の adapter は現状 `order = 4` を使いますが、コア自体は可変次数を受けられます。

### 3.2 `fmm_plan_type`

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

### 3.3 `fmm_state_type`

電荷に依存して毎回更新されるデータです。

- `src_q(n)`
- `multipole(ncoef, nnode)`
- `local(ncoef, n_target_nodes)`
- `multipole_active(nnode)`
- `local_active(n_target_nodes)`

`multipole` は source tree ノードごとの多重極係数、`local` は target tree ノードごとの局所展開係数です。
`*_active` は zero-node を早く飛ばすための 0/1 フラグです。

## 4. 数学的定義

### 4.1 kernel

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

近傍 direct 和でも遠方展開でも、同じ `G_\epsilon` を使います。

### 4.2 多重指数

多重指数 `alpha = (alpha_x, alpha_y, alpha_z)` を使います。

$$
|\alpha| = \alpha_x + \alpha_y + \alpha_z
$$

$$
\alpha! = \alpha_x! \, \alpha_y! \, \alpha_z!
$$

$$
\mathbf{r}^\alpha = r_x^{\alpha_x} r_y^{\alpha_y} r_z^{\alpha_z}
$$

### 4.3 P2M

node center を `c` とすると、葉ノードの multipole 係数は

$$
M_\alpha(c) = \sum_{j \in \text{leaf}} q_j
\frac{(\mathbf{x}_j - \mathbf{c})^\alpha}{\alpha!}
$$

で定義します。

### 4.4 M2M

子ノード中心 `c_child` の係数を親中心 `c_parent` に平行移動して集約します。
`d = c_child - c_parent` とすると

$$
M_\beta(c_{\mathrm{parent}}) =
\sum_{\alpha \le \beta}
M_\alpha(c_{\mathrm{child}})
\frac{\mathbf{d}^{\beta-\alpha}}{(\beta-\alpha)!}
$$

現行実装では `beta - alpha` に対応する index と
`d^(beta-alpha) / (beta-alpha)!` を `build_plan` 時に前計算します。

### 4.5 M2L

source node 中心 `c_s`、target node 中心 `c_t` に対して
`R = c_t - c_s` とします。

局所展開係数は

$$
L_\alpha(c_t) \mathrel{+}=
\sum_\beta (-1)^{|\beta|}
M_\beta(c_s)
D^{\alpha+\beta} G_\epsilon(R)
$$

で更新します。

ここで `D^gamma` は multi-index 微分です。
現行実装は `D^{alpha+beta} G_\epsilon(R)` を `m2l_deriv(:, pair)` として pair ごとに前計算します。

### 4.6 L2L

親中心 `c_parent` の局所展開を子中心 `c_child` へ平行移動します。
`d = c_child - c_parent` とすると

$$
L_\alpha(c_{\mathrm{child}}) \mathrel{+}=
\sum_{\gamma \ge \alpha}
L_\gamma(c_{\mathrm{parent}})
\frac{\mathbf{d}^{\gamma-\alpha}}{(\gamma-\alpha)!}
$$

これも `build_plan` 時に shift monomial を前計算します。

### 4.7 L2P

評価点 `x` が属する target leaf の中心を `c_leaf` とし、
`dr = x - c_leaf` とすると

$$
E_k(x) =
- \sum_{|\alpha| < p}
L_{\alpha + e_k}(c_{\mathrm{leaf}})
\frac{dr^\alpha}{\alpha!}
$$

で電場を評価します。
ここで `e_k` は軸 `k` の単位 multi-index です。

## 5. `build_plan` のアルゴリズム

`build_plan` は幾何依存処理だけを行います。

### 5.1 source tree

source 座標の bounding box を再帰的に 8 分割して octree を作ります。
停止条件は次のどちらかです。

- source 数 `<= leaf_max`
- bbox が十分に小さく、これ以上分割しても意味がない

### 5.2 target topology

target 側は 2 通りあります。

- `target_box` が無効:
  source tree の葉をそのまま target leaf として使う
- `target_box` が有効:
  box 全体を覆う別 target tree を作る

`periodic2` では target point を box 内に wrap してから target leaf を探します。

### 5.3 near/far と M2L pair cache

各 target leaf について source tree を再帰走査し、
near node と far node を作ります。

well-separated 判定は

$$
(r_s + r_t)^2 < \theta_{\mathrm{eff}}^2 \, \|\mathbf{d}\|^2
$$

です。

- `r_s`: source node 半径
- `r_t`: target node 半径
- `d`: node center 間ベクトル
- `theta_eff = theta` for `free`
- `theta_eff = theta / 2` for `periodic2`

`periodic2` では `d` に minimum-image 補正を入れます。

その後、dual-tree 再帰で M2L pair cache を作り、
target node ごとの索引配列も準備します。

### 5.4 build 時の前計算

`build_plan` の最後で、refresh ごとに変わらない量を前計算します。

- `source_parent_of`
- `parent_of`
- `source_p2m_basis`
- `m2m_term_count`, `m2m_alpha_list`, `m2m_delta_list`
- `l2l_term_count`, `l2l_gamma_list`, `l2l_delta_list`
- `source_shift_monomial`
- `target_shift_monomial`
- `m2l_deriv`

この前計算により `update_state` は charge-dependent な加算だけに近づきます。

### 5.5 擬似コード

```text
build_plan(src_pos, options):
  initialize_basis_tables(order)
  build_source_tree(src_pos)
  build_target_topology(target_box)
  build_interactions()
  precompute_translation_operators()
  precompute_m2l_derivatives()
```

## 6. `update_state` のアルゴリズム

`update_state` は legacy 実装の refresh に相当する処理です。
source 座標は変わらず、`src_q` だけが変わる前提です。

### 6.1 処理順

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

### 6.2 OpenMP 並列化

現行実装では、次の箇所に OpenMP を入れています。

- `update_state` 全体を 1 つの parallel region で囲み、その内側で `src_q` コピーと active flag 初期化
- `P2M`: source leaf ごとのループ
- `M2M`: 同一 depth の node ループ
- `M2L`: target node ごとのループ
- `L2L`: 同一 depth の node ループ
- `build_plan` 時の translation / M2L 微分前計算

各ループは「1 node 1 thread」になりやすいように書いてあり、
共有配列への更新は node 単位で独立させています。

### 6.3 実装上の最適化

`update_state` では次の無駄を避けています。

- `beta - alpha` の multi-index 差分を毎回計算しない
- 親子 center 差分のべき乗を毎回作り直さない
- `P2M` の monomial 基底を source ごとに build 時に前計算する
- `M2M/L2L` の有効な `(alpha, delta)` 項だけを圧縮して持ち、無効項を走査しない
- `M2L` では source node の active flag を見て zero-node を pair 単位で早期 skip する
- `M2L` で target node 列へ細かく何度も書かず、thread-local な `local_acc` にためてから戻す
- `P2M` で target leaf ではなく source leaf 専用 index を使う

## 7. `eval_point(s)` のアルゴリズム

評価時の処理は次の通りです。

```text
eval_point(r):
  if periodic2:
    wrap r into target box

  leaf = locate_target_leaf(r)
  if leaf not found:
    use direct sum over all sources
    add ewald-like correction if needed
    return

  evaluate local expansion at leaf center
  add near direct interactions
  add ewald-like far correction if needed
```

### 7.1 葉の特定

- target tree があるときはその葉を使う
- target tree が無いときは source 葉を使う
- `periodic2` では評価点を box に wrap してから探索する

### 7.2 近傍 direct

near list に入った source node については direct 和を取ります。
`periodic2` では `[-N, N] x [-N, N]` の画像シフトを陽に回します。

### 7.3 box 外 fallback

dual-target tree を使う場合、評価点が target box の外に出ることがあります。
そのときは target leaf を持たないので、全 source に対する direct 和へ fallback します。

## 8. `periodic2` と `ewald_like`

### 8.1 `periodic2`

`periodic2` は「ちょうど 2 軸だけ周期、残り 1 軸は開放」です。

近傍画像和は

$$
i, j \in [-N, N]
$$

の有限画像を陽に足します。

M2L でも同じ画像シフト集合を使い、各 pair の derivative を画像和で前計算します。

### 8.2 `ewald_like`

`ewald_like` は厳密な Ewald 和ではありません。
近傍有限画像の外側を、source node を 1 個の代表点電荷として近似し、
screened kernel で補正します。

代表点は node の monopole と一次モーメントから復元します。

$$
q = M_{(0,0,0)}
$$

$$
\mathbf{x}_{\mathrm{rep}}
= \mathbf{c}_{\mathrm{node}}
  + \frac{1}{q}
    \begin{bmatrix}
      M_{(1,0,0)} \\
      M_{(0,1,0)} \\
      M_{(0,0,1)}
    \end{bmatrix}
$$

screened point charge の補正電場は

$$
\mathbf{E}_{\mathrm{screen}} =
q \left(
\frac{\operatorname{erfc}(\alpha r)}{r^3}
+ \frac{2 \alpha}{\sqrt{\pi}}
\frac{e^{-\alpha^2 r^2}}{r^2}
\right) \mathbf{d}
$$

です。

- `d = x - x_rep`
- `r = ||d||`

補正対象の画像は、近傍層 `N` の外側から `N + L` 層までです。

## 9. 計算量の見方

固定次数 `p`、bounded interaction list を仮定すると、実用上は次のように見てよいです。

- `build_plan`: `O(n log n)` に近い
- `update_state`: `O(n)` に近い
- `eval_point`: `O(log n + n_near * n_img^2)` に近い
- `eval_points`: 上記の点評価を target ごとに並列実行

厳密な定数因子は次に強く依存します。

- `order`
- `theta`
- `leaf_max`
- `periodic_image_layers`
- `periodic_ewald_layers`
- target tree の有無

## 10. 現行実装の制約

この FMM コアは汎用 kernel FMM ではありません。

- kernel は Coulomb 固定
- simulator adapter の既定次数は `order = 4`
- source 座標は `build_plan` 後に不変とみなす
- 対応境界は `free` と `periodic2`
- `periodic2` は正確に 2 周期軸が必要
- far correction は `none` と `ewald_like` のみ
- `eval_point(s)` の返り値には `k_coulomb` を含めない

## 11. 実装との対応

主な対応箇所:

- API と型定義:
  `fmm_options_type`, `fmm_plan_type`, `fmm_state_type`
- plan 構築:
  `build_plan`
- charge refresh:
  `update_state`
- upward/downward:
  `p2m_leaf_moments`, `m2m_upward_pass`, `m2l_accumulate`, `l2l_downward_pass`
- 評価:
  `eval_point`, `eval_points`
- periodic2 補助:
  `add_point_charge_images_field`, `add_periodic2_ewald_like_correction`

設計上の責務分担は次の通りです。

- コア:
  幾何前処理、展開係数更新、近傍 direct、点評価
- BEACH adapter:
  `mesh_type` から `src_pos` を作る、`q_elem` を `src_q` へ流す、`k_coulomb` を最後に掛ける

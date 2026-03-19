# Fortran Coulomb FMM コア仕様

この文書は、現行 Fortran 実装の Coulomb FMM コア
[`src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90`](../src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90)
と、その実装を分割した関連ファイル群の仕様とアルゴリズムをまとめたものです。

- 公開 API / 境界: `src/physics/field_solver/fmm/api/`
- 内部共通実装: `src/physics/field_solver/fmm/internal/common/`
- tree / plan 実装: `src/physics/field_solver/fmm/internal/tree/`
- state / eval 実装: `src/physics/field_solver/fmm/internal/runtime/`
- periodic2 実装: `src/physics/field_solver/fmm/internal/periodic/`

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

### 2.1 BEACH adapter での使い方

BEACH の field solver adapter は、メッシュ要素重心を `src_pos` としてこのコアへ渡します。

- 初期化時は `build_plan` の直後に `update_state` を実行します。
- その後の refresh では、メッシュ幾何が変わらない通常運用を前提に既存 `plan` を再利用し、`src_q` 更新として `update_state` だけを呼びます。
- `build_plan` と legacy tree metadata の同期をやり直すのは、plan 未構築時・source 数変更時・要素数 0 件で plan/state を破棄したときだけです。

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
- `periodic_far_correction`: core が受ける値は `none`, `m2l_root_trunc`, `m2l_root_oracle`
- `periodic_ewald_alpha`, `periodic_ewald_layers`: `m2l_root_oracle` の build-time Ewald fit で使う分解パラメータと打切り深さ。`m2l_root_trunc` では `periodic_ewald_layers` を truncated root operator に畳み込む外側画像殻幅として使う
- `target_box_min/max`: dual-target tree を作るときの box

BEACH の adapter は現状 `order = 4` を使いますが、コア自体は可変次数を受けられます。
legacy な `m2l_root` は BEACH adapter 側で `m2l_root_trunc` に正規化されます。

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

近傍 direct 和でも遠方展開でも、同じ $G_\epsilon$ を使います。

### 4.2 多重指数

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

### 4.3 P2M

node center を $c$ とすると、葉ノードの multipole 係数は

$$
M_\alpha(c) = \sum_{j \in \text{leaf}} q_j
\frac{(\mathbf{x}_j - \mathbf{c})^\alpha}{\alpha!}
$$

で定義します。

### 4.4 M2M

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

### 4.5 M2L

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

### 4.6 L2L

親中心 $c_{\mathrm{parent}}$ の局所展開を子中心 $c_{\mathrm{child}}$ へ平行移動します。
$\mathbf{d} = c_{\mathrm{child}} - c_{\mathrm{parent}}$ とすると

$$
L_\alpha(c_{\mathrm{child}}) \mathrel{+}=
\sum_{\gamma \ge \alpha}
L_\gamma(c_{\mathrm{parent}})
\frac{\mathbf{d}^{\gamma-\alpha}}{(\gamma-\alpha)!}
$$

これも `build_plan` 時に shift monomial を前計算します。

### 4.7 L2P

評価点 $x$ が属する target leaf の中心を $c_{\mathrm{leaf}}$ とし、
$\mathbf{dr} = x - c_{\mathrm{leaf}}$ とすると

$$
E_k(x) = - \sum_{|\alpha| < p} L_{\alpha + e_k}(c_{\mathrm{leaf}}) \frac{\mathbf{dr}^\alpha}{\alpha!}
$$

で電場を評価します。
ここで $e_k$ は軸 $k$ の単位 multi-index です。

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

- $r_s$: source node 半径
- $r_t$: target node 半径
- $\mathbf{d}$: node center 間ベクトル
- $\theta_{\mathrm{eff}} = \theta$ for `free` と `periodic2`

`periodic2` では $\mathbf{d}$ に minimum-image 補正を入れます。

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
- `shift_axis1`, `shift_axis2`
- `periodic_ewald`
- `periodic_root_operator`
- `m2l_deriv`

この前計算により `update_state` は charge-dependent な加算だけに近づきます。

### 5.5 擬似コード

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

- $\beta - \alpha$ の multi-index 差分を毎回計算しない
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
  if plan is not built or state is not ready:
    return zero vector

  if periodic2:
    wrap r into target box

  leaf = locate_target_leaf(r)
  if leaf not found or leaf is not mapped to a leaf slot:
    use direct sum over all sources
    if periodic2 and far correction is trunc:
      add truncated far-image correction
    if periodic2 and far correction is oracle:
      add exact periodic Ewald correction
    return

  evaluate local expansion at leaf center
  add near direct interactions
  root local already carries periodic root correction when enabled
```

### 7.1 葉の特定

- `periodic2` では評価点を target box 内へ wrap してから探索する
- target tree があるときは target tree の葉を使う
- target tree が無いときは source tree の葉を使う
- leaf lookup に失敗した場合、あるいは leaf が tree の葉 slot に写像できない場合は direct fallback に入る

### 7.2 近傍 direct

near list に入った source index については direct 和を取ります。
`periodic2` では `[-N, N] x [-N, N]` の画像シフトを陽に回します。
fallback でも同じ direct kernel を使いますが、`periodic2` の far correction は mode に応じて別途加算されます。

### 7.3 box 外 fallback

dual-target tree を使う場合、評価点が target box の外に出ることがあります。
そのときは target leaf を持たないので、全 source に対する direct 和へ fallback します。
`m2l_root_trunc` では outer image shell の打切り和を direct fallback へ足し、`m2l_root_oracle` では build-time Ewald fit で得た exact periodic correction を direct fallback へ足します。

### 7.4 root 補正の位置

`m2l_root_trunc` と `m2l_root_oracle` の root 補正は `update_state` 側で `state%local(:, root)` に注入されます。
したがって通常の leaf 評価では、`eval_point(s)` は root 補正を再計算せず、`state` に載っている local 展開をそのまま使います。

## 8. `periodic2` と遠方補正

### 8.1 `periodic2`

`periodic2` は「ちょうど 2 軸だけ周期、残り 1 軸は開放」です。

近傍画像和は

$$
i, j \in [-N, N]
$$

の有限画像を陽に足します。

M2L でも同じ画像シフト集合を使い、各 pair の derivative を画像和で前計算します。

### 8.2 `m2l_root_trunc`

`m2l_root_trunc` は free-space の P2M / M2M / M2L / L2L / L2P をそのまま使い、近傍画像殻の外側だけを truncated far-image root operator として前計算して root local へ足すモードです。

- `periodic_image_layers = N`: explicit に残す近傍画像殻
- `periodic_ewald_layers = L`: truncated root operator に畳み込む外側画像殻 `N < max(|i|,|j|) <= N+L`
- runtime では `local(:, root) += T_per_root * multipole(:, root)` を加えた後に既存 L2L を実行
- 現行実装は root operator の monopole column を 0 にしており、constant/gauge mode を外した field-only closure になっている
- BEACH adapter では `sim.field_bc_mode="periodic2"` のとき `field_periodic_far_correction="none"` を後方互換として受理し、内部では `m2l_root_trunc` に正規化する。legacy な `m2l_root` も adapter 側で `m2l_root_trunc` に正規化される

この演算子は

$$
T_{\alpha\beta}^{\mathrm{root,trunc}}
\approx
(-1)^{|\beta|}
\sum_{(i,j)\in\Lambda_{\mathrm{outer}}}
D^{\alpha+\beta}\frac{1}{|R_{ij}|}
$$

という有限外側画像殻の打切り和で、exact Ewald residual ではありません。

### 8.3 `m2l_root_oracle`

`m2l_root_oracle` は build 時だけ exact periodic Ewald residual field を教師にして、root multipole から root local への演算子を proxy/check 点で数値同定するモードです。

- `periodic_image_layers = N`: runtime で explicit に残す近傍画像殻
- `periodic_ewald_alpha = alpha`: exact oracle に使う Ewald 分解パラメータ。`<= 0` なら自動決定
- `periodic_ewald_layers = L`: build-time oracle の real-space outer shell `N < max(|i|,|j|) <= N+L` と reciprocal cutoff `|h_1|, |h_2| <= L`
- build では exact periodic Ewald residual field を check points で評価し、field residual を least-squares fit して root local 演算子を作る
- runtime では `local(:, root) += T_root_oracle * multipole(:, root)` を加えるだけなので eval path に Ewald 本体は入らない
- tree 外 fallback では direct sum に exact periodic correction を直接足して、target box 外でも periodic residual を落とさない
- fit は potential ではなく field を使い、local の定数 potential mode は 0 に固定する

## 9. 計算量の見方

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

## 10. 現行実装の制約

この FMM コアは汎用 kernel FMM ではありません。

- kernel は Coulomb 固定
- simulator adapter の既定次数は `order = 4`
- source 座標は `build_plan` 後に不変とみなす
- 対応境界は `free` と `periodic2`
- `periodic2` は正確に 2 周期軸が必要
- far correction は `none`, `m2l_root_trunc`, `m2l_root_oracle`
- `m2l_root` は BEACH adapter が `m2l_root_trunc` に正規化する legacy 入力
- `eval_point(s)` の返り値には `k_coulomb` を含めない

## 11. 実装との対応

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
  `has_valid_target_box`, `use_periodic2_m2l_root_trunc`, `use_periodic2_m2l_root_oracle`,
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

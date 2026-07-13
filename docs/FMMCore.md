title: Coulomb FMM コア詳細

Lang: [日本語](FMMCore.md) | [English](FMMCore.en.md)

# Coulomb FMM コア詳細

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

BEACH の field solver adapter は、source model に応じて異なる幾何をコアへ渡します。

| source model | plan 構築 | `src_q(i)` の意味 |
|---|---|---|
| `point` | 要素重心を `src_pos` として `build_plan` に渡す | 重心に置く点電荷 |
| `triangle_p0` | 3 頂点を `build_panel_plan` に渡す | 三角形全体の総電荷。面密度は `src_q(i)/area(i)` |

- 初期化時は `build_plan` または `build_panel_plan` の直後に `update_state` を実行します。
- その後の refresh では、メッシュ幾何が変わらない通常運用を前提に既存 `plan` を再利用し、`src_q` 更新として `update_state` だけを呼びます。
- `build_plan` と legacy tree metadata の同期をやり直すのは、plan 未構築時・source 数変更時・要素数 0 件で plan/state を破棄したときだけです。

### 3. データ構造

#### 3.1 `fmm_options_type`

主な内部オプション:

- `theta`: well-separated 判定用パラメータ
- `leaf_max`: source octree の葉に許す最大 source 数
- `order`: Cartesian 展開次数
- `softening`: `point` source の softened Coulomb kernel に使う `epsilon`。`triangle_p0` では 0 が必須
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

#### 4.1 source kernel

コアには `point` と `triangle_p0` の 2 つの source kernel があります。

##### point source

`point` は要素重心に置いた点電荷を、softening 付き Coulomb kernel で評価します。

$$
G_\epsilon(\mathbf{r}) = \frac{1}{\sqrt{\lVert\mathbf{r}\rVert^2 + \epsilon^2}}
$$

$$
\phi(\mathbf{x}) = \sum_j q_j \, G_\epsilon(\mathbf{x} - \mathbf{x}_j)
$$

$$
\mathbf{E}(\mathbf{x}) = - \nabla \phi(\mathbf{x})
$$

近傍 direct 和と遠方展開は、同じ $G_\epsilon$ に基づきます。

##### P0 triangle source

`triangle_p0` は、`q_i` を三角形 $T_i$ の総電荷、$A_i$ をその面積とし、
$\sigma_i=q_i/A_i$ を三角形上の一定面密度として扱います。

$$
\phi_i(\mathbf{x}) =
\frac{q_i}{A_i}\int_{T_i}
\frac{1}{\lVert\mathbf{x}-\mathbf{y}\rVert}\,dA_{\mathbf{y}}
$$

近傍 direct 評価は、辺の対数項と立体角を使う解析的 P0 panel kernel です。
遠方 P2M には、tree node の中心 $\mathbf{c}$ に対する三角形上の monomial の
厳密面積平均を使います。

$$
M_{i,\alpha} =
q_i\left\langle(\mathbf{y}-\mathbf{c})^\alpha\right\rangle_{T_i}
= \frac{q_i}{A_i}\int_{T_i}
(\mathbf{y}-\mathbf{c})^\alpha\,dA_{\mathbf{y}}
$$

したがって area weighting は panel 積分と P2M 基底に含まれています。
`q_i` はすでに要素総電荷なので、`q_i` に $A_i$ をもう一度掛けません。
以後の M2M/M2L/L2L は、この panel moment に対する非 softening の
Coulomb/Laplace 展開です。`triangle_p0` では `softening=0` を強制し、
softened point source へ fallback しません。

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
(r_s + r_t)^2 < \theta_{\mathrm{eff}}^2 \, \lVert\mathbf{d}\rVert^2
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

Ewald2Pはruntime particle kernelではなく、`m2l_root_oracle`または`cached_kneq0`を作る
build-time teacherです。`triangle_p0`のcached経路でもteacherはproxy point chargeへ適用され、
結果をroot multipoleからlocal展開へのoperatorとしてfitします。実三角形の近傍評価は引き続き
解析panel kernel、遠方source表現はtriangle-averaged P2Mです。

$\alpha$はreal-spaceとreciprocal-spaceの収束を配分する数値パラメータで、Debye遮蔽を表しません。
直感的な分割とruntimeまでの流れは
[periodic2 zero modeと外部プラズマ](PeriodicZeroModeOuterPlasma.md#22-ewald2p-referenceとは)を参照してください。

##### 8.2.1 記法

周期軸を `a_1, a_2`、開放軸を `f` とします。
周期長、セル面積、画像集合、逆格子集合を次のように置きます。

$$
L_1 = \operatorname{periodic\_len}(1),\qquad
L_2 = \operatorname{periodic\_len}(2),\qquad
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

と書けます。ソース位置を $\mathbf s$、評価点を $\mathbf r$ とし、

$$
\mathbf R_{ij} = \mathbf r - \mathbf s - \mathbf L_{ij},\qquad
R_{ij} = \lVert\mathbf R_{ij}\rVert,\qquad
z = (\mathbf r - \mathbf s)\cdot \mathbf e_f
$$

を導入します。以下では $\alpha =$ `field_periodic_ewald_alpha`、$\epsilon =$ `softening` とします。

##### 8.2.2 実空間項

`add_screened_point_charge` が実装している screened Coulomb field は

$$
\mathbf E_\alpha(\mathbf R) =
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
\mathbf E_\epsilon(\mathbf R) =
q\,\frac{\mathbf R}{(R^2+\epsilon^2)^{3/2}}
$$

で、通常の runtime direct path と同じ softening を使います。

実装上の real-space 補正は

$$
\mathbf E_{\mathrm{real}} =
\sum_{(i,j)\in\mathcal I_{N+L}} \mathbf E_\alpha(\mathbf R_{ij})
{}- \sum_{(i,j)\in\mathcal I_N} \mathbf E_\epsilon(\mathbf R_{ij})
$$

です。実装では `r2 <= tiny(1.0d0)` の項はスキップするため、self interaction は入らない扱いです。`add_periodic2_exact_ewald_correction_single_source` に direct fallback の $\sum_{(i,j)\in\mathcal I_N}\mathbf E_\epsilon$ を足すと、inner image の softened 部分が打ち消され、outer shell 側は screened 形に置き換わります。

##### 8.2.3 逆空間項

`add_exact_periodic2_reciprocal_space_correction` が使う逆空間項は、$(m,n)\neq(0,0)$ に対して

$$
\theta_{mn} = \mathbf k_{mn}\cdot(\mathbf r-\mathbf s),\qquad
k_{mn} = \lVert\mathbf k_{mn}\rVert
$$

$$
G^\pm_{mn}(z) =
e^{\pm k_{mn} z}\operatorname{erfc}\!\left(\frac{k_{mn}}{2\alpha}\pm \alpha z\right)
$$

を定義すると

$$
\mathbf E_{\mathrm{rec}} =
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
\mathbf E_0 =
q\,\frac{2\pi}{A}\operatorname{erf}(\alpha z)\,\mathbf e_f
$$

です。single-source の oracle では `k=0` の電場寄与としてこの形を保持します。

##### 8.2.5 実装される補正

以上をまとめると、`add_periodic2_exact_ewald_correction_single_source` が 1 粒子分に加える補正は

$$
\mathbf E_{\mathrm{corr}} =
\mathbf E_{\mathrm{real}}
{}+ \mathbf E_{\mathrm{rec}}
{}+ \mathbf E_0
$$

です。`add_periodic2_exact_ewald_correction_all_sources` はまずこれを全ソースに対して総和します。

##### 8.2.6 `charged_walls` total-charge 補正

非中性 slab の `charged_walls` closure に合わせて、`add_periodic2_exact_ewald_correction_all_sources` は全ソース和のあとに total-charge 補正

$$
\mathbf E_{\mathrm{walls}}(z) =
\begin{cases}
\frac{2\pi Q_{\mathrm{tot}}}{A}\,\mathbf e_f, & z < z_{\mathrm{low}}, \\
0, & z_{\mathrm{low}} \le z \le z_{\mathrm{high}}, \\
-\frac{2\pi Q_{\mathrm{tot}}}{A}\,\mathbf e_f, & z > z_{\mathrm{high}}
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
\mathbf E_{\mathrm{fallback}} =
\sum_{(i,j)\in\mathcal I_N} \mathbf E_\epsilon(\mathbf R_{ij})
{}+
\mathbf E_{\mathrm{corr}}
{}+
\mathbf E_{\mathrm{walls}}
$$

です。`m2l_root_oracle` の build-time fit では check points が target box 内にあるため `\mathbf E_{\mathrm{walls}} = 0` となり、teacher には single-source の `\mathbf E_{\mathrm{corr}}` だけが使われます。`periodic_root_operator` 側では定数 potential mode を使わないため、monopole column は 0 に固定されます。

##### 8.2.7 `m2l_root_oracle`

`m2l_root_oracle` は、この Ewald2P 補正を teacher にして root multipole から root local への演算子を proxy/check 点で fit する明示 opt-in の高コスト診断モードです。無限周期の通常運用には `cached_kneq0` を使い、有限画像との互換性のため既定値だけは `none` に維持します。

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
- far correction は `none`（既定）, `auto`, `m2l_root_oracle`, `cached_kneq0`（`auto` は `none`、root oracle は診断、cached backend は production 非零モード）
- `eval_point(s)` の返り値には `k_coulomb` を含めない

### 10.1 cached periodic nonzero operator

#### 何を高速化するoperatorか

`periodic2`では、primary cellの外に同じ電荷分布がx/y方向へ無限に続きます。
通常のFMMは近傍画像`[-N,N]^2`までを高速に評価できますが、その外側にある無限個の画像が作る
滑らかな遠方場は残ります。全particle評価でEwald和を実行する代わりに、`cached_kneq0`は
この**遠方の差分だけ**をFMM local展開へ変換する線形operatorとして事前計算します。

| 入力 | cached operator | 出力 |
| --- | --- | --- |
| 現在電荷から作ったroot multipole | geometry固定の行列を適用 | target anchorごとの遠方local展開 |

cacheが保存するのは電場値、粒子位置、電荷履歴ではありません。固定geometryに対する
「source multipoleを遠方local展開へ写す行列」なので、batchごとに電荷が変わっても再利用できます。

#### 1回のfield評価で何を足すか

処理を順番に書くと次のようになります。

| 順序 | 成分 | 目的 |
| ---: | --- | --- |
| 1 | primary cell + 有限近傍画像 | singular/near fieldを通常のFMMとdirectで評価 |
| 2 | cached Ewald residual | 有限画像の外側にある滑らかな無限周期遠方場を補う |
| 3 | cached teacherに含まれた対称`k=0`を減算 | 非零モードbackendを`k!=0`だけにする |
| 4 | snapshotが物理的`k=0`を加算 | `symmetric_vacuum`、`e_bottom_zero`、outer plasmaを反映 |

したがって、`cached_kneq0`単体は全電場ではありません。3までが非零モードbackendの責任、
4は`electrostatic_snapshot`の責任です。`exclude_k0`は平均場を捨てる設定ではなく、
平均場を別のboundary providerが1回だけ加えるための重複防止規則です。

#### 数式との対応

上の1--3をまとめたruntime非零モードkernelは

$$
K_{k\ne0}
= K_\mathrm{shell}(N)
+ R_\mathrm{Ewald}^{\mathrm{full}}
- K_0^\mathrm{sym}
$$

| 項 | 作成時期 | 役割 |
| --- | --- | --- |
| $K_\mathrm{shell}(N)$ | 通常のFMM plan/runtime | primary cellと有限近傍画像 |
| $R_\mathrm{Ewald}^{\mathrm{full}}$ | cache cold build | full-periodic Ewald解と有限画像和の滑らかな差 |
| $K_0^\mathrm{sym}$ | charge state更新時 | cached full-periodic kernelから対称$k=0$を除く |

$K_0^\mathrm{sym}$ はsource高さの区分多項式累積stateを二分探索し、1評価あたり $O(\log n)$ です。
最終的なsurface fieldは

$$
K_\mathrm{surface}=K_{k\ne0}+K_0^\mathrm{physical}
$$

です。$K_0^\mathrm{physical}$のtriangle-height積分、lower boundary closure、outer-plasma接続は
[periodic2 zero modeと外部プラズマ](PeriodicZeroModeOuterPlasma.md)で説明します。

field fitとpotentialの定数modeは単位が異なるため、同じleast-squares列へ混ぜません。
potential gaugeは平均residualから別に固定します。

#### cache lifecycle

| 段階 | 動作 |
| --- | --- |
| identity作成 | geometry、target topology、order、周期長、画像層、generator/build versionからfingerprintを作る |
| warm read | version、fingerprint、shape、checksumが一致したoperatorだけを受理する |
| miss/corruption | filesystem lockを取得し、operatorを再生成する |
| publish | 同じdirectoryの`.tmp`をcloseしてからatomic renameする |
| checkpoint | operator本体は保存しない。再生成可能なcacheとして扱う |

この手順により、別jobは同じlock fileで直列化され、readerは部分書込みをcache hitとして受理しません。

#### MPI/OpenMPによるcold build

| 単位 | 担当 |
| --- | --- |
| cache I/Oとlock | MPI rank 0のみ |
| target operator slice | MPI rankへ均等分配。rank間差は最大1 target |
| target内のproxy列 | OpenMPで並列評価 |
| 正則化QR | targetごとに1回作り、全proxy RHSで再利用 |
| operator集約 | `MPI_Allreduce(SUM)`で全rankへ構築 |

warm runのfield evaluationとcharge refreshには、all-source Ewald和もoperator再fitもありません。

#### cold buildとwarm runの違い

| 経路 | Ewald teacher | QR fit | particle evaluation |
| --- | --- | --- | --- |
| 初回cache miss | 実行 | 実行してoperatorをpublish | build完了後に開始 |
| cache hit | 実行しない | 実行しない | 読み込んだoperatorを使用 |
| batch charge refresh | 実行しない | 実行しない | 新しいmultipoleへ同じoperatorを適用 |

「cold buildが重い」のは無限周期場を毎batch計算しているからではなく、再利用可能な行列を
初回だけ生成しているためです。warm runのhot pathには全source Ewald和は入りません。

#### SysA測定値

測定条件は2026-07-12の旧レゴリス入力、order 4、64 target、280 proxy、840 checkです。
時間はcache公開までではなく、表に明記した範囲の実測値です。

| 構成 | 実測時間 | 範囲 |
| --- | ---: | --- |
| 旧root-only、1 rank x 1 thread | 31分24秒 | cold operator build |
| QR再利用後、1 rank x 1 thread | 約25分45秒 | operator公開まで |
| 1 rank x 112 threads | 47.0秒 | cache prime + batch 1 |
| 2 ranks x 112 threads | 36.7秒 | cache prime + batch 1 |
| 4 ranks x 112 threads | 31.5秒 | cache prime + batch 1 |
| 6 ranks x 112 threads | 30.3秒 | cache prime + batch 1 |

全並列構成のchecksumは一致し、旧operatorとの差はFrobenius相対値 `1.73e-15` でした。
これは上記fixtureの測定結果であり、任意のgeometryに対する時間保証ではありません。

#### 運用指針

| 状況 | 推奨 |
| --- | --- |
| 専用cache prime | 1 rank x 112 threadsをcore効率とqueue footprintの基準にする |
| production allocationが既にある | 既存rankで生成してよい。6 x 112なら同fixtureで約30秒 |
| cold buildだけのために増員 | 4--6 rankの追加利得は小さい |
| 1 core | operator公開後の粒子batchでメモリ不足となった測定例があり、運用対象外 |
| warm cache | fingerprintとchecksumを確認し、そのまま再利用 |

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

title: Coulomb FMMコア詳細

Lang: [日本語](FMMCore.md) | [English](FMMCore.en.md)

# Coulomb FMMコア詳細

Coulomb FMMコアは、source geometryから作る`plan`と、電荷から作る`state`を分離し、
多重極展開と近傍Direct和から電場を評価します。ここでは、現行Fortran実装の低水準開発API、データ構造、
展開式、periodic2補正を、処理順に沿って説明します。

通常のFMMを利用するための数式と計算フローは[FMM](FMM.html)、periodic root operatorの構成と運用は
[periodic2遠方補正](PeriodicFarCorrection.html)にまとめています。このページはFortran内部配列と実装手順の
詳細を扱います。

APIの詳細は[`bem_coulomb_fmm_core` module page](../module/bem_coulomb_fmm_core.html)でも確認できます。

- 低水準開発 API / 境界: `src/physics/field_solver/fmm/api/`
- 内部共通実装: `src/physics/field_solver/fmm/internal/common/`
- tree / plan 実装: `src/physics/field_solver/fmm/internal/tree/`
- state / eval 実装: `src/physics/field_solver/fmm/internal/runtime/`
- periodic2 実装: `src/physics/field_solver/fmm/internal/periodic/`

コアの内部APIはsimulatorに依存せず、`mesh_type`や`sim_config`を直接`use`しません。
BEACHのfield solver adapterが、simulator側のデータを変換してこのコアを呼び出します。

## 1. 目的

このFMMコアは、固定されたsource geometryと可変電荷`src_q(n)`から、多数の評価点における
Coulomb電場を計算します。低水準の検証・periodic operator構築用にはmonopole位置
`src_pos(3,n)`を受けるgeneric planがあり、BEACHの表面電荷には三角形3頂点を受けるpanel planを
使います。generic planはBEACHの表面source modelではなく、TOMLやC/Python APIからは選べません。

計算の前提と責務は次のとおりです。

- kernel は 3D Coulomb のみ
- source 幾何と電荷更新を分離する
- `free` と `periodic2` のみを対象にする
- 近傍 direct 和もコア内部に含める
- simulator からは配列 API だけが見えるようにする

## 2. 低水準開発 API

### 2.1 Fortran API

コアが提供する主な手続きは次のとおりです。

```fortran
call build_plan(plan, src_pos, options)             ! generic monopole plan
call build_panel_plan(plan, v0, v1, v2, options)   ! BEACH surface plan
call update_state(plan, state, src_q)
call eval_points(plan, state, target_pos, e)
call eval_point(plan, state, r, e)
```

入力・出力の意味は次の通りです。

- `src_pos(3,n)`:
  低水準generic planのsource点座標。`build_plan`後は固定とみなす。
- `v0(3,n)`, `v1(3,n)`, `v2(3,n)`:
  panel planの三角形頂点。`build_panel_plan`後は固定とみなす。
- `src_q(n)`:
  source 点の電荷。`update_state` ごとに更新できる。
- `target_pos(3,m)` または `r(3)`:
  評価点。
- `e(3,m)` または `e(3)`:
  電場ベクトル。

呼び出し側では、次の点に注意が必要です。

- コアが返す電場には`k_coulomb`を掛けていません。BEACHのadapterが最後に掛けます。
- `build_plan` / `build_panel_plan`は幾何依存処理、`update_state`は電荷依存処理です。
- `eval_point(s)` は `plan` と `state` が ready な前提です。

### 2.2 C ABI / Python 連携

`src/physics/field_solver/bem_field_kernel_c.f90`は、Fortran APIを`iso_c_binding`の
opaque handle APIとして公開します。`make build-kernel`を実行すると、共有ライブラリ
`build/libbeach_field_kernel.so`が生成されます。

主な C ABI:

```text
beach_kernel_get_abi_version(major, minor)
beach_kernel_get_build_info(buffer, capacity, length)
beach_kernel_create(handle)
beach_kernel_destroy(handle)
beach_kernel_build(handle, vertex0_xyz, vertex1_xyz, vertex2_xyz, options...)
beach_kernel_update_charges(handle, src_q)
beach_kernel_eval_e(handle, target_pos, e)
beach_kernel_eval_phi(handle, target_pos, phi)
beach_kernel_force_on_charges(handle, target_pos, target_q, origin, force, torque)
```

公開ヘッダは Python package 内の `beach/include/beach_field_kernel.h` にあります。現行 ABI は
`2.0` です。C の呼び出し側は、ほかの関数を使う前に `beach_kernel_get_abi_version` を呼び、
major が一致し、library の minor が必要な minor 以上であることを確認してください。座標と
vector の配列は `values[3 * point_index + component]` の順で格納します。status code、
periodic far-correction code、handle の所有権は公開ヘッダに定義しています。

Pythonでは、`beach.fortran_results.kernel.FieldKernel`がこのABIを`ctypes`で呼び出します。
Python wrapper はversion問い合わせを提供するlibraryの互換性を読込時に検査します。
問い合わせsymbolがない旧libraryも拒否し、ABI v2以上を明示したlibraryだけを受理します。
`calc_object_forces_kernel`は、対象object自身のsource電荷をゼロにして
`sum(q_i E_not_self(r_i))`を評価します。これにより、自己力を除外しながら、
`periodic2 + cached_kneq0`を含むfield kernelをそのまま利用できます。

`Beach.scene()` / `BeachScene`は、Python側でobjectの剛体移動や回転を一時的に適用し、
変換後の三角形3頂点を同じABIへ渡します。剛体変換の補助処理には既定でNumPyを使い、
任意依存のNumba backendも選べます。電場の評価は、どちらの場合もFortran kernelが担当します。

### 2.3 BEACH adapterでの使い方

BEACH の field solver adapter は、各三角形の 3 頂点を `build_panel_plan` に渡します。
`src_q(i)` は三角形全体の総電荷で、面密度は `src_q(i)/area(i)` です。

- 初期化時は `build_panel_plan` の直後に `update_state` を実行します。
- その後のrefreshでは、mesh geometryが変わらない限り既存の`plan`を再利用し、`src_q`を渡して
  `update_state`だけを呼びます。
- plan を再構築するのは、planが未構築の場合、source数が変わった場合、
  または要素数0としてplan/stateを破棄した場合だけです。

## 3. データ構造

### 3.1 `fmm_options_type`

主な内部オプション:

- `theta`: well-separated 判定用パラメータ
- `leaf_max`: source octree の葉に許す最大 source 数
- `order`: Cartesian 展開次数（1以上）
- `softening`: 低水準の汎用monopole planだけが使う内部値。BEACH adapterは常に0を渡す
- `use_periodic2`: 2 周期軸モードの有効化
- `periodic_axes(2)`, `periodic_len(2)`: 周期軸と周期長
- `periodic_image_layers`: 近傍画像和の層数 `N`
- `periodic_far_correction`: coreが受ける値は`auto`, `none`, `cached_kneq0`。`periodic2`有効時の
  `auto`は互換性のため`none`へ正規化される
- `periodic_ewald_alpha`, `periodic_ewald_layers`: `cached_kneq0`のbuild-time Ewald fitで使う
  分解パラメータと打切り深さ
- `target_box_min/max`: dual-target tree を作るときの box

BEACH の adapter は現状 `order = 4` を使いますが、コア自体は1以上の可変次数を受けられます。
`order = 0` は電場のfar/local展開を持てないため、plan構築時に拒否します。
`periodic2` の `auto` は `none` に正規化されます。`cached_kneq0` は遠方補正を明示的に有効化します。

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

`multipole`はsource tree nodeごとの多重極係数、`local`はtarget tree nodeごとの局所展開係数です。
`*_active`はzero-nodeの計算を省略するための0/1 flagです。

## 4. 数学的定義

### 4.1 source kernel

BEACH runtimeのsource kernelはP0 triangleに固定されています。低水準FMM coreには
`build_plan`による汎用monopole planも内部API・回帰試験用として残っていますが、
BEACH field solver、公開C ABI、Python APIから選択する経路はありません。

runtimeの`q_i`を三角形$T_i$の総電荷、$A_i$をその面積とし、
$\sigma_i=q_i/A_i$ を三角形上の一定面密度として扱います。

$$
\phi_i(\mathbf{x}) =
\frac{q_i}{A_i}\int_{T_i}
\frac{1}{\lVert\mathbf{x}-\mathbf{y}\rVert}\,dA_{\mathbf{y}}
$$

近傍Direct評価には、辺の対数項と立体角を使う解析的P0 panel kernelを使います。
遠方P2Mには、tree nodeの中心$\mathbf{c}$に対する三角形上のmonomialの厳密な面積平均を使います。

$$
M_{i,\alpha} =
q_i\left\langle(\mathbf{y}-\mathbf{c})^\alpha\right\rangle_{T_i}
= \frac{q_i}{A_i}\int_{T_i}
(\mathbf{y}-\mathbf{c})^\alpha\,dA_{\mathbf{y}}
$$

area weightingはpanel積分とP2M基底に含まれています。`q_i`は要素の総電荷なので、
$A_i$を重ねて掛ける必要はありません。以後のM2M/M2L/L2Lには、このpanel momentに対する
Coulomb/Laplace展開を使います。近傍評価は解析的 panel kernel で処理します。

内部`build_plan`は、必要な低水準試験でだけ
$G_\epsilon(\mathbf r)=1/\sqrt{\lVert\mathbf r\rVert^2+\epsilon^2}$を使います。
これは削除済みの公開source modelを復活させるものではなく、`beach.toml`から設定できません。

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

node centerを$c$とすると、leaf nodeのmultipole係数は

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
D^{\alpha+\beta} G(R)
$$

で更新します。

ここで $D^\gamma$ は multi-index 微分です。
現行実装では、$D^{\alpha+\beta} G(R)$をpairごとに前計算し、
`m2l_deriv(:, pair)`へ保存します。

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

### 5.3 near/farとM2L pair cache

各target leafについてsource treeを再帰的に走査し、near nodeとfar nodeへ分類します。

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

その後、dual-tree再帰によってM2L pair cacheを作り、target nodeごとの索引配列を準備します。

### 5.4 build時の前計算

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

この前計算により、`update_state`では主に電荷に依存する加算だけを実行します。

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

`update_state`はlegacy実装のrefreshに相当します。source座標は固定され、`src_q`だけが変わることを
前提とします。

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

### 6.2 OpenMP並列化

現行実装では、次の処理をOpenMPで並列化しています。

- `update_state`全体を1つのparallel regionで囲み、その内側で`src_q`のコピーとactive flagの初期化を行う
- `P2M`: source leaf ごとのループ
- `M2M`: 同一 depth の node ループ
- `M2L`: target node ごとのループ
- `L2L`: 同一 depth の node ループ
- `build_plan` 時の translation / M2L 微分前計算

各ループは1 nodeを1 threadが担当しやすい構造とし、共有配列の更新をnode単位で独立させています。

### 6.3 実装上の最適化

`update_state` では次の無駄を避けています。

- $\beta - \alpha$ の multi-index 差分を毎回計算しない
- 親子 center 差分のべき乗を毎回作り直さない
- `P2M` の monomial 基底を source ごとに build 時に前計算する
- `M2M/L2L` の有効な `(alpha, delta)` 項だけを圧縮して持ち、無効項を走査しない
- `M2L`ではsource nodeのactive flagを調べ、zero-nodeのpairを早期にskipする
- `M2L`ではtarget node列を繰り返し更新せず、thread-localな`local_acc`へ蓄積してから書き戻す
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
    return

  evaluate local expansion at leaf center
  add near direct interactions
  root local already carries periodic root correction when enabled
```

### 7.1 leafの特定

- `periodic2` では評価点を target box 内へ wrap してから探索する
- target tree があるときは target tree の葉を使う
- target tree が無いときは source tree の葉を使う
- leaf lookupに失敗した場合、またはleafをtreeのleaf slotへ写像できない場合はDirect fallbackに入る

### 7.2 近傍Direct和

near listに入ったsource indexについてDirect和を計算します。`periodic2`では、
`[-N, N] x [-N, N]`の画像シフトも陽に加算します。fallbackでも同じDirect kernelを使います。

### 7.3 box外fallback

dual-target treeを使う場合、target boxの外にある評価点にはtarget leafがありません。この場合は、
全sourceのDirect和へfallbackします。`cached_kneq0`は固定target topologyを前提とするため、
target box外の評価をrejectします。

### 7.4 root補正の位置

`cached_kneq0`のroot補正は、`update_state`がtarget anchorのlocal展開へ注入します。
通常のleaf評価ではroot補正を再計算せず、`state`に保存されたlocal展開をそのまま使います。

## 8. `periodic2` と遠方補正

この節はFMM core内部から見た数式とfallbackを記録します。設定の選び方、operator fit、cache lifecycle、
`k=0` ownershipを通した説明は[periodic2遠方補正](PeriodicFarCorrection.html)に分離しています。

### 8.1 `periodic2`

`periodic2`では、2軸を周期境界、残りの1軸を開放境界とします。

近傍画像和は

$$
i, j \in [-N, N]
$$

の有限画像を陽に足します。

M2L でも同じ画像シフト集合を使い、各 pair の derivative を画像和で前計算します。

### 8.2 `periodic2` Ewald（Ewald2P）補正

`bem_coulomb_fmm_periodic_ewald.f90`は、2周期・1開放のCoulomb fieldに対するEwald型の補正を
実装します。ここでいう`exact`は、コードが実際に評価する有限和を指し、理論上の無限和そのものではありません。
real-spaceとreciprocal-spaceの打切り深さは、`field_periodic_image_layers = N`と
`field_periodic_ewald_layers = L`で決まります。この有限和をbuild-time oracleとして使います。

Ewald2Pはruntime particle kernelではなく、`cached_kneq0`を作るための
build-time teacherです。teacherはproxy monopoleに適用されます。
その結果を、root multipoleからlocal展開へのoperatorとしてfitします。実三角形の近傍評価には解析panel kernel、
遠方source表現にはtriangle-averaged P2Mを引き続き使います。

$\alpha$はreal-spaceとreciprocal-spaceの収束を配分する数値パラメータで、Debye遮蔽を表しません。
この分割をruntimeの場評価へ組み込む流れは、periodic2静電場の説明にまとめています。[<sup>1</sup>](PeriodicElectrostatics.md#ewald2p-teacher)

#### 8.2.1 記法

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

を導入します。以下では $\alpha =$ `field_periodic_ewald_alpha` とします。

#### 8.2.2 実空間項

内部helperが実装しているscreened Coulomb fieldは

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

低水準Ewald helperが内側画像和から差し引くdirect fieldは

$$
\mathbf E_\epsilon(\mathbf R) =
q\,\frac{\mathbf R}{(R^2+\epsilon^2)^{3/2}}
$$

です。BEACHのP0 panel経路では$\epsilon=0$なので、
$\mathbf E_\epsilon=\mathbf E_0=q\mathbf R/R^3$です。

実装上の real-space 補正は

$$
\mathbf E_{\mathrm{real}} =
\sum_{(i,j)\in\mathcal I_{N+L}} \mathbf E_\alpha(\mathbf R_{ij})
{}- \sum_{(i,j)\in\mathcal I_N} \mathbf E_\epsilon(\mathbf R_{ij})
$$

です。実装では`r2 <= tiny(1.0d0)`の項をskipするため、self interactionは含まれません。
`add_periodic2_exact_ewald_correction_single_source`にDirect fallbackの
$\sum_{(i,j)\in\mathcal I_N}\mathbf E_\epsilon$を加えると、inner imageのdirect部分が打ち消され、
outer shell側がscreened形式に置き換わります。

#### 8.2.3 逆空間項

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

#### 8.2.4 `k=0`項

`add_exact_periodic2_k0_correction` が実装しているゼロモード補正は

$$
\mathbf E_0 =
q\,\frac{2\pi}{A}\operatorname{erf}(\alpha z)\,\mathbf e_f
$$

です。single-source の Ewald teacher では `k=0` の電場寄与としてこの形を保持します。

#### 8.2.5 実装される補正

以上をまとめると、`add_periodic2_exact_ewald_correction_single_source` が 1 粒子分に加える補正は

$$
\mathbf E_{\mathrm{corr}} =
\mathbf E_{\mathrm{real}}
{}+ \mathbf E_{\mathrm{rec}}
{}+ \mathbf E_0
$$

です。このsingle-source補正をproxy/check点で評価し、cached operatorを生成します。

`field_periodic_ewald_alpha` が `<= 0` の場合、`resolve_periodic2_ewald_alpha` は

$$
\alpha = \frac{1.2}{(N+1)\min(L_1,L_2)}
$$

を自動選択します。`min(L_1,L_2)\le 0` なら `alpha = 0` としてoperator生成を無効化します。
また内部では `kmax = max(1, field_periodic_ewald_layers)` として逆空間の有限和を作ります。

`cached_kneq0`のcold buildはcheck pointsでEwald residualの電場と電位を評価し、root multipoleから
target localへのoperatorをfitします。電場fitで決まらない定数potential係数はpotential residualから
別にfitし、versioned cacheへ保存します。`m2l_root_oracle` は削除済みで、指定時はrejectします。
無限周期のproduction計算には`cached_kneq0`を使います。

## 9. 計算量の見方

展開次数$p$を固定し、interaction listの長さに上限があると仮定すると、計算量の目安は次のとおりです。

- `build_plan`: $O(n \log n)$ に近い
- `update_state`: $O(n)$ に近い
- `eval_point`: $O(\log n + n_{\mathrm{near}} \, n_{\mathrm{img}}^2)$ に近い
- `eval_points`: 上記の点評価を target ごとに並列実行

実際の定数係数は、次の設定に強く依存します。

- `order`
- `theta`
- `leaf_max`
- `periodic_image_layers`
- `periodic_ewald_layers`
- target tree の有無

## 10. 現行実装の制約

このFMMコアは、次の条件に特化した実装です。

- kernel は Coulomb 固定
- simulator adapter の既定次数は `order = 4`
- source 座標は `build_plan` 後に不変とみなす
- 対応境界は `free` と `periodic2`
- `periodic2` は正確に 2 周期軸が必要
- far correctionは`none`（既定）、`auto`、`cached_kneq0`。
  `auto`は`none`として動作し、cached backendはproductionの非零modeに使う
- `eval_point(s)` の返り値には `k_coulomb` を含めない

### 10.1 cached periodic nonzero operator

#### 何を高速化するoperatorか

`periodic2`では、primary cellの外に同じ電荷分布がx/y方向へ無限に続きます。
通常のFMMは近傍画像`[-N,N]^2`を高速に評価できますが、その外側にある無限個の画像が作る
滑らかな遠方場は含みません。`cached_kneq0`は、この**遠方の差分だけ**をFMM local展開へ変換する
線形operatorとして事前計算します。これにより、particleを評価するたびにEwald和を計算する必要がなくなります。

| 入力 | cached operator | 出力 |
| --- | --- | --- |
| 現在の電荷から作ったroot multipole | geometry固定の行列を適用 | target anchorごとの遠方local展開 |

cacheが保存するのは、固定geometryに対してsource multipoleを遠方local展開へ写す行列です。
電場値、粒子位置、電荷履歴は保存しません。このため、batchごとに電荷が変わっても同じ行列を再利用できます。

#### 1回のfield評価で何を足すか

1回のfield評価では、次の順に各成分を加算します。

| 順序 | 成分 | 目的 |
| ---: | --- | --- |
| 1 | primary cell + 有限近傍画像 | singular/near fieldを通常のFMMとdirectで評価 |
| 2 | cached Ewald residual | 有限画像の外側にある滑らかな無限周期遠方場を補う |
| 3 | cached teacherに含まれた対称`k=0`を減算 | 非零モードbackendを`k!=0`だけにする |
| 4 | 場の合成時に物理的`k=0`を加算 | `symmetric_vacuum`または`e_bottom_zero`を反映 |

`cached_kneq0`単体では全電場を構成しません。手順3までを非零mode backendが担当し、手順4を
`electrostatic_snapshot`が担当します。`exclude_k0`は平均場を除外する設定ではありません。
別のboundary providerが平均場を1回だけ加えるための重複防止規則です。

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

$K_0^\mathrm{sym}$は、source高さから作った区分多項式の累積stateを二分探索して評価します。
1評価あたりの計算量は$O(\log n)$です。
最終的なsurface fieldは

$$
K_\mathrm{surface}=K_{k\ne0}+K_0^\mathrm{physical}
$$

です。$K_0^\mathrm{physical}$のtriangle-height積分と下側境界条件は
[periodic2静電場](PeriodicElectrostatics.md)で説明します。

field fitとpotentialの定数modeは単位が異なるため、同じleast-squaresの列には含めません。
potential gaugeは、平均residualから別途固定します。

#### cache lifecycle

| 段階 | 動作 |
| --- | --- |
| identity作成 | geometry、target topology、order、周期長、画像層、generator/build versionからfingerprintを作る |
| warm read | version、fingerprint、shape、checksumが一致したoperatorだけを受理する |
| miss/corruption | filesystem lockを取得し、operatorを再生成する |
| publish | 同じdirectoryの`.tmp`をcloseしてからatomic renameする |
| checkpoint | operator本体は保存しない。再生成可能なcacheとして扱う |

同時に起動したjobは、同じlock fileによってcache生成を直列化します。readerは、書込み途中のfileを
cache hitとして受理しません。

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

cold buildでは、再利用可能なoperator行列を初回だけ生成します。batchごとに無限周期場を再構築するわけでは
ありません。warm runのhot pathでは、全sourceのEwald和を計算しません。

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

全並列構成でchecksumが一致し、旧operatorとの差はFrobenius相対値`1.73e-15`でした。
実測時間は上記fixtureに対する値であり、geometryによって変わります。

#### 運用指針

| 状況 | 推奨 |
| --- | --- |
| 専用cache prime | 1 rank x 112 threadsをcore効率とqueue footprintの基準にする |
| production allocationが既にある | 既存rankで生成してよい。6 x 112なら同fixtureで約30秒 |
| cold buildだけのために増員 | 4--6 rankの追加利得は小さい |
| 1 core | operator公開後の粒子batchでメモリ不足となった測定例があり、運用対象外 |
| warm cache | fingerprintとchecksumを確認し、そのまま再利用 |

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
  `build_plan`, `build_panel_plan`
  （`src/physics/field_solver/fmm/internal/tree/bem_coulomb_fmm_plan_ops.f90`）
- charge refresh:
  `update_state`, `p2m_leaf_moments`, `m2m_upward_pass`, `m2l_accumulate`, `l2l_downward_pass`
  （`src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_state_ops.f90`）
- 評価:
  `eval_point`, `eval_points`
  （`src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90`）
- periodic2 補助:
  `has_valid_target_box`, `use_periodic2_cached_kneq0`,
  `use_periodic2_root_operator`, `build_periodic_shift_values`, `add_point_charge_images_field`,
  `wrap_periodic2_point`, `apply_periodic2_minimum_image`, `distance_to_source_bbox`,
  `distance_to_source_bbox_periodic`
  （`src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic.f90`）
- periodic2 Ewald/oracle:
  `resolve_periodic2_ewald_alpha`, `precompute_periodic2_ewald_data`,
  `add_periodic2_exact_ewald_correction_single_source`
  （`src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_ewald.f90`）
- periodic2 root operator:
  `precompute_periodic_root_operator`
  （`src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90`）
- BEACH adapter:
  `src/physics/field_solver/bem_field_solver_config.f90`,
  `src/physics/field_solver/bem_field_solver_tree.f90`,
  `src/physics/field_solver/bem_field_solver_eval.f90`

コアとBEACH adapterの責務は次のように分かれます。

- コア:
  幾何前処理、展開係数更新、近傍 direct、点評価
- BEACH adapter:
  `mesh_type`の三角形3頂点で`build_panel_plan`を構築し、`q_elem`を`src_q`へ渡し、
  最後に`k_coulomb`を掛ける

---

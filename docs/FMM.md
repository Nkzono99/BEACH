title: FMM

Lang: [日本語](FMM.md) | [English](FMM.en.md)

# FMM

Fast Multipole Method（FMM）は、遠方にある多数のsourceの寄与を展開係数へまとめます。固定geometryから
作った計算planを再利用し、batchごとに電荷に依存する係数だけを更新します。各粒子位置では遠方展開と
近傍Direct和を合成するため、多数の境界要素を繰り返し評価する計算のコストを削減できます。

| 向いている計算 | 先に確認すること |
| --- | --- |
| 要素数が多く、1 batchで多数の粒子stepを評価する | Directとの差、release buildでの実測時間 |
| triangle P0を大規模meshで使う | panel Directとの差とmesh細分化 |
| `periodic2`場を使う | image和、nonzero/zero mode、outer modelの構成 |

```toml
[sim]
field_solver = "fmm"
field_bc_mode = "free"
use_box = true
box_min = [-0.5, -0.5, -0.1]
box_max = [ 0.5,  0.5,  1.0]
```

`use_box=true`にすると、粒子が動く領域を覆うtarget treeが作られます。target treeの外にある評価点は、
全sourceのDirect和へfallbackします。

## 遠方展開と近傍Direct和を合成する

FMMはsourceとtargetをoctreeへ分け、近傍と遠方を分離します。

```text
source charge
   │ P2M: leafのmultipoleを作る
   ▼
source multipoles
   │ M2M: childからparentへ集約
   ▼
far interactions
   │ M2L: target nodeのlocal expansionへ変換
   ▼
target locals
   │ L2L: parentからchildへ伝播
   ▼
L2P at target + near Direct interactions
```

遠方相互作用はCartesian多重極・局所展開で近似し、near listのsourceは選択したpointまたはpanel kernelで
直接評価します。現行simulator adapterの展開次数は`order = 4`固定で、入力からは変更できません。

### Cartesian展開の記法

多重指数を$\alpha=(\alpha_x,\alpha_y,\alpha_z)$とし、次を使います。

$$
|\alpha|=\alpha_x+\alpha_y+\alpha_z,qquad
\alpha!=\alpha_x!\alpha_y!\alpha_z!,qquad
\mathbf r^\alpha=r_x^{\alpha_x}r_y^{\alpha_y}r_z^{\alpha_z}.
$$

point sourceのkernelは

$$
G_\epsilon(\mathbf r)=\frac{1}{\sqrt{\lVert\mathbf r\rVert^2+\epsilon^2}},qquad
\phi(\mathbf x)=\sum_j q_jG_\epsilon(\mathbf x-\mathbf x_j),qquad
\mathbf E=-\nabla\phi
$$

です。以下の係数にはCoulomb定数を含めず、BEACH adapterが評価後に単位系とともに掛けます。

### P2MとM2Mでsourceを集約する

leaf中心を$\mathbf c$とすると、point sourceのmultipole係数は

$$
M_\alpha(\mathbf c)=
\sum_{j\in\mathrm{leaf}}q_j
\frac{(\mathbf x_j-\mathbf c)^\alpha}{\alpha!}
$$

です。triangle P0では点位置のmonomialの代わりに、三角形上の厳密な面積平均を使います。

$$
M_{i,\alpha}=\frac{q_i}{A_i}\int_{T_i}
\frac{(\mathbf y-\mathbf c)^\alpha}{\alpha!}\,dA_{\mathbf y}
$$

子中心$\mathbf c_c$から親中心$\mathbf c_p$へは、$\mathbf d=\mathbf c_c-\mathbf c_p$として

$$
M_\beta(\mathbf c_p)=
\sum_{\alpha\le\beta}M_\alpha(\mathbf c_c)
\frac{\mathbf d^{\beta-\alpha}}{(\beta-\alpha)!}
$$

を全childについて足します。P2M basisとM2Mのshift monomialはgeometryだけで決まるため、plan構築時に
前計算します。

### M2Lでmultipoleをtarget localへ変換する

source node中心$\mathbf c_s$、target node中心$\mathbf c_t$、$\mathbf R=\mathbf c_t-\mathbf c_s$に対して、
well-separatedなnode pairの寄与は

$$
L_\alpha(\mathbf c_t)\mathrel{+}=
\sum_\beta(-1)^{|\beta|}M_\beta(\mathbf c_s)
D^{\alpha+\beta}G_\epsilon(\mathbf R)
$$

です。$D^\gamma$はmulti-index微分です。どのsource/target nodeを組み合わせるかはplanのM2L pair cacheに、
kernel微分$D^{\alpha+\beta}G_\epsilon(\mathbf R)$は`m2l_deriv`に保存されます。したがってbatchごとの
`update_state`では、現在のmultipole係数との積和が主な処理になります。

### L2LとL2Pで評価点まで伝える

親localを子中心へ移すときは、$\mathbf d=\mathbf c_c-\mathbf c_p$として

$$
L_\alpha(\mathbf c_c)\mathrel{+}=
\sum_{\gamma\ge\alpha}L_\gamma(\mathbf c_p)
\frac{\mathbf d^{\gamma-\alpha}}{(\gamma-\alpha)!}
$$

を使います。評価点$\mathbf x$が属するleaf中心を$\mathbf c_l$、$\mathbf r=\mathbf x-\mathbf c_l$とすると、

$$
\phi_\mathrm{far}(\mathbf x)=
\sum_\alpha L_\alpha(\mathbf c_l)\frac{\mathbf r^\alpha}{\alpha!},
$$

$$
E_{\mathrm{far},k}(\mathbf x)=
-\sum_\alpha L_{\alpha+\mathbf e_k}(\mathbf c_l)
\frac{\mathbf r^\alpha}{\alpha!}.
$$

最後に同じleafのnear listを元のpoint/panel kernelでDirect評価して加えます。M2Lだけで全相互作用を
評価しているわけではなく、FMMの結果は常に「far local展開 + near Direct和」です。

## geometryと電荷更新を分けて再利用する

固定meshで複数のbatchを計算するため、FMMは幾何に依存する`plan`と電荷に依存する`state`を分けます。

| データ | 主な内容 | 更新時期 |
| --- | --- | --- |
| plan | source tree、target tree、near/far list、P2M basis、translation operator | 初期化時。geometry、要素数、主要optionsが変われば再構築 |
| state | 現在の`q_elem`、multipole係数、local係数 | batch開始時に場を更新するとき |

`build_plan`はsource treeとinteraction listを作り、geometryだけで決まる量を前計算します。
`update_state`は現在の要素電荷からP2M、M2M、M2L、L2Lを実行します。各粒子位置では、属するtarget leafの
local expansionとnear Direct和だけを評価します。

同じbatchの途中で粒子ごとにstateを更新することはありません。表面への堆積電荷はbatch末尾でcommitされ、
次batchの場を更新するときにstateへ反映されます。

## 粒子が動く領域をtarget treeで覆う

FMM planはtarget側を次のどちらかで構成します。

- `use_box=false`: source treeのleafをtarget leafとしても使う
- `use_box=true`: `box_min`から`box_max`を覆う独立したtarget treeを作る

粒子の評価領域がsourceのbounding boxより広い場合、前者ではtree外の点がDirect fallbackになりやすくなります。
通常は粒子が到達し得る領域をboxで覆います。free境界でもboxを設定でき、box面の粒子境界条件とtarget treeの
範囲に同じ幾何を使えます。

target leafが見つからない点では、全sourceの寄与を直接加算します。このfallbackが頻発すると、計算量は
Direct相当になります。box外へ出る粒子が多い場合は、境界通過が場評価より先に処理されていることと、
boxが意図した領域を覆っていることを確認してください。

<a id="source-kernel"></a>

## source離散化をnearとfarで揃える

### 重心点電荷

point modeは各要素重心をsource位置とし、near interactionには`sim.softening`を含むpoint kernelを使います。
far interactionは同じ電荷のmultipole展開です。内部座標ではsofteningも$L_0$で割って正規化されます。

### 三角形上の一定面電荷

triangle P0 modeは三角形全体をsource geometryとして保持します。source treeのboxには、重心だけでなく
三頂点も含まれます。P2Mには三角形上のmonomialの厳密な面平均を使い、near interactionとtree外fallbackには
[Direct solver](DirectSolver.html#triangle-p0)と同じ解析panel積分を使います。

near interactionもpoint sourceへ置き換えません。このmodeには`softening=0`、非退化三角形、insulator表面、
解決済みvacuum sideが必要です。現行FMM coreのpanel電場は、評価点が厳密にpanel面上にある場合に
principal-value traceを使います。面上の場を直接比較する検証では、Directのvacuum-side traceとの違いを
考慮してください。粒子追跡では、panel面に達する前の軌道線分との交差として衝突を処理します。

設定例は[`examples/panel_fmm.toml`](../examples/panel_fmm.toml)にあります。

## 近傍と遠方の分け方で精度を調整する

source node半径$r_s$、target node半径$r_t$、中心間距離$d$に対して、概念的には

$$
(r_s+r_t)^2 < \theta^2 d^2
$$

を満たすnode対をwell-separatedと判定し、遠方展開で評価します。それ以外は木を細分化し、最終的に
near Directで評価します。

| key | 影響 |
| --- | --- |
| `tree_theta` | 小さいほどfar判定が厳しく、一般に高精度・低速 |
| `tree_leaf_max` | leafのsource数。near Direct量、木の深さ、memoryのバランスを変える |
| `field_normalization` | 座標の数値スケール。物理単位やモデルは変えない |
| `softening` | point kernelだけに適用。triangle P0では0が必須 |

`tree_theta`と`tree_leaf_max`を入力に書かなければ、明示`fmm`でも要素数に応じて
`(0.40, 12)`、`(0.50, 16)`、`(0.58, 20)`、`(0.65, 24)`の順に自動設定します。区切りは
`1500`、`10000`、`50000`要素です。これらは誤差保証ではなく初期値です。

展開次数は現在4固定です。FMM近似の感度は`tree_theta`、`tree_leaf_max`、Directとの比較で確認します。
source離散化の収束は、これとは別にmesh細分化で確認します。

## 電場・電位と履歴出力を評価する

通常のtargetでは、電場も電位もlocal expansionとnear Direct和から評価します。要素中心の
`potential_history.csv`では同じFMM評価後にpoint kernelの有限自己項を補います。triangle P0は解析panel自己積分を
near kernelに含めます。[<sup>1</sup>](DirectSolver.html#要素中心の電位出力)

`potential_history`を書き出す時点では、最新の要素電荷でstateをrefreshします。そのため履歴を有効にすると、
通常のbatch field更新に加えて、stateのrefreshと全要素targetの評価が発生する場合があります。

## periodic2の遠方補正はlocal展開へ接続する

`periodic2`でも通常のP2MからL2Pまでは変わりません。有限画像の外側にある滑らかな周期遠方場だけを、
root multipoleからtarget localへの追加operatorとして`M2L`後・`L2L`前に注入します。そのため実装はFMMの
`plan`、`state`、local展開を共有しますが、通常のtree M2Lとは異なる計算段階として分離されています。

`none`とproduction用`cached_kneq0`の違い、Ewald2P teacher、cache、`k=0`の
ownershipは[periodic2遠方補正](PeriodicFarCorrection.html)に分けました。場全体の成分構成は
[periodic2静電場](PeriodicElectrostatics.html)、外部領域との結合は[外部プラズマモデル](OuterPlasmaModels.html)で
説明します。

小規模のsplit referenceは、`field_solver="direct"`とpanel spectral backendを組み合わせる別経路です。
対応する構成は[periodic2無限周期＋outer plasma構成](InfinitePeriodicOuterConfiguration.html)にまとめています。

## 固定費と粒子評価時間を分けて測る

固定次数でinteraction listが過度に増えない場合、目安は次のとおりです。

| 処理 | 目安 |
| --- | --- |
| plan構築 | $O(N\log N)$に近い。一度だけだがtranslationやperiodic operatorの前計算を含む |
| state更新 | $O(N)$に近い |
| 1点評価 | $O(\log N+N_{\mathrm{near}})$に近い |
| 多点評価 | targetごとの評価をOpenMPで並列化 |

実際の定数係数とmemory使用量は、展開係数数、source/target node数、M2L pair数、near list、periodic画像数に依存します。
小ケースではplan構築とstate更新の固定費がDirectを上回ることがあります。速度比較はdebug buildではなく
release profileで行い、初期化時間、batch更新時間、粒子評価時間を分けて見ます。

## Direct比較から帯電結果まで収束させる

推奨する確認順は次のとおりです。

1. 同じmesh、kernel、softening、正規化の小ケースをDirectとFMMで実行する。
2. 表面近傍、遠方、強い電荷勾配、正負電荷が相殺する領域で電場と電位を比較する。
3. `tree_theta`を小さくし、`tree_leaf_max`も変えて観測量の安定性を確認する。
4. source meshを細かくし、FMM近似とは独立にkernel離散化の収束を確認する。
5. 粒子の命中要素、escape数、batch後の`q_elem`、電荷収支を比較する。
6. periodic2では有限画像層、far correction、zero mode、cache cold/warmの一致を個別に確認する。

相対誤差は基準場がほぼ0の点で不安定になるため、絶対誤差または代表場で正規化した誤差も併記します。
FMMの一点誤差が小さくても、粗いpanel meshや大きい粒子時間刻みの誤差は小さくなりません。

## 現行実装で選べる範囲

- kernelはCoulomb固定
- simulator adapterの展開次数は4固定
- source geometryはplan構築後に固定
- 対応する場境界はfreeとperiodic2
- triangle P0はsofteningなし、insulator-onlyのPhase 1
- tree外targetはDirect fallback
- periodic zero modeとouter responseはFMM core単独では完結しない

係数配列、interaction cache、translation前計算、parallel loopなどの内部仕様は
[Coulomb FMMコア内部実装](FMMCore.html)に分けています。periodic operatorの内部仕様は
[periodic2遠方補正](PeriodicFarCorrection.html)から辿れます。

## Code reference

- BEACH adapterと固定order: [`bem_field_solver_config.f90`](../src/physics/field_solver/bem_field_solver_config.f90)
- plan/stateの更新: [`bem_field_solver_tree.f90`](../src/physics/field_solver/bem_field_solver_tree.f90)
- FMM public API: [`bem_coulomb_fmm_core.f90`](../src/physics/field_solver/fmm/api/bem_coulomb_fmm_core.f90)
- planとinteraction構築: [`bem_coulomb_fmm_plan_ops.f90`](../src/physics/field_solver/fmm/internal/tree/bem_coulomb_fmm_plan_ops.f90)
- stateのupward/downward pass: [`bem_coulomb_fmm_state_ops.f90`](../src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_state_ops.f90)
- target評価とfallback: [`bem_coulomb_fmm_eval_ops.f90`](../src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90)
- 主なsolver回帰テスト: [`test_dynamics_fmm.f90`](../tests/fortran/test_dynamics_fmm.f90)、[`test_dynamics_panel_fmm.f90`](../tests/fortran/test_dynamics_panel_fmm.f90)

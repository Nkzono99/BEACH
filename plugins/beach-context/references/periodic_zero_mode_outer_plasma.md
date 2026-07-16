title: periodic2 zero modeと外部プラズマ

# periodic2 zero modeと外部プラズマ

この文書は、2軸周期・1軸開放の静電場を、横方向に変化する`k!=0`成分、
平面平均された`k=0`成分、外部プラズマへどう分けて計算するかを説明します。
FMM自体の展開式は[FMMCore](FMMCore.md)、入力の選び方は
[FieldSolvers](FieldSolvers.md)、物理モデルの適用範囲はADR 0001とADR 0002を参照してください。

## 1. 最初に押さえる全体像

周期軸をx/y、開放軸をzとします。表面電荷が作る場を横方向Fourier modeに分けると、
役割は次のようになります。

| 成分 | 何を表すか | production実装 |
| --- | --- | --- |
| `k!=0` | x/y方向の局所的な電荷むら | FMM有限画像和 + `cached_kneq0`遠方operator |
| surface `k=0` | 各高さより下にある総電荷が作る平面平均場 | triangle-height累積多項式 |
| outer profile | interfaceから無限遠までのplasma応答 | `kinetic_1d`または`unified_linear_response` |

`cached_kneq0`は全電場を返すsolverではありません。意図的にsurface `k=0`を除いた
非零モードだけを返し、`electrostatic_snapshot`が物理的なzero modeを一度だけ合成します。
この分離により、横方向の無限周期補正と、開放方向の境界条件・シースモデルを独立に選べます。

## 2. `cached_kneq0`が保存するもの

### 2.1 有限画像和だけでは足りない

通常のperiodic FMMは、primary cellと設定した近傍画像

$$
(i,j)\in[-N,N]^2
$$

を陽に足します。これは近距離場を正確に扱えますが、範囲外に無限に続く周期画像の
滑らかな遠方場が欠けます。毎回Ewald和を全source・全particleへ評価すれば補えますが、
particle hot pathには重すぎます。

### 2.2 Ewald2P referenceとは

Coulombの$1/r$を周期画像へそのまま足すと収束が遅く、非中性slabでは平均場の境界条件も
明示しないと一意になりません。Ewald法は数値パラメータ$\alpha$を使って

$$
\frac{1}{r}=\frac{\operatorname{erfc}(\alpha r)}{r}
+\frac{\operatorname{erf}(\alpha r)}{r}
$$

と分けます。第1項は実空間で急速に減衰し、第2項は滑らかなので逆空間で少数modeへ展開できます。
BEACHのEwald2Pではx/yだけを逆格子へ変換し、zは開放座標のまま保持します。

| 部分 | BEACHでの評価 |
| --- | --- |
| real space | screened Coulombを有限画像範囲で和を取る |
| reciprocal `k!=0` | x/y逆格子modeを有限層まで和を取る |
| `k=0` | 開放z方向の平均場として別項に保持する |

$\alpha$は実空間と逆空間へ仕事を振り分ける**数値分割パラメータ**であり、Debye長や物理的な
遮蔽率ではありません。十分に収束した結果は$\alpha$へ依存しません。実装で`exact Ewald`と呼ぶ値は、
設定したreal/reciprocal cutoffでコードが評価する高精度teacherであり、無限和の解析的厳密値ではありません。

`cached_kneq0`のcold buildは、このEwald2P teacherから有限画像shellを引いた

$$
R_\mathrm{Ewald}^{\mathrm{full}}
=K_\mathrm{Ewald2P}^{\mathrm{truncated}}-K_\mathrm{shell}(N)
$$

だけをproxy/check pointでsampleします。Ewald2Pはoperatorを作る教師として使われ、warm runtimeの
particle evaluationには入りません。real-space項、逆空間項、zero-mode項の実装式は
[FMMCore 8.2](FMMCore.md#82-periodic2-ewaldewald2p補正)にあります。

### 2.3 遠方補正を線形operatorにする

幾何、周期長、FMM orderとtarget topologyが固定なら、source multipoleから遠方local展開への
写像は線形です。cold buildはproxy sourceとcheck pointを使い、Ewald referenceと有限画像和の差を
再現する行列を一度だけfitします。cacheに入るのは電荷stateや粒子位置ではなく、この行列です。

$$
\mathbf L_t^{\mathrm{far}}=\mathbf A_t\mathbf M_{\mathrm{root}}
$$

ここで$\mathbf M_{\mathrm{root}}$は現在電荷から作るroot multipole、
$\mathbf A_t$はtarget anchorごとのcached operator、$\mathbf L_t^{\mathrm{far}}$はlocal展開です。
電荷が変わっても$\mathbf A_t$は再利用でき、`update_state`では行列を現在のmultipoleへ適用するだけです。

### 2.4 なぜ`k=0`を引くのか

Ewald teacherのfull-periodic residualには、便宜的な対称vacuum `k=0`も含まれます。
しかし物理的なzero modeは`lower_boundary_model`やouter-plasma closureで決める必要があります。
そこでFMM core内部では、同じsource stateから対称`k=0`を再構成してcached結果から引きます。

$$
K_{k\ne0}=K_{\mathrm{shell}}+R_{\mathrm{Ewald}}^{\mathrm{full}}-K_0^{\mathrm{sym}}
$$

その後、場の合成処理が選択された境界条件の$K_0^{\mathrm{physical}}$を加えます。

$$
K_{\mathrm{surface}}=K_{k\ne0}+K_0^{\mathrm{physical}}
$$

`exclude_k0`は「平均場を無視する」という意味ではなく、**非零モードbackend内で重複加算しない**
というownership規則です。

### 2.5 cold buildとwarm run

| 時期 | 実行する処理 | 実行しない処理 |
| --- | --- | --- |
| cache miss | proxy/check配置、Ewald teacher評価、QR fit、checksum付きatomic publish | particle追跡 |
| cache hit | fingerprint・shape・checksum検証、operator読込 | Ewald再評価、再fit |
| charge refresh | P2M/M2M、cached行列適用、zero-mode state更新 | operator再生成 |
| particle eval | local展開、near direct、cached対称`k=0`減算、物理`k=0`加算 | all-source Ewald和 |

cache identityにはgeometry、periodic length、order、画像/Ewald層、source kernel、target topology、
generator/build versionなどが含まれます。一致しないcacheを近似的に流用しません。

## 3. surface `k=0`の数値アルゴリズム

### 3.1 三角形を高さ方向の累積電荷へ投影する

三角形$i$の総電荷を$q_i$とし、一様P0面密度を仮定します。
$F_i(z)$を「三角形面積のうち高さがz以下にある割合」とすると、平面平均された累積電荷は

$$
C(z)=\sum_i q_iF_i(z)
$$

です。$F_i$は3頂点の高さの間で区分二次関数になります。plan構築時に全頂点高さをsortして
breakpointを作り、各三角形が各区間へ加える二次係数を保存します。水平三角形はsheet chargeとして
別に保持し、面上評価ではminus trace、plus trace、principal valueを区別します。

### 3.2 charge refresh

batchで$q_i$が変わると、区間差分へ$q_i$倍した係数を加え、prefix sumで
$C(z)=a_0+a_1z+a_2z^2$の区間係数を作ります。さらに区間積分をbreakpointごとに前計算します。
geometry planは再構築しません。

### 3.3 fieldとpotential評価

周期セル面積を$A=L_xL_y$、下側遠方場を$E_{\mathrm{bottom}}$とすると、Gauss則から

$$
E_0(z)=E_{\mathrm{bottom}}+\frac{C(z)}{\epsilon_0A}
$$

です。電位はgauge点$(z_g,\phi_g)$から同じ累積電荷を積分して

$$
\phi_0(z)=\phi_g-E_{\mathrm{bottom}}(z-z_g)
-\frac{1}{\epsilon_0A}\int_{z_g}^{z}C(\zeta)\,d\zeta
$$

と評価します。区間探索は二分探索、区間内は二次式と三次primitiveなので1点あたり$O(\log N_z)$です。

### 3.4 lower boundary closure

総表面電荷を$Q=\sum_iq_i$とすると、現在の選択肢は次の2つです。

| model | $E_{\mathrm{bottom}}$ | $E_{\mathrm{top}}$ | 意味 |
| --- | ---: | ---: | --- |
| `symmetric_vacuum` | $-Q/(2\epsilon_0A)$ | $+Q/(2\epsilon_0A)$ | 上下に同じ真空半空間がある無外場closure |
| `e_bottom_zero` | $0$ | $Q/(\epsilon_0A)$ | 下側電束を0へ固定する旧計算再現用closure |

これは誘電体内部の遮蔽を解くモデルではありません。`symmetric_vacuum`は追加のinterface位置や
誘電率を導入しない最小の対称closure、`e_bottom_zero`は物理既定ではなくlegacy optionです。

## 4. outer-plasma modelとの接続

### 4.1 split modelとunified modelは違う

自己整合な外部シースの標準はsplit `kinetic_1d`です。`unified_linear_response`はその高精度版ではなく、
split windowを置けず、かつ線形応答で十分なrough surfaceに限る高度なscreening modelです。

split modelでは、局所領域と外部領域はinterfaceで接続されます。

| 領域 | 使用する場 |
| --- | --- |
| meshからinterfaceまで | `k!=0` surface field + surface `k=0` |
| interfaceより外側 | 1D outer profile。split contractでは横方向modeが十分減衰していることを要求 |

したがってouter profileを局所領域の全点へ加算しません。interfaceで電位と法線電場を一致させ、
外部へ出た粒子だけが同じprofileをreturn/escape判定に使います。

`unified_linear_response`は別経路です。surface projectionからfar boundaryまで一つのzero-mode
Poisson gridを解き、nonzero modeもplasma tailへ連続接続します。particle interfaceはfield境界ではなく
ownership面だけです。

流入粒子の加減速、外向き粒子のturning/return、Zhao系注入補正との違いは
[外部シースとreservoir粒子境界](sheath_reservoir_boundary.md)で説明します。

### 4.2 model別の数値処理

| `outer_plasma.model` | 位置付け | zero-mode処理 | 主な適用範囲 |
| --- | --- | --- | --- |
| `none` | 外部シースなし | surface `k=0`だけをboundary closureで評価 | 外部plasmaを置かない |
| `linear_debye` | 簡易・reference | $\Delta\phi\exp[-(z-z_I)/\lambda_D]$ | 小振幅split reference |
| `kinetic_1d` | **標準・推奨** | VDF closureを含む非線形1D Poisson solve | 単調・無衝突・非磁化sheath |
| `unified_linear_response` | 高度・限定用途 | rough surfaceとplasma sourceを同じ1D Poisson gridへ投入 | 線形応答でsplit windowがない場合 |

`cached_kneq0`のproduction経路が受理するouter modelは現在`none`、`kinetic_1d`、
`unified_linear_response`です。`linear_debye`は`panel_spectral_reference`の小規模検証経路です。

### 4.3 `linear_debye`

interfaceを$z_I$、$\Delta\phi=\phi_I-\phi_\infty$とすると

$$
\phi(z)=\phi_\infty+\Delta\phi e^{-(z-z_I)/\lambda_D},\qquad
E(z)=\frac{\Delta\phi}{\lambda_D}e^{-(z-z_I)/\lambda_D}
$$

を使います。$|\Delta\phi|/V_T$が設定した線形性上限を超える場合は適用外です。

### 4.4 `kinetic_1d`

interface fieldをsurface zero modeから受け取り、ambient electron half-Maxwellian、cold drifting ion、
任意のphotoelectron fluxから$\rho(\phi)$を構成し、

$$
\frac{d^2\phi}{dz^2}=-\frac{\rho(\phi)}{\epsilon_0}
$$

を解きます。格子は伸長可能な1D grid、内部点はconservative finite-volume residual、遠方は

$$
\phi'(L)+\frac{\phi(L)}{\lambda_{\mathrm{tail}}}=0
$$

です。解析density derivativeからbordered-tridiagonal Jacobianを組み、1 Newton stepは$O(N_z)$です。
backtracking、pseudo-transient、interface-field continuationは収束経路だけを改善し、最終受理には
元のPoisson residual、単調分枝、ion accessibility、Bohm条件をすべて要求します。MPIではrootが解き、
statusとprofileをbroadcastします。失敗時に別sheathや前batch解へfallbackしません。

### 4.5 `unified_linear_response`

高さごとのplasma accessible fractionを$f_{\mathrm{access}}(z)$、
$\kappa=1/\lambda_D$とすると、surfaceの平面平均sourceと

$$
\rho_{\mathrm{plasma}}(z)=-\epsilon_0f_{\mathrm{access}}(z)\kappa^2\phi(z)
$$

を同じ非一様1D gridへ入れます。bottom fieldとfar Robin条件を含むtridiagonal Poisson系を解きます。
非零モードはresponse startより上で$\alpha=\sqrt{k^2+\kappa^2}$のtailへ接続し、電位、法線場、
接線場を連続にします。linearity、height-field geometry、accessible-area収束、mode truncationの
いずれかがcontract外なら停止します。

## 5. batchごとの処理順

| 順序 | 処理 |
| ---: | --- |
| 1 | 新しい`q_elem`からFMM multipoleとsurface zero-mode stateを更新 |
| 2 | 必要なbatchでouter profileを解く。前profileは同一identity時だけ初期値に利用 |
| 3 | outer interface potentialに合うようzero-mode potential gaugeを設定 |
| 4 | Gauss residualとinterface診断を更新 |
| 5 | 粒子評価ではnonzero、zero、prescribed fieldをownership規則どおり1回ずつ合成 |

field evaluationは同じbatch中に固定された場を使います。個々の粒子命中のたびにoperatorや
outer profileを更新せず、batch commit後にまとめて更新します。

## 6. 診断と停止条件

最低限確認する量は、`interface_potential`、`interface_field`、`gauss_residual`、
`outer_integrated_charge`、nonlinear residual、linearity ratio、cache fingerprint/hitです。

Gauss診断では上側へ出るsurface電束を

$$
Q_{\mathrm{upper\ flux}}=Q+\epsilon_0AE_{\mathrm{bottom}}
$$

とし、outer integrated chargeとの和を残差にします。数値収束だけでなく、Bohm条件、単調性、
線形性など物理的な適用条件を満たさない場合もfail closedです。

## 7. 実装対応

| 処理 | Fortran実装 |
| --- | --- |
| cached nonzero operator | `fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90` |
| cached symmetric `k=0` subtraction | `fmm/internal/runtime/bem_coulomb_fmm_eval_ops.f90` |
| surface zero-mode plan/state | `periodic_zero_mode/bem_periodic_zero_mode_plan.f90` |
| zero-mode evaluation | `periodic_zero_mode/bem_periodic_zero_mode_eval.f90` |
| field ownership/composition | `bem_electrostatic_snapshot.f90` |
| linear outer model | `outer_plasma/bem_outer_plasma_linear.f90` |
| kinetic outer model | `outer_plasma/bem_outer_plasma_kinetic.f90` |
| unified linear model | `outer_plasma/bem_outer_plasma_unified.f90` |

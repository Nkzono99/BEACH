title: periodic2遠方補正

# periodic2遠方補正

`periodic2` の FMM は primary cell と有限画像を通常の tree FMM で評価し、その外側の滑らかな無限周期遠方場を
追加 operator で補います。このページは Ewald2P teacher、`cached_kneq0`、cache、FMM state への接続の正本です。

## 有限画像モデルと無限周期近似を選ぶ

```toml
[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0e-4, 1.0e-4, 3.2e-4]
periodic_axes = ["x", "y"]

[fields]
boundary = "periodic2"

[fields.solver]
method = "fmm"

[fields.periodic]
backend = "cached_kneq0"
image_layers = 1
cache_dir = ".beach_cache/periodic2"
generation_tolerance = 1.0e-8
```

`[domain]` が固定 target topology と周期軸を所有し、`fields.boundary` が periodic2 場を選びます。
`cached_kneq0` は x/y 周期・z 非周期の domain box 内 target だけを受理します。

| `fields.periodic.backend` | 計算内容 | 用途 |
| --- | --- | --- |
| `finite_images` | primary + 有限画像だけ | 有限画像model、収束比較 |
| `cached_kneq0` | versioned operatorを生成・再利用 | 無限周期nonzero modeのproduction経路 |
| `panel_spectral_reference` | 非零modeを直接のFourier和で評価（direct solver） | 小規模な検証の参照 |

`cached_kneq0` は z 方向の平均場 `k=0` を決めません。[periodic2 静電場](../models/PeriodicElectrostatics.md)の
場の合成処理が physical `k=0` を一度だけ加えます。このoperator自体は外部plasma/sheathを解かず、外部シース closure の
電位基準はzero-mode構築とは別の層で与えます。
`m2l_root_oracle`は削除済みで、設定すると`cached_kneq0`を案内して起動時にrejectされます。

## 有限画像shellを通常のFMMで評価する

周期軸を1、2とし、周期長を$L_1,L_2$とします。画像層$N$では

$$
\mathcal I_N=\{(i,j)\mid -N\le i,j\le N\}
$$

のshift

$$
\mathbf s_{ij}=iL_1\mathbf e_1+jL_2\mathbf e_2
$$

をsourceへ加えます。near listは固定 P0 triangle kernelでDirect評価し、well-separatedなnode pairは
各画像shiftを含めたM2L derivativeで評価します。この部分を$K_\mathrm{shell}(N)$と書きます。

`fields.periodic.backend="finite_images"`では$K_\mathrm{shell}(N)$が全結果です。画像層の外側は暗黙に近似されるのではなく、
寄与しません。

## Ewald2Pをruntime kernelではなくteacherとして使う

2軸周期・1軸開放のCoulomb和を作るため、kernelを数値parameter$\alpha$で分けます。

$$
\frac{1}{r}=
\frac{\operatorname{erfc}(\alpha r)}{r}
+\frac{\operatorname{erf}(\alpha r)}{r}
$$

第1項をreal-space画像和、第2項を周期2軸のreciprocal modeとして評価し、開放軸は実座標のまま残します。

| 成分 | 数値的な役割 |
| --- | --- |
| real space | 短距離側をscreened Coulomb和として収束させる |
| reciprocal `k\ne0` | 横方向に変化する滑らかなmodeを足す |
| `k=0` | 平面平均場。物理的なz境界条件と分離して扱う |

$\alpha$はreal/reciprocal計算量の配分parameterで、Debye遮蔽ではありません。実装上のteacherは設定した
`fields.periodic.ewald_layers`で打ち切った高精度有限和です。理論上の無限和を毎回直接評価するわけではありません。

teacherから有限画像分を引いた滑らかな残差を

$$
R_\mathrm{Ewald}^{\mathrm{full}}
=K_\mathrm{Ewald2P}^{\mathrm{truncated}}-K_\mathrm{shell}(N)
$$

とします。Ewald評価はoperatorを作るcold pathだけで使い、通常の粒子評価では実行しません。

## source節点のmultipoleからtarget localへのoperatorをfitする

固定geometryでは、遠方残差はsource電荷に対して線形です。source treeの深さ2の節点$s$のmultipole係数を
$\mathbf M_s$、target treeの深さ2のanchor node$t$のfar local係数を$\mathbf L_t^\mathrm{far}$とすると、

$$
\mathbf L_t^\mathrm{far}=\sum_s\mathbf A_{t,s}\mathbf M_s
$$

と書けます。線源全体を根のmultipole 1個で見ると、粒子層の高さが除いた画像までの距離（約$2L$）に近いため
展開がほとんど収束せず、帯電した三層レゴリスでは非零モードに数Vの面平均が漏れていました。深さ2の節点に
分けると節点が小さくなり、評価側のanchorも小さくなるので、同じ展開次数で誤差が大きく下がります。
cold buildは次の順に$\mathbf A_{t,s}$を作ります。

1. source節点ごとに節点を囲むproxy点を置き、proxy電荷からその節点のmultipoleへの行列と最小norm擬似逆を作る。
2. target anchorごとにcheck点を置き、check点の電場をlocal展開で再現する行列を正則化QRで分解する（source節点で共有）。
3. source節点ごとに、各proxy sourceのEwald residual（下の面平均の補正を含む）をcheck点で評価し、local係数を解く。
4. proxy-to-local写像と擬似逆を合成する。
5. geometry、source節点、target anchor、設定からfingerprintを作り、operatorとchecksumをcacheへ公開する。

組の数はtarget anchor数×source節点数（三層レゴリスで最大64×64）ですが、実行時の追加は小さな行列・vector積だけです。
初回の作成は長くなり（三層レゴリスで28スレッド1プロセス約27分、MPIではtarget anchorをrankに分ける）、
cacheに保存して再利用します。

`cached_kneq0`はこの線形写像をversioned cacheへ保存し、fingerprint、shape、checksumが一致するwarm runで
再利用します。proxy pointから作ったoperatorを、P0 triangle sourceの
triangle-averaged P2M係数へ適用します。

電場だけのfitではlocal展開の定数potential係数を決められないため、同じcheck点のpotential residualから
定数係数を別にfitします。これは電場のfitとpotential gaugeを同じ単位のleast-squares列へ混ぜないためです。

proxy 行列の各 multipole 行は長さの次数が異なるため、擬似逆の QR 分解前に行を正規化し、
右辺にも同じ尺度を適用します。これにより µm 単位の形状でも高次 moment が絶対しきい値で失われません。
cache の fingerprint は生成器の版を含むので、生成方法を変えた版は古い cache を再利用しません。

物理的な平均場は電荷源内部の高さ分布を保ち、電荷源より上では累積電荷を厳密に総電荷 $Q$ とします。
これにより、傾いた三角形の累積多項式を加減算した丸め残差が真空域へ外挿されることを防ぎます。
これは高さ別の平均場を消す処理ではありません。

## 通常M2Lの後、L2Lの前に補正を加える

batchごとの`update_state`は次の順で進みます。

```text
P2M
  -> M2M
  -> ordinary tree M2L
  -> L_far(t) += sum_s A_ts M_s  periodic far correction
  -> L2L
  -> L2P + near Direct
```

実装上は各target anchor nodeについて

$$
\mathbf L_t\leftarrow\mathbf L_t+
\sum_s\mathbf A_{t,s}\mathbf M_s
$$

を行います。その後は通常のL2Lが補正済みlocalをleafへ伝え、L2Pが粒子位置で評価します。したがってhot pathで
増える処理は行列・vector積と通常のlocal展開評価であり、all-source Ewald和ではありません。

tree の M2L pair cache と `m2l_deriv` は変更せず、補正を ordinary M2L と L2L の間に追加します。

## cached operatorから対称`k=0`を除く

Ewald teacher の full residual には fit を定義する対称 vacuum `k=0` が含まれます。有限画像殻
$K_\mathrm{shell}(N)$ をセル上で面平均すると、無限面電荷の平均場 $K_0^\mathrm{sym}$（高さ $h$ の単位電荷あたり
$-2\pi|h|/A$）と滑らかな残り $s(h)$ の和になります。$s$ は $1/r$ を画像ブロック $[c-(N+\tfrac12)L,\,c+(N+\tfrac12)L]^2$
上で積分した閉じた式から $K_0^\mathrm{sym}$ を引いたもので、$h=0$ でも解析的です（$c$ は線源の木の根の中心）。
残差の面平均はちょうど $-s$ なので、fit には面平均 0 の $R+s$ だけを渡し、$s$ は線源ごとの高さの表から
$S(z)=\sum_i q_i s_i(z)$ として実行時に引きます。FMM 側の nonzero kernel は

$$
K_{k\ne0}=K_\mathrm{shell}(N)
+\mathcal F\left[R_\mathrm{Ewald}^{\mathrm{full}}+s\right]-K_0^\mathrm{sym}-S
$$

で（$\mathcal F$ は上の operator による fit）、場の合成時に物理的な平均場を一度だけ加えます。
$S$ の高さの表は間隔 $L/32$ の 3 次 Hermite 補間で、電場は同じ多項式の微分なので電位の勾配と一致します。

$$
K_\mathrm{surface}=K_{k\ne0}+K_0^\mathrm{physical}
$$

`fields.periodic.backend="cached_kneq0"` は非零モードだけを評価し、physical zero mode を別に加えることで平均場の二重加算を
防ぎます。境界モデルと Gauss 則は [periodic2 静電場](../models/PeriodicElectrostatics.md)を参照してください。

## cacheと並列化の境界を理解する

| 段階 | 主な処理 |
| --- | --- |
| cold miss | Ewald teacher、QR fit、operator生成、checksum付きatomic publish |
| warm hit | fingerprint・shape・checksum検証、operator読込 |
| charge refresh | P2M/M2M/M2L、保存済みoperator適用、zero-mode state更新 |
| particle evaluation | local展開 + near Direct。Ewald再評価なし |

cache I/OとlockはMPI rootが担当します。cold buildではtarget operator sliceをMPI rankへ分配し、target内のproxy列を
OpenMPで評価し、完成したoperatorをrank間で集約します。charge履歴や粒子位置はcacheに含めません。

`cached_kneq0`はconfigured target box内の評価を前提とし、box外targetをDirect fallbackしません。これはcached
operatorが固定target topologyに対するproduction kernelだからです。

上端 H 近傍の局所展開誤差は、[上部真空の Fourier 切替](../models/PeriodicElectrostatics.md#上部真空の-fourier-評価への切替)で回避できます。
`fields.periodic.upper_fourier.switch_height="auto"` は上部真空の非零成分を Fourier 評価へ切り替えます。

## 収束とcache診断を確認する

1. `fields.periodic.image_layers`を増やし、near shellとoperatorの分割に依存しないことを確認する。
2. `fields.periodic.ewald_layers`と必要なら`fields.periodic.ewald_alpha`を変え、teacherの打切りを確認する。
3. cold buildとwarm cacheで電場・電位・主要観測量が丸め誤差範囲で一致することを確認する。
4. `periodic2_cache_fingerprint`、`periodic2_cache_hit`、`periodic2_operator_build_count`を確認する。
5. nonzero modeだけでなく、選択したphysical `k=0`とGauss residualも別に確認する。
6. 非零modeの面平均が0であることを確認する。外部シースと接続する場合は、`top_reference_history.csv`の上端面の
   平均電位が外部根の $\phi_0$ に一致する。`make test-fortran-far-correction`は一様面電荷で同じ性質を検査する。

### cold build だけを先に済ませる

operator の作成は target anchor を MPI rank に分配するので、本番 run より多い rank で作ると短くなる
（三層レゴリスで 28 コア合計の 1 プロセス × 28 スレッドと 4 rank × 7 スレッドはどちらも約 30 分。anchor は最大 64 個）。
`beach --prepare-periodic-cache beach.toml` は設定とメッシュを読み、operator を作成（または cache から読込）して
publish したら終了する。結果ファイルや checkpoint は書かず、`[run.restart]` も無視する。

```bash
# 例: 16 rank × 7 スレッドで cache だけ作り、本番は 4 rank で走らせる
srun -n 16 beach --prepare-periodic-cache beach.toml
srun -n 4 beach beach.toml
```

標準出力に `periodic2_cache_hit`、`periodic2_operator_build_count`、`periodic2_cache_fingerprint`、`periodic2_cache_path` を出す。
同じ fingerprint になる条件（box、tree の節点、展開次数、Ewald 設定、kernel）は上の fingerprint の節と同じ。

長さの単位に依存しないことやcacheの再利用だけでは、弱い外部シース場と比べられる精度は保証されません。帯電した
三層レゴリスの電荷（画像層1、展開次数4）で、層の上の真空域を直接のFourier和と比べると、非零modeの面平均は0.02 V以下、
横方向の電位誤差は0.01〜0.03 V、電場誤差は約2×10³ V/mでした。

有限画像和の意味と、無限周期・面平均成分の設定は[periodic2 の静電場](../models/PeriodicElectrostatics.md)にあります。

## Code reference

- mode選択とperiodic helper: [`bem_coulomb_fmm_periodic.f90`](../../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic.f90)
- Ewald2P teacher: [`bem_coulomb_fmm_periodic_ewald.f90`](../../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_ewald.f90)
- operator生成（source節点とtarget anchorの組）: [`bem_coulomb_fmm_periodic_root_ops.f90`](../../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90)
- 有限画像殻の面平均の滑らかな残り $s$ と高さの表: [`bem_periodic_shell_mean.f90`](../../src/physics/field_solver/periodic/bem_periodic_shell_mean.f90)
- cache I/O: [`bem_coulomb_fmm_periodic_cache.f90`](../../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_cache.f90)
- stateへのoperator適用: [`bem_coulomb_fmm_state_ops.f90`](../../src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_state_ops.f90)

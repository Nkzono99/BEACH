title: periodic2遠方補正

Lang: [日本語](PeriodicFarCorrection.md) | [English](PeriodicFarCorrection.en.md)

# periodic2遠方補正

`periodic2` の FMM は primary cell と有限画像を通常の tree FMM で評価し、その外側の滑らかな無限周期遠方場を
追加 operator で補います。このページは Ewald2P teacher、`cached_kneq0`、cache、FMM state への接続の正本です。

## 有限画像モデルと無限周期近似を選ぶ

```toml
[sim]
field_solver = "fmm"
field_periodic_image_layers = 1
field_periodic_far_correction = "cached_kneq0"
field_periodic_cache_dir = ".beach_cache/periodic2"
field_periodic_generation_tolerance = 1.0e-8

[domain]
box_origin = [0.0, 0.0, 0.0]
box_size = [1.0, 1.0, 1.0]
periodic_axes = ["x", "y"]

[field_boundary]
mode = "periodic2"
```

`[domain]` が固定 target topology と周期軸を所有し、`[field_boundary]` が periodic2 場を選びます。
現行 `cached_kneq0` は x/y 周期・z 非周期の domain box 内 target だけを受理します。

| `field_periodic_far_correction` | 計算内容 | 用途 |
| --- | --- | --- |
| `none` | primary + 有限画像だけ | 有限画像model、収束比較 |
| `auto` | 現行では`none`へ正規化 | 互換性のみ |
| `cached_kneq0` | versioned operatorを生成・再利用 | 無限周期nonzero modeのproduction経路 |

`cached_kneq0` は z 方向の平均場 `k=0` を決めません。[periodic2 静電場](PeriodicElectrostatics.html)の
場の合成処理が physical `k=0` を一度だけ加えます。self-consistent outer-plasma/sheath model は非対応です。
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

`field_periodic_far_correction="none"`では$K_\mathrm{shell}(N)$が全結果です。画像層の外側は暗黙に近似されるのではなく、
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
`field_periodic_ewald_layers`で打ち切った高精度有限和です。理論上の無限和を毎回直接評価するわけではありません。

teacherから有限画像分を引いた滑らかな残差を

$$
R_\mathrm{Ewald}^{\mathrm{full}}
=K_\mathrm{Ewald2P}^{\mathrm{truncated}}-K_\mathrm{shell}(N)
$$

とします。Ewald評価はoperatorを作るcold pathだけで使い、通常の粒子評価では実行しません。

## root multipoleからtarget localへのoperatorをfitする

固定geometryでは、遠方残差はsource電荷に対して線形です。source rootのmultipole係数を
$\mathbf M_\mathrm{root}$、target anchor node$t$のfar local係数を$\mathbf L_t^\mathrm{far}$とすると、

$$
\mathbf L_t^\mathrm{far}=\mathbf A_t\mathbf M_\mathrm{root}
$$

と書けます。cold buildは次の順に$\mathbf A_t$を作ります。

1. source rootを囲むproxy点を置き、proxy電荷からroot multipoleへの行列を作る。
2. target anchorごとにcheck点を置き、各proxy sourceのEwald residualを評価する。
3. check点の電場をlocal展開で再現する係数を正則化QRで解く。
4. proxy-to-local写像とproxy-to-multipoleの最小norm擬似逆を合成する。
5. geometryと設定からfingerprintを作り、operatorとchecksumをcacheへ公開する。

`cached_kneq0`はこの線形写像をversioned cacheへ保存し、fingerprint、shape、checksumが一致するwarm runで
再利用します。proxy pointから作ったoperatorを、P0 triangle sourceの
triangle-averaged P2M係数へ適用します。

電場だけのfitではlocal展開の定数potential係数を決められないため、同じcheck点のpotential residualから
定数係数を別にfitします。これは電場のfitとpotential gaugeを同じ単位のleast-squares列へ混ぜないためです。

## 通常M2Lの後、L2Lの前に補正を加える

batchごとの`update_state`は次の順で進みます。

```text
P2M
  -> M2M
  -> ordinary tree M2L
  -> L_far(t) += A_t M_root      periodic far correction
  -> L2L
  -> L2P + near Direct
```

実装上は各target anchor nodeについて

$$
\mathbf L_t\leftarrow\mathbf L_t+
\mathbf A_t\mathbf M_\mathrm{root}
$$

を行います。その後は通常のL2Lが補正済みlocalをleafへ伝え、L2Pが粒子位置で評価します。したがってhot pathで
増える処理は行列・vector積と通常のlocal展開評価であり、all-source Ewald和ではありません。

tree の M2L pair cache と `m2l_deriv` は変更せず、補正を ordinary M2L と L2L の間に追加します。

## cached operatorから対称`k=0`を除く

Ewald teacher の full residual には fit を定義する対称 vacuum `k=0` が含まれます。実際の平均場は
`lower_boundary_model` で決まるため、FMM 側の nonzero kernel を

$$
K_{k\ne0}=K_\mathrm{shell}(N)
+R_\mathrm{Ewald}^{\mathrm{full}}-K_0^\mathrm{sym}
$$

とし、場の合成時に物理的な平均場を一度だけ加えます。

$$
K_\mathrm{surface}=K_{k\ne0}+K_0^\mathrm{physical}
$$

`zero_mode_policy="exclude_k0"` は平均場を捨てる指定ではなく、FMM backend と physical zero mode の二重加算を
防ぐ ownership 規則です。境界モデルと Gauss 則は [periodic2 静電場](PeriodicElectrostatics.html)を参照してください。

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

## 収束とcache診断を確認する

1. `field_periodic_image_layers`を増やし、near shellとoperatorの分割に依存しないことを確認する。
2. `field_periodic_ewald_layers`と必要なら`field_periodic_ewald_alpha`を変え、teacherの打切りを確認する。
3. cold buildとwarm cacheで電場・電位・主要観測量が丸め誤差範囲で一致することを確認する。
4. `periodic2_cache_fingerprint`、`periodic2_cache_hit`、`periodic2_operator_build_count`を確認する。
5. nonzero modeだけでなく、選択したphysical `k=0`とGauss residualも別に確認する。

有限画像和の使い方は[有限画像構成](FinitePeriodicConfiguration.html)、無限周期とzero modeの
設定は[periodic2静電場](PeriodicElectrostatics.html)にあります。

## Code reference

- mode選択とperiodic helper: [`bem_coulomb_fmm_periodic.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic.f90)
- Ewald2P teacher: [`bem_coulomb_fmm_periodic_ewald.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_ewald.f90)
- root operator生成: [`bem_coulomb_fmm_periodic_root_ops.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_root_ops.f90)
- cache I/O: [`bem_coulomb_fmm_periodic_cache.f90`](../src/physics/field_solver/fmm/internal/periodic/bem_coulomb_fmm_periodic_cache.f90)
- stateへのoperator適用: [`bem_coulomb_fmm_state_ops.f90`](../src/physics/field_solver/fmm/internal/runtime/bem_coulomb_fmm_state_ops.f90)

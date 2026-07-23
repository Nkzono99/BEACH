title: 外部シースとreservoir粒子境界

# 外部シースとreservoir粒子境界

この文書は、`reservoir_face`からの流入、電位差による加減速、流入不能粒子、外向き粒子の
反射・return・escape、およびZhao系シースを説明します。公開入力は
`external_boundary.field`、`external_boundary.particles`、`external_boundary.ordinary_open`へ分けます。
以下に現れる`outer_plasma`、`coupling`、`sim`のselectorは、明記した旧構文例を除き、
facadeから導出される実行時contractです。zero modeとouter Poisson solverは
[periodic2 zero modeと外部プラズマ](periodic_zero_mode_outer_plasma.md)を参照してください。

## 1. 「シース」と呼ばれる機能の区別

BEACHには目的の異なる複数のreduced modelがあります。

| 機能 | 公開入力 | 変更するもの | 空間的な$E(z)$をparticle pusherへ与えるか |
| --- | --- | --- | --- |
| scalar流入障壁 | `particles.inflow_model="infinity_barrier"` | 流入flux、法線速度cutoff、face速度 | いいえ |
| legacy Zhao流入 | `inflow_model="legacy_sheath"` + `legacy_sheath_model="zhao_*"` | reservoir密度/VDF cutoff、ion drift、photoemission current | いいえ |
| 光電子なし浮遊流入 | `legacy_sheath_model="floating_no_photo"` | electron cutoff | いいえ |
| kinetic outer sheath | `field.model="kinetic_1d"` | **標準:** 自己整合1D電位・電場・密度profile | outer領域で使用 |
| unified linear response | `field.model="unified_linear_response"` | **高度:** localからfarまでのrough-surface線形screening | はい |
| 通常open面の障壁 | `ordinary_open.model="potential_barrier"` | 外向き粒子のescape/反射 | いいえ |

`legacy_sheath_model="zhao_*"`は既存speciesの注入分布を事前補正する静的モデルで、
simulation中のsurface chargeから外部Poisson問題を毎batch解きません。一方、
`field.model="kinetic_1d"` + `kinetic_closure="zhao_charge_driven"`はsurface chargeと
光電子populationから外部profileを更新する別の自己整合経路です。この2種類を混同したり重ねたりしません。

## 2. 共通のエネルギー規約

無限遠の電位を$\phi_\infty$、reservoir faceまたはinterfaceの電位を$\phi_f$、粒子電荷を$q$、
質量を$m$とします。1D静電写像では接線速度を変えず、法線エネルギーを保存します。

$$
\frac12m v_{n,f}^2+q\phi_f
=\frac12m v_{n,\infty}^2+q\phi_\infty
$$

したがって

$$
v_{n,f}^2=v_{n,\infty}^2-B,\qquad
B=\frac{2q(\phi_f-\phi_\infty)}{m}
$$

です。

| $B$ | 結果 |
| ---: | --- |
| $B>0$ | faceへ向かう途中で減速。$v_{n,\infty}<\sqrt B$の粒子はfaceへ到達しない |
| $B=0$ | 法線速度を変更しない |
| $B<0$ | faceへ向かう途中で加速し、$v_{n,f}>v_{n,\infty}$ |

到達不能な無限遠粒子はsimulation particleとして生成されません。これは生成済み粒子をbox面で
鏡面反射する処理とは異なり、上流VDFのaccessibility cutoffです。

## 3. `reservoir_face`の基本サンプリング

`reservoir_face`は注入面の内向き法線を$\mathbf n$として、drifting Maxwellianから
flux-weighted法線成分をsampleします。上流法線速度の下限を$v_{\min}$とすると

$$
\Gamma_{\mathrm{in}}=
\int_{v_n\ge v_{\min}}v_n f(\mathbf v)\,d^3v
$$

からbatchのmacro粒子数を決めます。電位障壁がある場合は

$$
v_{\min}=\sqrt{\max(B,0)}
$$

を使い、到達可能な上流速度だけをsampleした後、前節の式でface速度へ写像します。
接線成分は保持されます。位置はface矩形で一様sampleし、短いvelocity jitterの後に周期軸をwrap、
非周期軸をbox内へclampします。

この順序は、単にface Maxwellianをcutoffするだけではありません。`infinity_barrier`とouter-profile
inflowでは、上流fluxの選別と法線速度の加減速を両方行います。

## 4. face平均電位による `infinity_barrier`

### 4.1 face電位

各batch冒頭で更新した電場・電位を使い、注入口矩形を
`injection_face_phi_grid_n x injection_face_phi_grid_n`のcell-centered点で評価します。

$$
\bar\phi_f=\frac{1}{N^2}\sum_{a,b}\phi(\mathbf x_{ab})
$$

point/`triangle_p0` kernel、periodic2、physical `k=0`、outer state、`e0`はその場と同じ規約です。
この**面平均scalar電位**から$B=2q(\bar\phi_f-\phi_\infty)/m$を作ります。
同じ格子評価から母標準偏差・最小・最大も集計し、Maxwellian reservoirで
$|q|\sigma_\phi>0.1(k_BT+m u_n^2/2)$の場合はMPI rootが初回と最終batchに警告します。
この診断で電位評価回数は増えません。

### 4.2 何を表していないか

このmodelはfaceまでの電位差だけを使い、途中の$E(z)$、turning位置、flight time、空間電荷を解きません。
横方向にface電位が大きく変化していても平均値一つへ縮約します。したがって小さな注入口または
ほぼ一様なfaceのreduced modelであり、rough-surface sheathのproduction既定にはしません。

## 5. outer profileによる流入

`field.model="linear_debye"`または`"kinetic_1d"`と
`particles.mode="same_batch"`を選ぶと、z-high `reservoir_face` speciesを無限遠VDFとして解釈します。
この構成では`inflow_model="auto"`が同じ1D profileへ流入ownershipを渡し、内部の
`electrostatic_1d_instant_return` transferをBEACHが導出します。

- `kinetic_1d`は各batchで更新済みの`phi_interface-phi_infinity`から$B$を計算します。
- linear-Debye referenceはzero-mode surface chargeからinterface電位差を作ります。
- 到達可能な上流粒子をsampleし、同じエネルギー式でinterface速度へ写像します。

`kinetic_1d`のPoisson solveではsurface zero modeが与える`interface_field`をNeumann条件に使います。
ただし流入速度の写像に直接使うのは$E_I$ではなく、解いたprofileの電位差$\phi_I-\phi_\infty$です。
電場はprofileを決め、電位差が粒子のエネルギー変化を決めます。

## 6. 外向き粒子のescapeとreturn

### 6.1 `ordinary_open.model="potential_barrier"`

open face通過点の電位を$\phi_b$、外向き法線速度を$v_n$とすると、無限遠へ必要な障壁は

$$
U_b=q(\phi_\infty-\phi_b)
$$

です。$U_b>\tfrac12mv_n^2$なら法線速度だけを反転し、それ以外はescapeとします。
cornerで複数open faceを同時に横切る一般化はせずfail closedです。通過点電位は粒子運動と同じsnapshot規約で
局所`e0`電位も含みます。一様電場には有限な無限遠電位がないため、併用時の`phi_infty`は有効なreservoir基準として整合させます。

### 6.2 linear-Debye instant return

interfaceを外向きに通過した粒子について

$$
v_{n,\infty}^2=v_{n,I}^2+\frac{2q(\phi_I-\phi_\infty)}{m}
$$

を評価します。非負ならescape、負なら指数Debye profile内にturning pointがあります。
解析式で往復時間を求め、return時に法線速度を反転します。接線位置は$\mathbf v_t\tau_\mathrm{outer}$だけ
進め、x/yをperiodic wrapします。

### 6.3 kinetic-profile return

`field.model="kinetic_1d"` + `particles.mode="same_batch"`は同じ判定を、
実際に収束した離散$\phi(z)$で行います。正規化後のreturn IDは
`kinetic_1d_profile_return`です。
各区分で

$$
v_n^2(z)=v_{n,I}^2+\frac{2q[\phi_I-\phi(z)]}{m}
$$

を評価し、符号が変わる最初の区分でturning pointを補間します。grid上端まで到達した後は
far Robin exponential tailを解析積分します。往復時間が
`field_evolution_timescale`に対して長すぎる場合、frozen-field近似を受理せず停止します。

これはouter領域を小さい時間刻みで順次進めるtrajectory integratorではありません。処理は次の順です。

1. interface通過時の同時刻位置・速度を受け取る。
2. $v_{n,\infty}^2$が非負なら、その時点でinfinity escapeに分類する。
3. 負なら、保存エネルギーから離散profile各点の$v_n^2(z)$を走査する。
4. 最初に$v_n^2\le0$となる区間でturning pointを補間する。
5. 区分線形電位上の片道時間を
   $$
   \Delta t=\frac{2\Delta z}{\sqrt{v_{n,a}^2}+\sqrt{v_{n,b}^2}}
   $$
   で積算し、必要ならRobin tailの時間を解析的に加える。
6. 往復時間$\tau_\mathrm{outer}$後のinterface復帰状態を直接構成する。

復帰時は法線速度だけを反転し、接線速度は維持します。接線位置を
$\mathbf v_t\tau_\mathrm{outer}$だけ進めてx/yをperiodic wrapし、interface crossingを検出したlocal stepの
remaining `dt`だけを通常stepperで再積分します。outer flightをglobal simulation timeへ加算しません。

したがってこのmodelは「interfaceで速度とscalar障壁だけを比較する即時鏡面反射」より詳細ですが、
「外部軌道を時間刻みで明示追跡する」ものではありません。後者が必要な場合は
`field.model="unified_linear_response"` + `particles.mode="same_batch"`を使います。

### 6.4 unified 3D orbit

この組合せから導出される`electrostatic_3d_explicit_orbit`は、zero modeとscreened nonzero tailを
合成した3D電場中を固定刻みvelocity-Verletで追跡します。ownership面へ戻ればlocalへ返し、
far planeを外向きに横切ればescapeです。エネルギー誤差、flight time、frozen-field ratio、
step上限を検査します。

## 7. legacy Zhaoシース注入補正

### 7.1 解いているもの

legacy Zhaoはsolar-wind electron、ion、photoelectronの解析密度closureと電流条件から、
surface電位$\phi_0$、必要な場合の電位極小$\phi_m$、有効solar-wind electron密度
$n_{\mathrm{swe},\infty}$を非線形root solveで決めます。太陽高度$\alpha$は

$$
n_{\mathrm{phe},0}=n_{\mathrm{phe,ref}}\sin\alpha
$$

としてphotoelectron sourceへ入ります。

| branch | 実装上の電位・population構造 |
| --- | --- |
| A | $\phi_0>0$かつ$\phi_m<0$の非単調profile。下側にcaptured photoelectron、上側にreflected solar-wind electron |
| B | $\phi_0>0$の単調profile。photoelectron captureを含み、solar-wind electron reflectionは持たない |
| C | $\phi_0<0$の単調profile。solar-wind electron reflectionを含み、photoelectron cutoffは0 |

`zhao_auto`は低太陽高度`alpha < 20 deg`でC, A, B、それ以外でA, B, Cの順に解を試します。
指定branchのrootが得られなければ別の物理modelへfallbackせず停止します。

### 7.2 speciesへ適用する量

| species | 上書き |
| --- | --- |
| 最初の負`reservoir_face` | 有効密度、branch依存の法線cutoff |
| 最初の正`reservoir_face` | `sheath_reference_coordinate`指定時に局所ion密度、冷たいbeam法線速度 |
| 最初の負`photo_raycast` | 自由photoelectron電流密度、branch依存の放出cutoff、normal drift 0 |

`sheath_reference_coordinate`を指定すると、その平面をsheath座標$z_s=0$とし、共有reservoir faceまでの
距離でZhao 1D profileをsampleします。局所$\phi(z_s)$からelectron free/reflected population、
photoelectron free/captured population、ion速度を再構成します。この経路では既に局所VDFを構成しているため、
汎用のbarrier energy shiftをもう一度適用しません。

### 7.3 重要な非保証

legacy Zhao profileの$E(z)$はparticle pusherで用いる電場へ加算されません。任意のsurface geometryや
時間変化する`q_elem`に自己整合させる外部場でもありません。これは文献モデルに基づく
injection/photoemission closureです。軌道と同じ外部電位profileを必要とする場合は
`field.model="kinetic_1d"`を使います。そこで`kinetic_closure="zhao_charge_driven"`を選ぶと、
同じZhao解析densityを外部Poisson closureへ使い、surface chargeと光電子populationから毎batch更新します。
split windowを置けないrough surfaceで線形性条件を満たす場合だけ、高度な
`field.model="unified_linear_response"`を選びます。

## 8. `floating_no_photo`とtracked photoelectron

`floating_no_photo`はelectron Maxwellianのcutoff後fluxとion inflow fluxが一致する負の浮遊電位を
二分法で解き、electron reservoirへcutoffを適用します。空間的なsheath profileは作りません。

光電子は放出時の重みを保ったまま通常粒子として追跡します。通常open面では
`ordinary_open.model="potential_barrier"`、外部sheathでは`field.model`と
`particles.mode`から導出されるinterface transportがreturn / escapeを決めます。

## 9. 互換性と推奨用途

| 目的 | 推奨構成 |
| --- | --- |
| **標準:** 無限周期レゴリス + 自己整合1D sheath | `cached_kneq0` + `field.model="kinetic_1d"` + `particles.mode="same_batch"` |
| **高度:** rough surfaceでsplit windowがない線形検証 | `field.model="unified_linear_response"` + `particles.mode="same_batch"` |
| 有限画像のface scalar障壁 | `particles.inflow_model="infinity_barrier"` |
| legacy Zhao文献closureとの比較 | `inflow_model="legacy_sheath"` + `legacy_sheath_model="zhao_*"` |
| charge-driven Zhao外部シース | `field.model="kinetic_1d"` + `kinetic_closure="zhao_charge_driven"` |
| 光電子なしの簡易電流釣合い | `legacy_sheath_model="floating_no_photo"` |

tracked 1D profileは`inflow_model="auto"`を要求し、scalar barrierやlegacy Zhao流入を重ねません。
同じ電位差やcutoffを二重に適用しないためです。1D returnと3D explicit orbitは現在`b0=0`を要求します。

## 10. 実装対応

| 処理 | Fortran実装 |
| --- | --- |
| reservoir flux・速度sample | `src/particles/bem_injection.f90` |
| face平均電位・energy shift設定 | `src/config/bem_app_config_runtime.f90` |
| Zhao runtime coupling | `src/physics/sheath/bem_sheath_runtime.f90` |
| Zhao core/root/profile | `src/physics/sheath/bem_sheath_model_core.f90` |
| 1D outer return | `src/physics/outer_plasma/bem_outer_plasma_interface.f90` |
| unified 3D orbit | `src/physics/outer_plasma/bem_outer_plasma_orbit.f90` |
| open-face scalar reflection | `src/runtime/simulator/bem_particle_stepper.f90` |

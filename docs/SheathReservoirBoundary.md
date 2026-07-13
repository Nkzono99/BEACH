title: 外部シースとreservoir粒子境界

Lang: [日本語](SheathReservoirBoundary.md) | [English](SheathReservoirBoundary.en.md)

# 外部シースとreservoir粒子境界

この文書は、`reservoir_face`からの流入、電位差による加減速、流入不能粒子、外向き粒子の
反射・return・escape、およびZhao系シース注入補正を説明します。zero modeとouter Poisson solverは
[periodic2 zero modeと外部プラズマ](PeriodicZeroModeOuterPlasma.md)を参照してください。

## 1. 「シース」と呼ばれる機能の区別

BEACHには目的の異なる複数のreduced modelがあります。

| 機能 | 入力 | 変更するもの | 空間的な$E(z)$をparticle pusherへ与えるか |
| --- | --- | --- | --- |
| `reservoir_potential_model="infinity_barrier"` | 注入面平均電位と`phi_infty` | 流入flux、法線速度cutoff、face速度 | いいえ |
| `sheath_injection_model="zhao_*"` | 背景species、太陽高度、光電子基準密度 | reservoir密度/VDF cutoff、ion drift、photoemission current | いいえ |
| `floating_no_photo` | electron/ion流入flux | electron cutoff | いいえ |
| `outer_plasma.model="kinetic_1d"` | 無限遠VDF、surface zero-mode field | 自己整合1D電位・電場・密度profile | outer領域で使用 |
| `unified_linear_response` | surface field、accessible area、線形plasma応答 | localからfarまでのzero/nonzero field | はい |
| `open_boundary_model="potential_barrier"` | 境界通過点電位と`phi_infty` | 外向き粒子のescape/反射 | いいえ |

Zhao系を`kinetic_1d`の別解法として扱ってはいけません。Zhao系は既存speciesの注入分布を
事前補正するモデルで、simulation中のsurface chargeから外部Poisson問題を毎batch解くものではありません。

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

## 4. face平均電位によるlegacy `infinity_barrier`

### 4.1 face電位

各batch冒頭のrefresh済みelectrostatic snapshotを使い、注入口矩形を
`injection_face_phi_grid_n x injection_face_phi_grid_n`のcell-centered点で評価します。

$$
\bar\phi_f=\frac{1}{N^2}\sum_{a,b}\phi(\mathbf x_{ab})
$$

point/`triangle_p0` kernel、periodic2、physical `k=0`、outer state、`e0`はsnapshotと同じ規約です。
この**面平均scalar電位**から$B=2q(\bar\phi_f-\phi_\infty)/m$を作ります。

### 4.2 何を表していないか

このmodelはfaceまでの電位差だけを使い、途中の$E(z)$、turning位置、flight time、空間電荷を解きません。
横方向にface電位が大きく変化していても平均値一つへ縮約します。したがって小さな注入口または
ほぼ一様なfaceのlegacy/reduced modelであり、rough-surface sheathのproduction既定にはしません。

## 5. outer profileによる流入

`particle_transfer_mode="electrostatic_1d_instant_return"`では、z-high `reservoir_face` speciesを
無限遠VDFとして解釈します。

- `kinetic_1d`は各batchで更新済みの`phi_interface-phi_infinity`から$B$を計算します。
- linear-Debye referenceはzero-mode surface chargeからinterface電位差を作ります。
- 到達可能な上流粒子をsampleし、同じエネルギー式でinterface速度へ写像します。

`kinetic_1d`のPoisson solveではsurface zero modeが与える`interface_field`をNeumann条件に使います。
ただし流入速度の写像に直接使うのは$E_I$ではなく、解いたprofileの電位差$\phi_I-\phi_\infty$です。
電場はprofileを決め、電位差が粒子のエネルギー変化を決めます。

## 6. 外向き粒子のescapeとreturn

### 6.1 legacy `open_boundary_model="potential_barrier"`

open face通過点の電位を$\phi_b$、外向き法線速度を$v_n$とすると、無限遠へ必要な障壁は

$$
U_b=q(\phi_\infty-\phi_b)
$$

です。$U_b>\tfrac12mv_n^2$なら法線速度だけを反転し、それ以外はescapeとします。
cornerで複数open faceを同時に横切る一般化はせず、legacy/experimental modelとしてfail closedです。
一様`e0`は無限遠電位を定義できないため、このlegacy outflow barrierには含めません。

### 6.2 linear-Debye instant return

interfaceを外向きに通過した粒子について

$$
v_{n,\infty}^2=v_{n,I}^2+\frac{2q(\phi_I-\phi_\infty)}{m}
$$

を評価します。非負ならescape、負なら指数Debye profile内にturning pointがあります。
解析式で往復時間を求め、return時に法線速度を反転します。接線位置は$\mathbf v_t\tau_\mathrm{outer}$だけ
進め、x/yをperiodic wrapします。

### 6.3 kinetic-profile return

`kinetic_1d_profile_return`は同じ判定を、実際に収束した離散$\phi(z)$で行います。
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
6. 往復時間$\tau_\mathrm{outer}$後に相当するinterface復帰状態を直接構成する。

復帰時は法線速度だけを反転し、接線速度は維持します。接線位置を
$\mathbf v_t\tau_\mathrm{outer}$だけ進めてx/yをperiodic wrapし、interface crossingを検出したlocal stepの
remaining `dt`だけを通常stepperで再積分します。outer flightをglobal simulation timeへ加算しません。

したがってこのmodelは「interfaceで速度とscalar障壁だけを比較する即時鏡面反射」より詳細ですが、
「外部軌道を時間刻みで明示追跡する」ものではありません。後者が必要な場合は
`unified_linear_response + electrostatic_3d_explicit_orbit`を使います。

#### 即時帰還を選ぶ理由と適用範囲

`instant_return`の「即時」は、outer flightを粒子状態の写像には使う一方、global simulation timeを
進めないという意味です。turning粒子はinterfaceを出た時刻と同じsimulation時刻に戻り、
outward/returned chargeも同じbatchに計上されます。`dt`や`batch_duration`だけ待ってから戻す処理ではありません。

これは定常または準定常sheathを粒子領域から消去したreduced closureです。静電・無衝突の定常profileでは、
粒子がescapeするか戻るかは全エネルギーで決まり、十分に定常化した表面の平均return currentは
個々の往復時間に依存しません。往復時間によるouter粒子の滞在密度は`kinetic_1d`の
outgoing/returning density closureがPoisson solveへ含めます。このため、定常電位、長時間平均の電流収支、
定常化後の離脱力を主目的とする計算では即時帰還を標準とします。

一方、UV照射開始、plasma条件の急変、短時間pulseなど、往復時間と同程度以下で場や放出電流が変わる過渡では、
実際のreturn currentは過去のoutgoing currentに依存します。現行modelはその遅延、outer領域に一時滞在する
正味電荷、遅延に伴うovershootや振動を再現しません。したがって過渡電流や立ち上がり時間の評価には使わず、
初期過渡を除いた準定常結果を評価してください。

準定常性は

$$
\epsilon_\mathrm{ad}=\frac{\tau_\mathrm{outer}}{\tau_\mathrm{field}}
$$

で判定します。`field_evolution_timescale`が$\tau_\mathrm{field}$、
`max_frozen_field_ratio`が許容する$\epsilon_\mathrm{ad}$の上限です。これは数値刻み`dt`ではなく、
表面電位やouter profileが変化する物理時間との比較です。また
$\tau_\mathrm{outer}/\mathrm{batch\_duration}\gtrsim1$ならbatch単位のreturn currentは時間的に正しくありません。
その場合でも$\epsilon_\mathrm{ad}\ll1$で十分に定常化した後の長時間平均は利用できますが、batch履歴を
物理的な過渡応答として解釈してはいけません。persistent delayed-return queueは未実装であり、
`outer_queue_enabled=true`は拒否します。

### 6.4 unified 3D orbit

`electrostatic_3d_explicit_orbit`はzero modeとscreened nonzero tailを合成した3D電場中を
固定刻みvelocity-Verletで追跡します。ownership面へ戻ればlocalへ返し、far planeを外向きに
横切ればescapeです。エネルギー誤差、flight time、frozen-field ratio、step上限を検査します。

## 7. Zhao系シース注入補正

### 7.1 解いているもの

Zhao系はsolar-wind electron、ion、photoelectronの解析密度closureと電流条件から、
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

Zhao profileの$E(z)$はparticle pusherのfield snapshotへ加算されません。任意のsurface geometryや
時間変化する`q_elem`に自己整合させる外部場でもありません。Zhao系は文献モデルに基づく
injection/photoemission closureであり、軌道と同じ外部電位profileを必要とするproduction計算では
`kinetic_1d`または`unified_linear_response`を使います。

## 8. `floating_no_photo`とphotoelectron reduced closure

`floating_no_photo`はelectron Maxwellianのcutoff後fluxとion inflow fluxが一致する負の浮遊電位を
二分法で解き、electron reservoirへcutoffを適用します。空間的なsheath profileは作りません。

`photo_escape_model="boltzmann_cutoff"`は放出元自己項を除いた局所電位からescape fraction

$$
f_\mathrm{escape}=\exp\left[-\frac{|q|\max(\phi_\mathrm{emit}-\phi_\infty,0)}{k_BT_\mathrm{PE}}\right]
$$

だけを求めます。returning particleの軌道・再吸収位置・遅延は追跡しない即時reduced closureです。

## 9. 互換性と推奨用途

| 目的 | 推奨構成 |
| --- | --- |
| 無限周期レゴリス + 自己整合1D sheath | `cached_kneq0` + `kinetic_1d` + `kinetic_1d_profile_return` |
| rough surfaceでsplit windowがない線形検証 | `unified_linear_response` + explicit 3D orbit |
| 過去のface scalar障壁を再現 | `infinity_barrier` |
| Zhao文献closureとの比較 | `zhao_*`を単独で使用 |
| 光電子なしの簡易電流釣合い | `floating_no_photo` |

`kinetic_1d_profile_return`は`reservoir_potential_model`およびZhao系との併用を拒否します。
同じ電位差やcutoffを二重に適用しないためです。1D instant returnと3D explicit orbitは現在`b0=0`を要求します。

## 10. 実装対応

| 処理 | Fortran実装 |
| --- | --- |
| reservoir flux・速度sample | `src/particles/bem_injection.f90` |
| face平均電位・energy shift設定 | `src/config/bem_app_config_runtime.f90` |
| Zhao runtime coupling | `src/physics/sheath/bem_sheath_runtime.f90` |
| Zhao core/root/profile | `src/physics/sheath/bem_sheath_model_core.f90` |
| 1D outer return | `src/physics/outer_plasma/bem_outer_plasma_interface.f90` |
| unified 3D orbit | `src/physics/outer_plasma/bem_outer_plasma_orbit.f90` |
| legacy open-face reflection | `src/runtime/simulator/bem_particle_stepper.f90` |

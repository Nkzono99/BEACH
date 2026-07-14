title: シース注入closure

Lang: [日本語](SheathInjectionClosures.md) | [English](SheathInjectionClosures.en.md)

# シース注入closure

`sim.sheath_injection_model`は、解析的なsheath closureから`reservoir_face`と`photo_raycast`の密度、drift、
cutoffを決めます。補正結果はsource samplingへ適用し、生成後の粒子は、別に構成されてbatch内で固定された電場中を進みます。

| model | 解く量 | particle sourceへ反映する量 |
| --- | --- | --- |
| `none` | なし | 設定したVDFをそのまま使用 |
| `floating_no_photo` | 光電子なしの負の浮遊電位 | electronの法線cutoff |
| `zhao_a/b/c` | 指定したZhao branchの解析closure | electron/ion/photoelectronの密度、drift、cutoff、放出電流 |
| `zhao_auto` | 太陽高度に応じてZhao branchを探索 | 収束したbranchと同じ補正 |

自己整合なouter Poisson profileは[kinetic 1D外部プラズマ](KineticOuterPlasma.html)、rough surfaceを含む
線形fieldは[unified linear response](UnifiedLinearResponse.html)で説明します。

## sourceの役割からspeciesを選ぶ

有効なspeciesを設定順に走査し、次を使います。

- 最初の負電荷`reservoir_face`: solar-wind/ambient electron。
- 最初の正電荷`reservoir_face`: ion。
- Zhao系では最初の負電荷`photo_raycast`: photoelectron。

electronとionのreservoirは、同じ`inject_face`を使う必要があります。inward driftは各sourceの
`drift_velocity`から取得するか、`sheath_electron_drift_mode`と`sheath_ion_drift_mode`で決定します。結果は正の内向き速度でなければなりません。
Zhao系は正のelectron/photoelectron温度、ion密度、photoelectron基準密度、ion drift、粒子質量を要求します。

## `floating_no_photo`でelectron/ion流入を釣り合わせる

ionの流入flux$\Gamma_i$を通常のflux-weighted drifting Maxwellianから計算します。負のtrial電位$\phi_0\le0$に対する
electron cutoffは

$$
v_{e,\min}(\phi_0)
=\sqrt{\frac{2e\max(0,-\phi_0)}{m_e}}
$$

です。cutoff後のelectron fluxとion fluxの差

$$
F(\phi_0)=\Gamma_e(v_n\ge v_{e,\min})-\Gamma_i
$$

が0になる負の$\phi_0$を二分法で解きます。初期下限はelectron温度[eV]の128倍を負側に取り、解を挟み込めるまで
広げます。反復回数の上限は80回です。正のion fluxがない場合や、負の解を挟み込めない場合は停止します。

得た$\phi_0$は最初のelectron reservoirの$v_{\min}$へ変換されます。ion VDFは変えず、光電子も扱いません。
空間的な$\phi(z)$、$E(z)$、turning point、flight timeは作らない簡易current-balance closureです。

## Zhao closureの無次元量を作る

solar elevationを$\alpha$、reference photoelectron densityを$n_{\mathrm{phe,ref}}$とすると、surface source densityは

$$
n_{\mathrm{phe},0}=n_{\mathrm{phe,ref}}\sin\alpha
$$

です。実装はさらに

$$
v_{\mathrm{swe,th}}=\sqrt{\frac{2eT_\mathrm{swe}}{m_e}},
\qquad
v_{\mathrm{phe,th}}=\sqrt{\frac{2eT_\mathrm{phe}}{m_e}},
$$

$$
c_s=\sqrt{\frac{eT_\mathrm{swe}}{m_i}},
\quad
M=\frac{u_i}{c_s},
\quad
u=\frac{u_e}{v_\mathrm{swe,th}},
\quad
\tau=\frac{T_\mathrm{swe}}{T_\mathrm{phe}}
$$

を作り、解析densityとcurrent条件の非線形rootを解きます。未知量はsurface電位$\phi_0$、必要なbranchでの
potential minimum $\phi_m$、有効solar-wind electron密度$n_{\mathrm{swe},\infty}$です。

## 太陽高度とrootからZhao branchを決める

| branch | potentialとpopulation構造 |
| --- | --- |
| A | $\phi_0>0$、$\phi_m<0$の非単調profile。minimum下側にcaptured photoelectron、上側にreflected solar-wind electron |
| B | $\phi_0>0$の単調profile。photoelectron captureを含み、solar-wind electron reflectionなし |
| C | $\phi_0<0$の単調profile。solar-wind electron reflectionを含み、photoelectron cutoffなし |

`zhao_auto`は$\alpha<20^\circ$でC→A→B、それ以外でA→B→Cの順にrootを試します。この探索は、
Zhao family内でbranchを選ぶためのものです。すべて失敗すれば停止し、`floating_no_photo`やouter Poisson modelへ切り替えません。
`zhao_a/b/c`を明示した場合は指定branchだけを解き、失敗時に別branchへ移りません。

## closureの解を各sourceへ反映する

### ambient electron

有効密度を$n_{\mathrm{swe},\infty}$へ置き換え、branchに応じて

$$
v_{e,\min}=\sqrt{\frac{2e\max(0,\Delta\phi_e)}{m_e}},
\qquad
\Delta\phi_e=
\begin{cases}
-\phi_m & A,\\
0 & B,\\
-\phi_0 & C
\end{cases}
$$

を使います。

### photoelectron

放出法線cutoffは

$$
v_{pe,\min}=\sqrt{\frac{2e\max(0,\Delta\phi_{pe})}{m_{pe}}},
\qquad
\Delta\phi_{pe}=
\begin{cases}
\phi_0-\phi_m & A,\\
\phi_0 & B,\\
0 & C
\end{cases}
$$

です。normal driftを0へ上書きし、自由photoelectron電流密度を

$$
J_{pe}=\frac{|q_{pe}|n_{\mathrm{phe},0}v_{\mathrm{phe,th}}}{2\sqrt\pi}
\begin{cases}
\exp[(\phi_m-\phi_0)/T_\mathrm{phe}] & A,\\
\exp[-\phi_0/T_\mathrm{phe}] & B,\\
1 & C
\end{cases}
$$

へ置き換えます。$T$と$\phi$はこの指数ではeV/Vの対応する規約です。この電流密度からrayのmacro重みを決める方法は
[光電子の放出とライフサイクル](PhotoelectronEmission.html)で説明します。

### ionとlocal reference plane

`sheath_reference_coordinate`を指定すると、その平面をsheath座標$z_s=0$とします。この平面から共有reservoir faceまでの
距離において、Zhao 1D profileをsampleします。局所$\phi(z_s)$からfree/reflected electron、free/captured photoelectron、ionの速度と
密度を再構成し、electron density/cutoffとion cold-beam法線速度に反映します。

この経路はすでにlocal VDFを構成しているため、汎用の`reservoir_potential_model`によるbarrier energy shiftを
重ねません。

## Zhao closureを使う範囲

Zhao profileから得る量は、文献closureに基づくsource VDFの事前補正です。Zhao rootはBoris pusherで用いる電場とは独立しており、
batchごとの`q_elem`からは更新されません。したがって、任意の3D surface geometryと自己整合な外部場を必要とする構成には対応しません。

sourceと外向き粒子が同じ自己整合potential profileを共有する計算には、
`kinetic_1d + kinetic_1d_profile_return`または検証済みのunified構成を使います。

## 一つのsource補正だけを適用する

- `sheath_injection_model`は`reservoir_potential_model`との併用を拒否します。
- `velocity_distribution="grid"`のreservoirは現行Zhao/floating補正と併用できません。
- `kinetic_1d_profile_return`はZhao系と`reservoir_potential_model`を拒否します。
- Zhaoは負electron、正ion、負photoelectronの3 speciesを要求します。
- `floating_no_photo`は負electronと正ionだけを使い、photoelectron sourceを補正しません。

`photo_escape_model="boltzmann_cutoff"`はZhaoとは別の局所reduced closureです。式とtracked returnとの排他関係は
[光電子の放出とライフサイクル](PhotoelectronEmission.html)にまとめています。

## Code reference

- floating root、Zhao density/current/root solve: [`bem_sheath_model_core.f90`](../src/physics/sheath/bem_sheath_model_core.f90)
- species検出とruntime override: [`bem_sheath_runtime.f90`](../src/physics/sheath/bem_sheath_runtime.f90)
- source samplingへの適用: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- 入力と組合せの検証: [`bem_app_config_parser.f90`](../src/config/app_config_parser/bem_app_config_parser.f90)

title: reservoir注入

Lang: [日本語](ReservoirInjection.md) | [English](ReservoirInjection.en.md)

# reservoir注入

`source_mode="reservoir_face"`は、box外の分布関数から、指定面を通過して流入するmacro粒子を生成します。
上流分布のうちfaceへ到達できる粒子の流入量を積分し、同じ電位差で法線速度をface上へ写してから
位置と速度をsampleします。

## 開口の位置と流入面積を決める

`inject_face`は`x_low`、`x_high`、`y_low`、`y_high`、`z_low`、`z_high`のいずれかです。
`pos_low`と`pos_high`はそのbox面上の矩形開口を定めます。面の内向き単位法線を$\mathbf n$、
開口面積を$A$とします。

`reservoir_face`は`sim.use_box=true`、`sim.batch_duration>0`を要求します。開口の法線座標はbox面と一致し、
2つの接線方向は正の面積を持つ必要があります。

## Maxwell分布を流入fluxで重み付けする

上流のdrifting Maxwell分布で、法線driftを$u_n=\mathbf u\cdot\mathbf n$、
熱速度標準偏差を$\sigma=\sqrt{k_\mathrm{B}T/m}$とします。上流法線速度の下限が$v_{\min}$なら、
片側流入fluxは

$$
\Gamma_\mathrm{in}
=n\int_{v_n\ge v_{\min}}v_n f_n(v_n)\,dv_n
=n\left[u_n\{1-\Phi(a)\}+\sigma\varphi(a)\right],
\qquad
a=\frac{v_{\min}-u_n}{\sigma},
$$

です。$\Phi$と$\varphi$は標準正規分布のCDFと密度です。温度0では、$u_n\ge v_{\min}$のとき
$\Gamma_\mathrm{in}=nu_n$、それ以外は0です。

接線速度は通常のGaussian、法線速度は$v_n f_n(v_n)$に比例するflux-weighted half-range分布からsampleします。
この違いが、体積中のMaxwell分布をそのままfaceへ置く方法との重要な差です。

## 端数を持ち越して長時間のfluxを保つ

実粒子の期待流入数とmacro粒子の期待数は

$$
N_\mathrm{phys}=\Gamma_\mathrm{in}A\,\Delta t_\mathrm{batch},
\qquad
N_\mathrm{macro,expected}=\frac{N_\mathrm{phys}}{w}
$$

です。整数化で捨てた端数を$r$として次batchへ持ち越し、

$$
N_\mathrm{macro}=\left\lfloor r+N_\mathrm{macro,expected}\right\rfloor,
\qquad
r\leftarrow r+N_\mathrm{macro,expected}-N_\mathrm{macro}
$$

とします。したがって各batchの個数は揺れても、多数のbatchを通じた総流入量は期待fluxに一致します。

重みは正の`w_particle`を直接与えるか、`target_macro_particles_per_batch`から逆算します。
後者はsampling数を制御する設定であり、物理fluxを変更しません。`-1`はspecies 1で解決した重みを共有します。
両方を同時に指定することはできません。

## velocity gridの値をfluxへ変換する

`velocity_distribution="grid"`では、CSVの速度点と非負の値$f$からsampleします。

| `velocity_grid_pdf_kind` | CSVの$f$の解釈 | sampling重み |
| --- | --- | --- |
| `phase_space` | 位相空間分布 | $\max(v_n,0)f$ |
| `flux_weighted` | すでに流入粒子へ重み付け済み | $f$ |

いずれも$v_n\ge v_{\min}$かつ$v_n>0$の点だけを使います。`velocity_grid_sampling="rectilinear"`は完全な直交gridを補間し、
`discrete`は入力点から離散的にsampleします。`auto`は、可能な場合に`rectilinear`を使います。

grid分布では粒子数を`particle_flux_m2_s`または`current_density_a_m2`から決めます。密度・温度は使わず、
現行実装はZhao注入補正との併用を拒否します。

## 一つの電位差で到達条件とface速度を決める

無限遠または上流の電位を$\phi_\infty$、face/interface電位を$\phi_f$とします。1D静電写像は接線速度を保ち、
法線エネルギーを保存します。

$$
\frac12m v_{n,f}^2+q\phi_f
=\frac12m v_{n,\infty}^2+q\phi_\infty,
$$

$$
v_{n,f}^2=v_{n,\infty}^2-B,
\qquad
B=\frac{2q(\phi_f-\phi_\infty)}{m}.
$$

| $B$ | accessibilityとface速度 |
| ---: | --- |
| $B>0$ | $v_{n,\infty}\ge\sqrt B$だけが到達し、faceまでに減速 |
| $B=0$ | 法線速度を変更しない |
| $B<0$ | すべての流入粒子が到達でき、faceまでに加速 |

上流VDFは$v_{\min}=\sqrt{\max(B,0)}$で到達可能な軌道に切り、採用した速度を同じ$B$でface速度へ
写します。これにより、生成する粒子数とface到達時の速度が一つの電位差から決まります。

## 電位差を供給する外部モデルを選ぶ

| 構成 | $\phi_f-\phi_\infty$の取得 | 空間的な外部軌道 |
| --- | --- | --- |
| 補正なし | 0 | なし |
| `reservoir_potential_model="infinity_barrier"` | 注入開口の平均電位と`sim.phi_infty` | 解かない |
| split `linear_debye` | surface zero modeから作るinterface電位差 | return時は解析profile |
| split `kinetic_1d` | 収束したouter profileのinterface電位 | return時は離散profile |
| Zhao系注入closure | branchから局所VDFとcutoffを構成 | pusherへ$E(z)$を加えない |

legacy `infinity_barrier`は、batch開始時に構成した電場・電位を開口内の
`injection_face_phi_grid_n`四方のcell-centered点で評価し、scalar平均$\bar\phi_f$を使います。
point/triangle kernel、periodic場、zero mode、outer state、`sim.e0`を含む場と同じ電位規約ですが、
faceまでの途中の$E(z)$、turning位置、flight time、空間電荷は解きません。

split outer modelではz-highの`reservoir_face`を無限遠VDFとして解釈します。`kinetic_1d`ではinterface電場を
Poisson境界条件に使いますが、粒子速度を変える量は解から得た$\phi_I-\phi_\infty$です。[外部プラズマモデル](OuterPlasmaModels.html)で、
modelの選択とreturn側への同じprofileの適用を説明します。

Zhao branchと`floating_no_photo`がsource密度・cutoff・driftを上書きする規則は
[シース注入closure](SheathInjectionClosures.html)にまとめています。

## 生成位置の位相をjitterで散らす

位置は開口矩形内で一様sampleし、法線方向へ$10^{-12}$ mだけbox内部へ移します。
`position_jitter_dt>0`なら各粒子に$[0,\texttt{position\_jitter\_dt})$の一様な仮想時間を与え、
$\mathbf x\leftarrow\mathbf x+\mathbf v\tau$とします。その後、両面periodicな軸はprimary boxへwrapし、
その他の軸はbox内へclampします。

jitterは同じbatchに全粒子が厳密に同じface位置から始まる人工的な位相を散らします。変更するのは
生成位置だけで、global simulation timeと個々の追跡可能時間は保ちます。

## globalな端数をMPIとrestartで引き継ぐ

MPIではglobal個数と端数をrootで決めてから、各rankへ分配します。globalな端数は`macro_residuals.csv`へ保存し、
restart時に全rankへbroadcastします。確認するのはspecies別の注入個数と注入電荷、解決後のmacro重み、端数、
適用した$v_{\min}$と電位差、batch平均の吸収電流とescape電流です。

`batch_duration`を変えると1 batchの期待粒子数と場の更新間隔が同時に変わります。
[batch幅と安定性](BatchDurationStability.html)に、物理的な定常値の確認方法を示しています。

## Code reference

- flux積分、macro個数、Maxwell/grid sampling: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- barrierとouter profileのvelocity補正: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- geometryと入力組合せの検証: [`bem_app_config_parser_validate.f90`](../src/config/app_config_parser/bem_app_config_parser_validate.f90)
- Zhao species補正: [`bem_sheath_runtime.f90`](../src/physics/sheath/bem_sheath_runtime.f90)
- MPI-global macro残差: [`bem_restart.f90`](../src/runtime/bem_restart.f90)

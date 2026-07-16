title: kinetic 1D外部プラズマ

Lang: [日本語](KineticOuterPlasma.md) | [English](KineticOuterPlasma.en.md)

# kinetic 1D外部プラズマ

`outer_plasma.model="kinetic_1d"`は、periodic表面の上にある外部プラズマを、平面平均された1Dの静電・無衝突・
非磁化profileとして解きます。表面電荷が与えるinterface電場と無限遠の粒子VDFを接続し、interface電位、密度、
電流を自己整合に求めます。

収束した外部profileは双方向に使われます。無限遠VDFからinterfaceへ入る粒子の写像は
[`reservoir_face` の流入量と速度サンプリング](ReservoirInjection.html)、interfaceから出る粒子のescape/return写像は
[粒子のescapeとreturn](ParticleEscapeReturn.html)に対応します。

## local領域とouter領域をinterfaceで分ける

local粒子領域とouter領域はz-highのownership interface $z=z_I$で接続されます。

| 領域 | 使用する場 |
| --- | --- |
| meshからinterfaceまで | periodic $k\ne0$ surface fieldとsurface $k=0$ field |
| interfaceより外側 | `kinetic_1d`の平面平均profile |

outer profileをlocal領域の全点へ重ねません。このsplit構成では、interfaceまでに横方向modeが十分減衰すると仮定し、
interfaceで電位と法線電場を接続します。roughnessとplasma responseを同じ領域で解く必要がある場合は
[unified linear response](UnifiedLinearResponse.html)を使います。

## interface電場と無限遠VDFから電位を決める

未知量は、interfaceから無限遠reservoirまでの平面平均電位$\phi(z)$です。surface zero modeから得た
interface電場$E_I$をNeumann条件

$$
-\phi'(z_I)=E_I
$$

として与え、粒子VDFから$\rho(\phi)$を構成して

$$
\frac{d^2\phi}{dz^2}=-\frac{\rho(\phi)}{\epsilon_0}
$$

を解きます。無限遠電位はgaugeとして$\phi_\infty=0$です。interface電位$\phi_I=\phi(z_I)$は入力ではなく、
interface電場、VDFに基づく密度モデル、far条件を同時に満たす解から決まります。

Newton法は、Poisson方程式と境界条件からなる非線形方程式を解きます。batch間の表面電荷更新によって$E_I$が変わると、
$\phi_I$と各speciesの電流も変わります。electron、ion、photoelectron、外部回路のcurrent densityは、収束profileから得られる
診断量です。

## VDFを電位依存の電荷密度へ写す

z-highの最初の負・正`reservoir_face` speciesを、それぞれ無限遠ambient electronとionとして使います。
`photoelectron_closure="kinetic_mean"`では、最初の負電荷`photo_raycast` speciesの温度と放出電流密度を
平面平均photoelectron sourceとして加えます。

| population | 与える量 | outer密度の構成 |
| --- | --- | --- |
| ambient electron | $n_{e,\infty},T_e,q_e,m_e$ | half-Maxwellianを全エネルギー保存で写像し、吸収軌道とpotential-reflected軌道を含める |
| ion | $n_{i,\infty},T_i,q_i,m_i,u_{i,\infty}$ | cold beamをエネルギー保存とflux保存で写像 |
| photoelectron | $T_{pe},q_{pe},m_{pe},\Gamma_{pe,0}$ | surface half-Maxwellianのoutgoingとturning後のreturning平均population |

cold ion密度モデルは

$$
u_i(z)=\sqrt{u_{i,\infty}^2-\frac{2q_i\phi(z)}{m_i}},
\qquad
n_i(z)=n_{i,\infty}\frac{u_{i,\infty}}{u_i(z)}
$$

です。平方根の中が負になるprofileは、ionがinterfaceへ到達できないため拒否します。

photoelectronの無限遠escape率は

$$
f_{pe,\mathrm{esc}}
=\exp\left[-\frac{\max\{0,q_{pe}(\phi_\infty-\phi_I)\}}{T_{pe}}\right]
$$

で、残りがreturn populationです。ここで$T$は実装内部ではJ単位です。Poisson sourceは

$$
\rho(\phi)=q_en_e(\phi)+q_in_i(\phi)+q_{pe}n_{pe}(\phi)
$$

です。各解析密度モデルは密度に加え、Newton Jacobian用の$\partial n_s/\partial\phi$と必要な
$\partial n_s/\partial\phi_I$を返します。

`kinetic_mean`のoutgoing/returning densityは定常outer空間電荷の平均密度モデルです。tracked粒子の表面depositを置き換えず、
統計的return currentを別途表面へ加えません。[光電子の放出とライフサイクル](PhotoelectronEmission.html)で、
source電荷とtracked再吸収の収支を説明します。

## 有限gridをRobin tailで無限遠へ接続する

内部点は伸長可能なnonuniform 1D grid上のconservative finite-volume residualです。現行runtime値は次のとおりで、
`debye_length`以外はまだ個別のinput parameterではありません。

| 項目 | 現行値 |
| --- | ---: |
| grid点数 | 128 |
| 計算領域長 | $10\lambda_D$ |
| grid stretch | 2 |
| Newton最大反復 | 40 |
| residual tolerance | $10^{-8}$ |

有限grid上端$L$より先を$\lambda_\mathrm{tail}=\lambda_D$の指数tailとみなし、

$$
\phi'(L)+\frac{\phi(L)}{\lambda_\mathrm{tail}}=0
$$

というRobin条件を課し、無限遠gaugeへ指数緩和させます。grid外のtailは、return particleのflight-time積分にも使います。

## continuation付きNewton法で物理解を追う

解析density derivativeからbordered-tridiagonal Jacobianを作るため、1 Newton stepは$O(N_z)$です。
$\phi_I=\phi_1$への内部densityの依存が、通常のtridiagonal stencilへborder列を加えます。

1. Newton stepを計算する。
2. backtrackingでtrial stateをsupported monotone branch内に制限する。
3. 通常Newtonが停滞すればpseudo-transient regularizationを使う。
4. 前batch profileがある場合、以前のinterface fieldから現在値までcontinuationし、失敗した刻みを半分にする。

regularizationとcontinuationは収束経路だけを変えます。最終判定は必ず元の離散Poisson残差へ戻ります。

## supported sheath branchだけを受理する

| 条件 | 必要な内容 |
| --- | --- |
| original Poisson residual | regularized式ではなく元の残差がtolerance以下 |
| monotone branch | supported electron-repelling profileから外れない |
| ion accessibility | 全grid点で$u_i^2(z)>0$ |
| kinetic Bohm entry | $u_{i,\infty}\ge\sqrt{(T_e+\gamma_iT_i)/m_i}$ |
| infinity quasineutrality | $q_en_{e,\infty}+q_in_{i,\infty}\simeq0$ |

virtual cathodeを持つ非単調profile、trapped population、sub-Bohm ion inflowは適用外です。数値的に収束しても、
これらの条件を満たさない解は受理しません。`not_applicable`、`no_physical_solution`、`numerical_failure`のいずれかのstatusを付けて停止します。

## profileを指定strideで更新し全rankへ共有する

`outer_update_stride`で選ばれたbatchに、commit済みsurface chargeからinterface fieldを再構築してprofileを更新します。
同一model identity・同一gridの前profileだけをNewton初期値にできます。skipするbatchでは以前のouter stateを使いますが、
surface側の電場・電位は現在のcommit済み電荷から更新されます。

MPIではroot rankが1D solveを行い、status、profile、current diagnosticsを全rankへbroadcastします。同じbatchの粒子は
更新後の電場・電位を共有します。これらはbatch中に固定され、粒子1個のhitごとにはouter solveを実行しません。

## residual・電流・帯電を一緒に収束させる

収束した$z,\phi,E,\rho$は`outer_plasma_profile.csv`へ保存し、restart時の初期値にも使います。
最低限確認するのは、`interface_potential`、`interface_field`、`outer_integrated_charge`、species別とtotalのcurrent density、
Newton反復数、original nonlinear residualです。

さらに`debye_length`とinterface位置、必要ならsource samplingを変え、interface電位、電流、表面帯電が収束することを
確認します。returnを有効にした場合は[粒子のescapeとreturn](ParticleEscapeReturn.html)のflight time、
frozen-field ratio、準定常性も検査します。

## Code reference

- VDFに基づく密度モデルとnonlinear Poisson solve: [`bem_outer_plasma_kinetic.f90`](../src/physics/outer_plasma/bem_outer_plasma_kinetic.f90)
- runtime speciesからsolver optionを構成: [`bem_outer_plasma_kinetic_runtime.f90`](../src/runtime/bem_outer_plasma_kinetic_runtime.f90)
- surface fieldとの接続とMPI collective solve: [`bem_electrostatic_snapshot.f90`](../src/physics/bem_electrostatic_snapshot.f90)
- profile output: [`bem_output_writer.f90`](../src/runtime/bem_output_writer.f90)
- profile restart: [`bem_restart.f90`](../src/runtime/bem_restart.f90)

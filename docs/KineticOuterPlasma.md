title: kinetic 1D外部プラズマ

Lang: [日本語](KineticOuterPlasma.md) | [English](KineticOuterPlasma.en.md)

# kinetic 1D外部プラズマ

`outer_plasma.model="kinetic_1d"`は、periodic表面の上にある外部プラズマを、平面平均された1Dの静電・無衝突・
非磁化profileとして解きます。表面電荷が与えるinterface電場と無限遠の粒子VDFを接続し、interface電位、密度、
電流を自己整合に求めます。

BEACHでは、外部reservoirと平均シースを結ぶ標準モデルとして`kinetic_1d`を推奨します。対応する
`return_model`と`particle_transfer_mode`を組み合わせると、同じprofileを粒子の流入とescape/returnにも使います。
外部シースが必要なケースは、rough surface近傍の横方向screeningを解く明確な要件がない限り、このモデルから構成します。

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
interfaceで電位と法線電場を接続します。roughnessとplasma responseを同じ領域で線形に解く必要がある場合だけ、
[高度な粗面線形screening](UnifiedLinearResponse.html)を使います。

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
`photoelectron_density_model="kinetic_mean"`では、最初の負電荷`photo_raycast` speciesの温度と放出電流密度を
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

## Zhao populationを蓄積電荷へ接続する

`kinetic_closure="zhao_charge_driven"`は、free/reflected ambient electron、free/captured photoelectron、cold ionから
なるZhao populationを使い、現在の蓄積表面電荷が与える$E_I$を満たすprofileを解きます。
`zhao_branch="auto"`または`a`、`b`、`c`でbranchを選びます。既定の
`kinetic_closure="absorbing_maxwellian"`は従来の有限grid Poisson modelです。
Zhao closureでは、enabledな負電荷z-high `reservoir_face` ambient electronと正電荷z-high `reservoir_face` ionを
それぞれちょうど1つ要求します。`photoelectron_source_scale>0`では負電荷`photo_raycast` photoelectronも
ちょうど1つ要求し、`photoelectron_source_scale=0`では逆にenabledな光電子sourceを拒否します。
現行実装は`sim.sheath_electron_drift_mode="normal"`と`sim.sheath_ion_drift_mode="normal"`だけを受理し、
`photo_raycast.normal_drift_speed=0`およびcold-ion近似$T_i\le0.1T_e$を要求します。

charge-driven Zhaoでは、現在のinterface電場を指定するため、旧Zhaoのzero-current式をrootへ課しません。
Type B/Cでは無限遠準中性と

$$
2\int_{\psi_0}^{0}\hat\rho(\psi)\,d\psi=\hat E_I^2,
$$

Type Aでは無限遠準中性、upper branchのfar-field条件、および

$$
-2\int_{\psi_m}^{\psi_0}\hat\rho_{\mathrm{lower}}(\psi)\,d\psi=\hat E_I^2
$$

を解きます。ここで$\psi=\phi/T_{pe}$、$\hat E_I=E_I\lambda_{D,pe}/T_{pe}$です。旧zero-current式は
species別とtotalのcurrent-density診断へ変換され、帯電途中のtotal currentは非零でも構いません。
過渡population scale $\eta$はphotoelectron密度、無限遠準中性、Sagdeev項をscaleしますが、current診断の
raw photoelectron emission-current項はscaleせず、full tracked sourceのまま使います。
初期の$E_I=0$で強い光電子populationを含む準中性rootが存在しない場合だけ、outer populationが未形成の
ambient-only平坦stateをbranch `0`として使います。最初のtracked currentが電荷を作った後は通常のZhao rootを解き、
非零電場でこのbootstrapへfallbackしません。

UVなしは`outer_plasma.photoelectron_source_scale=0`で選びます。このときZhaoの光電子密度・raw emission currentは
厳密に0で、正規化密度と温度にはambient $n_\infty,T_e$を使います。設定されたambient electron densityは
準中性領域の総密度としてion densityとの準中性を検査し、solverがfree/reflected VDFに必要な入射electron densityを
解きます。reservoir macro数と速度cutoffにもその導出密度と同じprofile写像を使います。平坦な$E_I=0$はType B/Cの
接続点で、負の$E_I$ではType Cへ入ります。

full photoelectron populationの準定常A/B/C branchは、必ずしも$E_I=0$へ連続しません。要求電場がbranchの
可解域外なら、既定の`outer_queue_enabled=false`では`no_physical_solution`で停止します。

`outer_queue_enabled=true`では、追跡中の光電子がouter領域に滞在する間、そのmacro粒子重みをqueue inventoryとして保持します。
全rankの光電子数を水平面積で割ったcolumnをtargetとし、$0\le z\le10\lambda_{D,pe}$で
$n_{pe,f}+n_{pe,c}$を積分したZhao columnが一致するようにpopulation scale $\eta$を解きます。
`outer_photoelectron_population_fraction`という出力名ですが、$\eta$は確率ではなく定常reference populationに対する
occupancy scaleです。solverは$\eta=0$から連続する$0\le\eta\le16$の物理解を探索し、1を超える一時的overshootも
許します。`[0,1]`へclampせず、targetを無視したfull-population解やbranch `0`、disconnected branchへjump/fallbackしません。
queue modeは`zhao_branch="auto"`を要求し、縮退条件を満たす連続的なA/B等のbranch遷移だけを許します。
現在のbisectionはcolumnが$\eta$とともに単調増加するpathだけを受理し、foldを含む連続pathは未対応として停止します。
targetへ至る連結・単調解がない場合は`no_physical_solution`で停止します。

### Zhao continuation失敗を診断する

Zhao continuationがfail closedすると、MPI rootはgeneric errorの前に
`BEACH zhao-continuation`で始まる5行をstderrへ出力します。`call_stage`は`pre_batch`または
`post_enqueue`、`batch`は失敗したbatchです。残りの行にはsolver stage、reason code、underlying/return status、
`attempt`、`attempted_step`、target field・$\eta$・column、直前root、拒否candidate、root residual、
`normalized_potential_jump`、`log_density_jump`、`normalized_root_jump`が入ります。
branchを特定できない場合は`from_branch=-`または`candidate_branch=-`と出力します。実数は3桁指数を含む
full-range scientific formatで記録し、rootが5行をflushした後に全MPI rankが停止します。

公開Fortran procedure `trace_zhao_branch_atlas`は、指定したA/B/C branchの連結curveをpseudo-arclengthで追跡する
診断APIです。production continuation、root選択、fallbackには接続していません。有限density floorやpoint数への到達は
`search_limit`であり、別branchや退化接続を調べずに全Zhao定常解の`target_unreachable`を意味しません。

公開Fortran procedure `diagnose_zhao_ab_degeneracy`は、Type Bの密度ゼロ端を
$q=\sqrt{-\phi_m/T_{pe}}$で診断します。Type Aの準中性curveに沿うfar-field residualの$q^3$係数、
Type A/Bのinterface field積分差、入力Type B rootの準中性・field residual、および有限$q$ probeを
`zhao_ab_degeneracy_diagnostics_type`へ保存します。`regular_connection_conditions_met`は接続の必要条件だけを表し、
独立したType A componentの非存在や動的安定性を保証しません。このAPIもproduction branch選択を変更しません。

公開Fortran procedure `trace_zhao_field_column_homotopy`は、前状態$(E_0,N_0)$とtarget
$(E_1,N_1)$を直線補間し、photoelectron sourceが非零の固定Type B/C branch上でZhao残差と有限長column残差を同時に
pseudo-arclength追跡する診断APIです。accepted point間で$\lambda=1$を挟むと、$\lambda=1$を固定したcorrectorで
targetへ厳密に着地します。`target_reached`はこのcorrectorが収束した場合だけtrueとなり、
`homotopy_fold_detected`は接線の$\lambda$成分の反転を記録します。異なるbranchへの切替、production root選択、
fallbackには接続していません。
非単調Type Aの5座標correctorは現段階では数値的に十分条件付けされていないため明示的に拒否し、
columnが恒等的に0となるno-photo Type Cも退化系として拒否します。Type Aは既存の固定電場atlasとA/B端診断を使います。

強UV fixtureでは、forward Type B atlasを密度floorまで、reverse atlasを$\eta$下限まで追跡しても、
$E_I=0.9072962759\,\mathrm{V/m}$、$L=10\lambda_{D,pe}$のtarget columnを挟みませんでした。密度ゼロlimitでは
準中性curve上の$q^3$係数が約$2.9\times10^{-2}$で0ではないため、そのlimitへ連続するregular local Type A tangentの
必要条件を満たしません。
従って実行時rootを精密化したType B componentではtargetが不可達です。この結論は別の$L$やdisconnected Type A rootを
排除せず、production solverは引き続きfail closedします。数値、適用範囲、dynamic outerへ進む判断は
[ADR 0005](adr/0005-zhao-continuation-and-dynamic-outer.md)に記録しています。

batch 15の正常終了を再生して得たType B stateから、停止したbatch 16 targetまで上記の直線homotopyも追跡しました。
$E_0=0.8424570666\,\mathrm{V/m}$、$N_0=9.3202065681\times10^7\,\mathrm{m^{-2}}$から開始したforward curveは、
$\lambda\simeq0.33179$、$E_I\simeq0.86397\,\mathrm{V/m}$で有限ambient density floorへ達しました。
floorを$\log(n_{e,\infty}/n_{pe,ref})=-27$から$-24$へ変えても終端$\lambda$の差は$10^{-5}$未満で、
density-zero limitへの収束と整合します。全accepted pointでroot/column残差は$5\times10^{-10}$以下、
row-rank indicatorは有限かつ数値拒否閾値$10^{-12}$を上回りました。
target $E_1=0.9072962759\,\mathrm{V/m}$、$N_1=9.9455765203\times10^7\,\mathrm{m^{-2}}$には届かず、
$\lambda$ foldもありません。この端でも準中性$q^3$係数は約$2.9\times10^{-2}$であり、
regular local Type A tangentの必要条件を満たしません。
従って、batch 15 rootから$\lambda$増加方向へ進む固定branchの同時$(E_I,N_{pe})$準定常pathは、target前に有限floorで終端します。
逆向きの大域curve、disconnected component、異なる$L$はこの診断の対象外です。

同じ$0\le z\le10\lambda_{D,pe}$をqueue粒子の有限control volumeとし、
この範囲内でturningしなければ$L=10\lambda_{D,pe}$で外部reservoirへ吸収/escapeします。queue modeでは$L$外の
Robin tailを使ったreturn判定を行いません。各eventでは、`t_due`から最初のbatch-start pollまでの量子化遅延
`delta_poll`と、batch内crossing時刻のmidpoint近似誤差上限`batch_duration/2`を含む
`tau_outer + delta_poll + batch_duration/2`へ
`max_frozen_field_ratio * field_evolution_timescale`の上限を課します。`batch_duration`にも同じ上限を設定時に
要求します。時間離散化、batch末corrector、checkpointの規約は
[粒子のescapeとreturn](ParticleEscapeReturn.html#zhao-過渡closureでouter-flightをqueueする)を参照してください。

[`periodic2_zhao_charge_driven_outer.toml`](../examples/periodic2_zhao_charge_driven_outer.toml)は過渡を解かず、
1 coarse batchで既知のType A可解域へ入るbranch-entry smokeです。
[`periodic2_zhao_no_photo_outer.toml`](../examples/periodic2_zhao_no_photo_outer.toml)は同じambient入力だけで
flat stateからType Cへ帯電するno-photo smokeです。
[`periodic2_zhao_transient_outer.toml`](../examples/periodic2_zhao_transient_outer.toml)はstrong-UV queue closureの
frozen-field guard fixtureです。記載した物理timescaleでは長いouter flightをfail-closedで停止させることが期待値であり、
成功RUNではありません。return-currentを解釈するには、別途根拠を与えた遅いfield timescaleを用い、`batch_duration`・粒子数・
control-volume長の収束確認が必要です。

この初版ではz-high interfaceをZhaoの有効放出面として扱います。`photoelectron_source_scale>0`では
`sim.sheath_alpha_deg`と`sim.sheath_photoelectron_ref_density_cm3`を再利用し、最初の負電荷
`photo_raycast` speciesから質量と温度を得ます。
Zhao solverはこの温度$T_{pe}$を電位scale、温度と基準密度から導出したphotoelectron Debye長
$\lambda_{D,pe}$を長さscaleとして使い、収束stateの`outer_debye_length_m`へ導出値を出力します。

`outer_plasma.debye_length`と`outer_plasma.thermal_voltage`はZhao rootやprofileの物理scaleには使いません。
ただし現時点ではsplit-interface適用性診断のreference inputとして残り、前者は`interface_eta_gap`と
local-charge推定、後者は横方向の`interface_eta_phi_kneq0`と`interface_eta_field_kneq0`のscaleに使われます。
Zhaoの収束はこれらを変えて判定せず、profile grid、有効interface位置、tracked粒子数、`dt`やbatch解像度を変えて確認します。
現行のprofile gridはruntime固定の128点であり、productionのgrid収束調査には点数の入力化が必要です。

tracked sourceの`photo_raycast.emit_current_density_a_m2`は、有効平面での解析raw source

$$
J_{pe,\mathrm{raw}}=s_{UV}
\frac{|q_{pe}|n_{\mathrm{ref}}\sin(\alpha)v_{\mathrm{th},pe}}{2\sqrt{\pi}},
\qquad
v_{\mathrm{th},pe}=\sqrt{\frac{2T_{pe}}{m_{pe}}}
$$

と1%以内で一致させます。ここで$s_{UV}$は`photoelectron_source_scale`です。この速度式の$T_{pe}$はJへ換算した
熱エネルギーです。一致しない設定はruntimeで拒否します。
analytic raw currentはtracked sourceの整合性検査とcurrent-density診断に使いますが、root、surface charge、ledgerへ
別途加えません。表面電荷とledgerはtracked放出・再吸収だけで更新します。$\eta$もcurrent診断のraw photoelectron
emission-current項をscaleしません。
legacy Zhao `sheath_injection_model`、`reservoir_potential_model`、
`photoelectron_density_model="kinetic_mean"`との併用も拒否します。

この有効平面近似は、tracked rayの方向分布やrough surfaceからinterfaceへ到達したVDFをZhao outer populationへ
自己無撞着に接続しません。`ray_direction`は照射rayによる放出面sampling、$\alpha$は解析Zhao sourceを決める独立の入力です。
適用範囲と一般化に必要な境界VDFについては[ADR 0003](adr/0003-zhao-charge-driven-outer-closure.md)、
過渡queueの決定は[ADR 0004](adr/0004-zhao-transient-photoelectron-column-queue.md)、continuation診断と
dynamic outerへの段階移行は[ADR 0005](adr/0005-zhao-continuation-and-dynamic-outer.md)に記録しています。

## `absorbing_maxwellian`の有限gridをRobin tailで無限遠へ接続する

この節と次のNewton法は既定の`absorbing_maxwellian`に固有です。Zhaoは前節のSagdeev条件からrootとprofileを構成します。
`absorbing_maxwellian`の内部点は伸長可能なnonuniform 1D grid上のconservative finite-volume residualです。現行runtime値は次のとおりで、
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

## continuation付きNewton法で物理解を追う（`absorbing_maxwellian`）

解析density derivativeからbordered-tridiagonal Jacobianを作るため、1 Newton stepは$O(N_z)$です。
$\phi_I=\phi_1$への内部densityの依存が、通常のtridiagonal stencilへborder列を加えます。

1. Newton stepを計算する。
2. backtrackingでtrial stateをsupported monotone branch内に制限する。
3. 通常Newtonが停滞すればpseudo-transient regularizationを使う。
4. 前batch profileがある場合、以前のinterface fieldから現在値までcontinuationし、失敗した刻みを半分にする。

regularizationとcontinuationは収束経路だけを変えます。最終判定は必ず元の離散Poisson残差へ戻ります。

## closureごとのsupported sheath branchだけを受理する

| 条件 | 必要な内容 |
| --- | --- |
| original Poisson residual | regularized式ではなく元の残差がtolerance以下 |
| branch | `absorbing_maxwellian`はsupported monotone branch、Zhaoは要求したA/B/Cの符号とpopulation条件 |
| ion accessibility | 全grid点で$u_i^2(z)>0$ |
| kinetic Bohm entry | $u_{i,\infty}\ge\sqrt{(T_e+\gamma_iT_i)/m_i}$ |
| infinity quasineutrality | $q_en_{e,\infty}+q_in_{i,\infty}\simeq0$ |

`absorbing_maxwellian`ではvirtual cathodeを持つ非単調profileとtrapped populationは適用外です。
`zhao_charge_driven`はType Aの単一potential minimumと、Zhao式に含まれるreflected/captured populationだけを扱います。
sub-Bohm ion inflowや各closureの条件を外れた解は受理せず、`not_applicable`、`no_physical_solution`、
`numerical_failure`のいずれかのstatusを付けて停止します。

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

`absorbing_maxwellian`では`debye_length`、interface位置、必要ならsource samplingを変え、interface電位、電流、
表面帯電が収束することを確認します。Zhaoでは導出$\lambda_{D,pe}$を使うため、profile grid、有効interface位置、
tracked粒子数、時間解像度を収束軸にします。returnを有効にした場合は
[粒子のescapeとreturn](ParticleEscapeReturn.html)のflight time、frozen-field ratio、準定常性も検査します。

## Code reference

- VDFに基づく密度モデルとnonlinear Poisson solve: [`bem_outer_plasma_kinetic.f90`](../src/physics/outer_plasma/bem_outer_plasma_kinetic.f90)
- charge-driven Zhao rootと非単調profile: [`bem_outer_plasma_zhao.f90`](../src/physics/outer_plasma/bem_outer_plasma_zhao.f90)
- runtime speciesからsolver optionを構成: [`bem_outer_plasma_kinetic_runtime.f90`](../src/runtime/bem_outer_plasma_kinetic_runtime.f90)
- surface fieldとの接続とMPI collective solve: [`bem_electrostatic_snapshot.f90`](../src/physics/bem_electrostatic_snapshot.f90)
- profile output: [`bem_output_writer.f90`](../src/runtime/bem_output_writer.f90)
- profile restart: [`bem_restart.f90`](../src/runtime/bem_restart.f90)

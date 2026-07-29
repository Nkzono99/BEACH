title: kinetic 1D外部プラズマ

Lang: [日本語](KineticOuterPlasma.md) | [English](KineticOuterPlasma.en.md)

# kinetic 1D外部プラズマ

`external_boundary.field.model="kinetic_1d"`は、periodic表面の上にある外部プラズマを、平面平均された1Dの静電・無衝突・
非磁化profileとして解きます。表面電荷が与えるinterface電場と無限遠の粒子VDFを接続し、interface電位、密度、
電流を自己整合に求めます。

BEACHでは、外部reservoirと平均シースを結ぶ標準モデルとして`kinetic_1d`を推奨します。
`external_boundary.particles.mode="same_batch"`を選ぶと、同じprofileを粒子の流入と通常のescape/returnにも使います。
後述の `implicit_mean` 光電子では、最初の local trace を interface crossing で止め、更新後の平均場が決まるまで
crossing を保留します。`ambient_linear_debye` は解析 Maxwellian の backward-Euler 平均更新、
`zhao_charge_driven` は実測 interface energy CDF の一般化電荷根を使う兄弟経路です。平均更新後はどちらも
kinetic 1D profile mapper を使いますが、平均 escape / return 量の決め方と ray weight は共有しません。
外部シースが必要なケースは、rough surface近傍の横方向screeningを解く明確な要件がない限り、このモデルから構成します。
その要件があるケースは現行モデルでは未対応です。削除した旧近似と再設計条件は
[ADR 0010](adr/0010-remove-unified-linear-response.md)に記録しています。

収束した外部profileは通常の粒子移送で双方向に使われます。無限遠VDFからinterfaceへ入る粒子の写像は
[`reservoir_face` の流入量と速度サンプリング](ReservoirInjection.html)、interfaceから出る粒子のescape/return写像は
[粒子のescapeとreturn](ParticleEscapeReturn.html)に対応します。本ページでは、`implicit_mean` が
この共通写像を 1 回の軌道 sample として使い、速い平均帯電と遅い局所再分配を分ける方法を説明します。

## local領域とouter領域をinterfaceで分ける

local粒子領域とouter領域はz-highのownership interface $z=z_I$で接続されます。

| 領域 | 使用する場 |
| --- | --- |
| meshからinterfaceまで | periodic $k\ne0$ surface fieldとsurface $k=0$ field |
| interfaceより外側 | `kinetic_1d`の平面平均profile |

outer profileをlocal領域の全点へ重ねません。このsplit構成では、interfaceまでに横方向modeが十分減衰すると仮定し、
interfaceで電位と法線電場を接続します。

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

enabledな負電荷・正電荷z-high `reservoir_face` speciesをそれぞれちょうど1つ要求し、無限遠ambient electronと
ionとして使います。複数候補から先頭を暗黙に選びません。
`photoelectron_density_model="kinetic_mean"`でも、enabledな負電荷`photo_raycast` speciesをちょうど1つ要求し、
その温度と放出電流密度を平面平均photoelectron sourceとして加えます。

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

## ambient 線形 Debye 応答と tracked 光電子を分離する

`kinetic_closure="ambient_linear_debye"` は、現在の surface zero mode が与える interface 電場 $E_I$ から

$$
\phi_I=\lambda_D E_I,\qquad
\phi(z)=\phi_I\exp(-z/\lambda_D),\qquad
\rho_{\mathrm{amb}}(z)=-\frac{\epsilon_0}{\lambda_D^2}\phi(z)
$$

を解析的に構成します。空間 Poisson root や Newton 反復は行いません。無限遠ゲージは
$\phi(\infty)=0$ であり、この同じ profile を ambient reservoir の到達判定・速度写像と、
通常の z-high return / escape 判定、および mean solve 後の `implicit_mean` 光電子判定に使います。
診断用 ambient 流入電流は、各 `reservoir_face` species の
drifting Maxwellian を同じ interface barrier で切った片側流束から計算します。

この closure では公開 TOML に `photoelectron_density_model` を書きません。facade が内部状態を `none` に解決し、
`"none"` を含む key の明記自体を拒否します。enabled な `photo_raycast` source は拒否せず、光電子を 3D 領域内で
放出し、再吸収または z-high interface の外向き crossing まで追跡します。
光電子の平均密度と outer space charge は closure に含めません。

`external_boundary.particles.mode="same_batch"` と enabled な負電荷 `photo_raycast` species もある場合、
BEACH は内部の `coupling.update_mode="implicit_mean"` を自動選択します。公開 TOML に update mode は書きません。
この multirate 更新では、連続 Maxwellian closure で同じ batch 内に決める速い平均帯電と、
有限個の tracked ray から次の batch へ渡す遅い局所電荷再分配を分けます。

| 成分 | 更新方法 |
| --- | --- |
| 局所 $k\ne0$ | batch-start の非零 mode operator を固定し、最初の 3D trace と帰還後の local trace で共有 |
| 速い平均 mode | backward-Euler mean solver が総電荷、$\phi_I$、連続 Maxwellian escape / return 率を決定 |
| energy-resolved ray | 解いた mean profile 上で各 crossing を 1 回だけ追跡し、軌道と終端着地点を sample |
| 遅い局所再分配 | 解析 return 総電荷を再吸収 sample の source / destination 分布へ正規化して 1 回だけ commit |

解析 escape / return へ置換する対象は、z-high を外向きに crossing して deferred された成分だけです。
interface 到達前の局所再吸収、z-low など別の open 面からの escape、unresolved discard は最初の trace の
tracked 値を保持します。以下の mean solve と return 正規化は、それらを上書きしません。

光電子の最初の local trace は z-high interface を横切った時点で終了し、full macro weight、放出元 element、
crossing の位置・速度を mean 更新へ渡します。その batch の crossing charge から
$J_{pe,\mathrm{out}}^n$ を測り、最初の trace で stage された全 surface charge delta から
$J_{\mathrm{pending}}^n$ を測り、

$$
J_{\mathrm{other}}^n
=J_{\mathrm{pending}}^n-J_{e,\mathrm{tracked}}^n-J_{i,\mathrm{tracked}}^n-J_{pe,\mathrm{out}}^n
$$

と定義して、平面平均電位を

$$
C_A\frac{\phi_I^{n+1}-\phi_I^n}{\Delta t}
=J_{e,\mathrm{tracked}}^n+J_{i,\mathrm{tracked}}^n
+J_{\mathrm{other}}^n
+J_{pe,\mathrm{out}}^n f_{\mathrm{esc}}(\phi_I^{n+1})+J_{\mathrm{ext}}
$$

を満たすように進めます。`e_bottom_zero` では $C_A=\epsilon_0/\lambda_D$、
`symmetric_vacuum` では $C_A=2\epsilon_0/\lambda_D$ です。ambient 電流は、その batch で
実際に表面へ吸収された一意な z-high `reservoir_face` 電子・イオン種から測定します。それ以外の
追加 species、他面流入、別 open 面 escape、未解決粒子に伴う実測 surface change は
$J_{\mathrm{other}}^n$ として transaction に残ります。`emit_current_density_a_m2` は tracked 3D 放出の
粒子 weight を決めますが、その値を独立な PE 平均 source として再利用しません。平均 source の振幅には実測
$J_{pe,\mathrm{out}}^n$ を使います。設定温度が決める Maxwellian barrier 通過率 $f_{\mathrm{esc}}$ を
$\phi_I^{n+1}$ で評価するため、平均的な引戻し電位を同じ step 内で反映します。
この backward-Euler 解を、batch 末の平均総電荷と $\phi_I$ の一意な基準にします。解析的な return / escape 電荷は

$$
Q_{\mathrm{ret}}^{\mathrm{ana}}
=A\Delta t\,J_{pe,\mathrm{out}}^n(1-f_{\mathrm{esc}}),
\qquad
Q_{\mathrm{esc}}^{\mathrm{ana}}
=A\Delta t\,J_{pe,\mathrm{out}}^n f_{\mathrm{esc}}
$$

です。どちらも非負の電荷絶対量として表しています。有限 ray の離散分類でこの総量を置き換えません。
ray 追跡用の element 電荷では、解析 return 率を各 crossing の source countercharge に掛け、その分を
放出元で一時的に中和します。これは mean solver が決めた総電荷を $k=0$ snapshot へ写すための一時分布であり、
物理的な再吸収先ではありません。

各 crossing は full macro weight のまま、mean solve 後の profile 上の法線 energy

$$
v_n^2(z)
=v_{n,I}^2+\frac{2q[\phi_I-\phi(z)]}{m}
$$

を離散 profile と far tail に沿って検査します。$v_n^2$ が途中で 0 になれば turning return、
profile を通過できれば escape です。return では同じ mapper が outer flight time を積分し、
$\mathbf v_t\tau_{\mathrm{outer}}$ の横方向変位を周期 wrap して interface 直下へ戻します。
その packet を batch-start の $k\ne0$ と mean solve 後の $k=0$ の合成場で local 3D 領域へ再追跡します。
return 後に z-high を再通過した場合も同じ driver と profile mapperへ戻し、最終的に surface へ吸収されるか
無限遠へ escape するまで処理します。
return 後に z-high 以外の local open face へ出た shadow は、解析的な上向き escape と同一視せず fail closed で停止します。

この deferred packet は解析 Maxwellian closure が決めた return 総量の軌道・着地点を標本化する準定常 shadow です。
mapper は各 shadow の outer flight time と frozen-field ratio を計算して診断へ残しますが、ratio が
`max_frozen_field_ratio` を超えても停止しません。通常の `same_batch` 粒子と ambient species は shadow ではなく、
同じ上限を超えると従来どおり fail closed で停止します。

`summary.txt` はこの shadow について、今回の実行で最後に完了した batch の analytic weight 適用後の
return excursion だけを 2 つの値へ集約します。restart 後に 1 batch も進めない no-op 実行では、復元すべき
solver state ではなく run-local な派生診断なので、この 2 key は出力しません。
return excursion $j$ の正の電荷 magnitude を $W_j>0$、outer 往復時間を $\tau_j$、
水平面積を $A$、batch 幅を $\Delta t$ とすると、

$$
\bar{\tau}_{\mathrm{ret},Q}
=\frac{\sum_j W_j\tau_j}{\sum_j W_j},
\qquad
\widehat{\sigma}_{\mathrm{PE,ret}}
=\frac{\sum_j W_j\tau_j}{A\Delta t}
=J_{\mathrm{ret,gross}}^{(|Q|)}\bar{\tau}_{\mathrm{ret},Q}.
$$

`implicit_mean_last_returned_outer_flight_time_mean_s` は電荷 magnitude で重み付けした
$\bar{\tau}_{\mathrm{ret},Q}$、`implicit_mean_last_estimated_returning_photoelectron_column_charge_per_area_C_m2`
は Little の法則による $\widehat{\sigma}_{\mathrm{PE,ret}}$ の準定常 shadow 推定値です。どちらも非負で、
return excursion がなければ 0 です。後者は実在する queue や ledger の cloud stock ではなく、escape と分類された
outcome 自体には滞在時間を足しません。ただし、その ray が最終的に escape する前に完了した return excursion は
集計します。また、全 batch の累積値や最大値ではありません。

`outer_integrated_charge_per_area_C_m2` が非零なら、

$$
\chi_{\mathrm{PE,shadow}}
=\frac{\widehat{\sigma}_{\mathrm{PE,ret}}}
       {\left|\texttt{outer\_integrated\_charge\_per\_area\_C\_m2}\right|}
$$

は、省略した returning-photoelectron shadow column の magnitude を、1D outer profile の積分電荷 magnitude と
比較する尺度です。これはモデル制限の解釈用であり、BEACH が課す受理閾値ではありません。

軌道追跡は mean solve 後に 1 回だけ行い、その分類を mean solve へ戻しません。最終的に surface へ再吸収された
sample の集合を $\mathcal R$ とし、

$$
W_{\mathcal R}=\sum_{j\in\mathcal R}|q_j|w_j,
\qquad
s_{\mathrm{ret}}=\frac{Q_{\mathrm{ret}}^{\mathrm{ana}}}{W_{\mathcal R}}
$$

で共通 scale を求めます。各 sample の source leg、すなわち放出元 countercharge の一時中和量と、
destination leg である実 hit deposit の両方へ $s_{\mathrm{ret}}$ を掛けます。このため transaction は零和のまま、
return 総電荷が解析値と一致します。最終 state は pending に含まれる放出 countercharge を保持し、
正規化した destination deposit を加えます。
sampled escape fraction は軌道統計の診断量であり、平均総電荷を決めません。

この 1-pass 構成では、平均総電荷と引戻し電位が同じ batch 内で速く応答し、sampled reabsorption pattern による
局所電荷分布は commit 後の batch から $k\ne0$ へ反映されます。離散分類を同じ batch の mean solve へ戻さないため、
separatrix 近傍の少数 ray が return / escape を切り替える 2-cycle を作りません。

この経路は `outer_update_stride=1`、正の `batch_duration`、ちょうど 1 つの負電荷
`photo_raycast` species、`deposit_opposite_charge_on_emit=true`、および
`photo_raycast.normal_drift_speed=0` を要求します。解析 Maxwellian escape 率は放出法線 drift を含まないため、
非零の値は fail closed で拒否します。
`photoelectron_density_model` は公開 TOML では省略し、内部で `none` に解決されます。update mode や
return kernel の公開 key は追加せず、outer queue とは併用できません。この経路は mesh mode に固有の制約を追加しません。
mean solver が収束しない、解析 return 電荷が正なのに再吸収 sample がない、transaction の電荷が釣り合わない、
または許可された local trace 内に吸収・escape の終端へ到達しない場合は fail closed で run を停止します。

適用範囲は ambient plasma の線形応答が支配的な弱光電子放出です。光電子 space charge による
virtual cathode、space-charge-limited / inverse sheath、trapped population、非単調 profile は扱えません。
光電子の nonlinear mean density や outer cloud inventory もこの経路では時間発展させません。特に、未帯電状態からの
UV turn-on と遅延 return current を解く transient model ではありません。その用途には
`ambient_linear_debye` の時間発展へ対応する delayed inventory / queue を別途設計する必要があります。
`summary.txt` の `coupling_update_mode=implicit_mean` と、標準出力の
`BEACH implicit-mean` 行にある $\phi_I$、species 電流、`J_other_A_m2`、
`transaction_residual_C`、`mean_solver_iterations`、`sample_escape_fraction`、`return_weight_scale` を確認します。
charge ledger の `escaped_to_infinity_C` はtrackedな別open面 escapeと解析的なz-high deferred escapeの和、
`absorbed_on_surface_C` はtrackedなlocal absorptionと解析的なz-high return後absorptionの和です。
`discarded_unresolved_C` はtracked値のままです。対応する整数 count と
`sample_escape_fraction` は ray sample の最終分類を表します。
`interface_outward_gross_C` は初回 crossing とその後の recross、
`interface_returned_gross_C` は正規化 weight で local 領域へ戻した return event の電荷を記録します。
z-high deferred成分の符号付き解析escape電荷を
$Q_{\mathrm{esc,z-high}}^{\mathrm{analytic}}$とすると、両者は
`interface_returned_gross_C` = `interface_outward_gross_C` -
$Q_{\mathrm{esc,z-high}}^{\mathrm{analytic}}$を満たします。別open面からのtracked escapeも含む
`escaped_to_infinity_C`総量を、この式へ代入しません。

ambient electron/ion reservoir を持たず、`external_boundary.field.model="none"` で z-high を escape させる
UV-only ケースは、この closure の代用ではありません。局所放出・再吸収を調べる有限 box 過渡 control として扱い、
無限遠で閉じた準中性シースまたは定常解とは解釈しません。

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

charge-driven Zhaoでは、現在のinterface電場を指定するため、定常Zhaoのzero-current式をrootへ課しません。
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

UVなしは`external_boundary.field.photoelectron_source_scale=0`で選びます。このときZhaoの光電子密度・raw emission currentは
厳密に0で、正規化密度と温度にはambient $n_\infty,T_e$を使います。設定されたambient electron densityは
準中性領域の総密度としてion densityとの準中性を検査し、solverがfree/reflected VDFに必要な入射electron densityを
解きます。reservoir macro数と速度cutoffにもその導出密度と同じprofile写像を使います。平坦な$E_I=0$はType B/Cの
接続点で、負の$E_I$ではType Cへ入ります。

full photoelectron populationの準定常A/B/C branchは、必ずしも$E_I=0$へ連続しません。要求電場がbranchの
可解域外なら、`external_boundary.particles.mode="same_batch"`では`no_physical_solution`で停止します。

### 強 PE を実測 interface energy CDF で閉じる

photoemitting `zhao_charge_driven`を`same_batch`と組み合わせると、BEACHは内部の`implicit_mean`を
自動選択します。この経路は、前節のambient backward-Euler平均更新とは別のsolverです。
最初からcommit済みA/B/C branchを必要とするため、
`external_boundary.particles.steady_start_mode="zhao_floating"`と、その一様電荷を置く
`steady_start_mesh_id`も必須です。

最初の3D traceではbatch-startの場を固定し、光電子がz-highを最初に外向き通過したときの法線energyと
macro charge magnitudeを採取します。

$$
{\cal E}_{n,r}=\frac{1}{2}m_{pe}v_{z,r}^2,\qquad w_r=|\Delta q_r|
$$

各rankの$({\cal E}_{n,r},w_r)$はMPI rootへ集めます。rootはenergyを降順にstable sortし、
浮動小数点比較で完全に等しいenergyだけをgroup化してcharge weightを加算します。設定bin、histogram補間、
smoothingは使わないため、macro粒子ごとのweightを保ったexact empirical CDFです。

最初のtraceでstageされた全surface chargeから、interfaceへ出た全光電子source chargeを除いた値を
$Q_{\rm base}$とします。候補となるcell総電荷$Q$にはlower boundaryに応じて

$$
Q=\beta\epsilon_0AE_I,\qquad
\beta=
\begin{cases}
1,&\texttt{e\_bottom\_zero},\\
2,&\texttt{symmetric\_vacuum}
\end{cases}
$$

を対応させ、同じZhao branch上で1D profileを解きます。光電子の脱出に必要なenergyは、interfaceと無限遠の
端点差だけではなく、profile全体から

$$
B(Q)=\max\left[
0,\ \max_z q_{pe}\{\phi(z;Q)-\phi_I(Q)\},\
q_{pe}\{\phi_\infty-\phi_I(Q)\}
\right]
$$

と求めます。$q_{pe}<0$なので、Type Aのpotential minimumが作るvirtual-cathode barrierもこの最大値へ入ります。

実測CDFから、${\cal E}_{n,r}>B$のgroupをescape、equalityをturning / returnとして
$C_{\rm esc}(B)$を構成します。解く式は

$$
Q=Q_{\rm base}+C_{\rm esc}[B(Q)].
$$

barrierが1つのenergy groupを横切る場合だけ、そのgroupのescape weight
$0\le\theta\le1$を未知量として$B(Q)={\cal E}_k$を満たすfractional marginal rootを使います。
これにより、同じenergyを持つmacro粒子groupを恣意的に別groupへ分けません。

group数を$M$、$Q_k=Q_{\rm base}+C_k$とします。solverは最終source上の共通rootから
$Q_0$と$Q_M$へconnected pathを追跡し、全候補charge区間でbarrierの数値的な非減少性を検査します。
その上で

$$
P_k=
\begin{cases}
[B(Q_k)\ge{\cal E}_{k+1}], & 0\le k<M,\\
[B(Q_M)\ge{\cal E}_M], & k=M
\end{cases}
$$

のfirst-true indexを二分探索します。pure root選択のconnected candidate solve数は$O(\log M)$です。
$P_0$がtrueならall-return、$P_M$がfalseならall-escapeです。内部のfirst-true index $k$では、
$B(Q_k)<{\cal E}_k$ならpure root、そうでなければ$[Q_{k-1},Q_k]$でfractional marginal rootを解きます。
marginal二分はcharge幅まで続け、turning-point equalityをreturn側に置くためupper endpointを採用します。
各候補のZhao rootは独立なNewton solveで選び直さず、共通rootから連結parameter pathを追跡して求めます。

$$
\mathbf G\!\left(\mathbf y;E_I(\lambda),n_{pe,0}(\lambda)\right)=\mathbf 0,
\qquad
E_I(\lambda)=(1-\lambda)E_{I,0}+\lambda E_{I,1},
\qquad
n_{pe,0}(\lambda)=(1-\lambda)n_{pe,0}^{(0)}+\lambda n_{pe,0}^{(1)}.
$$

$\mathbf y$は固定したA/B/C branchのlog-encoded root変数です。Type A/Bは$E_I>0$、Type Cは$E_I<0$を要求し、
このchartでregularでない$E_I=0$も拒否します。adaptive pseudo-arclength predictor/correctorは、局所補正距離、
root jump、接線方向、Jacobian rankを検査しながら前のrootから進みます。$\lambda$方向の接線が消失・反転した場合や
固定parameter Jacobianのrankが失われた場合はtarget前のfoldとしてfail closedし、targetは$\lambda=1$の
event correctorで着地します。これは近傍rootへの数値的jumpを制限するlocal continuation guardであり、
任意に近い別sheetが存在しないことの数理証明ではありません。この内部経路に新しいTOML keyはありません。

最終sourceを固定したcharge候補では、各adaptive accepted pointで$B$のprescribed chargeに対する接線勾配を評価し、
負の勾配または逆向きのbarrier secantを検出すると停止します。endpoint、order-statisticのmidpoint、
marginal二分でもbarrierのcharge順を検査します。これはaccepted point間の連続区間をinterval arithmeticで
包含する数理的証明ではありません。最後に実測CDFからescape chargeを再計算し、
$Q-Q_{\rm base}-C_{\rm esc}=0$を要求します。rank低下、fold、barrier勾配反転、order predicateのbracket不整合、
またはmarginal energyをbracketできない場合は、別branchへfallbackせず停止します。

surface emissionとinterface source normalizationは別の量です。
`photo_raycast.emit_current_density_a_m2`はrough surfaceから放出するray数とmacro weightを決めます。
Zhao density closureへ渡す各batchのsource amplitudeは、実際にinterfaceへ到達したchargeから

$$
J_{pe,I}^{\rm meas}=\frac{\sum_r w_r}{A\Delta t},\qquad
s_{\rm eff}=
\frac{J_{pe,I}^{\rm meas}}
{|q_{pe}|n_{\rm ref}\sin(\alpha)\sqrt{2T_{pe}/m_{pe}}/(2\sqrt{\pi})}
$$

と解きます。設定した`photoelectron_source_scale`は`zhao_floating`の初期branch anchorを構成しますが、
同じbatchのinterface sourceとして再利用しません。受理した$s_{\rm eff}$はresolved source scaleとして
outer stateとcheckpointへ保存されます。

実測sourceが前batch値と異なる場合は、上記の連結pathで前batch rootから最終$s_{\rm eff}$上の共通anchorへ進みます。
その後の全charge候補はsourceを最終値に固定し、共通anchorからcommit済みA/B/C root sheet上を追います。
branch labelの変更、bootstrap branch、disconnected root、別closureへのfallbackは許しません。
各continuation stepでは$\phi_0$、$\phi_{\min}$、無限遠electron densityの無次元変化を0.25以下に制限します。

同じbatch-start場から得たfrozen cohortを更新後のprofileへ使うため、最終状態にも

$$
\frac{|q_e\Delta\phi_I|}{T_e}\le0.25,\quad
\frac{|\Delta B_{e,\mathrm{in}}|}{T_e}\le0.25,\quad
\frac{|q_i\Delta\phi_I|}{E_{i,n}}\le0.25,\quad
\left|\log\frac{n_{e,\infty}^{\mathrm{new}}}{n_{e,\infty}^{\mathrm{old}}}\right|\le0.25,\quad
|\log(s_{\rm eff}/s_{\rm old})|\le0.25,\quad
\frac{\lambda_{D,pe}|\Delta E_I|}{T_{pe}/|q_{pe}|}\le0.25,\quad
\frac{|\Delta B|}{T_{pe}}\le0.25
$$

を要求します。超過時にrayを取り直したり旧profileへ戻したりせずfail closedで停止します。
ここで$B_{e,\mathrm{in}}=\max(0,q_e\phi_I,q_e\phi_{\min})$、
$E_{i,n}=\max(T_i,m_i u_{i,n}^2/2)$です。最初の4項は、無限遠からinterfaceへのambient reservoir写像で
凍結される絶対電位差、Type Aを含む電子cutoff、ion drift energy、Zhaoが解く電子流入密度をそれぞれ監視します。
一方、光電子はprofileへ共通の電位移動ではなく、profile内のbarrier $B$とfield変形へ応答するため、
最後の2項は$T_{pe}$と$\lambda_{D,pe}$で測ります。同じ$\Delta\phi_I$を$T_{pe}$でも重複判定しません。

このCDFはbatchごとの生の実測標本なので、収束判定では総ray数だけでなく、interface到達数と
barrierより上のescape tailに残る有効標本数を確認します。source、field、barrierの変化が小さいまま
ambient-potential trust guardだけが統計的に失敗する場合も、guardを緩めず`rays_per_batch`を増やし、
少なくとも2倍のray数で$\phi_I$、escape率、表面総電荷が一致することを確認します。
ray数を増やしても同じguard超過が収束して残る場合は、統計誤差ではなく1 batchの物理更新幅が大きすぎます。
guardは緩めず、`batch_duration`と対応するbatch当たり放出量を小さくして時間刻み収束を確認します。

受理後は各rayのenergy groupが決めたescape / return weightを、放出元elementと実際の再吸収elementへ
そのまま使います。ambient経路の共通`return_weight_scale`は適用しません。各rayは原則1回のoutward crossingと
最大1回のreturnを持つ必要があり、return weightを持つrayが吸収されない場合、escape rayのterminal outcomeが
一致しない場合、または再越境した場合は、他rayへweightを付け替えず停止します。

実測するのはinterface source amplitudeと法線energy CDFです。Zhaoのfree/reflected ambient electron、
free/captured photoelectron、cold ionのdensity closureは解析式のままです。接線速度分布や任意のtrapped VDF、
衝突、磁化、時間依存Vlasov--Poisson、PICを解くmodelではありません。Type Aで扱うのも単一の
virtual-cathode minimumとZhao式に含まれるreflected / captured populationだけです。

現行MPI実装は全interface sampleをrootへ`MPI_Gatherv`してからgroup化とZhao solveを行い、受理したstepと
profileをbroadcastします。rootの一時memoryとgather通信量はinterfaceを横切る光電子数に比例します。
ray数を増やす収束調査では、このroot-local制約も見積もってください。

### 定常研究を Zhao 零電流根から開始する

強UVの定常観測量が目的で、未帯電状態からの立上がり過渡を研究対象にしない場合は、
`external_boundary.particles.steady_start_mode="zhao_floating"`を選べます。このmodeは最初batch前に、設定済みの無限遠準中性条件、
温度、drift、UV sourceを使って Zhao 零電流定常根を解きます。得られた根から
`phi(infinity)=0`のprofileを作り、そのinterface電場$E_I$に必要な一様平面電荷を初期値にします。

水平面積を$A$とすると、電荷はzero-modeの下側境界条件に応じて

$$
Q_{seed}=2\epsilon_0AE_I
\quad\text{(\texttt{symmetric\_vacuum})},
$$

または

$$
Q_{seed}=\epsilon_0AE_I
\quad\text{(\texttt{e\_bottom\_zero})}
$$

です。`steady_start_mesh_id`で選んだ水平平面の三角形に、この電荷を面積比で配ります。選択平面は同一高さで
periodic cellの水平面積全体を覆い、outer interfaceより下にある必要があります。現在は非重複・無欠損tilingを
構築時に保証する`mesh.mode="template"`だけを受け入れ、任意OBJ平面は拒否します。他のmeshは電荷0のままなので、
plane + sphereでplaneをmesh 1とすると、planeだけをseedしsphereは中性で開始します。

定常根から構成した同一profileを、最初のouter state、無限遠reservoirからinterfaceへの密度・速度写像、
およびinterface外のinstant return / escapeに使います。零電流条件を課すのは初期状態の構成時だけで、
後続batchのcharge-driven更新は従来どおり現在の表面電荷が決める$E_I$を解き、currentを診断します。
analytic currentを表面電荷に加えないため、tracked放出・流入・再吸収と二重計上しません。

この初期化はphysical transientではありません。`external_boundary.particles.mode="same_batch"`のinstant経路だけで使い、
新規実行で既存電荷を上書きしません。同一configで`output.resume=true`とした場合はcheckpointのmesh電荷と完全なouter stateを
復元し、零電流根から再seedしません。queue過渡closureとの併用は拒否します。UV立上がり、return-current遅延、
cloud inventoryの過渡が問いなら、引き続きqueueまたは動的outer modelが必要です。定常論文でも、warm startだけを
一意性や動的安定性の根拠にせず、独立に緩和させた状態または摂動seedから同じ定常量へ返るかを確認します。

この用途とqueue過渡closureとの役割分担は
[ADR 0006](adr/0006-zhao-stationary-warm-start.md)に記録しています。

`external_boundary.particles.mode="zhao_queue"`では、追跡中の光電子がouter領域に滞在する間、
そのmacro粒子重みをqueue inventoryとして保持します。
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

[`periodic2_zhao_charge_driven_outer.toml`](../examples/periodic2_zhao_charge_driven_outer.toml)は、
`zhao_floating`の定常Type A根をanchorとして、small batchで上記の実測CDFによる強PE更新を通すsmokeです。
正常終了は実装経路の確認であり、ray数、batch幅、interface位置に対する物理収束を示しません。
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

`external_boundary.field.debye_length`と`external_boundary.field.thermal_voltage`はZhao rootやprofileの物理scaleには使いません。
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

と、固定sourceを使う`implicit_mean`以外のZhao経路では1%以内で一致させます。ここで$s_{UV}$は
`photoelectron_source_scale`です。この速度式の$T_{pe}$はJへ換算した熱エネルギーです。
強PE `implicit_mean`ではsurface currentとの1%一致を課さず、実測interface currentから$s_{\rm eff}$を解きます。
analytic raw currentはtracked sourceの整合性検査とcurrent-density診断に使いますが、root、surface charge、ledgerへ
別途加えません。表面電荷とledgerはtracked放出・再吸収だけで更新します。$\eta$もcurrent診断のraw photoelectron
emission-current項をscaleしません。
`external_boundary.particles.inflow_model="infinity_barrier"`、
`external_boundary.field.photoelectron_density_model="kinetic_mean"`との併用も拒否します。

固定source Zhaoの有効平面近似は、tracked rayの方向分布やrough surfaceからinterfaceへ到達したVDFを
Zhao outer populationへ自己無撞着に接続しません。強PE `implicit_mean`も、実測するのはsource amplitudeと
法線energy CDFだけで、Zhao density populationや接線VDFは解析modelのままです。`ray_direction`は照射rayによる
放出面sampling、$\alpha$は解析Zhao density closureを決める独立の入力です。
適用範囲と一般化に必要な境界VDFについては[ADR 0003](adr/0003-zhao-charge-driven-outer-closure.md)、
過渡queueの決定は[ADR 0004](adr/0004-zhao-transient-photoelectron-column-queue.md)、continuation診断と
dynamic outerへの段階移行は[ADR 0005](adr/0005-zhao-continuation-and-dynamic-outer.md)、定常warm startは
[ADR 0006](adr/0006-zhao-stationary-warm-start.md)に記録しています。

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
個別 profile return を使う species について、
[粒子のescapeとreturn](ParticleEscapeReturn.html)の flight time、frozen-field ratio、準定常性も検査します。
ambient `implicit_mean`では、これらに加えて`mean_solver_iterations`、transaction residual、
sample escape fraction、共通return weight scaleのray数依存性を検査します。強PE empirical Zhaoでは、
branch、$\phi_{\min}$、$n_{e,\infty}$、full-profile barrier、resolved source scale、marginal energy / fraction、
empirical charge residual、recross fraction、terminal mismatch fractionを検査します。
どちらの`implicit_mean`でもfrozen-field ratioはshadow軌道の時間scaleを示す診断であり、
通常のsame-batch粒子と異なって上限超過だけではrunを停止しません。

## Code reference

- VDFに基づく密度モデルとnonlinear Poisson solve: [`bem_outer_plasma_kinetic.f90`](../src/physics/outer_plasma/bem_outer_plasma_kinetic.f90)
- charge-driven Zhao rootと非単調profile: [`bem_outer_plasma_zhao.f90`](../src/physics/outer_plasma/bem_outer_plasma_zhao.f90)
- runtime speciesからsolver optionを構成: [`bem_outer_plasma_kinetic_runtime.f90`](../src/runtime/bem_outer_plasma_kinetic_runtime.f90)
- ambient mean の backward-Euler 更新: [`bem_mean_charge_integrator.f90`](../src/physics/outer_plasma/bem_mean_charge_integrator.f90)、[`bem_dynamic_k0_mean.f90`](../src/runtime/coupling/bem_dynamic_k0_mean.f90)
- 強 PE の実測 energy CDF と非線形 Zhao 更新: [`bem_dynamic_k0_zhao.f90`](../src/runtime/coupling/bem_dynamic_k0_zhao.f90)
- one-pass ray と正規化 transaction: [`bem_mean_charge_transaction.f90`](../src/runtime/coupling/bem_mean_charge_transaction.f90)、[`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- 定常根と平面電荷のwarm start: [`bem_zhao_steady_start.f90`](../src/runtime/coupling/bem_zhao_steady_start.f90)
- surface fieldとの接続とMPI collective solve: [`bem_electrostatic_snapshot.f90`](../src/physics/bem_electrostatic_snapshot.f90)
- profile output: [`bem_output_writer.f90`](../src/runtime/bem_output_writer.f90)
- profile restart: [`bem_restart.f90`](../src/runtime/bem_restart.f90)

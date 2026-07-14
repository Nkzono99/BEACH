title: 粒子のescapeとreturn

Lang: [日本語](ParticleEscapeReturn.md) | [English](ParticleEscapeReturn.en.md)

# 粒子のescapeとreturn

open境界へ到達した粒子は、単純に除去するか、scalar障壁で反射するか、外部プラズマ領域へ渡して
return/escapeを判定します。この選択は粒子sourceに依存しません。reservoir粒子、光電子、`volume_seed`粒子は、
同じ状態で同じ面を横切れば、同じ境界処理を受けます。

## 境界通過時の状態を外部modelへ渡す

[粒子の衝突・境界イベント](ParticleEvents.html)は、mesh hitとbox面の通過のうち、軌道上で最初に起きるものを選びます。
通常のopen面では、交差位置と速度を使ってescapeを判定します。`particle_transfer_mode`が有効なz-high面だけは、
交差情報をouter modelへ渡します。

交差情報には次が含まれます。

- interface上の位置と外向き速度。
- local Boris step内の交差時刻。
- 交差後に残る`dt_remaining`。

outer modelがreturnと判定した場合、粒子をinterfaceの直内側へ戻します。その後、`dt_remaining`の間だけ、
通常のBoris更新と衝突・境界処理をやり直します。outer flight timeは、このlocal step remainderとは別の診断量です。

## outer stateを持たないopen面では粒子を除去する

`open_boundary_model="escape"`かつouter transfer対象でないopen面では、粒子をその場で除去します。
macro電荷$q w$をspecies別`escaped_to_infinity`へ計上し、表面電荷`q_elem`は変更しません。

複数面を同時に横切るcornerでも、openを含む通常の境界規則が決定論的に適用されます。reflect/periodic面との
組合せと残りstepの再積分は[粒子の衝突・境界イベント](ParticleEvents.html)で説明します。

## scalar barrierでは通過点の法線energyだけを比較する

`open_boundary_model="potential_barrier"`は、open面の通過点電位$\phi_b$と`sim.phi_infty`だけを使うreduced modelです。
外向き法線速度を$v_n>0$とすると、無限遠へ進むためのpotential barrierは

$$
U_b=q(\phi_\infty-\phi_b)
$$

です。

$$
\frac12 m v_n^2<U_b\quad\text{かつ}\quad U_b>0
$$

なら法線速度を反転し、残りstepを追跡します。それ以外はescapeです。接線速度は変えません。

このmodelは、外部状態として通過点のscalar potentialだけを保持します。判定結果はreflectまたはescapeです。
open面の外側にある$E(\mathbf x)$、turning位置、flight time、空間電荷は扱いません。

複数のopen面を同時に横切るcornerには対応しておらず、`unsupported_barrier_corner`で停止します。
通過点電位は粒子運動と同じsnapshot規約で評価するため`sim.e0`の局所電位も含みます。一様電場には有限な
無限遠電位がないため、`sim.e0` と併用する場合の `phi_infty` は有効なreservoir基準電位としてユーザが
整合させる必要があります。

## z-highを粒子ownershipのinterfaceにする

`coupling.particle_transfer_mode`は、z-high open面を粒子ownershipのinterfaceとして外部領域へ接続します。

| `particle_transfer_mode` | 対応するfield/return構成 | 粒子処理 |
| --- | --- | --- |
| `none` | outer transferなし | 通常のopen境界 |
| `electrostatic_1d_instant_return` | `linear_debye`または`kinetic_1d` | 保存エネルギーからescape/returnを直接写像 |
| `electrostatic_3d_explicit_orbit` | `unified_linear_response` | batch内で固定された3D電場で外部軌道を時間積分 |

1D/3D transferはいずれも、z-highがopenであることと`sim.b0=0`を要求します。現行のouter particle modelは
外部領域での磁場軌道を扱いません。return位置のx/yはprimary periodic cellへwrapされます。

## 一つのenergy式で流入と流出を結ぶ

interface電位を$\phi_I$、無限遠を$\phi_\infty$、interfaceでの外向き法線速度を$v_{n,I}$とすると、

$$
v_{n,\infty}^2
=v_{n,I}^2+\frac{2q(\phi_I-\phi_\infty)}{m}
$$

です。$v_{n,\infty}^2\ge0$なら無限遠へ到達できるためescape、負ならouter profile内にturning pointがあるため
returnです。これは[reservoir注入](ReservoirInjection.html)で無限遠からinterfaceへ写す式の逆向きです。
両方向で同じ$\phi_I-\phi_\infty$を使うことが、流入と流出を一つの外部modelへ結び付けます。

## linear Debye profileからreturn時間を求める

`outer_plasma.model="linear_debye"`では、interfaceから外側の電位差をDebye長$\lambda_D$の指数profileで表します。
escape不能粒子について

$$
D=-v_{n,\infty}^2>0
$$

とすると、実装が使う往復時間は

$$
\tau_\mathrm{outer}
=\frac{4\lambda_D}{\sqrt D}
\tan^{-1}\left(\frac{v_{n,I}}{\sqrt D}\right)
$$

です。return時は法線速度だけを$-v_{n,I}$へ反転し、接線速度は保持します。接線位置を
$\mathbf v_t\tau_\mathrm{outer}$だけ進め、x/yをperiodic wrapします。

この解析的な状態写像は、保存エネルギーから外部滞在時間と接線変位を求めてreturn状態を作ります。

## 離散kinetic profile上でturning pointを探す

`outer_plasma.model="kinetic_1d"`と`return_model="kinetic_1d_profile_return"`は、batch開始時に収束した離散
$\phi(z)$を使います。各grid点で

$$
v_n^2(z)=v_{n,I}^2+\frac{2q[\phi_I-\phi(z)]}{m}
$$

を評価し、最初に符号が変わる区間でturning pointを線形補間します。正の速度を持つ区間$a\to b$の片道時間は

$$
\Delta t=\frac{2\Delta z}{\sqrt{v_{n,a}^2}+\sqrt{v_{n,b}^2}}
$$

で積算します。turning区間では到達fractionまでを積分し、grid上端より先にturning pointがある場合は
far Robin exponential tailを解析積分します。往復時間を得た後のreturn位置・速度はlinear Debyeの場合と同じです。

profileは有限値、単調、interface点との整合を検査します。非単調profile、物理的turning pointをbracketしない区間、
非正のRobin tailなどはinvalid modelとして停止します。

## unified 3D field中で外部軌道を進める

`unified_linear_response + electrostatic_3d_explicit_orbit`は、zero modeとscreened nonzero modeを合成した
batch内で固定された電場で外部粒子を追跡します。固定刻み`outer_orbit_dt`のvelocity-Verlet更新は

$$
\mathbf v^{n+1/2}=\mathbf v^n+\frac{q\mathbf E(\mathbf x^n)}{2m}\Delta t_o,
$$

$$
\mathbf x^{n+1}=\mathbf x^n+\mathbf v^{n+1/2}\Delta t_o,
\qquad
\mathbf v^{n+1}=\mathbf v^{n+1/2}+\frac{q\mathbf E(\mathbf x^{n+1})}{2m}\Delta t_o
$$

です。ownership interfaceを内向きに再通過すればreturn、unified grid上端のfar planeを外向きに通過すればescapeです。
境界到達位置と速度は最後のouter step内で線形補間します。

`outer_orbit_max_steps`までにどちらの境界にも到達しない場合は、persistent queueが必要になるため停止します。
また、初期状態とreturn/escape境界へ到達したときの全エネルギー

$$
\mathcal E=\frac12m|\mathbf v|^2+q\phi(\mathbf x)
$$

の相対誤差が`outer_orbit_energy_tolerance`を越えても停止します。`outer_orbit_dt`はreturn/escape率、flight time、
energy errorに対して収束確認します。

## outer flightをglobal timeへ加えない近似

1D returnの「instant」は、outer flightを状態写像には使う一方、global simulation timeを進めないという意味です。
3D explicit orbitも現行実装では同じくouter flightをglobal clockへ加えません。粒子はinterfaceを出たlocal stepと
同じsimulation時刻へ戻り、outward/returned chargeは同じbatchに計上されます。

この近似の対象は、定常または準定常outer plasmaです。エネルギーに基づくescape/returnの分類と長時間平均電流を
保ちつつ、外部粒子stockを状態から除去できます。UV照射の開始、plasma条件の急変、短pulseへの過渡応答を扱うには、
過去のoutgoing currentを保持するdelayed-return queueが必要です。

## outer flight中にfieldを固定できる範囲を判定する

outer flight中に場を固定できるかを

$$
\epsilon_\mathrm{ad}
=\frac{\tau_\mathrm{outer}}{\texttt{field\_evolution\_timescale}}
$$

で評価し、`max_frozen_field_ratio`以下を要求します。分母の`field_evolution_timescale`には、表面電位やouter profileが
変化する物理時間を与えます。

$\tau_\mathrm{outer}/\mathrm{batch\_duration}\gtrsim1$の場合、十分に定常化した長時間平均は使えますが、
batchごとのreturn currentは正しい時間変化を表しません。

条件を越えた粒子を後のbatchへ保持するpersistent delayed-return queueは未実装です。
`outer_queue_enabled=true`は設定検証で拒否され、queueが必要な軌道は停止します。

## interface通過量と最終escapeを分けて記録する

species別に`interface_outward_gross`、`interface_returned_gross`、`escaped_to_infinity`を区別します。さらに、最大
`outer_flight_time`、frozen-field ratio、3D orbitのenergy errorを出力します。gross outwardとreturnedは、どちらもinterfaceの通過量です。
その差が正味escapeと一致するのは、そのspeciesのtransfer対象と電荷収支の集計期間が一致する場合に限ります。
各項目と`charge_ledger.csv`の読み方は[出力の読み方](OutputGuide.html)にまとめています。

## Code reference

- box境界処理とscalar barrier: [`bem_particle_stepper.f90`](../src/runtime/simulator/bem_particle_stepper.f90)
- linear/kinetic 1D return: [`bem_outer_plasma_interface.f90`](../src/physics/outer_plasma/bem_outer_plasma_interface.f90)
- unified 3D orbit: [`bem_outer_plasma_orbit.f90`](../src/physics/outer_plasma/bem_outer_plasma_orbit.f90)
- interface transferとdiagnostic集計: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- return modelの組合せ検証: [`bem_physics_config_types.f90`](../src/config/bem_physics_config_types.f90)

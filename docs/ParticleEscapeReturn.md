title: 粒子のescapeとreturn

Lang: [日本語](ParticleEscapeReturn.md) | [English](ParticleEscapeReturn.en.md)

# 粒子のescapeとreturn

open境界へ到達した粒子の扱いは、次の5つのmodelに整理できます。

1. `escape`: 粒子をその場で除去する。
2. `potential_barrier`: 通過点のscalar電位から反射またはescapeを判定する。
3. `linear_debye`: 解析的な1D Debye profileでescape/returnを写像する。
4. `kinetic_1d`: 自己整合に解いた離散1D profileでescape/returnを写像する。
5. `unified_linear_response`: 外部3D場中の軌道を積分してescape/returnを判定する。

この一覧は境界処理アルゴリズムの比較であり、運用上の推奨順位ではありません。自己整合な外部シースの標準経路は
`kinetic_1d`と対応する1D transferです。`unified_linear_response`と3D explicit orbitは、rough surfaceの線形screeningと
3D外部軌道の両方が必要な場合の高度な経路です。

ただし、これらは同じ設定キーの5択ではありません。`escape`と`potential_barrier`は通常のopen面に使う
`sim.open_boundary_model`です。残り3つはz-highを外部領域へのownership interfaceにしたときの
`outer_plasma.model`であり、対応する`return_model`と`particle_transfer_mode`を組み合わせます。

```text
open面を粒子が通過
├─ outer transferの対象外
│  ├─ escape                    無条件に除去
│  └─ potential_barrier         scalar障壁で反射またはescape
└─ z-highをouter transfer
   ├─ linear_debye              解析的1D return
   ├─ kinetic_1d                離散1D profile return
   └─ unified_linear_response   明示的3D outer orbit
```

## 5つのmodelを比較する

| model | 外部状態 | 粒子の判定方法 | return時間 | 主な用途 |
| --- | --- | --- | --- | --- |
| `escape` | なし | 常にescape | なし | 単純な有限box |
| `potential_barrier` | 通過点のscalar電位 | 法線energyと障壁を比較 | なし | 低コストな局所反射 |
| `linear_debye` | 解析的な指数1D profile | 保存energy | 解析式 | 簡易な準定常1D outer plasma |
| `kinetic_1d` | 収束した離散1D profile | 保存energyとturning point探索 | profileを積分 | **標準:** 自己整合な平均sheath |
| `unified_linear_response` | zero/nonzero modeを含む3D場 | 外部軌道を時間積分 | 軌道から計測 | **高度:** rough surfaceを含む線形3D応答 |

どのmodelでも粒子sourceには依存しません。reservoir粒子、光電子、`volume_seed`粒子は、同じ状態で同じ面を
横切れば同じ境界処理を受けます。外部場自体の選び方は
[境界・外部領域の構成を選ぶ](OuterPlasmaModels.html)にまとめています。

## 1. `escape`: open面で粒子を除去する

最も単純なmodelです。outer transfer対象でないopen面に対して次を指定します。

```toml
[sim]
open_boundary_model = "escape"
```

粒子は境界通過位置で直ちに除去されます。macro電荷$qw$はspecies別`escaped_to_infinity`へ計上され、
表面電荷`q_elem`は変更されません。外部場、turning point、flight time、return状態は持ちません。

複数面を同時に横切るcornerでも、openを含む通常の境界規則が決定論的に適用されます。reflect/periodic面との
組合せと残りstepの再積分は[粒子の衝突・境界イベント](ParticleEvents.html)で説明します。

## 2. `potential_barrier`: scalar障壁で反射を判定する

通過点のscalar電位だけを使うreduced modelです。外部の空間profileを構築せずに、局所的なenergy条件で
reflectまたはescapeを選びます。

```toml
[sim]
open_boundary_model = "potential_barrier"
phi_infty = 0.0
```

### 通過点の法線energyだけを比較する

open面の通過点電位を$\phi_b$、外向き法線速度を$v_n>0$とすると、無限遠へ進むためのpotential barrierは

$$
U_b=q(\phi_\infty-\phi_b)
$$

です。

$$
\frac12 m v_n^2<U_b\quad\text{かつ}\quad U_b>0
$$

なら法線速度を反転し、残りstepを追跡します。それ以外はescapeです。接線速度は変えません。

### このmodelが表さないもの

このmodelが保持する外部状態は通過点のscalar potentialだけです。open面の外側にある$E(\mathbf x)$、
turning位置、flight time、空間電荷は扱いません。複数のopen面を同時に横切るcornerにも対応しておらず、
`unsupported_barrier_corner`で停止します。

通過点電位は粒子運動と同じsnapshot規約で評価するため`sim.e0`の局所電位も含みます。一様電場には有限な
無限遠電位がないため、`sim.e0`と併用する場合の`phi_infty`は有効なreservoir基準電位としてユーザが
整合させる必要があります。

## 3. `linear_debye`: 解析的な1D profileでreturnを写像する

interfaceから外側の電位差をDebye長$\lambda_D$の指数profileで表すmodelです。z-highをownership interfaceにし、
解析式からescape/return、往復時間、接線変位を一度に求めます。

```toml
[outer_plasma]
model = "linear_debye"
return_model = "electrostatic_1d_instant_return"

[coupling]
particle_transfer_mode = "electrostatic_1d_instant_return"
```

### 保存energyでescapeとreturnを分ける

interface電位を$\phi_I$、無限遠を$\phi_\infty$、interfaceでの外向き法線速度を$v_{n,I}$とすると、

$$
v_{n,\infty}^2
=v_{n,I}^2+\frac{2q(\phi_I-\phi_\infty)}{m}
$$

です。$v_{n,\infty}^2\ge0$なら無限遠へ到達できるためescape、負なら指数profile内にturning pointがあるため
returnです。これは[`reservoir_face` の流入量と速度サンプリング](ReservoirInjection.html)で、無限遠からinterfaceへ写す式の逆向きです。

### linear Debye profileからreturn時間を求める

escape不能粒子について

$$
D=-v_{n,\infty}^2>0
$$

とすると、往復時間は

$$
\tau_\mathrm{outer}
=\frac{4\lambda_D}{\sqrt D}
\tan^{-1}\left(\frac{v_{n,I}}{\sqrt D}\right)
$$

です。return時は法線速度だけを$-v_{n,I}$へ反転し、接線速度は保持します。接線位置を
$\mathbf v_t\tau_\mathrm{outer}$だけ進め、x/yをprimary periodic cellへwrapします。

このmodelは低コストですが、外部電位を指数profileへ固定する近似です。自己整合なVDFや非線形sheathが必要なら
`kinetic_1d`を使います。

## 4. `kinetic_1d`: 離散sheath profileでreturnを求める

ambient electron/ionのVDFから非線形1D Poisson問題を解き、batch開始時に収束した$\phi(z)$を流入と流出の
両方へ使うmodelです。

```toml
[outer_plasma]
model = "kinetic_1d"
return_model = "kinetic_1d_profile_return"

[coupling]
particle_transfer_mode = "electrostatic_1d_instant_return"
```

### 離散kinetic profile上でturning pointを探す

各grid点で

$$
v_n^2(z)=v_{n,I}^2+\frac{2q[\phi_I-\phi(z)]}{m}
$$

を評価します。interface から profile を外向きに走査し、最初に$v_n^2$の符号が変わる区間があれば turning point を
線形補間して return とします。離散 profile と far Robin tail の全経路を通過できる場合だけ escape と判定します。

正の速度を持つ区間$a\to b$の片道時間は

$$
\Delta t=\frac{2\Delta z}{\sqrt{v_{n,a}^2}+\sqrt{v_{n,b}^2}}
$$

で積算します。turning区間では到達fractionまでを積分し、grid上端より先にturning pointがある場合は
far Robin exponential tailを解析積分します。往復後の速度反転、接線変位、periodic wrapは`linear_debye`と同じです。

Zhao closureの現実装は、このgrid外積分にもphotoelectron reference Debye長をtail scaleとして使います。真の漸近scale
`abs(phi_end/E_end)`とは一般に一致しないため、separatrix近傍のflight timeは暫定値です。production利用前にこの差を
評価し、必要ならouter stateへ独立なtail scaleを追加します。

### 非単調 Type A profile の障壁を全経路で調べる

`kinetic_closure="zhao_charge_driven"` の Type A profile は途中に potential minimum を持つため、interface と
無限遠の電位差だけでは粒子の到達可能性を決められません。外向き粒子は profile の全区間を z の順に走査し、
最初に到達不能になる区間を turning point とします。

無限遠からの流入では、無限遠での法線速度 $v_{n,\infty}$ に対して全 grid 点で

$$
v_n^2(z)=v_{n,\infty}^2+\frac{2q[\phi_\infty-\phi(z)]}{m}
$$

を検査します。どこか 1 点でも負になれば interface へ到達不能とします。これは endpoint の電位差ではなく、
経路上で最も厳しい potential-energy barrier を流入 cutoff に使うことに相当します。

### profileの物理解を検査する

profileは有限値、z座標の狭義単調増加、interface点との整合に加え、closureとresolved branchごとの電位shapeを検査します。

| closure / branch | 受理する電位shape |
| --- | --- |
| `absorbing_maxwellian` | 全gridで単調増加または単調減少 |
| Zhao `A` | interior minimumまで非増加、その後は非減少 |
| Zhao `B` | 全gridで非増加 |
| Zhao `C` | 全gridで非減少 |
| Zhao `0` | flat bootstrap |

この条件外の非単調profileは受理しません。
物理的 turning point を bracket しない区間や、必要な Robin tail が非正の場合は invalid model として停止します。
自己整合な平均 sheath を扱える一方、平面平均された
静電・無衝突・非磁化1D profileという仮定を持ちます。詳しくは
[外部場: kinetic 1D](KineticOuterPlasma.html)で説明します。

## 5. `unified_linear_response`: 外部3D軌道を積分する

rough surfaceからfar planeまでのzero modeとscreened nonzero modeを一つの線形応答場として構築するmodelです。
`unified_linear_response`だけでは粒子returnは有効になりません。次の3D orbit設定を組み合わせたときに、
外部粒子を明示的に追跡します。

```toml
[outer_plasma]
model = "unified_linear_response"
return_model = "electrostatic_3d_explicit_orbit"

[coupling]
particle_transfer_mode = "electrostatic_3d_explicit_orbit"
```

### unified 3D field中で外部軌道を進める

batch内で固定された電場に対し、固定刻み`outer_orbit_dt`のvelocity-Verlet更新を使います。

$$
\mathbf v^{n+1/2}=\mathbf v^n+\frac{q\mathbf E(\mathbf x^n)}{2m}\Delta t_o,
$$

$$
\mathbf x^{n+1}=\mathbf x^n+\mathbf v^{n+1/2}\Delta t_o,
\qquad
\mathbf v^{n+1}=\mathbf v^{n+1/2}+\frac{q\mathbf E(\mathbf x^{n+1})}{2m}\Delta t_o
$$

です。ownership interfaceを内向きに再通過すればreturn、unified grid上端のfar planeを外向きに通過すれば
escapeです。境界到達位置と速度は最後のouter step内で線形補間します。

### 軌道の収束とenergyを検査する

`outer_orbit_max_steps`までにどちらの境界にも到達しない場合は、persistent queueが必要になるため停止します。
また、初期状態とreturn/escape境界へ到達したときの全エネルギー

$$
\mathcal E=\frac12m|\mathbf v|^2+q\phi(\mathbf x)
$$

の相対誤差が`outer_orbit_energy_tolerance`を越えても停止します。`outer_orbit_dt`はreturn/escape率、flight time、
energy errorに対して収束確認します。場の構成と適用範囲は
[外部場: unified linear response](UnifiedLinearResponse.html)で説明します。

## outer transferに共通する処理

### z-highを粒子ownershipのinterfaceにする

`coupling.particle_transfer_mode`が有効なz-high open面だけが、交差情報をouter modelへ渡します。

| `particle_transfer_mode` | 対応するmodel | 粒子処理 |
| --- | --- | --- |
| `none` | outer transferなし | 通常のopen境界 |
| `electrostatic_1d_instant_return` | `linear_debye`または`kinetic_1d` | 保存energyからescape/returnを直接写像 |
| `electrostatic_3d_explicit_orbit` | `unified_linear_response` | batch内で固定された3D場で外部軌道を時間積分 |

1D/3D transferはいずれもz-highがopenであることと`sim.b0=0`を要求します。現行modelは外部領域での磁場軌道を
扱いません。

### 境界通過時の状態を外部modelへ渡す

[粒子の衝突・境界イベント](ParticleEvents.html)は、mesh hitとbox面通過のうち軌道上で最初のeventを選びます。
outer modelへ渡す交差情報は次のとおりです。

- interface上の位置と外向き速度。
- local Boris step内の交差時刻。
- 交差後に残る`dt_remaining`。

returnした粒子はinterfaceの直内側へ戻し、`dt_remaining`の間だけ通常のBoris更新とevent処理をやり直します。
outer flight timeは、このlocal step remainderとは別の診断量です。

### outer flightをglobal timeへ加えない近似

1D returnの「instant」は、outer flightを状態写像には使う一方、global simulation timeを進めないという意味です。
3D explicit orbitも現行実装では同じ規約です。粒子はinterfaceを出たlocal stepと同じsimulation時刻へ戻り、
outward/returned chargeは同じbatchに計上されます。

この近似は定常または準定常outer plasmaを対象にします。UV照射の開始、plasma条件の急変、短pulseへの過渡応答には、
過去のoutgoing currentを保持するdelayed-return queueが必要です。

### outer flight中にfieldを固定できるか検査する

snapshotの有効性を

$$
\epsilon_\mathrm{ad}
=\frac{\tau_\mathrm{outer}}{\texttt{field\_evolution\_timescale}}
$$

で評価し、`max_frozen_field_ratio`以下を要求します。$\tau_\mathrm{outer}/\mathrm{batch\_duration}\gtrsim1$では、
十分に定常化した長時間平均は使えますが、batchごとのreturn currentは正しい時間変化を表しません。

persistent delayed-return queueは未実装です。`outer_queue_enabled=true`は設定検証で拒否され、queueが必要な軌道は
停止します。

## 診断量でmodelの結果を確認する

species別に`interface_outward_gross`、`interface_returned_gross`、`escaped_to_infinity`を区別します。さらに、最大
`outer_flight_time`、frozen-field ratio、3D orbitのenergy errorを出力します。

gross outwardとreturnedはどちらもinterface通過量です。その差が正味escapeと一致するのは、そのspeciesの
transfer対象と電荷収支の集計期間が一致する場合に限ります。各項目と`charge_ledger.csv`の読み方は
[出力ファイルを調べる](OutputGuide.html)にまとめています。

## Code reference

- `escape`と`potential_barrier`: [`bem_particle_stepper.f90`](../src/runtime/simulator/bem_particle_stepper.f90)
- `linear_debye`と`kinetic_1d`: [`bem_outer_plasma_interface.f90`](../src/physics/outer_plasma/bem_outer_plasma_interface.f90)
- `unified_linear_response`の3D軌道: [`bem_outer_plasma_orbit.f90`](../src/physics/outer_plasma/bem_outer_plasma_orbit.f90)
- interface transferとdiagnostic集計: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- model組合せの検証: [`bem_physics_config_types.f90`](../src/config/bem_physics_config_types.f90)

title: 粒子のescapeとreturn

Lang: [日本語](ParticleEscapeReturn.md) | [English](ParticleEscapeReturn.en.md)

# 粒子のescapeとreturn

粒子が open 面へ到達したときは、最初に「z-high の外部領域がその面を所有するか」を判定します。
設定も同じ責務に合わせて2段に分けます。

1. outer が所有しない面は `external_boundary.ordinary_open.model` で処理する。
2. `particles.mode="same_batch"` または `"zhao_queue"` が所有する z-high は、
   `external_boundary.field.model` に対応する外部軌道で処理する。

```text
open面を粒子が通過
├─ outerが所有しない面
│  ├─ escape                    無条件に除去
│  └─ potential_barrier         scalar障壁で反射またはescape
└─ particles.modeが所有するz-high
   └─ kinetic_1d                離散1D profile return
```

通常 open 面の処理と kinetic outer transfer は別の責務です。通常 open 面には次の2択があります。

| `ordinary_open.model` | 外部状態 | 粒子の判定方法 | 主な用途 |
| --- | --- | --- | --- |
| `escape` | なし | 常にescape | 単純な有限box |
| `potential_barrier` | 通過点のscalar電位 | 法線energyと障壁を比較 | 低コストな局所反射 |

z-high を outer が所有する場合は `particles.mode` で追跡方法を選び、具体的な軌道は field から導出します。

| `field.model` | 外部状態 | 粒子の判定方法 | return時間 | 主な用途 |
| --- | --- | --- | --- | --- |
| `kinetic_1d` | 収束した離散1D profile | 保存energyとturning point探索 | profileを積分 | **標準:** 自己整合な平均sheath |

どの境界処理も粒子sourceには依存しません。reservoir粒子、光電子、`volume_seed`粒子は、同じ状態で同じ面を
横切れば同じ境界処理を受けます。外部場自体の選び方は
[境界・外部領域の構成を選ぶ](OuterPlasmaModels.html)にまとめています。

## 1. `escape`: open面で粒子を除去する

最も単純なmodelです。outer transfer対象でないopen面に対して次を指定します。

```toml
[external_boundary.field]
model = "none"

[external_boundary.particles]
mode = "local_source"

[external_boundary.ordinary_open]
model = "escape"
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
phi_infty = 0.0

[external_boundary.field]
model = "none"

[external_boundary.particles]
mode = "local_source"

[external_boundary.ordinary_open]
model = "potential_barrier"
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
`ambiguous_open_corner`で停止します。

通過点電位は粒子運動と同じsnapshot規約で評価するため`sim.e0`の局所電位も含みます。一様電場には有限な
無限遠電位がないため、`sim.e0`と併用する場合の`phi_infty`は有効なreservoir基準電位としてユーザが
整合させる必要があります。

## 3. `kinetic_1d`: 離散sheath profileでreturnを求める

ambient electron/ionのVDFから非線形1D Poisson問題を解き、batch開始時に収束した$\phi(z)$を流入と流出の
両方へ使うmodelです。

```toml
[external_boundary.field]
model = "kinetic_1d"
debye_length = 0.2
thermal_voltage = 2.0

[external_boundary.particles]
mode = "same_batch"
field_evolution_timescale = 1.0
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
far Robin exponential tailを解析積分します。往復後は法線速度を反転し、接線変位を加えてperiodic wrapします。

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
Robin tailを使うのはinstant経路であり、Zhao queue経路は後述する有限$L$境界で終了します。
自己整合な平均 sheath を扱える一方、平面平均された
静電・無衝突・非磁化1D profileという仮定を持ちます。詳しくは
[外部場: kinetic 1D](KineticOuterPlasma.html)で説明します。

## outer transferに共通する処理

### z-highを粒子ownershipのinterfaceにする

`external_boundary.particles.mode`が所有するz-high open面だけが、交差情報をouter modelへ渡します。

| `particles.mode` | 対応するfield | 粒子処理 |
| --- | --- | --- |
| `local_source` | すべて | 通常のopen境界 |
| `same_batch` | `kinetic_1d` | 保存energyからescape/returnを写像 |
| `zhao_queue` | Zhao `kinetic_1d` | 1D結果を後続batchのeventへ遅延 |

1D transferはz-highがopenであることと`sim.b0=0`を要求します。現行modelは外部領域での磁場軌道を扱いません。

### 境界通過時の状態を外部modelへ渡す

[粒子の衝突・境界イベント](ParticleEvents.html)は、mesh hitとbox面通過のうち軌道上で最初のeventを選びます。
outer modelへ渡す交差情報は次のとおりです。

- interface上の位置と外向き速度。
- local Boris step内の交差時刻。
- 交差後に残る`dt_remaining`。

instant modeでreturnした粒子はinterfaceの直内側へ戻し、`dt_remaining`の間だけ通常のBoris更新とevent処理を
やり直します。その残り時間内にz-highへ再到達した場合もouter transferへ戻し、local stepが完了するまで
returnと再交差を反復します。1 local stepあたり8回を上限とし、9回目は状態をcommitせずfail closedにします。
queue modeでは構成済みのreturn位置・速度をevent recordへ保存し、dueとなったbatchのfresh sourceの
後ろへ追加してlocal-domain粒子として再開します。outer flight timeは、local step remainderとは別の診断量です。

z-highの二次軌道補正後も周期または反射面を同時に横切るcornerでは、横方向の周期wrapまたは反射を先に合成してから
outer modelへ渡します。補正によってz-highと横面の先後・同時関係がどちら向きに変わっても作用順序を推測しません。
z-highと別のopen面を同時に横切ってownerが一意でない場合と同様にfail closedにします。

## outer flightをglobal timeへ加えない近似

1D returnの「instant」は、outer flightを状態写像には使う一方、global simulation timeを進めないという意味です。
粒子はinterfaceを出たlocal stepと同じsimulation時刻へ戻り、outward/returned chargeは同じbatchに計上されます。

この近似は定常または準定常outer plasmaを対象にします。UV照射の開始、plasma条件の急変、短pulseへの過渡応答には、
次節のdelayed-return queueを使います。

### 定常 Zhao profileからinstant経路を開始する

`steady_start_mode="zhao_floating"`は、Zhao 零電流定常根とその$E_I$に対応する平面電荷を新規実行の初期値にします。
この初期化はreturn algorithmを変更しません。`phi(infinity)=0`へ接続した同一kinetic profileを、無限遠reservoirからの
流入barrierとinterface外のinstant escape/returnの両方に使います。resumeではcheckpointのprofileとmesh電荷を復元し、
定常根から再seedしません。このmodeは未帯電状態からの過渡を短絡するため、return-current遅延の評価には使いません。
設定と電荷式は[kinetic 1D外部プラズマ](KineticOuterPlasma.html#定常研究を-zhao-零電流根から開始する)を参照してください。

### Zhao 過渡closureでouter flightをqueueする

強いUV照射の開始など、outer flightの遅延をbatch履歴へ反映する場合は次のZhao構成を使います。

```toml
[sim]
batch_duration = 2.5e-7

[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "zhao_charge_driven"
debye_length = 10.5
thermal_voltage = 10.0

[external_boundary.particles]
mode = "zhao_queue"
field_evolution_timescale = 2.0e-5
max_frozen_field_ratio = 0.2
```

`zhao_queue`は`zhao_charge_driven`専用で、`absorbing_maxwellian`との組合せは拒否します。
さらに`batch_duration <= max_frozen_field_ratio * field_evolution_timescale`を要求します。
`particles.mode="zhao_queue"`は、tracked粒子のouter flightとZhao光電子populationを一つの保存inventoryで接続します。
各rankはeventをlocal queueに保持し、Zhao closureへの入力として光電子のmacro粒子数をMPI全体で合計します。
batch開始時にdue eventを除いた
光電子column targetは

$$
N_{pe,q}=\frac{1}{A_{xy}}
\sum_{j\in\text{queued photoelectron}}w_j
$$

です。Zhao solverは有限control volume $L=10\lambda_{D,pe}$について

$$
N_{pe,Zhao}(\eta)=
\int_0^L\left[n_{pe,f}(z;\eta)+n_{pe,c}(z;\eta)\right]dz
=N_{pe,q}
$$

を満たすpopulation scale $\eta$を解きます。`outer_photoelectron_population_fraction`という出力名ですが、$\eta$は
確率ではなく定常reference populationに対するoccupancy scaleです。$\eta=0$から連続する物理解を$0\le\eta\le16$で
探索し、1を超える一時的overshootも許します。`[0,1]`へclampせず、targetを無視したfull-population解や
disconnected branchへjump/fallbackしません。queue modeは`zhao_branch="auto"`を要求し、縮退条件を満たす連続的なA/B等の
branch遷移だけを許します。現在はcolumnが$\eta$とともに単調増加するpathだけを受理し、foldは未対応です。targetへ到達する
連結・単調pathがなければ`no_physical_solution`で停止します。targetが0なら$\eta=0$を厳密に使います。
$\eta$はphotoelectron密度、無限遠準中性、Sagdeev項をscaleしますが、current診断のraw photoelectron
emission-current項はscaleしません。このanalytic raw currentはtracked sourceの整合性検査とcurrent-density診断に使い、
root、surface charge、ledgerへ別途加えません。

同じ$0\le z\le L$がqueueの粒子ownership領域です。離散profile上で$L$より手前にturning pointがあればreturn event、
$L$へ到達すれば外部reservoirに吸収されたescape eventとしてqueueへ入れます。queue modeは$L$外のRobin tailを延長して
returnを判定しません。

1 batchの時刻更新は次の順です。

1. batch $b$の開始時刻$t_b=(b-1)\Delta t_b$までにdueとなったrank-local eventをpopする。
2. pop後のglobal光電子inventoryから$N_{pe,q}$を計算し、Zhao profileと$\eta$を更新する。
3. fresh sourceを生成し、due return粒子を追加してlocal particle loopを実行する。due escapeはこのbatchの
   `escaped_to_infinity`へ計上する。
4. z-highを外向きに通過した粒子について、現在のprofileでinterface returnまたは$L$へのreservoir escapeと
   $\tau_{outer}$を計算し、
   batch中央を通過時刻とする
   $t_{due}=(b-\tfrac12)\Delta t_b+\tau_{outer}$でenqueueする。
5. surface chargeをcommitし、post-enqueue inventoryでZhao profileと$\eta$をcorrector更新する。この状態を次batchの
   continuation seedとcheckpointに使うため、straight runとsplit-resumeは同じbatchごとの演算列を通る。

BEACHのbatch内粒子は共通の物理時刻で同期していないため、crossing時刻にはbatch中央を使います。due eventはbatch開始時
だけreleaseされるため、return/escape時刻も`batch_duration`単位に量子化されます。enqueue時に決めたterminal状態と
due時刻は、その後にouter fieldが変わっても再積分しません。このclosureはflight delayとouter光電子columnを表しますが、
時間依存Vlasov–Poisson解、outer粒子間衝突、energy-resolved cloud evolutionではありません。

`batch_duration`を半分、`batch_count`を2倍にして同じ終了時刻を比較し、少なくとも$\eta$、column residual、return/escape
current、表面電荷、離脱力の収束を確認します。profileは固定128点の$0\le z\le10\lambda_{D,pe}$へ再標本化するため、
productionのcolumn grid収束にはgrid点数の入力化も必要です。tracked粒子数、horizontal area、effective interface位置も
独立に確認してください。

queue stateはserialでは`outer_event_queue.csv`、MPIでは各rankの`outer_event_queue_rankNNNNN.csv`へ保存します。
active eventのphase-space状態、terminal outcome、due時刻、`next_event_id`を保持します。resumeではqueue有効時に全rank分を
必須とします。queue-file schema 2はrank-local内容の`local_fingerprint`を持ち、summaryの
`outer_queue_fingerprint`は全rankのqueue内容と順序を束縛します。schema、rank、world size、完了batchに加えてglobal
event count、signed charge、fingerprintが一致しないcheckpointをfail closedで拒否します。

### outer flight中にfieldを固定できるか検査する

通常の instant 粒子では snapshot の有効性を

$$
\epsilon_\mathrm{ad}
=\frac{\tau_\mathrm{outer}}{\texttt{field\_evolution\_timescale}}
$$

で評価し、`max_frozen_field_ratio`以下を要求します。`implicit_mean` の deferred PE は、scalar closure が決めた
return 総量の軌道・着地点を標本化する準定常 shadow です。同じ $\epsilon_\mathrm{ad}$ を診断へ残しますが、
上限超過を停止条件にしません。通常の `same_batch` 粒子と ambient species の fail-closed 判定は変わりません。

`implicit_mean` のうち Zhao closure は、ambient linear-Debye closure の一様な
`return_weight_scale`を使いません。batch開始場でz-highへ初回到達した光電子から、法線運動energyと符号付きmacro電荷を
集めます。そのcharge-weighted実測CDFと、Zhao profile全体から得るvirtual cathodeを含む障壁を用いて、
$Q(\Phi_I)$を非線形に解きます。各energy群のreturn重みをその解から直接決め、同じ受理profileでouter軌道と
局所領域への帰還先を1回だけ追跡します。

この写像では、各deferred rayの外向きinterface通過が1回、return通過が高々1回であることを要求します。
正のreturn重みを持つrayはsurfaceへ再吸収されなければならず、return重みが0の非境界rayはreservoir escapeで
終端しなければなりません。再越境や終端channel不一致が有意なchargeを持つ場合、他rayへ重みを付け替えずfail closedで
停止します。障壁と一致するmarginal energy群だけは、1本の代表rayをreturn側に置いたままescape fraction $f$ と
return fraction $1-f$へ統計的に分割します。これは粒子軌道を途中で分岐させる操作ではありません。

$\tau_\mathrm{outer}/\mathrm{batch\_duration}\gtrsim1$では、
十分に定常化した長時間平均は使えますが、batchごとのreturn currentは正しい時間変化を表しません。
特に `implicit_mean` は UV turn-on の delayed return current や batch 間 outer inventory を解きません。
それらには、この closure に対応する delayed inventory / queue が別途必要です。

Zhao queue modeでも`field_evolution_timescale`と`max_frozen_field_ratio`は正値を要求します。eventはbatch開始時だけ
releaseされるため、$t_{due}$から最初のbatch-start pollまでの量子化遅延$\delta_{poll}$と、batch内crossing時刻の
midpoint近似誤差上限$\Delta t_b/2$を含め、

$$
\frac{\tau_{outer}+\delta_{poll}+\Delta t_b/2}{\texttt{field\_evolution\_timescale}}
\le\texttt{max\_frozen\_field\_ratio}
$$

を要求します。超過eventはdiscardやinstant fallbackをせず、enqueue前に停止します。さらに1 batchの幅についても

$$
\texttt{batch\_duration}
\le
\texttt{max\_frozen\_field\_ratio}\,
\texttt{field\_evolution\_timescale}
$$

を設定検証で要求し、違反時は停止します。

## 診断量でmodelの結果を確認する

species別に`interface_outward_gross`、`interface_returned_gross`、`escaped_to_infinity`を区別します。さらに、最大
`outer_flight_time`、frozen-field ratio、kinetic 1D return / escape写像の
`max_outer_energy_relative_error`を出力します。最後の値は法線運動エネルギーと静電エネルギーの保存残差を規格化した
診断です。`implicit_mean` PE shadow の ratio は設定上限を超え得ますが、それ自体は停止理由ではありません。

Zhao queue modeでは`summary.txt`の`outer_photoelectron_population_fraction`、
`outer_photoelectron_column_per_area_m2`、`outer_photoelectron_column_target_per_area_m2`、
`outer_photoelectron_column_residual_per_area_m2`、`outer_queue_event_count`、`outer_queue_signed_charge_C`を確認します。
`charge_ledger_outer_flight_charge_before_C`と`charge_ledger_outer_flight_charge_after_C`はqueue stockを電荷保存残差へ
組み込みます。

gross outwardとreturnedはどちらもinterface通過量です。その差が正味escapeと一致するのは、そのspeciesの
transfer対象と電荷収支の集計期間が一致する場合に限ります。各項目と`charge_ledger.csv`の読み方は
[出力ファイルを調べる](OutputGuide.html)にまとめています。

## Code reference

- `escape`と`potential_barrier`: [`bem_particle_stepper.f90`](../src/runtime/simulator/bem_particle_stepper.f90)
- `kinetic_1d`: [`bem_outer_plasma_interface.f90`](../src/physics/outer_plasma/bem_outer_plasma_interface.f90)
- delayed event queue: [`bem_outer_event_queue.f90`](../src/runtime/coupling/bem_outer_event_queue.f90)
- queue checkpoint: [`bem_outer_event_queue_io.f90`](../src/runtime/coupling/bem_outer_event_queue_io.f90)
- interface transferとdiagnostic集計: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- model組合せの検証: [`bem_physics_config_types.f90`](../src/config/bem_physics_config_types.f90)

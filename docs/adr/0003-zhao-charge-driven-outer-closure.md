# ADR 0003: charge-driven Zhao outer closure

- Status: experimental
- Date: 2026-07-21

## Context

既存の Zhao 系 `sheath_injection_model` は、定常電流バランスから求めた解析 profile で source VDF を
事前補正する。Boris pusher が使う場や蓄積表面電荷からは更新されない。一方、periodic2 の
`kinetic_1d` は蓄積電荷が作る z-high interface 電場を境界条件として outer profile を更新するが、
従来の half-Maxwellian closure は Type A の virtual cathode を表せない。

UV 照射開始直後や帯電途中に毎 batch `J=0` を課すと、解くべき表面電荷の時間発展を先に消してしまう。
必要なのは、現在の interface 電場を満たす準定常 profile と、その profile で評価した非零でもよい電流である。

## Decision

`outer_plasma.model="kinetic_1d"`へ選択可能な
`kinetic_closure="zhao_charge_driven"`を追加する。`zhao_branch`は`auto`、`a`、`b`、`c`を取る。
既定の`absorbing_maxwellian`は変更しない。legacy Zhao injection、reservoir potential correction、
`photoelectron_density_model="kinetic_mean"`との重複適用は拒否する。

UV source強度はouter-cloud occupancy $\eta$と分離し、`outer_plasma.photoelectron_source_scale=s_UV`で指定する。
既定値は1で、`s_UV=0`をno-photo経路とする。no-photoでは光電子species、$T_{pe}$、$n_{pe,ref}$を要求せず、
ambient $T_e$と準中性密度$n_\infty$をZhao式の電位・密度scaleへ使う。$E_I=0$は解析的な平坦Type B/C接続点、
$E_I<0$はType Cである。surface chargeが決める電場をroot条件とし、zero-currentは定常性の診断に残す。

Zhaoの無次元電位と電場を

```text
psi = phi / T_pe[eV]
E_hat = E_I lambda_D,pe / T_pe[V]
```

とする。Type B/Cでは無限遠準中性条件と

```text
2 integral_(psi_0)^0 rho_hat(psi) dpsi = E_hat^2
```

を解く。Type Aでは無限遠準中性条件、upper branchの`E(infinity)=0`、およびlower branchの

```text
-2 integral_(psi_m)^(psi_0) rho_hat_lower(psi) dpsi = E_hat^2
```

を解く。旧 Zhao の zero-current 式は root residual に使わず、electron、ion、photoelectron、total の
current-density診断へ変換する。従って帯電途中のtotal currentは非零でよい。定常化すれば旧 Zhao rootを
包含する状態として`J -> 0`へ近づける。

初期表面電荷がゼロのときは`E_I=0`である。正のambient electron密度で無限遠準中性を満たせる場合、`auto`は
`phi_0=phi_m=0`の退化Type Bを解析的に構成する。強い光電子源のためその密度が負になる場合は、outer
photoelectron populationがまだ形成されていない初期瞬間として、平坦なambient-only bootstrap stateを1回だけ使う。
このstateはZhao準定常rootではなくresolved branchを`0`と記録する。どちらの開始状態でも電流診断は一般に非零であり、
最初のbatchからtracked emission/depositionが表面電荷を更新する。非零$E_I$になった後はbootstrapへfallbackしない。

full photoelectron populationの準定常branchは一般に$E_I=0$へ連続しない。flat-plane fixtureの条件では、正電場の
Type Bは約$1.339\,\mathrm{V/m}$以上、Type Aは約$1.449$--$1.752\,\mathrm{V/m}$にだけ存在し、
$0<E_I<1.339\,\mathrm{V/m}$にはA/B解がない。これはNewton seed不足ではなく、無限遠準中性から要求される正の
ambient electron密度とSagdeev積分の可解域である。従って細かいbatchで未帯電状態から追跡するには、outer
photoelectron cloud occupancyなどを時間発展させる別の過渡closureが必要になる。現実装はこのgapを定常解で補間せず、
`outer_queue_enabled=false`では`no_physical_solution`として停止する。後から追加した有限column inventoryとevent queueによる
過渡closureは[ADR 0004](0004-zhao-transient-photoelectron-column-queue.md)で定義する。

Type A profileは途中にpotential minimumを持つ。粒子のescape/returnと無限遠からの流入はendpointの
電位差だけでは分類せず、離散profile全体を走査して最初のturning pointまたは最大障壁を使う。

## Interface interpretation

初版は z-high ownership interface を Zhao の有効放出面、かつType A lower sideとして扱う。このためType Aは
`E_I>0`、Type Bは`E_I>=0`、Type Cは`E_I<0`を要求する。interfaceが実際の光電子放出面から有限距離上にあり、
Type A minimumより上側に来る一般の場合は、この境界条件だけでは扱わない。

光電子ありでは`sheath_photoelectron_ref_density_cm3`と`sheath_alpha_deg`を既存 Zhao 入力として再利用する。solver profileの
物理scaleは、photoelectron温度$T_{pe}$と、$T_{pe}$および$n_{ref}$から導出した$\lambda_{D,pe}$である。
`outer_plasma.debye_length`と`outer_plasma.thermal_voltage`はZhao root/profileには入らない。ただし現時点では、
`interface_eta_gap`、横方向phi/field比、local-charge推定などsplit-interface適用性診断のreference inputとして残す。
ambient electron、ionはそれぞれenabledな負電荷z-high `reservoir_face`、正電荷z-high `reservoir_face` species
ちょうど1つから構成し、準中性密度を要求する。`s_UV>0`では負電荷`photo_raycast` speciesもちょうど1つ要求し、
`s_UV=0`ではenabledな光電子sourceとqueueを拒否する。
electron/ionは`sheath_electron_drift_mode="normal"`と`sheath_ion_drift_mode="normal"`だけを受理する。
photoelectronは`normal_drift_speed=0`、ionはcold-ion近似$T_i\le0.1T_e$を要求し、外れる状態を`not_applicable`とする。

tracked `photo_raycast`のraw source currentは、有効平面の解析値

```text
v_th,pe = sqrt(2 T_pe[J] / m_pe)
J_pe,raw = s_UV |q_pe| n_ref sin(alpha) v_th,pe / (2 sqrt(pi))
```

と1%以内で一致させ、一致しない設定をruntimeで拒否する。解析currentはrootの電流診断へだけ使い、表面電荷へ
二重加算しない。tracked `photo_raycast`粒子の放出・再吸収だけが従来どおり表面電荷を更新する。

初版はray方向やrough surfaceからinterfaceへ到達するVDFをZhao outer populationへ自己無撞着に接続しない。
`ray_direction`は照射rayによる放出面sampling、`sheath_alpha_deg`は解析Zhao sourceを決める独立の入力である。
従ってflat-plane fixtureでvertical rayを選んでも、$\alpha=0$を意味しない。

## Numerical and failure policy

- 符号制約を保つlog変数で2または3変数のNewton solveを行う。
- Sagdeev積分、準中性、Type A far-field残差を無次元化する。
- `auto`は要求電場の符号と太陽高度を使ってbranchを探索し、branchが存在しなければ停止する。
- 解なしは`no_physical_solution`、反復失敗は`numerical_failure`とし、旧注入モデルや線形modelへfallbackしない。
- profile、Gauss則、current density、解branchを出力し、設定はmodel fingerprintへ含める。

## Consequences

蓄積電荷から求めた電場と、Zhaoのfree/reflected/captured populationを同じouter stateで使える。Type Aを
含む準定常particle returnも同一profileで判定できる。一方、これは平面・無衝突・非磁化の有効面modelであり、
rough surfaceからinterfaceへ実際に通過した光電子のVDFを境界条件にするものではない。
またgrid外return時間は現状$\lambda_{D,pe}$のRobin tailを仮定し、Zhao endpointの`abs(phi/E)`を独立には保持しないため、
separatrix近傍のflight timeには系統誤差が残る。

production適用前には、有効interface位置、profile grid、tracked粒子数、`dt`/batch解像度の収束、outer flight timeと
帯電時間の比、tracked crossing histogramから得たflux/VDFとの比較が必要である。Zhao profileの収束確認に
`debye_length`や`thermal_voltage`を使わない。今回のflat-plane pilotは実装smokeであり、一般rough-surface closureの
妥当性を示す証拠にはしない。
初期の付属flat-plane exampleは配線検証のため、1 coarse batchで既知のType A可解域へ着地するようweightを選んだ
branch-entry smokeだった。現在のexampleは[ADR 0006](0006-zhao-stationary-warm-start.md)の明示的な定常warm startに置き換え、
通常のsmall batchを使う。どちらも未帯電過渡の時間積分やbatch収束を表さない。

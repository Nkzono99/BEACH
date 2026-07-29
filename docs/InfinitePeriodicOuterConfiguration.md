title: periodic2無限周期＋outer plasma構成

Lang: [日本語](InfinitePeriodicOuterConfiguration.md) | [English](InfinitePeriodicOuterConfiguration.en.md)

# periodic2無限周期＋outer plasma構成

この構成は、x/y無限周期のsurface fieldとz方向の外部plasmaモデルを一つの電場・電位として組み立てます。
near image、Ewaldから生成したfar `k\ne0` operator、物理`k=0`、outer responseを重複なく合成し、同じouter potentialを
reservoir inflowとparticle returnに使います。

対応する外部シース構成は、species VDF または ambient 線形応答から平均シースを作る split `kinetic_1d` です。

## 各field成分に一つのownerを割り当てる

| component | owner |
| --- | --- |
| primary/near image | FMMまたはspectral base field |
| infinite-periodic far `k\ne0` | `cached_kneq0`、小規模検証では`panel_spectral_reference` |
| surface `k=0` | periodic zero-mode plan/state |
| outer mean response | `kinetic_1d` |
| reservoir velocity map | outer interface potential差 |
| outward escape/return | 同じouter profile |

periodic operatorと`k=0`の式は[periodic2静電場](PeriodicElectrostatics.html)で説明します。

## split kineticの責務を確認する

| 構成 | mean plasma | nonzero mode | particle transfer |
| --- | --- | --- | --- |
| split kinetic | VDFに基づく密度モデルを使う非線形`kinetic_1d` | surface側でinterfaceまで減衰すると仮定 | `particles.mode="local_source"`、または`"same_batch"`の1D profile |

split kineticはspecies別VDF、Bohm entry、photoelectron mean densityを扱えます。一方で、surfaceとinterfaceの間には、
local surface fieldだけを扱うsplit windowを仮定します。[kinetic 1D外部プラズマ](KineticOuterPlasma.html)に
この分割と非線形solveをまとめています。

## batch開始時に確定した場を流入とreturnで共有する

1. commit済み`q_elem`からFMM/source multipoleとsurface zero modeをrefreshする。
2. cached far operatorをcurrent root multipoleへ適用する。
3. 必要なstrideでinterface fieldからouter profileを解く。
4. potential gauge、Gauss residual、interface diagnosticsを更新する。
5. outer stateを使ってz-high reservoirのglobal粒子数とinterface速度を決める。
6. photo raycastを行い、source reaction chargeをbatch差分へ記録する。
7. 全粒子をbatch内で固定された場の中で追跡する。
8. z-high outward crossingを同じouter stateでescape/returnへ写像する。
9. surface absorptionとemission差分をMPI all-reduceし、batch末尾にcommitする。

operator や outer Poisson problem は、particle hit ごとに解き直しません。通常の explicit 経路では source と return が
同じ batch-start outer state を共有し、surface charge は次の batch から field に反映されます。
後述する `ambient_linear_debye` の陰的平均経路だけは、最初の粒子 trace 後も $k\ne0$ を batch-start 値に保ちながら、
連続 Maxwellian mean solve で $k=0$ と outer state を更新し、energy-resolved ray を 1 回だけ局所再分配へ使います。

## 同じenergy式をkinetic流入とreturnに使う

z-highの最初の負・正`reservoir_face` speciesを無限遠electron/ion VDFとして使います。収束profileの
$\phi_I-\phi_\infty$により

- inflowでは到達可能な$v_{n,\infty}$を選び、$v_{n,I}$へ加減速する。
- outflowでは$v_{n,I}$から$v_{n,\infty}^2$を計算し、escape/turning returnを判定する。

tracked split kineticでは`external_boundary.particles.inflow_model="auto"`が同じprofileへ流入を委ねます。
このため、`inflow_model="infinity_barrier"`を重ねません。
後述の `implicit_mean` では、光電子の interface crossing を平均場の更新まで保留した後、同じ個別 outflow 写像へ渡します。
[`reservoir_face` の流入量と速度サンプリング](ReservoirInjection.html)が流入側、[粒子のescapeとreturn](ParticleEscapeReturn.html)が
流出側の写像を説明します。

## 平均outer密度とtracked光電子を別々に選ぶ

photoelectronの扱いはouter field modelとparticle returnを分けて選びます。

| 選択 | outer density | tracked particle |
| --- | --- | --- |
| `none`（`ambient_linear_debye` では内部 resolved 値、公開 key なし） | photoelectron mean densityなし | 通常のsource/軌道。`ambient_linear_debye + same_batch`ではinterfaceまで追跡 |
| `photoelectron_density_model="kinetic_mean"` + `particles.mode="local_source"` | 定常outgoing/returning mean density | z-highでは通常open処理 |
| `photoelectron_density_model="kinetic_mean"` + `particles.mode="same_batch"` | 同じmean density | 個々のinterface crossingもprofileでreturn/escape |

tracked returnでは`deposit_opposite_charge_on_emit=true`を要求します。平均密度モデルは
tracked surface depositを置き換えず、統計的return chargeを追加depositしません。
[光電子の放出とライフサイクル](PhotoelectronEmission.html)に、放出から再吸収までの電荷収支をまとめています。

## ambient 線形応答では局所 $k\ne0$ と平均 $k=0$ を別々に進める

`kinetic_closure="ambient_linear_debye"`、`particles.mode="same_batch"`、enabled な負電荷
`photo_raycast` species の組合せは、内部で `implicit_mean` に自動解決されます。公開設定に update mode は追加しません。

1. batch-start field で既存の 3D 粒子 trace を行い、interface 到達前の局所的な放出・再吸収と
   ambient 吸収を element へ stage し、光電子の z-high crossing を full macro weight と放出元情報付きで保留する。
2. 最初の trace が stage した全 surface charge delta を実測し、ambient electron / ion と PE 外向き成分を
   分離して残りを `J_other` として保存する。解析 Maxwellian backward-Euler solver で mean 総電荷、$\phi_I$、
   連続 escape / return 率を決める。
3. 解析 return 率を各 crossing の source countercharge へ掛けて発生元で一時中和し、
   mean solve の総電荷を持つ element 分布から $k=0$ と outer profile だけを refresh する。
4. 各 crossing を full weight のまま共通 kinetic 1D profile mapper へ渡し、法線 energy から return / escape を判定する。
   return は outer flight と横方向変位を含めて local 3D 領域へ戻し、再び z-high を横切れば同じ mapper へ渡す。
5. 解析 return 総電荷を、最終的に再吸収された ray sample の source leg（放出元での一時中和）と
   実 hit destination leg へ同じ scale で正規化配分し、零和 transaction を作る。
6. 正規化した実 deposit を通常の deposit と一度だけ commit し、更新後の $k=0$ と outer profile を次の batch へ渡す。

$k\ne0$ operator は batch-start 値に固定します。mean 総電荷と引戻し電位は同じ batch 内で更新する一方、
sampled 局所再分配は commit 後の batch から $k\ne0$ へ反映します。離散 ray 分類を mean solve へ戻さないため、
separatrix 近傍の ray による return / escape の 2-cycle を避けます。transaction の非零残差はエラーです。
`emit_current_density_a_m2` は tracked 3D 放出の weight を決めますが、PE 平均 source の振幅は実際の
interface crossing から測ります。標準出力の `J_other_A_m2` は、追加 species や他面・未解決 outcome を含む
残りの実測表面電流です。解析 Maxwellian closure が mean の return / escape 総電荷を決め、energy-resolved ray は
軌道と着地点の分布を sample します。`BEACH implicit-mean` 行は `transaction_residual_C`、
`mean_solver_iterations`、`sample_escape_fraction`、`return_weight_scale` を出力します。

この経路は `outer_update_stride=1`、正の `batch_duration`、`deposit_opposite_charge_on_emit=true` を要求し、
ちょうど 1 つの負電荷 `photo_raycast` species と `photo_raycast.normal_drift_speed=0` を要求します。
解析 Maxwellian escape 率は放出法線 drift を含みません。この経路は mesh mode に固有の制約を追加しません。
`photoelectron_density_model` は公開 TOML では省略し、内部で `none` に解決されます。
update mode や return kernel の公開設定は追加しません。
mean solver が収束しない、解析 return が正なのに再吸収 sample がない、transaction の電荷が釣り合わない、
または保留 ray が許可された trace 内に吸収・escape へ終端しない場合は fail closed で停止します。
保留した PE ray は準定常 shadow であり、outer flight time と frozen-field ratio は診断へ残しますが、
ratio 超過を停止条件にしません。
この ambient 線形経路では nonlinear photoelectron density と outer cloud inventory は未対応です。
物理式と診断は
[kinetic 1D 外部プラズマ](KineticOuterPlasma.html#ambient-線形-debye-応答と-tracked-光電子を分離する)
にまとめています。

ambient electron/ion reservoir がない UV-only 構成は、この無限遠 closure を満たしません。
`field.model="none"` と z-high escape を使う場合は有限 box の過渡 control であり、準中性な定常 outer sheath と
解釈しません。

## 強 PE では実測 energy CDF と非単調 Zhao profile を結ぶ

`kinetic_closure="zhao_charge_driven"`、`particles.mode="same_batch"`、負電荷 `photo_raycast`、
`steady_start_mode="zhao_floating"` の組合せも内部 `implicit_mean` を使います。公開 update-mode key は増えません。
この経路では ambient 線形 closure の Maxwellian return 率を使わず、各 rank の z-high crossing から
法線運動 energy と正の macro 電荷量を rank 0 へ集めます。charge-weighted empirical CDF と Zhao profile の
全経路 barrier $B(Q)$を使い、

$$
Q=Q_{\mathrm{base}}+Q_{\mathrm{escape}}\!\left[B(Q)\right]
$$

を解きます。同じ energy の macro 粒子は一群として扱い、separatrix が群の内部を通る場合だけ、その一群へ
fractional escape weight を与えます。解いた weight は MPI 全体へ broadcast し、各 crossing の source leg と
return後の実 hit destinationへ同じ零和 transactionとして適用します。`emit_current_density_a_m2` は3D表面放出量を
決め、Zhao mean source scale は実測 interface crossing current から別に解決します。

各候補 profile は、前の commit 済み root から

$$
G_b\!\left(y;E_I(\lambda),n_{\mathrm{pe0}}(\lambda)\right)=0
$$

を pseudo-arclength で追う同一 connected path 上に限定します。source scale と field の共通 anchorを作った後、
fixed-source sliceでは tangentから $dB/dQ\ge0$ を適応的に検査します。$\lambda$ tangentがゼロへ近づく fold、
fixed-parameter Jacobianのrank loss、barrierの低下、branch/domain endpointではfallbackせず停止します。
Type A/Bは$E_I>0$、Type Cは$E_I<0$を要求し、targetは通常のfixed-field Newtonではなく
$\lambda=1$のevent correctorで着地します。共通rootから$Q_0$と$Q_M$へ進むpathで全charge区間を検査し、
単調なorder predicateのfirst-true indexを二分探索します。pure root選択は$O(\log M)$です。
このlocal continuation guardは有限精度の数値判定であり、interval arithmeticによる大域的な数理証明では
ありません。最終候補では実測CDFの電荷残差も再検査し、order predicateまたはmarginal bracketが不整合なら
停止します。

$k\ne0$ operatorは最初の3D trace中はbatch-start値に固定され、解いたmean chargeとZhao outer stateを同じbatchで
更新します。個別rayの局所再分配が$k\ne0$へ反映されるのはcommit後です。保留rayはoutward crossingを1回、
returnを最大1回だけ許し、再越境や解いたweightとterminal outcomeの不一致はfail closedです。この準定常経路は
outer flight delayやbatch間outer inventoryを解きません。それらが必要な場合は`zhao_queue`を使います。

## productionではcached nonzero operatorを使う

FMM productionでは次のownershipを明示します。

```toml
[sim]
field_solver = "fmm"
field_bc_mode = "periodic2"
field_periodic_image_layers = 1
field_periodic_far_correction = "cached_kneq0"
field_periodic_ewald_layers = 4

[periodic2]
nonzero_mode_backend = "cached_kneq0"
zero_mode_policy = "exclude_k0"
lower_boundary_model = "symmetric_vacuum"
```

このblockだけではouter plasmaを有効にしません。`[external_boundary.field]`と
`[external_boundary.particles]`でkinetic構成を追加します。
最小のcache例は[`examples/periodic2_cached_panel.toml`](../examples/periodic2_cached_panel.toml)です。

小規模のreference計算にはDirect + `panel_spectral_reference`を使います。
[`examples/periodic2_kinetic_outer.toml`](../examples/periodic2_kinetic_outer.toml)は標準kinetic構成の小規模contract fixture、
[`examples/periodic2_zhao_transient_outer.toml`](../examples/periodic2_zhao_transient_outer.toml)はZhao過渡queueで長いflightを
物理timescale guardが拒否するexpected-fail fixtureです。
大規模productionでは、計算規模に合わせたbackendを別途選びます。

## 同じ物理効果を二重に加えない

次の組合せは、同じ成分や粒子処理を重複させます。一つの構成内では併用しません。

- `cached_kneq0`内部のsymmetric `k=0`と、最終的な場へ加えるphysical `k=0`。
- `kinetic_1d`のinterface potential mapと、有限画像 `infinity_barrier`。
- profile returnと、有限画像open `potential_barrier`のz-high処理。
- `photoelectron_density_model="kinetic_mean"`のouter densityと、surfaceへの架空の統計return deposit。

設定validationは主要なunsupported combinationをfail closedにしますが、数値収束と物理的な適用範囲はユーザーが確認します。

## outer flight中は選択したfieldを固定する

1D instant return は outer flight time を計算しますが、global simulation time には加えません。
通常の explicit 経路は batch-start field を flight 中に固定し、flight time と `field_evolution_timescale` の比が
`max_frozen_field_ratio` を超えると fail closed で停止します。`implicit_mean` 光電子は mean solve 後の $k=0$ と
batch-start $k\ne0$ を固定した準定常 shadow sampling です。flight time と ratio は同じ診断へ記録しますが、
shadow の ratio 超過では停止しません。通常の `same_batch` 粒子と ambient species の上限は変わりません。
対応する`kinetic_1d` + `zhao_charge_driven`構成では、1D eventをbatch間queueへ保存し、due batchでreturn/escapeを
計上できます。queueの外部領域は$L=10\lambda_{D,pe}$までで、そこに達した粒子はreservoirへescapeし、$L$外のRobin tailで
returnを判定しません。eventのterminal状態はenqueue時のfieldで決め、その後のfieldでは再積分しません。
各eventでは`tau_outer`、次のbatch-start pollまでの量子化遅延、midpoint crossing時刻誤差上限の合計へ
`max_frozen_field_ratio * field_evolution_timescale`の上限を課し、`batch_duration`にも同じ上限を設定時に要求します。
詳細は
[粒子のescapeとreturn](ParticleEscapeReturn.html#zhao-過渡closureでouter-flightをqueueする)を参照してください。

`implicit_mean` は UV turn-on の遅延 return current や batch 間の outer inventory を解きません。
それらが必要なら、この closure の時間発展に対応する delayed inventory / queue が別途必要です。

## componentごとに収束と収支を確認する

| 対象 | 変えるparameter | 比較量 |
| --- | --- | --- |
| near/far periodic | image layer、Ewald layer、cache cold/warm | $\phi,\mathbf E$、force、operator residual |
| physical `k=0` | 下側境界条件、height/grid refinement | Gauss residual、interface field |
| kinetic profile | Debye長、source sampling、outer stride | $\phi_I$、current、nonlinear residual |
| reservoir | macro target、batch duration | inflow current、macro residual |
| photoelectron | ray数、`dt`、outer return | emission、reabsorption、escape/return charge |
| implicit mean $k=0$ | `batch_duration`、ray数、macro粒子数、下側境界条件 | $\phi_I$、species電流、`mean_solver_iterations`、`transaction_residual_C`、`sample_escape_fraction`、`return_weight_scale` |
| Zhao過渡queue | `batch_duration`、time-scale guard、ray数、面積、interface位置 | $\eta$、column residual、return/escape current、force |
| batch coupling | `batch_duration` | steady surface chargeとcurrent balance |

`summary.txt`、`outer_plasma_profile.csv`、`charge_ledger.csv`、`charges.csv`を合わせて確認します。fieldだけ、particle countだけ、
またはsurface chargeだけの単独確認ではcomponent ownershipの二重加算やcharge lossを見落とします。

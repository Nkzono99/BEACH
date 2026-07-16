title: periodic2無限周期＋outer plasma構成

Lang: [日本語](InfinitePeriodicOuterConfiguration.md) | [English](InfinitePeriodicOuterConfiguration.en.md)

# periodic2無限周期＋outer plasma構成

この構成は、x/y無限周期のsurface fieldとz方向の外部plasmaモデルを一つの電場・電位として組み立てます。
near image、Ewaldから生成したfar `k\ne0` operator、物理`k=0`、outer responseを重複なく合成し、同じouter potentialを
reservoir inflowとparticle returnに使います。

## 各field成分に一つのownerを割り当てる

| component | owner |
| --- | --- |
| primary/near image | FMMまたはspectral base field |
| infinite-periodic far `k\ne0` | `cached_kneq0`、小規模検証では`panel_spectral_reference` |
| surface `k=0` | periodic zero-mode plan/state |
| outer mean response | `kinetic_1d`またはunified zero-mode solve |
| outer nonzero response | unified modelのreflection/transmission correction |
| reservoir velocity map | outer interface potential差 |
| outward escape/return | 同じouter profileまたは同じunified 3D field |

periodic operatorと`k=0`の式は[periodic2静電場](PeriodicElectrostatics.html)で説明します。

## nonlinear splitとlinear unifiedを選ぶ

| 構成 | mean plasma | nonzero mode | particle transfer |
| --- | --- | --- | --- |
| split kinetic | VDFに基づく密度モデルを使う非線形`kinetic_1d` | surface側でinterfaceまで減衰すると仮定 | `kinetic_1d_profile_return` |
| unified linear | accessible fraction付き線形Poisson | response startでscreened modeへ接続 | なし、または3D explicit orbit |

split kineticはspecies別VDF、Bohm entry、photoelectron mean densityを扱えます。一方で、surfaceとinterfaceの間には、
local surface fieldだけを扱うsplit windowを仮定します。[kinetic 1D外部プラズマ](KineticOuterPlasma.html)に
この分割と非線形solveをまとめています。

unified linearはroughness範囲からplasma responseを入れられますが、線形Debye応答モデルでありspecies VDFやBohm条件を
解きません。[unified linear response](UnifiedLinearResponse.html)で適用範囲とfield solveを説明します。

## batch開始時に確定した場を流入とreturnで共有する

1. commit済み`q_elem`からFMM/source multipoleとsurface zero modeをrefreshする。
2. cached far operatorをcurrent root multipoleへ適用する。
3. split構成では必要なstrideでinterface fieldからouter profileを解く。unified構成では場を更新するたびに解く。
4. potential gauge、Gauss residual、interface/linearity diagnosticsを更新する。
5. outer stateを使ってz-high reservoirのglobal粒子数とinterface速度を決める。
6. photo raycastを行い、source reaction chargeをbatch差分へ記録する。
7. 全粒子をbatch内で固定された場の中で追跡する。
8. z-high outward crossingを同じouter stateでescape/returnへ写像する。
9. surface absorptionとemission差分をMPI all-reduceし、batch末尾にcommitする。

operatorやouter Poisson problemは、particle hitごとに解き直しません。sourceとreturnは同じbatch-start outer stateを共有し、
surface chargeは次のbatchからfieldに反映されます。

## 同じenergy式をkinetic流入とreturnに使う

z-highの最初の負・正`reservoir_face` speciesを無限遠electron/ion VDFとして使います。収束profileの
$\phi_I-\phi_\infty$により

- inflowでは到達可能な$v_{n,\infty}$を選び、$v_{n,I}$へ加減速する。
- outflowでは$v_{n,I}$から$v_{n,\infty}^2$を計算し、escape/turning returnを判定する。

同じenergy equationを逆向きに使うため、`reservoir_potential_model`やZhao cutoffを重ねません。
[`reservoir_face` の流入量と速度サンプリング](ReservoirInjection.html)が流入側、[粒子のescapeとreturn](ParticleEscapeReturn.html)が
流出側の写像を説明します。

## 平均outer密度とtracked光電子を別々に選ぶ

photoelectronの扱いはouter field modelとparticle returnを分けて選びます。

| 選択 | outer density | tracked particle |
| --- | --- | --- |
| `photoelectron_closure="none"` | photoelectron mean densityなし | 通常のsource/軌道だけ |
| `kinetic_mean` + transferなし | 定常outgoing/returning mean density | z-highでは通常open処理 |
| `kinetic_mean` + profile return | 同じmean density | 個々のinterface crossingもprofileでreturn/escape |
| unified + explicit orbit | 光電子の平均密度モデルなし | 個々の3D外部軌道 |

tracked returnでは`deposit_opposite_charge_on_emit=true`を要求し、legacy `photo_escape_model`を併用しません。平均密度モデルは
tracked surface depositを置き換えず、統計的return chargeを追加depositしません。
[光電子の放出とライフサイクル](PhotoelectronEmission.html)に、放出から再吸収までの電荷収支をまとめています。

## productionではcached nonzero operatorを使う

FMM productionでは次のownershipを明示します。

```toml
[sim]
field_solver = "fmm"
field_bc_mode = "periodic2"
field_periodic_image_layers = 1
field_periodic_far_correction = "cached_kneq0"
field_periodic_ewald_layers = 4

[field]
element_kernel = "triangle_p0"

[periodic2]
nonzero_mode_backend = "cached_kneq0"
zero_mode_policy = "exclude_k0"
lower_boundary_model = "symmetric_vacuum"
```

このblockだけではouter plasmaを有効にしません。`[outer_plasma]`と`[coupling]`でkineticまたはunified構成を追加します。
最小のcache例は[`examples/periodic2_cached_panel.toml`](../examples/periodic2_cached_panel.toml)です。

小規模のreference計算にはDirect + `panel_spectral_reference`を使います。
[`examples/periodic2_kinetic_outer.toml`](../examples/periodic2_kinetic_outer.toml)と
[`examples/periodic2_unified_linear_response.toml`](../examples/periodic2_unified_linear_response.toml)は、各modelの前提を確認するための例です。
大規模productionでは、計算規模に合わせたbackendを別途選びます。

## 同じ物理効果を二重に加えない

次の組合せは、同じ成分や粒子処理を重複させます。一つの構成内では併用しません。

- `cached_kneq0`内部のsymmetric `k=0`と、最終的な場へ加えるphysical `k=0`。
- `kinetic_1d`のinterface potential mapと、有限画像 `infinity_barrier`。
- profile returnと、有限画像open `potential_barrier`のz-high処理。
- tracked photoelectron returnと`boltzmann_cutoff`。
- `kinetic_mean` outer densityと、surfaceへの架空の統計return deposit。
- unified base nonzero fieldと、reflection/transmission後のincident mode。

設定validationは主要なunsupported combinationをfail closedにしますが、数値収束と物理的な適用範囲はユーザーが確認します。

## outer flight中はbatch-start fieldを固定する

1D instant returnと3D explicit orbitはouter flight timeを計算しますが、global simulation timeには加えません。
flight timeと`field_evolution_timescale`の比は、`max_frozen_field_ratio`以下に保ちます。persistent delayed-return queueは未実装です。
steady/quasisteadyでないreturn currentには、[粒子のescapeとreturn](ParticleEscapeReturn.html)で示す制約が加わります。

## componentごとに収束と収支を確認する

| 対象 | 変えるparameter | 比較量 |
| --- | --- | --- |
| near/far periodic | image layer、Ewald layer、cache cold/warm | $\phi,\mathbf E$、force、operator residual |
| physical `k=0` | 下側境界条件、height/grid refinement | Gauss residual、interface field |
| kinetic profile | Debye長、source sampling、outer stride | $\phi_I$、current、nonlinear residual |
| unified profile | `unified_grid_points`、height sampling、mode layer | linearity、accessible fraction、Gauss residual |
| reservoir | macro target、batch duration | inflow current、macro residual |
| photoelectron | ray数、`dt`、outer return | emission、reabsorption、escape/return charge |
| batch coupling | `batch_duration` | steady surface chargeとcurrent balance |

`summary.txt`、`outer_plasma_profile.csv`、`charge_ledger.csv`、`charges.csv`を合わせて確認します。fieldだけ、particle countだけ、
またはsurface chargeだけの単独確認ではcomponent ownershipの二重加算やcharge lossを見落とします。

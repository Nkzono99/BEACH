title: 境界・外部領域の構成を選ぶ

Lang: [日本語](OuterPlasmaModels.md) | [English](OuterPlasmaModels.en.md)

# 境界・外部領域の構成を選ぶ

有限な粒子計算boxと外部reservoirをつなぐ処理は、粒子生成、流入補正、流出境界、外部場、box外粒子の
5段階に分かれます。これらは一つの`model`を選ぶ設定ではなく、必要な段階を矛盾しないように組み合わせる設定です。

外部reservoirと表面帯電を自己整合に結ぶ通常のproduction計算では、`kinetic_1d`を標準の外部シースモデルとして
使います。`unified_linear_response`はその上位版ではなく、rough surface近傍の横方向fieldを線形screeningする必要が
ある場合だけ選ぶ高度なモデルです。外部シースを使わないケースとの互換性を保つため、設定上の暗黙の既定は`none`の
ままです。

## 5つの設定段階を区別する

| 段階 | 主な設定 | 役割 |
| --- | --- | --- |
| 粒子生成 | `particles.species[].source_mode` | `reservoir_face`、`photo_raycast`などからmacro粒子を作る |
| 流入補正 | `sim.reservoir_potential_model`または`sim.sheath_injection_model` | 上流VDFの到達条件、密度、drift、cutoffを決める |
| 流出境界 | `sim.open_boundary_model` | open面で無条件escapeまたはscalar障壁による反射を選ぶ |
| 外部場 | `outer_plasma.model` | box外の1D profileまたはsurfaceからfarまでの応答場を構成する |
| box外粒子 | `coupling.particle_transfer_mode` | interfaceから出た粒子を1D写像または3D軌道へ渡す |

`reservoir_face`は粒子源です。`infinity_barrier`やシース流入補正は、その粒子源へ渡すVDFを補正します。
`potential_barrier`は逆向きにopen面へ出た粒子を処理します。`kinetic_1d`と`unified_linear_response`は、
scalar補正より多くの空間情報を持つ外部場です。

## 境界のownerを一つの実行時契約へ解決する

公開TOMLのキーは互換性のため上の5段階に分かれたままです。BEACHは流入補正、通常open面、
z-high輸送、queueのownerを一つの内部契約へ解決します。粒子loopは実行開始時に解決した契約を
read-onlyで共有し、batch injectionもhot loop外で同じresolverを使います。

1D transferを有効にした`linear_debye`と`kinetic_1d`は、returnだけでなくz-highの流入写像も同じprofileで
所有します。そのため`infinity_barrier`やlegacy `sheath_injection_model`との重ね合わせは設定時に拒否します。
`open_boundary_model`は別責務であり、outerが所有しない通常open面には引き続き適用できます。
この正規化は内部実装だけの変更で、既存TOMLを新しいbundleへ書き換える必要はありません。

## 流入側は粒子生成とVDF補正を分ける

```text
上流VDF
   │
   ├─ 補正なし
   ├─ infinity_barrier         face平均電位でenergy map
   ├─ sheath_injection_model  解析モデルで密度・drift・cutoffを補正
   └─ kinetic_1d profile      解いた外部電位でenergy map
             ↓
       reservoir_face
       fluxを積分して粒子を生成
```

`source_mode="reservoir_face"`は、与えられた上流VDFをflux-weighted分布へ変換し、box面上に粒子を生成します。
このsampling自体はシースを解きません。上図はreservoir流入の経路です。Zhao系モデルはこれに加えて、
対応する`photo_raycast` sourceの密度、cutoff、放出電流も補正します。

| 流入モデル | 入力から決めるもの | 空間profile |
| --- | --- | --- |
| 補正なし | 設定したVDFをface上の分布として使用 | なし |
| `infinity_barrier` | face平均電位と`phi_infty`から到達cutoff・face速度を計算 | なし |
| `floating_no_photo` | electron/ion流入が釣り合うelectron cutoff | なし |
| `zhao_*` | 解析branchからelectron/ion/photoelectronのVDFを補正 | 解析モデルだけ |
| `kinetic_1d` | 収束した外部Poisson profileからinterfaceへの流入を写像 | 離散1D profile |

`floating_no_photo`と`zhao_*`は`sim.sheath_injection_model`で選びます。これらはsource VDFの事前補正であり、
生成後の粒子を進める空間電場ではありません。詳しいsamplingは
[`reservoir_face` の流入量と速度サンプリング](ReservoirInjection.html)、
解析モデルは[流入VDFのシース補正](SheathInjectionClosures.html)にあります。

## 流出側はopen境界とouter transferを分ける

open面へ到達した粒子には、通常のopen境界規則または外部領域へのtransferを適用します。

| 流出モデル | 境界通過後の処理 | box外の状態 |
| --- | --- | --- |
| `open_boundary_model="escape"` | その場で粒子を除去 | なし |
| `open_boundary_model="potential_barrier"` | 通過点電位と法線運動energyから反射・escapeを判定 | scalar電位だけ |
| `electrostatic_1d_instant_return` | 1D profileからescape/turning point/往復時間を計算。対応するZhao構成では有限$10\lambda_{D,pe}$領域のeventをbatch間queueへ遅延可能 | 解析または離散1D profile、またはrank-local event queue |
| `electrostatic_3d_explicit_orbit` | batch内で固定された外部3D場中を時間積分 | zero/nonzero 3D field |

`potential_barrier`はsourceの種類に依存しません。reservoir粒子、光電子、volume seed粒子が同じ状態で同じopen面を
横切れば、同じ判定を受けます。z-highをouter ownership interfaceとしてtransferする構成では、その面の
escape/returnは`open_boundary_model`ではなく対応するouter modelが担当します。

判定式、return時間、Zhao過渡queue、3D軌道、準定常近似は
[open境界・escape・return](ParticleEscapeReturn.html)で説明します。

## 外部場の複雑さからmodelを選ぶ

まず`kinetic_1d`と対応するreturn/transfer設定で、無限遠VDF、interfaceへの流入、平均シース、escape/returnを
一つのprofileで閉じられるかを検討します。surfaceのroughness範囲とplasma responseが重なり、interfaceまでに横方向modeが減衰するsplit windowを
置けず、かつ線形応答で十分な場合に限って`unified_linear_response`を選びます。

| `outer_plasma.model` | 位置付け | 解くもの | 適した用途 |
| --- | --- | --- | --- |
| `none` | 外部シースなし | 外部場なし | 単純open境界、scalar barrier、解析的な注入補正 |
| `linear_debye` | 簡易・参照 | interfaceからの指数zero-mode応答 | 簡易1D instant return |
| `kinetic_1d` | **標準・推奨** | VDFに基づく密度モデルを含む非線形1D Poisson profile | 対応するtransferと組み合わせた自己整合な流入・電流・escape/return |
| `unified_linear_response` | 高度・限定用途 | rough surfaceからfarまでの線形zero/nonzero screening | 真空split windowを置けないrough surface |

`kinetic_1d`はambient electron/ionの無限遠VDF、Bohm条件、Poisson residual、無限遠準中性を満たす単調分枝を
解きます。`unified_linear_response`はspecies別VDFや電流balanceを解かず、Debye–Hückel型の線形応答と
plasma-accessible areaを場へ組み込みます。

標準構成は[外部シース: kinetic 1D](KineticOuterPlasma.html)にまとめています。roughnessとplasma応答が同じ領域に
重なるが線形近似でよい場合だけ、[高度な粗面線形screening](UnifiedLinearResponse.html)を使います。

## 典型的な構成を選ぶ

| 目的 | 流入 | 流出 | 外部場・transfer |
| --- | --- | --- | --- |
| 単純な有限box | `reservoir_face`、補正なし | `escape` | なし |
| 有限画像＋scalar barrier | `infinity_barrier` | `potential_barrier` | なし |
| Zhao文献モデル | `sheath_injection_model="zhao_*"` | 通常のopen境界 | なし |
| 自己整合1D sheath | kinetic profileから流入写像 | kinetic profile return | `kinetic_1d` + 1D instant transfer |
| 強UV立上がり | kinetic Zhao profileから流入写像 | due batchでprofile return/escape | `kinetic_1d` + `zhao_charge_driven` + transient queue |
| rough surface線形応答 | faceで指定したsource | openまたは3D outer orbit | `unified_linear_response` |

完成した設定例は[有限画像構成](FinitePeriodicConfiguration.html)と
[無限周期＋outer plasma構成](InfinitePeriodicOuterConfiguration.html)にあります。
Zhao過渡queueのfrozen-field guard fixtureは
[`periodic2_zhao_transient_outer.toml`](../examples/periodic2_zhao_transient_outer.toml)です。これは記載した物理timescaleで
長いflightをfail-closedに停止する負例であり、成功した物理検証RUNではありません。
このqueueでは$L=10\lambda_{D,pe}$到達をreservoir escapeとし、
各eventの`tau_outer`、次のbatch-start pollまでの遅延、midpoint crossing時刻誤差上限の合計へ
`max_frozen_field_ratio * field_evolution_timescale`の上限を課し、`batch_duration`にも同じ上限を要求します。

## 同じ役割の補正を重ねない

- `sheath_injection_model`と`reservoir_potential_model`は同時に使わない。
- `kinetic_1d`のprofile returnでは流入にも同じprofileを使い、legacy `sheath_injection_model="zhao_*"`や
  `infinity_barrier`を重ねない。Zhao populationを使う場合は`kinetic_closure="zhao_charge_driven"`を選ぶ。
- `unified_linear_response`はsource VDFや浮遊電流を決めない。reservoirを使う場合はface上の分布を別に定義する。
- outer transfer対象のz-high面と、通常の`potential_barrier`によるscalar反射を同じ処理として解釈しない。
- modelが失敗したときに、別の簡略化モデルへsilent fallbackしない。

profile statusだけでなく、Poisson residual、境界条件、電荷積分、current balance、grid refinement、
frozen-field ratioを、使用したmodelに応じて確認します。共通の手順は
[計算結果の妥当性確認](ValidationGuide.html)にまとめています。

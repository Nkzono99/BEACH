title: 外部境界を構成する

Lang: [日本語](OuterPlasmaModels.md) | [English](OuterPlasmaModels.en.md)

# 外部境界を構成する

このページでは、外部 plasma 応答の場、z-high を横切る粒子、その他の open 面を設定します。
通常の入力では `[external_boundary]` 配下だけを編集します。BEACH が実装用の return・transfer・queue
設定へ変換するため、それらの内部名を組み合わせる必要はありません。

## 最初に 3 つの責務を選ぶ

```text
[external_boundary]
├── [external_boundary.field]          外部 plasma 応答の場
├── [external_boundary.particles]      z-high の粒子と reservoir 流入
└── [external_boundary.ordinary_open]  outer が所有しない open 面
```

`field` と `particles` は必須です。`ordinary_open` は任意で、省略時は `escape` です。

標準の自己整合 1D 外部シースは次の構成です。

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "absorbing_maxwellian"
debye_length = 0.2
thermal_voltage = 2.0

[external_boundary.particles]
mode = "same_batch"
field_evolution_timescale = 1.0
```

この入力から、BEACH は kinetic profile による流入写像、1D return、同一 batch 内の復帰を選びます。
`interface_z` は `sim.box_max[2]`、更新方式は `explicit` に決まるため入力しません。
通常 open 面を単純に escape させる場合、`ordinary_open` は上のように省略します。

すべての open 面を単純に escape させ、流入補正や外部場も使わない場合は `[external_boundary]` 自体を省略します。
それ以外は、目的に最も近い行から構成を選んでください。

| 目的 | `field` | `particles.mode` | `inflow_model` |
| --- | --- | --- | --- |
| 外部場なし、または scalar barrier | `none` | `local_source` | `source_vdf` / `infinity_barrier` |
| 標準の自己整合 1D シース | `kinetic_1d` + `absorbing_maxwellian` | `same_batch` | `auto` |
| 蓄積電荷で閉じる Zhao 定常・準定常シース | `kinetic_1d` + `zhao_charge_driven` | `same_batch` | `auto` |
| Zhao の外部飛行遅延を含む過渡 | `kinetic_1d` + `zhao_charge_driven` | `zhao_queue` | `auto` |

通常流出に scalar barrier が必要な場合だけ次を追加します。

```toml
[external_boundary.ordinary_open]
model = "potential_barrier"
```

## 外部場を選ぶ

`external_boundary.field.model` は外部 plasma 応答の場だけを選びます。`kinetic_1d` は
interface 外を解きます。粒子を外部領域へ渡すかどうかは
`external_boundary.particles.mode` で別に選びます。

| `field.model` | 位置付け | 主な用途 |
| --- | --- | --- |
| `none` | 外部場なし | 通常 open、設定済み source VDF、scalar barrier |
| `kinetic_1d` | 対応する自己整合 1D シース | reservoir VDF、平均シース、流入、return を同じ profile で閉じる |

field-only も有効な構成です。

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "absorbing_maxwellian"
debye_length = 0.2
thermal_voltage = 10.0

[external_boundary.particles]
mode = "local_source"
inflow_model = "source_vdf"
```

## 接続面の粒子を選ぶ

| `particles.mode` | z-high を出た粒子 | 用途 |
| --- | --- | --- |
| `local_source` | `ordinary_open` で処理し、外部軌道を保持しない | 外部場なし、field-only、scalar 流入 |
| `same_batch` | kinetic 1D profile による return を計算 | 定常・準定常の tracked return |
| `zhao_queue` | Zhao 光電子の return/escape event を due batch まで保持 | 強 UV 立上がりなどの遅延電流 |

`particles.mode` が選ぶのは z-high 粒子の外部輸送だけです。reservoir からの流入 VDF は
後述の `inflow_model` で独立に選びます。

`same_batch` は `field.model="kinetic_1d"` と組み合わせ、離散 kinetic 1D profile 上で
return / escape を判定します。

`zhao_queue` は汎用的な遅延輸送ではありません。`kinetic_1d` と
`kinetic_closure="zhao_charge_driven"` の組合せだけで使用できます。

## reservoir 流入を選ぶ

`external_boundary.particles.inflow_model` は `reservoir_face` に渡す上流 VDF の扱いです。

| `inflow_model` | 動作 |
| --- | --- |
| `auto` | 1D tracked model では field profile、それ以外では設定した source VDF を使用 |
| `source_vdf` | 設定した VDF を face 上の分布として使用 |
| `infinity_barrier` | face 平均電位と `sim.phi_infty` から到達速度を補正 |

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "zhao_charge_driven"
# split-interface 診断の基準量。Zhao profile の scale ではない
debye_length = 10.5
thermal_voltage = 10.0

[external_boundary.particles]
mode = "same_batch" # 過渡 queue なら "zhao_queue"
field_evolution_timescale = 1.0
```

Zhao profile の物理的な長さ・電位 scale は、光電子ありでは $T_{pe}$ と基準密度、
光電子なしでは ambient $T_e$ と $n_\infty$ から導出します。`debye_length` と
`thermal_voltage` は Zhao root/profile を変えませんが、split-interface の gap、
lateral field、local-charge 診断の基準量として現在は必須です。
`mode="zhao_queue"` へ切り替える場合は、正の `sim.batch_duration` と光電子 source も必要です。
update stride は queue が内部で 1 に固定するため追加しません。

`kinetic_1d` の `same_batch` と `zhao_queue` では、同じ 1D profile が流入を所有します。
この場合、`inflow_model` は `auto` だけを許し、`infinity_barrier` との二重補正を拒否します。

## 通常 open 面を選ぶ

`external_boundary.ordinary_open.model` は outer model が所有しない open 面へ適用します。

| `ordinary_open.model` | 動作 |
| --- | --- |
| `escape` | 粒子を永久流出として除去 |
| `potential_barrier` | 局所電位と法線運動 energy から反射または escape を判定 |

`potential_barrier` は reservoir 流入の `infinity_barrier` と別の責務です。一方だけでも、両方でも使えます。
tracked outer model が z-high を所有する場合も、`ordinary_open` はその他の open 面に残ります。

## 設定エラーとして扱う組合せ

BEACH は矛盾する値を補正したり、別モデルへ silent fallback したりしません。

- `field.model="none"` と `same_batch` / `zhao_queue`
- kinetic の tracked 1D profile と local inflow 補正
- Zhao 以外の closure と `zhao_queue`
- `zhao_queue` と手動 branch、steady start
- 選んだ field / particle mode で効果のない key。既定値を明示しただけでも no-op として拒否
- 対応しない磁場、species、periodic2、zero-mode、時間スケール
- `[external_boundary]` と旧 `[outer_plasma]` / `[coupling]` の混在

物理入力や数値 guard は自動推定しません。species、Debye 長、温度、field timescale、
periodic2 backend などは使用する model に応じて明示し、矛盾時はエラーにします。
`zhao_queue` の update stride は内部で 1 に固定するため、公開入力へ書きません。

## 解決された契約を確認する

`summary.txt` には authoring 入力ではなく、実行に使った解決済みの流入、通常 open 面、interface transport、
粒子 mode を記録します。

旧 `[sim]` / `[outer_plasma]` / `[coupling]` からの移行は
[外部境界設定の移行](BoundaryConfigurationMigration.html)を使ってください。個々のモデルの物理と妥当性確認は
[外部シース: kinetic 1D](KineticOuterPlasma.html)、
[open 境界・escape・return](ParticleEscapeReturn.html)に分けています。

title: 外部境界設定を移行する

Lang: [日本語](BoundaryConfigurationMigration.md) | [English](BoundaryConfigurationMigration.en.md)

# 外部境界設定を移行する

このページは、旧 `[sim]` / `[outer_plasma]` / `[coupling]` の外部境界キーを
`[external_boundary]` authoring facade へ移すための対応表です。旧構文は互換入力として読み込めますが、
新しい example とガイドは facade を正本とします。

## 移行ルール

- `[external_boundary]` と旧 `[outer_plasma]` / `[coupling]` は同じファイルに書かない。
- facade 使用時は `sim.reservoir_potential_model`、`sim.sheath_injection_model`、
  `sim.open_boundary_model` を書かない。
- `sim.phi_infty`、legacy sheath の物理パラメータ、box、species、periodic2、solver は引き続き併用できる。
- 旧構文と新構文は同じ runtime contract へ解決され、物理 fingerprint に authoring 形式の違いを入れない。

同じ値を両側へ書いても merge しません。値の一致を前提にした last-wins は設定ミスを隠すため、混在はエラーです。

## キーの対応

| 旧設定 | 新設定 |
| --- | --- |
| `outer_plasma.model` | `external_boundary.field.model` |
| `outer_plasma.kinetic_closure` | `external_boundary.field.kinetic_closure` |
| その他の outer 物理・診断キー | `external_boundary.field` の同名キー |
| `coupling` の時間スケール・orbit・steady-start キー | `external_boundary.particles` の同名キー |
| `sim.open_boundary_model` | `external_boundary.ordinary_open.model` |
| `sim.reservoir_potential_model="infinity_barrier"` | `particles.inflow_model="infinity_barrier"` |
| `sim.sheath_injection_model=...` | `inflow_model="legacy_sheath"` と `legacy_sheath_model=...` |
| `outer_plasma.return_model` | 削除。field と particle mode から導出 |
| `coupling.particle_transfer_mode` | 削除。`particles.mode` から導出 |
| `coupling.outer_queue_enabled` | 削除。`particles.mode="zhao_queue"` で導出 |
| `coupling.update_mode` | 削除。現行は `explicit` 固定 |
| `outer_plasma.interface_z` | 削除。`sim.box_max[2]` から導出 |

## scalar 流入と通常 open 面

旧設定:

```toml
[sim]
reservoir_potential_model = "infinity_barrier"
open_boundary_model = "potential_barrier"
phi_infty = 0.0
```

新設定:

```toml
[sim]
phi_infty = 0.0

[external_boundary.field]
model = "none"

[external_boundary.particles]
mode = "local_source"
inflow_model = "infinity_barrier"

[external_boundary.ordinary_open]
model = "potential_barrier"
```

`infinity_barrier` は流入、`potential_barrier` は通常 open 面を担当します。一方だけを有効にする構成も維持されます。

## legacy Zhao 流入補正

旧設定:

```toml
[sim]
sheath_injection_model = "zhao_auto"
open_boundary_model = "escape"
sheath_alpha_deg = 60.0
sheath_photoelectron_ref_density_cm3 = 64.0
```

新設定:

```toml
[sim]
sheath_alpha_deg = 60.0
sheath_photoelectron_ref_density_cm3 = 64.0

[external_boundary.field]
model = "none"

[external_boundary.particles]
mode = "local_source"
inflow_model = "legacy_sheath"
legacy_sheath_model = "zhao_auto"

[external_boundary.ordinary_open]
model = "escape"
```

これは静的 source-VDF 補正です。`kinetic_closure="zhao_charge_driven"` へ自動変換しません。

## kinetic 1D tracked return

旧設定:

```toml
[outer_plasma]
model = "kinetic_1d"
kinetic_closure = "absorbing_maxwellian"
return_model = "kinetic_1d_profile_return"
interface_z = 1.0
debye_length = 0.2
thermal_voltage = 2.0

[coupling]
update_mode = "explicit"
particle_transfer_mode = "electrostatic_1d_instant_return"
field_evolution_timescale = 1.0
max_frozen_field_ratio = 0.1
outer_queue_enabled = false
```

新設定:

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "absorbing_maxwellian"
debye_length = 0.2
thermal_voltage = 2.0

[external_boundary.particles]
mode = "same_batch"
inflow_model = "auto"
field_evolution_timescale = 1.0
max_frozen_field_ratio = 0.1

[external_boundary.ordinary_open]
model = "escape"
```

linear field-only は `field.model="linear_debye"` と `particles.mode="local_source"`、
linear 1D return は `particles.mode="same_batch"` にします。`kinetic_1d` の field-only も同じ区別です。

## Zhao queue

旧 `outer_queue_enabled=true` は次へ移します。

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "zhao_charge_driven"
debye_length = 10.5
thermal_voltage = 10.0

[external_boundary.particles]
mode = "zhao_queue"
inflow_model = "auto"
field_evolution_timescale = 2.0e-5
max_frozen_field_ratio = 0.2
```

queue は `zhao_branch="auto"`、更新 stride 1、非 histogram を要求します。これらの固定値を繰り返し入力せず、
矛盾する明示値はエラーにします。

## unified 3D orbit

旧 explicit 3D return/transfer の対は、次の `same_batch` へ移します。

```toml
[external_boundary.field]
model = "unified_linear_response"
debye_length = 0.2
thermal_voltage = 10.0
unified_grid_points = 129

[external_boundary.particles]
mode = "same_batch"
inflow_model = "source_vdf"
field_evolution_timescale = 1.0e-4
max_frozen_field_ratio = 0.1
outer_orbit_dt = 1.0e-9
outer_orbit_max_steps = 10000
outer_orbit_energy_tolerance = 1.0e-3
```

unified 3D orbit は流入を所有しないため、`source_vdf`、`infinity_barrier`、legacy sheath を別に選べます。
field-only は `particles.mode="local_source"` を使います。

移行後は `beachx lint path/to/beach.toml` を実行し、実行後の `summary.txt` で解決された
inflow、ordinary open、interface transport、particle mode を確認してください。

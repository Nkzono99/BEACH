title: 外部境界設定を移行する

Lang: [日本語](BoundaryConfigurationMigration.md) | [English](BoundaryConfigurationMigration.en.md)

# 外部境界設定を移行する

このページは、旧 `[outer_plasma]` / `[coupling]` のうち現行モデルに対応する設定を
`[external_boundary]` authoring facade へ移す方法と、削除された外部境界モデルを含む入力の扱いを示します。
新しい入力では `[external_boundary]` を正本としてください。

## 現行の選択肢

現行の公開設定は次に限定されています。

| 責務 | 選択肢 |
| --- | --- |
| `external_boundary.field.model` | `none` / `kinetic_1d` |
| `external_boundary.particles.mode` | `local_source` / `same_batch` / `zhao_queue` |
| `external_boundary.particles.inflow_model` | `auto` / `source_vdf` / `infinity_barrier` |
| `external_boundary.ordinary_open.model` | `escape` / `potential_barrier` |

`potential_barrier` は通常 open 面からの流出、`infinity_barrier` は reservoir からの流入を担当します。
両者は独立して選べます。`kinetic_1d` の tracked mode では同じ kinetic profile が流入と return を所有するため、
`inflow_model="auto"` を使います。

## 削除された設定

解析的な線形 Debye 外部場、`unified_linear_response`、静的な sheath source 補正、
その線形モデル専用の光電子 histogram は削除されました。
該当する facade 値と旧 runtime selector はエラーになり、別モデルへ自動変換されません。

| 削除された入力・artifact | 現行の扱い |
| --- | --- |
| `external_boundary.field.model="linear_debye"` / `outer_plasma.model="linear_debye"` | 設定エラー。別の field model への alias はない |
| `external_boundary.field.model="unified_linear_response"` / `outer_plasma.model="unified_linear_response"` | 設定エラー。rough surface 3D screening に等価な現行モデルはない |
| `unified_grid_points`、`accessible_fraction_tolerance`、`max_linearity_ratio` | 設定エラー。unified field 専用の離散化・適用性設定は削除 |
| `electrostatic_3d_explicit_orbit` と `outer_orbit_*` | 設定エラー。専用 3D outer transfer とその積分設定は削除 |
| `external_boundary.particles.inflow_model="legacy_sheath"`、`legacy_sheath_model`、`sim.sheath_injection_model` | 設定エラー。静的 source 補正は実行しない |
| `external_boundary.field.infinity_potential` / `outer_plasma.infinity_potential` | 設定エラー。kinetic の無限遠 gauge は内部で 0 に固定 |
| `photoelectron_histogram_*` と `photoelectron_histogram.csv` | 設定、出力、checkpoint state を削除 |

移行先は「等価な置換」ではなく、計算目的に応じて選ぶ新しい物理モデルです。

| 以前の目的 | 検討する現行構成 | 重要な違い |
| --- | --- | --- |
| 簡易な 1D field / return | `kinetic_1d` + `local_source` または `same_batch` | species VDF と非線形 Poisson profile を解くため、線形解析解とは一致しない |
| rough surface を含む 3D screening | 直接の置換なし | 旧モデルが必要だったケースは、物理契約と検証問題から再設計する |
| source VDF を補正しない流入 | `field.model="none"` + `inflow_model="source_vdf"` | 設定した VDF を face 分布としてそのまま使う |
| 明示した無限遠電位による流入障壁 | `field.model="none"` + `inflow_model="infinity_barrier"` | 静的な電流釣合い closure ではなく scalar energy barrier |
| 蓄積電荷で閉じる Zhao sheath | `kinetic_1d` + `kinetic_closure="zhao_charge_driven"` | surface charge から profile を batch ごとに更新し、流入と return を同じ profile で閉じる |
| 線形モデル専用の光電子 histogram | 直接の置換なし | 必要な分布は出力粒子・event から別途解析する |

BEACH 1.5 / 1.6 の設定で削除値を使っている場合は、現行 parser でエラーになります。
削除されたモデル、または非デフォルト値の削除機能を有効にして作成した checkpoint は、現行 runtime contract と
model fingerprint が異なるため再開できません。`unified_linear_response` の checkpoint に移行経路はありません。
削除機能が旧デフォルト値のままだった `kinetic_1d` checkpoint-v4 は、設定から廃止キーを除けば fingerprint 互換です。
新しい構成は初期状態から実行してください。
[ADR 0010](adr/0010-remove-unified-linear-response.md)に削除理由と再設計条件を記録しています。

## 現行 raw 設定を facade へ移す

- `[external_boundary]` と旧 `[outer_plasma]` / `[coupling]` は同じファイルに書きません。
- `sim.open_boundary_model` と `sim.reservoir_potential_model` は facade 側へ移します。
- box、species、periodic2、solver、`sim.phi_infty` などの物理・数値設定はそのまま併用できます。
- 同じ値を両構文へ書いても merge しません。混在は設定エラーです。

| 旧設定 | facade |
| --- | --- |
| `outer_plasma.model` | `external_boundary.field.model` |
| `outer_plasma.kinetic_closure` | `external_boundary.field.kinetic_closure` |
| その他の現行 outer 物理・診断キー | `external_boundary.field` の同名キー |
| `coupling` の時間スケール・steady-start キー | `external_boundary.particles` の同名キー |
| `sim.open_boundary_model` | `external_boundary.ordinary_open.model` |
| `sim.reservoir_potential_model="infinity_barrier"` | `external_boundary.particles.inflow_model="infinity_barrier"` |
| `outer_plasma.return_model` | 記述しない。field と particle mode から導出 |
| `coupling.particle_transfer_mode` | 記述しない。`particles.mode` から導出 |
| `coupling.outer_queue_enabled` | 記述しない。`particles.mode="zhao_queue"` で導出 |
| `coupling.update_mode` | 記述しない。通常は `explicit`、`ambient_linear_debye + same_batch + photo_raycast` は内部の `implicit_mean` を自動導出 |
| `outer_plasma.interface_z` | 記述しない。`sim.box_max[2]` から導出 |

## scalar 流入と通常 open 面

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

`infinity_barrier` は流入、`potential_barrier` は通常 open 面を担当します。一方だけを有効にしても構いません。

## kinetic 1D tracked return

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

field-only なら `particles.mode="local_source"`、kinetic profile で return / escape まで閉じるなら
`particles.mode="same_batch"` を使います。

## Zhao queue

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

queue は `zhao_branch="auto"`、更新 stride 1 を要求します。固定値を繰り返し入力せず、
矛盾する明示値はエラーにします。

移行後は `beachx lint path/to/beach.toml` を実行し、実行後の `summary.txt` で解決された
inflow、ordinary open、interface transport、particle mode を確認してください。

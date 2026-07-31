title: 粒子源と境界流入

Lang: [日本語](ParticleSourcesBoundaries.md) | [English](ParticleSourcesBoundaries.en.md)

# 粒子源と境界流入

`source_mode` はシミュレーション領域内で粒子を生成する方法を選び、
`[particles.species.boundary_inflow]` は領域外から境界を横切る流入を加えます。2 つは異なる責務です。

## 生成位置に応じて設定を選ぶ

| 設定 | 粒子数を決める量 | 生成位置 | 適した用途 |
| --- | --- | --- | --- |
| `source_mode="volume_seed"` | `npcls_per_step` | `pos_low` から `pos_high` の体積 | 初期粒子、軌道試験、指定個数の生成 |
| `source_mode="plane_source"` | 流入流束、矩形面積、`batch_duration` | ボックス内部の axis-aligned 矩形面 | 内部に明示した一方向の面 source |
| `source_mode="photo_raycast"` | 電流密度、投影面積、`batch_duration`、光線数 | 光線が最初に命中した表面 | 光照射による表面放出 |
| `[particles.species.boundary_inflow]` | reservoir 流束、ボックス面積、`batch_duration` | 非周期のシミュレーション境界 | 外部 plasma reservoir からの連続流入 |

`boundary_inflow` は `source_mode` ではありません。純粋な境界流入 species は
`source_mode="volume_seed"` と `npcls_per_step=0` を使い、同じ species に初期粒子も必要なら
Maxwell 分布に限り `npcls_per_step>0` にします。速度 grid の境界流入では正の `npcls_per_step` を指定できません。
初版では、流束駆動の `plane_source` または旧 `reservoir_face` と
`boundary_inflow` を同じ species で併用できません。

## `volume_seed` で指定個数の粒子を作る

`volume_seed` は、各バッチに `npcls_per_step` 個の粒子を生成します。位置は直方体
`[pos_low, pos_high]` 内の一様分布、速度はドリフト付き Maxwell 分布です。
`thermal_speed` を指定した場合は、温度から求めた熱速度より優先します。

この source は物理的な面流束から粒子数を計算しません。`npcls_per_step=0` は、同じ species の
`boundary_inflow` だけを使う場合に有効です。正の値は Maxwell 境界流入と併用できますが、速度 grid とは併用できません。

## `plane_source` で内部矩形面から放出する

`plane_source` は、`pos_low` と `pos_high` が定めるボックス内部の axis-aligned 矩形面から
`source_normal` 方向へ粒子を生成します。

```toml
[[particles.species]]
source_mode = "plane_source"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.1093837139e-31
target_macro_particles_per_batch = 300
pos_low = [0.2, 0.2, 2.0]
pos_high = [0.8, 0.8, 2.0]
source_normal = [0.0, 0.0, -1.0]
```

`pos_low` と `pos_high` はちょうど 1 軸で同じ値を持ち、残る 2 軸は正の長さを持つ必要があります。
法線座標はボックス境界より厳密に内側、接線範囲はボックス内に置きます。`source_normal` は
zero-thickness 軸に沿う非ゼロベクトルです。
設定例では正負の単位ベクトルを推奨し、実装は `[2,0,0]` のような大きさを内部で正規化します。

粒子数は Maxwell reservoir の density・temperature、または速度グリッドの指定流束に、矩形面積と
`batch_duration` を掛けて決めます。法線速度は `source_normal` 方向の片側流束分布から生成します。
`[reservoir]` の `infinity_barrier`、`phi_infty`、`face_potential_grid_n` は内部平面には適用しません。

## `boundary_inflow` で外部 reservoir から流入させる

`[particles.species.boundary_inflow]` の 6 面キーへ `"reservoir"` を指定します。

```toml
[[particles.species]]
species_key = "solar_wind_electron"
source_mode = "volume_seed"
npcls_per_step = 0
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.1093837139e-31
target_macro_particles_per_batch = 300
drift_velocity = [0.0, 0.0, -4.0e5]

[particles.species.boundary_inflow]
z_high = "reservoir"
```

選択したボックス面全体から内向きに生成します。複数の非周期面を同時に指定できますが、
周期面は指定できません。外向き粒子の作用は `[particle_boundary]` または
`[particles.species.boundary]` で独立に指定し、流入面の有効作用は `open` にします。

流束、マクロ粒子数、面ごとの端数、`infinity_barrier` は
[シミュレーション境界からの reservoir 流入](ReservoirInjection.html)にまとめています。
このモデルはボックス外の軌道や自己無撞着 sheath を解きません。

## `photo_raycast` で照射面から粒子を放出する

`photo_raycast` は、ボックス面上の照射開口から光線を発射し、ボックス境界条件に従って進めます。
ボックス内で最初に命中した要素から粒子を放出し、要素法線に対する流束重み付き Maxwell 分布から速度を生成します。

放出元へ置く逆符号電荷、再吸収、closed PE は
[光電子の放出とライフサイクル](PhotoelectronEmission.html)で説明します。

## 旧 `reservoir_face` は互換入力

`source_mode="reservoir_face"` は既存ケースの挙動を維持する deprecated 入力です。`inject_face` と
`pos_low` / `pos_high` でボックス面上の開口を定め、旧方式の流入を生成します。BEACH はこれを
`boundary_inflow` や `plane_source` へ暗黙変換しません。

新しい外部 plasma ケースには `boundary_inflow`、内部の明示的な矩形面には `plane_source` を使ってください。
旧開口を境界全面の流入へ移す場合は、面積の変化により物理流入量が変わる点を確認します。

## 生成後は同じ粒子追跡へ入る

生成方法が異なっても、粒子は同じ追跡処理へ入ります。

| バッチ内の結果 | 処理 |
| --- | --- |
| メッシュへ吸収 | 命中要素へマクロ粒子電荷を堆積 |
| open 面で escape | 粒子を除去し、species 別の escape 量へ計上 |
| `reflect` / `redistributed_reflect` / periodic 面 | 同じ粒子で残り step を再積分 |
| `max_step` まで生存 | 未解決粒子としてバッチ末尾で破棄・計上 |

非周期面の global 作用は `[particle_boundary]`、species override は `[particles.species.boundary]`、
周期軸は `domain.periodic_axes` で指定します。粒子前進は [Boris 粒子更新](BorisPusher.html)、
mesh と box の event 順序は[粒子の衝突・境界イベント](ParticleEvents.html)、
open 面処理は[粒子の escape と局所 return](ParticleEscapeReturn.html)で説明します。

## MPI と再開で流入量を保つ

reservoir 流入と `plane_source` の全 rank 合計粒子数と端数は root rank で決め、各 rank へ分配します。
`photo_raycast` の `rays_per_batch` も全 rank の合計です。期待される流入量や放出量は MPI rank 数に
依存しませんが、乱数列と個々の粒子軌道は rank 数によって変わる場合があります。

境界流入の端数は species と face の組ごとに checkpoint へ保存し、再開時に復元します。端数処理と
確認項目は [reservoir 流入](ReservoirInjection.html)を参照してください。

## Code reference

- 粒子分布と光線追跡: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- source と境界流入のバッチ生成: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- 入力の検証: [`bem_app_config_parser_validate.f90`](../src/config/app_config_parser/bem_app_config_parser_validate.f90)
- 電荷収支とバッチ追跡: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- マクロ粒子端数の checkpoint: [`bem_restart.f90`](../src/runtime/bem_restart.f90)

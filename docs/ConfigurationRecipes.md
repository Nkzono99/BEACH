title: シミュレーションケースを設計する

Lang: [日本語](ConfigurationRecipes.md) | [English](ConfigurationRecipes.en.md)

# シミュレーションケースを設計する

公式入門ケースを基に、メッシュ、粒子源、境界条件を1項目ずつ置き換えるためのタスクガイドです。
全キーの定義は[入力パラメータリファレンス](Parameters.html)、設定ファイルの生成と検査は
[`beach.toml`を作成・検証する](Configuration.html)を参照してください。

## 共通手順

**前提:** BEACHをインストールし、作業ディレクトリを用意します。

**操作:**

```bash
beachx config init beach.toml
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

**期待する出力:** `lint`が設定を受理し、実行後に`outputs/latest/summary.txt`と`charges.csv`が作成されます。

**解釈:** 正常終了は設定と実行経路が動作したことを示します。数値収束や物理的妥当性は
[結果を検証する](ValidationGuide.html)に従って別途確認します。

**次の選択:**

| 目的 | 変更する場所 |
| --- | --- |
| 組み込み形状を使う | `[mesh]`, `[[mesh.templates]]` |
| OBJ形状を使う | `[mesh]` |
| 初期粒子、reservoir流入、光電子を選ぶ | `[[particles.species]]` |
| box、場境界、粒子境界を選ぶ | `[domain]`, `[field_boundary]`, `[particle_boundary]`, `[reservoir]` |
| 時系列を保存する | `[output]` |
| checkpointから続ける | `[output]` |

変更は1種類ずつ行い、そのたびに`beachx lint beach.toml`を実行してください。

## 組み込みメッシュを追加する

**前提:** `[mesh].mode="template"`を使います。各`[[mesh.templates]]`には別の`mesh_id`が割り当てられます。

**操作:** 既存のtemplateの後ろに形状を追加します。

| `kind` | 形状 | 寸法 | 解像度 |
| --- | --- | --- | --- |
| `plane` | XY長方形平面 | `size_x`, `size_y` | `nx`, `ny` |
| `plate_hole`, `plane_hole` | 円形穴付き平面 | `size_x`, `size_y`, `radius` | `n_theta`, `n_r` |
| `disk` | 円板 | `radius` | `n_theta`, `n_r` |
| `annulus` | 同心リング | `radius`, `inner_radius` | `n_theta`, `n_r` |
| `box` | 閉じた直方体表面 | `size` | `nx`, `ny`, `nz` |
| `cylinder` | z軸方向の円柱 | `radius`, `height`, `cap` | `n_theta`, `n_z` |
| `sphere` | 球面 | `radius` | `n_lon`, `n_lat` |

```toml
[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
enabled = true
surface_model = "insulator"
size_x = 1.0
size_y = 1.0
nx = 20
ny = 20
center = [0.5, 0.5, 0.02]

[[mesh.templates]]
kind = "sphere"
enabled = true
surface_model = "insulator"
center = [0.5, 0.5, 0.55]
radius = 0.12
n_lon = 24
n_lat = 12
```

**期待する出力:** `mesh_triangles.csv`に両方の`mesh_id`が現れ、`mesh_sources.csv`で元のtemplateを確認できます。

**解釈:** 分割数を増やすと場計算と衝突判定のコストも増えます。まず粗いメッシュで配置と衝突を確認します。
通常の帯電計算は`surface_model="insulator"`を使います。`conductor`は`field_boundary.mode="free"`だけに対応し、
`dielectric`の`epsilon_r`は現状ではメタデータであり、分極を解きません。

**次の選択:** 形状と要素数の制約は
[組み込み形状リファレンス](Parameters.html#meshtemplates-組み込み形状)にまとめています。

## OBJメッシュへ置き換える

**前提:** OBJファイルは三角形面を含み、実行時に読み取れる場所へ置きます。

**操作:**

```toml
[mesh]
mode = "obj"
obj_path = "mesh/object.obj"
surface_model = "insulator"
obj_scale = 1.0
obj_rotation = [0.0, 0.0, 0.0]
obj_offset = [0.0, 0.0, 0.0]
```

**期待する出力:** OBJ全体が1つの`mesh_id`として`mesh_triangles.csv`へ出力されます。

**解釈:** 変換順はscale → rotation → offsetです。

**次の選択:** 独立したobjectを別の`mesh_id`で解析する場合は、複数の組み込みtemplateを優先します。

## 粒子源を選ぶ

**前提:** `[[particles.species]]`の`source_mode`を1つ選びます。modeを変更するときは、別mode専用のキーを
残さずentry全体を置き換えます。

| `source_mode` | 用途 | 必須となる主なキー |
| --- | --- | --- |
| `volume_seed` | 箱内に初期粒子を置く小テスト | `npcls_per_step`, `pos_low`, `pos_high` |
| `plane_source` | 内部矩形面から一方向fluxを与える | `pos_low`, `pos_high`, `source_normal`, densityまたはgrid flux |
| `reservoir_face`（非推奨） | 旧形式の局所reservoir面 | `inject_face`, `pos_low`, `pos_high` |
| `photo_raycast` | 表面から光電子を放出する | `rays_per_batch`, `emit_current_density_a_m2`, `ray_direction` |

外部reservoirからの流入は`source_mode`ではなく`[particles.species.boundary_inflow]`で追加します。

### `volume_seed`

```toml
[[particles.species]]
source_mode = "volume_seed"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
npcls_per_step = 100
w_particle = 1.0
pos_low = [0.1, 0.1, 0.6]
pos_high = [0.9, 0.9, 0.9]
drift_velocity = [0.0, 0.0, -1.0e6]
temperature_k = 0.0
```

**期待する出力:** 各batchに`npcls_per_step`個のmacro粒子を生成します。

**解釈:** 物理的な面流束から個数を決めるsourceではありません。

### 境界reservoir流入

**前提:** `[domain]`と正の`sim.batch_duration`が必要です。流入障壁は`[reservoir]`で選びます。

```toml
[reservoir]
inflow_model = "source_vdf"
phi_infty = 0.0
face_potential_grid_n = 3
```

```toml
[[particles.species]]
source_mode = "volume_seed"
npcls_per_step = 0
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
target_macro_particles_per_batch = 300
drift_velocity = [0.0, 0.0, -4.0e5]

[particles.species.boundary_inflow]
z_high = "reservoir"
```

**期待する出力:** 外部reservoirの流束から粒子重みが決まり、z-high面全体から流入します。

**解釈:** `target_macro_particles_per_batch`は1 batchあたりの計算粒子数を固定します。重みを直接指定する場合は
代わりに`w_particle`を使います。`temperature_k`と`temperature_ev`は同時指定できません。

**次の選択:** 流束、重み、速度分布、`infinity_barrier`は
[シミュレーション境界からのreservoir流入](ReservoirInjection.html)で確認します。

### `plane_source`

**前提:** `[domain]`、正の`sim.batch_duration`、box内部のaxis-aligned矩形面が必要です。

```toml
[[particles.species]]
source_mode = "plane_source"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
target_macro_particles_per_batch = 300
pos_low = [0.2, 0.2, 2.0]
pos_high = [0.8, 0.8, 2.0]
source_normal = [0.0, 0.0, -1.0]
```

**期待する出力:** z=2.0の矩形面から-z方向へ一方向fluxを生成します。

**解釈:** 外部境界ではないため、`reservoir.inflow_model="infinity_barrier"`と`phi_infty`は適用しません。
旧`source_mode = "reservoir_face"`は互換入力として残りますが、新しいcaseでは用途に応じて
`boundary_inflow`または`plane_source`を選びます。

### `photo_raycast`

**前提:** 正の`sim.batch_duration`が必要です。

```toml
[[particles.species]]
source_mode = "photo_raycast"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
temperature_ev = 2.2
emit_current_density_a_m2 = 4.5e-6
rays_per_batch = 1000
deposit_opposite_charge_on_emit = true
inject_face = "z_high"
inject_region_mode = "face_fraction"
uv_low = [0.0, 0.0]
uv_high = [1.0, 1.0]
ray_direction = [0.0, 0.0, -1.0]
```

**期待する出力:** 最初にメッシュへ命中したrayから光電子を生成します。

**解釈:** `rays_per_batch`は照射ray数であり、生成粒子数は命中率で変わります。
`deposit_opposite_charge_on_emit=true`は放出元へ逆符号の電荷を残します。

**次の選択:** 放出、再吸収、closed PEの`neutral_return`は
[光電子の放出とライフサイクル](PhotoelectronEmission.html)にまとめています。

closed PEとして注入面をspecies単位で反射する場合:

```toml
[particles.species.boundary]
z_high = "reflect"
```

このtableは直前の`[[particles.species]]`に属します。`surface_charge_closure="neutral_return"`を使う場合は、
effectiveな`inject_face`境界が`reflect`または`redistributed_reflect`である必要があります。
例は接線位置を維持する`reflect`です。return位置を面内で一様再配置する選択肢は
詳しくは[光電子の放出とライフサイクル](PhotoelectronEmission.html)にあります。

## 2軸周期境界を使う

**前提:** `[domain]`でbox geometryとx/y周期を指定し、`field_solver="fmm"`を使います。
周期性は`domain.periodic_axes`だけで指定します。

**操作:**

```toml
[sim]
field_solver = "fmm"

[domain]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 4.0]
periodic_axes = ["x", "y"]

[field_boundary]
mode = "periodic2"

[particle_boundary]
z_low = "open"
z_high = "open"
ordinary_open_model = "escape"
```

**期待する出力:** x/yはfield、collision、raycastで同じ周期topologyを使い、z面は指定した粒子作用を使います。

**解釈:** `[particle_boundary]`に`periodic`は指定できず、周期面の作用は`[domain]`から決まります。

**次の選択:** operatorの選択は`[periodic2]`で行います。境界reservoir + closed PEの基準構成は
[periodic2有限画像構成](FinitePeriodicConfiguration.html)と
[`examples/periodic2_closed_photoelectron.toml`](../examples/periodic2_closed_photoelectron.toml)を使います。

## 履歴を保存する

**操作:**

```toml
[output]
dir = "outputs/latest"
history_stride = 1
write_mesh_potential = true
write_potential_history = true
```

**期待する出力:** `charge_history.csv`と`potential_history.csv`が作成されます。

**解釈:** 電位履歴は各保存点で再評価するため、要素数が多い場合は高コストです。

**次の選択:** まず電荷履歴だけを使い、相対電位が必要な場合に`write_potential_history`を有効にします。

## checkpointから再開する

**前提:** 読み込み元に有効なcheckpoint一式が必要です。

**操作:**

```toml
[output]
dir = "outputs/continuation"
resume = true
restart_from = "../previous/outputs/latest"
```

**期待する出力:** checkpointを読み込み、新しい出力を`outputs/continuation`へ保存します。

**解釈:** `sim.batch_count`は累積到達batch数です。checkpointが`batches=100`で
`batch_count=150`なら、追加実行は50 batchです。

**次の選択:** 同じdirectoryへ続けて書く場合は`restart_from`を省略し、`dir`をcheckpointのdirectoryにします。

title: シミュレーションケースを設計する

Lang: [日本語](ConfigurationRecipes.md) | [English](ConfigurationRecipes.en.md)

# シミュレーションケースを設計する

典型的なケースは、公式入門ケースの `beach.toml` を基に、メッシュと粒子源を目的に合わせて置き換えて作ります。
このページでは物理的な構成の選び方を扱います。全キーの定義は[入力パラメータリファレンス](Parameters.html)、
ファイルの生成や検査は[`beach.toml`を作成・検証する](Configuration.html)にまとめています。

以下は[公式入門ケース](Tutorial.html)を基準にした差分例です。
断片だけでは実行できない場合があります。

## 公式の実行手順

```bash
beachx config init
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

`beach.toml`はそのままFortran実行系に渡せます。box基準の座標・配置パラメータも読み込み時に実座標へ変換されます。

## レシピ一覧

| レシピ | 変更する主な場所 | 使う場面 |
| --- | --- | --- |
| 組み込みメッシュを追加 | `[mesh]`, `[[mesh.templates]]` | 平面、球、箱、円柱などを組み合わせる |
| OBJ メッシュ | `[mesh]` | 外部形状を使う |
| 粒子注入を選ぶ | `[[particles.species]]` | 流入、初期粒子、光電子放出を切り替える |
| 2軸周期境界（有限画像和） | `[sim]` | 指定した範囲の周期画像を使う |
| 高度な外部シース連成 | `[periodic2]`, `[external_boundary]` | 無限周期、`kinetic_1d`、UV 光電子を結合する |
| 履歴出力 | `[output]` | 時間発展を可視化する |
| 再開実行 | `[output]` | checkpoint から続ける |

## 組み込みメッシュを追加する

`[mesh].mode="template"` では、`[[mesh.templates]]` を 1 件書くたびに形状を追加できます。
有効な template には別々の `mesh_id` が割り当てられます。公式入門ケースの平面を残したまま球体を追加する場合は、
既存の `[[mesh.templates]]` の後ろへ新しい entry を追記します。

利用できる形状と、寸法・解像度を決める主なキーは次のとおりです。

| `kind` | 形状 | 寸法 | 解像度 |
| --- | --- | --- | --- |
| `plane` | XY 長方形平面 | `size_x`, `size_y` | `nx`, `ny` |
| `plate_hole`, `plane_hole` | 円形穴付き長方形平面 | `size_x`, `size_y`, `radius` | `n_theta`, `n_r` |
| `disk` | 円板 | `radius` | `n_theta`, `n_r` |
| `annulus` | 同心リング | `radius`, `inner_radius` | `n_theta`, `n_r` |
| `box` | 閉じた直方体表面 | `size = [sx, sy, sz]` | `nx`, `ny`, `nz` |
| `cylinder` | z 軸方向の円柱 | `radius`, `height`, `cap` | `n_theta`, `n_z` |
| `sphere` | 球面 | `radius` | `n_lon`, `n_lat` |

平面の例:

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
```

この平面へ球体を追加する例:

```toml
[[mesh.templates]]
kind = "sphere"
enabled = true
surface_model = "insulator"
center = [0.5, 0.5, 0.55]
radius = 0.12
n_lon = 24
n_lat = 12
```

`center` は形状の中心です。分割数を増やすと要素数と場計算・衝突判定のコストが増えるため、最初は小さい分割数で
配置と衝突を確認してから解像度を上げます。各形状の要素数と制約は
[組み込み形状リファレンス](Parameters.html#meshtemplates-組み込み形状)にまとめています。

通常の帯電計算には `surface_model="insulator"` を使います。`conductor` は `field_bc_mode="free"` だけに対応し、
`periodic2` とは併用できません。`dielectric` の `epsilon_r` は現行実装ではメタデータであり、誘電体分極を解きません。

### OBJ ファイルを使う

組み込み形状では表せない形状は、`[mesh]` を OBJ mode に置き換えます。

```toml
[mesh]
mode = "obj"
obj_path = "mesh/object.obj"
surface_model = "insulator"
obj_scale = 1.0
obj_rotation = [0.0, 0.0, 0.0]
obj_offset = [0.0, 0.0, 0.0]
```

変換は scale → rotation → offset の順です。1 つの OBJ ファイル全体が 1 つの `mesh_id` になります。
複数の独立した object を分けて扱う場合は、組み込み template を複数定義する方法を優先してください。

## 粒子注入を選ぶ

`[[particles.species]]` を 1 件書くたびに粒子種を追加します。公式入門ケースの電子源を変更する場合は既存 entry を編集し、
電子とイオンを同時に流入させる場合などは entry を追加します。
`source_mode` を切り替えるときは、既存 entry 全体を対応する例で置き換えてください。別 mode の
`npcls_per_step`、`w_particle`、`number_density_*` などを残すと validation error になります。

| `source_mode` | 使う場面 | 主なキー |
| --- | --- | --- |
| `volume_seed` | 箱内に初期粒子を置く小テスト | `npcls_per_step`, `pos_low`, `pos_high` |
| `reservoir_face` | 面から Maxwellian などの流入を与える通常ケース | `number_density_cm3`, `temperature_ev`, `inject_face`, `target_macro_particles_per_batch` |
| `photo_raycast` | 表面からの光電子放出を raycast で扱う | `rays_per_batch`, `emit_current_density_a_m2`, `ray_direction` |

`volume_seed` の例:

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

`volume_seed` は各 batch に `npcls_per_step` 個を生成します。物理的な面流束から個数を計算する方式ではありません。

`reservoir_face` または `photo_raycast` を使う前に、既存の `[sim]` へ対象時間 scale に合わせた正の
`batch_duration` を追加します。次は `batch_duration = 1.0e-5` とした場合の `reservoir_face` の最小例です。

```toml
[[particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
target_macro_particles_per_batch = 300
inject_face = "z_high"
inject_region_mode = "face_fraction"
uv_low = [0.0, 0.0]
uv_high = [1.0, 1.0]
drift_velocity = [0.0, 0.0, -4.0e5]
```

`reservoir_face` では `sim.use_box=true`、正の `sim.batch_duration`、`inject_face` が必要です。
`inject_region_mode="face_fraction"` の `uv_low` / `uv_high` は注入面内の割合で開口を指定します。
`target_macro_particles_per_batch` は 1 batch あたりの計算粒子数を固定し、粒子重みを流入量から解きます。
物理粒子数を重みで直接指定したい場合は `w_particle` を使います。
`temperature_k` と `temperature_ev` は同時に指定できません。

`photo_raycast` の例:

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

`photo_raycast` も正の `sim.batch_duration` を必要とします。`rays_per_batch` は照射 ray 数であり、実際に生成される
粒子数は最初のメッシュ命中率で決まります。`deposit_opposite_charge_on_emit=true` は放出元へ逆符号の電荷を残します。
放出、再吸収、open 面での処理は[光電子の放出とライフサイクル](PhotoelectronEmission.html)で説明します。

3 種類の粒子源の選択と生成後の共通処理は[粒子源の全体像](ParticleSourcesBoundaries.html)、
`reservoir_face` の流束・重み・速度分布は[`reservoir_face` の流入量と速度サンプリング](ReservoirInjection.html)にまとめています。

## 2軸周期境界を有限画像和で使う

`periodic2`は、3軸のうち2軸を周期軸とする場の境界条件です。このレシピでは、primary cellと
`field_periodic_image_layers`で指定した範囲の周期画像だけを足す有限画像和を使います。

```toml
[sim]
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 4.0]

bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"

field_solver = "fmm"
field_bc_mode = "periodic2"
field_periodic_image_layers = 1
field_periodic_far_correction = "none"
```

要件:

- `sim.use_box = true`
- ちょうど 2 軸が periodic
- 各周期軸の `box_max - box_min > 0`
- `sim.field_solver = "fmm"`
- `field_periodic_image_layers >= 0`

`field_periodic_image_layers = 1`なら、primary cellを含む$3\times3$ cellsをfield sourceとして足します。
`field_periodic_far_correction = "none"`なので、その外側の周期画像をEwald和やcached operatorで補いません。
画像層を増やし、電場、粒子flux、帯電分布などの目的量が変わらなくなることを確認してください。

これは無限個の周期画像を足した無限周期解ではありません。画像層の意味、z-highの局所太陽風reservoir、
光電子だけを閉じる`neutral_return`、top面基準電位、scalar potential barrierとの違いは
[periodic2有限画像構成](FinitePeriodicConfiguration.html)にまとめています。統合された実行例は
[`examples/periodic2_closed_photoelectron.toml`](../examples/periodic2_closed_photoelectron.toml)です。
無限周期operatorと外部sheathを使う場合は、次のレシピを選びます。

## 高度な外部シース連成は結合計算構成を使う

無限周期の `periodic2`、外部 `kinetic_1d` シース、reservoir の流入・帰還は、複数 section が同じ電位と
particle transfer を共有する高度な構成です。このページの断片を組み合わせず、
[periodic2 無限周期＋outer plasma の結合計算構成](InfinitePeriodicOuterConfiguration.html)を正本として設定してください。
標準の小規模 contract fixture は
[`examples/periodic2_kinetic_outer.toml`](../examples/periodic2_kinetic_outer.toml)です。

UV 光電子を外部シースの平均密度へ含める場合も、同じ結合計算構成の
「平均 outer 密度と tracked 光電子を別々に選ぶ」に従います。局所的な `photo_raycast` の設定、放出元電荷、
再吸収の確認は[光電子の放出とライフサイクル](PhotoelectronEmission.html)にまとめています。
`kinetic_mean`、tracked return、surface deposit は役割が異なるため、個別の TOML 断片として追加しません。

## 履歴を出す

```toml
[output]
dir = "outputs/latest"
history_stride = 1
write_mesh_potential = true
write_potential_history = true
```

`write_potential_history` は履歴ごとに電位を再評価するため、要素数が多いケースでは重くなります。
`sim.use_box=true`なら、同じbatchのz-high面平均を`top_reference_history.csv`にも記録します。
まず `charge_history.csv` だけで収束傾向を確認し、相対電位が必要なときに電位履歴を有効にしてください。

## 再開する

同じ出力ディレクトリから続ける場合:

```toml
[output]
dir = "outputs/latest"
resume = true
```

checkpoint 読み込み元と新しい出力先を分ける場合:

```toml
[output]
dir = "outputs/continuation"
resume = true
restart_from = "../previous/outputs/latest"
```

`sim.batch_count` は累積到達 batch 数です。checkpoint が `batches=100` で、新しい設定が `batch_count=150` なら追加で 50 batch だけ進みます。

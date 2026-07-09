title: 設定レシピ

Lang: [日本語](ConfigurationRecipes.md) | [English](ConfigurationRecipes.en.md)

# 設定レシピ

このページは、`beach.toml` をどう変えれば典型ケースを作れるかをまとめます。
全キーの定義は [入力パラメータリファレンス](Parameters.html)、高水準記法の詳細は [beachx config / 高水準記法ガイド](Configuration.html) を参照してください。

## 公式の実行手順

```bash
beachx config init
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

`beach.toml` はそのまま Fortran 実行系に渡せます。高水準記法も Fortran parser が読み込み時に正規化します。

## レシピ一覧

| レシピ | 変更する主な場所 | 使う場面 |
| --- | --- | --- |
| 平面メッシュで最小実行 | `[mesh]`, `[[mesh.templates]]` | まず動作確認する |
| 粒子注入を選ぶ | `[[particles.species]]` | 流入、初期粒子、光電子放出を切り替える |
| 2軸周期境界 | `[sim]` | 周期セルを使う |
| OBJ メッシュ | `[mesh]` | 外部形状を使う |
| 履歴出力 | `[output]` | 時間発展を可視化する |
| 再開実行 | `[output]` | checkpoint から続ける |

## 平面メッシュで最小実行

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

`nx` / `ny` を増やすと要素数が増えます。初回は小さく動かし、出力と粒子数を確認してから解像度を上げてください。

## 粒子注入を選ぶ

| `source_mode` | 使う場面 | 主なキー |
| --- | --- | --- |
| `volume_seed` | 箱内に初期粒子を置く小テスト | `npcls_per_step`, `pos_low`, `pos_high` |
| `reservoir_face` | 面から Maxwellian などの流入を与える通常ケース | `number_density_cm3`, `temperature_ev`, `inject_face`, `target_macro_particles_per_batch` |
| `photo_raycast` | 表面からの光電子放出を raycast で扱う | `rays_per_batch`, `emit_current_density_a_m2`, `ray_direction` |

`reservoir_face` の最小例:

```toml
[[particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
target_macro_particles_per_batch = 300
inject_face = "z_high"
pos_low = [0.0, 0.0, 4.0]
pos_high = [1.0, 1.0, 4.0]
drift_velocity = [0.0, 0.0, -4.0e5]
```

`target_macro_particles_per_batch` は 1 batch あたりの計算粒子数を固定し、粒子重みを流入量から解きます。
物理粒子数を重みで直接指定したい場合は `w_particle` を使います。

## 2軸周期境界

`periodic2` は 3 軸のうち 2 軸を周期軸にする場境界です。現行実装では `field_solver = "fmm"` が必要です。

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

## OBJ メッシュを使う

```toml
[mesh]
mode = "obj"
obj_path = "mesh/object.obj"
surface_model = "insulator"
```

OBJ では全要素に同じ `surface_model` が割り当てられます。
template mesh のように object ごとに surface model を分けたい場合は、複数 template を使うか、現行の制約を確認してください。

## 履歴を出す

```toml
[output]
dir = "outputs/latest"
history_stride = 1
write_mesh_potential = true
write_potential_history = true
```

`write_potential_history` は履歴ごとに電位を再評価するため、要素数が多いケースでは重くなります。
まず `charge_history.csv` だけで収束傾向を確認し、必要なときに電位履歴を有効にしてください。

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

title: 10分チュートリアル

Lang: [日本語](Tutorial.md) | [English](Tutorial.en.md)

# 10分チュートリアル

この公式入門ケースは、1個の電子を絶縁体平面へ向けて追跡し、吸収と表面電荷堆積までを確認します。
FMM、周期境界、光電子、導体・誘電体モデルは使いません。

## 1. 設定を作る

```bash
mkdir beach-tutorial
cd beach-tutorial
beachx config init beach.toml
beachx lint beach.toml
```

生成内容はリポジトリの
[`examples/tutorial_insulator.toml`](https://github.com/Nkzono99/BEACH/blob/main/examples/tutorial_insulator.toml)
と同一です。完全な設定は次のとおりです。

```toml
[sim]
dt = 1.0e-7
batch_count = 1
max_step = 10
softening = 1.0e-6
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 1.0]
bc_x_low = "reflect"
bc_x_high = "reflect"
bc_y_low = "reflect"
bc_y_high = "reflect"
bc_z_low = "open"
bc_z_high = "open"
rng_seed = 12345
open_boundary_model = "escape"
field_solver = "direct"
field_bc_mode = "free"
field_periodic_far_correction = "none"

[particles]
[[particles.species]]
source_mode = "volume_seed"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
w_particle = 1.0
npcls_per_step = 1
pos_low = [0.5, 0.5, 0.8]
pos_high = [0.5, 0.5, 0.8]
drift_velocity = [0.0, 0.0, -1.0e6]
temperature_k = 0.0

[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
enabled = true
surface_model = "insulator"
size_x = 1.0
size_y = 1.0
nx = 4
ny = 4
center = [0.5, 0.5, 0.2]

[output]
write_files = true
dir = "outputs/latest"
history_stride = 1
```

## 2. 実行する

```bash
beach beach.toml
beachx inspect outputs/latest
```

正常終了すると`outputs/latest/summary.txt`と`charges.csv`が生成されます。この決定論的ケースでは
`batches=1`、`processed_particles=1`となり、粒子は平面へ吸収されることを期待します。

```text
processed_particles=1
absorbed=1
batches=1
survived_max_step=0
```

## 3. 結果を見る

```bash
beachx inspect outputs/latest --save-mesh outputs/latest/charge.png
```

`charges.csv`では衝突要素だけに負電荷が堆積します。これは実行経路の確認用で、定常帯電状態ではありません。

![公式入門ケースの表面電荷密度](media/tutorial_insulator_charge.png)

## 4. 最初に変える値

| 目的 | キー |
| --- | --- |
| 粒子を増やす | `npcls_per_step` |
| 発射位置を広げる | `pos_low`, `pos_high` |
| 速度を変える | `drift_velocity`, `temperature_k` |
| 表面解像度を上げる | `nx`, `ny` |
| 長く計算する | `batch_count`, `max_step` |
| 履歴を間引く | `history_stride` |

次に[設定レシピ](ConfigurationRecipes.html)で流入・周期境界・OBJへ進めます。研究結果へ使う前に
[計算結果の妥当性確認](ValidationGuide.html)を実施してください。

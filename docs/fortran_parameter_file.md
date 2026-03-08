# Fortran パラメータファイル仕様（`beach.toml`）

この文書は、**現行実装（Fortran 実行系）で推奨する設定方法**をまとめたものです。  
「以前の書き方との互換」は最後に短く載せ、本文は初見向けに現在の推奨仕様だけを先に説明します。

## 1. 読み込みルール

- `beach .../config.toml` で明示指定できます。
- 引数なし実行時は、カレントディレクトリの `beach.toml` を自動で読み込みます。
- 開発時は `fpm run -- .../config.toml` でも同じ設定読込ルールです。
- 本実装は TOML のうち `key = value` とセクション/配列テーブルを使う軽量パーサです。
- 未知のセクション名やキー名はエラーになります。キーの打ち間違いを検出しやすい挙動です。

## 2. 初見向けの最小推奨例

まずは `reservoir_face`（物理流入ベース）を推奨します。

```toml
[sim]
dt = 2.0e-8
batch_duration_step = 60000.0
batch_count = 100
max_step = 10000
softening = 1.0e-6
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 10.0]
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
rng_seed = 12345

[[particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
target_macro_particles_per_batch = 5000
inject_face = "z_high"
pos_low = [0.0, 0.0, 10.0]
pos_high = [1.0, 1.0, 10.0]
drift_velocity = [0.0, 0.0, -4.0e5]

[[particles.species]]
source_mode = "reservoir_face"
number_density_cm3 = 5.0
temperature_ev = 10.0
q_particle = 1.602176634e-19
m_particle = 1.672482821616e-27
target_macro_particles_per_batch = -1
inject_face = "z_high"
pos_low = [0.0, 0.0, 10.0]
pos_high = [1.0, 1.0, 10.0]
drift_velocity = [0.0, 0.0, -4.0e5]

[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
enabled = true
size_x = 1.0
size_y = 1.0
nx = 20
ny = 20
center = [0.5, 0.5, 0.02]

[output]
write_files = true
dir = "outputs/latest"
history_stride = 1
```

## 3. セクションとキー

### 3.1 `[sim]`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `dt` | float | `1.0e-9` | 時間刻み [s] |
| `rng_seed` | int | `12345` | 乱数シード |
| `batch_count` | int | `1` | 1回の実行で処理するバッチ数 |
| `batch_duration` | float | `0.0` | 1バッチの物理時間 [s] |
| `batch_duration_step` | float | `0.0` | `batch_duration = dt * batch_duration_step` |
| `max_step` | int | `400` | 粒子あたり最大ステップ数 |
| `tol_rel` | float | `1.0e-8` | 相対変化量の監視値（停止条件には未使用） |
| `q_floor` | float | `1.0e-30` | `rel_change` 計算時の分母下限 |
| `softening` | float | `1.0e-6` | 電場計算の softening 長さ [m] |
| `use_hybrid` | bool | `true` | 予約キー（現行電場計算では未使用） |
| `r_switch_factor` | float | `3.0` | 予約キー（未使用） |
| `n_sub` | int | `2` | 予約キー（未使用） |
| `softening_factor` | float | `0.1` | 予約キー（未使用） |
| `b0` | float[3] | `[0,0,0]` | 一様磁場 [T] |
| `reservoir_potential_model` | string | `"none"` | `none` / `infinity_barrier` |
| `phi_infty` | float | `0.0` | 無限遠基準電位 [V] |
| `injection_face_phi_grid_n` | int | `3` | 注入面平均電位の格子分割数 `N x N` |
| `raycast_max_bounce` | int | `16` | `photo_raycast` レイ追跡の最大反射回数 |
| `use_box` | bool | `false` | ボックス境界を有効化 |
| `box_min` | float[3] | `[-1,-1,-1]` | ボックス下限 [m] |
| `box_max` | float[3] | `[1,1,1]` | ボックス上限 [m] |
| `bc_x_low` ほか | string | `"open"` | `open` / `reflect` / `periodic`（`open` は `outflow`,`escape` も可） |

`batch_duration` の解決ルール:

- `batch_duration` と `batch_duration_step` の同時指定はエラーです。
- `batch_duration_step` 指定時は `batch_duration = dt * batch_duration_step` に解決します。
- `reservoir_face` / `photo_raycast` を使う場合、解決後の `batch_duration > 0` が必須です。

重要な実行挙動:

- 実行ループは `batch_count` 分だけ進みます。
- `tol_rel` は出力監視値であり、現行実装では早期終了条件には使いません。

### 3.2 `[[particles.species]]`

`[[particles.species]]` は 1 件以上必須です。

### 共通キー

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `enabled` | bool | `true` | 種を有効化 |
| `source_mode` | string | `"volume_seed"` | `volume_seed` / `reservoir_face` / `photo_raycast` |
| `q_particle` | float | `-1.602176634e-19` | 粒子電荷 [C] |
| `m_particle` | float | `9.10938356e-31` | 粒子質量 [kg] |
| `pos_low` | float[3] | `[-0.4,-0.4,0.2]` | 位置下限 [m] |
| `pos_high` | float[3] | `[0.4,0.4,0.5]` | 位置上限 [m] |
| `drift_velocity` | float[3] | `[0,0,-8e5]` | ドリフト速度 [m/s] |
| `temperature_k` | float | `2.0e4` | 温度 [K] |
| `temperature_ev` | float | 未指定 | 温度 [eV]（`temperature_k` と排他） |
| `inject_face` | string | 未指定 | `x_low/x_high/y_low/y_high/z_low/z_high` |

### `source_mode = "volume_seed"`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `npcls_per_step` | int | `0` | 1バッチに生成するマクロ粒子数 |
| `w_particle` | float | `1.0` | マクロ粒子重み |

制約:

- 有効 species 全体で `npcls_per_step` 合計が 1 以上必要です（`reservoir_face` / `photo_raycast` を使わない場合）。
- `target_macro_particles_per_batch` は使用できません。

### `source_mode = "reservoir_face"`

| キー | 型 | 説明 |
|---|---|---|
| `number_density_cm3` / `number_density_m3` | float | 上流密度（どちらか片方必須） |
| `w_particle` | float | マクロ粒子重み（正値） |
| `target_macro_particles_per_batch` | int | `w_particle` 自動解決用（`>0` または `-1`） |

制約:

- `sim.use_box = true` 必須
- `sim.batch_duration > 0` 必須
- `inject_face` 必須
- `pos_low` / `pos_high` は指定 face 上になければエラー
- `w_particle` と `target_macro_particles_per_batch` は同時指定不可
- `target_macro_particles_per_batch = -1` は species[2] 以降のみ可（species[1] の `w_particle` を共有）

粒子数決定（1種あたり）:

- 期待マクロ粒子数 `n_macro_expected = gamma_in * A * batch_duration / w_particle`
- 実際の注入数は残差繰越つきで `floor(residual + n_macro_expected)`
- `target_macro_particles_per_batch > 0` のときは、その値になるよう `w_particle` を自動計算

### `source_mode = "photo_raycast"`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `emit_current_density_a_m2` | float | `0.0` | レイ垂直面基準の放出電流密度 [A/m^2] |
| `rays_per_batch` | int | `0` | 1バッチの発射レイ数 |
| `deposit_opposite_charge_on_emit` | bool | `false` | 放出元要素に逆符号電荷を堆積 |
| `normal_drift_speed` | float | `0.0` | 放出法線方向ドリフト [m/s] |
| `ray_direction` | float[3] | 未指定時は注入面内向き法線 | レイ方向 |

制約:

- `sim.use_box = true` 必須
- `sim.batch_duration > 0` 必須
- `emit_current_density_a_m2 > 0`, `rays_per_batch > 0` 必須
- `inject_face` 必須
- `q_particle` は非ゼロ、`m_particle > 0` 必須
- `ray_direction` 指定時は正規化可能で、注入面内向き法線への内積が正であること
- `npcls_per_step` / `number_density_*` / `w_particle` / `target_macro_particles_per_batch` は使用不可

`photo_raycast` の重み:

- `w_hit = J_perp * A_perp * batch_duration / (|q_particle| * rays_per_batch)`
- 実際の放出数はレイの命中率で決まるため、バッチごとの生成粒子数は `rays_per_batch` 以下になります。

### 3.3 `[mesh]`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `mode` | string | `"auto"` | `auto` / `obj` / `template` |
| `obj_path` | string | `"examples/simple_plate.obj"` | OBJ パス |

`mode = "auto"` のときは `obj_path` が存在すれば OBJ、なければ template を使います。

### 3.4 `[[mesh.templates]]`

共通キー:

- `enabled` (bool)
- `kind` (`plane` / `box` / `cylinder` / `sphere`)
- `center` (float[3])

`kind` ごとの主要キー:

- `plane`: `size_x`, `size_y`, `nx`, `ny`
- `box`: `size`, `nx`, `ny`, `nz`
- `cylinder`: `radius`, `height`, `n_theta`, `n_z`, `cap`
- `sphere`: `radius`, `n_lon`, `n_lat`

注意:

- `[[mesh.templates]]` を書いた場合、実際に使うテンプレート数は「定義した件数」で解決されます。

### 3.5 `[output]`

| キー | 型 | 既定値 | 説明 |
|---|---|---:|---|
| `write_files` | bool | `true` | ファイル出力の有効/無効 |
| `dir` | string | `"outputs/latest"` | 出力先ディレクトリ |
| `history_stride` | int | `1` | `charge_history.csv` の出力間隔（バッチ単位） |
| `resume` | bool | `false` | 既存チェックポイントから再開 |

出力ファイル:

- `summary.txt`
- `charges.csv`
- `mesh_triangles.csv`
- `mesh_sources.csv`
- `charge_history.csv`（`history_stride > 0` のとき）
- `rng_state.txt`
- `macro_residuals.csv`

`mesh_triangles.csv` には `mesh_id` 列が追加され、`mesh_sources.csv` で `mesh_id` ごとの元メッシュ種別と要素数を確認できます。

MPI実行（`world_size > 1`）では乱数状態・残差はrank別ファイルになります。

- `rng_state_rank00000.txt`, `rng_state_rank00001.txt`, ...
- `macro_residuals_rank00000.csv`, `macro_residuals_rank00001.csv`, ...

また、`summary.txt` には `mpi_world_size` が記録されます。

`resume = true` の要件:

- `write_files = true` 必須
- `output.dir` に `summary.txt` / `charges.csv` / `rng_state.txt` が必要
- `macro_residuals.csv` は存在すれば読み込みます

MPI実行での追加要件:

- `summary.txt` の `mpi_world_size` が現在のrank数と一致している必要があります
- 各rankに対応する `rng_state_rankNNNNN.txt` が必要です
- `macro_residuals_rankNNNNN.csv` は存在すれば読み込みます

## 4. キー検証

- 未知のセクション名・キー名はすべてエラーになります。
- `v0.3.0` 以降は旧キー互換を持たず、旧名は「未知キー」として扱います。
- `[particles]` は `[[particles.species]]` のコンテナとしてのみ使い、`[particles]` 直下に `key = value` は書かないでください。

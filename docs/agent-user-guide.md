# BEACH Agent User Guide

> AI Agent が BEACH シミュレーションを操作するためのリファレンスガイド。
> CLAUDE.md から `@import docs/agent-user-guide.md` で読み込むことを想定。

---

## 概要

BEACH (BEM + Accumulated CHarge) は、絶縁体表面への帯電蓄積をシミュレーションする境界要素法+粒子追跡ハイブリッドシミュレータである。

- **Fortran コア**: 粒子力学・電場ソルバー・衝突判定・電荷堆積
- **Python レイヤー**: 設定管理・後処理・可視化
- **バージョン**: 1.2.0

---

## クイックスタート

### 最小実行手順

```bash
# 1. インストール
pip install beach-bem

# 2. 設定ファイルを作成
beachx config init

# 3. 高水準記法を最終キーへレンダリング
beachx config render

# 4. シミュレーション実行
beach beach.toml

# 5. 結果確認
beachx inspect outputs/latest
```

## ビルドと実行

### 必要環境

| 項目 | 要件 |
|------|------|
| Fortran コンパイラ | gfortran (または互換) |
| fpm | v0.10+ |
| Python | 3.10+ |
| 主要 Python 依存 | matplotlib >= 3.8, numpy >= 1.24 |

### ビルドコマンド

```bash
make check                     # 開発用ビルド確認（version は dev 固定）
make build                     # git describe 付き version でビルド
make run CONFIG=beach.toml     # dev 固定 version で実行
make install-generic           # gfortran ポータブル
make install-camphor           # Intel コンパイラ最適化
make test                      # L1: Python + quick Fortran tests
```

`make check` / `make test` / `make run` は `BEACH_VERSION_MODE=dev` を使い、Fortran に渡す
version macro を安定させる。git hash が変わっても fpm の compile-flag hash が変わらないため、
開発中の差分コンパイルを再利用しやすい。git hash 付きの実行ファイルが必要な場合は
`make build VERSION_MODE=git` または `make install` を使う。

### 低レベル fpm 直接実行

```bash
fpm run --profile release --flag "-fopenmp" -- beach.toml
```

### MPI 並列実行

```bash
FPM_FC=mpiifort fpm run --profile release \
  --flag "-fpp -DUSE_MPI -qopenmp" \
  --runner "mpirun -n 4" -- beach.toml
```

### テスト

```bash
make test-l0      # L0: static/schema/build check
make test         # L1: normal development loop
make test-l2      # L2: contract/integration
make test-l3      # L3: heavy/release gate
make test-heavy   # heavy Fortran targets only
make test-full    # unfiltered fpm test
make test-mpi     # MPI テスト
pytest -q         # Python テストのみ
```

`make test` は L1 の alias で、通常の AI/開発内側ループではここまでを基本にする。
FMM 系の長時間 target は `make test-l3` / `make test-heavy` / `make test-fortran-heavy` / `make test-full` で明示実行する。
個別 target は `FPM_ACTION=test ./build.sh --target <name>` で確認できる。

---

## 設定パラメータ

### 設定ワークフロー

1. **beach.toml**: 通常の編集対象で、Fortran 実行ファイルが直接読む設定
2. **beachx lint**: TOML parse、JSON Schema、高水準記法、既知制約を検証
3. **beachx config render**: 高水準記法を `box_min` / `box_max` / `center` などの最終キーへ展開

### [sim] セクション — シミュレーション基本

| パラメータ | 型 | デフォルト | 説明 |
|------------|------|-----------|------|
| `dt` | float | 1.0e-9 | タイムステップ [s] |
| `rng_seed` | int | 12345 | 乱数シード |
| `batch_count` | int | 1 | 通常実行ではバッチ数。resume 時は累積の到達バッチ数 |
| `batch_duration` | float | — | バッチ持続時間 [s] (`batch_duration_step` と排他) |
| `batch_duration_step` | float | — | `batch_duration = dt * batch_duration_step` として解決 |
| `max_step` | int | 400 | 粒子あたり最大積分ステップ数 |
| `tol_rel` | float | 1.0e-8 | 監視メトリクス (早期停止条件ではない) |
| `q_floor` | float | 1.0e-30 | 相対変化計算の分母フロア |
| `softening` | float | 1.0e-6 | 電場ソフトニング長 [m] |

### 電場ソルバーパラメータ

| パラメータ | 型 | デフォルト | 選択肢 | 説明 |
|------------|------|-----------|--------|------|
| `field_solver` | string | "auto" | direct, treecode, fmm, auto | 電場評価手法 |
| `field_bc_mode` | string | "free" | free, periodic2 | 境界条件 (periodic2 は fmm 必須) |
| `field_periodic_image_layers` | int | 1 | >= 0 | periodic2 のイメージシェル層数 |
| `field_periodic_far_correction` | string | "auto" | auto, none, m2l_root_oracle | 遠方補正 |
| `field_periodic_ewald_alpha` | float | 0.0 | >= 0 | Ewald 分割パラメータ (0=自動) |
| `field_periodic_ewald_layers` | int | 4 | >= 0 | Ewald シェル深度 |
| `tree_theta` | float | 0.5 | (0, 1] | ツリー法 MAC パラメータ |
| `tree_leaf_max` | int | 16 | >= 1 | リーフノードあたり最大要素数 |
| `tree_min_nelem` | int | 256 | >= 1 | auto → treecode 切替閾値 |

**auto 推定テーブル** (tree_theta / tree_leaf_max 未指定時):

| 要素数 | theta | leaf_max |
|--------|-------|----------|
| < 1500 | 0.40 | 12 |
| 1500–9999 | 0.50 | 16 |
| 10000–49999 | 0.58 | 20 |
| 50000+ | 0.65 | 24 |

### 磁場・ポテンシャル

| パラメータ | 型 | デフォルト | 説明 |
|------------|------|-----------|------|
| `b0` | float[3] | [0, 0, 0] | 一様磁場 [T] |
| `reservoir_potential_model` | string | "none" | none, infinity_barrier |
| `phi_infty` | float | 0.0 | 無限遠参照電位 [V] |
| `injection_face_phi_grid_n` | int | 3 | 注入面ポテンシャルグリッド解像度 NxN |
| `raycast_max_bounce` | int | 16 | photo_raycast のレイ反射最大回数 |

### シースモデル

| パラメータ | 型 | デフォルト | 説明 |
|------------|------|-----------|------|
| `sheath_injection_model` | string | "none" | none, zhao_auto, zhao_a, zhao_b, zhao_c, floating_no_photo |
| `sheath_alpha_deg` | float | 60.0 | 太陽仰角 [deg] |
| `sheath_photoelectron_ref_density_cm3` | float | 64.0 | 参照光電子密度 [cm^-3] |
| `sheath_reference_coordinate` | float | — | シース参照面位置 [m] |
| `sheath_electron_drift_mode` | string | "normal" | normal, full |
| `sheath_ion_drift_mode` | string | "normal" | normal, full |

### ボックス境界条件

| パラメータ | 型 | デフォルト | 説明 |
|------------|------|-----------|------|
| `use_box` | bool | false | ボックス境界の有効化 |
| `box_min` | float[3] | [-1, -1, -1] | 下端 [m] |
| `box_max` | float[3] | [1, 1, 1] | 上端 [m] |
| `bc_{x,y,z}_{low,high}` | string | "open" | open, reflect, periodic |

**periodic2 制約**: ちょうど 2 軸が periodic、残り 1 軸は open/reflect

### [[particles.species]] セクション — 粒子種 (最低 1 種必須)

#### 共通パラメータ

| パラメータ | 型 | デフォルト | 説明 |
|------------|------|-----------|------|
| `enabled` | bool | true | 種の有効化 |
| `source_mode` | string | "volume_seed" | volume_seed, reservoir_face, photo_raycast |
| `q_particle` | float | -1.602e-19 | 粒子電荷 [C] |
| `m_particle` | float | 9.109e-31 | 粒子質量 [kg] |
| `pos_low` | float[3] | [-0.4, -0.4, 0.2] | 生成位置下限 [m] |
| `pos_high` | float[3] | [0.4, 0.4, 0.5] | 生成位置上限 [m] |
| `drift_velocity` | float[3] | [0, 0, -8e5] | ドリフト速度 [m/s] |
| `temperature_k` | float | 20000 | 温度 [K] (`temperature_ev` と排他) |
| `temperature_ev` | float | — | 温度 [eV] |

#### volume_seed モード

| パラメータ | 型 | デフォルト | 説明 |
|------------|------|-----------|------|
| `npcls_per_step` | int | 0 | バッチあたり生成粒子数 |
| `w_particle` | float | 1.0 | マクロ粒子重み |

**制約**: 全種の `npcls_per_step` 合計 >= 1

#### reservoir_face モード (物理フラックス注入)

| パラメータ | 型 | デフォルト | 説明 |
|------------|------|-----------|------|
| `number_density_cm3` | float | — | 上流密度 [cm^-3] (`number_density_m3` と排他) |
| `number_density_m3` | float | — | 上流密度 [m^-3] |
| `w_particle` | float | — | マクロ粒子重み (`target_macro_particles_per_batch` と排他) |
| `target_macro_particles_per_batch` | int | — | バッチあたり目標マクロ粒子数 (-1 で species[1] の重みを再利用) |
| `inject_face` | string | 必須 | x_low, x_high, y_low, y_high, z_low, z_high |

**制約**: `use_box = true` かつ `batch_duration > 0` が必要。`pos_low/pos_high` は指定面上に配置。

#### photo_raycast モード (光電子放出)

| パラメータ | 型 | デフォルト | 説明 |
|------------|------|-----------|------|
| `emit_current_density_a_m2` | float | 必須 | 放出電流密度 [A/m^2] |
| `rays_per_batch` | int | 必須 | バッチあたりレイ数 |
| `deposit_opposite_charge_on_emit` | bool | false | 放出元要素に逆符号電荷を堆積 |
| `normal_drift_speed` | float | 0.0 | 法線方向ドリフト速度 [m/s] |
| `ray_direction` | float[3] | 内向き法線 | レイ方向ベクトル |

### [mesh] セクション — ジオメトリ

| パラメータ | 型 | デフォルト | 説明 |
|------------|------|-----------|------|
| `mode` | string | "auto" | auto, obj, template |
| `obj_path` | string | examples/simple_plate.obj | OBJ ファイルパス |
| `surface_model` | string | "insulator" | OBJ メッシュ全体の表面モデル (`insulator`, `conductor`, `dielectric`) |
| `epsilon_r` | float | 1.0 | OBJ メッシュ全体の相対誘電率 (`>= 1`) |
| `obj_scale` | float | 1.0 | スケーリング係数 |
| `obj_rotation` | float[3] | [0, 0, 0] | 回転角 [deg] (外因性 x→y→z) |
| `obj_offset` | float[3] | [0, 0, 0] | 平行移動 [m] |

**変換順序**: scale → rotate → offset

#### [[mesh.templates]] — 手続き的メッシュ生成

共通: `enabled` (bool), `kind` (enum), `surface_model` (enum), `epsilon_r` (float), `center` (float[3])

`conductor` は mesh_id ごとの浮遊導体として等電位再配分されます。現行実装では `sim.field_bc_mode = "free"` のみ対応します。OBJ 入力はファイル全体を `mesh_id = 1` として読むため、1つの OBJ 内の離れた conductor 部品も同じ浮遊導体として扱われます。独立導体として扱う場合はテンプレート入力などで mesh_id を分けてください。`dielectric` は現行では object ごとの `epsilon_r` を保持するメタデータで、誘電体分極の物理分岐は今後の拡張点です。

| kind | 主要パラメータ |
|------|---------------|
| `plane` | `size_x`, `size_y`, `nx`, `ny` |
| `plate_hole` | `size_x`, `size_y`, `radius`, `n_theta`, `n_r` |
| `disk` | `radius`, `n_theta`, `n_r` |
| `annulus` | `radius`, `inner_radius`, `n_theta`, `n_r` |
| `box` | `size` (float[3]), `nx`, `ny`, `nz` |
| `cylinder` | `radius`, `height`, `n_theta`, `n_z`, `cap`, `cap_top`, `cap_bottom` |
| `sphere` | `radius`, `n_lon`, `n_lat` |

### [output] セクション — ファイル出力

| パラメータ | 型 | デフォルト | 説明 |
|------------|------|-----------|------|
| `write_files` | bool | true | ファイル出力の有効化 |
| `write_mesh_potential` | bool | false | 要素ポテンシャルを mesh_potential.csv に出力 |
| `write_potential_history` | bool | false | ポテンシャル履歴を出力 |
| `dir` | string | "outputs/latest" | 出力ディレクトリ |
| `history_stride` | int | 1 | 履歴出力間隔 [バッチ] (0 で無効化) |
| `resume` | bool | false | チェックポイントから再開 |

---

## 出力ファイル形式

出力先: `output.dir` で指定したディレクトリ (デフォルト `outputs/latest/`)

### 必須出力ファイル

| ファイル | 形式 | 内容 |
|----------|------|------|
| `summary.txt` | テキスト (Key-Value) | 実行メタデータ・統計情報 |
| `charges.csv` | CSV: `elem_idx, charge_C` | 最終要素電荷 |
| `mesh_triangles.csv` | CSV: `elem_idx, v0x, v0y, v0z, v1x, v1y, v1z, v2x, v2y, v2z, charge_C, mesh_id` | 三角形頂点・電荷・mesh_id |
| `mesh_sources.csv` | CSV: `mesh_id, source_kind, template_kind, surface_model, epsilon_r, elem_count` | メッシュソースメタデータ |
| `rng_state.txt` | テキスト | 乱数状態 (リジューム用) |

### オプション出力ファイル

| ファイル | 条件 | 形式 |
|----------|------|------|
| `charge_history.csv` | `history_stride > 0` | CSV: `batch, elem_idx, charge_C` |
| `potential_history.csv` | `write_potential_history = true` かつ `history_stride > 0` | CSV: `batch, elem_idx, potential_V` |
| `mesh_potential.csv` | `write_mesh_potential = true` | CSV: `elem_idx, potential_V` |
| `macro_residuals.csv` | reservoir_face 使用時 | CSV: 注入残差状態 |
| `performance_profile.csv` | `BEACH_PROFILE=1` 環境変数設定時 | CSV: 各領域の計測時間 |

### MPI 実行時の追加ファイル

- `rng_state_rank00000.txt`, `rng_state_rank00001.txt`, ...
- `macro_residuals_rank00000.csv`, `macro_residuals_rank00001.csv`, ...

---

## Python CLI (beachx)

### config — 設定管理

```bash
beachx lint [beach.toml]                           # schema と意味制約をまとめて検査
beachx config init [beach.toml]                    # beach.toml を新規作成
beachx config render [beach.toml]                  # 高水準記法を最終キーへレンダリング
beachx config validate [beach.toml]                # 高水準記法と意味制約の検証
beachx config diff left.toml right.toml            # 設定比較
```

### inspect — 結果確認

```bash
beachx inspect [output_dir]                        # サマリー表示
beachx inspect outputs/latest --show               # 表示
beachx inspect outputs/latest --save-bar charges.png
beachx inspect outputs/latest --save-mesh charges_mesh.png
beachx inspect outputs/latest --save-potential-mesh potential_mesh.png
beachx inspect outputs/latest --apply-periodic2-mesh
beachx inspect outputs/latest --periodic2-repeat 1
```

### animate — アニメーション作成

```bash
beachx animate [output_dir]                        # 履歴アニメーション
beachx animate outputs/latest --quantity charge     # 電荷 (デフォルト)
beachx animate outputs/latest --quantity potential   # ポテンシャル
beachx animate outputs/latest --save-gif charge.gif
beachx animate outputs/latest --total-frames 200
beachx animate outputs/latest --frame-stride 2 --fps 15
```

### slices — ポテンシャル断面図

```bash
beachx slices [output_dir]
beachx slices outputs/latest --grid-n 200
beachx slices outputs/latest --xy-z 0.5 --yz-x 0.5 --xz-y 0.5
beachx slices outputs/latest --vmin -20 --vmax 20
beachx slices outputs/latest --save slices.png
```

### coulomb — クーロン力行列

```bash
beachx coulomb [output_dir]
beachx coulomb outputs/latest --component z
beachx coulomb outputs/latest --target-kinds sphere
beachx coulomb outputs/latest --save forces.png
```

### mobility — クーロン移動度解析

```bash
beachx mobility [output_dir]
beachx mobility outputs/latest --density-kg-m3 2500
beachx mobility outputs/latest --mu-static 0.4
beachx mobility outputs/latest --save-csv mobility.csv
```

### workload — ワークロード推定

```bash
beachx workload beach.toml
beachx workload beach.toml --threads 8
```

---

## Python API

```python
from beach import Beach

run = Beach("outputs/latest")

# 電荷分布の可視化
fig, ax = run.plot_charges(step=-1)
fig, ax = run.plot_mesh()
fig, ax = run.plot_mesh_source_boxplot(quantity="charge", step=-1)

# ポテンシャル解析
fig, ax = run.plot_potential()
slices = run.compute_potential_slices(grid_n=200)
fig, axes = run.plot_potential_slices(slices)

# クーロン力解析
fig, ax = run.plot_coulomb_force_matrix(component="z")
result = run.analyze_coulomb_mobility(density_kg_m3=2500, mu_static=0.4)

# 電場・電界線
field = run.compute_electric_field_points(points)
lines = run.trace_field_lines(seed_points, direction="forward")
fig = run.plot_field_lines_3d(lines)

# アニメーション
run.animate_mesh(quantity="charge", save_path="charge.gif")
```

### 主要関数一覧

| 関数 | 用途 |
|------|------|
| `plot_mesh()` | 3D メッシュ可視化 |
| `plot_potential()` | 3D ポテンシャルメッシュ |
| `plot_charges()` | 要素電荷の棒グラフ |
| `plot_mesh_source_boxplot()` | メッシュソース別箱ひげ図 |
| `animate_mesh()` | 履歴からのメッシュアニメーション |
| `compute_potential_mesh()` | 要素ポテンシャル再構成 |
| `compute_potential_points()` | 任意点でのポテンシャル評価 |
| `compute_potential_slices()` | 断面データ計算 |
| `compute_electric_field_points()` | 電場ベクトル評価 |
| `plot_potential_slices()` | 2D 断面プロット |
| `plot_coulomb_force_matrix()` | 力行列ヒートマップ |
| `analyze_coulomb_mobility()` | 移動度解析 |
| `trace_field_lines()` | 電界線積分 |
| `plot_field_lines_3d()` | 電界線 + メッシュ描画 |

---

## `beachx config` の高水準記法

`beachx config render` は `beach.toml` 内の補助キーを、Fortran 実行系が読む最終キーへ展開する。

- `sim.box_origin` + `sim.box_size` -> `sim.box_min` / `sim.box_max`
- `inject_region_mode = "face_fraction"` + `uv_low` / `uv_high` -> `pos_low` / `pos_high`
- `mesh.templates` の `placement_mode = "box_anchor"` -> `center`
- `mesh.groups.*` の `scale_from` / `placement_mode` -> template ごとの実寸・実座標

`schema_version`、`use_presets`、`override`、`base_case` は現行の direct `beach.toml` では使わない。

---

## 設定例

### 最小構成 (volume_seed)

```toml
[sim]
dt = 1.0e-8
batch_count = 10
max_step = 500
use_box = true
box_min = [0, 0, 0]
box_max = [1, 1, 1]

[[particles.species]]
source_mode = "volume_seed"
npcls_per_step = 100

[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
size_x = 1.0
size_y = 1.0
nx = 10
ny = 10

[output]
dir = "outputs/test"
```

### 周期境界 + volume seed

```toml
[sim]
batch_count = 200
field_solver = "fmm"
field_bc_mode = "periodic2"
use_box = true
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 10.0]
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"

[[particles.species]]
source_mode = "volume_seed"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
npcls_per_step = 10

[[particles.species]]
source_mode = "volume_seed"
q_particle = 1.602176634e-19
m_particle = 1.672482821616e-27
npcls_per_step = 10
```

### 完全な後処理パイプライン

```bash
beach beach.toml
beachx inspect outputs/latest --save-mesh charges.png
beachx animate outputs/latest --quantity charge --save-gif charge.gif
beachx slices outputs/latest --save slices.png
beachx coulomb outputs/latest --save forces.png
beachx mobility outputs/latest --density-kg-m3 2500 --mu-static 0.4
```

---

## パラメータ制約チェックリスト

| 機能 | 必須条件 |
|------|---------|
| `reservoir_face` | `use_box=true`, `batch_duration>0`, `inject_face` 指定 |
| `photo_raycast` | `use_box=true`, `batch_duration>0`, `emit_current_density_a_m2>0`, `rays_per_batch>=1` |
| `periodic2` | `field_solver=fmm`, ちょうど 2 軸が periodic, `use_box=true` |
| シースモデル | `reservoir_potential_model = "none"` と互換 |
| リジューム | `write_files=true`, チェックポイントファイル存在, MPI サイズ一致 |
| 性能プロファイル | 環境変数 `BEACH_PROFILE=1` |
| MPI 実行 | `-DUSE_MPI` でコンパイル, MPI コンパイララッパー使用 |

---

## プロジェクト構造

```
BEACH/
├── app/main.f90              # Fortran エントリポイント
├── src/                      # Fortran ライブラリモジュール
│   ├── config/               #   設定パーサー
│   ├── core/                 #   シミュレータメインループ
│   ├── mesh/                 #   三角形メッシュ管理
│   ├── particles/            #   粒子力学 (Boris pusher)
│   ├── physics/              #   電場ソルバー (direct, treecode, FMM)
│   └── runtime/              #   出力・リスタート
├── beach/                    # Python パッケージ
│   ├── cli/                  #   CLI サブコマンド
│   ├── config/               #   設定レンダリング・検証
│   └── fortran_results/      #   後処理・可視化
├── schemas/                  # JSON Schema (バリデーション用)
│   └── beach.schema.json     #   beach.toml スキーマ
├── examples/                 # サンプル設定・スクリプト
├── tests/                    # テストスイート
├── docs/                     # ドキュメント
├── fpm.toml                  # Fortran パッケージマニフェスト
├── pyproject.toml            # Python パッケージメタデータ
├── Makefile                  # ビルド自動化
├── SPEC.md                   # Fortran 実装仕様
└── AGENTS.md                 # マルチエージェント設定
```

---

## 詳細ドキュメントへの参照

| ドキュメント | 内容 |
|-------------|------|
| `SPEC.md` | Fortran 実装仕様 (権威的) |
| `docs/fortran_parameter_file.md` | パラメータ詳細仕様 |
| `docs/fortran_workflow.md` | 実行ワークフロー・I/O |
| `docs/fortran_fmm_core.md` | FMM 数学・Ewald |
| `docs/config_workflow.md` | `beachx config` と高水準記法 |
| `docs/python_postprocess_api.md` | Python API リファレンス |
| `schemas/beach.schema.json` | IDE バリデーション用 JSON Schema |

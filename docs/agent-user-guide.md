# BEACH Agent User Guide

> AI Agent が BEACH シミュレーションを操作するためのリファレンスガイド。
> CLAUDE.md から `@import docs/agent-user-guide.md` で読み込むことを想定。

---

## 概要

BEACH (BEM + Accumulated CHarge) は、絶縁体表面への帯電蓄積をシミュレーションする境界要素法+粒子追跡ハイブリッドシミュレータである。

- **Fortran コア**: 粒子力学・電場ソルバー・衝突判定・電荷堆積
- **Python レイヤー**: 設定管理・後処理・可視化
- **バージョン**: 1.0.0

---

## クイックスタート

### 最小実行手順

```bash
# 1. インストール
pip install beach-bem

# 2. 設定ファイルを作成
beachx config init case.toml

# 3. beach.toml にレンダリング
beachx config render case.toml

# 4. シミュレーション実行
beach beach.toml

# 5. 結果確認
beachx inspect outputs/latest
```

### プリセットを使った設定作成

```bash
beachx config init case.toml \
  --preset sim/periodic2_fmm \
  --preset species/solarwind_electron \
  --preset species/solarwind_ion \
  --preset mesh/plane_basic

beachx config render case.toml
beach beach.toml
```

---

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
make                           # auto-detect コンパイラでビルド
make install-generic           # gfortran ポータブル
make install-camphor           # Intel コンパイラ最適化
```

### 開発者向け直接実行

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
fpm test          # Fortran テスト
make test         # OpenMP 付き
make test-mpi     # MPI テスト
pytest -q         # Python テスト
```

---

## 設定パラメータ

### 設定ワークフロー

1. **case.toml**: プリセット参照と上書きを含む高レベル設定
2. **プリセット**: 再利用可能な TOML フラグメント
3. **beach.toml**: マージ済み最終設定 (Fortran 実行ファイルが読む)

### [sim] セクション — シミュレーション基本

| パラメータ | 型 | デフォルト | 説明 |
|------------|------|-----------|------|
| `dt` | float | 1.0e-9 | タイムステップ [s] |
| `rng_seed` | int | 12345 | 乱数シード |
| `batch_count` | int | 1 | バッチ数 |
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
| `obj_scale` | float | 1.0 | スケーリング係数 |
| `obj_rotation` | float[3] | [0, 0, 0] | 回転角 [deg] (外因性 x→y→z) |
| `obj_offset` | float[3] | [0, 0, 0] | 平行移動 [m] |

**変換順序**: scale → rotate → offset

#### [[mesh.templates]] — 手続き的メッシュ生成

共通: `enabled` (bool), `kind` (enum), `center` (float[3])

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
| `mesh_triangles.csv` | CSV: `elem_idx, vertex_i, vertex_j, vertex_k, center_x, center_y, center_z, mesh_id` | 三角形メッシュ接続性 |
| `mesh_sources.csv` | CSV: `mesh_id, source_type, element_count` | メッシュソースメタデータ |
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
beachx config init [case.toml]                    # case.toml を新規作成
beachx config render [case.toml]                   # beach.toml にレンダリング
beachx config validate [case.toml]                 # バリデーション
beachx config diff left.toml right.toml            # 設定比較
beachx config save <name> [case.toml]              # テンプレートとして保存
beachx config list-saved                           # 保存済み一覧
```

### preset — プリセット管理

```bash
beachx preset list                                 # 全プリセット一覧
beachx preset show sim/periodic2_fmm               # 内容表示
beachx preset new sim/lab/custom                   # 新規作成
beachx preset new sim/lab/custom --from sim/periodic2_fmm  # クローン
beachx preset save <name> --section sim            # 設定から抽出
beachx preset path sim/periodic2_fmm               # ファイルパス表示
beachx preset edit sim/periodic2_fmm               # エディタで開く
beachx preset validate sim/periodic2_fmm           # バリデーション
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

## プリセットシステム

### 検索順序

1. `.beachx/presets/` (カレント〜親ディレクトリ)
2. `~/.config/beachx/presets/` (ユーザーレベル)
3. パッケージ同梱プリセット (組み込み)

### 命名規則

`category/subcategory/name` 形式。カテゴリ: `sim`, `species`, `mesh`, `output`

### マージルール

- `use_presets` リスト順に上から適用
- テーブル: ディープマージ
- スカラー: 後勝ち (last-one-wins)
- 配列: 追加 (`particles.species` と `mesh.templates` は ID によるディープマージ可)

### 組み込みプリセット

| プリセット | 概要 |
|-----------|------|
| `sim/periodic2_fmm` | 周期境界 + FMM ソルバー |
| `species/solarwind_electron` | 太陽風電子 (reservoir_face) |
| `species/solarwind_ion` | 太陽風イオン (reservoir_face) |
| `mesh/plane_basic` | 基本平面メッシュ (20x20) |
| `output/standard` | 標準出力設定 |

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

### 周期境界 + 太陽風 (case.toml)

```toml
use_presets = [
  "sim/periodic2_fmm",
  "species/solarwind_electron",
  "species/solarwind_ion",
  "mesh/plane_basic",
  "output/standard",
]

[sim]
batch_count = 200
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
│   ├── config/               #   設定合成システム
│   │   └── presets/          #   組み込みプリセット
│   └── fortran_results/      #   後処理・可視化
├── schemas/                  # JSON Schema (バリデーション用)
│   ├── beach.schema.json     #   beach.toml スキーマ
│   ├── beach.case.schema.json    # case.toml スキーマ
│   └── beach.preset.schema.json  # プリセットスキーマ
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
| `docs/config_workflow.md` | プリセット合成パターン |
| `docs/python_postprocess_api.md` | Python API リファレンス |
| `schemas/beach.schema.json` | IDE バリデーション用 JSON Schema |

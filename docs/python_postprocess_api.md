title: Python 後処理 API リファレンス

# Python 後処理 API リファレンス

BEACH の Python パッケージ (`beach`) は、Fortran シミュレーション結果の読み込み・解析・可視化を担う後処理レイヤです。
Fortran 実行系が出力するファイル群（`summary.txt`, `charges.csv`, `mesh_triangles.csv` 等）を読み込み、電位再構成、Coulomb 力計算、電場計算、電気力線追跡、3D 可視化を Python 側で行います。

## 1. パッケージ構成

| モジュール | 役割 |
|---|---|
| `beach.fortran_results.io` | 出力ディレクトリの読み込み (`load_fortran_result`, `list_fortran_runs`) |
| `beach.fortran_results.facade` | 高水準ファサード `Beach` クラス |
| `beach.fortran_results.potential` | 電位再構成 (`compute_potential_mesh`, `compute_potential_points`, `compute_potential_slices`) |
| `beach.fortran_results.coulomb` | Coulomb 力/トルク計算 (`calc_coulomb`) |
| `beach.fortran_results.field_lines` | 電場計算・電気力線追跡・3D 描画 (`compute_electric_field_points`, `trace_field_lines`, `plot_field_lines_3d`) |
| `beach.fortran_results.mobility` | Coulomb mobility 解析 (`analyze_coulomb_mobility`) |
| `beach.fortran_results.plotting` | 各種プロット (`plot_charge_mesh`, `plot_charges`, `plot_potential_mesh` 等) |
| `beach.fortran_results.animation` | 履歴アニメーション (`animate_history_mesh`) |
| `beach.fortran_results.history` | `charge_history.csv` のバッチステップ別アクセス (`FortranChargeHistory`) |
| `beach.fortran_results.types` | 公開データ型 (`FortranRunResult`, `CoulombInteraction` 等) |
| `beach.fortran_results.constants` | 物理定数 (`K_COULOMB`) |

すべての公開シンボルは `beach` トップレベルおよび `beach.fortran_results` からインポートできます。

```python
from beach import Beach, calc_coulomb, compute_electric_field_points, trace_field_lines
```

## 2. `Beach` ファサードクラス

出力ディレクトリを 1 つ束ねて、主要な解析・可視化メソッドを提供する高水準インターフェースです。

```python
b = Beach("outputs/latest")
```

### 2.1 コンストラクタ

| パラメータ | 型 | デフォルト | 説明 |
|---|---|---|---|
| `output_dir` | `str \| Path` | `"outputs/latest"` | Fortran 出力ディレクトリ |

### 2.2 プロパティ

| 名前 | 戻り型 | 説明 |
|---|---|---|
| `result` | `FortranRunResult` | 読み込み済み結果（遅延ロード）|
| `mesh_ids` | `tuple[int, ...]` | 利用可能な mesh ID 一覧 |

### 2.3 メソッド一覧

| メソッド | 委譲先 | 概要 |
|---|---|---|
| `reload()` | `load_fortran_result` | ディスクから再読み込み |
| `get_mesh(*mesh_ids, step)` | 内部 | mesh ID で `MeshSelection` を取得 |
| `get_mesh_charge(*mesh_ids, step)` | 内部 | mesh ID で要素電荷配列を取得 |
| `calc_coulomb(target, source, ...)` | `calc_coulomb` | Coulomb 力/トルク計算 |
| `analyze_coulomb_mobility(...)` | `analyze_coulomb_mobility` | オブジェクト別 mobility 解析 |
| `compute_potential(...)` | `compute_potential_mesh` | 重心での電位再構成 |
| `compute_potential_points(points, ...)` | `compute_potential_points` | 任意点での電位 |
| `compute_potential_slices(...)` | `compute_potential_slices` | XY/YZ/XZ 断面の電位 |
| `compute_electric_field(points, ...)` | `compute_electric_field_points` | 任意点での電場ベクトル |
| `trace_field_lines(seed_points, ...)` | `trace_field_lines` | 電気力線の RK4 追跡 |
| `plot_mesh(...)` | `plot_charge_mesh` | 電荷密度の 3D メッシュ描画 |
| `plot_potential(...)` | `plot_potential_mesh` | 電位の 3D メッシュ描画 |
| `plot_potential_slices(...)` | `plot_potential_slices` | 電位断面の描画 |
| `plot_field_lines(seed_points, ...)` | `plot_field_lines_3d` | 電気力線の 3D 描画 |
| `plot_bar()` | `plot_charges` | 要素電荷棒グラフ |
| `plot_mesh_source_boxplot(...)` | `plot_mesh_source_boxplot` | mesh source 別箱ひげ図 |
| `plot_coulomb_force_matrix(...)` | `plot_coulomb_force_matrix` | Coulomb 力行列プロット |
| `animate_mesh(...)` | `animate_history_mesh` | 電荷/電位履歴アニメーション |

## 3. 結果読み込み

### 3.1 `load_fortran_result(directory)`

Fortran 出力ディレクトリを読み込み、`FortranRunResult` を返します。

```python
from beach import load_fortran_result

result = load_fortran_result("outputs/latest")
print(f"要素数: {result.mesh_nelem}, バッチ数: {result.batches}")
print(f"吸収: {result.absorbed}, 脱出: {result.escaped}")
```

必須ファイル: `summary.txt`, `charges.csv`
オプションファイル: `mesh_triangles.csv`, `mesh_sources.csv`, `charge_history.csv`, `mesh_potential.csv`

### 3.2 `FortranRunResult` 型

| フィールド | 型 | 説明 |
|---|---|---|
| `directory` | `Path` | 出力ディレクトリパス |
| `mesh_nelem` | `int` | メッシュ要素数 |
| `processed_particles` | `int` | 処理済み粒子数 |
| `absorbed` | `int` | 吸収粒子数 |
| `escaped` | `int` | 脱出粒子数 |
| `batches` | `int` | 処理済みバッチ数 |
| `escaped_boundary` | `int` | 境界脱出粒子数 |
| `survived_max_step` | `int` | max_step 到達粒子数 |
| `last_rel_change` | `float` | 最終相対電荷変化量 |
| `charges` | `ndarray (mesh_nelem,)` | 要素電荷配列 [C] |
| `triangles` | `ndarray (mesh_nelem, 3, 3) \| None` | 三角形頂点座標 [m] |
| `mesh_ids` | `ndarray (mesh_nelem,) \| None` | 要素 mesh ID |
| `mesh_sources` | `dict[int, MeshSource] \| None` | mesh 種別メタデータ |
| `mesh_potential_v` | `ndarray (mesh_nelem,) \| None` | Fortran 出力の重心電位 [V] |
| `history` | `FortranChargeHistory \| None` | 電荷履歴アクセサ |

## 4. 電位再構成

### 4.1 `compute_potential_mesh(result, *, softening, self_term, periodic2, reference_point)`

三角形重心での電位を再構成します。Fortran が `mesh_potential.csv` を出力済みで条件が一致する場合はそちらを優先します。

| パラメータ | 型 | デフォルト | 単位 | 説明 |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (必須) | - | 結果オブジェクト |
| `softening` | `float \| None` | `None` | m | `None` で `sim.softening` を自動参照 |
| `self_term` | `str` | `"auto"` | - | 自己相互作用: `"auto"` / `"area_equivalent"` / `"exclude"` / `"softened_point"` |
| `periodic2` | `Mapping \| None` | `None` | - | 2 軸周期設定（後述）。`None` で自動判定 |
| `reference_point` | `Iterable[float] \| str \| None` | `None` | m | 基準電位点。`"species1_injection_center"` で species 1 注入面中心 |

戻り値: `ndarray (mesh_nelem,)` [V]

### 4.2 `compute_potential_points(result, points, *, softening, chunk_size, periodic2, reference_point)`

任意 3D 点での電位を計算します。

| パラメータ | 型 | デフォルト | 単位 | 説明 |
|---|---|---|---|---|
| `points` | `ndarray (n_points, 3)` | (必須) | m | サンプリング点座標 |
| `softening` | `float \| None` | `None` | m | `None` で自動 |
| `chunk_size` | `int` | `2048` | - | チャンク分割数 |
| `periodic2` | `Mapping \| None` | `None` | - | `None` で自動判定 |
| `reference_point` | `Iterable[float] \| str \| None` | `None` | m | 基準電位点 |

戻り値: `ndarray (n_points,)` [V]

### 4.3 `compute_potential_slices(result, *, box_min, box_max, grid_n, xy_z, yz_x, xz_y, ...)`

XY/YZ/XZ 平面上の電位断面を計算します。

戻り値: `dict[str, PotentialSlice2D]` (キー: `"xy"`, `"yz"`, `"xz"`)

## 5. Coulomb 力/トルク計算

### 5.1 `calc_coulomb(result, target, source, *, step, softening, torque_origin, periodic2)`

target メッシュグループが source から受ける Coulomb 力とトルクを計算します。

| パラメータ | 型 | デフォルト | 単位 | 説明 |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (必須) | - | 結果オブジェクト |
| `target` | `int \| MeshSelection \| Iterable` | (必須) | - | ターゲットメッシュグループ (group A) |
| `source` | `int \| MeshSelection \| Iterable` | (必須) | - | ソースメッシュグループ (group B) |
| `step` | `int \| None` | `-1` | - | 履歴ステップ。`-1` で最新、`None` で最終電荷 |
| `softening` | `float` | `0.0` | m | ソフトニング長 |
| `torque_origin` | `str` | `"target_center"` | - | トルク基準点: `"target_center"` / `"source_center"` / `"origin"` |
| `periodic2` | `Mapping \| None` | `None` | - | 2 軸周期境界設定。`None` で自動判定 |

戻り値: `CoulombInteraction`

### periodic2 パラメータによる周期クーロン和

`periodic2` が指定された場合、ソース電荷の画像シェル `ix in [-nimg, nimg], iy in [-nimg, nimg]` を 2 軸周期方向に生成し、最近接セル和としてクーロン力/トルクを計算します。

`periodic2=None`（デフォルト）の場合は、出力ディレクトリ近傍の `beach.toml` を探索し、`sim.field_bc_mode="periodic2"` が設定されていれば自動的に周期設定を適用します。これは `compute_potential_mesh` 等の他の関数と同じ自動判定ロジックです。

### 5.2 `CoulombInteraction` 型

| フィールド | 型 | 単位 | 説明 |
|---|---|---|---|
| `group_a_mesh_ids` | `tuple[int, ...]` | - | ターゲット mesh ID |
| `group_b_mesh_ids` | `tuple[int, ...]` | - | ソース mesh ID |
| `step` | `int \| None` | - | 使用した履歴ステップ |
| `softening` | `float` | m | 使用したソフトニング長 |
| `torque_origin_m` | `ndarray (3,)` | m | トルク基準点 |
| `force_on_a_N` | `ndarray (3,)` | N | group A に作用する正味力 |
| `force_on_b_N` | `ndarray (3,)` | N | group B に作用する正味力 |
| `torque_on_a_Nm` | `ndarray (3,)` | N m | group A に作用する正味トルク |
| `torque_on_b_Nm` | `ndarray (3,)` | N m | group B に作用する正味トルク |
| `mean_force_on_a_per_element_N` | `ndarray (3,)` | N | ターゲット要素あたり平均力 |
| `mean_torque_on_a_per_element_Nm` | `ndarray (3,)` | N m | ターゲット要素あたり平均トルク |

### 使用例

```python
from beach import Beach

b = Beach("outputs/latest")

# mesh_id=0 が mesh_id=1 から受ける Coulomb 力
interaction = b.calc_coulomb(target=0, source=1)
print(f"Force on target: {interaction.force_on_a_N} [N]")
print(f"Torque on target: {interaction.torque_on_a_Nm} [N m]")

# periodic2 を明示指定する場合
interaction_p = b.calc_coulomb(
    target=0, source=1,
    periodic2={"axes": [0, 1], "lengths": [0.01, 0.01], "image_layers": 2},
)
```

## 6. 電場計算

### 6.1 `compute_electric_field_points(result, points, *, softening, chunk_size, periodic2)`

任意 3D 点での電場ベクトルを、表面電荷からクーロン則で直接計算します。

計算式: `E(r) = K * sum_j q_j * (r - r_j) / |r - r_j|^3`

ここで `r_j` は三角形要素 `j` の重心、`q_j` は要素電荷です。

| パラメータ | 型 | デフォルト | 単位 | 説明 |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (必須) | - | 結果オブジェクト |
| `points` | `ndarray (n_points, 3)` | (必須) | m | サンプリング点座標 |
| `softening` | `float \| None` | `None` | m | ソフトニング長。`None` で `sim.softening` を自動参照 |
| `chunk_size` | `int` | `2048` | - | チャンクサイズ |
| `periodic2` | `Mapping \| None` | `None` | - | 2 軸周期設定。`None` で自動判定 |

戻り値: `ndarray (n_points, 3)` [V/m]

periodic2 モードでは、サンプリング点を周期セルに wrap した上で、ソース電荷の画像シェルを `ix in [-nimg, nimg], iy in [-nimg, nimg]` で重畳します。

### 使用例

```python
import numpy as np
from beach import Beach

b = Beach("outputs/latest")

# グリッド点での電場計算
x = np.linspace(0.0, 0.01, 50)
y = np.linspace(0.0, 0.01, 50)
xx, yy = np.meshgrid(x, y)
zz = np.full_like(xx, 0.005)
points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])

efield = b.compute_electric_field(points)
print(f"Shape: {efield.shape}")   # (2500, 3)
print(f"E-field [V/m]: {efield[0]}")
```

## 7. 電気力線追跡

### 7.1 `trace_field_lines(result, seed_points, *, ds, max_steps, softening, periodic2, direction, box_min, box_max)`

シード点から電場方向（または逆方向）に RK4 積分で電気力線を追跡します。

| パラメータ | 型 | デフォルト | 単位 | 説明 |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (必須) | - | 結果オブジェクト |
| `seed_points` | `ndarray (n_seeds, 3)` | (必須) | m | 力線の開始点座標 |
| `ds` | `float \| None` | `None` | m | 積分ステップサイズ。`None` でメッシュ平均辺長 x 0.5 から自動設定 |
| `max_steps` | `int` | `500` | - | 各方向の最大積分ステップ数 |
| `softening` | `float \| None` | `None` | m | ソフトニング長。`None` で自動 |
| `periodic2` | `Mapping \| None` | `None` | - | 2 軸周期設定。`None` で自動判定 |
| `direction` | `str` | `"both"` | - | 追跡方向: `"forward"` (電場方向) / `"backward"` (逆方向) / `"both"` (両方向) |
| `box_min` | `Iterable[float] \| None` | `None` | m | 境界ボックス下限。力線がこの外に出たら打ち切り |
| `box_max` | `Iterable[float] \| None` | `None` | m | 境界ボックス上限 |

戻り値: `list[ndarray]` -- 各要素は shape `(n_points_i, 3)` の力線座標 [m]

#### 追跡アルゴリズム

- 4 次 Runge-Kutta 法 (RK4) で電場の単位ベクトル方向に `ds` ずつ進行
- 各 RK4 ステージで電場ノルムが `1e-30` 未満になったら打ち切り
- `direction="both"` の場合は forward と backward を接合（シード点の重複を除去）
- `box_min` / `box_max` を超えた時点で打ち切り

#### 制約事項

- Python 側での直接和計算であり、大規模メッシュ（数万要素以上）では計算時間が長くなります
- Fortran の treecode/fmm は使用しません。Python 側の電場計算は要素重心の点電荷による direct 和のみです
- periodic2 では explicit image shell のみを再構成し、oracle residual（Ewald 遠方補正）は Python 側では再現しません

### 使用例

```python
import numpy as np
from beach import Beach

b = Beach("outputs/latest")

# シード点を手動指定
seeds = np.array([
    [0.005, 0.005, 0.008],
    [0.005, 0.005, 0.002],
])

lines = b.trace_field_lines(seeds, max_steps=1000)
print(f"力線数: {len(lines)}")
for i, line in enumerate(lines):
    print(f"  力線 {i}: {line.shape[0]} 点")
```

### 7.2 `plot_field_lines_3d(result, seed_points, *, ...)`

電気力線を 3D 描画し、オプションでメッシュ表面を電荷密度で着色してオーバーレイします。

| パラメータ | 型 | デフォルト | 説明 |
|---|---|---|---|
| `result` | `FortranRunResult \| object` | (必須) | 結果オブジェクト |
| `seed_points` | `ndarray (n_seeds, 3)` | (必須) | シード点 [m] |
| `ds` | `float \| None` | `None` | 積分ステップサイズ [m] |
| `max_steps` | `int` | `500` | 最大ステップ数 |
| `softening` | `float \| None` | `None` | ソフトニング長 [m] |
| `periodic2` | `Mapping \| None` | `None` | 周期設定 |
| `direction` | `str` | `"both"` | 追跡方向 |
| `box_min` | `Iterable[float] \| None` | `None` | 境界ボックス下限 [m] |
| `box_max` | `Iterable[float] \| None` | `None` | 境界ボックス上限 [m] |
| `show_mesh` | `bool` | `True` | メッシュオーバーレイの表示 |
| `mesh_alpha` | `float` | `0.25` | メッシュ面の透明度 |
| `mesh_cmap` | `str` | `"coolwarm"` | メッシュの面電荷密度カラーマップ |
| `line_color` | `str \| None` | `None` | 力線の固定色。`None` で `line_cmap` による色分け |
| `line_cmap` | `str` | `"plasma"` | 力線のカラーマップ（`line_color=None` 時） |
| `line_width` | `float` | `1.2` | 力線の線幅 |
| `view_elev` | `float` | `24.0` | 仰角 [deg] |
| `view_azim` | `float` | `-58.0` | 方位角 [deg] |
| `title` | `str` | `"Electric field lines"` | プロットタイトル |
| `figsize` | `tuple[float, float]` | `(9, 7)` | Figure サイズ [inch] |

戻り値: `(figure, axes)` -- matplotlib の Figure / Axes3D

#### 描画内容

- 各力線は線として描画。`line_color=None` の場合は力線ごとに `line_cmap` で色分け
- 力線の中間点に方向矢印（quiver）を描画
- シード点を赤色の散布点として描画
- `show_mesh=True` の場合、三角形メッシュを面電荷密度 `q / A` で着色し半透明で重ねる

### 使用例

```python
import numpy as np
from beach import Beach

b = Beach("outputs/latest")

seeds = np.array([
    [0.005, 0.005, 0.008],
    [0.003, 0.007, 0.006],
    [0.007, 0.003, 0.006],
])

fig, ax = b.plot_field_lines(
    seeds,
    max_steps=800,
    direction="both",
    show_mesh=True,
    mesh_alpha=0.3,
    line_cmap="viridis",
    view_elev=30,
    view_azim=-45,
)
fig.savefig("field_lines.png", dpi=150)
```

## 8. `periodic2` パラメータ仕様

電位再構成、Coulomb 力計算、電場計算、電気力線追跡のすべてで共通の `periodic2` パラメータが利用できます。

### 8.1 自動判定（`periodic2=None`）

デフォルトの `None` では、出力ディレクトリ近傍の `beach.toml` を探索し、`sim.field_bc_mode="periodic2"` が設定されている場合に自動的に周期境界設定を適用します。
設定ファイルが見つからない場合や `field_bc_mode` が `periodic2` でない場合は自由空間モードで計算します。

### 8.2 明示指定

`periodic2` は以下のキーを持つ `Mapping` で指定します。

| キー | 型 | 必須 | デフォルト | 説明 |
|---|---|---|---|---|
| `axes` | `list[int]` (長さ 2) | 必須 | - | 周期軸の 0-based インデックス (例: `[0, 1]` は x, y 軸) |
| `lengths` | `list[float]` (長さ 2) | 必須 | - | 各周期軸のボックス長 [m]。正の値 |
| `origins` | `list[float]` (長さ 2) | - | `[0.0, 0.0]` | 各周期軸のボックス原点 [m] |
| `box_min` | `list[float]` (長さ 3) | - | - | `origins` の代替。3D ボックス下限から周期軸の原点を抽出 |
| `image_layers` | `int` | - | `1` | 画像シェルの層数。各周期軸で `[-N, N]` を評価 |
| `far_correction` | `str` | - | `"auto"` | `"auto"` / `"none"` / `"m2l_root_oracle"`。Python 側は設定互換のために保持するが、oracle residual 自体は再現しない |
| `ewald_alpha` | `float` | - | `0.0` | Ewald 分解パラメータ（予約） |
| `ewald_layers` | `int` | - | `4` | Ewald 打切り深さ（予約） |

`origins` と `box_min` が両方指定された場合は `origins` が優先されます。

### 使用例

```python
p2 = {
    "axes": [0, 1],
    "lengths": [0.01, 0.01],
    "image_layers": 2,
}
potential = b.compute_potential(periodic2=p2)
efield = b.compute_electric_field(points, periodic2=p2)
interaction = b.calc_coulomb(target=0, source=1, periodic2=p2)
lines = b.trace_field_lines(seeds, periodic2=p2)
```

### 8.3 Python 側の periodic2 実装の制限

- Python 側では explicit image shell による直接和のみで周期和を再構成します
- Fortran 側 FMM の `m2l_root_oracle` による Ewald 遠方補正は Python 側では再現されません。`far_correction` は設定の互換性のために保持されますが、計算に影響しません
- 大きな `image_layers` を指定するほど精度は向上しますが、計算量は `(2*N+1)^2` 倍に増加します

## 9. Coulomb mobility 解析

### 9.1 `analyze_coulomb_mobility(result, *, step, softening, config_path, gravity, support_normal, ...)`

オブジェクト単位で Coulomb 力による滑り・転がり・浮上の傾向を解析します。

戻り値: `CoulombMobilityAnalysis` (`.records` に `CoulombMobilityRecord` のタプルを格納)

## 10. 可視化関数

### 10.1 電荷/電位メッシュ描画

```python
fig, ax = b.plot_mesh(cmap="coolwarm")                        # 電荷密度
fig, ax = b.plot_potential(reference_point="species1_injection_center")  # 電位
```

### 10.2 電位断面

```python
fig, axes = b.plot_potential_slices(
    box_min=[0, 0, 0], box_max=[0.01, 0.01, 0.01],
    xy_z=0.005,
)
```

### 10.3 履歴アニメーション

```python
gif_path = b.animate_mesh("charge_animation.gif", quantity="charge", fps=10)
```

### 10.4 Coulomb 力行列

```python
fig, ax = b.plot_coulomb_force_matrix(component="z")
```

## 11. CLI コマンド

| コマンド | 説明 |
|---|---|
| `beach-inspect <output_dir>` | 出力ディレクトリの要約表示 |
| `beach-animate-history <output_dir>` | 電荷/電位履歴のアニメーション GIF 生成 |
| `beach-estimate-workload <config.toml>` | ワークロード見積もり |
| `beach-plot-potential-slices <output_dir>` | 電位断面の描画 |
| `beach-plot-performance-profile <output_dir>` | パフォーマンスプロファイルの描画 |
| `beach-plot-coulomb-force-matrix <output_dir>` | Coulomb 力行列の描画 |

## 12. 物理定数

| シンボル | 値 | 単位 | 説明 |
|---|---|---|---|
| `K_COULOMB` | `8.9875517923e9` | N m^2 / C^2 | クーロン定数 |

title: Python 後処理 API リファレンス

Lang: [日本語](PythonPostprocessAPI.md) | [English](PythonPostprocessAPI.en.md)

# Python 後処理 API リファレンス

BEACHのPython package (`beach`) は、Fortranシミュレーションの結果を読み込み、解析・可視化するための後処理レイヤです。
`summary.txt`、`charges.csv`、`mesh_triangles.csv`などを読み込み、電位の再構成、Coulomb力と電場の計算、
電気力線の追跡、3D可視化を行います。

最初の図を作る手順は [後処理チュートリアル](PostprocessTutorial.html) にまとめています。

## 1. パッケージ構成

| モジュール | 役割 |
|---|---|
| `beach.fortran_results.io` | 出力ディレクトリの読み込み (`load_fortran_result`, `list_fortran_runs`) |
| `beach.fortran_results.facade` | 高水準ファサード `Beach` クラス |
| `beach.fortran_results.potential` | 電位再構成 (`compute_potential_mesh`, `compute_potential_points`, `compute_potential_slices`) |
| `beach.fortran_results.coulomb` | Coulomb 力/トルク計算 (`calc_coulomb`) |
| `beach.fortran_results.kernel` | Fortran FMM field kernel の共有ライブラリ呼び出し (`FieldKernel`, `calc_object_forces_kernel`) |
| `beach.fortran_results.object_interaction` | 凍結した source 電荷に対する object の力・トルク・鉛直経路 (`ObjectInteractionSnapshot`, `ObjectProbe`) |
| `beach.fortran_results.detachment` | 経路仕事、付着、重力、速度、from-rest barrier の immutable 結果型 |
| `beach.fortran_results.periodic_force_oracle` | 有限周期画像 shell と `E_bottom=0` 境界条件の収束診断 |
| `beach.fortran_results.scene` | object の一時移動・回転と編集後 scene の field-kernel 評価 |
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

`Beach("outputs/latest", config_path="path/to/beach.toml")`のように、設定ファイルを明示できます。
`config_path=None`の場合は、`output_dir/beach.toml`、親ディレクトリ、祖父ディレクトリの順に自動探索します。
objectやkernelの設定を参照する解析では、この`beach.toml`からobject kind/order、periodic2、treeパラメータを取得します。
全電場を再計算する高水準APIは、境界・一様場・外部場を推測しないため、この設定ファイルを必須とします。
simulatorが保存した`mesh_potential.csv`をそのまま読む経路では再計算しないため不要です。

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
| `calc_object_forces_kernel(...)` | `calc_object_forces_kernel` | Fortran field kernel による object 別合力計算 |
| `object_interaction_snapshot(...)` | `ObjectInteractionSnapshot.from_result` | 凍結 source に対する primary-only self exclusion の力・離脱経路解析 |
| `scene(step, ...)` | `BeachScene.from_result` | object を一時的に移動・回転する what-if scene |
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
| `mesh_sources` | `dict[int, MeshSource] \| None` | mesh 種別・surface model・epsilon_r メタデータ |
| `mesh_potential_v` | `ndarray (mesh_nelem,) \| None` | Fortran 出力の重心電位 [V] |
| `history` | `FortranChargeHistory \| None` | 電荷履歴アクセサ |

### 3.3 `FortranChargeHistory`

Fortranの`charge_history.csv`は、batchごとに全要素を記録するdense snapshotです。
各batchには、`elem_idx=1..mesh_nelem`の要素が1回ずつ含まれている必要があります。
要素の行順は問いませんが、batch groupは厳密な昇順でなければなりません。また、各batch内の
`processed_particles`と`rel_change`は全行で一致し、`charge_C`と`rel_change`は有限値である必要があります。

`load_fortran_result(...)`と履歴accessorの構築時には、履歴全体をすぐに読み込みません。
`batch_indices`などの履歴property、`get_step(...)`、`as_array()`のいずれかを最初に参照したときに、
CSV全体のbyte-offset indexを作りながら全batchを検証します。
欠損、重複、範囲外 index、非有限値、metadata 不一致、batch の逆行がある場合は、
batch と破損内容を含む `ValueError` を送出します。欠損要素を物理電荷 `0 C`
として補完することはありません。

検証後も、各batchの電荷vectorは要求されたときに読み込みます。`as_array()`を呼ぶまで、
履歴全体のmatrixは作りません。`FortranChargeHistory.from_arrays(...)`は、呼び出し元がdense matrixを渡す
trusted in-memory経路として動作します。

## 4. 電位再構成

### 4.1 `compute_potential_mesh(result, *, periodic2, reference_point, config_path, library_path)`

native field kernelのP0 triangle panel評価で、三角形重心の電位を計算します。Fortran が
`mesh_potential.csv` を出力済みで条件が一致する場合はそちらを優先します。

| パラメータ | 型 | デフォルト | 単位 | 説明 |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (必須) | - | 結果オブジェクト |
| `periodic2` | `Mapping \| None` | `None` | - | 2 軸周期設定（後述）。`None` で自動判定 |
| `reference_point` | `Iterable[float] \| str \| None` | `None` | m | 基準電位点。`"species1_injection_center"` で species 1 注入面中心 |
| `config_path` | `str \| Path \| None` | `None` | - | `beach.toml`の明示パス |
| `library_path` | `str \| Path \| None` | `None` | - | native field kernel共有libraryの明示パス |

戻り値: `ndarray (mesh_nelem,)` [V]

### 4.2 `compute_potential_points(result, points, *, chunk_size, periodic2, reference_point, config_path, library_path)`

任意 3D 点での電位を計算します。

| パラメータ | 型 | デフォルト | 単位 | 説明 |
|---|---|---|---|---|
| `points` | `ndarray (n_points, 3)` | (必須) | m | サンプリング点座標 |
| `chunk_size` | `int` | `2048` | - | チャンク分割数 |
| `periodic2` | `Mapping \| None` | `None` | - | `None` で自動判定 |
| `reference_point` | `Iterable[float] \| str \| None` | `None` | m | 基準電位点 |
| `config_path` | `str \| Path \| None` | `None` | - | `beach.toml`の明示パス |
| `library_path` | `str \| Path \| None` | `None` | - | native field kernel共有libraryの明示パス |

戻り値: `ndarray (n_points,)` [V]

### 4.3 `compute_potential_slices(result, *, box_min, box_max, grid_n, xy_z, yz_x, xz_y, ...)`

XY/YZ/XZ 平面上の電位断面を計算します。

戻り値: `dict[str, PotentialSlice2D]` (キー: `"xy"`, `"yz"`, `"xz"`)

## 5. Coulomb 力/トルク計算

### 5.1 `calc_coulomb(result, target, source, *, step, torque_origin, periodic2, quadrature_order, config_path, library_path)`

target メッシュグループが source から受ける Coulomb 力とトルクを計算します。

| パラメータ | 型 | デフォルト | 単位 | 説明 |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (必須) | - | 結果オブジェクト |
| `target` | `int \| MeshSelection \| Iterable` | (必須) | - | ターゲットメッシュグループ (group A) |
| `source` | `int \| MeshSelection \| Iterable` | (必須) | - | ソースメッシュグループ (group B) |
| `step` | `int \| None` | `-1` | - | 履歴ステップ。`-1` で最新、`None` で最終電荷 |
| `torque_origin` | `str` | `"target_center"` | - | トルク基準点: `"target_center"` / `"source_center"` / `"origin"` |
| `periodic2` | `Mapping \| None` | `None` | - | 2 軸周期境界設定。`None` で自動判定 |
| `quadrature_order` | `int` | `7` | - | target panelのGauss-Duffy積分次数 |
| `config_path` | `str \| Path \| None` | `None` | - | `beach.toml`の明示パス |
| `library_path` | `str \| Path \| None` | `None` | - | native field kernel共有libraryの明示パス |

戻り値: `CoulombInteraction`

### periodic2 パラメータによる周期クーロン和

`periodic2` が指定された場合、native field kernelをその2軸周期設定で構成して力とトルクを計算します。

`periodic2=None`（デフォルト）の場合は、出力ディレクトリ近傍の `beach.toml` を探索し、`sim.field_bc_mode="periodic2"` が設定されていれば自動的に周期設定を適用します。これは `compute_potential_mesh` 等の他の関数と同じ自動判定ロジックです。

### 5.2 `CoulombInteraction` 型

| フィールド | 型 | 単位 | 説明 |
|---|---|---|---|
| `group_a_mesh_ids` | `tuple[int, ...]` | - | ターゲット mesh ID |
| `group_b_mesh_ids` | `tuple[int, ...]` | - | ソース mesh ID |
| `step` | `int \| None` | - | 使用した履歴ステップ |
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

### 6.1 `compute_electric_field_points(result, points, *, chunk_size, periodic2, config_path, library_path)`

任意 3D 点での電場ベクトルを、native field kernelのP0 triangle panel評価で計算します。

| パラメータ | 型 | デフォルト | 単位 | 説明 |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (必須) | - | 結果オブジェクト |
| `points` | `ndarray (n_points, 3)` | (必須) | m | サンプリング点座標 |
| `chunk_size` | `int` | `2048` | - | チャンクサイズ |
| `periodic2` | `Mapping \| None` | `None` | - | 2 軸周期設定。`None` で自動判定 |
| `config_path` | `str \| Path \| None` | `None` | - | `beach.toml`の明示パス |
| `library_path` | `str \| Path \| None` | `None` | - | native field kernel共有libraryの明示パス |

戻り値: `ndarray (n_points, 3)` [V/m]

periodic2 モードでは、設定に応じて有限画像またはcached非零modeをnative kernelが評価します。

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

### 7.1 `trace_field_lines(result, seed_points, *, ds, max_steps, periodic2, direction, box_min, box_max, config_path, library_path)`

シード点から電場方向（または逆方向）に RK4 積分で電気力線を追跡します。

| パラメータ | 型 | デフォルト | 単位 | 説明 |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (必須) | - | 結果オブジェクト |
| `seed_points` | `ndarray (n_seeds, 3)` | (必須) | m | 力線の開始点座標 |
| `ds` | `float \| None` | `None` | m | 積分ステップサイズ。`None` でメッシュ平均辺長 x 0.5 から自動設定 |
| `max_steps` | `int` | `500` | - | 各方向の最大積分ステップ数 |
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

- native field kernel共有libraryが必要です
- RK4の各stageでfield kernelを評価するため、seed数と`max_steps`に比例して計算時間が増えます
- mesh衝突による自動停止は行わず、field閾値、`max_steps`、box外への退出で停止します

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
設定ファイルが見つからない場合、全電場を再計算する高水準APIは自由空間へ黙ってfallbackせず停止します。
設定に`field_bc_mode="free"`が明記されている場合だけ自由空間として再計算します。

### 8.2 明示指定

`periodic2` は以下のキーを持つ `Mapping` で指定します。

| キー | 型 | 必須 | デフォルト | 説明 |
|---|---|---|---|---|
| `axes` | `list[int]` (長さ 2) | 必須 | - | 周期軸の 0-based インデックス (例: `[0, 1]` は x, y 軸) |
| `lengths` | `list[float]` (長さ 2) | 必須 | - | 各周期軸のボックス長 [m]。正の値 |
| `origins` | `list[float]` (長さ 2) | - | `[0.0, 0.0]` | 各周期軸のボックス原点 [m] |
| `box_min` | `list[float]` (長さ 3) | - | - | `origins` の代替。3D ボックス下限から周期軸の原点を抽出 |
| `image_layers` | `int` | - | `1` | 画像シェルの層数。各周期軸で `[-N, N]` を評価 |
| `far_correction` | `str` | - | `"none"` | 明示指定では `"auto"` / `"none"`。`auto` は互換用に `"none"` として扱う |
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
- `cached_kneq0` native operatorは$k\ne0$成分だけを返します。一般の電位・電場・力APIはこれを全電場として返さず停止します。保存済み`mesh_potential.csv`、または物理的zero modeを明示合成する`ObjectInteractionSnapshot`を使ってください
- 削除済み `m2l_root_oracle` は、`periodic2=None` で出力近傍の過去メタデータを自動検出する場合だけ受理され、有限 image shell を表す `none` に正規化されます。明示した `periodic2` や設定ファイルでは拒否されます。Python 側で Ewald 遠方補正は再現しません
- 大きな `image_layers` を指定するほど精度は向上しますが、計算量は `(2*N+1)^2` 倍に増加します

## 9. Coulomb mobility 解析

### 9.1 `analyze_coulomb_mobility(result, *, step, config_path, library_path, gravity, support_normal, ...)`

オブジェクト単位で Coulomb 力による滑り・転がり・浮上の傾向を解析します。

戻り値: `CoulombMobilityAnalysis` (`.records` に `CoulombMobilityRecord` のタプルを格納)

## 10. Fortran field kernel 連携

### 10.1 `FieldKernel`

`make build-kernel` で生成した `build/libbeach_field_kernel.so` を `ctypes` 経由で読み込み、シミュレーションと同じFortran field coreでP0 triangle panelの電場・電位を評価します。`config_path` または自動探索された `beach.toml` から、periodic2、tree 設定を読みます。

```python
from beach import Beach, FieldKernel

run = Beach("outputs/latest")
with FieldKernel.from_result(run) as kernel:
    e = kernel.eval_e([[0.0, 0.0, 0.01]])
    phi = kernel.eval_phi([[0.0, 0.0, 0.01]])
```

共有ライブラリを別パスに置く場合は `library_path=` または環境変数 `BEACH_FIELD_KERNEL_LIB` を指定します。

`eval_e()` / `eval_phi()`は、構築時に選んだfree / finite periodic / cached periodic設定を使います。
ただしcached planの返り値は`cached_kneq0`の$k\ne0$成分と`sim.e0`だけで、物理的zero modeや
active outer stateを含むsimulator全電場ではありません。`FieldKernel`はこの成分を明示的に扱う低水準APIです。
`eval_e_direct()` / `eval_phi_direct()`は、同じsource geometryとchargeを非周期のexact direct法で評価する診断APIです。
periodic planでは使用できず、uniform `sim.e0`も加算しません。FMMの精度確認とprimary free-space subtractionの
小規模oracleに使います。

### 10.2 `calc_object_forces_kernel(result, ...)`

各objectについて、自身のsource電荷を0にしたうえで、`sum(q_i E_not_self(r_i))`とトルクを計算します。
self policyは`exclude_target_lattice`です。周期計算では、targetのprimary sourceと、target object自身の周期画像を除外します。
この高水準APIは未合成の`cached_kneq0`を拒否します。

object自身の周期画像を保持したまま離脱力を評価する場合は、後述の
`ObjectInteractionSnapshot.object_probe(...)`を使います。このAPIでは、self policyが
`exclude_primary_keep_images`に固定されています。

```python
from beach import Beach

run = Beach("outputs/latest")
records = run.calc_object_forces_kernel()
for record in records:
    print(record.mesh_id, record.total_charge_C, record.force_N)
```

### 10.3 `BeachScene`

`Beach.scene()`は、出力済みの帯電meshをPython側で一時的に編集するwhat-if viewです。
`move` / `rotate`は、元のsceneを変更せずに新しいsceneを返します。各要素の電荷はその要素に保持したまま、
重心と頂点だけを剛体変換します。編集後に`calc_object_forces_kernel`を呼ぶと、変換後のgeometryを
source/targetとしてFortran field kernelに渡します。

```python
from beach import Beach

run = Beach("outputs/latest", config_path="examples/beach.toml")
scene = run.scene()
moved = scene.move(2, by=[1.0e-3, 0.0, 0.0]).rotate(
    2,
    axis=[0.0, 0.0, 1.0],
    angle_deg=15.0,
)
records = moved.calc_object_forces_kernel(target_mesh_ids=[2])
print(records[0].force_N, records[0].torque_Nm)
```

Python側の剛体変換は、既定でNumPyを使います。Numbaを使う場合は、任意依存の`pip install ".[accel]"`を導入し、
`run.scene(transform_backend="numba")`を指定します。FMM、periodic2、遠方補正を含むfield評価はFortran kernel側で実行します。

### 10.4 `ObjectInteractionSnapshot` と凍結 source 経路

このAPIは、保存済みの全source geometryとchargeを一度凍結します。その凍結したsourceに対して、
選択したcentral-cell objectだけを独立したtarget probeとして動かします。self policyは`exclude_primary_keep_images`に固定され、
central-cell targetのprimary free-space寄与だけを差し引きます。target自身の周期画像、他objectの全画像、uniform fieldは残します。

```python
import numpy as np
from beach import AdhesionProfile, Beach

run = Beach("outputs/latest", config_path="beach.toml")
with run.object_interaction_snapshot(
    periodic_model="infinite_physical",
    cache_dir=".beach_cache/periodic2",
) as snapshot:
    probe = snapshot.object_probe(6)
    wrench = probe.wrench()
    path = probe.vertical_path(np.linspace(0.0, 2.0e-4, 65))

release = path.evaluate_release(
    mass_kg=2.0e-12,
    gravity_m_s2=9.80665,
    adhesion=AdhesionProfile.finite_range_constant(
        force_N=1.0e-10,
        range_m=2.0e-6,
    ),
)
print(wrench.force_N, wrench.torque_Nm)
print(path.status, path.work_relative_mismatch)
print(release.barrier_free_from_rest, release.endpoint_speed_m_s)
```

`periodic_model` の意味は次のとおりです。

| 値 | 場の定義 |
|---|---|
| `"configured"` | run の `beach.toml` をそのまま使用する。free、`far_correction="none"` の有限画像、または cached periodic のいずれにもなり得る |
| `"infinite_physical"` | x/y periodic run に対して cached `k != 0` と、`[periodic2].lower_boundary_model`に従う物理的な`k = 0` modeを組み合わせる。cache の生成・再利用条件を満たす必要がある |

完全な`beach.toml`と`sim.box_min` / `sim.box_max`が必要です。
`external_boundary.field.model`は現行仕様どおり`"none"`でなければなりません。

x/y periodic meshがcell seamをまたぐ場合、snapshot全体のsource geometryは、simulationとcache identityに一致する
saved表現のまま保持します。一方、`object_probe()`は選択したmeshだけを周期的に連結したbranchにunwrapします。
このtarget geometryをquadrature、central primaryの除外、剛体変換、面積重心、bounding radiusの計算で共通に使います。
`probe.target_geometry_representation`、
`target_triangles_m`、`geometric_area_centroid_m`、`vertex_bounding_center_m`、
`vertex_bounding_radius_m` でこの幾何を監査できます。

`ObjectWrench.components` は次の物理成分を保持します。

| key | 内容 |
|---|---|
| `other_objects_all_images` | target 以外の object とその周期画像 |
| `target_periodic_images` | target 自身の周期画像。central primary は除外済み |
| `external_uniform` | 設定された一様外部電場 |
| `total_external` | 上の3成分の和で、`ObjectWrench.force_N` / `torque_Nm` に一致 |

`numerical_metadata`のkernel/zero-mode内訳は、同じtotalを数値的に分解したものです。
torqueは基準点に依存するため、`ObjectWrench.torque_origin_m`と
`numerical_metadata["torque_origin_policy"]`をforceと一緒に保存してください。
`vertical_path()`では、各高さで使った基準点を`numerical_metadata["torque_origin_m"]`に保存します。

target integration は既定でorder 7のGauss-Duffy面積積分です。
`target_integration="centroid_compatibility"` は過去の重心積分結果との比較用で、
`auto`には選ばれません。

surface 上の機械的な合力には zero-mode の principal-value (PV) trace を使います。
simulator の粒子 pusher は表面の片側値 `zero_mode_trace_plus` を使うため、両者は目的が
異なります。cached 結果の `numerical_metadata["cached_kneq0_trace_correction"]` は、
native trace を PV 分解へ写すために **すでに `periodic_kneq0` へ含めた診断値** です。
この値を `force_N` や `periodic_kneq0` に再加算してはいけません。

`vertical_path()`では、source geometry/chargeを初期位置に固定し、target quadratureだけを`+z`方向に平行移動します。
全triangle頂点が非周期方向のbox/interface内にある経路だけを評価します。凍結した外部場に対するpotential energyは
`U_env = sum_i(q_i phi_env(r_i))`であり、係数`1/2`は付けません。

返り値では、
`electrostatic_work_J = integral(F_z dh)` と
`potential_difference_work_J = U_env(0)-U_env(h)` の不一致、適応細分化、
`status`を必ず確認してください。

`AdhesionProfile.finite_range_constant(F, d)`は、`0 <= h < d`で働く抵抗力と、`F*d`で飽和する付着仕事を表します。
`evaluate_release()`は、静電力の仕事から重力、付着、任意の散逸を差し引き、連続する区分経路全体を調べます。
`endpoint_positive=True`は終点での利用可能エネルギーのみを示し、途中のbarrierを通過できるかどうかは示しません。
静止状態から開始する場合は、`barrier_free_from_rest`と`first_inaccessible_displacement_m`で判定します。
速度配列は、利用可能エネルギーが負になる区間を0にclampします。

非中性の x/y 無限周期 cell では `E_bottom=0` zero mode により遠方の一定場と線形電位が
残り得ます。その場合、有限 box 内の endpoint work/speed は計算できますが、有限高さの
結果を無限遠への escape energy や終端速度と解釈できません。

### 10.5 有限画像 shell oracle

`finite_shell_wrench()`は、Fortran native finite-image kernel (`far_correction="none"`) で次の2つを評価します。

- raw symmetric shell
- 解析的な`Q_cell/(2 epsilon0 A_xy) e_z`を加えた`E_bottom=0`境界条件

返り値の`selected`は、要求した境界条件です。`symmetric`と`e_bottom_zero`も結果recordに残るため、境界条件間の差を確認できます。
sourceにはnative canonical-unwrapped表現を使います。Python側で周期画像を追加生成したり、移動後のtargetをprimary cellにwrapしたりしません。

```python
from beach import finite_shell_convergence, finite_shell_wrench

row = finite_shell_wrench(
    snapshot, probe, transform=None, image_layers=1, closure="e_bottom_zero"
)
shells = finite_shell_convergence(
    snapshot, probe, path.displacement_m, max_layers=12
)
```

`finite_shell_convergence()`は、`force_tail_proxy_N`と`work_tail_proxy_J`を使ってshell収束を判定します。
隣接shell間の増分が2回連続で小さい場合でもfalse convergenceが確認されたため、この増分だけでは判定しません。
snapshotが`infinite_physical`の場合は、`reference_model="infinite_physical"`に対する
`reference_force_error_N`、`reference_work_error_J`、`reference_converged`も判定に使います。
`increment_converged`は、tail proxyと物理referenceの両方を組み合わせたgateです。2回連続で真になったときに、
corrected pathを選択します。`status="not_converged"`の場合、`selected_image_layers`と`selected_path`は`None`です。

cached 無限周期場の opt-in oracle には、(1) 一様非中性 triangle plane の
`E_bottom=0` 解析解 (below `0`、above `sigma/epsilon0`、surface PV
`sigma/(2 epsilon0)`、全 object 力 `Q^2/(2 epsilon0 A)`) と、(2) 中性
`sigma_0 cos(kx)` sheet の `k != 0` 解 (`exp(-k |z-z0|)` 減衰) および triangle mesh
refinement を使います。これらは cache 生成を伴う opt-in 検証で、通常の軽量 test ではありません。

`ObjectForcePath.status="converged"`と`DetachmentResult.numerically_qualified=True`は、指定した経路積分toleranceを満たしたことを示します。
離脱について物理的な結論を得るには、mesh、quadrature、FMM/cache、shell、経路上端、charge snapshot、seedへの依存性も確認します。

## 11. 可視化関数

### 11.1 電荷/電位メッシュ描画

```python
fig, ax = b.plot_mesh(cmap="coolwarm")                        # 電荷密度
fig, ax = b.plot_potential(reference_point="species1_injection_center")  # 電位
```

### 11.2 電位断面

```python
fig, axes = b.plot_potential_slices(
    box_min=[0, 0, 0], box_max=[0.01, 0.01, 0.01],
    xy_z=0.005,
)
```

### 11.3 履歴アニメーション

```python
gif_path = b.animate_mesh("charge_animation.gif", quantity="charge", fps=10)
```

### 11.4 Coulomb 力行列

```python
fig, ax = b.plot_coulomb_force_matrix(component="z")
```

## 12. CLI コマンド

### 12.1 統一 CLI (`beachx`)

v1.0.0 以降は `beachx` 統一 CLI を推奨します。

| コマンド | 説明 |
|---|---|
| `beachx inspect <output_dir>` | 出力ディレクトリの要約表示 |
| `beachx animate <output_dir>` | 電荷/電位履歴のアニメーション GIF 生成 |
| `beachx workload <config.toml>` | ワークロード見積もり |
| `beachx slices <output_dir>` | 電位断面の描画 |
| `beachx profile <output_dir>` | パフォーマンスプロファイルの描画 |
| `beachx coulomb <output_dir>` | Coulomb 力行列の描画 |
| `beachx mobility <output_dir>` | Coulomb mobility 解析 |
| `beachx kernel-forces <output_dir>` | Fortran field kernel による object 別合力 CSV 出力 |
| `beachx object-detachment <output_dir>` | 周期画像を保持した凍結 source の wrench・離脱経路・仕事・速度解析 |
| `beachx lint <config.toml>` | TOML / JSON Schema / BEACH 制約の設定検査 |
| `beachx config validate <config.toml>` | 設定ファイルのバリデーション |
| `beachx model close-pack` | 密充填モデルの生成 |

`beachx inspect`は、`mesh_potential.csv`がある場合に、その事前計算済み配列から
`potential_min` / `potential_max`を表示します。ファイルがなければ電位は再構成せず、この2行を省略します。
破損、非有限値、要素数の不一致がある`mesh_potential.csv`は、入力errorとして扱います。

`--recompute-potential`を指定すると、`mesh_potential.csv`の有無にかかわらず、
`Beach.compute_potential`で電位を再計算して要約に使います。この計算はmesh規模によって$O(N^2)$になり得ます。

`--save-potential-mesh`と`--show`もpotential plotを明示的に要求するため、`--recompute-potential`なしでも電位を計算する場合があります。
`--recompute-potential`とpotential plotを同時に指定すると、要約用とplot用の計算が別々に実行され、電位評価が重複する場合があります。

### 12.2 旧 CLI（非推奨）

以下の旧エントリポイントは後方互換のため残されていますが、将来のバージョンで削除される可能性があります。

| コマンド | 説明 |
|---|---|
| `beach-inspect <output_dir>` | 出力ディレクトリの要約表示 |
| `beach-animate-history <output_dir>` | 電荷/電位履歴のアニメーション GIF 生成 |
| `beach-estimate-workload <config.toml>` | ワークロード見積もり |
| `beach-plot-potential-slices <output_dir>` | 電位断面の描画 |
| `beach-plot-performance-profile <output_dir>` | パフォーマンスプロファイルの描画 |
| `beach-plot-coulomb-force-matrix <output_dir>` | Coulomb 力行列の描画 |

`beach-inspect` にも `beachx inspect` と同じ precomputed-only の通常要約と
`--recompute-potential` の明示的な再計算規則が適用されます。

## 13. 物理定数

| シンボル | 値 | 単位 | 説明 |
|---|---|---|---|
| `K_COULOMB` | `8.9875517923e9` | N m^2 / C^2 | クーロン定数 |

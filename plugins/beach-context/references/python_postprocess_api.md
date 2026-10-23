title: Python 後処理 API リファレンス

# Python 後処理 API リファレンス

BEACHのPython package (`beach`) は、Fortranシミュレーションの結果を読み込み、解析・可視化するための後処理レイヤです。
`summary.txt`、`charges.csv`、`mesh_triangles.csv`などを読み込み、電位の再構成、Coulomb力と電場の計算、
電気力線の追跡、3D可視化を行います。

最初の図を作る手順は [後処理チュートリアル](../guide/PostprocessTutorial.md) にまとめています。

## 1. パッケージ構成

公開 API は `beach` トップレベルの名前です。サブパッケージは依存の層ごとに分かれ、下の層は上の層を import しません。

| サブパッケージ | 依存 | 役割 |
|---|---|---|
| `beach.io` | numpy | 出力ディレクトリの読み込み (`results`: `load_fortran_result`, `list_fortran_runs`)、`summary.txt` の parser (`summary`)、`charge_history.csv` のバッチステップ別アクセス (`history`: `FortranChargeHistory`)、公開データ型 (`types`: `FortranRunResult`, `CoulombInteraction` 等)、物理定数 (`constants`: `K_COULOMB`) |
| `beach.native` | numpy, ctypes | Fortran native field kernel の共有ライブラリ呼び出し (`kernel`: `FieldKernel`)、periodic2 plan の設定型 (`periodic2`)、physical periodic zero mode (`periodic_zero_mode`: `PeriodicZeroMode`) |
| `beach.analysis` | io, native | run からの kernel 構築 (`field_kernel`: `calc_object_forces_kernel`, `field_kernel_options_from_result`)、電位再構成 (`potential`)、電場計算・電気力線追跡 (`field_lines`)、Coulomb 力/トルク (`coulomb`)、Coulomb mobility (`mobility`)、凍結 source に対する object の力・離脱経路 (`object_interaction`, `detachment`)、有限周期画像 shell の収束診断 (`periodic_force_oracle`)、object の一時移動・回転 (`scene`) |
| `beach.plot` | analysis, matplotlib | 電荷/電位の 3D mesh (`mesh`)、電位断面 (`slices`)、Coulomb 力行列 (`matrix`)、mesh source 別箱ひげ図 (`boxplot`)、電気力線 (`field_lines`)、履歴アニメーション (`animation`) |
| `beach.facade` | 全層 | 高水準ファサード `Beach` クラス |
| `beach.config` | numpy 不要 | `beach.toml` の読み込み・正規化・schema 検査 |
| `beach.tooling` | numpy 不要 | ワークロード見積もり (`workload`)、close-pack モデル生成 (`closepack_config`)、runops adapter (`runops_adapter`) |

`beach.plot` は `beach.native` を直接 import せず、kernel が必要な描画は `beach.analysis` の関数を経由します。
`beach.config` と `beach.tooling` は numpy なしで import できます（runops は adapter を numpy の無い環境でも読み込みます）。
この依存方向は `tests/python/test_import_layers.py` が確認します。

すべての公開シンボルは `beach` トップレベルからインポートできます。

```python
from beach import Beach, calc_coulomb, compute_electric_field_points, trace_field_lines
```

## 2. `Beach` セッション

`Beach` は出力ディレクトリ 1 つ分のセッションです。読み込んだ結果、`beach.toml` を解決する `RunContext`、
その run から構築した native field kernel をまとめて持ち、解析・可視化メソッドを提供します。

```python
from beach import Beach

run = Beach("outputs/latest")
```

`Beach("outputs/latest", config_path="path/to/beach.toml")`のように、設定ファイルを明示できます。
`config_path=None`の場合は、`output_dir/beach.toml`、親ディレクトリ、祖父ディレクトリの順に自動探索します。
object 設定を参照する解析では、この `beach.toml` から object kind/order などを取得します。
新しい出力では、全電場再構成に必要な境界・box・一様場・periodic2・実際のsolver・tree/FMMパラメータを
`summary.txt` の `field_reconstruction_*` receipt から取得します。この receipt がある場合、近傍または
`config_path` で指定した別の `beach.toml` は field 設定を上書きしません。receipt のない旧出力だけは、
境界や外部場を推測しないため、実行時と同じ `beach.toml` が必要です。
simulatorが保存した`mesh_potential.csv`をそのまま読む経路では再計算しないため不要です。

### 2.1 コンストラクタ

| パラメータ | 型 | デフォルト | 説明 |
|---|---|---|---|
| `output_dir` | `str \| Path` | `"outputs/latest"` | Fortran 出力ディレクトリ |
| `config_path` | `str \| Path \| None` | `None` | `beach.toml` の明示パス。`None` で近傍を自動探索 |
| `library_path` | `str \| Path \| None` | `None` | native field kernel 共有 library。`None` で `BEACH_FIELD_KERNEL_LIB` または既定の探索 |

### 2.2 プロパティ

| 名前 | 戻り型 | 説明 |
|---|---|---|
| `result` | `FortranRunResult` | 読み込み済み結果（遅延ロード）|
| `context` | `RunContext` | 結果と、解決した `beach.toml`（`context.config_path`）。config 解決の唯一の持ち主 |
| `mesh_ids` | `tuple[int, ...]` | 利用可能な mesh ID 一覧 |
| `library_path` | `Path \| None` | セッションの native kernel library |

`Beach` を第 1 引数に受けたモジュール関数（`compute_potential_mesh(run)` など）も `run.context` を使うため、
`Beach(..., config_path=...)` の設定がそのまま引き継がれます。

### 2.3 メソッド一覧

「重さ」は native field kernel を構築するか（重い）、保存済みの値だけを使うか（軽い）です。`kernel=` 引数を持つ
重いメソッドはセッションの kernel cache（2.4 節）を使います。`analyze_coulomb_mobility` と
`object_interaction_snapshot` は従来どおり専用の kernel を構築します。

| メソッド | 委譲先 | 重さ | 概要 |
|---|---|---|---|
| `kernel(step, periodic2, theta, leaf_max, order)` | `field_kernel_from_result` | 重い（初回のみ） | セッションが所有する cached `FieldKernel` |
| `close()` | 内部 | - | cache 中の kernel をすべて閉じる。`with Beach(...) as run:` でも閉じる |
| `reload()` | `load_fortran_result` | 軽い | kernel を閉じてディスクから再読み込み |
| `get_mesh(*mesh_ids, step)` | 内部 | 軽い | mesh ID で `MeshSelection` を取得 |
| `get_mesh_charge(*mesh_ids, step)` | 内部 | 軽い | mesh ID で要素電荷配列を取得 |
| `calc_coulomb(target, source, ...)` | `calc_coulomb` | 重い | Coulomb 力/トルク計算 |
| `calc_object_forces_kernel(...)` | `calc_object_forces_kernel` | 重い | Fortran field kernel による object 別合力計算 |
| `object_interaction_snapshot(...)` | `ObjectInteractionSnapshot.from_result` | 重い（専用 kernel） | 凍結 source に対する primary-only self exclusion の力・離脱経路解析 |
| `scene(step, ...)` | `BeachScene.from_result` | 軽い | object を一時的に移動・回転する what-if scene |
| `analyze_coulomb_mobility(...)` | `analyze_coulomb_mobility` | 重い | オブジェクト別 mobility 解析 |
| `compute_potential(...)` | `compute_potential_mesh` | 保存値なら軽い | 重心での電位 |
| `compute_potential_points(points, ...)` | `compute_potential_points` | 重い | 任意点での電位 |
| `compute_potential_slices(...)` | `compute_potential_slices` | 重い | XY/YZ/XZ 断面の電位 |
| `compute_electric_field(points, ...)` | `compute_electric_field_points` | 重い | 任意点での電場ベクトル |
| `trace_field_lines(seed_points, ...)` | `trace_field_lines` | 重い | 電気力線の RK4 追跡 |
| `plot_mesh(...)` | `plot_charge_mesh` | 軽い | 電荷密度の 3D メッシュ描画 |
| `plot_potential(...)` | `plot_potential_mesh` | 保存値なら軽い | 電位の 3D メッシュ描画 |
| `plot_potential_slices(...)` | `plot_potential_slices` | 重い | 電位断面の描画 |
| `plot_field_lines(seed_points, ...)` | `plot_field_lines_3d` | 重い | 電気力線の 3D 描画 |
| `plot_bar()` | `plot_charges` | 軽い | 要素電荷棒グラフ |
| `plot_mesh_source_boxplot(...)` | `plot_mesh_source_boxplot` | 保存値の無い電位なら重い | mesh source 別箱ひげ図 |
| `plot_coulomb_force_matrix(...)` | `plot_coulomb_force_matrix` | 重い | Coulomb 力行列プロット |
| `animate_mesh(...)` | `animate_history_mesh` | `quantity="potential"` なら重い | 電荷/電位履歴アニメーション |

委譲メソッドは委譲先と同じ引数・既定値・docstring を持ちます（`help(Beach.plot_potential)` で確認できます）。
既定値は委譲先の関数だけが持ちます。メソッドは `result` に自分自身を渡し、`library_path` を省略すると
セッションの `library_path` を使います。

### 2.4 kernel の使い回し

`run.kernel(...)` は、run の電荷を載せた native kernel を引数ごとに 1 回だけ構築して cache します。鍵は電荷
（`step`）、`periodic2`、`theta`、`leaf_max`、`order`、セッションの `library_path` です。`step=-1`（最新の履歴
snapshot）と `step=None`（`charges.csv`）は、電荷が一致すれば同じ kernel を共有します。kernel はセッションが
所有するので、呼び出し側で閉じません。`close()`、`reload()`、`with` ブロックの終了で閉じます。

重いメソッドは、同じ plan の cached kernel を `kernel=` として委譲先に渡します。電位系（`compute_potential*`、
`plot_potential*`）は `charges.csv` の電荷（`step=None`）、電場・力・アニメーション系は `step=-1` の kernel を
使います。`kernel=`、`config_path=`、セッションと異なる `library_path=` を明示した呼び出しは cache を使わず、
その呼び出しで kernel を構築して閉じます。保存済み `mesh_potential.csv` を使う `compute_potential()` /
`plot_potential()`、`quantity="charge"` のアニメーション、`workers > 1` の力行列は kernel を構築しません。

```python
import numpy as np
from beach import Beach

seeds = np.array([[0.005, 0.005, 0.008]])
with Beach("outputs/latest", library_path="build/libbeach_field_kernel.so") as run:
    fig, _ = run.plot_potential(source="kernel")    # kernel を構築する
    fig, _ = run.plot_potential_slices(box_min=[0, 0, 0], box_max=[0.01, 0.01, 0.01])  # 同じ kernel
    lines = run.trace_field_lines(seeds)            # 最新履歴が charges.csv と同じなら同じ kernel
    records = run.calc_object_forces_kernel()       # 電荷を一時的に差し替えて戻す
# ここで kernel を閉じる
```

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
オプションファイル: `mesh_triangles.csv`, `mesh_sources.csv`, `charge_history.csv`, `mesh_potential.csv`,
`matching_plane_history.csv`, `secondary_emission_history.csv`

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
| `matching_plane_state` | `MatchingPlaneState \| None` | `zhao_stationary`の外部根を更新したrunの最後の外部状態。summaryでstateが無効なら`None` |
| `matching_plane_history` | `tuple[MatchingPlaneHistoryEntry, ...] \| None` | 外部状態の履歴。stateが無効またはCSVがなければ`None` |
| `secondary_emission_history` | `tuple[SecondaryEmissionHistoryEntry, ...] \| None` | SEE の入射種・受理バッチごとの履歴。CSV なしなら `None` |

`MatchingPlaneState`はセル平均電荷$D_H$、壁電位$\Phi_H=\phi_0$、流入束、外部の放出源、PE return / escape、
採用した外部根の数、直前の更新での相対変化を名前付きfieldで保持します。各fieldの意味は
[出力形式リファレンス](OutputReference.md#zhao_stationary)を参照してください。`MatchingPlaneHistoryEntry`は`batch`、`simulated_time_s`と、その時点の
`state: MatchingPlaneState`を保持します。`matching_plane_history.csv`がheaderだけなら履歴は空tuple `()`です。
このため、`None`は「履歴ファイルがない、または有効な外部状態ではない」、空tupleは
「有効なstateだが記録行がない」と区別できます。

`MatchingPlaneState.secondary_share` はそのバッチで SEE の脱出に充てた外部根の脱出束の割合、`secondary_budget_residual_a_m2` は差し引けなかった残差です。
`SecondaryEmissionHistoryEntry` は `batch, incident_species, hits, emitted_count, emitted_charge_c, mean_yield, generation_capped_count` を持ちます。

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

`compute_potential_mesh` は保存済み `mesh_potential.csv` を使える場合だけ軽く、それ以外の電位関数は native
field kernel を構築します（重い）。どの関数も `kernel=` で構築済みの kernel を受け取れます（10.1 節）。

`reference_point` の既定は、計算（`compute_potential_mesh` / `compute_potential_points` / `compute_potential_slices`）と
描画（`plot_potential_mesh` / `plot_potential_slices` / `animate_history_mesh`、`Beach` の同名メソッド）で共通の
`None`（絶対電位）です。species 1 の注入面中心を基準にした電位差は `reference_point="species1_injection_center"`
で明示します。

### 4.1 `compute_potential_mesh(result, *, periodic2, reference_point, config_path, library_path, source, kernel)`

三角形重心の電位を返します。保存済みの `mesh_potential.csv` を返すか、native field kernel の P0 triangle panel
評価で計算するかを `source` で選びます。

| `source` | 評価 |
|---|---|
| `"auto"`（既定） | `mesh_potential.csv` があり、`periodic2`・`reference_point`・`kernel` のどれも指定していなければ保存値（軽い）。それ以外は kernel（重い） |
| `"saved"` | 保存値だけを返す。`mesh_potential.csv` が無い場合と、`periodic2`・`reference_point`・`kernel` との併用は `ValueError` |
| `"kernel"` | 保存値があっても kernel で再計算する |

| パラメータ | 型 | デフォルト | 単位 | 説明 |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (必須) | - | 結果オブジェクト |
| `periodic2` | `Mapping \| None` | `None` | - | 2 軸周期設定（後述）。`None` で自動判定 |
| `reference_point` | `Iterable[float] \| str \| None` | `None` | m | 基準電位点。`"species1_injection_center"` で species 1 注入面中心 |
| `config_path` | `str \| Path \| None` | `None` | - | `beach.toml`の明示パス |
| `library_path` | `str \| Path \| None` | `None` | - | native field kernel共有libraryの明示パス |
| `source` | `str` | `"auto"` | - | `"auto"` / `"saved"` / `"kernel"`（上表） |
| `kernel` | `FieldKernel \| None` | `None` | - | この run から構築した kernel（10.1 節）。渡すとその plan と電荷で評価し、閉じない。`periodic2` と `library_path` は使わない |

戻り値: `ndarray (mesh_nelem,)` [V]

### 4.2 `compute_potential_points(result, points, *, chunk_size, periodic2, reference_point, config_path, library_path)`

任意 3D 点での電位を native field kernel で計算します（重い）。

| パラメータ | 型 | デフォルト | 単位 | 説明 |
|---|---|---|---|---|
| `points` | `ndarray (n_points, 3)` | (必須) | m | サンプリング点座標 |
| `chunk_size` | `int` | `2048` | - | チャンク分割数 |
| `periodic2` | `Mapping \| None` | `None` | - | `None` で自動判定 |
| `reference_point` | `Iterable[float] \| str \| None` | `None` | m | 基準電位点 |
| `config_path` | `str \| Path \| None` | `None` | - | `beach.toml`の明示パス |
| `library_path` | `str \| Path \| None` | `None` | - | native field kernel共有libraryの明示パス |
| `kernel` | `FieldKernel \| None` | `None` | - | この run から構築した kernel（10.1 節）。渡すとその plan と電荷で評価し、閉じない。`periodic2` と `library_path` は使わない |

戻り値: `ndarray (n_points,)` [V]

### 4.3 `compute_potential_slices(result, *, box_min, box_max, grid_n, xy_z, yz_x, xz_y, ...)`

XY/YZ/XZ 平面上の電位断面を native field kernel で計算します（重い）。断面の全点を 1 つの kernel で評価し、
`kernel=` を含む残りの引数は `compute_potential_points` と同じです。

戻り値: `dict[str, PotentialSlice2D]` (キー: `"xy"`, `"yz"`, `"xz"`)

## 5. Coulomb 力/トルク計算

### 5.1 `calc_coulomb(result, target, source, *, step, torque_origin, periodic2, quadrature_order, config_path, library_path)`

target メッシュグループが source から受ける Coulomb 力とトルクを native field kernel で計算します（重い）。

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
| `kernel` | `FieldKernel \| None` | `None` | - | この run の全 mesh で構築した kernel。呼び出し中だけ source 以外の電荷を 0 にし、終了後に元の電荷へ戻す |

戻り値: `CoulombInteraction`

`kernel=None` では source group の panel だけで kernel を構築します。`kernel=` を渡すと、その kernel の plan
（全 mesh）で source group の電荷だけを評価し、一様外部場の寄与は結果から除きます。両者は FMM の打ち切り誤差の
範囲で一致します。

### periodic2 パラメータによる周期クーロン和

`periodic2` が指定された場合、native field kernelをその2軸周期設定で構成して力とトルクを計算します。

`periodic2=None`（デフォルト）の場合は、出力ディレクトリ近傍の `beach.toml` を探索します。
`fields.boundary="periodic2"` なら、`domain.periodic_axes` と `domain.box_min` / `domain.box_max` から周期設定を
自動構築します。これは `compute_potential_mesh` などの他の関数と共通の自動判定です。

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

任意 3D 点での電場ベクトルを、native field kernelのP0 triangle panel評価で計算します（重い）。

| パラメータ | 型 | デフォルト | 単位 | 説明 |
|---|---|---|---|---|
| `result` | `FortranRunResult \| object` | (必須) | - | 結果オブジェクト |
| `points` | `ndarray (n_points, 3)` | (必須) | m | サンプリング点座標 |
| `chunk_size` | `int` | `2048` | - | チャンクサイズ |
| `periodic2` | `Mapping \| None` | `None` | - | 2 軸周期設定。`None` で自動判定 |
| `config_path` | `str \| Path \| None` | `None` | - | `beach.toml`の明示パス |
| `library_path` | `str \| Path \| None` | `None` | - | native field kernel共有libraryの明示パス |
| `kernel` | `FieldKernel \| None` | `None` | - | この run から構築した kernel（10.1 節）。渡すとその plan と電荷で評価し、閉じない。`periodic2` と `library_path` は使わない |

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

シード点から電場方向（または逆方向）に RK4 積分で電気力線を追跡します（重い: native field kernel）。

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
| `kernel` | `FieldKernel \| None` | `None` | - | この run から構築した kernel（10.1 節）。渡すとその plan と電荷で評価し、閉じない。`periodic2` と `library_path` は使わない |

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

電気力線を 3D 描画し、オプションでメッシュ表面を電荷密度で着色してオーバーレイします（重い: 追跡に native field
kernel を使う）。`kernel=` は `trace_field_lines` と同じです。

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

保存済み出力の `receipt.as_legacy_sim_mapping()` は 1.x 出力契約の内部 reader 設定を返し、この名前は公開 TOML の入力キーではありません。

### 8.1 自動判定（`periodic2=None`）

デフォルトの `None` では、まず出力の `field_reconstruction_*` receipt から、実行時に解決済みの
field boundary、周期軸、長さ、原点を復元します。receipt のない旧出力では、出力ディレクトリ近傍の
`beach.toml` を探索し、`fields.boundary="periodic2"` の場合は `domain.periodic_axes` と
`domain.box_min` / `domain.box_max` から自動解決します。旧出力で設定ファイルが見つからない場合、
全電場を再計算する高水準APIは自由空間へ黙ってfallbackせず停止します。field boundary が `free` に
解決される場合は自由空間として再計算します。
`panel_spectral_reference`のsplit-periodic receiptはfinite-image native kernelへ置換せず、全電場を再計算する
高水準APIで停止します。この場合はsimulator保存済み電位または専用のperiodic reference検証を使ってください。
receiptの実際のsolverが`direct`なら非周期exact-direct kernel、`fmm`ならreceiptのFMM展開次数を使います。
`treecode`はFMMへ置換せず、高水準の電場・電位・力APIで停止します。明示的な`periodic2`を渡しても
`treecode`をFMMへすり替えません。

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

`make build-kernel` で生成した `build/libbeach_field_kernel.so` を `ctypes` 経由で読み込み、シミュレーションと同じFortran field coreでP0 triangle panelの電場・電位を評価します。新しい出力では `field_reconstruction_*` receipt、旧出力では `config_path` または自動探索された `beach.toml` から solver、periodic2、tree 設定を読みます。

```python
from beach import Beach

with Beach("outputs/latest") as run:
    kernel = run.kernel()          # セッションが所有し、with の終了で閉じる
    e = kernel.eval_e([[0.0, 0.0, 0.01]])
    phi = kernel.eval_phi([[0.0, 0.0, 0.01]])
```

`FortranRunResult` から直接作る場合は `beach.analysis.field_kernel.field_kernel_from_result(result, ...)` を使います。
こちらは呼び出し側が所有するので、`with field_kernel_from_result(result) as kernel:` のように閉じます。

#### `kernel=` 引数の規約

電位・電場・力線・Coulomb 力・object 合力の関数と、その描画関数（`plot_potential_mesh`、`plot_potential_slices`、
`plot_field_lines_3d`、`plot_coulomb_force_matrix`、`animate_history_mesh(quantity="potential")`）は
`kernel: FieldKernel | None = None` を受け取ります。

- `None` では、従来どおり呼び出しごとに kernel を構築し、終了時に閉じます。
- 渡した kernel は閉じません。`Beach.kernel()` または `field_kernel_from_result()` で同じ run から構築したものを
  渡します。panel 数が run と異なる kernel と、全電場を表さない plan（未合成の `cached_kneq0` など）は
  `ValueError` です。
- 渡した kernel はその plan のまま評価します。`periodic2`、`theta`、`leaf_max`、`order`、`library_path` は
  kernel の構築に使わないため効きません。描画関数の `periodic2` は、表示 mesh の wrap と複製にだけ使います。
- `step` を持たない関数（電位・電場・力線）は kernel の電荷で評価します。`step` を持つ関数（`calc_coulomb`、
  `calc_object_forces_kernel`）とアニメーションは、呼び出し中だけ kernel の電荷を差し替え、終了時に元の電荷へ
  戻します。
- Coulomb 力行列（`plot_coulomb_force_matrix`）は worker process 間で kernel を共有できないため、
  `workers == 1` のときだけ渡された kernel を使います。`workers > 1` では各 entry が source-only kernel を構築します。

共有ライブラリを別パスに置く場合は `library_path=` または環境変数 `BEACH_FIELD_KERNEL_LIB` を指定します。

`eval_e()` / `eval_phi()` は、構築時に選んだ free / finite periodic / cached periodic 設定を使います。
自動構築ではresolved `direct`をexact-directへ、resolved `fmm`をFMMへ対応させ、resolved `treecode`は
再現不能として停止します。FMMの`order`はreceiptの展開次数がデフォルトで、明示`order=`だけが上書きします。
resolved `direct`では`eval_e()`、`eval_phi()`、`force_on_charges()`のすべてがexact-directを使い、uniform
`fields.external.electric_v_m`はそれぞれ一度だけ合成します。
ただし cached plan の返り値は `cached_kneq0` の $k\ne0$ 成分と `fields.external.electric_v_m` だけで、physical zero mode を含む
simulator 全電場ではありません。`FieldKernel` はこの成分を明示的に扱う低水準 API です。
`eval_e_direct()` / `eval_phi_direct()`は、同じsource geometryとchargeを非周期のexact direct法で評価する診断APIです。
periodic planでは使用できず、uniform `fields.external.electric_v_m`も加算しません。FMMの精度確認とprimary free-space subtractionの
小規模oracleに使います。

### 10.2 `calc_object_forces_kernel(result, ...)`

各objectについて、自身のsource電荷を0にしたうえで、`sum(q_i E_not_self(r_i))`とトルクを計算します（重い）。
`kernel=` を渡すと、object ごとの電荷の差し替えをその kernel で行い、終了時に元の電荷へ戻します。
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

`FieldKernelOptions.periodic_cache_dir=None` は、初期化時に `BEACH_CACHE_DIR`、`XDG_CACHE_HOME`、`HOME` の順で
共通の保存先を解決します。明示した文字列と結果ファイルの再構成情報は優先されます。

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
| `"configured"` | run の `beach.toml` をそのまま使用する。free、`fields.periodic.backend="finite_images"` の有限画像、または cached periodic のいずれにもなり得る |
| `"infinite_physical"` | x/y periodic run に対して cached `k != 0` と、`fields.periodic.lower_boundary_model`に従う物理的な`k = 0` modeを組み合わせる。cache の生成・再利用条件を満たす必要がある |

完全な `beach.toml` と `domain.box_min` / `domain.box_max` が必要です。

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
fig, ax = b.plot_potential()                                   # 絶対電位（mesh_potential.csv があれば保存値）
fig, ax = b.plot_potential(source="kernel")                    # 保存値があっても kernel で再計算
fig, ax = b.plot_potential(reference_point="species1_injection_center")  # species 1 注入面中心との電位差
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
| `beachx model close-pack` | 密充填モデルの生成 |

`beachx inspect`は、`mesh_potential.csv`がある場合に、その事前計算済み配列から
`potential_min` / `potential_max`を表示します。ファイルがなければ電位は再構成せず、この2行を省略します。
破損、非有限値、要素数の不一致がある`mesh_potential.csv`は、入力errorとして扱います。

`--recompute-potential`を指定すると、`mesh_potential.csv`の有無にかかわらず、
native field kernelで絶対電位を再計算して要約に使います。この計算はmesh規模によって$O(N^2)$になり得ます。

`--save-potential-mesh`と`--show`のpotential plotは絶対電位を描きます。`mesh_potential.csv`があればその値を使い、
なければ`--recompute-potential`なしでも電位を計算します。`mesh_potential.csv`がない出力で`--recompute-potential`と
potential plotを同時に指定すると、要約用とplot用の計算が別々に実行され、電位評価が重複します。

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

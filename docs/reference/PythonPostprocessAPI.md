title: Python 後処理 API

# Python 後処理 API

BEACH の Python API で、保存結果の読み込みから電位・電場・物体力の解析へ進めます。
通常は `Beach` をセッションとして使い、読み込みだけなら `load_run` を使います。
全クラス・関数の引数、戻り値、単位は、実装から生成する
[Python API リファレンス](https://nkzono99.github.io/BEACH/python/)で確認できます。
最初の図を作る手順は[後処理チュートリアル](../guide/PostprocessTutorial.md)にあります。

## import と依存の層

通常の後処理の公開名は `beach` から import します。設定と作業支援はそれぞれのサブパッケージを使います。

| import 先 | 用途 | 主な依存 |
|---|---|---|
| `beach` | `Beach`、`load_run`、`list_runs`、`RunResult`、`BatchHistory`、解析・描画関数、`K_COULOMB` | 名前を参照すると必要な層だけ読み込む |
| `beach.io` | `results` の reader、`history`、`types`、`summary` | numpy。`summary` は numpy 不要 |
| `beach.native` | `FieldKernel` と `PeriodicZeroMode` の共有ライブラリ呼び出し | numpy、ctypes |
| `beach.analysis` | 電位、電場、物体間力、離脱、剛体変換 | io、native |
| `beach.plot` | mesh、断面、力行列、履歴の描画 | analysis、matplotlib |
| `beach.config` | 設定の読み込み、正規化、schema 検査 | numpy 不要 |
| `beach.tooling` | 負荷見積もり、密充填モデル生成、runops adapter | numpy 不要 |

## `RunResult` を読み込む

```python
from beach import load_run, list_runs

result = load_run("outputs/latest")
print(result.summary.batches, result.summary.processed_particles)
print(result.charges.sum())                       # 最終電荷 [C]
if result.mesh.triangles is not None:
    print(result.mesh.triangles.shape)             # (nelem, 3, 3)
print(list_runs("outputs"))
```

`RunResult` は出力を用途別に保持します。属性の型と全項目は
[`RunResult` の生成リファレンス](https://nkzono99.github.io/BEACH/python/public.html#beach.RunResult)にあります。

| 属性 | 内容 |
|---|---|
| `summary` | 粒子数、バッチ数、相対変化、計算時間などの基本統計 |
| `mesh` | `triangles`、`mesh_ids`、`mesh_sources`、保存済みの `mesh_potential_v` |
| `charges` | `charges.csv` の最終要素電荷 [C] |
| `history` | `charge` と `potential` の `BatchHistory` |
| `receipts` | 場の再構成情報、外部シースの最終状態、fingerprint、checkpoint schema |
| `diagnostics` | 周期 cache、境界イベント、適応バッチの診断 |
| `field` | 場の source model と kernel 識別子 |
| `ledger` | 電荷収支の履歴 |
| `matching_plane_history`、`secondary_emission_history` | 外部シースと二次電子の履歴 |
| `top_reference_history`、`fixed_current_history`、`matching_plane_spectrum_history` | 上端基準、固定電流、上端 Fourier 成分の履歴 |
| `performance_profile` | 計測区間ごとの処理時間 |

任意のファイルがないとき、その配列・履歴は `None` です。行形式の履歴は、ファイルがあって記録行がない場合は
空 tuple `()` です。無効な外部シース状態など、各 reader 固有の扱いは生成リファレンスを参照してください。
CSV 列の定義と単位は[出力形式リファレンス](OutputReference.md)が正本です。

### `BatchHistory` で保存バッチを選ぶ

```python
history = result.history.charge
if history is not None and history.has_data:
    print(history.batch_indices)
    latest = history.get_step(-1)                  # 最後に保存されたバッチ
    first = history.get_step(int(history.batch_indices[0]))

potential_history = result.history.potential
if potential_history is not None and potential_history.has_data:
    latest_potential_v = potential_history.get_step(-1)
```

整数の `step` は CSV のバッチ番号、`-1` は最後の保存バッチです。解析関数の `step=None` は最終電荷を選びます。
履歴は最初の参照時に CSV 全体を検証して索引を作り、各バッチの値は要求時に読み込みます。
`as_array()` は全履歴を配列化するため、その分のメモリを使います。
要素の欠損や重複、不正な順序、非有限値は `ValueError` になり、電荷 `0 C` では補完しません。

## `Beach` セッションで解析する

```python
from beach import Beach

with Beach("outputs/latest", config_path="beach.toml") as run:
    print(run.result.summary.absorbed)
    fig, ax = run.plot_mesh()
    fig.savefig("charge.png", dpi=150)
    potential_fig, _ = run.plot_potential()
    potential_fig.savefig("potential.png", dpi=150)
```

`result` は最初の参照時に読み込みます。`context` が設定の解決を保持し、`kernel()` は同じ電荷・周期条件・
FMM 設定の native kernel を再利用します。転送メソッドを別の `step` で呼ぶときは、同じ周期条件・FMM 設定の
kernel が開いていればその step の電荷を一時的に読み込んで使い、呼び出し後に元の電荷へ戻します。
多数の step を解析しても kernel は設定ごとに 1 つです。`close()` と `with` の終了で kernel を閉じ、
`reload()` は kernel を閉じて出力を読み直します。セッションから借りた kernel を自分で閉じないでください。

関数を直接呼ぶ場合も同じ結果と引数を使います。例えば `compute_potential_mesh(run)` と
`run.compute_potential()` は同じ解析を行います。`Beach` の転送メソッドの引数と既定値は関数の定義から引き継ぎます。

### native kernel を用意する

保存済みの電荷・電位を読むだけなら共有ライブラリは不要です。任意点の電位・電場、力、電気力線を計算する場合は、
Fortran compiler と `fpm` が使える source checkout で次を実行します。

```bash
make -C /path/to/BEACH build-kernel
export BEACH_FIELD_KERNEL_LIB=/path/to/BEACH/build/libbeach_field_kernel.so
```

`library_path=` で場所を指定することもできます。`Beach.kernel()` の作成後は同じセッションで複数の解析に使えます。
`beach.analysis.field_kernel.field_kernel_from_result` で直接作る場合は、呼び出し側が `with` または `close()` で閉じます。

`FieldKernelOptions` に `box_min` と `box_max` を指定する場合は、両方とも有限な 3 成分の座標にします。
周期計算でも、長さが違う配列や片方だけの指定は、共有ライブラリへ渡す前に `ValueError` になります。

`kernel=` を受け付ける関数に渡した kernel は、その plan と電荷を使い、閉じません。plan を決める
`periodic2`、`theta`、`leaf_max`、`order`、`library_path` は kernel の作成時に指定してください。
`step` を使う力の計算やアニメーションは一時的に電荷を入れ替え、終了時に戻します。
並列の力行列では process 間で kernel を共有できず、`workers=1` の場合だけ明示 kernel を使います。

## 電位と電場

すべての電位の計算・描画で、`reference_point=None`（絶対電位）が既定です。
座標 [m] または `"species1_injection_center"` を渡すと、その点との差を求めます。

| `source` | mesh 電位の取得方法 |
|---|---|
| `"auto"` | `mesh_potential.csv` があり、`periodic2`・`reference_point`・`kernel` を指定しなければ保存値。それ以外は再計算 |
| `"saved"` | 保存値のみ。ファイルがない場合や上記の再計算用引数との併用は `ValueError` |
| `"kernel"` | native kernel で再計算 |

```python
import numpy as np
from beach import Beach

points = np.array([[0.0, 0.0, 0.01]])               # 評価点 [m]
with Beach("outputs/latest", config_path="beach.toml") as run:
    phi_v = run.compute_potential_points(points)
    e_v_m = run.compute_electric_field(points)
    lines_m = run.trace_field_lines(points, max_steps=100)
```

電気力線は場方向の RK4 積分であり、粒子軌道ではありません。mesh との衝突では停止せず、場の閾値、
最大ステップ数、box 外への退出で打ち切ります。引数と描画関数は
[生成 API](https://nkzono99.github.io/BEACH/python/public.html#beach.trace_field_lines)にあります。

### `periodic2` と再構成できる場

`periodic2=None` は保存された `field_reconstruction_*` 情報を使い、情報のない出力では近傍の `beach.toml` を探します。
有限画像を明示する例は `{"axes": [0, 1], "lengths": [0.01, 0.01], "image_layers": 2}` です。
軸は 0-based、長さと原点は m、画像数は各軸で `[-N, N]` です。全設定型は
[`Periodic2Config`](https://nkzono99.github.io/BEACH/python/advanced.html#beach.native.periodic2.Periodic2Config)にあります。

実際の solver が `direct` なら非周期の直接評価、`fmm` なら記録された展開次数を使います。
`treecode`、`panel_spectral_reference`、未合成の `cached_kneq0` を別の全電場へ置き換えることはありません。
一般の電位・電場・力 API で再現できない場合は保存値を使い、cached 周期場の物体力には次節の API を使います。
場境界の定義は[場 solver を選ぶ](../guide/FieldSolvers.md)、周期成分は[周期静電場](../models/PeriodicElectrostatics.md)を参照してください。

## 物体の力と離脱

`calc_coulomb` は選択した source から target への力・トルク、`calc_object_forces_kernel` は物体別合力を返します。
後者は target 自身の周期画像も除外します。target の primary 自己場だけを除き、自身の周期画像を残す場合は
`ObjectInteractionSnapshot.object_probe` を使います。`periodic_model="infinite_physical"` は x/y 周期の
`cached_kneq0` と物理的な平均場を合成します。

凍結した source に対して target を動かす経路、付着、仕事、離脱判定、剛体変換の手順は
[物体の力と離脱を調べる](../guide/ObjectForcesDetachment.md)にあります。
`ObjectWrench.components` は他物体の全画像 (`other_objects_all_images`)、自身の周期画像 (`target_periodic_images`)、
一様外部場 (`external_uniform`)、合計 (`total_external`) を保持します。トルクと一緒に `torque_origin_m` と
`numerical_metadata["torque_origin_policy"]` を保存してください。

表面の機械的な合力は平均極限値（principal value）を使います。診断値
`cached_kneq0_trace_correction` はすでに `periodic_kneq0` に含まれており、合力に再加算しません。
鉛直経路では source の形状と電荷を固定し、target だけを移動します。外部場に対するエネルギーは
`sum(q_i * phi_env(r_i))` で、係数 `1/2` は付けません。

`Beach.scene()` は保存電荷を保ったまま物体を移動・回転し、新しい配置で力を比較するための機能です。
`move` と `rotate` は新しい scene を返し、元の保存結果を書き換えません。

有限画像数を増やす `finite_shell_convergence` は、遠方の残りの寄与を見積もる `force_tail_proxy_N` と
`work_tail_proxy_J` で判定します。`infinite_physical` の比較では基準解との差も使い、2 回連続の合格を要求します。
隣接する画像数での差だけでは判定しません。`status="not_converged"` の場合は
`selected_image_layers` と `selected_path` が `None` です。

`ObjectForcePath.status` と `DetachmentResult.numerically_qualified` は数値的な判定であり、
物理的な離脱の結論には mesh、積分次数、場、経路上端、電荷の選び方への依存性も確認します。

## CLI から使う

`beachx --help` と `beachx <command> --help` がコマンドと全オプションの一覧です。
読み込み・図は[後処理チュートリアル](../guide/PostprocessTutorial.md)、物体解析は
[力と離脱の手順](../guide/ObjectForcesDetachment.md)を参照してください。

| オプション | 意味 |
|---|---|
| `--save PATH` | animate、coulomb、mobility、kernel-forces、slices、profile の 1 ファイルの保存先 |
| `--save-dir DIR` | inspect と object-detachment の複数ファイルの保存先 |
| `--json` | inspect、mobility、kernel-forces、coulomb、object-detachment、workload の標準出力を JSON 1 オブジェクトにする |
| `--config PATH`、`--library PATH` | 設定と native kernel の場所 |
| `--step final\|BATCH` | 最終電荷または保存バッチ。`-1` は最新の保存バッチ |

標準出力は `key=value` 行で、表は 1 行 1 レコードです。ファイルの保存先は `saved=<path>`、複数なら
`saved_<name>=<path>` です。機械処理には `--json` を使います。非有限値は文字列 `"inf"`、`"-inf"`、`"nan"` です。
`inspect` は保存済み mesh 電位がある場合だけ電位を要約し、`--recompute-potential` で再計算します。
`coulomb` の力の報告は選択した target ごとに並びます。source 別の値は source の選択順で、
合力はその全 source の和です。API が返す行列は source 行 × target 列です。

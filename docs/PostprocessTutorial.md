title: 後処理チュートリアル

Lang: [日本語](PostprocessTutorial.md) | [English](PostprocessTutorial.en.md)

# 後処理チュートリアル

BEACH の後処理は、実行概要を CLI で確認し、電荷履歴を動画にし、必要なときだけ Python で独自解析へ進むのが
基本です。通常の BEACH Python package だけで、作成済みの公式入門ケースの完走確認、電荷分布、保存済みの最終絶対電位、
保存 batch の選択まで実行できます。電位履歴の再構成と電位断面は、native field kernel を用意した追加経路です。

**前提:** [10 分チュートリアル](Tutorial.html)を完了し、現在のディレクトリに、その実行で使った
`beach.toml` と `outputs/tutorial` があります。repository の clone、repository root、
`examples/tutorial_insulator.toml` は必要ありません。

`outputs/tutorial` には少なくとも `summary.txt`、`charges.csv`、`mesh_triangles.csv`、
`mesh_sources.csv`、`charge_history.csv`、`mesh_potential.csv`、`potential_history.csv` が必要です。
不足する場合は[10 分チュートリアル](Tutorial.html)へ戻って実行し直してください。file の意味は
[出力ファイルを調べる](OutputGuide.html)、全 class・関数は
[Python 後処理 API リファレンス](PythonPostprocessAPI.html)を参照してください。

## 1. `beachx inspect` で完走と分布を確認する

最初に、図を作らず実行概要を読みます。

```bash
beachx inspect outputs/tutorial
```

公式入門ケースでは、少なくとも次を確認します。

- `directory=outputs/tutorial`
- `processed_particles=4000` と `batches=20`
- `absorbed`、`escaped`、`charge_sum` が表示される
- 保存済み `mesh_potential.csv` の `potential_min` と `potential_max` が表示される
- `charge_history_shape` と保存された batch index が表示される
- `mesh_ids` と `mesh_source` が表示される

`rng_seed=12345`、1 OpenMP thread での参照値は
[出力ファイルを調べる](OutputGuide.html#公式入門ケースの結果を確認する)にあります。thread 数が異なると、
乱数列に依存する吸収・escape 件数が変わる場合があります。

次に、同じ出力から基本図を保存します。

```bash
beachx inspect outputs/tutorial \
  --save-bar outputs/tutorial/charges_bar.png \
  --save-mesh outputs/tutorial/charge.png
```

正常終了すると、要素電荷の bar chart と、表面電荷密度で着色した mesh が指定先に保存されます。
画像の生成は後処理の成功を示しますが、数値収束や物理的妥当性の証明ではありません。

## 2. `beachx animate` で batch 間の変化を見る

公式入門ケースは `history_stride=1` なので、20 個の電荷 snapshot があります。標準配布だけで使える
表面電荷密度の GIF を作ります。

```bash
beachx animate outputs/tutorial \
  --quantity charge \
  --save-gif outputs/tutorial/charge_history.gif \
  --total-frames 20
```

正常終了すると、`saved_gif=outputs/tutorial/charge_history.gif`、`snapshots=20`、
`rendered_frames=20` が表示されます。

`animate` は `charge_history.csv` の保存済み snapshot を使います。`potential` mode では、その電荷と
場の再構成情報に加えて native field kernel が必要です。履歴がない別の run は、`output.history_stride` を
正にして simulation を再実行する必要があります。

## 3. Python API で図と保存 batch を選ぶ

CLI の確認が通ったら、独自の図や集計へ進みます。`config_path` には、現在のディレクトリで実行に使った
`beach.toml` を指定します。

```python
from beach import Beach

run = Beach(
    "outputs/tutorial",
    config_path="beach.toml",
)

print("absorbed:", run.result.absorbed)
print("escaped:", run.result.escaped)
print("batches:", run.result.batches)
print("mesh ids:", run.mesh_ids)

bar_fig, _ = run.plot_bar()
bar_fig.savefig("outputs/tutorial/python_charges_bar.png", dpi=150)

charge_fig, _ = run.plot_mesh()
charge_fig.savefig("outputs/tutorial/python_charge.png", dpi=150)

potential_fig, _ = run.plot_potential(reference_point=None)
potential_fig.savefig("outputs/tutorial/python_potential.png", dpi=150)
```

`reference_point=None` は `mesh_potential.csv` に保存された最終絶対電位を明示的に選ぶため、このケースでは
native field kernel は不要です。保存値と異なる基準で再評価する場合は、次節の kernel が必要です。

### 特定の mesh と batch を選ぶ

`beachx inspect` または `mesh_sources.csv` で `mesh_id` を確認してから選択します。`step=None` は
`charges.csv` の最終電荷、整数の `step` は `charge_history.csv` に保存された batch index を表します。

```python
mesh_id = run.mesh_ids[0]

final_mesh = run.get_mesh(mesh_id, step=None)
final_charge = run.get_mesh_charge(mesh_id, step=None)
print(final_mesh.triangles.shape, final_charge.shape)

# 公式ケースは history_stride=1 なので batch 10 が保存されています。
mesh_at_10 = run.get_mesh(mesh_id, step=10)
charge_at_10 = run.get_mesh_charge(mesh_id, step=10)
print(mesh_at_10.triangles.shape, charge_at_10.shape)
```

## 4. native field kernel で電位履歴と断面を再構成する

通常の Python package install には `libbeach_field_kernel.so` が含まれません。電位履歴を再構成する
`animate --quantity potential` と、任意点を評価する `slices` を使う場合だけ、Fortran compiler と `fpm` が
使える BEACH source checkout で library を build します。`/path/to/BEACH` は実際の checkout へ置き換えます。

```bash
make -C /path/to/BEACH build-kernel
export BEACH_FIELD_KERNEL_LIB=/path/to/BEACH/build/libbeach_field_kernel.so
```

環境変数を設定した shell で、まず電位履歴を作ります。

```bash
beachx animate outputs/tutorial \
  --quantity potential \
  --save-gif outputs/tutorial/potential_history.gif \
  --total-frames 20
```

正常終了すると、`saved_gif=outputs/tutorial/potential_history.gif`、`snapshots=20`、
`rendered_frames=20` が表示されます。次に、simulation box 中央の XY / YZ / XZ 断面を作ります。

```bash
beachx slices outputs/tutorial \
  --config beach.toml \
  --grid-n 200 \
  --save outputs/tutorial/potential_slices.png
```

正常終了すると `saved_potential_slices=outputs/tutorial/potential_slices.png` が表示されます。これらは保存電荷から
場を再構成する解析であり、mesh や solver に対する収束確認の代わりにはなりません。研究結果として採用する前に
[計算結果の妥当性確認](ValidationGuide.html)に従ってください。

## 5. 複数物体の解析へ進む

公式入門ケースは単一の固定平面なので、物体間の力や離脱には使えません。別の複数物体 run がある場合は、
[物体の力と離脱を調べる](ObjectForcesDetachment.html)へ進んでください。

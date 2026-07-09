title: 後処理チュートリアル

Lang: [日本語](PostprocessTutorial.md) | [English](PostprocessTutorial.en.md)

# 後処理チュートリアル

このページは、最初の実行結果を CLI と Python で確認するための短い手順です。
出力ファイルの意味は [出力の読み方](OutputGuide.html)、API の詳細は [Python 後処理 API リファレンス](PythonPostprocessAPI.html) を参照してください。

## CLI で概要を確認する

```bash
beachx inspect outputs/latest
```

画像を保存する場合:

```bash
beachx inspect outputs/latest \
  --save-bar outputs/latest/charges_bar.png \
  --save-mesh outputs/latest/charges_mesh.png \
  --save-potential-mesh outputs/latest/potential_mesh.png
```

`sim.field_bc_mode="periodic2"` の mesh を周期セルに寄せて描く場合:

```bash
beachx inspect outputs/latest \
  --save-mesh outputs/latest/charges_mesh_periodic.png \
  --apply-periodic2-mesh
```

## 履歴アニメーションを作る

`charge_history.csv` がある場合:

```bash
beachx animate outputs/latest \
  --quantity charge \
  --save-gif outputs/latest/charge_history.gif \
  --total-frames 200
```

`potential_history.csv` を出している場合は `--quantity potential` も使えます。

## Python API で読む

```python
from beach import Beach

b = Beach("outputs/latest")
print("absorbed:", b.result.absorbed)
print("escaped:", b.result.escaped)
print("batches:", b.result.batches)

b.plot_bar()
b.plot_mesh()
b.plot_potential()
```

設定ファイルが出力ディレクトリ近傍にない場合は明示します。

```python
b = Beach("outputs/latest", config_path="beach.toml")
```

## 特定 mesh だけを見る

`mesh_sources.csv` で `mesh_id` を確認してから、対象 mesh を選べます。

```python
from beach import Beach

b = Beach("outputs/latest")
mesh1 = b.get_mesh(1)
charge1 = b.get_mesh_charge(1)
print(mesh1.centers.shape, charge1.shape)
```

履歴がある場合は、batch index を指定できます。

```python
mesh1_step10 = b.get_mesh(1, step=10)
charge1_step10 = b.get_mesh_charge(1, step=10)
```

## 断面・力・移動しやすさ

少し進んだ解析では次の CLI を使います。

```bash
beachx slices outputs/latest \
  --grid-n 200 \
  --save outputs/latest/potential_slices.png

beachx coulomb outputs/latest \
  --component z \
  --save outputs/latest/coulomb_force_z.png

beachx mobility outputs/latest \
  --density-kg-m3 2500 \
  --mu-static 0.4 \
  --save-csv outputs/latest/mobility_summary.csv
```

これらは近傍の `beach.toml` から geometry や periodic2 設定を自動解決します。
見つからない場合は、対応する CLI の `--config` を指定してください。

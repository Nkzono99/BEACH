title: 後処理チュートリアル

Lang: [日本語](PostprocessTutorial.md) | [English](PostprocessTutorial.en.md)

# 後処理チュートリアル

最初にCLIで実行結果の全体を確認し、次にPythonで分布と履歴を読み込みます。
出力ファイルの意味は [出力の読み方](OutputGuide.html)、API の詳細は [Python 後処理 API リファレンス](PythonPostprocessAPI.html) から確認できます。

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

設定ファイルが出力directoryの近傍にない場合は、`config_path`で指定します。

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

## 周期画像を含む object の離脱力を見る

保存された charge snapshot を固定し、mesh 6 だけを上向きに動かす例です。

```bash
beachx object-detachment outputs/latest \
  --config beach.toml \
  --target-mesh-id 6 \
  --periodic-model infinite-physical \
  --z-max-m 2.0e-4 \
  --z-points 65 \
  --mass-kg 2.0e-12 \
  --gravity-m-s2 9.80665 \
  --adhesion-force-n 1.0e-10 \
  --adhesion-range-m 2.0e-6 \
  --output-dir outputs/latest/object_detachment
```

`object-detachment` の CLI 既定重力は月面の `1.62 m/s^2` です。この例は地上重力を
仮定して `9.80665 m/s^2` を明示しています。対象環境に合わせて変更してください。

`configured`はrunのfinite/cached設定をそのまま使います。`infinite-physical`は、
x/y periodic runのcached `k != 0`と物理的な`k = 0` modeを使います。

targetのcentral-cell primary自己場だけを除外し、target自身の周期画像が作る力は残します。
`instantaneous_wrench.csv`、`path.csv`、`summary.json`、`report.md`が生成されます。

同じ解析を Python から行う完全な例は
[`examples/analyze_periodic_object_detachment.py`](https://github.com/Nkzono99/BEACH/blob/main/examples/analyze_periodic_object_detachment.py)
です。最小の API 呼び出しは次のとおりです。

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
    adhesion=AdhesionProfile.none(),
)
```

正常終了し、CSV/JSONが生成されれば、解析処理自体は完了しています。離脱の物理的な妥当性は、
`path.status`、仕事と電位差の不一致、mesh/quadrature、finite shellまたはperiodic cache、経路上端、
charge snapshot、stochastic seedへの依存性から評価します。[<sup>1</sup>](ValidationGuide.html)
非中性周期cellで得た有限高さのspeedは、無限遠へのescape speedではありません。

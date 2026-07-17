title: 後処理チュートリアル

Lang: [日本語](PostprocessTutorial.md) | [English](PostprocessTutorial.en.md)

# 後処理チュートリアル

BEACHの後処理は、Python package `beach`を直接使う方法と、そのpackageで定型処理を行う`beachx`があります。
このページではPython APIを先に示し、その後に用意済みの可視化・解析commandを紹介します。全class・関数は
[Python後処理APIリファレンス](PythonPostprocessAPI.html)、出力 file の意味は[出力ファイルを調べる](OutputGuide.html)から確認できます。

## Python APIで後処理する

### 結果を読み、基本図を作る

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

### 特定meshだけを見る

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

### 周期画像を含むobjectの離脱力を調べる

保存されたcharge snapshotを固定し、mesh 6だけを上向きに動かす最小例です。

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

完全な例は
[`examples/analyze_periodic_object_detachment.py`](https://github.com/Nkzono99/BEACH/blob/main/examples/analyze_periodic_object_detachment.py)
にあります。

## `beachx`で用意済みの可視化・解析を使う

`beachx`は、Python APIの代表的な後処理をcommandとしてまとめた入口です。定型図やCSVをすぐ作る場合はこちらを使います。

### 実行結果の概要と基本図

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

`sim.field_bc_mode="periodic2"`のmeshを周期cellに寄せて描く場合:

```bash
beachx inspect outputs/latest \
  --save-mesh outputs/latest/charges_mesh_periodic.png \
  --apply-periodic2-mesh
```

### 履歴アニメーション

`charge_history.csv`がある場合:

```bash
beachx animate outputs/latest \
  --quantity charge \
  --save-gif outputs/latest/charge_history.gif \
  --total-frames 200
```

`potential_history.csv`を出している場合は`--quantity potential`も使えます。

### 断面・力・移動しやすさ

少し進んだ解析では次のcommandを使います。

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

これらは近傍の`beach.toml`からgeometryやperiodic2設定を自動解決します。
見つからない場合は、対応するcommandの`--config`を指定してください。

### 周期画像を含むobjectの離脱力

上のPython API例と同じ流れを、CSV/JSON/report出力までまとめて実行します。

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

`object-detachment`の既定重力は月面の`1.62 m/s^2`です。この例は地上重力を仮定して
`9.80665 m/s^2`を明示しています。対象環境に合わせて変更してください。

`configured`はrunのfinite/cached設定をそのまま使います。`infinite-physical`は、
x/y periodic runのcached `k != 0`と物理的な`k = 0` modeを使います。

targetのcentral-cell primary自己場だけを除外し、target自身の周期画像が作る力は残します。
`instantaneous_wrench.csv`、`path.csv`、`summary.json`、`report.md`が生成されます。

正常終了し、CSV/JSONが生成されれば、解析処理自体は完了しています。離脱の物理的な妥当性は、
`path.status`、仕事と電位差の不一致、mesh/quadrature、finite shellまたはperiodic cache、経路上端、
charge snapshot、stochastic seedへの依存性から評価します。[<sup>1</sup>](ValidationGuide.html)
非中性周期cellで得た有限高さのspeedは、無限遠へのescape speedではありません。

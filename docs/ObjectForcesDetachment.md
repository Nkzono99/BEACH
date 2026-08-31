title: 物体の力と離脱を調べる

Lang: [日本語](ObjectForcesDetachment.md) | [English](ObjectForcesDetachment.en.md)

# 物体の力と離脱を調べる

この task guide は、保存済みの複数物体 run から、定位置での物体間力と移動指標を確認し、保存電荷を
固定した鉛直離脱経路を評価する手順です。最後まで進むと、定位置の力、経路に沿った仕事、重力と付着を
含む離脱判定を review 用 artifact として保存できます。

**公式入門ケースの `outputs/tutorial` には適用できません。** 公式ケースは単一の固定平面だけを含むため、
物体間力、移動しやすさ、離脱を定義できません。別の run に、複数の `mesh_id` と、その run に実際に使った
完全な `beach.toml` が必要です。

**native field kernel も必要です。** 通常の Python package install は `libbeach_field_kernel.so` を
含まないため、`inspect` より後の object force command と Python 解析はそのままでは実行できません。
Fortran compiler と `fpm` が使える BEACH source checkout で library を build し、実際の path を設定します。

```bash
make -C /path/to/BEACH build-kernel
export BEACH_FIELD_KERNEL_LIB=/path/to/BEACH/build/libbeach_field_kernel.so
```

以下では、出力を `outputs/latest`、設定を現在のディレクトリにある `beach.toml`、移動対象を `mesh_id=6` と
仮定します。`/path/to/BEACH` とこれらの値は、すべて自分の環境と run に置き換えてください。

## 1. run と移動対象を確認する

先に出力を読み、実在する物体だけを選びます。

```bash
beachx inspect outputs/latest
```

`mesh_ids` に複数の ID があり、`mesh_source` または `mesh_sources.csv` から各形状を説明できることを確認します。
対象 ID `6` がなければ、以降の `6` を実在する ID へ変更します。ここで複数物体を確認できない run は、この解析に
進めません。

x/y 周期 run の mesh を primary cell 内へ表示したい場合だけ、次を使います。

```bash
beachx inspect outputs/latest \
  --save-mesh outputs/latest/charges_mesh_periodic.png \
  --apply-periodic2-mesh
```

## 2. 保存位置での力と移動指標を確認する

まず、物体を動かさずに3種類の解析を実行します。

```bash
beachx coulomb outputs/latest \
  --config beach.toml \
  --component z \
  --save outputs/latest/coulomb_force_z.png

beachx kernel-forces outputs/latest \
  --config beach.toml \
  --target-mesh-ids 6 \
  --save-csv outputs/latest/object_forces_kernel.csv

beachx mobility outputs/latest \
  --config beach.toml \
  --density-kg-m3 2500 \
  --mu-static 0.4 \
  --save-csv outputs/latest/mobility_summary.csv
```

`coulomb` は物体間の成分別 force matrix、`kernel-forces` は field kernel による物体別合力、`mobility` は
重力・支持面・摩擦係数を加えた lift / slide / roll 指標を出します。出力 file ができたことは解析の完了を
示しますが、運動開始を証明しません。密度、摩擦、支持面、mesh、field model の前提を先に検証してください。

## 3. 保存電荷を固定した離脱 artifact を作る

離脱解析は、保存済みの source geometry と電荷を初期位置に固定し、選んだ central-cell target だけを鉛直に
動かします。物体が動く間の再帯電や周囲物体の移動を再計算する力学 simulation ではありません。

場の定義を先に選びます。

| `--periodic-model` | 使用条件 |
| --- | --- |
| `configured` | run の free / finite / cached 設定をそのまま使う |
| `infinite-physical` | 対応する x/y periodic run で、cached `k != 0` と物理的な `k = 0` mode を合成する |

次の例は `infinite-physical` を使うため、対応する periodic run と cache が必要です。

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

正常終了すると、`instantaneous_wrench.csv`、`path.csv`、`summary.json`、`report.md` が生成されます。
`object-detachment` の既定重力は月面の `1.62 m/s^2` で、この例は地上重力を明示しています。付着を使う場合は、
`--adhesion-force-n` と `--adhesion-range-m` を必ず組にします。

`configured` は run の有限画像または cached 設定を使います。`infinite-physical` は cached `k != 0` と物理的な
`k = 0` mode を使います。target の central-cell primary 自己場だけを除外し、target 自身の周期画像が作る力は
残します。

## 4. Python で離脱経路を拡張する

独自の変位列や付着モデルを使う場合は Python API に進みます。次の例は同じ保存電荷を固定し、mesh 6 だけを
上向きに動かします。

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

print("force [N]:", wrench.force_N)
print("path status:", path.status)
print("barrier free:", release.barrier_free_from_rest)
```

完全な例は
[`examples/analyze_periodic_object_detachment.py`](https://github.com/Nkzono99/BEACH/blob/main/examples/analyze_periodic_object_detachment.py)
にあります。全 class、結果 field、`configured` と `infinite_physical` の詳細は
[Python 後処理 API リファレンス](PythonPostprocessAPI.html#104-objectinteractionsnapshot-と凍結-source-経路)を
参照してください。

## 5. 完了と妥当性を分けて判定する

zero exit status と artifact の生成は、解析処理が完了したことだけを示します。少なくとも次を確認します。

- `path.status` と、力の積分仕事・電位差から求めた仕事の不一致
- mesh と quadrature、有限画像 shell または periodic cache への依存性
- 経路上端と使用した charge snapshot への依存性
- stochastic seed を変えたときの力と離脱判定の変動
- `mobility` に与えた密度、摩擦、支持面と、離脱に与えた質量、重力、付着の妥当性

[計算結果の妥当性確認](ValidationGuide.html)に従って許容値を決めてください。非中性 periodic cell の有限高さで
得た speed は、無限遠への escape speed ではありません。

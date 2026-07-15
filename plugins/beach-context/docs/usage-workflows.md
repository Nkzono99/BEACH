# BEACH usage workflows

## 新規 case

```bash
mkdir beach-tutorial
cd beach-tutorial
beachx config init
beachx lint beach.toml
beach beach.toml
```

`beach.toml` は日常編集用であり、Fortran runtime が直接読む設定ファイルでもあります。高水準記法を使った場合は Fortran parser が読み込み時に正規化します。
生成内容は`references/examples/tutorial_insulator.toml`と同じで、x/y periodic、
`field_solver="fmm"`、`field_periodic_image_layers=1`、
`field_periodic_far_correction="none"`の$3\times3$有限画像和です。

## runtime config の直接実行

```bash
beach examples/beach.toml
```

引数なしの場合、カレントディレクトリの `beach.toml` を読みます。開発 checkout では `make` 後に `fpm run --profile release --flag "-fopenmp" -- examples/beach.toml` も使えます。

## 出力確認

```bash
beachx inspect outputs/latest
beachx animate outputs/latest --quantity charge --save-gif outputs/latest/charge_history.gif
beachx slices outputs/latest --grid-n 200 --save outputs/latest/potential_slices.png
beachx coulomb outputs/latest --component z --save outputs/latest/coulomb_force_z.png
beachx mobility outputs/latest --density-kg-m3 2500 --mu-static 0.4 --save-csv outputs/latest/mobility_summary.csv
```

## Python 解析

```python
from beach import Beach

run = Beach("outputs/latest")
print(run.result.absorbed, run.result.last_rel_change)
run.plot_mesh()
run.plot_potential()
```

## workload estimate

```bash
beachx estimate-workload beach.toml --threads 8
```

`reservoir_face` や `photo_raycast` では `batch_duration` / `batch_duration_step`、macro particle 数、mesh 要素数、history 出力が計算量と出力量を強く左右します。

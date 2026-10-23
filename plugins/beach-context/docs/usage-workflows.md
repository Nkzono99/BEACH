# BEACH usage workflows

## 新規 case

```bash
mkdir run_periodic2
cd run_periodic2
beachx config init
beachx lint beach.toml
beach beach.toml
```

`beach.toml` は日常編集用であり、Fortran runtime が直接読む設定ファイルでもあります。高水準記法を使った場合は Fortran parser が読み込み時に正規化します。

## runtime config の直接実行

```bash
beach examples/beach.toml
```

引数なしの場合、カレントディレクトリの `beach.toml` を読みます。開発 checkout では `make` 後に `fpm run --target beach --profile release --flag "-fopenmp" -- examples/beach.toml` も使えます。

## 出力確認

```bash
beachx inspect outputs/latest
beachx animate outputs/latest --quantity charge --save outputs/latest/charge_history.gif
beachx slices outputs/latest --grid-n 200 --save outputs/latest/potential_slices.png
beachx coulomb outputs/latest --component z --save outputs/latest/coulomb_force_z.png
beachx mobility outputs/latest --density-kg-m3 2500 --mu-static 0.4 --save outputs/latest/mobility_summary.csv
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
beachx workload beach.toml --threads 8
```

`inflow` や `photo_raycast` では `run.batch.duration_s` / `run.batch.duration_steps`、macro particle 数、mesh 要素数、history 出力が計算量と出力量を強く左右します。

## periodic2 の遠方演算子 cache を別 job で作る

`fields.periodic.backend="cached_kneq0"` は、初めての box・展開次数・Ewald 設定の実行で、初期化中に遠方演算子を作ります。
規模によっては数分〜数十分かかり、本番 job の時間枠と rank 数（再開時に固定）を消費します。cache は別 job で並列に作ります。

```bash
export BEACH_CACHE_DIR=/path/to/shared/beach-cache   # 計算ノードから見える共有 FS。全 job で同じ場所
srun -n 16 beach --prepare-periodic-cache beach.toml # cache job: 本番より多い rank で生成だけ
sbatch --dependency=afterok:<cache job id> run.sh    # 本番 job は cache job の後
```

- 本番の `summary.txt` で `periodic2_cache_hit=T` と `periodic2_operator_build_count=0` を確認する。
- 同じ box・設定なら物体集合が違っても同じ cache を使い、足りないセルだけが追加される。物体集合を変える sweep では、本番の前に全物体集合の `--prepare-periodic-cache` を 1 つの job で順に流す。
- 同じ cache への生成は lock で直列化される。cache のない本番 job を並べて投げると、待つ job の時間枠を浪費する。
- box の寸法、`fields.periodic.image_layers`、`ewald_layers`、展開次数を変えると別の cache になる。

# Fortran中心ワークフロー

このリポジトリでは、v0.x 以降の実行系を **Fortran（fpm）中心** とし、Python は後処理と可視化に利用します。

## 1. 実行フロー

1. Fortran 実行ファイル（`app/main.f90`）を `fpm` で実行
2. 計算結果を `outputs/...` に保存（`summary.txt`, `charges.csv`, `mesh_triangles.csv`, `rng_state.txt`）し、履歴は `charge_history.csv` に逐次追記
3. Python スクリプトで結果を読み込み、統計確認や可視化を実施

## 2. Fortran 実行

```bash
fpm run --profile release --flag "-fopenmp"
```

設定ファイル指定:

```bash
fpm run --profile release --flag "-fopenmp" -- examples/fortran_config.toml
```

OpenMP スレッド数指定:

```bash
OMP_NUM_THREADS=8 fpm run --profile release --flag "-fopenmp" -- examples/fortran_config.toml
```

再開実行:

```toml
[output]
dir = "outputs/latest"
resume = true
```

同じ設定で再度実行すると、`outputs/latest` に保存された前回の状態から `sim.batch_count` 分だけ続けて計算します。

## 3. Python 後処理

依存導入:

```bash
python -m pip install -U pip setuptools wheel
python -m pip install -e . --no-build-isolation
```

推奨API（Python）:

```python
from beach import Beach

beach = Beach("outputs/latest")
print(beach.result.absorbed, beach.result.escaped)

beach.plot_bar()
beach.plot_mesh()
beach.plot_potential()
beach.animate_mesh("outputs/latest/charge_history.gif", quantity="charge")
```

`output_path=None` を指定すると、`beach.animate_mesh()` は保存せず `FuncAnimation` を返します。

既存の関数API（`plot_charges`, `plot_charge_mesh`, `compute_potential_mesh`, `plot_potential_mesh`, `animate_history_mesh`）も継続利用可能で、`Beach` または `FortranRunResult` のどちらを渡しても動作します。

出力確認（CLI）:

```bash
python examples/inspect_fortran_output.py outputs/latest \
  --save-bar outputs/latest/charges_bar.png \
  --save-mesh outputs/latest/charges_mesh.png \
  --save-potential-mesh outputs/latest/potential_mesh.png \
  --potential-self-term area-equivalent

python examples/animate_fortran_history.py outputs/latest \
  --quantity charge \
  --save-gif outputs/latest/charge_history.gif

python examples/animate_fortran_history.py outputs/latest \
  --quantity potential \
  --save-gif outputs/latest/potential_history.gif \
  --potential-self-term area-equivalent
```

Python 後処理では、`charge mesh` に加えて `potential mesh` も生成できます。これは Fortran を再実行して電位を解き直すのではなく、出力済みの要素電荷から電位を再構成する近似です。既定では他要素の寄与を各要素重心の点電荷近似で扱い、自己項のみ要素面積に基づく有限値（面積等価円板近似）で評価します。そのため、Fortran 本体の電場計算と数値的一致を保証するものではありません。旧来の `1 / softening` 自己項を使いたい場合は、`--potential-self-term softened-point --potential-softening 1.0e-6` または `Beach.compute_potential(..., softening=1.0e-6, self_term="softened_point")`（または `compute_potential_mesh(...)`）を指定してください。`--potential-softening` の既定値は `0.0` で、`area-equivalent` / `exclude` では主に要素間カーネルの平滑化にのみ使われます。

## 4. 役割分担の目安

- **Fortran**
  - 粒子追跡の本計算
  - メッシュ生成/読込（テンプレート, OBJ）
  - 結果ファイル出力
- **Python**
  - 結果ファイルの読込・集計
  - 可視化（バー/3D 表示、電位メッシュ再構成）
  - 実験的なデータ解析スクリプト

## 5. 運用ルール（推奨）

- 新しい物理モデルや計算ロジックは、まず Fortran 側へ実装。
- 出力フォーマットを変更した場合は、Python の `beach/fortran_results.py` と `examples/inspect_fortran_output.py` も同時更新。
- パラメータ追加時は `docs/fortran_parameter_file.md` を更新。

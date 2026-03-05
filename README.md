# BEACH（BEM + Accumulated CHarge）

本リポジトリは、**BEACH（BEM + Accumulated CHarge）** として開発している、BEM（境界要素法）ベースの表面帯電＋テスト粒子シミュレータです。  
**メイン計算系は Fortran（fpm）** とし、**Python は可視化・後処理・補助解析**に使う構成です。

## 開発・運用方針（v0.x）

- 電場計算・粒子前進・衝突判定・電荷蓄積などのコア処理は Fortran 側で実行
- Python 側は結果読込、可視化、検証用ユーティリティを担当
- Fortran 実行結果は `summary.txt` / `charges.csv` / `mesh_triangles.csv` を出力し、`charge_history.csv` はバッチ間隔指定で逐次書き出し

## ディレクトリ構成

- [`src/`](src/), [`app/`](app/)：Fortran 本体（fpm プロジェクト）
- [`examples/`](examples/)：Fortran 設定例、Python 後処理例
- [`beach/`](beach/)：Python ライブラリ（結果読込・可視化などの後処理）
- [`docs/`](docs/)：運用・仕様ドキュメント

## セットアップ

### 1) Fortran 実行（メイン）

```bash
fpm run --profile release --flag "-fopenmp"
```

設定ファイルを指定して実行する場合：

```bash
fpm run --profile release --flag "-fopenmp" -- examples/fortran_config.toml
```

### 2) Python 環境（後処理・可視化）

```bash
python -m pip install -U pip setuptools wheel
python -m pip install -e . --no-build-isolation
```

## ドキュメント

- Fortran 中心ワークフロー: [`docs/fortran_workflow.md`](docs/fortran_workflow.md)
- Fortran パラメータファイル仕様: [`docs/fortran_parameter_file.md`](docs/fortran_parameter_file.md)
- シミュレータ仕様（v0.1）: [`SPEC.md`](SPEC.md)

## Python 後処理の利用例

`Beach` クラスを使うのが推奨導線です。

```python
from beach import (
    Beach,
    animate_history_mesh,
    compute_potential_mesh,
    plot_charge_mesh,
    plot_charges,
)

beach = Beach("outputs/latest")
result = beach.result  # FortranRunResult（遅延ロード）
print(result.absorbed, result.escaped, result.charges.sum())
print(result.charge_history.shape if result.charge_history is not None else 0)

fig_bar, ax_bar = beach.plot_bar()
fig_mesh, ax_mesh = beach.plot_mesh()
fig_phi, ax_phi = beach.plot_potential()
written = beach.animate_mesh("outputs/latest/charge_history.gif", quantity="charge")
print(written)

# output_path=None なら保存せず FuncAnimation を返す
anim = beach.animate_mesh(quantity="potential")

# 既存の関数APIも継続利用可能（Beach/FortranRunResult 両対応）
potential = compute_potential_mesh(beach, softening=0.0, self_term="area_equivalent")
print(potential.min(), potential.max())
fig_bar2, ax_bar2 = plot_charges(beach)
fig_mesh2, ax_mesh2 = plot_charge_mesh(beach)
written2 = animate_history_mesh(beach, "outputs/latest/potential_history.gif", quantity="potential")
```

`Beach.compute_potential()` / `compute_potential_mesh()` は、Fortran が出力した要素電荷を各要素重心ベースで再構成する後処理近似です。既定では、他要素の寄与は重心点電荷近似のまま、自己項のみ要素面積に基づく有限値（面積等価円板近似）で評価します。これは Fortran 本体の電場計算を厳密に再現するものではありません。旧来の `1 / softening` 自己項が必要な場合は、`beach.compute_potential(softening=1.0e-6, self_term="softened_point")` のように明示指定してください。

CLI での確認例：

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

## 参考ファイル

- [`examples/fortran_config.toml`](examples/fortran_config.toml)：複数テンプレート合成を含む設定サンプル

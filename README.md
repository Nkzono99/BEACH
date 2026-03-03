# BEACH（BEM + Accumulated CHarge）

本リポジトリは、**BEACH（BEM + Accumulated CHarge）** として開発している、BEM（境界要素法）ベースの表面帯電＋テスト粒子シミュレータです。  
**メイン計算系は Fortran（fpm）** とし、**Python は可視化・後処理・補助解析**に使う構成です。

## 開発・運用方針（v0.x）

- 電場計算・粒子前進・衝突判定・電荷蓄積などのコア処理は Fortran 側で実行
- Python 側は結果読込、可視化、検証用ユーティリティを担当
- Fortran 実行結果は `summary.txt` / `charges.csv` / `mesh_triangles.csv` を出力し、`charge_history.csv` はバッチ間隔指定で逐次書き出し

## ディレクトリ構成

- `src/`, `app/`：Fortran 本体（fpm プロジェクト）
- `examples/`：Fortran 設定例、Python 後処理例
- `beach/`：Python ライブラリ（後処理 + 既存プロトタイプ資産）
- `docs/`：運用・仕様ドキュメント

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

- Fortran 中心ワークフロー: `docs/fortran_workflow.md`
- Fortran パラメータファイル仕様: `docs/fortran_parameter_file.md`
- シミュレータ仕様（v0.1）: `SPEC.md`

## Python 後処理の利用例

Fortran 出力を Python で読み込む例：

```python
from beach import load_fortran_result, plot_charges, plot_charge_mesh

result = load_fortran_result("outputs/latest")
print(result.absorbed, result.escaped, result.charges.sum())
print(result.charge_history.shape if result.charge_history is not None else None)

fig_bar, ax_bar = plot_charges(result)
fig_mesh, ax_mesh = plot_charge_mesh(result)
```

CLI での確認例：

```bash
python examples/inspect_fortran_output.py outputs/latest \
  --save-bar outputs/latest/charges_bar.png \
  --save-mesh outputs/latest/charges_mesh.png
```

## 参考ファイル

- `examples/fortran_config.toml`：複数テンプレート合成を含む設定サンプル
- `examples/visualize_bem_list_3d.py`：CSV からの 3D 可視化サンプル

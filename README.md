# bemtracer

BEM（Boundary Element Method）ベースの表面帯電＋テスト粒子シミュレータです。**メイン実装は Fortran（fpm）**、Python は可視化・後処理・補助スクリプトに利用する構成へ整理しています。

## 方針（v0.x）

- コアシミュレーション（電場計算・粒子前進・衝突・電荷蓄積）は Fortran。
- Python 実装は互換リファレンスおよび解析ユーティリティとして維持。
- 出力は Fortran 側で `summary.txt` / `charges.csv` / `mesh_triangles.csv` を生成し、Python で読み込み・可視化。

## リポジトリ構成

- `src/`, `app/`: Fortran 本体（fpm）
- `examples/`: Fortran 設定例と Python 後処理例
- `bemtracer/`: Python ライブラリ（後処理 + 既存プロトタイプ）
- `docs/`: 設計・運用ドキュメント

## セットアップ

### Fortran 実行環境（メイン）

```bash
fpm run --profile release --flag "-fopenmp"
```

設定ファイルを渡す場合:

```bash
fpm run --profile release --flag "-fopenmp" -- examples/fortran_config.toml
```

### Python 環境（後処理・補助）

```bash
python -m pip install -U pip setuptools wheel
python -m pip install -e . --no-build-isolation
```

## ドキュメント

- Fortran ワークフロー: `docs/fortran_workflow.md`
- パラメータファイル仕様: `docs/fortran_parameter_file.md`
- シミュレータ仕様（v0.1）: `SPEC.md`

## Python 後処理の例

Fortran 出力の読み込み:

```python
from bemtracer import load_fortran_result, plot_charges, plot_charge_mesh

result = load_fortran_result("outputs/latest")
print(result.absorbed, result.escaped, result.charges.sum())

fig_bar, ax_bar = plot_charges(result)
fig_mesh, ax_mesh = plot_charge_mesh(result)
```

CLI 例:

```bash
python examples/inspect_fortran_output.py outputs/latest \
  --save-bar outputs/latest/charges_bar.png \
  --save-mesh outputs/latest/charges_mesh.png
```

## 参考

- `examples/fortran_config.toml`: 複数テンプレート合成を含む設定例
- `examples/visualize_bem_list_3d.py`: CSV ベースの 3D 可視化例

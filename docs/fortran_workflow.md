# Fortran中心ワークフロー

このリポジトリでは、v0.x 以降の実行系を **Fortran（fpm）中心** とし、Python は後処理と可視化に利用します。

## 1. 実行フロー

1. Fortran 実行ファイル（`app/main.f90`）を `fpm` で実行
2. 計算結果を `outputs/...` に保存（`summary.txt`, `charges.csv`, `mesh_triangles.csv`）し、履歴は `charge_history.csv` に逐次追記
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

## 3. Python 後処理

依存導入:

```bash
python -m pip install -U pip setuptools wheel
python -m pip install -e . --no-build-isolation
```

出力確認（CLI）:

```bash
python examples/inspect_fortran_output.py outputs/latest \
  --save-bar outputs/latest/charges_bar.png \
  --save-mesh outputs/latest/charges_mesh.png
```

## 4. 役割分担の目安

- **Fortran**
  - 粒子追跡の本計算
  - メッシュ生成/読込（テンプレート, OBJ）
  - 結果ファイル出力
- **Python**
  - 結果ファイルの読込・集計
  - 可視化（バー/3D 表示）
  - 実験的なデータ解析スクリプト

## 5. 運用ルール（推奨）

- 新しい物理モデルや計算ロジックは、まず Fortran 側へ実装。
- 出力フォーマットを変更した場合は、Python の `bemtracer/fortran_results.py` と `examples/inspect_fortran_output.py` も同時更新。
- パラメータ追加時は `docs/fortran_parameter_file.md` を更新。

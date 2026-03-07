# BEACH（BEM + Accumulated CHarge）

BEACH は、**BEM（境界要素法）ベースの表面帯電 + テスト粒子追跡シミュレータ**です。  
計算本体は Fortran（`fpm` 管理）、Python は後処理・可視化を担当します。

v0.x の主対象は **insulator accumulation（絶縁体への電荷蓄積）** です。

## 推奨ワークフロー

このリポジトリでは、**インストール済みの `beach` コマンドで実行する運用**を推奨します。  
`fpm run` は開発時の直接実行手段として残しています。

## 1. セットアップ

### 1.1 前提ツール

```bash
gfortran --version
fpm --version
python --version
```

### 1.2 Fortran 実行系のインストール（推奨）

```bash
make
```

`make` のデフォルトターゲットは `install` で、`BUILD_PROFILE=auto`（ホスト自動判定）で導入します。
環境に応じてプロファイルを明示することも可能です。

```bash
make install-generic
make install-camphor
```

`install.sh` の既定インストール先は `~/.local/bin/beach` です。必要に応じて PATH を通してください。

```bash
export PATH="$HOME/.local/bin:$PATH"
```

### 1.3 Python 後処理

```bash
python -m pip install -U pip setuptools wheel
python -m pip install -e . --no-build-isolation
```

## 2. 実行

### 2.1 推奨: `beach` コマンド

```bash
beach examples/beach.toml
```

引数なしの場合は、カレントディレクトリの `beach.toml` を自動読込します。

### 2.2 出力確認

```bash
ls outputs/latest
beach-inspect outputs/latest
```

主な出力:

- `summary.txt`
- `charges.csv`
- `mesh_triangles.csv`
- `charge_history.csv`（`history_stride > 0` のとき）

### 2.3 可視化

```bash
beach-inspect outputs/latest \
  --save-bar outputs/latest/charges_bar.png \
  --save-mesh outputs/latest/charges_mesh.png

beach-animate-history outputs/latest \
  --quantity charge \
  --save-gif outputs/latest/charge_history.gif \
  --total-frames 200
```

## 3. よく使うコマンド

### テスト

```bash
pytest -q
fpm test
```

### 負荷見積もり

```bash
beach-estimate-workload examples/beach.toml --threads 8
```

### 再開実行

```toml
[output]
dir = "outputs/latest"
resume = true
```

同じ `output.dir` で `beach` を再実行すると続きから計算します。

## 4. 開発者向け補足（`fpm`）

このプロジェクトは `fpm` で管理されています。開発時には直接実行も可能です。

```bash
fpm run --profile release --flag "-fopenmp" -- examples/beach.toml
```

MPI + OpenMP（`USE_MPI` 有効化）:

```bash
FPM_FC=mpiifort \
fpm run --profile release --flag "-fpp -DUSE_MPI -qopenmp" \
  --runner "mpirun -n 4" -- examples/beach.toml
```

## 5. プロジェクト構成

- [`src/`](src/), [`app/`](app/): Fortran 本体
- [`tests/fortran/`](tests/fortran/): Fortran テスト
- [`beach/`](beach/): Python 後処理ライブラリ
- [`tests/python/`](tests/python/): Python テスト
- [`examples/`](examples/): 設定例・補助スクリプト
- [`docs/`](docs/): 運用ドキュメント
- [`SPEC.md`](SPEC.md): 現行 Fortran 実装仕様

## 6. ドキュメント案内

- Fortran 実行フロー: [`docs/fortran_workflow.md`](docs/fortran_workflow.md)
- `beach.toml` 仕様: [`docs/fortran_parameter_file.md`](docs/fortran_parameter_file.md)
- 実装仕様（source of truth）: [`SPEC.md`](SPEC.md)
- 設定サンプル: [`examples/beach.toml`](examples/beach.toml)

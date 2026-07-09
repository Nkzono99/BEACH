# BEACH（BEM + Accumulated CHarge）

[![GitHub Pages](https://img.shields.io/website?url=https%3A%2F%2Fnkzono99.github.io%2FBEACH%2F&up_message=GitHub%20Pages&down_message=Pages%20down)](https://nkzono99.github.io/BEACH/)
[![Fortran Docs](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-docs.yml/badge.svg?branch=main)](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-docs.yml)
[![Fortran Format](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-format.yml/badge.svg?branch=main)](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-format.yml)
[![PyPI version](https://img.shields.io/pypi/v/beach-bem)](https://pypi.org/project/beach-bem/)
[![Python versions](https://img.shields.io/pypi/pyversions/beach-bem)](https://pypi.org/project/beach-bem/)
[![License: Apache 2.0](https://img.shields.io/badge/license-Apache%202.0-green.svg)](LICENSE)

BEACH は、**境界要素法（BEM）による表面電場計算**と
**テスト粒子追跡**を組み合わせた表面帯電シミュレータです。

- `beach`: Fortran 製のシミュレーション実行バイナリ
- `beachx`: 設定検査、可視化、ワークロード見積もりなどの Python CLI
- `beach/`: Fortran 出力を読むための Python ライブラリ

現行の主対象は **絶縁体表面への電荷蓄積（insulator accumulation）** です。
各 batch で電場計算、Boris pusher による粒子前進、メッシュ衝突判定、表面への電荷堆積を行い、表面電位の時間発展を出力します。

<div align="center">
  <img src="docs/images/potential_history.gif" alt="帯電シミュレーションの電位変化" width="80%">
  <p><i>電子ビーム照射下での絶縁体メッシュ上の電位分布の時間発展</i></p>
  <sub>3D model: <a href="https://www.turbosquid.com/ja/3d-models/rubber-duck-pbr-game-ready-model-2001526">Rubber Duck PBR Game Ready</a> (TurboSquid)</sub>
</div>

## 必要なもの

Python パッケージのインストール時に Fortran 実行バイナリもビルドします。
事前に次のコマンドが使える状態にしてください。

```bash
make --version
gfortran --version
fpm --version
python --version
```

`gfortran` 以外の Fortran compiler を使う場合は、環境に合わせて `FC` / `FPM_FC` などを設定してください。

## クイックスタート

通常利用では PyPI から導入します。

```bash
python -m pip install -U pip setuptools wheel
python -m pip install beach-bem
```

ユーザーインストール時に `beach` / `beachx` が見つからない場合は PATH を確認してください。

```bash
export PATH="$HOME/.local/bin:$PATH"
```

小さな設定を作って実行します。

```bash
mkdir run_periodic2
cd run_periodic2

beachx config init
beachx lint beach.toml
beach beach.toml
```

結果は既定で `outputs/latest/` に出力されます。

```bash
beachx inspect outputs/latest
beachx animate outputs/latest --quantity potential --save-gif potential_history.gif
```

開発版を直接試す場合:

```bash
python -m pip install "git+https://github.com/Nkzono99/BEACH.git"
```

## 設定ファイル

実行時に `beach` が読むのは `beach.toml` です。まずは `beachx config init` が生成したファイルを編集し、実行前に `beachx lint beach.toml` で検査してください。

`beachx config render` は、`box_origin` / `box_size`、`inject_region_mode = "face_fraction"`、`placement_mode`、`mesh.groups` などの高水準記法を、Fortran 実行系が読む最終キーへ展開します。
出力先を指定しない場合は入力ファイルを上書きするため、内容確認には `--stdout` か `--output` を使ってください。

```bash
beachx config render --stdout
beachx config render beach.toml --output beach.rendered.toml
```

詳細:

- 設定作成・検証: [`docs/Configuration.md`](docs/Configuration.md)
- `beach.toml` 仕様: [`docs/Parameters.md`](docs/Parameters.md)
- 実行・開発・テスト: [`docs/Workflow.md`](docs/Workflow.md)

## よく使うコマンド

| コマンド | 用途 |
| --- | --- |
| `beach beach.toml` | Fortran シミュレーションを実行 |
| `beachx config init [path]` | 小さな実行可能設定を作成 |
| `beachx lint beach.toml` | TOML / JSON Schema / BEACH 制約を検査 |
| `beachx config render beach.toml` | 高水準記法を最終 `beach.toml` キーへ展開 |
| `beachx inspect outputs/latest` | 出力ディレクトリの概要表示 |
| `beachx animate outputs/latest` | 電荷・電位履歴のアニメーション生成 |
| `beachx workload beach.toml --threads 8` | 実行前のワークロード見積もり |
| `beachx model close-pack` | 密充填球モデルの設定生成 |

旧 CLI 名（`beach-inspect`、`beach-animate-history`、`beach-estimate-workload` など）は互換用 alias です。新しい利用者には `beachx ...` を推奨します。

## ドキュメント

- 公開ドキュメント: [GitHub Pages](https://nkzono99.github.io/BEACH/)
- ドキュメント一覧: [`docs/index.md`](docs/index.md)
- アルゴリズム解説: [`docs/Algorithms.md`](docs/Algorithms.md)
- Python 後処理 API / CLI: [`docs/PythonPostprocessAPI.md`](docs/PythonPostprocessAPI.md)
- Fortran 実装仕様: [`SPEC.md`](SPEC.md)
- 変更履歴: [`CHANGELOG.md`](CHANGELOG.md)

## リポジトリ構成

- [`src/`](src/), [`app/`](app/): Fortran 本体
- [`beach/`](beach/): Python 後処理ライブラリと `beachx`
- [`examples/`](examples/): 設定例・補助スクリプト
- [`tests/fortran/`](tests/fortran/), [`tests/python/`](tests/python/): テスト
- [`docs/`](docs/): 利用者・開発者向けドキュメント
- [`schemas/`](schemas/), [`beach/config/schemas/`](beach/config/schemas/): `beach.toml` JSON Schema

## 開発者向け

ローカル checkout を編集する場合:

```bash
python -m pip install -e . --no-build-isolation
make check
```

Python 側の確認:

```bash
pytest -q
ruff check .
```

Fortran を含むテストは時間がかかるため、目的に応じて tiered target を選んでください。

```bash
make test-l0
make test-l1
make test-l2
make test-l3
```

KUDPC のログインノードでは `make test*` / `fpm test` / 長時間の Fortran 実行を直接走らせず、計算ノードへ投入してください。詳細は [`docs/Workflow.md`](docs/Workflow.md) を参照してください。

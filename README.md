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
各バッチで電場計算、Boris pusher による粒子前進、メッシュ衝突判定、表面への電荷堆積を行い、表面電位の時間発展を出力します。

<div align="center">
  <img src="https://nkzono99.github.io/BEACH/images/potential_history.gif" alt="帯電シミュレーションの電位変化" width="80%">
  <p><i>電子ビーム照射下での絶縁体メッシュ上の電位分布の時間発展</i></p>
  <sub>3D model: <a href="https://www.turbosquid.com/ja/3d-models/rubber-duck-pbr-game-ready-model-2001526">Rubber Duck PBR Game Ready</a> (TurboSquid)</sub>
</div>

## 必要なもの

Python パッケージのインストール時に、Fortran 実行バイナリもビルドします。通常の
`pip install` では `fpm` も隔離ビルド環境へ自動導入されるため、事前に必要なのは `make` と
Fortran コンパイラです。

```bash
make --version
gfortran --version
python --version
```

`gfortran` 以外の Fortran コンパイラを使う場合は、環境に合わせて `FC` / `FPM_FC` などを設定してください。

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
mkdir beach-tutorial
cd beach-tutorial

beachx config init
beachx lint beach.toml
beach beach.toml
```

この入門ケースでは、1 個の電子が平面へ吸収され、衝突した要素へ負電荷が堆積します。
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

`beach` は実行時に `beach.toml` を読みます。まずは `beachx config init` で生成した設定を編集し、
実行前に `beachx lint beach.toml` で検査してください。個別の設定項目と制約は、次のガイドにまとめています。

- 設定の作成と検証: [設定を編集する](https://nkzono99.github.io/BEACH/configuration.html)
- 全設定項目の仕様: [`beach.toml` パラメータ仕様](https://nkzono99.github.io/BEACH/parameters.html)

## よく使うコマンド

| コマンド | 用途 |
| --- | --- |
| `beach beach.toml` | Fortran シミュレーションを実行 |
| `beachx config init [path]` | 小さな実行可能設定を作成 |
| `beachx lint beach.toml` | TOML / JSON Schema / BEACH 制約を検査 |
| `beachx inspect outputs/latest` | 出力ディレクトリの概要表示 |
| `beachx animate outputs/latest` | 電荷・電位履歴のアニメーション生成 |
| `beachx workload beach.toml --threads 8` | 実行前のワークロード見積もり |
| `beachx model close-pack` | 密充填球モデルの設定生成 |

旧 CLI 名（`beach-inspect`、`beach-animate-history`、`beach-estimate-workload` など）は互換用の別名です。新しい利用者には `beachx ...` を推奨します。

## ドキュメント

利用者向け・開発者向けドキュメントは [BEACH ドキュメント](https://nkzono99.github.io/BEACH/) に集約しています。

- 初めて実行する: [インストール](https://nkzono99.github.io/BEACH/installation.html) → [10 分チュートリアル](https://nkzono99.github.io/BEACH/tutorial.html) → [出力の読み方](https://nkzono99.github.io/BEACH/output-guide.html)
- 研究ケースを作る: [設定レシピ](https://nkzono99.github.io/BEACH/configuration-recipes.html) → [計算結果の妥当性確認](https://nkzono99.github.io/BEACH/validation-guide.html)
- open境界と局所returnを設定する: [粒子のescapeと局所return](https://nkzono99.github.io/BEACH/particle-escape-return.html)。
- 数値モデルを確認する: [アルゴリズム概要](https://nkzono99.github.io/BEACH/algorithms.html) / [場ソルバー](https://nkzono99.github.io/BEACH/field-solvers.html) / [粒子更新](https://nkzono99.github.io/BEACH/particle-tracking-collision.html) / [FMM](https://nkzono99.github.io/BEACH/fmm-core.html)
- 開発する: [Fortran 中心ワークフロー](https://nkzono99.github.io/BEACH/workflow.html) / [Fortran API](https://nkzono99.github.io/BEACH/fortran/) / [Python 後処理 API](https://nkzono99.github.io/BEACH/python-postprocess-api.html)
- 問題を調べる: [トラブルシューティング](https://nkzono99.github.io/BEACH/troubleshooting.html) / [`batch_duration` の安定性](https://nkzono99.github.io/BEACH/batch-duration-stability.html)
- 英語版を読む: [BEACH Documentation](https://nkzono99.github.io/BEACH/en.html)

Pages ビルド:

- `make docs-deps`: Starlight 依存を `npm ci` でインストール
- `make docs-starlight`: Starlight 本体のみを生成
- `make docs-pages`: Starlight 本体と FORD サブサイトを `build/starlight-site/` に生成

## リポジトリ構成

- [`src/`](src/), [`app/`](app/): Fortran 本体
- [`beach/`](beach/): Python 後処理ライブラリと `beachx`
- [`docs-site/`](docs-site/): GitHub Pages 用 Starlight サイト
- [`examples/`](examples/): 設定例・補助スクリプト
- [`tests/fortran/`](tests/fortran/), [`tests/python/`](tests/python/): テスト
- [`docs/`](docs/): 利用者・開発者向けドキュメント
- [`schemas/`](schemas/), [`beach/config/schemas/`](beach/config/schemas/): `beach.toml` JSON Schema

## 開発者向け

ローカルチェックアウトを編集する場合は、次の最小構成で準備できます。

```bash
python -m pip install fpm
python -m pip install -e . --no-build-isolation
make check
```

`--no-build-isolation` や `make` を直接使う開発では、`fpm` を PATH 上へ用意してください。

Python 側の確認:

```bash
pytest -q
ruff check .
```

Fortran を含むテストは時間がかかるため、目的に応じて階層別のターゲットを選んでください。

```bash
make test-l0
make test-l1
make test-l2
make test-l3
```

KUDPC のログインノードでは `make test*` / `fpm test` / 長時間の Fortran 実行を直接走らせず、計算ノードへ投入してください。詳細は [Fortran 中心ワークフロー](https://nkzono99.github.io/BEACH/workflow.html) を参照してください。

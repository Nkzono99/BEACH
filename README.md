# BEACH（BEM + Accumulated CHarge）

[![GitHub Pages](https://img.shields.io/website?url=https%3A%2F%2Fnkzono99.github.io%2FBEACH%2F&up_message=GitHub%20Pages&down_message=Pages%20down)](https://nkzono99.github.io/BEACH/)
[![Fortran Docs](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-docs.yml/badge.svg?branch=main)](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-docs.yml)
[![Fortran Format](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-format.yml/badge.svg?branch=main)](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-format.yml)
[![PyPI version](https://img.shields.io/pypi/v/beach-bem)](https://pypi.org/project/beach-bem/)
[![Python versions](https://img.shields.io/pypi/pyversions/beach-bem)](https://pypi.org/project/beach-bem/)
[![License: Apache 2.0](https://img.shields.io/badge/license-Apache%202.0-green.svg)](LICENSE)

BEACH は、**BEM（境界要素法）ベースの表面帯電 + テスト粒子追跡シミュレータ**です。
計算本体は Fortran（`fpm` 管理）、Python は後処理・可視化を担当します。

現行の主対象は **insulator accumulation（絶縁体への電荷蓄積）** です。BEM による電荷堆積とテスト粒子追跡を batch ごとに繰り返し、表面電位の時間発展を計算します。

<div align="center">
  <img src="docs/images/potential_history.gif" alt="帯電シミュレーションの電位変化" width="80%">
  <p><i>電子ビーム照射下での絶縁体メッシュ上の電位分布の時間発展</i></p>
  <sub>3D model: <a href="https://www.turbosquid.com/ja/3d-models/rubber-duck-pbr-game-ready-model-2001526">Rubber Duck PBR Game Ready</a> (TurboSquid)</sub>
</div>

## クイックスタート

必要なツール:

```bash
make --version
gfortran --version
fpm --version
python --version
```

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
beachx config render
beach beach.toml
beachx inspect outputs/latest
beachx animate outputs/latest --quantity potential --save-gif potential_history.gif
```

開発版を直接試す場合:

```bash
python -m pip install "git+https://github.com/Nkzono99/BEACH.git"
```

## ドキュメント

- 公開ドキュメント: [GitHub Pages](https://nkzono99.github.io/BEACH/)
- ドキュメント一覧: [`docs/index.md`](docs/index.md)
- 実行・開発・テスト: [`docs/fortran_workflow.md`](docs/fortran_workflow.md)
- 設定作成・検証: [`docs/config_workflow.md`](docs/config_workflow.md)
- `beach.toml` 仕様: [`docs/fortran_parameter_file.md`](docs/fortran_parameter_file.md)
- Python 後処理 API / CLI: [`docs/python_postprocess_api.md`](docs/python_postprocess_api.md)
- Fortran 実装仕様: [`SPEC.md`](SPEC.md)

## リポジトリ構成

- [`src/`](src/), [`app/`](app/): Fortran 本体
- [`beach/`](beach/): Python 後処理ライブラリ
- [`examples/`](examples/): 設定例・補助スクリプト
- [`tests/fortran/`](tests/fortran/), [`tests/python/`](tests/python/): テスト
- [`docs/`](docs/): 運用ドキュメント
- [`schemas/`](schemas/): JSON Schema

## 開発者向け最短確認

```bash
python -m pip install -e . --no-build-isolation
make check
```

KUDPC のログインノードでは `make test*` / `fpm test` を直接実行せず、計算ノードへ投入してください。詳細は [`docs/fortran_workflow.md`](docs/fortran_workflow.md) を参照してください。

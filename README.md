# BEACH（BEM + Accumulated CHarge）

[![GitHub Pages](https://img.shields.io/website?url=https%3A%2F%2Fnkzono99.github.io%2FBEACH%2F&up_message=GitHub%20Pages&down_message=Pages%20down)](https://nkzono99.github.io/BEACH/)
[![Fortran Docs](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-docs.yml/badge.svg?branch=main)](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-docs.yml)
[![Fortran Format](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-format.yml/badge.svg?branch=main)](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-format.yml)
[![PyPI version](https://img.shields.io/pypi/v/beach-bem)](https://pypi.org/project/beach-bem/)
[![Python versions](https://img.shields.io/pypi/pyversions/beach-bem)](https://pypi.org/project/beach-bem/)
[![License: Apache 2.0](https://img.shields.io/badge/license-Apache%202.0-green.svg)](LICENSE)

BEACH は、三角形メッシュで表した物体表面の帯電を計算するシミュレータです。
表面電荷が作る電場を境界要素法（BEM）で評価し、その場の中で荷電粒子を追跡します。
表面に到達した粒子の電荷を蓄積し、次のバッチの電場へ反映することで、帯電の時間発展を求めます。

<div align="center">
  <img src="https://nkzono99.github.io/BEACH/images/potential_history.gif" alt="帯電シミュレーションの電位変化" width="80%">
  <p><i>電子ビーム照射下での絶縁体メッシュ上の電位分布の時間発展</i></p>
  <sub>3D model: <a href="https://www.turbosquid.com/ja/3d-models/rubber-duck-pbr-game-ready-model-2001526">Rubber Duck PBR Game Ready</a> (TurboSquid)</sub>
</div>

## 特長

- 表面電荷の電場を direct 法または高速多重極法（FMM）で評価
- Boris 法による荷電粒子追跡と、表面での電荷蓄積
- x・y 方向の周期セルと、外部シースモデルとの接続
- 光電子放出、二次電子放出、浮遊導体
- OBJ ファイルと組み込み形状を組み合わせたメッシュ入力
- MPI・OpenMP による並列計算と、checkpoint からの再開
- 設定検査、出力の確認・可視化を行う Python ツール `beachx`

## インストールと最初の実行

Linux、Python 3.10 以上、`git`、`make`、Fortran コンパイラが必要です。

```bash
python -m pip install "git+https://github.com/Nkzono99/BEACH.git"

mkdir beach-tutorial
cd beach-tutorial
beachx config init beach.toml
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/tutorial
```

`beachx lint` は設定ファイルのキー・値・組合せを実行前に検査します。
各手順の意味と出力の読み方は [入門チュートリアル](https://nkzono99.github.io/BEACH/tutorial.html) で説明しています。

## 目的別の入口

| 目的 | 読む順序 |
| --- | --- |
| 初めて使う | [インストール](https://nkzono99.github.io/BEACH/installation.html) → [入門チュートリアル](https://nkzono99.github.io/BEACH/tutorial.html) → [最初の出力を読む](https://nkzono99.github.io/BEACH/output-guide.html) |
| 研究ケースを作る | [ケース設計](https://nkzono99.github.io/BEACH/configuration-recipes.html) → [設定](https://nkzono99.github.io/BEACH/configuration.html) → [実行・再開](https://nkzono99.github.io/BEACH/execution.html) → [妥当性確認](https://nkzono99.github.io/BEACH/validation-guide.html) |
| 物理・数値モデルを理解する | [BEACH の計算サイクル](https://nkzono99.github.io/BEACH/algorithms.html) → [表面はどう帯電するか](https://nkzono99.github.io/BEACH/surface-models.html) → [場 solver の選び方](https://nkzono99.github.io/BEACH/field-solvers.html) |
| ソースを変更する | [アーキテクチャ](https://nkzono99.github.io/BEACH/architecture.html) → [開発とテスト](https://nkzono99.github.io/BEACH/workflow.html) |
| キーや API を調べる | [入力パラメータ](https://nkzono99.github.io/BEACH/parameters.html) / [Python API](https://nkzono99.github.io/BEACH/python-postprocess-api.html) / [Fortran API](https://nkzono99.github.io/BEACH/fortran/) |

ドキュメント全体の案内は [BEACH Documentation](https://nkzono99.github.io/BEACH/) にあります。

## 提供するコマンドとライブラリ

| 名前 | 役割 |
| --- | --- |
| `beach` | シミュレーション本体（Fortran） |
| `beachx` | 設定検査、出力確認、可視化、負荷見積もり |
| `beach` Python パッケージ | 出力の読み込みと解析 |

コマンドの全オプションは `beach --help`、`beachx --help` および `beachx <command> --help` で確認できます。

## 開発

```bash
python -m pip install fpm
python -m pip install -e . --no-build-isolation
make check
```

テストの種類と実行方法は [開発とテスト](https://nkzono99.github.io/BEACH/workflow.html) を参照してください。

## ライセンス

[Apache License 2.0](LICENSE)

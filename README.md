# BEACH（BEM + Accumulated CHarge）

[![GitHub Pages](https://img.shields.io/website?url=https%3A%2F%2Fnkzono99.github.io%2FBEACH%2F&up_message=GitHub%20Pages&down_message=Pages%20down)](https://nkzono99.github.io/BEACH/)
[![Fortran Docs](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-docs.yml/badge.svg?branch=main)](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-docs.yml)
[![Fortran Format](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-format.yml/badge.svg?branch=main)](https://github.com/Nkzono99/BEACH/actions/workflows/fortran-format.yml)
[![PyPI version](https://img.shields.io/pypi/v/beach-bem)](https://pypi.org/project/beach-bem/)
[![Python versions](https://img.shields.io/pypi/pyversions/beach-bem)](https://pypi.org/project/beach-bem/)
[![License: Apache 2.0](https://img.shields.io/badge/license-Apache%202.0-green.svg)](LICENSE)

BEACH は、三角形表面上の電荷が作る電場を境界要素法（BEM）で評価し、その場の中を動く
荷電粒子を追跡する表面帯電シミュレータです。粒子が表面へ吸収されるたびに電荷を蓄積し、
次のバッチの電場へ反映します。

現行版の中心は絶縁体表面の帯電です。体積プラズマを PIC として解くコードではありません。
外部シースは、必要な場合だけ計算領域上端の matching plane を介した準定常の境界応答として接続します。
詳しい適用範囲と更新順は [BEACH の計算サイクル](https://nkzono99.github.io/BEACH/algorithms.html) を参照してください。

<div align="center">
  <img src="https://nkzono99.github.io/BEACH/images/potential_history.gif" alt="帯電シミュレーションの電位変化" width="80%">
  <p><i>電子ビーム照射下での絶縁体メッシュ上の電位分布の時間発展</i></p>
  <sub>3D model: <a href="https://www.turbosquid.com/ja/3d-models/rubber-duck-pbr-game-ready-model-2001526">Rubber Duck PBR Game Ready</a> (TurboSquid)</sub>
</div>

## 最初の実行

Linux、Python 3.10 以上、`git`、`make`、Fortran コンパイラが必要です。
HPC のログインノードではシミュレーションを直接実行せず、各サイトのジョブスケジューラから
計算ノードへ投入してください。以下の `beach` は、ローカル環境または計算ノードの割当内で実行します。

```bash
python -m pip install "git+https://github.com/Nkzono99/BEACH.git"

mkdir beach-tutorial
cd beach-tutorial
beachx config init beach.toml
beachx lint beach.toml
OMP_NUM_THREADS=1 beach beach.toml
beachx inspect outputs/tutorial
```

この公式入門ケースは、毎 batch 200 個のマクロ電子を 20 batch 入射し、蓄積した負電荷が
後続電子を反射し始めるところまでを確認します。設定の意味、期待出力、図の読み方は
[10 分チュートリアル](https://nkzono99.github.io/BEACH/tutorial.html) にまとめています。
このサイトは `main` branch の挙動を説明します。PyPI の安定版は同じ入門設定をまだ含まない場合があるため、
上の最初の実行では GitHub 版を使用します。

## 目的別の入口

| 目的 | 読む順序 |
| --- | --- |
| 初めて使う | [インストール](https://nkzono99.github.io/BEACH/installation.html) → [10 分チュートリアル](https://nkzono99.github.io/BEACH/tutorial.html) → [最初の出力を読む](https://nkzono99.github.io/BEACH/output-guide.html) |
| 研究ケースを作る | [ケース設計](https://nkzono99.github.io/BEACH/configuration-recipes.html) → [設定](https://nkzono99.github.io/BEACH/configuration.html) → [実行・再開](https://nkzono99.github.io/BEACH/execution.html) → [妥当性確認](https://nkzono99.github.io/BEACH/validation-guide.html) |
| 物理・数値モデルを理解する | [BEACH の計算サイクル](https://nkzono99.github.io/BEACH/algorithms.html) → [表面はどう帯電するか](https://nkzono99.github.io/BEACH/surface-models.html) → [場 solver の選び方](https://nkzono99.github.io/BEACH/field-solvers.html) |
| ソースを変更する | [アーキテクチャ](https://nkzono99.github.io/BEACH/architecture.html) → [開発とテスト](https://nkzono99.github.io/BEACH/workflow.html) |
| キーや API を調べる | [入力パラメータ](https://nkzono99.github.io/BEACH/parameters.html) / [Python API](https://nkzono99.github.io/BEACH/python-postprocess-api.html) / [Fortran API](https://nkzono99.github.io/BEACH/fortran/) |

ドキュメント全体の案内は [BEACH Documentation](https://nkzono99.github.io/BEACH/) にあります。

## 提供するコマンドとライブラリ

| 名前 | 役割 |
| --- | --- |
| `beach` | Fortran シミュレーション本体 |
| `beach-zhao-response` | matching-plane 用 Zhao 応答表の生成（現在の source build） |
| `beachx` | 設定検査、出力確認、可視化、負荷見積もり |
| `beach` Python package | 出力の読込と独自解析 |

コマンドの全オプションは `beach --help`、`beachx --help` および `beachx <command> --help` で確認できます。
`beach-zhao-response` が必要で公開 package に含まれない場合は、開発版をインストールしてください。

## 開発

```bash
python -m pip install fpm
python -m pip install -e . --no-build-isolation
make check
```

変更範囲に応じたテスト tier、MPI/OpenMP、KUDPC での実行方法は
[開発とテスト](https://nkzono99.github.io/BEACH/workflow.html) を参照してください。

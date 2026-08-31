title: インストール

Lang: [日本語](Installation.md) | [English](Installation.en.md)

# インストール

BEACH の Python package、`beachx`、Fortran 実行バイナリ `beach` は同時にインストールされます。通常の
`pip install` では、Fortran Package Manager (`fpm`) も隔離された build 環境に自動導入されます。

初回は次の前提を確認してください。

## 動作要件

| 項目 | 要件 |
| --- | --- |
| OS | 主にLinuxで動作確認済み。HPCでは各サイトのcompiler/MPI環境に従う |
| Python | 3.10以上 |
| source取得 | `git` |
| build tool | `make` |
| Fortran | `gfortran`、Intel Fortranなど |

```bash
python --version
git --version
make --version
gfortran --version
```

## このドキュメントと一致する版をインストール

このサイトは GitHub の `main` branch にある設定と出力を説明します。20 batch の公式チュートリアルを
同じ期待値で実行するには、現行 source をインストールします。

```bash
python -m pip install -U pip setuptools wheel
python -m pip install "git+https://github.com/Nkzono99/BEACH.git"
```

確認:

```bash
beach --version
beachx --help
```

`beach` または `beachx` が見つからない場合:

```bash
python -m site --user-base
export PATH="$HOME/.local/bin:$PATH"
```

表示されたuser base配下の`bin`を、shellの設定に追加してください。

## PyPI の安定版

```bash
python -m pip install beach-bem
```

PyPI 版は公開時点で固定された安定版です。公開版の `beachx config init` が作る入門設定は、このサイトの
`main` branch より古い場合があります。その場合は、このサイトの 20 batch 期待値と混ぜず、上の GitHub 版を
使ってください。

## checkout を編集しながら使う

checkoutを編集しながら使う場合:

```bash
python -m pip install -e . --no-build-isolation
```

`--no-build-isolation`は、`pyproject.toml`に定義されたbuild dependencyの自動導入を無効にします。
この方法や`make`を直接使う場合は、先に`fpm`をPATH上に用意してください。

```bash
python -m pip install fpm
fpm --version
```

matching-plane の応答表を作る補助コマンド `beach-zhao-response` は、現在の source build には含まれますが、
PyPI 配布には含まれない場合があります。必要な場合は GitHub 版をインストールし、
[matching-plane 準定常連成](MatchingPlaneCoupling.html#table-backend用の応答表を作る)へ進んでください。

## 更新と削除

```bash
python -m pip install -U beach-bem
python -m pip uninstall beach-bem
```

compilerや開発時のfpmが見つからない場合は、先にOS/HPCのpackage/moduleを設定してください。
失敗例は[トラブルシューティング](Troubleshooting.html)にまとめています。

次は[10 分チュートリアル](Tutorial.html)へ進みます。

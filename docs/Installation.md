title: インストール

Lang: [日本語](Installation.md) | [English](Installation.en.md)

# インストール

BEACH の Python package と Fortran 実行バイナリは同時にインストールされます。通常の
`pip install`では、Fortran Package Manager (`fpm`)も隔離build環境へ自動的に導入されます。

初回は次の前提を確認してください。

## 動作要件

| 項目 | 要件 |
| --- | --- |
| OS | Linux を主に確認。HPCでは各サイトのcompiler/MPI環境に従う |
| Python | 3.10以上 |
| build tool | `make` |
| Fortran | `gfortran`、Intel Fortranなど |

```bash
python --version
make --version
gfortran --version
```

## PyPIからインストール

```bash
python -m pip install -U pip setuptools wheel
python -m pip install beach-bem
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

表示されたuser baseの`bin`を恒久的にshell設定へ追加してください。

## 開発版

```bash
python -m pip install "git+https://github.com/Nkzono99/BEACH.git"
```

checkoutを編集する場合:

```bash
python -m pip install -e . --no-build-isolation
```

`--no-build-isolation`は`pyproject.toml`に定義されたbuild dependencyの自動導入を無効にします。
この方法や`make`を直接使う開発では、先に`fpm`をPATH上へ用意してください。

```bash
python -m pip install fpm
fpm --version
```

## 更新と削除

```bash
python -m pip install -U beach-bem
python -m pip uninstall beach-bem
```

compilerまたは開発時のfpmが見つからない場合は、先にOS/HPCのpackage/moduleを設定してください。
失敗例は[トラブルシューティング](Troubleshooting.html)にまとめています。

次は[10分チュートリアル](Tutorial.html)へ進みます。

title: beach.tomlを作成・検証する

Lang: [日本語](Configuration.md) | [English](Configuration.en.md)

# `beach.toml`を作成・検証する

設定ファイルは `run / domain / mesh / particles / fields / sheath / output` の7グループ形式だけを受け付けます。完全な入力例は `examples/` を参照してください。

この文書は、直接編集する `beach.toml` と `beachx config` の使い方をまとめたものです。
メッシュ、粒子源、境界条件を選んで物理的な構成を組み立てる手順は
[シミュレーションケースを設計する](ConfigurationRecipes.html)にまとめています。

- Fortran 実行系 `beach` は `beach.toml` を直接読みます。
- `beachx config init` は、公式チュートリアル（プラズマ中の金属板の浮遊電位）の設定を作ります。
- 全キーと、読み込み時に値を計算する座標・配置パラメータは[入力パラメータリファレンス](Parameters.html)にまとめています。

## 1. 基本フロー

```bash
mkdir beach-tutorial
cd beach-tutorial

beachx config init beach.toml
$EDITOR beach.toml
beachx lint beach.toml
beach beach.toml
```

実行前には `beachx lint` を通し、`status=ok` を確認してから `beach` を実行します。

`box_origin` / `box_size`、`source.region_mode`、`mesh.groups`なども通常のTOML keyとして直接書けます。
これらがどの座標・寸法を計算し、明示値を置き換えるかは[座標・配置の補助パラメータ](Parameters.html#座標配置の補助パラメータ)を確認してください。

## 2. コマンド

### 2.1 `init`

新しい `beach.toml` を作ります。既に存在する場合は失敗します。

```bash
beachx config init
beachx config init run.toml
beachx config init --force
```

生成内容は[`examples/tutorial_floating_plate.toml`](https://github.com/Nkzono99/BEACH/blob/main/examples/tutorial_floating_plate.toml)と
同一です。密度 1 cm⁻³、温度 1 eV の電子と陽子を上端面から流入させ、x/y 周期の浮遊導体板が
電子電流と陽子電流の釣り合う浮遊電位まで帯電する 5 s（10000 batch）を追う公式入門ケースです。
`[fields] boundary="periodic2"`、`fields.solver.method="fmm"`、境界 reservoir 流入を使い、
光電子と二次電子は含みません。

### 2.2 `lint`

TOML の読込、同梱の JSON Schema によるキー・型・列挙値・値域の検査、座標・配置の正規化を行い、
続けて `beach --check-config` で Fortran の preflight (組合せ制約と派生値) を実行します。
成功時は `checks=toml,schema,preflight` と `status=ok` を表示します。preflight が拒否した場合は Fortran の診断をそのまま表示します。

```bash
beachx lint beach.toml
beachx lint run.toml --schema schemas/beach.schema.json
beachx lint run.toml --beach /path/to/beach
```

`--schema` は同梱の BEACH スキーマに追加の制約を課します。指定したスキーマは正規化前後の設定に適用され、
通常の BEACH 検証を無効化したり、その制約を緩めたりすることはできません。
`beach` 実行ファイルは `--beach`、環境変数 `BEACH_CONFIG_CHECK_EXE`、`PATH` の順に探します。見つからない場合は失敗します。

### 2.3 `beach --check-config`

`lint` が内部で呼ぶ Fortran の preflight です。通常実行と同じ既定値・読取・正規化と、組合せ制約・派生値の検査を行います。
検査対象の path は必須です。通常利用では `beachx lint` だけを実行すれば十分です。

```bash
beach --check-config beach.toml
```

成功時は終了コード 0 と次の表示を返します。検出した入力エラーは非ゼロの終了コードで報告します。

```text
config=beach.toml
checks=toml,semantic
status=ok
```

この検査ではシミュレーションを開始せず、結果ファイルも作りません。
OBJ、応答表、checkpoint などの外部データの内容や、実行後の数値・物理的妥当性は検証しません。
外部ファイルの読み込みとモデル初期化は通常の `beach beach.toml` で続けて確認します。

### 2.4 `diff`

2つの設定を意味的に比較します。既定では座標・配置パラメータを実座標と実寸へ変換してから比較します。

```bash
beachx config diff left.toml right.toml
beachx config diff --raw left.toml right.toml
```

## 3. スキーマ

`beach.toml` の先頭に `#:schema` directive を置くと、VS Code の Even Better TOML / Taplo などで補完や型検証を使えます。

```toml
#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.schema.json
```

ローカル checkout の schema を使う場合:

```toml
#:schema ../schemas/beach.schema.json
```

BEACH の Fortran パーサは「最初のセクションより前の `key = value`」を受け付けないため、`"$schema" = "..."` ではなくコメント directive を使ってください。

## 4. よくある失敗

### 4.1 top-level keyを置く位置が正しくない

設定は[公開 TOML セクション](Parameters.html#toml-の階層とセクション一覧)の下へ書きます。
最初のセクションより前に通常キーを置いたり、未知の top-level セクションを追加したりすると validation または Fortran 読み込みで失敗します。

### 4.2 同じ座標を2通りで指定する

`box_origin` / `box_size`と`box_min` / `box_max`のように、同じ座標を2通りで書くと検証で失敗します。
ただし`size_mode="box_fraction"`とgroup scaleは、対応する寸法を計算値で置き換える仕様です。対象キーは
[入力パラメータリファレンス](Parameters.html#座標配置の補助パラメータ)に明記しています。

### 4.3 実行前に設定を検査する

実行前には `beachx lint beach.toml` で設定を検証し、成功後に `beach beach.toml` を実行します。
lint は同梱 schema の検査の後に `beach --check-config` を呼ぶので、組合せ制約の診断も lint で得られます。
数値は有限値、整数項目は整数として記述してください。Fortran が格納できる長さを超えた文字列は、
切り詰めずにエラーとして報告します。

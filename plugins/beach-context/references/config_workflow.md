title: beach.tomlを作成・検証する

Lang: [日本語](Configuration.md) | [English](Configuration.en.md)

# `beach.toml`を作成・検証する

この文書は、直接編集する `beach.toml` と `beachx config` の使い方をまとめたものです。
メッシュ、粒子源、境界条件を選んで物理的な構成を組み立てる手順は
[シミュレーションケースを設計する](ConfigurationRecipes.html)にまとめています。

- Fortran 実行系 `beach` は `beach.toml` を直接読みます。
- `beachx config init` は、多数粒子・20 batch の公式チュートリアル設定を作ります。
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

`box_origin` / `box_size`、`inject_region_mode`、`mesh.groups`なども通常のTOML keyとして直接書けます。
これらがどの座標・寸法を計算し、明示値を置き換えるかは[座標・配置の補助パラメータ](Parameters.html#座標配置の補助パラメータ)を確認してください。

## 2. コマンド

### 2.1 `init`

新しい `beach.toml` を作ります。既に存在する場合は失敗します。

```bash
beachx config init
beachx config init run.toml
beachx config init --force
```

生成内容は[`examples/tutorial_insulator.toml`](https://github.com/Nkzono99/BEACH/blob/main/examples/tutorial_insulator.toml)と
同一です。`volume_seed` から毎 batch 200 個のマクロ電子を絶縁体平面へ入射し、20 batch にわたる
電荷分布と後続粒子への feedback を確認する公式入門ケースです。理解しやすい
`field_solver="direct"` と `[field_boundary] mode="free"` を使い、周期境界、ion species、
`photo_raycast` は含みません。

### 2.2 `lint`

TOML、同梱の JSON Schema、座標・配置パラメータの組合せ、BEACH の既知制約をまとめて検証します。
成功時は `checks=toml,schema,semantic` と `status=ok` を表示します。

```bash
beachx lint beach.toml
beachx lint run.toml --schema schemas/beach.schema.json
```

`--schema` は同梱の BEACH スキーマに追加の制約を課します。指定したスキーマは正規化前後の設定に適用され、
通常の BEACH 検証を無効化したり、その制約を緩めたりすることはできません。

### 2.3 `validate`

`beach.toml` を読み、`lint` と同じ同梱スキーマ、座標・配置の正規化、意味的制約を検証します。
成功時は設定 path と `status=ok` を表示します。追加のスキーマ制約やエラー表示件数を指定する場合は `lint` を使います。

```bash
beachx config validate
beachx config validate run.toml
```

### 2.4 `beach --check-config`

開発・診断用に、Fortran 実行系の読み込み・正規化と、実行に必要な最小条件だけを確認します。検査対象の path は必須です。
通常利用では `beachx lint` の後に `beach` を実行すればよく、このコマンドの追加実行は必要ありません。

```bash
beach --check-config beach.toml
```

成功時は終了コード 0 と次の表示を返します。検出した入力エラーは非ゼロの終了コードで報告します。

```text
config=beach.toml
checks=toml,semantic
status=ok
```

この `status=ok` は読取・正規化・残されたモデル成立条件の検査が通ったことだけを示します。
列挙値、値域、無効な機能への指定などの詳細な診断は `beachx lint` が担当するため、実行前 lint の代替にはなりません。

この検査ではシミュレーションを開始せず、結果ファイルも作りません。
OBJ、応答表、checkpoint などの外部データの内容や、実行後の数値・物理的妥当性は検証しません。
外部ファイルの読み込みとモデル初期化は通常の `beach beach.toml` で続けて確認します。

### 2.5 `diff`

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
Fortran 側の設定読込を単独で調べる場合は、開発・診断用の `beach --check-config beach.toml` を使えます。
数値は有限値、整数項目は整数として記述してください。Fortran が格納できる長さを超えた文字列は、
切り詰めずにエラーとして報告します。

title: 設定を編集する

Lang: [日本語](Configuration.md) | [English](Configuration.en.md)

# 設定を編集する

この文書は、直接編集する `beach.toml` と `beachx config` の使い方をまとめたものです。

- Fortran 実行系 `beach` は `beach.toml` を直接読みます。
- `beachx config init` は、最小限の実行可能な `beach.toml` を作ります。
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
同一です。`volume_seed` から 1 個の電子を絶縁体平面へ向けて追跡する、`field_solver="fmm"`、
`field_bc_mode="periodic2"` の公式入門ケースです。x/y 軸を周期境界とし、
`field_periodic_image_layers=1`、`field_periodic_far_correction="none"` による $3\times3$ cell の有限画像和を使います。
無限周期補正、ion species、`photo_raycast` は含みません。

### 2.2 `lint`

TOML parse、JSON Schema、座標・配置パラメータの組合せ、BEACHの既知制約をまとめて検証します。

```bash
beachx lint beach.toml
beachx lint run.toml --schema schemas/beach.schema.json
```

### 2.3 `validate`

`beach.toml`を読み、座標・配置パラメータの組合せと最終設定の既知制約を検証します。
JSON Schema も含めて確認したい場合は `beachx lint` を使います。

```bash
beachx config validate
beachx config validate run.toml
```

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

実行時の設定は `sim`、`particles`、`mesh`、`output` の下へ書きます。
最初のセクションより前に通常キーを置いたり、未知の top-level セクションを追加したりすると validation または Fortran 読み込みで失敗します。

### 4.2 同じ座標を2通りで指定する

`box_origin` / `box_size`と`box_min` / `box_max`のように、同じ座標を2通りで書くと検証で失敗します。
ただし`size_mode="box_fraction"`とgroup scaleは、対応する寸法を計算値で置き換える仕様です。対象キーは
[入力パラメータリファレンス](Parameters.html#座標配置の補助パラメータ)に明記しています。

### 4.3 実行前に設定を検査する

`beachx lint beach.toml`でTOML parse、schema、座標・配置パラメータの組合せ、既知のBEACH制約をまとめて確認できます。

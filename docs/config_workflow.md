title: beachx config / 高水準記法ガイド

# `beachx config` / 高水準記法ガイド

この文書は、直接編集する `beach.toml` と `beachx config` の使い方をまとめたものです。

- Fortran 実行系 `beach` が読むのは `beach.toml` だけです。
- `beachx config init` は小さく実行可能な `beach.toml` を作ります。
- `beachx config render` は同じ `beach.toml` 内の高水準記法を最終キーへ展開します。
- 最終キーの仕様は [Fortran パラメータファイル仕様](fortran_parameter_file.html) を参照してください。

## 1. 基本フロー

```bash
mkdir run_periodic2
cd run_periodic2

beachx config init
$EDITOR beach.toml
beachx config validate
beachx config render
beach beach.toml
```

`render` は既定で入力ファイルを上書きします。内容を確認したい場合は `--stdout` を使います。

```bash
beachx config render --stdout
beachx config render beach.toml --output rendered/beach.toml
```

## 2. コマンド

### 2.1 `init`

新しい `beach.toml` を作ります。既に存在する場合は失敗します。

```bash
beachx config init
beachx config init run.toml
beachx config init --force
```

初期値は、周期 2 軸 FMM、volume seed の電子・イオン、平面メッシュ、標準出力設定を含む小さな確認用設定です。

### 2.2 `validate`

`beach.toml` を読み、高水準記法の整合性と最終設定の既知制約を検証します。

```bash
beachx config validate
beachx config validate run.toml
```

### 2.3 `render`

高水準記法を Fortran 実行系が読む最終キーに展開します。

```bash
beachx config render
beachx config render run.toml --output rendered/beach.toml
```

### 2.4 `diff`

2 つの設定を意味的に比較します。既定では高水準記法を展開してから比較します。

```bash
beachx config diff left.toml right.toml
beachx config diff --raw left.toml right.toml
```

## 3. 高水準記法

高水準記法は、研究者が意図しやすい座標指定を `beach.toml` に直接書くための補助です。`beachx config render` 後は通常の数値キーに変換されます。

### 3.1 計算箱

`sim.box_origin` と `sim.box_size` を使うと、`box_min` / `box_max` を直接計算できます。

```toml
[sim]
use_box = true
box_origin = [0.0, 0.0, -1.0]
box_size = [1.0, 1.0, 2.0]
```

render 後:

```toml
[sim]
use_box = true
box_min = [0.0, 0.0, -1.0]
box_max = [1.0, 1.0, 1.0]
```

### 3.2 注入領域

`reservoir_face` と `photo_raycast` では、面内の割合で注入領域を指定できます。

```toml
[[particles.species]]
source_mode = "reservoir_face"
inject_face = "z_high"
inject_region_mode = "face_fraction"
uv_low = [0.25, 0.25]
uv_high = [0.75, 0.75]
```

render 後は `pos_low` / `pos_high` に展開されます。

### 3.3 メッシュ配置

`mesh.templates` では、計算箱に対するアンカー指定が使えます。

```toml
[[mesh.templates]]
kind = "plane"
size_mode = "box_fraction"
size_frac = [1.0, 1.0]
placement_mode = "box_anchor"
anchor = "z_low_face_center"
offset_frac = [0.0, 0.0, 0.02]
nx = 20
ny = 20
```

render 後は `size_x` / `size_y` / `center` などへ展開されます。

### 3.4 グループ配置

`mesh.groups` は、複数 template に共通の原点やスケールを与えるための table です。

```toml
[mesh.groups.cavity_unit]
origin_mode = "box_anchor"
anchor = "box_center"
scale_from = "box_x"
scale = 0.5

[[mesh.templates]]
group = "cavity_unit"
kind = "sphere"
radius = 0.2
center = [0.0, 0.0, 0.0]
```

render 後、`mesh.groups`、`group`、`scale_from` などの高水準キーは消え、template ごとの実座標と実寸に変換されます。

## 4. スキーマ

`beach.toml` の先頭に `#:schema` directive を置くと、VS Code の Even Better TOML / Taplo などで補完や型検証を使えます。

```toml
#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.schema.json
```

ローカル checkout の schema を使う場合:

```toml
#:schema ../schemas/beach.schema.json
```

BEACH の Fortran パーサは「最初のセクションより前の `key = value`」を受け付けないため、`"$schema" = "..."` ではなくコメント directive を使ってください。

## 5. よくある失敗

### 5.1 reserved top-level key

`schema_version`、`use_presets`、`override`、`base_case` は旧設定レイヤのキーで、現行の direct `beach.toml` では使いません。`sim`、`particles`、`mesh`、`output` の下へ直接設定を書いてください。

### 5.2 render 後に高水準キーが消える

正常です。`box_origin`、`box_size`、`inject_region_mode`、`mesh.groups` などは、`beach` が直接読む最終キーではありません。

### 5.3 既存ファイルを上書きしたくない

`render` に `--stdout` または `--output` を指定してください。

```bash
beachx config render --stdout
beachx config render beach.toml --output beach.rendered.toml
```

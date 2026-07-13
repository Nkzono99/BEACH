title: beachx config / 高水準記法ガイド

Lang: [日本語](Configuration.md) | [English](Configuration.en.md)

# `beachx config` / 高水準記法ガイド

この文書は、直接編集する `beach.toml` と `beachx config` の使い方をまとめたものです。

- Fortran 実行系 `beach` は `beach.toml` を直接読み、高水準記法を読み込み時に解決します。
- `beachx config init` は、最小限の実行可能な `beach.toml` を作ります。
- 最終キーの仕様は [Fortran パラメータファイル仕様](Parameters.html) を参照してください。

## 1. 基本フロー

```bash
mkdir run_periodic2
cd run_periodic2

beachx config init
$EDITOR beach.toml
beachx lint beach.toml
beach beach.toml
```

高水準記法を使っても、別ファイルへ展開する必要はありません。Fortran parserは、`box_origin` / `box_size`、`inject_region_mode`、`mesh.groups`などを実行時キーへ正規化し、その後に検証します。

## 2. コマンド

### 2.1 `init`

新しい `beach.toml` を作ります。既に存在する場合は失敗します。

```bash
beachx config init
beachx config init run.toml
beachx config init --force
```

初期値として、周期2軸FMM、volume seedの電子・イオン、`photo_raycast`による電子放出、平面mesh、標準出力を含む小規模な確認用設定を生成します。

### 2.2 `lint`

TOML parse、JSON Schema、高水準記法、BEACH の既知制約をまとめて検証します。

```bash
beachx lint beach.toml
beachx lint run.toml --schema schemas/beach.schema.json
```

### 2.3 `validate`

`beach.toml` を読み、高水準記法の整合性と最終設定の既知制約を検証します。
JSON Schema も含めて確認したい場合は `beachx lint` を使います。

```bash
beachx config validate
beachx config validate run.toml
```

### 2.4 `diff`

2 つの設定を意味的に比較します。既定では高水準記法を正規化してから比較します。

```bash
beachx config diff left.toml right.toml
beachx config diff --raw left.toml right.toml
```

## 3. 高水準記法

高水準記法は、研究者が意図しやすい座標指定を `beach.toml` に直接書くための補助です。Fortran parser が読み込み時に通常の数値キーへ正規化します。

### 3.1 計算箱

`sim.box_origin` と `sim.box_size` を使うと、`box_min` / `box_max` を直接計算できます。

```toml
[sim]
use_box = true
box_origin = [0.0, 0.0, -1.0]
box_size = [1.0, 1.0, 2.0]
```

読み込み時の解決結果:

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

読み込み時に `pos_low` / `pos_high` に解決されます。

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

読み込み時に `size_x` / `size_y` / `center` などへ解決されます。

### 3.4 グループ配置

`mesh.groups`は、複数のtemplateに共通の原点やスケールを与えるtableです。

```toml
[mesh.groups.cavity_unit]
placement_mode = "box_anchor"
anchor = "box_center"
scale_from = "box_x"
scale_factor = 0.5

[[mesh.templates]]
group = "cavity_unit"
kind = "sphere"
radius = 0.2
center = [0.0, 0.0, 0.0]
```

読み込み時に、`mesh.groups`、`group`、`scale_from`などを使って、templateごとの実座標と実寸を決定します。

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

### 5.1 top-level keyを置く位置が正しくない

実行時の設定は `sim`、`particles`、`mesh`、`output` の下へ書きます。
最初のセクションより前に通常キーを置いたり、未知の top-level セクションを追加したりすると validation または Fortran 読み込みで失敗します。

### 5.2 高水準キーと実行時キーを混ぜる

`box_origin` / `box_size` と `box_min` / `box_max` のように、同じ値を表す高水準キーと実行時キーを同時に書くと検証で失敗します。どちらか一方に揃えてください。

### 5.3 実行前に設定を検査する

`beachx lint beach.toml` で TOML parse、schema、高水準記法の整合性、既知の BEACH 制約をまとめて確認できます。

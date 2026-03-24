title: beachx config / preset / 高水準記法ガイド

# `beachx config` / preset / 高水準記法ガイド

この文書は、`case.toml` から `beach.toml` を生成する **preset 合成レイヤ** の使い方をまとめたものです。

- Fortran 実行系 `beach` が読むのは最終的な `beach.toml` だけです。
- 日常的な編集対象は `case.toml` と preset TOML です。
- 高水準記法は `case.toml` / preset でだけ使え、`beachx config render` が数値キーへ展開します。

`beach.toml` 自体のキー仕様は [Fortran パラメータファイル仕様](fortran_parameter_file.html) を参照してください。

## 1. 役割分担

| ファイル / 置き場 | 役割 | 通常の編集対象 |
|---|---|---|
| `case.toml` | 1 つの計算ケースを定義する入口。`use_presets` と `override` を持つ | はい |
| preset (`sim/...`, `mesh/...` など) | 再利用可能な設定断片 | はい |
| `beach.toml` | `beach` が直接読む最終設定 | いいえ。通常は生成物 |
| `~/.config/beachx/cases/*.toml` | `beachx config save` で保存する case 雛形 | 必要に応じて |

基本の流れは次です。

1. `beachx config init` で `case.toml` を作る
2. `case.toml` と必要な preset を編集する
3. `beachx config validate` で合成前提を確認する
4. `beachx config render` で `beach.toml` を生成する
5. `beach beach.toml` で実行する

## 2. 最短の使い方

```bash
mkdir run_periodic2
cd run_periodic2

beachx config init
$EDITOR case.toml
beachx config validate
beachx config render
beach beach.toml
```

`beachx config init` の既定値は、次の built-in preset を組み合わせた `case.toml` を作ります。

- `sim/periodic2_fmm`
- `species/solarwind_electron`
- `species/solarwind_ion`
- `mesh/plane_basic`
- `output/standard`

サンプルは
[examples/periodic2_basic/case.toml](https://github.com/Nkzono99/BEACH/blob/main/examples/periodic2_basic/case.toml)
にあります。

## 3. `case.toml` の構造

`case.toml` は次の 4 要素だけを持つ設計です。

- `schema_version`
- `title`
- `use_presets`
- `override`

例:

```toml
#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.case.schema.json

schema_version = 1
title = "Periodic2 base case"
use_presets = [
  "sim/periodic2_fmm",
  "species/solarwind_electron",
  "species/solarwind_ion",
  "mesh/plane_basic",
  "output/standard",
]

[override.output]
dir = "outputs/periodic2_debug"

[override.sim]
batch_count = 10
```

`override` 以下は最終 `beach.toml` の構造に対応します。たとえば次のように書けます。

- `[override.sim]`
- `[[override.particles.species]]`
- `[override.mesh]`
- `[[override.mesh.templates]]`
- `[override.output]`

## 4. merge ルール

合成順序は常に次です。

1. `use_presets` を上から順に適用
2. `override` を最後に適用

merge の意味は次のとおりです。

- table は deep merge です
- `particles.species` と `mesh.templates` は array-of-tables として append します
- 通常の配列と scalar は後勝ちで置換します
- `mesh.groups.<name>` は preset 間で同名定義を許しません
- `override.mesh.groups.<name>` による最終上書きは許可します

つまり、species や template は preset を足し合わせる使い方に向いています。一方で `sim.box_min` や `output.dir` のような単一値は後から置き換わります。

## 5. `beachx config` コマンド

### 5.1 `init`

新しい `case.toml` を作ります。

```bash
beachx config init
beachx config init --title "Close-pack debug"
beachx config init --preset sim/periodic2_fmm --preset output/standard
```

保存済み雛形から始めることもできます。

```bash
beachx config init --from baseline_periodic2
```

### 5.2 `validate`

`case.toml` を読み、preset 解決・高水準記法展開・最終 `beach.toml` 制約まで含めて検証します。

```bash
beachx config validate
```

### 5.3 `render`

`case.toml` を最終 `beach.toml` へ変換します。

```bash
beachx config render
beachx config render --stdout
beachx config render custom_case.toml --output out/beach.toml
```

### 5.4 `diff`

2 つの case あるいは rendered config の差を比較します。

```bash
beachx config diff left.toml right.toml
beachx config diff --rendered case_a.toml case_b.toml
```

`--rendered` を付けると、preset 展開後の最終差分を見られます。

### 5.5 `save` / `list-saved`

現在の `case.toml` を再利用用に保存できます。

```bash
beachx config save baseline_periodic2
beachx config list-saved
beachx config init --from baseline_periodic2
```

保存先は `~/.config/beachx/cases/` です。

## 6. preset の使い方

### 6.1 preset の探索順

同名 preset は次の優先順位で解決されます。

1. `case.toml` または現在位置から見て最も近い `.beachx/presets/`
2. その親ディレクトリの `.beachx/presets/`
3. さらに親方向の `.beachx/presets/`
4. `~/.config/beachx/presets/`
5. package 同梱の built-in preset

たとえば `sim/periodic2_fmm` は次のいずれかから探されます。

- `./.beachx/presets/sim/periodic2_fmm.toml`
- `../.beachx/presets/sim/periodic2_fmm.toml`
- `../../.beachx/presets/sim/periodic2_fmm.toml`
- `~/.config/beachx/presets/sim/periodic2_fmm.toml`
- `beach/config/presets/sim/periodic2_fmm.toml`

複数の project-local preset が見つかった場合は、最短パスになるもの、つまり最も近い親ディレクトリ側の preset が使われます。

### 6.2 命名規則

- 名前は拡張子なしです
- `sim/...`, `species/...`, `mesh/...`, `output/...` を基本カテゴリにします
- `.` や `..` を含むパス、絶対パス、`.toml` 付き指定は使えません

### 6.3 preset に書く内容

preset は `beach.toml` とほぼ同じ構造の TOML 断片です。ただしトップレベルは必要な table だけを書きます。

例:

```toml
#:schema https://raw.githubusercontent.com/Nkzono99/BEACH/main/schemas/beach.preset.schema.json

[mesh]
mode = "template"

[[mesh.templates]]
kind = "plane"
size_x = 1.0
size_y = 1.0
nx = 20
ny = 20
center = [0.5, 0.5, 0.02]
```

`schema_version` / `use_presets` / `override` / `base_case` のような `case.toml` 用キーは preset では使えません。

### 6.4 `beachx preset` コマンド

新規作成:

```bash
beachx preset new sim/lab/periodic2_fast
beachx preset new sim/lab/periodic2_fast --from sim/periodic2_fmm
beachx preset new output/project/debug --local
```

- 既定の保存先は `~/.config/beachx/presets/`
- `--local` を付けると `./.beachx/presets/`
- built-in preset を直接編集するのではなく、`--from` で複製して使う運用を想定しています

一覧・表示・パス確認:

```bash
beachx preset list
beachx preset show sim/lab/periodic2_fast
beachx preset path sim/lab/periodic2_fast
beachx preset validate sim/lab/periodic2_fast
```

編集:

```bash
beachx preset edit sim/lab/periodic2_fast
```

`$VISUAL` または `$EDITOR` で user/project-local preset を開きます。built-in preset は read-only 扱いです。

抽出保存:

```bash
beachx preset save sim/lab/current_run --section sim
beachx preset save species/lab/ion --section species --index 2 --rendered
beachx preset save mesh/lab/plate --section mesh.templates --index 1
```

- `--section sim|mesh|output|species|mesh.templates`
- `species` と `mesh.templates` は `--index` が必須です
- `--rendered` を付けると入力を `beach.toml` として扱います

## 7. 高水準記法

### 7.1 基本方針

- 高水準記法を使えるのは `case.toml` と preset だけです
- `beach.toml` には最終的な数値キーだけを残します
- Fortran 側 `beach` は高水準記法を解釈しません

### 7.2 `[sim]`: box 基準指定

```toml
[override.sim]
box_origin = [0.0, 0.0, 0.0]
box_size = [1.0, 1.0, 10.0]
```

これは render 時に次へ展開されます。

```toml
[sim]
box_min = [0.0, 0.0, 0.0]
box_max = [1.0, 1.0, 10.0]
```

同じ fragment の中で `box_origin` と `box_min`、または `box_size` と `box_max` を同時に書くことはできません。

### 7.3 `[[particles.species]]`: face 上の相対注入領域

`reservoir_face` / `photo_raycast` では、面上の矩形領域を `uv` で相対指定できます。

```toml
[[override.particles.species]]
source_mode = "reservoir_face"
inject_face = "z_high"
inject_region_mode = "face_fraction"
uv_low = [0.25, 0.25]
uv_high = [0.75, 0.50]
```

これは inject face 上の `[0,1] x [0,1]` 座標として解釈され、render 時に `pos_low` / `pos_high` へ変換されます。

制約:

- `inject_region_mode` は `reservoir_face` / `photo_raycast` でのみ使用可能です
- `face_fraction` では `uv_low` / `uv_high` が必要です
- `uv_*` は各成分が `[0,1]` に入っている必要があります
- `uv_low <= uv_high` が必要です
- `face_fraction` と `pos_low` / `pos_high` の併用はできません

### 7.4 `[[mesh.templates]]`: anchor + offset

direct template では、`center` を直接書かずに box 基準配置できます。

```toml
[[override.mesh.templates]]
kind = "plane"
placement_mode = "box_anchor"
anchor = "z_low_face_center"
offset = [0.0, 0.0, 0.02]
size_mode = "box_fraction"
size_frac = [1.0, 0.5]
nx = 40
ny = 20
```

利用できる `anchor` は次です。

- `box_center`
- `x_low_face_center`
- `x_high_face_center`
- `y_low_face_center`
- `y_high_face_center`
- `z_low_face_center`
- `z_high_face_center`

`offset` は実長さ、`offset_frac` は box サイズ比率です。同時指定はできません。

### 7.5 `size_mode = "box_fraction"`

box 追従の寸法指定は、現在次の shape でサポートしています。

- `plane`, `plane_hole`, `plate_hole`: `size_frac = [fx, fy]`
- `box`: `size_frac = [fx, fy, fz]`
- `sphere`: `size_frac = f` から `radius = f * min(box_size)`
- `cylinder`: `size_frac = [fr, fh]`

`disk` と `annulus` は現状 `size_mode = "box_fraction"` 非対応です。

### 7.6 `mesh.groups.*`: 複数 template をまとめて配置

複数 template を 1 つの group として定義し、group 単位で平行移動と一様スケーリングをかけられます。

```toml
[override.mesh.groups.cavity_unit]
placement_mode = "box_anchor"
anchor = "box_center"
scale_from = "box_x"
scale_factor = 0.25

[[override.mesh.templates]]
group = "cavity_unit"
kind = "plane"
center_local = [0.0, 0.0, -1.0]
size_x = 4.0
size_y = 2.0
nx = 8
ny = 4
```

group で使えるキー:

- `placement_mode`
- `anchor`
- `offset`
- `offset_frac`
- `scale`
- `scale_from`
- `scale_factor`

`scale_from` の候補:

- `box_x`
- `box_y`
- `box_z`
- `box_min_xy`
- `box_max_xy`
- `box_min_xyz`
- `box_max_xyz`

grouped template では `center_local` を使い、render 時に `center` へ変換します。group に属する template では `center` や個別 `placement_mode` は使えません。

### 7.7 render 後に残るキー

最終 `beach.toml` に残るのは既存の低水準キーだけです。

- `box_min`, `box_max`
- `pos_low`, `pos_high`
- `center`
- `size_x`, `size_y`, `size`
- `radius`, `inner_radius`, `height`

高水準キーは render 後に消し込まれます。

## 8. schema と editor 補完

生成されるファイルには次の schema directive が付きます。

- `case.toml` → `beach.case.schema.json`
- preset TOML → `beach.preset.schema.json`
- `beach.toml` → `beach.schema.json`

VS Code の Even Better TOML / Taplo では、この `#:schema ...` により補完と検証が効きます。

## 9. よくある詰まりどころ

### 9.1 preset が見つからない

まず次を順に確認してください。

1. 名前に `.toml` を付けていないか
2. 現在位置や `case.toml` の親方向に `.beachx/presets/` があるか
3. `~/.config/beachx/presets/` に保存したつもりの preset があるか
4. built-in preset 名を打ち間違えていないか

### 9.2 built-in preset を編集したい

built-in preset は package 同梱物なので、直接編集せず複製して使います。

```bash
beachx preset new sim/lab/periodic2_fast --from sim/periodic2_fmm
beachx preset edit sim/lab/periodic2_fast
```

### 9.3 `beach.toml` に高水準キーを書いてはいけないのか

はい。`beach.toml` は最終正規形に保つ方針です。高水準記法は `case.toml` / preset で使い、render で数値へ解決してください。

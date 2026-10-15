# periodic2 照射セルへの責務の絞り込み

作成日: 2026-10-13。基準: `main` の 2.0.0 直後。状態: 採用、実装済み (2026-10-15)。

BEACH は汎用の表面帯電計算を想定して設計されたが、実際の利用と開発は 2 軸周期セル (`periodic2`) に
外部 1-D シースを接続した照射実験に集中している。本書は、サポートする構成をこの照射実験に絞り、
汎用性のために残っている入力形式・機能・検証経路を整理する判断を記録する。利用者はまだおらず、
公開契約の互換は保たない。

## 1. 判断

| 項目 | 判断 |
| --- | --- |
| 製品構成 | `periodic2` 照射セル 1 種類だけをサポートする (2 節) |
| 参照構成 | free-space の direct / FMM と入門ケースは検証用の参照として残す。開発対象にしない |
| 設定の検証 | 意味的制約の正本は Fortran の preflight とする。`beachx lint` は TOML・JSON Schema の検査に加えて `beach --check-config` を呼び、Python 側の意味検証は削除する |
| 入力形式 | grouped 7 グループ形式だけを残し、flat な旧形式とその変換層を削除する |
| mesh | `mesh.mode` と `[mesh.obj]` を廃止し、OBJ と template を同じ `[[mesh.parts]]` の要素として並べられるようにする |
| 削除する機能 | `velocity_grid` 分布、`reservoir_face`、`plane_source`、`potential_barrier` |
| conductor | 残す。等電位条件を場 snapshot の電位で立て、periodic2 でも使えるようにする |
| 前提の集約 | 周期軸、壁の面 (`z_high`)、下端モデルといった照射セルの前提を契約型に集め、runtime の各所で直接参照しない |

汎用化を将来やり直す場合は、kernel 層 (mesh、panel kernel、pusher、collision、direct / FMM) を汎用化し直すのではなく、
照射セルと並ぶ別の experiment 層を追加する。そのために kernel 層へ照射セルの前提を持ち込まない。

## 2. サポートする構成

### 2.1 製品構成: periodic2 照射セル

| 要素 | 固定する値 |
| --- | --- |
| 領域 | `[domain]` の box。`periodic_axes = ["x", "y"]` |
| 場 | `fields.boundary = "periodic2"`、FMM、`cached_kneq0`、split zero mode (`exclude_k0`)。direct は `panel_spectral_reference` の検証にだけ使う |
| 表面 | insulator accumulation、`mesh_id` ごとの浮遊 conductor |
| 粒子境界 | z は `open`。`reflect` / `redistributed_reflect` は closed photoelectron にだけ使う |
| 粒子源 | `volume_seed`、boundary reservoir 流入、`photo_raycast` |
| 流入の写像 | `source_vdf` / `infinity_barrier`、Zhao 零電流根による kinetic 流入 |
| 表面電流 closure | `explicit`、`fixed_current`、`neutral_return`、`zhao_stationary` と outflow refresh |
| 時間幅 | 固定 `batch_duration` と adaptive nonzero-mode step |

### 2.2 参照構成

| ケース | 目的 |
| --- | --- |
| free-space + direct + `volume_seed` + plane template (入門ケース) | 10 分チュートリアル。帯電の基本挙動 |
| free-space + direct / FMM | panel kernel、FMM、有限 image shell と無限周期解の比較の参照 |

参照構成はビルドとテストを維持する。新しい機能は追加しない。

### 2.3 範囲外

上記以外の組合せは lint と Fortran の両方で拒否する。`SPEC.md` 2.3 の未実装項目はそのまま範囲外とする。

## 3. 削除する機能と理由

| 機能 | 現状 | 理由 |
| --- | --- | --- |
| `velocity_distribution = "grid"` と CSV 速度分布 | examples 0、Fortran test 1、608 行 | 当面の照射実験で使わない。再度必要になったときに設計し直す |
| `reservoir_face` | deprecated 済み。examples 0 | boundary reservoir 流入で代替済み |
| `plane_source` | examples 0 | boundary reservoir 流入と重複 |
| `ordinary_open_model = "potential_barrier"` | examples 0 | Zhao kinetic 流入の障壁が製品構成の経路 |
| flat な旧 TOML 形式と `beachx config migrate` | 変換層が Fortran 約 700 行、Python 約 300 行 | 利用者がいない。grouped 形式を直接読む |

## 4. 設定検証の単一ソース化

- `schemas/beach.schema.json` はキー名、型、列挙値、値域の正本である。Fortran の読取実装とのキー集合一致は
  `make schema-check` が確認する。
- 組合せの制約と派生値の正本は Fortran の preflight である。Python に同じ制約を書かない。
- `beachx lint` は TOML の読込、JSON Schema の検査、`beach --check-config` の順に実行し、Fortran の診断を
  そのまま表示する。`beach` 実行ファイルは `--beach`、環境変数 `BEACH_CONFIG_CHECK_EXE`、`PATH` の順に探す。
- Python の `load_config_file` と `normalize_config_document` は schema と authoring の正規化だけを行う。
  負荷見積もり (`beachx workload`) が必要とする正規化は Python に残す。
- 契約テストは「lint の採否」と「Fortran の採否」が一致することではなく、lint が Fortran の採否を
  そのまま返すことを確認する。

## 5. mesh の構成

`[[mesh.parts]]` の各要素は必須の `kind` で形を選ぶ。`obj` は必須の `path` と任意の `scale`、`rotation`、`offset` を持つ。
template kind (`plane`、`plate_hole` / `plane_hole`、`disk`、`annulus`、`box`、`cylinder`、`sphere`) は各形状のキーを持つ。
全 kind で `enabled` (既定 true)、`surface_side`、`surface_model` を指定できる。
`mesh.groups`、`group` と placement / 相対寸法のキーは template kind 専用とし、template kind の `path`、`scale`、`rotation` も拒否する。
有効な要素をリスト順に 1 から連番の `mesh_id` へ割り当て、OBJ 内の全三角形はその要素の番号を持つ。
conductor はこの番号ごとに 1 つの浮遊導体 group を構成する。有効な要素が 0 件のときは template 0 件と同じく mesh 構築できない。
`mesh.mode`、`[mesh.obj]`、`[[mesh.templates]]` は互換経路なしで拒否する。

## 6. 実施

| 工程 | 内容 | commit |
| --- | --- | --- |
| 0〜2 | 本書と `SPEC.md` の範囲宣言、検証の単一ソース化、機能削除、conductor の periodic2 対応 | `2c5d6c1` |
| 3 | grouped 形式だけを受け付け、文書と Fortran の診断を公開キー名に統一 | `dae2f98`、`46e9394`、`c51133b` |
| 4 | `[[mesh.parts]]` | `141b0b1` |
| 5 | 照射セルの前提を `irradiation_cell_contract_type` に集約し、直書きの再発防止テストを追加 | `8d2bb33` |

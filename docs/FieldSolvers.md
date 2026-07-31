title: 場の評価

Lang: [日本語](FieldSolvers.md) | [English](FieldSolvers.en.md)

# 場の評価

BEACHは各batchの開始時点の要素電荷`q_elem`から電場を作り、そのbatchで追跡する粒子に同じ場を使います。
粒子が表面へ運んだ電荷はbatch末尾でまとめて反映されるため、場が変わるのは次のbatchです。

要素電荷は常に三角形上の一定面密度で表します。利用者が選ぶのは、その相互作用を計算するsolverと
field boundaryです。

| 選択 | 役割 | 値 |
| --- | --- | --- |
| 要素kernel | 1要素の電荷分布 | `triangle_p0`（固定） |
| solver | 多数のsourceをどう足すか | `direct` / `treecode` / `fmm` / `auto` |
| field boundary | 周期画像や遠方場をどう含めるか | `free` / `periodic2` |

## 問題規模からsolverを選ぶ

| solver | 主な用途 | 場境界 | 近似 |
| --- | --- | --- | --- |
| [Direct](DirectSolver.html) | 小規模計算、基準解、split reference | free、条件付きperiodic2 | 解析panel kernelを全要素について直接評価 |
| [Treecode](Treecode.html) | 中規模のfree-space計算 | free | 遠方nodeをmonopole、近傍leafを解析panel kernelで評価 |
| [FMM](FMM.html) | 大規模計算、多数の評価点 | free、periodic2 | 遠方相互作用を多重極・局所展開で近似 |
| `auto` | free境界で要素数に応じて選択 | free | DirectまたはFMM |

`auto`は`nelem < tree_min_nelem`ならDirect、それ以上ならFMMを使います。既定のしきい値は`256`です。
solver間の速度差は要素数だけでなく、
粒子数、step数、評価点の分布にも依存します。実際の計算条件に近い小規模ケースで測定してください。

## solverと場境界の互換表

この表をsolverと場境界の互換性に関する正本とします。

| solver | `free` | `periodic2` |
| --- | --- | --- |
| `direct` | 対応 | split referenceのみ。`periodic2.nonzero_mode_backend="panel_spectral_reference"`、`zero_mode_policy="exclude_k0"`、対応するlower-boundary modelが必須 |
| `treecode` | 対応 | 非対応 |
| `fmm` | 対応 | 対応。無限周期productionでは`cached_kneq0`を使用 |
| `auto` | 対応 | 非対応 |

`periodic2` では `[domain]` の box、`periodic_axes=["x", "y"]`、非周期 z 軸が必要です。
Direct split reference は小規模な基準解・検証用であり、通常の periodic2 production 経路は FMM です。

## domain topology と field boundary を分ける

```toml
[domain]
box_origin = [0.0, 0.0, 0.0]
box_size = [1.0, 1.0, 1.0]
periodic_axes = ["x", "y"]

[field_boundary]
mode = "periodic2"
```

`[domain]` は計算 cell と periodic topology、`[field_boundary]` は場の closure を定めます。
`mode="free"` は自由空間場、`mode="periodic2"` は domain の x/y 周期を使う 2 軸周期場です。現行
`periodic2` は x/y 周期・z 非周期だけを受理します。粒子種別の境界設定は非周期面の `open` / `reflect` /
`redistributed_reflect` を選べますが、周期軸を追加・削除・上書きできません。

## triangle P0で要素電荷を離散化する

BEACHは要素の総電荷を三角形上の一定面密度として扱います。`triangle_p0`は暗黙の唯一の要素kernelであり、
`[field]` tableでは選択しません。有限で非退化な三角形と、各要素で解決済みのvacuum sideが必要です。
Treecodeは近傍leafを解析panel核、遠方nodeをmonopoleで
電場・電位とも評価します。[Direct](DirectSolver.html#triangle-p0)、[Treecode](Treecode.html)、
[FMM](FMM.html#source-kernel)で、それぞれのtriangle P0評価を説明します。

旧`[field]` tableと`sim.softening`は削除済みです。残した入力はunknown table / keyとして停止し、
別の要素modelへ読み替えません。

## 長さを正規化して数値スケールを整える

`sim.field_normalization`は内部座標の代表長さ$L_0$を決めます。入力と出力の単位は変わりません。

| 値 | $L_0$ | 原点$\mathbf{x}_0$ |
| --- | --- | --- |
| `si` | 1 m | 0 |
| `length` | `field_length_scale` | 0 |
| `box` | `[domain]` の3辺の最大値 | `domain.box_min` |
| `mesh` | mesh bounding boxの最大幅 | mesh bounding boxの最小点 |

内部では

$$
\mathbf{x}'=\frac{\mathbf{x}-\mathbf{x}_0}{L_0}
$$

として評価し、電場には$k_c/L_0^2$、電位には$k_c/L_0$を掛けてSIへ戻します。
`box` には正の幅を持つ `[domain]` が必要です。空 mesh で `mesh` を選んだ場合だけ
`field_length_scale` を使います。

## periodic2ではsolverと境界成分を組み合わせる

`periodic2` は solver だけでは決まらず、有限画像または無限画像の nonzero mode と physical zero mode を
組み合わせます。production は FMM、小規模 split reference は Direct の panel spectral backend を使います。
成分構成は [periodic2 場計算](PeriodicElectrostatics.html)を参照してください。

## solver誤差とmesh離散化誤差を分けて測る

新しいmeshでは、まず小さな同一ケースでDirectとの差を調べます。その後、source meshを細分化して
離散化誤差を、solver設定を変えて近似誤差を分けて確認します。場の一点比較だけでなく、吸収位置、蓄積電荷、
保存量など最終的に使う観測量も比較してください。[計算結果の妥当性確認](ValidationGuide.html)に
収束確認の手順をまとめています。

## Code reference

- solverの初期化と自動選択: [`bem_field_solver_config.f90`](../src/physics/field_solver/bem_field_solver_config.f90)
- 電場・電位の評価: [`bem_field_solver_eval.f90`](../src/physics/field_solver/bem_field_solver_eval.f90)
- Treecode/FMMの電荷更新: [`bem_field_solver_tree.f90`](../src/physics/field_solver/bem_field_solver_tree.f90)
- 設定値: [設定パラメータ](Parameters.html#場ソルバ)

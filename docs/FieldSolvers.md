title: 場の評価

Lang: [日本語](FieldSolvers.md) | [English](FieldSolvers.en.md)

# 場の評価

BEACHは各batchの開始時点の要素電荷`q_elem`から電場を作り、そのbatchで追跡する粒子に同じ場を使います。
粒子が表面へ運んだ電荷はbatch末尾でまとめて反映されるため、場が変わるのは次のbatchです。

場の評価では、まず要素電荷をどう表すかを選び、次にその相互作用を計算する方式を選びます。

| 選択 | 役割 | 値 |
| --- | --- | --- |
| source kernel | 1要素の電荷分布 | `point` / `triangle_p0` |
| solver | 多数のsourceをどう足すか | `direct` / `treecode` / `fmm` / `auto` |
| field boundary | 周期画像や遠方場をどう含めるか | `free` / `periodic2` |

## 問題規模とkernelからsolverを選ぶ

| solver | 主な用途 | source kernel | 場境界 | 近似 |
| --- | --- | --- | --- | --- |
| [Direct](DirectSolver.html) | 小規模計算、基準解 | point、triangle P0 | free | 選んだkernelを全要素について直接評価 |
| [Treecode](Treecode.html) | 中規模のpoint source | point | free | 遠方nodeをmonopoleで近似 |
| [FMM](FMM.html) | 大規模計算、多数の評価点 | point、triangle P0 | free、periodic2 | 遠方相互作用を多重極・局所展開で近似 |
| `auto` | free境界で要素数に応じて選択 | point、triangle P0 | free | pointはDirect/Treecode、triangle P0はDirect/FMM |

`auto`は`nelem < tree_min_nelem`ならDirectを使います。それ以上では、point sourceにTreecode、
triangle P0にFMMを選びます。既定のしきい値は`256`です。solver間の速度差は要素数だけでなく、
粒子数、step数、評価点の分布にも依存します。実際の計算条件に近い小規模ケースで測定してください。

## source kernelで要素電荷の離散化を決める

### Point

`field.element_kernel="point"`は、要素の総電荷を三角形重心の点電荷として扱います。
`sim.softening`で重心近傍の特異性を緩和できます。既存ケースとの互換既定です。

### Triangle P0

`field.element_kernel="triangle_p0"`は、要素の総電荷を三角形上の一定面密度として扱います。
近傍場と自己電位を三角形の解析kernelで扱えるため、重心点電荷とは異なる離散化です。

triangle P0では`sim.softening=0`、有限で非退化な三角形、各要素のvacuum sideが必要です。
現行Phase 1はinsulator表面だけに対応し、Treecodeには対応しません。詳細は
[Direct](DirectSolver.html#triangle-p0)と[FMM](FMM.html#source-kernel)を参照してください。

## 長さを正規化して数値スケールを整える

`sim.field_normalization`は内部座標の代表長さ$L_0$を決めます。入力と出力の単位は変わりません。

| 値 | $L_0$ | 原点$\mathbf{x}_0$ |
| --- | --- | --- |
| `si` | 1 m | 0 |
| `length` | `field_length_scale` | 0 |
| `box` | boxの3辺の最大値 | `box_min` |
| `mesh` | mesh bounding boxの最大幅 | mesh bounding boxの最小点 |

内部では

$$
\mathbf{x}'=\frac{\mathbf{x}-\mathbf{x}_0}{L_0},
\qquad
\epsilon'=\frac{\epsilon}{L_0}
$$

として評価し、電場には$k_c/L_0^2$、電位には$k_c/L_0$を掛けてSIへ戻します。
`box`には`use_box=true`と正のbox幅が必要です。空meshで`mesh`を選んだ場合だけ
`field_length_scale`を使います。

## periodic2ではsolverと境界成分を組み合わせる

`periodic2`はsolverを1つ選ぶだけでは決まりません。有限画像和または無限画像和、非zero mode、zero mode、
外部シース、reservoir粒子の加減速、photoelectronのescape/returnを組み合わせる計算構成です。

legacyの`sim.field_bc_mode="periodic2"`経路はFMMを使います。小規模検証用のsplit referenceは、Directの
panel spectral backendを使います。それぞれの成分構成は、
[periodic2場計算](PeriodicElectrostatics.html)と[外部プラズマモデル](OuterPlasmaModels.html)で説明します。

## solver誤差とsource離散化誤差を分けて測る

新しいmeshやkernelでは、まず小さな同一ケースでDirectとの差を調べます。その後、source meshを細分化して
離散化誤差を、solver設定を変えて近似誤差を分けて確認します。場の一点比較だけでなく、吸収位置、蓄積電荷、
保存量など最終的に使う観測量も比較してください。詳しい手順は
[計算結果の妥当性確認](ValidationGuide.html)を参照してください。

## Code reference

- solverの初期化と自動選択: [`bem_field_solver_config.f90`](../src/physics/field_solver/bem_field_solver_config.f90)
- 電場・電位の評価: [`bem_field_solver_eval.f90`](../src/physics/field_solver/bem_field_solver_eval.f90)
- Treecode/FMMの電荷更新: [`bem_field_solver_tree.f90`](../src/physics/field_solver/bem_field_solver_tree.f90)
- 設定値: [設定パラメータ](Parameters.html#場ソルバ)

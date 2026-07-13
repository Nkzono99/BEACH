title: 表面帯電モデル

Lang: [日本語](SurfaceModels.md) | [English](SurfaceModels.en.md)

# 表面帯電モデル

このページでは、粒子が三角形表面へ衝突した後の処理と、要素電荷が次のbatchへどう引き継がれるかを
説明します。現行版の中心モデルは、吸収された粒子の電荷を絶縁体表面へ蓄積するモデルです。

## 共通の相互作用

粒子と表面の相互作用は既定で吸収です。粒子`p`が要素`i`へ最初に衝突すると、その粒子は追跡から
除かれ、重みを含む電荷をbatch局所差分へ加えます。

$$
\Delta q_i \mathrel{+}=q_p w_p
$$

同じbatchの粒子をすべて処理した後、`delta_q_elem`を`q_elem`へcommitします。したがって、同じbatchで
堆積した電荷は次batchから場へ寄与します。

## surface model

| `surface_model` | 現行の処理 | 注意 |
| --- | --- | --- |
| `insulator` | 命中要素に電荷を保持 | 現行版の主対象 |
| `conductor` | object総電荷を保ち、要素電位が等しくなるよう再配分 | `field_bc_mode="free"`のみ |
| `dielectric` | 電荷を保持し、`epsilon_r`をmetadataとして出力 | 誘電分極は未実装 |

### insulator

`q_elem`は要素総電荷[C]です。`field.element_kernel="triangle_p0"`を選んだ場合も保存量は総電荷のままで、
場評価時だけ面積で割って一定面密度へ変換します。表面内の横方向伝導や漏洩は現行モデルに含みません。

### floating conductor

`conductor`要素は`mesh_id`ごとに1つの浮遊導体として扱います。粒子電荷をcommitした後、objectの
総電荷を保存しながら各要素の電位を等しくする線形系を解き、要素電荷を再配分します。

これは接地電位を与える導体境界ではありません。またperiodic/outer modelとの組み合わせには制約があるため、
[設定パラメータ](Parameters.html)のsupport matrixを確認してください。

### dielectric metadata

`dielectric`と`epsilon_r`は、現行版では形状と材料情報を出力へ残すためのmetadataです。誘電率境界条件、
束密度の連続条件、分極電荷の自己整合計算は行いません。誘電分極を含む結果として解釈しないでください。

## 光電子放出の電荷ledger

`photo_raycast`で`deposit_opposite_charge_on_emit=true`の場合、放出粒子と反対符号の電荷を放出元要素へ
`photo_emission_dq`として加えます。後続の再衝突による電荷とは別に集計し、batch commit時に統合します。

放出後に戻る光電子を個別追跡するか、reduced closureで即時中和するかは粒子境界／外部モデルの選択です。

## 確認する出力

- `charges.csv`: commit後の要素電荷
- `charge_history.csv`: batchごとの電荷履歴
- `summary.txt`: 吸収数と電荷ledger
- `mesh_sources.csv`: `mesh_id`とsurface model

電荷収支とmesh依存性は[計算結果の妥当性確認](ValidationGuide.html)で確認してください。実装式を含む旧統合詳細は
[粒子追跡と表面電荷蓄積](ParticleChargeLoop.html#11-電荷堆積と表面モデル)に残しています。

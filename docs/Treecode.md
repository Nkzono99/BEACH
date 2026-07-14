title: Treecode

Lang: [日本語](Treecode.md) | [English](Treecode.en.md)

# Treecode

Treecodeはpoint sourceをoctreeへまとめ、遠方のsource群を1つのmonopoleとして評価します。近傍nodeと、
正負電荷が相殺するnodeはDirect和へ戻します。そのため、遠方nodeをまとめられる割合が精度と速度を左右します。
中規模のfree境界計算で、Directよりsource評価数を減らすためのsolverです。

| 特性 | 内容 |
| --- | --- |
| 対応kernel | pointのみ |
| 場境界 | freeのみ |
| 遠方近似 | nodeの総電荷と電荷中心によるmonopole |
| 近傍評価 | leaf内のsoftened point chargeをDirect和 |
| 電位 | 現行adapterではTreecode化せずDirect評価 |

```toml
[sim]
field_solver = "treecode"
field_bc_mode = "free"
softening = 1.0e-6
```

## 固定geometryからoctreeを一度作る

初期化時に各要素の三角形重心をsource位置として、次の手順でoctreeを作ります。

1. node内のsource重心を囲むaxis-aligned bounding boxを求める。
2. box中心で8つのoctantに分ける。
3. source数が`tree_leaf_max`以下になるまで再帰する。

重心が数値的に同じ位置へ集中している場合や、分割後も全sourceが同じoctantへ残る場合は、それ以上
分割できません。この場合は`tree_leaf_max`を超えていても、そのnodeをleafにします。mesh geometryが
固定されている間は、木のtopologyを再利用します。

## batchごとにnodeの電荷momentを更新する

表面電荷が変わるたびに木を作り直すのではなく、各nodeの電荷momentだけをleafからrootへ集計します。
node $n$について

$$
Q_n=\sum_{i\in n}q_i,
\qquad
A_n=\sum_{i\in n}|q_i|,
$$

$$
\mathbf{c}_{Q,n}=
\begin{cases}
Q_n^{-1}\sum_{i\in n}q_i\mathbf{c}_i, & |Q_n|>\mathrm{tiny},\\
\mathbf{c}_{n}, & \text{otherwise}
\end{cases}
$$

を更新します。$\mathbf{c}_n$はnodeの幾何中心です。このrefreshによって、batch中に固定する電場・電位を
構成します。batch末尾でcommitされた要素電荷は、次のrefreshから使われます。

## 距離に応じてnodeをまとめる

leafでは全sourceをDirectに評価します。内部nodeでは、node半径$R$、node中心から評価点までの距離$d$、
`tree_theta`を使い、

$$
R < \theta(d-R)
$$

を満たすnodeだけを遠方候補とします。評価点がnodeのbounding sphere内にある場合は必ず子へ降ります。
$\theta$を小さくすると多くのnodeを展開するため、計算は遅くなりますが精度は上がります。
$\theta$を大きくすると多くのnodeをまとめるため、速くなりますが評価は粗くなります。

遠方nodeは、$Q_n$を$\mathbf{c}_{Q,n}$に置いたsoftened point chargeとして加算します。
Treecodeは評価点ごとにsource treeを走査し、FMMのようなtarget側のlocal expansionは作りません。

## 相殺するnodeはDirect和まで展開する

正負電荷が相殺するnodeでは、$Q_n$が小さくなり、電荷中心$\mathbf{c}_{Q,n}$がnodeの外へ大きく動くことがあります。
この状態を単一monopoleで近似すると不安定なため、BEACHは幾何条件に加えて

$$
|A_n-|Q_n||
\le 64\,\epsilon_{\mathrm{mach}}\max(A_n,|Q_n|)
$$

を満たすnodeだけを受理します。このtoleranceは同符号電荷の集計に伴う丸め差だけを許します。
実質的にmixed-signのnodeは遠方でも子へ降り、最終的にleafのDirect和で評価されます。

この判定は相殺時の精度を保つ一方で、正負電荷が細かく混在する分布ではTreecodeの高速化効果を弱めます。
そのようなケースではDirectとの誤差だけでなく実測時間も確認し、必要ならFMMを比較してください。

## 精度と速度を決める設定

| key | 役割 | 制約 |
| --- | --- | --- |
| `tree_theta` | 遠方nodeを受理する幾何条件 | $0 < \theta \le 1$ |
| `tree_leaf_max` | leafに保持するsource数の目安 | 1以上 |
| `tree_min_nelem` | `field_solver="auto"`の切替しきい値 | 1以上 |
| `softening` | leaf Direct和とmonopoleのsoftening [m] | 0以上 |

`tree_theta`と`tree_leaf_max`を入力に明示しない場合、明示`treecode`を含めて要素数から次の値を選びます。

| `nelem` | `tree_theta` | `tree_leaf_max` |
| ---: | ---: | ---: |
| `< 1500` | 0.40 | 12 |
| `1500`–`9999` | 0.50 | 16 |
| `10000`–`49999` | 0.58 | 20 |
| `>= 50000` | 0.65 | 24 |

これは速度と精度の初期値であり、ケース固有の誤差保証ではありません。入力にどちらか一方だけを明示した場合は、
その値だけを上書きし、もう一方には上表の値を使います。

## 現行実装は電場評価を高速化する

現行のTreecode経路が高速化するのは粒子位置での電場です。`eval_potential`はpoint sourceのDirect和を使い、
全要素中心の`potential_history.csv`も$O(N^2)$のDirect評価です。電位履歴を有効にした性能を見積もるときは、
このコストを別に含めてください。

木の構築は初期化時、node momentのrefreshは電荷更新後、木走査は各評価点で行われます。理想的な分布では
1点あたりの評価がDirectより大幅に少なくなりますが、厳密な計算量は木の偏り、$\theta$、leaf size、
mixed-sign nodeの割合に依存します。

## Direct比較で近似誤差を測る

まず同じpoint kernel、softening、正規化、meshを使うDirect結果と比較します。

1. 場が強い領域、表面近傍、遠方、電荷相殺領域に評価点を置く。
2. `tree_theta`を小さくして結果がDirectへ近づくことを確認する。
3. `tree_leaf_max`を変え、結果と実行時間の感度を測る。
4. 一点の場だけでなく、粒子の命中要素とbatch後の`q_elem`も比較する。
5. release buildで本計算に近い粒子数とstep数を使って速度を測る。

回帰テストは同符号nodeがmonopole経路を使うこと、強く相殺するmixed-sign nodeがDirect精度を保つこと、
長さ正規化後もSIのDirect場と一致することを個別に確認しています。

## Code reference

- octree構築とmoment refresh: [`bem_field_solver_tree.f90`](../src/physics/field_solver/bem_field_solver_tree.f90)
- MAC、mixed-sign guard、木走査: [`bem_field_solver_eval.f90`](../src/physics/field_solver/bem_field_solver_eval.f90)
- parameterの選択: [`bem_field_solver_config.f90`](../src/physics/field_solver/bem_field_solver_config.f90)
- 精度・相殺・正規化テスト: [`test_dynamics_field_solver.f90`](../tests/fortran/test_dynamics_field_solver.f90)
- batch結果のDirect等価性: [`test_simulator.f90`](../tests/fortran/test_simulator.f90)

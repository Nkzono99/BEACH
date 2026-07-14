title: Boris粒子更新

Lang: [日本語](BorisPusher.md) | [English](BorisPusher.en.md)

# Boris粒子更新

BEACHは1つの粒子stepを、候補軌道の構築と、その軌道上で最初に起きる衝突・境界処理の二段階で進めます。
候補軌道は、予測中点での電場評価、Boris法による速度更新、台形則による位置更新から成ります。
公開状態とcheckpointが保持する$(\mathbf{x},\mathbf{v})$は、stepの前後で同じ時刻に揃っています。

| 状態 | 時刻 |
| --- | --- |
| 入力 | $\mathbf{x}^n,\mathbf{v}^n$ |
| 電場標本 | 予測位置$\mathbf{x}_\mathrm{mid}$ |
| 出力候補 | $\mathbf{x}^{n+1},\mathbf{v}^{n+1}$ |
| 磁場 | 一様な`sim.b0` |

## 予測中点で電場を評価する

step幅を$\Delta t$として、現在速度から予測中点を作ります。

$$
\mathbf{x}_\mathrm{mid}=
\mathbf{x}^n+\frac{1}{2}\mathbf{v}^n\Delta t.
$$

この位置でfield snapshotを1回評価します。

$$
\mathbf{E}_\mathrm{mid}=
\mathbf{E}_\mathrm{snapshot}(\mathbf{x}_\mathrm{mid}).
$$

snapshotは、要素電荷が作る場、`sim.e0`、選択したperiodic zero mode、outer profileを一度ずつ合成したものです。
同じbatchの粒子は、batch開始時の要素電荷から作った同じsnapshotを使います。[<sup>1</sup>](FieldSolvers.html)

`use_box=true`では、場の評価位置だけをsolverの有効領域へ写します。

| 軸 | $\mathbf{x}_\mathrm{mid}$の処理 |
| --- | --- |
| low/highの両方がperiodic | primary boxへ`modulo`でwrap |
| その他 | `box_min`から`box_max`へclamp |

候補軌道の端点は物理座標のまま保持します。これにより、後段で三角形との衝突とbox境界の通過を
時間順に比較できます。

## Boris法で速度を進める

$q$を粒子電荷、$m$を質量、$\mathbf{B}=\texttt{sim.b0}$として、電場による半step加速を行います。

$$
\mathbf{v}^-=
\mathbf{v}^n+\frac{q}{m}\mathbf{E}_\mathrm{mid}\frac{\Delta t}{2}.
$$

次に磁場回転を

$$
\mathbf{t}=\frac{q}{m}\mathbf{B}\frac{\Delta t}{2},
\qquad
\mathbf{s}=\frac{2\mathbf{t}}{1+\lVert\mathbf{t}\rVert^2},
$$

$$
\mathbf{v}'=\mathbf{v}^-+\mathbf{v}^-\times\mathbf{t},
\qquad
\mathbf{v}^+=\mathbf{v}^-+\mathbf{v}'\times\mathbf{s}
$$

で計算し、残りの電場半stepを加えます。

$$
\mathbf{v}^{n+1}=
\mathbf{v}^++\frac{q}{m}\mathbf{E}_\mathrm{mid}\frac{\Delta t}{2}.
$$

実装では$t^2=\mathbf{t}\cdot\mathbf{t}$を使う標準的なBoris回転を用います。

## 台形則で候補位置を作る

速度更新後、同時刻の入出力速度を使う台形則で候補位置を作ります。

$$
\mathbf{x}^{n+1}=
\mathbf{x}^n+\frac{1}{2}
\left(\mathbf{v}^n+\mathbf{v}^{n+1}\right)\Delta t.
$$

一様電場だけなら、一定加速度から求めた速度と変位に一致します。空間的に滑らかに変化する静電場では、
予測中点での場評価と組み合わせることで、位置と速度がともに二次収束することを回帰テストで確認しています。

## 衝突・境界でstepを分割する

上の式で得られる$(\mathbf{x}^{n+1},\mathbf{v}^{n+1})$は、step全体を進めた候補状態です。
候補軌道上に三角形との衝突やbox境界との交差がなければ、この状態を確定します。

三角形との衝突と通常のbox面交差は、$\mathbf{x}^n$と$\mathbf{x}^{n+1}$を結ぶ
**step内の軌道線分**上で発生順を決めます。三角形との衝突が先なら、その位置で粒子を吸収します。
反射面または周期境界面との交差が先なら、交差時刻の位置と速度を補間し、残り時間について新たな
予測中点場とBoris候補を計算します。

z-high outer interfaceでは、Boris更新の両端と整合する二次軌道をcoupling経路で構成し、交差時刻を
再評価します。[<sup>2</sup>](ParticleEvents.html)

## Boris更新が保つ性質

| 条件 | 性質 |
| --- | --- |
| 電場なし・一様磁場 | $\lVert\mathbf{v}\rVert$を丸め誤差まで保存 |
| 一定の外部E/B | 速度更新は$\Delta t$の符号反転で時間反転可能 |
| 滑らかな外部場 | 現行の同時刻更新は位置・速度とも二次精度 |
| 三角形メッシュへの衝突、開放境界からの離脱 | 衝突・境界位置で追跡を終える非可逆過程 |
| batch間で変わる自己無撞着場 | batch内ではsnapshotを固定し、commit後に次の場へ更新 |

## 軌道と帯電結果から`dt`を決める

`dt`は粒子の回転、場の空間変化、衝突対象の幾何を解像できる大きさにします。

- cyclotron角$|q|\lVert\mathbf{B}\rVert\Delta t/m$
- 電場が空間的に変化する長さを粒子が通過する時間
- 三角形やbox近傍を横切る時間
- 外部interface profileが急に変わる領域の通過時間

三角形との交差はstep内の軌道線分で判定します。1 step中に軌道が強く曲がる場合や細かい形状を通過する場合は、
時間刻みを小さくして衝突位置を収束させます。`dt`、`dt/2`、必要なら`dt/4`で、軌道、命中要素、
吸収・escape統計、batch後の要素電荷を比較してください。

`max_step`を変えずに`dt`を半分にすると、追跡できる物理時間も半分になります。同じ物理時間まで比較する場合は、
`max_step`も調整します。

## Code reference

- Boris速度回転と台形位置: [`bem_pusher.f90`](../src/physics/bem_pusher.f90)
- 予測中点場とremainder再積分: [`bem_particle_stepper.f90`](../src/runtime/simulator/bem_particle_stepper.f90)
- E/B、pure-B保存、時間反転、二次収束テスト: [`test_dynamics_basic.f90`](../tests/fortran/test_dynamics_basic.f90)
- charged-meshとfield sampleテスト: [`test_particle_stepper.f90`](../tests/fortran/test_particle_stepper.f90)

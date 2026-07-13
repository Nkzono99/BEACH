title: 粒子の衝突・境界イベント

Lang: [日本語](ParticleEvents.md) | [English](ParticleEvents.en.md)

# 粒子の衝突・境界イベント

BEACHは、現在の軌道線分で最初に起きる衝突・境界イベントまで粒子を進めます。三角形へ衝突すれば
その位置で吸収し、反射面または周期境界面へ達すれば、残り時間の軌道を作り直します。この処理を
繰り返して1つの粒子stepを完了します。

## 最初に起きるものだけを確定する

現在状態$(\mathbf{x}_0,\mathbf{v}_0)$と[Boris粒子更新](BorisPusher.html)で得た候補
$(\mathbf{x}_1,\mathbf{v}_1)$に対し、次の順で処理します。

1. 軌道線分$\mathbf{x}(t)=\mathbf{x}_0+t(\mathbf{x}_1-\mathbf{x}_0)$の最初のmesh hitを照会する。
2. 候補終点がbox内なら、mesh hitまたは候補終点を確定する。
3. box外なら、最初のbox面fractionとmesh hitの$t$を比較する。
4. meshが同時刻または先なら吸収、box面が先ならその面の作用を適用する。
5. reflectまたはperiodicで粒子が生き残れば、残り時間の候補を作り直す。

mesh hitとbox面の差が$64\epsilon_\mathrm{mach}\max(1,|t|)$以内ならmeshを先として扱います。
同一stepで複数の衝突・境界イベントがあり得ても、常に現在の軌道線分で最初に起きるものまで進んでから
残り時間を再積分します。

| 最初の衝突・境界イベント | 確定する処理 |
| --- | --- |
| mesh | hit位置と要素indexを返し、粒子を吸収 |
| open face | 境界との交差位置でescape、またはouter interfaceへ渡す |
| reflect face | 法線速度を反転し、box内側から残り時間を進める |
| periodic face | 反対側faceの内側へ移し、速度を変えず残り時間を進める |

## 軌道線分と三角形の交点を求める

### 候補三角形を絞る

mesh初期化時に、各三角形のaxis-aligned bounding box（AABB）を作ります。

| 要素数 | 候補探索 |
| ---: | --- |
| 64未満 | 全要素のAABBを線形に調べる |
| 64以上 | 一様gridと3D-DDAで軌道線分が通るcellだけを調べる |

一様gridのcell幅は、1 cellあたり8要素を目安に決めます。cell数は1軸あたり最大128です。各三角形は、
そのAABBと重なるcellへCSR形式で登録します。これらは現行実装の固定値で、入力parameterではありません。

queryでは、まず軌道線分とgrid AABBの交差区間を求め、3D-DDAで通過cellを順にたどります。同じ三角形が
複数cellに登録されていても、保持するのは最小の交差parameterだけです。DDAの反復上限は
`nx + ny + nz + 3`です。cell index、増分、parameterの異常によって有限回で進めない場合は、
`collision_query_grid_stalled`を返します。

### 最初の交点を確定する

各候補三角形にはMöller–Trumbore法を使います。三角形を

$$
\mathbf{r}(u,v)=\mathbf{v}_0
+u(\mathbf{v}_1-\mathbf{v}_0)
+v(\mathbf{v}_2-\mathbf{v}_0)
$$

と書き、次をすべて満たす交点だけを採用します。

$$
0\le u\le1,
\qquad
0\le v,
\qquad
u+v\le1,
\qquad
0\le t\le1.
$$

determinantの大きさが、軌道線分の長さと二辺長の積に対して
$64\epsilon_\mathrm{mach}$以下なら、退化またはほぼ平行として除外します。determinantの符号では除外しないため、
衝突は三角形の両側から検出します。triangle windingは衝突の表裏を決めず、場のvacuum-side traceなど別の
用途に使われます。

複数三角形に命中した場合は最小$t$を採用します。交点は
$\mathbf{h}=\mathbf{x}_0+t(\mathbf{x}_1-\mathbf{x}_0)$です。

## 周期画像上の衝突をprimary cellへ対応づける

periodic2でも、meshが保持するのはprimary cellのbase要素だけです。軌道線分のAABBとcanonical mesh AABBを
重ねるために必要なimage shiftを、2つの周期軸について列挙します。周期長を$L$、対象軸の軌道線分の範囲を
$[p_\min,p_\max]$、mesh範囲を$[m_\min,m_\max]$とすると、

$$
n_\min=\left\lceil
\frac{p_\min-m_\max-\mathrm{tol}}{L}
\right\rceil,
\qquad
n_\max=\left\lfloor
\frac{p_\max-m_\min+\mathrm{tol}}{L}
\right\rfloor.
$$

各imageについて軌道線分を$-nL$だけbase mesh側へ写し、通常の交差判定を行います。hitは次を保持します。

| 値 | 意味 |
| --- | --- |
| `hit%pos` | 実際に命中したperiodic image上の物理座標 |
| `hit%pos_wrapped` | primary cellへwrapした座標 |
| `hit%image_shift` | 2周期軸のimage index |
| `hit%elem_idx` | base meshの要素index |

候補$t$が`1e-12`の相対tolerance内で一致する場合は、要素index、第1 image index、第2 image indexの順に
比較して一意に選びます。1軸のimage数または2軸の直積が4096を超えるqueryは、
`collision_query_image_limit`で停止します。

ここで列挙するimage範囲は粒子の軌道線分がmeshへ当たり得る範囲です。電場の`field_periodic_image_layers`とは
目的も決め方も異なります。場の画像和は[periodic2場計算](PeriodicElectrostatics.html)を参照してください。

## box面ごとの作用を適用する

`use_box=true`では、各軸のlow/high faceに`open`、`reflect`、`periodic`を設定します。候補軌道線分が複数faceへ
同時に達するcorner/edgeでは、最小fractionからmachine-epsilon tolerance内のfaceを1つのmaskへまとめます。

### open、reflect、periodic

`open_boundary_model="escape"`では、maskにopen faceが1つでも含まれていれば粒子を消滅させ、
`escaped_boundary`へ数えます。reflectだけのfaceでは、該当する速度成分を反転します。periodic faceでは
粒子を反対側へ移し、速度は変えません。境界処理後も生存する粒子は、`nearest`を使ってfaceから
1 floating-point値だけboxの内側へ置きます。

cornerでreflectとperiodicが組み合わさっても、face maskへまとめてから各軸へ作用するため、軸の走査順序に
依存しない結果になります。

### legacy potential barrier

`open_boundary_model="potential_barrier"`は、単一open faceから外向きに出る粒子について法線運動エネルギー

$$
K_n=\frac{1}{2}m v_\mathrm{out}^2
$$

と電位障壁

$$
\Delta U=q\left(\phi_\infty-\phi_\mathrm{boundary}\right)
$$

を比較します。$v_\mathrm{out}>0$、$\Delta U>0$かつ$K_n<\Delta U$ならreflectし、それ以外はescapeします。
これはlegacyの単一face reduced modelです。複数open faceが同時に関係するcornerは一般化せず、
`particle_step_unsupported_barrier_corner`で停止します。

reservoir粒子の加減速、outer sheath、photoelectron returnは、衝突位置を決めた後に適用する物理モデルです。
詳細は[粒子源](ParticleSourcesBoundaries.html)と[外部プラズマモデル](OuterPlasmaModels.html)を参照してください。

## 境界通過後の残り時間を進める

box境界イベント時の位置と速度は、候補軌道の入出力状態を交差fractionで補間します。そのうえで、
faceに垂直な座標をbox面の値へ正確に揃えます。粒子が生存する境界作用を適用した後、

$$
\Delta t_\mathrm{remain}=(1-t_\mathrm{event})\Delta t_\mathrm{segment}
$$

について、新しい予測中点場とBoris候補を計算します。新たな軌道線分でも、meshとの衝突とbox境界通過を
比較し、最初に起きるものを調べます。場も再評価するため、反射前のfield sampleは使い回しません。

1 particle step内で処理できるbox境界イベントは最大8回です。8回までは正しく処理し、9回目が必要なら
`particle_step_multiple_box_events`を返して未完成stateをcommitしません。これは非常に大きい`dt`、狭いbox、
または高速粒子を検出する安全上限でもあります。

## z-highから外部モデルへ粒子を渡す

particle transferが有効なとき、z-high open faceに達した粒子はその場でescapeさせず、interface crossingとして
callerへ返します。このfaceでは、Boris更新の入出力位置と速度に整合する二次軌道を使って交差時刻を再評価します。
payloadには、face、step全体に対するfraction、位置、速度、残り時間が入ります。

outer modelがlocal returnを返した場合は、戻り位置・速度から残り時間を再び通常のparticle stepで進めます。
infinityへescapeした場合は粒子を消滅させます。outer側の加減速やreturn条件は
[外部プラズマモデル](OuterPlasmaModels.html)を参照してください。

## 判定を完了できなければ停止する

collision queryは、必要な候補をすべて調べた場合だけ`ok`です。

| status | code | 意味 |
| --- | ---: | --- |
| `collision_query_ok` | 0 | query完了 |
| `collision_query_image_limit` | 1 | periodic image列挙が4096上限を超過 |
| `collision_query_index_range` | 2 | image boundが非有限または整数範囲外 |
| `collision_query_invalid_segment` | 3 | 軌道線分の端点が非有限 |
| `collision_query_grid_stalled` | 4 | grid geometry不正またはDDAが進行しない |
| `particle_step_invalid_boundary` | 1001 | particle、box、衝突・境界geometryが不正 |
| `particle_step_multiple_box_events` | 1002 | 1 stepで9回目のbox境界イベントが必要 |
| `particle_step_unsupported_barrier_corner` | 1003 | legacy barrierで複数open faceが同時に発生 |

これらの異常を「命中なし」とみなすと、粒子が表面を通り抜けて電荷ledgerを壊します。そのため、通常の
追跡はfail closedです。OpenMP内では、particle/step番号が最小の失敗情報を選びます。MPI実行では失敗rankと
位置・速度を全rankで共有し、同じbatch/rank/particle/step/statusを報告して停止します。photo raycastも、
species/ray/bounceについて同じ方針を使います。

## 衝突位置と帯電結果を収束させる

1. `dt`を半分にして軌道線分が衝突する要素と位置が安定するか確認する。
2. 薄い面、edge/corner近傍、三角形とほぼ平行な軌道線分を含む小ケースを作る。
3. reflect後とperiodic wrap後の残り時間内にmeshへ当たる軌道を確認する。
4. periodic seamで`pos`と`pos_wrapped`が意図した要素へ対応するか確認する。
5. `multiple_box_events`が出る場合、上限を回避するために`dt`を小さくする。
6. mesh refinementで最初のhitと最終的な表面電荷分布が収束するか確認する。

## Code reference

- box境界通過とmesh衝突の順序、残り時間の再積分: [`bem_particle_stepper.f90`](../src/runtime/simulator/bem_particle_stepper.f90)
- 三角形、grid DDA、periodic image: [`bem_collision.f90`](../src/physics/bem_collision.f90)
- box face検出と作用: [`bem_boundary.f90`](../src/physics/bem_boundary.f90)
- collision grid構築: [`bem_mesh.f90`](../src/mesh/bem_mesh.f90)
- 衝突・境界イベントの順序と8回上限のテスト: [`test_particle_stepper.f90`](../tests/fortran/test_particle_stepper.f90)
- 境界単体テスト: [`test_boundary.f90`](../tests/fortran/test_boundary.f90)
- 衝突・periodic・fail-closedテスト: [`test_dynamics_basic.f90`](../tests/fortran/test_dynamics_basic.f90)

title: 粒子の衝突・境界イベント

Lang: [日本語](ParticleEvents.md) | [English](ParticleEvents.en.md)

# 粒子の衝突・境界イベント

BEACH は、現在の軌道線分で最初に起きる mesh 衝突または box 境界イベントまで粒子を進めます。
このページは、交点探索、イベントの順序、境界通過後の残り時間、fail-closed status の参照です。

## 最初に起きるものだけを確定する

現在状態$(\mathbf{x}_0,\mathbf{v}_0)$と[Boris粒子更新](BorisPusher.html)で得た候補
$(\mathbf{x}_1,\mathbf{v}_1)$に対し、次の順で処理します。

1. 軌道線分$\mathbf{x}(t)=\mathbf{x}_0+t(\mathbf{x}_1-\mathbf{x}_0)$の最初のmesh hitを照会する。
2. 候補終点がbox内なら、mesh hitまたは候補終点を確定する。
3. box外なら、最初のbox面fractionとmesh hitの$t$を比較する。
4. meshが同時刻または先なら吸収、box面が先ならその面の作用を適用する。
5. `reflect`、`redistributed_reflect`、periodicのいずれかで粒子が生き残れば、残り時間の候補を作り直す。

mesh hitとbox面の差が$64\epsilon_\mathrm{mach}\max(1,|t|)$以内ならmeshを先として扱います。
同一stepで複数の衝突・境界イベントがあり得ても、常に現在の軌道線分で最初に起きるものまで進んでから
残り時間を再積分します。

| 最初の衝突・境界イベント | 確定する処理 |
| --- | --- |
| mesh | hit位置と要素indexを返し、粒子を吸収 |
| open face | 境界との交差位置で `particle_boundary.ordinary_open_model` を適用 |
| reflect face | 法線速度を反転し、box内側から残り時間を進める |
| `redistributed_reflect` face | 法線速度を反転し、面内位置を一様再配置して残り時間を進める |
| periodic face | 反対側faceの内側へ移し、速度を変えず残り時間を進める |

## 軌道線分と三角形の交点を求める

### 候補三角形を絞る

mesh初期化時に、各三角形のaxis-aligned bounding box（AABB）を作ります。

| 要素数 | 候補探索 |
| ---: | --- |
| 64未満 | 全要素のAABBを線形に調べる |
| 64以上 | 一様gridと3D-DDAで軌道線分が通るcellだけを調べる |

一様 grid は 1 cell あたり 8 要素を目安とし、各軸を最大 128 cell にします。各三角形は AABB と重なる
cell へ CSR 形式で登録します。これらは入力 parameter ではなく、現行実装の固定値です。

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
$[p_{\min},p_{\max}]$、mesh範囲を$[m_{\min},m_{\max}]$とすると、

$$
n_{\min}=\left\lceil
\frac{p_{\min}-m_{\max}-\mathrm{tol}}{L}
\right\rceil,
\qquad
n_{\max}=\left\lfloor
\frac{p_{\max}-m_{\min}+\mathrm{tol}}{L}
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
目的も決め方も異なります。[periodic2場計算](PeriodicElectrostatics.html)では、場の画像和を独立に定義します。

## box面ごとの作用を適用する

`domain.periodic_axes` に含む軸の両面は periodic です。非周期面は `[particle_boundary]` の
`x_low`、`x_high`、`y_low`、`y_high`、`z_low`、`z_high` で `open`、`reflect`、
`redistributed_reflect` のいずれかを選びます。
粒子境界へ `periodic` は指定できません。候補軌道線分が複数 face へ同時に達する corner/edge では、
最小 fraction から machine-epsilon tolerance 内の face を 1 つの mask へまとめます。

### open、reflect、redistributed_reflect、periodic

`particle_boundary.ordinary_open_model="escape"` では、mask に open face が 1 つでもあれば粒子を消滅させ、
`escaped_boundary` へ数えます。reflect系の face は該当する速度成分を反転し、periodic face は速度を変えずに
反対側へ移します。生存粒子のevent軸座標はbox scaleに応じた内側guardへ置きます。

通常の `reflect` はevent位置の接線成分を変更しません。`redistributed_reflect` は速度には同じ反射を適用し、
位置だけを再配置します。単一面eventでは、面内2軸をそれぞれbox spanの両端guardを除く範囲から一様に再標本化します。
edge / cornerの同時eventでmaskに`redistributed_reflect`が含まれる場合は、maskに含まれない軸だけを
一様再標本化します。event軸は再標本化せず、対応する面の内側guardへ置きます。このためcornerでは
再配置する軸はなく、edgeでは残る1軸だけを再配置します。

面内標本は共有乱数streamではなく、`sim.rng_seed`とbatch、rank、particle、step、event、axisのcounterから
生成します。同じMPI layoutとseedではOpenMP scheduleに依存せず再現できますが、rank分割を変えた同一軌道の
一致は保証しません。

`[particles.species.boundary]` は同じ 6 面を species ごとに `inherit`、`open`、`reflect`、
`redistributed_reflect` で上書きします。
`inherit` は global の `[particle_boundary]` を使います。`domain.periodic_axes` の面は topology なので、
global 設定でも species 設定でも上書きできません。closed PE の組合せは
[光電子の放出とライフサイクル](PhotoelectronEmission.html)にあります。

cornerでreflect系とperiodicが組み合わさっても、face maskへまとめてから各軸へ作用するため、軸の走査順序に
依存しない結果になります。

### potential barrier

`particle_boundary.ordinary_open_model="potential_barrier"` は、単一 open face から外向きに出る粒子について法線運動エネルギー

$$
K_n=\frac{1}{2}m v_\mathrm{out}^2
$$

と電位障壁

$$
\Delta U=q\left(\phi_\infty-\phi_\mathrm{boundary}\right)
$$

を比較します。$v_\mathrm{out}>0$、$\Delta U>0$かつ$K_n<\Delta U$ならreflectし、それ以外はescapeします。
これは単一faceのreduced modelです。複数open faceが同時に関係するcornerは一般化せず、
`particle_step_ambiguous_open_corner`で停止します。

reservoir 粒子の流入補正と closed PE はこのページの対象外です。[reservoir 注入](ReservoirInjection.html)と
[粒子の escape と局所 return](ParticleEscapeReturn.html)を参照してください。

## 境界通過後の残り時間を進める

box境界イベント位置は候補chordの交差fractionで求め、event maskに含まれる各軸を対応するbox面へ正確に
揃えます。event速度の向きはchord接線に合わせ、速さは予測中点電場がevent位置までに行う離散workから求めます。

$$
\lVert\mathbf{v}_\mathrm{event}\rVert^2=\lVert\mathbf{v}_0\rVert^2+
2(q/m)\mathbf{E}_\mathrm{mid}\cdot(\mathbf{x}_\mathrm{event}-\mathbf{x}_0)
$$

これにより、強い磁場で部分stepのBoris速度だけが内向きへ反転していても、外向きchord crossingへ内向き速度を
適用しません。純磁場ではevent速度のノルムを保存します。非正・非有限の速さ、または零長chordは
未解決eventとしてfail closedです。chord交差fractionは残り時間のfractionにも使いますが、加速または
磁場旋回中の真の交差時刻ではありません。この時間離散化近似が境界作用やpotential-barrier判定へ与える影響は、
`dt`、`dt/2`、必要なら`dt/4`で収束を確認します。粒子が生存する境界作用を適用した後、
reflect系/periodic後のevent軸座標はbox座標とspanに応じたguard幅だけ面内へ置きます。原点側でsubnormalに
なる1 ULP offsetを避け、残り軌道による次のevent fractionが0へunderflowする境界chatterを防ぎます。

$$
\Delta t_\mathrm{remain}=(1-t_\mathrm{event})\Delta t_\mathrm{segment}
$$

について、新しい予測中点場とBoris候補を計算します。新たな軌道線分でも、meshとの衝突とbox境界通過を
比較し、最初に起きるものを調べます。場も再評価するため、反射前のfield sampleは使い回しません。

1 回の local continuation では box 境界イベントを最大 8 回処理します。9 回目が必要なら
`particle_step_multiple_box_events` を返し、未完成 state を commit しません。これは大きすぎる `dt`、
狭い box、または高速粒子を検出する安全上限です。guard 幅を変えてもこの上限は変わりません。

`multiple_box_events_retry_backend="upper_panel_fourier"` を指定した `cached_kneq0` 構成では、この status が
出た 1 step だけを元の位置・速度から再試行します。再試行の非零モード場は triangle P0 電荷を
`periodic2.reference_mode_layers` まで Fourier 展開し、全 mesh 頂点の最大 z より上で成立する指数減衰形へ
因子化します。各評価では同じ periodic zero mode と `sim.e0` を 1 回だけ加えます。再試行中の場評価点が
1 点でもこの上部真空域に入らない場合は失敗とします。potential-barrier event の境界電位も同じ展開と
potential gauge で評価し、その評価点にも同じ成立域を要求します。成立域外の場合、または再試行後も
event を完了できない場合は、元の
`multiple_box_events` を保持して policy へ渡します。これは外部プラズマやシースを追加するモデルではなく、
通常粒子の FMM backend も変更しません。
幾何応答は snapshot 初期化時に 1 回だけ構築し、現在の面電荷との積は再試行の電場・電位評価時だけ計算します。
モード数と幾何応答メモリは `reference_mode_layers` の 2 乗で増えるため、この値は seam 誤差と再試行結果の
収束を確認して選びます。

既定の`multiple_box_events_policy="abort"`はこの時点でRUNをfail closedにします。有限画像和の定性的な
感度確認など、明示的に`"soft_discard"`を選んだ場合だけ、該当macro particleを消滅させます。標準エラーには
batchごとの全rank合計件数と絶対macro chargeだけを記録します。同じ集計値は`summary.txt`、restart、
charge ledgerにも残り、設定した累積上限のどちらかを超えるとRUNを停止します。
これは局所的な数値回避策であり、物理境界条件そのものの代替ではありません。
再試行の試行数と解決数は、policy にかかわらず `summary.txt` の
`multiple_box_events_retry_attempted` と `multiple_box_events_retry_resolved` に残ります。

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
| `particle_step_ambiguous_open_corner` | 1003 | potential barrier で複数 open face が同時に発生 |

旧名`particle_step_unsupported_barrier_corner`は、code 1003の互換aliasとして残します。

これらを「命中なし」とみなすと粒子が表面を通り抜けるため、追跡は fail closed です。OpenMP では
particle/step 番号が最小の失敗を選び、MPI では失敗 rank と状態を共有して全 rank が同じ
batch/rank/particle/step/status を報告します。photo raycast も species/ray/bounce について同じ方針です。

`grid_stalled` の内部原因は `BEACH_COLLISION_DIAGNOSTICS=1` で確認できます。失敗した DDA 分岐名、
`p0` / `p1`、grid 範囲、cell index、`t_cur` / `t_next` / `t_delta` などを標準エラーへ出します。
この環境変数は診断出力だけを有効にし、衝突判定を変更しません。

## 衝突位置と帯電結果を収束させる

1. `dt`を半分にして軌道線分が衝突する要素と位置が安定するか確認する。
2. 薄い面、edge/corner近傍、三角形とほぼ平行な軌道線分を含む小ケースを作る。
3. reflect後とperiodic wrap後の残り時間内にmeshへ当たる軌道を確認する。
4. periodic seamで`pos`と`pos_wrapped`が意図した要素へ対応するか確認する。
5. `multiple_box_events`が出る場合、まず`dt`を小さくする。`upper_panel_fourier`を使う場合は解決率と成立域、
   soft discardを使う場合はdiscard率と絶対電荷が結論を左右しないことを確認する。
6. mesh refinementで最初のhitと最終的な表面電荷分布が収束するか確認する。

## Code reference

- box境界通過とmesh衝突の順序、残り時間の再積分: [`bem_particle_stepper.f90`](../src/runtime/simulator/bem_particle_stepper.f90)
- 三角形、grid DDA、periodic image: [`bem_collision.f90`](../src/physics/bem_collision.f90)
- box face検出と作用: [`bem_boundary.f90`](../src/physics/bem_boundary.f90)
- collision grid構築: [`bem_mesh.f90`](../src/mesh/bem_mesh.f90)
- 衝突・境界イベントの順序と8回上限のテスト: [`test_particle_stepper.f90`](../tests/fortran/test_particle_stepper.f90)
- 境界単体テスト: [`test_boundary.f90`](../tests/fortran/test_boundary.f90)
- 衝突・periodic・fail-closedテスト: [`test_dynamics_basic.f90`](../tests/fortran/test_dynamics_basic.f90)

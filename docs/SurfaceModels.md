title: 表面電荷更新

Lang: [日本語](SurfaceModels.md) | [English](SurfaceModels.en.md)

# 表面電荷更新

表面に吸収された粒子の電荷は、batch中は差分として保持します。表面から粒子を放出するsourceを使う場合は、
放出に伴う表面側の反作用電荷も同じ差分へ加えます。batch末尾にこれらを`q_elem`へまとめてcommitし、
commit後のsurface model処理を反映した電荷を次のbatchのfield sourceにします。

## batch内の更新順

1. batch開始時の`q_elem`から、batch中に固定する電場・電位を構成する。
2. 各粒子の最初のmesh hitへ$q_pw_p$をthread-localに加える。
3. 表面放出sourceの反作用電荷を`photo_emission_dq`へ加える。
4. OpenMP threadの差分を足し、MPI all-reduceでglobal `dq`を作る。
5. `q_elem <- q_elem + dq`を実行する。
6. conductorがあれば、object総電荷を保って等電位化する。
7. commit前後の正味差分と`tol_rel` metricを計算する。

同じbatch内の後続粒子は、手順2や3で生じた電荷を場として見ません。電荷更新がfieldへ現れるのは次batch開始時の
場の更新です。このlagが`batch_duration`依存性を作るため、[batch幅と安定性](BatchDurationStability.html)で
収束確認します。

## 保存量と符号

`q_elem(i)`は要素$i$の総電荷[C]です。暗黙のP0 panel離散化でも保存量は面密度ではなく、場評価時だけ

$$
\sigma_i=\frac{q_i}{A_i}
$$

へ変換します。macro粒子$p$が要素$i$へ吸収されたときの差分は

$$
\Delta q_i\mathrel{+}=q_pw_p
$$

です。electronなら負、正ionなら正の電荷を堆積します。粒子は吸収後に追跡から除かれます。

collisionに使うordered triangleと、fieldのone-sided traceに使う`elem_vacuum_sign`は同じmesh geometryから作ります。
表面modelのためにtriangle winding自体を書き換えません。

## surface model

| `surface_model` | commit後の処理 | 現行の位置付け |
| --- | --- | --- |
| `insulator` | hit要素に電荷を保持 | v1.0の中心model |
| `conductor` | `mesh_id`ごとに総電荷を保存して等電位化 | `field_bc_mode="free"`のみ |
| `dielectric` | 電荷を保持し`epsilon_r`をmetadata出力 | polarization未実装 |

## insulator accumulation

`insulator`はcommit後の再配分を行いません。吸収または放出で変化した各要素の電荷をその要素に保持します。

現行modelは、表面内の横方向伝導、bulkへの漏洩、有限抵抗によるrelaxation、二次電子放出、specular/diffuse反射を扱いません。
v1.0のinteractionはabsorptionが基本であり、これらの効果は結果に含まれません。

## floating conductor

`surface_model="conductor"`要素は`mesh_id`ごとに一つのfloating conductor objectを作ります。grounded potentialを
与える境界ではありません。粒子差分をいったんcommitした後、objectごとの総電荷$Q_g^\mathrm{before}$を保存しながら、
同じobject内の要素重心potentialを等しくします。

未知量は、全conductor要素の電荷$q_j$と各groupのscaled equipotential $V_g$です。要素$i$がgroup $g(i)$なら

$$
\sum_jA_{ij}q_j-V_{g(i)}=-\phi_i^\mathrm{fixed}
$$

を課します。source要素$j$を単位総電荷のP0 triangleとしたpotential coefficientは

$$
A_{ij}=\frac{1}{A_j}\int_{T_j}
\frac{1}{|\mathbf c_i-\mathbf y|}\,dA_{\mathbf y}
$$

です。$\mathbf c_i$はtarget要素重心、$A_j$と$T_j$はsource要素の面積と三角形です。
自己項を含め、解析的P0 panel potentialをprincipal-value側規約で評価します。
$\phi_i^\mathrm{fixed}$はnon-conductor電荷と一様外部fieldが作るpotentialを`k_coulomb`で割った量です。

さらに、各groupに総電荷保存

$$
\sum_{i\in g}q_i=Q_g^\mathrm{before}
$$

を加え、$N_\mathrm{cond}+N_\mathrm{group}$次のdense square systemを部分pivot付きGauss消去で解きます。
解いた$q_i$でconductor要素だけを置換します。したがってconductor relaxationはobject間の電荷を移さず、
non-conductor要素も変更しません。

このmodelは、要素重心collocationとP0 triangle influence matrixを使うfloating conductorです。
periodic/outer fieldとは併用できず、現行実装は`field_bc_mode="free"`だけを受理します。
要素細分化に対するobject potential・電荷分布の収束を確認してください。

## dielectric metadata

`surface_model="dielectric"`と`epsilon_r`は、現行版ではgeometry/material identityを出力へ残すmetadataです。
誘電率interface条件、法線$\mathbf D$のjump、polarization charge、内部fieldを解きません。
`epsilon_r`を指定しても場や電荷更新はinsulatorから変わらないため、誘電分極を含む結果として解釈しないでください。

## OpenMPとMPI commit

particle loopは`dq_thread(nelem,nthreads)`へ吸収電荷を集計します。要素ごとのatomic updateを避け、loop終了後にthread軸を
sumします。光電子の`photo_emission_dq`を加えたlocal差分をMPI all-reduceし、全rankが同じglobal `q_elem`を持ちます。

conductor relaxationはall-reduce後の同じmesh stateへ各rankで決定論的に適用します。collision queryやphoto raycastが
不完全statusを返したbatchはcommitへ進まず、部分的な粒子配列や放出差分を使いません。

## `tol_rel`

conductor relaxationまで含むcommit前後の正味差分を

$$
\Delta\mathbf q=\mathbf q^\mathrm{after}-\mathbf q^\mathrm{before}
$$

として

$$
\mathrm{tol\_rel\ metric}
=\frac{\|\Delta\mathbf q\|_2}{\max(\|\mathbf q^\mathrm{after}\|_2,q_\mathrm{floor})}
$$

を計算します。現行仕様では、`tol_rel`はmonitor/output metricであり、early-stop条件ではありません。

光電子放出に伴う反作用電荷の符号と放出粒子の追跡は
[光電子の放出とライフサイクル](PhotoelectronEmission.html)、粒子種別の電荷収支、履歴、最終出力、再開用fileは
[出力ファイルを調べる](OutputGuide.html)で説明します。

## Code reference

- particle absorptionとbatch commit: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- conductor relaxation facade: [`bem_surface_models.f90`](../src/physics/bem_surface_models.f90)
- floating-conductor solver: [`bem_surface_models_conductor.f90`](../src/physics/bem_surface_models_conductor.f90)
- batch statistics: [`bem_simulator_stats.f90`](../src/runtime/simulator/bem_simulator_stats.f90)

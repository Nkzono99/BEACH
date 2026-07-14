title: 表面電荷更新

Lang: [日本語](SurfaceModels.md) | [English](SurfaceModels.en.md)

# 表面電荷更新

表面に吸収された粒子の電荷と、粒子放出に伴う表面側の反作用電荷は、batch中は差分として保持します。
batch末尾にこれらを`q_elem`へまとめてcommitします。commit後のsurface model処理を反映した電荷が、次のbatchのfield sourceになります。

## batch内の更新順

1. batch開始時の`q_elem`から、batch中に固定する電場・電位を構成する。
2. 各粒子の最初のmesh hitへ$q_pw_p$をthread-localに加える。
3. 光電子放出元の反作用電荷を`photo_emission_dq`へ加える。
4. OpenMP threadの差分を足し、MPI all-reduceでglobal `dq`を作る。
5. `q_elem <- q_elem + dq`を実行する。
6. conductorがあれば、object総電荷を保って等電位化する。
7. commit前後の正味差分と`tol_rel` metricを計算する。

同じbatch内の後続粒子は、手順2や3で生じた電荷を場として見ません。電荷更新がfieldへ現れるのは次batch開始時の
場の更新です。このlagが`batch_duration`依存性を作るため、[batch幅と安定性](BatchDurationStability.html)で
収束確認します。

## 保存量と符号

`q_elem(i)`は要素$i$の総電荷[C]です。`triangle_p0` kernelでも保存量は面密度ではなく、場評価時だけ

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

## 光電子放出

`photo_raycast`で`deposit_opposite_charge_on_emit=true`なら、放出元$i$へ

$$
\Delta q_{i,\mathrm{emit}}=-q_pw_p
$$

を加えます。これは後続collisionとは別の`photo_emission_dq`へ集計し、同じbatch commitに統合します。

放出粒子が要素$j$へ戻ると、通常の吸収として$+q_pw_p$を$j$へ加えます。同じ要素へのreturnは放出電荷を相殺し、
別要素へのreturnは表面間の電荷移送になります。ray weight、reduced escape、tracked outer returnの関係は
[光電子の放出とライフサイクル](PhotoelectronEmission.html)で説明します。

## floating conductor

`surface_model="conductor"`要素は`mesh_id`ごとに一つのfloating conductor objectを作ります。grounded potentialを
与える境界ではありません。粒子差分をいったんcommitした後、objectごとの総電荷$Q_g^\mathrm{before}$を保存しながら、
同じobject内の重心potentialを等しくします。

未知量は、全conductor要素の電荷$q_j$と各groupのscaled equipotential $V_g$です。要素$i$がgroup $g(i)$なら

$$
\sum_jA_{ij}q_j-V_{g(i)}=-\phi_i^\mathrm{fixed}
$$

を課します。potential coefficientは

$$
A_{ij}=
\begin{cases}
1/\epsilon, & i=j,\ \epsilon>0,\\
2\sqrt\pi/h_i, & i=j,\ \epsilon=0,\\
1/\sqrt{|\mathbf c_i-\mathbf c_j|^2+\epsilon^2}, & i\ne j
\end{cases}
$$

です。$\mathbf c_i$は要素重心、$h_i$は要素長さscale、$\epsilon$はfield softeningです。
$\phi_i^\mathrm{fixed}$はnon-conductor電荷と一様外部fieldが作るpotentialを`k_coulomb`で割った量です。

さらに、各groupに総電荷保存

$$
\sum_{i\in g}q_i=Q_g^\mathrm{before}
$$

を加え、$N_\mathrm{cond}+N_\mathrm{group}$次のdense square systemを部分pivot付きGauss消去で解きます。
解いた$q_i$でconductor要素だけを置換します。したがってconductor relaxationはobject間の電荷を移さず、
non-conductor要素も変更しません。

このmodelは、centroid point-potential coefficientを使う簡易floating conductorです。triangle P0による厳密なconductor BEM境界積分ではありません。
periodic/outer fieldとは併用できず、現行実装は`field_bc_mode="free"`だけを受理します。
要素細分化とsofteningに対するobject potential・電荷分布の収束を確認してください。

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

## charge ledger

`charge_ledger.csv`はspeciesごとにsigned chargeとcountを区別します。

| flux | 意味 |
| --- | --- |
| `injected_from_remote` | `volume_seed`/`reservoir_face`から入った電荷 |
| `emitted_from_surface` | `photo_raycast`でsurfaceから出たtracked電荷 |
| `absorbed_on_surface` | meshへ吸収された電荷 |
| `escaped_to_infinity` | open/outer modelでescapeした電荷 |
| `discarded_unresolved` | `max_step`で未解決のまま破棄した電荷 |
| `interface_outward_gross` / `returned_gross` | outer ownership面の往復量 |

transactional residualは、surface、local flight、outer flight、unresolved stockの前後差と、remote injection、escape、discardを合わせて計算します。
surface emissionとabsorptionはsurface/flight stock間の内部移送です。そのため、residual式で独立したexternal sourceとしては数えません。

species間で正負が相殺し得るresidualとは別に、`discarded_unresolved_abs`は
$\sum_s|Q_{s,\mathrm{discard}}|$を確認します。residualが小さくても未解決discardが大きい計算は受理しません。

## history、final output、restart

`history_stride>0`なら

$$
(\texttt{stats.batches}-1)\bmod\texttt{history\_stride}=0
$$

のbatchで`charge_history.csv`を書きます。batch 1は常に対象です。`write_potential_history=true`なら同じstrideで現在の
`q_elem`から電場・電位を更新し、要素重心の`potential_history.csv`も書きます。

主な確認fileは`charges.csv`、`charge_history.csv`、`charge_ledger.csv`、`summary.txt`、`mesh_sources.csv`です。
restart時には、mesh要素数、MPI world size、stats/chargeが有限値であること、ordered mesh/species/model fingerprint、ledger stockを検証します。
必要なouter profileも検証対象です。`output.resume=true`で必須checkpointが欠けている場合は、新規runへfallbackせず停止します。

## Code reference

- particle absorptionとbatch commit: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- conductor relaxation: [`bem_surface_models.f90`](../src/physics/bem_surface_models.f90)
- transactional charge ledger: [`bem_charge_ledger.f90`](../src/runtime/coupling/bem_charge_ledger.f90)
- batch statistics: [`bem_simulator_stats.f90`](../src/runtime/simulator/bem_simulator_stats.f90)
- history output: [`bem_simulator_io.f90`](../src/runtime/simulator/bem_simulator_io.f90)
- final output: [`bem_output_writer.f90`](../src/runtime/bem_output_writer.f90)
- restart validation: [`bem_restart.f90`](../src/runtime/bem_restart.f90)

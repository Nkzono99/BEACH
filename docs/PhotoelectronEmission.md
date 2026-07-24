title: 光電子の放出とライフサイクル

Lang: [日本語](PhotoelectronEmission.md) | [English](PhotoelectronEmission.en.md)

# 光電子の放出とライフサイクル

`source_mode="photo_raycast"`は、照射rayが最初に命中した表面から粒子を放出します。放出位置はmeshとraycastで決め、
放出に伴う表面電荷はsurface stateで保存します。生成後の光電子は通常の粒子と同じであり、batch内で固定された電場、
Boris更新、衝突判定、box境界を使います。

## 放出から再吸収までを同じbatchで追う

1. box面上の照射開口からrayを発射する。
2. box境界条件を適用しながら最初のmesh hitを探す。
3. 命中要素のplasma側法線から光電子を生成する。
4. 必要なら放出元要素へ$-qw$を記録する。
5. 通常粒子として追跡し、box境界へ達した後は共通のescape/return処理へ渡す。
6. 放出電荷と吸収電荷をbatch末尾に表面へcommitする。

放出と再吸収は同じbatchで起こり得ますが、途中で電場・電位を更新しません。したがって、放出が作る
正味表面電荷は次batchから場へ反映されます。

## 照射rayで放出面を決める

`inject_face`と`pos_low`/`pos_high`がbox面上の矩形開口を定めます。`ray_direction`は開口からbox内へ向く必要があり、
省略時は面の内向き法線です。開口面積を$A$、内向き法線を$\mathbf n_\mathrm{in}$、正規化したray方向を
$\hat{\mathbf d}$とすると、rayに垂直な投影面積は

$$
A_\mathrm{proj}=A\left|\hat{\mathbf d}\cdot\mathbf n_\mathrm{in}\right|
$$

です。rayの始点は開口内で一様sampleします。

rayは現在位置から次のbox面までのsegmentごとに最初のtriangle hitを探します。

| 到達したbox面 | rayの処理 |
| --- | --- |
| `open` | box外へ出て、そのrayは放出なしで終了 |
| `reflect` | 方向の法線成分を反転して継続 |
| `periodic` | 反対側へwrapして継続 |

`field_bc_mode="periodic2"`ではperiodic imageを含む最初のhitを探し、放出位置をprimary cellへwrapします。
hit要素は計算box内にある要素に限定します。`raycast_max_bounce`を越えてもhitしないrayは、粒子を生成しません。
衝突queryがimage上限、index範囲、DDA停止などの不完全statusを返した場合は、batchを停止します。

## ray 1本へ放出電流を割り当てる

放出電流密度を$J_\mathrm{emit}>0$、実粒子電荷を$q\ne0$、全rank合計のray数を$N_\mathrm{ray}$とすると、
hitしたrayが生成するmacro粒子重みは

$$
w_\mathrm{hit}
=\frac{J_\mathrm{emit}A_\mathrm{proj}\,\Delta t_\mathrm{batch}}
{|q|N_\mathrm{ray}}
$$

です。missしたrayは粒子を作らないため、遮蔽や見かけ面積はhit率として放出量へ入ります。
`w_particle`や`target_macro_particles_per_batch`は`photo_raycast`には指定しません。

`rays_per_batch`は物理放出量ではなくray積分のsampling数です。値を増やすと$w_\mathrm{hit}$が小さくなり、
照射可視率と放出位置のMonte Carlo noiseを減らせます。結果はray数に対して収束確認します。

## 表面法線から放出状態を作る

triangleの格納法線がray進行方向を向いている場合は反転し、入射rayと反対側を向く放出法線$\mathbf n_s$を作ります。
位置はhit点から$\mathbf n_s$方向へ$10^{-12}$ mずらし、直後に同じ要素へ再衝突することを避けます。

温度から$\sigma=\sqrt{k_\mathrm{B}T/m}$を求め、局所基底
$(\mathbf n_s,\mathbf t_1,\mathbf t_2)$で速度をsampleします。

- 法線速度は、drift `normal_drift_speed`を持つflux-weighted half-range Maxwell分布。
- 接線2成分は平均0、標準偏差$\sigma$のGaussian。
- Gaussian samplingは$6\sigma$で切る。

法線速度は正なので、生成直後の粒子は照射側へ表面から離れます。その後に戻るかescapeするかは、tracked orbitと
共通のbox境界またはouter transferが決めます。

## 放出・再吸収・escapeの電荷収支を確認する

`deposit_opposite_charge_on_emit=true`なら、放出元要素$i$へ

$$
\Delta q_{i,\mathrm{emit}}=-q w
$$

を加えます。電子では$q<0$なので表面には正電荷が残ります。この差分は`photo_emission_dq`として衝突堆積と別に集計し、
MPI all-reduce後に同じbatch commitへ加えます。

放出粒子が要素$j$へ再吸収されると、通常の吸収として$+qw$を$j$へ堆積します。同じ要素へ戻れば放出と吸収が相殺し、
別の要素へ戻れば表面内の正味電荷移送になります。現行のinsulator modelはその後の表面伝導を行いません。

## 生成後は共通のescape / return処理を使う

ray hitで生成する光電子の重みは常に$w_\mathrm{hit}$です。放出時にescape率を掛けて粒子重みを減らす設定はありません。
表面へ戻った粒子は通常の衝突として再吸収し、open面へ達した粒子にはreservoir粒子や`volume_seed`粒子と同じ
`external_boundary.ordinary_open`またはouter particle modeを適用します。

有限boxで外部領域を解かない場合は`external_boundary.ordinary_open.model="potential_barrier"`が
通過点の電位と法線運動エネルギーから反射またはescapeを決めます。自己整合な外部sheathを使う場合は
`external_boundary.field.model`と`external_boundary.particles.mode`が正本です。
scalar barrier、1D outer profile return、unified 3D explicit orbitの違いは
[粒子のescapeとreturn](ParticleEscapeReturn.html)で説明します。

tracked outer transferを使う`photo_raycast` speciesでは、放出と帰還の電荷収支を閉じるため
`deposit_opposite_charge_on_emit=true`を指定します。

## 光電子をouter plasmaの平均密度へ含める

`external_boundary.field.photoelectron_density_model="kinetic_mean"`は、最初の負電荷`photo_raycast` speciesの温度と放出電流密度を、
平面平均sourceとして1D Poisson密度モデルへ加えます。outgoing populationと、turning後のreturning populationが、
outer領域の空間電荷に寄与します。この平均密度モデルは、個々のtracked粒子の表面吸収を置き換えません。
統計的なreturn電荷を、別途表面へdepositすることもありません。

生成後のtracked光電子をz-high interfaceからouter領域へ渡す場合も、粒子sourceに依存しない共通の
escape/return処理を使います。外部flightをglobal timeへ加えない準定常近似と3D explicit orbitは
[粒子のescapeとreturn](ParticleEscapeReturn.html)、対応する場の作り方は[外部プラズマモデル](OuterPlasmaModels.html)で説明します。

`external_boundary.field.kinetic_closure="zhao_charge_driven"`は、蓄積電荷が決めるinterface電場を保つ
自己整合な1D outer profileを構築します。`external_boundary.particles.mode="zhao_queue"`では、さらに
tracked outer inventoryから光電子populationを閉じます。`mode="same_batch"`または`"zhao_queue"`では、このprofileが
z-high interfaceでの流入とreturn / escapeを決めます。

## 光電子放出の収束を確認する

`rays_per_batch`を増やし、hit率、放出電流、帯電分布が収束することを確認します。再吸収位置も評価する場合は、
`dt`を小さくして結果が変わらないことを確認します。放出、吸収、escapeを含む粒子種別の電荷収支と、
outer return 固有の診断値は[出力ファイルを調べる](OutputGuide.html)で確認します。

## Code reference

- ray伝播、hit、放出速度と重み: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- 放出電荷差分とsource生成: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- kinetic mean光電子密度モデル: [`bem_outer_plasma_kinetic.f90`](../src/physics/outer_plasma/bem_outer_plasma_kinetic.f90)
- kinetic meanのruntime組立: [`bem_outer_plasma_kinetic_runtime.f90`](../src/runtime/bem_outer_plasma_kinetic_runtime.f90)

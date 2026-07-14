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
5. 通常粒子として追跡し、再吸収、無限遠escape、outer returnを処理する。
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
- Zhao closureなどが$v_{\min}$を与える場合、法線速度は$v_n\ge v_{\min}$。
- Gaussian samplingは$6\sigma$で切る。

法線速度は正なので、生成直後の粒子は照射側へ表面から離れます。その後に戻るかescapeするかは、tracked orbitまたは
選択したreduced closureが決めます。

## 放出と再吸収を表面電荷へ記録する

`deposit_opposite_charge_on_emit=true`なら、放出元要素$i$へ

$$
\Delta q_{i,\mathrm{emit}}=-q w
$$

を加えます。電子では$q<0$なので表面には正電荷が残ります。この差分は`photo_emission_dq`として衝突堆積と別に集計し、
MPI all-reduce後に同じbatch commitへ加えます。

放出粒子が要素$j$へ再吸収されると、通常の吸収として$+qw$を$j$へ堆積します。同じ要素へ戻れば放出と吸収が相殺し、
別の要素へ戻れば表面内の正味電荷移送になります。現行のinsulator modelはその後の表面伝導を行いません。

## box外の挙動を一つのモデルで表す

box外へ向かう光電子には、次のいずれか一つの近似を適用します。選択した構成によって、生成するtracked重みと、
return位置やflight timeを表現する範囲が決まります。

| 構成 | 生成するtracked重み | returnの表現 |
| --- | --- | --- |
| `photo_escape_model="none"`、outer transferなし | $w_\mathrm{hit}$ | box内を追跡し、open面では通常のescape |
| `photo_escape_model="boltzmann_cutoff"` | $f_\mathrm{esc}w_\mathrm{hit}$ | 非escape成分を生成しない即時reduced closure |
| 1D outer profile return | $w_\mathrm{hit}$ | interfaceを出た個々の粒子を保存エネルギーでreturn/escape判定 |
| unified 3D explicit orbit | $w_\mathrm{hit}$ | zero/nonzero modeを含む3D outer場で個別軌道を追跡 |

legacy Boltzmann cutoffは、放出元のprimary self項を除いた局所電位$\phi_\mathrm{emit}$から

$$
f_\mathrm{esc}=
\exp\left[-\frac{|q|\max(\phi_\mathrm{emit}-\phi_\infty,0)}{k_\mathrm{B}T_\mathrm{PE}}\right]
$$

を求め、粒子重みをこの係数倍します。$T_\mathrm{PE}=0$で正の障壁がある場合、この係数は0です。
returning粒子の軌道、再吸収位置、flight timeは扱いません。非escape成分を戻す要素も求めません。
放出元へ置く逆符号電荷には、減衰後の粒子と同じ重みを使います。

個別returnを使う構成は`deposit_opposite_charge_on_emit=true`を要求し、legacy `photo_escape_model`とは併用しません。
escape cutoffとtracked returnは、どちらか一方だけを適用します。

## outer plasmaの平均密度と個別軌道を結び付ける

`outer_plasma.photoelectron_closure="kinetic_mean"`は、最初の負電荷`photo_raycast` speciesの温度と放出電流密度を、
平面平均sourceとして1D Poisson密度closureへ加えます。outgoing populationと、turning後のreturning populationが、
outer領域の空間電荷に寄与します。この平均closureは、個々のtracked粒子の表面吸収を置き換えません。
統計的なreturn電荷を、別途表面へdepositすることもありません。

`individual_return`または`kinetic_mean + kinetic_1d_profile_return`では、z-high interfaceを横切るtracked光電子を
保存エネルギーでescape/returnへ分類します。interfaceのoutward/returned chargeと放出速度histogramはdiagnosticです。
外部flightをglobal timeへ加えない準定常近似と3D explicit orbitは
[粒子のescapeとreturn](ParticleEscapeReturn.html)、対応する場の作り方は
[外部プラズマモデル](OuterPlasmaModels.html)で説明します。

Zhao系は、branchに応じて放出電流密度、法線cutoff、driftを与える注入closureです。tracked粒子は
Zhao profileの$E(z)$ではなく、通常の粒子追跡で使う、batch内で固定された電場中を進みます。

## 放出・再吸収・escapeの電荷収支を確認する

粒子種別の電荷収支では、surfaceからの放出、surfaceへの吸収、無限遠escape、未解決破棄を区別します。
outer interfaceのoutward/returned gross chargeも、それぞれ記録します。少なくとも次の量を比較します。

- `rays_per_batch`を増やしたときのhit率、放出電流、帯電分布。
- `dt`を小さくしたときの再吸収位置とescape/return率。
- 表面の放出$-qw$、tracked粒子の吸収/escape、batch間surface stockを合わせた電荷収支。
- outer return使用時の最大flight time、frozen-field ratio、energy error。

## Code reference

- ray伝播、hit、放出速度と重み: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- reduced escape係数と放出電荷差分: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- tracked-return互換性検証: [`bem_app_config_parser.f90`](../src/config/app_config_parser/bem_app_config_parser.f90)
- 通常追跡、outer transfer、電荷収支集計: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- kinetic mean photoelectron closure: [`bem_outer_plasma_photoelectron.f90`](../src/physics/outer_plasma/bem_outer_plasma_photoelectron.f90)

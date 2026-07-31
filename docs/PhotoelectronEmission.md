title: 光電子の放出とライフサイクル

Lang: [日本語](PhotoelectronEmission.md) | [English](PhotoelectronEmission.en.md)

# 光電子の放出とライフサイクル

`source_mode="photo_raycast"` は照射 ray が最初に命中した表面から粒子を放出します。このページは
raycast、放出量・速度、放出電荷、closed PE の `neutral_return` を説明します。生成後は通常粒子と同じ
固定場、Boris 更新、衝突判定、box 境界を使います。

## 放出から再吸収までを同じbatchで追う

1. box面上の照射開口からrayを発射する。
2. box境界条件を適用しながら最初のmesh hitを探す。
3. 命中要素のplasma側法線から光電子を生成する。
4. 必要なら放出元要素へ$-qw$を記録する。
5. 通常粒子として追跡し、box境界へ達した後は共通のescape/return処理へ渡す。
6. 放出電荷と吸収電荷をbatch末尾に表面へcommitする。

放出と再吸収は同じ batch で起こり得ますが、場は途中で更新しません。正味表面電荷は次 batch から場へ反映されます。

## 照射rayで放出面を決める

`inject_face`と`pos_low`/`pos_high`がbox面上の矩形開口を定めます。`ray_direction`は開口からbox内へ向く必要があり、
省略時は面の内向き法線です。開口面積を$A$、内向き法線を$\mathbf n_\mathrm{in}$、正規化したray方向を
$\hat{\mathbf d}$とすると、rayに垂直な投影面積は

$$
A_\mathrm{proj}=A\left|\hat{\mathbf d}\cdot\mathbf n_\mathrm{in}\right|
$$

です。rayの始点は開口内で一様sampleします。

rayは現在位置から次のbox面までのsegmentごとに最初のtriangle hitを探します。

| 到達した box 面 | ray の処理 |
| --- | --- |
| 非周期面 | box 外へ出て、その ray は放出なしで終了 |
| `domain.periodic_axes` の面 | 反対側へ wrap して継続 |

`field_boundary.mode="periodic2"` では periodic image を含む最初の hit を探し、放出位置を primary cell へ wrap します。
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

`rays_per_batch` は物理放出量ではなく ray 積分の sample 数です。増やすと $w_\mathrm{hit}$ が小さくなり、
照射可視率と放出位置の Monte Carlo noise を減らせます。結果は ray 数に対して収束確認します。

## 表面法線から放出状態を作る

triangleの格納法線がray進行方向を向いている場合は反転し、入射rayと反対側を向く放出法線$\mathbf n_s$を作ります。
位置はhit点から$\mathbf n_s$方向へ$10^{-12}$ mずらし、直後に同じ要素へ再衝突することを避けます。

温度から$\sigma=\sqrt{k_\mathrm{B}T/m}$を求め、局所基底
$(\mathbf n_s,\mathbf t_1,\mathbf t_2)$で速度をsampleします。

- 法線速度は、drift `normal_drift_speed`を持つflux-weighted half-range Maxwell分布。
- 接線2成分は平均0、標準偏差$\sigma$のGaussian。
- Gaussian samplingは$6\sigma$で切る。

法線速度は正なので、生成直後の粒子は照射側へ表面から離れます。その後の再吸収、escape、局所反射は
tracked orbit と共通の box 境界が決めます。

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

ray hit で生成する光電子の重みは常に $w_\mathrm{hit}$ です。放出時に escape 率を掛ける設定はありません。
表面へ戻れば通常衝突として再吸収し、open 面へ達すれば他の source と同じ
`particle_boundary.ordinary_open_model` を適用します。

### 光電子だけを注入面で閉じる

光電子の表面内再分配を閉じた軌道で評価する場合は、光電子 species へ局所反射を指定します。

```toml
[particle_boundary]
z_high = "open"
ordinary_open_model = "escape"

[[particles.species]]
species_key = "photoelectron"
source_mode = "photo_raycast"
inject_face = "z_high"
deposit_opposite_charge_on_emit = true
surface_charge_closure = "neutral_return"

[particles.species.boundary]
z_high = "reflect"
```

この例は z-high 通過時に法線速度だけを反転し、接線速度を保存して残り step を再積分します。
`[particles.species.boundary]` の既定 `inherit` を使う ambient species は global の open 契約に従います。
species 境界は 6 面すべてに `inherit`、`open`、`reflect`、`redistributed_reflect` を指定でき、closed PE では
`inject_face` と同じ面の有効作用が `reflect` または `redistributed_reflect` でなければなりません。
`domain.periodic_axes` の面は上書きできません。

通常の `reflect` は接線位置も維持します。境界から戻る光電子の面内位置を一様化する場合だけ、上のbaseline値を
次のように置き換えます。

```toml
[particles.species.boundary]
z_high = "redistributed_reflect"
```

`redistributed_reflect` は速度には通常の反射を適用し、単一面では面内2軸の位置だけをbox spanの両端guardを
除く範囲から一様再標本化します。
これは [Zimmerman et al. (2016)](https://doi.org/10.1002/2016JE005049) のtop-boundary PE returnで用いられた
水平位置randomizationを、任意の非周期面と同時face eventへ一般化した選択肢です。自己整合な外部sheathを追加する
ものではありません。同時eventの規則は[粒子の衝突・境界イベント](ParticleEvents.html)を参照してください。

species 境界の反射だけなら軌道を閉じるだけで、`max_step` までに戻らない粒子は未解決のままです。
`surface_charge_closure="neutral_return"`を加えると、1 batchの光電子放出電荷$S<0$と解決済み吸収電荷$R<0$を
MPI全体で測り、各帰還先depositを$S/R$倍します。放出元の反作用電荷と合わせた光電子の表面総電荷増分は
厳密に0になり、未帰還粒子は解決済み帰還先と同じ分布を持つと近似されます。

rawの`absorbed_on_surface_C`と`discarded_unresolved_C`は置き換えません。補正量、
`neutral_return_weight_scale`、`neutral_return_unresolved_fraction`を`charge_ledger.csv`へ別に記録します。
放出があるのに解決済み帰還がない場合、実escape、`soft_discard`、符号不整合は停止します。
未帰還率が5%を超える場合も、このclosureの固定適用範囲外として補正せず停止します。

完全反射は有限 box の注入面に人工的な鏡を置く試験条件であり、自己整合な sheath や準中性性を解きません。
`neutral_return` も未帰還軌道を解かず、正味光電子電流を 0 とする統計的 closure です。
`abs(weight_scale-1)` と未帰還率が十分小さくなるよう `max_step`、`dt`、ray 数、batch 幅を収束させます。
統合した場・流入・電位基準は
[periodic2有限画像構成](FinitePeriodicConfiguration.html)にあります。

通常の open 面に使う `particle_boundary.ordinary_open_model` と closed PE の使い分けは
通常の open 面との使い分けは[粒子の escape と局所 return](ParticleEscapeReturn.html)にまとめています。

## 光電子放出の収束を確認する

`rays_per_batch`を増やし、hit率、放出電流、帯電分布が収束することを確認します。再吸収位置も評価する場合は、
`dt`を小さくして結果が変わらないことを確認します。放出、吸収、escapeを含む粒子種別の電荷収支と、
closed PEの補正量と未解決率は[出力ファイルを調べる](OutputGuide.html)で確認します。

## Code reference

- ray伝播、hit、放出速度と重み: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- 放出電荷差分とsource生成: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)

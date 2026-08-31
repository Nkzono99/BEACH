title: periodic2有限画像構成

Lang: [日本語](FinitePeriodicConfiguration.md) | [English](FinitePeriodicConfiguration.en.md)

# periodic2有限画像構成

このページは、境界 reservoir + closed PE によってレゴリス帯電と光電子の表面内再分配を調べる
`periodic2` 構成の正本です。太陽風 VDF は z-high 面で定義し、光電子だけを閉じ、電位は z-high 面平均からの
差として読みます。このbaselineではmatching-plane outer-sheath closureを有効にしません。

## まず二つの構成を区別する

| 構成 | 太陽風流入 | 光電子 | 電位基準 | 用途 |
| --- | --- | --- | --- | --- |
| **境界reservoir + closed PE** | z-highのVDFを補正せず使用 | z-high反射 + `neutral_return` | z-high面平均 | このページの基準構成。表面内再分配とbatch感度 |
| scalar barrier | `infinity_barrier` | 共通の`potential_barrier` | `phi_infty` | 一つのscalar障壁による比較 |

同じrunでこれらを重ねません。特に、closed PE構成には`infinity_barrier`、
`potential_barrier`を追加しません。

## 推奨する統合構成

完全な実行例は
[`examples/periodic2_closed_photoelectron.toml`](../examples/periodic2_closed_photoelectron.toml) です。

```toml
[sim]
field_solver = "fmm"
field_periodic_image_layers = 1
field_periodic_far_correction = "none"

[domain]
box_origin = [0.0, 0.0, 0.0]
box_size = [1.0e-4, 1.0e-4, 1.0e-3]
periodic_axes = ["x", "y"]

[field_boundary]
mode = "periodic2"

[particle_boundary]
z_low = "open"
z_high = "open"
ordinary_open_model = "escape"

[reservoir]
inflow_model = "source_vdf"
face_potential_grid_n = 5
```

`[domain]` が cell と周期 topology を、`[field_boundary]` が場の closure を定めます。現行 `periodic2` は
`periodic_axes=["x", "y"]` と非周期 z 軸の組合せだけを受理します。周期性は全粒子種に共通で、species や
`[particle_boundary]` から上書きできません。
旧 `[external_boundary]` は削除済みです。外部 reservoir 条件は `[reservoir]`、species 別流入面は
`[particles.species.boundary_inflow]`、非周期面の外向き粒子作用は `[particle_boundary]` に指定します。

太陽風 electron/ion には z-high の `boundary_inflow="reservoir"` を使います。密度、温度、drift は
**z-high 面の局所 boundary VDF**であり、表面電位による到達率・速度補正は行いません。

光電子speciesだけに閉じた表面電荷closureを指定します。

```toml
[[particles.species]]
species_key = "photoelectron"
source_mode = "photo_raycast"
deposit_opposite_charge_on_emit = true
surface_charge_closure = "neutral_return"

[particles.species.boundary]
z_high = "reflect"
```

species 境界の `z_high="reflect"` は非周期 z-high 通過時の法線速度だけを反転します。太陽風 species は既定の
`"inherit"` なので、`particle_boundary.z_high="open"` と `ordinary_open_model="escape"` に従います。

このbaselineは境界eventの面内位置を維持します。topからのreturn位置をx-y面内で一様化する感度比較では、
`z_high="redistributed_reflect"`へ置き換えられます。速度作用は同じで、位置だけを再標本化します。
これはclosed PEの別のreturn-destination modelであり、自己整合sheathではありません。背景と同時eventの規則は
[光電子の放出とライフサイクル](PhotoelectronEmission.html)を参照してください。

## 1 batchで行うこと

1. batch 開始時の表面電荷から有限画像を含む場 snapshot を作る。
2. z-high の局所 VDF から太陽風を、ray hit から光電子を注入する。
3. 固定した snapshot で全粒子を追跡し、光電子だけを z-high で反射する。
4. 解決済み帰還先分布で未帰還光電子を統計的に閉じ、全電荷差分を commit する。
5. 必要なら同じ commit 後の電荷から要素電位と z-high 面平均を出力する。

batch中に場は更新しません。したがって`batch_duration`は物理時間幅であると同時にexplicit charge updateの
幅です。closed PEは光電子の正味電流による発散を除きますが、太陽風帯電のbatch幅依存性は残ります。
[<sup>1</sup>](FieldSolvers.html)

## `neutral_return`が閉じる量

負電荷の光電子speciesについて、1 batchのsigned macro chargeを

$$
S=\sum_{\mathrm{emitted}}qw<0,\qquad
R=\sum_{\mathrm{resolved\ absorption}}qw<0
$$

とします。`neutral_return`は実測した帰還先depositを

$$
s_{\mathrm{return}}=\frac{S}{R}
$$

倍します。放出元へ残す反作用電荷は$-S$なので、

$$
(-S)+s_{\mathrm{return}}R=0
$$

となります。つまり、このspeciesは表面総電荷を変えず、放出元から解決済み帰還先への分布だけを作ります。
太陽風による総電荷更新は拘束しません。

これは未帰還粒子の帰還先分布が、同じbatchで帰還した粒子と同じだという統計的closureです。
`s_return`が1から大きく離れるrunは、未解決軌道ではなくclosure仮定が局所分布を決めています。
この近似は少数の長寿命粒子だけを閉じる固定契約であり、未帰還率が5%を超える場合は補正せず停止します。
上限を緩める設定は持たず、`max_step`、`dt`、box高さを見直して5%以下へ収束させます。

`neutral_return`は「全ての$k_\parallel=0$構造を数式的に除去する」操作ではありません。表面総電荷の
monopole増分を0にしますが、異なる高さの面へ移った電荷は平面平均された鉛直dipoleを作り得ます。

次の場合は batch を受理せず停止します。

- 放出があるのに解決済み帰還電荷が0。
- 光電子の実 escape または `soft_discard`。
- 非有限値または符号不整合。

`max_step`まで残った粒子だけが統計補正の対象です。raw未帰還電荷、補正量、係数、未帰還率は
`charge_ledger.csv`へ別々に残します。

## z-high面を電位基準として読む

`[domain]` があり、`output.write_files=true`、`output.write_potential_history=true`、
`output.history_stride>0` なら、`potential_history.csv` と同じ batch に `top_reference_history.csv` を出力します。

$$
\bar\phi_{\mathrm{top}}
=\frac{1}{N^2}\sum_{a,b}\phi(x_a,y_b,z_{\mathrm{high}}^-),\qquad
\phi_{\mathrm{rel}}(\mathbf r)=\phi(\mathbf r)-\bar\phi_{\mathrm{top}}.
$$

$N$ は `reservoir.face_potential_grid_n` です。標本点は全 periodic cell 上面の cell center です。
二つの履歴を`batch`でjoinし、各要素の`potential_V`から同じbatchの
`potential_mean_V`を引きます。

```toml
[output]
history_stride = 1000
write_mesh_potential = true
write_potential_history = true
```

この面平均は無限遠電位でもプラズマ電位でもなく、太陽風流入へ feedback しない診断値です。定数 gauge は
除けますが、box 高さ、zero mode、有限画像打切りへの依存性は残ります。`potential_std_V` や min/max の幅が
大きい場合、z-high を一つの reservoir 面とみなす近似は弱くなります。

`write_mesh_potential=true`または`write_potential_history=true`なら、最終状態のz-high統計を
`summary.txt`にも記録します。最終`mesh_potential.csv`へ相対化するときはこの値を使えます。履歴strideに
最終batchが含まれない場合でも、最終mesh電位へ別batchのtop値を流用しません。

## 有限画像の意味

`field_periodic_image_layers=N`は、primary cellの周囲$N$層までをfield sourceへ加えます。

| $N$ | 含むcell |
| --- | --- |
| 0 | primaryだけ、$1\times1$ |
| 1 | 周囲1層、$3\times3$ |
| 2 | 周囲2層、$5\times5$ |

`field_periodic_far_correction="none"`では、その外側をEwald和やcached operatorで補いません。
top-relative potentialにしても無限周期解にはなりません。画像層を増やし、目的量が変わらなくなることを確認します。
無限周期の nonzero mode が必要なら `field_periodic_far_correction="cached_kneq0"` を使い、
[periodic2静電場](PeriodicElectrostatics.html)の収束手順を確認します。

## 受理前に確認する順序

1. `abs(neutral_return_weight_scale-1)`と`neutral_return_unresolved_fraction`が十分小さいことを確認する。
2. `dt`を下げるときも`max_step*dt`を維持または増加させ、`rays_per_batch`も増やして帰還先分布を収束させる。
3. `batch_duration`を少なくとも$T,T/2,T/4$と変え、電荷・相対電位・力の履歴を比較する。
4. z-high を上下し、`reservoir.face_potential_grid_n` を増やして top 平均とばらつきを確認する。
5. image layerを$N,N+1,N+2$と増やす。
6. 太陽風のmacro粒子数と乱数seedに対する統計誤差を確認する。`redistributed_reflect`を使う場合は、
   return位置の追加samplingについても粒子数とseedを変えて確認する。

完全反射は人工的な上端 mirror であり、自己整合 sheath や準中性解ではありません。この構成が示すのは、
指定した境界太陽風 flux の下で正味光電子電流を 0 としたときの表面内再分配です。

## scalar barrierを比較に使う場合

scalar barrier 比較では closed PE 設定を外し、`reservoir.inflow_model="infinity_barrier"`、
`reservoir.phi_infty`、`particle_boundary.ordinary_open_model="potential_barrier"` を組み合わせます。これは face 平均
scalar で上流 VDF を補正し、各 open 通過点の energy で reflect/escape を決める比較モデルです。設定、式、制約は
[境界から粒子を流入させる](ReservoirInjection.html)と
[粒子の escape と return](ParticleEscapeReturn.html)にまとめています。

title: periodic2有限画像構成

Lang: [日本語](FinitePeriodicConfiguration.md) | [English](FinitePeriodicConfiguration.en.md)

# periodic2有限画像構成

このページは、外部1Dシースを解かずに局所的なレゴリス帯電と光電子による表面内再分配を調べる
`periodic2` 構成の正本です。太陽風VDFはz-high面で定義し、光電子だけを閉じます。電位は
無限遠ではなくz-high面平均に対する差として読みます。

この構成の狙いは、外部シースを近似的に再現することではありません。太陽風が作る正味表面電荷と、
光電子が作る局所的な電荷移送を、少ない仮定で安定に分けることです。

## まず三つの構成を区別する

| 構成 | 太陽風流入 | 光電子 | 電位基準 | 用途 |
| --- | --- | --- | --- | --- |
| **局所reservoir + closed PE** | z-highのVDFを補正せず使用 | z-high反射 + `neutral_return` | z-high面平均 | このページの基準構成。局所再分配とbatch感度 |
| scalar barrier | `infinity_barrier` | 共通の`potential_barrier` | `phi_infty` | 一つのscalar障壁による比較 |
| 無限周期 + outer plasma | 1D profileから写像 | profileでreturn/escape | $\phi_\infty=0$ | 自己整合な平均シースが必要な場合 |

同じrunでこれらを重ねません。特に、closed PE構成には`infinity_barrier`、
`potential_barrier`、`kinetic_1d` particle transferを追加しません。

## 推奨する統合構成

完全な実行例は
[`examples/periodic2_closed_photoelectron.toml`](../examples/periodic2_closed_photoelectron.toml)
です。中心となる設定は次です。

```toml
[sim]
field_solver = "fmm"
field_bc_mode = "periodic2"
field_periodic_image_layers = 1
field_periodic_far_correction = "none"
use_box = true
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
injection_face_phi_grid_n = 5

[external_boundary.field]
model = "none"

[external_boundary.particles]
mode = "local_source"
inflow_model = "source_vdf"

[external_boundary.ordinary_open]
model = "escape"
```

太陽風electron/ionには通常の`reservoir_face`を使います。ここで指定する密度、温度、driftは
**z-high面の局所boundary VDF**です。無限遠VDFではなく、表面電位による到達率・速度補正も行いません。

光電子speciesだけに閉じた表面電荷closureを指定します。

```toml
[[particles.species]]
species_key = "photoelectron"
source_mode = "photo_raycast"
deposit_opposite_charge_on_emit = true
z_high_boundary = "reflect"
surface_charge_closure = "neutral_return"
```

`z_high_boundary="reflect"`はz-high通過時の法線速度だけを反転します。太陽風speciesは既定の
`"inherit"`なので、global open境界からescapeします。

## 1 batchで行うこと

1. batch開始時の表面電荷から有限画像を含む電場・電位の空間分布を作る。
2. z-highの局所VDFから太陽風を注入する。
3. ray hitから光電子を放出し、放出元へ逆符号の表面電荷を記録する。
4. 全粒子をbatch開始時に固定した電場・電位の空間分布で追跡する。光電子だけはz-highで反射する。
5. 解決済み光電子の帰還先分布を使い、未帰還分の表面電荷を統計的に閉じる。
6. 太陽風吸収、光電子の放出・帰還をまとめて表面へcommitする。
7. 必要なら要素電位とz-high面平均電位を同じcommit後の表面電荷から出力する。

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

次の場合は結果を作らず停止します。

- 放出があるのに解決済み帰還電荷が0。
- 光電子がopen面から実際にescape。
- `soft_discard`、非有限値、または符号不整合。
- outer particle transferや`implicit_mean`との併用。

`max_step`まで残った粒子だけが統計補正の対象です。raw未帰還電荷、補正量、係数、未帰還率は
`charge_ledger.csv`へ別々に残します。

## z-high面を電位基準として読む

`output.write_files=true`、`output.write_potential_history=true`、`output.history_stride>0`、
`sim.use_box=true`なら、
`potential_history.csv`と同じbatchに`top_reference_history.csv`を出力します。

$$
\bar\phi_{\mathrm{top}}
=\frac{1}{N^2}\sum_{a,b}\phi(x_a,y_b,z_{\mathrm{high}}^-),\qquad
\phi_{\mathrm{rel}}(\mathbf r)=\phi(\mathbf r)-\bar\phi_{\mathrm{top}}.
$$

$N$は`sim.injection_face_phi_grid_n`です。標本点は全periodic cell上面のcell centerです。
二つの履歴を`batch`でjoinし、各要素の`potential_V`から同じbatchの
`potential_mean_V`を引きます。

```toml
[output]
history_stride = 1000
write_mesh_potential = true
write_potential_history = true
```

この面平均は無限遠電位でもプラズマ電位でもなく、太陽風流入へfeedbackしない診断値です。定数gaugeは消せますが、
box高さ、zero mode、有限画像打切りへの依存性は消しません。`potential_std_V`やmin/maxの幅が大きい場合は、
z-highを一つのreservoir面として解釈する近似が弱くなっています。

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
無限周期のnonzero modeと物理的zero modeが必要なら
[periodic2無限周期＋outer plasma構成](InfinitePeriodicOuterConfiguration.html)を使います。

## 受理前に確認する順序

1. `abs(neutral_return_weight_scale-1)`と`neutral_return_unresolved_fraction`が十分小さいことを確認する。
2. `dt`を下げるときも`max_step*dt`を維持または増加させ、`rays_per_batch`も増やして帰還先分布を収束させる。
3. `batch_duration`を少なくとも$T,T/2,T/4$と変え、電荷・相対電位・力の履歴を比較する。
4. z-highを上下し、`injection_face_phi_grid_n`を増やしてtop平均とばらつきを確認する。
5. image layerを$N,N+1,N+2$と増やす。
6. 太陽風のmacro粒子数と乱数seedに対する統計誤差を確認する。

完全反射は人工的な上端mirrorであり、自己整合sheathや準中性解ではありません。この構成から言えるのは、
「指定した局所太陽風fluxの下で、正味光電子電流を0としたときに表面内再分配が作る場」です。

旧outputに残る`sim.softening`や旧`[field]` tableは現行入力では削除済みです。過去の0p5x設定を使う場合も、
現行例を基に移植してください。

## scalar barrierを比較に使う場合

scalar barrier比較ではclosed PE設定を外し、次を使います。

```toml
[external_boundary.particles]
mode = "local_source"
inflow_model = "infinity_barrier"

[external_boundary.ordinary_open]
model = "potential_barrier"

[sim]
phi_infty = 0.0
```

これはface平均scalarで上流VDFを補正し、各open通過点のenergyでreflect/escapeを決めます。一つの連続した
sheath profileを往復するモデルではありません。式と制約は
[reservoir流入](ReservoirInjection.html)と[粒子のescapeとreturn](ParticleEscapeReturn.html)を参照してください。

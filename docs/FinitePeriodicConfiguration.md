title: periodic2有限画像構成

Lang: [日本語](FinitePeriodicConfiguration.md) | [English](FinitePeriodicConfiguration.en.md)

# periodic2有限画像構成

この構成では、x/y周期boxのsurface fieldをprimary cellと指定した画像層の和として定義します。外部plasmaの
空間profileは状態に持たず、流入・流出・光電子には必要に応じてscalar potentialによるreduced closureを組み合わせます。

## 場・流入・流出を有限画像モデルで揃える

| 処理 | 選択 |
| --- | --- |
| surface field | `field_bc_mode="periodic2"`、finite image sum、far correctionなし |
| source | `volume_seed`、`reservoir_face`、`photo_raycast` |
| reservoir補正 | なし、または `infinity_barrier` |
| open outflow | 無条件`escape`、または `potential_barrier` |
| photoelectron | 通常追跡、または`boltzmann_cutoff` reduced escape |
| outer Poisson/profile | なし |

有限画像和で使うkernel、Ewald far correctionとの違い、衝突判定での周期画像探索は、各成分の説明に
まとめています。[<sup>1</sup>](PeriodicElectrostatics.html)

## 画像層が場の物理範囲を決める

画像層$N$なら$(2N+1)^2$ cellのsourceを加えます。これは「無限周期和を近似する計算」に使えますが、設定した$N$の
結果そのものは有限画像modelです。非中性cellでは画像層を増やしてもpotential gaugeやz方向far boundaryが自動的に決まる
わけではありません。

`field_periodic_far_correction="none"`または互換`auto`は、Ewald/cached far operatorを追加しません。この設定を使うには、
場の目的量、粒子軌道、吸収位置、帯電分布についてimage-layer refinementを行い、必要な精度への収束を確認します。

## reservoir流入をface平均potentialで補正する

補正なしでは、設定した上流VDFをface上のflux-weighted分布としてsampleし、potentialによるcutoffや加減速は行いません。

`reservoir_potential_model="infinity_barrier"`を選ぶと、batch開始時に有限画像和で構成した場から注入開口の平均potential
$\bar\phi_f$を求め、`phi_infty`との差で

- 上流VDFの到達可能な法線速度を選ぶ。
- 同じpotential差でface法線速度へenergy mapする。

という2処理を行います。[<sup>2</sup>](ReservoirInjection.html)

これはface平均scalarだけを使うため、途中の$E(z)$、turning position、flight time、space chargeを持ちません。画像層を変えると
$\bar\phi_f$も変わり得るので、粒子fluxだけでなくface potentialのimage convergenceも確認します。

平均に使う同じ `N x N` sampleから電位の母標準偏差・最小・最大も集計するため、診断のための追加の
電位評価はありません。Maxwellian reservoirで局所電位の母標準偏差に対応するenergyが熱・法線driftの
特徴energyの10%を超えると、MPI rootは初回と最終batchに面平均近似の警告を出します。

## open面でescapeまたはscalar反射を選ぶ

`open_boundary_model="escape"`はopen面を横切った粒子を除去します。`potential_barrier`は通過点potentialと`phi_infty`を
比較し、法線energyが不足する粒子だけをreflectします。

reservoir側の`infinity_barrier`とoutflow側の`potential_barrier`は、同じenergy conventionを使えます。ただし、
前者はface平均で上流VDFを絞り込み、後者は個々の通過点で生成済みの粒子をreflectします。両者は別の処理であり、
一つの空間的なsheath profileを順方向と逆方向に写すものではありません。[<sup>3</sup>](ParticleEscapeReturn.html)

## 光電子をbox内追跡またはreduced cutoffで扱う

`photo_escape_model="none"`なら、ray hitから生成した光電子を通常粒子として追跡し、再吸収またはopen escapeまで進めます。
box外へ出た後のreturnは表現しません。

`boltzmann_cutoff`は放出元のlocal potentialからescape fractionを求め、tracked macro重みを減らします。非escape成分の軌道と
再吸収位置は生成しないreduced closureです。source reaction chargeにも減衰後の同じ重みを使います。

ray visibility、surface charge、tracked/reduced escapeは一つのライフサイクルとして扱います。[<sup>4</sup>](PhotoelectronEmission.html)

## 有限画像構成を設定する

```toml
[sim]
field_bc_mode = "periodic2"
field_solver = "fmm"
field_periodic_image_layers = 1
field_periodic_far_correction = "none"
use_box = true
bc_x_low = "periodic"
bc_x_high = "periodic"
bc_y_low = "periodic"
bc_y_high = "periodic"
bc_z_low = "open"
bc_z_high = "open"
open_boundary_model = "potential_barrier"
reservoir_potential_model = "infinity_barrier"
phi_infty = 0.0
injection_face_phi_grid_n = 5
```

`infinity_barrier`を使う場合は`phi_infty`と`injection_face_phi_grid_n`を明示します。光電子のreduced cutoffにも同じ
`phi_infty` gaugeが必要です。[<sup>5</sup>](Parameters.html)

実行例[`examples/periodic2_basic/beach.toml`](../examples/periodic2_basic/beach.toml)は、数値経路を確認する最小構成です。
本計算の画像層とsampling数は、後述の収束確認から決めます。

## 有限画像モデルを選ぶ範囲

適している用途:

- image layerを明示して行う小規模比較。
- free/periodic境界の回帰試験。
- scalar barrierを使った有限画像reservoir比較。
- Ewald/cached構成に対するnear-image reference。

適していない用途:

- layer refinementなしで無限周期surfaceを主張する計算。
- 自己整合な外部space chargeやsheath potential profileが必要な計算。
- outer flight time、delayed return、3D external orbitが重要な過渡。

## 画像層から帯電分布まで収束させる

1. `field_periodic_image_layers`を$N,N+1,N+2$と増やす。
2. surface/particle位置の$\phi,\mathbf E$だけでなく、reservoir flux、absorption/escape率、帯電分布を比較する。
3. `infinity_barrier`使用時は`injection_face_phi_grid_n`も増やす。
4. `photo_raycast`使用時は`rays_per_batch`、粒子追跡は`dt`も独立に収束させる。
5. 電荷収支と`discarded_unresolved`を確認する。

収束が遅い、または非中性cellのfar fieldを物理的に閉じる必要がある場合は
[periodic2無限周期＋outer plasma構成](InfinitePeriodicOuterConfiguration.html)へ移ります。

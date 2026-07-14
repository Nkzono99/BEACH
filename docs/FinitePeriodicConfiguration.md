title: periodic2有限画像構成

Lang: [日本語](FinitePeriodicConfiguration.md) | [English](FinitePeriodicConfiguration.en.md)

# periodic2有限画像構成

この構成では、x/y周期boxのsurface fieldをprimary cellと指定した画像層の和として定義します。外部plasmaの
空間profileは状態に持たず、流入・流出・光電子には必要に応じてscalar potentialによる簡略化モデルを組み合わせます。

## このページで採用する典型構成

| 処理 | 選択 |
| --- | --- |
| surface field | `field_bc_mode="periodic2"`、finite image sum、far correctionなし |
| source | `reservoir_face`と`photo_raycast`。必要なら初期粒子に`volume_seed` |
| reservoir補正 | `infinity_barrier`で上流VDFをfaceへenergy map |
| open outflow | `potential_barrier`で通過点ごとにreflect/escape判定 |
| photoelectron | `photo_escape_model="none"`でbox内を通常追跡 |
| outer Poisson/profile | なし |

有限画像和で使うkernel、Ewald far correctionとの違い、衝突判定での周期画像探索は、各成分の説明に
まとめています。[<sup>1</sup>](PeriodicElectrostatics.html)

## この構成を有効にするparameter

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

| parameter | 典型値 | この構成での意味 |
| --- | --- | --- |
| `field_bc_mode` | `"periodic2"` | x/y周期の画像をfield sourceへ加える |
| `field_periodic_image_layers` | `1`から収束確認 | 加える有限画像shellの範囲を決める |
| `field_periodic_far_correction` | `"none"` | Ewald/cached far operatorを加えず、有限画像modelに固定する |
| `reservoir_potential_model` | `"infinity_barrier"` | face平均電位で上流VDFの到達条件とface速度を決める |
| `open_boundary_model` | `"potential_barrier"` | 個々の通過点のenergyでreflect/escapeを決める |
| `phi_infty` | 物理gaugeに合わせる | reservoirとopen outflowが共有する基準電位 |
| `injection_face_phi_grid_n` | `5`から収束確認 | face平均とばらつきを求める各周期軸の標本数 |

このページの典型構成では、reservoir流入とopen流出の基準電位を`phi_infty`で揃えます。

光電子sourceには、次の組合せを使います。

```toml
[[particles.species]]
source_mode = "photo_raycast"
emit_current_density_a_m2 = 2.0e-4
rays_per_batch = 500
deposit_opposite_charge_on_emit = true
photo_escape_model = "none"
q_particle = -1.602176634e-19
m_particle = 9.1093837139e-31
temperature_ev = 1.5
normal_drift_speed = 1.0e5
```

`photo_escape_model="none"`は放出重みをreduced cutoffで減らさず、生成した光電子をbox内で通常粒子として追跡する指定です。
表面へ戻れば吸収され、open面へ達すれば`potential_barrier`でreflectまたはescapeします。
`deposit_opposite_charge_on_emit=true`は放出元へ逆符号の反作用電荷を加えます。各parameterの制約は
[入力パラメータリファレンス](Parameters.html)にまとめています。

実行例[`examples/periodic2_basic/beach.toml`](../examples/periodic2_basic/beach.toml)はfield経路だけを確認する最小構成です。
reservoirと光電子には上のparameter fragmentを追加します。本計算の画像層とsampling数は、後述の収束確認から決めます。

## 何セル先の周期画像まで場に含めるか

x/y周期では、primary cellの左右・前後に同じcellが繰り返されます。この構成は、その無限個のcopyをすべて足すのではなく、
primary cellから指定した層までのcopyだけをfield sourceに含めます。

| `field_periodic_image_layers` | 場に含めるcell |
| --- | --- |
| `0` | primary cellだけ（$1\times1$） |
| `1` | 周囲1層まで（$3\times3=9$ cells） |
| `2` | 周囲2層まで（$5\times5=25$ cells） |

したがって画像層$N$は、mesh分割の細かさではなく「何cell先の電荷まで相互作用に含めるか」を決めます。
指定範囲より外側の周期画像は、この構成の場には寄与しません。

画像層を増やしても電場、粒子の吸収・escape率、最終的な帯電分布がほとんど変わらなくなれば、目的量に対して
十分な範囲を含めたと判断できます。まだ変化するなら、その画像層では遠方cellの影響を切り捨てすぎています。

`field_periodic_far_correction="none"`と互換aliasの`"auto"`は、範囲外のcellをEwald和やcached operatorで補いません。
また、cellの正味電荷が0でない場合は、画像層を増やすだけでは無限周期系のpotential基準やz方向の遠方境界は決まりません。
そのような物理的な遠方境界モデルが必要なら、有限画像層を増やし続けるのではなく無限周期＋outer plasma構成を使います。

## reservoir流入をface平均potentialで補正する

`reservoir_potential_model="infinity_barrier"`は、batch開始時に有限画像和で構成した場から注入開口の平均potential
$\bar\phi_f$を求め、`phi_infty`との差で

- 上流VDFの到達可能な法線速度を選ぶ。
- 同じpotential差でface法線速度へenergy mapする。

という2処理を行います。[<sup>2</sup>](ReservoirInjection.html)

これはface平均scalarだけを使うため、途中の$E(z)$、turning position、flight time、space chargeを持ちません。画像層を変えると
$\bar\phi_f$も変わり得るので、粒子fluxだけでなくface potentialのimage convergenceも確認します。

平均に使う同じ `N x N` sampleから電位の母標準偏差・最小・最大も集計するため、診断のための追加の
電位評価はありません。Maxwellian reservoirで局所電位の母標準偏差に対応するenergyが熱・法線driftの
特徴energyの10%を超えると、MPI rootは初回と最終batchに面平均近似の警告を出します。

## open面で通過点ごとのenergyを判定する

`open_boundary_model="potential_barrier"`は通過点potentialと`phi_infty`を比較し、法線energyが不足する粒子だけを
reflectし、障壁を越えられる粒子をescapeとして除去します。

reservoir側の`infinity_barrier`とoutflow側の`potential_barrier`は、同じenergy conventionを使えます。ただし、
前者はface平均で上流VDFを絞り込み、後者は個々の通過点で生成済みの粒子をreflectします。両者は別の処理であり、
一つの空間的なsheath profileを順方向と逆方向に写すものではありません。[<sup>3</sup>](ParticleEscapeReturn.html)

## 光電子をbox内で通常追跡する

`photo_escape_model="none"`では、ray hitから生成した光電子を通常粒子として追跡し、再吸収またはopen面まで進めます。
このページの構成では`boltzmann_cutoff`を重ねません。open面に達した光電子にも他の粒子と同じ`potential_barrier`を適用します。
box外へescapeした後のreturn軌道は表現しません。

ray visibility、surface charge、tracked/reduced escapeは一つのライフサイクルとして扱います。[<sup>4</sup>](PhotoelectronEmission.html)

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

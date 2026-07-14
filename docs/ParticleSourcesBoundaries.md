title: 粒子源

Lang: [日本語](ParticleSourcesBoundaries.md) | [English](ParticleSourcesBoundaries.en.md)

# 粒子源

各batchの粒子追跡は、`[[particles.species]]`ごとに新しい粒子群を作るところから始まります。
生成後はsourceの種類にかかわらず同じ粒子状態へ入り、[Boris粒子更新](BorisPusher.html)と
[粒子の衝突・境界イベント](ParticleEvents.html)で追跡されます。

| `source_mode` | 粒子数を決める量 | 生成位置 | 主な用途 |
| --- | --- | --- | --- |
| `volume_seed` | `npcls_per_step` | `pos_low`〜`pos_high` | 初期粒子、軌道試験 |
| `reservoir_face` | 流入flux、面積、`batch_duration` | 指定したbox面 | 外部reservoirからの連続流入 |
| `photo_raycast` | 電流密度、投影面積、`batch_duration`、ray数 | rayが最初に命中した表面 | 光照射による表面放出 |

## snapshot更新後にそのbatchの粒子を作る

粒子源は、batch開始時に場と外部プラズマのsnapshotを更新した後で評価されます。この順序により、
reservoirの速度補正と光電子のreduced escape率は、前batchまでにcommitされた表面電荷を見ることができます。
生成された粒子は同じsnapshot中を進み、吸収、escape、`max_step`到達のいずれかまで追跡されます。

粒子が生成された瞬間に表面電荷を変えるのは、`photo_raycast`で放出元へ逆符号電荷を置く場合だけです。
その差分も追跡中の吸収電荷と同様にbatch末尾でcommitされ、同じbatchの場は変えません。
[計算モデルの全体像](Algorithms.html)に、snapshot更新からcommitまでの順序を示しています。

## `volume_seed`で指定個数の初期粒子を作る

`volume_seed`は各batchに`npcls_per_step`個を生成します。位置は直方体
`[pos_low, pos_high]`内の一様分布、速度は

$$
\mathbf v=\mathbf u+\sigma\mathbf Z,
\qquad
\sigma=\sqrt{\frac{k_\mathrm{B}T}{m}}
$$

というdrift付きMaxwell分布です。`thermal_speed`を指定した場合は温度から求めた$\sigma$より優先します。
標準正規変量は各成分で$6\sigma$に切られます。

このsourceの注入量は`npcls_per_step`で直接決まります。物理fluxに基づく連続流入には`reservoir_face`を使います。

## `reservoir_face`で連続流入fluxを作る

`reservoir_face`はbox外に与えた密度・温度・driftまたは速度gridから、指定面を内向きに横切るfluxを求めます。
生成する粒子数は物理fluxから計算し、法線速度はflux-weighted分布からsampleします。上流とfaceの間に電位差がある構成では、
到達可能な上流粒子の選別とface速度へのエネルギー写像を同じ電位差で行います。

式、速度grid、macro粒子端数、`infinity_barrier`、outer profileとの接続は
[reservoir注入](ReservoirInjection.html)で説明します。

## `photo_raycast`で照射面から粒子を放出する

`photo_raycast`はbox面上の照射開口からrayを発射し、box境界条件に従ってrayを進め、最初に命中した
in-box要素から粒子を放出します。放出速度は要素法線に対するflux-weighted Maxwell分布です。

放出後の光電子は共通の粒子状態へ入り、ほかの粒子と同じ場、Boris更新、mesh衝突、box境界を使います。
放出元電荷、reduced escape closure、outer sheathでのreturnまでを
[光電子の放出とライフサイクル](PhotoelectronEmission.html)にまとめています。

## 生成後は共通の粒子状態と電荷ledgerへ入る

生成後に保持する主な量は位置$\mathbf x$、速度$\mathbf v$、実粒子1個の電荷$q$と質量$m$、macro粒子重み$w$、
species IDです。tracked macro粒子1個の電荷は$q w$です。表面へ吸収された場合も、charge ledgerへの計上にはこの値を使います。

| batch結果 | 処理 |
| --- | --- |
| meshへ吸収 | 命中要素へ$q w$を堆積 |
| open面から無限遠へescape | 粒子を除去し、species別escapeへ計上 |
| outer領域からreturn | 同じ粒子をinterfaceへ戻し、残りstepを再積分 |
| `max_step`まで生存 | unresolvedとしてbatch末尾で破棄・計上 |

mesh衝突とbox境界の順序は[粒子の衝突・境界イベント](ParticleEvents.html)、外部領域のfield modelは
[外部プラズマモデル](OuterPlasmaModels.html)、粒子のescape/return写像は
[粒子のescapeとreturn](ParticleEscapeReturn.html)、表面への電荷commitは[表面電荷更新](SurfaceModels.html)で説明します。

## global生成量をMPIとrestartで保つ

`reservoir_face`のglobal生成数とmacro粒子端数はrootで一度だけ決めてrankへ分配します。そのためMPI world sizeを
変えても期待流入量は変わりません。端数はspeciesごとの`macro_residual`として
`macro_residuals.csv`へ保存され、再開時に復元されます。

`photo_raycast`の`rays_per_batch`も全rank合計です。各rayのmacro重みはglobal ray数で割り、放出・吸収・escapeの
charge ledgerはMPI all-reduce後のglobal値として出力します。期待流入量はworld sizeに依存しませんが、乱数列と
個々の粒子軌道はworld sizeによって変わり得ます。

## Code reference

- 粒子分布とraycast: [`bem_injection.f90`](../src/particles/bem_injection.f90)
- sourceごとのbatch生成とmacro残差: [`bem_app_config_runtime.f90`](../src/config/bem_app_config_runtime.f90)
- source入力の検証: [`bem_app_config_parser_validate.f90`](../src/config/app_config_parser/bem_app_config_parser_validate.f90)
- charge ledgerとbatch追跡: [`bem_simulator_loop.f90`](../src/runtime/simulator/bem_simulator_loop.f90)
- macro残差のcheckpoint: [`bem_restart.f90`](../src/runtime/bem_restart.f90)

title: 外部境界を構成する

Lang: [日本語](OuterPlasmaModels.md) | [English](OuterPlasmaModels.en.md)

# 外部境界を構成する

このページでは、外部 plasma 応答の場、z-high を横切る粒子、その他の open 面を設定します。
通常の入力では `[external_boundary]` 配下だけを編集します。BEACH が実装用の return・transfer・queue
設定へ変換するため、それらの内部名を組み合わせる必要はありません。

## 最初に 3 つの責務を選ぶ

```text
[external_boundary]
├── [external_boundary.field]          外部 plasma 応答の場
├── [external_boundary.particles]      z-high の粒子と reservoir 流入
└── [external_boundary.ordinary_open]  outer が所有しない open 面
```

`field` と `particles` は必須です。`ordinary_open` は任意で、省略時は `escape` です。

標準の自己整合 1D 外部シースは次の構成です。

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "absorbing_maxwellian"
debye_length = 0.2
thermal_voltage = 2.0

[external_boundary.particles]
mode = "same_batch"
field_evolution_timescale = 1.0
```

この入力から、BEACH は kinetic profile による流入写像、1D return、同一 batch 内の復帰を選びます。
`interface_z` は `sim.box_max[2]` から導出されるため入力しません。更新方式も公開設定ではなく、
通常は `explicit`、後述の ambient 線形応答または強 PE Zhao と tracked 光電子の組合せでは
内部の `implicit_mean` に決まります。
通常 open 面を単純に escape させる場合、`ordinary_open` は上のように省略します。

すべての open 面を単純に escape させ、流入補正や外部場も使わない場合は `[external_boundary]` 自体を省略します。
それ以外は、目的に最も近い行から構成を選んでください。

| 目的 | `field` | `particles.mode` | `inflow_model` |
| --- | --- | --- | --- |
| 外部場なし、または scalar barrier | `none` | `local_source` | `source_vdf` / `infinity_barrier` |
| 標準の自己整合 1D シース | `kinetic_1d` + `absorbing_maxwellian` | `same_batch` | `auto` |
| 弱い光電子放出の局所分布と平均帯電を分離する線形応答 | `kinetic_1d` + `ambient_linear_debye` | `same_batch` | `auto` |
| virtual cathode を含む強 PE の実測 CDF 準定常シース | `kinetic_1d` + `zhao_charge_driven` | `same_batch` | `auto` |
| Zhao の外部飛行遅延を含む過渡 | `kinetic_1d` + `zhao_charge_driven` | `zhao_queue` | `auto` |

通常流出に scalar barrier が必要な場合だけ次を追加します。

```toml
[external_boundary.ordinary_open]
model = "potential_barrier"
```

## 外部場を選ぶ

`external_boundary.field.model` は外部 plasma 応答の場だけを選びます。`kinetic_1d` は
interface 外を解きます。粒子を外部領域へ渡すかどうかは
`external_boundary.particles.mode` で別に選びます。

| `field.model` | 位置付け | 主な用途 |
| --- | --- | --- |
| `none` | 外部場なし | 通常 open、設定済み source VDF、scalar barrier |
| `kinetic_1d` | 対応する自己整合 1D シース | reservoir VDF、平均シース、流入、return を同じ profile で閉じる |

field-only も有効な構成です。

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "absorbing_maxwellian"
debye_length = 0.2
thermal_voltage = 10.0

[external_boundary.particles]
mode = "local_source"
inflow_model = "source_vdf"
```

弱い光電子放出を線形近似で扱う場合は、次の closure を選べます。

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "ambient_linear_debye"
debye_length = 0.2
thermal_voltage = 2.0

[external_boundary.particles]
mode = "same_batch"
field_evolution_timescale = 1.0
```

この closure は ambient plasma だけを
$\phi(z)=\lambda_D E_I\exp(-z/\lambda_D)$ として応答させます。
`photoelectron_density_model` は公開 TOML に書きません。facade が内部値を `none` に解決し、
`"none"` を含む key の明記を拒否します。
enabled な負電荷 `photo_raycast` species があると、BEACH は公開 key を追加せず、内部の
`implicit_mean` 更新を自動選択します。最初の 3D trace では batch-start の $k\ne0$ を固定し、
放出、interface 到達前の再吸収、ambient 吸収を element ごとに記録します。光電子の z-high crossing は
mean solve まで保留し、その実測外向き電流と tracked ambient 吸収電流から
解析 Maxwellian backward-Euler solver で総電荷、$\phi_I$、連続 escape / return 率を決めます。
`emit_current_density_a_m2` は tracked 3D 放出の weight を決めますが、独立な PE 平均 source として再利用しません。

mean solve 後の kinetic 1D profile で、各 full-weight crossing ray を 1 回だけ energy-resolved 追跡します。
return では outer flight time と横方向変位を計算し、mean $k=0$ と batch-start $k\ne0$ の合成場で終端まで追跡します。
return 後に z-high を再通過した ray も同じ profile mapper へ戻します。

解析 return 総電荷は、最終的に再吸収された ray sample の source leg（放出元 countercharge の一時中和）と
実 hit destination leg へ
同じ `return_weight_scale` で正規化配分します。transaction は零和を保ち、pending の放出 countercharge を
保持したまま destination deposit を一度だけ commit します。
平均総電荷と引戻し電位は同じ batch 内で更新し、sampled 局所分布は次の batch の $k\ne0$ へ反映するため、
離散 ray 分類を mean solve へ反復して戻しません。

mean solver の失敗、解析 return が正なのに再吸収 sample がない場合、電荷不一致、
または許可された trace 内で吸収・escape へ終端しない ray は fail closed で停止します。deferred PE ray は
準定常 shadow なので、flight time と frozen-field ratio は診断しますが、ratio 超過では停止しません。
通常の `same_batch` 粒子と ambient species は上限超過で fail closed のままです。光電子は
1D 平均密度と outer space charge には入らず、nonlinear photoelectron density と outer cloud inventory は未対応です。
ledger の escape / reabsorption 電荷は解析 closure の総量、整数 count は ray sample の最終分類を表します。
標準出力は `transaction_residual_C`、`mean_solver_iterations`、`sample_escape_fraction`、
`return_weight_scale` を記録します。

この構成は `outer_update_stride=1`、正の `batch_duration`、
`deposit_opposite_charge_on_emit=true`、ちょうど 1 つの負電荷 `photo_raycast` species、
`photo_raycast.normal_drift_speed=0` を要求します。解析 Maxwellian escape 率は放出法線 drift を含みません。
この経路は mesh mode に固有の制約を追加しません。

この準定常 shadow sampling は UV turn-on の遅延 return current を時間発展させません。その用途には
`ambient_linear_debye` に対応する delayed inventory / queue を別途設計する必要があります。
virtual cathode、space-charge-limited / inverse sheath、trapped population を表せないため、
光電子空間電荷が ambient 線形応答に対して十分小さい場合だけ使います。
完全な例は
[`periodic2_ambient_linear_photo_outer.toml`](../examples/periodic2_ambient_linear_photo_outer.toml)
です。

### 強 PE では実測 CDF と Zhao profile を接続する

photoemitting `zhao_charge_driven + same_batch`では、BEACHは別の`implicit_mean`経路を自動選択します。
ambient経路の解析Maxwellian backward Eulerや共通`return_weight_scale`は使いません。
実測したz-high crossingの法線energy ${\cal E}_n=m_{pe}v_z^2/2$をcharge weight付きでexact group化し、
Type Aのvirtual cathodeを含むprofile全体のbarrier $B(Q)$と

$$
Q=Q_{\rm base}+C_{\rm esc}[B(Q)]
$$

を同時に満たす一般化電荷根を解きます。同一energy groupがbarrier上に来た場合だけfractional escape weightを使います。

surfaceの`emit_current_density_a_m2`は放出rayとmacro weightを決めます。Zhaoへ渡す各batchのsource amplitudeは、
surface設定値ではなく、実際にinterfaceへ到達した総chargeを$A\Delta t$で割った実測電流から解きます。
したがってrough surfaceでの再吸収やinterface transmissionがsource normalizationへ反映されます。

solverはcommit済みA/B/C rootを始点として
$\mathbf G(\mathbf y;E_I(\lambda),n_{pe,0}(\lambda))=\mathbf 0$をadaptive pseudo-arclengthで追跡します。
Type A/Bは$E_I>0$、Type Cは$E_I<0$を要求します。parameter foldではfail closedし、
$\lambda=1$のevent correctorでtargetへ着地します。
最終sourceを固定した候補pathでは、各adaptive accepted pointでbarrierのprescribed chargeに対する接線勾配を検査し、
負勾配や逆向きsecantでは停止します。共通rootから$Q_0$と$Q_M$へ進む2本のpathで全charge区間を検査した後、
単調なorder predicateのfirst-true indexを二分探索するため、pure root選択は$O(\log M)$です。
最後に実測CDFの電荷残差も再検査します。local continuation guardは
accepted point間をinterval arithmeticで包含する数理的証明ではありません。branch変更、rank低下、fold、
barrier勾配反転、order predicateのbracket不整合、marginal energyのbracket不成立、またはfrozen cohortの
0.25 trust region違反では、別branchや旧profileへfallbackせず停止します。この内部solverのために追加する
TOML keyはありません。

実測するのはinterface source amplitudeと法線energy CDFだけです。Zhao density populationは解析closureのままで、
任意のtrapped VDF、衝突、磁化、PICは解きません。MPIでは全interface sampleをrootへ集めるため、
root memoryとgather通信量はcrossing ray数に比例します。数式、trust region、ray終端条件は
[外部シース: kinetic 1D](KineticOuterPlasma.html#強-pe-を実測-interface-energy-cdf-で閉じる)を参照してください。

## 接続面の粒子を選ぶ

| `particles.mode` | z-high を出た粒子 | 用途 |
| --- | --- | --- |
| `local_source` | `ordinary_open` で処理し、外部軌道を保持しない | 外部場なし、field-only、scalar 流入 |
| `same_batch` | kinetic 1D profile による return を計算 | 定常・準定常の tracked return |
| `zhao_queue` | Zhao 光電子の return/escape event を due batch まで保持 | 強 UV 立上がりなどの遅延電流 |

`particles.mode` が選ぶのは z-high 粒子の外部輸送だけです。reservoir からの流入 VDF は
後述の `inflow_model` で独立に選びます。

`same_batch` は `field.model="kinetic_1d"` と組み合わせ、離散 kinetic 1D profile 上で
return / escape を判定します。上記の `implicit_mean` 光電子は interface crossing を一度保留して
平均場を更新してから、同じ profile mapper で粒子ごとに判定します。ambient と強 PE Zhao は平均場の
解法とray weightが異なります。他 species の same-batch return は変わりません。

`zhao_queue` は汎用的な遅延輸送ではありません。`kinetic_1d` と
`kinetic_closure="zhao_charge_driven"` の組合せだけで使用できます。

## reservoir 流入を選ぶ

`external_boundary.particles.inflow_model` は `reservoir_face` に渡す上流 VDF の扱いです。

| `inflow_model` | 動作 |
| --- | --- |
| `auto` | 1D tracked model では field profile、それ以外では設定した source VDF を使用 |
| `source_vdf` | 設定した VDF を face 上の分布として使用 |
| `infinity_barrier` | face 平均電位と `sim.phi_infty` から到達速度を補正 |

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "zhao_charge_driven"
# split-interface 診断の基準量。Zhao profile の scale ではない
debye_length = 10.5
thermal_voltage = 10.0

[external_boundary.particles]
mode = "same_batch"
steady_start_mode = "zhao_floating"
steady_start_mesh_id = 1
field_evolution_timescale = 1.0
```

Zhao profile の物理的な長さ・電位 scale は、光電子ありでは $T_{pe}$ と基準密度、
光電子なしでは ambient $T_e$ と $n_\infty$ から導出します。`debye_length` と
`thermal_voltage` は Zhao root/profile を変えませんが、split-interface の gap、
lateral field、local-charge 診断の基準量として現在は必須です。
強 PE `same_batch`は一様seed chargeを置ける、periodic cell全面を覆う水平な
`mesh.mode="template"` planeを`steady_start_mesh_id`で選ぶ必要があります。
`mode="zhao_queue"` へ切り替える場合は、`steady_start_mode`と`steady_start_mesh_id`を削除し、
正の `sim.batch_duration` と光電子 sourceを指定します。
update stride は queue が内部で 1 に固定するため追加しません。

`kinetic_1d` の `same_batch` と `zhao_queue` では、同じ 1D profile が流入を所有します。
この場合、`inflow_model` は `auto` だけを許し、`infinity_barrier` との二重補正を拒否します。

## 通常 open 面を選ぶ

`external_boundary.ordinary_open.model` は outer model が所有しない open 面へ適用します。

| `ordinary_open.model` | 動作 |
| --- | --- |
| `escape` | 粒子を永久流出として除去 |
| `potential_barrier` | 局所電位と法線運動 energy から反射または escape を判定 |

`potential_barrier` は reservoir 流入の `infinity_barrier` と別の責務です。一方だけでも、両方でも使えます。
tracked outer model が z-high を所有する場合も、`ordinary_open` はその他の open 面に残ります。

## 設定エラーとして扱う組合せ

BEACH は矛盾する値を補正したり、別モデルへ silent fallback したりしません。

- `field.model="none"` と `same_batch` / `zhao_queue`
- kinetic の tracked 1D profile と local inflow 補正
- Zhao 以外の closure と `zhao_queue`
- `zhao_queue` と手動 branch、steady start
- photoemitting Zhao の `same_batch` で `zhao_floating` steady start がない構成
- 選んだ field / particle mode で効果のない key。既定値を明示しただけでも no-op として拒否
- 対応しない磁場、species、periodic2、zero-mode、時間スケール
- `[external_boundary]` と旧 `[outer_plasma]` / `[coupling]` の混在

物理入力や数値 guard は自動推定しません。species、Debye 長、温度、field timescale、
periodic2 backend などは使用する model に応じて明示し、矛盾時はエラーにします。
`zhao_queue` の update stride は内部で 1 に固定するため、公開入力へ書きません。

## 解決された契約を確認する

`summary.txt` には authoring 入力ではなく、実行に使った解決済みの流入、通常 open 面、interface transport、
粒子 mode を記録します。

旧 `[sim]` / `[outer_plasma]` / `[coupling]` からの移行は
[外部境界設定の移行](BoundaryConfigurationMigration.html)を使ってください。個々のモデルの物理と妥当性確認は
[外部シース: kinetic 1D](KineticOuterPlasma.html)、
[open 境界・escape・return](ParticleEscapeReturn.html)に分けています。

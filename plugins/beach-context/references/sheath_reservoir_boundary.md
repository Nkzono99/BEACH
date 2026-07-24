# 外部シースと reservoir 粒子境界

この reference は、BEACH の外部場、reservoir 流入、z-high particle transfer、通常 open 面を
物理的な責務から選ぶためのものです。公開入力では `[external_boundary]` だけを編集します。

## 1. 公開設定の責務

```text
[external_boundary]
├── [external_boundary.field]          外部 plasma 応答の場
├── [external_boundary.particles]      z-high particle transfer と reservoir 流入
└── [external_boundary.ordinary_open]  outer が所有しない open 面
```

| 責務 | 選択肢 |
| --- | --- |
| `field.model` | `none` / `kinetic_1d` |
| `particles.mode` | `local_source` / `same_batch` / `zhao_queue` |
| `particles.inflow_model` | `auto` / `source_vdf` / `infinity_barrier` |
| `ordinary_open.model` | `escape` / `potential_barrier` |

`potential_barrier` は流出、`infinity_barrier` は流入を担当します。名前が似ていますが独立した設定です。

## 2. 目的から構成を選ぶ

| 目的 | field | particles | inflow |
| --- | --- | --- | --- |
| 外部場なし、設定した VDF をそのまま注入 | `none` | `local_source` | `source_vdf` |
| 外部場なし、無限遠から face までの scalar 障壁を適用 | `none` | `local_source` | `infinity_barrier` |
| 標準の自己整合 1D sheath | `kinetic_1d + absorbing_maxwellian` | `same_batch` | `auto` |
| 蓄積電荷で閉じる Zhao sheath | `kinetic_1d + zhao_charge_driven` | `same_batch` | `auto` |
| Zhao outer flight を後続 batch へ遅延 | `kinetic_1d + zhao_charge_driven` | `zhao_queue` | `auto` |

自己整合な外部シースには `kinetic_1d` を使います。

## 3. `source_vdf` と `infinity_barrier`

`source_vdf` は設定した `reservoir_face` VDF を注入面上の分布として使います。

`infinity_barrier` は開口内の平均電位 $\phi_f$ と `sim.phi_infty` から

$$
B=\frac{2q(\phi_f-\phi_\infty)}{m}
$$

を作ります。$B>0$ では $v_{n,\infty}\ge\sqrt B$ の粒子だけが face へ到達し、

$$
v_{n,f}^2=v_{n,\infty}^2-B
$$

で到達速度を求めます。$B<0$ では全流入粒子が加速します。このモデルは途中の $E(z)$、
turning point、flight time、space charge を解きません。

## 4. tracked `kinetic_1d` は流入と return を所有する

`kinetic_1d + same_batch` では、無限遠 electron / ion VDF から解いた同じ離散 profile を
reservoir 流入と z-high return / escape の両方に使います。

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "absorbing_maxwellian"
debye_length = 0.2
thermal_voltage = 2.0

[external_boundary.particles]
mode = "same_batch"
inflow_model = "auto"
field_evolution_timescale = 1.0
```

この構成へ `infinity_barrier` を重ねると同じエネルギー障壁を二重に補正するためエラーです。
無限遠 gauge は内部で `phi(infinity)=0` に固定され、公開・互換入力では変更できません。

Type A のような非単調 profile では endpoint の電位差だけで判定せず、全 profile 上の最初の
energy barrier を流入・流出の両方向で走査します。

## 5. Zhao charge-driven closure

蓄積電荷で閉じる Zhao sheath は `kinetic_1d` の closure として選びます。

```toml
[external_boundary.field]
model = "kinetic_1d"
kinetic_closure = "zhao_charge_driven"
debye_length = 10.5
thermal_voltage = 10.0

[external_boundary.particles]
mode = "same_batch"
inflow_model = "auto"
field_evolution_timescale = 2.0e-5
```

surface zero mode が与える interface field を境界条件として Zhao A/B/C profile を更新します。
`same_batch` は outer flight を状態写像へ使いますが global time へ加えません。
遅延電流が必要なら `mode="zhao_queue"` を使い、有限 control volume 内の tracked inventory と
Zhao photoelectron column を predictor/corrector で閉じます。

Zhao closure は background electron / ion species、必要に応じた `photo_raycast` species、
温度・密度・drift、`sim.sheath_alpha_deg`、`sim.sheath_photoelectron_ref_density_cm3` を参照します。
`debye_length` と `thermal_voltage` は split-interface 診断の reference scale であり、
Zhao root/profile の物理 scale ではありません。

## 6. 通常 open 面

`ordinary_open.model="escape"` は粒子を永久流出として除去します。
`potential_barrier` は境界通過点の電位 $\phi_b$ と法線速度 $v_n$ から

$$
U_b=q(\phi_\infty-\phi_b)
$$

を計算し、$U_b>0$ かつ $m v_n^2/2<U_b$ なら法線速度を反転します。それ以外は escape です。
outer が z-high を所有しても、その他の open 面にはこの規則が残ります。

## 7. fail-closed 条件

- `field.model="none"` と `same_batch` / `zhao_queue`
- tracked `kinetic_1d` と `inflow_model!="auto"`
- Zhao 以外の closure と `zhao_queue`
- unsupported な species、磁場、periodic2、zero mode、時間 scale
- 選択した model / mode で効果のない key
- `[external_boundary]` と旧 `[outer_plasma]` / `[coupling]` の混在

BEACH は矛盾した入力を別モデルへ変換せず停止します。旧モデルから移す場合も数値的・物理的な等価性を
仮定せず、[外部境界設定の移行](../../../docs/BoundaryConfigurationMigration.md)に従って計算目的を選び直します。

title: 粒子のescapeと局所return

Lang: [日本語](ParticleEscapeReturn.md) | [English](ParticleEscapeReturn.en.md)

# 粒子のescapeと局所return

BEACH は box 外の plasma を時間発展させません。open 面では
`external_boundary.ordinary_open.model="escape" | "potential_barrier"` を局所的に適用します。この処理は
source に依存せず、`reservoir_face`、`photo_raycast`、`volume_seed` に共通です。closed PE だけは
species 固有の z-high 反射と charge closure を組み合わせます。

## 1. `escape`: open面で粒子を除去する

```toml
[external_boundary.field]
model = "none"

[external_boundary.particles]
mode = "local_source"

[external_boundary.ordinary_open]
model = "escape"
```

粒子は境界通過位置で除去されます。macro電荷$qw$はspecies別の`escaped_to_infinity`へ計上され、
表面電荷`q_elem`は変更されません。

複数面を同時に横切る場合の規則と、reflect/periodic面を処理した後の残りstepの再積分は
[粒子の衝突・境界イベント](ParticleEvents.html)で説明します。

## 2. `potential_barrier`: scalar障壁で反射を判定する

`potential_barrier` は通過点の scalar 電位だけを使う reduced model です。

```toml
[sim]
phi_infty = 0.0

[external_boundary.field]
model = "none"

[external_boundary.particles]
mode = "local_source"

[external_boundary.ordinary_open]
model = "potential_barrier"
```

open面の通過点電位を$\phi_b$、外向き法線速度を$v_n>0$とすると、無限遠へ進むための障壁は

$$
U_b=q(\phi_\infty-\phi_b)
$$

です。

$$
\frac12 m v_n^2<U_b\quad\text{かつ}\quad U_b>0
$$

なら法線速度を反転し、残りstepを追跡します。それ以外はescapeです。接線速度は変えません。

open 面外側の電場、turning 位置、flight time、空間電荷は保持しません。複数 open 面を同時に横切る corner は
`ambiguous_open_corner` として停止します。

通過点電位は粒子運動と同じbatch開始時固定場で評価され、`sim.e0`の局所電位も含みます。
一様電場には有限な無限遠電位がないため、`sim.e0`と併用する場合は`phi_infty`を有効な
reservoir基準電位として整合させてください。

## 閉じた光電子のreturn

光電子をbox内で閉じる場合は、該当する負の`photo_raycast` speciesに
`z_high_boundary="reflect"`と`surface_charge_closure="neutral_return"`を指定します。
反射した光電子は同じ追跡粒子としてbox内に残り、batch末尾で放出counterchargeを補正します。
設定と収支条件は[光電子放出](PhotoelectronEmission.html)を参照してください。

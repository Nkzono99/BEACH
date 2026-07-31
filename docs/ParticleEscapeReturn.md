title: 粒子のescapeと局所return

Lang: [日本語](ParticleEscapeReturn.md) | [English](ParticleEscapeReturn.en.md)

# 粒子のescapeと局所return

BEACH は box 外の plasma を時間発展させません。`[particle_boundary]` で非周期面を `open`、`reflect`、
`redistributed_reflect` のいずれかにし、open 面へ `ordinary_open_model="escape" | "potential_barrier"` を適用します。この処理は
生成方法に依存せず、boundary reservoir流入、`plane_source`、`photo_raycast`、`volume_seed`、
deprecatedな`reservoir_face`に共通です。周期性は
`domain.periodic_axes` で定め、粒子境界へ `periodic` は指定しません。
完全修飾キーは `particle_boundary.ordinary_open_model` です。

## 1. `escape`: open面で粒子を除去する

```toml
[particle_boundary]
x_low = "open"
x_high = "open"
y_low = "open"
y_high = "open"
z_low = "open"
z_high = "open"
ordinary_open_model = "escape"
```

粒子は境界通過位置で除去されます。macro電荷$qw$はspecies別の`escaped_to_infinity`へ計上され、
表面電荷`q_elem`は変更されません。

複数面を同時に横切る場合の規則と、reflect/periodic面を処理した後の残りstepの再積分は
[粒子の衝突・境界イベント](ParticleEvents.html)で説明します。

## 2. `potential_barrier`: scalar障壁で反射を判定する

`potential_barrier` は通過点の scalar 電位だけを使う reduced model です。

```toml
[particle_boundary]
ordinary_open_model = "potential_barrier"

[reservoir]
phi_infty = 0.0
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
一様電場には有限な無限遠電位がないため、`sim.e0` と併用する場合は `reservoir.phi_infty` を有効な
reservoir 基準電位として整合させてください。

## 閉じた光電子の return

closed PE では、負の `photo_raycast` species の `[particles.species.boundary]` で `inject_face` と同じ面を
`reflect` または `redistributed_reflect` にし、`surface_charge_closure="neutral_return"` を指定します。
どちらも法線速度を反転して接線速度を維持しますが、return位置の扱いが異なります。

| action | return位置 |
| --- | --- |
| `reflect` | 境界eventの接線位置を維持 |
| `redistributed_reflect` | 単一面では面内2軸をbox spanの両端guardを除く範囲から一様再配置 |

edge / cornerの同時eventでは、`redistributed_reflect`はevent maskに含まれない軸だけを再配置します。
詳細は[粒子の衝突・境界イベント](ParticleEvents.html)を参照してください。global 面を open のままにすれば、
`inherit` の ambient species は通常の open 契約に従います。設定と収支条件は
[光電子放出](PhotoelectronEmission.html)にまとめています。

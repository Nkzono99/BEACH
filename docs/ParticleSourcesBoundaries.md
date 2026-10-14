title: 粒子をどこから入れるか

Lang: [日本語](ParticleSourcesBoundaries.md) | [English](ParticleSourcesBoundaries.en.md)

# 粒子をどこから入れるか

このページは、粒子を生成する位置と物理的な供給方法から、使う粒子源を選ぶための案内です。
領域内に粒子を置く方法は `source_mode`、領域外から box 境界を横切らせる方法は
`[particles.species.boundary_inflow]` が担当します。

> **選び方の要点:** 指定個数を置くなら `volume_seed`、外部 plasma reservoir から入れるなら
> `boundary_inflow`、照射面から光電子を放出するなら `photo_raycast` を使います。

読了後には、粒子源を一つ選び、詳細設定を読むべき専用ページまたはリファレンスへ進めます。

## 生成位置から選ぶ

| 目的 | 選ぶ設定 | 粒子数を決める量 | 生成位置 |
| --- | --- | --- | --- |
| 初期粒子や軌道試験を指定個数で作る | `source_mode="volume_seed"` | `npcls_per_step` | `pos_low` と `pos_high` の間 |
| box 外の plasma から流入させる | `boundary_inflow` | reservoir 流束 × box 面積 × `batch_duration` | 非周期の box 面 |
| 光が当たった表面から放出する | `source_mode="photo_raycast"` | 電流密度、投影面積、`batch_duration` | ray が最初に命中した表面 |

境界流入だけを使う species は、物理パラメータと `boundary_inflow` を指定します。
`source_mode` と `npcls_per_step` は省略し、体積生成などを使う場合に明示します。
省略時は体積生成を行いません。既定値と併用制約は[入力パラメータ](Parameters.html)を参照してください。

## 指定個数を置く: `volume_seed`

`volume_seed` は、各 batch に `npcls_per_step` 個の粒子を作ります。位置は `pos_low` と `pos_high` が
囲む直方体から、速度は設定したドリフト付き Maxwell 分布から標本化します。

物理的な面流束から個数を求める source ではありません。そのため `batch_duration=0` でも実行できますが、
その batch に物理秒を割り当てることはできません。公式の
[10 分チュートリアル](Tutorial.html)は、この方法で batch 間の表面帯電 feedback を見せます。

## 外部 reservoir から入れる: `boundary_inflow`

`boundary_inflow` は、非周期の box 面全体を横切る外部 plasma 流入です。選択した面は、粒子が外向きに
到達したときの作用を別途 `open` にします。通常の開放面では外向き粒子が脱出し、
[外部シースでの帰還](ParticleEscapeReturn.html#外部シースの障壁)は Zhao の kinetic 障壁が判定します。

```toml
[[particles.species]]
# density, temperature, charge, mass, and macro-particle weight

[particles.species.boundary_inflow]
z_high = "reservoir"
```

Maxwell 分布、面ごとの流束と端数、無限遠の plasma 電位から境界までの補正は
[境界から粒子を流入させる](ReservoirInjection.html)で説明します。このモデルは box 外の軌道や
自己無撞着な外部 sheath を解きません。外部 sheath 応答そのものが必要なら
[Zhao closure](ZhaoStationaryClosure.html)を検討してください。

## 照射面から放出する: `photo_raycast`

`photo_raycast` は box 面上の照射開口から ray を飛ばし、最初に命中した三角形から粒子を放出します。
放出元へ残す逆符号電荷、再吸収、closed photoelectron の return は
[光電子の放出とライフサイクル](PhotoelectronEmission.html)で扱います。

## 生成後の追跡は共通

どの方法で生成した粒子も、同じ粒子追跡へ入ります。

| 結果 | BEACH の処理 |
| --- | --- |
| 三角形へ吸収 | 命中要素の電荷差分へ加える |
| open 面から脱出 | 粒子を除去し、species 別 escape として数える |
| 反射または周期境界を通過 | 同じ粒子の残り時間を追跡する |
| `sim.max_step` まで生存 | 未解決 survivor として数え、batch 末尾で破棄する |

box 境界を出る粒子の扱いは[粒子の escape と return](ParticleEscapeReturn.html)、粒子 event の
判定順は[粒子の衝突・境界 event](ParticleEvents.html)を参照してください。

## 次に読むページ

- 形状、粒子源、境界、solver を一つのケースにまとめる: [ケースを設計する](ConfigurationRecipes.html)
- 境界流束と `phi_infty` を設定する: [境界から粒子を流入させる](ReservoirInjection.html)
- 表面衝突後の電荷を扱う: [表面はどう帯電するか](SurfaceModels.html)
- 完全な key と制約を検索する: [入力パラメータ](Parameters.html)

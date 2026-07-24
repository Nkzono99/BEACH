# ADR 0009: 線形 Debye 外部場と静的 sheath 流入補正の hard removal

- Status: accepted
- Date: 2026-07-24
- Supersedes: ADR 0008 の field / inflow 語彙と、該当 raw selector の互換維持

> この ADR が維持対象とした `unified_linear_response` は
> [ADR 0010](0010-remove-unified-linear-response.md)で削除した。以下の該当箇所は当時の移行判断を示す履歴である。

## Context

ADR 0007 と ADR 0008 は、外部境界を field、z-high particle lifecycle、reservoir 流入、
通常 open 面へ分離した。この責務分離は有効だが、公開語彙には物理的な閉じ方が異なる
`linear_debye`、`kinetic_1d`、`unified_linear_response` と、field から独立した
`legacy_sheath` source 補正が同時に残った。

`linear_debye` は解析的な指数 1D profile と専用 return、`legacy_sheath` は
surface charge から独立した静的 source-VDF 補正である。一方、標準の `kinetic_1d` は
species VDF と surface charge から nonlinear profile を更新し、tracked mode では同じ profile が
流入と return を所有する。名前を alias しても同じ計算にはならない。

線形 return 専用の光電子 histogram、任意の `infinity_potential`、静的 sheath runtime wrapper まで
互換維持すると、公開 schema、Fortran / Python parser、restart、出力、テストへ専用分岐が残り、
利用者が設定から物理 ownership を判断しにくい。

## Decision

公開設定を次へ限定する。

| 責務 | 値 |
| --- | --- |
| `external_boundary.field.model` | `none` / `kinetic_1d` / `unified_linear_response` |
| `external_boundary.particles.mode` | `local_source` / `same_batch` / `zhao_queue` |
| `external_boundary.particles.inflow_model` | `auto` / `source_vdf` / `infinity_barrier` |
| `external_boundary.ordinary_open.model` | `escape` / `potential_barrier` |

次を deprecation alias なしで削除する。

- public / raw `linear_debye` field model と解析 1D return
- public `legacy_sheath` / `legacy_sheath_model` と raw `sim.sheath_injection_model`
- static source-plan override と専用 runtime wrapper
- 線形 return 専用の photoelectron histogram、checkpoint、summary 項目
- public / raw `infinity_potential` key

`kinetic_1d` と `unified_linear_response` の無限遠 gauge は内部で
`phi(infinity)=0` に固定する。`summary.txt` や outer-state checkpoint に残る
無限遠電位は状態診断であり、authoring key ではない。

ADR 0007 の runtime responsibility、同一 local step の残り時間処理、通常 open 面、
`infinity_barrier` と `potential_barrier` の独立性は維持する。新しい strategy registry、
任意面 interface、汎用 preset 層は追加しない。

Zhao の density、charge、root、branch に関する共有数理 core は残す。
`kinetic_1d + zhao_charge_driven` と `steady_start_mode="zhao_floating"` が同じ数理を使うため、
静的 source 補正の削除を Zhao core 全体の削除へ広げない。

## Migration

削除値を別モデルへ自動変換しない。次は目的から選ぶ候補であり、数値的・物理的に等価な置換ではない。

| 以前の目的 | 検討する現行構成 | 非等価である理由 |
| --- | --- | --- |
| 解析的な 1D field / same-batch return | `kinetic_1d` + `local_source` / `same_batch` | species VDF と nonlinear Poisson profile、zero gauge を使う |
| rough surface を含む線形 3D screening | `unified_linear_response` | 3D field response であり、1D sheath を解かない |
| source VDF を補正しない流入 | `field.model="none"` + `inflow_model="source_vdf"` | 設定した face VDF をそのまま使う |
| 明示電位による scalar 流入障壁 | `field.model="none"` + `inflow_model="infinity_barrier"` | 静的な zero-current closure ではない |
| surface charge で閉じる Zhao sheath | `kinetic_1d + zhao_charge_driven` | profile を batch ごとに更新し、流入と return を同じ profile で所有する |
| 線形 return 専用 histogram | 直接の置換なし | 必要なら particle / event 出力から独立に解析する |

BEACH 1.5 / 1.6 で作成した設定のうち削除値を含むものは、現行 parser でエラーになる。
該当モデルの checkpoint は runtime contract、model fingerprint、保存 state が異なるため再開しない。
移行後の構成は初期状態から実行する。過去出力を解析資料として保持することと、
その checkpoint から計算を継続することは区別する。

## Consequences

標準 1D sheath は `kinetic_1d`、高度な rough-surface 線形 screening は
`unified_linear_response`、外部場なしは `none` となる。流入は `auto`、
設定済み VDF、scalar barrier の 3 択になり、field model と静的 source closure の
競合表を公開設定から除去できる。

parser、JSON Schema、examples、docs、model fingerprint、restart、summary、MPI state、
Fortran / Python tests は同じ hard-removal boundary へ更新する。解析 linear modelを
検証 oracle として残す場合も公開 TOML から到達できない内部名・内部テストに限定する。

## Validation

1. Fortran parser、Python validator、schema の field / inflow enum が一致すること。
2. 削除値と raw selector が未知値または removed-value error で fail closed になること。
3. `kinetic_1d` の field-only、same-batch return、Zhao queue と、unified 3D orbit の coverage を維持すること。
4. generic return、restart、output、MPI test を削除モデルの fixture に依存させないこと。
5. shared Zhao core の利用箇所を監査し、kinetic Zhao と stationary warm start を維持すること。
6. canonical docs / examples / plugin reference が現行語彙だけを案内すること。

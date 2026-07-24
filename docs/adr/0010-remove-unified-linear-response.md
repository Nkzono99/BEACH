# ADR 0010: `unified_linear_response` の hard removal

- Status: accepted
- Date: 2026-07-24
- Supersedes: [ADR 0002](0002-unified-periodic-outer-domain.md)

## Context

`unified_linear_response` は、rough surface の高さ範囲から遠方 plasma までを一つの線形応答領域として扱うために
導入した。平面平均では accessible fraction を掛けた Debye--Hückel 応答、非零 Fourier mode では
screened reflection/transmission tail を用い、`same_batch` では合成 3D 静電場中の明示的 outer orbit も提供した。

しかし、この機能は設定、snapshot、MPI、checkpoint、fingerprint、出力、schema、例、検証へ広く結合していた。
その一方で、公開例は平面形状に限られ、主用途である rough surface に対する物理的妥当性、数値収束、適用限界を
十分に検証できていなかった。`kinetic_1d` と異なる近似を使うにもかかわらず、同じ field model 選択肢へ並べたため、
設定上も「より高精度な sheath model」と誤解しやすかった。

## Decision

`unified_linear_response` を互換 alias や feature flag を残さず削除する。

削除対象は、次を一つの機能境界として扱う。

- accessible-fraction 付き zero-mode 線形応答
- screened nonzero-mode tail
- unified field 用の surface-to-far snapshot 更新
- `electrostatic_3d_explicit_orbit` による 3D outer return/transfer
- unified 専用の設定、診断、出力、checkpoint state、例、テスト、利用者向け文書

外部場として残す自己整合 model は `kinetic_1d` とする。Zhao closure、`source_vdf` /
`infinity_barrier`、`ordinary_open.escape` / `potential_barrier` は独立した責務として維持する。

## Consequences

これは breaking change である。`unified_linear_response` を指定した旧設定は拒否され、別 model へ自動変換しない。
unified state を含む checkpoint は再開できない。roughness と plasma 応答が重なり、妥当な split window を置けない
ケースには、現行実装内に物理的に等価な代替はない。

削除により、標準設定は `none` または `kinetic_1d` に絞られる。`same_batch` は kinetic 1D profile 上の
return/transfer、`zhao_queue` は Zhao の persistent queue だけを表すため、field と粒子 lifecycle の対応が
明確になる。

## Reintroduction criteria

rough-surface plasma response が必要になった場合は、旧実装をそのまま復元せず、少なくとも次を満たす設計から
再開する。

1. 対象とする 3D geometry、plasma closure、particle handoff の物理契約を先に定義する。
2. rough surface を含む代表問題と、解析解または独立 solver による oracle を用意する。
3. zero/nonzero mode、height sampling、panel quadrature、outer orbit の収束条件を検証する。
4. field response、particle return、reservoir inflow の責務を別 interface として実装する。
5. checkpoint、MPI、出力、設定 schema への結合を必要最小限にし、標準の `kinetic_1d` 経路へ影響させない。

ADR 0002 は当時の設計意図と失われた能力を記録する履歴として残す。

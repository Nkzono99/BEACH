# ADR 0008: 外部境界のauthoring facade

- Status: accepted
- Date: 2026-07-23
- Supersedes: ADR 0007の公開TOMLを変更しない判断

> `linear_debye`、静的 sheath 流入補正、専用 histogram の公開・互換設定は
> [ADR 0009](0009-remove-linear-debye-and-legacy-sheath.md)で削除した。3責務の facade と
> field / particle mode から runtime contract を導出する判断は引き続き有効である。

## Context

外部境界の実行時処理はADR 0007で、流入写像、通常open面、interface輸送、delayed queueの
read-only contractへ整理された。一方、公開TOMLには実装過程で追加された
`reservoir_potential_model`、`sheath_injection_model`、`open_boundary_model`、outer `return_model`、
coupling `particle_transfer_mode`、queue flagが残った。

これらを独立な選択肢として見せると、利用者はfield modelが決めるreturn/transferの内部IDまで再指定し、
許されない組合せを説明から復元しなければならない。完全な一括presetへまとめると、field-only、
scalar inflowと通常open面の独立性、unified 3D orbitでのlocal inflowを失う。

## Decision

通常利用者向けに、z-high固定の`[external_boundary]` authoring facadeを追加する。

```text
external_boundary
├── field
├── particles
└── ordinary_open
```

- `field`は外部場と、そのmodel固有の物理・診断パラメータを所有する。
- `particles`はz-high粒子のlifecycle、reservoir流入、時間・orbit・steady-start guardを所有する。
- `ordinary_open`はouterが所有しないopen面の規則だけを所有する。

`particles.mode`は`local_source | same_batch | zhao_queue`とする。`same_batch`の実装はfieldから導出する。

| field | `same_batch` |
| --- | --- |
| `linear_debye` | linear 1D return/transfer |
| `kinetic_1d` | kinetic-profile returnと1D transfer |
| `unified_linear_response` | explicit 3D return/transfer |

`zhao_queue`は`kinetic_1d + zhao_charge_driven`だけに限定する。汎用delayed transportとして抽象化しない。

`particles.inflow_model="auto"`は、tracked linear/kinetic 1Dでは同じprofileへ流入ownershipを渡し、
それ以外では設定済みsource VDFを使う。linear/kinetic 1D ownerにscalarまたはlegacy補正を重ねる入力は拒否する。
unified 3D orbitは流入を所有しないため、local inflowを独立に選べる。

次はfacadeから既存runtime設定へ導出し、公開facadeでは受け付けない。

- `outer_plasma.return_model`
- `coupling.particle_transfer_mode`
- `coupling.outer_queue_enabled`
- `coupling.update_mode="explicit"`
- `outer_plasma.interface_z=sim.box_max[2]`
- reservoir/sheath selectorと通常open selectorの内部配置

facadeは既存の`app_config`とADR 0007の`external_boundary_contract_type`へloweringする。simulation loop、
particle stepper、outer field、return計算、queue実装は変更しない。

旧`[outer_plasma]`、`[coupling]`、関連する`[sim]` selectorはcompatibility inputとして残す。
facadeと旧selectorを同じdocumentへ書いた場合は、値が一致していても拒否する。両構文を
同じresolved runtime値へ正規化し、authoring形式はmodel fingerprintへ含めない。

`summary.txt`にはresolved inflow、ordinary open、interface transport、particle modeを記録し、
利用者が派生結果を監査できるようにする。

## Rejected alternatives

- 全設定を一つの`mode` presetへまとめる。`infinity_barrier`と`potential_barrier`、field-onlyとtracked、
  unified orbitとlocal inflowは実際に独立なので、mode数の増加または機能欠落を招く。
- 旧raw selectorをadvanced overrideとしてfacadeと混在可能にする。優先順位と二重ownershipを再導入する。
- 任意面interface、抽象基底class、strategy registryを同時に導入する。現行ownership面はz-highだけであり、
  authoring問題の解決に不要である。
- field modelからspecies、periodic backend、lower zero-mode boundary、磁場を自動変更する。これらは物理・数値条件
  なので、矛盾をerrorとして報告する。

## Consequences

利用者は物理的なfield、particle lifecycle、通常open面を選び、内部return/transfer対を入力しない。
モデル固有パラメータは責務のsubtable内に限定され、schemaの補完候補も狭くなる。

runtime contractとfingerprintは維持されるため、旧設定から作ったcheckpointを同じresolved facade設定で再開できる。
parser、Python normalizer、JSON Schemaの3経路は同じmappingを持つ必要がある。旧構文は移行期間中もテストし、
canonical examplesとdocsはfacadeだけを使う。

## Validation

1. facadeの各valid構成が旧構文と同じruntime値とfingerprintへ解決されること。
2. field-onlyのlinear/kinetic/unifiedと、tracked 1D/3Dを別々に保持すること。
3. scalar、legacy、profile inflowのownership競合を拒否すること。
4. facadeと旧selectorの混在を、旧値がdefault値でも拒否すること。
5. Zhao queueのmodel、closure、branch、steady-start、histogram、時間guardを維持すること。
6. Fortran parser、Python normalizer、schema、全canonical examplesが同じ構文を受理すること。

title: 物理再設計の完遂監査

Lang: [日本語](PhysicsRedesignCompletionAudit.md) | [English](PhysicsRedesignCompletionAudit.en.md)

# 物理再設計の完遂監査

このページは、periodic2、P0 panel、outer plasma、粒子 handoff の再設計について、実装と検証根拠を
phase ごとに追跡する記録です。設計の source of truth は
`_handoff/BEACH_periodic2_sheath_panel_redesign_implementation_plan.md` です。`_handoff` は開発時の
handoff 領域で公開サイトには同期しないため、運用上必要な結論をここに固定します。

## 判定の読み方

- **完了**: phase の production path、fail-closed 条件、対応 test tier が存在する
- **scope 外**: 設計で初めから対象外とした物理で、未実装 phase ではない
- 定量利用の可否は case ごとの収束確認を含む。[計算結果の妥当性確認](ValidationGuide.html)に従う
- test 同士の一致だけでなく、解析解、保存則、収束、MPI/restart 決定性を根拠にする

## Phase 0–9 evidence matrix

| Phase | 判定 | production owner | 主な contract test / gate |
| --- | --- | --- | --- |
| 0 契約・ledger・移行基盤 | 完了 | typed physics config、`bem_charge_ledger`、model fingerprint、versioned restart | `test_charge_ledger`、`test_model_fingerprint`、`test_restart`、`test_simulator`、Python result tests |
| 1 direct P0 panel と surface side | 完了 | `bem_panel_geometry`、`bem_panel_kernel`、`bem_panel_self_terms`、element vacuum side | `test_panel_kernel` の off-surface oracle、potential continuity、`sigma/eps0` jump、self term。`test_panel_moments` の剛体変換・頂点順序・scaling |
| 2 `k0` 分離と linear outer | 完了 | `bem_periodic_zero_mode_*`、panel spectral `k!=0` reference、`electrostatic_snapshot`、linear Debye state | `test_periodic_zero_mode`、`test_outer_plasma_linear`、`test_electrostatic_snapshot`、far-correction oracle |
| 3 ambient/interface transaction | 完了 | earliest mesh/box event、typed crossing/outcome、reservoir flux mapping、outer coupler | `test_outer_plasma_interface`、`test_interface_particle_buffer`、`test_particle_stepper`、`test_simulator`、`test_mpi_hybrid` の global reservoir count |
| 4 photoelectron transfer | 完了 | individual-return histogram、previous-batch ownership、signed charge ledger | `test_outer_plasma_photoelectron`、`test_simulator` の emission/return transaction、MPI global histogram |
| 5 full panel FMM | 完了 | panel-aware topology、near subtract/add、exact panel P2M、public C kernel contract | `test_panel_geometry_near`、`test_panel_near_correction`、`test_coulomb_fmm_core_panel`、`test_dynamics_panel_fmm`、L2 C/kernel contract |
| 6 production infinite-periodic operator | 完了 | versioned `K_periodic,k!=0 - K_shell` cache、checksum/fingerprint、root build と broadcast | `test_periodic2_operator_cache`、`test_periodic2_infinite_operator`、`test_periodic2_cached_snapshot`、2-rank `test_periodic2_operator_cache_mpi` |
| 7 nonlinear kinetic outer sheath | 完了 | stretched-grid Poisson/Robin kinetic solver、branch/applicability status、root-only collective solve、restartable profile | `test_outer_plasma_kinetic_core`、`test_outer_plasma_kinetic`、`test_outer_plasma_kinetic_runtime`、`test_restart`、2-rank kinetic snapshot broadcast |
| 8 unified outer domain | 完了 | local mean/accessibility、screened nonzero tail、unified zero mode、explicit 3D electrostatic outer orbit | `test_outer_plasma_local_mean`、`test_periodic2_nonzero_tail`、`test_electrostatic_unified`、`test_outer_plasma_orbit`、simulator interface-height surface-charge invariance |
| 9 production promotion | 完了 | portable L2、HPC L3/far/MPI gates、RSS budget、convergence artifact、日英 docs/schema/example 同期 | `make test-physics-release` の manifest と `convergence.csv`、Starlight build、docs sync check |

元の実装計画は Phase 0–8 を定義し、section 9 に production release 条件を定義しています。この監査では
その最終 promotion を Phase 9 と呼びます。新しい物理 phase を追加したという意味ではありません。

## 保存・再開契約

新規 checkpoint は schema 3 で書き、kinetic/unified outer state の次の値を復元します。

- `z`、potential、electric field、charge density の全 profile
- solver status、nonlinear iteration、residual
- interface/infinity potential、interface field、Debye length、linearity metric
- integrated charge/area と electron、ion、photoelectron、total current

schema 2 の三列 profile は read-only migration として受理しますが、完全な held state とはみなしません。
初期 guess として読み、次の kinetic refresh を強制します。schema 3 の欠損 state は fail closed です。

## Interface と保存性

`outer_plasma.interface_z == sim.box_max[3]` は field domain の切断面ではなく、粒子所有権の handoff 面です。
production simulator test は同じ束縛軌道を二つの interface 高さで実行します。低い面では explicit outer
orbit へ渡して戻し、高い面では local domain で反転させます。両 case で次を要求します。

- absorbed/escaped count が一致する
- outward/returned transaction が各 ownership 規約で閉じる
- 最終総表面電荷が一致する
- charge-ledger residual が丸め誤差内である

粒子 pusher は同時刻 `(x^n,v^n)` の midpoint Boris contract です。純磁場速度保存、time reversal、
`dt` 半減時の二次収束を release convergence artifact で確認します。通常の box 内 hot path は追加 event
loop を通らず、過大 `dt` による複数 box event は fail closed とします。

## Release evidence

正式な入口は次です。

```bash
make test-l2
make test-physics-release
```

HPC gate は L3、far-correction、2-rank MPI hybrid、2-rank MPI cache concurrency を逐次実行します。
manifest の最終 `status=passed`、各 gate の `status=passed`、既定 8 GiB 未満の最大 RSS、次の六つの
convergence category をすべて要求します。

1. `boris_dt`
2. `panel_fmm_order`
3. `rough_panel_mesh`
4. `rough_outer_grid`
5. `rough_accessibility`
6. `outer_orbit_dt`

実行日時、commit、compiler、Slurm job、elapsed time、RSS は `build/physics-release/manifest.txt` に、
数値行は `build/physics-release/convergence.csv` に記録されます。生成物は release ごとに再作成し、
repository には固定値として commit しません。

## Superseded review stages

`docs/superpowers/specs/2026-07-10-review-remediation-design.md` は重要な先行レビューですが、後発の物理
再設計計画と矛盾する箇所は superseded です。特に次を旧文書どおりの完了条件にはしません。

- point-source kernel を最終基準にする案: continuous triangle P0 panel contract に置換
- charged-walls を非中性 periodic2 の一般 closure にする案: explicit zero-mode と outer provider に置換
- analytic periodic M2L を唯一の runtime backend にする案: versioned cached operator を production 採用
- Zhao closure を全 outer model の終着点にする案: legacy injection model として分離

一方、literal output path、checkpoint file の temp-write/atomic-rename、strict history、MPI global count、
共通 field/potential snapshot、root-only cache 生成など非競合の要件は現行設計に継承されています。
旧文書が提案した generation directory/current manifest 方式や Zhao prescribed profile を含む Stage 3–6
全項目の完遂を、このページは主張しません。旧 Stage 群は physical Phase 0–9 とは別の roadmap です。

## 明示的な scope 外

次は silent approximation を入れず、別の物理設計が必要な機能です。

- 誘電分極を自己無撞着に解く dielectric boundary condition
- 3D volume PIC/Poisson space charge と nonlinear lateral outer coupling
- magnetic field を含む outer orbit。現行 explicit outer orbit は electrostatic、`b0=0` の契約
- photoelectron statistical return の着地点分布・遅延・persistent queue
- panel-aware treecode と softened continuous-panel kernel の一般定義
- substrate の `E_bottom=0` を越える電気境界条件

これらを必要とする case は対応 model へ fallback せず validation error または not-applicable になります。

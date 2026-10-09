title: Glossary

Lang: [English](Glossary.en.md) | [日本語](Glossary.md)

# Glossary

This page lists the terms used in the BEACH documentation, their Japanese equivalents, and the configuration keys or
output names they correspond to. The documentation uses these terms in prose; use the identifier column when you search a
configuration or an output. Configuration keys use the grouped (seven-group) names.

## How a run progresses

| Term | Japanese | Identifier | Meaning |
|---|---|---|---|
| batch | バッチ | `run.batch_count` | One cycle that tracks particles in a frozen field and commits the charge increments to the surface |
| batch duration | バッチ幅 | `run.batch.duration_s` | Physical time represented by one batch. It sets the injected amount for time-proportional sources and the interval between surface-charge updates |
| commit | 確定反映 | — | Adding the collected charge increments to the surface charge once, at the end of a batch |
| time step | 時間刻み | `particles.tracking.dt_s` | Time by which a particle trajectory advances in one step |
| adaptive progression | 適応進行 | `run.batch.adaptive.*` | Retrying a batch with a shorter duration so that the nonzero-mode potential change stays below a limit |
| restart | 再開 | `run.restart` | Continuing a run from the state stored in a checkpoint |

## Particles

| Term | Japanese | Identifier | Meaning |
|---|---|---|---|
| species | 粒子種 | `particles.species[].species_key` | Particles sharing charge, mass, distribution, and source |
| macro-particle | マクロ粒子 | — | Computational particle that represents many physical particles |
| weight | 重み | `sampling.weight` | Number of physical particles represented by one macro-particle |
| particle source | 粒子源 | `source.mode` | How particles are created (`volume_seed`, `plane_source`, `photo_raycast`) |
| boundary inflow | 境界流入 | `particles.species[].inflow` | Injecting particles from an external plasma through a non-periodic box face |
| reservoir | 外部プラズマ（reservoir） | `particles.reservoir` | Plasma with a known distribution assumed to exist outside the box |
| inflow mapping | 流入の写像 | `particles.reservoir.inflow_model` | How the reservoir distribution is turned into the distribution on the inflow face (`source_vdf`, `infinity_barrier`) |
| absorption | 吸収 | `absorbed_*` | A particle hits a surface triangle, leaves its charge, and is removed |
| escape | 脱出 | `escaped_*` | A particle leaves through an `open` face and does not come back |
| reflection | 反射 | `reflect`, `redistributed_reflect` | Reversing the normal velocity at a box face and continuing the trajectory |
| return | 帰還 | — | A particle emitted from a surface is absorbed by a surface again |
| unresolved particle | 未解決粒子 | `discarded_unresolved_*` | A particle that is neither absorbed nor escaped within the per-particle step limit. It is discarded at the end of the batch |
| photoelectron (PE) | 光電子 | `source.mode = "photo_raycast"` | Electron emitted from an illuminated surface |

## Surfaces and charge

| Term | Japanese | Identifier | Meaning |
|---|---|---|---|
| surface model | 表面モデル | `surface_model` | How absorbed charge is kept on a surface (`insulator`, `conductor`) |
| reaction charge | 反作用電荷 | `source.deposit_opposite_charge_on_emit` | Charge of the opposite sign left on the emitting element |
| charge closure | 電荷の閉じ方 | `charging.closure` | What a species' charge is matched to in addition to the tracked result (`neutral_return`, `fixed_current`) |
| closed photoelectrons | 閉じた光電子 | `charging.closure = "neutral_return"` | Setup that reflects photoelectrons at the top face and matches the net photoelectron current to zero |
| target current / charge | 目標電流・目標電荷 | `charging.target_*_current_a`, `fixed_*_target_charge_C` | Current matched by a fixed-current closure, and that current multiplied by the batch duration |
| weight scale | 倍率 | `*_weight_scale`, `target_over_tracked` | Ratio of the target charge to the tracked charge. The tracked distribution is scaled by this ratio before it is deposited |
| charge ledger | 電荷台帳 | `charge_ledger.csv` | Per-species totals of injected, absorbed, escaped, and corrected charge and counts |

## Outer sheath

| Term | Japanese | Identifier | Meaning |
|---|---|---|---|
| sheath closure | 外部シースとの接続 | `sheath.closure = "zero_current"` | Assuming a planar stationary sheath outside the box and taking the top-face currents, barriers, and potential reference from its zero-current root |
| zero-current root | 零電流根 | `sheath.zhao.branch` | Stationary sheath solution with zero net current to the surface (Zhao Type A / B / C) |
| wall potential ($\phi_0$) | 壁電位 | `phi_H_V` in `matching_plane_history.csv` | Potential of the surface seen by the outer sheath (the BEACH top face), measured from the upstream plasma at 0 V |
| potential minimum ($\phi_m$) | 電位極小 | — | Minimum of the potential between the upstream plasma and the wall in Type A. It acts as a barrier for electrons and photoelectrons |
| upstream electron density ($n_{e,\infty}$) | 上流電子密度 | — | Density of the upstream electron Maxwellian given by the sheath root |
| outflow refresh | 外部根の更新 | `sheath.coupling.outflow_refresh_batches` | Re-solving the sheath root periodically from the photoelectrons observed to leave through the top face |
| top face (z-high) | 上端面 | `z_high` | Upper face of the box, where the outer sheath is connected |
| potential barrier | 電位障壁 | `particles.boundary.open_model = "potential_barrier"` | Test that classifies a particle leaving an open face as reflected or escaped from its potential difference to the upstream plasma |

## Periodic fields

| Term | Japanese | Identifier | Meaning |
|---|---|---|---|
| field solver | 場ソルバ | `fields.solver.method` | How the field is computed from surface charge (`direct`, `treecode`, `fmm`) |
| field boundary | 場境界 | `fields.boundary` | Space in which the field is solved (`free`, `periodic2`) |
| finite images | 有限画像 | `fields.periodic.image_layers` | Layers of periodic copies placed explicitly around the cell |
| nonzero mode ($k\ne0$) | 面内変動成分 | `fields.periodic.backend` | Field component that varies in x/y |
| zero mode ($k=0$) | 面平均成分 | `fields.periodic.lower_boundary_model` | Field component averaged over x/y. It is set by the total charge below each height and the lower boundary condition |
| far correction | 遠方補正 | `fields.periodic.backend = "cached_kneq0"` | Adding the nonzero-mode field of periodic copies beyond the finite images with a precomputed operator |

## Outputs

| Term | Japanese | Identifier | Meaning |
|---|---|---|---|
| run summary | 実行記録 | `summary.txt` | Run statistics and resolved settings written as `key=value` |
| history | 履歴 | `*_history.csv` | Time series written per batch |

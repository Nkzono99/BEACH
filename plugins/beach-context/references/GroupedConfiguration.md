title: BEACH 設定グループと移行

Lang: [日本語](GroupedConfiguration.md) | [English](GroupedConfiguration.en.md)

# 設定グループと移行

新しい `beach.toml` は、物理条件と実行・数値設定を次の7グループで記述します。
`beachx config init` と密充填メッシュの生成コマンドも、この形式を出力します。
Fortran の `beach` と Python の `beachx lint` は同じ形式を読み取ります。

| グループ | 設定するもの |
| --- | --- |
| `run` | batch数、乱数seed、batchの時間幅、適応更新、再開 |
| `domain` | boxの位置・大きさ、周期軸 |
| `mesh` | 表面形状、OBJ、テンプレート、配置グループ |
| `particles` | 粒子境界、reservoir、追跡刻み、粒子種、供給・サンプリング |
| `fields` | 電場境界、外部E/B、solver、周期場backendとキャッシュ |
| `sheath` | 外部シースの零電流closure、枝、species role、光電子源、外部根の更新 |
| `output` | 出力先、履歴、checkpoint、診断 |

旧形式の公開済み1.6入力は1.xの間読み取れます。旧形式の読み取りは2.0で削除する予定です。
ひとつのファイルで旧形式と新形式を混ぜることはできません。
内部の物理設定・既定値・SI量の意味は維持します。密度・温度をシース側へ重複指定せず、粒子種の名前で参照します。
Pythonの読み込み関数が返す正規化済みruntime辞書は既存の形式を維持します。

既存ファイルの移行は、別ファイルへ出力してください。変換前後を検証し、既存の出力先を上書きしません。
キャッシュやOBJなどのパスの基準は従来どおりです。実行のworking directoryを維持してください。

```bash
beachx config migrate old.toml grouped.toml
beachx lint grouped.toml
beach --check-config grouped.toml
```

## 境界流入だけの電子

旧設定では `source_mode="volume_seed"` と `npcls_per_step=0` が体積からの供給を止めていました。
新設定では `source` を省略し、流入面を `inflow` に指定します。`sampling` は粒子供給量の数値表現を管理します。

旧形式:

```toml
[[particles.species]]
species_key = "electron"
q_particle = -1.602176634e-19
m_particle = 9.10938356e-31
source_mode = "volume_seed"
npcls_per_step = 0
number_density_cm3 = 5.0
temperature_ev = 10.0
target_macro_particles_per_batch = 2048
[particles.species.boundary_inflow]
z_high = "reservoir"
```

新形式:

```toml
[[particles.species]]
species_key = "electron"
charge_c = -1.602176634e-19
mass_kg = 9.10938356e-31
[particles.species.distribution]
number_density_cm3 = 5.0
temperature_ev = 10.0
[particles.species.sampling]
target_macro_particles_per_batch = 2048
[particles.species.inflow]
z_high = "reservoir"
```

`source` を置くときは `mode` が必須です。`volume_seed`、`plane_source`、`photo_raycast` の物理モデルは従来と同じです。
旧 `reservoir_face` は移行時に境界流入へ置き換えません。有限の供給面や障壁モデルが同じとは限らないためです。
`volume_macro_particles_per_batch` は実際に毎batch発生させる数です。
`target_macro_particles_per_batch` は重みを決める目標値で、実現数の保証ではありません。
`-1` による先頭種の重みの共有も維持します。

## batch時間幅と粒子追跡の刻み

旧形式:

```toml
[sim]
batch_count = 50000
batch_duration = 2.0
dt = 2.0e-12
max_step = 100000
[output]
write_potential_history = true
history_stride = 1
checkpoint_stride = 50
```

新形式:

```toml
[run]
batch_count = 50000
[run.batch]
duration_s = 2.0
[particles.tracking]
dt_s = 2.0e-12
max_steps_per_particle = 100000
[output.history]
potential = true
stride_batches = 1
[output.checkpoint]
stride_batches = 50
```

`run.batch.duration_steps` は粒子追跡の `dt_s` を掛ける従来の倍率です。`duration_s` と同時には指定しません。
batchの蓄積時間と粒子軌道の刻みは異なる役割を持ちます。
再開は `[run.restart]` の存在で指定します。`from="outputs/previous"` を省略すると、
`output.dir` にあるcheckpointを従来どおり探索します。`from` はcheckpoint一式を含むディレクトリです。再開時の `batch_count` は累積の到達目標です。
`output.diagnostics.charge_rel_change_threshold` は診断値で、早期終了条件ではありません。

## 周期場と外部シース

旧形式:

```toml
[sim]
field_solver = "fmm"
field_periodic_far_correction = "cached_kneq0"
field_periodic_cache_dir = ".beach_cache/periodic2"
[field_boundary]
mode = "periodic2"
[periodic2]
nonzero_mode_backend = "cached_kneq0"
zero_mode_policy = "exclude_k0"
lower_boundary_model = "e_bottom_zero"
[surface_current_model]
model = "zhao_stationary"
zhao_branch = "auto"
electron_species = "solar_wind_electron"
ion_species = "solar_wind_ion"
photoelectron_species = "photoelectron"
solar_elevation_deg = 60.0
photoelectron_ref_density_m3 = 6.4e7
outflow_refresh_batches = 50
```

新形式:

```toml
[fields]
boundary = "periodic2"
[fields.solver]
method = "fmm"
[fields.periodic]
backend = "cached_kneq0"
lower_boundary_model = "e_bottom_zero"
cache_dir = ".beach_cache/periodic2"
[sheath]
closure = "zero_current"
[sheath.zhao]
branch = "auto"
[sheath.species]
electron = "solar_wind_electron"
ion = "solar_wind_ion"
photoelectron = "photoelectron"
[sheath.photoelectrons]
solar_elevation_deg = 60.0
ref_density_m3 = 6.4e7
[sheath.coupling]
outflow_refresh_batches = 50
```

`closure="zero_current"` は外部 1-D シースの Zhao 零電流根 $J_z=0$ です。`branch` は `auto/a/b/c` です。
`sheath.photoelectrons` は表面放出を決める太陽高度・基準密度・倍率、`sheath.coupling.outflow_refresh_batches` は
観測 PE 流出から外部根を解き直す accepted batch 間隔です。

周期場backendは `cached_kneq0`、`panel_spectral_reference`、`finite_images` の3つです。
前の2つは非零モードと零モードを分離し、`exclude_k0` を内部で導出します。
外部根の更新には split backend（`cached_kneq0` または `panel_spectral_reference`）を使います。
`finite_images` は従来の `none/auto` の有限image計算に対応し、無限周期解を表しません。
このbackendには分離モデルのlower boundary/referenceを指定できません。
`domain.periodic_axes` と場の境界closureは別の設定です。周期軸だけから `fields.boundary` を推測しません。
`auto_fmm_min_elements` はFMMへ切り替える要素数の下限です。

`cache_dir` は非零モードoperatorのディスク保存先です。初回は幾何・数値条件から生成し、
一致するcacheがある次回実行では再利用します。各batchで作り直す粒子のフラックス表ではありません。
省略時はworking directoryの `.beach_cache/periodic2` を使います。

完全な新形式の例は [通常の境界流入](../examples/beach.toml)、
[チュートリアル](../examples/tutorial_insulator.toml)、
[定常Zhao](../examples/grouped/zero_current.toml)、
[外部根の更新](../examples/grouped/zero_current_refresh.toml) を参照してください。
物理的な成立条件と既定値の詳細は[旧形式のパラメータ詳細](Parameters.html)を併読してください。
以下の対応表で新しい記述位置を確認できます。

## 全設定の記述位置

`domain`、`mesh.mode`、`mesh.groups.<name>`、`mesh.templates` とspeciesの `boundary` は同じ構造です。
各種の速度分布、供給、サンプリング、帯電closureは一つずつ保持します。

| 旧形式のキー | 新形式のキー |
| --- | --- |
| `sim.dt` | `particles.tracking.dt_s` |
| `sim.max_step` | `particles.tracking.max_steps_per_particle` |
| `sim.rng_seed` | `run.rng_seed` |
| `sim.batch_count` | `run.batch_count` |
| `sim.batch_duration` | `run.batch.duration_s` |
| `sim.batch_duration_step` | `run.batch.duration_steps` |
| `sim.tol_rel` | `output.diagnostics.charge_rel_change_threshold` |
| `sim.q_floor` | `output.diagnostics.charge_floor_c` |
| `sim.field_solver` | `fields.solver.method` |
| `sim.field_normalization` | `fields.solver.normalization` |
| `sim.field_length_scale` | `fields.solver.length_scale_m` |
| `sim.tree_theta` | `fields.solver.tree.theta` |
| `sim.tree_leaf_max` | `fields.solver.tree.leaf_max` |
| `sim.tree_min_nelem` | `fields.solver.auto_fmm_min_elements` |
| `sim.field_periodic_image_layers` | `fields.periodic.image_layers` |
| `sim.field_periodic_far_correction` | `fields.periodic.backend` |
| `sim.field_periodic_ewald_alpha` | `fields.periodic.ewald_alpha` |
| `sim.field_periodic_ewald_layers` | `fields.periodic.ewald_layers` |
| `sim.field_periodic_cache_dir` | `fields.periodic.cache_dir` |
| `sim.field_periodic_generation_tolerance` | `fields.periodic.generation_tolerance` |
| `sim.e0` | `fields.external.electric_v_m` |
| `sim.e0_abs` | `fields.external.electric_magnitude_v_m` |
| `sim.e0_phi_xy_deg` | `fields.external.electric_azimuth_deg` |
| `sim.e0_phi_z_deg` | `fields.external.electric_elevation_deg` |
| `sim.b0` | `fields.external.magnetic_t` |
| `sim.multiple_box_events_policy` | `particles.tracking.events.policy` |
| `sim.multiple_box_events_retry_backend` | `particles.tracking.events.retry_backend` |
| `sim.multiple_box_events_soft_discard_count_grace` | `particles.tracking.events.discard_count_grace` |
| `sim.multiple_box_events_soft_discard_fraction_limit` | `particles.tracking.events.discard_fraction_limit` |
| `sim.multiple_box_events_soft_discard_abs_charge_limit` | `particles.tracking.events.discard_charge_warning_c` |
| `sim.raycast_max_bounce` | `particles.raycast.max_bounces` |
| `periodic2.lower_boundary_model` | `fields.periodic.lower_boundary_model` |
| `periodic2.reference_mode_layers` | `fields.periodic.reference.mode_layers` |
| `periodic2.panel_quadrature_order` | `fields.periodic.reference.panel_quadrature_order` |
| `periodic2.max_nonzero_mode_potential_step` | `run.batch.adaptive.max_nonzero_mode_potential_step_v` |
| `surface_current_model.model` | `sheath.closure` |
| `surface_current_model.zhao_branch` | `sheath.zhao.branch` |
| `surface_current_model.electron_species` | `sheath.species.electron` |
| `surface_current_model.ion_species` | `sheath.species.ion` |
| `surface_current_model.photoelectron_species` | `sheath.species.photoelectron` |
| `surface_current_model.solar_elevation_deg` | `sheath.photoelectrons.solar_elevation_deg` |
| `surface_current_model.photoelectron_ref_density_m3` | `sheath.photoelectrons.ref_density_m3` |
| `surface_current_model.photoelectron_source_scale` | `sheath.photoelectrons.source_scale` |
| `surface_current_model.reference_area_m2` | `sheath.reference_area_m2` |
| `surface_current_model.outflow_refresh_batches` | `sheath.coupling.outflow_refresh_batches` |
| `output.write_files` | `output.enabled` |
| `output.dir` | `output.dir` |
| `output.write_mesh_potential` | `output.final.mesh_potential` |
| `output.write_potential_history` | `output.history.potential` |
| `output.history_stride` | `output.history.stride_batches` |
| `output.checkpoint_stride` | `output.checkpoint.stride_batches` |
| `output.restart_from` | `run.restart.from` |
| `field_boundary.mode` | `fields.boundary` |
| `particle_boundary.x_low` | `particles.boundary.x_low` |
| `particle_boundary.x_high` | `particles.boundary.x_high` |
| `particle_boundary.y_low` | `particles.boundary.y_low` |
| `particle_boundary.y_high` | `particles.boundary.y_high` |
| `particle_boundary.z_low` | `particles.boundary.z_low` |
| `particle_boundary.z_high` | `particles.boundary.z_high` |
| `particle_boundary.ordinary_open_model` | `particles.boundary.open_model` |
| `reservoir.inflow_model` | `particles.reservoir.inflow_model` |
| `reservoir.phi_infty` | `particles.reservoir.phi_infty_v` |
| `reservoir.face_potential_grid_n` | `particles.reservoir.face_potential_grid_n` |
| `mesh.obj_path` | `mesh.obj.path` |
| `mesh.obj_scale` | `mesh.obj.scale` |
| `mesh.obj_rotation` | `mesh.obj.rotation` |
| `mesh.obj_offset` | `mesh.obj.offset` |
| `mesh.surface_model` | `mesh.obj.surface_model` |
| `mesh.surface_side` | `mesh.obj.surface_side` |
| `particles.species[].species_key` | `particles.species[].species_key` |
| `particles.species[].enabled` | `particles.species[].enabled` |
| `particles.species[].q_particle` | `particles.species[].charge_c` |
| `particles.species[].m_particle` | `particles.species[].mass_kg` |
| `particles.species[].velocity_distribution` | `particles.species[].distribution.model` |
| `particles.species[].number_density_cm3` | `particles.species[].distribution.number_density_cm3` |
| `particles.species[].number_density_m3` | `particles.species[].distribution.number_density_m3` |
| `particles.species[].temperature_ev` | `particles.species[].distribution.temperature_ev` |
| `particles.species[].temperature_k` | `particles.species[].distribution.temperature_k` |
| `particles.species[].drift_velocity` | `particles.species[].distribution.drift_velocity_m_s` |
| `particles.species[].velocity_grid_path` | `particles.species[].distribution.grid_path` |
| `particles.species[].velocity_grid_pdf_kind` | `particles.species[].distribution.grid_pdf_kind` |
| `particles.species[].particle_flux_m2_s` | `particles.species[].distribution.particle_flux_m2_s` |
| `particles.species[].current_density_a_m2` | `particles.species[].distribution.current_density_a_m2` |
| `particles.species[].source_mode` | `particles.species[].source.mode` |
| `particles.species[].pos_low` | `particles.species[].source.pos_low` |
| `particles.species[].pos_high` | `particles.species[].source.pos_high` |
| `particles.species[].source_normal` | `particles.species[].source.normal` |
| `particles.species[].inject_face` | `particles.species[].source.inject_face` |
| `particles.species[].ray_direction` | `particles.species[].source.ray_direction` |
| `particles.species[].inject_region_mode` | `particles.species[].source.region_mode` |
| `particles.species[].uv_low` | `particles.species[].source.uv_low` |
| `particles.species[].uv_high` | `particles.species[].source.uv_high` |
| `particles.species[].emit_current_density_a_m2` | `particles.species[].source.emit_current_density_a_m2` |
| `particles.species[].deposit_opposite_charge_on_emit` | `particles.species[].source.deposit_opposite_charge_on_emit` |
| `particles.species[].normal_drift_speed` | `particles.species[].source.normal_drift_speed_m_s` |
| `particles.species[].npcls_per_step` | `particles.species[].sampling.volume_macro_particles_per_batch` |
| `particles.species[].w_particle` | `particles.species[].sampling.weight` |
| `particles.species[].target_macro_particles_per_batch` | `particles.species[].sampling.target_macro_particles_per_batch` |
| `particles.species[].rays_per_batch` | `particles.species[].sampling.rays_per_batch` |
| `particles.species[].velocity_grid_sampling` | `particles.species[].sampling.velocity_grid_sampling` |
| `particles.species[].surface_charge_closure` | `particles.species[].charging.closure` |
| `particles.species[].target_absorbed_current_a` | `particles.species[].charging.target_absorbed_current_a` |
| `particles.species[].target_emission_current_a` | `particles.species[].charging.target_emission_current_a` |
| `particles.species[].boundary` | `particles.species[].boundary` |
| `particles.species[].boundary_inflow` | `particles.species[].inflow` |
| `output.resume` | `run.restart` |
| `periodic2.nonzero_mode_backend` | `fields.periodic.backend` |
| `periodic2.zero_mode_policy` | `exclude_k0` (derived) |

title: BEACH Configuration Groups and Migration

Lang: [日本語](GroupedConfiguration.md) | [English](GroupedConfiguration.en.md)

# Configuration groups and migration

New `beach.toml` files use seven groups. `beachx config init` and the close-packed
mesh generator emit this layout. Both the Fortran runner and Python lint read it.

| Group | Responsibility |
| --- | --- |
| `run` | Batch target, random seed, batch clock, adaptive updates, restart |
| `domain` | Box geometry and periodic topology |
| `mesh` | Surfaces, OBJ input, templates, placement groups |
| `particles` | Boundary actions, reservoir, tracking, species, sources and sampling |
| `fields` | Field boundary, imposed E/B, solver, periodic backend and cache |
| `sheath` | External zero-current closure, branch, species roles, photoelectron source and outer-root refresh |
| `output` | Files, history, checkpoints and diagnostics |

Released flat 1.6 input remains readable throughout 1.x; remove that reader at 2.0.
Do not mix layouts in one file. Physical settings, defaults and units keep their
existing meaning. Sheath plasma properties reference species by name, rather than
duplicating density and temperature. Python loading still returns the existing
normalized runtime dictionary layout.

Migration validates input and output and refuses to overwrite a destination file.
Path bases, including cache and OBJ paths, remain unchanged; keep the simulation
working directory when comparing runs.

```bash
beachx config migrate old.toml grouped.toml
beachx lint grouped.toml
beach --check-config grouped.toml
```

## Electrons supplied only by boundary inflow

Old `source_mode="volume_seed"` with `npcls_per_step=0` disabled volume supply.
Omit `source` in the new layout and declare the inflow faces under `inflow`.
`sampling` controls the numerical representation of particle supply.

Old:

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

New:

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

A declared `source` requires `mode`. `volume_seed`, `plane_source`, and `photo_raycast`
retain their physical models. Migration keeps deprecated `reservoir_face` as an
explicit source; its finite aperture and barrier are not necessarily equivalent to
boundary inflow. `volume_macro_particles_per_batch` is an actual per-batch count.
`target_macro_particles_per_batch` sets particle weight and does not guarantee an
exact realized count. The existing `-1` option to share the first species weight remains.

## Batch duration and tracking time step

Old:

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

New:

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

`run.batch.duration_steps` is the existing multiplier of tracking `dt_s`; do not
combine it with `duration_s`. Charging batch duration and trajectory time step
serve different purposes. Presence of `[run.restart]` enables resume. An optional
`from="outputs/previous"` selects a checkpoint directory; without it, search `output.dir` as
before. `batch_count` remains a cumulative target during resume.
`output.diagnostics.charge_rel_change_threshold` is a diagnostic, not an early stop.

## Periodic fields and an outer sheath

Old:

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

New:

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

`closure="zero_current"` selects the Zhao zero-current root $J_z=0$ of an outer 1-D sheath; `branch` is `auto/a/b/c`.
`sheath.photoelectrons` holds the solar elevation, reference density, and scale that set the surface emission, and
`sheath.coupling.outflow_refresh_batches` is the accepted-batch interval for re-solving the outer root from the observed
PE outflow.

Periodic backends are `cached_kneq0`, `panel_spectral_reference`, and `finite_images`.
The first two split nonzero and zero modes and derive `exclude_k0` internally.
Use a split backend (`cached_kneq0` or `panel_spectral_reference`) for the outer-root refresh.
`finite_images` names the old `none/auto` finite image computation; it is not an
infinite periodic solution and cannot accept split lower-boundary/reference settings.
Periodic topology does not implicitly select `fields.boundary`.
`auto_fmm_min_elements` is the minimum element count for automatic FMM selection.

`cache_dir` stores the geometry-dependent nonzero-mode operator on disk. Generate
it on the first use and reuse a compatible cache on subsequent runs. This is not
a particle-flux table rebuilt each batch. The default is `.beach_cache/periodic2`
relative to the working directory.

Complete examples: [boundary inflow](../examples/beach.toml),
[tutorial](../examples/tutorial_insulator.toml),
[stationary Zhao](../examples/grouped/zero_current.toml),
[outer-root refresh](../examples/grouped/zero_current_refresh.toml).
The [flat parameter reference](Parameters.en.html) retains physical constraints
and default values; use the mapping below to find their new locations.

## Complete key mapping

`domain`, `mesh.mode`, `mesh.groups.<name>`, `mesh.templates`, and species
`boundary` keep their structure. Each species retains one distribution, source,
sampling configuration and charging closure.

| Flat key | Grouped key |
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
| `surface_current_model.zhao_upstream_band_tolerance` | `sheath.zhao.upstream_band_tolerance` |
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

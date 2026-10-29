!> 現行の局所 reservoir / closed PE 設定をFortran parserで検証する。
program test_app_config_parser
  use bem_kinds, only: dp, i32
  use bem_types, only: bc_open, bc_reflect, bc_redistributed_reflect
  use bem_app_config_types, only: particle_inflow_reservoir
  use bem_app_config, only: app_config, default_app_config, load_app_config, &
                            particles_per_batch_from_config
  use bem_config_helpers, only: resolve_particle_boundaries
  use bem_external_boundary_contract, only: external_inflow_none, external_inflow_scalar_barrier
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, delete_file_if_exists, run_self_probe, &
                          file_contains
  implicit none

  type(app_config) :: cfg
  integer(i32) :: effective_boundary_low(3), effective_boundary_high(3)
  integer :: i
  integer :: probe_exit
  integer(i32) :: batch_n
  real(dp) :: electron_temperature_ev
  character(len=64) :: run_mode
  character(len=512) :: probe_config_path
  character(len=*), parameter :: zhao_magnetized_path = 'test_zhao_magnetized_tmp.toml'
  character(len=*), parameter :: zhao_generic_barrier_path = 'test_zhao_generic_barrier_tmp.toml'
  character(len=*), parameter :: zhao_no_photo_branch_path = 'test_zhao_no_photo_branch_tmp.toml'
  character(len=*), parameter :: zhao_refresh_variant_path = 'test_zhao_refresh_variant_tmp.toml'
  character(len=*), parameter :: sim_variant_path = 'test_sim_variant_tmp.toml'
  character(len=*), parameter :: fixed_current_variant_path = 'test_fixed_current_variant_tmp.toml'
  character(len=*), parameter :: input_contract_path = 'test_input_contract_tmp.toml'
  character(len=*), parameter :: config_failure_path = 'test_zhao_config_failure_tmp.log'

  call get_command_argument(1, run_mode)
  if (trim(run_mode) == '--batch-count-overflow-probe') then
    call default_app_config(cfg)
    cfg%n_particle_species = 3_i32
    cfg%particle_species(1:3)%enabled = .true.
    cfg%particle_species(1:3)%npcls_per_step = huge(0_i32)
    batch_n = particles_per_batch_from_config(cfg)
    error stop 'overflow probe unexpectedly completed'
  end if
  if (trim(run_mode) == '--config-failure-probe') then
    call get_command_argument(2, probe_config_path)
    call default_app_config(cfg)
    call load_app_config(trim(probe_config_path), cfg)
    error stop 'invalid config probe unexpectedly completed'
  end if

  call test_init(50)

  call test_begin('volume_batch_count_integer_boundary')
  call default_app_config(cfg)
  cfg%n_particle_species = 3_i32
  cfg%particle_species(1:3)%enabled = .true.
  cfg%particle_species(1)%npcls_per_step = huge(0_i32) - 1_i32
  cfg%particle_species(2)%npcls_per_step = 1_i32
  cfg%particle_species(3)%npcls_per_step = 0_i32
  call assert_equal_i32(particles_per_batch_from_config(cfg), huge(0_i32), 'largest representable batch count')
  cfg%particle_species(3)%npcls_per_step = huge(0_i32)
  cfg%particle_species(3)%enabled = .false.
  call assert_equal_i32(particles_per_batch_from_config(cfg), huge(0_i32), 'disabled species do not count')
  call run_self_probe('--batch-count-overflow-probe', config_failure_path, probe_exit)
  call assert_true(probe_exit /= 0, 'overflow must be rejected before allocation')
  call assert_true(file_contains(config_failure_path, 'Total volume particles per batch exceed'), 'overflow diagnostic')
  call delete_file_if_exists(config_failure_path)
  call test_end()

  call test_begin('grouped_input_preserves_clock_fields_and_boundary_inflow')
  call write_grouped_config(input_contract_path, '')
  call default_app_config(cfg)
  call load_app_config(input_contract_path, cfg)
  call assert_close_dp(cfg%sim%dt, 1.0e-9_dp, 1.0e-20_dp, 'tracking clock')
  call assert_close_dp(cfg%sim%batch_duration, 2.0e-6_dp, 1.0e-16_dp, 'batch clock')
  call assert_equal_i32(cfg%sim%batch_count, 12_i32, 'batch count')
  call assert_close_dp(cfg%sim%e0(3), -5.0_dp, 1.0e-12_dp, 'external E')
  call assert_equal_i32(cfg%particle_species(1)%npcls_per_step, 0_i32, 'absent source has no volume particles')
  call assert_equal_i32(cfg%particle_species(1)%boundary_inflow_high(3), particle_inflow_reservoir, 'inflow')
  call test_end()

  call test_begin('grouped_input_restart')
  call write_grouped_config(input_contract_path, '[run.restart]')
  call default_app_config(cfg)
  call load_app_config(input_contract_path, cfg)
  call assert_true(cfg%resume_output, 'restart table presence enables resume')
  call test_end()

  call test_begin('mesh_parts_public_contract')
  call write_grouped_config(input_contract_path, &
                            '[[mesh.parts]]'//new_line('a')//'kind="obj"'//new_line('a')//'path="fixture.obj"'//new_line('a')// &
                            'surface_side="normal_plus"'//new_line('a')//'group="g"')
  call assert_config_rejected(input_contract_path, 'mesh.parts[1].group is only valid for template kinds')
  call write_grouped_config(input_contract_path, &
                            '[[mesh.parts]]'//new_line('a')//'kind="plane"'//new_line('a')//'path="fixture.obj"')
  call assert_config_rejected(input_contract_path, 'mesh.parts[1].path is only valid for kind="obj"')
  call write_grouped_config(input_contract_path, '[[mesh.parts]]'//new_line('a')//'kind="obj"')
  call assert_config_rejected(input_contract_path, 'mesh.parts[1].path is required')
  call write_grouped_config(input_contract_path, '[[mesh.parts]]'//new_line('a')//'enabled=false')
  call assert_config_rejected(input_contract_path, 'mesh.parts[1].kind is required')
  call write_grouped_config(input_contract_path, '[mesh]'//new_line('a')//'mode="template"')
  call assert_config_rejected(input_contract_path, 'Unknown or mixed-layout key: mesh.mode')
  call write_grouped_config(input_contract_path, '[mesh.obj]'//new_line('a')//'path="fixture.obj"')
  call assert_config_rejected(input_contract_path, 'Unknown or mixed-layout table: mesh.obj')
  call write_grouped_config(input_contract_path, '[[mesh.templates]]'//new_line('a')//'kind="plane"')
  call assert_config_rejected(input_contract_path, 'Unknown or mixed-layout key: mesh.templates')
  call test_end()

  call test_begin('flat_input_is_rejected')
  call write_flat_config(input_contract_path)
  call assert_config_rejected(input_contract_path, 'flat beach.toml layout is not supported')
  call test_end()

  call test_begin('grouped_input_rejects_missing_source_mode')
  call write_grouped_config(input_contract_path, '[particles.species.source]')
  call assert_config_rejected(input_contract_path, 'particles.species.source requires mode')
  call test_end()

  call test_begin('grouped_input_rejects_unknown_empty_table')
  call write_grouped_config(input_contract_path, '[fields.unknown]')
  call assert_config_rejected(input_contract_path, 'Unknown or mixed-layout table: fields.unknown')
  call test_end()

  call test_begin('preflight_rejects_duplicate_species_key')
  call write_grouped_config(input_contract_path, '[[particles.species]]'//new_line('a')//'species_key = "electron"')
  call assert_config_rejected(input_contract_path, 'particles.species[2].species_key values must be unique')
  call test_end()

  call test_begin('preflight_rejects_inflow_through_reflecting_face')
  call write_grouped_config(input_contract_path, '[particles.species.boundary]'//new_line('a')//'z_high = "reflect"')
  call assert_config_rejected(input_contract_path, 'requires an open particle boundary action')
  call test_end()

  call test_begin('preflight_rejects_absorbed_target_with_the_wrong_sign')
  call write_grouped_config( &
    input_contract_path, &
    '[particles.species.charging]'//new_line('a')//'closure = "fixed_current"'//new_line('a')// &
    'target_absorbed_current_a = 1.0e-9' &
    )
  call assert_config_rejected(input_contract_path, 'same sign as particles.species[1].charge_c')
  call test_end()

  call test_begin('preflight_rejects_target_current_with_explicit_closure')
  call write_grouped_config( &
    input_contract_path, '[particles.species.charging]'//new_line('a')//'target_absorbed_current_a = -1.0e-9' &
    )
  call assert_config_rejected(input_contract_path, 'require particles.species[1].charging.closure="fixed_current"')
  call test_end()

  call test_begin('preflight_rejects_particle_boundary_on_a_periodic_face')
  call write_periodic_grouped_config(input_contract_path, '[particles.boundary]'//new_line('a')//'x_low = "reflect"')
  call assert_config_rejected(input_contract_path, 'cannot override a periodic domain face')
  call delete_file_if_exists(input_contract_path)
  call test_end()

  call test_begin('config_checks_integer_range_before_conversion')
  call write_input_contract_config(input_contract_path, '1.0e6', '', 'outputs/test', '2147483648')
  call assert_config_rejected(input_contract_path, 'rng_seed must fit a 32-bit signed integer')
  call write_input_contract_config(input_contract_path, '1.0e6', '', 'outputs/test', '-2147483649')
  call assert_config_rejected(input_contract_path, 'rng_seed must fit a 32-bit signed integer')
  call write_input_contract_config(input_contract_path, '1.0e6', '', 'outputs/test', '2147483647')
  call default_app_config(cfg)
  call load_app_config(input_contract_path, cfg)
  call assert_equal_i32(cfg%sim%rng_seed, huge(0_i32), 'maximum integer input must remain exact')
  call write_input_contract_config(input_contract_path, '1.0e6', '', 'outputs/test', '-2147483648')
  call default_app_config(cfg)
  call load_app_config(input_contract_path, cfg)
  call assert_equal_i32(cfg%sim%rng_seed, -huge(0_i32) - 1_i32, 'minimum integer input must remain exact')
  call delete_file_if_exists(input_contract_path)
  call test_end()

  call test_begin('config_rejects_nonfinite_density_with_explicit_weight')
  call write_input_contract_config(input_contract_path, 'nan', '', 'outputs/test')
  call assert_config_rejected(input_contract_path, 'number_density_m3 must be finite')
  call write_input_contract_config(input_contract_path, 'inf', '', 'outputs/test')
  call assert_config_rejected(input_contract_path, 'number_density_m3 must be finite')
  call delete_file_if_exists(input_contract_path)
  call test_end()

  call test_begin('config_rejects_nonfinite_array_element')
  call write_input_contract_config( &
    input_contract_path, '1.0e6', &
    '[particles.species.distribution]'//new_line('a')//'drift_velocity_m_s = [0.0, 0.0, nan]', 'outputs/test' &
    )
  call assert_config_rejected(input_contract_path, &
                              'particles.species[1].distribution.drift_velocity_m_s must contain finite values')
  call delete_file_if_exists(input_contract_path)
  call test_end()

  call test_begin('config_rejects_truncated_paths_and_species_names')
  call write_input_contract_config(input_contract_path, '1.0e6', '', repeat('p', 257))
  call assert_config_rejected(input_contract_path, 'output.dir is too long')
  call write_input_contract_config( &
    input_contract_path, '1.0e6', 'species_key = "'//repeat('s', 65)//'"', 'outputs/test' &
    )
  call assert_config_rejected(input_contract_path, 'species_key is too long')
  call delete_file_if_exists(input_contract_path)
  call test_end()

  call test_begin('config_preserves_full_length_output_path')
  call write_input_contract_config(input_contract_path, '1.0e6', '', repeat('p', 256))
  call default_app_config(cfg)
  call load_app_config(input_contract_path, cfg)
  call assert_true(cfg%output_dir == repeat('p', 256), 'full-length path must not be shortened')
  call delete_file_if_exists(input_contract_path)
  call test_end()

  call test_begin('disabled_species_skip_semantic_preflight')
  call write_input_contract_config( &
    input_contract_path, '1.0e6', '[[particles.species]]'//new_line('a')// &
    'enabled = false'//new_line('a')//'mass_kg = -1.0'//new_line('a')// &
    '[particles.species.source]'//new_line('a')//'mode = "unused"', &
    'outputs/test' &
    )
  call default_app_config(cfg)
  call load_app_config(input_contract_path, cfg)
  call assert_equal_i32(cfg%n_particle_species, 2_i32, 'disabled species must remain present')
  call assert_true(.not. cfg%particle_species(2)%enabled, 'disabled flag must be preserved')
  call delete_file_if_exists(input_contract_path)
  call test_end()

  call test_begin('preflight_rejects_invalid_volume_seed_values')
  call write_volume_seed_config(input_contract_path, mass='0.0')
  call assert_config_rejected(input_contract_path, 'particles.species[1].mass_kg must be finite and > 0')
  call write_volume_seed_config(input_contract_path, weight='-1.0')
  call assert_config_rejected(input_contract_path, 'particles.species[1].sampling.weight must be finite and > 0')
  call write_volume_seed_config(input_contract_path, temperature='-1.0')
  call assert_config_rejected(input_contract_path, &
                              'particles.species[1].distribution must resolve to a finite temperature >= 0')
  call write_volume_seed_config(input_contract_path, pos_high='[0.6, 0.1, 0.6]')
  call assert_config_rejected(input_contract_path, &
                              'particles.species[1].source.pos_high must be >= particles.species[1].source.pos_low')
  ! reservoir 流入だけの species も同じ質量検査を通る。
  call write_input_contract_config( &
    input_contract_path, '1.0e6', 'species_key = "electron"'//new_line('a')//'mass_kg = 0.0', 'outputs/test' &
    )
  call assert_config_rejected(input_contract_path, 'particles.species[1].mass_kg must be finite and > 0')
  call delete_file_if_exists(input_contract_path)
  call test_end()

  call test_begin('preflight_accepts_cold_point_volume_seed')
  call write_volume_seed_config(input_contract_path, temperature='0.0', pos_high='[0.2, 0.2, 0.2]')
  call default_app_config(cfg)
  call load_app_config(input_contract_path, cfg)
  call assert_close_dp(cfg%particle_species(1)%temperature_k, 0.0_dp, 0.0_dp, 'zero temperature must be accepted')
  call assert_true(all(cfg%particle_species(1)%pos_high == cfg%particle_species(1)%pos_low), &
                   'point source bounds must be accepted unchanged')
  call delete_file_if_exists(input_contract_path)
  call test_end()

  call test_begin('default_config')
  call default_app_config(cfg)
  call assert_true(trim(cfg%sim%field_solver) == 'auto', 'default field solver mismatch')
  call assert_true(trim(cfg%sim%field_bc_mode) == 'free', 'default field boundary mismatch')
  call assert_true( &
    trim(cfg%sim%field_periodic_far_correction) == 'none', &
    'default periodic far correction mismatch' &
    )
  call assert_true(trim(cfg%sim%reservoir_potential_model) == 'none', 'default inflow model mismatch')
  call assert_equal_i32(cfg%boundary_contract%inflow_map, external_inflow_none, 'default inflow contract mismatch')
  call assert_true(trim(cfg%sim%multiple_box_events_retry_backend) == 'none', 'default retry backend mismatch')
  call assert_equal_i32(cfg%surface_current%outflow_refresh_batches, 0_i32, 'default outer-root refresh mismatch')
  call assert_equal_i32( &
    cfg%sim%multiple_box_events_soft_discard_count_grace, 1000_i32, &
    'default soft-discard count grace mismatch' &
    )
  call assert_close_dp( &
    cfg%sim%multiple_box_events_soft_discard_fraction_limit, 1.0e-6_dp, 1.0e-18_dp, &
    'default soft-discard fraction limit mismatch' &
    )
  call assert_equal_i32(cfg%checkpoint_stride, 0_i32, 'default checkpoint stride mismatch')
  call test_end()

  call test_begin('parses_soft_discard_count_grace')
  call write_sim_variant(sim_variant_path, 'discard_count_grace = 7')
  call default_app_config(cfg)
  call load_app_config(sim_variant_path, cfg)
  call assert_equal_i32( &
    cfg%sim%multiple_box_events_soft_discard_count_grace, 7_i32, &
    'soft-discard count grace parse mismatch' &
    )
  call delete_file_if_exists(sim_variant_path)
  call test_end()

  call test_begin('parses_soft_discard_fraction_limit')
  call write_sim_variant(sim_variant_path, 'discard_fraction_limit = 2.5e-5')
  call default_app_config(cfg)
  call load_app_config(sim_variant_path, cfg)
  call assert_close_dp( &
    cfg%sim%multiple_box_events_soft_discard_fraction_limit, 2.5e-5_dp, 1.0e-17_dp, &
    'soft-discard fraction limit parse mismatch' &
    )
  call delete_file_if_exists(sim_variant_path)
  call test_end()

  call test_begin('upper_fourier_parses_off')
  call write_first_line_variant(zhao_refresh_variant_path, 'lower_boundary_model =', &
                                'lower_boundary_model = "e_bottom_zero"'//new_line('a')// &
                                'upper_fourier.switch_height = "off"'//new_line('a')//'upper_fourier.mode_layers = 7')
  call default_app_config(cfg)
  call load_app_config(zhao_refresh_variant_path, cfg)
  call assert_true(cfg%periodic2%upper_fourier_switch == 'off', 'switch selection parsing')
  call assert_equal_i32(cfg%periodic2%upper_fourier_mode_layers, 7_i32, 'switch mode layers parsing')
  call delete_file_if_exists(zhao_refresh_variant_path)
  call test_end()

  call test_begin('upper_fourier_parses_auto')
  call write_first_line_variant(zhao_refresh_variant_path, 'lower_boundary_model =', &
                                'lower_boundary_model = "e_bottom_zero"'//new_line('a')// &
                                'upper_fourier.switch_height = "auto"'//new_line('a')//'upper_fourier.mode_layers = 7')
  call default_app_config(cfg)
  call load_app_config(zhao_refresh_variant_path, cfg)
  call assert_true(cfg%periodic2%upper_fourier_switch == 'auto', 'switch selection parsing')
  call assert_equal_i32(cfg%periodic2%upper_fourier_mode_layers, 7_i32, 'switch mode layers parsing')
  call delete_file_if_exists(zhao_refresh_variant_path)
  call test_end()

  call test_begin('resolved_inflow_contract')
  call write_grouped_config(input_contract_path, &
                            '[particles.reservoir]'//new_line('a')//'inflow_model="infinity_barrier"')
  call default_app_config(cfg)
  call load_app_config(input_contract_path, cfg)
  call assert_equal_i32( &
    cfg%boundary_contract%inflow_map, external_inflow_scalar_barrier, 'loaded inflow contract mismatch' &
    )
  call delete_file_if_exists(input_contract_path)
  call test_end()

  call test_begin('upper_fourier_parses_height')
  call write_first_line_variant(zhao_refresh_variant_path, 'lower_boundary_model =', &
                                'lower_boundary_model = "e_bottom_zero"'//new_line('a')// &
                                'upper_fourier.switch_height = 0.0005'//new_line('a')//'upper_fourier.mode_layers = 7')
  call default_app_config(cfg)
  call load_app_config(zhao_refresh_variant_path, cfg)
  call assert_true(cfg%periodic2%upper_fourier_switch == 'height', 'switch selection parsing')
  call assert_equal_i32(cfg%periodic2%upper_fourier_mode_layers, 7_i32, 'switch mode layers parsing')
  call assert_close_dp(cfg%periodic2%upper_fourier_switch_height, 0.0005_dp, 0.0_dp, 'numeric switch height')
  call delete_file_if_exists(zhao_refresh_variant_path)
  call test_end()

  call test_begin('upper_fourier_rejects_invalid_configuration')
  call write_grouped_config(input_contract_path, &
                            '[fields.periodic]'//new_line('a')//'backend="cached_kneq0"'//new_line('a')// &
                            'lower_boundary_model="e_top_from_sheath"')
  call assert_config_rejected(input_contract_path, 'e_top_from_sheath requires zhao_stationary')
  call delete_file_if_exists(input_contract_path)
  call write_first_line_variant(zhao_refresh_variant_path, 'lower_boundary_model =', &
                                'lower_boundary_model = "e_bottom_zero"'//new_line('a')// &
                                'upper_fourier.switch_height = "invalid"')
  call assert_config_rejected(zhao_refresh_variant_path, 'must be off, auto, or a height in metres')
  call write_first_line_variant(zhao_refresh_variant_path, 'lower_boundary_model =', &
                                'lower_boundary_model = "e_bottom_zero"'//new_line('a')//'upper_fourier.mode_layers = 33')
  call assert_config_rejected(zhao_refresh_variant_path, 'mode_layers must be between 0 and 32')
  call delete_file_if_exists(zhao_refresh_variant_path)
  call test_end()

  ! 無効値 0 の adaptive 上限は split 設定ではない。明示した finite_images を別 backend へ切り替えない。
  call test_begin('disabled_adaptive_limit_keeps_finite_images')
  call write_periodic_grouped_config(input_contract_path, &
                                     '[fields.periodic]'//new_line('a')//'backend = "finite_images"'//new_line('a')// &
                                     '[run.batch.adaptive]'//new_line('a')//'max_nonzero_mode_potential_step_v = 0.0')
  call default_app_config(cfg)
  call load_app_config(input_contract_path, cfg)
  call assert_true(trim(cfg%periodic2%nonzero_mode_backend) == 'legacy_finite_images', 'finite_images backend kept')
  call assert_true(trim(cfg%periodic2%lower_boundary_model) == 'legacy_implicit', 'finite_images lower boundary kept')
  call assert_close_dp(cfg%periodic2%max_nonzero_mode_potential_step, 0.0_dp, 0.0_dp, 'disabled adaptive limit')
  call write_periodic_grouped_config(input_contract_path, &
                                     '[fields.periodic]'//new_line('a')//'backend = "finite_images"'//new_line('a')// &
                                     '[run.batch.adaptive]'//new_line('a')//'max_nonzero_mode_potential_step_v = 0.0', &
                                     method='direct')
  call assert_config_rejected(input_contract_path, 'direct requires its exact free-space kernel')
  call write_periodic_grouped_config(input_contract_path, &
                                     '[fields.periodic]'//new_line('a')//'backend = "finite_images"'//new_line('a')// &
                                     '[run.batch.adaptive]'//new_line('a')//'max_nonzero_mode_potential_step_v = 0.1')
  call assert_config_rejected(input_contract_path, 'requires fields.periodic.backend="cached_kneq0"')
  call delete_file_if_exists(input_contract_path)
  call test_end()

  call test_begin('split_periodic_accepts_upper_fourier_retry')
  call write_sim_variant(sim_variant_path, 'retry_backend = "upper_panel_fourier"')
  call default_app_config(cfg)
  call load_app_config(sim_variant_path, cfg)
  call assert_true( &
    trim(cfg%sim%multiple_box_events_retry_backend) == 'upper_panel_fourier', &
    'split periodic case must accept upper Fourier retry with abort fallback' &
    )
  call delete_file_if_exists(sim_variant_path)
  call test_end()

  call test_begin('zhao_rejects_magnetized_closure')
  call write_zhao_variant(zhao_magnetized_path, 'magnetic_t = [0.0, 0.0, 1.0e-9]', .false.)
  call assert_config_rejected(zhao_magnetized_path, 'sheath.closure="zero_current" requires fields.external.magnetic_t=[0,0,0]')
  call delete_file_if_exists(zhao_magnetized_path)
  call test_end()

  call test_begin('zhao_rejects_generic_reservoir_barrier')
  call write_zhao_variant(zhao_generic_barrier_path, '', .true.)
  call assert_config_rejected(zhao_generic_barrier_path, 'cannot be combined with the generic reservoir potential model')
  call delete_file_if_exists(zhao_generic_barrier_path)
  call test_end()

  call test_begin('zhao_fixed_current_config')
  call default_app_config(cfg)
  call load_app_config('examples/zero_current.toml', cfg)
  call assert_true(trim(cfg%surface_current%model) == 'zhao_stationary', 'Zhao current model mismatch')
  call assert_true(trim(cfg%surface_current%zhao_branch) == 'auto', 'Zhao branch mismatch')
  call assert_true( &
    all([(trim(cfg%particle_species(i)%surface_charge_closure) == 'fixed_current', i=1, 3)]), &
    'Zhao species must use fixed_current' &
    )
  call resolve_particle_boundaries( &
    cfg%sim, cfg%particle_boundary_low, cfg%particle_boundary_high, cfg%particle_species(3), &
    effective_boundary_low, effective_boundary_high &
    )
  call assert_equal_i32(effective_boundary_high(3), bc_open, 'Zhao PE z-high boundary must be open')
  call test_end()

  call test_begin('zhao_accepts_zero_electron_drift')
  call write_first_line_variant(zhao_refresh_variant_path, 'drift_velocity_m_s =', 'drift_velocity_m_s = [0.0, 0.0, 0.0]')
  call default_app_config(cfg)
  call load_app_config(zhao_refresh_variant_path, cfg)
  call assert_close_dp( &
    cfg%particle_species(1)%drift_velocity(3), 0.0_dp, 0.0_dp, 'Zhao electron drift must accept zero' &
    )
  call delete_file_if_exists(zhao_refresh_variant_path)
  call test_end()

  call test_begin('reflected_electron_source_switch')
  call write_coupling_variant(zhao_refresh_variant_path, 'outflow_refresh_batches = 1'//new_line('a')// &
                              'include_reflected_electrons = true')
  call default_app_config(cfg)
  call load_app_config(zhao_refresh_variant_path, cfg)
  call assert_true(cfg%surface_current%include_reflected_electrons, 'reflected electron switch must parse')
  call default_app_config(cfg)
  call load_app_config('examples/zero_current_refresh.toml', cfg)
  call assert_true(cfg%surface_current%include_reflected_electrons .and. &
                   .not. cfg%surface_current%has_include_reflected_electrons, 'reflected electrons are on by default')
  ! 既定の true は固定根（更新なし）でも受け付け、効かないだけ。明示した true は更新を要求する。
  call default_app_config(cfg)
  call load_app_config('examples/zero_current.toml', cfg)
  call assert_true(cfg%surface_current%include_reflected_electrons .and. &
                   cfg%surface_current%outflow_refresh_batches == 0_i32, 'default switch with a fixed root is accepted')
  call write_coupling_variant(zhao_refresh_variant_path, 'include_reflected_electrons = true')
  call assert_config_rejected(zhao_refresh_variant_path, 'include_reflected_electrons requires')
  call write_coupling_variant(zhao_refresh_variant_path, 'outflow_refresh_batches = 1'//new_line('a')// &
                              'include_reflected_electrons = false', 'examples/zero_current_no_photo.toml')
  call assert_config_rejected(zhao_refresh_variant_path, 'requires a surface-emitted electron group')
  ! PE も SEE もない構成でも、反射電子を外部源に数えれば更新できる。数えなければ源がなく設定エラー。
  call write_coupling_variant(zhao_refresh_variant_path, 'outflow_refresh_batches = 1'//new_line('a')// &
                              'include_reflected_electrons = true', 'examples/zero_current_no_photo.toml')
  call default_app_config(cfg)
  call load_app_config(zhao_refresh_variant_path, cfg)
  call assert_true(cfg%surface_current%include_reflected_electrons .and. &
                   cfg%surface_current%outflow_refresh_batches == 1_i32, 'reflected-only refresh accepted')
  ! 既定の true なら、PE も SEE もない構成でも更新だけを書けば反射電子を源にして更新できる。
  call write_coupling_variant(zhao_refresh_variant_path, 'outflow_refresh_batches = 1', 'examples/zero_current_no_photo.toml')
  call default_app_config(cfg)
  call load_app_config(zhao_refresh_variant_path, cfg)
  call assert_true(cfg%surface_current%include_reflected_electrons, 'no-emitter refresh uses the default reflected source')
  call delete_file_if_exists(zhao_refresh_variant_path)
  call test_end()

  call test_begin('zhao_rejects_electron_density_unlike_solar_wind')
  call write_first_line_variant(zhao_refresh_variant_path, 'number_density_cm3 =', 'number_density_cm3 = 5.0')
  call assert_config_rejected(zhao_refresh_variant_path, 'share the solar-wind number density')
  call delete_file_if_exists(zhao_refresh_variant_path)
  call test_end()

  call test_begin('zhao_outflow_refresh_config')
  call default_app_config(cfg)
  call load_app_config('examples/zero_current_refresh.toml', cfg)
  call assert_true(trim(cfg%surface_current%model) == 'zhao_stationary', 'refresh Zhao model mismatch')
  call assert_equal_i32(cfg%surface_current%outflow_refresh_batches, 2_i32, 'outflow refresh interval mismatch')
  call assert_true(cfg%periodic2%upper_fourier_switch == 'auto', 'grouped Fourier auto parsing')
  call assert_true(cfg%periodic2%lower_boundary_model == 'e_top_from_sheath', 'sheath Neumann boundary parsing')
  call default_app_config(cfg)
  call load_app_config('examples/zero_current.toml', cfg)
  call assert_equal_i32(cfg%surface_current%outflow_refresh_batches, 0_i32, 'fixed Zhao root must not refresh')
  call test_end()

  ! 上端面の平均電位を壁電位 phi0 に固定できない場では、流入の写像と外部障壁が約 0 V の面で判定される。
  call test_begin('zhao_requires_split_zero_mode')
  call write_first_line_variant(zhao_refresh_variant_path, 'backend =', &
                                'backend = "finite_images"')
  call assert_config_rejected(zhao_refresh_variant_path, 'requires an explicit split-zero-mode')
  call write_first_line_variant(zhao_refresh_variant_path, 'lower_boundary_model =', '# lower boundary omitted')
  call assert_config_rejected(zhao_refresh_variant_path, 'requires an explicit split-zero-mode')
  call delete_file_if_exists(zhao_refresh_variant_path)
  call test_end()

  call test_begin('zhao_upstream_band_tolerance_config')
  call default_app_config(cfg)
  call load_app_config('examples/zero_current.toml', cfg)
  call assert_close_dp(cfg%surface_current%zhao_upstream_band_tolerance, 0.0_dp, 0.0_dp, &
                       'example must leave the width tolerance off')
  call write_first_line_variant(zhao_refresh_variant_path, 'upstream_field_floor_v_m =', &
                                'upstream_band_tolerance = 0.2')
  call default_app_config(cfg)
  call load_app_config(zhao_refresh_variant_path, cfg)
  call assert_close_dp(cfg%surface_current%zhao_upstream_band_tolerance, 0.2_dp, 0.0_dp, &
                       'upstream band tolerance mismatch')
  call write_first_line_variant(zhao_refresh_variant_path, 'upstream_field_floor_v_m =', &
                                'upstream_band_tolerance = 1.0')
  call assert_config_rejected(zhao_refresh_variant_path, 'sheath.zhao.upstream_band_tolerance must be >= 0 and < 1')
  call delete_file_if_exists(zhao_refresh_variant_path)
  call test_end()

  call test_begin('zhao_upstream_field_floor_config')
  call write_first_line_variant(zhao_refresh_variant_path, 'upstream_field_floor_v_m =', &
                                'upstream_field_floor_v_m = 5.0e-3')
  call default_app_config(cfg)
  call load_app_config(zhao_refresh_variant_path, cfg)
  call assert_close_dp(cfg%surface_current%zhao_upstream_field_floor_v_m, 5.0e-3_dp, 0.0_dp, 'upstream field floor mismatch')
  call write_first_line_variant(zhao_refresh_variant_path, 'upstream_field_floor_v_m =', &
                                'upstream_field_floor_v_m = -1.0e-3')
  call assert_config_rejected(zhao_refresh_variant_path, 'sheath.zhao.upstream_field_floor_v_m must be finite and >= 0')
  call delete_file_if_exists(zhao_refresh_variant_path)
  call test_end()

  call test_begin('zhao_ion_pressure_factor_config')
  call default_app_config(cfg)
  call load_app_config('examples/zero_current.toml', cfg)
  call assert_close_dp(cfg%surface_current%zhao_ion_pressure_factor, 1.0_dp, 0.0_dp, &
                       'default ion pressure factor must be isothermal')
  call write_first_line_variant(zhao_refresh_variant_path, 'upstream_field_floor_v_m =', &
                                'ion_pressure_factor = 3.0')
  call default_app_config(cfg)
  call load_app_config(zhao_refresh_variant_path, cfg)
  call assert_close_dp(cfg%surface_current%zhao_ion_pressure_factor, 3.0_dp, 0.0_dp, 'ion pressure factor mismatch')
  call write_first_line_variant(zhao_refresh_variant_path, 'upstream_field_floor_v_m =', &
                                'ion_pressure_factor = 0.0')
  call assert_config_rejected(zhao_refresh_variant_path, 'sheath.zhao.ion_pressure_factor must be finite and > 0')
  ! 太陽風イオン (約 0.86 keV) は温かくても超音速。入口エネルギーを下回る温度は上流へ接続しない。
  call write_first_line_variant(zhao_refresh_variant_path, 'temperature_ev = 0.1', 'temperature_ev = 2000.0')
  call assert_config_rejected(zhao_refresh_variant_path, 'Zhao ion entry must be supersonic')
  call delete_file_if_exists(zhao_refresh_variant_path)
  call test_end()

  call test_begin('zhao_no_photo_fixed_current_config')
  call default_app_config(cfg)
  call load_app_config('examples/zero_current_no_photo.toml', cfg)
  call assert_true(trim(cfg%surface_current%model) == 'zhao_stationary', 'no-PE Zhao current model mismatch')
  call assert_close_dp( &
    cfg%surface_current%photoelectron_source_scale, 0.0_dp, 0.0_dp, 'no-PE Zhao source scale mismatch' &
    )
  call assert_true(len_trim(cfg%surface_current%photoelectron_species) == 0, 'no-PE Zhao must omit PE species')
  call assert_equal_i32(cfg%n_particle_species, 2_i32, 'no-PE Zhao species count mismatch')
  call assert_true( &
    all([(trim(cfg%particle_species(i)%surface_charge_closure) == 'fixed_current', i=1, 2)]), &
    'no-PE Zhao ambient species must use fixed_current' &
    )
  call write_no_photo_zhao_variant(zhao_no_photo_branch_path, 'c')
  call default_app_config(cfg)
  call load_app_config(zhao_no_photo_branch_path, cfg)
  call assert_true(trim(cfg%surface_current%zhao_branch) == 'c', 'no-PE Zhao must accept explicit Type C')
  call delete_file_if_exists(zhao_no_photo_branch_path)
  call test_end()

  call test_begin('zhao_no_photo_rejects_non_c_branch')
  call write_no_photo_zhao_variant(zhao_no_photo_branch_path, 'a')
  call assert_config_rejected(zhao_no_photo_branch_path, 'sheath.zhao.branch="auto" or "c"')
  call delete_file_if_exists(zhao_no_photo_branch_path)
  call test_end()

  call test_begin('fixed_current_config')
  call default_app_config(cfg)
  call load_app_config('tests/fortran/fixed_current.toml', cfg)
  call assert_true( &
    trim(cfg%particle_species(1)%surface_charge_closure) == 'fixed_current', &
    'fixed-current closure mismatch' &
    )
  call assert_true(cfg%particle_species(1)%has_target_absorbed_current_a, 'fixed absorbed target presence mismatch')
  call assert_close_dp( &
    cfg%particle_species(1)%target_absorbed_current_a, -2.0_dp, 1.0e-15_dp, &
    'fixed absorbed target mismatch' &
    )
  call write_fixed_absorbed_variant('duration_steps = 2.0', '-2.0')
  call default_app_config(cfg)
  call load_app_config(fixed_current_variant_path, cfg)
  call assert_true(cfg%sim%batch_duration > 0.0_dp, 'fixed-current step duration was not resolved')
  call delete_file_if_exists(fixed_current_variant_path)
  call write_fixed_emission_variant('3.0e-6')
  call default_app_config(cfg)
  call load_app_config(fixed_current_variant_path, cfg)
  call assert_true( &
    cfg%particle_species(3)%has_target_emission_current_a, &
    'fixed emission target presence mismatch' &
    )
  call assert_close_dp( &
    cfg%particle_species(3)%target_emission_current_a, 3.0e-6_dp, 1.0e-18_dp, &
    'fixed emission target mismatch' &
    )
  call delete_file_if_exists(fixed_current_variant_path)
  call test_end()

  call test_begin('tutorial_config')
  call default_app_config(cfg)
  call load_app_config('examples/tutorial_electron_beam.toml', cfg)
  call assert_true(trim(cfg%sim%field_solver) == 'direct', 'tutorial field solver mismatch')
  call assert_true(trim(cfg%sim%field_bc_mode) == 'free', 'tutorial field boundary mismatch')
  call assert_equal_i32(cfg%n_particle_species, 1_i32, 'tutorial species count mismatch')
  call assert_true(cfg%particle_species(1)%npcls_per_step > 0_i32, 'tutorial must seed a volume beam')
  call assert_equal_i32(particles_per_batch_from_config(cfg), cfg%particle_species(1)%npcls_per_step, &
                        'tutorial batch must consist of its single volume species')
  call test_end()

  call test_begin('closed_photoelectron_config')
  call default_app_config(cfg)
  call load_app_config('examples/closed_photoelectron.toml', cfg)
  call assert_equal_i32(cfg%n_particle_species, 3_i32, 'closed case species count mismatch')
  call assert_true(trim(cfg%sim%reservoir_potential_model) == 'none', 'closed case must use source VDF inflow')
  call assert_true(trim(cfg%particle_species(1)%source_mode) == 'volume_seed', 'electron source mode mismatch')
  call assert_equal_i32( &
    cfg%particle_species(1)%boundary_inflow_high(3), particle_inflow_reservoir, &
    'electron z-high boundary inflow mismatch' &
    )
  call assert_true(trim(cfg%particle_species(2)%source_mode) == 'volume_seed', 'ion source mode mismatch')
  call assert_equal_i32( &
    cfg%particle_species(2)%boundary_inflow_high(3), particle_inflow_reservoir, &
    'ion z-high boundary inflow mismatch' &
    )
  call assert_true(trim(cfg%particle_species(3)%source_mode) == 'photo_raycast', 'photoelectron source mode mismatch')
  call assert_equal_i32(cfg%particle_species(3)%boundary_high(3), bc_reflect, 'photoelectron top boundary mismatch')
  call assert_true( &
    trim(cfg%particle_species(3)%surface_charge_closure) == 'neutral_return', &
    'photoelectron closure mismatch' &
    )
  ! The photoelectron distribution table belongs to the third species only.
  electron_temperature_ev = cfg%particle_species(1)%temperature_ev
  call write_without_line_variant(input_contract_path, 'examples/closed_photoelectron.toml', 'temperature_ev = 1.5', &
                                  'temperature_ev = 1.75')
  call default_app_config(cfg)
  call load_app_config(input_contract_path, cfg)
  call assert_close_dp(cfg%particle_species(3)%temperature_ev, 1.75_dp, 0.0_dp, 'photoelectron temperature mismatch')
  call assert_close_dp(cfg%particle_species(1)%temperature_ev, electron_temperature_ev, 0.0_dp, &
                       'the photoelectron table must not change the electron temperature')
  call delete_file_if_exists(input_contract_path)
  call test_end()

  call test_begin('all_particle_boundary_faces')
  call default_app_config(cfg)
  call load_app_config('tests/fortran/particle_boundary_faces.toml', cfg)
  call assert_true(all(cfg%particle_boundary_low == [bc_open, bc_reflect, bc_open]), 'global low faces mismatch')
  call assert_true( &
    all(cfg%particle_boundary_high == [bc_reflect, bc_open, bc_redistributed_reflect]), &
    'global high faces mismatch' &
    )
  call assert_equal_i32(cfg%particle_species(1)%boundary_low(1), bc_reflect, 'species x-low mismatch')
  call assert_equal_i32(cfg%particle_species(1)%boundary_high(1), bc_open, 'species x-high mismatch')
  call assert_equal_i32(cfg%particle_species(1)%boundary_high(2), bc_reflect, 'species y-high mismatch')
  call assert_equal_i32(cfg%particle_species(1)%boundary_low(3), bc_reflect, 'species z-low mismatch')
  call resolve_particle_boundaries( &
    cfg%sim, cfg%particle_boundary_low, cfg%particle_boundary_high, cfg%particle_species(1), &
    effective_boundary_low, effective_boundary_high &
    )
  call assert_true(all(effective_boundary_low == bc_reflect), 'effective low faces mismatch')
  call assert_true( &
    all(effective_boundary_high == [bc_open, bc_reflect, bc_redistributed_reflect]), &
    'effective high faces mismatch' &
    )
  call assert_equal_i32(cfg%checkpoint_stride, 2_i32, 'checkpoint stride mismatch')
  call test_end()

  call test_begin('secondary_emission_keys_and_target_resolution')
  call default_app_config(cfg)
  call load_app_config('examples/zero_current_secondary.toml', cfg)
  call assert_true(cfg%particle_species(1)%secondary_emission_enabled, 'SEE enabled')
  call assert_equal_i32(cfg%particle_species(1)%secondary_emission_target_idx, 4_i32, 'SEE target species')
  call assert_close_dp(cfg%particle_species(1)%secondary_emission_yield_max, 1.5_dp, 0.0_dp, 'Whipple yield_max')
  call assert_close_dp(cfg%particle_species(1)%secondary_emission_energy_max_ev, 500.0_dp, 0.0_dp, 'Whipple energy_max')
  call assert_equal_i32(cfg%secondary_emission_max_generations, 3_i32, 'default generation cap')
  call default_app_config(cfg)
  call load_app_config('examples/zero_current_no_photo_secondary.toml', cfg)
  call assert_true(len_trim(cfg%surface_current%secondary_electron_species) > 0, 'SEE-only sheath role')
  ! sheath.species.secondary_electron を省いても、SEE 種が 1 つなら自動でその種が SEE の役になる。
  call write_without_line_variant(input_contract_path, 'examples/zero_current_no_photo_secondary.toml', &
                                  'secondary_electron =')
  call default_app_config(cfg)
  call load_app_config(input_contract_path, cfg)
  call assert_true(trim(cfg%surface_current%secondary_electron_species) == 'secondary_electron', &
                   'a lone secondary_emission species takes the sheath SEE role')
  call delete_file_if_exists(input_contract_path)
  call assert_true(cfg%surface_current%outflow_refresh_batches > 0_i32, 'SEE-only refresh accepted')
  call test_end()

  call test_begin('secondary_source_rejects_independent_sampling_and_current_closures')
  call write_secondary_variant(input_contract_path, '[particles.species.sampling]'//new_line('a')//'weight = 1.0')
  call assert_config_rejected(input_contract_path, 'must omit sampling')
  call write_secondary_variant(input_contract_path, &
                               '[particles.species.charging]'//new_line('a')//'closure = "fixed_current"')
  call assert_config_rejected(input_contract_path, 'must omit particles.species[4].charging.closure')
  call write_secondary_variant(input_contract_path, &
                               '[particles.species.charging]'//new_line('a')//'closure = "neutral_return"')
  call assert_config_rejected(input_contract_path, 'must omit particles.species[4].charging.closure')
  call test_end()

  call test_begin('secondary_yield_requires_valid_parameters_and_enabled_target')
  call write_secondary_variant(input_contract_path, '', target='missing')
  call assert_config_rejected(input_contract_path, 'target_species must name an enabled')
  call write_secondary_variant(input_contract_path, '', yield_max=-1.0_dp)
  call assert_config_rejected(input_contract_path, 'whipple SEE requires')
  call delete_file_if_exists(input_contract_path)
  call test_end()

  call test_summary()

contains

  subroutine write_grouped_config(path, extra)
    character(len=*), intent(in) :: path, extra
    integer :: unit
    open (newunit=unit, file=path, status='replace', action='write')
    write (unit, '(a)') '[run]'
    write (unit, '(a)') 'batch_count = 12'
    write (unit, '(a)') '[run.batch]'
    write (unit, '(a)') 'duration_s = 2.0e-6'
    write (unit, '(a)') '[particles.tracking]'
    write (unit, '(a)') 'dt_s = 1.0e-9'
    write (unit, '(a)') '[domain]'
    write (unit, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (unit, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    write (unit, '(a)') '[[particles.species]]'
    write (unit, '(a)') 'species_key = "electron"'
    write (unit, '(a)') '[particles.species.distribution]'
    write (unit, '(a)') 'number_density_m3 = 1.0e6'
    write (unit, '(a)') 'temperature_ev = 1.0'
    write (unit, '(a)') '[particles.species.sampling]'
    write (unit, '(a)') 'target_macro_particles_per_batch = 100'
    write (unit, '(a)') '[particles.species.inflow]'
    write (unit, '(a)') 'z_high = "reservoir"'
    write (unit, '(a)') '[fields.external]'
    write (unit, '(a)') 'electric_v_m = [0.0, 0.0, -5.0]'
    write (unit, '(a)') '[output]'
    write (unit, '(a)') 'enabled = false'
    if (len_trim(extra) > 0) write (unit, '(a)') extra
    close (unit)
  end subroutine write_grouped_config

  !> 2 軸周期セルの最小構成。`extra` は末尾に追加する。`method` は場の solver（既定 fmm）。
  subroutine write_periodic_grouped_config(path, extra, method)
    character(len=*), intent(in) :: path, extra
    character(len=*), intent(in), optional :: method
    integer :: unit
    open (newunit=unit, file=path, status='replace', action='write')
    write (unit, '(a)') '[run]'
    write (unit, '(a)') 'batch_count = 2'
    write (unit, '(a)') '[run.batch]'
    write (unit, '(a)') 'duration_s = 1.0e-6'
    write (unit, '(a)') '[domain]'
    write (unit, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (unit, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    write (unit, '(a)') 'periodic_axes = ["x", "y"]'
    write (unit, '(a)') '[[particles.species]]'
    write (unit, '(a)') '[particles.species.source]'
    write (unit, '(a)') 'mode = "volume_seed"'
    write (unit, '(a)') '[particles.species.sampling]'
    write (unit, '(a)') 'volume_macro_particles_per_batch = 10'
    write (unit, '(a)') '[[mesh.parts]]'
    write (unit, '(a)') 'kind = "plane"'
    write (unit, '(a)') 'surface_side = "normal_plus"'
    write (unit, '(a)') 'center = [0.5, 0.5, 0.2]'
    write (unit, '(a)') '[fields]'
    write (unit, '(a)') 'boundary = "periodic2"'
    write (unit, '(a)') '[fields.solver]'
    if (present(method)) then
      write (unit, '(a)') 'method = "'//method//'"'
    else
      write (unit, '(a)') 'method = "fmm"'
    end if
    write (unit, '(a)') '[output]'
    write (unit, '(a)') 'enabled = false'
    if (len_trim(extra) > 0) write (unit, '(a)') extra
    close (unit)
  end subroutine write_periodic_grouped_config

  !> reservoir 流入を持たない volume_seed 1 種の構成。省略した値は有効な既定値を書く。
  subroutine write_volume_seed_config(path, mass, weight, temperature, pos_high)
    character(len=*), intent(in) :: path
    character(len=*), intent(in), optional :: mass, weight, temperature, pos_high
    integer :: unit
    open (newunit=unit, file=path, status='replace', action='write')
    write (unit, '(a)') '[particles.tracking]'
    write (unit, '(a)') 'dt_s = 1.0e-9'
    write (unit, '(a)') '[domain]'
    write (unit, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (unit, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    write (unit, '(a)') '[[particles.species]]'
    if (present(mass)) write (unit, '(a)') 'mass_kg = '//mass
    if (present(temperature)) then
      write (unit, '(a)') '[particles.species.distribution]'
      write (unit, '(a)') 'temperature_k = '//temperature
    end if
    write (unit, '(a)') '[particles.species.source]'
    write (unit, '(a)') 'mode = "volume_seed"'
    write (unit, '(a)') 'pos_low = [0.2, 0.2, 0.2]'
    if (present(pos_high)) then
      write (unit, '(a)') 'pos_high = '//pos_high
    else
      write (unit, '(a)') 'pos_high = [0.8, 0.8, 0.8]'
    end if
    write (unit, '(a)') '[particles.species.sampling]'
    write (unit, '(a)') 'volume_macro_particles_per_batch = 4'
    if (present(weight)) write (unit, '(a)') 'weight = '//weight
    write (unit, '(a)') '[output]'
    write (unit, '(a)') 'enabled = false'
    close (unit)
  end subroutine write_volume_seed_config

  !> 旧flat入力そのものを拒否する公開契約のfixture。
  subroutine write_flat_config(path)
    character(len=*), intent(in) :: path
    integer :: unit
    open (newunit=unit, file=path, status='replace', action='write')
    write (unit, '(a)') '[sim]'
    write (unit, '(a)') 'dt = 1.0e-9'
    write (unit, '(a)') '[[particles.species]]'
    write (unit, '(a)') 'npcls_per_step = 1'
    close (unit)
  end subroutine write_flat_config

  subroutine write_input_contract_config(path, density, species_extra, output_directory, integer_seed)
    character(len=*), intent(in) :: path, density, species_extra, output_directory
    character(len=*), intent(in), optional :: integer_seed
    integer :: unit

    open (newunit=unit, file=path, status='replace', action='write')
    write (unit, '(a)') '[run]'
    if (present(integer_seed)) write (unit, '(a)') 'rng_seed = '//integer_seed
    write (unit, '(a)') '[run.batch]'
    write (unit, '(a)') 'duration_s = 1.0e-6'
    write (unit, '(a)') '[particles.tracking]'
    write (unit, '(a)') 'dt_s = 1.0e-9'
    write (unit, '(a)') '[domain]'
    write (unit, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (unit, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    write (unit, '(a)') '[[particles.species]]'
    if (index(species_extra, 'species_key') == 1) write (unit, '(a)') species_extra
    write (unit, '(a)') '[particles.species.sampling]'
    write (unit, '(a)') 'weight = 1.0e6'
    write (unit, '(a)') '[particles.species.inflow]'
    write (unit, '(a)') 'z_high = "reservoir"'
    if (index(species_extra, '[particles.species.distribution]') /= 1) then
      write (unit, '(a)') '[particles.species.distribution]'
    end if
    if (index(species_extra, '[particles.species.distribution]') == 1) write (unit, '(a)') species_extra
    write (unit, '(a)') 'number_density_m3 = '//trim(density)
    write (unit, '(a)') 'temperature_ev = 1.0'
    if (index(species_extra, '[[') == 1) write (unit, '(a)') species_extra
    write (unit, '(a)') '[output]'
    write (unit, '(a)') 'enabled = false'
    write (unit, '(a)') 'dir = "'//output_directory//'"'
    close (unit)
  end subroutine write_input_contract_config

  subroutine write_zhao_variant(path, sim_line, replace_reservoir)
    character(len=*), intent(in) :: path, sim_line
    logical, intent(in) :: replace_reservoir
    character(len=1024) :: line
    integer :: source_unit, output_unit, ios
    logical :: inserted_sim, replaced_reservoir

    inserted_sim = len_trim(sim_line) == 0
    replaced_reservoir = .not. replace_reservoir
    open (newunit=source_unit, file='examples/zero_current.toml', &
          status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open Zhao example fixture'
    open (newunit=output_unit, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to create Zhao invalid-config fixture'
    do
      read (source_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (replace_reservoir .and. trim(line) == 'inflow_model = "source_vdf"') then
        write (output_unit, '(a)') 'inflow_model = "infinity_barrier"'
        replaced_reservoir = .true.
      else
        if (.not. inserted_sim .and. trim(line) == '[sheath]') then
          write (output_unit, '(a)') '[fields.external]'
          write (output_unit, '(a)') trim(sim_line)
          inserted_sim = .true.
        end if
        write (output_unit, '(a)') trim(line)
      end if
    end do
    close (source_unit)
    close (output_unit)
    if (.not. inserted_sim .or. .not. replaced_reservoir) then
      error stop 'failed to specialize Zhao invalid-config fixture'
    end if
  end subroutine write_zhao_variant

  subroutine write_no_photo_zhao_variant(path, zhao_branch)
    character(len=*), intent(in) :: path, zhao_branch
    character(len=1024) :: line
    integer :: source_unit, output_unit, ios
    logical :: replaced_branch

    replaced_branch = .false.
    open (newunit=source_unit, file='examples/zero_current_no_photo.toml', &
          status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open no-PE Zhao example fixture'
    open (newunit=output_unit, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to create no-PE Zhao invalid-config fixture'
    do
      read (source_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (trim(line) == 'branch = "auto"') then
        write (output_unit, '(a)') 'branch = "'//trim(zhao_branch)//'"'
        replaced_branch = .true.
      else
        write (output_unit, '(a)') trim(line)
      end if
    end do
    close (source_unit)
    close (output_unit)
    if (.not. replaced_branch) then
      error stop 'failed to specialize no-PE Zhao invalid-config fixture'
    end if
  end subroutine write_no_photo_zhao_variant

  !> Copy the Zhao example and replace the first line starting with `key`, which belongs to the electron species.
  subroutine write_first_line_variant(path, key, replacement)
    character(len=*), intent(in) :: path, key, replacement
    character(len=1024) :: line
    integer :: source_unit, target_unit, ios
    logical :: replaced, skipping_array

    replaced = .false.
    skipping_array = .false.
    open (newunit=source_unit, file='examples/zero_current.toml', status='old', action='read', &
          iostat=ios)
    if (ios /= 0) error stop 'failed to open Zhao variant source'
    open (newunit=target_unit, file=path, status='replace', action='write')
    do
      read (source_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (skipping_array) then
        if (index(line, ']') > 0) skipping_array = .false.
        cycle
      end if
      if (replacement == 'backend = "finite_images"' .and. index(line, 'lower_boundary_model =') == 1) cycle
      if (.not. replaced .and. index(line, key) == 1) then
        skipping_array = index(line, '[') > 0 .and. index(line, ']') == 0
        write (target_unit, '(A)') replacement
        replaced = .true.
      else
        write (target_unit, '(A)') trim(line)
      end if
    end do
    close (source_unit)
    close (target_unit)
    if (.not. replaced) error stop 'Zhao variant found no line to replace'
  end subroutine write_first_line_variant

  !> Copy the split-periodic Zhao example and append one tracking event setting.
  subroutine write_sim_variant(path, sim_extra)
    character(len=*), intent(in) :: path, sim_extra
    character(len=1024) :: line
    integer :: source_unit, target_unit, ios

    open (newunit=source_unit, file='examples/zero_current_refresh.toml', status='old', action='read', &
          iostat=ios)
    if (ios /= 0) error stop 'failed to open split-periodic config fixture'
    open (newunit=target_unit, file=path, status='replace', action='write')
    do
      read (source_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      write (target_unit, '(A)') trim(line)
    end do
    write (target_unit, '(A)') '[particles.tracking.events]'
    write (target_unit, '(A)') trim(sim_extra)
    close (source_unit)
    close (target_unit)
  end subroutine write_sim_variant

  subroutine write_fixed_absorbed_variant(duration_entry, target_current)
    character(len=*), intent(in) :: duration_entry, target_current
    integer :: output_unit, ios

    open (newunit=output_unit, file=fixed_current_variant_path, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to create fixed-current config fixture'
    write (output_unit, '(a)') '[run.batch]'
    write (output_unit, '(a)') trim(duration_entry)
    write (output_unit, '(a)') ''
    write (output_unit, '(a)') '[particles]'
    write (output_unit, '(a)') '[[particles.species]]'
    write (output_unit, '(a)') 'charge_c = -1.0'
    write (output_unit, '(a)') '[particles.species.source]'
    write (output_unit, '(a)') 'mode = "volume_seed"'
    write (output_unit, '(a)') '[particles.species.sampling]'
    write (output_unit, '(a)') 'volume_macro_particles_per_batch = 1'
    write (output_unit, '(a)') '[particles.species.charging]'
    write (output_unit, '(a)') 'closure = "fixed_current"'
    if (len_trim(target_current) > 0) then
      write (output_unit, '(a)') 'target_absorbed_current_a = '//trim(target_current)
    end if
    close (output_unit)
  end subroutine write_fixed_absorbed_variant

  subroutine write_fixed_emission_variant(target_current)
    character(len=*), intent(in) :: target_current
    character(len=1024) :: line
    integer :: source_unit, output_unit, ios
    logical :: replaced

    replaced = .false.
    open (newunit=source_unit, file='examples/closed_photoelectron.toml', &
          status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open closed-photoelectron config fixture'
    open (newunit=output_unit, file=fixed_current_variant_path, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to create fixed-emission config fixture'
    do
      read (source_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (.not. replaced .and. trim(line) == 'closure = "neutral_return"') then
        write (output_unit, '(a)') 'closure = "fixed_current"'
        write (output_unit, '(a)') 'target_emission_current_a = '//trim(target_current)
        replaced = .true.
      else
        write (output_unit, '(a)') trim(line)
      end if
    end do
    close (source_unit)
    close (output_unit)
    if (.not. replaced) error stop 'failed to specialize fixed-emission config fixture'
  end subroutine write_fixed_emission_variant

  !> Copy a config and drop the first line that starts with the given key, or write `replacement` in its place.
  subroutine write_without_line_variant(path, source_path, key, replacement)
    character(len=*), intent(in) :: path, source_path, key
    character(len=*), intent(in), optional :: replacement
    character(len=1024) :: line
    integer :: source_unit, target_unit, ios
    logical :: dropped

    dropped = .false.
    open (newunit=source_unit, file=source_path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open variant source'
    open (newunit=target_unit, file=path, status='replace', action='write')
    do
      read (source_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (.not. dropped .and. index(adjustl(line), key) == 1) then
        dropped = .true.
        if (present(replacement)) write (target_unit, '(A)') replacement
        cycle
      end if
      write (target_unit, '(A)') trim(line)
    end do
    close (source_unit)
    close (target_unit)
    if (.not. dropped) error stop 'variant key was not found'
  end subroutine write_without_line_variant

  !> Copy the fixed-root Zhao example and append a [sheath.coupling] table with the given lines.
  subroutine write_coupling_variant(path, coupling_lines, source_path)
    character(len=*), intent(in) :: path, coupling_lines
    character(len=*), intent(in), optional :: source_path
    character(len=1024) :: line
    character(len=:), allocatable :: source
    integer :: source_unit, target_unit, ios

    source = 'examples/zero_current.toml'
    if (present(source_path)) source = source_path
    open (newunit=source_unit, file=source, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open Zhao coupling variant source'
    open (newunit=target_unit, file=path, status='replace', action='write')
    do
      read (source_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      write (target_unit, '(A)') trim(line)
    end do
    close (source_unit)
    write (target_unit, '(A)') '[sheath.coupling]'
    write (target_unit, '(A)') coupling_lines
    close (target_unit)
  end subroutine write_coupling_variant

  subroutine write_secondary_variant(path, extra, target, yield_max)
    character(len=*), intent(in) :: path, extra
    character(len=*), intent(in), optional :: target
    real(dp), intent(in), optional :: yield_max
    integer :: input, output, ios
    character(len=512) :: line
    open (newunit=input, file='examples/zero_current_secondary.toml', status='old', action='read')
    open (newunit=output, file=path, status='replace', action='write')
    do
      read (input, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (present(target) .and. index(line, 'target_species =') > 0) then
        line = 'target_species = "'//trim(target)//'"'
      end if
      if (present(yield_max) .and. index(line, 'yield_max =') > 0) then
        write (line, '(a,es24.16)') 'yield_max = ', yield_max
      end if
      write (output, '(a)') trim(line)
    end do
    write (output, '(a)') extra
    close (input)
    close (output)
  end subroutine write_secondary_variant

  subroutine assert_config_rejected(path, expected_fragment)
    character(len=*), intent(in) :: path, expected_fragment
    integer :: exit_status

    call run_self_probe('--config-failure-probe "'//trim(path)//'"', config_failure_path, exit_status)
    call assert_true(exit_status /= 0, 'invalid config must be rejected')
    call assert_true(file_contains(config_failure_path, trim(expected_fragment)), &
                     'config failure message mismatch: '//trim(expected_fragment))
    call delete_file_if_exists(config_failure_path)
  end subroutine assert_config_rejected
end program test_app_config_parser

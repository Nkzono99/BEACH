!> 現行の局所 reservoir / closed PE 設定をFortran parserで検証する。
program test_app_config_parser
  use bem_kinds, only: dp, i32
  use bem_types, only: bc_open, bc_reflect, bc_redistributed_reflect
  use bem_app_config_types, only: particle_inflow_reservoir
  use bem_app_config, only: app_config, default_app_config, load_app_config, &
                            particles_per_batch_from_config
  use bem_config_helpers, only: resolve_particle_boundaries
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, delete_file_if_exists
  implicit none

  type(app_config) :: cfg
  integer(i32) :: effective_boundary_low(3), effective_boundary_high(3)
  integer :: i
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
  if (trim(run_mode) == '--config-failure-probe') then
    call get_command_argument(2, probe_config_path)
    call default_app_config(cfg)
    call load_app_config(trim(probe_config_path), cfg)
    error stop 'invalid config probe unexpectedly completed'
  end if

  call test_init(29)

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

  call test_begin('grouped_input_rejects_mixed_flat_input')
  call write_grouped_config(input_contract_path, '[sim]'//new_line('a')//'dt = 1.0e-9')
  call assert_config_rejected(input_contract_path, 'Unknown or mixed-layout table: sim')
  call test_end()

  call test_begin('grouped_input_rejects_missing_source_mode')
  call write_grouped_config(input_contract_path, '[particles.species.source]')
  call assert_config_rejected(input_contract_path, 'particles.species.source requires mode')
  call test_end()

  call test_begin('grouped_input_rejects_unknown_empty_table')
  call write_grouped_config(input_contract_path, '[fields.unknown]')
  call assert_config_rejected(input_contract_path, 'Unknown or mixed-layout table: fields.unknown')
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
    input_contract_path, '1.0e6', 'drift_velocity = [0.0, 0.0, nan]', 'outputs/test' &
    )
  call assert_config_rejected(input_contract_path, 'drift_velocity must contain finite values')
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
    'enabled = false'//new_line('a')//'m_particle = -1.0'//new_line('a')//'source_mode = "unused"', &
    'outputs/test' &
    )
  call default_app_config(cfg)
  call load_app_config(input_contract_path, cfg)
  call assert_equal_i32(cfg%n_particle_species, 2_i32, 'disabled species must remain present')
  call assert_true(.not. cfg%particle_species(2)%enabled, 'disabled flag must be preserved')
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
  call assert_true(trim(cfg%sim%open_boundary_model) == 'escape', 'default open model mismatch')
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
  call write_sim_variant(sim_variant_path, 'multiple_box_events_soft_discard_count_grace = 7')
  call default_app_config(cfg)
  call load_app_config(sim_variant_path, cfg)
  call assert_equal_i32( &
    cfg%sim%multiple_box_events_soft_discard_count_grace, 7_i32, &
    'soft-discard count grace parse mismatch' &
    )
  call delete_file_if_exists(sim_variant_path)
  call test_end()

  call test_begin('parses_soft_discard_fraction_limit')
  call write_sim_variant(sim_variant_path, 'multiple_box_events_soft_discard_fraction_limit = 2.5e-5')
  call default_app_config(cfg)
  call load_app_config(sim_variant_path, cfg)
  call assert_close_dp( &
    cfg%sim%multiple_box_events_soft_discard_fraction_limit, 2.5e-5_dp, 1.0e-17_dp, &
    'soft-discard fraction limit parse mismatch' &
    )
  call delete_file_if_exists(sim_variant_path)
  call test_end()

  call test_begin('split_periodic_accepts_upper_fourier_retry')
  call write_sim_variant(sim_variant_path, 'multiple_box_events_retry_backend = "upper_panel_fourier"')
  call default_app_config(cfg)
  call load_app_config(sim_variant_path, cfg)
  call assert_true( &
    trim(cfg%sim%multiple_box_events_retry_backend) == 'upper_panel_fourier', &
    'split periodic case must accept upper Fourier retry with abort fallback' &
    )
  call delete_file_if_exists(sim_variant_path)
  call test_end()

  call test_begin('zhao_rejects_magnetized_closure')
  call write_zhao_variant(zhao_magnetized_path, 'b0 = [0.0, 0.0, 1.0e-9]', .false.)
  call assert_config_rejected(zhao_magnetized_path, 'requires sim.b0=[0,0,0]')
  call delete_file_if_exists(zhao_magnetized_path)
  call test_end()

  call test_begin('zhao_rejects_generic_reservoir_barrier')
  call write_zhao_variant(zhao_generic_barrier_path, '', .true.)
  call assert_config_rejected(zhao_generic_barrier_path, 'cannot be combined with the generic reservoir potential model')
  call delete_file_if_exists(zhao_generic_barrier_path)
  call test_end()

  call test_begin('zhao_fixed_current_config')
  call default_app_config(cfg)
  call load_app_config('examples/periodic2_zhao_fixed_current.toml', cfg)
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
  call write_first_line_variant(zhao_refresh_variant_path, 'drift_velocity =', 'drift_velocity = [0.0, 0.0, 0.0]')
  call default_app_config(cfg)
  call load_app_config(zhao_refresh_variant_path, cfg)
  call assert_close_dp( &
    cfg%particle_species(1)%drift_velocity(3), 0.0_dp, 0.0_dp, 'Zhao electron drift must accept zero' &
    )
  call delete_file_if_exists(zhao_refresh_variant_path)
  call test_end()

  call test_begin('zhao_rejects_electron_density_unlike_solar_wind')
  call write_first_line_variant(zhao_refresh_variant_path, 'number_density_cm3 =', 'number_density_cm3 = 5.0')
  call assert_config_rejected(zhao_refresh_variant_path, 'share the solar-wind number density')
  call delete_file_if_exists(zhao_refresh_variant_path)
  call test_end()

  call test_begin('zhao_outflow_refresh_config')
  call default_app_config(cfg)
  call load_app_config('examples/periodic2_zhao_outflow_refresh.toml', cfg)
  call assert_true(trim(cfg%surface_current%model) == 'zhao_stationary', 'refresh Zhao model mismatch')
  call assert_equal_i32(cfg%surface_current%outflow_refresh_batches, 2_i32, 'outflow refresh interval mismatch')
  call default_app_config(cfg)
  call load_app_config('examples/periodic2_zhao_fixed_current.toml', cfg)
  call assert_equal_i32(cfg%surface_current%outflow_refresh_batches, 0_i32, 'fixed Zhao root must not refresh')
  call test_end()

  call test_begin('zhao_outflow_refresh_requires_split_zero_mode')
  call write_refresh_variant('examples/periodic2_zhao_fixed_current.toml', zhao_refresh_variant_path)
  call assert_config_rejected(zhao_refresh_variant_path, 'requires an explicit split-zero-mode')
  call delete_file_if_exists(zhao_refresh_variant_path)
  call test_end()

  call test_begin('zhao_upstream_band_tolerance_config')
  call write_refresh_variant('examples/periodic2_zhao_fixed_current.toml', zhao_refresh_variant_path, &
                             'zhao_upstream_band_tolerance = 0.2')
  call default_app_config(cfg)
  call load_app_config(zhao_refresh_variant_path, cfg)
  call assert_close_dp(cfg%surface_current%zhao_upstream_band_tolerance, 0.2_dp, 0.0_dp, &
                       'upstream band tolerance mismatch')
  call write_refresh_variant('examples/periodic2_zhao_fixed_current.toml', zhao_refresh_variant_path, &
                             'zhao_upstream_band_tolerance = 1.0')
  call assert_config_rejected(zhao_refresh_variant_path, 'zhao_upstream_band_tolerance must be >= 0 and < 1')
  call delete_file_if_exists(zhao_refresh_variant_path)
  call test_end()

  call test_begin('zhao_no_photo_fixed_current_config')
  call default_app_config(cfg)
  call load_app_config('examples/periodic2_zhao_no_photo_fixed_current.toml', cfg)
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
  call assert_config_rejected(zhao_no_photo_branch_path, 'zhao_branch="auto" or "c"')
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
  call write_fixed_absorbed_variant('batch_duration_step = 2.0', '-2.0')
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
  call load_app_config('examples/tutorial_insulator.toml', cfg)
  call assert_true(trim(cfg%sim%field_solver) == 'direct', 'tutorial field solver mismatch')
  call assert_true(trim(cfg%sim%field_bc_mode) == 'free', 'tutorial field boundary mismatch')
  call assert_equal_i32(cfg%n_particle_species, 1_i32, 'tutorial species count mismatch')
  call assert_equal_i32(particles_per_batch_from_config(cfg), 200_i32, 'tutorial batch particle count mismatch')
  call test_end()

  call test_begin('closed_photoelectron_config')
  call default_app_config(cfg)
  call load_app_config('examples/periodic2_closed_photoelectron.toml', cfg)
  call assert_equal_i32(cfg%n_particle_species, 3_i32, 'closed case species count mismatch')
  call assert_true(trim(cfg%sim%reservoir_potential_model) == 'none', 'closed case must use source VDF inflow')
  call assert_true(trim(cfg%sim%open_boundary_model) == 'escape', 'closed case must use ordinary escape')
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
  call assert_close_dp(cfg%particle_species(3)%temperature_ev, 1.5_dp, 1.0e-15_dp, 'photoelectron temperature mismatch')
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

  subroutine write_input_contract_config(path, density, species_extra, output_directory, integer_seed)
    character(len=*), intent(in) :: path, density, species_extra, output_directory
    character(len=*), intent(in), optional :: integer_seed
    integer :: unit

    open (newunit=unit, file=path, status='replace', action='write')
    write (unit, '(a)') '[sim]'
    write (unit, '(a)') 'dt = 1.0e-9'
    write (unit, '(a)') 'batch_duration = 1.0e-6'
    if (present(integer_seed)) write (unit, '(a)') 'rng_seed = '//integer_seed
    write (unit, '(a)') '[domain]'
    write (unit, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (unit, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    write (unit, '(a)') '[[particles.species]]'
    write (unit, '(a)') 'source_mode = "reservoir_face"'
    write (unit, '(a)') 'number_density_m3 = '//trim(density)
    write (unit, '(a)') 'temperature_ev = 1.0'
    write (unit, '(a)') 'w_particle = 1.0e6'
    write (unit, '(a)') 'inject_face = "z_high"'
    write (unit, '(a)') 'pos_low = [0.0, 0.0, 1.0]'
    write (unit, '(a)') 'pos_high = [1.0, 1.0, 1.0]'
    if (len_trim(species_extra) > 0) write (unit, '(a)') species_extra
    write (unit, '(a)') '[output]'
    write (unit, '(a)') 'write_files = false'
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
    open (newunit=source_unit, file='examples/periodic2_zhao_fixed_current.toml', &
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
        write (output_unit, '(a)') trim(line)
      end if
      if (.not. inserted_sim .and. trim(line) == '[sim]') then
        write (output_unit, '(a)') trim(sim_line)
        inserted_sim = .true.
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
    open (newunit=source_unit, file='examples/periodic2_zhao_no_photo_fixed_current.toml', &
          status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open no-PE Zhao example fixture'
    open (newunit=output_unit, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to create no-PE Zhao invalid-config fixture'
    do
      read (source_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (trim(line) == 'zhao_branch = "auto"') then
        write (output_unit, '(a)') 'zhao_branch = "'//trim(zhao_branch)//'"'
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
    logical :: replaced

    replaced = .false.
    open (newunit=source_unit, file='examples/periodic2_zhao_fixed_current.toml', status='old', action='read', &
          iostat=ios)
    if (ios /= 0) error stop 'failed to open Zhao variant source'
    open (newunit=target_unit, file=path, status='replace', action='write')
    do
      read (source_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (.not. replaced .and. index(line, key) == 1) then
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

  !> Copy the split-periodic Zhao example and append one [sim] setting.
  subroutine write_sim_variant(path, sim_extra)
    character(len=*), intent(in) :: path, sim_extra
    character(len=1024) :: line
    integer :: source_unit, target_unit, ios

    open (newunit=source_unit, file='examples/periodic2_zhao_outflow_refresh.toml', status='old', action='read', &
          iostat=ios)
    if (ios /= 0) error stop 'failed to open split-periodic config fixture'
    open (newunit=target_unit, file=path, status='replace', action='write')
    do
      read (source_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      write (target_unit, '(A)') trim(line)
      if (trim(line) == '[sim]') write (target_unit, '(A)') trim(sim_extra)
    end do
    close (source_unit)
    close (target_unit)
  end subroutine write_sim_variant

  subroutine write_fixed_absorbed_variant(duration_entry, target_current)
    character(len=*), intent(in) :: duration_entry, target_current
    integer :: output_unit, ios

    open (newunit=output_unit, file=fixed_current_variant_path, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to create fixed-current config fixture'
    write (output_unit, '(a)') '[sim]'
    write (output_unit, '(a)') trim(duration_entry)
    write (output_unit, '(a)') ''
    write (output_unit, '(a)') '[particles]'
    write (output_unit, '(a)') '[[particles.species]]'
    write (output_unit, '(a)') 'npcls_per_step = 1'
    write (output_unit, '(a)') 'q_particle = -1.0'
    write (output_unit, '(a)') 'surface_charge_closure = "fixed_current"'
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
    open (newunit=source_unit, file='examples/periodic2_closed_photoelectron.toml', &
          status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open closed-photoelectron config fixture'
    open (newunit=output_unit, file=fixed_current_variant_path, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to create fixed-emission config fixture'
    do
      read (source_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (.not. replaced .and. trim(line) == 'surface_charge_closure = "neutral_return"') then
        write (output_unit, '(a)') 'surface_charge_closure = "fixed_current"'
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

  !> Copy a Zhao example and enable the outflow refresh in its surface-current table.
  subroutine write_refresh_variant(source_path, path, inserted)
    character(len=*), intent(in) :: source_path, path
    character(len=*), intent(in), optional :: inserted
    character(len=1024) :: line
    character(len=128) :: added
    integer :: source_unit, target_unit, ios

    added = 'outflow_refresh_batches = 1'
    if (present(inserted)) added = inserted

    open (newunit=source_unit, file=source_path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open Zhao refresh variant source'
    open (newunit=target_unit, file=path, status='replace', action='write')
    do
      read (source_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      write (target_unit, '(A)') trim(line)
      if (trim(line) == '[surface_current_model]') write (target_unit, '(A)') trim(added)
    end do
    close (source_unit)
    close (target_unit)
  end subroutine write_refresh_variant

  subroutine assert_config_rejected(path, expected_fragment)
    character(len=*), intent(in) :: path, expected_fragment
    character(len=1024) :: executable_path, command, child_line
    integer :: child_exit_status, child_cmd_status, child_unit, child_ios
    logical :: saw_expected

    call get_command_argument(0, executable_path)
    call delete_file_if_exists(config_failure_path)
    command = '"'//trim(executable_path)//'" --config-failure-probe "'//trim(path)//'" > '// &
              config_failure_path//' 2>&1'
    call execute_command_line(trim(command), wait=.true., exitstat=child_exit_status, cmdstat=child_cmd_status)
    call assert_equal_i32(int(child_cmd_status, i32), 0_i32, 'config failure probe command status mismatch')
    call assert_true(child_exit_status /= 0, 'invalid config must be rejected')

    saw_expected = .false.
    open (newunit=child_unit, file=config_failure_path, status='old', action='read', iostat=child_ios)
    if (child_ios /= 0) error stop 'failed to read config failure probe output'
    do
      read (child_unit, '(A)', iostat=child_ios) child_line
      if (child_ios /= 0) exit
      saw_expected = saw_expected .or. index(child_line, trim(expected_fragment)) > 0
    end do
    close (child_unit)
    call assert_true(saw_expected, 'config failure message mismatch: '//trim(expected_fragment))
    call delete_file_if_exists(config_failure_path)
  end subroutine assert_config_rejected
end program test_app_config_parser

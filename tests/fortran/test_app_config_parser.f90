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
  character(len=*), parameter :: zhao_no_photo_stale_path = 'test_zhao_no_photo_stale_tmp.toml'
  character(len=*), parameter :: zhao_no_photo_branch_path = 'test_zhao_no_photo_branch_tmp.toml'
  character(len=*), parameter :: matching_variant_path = 'test_matching_plane_variant_tmp.toml'
  character(len=*), parameter :: fixed_current_emission_path = 'test_fixed_current_emission_tmp.toml'
  character(len=*), parameter :: matching_absolute_response_path = '/tmp/beach_matching_response.csv'
  character(len=*), parameter :: config_failure_path = 'test_zhao_config_failure_tmp.log'

  call get_command_argument(1, run_mode)
  if (trim(run_mode) == '--config-failure-probe') then
    call get_command_argument(2, probe_config_path)
    call default_app_config(cfg)
    call load_app_config(trim(probe_config_path), cfg)
    error stop 'invalid config probe unexpectedly completed'
  end if

  call test_init(36)

  call test_begin('default_config')
  call default_app_config(cfg)
  call assert_true(trim(cfg%sim%field_solver) == 'auto', 'default field solver mismatch')
  call assert_true(trim(cfg%sim%field_bc_mode) == 'free', 'default field boundary mismatch')
  call assert_true(trim(cfg%sim%reservoir_potential_model) == 'none', 'default inflow model mismatch')
  call assert_true(trim(cfg%sim%open_boundary_model) == 'escape', 'default open model mismatch')
  call assert_true(trim(cfg%sim%multiple_box_events_retry_backend) == 'none', 'default retry backend mismatch')
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
  call write_matching_variant( &
    matching_variant_path, 'matching_plane_quasistatic', '', '', &
    'multiple_box_events_soft_discard_count_grace = 7', '', '' &
    )
  call default_app_config(cfg)
  call load_app_config(matching_variant_path, cfg)
  call assert_equal_i32( &
    cfg%sim%multiple_box_events_soft_discard_count_grace, 7_i32, &
    'soft-discard count grace parse mismatch' &
    )
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('parses_soft_discard_fraction_limit')
  call write_matching_variant( &
    matching_variant_path, 'matching_plane_quasistatic', '', '', &
    'multiple_box_events_soft_discard_fraction_limit = 2.5e-5', '', '' &
    )
  call default_app_config(cfg)
  call load_app_config(matching_variant_path, cfg)
  call assert_close_dp( &
    cfg%sim%multiple_box_events_soft_discard_fraction_limit, 2.5e-5_dp, 1.0e-17_dp, &
    'soft-discard fraction limit parse mismatch' &
    )
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_accepts_upper_fourier_retry')
  call write_matching_variant( &
    matching_variant_path, 'matching_plane_quasistatic', '', '', &
    'multiple_box_events_retry_backend = "upper_panel_fourier"', '', '' &
    )
  call default_app_config(cfg)
  call load_app_config(matching_variant_path, cfg)
  call assert_true( &
    trim(cfg%sim%multiple_box_events_retry_backend) == 'upper_panel_fourier', &
    'matching-plane must accept upper Fourier retry with abort fallback' &
    )
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('rejects_unknown_multiple_box_retry_backend')
  call write_matching_variant( &
    matching_variant_path, 'matching_plane_quasistatic', '', '', &
    'multiple_box_events_retry_backend = "unknown"', '', '' &
    )
  call assert_config_rejected(matching_variant_path, 'must be "none" or "upper_panel_fourier"')
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_no_photo_config')
  call default_app_config(cfg)
  call load_app_config('tests/fortran/matching_plane_no_photo.toml', cfg)
  call assert_true( &
    trim(cfg%surface_current%model) == 'matching_plane_quasistatic', 'no-PE matching-plane model mismatch' &
    )
  call assert_true(len_trim(cfg%surface_current%photoelectron_species) == 0, 'no-PE matching plane must omit PE role')
  call assert_equal_i32(cfg%n_particle_species, 2_i32, 'no-PE matching-plane species count mismatch')
  call test_end()

  call test_begin('matching_plane_no_photo_online_config')
  call write_matching_no_photo_online(matching_variant_path)
  call default_app_config(cfg)
  call load_app_config(matching_variant_path, cfg)
  call assert_true(trim(cfg%surface_current%response_backend) == 'zhao_online', 'no-PE online backend mismatch')
  call assert_true(len_trim(cfg%surface_current%photoelectron_species) == 0, 'no-PE online matching has PE role')
  call assert_equal_i32(cfg%n_particle_species, 2_i32, 'no-PE online matching species count mismatch')
  call delete_file_if_exists(matching_variant_path)
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

  call test_begin('zhao_rejects_matching_response_backend')
  call write_zhao_variant(zhao_generic_barrier_path, '', .false., 'response_backend = "zhao_online"')
  call assert_config_rejected(zhao_generic_barrier_path, 'cannot use matching-plane-specific settings')
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
  call write_no_photo_zhao_variant(zhao_no_photo_branch_path, 'c', .false.)
  call default_app_config(cfg)
  call load_app_config(zhao_no_photo_branch_path, cfg)
  call assert_true(trim(cfg%surface_current%zhao_branch) == 'c', 'no-PE Zhao must accept explicit Type C')
  call delete_file_if_exists(zhao_no_photo_branch_path)
  call test_end()

  call test_begin('zhao_no_photo_rejects_explicit_pe_key')
  call write_no_photo_zhao_variant(zhao_no_photo_stale_path, 'auto', .true.)
  call assert_config_rejected(zhao_no_photo_stale_path, 'requires omitting all photoelectron-specific Zhao settings')
  call delete_file_if_exists(zhao_no_photo_stale_path)
  call test_end()

  call test_begin('zhao_no_photo_rejects_non_c_branch')
  call write_no_photo_zhao_variant(zhao_no_photo_branch_path, 'a', .false.)
  call assert_config_rejected(zhao_no_photo_branch_path, 'zhao_branch="auto" or "c"')
  call delete_file_if_exists(zhao_no_photo_branch_path)
  call test_end()

  call test_begin('matching_plane_config_defaults_and_relative_path')
  call default_app_config(cfg)
  call load_app_config('tests/fortran/matching_plane_quasistatic.toml', cfg)
  call assert_true( &
    trim(cfg%surface_current%model) == 'matching_plane_quasistatic', 'matching-plane model mismatch' &
    )
  call assert_true(trim(cfg%surface_current%response_backend) == 'table', 'matching backend default mismatch')
  call assert_true(.not. cfg%surface_current%has_response_backend, 'implicit table backend presence mismatch')
  call assert_true(.not. cfg%surface_current%implicit_zero_mode, 'implicit zero mode must default to false')
  call assert_true( &
    trim(cfg%surface_current%response_table_path) == &
    'tests/fortran/data/matching_response_table.csv', 'matching-plane relative response path mismatch' &
    )
  call assert_close_dp(cfg%surface_current%coupling_rtol, 1.0e-4_dp, 0.0_dp, 'matching coupling rtol mismatch')
  call assert_equal_i32(cfg%surface_current%coupling_max_iterations, 20_i32, 'matching iteration limit mismatch')
  call assert_close_dp(cfg%surface_current%coupling_relaxation, 0.5_dp, 0.0_dp, 'matching relaxation mismatch')
  call assert_true(all(cfg%surface_current%coupling_atol == 0.0_dp), 'default matching coupling atol mismatch')
  call assert_true( &
    all([(trim(cfg%particle_species(i)%surface_charge_closure) == 'explicit', i=1, 3)]), &
    'matching-plane species must retain explicit surface charge closure' &
    )
  call test_end()

  call test_begin('matching_plane_accepts_implicit_zero_mode')
  call write_matching_variant( &
    matching_variant_path, 'matching_plane_quasistatic', '', &
    'implicit_zero_mode = true', '', '', '' &
    )
  call default_app_config(cfg)
  call load_app_config(matching_variant_path, cfg)
  call assert_true(cfg%surface_current%implicit_zero_mode, 'implicit zero-mode setting mismatch')
  call assert_true(cfg%surface_current%has_implicit_zero_mode, 'implicit zero-mode presence mismatch')
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_retains_absolute_response_path')
  call write_matching_variant( &
    matching_variant_path, 'matching_plane_quasistatic', matching_absolute_response_path, &
    'coupling_rtol = 0.25', '', '', '' &
    )
  call default_app_config(cfg)
  call load_app_config(matching_variant_path, cfg)
  call assert_true( &
    trim(cfg%surface_current%response_table_path) == matching_absolute_response_path, &
    'matching-plane absolute response path mismatch' &
    )
  call assert_close_dp(cfg%surface_current%coupling_rtol, 0.25_dp, 0.0_dp, 'explicit matching rtol mismatch')
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_accepts_component_absolute_tolerance')
  call write_matching_variant( &
    matching_variant_path, 'matching_plane_quasistatic', '', &
    'coupling_atol = [0.0, 0.05, 0.0, 0.0]', '', '', '' &
    )
  call default_app_config(cfg)
  call load_app_config(matching_variant_path, cfg)
  call assert_true( &
    all(cfg%surface_current%coupling_atol == [0.0_dp, 0.05_dp, 0.0_dp, 0.0_dp]), &
    'explicit matching coupling atol mismatch' &
    )
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_rejects_negative_component_absolute_tolerance')
  call write_matching_variant( &
    matching_variant_path, 'matching_plane_quasistatic', '', &
    'coupling_atol = [0.0, -0.05, 0.0, 0.0]', '', '', '' &
    )
  call assert_config_rejected(matching_variant_path, 'coupling_atol entries must be finite and >= 0')
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_accepts_zhao_online_backend')
  call write_matching_online_variant(matching_variant_path, 'b', .false., '', '')
  call default_app_config(cfg)
  call load_app_config(matching_variant_path, cfg)
  call assert_true(trim(cfg%surface_current%response_backend) == 'zhao_online', 'online backend mismatch')
  call assert_true(cfg%surface_current%has_response_backend, 'online backend presence mismatch')
  call assert_true(trim(cfg%surface_current%zhao_branch) == 'b', 'online Zhao branch mismatch')
  call assert_true(.not. cfg%surface_current%has_response_table_path, 'online backend must omit response table')
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_rejects_empty_photoelectron_role')
  call write_matching_online_variant( &
    matching_variant_path, 'auto', .false., 'photoelectron_species = "photoelectron"', &
    'photoelectron_species = ""' &
    )
  call assert_config_rejected(matching_variant_path, 'photoelectron_species must be a non-empty string')
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_online_rejects_response_table')
  call write_matching_online_variant(matching_variant_path, 'auto', .true., '', '')
  call assert_config_rejected(matching_variant_path, 'zhao_online" cannot use response_table_path')
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_online_rejects_atol_on_inactive_axis')
  call write_matching_online_variant( &
    matching_variant_path, 'auto', .false., 'photoelectron_species = "photoelectron"', &
    'photoelectron_species = "photoelectron"'//new_line('a')// &
    'coupling_atol = [0.0, 0.0, 1.0, 0.0]' &
    )
  call assert_config_rejected(matching_variant_path, 'coupling_atol must be zero on inactive ambient-outward axes')
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_online_rejects_nonunit_charge')
  call write_matching_online_variant( &
    matching_variant_path, 'auto', .false., 'q_particle = -1.602176634e-19', &
    'q_particle = -3.204353268e-19' &
    )
  call assert_config_rejected(matching_variant_path, 'requires singly charged role species')
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_rejects_zhao_settings')
  call write_matching_variant( &
    matching_variant_path, 'matching_plane_quasistatic', '', 'zhao_branch = "auto"', '', '', '' &
    )
  call assert_config_rejected(matching_variant_path, 'cannot use Zhao-specific settings')
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('none_rejects_model_specific_settings')
  call write_matching_variant(matching_variant_path, 'none', '', '', '', '', '')
  call assert_config_rejected(matching_variant_path, 'model="none" cannot use Zhao or matching-plane settings')
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_rejects_external_field')
  call write_matching_variant( &
    matching_variant_path, 'matching_plane_quasistatic', '', '', 'e0 = [0.0, 0.0, 1.0]', '', '' &
    )
  call assert_config_rejected(matching_variant_path, 'requires sim.e0=sim.b0=[0,0,0]')
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_rejects_nonopen_z_low')
  call write_matching_variant( &
    matching_variant_path, 'matching_plane_quasistatic', '', '', '', 'reflect', '' &
    )
  call assert_config_rejected(matching_variant_path, 'z-low/z-high open particle boundaries')
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_rejects_ambient_volume_reseeding')
  call write_matching_variant( &
    matching_variant_path, 'matching_plane_quasistatic', '', '', '', '', '1' &
    )
  call assert_config_rejected(matching_variant_path, 'npcls_per_step=0')
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_accepts_bounded_soft_discard')
  call write_matching_variant( &
    matching_variant_path, 'matching_plane_quasistatic', '', '', &
    'multiple_box_events_policy = "soft_discard"', '', '' &
    )
  call default_app_config(cfg)
  call load_app_config(matching_variant_path, cfg)
  call assert_true( &
    trim(cfg%sim%multiple_box_events_policy) == 'soft_discard', &
    'matching-plane must accept bounded soft discard' &
    )
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_rejects_extra_fixed_current_species')
  call write_matching_variant( &
    matching_variant_path, 'matching_plane_quasistatic', '', '', '', '', '' &
    )
  call append_matching_fixed_current_species(matching_variant_path)
  call assert_config_rejected(matching_variant_path, 'any enabled species')
  call delete_file_if_exists(matching_variant_path)
  call test_end()

  call test_begin('matching_plane_rejects_extra_enabled_species')
  call write_matching_variant( &
    matching_variant_path, 'matching_plane_quasistatic', '', '', '', '', '' &
    )
  call append_matching_explicit_species(matching_variant_path)
  call assert_config_rejected(matching_variant_path, 'exactly its enabled')
  call delete_file_if_exists(matching_variant_path)
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
  call write_fixed_emission_variant('3.0e-6')
  call default_app_config(cfg)
  call load_app_config(fixed_current_emission_path, cfg)
  call assert_true( &
    cfg%particle_species(3)%has_target_emission_current_a, &
    'fixed emission target presence mismatch' &
    )
  call assert_close_dp( &
    cfg%particle_species(3)%target_emission_current_a, 3.0e-6_dp, 1.0e-18_dp, &
    'fixed emission target mismatch' &
    )
  call delete_file_if_exists(fixed_current_emission_path)
  call write_fixed_emission_variant('-3.0e-6')
  call assert_config_rejected(fixed_current_emission_path, 'target_emission_current_a sign must oppose q_particle')
  call delete_file_if_exists(fixed_current_emission_path)
  call test_end()

  call test_begin('tutorial_config')
  call default_app_config(cfg)
  call load_app_config('examples/tutorial_insulator.toml', cfg)
  call assert_true(trim(cfg%sim%field_solver) == 'fmm', 'tutorial field solver mismatch')
  call assert_true(trim(cfg%sim%field_bc_mode) == 'periodic2', 'tutorial field boundary mismatch')
  call assert_equal_i32(cfg%n_particle_species, 1_i32, 'tutorial species count mismatch')
  call assert_equal_i32(particles_per_batch_from_config(cfg), 1_i32, 'tutorial batch particle count mismatch')
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

  subroutine write_zhao_variant(path, sim_line, replace_reservoir, surface_line)
    character(len=*), intent(in) :: path, sim_line
    logical, intent(in) :: replace_reservoir
    character(len=*), intent(in), optional :: surface_line
    character(len=1024) :: line
    integer :: source_unit, output_unit, ios
    logical :: inserted_sim, inserted_surface, replaced_reservoir

    inserted_sim = len_trim(sim_line) == 0
    inserted_surface = .not. present(surface_line)
    if (present(surface_line)) inserted_surface = len_trim(surface_line) == 0
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
      if (present(surface_line)) then
        if (.not. inserted_surface .and. trim(line) == 'model = "zhao_stationary"') then
          write (output_unit, '(a)') trim(surface_line)
          inserted_surface = .true.
        end if
      end if
    end do
    close (source_unit)
    close (output_unit)
    if (.not. inserted_sim .or. .not. inserted_surface .or. .not. replaced_reservoir) then
      error stop 'failed to specialize Zhao invalid-config fixture'
    end if
  end subroutine write_zhao_variant

  subroutine write_no_photo_zhao_variant(path, zhao_branch, include_stale_photo_setting)
    character(len=*), intent(in) :: path, zhao_branch
    logical, intent(in) :: include_stale_photo_setting
    character(len=1024) :: line
    integer :: source_unit, output_unit, ios
    logical :: replaced_branch, inserted_stale_photo_setting

    replaced_branch = .false.
    inserted_stale_photo_setting = .not. include_stale_photo_setting
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
      if (include_stale_photo_setting .and. trim(line) == 'photoelectron_source_scale = 0.0') then
        write (output_unit, '(a)') 'solar_elevation_deg = 0.0'
        inserted_stale_photo_setting = .true.
      end if
    end do
    close (source_unit)
    close (output_unit)
    if (.not. replaced_branch .or. .not. inserted_stale_photo_setting) then
      error stop 'failed to specialize no-PE Zhao invalid-config fixture'
    end if
  end subroutine write_no_photo_zhao_variant

  subroutine write_matching_variant( &
    path, model_name, response_path, surface_extra, sim_extra, z_low_action, electron_npcls &
    )
    character(len=*), intent(in) :: path, model_name, response_path, surface_extra, sim_extra, z_low_action
    character(len=*), intent(in) :: electron_npcls
    character(len=1024) :: line
    integer :: source_unit, output_unit, ios
    logical :: replaced_model, replaced_response, inserted_surface_extra, inserted_sim_extra, replaced_z_low
    logical :: replaced_electron_npcls

    replaced_model = trim(model_name) == 'matching_plane_quasistatic'
    replaced_response = len_trim(response_path) == 0
    inserted_surface_extra = len_trim(surface_extra) == 0
    inserted_sim_extra = len_trim(sim_extra) == 0
    replaced_z_low = len_trim(z_low_action) == 0
    replaced_electron_npcls = len_trim(electron_npcls) == 0
    open (newunit=source_unit, file='tests/fortran/matching_plane_quasistatic.toml', &
          status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open matching-plane config fixture'
    open (newunit=output_unit, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to create matching-plane config variant'
    do
      read (source_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (.not. replaced_model .and. trim(line) == 'model = "matching_plane_quasistatic"') then
        write (output_unit, '(a)') 'model = "'//trim(model_name)//'"'
        replaced_model = .true.
      else if (len_trim(response_path) > 0 .and. &
               trim(line) == 'response_table_path = "data/matching_response_table.csv"') then
        write (output_unit, '(a)') 'response_table_path = "'//trim(response_path)//'"'
        replaced_response = .true.
      else if (.not. replaced_z_low .and. trim(line) == 'z_low = "open"') then
        write (output_unit, '(a)') 'z_low = "'//trim(z_low_action)//'"'
        replaced_z_low = .true.
      else if (.not. replaced_electron_npcls .and. trim(line) == 'npcls_per_step = 0') then
        write (output_unit, '(a)') 'npcls_per_step = '//trim(electron_npcls)
        replaced_electron_npcls = .true.
      else
        write (output_unit, '(a)') trim(line)
      end if
      if (.not. inserted_surface_extra .and. trim(line) == 'model = "matching_plane_quasistatic"') then
        write (output_unit, '(a)') trim(surface_extra)
        inserted_surface_extra = .true.
      end if
      if (.not. inserted_sim_extra .and. trim(line) == '[sim]') then
        write (output_unit, '(a)') trim(sim_extra)
        inserted_sim_extra = .true.
      end if
    end do
    close (source_unit)
    close (output_unit)
    if (.not. replaced_model .or. .not. replaced_response .or. .not. inserted_surface_extra .or. &
        .not. inserted_sim_extra .or. .not. replaced_z_low .or. .not. replaced_electron_npcls) then
      error stop 'failed to specialize matching-plane config fixture'
    end if
  end subroutine write_matching_variant

  subroutine write_matching_online_variant(path, zhao_branch, keep_response_path, replacement_from, replacement_to)
    character(len=*), intent(in) :: path, zhao_branch, replacement_from, replacement_to
    logical, intent(in) :: keep_response_path
    character(len=1024) :: line
    integer :: source_unit, output_unit, ios
    logical :: inserted_backend, inserted_branch, handled_response, replaced_value

    inserted_backend = .false.
    inserted_branch = len_trim(zhao_branch) == 0
    handled_response = .false.
    replaced_value = len_trim(replacement_from) == 0
    open (newunit=source_unit, file='tests/fortran/matching_plane_quasistatic.toml', &
          status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open matching-plane config fixture'
    open (newunit=output_unit, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to create online matching-plane config variant'
    do
      read (source_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (trim(line) == 'model = "matching_plane_quasistatic"') then
        write (output_unit, '(a)') trim(line)
        write (output_unit, '(a)') 'response_backend = "zhao_online"'
        inserted_backend = .true.
        if (.not. inserted_branch) then
          write (output_unit, '(a)') 'zhao_branch = "'//trim(zhao_branch)//'"'
          inserted_branch = .true.
        end if
      else if (trim(line) == 'response_table_path = "data/matching_response_table.csv"') then
        handled_response = .true.
        if (keep_response_path) write (output_unit, '(a)') trim(line)
      else if (.not. replaced_value .and. trim(line) == trim(replacement_from)) then
        write (output_unit, '(a)') trim(replacement_to)
        replaced_value = .true.
      else
        write (output_unit, '(a)') trim(line)
      end if
    end do
    close (source_unit)
    close (output_unit)
    if (.not. inserted_backend .or. .not. inserted_branch .or. .not. handled_response .or. &
        .not. replaced_value) then
      error stop 'failed to specialize online matching-plane config fixture'
    end if
  end subroutine write_matching_online_variant

  subroutine write_matching_no_photo_online(path)
    character(len=*), intent(in) :: path
    character(len=1024) :: line
    integer :: source_unit, output_unit, ios

    open (newunit=source_unit, file='tests/fortran/matching_plane_no_photo.toml', &
          status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open no-PE matching-plane config fixture'
    open (newunit=output_unit, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to create no-PE online matching-plane config fixture'
    do
      read (source_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (trim(line) == 'model = "matching_plane_quasistatic"') then
        write (output_unit, '(a)') trim(line)
        write (output_unit, '(a)') 'response_backend = "zhao_online"'
        write (output_unit, '(a)') 'zhao_branch = "auto"'
      else if (trim(line) /= 'response_table_path = "data/matching_response_table.csv"') then
        write (output_unit, '(a)') trim(line)
      end if
    end do
    close (source_unit)
    close (output_unit)
  end subroutine write_matching_no_photo_online

  subroutine append_matching_fixed_current_species(path)
    character(len=*), intent(in) :: path
    integer :: output_unit, ios

    open (newunit=output_unit, file=trim(path), status='old', position='append', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to append matching-plane fixed-current species fixture'
    write (output_unit, '(a)') ''
    write (output_unit, '(a)') '[[particles.species]]'
    write (output_unit, '(a)') 'species_key = "diagnostic_electron"'
    write (output_unit, '(a)') 'source_mode = "volume_seed"'
    write (output_unit, '(a)') 'npcls_per_step = 0'
    write (output_unit, '(a)') 'q_particle = -1.602176634e-19'
    write (output_unit, '(a)') 'm_particle = 9.1093837139e-31'
    write (output_unit, '(a)') 'surface_charge_closure = "fixed_current"'
    write (output_unit, '(a)') 'target_absorbed_current_a = -1.0e-9'
    close (output_unit)
  end subroutine append_matching_fixed_current_species

  subroutine append_matching_explicit_species(path)
    character(len=*), intent(in) :: path
    integer :: output_unit, ios

    open (newunit=output_unit, file=trim(path), status='old', position='append', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to append matching-plane explicit species fixture'
    write (output_unit, '(a)') ''
    write (output_unit, '(a)') '[[particles.species]]'
    write (output_unit, '(a)') 'species_key = "diagnostic_electron"'
    write (output_unit, '(a)') 'source_mode = "volume_seed"'
    write (output_unit, '(a)') 'npcls_per_step = 0'
    write (output_unit, '(a)') 'q_particle = -1.602176634e-19'
    write (output_unit, '(a)') 'm_particle = 9.1093837139e-31'
    write (output_unit, '(a)') 'surface_charge_closure = "explicit"'
    close (output_unit)
  end subroutine append_matching_explicit_species

  subroutine write_fixed_emission_variant(target_current)
    character(len=*), intent(in) :: target_current
    character(len=1024) :: line
    integer :: source_unit, output_unit, ios
    logical :: replaced

    replaced = .false.
    open (newunit=source_unit, file='examples/periodic2_closed_photoelectron.toml', &
          status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open closed-photoelectron config fixture'
    open (newunit=output_unit, file=fixed_current_emission_path, status='replace', action='write', iostat=ios)
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

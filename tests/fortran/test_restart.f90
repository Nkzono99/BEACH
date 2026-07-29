!> 再開用チェックポイントの保存/復元を検証するテスト。
program test_restart
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: pi, qe
  use bem_mesh, only: init_mesh, prepare_periodic2_collision_mesh
  use bem_restart, only: load_restart_checkpoint, write_rng_state_file, write_macro_residuals_file, &
                         restart_rng_state_path, restart_macro_residual_path, validate_restart_contract, &
                         restart_contract_ok, restart_contract_mismatch
  use bem_output_writer, only: write_result_files
  use bem_simulator, only: run_absorption_insulator
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, &
                            seed_particles_from_config, build_mesh_from_config
  use bem_charge_ledger, only: charge_ledger_type
  use bem_electrostatic_snapshot, only: electrostatic_diagnostics_type, electrostatic_restart_state_type
  use bem_filesystem, only: atomic_rename, filesystem_success
  use bem_types, only: mesh_type, sim_stats, injection_state, bc_open, bc_periodic
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_equal_i64, assert_close_dp, assert_allclose_1d, &
                          delete_file_if_exists, ensure_directory, remove_empty_directory
  implicit none

  type(mesh_type) :: mesh
  type(sim_stats) :: stats
  type(injection_state) :: state
  type(app_config) :: cfg, cfg_changed, implicit_cfg
  type(charge_ledger_type) :: ledger, restored_ledger
  type(electrostatic_diagnostics_type) :: electrostatic_diagnostics
  type(electrostatic_restart_state_type) :: electrostatic_state
  logical :: has_restart, exists
  logical :: saw_missing_scale_error
  integer(i32) :: contract_status
  integer :: child_exit_status, child_cmd_status, probe_unit, probe_ios
  character(len=256) :: contract_message
  character(len=512) :: profile_header
  character(len=512) :: probe_line
  integer :: profile_unit, profile_ios
  character(len=1024) :: rng_rank_path, residual_global_path
  character(len=1024) :: executable_path, probe_command
  character(len=64) :: run_mode
  character(len=*), parameter :: out_dir = 'test_restart_tmp'
  character(len=*), parameter :: implicit_out_dir = 'test_restart_implicit_zhao_tmp'
  character(len=*), parameter :: dynamic_split_out_dir = 'test_restart_dynamic_zhao_split_tmp'
  character(len=*), parameter :: missing_scale_probe_log = 'test_restart_missing_zhao_scale_probe.log'

  call get_command_argument(1, run_mode)
  if (trim(run_mode) == 'missing_implicit_zhao_source_scale_probe') then
    call run_missing_implicit_zhao_source_scale_probe()
    error stop 'missing implicit-Zhao source-scale probe unexpectedly completed'
  end if

  call delete_file_if_exists(out_dir//'/summary.txt')
  call delete_file_if_exists(out_dir//'/charges.csv')
  call delete_file_if_exists(out_dir//'/rng_state.txt')
  call delete_file_if_exists(out_dir//'/macro_residuals.csv')
  call delete_file_if_exists(out_dir//'/rng_state_rank00001.txt')
  call delete_file_if_exists(out_dir//'/macro_residuals_rank00001.csv')
  call delete_file_if_exists(out_dir//'/outer_plasma_profile.csv')
  call delete_file_if_exists(out_dir//'/mesh_triangles.csv')
  call delete_file_if_exists(out_dir//'/mesh_sources.csv')
  call remove_empty_directory(out_dir)
  call cleanup_restart_fixture(implicit_out_dir)
  call cleanup_restart_fixture(dynamic_split_out_dir)
  call delete_file_if_exists(missing_scale_probe_log)

  call test_init(9)

  call test_begin('missing_checkpoint')
  call build_test_mesh(mesh)
  call load_restart_checkpoint('test_restart_missing', mesh, stats, has_restart)
  call assert_true(.not. has_restart, 'missing checkpoint should not be reported as restart')
  call assert_equal_i32(stats%batches, 0_i32, 'missing checkpoint should keep stats at defaults')
  call test_end()

  call test_begin('kinetic_outer_profile_checkpoint')
  call default_app_config(cfg)
  cfg%outer_plasma%model = 'kinetic_1d'
  cfg%outer_plasma%kinetic_closure = 'zhao_charge_driven'
  cfg%outer_plasma%zhao_branch = 'auto'
  cfg%outer_plasma%photoelectron_density_model = 'none'
  cfg%outer_plasma%return_model = 'kinetic_1d_profile_return'
  cfg%outer_plasma%debye_length = 1.0_dp
  cfg%outer_plasma%thermal_voltage = 2.0_dp
  cfg%sim%batch_duration = 1.0e-6_dp
  cfg%coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
  cfg%coupling%field_evolution_timescale = 1.0_dp
  cfg%coupling%outer_queue_enabled = .true.
  stats = sim_stats()
  stats%batches = 2_i32
  stats%last_rel_change = 0.0_dp
  electrostatic_diagnostics%split_periodic_active = .true.
  electrostatic_diagnostics%interface_potential = -0.25_dp
  electrostatic_diagnostics%last_outer_update_batch = 2_i32
  electrostatic_diagnostics%outer_applicability_status = 0_i32
  electrostatic_diagnostics%outer_nonlinear_iterations = 7_i32
  electrostatic_diagnostics%outer_nonlinear_residual = 2.5e-11_dp
  electrostatic_diagnostics%outer_integrated_charge_per_area = -4.0e-12_dp
  electrostatic_diagnostics%outer_electron_current_density = -2.0e-6_dp
  electrostatic_diagnostics%outer_ion_current_density = 1.5e-6_dp
  electrostatic_diagnostics%outer_photoelectron_current_density = 0.25e-6_dp
  electrostatic_diagnostics%outer_total_current_density = -0.25e-6_dp
  electrostatic_diagnostics%outer_kinetic_closure = 'zhao_charge_driven'
  electrostatic_diagnostics%outer_zhao_branch = 'A'
  electrostatic_diagnostics%outer_zhao_phi0 = 1.25_dp
  electrostatic_diagnostics%outer_zhao_phi_minimum = -0.4_dp
  electrostatic_diagnostics%outer_zhao_electron_density_infinity = 7.5e6_dp
  electrostatic_diagnostics%outer_photoelectron_source_scale = 0.625_dp
  electrostatic_diagnostics%outer_photoelectron_population_fraction = 0.35_dp
  electrostatic_diagnostics%outer_photoelectron_column_per_area = 2.0e9_dp
  electrostatic_diagnostics%outer_photoelectron_column_target_per_area = 2.0e9_dp
  electrostatic_diagnostics%outer_photoelectron_column_residual_per_area = 1.0e-3_dp
  electrostatic_diagnostics%outer_queue_event_count = 3_i64
  electrostatic_diagnostics%outer_queue_signed_charge = -4.5e-15_dp
  electrostatic_diagnostics%outer_queue_fingerprint = '0123456789ABCDEF'
  electrostatic_diagnostics%max_outer_flight_time = 1.25_dp
  electrostatic_diagnostics%max_frozen_field_ratio = 0.125_dp
  electrostatic_diagnostics%max_outer_energy_relative_error = 2.0e-12_dp
  electrostatic_diagnostics%outer_profile_z = [1.0_dp, 1.5_dp, 2.0_dp]
  electrostatic_diagnostics%outer_profile_potential = [-0.25_dp, -0.1_dp, -0.02_dp]
  electrostatic_diagnostics%outer_profile_field = [0.3_dp, 0.23_dp, 0.16_dp]
  electrostatic_diagnostics%outer_profile_charge_density = [-1.0e-12_dp, -0.4e-12_dp, -0.08e-12_dp]
  call write_result_files(out_dir, mesh, stats, cfg, electrostatic_diagnostics=electrostatic_diagnostics)
  call write_rng_state_file(out_dir)
  open (newunit=profile_unit, file=out_dir//'/outer_plasma_profile.csv', status='old', action='read', iostat=profile_ios)
  call assert_equal_i32(int(profile_ios, i32), 0_i32, 'schema v4 outer profile should open')
  read (profile_unit, '(A)', iostat=profile_ios) profile_header
  close (profile_unit)
  call assert_true( &
    trim(profile_header) == 'point,z_m,potential_V,field_V_m,charge_density_C_m3', &
    'schema v4 outer profile header mismatch' &
    )
  call load_restart_checkpoint(out_dir, mesh, stats, has_restart, app=cfg, electrostatic_state=electrostatic_state)
  call assert_true(has_restart, 'kinetic profile checkpoint should load')
  call assert_true(electrostatic_state%outer_ready, 'kinetic outer state should be ready')
  call assert_allclose_1d( &
    electrostatic_state%outer_profile_z, [1.0_dp, 1.5_dp, 2.0_dp], 1.0e-15_dp, &
    'kinetic restart grid mismatch' &
    )
  call assert_allclose_1d( &
    electrostatic_state%outer_profile_potential, [-0.25_dp, -0.1_dp, -0.02_dp], 1.0e-15_dp, &
    'kinetic restart potential mismatch' &
    )
  call assert_allclose_1d( &
    electrostatic_state%outer_profile_field, [0.3_dp, 0.23_dp, 0.16_dp], 1.0e-15_dp, &
    'kinetic restart field mismatch' &
    )
  call assert_allclose_1d( &
    electrostatic_state%outer_profile_charge_density, [-1.0e-12_dp, -0.4e-12_dp, -0.08e-12_dp], 1.0e-24_dp, &
    'kinetic restart charge-density mismatch' &
    )
  call assert_equal_i32(electrostatic_state%outer_nonlinear_iterations, 7_i32, &
                        'kinetic restart iteration mismatch')
  call assert_close_dp(electrostatic_state%outer_nonlinear_residual, 2.5e-11_dp, 1.0e-24_dp, &
                       'kinetic restart residual mismatch')
  call assert_close_dp(electrostatic_state%outer_integrated_charge_per_area, -4.0e-12_dp, 1.0e-24_dp, &
                       'kinetic restart integrated charge mismatch')
  call assert_close_dp(electrostatic_state%outer_electron_current_density, -2.0e-6_dp, 1.0e-18_dp, &
                       'kinetic restart electron current mismatch')
  call assert_close_dp(electrostatic_state%outer_ion_current_density, 1.5e-6_dp, 1.0e-18_dp, &
                       'kinetic restart ion current mismatch')
  call assert_close_dp(electrostatic_state%outer_photoelectron_current_density, 0.25e-6_dp, 1.0e-18_dp, &
                       'kinetic restart photoelectron current mismatch')
  call assert_close_dp(electrostatic_state%outer_total_current_density, -0.25e-6_dp, 1.0e-18_dp, &
                       'kinetic restart total current mismatch')
  call assert_true(electrostatic_state%outer_zhao_branch == 'A', 'Zhao restart branch mismatch')
  call assert_close_dp(electrostatic_state%outer_zhao_phi0, 1.25_dp, 1.0e-15_dp, &
                       'Zhao restart phi0 mismatch')
  call assert_close_dp(electrostatic_state%outer_zhao_phi_minimum, -0.4_dp, 1.0e-15_dp, &
                       'Zhao restart minimum mismatch')
  call assert_close_dp(electrostatic_state%outer_zhao_electron_density_infinity, 7.5e6_dp, 1.0e-8_dp, &
                       'Zhao restart electron density mismatch')
  call assert_true(electrostatic_state%outer_zhao_source_scale_complete, &
                   'resolved Zhao source-scale restart state should be complete')
  call assert_close_dp(electrostatic_state%outer_photoelectron_source_scale, 0.625_dp, 1.0e-15_dp, &
                       'resolved Zhao source-scale restart mismatch')
  call assert_true(electrostatic_state%outer_zhao_transient_state_complete, &
                   'transient Zhao restart state should be complete')
  call assert_close_dp(electrostatic_state%outer_photoelectron_population_fraction, 0.35_dp, 1.0e-15_dp, &
                       'transient Zhao population mismatch')
  call assert_close_dp(electrostatic_state%outer_photoelectron_column_per_area, 2.0e9_dp, 1.0e-6_dp, &
                       'transient Zhao column mismatch')
  call assert_close_dp(electrostatic_state%outer_photoelectron_column_target_per_area, 2.0e9_dp, 1.0e-6_dp, &
                       'transient Zhao target column mismatch')
  call assert_close_dp(electrostatic_state%outer_photoelectron_column_residual_per_area, 1.0e-3_dp, 1.0e-15_dp, &
                       'transient Zhao column residual mismatch')
  call assert_true(electrostatic_state%outer_queue_inventory_complete, &
                   'outer queue restart inventory should be complete')
  call assert_equal_i64(electrostatic_state%outer_queue_event_count, 3_i64, &
                        'outer queue restart event count mismatch')
  call assert_close_dp(electrostatic_state%outer_queue_signed_charge, -4.5e-15_dp, 1.0e-30_dp, &
                       'outer queue restart signed charge mismatch')
  call assert_true(electrostatic_state%outer_queue_fingerprint == '0123456789ABCDEF', &
                   'outer queue restart fingerprint mismatch')
  call assert_true(electrostatic_state%outer_max_diagnostics_complete, &
                   'schema v4 cumulative outer diagnostics should be complete')
  call assert_close_dp(electrostatic_state%max_outer_flight_time, 1.25_dp, 1.0e-15_dp, &
                       'cumulative outer flight-time mismatch')
  call assert_close_dp(electrostatic_state%max_frozen_field_ratio, 0.125_dp, 1.0e-15_dp, &
                       'cumulative frozen-field ratio mismatch')
  call assert_close_dp(electrostatic_state%max_outer_energy_relative_error, 2.0e-12_dp, 1.0e-24_dp, &
                       'cumulative outer energy error mismatch')
  call test_end()

  call test_begin('implicit_Zhao_restart_requires_resolved_source_scale')
  call build_test_mesh(mesh)
  call configure_implicit_zhao_restart(implicit_cfg)
  stats = sim_stats()
  stats%batches = 1_i32
  stats%last_rel_change = 0.0_dp
  electrostatic_diagnostics = electrostatic_diagnostics_type()
  electrostatic_diagnostics%split_periodic_active = .true.
  electrostatic_diagnostics%interface_potential = 0.75_dp
  electrostatic_diagnostics%interface_field = 0.2_dp
  electrostatic_diagnostics%last_outer_update_batch = 1_i32
  electrostatic_diagnostics%outer_applicability_status = 0_i32
  electrostatic_diagnostics%outer_nonlinear_iterations = 4_i32
  electrostatic_diagnostics%outer_nonlinear_residual = 1.0e-12_dp
  electrostatic_diagnostics%outer_debye_length = 1.0_dp
  electrostatic_diagnostics%outer_integrated_charge_per_area = -1.0e-12_dp
  electrostatic_diagnostics%outer_kinetic_closure = 'zhao_charge_driven'
  electrostatic_diagnostics%outer_zhao_branch = 'A'
  electrostatic_diagnostics%outer_zhao_phi0 = 0.75_dp
  electrostatic_diagnostics%outer_zhao_phi_minimum = -0.25_dp
  electrostatic_diagnostics%outer_zhao_electron_density_infinity = 8.7e6_dp
  electrostatic_diagnostics%outer_photoelectron_source_scale = 0.42_dp
  electrostatic_diagnostics%outer_profile_z = [1.0_dp, 1.5_dp, 2.0_dp]
  electrostatic_diagnostics%outer_profile_potential = [0.75_dp, -0.25_dp, 0.0_dp]
  electrostatic_diagnostics%outer_profile_field = [0.2_dp, 0.0_dp, 0.0_dp]
  electrostatic_diagnostics%outer_profile_charge_density = [-1.0e-12_dp, -0.5e-12_dp, 0.0_dp]
  call write_result_files( &
    implicit_out_dir, mesh, stats, implicit_cfg, electrostatic_diagnostics=electrostatic_diagnostics &
    )
  call write_rng_state_file(implicit_out_dir)
  call load_restart_checkpoint( &
    implicit_out_dir, mesh, stats, has_restart, app=implicit_cfg, electrostatic_state=electrostatic_state &
    )
  call assert_true(has_restart .and. electrostatic_state%outer_zhao_source_scale_complete, &
                   'implicit Zhao checkpoint did not restore its resolved source scale')
  call assert_close_dp(electrostatic_state%outer_photoelectron_source_scale, 0.42_dp, 1.0e-15_dp, &
                       'implicit Zhao resolved source scale changed on restart')

  call remove_summary_key( &
    implicit_out_dir, 'outer_plasma_photoelectron_source_scale_resolved' &
    )
  call get_command_argument(0, executable_path)
  probe_command = trim(executable_path)//' missing_implicit_zhao_source_scale_probe > '// &
                  missing_scale_probe_log//' 2>&1'
  call execute_command_line( &
    trim(probe_command), wait=.true., exitstat=child_exit_status, cmdstat=child_cmd_status &
    )
  call assert_equal_i32(int(child_cmd_status, i32), 0_i32, &
                        'failed to launch missing implicit-Zhao source-scale probe')
  call assert_true(child_exit_status /= 0, &
                   'implicit Zhao restart without a resolved source scale must fail closed')
  saw_missing_scale_error = .false.
  open (newunit=probe_unit, file=missing_scale_probe_log, status='old', action='read', iostat=probe_ios)
  if (probe_ios /= 0) error stop 'failed to read missing implicit-Zhao source-scale probe output'
  do
    read (probe_unit, '(A)', iostat=probe_ios) probe_line
    if (probe_ios /= 0) exit
    saw_missing_scale_error = saw_missing_scale_error .or. &
                              index(probe_line, 'resolved implicit-Zhao photoelectron source scale') > 0
  end do
  close (probe_unit)
  call assert_true(saw_missing_scale_error, &
                   'implicit Zhao restart rejection lost its resolved-source-scale reason')
  call cleanup_restart_fixture(implicit_out_dir)
  call delete_file_if_exists(missing_scale_probe_log)
  call test_end()

  call test_begin('dynamic_Zhao_two_batches_match_checkpoint_split')
  call test_dynamic_zhao_split_restart()
  call test_end()

  call test_begin('schema_v2_outer_profile_migrates_by_forcing_refresh')
  call downgrade_outer_checkpoint_to_v2(out_dir)
  call load_restart_checkpoint(out_dir, mesh, stats, has_restart, app=cfg, electrostatic_state=electrostatic_state)
  call assert_true(has_restart, 'schema v2 kinetic profile checkpoint should load')
  call assert_equal_i32(electrostatic_state%checkpoint_schema_version, 2_i32, &
                        'schema v2 version should be preserved')
  call assert_true(electrostatic_state%outer_ready, 'schema v2 outer state should remain available as an initial guess')
  call assert_true(.not. electrostatic_state%outer_profile_complete, &
                   'schema v2 profile must force a complete root refresh')
  call assert_true(.not. electrostatic_state%outer_max_diagnostics_complete, &
                   'schema v2 must not claim cumulative outer-diagnostic restart coverage')
  call assert_allclose_1d( &
    electrostatic_state%outer_profile_potential, [-0.25_dp, -0.1_dp, -0.02_dp], 1.0e-15_dp, &
    'schema v2 migration potential mismatch' &
    )
  call test_end()

  call test_begin('schema_v2_contract_and_ledger_restore')
  call default_app_config(cfg)
  cfg%n_particle_species = 2_i32
  cfg%particle_species(1) = species_from_defaults()
  cfg%particle_species(1)%species_key = 'electron'
  cfg%particle_species(2) = species_from_defaults()
  cfg%particle_species(2)%species_key = 'ion'
  stats%processed_particles = 10_i64
  stats%absorbed = 7_i64
  stats%escaped = 3_i64
  stats%batches = 2_i32
  stats%last_rel_change = 1.0e-3_dp
  call ledger%init(2_i32)
  call ledger%reset(2_i32)
  ledger%surface_charge_before = 1.0_dp
  ledger%surface_charge_after = 2.0_dp
  ledger%injected_from_remote = [-3.0_dp, 4.0_dp]
  ledger%injected_count = [3_i64, 4_i64]
  mesh%q_elem = [1.0e-12_dp, -2.0e-12_dp]
  call write_result_files(out_dir, mesh, stats, cfg, charge_ledger=ledger)
  call write_rng_state_file(out_dir)
  call write_macro_residuals_file(out_dir, state)

  call build_test_mesh(mesh)
  call load_restart_checkpoint( &
    out_dir, mesh, stats, has_restart, state, app=cfg, charge_ledger=restored_ledger &
    )
  call assert_true(has_restart, 'schema v2 checkpoint should load')
  call assert_equal_i32(restored_ledger%batch_count, 2_i32, 'ledger batch count mismatch')
  call assert_close_dp(restored_ledger%surface_charge_before, 1.0_dp, 1.0e-12_dp, 'ledger stock mismatch')
  call assert_allclose_1d( &
    restored_ledger%injected_from_remote, [-3.0_dp, 4.0_dp], 1.0e-12_dp, 'ledger flux restore mismatch' &
    )

  cfg_changed = cfg
  cfg_changed%sim%dt = 2.0_dp*cfg%sim%dt
  call validate_restart_contract( &
    out_dir//'/summary.txt', mesh, cfg_changed, contract_status, contract_message &
    )
  call assert_equal_i32(contract_status, restart_contract_mismatch, 'model mismatch should be rejected')
  call validate_restart_contract(out_dir//'/summary.txt', mesh, cfg, contract_status, contract_message)
  call assert_equal_i32(contract_status, restart_contract_ok, 'matching contract should be accepted')
  call test_end()

  call test_begin('write_checkpoint')
  call ensure_directory(out_dir)
  call write_summary_fixture(out_dir)
  call write_charges_fixture(out_dir)
  call write_rng_state_file(out_dir)
  allocate (state%macro_residual(2))
  state%macro_residual = [0.25d0, 0.75d0]
  call write_macro_residuals_file(out_dir, state)
  inquire (file=trim(out_dir)//'/rng_state.txt', exist=exists)
  call assert_true(exists, 'rng_state.txt should be created')
  inquire (file=trim(out_dir)//'/macro_residuals.csv', exist=exists)
  call assert_true(exists, 'macro_residuals.csv should be created')
  call test_end()

  call test_begin('ranked_paths')
  rng_rank_path = restart_rng_state_path(out_dir, mpi_rank=1_i32, mpi_size=4_i32)
  residual_global_path = restart_macro_residual_path(out_dir, mpi_rank=1_i32, mpi_size=4_i32)
  call assert_true(trim(rng_rank_path) == trim(out_dir)//'/rng_state_rank00001.txt', 'ranked rng path mismatch')
  call assert_true( &
    trim(residual_global_path) == trim(out_dir)//'/macro_residuals.csv', 'MPI residual path must remain global' &
    )

  call write_rng_state_file(out_dir, mpi_rank=1_i32, mpi_size=4_i32)
  state%macro_residual = [0.5d0, 0.5d0]
  call write_macro_residuals_file(out_dir, state, mpi_rank=1_i32, mpi_size=4_i32)
  inquire (file=trim(rng_rank_path), exist=exists)
  call assert_true(exists, 'rng_state_rank00001.txt should be created')
  inquire (file=trim(out_dir)//'/macro_residuals_rank00001.csv', exist=exists)
  call assert_true(.not. exists, 'non-root must not create a ranked macro residual file')
  call test_end()

  call test_begin('load_checkpoint')
  call build_test_mesh(mesh)
  mesh%q_elem = [3.0d0, 4.0d0]
  state%macro_residual = 0.0d0
  call load_restart_checkpoint(out_dir, mesh, stats, has_restart, state)
  call assert_true(has_restart, 'complete checkpoint should be detected')
  call assert_equal_i64(stats%processed_particles, 10_i64, 'processed_particles mismatch')
  call assert_equal_i64(stats%absorbed, 7_i64, 'absorbed mismatch')
  call assert_equal_i64(stats%escaped, 3_i64, 'escaped mismatch')
  call assert_equal_i32(stats%batches, 2_i32, 'batches mismatch')
  call assert_equal_i64(stats%escaped_boundary, 1_i64, 'escaped_boundary mismatch')
  call assert_equal_i64(stats%survived_max_step, 2_i64, 'survived_max_step mismatch')
  call assert_equal_i64( &
    stats%multiple_box_events_soft_discarded, 4_i64, 'soft discarded count mismatch' &
    )
  call assert_equal_i32( &
    stats%adaptive_nonzero_mode_omp_threads, 6_i32, 'adaptive OpenMP thread count mismatch' &
    )
  call assert_close_dp( &
    stats%multiple_box_events_soft_discarded_abs_charge, 2.5d-15, 1.0d-27, &
    'soft discarded absolute charge mismatch' &
    )
  call assert_close_dp(stats%last_rel_change, 1.0d-3, 1.0d-12, 'last_rel_change mismatch')
  call assert_allclose_1d(mesh%q_elem, [1.0d-12, -2.0d-12], 1.0d-24, 'charge restore mismatch')
  call assert_allclose_1d(state%macro_residual, [0.25d0, 0.75d0], 1.0d-12, 'macro residual restore mismatch')
  call test_end()

  call delete_file_if_exists(out_dir//'/summary.txt')
  call delete_file_if_exists(out_dir//'/charges.csv')
  call delete_file_if_exists(out_dir//'/rng_state.txt')
  call delete_file_if_exists(out_dir//'/macro_residuals.csv')
  call delete_file_if_exists(out_dir//'/rng_state_rank00001.txt')
  call delete_file_if_exists(out_dir//'/macro_residuals_rank00001.csv')
  call delete_file_if_exists(out_dir//'/mesh_triangles.csv')
  call delete_file_if_exists(out_dir//'/mesh_sources.csv')
  call delete_file_if_exists(out_dir//'/charge_ledger.csv')
  call delete_file_if_exists(out_dir//'/outer_plasma_profile.csv')
  call remove_empty_directory(out_dir)
  call cleanup_restart_fixture(implicit_out_dir)
  call cleanup_restart_fixture(dynamic_split_out_dir)
  call delete_file_if_exists(missing_scale_probe_log)

  call test_summary()

contains

  subroutine test_dynamic_zhao_split_restart()
    type(mesh_type) :: straight_mesh, first_mesh, resumed_mesh
    type(app_config) :: straight_cfg, first_cfg, resumed_cfg
    type(sim_stats) :: straight_stats, first_stats, loaded_stats, resumed_stats
    type(injection_state) :: straight_injection, first_injection, resumed_injection
    type(charge_ledger_type) :: straight_ledger, first_ledger, resumed_ledger
    type(electrostatic_diagnostics_type) :: straight_diagnostics, first_diagnostics
    type(electrostatic_restart_state_type) :: straight_outer, first_outer, resumed_outer
    integer, allocatable :: straight_rng(:), resumed_rng(:)
    integer :: rng_size
    logical :: loaded
    real(dp) :: charge_tolerance, potential_tolerance, field_tolerance
    real(dp) :: coordinate_tolerance, density_tolerance, ledger_tolerance

    call configure_dynamic_type_a_fixture(straight_cfg)
    straight_cfg%sim%batch_count = 2_i32
    call build_mesh_from_config(straight_cfg, straight_mesh)
    call prepare_periodic2_collision_mesh(straight_mesh, straight_cfg%sim)
    allocate (straight_injection%macro_residual(straight_cfg%n_particle_species))
    straight_injection%macro_residual = 0.0_dp
    call seed_particles_from_config(straight_cfg)
    call run_absorption_insulator( &
      straight_mesh, straight_cfg, straight_stats, inject_state=straight_injection, &
      charge_ledger=straight_ledger, electrostatic_diagnostics=straight_diagnostics, &
      electrostatic_restart_state=straight_outer &
      )
    call random_seed(size=rng_size)
    allocate (straight_rng(rng_size), resumed_rng(rng_size))
    call random_seed(get=straight_rng)

    call assert_equal_i32(straight_stats%batches, 2_i32, &
                          'straight dynamic-Zhao fixture did not complete two batches')
    call assert_true(straight_outer%outer_ready .and. straight_outer%outer_profile_complete, &
                     'straight dynamic-Zhao run did not export a complete outer profile')
    call assert_true(straight_outer%outer_zhao_state_complete .and. &
                     straight_outer%outer_zhao_source_scale_complete, &
                     'straight dynamic-Zhao run did not export a complete resolved Zhao state')
    call assert_true(straight_outer%outer_zhao_branch == 'A', &
                     'strong-PE steady start did not stay on Zhao Type A')
    call assert_true(straight_ledger%emitted_count(3) > 0_i64 .and. &
                     straight_ledger%interface_outward_gross(3) < 0.0_dp, &
                     'strong-PE fixture did not exercise measured photoelectron interface transport')
    call assert_true( &
      minval(straight_outer%outer_profile_potential) < &
      min(straight_outer%outer_interface_potential, straight_outer%outer_infinity_potential), &
      'strong-PE Zhao Type A profile lost its virtual cathode' &
      )

    first_cfg = straight_cfg
    first_cfg%sim%batch_count = 1_i32
    call build_mesh_from_config(first_cfg, first_mesh)
    call prepare_periodic2_collision_mesh(first_mesh, first_cfg%sim)
    allocate (first_injection%macro_residual(first_cfg%n_particle_species))
    first_injection%macro_residual = 0.0_dp
    call seed_particles_from_config(first_cfg)
    call run_absorption_insulator( &
      first_mesh, first_cfg, first_stats, inject_state=first_injection, &
      charge_ledger=first_ledger, electrostatic_diagnostics=first_diagnostics, &
      electrostatic_restart_state=first_outer &
      )
    call assert_equal_i32(first_stats%batches, 1_i32, &
                          'split dynamic-Zhao fixture did not complete its first batch')
    call assert_true(first_outer%outer_zhao_branch == 'A', &
                     'split dynamic-Zhao first batch left Zhao Type A')
    call assert_true(first_outer%outer_zhao_state_complete .and. &
                     first_outer%outer_zhao_source_scale_complete, &
                     'split dynamic-Zhao first batch did not export a complete resolved Zhao state')

    call write_result_files( &
      dynamic_split_out_dir, first_mesh, first_stats, first_cfg, &
      charge_ledger=first_ledger, electrostatic_diagnostics=first_diagnostics &
      )
    call write_rng_state_file(dynamic_split_out_dir)
    call write_macro_residuals_file(dynamic_split_out_dir, first_injection)

    resumed_cfg = straight_cfg
    call build_mesh_from_config(resumed_cfg, resumed_mesh)
    call prepare_periodic2_collision_mesh(resumed_mesh, resumed_cfg%sim)
    allocate (resumed_injection%macro_residual(resumed_cfg%n_particle_species))
    resumed_injection%macro_residual = 0.0_dp
    call seed_particles_from_config(resumed_cfg)
    call load_restart_checkpoint( &
      dynamic_split_out_dir, resumed_mesh, loaded_stats, loaded, resumed_injection, &
      app=resumed_cfg, charge_ledger=resumed_ledger, electrostatic_state=resumed_outer &
      )
    call assert_true(loaded, 'dynamic-Zhao split checkpoint was not restored')
    call assert_equal_i32(loaded_stats%batches, 1_i32, &
                          'dynamic-Zhao split checkpoint restored the wrong batch count')
    call assert_true(resumed_outer%outer_zhao_branch == 'A', &
                     'dynamic-Zhao split checkpoint lost its committed Type-A branch')
    call assert_true(resumed_outer%outer_profile_complete .and. &
                     resumed_outer%outer_zhao_source_scale_complete, &
                     'dynamic-Zhao split checkpoint lost its complete outer state')
    call assert_close_dp( &
      resumed_outer%outer_photoelectron_source_scale, first_outer%outer_photoelectron_source_scale, &
      1.0e-14_dp*max(1.0_dp, abs(first_outer%outer_photoelectron_source_scale)), &
      'dynamic-Zhao split checkpoint changed its measured source scale' &
      )

    call run_absorption_insulator( &
      resumed_mesh, resumed_cfg, resumed_stats, initial_stats=loaded_stats, &
      inject_state=resumed_injection, charge_ledger=resumed_ledger, &
      electrostatic_restart_state=resumed_outer &
      )
    call random_seed(get=resumed_rng)

    charge_tolerance = 1.0e-10_dp*max( &
                       maxval(abs(straight_mesh%q_elem)), &
                       abs(straight_ledger%surface_charge_after), tiny(1.0_dp) &
                       )
    potential_tolerance = 1.0e-9_dp*max( &
                          maxval(abs(straight_outer%outer_profile_potential)), 1.0_dp &
                          )
    field_tolerance = 1.0e-9_dp*max( &
                      maxval(abs(straight_outer%outer_profile_field)), 1.0_dp &
                      )
    coordinate_tolerance = 1.0e-12_dp*max( &
                           maxval(abs(straight_outer%outer_profile_z)), 1.0_dp &
                           )
    density_tolerance = 1.0e-9_dp*max( &
                        maxval(abs(straight_outer%outer_profile_charge_density)), tiny(1.0_dp) &
                        )
    ledger_tolerance = 1.0e-10_dp*max( &
                       abs(straight_ledger%surface_charge_before), &
                       abs(straight_ledger%surface_charge_after), &
                       maxval(abs(straight_ledger%emitted_from_surface)), tiny(1.0_dp) &
                       )

    call assert_equal_i32(resumed_stats%batches, straight_stats%batches, &
                          'dynamic-Zhao split run completed a different batch count')
    call assert_equal_i64(resumed_stats%processed_particles, straight_stats%processed_particles, &
                          'dynamic-Zhao split run changed processed-particle count')
    call assert_equal_i64(resumed_stats%absorbed, straight_stats%absorbed, &
                          'dynamic-Zhao split run changed absorbed-particle count')
    call assert_equal_i64(resumed_stats%escaped, straight_stats%escaped, &
                          'dynamic-Zhao split run changed escaped-particle count')
    call assert_equal_i64(resumed_stats%escaped_boundary, straight_stats%escaped_boundary, &
                          'dynamic-Zhao split run changed boundary-escape count')
    call assert_equal_i64(resumed_stats%survived_max_step, straight_stats%survived_max_step, &
                          'dynamic-Zhao split run changed max-step survivor count')
    call assert_equal_i64( &
      resumed_stats%multiple_box_events_soft_discarded, &
      straight_stats%multiple_box_events_soft_discarded, &
      'dynamic-Zhao split run changed soft-discard count' &
      )
    call assert_allclose_1d( &
      resumed_mesh%q_elem, straight_mesh%q_elem, charge_tolerance, &
      'dynamic-Zhao split run changed element charge' &
      )
    call assert_true(resumed_outer%outer_zhao_branch == straight_outer%outer_zhao_branch, &
                     'dynamic-Zhao split run changed the committed branch')
    call assert_close_dp(resumed_outer%outer_interface_potential, straight_outer%outer_interface_potential, &
                         potential_tolerance, 'dynamic-Zhao split run changed interface potential')
    call assert_close_dp(resumed_outer%outer_interface_field, straight_outer%outer_interface_field, &
                         field_tolerance, 'dynamic-Zhao split run changed interface field')
    call assert_close_dp( &
      resumed_outer%outer_photoelectron_source_scale, straight_outer%outer_photoelectron_source_scale, &
      1.0e-10_dp*max(1.0_dp, abs(straight_outer%outer_photoelectron_source_scale)), &
      'dynamic-Zhao split run changed the measured source scale' &
      )
    call assert_close_dp(resumed_outer%outer_zhao_phi0, straight_outer%outer_zhao_phi0, &
                         potential_tolerance, 'dynamic-Zhao split run changed phi0')
    call assert_close_dp( &
      resumed_outer%outer_zhao_phi_minimum, straight_outer%outer_zhao_phi_minimum, &
      potential_tolerance, 'dynamic-Zhao split run changed the virtual-cathode potential' &
      )
    call assert_close_dp( &
      resumed_outer%outer_zhao_electron_density_infinity, &
      straight_outer%outer_zhao_electron_density_infinity, &
      1.0e-10_dp*max(1.0_dp, abs(straight_outer%outer_zhao_electron_density_infinity)), &
      'dynamic-Zhao split run changed the resolved infinity electron density' &
      )
    call assert_close_dp( &
      resumed_outer%outer_infinity_potential, straight_outer%outer_infinity_potential, &
      potential_tolerance, 'dynamic-Zhao split run changed the infinity gauge' &
      )
    call assert_close_dp( &
      resumed_outer%outer_integrated_charge_per_area, &
      straight_outer%outer_integrated_charge_per_area, &
      density_tolerance*max(1.0_dp, maxval(abs(straight_outer%outer_profile_z))), &
      'dynamic-Zhao split run changed the integrated outer charge' &
      )
    call assert_allclose_1d( &
      resumed_outer%outer_profile_z, straight_outer%outer_profile_z, coordinate_tolerance, &
      'dynamic-Zhao split run changed outer-profile coordinates' &
      )
    call assert_allclose_1d( &
      resumed_outer%outer_profile_potential, straight_outer%outer_profile_potential, &
      potential_tolerance, 'dynamic-Zhao split run changed outer potential' &
      )
    call assert_allclose_1d( &
      resumed_outer%outer_profile_field, straight_outer%outer_profile_field, &
      field_tolerance, 'dynamic-Zhao split run changed outer field' &
      )
    call assert_allclose_1d( &
      resumed_outer%outer_profile_charge_density, straight_outer%outer_profile_charge_density, &
      density_tolerance, 'dynamic-Zhao split run changed outer charge density' &
      )
    call assert_close_dp( &
      resumed_ledger%surface_charge_before, straight_ledger%surface_charge_before, &
      ledger_tolerance, 'dynamic-Zhao split run changed initial ledger stock' &
      )
    call assert_close_dp( &
      resumed_ledger%surface_charge_after, straight_ledger%surface_charge_after, &
      ledger_tolerance, 'dynamic-Zhao split run changed final ledger stock' &
      )
    call assert_equal_i32(resumed_ledger%batch_count, straight_ledger%batch_count, &
                          'dynamic-Zhao split run changed ledger batch count')
    call assert_close_dp( &
      resumed_ledger%local_flight_charge_before, straight_ledger%local_flight_charge_before, &
      ledger_tolerance, 'dynamic-Zhao split run changed initial local-flight stock' &
      )
    call assert_close_dp( &
      resumed_ledger%local_flight_charge_after, straight_ledger%local_flight_charge_after, &
      ledger_tolerance, 'dynamic-Zhao split run changed final local-flight stock' &
      )
    call assert_close_dp( &
      resumed_ledger%outer_flight_charge_before, straight_ledger%outer_flight_charge_before, &
      ledger_tolerance, 'dynamic-Zhao split run changed initial outer-flight stock' &
      )
    call assert_close_dp( &
      resumed_ledger%outer_flight_charge_after, straight_ledger%outer_flight_charge_after, &
      ledger_tolerance, 'dynamic-Zhao split run changed final outer-flight stock' &
      )
    call assert_close_dp( &
      resumed_ledger%unresolved_stock_before, straight_ledger%unresolved_stock_before, &
      ledger_tolerance, 'dynamic-Zhao split run changed initial unresolved stock' &
      )
    call assert_close_dp( &
      resumed_ledger%unresolved_stock_after, straight_ledger%unresolved_stock_after, &
      ledger_tolerance, 'dynamic-Zhao split run changed final unresolved stock' &
      )
    call assert_allclose_1d( &
      resumed_ledger%injected_from_remote, straight_ledger%injected_from_remote, &
      ledger_tolerance, 'dynamic-Zhao split run changed injected charge ledger' &
      )
    call assert_allclose_1d( &
      resumed_ledger%emitted_from_surface, straight_ledger%emitted_from_surface, &
      ledger_tolerance, 'dynamic-Zhao split run changed emitted charge ledger' &
      )
    call assert_allclose_1d( &
      resumed_ledger%absorbed_on_surface, straight_ledger%absorbed_on_surface, &
      ledger_tolerance, 'dynamic-Zhao split run changed absorbed charge ledger' &
      )
    call assert_allclose_1d( &
      resumed_ledger%escaped_to_infinity, straight_ledger%escaped_to_infinity, &
      ledger_tolerance, 'dynamic-Zhao split run changed escaped charge ledger' &
      )
    call assert_allclose_1d( &
      resumed_ledger%discarded_unresolved, straight_ledger%discarded_unresolved, &
      ledger_tolerance, 'dynamic-Zhao split run changed unresolved charge ledger' &
      )
    call assert_allclose_1d( &
      resumed_ledger%interface_outward_gross, straight_ledger%interface_outward_gross, &
      ledger_tolerance, 'dynamic-Zhao split run changed outward interface ledger' &
      )
    call assert_allclose_1d( &
      resumed_ledger%interface_returned_gross, straight_ledger%interface_returned_gross, &
      ledger_tolerance, 'dynamic-Zhao split run changed returned interface ledger' &
      )
    call assert_true(all(resumed_ledger%injected_count == straight_ledger%injected_count) .and. &
                     all(resumed_ledger%emitted_count == straight_ledger%emitted_count) .and. &
                     all(resumed_ledger%absorbed_count == straight_ledger%absorbed_count) .and. &
                     all(resumed_ledger%escaped_count == straight_ledger%escaped_count) .and. &
                     all(resumed_ledger%discarded_unresolved_count == &
                         straight_ledger%discarded_unresolved_count), &
                     'dynamic-Zhao split run changed charge-ledger event counts')
    call assert_true(all(resumed_rng == straight_rng), &
                     'dynamic-Zhao split run changed the final RNG state')
    call assert_allclose_1d( &
      resumed_injection%macro_residual, straight_injection%macro_residual, 1.0e-15_dp, &
      'dynamic-Zhao split run changed macro-particle residuals' &
      )

    call cleanup_restart_fixture(dynamic_split_out_dir)
  end subroutine test_dynamic_zhao_split_restart

  subroutine configure_dynamic_type_a_fixture(resolved)
    type(app_config), intent(out) :: resolved
    real(dp), parameter :: electron_mass = 9.1093837015e-31_dp
    real(dp), parameter :: proton_mass = 1.67262192369e-27_dp
    real(dp), parameter :: box_width = 9.899494936611664e-5_dp
    real(dp), parameter :: interface_z = 3.0e-5_dp
    real(dp), parameter :: support_z = 2.0e-6_dp
    real(dp), parameter :: ambient_density = 8.7e6_dp
    real(dp), parameter :: photoelectron_reference_density = 64.0e6_dp
    real(dp), parameter :: alpha_rad = 60.0_dp*pi/180.0_dp
    real(dp), parameter :: inward_drift = 4.68e5_dp*sin(alpha_rad)
    real(dp) :: photoelectron_current_density

    photoelectron_current_density = qe*photoelectron_reference_density*sin(alpha_rad)* &
                                    sqrt(2.0_dp*2.2_dp*qe/electron_mass)/(2.0_dp*sqrt(pi))

    call default_app_config(resolved)
    resolved%mesh_mode = 'template'
    resolved%sim%rng_seed = 97531_i32
    resolved%sim%dt = 2.0e-12_dp
    resolved%sim%batch_duration = 1.0e-7_dp
    resolved%sim%has_batch_duration = .true.
    resolved%sim%max_step = 4000_i32
    resolved%sim%q_floor = 1.0e-40_dp
    resolved%sim%field_solver = 'direct'
    resolved%sim%field_bc_mode = 'periodic2'
    resolved%sim%use_box = .true.
    resolved%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    resolved%sim%box_max = [box_width, box_width, interface_z]
    resolved%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    resolved%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    resolved%sim%reservoir_potential_model = 'none'
    resolved%sim%open_boundary_model = 'escape'
    resolved%sim%sheath_alpha_deg = 60.0_dp
    resolved%sim%sheath_photoelectron_ref_density_cm3 = &
      photoelectron_reference_density*1.0e-6_dp
    resolved%sim%sheath_electron_drift_mode = 'normal'
    resolved%sim%sheath_ion_drift_mode = 'normal'

    resolved%field%backend = 'direct'
    resolved%panel%kernel_id = 'triangle_p0_exact_direct'
    resolved%panel%surface_side_policy = 'per_element'
    resolved%periodic2%nonzero_mode_backend = 'panel_spectral_reference'
    resolved%periodic2%zero_mode_policy = 'exclude_k0'
    resolved%periodic2%lower_boundary_model = 'e_bottom_zero'
    resolved%periodic2%reference_mode_layers = 2_i32
    resolved%periodic2%panel_quadrature_order = 4_i32
    resolved%periodic2%interface_sample_n = 2_i32
    resolved%periodic2%interface_phi_tolerance = 1.0e6_dp
    resolved%periodic2%interface_field_tolerance = 1.0e6_dp

    resolved%outer_plasma%model = 'kinetic_1d'
    resolved%outer_plasma%kinetic_closure = 'zhao_charge_driven'
    resolved%outer_plasma%zhao_branch = 'auto'
    resolved%outer_plasma%photoelectron_source_scale = 1.0_dp
    resolved%outer_plasma%photoelectron_density_model = 'none'
    resolved%outer_plasma%return_model = 'kinetic_1d_profile_return'
    resolved%outer_plasma%interface_z = interface_z
    resolved%outer_plasma%debye_length = 1.0_dp
    resolved%outer_plasma%thermal_voltage = 10.0_dp
    resolved%coupling%update_mode = 'implicit_mean'
    resolved%coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
    resolved%coupling%steady_start_mode = 'zhao_floating'
    resolved%coupling%steady_start_mesh_id = 1_i32
    resolved%coupling%outer_update_stride = 1_i32
    resolved%coupling%field_evolution_timescale = 1.0_dp
    resolved%coupling%max_frozen_field_ratio = 1.0_dp
    resolved%coupling%outer_queue_enabled = .false.

    resolved%n_particle_species = 3_i32
    resolved%particle_species(1) = species_from_defaults()
    resolved%particle_species(1)%species_key = 'ambient_electron'
    resolved%particle_species(1)%source_mode = 'reservoir_face'
    resolved%particle_species(1)%inject_face = 'z_high'
    resolved%particle_species(1)%q_particle = -qe
    resolved%particle_species(1)%m_particle = electron_mass
    resolved%particle_species(1)%w_particle = 1.0e-3_dp
    resolved%particle_species(1)%has_w_particle = .true.
    resolved%particle_species(1)%number_density_m3 = ambient_density
    resolved%particle_species(1)%has_number_density_m3 = .true.
    resolved%particle_species(1)%temperature_ev = 12.0_dp
    resolved%particle_species(1)%has_temperature_ev = .true.
    resolved%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -inward_drift]
    resolved%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, interface_z]
    resolved%particle_species(1)%pos_high = [box_width, box_width, interface_z]

    resolved%particle_species(2) = species_from_defaults()
    resolved%particle_species(2)%species_key = 'ambient_proton'
    resolved%particle_species(2)%source_mode = 'reservoir_face'
    resolved%particle_species(2)%inject_face = 'z_high'
    resolved%particle_species(2)%q_particle = qe
    resolved%particle_species(2)%m_particle = proton_mass
    resolved%particle_species(2)%w_particle = 1.0e-3_dp
    resolved%particle_species(2)%has_w_particle = .true.
    resolved%particle_species(2)%number_density_m3 = ambient_density
    resolved%particle_species(2)%has_number_density_m3 = .true.
    resolved%particle_species(2)%temperature_ev = 0.1_dp
    resolved%particle_species(2)%has_temperature_ev = .true.
    resolved%particle_species(2)%drift_velocity = [0.0_dp, 0.0_dp, -inward_drift]
    resolved%particle_species(2)%pos_low = [0.0_dp, 0.0_dp, interface_z]
    resolved%particle_species(2)%pos_high = [box_width, box_width, interface_z]

    resolved%particle_species(3) = species_from_defaults()
    resolved%particle_species(3)%species_key = 'photoelectron'
    resolved%particle_species(3)%source_mode = 'photo_raycast'
    resolved%particle_species(3)%inject_face = 'z_high'
    resolved%particle_species(3)%q_particle = -qe
    resolved%particle_species(3)%m_particle = electron_mass
    resolved%particle_species(3)%temperature_ev = 2.2_dp
    resolved%particle_species(3)%has_temperature_ev = .true.
    resolved%particle_species(3)%emit_current_density_a_m2 = photoelectron_current_density
    resolved%particle_species(3)%rays_per_batch = 64_i32
    resolved%particle_species(3)%deposit_opposite_charge_on_emit = .true.
    resolved%particle_species(3)%has_deposit_opposite_charge_on_emit = .true.
    resolved%particle_species(3)%normal_drift_speed = 0.0_dp
    resolved%particle_species(3)%pos_low = [0.0_dp, 0.0_dp, interface_z]
    resolved%particle_species(3)%pos_high = [box_width, box_width, interface_z]
    resolved%particle_species(3)%ray_direction = [0.0_dp, 0.0_dp, -1.0_dp]
    resolved%particle_species(3)%has_ray_direction = .true.

    resolved%n_templates = 1_i32
    resolved%templates(1)%enabled = .true.
    resolved%templates(1)%kind = 'plane'
    resolved%templates(1)%surface_side_policy = 'normal_plus'
    resolved%templates(1)%size_x = box_width
    resolved%templates(1)%size_y = box_width
    resolved%templates(1)%nx = 2_i32
    resolved%templates(1)%ny = 2_i32
    resolved%templates(1)%center = [0.5_dp*box_width, 0.5_dp*box_width, support_z]
  end subroutine configure_dynamic_type_a_fixture

  subroutine configure_implicit_zhao_restart(resolved)
    type(app_config), intent(out) :: resolved

    call default_app_config(resolved)
    resolved%outer_plasma%model = 'kinetic_1d'
    resolved%outer_plasma%kinetic_closure = 'zhao_charge_driven'
    resolved%outer_plasma%zhao_branch = 'a'
    resolved%outer_plasma%photoelectron_density_model = 'none'
    resolved%outer_plasma%return_model = 'kinetic_1d_profile_return'
    resolved%outer_plasma%debye_length = 1.0_dp
    resolved%outer_plasma%thermal_voltage = 2.0_dp
    resolved%coupling%update_mode = 'implicit_mean'
    resolved%coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
    resolved%coupling%steady_start_mode = 'zhao_floating'
    resolved%coupling%outer_update_stride = 1_i32
    resolved%coupling%field_evolution_timescale = 1.0_dp
    resolved%coupling%outer_queue_enabled = .false.
  end subroutine configure_implicit_zhao_restart

  subroutine run_missing_implicit_zhao_source_scale_probe()
    type(mesh_type) :: probe_mesh
    type(sim_stats) :: probe_stats
    type(app_config) :: probe_cfg
    type(electrostatic_restart_state_type) :: probe_electrostatic_state
    logical :: probe_has_restart

    call build_test_mesh(probe_mesh)
    call configure_implicit_zhao_restart(probe_cfg)
    call load_restart_checkpoint( &
      implicit_out_dir, probe_mesh, probe_stats, probe_has_restart, app=probe_cfg, &
      electrostatic_state=probe_electrostatic_state &
      )
  end subroutine run_missing_implicit_zhao_source_scale_probe

  subroutine remove_summary_key(dir_path, removed_key)
    character(len=*), intent(in) :: dir_path, removed_key
    character(len=1024) :: summary_path, summary_tmp, line
    integer :: input_unit, output_unit, ios, rename_status

    summary_path = trim(dir_path)//'/summary.txt'
    summary_tmp = trim(summary_path)//'.tmp'
    open (newunit=input_unit, file=trim(summary_path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open summary while removing a restart key'
    open (newunit=output_unit, file=trim(summary_tmp), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to create summary without a restart key'
    do
      read (input_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (index(line, trim(removed_key)//'=') /= 1) write (output_unit, '(a)') trim(line)
    end do
    close (input_unit)
    close (output_unit)
    call atomic_rename(summary_tmp, summary_path, rename_status)
    if (rename_status /= filesystem_success) error stop 'failed to publish summary without a restart key'
  end subroutine remove_summary_key

  subroutine cleanup_restart_fixture(dir_path)
    character(len=*), intent(in) :: dir_path

    call delete_file_if_exists(trim(dir_path)//'/summary.txt')
    call delete_file_if_exists(trim(dir_path)//'/summary.txt.tmp')
    call delete_file_if_exists(trim(dir_path)//'/charges.csv')
    call delete_file_if_exists(trim(dir_path)//'/rng_state.txt')
    call delete_file_if_exists(trim(dir_path)//'/macro_residuals.csv')
    call delete_file_if_exists(trim(dir_path)//'/charge_ledger.csv')
    call delete_file_if_exists(trim(dir_path)//'/mesh_triangles.csv')
    call delete_file_if_exists(trim(dir_path)//'/mesh_sources.csv')
    call delete_file_if_exists(trim(dir_path)//'/outer_plasma_profile.csv')
    call remove_empty_directory(dir_path)
  end subroutine cleanup_restart_fixture

  subroutine downgrade_outer_checkpoint_to_v2(dir_path)
    character(len=*), intent(in) :: dir_path
    character(len=1024) :: summary_path, summary_tmp, profile_path, profile_tmp
    character(len=1024) :: line
    integer :: input_unit, output_unit, ios, rename_status
    integer(i32) :: point
    real(dp) :: z, potential, field, charge_density

    summary_path = trim(dir_path)//'/summary.txt'
    summary_tmp = trim(summary_path)//'.v2tmp'
    open (newunit=input_unit, file=trim(summary_path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open schema v4 summary for migration fixture'
    open (newunit=output_unit, file=trim(summary_tmp), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to create schema v2 summary fixture'
    do
      read (input_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (index(line, 'checkpoint_schema_version=') == 1) then
        write (output_unit, '(a)') 'checkpoint_schema_version=2'
      else
        write (output_unit, '(a)') trim(line)
      end if
    end do
    close (input_unit)
    close (output_unit)
    call atomic_rename(summary_tmp, summary_path, rename_status)
    if (rename_status /= filesystem_success) error stop 'failed to publish schema v2 summary fixture'

    profile_path = trim(dir_path)//'/outer_plasma_profile.csv'
    profile_tmp = trim(profile_path)//'.v2tmp'
    open (newunit=input_unit, file=trim(profile_path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open schema v4 profile for migration fixture'
    open (newunit=output_unit, file=trim(profile_tmp), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to create schema v2 profile fixture'
    read (input_unit, '(A)', iostat=ios) line
    if (ios /= 0) error stop 'schema v4 profile header is missing'
    write (output_unit, '(a)') 'point,z_m,potential_V'
    do
      read (input_unit, *, iostat=ios) point, z, potential, field, charge_density
      if (ios /= 0) exit
      write (output_unit, '(i0,2(a,es24.16))') point, ',', z, ',', potential
    end do
    close (input_unit)
    close (output_unit)
    call atomic_rename(profile_tmp, profile_path, rename_status)
    if (rename_status /= filesystem_success) error stop 'failed to publish schema v2 profile fixture'
  end subroutine downgrade_outer_checkpoint_to_v2

  !> 2要素メッシュを初期化する。
  !! @param[out] mesh 初期化済みテスト用メッシュ。
  subroutine build_test_mesh(mesh)
    type(mesh_type), intent(out) :: mesh
    real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2)

    v0(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    v1(:, 1) = [1.0d0, 0.0d0, 0.0d0]
    v2(:, 1) = [0.0d0, 1.0d0, 0.0d0]
    v0(:, 2) = [0.0d0, 0.0d0, 1.0d0]
    v1(:, 2) = [1.0d0, 0.0d0, 1.0d0]
    v2(:, 2) = [0.0d0, 1.0d0, 1.0d0]
    call init_mesh(mesh, v0, v1, v2)
  end subroutine build_test_mesh

  !> `summary.txt` のフィクスチャを書き出す。
  !! @param[in] dir_path 出力先ディレクトリ。
  subroutine write_summary_fixture(dir_path)
    character(len=*), intent(in) :: dir_path
    integer :: u, ios

    open (newunit=u, file=trim(dir_path)//'/summary.txt', status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open summary fixture'
    write (u, '(a)') 'mesh_nelem=2'
    write (u, '(a)') 'processed_particles=10'
    write (u, '(a)') 'absorbed=7'
    write (u, '(a)') 'escaped=3'
    write (u, '(a)') 'batches=2'
    write (u, '(a)') 'escaped_boundary=1'
    write (u, '(a)') 'survived_max_step=2'
    write (u, '(a)') 'multiple_box_events_soft_discarded=4'
    write (u, '(a)') 'multiple_box_events_soft_discarded_abs_charge_C=2.5e-15'
    write (u, '(a)') 'adaptive_nonzero_mode_omp_threads=6'
    write (u, '(a)') 'last_rel_change=1.0e-3'
    close (u)
  end subroutine write_summary_fixture

  !> `charges.csv` のフィクスチャを書き出す。
  !! @param[in] dir_path 出力先ディレクトリ。
  subroutine write_charges_fixture(dir_path)
    character(len=*), intent(in) :: dir_path
    integer :: u, ios

    open (newunit=u, file=trim(dir_path)//'/charges.csv', status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open charges fixture'
    write (u, '(a)') 'elem_idx,charge_C'
    write (u, '(a)') '1,1.0e-12'
    write (u, '(a)') '2,-2.0e-12'
    close (u)
  end subroutine write_charges_fixture

end program test_restart

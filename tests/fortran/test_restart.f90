!> 再開用チェックポイントの保存/復元を検証するテスト。
program test_restart
  use bem_kinds, only: dp, i32, i64
  use bem_mesh, only: init_mesh
  use bem_restart, only: load_restart_checkpoint, write_rng_state_file, write_macro_residuals_file, &
                         restart_rng_state_path, restart_macro_residual_path, validate_restart_contract, &
                         restart_contract_ok, restart_contract_mismatch
  use bem_output_writer, only: write_result_files
  use bem_app_config, only: app_config, default_app_config, species_from_defaults
  use bem_charge_ledger, only: charge_ledger_type
  use bem_electrostatic_snapshot, only: electrostatic_diagnostics_type, electrostatic_restart_state_type
  use bem_filesystem, only: atomic_rename, filesystem_success
  use bem_types, only: mesh_type, sim_stats, injection_state
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_equal_i64, assert_close_dp, assert_allclose_1d, &
                          delete_file_if_exists, ensure_directory, remove_empty_directory
  implicit none

  type(mesh_type) :: mesh
  type(sim_stats) :: stats
  type(injection_state) :: state
  type(app_config) :: cfg, cfg_changed
  type(charge_ledger_type) :: ledger, restored_ledger
  type(electrostatic_diagnostics_type) :: electrostatic_diagnostics
  type(electrostatic_restart_state_type) :: electrostatic_state
  logical :: has_restart, exists
  integer(i32) :: contract_status
  character(len=256) :: contract_message
  character(len=512) :: profile_header
  integer :: profile_unit, profile_ios
  character(len=1024) :: rng_rank_path, residual_global_path
  character(len=*), parameter :: out_dir = 'test_restart_tmp'

  call delete_file_if_exists(out_dir//'/summary.txt')
  call delete_file_if_exists(out_dir//'/charges.csv')
  call delete_file_if_exists(out_dir//'/rng_state.txt')
  call delete_file_if_exists(out_dir//'/macro_residuals.csv')
  call delete_file_if_exists(out_dir//'/rng_state_rank00001.txt')
  call delete_file_if_exists(out_dir//'/macro_residuals_rank00001.csv')
  call delete_file_if_exists(out_dir//'/outer_plasma_profile.csv')
  call remove_empty_directory(out_dir)

  call test_init(7)

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

  call test_summary()

contains

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

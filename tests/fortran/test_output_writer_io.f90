!> output writer の公開 I/O contract と resume 時の履歴追記を検証する。
program test_output_writer_io
  use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: qe
  use bem_mesh, only: init_mesh
  use bem_output_writer, only: ensure_output_dir, open_top_reference_history_writer, open_matching_plane_history_writer, &
                               write_result_files, write_top_reference_history_snapshot, &
                               open_fixed_current_history_writer, write_fixed_current_history_rows, &
                               open_matching_plane_spectrum_history_writer, write_matching_plane_spectrum_history_rows, &
                               write_matching_plane_history_snapshot, &
                               open_secondary_emission_history_writer, write_secondary_emission_history_rows
  use bem_outflow_spectrum, only: outflow_spectrum_bins, outflow_spectrum_edges_ev
  use bem_app_config, only: app_config, default_app_config, load_app_config
  use bem_external_boundary_contract, only: external_inflow_scalar_barrier
  use bem_types, only: mesh_type, sim_stats
  use bem_charge_ledger, only: charge_ledger_type
  use bem_electrostatic_snapshot, only: electrostatic_diagnostics_type
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, delete_file_if_exists, remove_empty_directory
  implicit none

  type(mesh_type) :: mesh
  type(app_config) :: cfg
  type(sim_stats) :: stats
  type(charge_ledger_type) :: ledger
  type(electrostatic_diagnostics_type) :: electrostatic_diagnostics
  logical :: exists, literal_created, marker_created, saw_integrator, saw_residual, saw_ledger_header
  logical :: saw_schema, saw_mesh_fp, saw_ledger_stock, saw_ledger_closure
  logical :: saw_build_schema, saw_build_version, saw_build_mode, saw_source_commit, saw_build_id
  logical :: saw_surface_current_model, saw_soft_discard_fraction
  logical :: saw_photoelectron_active_receipt
  logical :: saw_refresh_receipts(3), matching_history_opened
  logical :: saw_field_reconstruction(23), saw_auto_resolved_direct, saw_auto_resolved_fmm
  logical :: top_history_opened, saw_top_available, saw_top_definition, saw_top_last_batch, saw_top_mean
  integer :: literal_unit, ios, top_history_unit, matching_history_unit
  integer(i32) :: top_batch, top_sample_n
  real(dp) :: top_time, top_z, top_mean, top_std, top_min, top_max, switch_phi_rms, switch_ez_rms
  character(len=2048) :: line
  character(len=*), parameter :: out_dir_disabled = 'test_output_writer_io_disabled_tmp'
  character(len=*), parameter :: out_dir_ledger = 'test_output_writer_io_ledger_tmp'
  character(len=*), parameter :: out_dir_no_photo = 'test_output_writer_io_no_photo_tmp'
  character(len=*), parameter :: out_dir_matching = 'test_output_writer_io_matching_tmp'
  character(len=*), parameter :: literal_parent = 'test_output_writer_io_literal_tmp'
  character(len=*), parameter :: marker_path = 'test_output_writer_io_shell_marker_tmp'
  character(len=*), parameter :: literal_dir = &
                                 literal_parent//'/space $(touch '//marker_path//'); "double" ''single'''
  character(len=*), parameter :: expanded_dir = literal_parent//'/space ; double ''single'''
  character(len=*), parameter :: charge_ledger_csv_header = &
                                 'batch,species_idx,injected_from_remote_C,emitted_from_surface_C,'// &
                                 'absorbed_on_surface_C,escaped_to_infinity_C,discarded_unresolved_C,'// &
                                 'neutral_return_correction_C,neutral_return_weight_scale,'// &
                                 'neutral_return_unresolved_fraction,fixed_absorbed_target_charge_C,'// &
                                 'fixed_absorbed_weight_scale,fixed_emission_target_charge_C,'// &
                                 'fixed_emission_weight_scale,fixed_current_correction_C,'// &
                                 'fixed_absorbed_applied_charge_C,fixed_emission_applied_charge_C,'// &
                                 'fixed_escape_target_charge_C,fixed_escape_applied_charge_C,'// &
                                 'fixed_escape_correction_C,injected_count,emitted_count,absorbed_count,'// &
                                 'escaped_count,discarded_unresolved_count'

  interface
    integer(c_int) function c_rmdir(path) bind(C, name='rmdir')
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: path(*)
    end function c_rmdir
  end interface

  call test_init(10)

  call build_two_element_mesh(mesh)
  mesh%q_elem = [2.0d-12, -1.0d-12]

  call test_begin('secondary_history_and_legacy_matching_history_resume')
  call test_secondary_and_legacy_history()
  call test_end()

  call cleanup_output_dir(out_dir_disabled)
  call cleanup_output_dir(out_dir_ledger)
  call cleanup_output_dir(out_dir_no_photo)
  call cleanup_output_dir(out_dir_matching)

  call delete_file_if_exists(marker_path)
  call remove_test_directory(literal_dir)
  call remove_test_directory(expanded_dir)
  call remove_test_directory(literal_parent)

  call test_begin('output_directory_path_is_literal')
  call ensure_output_dir(literal_dir)
  open (newunit=literal_unit, file=literal_dir//'/probe', status='replace', action='write', iostat=ios)
  literal_created = (ios == 0)
  if (literal_created) close (literal_unit, status='delete')
  inquire (file=marker_path, exist=marker_created)

  if (marker_created) call delete_file_if_exists(marker_path)
  call remove_test_directory(literal_dir)
  call remove_test_directory(expanded_dir)
  call remove_test_directory(literal_parent)

  call assert_true(.not. marker_created, 'output directory path must not execute shell command substitution')
  call assert_true(literal_created, 'output directory path with shell metacharacters should be created literally')
  call test_end()

  call test_begin('resolved_auto_solver_and_disabled_mesh_potential')
  call default_app_config(cfg)
  stats = sim_stats()
  cfg%sim%batch_duration = 1.0_dp
  cfg%output_dir = out_dir_disabled
  cfg%write_mesh_potential = .false.
  call write_result_files(out_dir_disabled, mesh, stats, cfg)
  call scan_resolved_field_solver( &
    out_dir_disabled//'/summary.txt', 'direct', saw_auto_resolved_direct &
    )
  call assert_true(saw_auto_resolved_direct, 'auto field solver receipt should resolve a small mesh to direct')
  cfg%sim%tree_min_nelem = 2_i32
  call write_result_files(out_dir_disabled, mesh, stats, cfg)
  call scan_resolved_field_solver( &
    out_dir_disabled//'/summary.txt', 'fmm', saw_auto_resolved_fmm &
    )
  call assert_true(saw_auto_resolved_fmm, 'auto field solver receipt should record the resolved fmm backend')
  inquire (file=trim(out_dir_disabled)//'/mesh_potential.csv', exist=exists)
  call assert_true(.not. exists, 'mesh_potential.csv should not be written when output.write_mesh_potential=false')
  call test_end()

  call test_begin('resolved_local_boundary_receipt')
  call default_app_config(cfg)
  stats = sim_stats()
  cfg%sim%reservoir_potential_model = 'infinity_barrier'
  cfg%boundary_contract%inflow_map = external_inflow_scalar_barrier
  call write_result_files(out_dir_disabled, mesh, stats, cfg)
  call assert_resolved_boundary_summary( &
    out_dir_disabled//'/summary.txt', 'infinity_barrier' &
    )

  call test_end()

  call test_begin('top_reference_history_and_summary')
  call default_app_config(cfg)
  stats = sim_stats()
  cfg%sim%batch_duration = 1.0_dp
  cfg%sim%use_box = .true.
  cfg%output_dir = out_dir_disabled
  cfg%write_potential_history = .true.
  cfg%history_stride = 1_i32
  call open_top_reference_history_writer(cfg, .false., top_history_opened, top_history_unit)
  call assert_true(top_history_opened, 'top-reference history should open with potential history')
  call write_top_reference_history_snapshot( &
    top_history_unit, 3_i32, 1.5_dp, 2.0_dp, 5_i32, 4.0_dp, 0.25_dp, 3.5_dp, 4.5_dp, 0.125_dp, 2.5_dp &
    )
  close (top_history_unit)
  call open_top_reference_history_writer(cfg, .true., top_history_opened, top_history_unit)
  call assert_true(top_history_opened, 'top-reference history should reopen for resume')
  call write_top_reference_history_snapshot( &
    top_history_unit, 4_i32, 2.0_dp, 2.0_dp, 5_i32, 4.25_dp, 0.20_dp, 3.8_dp, 4.7_dp &
    )
  close (top_history_unit)

  open (newunit=literal_unit, file=out_dir_disabled//'/top_reference_history.csv', &
        status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'failed to open top-reference history fixture'
  read (literal_unit, '(A)', iostat=ios) line
  call assert_true( &
    ios == 0 .and. trim(line) == &
    'batch,simulated_time_s,z_high_m,sample_n,potential_mean_V,potential_std_V,potential_min_V,potential_max_V,'// &
    'switch_phi_rms_mismatch_V,switch_Ez_rms_mismatch_V_m', &
    'top-reference history header mismatch' &
    )
  read (literal_unit, *, iostat=ios) &
    top_batch, top_time, top_z, top_sample_n, top_mean, top_std, top_min, top_max, switch_phi_rms, switch_ez_rms
  call assert_true(ios == 0, 'failed to parse top-reference history row')
  call assert_equal_i32(top_batch, 3_i32, 'top-reference history batch mismatch')
  call assert_close_dp(switch_phi_rms, 0.125_dp, 0.0_dp, 'switch potential rms column')
  call assert_close_dp(switch_ez_rms, 2.5_dp, 0.0_dp, 'switch field rms column')
  call assert_equal_i32(top_sample_n, 5_i32, 'top-reference history sample count mismatch')
  call assert_true( &
    abs(top_time - 1.5_dp) < 1.0e-14_dp .and. abs(top_z - 2.0_dp) < 1.0e-14_dp .and. &
    abs(top_mean - 4.0_dp) < 1.0e-14_dp .and. abs(top_std - 0.25_dp) < 1.0e-14_dp .and. &
    abs(top_min - 3.5_dp) < 1.0e-14_dp .and. abs(top_max - 4.5_dp) < 1.0e-14_dp, &
    'top-reference history values mismatch' &
    )
  read (literal_unit, *, iostat=ios) &
    top_batch, top_time, top_z, top_sample_n, top_mean, top_std, top_min, top_max, switch_phi_rms, switch_ez_rms
  call assert_true(ios == 0, 'failed to parse resumed top-reference history row')
  call assert_equal_i32(top_batch, 4_i32, 'resumed top-reference history batch mismatch')
  call assert_close_dp(switch_phi_rms, 0.0_dp, 0.0_dp, 'off switch potential rms column')
  call assert_close_dp(switch_ez_rms, 0.0_dp, 0.0_dp, 'off switch field rms column')
  call assert_equal_i32(top_sample_n, 5_i32, 'resumed top-reference history sample count mismatch')
  call assert_true( &
    abs(top_time - 2.0_dp) < 1.0e-14_dp .and. abs(top_z - 2.0_dp) < 1.0e-14_dp .and. &
    abs(top_mean - 4.25_dp) < 1.0e-14_dp .and. abs(top_std - 0.20_dp) < 1.0e-14_dp .and. &
    abs(top_min - 3.8_dp) < 1.0e-14_dp .and. abs(top_max - 4.7_dp) < 1.0e-14_dp, &
    'resumed top-reference history values mismatch' &
    )
  read (literal_unit, '(A)', iostat=ios) line
  call assert_true(ios < 0, 'resumed top-reference history must not duplicate the header')
  close (literal_unit)

  electrostatic_diagnostics = electrostatic_diagnostics_type()
  electrostatic_diagnostics%top_reference_available = .true.
  electrostatic_diagnostics%top_reference_last_batch = 4_i32
  electrostatic_diagnostics%top_reference_simulated_time = 2.0_dp
  electrostatic_diagnostics%top_reference_z_high = 2.0_dp
  electrostatic_diagnostics%top_reference_sample_n = 5_i32
  electrostatic_diagnostics%top_reference_potential_mean = 4.25_dp
  electrostatic_diagnostics%top_reference_potential_std = 0.20_dp
  electrostatic_diagnostics%top_reference_potential_min = 3.8_dp
  electrostatic_diagnostics%top_reference_potential_max = 4.7_dp
  call write_result_files( &
    out_dir_disabled, mesh, stats, cfg, electrostatic_diagnostics=electrostatic_diagnostics &
    )
  call scan_top_reference_summary_fields( &
    out_dir_disabled//'/summary.txt', saw_top_available, saw_top_definition, saw_top_last_batch, saw_top_mean &
    )
  call assert_true( &
    saw_top_available .and. saw_top_definition .and. saw_top_last_batch .and. saw_top_mean, &
    'available top-reference diagnostics must be written to summary' &
    )
  call test_end()

  call test_begin('charge_ledger_and_model_metadata')
  call default_app_config(cfg)
  cfg%sim%use_box = .true.
  stats = sim_stats()
  call ledger%init(2_i32)
  call ledger%reset(1_i32)
  ledger%surface_charge_after = -3.0_dp
  ledger%local_flight_charge_before = -1.0_dp
  ledger%local_flight_charge_after = -2.0_dp
  ledger%injected_from_remote(1) = -3.0_dp
  ledger%absorbed_on_surface(1) = -3.0_dp
  ledger%neutral_return_correction(2) = -0.25_dp
  ledger%neutral_return_weight_scale(2) = 1.25_dp
  ledger%neutral_return_unresolved_fraction(2) = 0.2_dp
  ledger%injected_count(1) = 1
  ledger%absorbed_count(1) = 1
  cfg%output_dir = out_dir_ledger
  cfg%sim%field_solver = 'fmm'
  stats%processed_particles = 8
  stats%multiple_box_events_soft_discarded = 2
  call write_result_files(out_dir_ledger, mesh, stats, cfg, charge_ledger=ledger)

  saw_integrator = .false.
  saw_residual = .false.
  saw_schema = .false.
  saw_mesh_fp = .false.
  saw_ledger_stock = .false.
  saw_ledger_closure = .false.
  saw_build_schema = .false.
  saw_build_version = .false.
  saw_build_mode = .false.
  saw_source_commit = .false.
  saw_build_id = .false.
  saw_surface_current_model = .false.
  saw_soft_discard_fraction = .false.
  saw_field_reconstruction = .false.
  open (newunit=literal_unit, file=out_dir_ledger//'/summary.txt', status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'failed to open summary metadata fixture'
  do
    read (literal_unit, '(A)', iostat=ios) line
    if (ios /= 0) exit
    saw_integrator = saw_integrator .or. index(line, 'particle_time_centering=same_time_midpoint_boris') > 0
    saw_residual = saw_residual .or. &
                   summary_real_equals(line, 'charge_ledger_residual_C', -0.75_dp, 1.0e-14_dp)
    saw_schema = saw_schema .or. index(line, 'checkpoint_schema_version=11') > 0
    saw_mesh_fp = saw_mesh_fp .or. index(line, 'mesh_fingerprint=') > 0
    saw_ledger_stock = saw_ledger_stock .or. &
                       summary_real_equals( &
                       line, 'charge_ledger_local_flight_charge_before_C', -1.0_dp, 1.0e-14_dp &
                       )
    saw_ledger_closure = saw_ledger_closure .or. &
                         summary_real_equals( &
                         line, 'charge_ledger_neutral_return_correction_C', -0.25_dp, 1.0e-14_dp &
                         )
    saw_build_schema = saw_build_schema .or. index(line, 'build_info_schema_version=1') == 1
    saw_build_version = saw_build_version .or. index(line, 'build_version=') == 1
    saw_build_mode = saw_build_mode .or. index(line, 'build_version_mode=') == 1
    saw_source_commit = saw_source_commit .or. index(line, 'build_source_commit=') == 1
    saw_build_id = saw_build_id .or. index(line, 'build_id=') == 1
    saw_surface_current_model = saw_surface_current_model .or. index(line, 'surface_current_model=none') == 1
    saw_soft_discard_fraction = saw_soft_discard_fraction .or. &
                                summary_real_equals( &
                                line, 'multiple_box_events_soft_discard_fraction', 0.25_dp, 1.0e-14_dp &
                                )
    saw_field_reconstruction(1) = saw_field_reconstruction(1) .or. &
                                  trim(line) == 'field_reconstruction_schema_version=2'
    saw_field_reconstruction(2) = saw_field_reconstruction(2) .or. &
                                  trim(line) == 'field_reconstruction_resolved_field_solver=fmm'
    saw_field_reconstruction(3) = saw_field_reconstruction(3) .or. &
                                  trim(line) == 'field_reconstruction_fmm_expansion_order=4'
    saw_field_reconstruction(4) = saw_field_reconstruction(4) .or. &
                                  trim(line) == 'field_reconstruction_field_bc_mode=free'
    saw_field_reconstruction(5) = saw_field_reconstruction(5) .or. &
                                  index(line, 'field_reconstruction_tree_theta=') == 1
    saw_field_reconstruction(6) = saw_field_reconstruction(6) .or. &
                                  trim(line) == 'field_reconstruction_tree_leaf_max=12'
    saw_field_reconstruction(7) = saw_field_reconstruction(7) .or. &
                                  index(line, 'field_reconstruction_e0_V_m=') == 1
    saw_field_reconstruction(8) = saw_field_reconstruction(8) .or. &
                                  trim(line) == 'field_reconstruction_use_box=T'
    saw_field_reconstruction(9) = saw_field_reconstruction(9) .or. &
                                  index(line, 'field_reconstruction_box_min_m=') == 1
    saw_field_reconstruction(10) = saw_field_reconstruction(10) .or. &
                                   index(line, 'field_reconstruction_box_max_m=') == 1
    saw_field_reconstruction(11) = saw_field_reconstruction(11) .or. &
                                   index(line, 'field_reconstruction_boundary_low=') == 1
    saw_field_reconstruction(12) = saw_field_reconstruction(12) .or. &
                                   index(line, 'field_reconstruction_boundary_high=') == 1
    saw_field_reconstruction(13) = saw_field_reconstruction(13) .or. &
                                   trim(line) == 'field_reconstruction_periodic_image_layers=1'
    saw_field_reconstruction(14) = saw_field_reconstruction(14) .or. &
                                   trim(line) == 'field_reconstruction_periodic_far_correction=none'
    saw_field_reconstruction(15) = saw_field_reconstruction(15) .or. &
                                   trim(line) == 'field_reconstruction_periodic_nonzero_mode_backend=not_applicable'
    saw_field_reconstruction(16) = saw_field_reconstruction(16) .or. &
                                   trim(line) == 'field_reconstruction_periodic_zero_mode_policy=not_applicable'
    saw_field_reconstruction(17) = saw_field_reconstruction(17) .or. &
                                   trim(line) == 'field_reconstruction_periodic_lower_boundary_model=not_applicable'
    saw_field_reconstruction(18) = saw_field_reconstruction(18) .or. &
                                   trim(line) == 'field_reconstruction_periodic_reference_mode_layers=4'
    saw_field_reconstruction(19) = saw_field_reconstruction(19) .or. &
                                   trim(line) == 'field_reconstruction_periodic_panel_quadrature_order=12'
    saw_field_reconstruction(20) = saw_field_reconstruction(20) .or. &
                                   index(line, 'field_reconstruction_periodic_ewald_alpha=') == 1
    saw_field_reconstruction(21) = saw_field_reconstruction(21) .or. &
                                   trim(line) == 'field_reconstruction_periodic_ewald_layers=4'
    saw_field_reconstruction(22) = saw_field_reconstruction(22) .or. &
                                   trim(line) == &
                                   'field_reconstruction_periodic_cache_dir='//trim(cfg%sim%field_periodic_cache_dir)
    saw_field_reconstruction(23) = saw_field_reconstruction(23) .or. &
                                   index(line, 'field_reconstruction_periodic_generation_tolerance=') == 1
  end do
  close (literal_unit)
  open (newunit=literal_unit, file=out_dir_ledger//'/charge_ledger.csv', status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'failed to open charge ledger fixture'
  read (literal_unit, '(A)', iostat=ios) line
  close (literal_unit)
  saw_ledger_header = ios == 0 .and. trim(line) == charge_ledger_csv_header
  call assert_true(saw_integrator, 'summary should record the particle time-centering contract')
  call assert_true(saw_residual, 'summary should record the charge ledger residual')
  call assert_true(saw_schema, 'summary should record checkpoint schema v11')
  call assert_true(saw_mesh_fp, 'summary should identify the ordered mesh')
  call assert_true(saw_build_schema .and. saw_build_version .and. saw_build_mode .and. saw_source_commit .and. saw_build_id, &
                   'summary should record executable build origin')
  call assert_true(saw_ledger_stock, 'summary should record restartable charge stocks')
  call assert_true(saw_ledger_closure, 'summary should record the neutral-return correction')
  call assert_true(saw_surface_current_model, 'summary should record the surface-current model receipt')
  call assert_true(saw_soft_discard_fraction, 'summary should record the cumulative soft-discard fraction')
  call assert_true(all(saw_field_reconstruction), 'summary should record the field-reconstruction receipt')
  call assert_true(saw_ledger_header, 'charge ledger CSV header mismatch')
  call test_end()

  call test_begin('no_photo_surface_current_receipt')
  call default_app_config(cfg)
  stats = sim_stats()
  call load_app_config('examples/zero_current_no_photo.toml', cfg)
  cfg%output_dir = out_dir_no_photo
  call write_result_files(out_dir_no_photo, mesh, stats, cfg)
  call scan_no_photo_surface_current_receipt( &
    out_dir_no_photo//'/summary.txt', saw_photoelectron_active_receipt &
    )
  call assert_true( &
    saw_photoelectron_active_receipt, &
    'no-PE summary should record surface_current_model_photoelectron_active=F' &
    )
  call test_end()

  call test_begin('outflow_refresh_history_and_receipts')
  call default_app_config(cfg)
  stats = sim_stats()
  call load_app_config('examples/zero_current_refresh.toml', cfg)
  cfg%output_dir = out_dir_matching
  cfg%history_stride = 0_i32
  cfg%write_mesh_potential = .false.
  call ensure_output_dir(out_dir_matching)
  open (newunit=literal_unit, file=out_dir_matching//'/matching_plane_history.csv', &
        status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'failed to create stale outer-state history fixture'
  write (literal_unit, '(a)') 'stale'
  close (literal_unit)
  call open_matching_plane_history_writer(cfg, .false., matching_history_opened, matching_history_unit)
  call assert_true(.not. matching_history_opened, 'history_stride=0 must disable the outer-state history')
  inquire (file=out_dir_matching//'/matching_plane_history.csv', exist=exists)
  call assert_true(.not. exists, 'fresh run must remove a stale disabled outer-state history')
  cfg%history_stride = 1_i32
  call open_matching_plane_history_writer(cfg, .false., matching_history_opened, matching_history_unit)
  call assert_true(matching_history_opened, 'outflow refresh must open the outer-state history')
  if (matching_history_opened) close (matching_history_unit)

  call write_result_files(out_dir_matching, mesh, stats, cfg)
  call scan_summary_line( &
    out_dir_matching//'/summary.txt', 'surface_current_model_outflow_refresh_batches=2', saw_refresh_receipts(1) &
    )
  call scan_summary_line( &
    out_dir_matching//'/summary.txt', &
    'surface_current_model_dynamic_state_source=matching_plane_state_outflow_refresh', saw_refresh_receipts(2) &
    )
  call scan_summary_line( &
    out_dir_matching//'/summary.txt', 'surface_current_model_outer_return_position=cell_uniform', &
    saw_refresh_receipts(3) &
    )
  call assert_true(all(saw_refresh_receipts), 'summary should record the outflow refresh contract')
  call test_end()

  call test_begin('fixed_current_history_rows')
  call test_fixed_current_history_rows()
  call test_end()

  call test_begin('matching_plane_spectrum_history_and_summary')
  call test_matching_plane_spectrum_outputs()
  call test_end()

  call cleanup_output_dir(out_dir_disabled)
  call cleanup_output_dir(out_dir_ledger)
  call cleanup_output_dir(out_dir_no_photo)
  call cleanup_output_dir(out_dir_matching)

  call test_summary()

contains

  !> 固定電流historyはtargetを適用したchannelだけを窓合計で書き、追跡値0の比をnanにする。
  subroutine test_fixed_current_history_rows()
    logical :: opened, active(3, 3)
    integer :: unit_id, read_ios, row_count
    real(dp) :: tracked(3, 3), target(3, 3)
    character(len=512) :: rows(5)

    call default_app_config(cfg)
    call load_app_config('examples/zero_current_refresh.toml', cfg)
    cfg%output_dir = out_dir_matching
    cfg%history_stride = 0_i32
    call ensure_output_dir(out_dir_matching)
    open (newunit=unit_id, file=out_dir_matching//'/fixed_current_history.csv', status='replace', action='write')
    write (unit_id, '(a)') 'stale'
    close (unit_id)
    call open_fixed_current_history_writer(cfg, .false., opened, unit_id)
    call assert_true(.not. opened, 'history_stride=0 must disable the fixed-current history')
    inquire (file=out_dir_matching//'/fixed_current_history.csv', exist=exists)
    call assert_true(.not. exists, 'fresh run must remove a stale disabled fixed-current history')

    cfg%history_stride = 2_i32
    call open_fixed_current_history_writer(cfg, .false., opened, unit_id)
    call assert_true(opened, 'zhao_stationary must open the fixed-current history')
    if (.not. opened) return
    active = .false.
    active(1, 1:2) = .true.
    active(2:3, 3) = .true.
    tracked = 0.0_dp
    target = 0.0_dp
    tracked(1, 1) = -2.0_dp
    target(1, 1) = -3.0_dp
    tracked(1, 2) = 1.0e-12_dp
    target(1, 2) = 1.0e-12_dp
    tracked(2, 3) = 4.0e-12_dp
    target(2, 3) = 4.0e-12_dp
    target(3, 3) = -1.0e-13_dp
    call write_fixed_current_history_rows(unit_id, 3_i32, 6.0_dp, 2_i32, active, tracked, target)
    close (unit_id)

    open (newunit=unit_id, file=out_dir_matching//'/fixed_current_history.csv', status='old', action='read')
    row_count = 0
    do
      read (unit_id, '(a)', iostat=read_ios) line
      if (read_ios /= 0) exit
      row_count = row_count + 1
      if (row_count <= size(rows)) rows(row_count) = line
    end do
    close (unit_id)
    call assert_equal_i32(int(row_count, i32), 5_i32, 'header plus one row per active channel expected')
    if (row_count /= 5) return
    call assert_true( &
      trim(rows(1)) == 'batch,simulated_time_s,window_batches,species_idx,channel,tracked_charge_C,target_charge_C,'// &
      'target_over_tracked', 'fixed-current history header mismatch' &
      )
    call assert_true(index(rows(2), ',2,1,absorbed,') > 0, 'electron absorbed row mismatch')
    call assert_true(index(rows(2), ',1.5000000000000000E+00') > 0, 'electron target/tracked ratio mismatch')
    call assert_true(index(rows(3), ',2,2,absorbed,') > 0, 'ion absorbed row mismatch')
    call assert_true(index(rows(4), ',2,3,emission,') > 0, 'PE emission row mismatch')
    call assert_true( &
      index(rows(5), ',2,3,escape,') > 0 .and. rows(5) (len_trim(rows(5)) - 3:len_trim(rows(5))) == ',nan', &
      'PE escape row with no tracked charge must write nan' &
      )
  end subroutine test_fixed_current_history_rows

  !> 外部源スペクトルは更新ごとに正の bin だけを history へ書き、summary へは全 bin を 1 行で書く。
  subroutine test_secondary_and_legacy_history()
    integer(i64) :: counts(3, 2)
    real(dp) :: values(2, 2)
    logical :: opened
    integer :: unit_id, reader, comma, j
    character(len=2048) :: header, legacy_header, row

    call default_app_config(cfg)
    cfg%output_dir = out_dir_matching
    cfg%n_particle_species = 2_i32
    cfg%particle_species(1)%enabled = .true.
    cfg%particle_species(1)%species_key = 'incident,electron'
    cfg%particle_species(1)%secondary_emission_enabled = .true.
    cfg%history_stride = 0_i32
    call open_secondary_emission_history_writer(cfg, .false., opened, unit_id)
    call assert_true(opened, 'SEE history is independent of ordinary history stride')
    counts = 0_i64
    values = 0.0_dp
    counts(:, 1) = [1_i64, 2_i64, 3_i64]
    values(:, 1) = [-2.0_dp*qe, 1.5_dp]
    call write_secondary_emission_history_rows(unit_id, 1_i32, cfg, counts, values)
    close (unit_id)
    open (newunit=reader, file=out_dir_matching//'/secondary_emission_history.csv', status='old', action='read')
    read (reader, '(a)') header
    read (reader, '(a)') row
    close (reader)
    call assert_true(index(header, 'generation_capped_count') > 0, 'SEE history columns')
    call assert_true(index(row, '"incident,electron"') > 0, 'SEE history quotes incident species')

    cfg%history_stride = 1_i32
    cfg%surface_current%model = 'zhao_stationary'
    cfg%surface_current%outflow_refresh_batches = 1_i32
    call open_matching_plane_history_writer(cfg, .false., opened, unit_id)
    close (unit_id)
    open (newunit=reader, file=out_dir_matching//'/matching_plane_history.csv', status='old', action='read')
    read (reader, '(a)') header
    close (reader)
    comma = index(header, ',secondary_share')
    legacy_header = header(:comma - 1)
    open (newunit=unit_id, file=out_dir_matching//'/matching_plane_history.csv', status='replace', action='write')
    write (unit_id, '(a)') trim(legacy_header)
    write (unit_id, '(a)') '1,1,0,0,1,1,0,0,0,2,2,0,0,1,1,1,0'
    close (unit_id)
    call open_matching_plane_history_writer(cfg, .true., opened, unit_id)
    stats = sim_stats()
    stats%matching_plane%state_valid = .true.
    stats%matching_plane%secondary_share = 0.25_dp
    stats%matching_plane%secondary_budget_residual_a_m2 = -1.0e-8_dp
    call write_matching_plane_history_snapshot(unit_id, 2_i32, 2.0_dp, stats)
    close (unit_id)
    open (newunit=reader, file=out_dir_matching//'/matching_plane_history.csv', status='old', action='read')
    read (reader, '(a)') row
    call assert_true(trim(row) == trim(header), 'legacy header must upgrade to the 21-column contract')
    do j = 1, 2
      read (reader, '(a)') row
      comma = 0
      block
        integer :: k
        do k = 1, len_trim(row)
          if (row(k:k) == ',') comma = comma + 1
        end do
      end block
      call assert_equal_i32(int(comma, i32), 20_i32, 'legacy and new matching rows must have 21 columns')
    end do
    close (reader)
  end subroutine test_secondary_and_legacy_history

  subroutine test_matching_plane_spectrum_outputs()
    logical :: opened, found
    integer :: unit_id, read_ios, row_count, pos
    integer(i32) :: row_batch, row_bin
    real(dp) :: row_time, row_low, row_high, row_flux, edges(outflow_spectrum_bins + 1)
    real(dp) :: summary_spectrum(outflow_spectrum_bins)
    character(len=512) :: rows(3)
    character(len=8192) :: long_line

    call default_app_config(cfg)
    call load_app_config('examples/zero_current_refresh.toml', cfg)
    cfg%output_dir = out_dir_matching
    cfg%history_stride = 0_i32
    cfg%write_mesh_potential = .false.
    call ensure_output_dir(out_dir_matching)
    open (newunit=unit_id, file=out_dir_matching//'/matching_plane_spectrum_history.csv', status='replace', &
          action='write')
    write (unit_id, '(a)') 'stale'
    close (unit_id)
    call open_matching_plane_spectrum_history_writer(cfg, .false., opened, unit_id)
    call assert_true(.not. opened, 'history_stride=0 must disable the spectrum history')
    inquire (file=out_dir_matching//'/matching_plane_spectrum_history.csv', exist=exists)
    call assert_true(.not. exists, 'fresh run must remove a stale disabled spectrum history')

    cfg%history_stride = 5_i32
    call open_matching_plane_spectrum_history_writer(cfg, .false., opened, unit_id)
    call assert_true(opened, 'outflow refresh must open the spectrum history')
    if (.not. opened) return
    stats = sim_stats()
    stats%matching_plane%photoelectron_spectrum_m2_s(3) = 1.0e10_dp
    call write_matching_plane_spectrum_history_rows(unit_id, 1_i32, 2.0_dp, stats)
    stats%matching_plane%state_valid = .true.
    stats%matching_plane%photoelectron_spectrum_m2_s(100) = 2.0e10_dp
    call write_matching_plane_spectrum_history_rows(unit_id, 3_i32, 6.0_dp, stats)
    close (unit_id)

    open (newunit=unit_id, file=out_dir_matching//'/matching_plane_spectrum_history.csv', status='old', &
          action='read')
    row_count = 0
    do
      read (unit_id, '(a)', iostat=read_ios) line
      if (read_ios /= 0) exit
      row_count = row_count + 1
      if (row_count <= size(rows)) rows(row_count) = line
    end do
    close (unit_id)
    call assert_equal_i32( &
      int(row_count, i32), 3_i32, 'header plus one row per positive bin of the valid state expected' &
      )
    if (row_count /= 3) return
    call assert_true( &
      trim(rows(1)) == 'batch,simulated_time_s,bin,energy_low_eV,energy_high_eV,photoelectron_outward_flux_m2_s', &
      'spectrum history header mismatch' &
      )
    edges = outflow_spectrum_edges_ev()
    read (rows(3), *) row_batch, row_time, row_bin, row_low, row_high, row_flux
    call assert_equal_i32(row_batch, 3_i32, 'spectrum history batch mismatch')
    call assert_equal_i32(row_bin, 100_i32, 'spectrum history bin mismatch')
    call assert_true( &
      row_low == edges(100) .and. row_high == edges(101) .and. row_flux == 2.0e10_dp, &
      'spectrum history row must keep the bin edges and flux' &
      )

    call write_result_files(out_dir_matching, mesh, stats, cfg)
    found = .false.
    open (newunit=unit_id, file=out_dir_matching//'/summary.txt', status='old', action='read')
    do
      read (unit_id, '(a)', iostat=read_ios) long_line
      if (read_ios /= 0) exit
      pos = index(long_line, '=')
      if (pos <= 0) cycle
      if (long_line(:pos - 1) /= 'matching_plane_photoelectron_spectrum_m2_s') cycle
      read (long_line(pos + 1:), *, iostat=read_ios) summary_spectrum
      found = read_ios == 0
      exit
    end do
    close (unit_id)
    call assert_true(found, 'summary must write every spectrum bin on one line')
    call assert_true( &
      all(summary_spectrum == stats%matching_plane%photoelectron_spectrum_m2_s), &
      'summary spectrum must round-trip exactly' &
      )
  end subroutine test_matching_plane_spectrum_outputs

  subroutine scan_summary_line(summary_path, expected_line, found)
    character(len=*), intent(in) :: summary_path, expected_line
    logical, intent(out) :: found
    integer :: summary_unit, summary_ios
    character(len=2048) :: summary_line

    found = .false.
    open (newunit=summary_unit, file=trim(summary_path), status='old', action='read', iostat=summary_ios)
    if (summary_ios /= 0) error stop 'failed to open summary fixture'
    do
      read (summary_unit, '(A)', iostat=summary_ios) summary_line
      if (summary_ios /= 0) exit
      if (trim(summary_line) == trim(expected_line)) then
        found = .true.
        exit
      end if
    end do
    close (summary_unit)
  end subroutine scan_summary_line

  subroutine scan_resolved_field_solver(summary_path, expected_solver, found)
    character(len=*), intent(in) :: summary_path, expected_solver
    logical, intent(out) :: found
    integer :: summary_unit, summary_ios
    character(len=2048) :: summary_line

    found = .false.
    open (newunit=summary_unit, file=trim(summary_path), status='old', action='read', iostat=summary_ios)
    if (summary_ios /= 0) error stop 'failed to open field-solver receipt fixture'
    do
      read (summary_unit, '(A)', iostat=summary_ios) summary_line
      if (summary_ios /= 0) exit
      found = found .or. trim(summary_line) == &
        'field_reconstruction_resolved_field_solver='//trim(expected_solver)
    end do
    close (summary_unit)
  end subroutine scan_resolved_field_solver

  subroutine scan_no_photo_surface_current_receipt(summary_path, found)
    character(len=*), intent(in) :: summary_path
    logical, intent(out) :: found
    integer :: summary_unit, summary_ios
    character(len=2048) :: summary_line

    found = .false.
    open (newunit=summary_unit, file=trim(summary_path), status='old', action='read', iostat=summary_ios)
    if (summary_ios /= 0) error stop 'failed to open no-PE surface-current receipt fixture'
    do
      read (summary_unit, '(A)', iostat=summary_ios) summary_line
      if (summary_ios /= 0) exit
      found = found .or. trim(summary_line) == 'surface_current_model_photoelectron_active=F'
    end do
    close (summary_unit)
  end subroutine scan_no_photo_surface_current_receipt

  subroutine assert_resolved_boundary_summary(path, inflow_map)
    character(len=*), intent(in) :: path, inflow_map
    logical :: saw_inflow
    integer :: unit, read_status
    character(len=512) :: summary_line

    saw_inflow = .false.
    open (newunit=unit, file=path, status='old', action='read', iostat=read_status)
    if (read_status /= 0) error stop 'failed to open resolved local boundary summary fixture'
    do
      read (unit, '(A)', iostat=read_status) summary_line
      if (read_status /= 0) exit
      saw_inflow = saw_inflow .or. trim(summary_line) == 'reservoir_inflow_map='//trim(inflow_map)
    end do
    close (unit)
    call assert_true(saw_inflow, 'summary should record the resolved reservoir inflow map')
  end subroutine assert_resolved_boundary_summary

  subroutine scan_top_reference_summary_fields(path, saw_available, saw_definition, saw_last_batch, saw_mean)
    character(len=*), intent(in) :: path
    logical, intent(out) :: saw_available, saw_definition, saw_last_batch, saw_mean
    integer :: unit, read_status
    character(len=512) :: summary_line

    saw_available = .false.
    saw_definition = .false.
    saw_last_batch = .false.
    saw_mean = .false.
    open (newunit=unit, file=path, status='old', action='read', iostat=read_status)
    if (read_status /= 0) error stop 'failed to open top-reference summary fixture'
    do
      read (unit, '(A)', iostat=read_status) summary_line
      if (read_status /= 0) exit
      saw_available = saw_available .or. trim(summary_line) == 'top_reference_available=T'
      saw_definition = saw_definition .or. &
                       trim(summary_line) == 'top_reference_definition=box_z_high_plane_mean'
      saw_last_batch = saw_last_batch .or. trim(summary_line) == 'top_reference_last_batch=4'
      saw_mean = saw_mean .or. &
                 summary_real_equals( &
                 summary_line, 'top_reference_potential_mean_V', 4.25_dp, 1.0e-14_dp &
                 )
    end do
    close (unit)
  end subroutine scan_top_reference_summary_fields

  !> summary の対象 key を数値として読み、公開値と比較する。
  logical function summary_real_equals(summary_line, key, expected, tolerance) result(matches)
    character(len=*), intent(in) :: summary_line, key
    real(dp), intent(in) :: expected, tolerance
    real(dp) :: actual
    integer :: parse_status, value_start

    matches = .false.
    if (index(summary_line, trim(key)//'=') /= 1) return
    value_start = len_trim(key) + 2
    if (value_start > len_trim(summary_line)) return
    read (summary_line(value_start:), *, iostat=parse_status) actual
    if (parse_status /= 0) return
    matches = abs(actual - expected) <= tolerance
  end function summary_real_equals

  !> 2 要素メッシュを初期化する。
  subroutine build_two_element_mesh(mesh)
    type(mesh_type), intent(out) :: mesh
    real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2)

    v0(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    v1(:, 1) = [1.0d0, 0.0d0, 0.0d0]
    v2(:, 1) = [0.0d0, 1.0d0, 0.0d0]
    v0(:, 2) = [0.0d0, 0.0d0, 1.0d0]
    v1(:, 2) = [1.0d0, 0.0d0, 1.0d0]
    v2(:, 2) = [0.0d0, 1.0d0, 1.0d0]
    call init_mesh(mesh, v0, v1, v2)
  end subroutine build_two_element_mesh

  !> writer テストの一時出力を削除する。
  subroutine cleanup_output_dir(out_dir)
    character(len=*), intent(in) :: out_dir

    call delete_file_if_exists(out_dir//'/summary.txt')
    call delete_file_if_exists(out_dir//'/charges.csv')
    call delete_file_if_exists(out_dir//'/mesh_potential.csv')
    call delete_file_if_exists(out_dir//'/mesh_triangles.csv')
    call delete_file_if_exists(out_dir//'/mesh_sources.csv')
    call delete_file_if_exists(out_dir//'/charge_ledger.csv')
    call delete_file_if_exists(out_dir//'/checkpoint_complete.txt')
    call delete_file_if_exists(out_dir//'/checkpoint_complete.txt.tmp')
    call delete_file_if_exists(out_dir//'/top_reference_history.csv')
    call delete_file_if_exists(out_dir//'/matching_plane_history.csv')
    call delete_file_if_exists(out_dir//'/secondary_emission_history.csv')
    call delete_file_if_exists(out_dir//'/fixed_current_history.csv')
    call delete_file_if_exists(out_dir//'/matching_plane_spectrum_history.csv')
    call remove_empty_directory(out_dir)
  end subroutine cleanup_output_dir

  !> テストで作成した空ディレクトリを shell を使わず削除する。
  subroutine remove_test_directory(path)
    character(len=*), intent(in) :: path
    character(kind=c_char), allocatable :: c_path(:)
    integer :: i, n
    integer(c_int) :: status

    n = len_trim(path)
    allocate (c_path(n + 1))
    do i = 1, n
      c_path(i) = path(i:i)
    end do
    c_path(n + 1) = c_null_char
    status = c_rmdir(c_path)
  end subroutine remove_test_directory

end program test_output_writer_io

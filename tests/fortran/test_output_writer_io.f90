!> output writer の公開 I/O contract と resume 時の履歴追記を検証する。
program test_output_writer_io
  use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
  use bem_kinds, only: dp, i32
  use bem_mesh, only: init_mesh
  use bem_output_writer, only: ensure_output_dir, open_top_reference_history_writer, open_matching_plane_history_writer, &
                               write_result_files, write_top_reference_history_snapshot
  use bem_app_config, only: app_config, default_app_config, load_app_config
  use bem_types, only: mesh_type, sim_stats
  use bem_charge_ledger, only: charge_ledger_type
  use bem_electrostatic_snapshot, only: electrostatic_diagnostics_type
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, delete_file_if_exists, remove_empty_directory
  implicit none

  type(mesh_type) :: mesh
  type(app_config) :: cfg
  type(sim_stats) :: stats
  type(charge_ledger_type) :: ledger
  type(electrostatic_diagnostics_type) :: electrostatic_diagnostics
  logical :: exists, literal_created, marker_created, saw_integrator, saw_residual, saw_ledger_header
  logical :: saw_schema, saw_model_fp, saw_mesh_fp, saw_species_fp, saw_ledger_stock, saw_ledger_closure
  logical :: saw_build_schema, saw_build_version, saw_build_mode, saw_source_commit, saw_build_id
  logical :: saw_surface_current_model, saw_soft_discard_fraction
  logical :: saw_photoelectron_active_receipt
  logical :: saw_matching_receipts(13), matching_history_opened, saw_continuation_state
  logical :: saw_online_matching_receipts(9)
  logical :: saw_field_reconstruction(23), saw_auto_resolved_direct, saw_auto_resolved_fmm
  logical :: top_history_opened, saw_top_available, saw_top_definition, saw_top_last_batch, saw_top_mean
  integer :: literal_unit, ios, top_history_unit, matching_history_unit
  integer(i32) :: top_batch, top_sample_n
  real(dp) :: top_time, top_z, top_mean, top_std, top_min, top_max
  character(len=2048) :: line
  character(len=*), parameter :: out_dir_disabled = 'test_output_writer_io_disabled_tmp'
  character(len=*), parameter :: out_dir_ledger = 'test_output_writer_io_ledger_tmp'
  character(len=*), parameter :: out_dir_no_photo = 'test_output_writer_io_no_photo_tmp'
  character(len=*), parameter :: out_dir_matching = 'test_output_writer_io_matching_tmp'
  character(len=*), parameter :: out_dir_matching_online = 'test_output_writer_io_matching_online_tmp'
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

  call test_init(8)

  call build_two_element_mesh(mesh)
  mesh%q_elem = [2.0d-12, -1.0d-12]

  call cleanup_output_dir(out_dir_disabled)
  call cleanup_output_dir(out_dir_ledger)
  call cleanup_output_dir(out_dir_no_photo)
  call cleanup_output_dir(out_dir_matching)
  call cleanup_output_dir(out_dir_matching_online)

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
  cfg%sim%open_boundary_model = 'potential_barrier'
  call write_result_files(out_dir_disabled, mesh, stats, cfg)
  call assert_resolved_boundary_summary( &
    out_dir_disabled//'/summary.txt', 'infinity_barrier', 'potential_barrier' &
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
    top_history_unit, 3_i32, 1.5_dp, 2.0_dp, 5_i32, 4.0_dp, 0.25_dp, 3.5_dp, 4.5_dp &
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
    'batch,simulated_time_s,z_high_m,sample_n,potential_mean_V,potential_std_V,potential_min_V,potential_max_V', &
    'top-reference history header mismatch' &
    )
  read (literal_unit, *, iostat=ios) &
    top_batch, top_time, top_z, top_sample_n, top_mean, top_std, top_min, top_max
  call assert_true(ios == 0, 'failed to parse top-reference history row')
  call assert_equal_i32(top_batch, 3_i32, 'top-reference history batch mismatch')
  call assert_equal_i32(top_sample_n, 5_i32, 'top-reference history sample count mismatch')
  call assert_true( &
    abs(top_time - 1.5_dp) < 1.0e-14_dp .and. abs(top_z - 2.0_dp) < 1.0e-14_dp .and. &
    abs(top_mean - 4.0_dp) < 1.0e-14_dp .and. abs(top_std - 0.25_dp) < 1.0e-14_dp .and. &
    abs(top_min - 3.5_dp) < 1.0e-14_dp .and. abs(top_max - 4.5_dp) < 1.0e-14_dp, &
    'top-reference history values mismatch' &
    )
  read (literal_unit, *, iostat=ios) &
    top_batch, top_time, top_z, top_sample_n, top_mean, top_std, top_min, top_max
  call assert_true(ios == 0, 'failed to parse resumed top-reference history row')
  call assert_equal_i32(top_batch, 4_i32, 'resumed top-reference history batch mismatch')
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
  saw_model_fp = .false.
  saw_mesh_fp = .false.
  saw_species_fp = .false.
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
    saw_schema = saw_schema .or. index(line, 'checkpoint_schema_version=9') > 0
    saw_model_fp = saw_model_fp .or. index(line, 'model_fingerprint=') > 0
    saw_mesh_fp = saw_mesh_fp .or. index(line, 'mesh_fingerprint=') > 0
    saw_species_fp = saw_species_fp .or. index(line, 'species_fingerprint=') > 0
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
                                   trim(line) == 'field_reconstruction_periodic_cache_dir=.beach_cache/periodic2'
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
  call assert_true(saw_schema, 'summary should record checkpoint schema v9')
  call assert_true(saw_model_fp .and. saw_mesh_fp .and. saw_species_fp, 'summary should record restart fingerprints')
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
  call load_app_config('examples/periodic2_zhao_no_photo_fixed_current.toml', cfg)
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

  call test_begin('matching_plane_provenance_and_stale_history')
  call default_app_config(cfg)
  stats = sim_stats()
  call load_app_config('examples/periodic2_matching_plane_quasistatic.toml', cfg)
  cfg%output_dir = out_dir_matching
  cfg%history_stride = 0_i32
  cfg%write_mesh_potential = .false.
  call ensure_output_dir(out_dir_matching)
  open (newunit=literal_unit, file=out_dir_matching//'/matching_plane_history.csv', &
        status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'failed to create stale matching-plane history fixture'
  write (literal_unit, '(a)') 'stale'
  close (literal_unit)
  call open_matching_plane_history_writer(cfg, .false., matching_history_opened, matching_history_unit)
  call assert_true(.not. matching_history_opened, 'history_stride=0 must disable matching-plane history')
  inquire (file=out_dir_matching//'/matching_plane_history.csv', exist=exists)
  call assert_true(.not. exists, 'fresh run must remove stale disabled matching-plane history')

  call write_result_files(out_dir_matching, mesh, stats, cfg)
  call scan_matching_plane_receipts(out_dir_matching//'/summary.txt', saw_matching_receipts)
  call assert_true(all(saw_matching_receipts), 'summary should record matching-plane provenance and controls')
  call test_end()

  call test_begin('online_matching_plane_solver_receipt')
  call default_app_config(cfg)
  stats = sim_stats()
  call load_app_config('examples/periodic2_matching_plane_zhao_implicit.toml', cfg)
  cfg%output_dir = out_dir_matching_online
  cfg%write_mesh_potential = .false.
  call write_result_files(out_dir_matching_online, mesh, stats, cfg)
  call scan_online_matching_plane_receipts( &
    out_dir_matching_online//'/summary.txt', saw_online_matching_receipts &
    )
  call assert_true(all(saw_online_matching_receipts), 'summary should record online Zhao solver provenance')
  cfg%surface_current%zhao_branch = 'a'
  cfg%surface_current%zhao_root_selection = 'continuation'
  cfg%surface_current%implicit_zero_mode = .true.
  call write_result_files(out_dir_matching_online, mesh, stats, cfg)
  call scan_summary_line( &
    out_dir_matching_online//'/summary.txt', &
    'surface_current_model_outer_solver_state=accepted_endpoint_continuation_v1', &
    saw_continuation_state &
    )
  call assert_true(saw_continuation_state, 'summary should identify accepted-endpoint continuation state')
  call test_end()

  call cleanup_output_dir(out_dir_disabled)
  call cleanup_output_dir(out_dir_ledger)
  call cleanup_output_dir(out_dir_no_photo)
  call cleanup_output_dir(out_dir_matching)
  call cleanup_output_dir(out_dir_matching_online)

  call test_summary()

contains

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

  subroutine scan_matching_plane_receipts(summary_path, found)
    character(len=*), intent(in) :: summary_path
    logical, intent(out) :: found(13)
    integer :: summary_unit, summary_ios
    character(len=2048) :: summary_line

    found = .false.
    open (newunit=summary_unit, file=trim(summary_path), status='old', action='read', iostat=summary_ios)
    if (summary_ios /= 0) error stop 'failed to open matching-plane receipt fixture'
    do
      read (summary_unit, '(A)', iostat=summary_ios) summary_line
      if (summary_ios /= 0) exit
      found(1) = found(1) .or. trim(summary_line) == 'surface_current_model_response_backend=table'
      found(2) = found(2) .or. &
                 trim(summary_line) == &
                 'surface_current_model_response_table_path=examples/matching_plane_response_synthetic.csv'
      found(3) = found(3) .or. &
                 (index(summary_line, 'surface_current_model_response_content_fingerprint=') == 1 .and. &
                  len_trim(summary_line) == &
                  len('surface_current_model_response_content_fingerprint=') + 16)
      found(4) = found(4) .or. &
                 summary_real_equals( &
                 summary_line, 'surface_current_model_matching_plane_z_m', 1.0e-3_dp, 1.0e-16_dp &
                 )
      found(5) = found(5) .or. &
                 trim(summary_line) == 'surface_current_model_electron_species=solar_wind_electron'
      found(6) = found(6) .or. trim(summary_line) == 'surface_current_model_ion_species=solar_wind_ion'
      found(7) = found(7) .or. trim(summary_line) == 'surface_current_model_photoelectron_species=photoelectron'
      found(8) = found(8) .or. &
                 summary_real_equals( &
                 summary_line, 'surface_current_model_coupling_rtol', 1.0e-4_dp, 1.0e-16_dp &
                 )
      found(9) = found(9) .or. index(summary_line, 'surface_current_model_coupling_atol=') == 1
      found(10) = found(10) .or. trim(summary_line) == 'surface_current_model_coupling_max_iterations=20'
      found(11) = found(11) .or. &
                  summary_real_equals( &
                  summary_line, 'surface_current_model_coupling_relaxation', 0.5_dp, 1.0e-14_dp &
                  )
      found(12) = found(12) .or. &
                  trim(summary_line) == 'surface_current_model_dynamic_state_source=accepted_batch_fixed_point'
      found(13) = found(13) .or. trim(summary_line) == 'surface_current_model_implicit_zero_mode=F'
    end do
    close (summary_unit)
  end subroutine scan_matching_plane_receipts

  subroutine scan_online_matching_plane_receipts(summary_path, found)
    character(len=*), intent(in) :: summary_path
    logical, intent(out) :: found(9)
    integer :: summary_unit, summary_ios
    character(len=2048) :: summary_line

    found = .false.
    open (newunit=summary_unit, file=trim(summary_path), status='old', action='read', iostat=summary_ios)
    if (summary_ios /= 0) error stop 'failed to open online matching-plane receipt fixture'
    do
      read (summary_unit, '(A)', iostat=summary_ios) summary_line
      if (summary_ios /= 0) exit
      found(1) = found(1) .or. trim(summary_line) == 'surface_current_model_response_backend=zhao_online'
      found(2) = found(2) .or. &
                 trim(summary_line) == 'surface_current_model_response_contract=matching_plane_zhao_online_v1'
      found(3) = found(3) .or. trim(summary_line) == 'surface_current_model_zhao_branch=auto'
      found(4) = found(4) .or. &
                 trim(summary_line) == 'surface_current_model_zhao_root_selection=require_unique'
      found(5) = found(5) .or. &
                 trim(summary_line) == 'surface_current_model_outer_solver=charge_driven_finite_h_sagdeev'
      found(6) = found(6) .or. &
                 trim(summary_line) == &
                 'surface_current_model_photoelectron_closure=moment_matched_half_maxwellian'
      found(7) = found(7) .or. &
                 trim(summary_line) == 'surface_current_model_ambient_outward_feedback=transparent'
      found(8) = found(8) .or. trim(summary_line) == 'surface_current_model_outer_solver_state=stateless'
      found(9) = found(9) .or. trim(summary_line) == 'surface_current_model_implicit_zero_mode=T'
    end do
    close (summary_unit)
  end subroutine scan_online_matching_plane_receipts

  subroutine assert_resolved_boundary_summary(path, inflow_map, open_model)
    character(len=*), intent(in) :: path, inflow_map, open_model
    logical :: saw_inflow, saw_open
    integer :: unit, read_status
    character(len=512) :: summary_line

    saw_inflow = .false.
    saw_open = .false.
    open (newunit=unit, file=path, status='old', action='read', iostat=read_status)
    if (read_status /= 0) error stop 'failed to open resolved local boundary summary fixture'
    do
      read (unit, '(A)', iostat=read_status) summary_line
      if (read_status /= 0) exit
      saw_inflow = saw_inflow .or. trim(summary_line) == 'reservoir_inflow_map='//trim(inflow_map)
      saw_open = saw_open .or. trim(summary_line) == 'particle_ordinary_open_model='//trim(open_model)
    end do
    close (unit)
    call assert_true(saw_inflow, 'summary should record the resolved reservoir inflow map')
    call assert_true(saw_open, 'summary should record the resolved ordinary open model')
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

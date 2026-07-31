!> CSV 出力テスト: write_mesh_potential disabled 時の挙動検証。
program test_output_writer_io
  use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
  use bem_kinds, only: dp, i32
  use bem_mesh, only: init_mesh
  use bem_output_writer, only: ensure_output_dir, open_top_reference_history_writer, &
                               write_result_files, write_top_reference_history_snapshot
  use bem_app_config, only: app_config, default_app_config
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
  logical :: top_history_opened, saw_top_available, saw_top_definition, saw_top_last_batch, saw_top_mean
  integer :: literal_unit, ios, top_history_unit
  integer(i32) :: top_batch, top_sample_n
  real(dp) :: top_time, top_z, top_mean, top_std, top_min, top_max
  character(len=512) :: line
  character(len=*), parameter :: out_dir_disabled = 'test_output_writer_io_disabled_tmp'
  character(len=*), parameter :: out_dir_ledger = 'test_output_writer_io_ledger_tmp'
  character(len=*), parameter :: literal_parent = 'test_output_writer_io_literal_tmp'
  character(len=*), parameter :: marker_path = 'test_output_writer_io_shell_marker_tmp'
  character(len=*), parameter :: literal_dir = &
                                 literal_parent//'/space $(touch '//marker_path//'); "double" ''single'''
  character(len=*), parameter :: expanded_dir = literal_parent//'/space ; double ''single'''

  interface
    integer(c_int) function c_rmdir(path) bind(C, name='rmdir')
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: path(*)
    end function c_rmdir
  end interface

  call test_init(5)

  stats = sim_stats()

  call cleanup_output_dir(out_dir_disabled)
  call cleanup_output_dir(out_dir_ledger)

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

  call test_begin('write_mesh_potential_disabled')
  call build_two_element_mesh(mesh)
  mesh%q_elem = [2.0d-12, -1.0d-12]

  call default_app_config(cfg)
  cfg%sim%batch_duration = 1.0_dp
  cfg%output_dir = out_dir_disabled
  cfg%write_mesh_potential = .false.
  call write_result_files(out_dir_disabled, mesh, stats, cfg)
  inquire (file=trim(out_dir_disabled)//'/mesh_potential.csv', exist=exists)
  call assert_true(.not. exists, 'mesh_potential.csv should not be written when output.write_mesh_potential=false')
  call test_end()

  call test_begin('resolved_local_boundary_receipt')
  call default_app_config(cfg)
  cfg%sim%reservoir_potential_model = 'infinity_barrier'
  cfg%sim%open_boundary_model = 'potential_barrier'
  call write_result_files(out_dir_disabled, mesh, stats, cfg)
  call assert_resolved_boundary_summary( &
    out_dir_disabled//'/summary.txt', 'infinity_barrier', 'potential_barrier' &
    )

  call default_app_config(cfg)
  cfg%sim%batch_duration = 1.0_dp
  cfg%output_dir = out_dir_disabled
  cfg%write_mesh_potential = .false.
  call test_end()

  call test_begin('top_reference_history_and_summary')
  call default_app_config(cfg)
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
  call assert_equal_i32(top_batch, 4_i32, 'resumed top-reference history batch mismatch')
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
  open (newunit=literal_unit, file=out_dir_ledger//'/summary.txt', status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'failed to open summary metadata fixture'
  do
    read (literal_unit, '(A)', iostat=ios) line
    if (ios /= 0) exit
    saw_integrator = saw_integrator .or. index(line, 'particle_time_centering=same_time_midpoint_boris') > 0
    saw_residual = saw_residual .or. index(line, 'charge_ledger_residual_C=') > 0
    saw_schema = saw_schema .or. index(line, 'checkpoint_schema_version=6') > 0
    saw_model_fp = saw_model_fp .or. index(line, 'model_fingerprint=') > 0
    saw_mesh_fp = saw_mesh_fp .or. index(line, 'mesh_fingerprint=') > 0
    saw_species_fp = saw_species_fp .or. index(line, 'species_fingerprint=') > 0
    saw_ledger_stock = saw_ledger_stock .or. index(line, 'charge_ledger_local_flight_charge_before_C=') > 0
    saw_ledger_closure = saw_ledger_closure .or. &
                         index(line, 'charge_ledger_neutral_return_correction_C=') > 0
    saw_build_schema = saw_build_schema .or. index(line, 'build_info_schema_version=1') == 1
    saw_build_version = saw_build_version .or. index(line, 'build_version=') == 1
    saw_build_mode = saw_build_mode .or. index(line, 'build_version_mode=') == 1
    saw_source_commit = saw_source_commit .or. index(line, 'build_source_commit=') == 1
    saw_build_id = saw_build_id .or. index(line, 'build_id=') == 1
  end do
  close (literal_unit)
  open (newunit=literal_unit, file=out_dir_ledger//'/charge_ledger.csv', status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'failed to open charge ledger fixture'
  read (literal_unit, '(A)', iostat=ios) line
  close (literal_unit)
  saw_ledger_header = ios == 0 .and. index(line, 'species_idx') > 0 .and. &
                      index(line, 'discarded_unresolved_C') > 0 .and. &
                      index(line, 'neutral_return_correction_C') > 0 .and. &
                      index(line, 'neutral_return_weight_scale') > 0 .and. &
                      index(line, 'neutral_return_unresolved_fraction') > 0
  call assert_true(saw_integrator, 'summary should record the particle time-centering contract')
  call assert_true(saw_residual, 'summary should record the charge ledger residual')
  call assert_true(saw_schema, 'summary should record checkpoint schema v6')
  call assert_true(saw_model_fp .and. saw_mesh_fp .and. saw_species_fp, 'summary should record restart fingerprints')
  call assert_true(saw_build_schema .and. saw_build_version .and. saw_build_mode .and. saw_source_commit .and. saw_build_id, &
                   'summary should record executable build origin')
  call assert_true(saw_ledger_stock, 'summary should record restartable charge stocks')
  call assert_true(saw_ledger_closure, 'summary should record the neutral-return correction')
  call assert_true(saw_ledger_header, 'charge ledger CSV header mismatch')
  call test_end()

  call cleanup_output_dir(out_dir_disabled)
  call cleanup_output_dir(out_dir_ledger)

  call test_summary()

contains

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
      saw_mean = saw_mean .or. index(summary_line, 'top_reference_potential_mean_V=') == 1
    end do
    close (unit)
  end subroutine scan_top_reference_summary_fields

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
    call delete_file_if_exists(out_dir//'/top_reference_history.csv')
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

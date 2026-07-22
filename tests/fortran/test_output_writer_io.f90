!> CSV 出力テスト: write_mesh_potential disabled 時の挙動検証。
program test_output_writer_io
  use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
  use bem_kinds, only: dp, i32
  use bem_mesh, only: init_mesh
  use bem_output_writer, only: ensure_output_dir, write_result_files
  use bem_app_config, only: app_config, default_app_config
  use bem_types, only: mesh_type, sim_stats
  use bem_charge_ledger, only: charge_ledger_type
  use bem_electrostatic_snapshot, only: electrostatic_diagnostics_type
  use bem_outer_plasma_photoelectron, only: photoelectron_histogram_type, photoelectron_histogram_state_type
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, delete_file_if_exists, remove_empty_directory
  implicit none

  type(mesh_type) :: mesh
  type(app_config) :: cfg
  type(sim_stats) :: stats
  type(charge_ledger_type) :: ledger
  type(photoelectron_histogram_type) :: photo_batch
  type(photoelectron_histogram_state_type) :: photo_state
  type(electrostatic_diagnostics_type) :: electrostatic_diagnostics
  logical :: exists, literal_created, marker_created, saw_integrator, saw_residual, saw_ledger_header
  logical :: saw_schema, saw_model_fp, saw_mesh_fp, saw_species_fp, saw_ledger_stock, saw_photo_batch, saw_photo_flux
  logical :: saw_build_schema, saw_build_version, saw_build_mode, saw_source_commit, saw_build_id
  logical :: saw_queue_population, saw_queue_count, saw_queue_fingerprint
  logical :: saw_steady_start_mode, saw_steady_start_mesh_id
  integer :: literal_unit, ios
  character(len=512) :: line
  character(len=*), parameter :: out_dir_disabled = 'test_output_writer_io_disabled_tmp'
  character(len=*), parameter :: out_dir_ledger = 'test_output_writer_io_ledger_tmp'
  character(len=*), parameter :: out_dir_photo = 'test_output_writer_io_photo_tmp'
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
  call cleanup_output_dir(out_dir_photo)

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

  call test_begin('photoelectron_checkpoint_output')
  call photo_state%init(2_i32, 4.0_dp)
  call photo_state%begin_batch(photo_batch)
  call photo_batch%add(-1.0_dp, 2.0_dp, 3.0_dp, [1.0_dp, -2.0_dp, 1.0_dp])
  call photo_state%commit_batch(1_i32, photo_batch)
  stats%batches = 1_i32
  cfg%outer_plasma%photoelectron_histogram_enabled = .false.
  call write_result_files(out_dir_photo, mesh, stats, cfg, photoelectron_state=photo_state)
  inquire (file=out_dir_photo//'/photoelectron_histogram.csv', exist=exists)
  call assert_true(.not. exists, 'disabled photoelectron histogram must ignore a supplied ready state')
  saw_photo_batch = .false.
  open (newunit=literal_unit, file=out_dir_photo//'/summary.txt', status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'failed to open disabled photoelectron summary fixture'
  do
    read (literal_unit, '(A)', iostat=ios) line
    if (ios /= 0) exit
    saw_photo_batch = saw_photo_batch .or. index(line, 'photoelectron_last_completed_batch=') > 0
  end do
  close (literal_unit)
  call assert_true(.not. saw_photo_batch, 'disabled histogram must not add photoelectron summary fields')
  cfg%outer_plasma%photoelectron_histogram_enabled = .true.
  call write_result_files(out_dir_photo, mesh, stats, cfg, photoelectron_state=photo_state)
  inquire (file=out_dir_photo//'/photoelectron_histogram.csv', exist=exists)
  call assert_true(exists, 'photoelectron histogram checkpoint should be written')
  saw_photo_batch = .false.
  saw_photo_flux = .false.
  open (newunit=literal_unit, file=out_dir_photo//'/summary.txt', status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'failed to open photoelectron summary fixture'
  do
    read (literal_unit, '(A)', iostat=ios) line
    if (ios /= 0) exit
    saw_photo_batch = saw_photo_batch .or. index(line, 'photoelectron_last_completed_batch=1') > 0
    saw_photo_flux = saw_photo_flux .or. index(line, 'photoelectron_previous_signed_current_A=') > 0
  end do
  close (literal_unit)
  call assert_true(saw_photo_batch, 'summary should record photoelectron histogram batch ownership')
  call assert_true(saw_photo_flux, 'summary should record the outgoing photoelectron signed current')
  cfg%outer_plasma%photoelectron_histogram_enabled = .false.
  call test_end()

  call test_begin('queue_diagnostics_follow_enable_flag')
  electrostatic_diagnostics%outer_kinetic_closure = 'zhao_charge_driven'
  electrostatic_diagnostics%outer_zhao_branch = 'B'
  electrostatic_diagnostics%outer_photoelectron_population_fraction = 0.5_dp
  electrostatic_diagnostics%outer_queue_event_count = 2_i32
  electrostatic_diagnostics%outer_queue_fingerprint = 'FEDCBA9876543210'
  cfg%coupling%outer_queue_enabled = .false.
  call write_result_files( &
    out_dir_disabled, mesh, stats, cfg, electrostatic_diagnostics=electrostatic_diagnostics &
    )
  call scan_queue_summary_fields( &
    out_dir_disabled//'/summary.txt', saw_queue_population, saw_queue_count, saw_queue_fingerprint &
    )
  call assert_true(.not. saw_queue_population .and. .not. saw_queue_count .and. .not. saw_queue_fingerprint, &
                   'disabled queue must not write queue-specific diagnostics')
  cfg%coupling%outer_queue_enabled = .true.
  call write_result_files( &
    out_dir_disabled, mesh, stats, cfg, electrostatic_diagnostics=electrostatic_diagnostics &
    )
  call scan_queue_summary_fields( &
    out_dir_disabled//'/summary.txt', saw_queue_population, saw_queue_count, saw_queue_fingerprint &
    )
  call assert_true(saw_queue_population .and. saw_queue_count .and. saw_queue_fingerprint, &
                   'enabled queue must write queue-specific diagnostics')
  cfg%coupling%outer_queue_enabled = .false.
  call test_end()

  call test_begin('charge_ledger_and_model_metadata')
  call ledger%init(2_i32)
  call ledger%reset(1_i32)
  ledger%surface_charge_after = -3.0_dp
  ledger%local_flight_charge_before = -1.0_dp
  ledger%local_flight_charge_after = -2.0_dp
  ledger%injected_from_remote(1) = -3.0_dp
  ledger%absorbed_on_surface(1) = -3.0_dp
  ledger%injected_count(1) = 1
  ledger%absorbed_count(1) = 1
  cfg%output_dir = out_dir_ledger
  cfg%coupling%steady_start_mode = 'none'
  cfg%coupling%steady_start_mesh_id = 7_i32
  call write_result_files(out_dir_ledger, mesh, stats, cfg, charge_ledger=ledger)

  saw_integrator = .false.
  saw_residual = .false.
  saw_schema = .false.
  saw_model_fp = .false.
  saw_mesh_fp = .false.
  saw_species_fp = .false.
  saw_ledger_stock = .false.
  saw_build_schema = .false.
  saw_build_version = .false.
  saw_build_mode = .false.
  saw_source_commit = .false.
  saw_build_id = .false.
  saw_steady_start_mode = .false.
  saw_steady_start_mesh_id = .false.
  open (newunit=literal_unit, file=out_dir_ledger//'/summary.txt', status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'failed to open summary metadata fixture'
  do
    read (literal_unit, '(A)', iostat=ios) line
    if (ios /= 0) exit
    saw_integrator = saw_integrator .or. index(line, 'particle_time_centering=same_time_midpoint_boris') > 0
    saw_residual = saw_residual .or. index(line, 'charge_ledger_residual_C=') > 0
    saw_schema = saw_schema .or. index(line, 'checkpoint_schema_version=4') > 0
    saw_model_fp = saw_model_fp .or. index(line, 'model_fingerprint=') > 0
    saw_mesh_fp = saw_mesh_fp .or. index(line, 'mesh_fingerprint=') > 0
    saw_species_fp = saw_species_fp .or. index(line, 'species_fingerprint=') > 0
    saw_ledger_stock = saw_ledger_stock .or. index(line, 'charge_ledger_local_flight_charge_before_C=') > 0
    saw_build_schema = saw_build_schema .or. index(line, 'build_info_schema_version=1') == 1
    saw_build_version = saw_build_version .or. index(line, 'build_version=') == 1
    saw_build_mode = saw_build_mode .or. index(line, 'build_version_mode=') == 1
    saw_source_commit = saw_source_commit .or. index(line, 'build_source_commit=') == 1
    saw_build_id = saw_build_id .or. index(line, 'build_id=') == 1
    saw_steady_start_mode = saw_steady_start_mode .or. index(line, 'coupling_steady_start_mode=none') == 1
    saw_steady_start_mesh_id = saw_steady_start_mesh_id .or. index(line, 'coupling_steady_start_mesh_id=7') == 1
  end do
  close (literal_unit)
  open (newunit=literal_unit, file=out_dir_ledger//'/charge_ledger.csv', status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'failed to open charge ledger fixture'
  read (literal_unit, '(A)', iostat=ios) line
  close (literal_unit)
  saw_ledger_header = ios == 0 .and. index(line, 'species_idx') > 0 .and. index(line, 'discarded_unresolved_C') > 0
  call assert_true(saw_integrator, 'summary should record the particle time-centering contract')
  call assert_true(saw_residual, 'summary should record the charge ledger residual')
  call assert_true(saw_schema, 'summary should record checkpoint schema v4')
  call assert_true(saw_model_fp .and. saw_mesh_fp .and. saw_species_fp, 'summary should record restart fingerprints')
  call assert_true(saw_build_schema .and. saw_build_version .and. saw_build_mode .and. saw_source_commit .and. saw_build_id, &
                   'summary should record executable build origin')
  call assert_true(saw_steady_start_mode .and. saw_steady_start_mesh_id, &
                   'summary should record the steady-start configuration')
  call assert_true(saw_ledger_stock, 'summary should record restartable charge stocks')
  call assert_true(saw_ledger_header, 'charge ledger CSV header mismatch')
  call test_end()

  call cleanup_output_dir(out_dir_disabled)
  call cleanup_output_dir(out_dir_ledger)
  call cleanup_output_dir(out_dir_photo)

  call test_summary()

contains

  subroutine scan_queue_summary_fields(path, saw_population, saw_count, saw_fingerprint)
    character(len=*), intent(in) :: path
    logical, intent(out) :: saw_population, saw_count, saw_fingerprint
    integer :: unit, read_status
    character(len=512) :: summary_line

    saw_population = .false.
    saw_count = .false.
    saw_fingerprint = .false.
    open (newunit=unit, file=path, status='old', action='read', iostat=read_status)
    if (read_status /= 0) error stop 'failed to open queue summary fixture'
    do
      read (unit, '(A)', iostat=read_status) summary_line
      if (read_status /= 0) exit
      saw_population = saw_population .or. index(summary_line, 'outer_photoelectron_population_fraction=') == 1
      saw_count = saw_count .or. index(summary_line, 'outer_queue_event_count=') == 1
      saw_fingerprint = saw_fingerprint .or. index(summary_line, 'outer_queue_fingerprint=FEDCBA9876543210') == 1
    end do
    close (unit)
  end subroutine scan_queue_summary_fields

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
    call delete_file_if_exists(out_dir//'/photoelectron_histogram.csv')
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

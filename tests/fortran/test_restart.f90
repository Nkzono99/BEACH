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
  use bem_outer_plasma_photoelectron, only: photoelectron_histogram_type, photoelectron_coupling_state_type
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
  type(photoelectron_histogram_type) :: photo_batch
  type(photoelectron_coupling_state_type) :: photo_state, restored_photo_state
  logical :: has_restart, exists
  integer(i32) :: contract_status
  character(len=256) :: contract_message
  character(len=1024) :: rng_rank_path, residual_global_path
  character(len=*), parameter :: out_dir = 'test_restart_tmp'

  call delete_file_if_exists(out_dir//'/summary.txt')
  call delete_file_if_exists(out_dir//'/charges.csv')
  call delete_file_if_exists(out_dir//'/rng_state.txt')
  call delete_file_if_exists(out_dir//'/macro_residuals.csv')
  call delete_file_if_exists(out_dir//'/rng_state_rank00001.txt')
  call delete_file_if_exists(out_dir//'/macro_residuals_rank00001.csv')
  call delete_file_if_exists(out_dir//'/photoelectron_histogram.csv')
  call remove_empty_directory(out_dir)

  call test_init(5)

  call test_begin('missing_checkpoint')
  call build_test_mesh(mesh)
  call load_restart_checkpoint('test_restart_missing', mesh, stats, has_restart)
  call assert_true(.not. has_restart, 'missing checkpoint should not be reported as restart')
  call assert_equal_i32(stats%batches, 0_i32, 'missing checkpoint should keep stats at defaults')
  call test_end()

  call test_begin('schema_v2_contract_and_ledger_restore')
  call default_app_config(cfg)
  cfg%n_particle_species = 2_i32
  cfg%particle_species(1) = species_from_defaults()
  cfg%particle_species(1)%species_key = 'electron'
  cfg%particle_species(2) = species_from_defaults()
  cfg%particle_species(2)%species_key = 'ion'
  cfg%outer_plasma%photoelectron_closure = 'INDIVIDUAL_RETURN'
  cfg%outer_plasma%photoelectron_histogram_bins = 2_i32
  cfg%outer_plasma%photoelectron_histogram_energy_max = 4.0_dp
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
  call photo_state%init(2_i32, 4.0_dp)
  call photo_state%begin_batch(photo_batch)
  call photo_batch%add(-1.0_dp, 2.0_dp, 3.0_dp, [1.0_dp, -2.0_dp, 1.0_dp])
  call photo_state%commit_batch(1_i32, photo_batch)
  call photo_state%begin_batch(photo_batch)
  call photo_batch%add(-1.0_dp, 2.0_dp, 5.0_dp, [-1.0_dp, 1.0_dp, 2.0_dp])
  call photo_state%commit_batch(2_i32, photo_batch)
  call write_result_files(out_dir, mesh, stats, cfg, charge_ledger=ledger, photoelectron_state=photo_state)
  call write_rng_state_file(out_dir)
  call write_macro_residuals_file(out_dir, state)

  call build_test_mesh(mesh)
  call load_restart_checkpoint( &
    out_dir, mesh, stats, has_restart, state, app=cfg, charge_ledger=restored_ledger, &
    photoelectron_state=restored_photo_state &
    )
  call assert_true(has_restart, 'schema v2 checkpoint should load')
  call assert_equal_i32(restored_ledger%batch_count, 2_i32, 'ledger batch count mismatch')
  call assert_close_dp(restored_ledger%surface_charge_before, 1.0_dp, 1.0e-12_dp, 'ledger stock mismatch')
  call assert_allclose_1d( &
    restored_ledger%injected_from_remote, [-3.0_dp, 4.0_dp], 1.0e-12_dp, 'ledger flux restore mismatch' &
    )
  call assert_equal_i32(restored_photo_state%last_completed_batch, 2_i32, 'photoelectron batch restore mismatch')
  call assert_close_dp( &
    restored_photo_state%previous_batch%total_signed_charge(), -5.0_dp, 1.0e-12_dp, &
    'photoelectron previous-batch charge mismatch' &
    )
  call assert_close_dp( &
    restored_photo_state%cumulative%total_signed_charge(), -8.0_dp, 1.0e-12_dp, &
    'photoelectron cumulative charge mismatch' &
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
  call delete_file_if_exists(out_dir//'/photoelectron_histogram.csv')
  call remove_empty_directory(out_dir)

  call test_summary()

contains

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

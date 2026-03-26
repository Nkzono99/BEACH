!> 再開用チェックポイントの保存/復元を検証するテスト。
program test_restart
  use bem_kinds, only: dp, i32, i64
  use bem_mesh, only: init_mesh
  use bem_restart, only: load_restart_checkpoint, write_rng_state_file, write_macro_residuals_file, &
                         restart_rng_state_path, restart_macro_residual_path
  use bem_types, only: mesh_type, sim_stats, injection_state
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_equal_i64, assert_close_dp, assert_allclose_1d, &
                          delete_file_if_exists, ensure_directory, remove_empty_directory
  implicit none

  type(mesh_type) :: mesh
  type(sim_stats) :: stats
  type(injection_state) :: state
  logical :: has_restart, exists
  character(len=1024) :: rng_rank_path, residual_rank_path
  character(len=*), parameter :: out_dir = 'test_restart_tmp'

  call delete_file_if_exists(out_dir//'/summary.txt')
  call delete_file_if_exists(out_dir//'/charges.csv')
  call delete_file_if_exists(out_dir//'/rng_state.txt')
  call delete_file_if_exists(out_dir//'/macro_residuals.csv')
  call delete_file_if_exists(out_dir//'/rng_state_rank00001.txt')
  call delete_file_if_exists(out_dir//'/macro_residuals_rank00001.csv')
  call remove_empty_directory(out_dir)

  call test_init(4)

  call test_begin('missing_checkpoint')
  call build_test_mesh(mesh)
  call load_restart_checkpoint('test_restart_missing', mesh, stats, has_restart)
  call assert_true(.not. has_restart, 'missing checkpoint should not be reported as restart')
  call assert_equal_i32(stats%batches, 0_i32, 'missing checkpoint should keep stats at defaults')
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
  residual_rank_path = restart_macro_residual_path(out_dir, mpi_rank=1_i32, mpi_size=4_i32)
  call assert_true(trim(rng_rank_path) == trim(out_dir)//'/rng_state_rank00001.txt', 'ranked rng path mismatch')
  call assert_true( &
    trim(residual_rank_path) == trim(out_dir)//'/macro_residuals_rank00001.csv', 'ranked residual path mismatch' &
    )

  call write_rng_state_file(out_dir, mpi_rank=1_i32, mpi_size=4_i32)
  call write_macro_residuals_file(out_dir, state, mpi_rank=1_i32, mpi_size=4_i32)
  inquire (file=trim(rng_rank_path), exist=exists)
  call assert_true(exists, 'rng_state_rank00001.txt should be created')
  inquire (file=trim(residual_rank_path), exist=exists)
  call assert_true(exists, 'macro_residuals_rank00001.csv should be created')
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

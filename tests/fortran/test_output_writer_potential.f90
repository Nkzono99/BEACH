!> 事前計算済みメッシュ電位の CSV 出力契約を検証する。
program test_output_writer_potential
  use bem_kinds, only: dp, i32
  use bem_mesh, only: init_mesh
  use bem_output_writer, only: write_result_files
  use bem_app_config, only: app_config, default_app_config
  use bem_types, only: mesh_type, sim_stats
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, delete_file_if_exists, remove_empty_directory
  implicit none

  type(mesh_type) :: mesh
  type(app_config) :: cfg
  type(sim_stats) :: stats
  real(dp), allocatable :: values(:)
  real(dp), parameter :: expected_potential_v(2) = [1.234567890123456d3, -9.876543210987654d-7]
  character(len=*), parameter :: out_dir = 'test_output_writer_potential_tmp'

  call test_init(1)

  stats = sim_stats()

  call cleanup_output_dir(out_dir)

  call test_begin('precomputed_mesh_potential_is_written_verbatim')
  call build_two_element_mesh(mesh)
  mesh%q_elem = [2.0d-12, -1.0d-12]

  call default_app_config(cfg)
  cfg%output_dir = out_dir
  cfg%write_mesh_potential = .true.
  call write_result_files(out_dir, mesh, stats, cfg, mesh_potential_v=expected_potential_v)
  call read_potential_values(out_dir, values)
  call assert_equal_i32(int(size(values), i32), 2_i32, 'mesh_potential.csv row count mismatch')
  call assert_close_dp(values(1), expected_potential_v(1), 1.0d-12, 'mesh potential(1) mismatch')
  call assert_close_dp(values(2), expected_potential_v(2), 1.0d-20, 'mesh potential(2) mismatch')
  call test_end()

  call cleanup_output_dir(out_dir)

  call test_summary()

contains

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

  !> `mesh_potential.csv` を読み込む。
  subroutine read_potential_values(out_dir, values)
    character(len=*), intent(in) :: out_dir
    real(dp), allocatable, intent(out) :: values(:)
    character(len=256) :: header
    integer :: u, ios, file_row_count, idx, elem_idx
    logical :: exists

    inquire (file=trim(out_dir)//'/mesh_potential.csv', exist=exists)
    call assert_true(exists, 'mesh_potential.csv should exist')

    open (newunit=u, file=trim(out_dir)//'/mesh_potential.csv', status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to open mesh_potential.csv in test_output_writer_potential'

    read (u, '(a)', iostat=ios) header
    call assert_true(ios == 0, 'mesh_potential.csv header read failed')
    call assert_true(trim(header) == 'elem_idx,potential_V', 'mesh_potential.csv header mismatch')

    file_row_count = 0
    do
      read (u, *, iostat=ios)
      if (ios < 0) exit
      if (ios > 0) error stop 'failed to count mesh_potential.csv rows'
      file_row_count = file_row_count + 1
    end do

    rewind (u)
    read (u, '(a)', iostat=ios) header
    allocate (values(file_row_count))
    do idx = 1, file_row_count
      read (u, *, iostat=ios) elem_idx, values(idx)
      if (ios /= 0) error stop 'failed to parse mesh_potential.csv'
      call assert_equal_i32(int(elem_idx, i32), int(idx, i32), 'mesh_potential.csv element index mismatch')
    end do
    close (u)
  end subroutine read_potential_values

  !> writer テストの一時出力を削除する。
  subroutine cleanup_output_dir(out_dir)
    character(len=*), intent(in) :: out_dir

    call delete_file_if_exists(out_dir//'/summary.txt')
    call delete_file_if_exists(out_dir//'/charges.csv')
    call delete_file_if_exists(out_dir//'/mesh_potential.csv')
    call delete_file_if_exists(out_dir//'/checkpoint_complete.txt')
    call delete_file_if_exists(out_dir//'/checkpoint_complete.txt.tmp')
    call delete_file_if_exists(out_dir//'/mesh_triangles.csv')
    call delete_file_if_exists(out_dir//'/mesh_sources.csv')
    call remove_empty_directory(out_dir)
  end subroutine cleanup_output_dir

end program test_output_writer_potential

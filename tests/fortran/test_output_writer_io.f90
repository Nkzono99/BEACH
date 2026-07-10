!> CSV 出力テスト: write_mesh_potential disabled 時の挙動検証。
program test_output_writer_io
  use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
  use bem_kinds, only: dp, i32
  use bem_mesh, only: init_mesh
  use bem_output_writer, only: ensure_output_dir, write_result_files
  use bem_app_config, only: app_config, default_app_config
  use bem_types, only: mesh_type, sim_stats
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, delete_file_if_exists, remove_empty_directory
  implicit none

  type(mesh_type) :: mesh
  type(app_config) :: cfg
  type(sim_stats) :: stats
  logical :: exists, literal_created, marker_created
  integer :: literal_unit, ios
  character(len=*), parameter :: out_dir_disabled = 'test_output_writer_io_disabled_tmp'
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

  call test_init(2)

  stats = sim_stats()

  call cleanup_output_dir(out_dir_disabled)

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
  cfg%output_dir = out_dir_disabled
  cfg%write_mesh_potential = .false.
  call write_result_files(out_dir_disabled, mesh, stats, cfg)
  inquire (file=trim(out_dir_disabled)//'/mesh_potential.csv', exist=exists)
  call assert_true(.not. exists, 'mesh_potential.csv should not be written when output.write_mesh_potential=false')
  call test_end()

  call cleanup_output_dir(out_dir_disabled)

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

  !> writer テストの一時出力を削除する。
  subroutine cleanup_output_dir(out_dir)
    character(len=*), intent(in) :: out_dir

    call delete_file_if_exists(out_dir//'/summary.txt')
    call delete_file_if_exists(out_dir//'/charges.csv')
    call delete_file_if_exists(out_dir//'/mesh_potential.csv')
    call delete_file_if_exists(out_dir//'/mesh_triangles.csv')
    call delete_file_if_exists(out_dir//'/mesh_sources.csv')
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

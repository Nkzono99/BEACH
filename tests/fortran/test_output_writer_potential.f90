!> 電位計算 + CSV 出力テスト: free-space / periodic2 / FMM core のメッシュ電位検証。
program test_output_writer_potential
  use bem_kinds, only: dp, i32
  use bem_mesh, only: init_mesh
  use bem_output_writer, only: write_result_files
  use bem_app_config, only: app_config, default_app_config
  use bem_field_solver, only: field_solver_type
  use bem_panel_surface_sides, only: resolve_panel_surface_sides, panel_surface_side_ok
  use bem_physics_config_types, only: normalize_legacy_physics_config
  use bem_types, only: mesh_type, sim_stats, sim_config, bc_open, bc_periodic
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, delete_file_if_exists, remove_empty_directory
  implicit none

  type(mesh_type) :: mesh
  type(app_config) :: cfg
  type(sim_stats) :: stats
  type(field_solver_type) :: solver
  real(dp), allocatable :: values(:), potential_v(:)
  real(dp) :: free_reference(2), periodic_reference
  character(len=*), parameter :: out_dir_free = 'test_output_writer_pot_free_tmp'
  character(len=*), parameter :: out_dir_periodic = 'test_output_writer_pot_periodic_tmp'

  call test_init(4)

  stats = sim_stats()

  call cleanup_output_dir(out_dir_free)
  call cleanup_output_dir(out_dir_periodic)

  ! --- free-space mesh potential test ---
  call test_begin('free_space_mesh_potential')
  call build_two_element_mesh(mesh)
  mesh%q_elem = [2.0d-12, -1.0d-12]

  call default_app_config(cfg)
  cfg%output_dir = out_dir_free
  cfg%write_mesh_potential = .true.
  solver = field_solver_type()
  call solver%init(mesh, cfg%sim)
  call solver%refresh(mesh)
  allocate (potential_v(mesh%nelem))
  call solver%compute_mesh_potential(mesh, cfg%sim, potential_v)
  free_reference = potential_v
  call write_result_files(out_dir_free, mesh, stats, cfg, mesh_potential_v=potential_v)
  call read_potential_values(out_dir_free, values)
  call assert_equal_i32(int(size(values), i32), 2_i32, 'free-space mesh_potential.csv row count mismatch')
  call assert_close_dp(values(1), potential_v(1), 1.0d-12, 'free-space potential(1) mismatch')
  call assert_close_dp(values(2), potential_v(2), 1.0d-12, 'free-space potential(2) mismatch')
  deallocate (potential_v)
  call test_end()

  call test_begin('free_space_mesh_potential_length_normalized')
  call build_two_element_mesh(mesh)
  mesh%q_elem = [2.0d-12, -1.0d-12]
  call default_app_config(cfg)
  cfg%sim%field_normalization = 'length'
  cfg%sim%field_length_scale = 2.0d0
  call normalize_legacy_physics_config( &
    cfg%sim, cfg%field, cfg%periodic2, cfg%panel &
    )
  solver = field_solver_type()
  call solver%init(mesh, cfg%sim)
  call solver%refresh(mesh)
  allocate (potential_v(mesh%nelem))
  call solver%compute_mesh_potential(mesh, cfg%sim, potential_v)
  call assert_close_dp(potential_v(1), free_reference(1), 1.0d-9, 'normalized potential(1) mismatch')
  call assert_close_dp(potential_v(2), free_reference(2), 1.0d-9, 'normalized potential(2) mismatch')
  deallocate (potential_v)
  call test_end()

  ! --- periodic mesh potential test ---
  call test_begin('periodic_mesh_potential')
  call build_single_element_mesh(mesh)
  mesh%q_elem = [1.0d-12]
  call default_app_config(cfg)
  cfg%output_dir = out_dir_periodic
  cfg%write_mesh_potential = .true.
  cfg%sim%field_solver = 'fmm'
  cfg%sim%field_bc_mode = 'periodic2'
  cfg%sim%field_periodic_far_correction = 'none'
  cfg%sim%field_periodic_image_layers = 1_i32
  cfg%sim%field_periodic_ewald_layers = 1_i32
  cfg%sim%use_box = .true.
  cfg%sim%box_min = [0.0d0, 0.0d0, -1.0d0]
  cfg%sim%box_max = [1.0d0, 1.0d0, 1.0d0]
  cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
  call normalize_legacy_physics_config( &
    cfg%sim, cfg%field, cfg%periodic2, cfg%panel &
    )
  solver = field_solver_type()
  call solver%init(mesh, cfg%sim)
  call solver%refresh(mesh)
  allocate (potential_v(mesh%nelem))
  call solver%compute_mesh_potential(mesh, cfg%sim, potential_v)
  periodic_reference = potential_v(1)
  call write_result_files(out_dir_periodic, mesh, stats, cfg, mesh_potential_v=potential_v)
  call read_potential_values(out_dir_periodic, values)
  call assert_equal_i32(int(size(values), i32), 1_i32, 'periodic mesh_potential.csv row count mismatch')
  call assert_close_dp(values(1), periodic_reference, 5.0d-9, 'periodic potential mismatch')
  call test_end()

  ! --- FMM core mesh potential test ---
  call test_begin('fmm_core_mesh_potential')
  call test_fmm_core_mesh_potential(mesh, cfg%sim, periodic_reference, values)
  call test_end()
  deallocate (potential_v)

  call cleanup_output_dir(out_dir_free)
  call cleanup_output_dir(out_dir_periodic)

  call test_summary()

contains

  !> 2 要素メッシュを初期化する。
  subroutine build_two_element_mesh(mesh)
    type(mesh_type), intent(out) :: mesh
    real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2)
    integer(i32) :: status
    character(len=128) :: message

    v0(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    v1(:, 1) = [1.0d0, 0.0d0, 0.0d0]
    v2(:, 1) = [0.0d0, 1.0d0, 0.0d0]
    v0(:, 2) = [0.0d0, 0.0d0, 1.0d0]
    v1(:, 2) = [1.0d0, 0.0d0, 1.0d0]
    v2(:, 2) = [0.0d0, 1.0d0, 1.0d0]
    call init_mesh(mesh, v0, v1, v2)
    call resolve_panel_surface_sides(mesh, 'normal_plus', status, message)
    if (status /= panel_surface_side_ok) error stop 'two-element panel side setup failed: '//trim(message)
  end subroutine build_two_element_mesh

  !> 1 要素メッシュを初期化する。
  subroutine build_single_element_mesh(mesh)
    type(mesh_type), intent(out) :: mesh
    real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1)
    integer(i32) :: status
    character(len=128) :: message

    v0(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    v1(:, 1) = [1.0d0, 0.0d0, 0.0d0]
    v2(:, 1) = [0.0d0, 1.0d0, 0.0d0]
    call init_mesh(mesh, v0, v1, v2)
    call resolve_panel_surface_sides(mesh, 'normal_plus', status, message)
    if (status /= panel_surface_side_ok) error stop 'single-element panel side setup failed: '//trim(message)
  end subroutine build_single_element_mesh

  !> field_solver 経由で FMM メッシュ電位を計算し検証する。
  subroutine test_fmm_core_mesh_potential(mesh, sim, expected_phi, values)
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    real(dp), intent(in) :: expected_phi
    real(dp), allocatable, intent(out) :: values(:)
    type(field_solver_type) :: fmm_solver

    fmm_solver = field_solver_type()
    call fmm_solver%init(mesh, sim)
    call fmm_solver%refresh(mesh)
    allocate (values(mesh%nelem))
    call fmm_solver%compute_mesh_potential(mesh, sim, values)
    call assert_close_dp(values(1), expected_phi, 5.0d-9, 'periodic FMM mesh potential mismatch')
  end subroutine test_fmm_core_mesh_potential

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
    call delete_file_if_exists(out_dir//'/mesh_triangles.csv')
    call delete_file_if_exists(out_dir//'/mesh_sources.csv')
    call remove_empty_directory(out_dir)
  end subroutine cleanup_output_dir

end program test_output_writer_potential

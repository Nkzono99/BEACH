!> 最終出力 writer の CSV 生成と mesh 電位出力を検証するテスト。
program test_output_writer
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_mesh, only: init_mesh
  use bem_output_writer, only: write_result_files
  use bem_app_config, only: app_config, default_app_config
  use bem_field_solver, only: field_solver_type
  use bem_types, only: mesh_type, sim_stats, bc_open, bc_periodic
  use test_support, only: &
    assert_true, assert_equal_i32, assert_close_dp, delete_file_if_exists, remove_empty_directory
  implicit none

  type(mesh_type) :: mesh
  type(app_config) :: cfg
  type(sim_stats) :: stats
  type(field_solver_type) :: solver
  logical :: exists
  real(dp), parameter :: pi_dp = acos(-1.0d0)
  real(dp), allocatable :: values(:), potential_v(:)
  character(len=*), parameter :: out_dir_disabled = 'test_output_writer_disabled_tmp'
  character(len=*), parameter :: out_dir_free = 'test_output_writer_free_tmp'
  character(len=*), parameter :: out_dir_periodic = 'test_output_writer_periodic_tmp'

  stats = sim_stats()

  call cleanup_output_dir(out_dir_disabled)
  call cleanup_output_dir(out_dir_free)
  call cleanup_output_dir(out_dir_periodic)

  call build_two_element_mesh(mesh)
  mesh%q_elem = [2.0d-12, -1.0d-12]

  call default_app_config(cfg)
  cfg%output_dir = out_dir_disabled
  cfg%write_mesh_potential = .false.
  call write_result_files(out_dir_disabled, mesh, stats, cfg)
  inquire (file=trim(out_dir_disabled)//'/mesh_potential.csv', exist=exists)
  call assert_true(.not. exists, 'mesh_potential.csv should not be written when output.write_mesh_potential=false')

  call default_app_config(cfg)
  cfg%output_dir = out_dir_free
  cfg%write_mesh_potential = .true.
  cfg%sim%softening = 0.0d0
  solver = field_solver_type()
  call solver%init(mesh, cfg%sim)
  call solver%refresh(mesh)
  allocate (potential_v(mesh%nelem))
  call solver%compute_mesh_potential(mesh, cfg%sim, potential_v)
  call write_result_files(out_dir_free, mesh, stats, cfg, mesh_potential_v=potential_v)
  deallocate (potential_v)
  call read_potential_values(out_dir_free, values)
  call assert_equal_i32(int(size(values), i32), 2_i32, 'free-space mesh_potential.csv row count mismatch')
  call assert_close_dp(values(1), expected_free_potential_1(), 1.0d-9, 'free-space potential(1) mismatch')
  call assert_close_dp(values(2), expected_free_potential_2(), 1.0d-9, 'free-space potential(2) mismatch')

  call build_single_element_mesh(mesh)
  mesh%q_elem = [1.0d-12]
  call default_app_config(cfg)
  cfg%output_dir = out_dir_periodic
  cfg%write_mesh_potential = .true.
  cfg%sim%softening = 0.0d0
  cfg%sim%field_solver = 'fmm'
  cfg%sim%field_bc_mode = 'periodic2'
  cfg%sim%field_periodic_far_correction = 'm2l_root_oracle'
  cfg%sim%field_periodic_image_layers = 1_i32
  cfg%sim%field_periodic_ewald_layers = 1_i32
  cfg%sim%use_box = .true.
  cfg%sim%box_min = [0.0d0, 0.0d0, -1.0d0]
  cfg%sim%box_max = [1.0d0, 1.0d0, 1.0d0]
  cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
  solver = field_solver_type()
  call solver%init(mesh, cfg%sim)
  call solver%refresh(mesh)
  allocate (potential_v(mesh%nelem))
  call solver%compute_mesh_potential(mesh, cfg%sim, potential_v)
  call write_result_files(out_dir_periodic, mesh, stats, cfg, mesh_potential_v=potential_v)
  call read_potential_values(out_dir_periodic, values)
  call assert_equal_i32(int(size(values), i32), 1_i32, 'periodic mesh_potential.csv row count mismatch')
  call assert_close_dp(values(1), expected_periodic_potential(), 1.0d-9, 'periodic potential mismatch')

  call test_fmm_core_mesh_potential(mesh, cfg%sim, expected_periodic_potential(), values)
  deallocate (potential_v)

  call cleanup_output_dir(out_dir_disabled)
  call cleanup_output_dir(out_dir_free)
  call cleanup_output_dir(out_dir_periodic)

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

  !> 1 要素メッシュを初期化する。
  subroutine build_single_element_mesh(mesh)
    type(mesh_type), intent(out) :: mesh
    real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1)

    v0(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    v1(:, 1) = [1.0d0, 0.0d0, 0.0d0]
    v2(:, 1) = [0.0d0, 1.0d0, 0.0d0]
    call init_mesh(mesh, v0, v1, v2)
  end subroutine build_single_element_mesh

  !> field_solver 経由で FMM メッシュ電位を計算し検証する。
  subroutine test_fmm_core_mesh_potential(mesh, sim, expected_phi, values)
    use bem_types, only: sim_config
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
    call assert_close_dp(values(1), expected_phi, 1.0d-9, 'periodic FMM mesh potential mismatch')
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
    if (ios /= 0) error stop 'failed to open mesh_potential.csv in test_output_writer'

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

  !> free-space の 1 要素目期待電位。
  real(dp) function expected_free_potential_1()
    real(dp) :: self_coeff

    self_coeff = 2.0d0*sqrt(2.0d0*pi_dp)
    expected_free_potential_1 = k_coulomb*(self_coeff*2.0d-12 - 1.0d-12)
  end function expected_free_potential_1

  !> free-space の 2 要素目期待電位。
  real(dp) function expected_free_potential_2()
    real(dp) :: self_coeff

    self_coeff = 2.0d0*sqrt(2.0d0*pi_dp)
    expected_free_potential_2 = k_coulomb*(-self_coeff*1.0d-12 + 2.0d-12)
  end function expected_free_potential_2

  !> periodic2 image shell + oracle residual つき 1 要素期待電位。
  real(dp) function expected_periodic_potential()
    real(dp), parameter :: q = 1.0d-12
    real(dp), parameter :: alpha = 0.6d0
    real(dp) :: self_coeff

    self_coeff = 2.0d0*sqrt(2.0d0*pi_dp)
    expected_periodic_potential = k_coulomb*q*(self_coeff + direct_shell_sum() + oracle_real_space_correction(alpha) + &
                                               oracle_reciprocal_correction(alpha) + oracle_k0_correction(alpha))
  end function expected_periodic_potential

  real(dp) function direct_shell_sum()
    direct_shell_sum = 4.0d0 + 2.0d0*sqrt(2.0d0)
  end function direct_shell_sum

  real(dp) function oracle_real_space_correction(alpha)
    real(dp), intent(in) :: alpha
    integer :: ix, iy
    real(dp) :: r

    oracle_real_space_correction = 0.0d0
    do ix = -2, 2
      do iy = -2, 2
        if (ix == 0 .and. iy == 0) cycle
        r = sqrt(real(ix*ix + iy*iy, dp))
        oracle_real_space_correction = oracle_real_space_correction + erfc(alpha*r)/r
      end do
    end do
    do ix = -1, 1
      do iy = -1, 1
        if (ix == 0 .and. iy == 0) cycle
        r = sqrt(real(ix*ix + iy*iy, dp))
        oracle_real_space_correction = oracle_real_space_correction - 1.0d0/r
      end do
    end do
  end function oracle_real_space_correction

  real(dp) function oracle_reciprocal_correction(alpha)
    real(dp), intent(in) :: alpha
    integer :: h1, h2
    real(dp) :: kmag

    oracle_reciprocal_correction = 0.0d0
    do h1 = -1, 1
      do h2 = -1, 1
        if (h1 == 0 .and. h2 == 0) cycle
        kmag = 2.0d0*pi_dp*sqrt(real(h1*h1 + h2*h2, dp))
        oracle_reciprocal_correction = oracle_reciprocal_correction + (pi_dp/kmag)*2.0d0*erfc(0.5d0*kmag/alpha)
      end do
    end do
  end function oracle_reciprocal_correction

  real(dp) function oracle_k0_correction(alpha)
    real(dp), intent(in) :: alpha

    oracle_k0_correction = -2.0d0*pi_dp/(alpha*sqrt(pi_dp))
  end function oracle_k0_correction

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

end program test_output_writer

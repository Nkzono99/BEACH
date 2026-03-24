!> treecode/FMM ソルバモード自動選択とパラメータチューニングの検証テスト。
program test_dynamics_field_solver
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_types, only: mesh_type, sim_config, bc_open, bc_periodic
  use bem_templates, only: make_plane, make_sphere
  use bem_field, only: electric_field_at
  use bem_field_solver, only: field_solver_type
  use test_support, only: assert_true, assert_equal_i32, assert_close_dp
  implicit none

  call test_field_solver_auto_mode()
  call test_field_solver_explicit_mode_autotune()
  call test_treecode_field_accuracy()

contains

  subroutine test_field_solver_auto_mode()
    type(mesh_type) :: mesh_small, mesh_large
    type(mesh_type) :: mesh_mid, mesh_dense, mesh_huge
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim

    call make_plane(mesh_small, size_x=1.0d0, size_y=1.0d0, nx=1_i32, ny=1_i32, center=[0.0d0, 0.0d0, 0.0d0])
    sim = sim_config()
    sim%field_solver = 'AUTO'
    sim%tree_min_nelem = 64_i32
    call solver%init(mesh_small, sim)
    call assert_true(trim(solver%mode) == 'direct', 'auto solver should select direct for small meshes')

    call make_sphere(mesh_large, radius=0.6d0, n_lon=24_i32, n_lat=12_i32, center=[0.0d0, 0.0d0, 0.0d0])
    sim = sim_config()
    sim%field_solver = 'AUTO'
    sim%tree_min_nelem = 64_i32
    sim%tree_theta = 0.95d0
    sim%tree_leaf_max = 1_i32
    call solver%init(mesh_large, sim)
    call assert_true(trim(solver%mode) == 'treecode', 'auto solver should select treecode for dense meshes')
    call assert_close_dp(solver%theta, 0.40d0, 1.0d-12, 'auto theta mismatch for nelem<1500')
    call assert_equal_i32(solver%leaf_max, 12_i32, 'auto leaf_max mismatch for nelem<1500')

    call make_plane(mesh_mid, size_x=1.0d0, size_y=1.0d0, nx=30_i32, ny=30_i32, center=[0.0d0, 0.0d0, 0.0d0])
    sim = sim_config()
    sim%field_solver = 'AUTO'
    sim%tree_min_nelem = 1_i32
    call solver%init(mesh_mid, sim)
    call assert_close_dp(solver%theta, 0.50d0, 1.0d-12, 'auto theta mismatch for 1500<=nelem<10000')
    call assert_equal_i32(solver%leaf_max, 16_i32, 'auto leaf_max mismatch for 1500<=nelem<10000')

    call make_plane(mesh_dense, size_x=1.0d0, size_y=1.0d0, nx=80_i32, ny=80_i32, center=[0.0d0, 0.0d0, 0.0d0])
    sim = sim_config()
    sim%field_solver = 'AUTO'
    sim%tree_min_nelem = 1_i32
    call solver%init(mesh_dense, sim)
    call assert_close_dp(solver%theta, 0.58d0, 1.0d-12, 'auto theta mismatch for 10000<=nelem<50000')
    call assert_equal_i32(solver%leaf_max, 20_i32, 'auto leaf_max mismatch for 10000<=nelem<50000')

    call make_plane(mesh_huge, size_x=1.0d0, size_y=1.0d0, nx=160_i32, ny=160_i32, center=[0.0d0, 0.0d0, 0.0d0])
    sim = sim_config()
    sim%field_solver = 'AUTO'
    sim%tree_min_nelem = 1_i32
    call solver%init(mesh_huge, sim)
    call assert_close_dp(solver%theta, 0.65d0, 1.0d-12, 'auto theta mismatch for nelem>=50000')
    call assert_equal_i32(solver%leaf_max, 24_i32, 'auto leaf_max mismatch for nelem>=50000')
  end subroutine test_field_solver_auto_mode

  subroutine test_field_solver_explicit_mode_autotune()
    type(mesh_type) :: mesh_large
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim

    call make_sphere(mesh_large, radius=0.6d0, n_lon=24_i32, n_lat=12_i32, center=[0.0d0, 0.0d0, 0.0d0])

    sim = sim_config()
    sim%field_solver = 'treecode'
    sim%tree_theta = 0.95d0
    sim%tree_leaf_max = 3_i32
    call solver%init(mesh_large, sim)
    call assert_true(trim(solver%mode) == 'treecode', 'treecode mode mismatch')
    call assert_close_dp(solver%theta, 0.40d0, 1.0d-12, 'treecode autotune theta mismatch')
    call assert_equal_i32(solver%leaf_max, 12_i32, 'treecode autotune leaf_max mismatch')

    sim = sim_config()
    sim%field_solver = 'treecode'
    sim%tree_theta = 0.95d0
    sim%tree_leaf_max = 3_i32
    sim%has_tree_theta = .true.
    sim%has_tree_leaf_max = .true.
    call solver%init(mesh_large, sim)
    call assert_close_dp(solver%theta, 0.95d0, 1.0d-12, 'treecode explicit theta override mismatch')
    call assert_equal_i32(solver%leaf_max, 3_i32, 'treecode explicit leaf_max override mismatch')

    sim = sim_config()
    sim%field_solver = 'fmm'
    sim%tree_theta = 0.92d0
    sim%tree_leaf_max = 5_i32
    call solver%init(mesh_large, sim)
    call assert_true(trim(solver%mode) == 'fmm', 'fmm mode mismatch')
    call assert_close_dp(solver%theta, 0.40d0, 1.0d-12, 'fmm autotune theta mismatch')
    call assert_equal_i32(solver%leaf_max, 12_i32, 'fmm autotune leaf_max mismatch')

    sim = sim_config()
    sim%field_solver = 'fmm'
    sim%tree_theta = 0.92d0
    sim%tree_leaf_max = 5_i32
    sim%has_tree_theta = .true.
    sim%has_tree_leaf_max = .true.
    call solver%init(mesh_large, sim)
    call assert_close_dp(solver%theta, 0.92d0, 1.0d-12, 'fmm explicit theta override mismatch')
    call assert_equal_i32(solver%leaf_max, 5_i32, 'fmm explicit leaf_max override mismatch')
  end subroutine test_field_solver_explicit_mode_autotune

  subroutine test_treecode_field_accuracy()
    type(mesh_type) :: mesh_tree
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i, valid_count
    integer :: seed_size
    integer, allocatable :: seed(:)
    real(dp) :: r(3), e_direct(3), e_tree(3), max_rel_err
    real(dp) :: norm_direct, norm_diff, rel_err, norm_r

    call make_sphere(mesh_tree, radius=0.6d0, n_lon=24_i32, n_lat=12_i32, center=[0.0d0, 0.0d0, 0.0d0])
    mesh_tree%q_elem = 1.0d-12

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'treecode'
    sim%tree_theta = 0.10d0
    sim%tree_leaf_max = 6_i32
    sim%tree_min_nelem = 64_i32
    call solver%init(mesh_tree, sim)
    call solver%refresh(mesh_tree)

    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    do i = 1_i32, int(seed_size, i32)
      seed(i) = 314159 + 17*i
    end do
    call random_seed(put=seed)
    deallocate (seed)

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, 200_i32
      call random_number(r)
      r = 100.0d0*(r - 0.5d0)
      norm_r = sqrt(sum(r*r))
      if (norm_r < 20.0d0) then
        if (norm_r > 1.0d-12) then
          r = r*(20.0d0/norm_r)
        else
          r = [20.0d0, 0.0d0, 0.0d0]
        end if
      end if

      call electric_field_at(mesh_tree, r, sim%softening, e_direct)
      call solver%eval_e(mesh_tree, r, e_tree)

      norm_direct = sqrt(sum(e_direct*e_direct))
      if (norm_direct <= 1.0d-12) cycle
      norm_diff = sqrt(sum((e_tree - e_direct)*(e_tree - e_direct)))
      rel_err = norm_diff/norm_direct
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count >= 100_i32, 'treecode accuracy test has too few valid samples')
    call assert_true(max_rel_err <= 1.0d-3, 'treecode E relative error exceeds 1e-3')
  end subroutine test_treecode_field_accuracy

end program test_dynamics_field_solver

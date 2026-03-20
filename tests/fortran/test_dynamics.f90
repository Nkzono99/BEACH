!> 電場評価・Boris更新・衝突判定の基本動作を検証するテスト。
program test_dynamics
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_types, only: mesh_type, hit_info, sim_config, bc_open, bc_periodic
  use bem_mesh, only: init_mesh
  use bem_templates, only: make_plane, make_sphere
  use bem_field, only: electric_field_at
  use bem_field_solver, only: field_solver_type
  use bem_pusher, only: boris_push
  use bem_collision, only: find_first_hit, segment_triangle_intersect
  use test_support, only: assert_true, assert_equal_i32, assert_close_dp, assert_allclose_1d
  implicit none

  type(mesh_type) :: mesh_field, mesh_hit
  type(hit_info) :: hit
  real(dp) :: v0_field(3, 1), v1_field(3, 1), v2_field(3, 1), q0_field(1)
  real(dp) :: v0_hit(3, 2), v1_hit(3, 2), v2_hit(3, 2)
  real(dp) :: e(3), x_new(3), v_new(3), speed0, speed1
  real(dp) :: inv_r3, expected_ex

  v0_field(:, 1) = [1.0d0, 0.0d0, 0.0d0]
  v1_field(:, 1) = [0.0d0, 1.0d0, 0.0d0]
  v2_field(:, 1) = [-1.0d0, -1.0d0, 0.0d0]
  q0_field(1) = 2.0d-9
  call init_mesh(mesh_field, v0_field, v1_field, v2_field, q0=q0_field)

  call electric_field_at(mesh_field, [1.0d0, 0.0d0, 0.0d0], 0.5d0, e)
  inv_r3 = 1.0d0/(sqrt(1.25d0)*1.25d0)
  expected_ex = k_coulomb*q0_field(1)*inv_r3
  call assert_close_dp(e(1), expected_ex, abs(expected_ex)*1.0d-12, 'electric field Ex mismatch')
  call assert_close_dp(e(2), 0.0d0, 1.0d-20, 'electric field Ey should be zero')
  call assert_close_dp(e(3), 0.0d0, 1.0d-20, 'electric field Ez should be zero')

  call boris_push( &
    [0.0d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 0.0d0], &
    2.0d0, 1.0d0, 0.5d0, [1.0d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 0.0d0], &
    x_new, v_new &
    )
  call assert_allclose_1d(v_new, [1.0d0, 0.0d0, 0.0d0], 1.0d-12, 'boris (E only) velocity mismatch')
  call assert_allclose_1d(x_new, [0.5d0, 0.0d0, 0.0d0], 1.0d-12, 'boris (E only) position mismatch')

  call boris_push( &
    [0.0d0, 0.0d0, 0.0d0], [1.0d0, 2.0d0, -0.5d0], &
    1.0d0, 1.0d0, 0.1d0, [0.0d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 2.0d0], &
    x_new, v_new &
    )
  speed0 = sqrt(1.0d0 + 4.0d0 + 0.25d0)
  speed1 = sqrt(sum(v_new*v_new))
  call assert_close_dp(speed1, speed0, 1.0d-12, 'boris should preserve speed when E=0')

  v0_hit(:, 1) = [0.0d0, 0.0d0, 1.0d0]
  v1_hit(:, 1) = [1.0d0, 0.0d0, 1.0d0]
  v2_hit(:, 1) = [0.0d0, 1.0d0, 1.0d0]
  v0_hit(:, 2) = [0.0d0, 0.0d0, 0.0d0]
  v1_hit(:, 2) = [1.0d0, 0.0d0, 0.0d0]
  v2_hit(:, 2) = [0.0d0, 1.0d0, 0.0d0]
  call init_mesh(mesh_hit, v0_hit, v1_hit, v2_hit)

  call find_first_hit(mesh_hit, [0.2d0, 0.2d0, 2.0d0], [0.2d0, 0.2d0, -1.0d0], hit)
  call assert_true(hit%has_hit, 'segment should hit the mesh')
  call assert_equal_i32(hit%elem_idx, 1_i32, 'first hit should be the nearer triangle')
  call assert_close_dp(hit%t, 1.0d0/3.0d0, 1.0d-12, 'hit parameter mismatch')
  call assert_close_dp(hit%pos(3), 1.0d0, 1.0d-12, 'hit position mismatch')

  call find_first_hit(mesh_hit, [0.9d0, 0.9d0, 2.0d0], [0.9d0, 0.9d0, -1.0d0], hit)
  call assert_true(.not. hit%has_hit, 'segment outside triangle should not hit')

  call test_collision_grid_equivalence()
  call test_field_solver_auto_mode()
  call test_field_solver_explicit_mode_autotune()
  call test_treecode_field_accuracy()
  call test_fmm_field_accuracy()
  call test_fmm_free_box_dual_target_accuracy()
  call test_fmm_core_free_box_dual_target_accuracy()
  call test_fmm_periodic2_field_accuracy()
  call test_fmm_core_periodic2_field_accuracy()
  call test_fmm_periodic2_fallback_accuracy()
  call test_fmm_periodic2_image_layers_accuracy()
  call test_fmm_periodic2_m2l_cache_reuse()
  call test_fmm_periodic2_default_m2l_root_oracle_mode()

contains

  subroutine test_collision_grid_equivalence()
    type(mesh_type) :: mesh_grid
    integer :: i, seed_size
    integer, allocatable :: seed(:)
    real(dp) :: p0(3), p1(3)

    call make_plane(mesh_grid, size_x=1.0d0, size_y=1.0d0, nx=8_i32, ny=8_i32, center=[0.0d0, 0.0d0, 0.0d0])
    call assert_true(mesh_grid%use_collision_grid, 'collision grid should be enabled for dense mesh')

    call assert_hit_equivalent(mesh_grid, [1.8d0, 1.8d0, 1.0d0], [1.8d0, 1.8d0, -1.0d0], 'segment outside mesh AABB')
    call assert_hit_equivalent(mesh_grid, [-0.4d0, -0.4d0, 1.0d0], [0.4d0, -0.4d0, 1.0d0], 'axis parallel segment')
    call assert_hit_equivalent(mesh_grid, [0.1d0, 0.1d0, 0.0d0], [0.1d0, 0.1d0, 0.8d0], 'start point inside grid')
    call assert_hit_equivalent(mesh_grid, [-0.45d0, -0.35d0, 1.0d0], [0.45d0, 0.35d0, -1.0d0], 'multi-cell crossing')

    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    do i = 1, seed_size
      seed(i) = 12345 + 37*i
    end do
    call random_seed(put=seed)
    deallocate (seed)

    do i = 1, 200
      call random_number(p0)
      call random_number(p1)
      p0 = [1.4d0*(p0(1) - 0.5d0), 1.4d0*(p0(2) - 0.5d0), 3.0d0*(p0(3) - 0.5d0)]
      p1 = [1.4d0*(p1(1) - 0.5d0), 1.4d0*(p1(2) - 0.5d0), 3.0d0*(p1(3) - 0.5d0)]
      if (abs(p0(3)) < 1.0d-8) p0(3) = p0(3) + 1.0d-4
      if (abs(p1(3)) < 1.0d-8) p1(3) = p1(3) + 1.0d-4
      call assert_hit_equivalent(mesh_grid, p0, p1, 'random segment equivalence')
    end do
  end subroutine test_collision_grid_equivalence

  subroutine assert_hit_equivalent(mesh, p0, p1, label)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: p0(3), p1(3)
    character(len=*), intent(in) :: label
    type(hit_info) :: hit_fast, hit_ref

    call find_first_hit(mesh, p0, p1, hit_fast)
    call find_first_hit_bruteforce(mesh, p0, p1, hit_ref)

    call assert_true(hit_fast%has_hit .eqv. hit_ref%has_hit, trim(label)//': has_hit mismatch')
    if (.not. hit_ref%has_hit) return

    call assert_equal_i32(hit_fast%elem_idx, hit_ref%elem_idx, trim(label)//': elem_idx mismatch')
    call assert_close_dp(hit_fast%t, hit_ref%t, 1.0d-11, trim(label)//': t mismatch')
    call assert_allclose_1d(hit_fast%pos, hit_ref%pos, 1.0d-10, trim(label)//': position mismatch')
  end subroutine assert_hit_equivalent

  subroutine find_first_hit_bruteforce(mesh, p0, p1, hit)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: p0(3), p1(3)
    type(hit_info), intent(out) :: hit

    integer(i32) :: i
    logical :: ok
    real(dp) :: seg_min(3), seg_max(3), t, best_t, h(3)

    seg_min = min(p0, p1)
    seg_max = max(p0, p1)
    hit%has_hit = .false.
    hit%elem_idx = -1_i32
    hit%t = 0.0d0
    hit%pos = 0.0d0
    best_t = huge(1.0d0)

    do i = 1, mesh%nelem
      if (mesh%bb_max(1, i) < seg_min(1) .or. mesh%bb_min(1, i) > seg_max(1)) cycle
      if (mesh%bb_max(2, i) < seg_min(2) .or. mesh%bb_min(2, i) > seg_max(2)) cycle
      if (mesh%bb_max(3, i) < seg_min(3) .or. mesh%bb_min(3, i) > seg_max(3)) cycle
      call segment_triangle_intersect(p0, p1, mesh%v0(:, i), mesh%v1(:, i), mesh%v2(:, i), ok, t, h)
      if (.not. ok) cycle
      if (t < best_t) then
        best_t = t
        hit%has_hit = .true.
        hit%elem_idx = i
        hit%t = t
        hit%pos = h
      end if
    end do
  end subroutine find_first_hit_bruteforce

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

  subroutine test_fmm_field_accuracy()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i, valid_count
    integer :: seed_size
    integer, allocatable :: seed(:)
    real(dp) :: r(3), e_direct(3), e_fmm(3), max_rel_err
    real(dp) :: norm_direct, norm_diff, rel_err, norm_r

    call make_sphere(mesh_fmm, radius=0.6d0, n_lon=24_i32, n_lat=12_i32, center=[0.0d0, 0.0d0, 0.0d0])
    call assign_periodic_test_dipole_charges(mesh_fmm, 1.0d-12)

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'fmm'
    sim%tree_theta = 0.5d0
    sim%tree_leaf_max = 16_i32
    sim%tree_min_nelem = 64_i32
    call solver%init(mesh_fmm, sim)
    call solver%refresh(mesh_fmm)

    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    do i = 1_i32, int(seed_size, i32)
      seed(i) = 271828 + 29*i
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

      call electric_field_at(mesh_fmm, r, sim%softening, e_direct)
      call solver%eval_e(mesh_fmm, r, e_fmm)

      norm_direct = sqrt(sum(e_direct*e_direct))
      if (norm_direct <= 1.0d-12) cycle
      norm_diff = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))
      rel_err = norm_diff/norm_direct
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    write (*, '(A,I0,A,ES12.5)') 'test_fmm_field_accuracy: valid_count=', valid_count, ', max_rel_err=', max_rel_err
    call assert_true(valid_count >= 100_i32, 'fmm accuracy test has too few valid samples')
    call assert_true(max_rel_err <= 5.0d-3, 'fmm E relative error exceeds 5e-3')
  end subroutine test_fmm_field_accuracy

  subroutine test_fmm_free_box_dual_target_accuracy()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i, valid_count
    integer :: seed_size
    integer, allocatable :: seed(:)
    real(dp) :: r(3), e_direct(3), e_fmm(3), max_rel_err
    real(dp) :: norm_direct, norm_diff, rel_err, center_dist

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=18_i32, n_lat=10_i32, center=[0.5d0, 0.5d0, 0.0d0])
    mesh_fmm%q_elem = 1.0d-12

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'free'
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    call solver%init(mesh_fmm, sim)
    call solver%refresh(mesh_fmm)

    call assert_true(solver%target_tree_ready, 'fmm free/use_box should enable dual target tree')
    call assert_allclose_1d( &
      solver%target_node_center(:, 1_i32), 0.5d0*(sim%box_min + sim%box_max), 1.0d-12, &
      'fmm free/use_box target root center mismatch' &
      )
    call assert_allclose_1d( &
      solver%target_node_half_size(:, 1_i32), 0.5d0*(sim%box_max - sim%box_min), 1.0d-12, &
      'fmm free/use_box target root half-size mismatch' &
      )

    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    do i = 1_i32, int(seed_size, i32)
      seed(i) = 314159 + 31*i
    end do
    call random_seed(put=seed)
    deallocate (seed)

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, 200_i32
      call random_number(r)
      r = sim%box_min + r*(sim%box_max - sim%box_min)
      center_dist = sqrt((r(1) - 0.5d0)**2 + (r(2) - 0.5d0)**2 + r(3)**2)
      if (center_dist < 0.35d0) cycle

      call electric_field_at(mesh_fmm, r, sim%softening, e_direct)
      call solver%eval_e(mesh_fmm, r, e_fmm)

      norm_direct = sqrt(sum(e_direct*e_direct))
      if (norm_direct <= 1.0d-12) cycle
      norm_diff = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))
      rel_err = norm_diff/norm_direct
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    write (*, '(A,I0,A,ES12.5)') &
      'test_fmm_free_box_dual_target_accuracy: valid_count=', valid_count, ', max_rel_err=', max_rel_err
    call assert_true(valid_count >= 100_i32, 'fmm free/use_box test has too few valid samples')
    call assert_true(max_rel_err <= 2.5d-2, 'fmm free/use_box E relative error exceeds 2.5e-2')
  end subroutine test_fmm_free_box_dual_target_accuracy

  subroutine test_fmm_core_free_box_dual_target_accuracy()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i, valid_count
    integer :: seed_size
    integer, allocatable :: seed(:)
    real(dp) :: r(3), e_direct(3), e_fmm(3), max_rel_err
    real(dp) :: norm_direct, norm_diff, rel_err, center_dist

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=16_i32, n_lat=8_i32, center=[0.5d0, 0.5d0, 0.0d0])
    mesh_fmm%q_elem = 1.0d-12

    sim = sim_config()
    sim%softening = 0.0d0
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'free'
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    call solver%init(mesh_fmm, sim)
    call solver%refresh(mesh_fmm)

    call assert_true(solver%fmm_use_core, 'softening=0 free FMM should use the core path')
    call assert_true(solver%target_tree_ready, 'core free/use_box should expose target tree metadata')
    call assert_allclose_1d( &
      solver%target_node_center(:, 1_i32), 0.5d0*(sim%box_min + sim%box_max), 1.0d-12, &
      'core free/use_box target root center mismatch' &
      )
    call assert_allclose_1d( &
      solver%target_node_half_size(:, 1_i32), 0.5d0*(sim%box_max - sim%box_min), 1.0d-12, &
      'core free/use_box target root half-size mismatch' &
      )

    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    do i = 1_i32, int(seed_size, i32)
      seed(i) = 271829 + 17*i
    end do
    call random_seed(put=seed)
    deallocate (seed)

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, 120_i32
      call random_number(r)
      r = sim%box_min + r*(sim%box_max - sim%box_min)
      center_dist = sqrt((r(1) - 0.5d0)**2 + (r(2) - 0.5d0)**2 + r(3)**2)
      if (center_dist < 0.32d0) cycle

      call electric_field_at(mesh_fmm, r, 0.0d0, e_direct)
      call solver%eval_e(mesh_fmm, r, e_fmm)

      norm_direct = sqrt(sum(e_direct*e_direct))
      if (norm_direct <= 1.0d-12) cycle
      norm_diff = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))
      rel_err = norm_diff/norm_direct
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count >= 80_i32, 'core free/use_box test has too few valid samples')
    call assert_true(max_rel_err <= 5.0d-3, 'core free/use_box E relative error exceeds 5e-3')
  end subroutine test_fmm_core_free_box_dual_target_accuracy

  subroutine test_fmm_periodic2_field_accuracy()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i, valid_count, ref_layers
    real(dp) :: queries(3, 6)
    real(dp) :: r(3), e_direct(3), e_fmm(3), max_rel_err
    real(dp) :: norm_direct, norm_diff, rel_err

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=24_i32, n_lat=12_i32, center=[0.5d0, 0.5d0, 0.0d0])
    call assign_periodic_test_dipole_charges(mesh_fmm, 1.0d-12)

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'periodic2'
    sim%field_periodic_image_layers = 1_i32
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call solver%init(mesh_fmm, sim)
    call solver%refresh(mesh_fmm)
    ref_layers = sim%field_periodic_image_layers
    if (trim(solver%periodic_far_correction) == 'm2l_root_trunc') then
      ref_layers = ref_layers + solver%periodic_ewald_layers
    end if

    queries(:, 1) = [0.15d0, 0.15d0, -0.60d0]
    queries(:, 2) = [0.85d0, 0.20d0, -0.20d0]
    queries(:, 3) = [0.20d0, 0.80d0, 0.10d0]
    queries(:, 4) = [0.75d0, 0.75d0, 0.50d0]
    queries(:, 5) = [0.55d0, 0.35d0, -0.75d0]
    queries(:, 6) = [0.35d0, 0.60d0, 0.85d0]

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      r = queries(:, i)
      call electric_field_at_periodic2_images( &
        mesh_fmm, r, sim%softening, sim%box_min, sim%box_max, [1_i32, 2_i32], ref_layers, e_direct &
        )
      call solver%eval_e(mesh_fmm, r, e_fmm)

      norm_direct = sqrt(sum(e_direct*e_direct))
      if (norm_direct <= 1.0d-12) cycle
      norm_diff = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))
      rel_err = norm_diff/norm_direct
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    write (*, '(A,I0,A,ES12.5)') &
      'test_fmm_periodic2_field_accuracy(l2p): valid_count=', valid_count, ', max_rel_err=', max_rel_err
    call assert_true(valid_count == 6_i32, 'fmm periodic2 accuracy test lost valid samples')
    call assert_true(max_rel_err <= 2.0d-1, 'fmm periodic2 E relative error exceeds 2e-1')
  end subroutine test_fmm_periodic2_field_accuracy

  subroutine test_fmm_core_periodic2_field_accuracy()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i, valid_count, ref_layers
    real(dp) :: queries(3, 6)
    real(dp) :: r(3), e_direct(3), e_fmm(3), max_rel_err
    real(dp) :: norm_direct, norm_diff, rel_err

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=16_i32, n_lat=8_i32, center=[0.5d0, 0.5d0, 0.0d0])
    call assign_periodic_test_dipole_charges(mesh_fmm, 1.0d-12)

    sim = sim_config()
    sim%softening = 0.0d0
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'periodic2'
    sim%field_periodic_image_layers = 1_i32
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call solver%init(mesh_fmm, sim)
    call solver%refresh(mesh_fmm)

    call assert_true(solver%fmm_use_core, 'softening=0 periodic2 FMM should use the core path')
    call assert_true(solver%target_tree_ready, 'core periodic2 path should expose target tree metadata')
    call assert_true(solver%nleaf > 0_i32, 'core periodic2 path should expose leaf metadata')
    ref_layers = sim%field_periodic_image_layers
    if (trim(solver%periodic_far_correction) == 'm2l_root_trunc') then
      ref_layers = ref_layers + solver%periodic_ewald_layers
    end if

    queries(:, 1) = [0.15d0, 0.15d0, -0.60d0]
    queries(:, 2) = [0.85d0, 0.20d0, -0.20d0]
    queries(:, 3) = [0.20d0, 0.80d0, 0.10d0]
    queries(:, 4) = [0.75d0, 0.75d0, 0.50d0]
    queries(:, 5) = [0.55d0, 0.35d0, -0.75d0]
    queries(:, 6) = [0.35d0, 0.60d0, 0.85d0]

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      r = queries(:, i)
      call electric_field_at_periodic2_images( &
        mesh_fmm, r, 0.0d0, sim%box_min, sim%box_max, [1_i32, 2_i32], ref_layers, e_direct &
        )
      call solver%eval_e(mesh_fmm, r, e_fmm)

      norm_direct = sqrt(sum(e_direct*e_direct))
      if (norm_direct <= 1.0d-12) cycle
      norm_diff = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))
      rel_err = norm_diff/norm_direct
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    call assert_true(valid_count == 6_i32, 'core periodic2 accuracy test lost valid samples')
    call assert_true(max_rel_err <= 2.0d-1, 'core periodic2 E relative error exceeds 2e-1')
  end subroutine test_fmm_core_periodic2_field_accuracy

  subroutine test_fmm_periodic2_fallback_accuracy()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i, valid_count, ref_layers
    real(dp) :: queries(3, 6)
    real(dp) :: r(3), e_direct(3), e_fmm(3), max_rel_err
    real(dp) :: norm_direct, norm_diff, rel_err

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=24_i32, n_lat=12_i32, center=[0.5d0, 0.5d0, 0.0d0])
    call assign_periodic_test_dipole_charges(mesh_fmm, 1.0d-12)

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'periodic2'
    sim%field_periodic_image_layers = 1_i32
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call solver%init(mesh_fmm, sim)
    call solver%refresh(mesh_fmm)
    ref_layers = sim%field_periodic_image_layers
    if (trim(solver%periodic_far_correction) == 'm2l_root_trunc') then
      ref_layers = ref_layers + solver%periodic_ewald_layers
    end if

    queries(:, 1) = [0.15d0, 0.15d0, 1.10d0]
    queries(:, 2) = [0.85d0, 0.20d0, 1.25d0]
    queries(:, 3) = [0.20d0, 0.80d0, 1.40d0]
    queries(:, 4) = [0.75d0, 0.75d0, -1.15d0]
    queries(:, 5) = [0.55d0, 0.35d0, -1.30d0]
    queries(:, 6) = [0.35d0, 0.60d0, 1.55d0]

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      r = queries(:, i)
      call electric_field_at_periodic2_images( &
        mesh_fmm, r, sim%softening, sim%box_min, sim%box_max, [1_i32, 2_i32], ref_layers, e_direct &
        )
      call solver%eval_e(mesh_fmm, r, e_fmm)

      norm_direct = sqrt(sum(e_direct*e_direct))
      if (norm_direct <= 1.0d-12) cycle
      norm_diff = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))
      rel_err = norm_diff/norm_direct
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    write (*, '(A,I0,A,ES12.5)') &
      'test_fmm_periodic2_fallback_accuracy: valid_count=', valid_count, ', max_rel_err=', max_rel_err
    call assert_true(valid_count == 6_i32, 'fmm periodic2 fallback test lost valid samples')
    call assert_true(max_rel_err <= 2.0d-1, 'fmm periodic2 fallback E relative error exceeds 2e-1')
  end subroutine test_fmm_periodic2_fallback_accuracy

  subroutine test_fmm_periodic2_image_layers_accuracy()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i, valid_count, ref_layers
    real(dp) :: queries(3, 6)
    real(dp) :: r(3), e_direct(3), e_fmm(3), max_rel_err
    real(dp) :: norm_direct, norm_diff, rel_err

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=24_i32, n_lat=12_i32, center=[0.5d0, 0.5d0, 0.0d0])
    call assign_periodic_test_dipole_charges(mesh_fmm, 1.0d-12)

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'periodic2'
    sim%field_periodic_image_layers = 2_i32
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call solver%init(mesh_fmm, sim)
    call solver%refresh(mesh_fmm)
    ref_layers = sim%field_periodic_image_layers
    if (trim(solver%periodic_far_correction) == 'm2l_root_trunc') then
      ref_layers = ref_layers + solver%periodic_ewald_layers
    end if

    queries(:, 1) = [0.15d0, 0.15d0, -0.60d0]
    queries(:, 2) = [0.85d0, 0.20d0, -0.20d0]
    queries(:, 3) = [0.20d0, 0.80d0, 0.10d0]
    queries(:, 4) = [0.75d0, 0.75d0, 0.50d0]
    queries(:, 5) = [0.55d0, 0.35d0, -0.75d0]
    queries(:, 6) = [0.35d0, 0.60d0, 0.85d0]

    max_rel_err = 0.0d0
    valid_count = 0_i32
    do i = 1_i32, int(size(queries, 2), i32)
      r = queries(:, i)
      call electric_field_at_periodic2_images( &
        mesh_fmm, r, sim%softening, sim%box_min, sim%box_max, [1_i32, 2_i32], ref_layers, e_direct &
        )
      call solver%eval_e(mesh_fmm, r, e_fmm)

      norm_direct = sqrt(sum(e_direct*e_direct))
      if (norm_direct <= 1.0d-12) cycle
      norm_diff = sqrt(sum((e_fmm - e_direct)*(e_fmm - e_direct)))
      rel_err = norm_diff/norm_direct
      max_rel_err = max(max_rel_err, rel_err)
      valid_count = valid_count + 1_i32
    end do

    write (*, '(A,I0,A,ES12.5)') &
      'test_fmm_periodic2_image_layers_accuracy: valid_count=', valid_count, ', max_rel_err=', max_rel_err
    call assert_true(valid_count == 6_i32, 'fmm periodic2 image-layer test lost valid samples')
    call assert_true(max_rel_err <= 2.0d-1, 'fmm periodic2 image-layer E relative error exceeds 2e-1')
  end subroutine test_fmm_periodic2_image_layers_accuracy

  subroutine test_fmm_periodic2_m2l_cache_reuse()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver = field_solver_type()
    type(sim_config) :: sim

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=24_i32, n_lat=12_i32, center=[0.5d0, 0.5d0, 0.0d0])
    mesh_fmm%q_elem = 1.0d-12

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'periodic2'
    sim%field_periodic_image_layers = 1_i32
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call solver%init(mesh_fmm, sim)

    mesh_fmm%q_elem = 2.0d0*mesh_fmm%q_elem
    call solver%refresh(mesh_fmm)
  end subroutine test_fmm_periodic2_m2l_cache_reuse

  subroutine test_fmm_periodic2_default_m2l_root_oracle_mode()
    type(mesh_type) :: mesh_fmm
    type(field_solver_type) :: solver_default = field_solver_type()
    type(field_solver_type) :: solver_root = field_solver_type()
    type(sim_config) :: sim
    integer(i32) :: i
    real(dp) :: r(3), e_default(3), e_root(3)
    real(dp) :: max_delta_default_root

    call make_sphere(mesh_fmm, radius=0.2d0, n_lon=24_i32, n_lat=12_i32, center=[0.5d0, 0.5d0, 0.0d0])
    call assign_periodic_test_dipole_charges(mesh_fmm, 1.0d-12)

    sim = sim_config()
    sim%softening = 1.0d-4
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'periodic2'
    sim%field_periodic_image_layers = 1_i32
    sim%tree_min_nelem = 64_i32
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call solver_default%init(mesh_fmm, sim)
    call solver_default%refresh(mesh_fmm)

    sim%field_periodic_far_correction = 'm2l_root_oracle'
    sim%field_periodic_ewald_layers = 4_i32
    call solver_root%init(mesh_fmm, sim)
    call solver_root%refresh(mesh_fmm)

    call assert_true( &
      trim(solver_default%periodic_far_correction) == 'm2l_root_oracle', &
      'periodic2 default should normalize to m2l_root_oracle' &
      )
    call assert_true( &
      trim(solver_root%periodic_far_correction) == 'm2l_root_oracle', &
      'explicit periodic2 m2l_root_oracle should be preserved' &
      )
    call assert_true( &
      trim(solver_default%fmm_core_options%periodic_far_correction) == 'm2l_root_oracle', &
      'normalized periodic2 far correction should reach FMM core options' &
      )

    max_delta_default_root = 0.0d0
    do i = 1_i32, 4_i32
      r = [0.15d0 + 0.2d0*real(i - 1_i32, dp), 0.20d0 + 0.15d0*real(i - 1_i32, dp), -0.6d0 + 0.35d0*real(i - 1_i32, dp)]
      call solver_default%eval_e(mesh_fmm, r, e_default)
      call solver_root%eval_e(mesh_fmm, r, e_root)
      max_delta_default_root = max(max_delta_default_root, sqrt(sum((e_default - e_root)*(e_default - e_root))))
    end do

    call assert_true( &
      max_delta_default_root <= 1.0d-18, &
      'default periodic2 and explicit m2l_root_oracle should agree' &
      )
  end subroutine test_fmm_periodic2_default_m2l_root_oracle_mode

  subroutine assign_periodic_test_dipole_charges(mesh, scale)
    type(mesh_type), intent(inout) :: mesh
    real(dp), intent(in) :: scale
    integer(i32) :: elem_idx, npos, nneg
    real(dp) :: neg_scale

    npos = 0_i32
    nneg = 0_i32
    do elem_idx = 1_i32, mesh%nelem
      if (mesh%center_z(elem_idx) >= 0.0d0) then
        npos = npos + 1_i32
      else
        nneg = nneg + 1_i32
      end if
    end do

    if (npos <= 0_i32 .or. nneg <= 0_i32) then
      do elem_idx = 1_i32, mesh%nelem
        if (mod(elem_idx, 2_i32) == 0_i32) then
          mesh%q_elem(elem_idx) = -scale
        else
          mesh%q_elem(elem_idx) = scale
        end if
      end do
      return
    end if

    neg_scale = scale*real(npos, dp)/real(nneg, dp)
    do elem_idx = 1_i32, mesh%nelem
      if (mesh%center_z(elem_idx) >= 0.0d0) then
        mesh%q_elem(elem_idx) = scale
      else
        mesh%q_elem(elem_idx) = -neg_scale
      end if
    end do
  end subroutine assign_periodic_test_dipole_charges

  subroutine electric_field_at_periodic2_images(mesh, r, softening, box_min, box_max, periodic_axes, nimg, e)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: softening
    real(dp), intent(in) :: box_min(3), box_max(3)
    integer(i32), intent(in) :: periodic_axes(2)
    integer(i32), intent(in) :: nimg
    real(dp), intent(out) :: e(3)

    integer(i32) :: i, ix, iy
    integer(i32) :: axis1, axis2
    real(dp) :: l1, l2, soft2, shift1, shift2
    real(dp) :: src(3), d(3), r2, inv_r3
    real(dp) :: ex, ey, ez

    axis1 = periodic_axes(1)
    axis2 = periodic_axes(2)
    l1 = box_max(axis1) - box_min(axis1)
    l2 = box_max(axis2) - box_min(axis2)
    soft2 = softening*softening
    ex = 0.0d0
    ey = 0.0d0
    ez = 0.0d0

    do i = 1_i32, mesh%nelem
      do ix = -nimg, nimg
        do iy = -nimg, nimg
          src = [mesh%center_x(i), mesh%center_y(i), mesh%center_z(i)]
          shift1 = real(ix, dp)*l1
          shift2 = real(iy, dp)*l2
          src(axis1) = src(axis1) + shift1
          src(axis2) = src(axis2) + shift2
          d = r - src
          r2 = sum(d*d) + soft2
          inv_r3 = 1.0d0/(sqrt(r2)*r2)
          ex = ex + mesh%q_elem(i)*inv_r3*d(1)
          ey = ey + mesh%q_elem(i)*inv_r3*d(2)
          ez = ez + mesh%q_elem(i)*inv_r3*d(3)
        end do
      end do
    end do

    e(1) = k_coulomb*ex
    e(2) = k_coulomb*ey
    e(3) = k_coulomb*ez
  end subroutine electric_field_at_periodic2_images
end program test_dynamics

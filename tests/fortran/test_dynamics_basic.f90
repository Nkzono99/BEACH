!> 基本物理テスト: 電場評価・Boris更新・衝突判定・periodic2衝突の基礎検証。
program test_dynamics_basic
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_types, only: mesh_type, hit_info, sim_config, bc_open, bc_periodic
  use bem_mesh, only: init_mesh, prepare_periodic2_collision_mesh
  use bem_templates, only: make_plane
  use bem_field, only: electric_field_at
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
  call test_periodic2_collision_wrap()
  call test_periodic2_collision_multi_cell()
  call test_periodic2_collision_canonical_prepare()
  call test_periodic2_collision_grid_equivalence()

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

  subroutine test_periodic2_collision_wrap()
    type(mesh_type) :: mesh_periodic
    type(sim_config) :: sim
    type(hit_info) :: hit_periodic
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [0.1d0, 0.2d0, 0.0d0]
    tri_v1(:, 1) = [0.3d0, 0.2d0, 0.0d0]
    tri_v2(:, 1) = [0.1d0, 0.4d0, 0.0d0]
    call init_mesh(mesh_periodic, tri_v0, tri_v1, tri_v2)
    call init_periodic2_test_sim(sim)
    call prepare_periodic2_collision_mesh(mesh_periodic, sim)

    call find_first_hit(mesh_periodic, [1.2d0, 0.25d0, 1.0d0], [1.2d0, 0.25d0, -1.0d0], hit_periodic, sim=sim)
    call assert_true(hit_periodic%has_hit, 'periodic2 collision should hit opposite-side image')
    call assert_equal_i32(hit_periodic%elem_idx, 1_i32, 'periodic2 collision elem_idx mismatch')
    call assert_close_dp(hit_periodic%t, 0.5d0, 1.0d-12, 'periodic2 collision t mismatch')
    call assert_close_dp(hit_periodic%pos(1), 1.2d0, 1.0d-12, 'periodic2 collision unwrapped x mismatch')
    call assert_close_dp(hit_periodic%pos_wrapped(1), 0.2d0, 1.0d-12, 'periodic2 collision wrapped x mismatch')
    call assert_close_dp(hit_periodic%pos_wrapped(2), 0.25d0, 1.0d-12, 'periodic2 collision wrapped y mismatch')
    call assert_equal_i32(hit_periodic%image_shift(1), 1_i32, 'periodic2 collision x image shift mismatch')
    call assert_equal_i32(hit_periodic%image_shift(2), 0_i32, 'periodic2 collision y image shift mismatch')
  end subroutine test_periodic2_collision_wrap

  subroutine test_periodic2_collision_multi_cell()
    type(mesh_type) :: mesh_periodic
    type(sim_config) :: sim
    type(hit_info) :: hit_periodic
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [0.1d0, 0.2d0, 0.0d0]
    tri_v1(:, 1) = [0.3d0, 0.2d0, 0.0d0]
    tri_v2(:, 1) = [0.1d0, 0.4d0, 0.0d0]
    call init_mesh(mesh_periodic, tri_v0, tri_v1, tri_v2)
    call init_periodic2_test_sim(sim)
    call prepare_periodic2_collision_mesh(mesh_periodic, sim)

    call find_first_hit(mesh_periodic, [3.2d0, -0.75d0, 1.0d0], [3.2d0, -0.75d0, -1.0d0], hit_periodic, sim=sim)
    call assert_true(hit_periodic%has_hit, 'periodic2 collision should support multi-cell image search')
    call assert_close_dp(hit_periodic%pos_wrapped(1), 0.2d0, 1.0d-12, 'periodic2 multi-cell wrapped x mismatch')
    call assert_close_dp(hit_periodic%pos_wrapped(2), 0.25d0, 1.0d-12, 'periodic2 multi-cell wrapped y mismatch')
    call assert_equal_i32(hit_periodic%image_shift(1), 3_i32, 'periodic2 multi-cell x image shift mismatch')
    call assert_equal_i32(hit_periodic%image_shift(2), -1_i32, 'periodic2 multi-cell y image shift mismatch')
  end subroutine test_periodic2_collision_multi_cell

  subroutine test_periodic2_collision_canonical_prepare()
    type(mesh_type) :: mesh_periodic
    type(sim_config) :: sim
    type(hit_info) :: hit_periodic
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [0.1d0, 0.2d0, 0.0d0]
    tri_v1(:, 1) = [0.2d0, 0.2d0, 0.0d0]
    tri_v2(:, 1) = [-0.1d0, 0.4d0, 0.0d0]
    call init_mesh(mesh_periodic, tri_v0, tri_v1, tri_v2)
    call init_periodic2_test_sim(sim)
    call prepare_periodic2_collision_mesh(mesh_periodic, sim)

    call assert_true(mesh_periodic%periodic2_collision_ready, 'periodic2 collision prep flag should be set')
    call assert_close_dp(mesh_periodic%bb_min(1, 1), -0.1d0, 1.0d-12, 'canonical triangle bb_min mismatch')
    call assert_close_dp(mesh_periodic%bb_max(1, 1), 0.2d0, 1.0d-12, 'canonical triangle bb_max mismatch')

    call find_first_hit(mesh_periodic, [1.1d0, 0.25d0, 1.0d0], [1.1d0, 0.25d0, -1.0d0], hit_periodic, sim=sim)
    call assert_true(hit_periodic%has_hit, 'canonical periodic triangle should remain hittable')
    call assert_close_dp(hit_periodic%pos_wrapped(1), 0.1d0, 1.0d-12, 'canonical periodic triangle wrapped x mismatch')
  end subroutine test_periodic2_collision_canonical_prepare

  subroutine test_periodic2_collision_grid_equivalence()
    type(mesh_type) :: mesh_grid, mesh_linear
    type(sim_config) :: sim

    call make_plane(mesh_grid, size_x=1.0d0, size_y=1.0d0, nx=8_i32, ny=8_i32, center=[0.5d0, 0.5d0, 0.0d0])
    call init_periodic2_test_sim(sim)
    call prepare_periodic2_collision_mesh(mesh_grid, sim)
    call assert_true(mesh_grid%use_collision_grid, 'periodic2 dense mesh should use collision grid')

    mesh_linear = mesh_grid
    mesh_linear%use_collision_grid = .false.
    call assert_periodic_hit_equivalent( &
      mesh_grid, mesh_linear, sim, [2.2d0, -1.7d0, 1.0d0], [2.2d0, -1.7d0, -1.0d0], &
      'periodic2 grid/linear equivalence' &
      )
  end subroutine test_periodic2_collision_grid_equivalence

  subroutine assert_periodic_hit_equivalent(mesh_fast, mesh_ref, sim, p0, p1, label)
    type(mesh_type), intent(in) :: mesh_fast, mesh_ref
    type(sim_config), intent(in) :: sim
    real(dp), intent(in) :: p0(3), p1(3)
    character(len=*), intent(in) :: label
    type(hit_info) :: hit_fast, hit_ref

    call find_first_hit(mesh_fast, p0, p1, hit_fast, sim=sim)
    call find_first_hit(mesh_ref, p0, p1, hit_ref, sim=sim)

    call assert_true(hit_fast%has_hit .eqv. hit_ref%has_hit, trim(label)//': has_hit mismatch')
    if (.not. hit_ref%has_hit) return

    call assert_equal_i32(hit_fast%elem_idx, hit_ref%elem_idx, trim(label)//': elem_idx mismatch')
    call assert_close_dp(hit_fast%t, hit_ref%t, 1.0d-11, trim(label)//': t mismatch')
    call assert_allclose_1d(hit_fast%pos, hit_ref%pos, 1.0d-10, trim(label)//': position mismatch')
    call assert_allclose_1d(hit_fast%pos_wrapped, hit_ref%pos_wrapped, 1.0d-10, trim(label)//': wrapped mismatch')
    call assert_equal_i32(hit_fast%image_shift(1), hit_ref%image_shift(1), trim(label)//': shift(1) mismatch')
    call assert_equal_i32(hit_fast%image_shift(2), hit_ref%image_shift(2), trim(label)//': shift(2) mismatch')
  end subroutine assert_periodic_hit_equivalent

  subroutine init_periodic2_test_sim(sim)
    type(sim_config), intent(out) :: sim

    sim = sim_config()
    sim%field_solver = 'fmm'
    sim%field_bc_mode = 'periodic2'
    sim%use_box = .true.
    sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim%bc_high = [bc_periodic, bc_periodic, bc_open]
  end subroutine init_periodic2_test_sim

end program test_dynamics_basic

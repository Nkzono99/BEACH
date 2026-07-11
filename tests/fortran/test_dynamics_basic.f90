!> 基本物理テスト: 電場評価・Boris更新・衝突判定・periodic2衝突の基礎検証。
program test_dynamics_basic
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: k_coulomb
  use bem_types, only: mesh_type, hit_info, sim_config, bc_open, bc_periodic
  use bem_mesh, only: init_mesh, prepare_periodic2_collision_mesh
  use bem_templates, only: make_plane
  use bem_field, only: electric_field_at
  use bem_pusher, only: boris_push, boris_update_velocity
  use bem_collision, only: collision_query_image_limit, collision_query_index_range, collision_query_ok, &
                           compute_periodic_shift_bounds, find_first_hit, segment_triangle_intersect
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, assert_allclose_1d, delete_file_if_exists
  implicit none

  type(mesh_type) :: mesh_field, mesh_hit
  type(hit_info) :: hit
  real(dp) :: v0_field(3, 1), v1_field(3, 1), v2_field(3, 1), q0_field(1)
  real(dp) :: v0_hit(3, 2), v1_hit(3, 2), v2_hit(3, 2)
  real(dp) :: e(3), x_new(3), v_new(3), speed0, speed1
  real(dp) :: inv_r3, expected_ex
  character(len=*), parameter :: collision_failure_path = 'test_dynamics_basic_collision_failure_tmp.log'
  character(len=64) :: run_mode

  call get_command_argument(1, run_mode)
  if (trim(run_mode) == '--collision-status-omitted-probe') then
    call run_collision_status_omitted_probe()
    error stop 'collision status omitted probe unexpectedly completed'
  end if

  call test_init(21)

  call test_begin('electric_field_at')
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
  call test_end()

  call test_begin('boris_push_e_only')
  call boris_push( &
    [0.0d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 0.0d0], &
    2.0d0, 1.0d0, 0.5d0, [1.0d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 0.0d0], &
    x_new, v_new &
    )
  call assert_allclose_1d(v_new, [1.0d0, 0.0d0, 0.0d0], 1.0d-12, 'boris (E only) velocity mismatch')
  call assert_allclose_1d(x_new, [0.25d0, 0.0d0, 0.0d0], 1.0d-12, 'boris (E only) position mismatch')
  call test_end()

  call test_begin('boris_push_b_speed_preserve')
  call boris_push( &
    [0.0d0, 0.0d0, 0.0d0], [1.0d0, 2.0d0, -0.5d0], &
    1.0d0, 1.0d0, 0.1d0, [0.0d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 2.0d0], &
    x_new, v_new &
    )
  speed0 = sqrt(1.0d0 + 4.0d0 + 0.25d0)
  speed1 = sqrt(sum(v_new*v_new))
  call assert_close_dp(speed1, speed0, 1.0d-12, 'boris should preserve speed when E=0')
  call test_end()

  call test_begin('boris_push_pure_b_second_order')
  call test_boris_push_pure_b_second_order()
  call test_end()

  call test_begin('boris_update_velocity_e_only')
  call test_boris_update_velocity_e_only()
  call test_end()

  call test_begin('boris_update_velocity_pure_b_speed')
  call test_boris_update_velocity_pure_b_speed()
  call test_end()

  call test_begin('boris_update_velocity_time_reversal')
  call test_boris_update_velocity_time_reversal()
  call test_end()

  call test_begin('boris_push_velocity_helper_equivalence')
  call test_boris_push_velocity_helper_equivalence()
  call test_end()

  call test_begin('find_first_hit')
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
  call test_end()

  call test_begin('segment_triangle_intersect_small_scale')
  call test_segment_triangle_intersect_small_scale()
  call test_end()

  call test_begin('collision_grid_equivalence')
  call test_collision_grid_equivalence()
  call test_end()

  call test_begin('periodic2_collision_wrap')
  call test_periodic2_collision_wrap()
  call test_end()

  call test_begin('periodic2_collision_multi_cell')
  call test_periodic2_collision_multi_cell()
  call test_end()

  call test_begin('periodic2_collision_i32_upper_boundary')
  call test_periodic2_collision_i32_upper_boundary()
  call test_end()

  call test_begin('periodic2_collision_i32_lower_boundary')
  call test_periodic2_collision_i32_lower_boundary()
  call test_end()

  call test_begin('periodic_shift_bounds_legacy_signature')
  call test_periodic_shift_bounds_legacy_signature()
  call test_end()

  call test_begin('periodic2_collision_index_range_guard')
  call test_periodic2_collision_index_range_guard()
  call test_end()

  call test_begin('periodic2_collision_runaway_segment_guard')
  call test_periodic2_collision_runaway_segment_guard()
  call test_end()

  call test_begin('periodic2_collision_status_omitted_fails_closed')
  call test_periodic2_collision_status_omitted_fails_closed()
  call test_end()

  call test_begin('periodic2_collision_canonical_prepare')
  call test_periodic2_collision_canonical_prepare()
  call test_end()

  call test_begin('periodic2_collision_grid_equivalence')
  call test_periodic2_collision_grid_equivalence()
  call test_end()

  call test_summary()

contains

  subroutine test_boris_push_pure_b_second_order()
    integer(i32), parameter :: nstep_values(3) = [4_i32, 8_i32, 16_i32]
    real(dp), parameter :: final_time = 1.0d0, omega = 1.0d0
    integer(i32) :: refinement, step
    real(dp) :: dt, position(3), velocity(3), position_new(3), velocity_new(3)
    real(dp) :: position_exact(3), velocity_exact(3)
    real(dp) :: position_error(3), velocity_error(3)
    real(dp) :: position_ratio, velocity_ratio

    position_exact = [sin(omega*final_time)/omega, (cos(omega*final_time) - 1.0d0)/omega, 0.0d0]
    velocity_exact = [cos(omega*final_time), -sin(omega*final_time), 0.0d0]

    do refinement = 1_i32, 3_i32
      dt = final_time/real(nstep_values(refinement), dp)
      position = 0.0d0
      velocity = [1.0d0, 0.0d0, 0.0d0]
      do step = 1_i32, nstep_values(refinement)
        call boris_push( &
          position, velocity, 1.0d0, 1.0d0, dt, &
          [0.0d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 1.0d0], &
          position_new, velocity_new &
          )
        position = position_new
        velocity = velocity_new
      end do
      position_error(refinement) = sqrt(sum((position - position_exact)**2))
      velocity_error(refinement) = sqrt(sum((velocity - velocity_exact)**2))
    end do

    call assert_true(all(ieee_is_finite(position_error)), 'pure-B position errors must be finite')
    call assert_true(all(ieee_is_finite(velocity_error)), 'pure-B velocity errors must be finite')
    call assert_true(all(position_error > 0.0d0), 'pure-B position errors must be positive')
    call assert_true(all(velocity_error > 0.0d0), 'pure-B velocity errors must be positive')
    do refinement = 1_i32, 2_i32
      position_ratio = position_error(refinement)/position_error(refinement + 1_i32)
      velocity_ratio = velocity_error(refinement)/velocity_error(refinement + 1_i32)
      call assert_true( &
        position_ratio >= 3.2d0 .and. position_ratio <= 4.8d0, &
        'pure-B position error ratio must show second-order convergence' &
        )
      call assert_true( &
        velocity_ratio >= 3.2d0 .and. velocity_ratio <= 4.8d0, &
        'pure-B velocity error ratio must show second-order convergence' &
        )
    end do
    write (*, '(a,3(es16.8,a))') &
      'BEACH_CONVERGENCE,boris_dt,0.25_0.125_0.0625,', position_error(1), ',', position_error(2), ',', &
      position_error(3), ',second_order'
  end subroutine test_boris_push_pure_b_second_order

  subroutine test_boris_update_velocity_e_only()
    real(dp) :: velocity_new(3)

    call boris_update_velocity( &
      [0.0d0, 0.0d0, 0.0d0], 2.0d0, 1.0d0, 0.5d0, &
      [1.0d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 0.0d0], velocity_new &
      )

    call assert_allclose_1d( &
      velocity_new, [1.0d0, 0.0d0, 0.0d0], 0.0d0, &
      'boris velocity update (E only) mismatch' &
      )
  end subroutine test_boris_update_velocity_e_only

  subroutine test_boris_update_velocity_pure_b_speed()
    real(dp) :: velocity(3), velocity_new(3)

    velocity = [1.0d0, 0.0d0, 0.0d0]
    call boris_update_velocity( &
      velocity, 1.0d0, 1.0d0, 1.0d0, &
      [0.0d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 2.0d0], velocity_new &
      )

    call assert_allclose_1d( &
      velocity_new, [0.0d0, -1.0d0, 0.0d0], 0.0d0, &
      'boris velocity update (pure B) mismatch' &
      )
    call assert_close_dp( &
      sum(velocity_new*velocity_new), sum(velocity*velocity), 1.0d-12, &
      'boris velocity update should preserve squared speed when E=0' &
      )
  end subroutine test_boris_update_velocity_pure_b_speed

  subroutine test_boris_update_velocity_time_reversal()
    real(dp) :: velocity(3), velocity_forward(3), velocity_backward(3)

    velocity = [1.0d0, 0.0d0, 0.0d0]
    call boris_update_velocity( &
      velocity, 1.0d0, 1.0d0, 1.0d0, &
      [0.5d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 2.0d0], velocity_forward &
      )
    call assert_allclose_1d( &
      velocity_forward, [0.25d0, -1.25d0, 0.0d0], 0.0d0, &
      'forward Boris velocity update mismatch' &
      )

    call boris_update_velocity( &
      velocity_forward, 1.0d0, 1.0d0, -1.0d0, &
      [0.5d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 2.0d0], velocity_backward &
      )
    call assert_allclose_1d( &
      velocity_backward, velocity, 0.0d0, &
      'backward Boris velocity update should recover the initial velocity' &
      )
  end subroutine test_boris_update_velocity_time_reversal

  subroutine test_boris_push_velocity_helper_equivalence()
    integer(i64) :: push_bits(3), update_bits(3), expected_bits(3)
    real(dp) :: x(3), velocity(3), electric_field(3), magnetic_field(3)
    real(dp) :: x_push(3), velocity_push(3), velocity_update(3)

    x = [0.125_dp, -0.25_dp, 0.5_dp]
    velocity = [0.75_dp, -1.25_dp, 0.375_dp]
    electric_field = [0.5_dp, -0.25_dp, 0.75_dp]
    magnetic_field = [-0.375_dp, 0.625_dp, 0.25_dp]

    call boris_push( &
      x, velocity, -1.5_dp, 2.0_dp, 0.125_dp, electric_field, magnetic_field, &
      x_push, velocity_push &
      )
    call boris_update_velocity( &
      velocity, -1.5_dp, 2.0_dp, 0.125_dp, electric_field, magnetic_field, &
      velocity_update &
      )
    push_bits = transfer(velocity_push, push_bits)
    update_bits = transfer(velocity_update, update_bits)
    call assert_true(all(push_bits == update_bits), 'boris_push/helper velocity bits should match')

    call boris_push( &
      x, velocity, -1.5_dp, 2.0_dp, 0.0_dp, electric_field, magnetic_field, &
      x_push, velocity_push &
      )
    call boris_update_velocity( &
      velocity, -1.5_dp, 2.0_dp, 0.0_dp, electric_field, magnetic_field, &
      velocity_update &
      )
    expected_bits = transfer(velocity, expected_bits)
    push_bits = transfer(velocity_push, push_bits)
    update_bits = transfer(velocity_update, update_bits)
    call assert_true(all(push_bits == expected_bits), 'boris_push velocity should be identity when dt=0')
    call assert_true(all(update_bits == expected_bits), 'boris velocity helper should be identity when dt=0')
    expected_bits = transfer(x, expected_bits)
    push_bits = transfer(x_push, push_bits)
    call assert_true(all(push_bits == expected_bits), 'boris_push position should be identity when dt=0')
  end subroutine test_boris_push_velocity_helper_equivalence

  subroutine test_segment_triangle_intersect_small_scale()
    logical :: ok
    real(dp) :: p0(3), p1(3), tri0(3), tri1(3), tri2(3), t, h(3)
    real(dp), parameter :: scale = 9.899494936611664d-5

    tri0 = [0.0d0, 0.0d0, 0.0d0]*scale
    tri1 = [1.0d0, 0.0d0, 0.0d0]*scale
    tri2 = [0.0d0, 1.0d0, 0.0d0]*scale
    p0 = [0.25d0, 0.25d0, 0.02d0]*scale
    p1 = [0.25d0, 0.25d0, -0.02d0]*scale

    call segment_triangle_intersect(p0, p1, tri0, tri1, tri2, ok, t, h)

    call assert_true(ok, 'small-scale segment should hit triangle')
    call assert_close_dp(t, 0.5d0, 1.0d-12, 'small-scale hit parameter mismatch')
    call assert_allclose_1d(h, [0.25d0, 0.25d0, 0.0d0]*scale, 1.0d-16, 'small-scale hit position mismatch')
  end subroutine test_segment_triangle_intersect_small_scale

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

  subroutine test_periodic2_collision_runaway_segment_guard()
    type(mesh_type) :: mesh_periodic
    type(sim_config) :: sim
    type(hit_info) :: hit_periodic
    integer(i32) :: query_status
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [0.1d0, 0.2d0, 0.0d0]
    tri_v1(:, 1) = [0.3d0, 0.2d0, 0.0d0]
    tri_v2(:, 1) = [0.1d0, 0.4d0, 0.0d0]
    call init_mesh(mesh_periodic, tri_v0, tri_v1, tri_v2)
    call init_periodic2_test_sim(sim)
    call prepare_periodic2_collision_mesh(mesh_periodic, sim)

    call find_first_hit( &
      mesh_periodic, [0.2d0, 0.25d0, 1.0d0], [5000.2d0, 0.25d0, -1.0d0], hit_periodic, &
      sim=sim, status=query_status &
      )
    call assert_equal_i32( &
      query_status, collision_query_image_limit, &
      'periodic2 runaway segment should report the image enumeration limit' &
      )
  end subroutine test_periodic2_collision_runaway_segment_guard

  subroutine test_periodic2_collision_i32_upper_boundary()
    type(mesh_type) :: mesh_periodic
    type(sim_config) :: sim
    type(hit_info) :: hit_periodic
    integer(i32) :: query_status
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1), boundary_x

    tri_v0(:, 1) = [0.1d0, 0.2d0, 0.0d0]
    tri_v1(:, 1) = [0.3d0, 0.2d0, 0.0d0]
    tri_v2(:, 1) = [0.1d0, 0.4d0, 0.0d0]
    call init_mesh(mesh_periodic, tri_v0, tri_v1, tri_v2)
    call init_periodic2_test_sim(sim)
    call prepare_periodic2_collision_mesh(mesh_periodic, sim)

    boundary_x = real(huge(0_i32), dp) + 0.2d0
    call find_first_hit( &
      mesh_periodic, [boundary_x, 0.25d0, 1.0d0], [boundary_x, 0.25d0, -1.0d0], &
      hit_periodic, sim=sim, status=query_status &
      )
    call assert_equal_i32(query_status, collision_query_ok, 'i32 upper-bound collision status mismatch')
    call assert_true(hit_periodic%has_hit, 'i32 upper-bound periodic image should remain hittable')
    call assert_equal_i32( &
      hit_periodic%image_shift(1), huge(0_i32), 'i32 upper-bound image shift mismatch' &
      )
  end subroutine test_periodic2_collision_i32_upper_boundary

  subroutine test_periodic2_collision_i32_lower_boundary()
    type(mesh_type) :: mesh_periodic
    type(sim_config) :: sim
    type(hit_info) :: hit_periodic
    integer(i32) :: query_status
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1), boundary_x

    tri_v0(:, 1) = [0.1d0, 0.2d0, 0.0d0]
    tri_v1(:, 1) = [0.3d0, 0.2d0, 0.0d0]
    tri_v2(:, 1) = [0.1d0, 0.4d0, 0.0d0]
    call init_mesh(mesh_periodic, tri_v0, tri_v1, tri_v2)
    call init_periodic2_test_sim(sim)
    call prepare_periodic2_collision_mesh(mesh_periodic, sim)

    boundary_x = -real(huge(0_i32), dp) - 1.0d0 + 0.2d0
    call find_first_hit( &
      mesh_periodic, [boundary_x, 0.25d0, 1.0d0], [boundary_x, 0.25d0, -1.0d0], &
      hit_periodic, sim=sim, status=query_status &
      )
    call assert_equal_i32(query_status, collision_query_ok, 'i32 lower-bound collision status mismatch')
    call assert_true(hit_periodic%has_hit, 'i32 lower-bound periodic image should remain hittable')
    call assert_equal_i32( &
      hit_periodic%image_shift(1), -huge(0_i32) - 1_i32, 'i32 lower-bound image shift mismatch' &
      )
  end subroutine test_periodic2_collision_i32_lower_boundary

  subroutine test_periodic_shift_bounds_legacy_signature()
    type(mesh_type) :: mesh_periodic
    integer(i32) :: nmin, nmax
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [0.1d0, 0.2d0, 0.0d0]
    tri_v1(:, 1) = [0.3d0, 0.2d0, 0.0d0]
    tri_v2(:, 1) = [0.1d0, 0.4d0, 0.0d0]
    call init_mesh(mesh_periodic, tri_v0, tri_v1, tri_v2)

    call compute_periodic_shift_bounds( &
      mesh_periodic, [0.2d0, 0.25d0, 1.0d0], [0.2d0, 0.25d0, -1.0d0], 1_i32, 1.0d0, nmin, nmax &
      )
    call assert_equal_i32(nmin, 0_i32, 'legacy shift-bound nmin mismatch')
    call assert_equal_i32(nmax, 0_i32, 'legacy shift-bound nmax mismatch')
  end subroutine test_periodic_shift_bounds_legacy_signature

  subroutine test_periodic2_collision_index_range_guard()
    type(mesh_type) :: mesh_periodic
    type(sim_config) :: sim
    type(hit_info) :: hit_periodic
    integer(i32) :: query_status
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1), out_of_range_x

    tri_v0(:, 1) = [0.1d0, 0.2d0, 0.0d0]
    tri_v1(:, 1) = [0.3d0, 0.2d0, 0.0d0]
    tri_v2(:, 1) = [0.1d0, 0.4d0, 0.0d0]
    call init_mesh(mesh_periodic, tri_v0, tri_v1, tri_v2)
    call init_periodic2_test_sim(sim)
    call prepare_periodic2_collision_mesh(mesh_periodic, sim)

    out_of_range_x = real(huge(0_i32), dp) + 1024.0d0
    call find_first_hit( &
      mesh_periodic, [out_of_range_x, 0.25d0, 1.0d0], [out_of_range_x, 0.25d0, -1.0d0], &
      hit_periodic, sim=sim, status=query_status &
      )
    call assert_equal_i32( &
      query_status, collision_query_index_range, &
      'periodic2 collision should reject image indices outside the i32 range' &
      )
  end subroutine test_periodic2_collision_index_range_guard

  subroutine test_periodic2_collision_status_omitted_fails_closed()
    character(len=1024) :: executable_path, command, child_line
    integer :: child_exit_status, child_cmd_status, child_unit, child_ios
    logical :: saw_incomplete_message

    call get_command_argument(0, executable_path)
    call delete_file_if_exists(collision_failure_path)
    command = '"'//trim(executable_path)//'" --collision-status-omitted-probe > '//collision_failure_path//' 2>&1'
    call execute_command_line( &
      trim(command), wait=.true., exitstat=child_exit_status, cmdstat=child_cmd_status &
      )

    call assert_equal_i32(int(child_cmd_status, i32), 0_i32, 'collision status omitted probe command status mismatch')
    call assert_true(child_exit_status /= 0, 'collision status omitted probe should fail closed')

    saw_incomplete_message = .false.
    open (newunit=child_unit, file=collision_failure_path, status='old', action='read', iostat=child_ios)
    if (child_ios /= 0) error stop 'failed to read collision status omitted probe output'
    do
      read (child_unit, '(A)', iostat=child_ios) child_line
      if (child_ios /= 0) exit
      saw_incomplete_message = saw_incomplete_message .or. &
                               index(child_line, 'image enumeration limit exceeded') > 0
    end do
    close (child_unit)
    call delete_file_if_exists(collision_failure_path)

    call assert_true(saw_incomplete_message, 'status-omitted collision query should report incomplete work')
  end subroutine test_periodic2_collision_status_omitted_fails_closed

  subroutine run_collision_status_omitted_probe()
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

    call find_first_hit( &
      mesh_periodic, [0.2d0, 0.25d0, 1.0d0], [5000.2d0, 0.25d0, -1.0d0], hit_periodic, sim=sim &
      )
  end subroutine run_collision_status_omitted_probe

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

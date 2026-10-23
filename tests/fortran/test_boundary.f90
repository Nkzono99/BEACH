!> ボックス境界条件の各モードを検証するテスト。
program test_boundary
  use bem_kinds, only: dp, i32
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use bem_boundary, only: apply_ray_box_boundary, boundary_event_type, find_first_boundary_event, &
                          apply_escape_reflect_periodic_event
  use bem_status_codes, only: boundary_event_ok, boundary_event_invalid_geometry
  use bem_types, only: particle_domain_config, bc_open, bc_reflect, bc_periodic, bc_redistributed_reflect
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, assert_allclose_1d
  implicit none

  type(particle_domain_config) :: cfg
  real(dp) :: x(3), v(3), expected(3)
  logical :: alive

  call test_init(26)

  cfg = particle_domain_config()
  cfg%use_box = .true.
  cfg%box_min = [0.0d0, 0.0d0, 0.0d0]
  cfg%box_max = [1.0d0, 1.0d0, 1.0d0]

  ! photo raycast の ray 用 box 境界。粒子の境界処理（再配置反射を含む）は下の event API のテストが担う。
  call test_begin('ray_open_face_stops_ray')
  x = [1.1d0, 0.5d0, 0.5d0]
  v = [1.0d0, 2.0d0, 3.0d0]
  alive = .true.
  cfg%bc_high(1) = bc_open
  call apply_ray_box_boundary(cfg, x, v, alive)
  call assert_true(.not. alive, 'open face should stop the ray')
  call test_end()

  call test_begin('ray_reflect_face')
  x = [1.2d0, 0.5d0, 0.5d0]
  v = [1.0d0, 2.0d0, 3.0d0]
  alive = .true.
  cfg%bc_high(1) = bc_reflect
  call apply_ray_box_boundary(cfg, x, v, alive)
  expected = [0.8d0, 0.5d0, 0.5d0]
  call assert_true(alive, 'reflect face should keep the ray')
  call assert_allclose_1d(x, expected, 1.0d-10, 'ray reflect position mismatch')
  call assert_close_dp(v(1), -1.0d0, 1.0d-12, 'ray reflect direction mismatch')
  call assert_close_dp(v(2), 2.0d0, 1.0d-12, 'ray reflect should preserve the tangential direction')
  call test_end()

  call test_begin('ray_redistributed_reflect_face_stops_ray')
  ! ray は面内再配置の乱数を持たないので、再配置反射の面では止める。
  x = [1.2d0, 0.5d0, 0.5d0]
  v = [1.0d0, 2.0d0, 3.0d0]
  alive = .true.
  cfg%bc_high(1) = bc_redistributed_reflect
  call apply_ray_box_boundary(cfg, x, v, alive)
  call assert_true(.not. alive, 'redistributed reflect face should stop the ray')
  call test_end()

  call test_begin('ray_periodic_face')
  x = [-0.2d0, 0.5d0, 0.5d0]
  v = [1.0d0, 2.0d0, 3.0d0]
  alive = .true.
  cfg%bc_low(1) = bc_periodic
  call apply_ray_box_boundary(cfg, x, v, alive)
  expected = [0.8d0, 0.5d0, 0.5d0]
  call assert_true(alive, 'periodic face should keep the ray')
  call assert_allclose_1d(x, expected, 1.0d-10, 'ray periodic position mismatch')
  call assert_allclose_1d(v, [1.0d0, 2.0d0, 3.0d0], 1.0d-12, 'periodic face should preserve the direction')
  call test_end()

  call test_begin('ray_disabled_box')
  cfg = particle_domain_config()
  x = [1.2d0, 0.5d0, 0.5d0]
  v = [1.0d0, 2.0d0, 3.0d0]
  alive = .true.
  call apply_ray_box_boundary(cfg, x, v, alive)
  call assert_true(alive, 'disabled box should keep the ray')
  call assert_allclose_1d(x, [1.2d0, 0.5d0, 0.5d0], 1.0d-12, 'disabled box should preserve the position')
  call test_end()

  call test_begin('ray_degenerate_box')
  cfg = particle_domain_config()
  cfg%use_box = .true.
  cfg%box_min = [0.0d0, 0.0d0, 0.0d0]
  cfg%box_max = [0.0d0, 1.0d0, 1.0d0]
  x = [0.0d0, 0.5d0, 0.5d0]
  v = [1.0d0, 2.0d0, 3.0d0]
  alive = .true.
  call apply_ray_box_boundary(cfg, x, v, alive)
  call assert_true(.not. alive, 'degenerate box should stop the ray')
  call test_end()

  call test_begin('first_boundary_event_single_face')
  call test_first_boundary_event_single_face()
  call test_end()

  call test_begin('first_boundary_event_precedes_later_face')
  call test_first_boundary_event_precedes_later_face()
  call test_end()

  call test_begin('first_boundary_event_no_crossing')
  call test_first_boundary_event_no_crossing()
  call test_end()

  call test_begin('boundary_event_corner_reflect')
  call test_boundary_event_corner_reflect()
  call test_end()

  call test_begin('boundary_event_redistributed_reflect')
  call test_boundary_event_redistributed_reflect()
  call test_end()

  call test_begin('boundary_event_redistributed_reflect_requires_uniform')
  call test_boundary_event_redistributed_reflect_requires_uniform()
  call test_end()

  call test_begin('boundary_event_open_precedes_redistributed_reflect')
  call test_boundary_event_open_precedes_redistributed_reflect()
  call test_end()

  call test_begin('boundary_event_mixed_redistributed_reflect_periodic')
  call test_boundary_event_mixed_redistributed_reflect_periodic()
  call test_end()

  call test_begin('boundary_event_mixed_reflect_periodic')
  call test_boundary_event_mixed_reflect_periodic()
  call test_end()

  call test_begin('periodic_wrap_reentry_fraction_is_positive')
  call test_periodic_wrap_reentry_fraction_is_positive()
  call test_end()

  call test_begin('boundary_event_open_escape')
  call test_boundary_event_open_escape()
  call test_end()

  call test_begin('boundary_event_face_bc_override')
  call test_boundary_event_face_bc_override()
  call test_end()

  call test_begin('boundary_event_start_inward')
  call test_boundary_event_start_inward()
  call test_end()

  call test_begin('boundary_event_start_outward')
  call test_boundary_event_start_outward()
  call test_end()

  call test_begin('boundary_event_endpoint_on_face')
  call test_boundary_event_endpoint_on_face()
  call test_end()

  call test_begin('boundary_event_disabled_box')
  call test_boundary_event_disabled_box()
  call test_end()

  call test_begin('boundary_event_invalid_endpoint')
  call test_boundary_event_invalid_endpoint()
  call test_end()

  call test_begin('boundary_event_outside_start_is_invalid')
  call test_boundary_event_outside_start_is_invalid()
  call test_end()

  call test_begin('boundary_event_invalid_span')
  call test_boundary_event_invalid_span()
  call test_end()

  call test_begin('boundary_event_invalid_action_is_transactional')
  call test_boundary_event_invalid_action_is_transactional()
  call test_end()

  call test_summary()

contains

  subroutine init_event_box(event_cfg)
    type(particle_domain_config), intent(out) :: event_cfg

    event_cfg = particle_domain_config()
    event_cfg%use_box = .true.
    event_cfg%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    event_cfg%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
  end subroutine init_event_box

  subroutine test_first_boundary_event_single_face()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status

    call init_event_box(event_cfg)
    event_cfg%bc_high(1) = bc_reflect
    call find_first_boundary_event( &
      event_cfg, [0.25_dp, 0.5_dp, 0.5_dp], [1.25_dp, 0.5_dp, 0.5_dp], event, status &
      )

    call assert_equal_i32(status, boundary_event_ok, 'single-face event status mismatch')
    call assert_true(event%has_event, 'single-face crossing should produce an event')
    call assert_close_dp(event%fraction, 0.75_dp, 0.0_dp, 'single-face fraction mismatch')
    call assert_equal_i32(event%face_mask, 2_i32, 'single-face mask mismatch')
    call assert_true( &
      all(event%face_bc == [bc_open, bc_reflect, bc_open, bc_open, bc_open, bc_open]), &
      'captured face boundary conditions mismatch' &
      )
  end subroutine test_first_boundary_event_single_face

  subroutine test_first_boundary_event_precedes_later_face()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status

    call init_event_box(event_cfg)
    call find_first_boundary_event( &
      event_cfg, [0.25_dp, 0.25_dp, 0.5_dp], [1.25_dp, 1.75_dp, 0.5_dp], event, status &
      )

    call assert_equal_i32(status, boundary_event_ok, 'multiple-face event status mismatch')
    call assert_true(event%has_event, 'multiple-face chord should produce an event')
    call assert_close_dp(event%fraction, 0.5_dp, 0.0_dp, 'earliest face fraction mismatch')
    call assert_equal_i32(event%face_mask, 8_i32, 'later x-high face must not enter the first-event mask')
  end subroutine test_first_boundary_event_precedes_later_face

  subroutine test_first_boundary_event_no_crossing()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status

    call init_event_box(event_cfg)
    call find_first_boundary_event( &
      event_cfg, [0.25_dp, 0.5_dp, 0.5_dp], [0.75_dp, 0.25_dp, 0.5_dp], event, status &
      )
    call assert_equal_i32(status, boundary_event_ok, 'no-crossing event status mismatch')
    call assert_true(.not. event%has_event, 'inside chord should not produce an event')
    call assert_equal_i32(event%face_mask, 0_i32, 'no-crossing mask should be zero')
  end subroutine test_first_boundary_event_no_crossing

  subroutine test_boundary_event_corner_reflect()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status
    real(dp) :: event_x(3), event_v(3)
    logical :: event_alive, event_escaped

    call init_event_box(event_cfg)
    event_cfg%bc_high(1:2) = bc_reflect
    call find_first_boundary_event( &
      event_cfg, [0.5_dp, 0.5_dp, 0.5_dp], [1.5_dp, 1.5_dp, 0.5_dp], event, status &
      )
    call assert_equal_i32(status, boundary_event_ok, 'corner detection status mismatch')
    call assert_equal_i32(event%face_mask, 10_i32, 'corner event should contain x-high and y-high')

    event_x = [1.0_dp, 1.0_dp, 0.5_dp]
    event_v = [1.0_dp, 2.0_dp, 3.0_dp]
    event_alive = .true.
    call apply_escape_reflect_periodic_event( &
      event_cfg, event, event_x, event_v, event_alive, event_escaped, status &
      )
    call assert_equal_i32(status, boundary_event_ok, 'corner action status mismatch')
    call assert_true(event_alive .and. .not. event_escaped, 'reflect corner should keep particle alive')
    call assert_true(all(event_x(1:2) < 1.0_dp), 'corner reflection should move both coordinates inward')
    call assert_true(all(1.0_dp - event_x(1:2) >= 32.0_dp*epsilon(1.0_dp)), &
                     'corner reflection inset should exceed the event tolerance scale')
    call assert_allclose_1d(event_v, [-1.0_dp, -2.0_dp, 3.0_dp], 0.0_dp, 'corner reflection velocity mismatch')
  end subroutine test_boundary_event_corner_reflect

  subroutine test_boundary_event_redistributed_reflect()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status
    real(dp) :: event_x(3), event_v(3)
    logical :: event_alive, event_escaped

    call init_event_box(event_cfg)
    event_cfg%bc_high(3) = bc_redistributed_reflect
    call find_first_boundary_event( &
      event_cfg, [0.1_dp, 0.2_dp, 0.5_dp], [0.1_dp, 0.2_dp, 1.5_dp], event, status &
      )
    event_x = [0.1_dp, 0.2_dp, 1.0_dp]
    event_v = [4.0_dp, 5.0_dp, 6.0_dp]
    event_alive = .true.
    call apply_escape_reflect_periodic_event( &
      event_cfg, event, event_x, event_v, event_alive, event_escaped, status, &
      redistribution_uniform=[0.25_dp, 0.75_dp, 0.5_dp] &
      )

    call assert_equal_i32(status, boundary_event_ok, 'redistributed reflection status mismatch')
    call assert_true(event_alive .and. .not. event_escaped, 'redistributed reflection should survive')
    call assert_close_dp(event_x(1), 0.25_dp, 1.0e-12_dp, 'x tangential coordinate mismatch')
    call assert_close_dp(event_x(2), 0.75_dp, 1.0e-12_dp, 'y tangential coordinate mismatch')
    call assert_true(event_x(3) < 1.0_dp, 'normal coordinate should be inset from the face')
    call assert_allclose_1d( &
      event_v, [4.0_dp, 5.0_dp, -6.0_dp], 0.0_dp, &
      'redistributed reflection must only reverse normal velocity' &
      )
  end subroutine test_boundary_event_redistributed_reflect

  subroutine test_boundary_event_redistributed_reflect_requires_uniform()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status
    real(dp) :: event_x(3), event_v(3), x_before(3), v_before(3)
    logical :: event_alive, event_escaped

    call init_event_box(event_cfg)
    event_cfg%bc_high(3) = bc_redistributed_reflect
    call find_first_boundary_event( &
      event_cfg, [0.1_dp, 0.2_dp, 0.5_dp], [0.1_dp, 0.2_dp, 1.5_dp], event, status &
      )
    event_x = [0.1_dp, 0.2_dp, 1.0_dp]
    event_v = [4.0_dp, 5.0_dp, 6.0_dp]
    x_before = event_x
    v_before = event_v
    event_alive = .true.
    call apply_escape_reflect_periodic_event( &
      event_cfg, event, event_x, event_v, event_alive, event_escaped, status &
      )

    call assert_equal_i32( &
      status, boundary_event_invalid_geometry, &
      'redistributed reflection without supplied randomness must fail closed' &
      )
    call assert_allclose_1d(event_x, x_before, 0.0_dp, 'failed redistribution must preserve position')
    call assert_allclose_1d(event_v, v_before, 0.0_dp, 'failed redistribution must preserve velocity')
    call assert_true(event_alive, 'failed redistribution must preserve alive state')
  end subroutine test_boundary_event_redistributed_reflect_requires_uniform

  subroutine test_boundary_event_open_precedes_redistributed_reflect()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status
    real(dp) :: event_x(3), event_v(3)
    logical :: event_alive, event_escaped

    call init_event_box(event_cfg)
    event_cfg%bc_high(1) = bc_redistributed_reflect
    event_cfg%bc_high(2) = bc_open
    call find_first_boundary_event( &
      event_cfg, [0.5_dp, 0.5_dp, 0.5_dp], [1.5_dp, 1.5_dp, 0.5_dp], event, status &
      )
    event_x = [1.0_dp, 1.0_dp, 0.5_dp]
    event_v = [1.0_dp, 2.0_dp, 3.0_dp]
    event_alive = .true.
    call apply_escape_reflect_periodic_event( &
      event_cfg, event, event_x, event_v, event_alive, event_escaped, status &
      )

    call assert_equal_i32(status, boundary_event_ok, 'mixed open redistribution status mismatch')
    call assert_true(.not. event_alive .and. event_escaped, 'open face must take precedence without random input')
    call assert_allclose_1d( &
      event_x, [1.0_dp, 1.0_dp, 0.5_dp], 0.0_dp, &
      'open-precedence event must preserve its crossing point' &
      )
  end subroutine test_boundary_event_open_precedes_redistributed_reflect

  subroutine test_boundary_event_mixed_redistributed_reflect_periodic()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status
    real(dp) :: event_x(3), event_v(3)
    logical :: event_alive, event_escaped

    call init_event_box(event_cfg)
    event_cfg%bc_high(1) = bc_redistributed_reflect
    event_cfg%bc_low(2) = bc_periodic
    event_cfg%bc_high(2) = bc_periodic
    call find_first_boundary_event( &
      event_cfg, [0.5_dp, 0.5_dp, 0.5_dp], [1.5_dp, 1.5_dp, 0.5_dp], event, status &
      )
    event_x = [1.0_dp, 1.0_dp, 0.5_dp]
    event_v = [1.0_dp, 2.0_dp, 3.0_dp]
    event_alive = .true.
    call apply_escape_reflect_periodic_event( &
      event_cfg, event, event_x, event_v, event_alive, event_escaped, status, &
      redistribution_uniform=[0.25_dp, 0.5_dp, 0.75_dp] &
      )

    call assert_equal_i32(status, boundary_event_ok, 'mixed redistributed-periodic status mismatch')
    call assert_true(event_alive .and. .not. event_escaped, 'mixed redistributed-periodic event must survive')
    call assert_true(event_x(1) < 1.0_dp, 'redistributed-reflect normal coordinate must be inset')
    call assert_true(event_x(2) > 0.0_dp, 'periodic event coordinate must wrap inside the opposite face')
    call assert_close_dp(event_x(3), 0.75_dp, 1.0e-12_dp, 'only the non-event axis should be redistributed')
    call assert_allclose_1d( &
      event_v, [-1.0_dp, 2.0_dp, 3.0_dp], 0.0_dp, &
      'mixed redistributed-periodic velocity mismatch' &
      )
  end subroutine test_boundary_event_mixed_redistributed_reflect_periodic

  subroutine test_boundary_event_mixed_reflect_periodic()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status
    real(dp) :: event_x(3), event_v(3)
    logical :: event_alive, event_escaped

    call init_event_box(event_cfg)
    event_cfg%bc_high(1) = bc_reflect
    event_cfg%bc_high(2) = bc_periodic
    event_cfg%bc_low(2) = bc_periodic
    call find_first_boundary_event( &
      event_cfg, [0.5_dp, 0.5_dp, 0.5_dp], [1.5_dp, 1.5_dp, 0.5_dp], event, status &
      )
    event_x = [1.0_dp, 1.0_dp, 0.5_dp]
    event_v = [1.0_dp, 2.0_dp, 3.0_dp]
    event_alive = .true.
    call apply_escape_reflect_periodic_event( &
      event_cfg, event, event_x, event_v, event_alive, event_escaped, status &
      )

    call assert_equal_i32(status, boundary_event_ok, 'mixed action status mismatch')
    call assert_true(event_alive .and. .not. event_escaped, 'mixed action should keep particle alive')
    call assert_true(event_x(1) < 1.0_dp, 'reflect coordinate should be inside high face')
    call assert_true(event_x(2) > 0.0_dp, 'periodic coordinate should wrap inside low face')
    call assert_true(event_x(2) >= 32.0_dp*epsilon(1.0_dp), &
                     'periodic low-face inset should not be subnormal')
    call assert_allclose_1d(event_v, [-1.0_dp, 2.0_dp, 3.0_dp], 0.0_dp, 'mixed action velocity mismatch')
  end subroutine test_boundary_event_mixed_reflect_periodic

  subroutine test_periodic_wrap_reentry_fraction_is_positive()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event, reentry_event
    integer(i32) :: status
    real(dp) :: event_x(3), event_v(3)
    logical :: event_alive, event_escaped

    call init_event_box(event_cfg)
    event_cfg%box_max(1) = 1.0e-4_dp
    event_cfg%bc_low(1) = bc_periodic
    event_cfg%bc_high(1) = bc_periodic
    call find_first_boundary_event( &
      event_cfg, [0.9e-4_dp, 0.5_dp, 0.5_dp], [1.1e-4_dp, 0.5_dp, 0.5_dp], event, status &
      )
    event_x = [1.0e-4_dp, 0.5_dp, 0.5_dp]
    event_v = [1.0_dp, 0.0_dp, 0.0_dp]
    event_alive = .true.
    call apply_escape_reflect_periodic_event( &
      event_cfg, event, event_x, event_v, event_alive, event_escaped, status &
      )
    call assert_equal_i32(status, boundary_event_ok, 'periodic wrap status mismatch')
    call find_first_boundary_event( &
      event_cfg, event_x, [-1.0e-4_dp, 0.5_dp, 0.5_dp], reentry_event, status &
      )

    call assert_equal_i32(status, boundary_event_ok, 'periodic reentry status mismatch')
    call assert_true(reentry_event%has_event, 'outward reentry chord should produce an event')
    call assert_true(reentry_event%fraction > 0.0_dp, 'reentry fraction must not underflow to zero')
  end subroutine test_periodic_wrap_reentry_fraction_is_positive

  subroutine test_boundary_event_open_escape()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status
    real(dp) :: event_x(3), event_v(3)
    logical :: event_alive, event_escaped

    call init_event_box(event_cfg)
    call find_first_boundary_event( &
      event_cfg, [0.5_dp, 0.5_dp, 0.5_dp], [1.5_dp, 0.5_dp, 0.5_dp], event, status &
      )
    event_x = [1.0_dp, 0.5_dp, 0.5_dp]
    event_v = [1.0_dp, 2.0_dp, 3.0_dp]
    event_alive = .true.
    call apply_escape_reflect_periodic_event( &
      event_cfg, event, event_x, event_v, event_alive, event_escaped, status &
      )

    call assert_equal_i32(status, boundary_event_ok, 'open event action status mismatch')
    call assert_true(.not. event_alive, 'open event should terminate the particle')
    call assert_true(event_escaped, 'open event should report boundary escape')
    call assert_allclose_1d(event_x, [1.0_dp, 0.5_dp, 0.5_dp], 0.0_dp, 'open event should keep event point')
  end subroutine test_boundary_event_open_escape

  !> face 別境界条件を渡すと、box 設定の open 面を反射として扱える。event の face_bc は渡した値と照合する。
  subroutine test_boundary_event_face_bc_override()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status, face_bc(6)
    real(dp) :: event_x(3), event_v(3)
    logical :: event_alive, event_escaped

    call init_event_box(event_cfg)
    call find_first_boundary_event( &
      event_cfg, [0.5_dp, 0.5_dp, 0.5_dp], [1.5_dp, 0.5_dp, 0.5_dp], event, status &
      )
    face_bc = event%face_bc
    face_bc(2) = bc_reflect
    event%face_bc(2) = bc_reflect
    event_x = [1.0_dp, 0.5_dp, 0.5_dp]
    event_v = [1.0_dp, 2.0_dp, 3.0_dp]
    event_alive = .true.
    call apply_escape_reflect_periodic_event( &
      event_cfg, event, event_x, event_v, event_alive, event_escaped, status, face_bc=face_bc &
      )
    call assert_equal_i32(status, boundary_event_ok, 'face_bc override status mismatch')
    call assert_true(event_alive .and. .not. event_escaped, 'overridden reflect face should keep the particle')
    call assert_allclose_1d(event_v, [-1.0_dp, 2.0_dp, 3.0_dp], 0.0_dp, 'overridden face should reflect the normal velocity')
    call assert_true(event_x(1) < 1.0_dp, 'overridden reflect should move the particle inside')

    ! 渡した face 別境界条件と event が食い違えば、box 設定と食い違うときと同じく不正。
    event_x = [1.0_dp, 0.5_dp, 0.5_dp]
    event_v = [1.0_dp, 2.0_dp, 3.0_dp]
    event_alive = .true.
    call apply_escape_reflect_periodic_event( &
      event_cfg, event, event_x, event_v, event_alive, event_escaped, status &
      )
    call assert_equal_i32(status, boundary_event_invalid_geometry, 'event face_bc must match the box without override')
    call assert_true(event_alive, 'rejected event must not change the particle')
  end subroutine test_boundary_event_face_bc_override

  subroutine test_boundary_event_start_inward()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status

    call init_event_box(event_cfg)
    call find_first_boundary_event( &
      event_cfg, [1.0_dp, 0.5_dp, 0.5_dp], [0.5_dp, 0.5_dp, 0.5_dp], event, status &
      )
    call assert_equal_i32(status, boundary_event_ok, 'boundary-start inward status mismatch')
    call assert_true(.not. event%has_event, 'boundary-start inward chord should not create a zero-time event')
  end subroutine test_boundary_event_start_inward

  subroutine test_boundary_event_start_outward()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status

    call init_event_box(event_cfg)
    call find_first_boundary_event( &
      event_cfg, [1.0_dp, 0.5_dp, 0.5_dp], [1.25_dp, 0.5_dp, 0.5_dp], event, status &
      )
    call assert_equal_i32(status, boundary_event_ok, 'boundary-start outward status mismatch')
    call assert_true(event%has_event, 'boundary-start outward chord should create an event')
    call assert_close_dp(event%fraction, 0.0_dp, 0.0_dp, 'boundary-start outward fraction mismatch')
    call assert_equal_i32(event%face_mask, 2_i32, 'boundary-start outward mask mismatch')
  end subroutine test_boundary_event_start_outward

  subroutine test_boundary_event_endpoint_on_face()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status

    call init_event_box(event_cfg)
    call find_first_boundary_event( &
      event_cfg, [0.25_dp, 0.5_dp, 0.5_dp], [1.0_dp, 0.5_dp, 0.5_dp], event, status &
      )
    call assert_equal_i32(status, boundary_event_ok, 'endpoint-on-face status mismatch')
    call assert_true(event%has_event, 'endpoint exactly on an outward face should be an event')
    call assert_close_dp(event%fraction, 1.0_dp, 0.0_dp, 'endpoint-on-face fraction mismatch')
  end subroutine test_boundary_event_endpoint_on_face

  subroutine test_boundary_event_disabled_box()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status

    event_cfg = particle_domain_config()
    call find_first_boundary_event( &
      event_cfg, [0.0_dp, 0.0_dp, 0.0_dp], [10.0_dp, 0.0_dp, 0.0_dp], event, status &
      )
    call assert_equal_i32(status, boundary_event_ok, 'disabled-box event status mismatch')
    call assert_true(.not. event%has_event, 'disabled box should not create boundary events')
  end subroutine test_boundary_event_disabled_box

  subroutine test_boundary_event_invalid_endpoint()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status
    real(dp) :: invalid_x1(3)

    call init_event_box(event_cfg)
    invalid_x1 = [ieee_value(0.0_dp, ieee_quiet_nan), 0.5_dp, 0.5_dp]
    call find_first_boundary_event(event_cfg, [0.5_dp, 0.5_dp, 0.5_dp], invalid_x1, event, status)
    call assert_equal_i32(status, boundary_event_invalid_geometry, 'NaN endpoint should be invalid geometry')
    call assert_true(.not. event%has_event, 'invalid endpoint must not return an event')
  end subroutine test_boundary_event_invalid_endpoint

  subroutine test_boundary_event_outside_start_is_invalid()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status

    call init_event_box(event_cfg)
    call find_first_boundary_event( &
      event_cfg, [1.25_dp, 0.5_dp, 0.5_dp], [0.5_dp, 0.5_dp, 0.5_dp], event, status &
      )
    call assert_equal_i32(status, boundary_event_invalid_geometry, 'outside start should be invalid geometry')
    call assert_true(.not. event%has_event, 'outside start must not synthesize a reentry event')
  end subroutine test_boundary_event_outside_start_is_invalid

  subroutine test_boundary_event_invalid_span()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status

    call init_event_box(event_cfg)
    event_cfg%box_max(1) = event_cfg%box_min(1)
    call find_first_boundary_event( &
      event_cfg, [0.0_dp, 0.5_dp, 0.5_dp], [0.5_dp, 0.5_dp, 0.5_dp], event, status &
      )
    call assert_equal_i32(status, boundary_event_invalid_geometry, 'non-positive span should be invalid geometry')
  end subroutine test_boundary_event_invalid_span

  subroutine test_boundary_event_invalid_action_is_transactional()
    type(particle_domain_config) :: event_cfg
    type(boundary_event_type) :: event
    integer(i32) :: status
    real(dp) :: event_x(3), event_v(3), x_before(3), v_before(3)
    logical :: event_alive, event_escaped

    call init_event_box(event_cfg)
    event_cfg%bc_high(1) = bc_reflect
    call find_first_boundary_event( &
      event_cfg, [0.5_dp, 0.5_dp, 0.5_dp], [1.5_dp, 0.5_dp, 0.5_dp], event, status &
      )
    event_cfg%bc_high(1) = bc_periodic
    event_x = [1.0_dp, 0.5_dp, 0.5_dp]
    event_v = [1.0_dp, 2.0_dp, 3.0_dp]
    x_before = event_x
    v_before = event_v
    event_alive = .true.
    call apply_escape_reflect_periodic_event( &
      event_cfg, event, event_x, event_v, event_alive, event_escaped, status &
      )
    call assert_equal_i32(status, boundary_event_invalid_geometry, 'event/config mismatch should be invalid')
    call assert_true(all(event_x == x_before), 'invalid action must preserve position')
    call assert_true(all(event_v == v_before), 'invalid action must preserve velocity')
    call assert_true(event_alive, 'invalid action must preserve alive state')

  end subroutine test_boundary_event_invalid_action_is_transactional

end program test_boundary

!> ボックス境界条件の各モードを検証するテスト。
program test_boundary
  use bem_kinds, only: dp
  use bem_boundary, only: apply_box_boundary
  use bem_types, only: sim_config, bc_open, bc_reflect, bc_periodic
  use test_support, only: assert_true, assert_close_dp, assert_allclose_1d
  implicit none

  type(sim_config) :: cfg
  real(dp) :: x(3), v(3), expected(3)
  logical :: alive, escaped_boundary

  cfg = sim_config()
  cfg%use_box = .true.
  cfg%box_min = [0.0d0, 0.0d0, 0.0d0]
  cfg%box_max = [1.0d0, 1.0d0, 1.0d0]

  x = [1.1d0, 0.5d0, 0.5d0]
  v = [1.0d0, 2.0d0, 3.0d0]
  alive = .true.
  cfg%bc_high(1) = bc_open
  call apply_box_boundary(cfg, x, v, alive, escaped_boundary)
  call assert_true(.not. alive, 'open boundary should kill the particle')
  call assert_true(escaped_boundary, 'open boundary should mark escaped_boundary')

  x = [1.2d0, 0.5d0, 0.5d0]
  v = [1.0d0, 2.0d0, 3.0d0]
  alive = .true.
  cfg%bc_high(1) = bc_reflect
  call apply_box_boundary(cfg, x, v, alive, escaped_boundary)
  expected = [0.8d0, 0.5d0, 0.5d0]
  call assert_true(alive, 'reflect boundary should keep particle alive')
  call assert_true(.not. escaped_boundary, 'reflect boundary should not mark escaped')
  call assert_allclose_1d(x, expected, 1.0d-10, 'reflect position mismatch')
  call assert_close_dp(v(1), -1.0d0, 1.0d-12, 'reflect velocity mismatch')
  call assert_close_dp(v(2), 2.0d0, 1.0d-12, 'reflect should preserve tangential velocity')

  x = [-0.2d0, 0.5d0, 0.5d0]
  v = [1.0d0, 2.0d0, 3.0d0]
  alive = .true.
  cfg%bc_low(1) = bc_periodic
  call apply_box_boundary(cfg, x, v, alive, escaped_boundary)
  expected = [0.8d0, 0.5d0, 0.5d0]
  call assert_true(alive, 'periodic boundary should keep particle alive')
  call assert_true(.not. escaped_boundary, 'periodic boundary should not mark escaped')
  call assert_allclose_1d(x, expected, 1.0d-10, 'periodic position mismatch')
  call assert_allclose_1d(v, [1.0d0, 2.0d0, 3.0d0], 1.0d-12, 'periodic should preserve velocity')

  cfg = sim_config()
  x = [1.2d0, 0.5d0, 0.5d0]
  v = [1.0d0, 2.0d0, 3.0d0]
  alive = .true.
  call apply_box_boundary(cfg, x, v, alive, escaped_boundary)
  call assert_true(alive, 'disabled box should not change alive flag')
  call assert_true(.not. escaped_boundary, 'disabled box should not mark escaped')
  call assert_allclose_1d(x, [1.2d0, 0.5d0, 0.5d0], 1.0d-12, 'disabled box should preserve position')

  cfg = sim_config()
  cfg%use_box = .true.
  cfg%box_min = [0.0d0, 0.0d0, 0.0d0]
  cfg%box_max = [0.0d0, 1.0d0, 1.0d0]
  x = [0.0d0, 0.5d0, 0.5d0]
  v = [1.0d0, 2.0d0, 3.0d0]
  alive = .true.
  call apply_box_boundary(cfg, x, v, alive, escaped_boundary)
  call assert_true(.not. alive, 'degenerate box should stop the particle')
  call assert_true(escaped_boundary, 'degenerate box should mark escaped_boundary')

end program test_boundary

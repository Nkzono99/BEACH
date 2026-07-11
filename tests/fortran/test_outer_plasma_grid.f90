program test_outer_plasma_grid
  use bem_kinds, only: dp, i32
  use bem_outer_plasma_grid, only: outer_plasma_grid_type, init_outer_plasma_grid, interpolate_outer_profile
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp
  implicit none

  type(outer_plasma_grid_type) :: grid
  real(dp) :: value
  integer(i32) :: j

  call test_init(3)

  call test_begin('stretched grid has exact endpoints and positive cells')
  call init_outer_plasma_grid(9_i32, 12.0_dp, 2.0_dp, grid)
  call assert_close_dp(grid%z(1), 0.0_dp, 0.0_dp, 'grid origin mismatch')
  call assert_close_dp(grid%z(9), 12.0_dp, 1.0e-14_dp, 'grid endpoint mismatch')
  call assert_true(all(grid%dz > 0.0_dp), 'grid cells must be positive')
  do j = 2_i32, grid%n - 1_i32
    call assert_true(grid%dz(j) > grid%dz(j - 1_i32), 'positive stretch must increase cell width')
  end do
  call test_end()

  call test_begin('zero stretch is the uniform-grid limit')
  call init_outer_plasma_grid(5_i32, 4.0_dp, 0.0_dp, grid)
  call assert_close_dp(grid%z(2), 1.0_dp, 1.0e-15_dp, 'uniform node 2 mismatch')
  call assert_close_dp(grid%z(3), 2.0_dp, 1.0e-15_dp, 'uniform node 3 mismatch')
  call assert_close_dp(grid%z(4), 3.0_dp, 1.0e-15_dp, 'uniform node 4 mismatch')
  call test_end()

  call test_begin('profile interpolation is piecewise linear and bounded')
  call interpolate_outer_profile(grid, grid%z**2, 1.5_dp, value)
  call assert_close_dp(value, 2.5_dp, 1.0e-15_dp, 'linear interpolation mismatch')
  call interpolate_outer_profile(grid, grid%z**2, -1.0_dp, value)
  call assert_close_dp(value, 0.0_dp, 0.0_dp, 'left clamp mismatch')
  call interpolate_outer_profile(grid, grid%z**2, 8.0_dp, value)
  call assert_close_dp(value, 16.0_dp, 0.0_dp, 'right clamp mismatch')
  call test_end()

  call test_summary()
end program test_outer_plasma_grid

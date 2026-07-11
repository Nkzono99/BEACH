program test_outer_plasma_unified
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok, outer_plasma_not_applicable
  use bem_outer_plasma_unified, only: unified_outer_linear_options_type, solve_unified_outer_linear
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, assert_equal_i32
  implicit none

  type(unified_outer_linear_options_type) :: options
  type(outer_plasma_state_type) :: state, coarse_state
  real(dp), allocatable :: z(:), surface_field(:), accessible(:)
  integer(i32) :: status
  character(len=256) :: message

  call test_init(3)

  call test_begin('uniform screened half-space matches the exponential solution')
  call build_uniform_inputs(65_i32, 8.0_dp, z, surface_field, accessible)
  options%kappa = 1.0_dp
  options%tail_length = 1.0_dp
  options%bottom_field = 0.2_dp
  options%thermal_voltage = 1.0_dp
  options%max_linearity_ratio = 1.0_dp
  call solve_unified_outer_linear(z, surface_field, accessible, options, state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'uniform unified solve status mismatch')
  call assert_close_dp(state%potential(1), 0.2_dp, 1.0e-3_dp, 'screened interface potential mismatch')
  call assert_close_dp(state%field(1), 0.2_dp, 0.0_dp, 'screened interface field mismatch')
  call assert_close_dp(state%potential(17), 0.2_dp*exp(-2.0_dp), 2.0e-3_dp, 'screened profile mismatch')
  call test_end()

  call test_begin('surface projection and accessible plasma satisfy global Gauss closure')
  surface_field = 0.0_dp
  where (z >= 1.0_dp) surface_field = 0.15_dp
  accessible = 0.0_dp
  where (z > 1.0_dp) accessible = 1.0_dp
  options%bottom_field = 0.0_dp
  call solve_unified_outer_linear(z, surface_field, accessible, options, state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'surface unified solve status mismatch')
  call assert_close_dp( &
    state%integrated_charge_per_area + eps0*(surface_field(size(surface_field)) - surface_field(1)), &
    0.0_dp, 1.0e-24_dp, 'surface plus plasma charge must close at zero far and bottom fields' &
    )
  call test_end()

  call test_begin('grid refinement converges and linearity fails closed')
  call build_uniform_inputs(33_i32, 8.0_dp, z, surface_field, accessible)
  options%bottom_field = 0.2_dp
  call solve_unified_outer_linear(z, surface_field, accessible, options, coarse_state, status, message)
  call build_uniform_inputs(65_i32, 8.0_dp, z, surface_field, accessible)
  call solve_unified_outer_linear(z, surface_field, accessible, options, state, status, message)
  call assert_close_dp(state%potential(1), coarse_state%potential(1), 2.0e-3_dp, 'unified grid refinement mismatch')
  options%max_linearity_ratio = 0.1_dp
  call solve_unified_outer_linear(z, surface_field, accessible, options, state, status, message)
  call assert_equal_i32(status, outer_plasma_not_applicable, 'nonlinear rough response must fail closed')
  call assert_true(.not. state%ready, 'nonlinear rough state must not be ready')
  call test_end()

  call test_summary()

contains

  subroutine build_uniform_inputs(n, length, grid, projected_field, fraction)
    integer(i32), intent(in) :: n
    real(dp), intent(in) :: length
    real(dp), allocatable, intent(out) :: grid(:), projected_field(:), fraction(:)
    integer(i32) :: point

    allocate (grid(n), projected_field(n), fraction(n))
    do point = 1_i32, n
      grid(point) = length*real(point - 1_i32, dp)/real(n - 1_i32, dp)
    end do
    projected_field = 0.0_dp
    fraction = 1.0_dp
  end subroutine build_uniform_inputs

end program test_outer_plasma_unified

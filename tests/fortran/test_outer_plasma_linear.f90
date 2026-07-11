program test_outer_plasma_linear
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok, outer_plasma_not_applicable
  use bem_outer_plasma_linear, only: init_outer_plasma_linear, eval_outer_plasma_linear, &
                                     outer_plasma_integrated_charge_per_area, outer_plasma_gauss_residual_per_area
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_close_dp, assert_equal_i32
  implicit none

  type(outer_plasma_state_type) :: state
  real(dp) :: potential, field, charge_density, integrated_charge, residual
  integer(i32) :: status
  character(len=128) :: message

  call test_init(3)

  call test_begin('linear_debye_profile')
  call init_outer_plasma_linear( &
    interface_z=2.0_dp, interface_potential=3.0_dp, infinity_potential=1.0_dp, &
    debye_length=0.5_dp, linearity_ratio=0.2_dp, max_linearity_ratio=0.5_dp, &
    state=state, status=status, message=message &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'linear state status mismatch')
  call eval_outer_plasma_linear(state, 2.5_dp, potential, field, charge_density)
  call assert_close_dp(potential, 1.0_dp + 2.0_dp*exp(-1.0_dp), 1.0e-14_dp, 'potential mismatch')
  call assert_close_dp(field, 4.0_dp*exp(-1.0_dp), 1.0e-14_dp, 'field mismatch')
  call assert_close_dp(charge_density, -8.0_dp*eps0*exp(-1.0_dp), 1.0e-23_dp, 'charge density mismatch')
  call test_end()

  call test_begin('integrated_charge_and_gauss_residual')
  integrated_charge = outer_plasma_integrated_charge_per_area(state)
  residual = outer_plasma_gauss_residual_per_area(state)
  call assert_close_dp(integrated_charge, -4.0_dp*eps0, 1.0e-25_dp, 'integrated charge mismatch')
  call assert_close_dp(residual, 0.0_dp, 1.0e-25_dp, 'Gauss residual mismatch')
  call test_end()

  call test_begin('linearity_guard')
  call init_outer_plasma_linear( &
    interface_z=0.0_dp, interface_potential=5.0_dp, infinity_potential=0.0_dp, &
    debye_length=1.0_dp, linearity_ratio=1.5_dp, max_linearity_ratio=0.5_dp, &
    state=state, status=status, message=message &
    )
  call assert_equal_i32(status, outer_plasma_not_applicable, 'strong perturbation must be not_applicable')
  call test_end()

  call test_summary()
end program test_outer_plasma_linear

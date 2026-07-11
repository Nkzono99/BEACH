program test_outer_plasma_orbit
  use bem_kinds, only: dp, i32
  use bem_interface_types, only: interface_crossing_type, interface_particle_outcome_type, &
                                 interface_outcome_returned_local, interface_outcome_escaped_to_infinity
  use bem_outer_plasma_orbit, only: explicit_outer_orbit_options_type, trace_electrostatic_outer_orbit
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, assert_equal_i32
  implicit none

  type(interface_crossing_type) :: crossing
  type(interface_particle_outcome_type) :: outcome, coarse_outcome
  type(explicit_outer_orbit_options_type) :: options

  call test_init(3)
  crossing%has_crossing = .true.
  crossing%face_index = 6_i32
  crossing%position = [1.0_dp, 2.0_dp, 0.0_dp]
  crossing%velocity = [0.5_dp, -0.25_dp, 1.0_dp]
  options%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
  options%box_max = [10.0_dp, 10.0_dp, 0.0_dp]
  options%interface_z = 0.0_dp
  options%far_z = 10.0_dp
  options%dt = 0.1_dp
  options%max_steps = 200_i32
  options%field_evolution_timescale = 100.0_dp
  options%max_frozen_field_ratio = 1.0_dp
  options%energy_scale = 1.0_dp
  options%max_energy_relative_error = 1.0e-10_dp

  call test_begin('constant_3d_field_return_preserves_energy_and_lateral_shift')
  call trace_electrostatic_outer_orbit(constant_return_field, options, 1.0_dp, 1.0_dp, crossing, outcome)
  call assert_equal_i32(outcome%kind, interface_outcome_returned_local, 'constant field orbit must return')
  call assert_close_dp(outcome%outer_flight_time, 2.0_dp, 2.0e-13_dp, 'constant field flight time mismatch')
  call assert_close_dp(outcome%position(1), 2.4_dp, 2.0e-13_dp, 'constant field lateral shift mismatch')
  call assert_close_dp(outcome%position(2), 1.5_dp, 2.0e-13_dp, 'constant field y shift mismatch')
  call assert_close_dp(outcome%velocity(1), 0.9_dp, 2.0e-13_dp, 'constant field lateral acceleration mismatch')
  call assert_close_dp(outcome%velocity(3), -1.0_dp, 2.0e-13_dp, 'constant field return velocity mismatch')
  call assert_true(outcome%energy_relative_error <= 1.0e-12_dp, 'constant field energy must be conserved')
  call test_end()

  call test_begin('far_plane_crossing_escapes')
  crossing%position = [1.0_dp, 2.0_dp, 0.0_dp]
  crossing%velocity = [0.0_dp, 0.0_dp, 1.0_dp]
  options%far_z = 2.0_dp
  call trace_electrostatic_outer_orbit(zero_field, options, 1.0_dp, 1.0_dp, crossing, outcome)
  call assert_equal_i32(outcome%kind, interface_outcome_escaped_to_infinity, 'outward far crossing must escape')
  call assert_close_dp(outcome%outer_flight_time, 2.0_dp, 2.0e-13_dp, 'escape flight time mismatch')
  call test_end()

  call test_begin('nonlinear_field_energy_error_converges_with_dt')
  crossing%position = [1.0_dp, 2.0_dp, 1.0_dp]
  crossing%velocity = [0.0_dp, 0.0_dp, 0.5_dp]
  options%interface_z = 1.0_dp
  options%box_max(3) = 1.0_dp
  options%far_z = 10.0_dp
  options%max_energy_relative_error = 1.0_dp
  options%dt = 0.1_dp
  call trace_electrostatic_outer_orbit(harmonic_field, options, 1.0_dp, 1.0_dp, crossing, coarse_outcome)
  options%dt = 0.05_dp
  call trace_electrostatic_outer_orbit(harmonic_field, options, 1.0_dp, 1.0_dp, crossing, outcome)
  call assert_equal_i32(outcome%kind, interface_outcome_returned_local, 'harmonic orbit must return')
  call assert_true( &
    outcome%energy_relative_error < coarse_outcome%energy_relative_error, &
    'halving outer orbit dt must reduce the energy error' &
    )
  call test_end()
  call test_summary()

contains

  subroutine constant_return_field(position, potential, electric_field)
    real(dp), intent(in) :: position(3)
    real(dp), intent(out) :: potential, electric_field(3)

    potential = -0.2_dp*position(1) + position(3)
    electric_field = [0.2_dp, 0.0_dp, -1.0_dp]
  end subroutine constant_return_field

  subroutine zero_field(position, potential, electric_field)
    real(dp), intent(in) :: position(3)
    real(dp), intent(out) :: potential, electric_field(3)

    potential = 0.0_dp*position(1)
    electric_field = 0.0_dp
  end subroutine zero_field

  subroutine harmonic_field(position, potential, electric_field)
    real(dp), intent(in) :: position(3)
    real(dp), intent(out) :: potential, electric_field(3)

    potential = 0.5_dp*position(3)*position(3)
    electric_field = [0.0_dp, 0.0_dp, -position(3)]
  end subroutine harmonic_field

end program test_outer_plasma_orbit

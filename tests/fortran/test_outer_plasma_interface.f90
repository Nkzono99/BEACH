program test_outer_plasma_interface
  use bem_kinds, only: dp, i32
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok
  use bem_outer_plasma_linear, only: init_outer_plasma_linear
  use bem_interface_types, only: interface_crossing_type, interface_particle_outcome_type, &
                                 interface_outcome_returned_local, interface_outcome_escaped_to_infinity, &
                                 interface_outcome_queued_outer, interface_outcome_invalid_model
  use bem_outer_plasma_interface, only: map_infinity_normal_velocity_to_interface, map_outer_particle_linear_debye, &
                                        map_outer_particle_kinetic_profile
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, assert_equal_i32
  implicit none

  type(outer_plasma_state_type) :: state
  type(interface_crossing_type) :: crossing
  type(interface_particle_outcome_type) :: outcome
  integer(i32) :: status
  character(len=128) :: message
  logical :: accessible
  real(dp) :: mapped_speed

  call test_init(13)
  call init_outer_plasma_linear( &
    interface_z=1.0_dp, interface_potential=-1.0_dp, infinity_potential=0.0_dp, &
    debye_length=0.5_dp, linearity_ratio=0.1_dp, max_linearity_ratio=1.0_dp, &
    state=state, status=status, message=message &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'outer state initialization failed')

  call test_begin('infinity_to_interface_energy_map')
  call map_infinity_normal_velocity_to_interface(state, 1.0_dp, 1.0_dp, 1.0_dp, accessible, mapped_speed)
  call assert_true(accessible, 'attracted inflow should reach the interface')
  call assert_close_dp(mapped_speed, sqrt(3.0_dp), 1.0e-14_dp, 'inflow energy map mismatch')
  call map_infinity_normal_velocity_to_interface(state, -1.0_dp, 1.0_dp, 1.0_dp, accessible, mapped_speed)
  call assert_true(.not. accessible, 'sub-threshold repelled inflow must be rejected')
  call test_end()

  crossing%has_crossing = .true.
  crossing%face_index = 6_i32
  crossing%fraction = 0.25_dp
  crossing%position = [0.9_dp, 0.8_dp, 1.0_dp]
  crossing%velocity = [0.5_dp, 0.25_dp, 1.0_dp]
  crossing%dt_remaining = 0.75_dp

  call test_begin('kinetic_profile_constant_field_return')
  call set_kinetic_profile(state, [-1.0_dp, 0.0_dp], [1.0_dp, 2.0_dp])
  crossing%velocity = [0.5_dp, 0.25_dp, 1.0_dp]
  call map_outer_particle_kinetic_profile( &
    state, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 1.0_dp, 1.0_dp, crossing, &
    field_timescale=10.0_dp, max_frozen_field_ratio=1.0_dp, queue_enabled=.false., outcome=outcome &
    )
  call assert_equal_i32(outcome%kind, interface_outcome_returned_local, 'kinetic barrier must return the particle')
  call assert_close_dp(outcome%outer_flight_time, 2.0_dp, 1.0e-14_dp, 'kinetic flight time mismatch')
  call assert_close_dp(outcome%position(1), modulo(0.9_dp + 1.0_dp, 1.0_dp), &
                       1.0e-14_dp, 'kinetic lateral x shift mismatch')
  call assert_close_dp(outcome%position(2), modulo(0.8_dp + 0.5_dp, 1.0_dp), &
                       1.0e-14_dp, 'kinetic lateral y shift mismatch')
  call assert_close_dp(outcome%normal_energy_residual, 0.0_dp, 1.0e-14_dp, 'kinetic return energy mismatch')
  call test_end()

  call test_begin('kinetic_profile_robin_tail_return')
  call set_kinetic_profile(state, [-1.0_dp, -0.5_dp], [1.0_dp, 2.0_dp])
  state%infinity_potential = 0.0_dp
  state%debye_length = 1.0_dp
  crossing%velocity = [0.0_dp, 0.0_dp, 1.2_dp]
  call map_outer_particle_kinetic_profile( &
    state, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 1.0_dp, 1.0_dp, crossing, &
    field_timescale=100.0_dp, max_frozen_field_ratio=1.0_dp, queue_enabled=.false., outcome=outcome &
    )
  call assert_equal_i32(outcome%kind, interface_outcome_returned_local, &
                        'Robin-tail barrier must return the particle')
  call assert_close_dp( &
    outcome%outer_flight_time, &
    2.0_dp*(2.0_dp/(1.2_dp + sqrt(0.44_dp)) + &
            2.0_dp/sqrt(0.56_dp)*atan(sqrt(0.44_dp/0.56_dp))), &
    1.0e-13_dp, 'Robin-tail flight time mismatch' &
    )
  call test_end()

  call test_begin('kinetic_profile_escape')
  crossing%velocity(3) = 2.0_dp
  call map_outer_particle_kinetic_profile( &
    state, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 1.0_dp, 1.0_dp, crossing, &
    field_timescale=10.0_dp, max_frozen_field_ratio=1.0_dp, queue_enabled=.false., outcome=outcome &
    )
  call assert_equal_i32(outcome%kind, interface_outcome_escaped_to_infinity, &
                        'kinetic super-threshold particle must escape')
  call test_end()

  call test_begin('kinetic_profile_return_is_queued_with_target_state')
  call set_kinetic_profile(state, [-1.0_dp, 0.0_dp], [1.0_dp, 2.0_dp])
  crossing%velocity = [0.0_dp, 0.0_dp, 1.2_dp]
  call map_outer_particle_kinetic_profile( &
    state, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 1.0_dp, 1.0_dp, crossing, &
    field_timescale=100.0_dp, max_frozen_field_ratio=1.0_dp, queue_enabled=.true., outcome=outcome, &
    queue_poll_interval=1.0_dp &
    )
  call assert_equal_i32(outcome%kind, interface_outcome_queued_outer, 'kinetic return must enter persistent queue')
  call assert_equal_i32(outcome%queued_terminal_kind, interface_outcome_returned_local, &
                        'queued kinetic return must retain an explicit terminal kind')
  call assert_true(outcome%outer_flight_time > 0.0_dp, 'queued return must retain a positive due-time offset')
  call assert_true(outcome%velocity(3) < 0.0_dp, 'queued return target velocity must point into the local domain')
  call test_end()

  call test_begin('kinetic_robin_tail_return_exits_finite_queue_volume')
  call set_kinetic_profile(state, [-1.0_dp, -0.5_dp], [1.0_dp, 2.0_dp])
  state%infinity_potential = 0.0_dp
  state%debye_length = 1.0_dp
  crossing%velocity = [0.0_dp, 0.0_dp, 1.2_dp]
  call map_outer_particle_kinetic_profile( &
    state, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 1.0_dp, 1.0_dp, crossing, &
    field_timescale=100.0_dp, max_frozen_field_ratio=1.0_dp, queue_enabled=.true., outcome=outcome, &
    queue_poll_interval=1.0_dp &
    )
  call assert_equal_i32(outcome%kind, interface_outcome_queued_outer, &
                        'finite-volume departure must enter the persistent queue')
  call assert_equal_i32(outcome%queued_terminal_kind, interface_outcome_escaped_to_infinity, &
                        'a particle reaching L must leave the finite queue volume')
  call assert_close_dp(outcome%position(3), state%z(state%profile_n), 1.0e-14_dp, &
                       'finite-volume departure endpoint mismatch')
  call test_end()

  call test_begin('kinetic_profile_escape_is_queued_to_profile_end')
  crossing%velocity = [0.0_dp, 0.0_dp, 2.0_dp]
  call map_outer_particle_kinetic_profile( &
    state, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 1.0_dp, 1.0_dp, crossing, &
    field_timescale=10.0_dp, max_frozen_field_ratio=1.0_dp, queue_enabled=.true., outcome=outcome, &
    queue_poll_interval=1.0_dp &
    )
  call assert_equal_i32(outcome%kind, interface_outcome_queued_outer, 'kinetic escape must enter persistent queue')
  call assert_equal_i32(outcome%queued_terminal_kind, interface_outcome_escaped_to_infinity, &
                        'queued kinetic escape must retain an explicit terminal kind')
  call assert_close_dp(outcome%position(3), state%z(state%profile_n), 1.0e-14_dp, &
                       'queued escape endpoint must be the finite profile boundary')
  call assert_close_dp(outcome%velocity(3), sqrt(3.0_dp), 1.0e-14_dp, &
                       'queued escape target speed mismatch')
  call test_end()

  call test_begin('kinetic_queue_frozen_limit_includes_batch_poll_delay')
  call map_outer_particle_kinetic_profile( &
    state, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 1.0_dp, 1.0_dp, crossing, &
    field_timescale=10.0_dp, max_frozen_field_ratio=0.1_dp, queue_enabled=.true., outcome=outcome, &
    queue_poll_interval=1.0_dp &
    )
  call assert_true(outcome%outer_flight_time/10.0_dp < 0.1_dp, &
                   'fixture flight alone must fit the frozen-field bound')
  call assert_close_dp(outcome%frozen_field_ratio, 0.20_dp, 1.0e-14_dp, &
                       'frozen-field ratio must include poll delay and midpoint uncertainty')
  call assert_equal_i32(outcome%kind, interface_outcome_invalid_model, &
                        'queue poll quantization must not bypass the frozen-field limit')
  call test_end()

  call test_begin('kinetic_profile_nonmonotonic_barrier_returns_particle')
  call set_kinetic_profile(state, [0.0_dp, -1.0_dp, 0.0_dp], [1.0_dp, 2.0_dp, 3.0_dp])
  state%infinity_potential = 0.0_dp
  crossing%velocity(3) = 1.0_dp
  call map_outer_particle_kinetic_profile( &
    state, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], -1.0_dp, 1.0_dp, crossing, &
    field_timescale=10.0_dp, max_frozen_field_ratio=1.0_dp, queue_enabled=.false., outcome=outcome &
    )
  call assert_equal_i32(outcome%kind, interface_outcome_invalid_model, &
                        'default kinetic closure must reject a nonmonotonic profile')
  state%kinetic_closure = 'zhao_charge_driven'
  state%zhao_branch = 'A'
  call map_outer_particle_kinetic_profile( &
    state, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], -1.0_dp, 1.0_dp, crossing, &
    field_timescale=10.0_dp, max_frozen_field_ratio=1.0_dp, queue_enabled=.false., outcome=outcome &
    )
  call assert_equal_i32(outcome%kind, interface_outcome_returned_local, &
                        'an interior nonmonotonic barrier must return the particle')
  call assert_close_dp(outcome%outer_flight_time, 2.0_dp, 1.0e-14_dp, &
                       'nonmonotonic turning-point flight time mismatch')
  call test_end()

  call test_begin('kinetic_profile_nonmonotonic_barrier_filters_inflow')
  call map_infinity_normal_velocity_to_interface(state, -1.0_dp, 1.0_dp, 1.0_dp, accessible, mapped_speed)
  call assert_true(.not. accessible, 'an infinity particle below the interior barrier must be rejected')
  call map_infinity_normal_velocity_to_interface(state, -1.0_dp, 1.0_dp, 2.0_dp, accessible, mapped_speed)
  call assert_true(accessible, 'an infinity particle above the interior barrier must reach the interface')
  call assert_close_dp(mapped_speed, 2.0_dp, 1.0e-14_dp, 'nonmonotonic inflow energy map mismatch')
  call test_end()

  call init_outer_plasma_linear( &
    interface_z=1.0_dp, interface_potential=-1.0_dp, infinity_potential=0.0_dp, &
    debye_length=0.5_dp, linearity_ratio=0.1_dp, max_linearity_ratio=1.0_dp, &
    state=state, status=status, message=message &
    )
  crossing%velocity = [0.5_dp, 0.25_dp, 1.0_dp]

  call test_begin('instant_return_energy_time_and_shift')
  call map_outer_particle_linear_debye( &
    state, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 1.0_dp, 1.0_dp, crossing, &
    field_timescale=10.0_dp, max_frozen_field_ratio=1.0_dp, queue_enabled=.false., outcome=outcome &
    )
  call assert_equal_i32(outcome%kind, interface_outcome_returned_local, 'barrier must return the particle')
  call assert_close_dp(outcome%outer_flight_time, 0.5_dp*acos(-1.0_dp), 1.0e-14_dp, 'outer flight time mismatch')
  call assert_close_dp(outcome%position(1), modulo(0.9_dp + 0.25_dp*acos(-1.0_dp), 1.0_dp), &
                       1.0e-14_dp, 'lateral x shift mismatch')
  call assert_close_dp(outcome%position(2), modulo(0.8_dp + 0.125_dp*acos(-1.0_dp), 1.0_dp), &
                       1.0e-14_dp, 'lateral y shift mismatch')
  call assert_close_dp(outcome%velocity(3), -1.0_dp, 1.0e-14_dp, 'normal velocity reversal mismatch')
  call assert_close_dp(outcome%normal_energy_residual, 0.0_dp, 1.0e-14_dp, 'return energy residual mismatch')
  call test_end()

  call test_begin('escape_to_infinity')
  crossing%velocity(3) = 2.0_dp
  call map_outer_particle_linear_debye( &
    state, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 1.0_dp, 1.0_dp, crossing, &
    field_timescale=10.0_dp, max_frozen_field_ratio=1.0_dp, queue_enabled=.false., outcome=outcome &
    )
  call assert_equal_i32(outcome%kind, interface_outcome_escaped_to_infinity, 'super-threshold particle must escape')
  call test_end()

  call test_begin('frozen_field_violation_is_rejected_even_with_queue')
  crossing%velocity(3) = 1.0_dp
  call map_outer_particle_linear_debye( &
    state, [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 1.0_dp, 1.0_dp, crossing, &
    field_timescale=1.0_dp, max_frozen_field_ratio=0.1_dp, queue_enabled=.true., outcome=outcome &
    )
  call assert_equal_i32(outcome%kind, interface_outcome_invalid_model, &
                        'persistent queue must not bypass the frozen-field limit')
  call test_end()

  call test_summary()

contains

  subroutine set_kinetic_profile(value, potential, z)
    type(outer_plasma_state_type), intent(out) :: value
    real(dp), intent(in) :: potential(:), z(:)

    value%model = 'kinetic_1d'
    value%kinetic_closure = 'absorbing_maxwellian'
    value%ready = .true.
    value%interface_z = z(1)
    value%interface_potential = potential(1)
    value%infinity_potential = potential(size(potential))
    value%profile_n = int(size(potential), i32)
    allocate (value%potential(size(potential)), value%z(size(z)))
    value%potential = potential
    value%z = z
  end subroutine set_kinetic_profile
end program test_outer_plasma_interface

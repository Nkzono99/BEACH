module bem_outer_plasma_interface
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp
  use bem_outer_plasma_types, only: outer_plasma_state_type
  use bem_interface_types, only: interface_crossing_type, interface_particle_outcome_type, &
                                 interface_outcome_returned_local, interface_outcome_escaped_to_infinity, &
                                 interface_outcome_queued_outer, interface_outcome_invalid_model
  implicit none
  private

  public :: map_infinity_normal_velocity_to_interface
  public :: map_outer_particle_linear_debye

contains

  subroutine map_infinity_normal_velocity_to_interface(state, charge, mass, velocity_infinity, accessible, velocity_interface)
    type(outer_plasma_state_type), intent(in) :: state
    real(dp), intent(in) :: charge, mass, velocity_infinity
    logical, intent(out) :: accessible
    real(dp), intent(out) :: velocity_interface
    real(dp) :: speed_squared

    accessible = .false.
    velocity_interface = 0.0_dp
    if (.not. state%ready .or. mass <= 0.0_dp .or. velocity_infinity < 0.0_dp) return
    speed_squared = velocity_infinity*velocity_infinity - &
                    2.0_dp*charge*(state%interface_potential - state%infinity_potential)/mass
    if (.not. ieee_is_finite(speed_squared) .or. speed_squared < 0.0_dp) return
    velocity_interface = sqrt(speed_squared)
    accessible = .true.
  end subroutine map_infinity_normal_velocity_to_interface

  subroutine map_outer_particle_linear_debye( &
    state, box_min, box_max, charge, mass, crossing, field_timescale, max_frozen_field_ratio, &
    queue_enabled, outcome &
    )
    type(outer_plasma_state_type), intent(in) :: state
    real(dp), intent(in) :: box_min(3), box_max(3), charge, mass
    type(interface_crossing_type), intent(in) :: crossing
    real(dp), intent(in) :: field_timescale, max_frozen_field_ratio
    logical, intent(in) :: queue_enabled
    type(interface_particle_outcome_type), intent(out) :: outcome
    real(dp) :: normal_speed_squared, infinity_speed_squared, deficit, delta_potential
    real(dp) :: energy_before, energy_turn, inward_offset

    outcome = interface_particle_outcome_type()
    if (.not. state%ready .or. .not. crossing%has_crossing .or. crossing%face_index /= 6 .or. &
        mass <= 0.0_dp .or. crossing%velocity(3) <= 0.0_dp .or. &
        field_timescale <= 0.0_dp .or. max_frozen_field_ratio <= 0.0_dp) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'invalid z-high outer-interface crossing'
      return
    end if
    delta_potential = state%interface_potential - state%infinity_potential
    normal_speed_squared = crossing%velocity(3)*crossing%velocity(3)
    infinity_speed_squared = normal_speed_squared + 2.0_dp*charge*delta_potential/mass
    if (.not. ieee_is_finite(infinity_speed_squared)) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'non-finite outer normal energy'
      return
    end if
    if (infinity_speed_squared >= 0.0_dp) then
      outcome%kind = interface_outcome_escaped_to_infinity
      outcome%position = crossing%position
      outcome%velocity = crossing%velocity
      return
    end if

    deficit = -infinity_speed_squared
    outcome%outer_flight_time = 4.0_dp*state%debye_length/sqrt(deficit)* &
                                atan(sqrt(normal_speed_squared/deficit))
    outcome%frozen_field_ratio = outcome%outer_flight_time/field_timescale
    if (.not. ieee_is_finite(outcome%outer_flight_time) .or. &
        outcome%frozen_field_ratio > max_frozen_field_ratio) then
      if (queue_enabled) then
        outcome%kind = interface_outcome_queued_outer
        outcome%message = 'outer flight exceeds frozen-field limit'
      else
        outcome%kind = interface_outcome_invalid_model
        outcome%message = 'outer flight requires a persistent queue'
      end if
      return
    end if

    outcome%kind = interface_outcome_returned_local
    outcome%position = crossing%position + crossing%velocity*outcome%outer_flight_time
    outcome%position(1) = wrap_periodic(outcome%position(1), box_min(1), box_max(1))
    outcome%position(2) = wrap_periodic(outcome%position(2), box_min(2), box_max(2))
    inward_offset = 64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(box_max(3)))
    outcome%position(3) = box_max(3) - inward_offset
    outcome%velocity = crossing%velocity
    outcome%velocity(3) = -crossing%velocity(3)
    energy_before = 0.5_dp*mass*normal_speed_squared + charge*state%interface_potential
    energy_turn = charge*state%interface_potential + 0.5_dp*mass*outcome%velocity(3)**2
    outcome%normal_energy_residual = energy_turn - energy_before
  end subroutine map_outer_particle_linear_debye

  pure real(dp) function wrap_periodic(value, low, high) result(wrapped)
    real(dp), intent(in) :: value, low, high

    wrapped = low + modulo(value - low, high - low)
  end function wrap_periodic

end module bem_outer_plasma_interface

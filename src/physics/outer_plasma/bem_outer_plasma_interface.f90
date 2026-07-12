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
  public :: map_outer_particle_kinetic_profile

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

  subroutine map_outer_particle_kinetic_profile( &
    state, box_min, box_max, charge, mass, crossing, field_timescale, max_frozen_field_ratio, &
    queue_enabled, outcome &
    )
    type(outer_plasma_state_type), intent(in) :: state
    real(dp), intent(in) :: box_min(3), box_max(3), charge, mass
    type(interface_crossing_type), intent(in) :: crossing
    real(dp), intent(in) :: field_timescale, max_frozen_field_ratio
    logical, intent(in) :: queue_enabled
    type(interface_particle_outcome_type), intent(out) :: outcome
    real(dp) :: speed2, next_speed2, infinity_speed2, outward_time, segment_length, turning_fraction
    real(dp) :: deficit, tolerance
    integer :: point

    outcome = interface_particle_outcome_type()
    if (.not. valid_kinetic_profile_state(state)) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'invalid kinetic profile state at z-high interface'
      return
    end if
    if (.not. crossing%has_crossing .or. crossing%face_index /= 6) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'invalid kinetic z-high crossing metadata'
      return
    end if
    if (mass <= 0.0_dp .or. crossing%velocity(3) <= 0.0_dp) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'invalid kinetic z-high particle state'
      return
    end if
    if (field_timescale <= 0.0_dp .or. max_frozen_field_ratio <= 0.0_dp) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'invalid kinetic outer-interface timescale limits'
      return
    end if

    speed2 = crossing%velocity(3)**2
    infinity_speed2 = speed2 + &
                      2.0_dp*charge*(state%interface_potential - state%infinity_potential)/mass
    if (.not. ieee_is_finite(infinity_speed2)) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'non-finite kinetic outer normal energy'
      return
    end if
    tolerance = 256.0_dp*epsilon(1.0_dp)*max(1.0_dp, speed2, abs(infinity_speed2))
    if (infinity_speed2 >= -tolerance) then
      outcome%kind = interface_outcome_escaped_to_infinity
      outcome%position = crossing%position
      outcome%velocity = crossing%velocity
      return
    end if

    outward_time = 0.0_dp
    do point = 1, state%profile_n - 1
      next_speed2 = crossing%velocity(3)**2 + &
                    2.0_dp*charge*(state%interface_potential - state%potential(point + 1))/mass
      if (.not. ieee_is_finite(next_speed2)) then
        outcome%kind = interface_outcome_invalid_model
        outcome%message = 'non-finite kinetic profile normal energy'
        return
      end if
      segment_length = state%z(point + 1) - state%z(point)
      if (next_speed2 <= 0.0_dp) then
        if (speed2 <= 0.0_dp .or. next_speed2 >= speed2) then
          outcome%kind = interface_outcome_invalid_model
          outcome%message = 'kinetic profile does not bracket a physical turning point'
          return
        end if
        turning_fraction = speed2/(speed2 - next_speed2)
        outward_time = outward_time + 2.0_dp*segment_length*turning_fraction/sqrt(speed2)
        call finish_kinetic_return( &
          state, box_min, box_max, charge, mass, crossing, field_timescale, max_frozen_field_ratio, &
          queue_enabled, 2.0_dp*outward_time, outcome &
          )
        return
      end if
      outward_time = outward_time + 2.0_dp*segment_length/(sqrt(speed2) + sqrt(next_speed2))
      speed2 = next_speed2
    end do

    deficit = -infinity_speed2
    if (state%debye_length <= 0.0_dp .or. speed2 <= 0.0_dp) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'kinetic turning point requires a positive Robin tail'
      return
    end if
    outward_time = outward_time + &
                   2.0_dp*state%debye_length/sqrt(deficit)*atan(sqrt(speed2/deficit))
    call finish_kinetic_return( &
      state, box_min, box_max, charge, mass, crossing, field_timescale, max_frozen_field_ratio, &
      queue_enabled, 2.0_dp*outward_time, outcome &
      )
  end subroutine map_outer_particle_kinetic_profile

  logical function valid_kinetic_profile_state(state) result(valid)
    type(outer_plasma_state_type), intent(in) :: state
    real(dp) :: scale, tolerance
    logical :: nondecreasing, nonincreasing

    valid = state%ready .and. trim(state%model) == 'kinetic_1d' .and. state%profile_n >= 2 .and. &
            allocated(state%z) .and. allocated(state%potential)
    if (.not. valid) return
    valid = size(state%z) == state%profile_n .and. size(state%potential) == state%profile_n
    if (.not. valid) return
    valid = all(ieee_is_finite(state%z)) .and. all(ieee_is_finite(state%potential)) .and. &
            ieee_is_finite(state%interface_potential) .and. ieee_is_finite(state%infinity_potential) .and. &
            all(state%z(2:) > state%z(:state%profile_n - 1))
    if (.not. valid) return
    scale = max(1.0_dp, maxval(abs(state%potential)), abs(state%infinity_potential))
    tolerance = 256.0_dp*epsilon(1.0_dp)*scale
    valid = abs(state%z(1) - state%interface_z) <= tolerance .and. &
            abs(state%potential(1) - state%interface_potential) <= tolerance
    if (.not. valid) return
    nondecreasing = all(state%potential(2:) >= state%potential(:state%profile_n - 1) - tolerance) .and. &
                    state%infinity_potential >= state%potential(state%profile_n) - tolerance
    nonincreasing = all(state%potential(2:) <= state%potential(:state%profile_n - 1) + tolerance) .and. &
                    state%infinity_potential <= state%potential(state%profile_n) + tolerance
    valid = nondecreasing .or. nonincreasing
  end function valid_kinetic_profile_state

  subroutine finish_kinetic_return( &
    state, box_min, box_max, charge, mass, crossing, field_timescale, max_frozen_field_ratio, &
    queue_enabled, flight_time, outcome &
    )
    type(outer_plasma_state_type), intent(in) :: state
    real(dp), intent(in) :: box_min(3), box_max(3), charge, mass, field_timescale, max_frozen_field_ratio
    type(interface_crossing_type), intent(in) :: crossing
    logical, intent(in) :: queue_enabled
    real(dp), intent(in) :: flight_time
    type(interface_particle_outcome_type), intent(out) :: outcome
    real(dp) :: inward_offset, energy_before, energy_after, energy_scale

    outcome = interface_particle_outcome_type()
    outcome%outer_flight_time = flight_time
    outcome%frozen_field_ratio = flight_time/field_timescale
    if (.not. ieee_is_finite(flight_time) .or. flight_time <= 0.0_dp .or. &
        outcome%frozen_field_ratio > max_frozen_field_ratio) then
      if (queue_enabled) then
        outcome%kind = interface_outcome_queued_outer
        outcome%message = 'kinetic outer flight exceeds frozen-field limit'
      else
        outcome%kind = interface_outcome_invalid_model
        outcome%message = 'kinetic outer flight requires a persistent queue'
      end if
      return
    end if

    outcome%kind = interface_outcome_returned_local
    outcome%position = crossing%position + crossing%velocity*flight_time
    outcome%position(1) = wrap_periodic(outcome%position(1), box_min(1), box_max(1))
    outcome%position(2) = wrap_periodic(outcome%position(2), box_min(2), box_max(2))
    inward_offset = 64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(box_max(3)))
    outcome%position(3) = box_max(3) - inward_offset
    outcome%velocity = crossing%velocity
    outcome%velocity(3) = -crossing%velocity(3)
    energy_before = 0.5_dp*mass*crossing%velocity(3)**2 + charge*state%interface_potential
    energy_after = 0.5_dp*mass*outcome%velocity(3)**2 + charge*state%interface_potential
    outcome%normal_energy_residual = energy_after - energy_before
    energy_scale = max(abs(energy_before), 0.5_dp*mass*crossing%velocity(3)**2, tiny(1.0_dp))
    outcome%energy_relative_error = abs(outcome%normal_energy_residual)/energy_scale
  end subroutine finish_kinetic_return

  pure real(dp) function wrap_periodic(value, low, high) result(wrapped)
    real(dp), intent(in) :: value, low, high

    wrapped = low + modulo(value - low, high - low)
  end function wrap_periodic

end module bem_outer_plasma_interface

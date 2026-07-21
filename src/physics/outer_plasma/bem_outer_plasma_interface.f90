module bem_outer_plasma_interface
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp
  use bem_outer_plasma_types, only: outer_plasma_state_type
  use bem_interface_types, only: interface_crossing_type, interface_particle_outcome_type, &
                                 interface_outcome_returned_local, interface_outcome_escaped_to_infinity, &
                                 interface_outcome_queued_outer, interface_outcome_invalid_model
  use bem_string_utils, only: lower_ascii
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
    real(dp) :: speed_squared, profile_speed_squared, tolerance
    integer :: point

    accessible = .false.
    velocity_interface = 0.0_dp
    if (.not. state%ready .or. mass <= 0.0_dp .or. velocity_infinity < 0.0_dp) return
    speed_squared = velocity_infinity*velocity_infinity - &
                    2.0_dp*charge*(state%interface_potential - state%infinity_potential)/mass
    if (.not. ieee_is_finite(speed_squared) .or. speed_squared < 0.0_dp) return
    if (trim(state%model) == 'kinetic_1d') then
      if (.not. valid_kinetic_profile_state(state)) return
      tolerance = 256.0_dp*epsilon(1.0_dp)*max(1.0_dp, velocity_infinity*velocity_infinity, abs(speed_squared))
      do point = state%profile_n, 1, -1
        profile_speed_squared = velocity_infinity*velocity_infinity + &
                                2.0_dp*charge*(state%infinity_potential - state%potential(point))/mass
        if (.not. ieee_is_finite(profile_speed_squared) .or. profile_speed_squared < -tolerance) return
      end do
    end if
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
    if (.not. ieee_is_finite(outcome%outer_flight_time) .or. outcome%outer_flight_time <= 0.0_dp) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'invalid linear-Debye outer flight time'
      return
    end if

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
    if (outcome%frozen_field_ratio > max_frozen_field_ratio) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'linear-Debye outer flight violates the frozen-field limit'
    else if (queue_enabled) then
      outcome%kind = interface_outcome_queued_outer
      outcome%queued_terminal_kind = interface_outcome_returned_local
      outcome%message = 'linear-Debye return scheduled in persistent outer queue'
    else
      outcome%kind = interface_outcome_returned_local
    end if
  end subroutine map_outer_particle_linear_debye

  subroutine map_outer_particle_kinetic_profile( &
    state, box_min, box_max, charge, mass, crossing, field_timescale, max_frozen_field_ratio, &
    queue_enabled, outcome, queue_poll_interval &
    )
    type(outer_plasma_state_type), intent(in) :: state
    real(dp), intent(in) :: box_min(3), box_max(3), charge, mass
    type(interface_crossing_type), intent(in) :: crossing
    real(dp), intent(in) :: field_timescale, max_frozen_field_ratio
    logical, intent(in) :: queue_enabled
    type(interface_particle_outcome_type), intent(out) :: outcome
    real(dp), intent(in), optional :: queue_poll_interval
    real(dp) :: speed2, next_speed2, infinity_speed2, outward_time, segment_length, turning_fraction
    real(dp) :: deficit, tolerance, resolved_queue_poll_interval
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
    resolved_queue_poll_interval = 0.0_dp
    if (present(queue_poll_interval)) resolved_queue_poll_interval = queue_poll_interval
    if (queue_enabled .and. &
        (.not. ieee_is_finite(resolved_queue_poll_interval) .or. resolved_queue_poll_interval <= 0.0_dp)) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'persistent kinetic outer queue requires a positive finite poll interval'
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
          queue_enabled, resolved_queue_poll_interval, 2.0_dp*outward_time, outcome &
          )
        return
      end if
      outward_time = outward_time + 2.0_dp*segment_length/(sqrt(speed2) + sqrt(next_speed2))
      speed2 = next_speed2
    end do

    ! The persistent closure owns only the finite [interface,L] control volume.
    ! A queued particle that reaches L leaves that inventory even when the
    ! analytic Robin continuation would turn farther out.
    if (queue_enabled .or. infinity_speed2 >= -tolerance) then
      call finish_kinetic_escape( &
        state, box_min, box_max, charge, mass, crossing, field_timescale, max_frozen_field_ratio, &
        queue_enabled, resolved_queue_poll_interval, outward_time, speed2, outcome &
        )
      return
    end if

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
      queue_enabled, resolved_queue_poll_interval, 2.0_dp*outward_time, outcome &
      )
  end subroutine map_outer_particle_kinetic_profile

  logical function valid_kinetic_profile_state(state) result(valid)
    type(outer_plasma_state_type), intent(in) :: state
    real(dp) :: scale, tolerance
    integer :: minimum_index

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
    select case (trim(lower_ascii(state%kinetic_closure)))
    case ('zhao_charge_driven')
      select case (state%zhao_branch)
      case ('A')
        if (state%profile_n < 3) then
          valid = .false.
          return
        end if
        minimum_index = minloc(state%potential, dim=1)
        valid = minimum_index > 1 .and. minimum_index < state%profile_n
        if (.not. valid) return
        valid = all(state%potential(2:minimum_index) <= &
                    state%potential(1:minimum_index - 1) + tolerance) .and. &
          all(state%potential(minimum_index + 1:) >= &
              state%potential(minimum_index:state%profile_n - 1) - tolerance)
      case ('B')
        valid = all(state%potential(2:) <= state%potential(:state%profile_n - 1) + tolerance)
      case ('C')
        valid = all(state%potential(2:) >= state%potential(:state%profile_n - 1) - tolerance)
      case ('0')
        valid = maxval(abs(state%potential - state%potential(1))) <= tolerance
      case default
        valid = .false.
      end select
    case ('absorbing_maxwellian', 'none', '')
      valid = all(state%potential(2:) <= state%potential(:state%profile_n - 1) + tolerance) .or. &
              all(state%potential(2:) >= state%potential(:state%profile_n - 1) - tolerance)
    case default
      valid = .false.
    end select
  end function valid_kinetic_profile_state

  subroutine finish_kinetic_return( &
    state, box_min, box_max, charge, mass, crossing, field_timescale, max_frozen_field_ratio, &
    queue_enabled, queue_poll_interval, flight_time, outcome &
    )
    type(outer_plasma_state_type), intent(in) :: state
    real(dp), intent(in) :: box_min(3), box_max(3), charge, mass, field_timescale, max_frozen_field_ratio
    real(dp), intent(in) :: queue_poll_interval
    type(interface_crossing_type), intent(in) :: crossing
    logical, intent(in) :: queue_enabled
    real(dp), intent(in) :: flight_time
    type(interface_particle_outcome_type), intent(out) :: outcome
    real(dp) :: inward_offset, energy_before, energy_after, energy_scale

    outcome = interface_particle_outcome_type()
    outcome%outer_flight_time = flight_time
    if (.not. ieee_is_finite(flight_time) .or. flight_time <= 0.0_dp) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'invalid kinetic outer return flight time'
      return
    end if
    outcome%frozen_field_ratio = effective_frozen_duration( &
                                 flight_time, queue_enabled, queue_poll_interval &
                                 )/field_timescale

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
    if (outcome%frozen_field_ratio > max_frozen_field_ratio) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'kinetic outer return violates the frozen-field limit'
    else if (queue_enabled) then
      outcome%kind = interface_outcome_queued_outer
      outcome%queued_terminal_kind = interface_outcome_returned_local
      outcome%message = 'kinetic return scheduled in persistent outer queue'
    else
      outcome%kind = interface_outcome_returned_local
    end if
  end subroutine finish_kinetic_return

  subroutine finish_kinetic_escape( &
    state, box_min, box_max, charge, mass, crossing, field_timescale, max_frozen_field_ratio, &
    queue_enabled, queue_poll_interval, flight_time, final_speed_squared, outcome &
    )
    type(outer_plasma_state_type), intent(in) :: state
    real(dp), intent(in) :: box_min(3), box_max(3), charge, mass, field_timescale, max_frozen_field_ratio
    real(dp), intent(in) :: queue_poll_interval
    type(interface_crossing_type), intent(in) :: crossing
    logical, intent(in) :: queue_enabled
    real(dp), intent(in) :: flight_time, final_speed_squared
    type(interface_particle_outcome_type), intent(out) :: outcome
    real(dp) :: energy_before, energy_after, energy_scale

    outcome = interface_particle_outcome_type()
    if (.not. queue_enabled) then
      outcome%kind = interface_outcome_escaped_to_infinity
      outcome%position = crossing%position
      outcome%velocity = crossing%velocity
      return
    end if
    if (.not. ieee_is_finite(flight_time) .or. flight_time <= 0.0_dp .or. &
        .not. ieee_is_finite(final_speed_squared) .or. final_speed_squared < 0.0_dp) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'invalid kinetic outer escape flight time'
      return
    end if

    outcome%outer_flight_time = flight_time
    outcome%frozen_field_ratio = effective_frozen_duration( &
                                 flight_time, queue_enabled, queue_poll_interval &
                                 )/field_timescale
    outcome%position = crossing%position + crossing%velocity*flight_time
    outcome%position(1) = wrap_periodic(outcome%position(1), box_min(1), box_max(1))
    outcome%position(2) = wrap_periodic(outcome%position(2), box_min(2), box_max(2))
    outcome%position(3) = state%z(state%profile_n)
    outcome%velocity = crossing%velocity
    outcome%velocity(3) = sqrt(final_speed_squared)
    energy_before = 0.5_dp*mass*crossing%velocity(3)**2 + charge*state%interface_potential
    energy_after = 0.5_dp*mass*final_speed_squared + charge*state%potential(state%profile_n)
    outcome%normal_energy_residual = energy_after - energy_before
    energy_scale = max(abs(energy_before), 0.5_dp*mass*crossing%velocity(3)**2, tiny(1.0_dp))
    outcome%energy_relative_error = abs(outcome%normal_energy_residual)/energy_scale
    if (outcome%frozen_field_ratio > max_frozen_field_ratio) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'kinetic outer escape violates the frozen-field limit'
      return
    end if
    outcome%kind = interface_outcome_queued_outer
    outcome%queued_terminal_kind = interface_outcome_escaped_to_infinity
    outcome%message = 'kinetic escape scheduled to finite profile boundary'
  end subroutine finish_kinetic_escape

  !> Bound a midpoint-timestamped event through its crossing-time uncertainty and next queue poll.
  pure real(dp) function effective_frozen_duration(flight_time, queue_enabled, queue_poll_interval) result(duration)
    real(dp), intent(in) :: flight_time, queue_poll_interval
    logical, intent(in) :: queue_enabled
    real(dp) :: due_phase, phase_tolerance, quantization_delay

    duration = flight_time
    if (.not. queue_enabled) return

    due_phase = modulo(0.5_dp*queue_poll_interval + flight_time, queue_poll_interval)
    phase_tolerance = 64.0_dp*epsilon(1.0_dp)*max( &
                      queue_poll_interval, abs(flight_time), tiny(1.0_dp) &
                      )
    if (due_phase <= phase_tolerance .or. queue_poll_interval - due_phase <= phase_tolerance) then
      quantization_delay = 0.0_dp
    else
      quantization_delay = queue_poll_interval - due_phase
    end if
    duration = flight_time + quantization_delay + 0.5_dp*queue_poll_interval
  end function effective_frozen_duration

  pure real(dp) function wrap_periodic(value, low, high) result(wrapped)
    real(dp), intent(in) :: value, low, high

    wrapped = low + modulo(value - low, high - low)
  end function wrap_periodic

end module bem_outer_plasma_interface

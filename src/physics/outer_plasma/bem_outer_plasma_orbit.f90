!> 固定された静電場中で外部領域の粒子軌道を明示的に追跡する。
module bem_outer_plasma_orbit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_config
  use bem_physics_config_types, only: outer_plasma_config, coupling_config
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_interface_types, only: interface_crossing_type, interface_particle_outcome_type, &
                                 interface_outcome_returned_local, interface_outcome_escaped_to_infinity, &
                                 interface_outcome_invalid_model
  implicit none
  private

  type, public :: explicit_outer_orbit_options_type
    real(dp) :: box_min(3) = 0.0_dp
    real(dp) :: box_max(3) = 0.0_dp
    real(dp) :: interface_z = 0.0_dp
    real(dp) :: far_z = 0.0_dp
    real(dp) :: dt = 0.0_dp
    integer(i32) :: max_steps = 0_i32
    real(dp) :: field_evolution_timescale = 0.0_dp
    real(dp) :: max_frozen_field_ratio = 0.0_dp
    real(dp) :: energy_scale = 0.0_dp
    real(dp) :: max_energy_relative_error = 0.0_dp
  end type explicit_outer_orbit_options_type

  abstract interface
    subroutine electrostatic_field_evaluator(position, potential, electric_field)
      import dp
      real(dp), intent(in) :: position(3)
      real(dp), intent(out) :: potential, electric_field(3)
    end subroutine electrostatic_field_evaluator
  end interface

  public :: trace_electrostatic_outer_orbit
  public :: trace_unified_outer_particle

contains

  subroutine trace_unified_outer_particle(snapshot, mesh, sim, outer, coupling, charge, mass, crossing, outcome)
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(outer_plasma_config), intent(in) :: outer
    type(coupling_config), intent(in) :: coupling
    real(dp), intent(in) :: charge, mass
    type(interface_crossing_type), intent(in) :: crossing
    type(interface_particle_outcome_type), intent(out) :: outcome
    type(explicit_outer_orbit_options_type) :: options

    if (.not. snapshot%use_unified_outer .or. snapshot%unified_grid%n < 2_i32) then
      outcome = interface_particle_outcome_type()
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'explicit 3D outer orbit requires a unified electrostatic snapshot'
      return
    end if
    options%box_min = sim%box_min
    options%box_max = sim%box_max
    options%interface_z = outer%interface_z
    options%far_z = snapshot%unified_grid%z(snapshot%unified_grid%n)
    options%dt = coupling%outer_orbit_dt
    options%max_steps = coupling%outer_orbit_max_steps
    options%field_evolution_timescale = coupling%field_evolution_timescale
    options%max_frozen_field_ratio = coupling%max_frozen_field_ratio
    options%energy_scale = max(abs(charge)*outer%thermal_voltage, tiny(1.0_dp))
    options%max_energy_relative_error = coupling%outer_orbit_energy_tolerance
    call trace_electrostatic_outer_orbit(evaluate_snapshot, options, charge, mass, crossing, outcome)

  contains

    subroutine evaluate_snapshot(position, potential, electric_field)
      real(dp), intent(in) :: position(3)
      real(dp), intent(out) :: potential, electric_field(3)

      call snapshot%eval_local_phi(mesh, sim, position, potential)
      call snapshot%eval_local_e(mesh, position, electric_field)
    end subroutine evaluate_snapshot

  end subroutine trace_unified_outer_particle

  subroutine trace_electrostatic_outer_orbit(evaluate, options, charge, mass, crossing, outcome)
    procedure(electrostatic_field_evaluator) :: evaluate
    type(explicit_outer_orbit_options_type), intent(in) :: options
    real(dp), intent(in) :: charge, mass
    type(interface_crossing_type), intent(in) :: crossing
    type(interface_particle_outcome_type), intent(out) :: outcome
    real(dp) :: position(3), velocity(3), next_position(3), next_velocity(3), half_velocity(3)
    real(dp) :: electric_field(3), next_electric_field(3), potential, next_potential
    real(dp) :: acceleration(3), fraction, event_position(3), event_velocity(3)
    real(dp) :: initial_energy, event_energy, event_tolerance, inward_offset
    integer(i32) :: step

    outcome = interface_particle_outcome_type()
    if (.not. valid_inputs(options, charge, mass, crossing)) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'invalid explicit outer-orbit inputs'
      return
    end if

    position = crossing%position
    velocity = crossing%velocity
    call evaluate(position, potential, electric_field)
    if (.not. all(ieee_is_finite([potential, electric_field]))) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'non-finite outer field at interface'
      return
    end if
    initial_energy = 0.5_dp*mass*sum(velocity*velocity) + charge*potential
    event_tolerance = 64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(options%interface_z), abs(options%far_z))

    do step = 1_i32, options%max_steps
      acceleration = (charge/mass)*electric_field
      half_velocity = velocity + 0.5_dp*options%dt*acceleration
      next_position = position + options%dt*half_velocity
      call evaluate(next_position, next_potential, next_electric_field)
      if (.not. all(ieee_is_finite([next_position, next_potential, next_electric_field]))) then
        outcome%kind = interface_outcome_invalid_model
        outcome%message = 'non-finite explicit outer orbit'
        return
      end if
      next_velocity = half_velocity + 0.5_dp*options%dt*(charge/mass)*next_electric_field
      if (.not. all(ieee_is_finite(next_velocity))) then
        outcome%kind = interface_outcome_invalid_model
        outcome%message = 'non-finite explicit outer velocity'
        return
      end if

      if (next_position(3) <= options%interface_z + event_tolerance .and. next_velocity(3) < 0.0_dp) then
        fraction = crossing_fraction(position(3), next_position(3), options%interface_z, event_tolerance)
        event_position = position + fraction*(next_position - position)
        event_velocity = velocity + fraction*(next_velocity - velocity)
        call finish_event( &
          evaluate, options, charge, mass, initial_energy, step, fraction, event_position, event_velocity, &
          interface_outcome_returned_local, outcome &
          )
        if (outcome%kind /= interface_outcome_returned_local) return
        outcome%position(1) = wrap_periodic(outcome%position(1), options%box_min(1), options%box_max(1))
        outcome%position(2) = wrap_periodic(outcome%position(2), options%box_min(2), options%box_max(2))
        inward_offset = 64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(options%interface_z))
        outcome%position(3) = options%interface_z - inward_offset
        return
      end if

      if (next_position(3) >= options%far_z - event_tolerance .and. next_velocity(3) > 0.0_dp) then
        fraction = crossing_fraction(position(3), next_position(3), options%far_z, event_tolerance)
        event_position = position + fraction*(next_position - position)
        event_velocity = velocity + fraction*(next_velocity - velocity)
        call finish_event( &
          evaluate, options, charge, mass, initial_energy, step, fraction, event_position, event_velocity, &
          interface_outcome_escaped_to_infinity, outcome &
          )
        return
      end if

      position = next_position
      velocity = next_velocity
      potential = next_potential
      electric_field = next_electric_field
    end do

    outcome%kind = interface_outcome_invalid_model
    outcome%message = 'explicit outer orbit reached max_steps; persistent queue is required'
  end subroutine trace_electrostatic_outer_orbit

  subroutine finish_event( &
    evaluate, options, charge, mass, initial_energy, step, fraction, position, velocity, kind, outcome &
    )
    procedure(electrostatic_field_evaluator) :: evaluate
    type(explicit_outer_orbit_options_type), intent(in) :: options
    real(dp), intent(in) :: charge, mass, initial_energy, fraction, position(3), velocity(3)
    integer(i32), intent(in) :: step, kind
    type(interface_particle_outcome_type), intent(out) :: outcome
    real(dp) :: potential, electric_field(3), event_energy, scale

    call evaluate(position, potential, electric_field)
    if (.not. all(ieee_is_finite([potential, electric_field, position, velocity]))) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'non-finite explicit outer event'
      return
    end if
    event_energy = 0.5_dp*mass*sum(velocity*velocity) + charge*potential
    scale = max(abs(initial_energy), options%energy_scale, tiny(1.0_dp))
    outcome%energy_relative_error = abs(event_energy - initial_energy)/scale
    if (outcome%energy_relative_error > options%max_energy_relative_error) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'explicit outer orbit exceeds the energy error contract'
      return
    end if
    outcome%outer_flight_time = (real(step - 1_i32, dp) + fraction)*options%dt
    outcome%frozen_field_ratio = outcome%outer_flight_time/options%field_evolution_timescale
    if (outcome%frozen_field_ratio > options%max_frozen_field_ratio) then
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'explicit outer orbit requires a persistent queue'
      return
    end if
    outcome%kind = kind
    outcome%position = position
    outcome%velocity = velocity
  end subroutine finish_event

  logical function valid_inputs(options, charge, mass, crossing) result(valid)
    type(explicit_outer_orbit_options_type), intent(in) :: options
    real(dp), intent(in) :: charge, mass
    type(interface_crossing_type), intent(in) :: crossing

    valid = crossing%has_crossing .and. crossing%face_index == 6_i32 .and. crossing%velocity(3) > 0.0_dp .and. &
            mass > 0.0_dp .and. options%dt > 0.0_dp .and. options%max_steps > 0_i32 .and. &
            options%far_z > options%interface_z .and. all(options%box_max(1:2) > options%box_min(1:2)) .and. &
            options%field_evolution_timescale > 0.0_dp .and. options%max_frozen_field_ratio > 0.0_dp .and. &
            options%energy_scale > 0.0_dp .and. options%max_energy_relative_error > 0.0_dp .and. &
            all(ieee_is_finite([options%box_min, options%box_max, options%interface_z, options%far_z, options%dt, &
                                options%field_evolution_timescale, options%max_frozen_field_ratio, &
                                options%energy_scale, options%max_energy_relative_error, charge, mass, &
                                crossing%position, crossing%velocity]))
  end function valid_inputs

  pure real(dp) function crossing_fraction(z0, z1, plane, tolerance) result(fraction)
    real(dp), intent(in) :: z0, z1, plane, tolerance

    if (abs(z1 - plane) <= tolerance) then
      fraction = 1.0_dp
    else
      fraction = max(0.0_dp, min(1.0_dp, (plane - z0)/(z1 - z0)))
    end if
  end function crossing_fraction

  pure real(dp) function wrap_periodic(value, low, high) result(wrapped)
    real(dp), intent(in) :: value, low, high

    wrapped = low + modulo(value - low, high - low)
  end function wrap_periodic

end module bem_outer_plasma_orbit

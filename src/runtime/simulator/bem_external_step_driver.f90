!> local stepと外部輸送を、同一stepの残り時間が尽きるまで有限回接続する。
module bem_external_step_driver
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_config
  use bem_physics_config_types, only: coupling_config
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_particle_stepper, only: &
    particle_step_result, particle_step_ok, particle_step_invalid_external_model, particle_step_multiple_external_events, &
    advance_particle_step
  use bem_interface_types, only: &
    interface_crossing_type, interface_particle_outcome_type, interface_outcome_returned_local, &
    interface_outcome_escaped_to_infinity, interface_outcome_queued_outer, interface_outcome_invalid_model
  use bem_external_boundary_contract, only: &
    external_boundary_contract_type, external_transport_kinetic_1d
  use bem_outer_plasma_interface, only: map_outer_particle_kinetic_profile
  implicit none
  private

  integer(i32), parameter, public :: max_external_events_per_step = 8_i32

  type, public :: external_step_trace_type
    integer(i32) :: count = 0_i32
    type(interface_crossing_type) :: crossing(max_external_events_per_step)
    type(interface_particle_outcome_type) :: outcome(max_external_events_per_step)
  end type external_step_trace_type

  public :: continue_external_particle_step
  public :: external_trace_ends_at_infinity_escape

contains

  !> True only when the terminal open-boundary result came from the outer infinity map.
  pure function external_trace_ends_at_infinity_escape(trace) result(escaped_to_infinity)
    type(external_step_trace_type), intent(in) :: trace
    logical :: escaped_to_infinity

    escaped_to_infinity = trace%count > 0_i32
    if (escaped_to_infinity) then
      escaped_to_infinity = &
        trace%outcome(trace%count)%kind == interface_outcome_escaped_to_infinity
    end if
  end function external_trace_ends_at_infinity_escape

  subroutine continue_external_particle_step( &
    contract, snapshot, mesh, sim, coupling, bfield, charge, mass, batch_duration, initial_result, &
    final_result, trace, enforce_frozen_field_limit &
    )
    type(external_boundary_contract_type), intent(in) :: contract
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(coupling_config), intent(in) :: coupling
    real(dp), intent(in) :: bfield(3), charge, mass, batch_duration
    type(particle_step_result), intent(in) :: initial_result
    type(particle_step_result), intent(out) :: final_result
    type(external_step_trace_type), intent(out) :: trace
    logical, intent(in), optional :: enforce_frozen_field_limit

    type(particle_step_result) :: current_result, remainder_result
    type(interface_particle_outcome_type) :: outcome
    real(dp) :: dt_remaining, consumed_fraction
    integer(i32) :: event_index
    logical :: resolved_enforce_frozen_field_limit

    trace = external_step_trace_type()
    resolved_enforce_frozen_field_limit = .true.
    if (present(enforce_frozen_field_limit)) then
      resolved_enforce_frozen_field_limit = enforce_frozen_field_limit
    end if
    if (contract%queue_enabled) resolved_enforce_frozen_field_limit = .true.
    current_result = initial_result
    final_result = initial_result
    consumed_fraction = 0.0_dp
    do while (current_result%interface_crossing%has_crossing)
      if (trace%count >= max_external_events_per_step) then
        final_result = initial_result
        final_result%status = particle_step_multiple_external_events
        return
      end if
      event_index = trace%count + 1_i32
      trace%crossing(event_index) = current_result%interface_crossing
      if (event_index == 1_i32) then
        consumed_fraction = current_result%interface_crossing%fraction
      else
        consumed_fraction = consumed_fraction + &
                            (1.0_dp - consumed_fraction)*current_result%interface_crossing%fraction
        trace%crossing(event_index)%fraction = consumed_fraction
      end if
      call dispatch_external_interface_particle( &
        contract, snapshot, sim, coupling, charge, mass, current_result%interface_crossing, &
        batch_duration, resolved_enforce_frozen_field_limit, outcome &
        )
      trace%outcome(event_index) = outcome
      trace%count = event_index

      if (contract%queue_enabled) then
        select case (outcome%kind)
        case (interface_outcome_queued_outer)
          final_result = current_result
          final_result%interface_crossing%has_crossing = .false.
        case default
          final_result = current_result
          final_result%status = particle_step_invalid_external_model
        end select
        return
      end if

      select case (outcome%kind)
      case (interface_outcome_returned_local)
        dt_remaining = current_result%interface_crossing%dt_remaining
        if (dt_remaining > 0.0_dp) then
          call advance_particle_step( &
            mesh, sim, snapshot, bfield, outcome%position, outcome%velocity, charge, mass, dt_remaining, &
            remainder_result, boundary_contract=contract &
            )
          current_result = remainder_result
        else
          current_result%x = outcome%position
          current_result%v = outcome%velocity
          current_result%interface_crossing%has_crossing = .false.
        end if
      case (interface_outcome_escaped_to_infinity)
        current_result%x = outcome%position
        current_result%v = outcome%velocity
        current_result%escaped_boundary = .true.
        current_result%interface_crossing%has_crossing = .false.
      case default
        current_result%status = particle_step_invalid_external_model
        current_result%interface_crossing%has_crossing = .false.
      end select
      if (current_result%status /= particle_step_ok .or. current_result%absorbed .or. &
          current_result%escaped_boundary) exit
    end do
    final_result = current_result
  end subroutine continue_external_particle_step

  !> 解決済みtransport tagを既存の外部粒子modelへ配送する。
  subroutine dispatch_external_interface_particle( &
    contract, snapshot, sim, coupling, charge, mass, crossing, batch_duration, enforce_frozen_field_limit, outcome &
    )
    type(external_boundary_contract_type), intent(in) :: contract
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    type(sim_config), intent(in) :: sim
    type(coupling_config), intent(in) :: coupling
    real(dp), intent(in) :: charge, mass
    type(interface_crossing_type), intent(in) :: crossing
    real(dp), intent(in) :: batch_duration
    logical, intent(in) :: enforce_frozen_field_limit
    type(interface_particle_outcome_type), intent(out) :: outcome

    select case (contract%interface_transport)
    case (external_transport_kinetic_1d)
      call map_outer_particle_kinetic_profile( &
        snapshot%outer, sim%box_min, sim%box_max, charge, mass, crossing, &
        coupling%field_evolution_timescale, coupling%max_frozen_field_ratio, contract%queue_enabled, outcome, &
        queue_poll_interval=batch_duration, enforce_frozen_field_limit=enforce_frozen_field_limit &
        )
    case default
      outcome = interface_particle_outcome_type()
      outcome%kind = interface_outcome_invalid_model
      outcome%message = 'external boundary contract has no interface transport'
    end select
  end subroutine dispatch_external_interface_particle

end module bem_external_step_driver

!> outer return後の残りstepが再びinterfaceを越える経路を検証する。
program test_external_step_driver
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_config, bc_open, bc_periodic
  use bem_mesh, only: init_mesh
  use bem_panel_surface_sides, only: resolve_panel_surface_sides, panel_surface_side_ok
  use bem_physics_config_types, only: outer_plasma_config, coupling_config
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_particle_stepper, only: &
    advance_particle_step, particle_step_result, particle_step_ok, particle_step_multiple_external_events, &
    particle_step_invalid_external_model
  use bem_interface_types, only: &
    interface_outcome_returned_local, interface_outcome_escaped_to_infinity, interface_outcome_invalid_model
  use bem_external_boundary_contract, only: &
    external_boundary_contract_type, external_boundary_ok, resolve_external_boundary_contract
  use bem_external_step_driver, only: &
    external_step_trace_type, continue_external_particle_step, external_trace_ends_at_infinity_escape
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32, assert_close_dp
  implicit none

  type(mesh_type) :: mesh
  type(sim_config) :: sim
  type(outer_plasma_config) :: outer
  type(coupling_config) :: coupling
  type(electrostatic_snapshot_type) :: snapshot
  type(external_boundary_contract_type) :: contract
  type(particle_step_result) :: initial_result, final_result
  type(external_step_trace_type) :: trace
  real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1)
  integer(i32) :: status
  character(len=256) :: message

  call test_init(5)
  call test_begin('instant_return_remainder_can_reenter_transport')

  v0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
  v1(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
  v2(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
  call init_mesh(mesh, v0, v1, v2, q0=[0.0_dp])
  call resolve_panel_surface_sides(mesh, 'normal_plus', status, message)
  call assert_equal_i32(status, panel_surface_side_ok, 'panel side fixture is invalid')

  sim = sim_config()
  sim%field_solver = 'direct'
  sim%field_normalization = 'si'
  sim%use_box = .true.
  sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
  sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
  sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  sim%bc_high = [bc_periodic, bc_periodic, bc_open]
  sim%e0 = [0.0_dp, 0.0_dp, 4.0_dp]

  outer = outer_plasma_config()
  outer%model = 'kinetic_1d'
  outer%return_model = 'kinetic_1d_profile_return'
  outer%interface_z = 1.0_dp
  outer%debye_length = 0.1_dp
  outer%thermal_voltage = 10.0_dp
  coupling = coupling_config()
  coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
  coupling%field_evolution_timescale = 1.0e9_dp
  coupling%max_frozen_field_ratio = 0.1_dp
  call resolve_external_boundary_contract( &
    sim%reservoir_potential_model, sim%open_boundary_model, outer%model, &
    outer%kinetic_closure, outer%return_model, coupling%particle_transfer_mode, coupling%outer_queue_enabled, &
    contract, status, message &
    )
  call assert_equal_i32(status, external_boundary_ok, 'kinetic contract fixture is invalid')

  call snapshot%init(mesh, sim)
  call snapshot%refresh(mesh)
  snapshot%outer%model = 'kinetic_1d'
  snapshot%outer%ready = .true.
  snapshot%outer%interface_z = 1.0_dp
  snapshot%outer%interface_potential = -1.0_dp
  snapshot%outer%infinity_potential = 0.0_dp
  snapshot%outer%debye_length = 0.1_dp
  snapshot%outer%profile_n = 2_i32
  snapshot%outer%z = [1.0_dp, 1.1_dp]
  snapshot%outer%potential = [-1.0_dp, 0.0_dp]

  call advance_particle_step( &
    mesh, sim, snapshot, [0.0_dp, 0.0_dp, 0.0_dp], [0.2_dp, 0.2_dp, 0.9_dp], &
    [0.0_dp, 0.0_dp, 1.0_dp], 1.0_dp, 1.0_dp, 1.0_dp, initial_result, boundary_contract=contract &
    )
  call assert_true(initial_result%interface_crossing%has_crossing, 'initial z-high crossing is missing')

  call continue_external_particle_step( &
    contract, snapshot, mesh, sim, coupling, [0.0_dp, 0.0_dp, 0.0_dp], 1.0_dp, 1.0_dp, 1.0_dp, &
    initial_result, final_result, trace &
    )
  call assert_equal_i32(trace%count, 2_i32, 'remainder must produce exactly two external returns')
  call assert_equal_i32(trace%outcome(1)%kind, interface_outcome_returned_local, 'first outcome must return')
  call assert_equal_i32(trace%outcome(2)%kind, interface_outcome_returned_local, 'second outcome must return')
  call assert_true( &
    trace%crossing(2)%fraction > trace%crossing(1)%fraction .and. trace%crossing(2)%fraction <= 1.0_dp, &
    'trace crossing fractions must use the original local-step time base' &
    )
  call assert_close_dp( &
    trace%crossing(2)%fraction, 1.0_dp - trace%crossing(2)%dt_remaining, 1.0e-12_dp, &
    'second crossing fraction must be normalized by the original local-step duration' &
    )
  call assert_equal_i32(final_result%status, particle_step_ok, 'continued external step failed')
  call assert_true(.not. final_result%interface_crossing%has_crossing, 'continued step retained an interface event')
  call assert_true(.not. final_result%escaped_boundary, 'continued particle escaped unexpectedly')

  call test_end()

  call test_begin('ordinary_outer_model_failure_has_dedicated_status')
  coupling%field_evolution_timescale = 1.0e-12_dp
  call continue_external_particle_step( &
    contract, snapshot, mesh, sim, coupling, [0.0_dp, 0.0_dp, 0.0_dp], 1.0_dp, 1.0_dp, 1.0_dp, &
    initial_result, final_result, trace &
    )
  call assert_equal_i32(trace%count, 1_i32, 'invalid outer mapping must retain its attempted event')
  call assert_equal_i32(trace%outcome(1)%kind, interface_outcome_invalid_model, &
                        'strict frozen-field violation must remain an outer-model outcome')
  call assert_equal_i32(final_result%status, particle_step_invalid_external_model, &
                        'outer-model failure must not be reported as invalid boundary geometry')
  call assert_true(index(trace%outcome(1)%message, 'frozen-field limit') > 0, &
                   'outer-model trace lost its concrete failure reason')
  call test_end()

  call test_begin('shadow_outer_return_bypasses_only_frozen_stop')
  initial_result%interface_crossing%dt_remaining = 0.0_dp
  call continue_external_particle_step( &
    contract, snapshot, mesh, sim, coupling, [0.0_dp, 0.0_dp, 0.0_dp], 1.0_dp, 1.0_dp, 1.0_dp, &
    initial_result, final_result, trace, enforce_frozen_field_limit=.false. &
    )
  call assert_equal_i32(trace%count, 1_i32, 'shadow return must retain exactly one outer event')
  call assert_equal_i32(trace%outcome(1)%kind, interface_outcome_returned_local, &
                        'shadow transport must keep the energy-resolved return outcome')
  call assert_equal_i32(final_result%status, particle_step_ok, 'shadow return failed unexpectedly')
  call assert_true(trace%outcome(1)%frozen_field_ratio > coupling%max_frozen_field_ratio, &
                   'shadow return must preserve the frozen-field applicability diagnostic')
  coupling%field_evolution_timescale = 1.0e9_dp
  call test_end()

  call test_begin('local_open_escape_after_return_is_not_infinity_escape')
  sim%e0 = 0.0_dp
  call snapshot%init(mesh, sim)
  call snapshot%refresh(mesh)
  snapshot%outer%model = 'kinetic_1d'
  snapshot%outer%ready = .true.
  snapshot%outer%interface_z = 1.0_dp
  snapshot%outer%interface_potential = -1.0_dp
  snapshot%outer%infinity_potential = 0.0_dp
  snapshot%outer%debye_length = 0.1_dp
  snapshot%outer%profile_n = 2_i32
  snapshot%outer%z = [1.0_dp, 1.1_dp]
  snapshot%outer%potential = [-1.0_dp, 0.0_dp]
  call advance_particle_step( &
    mesh, sim, snapshot, [0.0_dp, 0.0_dp, 0.0_dp], [0.8_dp, 0.8_dp, 0.9_dp], &
    [0.0_dp, 0.0_dp, 1.0_dp], 1.0_dp, 1.0_dp, 1.0_dp, initial_result, boundary_contract=contract &
    )
  call assert_true(initial_result%interface_crossing%has_crossing, 'local-open fixture crossing is missing')
  initial_result%interface_crossing%dt_remaining = 2.0_dp
  call continue_external_particle_step( &
    contract, snapshot, mesh, sim, coupling, [0.0_dp, 0.0_dp, 0.0_dp], 1.0_dp, 1.0_dp, 1.0_dp, &
    initial_result, final_result, trace &
    )
  call assert_equal_i32(trace%count, 1_i32, 'local-open fixture must contain one outer return')
  call assert_equal_i32(trace%outcome(1)%kind, interface_outcome_returned_local, &
                        'local-open fixture must return from the outer region first')
  call assert_true(final_result%escaped_boundary, 'returned fixture did not leave through the local z-low face')
  call assert_true(.not. external_trace_ends_at_infinity_escape(trace), &
                   'local open-face escape was misclassified as outer infinity escape')

  snapshot%outer%interface_potential = 1.0_dp
  snapshot%outer%potential = [1.0_dp, 0.0_dp]
  call continue_external_particle_step( &
    contract, snapshot, mesh, sim, coupling, [0.0_dp, 0.0_dp, 0.0_dp], 1.0_dp, 1.0_dp, 1.0_dp, &
    initial_result, final_result, trace &
    )
  call assert_equal_i32(trace%outcome(1)%kind, interface_outcome_escaped_to_infinity, &
                        'infinity-escape control did not escape the outer profile')
  call assert_true(external_trace_ends_at_infinity_escape(trace), &
                   'outer infinity escape was not recognized')

  sim%e0 = [0.0_dp, 0.0_dp, 4.0_dp]
  call snapshot%init(mesh, sim)
  call snapshot%refresh(mesh)
  snapshot%outer%model = 'kinetic_1d'
  snapshot%outer%ready = .true.
  snapshot%outer%interface_z = 1.0_dp
  snapshot%outer%interface_potential = -1.0_dp
  snapshot%outer%infinity_potential = 0.0_dp
  snapshot%outer%debye_length = 0.1_dp
  snapshot%outer%profile_n = 2_i32
  snapshot%outer%z = [1.0_dp, 1.1_dp]
  snapshot%outer%potential = [-1.0_dp, 0.0_dp]
  call test_end()

  call test_begin('ninth_external_event_fails_with_dedicated_status')
  call advance_particle_step( &
    mesh, sim, snapshot, [0.0_dp, 0.0_dp, 0.0_dp], [0.2_dp, 0.2_dp, 0.9_dp], &
    [0.0_dp, 0.0_dp, 1.0_dp], 1.0_dp, 1.0_dp, 10.0_dp, initial_result, boundary_contract=contract &
    )
  call assert_true(initial_result%interface_crossing%has_crossing, 'cap fixture initial crossing is missing')
  call continue_external_particle_step( &
    contract, snapshot, mesh, sim, coupling, [0.0_dp, 0.0_dp, 0.0_dp], 1.0_dp, 1.0_dp, 1.0_dp, &
    initial_result, final_result, trace &
    )
  call assert_equal_i32(trace%count, 8_i32, 'external cap must retain exactly eight attempted outcomes')
  call assert_equal_i32( &
    final_result%status, particle_step_multiple_external_events, 'ninth external event needs its dedicated status' &
    )
  call assert_close_dp(final_result%x(3), initial_result%x(3), 0.0_dp, 'failed external step must not commit position')
  call assert_close_dp(final_result%v(3), initial_result%v(3), 0.0_dp, 'failed external step must not commit velocity')
  call test_end()

  call test_summary()

end program test_external_step_driver

program test_periodic_zero_mode
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_types, only: mesh_type
  use bem_mesh, only: init_mesh
  use bem_periodic_zero_mode_plan, only: periodic_zero_mode_plan_type, periodic_zero_mode_state_type, &
                                         build_periodic_zero_mode_plan, refresh_periodic_zero_mode_state, &
                                         symmetric_vacuum_bottom_field
  use bem_periodic_zero_mode_eval, only: zero_mode_trace_minus, zero_mode_trace_principal_value, &
                                         zero_mode_trace_plus, eval_periodic_zero_mode
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_close_dp, assert_equal_i32
  implicit none

  type(mesh_type) :: mesh
  type(periodic_zero_mode_plan_type) :: plan
  type(periodic_zero_mode_state_type) :: state
  real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1), field, potential
  real(dp), parameter :: charge = 3.0e-12_dp
  integer(i32) :: status
  character(len=128) :: message

  call test_init(4)

  call test_begin('horizontal_sheet_jump_and_gauge')
  v0(:, 1) = [0.0_dp, 0.0_dp, 0.5_dp]
  v1(:, 1) = [1.0_dp, 0.0_dp, 0.5_dp]
  v2(:, 1) = [0.0_dp, 1.0_dp, 0.5_dp]
  call init_mesh(mesh, v0, v1, v2, q0=[charge])
  call build_periodic_zero_mode_plan(mesh, 2.0_dp, plan, status, message)
  call assert_equal_i32(status, 0_i32, 'horizontal plan status mismatch')
  call refresh_periodic_zero_mode_state(plan, mesh%q_elem, 0.0_dp, 0.0_dp, 0.0_dp, state)
  call eval_periodic_zero_mode(plan, state, 0.5_dp, zero_mode_trace_minus, potential, field)
  call assert_close_dp(field, 0.0_dp, 1.0e-14_dp, 'sheet minus trace mismatch')
  call eval_periodic_zero_mode(plan, state, 0.5_dp, zero_mode_trace_principal_value, potential, field)
  call assert_close_dp(field, charge/(4.0_dp*eps0), 1.0e-12_dp, 'sheet PV trace mismatch')
  call eval_periodic_zero_mode(plan, state, 0.5_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(field, charge/(2.0_dp*eps0), 1.0e-12_dp, 'sheet plus trace mismatch')
  call eval_periodic_zero_mode(plan, state, 0.0_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(potential, 0.0_dp, 1.0e-14_dp, 'gauge potential mismatch')
  call test_end()

  call test_begin('inclined_triangle_exact_cumulative_area')
  v0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
  v1(:, 1) = [1.0_dp, 0.0_dp, 1.0_dp]
  v2(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp]
  call init_mesh(mesh, v0, v1, v2, q0=[charge])
  call build_periodic_zero_mode_plan(mesh, 1.0_dp, plan, status, message)
  call refresh_periodic_zero_mode_state(plan, mesh%q_elem, 0.0_dp, 0.0_dp, 0.0_dp, state)
  call eval_periodic_zero_mode(plan, state, 0.5_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(field, charge/(8.0_dp*eps0), 2.0e-12_dp, 'lower cumulative fraction mismatch')
  call eval_periodic_zero_mode(plan, state, 1.5_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(field, 7.0_dp*charge/(8.0_dp*eps0), 2.0e-12_dp, 'upper cumulative fraction mismatch')
  call test_end()

  call test_begin('potential_integrates_same_field')
  call eval_periodic_zero_mode(plan, state, 0.5_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(potential, -charge/(48.0_dp*eps0), 2.0e-12_dp, 'lower potential integral mismatch')
  call eval_periodic_zero_mode(plan, state, 2.5_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(potential, -1.5_dp*charge/eps0, 2.0e-12_dp, 'upper potential integral mismatch')
  call test_end()

  call test_begin('symmetric_vacuum_splits_far_field_equally')
  call refresh_periodic_zero_mode_state( &
    plan, mesh%q_elem, symmetric_vacuum_bottom_field(plan, mesh%q_elem), 0.0_dp, 0.0_dp, state &
    )
  call eval_periodic_zero_mode(plan, state, -1.0_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(field, -charge/(2.0_dp*eps0), 2.0e-12_dp, 'symmetric lower far field mismatch')
  call eval_periodic_zero_mode(plan, state, 3.0_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(field, charge/(2.0_dp*eps0), 2.0e-12_dp, 'symmetric upper far field mismatch')
  call test_end()

  call test_summary()
end program test_periodic_zero_mode

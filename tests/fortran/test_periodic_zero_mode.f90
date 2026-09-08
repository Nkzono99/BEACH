program test_periodic_zero_mode
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_types, only: mesh_type
  use bem_mesh, only: init_mesh
  use bem_periodic_zero_mode_plan, only: periodic_zero_mode_plan_type, periodic_zero_mode_state_type, &
                                         build_periodic_zero_mode_plan, refresh_periodic_zero_mode_state, &
                                         build_periodic_zero_mode_height_plan, &
                                         symmetric_vacuum_bottom_field, periodic_zero_mode_ok
  use bem_periodic_zero_mode_eval, only: zero_mode_trace_minus, zero_mode_trace_principal_value, &
                                         zero_mode_trace_plus, eval_periodic_zero_mode
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_close_dp, assert_equal_i32
  implicit none

  type(mesh_type) :: mesh
  type(periodic_zero_mode_plan_type) :: plan
  type(periodic_zero_mode_state_type) :: state
  real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2), field, potential
  real(dp), parameter :: charge = 3.0e-12_dp
  integer(i32) :: status
  character(len=128) :: message
  real(dp), parameter :: sheet_z(9) = [-100.0_dp, -4.0_dp, -1.0_dp, 0.0_dp, 0.5_dp, &
                                       2.0_dp, 5.0_dp, 20.0_dp, 60.0_dp]
  real(dp), parameter :: gauges(3) = [-120.0_dp, -0.3_dp, 80.0_dp]
  real(dp) :: heights(3, 9), charges(9), z, tolerance, expected_field, expected_potential
  integer(i32) :: i, j, gauge, offset, trace

  call test_init(3)

  call test_begin('horizontal_sheet_jumps_superposition_and_gauge')
  v0(:, 1) = [0.0_dp, 0.0_dp, 0.5_dp]
  v1(:, 1) = [1.0_dp, 0.0_dp, 0.5_dp]
  v2(:, 1) = [0.0_dp, 1.0_dp, 0.5_dp]
  v0(:, 2) = [1.0_dp, 1.0_dp, 1.5_dp]
  v1(:, 2) = [0.0_dp, 1.0_dp, 1.5_dp]
  v2(:, 2) = [1.0_dp, 0.0_dp, 1.5_dp]
  call init_mesh(mesh, v0, v1, v2, q0=[charge/3.0_dp, 2.0_dp*charge/3.0_dp])
  call build_periodic_zero_mode_plan(mesh, 2.0_dp, plan, status, message)
  call assert_equal_i32(status, periodic_zero_mode_ok, 'horizontal plan status mismatch')
  call refresh_periodic_zero_mode_state(plan, mesh%q_elem, 0.0_dp, 0.0_dp, 1.25_dp, state)
  call eval_periodic_zero_mode(plan, state, 0.5_dp, zero_mode_trace_minus, potential, field)
  call assert_close_dp(field, 0.0_dp, 1.0e-14_dp, 'sheet minus trace mismatch')
  call eval_periodic_zero_mode(plan, state, 0.5_dp, zero_mode_trace_principal_value, potential, field)
  call assert_close_dp(field, charge/(12.0_dp*eps0), 1.0e-12_dp, 'sheet PV trace mismatch')
  call eval_periodic_zero_mode(plan, state, 0.5_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(field, charge/(6.0_dp*eps0), 1.0e-12_dp, 'sheet plus trace mismatch')
  call eval_periodic_zero_mode(plan, state, 2.0_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(field, charge/(2.0_dp*eps0), 1.0e-12_dp, 'sheet superposition mismatch')
  call eval_periodic_zero_mode(plan, state, 0.0_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(potential, 1.25_dp, 1.0e-14_dp, 'gauge potential mismatch')
  call test_end()

  call test_begin('inclined_triangle_exact_field_potential_and_closure')
  ! Cyclic vertex order deliberately gives unsorted heights [2, 0, 1].
  v0(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp]
  v1(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
  v2(:, 1) = [1.0_dp, 0.0_dp, 1.0_dp]
  call init_mesh(mesh, v0(:, 1:1), v1(:, 1:1), v2(:, 1:1), q0=[charge])
  call build_periodic_zero_mode_plan(mesh, 1.0_dp, plan, status, message)
  call assert_equal_i32(status, periodic_zero_mode_ok, 'inclined plan status mismatch')
  call refresh_periodic_zero_mode_state(plan, mesh%q_elem, 0.0_dp, 0.0_dp, 0.0_dp, state)
  call eval_periodic_zero_mode(plan, state, 0.5_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(field, charge/(8.0_dp*eps0), 2.0e-12_dp, 'lower cumulative fraction mismatch')
  call eval_periodic_zero_mode(plan, state, 1.5_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(field, 7.0_dp*charge/(8.0_dp*eps0), 2.0e-12_dp, 'upper cumulative fraction mismatch')
  call eval_periodic_zero_mode(plan, state, 0.5_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(potential, -charge/(48.0_dp*eps0), 2.0e-12_dp, 'lower potential integral mismatch')
  call eval_periodic_zero_mode(plan, state, 2.5_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(potential, -1.5_dp*charge/eps0, 2.0e-12_dp, 'upper potential integral mismatch')
  call refresh_periodic_zero_mode_state( &
    plan, mesh%q_elem, symmetric_vacuum_bottom_field(plan, mesh%q_elem), 0.0_dp, 0.0_dp, state &
    )
  call eval_periodic_zero_mode(plan, state, -1.0_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(field, -charge/(2.0_dp*eps0), 2.0e-12_dp, 'symmetric lower far field mismatch')
  call eval_periodic_zero_mode(plan, state, 3.0_dp, zero_mode_trace_plus, potential, field)
  call assert_close_dp(field, charge/(2.0_dp*eps0), 2.0e-12_dp, 'symmetric upper far field mismatch')
  call test_end()

  call test_begin('signed_sheet_heights_traces_and_arbitrary_gauge')
  do i = 1, size(sheet_z)
    heights(:, i) = sheet_z(i)
    charges(i) = real((-1)**i, dp)*charge
  end do
  call build_periodic_zero_mode_height_plan(heights, 2.0_dp, plan, status, message)
  call assert_equal_i32(status, periodic_zero_mode_ok, 'signed sheet plan status mismatch')
  tolerance = 128.0_dp*epsilon(1.0_dp)*100.0_dp
  do gauge = 1, size(gauges)
    call refresh_periodic_zero_mode_state(plan, charges, 0.2_dp, gauges(gauge), 1.25_dp, state)
    do i = 1, size(sheet_z)
      do offset = -2, 2
        z = sheet_z(i) + real(offset, dp)*0.75_dp*tolerance
        do trace = zero_mode_trace_minus, zero_mode_trace_plus
          expected_field = 0.2_dp
          expected_potential = 1.25_dp - 0.2_dp*(z - gauges(gauge))
          do j = 1, size(sheet_z)
            if (abs(z - sheet_z(j)) <= tolerance) then
              expected_field = expected_field + 0.5_dp*real(1 + trace, dp)*charges(j)/(2.0_dp*eps0)
            else if (z > sheet_z(j)) then
              expected_field = expected_field + charges(j)/(2.0_dp*eps0)
            end if
            expected_potential = expected_potential - charges(j)/(2.0_dp*eps0)* &
                                 (max(z - sheet_z(j), 0.0_dp) - max(gauges(gauge) - sheet_z(j), 0.0_dp))
          end do
          call eval_periodic_zero_mode(plan, state, z, trace, potential, field)
          call assert_close_dp(field, expected_field, 1.0e-12_dp, 'signed sheet field mismatch')
          call assert_close_dp(potential, expected_potential, 1.0e-10_dp, 'signed sheet potential mismatch')
        end do
      end do
    end do
    call eval_periodic_zero_mode(plan, state, gauges(gauge), zero_mode_trace_plus, potential, field)
    call assert_close_dp(potential, 1.25_dp, 1.0e-14_dp, 'arbitrary gauge potential mismatch')
  end do
  call test_end()

  call test_summary()
end program test_periodic_zero_mode

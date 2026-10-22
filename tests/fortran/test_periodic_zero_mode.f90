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
  real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2), field, potential, reference_phi
  type(periodic_zero_mode_state_type) :: top_zero_state
  real(dp), parameter :: e_top = -2.7_dp, z_high = 3.0_dp
  real(dp), parameter :: charge = 3.0e-12_dp
  integer(i32) :: status
  character(len=128) :: message
  real(dp), parameter :: sheet_z(9) = [-100.0_dp, -4.0_dp, -1.0_dp, 0.0_dp, 0.5_dp, &
                                       2.0_dp, 5.0_dp, 20.0_dp, 60.0_dp]
  real(dp), parameter :: gauges(3) = [-120.0_dp, -0.3_dp, 80.0_dp]
  real(dp) :: heights(3, 9), charges(9), z, tolerance, expected_field, expected_potential
  integer(i32) :: i, j, gauge, offset, trace

  call test_init(6)

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

  call test_begin('charged_layers_with_prescribed_top_field')
  call refresh_periodic_zero_mode_state(plan, mesh%q_elem, -sum(mesh%q_elem)/(eps0*plan%area_xy), &
                                        z_high, 0.0_dp, top_zero_state)
  call refresh_periodic_zero_mode_state(plan, mesh%q_elem, e_top - sum(mesh%q_elem)/(eps0*plan%area_xy), &
                                        z_high, 0.0_dp, state)
  call eval_periodic_zero_mode(plan, state, z_high, zero_mode_trace_plus, potential, field)
  call assert_close_dp(field, e_top, abs(e_top)*1.0e-12_dp, 'Gauss top Neumann field')
  do i = 1, 5
    z = 0.4_dp*real(i, dp)
    call eval_periodic_zero_mode(plan, state, z, zero_mode_trace_plus, potential, field)
    call eval_periodic_zero_mode(plan, top_zero_state, z, zero_mode_trace_plus, reference_phi, field)
    call assert_close_dp(potential - reference_phi, e_top*(z_high - z), 1.0e-12_dp, 'top-field linear potential shift')
  end do
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

  call test_begin('thin_inclined_panels_have_constant_vacuum_field')
  do i = 1, size(charges)
    z = 1.0e-4_dp + real(i, dp)*1.3e-8_dp
    heights(:, i) = z + [0.0_dp, 0.3e-8_dp, 0.9e-8_dp]
    charges(i) = charge*real((-1)**i*i, dp)
  end do
  call build_periodic_zero_mode_height_plan(heights, 1.0e-8_dp, plan, status, message)
  call assert_equal_i32(status, periodic_zero_mode_ok, 'thin-panel plan status')
  expected_field = sum(charges)/(eps0*1.0e-8_dp)
  call refresh_periodic_zero_mode_state(plan, charges, 0.0_dp, 2.0e-4_dp, 0.0_dp, state)
  do i = 1, 3
    z = real(i + 2, dp)*1.0e-4_dp
    call eval_periodic_zero_mode(plan, state, z, zero_mode_trace_plus, potential, field)
    call assert_close_dp(field, expected_field, 8*epsilon(1.0_dp)*abs(expected_field), &
                         'vacuum Gauss field must use total charge')
    expected_potential = -expected_field*(z - 2.0e-4_dp)
    call assert_close_dp(potential, expected_potential, 32*epsilon(1.0_dp)*abs(expected_potential), &
                         'vacuum potential must remain linear')
  end do
  call test_end()

  call test_begin('first_moment_is_translation_invariant')
  call check_first_moment_translation()
  call test_end()

  call test_summary()

contains

  !> 一様面密度の panel を上端より上（下）から見ると、累積電荷の積分は q*(z-mean(頂点高さ)) になる。
  !! 薄い panel 群を z 方向へ平行移動しても、この一次モーメントと Gauss 則が保たれることを確認する。
  subroutine check_first_moment_translation()
    real(dp), parameter :: offsets(11) = [-10.0_dp, -1.0_dp, -0.1_dp, -1.0e-2_dp, -1.0e-3_dp, 0.0_dp, &
                                          1.0e-3_dp, 1.0e-2_dp, 0.1_dp, 1.0_dp, 10.0_dp]
    real(dp), parameter :: local_heights(3, 4) = reshape([ &
                                                         0.0_dp, 1.0e-6_dp, 2.0e-6_dp, &
                                                         2.5e-6_dp, 0.5e-6_dp, 1.2e-6_dp, &
                                                         1.5e-6_dp, 1.5e-6_dp, 1.5e-6_dp, &
                                                         3.0e-6_dp, 0.2e-6_dp, 0.9e-6_dp], [3, 4])
    real(dp), parameter :: panel_charge(4) = [1.0e-15_dp, -0.4e-15_dp, 0.25e-15_dp, 0.6e-15_dp]
    real(dp), parameter :: area = 1.0e-10_dp, e_bottom = 0.15_dp, phi_gauge = 0.3_dp
    real(dp) :: panel_heights(3, 4), z_below, z_above, z_inside, z_gauge, moment, field_above
    real(dp) :: expected, reference_inside(2), phi_inside
    integer(i32) :: k, j, gauge_case

    do k = 1, size(offsets)
      panel_heights = offsets(k) + local_heights
      call build_periodic_zero_mode_height_plan(panel_heights, area, plan, status, message)
      call assert_equal_i32(status, periodic_zero_mode_ok, 'translated plan status')
      z_below = offsets(k) - 1.0e-6_dp
      z_above = offsets(k) + 5.0e-6_dp
      z_inside = offsets(k) + 1.1e-6_dp
      field_above = e_bottom + sum(panel_charge)/(eps0*area)
      do gauge_case = 1, 2
        ! 1: 下側 gauge で上端より上を、2: 上側 gauge で下端より下を評価する。
        if (gauge_case == 1) then
          z_gauge = z_below
          z = z_above
        else
          z_gauge = z_above
          z = z_below
        end if
        ! sum_j q_j*(z_above - mean_j)。差は同じ binary64 座標の近い値どうしなので丸めない。
        moment = 0.0_dp
        do j = 1, size(panel_charge)
          moment = moment + panel_charge(j)*((z_above - panel_heights(1, j)) - &
                                             ((panel_heights(2, j) - panel_heights(1, j)) + &
                                              (panel_heights(3, j) - panel_heights(1, j)))/3.0_dp)
        end do
        call refresh_periodic_zero_mode_state(plan, panel_charge, e_bottom, z_gauge, phi_gauge, state)
        call eval_periodic_zero_mode(plan, state, z, zero_mode_trace_plus, potential, field)
        if (gauge_case == 1) then
          expected = phi_gauge - e_bottom*(z - z_gauge) - moment/(eps0*area)
          call assert_close_dp(field, field_above, 1.0e-12_dp*abs(field_above), 'translated Gauss field above')
        else
          ! 下端より下では累積電荷が 0 なので、gauge（z_above）からの積分は -moment になる。
          expected = phi_gauge + e_bottom*(z_gauge - z) + moment/(eps0*area)
          call assert_close_dp(field, e_bottom, 1.0e-12_dp*abs(field_above), 'translated Gauss field below')
        end if
        call assert_close_dp(potential, expected, 1.0e-10_dp*abs(moment)/(eps0*area), &
                             'translated first-moment potential')
        call eval_periodic_zero_mode(plan, state, z_gauge, zero_mode_trace_plus, potential, field)
        call assert_close_dp(potential, phi_gauge, 1.0e-12_dp, 'translated gauge potential')
        ! slab 内の値は gauge からの相対位置だけで決まる。原点 0 の値と比べる。
        call eval_periodic_zero_mode(plan, state, z_inside, zero_mode_trace_plus, phi_inside, field)
        if (k == 1) then
          reference_inside(gauge_case) = phi_inside
        else
          call assert_close_dp(phi_inside, reference_inside(gauge_case), &
                               1.0e-8_dp*abs(moment)/(eps0*area), 'translated in-slab potential')
        end if
      end do
    end do
  end subroutine check_first_moment_translation
end program test_periodic_zero_mode

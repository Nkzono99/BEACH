program test_panel_near_correction
  use bem_constants, only: k_coulomb
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, fmm_state_type, build_panel_plan, &
                                  update_state, eval_point, destroy_plan, destroy_state
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_kernel, only: panel_potential_field, panel_side_principal_value
  use test_support, only: test_begin, test_end, test_summary, assert_true, assert_equal_i32
  implicit none

  type(fmm_options_type) :: options
  type(fmm_plan_type) :: plan
  type(fmm_state_type) :: state
  type(panel_geometry_type) :: geometry
  real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1), q(1), target(3)
  real(dp) :: field(3), field_ref(3), field_point(3), potential, delta(3), r2
  real(dp) :: panel_error, point_error
  integer(i32) :: status

  call test_begin('panel_near_replaces_centroid_point_kernel')
  v0(:, 1) = [-4.0_dp, 0.0_dp, 0.0_dp]
  v1(:, 1) = [4.0_dp, 0.0_dp, 0.0_dp]
  v2(:, 1) = [0.0_dp, 0.2_dp, 0.0_dp]
  q = [3.0_dp]
  target = [0.0_dp, 0.05_dp, 0.08_dp]
  options%leaf_max = 1_i32
  options%order = 4_i32
  options%target_box_min = [-5.0_dp, -1.0_dp, -1.0_dp]
  options%target_box_max = [5.0_dp, 1.0_dp, 1.0_dp]
  call build_panel_plan(plan, v0, v1, v2, options)
  call update_state(plan, state, q)
  call eval_point(plan, state, target, field)
  call init_panel_geometry(v0(:, 1), v1(:, 1), v2(:, 1), geometry, status)
  call assert_equal_i32(status, panel_geometry_ok, 'panel fixture geometry status')
  call panel_potential_field(geometry, q(1), target, panel_side_principal_value, potential, field_ref)
  field_ref = field_ref/k_coulomb
  delta = target - geometry%centroid
  r2 = sum(delta*delta)
  field_point = q(1)*delta/(sqrt(r2)*r2)
  panel_error = sqrt(sum((field - field_ref)**2))
  point_error = sqrt(sum((field_point - field_ref)**2))
  call assert_true(panel_error <= 1.0e-12_dp*max(1.0_dp, sqrt(sum(field_ref*field_ref))), &
                   'panel near correction must match the analytic kernel')
  call assert_true(point_error > 1.0e-2_dp*sqrt(sum(field_ref*field_ref)), &
                   'fixture must expose a measurable centroid-point error')
  call destroy_state(state)
  call destroy_plan(plan)
  call test_end()
  call test_summary()
end program test_panel_near_correction

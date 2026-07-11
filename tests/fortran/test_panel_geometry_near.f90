program test_panel_geometry_near
  use bem_constants, only: k_coulomb
  use bem_kinds, only: dp, i32
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok, point_triangle_distance
  use bem_panel_kernel, only: panel_potential_field, panel_side_principal_value
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, fmm_state_type, build_panel_plan, update_state, &
                                  eval_point, eval_potential_point, destroy_plan, destroy_state
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp
  implicit none

  type(fmm_options_type) :: options
  type(fmm_plan_type) :: plan
  type(fmm_state_type) :: state
  type(panel_geometry_type) :: geometry
  real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2), expected, center(3), mean_x2
  real(dp) :: target(3), q(1), phi, phi_ref, field(3), field_ref(3)
  integer(i32) :: status, alpha_x, alpha_xx

  call test_init(5)

  call test_begin('point_to_triangle_distance')
  call init_panel_geometry( &
    [0.0_dp, 0.0_dp, 0.0_dp], [2.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 2.0_dp, 0.0_dp], geometry, status &
    )
  call assert_true(status == panel_geometry_ok, 'triangle fixture should be valid')
  call assert_close_dp(point_triangle_distance(geometry, [0.5_dp, 0.5_dp, 3.0_dp]), 3.0_dp, 1.0e-14_dp, &
                       'normal projection distance mismatch')
  expected = sqrt(2.0_dp)
  call assert_close_dp(point_triangle_distance(geometry, [3.0_dp, 1.0_dp, 0.0_dp]), expected, 1.0e-14_dp, &
                       'edge distance mismatch')
  call test_end()

  call test_begin('panel_tree_bbox_contains_all_vertices')
  v0(:, 1) = [-10.0_dp, 0.0_dp, 0.0_dp]
  v1(:, 1) = [10.0_dp, 0.0_dp, 0.0_dp]
  v2(:, 1) = [0.0_dp, 0.01_dp, 0.0_dp]
  v0(:, 2) = [0.0_dp, 2.0_dp, 0.0_dp]
  v1(:, 2) = [0.1_dp, 2.0_dp, 0.0_dp]
  v2(:, 2) = [0.0_dp, 2.1_dp, 0.0_dp]
  options%leaf_max = 1_i32
  options%order = 2_i32
  call build_panel_plan(plan, v0, v1, v2, options)
  call assert_true(all(v0 >= spread(plan%node_center(:, 1) - plan%node_half_size(:, 1), 2, 2)), &
                   'root lower bbox must contain v0')
  call assert_true(all(v1 <= spread(plan%node_center(:, 1) + plan%node_half_size(:, 1), 2, 2)), &
                   'root upper bbox must contain v1')
  call assert_true(plan%node_radius(1) >= 10.0_dp, 'root radius must include elongated panel vertices')
  call test_end()

  call test_begin('exact_panel_p2m_moments')
  call destroy_plan(plan)
  v0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
  v1(:, 1) = [2.0_dp, 0.0_dp, 0.0_dp]
  v2(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
  options%leaf_max = 8_i32
  options%order = 2_i32
  call build_panel_plan(plan, v0(:, :1), v1(:, :1), v2(:, :1), options)
  call init_panel_geometry(v0(:, 1), v1(:, 1), v2(:, 1), geometry, status)
  center = plan%node_center(:, 1)
  alpha_x = plan%alpha_map(1, 0, 0)
  alpha_xx = plan%alpha_map(2, 0, 0)
  mean_x2 = geometry%moment2(1, 1)/geometry%area - 2.0_dp*center(1)*geometry%centroid(1) + center(1)**2
  call assert_close_dp(plan%source_p2m_basis(1, 1), 1.0_dp, 1.0e-14_dp, 'panel monopole basis mismatch')
  call assert_close_dp(plan%source_p2m_basis(alpha_x, 1), geometry%centroid(1) - center(1), 1.0e-14_dp, &
                       'panel dipole basis mismatch')
  call assert_close_dp(plan%source_p2m_basis(alpha_xx, 1), 0.5_dp*mean_x2, 1.0e-14_dp, &
                       'panel quadrupole basis mismatch')
  call test_end()

  call test_begin('panel_near_uses_analytic_kernel')
  call destroy_plan(plan)
  options%target_box_min = [-1.0_dp, -1.0_dp, -1.0_dp]
  options%target_box_max = [1.0_dp, 1.0_dp, 1.0_dp]
  call build_panel_plan(plan, v0(:, :1), v1(:, :1), v2(:, :1), options)
  q = [2.5_dp]
  call update_state(plan, state, q)
  target = [0.25_dp, 0.25_dp, 0.5_dp]
  call eval_point(plan, state, target, field)
  call eval_potential_point(plan, state, target, phi)
  call panel_potential_field(geometry, q(1), target, panel_side_principal_value, phi_ref, field_ref)
  call assert_true(all(abs(field - field_ref/k_coulomb) < 1.0e-13_dp), 'near field must use analytic panel kernel')
  call assert_close_dp(phi, phi_ref/k_coulomb, 1.0e-13_dp, 'near potential must use analytic panel kernel')
  call test_end()

  call test_begin('panel_outside_box_uses_analytic_fallback')
  target = [2.0_dp, 1.5_dp, 0.75_dp]
  call eval_point(plan, state, target, field)
  call eval_potential_point(plan, state, target, phi)
  call panel_potential_field(geometry, q(1), target, panel_side_principal_value, phi_ref, field_ref)
  call assert_true(all(abs(field - field_ref/k_coulomb) < 1.0e-13_dp), 'fallback field must use analytic panel kernel')
  call assert_close_dp(phi, phi_ref/k_coulomb, 1.0e-13_dp, 'fallback potential must use analytic panel kernel')
  call test_end()

  call destroy_state(state)
  call destroy_plan(plan)
  call test_summary()
end program test_panel_geometry_near

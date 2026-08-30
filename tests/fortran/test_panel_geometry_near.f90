program test_panel_geometry_near
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use bem_constants, only: k_coulomb
  use bem_kinds, only: dp, i32
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok, panel_geometry_nonfinite, &
                                panel_geometry_degenerate, point_triangle_distance
  use bem_panel_kernel, only: panel_potential_field, panel_side_principal_value
  use bem_coulomb_fmm_core, only: fmm_options_type, fmm_plan_type, fmm_state_type, build_panel_plan, update_state, &
                                  eval_point, eval_potential_point, destroy_plan, destroy_state
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, assert_equal_i32
  implicit none

  type(fmm_options_type) :: options
  type(fmm_plan_type) :: plan
  type(fmm_state_type) :: state
  type(panel_geometry_type) :: geometry, transformed_geometry
  real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2), vertices(3, 6), center(3)
  real(dp) :: target(3), q(1), phi, phi_ref, field(3), field_ref(3)
  real(dp) :: shift(3), scale, max_vertex_radius, nan_value
  integer(i32) :: status, alpha_zero, alpha_x, alpha_xx
  real(dp), parameter :: triangle_mean_x = 2.0_dp/3.0_dp
  real(dp), parameter :: triangle_mean_x2 = 2.0_dp/3.0_dp

  call test_init(5)

  call test_begin('point_to_triangle_distance_regions_and_invariants')
  call init_panel_geometry( &
    [0.0_dp, 0.0_dp, 0.0_dp], [2.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 2.0_dp, 0.0_dp], geometry, status &
    )
  call assert_equal_i32(status, panel_geometry_ok, 'triangle fixture should be valid')
  call assert_close_dp(point_triangle_distance(geometry, [0.5_dp, 0.5_dp, 3.0_dp]), 3.0_dp, 1.0e-14_dp, &
                       'face-region distance mismatch')
  call assert_close_dp(point_triangle_distance(geometry, [0.5_dp, 0.5_dp, 0.0_dp]), 0.0_dp, 1.0e-14_dp, &
                       'coplanar interior distance mismatch')
  call assert_close_dp(point_triangle_distance(geometry, [1.0_dp, -1.0_dp, 0.0_dp]), 1.0_dp, 1.0e-14_dp, &
                       'v0-v1 edge-region distance mismatch')
  call assert_close_dp(point_triangle_distance(geometry, [-1.0_dp, 1.0_dp, 0.0_dp]), 1.0_dp, 1.0e-14_dp, &
                       'v0-v2 edge-region distance mismatch')
  call assert_close_dp(point_triangle_distance(geometry, [1.5_dp, 1.5_dp, 0.0_dp]), sqrt(0.5_dp), 1.0e-14_dp, &
                       'v1-v2 edge-region distance mismatch')
  call assert_close_dp(point_triangle_distance(geometry, [-1.0_dp, -1.0_dp, 0.0_dp]), sqrt(2.0_dp), 1.0e-14_dp, &
                       'v0 vertex-region distance mismatch')
  call assert_close_dp(point_triangle_distance(geometry, [3.0_dp, -1.0_dp, 0.0_dp]), sqrt(2.0_dp), 1.0e-14_dp, &
                       'v1 vertex-region distance mismatch')
  call assert_close_dp(point_triangle_distance(geometry, [-1.0_dp, 3.0_dp, 0.0_dp]), sqrt(2.0_dp), 1.0e-14_dp, &
                       'v2 vertex-region distance mismatch')

  shift = [7.0_dp, -4.0_dp, 2.0_dp]
  scale = 4.0_dp
  call init_panel_geometry(shift + scale*[0.0_dp, 0.0_dp, 0.0_dp], &
                           shift + scale*[0.0_dp, 2.0_dp, 0.0_dp], &
                           shift + scale*[2.0_dp, 0.0_dp, 0.0_dp], transformed_geometry, status)
  call assert_equal_i32(status, panel_geometry_ok, 'transformed reversed triangle fixture should be valid')
  call assert_close_dp( &
    point_triangle_distance(transformed_geometry, shift + scale*[1.5_dp, 1.5_dp, 0.25_dp]), &
    scale*0.75_dp, 1.0e-13_dp, 'distance must preserve orientation, translation, and uniform scaling' &
    )
  call test_end()

  call test_begin('invalid_triangle_statuses')
  call init_panel_geometry( &
    [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 0.0_dp, 0.0_dp], [2.0_dp, 0.0_dp, 0.0_dp], geometry, status &
    )
  call assert_equal_i32(status, panel_geometry_degenerate, 'collinear triangle must be rejected')
  nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
  call init_panel_geometry( &
    [nan_value, 0.0_dp, 0.0_dp], [1.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 1.0_dp, 0.0_dp], geometry, status &
    )
  call assert_equal_i32(status, panel_geometry_nonfinite, 'non-finite triangle must be rejected')
  call test_end()

  call test_begin('panel_tree_bbox_contains_all_vertices')
  v0(:, 1) = [-10.0_dp, 0.0_dp, 0.0_dp]
  v1(:, 1) = [10.0_dp, 0.0_dp, 0.0_dp]
  v2(:, 1) = [0.0_dp, 0.01_dp, 0.0_dp]
  v0(:, 2) = [0.0_dp, 2.0_dp, 0.0_dp]
  v1(:, 2) = [0.1_dp, 2.0_dp, 0.0_dp]
  v2(:, 2) = [0.0_dp, 2.1_dp, 0.0_dp]
  vertices(:, 1:2) = v0
  vertices(:, 3:4) = v1
  vertices(:, 5:6) = v2
  options = fmm_options_type()
  options%leaf_max = 1_i32
  options%order = 2_i32
  call build_panel_plan(plan, v0, v1, v2, options)
  call assert_true(all(minval(vertices, dim=2) >= plan%node_center(:, 1) - plan%node_half_size(:, 1)), &
                   'root lower bbox must contain every panel vertex')
  call assert_true(all(maxval(vertices, dim=2) <= plan%node_center(:, 1) + plan%node_half_size(:, 1)), &
                   'root upper bbox must contain every panel vertex')
  max_vertex_radius = maxval(sqrt(sum( &
                                  (vertices - spread(plan%node_center(:, 1), 2, size(vertices, 2)))**2, dim=1 &
                                  )))
  call assert_true( &
    plan%node_radius(1) + 1.0e-14_dp*max(1.0_dp, max_vertex_radius) >= max_vertex_radius, &
    'root radius must contain every panel vertex' &
    )
  call destroy_plan(plan)
  call test_end()

  call test_begin('exact_panel_p2m_moments')
  v0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
  v1(:, 1) = [2.0_dp, 0.0_dp, 0.0_dp]
  v2(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
  options = fmm_options_type()
  options%leaf_max = 8_i32
  options%order = 2_i32
  call build_panel_plan(plan, v0(:, :1), v1(:, :1), v2(:, :1), options)
  center = plan%node_center(:, 1)
  alpha_zero = plan%alpha_map(0, 0, 0)
  alpha_x = plan%alpha_map(1, 0, 0)
  alpha_xx = plan%alpha_map(2, 0, 0)
  call assert_close_dp(plan%source_p2m_basis(alpha_zero, 1), 1.0_dp, 1.0e-14_dp, 'panel monopole basis mismatch')
  call assert_close_dp(plan%source_p2m_basis(alpha_x, 1), triangle_mean_x - center(1), 1.0e-14_dp, &
                       'panel dipole basis mismatch')
  call assert_close_dp(plan%source_p2m_basis(alpha_xx, 1), &
                       0.5_dp*(triangle_mean_x2 - 2.0_dp*center(1)*triangle_mean_x + center(1)**2), 1.0e-14_dp, &
                       'panel quadrupole basis mismatch')
  call destroy_plan(plan)
  call test_end()

  call test_begin('panel_outside_box_uses_analytic_fallback')
  options = fmm_options_type()
  options%target_box_min = [-1.0_dp, -1.0_dp, -1.0_dp]
  options%target_box_max = [1.0_dp, 1.0_dp, 1.0_dp]
  call build_panel_plan(plan, v0(:, :1), v1(:, :1), v2(:, :1), options)
  call init_panel_geometry(v0(:, 1), v1(:, 1), v2(:, 1), geometry, status)
  q = [2.5_dp]
  call update_state(plan, state, q)
  target = [2.0_dp, 1.5_dp, 0.75_dp]
  call eval_point(plan, state, target, field)
  call eval_potential_point(plan, state, target, phi)
  call panel_potential_field(geometry, q(1), target, panel_side_principal_value, phi_ref, field_ref)
  call assert_true(all(abs(field - field_ref/k_coulomb) < 1.0e-13_dp), 'fallback field must use analytic panel kernel')
  call assert_close_dp(phi, phi_ref/k_coulomb, 1.0e-13_dp, 'fallback potential must use analytic panel kernel')
  call destroy_state(state)
  call destroy_plan(plan)
  call test_end()

  call test_summary()
end program test_panel_geometry_near

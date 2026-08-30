!> Triangle geometry, cached moments, cubature, and surface-side invariants.
program test_panel_moments
  use bem_kinds, only: dp, i32
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_quadrature, only: panel_quadrature_plan_type, build_panel_quadrature
  use bem_types, only: mesh_type
  use bem_mesh, only: init_mesh
  use bem_panel_surface_sides, only: resolve_panel_surface_sides, panel_surface_side_ok, panel_surface_side_open
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32, &
                          assert_close_dp, assert_allclose_1d
  implicit none

  type(panel_geometry_type) :: geom, cyclic, reversed, scaled
  type(panel_quadrature_plan_type) :: quad
  type(mesh_type) :: mesh
  real(dp) :: v0(3), v1(3), v2(3), cubature_value
  real(dp) :: tetra_v0(3, 4), tetra_v1(3, 4), tetra_v2(3, 4)
  real(dp), parameter :: unit_moment1(3) = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 0.0_dp]
  real(dp), parameter :: unit_moment2(9) = [ &
                         1.0_dp/12.0_dp, 1.0_dp/24.0_dp, 0.0_dp, &
                         1.0_dp/24.0_dp, 1.0_dp/12.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp &
                         ]
  character(len=128) :: message
  integer(i32) :: status
  integer :: x_power, y_power

  v0 = [0.0_dp, 0.0_dp, 0.0_dp]
  v1 = [1.0_dp, 0.0_dp, 0.0_dp]
  v2 = [0.0_dp, 1.0_dp, 0.0_dp]
  call test_init(6)

  call test_begin('unit_triangle_geometry_and_raw_moments')
  call init_panel_geometry(v0, v1, v2, geom, status)
  call assert_equal_i32(status, panel_geometry_ok, 'valid triangle status mismatch')
  call assert_close_dp(geom%area, 0.5_dp, 1.0e-15_dp, 'triangle area mismatch')
  call assert_allclose_1d(geom%centroid, [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.0_dp], 1.0e-15_dp, 'centroid mismatch')
  call assert_allclose_1d(geom%normal, [0.0_dp, 0.0_dp, 1.0_dp], 1.0e-15_dp, 'normal mismatch')
  call assert_close_dp(geom%moment0, 0.5_dp, 1.0e-15_dp, 'zeroth moment mismatch')
  call assert_allclose_1d(geom%moment1, unit_moment1, 1.0e-15_dp, 'first moment mismatch')
  call assert_allclose_1d(reshape(geom%moment2, [9]), unit_moment2, 1.0e-15_dp, 'second moment mismatch')
  call test_end()

  call test_begin('mesh_owns_panel_geometry')
  call init_mesh(mesh, reshape(v0, [3, 1]), reshape(v1, [3, 1]), reshape(v2, [3, 1]))
  call assert_close_dp(mesh%panel_area(1), 0.5_dp, 1.0e-15_dp, 'mesh panel area mismatch')
  call assert_allclose_1d(mesh%panel_moment1(:, 1), unit_moment1, 1.0e-15_dp, &
                          'mesh panel moment mismatch')
  call assert_allclose_1d(reshape(mesh%panel_moment2(:, :, 1), [9]), unit_moment2, 1.0e-15_dp, &
                          'mesh second moment mismatch')
  call assert_allclose_1d(mesh%panel_edge_length(:, 1), [1.0_dp, sqrt(2.0_dp), 1.0_dp], 1.0e-15_dp, &
                          'mesh edge cache mismatch')
  call assert_close_dp(sum(mesh%panel_quad_weight(:, 1)), 0.5_dp, 1.0e-15_dp, 'mesh quadrature weight mismatch')
  call test_end()

  call test_begin('open_surface_explicit_side')
  call init_mesh(mesh, reshape(v0, [3, 1]), reshape(v1, [3, 1]), reshape(v2, [3, 1]))
  call resolve_panel_surface_sides(mesh, 'normal_plus', status, message)
  call assert_equal_i32(status, panel_surface_side_ok, 'normal_plus side resolution failed')
  call assert_equal_i32(mesh%elem_vacuum_sign(1), 1_i32, 'vacuum sign mismatch')
  call assert_allclose_1d(mesh%vacuum_normals(:, 1), [0.0_dp, 0.0_dp, 1.0_dp], 1.0e-15_dp, &
                          'vacuum normal mismatch')
  call resolve_panel_surface_sides(mesh, 'outward_closed', status, message)
  call assert_equal_i32(status, panel_surface_side_open, 'open surface must reject outward_closed')
  call test_end()

  call test_begin('closed_surface_outward_side_is_winding_independent')
  tetra_v0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
  tetra_v1(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
  tetra_v2(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
  tetra_v0(:, 2) = [0.0_dp, 0.0_dp, 0.0_dp]
  tetra_v1(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
  tetra_v2(:, 2) = [0.0_dp, 0.0_dp, 1.0_dp]
  tetra_v0(:, 3) = [0.0_dp, 0.0_dp, 0.0_dp]
  tetra_v1(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp]
  tetra_v2(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
  tetra_v0(:, 4) = [1.0_dp, 0.0_dp, 0.0_dp]
  tetra_v1(:, 4) = [0.0_dp, 1.0_dp, 0.0_dp]
  tetra_v2(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
  call init_mesh(mesh, tetra_v0, tetra_v1, tetra_v2)
  call resolve_panel_surface_sides(mesh, 'outward_closed', status, message)
  call assert_equal_i32(status, panel_surface_side_ok, 'closed tetrahedron side resolution failed')
  call assert_true(all(mesh%elem_vacuum_sign == 1_i32), 'outward tetrahedron signs must be positive')
  call assert_true(all(sum(mesh%vacuum_normals*(mesh%centers - 0.25_dp), dim=1) > 0.0_dp), &
                   'vacuum normals must face away from tetrahedron interior')
  call init_mesh(mesh, tetra_v0, tetra_v2, tetra_v1)
  call resolve_panel_surface_sides(mesh, 'outward_closed', status, message)
  call assert_equal_i32(status, panel_surface_side_ok, 'inward-wound tetrahedron side resolution failed')
  call assert_true(all(mesh%elem_vacuum_sign == -1_i32), 'inward tetrahedron signs must be negative')
  call assert_true(all(sum(mesh%vacuum_normals*(mesh%centers - 0.25_dp), dim=1) > 0.0_dp), &
                   'inward-wound tetrahedron vacuum normals must face outward')
  call test_end()

  call test_begin('degree5_cubature_integrates_all_monomials')
  call build_panel_quadrature(geom, quad)
  do x_power = 0, 5
    do y_power = 0, 5 - x_power
      cubature_value = sum(quad%weight*quad%position(1, :)**x_power*quad%position(2, :)**y_power)
      call assert_close_dp(cubature_value, unit_triangle_monomial_integral(x_power, y_power), 1.0e-14_dp, &
                           'degree-5 cubature monomial mismatch')
    end do
  end do
  call test_end()

  call test_begin('moment_permutation_and_scaling_invariants')
  call init_panel_geometry(v1, v2, v0, cyclic, status)
  call assert_equal_i32(status, panel_geometry_ok, 'cyclic triangle status mismatch')
  call assert_close_dp(cyclic%area, geom%area, 1.0e-15_dp, 'cyclic area mismatch')
  call assert_allclose_1d(cyclic%moment1, geom%moment1, 1.0e-15_dp, 'cyclic first moment mismatch')
  call assert_allclose_1d(reshape(cyclic%moment2, [9]), reshape(geom%moment2, [9]), 1.0e-15_dp, 'cyclic moments mismatch')
  call assert_allclose_1d(cyclic%normal, geom%normal, 1.0e-15_dp, 'cyclic normal mismatch')
  call init_panel_geometry(v0, v2, v1, reversed, status)
  call assert_equal_i32(status, panel_geometry_ok, 'reversed triangle status mismatch')
  call assert_allclose_1d(reversed%normal, -geom%normal, 1.0e-15_dp, 'reversed normal mismatch')
  call assert_allclose_1d(reversed%moment1, geom%moment1, 1.0e-15_dp, 'reversed first moment mismatch')
  call assert_allclose_1d(reshape(reversed%moment2, [9]), reshape(geom%moment2, [9]), 1.0e-15_dp, 'reversed moments mismatch')
  call init_panel_geometry(2.0_dp*v0, 2.0_dp*v1, 2.0_dp*v2, scaled, status)
  call assert_equal_i32(status, panel_geometry_ok, 'scaled triangle status mismatch')
  call assert_close_dp(scaled%moment0, 4.0_dp*geom%moment0, 1.0e-14_dp, 'scaled zeroth moment mismatch')
  call assert_allclose_1d(scaled%moment1, 8.0_dp*geom%moment1, 1.0e-14_dp, 'scaled first moment mismatch')
  call assert_allclose_1d(reshape(scaled%moment2, [9]), 16.0_dp*reshape(geom%moment2, [9]), 1.0e-14_dp, &
                          'scaled second moment mismatch')
  call test_end()

  call test_summary()

contains

  pure real(dp) function unit_triangle_monomial_integral(px, py) result(integral)
    integer, intent(in) :: px, py

    integral = real(factorial(px)*factorial(py), dp)/real(factorial(px + py + 2), dp)
  end function unit_triangle_monomial_integral

  pure integer function factorial(n) result(value)
    integer, intent(in) :: n
    integer :: factor

    value = 1
    do factor = 2, n
      value = value*factor
    end do
  end function factorial

end program test_panel_moments

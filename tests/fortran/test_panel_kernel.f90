!> Analytic P0 triangle kernel against independent quadrature and side limits.
program test_panel_kernel
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_quadrature, only: panel_oracle_potential_field, panel_singular_potential_oracle
  use bem_panel_kernel, only: panel_potential_field, panel_side_principal_value, panel_side_normal_plus, &
                              panel_side_normal_minus
  use bem_panel_self_terms, only: panel_on_surface_integrals
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_equal_i32, assert_close_dp, &
                          assert_allclose_1d
  implicit none

  type(panel_geometry_type) :: geom, moved, small_geom, neighbor_geom
  real(dp) :: v0(3), v1(3), v2(3), target(3), phi, field(3), phi_ref, field_ref(3)
  real(dp) :: phi_plus, phi_minus, field_plus(3), field_minus(3), phi_self, self_ref
  real(dp) :: self_integral, self_field_integral(3)
  real(dp) :: rotation(3, 3), shift(3), moved_target(3), moved_phi, moved_field(3)
  real(dp) :: scale_value, offset_value, neighbor_field_plus(3), neighbor_field_minus(3)
  real(dp), parameter :: charge = 2.5e-12_dp
  integer(i32) :: status, scale_idx, offset_idx
  real(dp), parameter :: scales(5) = [1.0e-9_dp, 1.0e-6_dp, 1.0_dp, 1.0e3_dp, 1.0e6_dp]
  real(dp), parameter :: offsets(3) = [0.0_dp, 1.0e3_dp, 1.0e6_dp]

  v0 = [0.0_dp, 0.0_dp, 0.0_dp]
  v1 = [1.0_dp, 0.0_dp, 0.0_dp]
  v2 = [0.0_dp, 1.0_dp, 0.0_dp]
  call init_panel_geometry(v0, v1, v2, geom, status)

  call test_init(7)

  call test_begin('off_surface_analytic_matches_oracle')
  call assert_equal_i32(status, panel_geometry_ok, 'valid geometry status mismatch')
  target = [0.23_dp, 0.31_dp, 0.47_dp]
  call panel_potential_field(geom, charge, target, panel_side_principal_value, phi, field)
  call panel_oracle_potential_field(geom, charge, target, 32_i32, phi_ref, field_ref)
  call assert_close_dp(phi, phi_ref, 2.0e-11_dp*max(1.0_dp, abs(phi_ref)), 'off-surface potential mismatch')
  call assert_allclose_1d(field, field_ref, 2.0e-10_dp*max(1.0_dp, maxval(abs(field_ref))), &
                          'off-surface field mismatch')
  call test_end()

  call test_begin('geometry_status_is_scale_and_translation_invariant')
  do scale_idx = 1, size(scales)
    scale_value = scales(scale_idx)
    do offset_idx = 1, size(offsets)
      offset_value = offsets(offset_idx)
      shift = [offset_value, -0.5_dp*offset_value, 0.25_dp*offset_value]
      call init_panel_geometry( &
        shift, shift + [scale_value, 0.0_dp, 0.0_dp], shift + [0.0_dp, scale_value, 0.0_dp], moved, status &
        )
      call assert_equal_i32(status, panel_geometry_ok, 'translated/scaled triangle was classified as degenerate')
    end do
  end do
  call test_end()

  call test_begin('coplanar_neighbor_does_not_receive_normal_jump')
  scale_value = 1.0e-7_dp
  call init_panel_geometry( &
    [0.0_dp, 0.0_dp, 0.0_dp], [scale_value, 0.0_dp, 0.0_dp], [0.0_dp, scale_value, 0.0_dp], &
    small_geom, status &
    )
  call assert_equal_i32(status, panel_geometry_ok, 'small panel geometry status mismatch')
  call init_panel_geometry( &
    [scale_value, 0.0_dp, 0.0_dp], [scale_value, scale_value, 0.0_dp], &
    [0.0_dp, scale_value, 0.0_dp], neighbor_geom, status &
    )
  call assert_equal_i32(status, panel_geometry_ok, 'small neighbor geometry status mismatch')
  call panel_potential_field( &
    small_geom, charge, neighbor_geom%centroid, panel_side_normal_plus, phi_plus, neighbor_field_plus &
    )
  call panel_potential_field( &
    small_geom, charge, neighbor_geom%centroid, panel_side_normal_minus, phi_minus, neighbor_field_minus &
    )
  call assert_allclose_1d( &
    neighbor_field_plus, neighbor_field_minus, &
    1.0e-12_dp*max(1.0_dp, maxval(abs(neighbor_field_plus))), &
    'coplanar point outside the panel received a normal jump' &
    )
  call test_end()

  call test_begin('potential_is_continuous')
  target = geom%centroid + 1.0e-8_dp*geom%normal
  call panel_potential_field(geom, charge, target, panel_side_principal_value, phi_plus, field_plus)
  target = geom%centroid - 1.0e-8_dp*geom%normal
  call panel_potential_field(geom, charge, target, panel_side_principal_value, phi_minus, field_minus)
  call assert_close_dp(phi_plus, phi_minus, 1.0e-12_dp*max(1.0_dp, abs(phi_plus)), 'two-sided potential mismatch')
  call test_end()

  call test_begin('normal_field_jump')
  call panel_potential_field(geom, charge, geom%centroid, panel_side_normal_plus, phi_plus, field_plus)
  call panel_potential_field(geom, charge, geom%centroid, panel_side_normal_minus, phi_minus, field_minus)
  call assert_allclose_1d( &
    field_plus - field_minus, charge/(geom%area*eps0)*geom%normal, &
    1.0e-12_dp*max(1.0_dp, abs(charge/(geom%area*eps0))), 'normal field jump mismatch' &
    )
  call assert_close_dp(phi_plus, phi_minus, 1.0e-14_dp, 'one-sided potential mismatch')
  call test_end()

  call test_begin('centroid_self_potential')
  call panel_potential_field(geom, charge, geom%centroid, panel_side_principal_value, phi_self, field)
  call panel_on_surface_integrals(geom, geom%centroid, self_integral, self_field_integral)
  call panel_singular_potential_oracle(geom, charge, geom%centroid, 32_i32, self_ref)
  call assert_close_dp(phi_self, self_ref, 2.0e-11_dp*max(1.0_dp, abs(self_ref)), 'self potential mismatch')
  call assert_close_dp(phi_self, 8.9875517923e9_dp*charge/geom%area*self_integral, 1.0e-13_dp, &
                       'self-term owner potential mismatch')
  call assert_allclose_1d(field, 8.9875517923e9_dp*charge/geom%area*self_field_integral, 1.0e-13_dp, &
                          'self-term owner field mismatch')
  call test_end()

  call test_begin('rigid_transform_and_charge_scaling')
  rotation(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
  rotation(:, 2) = [-1.0_dp, 0.0_dp, 0.0_dp]
  rotation(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp]
  shift = [2.0_dp, -3.0_dp, 0.7_dp]
  call init_panel_geometry(matmul(rotation, v0) + shift, matmul(rotation, v1) + shift, &
                           matmul(rotation, v2) + shift, moved, status)
  target = [0.23_dp, 0.31_dp, 0.47_dp]
  moved_target = matmul(rotation, target) + shift
  call panel_potential_field(geom, charge, target, panel_side_principal_value, phi, field)
  call panel_potential_field(moved, 2.0_dp*charge, moved_target, panel_side_principal_value, moved_phi, moved_field)
  call assert_close_dp(moved_phi, 2.0_dp*phi, 1.0e-12_dp*max(1.0_dp, abs(phi)), 'transformed potential mismatch')
  call assert_allclose_1d(moved_field, 2.0_dp*matmul(rotation, field), &
                          1.0e-11_dp*max(1.0_dp, maxval(abs(field))), 'transformed field mismatch')
  call test_end()

  call test_summary()
end program test_panel_kernel

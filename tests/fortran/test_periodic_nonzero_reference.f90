program test_periodic_nonzero_reference
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type
  use bem_mesh, only: init_mesh
  use bem_coulomb_fmm_periodic_nonzero_reference, only: eval_periodic_nonzero_panel_reference
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_close_dp, assert_allclose_1d, &
                          assert_true
  implicit none

  type(mesh_type) :: mesh
  real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2), target(3), phi, field(3), periodic_phi, periodic_field(3)
  real(dp), parameter :: total_charge = 2.0e-12_dp

  call test_init(2)

  v0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
  v1(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
  v2(:, 1) = [1.0_dp, 1.0_dp, 0.0_dp]
  v0(:, 2) = [0.0_dp, 0.0_dp, 0.0_dp]
  v1(:, 2) = [1.0_dp, 1.0_dp, 0.0_dp]
  v2(:, 2) = [0.0_dp, 1.0_dp, 0.0_dp]
  call init_mesh(mesh, v0, v1, v2, q0=[0.5_dp*total_charge, 0.5_dp*total_charge])
  target = [0.31_dp, 0.27_dp, 0.5_dp]

  call test_begin('uniform_sheet_has_no_nonzero_mode')
  call eval_periodic_nonzero_panel_reference(mesh, target, 1.0_dp, 1.0_dp, 5_i32, 20_i32, phi, field)
  call assert_close_dp(phi, 0.0_dp, 2.0e-10_dp, 'uniform sheet nonzero potential mismatch')
  call assert_allclose_1d(field, [0.0_dp, 0.0_dp, 0.0_dp], 2.0e-9_dp, 'uniform sheet nonzero field mismatch')
  call test_end()

  call test_begin('nonneutral_reference_is_nonzero_and_periodic')
  mesh%q_elem = [total_charge, 0.0_dp]
  call eval_periodic_nonzero_panel_reference(mesh, target, 1.0_dp, 1.0_dp, 5_i32, 20_i32, phi, field)
  call assert_true(abs(phi) > 1.0e-8_dp, 'fixture must exercise a nonzero-mode potential')
  call assert_true(maxval(abs(field)) > 1.0e-8_dp, 'fixture must exercise a nonzero-mode field')
  call eval_periodic_nonzero_panel_reference( &
    mesh, target + [1.0_dp, -1.0_dp, 0.0_dp], 1.0_dp, 1.0_dp, 5_i32, 20_i32, &
    periodic_phi, periodic_field &
    )
  call assert_close_dp(periodic_phi, phi, 2.0e-12_dp, 'periodic potential mismatch')
  call assert_allclose_1d(periodic_field, field, 2.0e-11_dp, 'periodic field mismatch')
  call test_end()

  call test_summary()
end program test_periodic_nonzero_reference

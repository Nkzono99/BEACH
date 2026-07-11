program test_periodic2_nonzero_tail
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_periodic_nonzero_tail, only: &
    periodic_nonzero_tail_mode_type, init_periodic_nonzero_tail_mode, &
    eval_periodic_nonzero_tail_correction, periodic_nonzero_tail_linearity, &
    periodic_nonzero_tail_plan_type, build_periodic_nonzero_tail_plan, eval_periodic_nonzero_tail_plan
  use bem_types, only: mesh_type
  use bem_mesh, only: init_mesh
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_close_dp, assert_equal_i32
  implicit none

  type(periodic_nonzero_tail_mode_type) :: mode
  real(dp) :: potential_below, potential_above, field_below(3), field_above(3), potential, field(3), eta
  integer(i32) :: status
  type(periodic_nonzero_tail_plan_type) :: plan
  type(mesh_type) :: mesh
  real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1), translated_potential

  call test_init(4)

  call test_begin('screened correction is continuous at the handoff plane')
  call init_periodic_nonzero_tail_mode( &
    kx=2.0_dp, ky=0.0_dp, kappa=3.0_dp, handoff_z=1.0_dp, &
    incident_cos=1.0_dp, incident_sin=0.0_dp, mode=mode, status=status &
    )
  call assert_equal_i32(status, 0_i32, 'tail mode status mismatch')
  call eval_periodic_nonzero_tail_correction(mode, [0.0_dp, 0.0_dp, 1.0_dp - 1.0e-10_dp], &
                                             potential_below, field_below)
  call eval_periodic_nonzero_tail_correction(mode, [0.0_dp, 0.0_dp, 1.0_dp + 1.0e-10_dp], &
                                             potential_above, field_above)
  call assert_close_dp(potential_below, potential_above, 1.0e-9_dp, 'tail potential is discontinuous')
  call assert_close_dp(field_below(1), field_above(1), 1.0e-9_dp, 'tail tangential field is discontinuous')
  call assert_close_dp(field_below(3), field_above(3), 1.0e-9_dp, 'tail normal field is discontinuous')
  call test_end()

  call test_begin('panel mode plan is periodic and continuous')
  v0(:, 1) = [0.1_dp, 0.1_dp, 0.0_dp]
  v1(:, 1) = [0.4_dp, 0.1_dp, 0.0_dp]
  v2(:, 1) = [0.1_dp, 0.4_dp, 0.0_dp]
  call init_mesh(mesh, v0, v1, v2, q0=[2.0e-12_dp])
  call build_periodic_nonzero_tail_plan( &
    mesh, length_x=1.0_dp, length_y=1.0_dp, handoff_z=1.0_dp, kappa=2.0_dp, &
    mode_layers=1_i32, quadrature_order=4_i32, charge_abs=1.0_dp, temperature_j=1.0_dp, &
    plan=plan, status=status &
    )
  call assert_equal_i32(status, 0_i32, 'panel tail plan status mismatch')
  call assert_equal_i32(plan%nmodes, 4_i32, 'canonical tail mode count mismatch')
  call eval_periodic_nonzero_tail_plan(plan, [0.23_dp, 0.37_dp, 1.0_dp], potential, field)
  call eval_periodic_nonzero_tail_plan(plan, [1.23_dp, 0.37_dp, 1.0_dp], translated_potential, field_above)
  call assert_close_dp(translated_potential, potential, 1.0e-20_dp, 'tail plan periodic translation mismatch')
  call eval_periodic_nonzero_tail_plan(plan, [0.23_dp, 0.37_dp, 1.0_dp - 1.0e-10_dp], &
                                       potential_below, field_below)
  call eval_periodic_nonzero_tail_plan(plan, [0.23_dp, 0.37_dp, 1.0_dp + 1.0e-10_dp], &
                                       potential_above, field_above)
  call assert_close_dp(potential_below, potential_above, 1.0e-14_dp, 'panel tail plan potential discontinuity')
  call assert_close_dp(field_below(3), field_above(3), 1.0e-13_dp, 'panel tail plan field discontinuity')
  call test_end()

  call test_begin('zero susceptibility is exactly the free-space limit')
  call init_periodic_nonzero_tail_mode( &
    kx=2.0_dp, ky=1.0_dp, kappa=0.0_dp, handoff_z=1.0_dp, &
    incident_cos=3.0_dp, incident_sin=-2.0_dp, mode=mode, status=status &
    )
  call eval_periodic_nonzero_tail_correction(mode, [0.3_dp, 0.4_dp, 2.0_dp], potential, field)
  call assert_close_dp(potential, 0.0_dp, 1.0e-15_dp, 'vacuum-limit potential correction mismatch')
  call assert_close_dp(maxval(abs(field)), 0.0_dp, 1.0e-15_dp, 'vacuum-limit field correction mismatch')
  call test_end()

  call test_begin('linearity metric uses the transmitted mode amplitude')
  call init_periodic_nonzero_tail_mode( &
    kx=2.0_dp, ky=0.0_dp, kappa=3.0_dp, handoff_z=1.0_dp, &
    incident_cos=1.0_dp, incident_sin=0.0_dp, mode=mode, status=status &
    )
  eta = periodic_nonzero_tail_linearity(mode, charge_abs=2.0_dp, temperature_j=4.0_dp)
  call assert_close_dp(eta, 0.5_dp*mode%transmission, 1.0e-15_dp, 'tail linearity metric mismatch')
  call test_end()

  call test_summary()
end program test_periodic2_nonzero_tail

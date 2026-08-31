!> Matching-plane charge-driven Zhao 応答の物理分岐と入力契約を検証する。
program test_matching_plane_zhao
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, qe
  use bem_matching_plane_zhao, only: &
    matching_plane_zhao_model_type, matching_plane_zhao_diagnostics_type, &
    matching_plane_zhao_ok, matching_plane_zhao_invalid_argument, &
    matching_plane_zhao_no_physical_solution, matching_plane_zhao_numerical_failure, &
    matching_plane_zhao_ambiguous_solution
  use test_support, only: &
    test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32, &
    assert_close_dp, assert_allclose_1d
  implicit none

  real(dp), parameter :: electron_mass_kg = 9.1093837015e-31_dp
  real(dp), parameter :: proton_mass_kg = 1.67262192369e-27_dp
  real(dp), parameter :: ion_density_m3 = 8.7e6_dp
  real(dp), parameter :: electron_temperature_ev = 12.0_dp
  real(dp), parameter :: configured_photoelectron_temperature_ev = 12.0_dp
  real(dp), parameter :: drift_mps = 4.0529988897111727e5_dp
  real(dp), parameter :: type_a_photoelectron_temperature_ev = 2.2_dp
  real(dp), parameter :: type_a_phi0_v = 2.9712182827319435_dp
  real(dp), parameter :: type_a_phi_m_v = -0.8169121871620854_dp
  real(dp), parameter :: type_a_source_density_m3 = 5.5425625842204072e7_dp
  ! Type-A既知解を固定値化し、productionのrho積分でtest入力を再生成しない。
  real(dp), parameter :: type_a_input(5) = [ &
                         1.4187346568707933e-11_dp, 1.3754433596232731e13_dp, &
                         type_a_photoelectron_temperature_ev, 0.0_dp, 0.0_dp &
                         ]
  type(matching_plane_zhao_model_type) :: model
  type(matching_plane_zhao_diagnostics_type) :: diagnostics, inactive_diagnostics
  type(matching_plane_zhao_diagnostics_type) :: energy_a_diagnostics, energy_b_diagnostics
  real(dp) :: input(5), output(6), inactive_output(6), feedback_scales(4), energy_a_output(6)
  real(dp) :: electron_thermal_speed_mps, expected_electron_density_m3
  integer(i32) :: status
  character(len=512) :: message

  call test_init(8)

  call test_begin('zero_field_without_photoelectrons_is_degenerate_zhao_b')
  call initialize_model('auto', configured_photoelectron_temperature_ev)
  call assert_true(model%is_initialized(), 'Zhao matching model did not retain initialization')
  input = 0.0_dp
  call model%evaluate(input, output, status, message, diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'zero-field Zhao matching solve failed: '//trim(message))
  electron_thermal_speed_mps = sqrt(2.0_dp*qe*electron_temperature_ev/electron_mass_kg)
  expected_electron_density_m3 = &
    2.0_dp*ion_density_m3/(1.0_dp + erf(drift_mps/electron_thermal_speed_mps))
  call assert_true(diagnostics%branch == 'B', 'zero-field no-PE state must use degenerate Zhao-B')
  call assert_close_dp( &
    diagnostics%effective_photoelectron_temperature_ev, &
    configured_photoelectron_temperature_ev, 0.0_dp, &
    'zero-PE query did not use the configured photoelectron-temperature fallback' &
    )
  call assert_close_dp(output(1), 0.0_dp, 0.0_dp, 'zero-field matching potential changed')
  call assert_close_dp( &
    diagnostics%ambient_electron_density_m3, expected_electron_density_m3, &
    1.0e-8_dp*expected_electron_density_m3, 'zero-field ambient electron density mismatch' &
    )
  call assert_true(all(output(2:3) > 0.0_dp), 'zero-field ambient inward flux must be positive')
  call assert_true(all(output(4:6) == 0.0_dp), 'zero-field access/barrier potentials must vanish')
  call model%get_feedback_scales(feedback_scales, status, message)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'Zhao feedback scales failed: '//trim(message))
  call assert_true( &
    all(feedback_scales(1:2) > 0.0_dp) .and. all(feedback_scales(3:4) == 0.0_dp), &
    'Zhao feedback scales must disable ambient outward dependencies' &
    )
  call test_end()

  call test_begin('minimum_energy_selects_the_lower_energy_positive_branch')
  input = type_a_input
  call initialize_model('a', type_a_photoelectron_temperature_ev, 'minimum_energy')
  call model%evaluate(input, energy_a_output, status, message, energy_a_diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'minimum-energy Zhao-A solve failed: '//trim(message))
  call initialize_model('b', type_a_photoelectron_temperature_ev, 'minimum_energy')
  call model%evaluate(input, output, status, message, energy_b_diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'minimum-energy Zhao-B solve failed: '//trim(message))
  call assert_true( &
    energy_a_diagnostics%potential_energy_j_m2 < 0.0_dp .and. &
    energy_b_diagnostics%potential_energy_j_m2 < 0.0_dp, &
    'explicit Zhao roots did not report finite negative potential energies' &
    )
  call initialize_model('auto', type_a_photoelectron_temperature_ev, 'minimum_energy')
  call model%evaluate(input, output, status, message, diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'minimum-energy Zhao auto selection failed: '//trim(message))
  if (energy_a_diagnostics%potential_energy_j_m2 < energy_b_diagnostics%potential_energy_j_m2) then
    call assert_true(diagnostics%branch == 'A', 'minimum-energy auto selection did not choose Zhao-A')
  else
    call assert_true(diagnostics%branch == 'B', 'minimum-energy auto selection did not choose Zhao-B')
  end if
  call assert_close_dp( &
    diagnostics%potential_energy_j_m2, &
    min(energy_a_diagnostics%potential_energy_j_m2, energy_b_diagnostics%potential_energy_j_m2), &
    1.0e-6_dp*max( &
    abs(energy_a_diagnostics%potential_energy_j_m2), &
    abs(energy_b_diagnostics%potential_energy_j_m2) &
    ), &
    'minimum-energy auto selection returned the wrong energy' &
    )
  call test_end()

  call test_begin('explicit_zhao_b_solves_a_positive_field_profile')
  call initialize_model('b', configured_photoelectron_temperature_ev)
  input = 0.0_dp
  input(1) = 0.02_dp*eps0
  call model%evaluate(input, output, status, message, diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'positive-field Zhao-B solve failed: '//trim(message))
  call assert_true(diagnostics%branch == 'B' .and. output(1) > 0.0_dp, &
                   'positive-field Zhao-B response has the wrong branch or potential sign')
  call assert_true( &
    ieee_is_finite(diagnostics%minimum_field_squared_hat) .and. &
    diagnostics%minimum_field_squared_hat >= -1.0e-7_dp, &
    'positive-field Zhao-B path contains an imaginary-field interval' &
    )
  call assert_true(all(output(4:6) == 0.0_dp), 'Zhao-B access/barrier potentials must use the upstream gauge')
  call test_end()

  call test_begin('explicit_zhao_a_solves_a_positive_nonmonotonic_profile')
  call initialize_model('a', type_a_photoelectron_temperature_ev)
  call model%evaluate(type_a_input, output, status, message, diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'positive-field Zhao-A solve failed: '//trim(message))
  call assert_true(diagnostics%branch == 'A', 'positive-field Zhao-A response changed branch')
  call assert_close_dp(output(1), type_a_phi0_v, 5.0e-5_dp, 'positive-field Zhao-A potential mismatch')
  call assert_close_dp(output(4), type_a_phi_m_v, 5.0e-5_dp, 'Zhao-A electron access mismatch')
  call assert_close_dp(output(5), 0.0_dp, 0.0_dp, 'Zhao-A ion access must use the upstream gauge')
  call assert_close_dp(output(6), type_a_phi_m_v, 5.0e-5_dp, 'Zhao-A PE return barrier mismatch')
  call assert_close_dp( &
    diagnostics%effective_photoelectron_temperature_ev, type_a_photoelectron_temperature_ev, 0.0_dp, &
    'PE mean normal energy was not used as the half-Maxwellian temperature' &
    )
  call assert_close_dp( &
    diagnostics%photoelectron_source_density_m3, type_a_source_density_m3, &
    1.0e-12_dp*type_a_source_density_m3, &
    'PE outward flux was not mapped to the half-Maxwellian source density' &
    )
  call assert_true( &
    ieee_is_finite(diagnostics%minimum_field_squared_hat) .and. &
    diagnostics%minimum_field_squared_hat >= -1.0e-7_dp, &
    'positive-field Zhao-A path contains an imaginary-field interval' &
    )
  call test_end()

  call test_begin('auto_positive_field_fails_closed_when_uniqueness_is_uncertain')
  call initialize_model('auto', type_a_photoelectron_temperature_ev)
  call model%evaluate(type_a_input, output, status, message, diagnostics)
  call assert_true( &
    status == matching_plane_zhao_ambiguous_solution .or. &
    status == matching_plane_zhao_numerical_failure, &
    'auto positive-field Zhao query did not report a fail-closed uniqueness status' &
    )
  call assert_true(all(output == 0.0_dp), 'failed auto branch selection returned a partial response')
  call test_end()

  call test_begin('auto_negative_field_selects_zhao_c_and_ignores_ambient_outflow')
  call initialize_model('auto', configured_photoelectron_temperature_ev)
  input = 0.0_dp
  input(1) = -0.02_dp*eps0
  call model%evaluate(input, output, status, message, diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'nonzero-field Zhao matching solve failed: '//trim(message))
  call assert_true(diagnostics%branch == 'C', 'negative interface field must select Zhao-C')
  call assert_true( &
    ieee_is_finite(diagnostics%residual_norm) .and. diagnostics%residual_norm <= 1.0e-9_dp, &
    'nonzero-field Zhao charge residual did not converge' &
    )
  call assert_true( &
    ieee_is_finite(diagnostics%minimum_field_squared_hat) .and. &
    diagnostics%minimum_field_squared_hat >= -1.0e-7_dp, &
    'nonzero-field Zhao path contains an imaginary-field interval' &
    )
  call assert_true(output(1) < 0.0_dp .and. all(output(2:3) > 0.0_dp), &
                   'nonzero Zhao-C potential/flux response has the wrong sign')
  call assert_true(all(output(4:6) == 0.0_dp), 'Zhao-C access/barrier potentials must use the upstream gauge')
  input(4:5) = [1.0e30_dp, 2.0e30_dp]
  call model%evaluate(input, inactive_output, status, message, inactive_diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'inactive ambient-outflow query failed: '//trim(message))
  call assert_true(inactive_diagnostics%branch == 'C', 'inactive ambient outflow changed the Zhao branch')
  call assert_allclose_1d( &
    inactive_output, output, 0.0_dp, 'inactive ambient outflow changed the Zhao response' &
    )
  call test_end()

  call test_begin('explicit_branch_never_falls_back')
  call initialize_model('b', configured_photoelectron_temperature_ev)
  input = 0.0_dp
  input(1) = -0.02_dp*eps0
  call model%evaluate(input, output, status, message, diagnostics)
  call assert_equal_i32( &
    status, matching_plane_zhao_no_physical_solution, &
    'explicit Zhao-B silently fell back for a negative interface field' &
    )
  call assert_true(all(output == 0.0_dp), 'failed explicit branch returned a partial response')
  call test_end()

  call test_begin('positive_photoelectron_flux_requires_positive_mean_energy')
  call initialize_model('auto', configured_photoelectron_temperature_ev)
  input = 0.0_dp
  input(2) = 1.0e10_dp
  input(3) = 0.0_dp
  call model%evaluate(input, output, status, message, diagnostics)
  call assert_equal_i32( &
    status, matching_plane_zhao_invalid_argument, &
    'positive PE flux with zero mean normal energy was accepted' &
    )
  call assert_true(all(output == 0.0_dp), 'invalid PE energy returned a partial response')
  call test_end()

  call test_summary()

contains

  subroutine initialize_model(branch, photoelectron_temperature_ev, root_selection)
    character(len=*), intent(in) :: branch
    real(dp), intent(in) :: photoelectron_temperature_ev
    character(len=*), intent(in), optional :: root_selection
    character(len=16) :: selected_root_policy

    selected_root_policy = 'require_unique'
    if (present(root_selection)) selected_root_policy = root_selection

    call model%initialize( &
      branch, selected_root_policy, ion_density_m3, electron_temperature_ev, drift_mps, drift_mps, &
      proton_mass_kg, electron_mass_kg, photoelectron_temperature_ev, status, message &
      )
    call assert_equal_i32( &
      status, matching_plane_zhao_ok, 'Zhao matching initializer failed: '//trim(message) &
      )
  end subroutine initialize_model

end program test_matching_plane_zhao

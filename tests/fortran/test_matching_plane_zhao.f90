!> Matching-plane charge-driven Zhao 応答の入力契約と再現性を検証する。
program test_matching_plane_zhao
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, pi, qe
  use bem_sheath_model_core, only: zhao_params_type, build_zhao_params, evaluate_zhao_rho_hat
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
  real(dp), parameter :: type_a_ambient_density_m3 = 7.819215729579456e6_dp
  type(matching_plane_zhao_model_type) :: model
  type(matching_plane_zhao_diagnostics_type) :: diagnostics, repeated_diagnostics
  real(dp) :: input(5), type_a_input(5), output(6), repeated_output(6), feedback_scales(4)
  real(dp) :: electron_thermal_speed_mps, expected_electron_density_m3
  real(dp) :: expected_photoelectron_source_density_m3
  integer(i32) :: status
  character(len=512) :: message

  call test_init(8)

  call test_begin('zero_field_without_photoelectrons_is_degenerate_zhao_b')
  call model%initialize( &
    'auto', ion_density_m3, electron_temperature_ev, drift_mps, drift_mps, &
    proton_mass_kg, electron_mass_kg, configured_photoelectron_temperature_ev, status, message &
    )
  call assert_equal_i32(status, matching_plane_zhao_ok, 'Zhao matching initializer failed: '//trim(message))
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

  call test_begin('explicit_zhao_b_solves_a_positive_field_profile')
  call model%initialize( &
    'b', ion_density_m3, electron_temperature_ev, drift_mps, drift_mps, &
    proton_mass_kg, electron_mass_kg, configured_photoelectron_temperature_ev, status, message &
    )
  call assert_equal_i32(status, matching_plane_zhao_ok, 'explicit Zhao-B initializer failed')
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
  call make_type_a_matching_query(type_a_input)
  call model%initialize( &
    'a', ion_density_m3, electron_temperature_ev, drift_mps, drift_mps, &
    proton_mass_kg, electron_mass_kg, type_a_photoelectron_temperature_ev, status, message &
    )
  call assert_equal_i32(status, matching_plane_zhao_ok, 'explicit Zhao-A initializer failed')
  call model%evaluate(type_a_input, output, status, message, diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'positive-field Zhao-A solve failed: '//trim(message))
  call assert_true(diagnostics%branch == 'A', 'positive-field Zhao-A response changed branch')
  call assert_close_dp(output(1), type_a_phi0_v, 5.0e-5_dp, 'positive-field Zhao-A potential mismatch')
  call assert_true(output(4) < 0.0_dp .and. output(6) < 0.0_dp, &
                   'Zhao-A access and photoelectron barriers must expose the negative minimum')
  call assert_true( &
    ieee_is_finite(diagnostics%minimum_field_squared_hat) .and. &
    diagnostics%minimum_field_squared_hat >= -1.0e-7_dp, &
    'positive-field Zhao-A path contains an imaginary-field interval' &
    )
  call test_end()

  call test_begin('auto_positive_field_fails_closed_when_uniqueness_is_uncertain')
  call model%initialize( &
    'auto', ion_density_m3, electron_temperature_ev, drift_mps, drift_mps, &
    proton_mass_kg, electron_mass_kg, type_a_photoelectron_temperature_ev, status, message &
    )
  call assert_equal_i32(status, matching_plane_zhao_ok, 'auto Zhao initializer failed')
  call model%evaluate(type_a_input, output, status, message, diagnostics)
  call assert_true( &
    status == matching_plane_zhao_ambiguous_solution .or. &
    status == matching_plane_zhao_numerical_failure, &
    'auto positive-field Zhao query did not report a fail-closed uniqueness status' &
    )
  call assert_true(all(output == 0.0_dp), 'failed auto branch selection returned a partial response')
  call test_end()

  call test_begin('nonzero_charge_driven_root_is_finite_and_converged')
  input = 0.0_dp
  input(1) = -0.02_dp*eps0
  ! Ambient outward moments are intentionally inactive in Zhao v1.
  input(4:5) = [1.0e30_dp, 2.0e30_dp]
  call model%evaluate(input, output, status, message, diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'nonzero-field Zhao matching solve failed: '//trim(message))
  call assert_true(all(ieee_is_finite(output)), 'nonzero-field Zhao response is not finite')
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
  call test_end()

  call test_begin('explicit_branch_never_falls_back')
  call model%initialize( &
    'b', ion_density_m3, electron_temperature_ev, drift_mps, drift_mps, &
    proton_mass_kg, electron_mass_kg, configured_photoelectron_temperature_ev, status, message &
    )
  call assert_equal_i32(status, matching_plane_zhao_ok, 'explicit Zhao-B initializer failed')
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

  call test_begin('identical_queries_are_reproducible')
  call model%initialize( &
    'auto', ion_density_m3, electron_temperature_ev, drift_mps, drift_mps, &
    proton_mass_kg, electron_mass_kg, configured_photoelectron_temperature_ev, status, message &
    )
  call assert_equal_i32(status, matching_plane_zhao_ok, 'auto Zhao reinitializer failed')
  input = 0.0_dp
  input(2) = 1.0e10_dp
  input(3) = 3.0_dp
  call model%evaluate(input, output, status, message, diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'first deterministic Zhao query failed')
  call model%evaluate(input, repeated_output, status, message, repeated_diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'repeated deterministic Zhao query failed')
  call assert_allclose_1d(repeated_output, output, 0.0_dp, 'identical Zhao queries changed their response')
  call assert_true( &
    repeated_diagnostics%branch == diagnostics%branch .and. &
    repeated_diagnostics%nonlinear_iterations == diagnostics%nonlinear_iterations, &
    'identical Zhao queries changed root diagnostics' &
    )
  call assert_close_dp( &
    repeated_diagnostics%residual_norm, diagnostics%residual_norm, 0.0_dp, &
    'identical Zhao queries changed residual norm' &
    )
  expected_photoelectron_source_density_m3 = &
    2.0_dp*sqrt(acos(-1.0_dp))*input(2)/sqrt(2.0_dp*qe*input(3)/electron_mass_kg)
  call assert_close_dp( &
    diagnostics%effective_photoelectron_temperature_ev, input(3), 0.0_dp, &
    'PE mean normal energy was not used as the half-Maxwellian temperature' &
    )
  call assert_close_dp( &
    diagnostics%photoelectron_source_density_m3, expected_photoelectron_source_density_m3, &
    1.0e-12_dp*expected_photoelectron_source_density_m3, &
    'PE outward flux was not mapped to the half-Maxwellian source density' &
    )
  call test_end()

  call test_summary()

contains

  subroutine make_type_a_matching_query(query)
    real(dp), intent(out) :: query(5)

    integer, parameter :: integration_panels = 2048
    type(zhao_params_type) :: params
    real(dp) :: photoelectron_source_density_m3, density_hat
    real(dp) :: h, t, phi_hat, jacobian, rho_hat, weight, integral, field_squared_hat
    integer :: point

    photoelectron_source_density_m3 = 64.0e6_dp*sin(pi/3.0_dp)
    call build_zhao_params( &
      90.0_dp, ion_density_m3, ion_density_m3, electron_temperature_ev, &
      type_a_photoelectron_temperature_ev, drift_mps, drift_mps, proton_mass_kg, &
      electron_mass_kg, params, &
      photoelectron_source_scale=photoelectron_source_density_m3/ion_density_m3 &
      )
    density_hat = type_a_ambient_density_m3/params%n_phe_ref_m3
    h = 1.0_dp/real(integration_panels, dp)
    integral = 0.0_dp
    do point = 0, integration_panels
      t = real(point, dp)*h
      phi_hat = type_a_phi_m_v/params%t_phe_ev + &
                (type_a_phi0_v - type_a_phi_m_v)/params%t_phe_ev*sin(0.5_dp*pi*t)**2
      jacobian = (type_a_phi0_v - type_a_phi_m_v)/params%t_phe_ev*0.5_dp*pi*sin(pi*t)
      call evaluate_zhao_rho_hat( &
        params, 'A', 'lower', phi_hat, type_a_phi0_v/params%t_phe_ev, &
        type_a_phi_m_v/params%t_phe_ev, density_hat, rho_hat &
        )
      if (point == 0 .or. point == integration_panels) then
        weight = 1.0_dp
      else if (mod(point, 2) == 0) then
        weight = 2.0_dp
      else
        weight = 4.0_dp
      end if
      integral = integral + weight*rho_hat*jacobian
    end do
    integral = integral*h/3.0_dp
    field_squared_hat = -2.0_dp*integral
    call assert_true(field_squared_hat > 0.0_dp, 'Type-A regression fixture has no real interface field')

    query = 0.0_dp
    query(1) = eps0*(params%t_phe_ev/params%lambda_d_phe_ref_m)*sqrt(field_squared_hat)
    query(2) = photoelectron_source_density_m3*params%v_phe_th_mps/(2.0_dp*sqrt(pi))
    query(3) = type_a_photoelectron_temperature_ev
  end subroutine make_type_a_matching_query

end program test_matching_plane_zhao

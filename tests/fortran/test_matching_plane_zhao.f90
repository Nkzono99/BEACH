!> Public matching-plane contracts, independent roots, and physical scaling.
program test_matching_plane_zhao
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, qe, pi
  use bem_matching_plane_zhao, only: matching_plane_zhao_model_type, &
                                     matching_plane_zhao_diagnostics_type, matching_plane_zhao_root_seed_type, &
                                     matching_plane_zhao_ok, matching_plane_zhao_invalid_argument, &
                                     matching_plane_zhao_no_physical_solution, matching_plane_zhao_ambiguous_solution
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, &
                          assert_equal_i32, assert_close_dp, assert_allclose_1d
  implicit none

  real(dp), parameter :: me = 9.1093837015e-31_dp, mi = 1.67262192369e-27_dp
  real(dp), parameter :: ni = 8.7e6_dp, te = 12.0_dp, tpe = 2.2_dp
  real(dp), parameter :: vi = 4.0529988897111727e5_dp
  ! Independent velocity/Poisson quadrature, boundary_validation_20260923 case 5.
  real(dp), parameter :: phi_a = 3.3935508133919789_dp, minimum_a = -0.53407940085298677_dp
  real(dp), parameter :: ne_a = 9.4289391045459863e6_dp
  real(dp), parameter :: input_a(5) = [1.6_dp*eps0, 1.3754433596232731e13_dp, tpe, 0.0_dp, 0.0_dp]
  type(matching_plane_zhao_model_type) :: model
  type(matching_plane_zhao_diagnostics_type) :: diagnostics
  type(matching_plane_zhao_root_seed_type) :: seed, candidate, reconstructed
  real(dp) :: input(5), output(6), reference(6), scales(4), previous_phi
  real(dp) :: density_factor, temperature_factor, transition_vi, transition_flux
  integer(i32) :: status
  integer :: index
  character(len=512) :: message

  call test_init(11)
  call test_begin('zero_field_without_photoelectrons_is_flat_b')
  call initialize_model('auto', 'require_unique', vi)
  input = 0.0_dp
  call model%evaluate(input, output, status, message, diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'flat solve: '//trim(message))
  call assert_true(diagnostics%branch == 'B', 'flat branch')
  call assert_close_dp(output(1), 0.0_dp, 0.0_dp, 'flat potential')
  call assert_close_dp(diagnostics%ambient_electron_density_m3, &
                       2.0_dp*ni/(1.0_dp + erf(vi/sqrt(2.0_dp*qe*te/me))), 1.0e-7_dp, 'flat neutrality')
  call assert_true(all(output(4:6) == 0.0_dp), 'flat barriers')
  call model%get_feedback_scales(scales, status, message)
  call assert_true(all(scales(1:2) > 0.0_dp) .and. all(scales(3:4) == 0.0_dp), 'feedback dependencies')
  call initialize_model('auto', 'require_unique', 0.0_dp)
  call model%evaluate(input, output, status, message, diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'zero-drift flat solve: '//trim(message))
  call assert_true(diagnostics%branch == 'B', 'near-flat C candidates must coalesce with flat B')
  call test_end()

  call test_begin('zero_drift_type_a_matches_independent_orbit_root')
  call initialize_model('a', 'require_unique', 0.0_dp)
  call model%evaluate(input_a, output, status, message, diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'A solve: '//trim(message))
  call assert_close_dp(output(1), phi_a, 2.0e-5_dp, 'independent A potential')
  call assert_close_dp(output(4), minimum_a, 2.0e-5_dp, 'independent A minimum')
  call assert_close_dp(output(6), minimum_a, 2.0e-5_dp, 'A PE barrier')
  call assert_close_dp(diagnostics%ambient_electron_density_m3, ne_a, 10.0_dp, 'independent A amplitude')
  call assert_true(diagnostics%minimum_field_squared_hat >= -1.0e-7_dp, 'A real field path')
  call test_end()

  call test_begin('input_scales_preserve_dimensionless_type_a_root')
  do index = 1, 2
    density_factor = 10.0_dp**(4*index - 6)
    temperature_factor = 10.0_dp**(2*index - 3)
    call model%initialize('a', 'require_unique', ni*density_factor, te*temperature_factor, 0.0_dp, &
                          vi*sqrt(temperature_factor), mi, me, tpe*temperature_factor, status, message)
    input = input_a
    input(1) = input_a(1)*sqrt(density_factor*temperature_factor)
    input(2) = input_a(2)*density_factor*sqrt(temperature_factor)
    input(3) = input_a(3)*temperature_factor
    call model%evaluate(input, output, status, message, diagnostics)
    call assert_equal_i32(status, matching_plane_zhao_ok, 'scaled A solve: '//trim(message))
    call assert_close_dp(output(1)/temperature_factor, phi_a, 2.0e-5_dp, 'scaled potential')
    call assert_close_dp(output(4)/temperature_factor, minimum_a, 2.0e-5_dp, 'scaled minimum')
    call assert_close_dp(diagnostics%ambient_electron_density_m3/density_factor, ne_a, 10.0_dp, 'scaled density')
  end do
  call test_end()

  call test_begin('positive_drift_cannot_connect_a_or_c_to_field_free_infinity')
  call initialize_model('a', 'require_unique', vi)
  call model%evaluate(input_a, output, status, message)
  call assert_equal_i32(status, matching_plane_zhao_no_physical_solution, 'drifting A was accepted')
  call initialize_model('c', 'require_unique', vi)
  input = [-0.02_dp*eps0, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  call model%evaluate(input, output, status, message)
  call assert_equal_i32(status, matching_plane_zhao_no_physical_solution, 'drifting C was accepted')
  call test_end()

  call test_begin('zero_drift_type_b_matches_independent_orbit_root')
  call model%initialize('b', 'require_unique', ni, te, 0.0_dp, &
                        10.0_dp*sqrt(qe*te/mi), mi, me, 0.2_dp*te, status, message)
  input = [0.1_dp*sqrt(eps0*ni*qe*te), 0.3_dp*ni*sqrt(qe*te/me), 0.2_dp*te, 0.0_dp, 0.0_dp]
  call model%evaluate(input, output, status, message, diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'B solve: '//trim(message))
  call assert_close_dp(output(1), 0.0228867312916_dp*te, 1.0e-6_dp, 'independent B potential')
  call assert_close_dp(diagnostics%ambient_electron_density_m3/ni, 0.5003210848357_dp, 1.0e-6_dp, 'B amplitude')
  call assert_true(all(output(4:6) == 0.0_dp), 'B barriers')
  call test_end()

  call test_begin('accepted_continuation_seed_reuses_and_reconstructs_all_state')
  call initialize_model('a', 'continuation', 0.0_dp)
  call model%evaluate(input_a, reference, status, message, diagnostics, continuation_candidate=seed)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'continuation bootstrap: '//trim(message))
  call assert_true(seed%valid .and. seed%branch == 'A', 'bootstrap seed branch')
  input = input_a
  input(1) = input(1)*(1.0_dp + 1.0e-6_dp)
  call model%evaluate(input, output, status, message, diagnostics, seed, candidate)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'continuation: '//trim(message))
  call assert_true(diagnostics%continuation_used .and. .not. diagnostics%continuation_fallback_used, 'fast continuation')
  call model%reconstruct_seed(input, output, reconstructed, status, message)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'reconstruct A: '//trim(message))
  call assert_true(reconstructed%branch == 'A', 'reconstructed A branch')
  call assert_close_dp(reconstructed%ambient_electron_density_m3, candidate%ambient_electron_density_m3, &
                       1.0e-6_dp, 'reconstructed ambient density')
  seed = candidate
  seed%ambient_electron_density_m3 = seed%ambient_electron_density_m3*exp(10.0_dp)
  call model%evaluate(input, output, status, message, diagnostics, seed, candidate)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'distant continuation: '//trim(message))
  call assert_true(diagnostics%continuation_fallback_used, 'distant seed needs multistart')
  call assert_true(candidate%valid .and. diagnostics%continuation_root_jump > 0.25_dp, 'distant recovery receipt')
  call test_end()

  call test_begin('auto_bootstrap_does_not_hide_multiple_roots_with_energy_ranking')
  call initialize_model('auto', 'continuation', 0.0_dp)
  call model%evaluate(input_a, output, status, message, diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ambiguous_solution, 'ambiguous auto bootstrap')
  call assert_true(all(output == 0.0_dp), 'ambiguous bootstrap returned a response')
  call test_end()

  call test_begin('zero_field_nonflat_c_and_continuous_transition_to_negative_potential_a')
  transition_vi = 468.0e3_dp*sin(20.0_dp*pi/180.0_dp)
  transition_flux = 64.0e6_dp*sin(20.0_dp*pi/180.0_dp)*sqrt(2.0_dp*qe*tpe/me)/(2.0_dp*sqrt(pi))
  call model%initialize('auto', 'continuation', ni, te, 0.0_dp, transition_vi, mi, me, tpe, status, message)
  seed = matching_plane_zhao_root_seed_type()
  do index = -1, 1
    input = [real(index, dp)*0.01_dp*eps0, transition_flux, tpe, 0.0_dp, 0.0_dp]
    call model%evaluate(input, output, status, message, diagnostics, seed, candidate)
    call assert_equal_i32(status, matching_plane_zhao_ok, 'field transition: '//trim(message))
    call assert_true(output(1) < -5.0_dp, 'nonflat root was replaced by the flat state')
    if (index <= 0) call assert_true(diagnostics%branch == 'C', 'nonpositive field C')
    if (index > 0) call assert_true(diagnostics%branch == 'A', 'positive field negative-potential A')
    if (index > -1) call assert_true(abs(output(1) - previous_phi) < 0.01_dp, 'root transition discontinuity')
    if (index == 0) call assert_close_dp(output(1), -5.445785809762306_dp, 2.0e-5_dp, 'independent zero-field C')
    call model%reconstruct_seed(input, output, reconstructed, status, message)
    call assert_equal_i32(status, matching_plane_zhao_ok, 'transition seed reconstruction: '//trim(message))
    call assert_true(reconstructed%branch == candidate%branch, 'transition seed branch')
    call assert_close_dp(reconstructed%phi_m_v, candidate%phi_m_v, 1.0e-9_dp, 'global minimum in C/A seed')
    previous_phi = output(1)
    seed = reconstructed
  end do
  call test_end()

  call test_begin('negative_field_c_ignores_ambient_outward_flux_inputs')
  call initialize_model('auto', 'require_unique', 0.0_dp)
  input = [-0.02_dp*eps0, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  call model%evaluate(input, reference, status, message, diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'dark C: '//trim(message))
  call assert_true(diagnostics%branch == 'C' .and. reference(1) < 0.0_dp, 'negative field C')
  input(4:5) = [1.0e30_dp, 2.0e30_dp]
  call model%evaluate(input, output, status, message)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'inactive inputs: '//trim(message))
  call assert_allclose_1d(output, reference, 0.0_dp, 'inactive ambient outflow altered response')
  call test_end()

  call test_begin('explicit_branch_never_changes_branches')
  call initialize_model('b', 'continuation', 0.0_dp)
  call model%evaluate(input, output, status, message)
  call assert_equal_i32(status, matching_plane_zhao_no_physical_solution, 'negative field B was accepted')
  call assert_true(all(output == 0.0_dp), 'failed branch returned partial response')
  call test_end()

  call test_begin('positive_photoelectron_flux_requires_positive_mean_energy')
  call initialize_model('auto', 'require_unique', 0.0_dp)
  input = [0.0_dp, 1.0e10_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  call model%evaluate(input, output, status, message)
  call assert_equal_i32(status, matching_plane_zhao_invalid_argument, 'invalid PE mean energy accepted')
  call test_end()
  call test_summary()

contains

  subroutine initialize_model(branch, selection, electron_drift)
    character(len=*), intent(in) :: branch, selection
    real(dp), intent(in) :: electron_drift
    call model%initialize(branch, selection, ni, te, electron_drift, vi, mi, me, tpe, status, message)
    call assert_equal_i32(status, matching_plane_zhao_ok, 'initialize: '//trim(message))
  end subroutine initialize_model
end program test_matching_plane_zhao

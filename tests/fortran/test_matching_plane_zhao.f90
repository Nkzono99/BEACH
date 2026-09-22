!> Matching-plane charge-driven Zhao 応答の物理分岐と入力契約を検証する。
program test_matching_plane_zhao
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, qe, pi
  use bem_sheath_model_core, only: zhao_params_type, build_zhao_params, evaluate_zhao_density_hat
  use bem_matching_plane_zhao, only: &
    matching_plane_zhao_model_type, matching_plane_zhao_diagnostics_type, &
    matching_plane_zhao_root_seed_type, &
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
  ! Independent nested quadrature references. Type B uses SciPy adaptive
  ! quadrature of Zhao (2020), Eq. (3), with neutrality and the prescribed field.
  real(dp), parameter :: type_a_energy_reference_j_m2 = -1.2875334387049235e-11_dp
  real(dp), parameter :: type_b_energy_reference_j_m2 = -1.2958400777036397e-11_dp
  ! Type-A既知解を固定値化し、productionのrho積分でtest入力を再生成しない。
  real(dp), parameter :: type_a_input(5) = [ &
                         1.4187346568707933e-11_dp, 1.3754433596232731e13_dp, &
                         type_a_photoelectron_temperature_ev, 0.0_dp, 0.0_dp &
                         ]
  type(matching_plane_zhao_model_type) :: model
  type(matching_plane_zhao_diagnostics_type) :: diagnostics, inactive_diagnostics
  type(matching_plane_zhao_diagnostics_type) :: energy_a_diagnostics, energy_b_diagnostics
  type(matching_plane_zhao_diagnostics_type) :: threaded_energy_a_diagnostics
  type(matching_plane_zhao_root_seed_type) :: bootstrap_seed, continued_seed, reconstructed_seed
  type(matching_plane_zhao_root_seed_type) :: fallback_seed, fallback_candidate, distant_seed, rejected_seed
  real(dp) :: input(5), output(6), inactive_output(6), feedback_scales(4), energy_a_output(6)
  real(dp) :: threaded_energy_a_output(6), continuation_output(6), reconstructed_output(6)
  real(dp) :: electron_thermal_speed_mps, expected_electron_density_m3
  integer :: omp_threads_before
  integer(i32) :: status
  character(len=512) :: message

  call test_init(13)

  call test_begin('type_b_density_retains_the_original_velocity_cutoff')
  call assert_type_b_velocity_integral()
  call test_end()

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
  omp_threads_before = 1
!$ omp_threads_before = omp_get_max_threads()
!$ call omp_set_num_threads(1)
  call initialize_model('a', type_a_photoelectron_temperature_ev, 'minimum_energy')
  call model%evaluate(input, energy_a_output, status, message, energy_a_diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'minimum-energy Zhao-A solve failed: '//trim(message))
!$ call omp_set_num_threads(min(4, omp_threads_before))
  call model%evaluate(input, threaded_energy_a_output, status, message, threaded_energy_a_diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'threaded minimum-energy Zhao-A solve failed: '//trim(message))
!$ call omp_set_num_threads(omp_threads_before)
  call assert_allclose_1d( &
    threaded_energy_a_output, energy_a_output, 0.0_dp, &
    'OpenMP Zhao-A response changed from the serial response' &
    )
  call assert_true( &
    threaded_energy_a_diagnostics%branch == energy_a_diagnostics%branch, &
    'OpenMP Zhao-A branch changed from the serial branch' &
    )
  call assert_close_dp( &
    threaded_energy_a_diagnostics%potential_energy_j_m2, &
    energy_a_diagnostics%potential_energy_j_m2, 0.0_dp, &
    'OpenMP Zhao-A energy changed from the serial energy' &
    )
  call assert_close_dp( &
    threaded_energy_a_diagnostics%minimum_field_squared_hat, &
    energy_a_diagnostics%minimum_field_squared_hat, 0.0_dp, &
    'OpenMP Zhao-A profile changed from the serial profile' &
    )
  call initialize_model('b', type_a_photoelectron_temperature_ev, 'minimum_energy')
  call model%evaluate(input, output, status, message, energy_b_diagnostics)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'minimum-energy Zhao-B solve failed: '//trim(message))
  call assert_close_dp(output(1), 3.956627756217795_dp, 1.e-7_dp, 'independent Type B potential')
  call assert_close_dp(energy_b_diagnostics%ambient_electron_density_m3, 6742539.857612752_dp, &
                       1._dp, 'independent Type B ambient density')
  call assert_close_dp( &
    energy_a_diagnostics%potential_energy_j_m2, type_a_energy_reference_j_m2, &
    1.0e-7_dp*abs(type_a_energy_reference_j_m2), &
    'Zhao-A energy quadrature changed the independent reference energy' &
    )
  call assert_close_dp( &
    energy_b_diagnostics%potential_energy_j_m2, type_b_energy_reference_j_m2, &
    1.0e-7_dp*abs(type_b_energy_reference_j_m2), &
    'Zhao-B energy quadrature changed the independent reference energy' &
    )
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
  ! Independent handoff fixture: M=10, Tph/Te=0.2, G=0.3, E_H=0.1, zero electron drift.
  call model%initialize('b', 'require_unique', ion_density_m3, electron_temperature_ev, 0._dp, &
                        10*sqrt(qe*electron_temperature_ev/proton_mass_kg), proton_mass_kg, electron_mass_kg, &
                        0.2_dp*electron_temperature_ev, status, message)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'zero-drift B initialization')
  input = [0.1_dp*sqrt(eps0*ion_density_m3*qe*electron_temperature_ev), &
           0.3_dp*ion_density_m3*sqrt(qe*electron_temperature_ev/electron_mass_kg), &
           0.2_dp*electron_temperature_ev, 0._dp, 0._dp]
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
  call assert_close_dp(output(1), 0.0228867312916_dp*electron_temperature_ev, 1.e-7_dp, &
                       'zero-drift B differs from the independent orbit solver')
  call assert_close_dp(diagnostics%ambient_electron_density_m3/ion_density_m3, 0.5003210848357_dp, 1.e-7_dp, &
                       'zero-drift B amplitude differs from the independent orbit solver')
  call test_end()

  call test_begin('positive_field_without_photoelectrons_has_no_type_b_response')
  call initialize_model('b', configured_photoelectron_temperature_ev)
  input = 0._dp
  input(1) = 0.02_dp*eps0
  call model%evaluate(input, output, status, message, diagnostics)
  call assert_true(status == matching_plane_zhao_no_physical_solution .or. &
                   status == matching_plane_zhao_numerical_failure, 'dark positive-field B was accepted')
  call assert_true(all(output == 0._dp), 'dark positive-field B returned a partial response')
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

  call test_begin('type_a_continuation_reuses_and_reconstructs_the_accepted_root')
  call initialize_model('a', type_a_photoelectron_temperature_ev, 'continuation')
  call model%evaluate( &
    type_a_input, output, status, message, diagnostics, continuation_candidate=bootstrap_seed &
    )
  call assert_equal_i32(status, matching_plane_zhao_ok, 'Type-A continuation bootstrap failed: '//trim(message))
  call assert_true(bootstrap_seed%valid, 'Type-A continuation bootstrap did not return a root seed')
  call assert_true(.not. diagnostics%continuation_used, 'bootstrap was incorrectly reported as a continuation step')

  input = type_a_input
  input(1) = input(1)*(1.0_dp + 1.0e-6_dp)
  call model%evaluate( &
    input, continuation_output, status, message, diagnostics, &
    continuation_seed=bootstrap_seed, continuation_candidate=continued_seed &
    )
  call assert_equal_i32(status, matching_plane_zhao_ok, 'Type-A continuation step failed: '//trim(message))
  call assert_true(diagnostics%continuation_used, 'accepted seed was not used for Type-A continuation')
  call assert_true(.not. diagnostics%continuation_fallback_used, 'nearby Type-A root required a full multistart fallback')
  call assert_true( &
    diagnostics%continuation_root_jump >= 0.0_dp .and. diagnostics%continuation_root_jump <= 0.25_dp, &
    'accepted Type-A continuation root exceeded the bounded encoded distance' &
    )
  call assert_true(continued_seed%valid, 'accepted Type-A continuation root did not return the next seed')

  call model%reconstruct_seed(input, continuation_output, reconstructed_seed, status, message)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'Type-A restart seed reconstruction failed: '//trim(message))
  call assert_close_dp( &
    reconstructed_seed%ambient_electron_density_m3, continued_seed%ambient_electron_density_m3, &
    1.0e-12_dp*continued_seed%ambient_electron_density_m3, &
    'restart response did not reconstruct the accepted ambient electron density' &
    )
  call model%evaluate( &
    input, reconstructed_output, status, message, diagnostics, continuation_seed=reconstructed_seed &
    )
  call assert_equal_i32(status, matching_plane_zhao_ok, 'reconstructed Type-A seed was not reusable: '//trim(message))
  call assert_allclose_1d( &
    reconstructed_output, continuation_output, 1.0e-12_dp*maxval(abs(continuation_output)), &
    'reconstructed Type-A seed changed the matching response' &
    )
  call test_end()

  call test_begin('type_a_continuation_reacquires_a_nearby_root_after_local_newton_failure')
  input = [5.6749386274831732e-12_dp, 4.8140517586819559e12_dp, 1.1_dp, 0.0_dp, 0.0_dp]
  fallback_seed%valid = .true.
  fallback_seed%phi0_v = 4.1874671269993607e-1_dp
  fallback_seed%phi_m_v = -6.4058530638434064e-1_dp
  fallback_seed%ambient_electron_density_m3 = 8.5523788611962516e6_dp
  call model%evaluate( &
    input, output, status, message, diagnostics, &
    continuation_seed=fallback_seed, continuation_candidate=fallback_candidate &
    )
  call assert_equal_i32(status, matching_plane_zhao_ok, 'Type-A continuation fallback failed: '//trim(message))
  call assert_true(diagnostics%continuation_used, 'Type-A continuation did not use the supplied seed')
  call assert_true(diagnostics%continuation_fallback_used, 'local Newton failure did not trigger full multistart')
  call assert_close_dp(diagnostics%continuation_root_jump, 0.2_dp, 1.0e-10_dp, 'fallback root distance mismatch')
  call assert_true(fallback_candidate%valid, 'reacquired nearby Type-A root did not return a candidate seed')
  call assert_close_dp(output(1), 3.7889769433045589e-1_dp, 1.0e-10_dp, 'fallback matching potential mismatch')
  call assert_close_dp(output(4), -7.8241266005471111e-1_dp, 1.0e-10_dp, 'fallback potential minimum mismatch')
  call test_end()

  call test_begin('type_a_continuation_reacquires_a_unique_distant_root_with_multistart')
  input = type_a_input
  distant_seed = continued_seed
  distant_seed%ambient_electron_density_m3 = &
    distant_seed%ambient_electron_density_m3*exp(10.0_dp)
  call model%evaluate( &
    input, output, status, message, diagnostics, &
    continuation_seed=distant_seed, continuation_candidate=rejected_seed &
    )
  call assert_equal_i32(status, matching_plane_zhao_ok, 'Type-A continuation did not reacquire a distant root')
  call assert_true(diagnostics%continuation_fallback_used, 'distant root did not trigger full multistart recovery')
  call assert_true(rejected_seed%valid, 'reacquired distant root did not return a continuation seed')
  call assert_true( &
    diagnostics%continuation_root_jump > 0.25_dp, &
    'full multistart recovery did not record the large encoded root distance' &
    )
  call assert_close_dp(output(1), type_a_phi0_v, 5.0e-5_dp, 'distant fallback selected the wrong Type-A root')
  call assert_close_dp(output(4), type_a_phi_m_v, 5.0e-5_dp, 'distant fallback changed the Type-A minimum')
  call assert_close_dp(output(6), type_a_phi_m_v, 5.0e-5_dp, 'distant fallback changed the PE barrier')
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

  subroutine assert_type_b_velocity_integral()
    type(zhao_params_type) :: params
    real(dp), parameter :: potentials(6) = [0._dp, 0.01_dp, 0.2_dp, 1._dp, 9._dp, 1000._dp]
    real(dp), parameter :: amplitude = 0.7_dp
    integer, parameter :: panels = 20000
    real(dp) :: phi, cutoff, u, v, step, weight, integral, ion, free, reflected, photo, captured
    integer :: drift_index, potential_index, point

    call build_zhao_params(90._dp, ion_density_m3, ion_density_m3, 12._dp, 2.4_dp, 0._dp, &
                           100*sqrt(qe*12._dp/proton_mass_kg), proton_mass_kg, electron_mass_kg, params, &
                           photoelectron_source_scale=0._dp)
    do drift_index = 0, 2
      u = 0.5_dp*drift_index
      params%u = u
      do potential_index = 1, size(potentials)
        phi = potentials(potential_index)
        cutoff = sqrt(phi)
        step = 12._dp/(max(1._dp, (cutoff - u)/4)*panels)
        integral = 0
        ! Direct quadrature of the original local Maxwellian over v >= sqrt(phi).
        ! Velocity is in sqrt(2 Te/me) units; no erf/erfc is used by this oracle.
        do point = 0, panels
          v = cutoff + point*step
          weight = 2
          if (mod(point, 2) == 1) weight = 4
          if (point == 0 .or. point == panels) weight = 1
          integral = integral + weight*exp(phi - (v - u)**2)
        end do
        integral = amplitude*integral*step/(3*sqrt(pi))
        call evaluate_zhao_density_hat(params, 'B', 'monotonic', phi*params%tau, phi*params%tau, &
                                       phi*params%tau, amplitude, ion, free, reflected, photo, captured)
        call assert_true(ieee_is_finite(free), 'Type B electron density overflowed')
        call assert_close_dp(free, integral, 3.e-11_dp*integral, 'Type B velocity-domain integral mismatch')
        call assert_close_dp(reflected, 0._dp, 0._dp, 'Type B acquired reflected ambient electrons')
      end do
    end do
  end subroutine assert_type_b_velocity_integral

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

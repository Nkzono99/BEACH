program test_mean_charge_integrator
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, qe
  use bem_mean_charge_integrator, only: &
    mean_charge_step_ok, mean_charge_step_result_type, &
    planar_debye_capacitance_per_area, advance_mean_charge_backward_euler
  use bem_outer_plasma_kinetic, only: kinetic_outer_plasma_options_type
  use bem_dynamic_k0_mean, only: &
    dynamic_k0_ok, dynamic_k0_invalid, dynamic_k0_step_type, advance_dynamic_k0_mean
  use test_support, only: &
    test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, assert_equal_i32
  implicit none

  real(dp), parameter :: debye_length_m = 2.5_dp
  real(dp), parameter :: photo_current_a_m2 = 4.0_dp
  real(dp), parameter :: photo_temperature_v = 3.0_dp
  real(dp), parameter :: no_sink_a_m2 = 0.0_dp
  real(dp) :: capacitance_per_area, surface_charge_density, time_step, final_time
  real(dp) :: expected_potential, equilibrium_potential, explicit_potential
  real(dp) :: parameters(3)
  type(mean_charge_step_result_type) :: result
  type(kinetic_outer_plasma_options_type) :: kinetic_options
  type(dynamic_k0_step_type) :: k0_step
  integer(i32) :: step
  character(len=256) :: message

  call test_init(10)

  call test_begin('planar Debye capacitance maps mean charge to interface potential')
  capacitance_per_area = planar_debye_capacitance_per_area(debye_length_m)
  call assert_close_dp( &
    capacitance_per_area, eps0/debye_length_m, 1.0e-24_dp, &
    'planar Debye capacitance mismatch' &
    )
  call test_end()

  call test_begin('backward Euler solves the scalar current root')
  parameters = [6.0_dp, 2.0_dp, 1.0_dp]
  capacitance_per_area = 2.0_dp
  surface_charge_density = 0.0_dp
  call advance_mean_charge_backward_euler( &
    surface_charge_density, 5.0_dp, capacitance_per_area, &
    photo_with_constant_sink, parameters, result &
    )
  call assert_equal_i32(result%status, mean_charge_step_ok, 'backward Euler root solve failed')
  call assert_true(result%iterations > 0_i32, 'nonlinear backward Euler step did not iterate')
  call assert_true(result%potential_v > 0.0_dp, 'photoemission must raise the mean potential')
  equilibrium_potential = parameters(2)*log(parameters(1)/parameters(3))
  call assert_true( &
    result%potential_v < equilibrium_potential, &
    'implicit step crossed the stable photo-current equilibrium' &
    )
  call assert_true( &
    abs(result%backward_euler_residual_v) < 1.0e-10_dp, &
    'backward Euler voltage residual is too large' &
    )
  call assert_close_dp( &
    result%surface_charge_density_c_m2, 5.0_dp*result%current_density_a_m2, 2.0e-10_dp, &
    'backward Euler step does not conserve mean charge' &
    )
  call assert_close_dp( &
    result%surface_charge_density_c_m2, capacitance_per_area*result%potential_v, 1.0e-12_dp, &
    'mean charge and potential lost their capacitance relation' &
    )
  call test_end()

  call test_begin('analytic photo-current equilibrium is stationary')
  parameters = [6.0_dp, 2.0_dp, 1.0_dp]
  capacitance_per_area = 2.0_dp
  equilibrium_potential = parameters(2)*log(parameters(1)/parameters(3))
  surface_charge_density = capacitance_per_area*equilibrium_potential
  call advance_mean_charge_backward_euler( &
    surface_charge_density, 1.0e9_dp, capacitance_per_area, &
    photo_with_constant_sink, parameters, result &
    )
  call assert_equal_i32(result%status, mean_charge_step_ok, 'stationary mean charge step failed')
  call assert_close_dp( &
    result%potential_v, equilibrium_potential, 1.0e-11_dp, &
    'analytic photo-current equilibrium drifted' &
    )
  call test_end()

  call test_begin('UV-only mean potential follows the logarithmic transient')
  parameters = [photo_current_a_m2, photo_temperature_v, no_sink_a_m2]
  capacitance_per_area = 2.0_dp
  surface_charge_density = 0.0_dp
  final_time = 30.0_dp
  time_step = final_time/1000.0_dp
  do step = 1_i32, 1000_i32
    call advance_mean_charge_backward_euler( &
      surface_charge_density, time_step, capacitance_per_area, &
      photo_with_constant_sink, parameters, result &
      )
    call assert_equal_i32(result%status, mean_charge_step_ok, 'UV-only implicit step failed')
    surface_charge_density = result%surface_charge_density_c_m2
  end do
  expected_potential = photo_temperature_v*log( &
                       1.0_dp + photo_current_a_m2*final_time/ &
                       (capacitance_per_area*photo_temperature_v) &
                       )
  call assert_close_dp( &
    result%potential_v, expected_potential, 6.0e-3_dp, &
    'UV-only implicit solution lost the logarithmic charging transient' &
    )
  call test_end()

  call test_begin('implicit photo barrier avoids explicit batch overshoot')
  parameters = [1.0_dp, 1.0_dp, 0.0_dp]
  capacitance_per_area = 1.0_dp
  surface_charge_density = 0.0_dp
  time_step = 1000.0_dp
  explicit_potential = time_step*parameters(1)/capacitance_per_area
  call advance_mean_charge_backward_euler( &
    surface_charge_density, time_step, capacitance_per_area, &
    photo_with_constant_sink, parameters, result &
    )
  call assert_equal_i32(result%status, mean_charge_step_ok, 'large implicit photo step failed')
  call assert_close_dp( &
    result%explicit_predictor_potential_v, explicit_potential, 1.0e-12_dp, &
    'explicit predictor diagnostic mismatch' &
    )
  call assert_true(result%potential_v > 5.0_dp .and. result%potential_v < 6.0_dp, &
                   'implicit photo-barrier root left its analytic bracket')
  call assert_true( &
    result%explicit_predictor_potential_v > 100.0_dp*result%potential_v, &
    'fixture no longer demonstrates explicit batch overshoot' &
    )
  call test_end()

  call test_begin('runtime adapter closes planar surface current implicitly')
  kinetic_options = kinetic_outer_plasma_options_type()
  kinetic_options%kinetic_closure = 'ambient_linear_debye'
  kinetic_options%tail_length = 2.0_dp
  kinetic_options%electron_charge = -qe
  kinetic_options%electron_mass = 9.1093837139e-31_dp
  kinetic_options%electron_density_infinity = 1.0_dp
  kinetic_options%electron_temperature_j = 10.0_dp*qe
  kinetic_options%ion_charge = qe
  kinetic_options%ion_mass = 1.67262192595e-27_dp
  kinetic_options%ion_density_infinity = 1.0_dp
  kinetic_options%ion_temperature_j = 1.0_dp*qe
  kinetic_options%ion_drift_infinity = 1.0e5_dp
  kinetic_options%photoelectron_charge = -qe
  kinetic_options%photoelectron_mass = kinetic_options%electron_mass
  kinetic_options%photoelectron_temperature_j = 2.0_dp*qe
  kinetic_options%photoelectron_emission_flux = 1.0e10_dp
  call advance_dynamic_k0_mean( &
    kinetic_options, 'e_bottom_zero', 1.0_dp, 0.0_dp, 1.0e-3_dp, k0_step, message, &
    0.0_dp, 0.0_dp, -kinetic_options%photoelectron_charge*kinetic_options%photoelectron_emission_flux, 0.0_dp &
    )
  call assert_equal_i32(k0_step%status, dynamic_k0_ok, 'runtime dynamic k0 step failed: '//trim(message))
  call assert_close_dp( &
    k0_step%surface_charge_after_c - k0_step%surface_charge_before_c, &
    1.0e-3_dp*k0_step%total_current_density_a_m2, 1.0e-16_dp, &
    'runtime dynamic k0 step does not close its backward-Euler current' &
    )
  call assert_close_dp( &
    k0_step%interface_potential_after_v, &
    kinetic_options%tail_length*k0_step%interface_field_after_v_m, 1.0e-14_dp, &
    'runtime dynamic k0 potential and field are inconsistent' &
    )
  call test_end()

  call test_begin('runtime adapter keeps tracked ambient current explicit')
  call advance_dynamic_k0_mean( &
    kinetic_options, 'e_bottom_zero', 1.0_dp, 0.0_dp, 1.0e-3_dp, k0_step, message, &
    tracked_electron_current_density_a_m2=-2.0e-9_dp, &
    tracked_ion_current_density_a_m2=0.5e-9_dp, &
    tracked_photoelectron_outward_current_density_a_m2= &
    -kinetic_options%photoelectron_charge*kinetic_options%photoelectron_emission_flux, &
    tracked_additional_current_density_a_m2=0.0_dp &
    )
  call assert_equal_i32(k0_step%status, dynamic_k0_ok, 'tracked-current dynamic k0 step failed: '//trim(message))
  call assert_close_dp(k0_step%electron_current_density_a_m2, -2.0e-9_dp, 0.0_dp, &
                       'tracked electron current was replaced')
  call assert_close_dp(k0_step%ion_current_density_a_m2, 0.5e-9_dp, 0.0_dp, &
                       'tracked ion current was replaced')
  call assert_close_dp( &
    k0_step%surface_charge_after_c - k0_step%surface_charge_before_c, &
    1.0e-3_dp*k0_step%total_current_density_a_m2, 1.0e-16_dp, &
    'tracked-current IMEX step does not close its mean charge' &
    )
  call test_end()

  call test_begin('runtime adapter uses tracked photoelectron outward-current amplitude')
  call advance_dynamic_k0_mean( &
    kinetic_options, 'e_bottom_zero', 1.0_dp, 0.0_dp, 1.0e-3_dp, k0_step, message, &
    tracked_electron_current_density_a_m2=0.0_dp, &
    tracked_ion_current_density_a_m2=0.0_dp, &
    tracked_photoelectron_outward_current_density_a_m2=8.0e-9_dp, &
    tracked_additional_current_density_a_m2=0.0_dp &
    )
  call assert_equal_i32(k0_step%status, dynamic_k0_ok, &
                        'tracked-photoelectron dynamic k0 step failed: '//trim(message))
  call assert_close_dp( &
    k0_step%photoelectron_current_density_a_m2, &
    8.0e-9_dp*k0_step%photoelectron_escape_fraction, 1.0e-18_dp, &
    'tracked photoelectron outward current did not replace the configured emission amplitude' &
    )
  kinetic_options%photoelectron_charge = -0.5_dp*qe
  call advance_dynamic_k0_mean( &
    kinetic_options, 'e_bottom_zero', 1.0_dp, 0.0_dp, 1.0e-3_dp, k0_step, message, &
    tracked_electron_current_density_a_m2=0.0_dp, &
    tracked_ion_current_density_a_m2=0.0_dp, &
    tracked_photoelectron_outward_current_density_a_m2=8.0e-9_dp, &
    tracked_additional_current_density_a_m2=0.0_dp &
    )
  call assert_equal_i32(k0_step%status, dynamic_k0_ok, &
                        'tracked-photoelectron charge-independence step failed: '//trim(message))
  call assert_close_dp( &
    k0_step%photoelectron_current_density_a_m2, &
    8.0e-9_dp*k0_step%photoelectron_escape_fraction, 1.0e-18_dp, &
    'tracked photoelectron current amplitude depends on the ambient electron charge' &
    )
  kinetic_options%photoelectron_charge = -qe
  call advance_dynamic_k0_mean( &
    kinetic_options, 'e_bottom_zero', 1.0_dp, 0.0_dp, 1.0e-3_dp, k0_step, message, &
    tracked_electron_current_density_a_m2=0.0_dp, &
    tracked_ion_current_density_a_m2=0.0_dp, &
    tracked_photoelectron_outward_current_density_a_m2=0.0_dp, &
    tracked_additional_current_density_a_m2=0.0_dp &
    )
  call assert_equal_i32(k0_step%status, dynamic_k0_ok, &
                        'zero tracked-photoelectron dynamic k0 step failed: '//trim(message))
  call assert_close_dp(k0_step%photoelectron_current_density_a_m2, 0.0_dp, 0.0_dp, &
                       'zero tracked photoelectron outward current did not suppress the analytic source')
  call test_end()

  call test_begin('runtime adapter preserves aggregate non-closure current')
  call advance_dynamic_k0_mean( &
    kinetic_options, 'e_bottom_zero', 1.0_dp, 0.0_dp, 1.0e-3_dp, k0_step, message, &
    tracked_electron_current_density_a_m2=0.0_dp, &
    tracked_ion_current_density_a_m2=0.0_dp, &
    tracked_photoelectron_outward_current_density_a_m2=0.0_dp, &
    tracked_additional_current_density_a_m2=3.0e-9_dp &
    )
  call assert_equal_i32(k0_step%status, dynamic_k0_ok, &
                        'aggregate-current dynamic k0 step failed: '//trim(message))
  call assert_close_dp( &
    k0_step%additional_tracked_current_density_a_m2, 3.0e-9_dp, 0.0_dp, &
    'aggregate non-closure current was not retained' &
    )
  call assert_close_dp( &
    k0_step%total_current_density_a_m2, 3.0e-9_dp, 1.0e-20_dp, &
    'aggregate non-closure current was not included in total current' &
    )
  call assert_close_dp( &
    k0_step%surface_charge_after_c, 3.0e-12_dp, 1.0e-22_dp, &
    'aggregate non-closure current did not update mean surface charge' &
    )
  call test_end()

  call test_begin('runtime adapter rejects invalid tracked photoelectron current metadata')
  call advance_dynamic_k0_mean( &
    kinetic_options, 'e_bottom_zero', 1.0_dp, 0.0_dp, 1.0e-3_dp, k0_step, message, &
    tracked_electron_current_density_a_m2=0.0_dp, &
    tracked_ion_current_density_a_m2=0.0_dp, &
    tracked_photoelectron_outward_current_density_a_m2=-1.0_dp, &
    tracked_additional_current_density_a_m2=0.0_dp &
    )
  call assert_equal_i32(k0_step%status, dynamic_k0_invalid, &
                        'negative tracked photoelectron outward current was accepted')
  call advance_dynamic_k0_mean( &
    kinetic_options, 'e_bottom_zero', 1.0_dp, 0.0_dp, 1.0e-3_dp, k0_step, message, &
    tracked_electron_current_density_a_m2=0.0_dp, &
    tracked_ion_current_density_a_m2=0.0_dp, &
    tracked_photoelectron_outward_current_density_a_m2=ieee_value(0.0_dp, ieee_quiet_nan), &
    tracked_additional_current_density_a_m2=0.0_dp &
    )
  call assert_equal_i32(k0_step%status, dynamic_k0_invalid, &
                        'non-finite tracked photoelectron outward current was accepted')
  call advance_dynamic_k0_mean( &
    kinetic_options, 'e_bottom_zero', 1.0_dp, 0.0_dp, 1.0e-3_dp, k0_step, message, &
    tracked_electron_current_density_a_m2=0.0_dp, &
    tracked_ion_current_density_a_m2=0.0_dp, &
    tracked_photoelectron_outward_current_density_a_m2=0.0_dp, &
    tracked_additional_current_density_a_m2=ieee_value(0.0_dp, ieee_quiet_nan) &
    )
  call assert_equal_i32(k0_step%status, dynamic_k0_invalid, &
                        'non-finite aggregate tracked current was accepted')
  kinetic_options%photoelectron_charge = qe
  call advance_dynamic_k0_mean( &
    kinetic_options, 'e_bottom_zero', 1.0_dp, 0.0_dp, 1.0e-3_dp, k0_step, message, &
    tracked_electron_current_density_a_m2=0.0_dp, &
    tracked_ion_current_density_a_m2=0.0_dp, &
    tracked_photoelectron_outward_current_density_a_m2=1.0_dp, &
    tracked_additional_current_density_a_m2=0.0_dp &
    )
  call assert_equal_i32(k0_step%status, dynamic_k0_invalid, &
                        'tracked photoelectron outward current accepted nonnegative photoelectron charge')
  kinetic_options%photoelectron_charge = -qe
  call test_end()

  call test_summary()

contains

  pure function photo_with_constant_sink(potential_v, model_parameters) result(current_density_a_m2)
    real(dp), intent(in) :: potential_v
    real(dp), intent(in) :: model_parameters(:)
    real(dp) :: current_density_a_m2

    if (size(model_parameters) /= 3 .or. model_parameters(2) <= 0.0_dp) then
      current_density_a_m2 = 0.0_dp
      return
    end if
    current_density_a_m2 = model_parameters(1)* &
                           exp(-max(potential_v, 0.0_dp)/model_parameters(2)) - &
                           model_parameters(3)
  end function photo_with_constant_sink

end program test_mean_charge_integrator

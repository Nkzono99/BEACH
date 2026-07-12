program test_outer_plasma_kinetic_core
  use bem_kinds, only: dp, i32
  use bem_constants, only: qe
  use bem_outer_plasma_types, only: outer_plasma_ok, outer_plasma_no_physical_solution
  use bem_outer_plasma_kinetic, only: &
    eval_absorbing_maxwellian_density, eval_cold_drift_density, kinetic_bohm_speed, &
    eval_photoelectron_escape_return, eval_emitted_maxwellian_density, eval_kinetic_current_balance
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, assert_equal_i32
  implicit none

  real(dp), parameter :: ev = qe
  real(dp), parameter :: electron_mass = 9.1093837139e-31_dp
  real(dp) :: density, susceptibility, speed, escape_fraction, return_fraction
  real(dp) :: density_plus, density_minus, derivative_local, derivative_interface, derivative_fd
  real(dp) :: electron_current, ion_current, photo_current, total_current
  real(dp), parameter :: derivative_step = 1.0e-5_dp
  integer(i32) :: status

  call test_init(10)

  call test_begin('absorbing Maxwellian is normalized at infinity')
  call eval_absorbing_maxwellian_density( &
    phi=0.0_dp, phi_interface=-4.0_dp, charge=-ev, temperature_j=2.0_dp*ev, &
    density_infinity=3.0e6_dp, density=density, susceptibility=susceptibility, status=status &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'electron closure status mismatch')
  call assert_close_dp(density, 3.0e6_dp, 1.0e-8_dp, 'infinity density normalization mismatch')
  call test_end()

  call test_begin('absorbing Maxwellian analytic derivatives match finite differences')
  call eval_absorbing_maxwellian_density( &
    phi=-2.0_dp, phi_interface=-4.0_dp, charge=-ev, temperature_j=2.0_dp*ev, &
    density_infinity=3.0e6_dp, density=density, susceptibility=susceptibility, status=status, &
    derivative_interface=derivative_interface &
    )
  derivative_local = susceptibility
  call eval_absorbing_maxwellian_density( &
    phi=-2.0_dp + derivative_step, phi_interface=-4.0_dp, charge=-ev, temperature_j=2.0_dp*ev, &
    density_infinity=3.0e6_dp, density=density_plus, susceptibility=susceptibility, status=status &
    )
  call eval_absorbing_maxwellian_density( &
    phi=-2.0_dp - derivative_step, phi_interface=-4.0_dp, charge=-ev, temperature_j=2.0_dp*ev, &
    density_infinity=3.0e6_dp, density=density_minus, susceptibility=susceptibility, status=status &
    )
  derivative_fd = (density_plus - density_minus)/(2.0_dp*derivative_step)
  call assert_true( &
    abs(derivative_local - derivative_fd) <= 1.0e-8_dp*max(1.0_dp, abs(derivative_fd)), &
    'absorbing local derivative mismatch' &
    )
  call eval_absorbing_maxwellian_density( &
    phi=-2.0_dp, phi_interface=-4.0_dp + derivative_step, charge=-ev, temperature_j=2.0_dp*ev, &
    density_infinity=3.0e6_dp, density=density_plus, susceptibility=susceptibility, status=status &
    )
  call eval_absorbing_maxwellian_density( &
    phi=-2.0_dp, phi_interface=-4.0_dp - derivative_step, charge=-ev, temperature_j=2.0_dp*ev, &
    density_infinity=3.0e6_dp, density=density_minus, susceptibility=susceptibility, status=status &
    )
  derivative_fd = (density_plus - density_minus)/(2.0_dp*derivative_step)
  call assert_true( &
    abs(derivative_interface - derivative_fd) <= 1.0e-8_dp*max(1.0_dp, abs(derivative_fd)), &
    'absorbing interface derivative mismatch' &
    )
  call test_end()

  call test_begin('emitted Maxwellian counts outgoing and returning populations once')
  call eval_emitted_maxwellian_density( &
    phi=3.0_dp, phi_interface=3.0_dp, phi_infinity=0.0_dp, charge=-ev, mass=electron_mass, &
    temperature_j=1.5_dp*ev, emission_flux=1.0e10_dp, density=density, status=status &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'emitted density status mismatch')
  call assert_close_dp( &
    density, 1.0e10_dp*sqrt(acos(-1.0_dp)*electron_mass/(3.0_dp*ev))* &
    (1.0_dp + erf(sqrt(2.0_dp))), 1.0e-6_dp, 'interface emitted density mismatch' &
    )
  call eval_emitted_maxwellian_density( &
    phi=0.0_dp, phi_interface=3.0_dp, phi_infinity=0.0_dp, charge=-ev, mass=electron_mass, &
    temperature_j=1.5_dp*ev, emission_flux=1.0e10_dp, density=density, status=status &
    )
  call assert_close_dp( &
    density, 1.0e10_dp*sqrt(acos(-1.0_dp)*electron_mass/(3.0_dp*ev))*exp(-2.0_dp), &
    1.0e-6_dp, 'escaping emitted density mismatch' &
    )
  call test_end()

  call test_begin('emitted Maxwellian analytic derivatives match finite differences')
  call eval_emitted_maxwellian_density( &
    phi=1.5_dp, phi_interface=3.0_dp, phi_infinity=0.0_dp, charge=-ev, mass=electron_mass, &
    temperature_j=1.5_dp*ev, emission_flux=1.0e10_dp, density=density, status=status, &
    derivative_local=derivative_local, derivative_interface=derivative_interface &
    )
  call eval_emitted_maxwellian_density( &
    phi=1.5_dp + derivative_step, phi_interface=3.0_dp, phi_infinity=0.0_dp, charge=-ev, mass=electron_mass, &
    temperature_j=1.5_dp*ev, emission_flux=1.0e10_dp, density=density_plus, status=status &
    )
  call eval_emitted_maxwellian_density( &
    phi=1.5_dp - derivative_step, phi_interface=3.0_dp, phi_infinity=0.0_dp, charge=-ev, mass=electron_mass, &
    temperature_j=1.5_dp*ev, emission_flux=1.0e10_dp, density=density_minus, status=status &
    )
  derivative_fd = (density_plus - density_minus)/(2.0_dp*derivative_step)
  call assert_true( &
    abs(derivative_local - derivative_fd) <= 1.0e-8_dp*max(1.0_dp, abs(derivative_fd)), &
    'emitted local derivative mismatch' &
    )
  call eval_emitted_maxwellian_density( &
    phi=1.5_dp, phi_interface=3.0_dp + derivative_step, phi_infinity=0.0_dp, charge=-ev, mass=electron_mass, &
    temperature_j=1.5_dp*ev, emission_flux=1.0e10_dp, density=density_plus, status=status &
    )
  call eval_emitted_maxwellian_density( &
    phi=1.5_dp, phi_interface=3.0_dp - derivative_step, phi_infinity=0.0_dp, charge=-ev, mass=electron_mass, &
    temperature_j=1.5_dp*ev, emission_flux=1.0e10_dp, density=density_minus, status=status &
    )
  derivative_fd = (density_plus - density_minus)/(2.0_dp*derivative_step)
  call assert_true( &
    abs(derivative_interface - derivative_fd) <= 1.0e-8_dp*max(1.0_dp, abs(derivative_fd)), &
    'emitted interface derivative mismatch' &
    )
  call test_end()

  call test_begin('absorbing Maxwellian includes a finite loss cone')
  call eval_absorbing_maxwellian_density( &
    phi=-4.0_dp, phi_interface=-4.0_dp, charge=-ev, temperature_j=2.0_dp*ev, &
    density_infinity=3.0e6_dp, density=density, susceptibility=susceptibility, status=status &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'interface electron closure status mismatch')
  call assert_close_dp( &
    density, 3.0e6_dp*exp(-2.0_dp)/(1.0_dp + erf(sqrt(2.0_dp))), 1.0e-8_dp, &
    'absorbing-interface density mismatch' &
    )
  call test_end()

  call test_begin('cold drifting ion analytic derivative matches finite differences')
  call eval_cold_drift_density( &
    phi=-3.0_dp, charge=ev, mass=1836.0_dp*electron_mass, density_infinity=1.0e6_dp, &
    drift_infinity=2.0e4_dp, density=density, speed=speed, susceptibility=derivative_local, status=status &
    )
  call eval_cold_drift_density( &
    phi=-3.0_dp + derivative_step, charge=ev, mass=1836.0_dp*electron_mass, density_infinity=1.0e6_dp, &
    drift_infinity=2.0e4_dp, density=density_plus, speed=speed, susceptibility=susceptibility, status=status &
    )
  call eval_cold_drift_density( &
    phi=-3.0_dp - derivative_step, charge=ev, mass=1836.0_dp*electron_mass, density_infinity=1.0e6_dp, &
    drift_infinity=2.0e4_dp, density=density_minus, speed=speed, susceptibility=susceptibility, status=status &
    )
  derivative_fd = (density_plus - density_minus)/(2.0_dp*derivative_step)
  call assert_true( &
    abs(derivative_local - derivative_fd) <= 1.0e-8_dp*max(1.0_dp, abs(derivative_fd)), &
    'ion local derivative mismatch' &
    )
  call test_end()

  call test_begin('cold drifting ion conserves flux and energy')
  call eval_cold_drift_density( &
    phi=-3.0_dp, charge=ev, mass=1836.0_dp*electron_mass, density_infinity=1.0e6_dp, &
    drift_infinity=2.0e4_dp, density=density, speed=speed, susceptibility=susceptibility, status=status &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'ion closure status mismatch')
  call assert_close_dp(density*speed, 2.0e10_dp, 1.0e-4_dp, 'ion flux is not conserved')
  call assert_close_dp( &
    speed*speed, (2.0e4_dp)**2 - 2.0_dp*ev*(-3.0_dp)/(1836.0_dp*electron_mass), &
    1.0e-5_dp, 'ion energy is not conserved' &
    )
  call test_end()

  call test_begin('cold drifting ion rejects an inaccessible potential')
  call eval_cold_drift_density( &
    phi=100.0_dp, charge=ev, mass=1836.0_dp*electron_mass, density_infinity=1.0e6_dp, &
    drift_infinity=1.0e3_dp, density=density, speed=speed, susceptibility=susceptibility, status=status &
    )
  call assert_equal_i32(status, outer_plasma_no_physical_solution, 'inaccessible ion must fail physically')
  call test_end()

  call test_begin('Bohm entry and photoelectron return use energy barriers')
  speed = kinetic_bohm_speed(2.0_dp*ev, 0.1_dp*ev, 1836.0_dp*electron_mass, 1.0_dp)
  call assert_close_dp( &
    speed*speed, 2.1_dp*ev/(1836.0_dp*electron_mass), 1.0e-5_dp, 'Bohm speed mismatch' &
    )
  call eval_photoelectron_escape_return( &
    phi_interface=3.0_dp, phi_infinity=0.0_dp, charge=-ev, temperature_j=1.5_dp*ev, &
    escape_fraction=escape_fraction, return_fraction=return_fraction, status=status &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'photoelectron orbit status mismatch')
  call assert_close_dp(escape_fraction, exp(-2.0_dp), 1.0e-14_dp, 'photoelectron escape mismatch')
  call assert_close_dp(return_fraction, 1.0_dp - exp(-2.0_dp), 1.0e-14_dp, 'photoelectron return mismatch')
  call test_end()

  call test_begin('current ledger keeps signed components explicit')
  call eval_kinetic_current_balance( &
    electron_absorption_flux=4.0_dp, ion_absorption_flux=1.0_dp, photoelectron_emission_flux=2.0_dp, &
    photoelectron_escape_fraction=0.25_dp, electron_charge=-ev, ion_charge=ev, &
    external_current_density=0.5_dp*ev, electron_current=electron_current, ion_current=ion_current, &
    photoelectron_current=photo_current, total_current=total_current &
    )
  call assert_close_dp(electron_current, -4.0_dp*ev, 1.0e-32_dp, 'electron current mismatch')
  call assert_close_dp(ion_current, ev, 1.0e-32_dp, 'ion current mismatch')
  call assert_close_dp(photo_current, 0.5_dp*ev, 1.0e-32_dp, 'photoelectron escape current mismatch')
  call assert_close_dp(total_current, -2.0_dp*ev, 1.0e-32_dp, 'total current mismatch')
  call test_end()

  call test_summary()
end program test_outer_plasma_kinetic_core

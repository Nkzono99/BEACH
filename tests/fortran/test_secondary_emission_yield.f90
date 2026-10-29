!> SEE の収率、確率丸め、共有放出 sampler の物理的契約。
program test_secondary_emission_yield
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: qe
  use bem_secondary_emission, only: secondary_emission_yield, secondary_emission_count
  use bem_photoelectron_injection, only: sample_photo_emission_velocity
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_close_dp, assert_true
  implicit none
  integer(i32) :: i, count
  real(dp) :: delta, mean_count, v(3), mean_energy, tangent_variance, sigma
  real(dp), parameter :: mass = 9.1093837139e-31_dp
  ! E/E_max = 0.5 and 2 at normal and 60 degree incidence, then exponents 5e-6 and 2e-5 on either side
  ! of the small-exponent series switch at 1e-5.
  real(dp), parameter :: whipple_energy(5) = [250.0_dp, 1000.0_dp, 500.0_dp, 0.03214728_dp, 0.08976656_dp]
  real(dp), parameter :: whipple_cosine(5) = [1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp, 1.0_dp]
  real(dp), parameter :: whipple_reference(5) = [ &
                         1.2590423225330811_dp, 2.4786675734494112_dp, 2.2731648290899127_dp, &
                         2.4495394586735255e-4_dp, 6.839928032742621e-4_dp]

  call test_init(3)
  call test_begin('whipple_and_constant_yield')
  delta = secondary_emission_yield('whipple', 500.0_dp, 1.0_dp, 1.5_dp, 500.0_dp, 0.0_dp)
  call assert_close_dp(delta, 1.114_dp*1.5_dp*(1.0_dp - exp(-2.28_dp)), 1.0e-14_dp, 'normal incidence at E_max')
  call assert_close_dp(secondary_emission_yield('whipple', 0.0_dp, 1.0_dp, 1.5_dp, 500.0_dp, 0.0_dp), &
                       0.0_dp, 0.0_dp, 'zero-energy yield')
  call assert_close_dp(secondary_emission_yield('whipple', -1.0_dp, 1.0_dp, 1.5_dp, 500.0_dp, 0.0_dp), &
                       0.0_dp, 0.0_dp, 'negative-energy yield')
  delta = 1.114_dp*1.5_dp/0.1_dp*(1.0_dp - exp(-0.228_dp))
  call assert_close_dp(secondary_emission_yield('whipple', 500.0_dp, 0.0_dp, 1.5_dp, 500.0_dp, 0.0_dp), &
                       delta, 1.0e-14_dp, 'grazing-angle cosine must clamp to 0.1')
  call assert_close_dp(secondary_emission_yield('constant', 123.0_dp, 0.2_dp, 0.0_dp, 0.0_dp, 2.0_dp), &
                       2.0_dp, 0.0_dp, 'constant yield')
  ! Away from E = E_max the curve depends on (E/E_max)**1.35 and (E/E_max)**(-0.35). References are
  ! 1.114 y_max / c r**(-0.35) (-expm1(-2.28 c r**1.35)) evaluated independently in double precision.
  do i = 1, size(whipple_energy)
    delta = secondary_emission_yield('whipple', whipple_energy(i), whipple_cosine(i), 1.5_dp, 500.0_dp, 0.0_dp)
    ! 1e-10 leaves room for the eps/exponent cancellation of 1 - exp(-x) just above the switch.
    call assert_close_dp(delta, whipple_reference(i), 1.0e-10_dp*whipple_reference(i), &
                         'Whipple yield off E_max mismatch')
  end do
  call test_end()

  call test_begin('bernoulli_rounding_preserves_charge_scaled_expectation')
  mean_count = 0.0_dp
  do i = 1, 20000
    count = secondary_emission_count(0.65_dp, 2.0_dp, 73_i32, int([1, 0, i, 2, 0], i64))
    call assert_true(count == 1 .or. count == 2, 'rounding must stay between floor and ceiling')
    mean_count = mean_count + real(count, dp)
  end do
  call assert_close_dp(mean_count/20000.0_dp, 1.3_dp, 0.015_dp, 'Bernoulli charge-scaled expectation')
  call test_end()

  call test_begin('shared_half_maxwell_sampler_and_emission_charge_conservation')
  sigma = sqrt(2.0_dp*qe/mass)
  mean_energy = 0.0_dp
  tangent_variance = 0.0_dp
  do i = 1, 10000
    call sample_photo_emission_velocity(sigma, 0.0_dp, [0.0_dp, 0.0_dp, 1.0_dp], &
                                        [1.0_dp, 0.0_dp, 0.0_dp], [0.0_dp, 1.0_dp, 0.0_dp], v, &
                                        rng_seed=73_i32, rng_counter=int([2, 0, i, 1, 1], i64))
    call assert_true(v(3) >= 0.0_dp, 'normal velocity must be outward')
    mean_energy = mean_energy + 0.5_dp*mass*v(3)**2/qe
    tangent_variance = tangent_variance + sum(v(1:2)**2)/(2.0_dp*sigma**2)
  end do
  call assert_close_dp(mean_energy/10000.0_dp, 2.0_dp, 0.06_dp, 'flux-weighted normal kinetic energy mean equals T')
  call assert_close_dp(tangent_variance/10000.0_dp, 1.0_dp, 0.04_dp, 'Gaussian tangential variance')
  call test_end()
  call test_summary()
end program test_secondary_emission_yield

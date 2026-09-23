!> PE エネルギー分布の保存則、Maxwell 極限、独立積分との整合を検証する。
program test_pe_spectrum
  use bem_kinds, only: dp, i32
  use bem_constants, only: qe, pi
  use bem_pe_spectrum, only: pe_spectrum_type
  use bem_sheath_model_orbits, only: electron_density
  use bem_sheath_model_core, only: zhao_params_type, evaluate_zhao_rho_hat, integrate_zhao_rho, swe_free_current_term
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp
  implicit none

  real(dp), parameter :: electron_mass = 9.1093837015e-31_dp
  real(dp), parameter :: source_flux = 3.1e13_dp, temperature = 2.4_dp
  type(pe_spectrum_type) :: spectrum, other, saved, fine, low, high
  type(zhao_params_type) :: params
  real(dp) :: free, returning, expected_free, expected_returning, value, expected, barrier
  real(dp) :: first_mean, low_mean, high_mean, high_weight, previous_flux, previous_tail
  real(dp) :: phi_h, phi_min, phi, gap, reference_density, fine_error, coarse_error
  real(dp) :: drift, psi, orbit_barrier, local_flux, current_term, rho, plus, minus
  real(dp), parameter :: orbit_psi(4) = [0.25_dp, -0.03_dp, 0.25_dp, -0.3_dp]
  real(dp), parameter :: orbit_barriers(4) = [-0.08_dp, -0.08_dp, 0.0_dp, -0.5_dp]
  integer :: j, n, k
  logical :: upper

  call test_init(12)

  call test_begin('empty_and_reset_spectra_have_zero_flux_and_density')
  call assert_close_dp(spectrum%total_flux(), 0.0_dp, 0.0_dp, 'empty total flux')
  call assert_close_dp(spectrum%mean_energy(), 0.0_dp, 0.0_dp, 'empty mean-energy convention')
  call assert_close_dp(spectrum%tail_flux(-1.0_dp), 0.0_dp, 0.0_dp, 'empty tail')
  call spectrum%density(1.0_dp, 2.0_dp, -1.0_dp, electron_mass, free, returning)
  call assert_close_dp(free + returning, 0.0_dp, 0.0_dp, 'empty density')
  call assert_close_dp(spectrum%integrated_density(1.0_dp, 2.0_dp, -1.0_dp, electron_mass), &
                       0.0_dp, 0.0_dp, 'empty density primitive')
  call spectrum%add(2.0_dp, 7.0_dp)
  n = size(spectrum%flux)
  call spectrum%reset()
  call assert_true(size(spectrum%flux) == n, 'reset retains grid extent')
  call assert_close_dp(spectrum%total_flux(), 0.0_dp, 0.0_dp, 'reset clears weights')
  call assert_close_dp(spectrum%mean_energy(), 0.0_dp, 0.0_dp, 'zero-flux mean-energy convention')
  call spectrum%clear()
  call assert_true(.not. allocated(spectrum%flux), 'clear deallocates weights')
  call spectrum%bootstrap_maxwellian(0.0_dp, temperature)
  call assert_close_dp(spectrum%total_flux(), 0.0_dp, 0.0_dp, 'zero-flux Maxwell bootstrap')
  call test_end()

  call test_begin('sampling_preserves_all_weights_and_extends_high_energy_support')
  call spectrum%add(0.0_dp, 1.0_dp)
  call spectrum%add(1.0e-25_dp, 2.0_dp)
  value = spectrum%edge(7)
  call spectrum%add(value, 4.0_dp)
  call assert_close_dp(spectrum%flux(8), 4.0_dp, 0.0_dp, 'exact-edge sample enters its upper bin')
  call spectrum%add(1.0e12_dp, 8.0_dp)
  call assert_close_dp(spectrum%total_flux(), 15.0_dp, 0.0_dp, 'sampled flux conserved')
  call assert_close_dp(spectrum%tail_flux(1.0e6_dp), 8.0_dp, 0.0_dp, 'high-energy samples retained')
  call assert_true(spectrum%edge(size(spectrum%flux)) > 1.0e12_dp, 'energy support extends to the sample')
  call assert_close_dp(spectrum%tail_flux(-1.0_dp), 15.0_dp, 0.0_dp, 'negative barrier transmits all flux')
  call assert_close_dp(spectrum%tail_flux(spectrum%edge(size(spectrum%flux))), &
                       0.0_dp, 0.0_dp, 'barrier beyond support transmits no flux')
  ! The energy/scale ratio overflows, but its logarithm and the grid edge do not.
  call other%clear()
  other%energy_scale_ev = 1.0e-300_dp
  call other%add(1.0e300_dp, 3.0_dp)
  call assert_close_dp(other%total_flux(), 3.0_dp, 0.0_dp, 'extreme finite energy is not discarded')
  call assert_true(other%edge(size(other%flux)) > 1.0e300_dp, 'scaled edge avoids intermediate overflow')
  call other%clear()
  other%energy_scale_ev = spectrum%energy_scale_ev
  call test_end()

  call test_begin('combination_conserves_flux_tail_and_first_energy_moment')
  call spectrum%clear()
  call spectrum%add(0.7_dp, 11.0_dp)
  call spectrum%add(2.0_dp, 7.0_dp)
  call other%add(10.0_dp, 13.0_dp)
  previous_flux = spectrum%total_flux()
  previous_tail = spectrum%tail_flux(1.0_dp)
  first_mean = spectrum%mean_energy()
  expected = (0.3_dp*previous_flux*first_mean + 0.7_dp*other%total_flux()*other%mean_energy())/ &
             (0.3_dp*previous_flux + 0.7_dp*other%total_flux())
  call spectrum%combine(other, 0.3_dp, 0.7_dp)
  call assert_close_dp(spectrum%total_flux(), 0.3_dp*previous_flux + 0.7_dp*other%total_flux(), &
                       1.0e-13_dp, 'combined flux')
  call assert_close_dp(spectrum%tail_flux(1.0_dp), 0.3_dp*previous_tail + 0.7_dp*other%tail_flux(1.0_dp), &
                       1.0e-13_dp, 'combined tail')
  call assert_close_dp(spectrum%mean_energy(), expected, 1.0e-13_dp, 'combined energy moment')
  call assert_true(size(spectrum%flux) == size(other%flux), 'combination zero-pads shorter support')
  saved = spectrum
  call other%clear()
  call spectrum%combine(other, 1.0_dp, 1.0_dp)
  call assert_close_dp(spectrum%total_flux(), saved%total_flux(), 0.0_dp, 'empty addition is neutral')
  call test_end()

  call test_begin('maxwellian_bootstrap_converges_to_analytic_tail_and_mean')
  spectrum%energy_scale_ev = 0.25_dp*temperature
  call spectrum%bootstrap_maxwellian(source_flux, temperature)
  fine%energy_scale_ev = spectrum%energy_scale_ev
  fine%bins_per_decade = 4096_i32
  call fine%bootstrap_maxwellian(source_flux, temperature)
  call assert_true(fine%edge(size(fine%flux)) >= 40.0_dp*temperature, 'bootstrap resolves the thermal tail')
  call assert_close_dp(fine%total_flux(), source_flux, 1.0e-13_dp*source_flux, 'bootstrap flux normalization')
  call assert_close_dp(fine%mean_energy(), temperature, 3.0e-7_dp*temperature, 'Maxwell mean normal energy')
  coarse_error = 0.0_dp
  fine_error = 0.0_dp
  do j = 0, 12
    barrier = 0.73_dp*real(j, dp)
    expected = source_flux*exp(-barrier/temperature)
    coarse_error = max(coarse_error, abs(spectrum%tail_flux(barrier) - expected)/expected)
    fine_error = max(fine_error, abs(fine%tail_flux(barrier) - expected)/expected)
  end do
  call assert_true(fine_error < 1.0e-6_dp, 'Maxwell tail agrees with its independent exponential formula')
  call assert_true(fine_error < 1.0e-3_dp*coarse_error, 'grid refinement improves Maxwell tail accuracy')
  call test_end()

  call test_begin('density_matches_maxwell_velocity_integrals_on_both_sides')
  phi_h = 1.7_dp
  phi_min = -0.8_dp
  do j = 0, 6
    phi = phi_min + 0.5_dp*real(j, dp)
    call maxwell_density(phi, phi_h, phi_min, expected_free, expected_returning)
    call fine%density(phi, phi_h, phi_min, electron_mass, free, returning)
    call assert_close_dp(free, expected_free, 2.0e-5_dp*expected_free, 'free Maxwell density')
    call assert_close_dp(returning, expected_returning, &
                         2.0e-5_dp*max(expected_free, expected_returning), 'both reflected legs are included')
    call fine%density(phi, phi_h, phi_min, electron_mass, free, returning, upper_side=.true.)
    call assert_close_dp(free, expected_free, 2.0e-5_dp*expected_free, 'upper-side free Maxwell density')
    call assert_close_dp(returning, 0.0_dp, 0.0_dp, 'upper side has no reflected PE orbits')
  end do
  call maxwell_density(-1.0_dp, -3.0_dp, -3.0_dp, expected_free, expected_returning)
  call fine%density(-1.0_dp, -3.0_dp, -3.0_dp, electron_mass, free, returning, upper_side=.true.)
  call assert_close_dp(free, expected_free, 2.0e-5_dp*expected_free, 'monotonic C accelerates all outward PE')
  call assert_close_dp(returning, 0.0_dp, 0.0_dp, 'monotonic C has no reflected PE')
  call test_end()

  call test_begin('exact_density_primitive_matches_independent_maxwell_quadrature')
  do j = 0, 5
    upper = mod(j, 2) == 1
    phi_h = 1.7_dp
    phi_min = -0.8_dp
    phi = phi_min + 0.7_dp*real(1 + j/2, dp)
    value = fine%integrated_density(phi, phi_h, phi_min, electron_mass, upper_side=upper)
    expected = maxwell_integral_oracle(phi, phi_h, phi_min, upper)
    call assert_close_dp(value, expected, 2.0e-5_dp*expected, 'Maxwell density primitive versus quadrature')
  end do
  value = fine%integrated_density(-1.0_dp, -3.0_dp, -3.0_dp, electron_mass, upper_side=.true.)
  expected = maxwell_integral_oracle(-1.0_dp, -3.0_dp, -3.0_dp, .true.)
  call assert_close_dp(value, expected, 2.0e-5_dp*expected, 'monotonic C density primitive')
  ! A signed primitive remains meaningful when evaluating below the chosen reference.
  value = fine%integrated_density(-1.3_dp, 1.7_dp, -0.8_dp, electron_mass, upper_side=.true.)
  expected = maxwell_integral_oracle(-1.3_dp, 1.7_dp, -0.8_dp, .true.)
  call assert_close_dp(value, expected, 2.0e-5_dp*abs(expected), 'signed primitive below its reference')
  call test_end()

  call test_begin('density_primitive_is_stable_arbitrarily_close_to_the_minimum')
  phi_h = 1.7_dp
  phi_min = -0.8_dp
  call fine%density(phi_min, phi_h, phi_min, electron_mass, free, returning)
  reference_density = free
  call assert_close_dp(fine%integrated_density(phi_min, phi_h, phi_min, electron_mass), &
                       0.0_dp, 0.0_dp, 'primitive is exactly zero at its reference')
  do j = 4, 12, 2
    phi = phi_min + 10.0_dp**(-j)
    gap = phi - phi_min
    value = fine%integrated_density(phi, phi_h, phi_min, electron_mass, upper_side=.true.)
    ! For this branch n changes by O(sqrt(gap)); avoid asserting a flat profile.
    call assert_close_dp(value/gap, reference_density, &
                         (2.0_dp*sqrt(gap) + 2.0e-4_dp)*reference_density, 'small-gap primitive retains significance')
  end do
  call test_end()

  call test_begin('equal_flux_and_mean_energy_do_not_imply_equal_escape_flux')
  call spectrum%clear()
  spectrum%energy_scale_ev = 1.0_dp
  call spectrum%add(2.0_dp, 1.0_dp)
  call low%add(0.5_dp, 1.0_dp)
  call high%add(3.5_dp, 1.0_dp)
  first_mean = spectrum%mean_energy()
  low_mean = low%mean_energy()
  high_mean = high%mean_energy()
  high_weight = (first_mean - low_mean)/(high_mean - low_mean)
  call low%combine(high, 1.0_dp - high_weight, high_weight)
  call assert_close_dp(low%total_flux(), spectrum%total_flux(), 1.0e-14_dp, 'equal source flux')
  call assert_close_dp(low%mean_energy(), first_mean, 1.0e-14_dp, 'equal source mean energy')
  call assert_close_dp(spectrum%tail_flux(3.0_dp), 0.0_dp, 0.0_dp, 'narrow 2 eV source cannot cross 3 eV barrier')
  call assert_true(low%tail_flux(3.0_dp) > 0.45_dp .and. low%tail_flux(3.0_dp) < 0.55_dp, &
                   'bimodal source transmits approximately half the flux')
  call test_end()

  call test_begin('piecewise_density_primitive_differentiates_to_the_actual_spectrum')
  call low%combine(spectrum, 0.4_dp, 0.6_dp)
  phi_h = 1.7_dp
  phi_min = -0.8_dp
  do j = 0, 12
    phi = phi_min + 0.17_dp*real(j + 1, dp)
    gap = 1.0e-6_dp
    call low%density(phi, phi_h, phi_min, electron_mass, free, returning)
    value = (low%integrated_density(phi + gap, phi_h, phi_min, electron_mass) - &
             low%integrated_density(phi - gap, phi_h, phi_min, electron_mass))/(2.0_dp*gap)
    call assert_close_dp(value, free + returning, 2.0e-7_dp*max(free + returning, 1.0e-20_dp), &
                         'non-Maxwell primitive derivative equals its orbit density')
  end do
  call test_end()

  call test_begin('drifting_electron_density_and_flux_preserve_upstream_orbits')
  do k = 0, 2
    drift = 0.2_dp*real(k - 1, dp)
    do j = 1, size(orbit_psi)
      psi = orbit_psi(j)
      orbit_barrier = orbit_barriers(j)
      call electron_density(psi, orbit_barrier, drift, free, returning)
      expected = upstream_density_oracle(psi, orbit_barrier, drift)
      call assert_close_dp(free, expected, 3.0e-9_dp, 'passing density versus upstream-velocity integral')
      expected = upstream_reflected_density_oracle(psi, orbit_barrier, drift)
      call assert_close_dp(returning, expected, 3.0e-9_dp, 'both reflected legs versus upstream-velocity integral')
      local_flux = local_flux_oracle(psi, orbit_barrier, drift)
      params%t_swe_ev = 12.0_dp
      params%t_phe_ev = temperature
      params%v_swe_th_mps = sqrt(2.0_dp*qe*params%t_swe_ev/electron_mass)
      params%v_phe_th_mps = sqrt(2.0_dp*qe*params%t_phe_ev/electron_mass)
      params%v_d_electron_mps = drift*params%v_swe_th_mps
      current_term = params%v_phe_th_mps/(2.0_dp*sqrt(pi))* &
                     swe_free_current_term(params, 1.0_dp, sqrt(-orbit_barrier) - drift)
      call assert_close_dp(local_flux, current_term/params%v_swe_th_mps, 3.0e-10_dp, &
                           'local VDF flux equals the upstream current expression')
    end do
  end do
  call test_end()

  call test_begin('sagdeev_integral_uses_the_same_orbit_density_and_spectrum')
  params%n_swi_inf_m3 = 8.7e6_dp
  params%n_phe_ref_m3 = 8.7e6_dp
  params%n_phe0_m3 = 2.0_dp*sqrt(pi)*source_flux/params%v_phe_th_mps
  params%t_swe_ev = 12.0_dp
  params%t_phe_ev = temperature
  params%tau = params%t_swe_ev/params%t_phe_ev
  params%mach = 5.0_dp
  params%m_e_kg = electron_mass
  params%pe_spectrum = fine
  do k = 0, 1
    params%u = 0.2_dp*real(k, dp)
    do j = 0, 1
      phi_h = 1.7_dp/temperature
      phi_min = -0.8_dp/temperature
      phi = -0.3_dp/temperature
      gap = 1.0e-5_dp
      if (j == 0) then
        call evaluate_zhao_rho_hat(params, 'A', 'lower', phi, phi_h, phi_min, 0.9_dp, rho)
        plus = integrate_zhao_rho(params, 'A', 'lower', phi_min, phi + gap, phi_h, phi_min, 0.9_dp)
        minus = integrate_zhao_rho(params, 'A', 'lower', phi_min, phi - gap, phi_h, phi_min, 0.9_dp)
      else
        call evaluate_zhao_rho_hat(params, 'A', 'upper', phi, phi_h, phi_min, 0.9_dp, rho)
        plus = integrate_zhao_rho(params, 'A', 'upper', phi_min, phi + gap, phi_h, phi_min, 0.9_dp)
        minus = integrate_zhao_rho(params, 'A', 'upper', phi_min, phi - gap, phi_h, phi_min, 0.9_dp)
      end if
      call assert_close_dp((plus - minus)/(2.0_dp*gap), rho, 3.0e-6_dp*max(1.0_dp, abs(rho)), &
                           'Sagdeev derivative matches the actual total space-charge density')
    end do
  end do
  call test_end()

  call test_begin('positive_drift_algebraic_root_has_negative_upstream_field_squared')
  call params%pe_spectrum%clear()
  params%n_phe_ref_m3 = 64.0e6_dp
  params%n_phe0_m3 = 64.0e6_dp*sin(pi/3.0_dp)
  params%t_phe_ev = 2.2_dp
  params%t_swe_ev = 12.0_dp
  params%tau = params%t_swe_ev/params%t_phe_ev
  params%v_d_ion_mps = 468.0e3_dp*sin(pi/3.0_dp)
  params%mach = params%v_d_ion_mps/sqrt(qe*params%t_swe_ev/1.67262192369e-27_dp)
  params%u = params%v_d_ion_mps/sqrt(2.0_dp*qe*params%t_swe_ev/electron_mass)
  value = 2.0_dp*integrate_zhao_rho(params, 'A', 'upper', -0.01_dp, 0.0_dp, &
                                    2.6544268139403324_dp/params%t_phe_ev, &
                                    -1.1897238115913789_dp/params%t_phe_ev, 7923268.55824609_dp/params%n_phe_ref_m3)
  call assert_true(value < -1.0e-7_dp, 'independent upstream field obstruction survives without the profile guard')
  call test_end()

  call test_summary()

contains

  ! An independent Simpson integral in upstream speed with its density Jacobian.
  real(dp) function upstream_density_oracle(psi, barrier, drift) result(value)
    real(dp), intent(in) :: psi, barrier, drift
    integer, parameter :: panels = 32768
    real(dp) :: a, h, weight
    integer :: k

    h = 12.0_dp/real(panels, dp)
    value = 0.0_dp
    do k = 0, panels
      a = sqrt(-barrier) + h*real(k, dp)
      weight = 2.0_dp
      if (mod(k, 2) == 1) weight = 4.0_dp
      if (k == 0 .or. k == panels) weight = 1.0_dp
      value = value + weight*a*exp(-(a - drift)**2)/sqrt(a*a + psi)
    end do
    value = value*h/(3.0_dp*sqrt(pi))
  end function upstream_density_oracle

  real(dp) function upstream_reflected_density_oracle(psi, barrier, drift) result(value)
    real(dp), intent(in) :: psi, barrier, drift
    integer, parameter :: panels = 8192
    real(dp) :: t, a, minimum_squared, h, weight, jacobian
    integer :: k

    minimum_squared = max(0.0_dp, -psi)
    h = sqrt(max(0.0_dp, -barrier - minimum_squared))/real(panels, dp)
    value = 0.0_dp
    if (h == 0.0_dp) return
    ! a^2=a_min^2+t^2 removes the upstream density Jacobian's turning-point pole.
    do k = 0, panels
      t = h*real(k, dp)
      a = sqrt(minimum_squared + t*t)
      jacobian = 1.0_dp
      if (psi > 0.0_dp) jacobian = t/sqrt(t*t + psi)
      weight = 2.0_dp
      if (mod(k, 2) == 1) weight = 4.0_dp
      if (k == 0 .or. k == panels) weight = 1.0_dp
      value = value + weight*jacobian*exp(-(a - drift)**2)
    end do
    value = 2.0_dp*value*h/(3.0_dp*sqrt(pi))
  end function upstream_reflected_density_oracle

  ! Integrate w*f(w) at the local potential; the squared map resolves the cutoff.
  real(dp) function local_flux_oracle(psi, barrier, drift) result(value)
    real(dp), intent(in) :: psi, barrier, drift
    integer, parameter :: panels = 8192
    real(dp) :: t, w, upstream, h, weight
    integer :: k

    h = 4.0_dp/real(panels, dp)
    value = 0.0_dp
    do k = 0, panels
      t = h*real(k, dp)
      w = sqrt(psi - barrier) + t*t
      upstream = sqrt(max(0.0_dp, w*w - psi))
      weight = 2.0_dp
      if (mod(k, 2) == 1) weight = 4.0_dp
      if (k == 0 .or. k == panels) weight = 1.0_dp
      value = value + weight*w*exp(-(upstream - drift)**2)*2.0_dp*t
    end do
    value = value*h/(3.0_dp*sqrt(pi))
  end function local_flux_oracle

  ! Analytic velocity integrals of F(K)=Gamma/T*exp(-K/T), independent of bins.
  subroutine maxwell_density(phi, phi_h, phi_min, free, returning)
    real(dp), intent(in) :: phi, phi_h, phi_min
    real(dp), intent(out) :: free, returning
    real(dp) :: b, escape_energy, prefactor, start

    b = phi_h - phi
    escape_energy = max(phi_h - phi_min, 0.0_dp)
    prefactor = source_flux*sqrt(pi*electron_mass/(2.0_dp*qe*temperature))*exp(-b/temperature)
    free = prefactor*erfc(sqrt(max(escape_energy - b, 0.0_dp)/temperature))
    start = max(0.0_dp, b)
    returning = 0.0_dp
    if (escape_energy > start) returning = 2.0_dp*prefactor*( &
                                           erf(sqrt((escape_energy - b)/temperature)) - erf(sqrt((start - b)/temperature)))
  end subroutine maxwell_density

  ! Simpson quadrature in phi=phi_min+(phi-phi_min)*t^2 resolves the endpoint cusp.
  real(dp) function maxwell_integral_oracle(phi, phi_h, phi_min, upper) result(value)
    real(dp), intent(in) :: phi, phi_h, phi_min
    logical, intent(in) :: upper
    integer, parameter :: panels = 32768
    real(dp) :: t, local_phi, free, returning, weight
    integer :: k

    value = 0.0_dp
    do k = 0, panels
      t = real(k, dp)/real(panels, dp)
      local_phi = phi_min + (phi - phi_min)*t*t
      call maxwell_density(local_phi, phi_h, phi_min, free, returning)
      if (upper) returning = 0.0_dp
      weight = 2.0_dp
      if (mod(k, 2) == 1) weight = 4.0_dp
      if (k == 0 .or. k == panels) weight = 1.0_dp
      value = value + weight*(free + returning)*2.0_dp*(phi - phi_min)*t
    end do
    value = value/(3.0_dp*real(panels, dp))
  end function maxwell_integral_oracle

end program test_pe_spectrum

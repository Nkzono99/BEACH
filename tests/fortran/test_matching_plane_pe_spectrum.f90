!> Spectrum closure の根・中性条件・接続場を独立な Maxwell 混合分布で検証する。
program test_matching_plane_pe_spectrum
  use bem_kinds, only: dp, i32
  use bem_constants, only: qe, eps0, pi
  use bem_pe_spectrum, only: pe_spectrum_type
  use bem_matching_plane_zhao, only: matching_plane_zhao_model_type, matching_plane_zhao_diagnostics_type, &
                                     matching_plane_zhao_root_seed_type, matching_plane_zhao_ok
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, assert_equal_i32
  implicit none

  real(dp), parameter :: me = 9.1093837015e-31_dp, mi = 1.67262192369e-27_dp
  real(dp), parameter :: ni = 5.0e6_dp, te = 10.0_dp, mach = 10.0_dp
  real(dp), parameter :: flux_unit = ni*sqrt(qe*te/me), d_unit = sqrt(eps0*ni*qe*te)
  type(matching_plane_zhao_model_type) :: spectral, moments
  type(matching_plane_zhao_diagnostics_type) :: diagnostic, reference_diagnostic
  type(matching_plane_zhao_root_seed_type) :: seed
  type(pe_spectrum_type) :: spectrum, second, empty
  real(dp) :: input(5), response(6), reference(6), drift, g, field, fraction
  real(dp) :: free, returning, expected, h, minimum, amplitude, escape, g_components(2), t_components(2)
  real(dp) :: finite_ni, finite_te, finite_g, finite_mach, finite_u, finite_d_unit, finite_flux_unit
  integer(i32) :: status
  integer :: branch_index, drift_index
  character(len=1) :: branch
  character(len=512) :: message

  call test_init(7)
  spectrum%energy_scale_ev = 0.1_dp*te
  spectrum%bins_per_decade = 1024_i32
  second%energy_scale_ev = spectrum%energy_scale_ev
  second%bins_per_decade = spectrum%bins_per_decade

  call test_begin('maxwell_spectrum_reproduces_B_and_C_moment_roots_with_and_without_drift')
  do branch_index = 1, 2
    branch = 'B'
    g = 0.4_dp
    field = 0.5_dp
    if (branch_index == 2) then
      branch = 'C'
      g = 0.1_dp
      field = -0.5_dp
    end if
    do drift_index = 0, 1
      drift = 0.2_dp*real(drift_index, dp)
      call spectrum%bootstrap_maxwellian(g*flux_unit, 0.2_dp*te)
      call initialize(spectral, branch, drift)
      call initialize(moments, branch, drift)
      call spectral%set_photoelectron_spectrum(spectrum)
      input = [field*d_unit, spectrum%total_flux(), spectrum%mean_energy(), 0.0_dp, 0.0_dp]
      call spectral%evaluate(input, response, status, message, diagnostic)
      call assert_equal_i32(status, matching_plane_zhao_ok, 'spectral '//branch//' response: '//trim(message))
      if (status /= matching_plane_zhao_ok) cycle
      call moments%evaluate(input, reference, status, message, reference_diagnostic)
      call assert_equal_i32(status, matching_plane_zhao_ok, 'moment '//branch//' response: '//trim(message))
      if (status /= matching_plane_zhao_ok) cycle
      call assert_close_dp(response(1)/te, reference(1)/te, 2.0e-4_dp, 'Maxwell interface potential')
      call assert_close_dp(diagnostic%ambient_electron_density_m3/ni, &
                           reference_diagnostic%ambient_electron_density_m3/ni, 2.0e-4_dp, 'Maxwell ambient amplitude')
      call assert_close_dp(response(2)/flux_unit, reference(2)/flux_unit, 2.0e-4_dp, 'Maxwell electron inward flux')
      minimum = 0.0_dp
      if (branch == 'C') minimum = response(1)/te
      call check_orbit_closure(branch, response(1)/te, minimum, diagnostic%ambient_electron_density_m3/ni, &
                               field, mach, drift, [g, 0.0_dp], [0.2_dp, 1.0_dp])
    end do
  end do
  call test_end()

  call test_begin('zero_drift_A_matches_the_independent_source_kinetic_oracle')
  ! Archived independent solver, M=10, Tph/Te=.2, G=.4, EH=.5.
  ! The old closed-form A integral is singular at u=0, so it is not the oracle.
  call spectrum%bootstrap_maxwellian(0.4_dp*flux_unit, 0.2_dp*te)
  call initialize(spectral, 'A', 0.0_dp, 'continuation')
  call spectral%set_photoelectron_spectrum(spectrum)
  seed%valid = .true.
  seed%phi0_v = 0.1854122640629287_dp*te
  seed%phi_m_v = -0.01182341861232325_dp*te
  seed%ambient_electron_density_m3 = 1.204317385077186_dp*ni
  input = [0.5_dp*d_unit, spectrum%total_flux(), spectrum%mean_energy(), 0.0_dp, 0.0_dp]
  call spectral%evaluate(input, response, status, message, diagnostic, continuation_seed=seed)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'zero-drift spectral A: '//trim(message))
  if (status == matching_plane_zhao_ok) then
    call assert_close_dp(response(1)/te, 0.1854122640629287_dp, 2.0e-4_dp, 'independent A interface potential')
    call assert_close_dp(response(6)/te, -0.01182341861232325_dp, 2.0e-4_dp, 'independent A potential minimum')
    call assert_close_dp(diagnostic%ambient_electron_density_m3/ni, 1.204317385077186_dp, &
                         2.0e-4_dp, 'independent A ambient amplitude')
    escape = spectrum%tail_flux(response(1) - response(6))/flux_unit
    call assert_close_dp(escape, 0.1491997681377198_dp, 2.0e-4_dp, 'independent A escaping flux')
    call check_orbit_closure('A', response(1)/te, response(6)/te, diagnostic%ambient_electron_density_m3/ni, &
                             0.5_dp, mach, 0.0_dp, [0.4_dp, 0.0_dp], [0.2_dp, 1.0_dp])
  end if
  call test_end()

  call test_begin('finite_drift_A_preserves_the_maxwell_limit_and_actual_upper_integral')
  finite_ni = 8.7e6_dp
  finite_te = 12.0_dp
  finite_flux_unit = finite_ni*sqrt(qe*finite_te/me)
  finite_d_unit = sqrt(eps0*finite_ni*qe*finite_te)
  drift = 4.0529988897111727e5_dp
  finite_mach = drift/sqrt(qe*finite_te/mi)
  finite_u = drift/sqrt(2.0_dp*qe*finite_te/me)
  finite_g = 1.3754433596232731e13_dp/finite_flux_unit
  call spectrum%bootstrap_maxwellian(1.3754433596232731e13_dp, 2.2_dp)
  call moments%initialize('a', 'minimum_energy', finite_ni, finite_te, drift, drift, mi, me, 2.2_dp, status, message)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'finite-drift moment initialization')
  call moments%set_photoelectron_spectrum(empty)
  input = [1.4187346568707933e-11_dp, spectrum%total_flux(), 2.2_dp, 0.0_dp, 0.0_dp]
  call moments%evaluate(input, reference, status, message, reference_diagnostic)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'finite-drift moment A: '//trim(message))
  if (status == matching_plane_zhao_ok) then
    call spectral%initialize('a', 'continuation', finite_ni, finite_te, drift, drift, mi, me, 2.2_dp, status, message)
    call assert_equal_i32(status, matching_plane_zhao_ok, 'finite-drift spectral initialization')
    call spectral%set_photoelectron_spectrum(spectrum)
    seed%phi0_v = reference(1)
    seed%phi_m_v = reference(6)
    seed%ambient_electron_density_m3 = reference_diagnostic%ambient_electron_density_m3
    call spectral%evaluate(input, response, status, message, diagnostic, continuation_seed=seed)
    call assert_equal_i32(status, matching_plane_zhao_ok, 'finite-drift spectral A: '//trim(message))
    if (status == matching_plane_zhao_ok) then
      call assert_close_dp(response(1)/finite_te, reference(1)/finite_te, 2.0e-4_dp, 'finite-drift Maxwell A potential')
      call assert_close_dp(response(6)/finite_te, reference(6)/finite_te, 2.0e-4_dp, 'finite-drift Maxwell A minimum')
      call assert_close_dp(diagnostic%ambient_electron_density_m3/finite_ni, &
                           reference_diagnostic%ambient_electron_density_m3/finite_ni, &
                           2.0e-4_dp, 'finite-drift Maxwell A ambient amplitude')
      call check_orbit_closure('A', response(1)/finite_te, response(6)/finite_te, &
                               diagnostic%ambient_electron_density_m3/finite_ni, input(1)/finite_d_unit, finite_mach, &
                               finite_u, [finite_g, 0.0_dp], [2.2_dp/finite_te, 1.0_dp])
    end if
  end if
  call test_end()

  call test_begin('same_H_flux_and_mean_have_different_roots_and_escape_for_a_filtered_mixture')
  ! Exact continuum transport of 0.6*Maxwell(.04)+0.4*Maxwell(.4) through
  ! a retarding inner drop .04 Te; normalize the transmitted source to G=.4.
  fraction = 0.6_dp*exp(-1.0_dp)/(0.6_dp*exp(-1.0_dp) + 0.4_dp*exp(-0.1_dp))
  g_components = 0.4_dp*[fraction, 1.0_dp - fraction]
  t_components = [0.04_dp, 0.4_dp]
  call spectrum%bootstrap_maxwellian(g_components(1)*flux_unit, t_components(1)*te)
  call second%bootstrap_maxwellian(g_components(2)*flux_unit, t_components(2)*te)
  call spectrum%combine(second, 1.0_dp, 1.0_dp)
  call initialize(spectral, 'B', 0.0_dp)
  call initialize(moments, 'B', 0.0_dp)
  call spectral%set_photoelectron_spectrum(spectrum)
  input = [0.3_dp*d_unit, spectrum%total_flux(), spectrum%mean_energy(), 0.0_dp, 0.0_dp]
  call spectral%evaluate(input, response, status, message, diagnostic)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'mixture spectral B: '//trim(message))
  if (status == matching_plane_zhao_ok) then
    h = response(1)/te
    amplitude = diagnostic%ambient_electron_density_m3/ni
    call assert_close_dp(h, 0.0582_dp, 2.0e-4_dp, 'independent transported-mixture root')
    call check_orbit_closure('B', h, 0.0_dp, amplitude, 0.3_dp, mach, 0.0_dp, g_components, t_components)
    escape = spectrum%tail_flux(response(1))/flux_unit
    expected = sum(g_components*exp(-h/t_components))
    call assert_close_dp(escape, expected, 2.0e-5_dp, 'root barrier uses the same mixture escape integral')
    call moments%evaluate(input, reference, status, message, reference_diagnostic)
    call assert_equal_i32(status, matching_plane_zhao_ok, 'moment-matched B: '//trim(message))
    if (status == matching_plane_zhao_ok) then
      call assert_true(reference(1)/te - h > 0.03_dp, 'preserving shape changes the root despite identical moments')
      call assert_true(abs(input(2)*exp(-reference(1)/input(3))/flux_unit - escape) > 0.015_dp, &
                       'preserving shape changes net escape at the respective roots')
    end if
  end if
  call test_end()

  call test_begin('zero_field_nonmaxwell_neutrality_uses_density_instead_of_mean_energy')
  spectrum%flux = 0.25_dp*spectrum%flux
  call initialize(spectral, 'B', 0.2_dp)
  call spectral%set_photoelectron_spectrum(spectrum)
  input = [0.0_dp, spectrum%total_flux(), spectrum%mean_energy(), 0.0_dp, 0.0_dp]
  call spectrum%density(0.0_dp, 0.0_dp, 0.0_dp, me, free, returning, upper_side=.true.)
  expected = 2.0_dp*(ni - free)/(1.0_dp + erf(0.2_dp))
  call spectral%evaluate(input, response, status, message, diagnostic)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'zero-field spectral response: '//trim(message))
  if (status == matching_plane_zhao_ok) then
    call assert_close_dp(diagnostic%ambient_electron_density_m3, expected, 1.0e-10_dp*ni, &
                         'zero-field mixture satisfies actual density neutrality')
  end if
  call test_end()

  call test_begin('allocated_zero_spectrum_preserves_the_no_PE_C_solution')
  call spectrum%clear()
  allocate (spectrum%flux(0))
  call initialize(spectral, 'C', 0.2_dp)
  call initialize(moments, 'C', 0.2_dp)
  call spectral%set_photoelectron_spectrum(spectrum)
  input = [-0.3_dp*d_unit, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  call spectral%evaluate(input, response, status, message, diagnostic)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'zero PE spectral C: '//trim(message))
  call moments%evaluate(input, reference, status, message, reference_diagnostic)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'zero PE moment C: '//trim(message))
  call assert_close_dp(response(1), reference(1), 1.0e-7_dp*te, 'zero PE C response unchanged')
  call test_end()

  call test_begin('model_reinitialization_clears_previously_attached_spectrum')
  call spectrum%bootstrap_maxwellian(0.1_dp*flux_unit, 0.2_dp*te)
  call spectral%set_photoelectron_spectrum(spectrum)
  call spectral%initialize('b', 'require_unique', ni, te, 0.0_dp, mach*sqrt(qe*te/mi), &
                           mi, me, 0.2_dp*te, status, message)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'model reinitialization')
  input = 0.0_dp
  call spectral%evaluate(input, response, status, message, diagnostic)
  call assert_equal_i32(status, matching_plane_zhao_ok, 'reinitialized zero PE response: '//trim(message))
  call assert_close_dp(diagnostic%ambient_electron_density_m3, 2.0_dp*ni, 1.0e-8_dp*ni, &
                       'previously attached spectrum cannot leak into reinitialized model')
  call test_end()

  call test_summary()

contains

  subroutine initialize(model, selected_branch, u, policy)
    type(matching_plane_zhao_model_type), intent(inout) :: model
    character(len=*), intent(in) :: selected_branch
    real(dp), intent(in) :: u
    character(len=*), intent(in), optional :: policy
    character(len=16) :: selected_policy

    selected_policy = 'require_unique'
    if (present(policy)) selected_policy = policy
    call model%initialize(selected_branch, trim(selected_policy), ni, te, u*sqrt(2.0_dp*qe*te/me), &
                          mach*sqrt(qe*te/mi), mi, me, 0.2_dp*te, status, message)
    call assert_equal_i32(status, matching_plane_zhao_ok, 'model initialization: '//trim(message))
    call model%set_photoelectron_spectrum(empty)
  end subroutine initialize

  subroutine check_orbit_closure(branch, h, minimum, amplitude, field, mach, u, fluxes, temperatures)
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: h, minimum, amplitude, field, mach, u, fluxes(2), temperatures(2)
    real(dp) :: rho, integral, actual_field_squared, lowest, probe
    integer :: sample

    rho = oracle_rho(0.0_dp, branch, 'upper', h, minimum, amplitude, mach, u, fluxes, temperatures)
    call assert_close_dp(rho, 0.0_dp, 4.0e-5_dp, 'independent far-boundary neutrality')
    if (branch == 'A') then
      integral = rho_integral(minimum, h, branch, 'lower', h, minimum, amplitude, mach, u, fluxes, temperatures)
      actual_field_squared = -2.0_dp*integral
      integral = rho_integral(minimum, 0.0_dp, branch, 'upper', h, minimum, amplitude, mach, u, fluxes, temperatures)
      call assert_close_dp(integral, 0.0_dp, 4.0e-5_dp, 'independent upper-side zero-field closure')
    else
      integral = rho_integral(h, 0.0_dp, branch, 'monotonic', h, minimum, amplitude, mach, u, fluxes, temperatures)
      actual_field_squared = 2.0_dp*integral
    end if
    call assert_close_dp(actual_field_squared, field*field, 8.0e-5_dp, 'independent interface field closure')
    lowest = 0.0_dp
    do sample = 0, 24
      if (branch == 'A') then
        probe = minimum + (h - minimum)*real(sample, dp)/24.0_dp
        integral = rho_integral(minimum, probe, branch, 'lower', h, minimum, amplitude, mach, u, fluxes, temperatures)
        lowest = min(lowest, -2.0_dp*integral)
        probe = minimum*(1.0_dp - real(sample, dp)/24.0_dp)
        integral = rho_integral(minimum, probe, branch, 'upper', h, minimum, amplitude, mach, u, fluxes, temperatures)
        lowest = min(lowest, -2.0_dp*integral)
      else
        probe = h*(1.0_dp - real(sample, dp)/24.0_dp)
        integral = rho_integral(probe, 0.0_dp, branch, 'monotonic', h, minimum, amplitude, mach, u, fluxes, temperatures)
        lowest = min(lowest, 2.0_dp*integral)
      end if
    end do
    call assert_true(lowest > -8.0e-5_dp, 'independent sampled profile has real electric field')
  end subroutine check_orbit_closure

  ! Independent continuum oracle: two exponential flux laws, not histogram code.
  real(dp) function oracle_rho(phi, branch, side, h, minimum, amplitude, mach, u, fluxes, temperatures) result(rho)
    real(dp), intent(in) :: phi, h, minimum, amplitude, mach, u, fluxes(2), temperatures(2)
    character(len=1), intent(in) :: branch
    character(len=*), intent(in) :: side
    real(dp) :: cutoff, electrons, photoelectrons, prefactor, photo_cutoff
    integer :: component

    cutoff = sqrt(max(phi - minimum, 0.0_dp))
    electrons = 0.5_dp*amplitude*exp(phi)*erfc(cutoff - u)
    if (branch == 'C' .or. (branch == 'A' .and. side == 'upper')) &
      electrons = electrons + amplitude*exp(phi)*(erf(cutoff - u) + erf(u))
    photoelectrons = 0.0_dp
    do component = 1, 2
      prefactor = fluxes(component)*sqrt(pi/(2.0_dp*temperatures(component)))*exp((phi - h)/temperatures(component))
      photo_cutoff = sqrt(max(phi - minimum, 0.0_dp)/temperatures(component))
      photoelectrons = photoelectrons + prefactor*erfc(photo_cutoff)
      if (branch == 'B' .or. (branch == 'A' .and. side == 'lower')) &
        photoelectrons = photoelectrons + 2.0_dp*prefactor*erf(photo_cutoff)
    end do
    rho = 1.0_dp/sqrt(1.0_dp - 2.0_dp*phi/(mach*mach)) - electrons - photoelectrons
  end function oracle_rho

  real(dp) function rho_integral(lo, hi, branch, side, h, minimum, amplitude, mach, u, fluxes, temperatures) result(value)
    real(dp), intent(in) :: lo, hi, h, minimum, amplitude, mach, u, fluxes(2), temperatures(2)
    character(len=1), intent(in) :: branch
    character(len=*), intent(in) :: side
    integer, parameter :: panels = 512
    integer :: point
    real(dp) :: t, phi, weight, jacobian

    value = 0.0_dp
    do point = 0, panels
      t = real(point, dp)/real(panels, dp)
      phi = lo + (hi - lo)*t*t
      jacobian = 2.0_dp*(hi - lo)*t
      weight = 2.0_dp
      if (mod(point, 2) == 1) weight = 4.0_dp
      if (point == 0 .or. point == panels) weight = 1.0_dp
      value = value + weight*jacobian*oracle_rho(phi, branch, side, h, minimum, amplitude, mach, u, fluxes, temperatures)
    end do
    value = value/(3.0_dp*real(panels, dp))
  end function rho_integral

end program test_matching_plane_pe_spectrum

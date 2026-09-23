!> Spectrum closure の根・中性条件・接続場を独立な Maxwell 混合分布で検証する。
program test_matching_plane_pe_spectrum
  use bem_kinds, only: dp, i32
  use bem_constants, only: qe, eps0, pi
  use bem_pe_spectrum, only: pe_spectrum_type
  use bem_app_config, only: app_config, default_app_config, species_from_defaults
  use bem_matching_plane_response_provider, only: matching_plane_response_provider_type, &
                                                  matching_plane_provider_ok, matching_plane_provider_ambiguous_solution
  use bem_matching_plane_implicit, only: solve_matching_implicit_zero_mode
  use bem_mpi, only: mpi_context
  use bem_matching_plane_zhao, only: matching_plane_zhao_model_type, matching_plane_zhao_diagnostics_type, &
                                     matching_plane_zhao_root_seed_type, matching_plane_zhao_ok, &
                                     matching_plane_zhao_no_physical_solution
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
  integer(i32) :: status
  integer :: branch_index, drift_index
  character(len=1) :: branch
  character(len=512) :: message

  call test_init(8)
  spectrum%energy_scale_ev = 0.1_dp*te
  spectrum%bins_per_decade = 1024_i32
  second%energy_scale_ev = spectrum%energy_scale_ev
  second%bins_per_decade = spectrum%bins_per_decade

  call test_begin('maxwell_spectrum_reproduces_admissible_B_and_C_moment_roots')
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
      if (branch == 'C' .and. drift > 0.0_dp) then
        call assert_equal_i32(status, matching_plane_zhao_no_physical_solution, 'drifting C upstream obstruction')
        cycle
      end if
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
  seed%branch = 'A'
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

  call test_begin('finite_drift_A_rejects_the_nonphysical_neutral_infinity_connection')
  call spectrum%bootstrap_maxwellian(0.4_dp*flux_unit, 0.2_dp*te)
  call initialize(spectral, 'A', 0.2_dp)
  call initialize(moments, 'A', 0.2_dp)
  call spectral%set_photoelectron_spectrum(spectrum)
  input = [0.5_dp*d_unit, spectrum%total_flux(), spectrum%mean_energy(), 0.0_dp, 0.0_dp]
  call moments%evaluate(input, reference, status, message, reference_diagnostic)
  call assert_equal_i32(status, matching_plane_zhao_no_physical_solution, 'drifting moment A upstream obstruction')
  call spectral%evaluate(input, response, status, message, diagnostic)
  call assert_equal_i32(status, matching_plane_zhao_no_physical_solution, 'drifting spectral A upstream obstruction')
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
  call initialize(spectral, 'C', 0.0_dp)
  call initialize(moments, 'C', 0.0_dp)
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

  call test_begin('captured_spectrum_implicit_endpoints_preserve_ambiguity_and_continuation')
  call check_captured_implicit_endpoints()
  call test_end()

  call test_summary()

contains

  subroutine check_captured_implicit_endpoints()
    ! First-batch replay spectrum captured on 2026-09-23.  Independent SciPy
    ! upstream-velocity and Poisson quadrature found these three BE endpoints.
    real(dp), parameter :: root_phi(3) = [5.829194512583357_dp, 5.945909126905672_dp, 6.055444426356461_dp]
    real(dp), parameter :: root_density(3) = [4425184.213264361_dp, 4380830.616359554_dp, 4321700.853467314_dp]
    real(dp), parameter :: root_displacement(3) = [ &
                           1.769541040692815e-11_dp, 1.7696825059325093e-11_dp, 1.7704024051544765e-11_dp]
    real(dp), parameter :: captured_me = 9.1093837139e-31_dp
    type(app_config) :: cfg
    type(matching_plane_response_provider_type) :: provider
    type(matching_plane_zhao_root_seed_type) :: previous, candidate
    type(pe_spectrum_type) :: captured
    type(mpi_context) :: serial
    real(dp) :: feedback(4), output(6), displacement, lower_edge, upper_edge, net_current
    integer :: unit_id, ios, bin_index, root_index
    integer(i32) :: provider_status
    logical :: handled
    character(len=512) :: provider_message, header

    captured%energy_scale_ev = 2.19999999970613391_dp
    captured%bins_per_decade = 128_i32
    allocate (captured%flux(207))
    open (newunit=unit_id, file='tests/fixtures/matching_plane_pe_captured_spectrum.csv', &
          status='old', action='read', iostat=ios)
    call assert_true(ios == 0, 'captured spectrum fixture is available')
    if (ios /= 0) return
    read (unit_id, '(a)', iostat=ios) header
    call assert_true(ios == 0 .and. trim(header) == 'energy_low_ev,energy_high_ev,flux_m2_s', 'captured spectrum header')
    do bin_index = 1, size(captured%flux)
      read (unit_id, *, iostat=ios) lower_edge, upper_edge, captured%flux(bin_index)
      call assert_true(ios == 0, 'captured spectrum bin is readable')
      if (ios /= 0) then
        close (unit_id)
        return
      end if
      call assert_close_dp(lower_edge, captured%edge(bin_index - 1), 1.0e-12_dp, 'captured lower energy edge')
      call assert_close_dp(upper_edge, captured%edge(bin_index), 1.0e-12_dp, 'captured upper energy edge')
    end do
    close (unit_id)
    call assert_close_dp(captured%total_flux(), 1.6662817804987412e13_dp, 1.0_dp, 'captured outward flux')
    feedback = [captured%total_flux(), captured%mean_energy(), 0.0_dp, 0.0_dp]

    call default_app_config(cfg)
    cfg%surface_current%model = 'matching_plane_quasistatic'
    cfg%surface_current%response_backend = 'zhao_online'
    cfg%surface_current%zhao_branch = 'auto'
    cfg%surface_current%zhao_root_selection = 'continuation'
    cfg%surface_current%photoelectron_closure = 'energy_spectrum'
    cfg%surface_current%electron_species = 'electron'
    cfg%surface_current%ion_species = 'ion'
    cfg%surface_current%photoelectron_species = 'photoelectron'
    cfg%n_particle_species = 3_i32
    cfg%particle_species(1:3) = species_from_defaults()
    cfg%particle_species(1)%species_key = 'electron'
    cfg%particle_species(1)%q_particle = -qe
    cfg%particle_species(1)%m_particle = captured_me
    cfg%particle_species(1)%temperature_ev = 10.0_dp
    cfg%particle_species(1)%has_temperature_ev = .true.
    cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -4.0e5_dp]
    cfg%particle_species(2)%species_key = 'ion'
    cfg%particle_species(2)%q_particle = qe
    cfg%particle_species(2)%m_particle = mi
    cfg%particle_species(2)%number_density_m3 = 5.0e6_dp
    cfg%particle_species(2)%drift_velocity = [0.0_dp, 0.0_dp, -4.0e5_dp]
    cfg%particle_species(3)%species_key = 'photoelectron'
    cfg%particle_species(3)%q_particle = -qe
    cfg%particle_species(3)%m_particle = captured_me
    cfg%particle_species(3)%temperature_ev = captured%energy_scale_ev
    cfg%particle_species(3)%has_temperature_ev = .true.
    serial = mpi_context()
    call provider%initialize(cfg, serial, provider_status, provider_message)
    call assert_equal_i32(provider_status, matching_plane_provider_ok, 'captured provider initialization: '// &
                          trim(provider_message))
    if (provider_status /= matching_plane_provider_ok) return
    call provider%set_photoelectron_spectrum(captured)
    previous = matching_plane_zhao_root_seed_type()
    call provider%solve_implicit_endpoint(feedback, 0.0_dp, 2.0_dp, -qe, qe, .true., -qe, &
                                          previous, handled, displacement, output, candidate, provider_status, provider_message)
    call assert_true(handled, 'positive-drift online implicit endpoint has a dedicated solve')
    call assert_equal_i32(provider_status, matching_plane_provider_ambiguous_solution, &
                          'three physical endpoints without a seed remain ambiguous: '//trim(provider_message))
    call assert_true(.not. candidate%valid, 'ambiguous endpoint cannot publish a continuation seed')

    do root_index = 1, size(root_phi)
      previous%valid = .true.
      previous%branch = 'B'
      previous%phi0_v = root_phi(root_index) + 0.01_dp
      previous%phi_m_v = previous%phi0_v
      previous%ambient_electron_density_m3 = root_density(root_index)
      call provider%solve_implicit_endpoint(feedback, 0.0_dp, 2.0_dp, -qe, qe, .true., -qe, &
                                            previous, handled, displacement, output, candidate, provider_status, provider_message)
      call assert_equal_i32(provider_status, matching_plane_provider_ok, 'captured seeded endpoint: '//trim(provider_message))
      if (provider_status /= matching_plane_provider_ok) cycle
      call assert_true(handled .and. candidate%valid .and. candidate%branch == 'B', 'certified Type-B endpoint seed')
      call assert_close_dp(output(1), root_phi(root_index), 5.0e-6_dp, 'independent implicit interface potential')
      call assert_close_dp(candidate%ambient_electron_density_m3, root_density(root_index), 10.0_dp, &
                           'independent implicit electron normalization')
      call assert_close_dp(displacement, root_displacement(root_index), 5.0e-17_dp, 'independent implicit displacement')
      net_current = qe*(output(3) - output(2) + captured%tail_flux(output(1)))
      call assert_close_dp(displacement - 2.0_dp*net_current, 0.0_dp, 2.0e-17_dp, 'backward Euler charge balance')
    end do

    ! The previous replay seed must select the highest-potential endpoint, even
    ! when the whole physical interval fits between two old displacement probes.
    previous%phi0_v = 6.96756157119758868_dp
    previous%phi_m_v = previous%phi0_v
    previous%ambient_electron_density_m3 = 4.22844877100958023e6_dp
    call solve_matching_implicit_zero_mode( &
      provider, serial, 0.0_dp, 2.43269927022345257e-11_dp, 2.0_dp, .false., 0.0_dp, 0.0_dp, &
      8.42198694632710321e-12_dp, 0_i32, feedback, -qe, qe, .true., -qe, previous, candidate, displacement, output)
    call assert_close_dp(output(1), root_phi(3), 5.0e-6_dp, 'outer implicit solve preserves continuation selection')
    call assert_close_dp(displacement, root_displacement(3), 5.0e-17_dp, 'outer implicit displacement')
    call assert_true(candidate%valid, 'outer implicit solve publishes its certified seed')
  end subroutine check_captured_implicit_endpoints

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
    ! The drifted incoming distribution is specified at infinity.  Integrate
    ! in upstream velocity, independently of the production local quadrature.
    if (u /= 0.0_dp) electrons = amplitude*upstream_density(phi, minimum, u)
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

  real(dp) function upstream_density(phi, minimum, u) result(value)
    real(dp), intent(in) :: phi, minimum, u
    integer, parameter :: panels = 1024
    integer :: point
    real(dp) :: a, step, weight, factor

    step = 12.0_dp/real(panels, dp)
    value = 0.0_dp
    do point = 0, panels
      a = sqrt(max(-minimum, 0.0_dp)) + real(point, dp)*step
      factor = 1.0_dp
      if (a*a + phi > 0.0_dp) factor = a/sqrt(a*a + phi)
      weight = 2.0_dp
      if (mod(point, 2) == 1) weight = 4.0_dp
      if (point == 0 .or. point == panels) weight = 1.0_dp
      value = value + weight*factor*exp(-(a - u)**2)
    end do
    value = value*step/(3.0_dp*sqrt(pi))
  end function upstream_density

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

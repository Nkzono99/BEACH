program test_outer_plasma_zhao
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, qe
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok, &
                                    outer_plasma_invalid, outer_plasma_no_physical_solution, &
                                    outer_plasma_numerical_failure
  use bem_sheath_model_core, only: zhao_params_type, build_zhao_params, solve_zhao_unknowns, &
                                   evaluate_zhao_density_hat
  use bem_outer_plasma_zhao, only: zhao_charge_root_type, solve_zhao_charge_root, &
                                   evaluate_zhao_interface_field, build_zhao_outer_profile, &
                                   solve_outer_plasma_zhao_column
  use bem_outer_plasma_kinetic, only: kinetic_outer_plasma_options_type, solve_outer_plasma_kinetic
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, &
                          assert_close_dp, assert_equal_i32
  implicit none

  integer(i32), parameter :: profile_points = 257_i32
  real(dp), parameter :: electron_mass = 9.1093837015e-31_dp
  real(dp), parameter :: proton_mass = 1.67262192369e-27_dp
  real(dp), parameter :: population_gap_field = 9.2776117384997914e-2_dp
  real(dp), parameter :: runtime_gap_field = 4.1743975711591802e-1_dp
  real(dp), parameter :: runtime_previous_column = 7.22108e7_dp
  real(dp), parameter :: runtime_target_column = 7.56494e7_dp
  real(dp), parameter :: control_length_m = 2.0_dp
  real(dp), parameter :: feasible_population_fractions(3) = [0.0_dp, 0.1_dp, 0.25_dp]
  type(zhao_params_type) :: params
  type(kinetic_outer_plasma_options_type) :: kinetic_options
  type(zhao_charge_root_type) :: stationary, charge_root, perturbed_root, zero_root, reference_root
  type(zhao_charge_root_type) :: continuation_target_root, field_target_root
  type(outer_plasma_state_type) :: state, reference_state, continuation_target_state
  type(outer_plasma_state_type) :: field_target_state
  real(dp) :: equilibrium_field, perturbed_field, gauss_expected, gauss_scale, reference_column
  real(dp) :: continuation_target_column, field_target, field_target_column
  real(dp) :: manually_integrated_column, type_a_minimum_z, type_a_native_extent
  real(dp) :: phi0_v, phi_m_v, density_m3
  integer(i32) :: status
  character(len=1) :: branch
  character(len=256) :: message
  integer :: fraction_index

  call test_init(14)

  call test_begin('stationary Zhao-A root is recovered by its interface field')
  call configure_params(60.0_dp, params)
  call solve_zhao_unknowns('zhao_a', params, phi0_v, phi_m_v, density_m3, branch)
  stationary = zhao_charge_root_type( &
               branch=branch, phi0_v=phi0_v, phi_m_v=phi_m_v, n_swe_inf_m3=density_m3 &
               )
  call evaluate_zhao_interface_field(params, stationary, equilibrium_field, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'stationary Zhao-A field evaluation failed')
  call assert_true(equilibrium_field > 0.0_dp, 'stationary Zhao-A lower field must be positive')
  call solve_zhao_charge_root( &
    'zhao_a', params, equilibrium_field, charge_root, status, message, initial_root=stationary &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'charge-driven Zhao-A solve failed: '//trim(message))
  call assert_close_dp(charge_root%phi0_v, phi0_v, 2.0e-5_dp, 'Zhao-A phi0 stationary recovery mismatch')
  call assert_close_dp(charge_root%phi_m_v, phi_m_v, 2.0e-5_dp, 'Zhao-A phi_m stationary recovery mismatch')
  call assert_close_dp( &
    charge_root%n_swe_inf_m3, density_m3, 2.0e2_dp, 'Zhao-A density stationary recovery mismatch' &
    )
  call assert_true( &
    abs(charge_root%net_current_density_a_m2) < 1.0e-11_dp, &
    'stationary Zhao-A net current is not zero' &
    )
  call build_zhao_outer_profile(params, charge_root, profile_points, state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'Zhao-A profile reconstruction failed: '//trim(message))
  call assert_true(state%ready .and. state%profile_n == profile_points, 'Zhao-A profile state is not ready')
  call assert_true( &
    minval(state%potential) < 0.0_dp .and. state%potential(1) > 0.0_dp, &
    'Zhao-A profile must contain a negative minimum below its positive interface' &
    )
  gauss_expected = eps0*(state%field(state%profile_n) - state%interface_field)
  gauss_scale = max(abs(gauss_expected), abs(eps0*state%interface_field), 1.0e-20_dp)
  call assert_true( &
    abs(state%integrated_charge_per_area - gauss_expected) <= 5.0e-2_dp*gauss_scale, &
    'Zhao-A finite-profile Gauss closure mismatch' &
    )
  type_a_minimum_z = state%z(minloc(state%potential, dim=1))
  type_a_native_extent = state%z(state%profile_n)
  call build_zhao_outer_profile( &
    params, charge_root, profile_points, state, status, message, control_length_m=type_a_native_extent &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'finite Zhao-A profile reconstruction failed: '//trim(message))
  call assert_close_dp( &
    minval(state%potential), charge_root%phi_m_v, 1.0e-12_dp, &
    'finite Zhao-A profile must retain its resolved potential minimum' &
    )
  call build_zhao_outer_profile( &
    params, charge_root, profile_points, state, status, message, control_length_m=0.5_dp*type_a_minimum_z &
    )
  call assert_equal_i32(status, outer_plasma_no_physical_solution, &
                        'a finite Zhao-A profile ending before its minimum must fail closed')
  call test_end()

  call test_begin('field perturbation leaves a nonzero Zhao-A charging current')
  perturbed_field = 1.01_dp*equilibrium_field
  call solve_zhao_charge_root( &
    'zhao_a', params, perturbed_field, perturbed_root, status, message, initial_root=charge_root &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'perturbed Zhao-A solve failed: '//trim(message))
  call assert_true( &
    abs(perturbed_root%net_current_density_a_m2) > 1.0e-14_dp, &
    'perturbed Zhao-A state must retain a nonzero current diagnostic' &
    )
  call test_end()

  call test_begin('auto continuation preserves a compatible previous branch')
  call configure_params(60.0_dp, params)
  call solve_zhao_unknowns('zhao_b', params, phi0_v, phi_m_v, density_m3, branch)
  stationary = zhao_charge_root_type( &
               branch=branch, phi0_v=phi0_v, phi_m_v=phi_m_v, n_swe_inf_m3=density_m3 &
               )
  call evaluate_zhao_interface_field(params, stationary, equilibrium_field, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'stationary Zhao-B field evaluation failed')
  call solve_zhao_charge_root( &
    'auto', params, equilibrium_field, charge_root, status, message, initial_root=stationary &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'auto Zhao-B continuation failed: '//trim(message))
  call assert_true(charge_root%branch == 'B', 'auto continuation changed a compatible Zhao-B branch')
  call test_end()

  call test_begin('stationary Zhao-C root is recovered by its interface field')
  call configure_params(10.0_dp, params)
  call solve_zhao_unknowns('zhao_c', params, phi0_v, phi_m_v, density_m3, branch)
  stationary = zhao_charge_root_type( &
               branch=branch, phi0_v=phi0_v, phi_m_v=phi_m_v, n_swe_inf_m3=density_m3 &
               )
  call evaluate_zhao_interface_field(params, stationary, equilibrium_field, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'stationary Zhao-C field evaluation failed')
  call assert_true(equilibrium_field < 0.0_dp, 'stationary Zhao-C field must be negative')
  call solve_zhao_charge_root( &
    'zhao_c', params, equilibrium_field, charge_root, status, message, initial_root=stationary &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'charge-driven Zhao-C solve failed: '//trim(message))
  call assert_close_dp(charge_root%phi0_v, phi0_v, 2.0e-5_dp, 'Zhao-C phi0 stationary recovery mismatch')
  call assert_close_dp( &
    charge_root%n_swe_inf_m3, density_m3, 2.0e2_dp, 'Zhao-C density stationary recovery mismatch' &
    )
  call assert_true( &
    abs(charge_root%net_current_density_a_m2) < 1.0e-11_dp, &
    'stationary Zhao-C net current is not zero' &
    )
  call build_zhao_outer_profile(params, charge_root, profile_points, state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'Zhao-C profile reconstruction failed: '//trim(message))
  call assert_true(all(state%potential <= 1.0e-12_dp), 'Zhao-C profile must stay on its negative branch')
  gauss_expected = eps0*(state%field(state%profile_n) - state%interface_field)
  gauss_scale = max(abs(gauss_expected), abs(eps0*state%interface_field), 1.0e-20_dp)
  call assert_true( &
    abs(state%integrated_charge_per_area - gauss_expected) <= 5.0e-2_dp*gauss_scale, &
    'Zhao-C finite-profile Gauss closure mismatch' &
    )
  call test_end()

  call test_begin('branch-incompatible prescribed fields fail physically')
  call solve_zhao_charge_root('zhao_c', params, abs(equilibrium_field), charge_root, status, message)
  call assert_equal_i32( &
    status, outer_plasma_no_physical_solution, 'positive field must not select explicit Zhao-C' &
    )
  call configure_params(60.0_dp, params)
  call solve_zhao_charge_root('zhao_a', params, -abs(equilibrium_field), charge_root, status, message)
  call assert_equal_i32( &
    status, outer_plasma_no_physical_solution, 'negative field must not select explicit Zhao-A lower side' &
    )
  call test_end()

  call test_begin('partial photoelectron population closes the low-field gap')
  call configure_params(60.0_dp, params)
  do fraction_index = 1, size(feasible_population_fractions)
    params%photoelectron_population_fraction = feasible_population_fractions(fraction_index)
    call solve_zhao_charge_root( &
      'zhao_b', params, population_gap_field, charge_root, status, message &
      )
    call assert_equal_i32( &
      status, outer_plasma_ok, 'partial-population Zhao-B solve failed: '//trim(message) &
      )
    call assert_true(charge_root%branch == 'B', 'partial-population solve must remain on Zhao-B')
    call assert_close_dp( &
      charge_root%photoelectron_population_fraction, feasible_population_fractions(fraction_index), &
      1.0e-14_dp, 'Zhao root lost its photoelectron population fraction' &
      )
    call build_zhao_outer_profile(params, charge_root, profile_points, state, status, message)
    call assert_equal_i32( &
      status, outer_plasma_ok, 'partial-population profile reconstruction failed: '//trim(message) &
      )
    call assert_close_dp( &
      state%photoelectron_population_fraction, feasible_population_fractions(fraction_index), &
      1.0e-14_dp, 'Zhao state lost its photoelectron population fraction' &
      )
    call assert_true( &
      state%photoelectron_current_density > 0.0_dp, &
      'partial outer population must not scale the tracked full emission current' &
      )
  end do
  params%photoelectron_population_fraction = 0.0_dp
  call solve_zhao_charge_root('zhao_auto', params, 0.0_dp, zero_root, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'ambient-only zero-field Zhao state failed: '//trim(message))
  call assert_true(zero_root%branch == 'B', 'ambient-only zero-field state must not use transient branch 0')
  params%photoelectron_population_fraction = 1.0_dp
  call solve_zhao_charge_root('zhao_auto', params, population_gap_field, charge_root, status, message)
  call assert_equal_i32( &
    status, outer_plasma_numerical_failure, &
    'full-photoelectron population must preserve the unresolved Newton failure' &
    )
  call assert_true(charge_root%branch /= '0', 'nonzero field must never fall back to transient branch 0')
  call test_end()

  call test_begin('zero-field state provides a finite transient bootstrap')
  call solve_zhao_charge_root('zhao_auto', params, 0.0_dp, zero_root, status, message)
  call assert_equal_i32( &
    status, outer_plasma_no_physical_solution, &
    'zero-field Zhao bootstrap must require explicit permission' &
    )
  call assert_true(zero_root%branch /= '0', 'automatic Zhao solve created an implicit transient branch')
  call solve_zhao_charge_root( &
    'zhao_a', params, 0.0_dp, zero_root, status, message, allow_transient_bootstrap=.true. &
    )
  call assert_equal_i32( &
    status, outer_plasma_no_physical_solution, 'explicit Zhao-A was hijacked by transient bootstrap' &
    )
  call assert_true(zero_root%branch /= '0', 'explicit Zhao-A must never select transient branch 0')
  call solve_zhao_charge_root( &
    'zhao_b', params, 0.0_dp, zero_root, status, message, allow_transient_bootstrap=.true. &
    )
  call assert_equal_i32( &
    status, outer_plasma_no_physical_solution, 'explicit Zhao-B was hijacked by transient bootstrap' &
    )
  call assert_true(zero_root%branch /= '0', 'explicit Zhao-B must never select transient branch 0')
  call solve_zhao_charge_root( &
    'zhao_c', params, 0.0_dp, zero_root, status, message, allow_transient_bootstrap=.true. &
    )
  call assert_equal_i32( &
    status, outer_plasma_no_physical_solution, 'explicit Zhao-C was hijacked by transient bootstrap' &
    )
  call assert_true(zero_root%branch /= '0', 'explicit Zhao-C must never select transient branch 0')
  call solve_zhao_charge_root( &
    'zhao_auto', params, 0.0_dp, zero_root, status, message, &
    allow_transient_bootstrap=.true. &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'zero-field Zhao bootstrap failed: '//trim(message))
  call assert_true(zero_root%branch == '0', 'high-photoemission zero-field state should use transient bootstrap')
  call assert_true( &
    abs(zero_root%net_current_density_a_m2) > 0.0_dp, &
    'zero-field transient bootstrap must retain its charging current' &
    )
  call build_zhao_outer_profile(params, zero_root, profile_points, state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'zero-field profile build failed: '//trim(message))
  call assert_true(maxval(abs(state%potential)) == 0.0_dp, 'zero-field bootstrap potential must be flat')
  call assert_true(maxval(abs(state%charge_density)) == 0.0_dp, 'zero-field bootstrap charge must be flat')
  call solve_outer_plasma_zhao_column( &
    'auto', params, 0.0_dp, profile_points, control_length_m, 0.0_dp, &
    state, charge_root, status, message, initial_root=zero_root &
    )
  call assert_equal_i32( &
    status, outer_plasma_invalid, 'column closure must reject a transient branch-0 anchor' &
    )
  call assert_true(len_trim(message) > 0, 'branch-0 continuation rejection needs a diagnostic message')
  call test_end()

  call test_begin('resolved profile stores the free plus captured photoelectron column')
  call configure_params(60.0_dp, params)
  params%photoelectron_population_fraction = 0.1_dp
  call solve_zhao_charge_root( &
    'zhao_b', params, population_gap_field, reference_root, status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'column reference root failed: '//trim(message))
  call build_zhao_outer_profile( &
    params, reference_root, profile_points, reference_state, status, message, &
    control_length_m=control_length_m &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'finite-column profile failed: '//trim(message))
  call assert_close_dp( &
    reference_state%z(reference_state%profile_n), control_length_m, 1.0e-14_dp, &
    'finite-column profile does not end at the requested control length' &
    )
  manually_integrated_column = resolved_photoelectron_column(params, reference_root, reference_state)
  call assert_close_dp( &
    reference_state%photoelectron_column_per_area, manually_integrated_column, &
    1.0e-12_dp*max(1.0_dp, manually_integrated_column), &
    'stored photoelectron column is not the resolved-profile trapezoid integral' &
    )
  reference_column = reference_state%photoelectron_column_per_area
  call assert_true(reference_column > 0.0_dp, 'partial Zhao profile must contain photoelectrons')
  call test_end()

  call test_begin('column closure recovers eta from a finite-volume inventory')
  params%photoelectron_population_fraction = 1.0_dp
  call solve_outer_plasma_zhao_column( &
    'auto', params, population_gap_field, profile_points, control_length_m, reference_column, &
    state, charge_root, status, message, initial_root=reference_root &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'photoelectron-column closure failed: '//trim(message))
  call assert_close_dp( &
    charge_root%photoelectron_population_fraction, 0.1_dp, 1.0e-12_dp, &
    'photoelectron-column closure did not recover its reference eta' &
    )
  call assert_close_dp( &
    state%photoelectron_column_per_area, reference_column, &
    1.0e-8_dp*max(1.0_dp, reference_column), &
    'photoelectron-column closure did not reproduce the queue inventory' &
    )
  call assert_close_dp( &
    state%photoelectron_column_residual_per_area, 0.0_dp, &
    1.0e-8_dp*max(1.0_dp, reference_column), &
    'photoelectron-column closure residual is too large' &
    )
  call test_end()

  call test_begin('auto column closure follows the previous branch in both eta directions')
  params%photoelectron_population_fraction = 0.05_dp
  call solve_zhao_charge_root( &
    'zhao_b', params, population_gap_field, continuation_target_root, status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'continuation target root failed: '//trim(message))
  call build_zhao_outer_profile( &
    params, continuation_target_root, profile_points, continuation_target_state, status, message, &
    control_length_m=control_length_m &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'continuation target profile failed: '//trim(message))
  continuation_target_column = continuation_target_state%photoelectron_column_per_area
  params%photoelectron_population_fraction = 1.0_dp
  call solve_outer_plasma_zhao_column( &
    'auto', params, population_gap_field, profile_points, control_length_m, &
    continuation_target_column, state, charge_root, status, message, initial_root=reference_root &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'auto previous-root column solve failed: '//trim(message))
  call assert_true(charge_root%branch == 'B', 'auto column solve left the previous Zhao-B branch')
  call assert_close_dp( &
    charge_root%photoelectron_population_fraction, 0.05_dp, 1.0e-12_dp, &
    'auto column solve did not follow the previous root downward in eta' &
    )
  call assert_close_dp( &
    state%photoelectron_column_per_area, continuation_target_column, &
    1.0e-8_dp*max(1.0_dp, continuation_target_column), &
    'auto previous-root solve did not reproduce its target inventory' &
    )
  call solve_outer_plasma_zhao_column( &
    'zhao_b', params, population_gap_field, profile_points, control_length_m, &
    continuation_target_column, state, charge_root, status, message, initial_root=reference_root &
    )
  call assert_equal_i32( &
    status, outer_plasma_invalid, 'column closure must reject a fixed Zhao branch' &
    )
  call assert_true(len_trim(message) > 0, 'fixed-branch column rejection needs a diagnostic message')
  call solve_outer_plasma_zhao_column( &
    'auto', params, population_gap_field, profile_points, control_length_m, reference_column, &
    state, charge_root, status, message, initial_root=continuation_target_root &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'auto upward column solve failed: '//trim(message))
  call assert_true(charge_root%branch == 'B', 'upward auto column solve left Zhao-B')
  call assert_close_dp( &
    charge_root%photoelectron_population_fraction, 0.1_dp, 1.0e-12_dp, &
    'auto column solve did not follow the previous root upward in eta' &
    )
  call test_end()

  call test_begin('auto column closure homotopies the previous field before solving eta')
  params%photoelectron_population_fraction = 0.1_dp
  field_target = 1.01_dp*population_gap_field
  call solve_zhao_charge_root( &
    'zhao_b', params, field_target, field_target_root, status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'field-homotopy target root failed: '//trim(message))
  call build_zhao_outer_profile( &
    params, field_target_root, profile_points, field_target_state, status, message, &
    control_length_m=control_length_m &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'field-homotopy target profile failed: '//trim(message))
  field_target_column = field_target_state%photoelectron_column_per_area
  params%photoelectron_population_fraction = 1.0_dp
  call solve_outer_plasma_zhao_column( &
    'auto', params, field_target, profile_points, control_length_m, field_target_column, &
    state, charge_root, status, message, initial_root=reference_root &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'auto field homotopy failed: '//trim(message))
  call assert_true(charge_root%branch == 'B', 'field homotopy left the previous Zhao-B branch')
  call assert_close_dp( &
    charge_root%interface_field_v_m, field_target, 1.0e-12_dp, &
    'field homotopy did not reach the requested interface field' &
    )
  call assert_close_dp( &
    charge_root%photoelectron_population_fraction, 0.1_dp, 1.0e-12_dp, &
    'field homotopy changed eta for an already matched inventory' &
    )
  call test_end()

  call test_begin('kinetic options route queue column closure through continuation state')
  kinetic_options = kinetic_outer_plasma_options_type()
  kinetic_options%kinetic_closure = 'zhao_charge_driven'
  kinetic_options%zhao_branch = 'auto'
  kinetic_options%grid_points = profile_points
  kinetic_options%domain_length = control_length_m
  kinetic_options%interface_field = population_gap_field
  kinetic_options%electron_charge = -qe
  kinetic_options%electron_mass = electron_mass
  kinetic_options%electron_temperature_j = 12.0_dp*qe
  kinetic_options%electron_drift_infinity = 4.68e5_dp*sin(60.0_dp*pi_value()/180.0_dp)
  kinetic_options%ion_charge = qe
  kinetic_options%ion_mass = proton_mass
  kinetic_options%ion_density_infinity = 8.7e6_dp
  kinetic_options%ion_drift_infinity = kinetic_options%electron_drift_infinity
  kinetic_options%photoelectron_charge = -qe
  kinetic_options%photoelectron_mass = electron_mass
  kinetic_options%photoelectron_temperature_j = 2.2_dp*qe
  kinetic_options%photoelectron_reference_density = 64.0e6_dp
  kinetic_options%photoelectron_population_fraction = 1.0_dp
  kinetic_options%photoelectron_column_closure_enabled = .true.
  kinetic_options%photoelectron_column_target_m2 = reference_column
  kinetic_options%zhao_alpha_deg = 60.0_dp
  call solve_outer_plasma_kinetic( &
    kinetic_options, state, status, message, initial_state=reference_state &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'kinetic column-closure routing failed: '//trim(message))
  call assert_close_dp( &
    state%photoelectron_population_fraction, 0.1_dp, 1.0e-12_dp, &
    'kinetic continuation did not preserve the column-selected eta' &
    )
  call assert_close_dp( &
    state%z(state%profile_n), control_length_m, 1.0e-14_dp, &
    'kinetic column closure did not retain the finite control length' &
    )
  call test_end()

  call test_begin('column closure treats zero inventory exactly and rejects unreachable inventory')
  call solve_outer_plasma_zhao_column( &
    'auto', params, population_gap_field, profile_points, control_length_m, 0.0_dp, &
    state, charge_root, status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'zero-column closure failed: '//trim(message))
  call assert_close_dp( &
    charge_root%photoelectron_population_fraction, 0.0_dp, 0.0_dp, &
    'zero-column closure must use eta=0 without clamping' &
    )
  call assert_close_dp( &
    state%photoelectron_column_per_area, 0.0_dp, 0.0_dp, &
    'zero-column closure must contain no photoelectrons' &
    )
  call solve_outer_plasma_zhao_column( &
    'auto', params, population_gap_field, profile_points, control_length_m, 1.0e30_dp, &
    state, charge_root, status, message, initial_root=reference_root &
    )
  call assert_equal_i32( &
    status, outer_plasma_no_physical_solution, &
    'unreachable photoelectron inventory must report no physical solution' &
    )
  call assert_true(len_trim(message) > 0, 'unreachable photoelectron inventory needs a diagnostic message')
  call test_end()

  call test_begin('column closure brackets a target before a later Zhao path endpoint')
  call configure_params(60.0_dp, params)
  params%photoelectron_population_fraction = 1.0_dp
  call solve_outer_plasma_zhao_column( &
    'auto', params, runtime_gap_field, 128_i32, 10.0_dp*params%lambda_d_phe_ref_m, &
    runtime_previous_column, reference_state, reference_root, status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'runtime predecessor column solve failed: '//trim(message))
  call solve_outer_plasma_zhao_column( &
    'auto', params, runtime_gap_field, 128_i32, 10.0_dp*params%lambda_d_phe_ref_m, &
    runtime_target_column, state, charge_root, status, message, initial_root=reference_root &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'runtime target column solve failed: '//trim(message))
  call assert_true(charge_root%branch == 'B', 'runtime target solve left the connected Zhao-B branch')
  call assert_true( &
    charge_root%photoelectron_population_fraction > reference_root%photoelectron_population_fraction .and. &
    charge_root%photoelectron_population_fraction < 2.0_dp*reference_root%photoelectron_population_fraction, &
    'runtime target solve searched beyond the first eta bracket' &
    )
  call assert_close_dp( &
    state%photoelectron_column_per_area, runtime_target_column, &
    1.0e-8_dp*runtime_target_column, 'runtime target column residual is too large' &
    )
  call test_end()

  call test_summary()

contains

  subroutine configure_params(alpha_deg, p)
    real(dp), intent(in) :: alpha_deg
    type(zhao_params_type), intent(out) :: p

    real(dp) :: speed

    speed = 4.68e5_dp*sin(alpha_deg*pi_value()/180.0_dp)
    call build_zhao_params( &
      alpha_deg, 8.7e6_dp, 64.0e6_dp, 12.0_dp, 2.2_dp, speed, speed, &
      proton_mass, electron_mass, p &
      )
  end subroutine configure_params

  pure real(dp) function pi_value() result(value)
    value = acos(-1.0_dp)
  end function pi_value

  real(dp) function resolved_photoelectron_column(p, root, profile) result(column)
    type(zhao_params_type), intent(in) :: p
    type(zhao_charge_root_type), intent(in) :: root
    type(outer_plasma_state_type), intent(in) :: profile

    real(dp), allocatable :: density(:)
    real(dp) :: phi_hat, phi0_hat, phi_m_hat, ambient_density_hat
    real(dp) :: n_swi_hat, n_swe_f_hat, n_swe_r_hat, n_phe_f_hat, n_phe_c_hat
    integer :: point, minimum_point
    character(len=16) :: side

    allocate (density(profile%profile_n))
    phi0_hat = root%phi0_v/p%t_phe_ev
    phi_m_hat = root%phi_m_v/p%t_phe_ev
    ambient_density_hat = root%n_swe_inf_m3/p%n_phe_ref_m3
    minimum_point = minloc(profile%potential, dim=1)
    do point = 1, profile%profile_n
      if (root%branch == 'A') then
        if (point <= minimum_point) then
          side = 'lower'
        else
          side = 'upper'
        end if
      else
        side = 'monotonic'
      end if
      phi_hat = profile%potential(point)/p%t_phe_ev
      call evaluate_zhao_density_hat( &
        p, root%branch, side, phi_hat, phi0_hat, phi_m_hat, ambient_density_hat, &
        n_swi_hat, n_swe_f_hat, n_swe_r_hat, n_phe_f_hat, n_phe_c_hat &
        )
      density(point) = p%n_phe_ref_m3*(n_phe_f_hat + n_phe_c_hat)
    end do
    column = 0.0_dp
    do point = 1, profile%profile_n - 1
      column = column + 0.5_dp*(density(point) + density(point + 1))* &
               (profile%z(point + 1) - profile%z(point))
    end do
  end function resolved_photoelectron_column

end program test_outer_plasma_zhao

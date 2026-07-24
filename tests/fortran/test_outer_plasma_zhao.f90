program test_outer_plasma_zhao
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, qe
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok, &
                                    outer_plasma_invalid, outer_plasma_no_physical_solution, &
                                    outer_plasma_numerical_failure
  use bem_sheath_model_core, only: zhao_params_type, build_zhao_params, solve_zhao_unknowns, &
                                   evaluate_zhao_density_hat
  use bem_outer_plasma_zhao, only: zhao_charge_root_type, solve_zhao_charge_root, &
                                   evaluate_zhao_interface_field, build_zhao_outer_profile, &
                                   solve_outer_plasma_zhao_column, &
                                   zhao_continuation_diagnostics_type, &
                                   zhao_continuation_reason_guard_rejected, &
                                   zhao_continuation_reason_numerical_failure, &
                                   zhao_continuation_reason_search_limit, &
                                   zhao_continuation_reason_target_bracketed, &
                                   zhao_branch_atlas_options_type, zhao_branch_atlas_type, &
                                   zhao_field_column_homotopy_options_type, &
                                   zhao_field_column_homotopy_type, &
                                   zhao_ab_degeneracy_diagnostics_type, &
                                   trace_zhao_branch_atlas, write_zhao_continuation_diagnostics, &
                                   trace_zhao_field_column_homotopy, diagnose_zhao_ab_degeneracy
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
  type(zhao_charge_root_type) :: continuation_target_root, field_target_root, nonfinite_root
  type(zhao_continuation_diagnostics_type) :: continuation_diagnostics
  type(zhao_branch_atlas_options_type) :: atlas_options
  type(zhao_branch_atlas_type) :: branch_atlas, coarse_branch_atlas
  type(zhao_branch_atlas_type) :: reverse_branch_atlas, candidate_branch_atlas
  type(zhao_field_column_homotopy_options_type) :: homotopy_options
  type(zhao_field_column_homotopy_type) :: field_column_homotopy
  type(zhao_field_column_homotopy_type) :: coarse_field_column_homotopy
  type(zhao_ab_degeneracy_diagnostics_type) :: ab_diagnostics
  type(outer_plasma_state_type) :: state, reference_state, continuation_target_state
  type(outer_plasma_state_type) :: field_target_state
  real(dp) :: equilibrium_field, perturbed_field, gauss_expected, gauss_scale, reference_column
  real(dp) :: continuation_target_column, field_target, field_target_column
  real(dp) :: manually_integrated_column, type_a_minimum_z, type_a_native_extent
  real(dp) :: phi0_v, phi_m_v, density_m3, ab_limit_coefficient
  integer(i32) :: status
  character(len=1) :: branch
  character(len=256) :: message
  character(len=2048) :: diagnostic_lines(5)
  integer :: fraction_index, atlas_index, diagnostic_unit, diagnostic_ios

  call test_init(29)

  call test_begin('stationary Zhao-A root is recovered by its interface field')
  call configure_params(60.0_dp, params)
  call solve_zhao_unknowns('zhao_a', params, phi0_v, phi_m_v, density_m3, branch)
  call assert_true(branch == 'A', 'stationary Zhao-A reference selected the wrong branch')
  call assert_close_dp(phi0_v, 2.9712182827319435_dp, 5.0e-6_dp, &
                       'stationary Zhao-A reference phi0 changed')
  call assert_close_dp(phi_m_v, -0.8169121871620854_dp, 5.0e-6_dp, &
                       'stationary Zhao-A reference phi_m changed')
  call assert_close_dp(density_m3, 7.819215729579456e6_dp, 5.0e1_dp, &
                       'stationary Zhao-A reference density changed')
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

  call test_begin('continuation diagnostics use an unambiguous full-range schema')
  continuation_diagnostics = zhao_continuation_diagnostics_type()
  continuation_diagnostics%from_density_m3 = huge(1.0_dp)
  continuation_diagnostics%candidate_density_m3 = tiny(1.0_dp)
  open (newunit=diagnostic_unit, status='scratch', action='readwrite', iostat=diagnostic_ios)
  call assert_true(diagnostic_ios == 0, 'failed to open scratch stream for Zhao diagnostics')
  if (diagnostic_ios == 0) then
    call write_zhao_continuation_diagnostics( &
      diagnostic_unit, continuation_diagnostics, 'test_stage', 42_i32 &
      )
    rewind (diagnostic_unit)
    diagnostic_lines = ''
    do atlas_index = 1, size(diagnostic_lines)
      read (diagnostic_unit, '(A)', iostat=diagnostic_ios) diagnostic_lines(atlas_index)
      call assert_true(diagnostic_ios == 0, 'Zhao diagnostics did not emit exactly five readable records')
    end do
    close (diagnostic_unit)
    call assert_true( &
      index(diagnostic_lines(1), 'call_stage=test_stage') > 0 .and. &
      index(diagnostic_lines(1), 'batch=42') > 0, &
      'Zhao diagnostics lost call-stage or batch context' &
      )
    call assert_true( &
      index(diagnostic_lines(3), 'from_branch=-') > 0 .and. &
      index(diagnostic_lines(4), 'candidate_branch=-') > 0, &
      'Zhao diagnostics did not encode missing branches with a sentinel' &
      )
    do atlas_index = 1, size(diagnostic_lines)
      call assert_true( &
        index(diagnostic_lines(atlas_index), '*') == 0, &
        'Zhao diagnostics overflowed its full-range real format' &
        )
    end do
  end if
  call test_end()

  call test_begin('no-photo Zhao-C stationary root is recovered without photoelectron inputs')
  call build_zhao_params( &
    60.0_dp, 8.7e6_dp, 8.7e6_dp, 12.0_dp, 12.0_dp, 4.0529988897111727e5_dp, &
    4.0529988897111727e5_dp, proton_mass, electron_mass, params, &
    photoelectron_population_fraction=0.0_dp, photoelectron_source_scale=0.0_dp &
    )
  call assert_close_dp(params%n_phe0_m3, 0.0_dp, 0.0_dp, 'no-photo Zhao source density must vanish')
  call solve_zhao_unknowns('zhao_c', params, phi0_v, phi_m_v, density_m3, branch)
  call assert_true(branch == 'C', 'stationary no-photo Zhao root must select Type C')
  stationary = zhao_charge_root_type( &
               branch=branch, phi0_v=phi0_v, phi_m_v=phi_m_v, n_swe_inf_m3=density_m3 &
               )
  call evaluate_zhao_interface_field(params, stationary, equilibrium_field, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'no-photo Zhao-C field evaluation failed: '//trim(message))
  call assert_true(equilibrium_field < 0.0_dp, 'stationary no-photo Zhao-C field must be negative')
  call solve_zhao_charge_root( &
    'zhao_c', params, equilibrium_field, charge_root, status, message, initial_root=stationary &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'charge-driven no-photo Zhao-C solve failed: '//trim(message))
  call assert_close_dp(charge_root%phi0_v, phi0_v, 2.0e-5_dp, 'no-photo Zhao-C phi0 recovery mismatch')
  call assert_close_dp(charge_root%n_swe_inf_m3, density_m3, 2.0e2_dp, &
                       'no-photo Zhao-C density recovery mismatch')
  call assert_true(abs(charge_root%net_current_density_a_m2) < 1.0e-11_dp, &
                   'stationary no-photo Zhao-C net current is not zero')
  call test_end()

  call test_begin('field-column homotopy rejects a degenerate no-photo Zhao-C column')
  call trace_zhao_field_column_homotopy( &
    params, equilibrium_field, 0.0_dp, equilibrium_field, 0.0_dp, &
    profile_points, control_length_m, charge_root, field_column_homotopy, status, message &
    )
  call assert_equal_i32( &
    status, outer_plasma_invalid, 'no-photo Zhao-C homotopy did not reject its zero column' &
    )
  call assert_true( &
    trim(field_column_homotopy%termination_reason) == 'degenerate_zero_photoelectron_column', &
    'no-photo Zhao-C homotopy lost its degenerate-column classification' &
    )
  call test_end()

  call test_begin('field perturbation leaves a nonzero Zhao-A charging current')
  call configure_params(60.0_dp, params)
  call solve_zhao_unknowns('zhao_a', params, phi0_v, phi_m_v, density_m3, branch)
  stationary = zhao_charge_root_type( &
               branch=branch, phi0_v=phi0_v, phi_m_v=phi_m_v, n_swe_inf_m3=density_m3 &
               )
  call evaluate_zhao_interface_field(params, stationary, equilibrium_field, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'stationary Zhao-A field evaluation failed before perturbation')
  call solve_zhao_charge_root( &
    'zhao_a', params, equilibrium_field, charge_root, status, message, initial_root=stationary &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'stationary Zhao-A recovery failed before perturbation')
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
  field_target = 1.01_dp*equilibrium_field
  call solve_zhao_charge_root( &
    'zhao_c', params, field_target, field_target_root, status, message, initial_root=charge_root &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'nearby Zhao-C target solve failed: '//trim(message))
  call build_zhao_outer_profile( &
    params, field_target_root, profile_points, field_target_state, status, message, &
    control_length_m=state%z(state%profile_n) &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'nearby Zhao-C target profile failed: '//trim(message))
  field_target_column = field_target_state%photoelectron_column_per_area
  call assert_true( &
    abs(field_target - equilibrium_field) > 1.0e-3_dp*abs(equilibrium_field) .and. &
    abs(field_target_column - state%photoelectron_column_per_area) > &
    1.0e-4_dp*max(field_target_column, state%photoelectron_column_per_area), &
    'nearby Zhao-C fixture did not change both prescribed quantities' &
    )
  homotopy_options = zhao_field_column_homotopy_options_type()
  homotopy_options%max_points = 96_i32
  homotopy_options%homotopy_max = 1.0_dp
  call trace_zhao_field_column_homotopy( &
    params, equilibrium_field, state%photoelectron_column_per_area, &
    field_target, field_target_column, profile_points, state%z(state%profile_n), &
    charge_root, field_column_homotopy, status, message, options=homotopy_options &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'populated Zhao-C homotopy failed: '//trim(message))
  call assert_true( &
    field_column_homotopy%target_reached .and. field_column_homotopy%branch == 'C' .and. &
    field_column_homotopy%point_count > 2_i32 .and. &
    .not. field_column_homotopy%homotopy_fold_detected .and. &
    trim(field_column_homotopy%termination_reason) == 'target_reached', &
    'populated Zhao-C homotopy did not reach its rank-regular target' &
    )
  call assert_close_dp( &
    field_column_homotopy%points(field_column_homotopy%point_count)%homotopy_fraction, &
    1.0_dp, 1.0e-14_dp, 'populated Zhao-C homotopy did not land exactly at lambda=1' &
    )
  call assert_close_dp( &
    field_column_homotopy%points(field_column_homotopy%point_count)%root%interface_field_v_m, &
    field_target, 1.0e-12_dp, 'populated Zhao-C homotopy missed its target field' &
    )
  call assert_close_dp( &
    field_column_homotopy%points(field_column_homotopy%point_count)%root% &
    photoelectron_column_per_area, field_target_column, &
    1.0e-8_dp*field_target_column, 'populated Zhao-C homotopy missed its target column' &
    )
  call assert_true( &
    abs(field_column_homotopy%points(field_column_homotopy%point_count)%root%phi0_v - &
        field_target_root%phi0_v)/params%t_phe_ev < 5.0e-8_dp .and. &
    abs(log(field_column_homotopy%points(field_column_homotopy%point_count)%root% &
            n_swe_inf_m3/field_target_root%n_swe_inf_m3)) < 1.0e-8_dp .and. &
    abs(field_column_homotopy%points(field_column_homotopy%point_count)%root% &
        photoelectron_population_fraction - &
        field_target_root%photoelectron_population_fraction) < 1.0e-8_dp, &
    'populated Zhao-C homotopy did not recover its independent target root' &
    )
  do atlas_index = 1, int(field_column_homotopy%point_count)
    call assert_true( &
      field_column_homotopy%points(atlas_index)%root%branch == 'C' .and. &
      field_column_homotopy%points(atlas_index)%root%residual_norm <= &
      5.0_dp*homotopy_options%residual_tolerance .and. &
      abs(field_column_homotopy%points(atlas_index)%normalized_column_residual) <= &
      5.0_dp*homotopy_options%residual_tolerance .and. &
      ieee_is_finite(field_column_homotopy%points(atlas_index)%row_rank_indicator) .and. &
      field_column_homotopy%points(atlas_index)%row_rank_indicator > 1.0e-12_dp, &
      'populated Zhao-C homotopy accepted an unconverged or rank-deficient point' &
      )
  end do
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

  call test_begin('field-column homotopy lands on a smooth Zhao-B target')
  homotopy_options = zhao_field_column_homotopy_options_type()
  homotopy_options%max_points = 96_i32
  homotopy_options%homotopy_max = 1.0_dp
  call trace_zhao_field_column_homotopy( &
    params, population_gap_field, reference_column, field_target, field_target_column, &
    profile_points, control_length_m, reference_root, field_column_homotopy, status, message, &
    options=homotopy_options &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'smooth Zhao-B field-column homotopy failed: '//trim(message))
  call assert_true(field_column_homotopy%target_reached, 'smooth Zhao-B homotopy did not reach lambda=1')
  call assert_true( &
    .not. field_column_homotopy%homotopy_fold_detected .and. &
    field_column_homotopy%point_count > 2_i32 .and. &
    trim(field_column_homotopy%termination_reason) == 'target_reached', &
    'smooth Zhao-B homotopy introduced an unexpected lambda fold' &
    )
  call assert_close_dp( &
    field_column_homotopy%points(field_column_homotopy%point_count)%homotopy_fraction, &
    1.0_dp, 1.0e-14_dp, 'smooth Zhao-B homotopy did not land exactly at lambda=1' &
    )
  call assert_close_dp( &
    field_column_homotopy%points(field_column_homotopy%point_count)%root% &
    photoelectron_column_per_area, field_target_column, &
    1.0e-8_dp*field_target_column, 'smooth Zhao-B homotopy missed its target column' &
    )
  call assert_true( &
    abs(field_column_homotopy%points(field_column_homotopy%point_count)%root%phi0_v - &
        field_target_root%phi0_v)/params%t_phe_ev < 1.0e-8_dp .and. &
    abs(log(field_column_homotopy%points(field_column_homotopy%point_count)%root% &
            n_swe_inf_m3/field_target_root%n_swe_inf_m3)) < 1.0e-8_dp .and. &
    abs(field_column_homotopy%points(field_column_homotopy%point_count)%root% &
        photoelectron_population_fraction - &
        field_target_root%photoelectron_population_fraction) < 1.0e-8_dp, &
    'smooth Zhao-B homotopy did not recover its independent target root' &
    )
  call assert_true( &
    all(field_column_homotopy%points%root%residual_norm <= &
        5.0_dp*homotopy_options%residual_tolerance) .and. &
    all(abs(field_column_homotopy%points%normalized_column_residual) <= &
        5.0_dp*homotopy_options%residual_tolerance) .and. &
    all(ieee_is_finite(field_column_homotopy%points%row_rank_indicator)) .and. &
    all(field_column_homotopy%points%row_rank_indicator > 1.0e-12_dp), &
    'smooth Zhao-B homotopy accepted an unconverged or non-finite point' &
    )
  call test_end()

  call test_begin('field-column homotopy reports a bounded point search')
  homotopy_options = zhao_field_column_homotopy_options_type()
  homotopy_options%max_points = 2_i32
  call trace_zhao_field_column_homotopy( &
    params, population_gap_field, reference_column, field_target, field_target_column, &
    profile_points, control_length_m, reference_root, field_column_homotopy, status, message, &
    options=homotopy_options &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'bounded Zhao-B homotopy failed numerically')
  call assert_true( &
    .not. field_column_homotopy%target_reached .and. &
    field_column_homotopy%point_count == 2_i32 .and. &
    trim(field_column_homotopy%termination_reason) == 'point_limit', &
    'bounded Zhao-B homotopy lost its point-limit classification' &
    )
  homotopy_options = zhao_field_column_homotopy_options_type()
  homotopy_options%log_density_floor = -30.0_dp
  call trace_zhao_field_column_homotopy( &
    params, population_gap_field, reference_column, field_target, field_target_column, &
    profile_points, control_length_m, reference_root, field_column_homotopy, status, message, &
    options=homotopy_options &
    )
  call assert_equal_i32(status, outer_plasma_invalid, 'decoder-edge density floor was not rejected')
  call assert_true( &
    trim(field_column_homotopy%termination_reason) == 'invalid_options', &
    'decoder-edge density floor lost its invalid-option classification' &
    )
  atlas_options = zhao_branch_atlas_options_type()
  atlas_options%log_density_floor = -30.0_dp
  call trace_zhao_branch_atlas( &
    params, population_gap_field, profile_points, control_length_m, reference_root, &
    branch_atlas, status, message, options=atlas_options &
    )
  call assert_equal_i32(status, outer_plasma_invalid, 'atlas decoder-edge density floor was not rejected')
  call assert_true( &
    trim(branch_atlas%termination_reason) == 'invalid_options', &
    'atlas decoder-edge density floor lost its invalid-option classification' &
    )
  call test_end()

  call test_begin('field-column homotopy honors a requested eta lower bound')
  homotopy_options = zhao_field_column_homotopy_options_type()
  homotopy_options%max_points = 96_i32
  homotopy_options%eta_min = 7.5e-2_dp
  call trace_zhao_field_column_homotopy( &
    params, population_gap_field, reference_column, population_gap_field, &
    continuation_target_column, profile_points, control_length_m, reference_root, &
    field_column_homotopy, status, message, options=homotopy_options &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'eta-bounded Zhao-B homotopy failed numerically')
  call assert_true( &
    .not. field_column_homotopy%target_reached .and. &
    trim(field_column_homotopy%termination_reason) == 'eta_lower_search_limit' .and. &
    minval(field_column_homotopy%points%root%photoelectron_population_fraction) >= &
    homotopy_options%eta_min, &
    'eta-bounded Zhao-B homotopy crossed its requested lower bound' &
    )
  call test_end()

  call test_begin('column continuation classifies a non-finite root numerically')
  nonfinite_root = reference_root
  nonfinite_root%phi0_v = ieee_value(0.0_dp, ieee_quiet_nan)
  call solve_outer_plasma_zhao_column( &
    'auto', params, population_gap_field, profile_points, control_length_m, reference_column, &
    state, charge_root, status, message, initial_root=nonfinite_root, &
    diagnostics=continuation_diagnostics &
    )
  call assert_equal_i32( &
    status, outer_plasma_numerical_failure, &
    'non-finite continuation root was not classified numerically: '//trim(message) &
    )
  call assert_equal_i32( &
    continuation_diagnostics%reason_code, zhao_continuation_reason_numerical_failure, &
    'non-finite continuation root lost its numerical-failure reason' &
    )
  call assert_true( &
    trim(continuation_diagnostics%reason) == 'nonfinite_state', &
    'non-finite continuation root recorded the wrong reason detail' &
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
  call build_zhao_params( &
    60.0_dp, 8.7e6_dp, 8.7e6_dp, 12.0_dp, 12.0_dp, 4.0529988897111727e5_dp, &
    4.0529988897111727e5_dp, proton_mass, electron_mass, params, &
    photoelectron_population_fraction=0.0_dp, photoelectron_source_scale=0.0_dp &
    )
  call solve_zhao_unknowns('zhao_c', params, phi0_v, phi_m_v, density_m3, branch)
  stationary = zhao_charge_root_type( &
               branch=branch, phi0_v=phi0_v, phi_m_v=phi_m_v, n_swe_inf_m3=density_m3 &
               )
  call evaluate_zhao_interface_field(params, stationary, equilibrium_field, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'no-photo eta-limit field evaluation failed')
  continuation_diagnostics = zhao_continuation_diagnostics_type()
  call solve_outer_plasma_zhao_column( &
    'auto', params, equilibrium_field, profile_points, control_length_m, 1.0_dp, &
    state, charge_root, status, message, diagnostics=continuation_diagnostics &
    )
  call assert_equal_i32( &
    status, outer_plasma_no_physical_solution, &
    'no-photo positive column must stop at the finite eta search limit' &
    )
  call assert_equal_i32( &
    continuation_diagnostics%reason_code, zhao_continuation_reason_search_limit, &
    'finite eta search bound was misclassified as target-unreachable' &
    )
  call assert_true( &
    trim(continuation_diagnostics%reason) == 'eta_upper_search_limit', &
    'finite eta search bound recorded the wrong reason detail' &
    )
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

  call test_begin('strong-UV column failure records the near-zero-density Zhao-B jump')
  call build_zhao_params( &
    60.0_dp, 5.0e6_dp, 130.68910376476308e6_dp, 10.0_dp, 2.2_dp, &
    4.0e5_dp, 4.0e5_dp, 1.672482821616e-27_dp, 9.10938356e-31_dp, params, &
    photoelectron_population_fraction=1.0_dp, photoelectron_source_scale=1.0_dp &
    )
  reference_root = zhao_charge_root_type( &
                   branch='B', interface_side='monotonic', &
                   phi0_v=2.2999370729879152_dp, phi_m_v=2.2999370729879152_dp, &
                   n_swe_inf_m3=2.2336848555940270e-1_dp, &
                   photoelectron_population_fraction=2.5133484117801319e-1_dp, &
                   photoelectron_column_per_area=9.7614063582301259e7_dp, &
                   photoelectron_column_target_per_area=9.9455765203101724e7_dp, &
                   photoelectron_column_residual_per_area=-1.8417016208004653e6_dp, &
                   interface_field_v_m=9.0729627587860184e-1_dp, &
                   residual_norm=9.9958608057448828e-10_dp &
                   )
  call solve_outer_plasma_zhao_column( &
    'auto', params, 9.0729627587860184e-1_dp, 128_i32, &
    10.0_dp*params%lambda_d_phe_ref_m, 9.9455765203101724e7_dp, &
    state, charge_root, status, message, initial_root=reference_root, &
    diagnostics=continuation_diagnostics &
    )
  call assert_equal_i32( &
    status, outer_plasma_numerical_failure, &
    'strong-UV characterization changed its aggregated return status: '//trim(message) &
    )
  call assert_equal_i32( &
    continuation_diagnostics%reason_code, zhao_continuation_reason_guard_rejected, &
    'strong-UV characterization lost the same-branch guard decision' &
    )
  call assert_true( &
    trim(continuation_diagnostics%reason) == 'same_branch_jump_guard', &
    'strong-UV characterization recorded the wrong reason detail' &
    )
  call assert_true( &
    continuation_diagnostics%candidate_available .and. &
    continuation_diagnostics%from_branch == 'B' .and. &
    continuation_diagnostics%candidate_branch == 'B', &
    'strong-UV characterization did not preserve the Zhao-B root pair' &
    )
  call assert_close_dp( &
    continuation_diagnostics%from_density_m3, 2.2336848555940270e-1_dp, 1.0e-10_dp, &
    'strong-UV characterization changed the accepted ambient density' &
    )
  call assert_close_dp( &
    continuation_diagnostics%candidate_density_m3, 1.6511518550182400e-1_dp, 5.0e-4_dp, &
    'strong-UV characterization changed the rejected ambient density' &
    )
  call assert_true( &
    continuation_diagnostics%normalized_potential_jump < 1.0e-7_dp .and. &
    continuation_diagnostics%log_density_jump > 0.25_dp, &
    'strong-UV root jump is no longer isolated to the near-zero density coordinate' &
    )
  call assert_true( &
    reference_root%n_swe_inf_m3/params%n_phe_ref_m3 < 2.0e-9_dp, &
    'strong-UV predecessor is no longer at the ambient-density positivity boundary' &
    )
  atlas_options = zhao_branch_atlas_options_type()
  atlas_options%eta_direction = 1_i32
  call trace_zhao_branch_atlas( &
    params, reference_root%interface_field_v_m, 128_i32, &
    10.0_dp*params%lambda_d_phe_ref_m, reference_root, branch_atlas, status, message, &
    target_column_m2=9.9455765203101724e7_dp, options=atlas_options &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'strong-UV branch atlas failed: '//trim(message))
  call assert_equal_i32( &
    branch_atlas%termination_reason_code, zhao_continuation_reason_search_limit, &
    'strong-UV forward Zhao-B atlas did not report its finite density floor' &
    )
  call assert_true( &
    trim(branch_atlas%termination_reason) == 'ambient_density_floor_limit', &
    'strong-UV forward Zhao-B atlas reported the wrong search limit' &
    )
  call assert_true(branch_atlas%point_count > 20_i32, 'strong-UV branch atlas stopped prematurely')
  call assert_true(.not. branch_atlas%target_bracketed, 'strong-UV atlas falsely bracketed the target')
  call assert_true( &
    branch_atlas%maximum_column_m2 < 9.9455765203101724e7_dp - 1.8e6_dp, &
    'strong-UV forward Zhao-B tail approached the target unexpectedly' &
    )
  call assert_true( &
    branch_atlas%points(branch_atlas%point_count)%root%n_swe_inf_m3/params%n_phe_ref_m3 < 5.0e-12_dp, &
    'strong-UV atlas did not approach its configured ambient-density floor' &
    )
  call assert_true( &
    .not. branch_atlas%eta_fold_detected .and. .not. branch_atlas%column_fold_detected, &
    'strong-UV forward Zhao-B tail was incorrectly classified as a fold' &
    )
  call assert_true( &
    branch_atlas%seed_reanchored .and. &
    branch_atlas%seed_refinement_jump > 0.25_dp .and. &
    branch_atlas%seed_refinement_jump <= atlas_options%seed_refinement_limit, &
    'strong-UV atlas did not disclose refinement beyond the production same-root guard' &
    )
  do atlas_index = 1, int(branch_atlas%point_count)
    call assert_true( &
      branch_atlas%points(atlas_index)%root%branch == 'B' .and. &
      branch_atlas%points(atlas_index)%root%residual_norm <= 5.0e-12_dp, &
      'strong-UV atlas left its fixed converged Zhao-B curve' &
      )
  end do
  atlas_options%log_density_floor = -24.0_dp
  call trace_zhao_branch_atlas( &
    params, reference_root%interface_field_v_m, 128_i32, &
    10.0_dp*params%lambda_d_phe_ref_m, reference_root, coarse_branch_atlas, status, message, &
    target_column_m2=9.9455765203101724e7_dp, options=atlas_options &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'coarse-floor strong-UV atlas failed: '//trim(message))
  call assert_equal_i32( &
    coarse_branch_atlas%termination_reason_code, zhao_continuation_reason_search_limit, &
    'coarse-floor strong-UV atlas did not stop at its density floor' &
    )
  call assert_true( &
    trim(coarse_branch_atlas%termination_reason) == 'ambient_density_floor_limit', &
    'coarse-floor strong-UV atlas stopped at a different search limit' &
    )
  call assert_true( &
    coarse_branch_atlas%points(coarse_branch_atlas%point_count)%root%n_swe_inf_m3/ &
    params%n_phe_ref_m3 < 5.0e-11_dp, &
    'coarse-floor strong-UV atlas did not approach its configured density floor' &
    )
  call assert_true( &
    abs(branch_atlas%maximum_column_m2 - coarse_branch_atlas%maximum_column_m2) < 10.0_dp, &
    'strong-UV forward Zhao-B maximum column is not converged with density floor' &
    )
  call assert_true( &
    .not. coarse_branch_atlas%target_bracketed, &
    'coarse-floor strong-UV atlas falsely bracketed the target' &
    )
  call test_end()

  call test_begin('strong-UV Zhao-B endpoint has no regular local Type-A tangent')
  call diagnose_zhao_ab_degeneracy( &
    params, branch_atlas%points(branch_atlas%point_count)%root, &
    ab_diagnostics, status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'strong-UV A/B degeneracy diagnostic failed: '//trim(message))
  call assert_true( &
    trim(ab_diagnostics%classification) == 'no_regular_type_a_tangent', &
    'strong-UV A/B endpoint received the wrong local classification' &
    )
  call assert_true( &
    ab_diagnostics%density_zero_limit .and. &
    ab_diagnostics%b_root_field_condition_met .and. &
    ab_diagnostics%limiting_field_condition_met, &
    'strong-UV B endpoint did not approach its density-zero field limit' &
    )
  call assert_true( &
    .not. ab_diagnostics%far_field_tangent_condition_met .and. &
    .not. ab_diagnostics%regular_connection_conditions_met, &
    'strong-UV B endpoint falsely admitted a regular local Type-A connection' &
    )
  call assert_true( &
    ab_diagnostics%ambient_density_ratio < 5.0e-12_dp .and. &
    abs(ab_diagnostics%doubled_quasineutral_residual_hat) < 1.0e-10_dp .and. &
    abs(ab_diagnostics%b_field_squared_residual_hat) < 1.0e-10_dp, &
    'strong-UV B endpoint is no longer at its quasineutral density boundary' &
    )
  call assert_true( &
    ab_diagnostics%limiting_field_squared_jump_hat <= 0.0_dp .and. &
    ab_diagnostics%quasineutral_far_field_q3_coefficient > 2.8e-2_dp .and. &
    ab_diagnostics%quasineutral_far_field_q3_coefficient < 3.0e-2_dp, &
    'strong-UV Type-A far-field q^3 limit changed unexpectedly' &
    )
  call assert_true( &
    ab_diagnostics%probe_available .and. &
    ab_diagnostics%probe_ambient_density_ratio > 0.0_dp .and. &
    abs(ab_diagnostics%probe_quasineutral_residual) < 1.0e-12_dp .and. &
    abs(ab_diagnostics%probe_far_field_q3_coefficient - &
        ab_diagnostics%quasineutral_far_field_q3_coefficient) < 5.0e-4_dp, &
    'finite-q Type-A probe did not converge to its analytic far-field limit' &
    )
  call assert_true( &
    abs(ab_diagnostics%probe_field_squared_residual_hat) > 5.0e-5_dp, &
    'finite-q Type-A probe unexpectedly solved the fixed-field equation' &
    )
  ab_limit_coefficient = ab_diagnostics%quasineutral_far_field_q3_coefficient
  zero_root = branch_atlas%points(branch_atlas%point_count)%root
  zero_root%n_swe_inf_m3 = 0.0_dp
  zero_root%photoelectron_population_fraction = &
    2.0_dp*params%n_swi_inf_m3*exp(zero_root%phi0_v/params%t_phe_ev)/params%n_phe0_m3
  call diagnose_zhao_ab_degeneracy( &
    params, zero_root, ab_diagnostics, status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'exact density-zero A/B diagnostic failed: '//trim(message))
  call assert_true( &
    trim(ab_diagnostics%classification) == 'no_regular_type_a_tangent' .and. &
    ab_diagnostics%ambient_density_ratio == 0.0_dp .and. &
    abs(ab_diagnostics%quasineutral_far_field_q3_coefficient - ab_limit_coefficient) < 1.0e-6_dp, &
    'finite-density atlas endpoint did not converge to its analytic density-zero limit' &
    )
  call test_end()

  call test_begin('strong-UV reverse Zhao-B atlas does not reach the larger target column')
  atlas_options = zhao_branch_atlas_options_type()
  atlas_options%eta_direction = -1_i32
  atlas_options%max_points = 512_i32
  call trace_zhao_branch_atlas( &
    params, reference_root%interface_field_v_m, 128_i32, &
    10.0_dp*params%lambda_d_phe_ref_m, branch_atlas%points(1)%root, &
    reverse_branch_atlas, status, message, &
    target_column_m2=9.9455765203101724e7_dp, options=atlas_options &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'reverse strong-UV Zhao-B atlas failed: '//trim(message))
  call assert_true( &
    .not. reverse_branch_atlas%target_bracketed .and. &
    reverse_branch_atlas%maximum_column_m2 < 9.9455765203101724e7_dp - 1.8e6_dp, &
    'reverse strong-UV Zhao-B atlas unexpectedly reached the target column' &
    )
  call assert_true( &
    .not. reverse_branch_atlas%eta_fold_detected .and. &
    .not. reverse_branch_atlas%column_fold_detected, &
    'reverse strong-UV Zhao-B atlas introduced an unexpected fold' &
    )
  call assert_true( &
    reverse_branch_atlas%point_count > 1_i32 .and. &
    reverse_branch_atlas%points(reverse_branch_atlas%point_count)%root% &
    photoelectron_population_fraction < &
    reverse_branch_atlas%points(1)%root%photoelectron_population_fraction, &
    'reverse strong-UV Zhao-B atlas did not advance toward smaller eta' &
    )
  call assert_true( &
    reverse_branch_atlas%termination_reason_code == zhao_continuation_reason_search_limit .and. &
    trim(reverse_branch_atlas%termination_reason) == 'eta_lower_search_limit', &
    'reverse strong-UV Zhao-B atlas stopped before its eta lower limit' &
    )
  do atlas_index = 1, int(reverse_branch_atlas%point_count)
    call assert_true( &
      reverse_branch_atlas%points(atlas_index)%root%branch == 'B' .and. &
      reverse_branch_atlas%points(atlas_index)%root%residual_norm <= 5.0e-12_dp, &
      'reverse strong-UV atlas left its converged Zhao-B curve' &
      )
  end do
  call test_end()

  call test_begin('strong-UV density candidates refine onto the same Zhao-B root')
  continuation_target_root = zhao_charge_root_type( &
                             branch='B', interface_side='monotonic', &
                             phi0_v=2.2999370655044924_dp, phi_m_v=2.2999370655044924_dp, &
                             n_swe_inf_m3=1.6511518550182400e-1_dp, &
                             photoelectron_population_fraction=reference_root%photoelectron_population_fraction, &
                             interface_field_v_m=reference_root%interface_field_v_m &
                             )
  atlas_options = zhao_branch_atlas_options_type()
  atlas_options%max_points = 1_i32
  call trace_zhao_branch_atlas( &
    params, reference_root%interface_field_v_m, 128_i32, &
    10.0_dp*params%lambda_d_phe_ref_m, continuation_target_root, &
    candidate_branch_atlas, status, message, options=atlas_options &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'candidate Zhao-B seed refinement failed: '//trim(message))
  call assert_true( &
    branch_atlas%seed_reanchored .and. branch_atlas%seed_refinement_jump > 0.25_dp .and. &
    .not. candidate_branch_atlas%seed_reanchored .and. &
    candidate_branch_atlas%seed_refinement_jump < 5.0e-2_dp, &
    'strong-UV runtime and candidate seeds recorded the wrong refinement neighborhoods' &
    )
  call assert_true( &
    abs(candidate_branch_atlas%points(1)%root%phi0_v - &
        branch_atlas%points(1)%root%phi0_v)/params%t_phe_ev < 1.0e-8_dp .and. &
    abs(log(candidate_branch_atlas%points(1)%root%n_swe_inf_m3/ &
            branch_atlas%points(1)%root%n_swe_inf_m3)) < 1.0e-5_dp, &
    'strong-UV density candidates did not refine onto the same local Zhao-B root' &
    )
  call test_end()

  call test_begin('strong-UV batch transition traces the coupled Zhao-B manifold')
  ! Recovered from project_dust_release run R20260722-0005, whose child step
  ! 8184334.0 failed at batch 16,
  ! by the 15-batch replay job 8186918 using BEACH commit 5ccc7f9.  ADR 0005
  ! records the executable hash and the limitation that replay raw files were
  ! not retained.
  reference_root = zhao_charge_root_type( &
                   branch='B', interface_side='monotonic', &
                   phi0_v=2.1622086033018255_dp, phi_m_v=2.1622086033018255_dp, &
                   n_swe_inf_m3=1.3842809338228058e5_dp, &
                   photoelectron_population_fraction=2.3203997288209771e-1_dp, &
                   photoelectron_column_per_area=9.3202065343372747e7_dp, &
                   photoelectron_column_target_per_area=9.3202065681229725e7_dp, &
                   photoelectron_column_residual_per_area=-3.3785697817802429e-1_dp, &
                   interface_field_v_m=8.4245706656096919e-1_dp, &
                   residual_norm=1.4155343563970746e-15_dp &
                   )
  homotopy_options = zhao_field_column_homotopy_options_type()
  homotopy_options%log_density_floor = -27.0_dp
  call trace_zhao_field_column_homotopy( &
    params, reference_root%interface_field_v_m, &
    reference_root%photoelectron_column_target_per_area, 9.0729627587860184e-1_dp, &
    9.9455765203101724e7_dp, 128_i32, 10.0_dp*params%lambda_d_phe_ref_m, &
    reference_root, field_column_homotopy, status, message, options=homotopy_options &
    )
  call assert_equal_i32( &
    status, outer_plasma_ok, 'strong-UV coupled Zhao-B homotopy failed numerically: '//trim(message) &
    )
  call assert_true( &
    field_column_homotopy%point_count > 1_i32, &
    'strong-UV coupled Zhao-B homotopy did not leave its batch-15 seed' &
    )
  call assert_true( &
    .not. field_column_homotopy%target_reached .and. &
    .not. field_column_homotopy%homotopy_fold_detected, &
    'strong-UV coupled Zhao-B path unexpectedly reached or folded before its endpoint' &
    )
  call assert_true( &
    field_column_homotopy%termination_reason_code == zhao_continuation_reason_search_limit .and. &
    trim(field_column_homotopy%termination_reason) == 'ambient_density_floor_limit', &
    'strong-UV coupled Zhao-B path did not terminate at its finite density floor' &
    )
  call assert_true( &
    abs(field_column_homotopy%points(field_column_homotopy%point_count)%homotopy_fraction - &
        0.3317869678461613_dp) < 1.0e-5_dp .and. &
    field_column_homotopy%points(field_column_homotopy%point_count)%root%n_swe_inf_m3/ &
    params%n_phe_ref_m3 < 5.0e-12_dp, &
    'strong-UV coupled Zhao-B endpoint moved away from the one-third homotopy density limit' &
    )
  call assert_true( &
    log(field_column_homotopy%points(field_column_homotopy%point_count)%root%n_swe_inf_m3/ &
        params%n_phe_ref_m3) <= homotopy_options%log_density_floor .and. &
    log(field_column_homotopy%points(field_column_homotopy%point_count)%root%n_swe_inf_m3/ &
        params%n_phe_ref_m3) >= homotopy_options%log_density_floor - &
    field_column_homotopy%points(field_column_homotopy%point_count)%accepted_step - 1.0e-10_dp .and. &
    log(field_column_homotopy%points(field_column_homotopy%point_count - 1_i32)%root% &
        n_swe_inf_m3/params%n_phe_ref_m3) > homotopy_options%log_density_floor .and. &
    .not. field_column_homotopy%seed_reanchored .and. &
    field_column_homotopy%seed_refinement_jump < 1.0e-6_dp, &
    'strong-UV coupled Zhao-B trace did not stop at its first density-floor crossing' &
    )
  call assert_true( &
    abs(field_column_homotopy%points(field_column_homotopy%point_count)% &
        normalized_column_residual) < 1.0e-8_dp, &
    'strong-UV coupled Zhao-B endpoint does not close its prescribed intermediate column' &
    )
  do atlas_index = 1, int(field_column_homotopy%point_count)
    call assert_true( &
      field_column_homotopy%points(atlas_index)%root%branch == 'B' .and. &
      field_column_homotopy%points(atlas_index)%root%residual_norm <= &
      5.0_dp*homotopy_options%residual_tolerance .and. &
      abs(field_column_homotopy%points(atlas_index)%normalized_column_residual) <= &
      5.0_dp*homotopy_options%residual_tolerance .and. &
      ieee_is_finite(field_column_homotopy%points(atlas_index)%row_rank_indicator) .and. &
      field_column_homotopy%points(atlas_index)%row_rank_indicator > 1.0e-12_dp, &
      'strong-UV coupled Zhao-B trace accepted an unconverged or rank-deficient point' &
      )
  end do
  call diagnose_zhao_ab_degeneracy( &
    params, field_column_homotopy%points(field_column_homotopy%point_count)%root, &
    ab_diagnostics, status, message &
    )
  call assert_equal_i32( &
    status, outer_plasma_ok, 'strong-UV coupled endpoint A/B diagnostic failed: '//trim(message) &
    )
  call assert_true( &
    trim(ab_diagnostics%classification) == 'no_regular_type_a_tangent' .and. &
    .not. ab_diagnostics%regular_connection_conditions_met .and. &
    ab_diagnostics%quasineutral_far_field_q3_coefficient > 2.8e-2_dp .and. &
    ab_diagnostics%quasineutral_far_field_q3_coefficient < 3.0e-2_dp, &
    'strong-UV coupled endpoint falsely admitted a regular local Type-A connection' &
    )
  homotopy_options%log_density_floor = -24.0_dp
  call trace_zhao_field_column_homotopy( &
    params, reference_root%interface_field_v_m, &
    reference_root%photoelectron_column_target_per_area, 9.0729627587860184e-1_dp, &
    9.9455765203101724e7_dp, 128_i32, 10.0_dp*params%lambda_d_phe_ref_m, &
    reference_root, coarse_field_column_homotopy, status, message, options=homotopy_options &
    )
  call assert_equal_i32( &
    status, outer_plasma_ok, 'coarse-floor coupled Zhao-B homotopy failed numerically: '//trim(message) &
    )
  call assert_true( &
    coarse_field_column_homotopy%termination_reason_code == &
    zhao_continuation_reason_search_limit .and. &
    trim(coarse_field_column_homotopy%termination_reason) == 'ambient_density_floor_limit' .and. &
    .not. coarse_field_column_homotopy%target_reached .and. &
    .not. coarse_field_column_homotopy%homotopy_fold_detected .and. &
    abs(coarse_field_column_homotopy%points(coarse_field_column_homotopy%point_count)% &
        homotopy_fraction - &
        field_column_homotopy%points(field_column_homotopy%point_count)%homotopy_fraction) < 1.0e-5_dp, &
    'strong-UV coupled endpoint depends materially on the diagnostic density floor' &
    )
  call assert_true( &
    log(coarse_field_column_homotopy%points(coarse_field_column_homotopy%point_count)%root% &
        n_swe_inf_m3/params%n_phe_ref_m3) <= homotopy_options%log_density_floor .and. &
    log(coarse_field_column_homotopy%points(coarse_field_column_homotopy%point_count)%root% &
        n_swe_inf_m3/params%n_phe_ref_m3) >= homotopy_options%log_density_floor - &
    coarse_field_column_homotopy%points(coarse_field_column_homotopy%point_count)%accepted_step - &
    1.0e-10_dp .and. &
    log(coarse_field_column_homotopy%points( &
        coarse_field_column_homotopy%point_count - 1_i32)%root%n_swe_inf_m3/ &
        params%n_phe_ref_m3) > homotopy_options%log_density_floor, &
    'coarse-floor coupled Zhao-B trace did not stop at its first density-floor crossing' &
    )
  do atlas_index = 1, int(coarse_field_column_homotopy%point_count)
    call assert_true( &
      coarse_field_column_homotopy%points(atlas_index)%root%branch == 'B' .and. &
      coarse_field_column_homotopy%points(atlas_index)%root%residual_norm <= &
      5.0_dp*homotopy_options%residual_tolerance .and. &
      abs(coarse_field_column_homotopy%points(atlas_index)%normalized_column_residual) <= &
      5.0_dp*homotopy_options%residual_tolerance .and. &
      ieee_is_finite(coarse_field_column_homotopy%points(atlas_index)%row_rank_indicator) .and. &
      coarse_field_column_homotopy%points(atlas_index)%row_rank_indicator > 1.0e-12_dp, &
      'coarse-floor coupled Zhao-B trace accepted an unconverged or rank-deficient point' &
      )
  end do
  call test_end()

  call test_begin('branch atlas brackets a reachable target on a smooth Zhao-B curve')
  call configure_params(60.0_dp, params)
  params%photoelectron_population_fraction = 0.10_dp
  call solve_zhao_charge_root( &
    'zhao_b', params, population_gap_field, reference_root, status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'smooth Zhao-B atlas seed failed: '//trim(message))
  params%photoelectron_population_fraction = 0.14_dp
  call solve_zhao_charge_root( &
    'zhao_b', params, population_gap_field, continuation_target_root, status, message, &
    initial_root=reference_root &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'smooth Zhao-B atlas target root failed: '//trim(message))
  call build_zhao_outer_profile( &
    params, continuation_target_root, profile_points, continuation_target_state, status, message, &
    control_length_m=control_length_m &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'smooth Zhao-B atlas target profile failed: '//trim(message))
  atlas_options = zhao_branch_atlas_options_type()
  atlas_options%max_points = 32_i32
  atlas_options%initial_step = 1.0e-2_dp
  atlas_options%maximum_step = 2.5e-2_dp
  call trace_zhao_branch_atlas( &
    params, population_gap_field, profile_points, control_length_m, reference_root, &
    branch_atlas, status, message, &
    target_column_m2=continuation_target_state%photoelectron_column_per_area, options=atlas_options &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'smooth Zhao-B branch atlas failed: '//trim(message))
  call assert_equal_i32( &
    branch_atlas%termination_reason_code, zhao_continuation_reason_target_bracketed, &
    'smooth Zhao-B branch atlas did not bracket its reachable target' &
    )
  call assert_true(branch_atlas%target_bracketed, 'smooth Zhao-B atlas lost its target crossing')
  call assert_true(branch_atlas%point_count > 1_i32, 'smooth Zhao-B atlas did not advance')
  call test_end()

  call test_begin('branch atlas exercises the four-coordinate Zhao-A corrector')
  call configure_params(60.0_dp, params)
  call solve_zhao_unknowns('zhao_a', params, phi0_v, phi_m_v, density_m3, branch)
  reference_root = zhao_charge_root_type( &
                   branch=branch, phi0_v=phi0_v, phi_m_v=phi_m_v, n_swe_inf_m3=density_m3, &
                   photoelectron_population_fraction=params%photoelectron_population_fraction &
                   )
  call evaluate_zhao_interface_field(params, reference_root, equilibrium_field, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'Zhao-A atlas seed field failed: '//trim(message))
  reference_root%interface_field_v_m = equilibrium_field
  atlas_options = zhao_branch_atlas_options_type()
  atlas_options%max_points = 2_i32
  atlas_options%eta_direction = -1_i32
  atlas_options%initial_step = 1.0e-3_dp
  atlas_options%maximum_step = 1.0e-3_dp
  call trace_zhao_branch_atlas( &
    params, equilibrium_field, profile_points, type_a_native_extent, reference_root, &
    branch_atlas, status, message, options=atlas_options &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'Zhao-A branch atlas failed: '//trim(message))
  call assert_equal_i32(branch_atlas%point_count, 2_i32, 'Zhao-A branch atlas did not accept one step')
  call assert_true( &
    branch_atlas%points(2)%root%branch == 'A' .and. &
    branch_atlas%points(2)%root%residual_norm <= 1.1e-12_dp, &
    'Zhao-A four-coordinate corrector left its converged branch' &
    )
  call test_end()

  call test_begin('field-column homotopy rejects unsupported Zhao-A coordinates')
  homotopy_options = zhao_field_column_homotopy_options_type()
  call trace_zhao_field_column_homotopy( &
    params, equilibrium_field, branch_atlas%points(1)%root%photoelectron_column_per_area, &
    equilibrium_field, branch_atlas%points(1)%root%photoelectron_column_per_area, &
    profile_points, type_a_native_extent, branch_atlas%points(1)%root, &
    field_column_homotopy, status, message, options=homotopy_options &
    )
  call assert_equal_i32(status, outer_plasma_invalid, 'Zhao-A field-column homotopy was not rejected')
  call assert_true( &
    trim(field_column_homotopy%termination_reason) == 'unsupported_type_a_homotopy', &
    'Zhao-A field-column homotopy lost its unsupported-coordinate classification' &
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

program test_dynamic_k0_zhao
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, pi, qe
  use bem_dynamic_k0_mean, only: dynamic_k0_ok, dynamic_k0_step_type
  use bem_dynamic_k0_zhao, only: &
    measured_interface_energy_distribution_type, &
    build_measured_interface_energy_distribution, advance_dynamic_k0_zhao, &
    zhao_profile_barrier_energy, measured_sample_escape_fraction, &
    dynamic_zhao_no_physical_root
  use bem_outer_plasma_kinetic, only: &
    kinetic_outer_plasma_options_type, solve_outer_plasma_zhao_stationary, &
    continue_outer_plasma_zhao_connected, continue_outer_plasma_zhao_connected_root, &
    materialize_outer_plasma_zhao_root, zhao_charge_root_type, &
    zhao_charge_root_barrier_energy, zhao_connected_path_certificate_type
  use bem_outer_plasma_types, only: &
    outer_plasma_state_type, outer_plasma_ok, outer_plasma_no_physical_solution
  use test_support, only: &
    test_init, test_begin, test_end, test_summary, assert_true, &
    assert_close_dp, assert_equal_i32
  implicit none

  real(dp), parameter :: electron_mass = 9.1093837015e-31_dp
  real(dp), parameter :: proton_mass = 1.67262192369e-27_dp
  real(dp), parameter :: area_xy_m2 = (9.89949e-5_dp)**2
  real(dp), parameter :: time_step_s = 1.0e-8_dp
  type(kinetic_outer_plasma_options_type) :: options
  type(kinetic_outer_plasma_options_type) :: closepack_options
  type(kinetic_outer_plasma_options_type) :: connected_options, type_c_options, extreme_options
  type(outer_plasma_state_type) :: stationary, accepted, closepack_stationary
  type(outer_plasma_state_type) :: connected_state, root_state, reverse_state, type_c_stationary
  type(zhao_charge_root_type) :: connected_root
  type(zhao_connected_path_certificate_type) :: path_certificate
  type(measured_interface_energy_distribution_type) :: distribution
  type(dynamic_k0_step_type) :: step
  real(dp) :: barrier_j, barrier_ev, reference_current_density
  real(dp) :: emitted_charge_c, anchor_charge_c, base_charge_c
  real(dp) :: effective_source_scale, high_fraction, marginal_fraction, tiny_group_charge_c
  real(dp) :: target_field, analytic_barrier_j
  real(dp) :: sample_energy_j(4), sample_charge_c(4)
  real(dp), allocatable :: many_energy_j(:), many_charge_c(:)
  integer(i32) :: status, connected_status, sample_index
  character(len=256) :: message

  call test_init(22)

  call test_begin('measured interface samples are stably grouped by descending energy')
  sample_energy_j = [2.0_dp, 4.0_dp, 2.0_dp, 1.0_dp]*qe
  sample_charge_c = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]*qe
  call build_measured_interface_energy_distribution( &
    sample_energy_j, sample_charge_c, distribution, status, message &
    )
  call assert_equal_i32(status, dynamic_k0_ok, 'measured energy grouping failed: '//trim(message))
  call assert_equal_i32(distribution%group_count, 3_i32, 'equal measured energies were not grouped')
  call assert_close_dp(distribution%energy_j(1)/qe, 4.0_dp, 1.0e-14_dp, &
                       'measured energies were not sorted descending')
  call assert_close_dp(distribution%energy_j(2)/qe, 2.0_dp, 1.0e-14_dp, &
                       'middle measured energy changed')
  call assert_close_dp(distribution%group_charge_c(2)/qe, 4.0_dp, 1.0e-14_dp, &
                       'equal-energy macro charges were not accumulated')
  call assert_close_dp(distribution%cumulative_charge_c(3)/qe, 10.0_dp, 1.0e-13_dp, &
                       'measured cumulative charge is not closed')
  call test_end()

  call configure_type_a_options(options)
  call solve_outer_plasma_zhao_stationary(options, stationary, status, message)

  call test_begin('stationary Zhao-A profile retains the analytic virtual cathode barrier')
  call assert_equal_i32(status, outer_plasma_ok, 'stationary Zhao-A solve failed: '//trim(message))
  call assert_true(stationary%ready .and. stationary%zhao_branch == 'A', &
                   'stationary strong-PE anchor did not select Zhao-A')
  call assert_close_dp(stationary%interface_potential, 2.9712182827319435_dp, 5.0e-6_dp, &
                       'stationary Zhao-A interface potential changed')
  call assert_close_dp(minval(stationary%potential), -0.8169121871620854_dp, 5.0e-6_dp, &
                       'stationary Zhao-A virtual-cathode minimum changed')
  barrier_j = zhao_profile_barrier_energy(stationary, -qe)
  barrier_ev = barrier_j/qe
  call assert_close_dp(barrier_ev, 3.788130469894029_dp, 1.0e-5_dp, &
                       'Zhao-A profile barrier changed')
  high_fraction = exp(-barrier_j/options%photoelectron_temperature_j)
  call assert_close_dp(high_fraction, 0.17873026907185158_dp, 2.0e-8_dp, &
                       'analytic Zhao-A Maxwellian escape fraction changed')
  call test_end()

  call test_begin('connected Zhao-A field path lands exactly and reverses on one sheet')
  connected_options = options
  target_field = 1.01_dp*stationary%interface_field
  connected_options%interface_field = target_field
  call continue_outer_plasma_zhao_connected( &
    connected_options, stationary, .true., connected_state, path_certificate, &
    status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, &
                        'connected Zhao-A forward path failed: '//trim(message))
  call assert_true( &
    path_certificate%target_reached .and. .not. path_certificate%parameter_fold_detected .and. &
    .not. path_certificate%nonmonotone_barrier_detected .and. &
    path_certificate%target_bracketed .and. path_certificate%target_roundtrip_verified .and. &
    path_certificate%accepted_points > 1_i32, &
    'connected Zhao-A forward path did not produce a regular target certificate' &
    )
  call assert_close_dp(connected_state%interface_field, target_field, &
                       1.0e-12_dp*abs(target_field), &
                       'connected Zhao-A path did not land at its exact field')
  analytic_barrier_j = abs(options%photoelectron_charge)* &
                       (connected_state%zhao_phi0 - min( &
                        connected_state%zhao_phi_minimum, 0.0_dp &
                        ))
  call assert_close_dp(zhao_profile_barrier_energy( &
                       connected_state, options%photoelectron_charge), analytic_barrier_j, &
                       1.0e-11_dp*options%photoelectron_temperature_j, &
                       'connected Zhao-A profile and analytic barriers disagree')
  call continue_outer_plasma_zhao_connected_root( &
    connected_options, stationary, .true., connected_root, path_certificate, &
    status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, &
                        'root-only connected Zhao-A path failed: '//trim(message))
  call assert_true(path_certificate%target_reached .and. &
                   path_certificate%target_roundtrip_verified, &
                   'root-only connected Zhao-A path lost its target certificate')
  call assert_close_dp( &
    zhao_charge_root_barrier_energy(connected_root, options%photoelectron_charge), &
    analytic_barrier_j, 1.0e-11_dp*options%photoelectron_temperature_j, &
    'root-only connected Zhao-A barrier changed' &
    )
  call materialize_outer_plasma_zhao_root( &
    connected_options, connected_root, root_state, status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, &
                        'accepted Zhao-A root materialization failed: '//trim(message))
  call assert_close_dp( &
    zhao_profile_barrier_energy(root_state, options%photoelectron_charge), &
    zhao_charge_root_barrier_energy(connected_root, options%photoelectron_charge), &
    1.0e-11_dp*options%photoelectron_temperature_j, &
    'materialized Zhao-A profile disagrees with its exact root barrier' &
    )
  call assert_close_dp(maxval(abs(root_state%potential - connected_state%potential)), &
                       0.0_dp, 1.0e-12_dp*options%photoelectron_temperature_j/qe, &
                       'root-only and full Zhao-A profiles differ')
  connected_options%interface_field = stationary%interface_field
  call continue_outer_plasma_zhao_connected( &
    connected_options, connected_state, .true., reverse_state, path_certificate, &
    status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, &
                        'connected Zhao-A reverse path failed: '//trim(message))
  call assert_close_dp(reverse_state%zhao_phi0, stationary%zhao_phi0, 2.0e-8_dp, &
                       'connected Zhao-A round trip changed phi0')
  call assert_close_dp(reverse_state%zhao_phi_minimum, stationary%zhao_phi_minimum, &
                       2.0e-8_dp, 'connected Zhao-A round trip changed phi_min')
  call test_end()

  call test_begin('connected Zhao source and field homotopy reaches the same A sheet')
  connected_options = options
  connected_options%photoelectron_source_scale = 1.001_dp
  connected_options%interface_field = 1.001_dp*stationary%interface_field
  call continue_outer_plasma_zhao_connected( &
    connected_options, stationary, .false., connected_state, path_certificate, &
    status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, &
                        'connected Zhao source-field path failed: '//trim(message))
  call assert_true(path_certificate%target_reached .and. connected_state%zhao_branch == 'A', &
                   'connected source-field homotopy did not retain the Zhao-A sheet')
  call assert_true( &
    path_certificate%target_bracketed .and. path_certificate%target_roundtrip_verified, &
    'connected source-field homotopy did not certify its target event' &
    )
  call assert_close_dp(connected_state%photoelectron_source_scale, 1.001_dp, 1.0e-14_dp, &
                       'connected source-field homotopy lost its target source scale')
  call test_end()

  call test_begin('connected Zhao-A path stops before a field fold')
  connected_options = options
  connected_options%interface_field = 0.85_dp*stationary%interface_field
  call continue_outer_plasma_zhao_connected( &
    connected_options, stationary, .true., connected_state, path_certificate, &
    connected_status, message &
    )
  call assert_equal_i32(connected_status, outer_plasma_no_physical_solution, &
                        'connected Zhao-A path crossed its field fold')
  call assert_true( &
    path_certificate%parameter_fold_detected .and. .not. path_certificate%target_reached .and. &
    trim(path_certificate%reason) == 'parameter_fold_before_target', &
    'connected Zhao-A fold lost its fail-closed classification' &
    )
  call assert_true(.not. connected_state%ready .and. &
                   .not. allocated(connected_state%potential), &
                   'failed connected Zhao path exposed a partial profile')
  call assert_true( &
    path_certificate%minimum_row_rank_indicator > 1.0e-2_dp .and. &
    path_certificate%minimum_lambda_tangent <= 1.1e-6_dp .and. &
    path_certificate%final_residual <= 1.0e-11_dp, &
    'connected Zhao-A fold fixture no longer has a regular augmented curve' &
    )
  call test_end()

  call test_begin('connected Zhao path rejects a branch-incompatible target field sign')
  connected_options = options
  connected_options%interface_field = -abs(stationary%interface_field)
  call continue_outer_plasma_zhao_connected( &
    connected_options, stationary, .true., connected_state, path_certificate, &
    connected_status, message &
    )
  call assert_equal_i32(connected_status, outer_plasma_no_physical_solution, &
                        'connected Zhao-A path accepted a negative interface field')
  call assert_true( &
    .not. path_certificate%target_reached .and. &
    trim(path_certificate%reason) == 'branch_field_sign_mismatch', &
    'connected Zhao-A sign rejection lost its physical classification' &
    )
  connected_options%interface_field = 0.0_dp
  call continue_outer_plasma_zhao_connected( &
    connected_options, stationary, .true., connected_state, path_certificate, &
    connected_status, message &
    )
  call assert_equal_i32(connected_status, outer_plasma_no_physical_solution, &
                        'connected Zhao-A path accepted a zero interface field')
  call assert_true( &
    .not. path_certificate%target_reached .and. &
    trim(path_certificate%reason) == 'branch_field_sign_mismatch', &
    'connected Zhao-A zero-field rejection lost its physical classification' &
    )
  call test_end()

  reference_current_density = abs(options%photoelectron_charge)* &
                              options%photoelectron_reference_density* &
                              sin(options%zhao_alpha_deg*pi/180.0_dp)* &
                              sqrt(2.0_dp*options%photoelectron_temperature_j/ &
                                   options%photoelectron_mass)/(2.0_dp*sqrt(pi))
  emitted_charge_c = reference_current_density*area_xy_m2*time_step_s

  call test_begin('empirical Zhao update recovers an integer escape-set fixed point')
  call assert_close_dp(stationary%photoelectron_current_density/reference_current_density, &
                       high_fraction, 2.0e-8_dp, &
                       'stationary Zhao-A current does not match its analytic barrier fraction')
  anchor_charge_c = eps0*area_xy_m2*stationary%interface_field
  sample_energy_j = [barrier_j + 0.5_dp*qe, barrier_j - 0.5_dp*qe, &
                     barrier_j - 0.5_dp*qe, barrier_j - 0.5_dp*qe]
  sample_charge_c = emitted_charge_c*[high_fraction, &
                                      (1.0_dp - high_fraction)/3.0_dp, &
                                      (1.0_dp - high_fraction)/3.0_dp, &
                                      (1.0_dp - high_fraction)/3.0_dp]
  call build_measured_interface_energy_distribution( &
    sample_energy_j, sample_charge_c, distribution, status, message &
    )
  call assert_equal_i32(status, dynamic_k0_ok, 'integer escape-set distribution failed: '//trim(message))
  base_charge_c = anchor_charge_c - high_fraction*emitted_charge_c
  call advance_dynamic_k0_zhao( &
    options, stationary, 'e_bottom_zero', area_xy_m2, anchor_charge_c, &
    base_charge_c, time_step_s, distribution, step, accepted, effective_source_scale, message &
    )
  call assert_equal_i32(step%status, dynamic_k0_ok, 'integer Zhao fixed point failed: '//trim(message))
  call assert_close_dp(effective_source_scale, 1.0_dp, 2.0e-14_dp, &
                       'measured source normalization did not recover unity')
  call assert_close_dp(step%surface_charge_after_c, anchor_charge_c, &
                       1.0e-8_dp*emitted_charge_c, &
                       'integer Zhao fixed point did not recover the anchor charge')
  call assert_close_dp(step%photoelectron_escape_fraction, high_fraction, 1.0e-8_dp, &
                       'integer Zhao escape fraction changed')
  call assert_true(step%zhao_branch == 'A' .and. accepted%zhao_branch == 'A', &
                   'integer Zhao update left the committed branch')
  call assert_close_dp(step%photoelectron_barrier_energy_j/qe, barrier_ev, 2.0e-5_dp, &
                       'integer Zhao update changed the anchor barrier')
  call assert_close_dp(measured_sample_escape_fraction(sample_energy_j(1), step), &
                       1.0_dp, 0.0_dp, 'above-barrier measured ray did not escape')
  call assert_close_dp(measured_sample_escape_fraction(sample_energy_j(2), step), &
                       0.0_dp, 0.0_dp, 'below-barrier measured ray did not return')
  call assert_close_dp(accepted%photoelectron_current_density, &
                       high_fraction*reference_current_density, &
                       1.0e-8_dp*reference_current_density, &
                       'accepted Zhao current did not use the measured escaped charge')
  call assert_close_dp(accepted%total_current_density, &
                       accepted%electron_current_density + accepted%ion_current_density + &
                       accepted%photoelectron_current_density, &
                       1.0e-12_dp*reference_current_density, &
                       'accepted Zhao total current is inconsistent with its components')
  call test_end()

  call test_begin('large empirical CDF matches a pure breakpoint oracle in logarithmic solves')
  allocate (many_energy_j(512), many_charge_c(512))
  do sample_index = 1_i32, 512_i32
    many_energy_j(sample_index) = barrier_j + &
                                  (real(213_i32 - sample_index, dp) + 0.5_dp)*1.0e-3_dp*qe
  end do
  many_charge_c = emitted_charge_c/512.0_dp
  call build_measured_interface_energy_distribution( &
    many_energy_j, many_charge_c, distribution, status, message &
    )
  call assert_equal_i32(status, dynamic_k0_ok, 'large empirical CDF build failed: '//trim(message))
  call assert_equal_i32(distribution%group_count, 512_i32, 'unique empirical energies were merged')
  base_charge_c = anchor_charge_c - 213.0_dp*emitted_charge_c/512.0_dp
  call advance_dynamic_k0_zhao( &
    options, stationary, 'e_bottom_zero', area_xy_m2, anchor_charge_c, &
    base_charge_c, time_step_s, distribution, step, accepted, effective_source_scale, message &
    )
  call assert_equal_i32(step%status, dynamic_k0_ok, 'large empirical CDF root failed: '//trim(message))
  call assert_close_dp(step%surface_charge_after_c, anchor_charge_c, &
                       1.0e-8_dp*emitted_charge_c, &
                       'large empirical CDF did not recover the analytic charge')
  call assert_close_dp(step%photoelectron_escape_fraction, 213.0_dp/512.0_dp, 1.0e-12_dp, &
                       'large empirical CDF selected the wrong escape set')
  call assert_true(step%iterations >= 12_i32 .and. step%iterations <= 13_i32, &
                   'large empirical CDF did not use logarithmic order-statistic selection')
  deallocate (many_energy_j, many_charge_c)
  call test_end()

  call test_begin('empirical Zhao update resolves a fractional marginal macro group')
  anchor_charge_c = 2.0_dp*eps0*area_xy_m2*stationary%interface_field
  marginal_fraction = 0.35_dp
  sample_energy_j = barrier_j
  sample_charge_c = 0.25_dp*emitted_charge_c
  call build_measured_interface_energy_distribution( &
    sample_energy_j, sample_charge_c, distribution, status, message &
    )
  call assert_equal_i32(status, dynamic_k0_ok, 'marginal distribution build failed: '//trim(message))
  base_charge_c = anchor_charge_c - marginal_fraction*emitted_charge_c
  call advance_dynamic_k0_zhao( &
    options, stationary, 'symmetric_vacuum', area_xy_m2, anchor_charge_c, &
    base_charge_c, time_step_s, distribution, step, accepted, effective_source_scale, message &
    )
  call assert_equal_i32(step%status, dynamic_k0_ok, 'marginal Zhao fixed point failed: '//trim(message))
  call assert_close_dp(step%surface_charge_after_c, anchor_charge_c, &
                       5.0e-7_dp*emitted_charge_c, &
                       'marginal Zhao root did not recover the anchor charge')
  call assert_close_dp(step%marginal_photoelectron_escape_fraction, marginal_fraction, 5.0e-7_dp, &
                       'marginal Zhao macro fraction changed')
  call assert_close_dp(measured_sample_escape_fraction(barrier_j, step), marginal_fraction, &
                       5.0e-7_dp, 'marginal measured ray did not receive its solved fraction')
  call assert_close_dp(step%backward_euler_residual_charge_c, 0.0_dp, &
                       1.0e-12_dp*emitted_charge_c, &
                       'marginal Zhao transaction did not close charge')
  call test_end()

  call test_begin('tiny marginal energy group converges in charge as well as barrier')
  anchor_charge_c = eps0*area_xy_m2*stationary%interface_field
  tiny_group_charge_c = 1.0e-6_dp*emitted_charge_c
  sample_energy_j = [barrier_j + 0.5_dp*qe, barrier_j, &
                     barrier_j - 0.5_dp*qe, barrier_j - 0.5_dp*qe]
  sample_charge_c = [ &
                    high_fraction*emitted_charge_c, tiny_group_charge_c, &
                    0.5_dp*((1.0_dp - high_fraction)*emitted_charge_c - tiny_group_charge_c), &
                    0.5_dp*((1.0_dp - high_fraction)*emitted_charge_c - tiny_group_charge_c) &
                    ]
  call build_measured_interface_energy_distribution( &
    sample_energy_j, sample_charge_c, distribution, status, message &
    )
  call assert_equal_i32(status, dynamic_k0_ok, 'tiny-marginal distribution build failed: '//trim(message))
  base_charge_c = anchor_charge_c - sample_charge_c(1) - 0.25_dp*tiny_group_charge_c
  call advance_dynamic_k0_zhao( &
    options, stationary, 'e_bottom_zero', area_xy_m2, anchor_charge_c, &
    base_charge_c, time_step_s, distribution, step, accepted, effective_source_scale, message &
    )
  call assert_equal_i32(step%status, dynamic_k0_ok, 'tiny-marginal Zhao root failed: '//trim(message))
  call assert_close_dp(step%marginal_photoelectron_escape_fraction, 0.25_dp, 5.0e-3_dp, &
                       'tiny marginal group was hidden by the barrier tolerance')
  call assert_true(step%photoelectron_barrier_energy_j >= step%marginal_photoelectron_energy_j, &
                   'tiny marginal root did not retain the return-side profile')
  call assert_close_dp(step%surface_charge_after_c, anchor_charge_c, &
                       2.0_dp*max(1.0e-12_dp*emitted_charge_c, tiny_group_charge_c*2.0e-3_dp), &
                       'tiny marginal root did not recover the prescribed charge')
  call test_end()
  anchor_charge_c = 2.0_dp*eps0*area_xy_m2*stationary%interface_field

  call test_begin('all-return k=0 candidate obeys the turning-point equality convention')
  sample_energy_j = barrier_j
  sample_charge_c = 0.25_dp*emitted_charge_c
  call build_measured_interface_energy_distribution( &
    sample_energy_j, sample_charge_c, distribution, status, message &
    )
  call assert_equal_i32(status, dynamic_k0_ok, 'all-return distribution build failed: '//trim(message))
  base_charge_c = anchor_charge_c
  call advance_dynamic_k0_zhao( &
    options, stationary, 'symmetric_vacuum', area_xy_m2, anchor_charge_c, &
    base_charge_c, time_step_s, distribution, step, accepted, effective_source_scale, message &
    )
  call assert_equal_i32(step%status, dynamic_k0_ok, 'all-return Zhao root failed: '//trim(message))
  call assert_close_dp(step%photoelectron_escape_fraction, 0.0_dp, 0.0_dp, &
                       'turning-point equality must not be counted as ordinary escape')
  call assert_close_dp(measured_sample_escape_fraction(barrier_j, step), 0.0_dp, 0.0_dp, &
                       'sample mapper disagrees with the turning-point equality convention')
  call assert_close_dp(step%surface_charge_after_c, anchor_charge_c, &
                       1.0e-12_dp*emitted_charge_c, 'all-return k=0 charge changed')
  call assert_equal_i32(step%iterations, 3_i32, &
                        'all-return endpoint did not stop after the two global certificates')
  call test_end()

  call test_begin('unit marginal endpoint keeps its full escaping statistical weight')
  sample_energy_j = barrier_j
  sample_charge_c = 0.25_dp*emitted_charge_c
  call build_measured_interface_energy_distribution( &
    sample_energy_j, sample_charge_c, distribution, status, message &
    )
  call assert_equal_i32(status, dynamic_k0_ok, 'unit-marginal distribution build failed: '//trim(message))
  base_charge_c = anchor_charge_c - sum(sample_charge_c)
  call advance_dynamic_k0_zhao( &
    options, stationary, 'symmetric_vacuum', area_xy_m2, anchor_charge_c, &
    base_charge_c, time_step_s, distribution, step, accepted, effective_source_scale, message &
    )
  call assert_equal_i32(step%status, dynamic_k0_ok, 'unit-marginal Zhao root failed: '//trim(message))
  call assert_close_dp(step%marginal_photoelectron_escape_fraction, 1.0_dp, 5.0e-7_dp, &
                       'unit marginal did not retain full escape weight')
  call assert_close_dp(measured_sample_escape_fraction(barrier_j, step), 1.0_dp, 5.0e-7_dp, &
                       'unit marginal sample mapper lost its escape weight')
  call test_end()

  call test_begin('multi-group unit marginal endpoint preserves every escaping group')
  sample_energy_j = [ &
                    barrier_j + 0.5_dp*qe, barrier_j + 0.5_dp*qe, &
                    barrier_j, barrier_j &
                    ]
  sample_charge_c = [0.4999995_dp, 0.4999995_dp, 0.0000005_dp, 0.0000005_dp]* &
                    emitted_charge_c
  call build_measured_interface_energy_distribution( &
    sample_energy_j, sample_charge_c, distribution, status, message &
    )
  call assert_equal_i32(status, dynamic_k0_ok, &
                        'multi-group endpoint distribution build failed: '//trim(message))
  call assert_equal_i32(distribution%group_count, 2_i32, &
                        'multi-group endpoint fixture lost its final energy group')
  base_charge_c = anchor_charge_c - sum(sample_charge_c)
  call advance_dynamic_k0_zhao( &
    options, stationary, 'symmetric_vacuum', area_xy_m2, anchor_charge_c, &
    base_charge_c, time_step_s, distribution, step, accepted, effective_source_scale, message &
    )
  call assert_equal_i32(step%status, dynamic_k0_ok, &
                        'multi-group unit-marginal Zhao root failed: '//trim(message))
  call assert_close_dp(step%marginal_photoelectron_escape_fraction, 1.0_dp, 5.0e-7_dp, &
                       'final marginal group lost its full escape weight')
  call assert_close_dp(measured_sample_escape_fraction(barrier_j + 0.5_dp*qe, step), &
                       1.0_dp, 0.0_dp, 'higher energy group was not mapped to escape')
  call assert_close_dp(measured_sample_escape_fraction(barrier_j, step), 1.0_dp, 5.0e-7_dp, &
                       'final marginal group was not mapped to escape')
  call assert_close_dp(step%surface_charge_after_c, anchor_charge_c, &
                       1.0e-12_dp*emitted_charge_c, &
                       'multi-group unit-marginal charge did not close')
  call test_end()

  call test_begin('all-escape k=M candidate closes the opposite CDF endpoint')
  sample_energy_j = barrier_j + 0.5_dp*qe
  sample_charge_c = 0.25_dp*emitted_charge_c
  call build_measured_interface_energy_distribution( &
    sample_energy_j, sample_charge_c, distribution, status, message &
    )
  call assert_equal_i32(status, dynamic_k0_ok, 'all-escape distribution build failed: '//trim(message))
  base_charge_c = anchor_charge_c - emitted_charge_c
  call advance_dynamic_k0_zhao( &
    options, stationary, 'symmetric_vacuum', area_xy_m2, anchor_charge_c, &
    base_charge_c, time_step_s, distribution, step, accepted, effective_source_scale, message &
    )
  call assert_equal_i32(step%status, dynamic_k0_ok, 'all-escape Zhao root failed: '//trim(message))
  call assert_close_dp(step%photoelectron_escape_fraction, 1.0_dp, 0.0_dp, &
                       'all-escape CDF endpoint changed')
  call assert_close_dp(step%surface_charge_after_c, anchor_charge_c, &
                       1.0e-12_dp*emitted_charge_c, 'all-escape k=M charge changed')
  call assert_equal_i32(step%iterations, 3_i32, &
                        'all-escape endpoint did not stop after the two global certificates')
  call test_end()

  call test_begin('measured source-scale change stays on the connected Zhao-A root')
  sample_energy_j = 0.0_dp
  sample_charge_c = 0.25025_dp*emitted_charge_c
  call build_measured_interface_energy_distribution( &
    sample_energy_j, sample_charge_c, distribution, status, message &
    )
  call assert_equal_i32(status, dynamic_k0_ok, 'scaled-source distribution build failed: '//trim(message))
  base_charge_c = anchor_charge_c
  call advance_dynamic_k0_zhao( &
    options, stationary, 'symmetric_vacuum', area_xy_m2, anchor_charge_c, &
    base_charge_c, time_step_s, distribution, step, accepted, effective_source_scale, message &
    )
  call assert_equal_i32(step%status, dynamic_k0_ok, 'scaled-source Zhao continuation failed: '//trim(message))
  call assert_close_dp(effective_source_scale, 1.001_dp, 2.0e-14_dp, &
                       'measured source-scale continuation changed its normalization')
  call assert_true(step%zhao_branch == 'A' .and. accepted%zhao_branch == 'A', &
                   'measured source-scale continuation changed branch')
  call test_end()

  call test_begin('ambient and photoelectron frozen-cohort trust scales remain separated')
  sample_energy_j = 0.0_dp
  sample_charge_c = 0.2575_dp*emitted_charge_c
  call build_measured_interface_energy_distribution( &
    sample_energy_j, sample_charge_c, distribution, status, message &
    )
  call assert_equal_i32(status, dynamic_k0_ok, &
                        'split-scale distribution build failed: '//trim(message))
  base_charge_c = anchor_charge_c
  call advance_dynamic_k0_zhao( &
    options, stationary, 'symmetric_vacuum', area_xy_m2, anchor_charge_c, &
    base_charge_c, time_step_s, distribution, step, accepted, effective_source_scale, message &
    )
  call assert_equal_i32(step%status, dynamic_k0_ok, &
                        'ambient-scaled potential trust check failed: '//trim(message))
  call assert_close_dp(effective_source_scale, 1.03_dp, 2.0e-14_dp, &
                       'split-scale fixture changed its measured source normalization')
  call assert_true( &
    abs(step%interface_potential_after_v - stationary%interface_potential)/ &
    (options%photoelectron_temperature_j/abs(options%photoelectron_charge)) > 0.25_dp, &
    'split-scale fixture no longer exceeds the obsolete Tpe potential guard' &
    )
  call assert_true( &
    abs(step%interface_potential_after_v - stationary%interface_potential)/ &
    (options%electron_temperature_j/abs(options%electron_charge)) < 0.25_dp, &
    'split-scale fixture left the ambient-electron potential trust region' &
    )
  call test_end()

  call test_begin('large frozen-cohort source jump fails closed')
  sample_energy_j = 0.0_dp
  sample_charge_c = 0.325_dp*emitted_charge_c
  call build_measured_interface_energy_distribution( &
    sample_energy_j, sample_charge_c, distribution, status, message &
    )
  call assert_equal_i32(status, dynamic_k0_ok, 'large-jump distribution build failed: '//trim(message))
  call advance_dynamic_k0_zhao( &
    options, stationary, 'symmetric_vacuum', area_xy_m2, anchor_charge_c, &
    base_charge_c, time_step_s, distribution, step, accepted, effective_source_scale, message &
    )
  call assert_equal_i32(step%status, dynamic_zhao_no_physical_root, &
                        'large frozen-cohort source jump was not rejected')
  call assert_true(.not. accepted%ready, &
                   'failed Zhao update exposed a partially accepted outer state')
  call test_end()

  call test_begin('unrepresentable ambient ion energy fails closed before continuation')
  extreme_options = options
  extreme_options%ion_drift_infinity = huge(1.0_dp)
  call advance_dynamic_k0_zhao( &
    extreme_options, stationary, 'symmetric_vacuum', area_xy_m2, anchor_charge_c, &
    base_charge_c, time_step_s, distribution, step, accepted, effective_source_scale, message &
    )
  call assert_equal_i32(step%status, dynamic_zhao_no_physical_root, &
                        'overflowing ambient ion energy was not rejected')
  call assert_true(.not. accepted%ready, &
                   'overflowing ion-energy failure exposed a partially accepted outer state')
  call assert_true(index(message, 'ambient ion energy scale is not representable') > 0, &
                   'overflowing ion-energy failure lost its fail-closed diagnostic')
  call test_end()

  call test_begin('effective closepack source admits the explicit monotonic Zhao-B stationary root')
  call configure_closepack_options(closepack_options)
  call solve_outer_plasma_zhao_stationary(closepack_options, closepack_stationary, status, message)
  call assert_equal_i32(status, outer_plasma_ok, &
                        'effective closepack Zhao stationary solve failed: '//trim(message))
  call assert_true(closepack_stationary%ready .and. closepack_stationary%zhao_branch == 'B', &
                   'explicit effective closepack source did not resolve Zhao-B')
  call assert_close_dp(closepack_stationary%interface_potential, 5.855819944830117_dp, 5.0e-6_dp, &
                       'effective closepack Zhao-B potential changed')
  call assert_close_dp(closepack_stationary%zhao_electron_density_infinity, &
                       4.228122218797317e6_dp, 5.0_dp, &
                       'effective closepack Zhao-B ambient density changed')
  call test_end()

  call test_begin('connected closepack Zhao-B field path stays regular and monotone')
  connected_options = closepack_options
  connected_options%interface_field = 1.005_dp*closepack_stationary%interface_field
  call continue_outer_plasma_zhao_connected( &
    connected_options, closepack_stationary, .true., connected_state, path_certificate, &
    status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, &
                        'connected closepack Zhao-B path failed: '//trim(message))
  call assert_true( &
    path_certificate%target_reached .and. connected_state%zhao_branch == 'B' .and. &
    path_certificate%target_bracketed .and. path_certificate%target_roundtrip_verified .and. &
    path_certificate%minimum_lambda_tangent > 0.0_dp .and. &
    path_certificate%minimum_fixed_parameter_rank_indicator > 0.0_dp, &
    'connected closepack Zhao-B path lost its regular sheet certificate' &
    )
  analytic_barrier_j = abs(closepack_options%photoelectron_charge)* &
                       max(connected_state%zhao_phi0, 0.0_dp)
  call assert_close_dp(zhao_profile_barrier_energy( &
                       connected_state, closepack_options%photoelectron_charge), &
                       analytic_barrier_j, &
                       1.0e-11_dp*closepack_options%photoelectron_temperature_j, &
                       'connected Zhao-B profile and analytic barriers disagree')
  connected_root = zhao_charge_root_type( &
                   branch='B', phi0_v=connected_state%zhao_phi0, &
                   phi_m_v=connected_state%zhao_phi_minimum &
                   )
  call assert_close_dp( &
    zhao_charge_root_barrier_energy(connected_root, closepack_options%photoelectron_charge), &
    zhao_profile_barrier_energy(connected_state, closepack_options%photoelectron_charge), &
    1.0e-11_dp*closepack_options%photoelectron_temperature_j, &
    'connected Zhao-B exact root and profile barriers disagree' &
    )
  connected_options%interface_field = -abs(closepack_stationary%interface_field)
  call continue_outer_plasma_zhao_connected( &
    connected_options, closepack_stationary, .true., connected_state, path_certificate, &
    connected_status, message &
    )
  call assert_equal_i32(connected_status, outer_plasma_no_physical_solution, &
                        'connected Zhao-B path accepted a negative interface field')
  call assert_true(trim(path_certificate%reason) == 'branch_field_sign_mismatch', &
                   'connected Zhao-B sign rejection lost its physical classification')
  call test_end()

  call test_begin('connected Zhao-C path retains its zero outward electron barrier')
  call configure_type_c_options(type_c_options)
  call solve_outer_plasma_zhao_stationary(type_c_options, type_c_stationary, status, message)
  call assert_equal_i32(status, outer_plasma_ok, &
                        'stationary Zhao-C setup failed: '//trim(message))
  call assert_true(type_c_stationary%zhao_branch == 'C' .and. &
                   type_c_stationary%interface_field < 0.0_dp, &
                   'stationary Zhao-C setup selected the wrong branch')
  connected_options = type_c_options
  connected_options%interface_field = 1.005_dp*type_c_stationary%interface_field
  call continue_outer_plasma_zhao_connected( &
    connected_options, type_c_stationary, .true., connected_state, path_certificate, &
    status, message &
    )
  call assert_equal_i32(status, outer_plasma_ok, &
                        'connected Zhao-C path failed: '//trim(message))
  call assert_true( &
    path_certificate%target_reached .and. connected_state%zhao_branch == 'C' .and. &
    path_certificate%target_bracketed .and. path_certificate%target_roundtrip_verified, &
    'connected Zhao-C path left its committed sheet or skipped target certification' &
    )
  call assert_close_dp(zhao_profile_barrier_energy( &
                       connected_state, type_c_options%photoelectron_charge), &
                       0.0_dp, 1.0e-30_dp, &
                       'connected Zhao-C outward electron barrier is not zero')
  connected_root = zhao_charge_root_type( &
                   branch='C', phi0_v=connected_state%zhao_phi0, &
                   phi_m_v=connected_state%zhao_phi_minimum &
                   )
  call assert_close_dp( &
    zhao_charge_root_barrier_energy(connected_root, type_c_options%photoelectron_charge), &
    zhao_profile_barrier_energy(connected_state, type_c_options%photoelectron_charge), &
    1.0e-30_dp, 'connected Zhao-C exact root and profile barriers disagree' &
    )
  connected_options%interface_field = abs(type_c_stationary%interface_field)
  call continue_outer_plasma_zhao_connected( &
    connected_options, type_c_stationary, .true., connected_state, path_certificate, &
    connected_status, message &
    )
  call assert_equal_i32(connected_status, outer_plasma_no_physical_solution, &
                        'connected Zhao-C path accepted a positive interface field')
  call assert_true(trim(path_certificate%reason) == 'branch_field_sign_mismatch', &
                   'connected Zhao-C sign rejection lost its physical classification')
  call test_end()

  call test_begin('invalid public distributions fail before Zhao continuation')
  distribution = measured_interface_energy_distribution_type(ready=.true., group_count=1_i32)
  call advance_dynamic_k0_zhao( &
    options, stationary, 'symmetric_vacuum', area_xy_m2, anchor_charge_c, &
    anchor_charge_c, time_step_s, distribution, step, accepted, effective_source_scale, message &
    )
  call assert_true(step%status /= dynamic_k0_ok, 'incomplete measured distribution was accepted')
  sample_energy_j = barrier_j
  sample_charge_c = 0.25_dp*emitted_charge_c
  call build_measured_interface_energy_distribution( &
    sample_energy_j, sample_charge_c, distribution, status, message &
    )
  distribution%cumulative_charge_c(1) = 0.5_dp*distribution%group_charge_c(1)
  distribution%total_charge_c = distribution%cumulative_charge_c(1)
  call advance_dynamic_k0_zhao( &
    options, stationary, 'symmetric_vacuum', area_xy_m2, anchor_charge_c, &
    anchor_charge_c, time_step_s, distribution, step, accepted, effective_source_scale, message &
    )
  call assert_true(step%status /= dynamic_k0_ok, &
                   'cumulative/group-charge mismatch was accepted')
  call test_end()

  call test_summary()

contains

  subroutine configure_type_a_options(resolved)
    type(kinetic_outer_plasma_options_type), intent(out) :: resolved

    resolved = kinetic_outer_plasma_options_type()
    resolved%kinetic_closure = 'zhao_charge_driven'
    resolved%zhao_branch = 'a'
    resolved%grid_points = 257_i32
    resolved%domain_length = 2.0_dp
    resolved%electron_charge = -qe
    resolved%electron_mass = electron_mass
    resolved%electron_temperature_j = 12.0_dp*qe
    resolved%electron_drift_infinity = 4.68e5_dp*sin(60.0_dp*pi/180.0_dp)
    resolved%ion_charge = qe
    resolved%ion_mass = proton_mass
    resolved%ion_density_infinity = 8.7e6_dp
    resolved%ion_temperature_j = 0.1_dp*qe
    resolved%ion_drift_infinity = resolved%electron_drift_infinity
    resolved%photoelectron_charge = -qe
    resolved%photoelectron_mass = electron_mass
    resolved%photoelectron_temperature_j = 2.2_dp*qe
    resolved%photoelectron_reference_density = 64.0e6_dp
    resolved%photoelectron_population_fraction = 1.0_dp
    resolved%photoelectron_source_scale = 1.0_dp
    resolved%photoelectron_column_closure_enabled = .false.
    resolved%zhao_alpha_deg = 60.0_dp
  end subroutine configure_type_a_options

  subroutine configure_closepack_options(resolved)
    type(kinetic_outer_plasma_options_type), intent(out) :: resolved

    call configure_type_a_options(resolved)
    resolved%zhao_branch = 'b'
    resolved%electron_density_infinity = 5.0e6_dp
    resolved%electron_temperature_j = 10.0_dp*qe
    resolved%electron_drift_infinity = 4.0e5_dp
    resolved%ion_density_infinity = 5.0e6_dp
    resolved%ion_drift_infinity = 4.0e5_dp
    resolved%photoelectron_reference_density = 78.87291613928085e6_dp
  end subroutine configure_closepack_options

  subroutine configure_type_c_options(resolved)
    type(kinetic_outer_plasma_options_type), intent(out) :: resolved

    call configure_type_a_options(resolved)
    resolved%zhao_branch = 'c'
    resolved%zhao_alpha_deg = 10.0_dp
    resolved%electron_drift_infinity = 4.68e5_dp*sin(10.0_dp*pi/180.0_dp)
    resolved%ion_drift_infinity = resolved%electron_drift_infinity
  end subroutine configure_type_c_options

end program test_dynamic_k0_zhao

program test_outer_plasma_kinetic
  use bem_kinds, only: dp, i32
  use bem_constants, only: qe, eps0
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok, outer_plasma_no_physical_solution
  use bem_outer_plasma_kinetic, only: kinetic_outer_plasma_options_type, solve_outer_plasma_kinetic, &
                                      eval_kinetic_residual_jacobian_action
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, assert_equal_i32
  implicit none

  real(dp), parameter :: electron_mass = 9.1093837139e-31_dp
  type(kinetic_outer_plasma_options_type) :: options
  type(outer_plasma_state_type) :: state
  type(outer_plasma_state_type) :: ambient_state
  type(outer_plasma_state_type) :: coarse_state
  real(dp), allocatable :: direction(:), residual(:), jacobian_action(:), residual_plus(:), residual_minus(:)
  real(dp), parameter :: jacobian_step = 1.0e-5_dp
  integer(i32) :: status
  integer(i32) :: continuation_steps
  character(len=256) :: message

  call test_init(11)

  call test_begin('vacuum Neumann Robin problem matches its analytic solution')
  options = reference_options()
  options%electron_density_infinity = 0.0_dp
  options%ion_density_infinity = 0.0_dp
  options%interface_field = -0.25_dp
  options%domain_length = 3.0_dp
  options%tail_length = 2.0_dp
  call solve_outer_plasma_kinetic(options, state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'vacuum solve status mismatch')
  call assert_true(state%ready, 'vacuum state must be ready')
  call assert_close_dp(state%potential(1), -1.25_dp, 1.0e-10_dp, 'vacuum interface potential mismatch')
  call assert_close_dp(state%potential(state%profile_n), -0.5_dp, 1.0e-10_dp, 'vacuum far potential mismatch')
  call assert_true(state%nonlinear_residual < 1.0e-10_dp, 'vacuum residual is too large')
  call assert_close_dp( &
    state%integrated_charge_per_area, eps0*(state%field(state%profile_n) - state%interface_field), &
    1.0e-20_dp, 'finite-domain Gauss closure mismatch' &
    )
  call test_end()

  call test_begin('analytic kinetic Jacobian action matches residual differences')
  options = reference_options()
  options%interface_field = -0.02_dp
  call solve_outer_plasma_kinetic(options, state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'Jacobian source solve failed')
  allocate (direction(state%profile_n), residual(state%profile_n), jacobian_action(state%profile_n), &
            residual_plus(state%profile_n), residual_minus(state%profile_n))
  direction = 0.1_dp*state%potential
  call eval_kinetic_residual_jacobian_action( &
    options, state%potential, direction, residual, jacobian_action, status &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'analytic Jacobian action failed')
  call eval_kinetic_residual_jacobian_action( &
    options, state%potential + jacobian_step*direction, direction, residual_plus, jacobian_action, status &
    )
  call eval_kinetic_residual_jacobian_action( &
    options, state%potential - jacobian_step*direction, direction, residual_minus, jacobian_action, status &
    )
  call eval_kinetic_residual_jacobian_action( &
    options, state%potential, direction, residual, jacobian_action, status &
    )
  call assert_true( &
    maxval(abs(jacobian_action - (residual_plus - residual_minus)/(2.0_dp*jacobian_step))) <= &
    2.0e-6_dp*max(1.0_dp, maxval(abs(jacobian_action))), &
    'analytic Jacobian action mismatch' &
    )
  deallocate (direction, residual, jacobian_action, residual_plus, residual_minus)
  call test_end()

  call test_begin('nonlinear profile is stable under grid refinement')
  options = reference_options()
  options%interface_field = -0.02_dp
  call solve_outer_plasma_kinetic(options, coarse_state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'coarse kinetic solve status mismatch')
  options%grid_points = 65_i32
  call solve_outer_plasma_kinetic(options, state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'refined kinetic solve status mismatch')
  call assert_close_dp( &
    state%interface_potential, coarse_state%interface_potential, 2.0e-3_dp, &
    'kinetic interface potential did not converge with grid refinement' &
    )
  call assert_close_dp( &
    state%integrated_charge_per_area, eps0*(state%field(state%profile_n) - state%interface_field), &
    5.0e-15_dp, 'refined nonlinear Gauss closure mismatch' &
    )
  call test_end()

  call test_begin('lunar warm start crosses the former runtime failure field')
  options = reference_options()
  options%grid_points = 128_i32
  options%domain_length = 105.132_dp
  options%tail_length = 10.5132_dp
  options%electron_density_infinity = 5.0e6_dp
  options%electron_temperature_j = 10.0_dp*qe
  options%ion_density_infinity = 5.0e6_dp
  options%ion_temperature_j = 10.0_dp*qe
  options%ion_mass = 1.672482821616e-27_dp
  options%ion_drift_infinity = 4.0e5_dp
  options%max_iterations = 40_i32
  options%residual_tolerance = 1.0e-8_dp
  options%interface_field = 0.0_dp
  call solve_outer_plasma_kinetic(options, ambient_state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'lunar continuation source solve failed: '//trim(message))
  options%interface_field = -0.70_dp
  call solve_outer_plasma_kinetic( &
    options, state, status, message, initial_potential=ambient_state%potential, &
    continuation_steps=continuation_steps &
    )
  call assert_equal_i32(status, outer_plasma_ok, 'large lunar field change failed: '//trim(message))
  call assert_true(continuation_steps > 1_i32, 'large lunar field change did not use continuation')
  ambient_state = state
  options%interface_field = -0.72898324579369622_dp
  call solve_outer_plasma_kinetic(options, state, status, message, initial_potential=ambient_state%potential)
  call assert_equal_i32(status, outer_plasma_ok, 'former runtime failure field did not converge: '//trim(message))
  call assert_true(state%nonlinear_residual <= options%residual_tolerance, 'runtime failure field residual is too large')
  call assert_true( &
    all(state%potential(2:) >= state%potential(:state%profile_n - 1_i32)), &
    'runtime failure field left the monotonic branch' &
    )
  call test_end()

  call test_begin('strong photoelectron space charge changes the positive monotonic branch')
  options = reference_options()
  options%interface_field = 0.05_dp
  options%ion_drift_infinity = 4.0_dp*sqrt(options%electron_temperature_j/options%ion_mass)
  call solve_outer_plasma_kinetic(options, ambient_state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'ambient positive-branch solve status mismatch')
  options%photoelectron_charge = -qe
  options%photoelectron_mass = electron_mass
  options%photoelectron_temperature_j = 1.5_dp*qe
  options%photoelectron_emission_flux = 5.0e10_dp
  call solve_outer_plasma_kinetic(options, state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'strong-photoelectron solve status mismatch: '//trim(message))
  call assert_true(state%interface_potential > 0.0_dp, 'positive monotonic branch was not selected')
  call assert_true( &
    abs(state%interface_potential - ambient_state%interface_potential) > 1.0e-6_dp, &
    'strong photoelectron density did not affect the sheath profile' &
    )
  call assert_true(state%nonlinear_residual < options%residual_tolerance, 'strong-photoelectron residual is too large')
  call test_end()

  call test_begin('photoelectron space charge starts from a neutral surface field')
  options = reference_options()
  options%interface_field = 0.0_dp
  options%photoelectron_charge = -qe
  options%photoelectron_mass = electron_mass
  options%photoelectron_temperature_j = 1.5_dp*qe
  options%photoelectron_emission_flux = 5.0e10_dp
  call solve_outer_plasma_kinetic(options, state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'zero-field photoelectron solve status mismatch: '//trim(message))
  call assert_true(state%interface_potential < 0.0_dp, 'photoelectron space charge must lower the interface potential')
  call assert_true( &
    all(state%potential(2:) >= state%potential(:state%profile_n - 1_i32)), &
    'zero-field photoelectron profile must remain monotonic' &
    )
  call test_end()

  call test_begin('quasineutral zero field is an exact kinetic equilibrium')
  options = reference_options()
  options%interface_field = 0.0_dp
  call solve_outer_plasma_kinetic(options, state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'quasineutral solve status mismatch')
  call assert_true(maxval(abs(state%potential)) < 1.0e-10_dp, 'quasineutral potential must vanish')
  call assert_true(maxval(abs(state%charge_density)) < 1.0e-12_dp, 'quasineutral charge must vanish')
  call test_end()

  call test_begin('warm start follows a changed interface field')
  options = reference_options()
  options%interface_field = 0.0_dp
  call solve_outer_plasma_kinetic(options, ambient_state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'zero-field warm-start source solve failed')
  options%interface_field = -0.002_dp
  call solve_outer_plasma_kinetic(options, state, status, message, initial_potential=ambient_state%potential)
  call assert_equal_i32(status, outer_plasma_ok, 'changed-field warm start failed: '//trim(message))
  call assert_true(state%ready, 'changed-field warm start must produce a ready state')
  call test_end()

  call test_begin('warm start crosses the monotonic branch sign')
  options = reference_options()
  options%interface_field = 0.05_dp
  options%ion_drift_infinity = 4.0_dp*sqrt(options%electron_temperature_j/options%ion_mass)
  call solve_outer_plasma_kinetic(options, ambient_state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'positive-branch warm-start source solve failed')
  options%interface_field = -0.02_dp
  call solve_outer_plasma_kinetic(options, state, status, message, initial_potential=ambient_state%potential)
  call assert_equal_i32(status, outer_plasma_ok, 'branch-crossing warm start failed: '//trim(message))
  call assert_true(state%ready, 'branch-crossing warm start must produce a ready state')
  call test_end()

  call test_begin('weak lunar-field warm start keeps a valid Jacobian')
  options = reference_options()
  options%grid_points = 128_i32
  options%domain_length = 105.132_dp
  options%tail_length = 10.5132_dp
  options%electron_density_infinity = 5.0e6_dp
  options%electron_temperature_j = 10.0_dp*qe
  options%ion_density_infinity = 5.0e6_dp
  options%ion_temperature_j = 10.0_dp*qe
  options%ion_mass = 1.672482821616e-27_dp
  options%ion_drift_infinity = 4.0e5_dp
  options%max_iterations = 40_i32
  options%residual_tolerance = 1.0e-8_dp
  options%interface_field = 0.0_dp
  call solve_outer_plasma_kinetic(options, ambient_state, status, message)
  call assert_equal_i32(status, outer_plasma_ok, 'lunar zero-field source solve failed')
  options%interface_field = -1.9340193196540115e-3_dp
  call solve_outer_plasma_kinetic(options, state, status, message, initial_potential=ambient_state%potential)
  call assert_equal_i32(status, outer_plasma_ok, 'weak lunar-field warm start failed: '//trim(message))
  call assert_true(state%ready, 'weak lunar-field warm start must produce a ready state')
  call test_end()

  call test_begin('sub-Bohm ion entry has no physical solution')
  options = reference_options()
  options%ion_drift_infinity = 0.5_dp*sqrt( &
                               (options%electron_temperature_j + options%ion_temperature_j)/options%ion_mass &
                               )
  call solve_outer_plasma_kinetic(options, state, status, message)
  call assert_equal_i32(status, outer_plasma_no_physical_solution, 'sub-Bohm entry must fail physically')
  call assert_true(.not. state%ready, 'sub-Bohm state must not be ready')
  call test_end()

  call test_summary()

contains

  function reference_options() result(value)
    type(kinetic_outer_plasma_options_type) :: value

    value%grid_points = 33_i32
    value%domain_length = 8.0_dp
    value%grid_stretch = 1.5_dp
    value%tail_length = 1.0_dp
    value%interface_field = 0.0_dp
    value%electron_charge = -qe
    value%electron_mass = electron_mass
    value%electron_density_infinity = 1.0e6_dp
    value%electron_temperature_j = 2.0_dp*qe
    value%ion_charge = qe
    value%ion_mass = 1836.0_dp*electron_mass
    value%ion_density_infinity = 1.0e6_dp
    value%ion_temperature_j = 0.0_dp
    value%ion_gamma = 1.0_dp
    value%ion_drift_infinity = 2.0_dp*sqrt(value%electron_temperature_j/value%ion_mass)
    value%max_iterations = 30_i32
    value%residual_tolerance = 1.0e-10_dp
  end function reference_options

end program test_outer_plasma_kinetic

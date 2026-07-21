program test_outer_plasma_zhao
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok, &
                                    outer_plasma_no_physical_solution
  use bem_sheath_model_core, only: zhao_params_type, build_zhao_params, solve_zhao_unknowns
  use bem_outer_plasma_zhao, only: zhao_charge_root_type, solve_zhao_charge_root, &
                                   evaluate_zhao_interface_field, build_zhao_outer_profile
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, &
                          assert_close_dp, assert_equal_i32
  implicit none

  integer(i32), parameter :: profile_points = 257_i32
  real(dp), parameter :: electron_mass = 9.1093837015e-31_dp
  real(dp), parameter :: proton_mass = 1.67262192369e-27_dp
  type(zhao_params_type) :: params
  type(zhao_charge_root_type) :: stationary, charge_root, perturbed_root, zero_root
  type(outer_plasma_state_type) :: state
  real(dp) :: equilibrium_field, perturbed_field, gauss_expected, gauss_scale
  real(dp) :: phi0_v, phi_m_v, density_m3
  integer(i32) :: status
  character(len=1) :: branch
  character(len=256) :: message

  call test_init(6)

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

  call test_begin('zero-field state provides a finite transient bootstrap')
  call solve_zhao_charge_root('zhao_auto', params, 0.0_dp, zero_root, status, message)
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

end program test_outer_plasma_zhao

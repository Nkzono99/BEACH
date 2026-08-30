!> 場・periodic2・panel設定の正規化とfail-closed契約を検証する。
program test_physics_config_types
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_config, bc_open, bc_periodic
  use bem_physics_config_types, only: &
    field_physics_config, periodic2_physics_config, panel_kernel_config, &
    physics_config_ok, physics_config_invalid_combination, physics_config_unavailable, &
    normalize_legacy_physics_config, validate_active_physics_config
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32
  implicit none

  type(sim_config) :: sim
  type(field_physics_config) :: field
  type(periodic2_physics_config) :: periodic2
  type(panel_kernel_config) :: panel
  integer(i32) :: status
  character(len=256) :: message

  call test_init(6)

  call test_begin('free_field_normalization')
  sim = sim_config()
  sim%field_solver = 'treecode'
  sim%field_normalization = 'mesh'
  sim%field_bc_mode = 'free'
  call normalize_legacy_physics_config(sim, field, periodic2, panel)
  call assert_true(trim(field%backend) == 'treecode', 'free backend mismatch')
  call assert_true(trim(field%normalization) == 'mesh', 'free normalization mismatch')
  call assert_true(trim(panel%kernel_id) == 'triangle_p0_exact_tree_near', 'tree panel kernel mismatch')
  call assert_true(trim(periodic2%nonzero_mode_backend) == 'not_applicable', 'free periodic backend mismatch')
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_ok, 'free triangle config should be valid')
  call test_end()

  call test_begin('upper_panel_fourier_retry_contract')
  call configure_xy_periodic(sim, 'fmm', 'cached_kneq0')
  sim%multiple_box_events_retry_backend = 'upper_panel_fourier'
  call normalize_legacy_physics_config(sim, field, periodic2, panel)
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_ok, 'upper Fourier retry should accept cached kneq0')
  periodic2%nonzero_mode_backend = 'panel_spectral_reference'
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'upper Fourier retry must reject other backends')
  call test_end()

  call test_begin('periodic_far_correction_normalization')
  call configure_xy_periodic(sim, 'fmm', 'none')
  call normalize_legacy_physics_config(sim, field, periodic2, panel)
  call assert_true(trim(periodic2%nonzero_mode_backend) == 'legacy_finite_images', 'finite backend mismatch')
  call assert_true(trim(periodic2%zero_mode_policy) == 'legacy_not_decomposed', 'finite zero mode mismatch')
  call assert_true(trim(periodic2%lower_boundary_model) == 'legacy_implicit', 'finite lower boundary mismatch')
  sim%bc_low = [bc_periodic, bc_open, bc_periodic]
  sim%bc_high = [bc_periodic, bc_open, bc_periodic]
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_unavailable, 'legacy periodic2 requires x/y periodic axes')
  sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  sim%bc_high = [bc_periodic, bc_periodic, bc_open]
  sim%field_periodic_far_correction = 'auto'
  call normalize_legacy_physics_config(sim, field, periodic2, panel)
  call assert_true(trim(periodic2%nonzero_mode_backend) == 'legacy_finite_images', 'auto must select finite images')
  sim%field_periodic_far_correction = 'unknown'
  call normalize_legacy_physics_config(sim, field, periodic2, panel)
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'unknown far correction must fail closed')
  call test_end()

  call test_begin('cached_kneq0_contract')
  call configure_xy_periodic(sim, 'fmm', 'cached_kneq0')
  call normalize_legacy_physics_config(sim, field, periodic2, panel)
  call assert_true(trim(periodic2%nonzero_mode_backend) == 'cached_kneq0', 'cached backend mismatch')
  call assert_true(trim(periodic2%zero_mode_policy) == 'exclude_k0', 'cached zero mode mismatch')
  call assert_true(trim(periodic2%lower_boundary_model) == 'e_bottom_zero', 'cached lower boundary mismatch')
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_ok, 'cached kneq0 config should be valid')
  sim%bc_low = [bc_periodic, bc_open, bc_periodic]
  sim%bc_high = [bc_periodic, bc_open, bc_periodic]
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_unavailable, 'cached kneq0 requires x/y periodic axes')
  sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  sim%bc_high = [bc_periodic, bc_periodic, bc_open]
  periodic2%zero_mode_policy = 'legacy_not_decomposed'
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'cached kneq0 requires split zero mode')
  periodic2%zero_mode_policy = 'exclude_k0'
  periodic2%lower_boundary_model = 'symmetric_vacuum'
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_ok, 'symmetric vacuum should be valid')
  periodic2%lower_boundary_model = 'unknown'
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'unknown lower model must fail closed')
  call test_end()

  call test_begin('panel_spectral_reference_contract')
  call configure_xy_periodic(sim, 'direct', 'none')
  call normalize_legacy_physics_config(sim, field, periodic2, panel)
  periodic2%nonzero_mode_backend = 'panel_spectral_reference'
  periodic2%zero_mode_policy = 'exclude_k0'
  periodic2%lower_boundary_model = 'e_bottom_zero'
  periodic2%max_nonzero_mode_potential_step = 0.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_ok, 'panel spectral reference should remain available')
  periodic2%max_nonzero_mode_potential_step = 1.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'panel spectral adaptive limit must fail closed')
  call test_end()

  call test_begin('adaptive_limit_must_be_finite')
  call configure_xy_periodic(sim, 'fmm', 'cached_kneq0')
  call normalize_legacy_physics_config(sim, field, periodic2, panel)
  periodic2%max_nonzero_mode_potential_step = ieee_value(0.0_dp, ieee_quiet_nan)
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'NaN adaptive limit must fail closed')
  periodic2%max_nonzero_mode_potential_step = -1.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'negative adaptive limit must fail closed')
  call test_end()

  call test_summary()

contains

  subroutine configure_xy_periodic(config, solver, far_correction)
    type(sim_config), intent(out) :: config
    character(len=*), intent(in) :: solver, far_correction

    config = sim_config()
    config%field_solver = solver
    config%field_bc_mode = 'periodic2'
    config%field_periodic_far_correction = far_correction
    config%use_box = .true.
    config%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    config%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    config%bc_low = [bc_periodic, bc_periodic, bc_open]
    config%bc_high = [bc_periodic, bc_periodic, bc_open]
  end subroutine configure_xy_periodic

end program test_physics_config_types

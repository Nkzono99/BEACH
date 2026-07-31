!> 場・periodic2・panel設定の正規化とfail-closed契約を検証する。
program test_physics_config_types
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_config, bc_open, bc_periodic
  use bem_physics_config_types, only: &
    field_physics_config, periodic2_physics_config, panel_kernel_config, &
    physics_config_ok, physics_config_invalid_combination, normalize_legacy_physics_config, &
    validate_active_physics_config
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32
  implicit none

  type(sim_config) :: sim
  type(field_physics_config) :: field
  type(periodic2_physics_config) :: periodic2
  type(panel_kernel_config) :: panel
  integer(i32) :: status
  character(len=256) :: message

  call test_init(5)

  call test_begin('free_field_normalization')
  sim = sim_config()
  sim%field_solver = 'treecode'
  sim%field_normalization = 'mesh'
  sim%field_bc_mode = 'free'
  call normalize_legacy_physics_config(sim, field, periodic2, panel)
  call assert_true(trim(field%backend) == 'treecode', 'free backend mismatch')
  call assert_true(trim(field%normalization) == 'mesh', 'free normalization mismatch')
  call assert_true(trim(periodic2%nonzero_mode_backend) == 'not_applicable', 'free periodic backend mismatch')
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_ok, 'free triangle config should be valid')
  call test_end()

  call test_begin('finite_periodic_normalization')
  sim = sim_config()
  sim%field_solver = 'fmm'
  sim%field_bc_mode = 'periodic2'
  sim%field_periodic_far_correction = 'none'
  call normalize_legacy_physics_config(sim, field, periodic2, panel)
  call assert_true(trim(periodic2%nonzero_mode_backend) == 'legacy_finite_images', 'finite backend mismatch')
  call assert_true(trim(periodic2%zero_mode_policy) == 'legacy_not_decomposed', 'finite zero mode mismatch')
  call test_end()

  call test_begin('cached_kneq0_contract')
  sim%field_periodic_far_correction = 'cached_kneq0'
  sim%use_box = .true.
  sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
  sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
  sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  sim%bc_high = [bc_periodic, bc_periodic, bc_open]
  call normalize_legacy_physics_config(sim, field, periodic2, panel)
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_ok, 'cached kneq0 config should be valid')
  periodic2%lower_boundary_model = 'symmetric_vacuum'
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_ok, 'symmetric vacuum should be valid')
  periodic2%lower_boundary_model = 'unknown'
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'unknown lower model must fail closed')
  call test_end()

  call test_begin('panel_spectral_reference_contract')
  sim%field_solver = 'direct'
  field%backend = 'direct'
  periodic2%nonzero_mode_backend = 'panel_spectral_reference'
  periodic2%zero_mode_policy = 'exclude_k0'
  periodic2%lower_boundary_model = 'e_bottom_zero'
  periodic2%max_nonzero_mode_potential_step = 0.0_dp
  panel%kernel_id = 'triangle_p0_exact_direct'
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_ok, 'panel spectral reference should remain available')
  periodic2%max_nonzero_mode_potential_step = 1.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'panel spectral adaptive limit must fail closed')
  call test_end()

  call test_begin('adaptive_limit_must_be_finite')
  sim%field_solver = 'fmm'
  field%backend = 'fmm'
  periodic2%nonzero_mode_backend = 'cached_kneq0'
  panel%kernel_id = 'triangle_p0_exact_p2m_near'
  periodic2%lower_boundary_model = 'e_bottom_zero'
  periodic2%max_nonzero_mode_potential_step = ieee_value(0.0_dp, ieee_quiet_nan)
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'NaN adaptive limit must fail closed')
  periodic2%max_nonzero_mode_potential_step = -1.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'negative adaptive limit must fail closed')
  call test_end()

  call test_summary()
end program test_physics_config_types

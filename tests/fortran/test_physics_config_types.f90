!> 物理モデル設定の legacy normalization と fail-closed 契約を検証する。
program test_physics_config_types
  use bem_kinds, only: i32
  use bem_types, only: sim_config
  use bem_physics_config_types, only: &
    field_physics_config, periodic2_physics_config, panel_kernel_config, outer_plasma_config, coupling_config, &
    physics_config_ok, physics_config_invalid_combination, physics_config_unavailable, &
    normalize_legacy_physics_config, validate_phase0_physics_config, validate_phase1_panel_config
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32
  implicit none

  type(sim_config) :: sim
  type(field_physics_config) :: field
  type(periodic2_physics_config) :: periodic2
  type(panel_kernel_config) :: panel
  type(outer_plasma_config) :: outer
  type(coupling_config) :: coupling
  integer(i32) :: status
  character(len=256) :: message

  call test_init(7)

  call test_begin('free_legacy_normalization')
  sim = sim_config()
  sim%field_solver = 'treecode'
  sim%field_normalization = 'mesh'
  sim%field_bc_mode = 'free'
  call normalize_legacy_physics_config(sim, field, periodic2, panel, outer, coupling)
  call assert_true(trim(field%backend) == 'treecode', 'free legacy backend mismatch')
  call assert_true(trim(field%normalization) == 'mesh', 'free legacy normalization mismatch')
  call assert_true(trim(panel%source_model) == 'point', 'legacy source model must remain point')
  call assert_true(trim(periodic2%nonzero_mode_backend) == 'not_applicable', 'free nonzero backend mismatch')
  call assert_true(trim(periodic2%zero_mode_policy) == 'not_applicable', 'free zero mode mismatch')
  call validate_phase0_physics_config(field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'free legacy config should be valid')
  call test_end()

  call test_begin('triangle_panel_direct_free_available')
  sim = sim_config()
  sim%field_solver = 'direct'
  sim%field_bc_mode = 'free'
  sim%softening = 0.0d0
  panel = panel_kernel_config()
  panel%source_model = 'triangle_p0'
  panel%kernel_id = 'triangle_p0_exact_direct'
  call validate_phase1_panel_config(sim, panel, status, message)
  call assert_equal_i32(status, physics_config_ok, 'direct free triangle panel should be available')
  call test_end()

  call test_begin('triangle_panel_unsupported_solver_rejected')
  sim%field_solver = 'auto'
  call validate_phase1_panel_config(sim, panel, status, message)
  call assert_equal_i32(status, physics_config_unavailable, 'triangle panel auto solver must be rejected')
  sim%field_solver = 'direct'
  sim%softening = 1.0e-6
  call validate_phase1_panel_config(sim, panel, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'triangle panel softening must be rejected')
  call test_end()

  call test_begin('finite_image_legacy_normalization')
  sim = sim_config()
  sim%field_solver = 'fmm'
  sim%field_bc_mode = 'periodic2'
  sim%field_periodic_far_correction = 'none'
  call normalize_legacy_physics_config(sim, field, periodic2, panel, outer, coupling)
  call assert_true(trim(periodic2%nonzero_mode_backend) == 'legacy_finite_images', 'finite image backend mismatch')
  call assert_true(trim(periodic2%zero_mode_policy) == 'legacy_not_decomposed', 'finite image zero policy mismatch')
  call validate_phase0_physics_config(field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'finite image legacy config should be valid')
  call test_end()

  call test_begin('root_oracle_legacy_normalization')
  sim%field_periodic_far_correction = 'm2l_root_oracle'
  call normalize_legacy_physics_config(sim, field, periodic2, panel, outer, coupling)
  call assert_true(trim(periodic2%nonzero_mode_backend) == 'legacy_root_oracle', 'root oracle backend mismatch')
  call assert_true(trim(periodic2%zero_mode_policy) == 'legacy_charged_walls', 'root oracle zero policy mismatch')
  call validate_phase0_physics_config(field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'root oracle legacy config should be valid')
  call test_end()

  call test_begin('mismatched_zero_policy_rejected')
  periodic2%zero_mode_policy = 'exclude_k0'
  call validate_phase0_physics_config(field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'mismatched zero policy must be rejected')
  call assert_true(len_trim(message) > 0, 'invalid combination should report a message')
  call test_end()

  call test_begin('future_panel_mode_unavailable')
  sim = sim_config()
  call normalize_legacy_physics_config(sim, field, periodic2, panel, outer, coupling)
  panel%source_model = 'triangle_p0'
  call validate_phase0_physics_config(field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_unavailable, 'triangle panel must remain gated in Phase 0')
  call test_end()

  call test_summary()
end program test_physics_config_types

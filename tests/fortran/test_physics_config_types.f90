!> 物理モデル設定の legacy normalization と fail-closed 契約を検証する。
program test_physics_config_types
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_config, bc_open, bc_periodic
  use bem_physics_config_types, only: &
    field_physics_config, periodic2_physics_config, panel_kernel_config, outer_plasma_config, coupling_config, &
    physics_config_ok, physics_config_invalid_combination, physics_config_unavailable, &
    normalize_legacy_physics_config, validate_phase0_physics_config, validate_phase1_panel_config, &
    validate_active_physics_config
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

  call test_init(14)

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

  call test_begin('cached_kinetic_outer_contract')
  sim = sim_config()
  sim%field_solver = 'fmm'
  sim%field_bc_mode = 'periodic2'
  sim%field_periodic_far_correction = 'cached_kneq0'
  sim%use_box = .true.
  sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
  sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
  sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  sim%bc_high = [bc_periodic, bc_periodic, bc_open]
  call normalize_legacy_physics_config(sim, field, periodic2, panel, outer, coupling)
  outer%model = 'kinetic_1d'
  outer%photoelectron_closure = 'kinetic_mean'
  outer%interface_z = 1.0_dp
  outer%debye_length = 0.2_dp
  outer%thermal_voltage = 2.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'cached kinetic_1d config should be valid')
  coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_unavailable, 'cached kinetic return must fail closed')
  call test_end()

  call test_begin('cached_kneq0_active_contract')
  sim = sim_config()
  sim%field_solver = 'fmm'
  sim%field_bc_mode = 'periodic2'
  sim%field_periodic_far_correction = 'cached_kneq0'
  sim%use_box = .true.
  sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
  sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
  sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  sim%bc_high = [bc_periodic, bc_periodic, bc_open]
  call normalize_legacy_physics_config(sim, field, periodic2, panel, outer, coupling)
  call assert_true(trim(periodic2%nonzero_mode_backend) == 'cached_kneq0', 'cached backend normalization mismatch')
  call assert_true(trim(periodic2%zero_mode_policy) == 'exclude_k0', 'cached zero policy normalization mismatch')
  call assert_true(trim(periodic2%lower_boundary_model) == 'e_bottom_zero', 'cached lower boundary mismatch')
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'cached_kneq0 active config should be valid')
  periodic2%zero_mode_policy = 'legacy_not_decomposed'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'cached_kneq0 must reject an unsplit zero mode')
  call test_end()

  call test_begin('instant_return_active_config')
  sim = sim_config()
  sim%field_solver = 'direct'
  sim%field_bc_mode = 'periodic2'
  sim%softening = 0.0_dp
  sim%use_box = .true.
  sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
  sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
  sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  sim%bc_high = [bc_periodic, bc_periodic, bc_open]
  field = field_physics_config(backend='direct')
  periodic2 = periodic2_physics_config( &
              nonzero_mode_backend='panel_spectral_reference', zero_mode_policy='exclude_k0', &
              lower_boundary_model='e_bottom_zero' &
              )
  panel = panel_kernel_config(source_model='triangle_p0', kernel_id='triangle_p0_exact_direct')
  outer = outer_plasma_config( &
          model='linear_debye', return_model='electrostatic_1d_instant_return', interface_z=1.0_dp, &
          debye_length=0.2_dp, thermal_voltage=10.0_dp &
          )
  coupling = coupling_config( &
             particle_transfer_mode='electrostatic_1d_instant_return', field_evolution_timescale=1.0_dp &
             )
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'instant return config should be valid')
  call test_end()

  call test_begin('individual_photoelectron_closure_contract')
  outer%photoelectron_closure = 'individual_return'
  outer%photoelectron_histogram_bins = 8_i32
  outer%photoelectron_histogram_energy_max = 12.0_dp
  outer%photoelectron_ambient_charge_scale = 4.0_dp
  outer%max_photoelectron_charge_ratio = 0.2_dp
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'individual photoelectron closure should be valid')
  outer%photoelectron_ambient_charge_scale = 0.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'missing photoelectron scale must be rejected')
  outer%photoelectron_ambient_charge_scale = 4.0_dp
  call test_end()

  call test_begin('statistical_photoelectron_closure_unavailable')
  outer%photoelectron_closure = 'statistical_return'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_unavailable, 'statistical photoelectron closure must fail closed')
  outer%photoelectron_closure = 'none'
  call test_end()

  call test_begin('persistent_outer_queue_unavailable')
  coupling%outer_queue_enabled = .true.
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_unavailable, 'persistent outer queue must fail closed')
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

  call test_begin('triangle_panel_fmm_available')
  sim = sim_config()
  sim%field_solver = 'fmm'
  sim%field_bc_mode = 'free'
  sim%softening = 0.0_dp
  panel = panel_kernel_config(source_model='triangle_p0', kernel_id='triangle_p0_exact_p2m_near')
  call validate_phase1_panel_config(sim, panel, status, message)
  call assert_equal_i32(status, physics_config_ok, 'triangle panel FMM should be available')
  sim%field_periodic_far_correction = 'm2l_root_oracle'
  call validate_phase1_panel_config(sim, panel, status, message)
  call assert_equal_i32(status, physics_config_unavailable, 'point-source root oracle must be rejected for panels')
  call test_end()

  call test_begin('triangle_panel_auto_and_unsupported_solver')
  sim%field_solver = 'auto'
  panel%kernel_id = 'triangle_p0_exact_auto'
  call validate_phase1_panel_config(sim, panel, status, message)
  call assert_equal_i32(status, physics_config_ok, 'triangle panel auto solver should be available')
  sim%field_solver = 'treecode'
  call validate_phase1_panel_config(sim, panel, status, message)
  call assert_equal_i32(status, physics_config_unavailable, 'triangle panel treecode must be rejected')
  sim%field_solver = 'direct'
  panel%kernel_id = 'triangle_p0_exact_direct'
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

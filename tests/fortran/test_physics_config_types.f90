!> 物理モデル設定の legacy normalization と fail-closed 契約を検証する。
program test_physics_config_types
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
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

  call test_init(15)

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
  outer%photoelectron_density_model = 'kinetic_mean'
  call validate_phase0_physics_config(field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'phase0 must reject kinetic_mean without kinetic_1d')
  outer%photoelectron_density_model = 'none'
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
  outer%photoelectron_density_model = 'kinetic_mean'
  outer%interface_z = 1.0_dp
  outer%debye_length = 0.2_dp
  outer%thermal_voltage = 2.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'cached kinetic_1d config should be valid')

  outer%model = 'unified_linear_response'
  outer%photoelectron_density_model = 'none'
  sim%softening = 0.0_dp
  panel = panel_kernel_config(source_model='triangle_p0', kernel_id='triangle_p0_exact_p2m_near')
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'cached unified linear response config should be valid')

  outer%unified_grid_points = 16_i32
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_true(status /= physics_config_ok, 'unified grid must reject fewer than 17 points')
  outer%unified_grid_points = 129_i32

  outer%accessible_fraction_tolerance = 0.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_true(status /= physics_config_ok, 'unified accessibility tolerance must be positive')
  outer%accessible_fraction_tolerance = 0.1_dp

  sim%e0(3) = 1.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_true(status /= physics_config_ok, 'unified linear response must reject a prescribed uniform field')
  sim%e0 = 0.0_dp

  coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
  outer%return_model = 'electrostatic_1d_instant_return'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_true(status /= physics_config_ok, 'unified linear response must reject scalar instant return')

  coupling%particle_transfer_mode = 'electrostatic_3d_explicit_orbit'
  coupling%field_evolution_timescale = 1.0_dp
  coupling%outer_orbit_dt = 1.0e-3_dp
  coupling%outer_orbit_max_steps = 1000_i32
  coupling%outer_orbit_energy_tolerance = 1.0e-4_dp
  outer%return_model = 'electrostatic_3d_explicit_orbit'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'unified explicit 3D outer orbit should be valid')
  sim%b0(3) = 1.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_unavailable, 'explicit 3D outer orbit must reject nonzero b0')
  sim%b0 = 0.0_dp
  coupling%particle_transfer_mode = 'none'
  outer%return_model = 'none'
  outer%model = 'kinetic_1d'
  outer%photoelectron_density_model = 'kinetic_mean'
  coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
  coupling%field_evolution_timescale = 1.0_dp
  outer%return_model = 'kinetic_1d_profile_return'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'cached kinetic profile return should be valid')
  outer%return_model = 'unsupported_return'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_true(status /= physics_config_ok, 'kinetic outer model must reject an unsupported return identifier')
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
  periodic2%lower_boundary_model = 'symmetric_vacuum'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'cached_kneq0 symmetric vacuum should be valid')
  periodic2%lower_boundary_model = 'unknown'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'cached_kneq0 must reject an unknown lower model')
  periodic2%lower_boundary_model = 'e_bottom_zero'
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
          model='kinetic_1d', return_model='kinetic_1d_profile_return', interface_z=1.0_dp, &
          debye_length=0.2_dp, thermal_voltage=10.0_dp &
          )
  coupling = coupling_config( &
             particle_transfer_mode='electrostatic_1d_instant_return', field_evolution_timescale=1.0_dp &
             )
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'instant return config should be valid')
  sim%open_boundary_model = 'potential_barrier'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'ordinary open barrier must remain valid outside z-high ownership')
  sim%reservoir_potential_model = 'infinity_barrier'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32( &
    status, physics_config_invalid_combination, 'kinetic transfer must reject a competing scalar inflow owner' &
    )
  sim%reservoir_potential_model = 'none'
  sim%open_boundary_model = 'escape'
  call test_end()

  call test_begin('unknown_photoelectron_density_model')
  outer%photoelectron_density_model = 'statistical_return'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'unknown photoelectron density model must fail closed')
  outer%photoelectron_density_model = 'none'
  call test_end()

  call test_begin('zhao_charge_driven_closure_contract')
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
  outer%kinetic_closure = 'zhao_charge_driven'
  outer%zhao_branch = 'c'
  outer%photoelectron_density_model = 'none'
  outer%interface_z = 1.0_dp
  outer%debye_length = 0.2_dp
  outer%thermal_voltage = 2.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'Zhao charge-driven closure should be valid')

  outer%model = 'unified_linear_response'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_true(status /= physics_config_ok, 'Zhao closure must require kinetic_1d')
  outer%model = 'kinetic_1d'
  outer%photoelectron_density_model = 'kinetic_mean'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_true(status /= physics_config_ok, 'Zhao closure must own its photoelectron population')
  outer%photoelectron_density_model = 'none'
  sim%reservoir_potential_model = 'infinity_barrier'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_true(status /= physics_config_ok, 'Zhao closure must reject a reservoir-potential correction')
  sim%reservoir_potential_model = 'none'
  sim%sheath_photoelectron_ref_density_cm3 = 0.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_true(status /= physics_config_ok, 'Zhao closure must require positive photoelectron reference density')
  outer%photoelectron_source_scale = 0.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'no-photo Zhao must not require a photoelectron reference density')
  outer%photoelectron_source_scale = -1.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_true(status /= physics_config_ok, 'negative photoelectron source scale must fail closed')
  outer%photoelectron_source_scale = 1.0_dp
  sim%sheath_photoelectron_ref_density_cm3 = 64.0_dp
  outer%zhao_branch = 'unknown'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_true(status /= physics_config_ok, 'unknown Zhao branch must fail closed')
  outer%zhao_branch = 'auto'
  outer%kinetic_closure = 'unknown'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_true(status /= physics_config_ok, 'unknown kinetic closure must fail closed')
  outer%kinetic_closure = 'absorbing_maxwellian'
  outer%zhao_branch = 'c'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_true(status /= physics_config_ok, 'explicit Zhao branch must require the Zhao closure')
  outer%kinetic_closure = 'zhao_charge_driven'
  outer%zhao_branch = 'auto'
  outer%return_model = 'kinetic_1d_profile_return'
  coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
  coupling%field_evolution_timescale = 1.0_dp
  call test_end()

  call test_begin('zhao_floating_steady_start_contract')
  coupling%steady_start_mode = 'zhao_floating'
  coupling%steady_start_mesh_id = 2_i32
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'Zhao floating steady start should be valid')
  coupling%steady_start_mesh_id = 0_i32
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'steady start must require a positive mesh ID')
  coupling%steady_start_mesh_id = 2_i32
  coupling%outer_queue_enabled = .true.
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'steady start must reject the transient queue')
  coupling%outer_queue_enabled = .false.
  outer%return_model = 'none'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'steady start must require kinetic profile return')
  outer%return_model = 'kinetic_1d_profile_return'
  periodic2%zero_mode_policy = 'legacy_not_decomposed'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'steady start must require the split zero mode')
  periodic2%zero_mode_policy = 'exclude_k0'
  coupling%steady_start_mode = 'unknown'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'unknown steady-start mode must fail closed')
  coupling%steady_start_mode = 'none'
  call test_end()

  call test_begin('persistent_outer_queue_contract')
  coupling%outer_queue_enabled = .true.
  outer%photoelectron_source_scale = 0.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_true(status /= physics_config_ok, 'persistent queue must require a photoelectron source')
  outer%photoelectron_source_scale = 1.0_dp
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'outer queue must require physical batch time')
  sim%batch_duration = 1.0e-6_dp
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_ok, 'bounded Zhao outer queue should be available')
  coupling%field_evolution_timescale = ieee_value(0.0_dp, ieee_quiet_nan)
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, &
                        'outer queue must reject a non-finite field-evolution timescale')
  coupling%field_evolution_timescale = 1.0_dp
  outer%zhao_branch = 'b'
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, &
                        'outer queue must require automatic branch continuation')
  outer%zhao_branch = 'auto'
  sim%batch_duration = 0.2_dp
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, &
                        'outer queue batch time must resolve field evolution')
  sim%batch_duration = 1.0e-6_dp
  coupling%outer_update_stride = 2_i32
  call validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'outer queue must refresh the closure every batch')
  coupling%outer_update_stride = 1_i32
  coupling%outer_queue_enabled = .false.
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
  call test_end()

  call test_begin('triangle_panel_auto_and_treecode')
  sim%field_solver = 'auto'
  panel%kernel_id = 'triangle_p0_exact_auto'
  call validate_phase1_panel_config(sim, panel, status, message)
  call assert_equal_i32(status, physics_config_ok, 'triangle panel auto solver should be available')
  sim%field_solver = 'treecode'
  panel%kernel_id = 'triangle_p0_exact_tree_near'
  call validate_phase1_panel_config(sim, panel, status, message)
  call assert_equal_i32(status, physics_config_ok, 'triangle panel treecode should be available')
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

  call test_begin('removed_root_oracle_normalization')
  sim%field_periodic_far_correction = 'm2l_root_oracle'
  call normalize_legacy_physics_config(sim, field, periodic2, panel, outer, coupling)
  call assert_true(trim(periodic2%nonzero_mode_backend) == 'invalid', 'removed root oracle must normalize to invalid')
  call assert_true(trim(periodic2%zero_mode_policy) == 'invalid', 'removed root oracle zero policy must be invalid')
  call validate_phase0_physics_config(field, periodic2, panel, outer, coupling, status, message)
  call assert_equal_i32(status, physics_config_invalid_combination, 'removed root oracle config must be invalid')
  call test_end()

  call test_begin('mismatched_zero_policy_rejected')
  periodic2%nonzero_mode_backend = 'cached_kneq0'
  periodic2%zero_mode_policy = 'legacy_not_decomposed'
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

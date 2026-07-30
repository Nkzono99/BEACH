!> model、ordered mesh、ordered species の fingerprint 感度と決定性を検証する。
program test_model_fingerprint
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type
  use bem_mesh, only: init_mesh
  use bem_app_config, only: app_config, default_app_config, species_from_defaults
  use bem_physics_config_types, only: normalize_legacy_physics_config
  use bem_model_fingerprint, only: model_fingerprint, mesh_fingerprint, species_fingerprint
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32
  implicit none

  type(mesh_type) :: mesh, mesh_changed
  type(app_config) :: cfg, cfg_changed
  real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2)
  character(len=16) :: fp_a, fp_b

  v0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
  v1(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
  v2(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
  v0(:, 2) = [0.0_dp, 0.0_dp, 1.0_dp]
  v1(:, 2) = [1.0_dp, 0.0_dp, 1.0_dp]
  v2(:, 2) = [0.0_dp, 1.0_dp, 1.0_dp]
  call init_mesh(mesh, v0, v1, v2)
  call default_app_config(cfg)
  cfg%n_particle_species = 2_i32
  cfg%particle_species(1) = species_from_defaults()
  cfg%particle_species(1)%species_key = 'electron'
  cfg%particle_species(1)%q_particle = -1.0_dp
  cfg%particle_species(1)%m_particle = 2.0_dp
  cfg%particle_species(1)%w_particle = 3.0_dp
  cfg%particle_species(2) = species_from_defaults()
  cfg%particle_species(2)%species_key = 'ion'
  cfg%particle_species(2)%q_particle = 1.0_dp
  cfg%particle_species(2)%m_particle = 4.0_dp
  cfg%particle_species(2)%w_particle = 5.0_dp

  call test_init(19)

  call test_begin('deterministic_fingerprint')
  fp_a = mesh_fingerprint(mesh)
  fp_b = mesh_fingerprint(mesh)
  call assert_true(fp_a == fp_b, 'mesh fingerprint must be deterministic')
  call assert_equal_i32(int(len_trim(fp_a), i32), 16_i32, 'fingerprint length mismatch')
  call test_end()

  call test_begin('checkpoint_v4_default_model_fingerprint')
  fp_a = model_fingerprint(cfg)
  call assert_true( &
    fp_a == '080673E71F479214', &
    'checkpoint-v4 default model fingerprint stream changed: got '//fp_a &
    )
  call test_end()

  call test_begin('steady_start_mode_change_detected')
  cfg_changed = cfg
  cfg_changed%coupling%steady_start_mode = 'zhao_floating'
  call assert_true( &
    model_fingerprint(cfg_changed) /= model_fingerprint(cfg), &
    'steady-start mode must alter fingerprint' &
    )
  call test_end()

  call test_begin('inactive_steady_start_mesh_id_ignored')
  cfg_changed = cfg
  cfg_changed%coupling%steady_start_mesh_id = 2_i32
  call assert_true( &
    model_fingerprint(cfg_changed) == model_fingerprint(cfg), &
    'inactive steady-start mesh ID must preserve the legacy fingerprint stream' &
    )
  call test_end()

  call test_begin('active_steady_start_mesh_id_change_detected')
  cfg_changed = cfg
  cfg_changed%coupling%steady_start_mode = 'zhao_floating'
  fp_a = model_fingerprint(cfg_changed)
  cfg_changed%coupling%steady_start_mesh_id = 2_i32
  call assert_true( &
    model_fingerprint(cfg_changed) /= fp_a, &
    'active steady-start mesh ID must alter fingerprint' &
    )
  call test_end()

  call test_begin('split_numeric_contract_change_detected')
  cfg_changed = cfg
  cfg_changed%periodic2%interface_sample_n = cfg%periodic2%interface_sample_n + 1_i32
  call assert_true( &
    model_fingerprint(cfg_changed) /= model_fingerprint(cfg), &
    'split periodic numeric contract must alter fingerprint' &
    )
  call test_end()

  call test_begin('photoelectron_density_model_change_detected')
  cfg_changed = cfg
  cfg_changed%outer_plasma%photoelectron_density_model = 'kinetic_mean'
  call assert_true( &
    model_fingerprint(cfg_changed) /= model_fingerprint(cfg), &
    'photoelectron density model must alter fingerprint' &
    )
  call test_end()

  call test_begin('kinetic_closure_change_detected')
  cfg_changed = cfg
  cfg_changed%outer_plasma%kinetic_closure = 'zhao_charge_driven'
  call assert_true( &
    model_fingerprint(cfg_changed) /= model_fingerprint(cfg), &
    'kinetic closure must alter fingerprint' &
    )
  call test_end()

  call test_begin('zhao_branch_change_detected')
  cfg_changed = cfg
  cfg_changed%outer_plasma%zhao_branch = 'c'
  call assert_true( &
    model_fingerprint(cfg_changed) /= model_fingerprint(cfg), &
    'Zhao branch must alter fingerprint' &
    )
  call test_end()

  call test_begin('photoelectron_source_scale_change_detected')
  cfg_changed = cfg
  cfg_changed%outer_plasma%photoelectron_source_scale = 0.0_dp
  call assert_true( &
    model_fingerprint(cfg_changed) /= model_fingerprint(cfg), &
    'photoelectron source scale must alter fingerprint' &
    )
  call test_end()

  call test_begin('mesh_vacuum_side_change_detected')
  mesh_changed = mesh
  mesh_changed%elem_vacuum_sign(1) = 1_i32
  call assert_true(mesh_fingerprint(mesh_changed) /= mesh_fingerprint(mesh), 'vacuum side must alter mesh fingerprint')
  call test_end()

  call test_begin('mesh_vertex_change_detected')
  mesh_changed = mesh
  mesh_changed%v2(3, 2) = mesh_changed%v2(3, 2) + 1.0e-9_dp
  call assert_true(mesh_fingerprint(mesh_changed) /= mesh_fingerprint(mesh), 'mesh vertex change must alter fingerprint')
  call test_end()

  call test_begin('species_contract_change_detected')
  cfg_changed = cfg
  cfg_changed%particle_species(2)%q_particle = 2.0_dp
  call assert_true( &
    species_fingerprint(cfg_changed) /= species_fingerprint(cfg), &
    'species contract change must alter fingerprint' &
    )
  call test_end()

  call test_begin('species_z_high_boundary_change_detected')
  cfg_changed = cfg
  cfg_changed%particle_species(1)%z_high_boundary = 'reflect'
  call assert_true( &
    species_fingerprint(cfg_changed) /= species_fingerprint(cfg), &
    'species z-high boundary override must alter fingerprint' &
    )
  call test_end()

  call test_begin('species_surface_charge_closure_change_detected')
  cfg_changed = cfg
  cfg_changed%particle_species(1)%surface_charge_closure = 'neutral_return'
  call assert_true( &
    species_fingerprint(cfg_changed) /= species_fingerprint(cfg), &
    'species surface-charge closure must alter fingerprint' &
    )
  call test_end()

  call test_begin('species_order_change_detected')
  cfg_changed = cfg
  cfg_changed%particle_species(1) = cfg%particle_species(2)
  cfg_changed%particle_species(2) = cfg%particle_species(1)
  call assert_true(species_fingerprint(cfg_changed) /= species_fingerprint(cfg), 'species order must alter fingerprint')
  call test_end()

  call test_begin('species_presence_flag_change_detected')
  cfg_changed = cfg
  cfg_changed%particle_species(1)%has_temperature_ev = .true.
  call assert_true( &
    species_fingerprint(cfg_changed) /= species_fingerprint(cfg), &
    'species presence flags must alter fingerprint' &
    )
  call test_end()

  call test_begin('model_backend_change_detected')
  cfg_changed = cfg
  cfg_changed%sim%field_solver = 'fmm'
  cfg_changed%sim%field_bc_mode = 'periodic2'
  cfg_changed%sim%field_periodic_far_correction = 'none'
  call normalize_legacy_physics_config( &
    cfg_changed%sim, cfg_changed%field, cfg_changed%periodic2, cfg_changed%panel, &
    cfg_changed%outer_plasma, cfg_changed%coupling &
    )
  call assert_true(model_fingerprint(cfg_changed) /= model_fingerprint(cfg), 'model backend change must alter fingerprint')
  call test_end()

  call test_begin('periodic_generation_tolerance_change_detected')
  cfg_changed = cfg
  cfg_changed%sim%field_periodic_generation_tolerance = &
    2.0_dp*cfg%sim%field_periodic_generation_tolerance
  call assert_true( &
    model_fingerprint(cfg_changed) /= model_fingerprint(cfg), &
    'periodic generation tolerance must alter the model fingerprint' &
    )
  call test_end()

  call test_summary()
end program test_model_fingerprint

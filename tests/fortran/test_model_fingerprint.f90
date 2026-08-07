!> model、ordered mesh、ordered speciesのfingerprint感度と決定性を検証する。
program test_model_fingerprint
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, bc_reflect, bc_redistributed_reflect
  use bem_mesh, only: init_mesh
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, particle_inflow_reservoir
  use bem_physics_config_types, only: normalize_legacy_physics_config
  use bem_model_fingerprint, only: model_fingerprint, mesh_fingerprint, species_fingerprint
  use bem_injection_velocity_grid, only: reset_velocity_grid_snapshot_cache
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32, &
                          delete_file_if_exists
  implicit none

  type(mesh_type) :: mesh, mesh_changed
  type(app_config) :: cfg, cfg_changed
  real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2)
  character(len=16) :: fp_a, fp_b
  character(len=*), parameter :: velocity_grid_path = 'test_model_fingerprint_velocity_grid.csv'
  integer :: u, ios

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

  call test_init(18)

  call test_begin('deterministic_fingerprint')
  fp_a = mesh_fingerprint(mesh)
  fp_b = mesh_fingerprint(mesh)
  call assert_true(fp_a == fp_b, 'mesh fingerprint must be deterministic')
  call assert_equal_i32(int(len_trim(fp_a), i32), 16_i32, 'fingerprint length mismatch')
  call test_end()

  call test_begin('mesh_vacuum_side_change_detected')
  mesh_changed = mesh
  mesh_changed%elem_vacuum_sign(1) = 1_i32
  call assert_true(mesh_fingerprint(mesh_changed) /= mesh_fingerprint(mesh), 'vacuum side must alter fingerprint')
  call test_end()

  call test_begin('mesh_vertex_change_detected')
  mesh_changed = mesh
  mesh_changed%v2(3, 2) = mesh_changed%v2(3, 2) + 1.0e-9_dp
  call assert_true(mesh_fingerprint(mesh_changed) /= mesh_fingerprint(mesh), 'vertex change must alter fingerprint')
  call test_end()

  call test_begin('species_contract_change_detected')
  cfg_changed = cfg
  cfg_changed%particle_species(2)%q_particle = 2.0_dp
  call assert_true(species_fingerprint(cfg_changed) /= species_fingerprint(cfg), 'charge change must alter fingerprint')
  call test_end()

  call test_begin('species_particle_boundary_change_detected')
  cfg_changed = cfg
  cfg_changed%particle_species(1)%boundary_low(1) = bc_reflect
  call assert_true(species_fingerprint(cfg_changed) /= species_fingerprint(cfg), 'particle boundary must alter fingerprint')
  call test_end()

  call test_begin('species_boundary_inflow_change_detected')
  cfg_changed = cfg
  cfg_changed%particle_species(1)%boundary_inflow_high(3) = particle_inflow_reservoir
  call assert_true(species_fingerprint(cfg_changed) /= species_fingerprint(cfg), 'boundary inflow must alter fingerprint')
  call test_end()

  call test_begin('plane_source_normal_change_detected')
  cfg_changed = cfg
  cfg_changed%particle_species(1)%source_mode = 'plane_source'
  cfg_changed%particle_species(1)%has_source_normal = .true.
  cfg_changed%particle_species(1)%source_normal = [1.0_dp, 0.0_dp, 0.0_dp]
  fp_a = species_fingerprint(cfg_changed)
  cfg_changed%particle_species(1)%source_normal = [-1.0_dp, 0.0_dp, 0.0_dp]
  fp_b = species_fingerprint(cfg_changed)
  call assert_true(fp_a /= fp_b, 'plane source normal must alter fingerprint')
  call test_end()

  call test_begin('global_particle_boundary_change_detected')
  cfg_changed = cfg
  cfg_changed%particle_boundary_high(3) = bc_reflect
  call assert_true(model_fingerprint(cfg_changed) /= model_fingerprint(cfg), 'global particle boundary must alter fingerprint')
  call test_end()

  call test_begin('redistributed_reflect_seed_change_detected')
  cfg_changed = cfg
  cfg_changed%sim%rng_seed = cfg_changed%sim%rng_seed + 1_i32
  call assert_true( &
    model_fingerprint(cfg_changed) == model_fingerprint(cfg), &
    'ordinary-boundary checkpoints must keep their existing RNG-state compatibility' &
    )
  cfg_changed = cfg
  cfg_changed%particle_species(1)%boundary_high(3) = bc_redistributed_reflect
  fp_a = model_fingerprint(cfg_changed)
  cfg_changed%sim%rng_seed = cfg_changed%sim%rng_seed + 1_i32
  fp_b = model_fingerprint(cfg_changed)
  call assert_true(fp_a /= fp_b, 'redistributed-reflect RNG seed must alter the model fingerprint')
  call test_end()

  call test_begin('surface_charge_closure_change_detected')
  cfg_changed = cfg
  cfg_changed%particle_species(1)%surface_charge_closure = 'neutral_return'
  call assert_true(species_fingerprint(cfg_changed) /= species_fingerprint(cfg), 'closure must alter fingerprint')
  call test_end()

  call test_begin('surface_current_model_change_detected')
  cfg_changed = cfg
  cfg_changed%surface_current%model = 'zhao_stationary'
  cfg_changed%surface_current%solar_elevation_deg = 60.0_dp
  call assert_true(model_fingerprint(cfg_changed) /= model_fingerprint(cfg), 'surface current model must alter fingerprint')
  call test_end()

  call test_begin('species_order_change_detected')
  cfg_changed = cfg
  cfg_changed%particle_species(1) = cfg%particle_species(2)
  cfg_changed%particle_species(2) = cfg%particle_species(1)
  call assert_true(species_fingerprint(cfg_changed) /= species_fingerprint(cfg), 'order must alter fingerprint')
  call test_end()

  call test_begin('species_presence_flag_change_detected')
  cfg_changed = cfg
  cfg_changed%particle_species(1)%has_temperature_ev = .true.
  call assert_true(species_fingerprint(cfg_changed) /= species_fingerprint(cfg), 'presence flag must alter fingerprint')
  call test_end()

  call test_begin('velocity_grid_snapshot_is_immutable_within_run')
  call delete_file_if_exists(velocity_grid_path)
  open (newunit=u, file=velocity_grid_path, status='replace', action='write', iostat=ios)
  call assert_equal_i32(int(ios, i32), 0_i32, 'failed to create velocity-grid fingerprint fixture')
  write (u, '(a)') 'vx,vy,vz,f'
  write (u, '(a)') '1,0,0,1'
  close (u)
  cfg_changed = cfg
  ! Match runtime dispatch: programmatically constructed configs may use mixed case.
  cfg_changed%particle_species(1)%velocity_distribution = 'GRID'
  cfg_changed%particle_species(1)%velocity_grid_path = velocity_grid_path
  fp_a = species_fingerprint(cfg_changed)
  open (newunit=u, file=velocity_grid_path, status='replace', action='write', iostat=ios)
  call assert_equal_i32(int(ios, i32), 0_i32, 'failed to replace velocity-grid fingerprint fixture')
  write (u, '(a)') 'vx,vy,vz,f'
  write (u, '(a)') '1,0,0,2'
  close (u)
  fp_b = species_fingerprint(cfg_changed)
  call assert_true(fp_a == fp_b, 'an active run must retain its initial velocity-grid snapshot')
  call reset_velocity_grid_snapshot_cache()
  fp_b = species_fingerprint(cfg_changed)
  call assert_true(fp_a /= fp_b, 'a new run must fingerprint the replacement velocity-grid contents')
  call reset_velocity_grid_snapshot_cache()
  call delete_file_if_exists(velocity_grid_path)
  call test_end()

  call test_begin('disabled_velocity_grid_does_not_require_file')
  cfg_changed = cfg
  cfg_changed%particle_species(1)%enabled = .false.
  cfg_changed%particle_species(1)%velocity_distribution = 'grid'
  cfg_changed%particle_species(1)%velocity_grid_path = 'missing-disabled-velocity-grid.csv'
  fp_a = species_fingerprint(cfg_changed)
  call assert_equal_i32(int(len_trim(fp_a), i32), 16_i32, 'disabled species fingerprint length mismatch')
  call test_end()

  call test_begin('model_backend_change_detected')
  cfg_changed = cfg
  cfg_changed%sim%field_solver = 'fmm'
  cfg_changed%sim%field_bc_mode = 'periodic2'
  cfg_changed%sim%field_periodic_far_correction = 'none'
  call normalize_legacy_physics_config( &
    cfg_changed%sim, cfg_changed%field, cfg_changed%periodic2, cfg_changed%panel &
    )
  call assert_true(model_fingerprint(cfg_changed) /= model_fingerprint(cfg), 'backend change must alter fingerprint')
  call test_end()

  call test_begin('tree_override_presence_change_detected')
  cfg_changed = cfg
  cfg_changed%sim%has_tree_theta = .true.
  call assert_true( &
    model_fingerprint(cfg_changed) /= model_fingerprint(cfg), &
    'explicit tree theta must differ from an equal-valued automatic setting' &
    )
  cfg_changed = cfg
  cfg_changed%sim%has_tree_leaf_max = .true.
  call assert_true( &
    model_fingerprint(cfg_changed) /= model_fingerprint(cfg), &
    'explicit tree leaf size must differ from an equal-valued automatic setting' &
    )
  call test_end()

  call test_begin('periodic_generation_tolerance_change_detected')
  cfg_changed = cfg
  cfg_changed%sim%field_periodic_generation_tolerance = 2.0_dp*cfg%sim%field_periodic_generation_tolerance
  call assert_true(model_fingerprint(cfg_changed) /= model_fingerprint(cfg), 'tolerance must alter fingerprint')
  call test_end()

  call test_summary()
end program test_model_fingerprint

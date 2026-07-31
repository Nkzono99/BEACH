!> 現行の局所 reservoir / closed PE 設定をFortran parserで検証する。
program test_app_config_parser
  use bem_kinds, only: dp, i32
  use bem_types, only: bc_open, bc_reflect, bc_redistributed_reflect
  use bem_app_config, only: app_config, default_app_config, load_app_config, &
                            particles_per_batch_from_config, total_particles_from_config
  use bem_config_helpers, only: resolve_particle_boundaries
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp
  implicit none

  type(app_config) :: cfg
  integer(i32) :: effective_boundary_low(3), effective_boundary_high(3)

  call test_init(4)

  call test_begin('default_config')
  call default_app_config(cfg)
  call assert_true(trim(cfg%sim%field_solver) == 'auto', 'default field solver mismatch')
  call assert_true(trim(cfg%sim%field_bc_mode) == 'free', 'default field boundary mismatch')
  call assert_true(trim(cfg%sim%reservoir_potential_model) == 'none', 'default inflow model mismatch')
  call assert_true(trim(cfg%sim%open_boundary_model) == 'escape', 'default open model mismatch')
  call test_end()

  call test_begin('tutorial_config')
  call default_app_config(cfg)
  call load_app_config('examples/tutorial_insulator.toml', cfg)
  call assert_true(trim(cfg%sim%field_solver) == 'fmm', 'tutorial field solver mismatch')
  call assert_true(trim(cfg%sim%field_bc_mode) == 'periodic2', 'tutorial field boundary mismatch')
  call assert_equal_i32(cfg%n_particle_species, 1_i32, 'tutorial species count mismatch')
  call assert_equal_i32(particles_per_batch_from_config(cfg), 1_i32, 'tutorial batch particle count mismatch')
  call assert_equal_i32(total_particles_from_config(cfg), cfg%sim%batch_count, 'tutorial total particle count mismatch')
  call test_end()

  call test_begin('closed_photoelectron_config')
  call default_app_config(cfg)
  call load_app_config('examples/periodic2_closed_photoelectron.toml', cfg)
  call assert_equal_i32(cfg%n_particle_species, 3_i32, 'closed case species count mismatch')
  call assert_true(trim(cfg%sim%reservoir_potential_model) == 'none', 'closed case must use source VDF inflow')
  call assert_true(trim(cfg%sim%open_boundary_model) == 'escape', 'closed case must use ordinary escape')
  call assert_true(trim(cfg%particle_species(3)%source_mode) == 'photo_raycast', 'photoelectron source mode mismatch')
  call assert_equal_i32(cfg%particle_species(3)%boundary_high(3), bc_reflect, 'photoelectron top boundary mismatch')
  call assert_true( &
    trim(cfg%particle_species(3)%surface_charge_closure) == 'neutral_return', &
    'photoelectron closure mismatch' &
    )
  call assert_close_dp(cfg%particle_species(3)%temperature_ev, 1.5_dp, 1.0e-15_dp, 'photoelectron temperature mismatch')
  call test_end()

  call test_begin('all_particle_boundary_faces')
  call default_app_config(cfg)
  call load_app_config('tests/fortran/particle_boundary_faces.toml', cfg)
  call assert_true(all(cfg%particle_boundary_low == [bc_open, bc_reflect, bc_open]), 'global low faces mismatch')
  call assert_true( &
    all(cfg%particle_boundary_high == [bc_reflect, bc_open, bc_redistributed_reflect]), &
    'global high faces mismatch' &
    )
  call assert_equal_i32(cfg%particle_species(1)%boundary_low(1), bc_reflect, 'species x-low mismatch')
  call assert_equal_i32(cfg%particle_species(1)%boundary_high(1), bc_open, 'species x-high mismatch')
  call assert_equal_i32(cfg%particle_species(1)%boundary_high(2), bc_reflect, 'species y-high mismatch')
  call assert_equal_i32(cfg%particle_species(1)%boundary_low(3), bc_reflect, 'species z-low mismatch')
  call resolve_particle_boundaries( &
    cfg%sim, cfg%particle_boundary_low, cfg%particle_boundary_high, cfg%particle_species(1), &
    effective_boundary_low, effective_boundary_high &
    )
  call assert_true(all(effective_boundary_low == bc_reflect), 'effective low faces mismatch')
  call assert_true( &
    all(effective_boundary_high == [bc_open, bc_reflect, bc_redistributed_reflect]), &
    'effective high faces mismatch' &
    )
  call test_end()

  call test_summary()
end program test_app_config_parser

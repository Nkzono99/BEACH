!> app_config の設定読込と派生値計算を検証するテスト。
program test_app_config_parser
  use bem_kinds, only: dp, i32
  use bem_types, only: bc_open, bc_reflect, bc_periodic
  use bem_app_config, only: app_config, default_app_config, load_app_config, &
                            particles_per_batch_from_config, total_particles_from_config
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, assert_allclose_1d, delete_file_if_exists
  implicit none

  type(app_config) :: cfg, photo_cfg, large_cfg, periodic_cfg, periodic_oracle_cfg, sheath_cfg, high_level_cfg
  type(app_config) :: multiline_cfg, toml_syntax_cfg, panel_cfg, split_cfg
  character(len=*), parameter :: cfg_path = 'test_app_config_parser_tmp.toml'
  character(len=*), parameter :: photo_cfg_path = 'test_app_config_parser_photo_tmp.toml'
  character(len=*), parameter :: large_cfg_path = 'test_app_config_parser_large_tmp.toml'
  character(len=*), parameter :: periodic_cfg_path = 'test_app_config_parser_periodic_tmp.toml'
  character(len=*), parameter :: periodic_oracle_cfg_path = 'test_app_config_parser_periodic_oracle_tmp.toml'
  character(len=*), parameter :: sheath_cfg_path = 'test_app_config_parser_sheath_tmp.toml'
  character(len=*), parameter :: high_level_cfg_path = 'test_app_config_parser_high_level_tmp.toml'
  character(len=*), parameter :: multiline_cfg_path = 'test_app_config_parser_multiline_tmp.toml'
  character(len=*), parameter :: toml_syntax_cfg_path = 'test_app_config_parser_toml_syntax_tmp.toml'
  character(len=*), parameter :: panel_cfg_path = 'test_app_config_parser_panel_tmp.toml'
  character(len=*), parameter :: split_cfg_path = 'test_app_config_parser_split_tmp.toml'

  call write_config_fixture(cfg_path)
  call write_photo_config_fixture(photo_cfg_path)
  call write_large_config_fixture(large_cfg_path)
  call write_periodic_config_fixture(periodic_cfg_path)
  call write_periodic_oracle_config_fixture(periodic_oracle_cfg_path)
  call write_sheath_config_fixture(sheath_cfg_path)
  call write_high_level_config_fixture(high_level_cfg_path)
  call write_multiline_config_fixture(multiline_cfg_path)
  call write_toml_syntax_config_fixture(toml_syntax_cfg_path)
  call write_panel_config_fixture(panel_cfg_path)
  call write_split_config_fixture(split_cfg_path)

  call test_init(11)

  call test_begin('defaults_and_basic_config')
  call default_app_config(cfg)
  call assert_true(trim(cfg%sim%field_solver) == 'auto', 'default field_solver mismatch')
  call assert_true(trim(cfg%sim%field_normalization) == 'si', 'default field_normalization mismatch')
  call assert_close_dp(cfg%sim%field_length_scale, 1.0d0, 1.0d-15, 'default field_length_scale mismatch')
  call assert_true(trim(cfg%mesh_surface_model) == 'insulator', 'default mesh surface_model mismatch')
  call assert_close_dp(cfg%mesh_epsilon_r, 1.0d0, 1.0d-15, 'default mesh epsilon_r mismatch')
  call assert_true(trim(cfg%sim%field_bc_mode) == 'free', 'default field_bc_mode mismatch')
  call assert_equal_i32(cfg%sim%field_periodic_image_layers, 1_i32, 'default field_periodic_image_layers mismatch')
  call assert_true(trim(cfg%sim%field_periodic_far_correction) == 'none', 'default field_periodic_far_correction mismatch')
  call assert_true(trim(cfg%field%backend) == 'auto', 'default typed field backend mismatch')
  call assert_true(trim(cfg%panel%source_model) == 'point', 'default typed source model mismatch')
  call assert_true( &
    trim(cfg%periodic2%nonzero_mode_backend) == 'not_applicable', &
    'default typed periodic backend mismatch' &
    )
  call assert_true(trim(cfg%outer_plasma%model) == 'none', 'default typed outer model mismatch')
  call assert_close_dp(cfg%sim%field_periodic_ewald_alpha, 0.0d0, 1.0d-15, 'default field_periodic_ewald_alpha mismatch')
  call assert_equal_i32(cfg%sim%field_periodic_ewald_layers, 4_i32, 'default field_periodic_ewald_layers mismatch')
  call assert_true( &
    trim(cfg%sim%field_periodic_cache_dir) == '.beach_cache/periodic2', &
    'default field_periodic_cache_dir mismatch' &
    )
  call assert_close_dp( &
    cfg%sim%field_periodic_generation_tolerance, 1.0e-8_dp, 1.0e-20_dp, &
    'default field_periodic_generation_tolerance mismatch' &
    )
  call assert_close_dp(cfg%sim%tree_theta, 0.5d0, 1.0d-15, 'default tree_theta mismatch')
  call assert_true(.not. cfg%sim%has_tree_theta, 'default has_tree_theta should be false')
  call assert_equal_i32(cfg%sim%tree_leaf_max, 16_i32, 'default tree_leaf_max mismatch')
  call assert_true(.not. cfg%sim%has_tree_leaf_max, 'default has_tree_leaf_max should be false')
  call assert_equal_i32(cfg%sim%tree_min_nelem, 256_i32, 'default tree_min_nelem mismatch')
  call assert_true(.not. cfg%write_mesh_potential, 'default write_mesh_potential should be false')
  call assert_true(len_trim(cfg%output_restart_from) == 0, 'default output_restart_from should be empty')
  call load_app_config(cfg_path, cfg)

  call assert_true(trim(cfg%mesh_mode) == 'template', 'mesh.mode was not parsed')
  call assert_true(trim(cfg%mesh_surface_model) == 'dielectric', 'mesh surface_model was not parsed')
  call assert_close_dp(cfg%mesh_epsilon_r, 2.5d0, 1.0d-15, 'mesh epsilon_r was not parsed')
  call assert_true(trim(cfg%templates(2)%kind) == 'sphere', 'second template kind mismatch')
  call assert_true(trim(cfg%templates(3)%kind) == 'annulus', 'third template kind mismatch')
  call assert_close_dp(cfg%templates(3)%radius, 0.3d0, 1.0d-12, 'annulus radius mismatch')
  call assert_close_dp(cfg%templates(3)%inner_radius, 0.1d0, 1.0d-12, 'annulus inner_radius mismatch')
  call assert_equal_i32(cfg%templates(3)%n_theta, 16_i32, 'annulus n_theta mismatch')
  call assert_equal_i32(cfg%templates(3)%n_r, 3_i32, 'annulus n_r mismatch')
  call assert_true(trim(cfg%templates(4)%kind) == 'cylinder', 'fourth template kind mismatch')
  call assert_true(cfg%templates(4)%has_cap_top, 'cylinder has_cap_top mismatch')
  call assert_true(cfg%templates(4)%has_cap_bottom, 'cylinder has_cap_bottom mismatch')
  call assert_true(cfg%templates(4)%cap_top, 'cylinder cap_top mismatch')
  call assert_true(.not. cfg%templates(4)%cap_bottom, 'cylinder cap_bottom mismatch')
  call assert_true(trim(cfg%templates(5)%kind) == 'plate_hole', 'fifth template kind mismatch')
  call assert_close_dp(cfg%templates(5)%size_x, 1.5d0, 1.0d-12, 'plate_hole size_x mismatch')
  call assert_close_dp(cfg%templates(5)%size_y, 0.8d0, 1.0d-12, 'plate_hole size_y mismatch')
  call assert_close_dp(cfg%templates(5)%radius, 0.2d0, 1.0d-12, 'plate_hole radius mismatch')
  call assert_equal_i32(cfg%templates(5)%n_theta, 20_i32, 'plate_hole n_theta mismatch')
  call assert_equal_i32(cfg%templates(5)%n_r, 2_i32, 'plate_hole n_r mismatch')
  call assert_equal_i32(cfg%n_particle_species, 2_i32, 'n_particle_species mismatch')
  call assert_true(trim(cfg%particle_species(1)%species_key) == 'electron', 'explicit species_key mismatch')
  call assert_true(trim(cfg%particle_species(2)%species_key) == 'species_2', 'generated species_key mismatch')
  call assert_equal_i32(particles_per_batch_from_config(cfg), 5_i32, 'per-batch particle count mismatch')
  call assert_equal_i32(total_particles_from_config(cfg), 15_i32, 'total particle count mismatch')
  call assert_equal_i32(cfg%n_particles, 15_i32, 'cached n_particles mismatch')
  call assert_equal_i32(cfg%sim%bc_low(1), bc_periodic, 'bc_x_low mismatch')
  call assert_equal_i32(cfg%sim%bc_high(2), bc_reflect, 'bc_y_high mismatch')
  call assert_equal_i32(cfg%sim%bc_low(3), bc_open, 'bc_z_low mismatch')
  call assert_true(trim(cfg%sim%reservoir_potential_model) == 'infinity_barrier', 'reservoir_potential_model mismatch')
  call assert_close_dp(cfg%sim%phi_infty, -2.0d0, 1.0d-12, 'phi_infty mismatch')
  call assert_true(trim(cfg%sim%open_boundary_model) == 'potential_barrier', 'open_boundary_model mismatch')
  call assert_equal_i32(cfg%sim%injection_face_phi_grid_n, 5_i32, 'injection_face_phi_grid_n mismatch')
  call assert_equal_i32(cfg%history_stride, 2_i32, 'history_stride mismatch')
  call assert_true(cfg%write_mesh_potential, 'write_mesh_potential mismatch')
  call assert_true(trim(cfg%output_restart_from) == 'outputs/parent', 'output.restart_from mismatch')
  call assert_close_dp(cfg%sim%dt, 2.5d-9, 1.0d-15, 'dt mismatch')
  call assert_true(trim(cfg%sim%field_solver) == 'fmm', 'field_solver mismatch')
  call assert_true(trim(cfg%sim%field_normalization) == 'length', 'field_normalization mismatch')
  call assert_close_dp(cfg%sim%field_length_scale, 2.0d-6, 1.0d-18, 'field_length_scale mismatch')
  call assert_true(trim(cfg%sim%field_bc_mode) == 'free', 'field_bc_mode mismatch')
  call assert_true(trim(cfg%templates(1)%surface_model) == 'insulator', 'default surface_model mismatch')
  call assert_true(trim(cfg%templates(2)%surface_model) == 'conductor', 'second template surface_model mismatch')
  call assert_true(trim(cfg%templates(3)%surface_model) == 'dielectric', 'third template surface_model mismatch')
  call assert_close_dp(cfg%templates(3)%epsilon_r, 3.9d0, 1.0d-15, 'third template epsilon_r mismatch')
  call assert_close_dp(cfg%sim%tree_theta, 0.35d0, 1.0d-15, 'tree_theta mismatch')
  call assert_true(cfg%sim%has_tree_theta, 'has_tree_theta mismatch')
  call assert_equal_i32(cfg%sim%tree_leaf_max, 12_i32, 'tree_leaf_max mismatch')
  call assert_true(cfg%sim%has_tree_leaf_max, 'has_tree_leaf_max mismatch')
  call assert_equal_i32(cfg%sim%tree_min_nelem, 1024_i32, 'tree_min_nelem mismatch')
  call test_end()

  call test_begin('split_periodic_outer_config')
  call default_app_config(split_cfg)
  call load_app_config(split_cfg_path, split_cfg)
  call assert_true( &
    trim(split_cfg%periodic2%nonzero_mode_backend) == 'panel_spectral_reference', &
    'split nonzero backend mismatch' &
    )
  call assert_true(trim(split_cfg%periodic2%zero_mode_policy) == 'exclude_k0', 'split zero policy mismatch')
  call assert_equal_i32(split_cfg%periodic2%reference_mode_layers, 5_i32, 'split mode layers mismatch')
  call assert_equal_i32(split_cfg%periodic2%panel_quadrature_order, 16_i32, 'split quadrature order mismatch')
  call assert_equal_i32(split_cfg%periodic2%interface_sample_n, 7_i32, 'split interface sample mismatch')
  call assert_true(trim(split_cfg%outer_plasma%model) == 'linear_debye', 'split outer model mismatch')
  call assert_close_dp(split_cfg%outer_plasma%interface_z, 1.0_dp, 1.0e-15_dp, 'split interface mismatch')
  call assert_close_dp(split_cfg%outer_plasma%debye_length, 0.2_dp, 1.0e-15_dp, 'split Debye length mismatch')
  call assert_equal_i32(split_cfg%outer_plasma%unified_grid_points, 65_i32, 'unified grid points mismatch')
  call assert_close_dp( &
    split_cfg%outer_plasma%accessible_fraction_tolerance, 0.04_dp, 1.0e-15_dp, &
    'accessible-fraction tolerance mismatch' &
    )
  call assert_close_dp(split_cfg%outer_plasma%max_gap_ratio, 4.0_dp, 1.0e-15_dp, 'split gap limit mismatch')
  call assert_close_dp( &
    split_cfg%outer_plasma%max_local_charge_ratio, 6.0_dp, 1.0e-15_dp, 'split local-charge limit mismatch' &
    )
  call assert_true( &
    trim(split_cfg%outer_plasma%photoelectron_closure) == 'individual_return', &
    'photoelectron closure mismatch' &
    )
  call assert_equal_i32( &
    split_cfg%outer_plasma%photoelectron_histogram_bins, 8_i32, 'photoelectron histogram bins mismatch' &
    )
  call assert_close_dp( &
    split_cfg%outer_plasma%photoelectron_histogram_energy_max, 12.0_dp, 1.0e-15_dp, &
    'photoelectron histogram energy mismatch' &
    )
  call assert_close_dp( &
    split_cfg%outer_plasma%photoelectron_ambient_charge_scale, 4.0_dp, 1.0e-15_dp, &
    'photoelectron ambient charge scale mismatch' &
    )
  call assert_close_dp( &
    split_cfg%outer_plasma%max_photoelectron_charge_ratio, 0.2_dp, 1.0e-15_dp, &
    'photoelectron charge ratio mismatch' &
    )
  call assert_true( &
    trim(split_cfg%coupling%particle_transfer_mode) == 'electrostatic_1d_instant_return', &
    'split transfer mode mismatch' &
    )
  call assert_close_dp( &
    split_cfg%coupling%field_evolution_timescale, 2.0_dp, 1.0e-15_dp, 'split field timescale mismatch' &
    )
  call assert_close_dp(split_cfg%coupling%outer_orbit_dt, 0.01_dp, 1.0e-15_dp, 'outer orbit dt mismatch')
  call assert_equal_i32(split_cfg%coupling%outer_orbit_max_steps, 400_i32, 'outer orbit max steps mismatch')
  call assert_close_dp( &
    split_cfg%coupling%outer_orbit_energy_tolerance, 2.0e-4_dp, 1.0e-15_dp, &
    'outer orbit energy tolerance mismatch' &
    )
  call test_end()

  call test_begin('triangle_panel_config')
  call default_app_config(panel_cfg)
  call load_app_config(panel_cfg_path, panel_cfg)
  call assert_true(trim(panel_cfg%panel%source_model) == 'triangle_p0', 'panel source model mismatch')
  call assert_true(trim(panel_cfg%panel%kernel_id) == 'triangle_p0_exact_direct', 'panel kernel id mismatch')
  call assert_true(trim(panel_cfg%templates(1)%surface_side_policy) == 'normal_plus', 'template side policy mismatch')
  call test_end()

  call test_begin('photo_raycast_config')
  call default_app_config(photo_cfg)
  call load_app_config(photo_cfg_path, photo_cfg)

  call assert_equal_i32(photo_cfg%n_particle_species, 1_i32, 'photo n_particle_species mismatch')
  call assert_true(trim(photo_cfg%particle_species(1)%source_mode) == 'photo_raycast', 'photo source_mode mismatch')
  call assert_close_dp(photo_cfg%particle_species(1)%emit_current_density_a_m2, 2.0d-3, 1.0d-15, 'photo emit_current mismatch')
  call assert_equal_i32(photo_cfg%particle_species(1)%rays_per_batch, 40_i32, 'photo rays_per_batch mismatch')
  call assert_close_dp(photo_cfg%particle_species(1)%normal_drift_speed, 1.5d5, 1.0d-12, 'photo normal_drift_speed mismatch')
  call assert_true( &
    photo_cfg%particle_species(1)%deposit_opposite_charge_on_emit, &
    'photo deposit_opposite_charge_on_emit mismatch' &
    )
  call assert_true( &
    trim(photo_cfg%particle_species(1)%photo_escape_model) == 'boltzmann_cutoff', &
    'photo_escape_model mismatch' &
    )
  call assert_allclose_1d( &
    photo_cfg%particle_species(1)%ray_direction, [0.0d0, 0.0d0, -1.0d0], 1.0d-12, 'photo ray_direction mismatch' &
    )
  call assert_equal_i32(photo_cfg%sim%raycast_max_bounce, 16_i32, 'photo default raycast_max_bounce mismatch')
  call assert_equal_i32(particles_per_batch_from_config(photo_cfg), 0_i32, 'photo per-batch particle count mismatch')
  call assert_equal_i32(total_particles_from_config(photo_cfg), 0_i32, 'photo total particle count mismatch')
  call assert_equal_i32(photo_cfg%n_particles, 0_i32, 'photo cached n_particles mismatch')
  call test_end()

  call test_begin('large_species_config')
  call default_app_config(large_cfg)
  call load_app_config(large_cfg_path, large_cfg)

  call assert_equal_i32(large_cfg%n_particle_species, 12_i32, 'large n_particle_species mismatch')
  call assert_equal_i32(large_cfg%n_templates, 12_i32, 'large n_templates mismatch')
  call assert_true(trim(large_cfg%templates(12)%kind) == 'sphere', 'large 12th template kind mismatch')
  call assert_equal_i32(large_cfg%particle_species(12)%npcls_per_step, 1_i32, 'large 12th species npcls mismatch')
  call assert_equal_i32(particles_per_batch_from_config(large_cfg), 12_i32, 'large per-batch particle count mismatch')
  call assert_equal_i32(total_particles_from_config(large_cfg), 24_i32, 'large total particle count mismatch')
  call test_end()

  call test_begin('periodic_config')
  call default_app_config(periodic_cfg)
  call load_app_config(periodic_cfg_path, periodic_cfg)
  call assert_true(trim(periodic_cfg%sim%field_solver) == 'fmm', 'periodic field_solver mismatch')
  call assert_true(trim(periodic_cfg%sim%field_bc_mode) == 'periodic2', 'periodic field_bc_mode mismatch')
  call assert_true(periodic_cfg%sim%use_box, 'periodic use_box mismatch')
  call assert_equal_i32(periodic_cfg%sim%bc_low(1), bc_periodic, 'periodic bc_x_low mismatch')
  call assert_equal_i32(periodic_cfg%sim%bc_high(1), bc_periodic, 'periodic bc_x_high mismatch')
  call assert_equal_i32(periodic_cfg%sim%bc_low(2), bc_periodic, 'periodic bc_y_low mismatch')
  call assert_equal_i32(periodic_cfg%sim%bc_high(2), bc_periodic, 'periodic bc_y_high mismatch')
  call assert_equal_i32(periodic_cfg%sim%field_periodic_image_layers, 2_i32, 'periodic field_periodic_image_layers mismatch')
  call assert_true( &
    trim(periodic_cfg%sim%field_periodic_far_correction) == 'cached_kneq0', &
    'periodic far correction mismatch' &
    )
  call assert_true( &
    trim(periodic_cfg%periodic2%nonzero_mode_backend) == 'cached_kneq0', &
    'periodic typed backend mismatch' &
    )
  call assert_true( &
    trim(periodic_cfg%periodic2%zero_mode_policy) == 'exclude_k0', &
    'periodic typed zero policy mismatch' &
    )
  call assert_true(trim(periodic_cfg%sim%field_periodic_cache_dir) == 'test-cache', 'periodic cache dir mismatch')
  call assert_close_dp( &
    periodic_cfg%sim%field_periodic_generation_tolerance, 2.5e-7_dp, 1.0e-18_dp, &
    'periodic generation tolerance mismatch' &
    )
  call assert_close_dp(periodic_cfg%sim%field_periodic_ewald_alpha, 1.5d0, 1.0d-12, 'periodic ewald alpha mismatch')
  call assert_equal_i32(periodic_cfg%sim%field_periodic_ewald_layers, 5_i32, 'periodic ewald layers mismatch')
  call test_end()

  call test_begin('periodic_oracle_config')
  call default_app_config(periodic_oracle_cfg)
  call load_app_config(periodic_oracle_cfg_path, periodic_oracle_cfg)
  call assert_true(trim(periodic_oracle_cfg%sim%field_bc_mode) == 'periodic2', 'periodic oracle field_bc_mode mismatch')
  call assert_true( &
    trim(periodic_oracle_cfg%sim%field_periodic_far_correction) == 'm2l_root_oracle', &
    'periodic oracle far correction mismatch' &
    )
  call assert_true( &
    trim(periodic_oracle_cfg%periodic2%nonzero_mode_backend) == 'legacy_root_oracle', &
    'periodic oracle typed backend mismatch' &
    )
  call assert_true( &
    trim(periodic_oracle_cfg%periodic2%zero_mode_policy) == 'legacy_charged_walls', &
    'periodic oracle typed zero policy mismatch' &
    )
  call assert_equal_i32( &
    periodic_oracle_cfg%sim%field_periodic_ewald_layers, 3_i32, 'periodic oracle ewald layers mismatch' &
    )
  call assert_close_dp( &
    periodic_oracle_cfg%sim%field_periodic_ewald_alpha, 0.0d0, 1.0d-15, 'periodic oracle ewald alpha mismatch' &
    )
  call test_end()

  call test_begin('sheath_config')
  call default_app_config(sheath_cfg)
  call load_app_config(sheath_cfg_path, sheath_cfg)
  call assert_true(trim(sheath_cfg%sim%sheath_injection_model) == 'zhao_auto', 'sheath model mismatch')
  call assert_close_dp(sheath_cfg%sim%sheath_alpha_deg, 42.0d0, 1.0d-12, 'sheath alpha mismatch')
  call assert_close_dp( &
    sheath_cfg%sim%sheath_photoelectron_ref_density_cm3, 48.0d0, 1.0d-12, 'sheath photoelectron density mismatch' &
    )
  call assert_true(sheath_cfg%sim%has_sheath_reference_coordinate, 'sheath reference coordinate flag mismatch')
  call assert_close_dp( &
    sheath_cfg%sim%sheath_reference_coordinate, 0.25d0, 1.0d-12, 'sheath reference coordinate mismatch' &
    )
  call assert_true(trim(sheath_cfg%sim%sheath_electron_drift_mode) == 'full', 'sheath electron drift mode mismatch')
  call assert_true(trim(sheath_cfg%sim%sheath_ion_drift_mode) == 'normal', 'sheath ion drift mode mismatch')
  call test_end()

  call test_begin('high_level_config_normalization')
  call default_app_config(high_level_cfg)
  call load_app_config(high_level_cfg_path, high_level_cfg)
  call assert_allclose_1d(high_level_cfg%sim%box_min, [1.0d0, 2.0d0, 3.0d0], 1.0d-12, 'high-level box_min mismatch')
  call assert_allclose_1d(high_level_cfg%sim%box_max, [3.0d0, 6.0d0, 9.0d0], 1.0d-12, 'high-level box_max mismatch')
  call assert_allclose_1d( &
    high_level_cfg%particle_species(1)%pos_low, [1.5d0, 4.0d0, 9.0d0], 1.0d-12, &
    'face_fraction pos_low mismatch' &
    )
  call assert_allclose_1d( &
    high_level_cfg%particle_species(1)%pos_high, [2.5d0, 6.0d0, 9.0d0], 1.0d-12, &
    'face_fraction pos_high mismatch' &
    )
  call assert_allclose_1d( &
    high_level_cfg%templates(1)%center, [3.2d0, 6.0d0, 9.0d0], 1.0d-12, &
    'grouped template center mismatch' &
    )
  call assert_close_dp(high_level_cfg%templates(1)%size_x, 2.0d0, 1.0d-12, 'grouped template size_x mismatch')
  call assert_close_dp(high_level_cfg%templates(1)%size_y, 4.0d0, 1.0d-12, 'grouped template size_y mismatch')
  call assert_allclose_1d( &
    high_level_cfg%templates(2)%center, [2.2d0, 4.0d0, 3.1d0], 1.0d-12, &
    'box_anchor template center mismatch' &
    )
  call assert_close_dp(high_level_cfg%templates(2)%radius, 0.5d0, 1.0d-12, 'box_fraction sphere radius mismatch')
  call test_end()

  call test_begin('multiline_array_config')
  call default_app_config(multiline_cfg)
  call load_app_config(multiline_cfg_path, multiline_cfg)
  call assert_true(multiline_cfg%sim%has_e0_vector, 'multiline e0 flag mismatch')
  call assert_allclose_1d(multiline_cfg%sim%e0, [1.0d0, 2.0d0, 3.0d0], 1.0d-12, 'multiline e0 mismatch')
  call assert_allclose_1d(multiline_cfg%sim%b0, [4.0d0, 5.0d0, 6.0d0], 1.0d-12, 'multiline b0 mismatch')
  call assert_allclose_1d(multiline_cfg%sim%box_min, [-1.0d0, -2.0d0, -3.0d0], 1.0d-12, 'box_min mismatch')
  call assert_allclose_1d(multiline_cfg%sim%box_max, [1.0d0, 2.0d0, 3.0d0], 1.0d-12, 'box_max mismatch')
  call assert_allclose_1d( &
    multiline_cfg%particle_species(1)%pos_low, [-0.2d0, -0.1d0, 0.0d0], 1.0d-12, 'pos_low mismatch' &
    )
  call assert_allclose_1d( &
    multiline_cfg%particle_species(1)%pos_high, [0.2d0, 0.1d0, 0.5d0], 1.0d-12, 'pos_high mismatch' &
    )
  call assert_allclose_1d( &
    multiline_cfg%particle_species(1)%drift_velocity, [7.0d0, 8.0d0, 9.0d0], 1.0d-12, 'drift_velocity mismatch' &
    )
  call assert_allclose_1d(multiline_cfg%templates(1)%center, [0.1d0, 0.2d0, 0.3d0], 1.0d-12, 'center mismatch')
  call assert_equal_i32(particles_per_batch_from_config(multiline_cfg), 2_i32, 'multiline per-batch mismatch')
  call test_end()

  call test_begin('toml_native_syntax_config')
  call default_app_config(toml_syntax_cfg)
  call load_app_config(toml_syntax_cfg_path, toml_syntax_cfg)
  call assert_equal_i32(toml_syntax_cfg%sim%batch_count, 1_i32, 'dotted batch_count mismatch')
  call assert_allclose_1d(toml_syntax_cfg%sim%box_min, [0.0d0, 0.0d0, 0.0d0], 1.0d-12, 'dotted box_min mismatch')
  call assert_allclose_1d(toml_syntax_cfg%sim%box_max, [1.0d0, 1.0d0, 1.0d0], 1.0d-12, 'dotted box_max mismatch')
  call assert_equal_i32(toml_syntax_cfg%n_particle_species, 1_i32, 'inline species count mismatch')
  call assert_equal_i32(toml_syntax_cfg%particle_species(1)%npcls_per_step, 3_i32, 'inline species npcls mismatch')
  call assert_allclose_1d( &
    toml_syntax_cfg%particle_species(1)%drift_velocity, [1.0d0, 2.0d0, 3.0d0], 1.0d-12, &
    'inline species drift_velocity mismatch' &
    )
  call assert_equal_i32(toml_syntax_cfg%n_templates, 1_i32, 'inline template count mismatch')
  call assert_true(trim(toml_syntax_cfg%templates(1)%kind) == 'sphere', 'inline template kind mismatch')
  call assert_allclose_1d(toml_syntax_cfg%templates(1)%center, [0.1d0, 0.2d0, 0.3d0], 1.0d-12, 'inline center mismatch')
  call assert_true(trim(toml_syntax_cfg%output_dir) == 'outputs/#literal', 'literal hash in string mismatch')
  call test_end()

  call delete_file_if_exists(cfg_path)
  call delete_file_if_exists(photo_cfg_path)
  call delete_file_if_exists(large_cfg_path)
  call delete_file_if_exists(periodic_cfg_path)
  call delete_file_if_exists(periodic_oracle_cfg_path)
  call delete_file_if_exists(sheath_cfg_path)
  call delete_file_if_exists(high_level_cfg_path)
  call delete_file_if_exists(multiline_cfg_path)
  call delete_file_if_exists(toml_syntax_cfg_path)
  call delete_file_if_exists(panel_cfg_path)
  call delete_file_if_exists(split_cfg_path)

  call test_summary()

contains

  !> テスト専用の一時設定ファイルを書き出す。
  !! @param[in] path 書き出し先TOMLファイルパス。
  subroutine write_config_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'dt = 2.5e-9'
    write (u, '(a)') 'batch_count = 3 # comment should be ignored'
    write (u, '(a)') 'reservoir_potential_model = "infinity_barrier"'
    write (u, '(a)') 'phi_infty = -2.0'
    write (u, '(a)') 'open_boundary_model = "potential_barrier"'
    write (u, '(a)') 'injection_face_phi_grid_n = 5'
    write (u, '(a)') 'field_solver = "fmm"'
    write (u, '(a)') 'field_normalization = "length"'
    write (u, '(a)') 'field_length_scale = 2.0e-6'
    write (u, '(a)') 'field_bc_mode = "free"'
    write (u, '(a)') 'tree_theta = 0.35'
    write (u, '(a)') 'tree_leaf_max = 12'
    write (u, '(a)') 'tree_min_nelem = 1024'
    write (u, '(a)') 'bc_x_low = "periodic"'
    write (u, '(a)') 'bc_y_high = "reflect"'
    write (u, '(a)') 'bc_z_low = "open"'
    write (u, '(a)') ''
    write (u, '(a)') '[particles]'
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'species_key = "electron"'
    write (u, '(a)') 'npcls_per_step = 4'
    write (u, '(a)') 'temperature_k = 10000.0'
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'npcls_per_step = 1'
    write (u, '(a)') 'drift_velocity = [1.0, 0.0, -2.0]'
    write (u, '(a)') ''
    write (u, '(a)') '[mesh]'
    write (u, '(a)') 'mode = "template"'
    write (u, '(a)') 'surface_model = "dielectric"'
    write (u, '(a)') 'epsilon_r = 2.5'
    write (u, '(a)') '[[mesh.templates]]'
    write (u, '(a)') 'kind = "plane"'
    write (u, '(a)') 'enabled = true'
    write (u, '(a)') '[[mesh.templates]]'
    write (u, '(a)') 'kind = "sphere"'
    write (u, '(a)') 'enabled = true'
    write (u, '(a)') 'surface_model = "conductor"'
    write (u, '(a)') 'radius = 0.25'
    write (u, '(a)') '[[mesh.templates]]'
    write (u, '(a)') 'kind = "annulus"'
    write (u, '(a)') 'enabled = true'
    write (u, '(a)') 'surface_model = "dielectric"'
    write (u, '(a)') 'epsilon_r = 3.9'
    write (u, '(a)') 'radius = 0.3'
    write (u, '(a)') 'inner_radius = 0.1'
    write (u, '(a)') 'n_theta = 16'
    write (u, '(a)') 'n_r = 3'
    write (u, '(a)') '[[mesh.templates]]'
    write (u, '(a)') 'kind = "cylinder"'
    write (u, '(a)') 'enabled = true'
    write (u, '(a)') 'cap_top = true'
    write (u, '(a)') 'cap_bottom = false'
    write (u, '(a)') '[[mesh.templates]]'
    write (u, '(a)') 'kind = "plate_hole"'
    write (u, '(a)') 'enabled = true'
    write (u, '(a)') 'size_x = 1.5'
    write (u, '(a)') 'size_y = 0.8'
    write (u, '(a)') 'radius = 0.2'
    write (u, '(a)') 'n_theta = 20'
    write (u, '(a)') 'n_r = 2'
    write (u, '(a)') ''
    write (u, '(a)') '[output]'
    write (u, '(a)') 'history_stride = 2'
    write (u, '(a)') 'write_mesh_potential = true'
    write (u, '(a)') 'resume = true'
    write (u, '(a)') 'restart_from = "outputs/parent"'

    close (u)
  end subroutine write_config_fixture

  !> テスト専用の photo_raycast 設定ファイルを書き出す。
  !! @param[in] path 書き出し先TOMLファイルパス。
  subroutine write_photo_config_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open photo config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'batch_count = 4'
    write (u, '(a)') 'batch_duration = 1.0e-7'
    write (u, '(a)') 'use_box = true'
    write (u, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (u, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'source_mode = "photo_raycast"'
    write (u, '(a)') 'emit_current_density_a_m2 = 2.0e-3'
    write (u, '(a)') 'rays_per_batch = 40'
    write (u, '(a)') 'deposit_opposite_charge_on_emit = true'
    write (u, '(a)') 'photo_escape_model = "boltzmann_cutoff"'
    write (u, '(a)') 'normal_drift_speed = 1.5e5'
    write (u, '(a)') 'q_particle = -1.0'
    write (u, '(a)') 'm_particle = 1.0'
    write (u, '(a)') 'temperature_k = 0.0'
    write (u, '(a)') 'inject_face = "z_high"'
    write (u, '(a)') 'pos_low = [0.0, 0.0, 1.0]'
    write (u, '(a)') 'pos_high = [1.0, 1.0, 1.0]'

    close (u)
  end subroutine write_photo_config_fixture

  !> 既定初期容量を超える species/template を含む設定ファイルを書き出す。
  !! @param[in] path 書き出し先TOMLファイルパス。
  subroutine write_large_config_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios, i

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open large config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'batch_count = 2'
    write (u, '(a)') ''
    write (u, '(a)') '[particles]'
    do i = 1, 12
      write (u, '(a)') '[[particles.species]]'
      write (u, '(a)') 'npcls_per_step = 1'
    end do
    write (u, '(a)') ''
    write (u, '(a)') '[mesh]'
    write (u, '(a)') 'mode = "template"'
    do i = 1, 12
      write (u, '(a)') '[[mesh.templates]]'
      write (u, '(a)') 'kind = "sphere"'
      write (u, '(a)') 'enabled = true'
      write (u, '(a)') 'radius = 0.1'
      write (u, '(a)') 'center = [0.0, 0.0, 0.2]'
    end do

    close (u)
  end subroutine write_large_config_fixture

  !> `field_bc_mode="periodic2"` で cached nonzero-mode operator の設定を書き出す。
  !! @param[in] path 書き出し先TOMLファイルパス。
  subroutine write_periodic_config_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open periodic config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'batch_count = 2'
    write (u, '(a)') 'field_solver = "fmm"'
    write (u, '(a)') 'field_bc_mode = "periodic2"'
    write (u, '(a)') 'use_box = true'
    write (u, '(a)') 'box_min = [0.0, 0.0, -1.0]'
    write (u, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    write (u, '(a)') 'bc_x_low = "periodic"'
    write (u, '(a)') 'bc_x_high = "periodic"'
    write (u, '(a)') 'bc_y_low = "periodic"'
    write (u, '(a)') 'bc_y_high = "periodic"'
    write (u, '(a)') 'bc_z_low = "open"'
    write (u, '(a)') 'bc_z_high = "open"'
    write (u, '(a)') 'field_periodic_image_layers = 2'
    write (u, '(a)') 'field_periodic_far_correction = "cached_kneq0"'
    write (u, '(a)') 'field_periodic_ewald_alpha = 1.5'
    write (u, '(a)') 'field_periodic_ewald_layers = 5'
    write (u, '(a)') 'field_periodic_cache_dir = "test-cache"'
    write (u, '(a)') 'field_periodic_generation_tolerance = 2.5e-7'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'npcls_per_step = 1'
    write (u, '(a)') ''
    write (u, '(a)') '[mesh]'
    write (u, '(a)') 'mode = "template"'
    write (u, '(a)') '[[mesh.templates]]'
    write (u, '(a)') 'kind = "plane"'
    write (u, '(a)') 'enabled = true'

    close (u)
  end subroutine write_periodic_config_fixture

  !> `m2l_root_oracle` が受理される periodic2 設定を書き出す。
  subroutine write_periodic_oracle_config_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open periodic oracle config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'batch_count = 2'
    write (u, '(a)') 'field_solver = "fmm"'
    write (u, '(a)') 'field_bc_mode = "periodic2"'
    write (u, '(a)') 'use_box = true'
    write (u, '(a)') 'box_min = [0.0, 0.0, -1.0]'
    write (u, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    write (u, '(a)') 'bc_x_low = "periodic"'
    write (u, '(a)') 'bc_x_high = "periodic"'
    write (u, '(a)') 'bc_y_low = "periodic"'
    write (u, '(a)') 'bc_y_high = "periodic"'
    write (u, '(a)') 'bc_z_low = "open"'
    write (u, '(a)') 'bc_z_high = "open"'
    write (u, '(a)') 'field_periodic_image_layers = 1'
    write (u, '(a)') 'field_periodic_far_correction = "m2l_root_oracle"'
    write (u, '(a)') 'field_periodic_ewald_layers = 3'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'npcls_per_step = 1'
    write (u, '(a)') ''
    write (u, '(a)') '[mesh]'
    write (u, '(a)') 'mode = "template"'
    write (u, '(a)') '[[mesh.templates]]'
    write (u, '(a)') 'kind = "plane"'
    write (u, '(a)') 'enabled = true'

    close (u)
  end subroutine write_periodic_oracle_config_fixture

  !> シース注入キーの受理を確認する設定を書き出す。
  subroutine write_sheath_config_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open sheath config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'batch_count = 1'
    write (u, '(a)') 'sheath_injection_model = "zhao_auto"'
    write (u, '(a)') 'sheath_alpha_deg = 42.0'
    write (u, '(a)') 'sheath_photoelectron_ref_density_cm3 = 48.0'
    write (u, '(a)') 'sheath_reference_coordinate = 0.25'
    write (u, '(a)') 'sheath_electron_drift_mode = "full"'
    write (u, '(a)') 'sheath_ion_drift_mode = "normal"'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'npcls_per_step = 1'
    write (u, '(a)') ''
    write (u, '(a)') '[mesh]'
    write (u, '(a)') 'mode = "template"'
    write (u, '(a)') '[[mesh.templates]]'
    write (u, '(a)') 'kind = "plane"'
    write (u, '(a)') 'enabled = true'

    close (u)
  end subroutine write_sheath_config_fixture

  !> 高水準 authoring キーを Fortran 側で正規化することを確認する設定を書き出す。
  subroutine write_high_level_config_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open high-level config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'batch_count = 1'
    write (u, '(a)') 'batch_duration = 1.0'
    write (u, '(a)') 'use_box = true'
    write (u, '(a)') 'box_origin = [1.0, 2.0, 3.0]'
    write (u, '(a)') 'box_size = [2.0, 4.0, 6.0]'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'source_mode = "reservoir_face"'
    write (u, '(a)') 'number_density_m3 = 1000.0'
    write (u, '(a)') 'temperature_k = 0.0'
    write (u, '(a)') 'q_particle = -1.0'
    write (u, '(a)') 'm_particle = 1.0'
    write (u, '(a)') 'w_particle = 100.0'
    write (u, '(a)') 'inject_face = "z_high"'
    write (u, '(a)') 'inject_region_mode = "face_fraction"'
    write (u, '(a)') 'uv_low = [0.25, 0.5]'
    write (u, '(a)') 'uv_high = [0.75, 1.0]'
    write (u, '(a)') 'drift_velocity = [0.0, 0.0, -1.0]'
    write (u, '(a)') ''
    write (u, '(a)') '[mesh]'
    write (u, '(a)') 'mode = "template"'
    write (u, '(a)') ''
    write (u, '(a)') '[mesh.groups.unit]'
    write (u, '(a)') 'placement_mode = "box_anchor"'
    write (u, '(a)') 'anchor = "box_center"'
    write (u, '(a)') 'offset_frac = [0.1, 0.0, 0.0]'
    write (u, '(a)') 'scale_from = "box_min_xy"'
    write (u, '(a)') 'scale_factor = 0.5'
    write (u, '(a)') ''
    write (u, '(a)') '[[mesh.templates]]'
    write (u, '(a)') 'kind = "plane"'
    write (u, '(a)') 'group = "unit"'
    write (u, '(a)') 'center_local = [1.0, 2.0, 3.0]'
    write (u, '(a)') 'size_x = 2.0'
    write (u, '(a)') 'size_y = 4.0'
    write (u, '(a)') 'nx = 4'
    write (u, '(a)') 'ny = 5'
    write (u, '(a)') ''
    write (u, '(a)') '[[mesh.templates]]'
    write (u, '(a)') 'kind = "sphere"'
    write (u, '(a)') 'placement_mode = "box_anchor"'
    write (u, '(a)') 'anchor = "z_low_face_center"'
    write (u, '(a)') 'offset = [0.2, 0.0, 0.1]'
    write (u, '(a)') 'size_mode = "box_fraction"'
    write (u, '(a)') 'size_frac = 0.25'
    write (u, '(a)') 'n_lon = 8'
    write (u, '(a)') 'n_lat = 6'

    close (u)
  end subroutine write_high_level_config_fixture

  !> 標準 TOML の複数行配列を含む設定を書き出す。
  subroutine write_multiline_config_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open multiline config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'batch_count = 1'
    write (u, '(a)') 'use_box = true'
    write (u, '(a)') 'e0 = ['
    write (u, '(a)') '  1.0,'
    write (u, '(a)') '  2.0,'
    write (u, '(a)') '  3.0,'
    write (u, '(a)') ']'
    write (u, '(a)') 'b0 = ['
    write (u, '(a)') '  4.0,'
    write (u, '(a)') '  5.0,'
    write (u, '(a)') '  6.0,'
    write (u, '(a)') ']'
    write (u, '(a)') 'box_min = ['
    write (u, '(a)') '  -1.0,'
    write (u, '(a)') '  -2.0,'
    write (u, '(a)') '  -3.0,'
    write (u, '(a)') ']'
    write (u, '(a)') 'box_max = ['
    write (u, '(a)') '  1.0,'
    write (u, '(a)') '  2.0,'
    write (u, '(a)') '  3.0,'
    write (u, '(a)') ']'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'npcls_per_step = 2'
    write (u, '(a)') 'pos_low = ['
    write (u, '(a)') '  -0.2,'
    write (u, '(a)') '  -0.1,'
    write (u, '(a)') '  0.0,'
    write (u, '(a)') ']'
    write (u, '(a)') 'pos_high = ['
    write (u, '(a)') '  0.2,'
    write (u, '(a)') '  0.1,'
    write (u, '(a)') '  0.5,'
    write (u, '(a)') ']'
    write (u, '(a)') 'drift_velocity = ['
    write (u, '(a)') '  7.0,'
    write (u, '(a)') '  8.0,'
    write (u, '(a)') '  9.0,'
    write (u, '(a)') ']'
    write (u, '(a)') ''
    write (u, '(a)') '[mesh]'
    write (u, '(a)') 'mode = "template"'
    write (u, '(a)') '[[mesh.templates]]'
    write (u, '(a)') 'kind = "sphere"'
    write (u, '(a)') 'enabled = true'
    write (u, '(a)') 'radius = 0.2'
    write (u, '(a)') 'center = ['
    write (u, '(a)') '  0.1,'
    write (u, '(a)') '  0.2,'
    write (u, '(a)') '  0.3,'
    write (u, '(a)') ']'

    close (u)
  end subroutine write_multiline_config_fixture

  !> dotted key と inline table 配列を含む標準 TOML 設定を書き出す。
  subroutine write_toml_syntax_config_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open TOML syntax config fixture'

    write (u, '(a)') 'sim.batch_count = 1'
    write (u, '(a)') 'sim.use_box = true'
    write (u, '(a)') 'sim.box_min = [0.0, 0.0, 0.0]'
    write (u, '(a)') 'sim.box_max = [1.0, 1.0, 1.0]'
    write (u, '(a)') ''
    write (u, '(a)') 'particles.species = ['
    write (u, '(a)') '  {'
    write (u, '(a)') '    npcls_per_step = 3,'
    write (u, '(a)') '    drift_velocity = [1.0, 2.0, 3.0],'
    write (u, '(a)') '  },'
    write (u, '(a)') ']'
    write (u, '(a)') ''
    write (u, '(a)') 'mesh.mode = "template"'
    write (u, '(a)') 'mesh.templates = ['
    write (u, '(a)') '  {'
    write (u, '(a)') '    kind = "sphere",'
    write (u, '(a)') '    enabled = true,'
    write (u, '(a)') '    radius = 0.2,'
    write (u, '(a)') '    center = [0.1, 0.2, 0.3],'
    write (u, '(a)') '  },'
    write (u, '(a)') ']'
    write (u, '(a)') ''
    write (u, '(a)') 'output.dir = "outputs/#literal"'

    close (u)
  end subroutine write_toml_syntax_config_fixture

  subroutine write_panel_config_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=path, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open panel config fixture'
    write (u, '(a)') '[sim]'
    write (u, '(a)') 'batch_count = 1'
    write (u, '(a)') 'dt = 1.0e-9'
    write (u, '(a)') 'max_step = 2'
    write (u, '(a)') 'field_solver = "direct"'
    write (u, '(a)') 'field_bc_mode = "free"'
    write (u, '(a)') 'softening = 0.0'
    write (u, '(a)') '[field]'
    write (u, '(a)') 'element_kernel = "triangle_p0"'
    write (u, '(a)') '[particles]'
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'npcls_per_step = 1'
    write (u, '(a)') '[mesh]'
    write (u, '(a)') 'mode = "template"'
    write (u, '(a)') '[[mesh.templates]]'
    write (u, '(a)') 'kind = "plane"'
    write (u, '(a)') 'surface_side = "normal_plus"'
    write (u, '(a)') 'nx = 1'
    write (u, '(a)') 'ny = 1'
    close (u)
  end subroutine write_panel_config_fixture

  subroutine write_split_config_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=path, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open split config fixture'
    write (u, '(a)') '[sim]'
    write (u, '(a)') 'batch_count = 1'
    write (u, '(a)') 'dt = 1.0e-9'
    write (u, '(a)') 'max_step = 2'
    write (u, '(a)') 'field_solver = "direct"'
    write (u, '(a)') 'field_bc_mode = "periodic2"'
    write (u, '(a)') 'softening = 0.0'
    write (u, '(a)') 'use_box = true'
    write (u, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (u, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    write (u, '(a)') 'bc_x_low = "periodic"'
    write (u, '(a)') 'bc_x_high = "periodic"'
    write (u, '(a)') 'bc_y_low = "periodic"'
    write (u, '(a)') 'bc_y_high = "periodic"'
    write (u, '(a)') 'bc_z_low = "open"'
    write (u, '(a)') 'bc_z_high = "open"'
    write (u, '(a)') '[field]'
    write (u, '(a)') 'element_kernel = "triangle_p0"'
    write (u, '(a)') '[periodic2]'
    write (u, '(a)') 'nonzero_mode_backend = "panel_spectral_reference"'
    write (u, '(a)') 'zero_mode_policy = "exclude_k0"'
    write (u, '(a)') 'lower_boundary_model = "e_bottom_zero"'
    write (u, '(a)') 'reference_mode_layers = 5'
    write (u, '(a)') 'panel_quadrature_order = 16'
    write (u, '(a)') 'interface_sample_n = 7'
    write (u, '(a)') '[outer_plasma]'
    write (u, '(a)') 'model = "linear_debye"'
    write (u, '(a)') 'photoelectron_closure = "individual_return"'
    write (u, '(a)') 'return_model = "electrostatic_1d_instant_return"'
    write (u, '(a)') 'interface_z = 1.0'
    write (u, '(a)') 'infinity_potential = 0.0'
    write (u, '(a)') 'debye_length = 0.2'
    write (u, '(a)') 'thermal_voltage = 10.0'
    write (u, '(a)') 'unified_grid_points = 65'
    write (u, '(a)') 'accessible_fraction_tolerance = 0.04'
    write (u, '(a)') 'max_linearity_ratio = 0.5'
    write (u, '(a)') 'max_gap_ratio = 4.0'
    write (u, '(a)') 'max_local_charge_ratio = 6.0'
    write (u, '(a)') 'photoelectron_histogram_bins = 8'
    write (u, '(a)') 'photoelectron_histogram_energy_max = 12.0'
    write (u, '(a)') 'photoelectron_ambient_charge_scale = 4.0'
    write (u, '(a)') 'max_photoelectron_charge_ratio = 0.2'
    write (u, '(a)') '[coupling]'
    write (u, '(a)') 'update_mode = "explicit"'
    write (u, '(a)') 'particle_transfer_mode = "electrostatic_1d_instant_return"'
    write (u, '(a)') 'outer_update_stride = 1'
    write (u, '(a)') 'field_evolution_timescale = 2.0'
    write (u, '(a)') 'max_frozen_field_ratio = 0.1'
    write (u, '(a)') 'outer_orbit_dt = 0.01'
    write (u, '(a)') 'outer_orbit_max_steps = 400'
    write (u, '(a)') 'outer_orbit_energy_tolerance = 2.0e-4'
    write (u, '(a)') 'outer_queue_enabled = false'
    write (u, '(a)') '[particles]'
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'npcls_per_step = 1'
    write (u, '(a)') '[mesh]'
    write (u, '(a)') 'mode = "template"'
    write (u, '(a)') '[[mesh.templates]]'
    write (u, '(a)') 'kind = "plane"'
    write (u, '(a)') 'surface_side = "normal_plus"'
    write (u, '(a)') 'center = [0.5, 0.5, 0.25]'
    close (u)
  end subroutine write_split_config_fixture

end program test_app_config_parser

!> app_config の設定読込と派生値計算を検証するテスト。
program test_app_config_parser
  use bem_kinds, only: dp, i32
  use bem_types, only: bc_open, bc_reflect, bc_periodic
  use bem_app_config, only: app_config, default_app_config, load_app_config, &
                            particles_per_batch_from_config, total_particles_from_config
  use bem_app_config_authoring, only: app_config_authoring, lower_external_boundary_authoring
  use bem_model_fingerprint, only: model_fingerprint
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, assert_allclose_1d, delete_file_if_exists
  implicit none

  type(app_config) :: cfg, photo_cfg, large_cfg, periodic_cfg, zhao_cfg, high_level_cfg
  type(app_config) :: multiline_cfg, toml_syntax_cfg, panel_cfg, panel_tree_cfg, split_cfg
  type(app_config) :: external_boundary_cfg, none_boundary_cfg
  type(app_config) :: implicit_mean_cfg, ambient_only_cfg, positive_photo_cfg, non_template_mean_cfg
  character(len=*), parameter :: cfg_path = 'test_app_config_parser_tmp.toml'
  character(len=*), parameter :: photo_cfg_path = 'test_app_config_parser_photo_tmp.toml'
  character(len=*), parameter :: large_cfg_path = 'test_app_config_parser_large_tmp.toml'
  character(len=*), parameter :: periodic_cfg_path = 'test_app_config_parser_periodic_tmp.toml'
  character(len=*), parameter :: periodic_oracle_cfg_path = 'test_app_config_parser_periodic_oracle_tmp.toml'
  character(len=*), parameter :: periodic_oracle_output_path = 'test_app_config_parser_periodic_oracle_tmp.out'
  character(len=*), parameter :: zhao_cfg_path = 'test_app_config_parser_zhao_tmp.toml'
  character(len=*), parameter :: high_level_cfg_path = 'test_app_config_parser_high_level_tmp.toml'
  character(len=*), parameter :: multiline_cfg_path = 'test_app_config_parser_multiline_tmp.toml'
  character(len=*), parameter :: toml_syntax_cfg_path = 'test_app_config_parser_toml_syntax_tmp.toml'
  character(len=*), parameter :: panel_cfg_path = 'test_app_config_parser_panel_tmp.toml'
  character(len=*), parameter :: panel_tree_cfg_path = 'test_app_config_parser_panel_tree_tmp.toml'
  character(len=*), parameter :: split_cfg_path = 'test_app_config_parser_split_tmp.toml'
  character(len=*), parameter :: external_boundary_cfg_path = 'test_app_config_parser_external_boundary_tmp.toml'
  character(len=*), parameter :: mixed_boundary_cfg_path = 'test_app_config_parser_mixed_boundary_tmp.toml'
  character(len=*), parameter :: mixed_boundary_output_path = 'test_app_config_parser_mixed_boundary_tmp.out'
  character(len=*), parameter :: noop_boundary_cfg_path = 'test_app_config_parser_noop_boundary_tmp.toml'
  character(len=*), parameter :: noop_boundary_output_path = 'test_app_config_parser_noop_boundary_tmp.out'
  character(len=*), parameter :: none_boundary_cfg_path = 'test_app_config_parser_none_boundary_tmp.toml'
  character(len=*), parameter :: implicit_mean_cfg_path = 'test_app_config_parser_implicit_mean_tmp.toml'
  character(len=*), parameter :: ambient_only_cfg_path = 'test_app_config_parser_ambient_only_tmp.toml'
  character(len=*), parameter :: positive_photo_cfg_path = 'test_app_config_parser_positive_photo_tmp.toml'
  character(len=*), parameter :: missing_deposit_cfg_path = 'test_app_config_parser_missing_deposit_tmp.toml'
  character(len=*), parameter :: missing_deposit_output_path = 'test_app_config_parser_missing_deposit_tmp.out'
  character(len=*), parameter :: non_template_mean_cfg_path = 'test_app_config_parser_non_template_mean_tmp.toml'
  character(len=*), parameter :: authored_implicit_cfg_path = 'test_app_config_parser_authored_implicit_tmp.toml'
  character(len=*), parameter :: authored_implicit_output_path = 'test_app_config_parser_authored_implicit_tmp.out'
  character(len=64) :: run_mode, probe_argument, noop_kind

  call get_command_argument(1, run_mode)
  if (trim(run_mode) == 'probe_authored_implicit_mean') then
    call write_authored_implicit_mean_config_fixture(authored_implicit_cfg_path)
    call default_app_config(cfg)
    call load_app_config(authored_implicit_cfg_path, cfg)
    error stop 'user-authored implicit_mean probe unexpectedly completed'
  end if
  if (trim(run_mode) == 'probe_missing_1d_photo_deposit') then
    call write_ambient_linear_debye_config_fixture( &
      missing_deposit_cfg_path, include_photo=.true., deposit_countercharge=.false. &
      )
    call default_app_config(cfg)
    call load_app_config(missing_deposit_cfg_path, cfg)
    error stop 'missing 1D photoelectron countercharge probe unexpectedly completed'
  end if
  if (trim(run_mode) == 'probe_mixed_external_boundary_owner') then
    call get_command_argument(2, probe_argument)
    call write_external_boundary_config_fixture(mixed_boundary_cfg_path, conflict_owner=trim(probe_argument))
    call default_app_config(cfg)
    call load_app_config(mixed_boundary_cfg_path, cfg)
    error stop 'mixed external-boundary owner probe unexpectedly completed'
  end if
  noop_kind = ''
  select case (trim(run_mode))
  case ('probe_none_field_option')
    noop_kind = 'field'
  case ('probe_none_local_coupling_option')
    noop_kind = 'coupling'
  case ('probe_kinetic_zhao_option')
    noop_kind = 'kinetic_zhao'
  case ('probe_local_time_guard')
    noop_kind = 'local_time_guard'
  case ('probe_queue_update_stride')
    noop_kind = 'queue_stride'
  case ('probe_kinetic_steady_start')
    noop_kind = 'kinetic_steady'
  end select
  if (len_trim(noop_kind) > 0) then
    call write_external_boundary_noop_fixture(noop_boundary_cfg_path, trim(noop_kind))
    call default_app_config(cfg)
    call load_app_config(noop_boundary_cfg_path, cfg)
    error stop 'external-boundary applicability probe unexpectedly completed'
  end if
  if (trim(run_mode) == 'probe_active_field_without_box') then
    call write_external_boundary_config_fixture(noop_boundary_cfg_path, omit_box=.true.)
    call default_app_config(cfg)
    call load_app_config(noop_boundary_cfg_path, cfg)
    error stop 'active field without explicit box probe unexpectedly completed'
  end if
  if (trim(run_mode) == 'probe_removed_root_oracle') then
    call write_periodic_oracle_config_fixture(periodic_oracle_cfg_path)
    call default_app_config(cfg)
    call load_app_config(periodic_oracle_cfg_path, cfg)
    error stop 'removed root oracle probe unexpectedly completed'
  end if
  call write_config_fixture(cfg_path)
  call write_photo_config_fixture(photo_cfg_path)
  call write_large_config_fixture(large_cfg_path)
  call write_periodic_config_fixture(periodic_cfg_path)
  call write_zhao_physics_config_fixture(zhao_cfg_path)
  call write_high_level_config_fixture(high_level_cfg_path)
  call write_multiline_config_fixture(multiline_cfg_path)
  call write_toml_syntax_config_fixture(toml_syntax_cfg_path)
  call write_panel_config_fixture(panel_cfg_path, 'direct')
  call write_panel_config_fixture(panel_tree_cfg_path, 'treecode')
  call write_split_config_fixture(split_cfg_path)
  call write_external_boundary_config_fixture(external_boundary_cfg_path)
  call write_external_boundary_noop_fixture(none_boundary_cfg_path, 'none')
  call write_ambient_linear_debye_config_fixture(implicit_mean_cfg_path, include_photo=.true., deposit_countercharge=.true.)
  call write_ambient_linear_debye_config_fixture(ambient_only_cfg_path, include_photo=.false., deposit_countercharge=.false.)
  call write_ambient_linear_debye_config_fixture( &
    positive_photo_cfg_path, include_photo=.true., deposit_countercharge=.true., &
    photo_charge=1.602176634e-19_dp &
    )
  call write_ambient_linear_debye_config_fixture( &
    non_template_mean_cfg_path, include_photo=.true., deposit_countercharge=.true., mesh_mode='obj' &
    )

  call test_init(20)

  call test_begin('defaults_and_basic_config')
  call default_app_config(cfg)
  call assert_true(trim(cfg%mesh_mode) == 'template', 'default mesh_mode mismatch')
  call assert_true(trim(cfg%sim%field_solver) == 'auto', 'default field_solver mismatch')
  call assert_true(trim(cfg%sim%field_normalization) == 'si', 'default field_normalization mismatch')
  call assert_close_dp(cfg%sim%field_length_scale, 1.0d0, 1.0d-15, 'default field_length_scale mismatch')
  call assert_true(trim(cfg%mesh_surface_model) == 'insulator', 'default mesh surface_model mismatch')
  call assert_close_dp(cfg%mesh_epsilon_r, 1.0d0, 1.0d-15, 'default mesh epsilon_r mismatch')
  call assert_true(trim(cfg%sim%field_bc_mode) == 'free', 'default field_bc_mode mismatch')
  call assert_equal_i32(cfg%sim%field_periodic_image_layers, 1_i32, 'default field_periodic_image_layers mismatch')
  call assert_true(trim(cfg%sim%field_periodic_far_correction) == 'none', 'default field_periodic_far_correction mismatch')
  call assert_true(trim(cfg%field%backend) == 'auto', 'default typed field backend mismatch')
  call assert_true( &
    trim(cfg%periodic2%nonzero_mode_backend) == 'not_applicable', &
    'default typed periodic backend mismatch' &
    )
  call assert_true(trim(cfg%outer_plasma%model) == 'none', 'default typed outer model mismatch')
  call assert_true( &
    trim(cfg%outer_plasma%kinetic_closure) == 'absorbing_maxwellian', &
    'default kinetic closure mismatch' &
    )
  call assert_true(trim(cfg%outer_plasma%zhao_branch) == 'auto', 'default Zhao branch mismatch')
  call assert_true(trim(cfg%coupling%steady_start_mode) == 'none', 'default steady-start mode mismatch')
  call assert_equal_i32(cfg%coupling%steady_start_mesh_id, 1_i32, 'default steady-start mesh ID mismatch')
  call assert_close_dp( &
    cfg%outer_plasma%photoelectron_source_scale, 1.0_dp, 0.0_dp, &
    'default photoelectron source scale mismatch' &
    )
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
  call assert_true(trim(cfg%sim%multiple_box_events_policy) == 'abort', 'default multiple box event policy mismatch')
  call assert_equal_i32( &
    cfg%sim%multiple_box_events_soft_discard_count_limit, 1000_i32, &
    'default multiple box event soft discard count limit mismatch' &
    )
  call assert_close_dp( &
    cfg%sim%multiple_box_events_soft_discard_abs_charge_limit, 1.0d-12, 1.0d-24, &
    'default multiple box event soft discard charge limit mismatch' &
    )
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
  call assert_true(trim(cfg%sim%multiple_box_events_policy) == 'soft_discard', 'multiple box event policy mismatch')
  call assert_equal_i32( &
    cfg%sim%multiple_box_events_soft_discard_count_limit, 100000_i32, &
    'multiple box event soft discard count limit mismatch' &
    )
  call assert_close_dp( &
    cfg%sim%multiple_box_events_soft_discard_abs_charge_limit, 1.0d-9, 1.0d-21, &
    'multiple box event soft discard charge limit mismatch' &
    )
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
  call assert_true( &
    trim(split_cfg%periodic2%lower_boundary_model) == 'symmetric_vacuum', &
    'split lower boundary mismatch' &
    )
  call assert_equal_i32(split_cfg%periodic2%reference_mode_layers, 5_i32, 'split mode layers mismatch')
  call assert_equal_i32(split_cfg%periodic2%panel_quadrature_order, 16_i32, 'split quadrature order mismatch')
  call assert_equal_i32(split_cfg%periodic2%interface_sample_n, 7_i32, 'split interface sample mismatch')
  call assert_true(trim(split_cfg%outer_plasma%model) == 'kinetic_1d', 'split outer model mismatch')
  call assert_true( &
    trim(split_cfg%outer_plasma%kinetic_closure) == 'absorbing_maxwellian', &
    'split default kinetic closure mismatch' &
    )
  call assert_true(trim(split_cfg%outer_plasma%zhao_branch) == 'auto', 'split default Zhao branch mismatch')
  call assert_close_dp( &
    split_cfg%outer_plasma%photoelectron_source_scale, 1.0_dp, 0.0_dp, &
    'split default photoelectron source scale mismatch' &
    )
  call assert_close_dp(split_cfg%outer_plasma%interface_z, 1.0_dp, 1.0e-15_dp, 'split interface mismatch')
  call assert_close_dp(split_cfg%outer_plasma%debye_length, 0.2_dp, 1.0e-15_dp, 'split Debye length mismatch')
  call assert_close_dp(split_cfg%outer_plasma%max_gap_ratio, 4.0_dp, 1.0e-15_dp, 'split gap limit mismatch')
  call assert_close_dp( &
    split_cfg%outer_plasma%max_local_charge_ratio, 6.0_dp, 1.0e-15_dp, 'split local-charge limit mismatch' &
    )
  call assert_true( &
    trim(split_cfg%outer_plasma%photoelectron_density_model) == 'none', &
    'default photoelectron density model mismatch' &
    )
  call assert_true( &
    trim(split_cfg%outer_plasma%return_model) == 'kinetic_1d_profile_return', &
    'split return model mismatch' &
    )
  call assert_true( &
    trim(split_cfg%coupling%particle_transfer_mode) == 'electrostatic_1d_instant_return', &
    'split transfer mode mismatch' &
    )
  call assert_true(trim(split_cfg%coupling%steady_start_mode) == 'none', 'split default steady-start mode mismatch')
  call assert_equal_i32(split_cfg%coupling%steady_start_mesh_id, 1_i32, 'split default steady-start mesh ID mismatch')
  call assert_close_dp( &
    split_cfg%coupling%field_evolution_timescale, 2.0_dp, 1.0e-15_dp, 'split field timescale mismatch' &
    )
  call test_end()

  call test_begin('external_boundary_facade_lowering')
  call default_app_config(external_boundary_cfg)
  call load_app_config(external_boundary_cfg_path, external_boundary_cfg)
  call assert_true( &
    trim(external_boundary_cfg%outer_plasma%model) == 'kinetic_1d', &
    'facade outer field model mismatch' &
    )
  call assert_close_dp( &
    external_boundary_cfg%outer_plasma%interface_z, 1.0_dp, 1.0e-15_dp, &
    'facade interface location mismatch' &
    )
  call assert_true( &
    trim(external_boundary_cfg%outer_plasma%return_model) == 'kinetic_1d_profile_return', &
    'facade return model mismatch' &
    )
  call assert_true( &
    trim(external_boundary_cfg%coupling%particle_transfer_mode) == 'electrostatic_1d_instant_return', &
    'facade transfer mode mismatch' &
    )
  call assert_true(trim(external_boundary_cfg%coupling%update_mode) == 'explicit', 'facade update mode mismatch')
  call assert_true(.not. external_boundary_cfg%coupling%outer_queue_enabled, 'facade queue mismatch')
  call assert_true( &
    trim(external_boundary_cfg%sim%reservoir_potential_model) == 'none', &
    'profile-owned facade inflow must disable the scalar correction' &
    )
  call assert_true( &
    trim(external_boundary_cfg%sim%open_boundary_model) == 'potential_barrier', &
    'facade ordinary-open model mismatch' &
    )
  call assert_close_dp( &
    external_boundary_cfg%coupling%field_evolution_timescale, 2.0_dp, 1.0e-15_dp, &
    'facade field timescale mismatch' &
    )
  call default_app_config(none_boundary_cfg)
  call load_app_config(none_boundary_cfg_path, none_boundary_cfg)
  call assert_true(trim(none_boundary_cfg%outer_plasma%model) == 'none', 'none facade field mismatch')
  call assert_close_dp( &
    none_boundary_cfg%outer_plasma%interface_z, 0.0_dp, 0.0_dp, &
    'none facade must preserve the runtime interface default' &
    )
  call assert_true(trim(none_boundary_cfg%outer_plasma%return_model) == 'none', 'none facade return mismatch')
  call assert_true(trim(none_boundary_cfg%coupling%particle_transfer_mode) == 'none', 'none facade transfer mismatch')
  call assert_true( &
    model_fingerprint(external_boundary_cfg) == model_fingerprint(split_cfg), &
    'facade and raw external-boundary syntax must have identical model fingerprints' &
    )
  call test_end()

  call test_begin('external_boundary_facade_mode_matrix')
  call assert_external_boundary_mode_matrix()
  call test_end()

  call test_begin('external_boundary_facade_rejects_raw_selector')
  call assert_mixed_external_boundary_owner_rejected()
  call test_end()

  call test_begin('external_boundary_facade_rejects_noop_options')
  call assert_external_boundary_noop_options_rejected()
  call test_end()

  call test_begin('triangle_panel_config')
  call default_app_config(panel_cfg)
  call load_app_config(panel_cfg_path, panel_cfg)
  call assert_true(trim(panel_cfg%panel%kernel_id) == 'triangle_p0_exact_direct', 'panel kernel id mismatch')
  call assert_true(trim(panel_cfg%templates(1)%surface_side_policy) == 'normal_plus', 'template side policy mismatch')
  call test_end()

  call test_begin('triangle_panel_treecode_config')
  call default_app_config(panel_tree_cfg)
  call load_app_config(panel_tree_cfg_path, panel_tree_cfg)
  call assert_true(trim(panel_tree_cfg%sim%field_solver) == 'treecode', 'panel tree solver mismatch')
  call assert_true( &
    trim(panel_tree_cfg%panel%kernel_id) == 'triangle_p0_exact_tree_near', &
    'panel tree kernel id mismatch' &
    )
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

  call test_begin('removed_root_oracle_config')
  call assert_removed_root_oracle_rejected()
  call test_end()

  call test_begin('zhao_physics_inputs')
  call default_app_config(zhao_cfg)
  call load_app_config(zhao_cfg_path, zhao_cfg)
  call assert_close_dp(zhao_cfg%sim%sheath_alpha_deg, 42.0d0, 1.0d-12, 'Zhao alpha mismatch')
  call assert_close_dp( &
    zhao_cfg%sim%sheath_photoelectron_ref_density_cm3, 48.0d0, 1.0d-12, &
    'Zhao photoelectron density mismatch' &
    )
  call assert_true(trim(zhao_cfg%sim%sheath_electron_drift_mode) == 'full', 'Zhao electron drift mode mismatch')
  call assert_true(trim(zhao_cfg%sim%sheath_ion_drift_mode) == 'normal', 'Zhao ion drift mode mismatch')
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

  call test_begin('tracked_1d_photo_requires_countercharge')
  call assert_missing_1d_photo_deposit_rejected()
  call test_end()

  call test_begin('ambient_linear_debye_derives_implicit_mean_only_with_photo')
  call default_app_config(implicit_mean_cfg)
  call load_app_config(implicit_mean_cfg_path, implicit_mean_cfg)
  call assert_true( &
    trim(implicit_mean_cfg%outer_plasma%kinetic_closure) == 'ambient_linear_debye', &
    'photo fixture kinetic closure mismatch' &
    )
  call assert_true( &
    trim(implicit_mean_cfg%coupling%update_mode) == 'implicit_mean', &
    'enabled photo_raycast must derive implicit_mean coupling' &
    )
  call default_app_config(ambient_only_cfg)
  call load_app_config(ambient_only_cfg_path, ambient_only_cfg)
  call assert_true( &
    trim(ambient_only_cfg%outer_plasma%kinetic_closure) == 'ambient_linear_debye', &
    'ambient-only fixture kinetic closure mismatch' &
    )
  call assert_true( &
    trim(ambient_only_cfg%coupling%update_mode) == 'explicit', &
    'ambient_linear_debye without photo_raycast must retain explicit coupling' &
    )
  call default_app_config(positive_photo_cfg)
  call load_app_config(positive_photo_cfg_path, positive_photo_cfg)
  call assert_true( &
    trim(positive_photo_cfg%coupling%update_mode) == 'explicit', &
    'positive photo_raycast must not derive the negative-photo implicit mean coupling' &
    )
  call test_end()

  call test_begin('legacy_coupling_rejects_authored_implicit_mean')
  call assert_authored_implicit_mean_rejected()
  call test_end()

  call test_begin('implicit_mean_accepts_non_template_source_provenance')
  call default_app_config(non_template_mean_cfg)
  call load_app_config(non_template_mean_cfg_path, non_template_mean_cfg)
  call assert_true( &
    trim(non_template_mean_cfg%mesh_mode) == 'obj', &
    'non-template implicit-mean fixture lost its mesh mode' &
    )
  call assert_true( &
    trim(non_template_mean_cfg%coupling%update_mode) == 'implicit_mean', &
    'source-provenance coupling must not require a generated support plane' &
    )
  call test_end()

  call delete_file_if_exists(cfg_path)
  call delete_file_if_exists(photo_cfg_path)
  call delete_file_if_exists(large_cfg_path)
  call delete_file_if_exists(periodic_cfg_path)
  call delete_file_if_exists(periodic_oracle_cfg_path)
  call delete_file_if_exists(periodic_oracle_output_path)
  call delete_file_if_exists(zhao_cfg_path)
  call delete_file_if_exists(high_level_cfg_path)
  call delete_file_if_exists(multiline_cfg_path)
  call delete_file_if_exists(toml_syntax_cfg_path)
  call delete_file_if_exists(panel_cfg_path)
  call delete_file_if_exists(panel_tree_cfg_path)
  call delete_file_if_exists(split_cfg_path)
  call delete_file_if_exists(external_boundary_cfg_path)
  call delete_file_if_exists(mixed_boundary_cfg_path)
  call delete_file_if_exists(mixed_boundary_output_path)
  call delete_file_if_exists(noop_boundary_cfg_path)
  call delete_file_if_exists(noop_boundary_output_path)
  call delete_file_if_exists(none_boundary_cfg_path)
  call delete_file_if_exists(implicit_mean_cfg_path)
  call delete_file_if_exists(ambient_only_cfg_path)
  call delete_file_if_exists(positive_photo_cfg_path)
  call delete_file_if_exists(missing_deposit_cfg_path)
  call delete_file_if_exists(missing_deposit_output_path)
  call delete_file_if_exists(non_template_mean_cfg_path)
  call delete_file_if_exists(authored_implicit_cfg_path)
  call delete_file_if_exists(authored_implicit_output_path)

  call test_summary()

contains

  subroutine assert_authored_implicit_mean_rejected()
    character(len=1024) :: executable_path, line
    character(len=4096) :: command
    integer :: child_exit_status, child_cmd_status, u, ios
    logical :: saw_requirement

    call get_command_argument(0, executable_path)
    command = '"'//trim(executable_path)//'" probe_authored_implicit_mean > "'// &
              authored_implicit_output_path//'" 2>&1'
    call execute_command_line( &
      trim(command), wait=.true., exitstat=child_exit_status, cmdstat=child_cmd_status &
      )
    call assert_equal_i32(int(child_cmd_status, i32), 0_i32, 'authored implicit_mean probe command status mismatch')
    call assert_true(child_exit_status /= 0, 'user-authored coupling.update_mode=implicit_mean must fail')
    saw_requirement = .false.
    open (newunit=u, file=authored_implicit_output_path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to read authored implicit_mean probe output'
    do
      read (u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      saw_requirement = saw_requirement .or. &
                        index(line, 'coupling.update_mode accepts only "explicit" in user input') > 0
    end do
    close (u)
    call assert_true(saw_requirement, 'authored implicit_mean probe must report the explicit-only input contract')
    call delete_file_if_exists(authored_implicit_cfg_path)
    call delete_file_if_exists(authored_implicit_output_path)
  end subroutine assert_authored_implicit_mean_rejected

  subroutine assert_removed_root_oracle_rejected()
    character(len=1024) :: executable_path, line
    character(len=4096) :: command
    integer :: child_exit_status, child_cmd_status, u, ios
    logical :: saw_removal

    call get_command_argument(0, executable_path)
    command = '"'//trim(executable_path)//'" probe_removed_root_oracle > "'// &
              periodic_oracle_output_path//'" 2>&1'
    call execute_command_line( &
      trim(command), wait=.true., exitstat=child_exit_status, cmdstat=child_cmd_status &
      )
    call assert_equal_i32(int(child_cmd_status, i32), 0_i32, 'removed root oracle probe command status mismatch')
    call assert_true(child_exit_status /= 0, 'removed root oracle configuration must fail')
    saw_removal = .false.
    open (newunit=u, file=periodic_oracle_output_path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to read removed root oracle probe output'
    do
      read (u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      saw_removal = saw_removal .or. index(line, '"m2l_root_oracle" was removed; use "cached_kneq0"') > 0
    end do
    close (u)
    call assert_true(saw_removal, 'removed root oracle probe must recommend cached_kneq0')
    call delete_file_if_exists(periodic_oracle_cfg_path)
    call delete_file_if_exists(periodic_oracle_output_path)
  end subroutine assert_removed_root_oracle_rejected

  subroutine assert_missing_1d_photo_deposit_rejected()
    character(len=1024) :: executable_path, line
    character(len=4096) :: command
    integer :: child_exit_status, child_cmd_status, u, ios
    logical :: saw_requirement

    call get_command_argument(0, executable_path)
    command = '"'//trim(executable_path)//'" probe_missing_1d_photo_deposit > "'// &
              missing_deposit_output_path//'" 2>&1'
    call execute_command_line( &
      trim(command), wait=.true., exitstat=child_exit_status, cmdstat=child_cmd_status &
      )
    call assert_equal_i32(int(child_cmd_status, i32), 0_i32, 'missing deposit probe command status mismatch')
    call assert_true(child_exit_status /= 0, 'tracked 1D photoelectron transfer without countercharge must fail')
    saw_requirement = .false.
    open (newunit=u, file=missing_deposit_output_path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to read missing deposit probe output'
    do
      read (u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      saw_requirement = saw_requirement .or. &
                        index(line, 'tracked outer transfer requires photo_raycast '// &
                              'deposit_opposite_charge_on_emit=true') > 0
    end do
    close (u)
    call assert_true(saw_requirement, 'missing deposit probe must report the tracked outer transfer requirement')
    call delete_file_if_exists(missing_deposit_cfg_path)
    call delete_file_if_exists(missing_deposit_output_path)
  end subroutine assert_missing_1d_photo_deposit_rejected

  subroutine assert_external_boundary_mode_matrix()
    type(app_config) :: mode_cfg
    type(app_config_authoring) :: mode_authoring

    call default_app_config(mode_cfg)
    mode_cfg%sim%box_max(3) = 2.0_dp
    call prepare_external_boundary_authoring(mode_authoring, 'kinetic_1d', 'same_batch', 'auto')
    mode_authoring%external_boundary%field%has_kinetic_closure = .true.
    mode_authoring%external_boundary%field%has_photoelectron_density_model = .true.
    mode_authoring%external_boundary%field%photoelectron_density_model = 'kinetic_mean'
    mode_authoring%external_boundary%particles%has_time_guard_option = .true.
    mode_authoring%external_boundary%particles%field_evolution_timescale = 2.0_dp
    call lower_external_boundary_authoring(mode_cfg, mode_authoring)
    call assert_true( &
      trim(mode_cfg%outer_plasma%return_model) == 'kinetic_1d_profile_return', &
      'kinetic same-batch return mapping mismatch' &
      )
    call assert_true( &
      trim(mode_cfg%coupling%particle_transfer_mode) == 'electrostatic_1d_instant_return', &
      'kinetic same-batch transfer mapping mismatch' &
      )
    call assert_true(.not. mode_cfg%coupling%outer_queue_enabled, 'kinetic same-batch queue mismatch')
    call assert_true( &
      trim(mode_cfg%sim%reservoir_potential_model) == 'none', &
      'kinetic profile-owned inflow must disable the scalar inflow correction' &
      )

    call default_app_config(mode_cfg)
    mode_cfg%sim%box_max(3) = 4.0_dp
    call prepare_external_boundary_authoring(mode_authoring, 'kinetic_1d', 'zhao_queue', 'auto')
    mode_authoring%external_boundary%field%kinetic_closure = 'zhao_charge_driven'
    mode_authoring%external_boundary%field%has_kinetic_closure = .true.
    mode_authoring%external_boundary%particles%has_time_guard_option = .true.
    mode_authoring%external_boundary%particles%field_evolution_timescale = 2.0_dp
    call lower_external_boundary_authoring(mode_cfg, mode_authoring)
    call assert_true( &
      trim(mode_cfg%outer_plasma%return_model) == 'kinetic_1d_profile_return', &
      'Zhao queue return mapping mismatch' &
      )
    call assert_true( &
      trim(mode_cfg%coupling%particle_transfer_mode) == 'electrostatic_1d_instant_return', &
      'Zhao queue transfer mapping mismatch' &
      )
    call assert_true(mode_cfg%coupling%outer_queue_enabled, 'Zhao queue flag mapping mismatch')

    call default_app_config(mode_cfg)
    mode_cfg%sim%box_max(3) = 5.0_dp
    call prepare_external_boundary_authoring(mode_authoring, 'kinetic_1d', 'same_batch', 'auto')
    mode_authoring%external_boundary%field%kinetic_closure = 'zhao_charge_driven'
    mode_authoring%external_boundary%field%has_kinetic_closure = .true.
    mode_authoring%external_boundary%particles%steady_start_mode = 'zhao_floating'
    mode_authoring%external_boundary%particles%steady_start_mesh_id = 2_i32
    mode_authoring%external_boundary%particles%has_steady_start_mode = .true.
    mode_authoring%external_boundary%particles%has_steady_start_mesh_id = .true.
    mode_authoring%external_boundary%particles%has_time_guard_option = .true.
    mode_authoring%external_boundary%particles%field_evolution_timescale = 2.0_dp
    call lower_external_boundary_authoring(mode_cfg, mode_authoring)
    call assert_true( &
      trim(mode_cfg%coupling%steady_start_mode) == 'zhao_floating', &
      'Zhao same-batch steady-start mapping mismatch' &
      )
    call assert_equal_i32(mode_cfg%coupling%steady_start_mesh_id, 2_i32, 'Zhao steady-start mesh ID mismatch')

    call default_app_config(mode_cfg)
    mode_cfg%sim%box_max(3) = 6.0_dp
    call prepare_external_boundary_authoring(mode_authoring, 'kinetic_1d', 'local_source', 'source_vdf')
    mode_authoring%external_boundary%particles%has_outer_update_stride = .true.
    mode_authoring%external_boundary%particles%outer_update_stride = 3_i32
    call lower_external_boundary_authoring(mode_cfg, mode_authoring)
    call assert_equal_i32(mode_cfg%coupling%outer_update_stride, 3_i32, 'kinetic local-source update stride mismatch')
    call assert_true( &
      trim(mode_cfg%sim%reservoir_potential_model) == 'none', &
      'source-vdf inflow must not enable the scalar inflow correction' &
      )
  end subroutine assert_external_boundary_mode_matrix

  subroutine prepare_external_boundary_authoring(authoring, field_model, particle_mode, inflow_model)
    type(app_config_authoring), intent(out) :: authoring
    character(len=*), intent(in) :: field_model, particle_mode, inflow_model

    authoring%external_boundary%present = .true.
    authoring%external_boundary%field%present = .true.
    authoring%external_boundary%field%has_model = .true.
    authoring%external_boundary%field%model = trim(field_model)
    authoring%external_boundary%particles%present = .true.
    authoring%external_boundary%particles%has_mode = .true.
    authoring%external_boundary%particles%mode = trim(particle_mode)
    authoring%external_boundary%particles%inflow_model = trim(inflow_model)
    if (trim(field_model) /= 'none') authoring%sim%has_box_max = .true.
  end subroutine prepare_external_boundary_authoring

  subroutine assert_mixed_external_boundary_owner_rejected()
    call assert_external_boundary_owner_rejected( &
      'outer_plasma', '[external_boundary] cannot be combined with legacy [outer_plasma]' &
      )
    call assert_external_boundary_owner_rejected( &
      'coupling', '[external_boundary] cannot be combined with legacy [coupling]' &
      )
    call assert_external_boundary_owner_rejected( &
      'reservoir', '[external_boundary] cannot be combined with sim.reservoir_potential_model' &
      )
    call assert_external_boundary_owner_rejected( &
      'open', '[external_boundary] cannot be combined with sim.open_boundary_model' &
      )
  end subroutine assert_mixed_external_boundary_owner_rejected

  subroutine assert_external_boundary_owner_rejected(owner, expected_message)
    character(len=*), intent(in) :: owner, expected_message
    character(len=1024) :: executable_path, line
    character(len=4096) :: command
    integer :: child_exit_status, child_cmd_status, u, ios
    logical :: saw_conflict

    call get_command_argument(0, executable_path)
    command = '"'//trim(executable_path)//'" probe_mixed_external_boundary_owner '//trim(owner)//' > "'// &
              mixed_boundary_output_path//'" 2>&1'
    call execute_command_line( &
      trim(command), wait=.true., exitstat=child_exit_status, cmdstat=child_cmd_status &
      )
    call assert_equal_i32(int(child_cmd_status, i32), 0_i32, 'mixed owner probe command status mismatch')
    call assert_true(child_exit_status /= 0, 'facade combined with a raw selector must fail')
    saw_conflict = .false.
    open (newunit=u, file=mixed_boundary_output_path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to read mixed owner probe output'
    do
      read (u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      saw_conflict = saw_conflict .or. index(line, trim(expected_message)) > 0
    end do
    close (u)
    call assert_true(saw_conflict, 'mixed owner probe must report the explicit owner conflict')
    call delete_file_if_exists(mixed_boundary_cfg_path)
    call delete_file_if_exists(mixed_boundary_output_path)
  end subroutine assert_external_boundary_owner_rejected

  subroutine assert_external_boundary_noop_options_rejected()
    call assert_external_boundary_probe_rejected( &
      'probe_none_field_option', 'model="none" does not accept additional field options' &
      )
    call assert_external_boundary_probe_rejected( &
      'probe_none_local_coupling_option', 'mode="local_source" does not accept coupling options' &
      )
    call assert_external_boundary_probe_rejected( &
      'probe_active_field_without_box', 'requires explicit sim.box_max or sim.box_origin/box_size' &
      )
    call assert_external_boundary_probe_rejected( &
      'probe_kinetic_zhao_option', 'Zhao field options require kinetic_1d' &
      )
    call assert_external_boundary_probe_rejected( &
      'probe_local_time_guard', 'time-guard options require particles.mode="same_batch" or "zhao_queue"' &
      )
    call assert_external_boundary_probe_rejected( &
      'probe_queue_update_stride', 'outer_update_stride requires a kinetic_1d field and a non-queue particle mode' &
      )
    call assert_external_boundary_probe_rejected( &
      'probe_kinetic_steady_start', 'Steady-start options require Zhao kinetic_1d' &
      )
  end subroutine assert_external_boundary_noop_options_rejected

  subroutine assert_external_boundary_probe_rejected(probe_mode, expected_message)
    character(len=*), intent(in) :: probe_mode, expected_message
    character(len=1024) :: executable_path, line
    character(len=4096) :: command
    integer :: child_exit_status, child_cmd_status, u, ios
    logical :: saw_expected_message

    call get_command_argument(0, executable_path)
    command = '"'//trim(executable_path)//'" '//trim(probe_mode)//' > "'// &
              noop_boundary_output_path//'" 2>&1'
    call execute_command_line( &
      trim(command), wait=.true., exitstat=child_exit_status, cmdstat=child_cmd_status &
      )
    call assert_equal_i32(int(child_cmd_status, i32), 0_i32, 'no-op option probe command status mismatch')
    call assert_true(child_exit_status /= 0, 'facade no-op option probe must fail')
    saw_expected_message = .false.
    open (newunit=u, file=noop_boundary_output_path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'failed to read no-op option probe output'
    do
      read (u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      saw_expected_message = saw_expected_message .or. index(line, trim(expected_message)) > 0
    end do
    close (u)
    call assert_true(saw_expected_message, 'no-op option probe must report the facade applicability conflict')
    call delete_file_if_exists(noop_boundary_cfg_path)
    call delete_file_if_exists(noop_boundary_output_path)
  end subroutine assert_external_boundary_probe_rejected

  !> ambient linear Debye のmean更新導出を検証する設定を書き出す。
  !! @param[in] path 書き出し先TOMLファイルパス。
  !! @param[in] include_photo photo_raycast speciesを含める場合はtrue。
  !! @param[in] deposit_countercharge 光電子放出時の逆電荷をdepositする場合はtrue。
  !! @param[in] mesh_mode optionalなmesh mode。省略時はtemplate。
  !! @param[in] photo_charge optionalなphoto_raycast粒子電荷。省略時は電子電荷。
  subroutine write_ambient_linear_debye_config_fixture( &
    path, include_photo, deposit_countercharge, mesh_mode, photo_charge &
    )
    character(len=*), intent(in) :: path
    logical, intent(in) :: include_photo, deposit_countercharge
    character(len=*), intent(in), optional :: mesh_mode
    real(dp), intent(in), optional :: photo_charge
    integer :: u, ios
    character(len=16) :: resolved_mesh_mode
    real(dp) :: resolved_photo_charge

    resolved_mesh_mode = 'template'
    if (present(mesh_mode)) resolved_mesh_mode = trim(mesh_mode)
    resolved_photo_charge = -1.602176634e-19_dp
    if (present(photo_charge)) resolved_photo_charge = photo_charge

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open ambient linear Debye config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'dt = 1.0e-9'
    write (u, '(a)') 'batch_duration = 1.0e-8'
    write (u, '(a)') 'batch_count = 1'
    write (u, '(a)') 'max_step = 4'
    write (u, '(a)') 'field_solver = "direct"'
    write (u, '(a)') 'field_bc_mode = "periodic2"'
    write (u, '(a)') 'use_box = true'
    write (u, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (u, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    write (u, '(a)') 'bc_x_low = "periodic"'
    write (u, '(a)') 'bc_x_high = "periodic"'
    write (u, '(a)') 'bc_y_low = "periodic"'
    write (u, '(a)') 'bc_y_high = "periodic"'
    write (u, '(a)') 'bc_z_low = "open"'
    write (u, '(a)') 'bc_z_high = "open"'
    write (u, '(a)') ''
    write (u, '(a)') '[periodic2]'
    write (u, '(a)') 'nonzero_mode_backend = "panel_spectral_reference"'
    write (u, '(a)') 'zero_mode_policy = "exclude_k0"'
    write (u, '(a)') 'lower_boundary_model = "symmetric_vacuum"'
    write (u, '(a)') 'reference_mode_layers = 2'
    write (u, '(a)') 'panel_quadrature_order = 4'
    write (u, '(a)') 'interface_sample_n = 2'
    write (u, '(a)') 'interface_phi_tolerance = 1.0e-3'
    write (u, '(a)') 'interface_field_tolerance = 1.0e-3'
    write (u, '(a)') ''
    write (u, '(a)') '[external_boundary.field]'
    write (u, '(a)') 'model = "kinetic_1d"'
    write (u, '(a)') 'kinetic_closure = "ambient_linear_debye"'
    write (u, '(a)') 'debye_length = 0.2'
    write (u, '(a)') 'thermal_voltage = 2.0'
    write (u, '(a)') ''
    write (u, '(a)') '[external_boundary.particles]'
    write (u, '(a)') 'mode = "same_batch"'
    write (u, '(a)') 'field_evolution_timescale = 1.0'
    write (u, '(a)') ''
    write (u, '(a)') '[particles]'
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'species_key = "ambient_electron"'
    write (u, '(a)') 'source_mode = "reservoir_face"'
    write (u, '(a)') 'number_density_m3 = 1.0e6'
    write (u, '(a)') 'target_macro_particles_per_batch = 1'
    write (u, '(a)') 'q_particle = -1.602176634e-19'
    write (u, '(a)') 'm_particle = 9.1093837139e-31'
    write (u, '(a)') 'inject_face = "z_high"'
    write (u, '(a)') 'pos_low = [0.0, 0.0, 1.0]'
    write (u, '(a)') 'pos_high = [1.0, 1.0, 1.0]'
    write (u, '(a)') 'drift_velocity = [0.0, 0.0, -1.0e5]'
    write (u, '(a)') 'temperature_ev = 2.0'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'species_key = "ambient_proton"'
    write (u, '(a)') 'source_mode = "reservoir_face"'
    write (u, '(a)') 'number_density_m3 = 1.0e6'
    write (u, '(a)') 'target_macro_particles_per_batch = 1'
    write (u, '(a)') 'q_particle = 1.602176634e-19'
    write (u, '(a)') 'm_particle = 1.67262192595e-27'
    write (u, '(a)') 'inject_face = "z_high"'
    write (u, '(a)') 'pos_low = [0.0, 0.0, 1.0]'
    write (u, '(a)') 'pos_high = [1.0, 1.0, 1.0]'
    write (u, '(a)') 'drift_velocity = [0.0, 0.0, -4.0e4]'
    write (u, '(a)') 'temperature_ev = 0.1'
    if (include_photo) then
      write (u, '(a)') ''
      write (u, '(a)') '[[particles.species]]'
      write (u, '(a)') 'species_key = "photoelectron"'
      write (u, '(a)') 'source_mode = "photo_raycast"'
      write (u, '(a)') 'emit_current_density_a_m2 = 1.0e-3'
      write (u, '(a)') 'rays_per_batch = 4'
      if (deposit_countercharge) write (u, '(a)') 'deposit_opposite_charge_on_emit = true'
      write (u, '(a,es24.16)') 'q_particle = ', resolved_photo_charge
      write (u, '(a)') 'm_particle = 9.1093837139e-31'
      write (u, '(a)') 'temperature_k = 0.0'
      write (u, '(a)') 'inject_face = "z_high"'
      write (u, '(a)') 'pos_low = [0.0, 0.0, 1.0]'
      write (u, '(a)') 'pos_high = [1.0, 1.0, 1.0]'
      write (u, '(a)') 'ray_direction = [0.0, 0.0, -1.0]'
    end if
    write (u, '(a)') ''
    write (u, '(a)') '[mesh]'
    write (u, '(a)') 'mode = "'//trim(resolved_mesh_mode)//'"'
    if (trim(resolved_mesh_mode) == 'template') then
      write (u, '(a)') '[[mesh.templates]]'
      write (u, '(a)') 'enabled = true'
      write (u, '(a)') 'kind = "plane"'
      write (u, '(a)') 'surface_model = "insulator"'
      write (u, '(a)') 'surface_side = "normal_plus"'
      write (u, '(a)') 'center = [0.5, 0.5, 0.25]'
      write (u, '(a)') 'size_x = 1.0'
      write (u, '(a)') 'size_y = 1.0'
      write (u, '(a)') 'nx = 1'
      write (u, '(a)') 'ny = 1'
    else
      write (u, '(a)') 'obj_path = "unused.obj"'
      write (u, '(a)') 'surface_side = "normal_plus"'
    end if
    write (u, '(a)') ''
    write (u, '(a)') '[output]'
    write (u, '(a)') 'write_files = false'

    close (u)
  end subroutine write_ambient_linear_debye_config_fixture

  subroutine write_authored_implicit_mean_config_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open authored implicit_mean config fixture'
    write (u, '(a)') '[coupling]'
    write (u, '(a)') 'update_mode = "implicit_mean"'
    close (u)
  end subroutine write_authored_implicit_mean_config_fixture

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
    write (u, '(a)') 'multiple_box_events_policy = "soft_discard"'
    write (u, '(a)') 'multiple_box_events_soft_discard_count_limit = 100000'
    write (u, '(a)') 'multiple_box_events_soft_discard_abs_charge_limit = 1.0e-9'
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

  !> 削除済み `m2l_root_oracle` を指定する periodic2 設定を書き出す。
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

  !> Zhao kinetic closure が共有する物理入力キーの受理を確認する設定を書き出す。
  subroutine write_zhao_physics_config_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open Zhao physics config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'batch_count = 1'
    write (u, '(a)') 'sheath_alpha_deg = 42.0'
    write (u, '(a)') 'sheath_photoelectron_ref_density_cm3 = 48.0'
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
  end subroutine write_zhao_physics_config_fixture

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

  subroutine write_panel_config_fixture(path, solver)
    character(len=*), intent(in) :: path, solver
    integer :: u, ios

    open (newunit=u, file=path, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open panel config fixture'
    write (u, '(a)') '[sim]'
    write (u, '(a)') 'batch_count = 1'
    write (u, '(a)') 'dt = 1.0e-9'
    write (u, '(a)') 'max_step = 2'
    write (u, '(a)') 'field_solver = "'//trim(solver)//'"'
    write (u, '(a)') 'field_bc_mode = "free"'
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
    write (u, '(a)') 'use_box = true'
    write (u, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (u, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    write (u, '(a)') 'bc_x_low = "periodic"'
    write (u, '(a)') 'bc_x_high = "periodic"'
    write (u, '(a)') 'bc_y_low = "periodic"'
    write (u, '(a)') 'bc_y_high = "periodic"'
    write (u, '(a)') 'bc_z_low = "open"'
    write (u, '(a)') 'bc_z_high = "open"'
    write (u, '(a)') 'open_boundary_model = "potential_barrier"'
    write (u, '(a)') '[periodic2]'
    write (u, '(a)') 'nonzero_mode_backend = "panel_spectral_reference"'
    write (u, '(a)') 'zero_mode_policy = "exclude_k0"'
    write (u, '(a)') 'lower_boundary_model = "symmetric_vacuum"'
    write (u, '(a)') 'reference_mode_layers = 5'
    write (u, '(a)') 'panel_quadrature_order = 16'
    write (u, '(a)') 'interface_sample_n = 7'
    write (u, '(a)') '[outer_plasma]'
    write (u, '(a)') 'model = "kinetic_1d"'
    write (u, '(a)') 'return_model = "kinetic_1d_profile_return"'
    write (u, '(a)') 'interface_z = 1.0'
    write (u, '(a)') 'debye_length = 0.2'
    write (u, '(a)') 'thermal_voltage = 10.0'
    write (u, '(a)') 'max_gap_ratio = 4.0'
    write (u, '(a)') 'max_local_charge_ratio = 6.0'
    write (u, '(a)') '[coupling]'
    write (u, '(a)') 'update_mode = "explicit"'
    write (u, '(a)') 'particle_transfer_mode = "electrostatic_1d_instant_return"'
    write (u, '(a)') 'outer_update_stride = 1'
    write (u, '(a)') 'field_evolution_timescale = 2.0'
    write (u, '(a)') 'max_frozen_field_ratio = 0.1'
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

  subroutine write_external_boundary_config_fixture(path, conflict_owner, omit_box)
    character(len=*), intent(in) :: path
    character(len=*), intent(in), optional :: conflict_owner
    logical, intent(in), optional :: omit_box
    integer :: u, ios
    logical :: omit_explicit_box
    character(len=32) :: owner

    open (newunit=u, file=path, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open external-boundary config fixture'
    owner = ''
    if (present(conflict_owner)) owner = trim(conflict_owner)
    omit_explicit_box = .false.
    if (present(omit_box)) omit_explicit_box = omit_box
    write (u, '(a)') '[sim]'
    write (u, '(a)') 'batch_count = 1'
    write (u, '(a)') 'dt = 1.0e-9'
    write (u, '(a)') 'max_step = 2'
    write (u, '(a)') 'field_solver = "direct"'
    write (u, '(a)') 'field_bc_mode = "periodic2"'
    write (u, '(a)') 'use_box = true'
    if (.not. omit_explicit_box) then
      write (u, '(a)') 'box_origin = [0.0, 0.0, 0.0]'
      write (u, '(a)') 'box_size = [1.0, 1.0, 1.0]'
    end if
    write (u, '(a)') 'bc_x_low = "periodic"'
    write (u, '(a)') 'bc_x_high = "periodic"'
    write (u, '(a)') 'bc_y_low = "periodic"'
    write (u, '(a)') 'bc_y_high = "periodic"'
    write (u, '(a)') 'bc_z_low = "open"'
    write (u, '(a)') 'bc_z_high = "open"'
    if (trim(owner) == 'reservoir') write (u, '(a)') 'reservoir_potential_model = "none"'
    if (trim(owner) == 'open') write (u, '(a)') 'open_boundary_model = "escape"'
    write (u, '(a)') '[periodic2]'
    write (u, '(a)') 'nonzero_mode_backend = "panel_spectral_reference"'
    write (u, '(a)') 'zero_mode_policy = "exclude_k0"'
    write (u, '(a)') 'lower_boundary_model = "symmetric_vacuum"'
    write (u, '(a)') 'reference_mode_layers = 5'
    write (u, '(a)') 'panel_quadrature_order = 16'
    write (u, '(a)') 'interface_sample_n = 7'
    if (trim(owner) == 'outer_plasma') write (u, '(a)') '[outer_plasma]'
    if (trim(owner) == 'coupling') write (u, '(a)') '[coupling]'
    write (u, '(a)') '[external_boundary.field]'
    write (u, '(a)') 'model = "kinetic_1d"'
    write (u, '(a)') 'debye_length = 0.2'
    write (u, '(a)') 'thermal_voltage = 10.0'
    write (u, '(a)') 'max_gap_ratio = 4.0'
    write (u, '(a)') 'max_local_charge_ratio = 6.0'
    write (u, '(a)') '[external_boundary.particles]'
    write (u, '(a)') 'mode = "same_batch"'
    write (u, '(a)') 'inflow_model = "auto"'
    write (u, '(a)') 'outer_update_stride = 1'
    write (u, '(a)') 'field_evolution_timescale = 2.0'
    write (u, '(a)') 'max_frozen_field_ratio = 0.1'
    write (u, '(a)') '[external_boundary.ordinary_open]'
    write (u, '(a)') 'model = "potential_barrier"'
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
  end subroutine write_external_boundary_config_fixture

  subroutine write_external_boundary_noop_fixture(path, conflict_kind)
    character(len=*), intent(in) :: path, conflict_kind
    integer :: u, ios
    character(len=32) :: field_model, particle_mode

    field_model = 'none'
    particle_mode = 'local_source'
    select case (trim(conflict_kind))
    case ('field', 'coupling', 'none')
      continue
    case ('kinetic_zhao', 'kinetic_steady')
      field_model = 'kinetic_1d'
      particle_mode = 'same_batch'
    case ('queue_stride')
      field_model = 'kinetic_1d'
      particle_mode = 'zhao_queue'
    case ('local_time_guard')
      field_model = 'kinetic_1d'
    case default
      error stop 'unknown external-boundary no-op fixture kind'
    end select

    open (newunit=u, file=path, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open external-boundary no-op fixture'
    if (trim(field_model) /= 'none') then
      write (u, '(a)') '[sim]'
      write (u, '(a)') 'box_origin = [0.0, 0.0, 0.0]'
      write (u, '(a)') 'box_size = [1.0, 1.0, 1.0]'
    end if
    write (u, '(a)') '[external_boundary.field]'
    write (u, '(a)') 'model = "'//trim(field_model)//'"'
    if (trim(field_model) /= 'none') then
      write (u, '(a)') 'debye_length = 1.0'
      write (u, '(a)') 'thermal_voltage = 1.0'
    end if
    if (trim(conflict_kind) == 'field') write (u, '(a)') 'debye_length = 1.0'
    if (trim(conflict_kind) == 'kinetic_zhao') write (u, '(a)') 'zhao_branch = "auto"'
    if (trim(conflict_kind) == 'queue_stride') write (u, '(a)') 'kinetic_closure = "zhao_charge_driven"'
    write (u, '(a)') '[external_boundary.particles]'
    write (u, '(a)') 'mode = "'//trim(particle_mode)//'"'
    if (trim(conflict_kind) == 'coupling') write (u, '(a)') 'outer_update_stride = 1'
    if (trim(conflict_kind) == 'local_time_guard') write (u, '(a)') 'field_evolution_timescale = 2.0'
    if (trim(conflict_kind) == 'queue_stride') write (u, '(a)') 'outer_update_stride = 1'
    if (trim(conflict_kind) == 'kinetic_steady') write (u, '(a)') 'steady_start_mode = "zhao_floating"'
    write (u, '(a)') '[particles]'
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'npcls_per_step = 1'
    close (u)
  end subroutine write_external_boundary_noop_fixture

end program test_app_config_parser

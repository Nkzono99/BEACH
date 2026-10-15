!> 設定の正規化と領域別 preflight の実行順を管理する。
submodule(bem_app_config_parser) bem_app_config_parser_finalize
  implicit none
contains

  module procedure finalize_loaded_config

  call lower_boundary_authoring(cfg, authoring)
  call normalize_high_level_config(cfg, authoring)
  call normalize_legacy_physics_config(cfg%sim, cfg%field, cfg%periodic2, cfg%panel)
  call apply_physics_authoring(cfg, authoring)
  cfg%sim%field_solver = lower_ascii(trim(cfg%sim%field_solver))
  cfg%sim%field_normalization = lower_ascii(trim(cfg%sim%field_normalization))
  cfg%sim%field_bc_mode = lower_ascii(trim(cfg%sim%field_bc_mode))
  cfg%sim%field_periodic_far_correction = lower_ascii(trim(cfg%sim%field_periodic_far_correction))

  call validate_simulation_config(cfg)
  call validate_particle_species_config(cfg)
  call validate_surface_current_model_config( &
    cfg, authoring%periodic2%present .and. &
    authoring%periodic2%has_nonzero_mode_backend .and. &
    authoring%periodic2%has_zero_mode_policy .and. &
    authoring%periodic2%has_lower_boundary_model &
    )
  call validate_field_config(cfg)
  end procedure finalize_loaded_config

  module procedure apply_physics_authoring

  if (authoring%periodic2%present) then
    cfg%periodic2%nonzero_mode_backend = authoring%periodic2%nonzero_mode_backend
    cfg%periodic2%zero_mode_policy = authoring%periodic2%zero_mode_policy
    cfg%periodic2%lower_boundary_model = authoring%periodic2%lower_boundary_model
    cfg%periodic2%reference_mode_layers = authoring%periodic2%reference_mode_layers
    cfg%periodic2%panel_quadrature_order = authoring%periodic2%panel_quadrature_order
    cfg%periodic2%max_nonzero_mode_potential_step = &
      authoring%periodic2%max_nonzero_mode_potential_step
  end if
  end procedure apply_physics_authoring

  module procedure validate_field_config
  integer(i32) :: status
  character(len=256) :: message
  type(field_physics_config) :: field_config
  type(panel_kernel_config) :: panel_config
  call derive_field_panel_config(cfg%sim, field_config, panel_config)
  call validate_active_physics_config(cfg%sim, field_config, cfg%periodic2, panel_config, status, message)
  if (status /= physics_config_ok) call stop_config_error(message)
  if (config_uses_surface_model(cfg, 'dielectric')) error stop 'dielectric surface model is not implemented.'
  call validate_field_parameters(cfg)
  call validate_adaptive_batch_config(cfg)
  call validate_mesh_config(cfg)
  end procedure validate_field_config

  !> solver 固有パラメータは、その solver が有効なときだけ検査する。
  subroutine validate_field_parameters(cfg)
    type(app_config), intent(in) :: cfg

    if (trim(cfg%sim%field_solver) /= 'direct') then
      if (.not. ieee_is_finite(cfg%sim%tree_theta) .or. cfg%sim%tree_theta <= 0.0_dp .or. cfg%sim%tree_theta > 1.0_dp) then
        error stop 'fields.solver.tree.theta must satisfy 0 < theta <= 1 for tree-capable solvers.'
      end if
      if (cfg%sim%tree_leaf_max < 1_i32) error stop 'fields.solver.tree.leaf_max must be >= 1 for tree-capable solvers.'
      if (cfg%sim%tree_min_nelem < 1_i32) error stop 'fields.solver.auto_fmm_min_elements must be >= 1 for tree-capable solvers.'
    end if
    if (trim(cfg%sim%field_bc_mode) /= 'periodic2') return
    if (cfg%sim%field_periodic_image_layers < 0_i32) then
      error stop 'fields.periodic.image_layers must be >= 0 for periodic2.'
    end if
    if (trim(cfg%sim%field_periodic_far_correction) /= 'cached_kneq0') return
    if (cfg%sim%field_periodic_image_layers < 1_i32) then
      error stop 'fields.periodic.image_layers must be >= 1 for cached_kneq0.'
    end if
    if (cfg%sim%field_periodic_ewald_layers < 1_i32) then
      error stop 'fields.periodic.ewald_layers must be >= 1 for cached_kneq0.'
    end if
    if (.not. ieee_is_finite(cfg%sim%field_periodic_ewald_alpha) .or. cfg%sim%field_periodic_ewald_alpha < 0.0_dp) then
      error stop 'fields.periodic.ewald_alpha must be finite and >= 0 for cached_kneq0.'
    end if
    if (.not. ieee_is_finite(cfg%sim%field_periodic_generation_tolerance) .or. &
        cfg%sim%field_periodic_generation_tolerance <= 0.0_dp) then
      error stop 'fields.periodic.generation_tolerance must be finite and > 0 for cached_kneq0.'
    end if
    if (len_trim(cfg%sim%field_periodic_cache_dir) == 0) then
      error stop 'fields.periodic.cache_dir must be a non-empty string for cached_kneq0.'
    end if
  end subroutine validate_field_parameters

  !> adaptive nonzero-mode step は、棄却した trial を巻き戻せる時間スケール付き粒子源だけを許す。
  subroutine validate_adaptive_batch_config(cfg)
    type(app_config), intent(in) :: cfg
    integer :: s
    logical :: flux_driven

    if (cfg%periodic2%max_nonzero_mode_potential_step <= 0.0_dp) return
    if (trim(lower_ascii(cfg%periodic2%nonzero_mode_backend)) /= 'cached_kneq0') then
      error stop 'run.batch.adaptive.max_nonzero_mode_potential_step_v requires fields.periodic.backend="cached_kneq0".'
    end if
    if (.not. ieee_is_finite(cfg%sim%batch_duration) .or. cfg%sim%batch_duration <= 0.0_dp) then
      error stop 'run.batch.adaptive.max_nonzero_mode_potential_step_v requires a finite positive run.batch.duration_s '// &
        'or resolved run.batch.duration_steps.'
    end if
    do s = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(s)%enabled) cycle
      flux_driven = any(cfg%particle_species(s)%boundary_inflow_low == particle_inflow_reservoir) .or. &
                    any(cfg%particle_species(s)%boundary_inflow_high == particle_inflow_reservoir)
      if (trim(cfg%particle_species(s)%source_mode) == 'volume_seed' .and. cfg%particle_species(s)%npcls_per_step > 0_i32) then
        error stop 'run.batch.adaptive.max_nonzero_mode_potential_step_v requires time-scaled reservoir or photo_raycast sources.'
      end if
      if (flux_driven .and. .not. cfg%particle_species(s)%has_target_macro_particles_per_batch) then
        error stop 'adaptive reservoir injection requires '//species_config_path(int(s))// &
          '.sampling.target_macro_particles_per_batch instead of fixed '//species_config_path(int(s))//'.sampling.weight.'
      end if
    end do
  end subroutine validate_adaptive_batch_config

  !> mesh の構築前に決まる幾何条件を検査する。
  subroutine validate_mesh_config(cfg)
    type(app_config), intent(in) :: cfg
    integer :: i

    select case (trim(cfg%mesh_mode))
    case ('obj', 'auto')
      if (trim(cfg%mesh_surface_side_policy) == 'not_specified') then
        error stop 'mesh.obj.surface_side is required for mesh.mode="obj" or "auto".'
      end if
    case ('template')
      continue
    case default
      error stop 'mesh.mode must be "auto", "obj", or "template".'
    end select
    do i = 1, cfg%n_templates
      if (.not. cfg%templates(i)%enabled) cycle
      select case (trim(lower_ascii(cfg%templates(i)%kind)))
      case ('plate_hole', 'plane_hole')
        if (cfg%templates(i)%radius >= 0.5_dp*min(cfg%templates(i)%size_x, cfg%templates(i)%size_y)) then
          error stop 'mesh.templates plate_hole radius must be smaller than half of min(size_x, size_y).'
        end if
      case ('annulus')
        if (cfg%templates(i)%inner_radius >= cfg%templates(i)%radius) then
          error stop 'mesh.templates annulus inner_radius must be smaller than radius.'
        end if
      end select
    end do
  end subroutine validate_mesh_config

end submodule bem_app_config_parser_finalize

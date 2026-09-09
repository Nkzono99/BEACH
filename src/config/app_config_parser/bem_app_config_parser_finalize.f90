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
  if (trim(cfg%sim%field_bc_mode) /= 'free' .and. config_uses_surface_model(cfg, 'conductor')) then
    error stop 'surface_model="conductor" currently requires field_boundary.mode="free".'
  end if
  end procedure validate_field_config

end submodule bem_app_config_parser_finalize

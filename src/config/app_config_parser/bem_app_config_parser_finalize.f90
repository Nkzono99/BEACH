!> 設定の正規化と領域別 preflight の実行順を管理する。
submodule(bem_app_config_parser) bem_app_config_parser_finalize
  implicit none
contains

  module procedure finalize_loaded_config

  call lower_boundary_authoring(cfg, authoring)
  call validate_mesh_part_keys(cfg, authoring)
  call normalize_high_level_config(cfg, authoring)
  call normalize_legacy_physics_config(cfg%sim, cfg%field, cfg%periodic2, cfg%panel)
  call apply_physics_authoring(cfg, authoring)
  cfg%sim%field_solver = lower_ascii(trim(cfg%sim%field_solver))
  cfg%sim%field_normalization = lower_ascii(trim(cfg%sim%field_normalization))
  cfg%sim%field_bc_mode = lower_ascii(trim(cfg%sim%field_bc_mode))
  cfg%sim%field_periodic_far_correction = lower_ascii(trim(cfg%sim%field_periodic_far_correction))

  call validate_simulation_config(cfg)
  call validate_particle_species_config(cfg)
  call resolve_sheath_secondary_role(cfg)
  call validate_surface_current_model_config( &
    cfg, authoring%periodic2%present .and. &
    authoring%periodic2%has_nonzero_mode_backend .and. &
    authoring%periodic2%has_zero_mode_policy .and. &
    authoring%periodic2%has_lower_boundary_model &
    )
  call validate_field_config(cfg)
  end procedure finalize_loaded_config

  !> zero_current の外部シースで SEE 種が 1 つだけなら、sheath.species.secondary_electron を省略してもその種を
  !! SEE の役に当てる。外部根の脱出・帰還は放出電子全体の値なので、登録しない SEE 種は表面収支に入らずセル総電荷を
  !! random walk させる。複数あって指定がなければ、どれを当てるか決められないので設定エラーにする。
  subroutine resolve_sheath_secondary_role(cfg)
    type(app_config), intent(inout) :: cfg
    integer(i32) :: s, found, count

    if (trim(lower_ascii(cfg%surface_current%model)) /= 'zhao_stationary') return
    if (len_trim(cfg%surface_current%secondary_electron_species) > 0) return
    count = 0_i32
    found = 0_i32
    do s = 1_i32, cfg%n_particle_species
      if (.not. cfg%particle_species(s)%enabled) cycle
      if (trim(lower_ascii(cfg%particle_species(s)%source_mode)) /= 'secondary_emission') cycle
      count = count + 1_i32
      found = s
    end do
    if (count == 0_i32) return
    if (count > 1_i32) then
      call stop_config_error('sheath.closure="zero_current" has several secondary_emission species; '// &
                             'set sheath.species.secondary_electron.')
    end if
    cfg%surface_current%secondary_electron_species = cfg%particle_species(found)%species_key
  end subroutine resolve_sheath_secondary_role

  module procedure apply_physics_authoring

  if (authoring%periodic2%present) then
    cfg%periodic2%nonzero_mode_backend = authoring%periodic2%nonzero_mode_backend
    cfg%periodic2%zero_mode_policy = authoring%periodic2%zero_mode_policy
    cfg%periodic2%lower_boundary_model = authoring%periodic2%lower_boundary_model
    cfg%periodic2%upper_fourier_switch = authoring%periodic2%upper_fourier_switch
    cfg%periodic2%upper_fourier_switch_height = authoring%periodic2%upper_fourier_switch_height
    cfg%periodic2%upper_fourier_mode_layers = authoring%periodic2%upper_fourier_mode_layers
    cfg%periodic2%reference_mode_layers = authoring%periodic2%reference_mode_layers
    cfg%periodic2%panel_quadrature_order = authoring%periodic2%panel_quadrature_order
  end if
  ! adaptive 上限は backend を選ばない。正の値と backend の組み合わせは validate_adaptive_batch_config が検査する。
  cfg%periodic2%max_nonzero_mode_potential_step = authoring%periodic2%max_nonzero_mode_potential_step
  end procedure apply_physics_authoring

  module procedure validate_field_config
  integer(i32) :: status
  character(len=256) :: message
  type(field_physics_config) :: field_config
  type(panel_kernel_config) :: panel_config
  call derive_field_panel_config(cfg%sim, field_config, panel_config)
  call validate_active_physics_config(cfg%sim, field_config, cfg%periodic2, panel_config, status, message)
  if (status /= physics_config_ok) call stop_config_error(message)
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

  !> kind と公開キーの組合せは authoring の座標変換前に一度だけ検査する。
  subroutine validate_mesh_part_keys(cfg, authoring)
    type(app_config), intent(in) :: cfg
    type(app_config_authoring), intent(in) :: authoring
    integer :: i
    character(len=:), allocatable :: path

    do i = 1, cfg%n_mesh_parts
      path = mesh_part_config_path(i)
      select case (trim(cfg%mesh_parts(i)%kind))
      case ('obj')
        if (len_trim(authoring%mesh_parts(i)%template_only_key) > 0) then
          error stop path//'.'//trim(authoring%mesh_parts(i)%template_only_key)//' is only valid for template kinds.'
        end if
        if (len_trim(cfg%mesh_parts(i)%path) == 0) error stop path//'.path is required for kind="obj".'
      case ('plane', 'plate_hole', 'plane_hole', 'disk', 'annulus', 'box', 'cylinder', 'sphere')
        if (len_trim(authoring%mesh_parts(i)%obj_only_key) > 0) then
          error stop path//'.'//trim(authoring%mesh_parts(i)%obj_only_key)//' is only valid for kind="obj".'
        end if
      case ('')
        error stop path//'.kind is required.'
      case default
        error stop path//'.kind is unknown.'
      end select
      if (.not. cfg%mesh_parts(i)%enabled) cycle
      select case (trim(cfg%mesh_parts(i)%surface_side_policy))
      case ('normal_plus', 'normal_minus', 'outward_closed')
        continue
      case default
        error stop path//'.surface_side must be normal_plus, normal_minus, or outward_closed.'
      end select
      select case (trim(cfg%mesh_parts(i)%surface_model))
      case ('insulator', 'conductor')
        continue
      case default
        error stop path//'.surface_model must be insulator or conductor.'
      end select
    end do
  end subroutine validate_mesh_part_keys

  !> mesh の構築前に決まる幾何条件を検査する。
  subroutine validate_mesh_config(cfg)
    type(app_config), intent(in) :: cfg
    integer :: i

    do i = 1, cfg%n_mesh_parts
      if (.not. cfg%mesh_parts(i)%enabled) cycle
      select case (trim(lower_ascii(cfg%mesh_parts(i)%kind)))
      case ('plate_hole', 'plane_hole')
        if (cfg%mesh_parts(i)%radius >= 0.5_dp*min(cfg%mesh_parts(i)%size_x, cfg%mesh_parts(i)%size_y)) then
          error stop mesh_part_config_path(i)//'.radius must be smaller than half of min(size_x, size_y).'
        end if
      case ('annulus')
        if (cfg%mesh_parts(i)%inner_radius >= cfg%mesh_parts(i)%radius) then
          error stop mesh_part_config_path(i)//'.inner_radius must be smaller than radius.'
        end if
      end select
    end do
  end subroutine validate_mesh_config

end submodule bem_app_config_parser_finalize

!> TOML設定ファイルを `toml-f` で読み込み、`app_config` へ反映する。
module bem_app_config_parser
  use bem_kinds, only: dp, i32
  use bem_types, only: bc_open, bc_reflect, bc_periodic, bc_redistributed_reflect
  use bem_app_config_types, only: &
    app_config, particle_species_spec, template_spec, max_templates, max_particle_species, particle_bc_inherit, &
    particle_inflow_reservoir, species_from_defaults
  use bem_physics_config_types, only: &
    field_physics_config, panel_kernel_config, normalize_legacy_physics_config, derive_field_panel_config, &
    validate_active_physics_config, physics_config_ok
  use bem_app_config_authoring, only: &
    app_config_authoring, particle_authoring_spec, template_authoring_spec, mesh_group_authoring_spec, &
    domain_authoring_spec, field_boundary_authoring_spec, particle_boundary_authoring_spec, reservoir_authoring_spec, &
    init_app_config_authoring, ensure_authoring_particle_capacity, ensure_authoring_template_capacity, &
    ensure_authoring_group_capacity, normalize_high_level_config, lower_boundary_authoring
  use bem_string_utils, only: lower_ascii
  use bem_config_toml, only: require_toml_success, get_toml_real, get_toml_int, get_toml_logical, &
                             get_toml_string, get_toml_real2, get_toml_real3, get_toml_real4, &
                             get_toml_real_scalar_or_array3, get_toml_particle_boundary_mode, &
                             get_toml_boundary_inflow_mode, stop_config_error
  use bem_injection_flux, only: compute_inflow_flux_from_drifting_maxwellian
  use bem_injection_geometry, only: compute_face_area_from_bounds, resolve_face_axes
  ! Intel 2023 で同名の TOML constructor と host association が衝突しないよう型名を分ける。
  use tomlf, only: config_toml_array => toml_array, config_toml_table => toml_table, &
                   toml_error, toml_key, toml_parse, toml_stat, get_value, toml_len => len
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  private :: finalize_loaded_config
  private :: validate_simulation_config
  private :: validate_mesh_config
  private :: validate_field_config
  private :: validate_particle_species_config
  private :: validate_source_workload
  private :: validate_particle_boundary_override
  private :: validate_particle_boundary_inflow
  private :: is_automatic_current_species
  private :: validate_surface_current_model_config
  private :: stop_config_error

  interface
    !> 読み込み済み設定を正規化し、派生値を確定して全体整合性を検証する。
    module subroutine finalize_loaded_config(cfg, authoring)
      type(app_config), intent(inout) :: cfg
      type(app_config_authoring), intent(in) :: authoring
    end subroutine finalize_loaded_config

    !> `sim.batch_duration` と `sim.batch_duration_step` の整合を検証して確定値を反映する。
    module subroutine resolve_batch_duration(cfg)
      type(app_config), intent(inout) :: cfg
    end subroutine resolve_batch_duration

    !> `sim.e0` または `sim.e0_abs` + angle 指定を内部ベクトルへ正規化する。
    module subroutine resolve_external_e_field(cfg)
      type(app_config), intent(inout) :: cfg
    end subroutine resolve_external_e_field

    !> `reservoir_face` 粒子種の入力値を検証し、必要なら `w_particle` を解決する。
    module subroutine validate_reservoir_species(cfg, species_idx)
      type(app_config), intent(inout) :: cfg
      integer, intent(in) :: species_idx
    end subroutine validate_reservoir_species

    !> box境界に結び付いた reservoir 流入設定を検証し、必要なら重みを解決する。
    module subroutine validate_boundary_inflow_species(cfg, species_idx)
      type(app_config), intent(inout) :: cfg
      integer, intent(in) :: species_idx
    end subroutine validate_boundary_inflow_species

    !> `plane_source` 粒子種の物理量と内部矩形面を検証し、必要なら重みを解決する。
    module subroutine validate_plane_source_species(cfg, species_idx)
      type(app_config), intent(inout) :: cfg
      integer, intent(in) :: species_idx
    end subroutine validate_plane_source_species

    !> `photo_raycast` 粒子種の入力値を検証し、発射方向などを正規化する。
    module subroutine validate_photo_raycast_species(cfg, species_idx)
      type(app_config), intent(inout) :: cfg
      integer, intent(in) :: species_idx
    end subroutine validate_photo_raycast_species

    module subroutine apply_periodic2_toml_table(table, authoring)
      type(config_toml_table), intent(inout) :: table
      type(app_config_authoring), intent(inout) :: authoring
    end subroutine apply_periodic2_toml_table

    module subroutine apply_domain_toml_table(table, domain)
      type(config_toml_table), intent(inout) :: table
      type(domain_authoring_spec), intent(inout) :: domain
    end subroutine apply_domain_toml_table

    module subroutine apply_field_boundary_toml_table(table, field)
      type(config_toml_table), intent(inout) :: table
      type(field_boundary_authoring_spec), intent(inout) :: field
    end subroutine apply_field_boundary_toml_table

    module subroutine apply_particle_boundary_toml_table(table, particles)
      type(config_toml_table), intent(inout) :: table
      type(particle_boundary_authoring_spec), intent(inout) :: particles
    end subroutine apply_particle_boundary_toml_table

    module subroutine apply_reservoir_toml_table(table, reservoir)
      type(config_toml_table), intent(inout) :: table
      type(reservoir_authoring_spec), intent(inout) :: reservoir
    end subroutine apply_reservoir_toml_table

    module subroutine apply_sim_toml_table(cfg, table)
      type(app_config), intent(inout) :: cfg
      type(config_toml_table), intent(inout) :: table
    end subroutine apply_sim_toml_table

    module subroutine apply_output_toml_table(cfg, table)
      type(app_config), intent(inout) :: cfg
      type(config_toml_table), intent(inout) :: table
    end subroutine apply_output_toml_table

    module subroutine ensure_particle_species_capacity(cfg, required_size)
      type(app_config), intent(inout) :: cfg
      integer, intent(in) :: required_size
    end subroutine ensure_particle_species_capacity

    module subroutine apply_particles_toml_table(cfg, table, authoring)
      type(app_config), intent(inout) :: cfg
      type(config_toml_table), intent(inout) :: table
      type(app_config_authoring), intent(inout) :: authoring
    end subroutine apply_particles_toml_table

    module subroutine read_particle_species_array(cfg, table, key, authoring)
      type(app_config), intent(inout) :: cfg
      type(config_toml_table), intent(inout) :: table
      type(toml_key), intent(in) :: key
      type(app_config_authoring), intent(inout) :: authoring
    end subroutine read_particle_species_array

    module subroutine apply_particles_species_toml_table(spec, table, auth)
      type(particle_species_spec), intent(inout) :: spec
      type(config_toml_table), intent(inout) :: table
      type(particle_authoring_spec), intent(inout) :: auth
    end subroutine apply_particles_species_toml_table

    module subroutine apply_species_boundary_toml_table(spec, table)
      type(particle_species_spec), intent(inout) :: spec
      type(config_toml_table), intent(inout) :: table
    end subroutine apply_species_boundary_toml_table

    module subroutine apply_species_boundary_inflow_toml_table(spec, table)
      type(particle_species_spec), intent(inout) :: spec
      type(config_toml_table), intent(inout) :: table
    end subroutine apply_species_boundary_inflow_toml_table

    module subroutine ensure_template_capacity(cfg, required_size)
      type(app_config), intent(inout) :: cfg
      integer, intent(in) :: required_size
    end subroutine ensure_template_capacity

    module subroutine apply_mesh_toml_table(cfg, table, authoring)
      type(app_config), intent(inout) :: cfg
      type(config_toml_table), intent(inout) :: table
      type(app_config_authoring), intent(inout) :: authoring
    end subroutine apply_mesh_toml_table

    module subroutine read_mesh_groups_table(table, key, authoring)
      type(config_toml_table), intent(inout) :: table
      type(toml_key), intent(in) :: key
      type(app_config_authoring), intent(inout) :: authoring
    end subroutine read_mesh_groups_table

    module subroutine apply_mesh_group_toml_table(group, table)
      type(mesh_group_authoring_spec), intent(inout) :: group
      type(config_toml_table), intent(inout) :: table
    end subroutine apply_mesh_group_toml_table

    module subroutine read_template_array(cfg, table, key, authoring)
      type(app_config), intent(inout) :: cfg
      type(config_toml_table), intent(inout) :: table
      type(toml_key), intent(in) :: key
      type(app_config_authoring), intent(inout) :: authoring
    end subroutine read_template_array

    module subroutine apply_template_toml_table(spec, table, auth)
      type(template_spec), intent(inout) :: spec
      type(config_toml_table), intent(inout) :: table
      type(template_authoring_spec), intent(inout) :: auth
    end subroutine apply_template_toml_table

    module subroutine apply_surface_current_model_toml_table(cfg, table)
      type(app_config), intent(inout) :: cfg
      type(config_toml_table), intent(inout) :: table
    end subroutine apply_surface_current_model_toml_table

    module subroutine validate_simulation_config(cfg)
      type(app_config), intent(inout) :: cfg
    end subroutine validate_simulation_config

    module subroutine validate_mesh_config(cfg)
      type(app_config), intent(inout) :: cfg
    end subroutine validate_mesh_config

    module subroutine validate_field_config(cfg)
      type(app_config), intent(inout) :: cfg
    end subroutine validate_field_config

    module subroutine validate_particle_species_config( &
      cfg, per_batch_particles, has_dynamic_source_species, has_enabled_volume_seed &
      )
      type(app_config), intent(inout) :: cfg
      integer(i32), intent(out) :: per_batch_particles
      logical, intent(out) :: has_dynamic_source_species, has_enabled_volume_seed
    end subroutine validate_particle_species_config

    module subroutine validate_source_workload(cfg, per_batch_particles, has_dynamic_source_species, has_enabled_volume_seed)
      type(app_config), intent(in) :: cfg
      integer(i32), intent(in) :: per_batch_particles
      logical, intent(in) :: has_dynamic_source_species, has_enabled_volume_seed
    end subroutine validate_source_workload

    module subroutine validate_particle_boundary_override(action, topology_action, context)
      integer(i32), intent(in) :: action, topology_action
      character(len=*), intent(in) :: context
    end subroutine validate_particle_boundary_override

    module subroutine validate_particle_boundary_inflow(inflow, effective_topology_action, effective_particle_action, context)
      integer(i32), intent(in) :: inflow, effective_topology_action, effective_particle_action
      character(len=*), intent(in) :: context
    end subroutine validate_particle_boundary_inflow

    module logical function is_automatic_current_species(cfg, species_idx) result(selected)
      type(app_config), intent(in) :: cfg
      integer, intent(in) :: species_idx
    end function is_automatic_current_species

    module subroutine validate_surface_current_model_config(cfg, periodic2_split_explicit)
      type(app_config), intent(in) :: cfg
      logical, intent(in) :: periodic2_split_explicit
    end subroutine validate_surface_current_model_config

    module subroutine apply_physics_authoring(cfg, authoring)
      type(app_config), intent(inout) :: cfg
      type(app_config_authoring), intent(in) :: authoring
    end subroutine apply_physics_authoring

  end interface

contains

  !> `.toml` 拡張子の設定ファイルを読み込み、既存値へ上書き適用する。
  !! @param[in] path 読み込む設定ファイルパス（`.toml` 必須）。
  !! @param[inout] cfg 読み込み結果で上書きするアプリ設定。
  subroutine load_app_config(path, cfg)
    character(len=*), intent(in) :: path
    type(app_config), intent(inout) :: cfg

    if (.not. has_suffix(lower_ascii(trim(path)), '.toml')) then
      error stop 'Only TOML config is supported. Please pass a .toml file.'
    end if
    call load_toml_config(path, cfg)
  end subroutine load_app_config

  !> 文字列が指定した接尾辞で終わるかを判定する。
  pure logical function has_suffix(s, suffix)
    character(len=*), intent(in) :: s
    character(len=*), intent(in) :: suffix
    integer :: ls, lf

    ls = len_trim(s)
    lf = len_trim(suffix)
    if (lf > ls) then
      has_suffix = .false.
    else
      has_suffix = (s(ls - lf + 1:ls) == suffix(1:lf))
    end if
  end function has_suffix

  !> TOML 文書を `toml-f` で解釈して設定へ反映する。
  !! 現在は `sim` / `mesh` / `output` / `[[mesh.templates]]` / `[[particles.species]]` を扱う。
  !! @param[in] path 読み込むTOMLファイルパス。
  !! @param[inout] cfg 読み込み結果で更新するアプリ設定。
  subroutine load_toml_config(path, cfg)
    character(len=*), intent(in) :: path
    type(app_config), intent(inout) :: cfg
    integer :: u, ios
    type(config_toml_table), allocatable :: document
    type(toml_error), allocatable :: parse_error
    type(app_config_authoring) :: authoring

    if (.not. allocated(cfg%templates)) then
      allocate (cfg%templates(max_templates))
      cfg%n_templates = 0_i32
    end if
    if (.not. allocated(cfg%particle_species)) then
      allocate (cfg%particle_species(max_particle_species))
      cfg%particle_species = particle_species_spec()
      cfg%n_particle_species = 0_i32
    end if
    call init_app_config_authoring(authoring, size(cfg%templates), size(cfg%particle_species))

    open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Could not open TOML file.'
    call toml_parse(document, u, parse_error)
    close (u)
    if (allocated(parse_error)) then
      error stop 'Failed to parse TOML config: '//parse_error%message
    end if
    if (.not. allocated(document)) error stop 'Failed to parse TOML config.'

    call apply_toml_document(cfg, document, authoring)
    call document%destroy
    call resolve_surface_current_model_path(path, cfg)
    call finalize_loaded_config(cfg, authoring)
  end subroutine load_toml_config

  !> response table の相対パスを設定ファイルの配置ディレクトリ基準へ解決する。
  subroutine resolve_surface_current_model_path(config_path, cfg)
    character(len=*), intent(in) :: config_path
    type(app_config), intent(inout) :: cfg
    character(len=:), allocatable :: resolved_path
    integer :: directory_end

    if (.not. cfg%surface_current%has_response_table_path) return
    if (len_trim(cfg%surface_current%response_table_path) == 0) return
    if (cfg%surface_current%response_table_path(1:1) == '/') return

    directory_end = scan(trim(config_path), '/', back=.true.)
    if (directory_end == 0) return
    resolved_path = config_path(:directory_end)//trim(cfg%surface_current%response_table_path)
    if (len(resolved_path) > len(cfg%surface_current%response_table_path)) then
      error stop 'Resolved surface_current_model.response_table_path is too long.'
    end if
    cfg%surface_current%response_table_path = resolved_path
  end subroutine resolve_surface_current_model_path

  !> `toml-f` のルートテーブルから既知セクションを読み込む。
  subroutine apply_toml_document(cfg, document, authoring)
    type(app_config), intent(inout) :: cfg
    type(config_toml_table), intent(inout) :: document
    type(app_config_authoring), intent(inout) :: authoring
    type(toml_key), allocatable :: keys(:)
    type(config_toml_table), pointer :: section
    integer :: ikey, stat
    character(len=:), allocatable :: key_name

    call document%get_keys(keys)
    do ikey = 1, size(keys)
      key_name = lower_ascii(trim(keys(ikey)%key))
      nullify (section)
      call get_value(document, keys(ikey), section, requested=.false., stat=stat)
      call require_toml_success(stat, '['//trim(keys(ikey)%key)//']')
      select case (trim(key_name))
      case ('sim')
        if (.not. associated(section)) error stop 'TOML section [sim] must be a table.'
        call apply_sim_toml_table(cfg, section)
      case ('domain')
        if (.not. associated(section)) error stop 'TOML section [domain] must be a table.'
        call apply_domain_toml_table(section, authoring%domain)
      case ('field_boundary')
        if (.not. associated(section)) error stop 'TOML section [field_boundary] must be a table.'
        call apply_field_boundary_toml_table(section, authoring%field_boundary)
      case ('particle_boundary')
        if (.not. associated(section)) error stop 'TOML section [particle_boundary] must be a table.'
        call apply_particle_boundary_toml_table(section, authoring%particle_boundary)
      case ('reservoir')
        if (.not. associated(section)) error stop 'TOML section [reservoir] must be a table.'
        call apply_reservoir_toml_table(section, authoring%reservoir)
      case ('surface_current_model')
        if (.not. associated(section)) error stop 'TOML section [surface_current_model] must be a table.'
        call apply_surface_current_model_toml_table(cfg, section)
      case ('particles')
        if (.not. associated(section)) error stop 'TOML section [particles] must be a table.'
        call apply_particles_toml_table(cfg, section, authoring)
      case ('periodic2')
        if (.not. associated(section)) error stop 'TOML section [periodic2] must be a table.'
        call apply_periodic2_toml_table(section, authoring)
      case ('mesh')
        if (.not. associated(section)) error stop 'TOML section [mesh] must be a table.'
        call apply_mesh_toml_table(cfg, section, authoring)
      case ('output')
        if (.not. associated(section)) error stop 'TOML section [output] must be a table.'
        call apply_output_toml_table(cfg, section)
      case default
        error stop 'Unknown TOML section or top-level key: '//trim(keys(ikey)%key)
      end select
    end do
  end subroutine apply_toml_document

  !> 現在のメッシュ入力設定が conductor 表面を生成し得るかを返す。
  logical function config_uses_conductor_surface_model(cfg) result(uses_conductor)
    type(app_config), intent(in) :: cfg
    character(len=16) :: mode
    logical :: has_obj
    integer :: i

    uses_conductor = .false.
    mode = trim(cfg%mesh_mode)
    select case (mode)
    case ('obj')
      uses_conductor = trim(cfg%mesh_surface_model) == 'conductor'
      return
    case ('auto')
      inquire (file=trim(cfg%obj_path), exist=has_obj)
      if (has_obj) then
        uses_conductor = trim(cfg%mesh_surface_model) == 'conductor'
        return
      end if
    end select
    do i = 1, cfg%n_templates
      if (.not. cfg%templates(i)%enabled) cycle
      if (trim(cfg%templates(i)%surface_model) == 'conductor') then
        uses_conductor = .true.
        return
      end if
    end do
  end function config_uses_conductor_surface_model

end module bem_app_config_parser

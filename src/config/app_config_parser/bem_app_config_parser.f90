!> TOML設定ファイルを `toml-f` で読み込み、`app_config` へ反映する。
module bem_app_config_parser
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
  use bem_types, only: bc_open, bc_reflect, bc_periodic
  use bem_app_config_types, only: &
    app_config, particle_species_spec, template_spec, max_templates, max_particle_species, species_from_defaults
  use bem_string_utils, only: lower_ascii
  use tomlf, only: toml_array, toml_error, toml_key, toml_parse, toml_stat, toml_table, get_value, toml_len => len
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  interface
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

    !> `photo_raycast` 粒子種の入力値を検証し、発射方向などを正規化する。
    module subroutine validate_photo_raycast_species(cfg, species_idx)
      type(app_config), intent(inout) :: cfg
      integer, intent(in) :: species_idx
    end subroutine validate_photo_raycast_species

    !> drifting Maxwellian に基づく片側流入束 `[1/m^2/s]` を返す。
    pure module function compute_inflow_flux_from_drifting_maxwellian( &
      number_density_m3, temperature_k, m_particle, drift_velocity, inward_normal &
      ) result(gamma_in)
      real(dp), intent(in) :: number_density_m3
      real(dp), intent(in) :: temperature_k
      real(dp), intent(in) :: m_particle
      real(dp), intent(in) :: drift_velocity(3)
      real(dp), intent(in) :: inward_normal(3)
      real(dp) :: gamma_in
    end function compute_inflow_flux_from_drifting_maxwellian

    !> 標準正規分布の確率密度関数値を返す。
    pure module function standard_normal_pdf(x) result(pdf)
      real(dp), intent(in) :: x
      real(dp) :: pdf
    end function standard_normal_pdf

    !> 標準正規分布の累積分布関数値を返す。
    pure module function standard_normal_cdf(x) result(cdf)
      real(dp), intent(in) :: x
      real(dp) :: cdf
    end function standard_normal_cdf

    !> 注入面上矩形開口の面積 `[m^2]` を計算する。
    pure module function compute_face_area_from_bounds(inject_face, pos_low, pos_high) result(area)
      character(len=*), intent(in) :: inject_face
      real(dp), intent(in) :: pos_low(3), pos_high(3)
      real(dp) :: area
    end function compute_face_area_from_bounds

    !> 注入面名から接線2軸インデックスを返す。
    pure module subroutine resolve_face_axes(inject_face, axis_t1, axis_t2)
      character(len=*), intent(in) :: inject_face
      integer, intent(out) :: axis_t1, axis_t2
    end subroutine resolve_face_axes

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
    integer :: u, ios, i, axis
    integer(i32) :: per_batch_particles
    integer(i32) :: n_periodic_axes
    logical :: has_dynamic_source_species
    type(toml_table), allocatable :: document
    type(toml_error), allocatable :: parse_error

    if (.not. allocated(cfg%templates)) then
      allocate (cfg%templates(max_templates))
      cfg%n_templates = 0_i32
    end if
    if (.not. allocated(cfg%particle_species)) then
      allocate (cfg%particle_species(max_particle_species))
      cfg%particle_species = particle_species_spec()
      cfg%n_particle_species = 0_i32
    end if

    open (newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Could not open TOML file.'
    call toml_parse(document, u, parse_error)
    close (u)
    if (allocated(parse_error)) then
      error stop 'Failed to parse TOML config: '//parse_error%message
    end if
    if (.not. allocated(document)) error stop 'Failed to parse TOML config.'

    call apply_toml_document(cfg, document)
    call document%destroy

    if (cfg%sim%batch_count <= 0_i32) error stop 'sim.batch_count must be > 0.'
    if (.not. ieee_is_finite(cfg%sim%dt) .or. cfg%sim%dt <= 0.0d0) then
      error stop 'sim.dt must be finite and > 0.'
    end if
    if (cfg%sim%max_step <= 0_i32) error stop 'sim.max_step must be > 0.'
    if (.not. ieee_is_finite(cfg%sim%tol_rel) .or. cfg%sim%tol_rel < 0.0d0) then
      error stop 'sim.tol_rel must be finite and >= 0.'
    end if
    if (.not. ieee_is_finite(cfg%sim%q_floor) .or. cfg%sim%q_floor <= 0.0d0) then
      error stop 'sim.q_floor must be finite and > 0.'
    end if
    if (.not. ieee_is_finite(cfg%sim%softening) .or. cfg%sim%softening < 0.0d0) then
      error stop 'sim.softening must be finite and >= 0.'
    end if
    if (.not. all(ieee_is_finite(cfg%sim%b0))) error stop 'sim.b0 must contain finite values.'
    if (cfg%sim%use_box) then
      if (.not. all(ieee_is_finite(cfg%sim%box_min)) .or. .not. all(ieee_is_finite(cfg%sim%box_max))) then
        error stop 'sim.box_min/box_max must contain finite values.'
      end if
      if (any(cfg%sim%box_max <= cfg%sim%box_min)) then
        error stop 'sim.box_max must be greater than sim.box_min on all axes when sim.use_box=true.'
      end if
    end if
    if (cfg%n_particle_species <= 0_i32) error stop 'At least one [[particles.species]] entry is required.'

    cfg%mesh_mode = lower_ascii(trim(cfg%mesh_mode))
    select case (trim(cfg%mesh_mode))
    case ('auto', 'obj', 'template')
      continue
    case default
      error stop 'mesh.mode must be "auto", "obj", or "template".'
    end select
    cfg%mesh_surface_model = lower_ascii(trim(cfg%mesh_surface_model))
    select case (trim(cfg%mesh_surface_model))
    case ('insulator', 'conductor', 'dielectric')
      continue
    case default
      error stop 'mesh.surface_model must be "insulator", "conductor", or "dielectric".'
    end select
    if (.not. ieee_is_finite(cfg%mesh_epsilon_r) .or. cfg%mesh_epsilon_r < 1.0d0) then
      error stop 'mesh.epsilon_r must be finite and >= 1.'
    end if
    do i = 1, cfg%n_templates
      cfg%templates(i)%surface_model = lower_ascii(trim(cfg%templates(i)%surface_model))
      select case (trim(cfg%templates(i)%surface_model))
      case ('insulator', 'conductor', 'dielectric')
        continue
      case default
        error stop 'mesh.templates.surface_model must be "insulator", "conductor", or "dielectric".'
      end select
      if (.not. ieee_is_finite(cfg%templates(i)%epsilon_r) .or. cfg%templates(i)%epsilon_r < 1.0d0) then
        error stop 'mesh.templates.epsilon_r must be finite and >= 1.'
      end if
    end do
    cfg%sim%field_solver = lower_ascii(trim(cfg%sim%field_solver))
    select case (trim(cfg%sim%field_solver))
    case ('direct', 'treecode', 'fmm', 'auto')
      continue
    case default
      error stop 'sim.field_solver must be "direct", "treecode", "fmm", or "auto".'
    end select
    cfg%sim%field_normalization = lower_ascii(trim(cfg%sim%field_normalization))
    select case (trim(cfg%sim%field_normalization))
    case ('si', 'box', 'mesh', 'length')
      continue
    case default
      error stop 'sim.field_normalization must be "si", "box", "mesh", or "length".'
    end select
    if (.not. ieee_is_finite(cfg%sim%field_length_scale) .or. cfg%sim%field_length_scale <= 0.0d0) then
      error stop 'sim.field_length_scale must be finite and > 0.'
    end if
    if (trim(cfg%sim%field_normalization) == 'length' .and. cfg%sim%field_length_scale <= 0.0d0) then
      error stop 'sim.field_normalization="length" requires sim.field_length_scale > 0.'
    end if
    cfg%sim%field_bc_mode = lower_ascii(trim(cfg%sim%field_bc_mode))
    select case (trim(cfg%sim%field_bc_mode))
    case ('free', 'periodic2')
      continue
    case default
      error stop 'sim.field_bc_mode must be "free" or "periodic2".'
    end select
    cfg%sim%field_periodic_far_correction = lower_ascii(trim(cfg%sim%field_periodic_far_correction))
    select case (trim(cfg%sim%field_periodic_far_correction))
    case ('auto')
      continue
    case ('none')
      continue
    case ('m2l_root_oracle')
      continue
    case default
      error stop 'sim.field_periodic_far_correction must be "auto", "none", '// &
        'or "m2l_root_oracle".'
    end select
    if (trim(cfg%sim%field_periodic_far_correction) == 'm2l_root_oracle') then
      if (trim(cfg%sim%field_solver) /= 'fmm' .or. trim(cfg%sim%field_bc_mode) /= 'periodic2') then
        error stop 'sim.field_periodic_far_correction requires sim.field_solver="fmm" and sim.field_bc_mode="periodic2".'
      end if
      if (cfg%sim%field_periodic_ewald_layers < 1_i32) then
        error stop 'sim.field_periodic_ewald_layers must be >= 1 when far correction is enabled.'
      end if
    end if
    if (cfg%sim%field_periodic_image_layers < 0_i32) then
      error stop 'sim.field_periodic_image_layers must be >= 0.'
    end if
    if (.not. ieee_is_finite(cfg%sim%field_periodic_ewald_alpha) .or. cfg%sim%field_periodic_ewald_alpha < 0.0d0) then
      error stop 'sim.field_periodic_ewald_alpha must be finite and >= 0.'
    end if
    if (cfg%sim%field_periodic_ewald_layers < 0_i32) then
      error stop 'sim.field_periodic_ewald_layers must be >= 0.'
    end if
    select case (trim(cfg%sim%field_solver))
    case ('direct', 'treecode', 'auto')
      if (trim(cfg%sim%field_bc_mode) /= 'free') then
        error stop 'sim.field_bc_mode must be "free" when sim.field_solver is "direct", "treecode", or "auto".'
      end if
    case ('fmm')
      if (trim(cfg%sim%field_bc_mode) == 'periodic2') then
        if (.not. cfg%sim%use_box) then
          error stop 'sim.field_bc_mode="periodic2" requires sim.use_box=true.'
        end if
        n_periodic_axes = 0_i32
        do axis = 1, 3
          if ((cfg%sim%bc_low(axis) == bc_periodic) .neqv. (cfg%sim%bc_high(axis) == bc_periodic)) then
            error stop 'periodic2 requires bc_low(axis)=bc_high(axis)=periodic for periodic axes.'
          end if
          if (cfg%sim%bc_low(axis) == bc_periodic) then
            n_periodic_axes = n_periodic_axes + 1_i32
            if (cfg%sim%box_max(axis) <= cfg%sim%box_min(axis)) then
              error stop 'periodic2 requires positive box length on periodic axes.'
            end if
          end if
        end do
        if (n_periodic_axes /= 2_i32) then
          error stop 'sim.field_bc_mode="periodic2" requires exactly two periodic axes.'
        end if
      end if
    end select
    if (trim(cfg%sim%field_bc_mode) /= 'free' .and. config_uses_conductor_surface_model(cfg)) then
      error stop 'surface_model="conductor" currently requires sim.field_bc_mode="free".'
    end if
    if (.not. ieee_is_finite(cfg%sim%tree_theta) .or. cfg%sim%tree_theta <= 0.0d0 .or. cfg%sim%tree_theta > 1.0d0) then
      error stop 'sim.tree_theta must be finite and satisfy 0 < theta <= 1.'
    end if
    if (cfg%sim%tree_leaf_max < 1_i32) then
      error stop 'sim.tree_leaf_max must be >= 1.'
    end if
    if (cfg%sim%tree_min_nelem < 1_i32) then
      error stop 'sim.tree_min_nelem must be >= 1.'
    end if
    call resolve_external_e_field(cfg)
    cfg%sim%reservoir_potential_model = lower_ascii(trim(cfg%sim%reservoir_potential_model))
    select case (trim(cfg%sim%reservoir_potential_model))
    case ('none', 'infinity_barrier')
      continue
    case default
      error stop 'sim.reservoir_potential_model must be "none" or "infinity_barrier".'
    end select
    if (cfg%sim%injection_face_phi_grid_n < 1_i32) then
      error stop 'sim.injection_face_phi_grid_n must be >= 1.'
    end if
    if (cfg%sim%raycast_max_bounce < 1_i32) then
      error stop 'sim.raycast_max_bounce must be >= 1.'
    end if
    if (.not. ieee_is_finite(cfg%sim%phi_infty)) then
      error stop 'sim.phi_infty must be finite.'
    end if
    cfg%sim%sheath_injection_model = lower_ascii(trim(cfg%sim%sheath_injection_model))
    select case (trim(cfg%sim%sheath_injection_model))
    case ('none', 'zhao_auto', 'zhao_a', 'zhao_b', 'zhao_c', 'floating_no_photo')
      continue
    case default
      error stop 'sim.sheath_injection_model must be "none", "zhao_auto", "zhao_a", "zhao_b", "zhao_c", or "floating_no_photo".'
    end select
    cfg%sim%sheath_electron_drift_mode = lower_ascii(trim(cfg%sim%sheath_electron_drift_mode))
    select case (trim(cfg%sim%sheath_electron_drift_mode))
    case ('normal', 'full')
      continue
    case default
      error stop 'sim.sheath_electron_drift_mode must be "normal" or "full".'
    end select
    cfg%sim%sheath_ion_drift_mode = lower_ascii(trim(cfg%sim%sheath_ion_drift_mode))
    select case (trim(cfg%sim%sheath_ion_drift_mode))
    case ('normal', 'full')
      continue
    case default
      error stop 'sim.sheath_ion_drift_mode must be "normal" or "full".'
    end select
    if (.not. ieee_is_finite(cfg%sim%sheath_alpha_deg) .or. cfg%sim%sheath_alpha_deg < 0.0d0 .or. &
        cfg%sim%sheath_alpha_deg > 90.0d0) then
      error stop 'sim.sheath_alpha_deg must be finite and satisfy 0 <= alpha <= 90.'
    end if
    if (index(trim(cfg%sim%sheath_injection_model), 'zhao_') == 1) then
      if (.not. ieee_is_finite(cfg%sim%sheath_photoelectron_ref_density_cm3) .or. &
          cfg%sim%sheath_photoelectron_ref_density_cm3 <= 0.0d0) then
        error stop 'sim.sheath_photoelectron_ref_density_cm3 must be > 0 for Zhao sheath injection.'
      end if
    end if
    if (cfg%sim%has_sheath_reference_coordinate) then
      if (.not. ieee_is_finite(cfg%sim%sheath_reference_coordinate)) then
        error stop 'sim.sheath_reference_coordinate must be finite.'
      end if
    end if
    if (trim(cfg%sim%sheath_injection_model) /= 'none' .and. trim(cfg%sim%reservoir_potential_model) /= 'none') then
      error stop 'sim.sheath_injection_model currently requires sim.reservoir_potential_model="none".'
    end if
    call resolve_batch_duration(cfg)
    per_batch_particles = 0_i32
    has_dynamic_source_species = .false.
    do i = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(i)%enabled) cycle

      cfg%particle_species(i)%source_mode = lower_ascii(trim(cfg%particle_species(i)%source_mode))
      cfg%particle_species(i)%velocity_distribution = lower_ascii(trim(cfg%particle_species(i)%velocity_distribution))
      cfg%particle_species(i)%velocity_grid_pdf_kind = lower_ascii(trim(cfg%particle_species(i)%velocity_grid_pdf_kind))
      cfg%particle_species(i)%velocity_grid_sampling = lower_ascii(trim(cfg%particle_species(i)%velocity_grid_sampling))
      if (.not. all(ieee_is_finite(cfg%particle_species(i)%pos_low)) .or. &
          .not. all(ieee_is_finite(cfg%particle_species(i)%pos_high))) then
        error stop 'particles.species.pos_low/pos_high must contain finite values.'
      end if
      if (.not. all(ieee_is_finite(cfg%particle_species(i)%drift_velocity))) then
        error stop 'particles.species.drift_velocity must contain finite values.'
      end if
      if (.not. ieee_is_finite(cfg%particle_species(i)%q_particle) .or. &
          abs(cfg%particle_species(i)%q_particle) <= 0.0d0) then
        error stop 'particles.species.q_particle must be finite and non-zero.'
      end if
      if (.not. ieee_is_finite(cfg%particle_species(i)%m_particle) .or. cfg%particle_species(i)%m_particle <= 0.0d0) then
        error stop 'particles.species.m_particle must be finite and > 0.'
      end if
      if (.not. ieee_is_finite(cfg%particle_species(i)%w_particle) .or. cfg%particle_species(i)%w_particle <= 0.0d0) then
        error stop 'particles.species.w_particle must be finite and > 0.'
      end if
      if (cfg%particle_species(i)%has_temperature_ev) then
        if (.not. ieee_is_finite(cfg%particle_species(i)%temperature_ev) .or. &
            cfg%particle_species(i)%temperature_ev < 0.0d0) then
          error stop 'particles.species.temperature_ev must be finite and >= 0.'
        end if
      end if
      if (cfg%particle_species(i)%has_temperature_k) then
        if (.not. ieee_is_finite(cfg%particle_species(i)%temperature_k) .or. &
            cfg%particle_species(i)%temperature_k < 0.0d0) then
          error stop 'particles.species.temperature_k must be finite and >= 0.'
        end if
      end if
      select case (trim(cfg%particle_species(i)%velocity_distribution))
      case ('maxwellian', 'grid')
        continue
      case default
        error stop 'particles.species.velocity_distribution must be "maxwellian" or "grid".'
      end select
      select case (trim(cfg%particle_species(i)%velocity_grid_pdf_kind))
      case ('phase_space', 'flux_weighted')
        continue
      case default
        error stop 'particles.species.velocity_grid_pdf_kind must be "phase_space" or "flux_weighted".'
      end select
      select case (trim(cfg%particle_species(i)%velocity_grid_sampling))
      case ('auto', 'rectilinear', 'discrete')
        continue
      case default
        error stop 'particles.species.velocity_grid_sampling must be "auto", "rectilinear", or "discrete".'
      end select
      select case (trim(cfg%particle_species(i)%source_mode))
      case ('volume_seed')
        if (cfg%particle_species(i)%npcls_per_step < 0_i32) then
          error stop 'particles.species.npcls_per_step must be >= 0.'
        end if
        if (trim(cfg%particle_species(i)%velocity_distribution) /= 'maxwellian' .or. &
            len_trim(cfg%particle_species(i)%velocity_grid_path) > 0 .or. &
            trim(cfg%particle_species(i)%velocity_grid_sampling) /= 'auto' .or. &
            cfg%particle_species(i)%has_particle_flux_m2_s .or. cfg%particle_species(i)%has_current_density_a_m2) then
          error stop 'velocity_distribution="grid" and flux keys are only valid for reservoir_face.'
        end if
        if (cfg%particle_species(i)%has_target_macro_particles_per_batch) then
          error stop 'target_macro_particles_per_batch is only valid for reservoir_face.'
        end if
        if (abs(cfg%particle_species(i)%emit_current_density_a_m2) > 0.0d0 .or. &
            cfg%particle_species(i)%rays_per_batch /= 0_i32 .or. cfg%particle_species(i)%has_ray_direction .or. &
            cfg%particle_species(i)%has_deposit_opposite_charge_on_emit) then
          error stop 'photo_raycast keys are only valid for source_mode="photo_raycast".'
        end if
        per_batch_particles = per_batch_particles + cfg%particle_species(i)%npcls_per_step
      case ('reservoir_face')
        has_dynamic_source_species = .true.
        call validate_reservoir_species(cfg, i)
      case ('photo_raycast')
        has_dynamic_source_species = .true.
        call validate_photo_raycast_species(cfg, i)
      case default
        error stop 'Unknown particles.species.source_mode.'
      end select
    end do

    if (per_batch_particles <= 0_i32 .and. .not. has_dynamic_source_species) then
      error stop 'At least one enabled [[particles.species]] entry must have npcls_per_step > 0.'
    end if
    cfg%n_particles = cfg%sim%batch_count*per_batch_particles
  end subroutine load_toml_config

  !> `toml-f` のルートテーブルから既知セクションを読み込む。
  subroutine apply_toml_document(cfg, document)
    type(app_config), intent(inout) :: cfg
    type(toml_table), intent(inout) :: document
    type(toml_key), allocatable :: keys(:)
    type(toml_table), pointer :: section
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
      case ('particles')
        if (.not. associated(section)) error stop 'TOML section [particles] must be a table.'
        call apply_particles_toml_table(cfg, section)
      case ('mesh')
        if (.not. associated(section)) error stop 'TOML section [mesh] must be a table.'
        call apply_mesh_toml_table(cfg, section)
      case ('output')
        if (.not. associated(section)) error stop 'TOML section [output] must be a table.'
        call apply_output_toml_table(cfg, section)
      case default
        error stop 'Unknown TOML section or top-level key: '//trim(keys(ikey)%key)
      end select
    end do
  end subroutine apply_toml_document

  !> `toml-f` の status を BEACH の停止メッセージへ変換する。
  subroutine require_toml_success(stat, context)
    integer, intent(in) :: stat
    character(len=*), intent(in) :: context

    if (stat /= toml_stat%success) then
      error stop 'Invalid TOML value for '//trim(context)//'.'
    end if
  end subroutine require_toml_success

  !> TOML の数値キーを倍精度実数として読み込む。
  subroutine get_toml_real(table, key, value, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    real(dp), intent(out) :: value
    character(len=*), intent(in) :: context
    integer :: stat

    call get_value(table, key, value, stat=stat)
    call require_toml_success(stat, context)
  end subroutine get_toml_real

  !> TOML の整数キーを `integer(i32)` として読み込む。
  subroutine get_toml_int(table, key, value, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    integer(i32), intent(out) :: value
    character(len=*), intent(in) :: context
    integer :: stat, tmp

    call get_value(table, key, tmp, stat=stat)
    call require_toml_success(stat, context)
    value = int(tmp, i32)
  end subroutine get_toml_int

  !> TOML の論理値キーを読み込む。
  subroutine get_toml_logical(table, key, value, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    logical, intent(out) :: value
    character(len=*), intent(in) :: context
    integer :: stat

    call get_value(table, key, value, stat=stat)
    call require_toml_success(stat, context)
  end subroutine get_toml_logical

  !> TOML の文字列キーを固定長文字列へ読み込む。
  subroutine get_toml_string(table, key, value, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    character(len=*), intent(out) :: value
    character(len=*), intent(in) :: context
    character(len=:), allocatable :: tmp
    integer :: stat

    call get_value(table, key, tmp, stat=stat)
    call require_toml_success(stat, context)
    if (.not. allocated(tmp)) error stop 'Invalid TOML value for '//trim(context)//'.'
    value = ''
    value = trim(tmp)
  end subroutine get_toml_string

  !> TOML の3成分数値配列キーを読み込む。
  subroutine get_toml_real3(table, key, value, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    real(dp), intent(out) :: value(3)
    character(len=*), intent(in) :: context
    type(toml_array), pointer :: array
    integer :: i, stat

    nullify (array)
    call get_value(table, key, array, stat=stat)
    call require_toml_success(stat, context)
    if (.not. associated(array)) error stop 'Invalid TOML value for '//trim(context)//'.'
    if (toml_len(array) /= 3) then
      error stop trim(context)//' must be an array of 3 numbers.'
    end if
    do i = 1, 3
      call get_value(array, i, value(i), stat=stat)
      call require_toml_success(stat, context)
    end do
  end subroutine get_toml_real3

  !> TOML の文字列キーを境界条件モードへ変換する。
  subroutine get_toml_boundary_mode(table, key, value, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    integer(i32), intent(out) :: value
    character(len=*), intent(in) :: context
    character(len=64) :: mode

    call get_toml_string(table, key, mode, context)
    select case (trim(lower_ascii(mode)))
    case ('open', 'outflow', 'escape')
      value = bc_open
    case ('reflect', 'reflection')
      value = bc_reflect
    case ('periodic')
      value = bc_periodic
    case default
      error stop 'Unknown boundary condition mode in [sim].'
    end select
  end subroutine get_toml_boundary_mode

  !> `[[mesh.templates]]` の読み込み数に応じてテンプレート配列容量を拡張する。
  !! @param[inout] cfg 容量拡張対象のアプリ設定。
  !! @param[in] required_size 必要最小要素数。
  subroutine ensure_template_capacity(cfg, required_size)
    type(app_config), intent(inout) :: cfg
    integer, intent(in) :: required_size
    type(template_spec), allocatable :: grown(:)
    integer :: old_capacity, new_capacity

    if (required_size <= 0) return
    if (allocated(cfg%templates)) then
      old_capacity = size(cfg%templates)
    else
      old_capacity = 0
    end if
    if (old_capacity >= required_size) return

    new_capacity = max(required_size, max(max_templates, max(1, 2*old_capacity)))
    allocate (grown(new_capacity))
    if (old_capacity > 0) grown(1:old_capacity) = cfg%templates(1:old_capacity)
    call move_alloc(grown, cfg%templates)
  end subroutine ensure_template_capacity

  !> `[[particles.species]]` の読み込み数に応じて粒子種配列容量を拡張する。
  !! @param[inout] cfg 容量拡張対象のアプリ設定。
  !! @param[in] required_size 必要最小要素数。
  subroutine ensure_particle_species_capacity(cfg, required_size)
    type(app_config), intent(inout) :: cfg
    integer, intent(in) :: required_size
    type(particle_species_spec), allocatable :: grown(:)
    integer :: old_capacity, new_capacity

    if (required_size <= 0) return
    if (allocated(cfg%particle_species)) then
      old_capacity = size(cfg%particle_species)
    else
      old_capacity = 0
    end if
    if (old_capacity >= required_size) return

    new_capacity = max(required_size, max(max_particle_species, max(1, 2*old_capacity)))
    allocate (grown(new_capacity))
    grown = particle_species_spec()
    if (old_capacity > 0) grown(1:old_capacity) = cfg%particle_species(1:old_capacity)
    call move_alloc(grown, cfg%particle_species)
  end subroutine ensure_particle_species_capacity

  !> `[sim]` TOML テーブルを `sim_config` へ適用する。
  subroutine apply_sim_toml_table(cfg, table)
    type(app_config), intent(inout) :: cfg
    type(toml_table), intent(inout) :: table
    type(toml_key), allocatable :: keys(:)
    integer :: ikey
    character(len=:), allocatable :: k

    call table%get_keys(keys)
    do ikey = 1, size(keys)
      k = lower_ascii(trim(keys(ikey)%key))
      select case (trim(k))
      case ('dt')
        call get_toml_real(table, keys(ikey), cfg%sim%dt, 'sim.dt')
      case ('rng_seed')
        call get_toml_int(table, keys(ikey), cfg%sim%rng_seed, 'sim.rng_seed')
      case ('batch_count')
        call get_toml_int(table, keys(ikey), cfg%sim%batch_count, 'sim.batch_count')
      case ('batch_duration')
        call get_toml_real(table, keys(ikey), cfg%sim%batch_duration, 'sim.batch_duration')
        cfg%sim%has_batch_duration = .true.
      case ('batch_duration_step')
        call get_toml_real(table, keys(ikey), cfg%sim%batch_duration_step, 'sim.batch_duration_step')
        cfg%sim%has_batch_duration_step = .true.
      case ('max_step')
        call get_toml_int(table, keys(ikey), cfg%sim%max_step, 'sim.max_step')
      case ('tol_rel')
        call get_toml_real(table, keys(ikey), cfg%sim%tol_rel, 'sim.tol_rel')
      case ('q_floor')
        call get_toml_real(table, keys(ikey), cfg%sim%q_floor, 'sim.q_floor')
      case ('softening')
        call get_toml_real(table, keys(ikey), cfg%sim%softening, 'sim.softening')
      case ('field_solver')
        call get_toml_string(table, keys(ikey), cfg%sim%field_solver, 'sim.field_solver')
        cfg%sim%field_solver = lower_ascii(trim(cfg%sim%field_solver))
      case ('field_normalization')
        call get_toml_string(table, keys(ikey), cfg%sim%field_normalization, 'sim.field_normalization')
        cfg%sim%field_normalization = lower_ascii(trim(cfg%sim%field_normalization))
      case ('field_length_scale')
        call get_toml_real(table, keys(ikey), cfg%sim%field_length_scale, 'sim.field_length_scale')
      case ('field_bc_mode')
        call get_toml_string(table, keys(ikey), cfg%sim%field_bc_mode, 'sim.field_bc_mode')
        cfg%sim%field_bc_mode = lower_ascii(trim(cfg%sim%field_bc_mode))
      case ('field_periodic_image_layers')
        call get_toml_int(table, keys(ikey), cfg%sim%field_periodic_image_layers, 'sim.field_periodic_image_layers')
      case ('field_periodic_far_correction')
        call get_toml_string( &
          table, keys(ikey), cfg%sim%field_periodic_far_correction, 'sim.field_periodic_far_correction' &
          )
        cfg%sim%field_periodic_far_correction = lower_ascii(trim(cfg%sim%field_periodic_far_correction))
      case ('field_periodic_ewald_alpha')
        call get_toml_real(table, keys(ikey), cfg%sim%field_periodic_ewald_alpha, 'sim.field_periodic_ewald_alpha')
      case ('field_periodic_ewald_layers')
        call get_toml_int(table, keys(ikey), cfg%sim%field_periodic_ewald_layers, 'sim.field_periodic_ewald_layers')
      case ('tree_theta')
        call get_toml_real(table, keys(ikey), cfg%sim%tree_theta, 'sim.tree_theta')
        cfg%sim%has_tree_theta = .true.
      case ('tree_leaf_max')
        call get_toml_int(table, keys(ikey), cfg%sim%tree_leaf_max, 'sim.tree_leaf_max')
        cfg%sim%has_tree_leaf_max = .true.
      case ('tree_min_nelem')
        call get_toml_int(table, keys(ikey), cfg%sim%tree_min_nelem, 'sim.tree_min_nelem')
      case ('e0')
        call get_toml_real3(table, keys(ikey), cfg%sim%e0, 'sim.e0')
        cfg%sim%has_e0_vector = .true.
      case ('e0_abs')
        call get_toml_real(table, keys(ikey), cfg%sim%e0_abs, 'sim.e0_abs')
        cfg%sim%has_e0_abs = .true.
      case ('e0_phi_xy_deg')
        call get_toml_real(table, keys(ikey), cfg%sim%e0_phi_xy_deg, 'sim.e0_phi_xy_deg')
        cfg%sim%has_e0_phi_xy_deg = .true.
      case ('e0_phi_z_deg')
        call get_toml_real(table, keys(ikey), cfg%sim%e0_phi_z_deg, 'sim.e0_phi_z_deg')
        cfg%sim%has_e0_phi_z_deg = .true.
      case ('b0')
        call get_toml_real3(table, keys(ikey), cfg%sim%b0, 'sim.b0')
      case ('reservoir_potential_model')
        call get_toml_string(table, keys(ikey), cfg%sim%reservoir_potential_model, 'sim.reservoir_potential_model')
        cfg%sim%reservoir_potential_model = lower_ascii(trim(cfg%sim%reservoir_potential_model))
      case ('phi_infty')
        call get_toml_real(table, keys(ikey), cfg%sim%phi_infty, 'sim.phi_infty')
      case ('injection_face_phi_grid_n')
        call get_toml_int(table, keys(ikey), cfg%sim%injection_face_phi_grid_n, 'sim.injection_face_phi_grid_n')
      case ('raycast_max_bounce')
        call get_toml_int(table, keys(ikey), cfg%sim%raycast_max_bounce, 'sim.raycast_max_bounce')
      case ('sheath_injection_model')
        call get_toml_string(table, keys(ikey), cfg%sim%sheath_injection_model, 'sim.sheath_injection_model')
        cfg%sim%sheath_injection_model = lower_ascii(trim(cfg%sim%sheath_injection_model))
      case ('sheath_alpha_deg')
        call get_toml_real(table, keys(ikey), cfg%sim%sheath_alpha_deg, 'sim.sheath_alpha_deg')
      case ('sheath_photoelectron_ref_density_cm3')
        call get_toml_real( &
          table, keys(ikey), cfg%sim%sheath_photoelectron_ref_density_cm3, &
          'sim.sheath_photoelectron_ref_density_cm3' &
          )
      case ('sheath_reference_coordinate')
        call get_toml_real(table, keys(ikey), cfg%sim%sheath_reference_coordinate, 'sim.sheath_reference_coordinate')
        cfg%sim%has_sheath_reference_coordinate = .true.
      case ('sheath_electron_drift_mode')
        call get_toml_string(table, keys(ikey), cfg%sim%sheath_electron_drift_mode, 'sim.sheath_electron_drift_mode')
        cfg%sim%sheath_electron_drift_mode = lower_ascii(trim(cfg%sim%sheath_electron_drift_mode))
      case ('sheath_ion_drift_mode')
        call get_toml_string(table, keys(ikey), cfg%sim%sheath_ion_drift_mode, 'sim.sheath_ion_drift_mode')
        cfg%sim%sheath_ion_drift_mode = lower_ascii(trim(cfg%sim%sheath_ion_drift_mode))
      case ('use_box')
        call get_toml_logical(table, keys(ikey), cfg%sim%use_box, 'sim.use_box')
      case ('box_min')
        call get_toml_real3(table, keys(ikey), cfg%sim%box_min, 'sim.box_min')
      case ('box_max')
        call get_toml_real3(table, keys(ikey), cfg%sim%box_max, 'sim.box_max')
      case ('bc_x_low')
        call get_toml_boundary_mode(table, keys(ikey), cfg%sim%bc_low(1), 'sim.bc_x_low')
      case ('bc_x_high')
        call get_toml_boundary_mode(table, keys(ikey), cfg%sim%bc_high(1), 'sim.bc_x_high')
      case ('bc_y_low')
        call get_toml_boundary_mode(table, keys(ikey), cfg%sim%bc_low(2), 'sim.bc_y_low')
      case ('bc_y_high')
        call get_toml_boundary_mode(table, keys(ikey), cfg%sim%bc_high(2), 'sim.bc_y_high')
      case ('bc_z_low')
        call get_toml_boundary_mode(table, keys(ikey), cfg%sim%bc_low(3), 'sim.bc_z_low')
      case ('bc_z_high')
        call get_toml_boundary_mode(table, keys(ikey), cfg%sim%bc_high(3), 'sim.bc_z_high')
      case default
        error stop 'Unknown key in [sim]: '//trim(keys(ikey)%key)
      end select
    end do
  end subroutine apply_sim_toml_table

  !> `[particles]` TOML テーブルを適用する。
  subroutine apply_particles_toml_table(cfg, table)
    type(app_config), intent(inout) :: cfg
    type(toml_table), intent(inout) :: table
    type(toml_key), allocatable :: keys(:)
    integer :: ikey
    character(len=:), allocatable :: k

    call table%get_keys(keys)
    do ikey = 1, size(keys)
      k = lower_ascii(trim(keys(ikey)%key))
      select case (trim(k))
      case ('species')
        call read_particle_species_array(cfg, table, keys(ikey))
      case default
        error stop 'Unknown key in [particles]: '//trim(keys(ikey)%key)
      end select
    end do
  end subroutine apply_particles_toml_table

  !> `[[particles.species]]` の配列テーブルを読み込む。
  subroutine read_particle_species_array(cfg, table, key)
    type(app_config), intent(inout) :: cfg
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    type(toml_array), pointer :: array
    type(toml_table), pointer :: child
    integer :: ispec, n, stat

    nullify (array)
    call get_value(table, key, array, stat=stat)
    call require_toml_success(stat, 'particles.species')
    if (.not. associated(array)) error stop 'particles.species must be an array of tables.'

    n = toml_len(array)
    call ensure_particle_species_capacity(cfg, n)
    if (n > 0) cfg%particle_species(1:n) = particle_species_spec()
    do ispec = 1, n
      nullify (child)
      call get_value(array, ispec, child, stat=stat)
      call require_toml_success(stat, 'particles.species entry')
      if (.not. associated(child)) error stop 'particles.species entries must be tables.'
      cfg%particle_species(ispec) = species_from_defaults()
      cfg%particle_species(ispec)%enabled = .true.
      call apply_particles_species_toml_table(cfg%particle_species(ispec), child)
    end do
    cfg%n_particle_species = int(n, i32)
  end subroutine read_particle_species_array

  !> `[[particles.species]]` の1要素を粒子種設定へ適用する。
  subroutine apply_particles_species_toml_table(spec, table)
    type(particle_species_spec), intent(inout) :: spec
    type(toml_table), intent(inout) :: table
    type(toml_key), allocatable :: keys(:)
    integer :: ikey
    character(len=:), allocatable :: k

    call table%get_keys(keys)
    do ikey = 1, size(keys)
      k = lower_ascii(trim(keys(ikey)%key))
      select case (trim(k))
      case ('enabled')
        call get_toml_logical(table, keys(ikey), spec%enabled, 'particles.species.enabled')
      case ('npcls_per_step')
        call get_toml_int(table, keys(ikey), spec%npcls_per_step, 'particles.species.npcls_per_step')
        spec%has_npcls_per_step = .true.
      case ('source_mode')
        call get_toml_string(table, keys(ikey), spec%source_mode, 'particles.species.source_mode')
        spec%source_mode = lower_ascii(trim(spec%source_mode))
      case ('number_density_cm3')
        call get_toml_real(table, keys(ikey), spec%number_density_cm3, 'particles.species.number_density_cm3')
        spec%has_number_density_cm3 = .true.
      case ('number_density_m3')
        call get_toml_real(table, keys(ikey), spec%number_density_m3, 'particles.species.number_density_m3')
        spec%has_number_density_m3 = .true.
      case ('q_particle')
        call get_toml_real(table, keys(ikey), spec%q_particle, 'particles.species.q_particle')
      case ('m_particle')
        call get_toml_real(table, keys(ikey), spec%m_particle, 'particles.species.m_particle')
      case ('w_particle')
        call get_toml_real(table, keys(ikey), spec%w_particle, 'particles.species.w_particle')
        spec%has_w_particle = .true.
      case ('target_macro_particles_per_batch')
        call get_toml_int( &
          table, keys(ikey), spec%target_macro_particles_per_batch, &
          'particles.species.target_macro_particles_per_batch' &
          )
        spec%has_target_macro_particles_per_batch = .true.
      case ('pos_low')
        call get_toml_real3(table, keys(ikey), spec%pos_low, 'particles.species.pos_low')
      case ('pos_high')
        call get_toml_real3(table, keys(ikey), spec%pos_high, 'particles.species.pos_high')
      case ('velocity_distribution')
        call get_toml_string(table, keys(ikey), spec%velocity_distribution, 'particles.species.velocity_distribution')
        spec%velocity_distribution = lower_ascii(trim(spec%velocity_distribution))
      case ('velocity_grid_path')
        call get_toml_string(table, keys(ikey), spec%velocity_grid_path, 'particles.species.velocity_grid_path')
      case ('velocity_grid_pdf_kind')
        call get_toml_string(table, keys(ikey), spec%velocity_grid_pdf_kind, 'particles.species.velocity_grid_pdf_kind')
        spec%velocity_grid_pdf_kind = lower_ascii(trim(spec%velocity_grid_pdf_kind))
      case ('velocity_grid_sampling')
        call get_toml_string(table, keys(ikey), spec%velocity_grid_sampling, 'particles.species.velocity_grid_sampling')
        spec%velocity_grid_sampling = lower_ascii(trim(spec%velocity_grid_sampling))
      case ('particle_flux_m2_s')
        call get_toml_real(table, keys(ikey), spec%particle_flux_m2_s, 'particles.species.particle_flux_m2_s')
        spec%has_particle_flux_m2_s = .true.
      case ('current_density_a_m2')
        call get_toml_real(table, keys(ikey), spec%current_density_a_m2, 'particles.species.current_density_a_m2')
        spec%has_current_density_a_m2 = .true.
      case ('drift_velocity')
        call get_toml_real3(table, keys(ikey), spec%drift_velocity, 'particles.species.drift_velocity')
      case ('temperature_k')
        call get_toml_real(table, keys(ikey), spec%temperature_k, 'particles.species.temperature_k')
        spec%has_temperature_k = .true.
      case ('temperature_ev')
        call get_toml_real(table, keys(ikey), spec%temperature_ev, 'particles.species.temperature_ev')
        spec%has_temperature_ev = .true.
      case ('emit_current_density_a_m2')
        call get_toml_real(table, keys(ikey), spec%emit_current_density_a_m2, 'particles.species.emit_current_density_a_m2')
      case ('rays_per_batch')
        call get_toml_int(table, keys(ikey), spec%rays_per_batch, 'particles.species.rays_per_batch')
      case ('deposit_opposite_charge_on_emit')
        call get_toml_logical( &
          table, keys(ikey), spec%deposit_opposite_charge_on_emit, &
          'particles.species.deposit_opposite_charge_on_emit' &
          )
        spec%has_deposit_opposite_charge_on_emit = .true.
      case ('normal_drift_speed')
        call get_toml_real(table, keys(ikey), spec%normal_drift_speed, 'particles.species.normal_drift_speed')
      case ('ray_direction')
        call get_toml_real3(table, keys(ikey), spec%ray_direction, 'particles.species.ray_direction')
        spec%has_ray_direction = .true.
      case ('inject_face')
        call get_toml_string(table, keys(ikey), spec%inject_face, 'particles.species.inject_face')
        spec%inject_face = lower_ascii(trim(spec%inject_face))
      case default
        error stop 'Unknown key in [[particles.species]]: '//trim(keys(ikey)%key)
      end select
    end do
  end subroutine apply_particles_species_toml_table

  !> `[mesh]` TOML テーブルをメッシュ入力設定へ適用する。
  subroutine apply_mesh_toml_table(cfg, table)
    type(app_config), intent(inout) :: cfg
    type(toml_table), intent(inout) :: table
    type(toml_key), allocatable :: keys(:)
    integer :: ikey
    character(len=:), allocatable :: k

    call table%get_keys(keys)
    do ikey = 1, size(keys)
      k = lower_ascii(trim(keys(ikey)%key))
      select case (trim(k))
      case ('mode')
        call get_toml_string(table, keys(ikey), cfg%mesh_mode, 'mesh.mode')
      case ('obj_path')
        call get_toml_string(table, keys(ikey), cfg%obj_path, 'mesh.obj_path')
      case ('surface_model')
        call get_toml_string(table, keys(ikey), cfg%mesh_surface_model, 'mesh.surface_model')
        cfg%mesh_surface_model = lower_ascii(trim(cfg%mesh_surface_model))
      case ('epsilon_r')
        call get_toml_real(table, keys(ikey), cfg%mesh_epsilon_r, 'mesh.epsilon_r')
      case ('obj_scale')
        call get_toml_real(table, keys(ikey), cfg%obj_scale, 'mesh.obj_scale')
      case ('obj_rotation')
        call get_toml_real3(table, keys(ikey), cfg%obj_rotation, 'mesh.obj_rotation')
      case ('obj_offset')
        call get_toml_real3(table, keys(ikey), cfg%obj_offset, 'mesh.obj_offset')
      case ('templates')
        call read_template_array(cfg, table, keys(ikey))
      case default
        error stop 'Unknown key in [mesh]: '//trim(keys(ikey)%key)
      end select
    end do
  end subroutine apply_mesh_toml_table

  !> `[[mesh.templates]]` の配列テーブルを読み込む。
  subroutine read_template_array(cfg, table, key)
    type(app_config), intent(inout) :: cfg
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    type(toml_array), pointer :: array
    type(toml_table), pointer :: child
    integer :: itemplate, n, stat

    nullify (array)
    call get_value(table, key, array, stat=stat)
    call require_toml_success(stat, 'mesh.templates')
    if (.not. associated(array)) error stop 'mesh.templates must be an array of tables.'

    n = toml_len(array)
    call ensure_template_capacity(cfg, n)
    if (n > 0) cfg%templates(1:n) = template_spec()
    do itemplate = 1, n
      nullify (child)
      call get_value(array, itemplate, child, stat=stat)
      call require_toml_success(stat, 'mesh.templates entry')
      if (.not. associated(child)) error stop 'mesh.templates entries must be tables.'
      cfg%templates(itemplate)%enabled = .true.
      call apply_template_toml_table(cfg%templates(itemplate), child)
    end do
    cfg%n_templates = int(n, i32)
  end subroutine read_template_array

  !> `[[mesh.templates]]` の1要素をテンプレート設定へ適用する。
  subroutine apply_template_toml_table(spec, table)
    type(template_spec), intent(inout) :: spec
    type(toml_table), intent(inout) :: table
    type(toml_key), allocatable :: keys(:)
    integer :: ikey
    character(len=:), allocatable :: k

    call table%get_keys(keys)
    do ikey = 1, size(keys)
      k = lower_ascii(trim(keys(ikey)%key))
      select case (trim(k))
      case ('enabled')
        call get_toml_logical(table, keys(ikey), spec%enabled, 'mesh.templates.enabled')
      case ('kind')
        call get_toml_string(table, keys(ikey), spec%kind, 'mesh.templates.kind')
      case ('surface_model')
        call get_toml_string(table, keys(ikey), spec%surface_model, 'mesh.templates.surface_model')
        spec%surface_model = lower_ascii(trim(spec%surface_model))
      case ('epsilon_r')
        call get_toml_real(table, keys(ikey), spec%epsilon_r, 'mesh.templates.epsilon_r')
      case ('center')
        call get_toml_real3(table, keys(ikey), spec%center, 'mesh.templates.center')
      case ('size_x')
        call get_toml_real(table, keys(ikey), spec%size_x, 'mesh.templates.size_x')
      case ('size_y')
        call get_toml_real(table, keys(ikey), spec%size_y, 'mesh.templates.size_y')
      case ('size')
        call get_toml_real3(table, keys(ikey), spec%size, 'mesh.templates.size')
      case ('nx')
        call get_toml_int(table, keys(ikey), spec%nx, 'mesh.templates.nx')
      case ('ny')
        call get_toml_int(table, keys(ikey), spec%ny, 'mesh.templates.ny')
      case ('nz')
        call get_toml_int(table, keys(ikey), spec%nz, 'mesh.templates.nz')
      case ('radius')
        call get_toml_real(table, keys(ikey), spec%radius, 'mesh.templates.radius')
      case ('inner_radius')
        call get_toml_real(table, keys(ikey), spec%inner_radius, 'mesh.templates.inner_radius')
      case ('height')
        call get_toml_real(table, keys(ikey), spec%height, 'mesh.templates.height')
      case ('n_theta')
        call get_toml_int(table, keys(ikey), spec%n_theta, 'mesh.templates.n_theta')
      case ('n_r')
        call get_toml_int(table, keys(ikey), spec%n_r, 'mesh.templates.n_r')
      case ('n_z')
        call get_toml_int(table, keys(ikey), spec%n_z, 'mesh.templates.n_z')
      case ('cap')
        call get_toml_logical(table, keys(ikey), spec%cap, 'mesh.templates.cap')
      case ('cap_top')
        call get_toml_logical(table, keys(ikey), spec%cap_top, 'mesh.templates.cap_top')
        spec%has_cap_top = .true.
      case ('cap_bottom')
        call get_toml_logical(table, keys(ikey), spec%cap_bottom, 'mesh.templates.cap_bottom')
        spec%has_cap_bottom = .true.
      case ('n_lon')
        call get_toml_int(table, keys(ikey), spec%n_lon, 'mesh.templates.n_lon')
      case ('n_lat')
        call get_toml_int(table, keys(ikey), spec%n_lat, 'mesh.templates.n_lat')
      case default
        error stop 'Unknown key in [[mesh.templates]]: '//trim(keys(ikey)%key)
      end select
    end do
  end subroutine apply_template_toml_table

  !> `[output]` TOML テーブルを出力制御設定へ適用する。
  subroutine apply_output_toml_table(cfg, table)
    type(app_config), intent(inout) :: cfg
    type(toml_table), intent(inout) :: table
    type(toml_key), allocatable :: keys(:)
    integer :: ikey
    character(len=:), allocatable :: k

    call table%get_keys(keys)
    do ikey = 1, size(keys)
      k = lower_ascii(trim(keys(ikey)%key))
      select case (trim(k))
      case ('write_files')
        call get_toml_logical(table, keys(ikey), cfg%write_output, 'output.write_files')
      case ('write_mesh_potential')
        call get_toml_logical(table, keys(ikey), cfg%write_mesh_potential, 'output.write_mesh_potential')
      case ('write_potential_history')
        call get_toml_logical(table, keys(ikey), cfg%write_potential_history, 'output.write_potential_history')
      case ('dir')
        call get_toml_string(table, keys(ikey), cfg%output_dir, 'output.dir')
      case ('history_stride')
        call get_toml_int(table, keys(ikey), cfg%history_stride, 'output.history_stride')
      case ('resume')
        call get_toml_logical(table, keys(ikey), cfg%resume_output, 'output.resume')
      case default
        error stop 'Unknown key in [output]: '//trim(keys(ikey)%key)
      end select
    end do
  end subroutine apply_output_toml_table

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

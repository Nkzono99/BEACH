!> TOML設定ファイルを `toml-f` で読み込み、`app_config` へ反映する。
module bem_app_config_parser
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
  use bem_types, only: bc_open, bc_reflect, bc_periodic
  use bem_app_config_types, only: &
    app_config, particle_species_spec, template_spec, max_templates, max_particle_species, species_from_defaults
  use bem_physics_config_types, only: normalize_legacy_physics_config, validate_active_physics_config, physics_config_ok
  use bem_app_config_authoring, only: &
    app_config_authoring, sim_authoring_spec, particle_authoring_spec, template_authoring_spec, mesh_group_authoring_spec, &
    init_app_config_authoring, ensure_authoring_particle_capacity, ensure_authoring_template_capacity, &
    ensure_authoring_group_capacity, normalize_high_level_config
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
    integer :: u, ios, i, j, axis
    integer(i32) :: per_batch_particles
    integer(i32) :: n_periodic_axes
    logical :: has_dynamic_source_species
    character(len=64) :: generated_species_key
    type(toml_table), allocatable :: document
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
    call normalize_high_level_config(cfg, authoring)

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
    if (len_trim(cfg%output_restart_from) > 0 .and. .not. cfg%resume_output) then
      error stop 'output.restart_from requires output.resume = true.'
    end if

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
        if (.not. (trim(cfg%sim%field_solver) == 'direct' .and. authoring%periodic2%present .and. &
                   trim(authoring%periodic2%nonzero_mode_backend) == 'panel_spectral_reference')) then
          error stop 'sim.field_bc_mode must be "free" for this direct/treecode/auto configuration.'
        end if
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
    cfg%sim%open_boundary_model = lower_ascii(trim(cfg%sim%open_boundary_model))
    select case (trim(cfg%sim%open_boundary_model))
    case ('escape', 'potential_barrier')
      continue
    case default
      error stop 'sim.open_boundary_model must be "escape" or "potential_barrier".'
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
      if (len_trim(cfg%particle_species(i)%species_key) == 0) then
        write (generated_species_key, '(a,i0)') 'species_', i
        cfg%particle_species(i)%species_key = trim(generated_species_key)
      end if
      do j = 1, i - 1
        if (trim(cfg%particle_species(i)%species_key) == trim(cfg%particle_species(j)%species_key)) then
          error stop 'particles.species.species_key values must be unique.'
        end if
      end do
      if (.not. cfg%particle_species(i)%enabled) cycle

      cfg%particle_species(i)%source_mode = lower_ascii(trim(cfg%particle_species(i)%source_mode))
      cfg%particle_species(i)%velocity_distribution = lower_ascii(trim(cfg%particle_species(i)%velocity_distribution))
      cfg%particle_species(i)%velocity_grid_pdf_kind = lower_ascii(trim(cfg%particle_species(i)%velocity_grid_pdf_kind))
      cfg%particle_species(i)%velocity_grid_sampling = lower_ascii(trim(cfg%particle_species(i)%velocity_grid_sampling))
      cfg%particle_species(i)%photo_escape_model = lower_ascii(trim(cfg%particle_species(i)%photo_escape_model))
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
            cfg%particle_species(i)%has_deposit_opposite_charge_on_emit .or. &
            cfg%particle_species(i)%has_photo_escape_model) then
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
    call normalize_legacy_physics_config( &
      cfg%sim, cfg%field, cfg%periodic2, cfg%panel, cfg%outer_plasma, cfg%coupling &
      )
    call apply_field_authoring(cfg, authoring)
    call apply_physics_authoring(cfg, authoring)
    if (trim(lower_ascii(cfg%outer_plasma%photoelectron_closure)) == 'individual_return') then
      do i = 1, cfg%n_particle_species
        if (.not. cfg%particle_species(i)%enabled) cycle
        if (trim(lower_ascii(cfg%particle_species(i)%source_mode)) /= 'photo_raycast') cycle
        if (.not. cfg%particle_species(i)%deposit_opposite_charge_on_emit) then
          error stop 'individual photoelectron return requires deposit_opposite_charge_on_emit=true.'
        end if
        if (cfg%particle_species(i)%has_photo_escape_model .and. &
            trim(lower_ascii(cfg%particle_species(i)%photo_escape_model)) /= 'none') then
          error stop 'individual photoelectron return cannot use the legacy photo_escape_model.'
        end if
      end do
    end if
  end subroutine load_toml_config

  !> `toml-f` のルートテーブルから既知セクションを読み込む。
  subroutine apply_toml_document(cfg, document, authoring)
    type(app_config), intent(inout) :: cfg
    type(toml_table), intent(inout) :: document
    type(app_config_authoring), intent(inout) :: authoring
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
        call apply_sim_toml_table(cfg, section, authoring%sim)
      case ('particles')
        if (.not. associated(section)) error stop 'TOML section [particles] must be a table.'
        call apply_particles_toml_table(cfg, section, authoring)
      case ('field')
        if (.not. associated(section)) error stop 'TOML section [field] must be a table.'
        call apply_field_toml_table(section, authoring)
      case ('periodic2')
        if (.not. associated(section)) error stop 'TOML section [periodic2] must be a table.'
        call apply_periodic2_toml_table(section, authoring)
      case ('outer_plasma')
        if (.not. associated(section)) error stop 'TOML section [outer_plasma] must be a table.'
        call apply_outer_plasma_toml_table(section, authoring)
      case ('coupling')
        if (.not. associated(section)) error stop 'TOML section [coupling] must be a table.'
        call apply_coupling_toml_table(section, authoring)
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

  subroutine apply_field_toml_table(table, authoring)
    type(toml_table), intent(inout) :: table
    type(app_config_authoring), intent(inout) :: authoring
    type(toml_key), allocatable :: keys(:)
    integer :: ikey
    character(len=:), allocatable :: k

    call table%get_keys(keys)
    do ikey = 1, size(keys)
      k = lower_ascii(trim(keys(ikey)%key))
      select case (trim(k))
      case ('element_kernel')
        call get_toml_string(table, keys(ikey), authoring%field%element_kernel, 'field.element_kernel')
        authoring%field%element_kernel = lower_ascii(trim(authoring%field%element_kernel))
        authoring%field%has_element_kernel = .true.
      case default
        error stop 'Unknown key in [field]: '//trim(keys(ikey)%key)
      end select
    end do
  end subroutine apply_field_toml_table

  subroutine apply_field_authoring(cfg, authoring)
    type(app_config), intent(inout) :: cfg
    type(app_config_authoring), intent(in) :: authoring
    if (.not. authoring%field%has_element_kernel) return
    select case (trim(authoring%field%element_kernel))
    case ('point')
      cfg%panel%source_model = 'point'
      cfg%panel%kernel_id = 'softened_point'
      cfg%panel%surface_side_policy = 'not_applicable'
    case ('triangle_p0')
      cfg%panel%source_model = 'triangle_p0'
      if (trim(lower_ascii(cfg%sim%field_solver)) == 'fmm') then
        cfg%panel%kernel_id = 'triangle_p0_exact_p2m_near'
      else if (trim(lower_ascii(cfg%sim%field_solver)) == 'auto') then
        cfg%panel%kernel_id = 'triangle_p0_exact_auto'
      else
        cfg%panel%kernel_id = 'triangle_p0_exact_direct'
      end if
      cfg%panel%surface_side_policy = 'per_element'
    case default
      error stop 'field.element_kernel must be "point" or "triangle_p0".'
    end select
  end subroutine apply_field_authoring

  subroutine apply_periodic2_toml_table(table, authoring)
    type(toml_table), intent(inout) :: table
    type(app_config_authoring), intent(inout) :: authoring
    type(toml_key), allocatable :: keys(:)
    integer :: ikey
    character(len=:), allocatable :: k

    authoring%periodic2%present = .true.
    call table%get_keys(keys)
    do ikey = 1, size(keys)
      k = lower_ascii(trim(keys(ikey)%key))
      select case (trim(k))
      case ('nonzero_mode_backend')
        call get_toml_string( &
          table, keys(ikey), authoring%periodic2%nonzero_mode_backend, 'periodic2.nonzero_mode_backend' &
          )
        authoring%periodic2%nonzero_mode_backend = lower_ascii(trim(authoring%periodic2%nonzero_mode_backend))
      case ('zero_mode_policy')
        call get_toml_string(table, keys(ikey), authoring%periodic2%zero_mode_policy, 'periodic2.zero_mode_policy')
        authoring%periodic2%zero_mode_policy = lower_ascii(trim(authoring%periodic2%zero_mode_policy))
      case ('lower_boundary_model')
        call get_toml_string( &
          table, keys(ikey), authoring%periodic2%lower_boundary_model, 'periodic2.lower_boundary_model' &
          )
        authoring%periodic2%lower_boundary_model = lower_ascii(trim(authoring%periodic2%lower_boundary_model))
      case ('reference_mode_layers')
        call get_toml_int( &
          table, keys(ikey), authoring%periodic2%reference_mode_layers, 'periodic2.reference_mode_layers' &
          )
      case ('panel_quadrature_order')
        call get_toml_int( &
          table, keys(ikey), authoring%periodic2%panel_quadrature_order, 'periodic2.panel_quadrature_order' &
          )
      case ('interface_sample_n')
        call get_toml_int(table, keys(ikey), authoring%periodic2%interface_sample_n, 'periodic2.interface_sample_n')
      case ('interface_phi_tolerance')
        call get_toml_real( &
          table, keys(ikey), authoring%periodic2%interface_phi_tolerance, 'periodic2.interface_phi_tolerance' &
          )
      case ('interface_field_tolerance')
        call get_toml_real( &
          table, keys(ikey), authoring%periodic2%interface_field_tolerance, 'periodic2.interface_field_tolerance' &
          )
      case default
        error stop 'Unknown key in [periodic2]: '//trim(keys(ikey)%key)
      end select
    end do
  end subroutine apply_periodic2_toml_table

  subroutine apply_outer_plasma_toml_table(table, authoring)
    type(toml_table), intent(inout) :: table
    type(app_config_authoring), intent(inout) :: authoring
    type(toml_key), allocatable :: keys(:)
    integer :: ikey
    character(len=:), allocatable :: k

    authoring%outer_plasma%present = .true.
    call table%get_keys(keys)
    do ikey = 1, size(keys)
      k = lower_ascii(trim(keys(ikey)%key))
      select case (trim(k))
      case ('model')
        call get_toml_string(table, keys(ikey), authoring%outer_plasma%model, 'outer_plasma.model')
        authoring%outer_plasma%model = lower_ascii(trim(authoring%outer_plasma%model))
      case ('photoelectron_closure')
        call get_toml_string( &
          table, keys(ikey), authoring%outer_plasma%photoelectron_closure, 'outer_plasma.photoelectron_closure' &
          )
        authoring%outer_plasma%photoelectron_closure = &
          lower_ascii(trim(authoring%outer_plasma%photoelectron_closure))
      case ('return_model')
        call get_toml_string(table, keys(ikey), authoring%outer_plasma%return_model, 'outer_plasma.return_model')
        authoring%outer_plasma%return_model = lower_ascii(trim(authoring%outer_plasma%return_model))
      case ('interface_z')
        call get_toml_real(table, keys(ikey), authoring%outer_plasma%interface_z, 'outer_plasma.interface_z')
      case ('infinity_potential')
        call get_toml_real( &
          table, keys(ikey), authoring%outer_plasma%infinity_potential, 'outer_plasma.infinity_potential' &
          )
      case ('debye_length')
        call get_toml_real(table, keys(ikey), authoring%outer_plasma%debye_length, 'outer_plasma.debye_length')
      case ('thermal_voltage')
        call get_toml_real(table, keys(ikey), authoring%outer_plasma%thermal_voltage, 'outer_plasma.thermal_voltage')
      case ('max_linearity_ratio')
        call get_toml_real( &
          table, keys(ikey), authoring%outer_plasma%max_linearity_ratio, 'outer_plasma.max_linearity_ratio' &
          )
      case ('max_gap_ratio')
        call get_toml_real(table, keys(ikey), authoring%outer_plasma%max_gap_ratio, 'outer_plasma.max_gap_ratio')
      case ('max_local_charge_ratio')
        call get_toml_real( &
          table, keys(ikey), authoring%outer_plasma%max_local_charge_ratio, 'outer_plasma.max_local_charge_ratio' &
          )
      case ('photoelectron_histogram_bins')
        call get_toml_int( &
          table, keys(ikey), authoring%outer_plasma%photoelectron_histogram_bins, &
          'outer_plasma.photoelectron_histogram_bins' &
          )
      case ('photoelectron_histogram_energy_max')
        call get_toml_real( &
          table, keys(ikey), authoring%outer_plasma%photoelectron_histogram_energy_max, &
          'outer_plasma.photoelectron_histogram_energy_max' &
          )
      case ('photoelectron_ambient_charge_scale')
        call get_toml_real( &
          table, keys(ikey), authoring%outer_plasma%photoelectron_ambient_charge_scale, &
          'outer_plasma.photoelectron_ambient_charge_scale' &
          )
      case ('max_photoelectron_charge_ratio')
        call get_toml_real( &
          table, keys(ikey), authoring%outer_plasma%max_photoelectron_charge_ratio, &
          'outer_plasma.max_photoelectron_charge_ratio' &
          )
      case default
        error stop 'Unknown key in [outer_plasma]: '//trim(keys(ikey)%key)
      end select
    end do
  end subroutine apply_outer_plasma_toml_table

  subroutine apply_coupling_toml_table(table, authoring)
    type(toml_table), intent(inout) :: table
    type(app_config_authoring), intent(inout) :: authoring
    type(toml_key), allocatable :: keys(:)
    integer :: ikey
    character(len=:), allocatable :: k

    authoring%coupling%present = .true.
    call table%get_keys(keys)
    do ikey = 1, size(keys)
      k = lower_ascii(trim(keys(ikey)%key))
      select case (trim(k))
      case ('update_mode')
        call get_toml_string(table, keys(ikey), authoring%coupling%update_mode, 'coupling.update_mode')
        authoring%coupling%update_mode = lower_ascii(trim(authoring%coupling%update_mode))
      case ('particle_transfer_mode')
        call get_toml_string( &
          table, keys(ikey), authoring%coupling%particle_transfer_mode, 'coupling.particle_transfer_mode' &
          )
        authoring%coupling%particle_transfer_mode = lower_ascii(trim(authoring%coupling%particle_transfer_mode))
      case ('outer_update_stride')
        call get_toml_int( &
          table, keys(ikey), authoring%coupling%outer_update_stride, 'coupling.outer_update_stride' &
          )
      case ('field_evolution_timescale')
        call get_toml_real( &
          table, keys(ikey), authoring%coupling%field_evolution_timescale, 'coupling.field_evolution_timescale' &
          )
      case ('max_frozen_field_ratio')
        call get_toml_real( &
          table, keys(ikey), authoring%coupling%max_frozen_field_ratio, 'coupling.max_frozen_field_ratio' &
          )
      case ('outer_queue_enabled')
        call get_toml_logical( &
          table, keys(ikey), authoring%coupling%outer_queue_enabled, 'coupling.outer_queue_enabled' &
          )
      case default
        error stop 'Unknown key in [coupling]: '//trim(keys(ikey)%key)
      end select
    end do
  end subroutine apply_coupling_toml_table

  subroutine apply_physics_authoring(cfg, authoring)
    type(app_config), intent(inout) :: cfg
    type(app_config_authoring), intent(in) :: authoring
    integer(i32) :: status
    character(len=256) :: message

    if (authoring%periodic2%present) then
      cfg%periodic2%nonzero_mode_backend = authoring%periodic2%nonzero_mode_backend
      cfg%periodic2%zero_mode_policy = authoring%periodic2%zero_mode_policy
      cfg%periodic2%lower_boundary_model = authoring%periodic2%lower_boundary_model
      cfg%periodic2%reference_mode_layers = authoring%periodic2%reference_mode_layers
      cfg%periodic2%panel_quadrature_order = authoring%periodic2%panel_quadrature_order
      cfg%periodic2%interface_sample_n = authoring%periodic2%interface_sample_n
      cfg%periodic2%interface_phi_tolerance = authoring%periodic2%interface_phi_tolerance
      cfg%periodic2%interface_field_tolerance = authoring%periodic2%interface_field_tolerance
    end if
    if (authoring%outer_plasma%present) then
      cfg%outer_plasma%model = authoring%outer_plasma%model
      cfg%outer_plasma%photoelectron_closure = authoring%outer_plasma%photoelectron_closure
      cfg%outer_plasma%return_model = authoring%outer_plasma%return_model
      cfg%outer_plasma%interface_z = authoring%outer_plasma%interface_z
      cfg%outer_plasma%infinity_potential = authoring%outer_plasma%infinity_potential
      cfg%outer_plasma%debye_length = authoring%outer_plasma%debye_length
      cfg%outer_plasma%thermal_voltage = authoring%outer_plasma%thermal_voltage
      cfg%outer_plasma%max_linearity_ratio = authoring%outer_plasma%max_linearity_ratio
      cfg%outer_plasma%max_gap_ratio = authoring%outer_plasma%max_gap_ratio
      cfg%outer_plasma%max_local_charge_ratio = authoring%outer_plasma%max_local_charge_ratio
      cfg%outer_plasma%photoelectron_histogram_bins = authoring%outer_plasma%photoelectron_histogram_bins
      cfg%outer_plasma%photoelectron_histogram_energy_max = &
        authoring%outer_plasma%photoelectron_histogram_energy_max
      cfg%outer_plasma%photoelectron_ambient_charge_scale = authoring%outer_plasma%photoelectron_ambient_charge_scale
      cfg%outer_plasma%max_photoelectron_charge_ratio = authoring%outer_plasma%max_photoelectron_charge_ratio
    end if
    if (authoring%coupling%present) then
      cfg%coupling%update_mode = authoring%coupling%update_mode
      cfg%coupling%particle_transfer_mode = authoring%coupling%particle_transfer_mode
      cfg%coupling%outer_update_stride = authoring%coupling%outer_update_stride
      cfg%coupling%field_evolution_timescale = authoring%coupling%field_evolution_timescale
      cfg%coupling%max_frozen_field_ratio = authoring%coupling%max_frozen_field_ratio
      cfg%coupling%outer_queue_enabled = authoring%coupling%outer_queue_enabled
    end if
    call validate_active_physics_config( &
      cfg%sim, cfg%field, cfg%periodic2, cfg%panel, cfg%outer_plasma, cfg%coupling, status, message &
      )
    if (status /= physics_config_ok) error stop trim(message)
  end subroutine apply_physics_authoring

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

  !> TOML の2成分数値配列キーを読み込む。
  subroutine get_toml_real2(table, key, value, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    real(dp), intent(out) :: value(2)
    character(len=*), intent(in) :: context
    type(toml_array), pointer :: array
    integer :: i, stat

    nullify (array)
    call get_value(table, key, array, stat=stat)
    call require_toml_success(stat, context)
    if (.not. associated(array)) error stop 'Invalid TOML value for '//trim(context)//'.'
    if (toml_len(array) /= 2) then
      error stop trim(context)//' must be an array of 2 numbers.'
    end if
    do i = 1, 2
      call get_value(array, i, value(i), stat=stat)
      call require_toml_success(stat, context)
    end do
  end subroutine get_toml_real2

  !> TOML の scalar または最大3成分数値配列キーを読み込む。
  subroutine get_toml_real_scalar_or_array3(table, key, value, value_len, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    real(dp), intent(out) :: value(3)
    integer(i32), intent(out) :: value_len
    character(len=*), intent(in) :: context
    type(toml_array), pointer :: array
    real(dp) :: scalar_value
    integer :: i, n, stat

    value = 0.0d0
    value_len = 0_i32
    call get_value(table, key, scalar_value, stat=stat)
    if (stat == toml_stat%success) then
      value(1) = scalar_value
      value_len = 1_i32
      return
    end if

    nullify (array)
    call get_value(table, key, array, stat=stat)
    call require_toml_success(stat, context)
    if (.not. associated(array)) error stop 'Invalid TOML value for '//trim(context)//'.'
    n = toml_len(array)
    if (n < 1 .or. n > 3) then
      error stop trim(context)//' must be a number or an array of 1 to 3 numbers.'
    end if
    do i = 1, n
      call get_value(array, i, value(i), stat=stat)
      call require_toml_success(stat, context)
    end do
    value_len = int(n, i32)
  end subroutine get_toml_real_scalar_or_array3

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
  subroutine apply_sim_toml_table(cfg, table, sim_auth)
    type(app_config), intent(inout) :: cfg
    type(toml_table), intent(inout) :: table
    type(sim_authoring_spec), intent(inout) :: sim_auth
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
      case ('open_boundary_model')
        call get_toml_string(table, keys(ikey), cfg%sim%open_boundary_model, 'sim.open_boundary_model')
        cfg%sim%open_boundary_model = lower_ascii(trim(cfg%sim%open_boundary_model))
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
        sim_auth%has_box_min = .true.
      case ('box_max')
        call get_toml_real3(table, keys(ikey), cfg%sim%box_max, 'sim.box_max')
        sim_auth%has_box_max = .true.
      case ('box_origin')
        call get_toml_real3(table, keys(ikey), sim_auth%box_origin, 'sim.box_origin')
        sim_auth%has_box_origin = .true.
      case ('box_size')
        call get_toml_real3(table, keys(ikey), sim_auth%box_size, 'sim.box_size')
        sim_auth%has_box_size = .true.
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
  subroutine apply_particles_toml_table(cfg, table, authoring)
    type(app_config), intent(inout) :: cfg
    type(toml_table), intent(inout) :: table
    type(app_config_authoring), intent(inout) :: authoring
    type(toml_key), allocatable :: keys(:)
    integer :: ikey
    character(len=:), allocatable :: k

    call table%get_keys(keys)
    do ikey = 1, size(keys)
      k = lower_ascii(trim(keys(ikey)%key))
      select case (trim(k))
      case ('species')
        call read_particle_species_array(cfg, table, keys(ikey), authoring)
      case default
        error stop 'Unknown key in [particles]: '//trim(keys(ikey)%key)
      end select
    end do
  end subroutine apply_particles_toml_table

  !> `[[particles.species]]` の配列テーブルを読み込む。
  subroutine read_particle_species_array(cfg, table, key, authoring)
    type(app_config), intent(inout) :: cfg
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    type(app_config_authoring), intent(inout) :: authoring
    type(toml_array), pointer :: array
    type(toml_table), pointer :: child
    integer :: ispec, n, stat

    nullify (array)
    call get_value(table, key, array, stat=stat)
    call require_toml_success(stat, 'particles.species')
    if (.not. associated(array)) error stop 'particles.species must be an array of tables.'

    n = toml_len(array)
    call ensure_particle_species_capacity(cfg, n)
    call ensure_authoring_particle_capacity(authoring, n)
    if (n > 0) cfg%particle_species(1:n) = particle_species_spec()
    if (n > 0) authoring%particle_species(1:n) = particle_authoring_spec()
    do ispec = 1, n
      nullify (child)
      call get_value(array, ispec, child, stat=stat)
      call require_toml_success(stat, 'particles.species entry')
      if (.not. associated(child)) error stop 'particles.species entries must be tables.'
      cfg%particle_species(ispec) = species_from_defaults()
      cfg%particle_species(ispec)%enabled = .true.
      call apply_particles_species_toml_table(cfg%particle_species(ispec), child, authoring%particle_species(ispec))
    end do
    cfg%n_particle_species = int(n, i32)
  end subroutine read_particle_species_array

  !> `[[particles.species]]` の1要素を粒子種設定へ適用する。
  subroutine apply_particles_species_toml_table(spec, table, auth)
    type(particle_species_spec), intent(inout) :: spec
    type(toml_table), intent(inout) :: table
    type(particle_authoring_spec), intent(inout) :: auth
    type(toml_key), allocatable :: keys(:)
    integer :: ikey
    character(len=:), allocatable :: k

    call table%get_keys(keys)
    do ikey = 1, size(keys)
      k = lower_ascii(trim(keys(ikey)%key))
      select case (trim(k))
      case ('species_key')
        call get_toml_string(table, keys(ikey), spec%species_key, 'particles.species.species_key')
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
        auth%has_pos_low = .true.
      case ('pos_high')
        call get_toml_real3(table, keys(ikey), spec%pos_high, 'particles.species.pos_high')
        auth%has_pos_high = .true.
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
      case ('photo_escape_model')
        call get_toml_string(table, keys(ikey), spec%photo_escape_model, 'particles.species.photo_escape_model')
        spec%photo_escape_model = lower_ascii(trim(spec%photo_escape_model))
        spec%has_photo_escape_model = .true.
      case ('normal_drift_speed')
        call get_toml_real(table, keys(ikey), spec%normal_drift_speed, 'particles.species.normal_drift_speed')
      case ('ray_direction')
        call get_toml_real3(table, keys(ikey), spec%ray_direction, 'particles.species.ray_direction')
        spec%has_ray_direction = .true.
      case ('inject_face')
        call get_toml_string(table, keys(ikey), spec%inject_face, 'particles.species.inject_face')
        spec%inject_face = lower_ascii(trim(spec%inject_face))
      case ('inject_region_mode')
        call get_toml_string(table, keys(ikey), auth%inject_region_mode, 'particles.species.inject_region_mode')
        auth%inject_region_mode = lower_ascii(trim(auth%inject_region_mode))
        auth%has_inject_region_mode = .true.
      case ('uv_low')
        call get_toml_real2(table, keys(ikey), auth%uv_low, 'particles.species.uv_low')
        auth%has_uv_low = .true.
      case ('uv_high')
        call get_toml_real2(table, keys(ikey), auth%uv_high, 'particles.species.uv_high')
        auth%has_uv_high = .true.
      case default
        error stop 'Unknown key in [[particles.species]]: '//trim(keys(ikey)%key)
      end select
    end do
  end subroutine apply_particles_species_toml_table

  !> `[mesh]` TOML テーブルをメッシュ入力設定へ適用する。
  subroutine apply_mesh_toml_table(cfg, table, authoring)
    type(app_config), intent(inout) :: cfg
    type(toml_table), intent(inout) :: table
    type(app_config_authoring), intent(inout) :: authoring
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
      case ('surface_side')
        call get_toml_string(table, keys(ikey), cfg%mesh_surface_side_policy, 'mesh.surface_side')
        cfg%mesh_surface_side_policy = lower_ascii(trim(cfg%mesh_surface_side_policy))
      case ('epsilon_r')
        call get_toml_real(table, keys(ikey), cfg%mesh_epsilon_r, 'mesh.epsilon_r')
      case ('obj_scale')
        call get_toml_real(table, keys(ikey), cfg%obj_scale, 'mesh.obj_scale')
      case ('obj_rotation')
        call get_toml_real3(table, keys(ikey), cfg%obj_rotation, 'mesh.obj_rotation')
      case ('obj_offset')
        call get_toml_real3(table, keys(ikey), cfg%obj_offset, 'mesh.obj_offset')
      case ('templates')
        call read_template_array(cfg, table, keys(ikey), authoring)
      case ('groups')
        call read_mesh_groups_table(table, keys(ikey), authoring)
      case default
        error stop 'Unknown key in [mesh]: '//trim(keys(ikey)%key)
      end select
    end do
  end subroutine apply_mesh_toml_table

  !> `[mesh.groups.*]` の table 群を読み込む。
  subroutine read_mesh_groups_table(table, key, authoring)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    type(app_config_authoring), intent(inout) :: authoring
    type(toml_table), pointer :: groups_table, group_table
    type(toml_key), allocatable :: group_keys(:)
    integer :: igroup, stat

    nullify (groups_table)
    call get_value(table, key, groups_table, stat=stat)
    call require_toml_success(stat, 'mesh.groups')
    if (.not. associated(groups_table)) error stop 'mesh.groups must be a table.'

    call groups_table%get_keys(group_keys)
    call ensure_authoring_group_capacity(authoring, size(group_keys))
    authoring%n_groups = int(size(group_keys), i32)
    if (size(group_keys) > 0) authoring%groups(1:size(group_keys)) = mesh_group_authoring_spec()
    do igroup = 1, size(group_keys)
      nullify (group_table)
      call get_value(groups_table, group_keys(igroup), group_table, stat=stat)
      call require_toml_success(stat, 'mesh.groups entry')
      if (.not. associated(group_table)) error stop 'mesh.groups entries must be tables.'
      authoring%groups(igroup)%name = trim(group_keys(igroup)%key)
      call apply_mesh_group_toml_table(authoring%groups(igroup), group_table)
    end do
  end subroutine read_mesh_groups_table

  !> 1つの `[mesh.groups.name]` table を読み込む。
  subroutine apply_mesh_group_toml_table(group, table)
    type(mesh_group_authoring_spec), intent(inout) :: group
    type(toml_table), intent(inout) :: table
    type(toml_key), allocatable :: keys(:)
    integer :: ikey
    character(len=:), allocatable :: k

    call table%get_keys(keys)
    do ikey = 1, size(keys)
      k = lower_ascii(trim(keys(ikey)%key))
      select case (trim(k))
      case ('placement_mode')
        call get_toml_string(table, keys(ikey), group%placement_mode, 'mesh.groups.placement_mode')
        group%placement_mode = lower_ascii(trim(group%placement_mode))
        group%has_placement_mode = .true.
      case ('anchor')
        call get_toml_string(table, keys(ikey), group%anchor, 'mesh.groups.anchor')
        group%anchor = lower_ascii(trim(group%anchor))
        group%has_anchor = .true.
      case ('offset')
        call get_toml_real3(table, keys(ikey), group%offset, 'mesh.groups.offset')
        group%has_offset = .true.
      case ('offset_frac')
        call get_toml_real3(table, keys(ikey), group%offset_frac, 'mesh.groups.offset_frac')
        group%has_offset_frac = .true.
      case ('scale')
        call get_toml_real(table, keys(ikey), group%scale, 'mesh.groups.scale')
        group%has_scale = .true.
      case ('scale_from')
        call get_toml_string(table, keys(ikey), group%scale_from, 'mesh.groups.scale_from')
        group%scale_from = lower_ascii(trim(group%scale_from))
        group%has_scale_from = .true.
      case ('scale_factor')
        call get_toml_real(table, keys(ikey), group%scale_factor, 'mesh.groups.scale_factor')
        group%has_scale_factor = .true.
      case default
        error stop 'Unknown key in [mesh.groups.*]: '//trim(keys(ikey)%key)
      end select
    end do
  end subroutine apply_mesh_group_toml_table

  !> `[[mesh.templates]]` の配列テーブルを読み込む。
  subroutine read_template_array(cfg, table, key, authoring)
    type(app_config), intent(inout) :: cfg
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    type(app_config_authoring), intent(inout) :: authoring
    type(toml_array), pointer :: array
    type(toml_table), pointer :: child
    integer :: itemplate, n, stat

    nullify (array)
    call get_value(table, key, array, stat=stat)
    call require_toml_success(stat, 'mesh.templates')
    if (.not. associated(array)) error stop 'mesh.templates must be an array of tables.'

    n = toml_len(array)
    call ensure_template_capacity(cfg, n)
    call ensure_authoring_template_capacity(authoring, n)
    if (n > 0) cfg%templates(1:n) = template_spec()
    if (n > 0) authoring%templates(1:n) = template_authoring_spec()
    do itemplate = 1, n
      nullify (child)
      call get_value(array, itemplate, child, stat=stat)
      call require_toml_success(stat, 'mesh.templates entry')
      if (.not. associated(child)) error stop 'mesh.templates entries must be tables.'
      cfg%templates(itemplate)%enabled = .true.
      call apply_template_toml_table(cfg%templates(itemplate), child, authoring%templates(itemplate))
    end do
    cfg%n_templates = int(n, i32)
  end subroutine read_template_array

  !> `[[mesh.templates]]` の1要素をテンプレート設定へ適用する。
  subroutine apply_template_toml_table(spec, table, auth)
    type(template_spec), intent(inout) :: spec
    type(toml_table), intent(inout) :: table
    type(template_authoring_spec), intent(inout) :: auth
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
      case ('surface_side')
        call get_toml_string(table, keys(ikey), spec%surface_side_policy, 'mesh.templates.surface_side')
        spec%surface_side_policy = lower_ascii(trim(spec%surface_side_policy))
      case ('epsilon_r')
        call get_toml_real(table, keys(ikey), spec%epsilon_r, 'mesh.templates.epsilon_r')
      case ('center')
        call get_toml_real3(table, keys(ikey), spec%center, 'mesh.templates.center')
        auth%has_center = .true.
      case ('size_x')
        call get_toml_real(table, keys(ikey), spec%size_x, 'mesh.templates.size_x')
        auth%has_size_x = .true.
      case ('size_y')
        call get_toml_real(table, keys(ikey), spec%size_y, 'mesh.templates.size_y')
        auth%has_size_y = .true.
      case ('size')
        call get_toml_real3(table, keys(ikey), spec%size, 'mesh.templates.size')
        auth%has_size = .true.
      case ('nx')
        call get_toml_int(table, keys(ikey), spec%nx, 'mesh.templates.nx')
      case ('ny')
        call get_toml_int(table, keys(ikey), spec%ny, 'mesh.templates.ny')
      case ('nz')
        call get_toml_int(table, keys(ikey), spec%nz, 'mesh.templates.nz')
      case ('radius')
        call get_toml_real(table, keys(ikey), spec%radius, 'mesh.templates.radius')
        auth%has_radius = .true.
      case ('inner_radius')
        call get_toml_real(table, keys(ikey), spec%inner_radius, 'mesh.templates.inner_radius')
        auth%has_inner_radius = .true.
      case ('height')
        call get_toml_real(table, keys(ikey), spec%height, 'mesh.templates.height')
        auth%has_height = .true.
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
      case ('group')
        call get_toml_string(table, keys(ikey), auth%group, 'mesh.templates.group')
        auth%has_group = .true.
      case ('center_local')
        call get_toml_real3(table, keys(ikey), auth%center_local, 'mesh.templates.center_local')
        auth%has_center_local = .true.
      case ('placement_mode')
        call get_toml_string(table, keys(ikey), auth%placement_mode, 'mesh.templates.placement_mode')
        auth%placement_mode = lower_ascii(trim(auth%placement_mode))
        auth%has_placement_mode = .true.
      case ('anchor')
        call get_toml_string(table, keys(ikey), auth%anchor, 'mesh.templates.anchor')
        auth%anchor = lower_ascii(trim(auth%anchor))
        auth%has_anchor = .true.
      case ('offset')
        call get_toml_real3(table, keys(ikey), auth%offset, 'mesh.templates.offset')
        auth%has_offset = .true.
      case ('offset_frac')
        call get_toml_real3(table, keys(ikey), auth%offset_frac, 'mesh.templates.offset_frac')
        auth%has_offset_frac = .true.
      case ('size_mode')
        call get_toml_string(table, keys(ikey), auth%size_mode, 'mesh.templates.size_mode')
        auth%size_mode = lower_ascii(trim(auth%size_mode))
        auth%has_size_mode = .true.
      case ('size_frac')
        call get_toml_real_scalar_or_array3( &
          table, keys(ikey), auth%size_frac, auth%size_frac_len, 'mesh.templates.size_frac' &
          )
        auth%has_size_frac = .true.
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
      case ('restart_from')
        call get_toml_string(table, keys(ikey), cfg%output_restart_from, 'output.restart_from')
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

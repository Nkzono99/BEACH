!> TOML風設定ファイルを `app_config` へ読み込む軽量パーサ。
module bem_app_config_parser
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
  use bem_types, only: bc_open, bc_reflect, bc_periodic
  use bem_app_config_types, only: &
    app_config, particle_species_spec, template_spec, max_templates, max_particle_species, species_from_defaults
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  interface
    !> `sim.batch_duration` と `sim.batch_duration_step` の整合を検証して確定値を反映する。
    module subroutine resolve_batch_duration(cfg)
      type(app_config), intent(inout) :: cfg
    end subroutine resolve_batch_duration

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

    !> 注入面名から法線軸と境界座標値を解決する。
    module subroutine resolve_inject_face(box_min, box_max, inject_face, axis, boundary_value)
      real(dp), intent(in) :: box_min(3), box_max(3)
      character(len=*), intent(in) :: inject_face
      integer, intent(out) :: axis
      real(dp), intent(out) :: boundary_value
    end subroutine resolve_inject_face

    !> 注入面名から内向き法線ベクトルを解決する。
    module subroutine resolve_inward_normal(inject_face, inward_normal)
      character(len=*), intent(in) :: inject_face
      real(dp), intent(out) :: inward_normal(3)
    end subroutine resolve_inward_normal

    !> 粒子種設定から実効数密度 `[1/m^3]` を返す。
    pure module function species_number_density_m3(spec) result(number_density_m3)
      type(particle_species_spec), intent(in) :: spec
      real(dp) :: number_density_m3
    end function species_number_density_m3

    !> 粒子種設定から実効温度 `[K]` を返す。
    pure module function species_temperature_k(spec) result(temperature_k)
      type(particle_species_spec), intent(in) :: spec
      real(dp) :: temperature_k
    end function species_temperature_k

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

    !> `key = value` 行を分割し、正規化キーと値文字列を返す。
    module subroutine split_key_value(line, key, value)
      character(len=*), intent(in) :: line
      character(len=*), intent(out) :: key
      character(len=*), intent(out) :: value
    end subroutine split_key_value

    !> 文字列表現を倍精度実数へ変換する。
    module subroutine parse_real(text, out)
      character(len=*), intent(in) :: text
      real(dp), intent(out) :: out
    end subroutine parse_real

    !> 文字列表現を `integer(i32)` へ変換する。
    module subroutine parse_int(text, out)
      character(len=*), intent(in) :: text
      integer(i32), intent(out) :: out
    end subroutine parse_int

    !> `true/false` 表現を論理値へ変換する。
    module subroutine parse_logical(text, out)
      character(len=*), intent(in) :: text
      logical, intent(out) :: out
    end subroutine parse_logical

    !> 引用符付き/裸の文字列表現を値へ変換する。
    module subroutine parse_string(text, out)
      character(len=*), intent(in) :: text
      character(len=*), intent(out) :: out
    end subroutine parse_string

    !> `[x,y,z]` 形式の3成分ベクトル文字列を配列へ変換する。
    module subroutine parse_real3(text, out)
      character(len=*), intent(in) :: text
      real(dp), intent(out) :: out(3)
    end subroutine parse_real3

    !> 境界条件モード文字列を内部定数へ変換する。
    module subroutine parse_boundary_mode(text, out)
      character(len=*), intent(in) :: text
      integer(i32), intent(out) :: out
    end subroutine parse_boundary_mode

    !> 行文字列から `#` 以降のコメント部分を除去する。
    pure module function strip_comment(line) result(out)
      character(len=*), intent(in) :: line
      character(len=len(line)) :: out
    end function strip_comment

    !> ASCII 英字だけを小文字化した文字列を返す。
    pure module function lower(s) result(o)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: o
    end function lower

    !> 文字列が指定接尾辞で終わるかを判定する。
    pure module function ends_with(s, suffix) result(ends_it)
      character(len=*), intent(in) :: s
      character(len=*), intent(in) :: suffix
      logical :: ends_it
    end function ends_with
  end interface

contains

  !> `.toml` 拡張子の設定ファイルを読み込み、既存値へ上書き適用する。
  !! @param[in] path 読み込む設定ファイルパス（`.toml` 必須）。
  !! @param[inout] cfg 読み込み結果で上書きするアプリ設定。
  subroutine load_app_config(path, cfg)
    character(len=*), intent(in) :: path
    type(app_config), intent(inout) :: cfg

    if (.not. ends_with(lower(trim(path)), '.toml')) then
      error stop 'Only TOML config is supported. Please pass a .toml file.'
    end if
    call load_toml_config(path, cfg)
  end subroutine load_app_config

  !> 最小限の TOML セクションと `key = value` を解釈して設定へ反映する。
  !! 現在は `sim` / `mesh` / `output` / `[[mesh.templates]]` / `[[particles.species]]` を扱う。
  !! @param[in] path 読み込むTOMLファイルパス。
  !! @param[inout] cfg 読み込み結果で更新するアプリ設定。
  subroutine load_toml_config(path, cfg)
    character(len=*), intent(in) :: path
    type(app_config), intent(inout) :: cfg
    integer :: u, ios, i, axis, t_idx, s_idx
    integer(i32) :: per_batch_particles
    integer(i32) :: n_periodic_axes
    logical :: has_dynamic_source_species
    character(len=512) :: raw, line, section

    t_idx = 0
    s_idx = 0
    section = ''
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

    do
      read (u, '(A)', iostat=ios) raw
      if (ios /= 0) exit
      line = strip_comment(trim(raw))
      if (len_trim(line) == 0) cycle

      if (line(1:1) == '[') then
        ! 配列テーブルは専用カウンタを進める。
        if (trim(line) == '[[mesh.templates]]') then
          t_idx = t_idx + 1
          call ensure_template_capacity(cfg, t_idx)
          if (int(t_idx, i32) > cfg%n_templates) cfg%n_templates = int(t_idx, i32)
          cfg%templates(t_idx)%enabled = .true.
          section = 'mesh.template'
        else if (trim(line) == '[[particles.species]]') then
          s_idx = s_idx + 1
          call ensure_particle_species_capacity(cfg, s_idx)
          if (int(s_idx, i32) > cfg%n_particle_species) cfg%n_particle_species = int(s_idx, i32)
          cfg%particle_species(s_idx) = species_from_defaults()
          cfg%particle_species(s_idx)%enabled = .true.
          section = 'particles.species'
        else
          section = lower(trim(adjustl(line(2:len_trim(line) - 1))))
        end if
        cycle
      end if

      select case (trim(section))
      case ('')
        error stop 'Found key-value pair before any TOML section.'
      case ('sim')
        call apply_sim_kv(cfg, line)
      case ('particles')
        call apply_particles_kv(line)
      case ('particles.species')
        call apply_particles_species_kv(cfg%particle_species(s_idx), line)
      case ('mesh')
        call apply_mesh_kv(cfg, line)
      case ('mesh.template')
        call apply_template_kv(cfg%templates(t_idx), line)
      case ('output')
        call apply_output_kv(cfg, line)
      case default
        error stop 'Unknown TOML section: ['//trim(section)//']'
      end select
    end do
    close (u)

    if (cfg%sim%batch_count <= 0_i32) error stop 'sim.batch_count must be > 0.'
    if (s_idx <= 0) error stop 'At least one [[particles.species]] entry is required.'
    if (t_idx > 0) cfg%n_templates = int(t_idx, i32)

    cfg%n_particle_species = int(s_idx, i32)
    cfg%sim%field_solver = lower(trim(cfg%sim%field_solver))
    select case (trim(cfg%sim%field_solver))
    case ('direct', 'treecode', 'fmm', 'auto')
      continue
    case default
      error stop 'sim.field_solver must be "direct", "treecode", "fmm", or "auto".'
    end select
    cfg%sim%field_bc_mode = lower(trim(cfg%sim%field_bc_mode))
    select case (trim(cfg%sim%field_bc_mode))
    case ('free', 'periodic2')
      continue
    case default
      error stop 'sim.field_bc_mode must be "free" or "periodic2".'
    end select
    cfg%sim%field_periodic_far_correction = lower(trim(cfg%sim%field_periodic_far_correction))
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
    if (.not. ieee_is_finite(cfg%sim%tree_theta) .or. cfg%sim%tree_theta <= 0.0d0 .or. cfg%sim%tree_theta > 1.0d0) then
      error stop 'sim.tree_theta must be finite and satisfy 0 < theta <= 1.'
    end if
    if (cfg%sim%tree_leaf_max < 1_i32) then
      error stop 'sim.tree_leaf_max must be >= 1.'
    end if
    if (cfg%sim%tree_min_nelem < 1_i32) then
      error stop 'sim.tree_min_nelem must be >= 1.'
    end if
    cfg%sim%reservoir_potential_model = lower(trim(cfg%sim%reservoir_potential_model))
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
    cfg%sim%sheath_injection_model = lower(trim(cfg%sim%sheath_injection_model))
    select case (trim(cfg%sim%sheath_injection_model))
    case ('none', 'zhao_auto', 'zhao_a', 'zhao_b', 'zhao_c', 'floating_no_photo')
      continue
    case default
      error stop 'sim.sheath_injection_model must be "none", "zhao_auto", "zhao_a", "zhao_b", "zhao_c", or "floating_no_photo".'
    end select
    cfg%sim%sheath_electron_drift_mode = lower(trim(cfg%sim%sheath_electron_drift_mode))
    select case (trim(cfg%sim%sheath_electron_drift_mode))
    case ('normal', 'full')
      continue
    case default
      error stop 'sim.sheath_electron_drift_mode must be "normal" or "full".'
    end select
    cfg%sim%sheath_ion_drift_mode = lower(trim(cfg%sim%sheath_ion_drift_mode))
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
    do i = 1, s_idx
      if (.not. cfg%particle_species(i)%enabled) cycle

      cfg%particle_species(i)%source_mode = lower(trim(cfg%particle_species(i)%source_mode))
      select case (trim(cfg%particle_species(i)%source_mode))
      case ('volume_seed')
        if (cfg%particle_species(i)%npcls_per_step < 0_i32) then
          error stop 'particles.species.npcls_per_step must be >= 0.'
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

  !> `[sim]` セクションのキーを `sim_config` へ適用する。
  !! @param[inout] cfg 更新対象のアプリ設定。
  !! @param[in] line `key = value` 形式の設定行。
  subroutine apply_sim_kv(cfg, line)
    type(app_config), intent(inout) :: cfg
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v

    call split_key_value(line, k, v)
    select case (trim(k))
    case ('dt')
      call parse_real(v, cfg%sim%dt)
    case ('rng_seed')
      call parse_int(v, cfg%sim%rng_seed)
    case ('batch_count')
      call parse_int(v, cfg%sim%batch_count)
    case ('batch_duration')
      call parse_real(v, cfg%sim%batch_duration)
      cfg%sim%has_batch_duration = .true.
    case ('batch_duration_step')
      call parse_real(v, cfg%sim%batch_duration_step)
      cfg%sim%has_batch_duration_step = .true.
    case ('max_step')
      call parse_int(v, cfg%sim%max_step)
    case ('tol_rel')
      call parse_real(v, cfg%sim%tol_rel)
    case ('q_floor')
      call parse_real(v, cfg%sim%q_floor)
    case ('softening')
      call parse_real(v, cfg%sim%softening)
    case ('field_solver')
      call parse_string(v, cfg%sim%field_solver)
      cfg%sim%field_solver = lower(trim(cfg%sim%field_solver))
    case ('field_bc_mode')
      call parse_string(v, cfg%sim%field_bc_mode)
      cfg%sim%field_bc_mode = lower(trim(cfg%sim%field_bc_mode))
    case ('field_periodic_image_layers')
      call parse_int(v, cfg%sim%field_periodic_image_layers)
    case ('field_periodic_far_correction')
      call parse_string(v, cfg%sim%field_periodic_far_correction)
      cfg%sim%field_periodic_far_correction = lower(trim(cfg%sim%field_periodic_far_correction))
    case ('field_periodic_ewald_alpha')
      call parse_real(v, cfg%sim%field_periodic_ewald_alpha)
    case ('field_periodic_ewald_layers')
      call parse_int(v, cfg%sim%field_periodic_ewald_layers)
    case ('tree_theta')
      call parse_real(v, cfg%sim%tree_theta)
      cfg%sim%has_tree_theta = .true.
    case ('tree_leaf_max')
      call parse_int(v, cfg%sim%tree_leaf_max)
      cfg%sim%has_tree_leaf_max = .true.
    case ('tree_min_nelem')
      call parse_int(v, cfg%sim%tree_min_nelem)
    case ('use_hybrid')
      call parse_logical(v, cfg%sim%use_hybrid)
    case ('r_switch_factor')
      call parse_real(v, cfg%sim%r_switch_factor)
    case ('n_sub')
      call parse_int(v, cfg%sim%n_sub)
    case ('softening_factor')
      call parse_real(v, cfg%sim%softening_factor)
    case ('b0')
      call parse_real3(v, cfg%sim%b0)
    case ('reservoir_potential_model')
      call parse_string(v, cfg%sim%reservoir_potential_model)
      cfg%sim%reservoir_potential_model = lower(trim(cfg%sim%reservoir_potential_model))
    case ('phi_infty')
      call parse_real(v, cfg%sim%phi_infty)
    case ('injection_face_phi_grid_n')
      call parse_int(v, cfg%sim%injection_face_phi_grid_n)
    case ('raycast_max_bounce')
      call parse_int(v, cfg%sim%raycast_max_bounce)
    case ('sheath_injection_model')
      call parse_string(v, cfg%sim%sheath_injection_model)
      cfg%sim%sheath_injection_model = lower(trim(cfg%sim%sheath_injection_model))
    case ('sheath_alpha_deg')
      call parse_real(v, cfg%sim%sheath_alpha_deg)
    case ('sheath_photoelectron_ref_density_cm3')
      call parse_real(v, cfg%sim%sheath_photoelectron_ref_density_cm3)
    case ('sheath_reference_coordinate')
      call parse_real(v, cfg%sim%sheath_reference_coordinate)
      cfg%sim%has_sheath_reference_coordinate = .true.
    case ('sheath_electron_drift_mode')
      call parse_string(v, cfg%sim%sheath_electron_drift_mode)
      cfg%sim%sheath_electron_drift_mode = lower(trim(cfg%sim%sheath_electron_drift_mode))
    case ('sheath_ion_drift_mode')
      call parse_string(v, cfg%sim%sheath_ion_drift_mode)
      cfg%sim%sheath_ion_drift_mode = lower(trim(cfg%sim%sheath_ion_drift_mode))
    case ('use_box')
      call parse_logical(v, cfg%sim%use_box)
    case ('box_min')
      call parse_real3(v, cfg%sim%box_min)
    case ('box_max')
      call parse_real3(v, cfg%sim%box_max)
    case ('bc_x_low')
      call parse_boundary_mode(v, cfg%sim%bc_low(1))
    case ('bc_x_high')
      call parse_boundary_mode(v, cfg%sim%bc_high(1))
    case ('bc_y_low')
      call parse_boundary_mode(v, cfg%sim%bc_low(2))
    case ('bc_y_high')
      call parse_boundary_mode(v, cfg%sim%bc_high(2))
    case ('bc_z_low')
      call parse_boundary_mode(v, cfg%sim%bc_low(3))
    case ('bc_z_high')
      call parse_boundary_mode(v, cfg%sim%bc_high(3))
    case default
      error stop 'Unknown key in [sim]: '//trim(k)
    end select
  end subroutine apply_sim_kv

  !> `[particles]` セクションのキーを検証する。
  !! @param[in] line `key = value` 形式の設定行。
  subroutine apply_particles_kv(line)
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v

    call split_key_value(line, k, v)
    error stop 'Unknown key in [particles]: '//trim(k)
  end subroutine apply_particles_kv

  !> `[[particles.species]]` のキーを粒子種設定へ適用する。
  !! @param[inout] spec 更新対象の粒子種設定。
  !! @param[in] line `key = value` 形式の設定行。
  subroutine apply_particles_species_kv(spec, line)
    type(particle_species_spec), intent(inout) :: spec
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v

    call split_key_value(line, k, v)
    select case (trim(k))
    case ('enabled')
      call parse_logical(v, spec%enabled)
    case ('npcls_per_step')
      call parse_int(v, spec%npcls_per_step)
      spec%has_npcls_per_step = .true.
    case ('source_mode')
      call parse_string(v, spec%source_mode)
      spec%source_mode = lower(trim(spec%source_mode))
    case ('number_density_cm3')
      call parse_real(v, spec%number_density_cm3)
      spec%has_number_density_cm3 = .true.
    case ('number_density_m3')
      call parse_real(v, spec%number_density_m3)
      spec%has_number_density_m3 = .true.
    case ('q_particle')
      call parse_real(v, spec%q_particle)
    case ('m_particle')
      call parse_real(v, spec%m_particle)
    case ('w_particle')
      call parse_real(v, spec%w_particle)
      spec%has_w_particle = .true.
    case ('target_macro_particles_per_batch')
      call parse_int(v, spec%target_macro_particles_per_batch)
      spec%has_target_macro_particles_per_batch = .true.
    case ('pos_low')
      call parse_real3(v, spec%pos_low)
    case ('pos_high')
      call parse_real3(v, spec%pos_high)
    case ('drift_velocity')
      call parse_real3(v, spec%drift_velocity)
    case ('temperature_k')
      call parse_real(v, spec%temperature_k)
      spec%has_temperature_k = .true.
    case ('temperature_ev')
      call parse_real(v, spec%temperature_ev)
      spec%has_temperature_ev = .true.
    case ('emit_current_density_a_m2')
      call parse_real(v, spec%emit_current_density_a_m2)
    case ('rays_per_batch')
      call parse_int(v, spec%rays_per_batch)
    case ('deposit_opposite_charge_on_emit')
      call parse_logical(v, spec%deposit_opposite_charge_on_emit)
      spec%has_deposit_opposite_charge_on_emit = .true.
    case ('normal_drift_speed')
      call parse_real(v, spec%normal_drift_speed)
    case ('ray_direction')
      call parse_real3(v, spec%ray_direction)
      spec%has_ray_direction = .true.
    case ('inject_face')
      call parse_string(v, spec%inject_face)
      spec%inject_face = lower(trim(spec%inject_face))
    case default
      error stop 'Unknown key in [[particles.species]]: '//trim(k)
    end select
  end subroutine apply_particles_species_kv

  !> `[mesh]` セクションのキーをメッシュ入力設定へ適用する。
  !! @param[inout] cfg 更新対象のアプリ設定。
  !! @param[in] line `key = value` 形式の設定行。
  subroutine apply_mesh_kv(cfg, line)
    type(app_config), intent(inout) :: cfg
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v

    call split_key_value(line, k, v)
    select case (trim(k))
    case ('mode')
      call parse_string(v, cfg%mesh_mode)
    case ('obj_path')
      call parse_string(v, cfg%obj_path)
    case default
      error stop 'Unknown key in [mesh]: '//trim(k)
    end select
  end subroutine apply_mesh_kv

  !> `[[mesh.templates]]` のキーをテンプレート設定へ適用する。
  !! @param[inout] spec 更新対象のテンプレート設定。
  !! @param[in] line `key = value` 形式の設定行。
  subroutine apply_template_kv(spec, line)
    type(template_spec), intent(inout) :: spec
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v

    call split_key_value(line, k, v)
    select case (trim(k))
    case ('enabled')
      call parse_logical(v, spec%enabled)
    case ('kind')
      call parse_string(v, spec%kind)
    case ('center')
      call parse_real3(v, spec%center)
    case ('size_x')
      call parse_real(v, spec%size_x)
    case ('size_y')
      call parse_real(v, spec%size_y)
    case ('size')
      call parse_real3(v, spec%size)
    case ('nx')
      call parse_int(v, spec%nx)
    case ('ny')
      call parse_int(v, spec%ny)
    case ('nz')
      call parse_int(v, spec%nz)
    case ('radius')
      call parse_real(v, spec%radius)
    case ('inner_radius')
      call parse_real(v, spec%inner_radius)
    case ('height')
      call parse_real(v, spec%height)
    case ('n_theta')
      call parse_int(v, spec%n_theta)
    case ('n_r')
      call parse_int(v, spec%n_r)
    case ('n_z')
      call parse_int(v, spec%n_z)
    case ('cap')
      call parse_logical(v, spec%cap)
    case ('cap_top')
      call parse_logical(v, spec%cap_top)
      spec%has_cap_top = .true.
    case ('cap_bottom')
      call parse_logical(v, spec%cap_bottom)
      spec%has_cap_bottom = .true.
    case ('n_lon')
      call parse_int(v, spec%n_lon)
    case ('n_lat')
      call parse_int(v, spec%n_lat)
    case default
      error stop 'Unknown key in [[mesh.templates]]: '//trim(k)
    end select
  end subroutine apply_template_kv

  !> `[output]` セクションのキーを出力制御設定へ適用する。
  !! @param[inout] cfg 更新対象のアプリ設定。
  !! @param[in] line `key = value` 形式の設定行。
  subroutine apply_output_kv(cfg, line)
    type(app_config), intent(inout) :: cfg
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v

    call split_key_value(line, k, v)
    select case (trim(k))
    case ('write_files')
      call parse_logical(v, cfg%write_output)
    case ('write_mesh_potential')
      call parse_logical(v, cfg%write_mesh_potential)
    case ('dir')
      call parse_string(v, cfg%output_dir)
    case ('history_stride')
      call parse_int(v, cfg%history_stride)
    case ('resume')
      call parse_logical(v, cfg%resume_output)
    case default
      error stop 'Unknown key in [output]: '//trim(k)
    end select
  end subroutine apply_output_kv

end module bem_app_config_parser

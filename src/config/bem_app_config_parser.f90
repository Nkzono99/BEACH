!> TOML風設定ファイルを `app_config` へ読み込む軽量パーサ。
module bem_app_config_parser
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
  use bem_types, only: bc_open, bc_reflect, bc_periodic
  use bem_app_config_types, only: &
    app_config, particle_species_spec, template_spec, max_templates, max_particle_species, species_from_defaults
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

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
    integer :: u, ios, i, t_idx, s_idx
    integer(i32) :: per_batch_particles
    logical :: has_dynamic_source_species
    character(len=512) :: raw, line, section

    t_idx = 0
    s_idx = 0
    section = ''

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
          if (t_idx > max_templates) error stop 'Too many mesh.templates entries.'
          cfg%templates(t_idx)%enabled = .true.
          section = 'mesh.template'
        else if (trim(line) == '[[particles.species]]') then
          s_idx = s_idx + 1
          if (s_idx > max_particle_species) error stop 'Too many particles.species entries.'
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
        error stop ('Unknown TOML section: [' // trim(section) // ']')
      end select
    end do
    close (u)

    if (cfg%sim%batch_count <= 0_i32) error stop 'sim.batch_count must be > 0.'
    if (s_idx <= 0) error stop 'At least one [[particles.species]] entry is required.'

    cfg%n_particle_species = s_idx
    cfg%sim%field_solver = lower(trim(cfg%sim%field_solver))
    select case (trim(cfg%sim%field_solver))
    case ('direct', 'treecode', 'auto')
      continue
    case default
      error stop 'sim.field_solver must be "direct", "treecode", or "auto".'
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
    cfg%n_particles = cfg%sim%batch_count * per_batch_particles
  end subroutine load_toml_config

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
    case ('tree_theta')
      call parse_real(v, cfg%sim%tree_theta)
    case ('tree_leaf_max')
      call parse_int(v, cfg%sim%tree_leaf_max)
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
      error stop ('Unknown key in [sim]: ' // trim(k))
    end select
  end subroutine apply_sim_kv

  !> `[particles]` セクションのキーを検証する。
  !! @param[in] line `key = value` 形式の設定行。
  subroutine apply_particles_kv(line)
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v

    call split_key_value(line, k, v)
    error stop ('Unknown key in [particles]: ' // trim(k))
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
      error stop ('Unknown key in [[particles.species]]: ' // trim(k))
    end select
  end subroutine apply_particles_species_kv

  !> `sim.batch_duration` / `sim.batch_duration_step` を解決して確定値へ反映する。
  !! @param[inout] cfg 検証・更新対象のアプリ設定。
  subroutine resolve_batch_duration(cfg)
    type(app_config), intent(inout) :: cfg
    real(dp) :: batch_duration

    if (cfg%sim%has_batch_duration .and. cfg%sim%has_batch_duration_step) then
      error stop 'sim.batch_duration and sim.batch_duration_step cannot be used together.'
    end if

    if (cfg%sim%has_batch_duration_step) then
      if (.not. ieee_is_finite(cfg%sim%batch_duration_step) .or. cfg%sim%batch_duration_step <= 0.0d0) then
        error stop 'sim.batch_duration_step must be > 0.'
      end if
      if (.not. ieee_is_finite(cfg%sim%dt) .or. cfg%sim%dt <= 0.0d0) then
        error stop 'sim.dt must be > 0 when sim.batch_duration_step is set.'
      end if
      batch_duration = cfg%sim%dt * cfg%sim%batch_duration_step
      if (.not. ieee_is_finite(batch_duration) .or. batch_duration <= 0.0d0) then
        error stop 'sim.batch_duration_step produced invalid sim.batch_duration.'
      end if
      cfg%sim%batch_duration = batch_duration
      cfg%sim%has_batch_duration = .true.
    else if (cfg%sim%has_batch_duration) then
      if (.not. ieee_is_finite(cfg%sim%batch_duration)) then
        error stop 'sim.batch_duration must be finite.'
      end if
    end if
  end subroutine resolve_batch_duration

  !> `reservoir_face` 用の必須項目と整合性を検証する。
  !! @param[inout] cfg 検証対象のアプリ設定。
  !! @param[in] species_idx 検証する粒子種のインデックス（1始まり）。
  subroutine validate_reservoir_species(cfg, species_idx)
    type(app_config), intent(inout) :: cfg
    integer, intent(in) :: species_idx

    integer :: axis, axis_t1, axis_t2
    real(dp) :: boundary_value, area
    real(dp) :: number_density_m3, temperature_k, gamma_in, w_particle
    real(dp) :: inward_normal(3)
    type(particle_species_spec) :: spec

    spec = cfg%particle_species(species_idx)

    if (spec%has_npcls_per_step) then
      error stop 'particles.species.npcls_per_step is auto-computed for reservoir_face.'
    end if
    if (abs(spec%emit_current_density_a_m2) > 0.0d0 .or. spec%rays_per_batch /= 0_i32 .or. &
        spec%has_ray_direction .or. spec%has_deposit_opposite_charge_on_emit) then
      error stop 'photo_raycast keys are not allowed for reservoir_face.'
    end if
    if (spec%has_w_particle .and. spec%has_target_macro_particles_per_batch) then
      error stop 'reservoir_face does not allow both w_particle and target_macro_particles_per_batch.'
    end if
    if (.not. spec%has_w_particle .and. .not. spec%has_target_macro_particles_per_batch) then
      error stop 'reservoir_face requires either w_particle or target_macro_particles_per_batch.'
    end if
    if (spec%has_w_particle) then
      if (spec%w_particle <= 0.0d0) error stop 'particles.species.w_particle must be > 0 for reservoir_face.'
    end if
    if (spec%has_target_macro_particles_per_batch) then
      if (spec%target_macro_particles_per_batch == 0_i32 .or. spec%target_macro_particles_per_batch < -1_i32) then
        error stop 'particles.species.target_macro_particles_per_batch must be > 0 or -1.'
      end if
      if (spec%target_macro_particles_per_batch == -1_i32) then
        if (species_idx == 1) then
          error stop 'particles.species[1].target_macro_particles_per_batch cannot be -1.'
        end if
        if (.not. cfg%particle_species(1)%enabled) then
          error stop 'target_macro_particles_per_batch=-1 requires particles.species[1] to be enabled.'
        end if
        if (trim(lower(cfg%particle_species(1)%source_mode)) /= 'reservoir_face') then
          error stop 'target_macro_particles_per_batch=-1 requires particles.species[1].source_mode="reservoir_face".'
        end if
        if (.not. cfg%particle_species(1)%has_w_particle .or. cfg%particle_species(1)%w_particle <= 0.0d0) then
          error stop 'target_macro_particles_per_batch=-1 requires species[1] to resolve a positive w_particle.'
        end if
      end if
    end if
    if (.not. cfg%sim%use_box) then
      error stop 'particles.species.source_mode="reservoir_face" requires sim.use_box = true.'
    end if
    if (cfg%sim%batch_duration <= 0.0d0) then
      error stop 'sim.batch_duration must be > 0 for reservoir_face.'
    end if
    if (spec%has_number_density_cm3 .and. spec%has_number_density_m3) then
      error stop 'Specify either number_density_cm3 or number_density_m3, not both.'
    end if
    if (.not. spec%has_number_density_cm3 .and. .not. spec%has_number_density_m3) then
      error stop 'reservoir_face requires number_density_cm3 or number_density_m3.'
    end if
    if (spec%has_number_density_cm3) then
      if (spec%number_density_cm3 <= 0.0d0) error stop 'number_density_cm3 must be > 0.'
    else
      if (spec%number_density_m3 <= 0.0d0) error stop 'number_density_m3 must be > 0.'
    end if
    if (spec%has_temperature_ev .and. spec%has_temperature_k) then
      error stop 'Specify either temperature_ev or temperature_k, not both.'
    end if
    if (spec%has_temperature_ev) then
      if (spec%temperature_ev < 0.0d0) error stop 'temperature_ev must be >= 0.'
    else if (spec%has_temperature_k) then
      if (spec%temperature_k < 0.0d0) error stop 'temperature_k must be >= 0.'
    else
      if (spec%temperature_k < 0.0d0) error stop 'temperature_k must be >= 0.'
    end if
    if (spec%m_particle <= 0.0d0) then
      error stop 'm_particle must be > 0.'
    end if

    call resolve_inject_face(cfg%sim%box_min, cfg%sim%box_max, spec%inject_face, axis, boundary_value)
    axis_t1 = modulo(axis, 3) + 1
    axis_t2 = modulo(axis + 1, 3) + 1
    if (abs(spec%pos_low(axis) - boundary_value) > 1.0d-12 .or. abs(spec%pos_high(axis) - boundary_value) > 1.0d-12) then
      error stop 'reservoir_face pos_low/pos_high must lie on the selected box face.'
    end if
    if (spec%pos_high(axis_t1) < spec%pos_low(axis_t1) .or. spec%pos_high(axis_t2) < spec%pos_low(axis_t2)) then
      error stop 'reservoir_face tangential bounds must satisfy pos_high >= pos_low.'
    end if
    if ((spec%pos_high(axis_t1) - spec%pos_low(axis_t1)) <= 0.0d0 .or. &
        (spec%pos_high(axis_t2) - spec%pos_low(axis_t2)) <= 0.0d0) then
      error stop 'reservoir_face opening area must be positive.'
    end if
    area = compute_face_area_from_bounds(spec%inject_face, spec%pos_low, spec%pos_high)
    if (area <= 0.0d0) then
      error stop 'reservoir_face opening area must be positive.'
    end if

    if (spec%has_target_macro_particles_per_batch) then
      if (spec%target_macro_particles_per_batch == -1_i32) then
        w_particle = cfg%particle_species(1)%w_particle
      else
        number_density_m3 = species_number_density_m3(spec)
        temperature_k = species_temperature_k(spec)
        call resolve_inward_normal(spec%inject_face, inward_normal)
        gamma_in = compute_inflow_flux_from_drifting_maxwellian( &
          number_density_m3, temperature_k, spec%m_particle, spec%drift_velocity, inward_normal &
        )
        w_particle = gamma_in * area * cfg%sim%batch_duration / real(spec%target_macro_particles_per_batch, dp)
      end if
      if (.not. ieee_is_finite(w_particle) .or. w_particle <= 0.0d0) then
        error stop 'target_macro_particles_per_batch produced invalid w_particle.'
      end if
      spec%w_particle = w_particle
      spec%has_w_particle = .true.
    end if

    cfg%particle_species(species_idx) = spec
  end subroutine validate_reservoir_species

  !> `photo_raycast` 用の必須項目と整合性を検証する。
  !! @param[inout] cfg 検証対象のアプリ設定。
  !! @param[in] species_idx 検証する粒子種のインデックス（1始まり）。
  subroutine validate_photo_raycast_species(cfg, species_idx)
    type(app_config), intent(inout) :: cfg
    integer, intent(in) :: species_idx

    integer :: axis, axis_t1, axis_t2
    real(dp) :: boundary_value, area, direction_norm, inward_dot
    real(dp) :: inward_normal(3)
    type(particle_species_spec) :: spec

    spec = cfg%particle_species(species_idx)

    if (spec%has_npcls_per_step) then
      error stop 'particles.species.npcls_per_step is not allowed for photo_raycast.'
    end if
    if (spec%has_number_density_cm3 .or. spec%has_number_density_m3) then
      error stop 'number_density_cm3/number_density_m3 are not allowed for photo_raycast.'
    end if
    if (spec%has_w_particle) then
      error stop 'w_particle is not allowed for photo_raycast.'
    end if
    if (spec%has_target_macro_particles_per_batch) then
      error stop 'target_macro_particles_per_batch is not allowed for photo_raycast.'
    end if
    if (.not. cfg%sim%use_box) then
      error stop 'particles.species.source_mode="photo_raycast" requires sim.use_box = true.'
    end if
    if (cfg%sim%batch_duration <= 0.0d0) then
      error stop 'sim.batch_duration must be > 0 for photo_raycast.'
    end if
    if (.not. ieee_is_finite(spec%emit_current_density_a_m2) .or. spec%emit_current_density_a_m2 <= 0.0d0) then
      error stop 'photo_raycast requires emit_current_density_a_m2 > 0.'
    end if
    if (spec%rays_per_batch <= 0_i32) then
      error stop 'photo_raycast requires rays_per_batch > 0.'
    end if
    if (.not. ieee_is_finite(spec%normal_drift_speed)) then
      error stop 'normal_drift_speed must be finite.'
    end if
    if (spec%has_temperature_ev .and. spec%has_temperature_k) then
      error stop 'Specify either temperature_ev or temperature_k, not both.'
    end if
    if (spec%has_temperature_ev) then
      if (.not. ieee_is_finite(spec%temperature_ev) .or. spec%temperature_ev < 0.0d0) then
        error stop 'temperature_ev must be finite and >= 0.'
      end if
    else
      if (.not. ieee_is_finite(spec%temperature_k) .or. spec%temperature_k < 0.0d0) then
        error stop 'temperature_k must be finite and >= 0.'
      end if
    end if
    if (.not. ieee_is_finite(spec%m_particle) .or. spec%m_particle <= 0.0d0) then
      error stop 'm_particle must be finite and > 0.'
    end if
    if (.not. ieee_is_finite(spec%q_particle) .or. abs(spec%q_particle) <= 0.0d0) then
      error stop 'q_particle must be finite and non-zero for photo_raycast.'
    end if

    call resolve_inject_face(cfg%sim%box_min, cfg%sim%box_max, spec%inject_face, axis, boundary_value)
    axis_t1 = modulo(axis, 3) + 1
    axis_t2 = modulo(axis + 1, 3) + 1
    if (abs(spec%pos_low(axis) - boundary_value) > 1.0d-12 .or. abs(spec%pos_high(axis) - boundary_value) > 1.0d-12) then
      error stop 'photo_raycast pos_low/pos_high must lie on the selected box face.'
    end if
    if (spec%pos_high(axis_t1) < spec%pos_low(axis_t1) .or. spec%pos_high(axis_t2) < spec%pos_low(axis_t2)) then
      error stop 'photo_raycast tangential bounds must satisfy pos_high >= pos_low.'
    end if
    if ((spec%pos_high(axis_t1) - spec%pos_low(axis_t1)) <= 0.0d0 .or. &
        (spec%pos_high(axis_t2) - spec%pos_low(axis_t2)) <= 0.0d0) then
      error stop 'photo_raycast opening area must be positive.'
    end if
    area = compute_face_area_from_bounds(spec%inject_face, spec%pos_low, spec%pos_high)
    if (.not. ieee_is_finite(area) .or. area <= 0.0d0) then
      error stop 'photo_raycast opening area must be positive.'
    end if

    call resolve_inward_normal(spec%inject_face, inward_normal)
    if (spec%has_ray_direction) then
      if (.not. all(ieee_is_finite(spec%ray_direction))) then
        error stop 'ray_direction must be finite.'
      end if
      direction_norm = sqrt(sum(spec%ray_direction * spec%ray_direction))
      if (direction_norm <= 0.0d0) then
        error stop 'ray_direction norm must be > 0.'
      end if
      spec%ray_direction = spec%ray_direction / direction_norm
    else
      spec%ray_direction = inward_normal
    end if
    inward_dot = dot_product(spec%ray_direction, inward_normal)
    if (.not. ieee_is_finite(inward_dot) .or. inward_dot <= 0.0d0) then
      error stop 'ray_direction must point inward from inject_face.'
    end if

    cfg%particle_species(species_idx) = spec
  end subroutine validate_photo_raycast_species

  !> 注入面名から法線軸と対応する境界座標を返す。
  !! @param[in] box_min ボックス下限座標 `(x,y,z)` [m]。
  !! @param[in] box_max ボックス上限座標 `(x,y,z)` [m]。
  !! @param[in] inject_face 注入面識別子（`x_low/x_high/y_low/y_high/z_low/z_high`）。
  !! @param[out] axis 注入面の法線軸インデックス（1:x, 2:y, 3:z）。
  !! @param[out] boundary_value 指定面の境界座標値 [m]。
  subroutine resolve_inject_face(box_min, box_max, inject_face, axis, boundary_value)
    real(dp), intent(in) :: box_min(3), box_max(3)
    character(len=*), intent(in) :: inject_face
    integer, intent(out) :: axis
    real(dp), intent(out) :: boundary_value

    select case (trim(lower(inject_face)))
    case ('x_low')
      axis = 1
      boundary_value = box_min(1)
    case ('x_high')
      axis = 1
      boundary_value = box_max(1)
    case ('y_low')
      axis = 2
      boundary_value = box_min(2)
    case ('y_high')
      axis = 2
      boundary_value = box_max(2)
    case ('z_low')
      axis = 3
      boundary_value = box_min(3)
    case ('z_high')
      axis = 3
      boundary_value = box_max(3)
    case default
      error stop 'Unknown particles.species.inject_face.'
    end select
  end subroutine resolve_inject_face

  !> 注入面名から内向き法線ベクトルを返す。
  !! @param[in] inject_face 注入面識別子（`x_low/x_high/y_low/y_high/z_low/z_high`）。
  !! @param[out] inward_normal 注入面の内向き単位法線ベクトル。
  subroutine resolve_inward_normal(inject_face, inward_normal)
    character(len=*), intent(in) :: inject_face
    real(dp), intent(out) :: inward_normal(3)

    inward_normal = 0.0d0
    select case (trim(lower(inject_face)))
    case ('x_low')
      inward_normal(1) = 1.0d0
    case ('x_high')
      inward_normal(1) = -1.0d0
    case ('y_low')
      inward_normal(2) = 1.0d0
    case ('y_high')
      inward_normal(2) = -1.0d0
    case ('z_low')
      inward_normal(3) = 1.0d0
    case ('z_high')
      inward_normal(3) = -1.0d0
    case default
      error stop 'Unknown particles.species.inject_face.'
    end select
  end subroutine resolve_inward_normal

  !> 粒子種設定から実効密度[m^-3]を返す。
  !! @param[in] spec 粒子種設定。
  !! @return number_density_m3 実効粒子数密度 [1/m^3]。
  pure real(dp) function species_number_density_m3(spec) result(number_density_m3)
    type(particle_species_spec), intent(in) :: spec

    number_density_m3 = spec%number_density_m3
    if (spec%has_number_density_cm3) number_density_m3 = spec%number_density_cm3 * 1.0d6
  end function species_number_density_m3

  !> 粒子種設定から実効温度[K]を返す。
  !! @param[in] spec 粒子種設定。
  !! @return temperature_k 実効温度 [K]。
  pure real(dp) function species_temperature_k(spec) result(temperature_k)
    type(particle_species_spec), intent(in) :: spec

    temperature_k = spec%temperature_k
    if (spec%has_temperature_ev) temperature_k = spec%temperature_ev * 1.160451812d4
  end function species_temperature_k

  !> drifting Maxwellian の片側流入束 [#/m^2/s] を返す。
  !! @param[in] number_density_m3 粒子数密度 [1/m^3]。
  !! @param[in] temperature_k 温度 [K]。
  !! @param[in] m_particle 粒子1個あたりの質量 [kg]。
  !! @param[in] drift_velocity ドリフト速度ベクトル `(vx,vy,vz)` [m/s]。
  !! @param[in] inward_normal 注入面の内向き単位法線ベクトル。
  !! @return gamma_in 片側流入束 [1/m^2/s]。
  pure real(dp) function compute_inflow_flux_from_drifting_maxwellian( &
    number_density_m3, temperature_k, m_particle, drift_velocity, inward_normal &
  ) result(gamma_in)
    real(dp), intent(in) :: number_density_m3
    real(dp), intent(in) :: temperature_k
    real(dp), intent(in) :: m_particle
    real(dp), intent(in) :: drift_velocity(3)
    real(dp), intent(in) :: inward_normal(3)
    real(dp) :: sigma, alpha, u_n

    u_n = dot_product(drift_velocity, inward_normal)
    sigma = sqrt(k_boltzmann * temperature_k / m_particle)
    if (sigma <= 0.0d0) then
      gamma_in = number_density_m3 * max(0.0d0, u_n)
      return
    end if

    alpha = u_n / sigma
    gamma_in = number_density_m3 * (sigma * standard_normal_pdf(alpha) + u_n * standard_normal_cdf(alpha))
  end function compute_inflow_flux_from_drifting_maxwellian

  !> 標準正規分布の PDF を返す。
  !! @param[in] x 評価点。
  !! @return pdf 標準正規分布の確率密度。
  pure real(dp) function standard_normal_pdf(x) result(pdf)
    real(dp), intent(in) :: x
    real(dp), parameter :: inv_sqrt_2pi = 3.98942280401432678d-1

    pdf = inv_sqrt_2pi * exp(-0.5d0 * x * x)
  end function standard_normal_pdf

  !> 標準正規分布の CDF を返す。
  !! @param[in] x 評価点。
  !! @return cdf 標準正規分布の累積分布値。
  pure real(dp) function standard_normal_cdf(x) result(cdf)
    real(dp), intent(in) :: x
    real(dp), parameter :: inv_sqrt_2 = 7.07106781186547524d-1

    cdf = 0.5d0 * (1.0d0 + erf(x * inv_sqrt_2))
  end function standard_normal_cdf

  !> 注入面上の矩形開口から有効面積[m^2]を返す。
  !! @param[in] inject_face 注入面識別子（`x_low/x_high/y_low/y_high/z_low/z_high`）。
  !! @param[in] pos_low 開口領域の下限座標 `(x,y,z)` [m]。
  !! @param[in] pos_high 開口領域の上限座標 `(x,y,z)` [m]。
  !! @return area 注入開口の有効面積 [m^2]。
  pure real(dp) function compute_face_area_from_bounds(inject_face, pos_low, pos_high) result(area)
    character(len=*), intent(in) :: inject_face
    real(dp), intent(in) :: pos_low(3), pos_high(3)
    integer :: axis_t1, axis_t2

    call resolve_face_axes(inject_face, axis_t1, axis_t2)
    area = (pos_high(axis_t1) - pos_low(axis_t1)) * (pos_high(axis_t2) - pos_low(axis_t2))
  end function compute_face_area_from_bounds

  !> 注入面名から接線2軸を返す。
  !! @param[in] inject_face 注入面識別子（`x_low/x_high/y_low/y_high/z_low/z_high`）。
  !! @param[out] axis_t1 注入面の第1接線軸インデックス（1:x, 2:y, 3:z）。
  !! @param[out] axis_t2 注入面の第2接線軸インデックス（1:x, 2:y, 3:z）。
  pure subroutine resolve_face_axes(inject_face, axis_t1, axis_t2)
    character(len=*), intent(in) :: inject_face
    integer, intent(out) :: axis_t1, axis_t2

    select case (trim(lower(inject_face)))
    case ('x_low', 'x_high')
      axis_t1 = 2
      axis_t2 = 3
    case ('y_low', 'y_high')
      axis_t1 = 3
      axis_t2 = 1
    case ('z_low', 'z_high')
      axis_t1 = 1
      axis_t2 = 2
    case default
      error stop 'Unknown particles.species.inject_face.'
    end select
  end subroutine resolve_face_axes

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
      error stop ('Unknown key in [mesh]: ' // trim(k))
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
    case ('height')
      call parse_real(v, spec%height)
    case ('n_theta')
      call parse_int(v, spec%n_theta)
    case ('n_z')
      call parse_int(v, spec%n_z)
    case ('cap')
      call parse_logical(v, spec%cap)
    case ('n_lon')
      call parse_int(v, spec%n_lon)
    case ('n_lat')
      call parse_int(v, spec%n_lat)
    case default
      error stop ('Unknown key in [[mesh.templates]]: ' // trim(k))
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
    case ('dir')
      call parse_string(v, cfg%output_dir)
    case ('history_stride')
      call parse_int(v, cfg%history_stride)
    case ('resume')
      call parse_logical(v, cfg%resume_output)
    case default
      error stop ('Unknown key in [output]: ' // trim(k))
    end select
  end subroutine apply_output_kv

  !> `key = value` 形式の1行を分割し、前後空白を取り除いて返す。
  !! @param[in] line 分割対象の入力行。
  !! @param[out] key 正規化したキー文字列（小文字化済み）。
  !! @param[out] value 前後空白を除去した値文字列。
  subroutine split_key_value(line, key, value)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: key
    character(len=*), intent(out) :: value
    integer :: p

    p = index(line, '=')
    if (p <= 0) then
      key = ''
      value = ''
      return
    end if
    key = lower(trim(adjustl(line(:p - 1))))
    value = trim(adjustl(line(p + 1:)))
  end subroutine split_key_value

  !> 文字列表現を倍精度実数へ変換する。
  !! @param[in] text 実数を表す文字列。
  !! @param[out] out 変換後の倍精度実数。
  subroutine parse_real(text, out)
    character(len=*), intent(in) :: text
    real(dp), intent(out) :: out

    read (text, *) out
  end subroutine parse_real

  !> 文字列表現を 32bit 整数へ変換する。
  !! @param[in] text 整数を表す文字列。
  !! @param[out] out 変換後の32bit整数。
  subroutine parse_int(text, out)
    character(len=*), intent(in) :: text
    integer(i32), intent(out) :: out

    read (text, *) out
  end subroutine parse_int

  !> `true/.true.` を真とみなして論理値へ変換する。
  !! @param[in] text 論理値文字列。
  !! @param[out] out 変換後の論理値。
  subroutine parse_logical(text, out)
    character(len=*), intent(in) :: text
    logical, intent(out) :: out
    character(len=32) :: t

    t = lower(trim(adjustl(text)))
    select case (t)
    case ('true', '.true.')
      out = .true.
    case ('false', '.false.')
      out = .false.
    case default
      error stop 'Logical value must be true/false (or .true./.false.).'
    end select
  end subroutine parse_logical

  !> 前後の二重引用符を外して文字列値へ変換する。
  !! @param[in] text 文字列リテラルまたは裸文字列。
  !! @param[out] out 引用符を除去した文字列値。
  subroutine parse_string(text, out)
    character(len=*), intent(in) :: text
    character(len=*), intent(out) :: out
    character(len=:), allocatable :: tmp

    tmp = trim(adjustl(text))
    if (len(tmp) >= 2 .and. tmp(1:1) == '"' .and. tmp(len(tmp):len(tmp)) == '"') then
      tmp = tmp(2:len(tmp) - 1)
    end if
    out = trim(tmp)
  end subroutine parse_string

  !> `[x, y, z]` 形式の3成分ベクトルを配列へ変換する。
  !! @param[in] text 3成分ベクトル文字列。
  !! @param[out] out 変換後の3成分実数ベクトル。
  subroutine parse_real3(text, out)
    character(len=*), intent(in) :: text
    real(dp), intent(out) :: out(3)
    character(len=256) :: t

    t = trim(adjustl(text))
    if (t(1:1) == '[') t = t(2:)
    if (t(len_trim(t):len_trim(t)) == ']') t = t(:len_trim(t) - 1)
    read (t, *) out(1), out(2), out(3)
  end subroutine parse_real3

  !> 境界条件の文字列を内部定数へ変換する。
  !! 未知のモードは実行継続せず停止する。
  !! @param[in] text 境界条件モード文字列。
  !! @param[out] out 変換後の境界条件定数。
  subroutine parse_boundary_mode(text, out)
    character(len=*), intent(in) :: text
    integer(i32), intent(out) :: out
    character(len=64) :: mode

    call parse_string(text, mode)
    select case (trim(lower(mode)))
    case ('open', 'outflow', 'escape')
      out = bc_open
    case ('reflect', 'reflection')
      out = bc_reflect
    case ('periodic')
      out = bc_periodic
    case default
      error stop 'Unknown boundary condition mode in [sim].'
    end select
  end subroutine parse_boundary_mode

  !> `#` 以降の行内コメントを除去する。
  !! @param[in] line コメント除去対象の入力行。
  !! @return out `#` 以降を取り除いた行文字列。
  pure function strip_comment(line) result(out)
    character(len=*), intent(in) :: line
    character(len=len(line)) :: out
    integer :: p

    p = index(line, '#')
    if (p > 0) then
      out = trim(line(:p - 1))
    else
      out = trim(line)
    end if
  end function strip_comment

  !> ASCII 英字だけを小文字化し、設定キー比較を簡単にする。
  !! @param[in] s 変換対象文字列。
  !! @return o 小文字化した文字列。
  pure function lower(s) result(o)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: o
    integer :: i, c

    o = s
    do i = 1, len(s)
      c = iachar(s(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) o(i:i) = achar(c + 32)
    end do
  end function lower

  !> 文字列 `s` が `suffix` で終わるかを返す。
  !! @param[in] s 判定対象文字列。
  !! @param[in] suffix 末尾一致を判定する接尾辞。
  !! @return ends_with `s` が `suffix` で終われば `.true.`。
  pure logical function ends_with(s, suffix)
    character(len=*), intent(in) :: s
    character(len=*), intent(in) :: suffix
    integer :: ls, lf

    ls = len_trim(s)
    lf = len_trim(suffix)
    if (lf > ls) then
      ends_with = .false.
    else
      ends_with = (s(ls - lf + 1:ls) == suffix(1:lf))
    end if
  end function ends_with

end module bem_app_config_parser

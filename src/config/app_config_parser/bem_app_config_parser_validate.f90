!> `bem_app_config_parser` の入力検証・物理量導出手続きを実装する submodule。
submodule(bem_app_config_parser) bem_app_config_parser_validate
  use bem_config_helpers, only: resolve_inject_face, resolve_inward_normal, species_number_density_m3, species_temperature_k
  implicit none
contains

  !> `batch_duration` 系パラメータの排他条件と値域を検証し、確定値を反映する。
  module procedure resolve_batch_duration
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
    batch_duration = cfg%sim%dt*cfg%sim%batch_duration_step
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
  end procedure resolve_batch_duration

  !> `reservoir_face` 粒子種の必須項目と幾何条件を検証する。
  module procedure validate_reservoir_species
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
      if (trim(lower_ascii(cfg%particle_species(1)%source_mode)) /= 'reservoir_face') then
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
      w_particle = gamma_in*area*cfg%sim%batch_duration/real(spec%target_macro_particles_per_batch, dp)
    end if
    if (.not. ieee_is_finite(w_particle) .or. w_particle <= 0.0d0) then
      error stop 'target_macro_particles_per_batch produced invalid w_particle.'
    end if
    spec%w_particle = w_particle
    spec%has_w_particle = .true.
  end if

  cfg%particle_species(species_idx) = spec
  end procedure validate_reservoir_species

  !> `photo_raycast` 粒子種の必須項目と幾何/方向条件を検証する。
  module procedure validate_photo_raycast_species
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
    direction_norm = sqrt(sum(spec%ray_direction*spec%ray_direction))
    if (direction_norm <= 0.0d0) then
      error stop 'ray_direction norm must be > 0.'
    end if
    spec%ray_direction = spec%ray_direction/direction_norm
  else
    spec%ray_direction = inward_normal
  end if
  inward_dot = dot_product(spec%ray_direction, inward_normal)
  if (.not. ieee_is_finite(inward_dot) .or. inward_dot <= 0.0d0) then
    error stop 'ray_direction must point inward from inject_face.'
  end if

  cfg%particle_species(species_idx) = spec
  end procedure validate_photo_raycast_species

  !> drifting Maxwellian の片側流入束 `[1/m^2/s]` を評価する。
  module procedure compute_inflow_flux_from_drifting_maxwellian
  real(dp) :: sigma, alpha, u_n

  u_n = dot_product(drift_velocity, inward_normal)
  sigma = sqrt(k_boltzmann*temperature_k/m_particle)
  if (sigma <= 0.0d0) then
    gamma_in = number_density_m3*max(0.0d0, u_n)
    return
  end if

  alpha = u_n/sigma
  gamma_in = number_density_m3*(sigma*standard_normal_pdf(alpha) + u_n*standard_normal_cdf(alpha))
  end procedure compute_inflow_flux_from_drifting_maxwellian

  !> 標準正規分布の PDF を評価する。
  module procedure standard_normal_pdf
  real(dp), parameter :: inv_sqrt_2pi = 3.98942280401432678d-1

  pdf = inv_sqrt_2pi*exp(-0.5d0*x*x)
  end procedure standard_normal_pdf

  !> 標準正規分布の CDF を評価する。
  module procedure standard_normal_cdf
  real(dp), parameter :: inv_sqrt_2 = 7.07106781186547524d-1

  cdf = 0.5d0*(1.0d0 + erf(x*inv_sqrt_2))
  end procedure standard_normal_cdf

  !> 注入面上開口矩形の有効面積 `[m^2]` を計算する。
  module procedure compute_face_area_from_bounds
  integer :: axis_t1, axis_t2

  call resolve_face_axes(inject_face, axis_t1, axis_t2)
  area = (pos_high(axis_t1) - pos_low(axis_t1))*(pos_high(axis_t2) - pos_low(axis_t2))
  end procedure compute_face_area_from_bounds

  !> 注入面識別子から接線2軸インデックスを返す。
  module procedure resolve_face_axes
  select case (trim(lower_ascii(inject_face)))
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
  end procedure resolve_face_axes

end submodule bem_app_config_parser_validate

!> 有効粒子種の境界・電流・粒子源の検証と workload 条件。
submodule(bem_app_config_parser) bem_app_config_parser_preflight_particles
  use bem_config_helpers, only: resolve_particle_boundaries, particle_boundary_action_for_face
  use bem_app_config_types, only: particle_inflow_none, particle_inflow_reservoir
  implicit none
contains

  module procedure validate_particle_species_config
  integer :: i, j, axis
  integer(i32) :: effective_boundary_low(3), effective_boundary_high(3), inject_face_boundary
  character(len=64) :: generated_species_key
  logical :: has_boundary_inflow
  per_batch_particles = 0_i32
  has_dynamic_source_species = .false.
  has_enabled_volume_seed = .false.
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
    cfg%particle_species(i)%surface_charge_closure = &
      lower_ascii(trim(cfg%particle_species(i)%surface_charge_closure))
    has_boundary_inflow = any(cfg%particle_species(i)%boundary_inflow_low /= particle_inflow_none) .or. &
                          any(cfg%particle_species(i)%boundary_inflow_high /= particle_inflow_none)
    do axis = 1, 3
      call validate_particle_boundary_override( &
        cfg%particle_species(i)%boundary_low(axis), cfg%sim%bc_low(axis), &
        'particles.species.boundary low face' &
        )
      call validate_particle_boundary_override( &
        cfg%particle_species(i)%boundary_high(axis), cfg%sim%bc_high(axis), &
        'particles.species.boundary high face' &
        )
    end do
    if (.not. cfg%sim%use_box .and. &
        (any(cfg%particle_species(i)%boundary_low /= particle_bc_inherit) .or. &
         any(cfg%particle_species(i)%boundary_high /= particle_bc_inherit))) then
      error stop 'particles.species.boundary requires a finite [domain].'
    end if
    call resolve_particle_boundaries( &
      cfg%sim, cfg%particle_boundary_low, cfg%particle_boundary_high, cfg%particle_species(i), &
      effective_boundary_low, effective_boundary_high &
      )
    do axis = 1, 3
      call validate_particle_boundary_inflow( &
        cfg%particle_species(i)%boundary_inflow_low(axis), cfg%sim%bc_low(axis), effective_boundary_low(axis), &
        'particles.species.boundary_inflow low face' &
        )
      call validate_particle_boundary_inflow( &
        cfg%particle_species(i)%boundary_inflow_high(axis), cfg%sim%bc_high(axis), effective_boundary_high(axis), &
        'particles.species.boundary_inflow high face' &
        )
    end do
    select case (trim(cfg%particle_species(i)%surface_charge_closure))
    case ('explicit')
      if (cfg%particle_species(i)%has_target_absorbed_current_a .or. &
          cfg%particle_species(i)%has_target_emission_current_a) then
        error stop 'target surface currents require surface_charge_closure="fixed_current".'
      end if
    case ('fixed_current')
      if (.not. ieee_is_finite(cfg%sim%batch_duration) .or. cfg%sim%batch_duration <= 0.0_dp) then
        error stop 'sim.batch_duration must be > 0 for fixed_current.'
      end if
      if (.not. cfg%particle_species(i)%has_target_absorbed_current_a .and. &
          .not. cfg%particle_species(i)%has_target_emission_current_a .and. &
          .not. is_automatic_current_species(cfg, i)) then
        error stop 'surface_charge_closure="fixed_current" requires at least one target current.'
      end if
      if (cfg%particle_species(i)%has_target_absorbed_current_a) then
        if (.not. ieee_is_finite(cfg%particle_species(i)%target_absorbed_current_a)) then
          error stop 'particles.species.target_absorbed_current_a must be finite.'
        end if
        if (cfg%particle_species(i)%target_absorbed_current_a*cfg%particle_species(i)%q_particle < 0.0_dp) then
          error stop 'target_absorbed_current_a sign must match q_particle.'
        end if
      end if
      if (cfg%particle_species(i)%has_target_emission_current_a) then
        if (trim(cfg%particle_species(i)%source_mode) /= 'photo_raycast' .or. &
            .not. cfg%particle_species(i)%deposit_opposite_charge_on_emit) then
          error stop 'target_emission_current_a requires photo_raycast with opposite-charge emission deposit.'
        end if
        if (.not. ieee_is_finite(cfg%particle_species(i)%target_emission_current_a)) then
          error stop 'particles.species.target_emission_current_a must be finite.'
        end if
        if (cfg%particle_species(i)%target_emission_current_a*cfg%particle_species(i)%q_particle > 0.0_dp) then
          error stop 'target_emission_current_a sign must oppose q_particle.'
        end if
      end if
    case ('neutral_return')
      if (trim(cfg%particle_species(i)%source_mode) /= 'photo_raycast' .or. &
          cfg%particle_species(i)%q_particle >= 0.0_dp) then
        error stop 'surface_charge_closure="neutral_return" requires a negative photo_raycast species.'
      end if
      if (.not. cfg%particle_species(i)%deposit_opposite_charge_on_emit) then
        error stop 'surface_charge_closure="neutral_return" requires deposit_opposite_charge_on_emit=true.'
      end if
      inject_face_boundary = particle_boundary_action_for_face( &
                             effective_boundary_low, effective_boundary_high, cfg%particle_species(i)%inject_face &
                             )
      if (inject_face_boundary /= bc_reflect .and. inject_face_boundary /= bc_redistributed_reflect) then
        error stop 'surface_charge_closure="neutral_return" requires a reflecting action on the species inject_face.'
      end if
    case default
      error stop 'particles.species.surface_charge_closure must be "explicit", "fixed_current", or "neutral_return".'
    end select
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
      if (cfg%particle_species(i)%has_source_normal) then
        error stop 'source_normal is only valid for source_mode="plane_source".'
      end if
      has_enabled_volume_seed = has_enabled_volume_seed .or. cfg%particle_species(i)%npcls_per_step > 0_i32
      if (cfg%particle_species(i)%npcls_per_step < 0_i32) then
        error stop 'particles.species.npcls_per_step must be >= 0.'
      end if
      if (.not. has_boundary_inflow .and. &
          (trim(cfg%particle_species(i)%velocity_distribution) /= 'maxwellian' .or. &
           len_trim(cfg%particle_species(i)%velocity_grid_path) > 0 .or. &
           trim(cfg%particle_species(i)%velocity_grid_sampling) /= 'auto' .or. &
           cfg%particle_species(i)%has_particle_flux_m2_s .or. cfg%particle_species(i)%has_current_density_a_m2)) then
        error stop 'velocity_distribution="grid" and flux keys are only valid for reservoir_face.'
      end if
      if (.not. has_boundary_inflow .and. cfg%particle_species(i)%has_target_macro_particles_per_batch) then
        error stop 'target_macro_particles_per_batch is only valid for reservoir_face.'
      end if
      if (abs(cfg%particle_species(i)%emit_current_density_a_m2) > 0.0d0 .or. &
          cfg%particle_species(i)%rays_per_batch /= 0_i32 .or. cfg%particle_species(i)%has_ray_direction .or. &
          cfg%particle_species(i)%has_deposit_opposite_charge_on_emit) then
        error stop 'photo_raycast keys are only valid for source_mode="photo_raycast".'
      end if
      per_batch_particles = per_batch_particles + cfg%particle_species(i)%npcls_per_step
      if (has_boundary_inflow) then
        has_dynamic_source_species = .true.
        call validate_boundary_inflow_species(cfg, i)
      end if
    case ('reservoir_face')
      if (cfg%particle_species(i)%has_source_normal) then
        error stop 'source_normal is only valid for source_mode="plane_source".'
      end if
      if (has_boundary_inflow) then
        error stop 'source_mode="reservoir_face" cannot be combined with boundary_inflow.'
      end if
      has_dynamic_source_species = .true.
      call validate_reservoir_species(cfg, i)
    case ('plane_source')
      if (has_boundary_inflow) then
        error stop 'source_mode="plane_source" cannot be combined with boundary_inflow.'
      end if
      has_dynamic_source_species = .true.
      call validate_plane_source_species(cfg, i)
    case ('photo_raycast')
      if (cfg%sim%raycast_max_bounce < 1_i32) then
        error stop 'sim.raycast_max_bounce must be >= 1 when photo_raycast is enabled.'
      end if
      if (cfg%particle_species(i)%has_source_normal) then
        error stop 'source_normal is only valid for source_mode="plane_source".'
      end if
      if (has_boundary_inflow) then
        error stop 'source_mode="photo_raycast" cannot be combined with boundary_inflow.'
      end if
      has_dynamic_source_species = .true.
      call validate_photo_raycast_species(cfg, i)
    case default
      error stop 'Unknown particles.species.source_mode.'
    end select
  end do
  end procedure validate_particle_species_config

  module procedure validate_source_workload
  integer :: i
  if (per_batch_particles <= 0_i32 .and. .not. has_dynamic_source_species) then
    error stop 'At least one enabled [[particles.species]] entry must have npcls_per_step > 0.'
  end if
  if (cfg%periodic2%max_nonzero_mode_potential_step > 0.0_dp) then
    if (trim(lower_ascii(cfg%periodic2%nonzero_mode_backend)) /= 'cached_kneq0') then
      error stop 'periodic2.max_nonzero_mode_potential_step requires nonzero_mode_backend="cached_kneq0".'
    end if
    if (.not. ieee_is_finite(cfg%sim%batch_duration) .or. cfg%sim%batch_duration <= 0.0_dp) then
      error stop 'periodic2.max_nonzero_mode_potential_step requires a positive sim.batch_duration.'
    end if
    if (has_enabled_volume_seed) then
      error stop 'periodic2.max_nonzero_mode_potential_step requires time-scaled reservoir_face/photo_raycast sources.'
    end if
    do i = 1_i32, cfg%n_particle_species
      if (.not. cfg%particle_species(i)%enabled) cycle
      if (trim(lower_ascii(cfg%particle_species(i)%source_mode)) /= 'reservoir_face' .and. &
          trim(lower_ascii(cfg%particle_species(i)%source_mode)) /= 'plane_source' .and. &
          .not. any(cfg%particle_species(i)%boundary_inflow_low == particle_inflow_reservoir) .and. &
          .not. any(cfg%particle_species(i)%boundary_inflow_high == particle_inflow_reservoir)) cycle
      if (.not. cfg%particle_species(i)%has_target_macro_particles_per_batch) then
        error stop 'adaptive flux-driven injection requires target_macro_particles_per_batch instead of fixed w_particle.'
      end if
    end do
  end if
  end procedure validate_source_workload

  module procedure validate_particle_boundary_override

  select case (action)
  case (particle_bc_inherit, bc_open, bc_reflect, bc_redistributed_reflect)
    continue
  case default
    call stop_config_error(trim(context)//' must be inherit, open, reflect, or redistributed_reflect.')
  end select
  if (topology_action == bc_periodic .and. action /= particle_bc_inherit) then
    call stop_config_error(trim(context)//' cannot override a periodic domain face.')
  end if
  end procedure validate_particle_boundary_override

  module procedure validate_particle_boundary_inflow

  select case (inflow)
  case (particle_inflow_none)
    return
  case (particle_inflow_reservoir)
    continue
  case default
    call stop_config_error(trim(context)//' must be none or reservoir.')
  end select
  if (effective_topology_action == bc_periodic) then
    call stop_config_error(trim(context)//' cannot inject through a periodic domain face.')
  end if
  if (effective_particle_action /= bc_open) then
    call stop_config_error(trim(context)//' requires the effective particle action to be open.')
  end if
  end procedure validate_particle_boundary_inflow

end submodule bem_app_config_parser_preflight_particles

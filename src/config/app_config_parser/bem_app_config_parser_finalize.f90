!> 読み込み済み TOML 設定の正規化・派生値確定・検証を実装する submodule。
submodule(bem_app_config_parser) bem_app_config_parser_finalize
  use bem_constants, only: qe
  use bem_config_helpers, only: resolve_particle_boundaries, particle_boundary_action_for_face, &
                                species_number_density_m3, species_temperature_k
  use bem_app_config_types, only: particle_inflow_none, particle_inflow_reservoir
  implicit none
contains

  module procedure finalize_loaded_config
  integer :: i, j, axis
  integer(i32) :: per_batch_particles, physics_status
  integer(i32) :: n_periodic_axes
  integer(i32) :: effective_boundary_low(3), effective_boundary_high(3), inject_face_boundary
  logical :: has_dynamic_source_species, has_enabled_volume_seed, adaptive_nonzero_mode, has_boundary_inflow
  character(len=64) :: generated_species_key
  character(len=256) :: physics_message
  type(field_physics_config) :: field_config
  type(panel_kernel_config) :: panel_config

  call lower_boundary_authoring(cfg, authoring)
  call normalize_high_level_config(cfg, authoring)
  call normalize_legacy_physics_config(cfg%sim, cfg%field, cfg%periodic2, cfg%panel)
  call apply_physics_authoring(cfg, authoring)

  if (.not. ieee_is_finite(cfg%periodic2%max_nonzero_mode_potential_step) .or. &
      cfg%periodic2%max_nonzero_mode_potential_step < 0.0_dp) then
    error stop 'periodic2.max_nonzero_mode_potential_step must be finite and >= 0.'
  end if
  adaptive_nonzero_mode = cfg%periodic2%max_nonzero_mode_potential_step > 0.0_dp

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
  if (.not. all(ieee_is_finite(cfg%sim%b0))) error stop 'sim.b0 must contain finite values.'
  if (cfg%sim%use_box) then
    if (.not. all(ieee_is_finite(cfg%sim%box_min)) .or. .not. all(ieee_is_finite(cfg%sim%box_max))) then
      error stop 'sim.box_min/box_max must contain finite values.'
    end if
    if (any(cfg%sim%box_max <= cfg%sim%box_min)) then
      error stop 'sim.box_max must be greater than sim.box_min on all axes when sim.use_box=true.'
    end if
  end if
  do axis = 1, 3
    call validate_particle_boundary_override( &
      cfg%particle_boundary_low(axis), cfg%sim%bc_low(axis), 'particle_boundary low face' &
      )
    call validate_particle_boundary_override( &
      cfg%particle_boundary_high(axis), cfg%sim%bc_high(axis), 'particle_boundary high face' &
      )
  end do
  if (.not. cfg%sim%use_box .and. &
      (any(cfg%particle_boundary_low /= particle_bc_inherit) .or. &
       any(cfg%particle_boundary_high /= particle_bc_inherit))) then
    error stop '[particle_boundary] requires a finite [domain].'
  end if
  if (cfg%n_particle_species <= 0_i32) error stop 'At least one [[particles.species]] entry is required.'
  if (len_trim(cfg%output_restart_from) > 0 .and. .not. cfg%resume_output) then
    error stop 'output.restart_from requires output.resume = true.'
  end if
  if (cfg%checkpoint_stride < 0_i32) then
    error stop 'output.checkpoint_stride must be >= 0.'
  end if
  if (cfg%checkpoint_stride > 0_i32 .and. .not. cfg%write_output) then
    error stop 'output.checkpoint_stride > 0 requires output.write_files = true.'
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
  case ('insulator', 'conductor')
    continue
  case ('dielectric')
    error stop 'mesh.surface_model="dielectric" is not implemented; use "insulator" for charge accumulation.'
  case default
    error stop 'mesh.surface_model must be "insulator" or "conductor".'
  end select
  do i = 1, cfg%n_templates
    cfg%templates(i)%surface_model = lower_ascii(trim(cfg%templates(i)%surface_model))
    select case (trim(cfg%templates(i)%surface_model))
    case ('insulator', 'conductor')
      continue
    case ('dielectric')
      error stop 'mesh.templates.surface_model="dielectric" is not implemented; use "insulator".'
    case default
      error stop 'mesh.templates.surface_model must be "insulator" or "conductor".'
    end select
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
    error stop 'sim.field_periodic_far_correction="m2l_root_oracle" was removed; use "cached_kneq0".'
  case ('cached_kneq0')
    continue
  case default
    error stop 'sim.field_periodic_far_correction must be "auto", "none", or "cached_kneq0".'
  end select
  if (trim(cfg%sim%field_periodic_far_correction) == 'cached_kneq0') then
    if (trim(cfg%sim%field_solver) /= 'fmm' .or. trim(cfg%sim%field_bc_mode) /= 'periodic2') then
      error stop 'sim.field_periodic_far_correction requires sim.field_solver="fmm" and sim.field_bc_mode="periodic2".'
    end if
    if (cfg%sim%field_periodic_ewald_layers < 1_i32) then
      error stop 'sim.field_periodic_ewald_layers must be >= 1 when far correction is enabled.'
    end if
  end if
  if (trim(cfg%sim%field_bc_mode) == 'periodic2') then
    if (cfg%sim%field_periodic_image_layers < 0_i32) then
      error stop 'sim.field_periodic_image_layers must be >= 0 for periodic2.'
    end if
  end if
  if (trim(cfg%sim%field_periodic_far_correction) == 'cached_kneq0') then
    if (cfg%sim%field_periodic_image_layers < 1_i32) then
      error stop 'cached_kneq0 requires sim.field_periodic_image_layers >= 1.'
    end if
    if (.not. ieee_is_finite(cfg%sim%field_periodic_ewald_alpha) .or. cfg%sim%field_periodic_ewald_alpha < 0.0d0) then
      error stop 'sim.field_periodic_ewald_alpha must be finite and >= 0 for cached_kneq0.'
    end if
    if (len_trim(cfg%sim%field_periodic_cache_dir) == 0) then
      error stop 'sim.field_periodic_cache_dir must not be empty for cached_kneq0.'
    end if
    if (.not. ieee_is_finite(cfg%sim%field_periodic_generation_tolerance) .or. &
        cfg%sim%field_periodic_generation_tolerance <= 0.0_dp) then
      error stop 'sim.field_periodic_generation_tolerance must be finite and > 0 for cached_kneq0.'
    end if
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
  if (trim(cfg%sim%field_solver) /= 'direct') then
    if (.not. ieee_is_finite(cfg%sim%tree_theta) .or. cfg%sim%tree_theta <= 0.0d0 .or. cfg%sim%tree_theta > 1.0d0) then
      error stop 'sim.tree_theta must be finite and satisfy 0 < theta <= 1 for tree-capable solvers.'
    end if
    if (cfg%sim%tree_leaf_max < 1_i32) then
      error stop 'sim.tree_leaf_max must be >= 1 for tree-capable solvers.'
    end if
    if (cfg%sim%tree_min_nelem < 1_i32) then
      error stop 'sim.tree_min_nelem must be >= 1 for tree-capable solvers.'
    end if
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
  cfg%sim%multiple_box_events_policy = lower_ascii(trim(cfg%sim%multiple_box_events_policy))
  select case (trim(cfg%sim%multiple_box_events_policy))
  case ('abort', 'soft_discard')
    continue
  case default
    error stop 'sim.multiple_box_events_policy must be "abort" or "soft_discard".'
  end select
  if (trim(cfg%sim%multiple_box_events_policy) == 'soft_discard') then
    if (cfg%sim%multiple_box_events_soft_discard_count_limit < 1_i32) then
      error stop 'sim.multiple_box_events_soft_discard_count_limit must be >= 1 for soft_discard.'
    end if
    if (.not. ieee_is_finite(cfg%sim%multiple_box_events_soft_discard_abs_charge_limit) .or. &
        cfg%sim%multiple_box_events_soft_discard_abs_charge_limit <= 0.0_dp) then
      error stop 'sim.multiple_box_events_soft_discard_abs_charge_limit must be finite and > 0 for soft_discard.'
    end if
  end if
  if (cfg%sim%injection_face_phi_grid_n < 1_i32) then
    error stop 'sim.injection_face_phi_grid_n must be >= 1.'
  end if
  if (.not. ieee_is_finite(cfg%sim%phi_infty)) then
    error stop 'sim.phi_infty must be finite.'
  end if
  call resolve_batch_duration(cfg)
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

  call validate_surface_current_model_config(cfg)

  if (per_batch_particles <= 0_i32 .and. .not. has_dynamic_source_species) then
    error stop 'At least one enabled [[particles.species]] entry must have npcls_per_step > 0.'
  end if
  if (adaptive_nonzero_mode) then
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
  call derive_field_panel_config(cfg%sim, field_config, panel_config)
  call validate_active_physics_config(cfg%sim, field_config, cfg%periodic2, panel_config, physics_status, physics_message)
  if (physics_status /= physics_config_ok) error stop trim(physics_message)
  end procedure finalize_loaded_config

  logical function is_automatic_current_species(cfg, species_idx) result(selected)
    type(app_config), intent(in) :: cfg
    integer, intent(in) :: species_idx

    selected = .false.
    if (trim(lower_ascii(cfg%surface_current%model)) == 'none') return
    selected = trim(cfg%particle_species(species_idx)%species_key) == trim(cfg%surface_current%electron_species) .or. &
               trim(cfg%particle_species(species_idx)%species_key) == trim(cfg%surface_current%ion_species) .or. &
               trim(cfg%particle_species(species_idx)%species_key) == trim(cfg%surface_current%photoelectron_species)
  end function is_automatic_current_species

  subroutine validate_surface_current_model_config(cfg)
    type(app_config), intent(in) :: cfg
    integer :: electron_idx, ion_idx, photo_idx
    integer(i32) :: effective_boundary_low(3), effective_boundary_high(3)
    real(dp) :: electron_temperature, ion_temperature
    logical :: photoelectron_active

    select case (trim(lower_ascii(cfg%surface_current%model)))
    case ('none')
      return
    case ('zhao_stationary')
      continue
    case default
      error stop 'surface_current_model.model must be "none" or "zhao_stationary".'
    end select
    if (any(cfg%sim%b0 /= 0.0_dp)) then
      error stop 'surface_current_model="zhao_stationary" requires sim.b0=[0,0,0] for its unmagnetized sheath closure.'
    end if
    if (trim(lower_ascii(cfg%sim%reservoir_potential_model)) /= 'none') then
      error stop 'Zhao kinetic inflow cannot be combined with the generic reservoir potential model.'
    end if
    select case (trim(lower_ascii(cfg%surface_current%zhao_branch)))
    case ('auto', 'a', 'b', 'c')
      continue
    case default
      error stop 'surface_current_model.zhao_branch must be "auto", "a", "b", or "c".'
    end select
    if (.not. ieee_is_finite(cfg%surface_current%photoelectron_source_scale) .or. &
        cfg%surface_current%photoelectron_source_scale < 0.0_dp) then
      error stop 'surface_current_model.photoelectron_source_scale must be finite and >= 0.'
    end if
    photoelectron_active = cfg%surface_current%photoelectron_source_scale > 0.0_dp
    if (.not. photoelectron_active .and. &
        trim(lower_ascii(cfg%surface_current%zhao_branch)) /= 'auto' .and. &
        trim(lower_ascii(cfg%surface_current%zhao_branch)) /= 'c') then
      error stop 'photoelectron_source_scale=0 requires surface_current_model.zhao_branch="auto" or "c".'
    end if
    if (photoelectron_active) then
      if (.not. ieee_is_finite(cfg%surface_current%solar_elevation_deg) .or. &
          cfg%surface_current%solar_elevation_deg <= 0.0_dp .or. &
          cfg%surface_current%solar_elevation_deg > 90.0_dp) then
        error stop 'surface_current_model.solar_elevation_deg must be finite and in (0, 90].'
      end if
      if (.not. ieee_is_finite(cfg%surface_current%photoelectron_ref_density_m3) .or. &
          cfg%surface_current%photoelectron_ref_density_m3 <= 0.0_dp) then
        error stop 'surface_current_model.photoelectron_ref_density_m3 must be finite and > 0.'
      end if
      if (len_trim(cfg%surface_current%photoelectron_species) == 0) then
        error stop 'positive photoelectron_source_scale requires surface_current_model.photoelectron_species.'
      end if
    else
      if (cfg%surface_current%has_photoelectron_species .or. &
          cfg%surface_current%has_solar_elevation_deg .or. &
          cfg%surface_current%has_photoelectron_ref_density_m3 .or. &
          len_trim(cfg%surface_current%photoelectron_species) > 0 .or. &
          cfg%surface_current%solar_elevation_deg /= 0.0_dp .or. &
          cfg%surface_current%photoelectron_ref_density_m3 /= 0.0_dp) then
        error stop 'photoelectron_source_scale=0 requires omitting all photoelectron-specific Zhao settings.'
      end if
    end if
    if (cfg%surface_current%has_reference_area_m2) then
      if (.not. ieee_is_finite(cfg%surface_current%reference_area_m2) .or. &
          cfg%surface_current%reference_area_m2 <= 0.0_dp) then
        error stop 'surface_current_model.reference_area_m2 must be finite and > 0.'
      end if
    else if (.not. cfg%sim%use_box .or. any(cfg%sim%box_max(1:2) <= cfg%sim%box_min(1:2))) then
      error stop 'surface_current_model requires reference_area_m2 or a finite x-y domain area.'
    end if

    electron_idx = find_species_index(cfg, cfg%surface_current%electron_species)
    ion_idx = find_species_index(cfg, cfg%surface_current%ion_species)
    photo_idx = 0
    if (photoelectron_active) photo_idx = find_species_index(cfg, cfg%surface_current%photoelectron_species)
    if (electron_idx == ion_idx .or. &
        (photoelectron_active .and. (electron_idx == photo_idx .or. ion_idx == photo_idx))) then
      error stop 'surface_current_model species references must be distinct.'
    end if
    call validate_automatic_current_species(cfg, electron_idx, 'electron')
    call validate_automatic_current_species(cfg, ion_idx, 'ion')
    if (photoelectron_active) call validate_automatic_current_species(cfg, photo_idx, 'photoelectron')
    if (cfg%particle_species(electron_idx)%q_particle >= 0.0_dp .or. &
        cfg%particle_species(ion_idx)%q_particle <= 0.0_dp) then
      error stop 'Zhao current species must be negative electron and positive ion.'
    end if
    if (abs(abs(cfg%particle_species(electron_idx)%q_particle) - qe) > 1.0e-6_dp*qe .or. &
        abs(abs(cfg%particle_species(ion_idx)%q_particle) - qe) > 1.0e-6_dp*qe) then
      error stop 'Zhao stationary current model currently requires singly charged electron and ion species.'
    end if
    if (photoelectron_active) then
      if (cfg%particle_species(photo_idx)%q_particle >= 0.0_dp) then
        error stop 'Zhao photoelectron species must have negative charge.'
      end if
      if (abs(abs(cfg%particle_species(photo_idx)%q_particle) - qe) > 1.0e-6_dp*qe) then
        error stop 'Zhao stationary current model currently requires singly charged photoelectrons.'
      end if
      if (abs(cfg%particle_species(photo_idx)%m_particle - cfg%particle_species(electron_idx)%m_particle) > &
          1.0e-6_dp*cfg%particle_species(electron_idx)%m_particle) then
        error stop 'Zhao stationary current model requires matching ambient-electron and photoelectron masses.'
      end if
      if (trim(cfg%particle_species(photo_idx)%source_mode) /= 'photo_raycast' .or. &
          .not. cfg%particle_species(photo_idx)%deposit_opposite_charge_on_emit) then
        error stop 'Zhao photoelectron current requires photo_raycast with opposite-charge emission deposit.'
      end if
    end if
    if (.not. is_z_high_reservoir(cfg%particle_species(electron_idx)) .or. &
        .not. is_z_high_reservoir(cfg%particle_species(ion_idx))) then
      error stop 'Zhao ambient electron and ion species require z-high reservoir inflow.'
    end if
    if (photoelectron_active) then
      if (trim(lower_ascii(cfg%particle_species(photo_idx)%inject_face)) /= 'z_high') then
        error stop 'Zhao photoelectron species requires inject_face="z_high".'
      end if
    end if
    call resolve_particle_boundaries( &
      cfg%sim, cfg%particle_boundary_low, cfg%particle_boundary_high, cfg%particle_species(electron_idx), &
      effective_boundary_low, effective_boundary_high &
      )
    if (effective_boundary_high(3) /= bc_open) then
      error stop 'Zhao ambient-electron kinetic closure requires an open z-high particle boundary.'
    end if
    call resolve_particle_boundaries( &
      cfg%sim, cfg%particle_boundary_low, cfg%particle_boundary_high, cfg%particle_species(ion_idx), &
      effective_boundary_low, effective_boundary_high &
      )
    if (effective_boundary_high(3) /= bc_open) then
      error stop 'Zhao ion kinetic closure requires an open z-high particle boundary.'
    end if
    if (photoelectron_active) then
      call resolve_particle_boundaries( &
        cfg%sim, cfg%particle_boundary_low, cfg%particle_boundary_high, cfg%particle_species(photo_idx), &
        effective_boundary_low, effective_boundary_high &
        )
      if (effective_boundary_high(3) /= bc_open) then
        error stop 'Zhao photoelectron kinetic closure requires an open z-high particle boundary.'
      end if
    end if
    if (-cfg%particle_species(electron_idx)%drift_velocity(3) <= 0.0_dp .or. &
        -cfg%particle_species(ion_idx)%drift_velocity(3) <= 0.0_dp) then
      error stop 'Zhao ambient species require positive inward drift at z-high.'
    end if
    if (species_number_density_m3(cfg%particle_species(ion_idx)) <= 0.0_dp) then
      error stop 'Zhao ion species requires a positive number density.'
    end if
    electron_temperature = species_temperature_k(cfg%particle_species(electron_idx))
    ion_temperature = species_temperature_k(cfg%particle_species(ion_idx))
    if (electron_temperature <= 0.0_dp) then
      error stop 'Zhao electron temperature must be positive.'
    end if
    if (photoelectron_active) then
      if (species_temperature_k(cfg%particle_species(photo_idx)) <= 0.0_dp) then
        error stop 'Zhao photoelectron temperature must be positive.'
      end if
    end if
    if (ion_temperature > 0.1_dp*electron_temperature) then
      error stop 'Zhao stationary current model requires cold ions with T_i <= 0.1 T_e.'
    end if
  end subroutine validate_surface_current_model_config

  integer function find_species_index(cfg, species_key) result(species_idx)
    type(app_config), intent(in) :: cfg
    character(len=*), intent(in) :: species_key
    integer :: idx

    species_idx = 0
    do idx = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(idx)%enabled) cycle
      if (trim(cfg%particle_species(idx)%species_key) /= trim(species_key)) cycle
      species_idx = idx
      return
    end do
    error stop 'surface_current_model references an unknown or disabled species: '//trim(species_key)
  end function find_species_index

  subroutine validate_automatic_current_species(cfg, species_idx, role)
    type(app_config), intent(in) :: cfg
    integer, intent(in) :: species_idx
    character(len=*), intent(in) :: role

    if (trim(cfg%particle_species(species_idx)%surface_charge_closure) /= 'fixed_current') then
      error stop 'surface_current_model '//trim(role)//' species requires surface_charge_closure="fixed_current".'
    end if
    if (cfg%particle_species(species_idx)%has_target_absorbed_current_a .or. &
        cfg%particle_species(species_idx)%has_target_emission_current_a) then
      error stop 'surface_current_model species cannot also specify manual target currents.'
    end if
  end subroutine validate_automatic_current_species

  logical function is_z_high_reservoir(spec) result(enabled)
    type(particle_species_spec), intent(in) :: spec

    enabled = spec%boundary_inflow_high(3) == particle_inflow_reservoir .or. &
              (trim(spec%source_mode) == 'reservoir_face' .and. trim(lower_ascii(spec%inject_face)) == 'z_high')
  end function is_z_high_reservoir

  subroutine validate_particle_boundary_override(action, topology_action, context)
    integer(i32), intent(in) :: action, topology_action
    character(len=*), intent(in) :: context

    select case (action)
    case (particle_bc_inherit, bc_open, bc_reflect, bc_redistributed_reflect)
      continue
    case default
      error stop trim(context)//' must be inherit, open, reflect, or redistributed_reflect.'
    end select
    if (topology_action == bc_periodic .and. action /= particle_bc_inherit) then
      error stop trim(context)//' cannot override a periodic domain face.'
    end if
  end subroutine validate_particle_boundary_override

  subroutine validate_particle_boundary_inflow(inflow, effective_topology_action, effective_particle_action, context)
    integer(i32), intent(in) :: inflow, effective_topology_action, effective_particle_action
    character(len=*), intent(in) :: context

    select case (inflow)
    case (particle_inflow_none)
      return
    case (particle_inflow_reservoir)
      continue
    case default
      error stop trim(context)//' must be none or reservoir.'
    end select
    if (effective_topology_action == bc_periodic) then
      error stop trim(context)//' cannot inject through a periodic domain face.'
    end if
    if (effective_particle_action /= bc_open) then
      error stop trim(context)//' requires the effective particle action to be open.'
    end if
  end subroutine validate_particle_boundary_inflow

end submodule bem_app_config_parser_finalize

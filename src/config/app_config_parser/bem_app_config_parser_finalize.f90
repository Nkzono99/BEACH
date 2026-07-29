!> 読み込み済み TOML 設定の正規化・派生値確定・検証を実装する submodule。
submodule(bem_app_config_parser) bem_app_config_parser_finalize
  implicit none
contains

  module procedure finalize_loaded_config
  integer :: i, j, axis
  integer(i32) :: per_batch_particles, physics_status, negative_photo_raycast_count
  integer(i32) :: implicit_photo_raycast_index
  integer(i32) :: n_periodic_axes
  logical :: has_dynamic_source_species
  character(len=64) :: generated_species_key
  character(len=256) :: physics_message

  call normalize_high_level_config(cfg, authoring)
  call normalize_legacy_physics_config( &
    cfg%sim, cfg%field, cfg%periodic2, cfg%panel, cfg%outer_plasma, cfg%coupling &
    )
  call apply_physics_authoring(cfg, authoring)
  call lower_external_boundary_authoring(cfg, authoring)

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
  if (cfg%sim%field_periodic_image_layers < 0_i32) then
    error stop 'sim.field_periodic_image_layers must be >= 0.'
  end if
  if (trim(cfg%sim%field_periodic_far_correction) == 'cached_kneq0' .and. &
      cfg%sim%field_periodic_image_layers < 1_i32) then
    error stop 'cached_kneq0 requires sim.field_periodic_image_layers >= 1.'
  end if
  if (.not. ieee_is_finite(cfg%sim%field_periodic_ewald_alpha) .or. cfg%sim%field_periodic_ewald_alpha < 0.0d0) then
    error stop 'sim.field_periodic_ewald_alpha must be finite and >= 0.'
  end if
  if (cfg%sim%field_periodic_ewald_layers < 0_i32) then
    error stop 'sim.field_periodic_ewald_layers must be >= 0.'
  end if
  if (len_trim(cfg%sim%field_periodic_cache_dir) == 0) then
    error stop 'sim.field_periodic_cache_dir must not be empty.'
  end if
  if (.not. ieee_is_finite(cfg%sim%field_periodic_generation_tolerance) .or. &
      cfg%sim%field_periodic_generation_tolerance <= 0.0_dp) then
    error stop 'sim.field_periodic_generation_tolerance must be finite and > 0.'
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
  cfg%sim%multiple_box_events_policy = lower_ascii(trim(cfg%sim%multiple_box_events_policy))
  select case (trim(cfg%sim%multiple_box_events_policy))
  case ('abort', 'soft_discard')
    continue
  case default
    error stop 'sim.multiple_box_events_policy must be "abort" or "soft_discard".'
  end select
  if (cfg%sim%multiple_box_events_soft_discard_count_limit < 1_i32) then
    error stop 'sim.multiple_box_events_soft_discard_count_limit must be >= 1.'
  end if
  if (.not. ieee_is_finite(cfg%sim%multiple_box_events_soft_discard_abs_charge_limit) .or. &
      cfg%sim%multiple_box_events_soft_discard_abs_charge_limit <= 0.0_dp) then
    error stop 'sim.multiple_box_events_soft_discard_abs_charge_limit must be finite and > 0.'
  end if
  if (cfg%sim%injection_face_phi_grid_n < 1_i32) then
    error stop 'sim.injection_face_phi_grid_n must be >= 1.'
  end if
  if (cfg%sim%raycast_max_bounce < 1_i32) then
    error stop 'sim.raycast_max_bounce must be >= 1.'
  end if
  if (.not. ieee_is_finite(cfg%sim%phi_infty)) then
    error stop 'sim.phi_infty must be finite.'
  end if
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
  if (.not. ieee_is_finite(cfg%sim%sheath_photoelectron_ref_density_cm3) .or. &
      cfg%sim%sheath_photoelectron_ref_density_cm3 < 0.0d0) then
    error stop 'sim.sheath_photoelectron_ref_density_cm3 must be finite and >= 0.'
  end if
  call resolve_batch_duration(cfg)
  per_batch_particles = 0_i32
  has_dynamic_source_species = .false.
  negative_photo_raycast_count = 0_i32
  implicit_photo_raycast_index = 0_i32
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
      if (cfg%particle_species(i)%q_particle < 0.0_dp) then
        negative_photo_raycast_count = negative_photo_raycast_count + 1_i32
        implicit_photo_raycast_index = i
      end if
      call validate_photo_raycast_species(cfg, i)
    case default
      error stop 'Unknown particles.species.source_mode.'
    end select
  end do

  if (per_batch_particles <= 0_i32 .and. .not. has_dynamic_source_species) then
    error stop 'At least one enabled [[particles.species]] entry must have npcls_per_step > 0.'
  end if
  cfg%n_particles = cfg%sim%batch_count*per_batch_particles
  ! kinetic_1dのsame-batch tracked photoelectronは、公開設定を増やさず
  ! mean k=0だけを陰的に更新するmultirate couplingへlowerする。
  if (negative_photo_raycast_count > 0_i32 .and. &
      trim(lower_ascii(cfg%outer_plasma%model)) == 'kinetic_1d' .and. &
      (trim(lower_ascii(cfg%outer_plasma%kinetic_closure)) == 'ambient_linear_debye' .or. &
       trim(lower_ascii(cfg%outer_plasma%kinetic_closure)) == 'zhao_charge_driven') .and. &
      trim(lower_ascii(cfg%coupling%particle_transfer_mode)) == 'electrostatic_1d_instant_return' .and. &
      .not. cfg%coupling%outer_queue_enabled) then
    if (negative_photo_raycast_count /= 1_i32) then
      error stop 'implicit_mean requires exactly one enabled negative photo_raycast species.'
    end if
    if (abs(cfg%particle_species(implicit_photo_raycast_index)%normal_drift_speed) > &
        64.0_dp*epsilon(1.0_dp)) then
      error stop 'implicit_mean Maxwellian photoelectron closure requires normal_drift_speed=0.'
    end if
    if (trim(lower_ascii(cfg%outer_plasma%kinetic_closure)) == 'zhao_charge_driven' .and. &
        trim(lower_ascii(cfg%coupling%steady_start_mode)) /= 'zhao_floating') then
      error stop 'strong-photoelectron implicit Zhao coupling requires steady_start_mode="zhao_floating".'
    end if
    cfg%coupling%update_mode = 'implicit_mean'
  end if
  call validate_active_physics_config( &
    cfg%sim, cfg%field, cfg%periodic2, cfg%panel, cfg%outer_plasma, cfg%coupling, physics_status, physics_message &
    )
  if (physics_status /= physics_config_ok) error stop trim(physics_message)
  if (trim(lower_ascii(cfg%coupling%particle_transfer_mode)) == 'electrostatic_1d_instant_return') then
    do i = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(i)%enabled) cycle
      if (trim(lower_ascii(cfg%particle_species(i)%source_mode)) /= 'photo_raycast') cycle
      if (.not. cfg%particle_species(i)%deposit_opposite_charge_on_emit) then
        error stop 'tracked outer transfer requires photo_raycast deposit_opposite_charge_on_emit=true.'
      end if
    end do
  end if
  end procedure finalize_loaded_config

end submodule bem_app_config_parser_finalize

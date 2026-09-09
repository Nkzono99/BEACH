!> 表面電流と matching-plane のモデル・粒子種・境界条件の preflight。
submodule(bem_app_config_parser) bem_app_config_parser_preflight_surface
  use bem_constants, only: qe
  use bem_config_helpers, only: resolve_particle_boundaries, species_number_density_m3, species_temperature_k
  use bem_app_config_types, only: particle_inflow_none, particle_inflow_reservoir
  implicit none
contains

  module procedure is_automatic_current_species

  selected = .false.
  if (trim(lower_ascii(cfg%surface_current%model)) /= 'zhao_stationary') return
  selected = trim(cfg%particle_species(species_idx)%species_key) == trim(cfg%surface_current%electron_species) .or. &
             trim(cfg%particle_species(species_idx)%species_key) == trim(cfg%surface_current%ion_species) .or. &
             trim(cfg%particle_species(species_idx)%species_key) == trim(cfg%surface_current%photoelectron_species)
  end procedure is_automatic_current_species
  module procedure validate_surface_current_model_config
  integer :: electron_idx, ion_idx, photo_idx
  integer(i32) :: effective_boundary_low(3), effective_boundary_high(3)
  real(dp) :: electron_temperature, ion_temperature
  logical :: photoelectron_active

  select case (trim(lower_ascii(cfg%surface_current%model)))
  case ('none')
    if (has_surface_current_model_settings(cfg)) then
      error stop 'surface_current_model.model="none" cannot use Zhao or matching-plane settings.'
    end if
    return
  case ('zhao_stationary')
    continue
  case ('matching_plane_quasistatic')
    call validate_matching_plane_config(cfg, periodic2_split_explicit)
    return
  case default
    error stop 'surface_current_model.model must be "none", "zhao_stationary", or "matching_plane_quasistatic".'
  end select
  if (cfg%surface_current%has_response_backend .or. &
      cfg%surface_current%has_response_table_path .or. &
      cfg%surface_current%has_zhao_root_selection .or. &
      cfg%surface_current%has_implicit_zero_mode .or. &
      cfg%surface_current%has_coupling_rtol .or. &
      cfg%surface_current%has_coupling_atol .or. &
      cfg%surface_current%has_coupling_max_iterations .or. &
      cfg%surface_current%has_coupling_relaxation .or. &
      trim(lower_ascii(cfg%surface_current%response_backend)) /= 'table' .or. &
      len_trim(cfg%surface_current%response_table_path) > 0 .or. &
      cfg%surface_current%implicit_zero_mode .or. &
      cfg%surface_current%coupling_rtol /= 1.0e-4_dp .or. &
      any(cfg%surface_current%coupling_atol /= 0.0_dp) .or. &
      cfg%surface_current%coupling_max_iterations /= 20_i32 .or. &
      cfg%surface_current%coupling_relaxation /= 0.5_dp) then
    error stop 'surface_current_model="zhao_stationary" cannot use matching-plane-specific settings.'
  end if
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
  end procedure validate_surface_current_model_config
  logical function has_surface_current_model_settings(cfg) result(has_settings)
    type(app_config), intent(in) :: cfg

    has_settings = cfg%surface_current%has_response_backend .or. &
                   cfg%surface_current%has_zhao_branch .or. &
                   cfg%surface_current%has_zhao_root_selection .or. &
                   cfg%surface_current%has_electron_species .or. &
                   cfg%surface_current%has_ion_species .or. &
                   cfg%surface_current%has_photoelectron_species .or. &
                   cfg%surface_current%has_solar_elevation_deg .or. &
                   cfg%surface_current%has_photoelectron_ref_density_m3 .or. &
                   cfg%surface_current%has_photoelectron_source_scale .or. &
                   cfg%surface_current%has_reference_area_m2 .or. &
                   cfg%surface_current%has_response_table_path .or. &
                   cfg%surface_current%has_implicit_zero_mode .or. &
                   cfg%surface_current%has_coupling_rtol .or. &
                   cfg%surface_current%has_coupling_atol .or. &
                   cfg%surface_current%has_coupling_max_iterations .or. &
                   cfg%surface_current%has_coupling_relaxation .or. &
                   len_trim(cfg%surface_current%electron_species) > 0 .or. &
                   len_trim(cfg%surface_current%ion_species) > 0 .or. &
                   len_trim(cfg%surface_current%photoelectron_species) > 0 .or. &
                   trim(lower_ascii(cfg%surface_current%response_backend)) /= 'table' .or. &
                   len_trim(cfg%surface_current%response_table_path) > 0 .or. &
                   cfg%surface_current%implicit_zero_mode .or. &
                   trim(lower_ascii(cfg%surface_current%zhao_branch)) /= 'auto' .or. &
                   trim(lower_ascii(cfg%surface_current%zhao_root_selection)) /= 'require_unique' .or. &
                   cfg%surface_current%solar_elevation_deg /= 0.0_dp .or. &
                   cfg%surface_current%photoelectron_ref_density_m3 /= 0.0_dp .or. &
                   cfg%surface_current%photoelectron_source_scale /= 1.0_dp .or. &
                   cfg%surface_current%reference_area_m2 /= 0.0_dp .or. &
                   cfg%surface_current%coupling_rtol /= 1.0e-4_dp .or. &
                   any(cfg%surface_current%coupling_atol /= 0.0_dp) .or. &
                   cfg%surface_current%coupling_max_iterations /= 20_i32 .or. &
                   cfg%surface_current%coupling_relaxation /= 0.5_dp
  end function has_surface_current_model_settings

  subroutine validate_matching_plane_config(cfg, periodic2_split_explicit)
    type(app_config), intent(in) :: cfg
    logical, intent(in) :: periodic2_split_explicit
    integer :: electron_idx, ion_idx, photo_idx, species_idx
    logical :: photoelectron_active

    if (cfg%surface_current%has_solar_elevation_deg .or. &
        cfg%surface_current%has_photoelectron_ref_density_m3 .or. &
        cfg%surface_current%has_photoelectron_source_scale .or. &
        cfg%surface_current%has_reference_area_m2 .or. &
        cfg%surface_current%solar_elevation_deg /= 0.0_dp .or. &
        cfg%surface_current%photoelectron_ref_density_m3 /= 0.0_dp .or. &
        cfg%surface_current%photoelectron_source_scale /= 1.0_dp .or. &
        cfg%surface_current%reference_area_m2 /= 0.0_dp) then
      error stop 'matching_plane_quasistatic cannot use stationary-Zhao source settings or reference_area_m2.'
    end if
    select case (trim(lower_ascii(cfg%surface_current%response_backend)))
    case ('table')
      if (cfg%surface_current%has_zhao_branch .or. cfg%surface_current%has_zhao_root_selection .or. &
          trim(lower_ascii(cfg%surface_current%zhao_branch)) /= 'auto' .or. &
          trim(lower_ascii(cfg%surface_current%zhao_root_selection)) /= 'require_unique') then
        error stop 'matching_plane_quasistatic response_backend="table" cannot use Zhao-specific settings.'
      end if
      if (.not. cfg%surface_current%has_response_table_path .or. &
          len_trim(cfg%surface_current%response_table_path) == 0) then
        error stop 'matching_plane_quasistatic response_backend="table" requires response_table_path.'
      end if
    case ('zhao_online')
      if (cfg%surface_current%has_response_table_path .or. &
          len_trim(cfg%surface_current%response_table_path) > 0) then
        error stop 'matching_plane_quasistatic response_backend="zhao_online" cannot use response_table_path.'
      end if
      select case (trim(lower_ascii(cfg%surface_current%zhao_branch)))
      case ('auto', 'a', 'b', 'c')
        continue
      case default
        error stop 'surface_current_model.zhao_branch must be "auto", "a", "b", or "c".'
      end select
      select case (trim(lower_ascii(cfg%surface_current%zhao_root_selection)))
      case ('require_unique', 'minimum_energy', 'continuation')
        continue
      case default
        error stop 'surface_current_model.zhao_root_selection must be "require_unique", "minimum_energy", '// &
          'or "continuation".'
      end select
      if (trim(lower_ascii(cfg%surface_current%zhao_root_selection)) == 'continuation') then
        if (trim(lower_ascii(cfg%surface_current%zhao_branch)) /= 'a' .or. &
            .not. cfg%surface_current%implicit_zero_mode) then
          error stop 'surface_current_model.zhao_root_selection="continuation" requires '// &
            'response_backend="zhao_online", zhao_branch="a", and implicit_zero_mode=true.'
        end if
      end if
    case default
      error stop 'surface_current_model.response_backend must be "table" or "zhao_online".'
    end select
    if (cfg%surface_current%implicit_zero_mode) then
      if (trim(lower_ascii(cfg%periodic2%lower_boundary_model)) /= 'e_bottom_zero') then
        error stop 'surface_current_model.implicit_zero_mode requires periodic2.lower_boundary_model="e_bottom_zero".'
      end if
    end if
    if (.not. ieee_is_finite(cfg%surface_current%coupling_rtol) .or. &
        cfg%surface_current%coupling_rtol <= 0.0_dp .or. cfg%surface_current%coupling_rtol > 1.0_dp) then
      error stop 'surface_current_model.coupling_rtol must be finite and in (0, 1].'
    end if
    if (any(.not. ieee_is_finite(cfg%surface_current%coupling_atol)) .or. &
        any(cfg%surface_current%coupling_atol < 0.0_dp)) then
      error stop 'surface_current_model.coupling_atol entries must be finite and >= 0.'
    end if
    if (cfg%surface_current%coupling_max_iterations < 1_i32) then
      error stop 'surface_current_model.coupling_max_iterations must be >= 1.'
    end if
    if (.not. ieee_is_finite(cfg%surface_current%coupling_relaxation) .or. &
        cfg%surface_current%coupling_relaxation <= 0.0_dp .or. &
        cfg%surface_current%coupling_relaxation > 1.0_dp) then
      error stop 'surface_current_model.coupling_relaxation must be finite and in (0, 1].'
    end if
    if (trim(lower_ascii(cfg%sim%field_bc_mode)) /= 'periodic2' .or. .not. cfg%sim%use_box) then
      error stop 'matching_plane_quasistatic requires a periodic2 [domain] box.'
    end if
    if (any(cfg%sim%bc_low(1:2) /= bc_periodic) .or. any(cfg%sim%bc_high(1:2) /= bc_periodic) .or. &
        cfg%sim%bc_low(3) /= bc_open .or. cfg%sim%bc_high(3) /= bc_open) then
      error stop 'matching_plane_quasistatic requires x/y periodic and z-open box topology.'
    end if
    if (.not. periodic2_split_explicit) then
      error stop 'matching_plane_quasistatic requires an explicit split-zero-mode [periodic2] table.'
    end if
    select case (trim(lower_ascii(cfg%periodic2%nonzero_mode_backend)))
    case ('cached_kneq0', 'panel_spectral_reference')
      continue
    case default
      error stop 'matching_plane_quasistatic requires a split periodic2 nonzero-mode backend.'
    end select
    if (trim(lower_ascii(cfg%periodic2%zero_mode_policy)) /= 'exclude_k0') then
      error stop 'matching_plane_quasistatic requires periodic2.zero_mode_policy="exclude_k0".'
    end if
    select case (trim(lower_ascii(cfg%periodic2%lower_boundary_model)))
    case ('e_bottom_zero', 'symmetric_vacuum')
      continue
    case default
      error stop 'matching_plane_quasistatic requires a supported periodic2 lower boundary model.'
    end select
    if (any(cfg%sim%e0 /= 0.0_dp) .or. any(cfg%sim%b0 /= 0.0_dp)) then
      error stop 'matching_plane_quasistatic requires sim.e0=sim.b0=[0,0,0].'
    end if
    if (trim(lower_ascii(cfg%sim%reservoir_potential_model)) /= 'none') then
      error stop 'matching_plane_quasistatic cannot use the generic reservoir potential model.'
    end if
    if (trim(lower_ascii(cfg%sim%open_boundary_model)) /= 'escape') then
      error stop 'matching_plane_quasistatic requires particle_boundary.ordinary_open_model="escape".'
    end if
    if (.not. cfg%surface_current%has_electron_species .or. &
        len_trim(cfg%surface_current%electron_species) == 0 .or. &
        .not. cfg%surface_current%has_ion_species .or. &
        len_trim(cfg%surface_current%ion_species) == 0) then
      error stop 'matching_plane_quasistatic requires electron and ion species roles.'
    end if
    if (cfg%surface_current%has_photoelectron_species .and. &
        len_trim(cfg%surface_current%photoelectron_species) == 0) then
      error stop 'matching_plane_quasistatic photoelectron_species must be a non-empty string when provided.'
    end if
    photoelectron_active = cfg%surface_current%has_photoelectron_species

    electron_idx = find_species_index(cfg, cfg%surface_current%electron_species)
    ion_idx = find_species_index(cfg, cfg%surface_current%ion_species)
    photo_idx = 0
    if (photoelectron_active) photo_idx = find_species_index(cfg, cfg%surface_current%photoelectron_species)
    if (electron_idx == ion_idx .or. &
        (photoelectron_active .and. (electron_idx == photo_idx .or. ion_idx == photo_idx))) then
      error stop 'matching_plane_quasistatic species references must be distinct.'
    end if
    call validate_matching_species(cfg, electron_idx, 'electron')
    call validate_matching_species(cfg, ion_idx, 'ion')
    if (photoelectron_active) call validate_matching_species(cfg, photo_idx, 'photoelectron')
    do species_idx = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(species_idx)%enabled) cycle
      if (trim(lower_ascii(cfg%particle_species(species_idx)%surface_charge_closure)) == 'fixed_current' .or. &
          cfg%particle_species(species_idx)%has_target_absorbed_current_a .or. &
          cfg%particle_species(species_idx)%has_target_emission_current_a) then
        error stop 'matching_plane_quasistatic cannot use manual fixed_current targets on any enabled species.'
      end if
    end do
    if (count(cfg%particle_species(1:cfg%n_particle_species)%enabled) /= merge(3, 2, photoelectron_active)) then
      error stop 'matching_plane_quasistatic requires exactly its enabled electron, ion, and optional photoelectron roles.'
    end if
    if (cfg%particle_species(electron_idx)%q_particle >= 0.0_dp .or. &
        cfg%particle_species(ion_idx)%q_particle <= 0.0_dp) then
      error stop 'matching_plane_quasistatic requires negative electron and positive ion species.'
    end if
    if (photoelectron_active) then
      if (cfg%particle_species(photo_idx)%q_particle >= 0.0_dp) then
        error stop 'matching_plane_quasistatic requires a negative photoelectron species.'
      end if
    end if
    call validate_matching_ambient_source(cfg, electron_idx, 'electron')
    call validate_matching_ambient_source(cfg, ion_idx, 'ion')
    if (photoelectron_active) then
      if (trim(lower_ascii(cfg%particle_species(photo_idx)%source_mode)) /= 'photo_raycast' .or. &
          .not. cfg%particle_species(photo_idx)%deposit_opposite_charge_on_emit .or. &
          trim(lower_ascii(cfg%particle_species(photo_idx)%inject_face)) /= 'z_high') then
        error stop 'matching_plane_quasistatic photoelectrons require opposite-deposit photo_raycast from z_high.'
      end if
    end if

    call validate_matching_species_boundaries(cfg, electron_idx, 'electron')
    call validate_matching_species_boundaries(cfg, ion_idx, 'ion')
    if (photoelectron_active) call validate_matching_species_boundaries(cfg, photo_idx, 'photoelectron')
    if (trim(lower_ascii(cfg%surface_current%response_backend)) == 'zhao_online') then
      call validate_matching_plane_zhao_online(cfg, electron_idx, ion_idx, photo_idx, photoelectron_active)
    end if
  end subroutine validate_matching_plane_config

  subroutine validate_matching_plane_zhao_online(cfg, electron_idx, ion_idx, photo_idx, photoelectron_active)
    type(app_config), intent(in) :: cfg
    integer, intent(in) :: electron_idx, ion_idx, photo_idx
    logical, intent(in) :: photoelectron_active
    real(dp) :: electron_mass, photoelectron_mass
    real(dp) :: electron_temperature, ion_temperature, photoelectron_temperature
    real(dp) :: ion_density

    if (any(cfg%surface_current%coupling_atol(3:4) > 0.0_dp)) then
      error stop 'matching_plane_quasistatic zhao_online coupling_atol must be zero on inactive ambient-outward axes.'
    end if

    if (.not. ieee_is_finite(cfg%particle_species(electron_idx)%q_particle) .or. &
        .not. ieee_is_finite(cfg%particle_species(ion_idx)%q_particle) .or. &
        abs(abs(cfg%particle_species(electron_idx)%q_particle) - qe) > 1.0e-6_dp*qe .or. &
        abs(abs(cfg%particle_species(ion_idx)%q_particle) - qe) > 1.0e-6_dp*qe) then
      error stop 'matching_plane_quasistatic zhao_online requires singly charged role species.'
    end if
    if (photoelectron_active) then
      if (.not. ieee_is_finite(cfg%particle_species(photo_idx)%q_particle) .or. &
          abs(abs(cfg%particle_species(photo_idx)%q_particle) - qe) > 1.0e-6_dp*qe) then
        error stop 'matching_plane_quasistatic zhao_online requires singly charged role species.'
      end if
    end if

    electron_mass = cfg%particle_species(electron_idx)%m_particle
    if (.not. ieee_is_finite(electron_mass) .or. electron_mass <= 0.0_dp) then
      error stop 'matching_plane_quasistatic zhao_online requires a positive electron mass.'
    end if
    if (photoelectron_active) then
      photoelectron_mass = cfg%particle_species(photo_idx)%m_particle
      if (.not. ieee_is_finite(photoelectron_mass) .or. photoelectron_mass <= 0.0_dp .or. &
          abs(photoelectron_mass - electron_mass) > 1.0e-6_dp*electron_mass) then
        error stop 'matching_plane_quasistatic zhao_online requires matching ambient-electron and photoelectron masses.'
      end if
    end if

    electron_temperature = species_temperature_k(cfg%particle_species(electron_idx))
    ion_temperature = species_temperature_k(cfg%particle_species(ion_idx))
    if (.not. ieee_is_finite(electron_temperature) .or. electron_temperature <= 0.0_dp) then
      error stop 'matching_plane_quasistatic zhao_online requires a positive electron temperature.'
    end if
    if (photoelectron_active) then
      photoelectron_temperature = species_temperature_k(cfg%particle_species(photo_idx))
      if (.not. ieee_is_finite(photoelectron_temperature) .or. photoelectron_temperature <= 0.0_dp) then
        error stop 'matching_plane_quasistatic zhao_online requires a positive photoelectron temperature.'
      end if
    end if
    if (.not. ieee_is_finite(ion_temperature) .or. ion_temperature < 0.0_dp .or. &
        ion_temperature > 0.1_dp*electron_temperature) then
      error stop 'matching_plane_quasistatic zhao_online requires cold ions with T_i <= 0.1 T_e.'
    end if

    if (.not. ieee_is_finite(cfg%particle_species(electron_idx)%drift_velocity(3)) .or. &
        .not. ieee_is_finite(cfg%particle_species(ion_idx)%drift_velocity(3)) .or. &
        cfg%particle_species(electron_idx)%drift_velocity(3) >= 0.0_dp .or. &
        cfg%particle_species(ion_idx)%drift_velocity(3) >= 0.0_dp) then
      error stop 'matching_plane_quasistatic zhao_online requires positive ambient inward drift at z-high.'
    end if

    ion_density = species_number_density_m3(cfg%particle_species(ion_idx))
    if (.not. ieee_is_finite(ion_density) .or. ion_density <= 0.0_dp) then
      error stop 'matching_plane_quasistatic zhao_online requires a positive ion number density.'
    end if
  end subroutine validate_matching_plane_zhao_online

  subroutine validate_matching_species(cfg, species_idx, role)
    type(app_config), intent(in) :: cfg
    integer, intent(in) :: species_idx
    character(len=*), intent(in) :: role

    if (trim(lower_ascii(cfg%particle_species(species_idx)%surface_charge_closure)) /= 'explicit') then
      call stop_config_error( &
        'matching_plane_quasistatic '//trim(role)//' species requires surface_charge_closure="explicit".' &
        )
    end if
  end subroutine validate_matching_species

  subroutine validate_matching_ambient_source(cfg, species_idx, role)
    type(app_config), intent(in) :: cfg
    integer, intent(in) :: species_idx
    character(len=*), intent(in) :: role

    if (trim(lower_ascii(cfg%particle_species(species_idx)%source_mode)) /= 'volume_seed' .or. &
        cfg%particle_species(species_idx)%npcls_per_step /= 0_i32 .or. &
        any(cfg%particle_species(species_idx)%boundary_inflow_low /= particle_inflow_none) .or. &
        any(cfg%particle_species(species_idx)%boundary_inflow_high(1:2) /= particle_inflow_none) .or. &
        cfg%particle_species(species_idx)%boundary_inflow_high(3) /= particle_inflow_reservoir) then
      call stop_config_error( &
        'matching_plane_quasistatic '//trim(role)// &
        ' species requires volume_seed, npcls_per_step=0, and only z-high boundary_inflow="reservoir".' &
        )
    end if
  end subroutine validate_matching_ambient_source

  subroutine validate_matching_species_boundaries(cfg, species_idx, role)
    type(app_config), intent(in) :: cfg
    integer, intent(in) :: species_idx
    character(len=*), intent(in) :: role
    integer(i32) :: effective_boundary_low(3), effective_boundary_high(3)

    call resolve_particle_boundaries( &
      cfg%sim, cfg%particle_boundary_low, cfg%particle_boundary_high, cfg%particle_species(species_idx), &
      effective_boundary_low, effective_boundary_high &
      )
    if (.not. all(effective_boundary_low == [bc_periodic, bc_periodic, bc_open]) .or. &
        .not. all(effective_boundary_high == [bc_periodic, bc_periodic, bc_open])) then
      call stop_config_error( &
        'matching_plane_quasistatic '//trim(role)// &
        ' species requires x/y periodic and z-low/z-high open particle boundaries.' &
        )
    end if
  end subroutine validate_matching_species_boundaries

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
    call stop_config_error('surface_current_model references an unknown or disabled species: '//trim(species_key))
  end function find_species_index

  subroutine validate_automatic_current_species(cfg, species_idx, role)
    type(app_config), intent(in) :: cfg
    integer, intent(in) :: species_idx
    character(len=*), intent(in) :: role

    if (trim(cfg%particle_species(species_idx)%surface_charge_closure) /= 'fixed_current') then
      call stop_config_error( &
        'surface_current_model '//trim(role)//' species requires surface_charge_closure="fixed_current".' &
        )
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

end submodule bem_app_config_parser_preflight_surface

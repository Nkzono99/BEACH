!> lint 済みの表面電流設定について、種の参照と物理モデルの成立条件を検査する。
submodule(bem_app_config_parser) bem_app_config_parser_preflight_surface
  use bem_constants, only: qe
  use bem_config_helpers, only: resolve_particle_boundaries, species_number_density_m3, species_temperature_k
  use bem_app_config_types, only: particle_inflow_reservoir
  implicit none
contains

  module procedure validate_surface_current_model_config
  integer :: electron_idx, ion_idx, photo_idx
  integer(i32) :: effective_boundary_low(3), effective_boundary_high(3)
  real(dp) :: electron_temperature, ion_temperature
  logical :: photoelectron_active

  if (cfg%surface_current%outflow_refresh_batches < 0_i32) then
    error stop 'surface_current_model.outflow_refresh_batches must be >= 0.'
  end if
  select case (trim(lower_ascii(cfg%surface_current%model)))
  case ('none')
    if (cfg%surface_current%outflow_refresh_batches /= 0_i32) then
      error stop 'surface_current_model.outflow_refresh_batches requires model="zhao_stationary".'
    end if
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
  photoelectron_active = cfg%surface_current%photoelectron_source_scale > 0.0_dp
  if (.not. photoelectron_active .and. &
      trim(lower_ascii(cfg%surface_current%zhao_branch)) /= 'auto' .and. &
      trim(lower_ascii(cfg%surface_current%zhao_branch)) /= 'c') then
    error stop 'photoelectron_source_scale=0 requires surface_current_model.zhao_branch="auto" or "c".'
  end if
  if (.not. cfg%surface_current%has_reference_area_m2 .and. &
      (.not. cfg%sim%use_box .or. any(cfg%sim%box_max(1:2) <= cfg%sim%box_min(1:2)))) then
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
  if (-cfg%particle_species(electron_idx)%drift_velocity(3) < 0.0_dp .or. &
      -cfg%particle_species(ion_idx)%drift_velocity(3) <= 0.0_dp) then
    error stop 'Zhao ambient species require nonnegative electron and positive ion inward drift at z-high.'
  end if
  if (species_number_density_m3(cfg%particle_species(ion_idx)) <= 0.0_dp) then
    error stop 'Zhao ion species requires a positive number density.'
  end if
  ! 電子とイオンのnumber densityは同じ無限遠の太陽風密度を表す。注入する電子の密度は外部根が決める。
  if (abs(species_number_density_m3(cfg%particle_species(electron_idx)) - &
          species_number_density_m3(cfg%particle_species(ion_idx))) > &
      1.0e-9_dp*species_number_density_m3(cfg%particle_species(ion_idx))) then
    error stop 'Zhao ambient electron and ion species must share the solar-wind number density; '// &
      'the injected electron reservoir density is derived from the Zhao root.'
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
  if (cfg%surface_current%outflow_refresh_batches > 0_i32) then
    call validate_zhao_outflow_refresh_config(cfg, photoelectron_active, periodic2_split_explicit)
  end if
  end procedure validate_surface_current_model_config

  !> 観測PE流出による外部根の更新は、z-high面の平均電位を外部の壁電位へ固定できる周期セルに限る。
  subroutine validate_zhao_outflow_refresh_config(cfg, photoelectron_active, periodic2_split_explicit)
    type(app_config), intent(in) :: cfg
    logical, intent(in) :: photoelectron_active, periodic2_split_explicit

    if (.not. photoelectron_active) then
      error stop 'surface_current_model.outflow_refresh_batches requires photoelectron_source_scale > 0.'
    end if
    if (trim(lower_ascii(cfg%sim%field_bc_mode)) /= 'periodic2' .or. .not. cfg%sim%use_box) then
      error stop 'surface_current_model.outflow_refresh_batches requires a periodic2 [domain] box.'
    end if
    if (any(cfg%sim%bc_low(1:2) /= bc_periodic) .or. any(cfg%sim%bc_high(1:2) /= bc_periodic)) then
      error stop 'surface_current_model.outflow_refresh_batches requires x/y periodic axes.'
    end if
    if (.not. periodic2_split_explicit) then
      error stop 'surface_current_model.outflow_refresh_batches requires an explicit split-zero-mode [periodic2] table.'
    end if
    select case (trim(lower_ascii(cfg%periodic2%nonzero_mode_backend)))
    case ('cached_kneq0', 'panel_spectral_reference')
      continue
    case default
      error stop 'surface_current_model.outflow_refresh_batches requires a split periodic2 nonzero-mode backend.'
    end select
    if (trim(lower_ascii(cfg%periodic2%zero_mode_policy)) /= 'exclude_k0') then
      error stop 'surface_current_model.outflow_refresh_batches requires periodic2.zero_mode_policy="exclude_k0".'
    end if
  end subroutine validate_zhao_outflow_refresh_config

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

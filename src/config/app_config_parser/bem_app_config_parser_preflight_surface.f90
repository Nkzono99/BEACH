!> lint 済みの表面電流設定について、種の参照と物理モデルの成立条件を検査する。
submodule(bem_app_config_parser) bem_app_config_parser_preflight_surface
  use, intrinsic :: iso_fortran_env, only: error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_constants, only: qe, k_boltzmann, pi
  use bem_config_helpers, only: resolve_particle_boundaries, species_number_density_m3, species_temperature_k
  use bem_app_config_types, only: particle_inflow_reservoir
  implicit none
contains

  module procedure validate_surface_current_model_config
  integer :: electron_idx, ion_idx, photo_idx, secondary_idx
  integer(i32) :: effective_boundary_low(3), effective_boundary_high(3)
  real(dp) :: electron_temperature, ion_temperature, ion_entry_energy_ev, ion_pressure_energy_ev
  logical :: photoelectron_active, secondary_active

  if (cfg%surface_current%outflow_refresh_batches < 0_i32) then
    error stop 'sheath.coupling.outflow_refresh_batches must be >= 0.'
  end if
  if (.not. (cfg%surface_current%zhao_upstream_band_tolerance >= 0.0_dp .and. &
             cfg%surface_current%zhao_upstream_band_tolerance < 1.0_dp)) then
    error stop 'sheath.zhao.upstream_band_tolerance must be >= 0 and < 1.'
  end if
  if (.not. (cfg%surface_current%zhao_upstream_field_floor_v_m >= 0.0_dp) .or. &
      .not. ieee_is_finite(cfg%surface_current%zhao_upstream_field_floor_v_m)) then
    error stop 'sheath.zhao.upstream_field_floor_v_m must be finite and >= 0.'
  end if
  if (.not. (cfg%surface_current%zhao_ion_pressure_factor > 0.0_dp) .or. &
      .not. ieee_is_finite(cfg%surface_current%zhao_ion_pressure_factor)) then
    error stop 'sheath.zhao.ion_pressure_factor must be finite and > 0.'
  end if
  select case (trim(lower_ascii(cfg%surface_current%model)))
  case ('none')
    if (cfg%surface_current%outflow_refresh_batches /= 0_i32) then
      error stop 'sheath.coupling.outflow_refresh_batches requires sheath.closure="zero_current".'
    end if
    if (cfg%surface_current%zhao_upstream_band_tolerance /= 0.0_dp) then
      error stop 'sheath.zhao.upstream_band_tolerance requires sheath.closure="zero_current".'
    end if
    if (cfg%surface_current%zhao_upstream_field_floor_v_m /= 0.0_dp) then
      error stop 'sheath.zhao.upstream_field_floor_v_m requires sheath.closure="zero_current".'
    end if
    if (cfg%surface_current%zhao_ion_pressure_factor /= 1.0_dp) then
      error stop 'sheath.zhao.ion_pressure_factor requires sheath.closure="zero_current".'
    end if
    return
  case ('zhao_stationary')
    continue
  case default
    error stop 'sheath.closure must be "none" or "zero_current".'
  end select
  if (any(cfg%sim%b0 /= 0.0_dp)) then
    error stop 'sheath.closure="zero_current" requires fields.external.magnetic_t=[0,0,0] for its unmagnetized sheath closure.'
  end if
  if (trim(lower_ascii(cfg%sim%reservoir_potential_model)) /= 'none') then
    error stop 'Zhao kinetic inflow cannot be combined with the generic reservoir potential model.'
  end if
  if (.not. ieee_is_finite(cfg%surface_current%photoelectron_source_scale) .or. &
      cfg%surface_current%photoelectron_source_scale < 0.0_dp) then
    error stop 'sheath.photoelectrons.source_scale must be finite and >= 0.'
  end if
  photoelectron_active = cfg%surface_current%photoelectron_source_scale > 0.0_dp
  secondary_active = len_trim(cfg%surface_current%secondary_electron_species) > 0
  if (.not. photoelectron_active .and. &
      trim(lower_ascii(cfg%surface_current%zhao_branch)) /= 'auto' .and. &
      trim(lower_ascii(cfg%surface_current%zhao_branch)) /= 'c') then
    error stop 'sheath.photoelectrons.source_scale=0 requires sheath.zhao.branch="auto" or "c".'
  end if
  if (photoelectron_active) then
    if (.not. ieee_is_finite(cfg%surface_current%solar_elevation_deg) .or. &
        cfg%surface_current%solar_elevation_deg <= 0.0_dp .or. cfg%surface_current%solar_elevation_deg > 90.0_dp) then
      error stop 'sheath.photoelectrons.solar_elevation_deg must satisfy 0 < elevation <= 90.'
    end if
    if (.not. ieee_is_finite(cfg%surface_current%photoelectron_ref_density_m3) .or. &
        cfg%surface_current%photoelectron_ref_density_m3 <= 0.0_dp) then
      error stop 'sheath.photoelectrons.ref_density_m3 must be finite and > 0.'
    end if
  else if (len_trim(cfg%surface_current%photoelectron_species) > 0 .or. &
           cfg%surface_current%solar_elevation_deg /= 0.0_dp .or. &
           cfg%surface_current%photoelectron_ref_density_m3 /= 0.0_dp) then
    error stop 'sheath.photoelectrons.source_scale=0 requires omitting all photoelectron-specific Zhao settings.'
  end if
  if (.not. cfg%surface_current%has_reference_area_m2 .and. &
      (.not. cfg%sim%use_box .or. any(cfg%sim%box_max(1:2) <= cfg%sim%box_min(1:2)))) then
    error stop 'sheath requires sheath.reference_area_m2 or a finite x-y domain area.'
  end if

  electron_idx = find_species_index(cfg, cfg%surface_current%electron_species)
  ion_idx = find_species_index(cfg, cfg%surface_current%ion_species)
  photo_idx = 0
  if (photoelectron_active) photo_idx = find_species_index(cfg, cfg%surface_current%photoelectron_species)
  if (electron_idx == ion_idx .or. &
      (photoelectron_active .and. (electron_idx == photo_idx .or. ion_idx == photo_idx))) then
    error stop 'sheath species references must be distinct.'
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
  if (secondary_active) then
    secondary_idx = find_species_index(cfg, cfg%surface_current%secondary_electron_species)
    if (trim(cfg%particle_species(secondary_idx)%source_mode) /= 'secondary_emission' .or. &
        abs(cfg%particle_species(secondary_idx)%q_particle + qe) > 1.0e-6_dp*qe .or. &
        abs(cfg%particle_species(secondary_idx)%m_particle - cfg%particle_species(electron_idx)%m_particle) > &
        1.0e-6_dp*cfg%particle_species(electron_idx)%m_particle) then
      call stop_config_error('Zhao secondary_electron requires a secondary_emission species '// &
                             'with charge -e and ambient electron mass.')
    end if
    call resolve_particle_boundaries(cfg%sim, cfg%particle_boundary_low, cfg%particle_boundary_high, &
                                     cfg%particle_species(secondary_idx), effective_boundary_low, effective_boundary_high)
    if (effective_boundary_high(3) /= bc_open) then
      call stop_config_error('Zhao secondary-electron kinetic closure requires an open z-high particle boundary.')
    end if
  end if
  ! 壁・法線の実行時表現は bem_irradiation_cell_contract。ここは製品の接続面を制限する preflight 規則。
  if (.not. is_z_high_reservoir(cfg%particle_species(electron_idx)) .or. &
      .not. is_z_high_reservoir(cfg%particle_species(ion_idx))) then
    error stop 'Zhao ambient electron and ion species require z-high reservoir inflow.'
  end if
  if (photoelectron_active) then
    if (trim(lower_ascii(cfg%particle_species(photo_idx)%inject_face)) /= 'z_high') then
      error stop species_config_path(photo_idx)//'.source.inject_face must be "z_high" for Zhao photoelectrons.'
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
  ! 外部根のイオン流体 (Sana & Mishra 2026) は、入口が音速より速い場合だけ上流へ接続する。
  ion_entry_energy_ev = 0.5_dp*cfg%particle_species(ion_idx)%m_particle* &
                        cfg%particle_species(ion_idx)%drift_velocity(3)**2/qe
  ion_pressure_energy_ev = cfg%surface_current%zhao_ion_pressure_factor*ion_temperature*k_boltzmann/qe
  if (.not. ieee_is_finite(ion_temperature) .or. ion_temperature < 0.0_dp) then
    error stop 'Zhao ion temperature must be finite and nonnegative.'
  end if
  if (.not. (2.0_dp*ion_entry_energy_ev > ion_pressure_energy_ev)) then
    error stop 'Zhao ion entry must be supersonic: m_i u_z^2/e must exceed sheath.zhao.ion_pressure_factor*T_i.'
  end if
  if (photoelectron_active) call warn_photoelectron_source_mismatch(cfg, photo_idx)
  call validate_zhao_wall_gauge_config(cfg, periodic2_split_explicit)
  if (cfg%surface_current%outflow_refresh_batches > 0_i32 .and. .not. (photoelectron_active .or. secondary_active)) then
    error stop 'sheath.coupling.outflow_refresh_batches requires a surface-emitted electron group.'
  end if
  end procedure validate_surface_current_model_config

  !> 表面放出の生の束と照射幾何が、外部シースへ渡す放出源 $e\,s\,n_{ref}\sin\alpha\,v_{th}/(2\sqrt\pi)$ と
  !! 一致しないときに警告する。固定電流の倍率が差を吸収するので停止はしない。
  !! `emit_current_density_a_m2` はレイ垂直面基準なので、表面放出の面密度は ray 仰角の sin を掛けた値になる。
  subroutine warn_photoelectron_source_mismatch(cfg, photo_idx)
    type(app_config), intent(in) :: cfg
    integer, intent(in) :: photo_idx
    real(dp), parameter :: current_tolerance = 0.05_dp, elevation_tolerance_deg = 1.0_dp
    real(dp) :: ray(3), ray_norm, ray_sin, ray_elevation_deg, thermal_speed, outer_current, surface_current

    ray = cfg%particle_species(photo_idx)%ray_direction
    ray_norm = sqrt(sum(ray*ray))
    if (.not. ieee_is_finite(ray_norm) .or. ray_norm <= 0.0_dp) then
      ! 省略時は inject_face の内向き法線（垂直入射）が validate 段で入る。
      ray_sin = 1.0_dp
    else
      ray_sin = abs(ray(3))/ray_norm
    end if
    ray_elevation_deg = asin(min(1.0_dp, ray_sin))*180.0_dp/pi
    thermal_speed = sqrt(2.0_dp*k_boltzmann*species_temperature_k(cfg%particle_species(photo_idx))/ &
                         cfg%particle_species(photo_idx)%m_particle)
    outer_current = qe*cfg%surface_current%photoelectron_source_scale*cfg%surface_current%photoelectron_ref_density_m3* &
                    sin(cfg%surface_current%solar_elevation_deg*pi/180.0_dp)*thermal_speed/(2.0_dp*sqrt(pi))
    surface_current = cfg%particle_species(photo_idx)%emit_current_density_a_m2*ray_sin
    if (.not. all(ieee_is_finite([outer_current, surface_current])) .or. outer_current <= 0.0_dp) return
    if (abs(surface_current - outer_current) > current_tolerance*outer_current) then
      write (error_unit, '(a,es12.4,a,es12.4,a)') &
        'WARNING: photo_raycast surface emission ', surface_current, &
        ' A/m^2 ('//species_config_path(photo_idx)// &
        '.source.emit_current_density_a_m2 x sin of the ray elevation) differs from the outer Zhao photoelectron source ', &
        outer_current, ' A/m^2 (e s n_ref sin(alpha) v_th/(2 sqrt(pi))); the fixed_current emission scale absorbs '// &
        'the difference, so set '//species_config_path(photo_idx)// &
        '.source.emit_current_density_a_m2 to the normal-incidence value and the ray to the solar elevation.'
      flush (error_unit)
    end if
    if (abs(ray_elevation_deg - cfg%surface_current%solar_elevation_deg) > elevation_tolerance_deg) then
      write (error_unit, '(a,f6.1,a,f6.1,a)') &
        'WARNING: photo_raycast ray elevation ', ray_elevation_deg, &
        ' deg differs from sheath.photoelectrons.solar_elevation_deg ', cfg%surface_current%solar_elevation_deg, &
        ' deg; the illumination inside the cell and the outer photoelectron source assume different sun directions.'
      flush (error_unit)
    end if
  end subroutine warn_photoelectron_source_mismatch

  !> 外部根は上端面の水平平均電位が壁電位 phi0 であることを前提に、流入の写像と外部障壁の判定を組む。
  !! 壁と周期軸の実行時参照は bem_irradiation_cell_contract に集約する。
  !! この電位を固定できるのは、零モードを分けて持つ周期セルだけである。
  subroutine validate_zhao_wall_gauge_config(cfg, periodic2_split_explicit)
    type(app_config), intent(in) :: cfg
    logical, intent(in) :: periodic2_split_explicit
    character(len=*), parameter :: prefix = 'sheath.closure="zero_current" '
    character(len=*), parameter :: reason = &
                                   ' The z-high mean potential is tied to the outer wall potential phi0 '// &
                                   'through the periodic zero mode.'

    if (trim(lower_ascii(cfg%sim%field_bc_mode)) /= 'periodic2' .or. .not. cfg%sim%use_box) then
      error stop prefix//'requires a periodic2 [domain] box.'//reason
    end if
    if (any(cfg%sim%bc_low(1:2) /= bc_periodic) .or. any(cfg%sim%bc_high(1:2) /= bc_periodic)) then
      error stop prefix//'requires x/y periodic axes.'//reason
    end if
    if (.not. periodic2_split_explicit) then
      error stop prefix//'requires an explicit split-zero-mode [fields.periodic] table.'//reason
    end if
    select case (trim(lower_ascii(cfg%periodic2%nonzero_mode_backend)))
    case ('cached_kneq0', 'panel_spectral_reference')
      continue
    case default
      error stop prefix//'requires a split periodic2 nonzero-mode backend '// &
        '("cached_kneq0" or "panel_spectral_reference").'//reason
    end select
    if (trim(lower_ascii(cfg%periodic2%zero_mode_policy)) /= 'exclude_k0') then
      error stop prefix//'requires the split zero mode selected by fields.periodic.backend.'//reason
    end if
  end subroutine validate_zhao_wall_gauge_config

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
    call stop_config_error('sheath references an unknown or disabled species: '//trim(species_key))
  end function find_species_index

  subroutine validate_automatic_current_species(cfg, species_idx, role)
    type(app_config), intent(in) :: cfg
    integer, intent(in) :: species_idx
    character(len=*), intent(in) :: role

    if (trim(cfg%particle_species(species_idx)%surface_charge_closure) /= 'fixed_current') then
      call stop_config_error( &
        'sheath.species.'//trim(role)//' requires '//species_config_path(species_idx)//'.charging.closure="fixed_current".' &
        )
    end if
    if (cfg%particle_species(species_idx)%has_target_absorbed_current_a .or. &
        cfg%particle_species(species_idx)%has_target_emission_current_a) then
      error stop 'sheath species cannot also specify manual target currents.'
    end if
  end subroutine validate_automatic_current_species

  logical function is_z_high_reservoir(spec) result(enabled)
    type(particle_species_spec), intent(in) :: spec

    enabled = spec%boundary_inflow_high(3) == particle_inflow_reservoir
  end function is_z_high_reservoir

end submodule bem_app_config_parser_preflight_surface

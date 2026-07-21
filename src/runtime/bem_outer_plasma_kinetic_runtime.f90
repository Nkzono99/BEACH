module bem_outer_plasma_kinetic_runtime
  use bem_kinds, only: dp, i32
  use bem_constants, only: qe, k_boltzmann, eps0, pi
  use bem_app_config_types, only: app_config, particle_species_spec
  use bem_outer_plasma_types, only: outer_plasma_ok, outer_plasma_not_applicable, outer_plasma_invalid
  use bem_outer_plasma_kinetic, only: kinetic_outer_plasma_options_type
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  real(dp), parameter :: electron_mass_si = 9.1093837139e-31_dp

  public :: resolve_kinetic_outer_options

contains

  subroutine resolve_kinetic_outer_options(app, interface_field, options, status, message)
    type(app_config), intent(in) :: app
    real(dp), intent(in) :: interface_field
    type(kinetic_outer_plasma_options_type), intent(out) :: options
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    integer(i32) :: species, electron_index, ion_index, photoelectron_index
    integer(i32) :: electron_count, ion_count, photoelectron_count
    logical :: use_zhao

    options = kinetic_outer_plasma_options_type()
    options%kinetic_closure = trim(lower_ascii(app%outer_plasma%kinetic_closure))
    options%zhao_branch = trim(lower_ascii(app%outer_plasma%zhao_branch))
    options%zhao_alpha_deg = app%sim%sheath_alpha_deg
    options%photoelectron_reference_density = 1.0e6_dp*app%sim%sheath_photoelectron_ref_density_cm3
    use_zhao = trim(options%kinetic_closure) == 'zhao_charge_driven'
    status = outer_plasma_invalid
    message = ''
    if (use_zhao .and. &
        (trim(lower_ascii(app%sim%sheath_electron_drift_mode)) /= 'normal' .or. &
         trim(lower_ascii(app%sim%sheath_ion_drift_mode)) /= 'normal')) then
      status = outer_plasma_not_applicable
      message = 'Zhao charge-driven closure currently requires normal electron and ion drift modes'
      return
    end if
    if (trim(lower_ascii(app%outer_plasma%model)) /= 'kinetic_1d' .or. &
        app%outer_plasma%debye_length <= 0.0_dp .or. app%outer_plasma%thermal_voltage <= 0.0_dp) then
      message = 'kinetic_1d requires positive debye_length and thermal_voltage'
      return
    end if
    electron_index = 0_i32
    ion_index = 0_i32
    photoelectron_index = 0_i32
    electron_count = 0_i32
    ion_count = 0_i32
    photoelectron_count = 0_i32
    do species = 1_i32, app%n_particle_species
      if (.not. app%particle_species(species)%enabled) cycle
      if ((trim(lower_ascii(app%outer_plasma%photoelectron_density_model)) == 'kinetic_mean' .or. use_zhao) .and. &
          trim(lower_ascii(app%particle_species(species)%source_mode)) == 'photo_raycast' .and. &
          app%particle_species(species)%q_particle < 0.0_dp) then
        photoelectron_count = photoelectron_count + 1_i32
        if (photoelectron_index == 0_i32) photoelectron_index = species
        cycle
      end if
      if (trim(lower_ascii(app%particle_species(species)%source_mode)) /= 'reservoir_face') cycle
      if (trim(lower_ascii(app%particle_species(species)%inject_face)) /= 'z_high') cycle
      if (app%particle_species(species)%q_particle < 0.0_dp) then
        electron_count = electron_count + 1_i32
        if (electron_index == 0_i32) electron_index = species
      end if
      if (app%particle_species(species)%q_particle > 0.0_dp) then
        ion_count = ion_count + 1_i32
        if (ion_index == 0_i32) ion_index = species
      end if
    end do
    if (use_zhao .and. &
        (electron_count /= 1_i32 .or. ion_count /= 1_i32 .or. photoelectron_count /= 1_i32)) then
      status = outer_plasma_not_applicable
      message = 'Zhao charge-driven closure requires exactly one ambient electron, ion, and photoelectron species'
      return
    end if
    if (electron_index == 0_i32 .or. ion_index == 0_i32) then
      status = outer_plasma_not_applicable
      message = 'kinetic_1d requires negative and positive z_high reservoir_face species'
      return
    end if
    if ((trim(lower_ascii(app%outer_plasma%photoelectron_density_model)) == 'kinetic_mean' .or. use_zhao) .and. &
        photoelectron_index == 0_i32) then
      status = outer_plasma_not_applicable
      message = 'selected kinetic closure requires a negative photo_raycast species'
      return
    end if

    options%grid_points = 128_i32
    options%domain_length = 10.0_dp*app%outer_plasma%debye_length
    options%grid_stretch = 2.0_dp
    options%tail_length = app%outer_plasma%debye_length
    options%interface_field = interface_field
    options%max_iterations = 40_i32
    options%residual_tolerance = 1.0e-8_dp
    call map_electron(app%particle_species(electron_index), options, status, message)
    if (status /= outer_plasma_ok) return
    if (use_zhao) then
      call map_ion( &
        app%particle_species(ion_index), options, status, message, &
        inward_drift=resolve_inward_drift(app%particle_species(ion_index), app%sim%sheath_ion_drift_mode) &
        )
    else
      call map_ion(app%particle_species(ion_index), options, status, message)
    end if
    if (status /= outer_plasma_ok) return
    if (use_zhao) then
      if (abs(options%ion_charge - qe) > 0.1_dp*qe) then
        status = outer_plasma_not_applicable
        message = 'Zhao charge-driven closure supports singly charged physical ions'
        return
      end if
      if (options%ion_temperature_j > 0.1_dp*options%electron_temperature_j) then
        status = outer_plasma_not_applicable
        message = 'Zhao charge-driven cold-ion closure requires ion temperature at most 10 percent of electron temperature'
        return
      end if
      options%electron_drift_infinity = resolve_inward_drift( &
                                        app%particle_species(electron_index), app%sim%sheath_electron_drift_mode &
                                        )
      options%ion_drift_infinity = resolve_inward_drift( &
                                   app%particle_species(ion_index), app%sim%sheath_ion_drift_mode &
                                   )
      if (options%electron_drift_infinity <= 0.0_dp .or. options%ion_drift_infinity <= 0.0_dp .or. &
          options%photoelectron_reference_density <= 0.0_dp) then
        status = outer_plasma_not_applicable
        message = 'Zhao charge-driven closure requires positive inward drifts and photoelectron reference density'
        return
      end if
    end if
    if (photoelectron_index > 0_i32) then
      if (use_zhao) then
        call map_zhao_photoelectron(app%particle_species(photoelectron_index), options, status, message)
        if (status == outer_plasma_ok) then
          options%tail_length = sqrt( &
                                eps0*options%photoelectron_temperature_j/ &
                                (options%photoelectron_reference_density*options%photoelectron_charge**2) &
                                )
          options%domain_length = 10.0_dp*options%tail_length
        end if
      else
        call map_photoelectron(app%particle_species(photoelectron_index), options, status, message)
      end if
    end if
  end subroutine resolve_kinetic_outer_options

  subroutine map_electron(species, options, status, message)
    type(particle_species_spec), intent(in) :: species
    type(kinetic_outer_plasma_options_type), intent(inout) :: options
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = outer_plasma_invalid
    message = ''
    options%electron_charge = species%q_particle
    options%electron_mass = species%m_particle
    options%electron_density_infinity = species_density_m3(species)
    options%electron_temperature_j = species_temperature_j(species)
    if (options%electron_charge >= 0.0_dp .or. options%electron_mass <= 0.0_dp .or. &
        options%electron_density_infinity <= 0.0_dp .or. options%electron_temperature_j <= 0.0_dp) then
      message = 'kinetic electron reservoir requires positive density, mass, and temperature'
      return
    end if
    if (abs(options%electron_mass - electron_mass_si) > 0.5_dp*electron_mass_si .or. &
        abs(abs(options%electron_charge) - qe) > 0.1_dp*qe) then
      status = outer_plasma_not_applicable
      message = 'Phase 7 ambient electron closure supports singly charged physical electrons'
      return
    end if
    status = outer_plasma_ok
  end subroutine map_electron

  subroutine map_ion(species, options, status, message, inward_drift)
    type(particle_species_spec), intent(in) :: species
    type(kinetic_outer_plasma_options_type), intent(inout) :: options
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    real(dp), intent(in), optional :: inward_drift

    status = outer_plasma_invalid
    message = ''
    options%ion_charge = species%q_particle
    options%ion_mass = species%m_particle
    options%ion_density_infinity = species_density_m3(species)
    options%ion_temperature_j = species_temperature_j(species)
    options%ion_gamma = 1.0_dp
    options%ion_drift_infinity = -species%drift_velocity(3)
    if (present(inward_drift)) options%ion_drift_infinity = inward_drift
    if (options%ion_charge <= 0.0_dp .or. options%ion_mass <= 0.0_dp .or. &
        options%ion_density_infinity <= 0.0_dp .or. options%ion_temperature_j < 0.0_dp .or. &
        options%ion_drift_infinity <= 0.0_dp) then
      message = 'kinetic ion reservoir requires inward drift and positive density and mass'
      return
    end if
    status = outer_plasma_ok
  end subroutine map_ion

  subroutine map_photoelectron(species, options, status, message)
    type(particle_species_spec), intent(in) :: species
    type(kinetic_outer_plasma_options_type), intent(inout) :: options
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = outer_plasma_invalid
    message = ''
    options%photoelectron_charge = species%q_particle
    options%photoelectron_mass = species%m_particle
    options%photoelectron_temperature_j = species_temperature_j(species)
    if (species%emit_current_density_a_m2 > 0.0_dp .and. species%q_particle /= 0.0_dp) then
      options%photoelectron_emission_flux = species%emit_current_density_a_m2/abs(species%q_particle)
    end if
    if (options%photoelectron_charge >= 0.0_dp .or. options%photoelectron_mass <= 0.0_dp .or. &
        options%photoelectron_temperature_j <= 0.0_dp .or. options%photoelectron_emission_flux <= 0.0_dp) then
      message = 'kinetic photoelectron source requires negative charge, temperature, and emitted current density'
      return
    end if
    status = outer_plasma_ok
  end subroutine map_photoelectron

  subroutine map_zhao_photoelectron(species, options, status, message)
    type(particle_species_spec), intent(in) :: species
    type(kinetic_outer_plasma_options_type), intent(inout) :: options
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp) :: expected_emit_current_density, photoelectron_thermal_speed

    status = outer_plasma_invalid
    message = ''
    options%photoelectron_charge = species%q_particle
    options%photoelectron_mass = species%m_particle
    options%photoelectron_temperature_j = species_temperature_j(species)
    if (options%photoelectron_charge >= 0.0_dp .or. options%photoelectron_mass <= 0.0_dp .or. &
        options%photoelectron_temperature_j <= 0.0_dp) then
      message = 'Zhao charge-driven closure requires a negative physical photoelectron species with positive temperature'
      return
    end if
    if (abs(options%photoelectron_mass - electron_mass_si) > 0.5_dp*electron_mass_si .or. &
        abs(abs(options%photoelectron_charge) - qe) > 0.1_dp*qe) then
      status = outer_plasma_not_applicable
      message = 'Zhao charge-driven closure supports singly charged physical photoelectrons'
      return
    end if
    if (abs(species%normal_drift_speed) > 64.0_dp*epsilon(1.0_dp)) then
      status = outer_plasma_not_applicable
      message = 'Zhao charge-driven closure requires zero photoelectron normal_drift_speed'
      return
    end if
    photoelectron_thermal_speed = sqrt( &
                                  2.0_dp*options%photoelectron_temperature_j/options%photoelectron_mass &
                                  )
    expected_emit_current_density = abs(options%photoelectron_charge)* &
                                    options%photoelectron_reference_density*sin(options%zhao_alpha_deg*pi/180.0_dp)* &
                                    photoelectron_thermal_speed/(2.0_dp*sqrt(pi))
    if (expected_emit_current_density <= 0.0_dp .or. species%emit_current_density_a_m2 <= 0.0_dp .or. &
        abs(species%emit_current_density_a_m2 - expected_emit_current_density) > &
        0.01_dp*expected_emit_current_density) then
      status = outer_plasma_not_applicable
      message = 'Zhao photo_raycast current must match n_ref, alpha, and photoelectron temperature within 1 percent'
      return
    end if
    status = outer_plasma_ok
  end subroutine map_zhao_photoelectron

  pure real(dp) function resolve_inward_drift(species, drift_mode) result(speed)
    type(particle_species_spec), intent(in) :: species
    character(len=*), intent(in) :: drift_mode

    select case (trim(lower_ascii(drift_mode)))
    case ('normal')
      speed = -species%drift_velocity(3)
    case ('full')
      speed = sqrt(sum(species%drift_velocity*species%drift_velocity))
    case default
      speed = 0.0_dp
    end select
  end function resolve_inward_drift

  pure real(dp) function species_density_m3(species) result(density)
    type(particle_species_spec), intent(in) :: species

    if (species%has_number_density_m3) then
      density = species%number_density_m3
    else if (species%has_number_density_cm3) then
      density = 1.0e6_dp*species%number_density_cm3
    else
      density = 0.0_dp
    end if
  end function species_density_m3

  pure real(dp) function species_temperature_j(species) result(temperature)
    type(particle_species_spec), intent(in) :: species

    if (species%has_temperature_ev) then
      temperature = qe*species%temperature_ev
    else
      temperature = k_boltzmann*species%temperature_k
    end if
  end function species_temperature_j

end module bem_outer_plasma_kinetic_runtime

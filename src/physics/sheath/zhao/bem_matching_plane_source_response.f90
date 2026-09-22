!> SI adapter for the source-connected density model; no legacy fallback.
submodule(bem_matching_plane_zhao) bem_matching_plane_source_response
  use bem_source_kinetic_sheath, only: solve_source_kinetic, source_kinetic_root, &
                                       source_kinetic_ok, source_kinetic_invalid, source_kinetic_excluded
  implicit none
contains
  module procedure evaluate_source_response
  type(source_kinetic_root) :: root
  real(dp) :: temperature, velocity, flux_scale, field_scale, mach
  integer :: i, selected, count

  temperature = diagnostics%effective_photoelectron_temperature_ev
  velocity = sqrt(qe*self%electron_temperature_ev/self%electron_mass_kg)
  flux_scale = self%ion_density_m3*velocity
  field_scale = sqrt(self%ion_density_m3*qe*self%electron_temperature_ev/eps0)
  mach = self%ion_drift_mps/sqrt(qe*self%electron_temperature_ev/self%ion_mass_kg)
  call solve_source_kinetic( &
    mach, temperature/self%electron_temperature_ev, input(2)/flux_scale, input(1)/(eps0*field_scale), &
    self%ion_mass_kg/self%electron_mass_kg, diagnostics%kinetic &
    )
  output = 0
  message = diagnostics%kinetic%message
  select case (diagnostics%kinetic%status)
  case (source_kinetic_invalid)
    status = matching_plane_zhao_invalid_argument
    return
  case (source_kinetic_excluded)
    status = matching_plane_zhao_no_physical_solution
    return
  case (source_kinetic_ok)
    continue
  case default
    status = matching_plane_zhao_numerical_failure
    return
  end select
  count = 0
  selected = 0
  do i = 1, size(diagnostics%kinetic%roots)
    if (self%branch_model /= 'auto') then
      if (trim(lower_ascii(diagnostics%kinetic%roots(i)%branch)) /= self%branch_model) cycle
    end if
    count = count + 1
    selected = i
  end do
  if (count == 0) then
    status = matching_plane_zhao_numerical_failure
    message = 'No root of the requested source_kinetic Type was detected; finite search remains unresolved.'
    return
  else if (count > 1) then
    status = matching_plane_zhao_ambiguous_solution
    message = 'Multiple source_kinetic roots detected; all candidates remain in the atlas diagnostics.'
    return
  end if
  root = diagnostics%kinetic%roots(selected)
  output(1) = root%phi_h*self%electron_temperature_ev
  output(2) = root%electron_flux*flux_scale
  output(3) = root%ion_flux*flux_scale
  if (root%branch == 'A' .or. root%branch == 'N') then
    output(4) = root%phi_min*self%electron_temperature_ev
    output(6) = output(4)
  end if
  diagnostics%branch = root%branch
  diagnostics%ambient_electron_density_m3 = root%amplitude*self%ion_density_m3
  diagnostics%residual_norm = max(abs(root%neutrality_residual), abs(root%outer_residual), &
                                  abs(root%field_squared_residual))
  diagnostics%minimum_field_squared_hat = root%minimum_field_squared
  status = matching_plane_zhao_ok
  end procedure
end submodule

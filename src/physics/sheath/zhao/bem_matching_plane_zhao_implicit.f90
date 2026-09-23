!> Type B の上流中性条件から Ne を消去し、電位で BE endpoint を直接探索する。
!! D(phi_H) が折り返す場合も、外側の D scan で狭い区間を飛び越さない。
submodule(bem_matching_plane_zhao) bem_matching_plane_zhao_implicit
  use bem_sheath_model_core, only: integrate_zhao_rho
  implicit none

  integer, parameter :: thermal_grid_intervals = 256
  integer, parameter :: bin_grid_intervals = 4
  integer, parameter :: scalar_max_iterations = 96

contains

  module procedure solve_matching_implicit_endpoint

  type(zhao_params_type) :: params
  type(zhao_matching_root_type), allocatable :: roots(:)
  type(zhao_matching_root_type) :: trial_root, selected_root
  real(dp), allocatable :: grid(:)
  real(dp) :: input(5), tpe, npe, ion_limit, phi_max, source_max, field_scale, displacement_scale
  real(dp) :: residual_tolerance, current_scale, phi_left, phi_right, f_left, f_right, trial_d, trial_response(6)
  real(dp) :: root_phi, root_residual, distance, nearest_distance, second_distance
  real(dp) :: old_value, left_edge, right_edge
  integer :: grid_count, capacity, bin_count, bin, quarter, index, previous, root_count, selected_index
  integer :: iterations
  integer(i32) :: profile_status
  logical :: left_valid, right_valid, root_valid, duplicate, have_seed, saw_profile_failure
  character(len=512) :: profile_message

  handled = self%branch_model == 'b' .or. &
            (self%branch_model == 'auto' .and. self%electron_drift_mps > 0.0_dp)
  displacement = 0.0_dp
  response = 0.0_dp
  candidate = matching_plane_zhao_root_seed_type()
  status = matching_plane_zhao_ok
  message = ''
  if (.not. handled) return
  status = matching_plane_zhao_invalid_argument
  if (.not. all(ieee_is_finite([displacement_before, duration, electron_charge, ion_charge, photoelectron_charge])) &
      .or. duration <= 0.0_dp) then
    message = 'matching-plane implicit endpoint inputs must be finite with a positive duration.'
    return
  end if
  input = [displacement_before, feedback]
  call prepare_matching_zhao_query(self, input, params, tpe, npe, status, message)
  if (status /= matching_plane_zhao_ok) return
  have_seed = self%root_selection == 'continuation' .and. seed%valid
  if (have_seed) then
    if (.not. all(ieee_is_finite([seed%phi0_v, seed%phi_m_v, seed%ambient_electron_density_m3])) .or. &
        seed%ambient_electron_density_m3 <= 0.0_dp) then
      status = matching_plane_zhao_invalid_argument
      message = 'matching-plane implicit endpoint continuation seed is invalid.'
      return
    end if
  end if

  field_scale = params%t_phe_ev/params%lambda_d_phe_ref_m
  displacement_scale = eps0*field_scale
  current_scale = abs(electron_charge)*params%n_swi_inf_m3*params%v_swe_th_mps + &
                  abs(ion_charge)*params%n_swi_inf_m3*params%v_d_ion_mps
  if (photoelectron_active) current_scale = current_scale + abs(photoelectron_charge)*feedback(1)
  ! Resolve BE cancellation using the larger of field-scale and floating-point
  ! current-roundoff tolerances; this does not change the prescribed duration.
  residual_tolerance = max(sqrt(epsilon(1.0_dp))*displacement_scale, &
                           128.0_dp*epsilon(1.0_dp)*duration*current_scale)
  ion_limit = 0.5_dp*params%t_swe_ev*params%mach**2
  source_max = 32.0_dp*params%t_phe_ev
  bin_count = 0
  if (allocated(params%pe_spectrum%flux)) then
    bin_count = size(params%pe_spectrum%flux)
    source_max = params%pe_spectrum%edge(bin_count)
  end if
  phi_max = max(source_max, 8.0_dp*params%t_swe_ev, &
                params%t_phe_ev*(displacement_before/displacement_scale)**2)
  if (have_seed) phi_max = max(phi_max, 2.0_dp*max(0.0_dp, seed%phi0_v) + params%t_phe_ev)
  phi_max = min(phi_max, ion_limit*(1.0_dp - 64.0_dp*epsilon(1.0_dp)))
  if (.not. ieee_is_finite(phi_max) .or. phi_max <= 0.0_dp) then
    status = matching_plane_zhao_numerical_failure
    message = 'matching-plane implicit endpoint potential search scale is invalid.'
    return
  end if

  capacity = thermal_grid_intervals + 2 + bin_grid_intervals*bin_count
  if (have_seed) capacity = capacity + 1
  allocate (grid(capacity), roots(capacity))
  grid_count = 0
  do index = 0, thermal_grid_intervals
    grid_count = grid_count + 1
    grid(grid_count) = phi_max*(real(index, dp)/real(thermal_grid_intervals, dp))**2
  end do
  do bin = 1, bin_count
    left_edge = params%pe_spectrum%edge(bin - 1)
    right_edge = params%pe_spectrum%edge(bin)
    do quarter = 1, bin_grid_intervals
      old_value = left_edge + (right_edge - left_edge)*real(quarter, dp)/real(bin_grid_intervals, dp)
      if (old_value > phi_max) cycle
      grid_count = grid_count + 1
      grid(grid_count) = old_value
    end do
  end do
  if (have_seed .and. seed%phi0_v > 0.0_dp .and. seed%phi0_v < phi_max) then
    grid_count = grid_count + 1
    grid(grid_count) = seed%phi0_v
  end if
  ! Insertion sort is bounded by the recorded spectrum size and keeps bin edges exact.
  do index = 2, grid_count
    old_value = grid(index)
    previous = index - 1
    do while (previous >= 1)
      if (grid(previous) <= old_value) exit
      grid(previous + 1) = grid(previous)
      previous = previous - 1
    end do
    grid(previous + 1) = old_value
  end do

  root_count = 0
  roots = zhao_matching_root_type()
  saw_profile_failure = .false.
  phi_left = grid(1)
  call evaluate_implicit_b_state(params, phi_left, displacement_before, duration, electron_charge, ion_charge, &
                                 photoelectron_active, photoelectron_charge, trial_root, trial_d, trial_response, &
                                 f_left, left_valid)
  do index = 1, grid_count
    phi_right = grid(index)
    if (index > 1 .and. phi_right <= phi_left) cycle
    call evaluate_implicit_b_state(params, phi_right, displacement_before, duration, electron_charge, ion_charge, &
                                   photoelectron_active, photoelectron_charge, trial_root, trial_d, trial_response, &
                                   f_right, right_valid)
    root_valid = .false.
    iterations = 0
    if (right_valid) then
      if (abs(f_right) <= residual_tolerance) then
        root_phi = phi_right
        root_valid = .true.
      else if (index > 1 .and. left_valid) then
        if ((f_left < 0.0_dp .and. f_right > 0.0_dp) .or. (f_left > 0.0_dp .and. f_right < 0.0_dp)) then
          call bisect_implicit_b_root(params, phi_left, phi_right, f_left, f_right, displacement_before, duration, &
                                      electron_charge, ion_charge, photoelectron_active, photoelectron_charge, &
                                      residual_tolerance, root_phi, iterations, root_valid)
        end if
      end if
    end if
    if (root_valid) then
      call evaluate_implicit_b_state(params, root_phi, displacement_before, duration, electron_charge, ion_charge, &
                                     photoelectron_active, photoelectron_charge, trial_root, trial_d, trial_response, &
                                     root_residual, root_valid)
      if (root_valid) then
        trial_root%residual_norm = abs(root_residual)/max(displacement_scale, duration*current_scale)
        trial_root%nonlinear_iterations = int(iterations, i32)
        call validate_matching_root_profile(params, trial_root, trial_d/displacement_scale, &
                                            profile_status, profile_message)
        if (profile_status == matching_plane_zhao_ok) then
          duplicate = .false.
          do previous = 1, root_count
            if (abs(trial_root%phi0_v - roots(previous)%phi0_v) <= &
                1.0e-6_dp*max(params%t_phe_ev, abs(trial_root%phi0_v), abs(roots(previous)%phi0_v))) then
              duplicate = .true.
              if (trial_root%residual_norm < roots(previous)%residual_norm) roots(previous) = trial_root
              exit
            end if
          end do
          if (.not. duplicate) then
            root_count = root_count + 1
            roots(root_count) = trial_root
          end if
        else if (profile_status == matching_plane_zhao_numerical_failure) then
          saw_profile_failure = .true.
        end if
      end if
    end if
    phi_left = phi_right
    f_left = f_right
    left_valid = right_valid
  end do

  status = matching_plane_zhao_numerical_failure
  if (saw_profile_failure) then
    message = 'a located implicit Zhao-B endpoint profile could not be certified numerically.'
    return
  end if
  if (root_count == 0) then
    message = 'finite potential search found no admissible implicit Zhao-B endpoint.'
    return
  end if
  selected_index = 1
  if (have_seed) then
    nearest_distance = huge(1.0_dp)
    second_distance = huge(1.0_dp)
    do index = 1, root_count
      distance = max(abs(roots(index)%phi0_v - seed%phi0_v)/params%t_phe_ev, &
                     abs(min(seed%phi_m_v, 0.0_dp))/params%t_phe_ev, &
                     abs(log(roots(index)%ambient_electron_density_m3/seed%ambient_electron_density_m3)))
      if (distance < nearest_distance) then
        second_distance = nearest_distance
        nearest_distance = distance
        selected_index = index
      else if (distance < second_distance) then
        second_distance = distance
      end if
    end do
    if (root_count > 1 .and. abs(second_distance - nearest_distance) <= &
        1.0e-6_dp*max(1.0_dp, nearest_distance)) then
      status = matching_plane_zhao_ambiguous_solution
      message = 'implicit Zhao-B continuation found indistinguishable nearest endpoints.'
      return
    end if
    selected_root = roots(selected_index)
  else if (self%root_selection == 'minimum_energy') then
    call select_minimum_energy_root(params, roots, root_count, selected_root, status, message)
    if (status /= matching_plane_zhao_ok) return
  else if (root_count > 1) then
    status = matching_plane_zhao_ambiguous_solution
    message = 'finite potential search found multiple physical implicit Zhao-B endpoints.'
    return
  else
    selected_root = roots(1)
  end if
  status = matching_plane_zhao_numerical_failure
  call evaluate_implicit_b_state(params, selected_root%phi0_v, displacement_before, duration, &
                                 electron_charge, ion_charge, photoelectron_active, photoelectron_charge, &
                                 trial_root, displacement, response, root_residual, root_valid)
  if (.not. root_valid) then
    message = 'selected implicit Zhao-B endpoint could not be reevaluated.'
    return
  end if
  candidate%branch = 'B'
  candidate%valid = .true.
  candidate%phi0_v = selected_root%phi0_v
  candidate%phi_m_v = selected_root%phi0_v
  candidate%ambient_electron_density_m3 = selected_root%ambient_electron_density_m3
  status = matching_plane_zhao_ok
  message = ''
  end procedure solve_matching_implicit_endpoint

  subroutine evaluate_implicit_b_state(params, phi, displacement_before, duration, electron_charge, ion_charge, &
                                       photoelectron_active, photoelectron_charge, root, displacement, response, residual, valid)
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: phi, displacement_before, duration, electron_charge, ion_charge, photoelectron_charge
    logical, intent(in) :: photoelectron_active
    type(zhao_matching_root_type), intent(out) :: root
    real(dp), intent(out) :: displacement, response(6), residual
    logical, intent(out) :: valid
    real(dp) :: ion, free, reflected, photo, captured, density_hat, phi_hat, e2_hat, escape_flux, current
    integer(i32) :: status
    character(len=512) :: message

    root = zhao_matching_root_type()
    root%branch = 'B'
    root%phi0_v = phi
    root%phi_m_v = phi
    displacement = 0.0_dp
    response = 0.0_dp
    residual = 0.0_dp
    valid = .false.
    phi_hat = phi/params%t_phe_ev
    call evaluate_zhao_density_hat(params, 'B', 'monotonic', 0.0_dp, phi_hat, phi_hat, &
                                   1.0_dp, ion, free, reflected, photo, captured)
    if (.not. all(ieee_is_finite([ion, free, reflected, photo, captured]))) return
    if (free + reflected <= 0.0_dp .or. ion <= photo + captured) return
    density_hat = (ion - photo - captured)/(free + reflected)
    root%ambient_electron_density_m3 = density_hat*params%n_phe_ref_m3
    e2_hat = 2.0_dp*integrate_zhao_rho(params, 'B', 'monotonic', phi_hat, 0.0_dp, phi_hat, phi_hat, density_hat)
    if (.not. ieee_is_finite(e2_hat) .or. e2_hat < 0.0_dp) return
    displacement = eps0*(params%t_phe_ev/params%lambda_d_phe_ref_m)*sqrt(e2_hat)
    call compose_matching_response(params, root, response, status, message)
    if (status /= matching_plane_zhao_ok) return
    current = electron_charge*response(2) + ion_charge*response(3)
    if (photoelectron_active) then
      if (allocated(params%pe_spectrum%flux)) then
        escape_flux = params%pe_spectrum%tail_flux(phi)
      else
        escape_flux = params%n_phe0_m3*params%v_phe_th_mps/(2.0_dp*sqrt(pi))*exp(-phi_hat)
      end if
      current = current - photoelectron_charge*escape_flux
    end if
    residual = displacement - displacement_before - duration*current
    valid = all(ieee_is_finite([displacement, residual]))
  end subroutine evaluate_implicit_b_state

  subroutine bisect_implicit_b_root(params, left, right, f_left, f_right, displacement_before, duration, &
                                    electron_charge, ion_charge, photoelectron_active, photoelectron_charge, &
                                    tolerance, phi, iterations, valid)
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: left, right, f_left, f_right, displacement_before, duration
    real(dp), intent(in) :: electron_charge, ion_charge, photoelectron_charge, tolerance
    logical, intent(in) :: photoelectron_active
    real(dp), intent(out) :: phi
    integer, intent(out) :: iterations
    logical, intent(out) :: valid
    type(zhao_matching_root_type) :: root
    real(dp) :: lo, hi, flo, fhi, fmid, displacement, response(6)

    lo = left
    hi = right
    flo = f_left
    fhi = f_right
    valid = .false.
    do iterations = 1, scalar_max_iterations
      phi = 0.5_dp*lo + 0.5_dp*hi
      call evaluate_implicit_b_state(params, phi, displacement_before, duration, electron_charge, ion_charge, &
                                     photoelectron_active, photoelectron_charge, root, displacement, response, fmid, valid)
      ! An invalid midpoint separates domains; never bridge it with a bracket.
      if (.not. valid) return
      if (abs(fmid) <= tolerance) return
      if ((flo < 0.0_dp .and. fmid > 0.0_dp) .or. (flo > 0.0_dp .and. fmid < 0.0_dp)) then
        hi = phi
        fhi = fmid
      else
        lo = phi
        flo = fmid
      end if
      if (hi - lo <= 64.0_dp*epsilon(1.0_dp)*max(params%t_phe_ev, abs(phi))) exit
    end do
    if (abs(flo) <= abs(fhi)) then
      phi = lo
      valid = abs(flo) <= tolerance
    else
      phi = hi
      valid = abs(fhi) <= tolerance
    end if
  end subroutine bisect_implicit_b_root

end submodule bem_matching_plane_zhao_implicit

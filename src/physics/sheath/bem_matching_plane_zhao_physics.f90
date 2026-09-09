!> Zhao query の物理量への変換、残差式、Sagdeev 積分、接続プロファイルと流入応答。
!! 根の探索順・選択方針を持たず、与えられた状態の物理量と成立条件を評価する。
submodule(bem_matching_plane_zhao) bem_matching_plane_zhao_physics
  use bem_sheath_model_core, only: evaluate_zhao_rho_hat, &
                                   zhao_residuals_type_a, zhao_residuals_type_b, zhao_residuals_type_c
  implicit none

  integer, parameter :: rho_quadrature_panels = 256
  integer, parameter :: energy_quadrature_panels = 128
  integer, parameter :: profile_validation_samples = 32
  real(dp), parameter :: profile_negative_tolerance = 1.0e-7_dp
  real(dp), parameter :: profile_endpoint_tolerance = 1.0e-5_dp

contains

  module procedure prepare_matching_zhao_query

  real(dp) :: photoelectron_flux_m2_s, photoelectron_thermal_speed_mps

  params = zhao_params_type()
  photoelectron_temperature_ev = 0.0_dp
  photoelectron_source_density_m3 = 0.0_dp
  status = matching_plane_zhao_invalid_argument
  message = ''
  if (.not. self%initialized) then
    message = 'matching-plane Zhao model is not initialized.'
    return
  end if
  if (.not. all(ieee_is_finite(input))) then
    message = 'matching-plane Zhao query must be finite.'
    return
  end if
  if (input(matching_plane_input_photoelectron_outward_flux) < 0.0_dp .or. &
      input(matching_plane_input_photoelectron_mean_normal_energy) < 0.0_dp .or. &
      input(matching_plane_input_electron_outward_flux) < 0.0_dp .or. &
      input(matching_plane_input_ion_outward_flux) < 0.0_dp) then
    message = 'matching-plane Zhao fluxes and photoelectron energy must be nonnegative.'
    return
  end if

  photoelectron_flux_m2_s = input(matching_plane_input_photoelectron_outward_flux)
  if (photoelectron_flux_m2_s > 0.0_dp) then
    photoelectron_temperature_ev = input(matching_plane_input_photoelectron_mean_normal_energy)
    if (photoelectron_temperature_ev <= 0.0_dp) then
      message = 'positive matching-plane photoelectron flux requires positive mean normal energy.'
      return
    end if
  else
    photoelectron_temperature_ev = self%configured_photoelectron_temperature_ev
  end if

  photoelectron_thermal_speed_mps = sqrt( &
                                    2.0_dp*qe*photoelectron_temperature_ev/self%electron_mass_kg &
                                    )
  if (.not. ieee_is_finite(photoelectron_thermal_speed_mps) .or. &
      photoelectron_thermal_speed_mps <= 0.0_dp) then
    message = 'matching-plane Zhao photoelectron thermal speed is invalid.'
    return
  end if
  photoelectron_source_density_m3 = &
    2.0_dp*sqrt(pi)*photoelectron_flux_m2_s/photoelectron_thermal_speed_mps
  if (.not. ieee_is_finite(photoelectron_source_density_m3) .or. &
      photoelectron_source_density_m3 < 0.0_dp) then
    message = 'matching-plane Zhao photoelectron moment map is invalid.'
    return
  end if

  call prepare_matching_zhao_params( &
    self, photoelectron_temperature_ev, photoelectron_source_density_m3, params, status, message &
    )
  end procedure prepare_matching_zhao_query

  !> `build_zhao_params` の error-stop API を通さず online query を正規化する。
  subroutine prepare_matching_zhao_params( &
    self, photoelectron_temperature_ev, photoelectron_source_density_m3, &
    params, status, message &
    )
    class(matching_plane_zhao_model_type), intent(in) :: self
    real(dp), intent(in) :: photoelectron_temperature_ev, photoelectron_source_density_m3
    type(zhao_params_type), intent(out) :: params
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    params = zhao_params_type()
    status = matching_plane_zhao_numerical_failure
    message = ''
    params%alpha_rad = 0.5_dp*pi
    params%n_swi_inf_m3 = self%ion_density_m3
    params%n_phe_ref_m3 = self%ion_density_m3
    params%n_phe0_m3 = photoelectron_source_density_m3
    params%photoelectron_population_fraction = 1.0_dp
    params%t_swe_ev = self%electron_temperature_ev
    params%t_phe_ev = photoelectron_temperature_ev
    params%v_d_electron_mps = self%electron_drift_mps
    params%v_d_ion_mps = self%ion_drift_mps
    params%m_i_kg = self%ion_mass_kg
    params%m_e_kg = self%electron_mass_kg
    params%v_swe_th_mps = sqrt(2.0_dp*qe*params%t_swe_ev/params%m_e_kg)
    params%v_phe_th_mps = sqrt(2.0_dp*qe*params%t_phe_ev/params%m_e_kg)
    params%cs_mps = sqrt(qe*params%t_swe_ev/params%m_i_kg)
    params%mach = params%v_d_ion_mps/params%cs_mps
    params%u = params%v_d_electron_mps/params%v_swe_th_mps
    params%tau = params%t_swe_ev/params%t_phe_ev
    params%lambda_d_phe_ref_m = sqrt( &
                                eps0*qe*params%t_phe_ev/(params%n_phe_ref_m3*qe*qe) &
                                )
    if (.not. all(ieee_is_finite([ &
                                 params%alpha_rad, params%n_swi_inf_m3, params%n_phe_ref_m3, params%n_phe0_m3, &
                                 params%photoelectron_population_fraction, params%t_swe_ev, params%t_phe_ev, &
                                 params%v_d_electron_mps, params%v_d_ion_mps, params%m_i_kg, params%m_e_kg, &
                                 params%v_swe_th_mps, params%v_phe_th_mps, params%cs_mps, params%mach, &
                                 params%u, params%tau, params%lambda_d_phe_ref_m &
                                 ])) .or. params%n_phe0_m3 < 0.0_dp .or. &
        min(params%v_swe_th_mps, params%v_phe_th_mps, params%cs_mps, &
            params%mach, params%tau, params%lambda_d_phe_ref_m) <= 0.0_dp) then
      params = zhao_params_type()
      message = 'matching-plane Zhao normalization produced invalid parameters.'
      return
    end if
    status = matching_plane_zhao_ok
  end subroutine prepare_matching_zhao_params

  module procedure encode_matching_unknowns

  y = 0.0_dp
  valid = density_m3 > 0.0_dp .and. params%n_phe_ref_m3 > 0.0_dp .and. params%t_phe_ev > 0.0_dp
  if (.not. valid) return
  select case (branch)
  case ('A')
    valid = phi0_v > 0.0_dp .and. phi_m_v < 0.0_dp
    if (.not. valid) return
    y(1) = log(phi0_v/params%t_phe_ev)
    y(2) = log(-phi_m_v/params%t_phe_ev)
    y(3) = log(density_m3/params%n_phe_ref_m3)
  case ('B')
    valid = phi0_v > 0.0_dp
    if (.not. valid) return
    y(1) = log(phi0_v/params%t_phe_ev)
    y(2) = log(density_m3/params%n_phe_ref_m3)
  case ('C')
    valid = phi0_v < 0.0_dp
    if (.not. valid) return
    y(1) = log(-phi0_v/params%t_phe_ev)
    y(2) = log(density_m3/params%n_phe_ref_m3)
  case default
    valid = .false.
  end select
  valid = valid .and. all(ieee_is_finite(y))
  end procedure encode_matching_unknowns

  module procedure decode_matching_unknowns

  phi0_v = 0.0_dp
  phi_m_v = 0.0_dp
  density_m3 = 0.0_dp
  valid = all(ieee_is_finite(y))
  if (.not. valid .or. y(1) < -50.0_dp .or. y(1) > log(200.0_dp)) then
    valid = .false.
    return
  end if
  select case (branch)
  case ('A')
    if (y(2) < -50.0_dp .or. y(2) > log(200.0_dp) .or. &
        y(3) < -30.0_dp .or. y(3) > log(1.0e6_dp)) then
      valid = .false.
      return
    end if
    phi0_v = params%t_phe_ev*exp(y(1))
    phi_m_v = -params%t_phe_ev*exp(y(2))
    density_m3 = params%n_phe_ref_m3*exp(y(3))
  case ('B')
    if (y(2) < -30.0_dp .or. y(2) > log(1.0e6_dp)) then
      valid = .false.
      return
    end if
    phi0_v = params%t_phe_ev*exp(y(1))
    phi_m_v = phi0_v
    density_m3 = params%n_phe_ref_m3*exp(y(2))
  case ('C')
    if (y(2) < -30.0_dp .or. y(2) > log(1.0e6_dp)) then
      valid = .false.
      return
    end if
    phi0_v = -params%t_phe_ev*exp(y(1))
    phi_m_v = phi0_v
    density_m3 = params%n_phe_ref_m3*exp(y(2))
  case default
    valid = .false.
  end select
  valid = valid .and. all(ieee_is_finite([phi0_v, phi_m_v, density_m3]))
  end procedure decode_matching_unknowns

  module procedure evaluate_charge_residual

  real(dp) :: phi0_v, phi_m_v, density_m3, phi0_hat, phi_m_hat, density_hat
  real(dp) :: raw(3), integral, field_squared, field_residual_scale
  real(dp) :: x3(3), x2(2)
  logical :: integral_ok

  residual = 0.0_dp
  call decode_matching_unknowns( &
    params, branch, y, phi0_v, phi_m_v, density_m3, valid &
    )
  if (.not. valid) return
  phi0_hat = phi0_v/params%t_phe_ev
  phi_m_hat = phi_m_v/params%t_phe_ev
  density_hat = density_m3/params%n_phe_ref_m3
  if (.not. ion_accessible(params, max(phi0_hat, 0.0_dp))) then
    valid = .false.
    return
  end if

  select case (branch)
  case ('A')
    x3 = [phi0_v, phi_m_v, density_m3]
    call zhao_residuals_type_a(params, x3, raw)
    call integrate_matching_rho_hat( &
      params, branch, 'lower', phi_m_hat, phi0_hat, phi0_hat, phi_m_hat, &
      density_hat, integral, integral_ok &
      )
    if (.not. integral_ok) then
      valid = .false.
      return
    end if
    field_squared = -2.0_dp*integral
    if (field_squared < -1.0e-10_dp) then
      valid = .false.
      return
    end if
    field_residual_scale = max(1.0_dp, target_field_hat*target_field_hat)
    residual(1) = raw(1)/params%n_phe_ref_m3
    residual(2) = (max(0.0_dp, field_squared) - target_field_hat*target_field_hat)/field_residual_scale
    residual(3) = raw(3)
  case ('B', 'C')
    x2 = [phi0_v, density_m3]
    if (branch == 'B') then
      call zhao_residuals_type_b(params, x2, raw(1:2))
    else
      call zhao_residuals_type_c(params, x2, raw(1:2))
    end if
    call integrate_matching_rho_hat( &
      params, branch, 'monotonic', phi0_hat, 0.0_dp, phi0_hat, phi_m_hat, &
      density_hat, integral, integral_ok &
      )
    if (.not. integral_ok) then
      valid = .false.
      return
    end if
    field_squared = 2.0_dp*integral
    if (field_squared < -1.0e-10_dp) then
      valid = .false.
      return
    end if
    field_residual_scale = max(1.0_dp, target_field_hat*target_field_hat)
    residual(1) = raw(1)/params%n_phe_ref_m3
    residual(2) = (max(0.0_dp, field_squared) - target_field_hat*target_field_hat)/field_residual_scale
  case default
    valid = .false.
    return
  end select
  valid = all(ieee_is_finite(residual))
  end procedure evaluate_charge_residual

  subroutine integrate_matching_rho_hat( &
    params, branch, side, lower_phi_hat, upper_phi_hat, phi0_hat, phi_m_hat, &
    density_hat, integral, success &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    character(len=*), intent(in) :: side
    real(dp), intent(in) :: lower_phi_hat, upper_phi_hat, phi0_hat, phi_m_hat, density_hat
    real(dp), intent(out) :: integral
    logical, intent(out) :: success

    real(dp) :: t, phi_hat, jacobian, rho_hat, summand, weight, h
    integer :: point

    integral = 0.0_dp
    success = .false.
    if (.not. all(ieee_is_finite([ &
                                 lower_phi_hat, upper_phi_hat, phi0_hat, phi_m_hat, density_hat &
                                 ])) .or. density_hat <= 0.0_dp) return
    h = 1.0_dp/real(rho_quadrature_panels, dp)
    do point = 0, rho_quadrature_panels
      t = real(point, dp)*h
      phi_hat = lower_phi_hat + (upper_phi_hat - lower_phi_hat)*sin(0.5_dp*pi*t)**2
      jacobian = (upper_phi_hat - lower_phi_hat)*0.5_dp*pi*sin(pi*t)
      if (.not. ion_accessible(params, phi_hat)) return
      call evaluate_zhao_rho_hat( &
        params, branch, side, phi_hat, phi0_hat, phi_m_hat, density_hat, rho_hat &
        )
      if (.not. ieee_is_finite(rho_hat)) return
      summand = rho_hat*jacobian
      if (point == 0 .or. point == rho_quadrature_panels) then
        weight = 1.0_dp
      else if (mod(point, 2) == 0) then
        weight = 2.0_dp
      else
        weight = 4.0_dp
      end if
      integral = integral + weight*summand
    end do
    integral = integral*h/3.0_dp
    success = ieee_is_finite(integral)
  end subroutine integrate_matching_rho_hat

  module procedure validate_matching_root_profile

  real(dp) :: phi0_hat, phi_m_hat, density_hat, phi_hat, fraction
  real(dp) :: integral, field_squared, interface_field_squared, upper_endpoint_field_squared
  real(dp) :: minimum_field_squared, field_squared_scale
  integer :: point
  logical :: integral_ok

  status = matching_plane_zhao_numerical_failure
  message = ''
  root%minimum_field_squared_hat = huge(1.0_dp)
  phi0_hat = root%phi0_v/params%t_phe_ev
  phi_m_hat = root%phi_m_v/params%t_phe_ev
  density_hat = root%ambient_electron_density_m3/params%n_phe_ref_m3
  field_squared_scale = max(1.0_dp, target_field_hat*target_field_hat)
  if (.not. all(ieee_is_finite([ &
                               phi0_hat, phi_m_hat, density_hat, field_squared_scale &
                               ])) .or. density_hat <= 0.0_dp) then
    message = 'matching-plane Zhao profile normalization is invalid.'
    return
  end if

  minimum_field_squared = huge(1.0_dp)
  interface_field_squared = huge(1.0_dp)
  upper_endpoint_field_squared = 0.0_dp
  select case (root%branch)
  case ('A')
    do point = 0, profile_validation_samples
      fraction = real(point, dp)/real(profile_validation_samples, dp)
      phi_hat = phi_m_hat + fraction*(phi0_hat - phi_m_hat)
      call integrate_matching_rho_hat( &
        params, 'A', 'lower', phi_m_hat, phi_hat, phi0_hat, phi_m_hat, &
        density_hat, integral, integral_ok &
        )
      if (.not. integral_ok) then
        message = 'matching-plane Zhao lower profile integration failed.'
        return
      end if
      field_squared = -2.0_dp*integral
      minimum_field_squared = min(minimum_field_squared, field_squared)
      if (point == profile_validation_samples) interface_field_squared = field_squared
    end do
    do point = 0, profile_validation_samples
      fraction = real(point, dp)/real(profile_validation_samples, dp)
      phi_hat = phi_m_hat + fraction*(0.0_dp - phi_m_hat)
      call integrate_matching_rho_hat( &
        params, 'A', 'upper', phi_m_hat, phi_hat, phi0_hat, phi_m_hat, &
        density_hat, integral, integral_ok &
        )
      if (.not. integral_ok) then
        message = 'matching-plane Zhao upper profile integration failed.'
        return
      end if
      field_squared = -2.0_dp*integral
      minimum_field_squared = min(minimum_field_squared, field_squared)
      if (point == profile_validation_samples) upper_endpoint_field_squared = field_squared
    end do
  case ('B', 'C')
    do point = 0, profile_validation_samples
      fraction = real(point, dp)/real(profile_validation_samples, dp)
      phi_hat = phi0_hat + fraction*(0.0_dp - phi0_hat)
      call integrate_matching_rho_hat( &
        params, root%branch, 'monotonic', phi_hat, 0.0_dp, phi0_hat, phi_m_hat, &
        density_hat, integral, integral_ok &
        )
      if (.not. integral_ok) then
        message = 'matching-plane Zhao monotonic profile integration failed.'
        return
      end if
      field_squared = 2.0_dp*integral
      minimum_field_squared = min(minimum_field_squared, field_squared)
      if (point == 0) interface_field_squared = field_squared
    end do
  case default
    message = 'matching-plane Zhao profile has an unknown branch.'
    return
  end select

  if (.not. all(ieee_is_finite([ &
                               minimum_field_squared, interface_field_squared, upper_endpoint_field_squared &
                               ]))) then
    message = 'matching-plane Zhao profile field is non-finite.'
    return
  end if
  root%minimum_field_squared_hat = minimum_field_squared
  if (minimum_field_squared < -profile_negative_tolerance*field_squared_scale) then
    status = matching_plane_zhao_no_physical_solution
    message = 'matching-plane Zhao profile requires an imaginary electric field.'
    return
  end if
  if (abs(interface_field_squared - target_field_hat*target_field_hat) > &
      profile_endpoint_tolerance*field_squared_scale) then
    message = 'matching-plane Zhao profile does not reproduce the interface field.'
    return
  end if
  if (root%branch == 'A' .and. &
      abs(upper_endpoint_field_squared) > profile_endpoint_tolerance*field_squared_scale) then
    message = 'matching-plane Zhao-A upper profile does not reach zero upstream field.'
    return
  end if
  status = matching_plane_zhao_ok
  end procedure validate_matching_root_profile

  module procedure evaluate_root_potential_energy

  real(dp) :: phi0_hat, phi_m_hat, density_hat, energy_hat, segment_energy_hat
  logical :: success

  status = matching_plane_zhao_numerical_failure
  message = ''
  root%potential_energy_j_m2 = huge(1.0_dp)
  phi0_hat = root%phi0_v/params%t_phe_ev
  phi_m_hat = root%phi_m_v/params%t_phe_ev
  density_hat = root%ambient_electron_density_m3/params%n_phe_ref_m3
  if (.not. all(ieee_is_finite([phi0_hat, phi_m_hat, density_hat])) .or. &
      density_hat <= 0.0_dp .or. params%lambda_d_phe_ref_m <= 0.0_dp) then
    message = 'matching-plane Zhao potential-energy normalization is invalid.'
    return
  end if

  energy_hat = 0.0_dp
  select case (root%branch)
  case ('A')
    call integrate_matching_field_energy_hat( &
      params, root%branch, 'lower', phi_m_hat, phi0_hat, phi0_hat, phi_m_hat, &
      density_hat, segment_energy_hat, success &
      )
    if (.not. success) then
      message = 'matching-plane Zhao-A lower potential-energy integral failed.'
      return
    end if
    energy_hat = energy_hat + segment_energy_hat
    call integrate_matching_field_energy_hat( &
      params, root%branch, 'upper', phi_m_hat, 0.0_dp, phi0_hat, phi_m_hat, &
      density_hat, segment_energy_hat, success &
      )
    if (.not. success) then
      message = 'matching-plane Zhao-A upper potential-energy integral failed.'
      return
    end if
    energy_hat = energy_hat + segment_energy_hat
  case ('B', 'C')
    call integrate_matching_field_energy_hat( &
      params, root%branch, 'monotonic', phi0_hat, 0.0_dp, phi0_hat, phi_m_hat, &
      density_hat, energy_hat, success &
      )
    if (.not. success) then
      message = 'matching-plane Zhao monotonic potential-energy integral failed.'
      return
    end if
  case default
    message = 'matching-plane Zhao potential-energy root has an unknown branch.'
    return
  end select
  root%potential_energy_j_m2 = -0.5_dp*eps0*params%t_phe_ev*params%t_phe_ev* &
                               energy_hat/params%lambda_d_phe_ref_m
  if (.not. ieee_is_finite(root%potential_energy_j_m2) .or. root%potential_energy_j_m2 > 0.0_dp) then
    root%potential_energy_j_m2 = huge(1.0_dp)
    message = 'matching-plane Zhao potential energy is invalid.'
    return
  end if
  status = matching_plane_zhao_ok
  end procedure evaluate_root_potential_energy

  subroutine integrate_matching_field_energy_hat( &
    params, branch, side, start_phi_hat, end_phi_hat, phi0_hat, phi_m_hat, &
    density_hat, energy_hat, success &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    character(len=*), intent(in) :: side
    real(dp), intent(in) :: start_phi_hat, end_phi_hat, phi0_hat, phi_m_hat, density_hat
    real(dp), intent(out) :: energy_hat
    logical, intent(out) :: success

    real(dp) :: t, phi_hat, jacobian, rho_integral, field_squared, summand, weight, h
    real(dp) :: field_squared_scale
    real(dp) :: energy_integrand(0:energy_quadrature_panels)
    integer :: point
    logical :: integral_ok, point_ok(0:energy_quadrature_panels)

    energy_hat = 0.0_dp
    success = .false.
    if (.not. all(ieee_is_finite([ &
                                 start_phi_hat, end_phi_hat, phi0_hat, phi_m_hat, density_hat &
                                 ])) .or. density_hat <= 0.0_dp) return
    h = 1.0_dp/real(energy_quadrature_panels, dp)
    field_squared_scale = max(1.0_dp, phi0_hat*phi0_hat, phi_m_hat*phi_m_hat)
    energy_integrand = 0.0_dp
    point_ok = .false.
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(params,branch,side,start_phi_hat,end_phi_hat,phi0_hat,phi_m_hat,density_hat,h) &
    !$omp shared(field_squared_scale,energy_integrand,point_ok) &
    !$omp private(point,t,phi_hat,jacobian,rho_integral,field_squared,integral_ok)
    do point = 0, energy_quadrature_panels
      t = real(point, dp)*h
      phi_hat = start_phi_hat + (end_phi_hat - start_phi_hat)*sin(0.5_dp*pi*t)**2
      jacobian = (end_phi_hat - start_phi_hat)*0.5_dp*pi*sin(pi*t)
      if (branch == 'A') then
        call integrate_matching_rho_hat( &
          params, branch, side, phi_m_hat, phi_hat, phi0_hat, phi_m_hat, &
          density_hat, rho_integral, integral_ok &
          )
        field_squared = -2.0_dp*rho_integral
      else
        call integrate_matching_rho_hat( &
          params, branch, side, phi_hat, 0.0_dp, phi0_hat, phi_m_hat, &
          density_hat, rho_integral, integral_ok &
          )
        field_squared = 2.0_dp*rho_integral
      end if
      if (.not. integral_ok) cycle
      if (field_squared < -profile_negative_tolerance*field_squared_scale) cycle
      energy_integrand(point) = sqrt(max(0.0_dp, field_squared))*abs(jacobian)
      point_ok(point) = .true.
    end do
    !$omp end parallel do
    if (.not. all(point_ok)) return

    ! Preserve the original Simpson accumulation order for reproducible root ranking.
    do point = 0, energy_quadrature_panels
      summand = energy_integrand(point)
      if (point == 0 .or. point == energy_quadrature_panels) then
        weight = 1.0_dp
      else if (mod(point, 2) == 0) then
        weight = 2.0_dp
      else
        weight = 4.0_dp
      end if
      energy_hat = energy_hat + weight*summand
    end do
    energy_hat = energy_hat*h/3.0_dp
    success = ieee_is_finite(energy_hat) .and. energy_hat >= 0.0_dp
  end subroutine integrate_matching_field_energy_hat

  module procedure compose_matching_response

  real(dp) :: electron_cutoff, electron_term, number_flux_scale
  real(dp) :: electron_access_potential_v, photoelectron_barrier_potential_v

  output = 0.0_dp
  status = matching_plane_zhao_numerical_failure
  message = ''
  select case (root%branch)
  case ('A')
    electron_cutoff = sqrt(max(0.0_dp, -root%phi_m_v/params%t_swe_ev)) - params%u
    electron_access_potential_v = root%phi_m_v
    photoelectron_barrier_potential_v = root%phi_m_v
  case ('B')
    electron_cutoff = -params%u
    electron_access_potential_v = 0.0_dp
    photoelectron_barrier_potential_v = 0.0_dp
  case ('C')
    electron_cutoff = sqrt(max(0.0_dp, -root%phi0_v/params%t_swe_ev)) - params%u
    electron_access_potential_v = 0.0_dp
    photoelectron_barrier_potential_v = 0.0_dp
  case default
    message = 'matching-plane Zhao root has an unknown branch.'
    return
  end select

  electron_term = swe_free_current_term( &
                  params, root%ambient_electron_density_m3, electron_cutoff &
                  )
  number_flux_scale = params%v_phe_th_mps/(2.0_dp*sqrt(pi))
  output(matching_plane_output_matching_potential) = root%phi0_v
  output(matching_plane_output_electron_inward_flux) = number_flux_scale*electron_term
  output(matching_plane_output_ion_inward_flux) = params%n_swi_inf_m3*params%v_d_ion_mps
  output(matching_plane_output_electron_access_potential) = electron_access_potential_v
  output(matching_plane_output_ion_access_potential) = 0.0_dp
  output(matching_plane_output_photoelectron_barrier_potential) = &
    photoelectron_barrier_potential_v
  if (.not. all(ieee_is_finite(output)) .or. &
      any(output(matching_plane_output_electron_inward_flux:matching_plane_output_ion_inward_flux) < 0.0_dp)) then
    output = 0.0_dp
    message = 'matching-plane Zhao response is non-finite or has a negative inward flux.'
    return
  end if
  status = matching_plane_zhao_ok
  end procedure compose_matching_response

  pure logical function ion_accessible(params, phi_hat) result(accessible)
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: phi_hat

    accessible = params%tau > 0.0_dp .and. params%mach > 0.0_dp .and. &
                 1.0_dp - 2.0_dp*phi_hat/(params%tau*params%mach*params%mach) > 0.0_dp
  end function ion_accessible

end submodule bem_matching_plane_zhao_physics

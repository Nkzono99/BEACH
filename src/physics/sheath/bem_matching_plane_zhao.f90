!> Matching-plane 用の stateless charge-driven Zhao 外部シース応答。
!!
!! 応答は整合面直下の D_z/eps0 を外部側 interface field とみなし、
!! 上流の準中性条件と Sagdeev 積分を満たす A/B/C branch を解く。
!! 定常表面電流の零電流条件は課さないため、BEACH 側の帯電過程を消去しない。
!! 平面・半無限の外部問題は z 方向に並進対称なので、絶対高度 H は数値式に
!! 入らず、runtime がこの応答を domain の z-high gauge へ結び付ける。
module bem_matching_plane_zhao
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, pi, qe
  use bem_matching_plane_response, only: &
    matching_plane_response_input_count, matching_plane_response_output_count, &
    matching_plane_input_displacement, matching_plane_input_photoelectron_outward_flux, &
    matching_plane_input_photoelectron_mean_normal_energy, &
    matching_plane_input_electron_outward_flux, matching_plane_input_ion_outward_flux, &
    matching_plane_output_matching_potential, matching_plane_output_electron_inward_flux, &
    matching_plane_output_ion_inward_flux, matching_plane_output_electron_access_potential, &
    matching_plane_output_ion_access_potential, matching_plane_output_photoelectron_barrier_potential
  use bem_sheath_model_core, only: &
    zhao_params_type, evaluate_zhao_rho_hat, &
    zhao_residuals_type_a, zhao_residuals_type_b, zhao_residuals_type_c, &
    swe_free_current_term
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  integer, parameter :: root_max_iterations = 60
  integer, parameter :: root_max_backtracks = 24
  integer, parameter :: rho_quadrature_panels = 256
  integer, parameter :: energy_quadrature_panels = 128
  integer, parameter :: profile_validation_samples = 32
  real(dp), parameter :: root_tolerance = 1.0e-9_dp
  real(dp), parameter :: zero_field_tolerance_hat = 1.0e-12_dp
  real(dp), parameter :: root_cluster_tolerance = 1.0e-6_dp
  real(dp), parameter :: profile_negative_tolerance = 1.0e-7_dp
  real(dp), parameter :: profile_endpoint_tolerance = 1.0e-5_dp
  real(dp), parameter :: energy_tie_tolerance = 1.0e-6_dp

  integer(i32), parameter, public :: matching_plane_zhao_ok = 0_i32
  integer(i32), parameter, public :: matching_plane_zhao_invalid_argument = 1_i32
  integer(i32), parameter, public :: matching_plane_zhao_no_physical_solution = 2_i32
  integer(i32), parameter, public :: matching_plane_zhao_numerical_failure = 3_i32
  integer(i32), parameter, public :: matching_plane_zhao_ambiguous_solution = 4_i32

  type, public :: matching_plane_zhao_diagnostics_type
    character(len=1) :: branch = ' '
    real(dp) :: interface_field_v_m = 0.0_dp
    real(dp) :: effective_photoelectron_temperature_ev = 0.0_dp
    real(dp) :: photoelectron_source_density_m3 = 0.0_dp
    real(dp) :: ambient_electron_density_m3 = 0.0_dp
    real(dp) :: residual_norm = huge(1.0_dp)
    real(dp) :: minimum_field_squared_hat = huge(1.0_dp)
    real(dp) :: potential_energy_j_m2 = huge(1.0_dp)
    integer(i32) :: nonlinear_iterations = 0_i32
  end type matching_plane_zhao_diagnostics_type

  type :: zhao_matching_root_type
    character(len=1) :: branch = ' '
    real(dp) :: phi0_v = 0.0_dp
    real(dp) :: phi_m_v = 0.0_dp
    real(dp) :: ambient_electron_density_m3 = 0.0_dp
    real(dp) :: residual_norm = huge(1.0_dp)
    real(dp) :: minimum_field_squared_hat = huge(1.0_dp)
    real(dp) :: potential_energy_j_m2 = huge(1.0_dp)
    integer(i32) :: nonlinear_iterations = 0_i32
  end type zhao_matching_root_type

  type, public :: matching_plane_zhao_model_type
    private
    logical :: initialized = .false.
    character(len=9) :: branch_model = 'auto'
    character(len=16) :: root_selection = 'require_unique'
    real(dp) :: ion_density_m3 = 0.0_dp
    real(dp) :: electron_temperature_ev = 0.0_dp
    real(dp) :: electron_drift_mps = 0.0_dp
    real(dp) :: ion_drift_mps = 0.0_dp
    real(dp) :: ion_mass_kg = 0.0_dp
    real(dp) :: electron_mass_kg = 0.0_dp
    real(dp) :: configured_photoelectron_temperature_ev = 0.0_dp
  contains
    procedure, public :: initialize => initialize_matching_plane_zhao
    procedure, public :: evaluate => evaluate_matching_plane_zhao
    procedure, public :: get_feedback_scales => get_matching_plane_zhao_feedback_scales
    procedure, public :: is_initialized => matching_plane_zhao_is_initialized
  end type matching_plane_zhao_model_type

contains

  subroutine initialize_matching_plane_zhao( &
    self, branch_model, root_selection, ion_density_m3, electron_temperature_ev, electron_drift_mps, &
    ion_drift_mps, ion_mass_kg, electron_mass_kg, configured_photoelectron_temperature_ev, &
    status, message &
    )
    class(matching_plane_zhao_model_type), intent(inout) :: self
    character(len=*), intent(in) :: branch_model
    character(len=*), intent(in) :: root_selection
    real(dp), intent(in) :: ion_density_m3, electron_temperature_ev, electron_drift_mps
    real(dp), intent(in) :: ion_drift_mps, ion_mass_kg, electron_mass_kg
    real(dp), intent(in) :: configured_photoelectron_temperature_ev
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    character(len=:), allocatable :: normalized_branch, normalized_root_selection

    self%initialized = .false.
    self%branch_model = 'auto'
    self%root_selection = 'require_unique'
    self%ion_density_m3 = 0.0_dp
    self%electron_temperature_ev = 0.0_dp
    self%electron_drift_mps = 0.0_dp
    self%ion_drift_mps = 0.0_dp
    self%ion_mass_kg = 0.0_dp
    self%electron_mass_kg = 0.0_dp
    self%configured_photoelectron_temperature_ev = 0.0_dp
    status = matching_plane_zhao_invalid_argument
    message = ''
    normalized_branch = trim(lower_ascii(branch_model))
    select case (normalized_branch)
    case ('auto', 'zhao_auto')
      self%branch_model = 'auto'
    case ('a', 'zhao_a')
      self%branch_model = 'a'
    case ('b', 'zhao_b')
      self%branch_model = 'b'
    case ('c', 'zhao_c')
      self%branch_model = 'c'
    case default
      message = 'matching-plane Zhao branch must be auto, a, b, or c.'
      return
    end select
    normalized_root_selection = trim(lower_ascii(root_selection))
    select case (normalized_root_selection)
    case ('require_unique', 'minimum_energy')
      self%root_selection = normalized_root_selection
    case default
      message = 'matching-plane Zhao root selection must be require_unique or minimum_energy.'
      return
    end select
    if (.not. all(ieee_is_finite([ &
                                 ion_density_m3, electron_temperature_ev, electron_drift_mps, ion_drift_mps, &
                                 ion_mass_kg, electron_mass_kg, configured_photoelectron_temperature_ev &
                                 ]))) then
      message = 'matching-plane Zhao initialization values must be finite.'
      return
    end if
    if (ion_density_m3 <= 0.0_dp .or. electron_temperature_ev <= 0.0_dp .or. &
        ion_drift_mps <= 0.0_dp .or. ion_mass_kg <= 0.0_dp .or. electron_mass_kg <= 0.0_dp .or. &
        configured_photoelectron_temperature_ev <= 0.0_dp) then
      message = 'matching-plane Zhao densities, temperatures, ion drift, and masses must be positive.'
      return
    end if

    self%ion_density_m3 = ion_density_m3
    self%electron_temperature_ev = electron_temperature_ev
    self%electron_drift_mps = electron_drift_mps
    self%ion_drift_mps = ion_drift_mps
    self%ion_mass_kg = ion_mass_kg
    self%electron_mass_kg = electron_mass_kg
    self%configured_photoelectron_temperature_ev = configured_photoelectron_temperature_ev
    self%initialized = .true.
    status = matching_plane_zhao_ok
  end subroutine initialize_matching_plane_zhao

  subroutine evaluate_matching_plane_zhao(self, input, output, status, message, diagnostics)
    class(matching_plane_zhao_model_type), intent(in) :: self
    real(dp), intent(in) :: input(matching_plane_response_input_count)
    real(dp), intent(out) :: output(matching_plane_response_output_count)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(matching_plane_zhao_diagnostics_type), intent(out), optional :: diagnostics

    type(zhao_params_type) :: params
    type(zhao_matching_root_type) :: root
    type(matching_plane_zhao_diagnostics_type) :: local_diagnostics
    real(dp) :: interface_field_v_m, photoelectron_flux_m2_s
    real(dp) :: photoelectron_temperature_ev, photoelectron_source_density_m3
    real(dp) :: photoelectron_thermal_speed_mps

    output = 0.0_dp
    status = matching_plane_zhao_invalid_argument
    message = ''
    local_diagnostics = matching_plane_zhao_diagnostics_type()
    if (.not. self%initialized) then
      message = 'matching-plane Zhao model is not initialized.'
      call assign_diagnostics(diagnostics, local_diagnostics)
      return
    end if
    if (.not. all(ieee_is_finite(input))) then
      message = 'matching-plane Zhao query must be finite.'
      call assign_diagnostics(diagnostics, local_diagnostics)
      return
    end if
    if (input(matching_plane_input_photoelectron_outward_flux) < 0.0_dp .or. &
        input(matching_plane_input_photoelectron_mean_normal_energy) < 0.0_dp .or. &
        input(matching_plane_input_electron_outward_flux) < 0.0_dp .or. &
        input(matching_plane_input_ion_outward_flux) < 0.0_dp) then
      message = 'matching-plane Zhao fluxes and photoelectron energy must be nonnegative.'
      call assign_diagnostics(diagnostics, local_diagnostics)
      return
    end if

    photoelectron_flux_m2_s = input(matching_plane_input_photoelectron_outward_flux)
    if (photoelectron_flux_m2_s > 0.0_dp) then
      photoelectron_temperature_ev = input(matching_plane_input_photoelectron_mean_normal_energy)
      if (photoelectron_temperature_ev <= 0.0_dp) then
        message = 'positive matching-plane photoelectron flux requires positive mean normal energy.'
        call assign_diagnostics(diagnostics, local_diagnostics)
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
      call assign_diagnostics(diagnostics, local_diagnostics)
      return
    end if
    photoelectron_source_density_m3 = &
      2.0_dp*sqrt(pi)*photoelectron_flux_m2_s/photoelectron_thermal_speed_mps
    if (.not. ieee_is_finite(photoelectron_source_density_m3) .or. &
        photoelectron_source_density_m3 < 0.0_dp) then
      message = 'matching-plane Zhao photoelectron moment map is invalid.'
      call assign_diagnostics(diagnostics, local_diagnostics)
      return
    end if

    call prepare_matching_zhao_params( &
      self, photoelectron_temperature_ev, photoelectron_source_density_m3, &
      params, status, message &
      )
    if (status /= matching_plane_zhao_ok) then
      call assign_diagnostics(diagnostics, local_diagnostics)
      return
    end if

    interface_field_v_m = input(matching_plane_input_displacement)/eps0
    local_diagnostics%interface_field_v_m = interface_field_v_m
    local_diagnostics%effective_photoelectron_temperature_ev = photoelectron_temperature_ev
    local_diagnostics%photoelectron_source_density_m3 = photoelectron_source_density_m3
    call solve_matching_root( &
      trim(self%branch_model), trim(self%root_selection), params, interface_field_v_m, root, status, message &
      )
    if (status /= matching_plane_zhao_ok) then
      call assign_diagnostics(diagnostics, local_diagnostics)
      return
    end if

    local_diagnostics%branch = root%branch
    local_diagnostics%ambient_electron_density_m3 = root%ambient_electron_density_m3
    local_diagnostics%residual_norm = root%residual_norm
    local_diagnostics%minimum_field_squared_hat = root%minimum_field_squared_hat
    local_diagnostics%potential_energy_j_m2 = root%potential_energy_j_m2
    local_diagnostics%nonlinear_iterations = root%nonlinear_iterations
    call compose_matching_response(params, root, output, status, message)
    if (status /= matching_plane_zhao_ok) then
      output = 0.0_dp
      call assign_diagnostics(diagnostics, local_diagnostics)
      return
    end if
    call assign_diagnostics(diagnostics, local_diagnostics)
  end subroutine evaluate_matching_plane_zhao

  !> `build_zhao_params` の error-stop APIを通さず、online queryをfail-closedに正規化する。
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

  subroutine get_matching_plane_zhao_feedback_scales(self, scales, status, message)
    class(matching_plane_zhao_model_type), intent(in) :: self
    real(dp), intent(out) :: scales(4)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp) :: configured_photoelectron_thermal_speed_mps

    scales = 0.0_dp
    status = matching_plane_zhao_invalid_argument
    message = ''
    if (.not. self%initialized) then
      message = 'matching-plane Zhao model is not initialized.'
      return
    end if
    configured_photoelectron_thermal_speed_mps = sqrt( &
                                                 2.0_dp*qe*self%configured_photoelectron_temperature_ev/self%electron_mass_kg &
                                                 )
    scales(1) = self%ion_density_m3*configured_photoelectron_thermal_speed_mps/(2.0_dp*sqrt(pi))
    scales(2) = self%configured_photoelectron_temperature_ev
    ! v1 Zhao moment closure does not reconstruct ambient outward VDFs.  Their
    ! query entries are accepted for the common 5-input ABI but are inactive.
    scales(3:4) = 0.0_dp
    if (.not. all(ieee_is_finite(scales)) .or. any(scales(1:2) <= 0.0_dp)) then
      scales = 0.0_dp
      status = matching_plane_zhao_numerical_failure
      message = 'matching-plane Zhao feedback scales are invalid.'
      return
    end if
    status = matching_plane_zhao_ok
  end subroutine get_matching_plane_zhao_feedback_scales

  pure logical function matching_plane_zhao_is_initialized(self) result(initialized)
    class(matching_plane_zhao_model_type), intent(in) :: self

    initialized = self%initialized
  end function matching_plane_zhao_is_initialized

  subroutine solve_matching_root(model, root_selection, params, interface_field_v_m, root, status, message)
    character(len=*), intent(in) :: model, root_selection
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: interface_field_v_m
    type(zhao_matching_root_type), intent(out) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    character(len=1) :: order(3), candidate
    type(zhao_matching_root_type) :: trial_root, successful_root, successful_roots(3)
    real(dp) :: field_scale, target_field_hat, degenerate_density_m3
    integer :: candidate_count, candidate_index, successful_count
    logical :: saw_numerical_failure, saw_ambiguous_solution

    root = zhao_matching_root_type()
    status = matching_plane_zhao_no_physical_solution
    message = ''
    field_scale = params%t_phe_ev/params%lambda_d_phe_ref_m
    target_field_hat = interface_field_v_m/field_scale
    if (.not. all(ieee_is_finite([field_scale, target_field_hat])) .or. field_scale <= 0.0_dp) then
      status = matching_plane_zhao_numerical_failure
      message = 'matching-plane Zhao field normalization is invalid.'
      return
    end if

    call matching_branch_order(model, target_field_hat, order, candidate_count, status, message)
    if (status /= matching_plane_zhao_ok) return
    if (abs(target_field_hat) <= zero_field_tolerance_hat) then
      if (trim(model) /= 'auto' .and. trim(model) /= 'b') then
        status = matching_plane_zhao_no_physical_solution
        message = 'requested Zhao branch does not contain the zero-field state.'
        return
      end if
      degenerate_density_m3 = ( &
                              2.0_dp*params%n_swi_inf_m3 - params%n_phe0_m3 &
                              )/(1.0_dp + erf(params%u))
      if (.not. ieee_is_finite(degenerate_density_m3) .or. degenerate_density_m3 <= 0.0_dp) then
        status = matching_plane_zhao_no_physical_solution
        message = 'zero-field Zhao-B state has no positive ambient electron density.'
        return
      end if
      root%branch = 'B'
      root%ambient_electron_density_m3 = degenerate_density_m3
      root%residual_norm = 0.0_dp
      root%minimum_field_squared_hat = 0.0_dp
      root%potential_energy_j_m2 = 0.0_dp
      root%nonlinear_iterations = 0_i32
      status = matching_plane_zhao_ok
      message = 'zero-field degenerate Zhao-B state'
      return
    end if

    saw_numerical_failure = .false.
    saw_ambiguous_solution = .false.
    successful_count = 0
    successful_root = zhao_matching_root_type()
    do candidate_index = 1, candidate_count
      candidate = order(candidate_index)
      call solve_one_matching_branch( &
        params, candidate, target_field_hat, root_selection, trial_root, status, message &
        )
      if (status == matching_plane_zhao_ok) then
        successful_count = successful_count + 1
        successful_roots(successful_count) = trial_root
        if (successful_count == 1) successful_root = trial_root
      end if
      if (status == matching_plane_zhao_numerical_failure) saw_numerical_failure = .true.
      if (status == matching_plane_zhao_ambiguous_solution) saw_ambiguous_solution = .true.
    end do
    if (saw_ambiguous_solution) then
      status = matching_plane_zhao_ambiguous_solution
      message = 'matching-plane Zhao branch search found multiple roots within one branch.'
      return
    end if
    if (successful_count == 1 .and. saw_numerical_failure .and. trim(model) == 'auto') then
      status = matching_plane_zhao_numerical_failure
      message = 'matching-plane Zhao auto selection could not certify a unique branch.'
      return
    else if (successful_count == 1) then
      root = successful_root
      status = matching_plane_zhao_ok
      message = ''
      return
    else if (successful_count > 1) then
      if (saw_numerical_failure .and. trim(model) == 'auto' .and. &
          trim(root_selection) == 'minimum_energy') then
        status = matching_plane_zhao_numerical_failure
        message = 'matching-plane Zhao minimum-energy selection could not certify every candidate branch.'
      else if (trim(root_selection) == 'minimum_energy') then
        call select_minimum_energy_root( &
          params, successful_roots, successful_count, root, status, message &
          )
      else
        status = matching_plane_zhao_ambiguous_solution
        message = 'matching-plane Zhao auto selection is ambiguous across multiple physical branches.'
      end if
      return
    end if
    if (saw_numerical_failure) then
      status = matching_plane_zhao_numerical_failure
      message = 'matching-plane Zhao branch search did not converge.'
    else
      status = matching_plane_zhao_no_physical_solution
      message = 'no Zhao branch satisfies the prescribed matching-plane field.'
    end if
  end subroutine solve_matching_root

  subroutine matching_branch_order(model, target_field_hat, order, count, status, message)
    character(len=*), intent(in) :: model
    real(dp), intent(in) :: target_field_hat
    character(len=1), intent(out) :: order(3)
    integer, intent(out) :: count
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    order = ' '
    count = 0
    status = matching_plane_zhao_ok
    message = ''
    select case (trim(model))
    case ('a')
      order(1) = 'A'
      count = 1
    case ('b')
      order(1) = 'B'
      count = 1
    case ('c')
      order(1) = 'C'
      count = 1
    case ('auto')
      if (target_field_hat > 0.0_dp) then
        order = ['A', 'B', 'C']
      else
        order = ['C', 'A', 'B']
      end if
      count = 3
    case default
      status = matching_plane_zhao_invalid_argument
      message = 'unknown matching-plane Zhao branch.'
    end select
  end subroutine matching_branch_order

  subroutine solve_one_matching_branch(params, branch, target_field_hat, root_selection, root, status, message)
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    character(len=*), intent(in) :: root_selection
    real(dp), intent(in) :: target_field_hat
    type(zhao_matching_root_type), intent(out) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp) :: guesses(3, 8), y(3), norm
    type(zhao_matching_root_type) :: candidate_root, unique_roots(8)
    integer :: guess_count, guess_index, iterations, unique_count, root_index
    integer(i32) :: profile_status
    logical :: success, compatible, duplicate_root
    logical :: saw_nonphysical_profile, saw_numerical_profile_failure
    character(len=512) :: profile_message

    root = zhao_matching_root_type()
    root%branch = branch
    status = matching_plane_zhao_no_physical_solution
    message = ''
    compatible = (branch == 'C' .and. target_field_hat < 0.0_dp) .or. &
                 ((branch == 'A' .or. branch == 'B') .and. target_field_hat > 0.0_dp)
    if (.not. compatible) then
      message = 'Zhao branch and matching-plane field signs are incompatible.'
      return
    end if
    call make_matching_branch_guesses(params, branch, guesses, guess_count)
    unique_count = 0
    saw_nonphysical_profile = .false.
    saw_numerical_profile_failure = .false.
    do guess_index = 1, guess_count
      call newton_matching_branch( &
        params, branch, target_field_hat, guesses(:, guess_index), y, norm, iterations, success &
        )
      if (.not. success) cycle
      candidate_root = zhao_matching_root_type()
      candidate_root%branch = branch
      call decode_matching_unknowns( &
        params, branch, y, candidate_root%phi0_v, candidate_root%phi_m_v, &
        candidate_root%ambient_electron_density_m3, success &
        )
      if (.not. success) cycle
      candidate_root%residual_norm = norm
      candidate_root%nonlinear_iterations = int(iterations, i32)
      call validate_matching_root_profile( &
        params, candidate_root, target_field_hat, profile_status, profile_message &
        )
      if (profile_status == matching_plane_zhao_no_physical_solution) then
        saw_nonphysical_profile = .true.
        cycle
      else if (profile_status /= matching_plane_zhao_ok) then
        saw_numerical_profile_failure = .true.
        cycle
      end if

      duplicate_root = .false.
      do root_index = 1, unique_count
        if (.not. matching_roots_equivalent(params, candidate_root, unique_roots(root_index))) cycle
        duplicate_root = .true.
        if (candidate_root%residual_norm < unique_roots(root_index)%residual_norm) then
          unique_roots(root_index) = candidate_root
        end if
        exit
      end do
      if (.not. duplicate_root) then
        unique_count = unique_count + 1
        unique_roots(unique_count) = candidate_root
      end if
    end do
    if (unique_count > 1) then
      if (trim(root_selection) == 'minimum_energy') then
        call select_minimum_energy_root(params, unique_roots, unique_count, root, status, message)
      else
        status = matching_plane_zhao_ambiguous_solution
        message = 'charge-driven Zhao solve found multiple roots in the requested branch.'
      end if
    else if (saw_numerical_profile_failure) then
      status = matching_plane_zhao_numerical_failure
      message = 'charge-driven Zhao root profile could not be certified numerically.'
    else if (unique_count == 1) then
      root = unique_roots(1)
      if (trim(root_selection) == 'minimum_energy') then
        call evaluate_root_potential_energy(params, root, status, message)
      else
        status = matching_plane_zhao_ok
        message = ''
      end if
    else if (saw_nonphysical_profile) then
      status = matching_plane_zhao_no_physical_solution
      message = 'charge-driven Zhao endpoint root has no real connecting field profile.'
    else
      status = matching_plane_zhao_numerical_failure
      message = 'charge-driven Zhao Newton solve did not converge.'
    end if
  end subroutine solve_one_matching_branch

  subroutine select_minimum_energy_root(params, roots, root_count, root, status, message)
    type(zhao_params_type), intent(in) :: params
    type(zhao_matching_root_type), intent(in) :: roots(:)
    integer, intent(in) :: root_count
    type(zhao_matching_root_type), intent(out) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    type(zhao_matching_root_type) :: candidates(size(roots))
    real(dp) :: energy_scale
    integer :: candidate_index, best_index

    root = zhao_matching_root_type()
    candidates = roots
    status = matching_plane_zhao_numerical_failure
    message = ''
    if (root_count < 1 .or. root_count > size(roots)) then
      message = 'matching-plane Zhao minimum-energy selection received an invalid candidate count.'
      return
    end if
    do candidate_index = 1, root_count
      call evaluate_root_potential_energy(params, candidates(candidate_index), status, message)
      if (status /= matching_plane_zhao_ok) return
    end do
    best_index = 1
    do candidate_index = 2, root_count
      if (candidates(candidate_index)%potential_energy_j_m2 < &
          candidates(best_index)%potential_energy_j_m2) best_index = candidate_index
    end do
    do candidate_index = 1, root_count
      if (candidate_index == best_index) cycle
      energy_scale = max( &
                     abs(candidates(best_index)%potential_energy_j_m2), &
                     abs(candidates(candidate_index)%potential_energy_j_m2), tiny(1.0_dp) &
                     )
      if (abs( &
          candidates(candidate_index)%potential_energy_j_m2 - &
          candidates(best_index)%potential_energy_j_m2 &
          ) <= energy_tie_tolerance*energy_scale) then
        status = matching_plane_zhao_ambiguous_solution
        message = 'matching-plane Zhao minimum-energy candidates are numerically tied.'
        return
      end if
    end do
    root = candidates(best_index)
    status = matching_plane_zhao_ok
    message = ''
  end subroutine select_minimum_energy_root

  subroutine evaluate_root_potential_energy(params, root, status, message)
    type(zhao_params_type), intent(in) :: params
    type(zhao_matching_root_type), intent(inout) :: root
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

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
  end subroutine evaluate_root_potential_energy

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
    integer :: point
    logical :: integral_ok

    energy_hat = 0.0_dp
    success = .false.
    if (.not. all(ieee_is_finite([ &
                                 start_phi_hat, end_phi_hat, phi0_hat, phi_m_hat, density_hat &
                                 ])) .or. density_hat <= 0.0_dp) return
    h = 1.0_dp/real(energy_quadrature_panels, dp)
    field_squared_scale = max(1.0_dp, phi0_hat*phi0_hat, phi_m_hat*phi_m_hat)
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
      if (.not. integral_ok .or. &
          field_squared < -profile_negative_tolerance*field_squared_scale) return
      summand = sqrt(max(0.0_dp, field_squared))*abs(jacobian)
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

  subroutine validate_matching_root_profile(params, root, target_field_hat, status, message)
    type(zhao_params_type), intent(in) :: params
    type(zhao_matching_root_type), intent(inout) :: root
    real(dp), intent(in) :: target_field_hat
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

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
  end subroutine validate_matching_root_profile

  pure logical function matching_roots_equivalent(params, first, second) result(equivalent)
    type(zhao_params_type), intent(in) :: params
    type(zhao_matching_root_type), intent(in) :: first, second

    real(dp) :: first_phi0_hat, second_phi0_hat, first_phi_m_hat, second_phi_m_hat
    real(dp) :: log_density_ratio

    equivalent = .false.
    if (first%branch /= second%branch) return
    if (min(first%ambient_electron_density_m3, second%ambient_electron_density_m3) <= 0.0_dp) return
    first_phi0_hat = first%phi0_v/params%t_phe_ev
    second_phi0_hat = second%phi0_v/params%t_phe_ev
    first_phi_m_hat = first%phi_m_v/params%t_phe_ev
    second_phi_m_hat = second%phi_m_v/params%t_phe_ev
    log_density_ratio = log(first%ambient_electron_density_m3/second%ambient_electron_density_m3)
    if (.not. all(ieee_is_finite([ &
                                 first_phi0_hat, second_phi0_hat, first_phi_m_hat, second_phi_m_hat, &
                                 log_density_ratio &
                                 ]))) return
    equivalent = abs(first_phi0_hat - second_phi0_hat) <= &
      root_cluster_tolerance*max(1.0_dp, abs(first_phi0_hat), abs(second_phi0_hat)) .and. &
      abs(first_phi_m_hat - second_phi_m_hat) <= &
      root_cluster_tolerance*max(1.0_dp, abs(first_phi_m_hat), abs(second_phi_m_hat)) .and. &
      abs(log_density_ratio) <= root_cluster_tolerance
  end function matching_roots_equivalent

  subroutine make_matching_branch_guesses(params, branch, guesses, count)
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(out) :: guesses(3, 8)
    integer, intent(out) :: count

    real(dp) :: density
    logical :: valid

    guesses = 0.0_dp
    select case (branch)
    case ('A')
      count = 5
      call encode_matching_unknowns(params, branch, 3.6_dp, -0.5_dp, 8.2e6_dp, guesses(:, 1), valid)
      call encode_matching_unknowns(params, branch, 2.8_dp, -0.3_dp, 8.0e6_dp, guesses(:, 2), valid)
      call encode_matching_unknowns(params, branch, 4.5_dp, -0.8_dp, 8.4e6_dp, guesses(:, 3), valid)
      call encode_matching_unknowns( &
        params, branch, 1.3_dp*params%t_phe_ev, -0.3_dp*params%t_phe_ev, &
        max(0.9_dp*params%n_swi_inf_m3, tiny(1.0_dp)), guesses(:, 4), valid &
        )
      call encode_matching_unknowns( &
        params, branch, 0.8_dp*params%t_phe_ev, -0.1_dp*params%t_phe_ev, &
        params%n_swi_inf_m3, guesses(:, 5), valid &
        )
    case ('B')
      count = 7
      call encode_matching_unknowns(params, branch, 1.3_dp, 1.3_dp, 7.0e6_dp, guesses(:, 1), valid)
      call encode_matching_unknowns(params, branch, 0.8_dp, 0.8_dp, 6.5e6_dp, guesses(:, 2), valid)
      call encode_matching_unknowns(params, branch, 2.0_dp, 2.0_dp, 7.8e6_dp, guesses(:, 3), valid)
      call encode_matching_unknowns( &
        params, branch, 0.6_dp*params%t_phe_ev, 0.6_dp*params%t_phe_ev, &
        params%n_swi_inf_m3, guesses(:, 4), valid &
        )
      call encode_matching_unknowns( &
        params, branch, 0.2_dp*params%t_phe_ev, 0.2_dp*params%t_phe_ev, &
        0.8_dp*params%n_swi_inf_m3, guesses(:, 5), valid &
        )
      density = max( &
                (2.0_dp*params%n_swi_inf_m3 - params%n_phe0_m3)/(1.0_dp + erf(params%u)), &
                0.5_dp*params%n_swi_inf_m3 &
                )
      call encode_matching_unknowns( &
        params, branch, 0.02_dp*params%t_phe_ev, 0.02_dp*params%t_phe_ev, &
        density, guesses(:, 6), valid &
        )
      call encode_matching_unknowns( &
        params, branch, 0.002_dp*params%t_phe_ev, 0.002_dp*params%t_phe_ev, &
        density, guesses(:, 7), valid &
        )
    case ('C')
      count = 7
      call encode_matching_unknowns(params, branch, -0.5_dp, -0.5_dp, 6.0e6_dp, guesses(:, 1), valid)
      call encode_matching_unknowns(params, branch, -2.0_dp, -2.0_dp, 7.0e6_dp, guesses(:, 2), valid)
      call encode_matching_unknowns(params, branch, -5.0_dp, -5.0_dp, 8.0e6_dp, guesses(:, 3), valid)
      call encode_matching_unknowns(params, branch, -10.0_dp, -10.0_dp, 8.2e6_dp, guesses(:, 4), valid)
      call encode_matching_unknowns(params, branch, -15.0_dp, -15.0_dp, 8.5e6_dp, guesses(:, 5), valid)
      call encode_matching_unknowns( &
        params, branch, -params%t_phe_ev, -params%t_phe_ev, &
        params%n_swi_inf_m3, guesses(:, 6), valid &
        )
      call encode_matching_unknowns( &
        params, branch, -3.0_dp*params%t_phe_ev, -3.0_dp*params%t_phe_ev, &
        0.9_dp*params%n_swi_inf_m3, guesses(:, 7), valid &
        )
    case default
      count = 0
    end select
  end subroutine make_matching_branch_guesses

  subroutine encode_matching_unknowns(params, branch, phi0_v, phi_m_v, density_m3, y, valid)
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: phi0_v, phi_m_v, density_m3
    real(dp), intent(out) :: y(3)
    logical, intent(out) :: valid

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
  end subroutine encode_matching_unknowns

  subroutine decode_matching_unknowns(params, branch, y, phi0_v, phi_m_v, density_m3, valid)
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: y(3)
    real(dp), intent(out) :: phi0_v, phi_m_v, density_m3
    logical, intent(out) :: valid

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
  end subroutine decode_matching_unknowns

  subroutine newton_matching_branch( &
    params, branch, target_field_hat, y0, y_out, final_norm, iterations, success &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: target_field_hat, y0(3)
    real(dp), intent(out) :: y_out(3), final_norm
    integer, intent(out) :: iterations
    logical, intent(out) :: success

    real(dp) :: y(3), f(3), jac(3, 3), delta(3), trial(3), trial_f(3)
    real(dp) :: norm, trial_norm, step
    integer :: n, iteration, backtrack
    logical :: valid, jacobian_ok, linear_ok, trial_valid

    n = merge(3, 2, branch == 'A')
    y = y0
    call evaluate_charge_residual(params, branch, target_field_hat, y, f, valid)
    if (.not. valid) then
      y_out = y
      final_norm = huge(1.0_dp)
      iterations = 0
      success = .false.
      return
    end if
    norm = maxval(abs(f(1:n)))
    success = .false.
    do iteration = 0, root_max_iterations
      if (norm <= root_tolerance) then
        success = .true.
        exit
      end if
      if (iteration == root_max_iterations) exit
      call matching_numerical_jacobian( &
        params, branch, target_field_hat, y, f, n, jac, jacobian_ok &
        )
      if (.not. jacobian_ok) exit
      call solve_matching_small_system(jac, -f, n, delta, linear_ok)
      if (.not. linear_ok) exit
      step = 1.0_dp
      do backtrack = 1, root_max_backtracks
        trial = y + step*delta
        call evaluate_charge_residual( &
          params, branch, target_field_hat, trial, trial_f, trial_valid &
          )
        if (trial_valid) then
          trial_norm = maxval(abs(trial_f(1:n)))
          if (trial_norm < norm) then
            y = trial
            f = trial_f
            norm = trial_norm
            exit
          end if
        end if
        step = 0.5_dp*step
      end do
      if (backtrack > root_max_backtracks) exit
    end do
    y_out = y
    final_norm = norm
    iterations = iteration
  end subroutine newton_matching_branch

  subroutine matching_numerical_jacobian( &
    params, branch, target_field_hat, y, f0, n, jac, success &
    )
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: target_field_hat, y(3), f0(3)
    integer, intent(in) :: n
    real(dp), intent(out) :: jac(3, 3)
    logical, intent(out) :: success

    real(dp) :: yp(3), ym(3), fp(3), fm(3), h
    integer :: column
    logical :: plus_valid, minus_valid

    jac = 0.0_dp
    success = .true.
    do column = 1, n
      h = epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp, abs(y(column)))
      yp = y
      ym = y
      yp(column) = yp(column) + h
      ym(column) = ym(column) - h
      call evaluate_charge_residual(params, branch, target_field_hat, yp, fp, plus_valid)
      call evaluate_charge_residual(params, branch, target_field_hat, ym, fm, minus_valid)
      if (plus_valid .and. minus_valid) then
        jac(1:n, column) = (fp(1:n) - fm(1:n))/(2.0_dp*h)
      else if (plus_valid) then
        jac(1:n, column) = (fp(1:n) - f0(1:n))/h
      else if (minus_valid) then
        jac(1:n, column) = (f0(1:n) - fm(1:n))/h
      else
        success = .false.
        return
      end if
    end do
    success = all(ieee_is_finite(jac(1:n, 1:n)))
  end subroutine matching_numerical_jacobian

  subroutine solve_matching_small_system(a_in, b_in, n, x, success)
    real(dp), intent(in) :: a_in(3, 3), b_in(3)
    integer, intent(in) :: n
    real(dp), intent(out) :: x(3)
    logical, intent(out) :: success

    real(dp) :: a(3, 3), b(3), factor, pivot_value, tmp
    integer :: i, j, k, pivot

    a = a_in
    b = b_in
    x = 0.0_dp
    success = .false.
    do k = 1, n
      pivot = k
      do i = k + 1, n
        if (abs(a(i, k)) > abs(a(pivot, k))) pivot = i
      end do
      if (.not. ieee_is_finite(a(pivot, k)) .or. abs(a(pivot, k)) <= 1.0e-14_dp) return
      if (pivot /= k) then
        do j = k, n
          tmp = a(k, j)
          a(k, j) = a(pivot, j)
          a(pivot, j) = tmp
        end do
        tmp = b(k)
        b(k) = b(pivot)
        b(pivot) = tmp
      end if
      pivot_value = a(k, k)
      do i = k + 1, n
        factor = a(i, k)/pivot_value
        a(i, k:n) = a(i, k:n) - factor*a(k, k:n)
        b(i) = b(i) - factor*b(k)
      end do
    end do
    do i = n, 1, -1
      x(i) = b(i)
      do j = i + 1, n
        x(i) = x(i) - a(i, j)*x(j)
      end do
      x(i) = x(i)/a(i, i)
    end do
    success = all(ieee_is_finite(x(1:n)))
  end subroutine solve_matching_small_system

  subroutine evaluate_charge_residual(params, branch, target_field_hat, y, residual, valid)
    type(zhao_params_type), intent(in) :: params
    character(len=1), intent(in) :: branch
    real(dp), intent(in) :: target_field_hat, y(3)
    real(dp), intent(out) :: residual(3)
    logical, intent(out) :: valid

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
  end subroutine evaluate_charge_residual

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

  subroutine compose_matching_response(params, root, output, status, message)
    type(zhao_params_type), intent(in) :: params
    type(zhao_matching_root_type), intent(in) :: root
    real(dp), intent(out) :: output(matching_plane_response_output_count)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

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
  end subroutine compose_matching_response

  pure logical function ion_accessible(params, phi_hat) result(accessible)
    type(zhao_params_type), intent(in) :: params
    real(dp), intent(in) :: phi_hat

    accessible = params%tau > 0.0_dp .and. params%mach > 0.0_dp .and. &
                 1.0_dp - 2.0_dp*phi_hat/(params%tau*params%mach*params%mach) > 0.0_dp
  end function ion_accessible

  subroutine assign_diagnostics(destination, source)
    type(matching_plane_zhao_diagnostics_type), intent(out), optional :: destination
    type(matching_plane_zhao_diagnostics_type), intent(in) :: source

    if (present(destination)) destination = source
  end subroutine assign_diagnostics

end module bem_matching_plane_zhao

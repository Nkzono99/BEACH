!> Matching-plane 用の charge-driven Zhao 外部シース応答。
!!
!! 応答は整合面直下の D_z/eps0 を外部側 interface field とみなし、
!! 上流の準中性条件と Sagdeev 積分を満たす A/B/C branch を解く。
!! 定常表面電流の零電流条件は課さないため、BEACH 側の帯電過程を消去しない。
!! 平面・半無限の外部問題は z 方向に並進対称なので、絶対高度 H は数値式に
!! 入らず、runtime がこの応答を domain の z-high gauge へ結び付ける。
!! 既定policyはqueryごとにstatelessで、opt-inのType A continuationだけが
!! 呼出側から渡されたaccepted endpoint rootをseedとして使う。
!! この module は公開型と評価の入口を持ち、物理式を physics、Newton 法を
!! numerics、根の選択と継続を roots の非公開 submodule へ委譲する。
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
    zhao_params_type, swe_free_current_term
  use bem_string_utils, only: lower_ascii
  implicit none
  private

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
    logical :: continuation_used = .false.
    logical :: continuation_fallback_used = .false.
    real(dp) :: continuation_root_jump = 0.0_dp
  end type matching_plane_zhao_diagnostics_type

  type, public :: matching_plane_zhao_root_seed_type
    logical :: valid = .false.
    real(dp) :: phi0_v = 0.0_dp
    real(dp) :: phi_m_v = 0.0_dp
    real(dp) :: ambient_electron_density_m3 = 0.0_dp
  end type matching_plane_zhao_root_seed_type

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
    procedure, public :: reconstruct_seed => reconstruct_matching_plane_zhao_type_a_seed
    procedure, public :: get_feedback_scales => get_matching_plane_zhao_feedback_scales
    procedure, public :: is_initialized => matching_plane_zhao_is_initialized
  end type matching_plane_zhao_model_type

  ! Private entry points shared by the three implementation submodules.
  interface
    module subroutine solve_matching_root(model, root_selection, params, interface_field_v_m, root, status, message)
      character(len=*), intent(in) :: model, root_selection
      type(zhao_params_type), intent(in) :: params
      real(dp), intent(in) :: interface_field_v_m
      type(zhao_matching_root_type), intent(out) :: root
      integer(i32), intent(out) :: status
      character(len=*), intent(out) :: message
    end subroutine solve_matching_root

    module subroutine solve_matching_type_a_continuation( &
      params, interface_field_v_m, seed, root, fallback_used, root_jump, status, message &
      )
      type(zhao_params_type), intent(in) :: params
      real(dp), intent(in) :: interface_field_v_m
      type(matching_plane_zhao_root_seed_type), intent(in) :: seed
      type(zhao_matching_root_type), intent(out) :: root
      logical, intent(out) :: fallback_used
      real(dp), intent(out) :: root_jump
      integer(i32), intent(out) :: status
      character(len=*), intent(out) :: message
    end subroutine solve_matching_type_a_continuation

    module subroutine make_matching_branch_guesses(params, branch, guesses, count)
      type(zhao_params_type), intent(in) :: params
      character(len=1), intent(in) :: branch
      real(dp), intent(out) :: guesses(3, 8)
      integer, intent(out) :: count
    end subroutine make_matching_branch_guesses

    module subroutine newton_matching_branch( &
      params, branch, target_field_hat, y0, y_out, final_norm, iterations, success &
      )
      type(zhao_params_type), intent(in) :: params
      character(len=1), intent(in) :: branch
      real(dp), intent(in) :: target_field_hat, y0(3)
      real(dp), intent(out) :: y_out(3), final_norm
      integer, intent(out) :: iterations
      logical, intent(out) :: success
    end subroutine newton_matching_branch

    module subroutine prepare_matching_zhao_query( &
      self, input, params, photoelectron_temperature_ev, photoelectron_source_density_m3, status, message &
      )
      class(matching_plane_zhao_model_type), intent(in) :: self
      real(dp), intent(in) :: input(matching_plane_response_input_count)
      type(zhao_params_type), intent(out) :: params
      real(dp), intent(out) :: photoelectron_temperature_ev, photoelectron_source_density_m3
      integer(i32), intent(out) :: status
      character(len=*), intent(out) :: message
    end subroutine prepare_matching_zhao_query

    module subroutine encode_matching_unknowns(params, branch, phi0_v, phi_m_v, density_m3, y, valid)
      type(zhao_params_type), intent(in) :: params
      character(len=1), intent(in) :: branch
      real(dp), intent(in) :: phi0_v, phi_m_v, density_m3
      real(dp), intent(out) :: y(3)
      logical, intent(out) :: valid
    end subroutine encode_matching_unknowns

    module subroutine decode_matching_unknowns(params, branch, y, phi0_v, phi_m_v, density_m3, valid)
      type(zhao_params_type), intent(in) :: params
      character(len=1), intent(in) :: branch
      real(dp), intent(in) :: y(3)
      real(dp), intent(out) :: phi0_v, phi_m_v, density_m3
      logical, intent(out) :: valid
    end subroutine decode_matching_unknowns

    module subroutine evaluate_charge_residual(params, branch, target_field_hat, y, residual, valid)
      type(zhao_params_type), intent(in) :: params
      character(len=1), intent(in) :: branch
      real(dp), intent(in) :: target_field_hat, y(3)
      real(dp), intent(out) :: residual(3)
      logical, intent(out) :: valid
    end subroutine evaluate_charge_residual

    module subroutine validate_matching_root_profile(params, root, target_field_hat, status, message)
      type(zhao_params_type), intent(in) :: params
      type(zhao_matching_root_type), intent(inout) :: root
      real(dp), intent(in) :: target_field_hat
      integer(i32), intent(out) :: status
      character(len=*), intent(out) :: message
    end subroutine validate_matching_root_profile

    module subroutine evaluate_root_potential_energy(params, root, status, message)
      type(zhao_params_type), intent(in) :: params
      type(zhao_matching_root_type), intent(inout) :: root
      integer(i32), intent(out) :: status
      character(len=*), intent(out) :: message
    end subroutine evaluate_root_potential_energy

    module subroutine compose_matching_response(params, root, output, status, message)
      type(zhao_params_type), intent(in) :: params
      type(zhao_matching_root_type), intent(in) :: root
      real(dp), intent(out) :: output(matching_plane_response_output_count)
      integer(i32), intent(out) :: status
      character(len=*), intent(out) :: message
    end subroutine compose_matching_response
  end interface

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
    case ('require_unique', 'minimum_energy', 'continuation')
      self%root_selection = normalized_root_selection
    case default
      message = 'matching-plane Zhao root selection must be require_unique, minimum_energy, or continuation.'
      return
    end select
    if (self%root_selection == 'continuation' .and. self%branch_model /= 'a') then
      message = 'matching-plane Zhao continuation requires an explicit Type-A branch.'
      return
    end if
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

  subroutine evaluate_matching_plane_zhao( &
    self, input, output, status, message, diagnostics, continuation_seed, continuation_candidate &
    )
    class(matching_plane_zhao_model_type), intent(in) :: self
    real(dp), intent(in) :: input(matching_plane_response_input_count)
    real(dp), intent(out) :: output(matching_plane_response_output_count)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(matching_plane_zhao_diagnostics_type), intent(out), optional :: diagnostics
    type(matching_plane_zhao_root_seed_type), intent(in), optional :: continuation_seed
    type(matching_plane_zhao_root_seed_type), intent(out), optional :: continuation_candidate

    type(zhao_params_type) :: params
    type(zhao_matching_root_type) :: root
    type(matching_plane_zhao_diagnostics_type) :: local_diagnostics
    real(dp) :: interface_field_v_m
    real(dp) :: photoelectron_temperature_ev, photoelectron_source_density_m3
    real(dp) :: continuation_root_jump
    logical :: continuation_fallback_used, have_continuation_seed

    output = 0.0_dp
    if (present(continuation_candidate)) continuation_candidate = matching_plane_zhao_root_seed_type()
    status = matching_plane_zhao_invalid_argument
    message = ''
    local_diagnostics = matching_plane_zhao_diagnostics_type()
    call prepare_matching_zhao_query( &
      self, input, params, photoelectron_temperature_ev, photoelectron_source_density_m3, status, message &
      )
    if (status /= matching_plane_zhao_ok) then
      call assign_diagnostics(diagnostics, local_diagnostics)
      return
    end if

    interface_field_v_m = input(matching_plane_input_displacement)/eps0
    local_diagnostics%interface_field_v_m = interface_field_v_m
    local_diagnostics%effective_photoelectron_temperature_ev = photoelectron_temperature_ev
    local_diagnostics%photoelectron_source_density_m3 = photoelectron_source_density_m3
    have_continuation_seed = .false.
    if (present(continuation_seed)) have_continuation_seed = continuation_seed%valid
    if (self%root_selection == 'continuation' .and. have_continuation_seed) then
      local_diagnostics%continuation_used = .true.
      call solve_matching_type_a_continuation( &
        params, interface_field_v_m, continuation_seed, root, continuation_fallback_used, &
        continuation_root_jump, status, message &
        )
      local_diagnostics%continuation_fallback_used = continuation_fallback_used
      local_diagnostics%continuation_root_jump = continuation_root_jump
    else if (self%root_selection == 'continuation') then
      call solve_matching_root( &
        'a', 'minimum_energy', params, interface_field_v_m, root, status, message &
        )
    else
      call solve_matching_root( &
        trim(self%branch_model), trim(self%root_selection), params, interface_field_v_m, root, status, message &
        )
    end if
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
    if (present(continuation_candidate) .and. root%branch == 'A') then
      continuation_candidate%valid = .true.
      continuation_candidate%phi0_v = root%phi0_v
      continuation_candidate%phi_m_v = root%phi_m_v
      continuation_candidate%ambient_electron_density_m3 = root%ambient_electron_density_m3
    end if
    call assign_diagnostics(diagnostics, local_diagnostics)
  end subroutine evaluate_matching_plane_zhao

  subroutine reconstruct_matching_plane_zhao_type_a_seed( &
    self, input, response, seed, status, message &
    )
    class(matching_plane_zhao_model_type), intent(in) :: self
    real(dp), intent(in) :: input(matching_plane_response_input_count)
    real(dp), intent(in) :: response(matching_plane_response_output_count)
    type(matching_plane_zhao_root_seed_type), intent(out) :: seed
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    type(zhao_params_type) :: params
    real(dp) :: photoelectron_temperature_ev, photoelectron_source_density_m3
    real(dp) :: electron_cutoff, unit_density_electron_term, number_flux_scale, flux_coefficient

    seed = matching_plane_zhao_root_seed_type()
    status = matching_plane_zhao_invalid_argument
    message = ''
    if (.not. self%initialized) then
      message = 'matching-plane Zhao model is not initialized.'
      return
    end if
    if (self%branch_model /= 'a') then
      message = 'matching-plane Zhao Type-A seed reconstruction requires an explicit Type-A model.'
      return
    end if
    if (.not. all(ieee_is_finite(response))) then
      message = 'matching-plane Zhao restart response must be finite.'
      return
    end if

    call prepare_matching_zhao_query( &
      self, input, params, photoelectron_temperature_ev, photoelectron_source_density_m3, status, message &
      )
    if (status /= matching_plane_zhao_ok) return

    seed%phi0_v = response(matching_plane_output_matching_potential)
    seed%phi_m_v = response(matching_plane_output_electron_access_potential)
    electron_cutoff = sqrt(max(0.0_dp, -seed%phi_m_v/params%t_swe_ev)) - params%u
    unit_density_electron_term = swe_free_current_term(params, 1.0_dp, electron_cutoff)
    number_flux_scale = params%v_phe_th_mps/(2.0_dp*sqrt(pi))
    flux_coefficient = number_flux_scale*unit_density_electron_term
    if (.not. all(ieee_is_finite([electron_cutoff, flux_coefficient])) .or. &
        flux_coefficient <= 0.0_dp .or. &
        response(matching_plane_output_electron_inward_flux) <= 0.0_dp) then
      seed = matching_plane_zhao_root_seed_type()
      message = 'matching-plane Zhao restart response cannot reconstruct a positive electron density.'
      return
    end if
    seed%ambient_electron_density_m3 = &
      response(matching_plane_output_electron_inward_flux)/flux_coefficient
    seed%valid = seed%phi0_v > 0.0_dp .and. seed%phi_m_v < 0.0_dp .and. &
                 ieee_is_finite(seed%ambient_electron_density_m3) .and. &
                 seed%ambient_electron_density_m3 > 0.0_dp
    if (.not. seed%valid) then
      seed = matching_plane_zhao_root_seed_type()
      message = 'matching-plane Zhao restart response is not a valid Type-A root seed.'
      return
    end if
    status = matching_plane_zhao_ok
  end subroutine reconstruct_matching_plane_zhao_type_a_seed

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

  subroutine assign_diagnostics(destination, source)
    type(matching_plane_zhao_diagnostics_type), intent(out), optional :: destination
    type(matching_plane_zhao_diagnostics_type), intent(in) :: source

    if (present(destination)) destination = source
  end subroutine assign_diagnostics

end module bem_matching_plane_zhao

!> Matching-plane 応答の table / online Zhao backend を同じ契約で提供する。
module bem_matching_plane_response_provider
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_app_config_types, only: app_config
  use bem_matching_plane_response, only: matching_plane_response_table_type, &
                                         matching_plane_response_input_count, &
                                         matching_plane_response_output_count, &
                                         matching_plane_response_ok
  use bem_matching_plane_zhao, only: matching_plane_zhao_model_type, matching_plane_zhao_root_seed_type, &
                                     matching_plane_zhao_ok, &
                                     matching_plane_zhao_diagnostics_type, &
                                     matching_plane_zhao_invalid_argument, &
                                     matching_plane_zhao_no_physical_solution, &
                                     matching_plane_zhao_ambiguous_solution
  use bem_mpi, only: mpi_context
  implicit none
  private

  integer(i32), parameter, public :: matching_plane_provider_ok = 0_i32
  integer(i32), parameter, public :: matching_plane_provider_invalid_argument = 1_i32
  integer(i32), parameter, public :: matching_plane_provider_load_failure = 2_i32
  integer(i32), parameter, public :: matching_plane_provider_no_physical_solution = 3_i32
  integer(i32), parameter, public :: matching_plane_provider_numerical_failure = 4_i32
  integer(i32), parameter, public :: matching_plane_provider_ambiguous_solution = 5_i32

  integer(i32), parameter :: provider_backend_none = 0_i32
  integer(i32), parameter :: provider_backend_table = 1_i32
  integer(i32), parameter :: provider_backend_zhao_online = 2_i32

  type, public :: matching_plane_response_provider_type
    private
    logical :: active = .false.
    integer(i32) :: backend = provider_backend_none
    real(dp) :: matching_plane_z_m = 0.0_dp
    real(dp) :: feedback_min(4) = 0.0_dp
    real(dp) :: feedback_max(4) = 0.0_dp
    real(dp) :: feedback_scale(4) = 0.0_dp
    logical :: feedback_bounded(4) = .false.
    real(dp) :: implicit_displacement_min = 0.0_dp
    real(dp) :: implicit_displacement_max = 0.0_dp
    real(dp) :: implicit_displacement_scale = 0.0_dp
    real(dp) :: implicit_feedback_reference(4) = 0.0_dp
    logical :: implicit_displacement_bounded = .false.
    logical :: implicit_zero_mode_supported = .false.
    character(len=16) :: content_fingerprint = ''
    type(matching_plane_response_table_type) :: table
    type(matching_plane_zhao_model_type) :: zhao
  contains
    procedure, public :: initialize => initialize_matching_plane_response_provider
    procedure, public :: evaluate => evaluate_matching_plane_response_provider
    procedure, public :: evaluate_local => evaluate_matching_plane_response_provider_local
    procedure, public :: reconstruct_continuation_seed_local => &
      reconstruct_matching_plane_continuation_seed_local
    procedure, public :: evaluate_zhao_local => evaluate_matching_plane_zhao_provider_local
    procedure, public :: validate_feedback => validate_matching_plane_provider_feedback
    procedure, public :: feedback_converged => matching_plane_provider_feedback_converged
    procedure, public :: feedback_residual => matching_plane_provider_feedback_residual
    procedure, public :: get_matching_plane_z => get_provider_matching_plane_z
    procedure, public :: get_feedback_scales => get_provider_feedback_scales
    procedure, public :: get_implicit_zero_mode_contract => get_provider_implicit_zero_mode_contract
    procedure, public :: get_backend_name => get_provider_backend_name
    procedure, public :: get_content_fingerprint => get_provider_content_fingerprint
    procedure, public :: is_active => matching_plane_provider_is_active
  end type matching_plane_response_provider_type

  interface
    !> 設定を検証済み backend snapshot へ変換する。
    module subroutine initialize_matching_plane_response_provider(self, cfg, mpi, status, message)
      class(matching_plane_response_provider_type), intent(inout) :: self
      type(app_config), intent(in) :: cfg
      type(mpi_context), intent(in) :: mpi
      integer(i32), intent(out) :: status
      character(len=*), intent(out) :: message
    end subroutine initialize_matching_plane_response_provider
    !> 全rankで同じ query を確認し、online Zhao はrootだけで解いてbroadcastする。
    module subroutine evaluate_matching_plane_response_provider(self, input, mpi, output, status, message)
      class(matching_plane_response_provider_type), intent(inout) :: self
      real(dp), intent(in) :: input(matching_plane_response_input_count)
      type(mpi_context), intent(in) :: mpi
      real(dp), intent(out) :: output(matching_plane_response_output_count)
      integer(i32), intent(out) :: status
      character(len=*), intent(out) :: message
    end subroutine evaluate_matching_plane_response_provider

    module subroutine map_zhao_failure(zhao_status, zhao_message, status, message)
      integer(i32), intent(in) :: zhao_status
      character(len=*), intent(in) :: zhao_message
      integer(i32), intent(out) :: status
      character(len=*), intent(out) :: message
    end subroutine map_zhao_failure

    module subroutine accept_provider(status, message)
      integer(i32), intent(out) :: status
      character(len=*), intent(out) :: message
    end subroutine accept_provider

    module subroutine reject_provider(code, text, status, message)
      integer(i32), intent(in) :: code
      character(len=*), intent(in) :: text
      integer(i32), intent(out) :: status
      character(len=*), intent(out) :: message
    end subroutine reject_provider
  end interface

contains

  !> MPI collectiveを伴わず、呼出rankだけで応答を評価する。
  !! implicit zero-mode root はMPI root上でこの入口を反復し、最終結果だけをbroadcastする。
  subroutine evaluate_matching_plane_response_provider_local( &
    self, input, output, status, message, continuation_seed, continuation_candidate &
    )
    class(matching_plane_response_provider_type), intent(inout) :: self
    real(dp), intent(in) :: input(matching_plane_response_input_count)
    real(dp), intent(out) :: output(matching_plane_response_output_count)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(matching_plane_zhao_root_seed_type), intent(in), optional :: continuation_seed
    type(matching_plane_zhao_root_seed_type), intent(out), optional :: continuation_candidate

    integer(i32) :: backend_status
    character(len=512) :: backend_message

    output = 0.0_dp
    if (present(continuation_candidate)) continuation_candidate = matching_plane_zhao_root_seed_type()
    call accept_provider(status, message)
    backend_message = ''
    if (.not. self%active .or. self%backend == provider_backend_none .or. &
        any(.not. ieee_is_finite(input)) .or. any(input(2:5) < 0.0_dp)) then
      call reject_provider( &
        matching_plane_provider_invalid_argument, &
        'matching-plane response provider or query is invalid.', status, message &
        )
      return
    end if

    select case (self%backend)
    case (provider_backend_table)
      call self%table%evaluate(input, output, backend_status, backend_message)
      if (backend_status /= matching_plane_response_ok) then
        call reject_provider( &
          matching_plane_provider_load_failure, trim(backend_message), status, message &
          )
        return
      end if
    case (provider_backend_zhao_online)
      call self%zhao%evaluate( &
        input, output, backend_status, backend_message, &
        continuation_seed=continuation_seed, continuation_candidate=continuation_candidate &
        )
      if (backend_status /= matching_plane_zhao_ok) then
        call map_zhao_failure(backend_status, backend_message, status, message)
        return
      end if
    case default
      call reject_provider( &
        matching_plane_provider_invalid_argument, &
        'matching-plane response provider is not initialized.', status, message &
        )
      return
    end select

    if (any(.not. ieee_is_finite(output)) .or. any(output(2:3) < 0.0_dp)) then
      output = 0.0_dp
      call reject_provider( &
        matching_plane_provider_numerical_failure, &
        'matching-plane response output is invalid.', status, message &
        )
    end if
  end subroutine evaluate_matching_plane_response_provider_local

  !> 保存済みの有限なonline Zhao応答からcontinuation seedを再構成する。
  subroutine reconstruct_matching_plane_continuation_seed_local( &
    self, input, output, seed, status, message &
    )
    class(matching_plane_response_provider_type), intent(in) :: self
    real(dp), intent(in) :: input(matching_plane_response_input_count)
    real(dp), intent(in) :: output(matching_plane_response_output_count)
    type(matching_plane_zhao_root_seed_type), intent(out) :: seed
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    integer(i32) :: backend_status
    character(len=512) :: backend_message

    seed = matching_plane_zhao_root_seed_type()
    call accept_provider(status, message)
    if (.not. self%active .or. self%backend /= provider_backend_zhao_online) then
      call reject_provider( &
        matching_plane_provider_invalid_argument, &
        'continuation seed reconstruction requires an initialized online Zhao provider.', status, message &
        )
      return
    end if

    call self%zhao%reconstruct_seed(input, output, seed, backend_status, backend_message)
    if (backend_status /= matching_plane_zhao_ok) then
      call map_zhao_failure(backend_status, backend_message, status, message)
    end if
  end subroutine reconstruct_matching_plane_continuation_seed_local

  !> Online Zhao backendをMPI collectiveなしで評価し、Zhao固有statusを保持する。
  !! Solvability atlas専用の入口であり、通常runtimeは共通evaluate APIを使う。
  subroutine evaluate_matching_plane_zhao_provider_local( &
    self, input, output, status, message, diagnostics &
    )
    class(matching_plane_response_provider_type), intent(in) :: self
    real(dp), intent(in) :: input(matching_plane_response_input_count)
    real(dp), intent(out) :: output(matching_plane_response_output_count)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    type(matching_plane_zhao_diagnostics_type), intent(out) :: diagnostics

    output = 0.0_dp
    diagnostics = matching_plane_zhao_diagnostics_type()
    status = matching_plane_zhao_invalid_argument
    message = ''
    if (.not. self%active .or. self%backend /= provider_backend_zhao_online) then
      message = 'matching-plane provider is not an initialized online Zhao backend.'
      return
    end if
    call self%zhao%evaluate(input, output, status, message, diagnostics)
  end subroutine evaluate_matching_plane_zhao_provider_local

  subroutine validate_matching_plane_provider_feedback(self, feedback, status, message)
    class(matching_plane_response_provider_type), intent(in) :: self
    real(dp), intent(in) :: feedback(4)
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    integer(i32) :: axis
    real(dp) :: tolerance

    call accept_provider(status, message)
    if (any(.not. ieee_is_finite(feedback)) .or. any(feedback < 0.0_dp)) then
      call reject_provider( &
        matching_plane_provider_invalid_argument, &
        'matching-plane observed feedback must be finite and nonnegative.', status, message &
        )
      return
    end if
    do axis = 1_i32, 4_i32
      if (.not. self%feedback_bounded(axis)) cycle
      tolerance = 64.0_dp*epsilon(1.0_dp)*max( &
                  1.0_dp, abs(self%feedback_min(axis)), abs(self%feedback_max(axis)) &
                  )
      if (feedback(axis) < self%feedback_min(axis) - tolerance .or. &
          feedback(axis) > self%feedback_max(axis) + tolerance) then
        call reject_provider( &
          matching_plane_provider_invalid_argument, &
          'matching-plane observed feedback is outside the response backend domain.', &
          status, message &
          )
        return
      end if
    end do
  end subroutine validate_matching_plane_provider_feedback

  logical function matching_plane_provider_feedback_converged( &
    self, previous, observed, relative_tolerance, absolute_tolerance &
    ) result(converged)
    class(matching_plane_response_provider_type), intent(in) :: self
    real(dp), intent(in) :: previous(4), observed(4)
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance(4)
    integer(i32) :: axis
    real(dp) :: threshold

    converged = .true.
    do axis = 1_i32, 4_i32
      if (self%feedback_scale(axis) <= 0.0_dp) cycle
      threshold = max( &
                  relative_tolerance*self%feedback_scale(axis), &
                  absolute_tolerance(axis) &
                  )
      if (abs(observed(axis) - previous(axis)) > threshold) then
        converged = .false.
        return
      end if
    end do
  end function matching_plane_provider_feedback_converged

  real(dp) function matching_plane_provider_feedback_residual( &
    self, previous, observed, relative_tolerance, absolute_tolerance &
    ) result(residual)
    class(matching_plane_response_provider_type), intent(in) :: self
    real(dp), intent(in) :: previous(4), observed(4)
    real(dp), intent(in), optional :: relative_tolerance, absolute_tolerance(4)
    integer(i32) :: axis
    real(dp) :: absolute_defect
    logical :: use_mixed_tolerance

    residual = 0.0_dp
    use_mixed_tolerance = present(relative_tolerance) .and. present(absolute_tolerance)
    if (present(relative_tolerance) .neqv. present(absolute_tolerance)) then
      error stop 'matching-plane residual requires both relative and absolute tolerances.'
    end if
    do axis = 1_i32, 4_i32
      if (self%feedback_scale(axis) <= 0.0_dp) cycle
      absolute_defect = abs(observed(axis) - previous(axis))
      if (use_mixed_tolerance) then
        if (absolute_tolerance(axis) > relative_tolerance*self%feedback_scale(axis)) then
          residual = max( &
                     residual, relative_tolerance*(absolute_defect/absolute_tolerance(axis)) &
                     )
        else
          residual = max(residual, absolute_defect/self%feedback_scale(axis))
        end if
      else
        residual = max(residual, absolute_defect/self%feedback_scale(axis))
      end if
    end do
  end function matching_plane_provider_feedback_residual

  subroutine get_provider_matching_plane_z(self, matching_plane_z_m, status, message)
    class(matching_plane_response_provider_type), intent(in) :: self
    real(dp), intent(out) :: matching_plane_z_m
    integer(i32), intent(out), optional :: status
    character(len=*), intent(out), optional :: message

    matching_plane_z_m = self%matching_plane_z_m
    if (present(status)) status = matching_plane_provider_ok
    if (present(message)) message = ''
    if (self%active .and. self%backend /= provider_backend_none) return
    matching_plane_z_m = 0.0_dp
    if (present(status)) status = matching_plane_provider_invalid_argument
    if (present(message)) message = 'matching-plane response provider is not initialized.'
  end subroutine get_provider_matching_plane_z

  subroutine get_provider_feedback_scales(self, scales)
    class(matching_plane_response_provider_type), intent(in) :: self
    real(dp), intent(out) :: scales(4)

    scales = self%feedback_scale
  end subroutine get_provider_feedback_scales

  subroutine get_provider_implicit_zero_mode_contract( &
    self, supported, displacement_bounded, displacement_min, displacement_max, displacement_scale, &
    feedback_reference &
    )
    class(matching_plane_response_provider_type), intent(in) :: self
    logical, intent(out) :: supported, displacement_bounded
    real(dp), intent(out) :: displacement_min, displacement_max, displacement_scale, feedback_reference(4)

    supported = self%active .and. self%implicit_zero_mode_supported
    displacement_bounded = self%implicit_displacement_bounded
    displacement_min = self%implicit_displacement_min
    displacement_max = self%implicit_displacement_max
    displacement_scale = self%implicit_displacement_scale
    feedback_reference = self%implicit_feedback_reference
  end subroutine get_provider_implicit_zero_mode_contract

  function get_provider_backend_name(self) result(name)
    class(matching_plane_response_provider_type), intent(in) :: self
    character(len=16) :: name

    select case (self%backend)
    case (provider_backend_table)
      name = 'table'
    case (provider_backend_zhao_online)
      name = 'zhao_online'
    case default
      name = 'none'
    end select
  end function get_provider_backend_name

  function get_provider_content_fingerprint(self) result(fingerprint)
    class(matching_plane_response_provider_type), intent(in) :: self
    character(len=16) :: fingerprint

    fingerprint = self%content_fingerprint
  end function get_provider_content_fingerprint

  logical function matching_plane_provider_is_active(self) result(active)
    class(matching_plane_response_provider_type), intent(in) :: self

    active = self%active
  end function matching_plane_provider_is_active

end module bem_matching_plane_response_provider

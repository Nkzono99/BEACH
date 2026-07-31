!> シミュレーション実行中に再利用するバッチ作業配列を管理する。
module bem_simulator_workspace
  use bem_kinds, only: dp, i32, i64
  implicit none
  private

  public :: simulator_batch_workspace_type

  type :: simulator_batch_workspace_type
    real(dp), allocatable :: dq_thread(:, :)
    real(dp), allocatable :: dq(:)
    real(dp), allocatable :: photo_emission_dq(:)
    logical, allocatable :: escaped_boundary_flag(:)
    logical, allocatable :: absorbed_flag(:)
    integer(i32), allocatable :: absorbed_element(:)
    logical, allocatable :: soft_discarded_boundary_flag(:)
    real(dp), allocatable :: q_before(:)
    real(dp), allocatable :: candidate_charge(:)
    logical :: charge_candidate_ready = .false.
    real(dp), allocatable :: neutral_return_charge_values(:)
    real(dp), allocatable :: neutral_return_emitted_charge(:)
    real(dp), allocatable :: neutral_return_absorbed_charge(:)
    real(dp), allocatable :: neutral_return_unresolved_charge(:)
    real(dp), allocatable :: neutral_return_weight_scale(:)
    real(dp), allocatable :: neutral_return_correction(:)
    real(dp), allocatable :: neutral_return_unresolved_fraction(:)
    integer(i64), allocatable :: neutral_return_terminal_counts(:)
    real(dp), allocatable :: ledger_charge_values(:)
    integer(i64), allocatable :: ledger_count_values(:)
  contains
    procedure :: init => init_simulator_batch_workspace
    procedure :: reset_before_injection
    procedure :: prepare_particle_flags
  end type simulator_batch_workspace_type

contains

  subroutine init_simulator_batch_workspace(self, nelem, nspecies, nthreads, candidate_charge_enabled)
    class(simulator_batch_workspace_type), intent(inout) :: self
    integer(i32), intent(in) :: nelem, nspecies, nthreads
    logical, intent(in), optional :: candidate_charge_enabled
    logical :: allocate_candidate

    if (nelem <= 0_i32) error stop 'simulator workspace requires nelem > 0.'
    if (nspecies <= 0_i32) error stop 'simulator workspace requires nspecies > 0.'
    if (nthreads <= 0_i32) error stop 'simulator workspace requires nthreads > 0.'
    if (allocated(self%dq_thread)) error stop 'simulator workspace is already initialized.'

    allocate (self%dq_thread(nelem, nthreads), self%dq(nelem), self%photo_emission_dq(nelem))
    allocate (self%q_before(nelem))
    allocate_candidate = .false.
    if (present(candidate_charge_enabled)) allocate_candidate = candidate_charge_enabled
    if (allocate_candidate) then
      allocate (self%candidate_charge(nelem))
    else
      allocate (self%candidate_charge(0))
    end if
    allocate (self%neutral_return_charge_values(6_i32*nspecies))
    allocate (self%neutral_return_terminal_counts(3_i32*nspecies))
    allocate ( &
      self%neutral_return_emitted_charge(nspecies), self%neutral_return_absorbed_charge(nspecies), &
      self%neutral_return_unresolved_charge(nspecies), self%neutral_return_weight_scale(nspecies), &
      self%neutral_return_correction(nspecies), self%neutral_return_unresolved_fraction(nspecies) &
      )
    allocate (self%ledger_charge_values(5_i32*nspecies), self%ledger_count_values(5_i32*nspecies))

    self%dq = 0.0_dp
    self%q_before = 0.0_dp
    self%candidate_charge = 0.0_dp
    self%ledger_charge_values = 0.0_dp
    self%ledger_count_values = 0_i64
    call self%reset_before_injection()
  end subroutine init_simulator_batch_workspace

  subroutine reset_before_injection(self)
    class(simulator_batch_workspace_type), intent(inout) :: self

    if (.not. allocated(self%dq_thread)) error stop 'simulator workspace is not initialized.'
    self%dq_thread = 0.0_dp
    self%photo_emission_dq = 0.0_dp
    self%candidate_charge = 0.0_dp
    self%charge_candidate_ready = .false.
    self%neutral_return_charge_values = 0.0_dp
    self%neutral_return_terminal_counts = 0_i64
    self%neutral_return_emitted_charge = 0.0_dp
    self%neutral_return_absorbed_charge = 0.0_dp
    self%neutral_return_unresolved_charge = 0.0_dp
    self%neutral_return_weight_scale = 1.0_dp
    self%neutral_return_correction = 0.0_dp
    self%neutral_return_unresolved_fraction = 0.0_dp
  end subroutine reset_before_injection

  subroutine prepare_particle_flags(self, particle_count)
    class(simulator_batch_workspace_type), intent(inout) :: self
    integer(i32), intent(in) :: particle_count

    if (particle_count < 0_i32) error stop 'simulator workspace particle count must be >= 0.'
    call ensure_logical_capacity(self%escaped_boundary_flag, particle_count)
    call ensure_logical_capacity(self%absorbed_flag, particle_count)
    call ensure_i32_capacity(self%absorbed_element, particle_count)
    call ensure_logical_capacity(self%soft_discarded_boundary_flag, particle_count)
    self%escaped_boundary_flag(:particle_count) = .false.
    self%absorbed_flag(:particle_count) = .false.
    self%absorbed_element(:particle_count) = -1_i32
    self%soft_discarded_boundary_flag(:particle_count) = .false.
  end subroutine prepare_particle_flags

  subroutine ensure_logical_capacity(values, required_size)
    logical, allocatable, intent(inout) :: values(:)
    integer(i32), intent(in) :: required_size
    logical, allocatable :: grown(:)

    if (allocated(values)) then
      if (size(values) >= required_size) return
    end if
    allocate (grown(required_size))
    call move_alloc(grown, values)
  end subroutine ensure_logical_capacity

  subroutine ensure_i32_capacity(values, required_size)
    integer(i32), allocatable, intent(inout) :: values(:)
    integer(i32), intent(in) :: required_size
    integer(i32), allocatable :: grown(:)

    if (allocated(values)) then
      if (size(values) >= required_size) return
    end if
    allocate (grown(required_size))
    call move_alloc(grown, values)
  end subroutine ensure_i32_capacity

end module bem_simulator_workspace

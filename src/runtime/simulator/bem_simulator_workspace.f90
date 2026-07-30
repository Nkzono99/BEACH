!> シミュレーション実行中に再利用するバッチ作業配列を管理する。
module bem_simulator_workspace
  use bem_kinds, only: dp, i32, i64
  use bem_interface_types, only: interface_crossing_type
  use bem_outer_event_queue, only: outer_event_record_type
  implicit none
  private

  public :: simulator_batch_workspace_type

  !> 1回の simulation run を通して再利用する作業領域。
  type :: simulator_batch_workspace_type
    real(dp), allocatable :: dq_thread(:, :)
    real(dp), allocatable :: dq(:)
    real(dp), allocatable :: photo_emission_dq(:)
    real(dp), allocatable :: interface_outward_thread(:, :)
    real(dp), allocatable :: interface_returned_thread(:, :)
    real(dp), allocatable :: interface_tau_max_thread(:)
    real(dp), allocatable :: interface_frozen_ratio_max_thread(:)
    real(dp), allocatable :: interface_energy_error_max_thread(:)
    logical, allocatable :: escaped_boundary_flag(:)
    logical, allocatable :: absorbed_flag(:)
    integer(i32), allocatable :: absorbed_element(:)
    logical, allocatable :: soft_discarded_boundary_flag(:)
    logical, allocatable :: queued_outer_flag(:)
    logical, allocatable :: deferred_mean_interface_flag(:)
    integer(i32), allocatable :: deferred_mean_interface_step(:)
    type(interface_crossing_type), allocatable :: deferred_mean_interface_crossing(:)
    integer(i32), allocatable :: deferred_mean_return_element(:)
    logical, allocatable :: deferred_mean_terminal_absorbed(:)
    logical, allocatable :: deferred_mean_terminal_escaped(:)
    type(outer_event_record_type), allocatable :: outer_event_staging(:)
    real(dp), allocatable :: q_before(:)
    real(dp), allocatable :: mean_pending_charge(:)
    real(dp), allocatable :: mean_deferred_source_charge(:)
    real(dp), allocatable :: mean_returned_destination_charge(:)
    real(dp), allocatable :: mean_candidate_charge(:)
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

  !> mesh/species/thread 数が一定の run-local 作業配列を一度だけ確保する。
  subroutine init_simulator_batch_workspace( &
    self, nelem, nspecies, nthreads, implicit_mean_enabled, candidate_charge_enabled &
    )
    class(simulator_batch_workspace_type), intent(inout) :: self
    integer(i32), intent(in) :: nelem, nspecies, nthreads
    logical, intent(in), optional :: implicit_mean_enabled
    logical, intent(in), optional :: candidate_charge_enabled
    logical :: allocate_mean_storage, allocate_candidate_storage

    if (nelem <= 0_i32) error stop 'simulator workspace requires nelem > 0.'
    if (nspecies <= 0_i32) error stop 'simulator workspace requires nspecies > 0.'
    if (nthreads <= 0_i32) error stop 'simulator workspace requires nthreads > 0.'
    if (allocated(self%dq_thread)) error stop 'simulator workspace is already initialized.'

    allocate (self%dq_thread(nelem, nthreads), self%dq(nelem), self%photo_emission_dq(nelem))
    allocate (self%interface_outward_thread(nspecies, nthreads))
    allocate (self%interface_returned_thread(nspecies, nthreads))
    allocate (self%interface_tau_max_thread(nthreads), self%interface_frozen_ratio_max_thread(nthreads))
    allocate (self%interface_energy_error_max_thread(nthreads))
    allocate (self%q_before(nelem))
    allocate (self%neutral_return_charge_values(6_i32*nspecies))
    allocate (self%neutral_return_terminal_counts(3_i32*nspecies))
    allocate ( &
      self%neutral_return_emitted_charge(nspecies), self%neutral_return_absorbed_charge(nspecies), &
      self%neutral_return_unresolved_charge(nspecies), self%neutral_return_weight_scale(nspecies), &
      self%neutral_return_correction(nspecies), self%neutral_return_unresolved_fraction(nspecies) &
      )
    allocate_mean_storage = .false.
    if (present(implicit_mean_enabled)) allocate_mean_storage = implicit_mean_enabled
    allocate_candidate_storage = allocate_mean_storage
    if (present(candidate_charge_enabled)) then
      allocate_candidate_storage = allocate_candidate_storage .or. candidate_charge_enabled
    end if
    if (allocate_candidate_storage) then
      allocate (self%mean_candidate_charge(nelem))
    else
      allocate (self%mean_candidate_charge(0))
    end if
    if (allocate_mean_storage) then
      allocate ( &
        self%mean_pending_charge(nelem), self%mean_deferred_source_charge(nelem), &
        self%mean_returned_destination_charge(nelem) &
        )
    else
      allocate ( &
        self%mean_pending_charge(0), self%mean_deferred_source_charge(0), &
        self%mean_returned_destination_charge(0) &
        )
    end if
    allocate (self%ledger_charge_values(7_i32*nspecies), self%ledger_count_values(5_i32*nspecies))

    self%dq = 0.0_dp
    self%q_before = 0.0_dp
    self%mean_pending_charge = 0.0_dp
    self%mean_deferred_source_charge = 0.0_dp
    self%mean_returned_destination_charge = 0.0_dp
    self%mean_candidate_charge = 0.0_dp
    self%charge_candidate_ready = .false.
    self%neutral_return_charge_values = 0.0_dp
    self%neutral_return_terminal_counts = 0_i64
    self%ledger_charge_values = 0.0_dp
    self%ledger_count_values = 0_i64
    call self%reset_before_injection()
  end subroutine init_simulator_batch_workspace

  !> 粒子生成前に、前バッチから持ち越してはいけない集計値を初期化する。
  subroutine reset_before_injection(self)
    class(simulator_batch_workspace_type), intent(inout) :: self

    if (.not. allocated(self%dq_thread)) error stop 'simulator workspace is not initialized.'
    self%dq_thread = 0.0_dp
    self%photo_emission_dq = 0.0_dp
    self%interface_outward_thread = 0.0_dp
    self%interface_returned_thread = 0.0_dp
    self%interface_tau_max_thread = 0.0_dp
    self%interface_frozen_ratio_max_thread = 0.0_dp
    self%interface_energy_error_max_thread = 0.0_dp
    self%mean_pending_charge = 0.0_dp
    self%mean_deferred_source_charge = 0.0_dp
    self%mean_returned_destination_charge = 0.0_dp
    self%mean_candidate_charge = 0.0_dp
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

  !> 可変粒子数に合わせて outcome flag を grow-only で確保し、有効範囲を初期化する。
  subroutine prepare_particle_flags(self, particle_count, outer_queue_enabled, implicit_mean_enabled)
    class(simulator_batch_workspace_type), intent(inout) :: self
    integer(i32), intent(in) :: particle_count
    logical, intent(in), optional :: outer_queue_enabled, implicit_mean_enabled
    logical :: prepare_outer_staging, prepare_mean_staging

    if (particle_count < 0_i32) error stop 'simulator workspace particle count must be >= 0.'
    call ensure_logical_capacity(self%escaped_boundary_flag, particle_count)
    call ensure_logical_capacity(self%absorbed_flag, particle_count)
    call ensure_i32_capacity(self%absorbed_element, particle_count)
    call ensure_logical_capacity(self%soft_discarded_boundary_flag, particle_count)
    call ensure_logical_capacity(self%queued_outer_flag, particle_count)
    prepare_mean_staging = .false.
    if (present(implicit_mean_enabled)) prepare_mean_staging = implicit_mean_enabled
    if (prepare_mean_staging) then
      call ensure_logical_capacity(self%deferred_mean_interface_flag, particle_count)
      call ensure_i32_capacity(self%deferred_mean_interface_step, particle_count)
      call ensure_interface_crossing_capacity(self%deferred_mean_interface_crossing, particle_count)
      call ensure_i32_capacity(self%deferred_mean_return_element, particle_count)
      call ensure_logical_capacity(self%deferred_mean_terminal_absorbed, particle_count)
      call ensure_logical_capacity(self%deferred_mean_terminal_escaped, particle_count)
    else
      if (.not. allocated(self%deferred_mean_interface_flag)) allocate (self%deferred_mean_interface_flag(0))
      if (.not. allocated(self%deferred_mean_interface_step)) allocate (self%deferred_mean_interface_step(0))
      if (.not. allocated(self%deferred_mean_interface_crossing)) allocate (self%deferred_mean_interface_crossing(0))
      if (.not. allocated(self%deferred_mean_return_element)) allocate (self%deferred_mean_return_element(0))
      if (.not. allocated(self%deferred_mean_terminal_absorbed)) allocate (self%deferred_mean_terminal_absorbed(0))
      if (.not. allocated(self%deferred_mean_terminal_escaped)) allocate (self%deferred_mean_terminal_escaped(0))
    end if
    prepare_outer_staging = .false.
    if (present(outer_queue_enabled)) prepare_outer_staging = outer_queue_enabled
    if (prepare_outer_staging) then
      call ensure_outer_event_capacity(self%outer_event_staging, particle_count)
    else if (.not. allocated(self%outer_event_staging)) then
      allocate (self%outer_event_staging(0))
    end if
    self%escaped_boundary_flag(:particle_count) = .false.
    self%absorbed_flag(:particle_count) = .false.
    self%absorbed_element(:particle_count) = -1_i32
    self%soft_discarded_boundary_flag(:particle_count) = .false.
    self%queued_outer_flag(:particle_count) = .false.
    if (prepare_mean_staging) then
      self%deferred_mean_interface_flag(:particle_count) = .false.
      self%deferred_mean_interface_step(:particle_count) = 0_i32
      self%deferred_mean_interface_crossing(:particle_count) = interface_crossing_type()
      self%deferred_mean_return_element(:particle_count) = -1_i32
      self%deferred_mean_terminal_absorbed(:particle_count) = .false.
      self%deferred_mean_terminal_escaped(:particle_count) = .false.
    end if
    if (prepare_outer_staging) self%outer_event_staging(:particle_count) = outer_event_record_type()
  end subroutine prepare_particle_flags

  !> logical 配列を既存容量を保ちながら必要時だけ拡張する。
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

  !> integer staging 配列を既存容量を保ちながら必要時だけ拡張する。
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

  !> interface crossing staging 配列を既存容量を保ちながら必要時だけ拡張する。
  subroutine ensure_interface_crossing_capacity(values, required_size)
    type(interface_crossing_type), allocatable, intent(inout) :: values(:)
    integer(i32), intent(in) :: required_size
    type(interface_crossing_type), allocatable :: grown(:)

    if (allocated(values)) then
      if (size(values) >= required_size) return
    end if
    allocate (grown(required_size))
    call move_alloc(grown, values)
  end subroutine ensure_interface_crossing_capacity

  !> queue event staging 配列を粒子数に合わせて grow-only で確保する。
  subroutine ensure_outer_event_capacity(values, required_size)
    type(outer_event_record_type), allocatable, intent(inout) :: values(:)
    integer(i32), intent(in) :: required_size
    type(outer_event_record_type), allocatable :: grown(:)

    if (allocated(values)) then
      if (size(values) >= required_size) return
    end if
    allocate (grown(required_size))
    call move_alloc(grown, values)
  end subroutine ensure_outer_event_capacity

end module bem_simulator_workspace

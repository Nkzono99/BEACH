!> シミュレーション実行中に再利用するバッチ作業配列を管理する。
module bem_simulator_workspace
  use bem_kinds, only: dp, i32, i64
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
    logical, allocatable :: soft_discarded_boundary_flag(:)
    real(dp), allocatable :: q_before(:)
    real(dp), allocatable :: ledger_charge_values(:)
    integer(i64), allocatable :: ledger_count_values(:)
  contains
    procedure :: init => init_simulator_batch_workspace
    procedure :: reset_before_injection
    procedure :: prepare_particle_flags
  end type simulator_batch_workspace_type

contains

  !> mesh/species/thread 数が一定の run-local 作業配列を一度だけ確保する。
  subroutine init_simulator_batch_workspace(self, nelem, nspecies, nthreads)
    class(simulator_batch_workspace_type), intent(inout) :: self
    integer(i32), intent(in) :: nelem, nspecies, nthreads

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
    allocate (self%ledger_charge_values(7_i32*nspecies), self%ledger_count_values(5_i32*nspecies))

    self%dq = 0.0_dp
    self%q_before = 0.0_dp
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
  end subroutine reset_before_injection

  !> 可変粒子数に合わせて outcome flag を grow-only で確保し、有効範囲を初期化する。
  subroutine prepare_particle_flags(self, particle_count)
    class(simulator_batch_workspace_type), intent(inout) :: self
    integer(i32), intent(in) :: particle_count

    if (particle_count < 0_i32) error stop 'simulator workspace particle count must be >= 0.'
    call ensure_logical_capacity(self%escaped_boundary_flag, particle_count)
    call ensure_logical_capacity(self%absorbed_flag, particle_count)
    call ensure_logical_capacity(self%soft_discarded_boundary_flag, particle_count)
    self%escaped_boundary_flag(:particle_count) = .false.
    self%absorbed_flag(:particle_count) = .false.
    self%soft_discarded_boundary_flag(:particle_count) = .false.
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

end module bem_simulator_workspace

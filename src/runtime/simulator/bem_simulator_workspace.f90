!> シミュレーション実行中に再利用するバッチ作業配列を管理する。
module bem_simulator_workspace
  use bem_kinds, only: dp, i32, i64
  use bem_outflow_spectrum, only: outflow_spectrum_bins
  implicit none
  private

  public :: simulator_batch_workspace_type

  !> fixed_current_charge_values(species, channel) の channel。旧 1 次元配列の (channel-1)*n+species と同じ並び。
  integer(i32), parameter, public :: fixed_charge_absorbed = 1_i32
  integer(i32), parameter, public :: fixed_charge_emission = 2_i32
  integer(i32), parameter, public :: fixed_charge_escaped = 3_i32
  integer(i32), parameter, public :: fixed_charge_barrier_return = 4_i32
  integer(i32), parameter :: fixed_charge_channels = 4_i32
  !> neutral_return_charge_values(species, channel) の channel。
  integer(i32), parameter, public :: neutral_charge_emitted = 1_i32
  integer(i32), parameter, public :: neutral_charge_absorbed = 2_i32
  integer(i32), parameter, public :: neutral_charge_unresolved = 3_i32
  integer(i32), parameter, public :: neutral_charge_escaped = 4_i32
  integer(i32), parameter, public :: neutral_charge_soft = 5_i32
  integer(i32), parameter, public :: neutral_charge_invalid = 6_i32
  integer(i32), parameter :: neutral_charge_channels = 6_i32
  !> neutral_return_terminal_counts(species, kind) の kind。
  integer(i32), parameter, public :: neutral_count_escaped = 1_i32
  integer(i32), parameter, public :: neutral_count_soft = 2_i32
  integer(i32), parameter, public :: neutral_count_invalid = 3_i32
  integer(i32), parameter :: neutral_count_kinds = 3_i32

  type :: simulator_batch_workspace_type
    real(dp), allocatable :: dq_thread(:, :)
    real(dp), allocatable :: dq(:)
    real(dp), allocatable :: photo_emission_dq(:, :)
    logical, allocatable :: escaped_boundary_flag(:)
    logical, allocatable :: absorbed_flag(:)
    integer(i32), allocatable :: absorbed_element(:)
    logical, allocatable :: soft_discarded_boundary_flag(:)
    real(dp), allocatable :: q_before(:)
    real(dp), allocatable :: candidate_charge(:)
    logical :: charge_candidate_ready = .false.
    !> neutral_return が rank 合計する電荷（species, neutral_charge_*）。
    real(dp), allocatable :: neutral_return_charge_values(:, :)
    real(dp), allocatable :: neutral_return_emitted_charge(:)
    real(dp), allocatable :: neutral_return_absorbed_charge(:)
    real(dp), allocatable :: neutral_return_unresolved_charge(:)
    real(dp), allocatable :: neutral_return_weight_scale(:)
    real(dp), allocatable :: neutral_return_correction(:)
    real(dp), allocatable :: neutral_return_unresolved_fraction(:)
    !> 固定電流 closure が rank 合計する raw 電荷（species, fixed_charge_*）。吸収、放出反作用、脱出、
    !! 外部障壁で戻った電荷（H 外向き通過のうち戻った分、帯電符号）。
    real(dp), allocatable :: fixed_current_charge_values(:, :)
    real(dp), allocatable :: fixed_absorbed_target_charge(:)
    real(dp), allocatable :: fixed_absorbed_weight_scale(:)
    real(dp), allocatable :: fixed_emission_target_charge(:)
    real(dp), allocatable :: fixed_emission_weight_scale(:)
    real(dp), allocatable :: fixed_escape_target_charge(:)
    real(dp), allocatable :: fixed_escape_correction(:)
    real(dp), allocatable :: fixed_current_correction(:)
    !> 脱出目標の調整で収支に組み込んだ電荷と、組み込めなかった電荷（どちらも帯電符号の大きさ）。
    !! partner を持つ種では差し引いた partner の脱出とその超過、放出目標が追跡した源の種では
    !! 収支に載せた根の脱出と、放出がなく載せられなかった根の脱出。
    real(dp), allocatable :: fixed_escape_closed_charge(:)
    real(dp), allocatable :: fixed_escape_unclosed_charge(:)
    !> この batch で固定電流 target を適用した channel。1=吸収、2=放出、3=脱出。
    logical, allocatable :: fixed_channel_active(:, :)
    !> neutral_return の非対応な終端の数（species, neutral_count_*）。
    integer(i64), allocatable :: neutral_return_terminal_counts(:, :)
    real(dp), allocatable :: ledger_charge_values(:)
    integer(i64), allocatable :: ledger_count_values(:)
    ! (outward number, outward normal energy [J], barrier return number, barrier escape number,
    !  outward number per normal-energy bin of bem_outflow_spectrum)
    real(dp), allocatable :: matching_plane_moments_thread(:, :, :)
    ! hits, emitted macro count, generation-capped count; charge [C], sum of yields.
    integer(i64), allocatable :: secondary_counts_thread(:, :, :)
    real(dp), allocatable :: secondary_values_thread(:, :, :)
  contains
    procedure :: init => init_simulator_batch_workspace
    procedure :: reset_before_injection
    procedure :: prepare_particle_flags
    procedure :: extend_particle_flags
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

    allocate (self%dq_thread(nelem, nthreads), self%dq(nelem), self%photo_emission_dq(nelem, nspecies))
    allocate (self%q_before(nelem))
    allocate_candidate = .false.
    if (present(candidate_charge_enabled)) allocate_candidate = candidate_charge_enabled
    if (allocate_candidate) then
      allocate (self%candidate_charge(nelem))
    else
      allocate (self%candidate_charge(0))
    end if
    allocate (self%neutral_return_charge_values(nspecies, neutral_charge_channels))
    allocate (self%neutral_return_terminal_counts(nspecies, neutral_count_kinds))
    allocate ( &
      self%neutral_return_emitted_charge(nspecies), self%neutral_return_absorbed_charge(nspecies), &
      self%neutral_return_unresolved_charge(nspecies), self%neutral_return_weight_scale(nspecies), &
      self%neutral_return_correction(nspecies), self%neutral_return_unresolved_fraction(nspecies) &
      )
    allocate (self%fixed_current_charge_values(nspecies, fixed_charge_channels))
    allocate ( &
      self%fixed_absorbed_target_charge(nspecies), self%fixed_absorbed_weight_scale(nspecies), &
      self%fixed_emission_target_charge(nspecies), self%fixed_emission_weight_scale(nspecies), &
      self%fixed_escape_target_charge(nspecies), &
      self%fixed_escape_correction(nspecies), &
      self%fixed_current_correction(nspecies), &
      self%fixed_escape_closed_charge(nspecies), self%fixed_escape_unclosed_charge(nspecies) &
      )
    allocate (self%fixed_channel_active(3_i32, nspecies))
    allocate (self%ledger_charge_values(5_i32*nspecies), self%ledger_count_values(5_i32*nspecies))
    ! 1:4 は H 通過モーメント、5 以降は外向き法線運動エネルギーの bin 別重み。
    allocate (self%matching_plane_moments_thread(4_i32 + outflow_spectrum_bins, nspecies, nthreads))

    allocate (self%secondary_counts_thread(3, nspecies, nthreads), self%secondary_values_thread(2, nspecies, nthreads))
    self%dq = 0.0_dp
    self%q_before = 0.0_dp
    self%candidate_charge = 0.0_dp
    self%ledger_charge_values = 0.0_dp
    self%ledger_count_values = 0_i64
    self%matching_plane_moments_thread = 0.0_dp
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
    self%fixed_current_charge_values = 0.0_dp
    self%fixed_absorbed_target_charge = 0.0_dp
    self%fixed_absorbed_weight_scale = 1.0_dp
    self%fixed_emission_target_charge = 0.0_dp
    self%fixed_emission_weight_scale = 1.0_dp
    self%fixed_escape_target_charge = 0.0_dp
    self%fixed_escape_correction = 0.0_dp
    self%fixed_current_correction = 0.0_dp
    self%fixed_escape_closed_charge = 0.0_dp
    self%fixed_escape_unclosed_charge = 0.0_dp
    self%fixed_channel_active = .false.
    self%matching_plane_moments_thread = 0.0_dp
    self%secondary_counts_thread = 0_i64
    self%secondary_values_thread = 0.0_dp
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

  !> Preserve terminal flags of earlier generations when appending secondary particles.
  subroutine extend_particle_flags(self, old_count, new_count)
    class(simulator_batch_workspace_type), intent(inout) :: self
    integer(i32), intent(in) :: old_count, new_count
    call ensure_logical_capacity(self%escaped_boundary_flag, new_count)
    call ensure_logical_capacity(self%absorbed_flag, new_count)
    call ensure_i32_capacity(self%absorbed_element, new_count)
    call ensure_logical_capacity(self%soft_discarded_boundary_flag, new_count)
    self%escaped_boundary_flag(old_count + 1:new_count) = .false.
    self%absorbed_flag(old_count + 1:new_count) = .false.
    self%absorbed_element(old_count + 1:new_count) = -1_i32
    self%soft_discarded_boundary_flag(old_count + 1:new_count) = .false.
  end subroutine extend_particle_flags

  subroutine ensure_logical_capacity(values, required_size)
    logical, allocatable, intent(inout) :: values(:)
    integer(i32), intent(in) :: required_size
    logical, allocatable :: grown(:)

    if (allocated(values)) then
      if (size(values) >= required_size) return
    end if
    allocate (grown(required_size))
    if (allocated(values)) grown(:size(values)) = values
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
    if (allocated(values)) grown(:size(values)) = values
    call move_alloc(grown, values)
  end subroutine ensure_i32_capacity

end module bem_simulator_workspace

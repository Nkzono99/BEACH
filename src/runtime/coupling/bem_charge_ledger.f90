!> batch 間の signed charge stock と移送 flux から電荷収支を集計する。
module bem_charge_ledger
  use bem_kinds, only: dp, i32, i64
  implicit none
  private

  type, public :: charge_ledger_type
    integer(i32) :: nspecies = 0_i32
    integer(i32) :: batch_count = 0_i32
    real(dp) :: surface_charge_before = 0.0_dp
    real(dp) :: surface_charge_after = 0.0_dp
    real(dp) :: local_flight_charge_before = 0.0_dp
    real(dp) :: local_flight_charge_after = 0.0_dp
    real(dp) :: unresolved_stock_before = 0.0_dp
    real(dp) :: unresolved_stock_after = 0.0_dp
    real(dp), allocatable :: injected_from_remote(:)
    real(dp), allocatable :: emitted_from_surface(:)
    real(dp), allocatable :: absorbed_on_surface(:)
    real(dp), allocatable :: escaped_to_infinity(:)
    real(dp), allocatable :: discarded_unresolved(:)
    real(dp), allocatable :: neutral_return_correction(:)
    real(dp), allocatable :: neutral_return_weight_scale(:)
    real(dp), allocatable :: neutral_return_unresolved_fraction(:)
    integer(i64), allocatable :: injected_count(:)
    integer(i64), allocatable :: emitted_count(:)
    integer(i64), allocatable :: absorbed_count(:)
    integer(i64), allocatable :: escaped_count(:)
    integer(i64), allocatable :: discarded_unresolved_count(:)
  contains
    procedure :: init => init_charge_ledger
    procedure :: reset => reset_charge_ledger
    procedure :: residual => charge_ledger_residual
    procedure :: discarded_unresolved_abs
    procedure :: has_unresolved_discard
  end type charge_ledger_type

  public :: accumulate_charge_ledger

contains

  !> species 数に合わせて ledger の flux/count 配列を確保する。
  subroutine init_charge_ledger(self, nspecies)
    class(charge_ledger_type), intent(inout) :: self
    integer(i32), intent(in) :: nspecies

    if (nspecies < 1_i32) error stop 'charge ledger nspecies must be >= 1.'
    if (allocated(self%injected_from_remote)) then
      if (size(self%injected_from_remote) /= nspecies) call release_charge_ledger_arrays(self)
    end if
    if (.not. allocated(self%injected_from_remote)) then
      allocate ( &
        self%injected_from_remote(nspecies), self%emitted_from_surface(nspecies), &
        self%absorbed_on_surface(nspecies), self%escaped_to_infinity(nspecies), &
        self%discarded_unresolved(nspecies), self%neutral_return_correction(nspecies), &
        self%neutral_return_weight_scale(nspecies), self%neutral_return_unresolved_fraction(nspecies), &
        self%injected_count(nspecies), &
        self%emitted_count(nspecies), self%absorbed_count(nspecies), self%escaped_count(nspecies), &
        self%discarded_unresolved_count(nspecies) &
        )
    end if
    self%nspecies = nspecies
    call self%reset(0_i32)
  end subroutine init_charge_ledger

  !> stock と flux を 0 に戻し、対象 batch 番号を設定する。
  subroutine reset_charge_ledger(self, batch_count)
    class(charge_ledger_type), intent(inout) :: self
    integer(i32), intent(in) :: batch_count

    if (.not. allocated(self%injected_from_remote)) error stop 'charge ledger must be initialized before reset.'
    if (batch_count < 0_i32) error stop 'charge ledger batch_count must be >= 0.'
    self%batch_count = batch_count
    self%surface_charge_before = 0.0_dp
    self%surface_charge_after = 0.0_dp
    self%local_flight_charge_before = 0.0_dp
    self%local_flight_charge_after = 0.0_dp
    self%unresolved_stock_before = 0.0_dp
    self%unresolved_stock_after = 0.0_dp
    self%injected_from_remote = 0.0_dp
    self%emitted_from_surface = 0.0_dp
    self%absorbed_on_surface = 0.0_dp
    self%escaped_to_infinity = 0.0_dp
    self%discarded_unresolved = 0.0_dp
    self%neutral_return_correction = 0.0_dp
    self%neutral_return_weight_scale = 1.0_dp
    self%neutral_return_unresolved_fraction = 0.0_dp
    self%injected_count = 0_i64
    self%emitted_count = 0_i64
    self%absorbed_count = 0_i64
    self%escaped_count = 0_i64
    self%discarded_unresolved_count = 0_i64
  end subroutine reset_charge_ledger

  !> stock 差分と remote flux から transactional conservation residual を返す。
  pure real(dp) function charge_ledger_residual(self) result(residual)
    class(charge_ledger_type), intent(in) :: self

    residual = (self%surface_charge_after - self%surface_charge_before) + &
               (self%local_flight_charge_after - self%local_flight_charge_before) + &
               (self%unresolved_stock_after - self%unresolved_stock_before) - &
               sum(self%injected_from_remote) + sum(self%escaped_to_infinity) + &
               sum(self%discarded_unresolved) - sum(self%neutral_return_correction)
  end function charge_ledger_residual

  !> species 間で相殺しない max-step discard charge の絶対値和を返す。
  pure real(dp) function discarded_unresolved_abs(self) result(total_abs)
    class(charge_ledger_type), intent(in) :: self

    total_abs = sum(abs(self%discarded_unresolved))
  end function discarded_unresolved_abs

  !> residual とは独立に discard charge/count が許容値を超えたか返す。
  pure logical function has_unresolved_discard(self, max_abs_charge, max_count) result(exceeded)
    class(charge_ledger_type), intent(in) :: self
    real(dp), intent(in) :: max_abs_charge
    integer(i64), intent(in) :: max_count

    exceeded = self%discarded_unresolved_abs() > max(0.0_dp, max_abs_charge) .or. &
               sum(self%discarded_unresolved_count) > max(0_i64, max_count)
  end function has_unresolved_discard

  !> 連続する batch ledger を、初期 stock と最終 stock を保って累積する。
  subroutine accumulate_charge_ledger(cumulative, batch)
    type(charge_ledger_type), intent(inout) :: cumulative
    type(charge_ledger_type), intent(in) :: batch
    logical :: first_batch
    integer(i32) :: species_idx

    if (batch%nspecies < 1_i32 .or. .not. allocated(batch%injected_from_remote)) then
      error stop 'cannot accumulate an uninitialized charge ledger.'
    end if
    first_batch = .not. allocated(cumulative%injected_from_remote) .or. cumulative%batch_count == 0_i32
    if (.not. allocated(cumulative%injected_from_remote)) then
      call cumulative%init(batch%nspecies)
    else if (cumulative%nspecies /= batch%nspecies) then
      error stop 'charge ledger species count mismatch during accumulation.'
    end if
    if (first_batch) then
      cumulative%surface_charge_before = batch%surface_charge_before
      cumulative%local_flight_charge_before = batch%local_flight_charge_before
      cumulative%unresolved_stock_before = batch%unresolved_stock_before
    end if
    cumulative%surface_charge_after = batch%surface_charge_after
    cumulative%local_flight_charge_after = batch%local_flight_charge_after
    cumulative%unresolved_stock_after = batch%unresolved_stock_after
    cumulative%batch_count = max(cumulative%batch_count, batch%batch_count)
    cumulative%injected_from_remote = cumulative%injected_from_remote + batch%injected_from_remote
    cumulative%emitted_from_surface = cumulative%emitted_from_surface + batch%emitted_from_surface
    cumulative%absorbed_on_surface = cumulative%absorbed_on_surface + batch%absorbed_on_surface
    cumulative%escaped_to_infinity = cumulative%escaped_to_infinity + batch%escaped_to_infinity
    cumulative%discarded_unresolved = cumulative%discarded_unresolved + batch%discarded_unresolved
    cumulative%neutral_return_correction = &
      cumulative%neutral_return_correction + batch%neutral_return_correction
    cumulative%neutral_return_weight_scale = 1.0_dp
    cumulative%neutral_return_unresolved_fraction = 0.0_dp
    do species_idx = 1_i32, cumulative%nspecies
      if (cumulative%neutral_return_correction(species_idx) == 0.0_dp) cycle
      if (cumulative%emitted_from_surface(species_idx) >= 0.0_dp .or. &
          cumulative%absorbed_on_surface(species_idx) >= 0.0_dp) then
        error stop 'cumulative neutral-return diagnostics require negative emitted and absorbed charge.'
      end if
      cumulative%neutral_return_weight_scale(species_idx) = &
        cumulative%emitted_from_surface(species_idx)/cumulative%absorbed_on_surface(species_idx)
      cumulative%neutral_return_unresolved_fraction(species_idx) = &
        cumulative%discarded_unresolved(species_idx)/cumulative%emitted_from_surface(species_idx)
    end do
    cumulative%injected_count = cumulative%injected_count + batch%injected_count
    cumulative%emitted_count = cumulative%emitted_count + batch%emitted_count
    cumulative%absorbed_count = cumulative%absorbed_count + batch%absorbed_count
    cumulative%escaped_count = cumulative%escaped_count + batch%escaped_count
    cumulative%discarded_unresolved_count = &
      cumulative%discarded_unresolved_count + batch%discarded_unresolved_count
  end subroutine accumulate_charge_ledger

  subroutine release_charge_ledger_arrays(self)
    type(charge_ledger_type), intent(inout) :: self

    deallocate ( &
      self%injected_from_remote, self%emitted_from_surface, self%absorbed_on_surface, &
      self%escaped_to_infinity, self%discarded_unresolved, self%neutral_return_correction, &
      self%neutral_return_weight_scale, self%neutral_return_unresolved_fraction, &
      self%injected_count, self%emitted_count, self%absorbed_count, &
      self%escaped_count, self%discarded_unresolved_count &
      )
  end subroutine release_charge_ledger_arrays

end module bem_charge_ledger

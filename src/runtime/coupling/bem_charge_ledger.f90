!> batch 間の signed charge stock と移送 flux から電荷収支を集計する。
module bem_charge_ledger
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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
    real(dp), allocatable :: fixed_absorbed_target_charge(:)
    real(dp), allocatable :: fixed_absorbed_weight_scale(:)
    real(dp), allocatable :: fixed_emission_target_charge(:)
    real(dp), allocatable :: fixed_emission_weight_scale(:)
    real(dp), allocatable :: fixed_escape_target_charge(:)
    real(dp), allocatable :: fixed_escape_correction(:)
    real(dp), allocatable :: fixed_current_correction(:)
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
  public :: checked_accumulate_charge
  public :: finite_charge_sum

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
        self%fixed_absorbed_target_charge(nspecies), self%fixed_absorbed_weight_scale(nspecies), &
        self%fixed_emission_target_charge(nspecies), self%fixed_emission_weight_scale(nspecies), &
        self%fixed_escape_target_charge(nspecies), &
        self%fixed_escape_correction(nspecies), &
        self%fixed_current_correction(nspecies), &
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
    self%fixed_absorbed_target_charge = 0.0_dp
    self%fixed_absorbed_weight_scale = 1.0_dp
    self%fixed_emission_target_charge = 0.0_dp
    self%fixed_emission_weight_scale = 1.0_dp
    self%fixed_escape_target_charge = 0.0_dp
    self%fixed_escape_correction = 0.0_dp
    self%fixed_current_correction = 0.0_dp
    self%injected_count = 0_i64
    self%emitted_count = 0_i64
    self%absorbed_count = 0_i64
    self%escaped_count = 0_i64
    self%discarded_unresolved_count = 0_i64
  end subroutine reset_charge_ledger

  !> stock 差分と remote flux から transactional conservation residual を返す。
  pure real(dp) function charge_ledger_residual(self) result(residual)
    class(charge_ledger_type), intent(in) :: self
    real(dp), allocatable :: terms(:)
    integer :: n, offset

    n = self%nspecies
    allocate (terms(6 + 5*n))
    terms(1:6) = [ &
                 self%surface_charge_after, -self%surface_charge_before, &
                 self%local_flight_charge_after, -self%local_flight_charge_before, &
                 self%unresolved_stock_after, -self%unresolved_stock_before &
                 ]
    offset = 6
    terms(offset + 1:offset + n) = -self%injected_from_remote
    offset = offset + n
    terms(offset + 1:offset + n) = self%escaped_to_infinity
    offset = offset + n
    terms(offset + 1:offset + n) = self%discarded_unresolved
    offset = offset + n
    terms(offset + 1:offset + n) = -self%neutral_return_correction
    offset = offset + n
    terms(offset + 1:offset + n) = -self%fixed_current_correction
    residual = finite_charge_sum(terms, 'charge ledger residual')
  end function charge_ledger_residual

  !> species 間で相殺しない max-step discard charge の絶対値和を返す。
  pure real(dp) function discarded_unresolved_abs(self) result(total_abs)
    class(charge_ledger_type), intent(in) :: self

    total_abs = finite_charge_sum(abs(self%discarded_unresolved), 'discarded unresolved absolute charge')
  end function discarded_unresolved_abs

  !> residual とは独立に discard charge/count が許容値を超えたか返す。
  pure logical function has_unresolved_discard(self, max_abs_charge, max_count) result(exceeded)
    class(charge_ledger_type), intent(in) :: self
    real(dp), intent(in) :: max_abs_charge
    integer(i64), intent(in) :: max_count

    exceeded = self%discarded_unresolved_abs() > max(0.0_dp, max_abs_charge) .or. &
               checked_nonnegative_count_sum(self%discarded_unresolved_count, 'discarded unresolved count') > &
               max(0_i64, max_count)
  end function has_unresolved_discard

  !> 連続する batch ledger を、初期 stock と最終 stock を保って累積する。
  subroutine accumulate_charge_ledger(cumulative, batch)
    type(charge_ledger_type), intent(inout) :: cumulative
    type(charge_ledger_type), intent(in) :: batch
    type(charge_ledger_type) :: candidate
    logical :: first_batch, absorbed_fixed_active, emission_fixed_active
    integer(i32) :: species_idx

    if (batch%nspecies < 1_i32 .or. .not. allocated(batch%injected_from_remote)) then
      error stop 'cannot accumulate an uninitialized charge ledger.'
    end if
    first_batch = .not. allocated(cumulative%injected_from_remote) .or. cumulative%batch_count == 0_i32
    candidate = cumulative
    if (.not. allocated(candidate%injected_from_remote)) then
      call candidate%init(batch%nspecies)
    else if (candidate%nspecies /= batch%nspecies) then
      error stop 'charge ledger species count mismatch during accumulation.'
    end if
    call validate_ledger_finite(batch, 'batch charge ledger')
    call validate_ledger_finite(candidate, 'cumulative charge ledger')
    if (first_batch) then
      candidate%surface_charge_before = batch%surface_charge_before
      candidate%local_flight_charge_before = batch%local_flight_charge_before
      candidate%unresolved_stock_before = batch%unresolved_stock_before
    end if
    candidate%surface_charge_after = batch%surface_charge_after
    candidate%local_flight_charge_after = batch%local_flight_charge_after
    candidate%unresolved_stock_after = batch%unresolved_stock_after
    candidate%batch_count = max(candidate%batch_count, batch%batch_count)
    call checked_accumulate_charge_array(candidate%injected_from_remote, batch%injected_from_remote, 'injected charge')
    call checked_accumulate_charge_array(candidate%emitted_from_surface, batch%emitted_from_surface, 'emitted charge')
    call checked_accumulate_charge_array(candidate%absorbed_on_surface, batch%absorbed_on_surface, 'absorbed charge')
    call checked_accumulate_charge_array(candidate%escaped_to_infinity, batch%escaped_to_infinity, 'escaped charge')
    call checked_accumulate_charge_array(candidate%discarded_unresolved, batch%discarded_unresolved, 'discarded charge')
    call checked_accumulate_charge_array( &
      candidate%neutral_return_correction, batch%neutral_return_correction, 'neutral-return correction' &
      )
    call checked_accumulate_charge_array( &
      candidate%fixed_absorbed_target_charge, batch%fixed_absorbed_target_charge, 'fixed absorbed target charge' &
      )
    call checked_accumulate_charge_array( &
      candidate%fixed_emission_target_charge, batch%fixed_emission_target_charge, 'fixed emission target charge' &
      )
    call checked_accumulate_charge_array( &
      candidate%fixed_escape_target_charge, batch%fixed_escape_target_charge, 'fixed escape target charge' &
      )
    call checked_accumulate_charge_array( &
      candidate%fixed_escape_correction, batch%fixed_escape_correction, 'fixed escape correction' &
      )
    call checked_accumulate_charge_array( &
      candidate%fixed_current_correction, batch%fixed_current_correction, 'fixed current correction' &
      )
    do species_idx = 1_i32, candidate%nspecies
      absorbed_fixed_active = candidate%fixed_absorbed_target_charge(species_idx) /= 0.0_dp .or. &
                              candidate%fixed_absorbed_weight_scale(species_idx) /= 1.0_dp .or. &
                              batch%fixed_absorbed_weight_scale(species_idx) /= 1.0_dp
      emission_fixed_active = candidate%fixed_emission_target_charge(species_idx) /= 0.0_dp .or. &
                              candidate%fixed_emission_weight_scale(species_idx) /= 1.0_dp .or. &
                              batch%fixed_emission_weight_scale(species_idx) /= 1.0_dp
      candidate%fixed_absorbed_weight_scale(species_idx) = 1.0_dp
      candidate%fixed_emission_weight_scale(species_idx) = 1.0_dp
      if (absorbed_fixed_active) then
        if (candidate%absorbed_on_surface(species_idx) /= 0.0_dp) then
          candidate%fixed_absorbed_weight_scale(species_idx) = checked_charge_ratio( &
                                                               candidate%fixed_absorbed_target_charge(species_idx), &
                                                               candidate%absorbed_on_surface(species_idx), &
                                                               'cumulative fixed absorbed weight scale' &
                                                               )
        end if
      end if
      if (emission_fixed_active) then
        if (candidate%emitted_from_surface(species_idx) /= 0.0_dp) then
          candidate%fixed_emission_weight_scale(species_idx) = checked_charge_ratio( &
                                                               -candidate%fixed_emission_target_charge(species_idx), &
                                                               candidate%emitted_from_surface(species_idx), &
                                                               'cumulative fixed emission weight scale' &
                                                               )
        end if
      end if
    end do
    candidate%neutral_return_weight_scale = 1.0_dp
    candidate%neutral_return_unresolved_fraction = 0.0_dp
    do species_idx = 1_i32, candidate%nspecies
      if (candidate%neutral_return_correction(species_idx) == 0.0_dp) cycle
      if (candidate%emitted_from_surface(species_idx) >= 0.0_dp .or. &
          candidate%absorbed_on_surface(species_idx) >= 0.0_dp) then
        error stop 'cumulative neutral-return diagnostics require negative emitted and absorbed charge.'
      end if
      candidate%neutral_return_weight_scale(species_idx) = checked_charge_ratio( &
                                                           candidate%emitted_from_surface(species_idx), &
                                                           candidate%absorbed_on_surface(species_idx), &
                                                           'cumulative neutral-return weight scale' &
                                                           )
      candidate%neutral_return_unresolved_fraction(species_idx) = checked_charge_ratio( &
                                                                  candidate%discarded_unresolved(species_idx), &
                                                                  candidate%emitted_from_surface(species_idx), &
                                                                  'cumulative neutral-return unresolved fraction' &
                                                                  )
    end do
    call checked_accumulate_count_array(candidate%injected_count, batch%injected_count, 'injected count')
    call checked_accumulate_count_array(candidate%emitted_count, batch%emitted_count, 'emitted count')
    call checked_accumulate_count_array(candidate%absorbed_count, batch%absorbed_count, 'absorbed count')
    call checked_accumulate_count_array(candidate%escaped_count, batch%escaped_count, 'escaped count')
    call checked_accumulate_count_array( &
      candidate%discarded_unresolved_count, batch%discarded_unresolved_count, 'discarded unresolved count' &
      )
    call validate_ledger_finite(candidate, 'accumulated charge ledger')
    cumulative = candidate
  end subroutine accumulate_charge_ledger

  !> 有限なcharge incrementをoverflowさせずscalar accumulatorへ加える。
  subroutine checked_accumulate_charge(accumulator, increment, context)
    real(dp), intent(inout) :: accumulator
    real(dp), intent(in) :: increment
    character(len=*), intent(in) :: context

    if (.not. ieee_is_finite(accumulator) .or. .not. ieee_is_finite(increment)) then
      error stop trim(context)//' contains a non-finite charge.'
    end if
    if (increment > 0.0_dp .and. accumulator > huge(accumulator) - increment) then
      error stop trim(context)//' charge accumulation overflowed.'
    end if
    if (increment < 0.0_dp .and. accumulator < -huge(accumulator) - increment) then
      error stop trim(context)//' charge accumulation overflowed.'
    end if
    accumulator = accumulator + increment
    if (.not. ieee_is_finite(accumulator)) error stop trim(context)//' charge accumulation is not finite.'
  end subroutine checked_accumulate_charge

  !> scale後のNeumaier和により、中間overflowを避けてsigned charge総和を返す。
  pure real(dp) function finite_charge_sum(values, context) result(total)
    real(dp), intent(in) :: values(:)
    character(len=*), intent(in) :: context
    real(dp) :: compensation, scaled_sum, scaled_value, scale_value, updated_sum
    integer :: i

    if (.not. all(ieee_is_finite(values))) error stop trim(context)//' contains a non-finite charge.'
    if (size(values) == 0) then
      total = 0.0_dp
      return
    end if
    scale_value = maxval(abs(values))
    if (scale_value == 0.0_dp) then
      total = 0.0_dp
      return
    end if
    scaled_sum = 0.0_dp
    compensation = 0.0_dp
    do i = 1, size(values)
      scaled_value = values(i)/scale_value
      updated_sum = scaled_sum + scaled_value
      if (abs(scaled_sum) >= abs(scaled_value)) then
        compensation = compensation + (scaled_sum - updated_sum) + scaled_value
      else
        compensation = compensation + (scaled_value - updated_sum) + scaled_sum
      end if
      scaled_sum = updated_sum
    end do
    scaled_sum = scaled_sum + compensation
    if (scale_value >= 1.0_dp .and. abs(scaled_sum) > huge(total)/scale_value) then
      error stop trim(context)//' charge sum overflowed.'
    end if
    total = scale_value*scaled_sum
    if (.not. ieee_is_finite(total)) error stop trim(context)//' charge sum is not finite.'
  end function finite_charge_sum

  subroutine checked_accumulate_charge_array(accumulator, increment, context)
    real(dp), intent(inout) :: accumulator(:)
    real(dp), intent(in) :: increment(:)
    character(len=*), intent(in) :: context
    integer :: i

    if (size(accumulator) /= size(increment)) error stop trim(context)//' array size mismatch.'
    do i = 1, size(accumulator)
      call checked_accumulate_charge(accumulator(i), increment(i), context)
    end do
  end subroutine checked_accumulate_charge_array

  subroutine checked_accumulate_count_array(accumulator, increment, context)
    integer(i64), intent(inout) :: accumulator(:)
    integer(i64), intent(in) :: increment(:)
    character(len=*), intent(in) :: context
    integer :: i

    if (size(accumulator) /= size(increment)) error stop trim(context)//' array size mismatch.'
    do i = 1, size(accumulator)
      if (accumulator(i) < 0_i64 .or. increment(i) < 0_i64) then
        error stop trim(context)//' must be nonnegative.'
      end if
      if (accumulator(i) > huge(accumulator(i)) - increment(i)) then
        error stop trim(context)//' accumulation overflowed.'
      end if
      accumulator(i) = accumulator(i) + increment(i)
    end do
  end subroutine checked_accumulate_count_array

  pure integer(i64) function checked_nonnegative_count_sum(values, context) result(total)
    integer(i64), intent(in) :: values(:)
    character(len=*), intent(in) :: context
    integer :: i

    total = 0_i64
    do i = 1, size(values)
      if (values(i) < 0_i64) error stop trim(context)//' must be nonnegative.'
      if (total > huge(total) - values(i)) error stop trim(context)//' sum overflowed.'
      total = total + values(i)
    end do
  end function checked_nonnegative_count_sum

  real(dp) function checked_charge_ratio(numerator, denominator, context) result(ratio)
    real(dp), intent(in) :: numerator, denominator
    character(len=*), intent(in) :: context

    if (.not. ieee_is_finite(numerator) .or. .not. ieee_is_finite(denominator) .or. denominator == 0.0_dp) then
      error stop trim(context)//' operands are invalid.'
    end if
    if (abs(denominator) < 1.0_dp .and. abs(numerator) > huge(ratio)*abs(denominator)) then
      error stop trim(context)//' overflowed.'
    end if
    ratio = numerator/denominator
    if (.not. ieee_is_finite(ratio)) error stop trim(context)//' is not finite.'
  end function checked_charge_ratio

  subroutine validate_ledger_finite(ledger, context)
    type(charge_ledger_type), intent(in) :: ledger
    character(len=*), intent(in) :: context

    if (ledger%batch_count < 0_i32) error stop trim(context)//' batch count must be nonnegative.'
    if (.not. all(ieee_is_finite([ &
                                 ledger%surface_charge_before, ledger%surface_charge_after, &
                                 ledger%local_flight_charge_before, ledger%local_flight_charge_after, &
                                 ledger%unresolved_stock_before, ledger%unresolved_stock_after &
                                 ]))) error stop trim(context)//' stocks must be finite.'
    if (.not. all(ieee_is_finite(ledger%injected_from_remote)) .or. &
        .not. all(ieee_is_finite(ledger%emitted_from_surface)) .or. &
        .not. all(ieee_is_finite(ledger%absorbed_on_surface)) .or. &
        .not. all(ieee_is_finite(ledger%escaped_to_infinity)) .or. &
        .not. all(ieee_is_finite(ledger%discarded_unresolved)) .or. &
        .not. all(ieee_is_finite(ledger%neutral_return_correction)) .or. &
        .not. all(ieee_is_finite(ledger%neutral_return_weight_scale)) .or. &
        .not. all(ieee_is_finite(ledger%neutral_return_unresolved_fraction)) .or. &
        .not. all(ieee_is_finite(ledger%fixed_absorbed_target_charge)) .or. &
        .not. all(ieee_is_finite(ledger%fixed_absorbed_weight_scale)) .or. &
        .not. all(ieee_is_finite(ledger%fixed_emission_target_charge)) .or. &
        .not. all(ieee_is_finite(ledger%fixed_emission_weight_scale)) .or. &
        .not. all(ieee_is_finite(ledger%fixed_escape_target_charge)) .or. &
        .not. all(ieee_is_finite(ledger%fixed_escape_correction)) .or. &
        .not. all(ieee_is_finite(ledger%fixed_current_correction))) then
      error stop trim(context)//' diagnostics must be finite.'
    end if
    if (any(ledger%injected_count < 0_i64) .or. any(ledger%emitted_count < 0_i64) .or. &
        any(ledger%absorbed_count < 0_i64) .or. any(ledger%escaped_count < 0_i64) .or. &
        any(ledger%discarded_unresolved_count < 0_i64)) then
      error stop trim(context)//' counts must be nonnegative.'
    end if
  end subroutine validate_ledger_finite

  subroutine release_charge_ledger_arrays(self)
    type(charge_ledger_type), intent(inout) :: self

    deallocate ( &
      self%injected_from_remote, self%emitted_from_surface, self%absorbed_on_surface, &
      self%escaped_to_infinity, self%discarded_unresolved, self%neutral_return_correction, &
      self%neutral_return_weight_scale, self%neutral_return_unresolved_fraction, &
      self%fixed_absorbed_target_charge, self%fixed_absorbed_weight_scale, &
      self%fixed_emission_target_charge, self%fixed_emission_weight_scale, &
      self%fixed_escape_target_charge, self%fixed_escape_correction, &
      self%fixed_current_correction, &
      self%injected_count, self%emitted_count, self%absorbed_count, &
      self%escaped_count, self%discarded_unresolved_count &
      )
  end subroutine release_charge_ledger_arrays

end module bem_charge_ledger

!> 帰還粒子の発生元と着地点を区別して、平均電荷更新用の候補を構成する。
module bem_mean_charge_transaction
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  implicit none
  private

  integer(i32), parameter, public :: mean_charge_transaction_ok = 0_i32
  integer(i32), parameter, public :: mean_charge_transaction_invalid = 1_i32
  integer(i32), parameter, public :: mean_charge_transaction_size_mismatch = 2_i32
  integer(i32), parameter, public :: mean_charge_transaction_nonfinite = 3_i32
  integer(i32), parameter, public :: mean_charge_transaction_charge_mismatch = 4_i32

  !> 1 batchの予測電荷と、帰還着地後の実電荷transaction。
  !!
  !! `predictor_charge_c` は未確定の帰還分を発生元で一時的に中和した値、
  !! `actual_charge_c` は同じ帰還分を実際の着地点へdepositした値である。
  type, public :: mean_charge_transaction_type
    logical :: ready = .false.
    integer(i32) :: element_count = 0_i32
    real(dp), allocatable :: predictor_charge_c(:)
    real(dp), allocatable :: actual_charge_c(:)
    real(dp) :: pending_total_charge_c = 0.0_dp
    real(dp) :: deferred_source_total_charge_c = 0.0_dp
    real(dp) :: returned_destination_total_charge_c = 0.0_dp
    real(dp) :: predictor_total_charge_c = 0.0_dp
    real(dp) :: actual_total_charge_c = 0.0_dp
    real(dp) :: source_return_residual_charge_c = 0.0_dp
    real(dp) :: predictor_actual_residual_charge_c = 0.0_dp
    real(dp) :: charge_tolerance_c = 0.0_dp
  end type mean_charge_transaction_type

  public :: build_mean_charge_transaction

contains

  !> pending、帰還発生元、帰還着地点からpredictorとactualを構成する。
  !!
  !! `deferred_source_charge_c` は帰還割合を適用済みの正電荷であり、
  !! `returned_destination_charge_c` は同じmacro chargeの負depositである。
  !! escape分は `pending_charge_c` に残す。両者の総量が打ち消し合わない場合は、
  !! 任意のsupport面へ残差を再分配せず、明示的な失敗として返す。
  subroutine build_mean_charge_transaction( &
    pending_charge_c, deferred_source_charge_c, returned_destination_charge_c, &
    transaction, status, message &
    )
    real(dp), intent(in) :: pending_charge_c(:)
    real(dp), intent(in) :: deferred_source_charge_c(:)
    real(dp), intent(in) :: returned_destination_charge_c(:)
    type(mean_charge_transaction_type), intent(out) :: transaction
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    integer(i32) :: element_count
    real(dp), allocatable :: predictor_charge_c(:), actual_charge_c(:)
    real(dp) :: pending_total, deferred_source_total, returned_destination_total
    real(dp) :: predictor_total, actual_total
    real(dp) :: source_return_residual, predictor_actual_residual
    real(dp) :: charge_scale, charge_tolerance

    transaction = mean_charge_transaction_type()
    status = mean_charge_transaction_invalid
    message = ''

    element_count = int(size(pending_charge_c), i32)
    if (element_count < 1_i32) then
      message = 'mean-charge transaction requires at least one element'
      return
    end if
    if (size(deferred_source_charge_c) /= size(pending_charge_c) .or. &
        size(returned_destination_charge_c) /= size(pending_charge_c)) then
      status = mean_charge_transaction_size_mismatch
      message = 'mean-charge transaction vectors must have identical sizes'
      return
    end if
    if (any(.not. ieee_is_finite(pending_charge_c)) .or. &
        any(.not. ieee_is_finite(deferred_source_charge_c)) .or. &
        any(.not. ieee_is_finite(returned_destination_charge_c))) then
      status = mean_charge_transaction_nonfinite
      message = 'mean-charge transaction received a non-finite charge'
      return
    end if
    if (any(deferred_source_charge_c < 0.0_dp)) then
      message = 'mean-charge deferred source charge must be nonnegative'
      return
    end if
    if (any(returned_destination_charge_c > 0.0_dp)) then
      message = 'mean-charge returned destination charge must be nonpositive'
      return
    end if

    allocate (predictor_charge_c(element_count), actual_charge_c(element_count))
    predictor_charge_c = pending_charge_c - deferred_source_charge_c
    actual_charge_c = pending_charge_c + returned_destination_charge_c
    if (any(.not. ieee_is_finite(predictor_charge_c)) .or. &
        any(.not. ieee_is_finite(actual_charge_c))) then
      status = mean_charge_transaction_nonfinite
      message = 'mean-charge transaction candidate overflowed'
      return
    end if

    pending_total = sum(pending_charge_c)
    deferred_source_total = sum(deferred_source_charge_c)
    returned_destination_total = sum(returned_destination_charge_c)
    predictor_total = sum(predictor_charge_c)
    actual_total = sum(actual_charge_c)
    source_return_residual = deferred_source_total + returned_destination_total
    predictor_actual_residual = actual_total - predictor_total
    charge_scale = max( &
                   sum(abs(pending_charge_c)), sum(abs(deferred_source_charge_c)), &
                   sum(abs(returned_destination_charge_c)), sum(abs(predictor_charge_c)), &
                   sum(abs(actual_charge_c)), tiny(1.0_dp) &
                   )
    if (.not. all(ieee_is_finite([ &
                                 pending_total, deferred_source_total, returned_destination_total, &
                                 predictor_total, actual_total, source_return_residual, &
                                 predictor_actual_residual, charge_scale &
                                 ]))) then
      status = mean_charge_transaction_nonfinite
      message = 'mean-charge transaction totals are not finite'
      return
    end if
    charge_tolerance = 256.0_dp*epsilon(1.0_dp)*charge_scale
    if (.not. ieee_is_finite(charge_tolerance)) then
      status = mean_charge_transaction_nonfinite
      message = 'mean-charge transaction tolerance is not finite'
      return
    end if
    if (abs(source_return_residual) > charge_tolerance .or. &
        abs(predictor_actual_residual) > charge_tolerance) then
      status = mean_charge_transaction_charge_mismatch
      message = 'mean-charge returned charge does not cancel its deferred source charge'
      return
    end if

    transaction%element_count = element_count
    call move_alloc(predictor_charge_c, transaction%predictor_charge_c)
    call move_alloc(actual_charge_c, transaction%actual_charge_c)
    transaction%pending_total_charge_c = pending_total
    transaction%deferred_source_total_charge_c = deferred_source_total
    transaction%returned_destination_total_charge_c = returned_destination_total
    transaction%predictor_total_charge_c = predictor_total
    transaction%actual_total_charge_c = actual_total
    transaction%source_return_residual_charge_c = source_return_residual
    transaction%predictor_actual_residual_charge_c = predictor_actual_residual
    transaction%charge_tolerance_c = charge_tolerance
    transaction%ready = .true.
    status = mean_charge_transaction_ok
    message = 'mean-charge transaction constructed'
  end subroutine build_mean_charge_transaction

end module bem_mean_charge_transaction

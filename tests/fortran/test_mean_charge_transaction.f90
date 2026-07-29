program test_mean_charge_transaction
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use bem_kinds, only: dp, i32
  use bem_mean_charge_transaction, only: &
    mean_charge_transaction_ok, mean_charge_transaction_invalid, &
    mean_charge_transaction_size_mismatch, mean_charge_transaction_nonfinite, &
    mean_charge_transaction_charge_mismatch, mean_charge_transaction_type, &
    build_mean_charge_transaction
  use test_support, only: &
    test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, &
    assert_equal_i32, assert_allclose_1d
  implicit none

  call test_init(5)

  call test_begin('source-aware predictor preserves escape charge without a support plane')
  call test_source_aware_candidates()
  call test_end()

  call test_begin('distributed source and destination charges preserve the total')
  call test_distributed_return()
  call test_end()

  call test_begin('vector size mismatch fails closed')
  call test_size_mismatch()
  call test_end()

  call test_begin('non-finite input fails closed')
  call test_nonfinite_input()
  call test_end()

  call test_begin('unbalanced return charge fails closed')
  call test_unbalanced_return()
  call test_end()

  call test_summary()

contains

  subroutine test_source_aware_candidates()
    type(mean_charge_transaction_type) :: transaction
    real(dp) :: pending(4), deferred_source(4), returned_destination(4)
    real(dp) :: expected_predictor(4), expected_actual(4)
    integer(i32) :: status
    character(len=256) :: message

    ! element 2 の放出反作用10 Cのうち6 Cが帰還し、4 Cがescapeする。
    pending = [2.0_dp, 10.0_dp, -1.0_dp, 0.5_dp]
    deferred_source = [0.0_dp, 6.0_dp, 0.0_dp, 0.0_dp]
    returned_destination = [0.0_dp, 0.0_dp, -2.0_dp, -4.0_dp]
    expected_predictor = [2.0_dp, 4.0_dp, -1.0_dp, 0.5_dp]
    expected_actual = [2.0_dp, 10.0_dp, -3.0_dp, -3.5_dp]

    call build_mean_charge_transaction( &
      pending, deferred_source, returned_destination, transaction, status, message &
      )

    call assert_equal_i32(status, mean_charge_transaction_ok, &
                          'source-aware transaction failed: '//trim(message))
    call assert_true(transaction%ready, 'successful transaction is not ready')
    call assert_equal_i32(transaction%element_count, 4_i32, 'transaction element count mismatch')
    call assert_allclose_1d(transaction%predictor_charge_c, expected_predictor, 0.0_dp, &
                            'predictor changed unrelated or escaping charge')
    call assert_allclose_1d(transaction%actual_charge_c, expected_actual, 0.0_dp, &
                            'actual transaction did not use real return destinations')
    call assert_close_dp(transaction%predictor_total_charge_c, 5.5_dp, 0.0_dp, &
                         'predictor total charge mismatch')
    call assert_close_dp(transaction%actual_total_charge_c, 5.5_dp, 0.0_dp, &
                         'actual total charge mismatch')
    call assert_close_dp(transaction%source_return_residual_charge_c, 0.0_dp, &
                         transaction%charge_tolerance_c, 'source/return invariant did not close')
    call assert_close_dp(transaction%predictor_actual_residual_charge_c, 0.0_dp, &
                         transaction%charge_tolerance_c, 'predictor/actual invariant did not close')
  end subroutine test_source_aware_candidates

  subroutine test_distributed_return()
    type(mean_charge_transaction_type) :: transaction
    real(dp) :: pending(6), deferred_source(6), returned_destination(6)
    integer(i32) :: status
    character(len=256) :: message

    pending = [0.25_dp, 3.0_dp, -0.5_dp, 5.0_dp, 0.0_dp, 0.75_dp]
    deferred_source = [0.0_dp, 1.25_dp, 0.0_dp, 2.75_dp, 0.0_dp, 0.0_dp]
    returned_destination = [-0.5_dp, 0.0_dp, -1.0_dp, 0.0_dp, -2.5_dp, 0.0_dp]

    call build_mean_charge_transaction( &
      pending, deferred_source, returned_destination, transaction, status, message &
      )

    call assert_equal_i32(status, mean_charge_transaction_ok, &
                          'distributed return transaction failed: '//trim(message))
    call assert_close_dp( &
      transaction%predictor_total_charge_c, transaction%actual_total_charge_c, &
      transaction%charge_tolerance_c, 'distributed return changed total charge' &
      )
    call assert_close_dp(transaction%predictor_charge_c(6), pending(6), 0.0_dp, &
                         'transaction changed an unrelated element')
  end subroutine test_distributed_return

  subroutine test_size_mismatch()
    type(mean_charge_transaction_type) :: transaction
    integer(i32) :: status
    character(len=256) :: message

    call build_mean_charge_transaction( &
      [1.0_dp, 2.0_dp, 3.0_dp], [0.0_dp, 1.0_dp], [0.0_dp, -1.0_dp, 0.0_dp], &
      transaction, status, message &
      )

    call assert_equal_i32(status, mean_charge_transaction_size_mismatch, &
                          'size mismatch returned the wrong status')
    call assert_true(.not. transaction%ready, 'size mismatch exposed a ready transaction')
    call assert_true(.not. allocated(transaction%predictor_charge_c), &
                     'size mismatch exposed a predictor candidate')
    call assert_true(.not. allocated(transaction%actual_charge_c), &
                     'size mismatch exposed an actual candidate')
  end subroutine test_size_mismatch

  subroutine test_nonfinite_input()
    type(mean_charge_transaction_type) :: transaction
    real(dp) :: pending(2)
    integer(i32) :: status
    character(len=256) :: message

    pending = [1.0_dp, ieee_value(0.0_dp, ieee_quiet_nan)]
    call build_mean_charge_transaction( &
      pending, [1.0_dp, 0.0_dp], [-1.0_dp, 0.0_dp], transaction, status, message &
      )

    call assert_equal_i32(status, mean_charge_transaction_nonfinite, &
                          'non-finite input returned the wrong status')
    call assert_true(.not. transaction%ready, 'non-finite input exposed a ready transaction')
    call assert_true(.not. allocated(transaction%predictor_charge_c), &
                     'non-finite input exposed a predictor candidate')
    call assert_true(.not. allocated(transaction%actual_charge_c), &
                     'non-finite input exposed an actual candidate')
  end subroutine test_nonfinite_input

  subroutine test_unbalanced_return()
    type(mean_charge_transaction_type) :: transaction
    integer(i32) :: status
    character(len=256) :: message

    call build_mean_charge_transaction( &
      [4.0_dp, 0.0_dp], [2.0_dp, 0.0_dp], [0.0_dp, -1.5_dp], &
      transaction, status, message &
      )

    call assert_equal_i32(status, mean_charge_transaction_charge_mismatch, &
                          'unbalanced return returned the wrong status')
    call assert_true(.not. transaction%ready, 'unbalanced return exposed a ready transaction')
    call assert_true(.not. allocated(transaction%predictor_charge_c), &
                     'unbalanced return exposed a predictor candidate')
    call assert_true(.not. allocated(transaction%actual_charge_c), &
                     'unbalanced return exposed an actual candidate')

    call build_mean_charge_transaction( &
      [1.0_dp], [-1.0_dp], [1.0_dp], transaction, status, message &
      )
    call assert_equal_i32(status, mean_charge_transaction_invalid, &
                          'invalid source sign returned the wrong status')
    call assert_true(.not. transaction%ready, 'invalid charge signs exposed a ready transaction')
  end subroutine test_unbalanced_return

end program test_mean_charge_transaction

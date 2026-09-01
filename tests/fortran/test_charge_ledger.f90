!> 電荷収支の符号規約、stock 差分、累積を検証する。
program test_charge_ledger
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use bem_kinds, only: dp, i32, i64
  use bem_charge_ledger, only: charge_ledger_type, accumulate_charge_ledger, finite_charge_sum
  use test_support, only: &
    test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32, assert_close_dp, assert_allclose_1d
  implicit none

  type(charge_ledger_type) :: ledger, batch2, cumulative
  character(len=64) :: run_mode

  call get_command_argument(1, run_mode)
  select case (trim(run_mode))
  case ('--real-overflow-probe')
    call run_real_overflow_probe()
    stop 0
  case ('--count-overflow-probe')
    call run_count_overflow_probe()
    stop 0
  case ('--ratio-overflow-probe')
    call run_ratio_overflow_probe()
    stop 0
  case ('--negative-count-probe')
    call run_negative_count_probe()
    stop 0
  case ('--nonfinite-sum-probe')
    call run_nonfinite_sum_probe()
    stop 0
  case ('--sum-overflow-probe')
    call run_sum_overflow_probe()
    stop 0
  end select

  call test_init(10)

  call test_begin('initialization_and_reset_defaults')
  call ledger%init(2_i32)
  call assert_equal_i32(ledger%nspecies, 2_i32, 'ledger species count mismatch')
  call assert_reset_defaults(ledger, 0_i32, 'initialized ledger')
  ledger%surface_charge_before = 2.0_dp
  ledger%surface_charge_after = 2.0_dp
  ledger%local_flight_charge_before = 2.0_dp
  ledger%local_flight_charge_after = 2.0_dp
  ledger%unresolved_stock_before = 2.0_dp
  ledger%unresolved_stock_after = 2.0_dp
  ledger%injected_from_remote = 2.0_dp
  ledger%emitted_from_surface = 2.0_dp
  ledger%absorbed_on_surface = 2.0_dp
  ledger%escaped_to_infinity = 2.0_dp
  ledger%discarded_unresolved = 2.0_dp
  ledger%neutral_return_correction = 2.0_dp
  ledger%neutral_return_weight_scale = 2.0_dp
  ledger%neutral_return_unresolved_fraction = 2.0_dp
  ledger%fixed_absorbed_target_charge = 2.0_dp
  ledger%fixed_absorbed_weight_scale = 2.0_dp
  ledger%fixed_emission_target_charge = 2.0_dp
  ledger%fixed_emission_weight_scale = 2.0_dp
  ledger%fixed_escape_target_charge = 2.0_dp
  ledger%fixed_escape_correction = 2.0_dp
  ledger%fixed_current_correction = 2.0_dp
  ledger%injected_count = 2_i64
  ledger%emitted_count = 2_i64
  ledger%absorbed_count = 2_i64
  ledger%escaped_count = 2_i64
  ledger%discarded_unresolved_count = 2_i64
  call ledger%reset(7_i32)
  call assert_reset_defaults(ledger, 7_i32, 'reset ledger')
  call test_end()

  call test_begin('cumulative_fixed_current_diagnostics')
  cumulative = charge_ledger_type()
  call ledger%init(2_i32)
  call ledger%reset(1_i32)
  ledger%absorbed_on_surface(1) = -2.0_dp
  ledger%fixed_absorbed_target_charge(1) = -1.0_dp
  ledger%fixed_absorbed_weight_scale(1) = 7.0_dp
  ledger%emitted_from_surface(2) = -4.0_dp
  ledger%fixed_emission_target_charge(2) = 8.0_dp
  ledger%fixed_emission_weight_scale(2) = 8.0_dp
  ledger%fixed_escape_target_charge(2) = -1.0_dp
  ledger%fixed_escape_correction(2) = -0.25_dp
  ledger%fixed_current_correction(1) = -0.5_dp
  call batch2%init(2_i32)
  call batch2%reset(2_i32)
  batch2%absorbed_on_surface(1) = -6.0_dp
  batch2%fixed_absorbed_target_charge(1) = -3.0_dp
  batch2%fixed_absorbed_weight_scale(1) = 9.0_dp
  batch2%emitted_from_surface(2) = -1.0_dp
  batch2%fixed_emission_target_charge(2) = 2.0_dp
  batch2%fixed_emission_weight_scale(2) = 10.0_dp
  batch2%fixed_escape_target_charge(2) = -2.0_dp
  batch2%fixed_escape_correction(2) = -0.5_dp
  batch2%fixed_current_correction(1) = -1.0_dp
  call accumulate_charge_ledger(cumulative, ledger)
  call accumulate_charge_ledger(cumulative, batch2)
  call assert_allclose_1d( &
    [ &
    cumulative%fixed_absorbed_weight_scale(1), cumulative%fixed_emission_weight_scale(2), &
    cumulative%fixed_escape_target_charge(2), cumulative%fixed_escape_correction(2), &
    cumulative%fixed_current_correction(1) &
    ], &
    [0.5_dp, 2.0_dp, -3.0_dp, -0.75_dp, -1.5_dp], 1.0d-14, &
    'cumulative fixed-current diagnostics mismatch' &
    )
  call test_end()

  call test_begin('residual_channel_signs_and_ownership')
  call ledger%init(2_i32)
  call ledger%reset(1_i32)
  ledger%surface_charge_before = 1.0_dp
  ledger%surface_charge_after = 4.0_dp
  ledger%local_flight_charge_before = 2.0_dp
  ledger%local_flight_charge_after = 7.0_dp
  ledger%unresolved_stock_before = -3.0_dp
  ledger%unresolved_stock_after = -2.0_dp
  ledger%injected_from_remote = [-4.0_dp, 6.0_dp]
  ledger%escaped_to_infinity = [-7.0_dp, 8.0_dp]
  ledger%discarded_unresolved = [-9.0_dp, 10.0_dp]
  ledger%neutral_return_correction = [-11.0_dp, 12.0_dp]
  ledger%fixed_current_correction = [-13.0_dp, 14.0_dp]
  call assert_close_dp(ledger%residual(), 7.0_dp, 1.0d-14, 'owned residual channel signs mismatch')
  ! Raw surface channels are diagnostics; their applied effects already appear in stocks/corrections.
  ledger%emitted_from_surface = [101.0_dp, 103.0_dp]
  ledger%absorbed_on_surface = [107.0_dp, 110.0_dp]
  call assert_close_dp(ledger%residual(), 7.0_dp, 1.0d-14, 'raw surface diagnostics entered residual')
  ledger%fixed_absorbed_target_charge = [17.0_dp, 19.0_dp]
  ledger%fixed_absorbed_weight_scale = [2.0_dp, 3.0_dp]
  ledger%fixed_emission_target_charge = [23.0_dp, 29.0_dp]
  ledger%fixed_emission_weight_scale = [4.0_dp, 5.0_dp]
  ledger%fixed_escape_target_charge = [31.0_dp, 37.0_dp]
  ledger%fixed_escape_correction = [41.0_dp, 43.0_dp]
  call assert_close_dp(ledger%residual(), 7.0_dp, 1.0d-14, 'fixed-current diagnostics entered residual')
  call test_end()

  call test_begin('discard_diagnostics_do_not_cancel_between_species')
  call ledger%init(2_i32)
  call ledger%reset(1_i32)
  ledger%injected_from_remote = [-1.0_dp, 2.0_dp]
  ledger%discarded_unresolved = [-1.0_dp, 2.0_dp]
  ledger%discarded_unresolved_count = [1_i64, 2_i64]
  call assert_close_dp(ledger%residual(), 0.0_dp, 0.0_dp, 'discarded charge residual mismatch')
  call assert_close_dp(ledger%discarded_unresolved_abs(), 3.0_dp, 0.0_dp, 'discarded absolute charge mismatch')
  call assert_true( &
    .not. ledger%has_unresolved_discard(3.0_dp, 3_i64), &
    'discard threshold equality must be accepted' &
    )
  call assert_true( &
    ledger%has_unresolved_discard(2.5_dp, 3_i64), &
    'discard absolute-charge threshold must not use signed cancellation' &
    )
  call assert_true( &
    ledger%has_unresolved_discard(huge(1.0_dp), 2_i64), &
    'discard count threshold must remain independent of charge residual' &
    )
  call test_end()

  call test_begin('cumulative_batches')
  cumulative = charge_ledger_type()
  call ledger%init(2_i32)
  call ledger%reset(1_i32)
  ledger%surface_charge_before = 0.0_dp
  ledger%surface_charge_after = -2.0_dp
  ledger%local_flight_charge_before = 1.0_dp
  ledger%local_flight_charge_after = 3.0_dp
  ledger%unresolved_stock_before = -4.0_dp
  ledger%unresolved_stock_after = -6.0_dp
  ledger%injected_from_remote(1) = -3.0_dp
  ledger%emitted_from_surface = [-3.0_dp, -5.0_dp]
  ledger%absorbed_on_surface = [-7.0_dp, -11.0_dp]
  ledger%discarded_unresolved(1) = 0.25_dp
  ledger%fixed_current_correction(1) = 0.5_dp
  ledger%injected_count = [2_i64, 3_i64]
  ledger%emitted_count = [4_i64, 5_i64]
  ledger%absorbed_count = [6_i64, 7_i64]
  ledger%escaped_count = [8_i64, 9_i64]
  ledger%discarded_unresolved_count = [10_i64, 11_i64]
  call batch2%init(2_i32)
  call batch2%reset(2_i32)
  batch2%surface_charge_before = -2.0_dp
  batch2%surface_charge_after = -1.0_dp
  batch2%local_flight_charge_before = 3.0_dp
  batch2%local_flight_charge_after = 2.0_dp
  batch2%unresolved_stock_before = -6.0_dp
  batch2%unresolved_stock_after = -5.0_dp
  batch2%escaped_to_infinity(2) = -1.0_dp
  batch2%emitted_from_surface = [-13.0_dp, -17.0_dp]
  batch2%absorbed_on_surface = [-19.0_dp, -23.0_dp]
  batch2%discarded_unresolved(2) = 0.75_dp
  batch2%fixed_current_correction(2) = 1.5_dp
  batch2%injected_count = [1_i64, 2_i64]
  batch2%emitted_count = [3_i64, 4_i64]
  batch2%absorbed_count = [5_i64, 6_i64]
  batch2%escaped_count = [7_i64, 8_i64]
  batch2%discarded_unresolved_count = [9_i64, 10_i64]
  call accumulate_charge_ledger(cumulative, ledger)
  call accumulate_charge_ledger(cumulative, batch2)
  call assert_equal_i32(cumulative%batch_count, 2_i32, 'cumulative batch count mismatch')
  call assert_allclose_1d( &
    [ &
    cumulative%surface_charge_before, cumulative%surface_charge_after, &
    cumulative%local_flight_charge_before, cumulative%local_flight_charge_after, &
    cumulative%unresolved_stock_before, cumulative%unresolved_stock_after &
    ], &
    [0.0_dp, -1.0_dp, 1.0_dp, 2.0_dp, -4.0_dp, -5.0_dp], 0.0_dp, &
    'cumulative stock endpoints mismatch' &
    )
  call assert_allclose_1d( &
    [ &
    cumulative%injected_from_remote, cumulative%emitted_from_surface, cumulative%absorbed_on_surface, &
    cumulative%escaped_to_infinity, cumulative%discarded_unresolved, &
    cumulative%fixed_current_correction &
    ], &
    [-3.0_dp, 0.0_dp, -16.0_dp, -22.0_dp, -26.0_dp, -34.0_dp, &
     0.0_dp, -1.0_dp, 0.25_dp, 0.75_dp, 0.5_dp, 1.5_dp], &
    0.0_dp, 'cumulative charge channels mismatch' &
    )
  call assert_true( &
    all(cumulative%injected_count == [3_i64, 5_i64]) .and. &
    all(cumulative%emitted_count == [7_i64, 9_i64]) .and. &
    all(cumulative%absorbed_count == [11_i64, 13_i64]) .and. &
    all(cumulative%escaped_count == [15_i64, 17_i64]) .and. &
    all(cumulative%discarded_unresolved_count == [19_i64, 21_i64]), &
    'cumulative count channels mismatch' &
    )
  call assert_close_dp(cumulative%residual(), 0.0_dp, 1.0d-14, 'cumulative residual mismatch')
  call test_end()

  call test_begin('cumulative_neutral_return_diagnostics')
  cumulative = charge_ledger_type()
  call ledger%init(2_i32)
  call ledger%reset(1_i32)
  ledger%emitted_from_surface(2) = -8.0_dp
  ledger%absorbed_on_surface(2) = -8.0_dp
  call batch2%init(2_i32)
  call batch2%reset(2_i32)
  batch2%emitted_from_surface(2) = -2.0_dp
  batch2%absorbed_on_surface(2) = -1.0_dp
  batch2%discarded_unresolved(2) = -1.0_dp
  batch2%neutral_return_correction(2) = -1.0_dp
  call accumulate_charge_ledger(cumulative, ledger)
  call accumulate_charge_ledger(cumulative, batch2)
  call assert_allclose_1d( &
    [ &
    cumulative%emitted_from_surface(2), cumulative%absorbed_on_surface(2), &
    cumulative%discarded_unresolved(2), cumulative%neutral_return_correction(2), &
    cumulative%neutral_return_weight_scale(2), cumulative%neutral_return_unresolved_fraction(2) &
    ], &
    [-10.0_dp, -9.0_dp, -1.0_dp, -1.0_dp, 10.0_dp/9.0_dp, 0.1_dp], 1.0d-14, &
    'cumulative neutral-return diagnostics mismatch' &
    )
  call assert_close_dp(cumulative%residual(), 0.0_dp, 1.0d-14, 'cumulative neutral-return residual mismatch')
  call test_end()

  call test_begin('cumulative_neutral_return_diagnostics_degrade_without_stopping')
  cumulative = charge_ledger_type()
  call ledger%init(1_i32)
  call ledger%reset(1_i32)
  ledger%emitted_from_surface(1) = 1.0_dp
  ledger%absorbed_on_surface(1) = -1.0_dp
  ledger%neutral_return_correction(1) = 0.5_dp
  call accumulate_charge_ledger(cumulative, ledger)
  call assert_close_dp( &
    cumulative%neutral_return_correction(1), 0.5_dp, 0.0_dp, &
    'neutral-return correction was lost when derived diagnostics were unavailable' &
    )
  call assert_close_dp( &
    cumulative%neutral_return_weight_scale(1), 1.0_dp, 0.0_dp, &
    'unavailable neutral-return weight scale did not retain its neutral default' &
    )
  call assert_close_dp( &
    cumulative%neutral_return_unresolved_fraction(1), 0.0_dp, 0.0_dp, &
    'unavailable neutral-return unresolved fraction did not retain its neutral default' &
    )
  call test_end()

  call test_begin('finite_charge_sum_handles_overflow_scale_and_cancellation')
  call assert_close_dp( &
    finite_charge_sum( &
    [0.75_dp*huge(1.0_dp), 0.75_dp*huge(1.0_dp), -0.75_dp*huge(1.0_dp), -0.75_dp*huge(1.0_dp)], &
    'scaled signed sum test' &
    ), &
    0.0_dp, 0.0_dp, 'scaled signed charge sum mismatch' &
    )
  call assert_close_dp( &
    finite_charge_sum([1.0e16_dp, 1.0_dp, -1.0e16_dp], 'compensated signed sum test'), &
    1.0_dp, 1.0e-12_dp, 'compensated signed charge sum mismatch' &
    )
  call assert_close_dp( &
    finite_charge_sum([1.0_dp, 1.0e16_dp, -1.0e16_dp], 'reverse compensated signed sum test'), &
    1.0_dp, 1.0e-12_dp, 'reverse compensated signed charge sum mismatch' &
    )
  call test_end()

  call test_begin('invalid_numeric_accumulation_fails_closed')
  call assert_probe_fails('--real-overflow-probe', 'real ledger overflow probe must fail')
  call assert_probe_fails('--count-overflow-probe', 'count ledger overflow probe must fail')
  call assert_probe_fails('--ratio-overflow-probe', 'derived ledger ratio overflow probe must fail')
  call assert_probe_fails('--nonfinite-sum-probe', 'nonfinite charge sum probe must fail')
  call assert_probe_fails('--sum-overflow-probe', 'charge sum overflow probe must fail')
  call test_end()

  call test_begin('negative_count_fails_at_accumulation_preflight')
  call assert_probe_fails('--negative-count-probe', 'negative ledger count probe must fail')
  call test_end()

  call test_summary()

contains

  subroutine assert_reset_defaults(value, expected_batch_count, label)
    type(charge_ledger_type), intent(in) :: value
    integer(i32), intent(in) :: expected_batch_count
    character(len=*), intent(in) :: label

    call assert_equal_i32(value%batch_count, expected_batch_count, trim(label)//' batch count mismatch')
    call assert_close_dp(value%residual(), 0.0_dp, 0.0_dp, trim(label)//' residual mismatch')
    call assert_true( &
      all([ &
          value%surface_charge_before, value%surface_charge_after, &
          value%local_flight_charge_before, value%local_flight_charge_after, &
          value%unresolved_stock_before, value%unresolved_stock_after, &
          value%injected_from_remote, value%emitted_from_surface, value%absorbed_on_surface, &
          value%escaped_to_infinity, value%discarded_unresolved, value%neutral_return_correction, &
          value%neutral_return_unresolved_fraction, value%fixed_absorbed_target_charge, &
          value%fixed_emission_target_charge, value%fixed_escape_target_charge, &
          value%fixed_escape_correction, value%fixed_current_correction &
          ] == 0.0_dp), &
      trim(label)//' zero-valued defaults mismatch' &
      )
    call assert_true( &
      all(value%neutral_return_weight_scale == 1.0_dp) .and. &
      all(value%fixed_absorbed_weight_scale == 1.0_dp) .and. &
      all(value%fixed_emission_weight_scale == 1.0_dp), &
      trim(label)//' unit scale defaults mismatch' &
      )
    call assert_true( &
      all(value%injected_count == 0_i64) .and. all(value%emitted_count == 0_i64) .and. &
      all(value%absorbed_count == 0_i64) .and. all(value%escaped_count == 0_i64) .and. &
      all(value%discarded_unresolved_count == 0_i64), &
      trim(label)//' count defaults mismatch' &
      )
  end subroutine assert_reset_defaults

  subroutine assert_probe_fails(mode, message)
    character(len=*), intent(in) :: mode, message
    character(len=1024) :: executable_path, command
    integer :: exit_status, command_status

    call get_command_argument(0, executable_path)
    command = '"'//trim(executable_path)//'" '//trim(mode)//' > /dev/null 2>&1'
    call execute_command_line(trim(command), wait=.true., exitstat=exit_status, cmdstat=command_status)
    call assert_equal_i32(int(command_status, i32), 0_i32, trim(message)//' command status')
    call assert_true(exit_status /= 0, message)
  end subroutine assert_probe_fails

  subroutine run_real_overflow_probe()
    type(charge_ledger_type) :: accumulated, increment

    call accumulated%init(1_i32)
    call accumulated%reset(1_i32)
    call increment%init(1_i32)
    call increment%reset(2_i32)
    increment%injected_from_remote(1) = 0.0096_dp*huge(1.0_dp)
    ! The subtraction rounds to the precheck boundary, while adding the same increment still overflows.
    accumulated%injected_from_remote(1) = huge(1.0_dp) - increment%injected_from_remote(1)
    call accumulate_charge_ledger(accumulated, increment)
  end subroutine run_real_overflow_probe

  subroutine run_count_overflow_probe()
    type(charge_ledger_type) :: accumulated, increment

    call accumulated%init(1_i32)
    call accumulated%reset(1_i32)
    accumulated%injected_count(1) = huge(0_i64)
    call increment%init(1_i32)
    call increment%reset(2_i32)
    increment%injected_count(1) = 1_i64
    call accumulate_charge_ledger(accumulated, increment)
  end subroutine run_count_overflow_probe

  subroutine run_ratio_overflow_probe()
    type(charge_ledger_type) :: accumulated, increment

    accumulated = charge_ledger_type()
    call increment%init(1_i32)
    call increment%reset(1_i32)
    increment%absorbed_on_surface(1) = -tiny(1.0_dp)
    increment%fixed_absorbed_target_charge(1) = -huge(1.0_dp)
    increment%fixed_absorbed_weight_scale(1) = 2.0_dp
    call accumulate_charge_ledger(accumulated, increment)
  end subroutine run_ratio_overflow_probe

  subroutine run_negative_count_probe()
    type(charge_ledger_type) :: invalid_batch, destination

    call invalid_batch%init(1_i32)
    call invalid_batch%reset(1_i32)
    invalid_batch%injected_count(1) = -1_i64
    call accumulate_charge_ledger(destination, invalid_batch)
  end subroutine run_negative_count_probe

  subroutine run_nonfinite_sum_probe()
    real(dp) :: ignored

    ignored = finite_charge_sum([ieee_value(0.0_dp, ieee_quiet_nan)], 'nonfinite sum probe')
    write (*, *) ignored
  end subroutine run_nonfinite_sum_probe

  subroutine run_sum_overflow_probe()
    real(dp) :: ignored

    ignored = finite_charge_sum([huge(1.0_dp), huge(1.0_dp)], 'sum overflow probe')
    write (*, *) ignored
  end subroutine run_sum_overflow_probe
end program test_charge_ledger

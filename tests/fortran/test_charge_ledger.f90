!> 電荷収支の符号規約、stock 差分、累積を検証する。
program test_charge_ledger
  use bem_kinds, only: dp, i32, i64
  use bem_charge_ledger, only: charge_ledger_type, accumulate_charge_ledger
  use test_support, only: &
    test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32, assert_equal_i64, assert_close_dp
  implicit none

  type(charge_ledger_type) :: ledger, batch2, cumulative

  call test_init(10)

  call test_begin('initialization')
  call ledger%init(2_i32)
  call ledger%reset(1_i32)
  call assert_equal_i32(ledger%nspecies, 2_i32, 'ledger species count mismatch')
  call assert_equal_i32(ledger%batch_count, 1_i32, 'ledger batch count mismatch')
  call assert_close_dp(ledger%residual(), 0.0_dp, 0.0_dp, 'empty ledger residual mismatch')
  call assert_equal_i64(sum(ledger%discarded_unresolved_count), 0_i64, 'empty discard count mismatch')
  call test_end()

  call test_begin('cumulative_fixed_current_diagnostics')
  cumulative = charge_ledger_type()
  call ledger%reset(1_i32)
  ledger%absorbed_on_surface(1) = -2.0_dp
  ledger%fixed_absorbed_target_charge(1) = -5.0_dp
  ledger%fixed_absorbed_applied_charge(1) = -5.0_dp
  ledger%fixed_absorbed_weight_scale(1) = 2.5_dp
  ledger%emitted_from_surface(2) = -4.0_dp
  ledger%fixed_emission_target_charge(2) = 2.0_dp
  ledger%fixed_emission_applied_charge(2) = 2.0_dp
  ledger%fixed_emission_weight_scale(2) = 0.5_dp
  ledger%fixed_escape_target_charge(2) = -1.0_dp
  ledger%fixed_escape_applied_charge(2) = -1.0_dp
  ledger%fixed_escape_correction(2) = -0.25_dp
  call batch2%init(2_i32)
  call batch2%reset(2_i32)
  batch2%absorbed_on_surface(1) = -6.0_dp
  batch2%fixed_absorbed_target_charge(1) = -3.0_dp
  batch2%fixed_absorbed_applied_charge(1) = -3.0_dp
  batch2%fixed_absorbed_weight_scale(1) = 0.5_dp
  batch2%emitted_from_surface(2) = -1.0_dp
  batch2%fixed_emission_target_charge(2) = 3.0_dp
  batch2%fixed_emission_applied_charge(2) = 3.0_dp
  batch2%fixed_emission_weight_scale(2) = 3.0_dp
  batch2%fixed_escape_target_charge(2) = -2.0_dp
  batch2%fixed_escape_applied_charge(2) = -2.0_dp
  batch2%fixed_escape_correction(2) = -0.5_dp
  call accumulate_charge_ledger(cumulative, ledger)
  call accumulate_charge_ledger(cumulative, batch2)
  call assert_close_dp( &
    cumulative%fixed_absorbed_weight_scale(1), 1.0_dp, 1.0d-14, &
    'cumulative fixed absorbed scale must derive from cumulative raw charge' &
    )
  call assert_close_dp( &
    cumulative%fixed_emission_weight_scale(2), 1.0_dp, 1.0d-14, &
    'cumulative fixed emission scale must derive from cumulative raw charge' &
    )
  call assert_close_dp( &
    cumulative%fixed_escape_target_charge(2), -3.0_dp, 1.0d-14, &
    'cumulative fixed escape target mismatch' &
    )
  call assert_close_dp( &
    cumulative%fixed_escape_applied_charge(2), -3.0_dp, 1.0d-14, &
    'cumulative fixed escape applied mismatch' &
    )
  call assert_close_dp( &
    cumulative%fixed_escape_correction(2), -0.75_dp, 1.0d-14, &
    'cumulative fixed escape correction mismatch' &
    )
  call test_end()

  call test_begin('fixed_current_correction_closes_external_surface_source')
  call ledger%reset(1_i32)
  ledger%surface_charge_after = -5.0_dp
  ledger%injected_from_remote(1) = -2.0_dp
  ledger%absorbed_on_surface(1) = -2.0_dp
  ledger%fixed_absorbed_target_charge(1) = -5.0_dp
  ledger%fixed_absorbed_weight_scale(1) = 2.5_dp
  ledger%fixed_current_correction(1) = -3.0_dp
  call assert_close_dp(ledger%residual(), 0.0_dp, 1.0d-14, 'fixed-current residual mismatch')
  call test_end()

  call test_begin('remote_injection_absorption')
  call ledger%reset(1_i32)
  ledger%surface_charge_before = 0.0_dp
  ledger%surface_charge_after = -3.0_dp
  ledger%injected_from_remote(1) = -3.0_dp
  ledger%absorbed_on_surface(1) = -3.0_dp
  call assert_close_dp(ledger%residual(), 0.0_dp, 1.0d-14, 'injection/absorption residual mismatch')
  call test_end()

  call test_begin('photoelectron_escape')
  call ledger%reset(1_i32)
  ledger%surface_charge_before = 0.0_dp
  ledger%surface_charge_after = 2.0_dp
  ledger%emitted_from_surface(2) = -2.0_dp
  ledger%escaped_to_infinity(2) = -2.0_dp
  call assert_close_dp(ledger%residual(), 0.0_dp, 1.0d-14, 'photoelectron escape residual mismatch')
  call test_end()

  call test_begin('persistent_local_stock')
  call ledger%reset(1_i32)
  ledger%local_flight_charge_before = 0.0_dp
  ledger%local_flight_charge_after = -4.0_dp
  ledger%injected_from_remote(1) = -4.0_dp
  call assert_close_dp(ledger%residual(), 0.0_dp, 1.0d-14, 'persistent local stock residual mismatch')
  call test_end()

  call test_begin('discard_is_separate_from_residual')
  call ledger%reset(1_i32)
  ledger%injected_from_remote(1) = -1.0_dp
  ledger%discarded_unresolved(1) = -1.0_dp
  ledger%discarded_unresolved_count(1) = 1_i64
  call assert_close_dp(ledger%residual(), 0.0_dp, 1.0d-14, 'discarded charge residual mismatch')
  call assert_close_dp(ledger%discarded_unresolved_abs(), 1.0_dp, 1.0d-14, 'discarded absolute charge mismatch')
  call assert_true(ledger%has_unresolved_discard(0.5_dp, 0_i64), 'discard threshold must be independent of residual')
  call test_end()

  call test_begin('neutral_return_preserves_raw_unresolved_diagnostic')
  call ledger%reset(1_i32)
  ledger%surface_charge_before = 0.0_dp
  ledger%surface_charge_after = 0.0_dp
  ledger%emitted_from_surface(2) = -10.0_dp
  ledger%absorbed_on_surface(2) = -9.0_dp
  ledger%discarded_unresolved(2) = -1.0_dp
  ledger%neutral_return_correction(2) = -1.0_dp
  ledger%neutral_return_weight_scale(2) = 10.0_dp/9.0_dp
  ledger%neutral_return_unresolved_fraction(2) = 0.1_dp
  call assert_close_dp(ledger%residual(), 0.0_dp, 1.0d-14, 'neutral-return residual mismatch')
  call assert_close_dp( &
    ledger%discarded_unresolved_abs(), 1.0_dp, 1.0d-14, &
    'neutral-return must preserve the raw unresolved diagnostic' &
    )
  call test_end()

  call test_begin('cumulative_batches')
  call ledger%reset(1_i32)
  ledger%surface_charge_before = 0.0_dp
  ledger%surface_charge_after = -2.0_dp
  ledger%injected_from_remote(1) = -2.0_dp
  ledger%absorbed_on_surface(1) = -2.0_dp
  call batch2%init(2_i32)
  call batch2%reset(2_i32)
  batch2%surface_charge_before = -2.0_dp
  batch2%surface_charge_after = -1.0_dp
  batch2%emitted_from_surface(2) = -1.0_dp
  batch2%escaped_to_infinity(2) = -1.0_dp
  call accumulate_charge_ledger(cumulative, ledger)
  call accumulate_charge_ledger(cumulative, batch2)
  call assert_equal_i32(cumulative%batch_count, 2_i32, 'cumulative batch count mismatch')
  call assert_close_dp(cumulative%surface_charge_before, 0.0_dp, 1.0d-14, 'cumulative surface before mismatch')
  call assert_close_dp(cumulative%surface_charge_after, -1.0_dp, 1.0d-14, 'cumulative surface after mismatch')
  call assert_close_dp(cumulative%residual(), 0.0_dp, 1.0d-14, 'cumulative residual mismatch')
  call test_end()

  call test_begin('cumulative_neutral_return_diagnostics')
  cumulative = charge_ledger_type()
  call ledger%reset(1_i32)
  ledger%emitted_from_surface(2) = -8.0_dp
  ledger%absorbed_on_surface(2) = -8.0_dp
  call batch2%reset(2_i32)
  batch2%emitted_from_surface(2) = -2.0_dp
  batch2%absorbed_on_surface(2) = -1.0_dp
  batch2%discarded_unresolved(2) = -1.0_dp
  batch2%neutral_return_correction(2) = -1.0_dp
  call accumulate_charge_ledger(cumulative, ledger)
  call accumulate_charge_ledger(cumulative, batch2)
  call assert_close_dp( &
    cumulative%neutral_return_correction(2), -1.0_dp, 1.0d-14, &
    'cumulative neutral-return correction mismatch' &
    )
  call assert_close_dp( &
    cumulative%neutral_return_weight_scale(2), 10.0_dp/9.0_dp, 1.0d-14, &
    'cumulative neutral-return scale must derive from cumulative raw charge' &
    )
  call assert_close_dp( &
    cumulative%neutral_return_unresolved_fraction(2), 0.1_dp, 1.0d-14, &
    'cumulative neutral-return fraction must derive from cumulative raw charge' &
    )
  call assert_close_dp(cumulative%residual(), 0.0_dp, 1.0d-14, 'cumulative neutral-return residual mismatch')
  call test_end()

  call test_summary()
end program test_charge_ledger

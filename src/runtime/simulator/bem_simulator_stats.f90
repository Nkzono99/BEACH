!> `bem_simulator` のバッチ集計・統計更新処理を実装する submodule。
submodule(bem_simulator) bem_simulator_stats
  implicit none
contains

  !> バッチ内粒子の吸着/脱出/生存内訳をカウントする。
  module procedure count_batch_outcomes
  integer(i32) :: i

  batch_counts = 0_i64
  soft_discarded_abs_charge = 0.0_dp
  batch_counts(1) = int(pcls_batch%n, i64)
  do i = 1, pcls_batch%n
    if (absorbed_flag(i)) then
      batch_counts(2) = batch_counts(2) + 1_i64
    else if (escaped_boundary_flag(i)) then
      batch_counts(3) = batch_counts(3) + 1_i64
      batch_counts(4) = batch_counts(4) + 1_i64
    else if (soft_discarded_boundary_flag(i)) then
      batch_counts(3) = batch_counts(3) + 1_i64
      batch_counts(6) = batch_counts(6) + 1_i64
      soft_discarded_abs_charge = soft_discarded_abs_charge + abs(pcls_batch%q(i)*pcls_batch%w(i))
    else if (pcls_batch%alive(i)) then
      batch_counts(3) = batch_counts(3) + 1_i64
      batch_counts(5) = batch_counts(5) + 1_i64
    end if
  end do
  end procedure count_batch_outcomes

  !> バッチ単位の集計値を累積統計 `sim_stats` に加算する。
  module procedure accumulate_batch_stats
  stats%batches = stats%batches + 1_i32
  stats%last_rel_change = rel
  stats%processed_particles = checked_add_nonnegative_i64( &
                              stats%processed_particles, batch_counts(1), 'processed_particles' &
                              )
  stats%absorbed = checked_add_nonnegative_i64(stats%absorbed, batch_counts(2), 'absorbed')
  stats%escaped = checked_add_nonnegative_i64(stats%escaped, batch_counts(3), 'escaped')
  stats%escaped_boundary = checked_add_nonnegative_i64( &
                           stats%escaped_boundary, batch_counts(4), 'escaped_boundary' &
                           )
  stats%survived_max_step = checked_add_nonnegative_i64( &
                            stats%survived_max_step, batch_counts(5), 'survived_max_step' &
                            )
  stats%multiple_box_events_retry_attempted = checked_add_nonnegative_i64( &
                                              stats%multiple_box_events_retry_attempted, retry_counts(1), &
                                              'multiple_box_events_retry_attempted' &
                                              )
  stats%multiple_box_events_retry_resolved = checked_add_nonnegative_i64( &
                                             stats%multiple_box_events_retry_resolved, retry_counts(2), &
                                             'multiple_box_events_retry_resolved' &
                                             )
  stats%multiple_box_events_soft_discarded = checked_add_nonnegative_i64( &
                                             stats%multiple_box_events_soft_discarded, batch_counts(6), &
                                             'multiple_box_events_soft_discarded' &
                                             )
  stats%multiple_box_events_soft_discarded_abs_charge = checked_add_nonnegative_real( &
                                                        stats%multiple_box_events_soft_discarded_abs_charge, &
                                                        soft_discarded_abs_charge, &
                                                        'multiple_box_events_soft_discarded_abs_charge' &
                                                        )
  end procedure accumulate_batch_stats

  !> 非負の64bit累積カウンタを符号反転させずに加算する。
  function checked_add_nonnegative_i64(accumulated, increment, field_name) result(total)
    integer(i64), intent(in) :: accumulated, increment
    character(len=*), intent(in) :: field_name
    integer(i64) :: total

    if (accumulated < 0_i64 .or. increment < 0_i64) then
      error stop 'simulation statistic must be nonnegative: '//trim(field_name)
    end if
    if (accumulated > huge(accumulated) - increment) then
      error stop 'simulation statistic overflow: '//trim(field_name)
    end if
    total = accumulated + increment
  end function checked_add_nonnegative_i64

  function checked_add_nonnegative_real(accumulated, increment, field_name) result(total)
    real(dp), intent(in) :: accumulated, increment
    character(len=*), intent(in) :: field_name
    real(dp) :: total

    if (.not. all(ieee_is_finite([accumulated, increment])) .or. &
        accumulated < 0.0_dp .or. increment < 0.0_dp) then
      error stop 'simulation statistic must be finite and nonnegative: '//trim(field_name)
    end if
    total = accumulated + increment
    if (.not. ieee_is_finite(total)) then
      error stop 'simulation statistic overflow: '//trim(field_name)
    end if
  end function checked_add_nonnegative_real

end submodule bem_simulator_stats

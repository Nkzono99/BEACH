!> `bem_simulator` のバッチ集計・統計更新処理を実装する submodule。
submodule(bem_simulator) bem_simulator_stats
  implicit none
contains

  !> バッチ内粒子の吸着/脱出/生存内訳をカウントする。
  module procedure count_batch_outcomes
  integer(i32) :: i

  batch_counts = 0_i32
  soft_discarded_abs_charge = 0.0_dp
  batch_counts(1) = pcls_batch%n
  do i = 1, pcls_batch%n
    if (queued_outer_flag(i)) then
      cycle
    else if (absorbed_flag(i)) then
      batch_counts(2) = batch_counts(2) + 1_i32
    else if (escaped_boundary_flag(i)) then
      batch_counts(3) = batch_counts(3) + 1_i32
      batch_counts(4) = batch_counts(4) + 1_i32
    else if (soft_discarded_boundary_flag(i)) then
      batch_counts(3) = batch_counts(3) + 1_i32
      batch_counts(6) = batch_counts(6) + 1_i32
      soft_discarded_abs_charge = soft_discarded_abs_charge + abs(pcls_batch%q(i)*pcls_batch%w(i))
    else if (pcls_batch%alive(i)) then
      batch_counts(3) = batch_counts(3) + 1_i32
      batch_counts(5) = batch_counts(5) + 1_i32
    end if
  end do
  end procedure count_batch_outcomes

  !> バッチ単位の集計値を累積統計 `sim_stats` に加算する。
  module procedure accumulate_batch_stats
  stats%batches = stats%batches + 1_i32
  stats%last_rel_change = rel
  stats%processed_particles = stats%processed_particles + batch_counts(1)
  stats%absorbed = stats%absorbed + batch_counts(2)
  stats%escaped = stats%escaped + batch_counts(3)
  stats%escaped_boundary = stats%escaped_boundary + batch_counts(4)
  stats%survived_max_step = stats%survived_max_step + batch_counts(5)
  stats%multiple_box_events_soft_discarded = &
    stats%multiple_box_events_soft_discarded + int(batch_counts(6), i64)
  stats%multiple_box_events_soft_discarded_abs_charge = &
    stats%multiple_box_events_soft_discarded_abs_charge + soft_discarded_abs_charge
  end procedure accumulate_batch_stats

end submodule bem_simulator_stats

submodule (bem_simulator) bem_simulator_stats
  implicit none
contains

  module procedure count_batch_outcomes
    integer(i32) :: i

    batch_counts = 0_i32
    batch_counts(1) = pcls_batch%n
    do i = 1, pcls_batch%n
      if (absorbed_flag(i)) then
        batch_counts(2) = batch_counts(2) + 1_i32
      else if (escaped_boundary_flag(i)) then
        batch_counts(3) = batch_counts(3) + 1_i32
        batch_counts(4) = batch_counts(4) + 1_i32
      else if (pcls_batch%alive(i)) then
        batch_counts(3) = batch_counts(3) + 1_i32
        batch_counts(5) = batch_counts(5) + 1_i32
      end if
    end do
  end procedure count_batch_outcomes

  module procedure accumulate_batch_stats
    stats%batches = stats%batches + 1_i32
    stats%last_rel_change = rel
    stats%processed_particles = stats%processed_particles + batch_counts(1)
    stats%absorbed = stats%absorbed + batch_counts(2)
    stats%escaped = stats%escaped + batch_counts(3)
    stats%escaped_boundary = stats%escaped_boundary + batch_counts(4)
    stats%survived_max_step = stats%survived_max_step + batch_counts(5)
    stats%field_time_s = stats%field_time_s + field_time_s
    stats%push_time_s = stats%push_time_s + push_time_s
    stats%collision_time_s = stats%collision_time_s + collision_time_s
  end procedure accumulate_batch_stats

end submodule bem_simulator_stats

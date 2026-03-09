submodule (bem_simulator) bem_simulator_io
  implicit none
contains

  module procedure print_batch_progress
    print '(a,i0,a,i0,a,es12.4,a)', &
      '---------- batch ', batch_idx, '/', final_batch_idx, ' rel_change=', rel_change, ' ----------'
    flush (output_unit)
  end procedure print_batch_progress

  module procedure maybe_write_history_snapshot
    if (.not. history_enabled) return
    if (mod(stats%batches - 1_i32, hist_stride) /= 0_i32) return
    call write_history_snapshot(hist_unit, stats%batches, stats%processed_particles, rel, q_elem)
  end procedure maybe_write_history_snapshot

  module procedure write_history_snapshot
    integer(i32) :: elem_idx

    do elem_idx = 1, size(q_elem)
      write(unit_id, '(i0,a,i0,a,es24.16,a,i0,a,es24.16)') batch_idx, ',', processed_particles, ',', &
        rel_change, ',', elem_idx, ',', q_elem(elem_idx)
    end do
  end procedure write_history_snapshot

  module procedure wall_time_seconds
    !$ time_s = omp_get_wtime()
    !$ return
    call cpu_time(time_s)
  end procedure wall_time_seconds

end submodule bem_simulator_io

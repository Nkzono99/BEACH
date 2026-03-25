!> `bem_simulator` の進捗表示と履歴出力を実装する submodule。
submodule(bem_simulator) bem_simulator_io
  implicit none
contains

  !> 現在バッチ番号と相対変化量を標準出力へ表示する。
  module procedure print_batch_progress
  print '(a,i0,a,i0,a,es12.4,a)', &
    '---------- batch ', batch_idx, '/', final_batch_idx, ' rel_change=', rel_change, ' ----------'
  flush (output_unit)
  end procedure print_batch_progress

  !> 履歴出力条件を満たすバッチだけ電荷スナップショットを書き出す。
  module procedure maybe_write_history_snapshot
  if (.not. history_enabled) return
  if (mod(stats%batches - 1_i32, hist_stride) /= 0_i32) return
  call write_history_snapshot(hist_unit, stats%batches, stats%processed_particles, rel, q_elem)
  end procedure maybe_write_history_snapshot

  !> 全要素電荷を履歴CSV形式で1バッチ分書き出す。
  module procedure write_history_snapshot
  integer(i32) :: elem_idx

  do elem_idx = 1, size(q_elem)
    write (unit_id, '(i0,a,i0,a,es24.16,a,i0,a,es24.16)') batch_idx, ',', processed_particles, ',', &
      rel_change, ',', elem_idx, ',', q_elem(elem_idx)
  end do
  end procedure write_history_snapshot

  !> 電位履歴出力条件を満たすバッチだけ電位スナップショットを書き出す。
  module procedure maybe_write_potential_history_snapshot
  if (.not. potential_history_enabled) return
  if (mod(stats%batches - 1_i32, hist_stride) /= 0_i32) return
  call field_solver%refresh(mesh)
  call field_solver%compute_mesh_potential(mesh, sim, potential_buf)
  call write_potential_history_snapshot(pot_hist_unit, stats%batches, potential_buf)
  end procedure maybe_write_potential_history_snapshot

  !> 全要素電位を電位履歴CSV形式で1バッチ分書き出す。
  module procedure write_potential_history_snapshot
  integer(i32) :: elem_idx

  do elem_idx = 1, size(potential_v)
    write (unit_id, '(i0,a,i0,a,es24.16)') batch_idx, ',', elem_idx, ',', potential_v(elem_idx)
  end do
  end procedure write_potential_history_snapshot

end submodule bem_simulator_io

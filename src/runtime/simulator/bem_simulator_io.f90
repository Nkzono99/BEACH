!> `bem_simulator` の進捗表示と履歴出力を実装する submodule。
submodule(bem_simulator) bem_simulator_io
  use bem_app_config_runtime, only: compute_z_high_box_potential_statistics
  use bem_output_writer, only: write_top_reference_history_snapshot
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
  real(dp) :: top_phi_mean, top_phi_std, top_phi_min, top_phi_max

  if (.not. potential_history_enabled) return
  if (mod(stats%batches - 1_i32, hist_stride) /= 0_i32) return
  call snapshot%refresh(mesh, update_outer=.false.)
  call snapshot%compute_mesh_potential(mesh, sim, potential_buf)
  call write_potential_history_snapshot(pot_hist_unit, stats%batches, potential_buf)
  if (top_reference_history_enabled) then
    call compute_z_high_box_potential_statistics( &
      mesh, sim, snapshot, top_phi_mean, top_phi_std, top_phi_min, top_phi_max &
      )
    call write_top_reference_history_snapshot( &
      top_reference_history_unit, stats%batches, stats%simulated_time, sim%box_max(3), &
      sim%injection_face_phi_grid_n, top_phi_mean, top_phi_std, top_phi_min, top_phi_max &
      )
  end if
  end procedure maybe_write_potential_history_snapshot

  !> 全要素電位を電位履歴CSV形式で1バッチ分書き出す。
  module procedure write_potential_history_snapshot
  integer(i32) :: elem_idx

  do elem_idx = 1, size(potential_v)
    write (unit_id, '(i0,a,i0,a,es24.16)') batch_idx, ',', elem_idx, ',', potential_v(elem_idx)
  end do
  end procedure write_potential_history_snapshot

end submodule bem_simulator_io

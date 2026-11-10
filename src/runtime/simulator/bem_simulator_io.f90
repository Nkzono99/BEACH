!> `bem_simulator` の進捗表示と履歴出力を実装する submodule。
submodule(bem_simulator) bem_simulator_io
  use bem_app_config_runtime, only: compute_wall_box_potential_statistics
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
  ! 最長行: i32 batch(11) + i64 粒子数(20) + es24.16 + i32 要素番号(11) + es24.16 + 区切り 4 = 94 文字。
  character(len=96), allocatable :: rows(:)

  allocate (rows(size(q_elem)))
  !$omp parallel do default(none) schedule(static) shared(rows, q_elem, batch_idx, processed_particles, rel_change)
  do elem_idx = 1, size(q_elem)
    write (rows(elem_idx), '(i0,a,i0,a,es24.16,a,i0,a,es24.16)') batch_idx, ',', processed_particles, ',', &
      rel_change, ',', elem_idx, ',', q_elem(elem_idx)
  end do
  !$omp end parallel do
  call write_rows_in_one_record(unit_id, rows)
  end procedure write_history_snapshot

  !> 要素電位と z-high 参照電位の history を 1 batch 分書き出す。場は呼ぶ前に現在の電荷で更新しておく。
  module procedure write_potential_history_outputs
  real(dp) :: top_phi_mean, top_phi_std, top_phi_min, top_phi_max, switch_phi_rms, switch_ez_rms

  call snapshot%compute_mesh_potential(mesh, sim, potential_buf)
  call write_potential_history_snapshot(pot_hist_unit, stats%batches, potential_buf)
  if (top_reference_history_enabled) then
    call compute_wall_box_potential_statistics( &
      mesh, sim, snapshot, top_phi_mean, top_phi_std, top_phi_min, top_phi_max, &
      cell=cell &
      )
    call snapshot%measure_switch_mismatch(mesh, sim, switch_phi_rms, switch_ez_rms)
    call write_top_reference_history_snapshot( &
      top_reference_history_unit, stats%batches, stats%simulated_time, cell%wall_coordinate(), &
      sim%injection_face_phi_grid_n, top_phi_mean, top_phi_std, top_phi_min, top_phi_max, switch_phi_rms, switch_ez_rms &
      )
  end if
  end procedure write_potential_history_outputs

  !> 全要素電位を電位履歴CSV形式で1バッチ分書き出す。
  module procedure write_potential_history_snapshot
  integer(i32) :: elem_idx
  ! 最長行: i32 batch(11) + i32 要素番号(11) + es24.16 + 区切り 2 = 48 文字。
  character(len=48), allocatable :: rows(:)

  allocate (rows(size(potential_v)))
  !$omp parallel do default(none) schedule(static) shared(rows, potential_v, batch_idx)
  do elem_idx = 1, size(potential_v)
    write (rows(elem_idx), '(i0,a,i0,a,es24.16)') batch_idx, ',', elem_idx, ',', potential_v(elem_idx)
  end do
  !$omp end parallel do
  call write_rows_in_one_record(unit_id, rows)
  end procedure write_potential_history_snapshot

  !> 行を改行でつないで 1 回の WRITE で書く。Intel Fortran の既定（無バッファ）では WRITE 1 回が
  !! 書き込み要求 1 回になり、要素ごとに書くと Lustre 上で 1 万行あたり 0.1 秒以上かかる。
  !! ファイルの内容は 1 行ずつ書いた場合と同じ。
  subroutine write_rows_in_one_record(unit_id, rows)
    integer, intent(in) :: unit_id
    character(len=*), intent(in) :: rows(:)
    character(len=:), allocatable :: text
    integer, allocatable :: row_len(:)
    integer :: row_idx, cursor

    if (size(rows) == 0) return
    allocate (row_len(size(rows)))
    do row_idx = 1, size(rows)
      row_len(row_idx) = len_trim(rows(row_idx))
    end do
    allocate (character(len=sum(row_len) + size(rows) - 1) :: text)
    cursor = 0
    do row_idx = 1, size(rows)
      if (row_idx > 1) then
        cursor = cursor + 1
        text(cursor:cursor) = new_line('a')
      end if
      ! 代入先が短いので、行は末尾の空白を落として入る。
      text(cursor + 1:cursor + row_len(row_idx)) = rows(row_idx)
      cursor = cursor + row_len(row_idx)
    end do
    write (unit_id, '(a)') text
  end subroutine write_rows_in_one_record

end submodule bem_simulator_io

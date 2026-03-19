!> 吸着(insulator)モデルのメインループを実行し、電荷堆積と統計更新を行うモジュール。
module bem_simulator
!$ use omp_lib
  use, intrinsic :: iso_fortran_env, only: output_unit
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_stats, mesh_type, particles_soa, hit_info, injection_state
  use bem_app_config, only: app_config, init_particle_batch_from_config
  use bem_field_solver, only: field_solver_type
  use bem_pusher, only: boris_push
  use bem_collision, only: find_first_hit
  use bem_boundary, only: apply_box_boundary
  use bem_mpi, only: mpi_context, mpi_is_root, mpi_allreduce_sum_real_dp_array, mpi_allreduce_max_real_dp_array, &
                     mpi_allreduce_sum_i32_array
  implicit none
  private

  public :: run_absorption_insulator

  interface
    !> 粒子をバッチ処理し、衝突時は要素へ電荷堆積、非衝突時は脱出として統計を更新する。
    module subroutine run_absorption_insulator(mesh, app, stats, history_unit, history_stride, initial_stats, inject_state, mpi)
      type(mesh_type), intent(inout) :: mesh
      type(app_config), intent(in) :: app
      type(sim_stats), intent(out) :: stats
      integer, intent(in), optional :: history_unit
      integer(i32), intent(in), optional :: history_stride
      type(sim_stats), intent(in), optional :: initial_stats
      type(injection_state), intent(inout), optional :: inject_state
      type(mpi_context), intent(in), optional :: mpi
    end subroutine run_absorption_insulator

    !> 1バッチ分の粒子群と作業配列を初期化する。
    module subroutine prepare_batch_state( &
      mesh, app, stats, local_batch_idx, batch_idx, dq_thread, pcls_batch, escaped_boundary_flag, absorbed_flag, &
      photo_emission_dq, mpi, inject_state &
      )
      type(mesh_type), intent(in) :: mesh
      type(app_config), intent(in) :: app
      type(sim_stats), intent(in) :: stats
      integer(i32), intent(in) :: local_batch_idx
      integer(i32), intent(out) :: batch_idx
      real(dp), intent(inout) :: dq_thread(:, :)
      type(particles_soa), intent(out) :: pcls_batch
      logical, allocatable, intent(out) :: escaped_boundary_flag(:)
      logical, allocatable, intent(out) :: absorbed_flag(:)
      real(dp), intent(out) :: photo_emission_dq(:)
      type(mpi_context), intent(in) :: mpi
      type(injection_state), intent(inout), optional :: inject_state
    end subroutine prepare_batch_state

    !> 1バッチぶんの粒子を前進させ、スレッド別に堆積電荷を集計する。
    module subroutine process_particle_batch( &
      mesh, app, field_solver, pcls_batch, dq_thread, escaped_boundary_flag, absorbed_flag, bfield, &
      field_time_s, push_time_s, collision_time_s &
      )
      type(mesh_type), intent(in) :: mesh
      type(app_config), intent(in) :: app
      type(field_solver_type), intent(inout) :: field_solver
      type(particles_soa), intent(inout) :: pcls_batch
      real(dp), intent(inout) :: dq_thread(:, :)
      logical, intent(inout) :: escaped_boundary_flag(:)
      logical, intent(inout) :: absorbed_flag(:)
      real(dp), intent(in) :: bfield(3)
      real(dp), intent(out) :: field_time_s, push_time_s, collision_time_s
    end subroutine process_particle_batch

    !> スレッド別に集計した電荷差分をメッシュへ反映し、相対変化量を返す。
    module subroutine commit_batch_charge(mesh, q_floor, dq_thread, photo_emission_dq, dq, rel, mpi)
      type(mesh_type), intent(inout) :: mesh
      real(dp), intent(in) :: q_floor
      real(dp), intent(in) :: dq_thread(:, :)
      real(dp), intent(in) :: photo_emission_dq(:)
      real(dp), intent(out) :: dq(:)
      real(dp), intent(out) :: rel
      type(mpi_context), intent(in) :: mpi
    end subroutine commit_batch_charge

    !> 今バッチの粒子処理結果を局所集計する。
    module subroutine count_batch_outcomes(pcls_batch, escaped_boundary_flag, absorbed_flag, batch_counts)
      type(particles_soa), intent(in) :: pcls_batch
      logical, intent(in) :: escaped_boundary_flag(:)
      logical, intent(in) :: absorbed_flag(:)
      integer(i32), intent(out) :: batch_counts(5)
    end subroutine count_batch_outcomes

    !> バッチ完了後の統計値を更新する。
    module subroutine accumulate_batch_stats(stats, batch_counts, rel, field_time_s, push_time_s, collision_time_s)
      type(sim_stats), intent(inout) :: stats
      integer(i32), intent(in) :: batch_counts(5)
      real(dp), intent(in) :: rel
      real(dp), intent(in) :: field_time_s, push_time_s, collision_time_s
    end subroutine accumulate_batch_stats

    !> ルートランクでバッチ完了時の進捗と相対変化量を標準出力へ表示する。
    module subroutine print_batch_progress(batch_idx, final_batch_idx, rel_change)
      integer(i32), intent(in) :: batch_idx
      integer(i32), intent(in) :: final_batch_idx
      real(dp), intent(in) :: rel_change
    end subroutine print_batch_progress

    !> 履歴出力が有効で、指定ストライドを満たす場合のみ電荷履歴を書き出す。
    module subroutine maybe_write_history_snapshot(history_enabled, hist_unit, hist_stride, stats, rel, q_elem)
      logical, intent(in) :: history_enabled
      integer, intent(in) :: hist_unit
      integer(i32), intent(in) :: hist_stride
      type(sim_stats), intent(in) :: stats
      real(dp), intent(in) :: rel
      real(dp), intent(in) :: q_elem(:)
    end subroutine maybe_write_history_snapshot

    !> 現時点の要素電荷を CSV 行群として書き出す。
    module subroutine write_history_snapshot(unit_id, batch_idx, processed_particles, rel_change, q_elem)
      integer, intent(in) :: unit_id
      integer(i32), intent(in) :: batch_idx
      integer(i32), intent(in) :: processed_particles
      real(dp), intent(in) :: rel_change
      real(dp), intent(in) :: q_elem(:)
    end subroutine write_history_snapshot

    !> OpenMP有効時は `omp_get_wtime`、それ以外は `cpu_time` を返す簡易タイマ。
    module function wall_time_seconds() result(time_s)
      real(dp) :: time_s
    end function wall_time_seconds
  end interface

end module bem_simulator

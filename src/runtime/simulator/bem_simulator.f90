!> 吸着(insulator)モデルのメインループを実行し、電荷堆積と統計更新を行うモジュール。
module bem_simulator
!$ use omp_lib
  use, intrinsic :: iso_fortran_env, only: output_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32, i64
  use bem_types, only: sim_stats, mesh_type, particles_soa, injection_state, sim_config, hit_info
  use bem_app_config, only: app_config, init_particle_batch_from_config
  use bem_app_config_runtime, only: particle_source_plan_type, build_particle_source_plan
  use bem_particles, only: append_particles
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type, electrostatic_diagnostics_type, &
                                        electrostatic_restart_state_type
  use bem_outer_coupler, only: outer_coupler_type
  use bem_particle_stepper, only: build_particle_step_candidate, resolve_particle_boundary_candidate, advance_particle_step, &
                                  particle_step_result, particle_step_invalid_boundary, &
                                  particle_step_multiple_box_events, particle_step_unsupported_barrier_corner
  use bem_collision, only: collision_query_grid_stalled, collision_query_image_limit, &
                           collision_query_index_range, collision_query_invalid_segment, collision_query_ok, find_first_hit
  use bem_surface_models, only: apply_surface_model_charge_relaxation
  use bem_charge_ledger, only: charge_ledger_type, accumulate_charge_ledger
  use bem_string_utils, only: lower_ascii
  use bem_interface_types, only: interface_particle_outcome_type, interface_outcome_returned_local, &
                                 interface_outcome_escaped_to_infinity, interface_outcome_queued_outer
  use bem_outer_event_queue, only: outer_event_queue_type, outer_event_record_type, &
                                   outer_event_outcome_return, outer_event_outcome_escape, &
                                   outer_event_queue_global_fingerprint
  use bem_outer_plasma_interface, only: map_outer_particle_linear_debye, map_outer_particle_kinetic_profile
  use bem_outer_plasma_orbit, only: trace_unified_outer_particle
  use bem_outer_plasma_photoelectron, only: photoelectron_histogram_type, photoelectron_histogram_state_type, &
                                            validate_photoelectron_linear_applicability, photoelectron_applicability_ok
  use bem_outer_plasma_kinetic, only: kinetic_outer_plasma_options_type
  use bem_outer_plasma_kinetic_runtime, only: resolve_kinetic_outer_options
  use bem_zhao_steady_start, only: initialize_zhao_floating_steady_start
  use bem_outer_plasma_types, only: outer_plasma_ok, outer_plasma_state_type
  use bem_simulator_workspace, only: simulator_batch_workspace_type
  use bem_mpi, only: mpi_context, mpi_is_root, mpi_allreduce_sum_real_dp_array, mpi_allreduce_sum_i32_array, &
                     mpi_allreduce_sum_real_dp_scalar, mpi_allreduce_sum_i32_scalar, mpi_allreduce_sum_i64_array, &
                     mpi_allreduce_max_real_dp_array, &
                     mpi_select_lowest_rank_i32_values
  implicit none
  private

  public :: run_absorption_insulator

  interface
    !> 粒子をバッチ処理し、衝突時は要素へ電荷堆積、非衝突時は脱出として統計を更新する。
    module subroutine run_absorption_insulator( &
      mesh, app, stats, history_unit, history_stride, initial_stats, inject_state, mpi, mesh_potential_v, &
      potential_history_unit, charge_ledger, electrostatic_diagnostics, electrostatic_restart_state, photoelectron_state, &
      outer_queue_state &
      )
      type(mesh_type), intent(inout) :: mesh
      type(app_config), intent(in) :: app
      type(sim_stats), intent(out) :: stats
      integer, intent(in), optional :: history_unit
      integer(i32), intent(in), optional :: history_stride
      type(sim_stats), intent(in), optional :: initial_stats
      type(injection_state), intent(inout), optional :: inject_state
      type(mpi_context), intent(in), optional :: mpi
      real(dp), allocatable, intent(out), optional :: mesh_potential_v(:)
      integer, intent(in), optional :: potential_history_unit
      type(charge_ledger_type), intent(inout), optional :: charge_ledger
      type(electrostatic_diagnostics_type), intent(out), optional :: electrostatic_diagnostics
      type(electrostatic_restart_state_type), intent(inout), optional :: electrostatic_restart_state
      type(photoelectron_histogram_state_type), intent(inout), optional :: photoelectron_state
      type(outer_event_queue_type), intent(inout), optional :: outer_queue_state
    end subroutine run_absorption_insulator

    !> 1バッチ分の粒子群と作業配列を初期化する。
    module subroutine prepare_batch_state( &
      mesh, app, source_plan, snapshot, stats, batch_idx, workspace, pcls_batch, mpi, outer_state, inject_state, &
      collision_failure_status, collision_failure_species, &
      collision_failure_ray, collision_failure_bounce &
      )
      type(mesh_type), intent(in) :: mesh
      type(app_config), intent(in) :: app
      type(particle_source_plan_type), intent(in) :: source_plan
      type(electrostatic_snapshot_type), intent(inout) :: snapshot
      type(sim_stats), intent(in) :: stats
      integer(i32), intent(out) :: batch_idx
      type(simulator_batch_workspace_type), intent(inout) :: workspace
      type(particles_soa), intent(out) :: pcls_batch
      type(mpi_context), intent(in) :: mpi
      type(outer_plasma_state_type), intent(in) :: outer_state
      type(injection_state), intent(inout), optional :: inject_state
      integer(i32), intent(out) :: collision_failure_status, collision_failure_species
      integer(i32), intent(out) :: collision_failure_ray, collision_failure_bounce
    end subroutine prepare_batch_state

    !> 1バッチぶんの粒子を前進させ、スレッド別に堆積電荷を集計する。
    module subroutine process_particle_batch( &
      mesh, app, snapshot, pcls_batch, dq_thread, escaped_boundary_flag, absorbed_flag, bfield, batch_idx, mpi_rank, &
      soft_discarded_boundary_flag, queued_outer_flag, outer_event_staging, &
      interface_outward_thread, interface_returned_thread, &
      collision_failure_status, collision_failure_particle, &
      collision_failure_step, collision_failure_x, collision_failure_v, interface_tau_max_thread, &
      interface_frozen_ratio_max_thread, interface_energy_error_max_thread, photoelectron_histogram_thread &
      )
      type(mesh_type), intent(in) :: mesh
      type(app_config), intent(in) :: app
      type(electrostatic_snapshot_type), intent(inout) :: snapshot
      type(particles_soa), intent(inout) :: pcls_batch
      real(dp), intent(inout) :: dq_thread(:, :)
      logical, intent(inout) :: escaped_boundary_flag(:)
      logical, intent(inout) :: absorbed_flag(:)
      logical, intent(inout) :: soft_discarded_boundary_flag(:)
      logical, intent(inout) :: queued_outer_flag(:)
      type(outer_event_record_type), intent(inout) :: outer_event_staging(:)
      real(dp), intent(in) :: bfield(3)
      integer(i32), intent(in) :: batch_idx
      integer(i32), intent(in) :: mpi_rank
      real(dp), intent(out) :: interface_outward_thread(:, :), interface_returned_thread(:, :)
      integer(i32), intent(out) :: collision_failure_status, collision_failure_particle, collision_failure_step
      real(dp), intent(out) :: collision_failure_x(3), collision_failure_v(3)
      real(dp), intent(out) :: interface_tau_max_thread(:), interface_frozen_ratio_max_thread(:)
      real(dp), intent(out) :: interface_energy_error_max_thread(:)
      type(photoelectron_histogram_type), intent(inout), optional :: photoelectron_histogram_thread(:)
    end subroutine process_particle_batch

    !> スレッド別に集計した電荷差分をメッシュへ反映し、相対変化量を返す。
    module subroutine commit_batch_charge( &
      mesh, q_floor, softening, external_e, field_bc_mode, workspace, rel, mpi &
      )
      type(mesh_type), intent(inout) :: mesh
      real(dp), intent(in) :: q_floor
      real(dp), intent(in) :: softening
      real(dp), intent(in) :: external_e(3)
      character(len=*), intent(in) :: field_bc_mode
      type(simulator_batch_workspace_type), intent(inout) :: workspace
      real(dp), intent(out) :: rel
      type(mpi_context), intent(in) :: mpi
    end subroutine commit_batch_charge

    !> 今バッチの粒子処理結果を局所集計する。
    module subroutine count_batch_outcomes( &
      pcls_batch, escaped_boundary_flag, absorbed_flag, soft_discarded_boundary_flag, &
      queued_outer_flag, batch_counts, soft_discarded_abs_charge &
      )
      type(particles_soa), intent(in) :: pcls_batch
      logical, intent(in) :: escaped_boundary_flag(:)
      logical, intent(in) :: absorbed_flag(:)
      logical, intent(in) :: soft_discarded_boundary_flag(:)
      logical, intent(in) :: queued_outer_flag(:)
      integer(i32), intent(out) :: batch_counts(6)
      real(dp), intent(out) :: soft_discarded_abs_charge
    end subroutine count_batch_outcomes

    !> バッチ完了後の統計値を更新する。
    module subroutine accumulate_batch_stats(stats, batch_counts, soft_discarded_abs_charge, rel)
      type(sim_stats), intent(inout) :: stats
      integer(i32), intent(in) :: batch_counts(6)
      real(dp), intent(in) :: soft_discarded_abs_charge
      real(dp), intent(in) :: rel
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
      integer(i64), intent(in) :: processed_particles
      real(dp), intent(in) :: rel_change
      real(dp), intent(in) :: q_elem(:)
    end subroutine write_history_snapshot

    !> 電位履歴出力条件を満たすバッチだけ電位スナップショットを書き出す。
    module subroutine maybe_write_potential_history_snapshot( &
      potential_history_enabled, pot_hist_unit, hist_stride, stats, &
      snapshot, mesh, sim, potential_buf &
      )
      logical, intent(in) :: potential_history_enabled
      integer, intent(in) :: pot_hist_unit
      integer(i32), intent(in) :: hist_stride
      type(sim_stats), intent(in) :: stats
      type(electrostatic_snapshot_type), intent(inout) :: snapshot
      type(mesh_type), intent(inout) :: mesh
      type(sim_config), intent(in) :: sim
      real(dp), intent(inout) :: potential_buf(:)
    end subroutine maybe_write_potential_history_snapshot

    !> 全要素電位を電位履歴CSV形式で1バッチ分書き出す。
    module subroutine write_potential_history_snapshot(unit_id, batch_idx, potential_v)
      integer, intent(in) :: unit_id
      integer(i32), intent(in) :: batch_idx
      real(dp), intent(in) :: potential_v(:)
    end subroutine write_potential_history_snapshot

  end interface

end module bem_simulator

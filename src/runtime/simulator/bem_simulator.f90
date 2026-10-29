!> 吸着(insulator)モデルのメインループを実行し、電荷堆積と統計更新を行う。
module bem_simulator
!$ use omp_lib
  use, intrinsic :: iso_fortran_env, only: output_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32, i64
  use bem_types, only: sim_stats, mesh_type, particles_soa, injection_state, sim_config, particle_domain_config, hit_info
  use bem_app_config, only: app_config, init_particle_batch_from_config
  use bem_physics_config_types, only: field_physics_config, panel_kernel_config, derive_field_panel_config
  use bem_irradiation_cell_contract, only: &
    irradiation_cell_contract_type, resolve_irradiation_cell_contract, irradiation_cell_ok
  use bem_config_helpers, only: resolve_particle_boundaries
  use bem_app_config_runtime, only: particle_source_plan_type, build_particle_source_plan, &
                                    batch_injection_layout_type, build_batch_injection_layout, generate_batch_slot
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type, electrostatic_diagnostics_type
  use bem_particle_stepper, only: build_particle_step_candidate, resolve_particle_boundary_candidate, &
                                  advance_particle_step_upper_panel_fourier, &
                                  particle_step_result
  use bem_status_codes, only: particle_step_invalid_boundary, particle_step_multiple_box_events, &
                              particle_step_ambiguous_open_corner, particle_status_name, collision_query_ok, &
                              particle_track_failure_type, photo_ray_failure_type
  use bem_collision, only: find_first_hit
  use bem_surface_models, only: floating_conductor_solver_type
  use bem_charge_ledger, only: charge_ledger_type, accumulate_charge_ledger
  use bem_string_utils, only: lower_ascii
  use bem_external_boundary_contract, only: external_boundary_contract_type
  use bem_simulator_workspace, only: simulator_batch_workspace_type
  use bem_surface_closure_contract, only: surface_closure_contract_type
  use bem_surface_current_model, only: evaluate_surface_closure
  use bem_run_outputs, only: run_outputs_type
  use bem_output_writer, only: write_matching_plane_history_snapshot, write_fixed_current_history_rows, &
                               write_matching_plane_spectrum_history_rows, write_secondary_emission_history_rows
  use bem_mpi, only: mpi_context, mpi_is_root, mpi_allreduce_sum_real_dp_array, mpi_allreduce_sum_real_dp_scalar, &
                     mpi_shared_counter_type, &
                     mpi_allreduce_sum_i32_scalar, mpi_allreduce_sum_i64_array, &
                     mpi_allreduce_min_i32_scalar, mpi_allreduce_max_i32_scalar, mpi_select_lowest_rank_i32_values
  implicit none
  private

  !> `batch_counts(:)` の各要素の意味。MPI 合計は配列のまま行うので、順序はここだけで決める。
  !! escaped は境界脱出・最大 step 到達・soft discard の合計（`sim_stats%escaped` と同じ定義）。
  integer(i32), parameter :: batch_count_processed = 1_i32
  integer(i32), parameter :: batch_count_absorbed = 2_i32
  integer(i32), parameter :: batch_count_escaped = 3_i32
  integer(i32), parameter :: batch_count_escaped_boundary = 4_i32
  integer(i32), parameter :: batch_count_survived_max_step = 5_i32
  integer(i32), parameter :: batch_count_soft_discarded = 6_i32
  integer(i32), parameter :: batch_count_kinds = 6_i32

  public :: run_absorption_insulator
  !> 固定電流 closure は contract との契約をテストから直接確認する。
  public :: apply_fixed_surface_current_closure

  interface
    module subroutine run_absorption_insulator( &
      mesh, app, stats, initial_stats, inject_state, mpi, charge_ledger, electrostatic_diagnostics, outputs &
      )
      type(mesh_type), intent(inout) :: mesh
      type(app_config), intent(in) :: app
      type(sim_stats), intent(out) :: stats
      type(sim_stats), intent(in), optional :: initial_stats
      type(injection_state), intent(inout), optional :: inject_state
      type(mpi_context), intent(in), optional :: mpi
      type(charge_ledger_type), intent(inout), optional :: charge_ledger
      type(electrostatic_diagnostics_type), intent(out), optional :: electrostatic_diagnostics
      !> history の出力先と終了時のメッシュ電位。省略時は何も書かない。
      type(run_outputs_type), intent(inout), optional :: outputs
    end subroutine run_absorption_insulator

    module subroutine prepare_batch_state( &
      mesh, app, source_plan, snapshot, stats, batch_idx, workspace, pcls_batch, injection_layout, mpi, inject_state, &
      photo_failure &
      )
      type(mesh_type), intent(in) :: mesh
      type(app_config), intent(in) :: app
      type(particle_source_plan_type), intent(in) :: source_plan
      type(electrostatic_snapshot_type), intent(inout) :: snapshot
      type(sim_stats), intent(in) :: stats
      integer(i32), intent(out) :: batch_idx
      type(simulator_batch_workspace_type), intent(inout) :: workspace
      type(particles_soa), intent(out) :: pcls_batch
      !> `shared` のときだけ構築される global slot layout。`static` では `ready=.false.`。
      type(batch_injection_layout_type), intent(out) :: injection_layout
      type(mpi_context), intent(in) :: mpi
      type(injection_state), intent(inout), optional :: inject_state
      !> この rank の photo raycast 注入で最初に失敗した ray。
      type(photo_ray_failure_type), intent(out) :: photo_failure
    end subroutine prepare_batch_state

    module subroutine process_particle_batch( &
      mesh, app, boundary_contract, current_model, snapshot, pcls_batch, workspace, bfield, batch_idx, mpi_rank, &
      actual_team_size, track_failure, retry_counts, substep_counts, &
      cell, injection_layout, pool_counter, pool_claims, pool_claim_seconds &
      )
      type(irradiation_cell_contract_type), intent(in) :: cell
      type(mesh_type), intent(in) :: mesh
      type(app_config), intent(in) :: app
      type(external_boundary_contract_type), intent(in) :: boundary_contract
      type(surface_closure_contract_type), intent(in) :: current_model
      type(electrostatic_snapshot_type), intent(inout) :: snapshot
      type(particles_soa), intent(inout) :: pcls_batch
      type(simulator_batch_workspace_type), intent(inout) :: workspace
      real(dp), intent(in) :: bfield(3)
      integer(i32), intent(in) :: batch_idx, mpi_rank
      integer(i32), intent(out) :: actual_team_size
      !> この rank で最初に失敗した粒子追跡。
      type(particle_track_failure_type), intent(out) :: track_failure
      integer(i64), intent(out) :: retry_counts(2)
      !> 半分にした dt step 数、その sub-step 総数、半分化の上限に達した step 数、dt より長く進めた push 数。
      integer(i64), intent(out) :: substep_counts(4)
      !> `shared` の slot layout と claim カウンタ。`static` では layout は未使用。
      type(batch_injection_layout_type), intent(in) :: injection_layout
      type(mpi_shared_counter_type), intent(inout) :: pool_counter
      !> この rank の claim 回数と claim の critical 区間で費やした秒数。
      integer(i64), intent(out) :: pool_claims
      real(dp), intent(out) :: pool_claim_seconds
    end subroutine process_particle_batch

    module subroutine commit_batch_charge(mesh, sim, snapshot, conductor, workspace, rel, mpi)
      type(mesh_type), intent(inout) :: mesh
      type(sim_config), intent(in) :: sim
      type(electrostatic_snapshot_type), intent(inout) :: snapshot
      type(floating_conductor_solver_type), intent(in) :: conductor
      type(simulator_batch_workspace_type), intent(inout) :: workspace
      real(dp), intent(out) :: rel
      type(mpi_context), intent(in) :: mpi
    end subroutine commit_batch_charge

    module subroutine prepare_adaptive_charge_candidate(mesh, sim, snapshot, conductor, workspace, mpi)
      type(mesh_type), intent(inout) :: mesh
      type(sim_config), intent(in) :: sim
      type(electrostatic_snapshot_type), intent(inout) :: snapshot
      type(floating_conductor_solver_type), intent(in) :: conductor
      type(simulator_batch_workspace_type), intent(inout) :: workspace
      type(mpi_context), intent(in) :: mpi
    end subroutine prepare_adaptive_charge_candidate

    module subroutine apply_neutral_return_surface_closure(app, pcls_batch, fresh_particle_count, workspace, mpi)
      type(app_config), intent(in) :: app
      type(particles_soa), intent(in) :: pcls_batch
      integer(i32), intent(in) :: fresh_particle_count
      type(simulator_batch_workspace_type), intent(inout) :: workspace
      type(mpi_context), intent(in) :: mpi
    end subroutine apply_neutral_return_surface_closure

    module subroutine apply_fixed_surface_current_closure( &
      app, current_model, pcls_batch, fresh_particle_count, workspace, mpi &
      )
      type(app_config), intent(in) :: app
      type(surface_closure_contract_type), intent(in) :: current_model
      type(particles_soa), intent(in) :: pcls_batch
      integer(i32), intent(in) :: fresh_particle_count
      type(simulator_batch_workspace_type), intent(inout) :: workspace
      type(mpi_context), intent(in) :: mpi
    end subroutine apply_fixed_surface_current_closure

    module subroutine record_batch_initial_charge(app, pcls_batch, fresh_particle_count, ledger)
      type(app_config), intent(in) :: app
      type(particles_soa), intent(in) :: pcls_batch
      integer(i32), intent(in) :: fresh_particle_count
      type(charge_ledger_type), intent(inout) :: ledger
    end subroutine record_batch_initial_charge

    module subroutine record_batch_outcome_charge( &
      pcls_batch, escaped_boundary_flag, absorbed_flag, soft_discarded_boundary_flag, ledger &
      )
      type(particles_soa), intent(in) :: pcls_batch
      logical, intent(in) :: escaped_boundary_flag(:), absorbed_flag(:), soft_discarded_boundary_flag(:)
      type(charge_ledger_type), intent(inout) :: ledger
    end subroutine record_batch_outcome_charge

    module subroutine reduce_charge_ledger_fluxes(ledger, mpi, workspace)
      type(charge_ledger_type), intent(inout) :: ledger
      type(mpi_context), intent(in) :: mpi
      type(simulator_batch_workspace_type), intent(inout) :: workspace
    end subroutine reduce_charge_ledger_fluxes

    module subroutine count_batch_outcomes( &
      pcls_batch, escaped_boundary_flag, absorbed_flag, soft_discarded_boundary_flag, &
      batch_counts, soft_discarded_abs_charge &
      )
      type(particles_soa), intent(in) :: pcls_batch
      logical, intent(in) :: escaped_boundary_flag(:), absorbed_flag(:), soft_discarded_boundary_flag(:)
      integer(i64), intent(out) :: batch_counts(batch_count_kinds)
      real(dp), intent(out) :: soft_discarded_abs_charge
    end subroutine count_batch_outcomes

    module subroutine accumulate_batch_stats(stats, batch_counts, soft_discarded_abs_charge, retry_counts, &
                                             substep_counts, rel)
      type(sim_stats), intent(inout) :: stats
      integer(i64), intent(in) :: batch_counts(batch_count_kinds)
      integer(i64), intent(in) :: retry_counts(2)
      integer(i64), intent(in) :: substep_counts(4)
      real(dp), intent(in) :: soft_discarded_abs_charge, rel
    end subroutine accumulate_batch_stats

    module subroutine accumulate_secondary_stats(stats, counts, values)
      type(sim_stats), intent(inout) :: stats
      integer(i64), intent(in) :: counts(:, :)
      real(dp), intent(in) :: values(:, :)
    end subroutine accumulate_secondary_stats

    module subroutine print_batch_progress(batch_idx, final_batch_idx, rel_change)
      integer(i32), intent(in) :: batch_idx, final_batch_idx
      real(dp), intent(in) :: rel_change
    end subroutine print_batch_progress

    module subroutine maybe_write_history_snapshot(history_enabled, hist_unit, hist_stride, stats, rel, q_elem)
      logical, intent(in) :: history_enabled
      integer, intent(in) :: hist_unit
      integer(i32), intent(in) :: hist_stride
      type(sim_stats), intent(in) :: stats
      real(dp), intent(in) :: rel, q_elem(:)
    end subroutine maybe_write_history_snapshot

    module subroutine write_history_snapshot(unit_id, batch_idx, processed_particles, rel_change, q_elem)
      integer, intent(in) :: unit_id
      integer(i32), intent(in) :: batch_idx
      integer(i64), intent(in) :: processed_particles
      real(dp), intent(in) :: rel_change, q_elem(:)
    end subroutine write_history_snapshot

    module subroutine maybe_write_potential_history_snapshot( &
      potential_history_enabled, pot_hist_unit, hist_stride, stats, snapshot, mesh, sim, potential_buf, &
      top_reference_history_enabled, top_reference_history_unit, &
      cell &
      )
      type(irradiation_cell_contract_type), intent(in) :: cell
      logical, intent(in) :: potential_history_enabled
      integer, intent(in) :: pot_hist_unit
      integer(i32), intent(in) :: hist_stride
      type(sim_stats), intent(in) :: stats
      type(electrostatic_snapshot_type), intent(inout) :: snapshot
      type(mesh_type), intent(inout) :: mesh
      type(sim_config), intent(in) :: sim
      real(dp), intent(inout) :: potential_buf(:)
      logical, intent(in) :: top_reference_history_enabled
      integer, intent(in) :: top_reference_history_unit
    end subroutine maybe_write_potential_history_snapshot

    module subroutine write_potential_history_snapshot(unit_id, batch_idx, potential_v)
      integer, intent(in) :: unit_id
      integer(i32), intent(in) :: batch_idx
      real(dp), intent(in) :: potential_v(:)
    end subroutine write_potential_history_snapshot
  end interface

end module bem_simulator

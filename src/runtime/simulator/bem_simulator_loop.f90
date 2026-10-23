!> `bem_simulator` の主ループと粒子処理計算を実装する submodule。
submodule(bem_simulator) bem_simulator_loop
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use bem_app_config_runtime, only: compute_wall_box_potential_statistics
  use bem_zhao_outflow_refresh, only: zhao_outflow_refresh_type
  use bem_types, only: adaptive_nonzero_mode_stats
  use bem_simulator_workspace, only: fixed_charge_absorbed, fixed_charge_emission, fixed_charge_escaped
  use bem_periodic_zero_mode_plan, only: periodic_zero_mode_state_type
  use bem_performance_profile, only: perf_region_batch_total, perf_region_begin, perf_region_commit_charge, &
                                     perf_region_count_outcomes, perf_region_end, perf_region_field_refresh, &
                                     perf_region_field_solver_init, perf_region_history_write, perf_region_mpi_reduce, &
                                     perf_region_particle_batch, perf_region_prepare_batch, perf_region_simulation_total, &
                                     perf_region_stats_update
  use bem_mpi, only: mpi_allreduce_max_real_dp_array
  use bem_periodic_checkpoint, only: maybe_write_periodic_checkpoint
  use bem_charge_ledger, only: finite_charge_sum
  implicit none
  integer(i32), parameter :: adaptive_max_halvings = 24_i32
contains

  module procedure run_absorption_insulator
  integer(i32) :: batch_idx, final_batch_idx, batch_count_this_run, local_batch_idx, nth, hist_stride
  integer(i32) :: particle_team_size
  type(particle_track_failure_type) :: track_failure
  type(photo_ray_failure_type) :: photo_failure
  integer(i32) :: trial_halvings, fresh_particle_count
  integer :: hist_unit, pot_hist_unit, top_ref_hist_unit, matching_hist_unit, fixed_hist_unit, spectrum_hist_unit
  integer :: secondary_hist_unit
  logical :: history_enabled, potential_history_enabled, top_reference_history_enabled
  logical :: ledger_enabled, adaptive_nonzero_mode, omp_dynamic_before
  logical :: replay_active, matching_history_enabled, outflow_refreshed, fixed_history_enabled
  logical :: spectrum_history_enabled, secondary_history_enabled
  integer(i32) :: fixed_window_batches
  ! 固定電流historyの窓合計。(channel, species)、channelは1=吸収、2=放出、3=脱出。
  real(dp), allocatable :: fixed_window_tracked(:, :), fixed_window_target(:, :)
  logical, allocatable :: fixed_window_active(:, :)
  real(dp), allocatable :: potential_buf(:), injection_residual_before(:), boundary_injection_residual_before(:, :)
  integer(i64) :: batch_counts(6), batch_retry_counts(2), batch_substep_counts(4)
  real(dp) :: bfield(3), rel, t0, sim_t0, batch_t0, batch_soft_discarded_abs_charge
  real(dp) :: trial_batch_duration, adaptive_potential_step
  type(particles_soa) :: pcls_batch
  type(mpi_context) :: mpi_ctx
  type(electrostatic_snapshot_type) :: snapshot
  type(floating_conductor_solver_type) :: conductor_solver
  type(electrostatic_diagnostics_type) :: committed_snapshot_diagnostics
  type(periodic_zero_mode_state_type) :: committed_zero_state
  type(charge_ledger_type) :: batch_ledger
  type(zhao_outflow_refresh_type) :: outflow_refresh
  type(simulator_batch_workspace_type) :: workspace
  type(particle_source_plan_type) :: source_plan
  type(batch_injection_layout_type) :: injection_layout
  type(mpi_shared_counter_type) :: pool_counter
  logical :: shared_pool
  integer(i64) :: batch_pool_claims
  real(dp) :: batch_pool_claim_seconds
  type(external_boundary_contract_type) :: boundary_contract
  type(surface_closure_contract_type) :: surface_closure
  type(irradiation_cell_contract_type) :: cell
  type(app_config) :: trial_app
  integer(i64) :: secondary_counts(3, app%n_particle_species)
  real(dp) :: secondary_values(2, app%n_particle_species)
  type(sim_stats) :: stats_candidate

  call initialize_run()
  call resolve_history_outputs()
  call initialize_fields()

  do local_batch_idx = 1_i32, batch_count_this_run
    call perf_region_begin(perf_region_batch_total, batch_t0)
    batch_idx = stats%batches + 1_i32
    call perf_region_begin(perf_region_field_refresh, t0)
    call snapshot%refresh(mesh, mpi_ctx)
    call perf_region_end(perf_region_field_refresh, t0)

    call run_accepted_trial()
    call reduce_batch_outcomes()
    call build_stats_candidate()
    if (ledger_enabled) call record_batch_ledger()
    call commit_accepted_batch()

    call perf_region_begin(perf_region_history_write, t0)
    if (mpi_is_root(mpi_ctx)) call write_batch_history()
    call perf_region_end(perf_region_history_write, t0)
    call maybe_write_periodic_checkpoint(app, mesh, stats, inject_state, mpi_ctx, charge_ledger)
    call perf_region_end(perf_region_batch_total, batch_t0)
  end do
  if (shared_pool) call pool_counter%destroy()
  call perf_region_end(perf_region_simulation_total, sim_t0)
!$ if (replay_active) call omp_set_dynamic(omp_dynamic_before)

  call finish_run_outputs()

contains

  !> 統計・charge ledger・境界契約を初期化し、replay する trial の OpenMP team size を確定する。
  subroutine initialize_run()
    integer(i32) :: team_size_min, team_size_max, boundary_status
    character(len=256) :: boundary_message

    stats = sim_stats()
    if (present(initial_stats)) stats = initial_stats
    mpi_ctx = mpi_context()
    if (present(mpi)) mpi_ctx = mpi
    ledger_enabled = present(charge_ledger)
    if (ledger_enabled) then
      if (.not. allocated(charge_ledger%injected_from_remote)) then
        call charge_ledger%init(app%n_particle_species)
      else if (charge_ledger%nspecies /= app%n_particle_species) then
        error stop 'charge ledger species count does not match app config.'
      end if
      call batch_ledger%init(app%n_particle_species)
    end if

    call resolve_external_boundary_contract( &
      app%sim%reservoir_potential_model, &
      boundary_contract, boundary_status, boundary_message &
      )
    if (boundary_status /= external_boundary_ok) error stop trim(boundary_message)
    call resolve_irradiation_cell_contract(app, cell, boundary_status, boundary_message)
    if (boundary_status /= irradiation_cell_ok) error stop trim(boundary_message)
    adaptive_nonzero_mode = app%periodic2%max_nonzero_mode_potential_step > 0.0_dp
    replay_active = adaptive_nonzero_mode
    omp_dynamic_before = .false.
!$  omp_dynamic_before = omp_get_dynamic()
    nth = 1_i32
    if (replay_active) then
!$    call omp_set_dynamic(.false.)
      !$omp parallel default(none) shared(nth)
      !$omp single
!$    nth = int(omp_get_num_threads(), i32)
      !$omp end single
      !$omp end parallel

      team_size_min = nth
      team_size_max = nth
      call mpi_allreduce_min_i32_scalar(mpi_ctx, team_size_min)
      call mpi_allreduce_max_i32_scalar(mpi_ctx, team_size_max)
      if (team_size_min /= team_size_max) then
        if (mpi_is_root(mpi_ctx)) then
          write (error_unit, '(a,i0,a,i0)') &
            'replayed-trial OpenMP team-size mismatch across MPI ranks: min=', team_size_min, ' max=', team_size_max
          flush (error_unit)
        end if
        error stop 'replayed-trial OpenMP team size must match across MPI ranks.'
      end if

      ! Rejected trials within this run still use one fixed team size.  A restart may
      ! use a different team size because checkpoint compatibility is numerical, not
      ! bitwise-replay compatibility with the previous process.
      if (adaptive_nonzero_mode) stats%adaptive_nonzero_mode%omp_threads = nth
    else
!$    nth = max(1_i32, int(omp_get_max_threads(), i32))
    end if
    if (.not. adaptive_nonzero_mode) then
      stats%adaptive_nonzero_mode = adaptive_nonzero_mode_stats()
    end if
    call validate_soft_discard_initial_state(app%sim, stats)
    call enforce_soft_discard_limits(app%sim, stats, stats%batches, 0_i64, 0_i64, 0.0_dp, mpi_ctx)
    call workspace%init(mesh%nelem, app%n_particle_species, nth, candidate_charge_enabled=adaptive_nonzero_mode)
    call evaluate_surface_closure(app, surface_closure, cell)
  end subroutine initialize_run

  !> 出力先の unit と間隔を取り出す。unit=-1 の history は書かない。間隔が 0 以下なら matching-plane・spectrum・
  !! 固定電流 history を止め、電荷・要素電位 history は毎 batch 扱いにする。
  subroutine resolve_history_outputs()
    integer(i32) :: history_stride

    hist_unit = -1
    pot_hist_unit = -1
    top_ref_hist_unit = -1
    matching_hist_unit = -1
    spectrum_hist_unit = -1
    fixed_hist_unit = -1
    secondary_hist_unit = -1
    history_stride = 1_i32
    if (present(outputs)) then
      hist_unit = outputs%history_unit
      pot_hist_unit = outputs%potential_history_unit
      top_ref_hist_unit = outputs%top_reference_history_unit
      matching_hist_unit = outputs%matching_plane_history_unit
      spectrum_hist_unit = outputs%matching_plane_spectrum_history_unit
      fixed_hist_unit = outputs%fixed_current_history_unit
      secondary_hist_unit = outputs%secondary_emission_history_unit
      history_stride = outputs%history_stride
      if (allocated(outputs%mesh_potential_v)) deallocate (outputs%mesh_potential_v)
    end if
    history_enabled = hist_unit /= -1
    potential_history_enabled = pot_hist_unit /= -1
    top_reference_history_enabled = top_ref_hist_unit /= -1
    matching_history_enabled = matching_hist_unit /= -1 .and. history_stride > 0_i32
    spectrum_history_enabled = spectrum_hist_unit /= -1 .and. history_stride > 0_i32
    fixed_history_enabled = fixed_hist_unit /= -1 .and. history_stride > 0_i32
    secondary_history_enabled = secondary_hist_unit /= -1
    hist_stride = max(1_i32, history_stride)
    allocate ( &
      fixed_window_tracked(3, app%n_particle_species), fixed_window_target(3, app%n_particle_species), &
      fixed_window_active(3, app%n_particle_species) &
      )
    call reset_fixed_current_window()
    if (potential_history_enabled) allocate (potential_buf(mesh%nelem))
  end subroutine resolve_history_outputs

  !> 再開 batch を検査し、場ソルバ・浮遊導体・外部根の更新を初期化して、replay と共有 pool の作業領域を用意する。
  subroutine initialize_fields()
    type(field_physics_config) :: field_config
    type(panel_kernel_config) :: panel_config

    bfield = app%sim%b0
    final_batch_idx = app%sim%batch_count
    if (stats%batches < 0_i32 .or. stats%batches > final_batch_idx) then
      error stop 'checkpoint batch count is outside run.batch_count.'
    end if
    if (.not. adaptive_nonzero_mode .and. stats%batches > 0_i32 .and. stats%simulated_time == 0.0_dp) then
      stats%simulated_time = real(stats%batches, dp)*app%sim%batch_duration
    end if
    batch_count_this_run = final_batch_idx - stats%batches

    call perf_region_begin(perf_region_simulation_total, sim_t0)
    call perf_region_begin(perf_region_field_solver_init, t0)
    call derive_field_panel_config(app%sim, field_config, panel_config)
    call snapshot%init(mesh, app%sim, field_config, app%periodic2, panel_config, mpi=mpi_ctx)
    call conductor_solver%init(mesh, app%sim, snapshot)
    call perf_region_end(perf_region_field_solver_init, t0)
    call outflow_refresh%initialize(app, mesh, stats%matching_plane, surface_closure, cell)
    call apply_surface_plane_gauge(cell, mesh, surface_closure, snapshot, mpi_ctx, .true.)
    call record_matching_fields(stats, snapshot, surface_closure, mesh)

    if (replay_active) then
      if (present(inject_state)) then
        if (allocated(inject_state%macro_residual)) allocate (injection_residual_before(size(inject_state%macro_residual)))
        if (allocated(inject_state%boundary_macro_residual)) then
          allocate ( &
            boundary_injection_residual_before( &
            size(inject_state%boundary_macro_residual, 1), size(inject_state%boundary_macro_residual, 2) &
            ) &
            )
        end if
      end if
    end if
    adaptive_potential_step = 0.0_dp

    shared_pool = trim(app%sim%particle_distribution) == 'shared'
    if (shared_pool) call pool_counter%init(mpi_ctx)
  end subroutine initialize_fields

  !> batch の trial を実行し、adaptive nonzero-mode の電位変化が上限に収まるまで duration を半分にして
  !! 注入残差と場の状態を巻き戻してやり直す。受理した trial の粒子と closure が workspace に残る。
  subroutine run_accepted_trial()
    integer(i32) :: species_idx, particle_team_size_min, particle_team_size_max
    real(dp) :: duration_ratio, adaptive_metric_values(1)
    logical :: trial_accepted

    trial_batch_duration = app%sim%batch_duration
    trial_halvings = 0_i32
    trial_accepted = .false.
    if (replay_active) then
      if (allocated(injection_residual_before)) injection_residual_before = inject_state%macro_residual
      if (allocated(boundary_injection_residual_before)) then
        boundary_injection_residual_before = inject_state%boundary_macro_residual
      end if
      committed_zero_state = snapshot%zero_state
      committed_snapshot_diagnostics = snapshot%diagnostics
    end if

    do while (.not. trial_accepted)
      if (replay_active) then
        if (allocated(injection_residual_before)) inject_state%macro_residual = injection_residual_before
        if (allocated(boundary_injection_residual_before)) then
          inject_state%boundary_macro_residual = boundary_injection_residual_before
        end if
        snapshot%zero_state = committed_zero_state
        snapshot%diagnostics = committed_snapshot_diagnostics
      end if
      trial_app = app
      trial_app%sim%batch_duration = trial_batch_duration
      duration_ratio = trial_batch_duration/app%sim%batch_duration
      do species_idx = 1_i32, trial_app%n_particle_species
        if (.not. trial_app%particle_species(species_idx)%enabled) cycle
        if (.not. trial_app%particle_species(species_idx)%has_target_macro_particles_per_batch) cycle
        if (.not. any(trial_app%particle_species(species_idx)%boundary_inflow_low /= 0_i32) .and. &
            .not. any(trial_app%particle_species(species_idx)%boundary_inflow_high /= 0_i32)) cycle
        trial_app%particle_species(species_idx)%w_particle = &
          app%particle_species(species_idx)%w_particle*duration_ratio
      end do

      call perf_region_begin(perf_region_prepare_batch, t0)
      call build_particle_source_plan( &
        trial_app, source_plan, mpi=mpi_ctx, &
        kinetic_inflow_active=surface_closure%has_inflow_kinetic_map, &
        kinetic_reservoir_potential_v=surface_closure%inflow_reservoir_potential_v, &
        kinetic_access_potential_v=surface_closure%inflow_access_potential_v, &
        kinetic_inflow_face=surface_closure%inflow_kinetic_face, &
        number_flux_override_active=surface_closure%has_inflow_number_flux, &
        number_flux_override_m2_s=surface_closure%inflow_number_flux_m2_s, &
        kinetic_reservoir_density_m3=surface_closure%inflow_reservoir_density_m3 &
        )
      call prepare_batch_state( &
        mesh, trial_app, source_plan, snapshot, stats, batch_idx, workspace, pcls_batch, injection_layout, mpi_ctx, &
        inject_state, photo_failure &
        )
      call perf_region_end(perf_region_prepare_batch, t0)
      fresh_particle_count = pcls_batch%n
      call stop_on_photo_collision_failure(photo_failure)

      call perf_region_begin(perf_region_particle_batch, t0)
      ! Every trial of a batch starts the shared pool from slot 1 again.
      if (shared_pool) call pool_counter%reset(mpi_ctx)
      call process_particle_batch( &
        mesh, trial_app, boundary_contract, surface_closure, snapshot, pcls_batch, workspace, &
        bfield, batch_idx, mpi_ctx%rank, particle_team_size, track_failure, &
        batch_retry_counts, batch_substep_counts, &
        cell=cell, injection_layout=injection_layout, pool_counter=pool_counter, pool_claims=batch_pool_claims, &
        pool_claim_seconds=batch_pool_claim_seconds &
        )
      fresh_particle_count = pcls_batch%n
      call perf_region_end(perf_region_particle_batch, t0)

      if (replay_active) then
        particle_team_size_min = particle_team_size
        particle_team_size_max = particle_team_size
        call mpi_allreduce_min_i32_scalar(mpi_ctx, particle_team_size_min)
        call mpi_allreduce_max_i32_scalar(mpi_ctx, particle_team_size_max)
        if (particle_team_size_min /= nth .or. particle_team_size_max /= nth) then
          if (mpi_is_root(mpi_ctx)) then
            write (error_unit, '(a,i0,a,i0,a,i0)') &
              'replayed-trial OpenMP team size changed after replay probe: expected=', nth, &
              ' min=', particle_team_size_min, ' max=', particle_team_size_max
            flush (error_unit)
          end if
          error stop 'replayed-trial OpenMP team size changed after the replay probe.'
        end if
      end if

      call stop_on_particle_step_failure(track_failure)

      call apply_neutral_return_surface_closure(trial_app, pcls_batch, fresh_particle_count, workspace, mpi_ctx)
      call apply_fixed_surface_current_closure( &
        trial_app, surface_closure, pcls_batch, fresh_particle_count, workspace, mpi_ctx &
        )

      if (adaptive_nonzero_mode) then
        call prepare_adaptive_charge_candidate(mesh, app%sim, snapshot, conductor_solver, workspace, mpi_ctx)
        adaptive_metric_values = 0.0_dp
        if (mpi_is_root(mpi_ctx)) then
          call snapshot%measure_kneq0_potential_step( &
            mesh, workspace%candidate_charge, adaptive_metric_values(1) &
            )
        end if
        call mpi_allreduce_max_real_dp_array(mpi_ctx, adaptive_metric_values)
        adaptive_potential_step = adaptive_metric_values(1)
        trial_accepted = adaptive_potential_step <= app%periodic2%max_nonzero_mode_potential_step
      else
        trial_accepted = .true.
      end if
      if (trial_accepted) exit

      workspace%charge_candidate_ready = .false.
      if (trial_halvings >= adaptive_max_halvings) then
        error stop 'adaptive nonzero-mode batch failed after 24 duration halvings.'
      end if
      trial_halvings = trial_halvings + 1_i32
      trial_batch_duration = scale(app%sim%batch_duration, -trial_halvings)
      if (.not. ieee_is_finite(trial_batch_duration) .or. trial_batch_duration <= 0.0_dp) then
        error stop 'adaptive nonzero-mode batch duration became invalid.'
      end if
    end do
  end subroutine run_accepted_trial

  !> photo ray の衝突判定が失敗した rank があれば、最小 rank の失敗を報告して止める。
  subroutine stop_on_photo_collision_failure(failure)
    type(photo_ray_failure_type), intent(in) :: failure
    integer(i32) :: failure_count, failure_rank, local_values(4), selected_values(4)

    failure_count = merge(1_i32, 0_i32, failure%status /= collision_query_ok)
    call mpi_allreduce_sum_i32_scalar(mpi_ctx, failure_count)
    if (failure_count == 0_i32) return
    local_values = [failure%species, failure%ray, failure%bounce, failure%status]
    call mpi_select_lowest_rank_i32_values( &
      mpi_ctx, failure%status /= collision_query_ok, local_values, failure_rank, selected_values &
      )
    call stop_for_photo_collision_failure( &
      batch_idx, failure_rank, selected_values(1), selected_values(2), selected_values(3), selected_values(4) &
      )
  end subroutine stop_on_photo_collision_failure

  !> 粒子 step が失敗した rank があれば、最小 rank の失敗の粒子・位置・速度を報告して止める。
  subroutine stop_on_particle_step_failure(failure)
    type(particle_track_failure_type), intent(in) :: failure
    integer(i32) :: failure_count, failure_rank, local_values(3), selected_values(3)
    real(dp) :: selected_state(6)

    failure_count = merge(1_i32, 0_i32, failure%status /= collision_query_ok)
    call mpi_allreduce_sum_i32_scalar(mpi_ctx, failure_count)
    if (failure_count == 0_i32) return
    local_values = [failure%status, failure%particle, failure%step]
    call mpi_select_lowest_rank_i32_values( &
      mpi_ctx, failure%status /= collision_query_ok, local_values, failure_rank, selected_values &
      )
    selected_state = 0.0_dp
    if (mpi_ctx%rank == failure_rank) then
      selected_state(1:3) = failure%x
      selected_state(4:6) = failure%v
    end if
    call mpi_allreduce_sum_real_dp_array(mpi_ctx, selected_state)
    call stop_for_collision_failure( &
      batch_idx, failure_rank, selected_values(1), selected_values(2), &
      selected_values(3), trial_app%sim%dt, selected_state(1:3), selected_state(4:6) &
      )
  end subroutine stop_on_particle_step_failure

  !> 受理した trial の粒子の終端を数えて rank 合計し、soft-discard の上限を検査する。
  subroutine reduce_batch_outcomes()
    integer(i64) :: pool_claim_counts(1)
    real(dp) :: pool_claim_values(1)

    call perf_region_begin(perf_region_count_outcomes, t0)
    call count_batch_outcomes( &
      pcls_batch, workspace%escaped_boundary_flag, workspace%absorbed_flag, &
      workspace%soft_discarded_boundary_flag, batch_counts, batch_soft_discarded_abs_charge &
      )
    call perf_region_end(perf_region_count_outcomes, t0)
    call perf_region_begin(perf_region_mpi_reduce, t0)
    call mpi_allreduce_sum_i64_array(mpi_ctx, batch_counts)
    call mpi_allreduce_sum_i64_array(mpi_ctx, batch_retry_counts)
    call mpi_allreduce_sum_i64_array(mpi_ctx, batch_substep_counts)
    call mpi_allreduce_sum_real_dp_scalar(mpi_ctx, batch_soft_discarded_abs_charge)
    if (shared_pool) then
      pool_claim_values = [batch_pool_claim_seconds]
      pool_claim_counts = [batch_pool_claims]
      call mpi_allreduce_sum_real_dp_array(mpi_ctx, pool_claim_values)
      call mpi_allreduce_sum_i64_array(mpi_ctx, pool_claim_counts)
      batch_pool_claim_seconds = pool_claim_values(1)
      batch_pool_claims = pool_claim_counts(1)
    end if
    call perf_region_end(perf_region_mpi_reduce, t0)
    call enforce_soft_discard_limits( &
      app%sim, stats, batch_idx, batch_counts(6), batch_counts(1), batch_soft_discarded_abs_charge, mpi_ctx &
      )
    call report_soft_discard_summary(batch_idx, batch_counts(6), batch_soft_discarded_abs_charge, mpi_ctx)
  end subroutine reduce_batch_outcomes

  !> 受理した batch の統計更新を stats_candidate に組み立てる。overflow の検査は全て通ってから commit する。
  subroutine build_stats_candidate()
    integer(i64) :: secondary_count_values(3*app%n_particle_species)
    real(dp) :: secondary_real_values(2*app%n_particle_species), projected_simulated_time

    projected_simulated_time = stats%simulated_time + trial_batch_duration
    if (.not. ieee_is_finite(stats%simulated_time) .or. stats%simulated_time < 0.0_dp .or. &
        .not. ieee_is_finite(projected_simulated_time)) then
      error stop 'simulation statistic overflow: simulated_time'
    end if
    ! Build the complete statistics update before mutating mesh/ledger state.  The
    ! checked add routines are the single overflow guard; commit only after all pass.
    secondary_counts = sum(workspace%secondary_counts_thread, dim=3)
    secondary_values = sum(workspace%secondary_values_thread, dim=3)
    secondary_count_values = reshape(secondary_counts, [size(secondary_count_values)])
    secondary_real_values = reshape(secondary_values, [size(secondary_real_values)])
    call mpi_allreduce_sum_i64_array(mpi_ctx, secondary_count_values)
    call mpi_allreduce_sum_real_dp_array(mpi_ctx, secondary_real_values)
    secondary_counts = reshape(secondary_count_values, shape(secondary_counts))
    secondary_values = reshape(secondary_real_values, shape(secondary_values))
    stats_candidate = stats
    call accumulate_secondary_stats(stats_candidate, secondary_counts, secondary_values)
    call accumulate_batch_stats( &
      stats_candidate, batch_counts, batch_soft_discarded_abs_charge, batch_retry_counts, batch_substep_counts, 0.0_dp &
      )
    stats_candidate%simulated_time = projected_simulated_time
    if (shared_pool) then
      stats_candidate%particle_pool_claims = stats%particle_pool_claims + batch_pool_claims
      stats_candidate%particle_pool_claim_seconds = stats%particle_pool_claim_seconds + batch_pool_claim_seconds
    end if
    if (adaptive_nonzero_mode) then
      stats_candidate%adaptive_nonzero_mode%last_batch_duration = trial_batch_duration
      stats_candidate%adaptive_nonzero_mode%last_potential_step = adaptive_potential_step
      stats_candidate%adaptive_nonzero_mode%rejected_trials = checked_add_adaptive_rejected_trials( &
                                                              stats%adaptive_nonzero_mode%rejected_trials, &
                                                              trial_halvings &
                                                              )
    end if
  end subroutine build_stats_candidate

  !> 表面電荷の commit 前に、この batch の注入・終端・closure 補正を batch ledger に記録する。
  subroutine record_batch_ledger()
    call batch_ledger%reset(batch_idx)
    batch_ledger%surface_charge_before = finite_charge_sum(mesh%q_elem, 'batch surface charge before commit')
    call record_batch_initial_charge(trial_app, pcls_batch, fresh_particle_count, batch_ledger)
    call record_batch_outcome_charge( &
      pcls_batch, workspace%escaped_boundary_flag, workspace%absorbed_flag, &
      workspace%soft_discarded_boundary_flag, batch_ledger &
      )
    call reduce_charge_ledger_fluxes(batch_ledger, mpi_ctx, workspace)
    call copy_closure_corrections_to_ledger(workspace, batch_ledger)
  end subroutine record_batch_ledger

  !> 表面電荷を commit し、外部根の更新・壁電位の固定・matching-plane receipt を進めて統計を確定する。
  subroutine commit_accepted_batch()
    call perf_region_begin(perf_region_commit_charge, t0)
    call commit_batch_charge(mesh, app%sim, snapshot, conductor_solver, workspace, rel, mpi_ctx)
    call perf_region_end(perf_region_commit_charge, t0)
    call outflow_refresh%commit_batch( &
      trial_app, mesh, mpi_ctx, workspace%matching_plane_moments_thread, &
      workspace%fixed_current_charge_values(:, fixed_charge_emission), &
      workspace%fixed_escape_closed_charge, workspace%fixed_escape_unclosed_charge, &
      batch_idx, stats_candidate%matching_plane, surface_closure, outflow_refreshed &
      )
    if (outflow_refreshed) call apply_surface_plane_gauge(cell, mesh, surface_closure, snapshot, mpi_ctx, .false.)
    call record_matching_fields(stats_candidate, snapshot, surface_closure, mesh)
    if (fixed_history_enabled) call accumulate_fixed_current_window()
    stats_candidate%last_rel_change = rel
    if (ledger_enabled) then
      batch_ledger%surface_charge_after = finite_charge_sum(mesh%q_elem, 'batch surface charge after commit')
      call accumulate_charge_ledger(charge_ledger, batch_ledger)
    end if

    call perf_region_begin(perf_region_stats_update, t0)
    stats = stats_candidate
    call perf_region_end(perf_region_stats_update, t0)
  end subroutine commit_accepted_batch

  !> root rank で batch の進捗と各 history を書く。
  subroutine write_batch_history()
    if (secondary_history_enabled) then
      call write_secondary_emission_history_rows(secondary_hist_unit, batch_idx, app, secondary_counts, secondary_values)
    end if
    call print_batch_progress(batch_idx, final_batch_idx, rel)
    call maybe_write_history_snapshot(history_enabled, hist_unit, hist_stride, stats, rel, mesh%q_elem)
    if (matching_history_enabled .and. mod(batch_idx - 1_i32, hist_stride) == 0_i32) then
      call write_matching_plane_history_snapshot(matching_hist_unit, batch_idx, stats%simulated_time, stats)
    end if
    ! 外部源は更新時にだけ変わるので、スペクトルは間隔によらず更新した batch で書く。
    if (spectrum_history_enabled .and. outflow_refreshed) then
      call write_matching_plane_spectrum_history_rows(spectrum_hist_unit, batch_idx, stats%simulated_time, stats)
    end if
    if (fixed_history_enabled .and. mod(batch_idx - 1_i32, hist_stride) == 0_i32) then
      call write_fixed_current_history_rows( &
        fixed_hist_unit, batch_idx, stats%simulated_time, fixed_window_batches, fixed_window_active, &
        fixed_window_tracked, fixed_window_target &
        )
      call reset_fixed_current_window()
    end if
    if (potential_history_enabled) then
      call maybe_write_potential_history_snapshot( &
        potential_history_enabled, pot_hist_unit, hist_stride, stats, snapshot, mesh, app%sim, potential_buf, &
        top_reference_history_enabled, top_ref_hist_unit, &
        cell=cell &
        )
    end if
  end subroutine write_batch_history

  !> 終了時のメッシュ電位と静電診断（z-high 参照電位の統計を含む）を返す。
  subroutine finish_run_outputs()
    real(dp) :: top_phi_mean, top_phi_std, top_phi_min, top_phi_max

    if (present(outputs)) then
      if (outputs%mesh_potential .and. mpi_is_root(mpi_ctx)) then
        call snapshot%refresh(mesh)
        allocate (outputs%mesh_potential_v(mesh%nelem))
        call snapshot%compute_mesh_potential(mesh, app%sim, outputs%mesh_potential_v)
      end if
    end if
    if (present(electrostatic_diagnostics)) then
      call snapshot%get_diagnostics(electrostatic_diagnostics)
      if (mpi_is_root(mpi_ctx) .and. app%sim%use_box) then
        call snapshot%refresh(mesh)
        call compute_wall_box_potential_statistics( &
          mesh, app%sim, snapshot, top_phi_mean, top_phi_std, top_phi_min, top_phi_max, &
          cell=cell &
          )
        electrostatic_diagnostics%top_reference_available = .true.
        electrostatic_diagnostics%top_reference_last_batch = stats%batches
        electrostatic_diagnostics%top_reference_simulated_time = stats%simulated_time
        electrostatic_diagnostics%top_reference_z_high = cell%wall_coordinate()
        electrostatic_diagnostics%top_reference_sample_n = app%sim%injection_face_phi_grid_n
        electrostatic_diagnostics%top_reference_potential_mean = top_phi_mean
        electrostatic_diagnostics%top_reference_potential_std = top_phi_std
        electrostatic_diagnostics%top_reference_potential_min = top_phi_min
        electrostatic_diagnostics%top_reference_potential_max = top_phi_max
      end if
    end if
  end subroutine finish_run_outputs

  subroutine reset_fixed_current_window()
    fixed_window_batches = 0_i32
    fixed_window_tracked = 0.0_dp
    fixed_window_target = 0.0_dp
    fixed_window_active = .false.
  end subroutine reset_fixed_current_window

  !> accepted batch の固定電流 closure の追跡電荷と目標電荷を history 窓へ足す。
  subroutine accumulate_fixed_current_window()
    fixed_window_tracked(1, :) = fixed_window_tracked(1, :) + &
                                 workspace%fixed_current_charge_values(:, fixed_charge_absorbed)
    fixed_window_tracked(2, :) = fixed_window_tracked(2, :) + &
                                 workspace%fixed_current_charge_values(:, fixed_charge_emission)
    fixed_window_tracked(3, :) = fixed_window_tracked(3, :) + &
                                 workspace%fixed_current_charge_values(:, fixed_charge_escaped)
    fixed_window_target(1, :) = fixed_window_target(1, :) + workspace%fixed_absorbed_target_charge
    fixed_window_target(2, :) = fixed_window_target(2, :) + workspace%fixed_emission_target_charge
    fixed_window_target(3, :) = fixed_window_target(3, :) + workspace%fixed_escape_target_charge
    fixed_window_active = fixed_window_active .or. workspace%fixed_channel_active
    fixed_window_batches = fixed_window_batches + 1_i32
  end subroutine accumulate_fixed_current_window
  end procedure run_absorption_insulator

  !> neutral_return と固定電流 closure がこの batch に当てた補正を ledger の receipt へ写す。
  subroutine copy_closure_corrections_to_ledger(workspace, ledger)
    type(simulator_batch_workspace_type), intent(in) :: workspace
    type(charge_ledger_type), intent(inout) :: ledger

    ledger%neutral_return_correction = workspace%neutral_return_correction
    ledger%neutral_return_weight_scale = workspace%neutral_return_weight_scale
    ledger%neutral_return_unresolved_fraction = workspace%neutral_return_unresolved_fraction
    ledger%fixed_absorbed_target_charge = workspace%fixed_absorbed_target_charge
    ledger%fixed_absorbed_weight_scale = workspace%fixed_absorbed_weight_scale
    ledger%fixed_emission_target_charge = workspace%fixed_emission_target_charge
    ledger%fixed_emission_weight_scale = workspace%fixed_emission_weight_scale
    ledger%fixed_escape_target_charge = workspace%fixed_escape_target_charge
    ledger%fixed_escape_correction = workspace%fixed_escape_correction
    ledger%fixed_current_correction = workspace%fixed_current_correction
  end subroutine copy_closure_corrections_to_ledger

  !> 外部シース解が与える壁電位を、z-high面の水平平均電位として周期k=0成分へ固定する。
  subroutine apply_surface_plane_gauge(cell, mesh, surface_closure, snapshot, mpi, warn_if_unavailable)
    type(irradiation_cell_contract_type), intent(in) :: cell
    type(mesh_type), intent(in) :: mesh
    type(surface_closure_contract_type), intent(in) :: surface_closure
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    type(mpi_context), intent(in) :: mpi
    logical, intent(in) :: warn_if_unavailable

    if (.not. surface_closure%has_plane_gauge) return
    if (.not. snapshot%use_zero_mode) then
      if (warn_if_unavailable .and. mpi_is_root(mpi)) then
        write (error_unit, '(a)') &
          'WARNING: the surface-current outer barrier is compared with the z-high potential, but this field '// &
          'backend has no periodic zero mode to tie that potential to the outer wall potential phi0.'
        flush (error_unit)
      end if
      return
    end if
    call snapshot%set_matching_plane_top_field(mesh, surface_closure%top_field_v_m)
    call snapshot%set_matching_plane_gauge(mesh, cell%wall_coordinate(), surface_closure%plane_gauge_potential_v)
  end subroutine apply_surface_plane_gauge

  !> D_H is the BEACH-side displacement, including the selected lower-boundary field.
  subroutine record_matching_fields(stats, snapshot, closure, mesh)
    use bem_constants, only: eps0
    type(sim_stats), intent(inout) :: stats
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    type(mesh_type), intent(in) :: mesh
    type(surface_closure_contract_type), intent(in) :: closure

    stats%matching_plane%top_field_v_m = closure%top_field_v_m
    if (.not. snapshot%use_zero_mode) return
    call snapshot%refresh_zero_mode(mesh)
    stats%matching_plane%top_cell_mean_v_m = snapshot%zero_state%e_bottom + &
                                             snapshot%zero_state%total_charge/(eps0*snapshot%zero_plan%area_xy)
    if (stats%matching_plane%state_valid) then
      stats%matching_plane%displacement_c_m2 = eps0*stats%matching_plane%top_cell_mean_v_m
    end if
  end subroutine record_matching_fields

  function checked_add_adaptive_rejected_trials(accumulated, trial_halvings) result(total)
    integer(i64), intent(in) :: accumulated
    integer(i32), intent(in) :: trial_halvings
    integer(i64) :: increment, total

    increment = int(trial_halvings, i64)
    if (accumulated < 0_i64 .or. increment < 0_i64) then
      error stop 'adaptive rejected-trial count must be nonnegative.'
    end if
    if (accumulated > huge(accumulated) - increment) then
      error stop 'adaptive rejected-trial count overflow.'
    end if
    total = accumulated + increment
  end function checked_add_adaptive_rejected_trials

  subroutine stop_for_collision_failure(batch_idx, rank, status, particle, step, dt, x, v)
    integer(i32), intent(in) :: batch_idx, rank, status, particle, step
    real(dp), intent(in) :: dt, x(3), v(3)
    write (error_unit, '(a,i0,a,i0,a,i0,a,i0,a,a,a,i0,a,es13.5,a,3es13.5,a,3es13.5)') &
      'particle step failed: batch=', batch_idx, ' particle=', particle, ' step=', step, &
      ' rank=', rank, ' status=', trim(particle_status_name(status)), ' code=', status, &
      ' dt=', dt, ' x=', x, ' v=', v
    flush (error_unit)
    error stop 'particle step failed.'
  end subroutine stop_for_collision_failure

  subroutine stop_for_photo_collision_failure(batch_idx, rank, species, ray, bounce, status)
    integer(i32), intent(in) :: batch_idx, rank, species, ray, bounce, status
    write (error_unit, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,a,a,i0)') &
      'photo ray collision failed: batch=', batch_idx, ' rank=', rank, ' species=', species, &
      ' ray=', ray, ' bounce=', bounce, ' status=', trim(particle_status_name(status)), ' code=', status
    flush (error_unit)
    error stop 'photo ray collision query failed.'
  end subroutine stop_for_photo_collision_failure

  !> accepted trial の soft-discard 集約だけを root へ出力する。
  subroutine report_soft_discard_summary(batch_idx, global_count, global_abs_charge, mpi)
    integer(i32), intent(in) :: batch_idx
    integer(i64), intent(in) :: global_count
    real(dp), intent(in) :: global_abs_charge
    type(mpi_context), intent(in) :: mpi

    if (global_count == 0_i64) return

    if (mpi_is_root(mpi)) then
      write (error_unit, '(a,i0,a,i0,a,es13.5)') &
        'multiple_box_events soft discard accepted: batch=', batch_idx, &
        ' global_count=', global_count, ' global_abs_charge_C=', global_abs_charge
      flush (error_unit)
    end if
  end subroutine report_soft_discard_summary

  subroutine validate_soft_discard_initial_state(sim, stats)
    type(sim_config), intent(in) :: sim
    type(sim_stats), intent(in) :: stats

    if (trim(lower_ascii(sim%multiple_box_events_policy)) /= 'soft_discard') return
    if (stats%processed_particles < 0_i64 .or. stats%multiple_box_events_soft_discarded < 0_i64 .or. &
        stats%multiple_box_events_soft_discarded > stats%processed_particles .or. &
        .not. ieee_is_finite(stats%multiple_box_events_soft_discarded_abs_charge) .or. &
        stats%multiple_box_events_soft_discarded_abs_charge < 0.0_dp) then
      error stop 'soft-discard initial statistics are inconsistent.'
    end if
  end subroutine validate_soft_discard_initial_state

  subroutine enforce_soft_discard_limits( &
    sim, stats, batch_idx, batch_count, batch_processed, batch_abs_charge, mpi &
    )
    type(sim_config), intent(in) :: sim
    type(sim_stats), intent(in) :: stats
    integer(i32), intent(in) :: batch_idx
    integer(i64), intent(in) :: batch_count, batch_processed
    real(dp), intent(in) :: batch_abs_charge
    type(mpi_context), intent(in) :: mpi
    integer(i64) :: projected_count, projected_processed, count_grace
    real(dp) :: projected_fraction, projected_abs_charge
    logical :: fraction_exceeded, charge_warning

    if (trim(lower_ascii(sim%multiple_box_events_policy)) /= 'soft_discard') return
    count_grace = int(sim%multiple_box_events_soft_discard_count_grace, i64)
    if (batch_count < 0_i64 .or. batch_processed < 0_i64) then
      error stop 'soft-discard batch counters must be nonnegative.'
    end if
    if (stats%multiple_box_events_soft_discarded > huge(projected_count) - batch_count) then
      projected_count = huge(projected_count)
    else
      projected_count = stats%multiple_box_events_soft_discarded + batch_count
    end if
    if (stats%processed_particles > huge(projected_processed) - batch_processed) then
      projected_processed = huge(projected_processed)
    else
      projected_processed = stats%processed_particles + batch_processed
    end if
    projected_fraction = 0.0_dp
    if (projected_processed > 0_i64) then
      projected_fraction = real(projected_count, dp)/real(projected_processed, dp)
    end if
    projected_abs_charge = stats%multiple_box_events_soft_discarded_abs_charge + batch_abs_charge
    if (.not. ieee_is_finite(projected_abs_charge)) then
      error stop 'multiple_box_events soft-discard cumulative absolute charge is not finite.'
    end if
    fraction_exceeded = projected_count > count_grace .and. &
                        projected_fraction > sim%multiple_box_events_soft_discard_fraction_limit
    charge_warning = &
      stats%multiple_box_events_soft_discarded_abs_charge <= &
      sim%multiple_box_events_soft_discard_abs_charge_limit .and. &
      projected_abs_charge > sim%multiple_box_events_soft_discard_abs_charge_limit

    if (charge_warning .and. mpi_is_root(mpi)) then
      write (error_unit, '(a,i0,a,es13.5,a,es13.5)') &
        'WARNING: multiple_box_events soft-discard cumulative absolute-charge threshold crossed: batch=', &
        batch_idx, ' abs_charge_C=', projected_abs_charge, &
        ' particles.tracking.events.discard_charge_warning_c=', sim%multiple_box_events_soft_discard_abs_charge_limit
      flush (error_unit)
    end if
    if (.not. fraction_exceeded) return

    if (mpi_is_root(mpi)) then
      write (error_unit, '(a,i0,a,i0,a,i0,a,i0,a,es13.5,a,es13.5,a,es13.5,a,es13.5)') &
        'multiple_box_events soft-discard cumulative fraction limit exceeded: batch=', batch_idx, &
        ' count=', projected_count, ' particles.tracking.events.discard_count_grace=', count_grace, &
        ' processed=', projected_processed, ' fraction=', projected_fraction, &
        ' particles.tracking.events.discard_fraction_limit=', sim%multiple_box_events_soft_discard_fraction_limit, &
        ' abs_charge_C=', projected_abs_charge, &
        ' particles.tracking.events.discard_charge_warning_c=', sim%multiple_box_events_soft_discard_abs_charge_limit
      flush (error_unit)
    end if
    error stop 'multiple_box_events soft-discard cumulative fraction limit exceeded.'
  end subroutine enforce_soft_discard_limits

end submodule bem_simulator_loop

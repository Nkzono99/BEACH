!> `bem_simulator` の主ループと粒子処理計算を実装する submodule。
submodule(bem_simulator) bem_simulator_loop
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use bem_performance_profile, only: perf_region_batch_total, perf_region_begin, perf_region_commit_charge, &
                                     perf_region_count_outcomes, perf_region_end, perf_region_field_refresh, &
                                     perf_region_field_solver_init, perf_region_history_write, perf_region_mpi_reduce, &
                                     perf_region_particle_batch, perf_region_prepare_batch, perf_region_simulation_total, &
                                     perf_region_stats_update
  implicit none
contains

  !> 吸着モデルのバッチループを実行し、電荷更新と統計集計を進める。
  module procedure run_absorption_insulator
  integer(i32) :: batch_idx, final_batch_idx, batch_count_this_run, local_batch_idx, nth, hist_stride
  integer(i32) :: fresh_particle_count, due_escape_count_total
  integer(i64) :: outer_queue_event_count_before, outer_queue_event_count_after_pop, outer_queue_event_count_after
  integer(i32) :: collision_failure_count, collision_failure_rank, collision_failure_status
  integer(i32) :: collision_failure_particle, collision_failure_step
  integer(i32) :: local_failure_values(3), selected_failure_values(3)
  integer(i32) :: photo_failure_count, photo_failure_rank, photo_failure_status, photo_failure_species
  integer(i32) :: photo_failure_ray, photo_failure_bounce
  integer(i32) :: photo_local_failure_values(4), photo_selected_failure_values(4)
  integer :: hist_unit
  logical :: history_enabled
  logical :: potential_history_enabled
  integer :: pot_hist_unit
  real(dp), allocatable :: potential_buf(:)
  type(photoelectron_histogram_type), allocatable :: photoelectron_histogram_thread(:)
  type(photoelectron_histogram_type) :: photoelectron_batch_histogram
  integer(i32) :: batch_counts(6)
  real(dp) :: bfield(3), rel, t0, sim_t0, batch_t0, batch_soft_discarded_abs_charge
  real(dp) :: collision_failure_x(3), collision_failure_v(3), selected_failure_state(6)
  real(dp) :: max_outer_flight_time, max_frozen_field_ratio, max_outer_energy_relative_error
  real(dp) :: outer_max_diagnostics(3)
  real(dp) :: outer_queue_charge_before, outer_queue_charge_after_pop, outer_queue_charge_after
  real(dp) :: outer_queue_photoelectron_number_before, outer_queue_photoelectron_number_after_pop
  real(dp) :: outer_queue_photoelectron_number_after
  real(dp) :: outer_queue_area_xy, outer_queue_photoelectron_column_target
  character(len=16) :: outer_queue_fingerprint
  type(particles_soa) :: pcls_batch
  type(outer_event_record_type), allocatable :: due_outer_events(:)
  real(dp), allocatable :: due_returned_charge(:), due_escaped_charge(:)
  integer(i64), allocatable :: due_escaped_count(:)
  type(mpi_context) :: mpi_ctx
  type(electrostatic_snapshot_type) :: snapshot
  type(outer_coupler_type) :: outer_coupler
  type(charge_ledger_type) :: batch_ledger
  type(simulator_batch_workspace_type) :: workspace
  type(particle_source_plan_type) :: source_plan
  logical :: ledger_enabled
  logical :: outer_queue_enabled
  logical :: photoelectron_histogram_enabled
  integer(i32) :: thread_index, photoelectron_status
  integer(i32) :: kinetic_status
  integer(i32) :: boundary_status
  character(len=256) :: kinetic_message
  character(len=256) :: boundary_message
  type(kinetic_outer_plasma_options_type) :: kinetic_options
  type(external_boundary_contract_type) :: boundary_contract
  type(outer_plasma_state_type) :: steady_start_state
  real(dp) :: steady_start_charge

  stats = sim_stats()
  if (present(initial_stats)) stats = initial_stats
  mpi_ctx = mpi_context()
  if (present(mpi)) mpi_ctx = mpi
  call resolve_external_boundary_contract( &
    app%sim%reservoir_potential_model, app%sim%sheath_injection_model, app%sim%open_boundary_model, &
    app%outer_plasma%model, app%outer_plasma%kinetic_closure, app%outer_plasma%return_model, &
    app%coupling%particle_transfer_mode, app%coupling%outer_queue_enabled, boundary_contract, boundary_status, &
    boundary_message &
    )
  if (boundary_status /= external_boundary_ok) error stop trim(boundary_message)
  ledger_enabled = present(charge_ledger)
  outer_queue_enabled = app%coupling%outer_queue_enabled
  if (outer_queue_enabled .and. .not. present(outer_queue_state)) then
    error stop 'outer_queue_enabled=true requires a persistent outer queue state.'
  end if
  outer_queue_area_xy = 1.0_dp
  if (outer_queue_enabled) then
    outer_queue_area_xy = (app%sim%box_max(1) - app%sim%box_min(1))* &
                          (app%sim%box_max(2) - app%sim%box_min(2))
    if (.not. ieee_is_finite(outer_queue_area_xy) .or. outer_queue_area_xy <= 0.0_dp) then
      error stop 'outer_queue_enabled=true requires a positive finite horizontal box area.'
    end if
  end if
  allocate (due_returned_charge(app%n_particle_species), due_escaped_charge(app%n_particle_species))
  allocate (due_escaped_count(app%n_particle_species))
  if (ledger_enabled) then
    if (.not. allocated(charge_ledger%injected_from_remote)) then
      call charge_ledger%init(app%n_particle_species)
    else if (charge_ledger%nspecies /= app%n_particle_species) then
      error stop 'charge ledger species count does not match app config.'
    end if
    call batch_ledger%init(app%n_particle_species)
  end if

  nth = 1
!$ nth = max(1, omp_get_max_threads())
  call workspace%init(mesh%nelem, app%n_particle_species, nth)
  max_outer_flight_time = 0.0_dp
  max_frozen_field_ratio = 0.0_dp
  max_outer_energy_relative_error = 0.0_dp
  if (present(electrostatic_restart_state)) then
    if (electrostatic_restart_state%outer_max_diagnostics_complete) then
      outer_max_diagnostics = [ &
                              electrostatic_restart_state%max_outer_flight_time, &
                              electrostatic_restart_state%max_frozen_field_ratio, &
                              electrostatic_restart_state%max_outer_energy_relative_error &
                              ]
      if (.not. all(ieee_is_finite(outer_max_diagnostics)) .or. any(outer_max_diagnostics < 0.0_dp)) then
        error stop 'Restart cumulative outer diagnostics must be finite and nonnegative.'
      end if
      max_outer_flight_time = outer_max_diagnostics(1)
      max_frozen_field_ratio = outer_max_diagnostics(2)
      max_outer_energy_relative_error = outer_max_diagnostics(3)
    end if
  end if
  outer_queue_fingerprint = ''
  photoelectron_histogram_enabled = app%outer_plasma%photoelectron_histogram_enabled
  if (photoelectron_histogram_enabled .and. .not. present(photoelectron_state)) then
    error stop 'photoelectron_histogram_enabled=true requires a photoelectron histogram state.'
  end if
  if (photoelectron_histogram_enabled) then
    allocate (photoelectron_histogram_thread(nth))
    if (.not. photoelectron_state%ready) then
      call photoelectron_state%init( &
        app%outer_plasma%photoelectron_histogram_bins, app%outer_plasma%photoelectron_histogram_energy_max &
        )
    end if
    if (photoelectron_state%last_completed_batch /= stats%batches) then
      error stop 'Photoelectron histogram checkpoint batch does not match simulation stats.'
    end if
    do thread_index = 1_i32, nth
      call photoelectron_histogram_thread(thread_index)%init( &
        app%outer_plasma%photoelectron_histogram_bins, app%outer_plasma%photoelectron_histogram_energy_max &
        )
    end do
  end if

  history_enabled = present(history_unit)
  hist_unit = 0
  if (history_enabled) hist_unit = history_unit
  hist_stride = 1_i32
  if (present(history_stride)) hist_stride = max(1_i32, history_stride)
  potential_history_enabled = present(potential_history_unit)
  pot_hist_unit = 0
  if (potential_history_enabled) then
    pot_hist_unit = potential_history_unit
    allocate (potential_buf(mesh%nelem))
  end if
  bfield = app%sim%b0
  final_batch_idx = app%sim%batch_count
  if (stats%batches < 0_i32) error stop 'Initial simulation batch count must be >= 0.'
  if (stats%batches > final_batch_idx) then
    error stop 'sim.batch_count must be >= completed checkpoint batches when resuming.'
  end if
  batch_count_this_run = final_batch_idx - stats%batches
  call perf_region_begin(perf_region_simulation_total, sim_t0)
  call perf_region_begin(perf_region_field_solver_init, t0)
  if (trim(lower_ascii(app%outer_plasma%model)) == 'kinetic_1d') then
    call resolve_kinetic_outer_options(app, 0.0_dp, kinetic_options, kinetic_status, kinetic_message)
    if (kinetic_status /= outer_plasma_ok) error stop 'kinetic outer options: '//trim(kinetic_message)
    call snapshot%init( &
      mesh, app%sim, app%field, app%periodic2, app%panel, app%outer_plasma, &
      kinetic_options=kinetic_options, mpi=mpi_ctx &
      )
  else
    call snapshot%init(mesh, app%sim, app%field, app%periodic2, app%panel, app%outer_plasma, mpi=mpi_ctx)
  end if
  if (present(electrostatic_restart_state)) then
    call snapshot%restore_outer_state(electrostatic_restart_state)
  end if
  if (trim(lower_ascii(app%coupling%steady_start_mode)) == 'zhao_floating') then
    if (stats%batches == 0_i32) then
      if (snapshot%outer%ready) then
        error stop 'Zhao floating steady start cannot overwrite an outer state on a fresh run.'
      end if
      call initialize_zhao_floating_steady_start( &
        mesh, app%mesh_mode, app%sim, app%periodic2, app%coupling, kinetic_options, steady_start_state, &
        steady_start_charge, kinetic_status, kinetic_message &
        )
      if (kinetic_status /= outer_plasma_ok) then
        error stop 'Zhao floating steady start: '//trim(kinetic_message)
      end if
      snapshot%outer = steady_start_state
      if (mpi_is_root(mpi_ctx)) then
        write (output_unit, '(a,a,a,es24.16,a,es24.16,a,i0)') &
          'zhao_steady_start_branch=', steady_start_state%zhao_branch, &
          ' interface_field_V_m=', steady_start_state%interface_field, &
          ' surface_charge_C=', steady_start_charge, &
          ' mesh_id=', app%coupling%steady_start_mesh_id
      end if
    else
      if (.not. snapshot%outer%ready) then
        error stop 'Zhao floating steady-start resume requires a complete restored outer state.'
      end if
      if (mpi_is_root(mpi_ctx)) then
        write (output_unit, '(a,i0)') 'zhao_steady_start_restored_after_batches=', stats%batches
      end if
    end if
  end if
  if (present(electrostatic_restart_state)) then
    call outer_coupler%init(app%coupling, electrostatic_restart_state%last_outer_update_batch)
  else
    call outer_coupler%init(app%coupling)
  end if
  call perf_region_end(perf_region_field_solver_init, t0)

  do local_batch_idx = 1, batch_count_this_run
    call perf_region_begin(perf_region_batch_total, batch_t0)
    if (photoelectron_histogram_enabled) then
      do thread_index = 1_i32, nth
        call photoelectron_histogram_thread(thread_index)%reset()
      end do
    end if

    batch_idx = stats%batches + 1_i32
    due_returned_charge = 0.0_dp
    due_escaped_charge = 0.0_dp
    due_escaped_count = 0_i64
    due_escape_count_total = 0_i32
    outer_queue_charge_before = 0.0_dp
    outer_queue_charge_after_pop = 0.0_dp
    outer_queue_charge_after = 0.0_dp
    outer_queue_photoelectron_number_before = 0.0_dp
    outer_queue_photoelectron_number_after_pop = 0.0_dp
    outer_queue_photoelectron_number_after = 0.0_dp
    outer_queue_photoelectron_column_target = 0.0_dp
    outer_queue_event_count_before = 0_i64
    outer_queue_event_count_after_pop = 0_i64
    outer_queue_event_count_after = 0_i64
    if (outer_queue_enabled) then
      call measure_outer_queue_state( &
        app, outer_queue_state, mpi_ctx, outer_queue_charge_before, outer_queue_photoelectron_number_before, &
        outer_queue_event_count_before &
        )
      call outer_queue_state%pop_due( &
        real(batch_idx - 1_i32, dp)*app%sim%batch_duration, due_outer_events &
        )
      call tally_due_outer_events( &
        app, due_outer_events, due_returned_charge, due_escaped_charge, due_escaped_count, due_escape_count_total &
        )
      ! This post-pop global inventory is the hand-off point for the transient Zhao column closure.
      call measure_outer_queue_state( &
        app, outer_queue_state, mpi_ctx, outer_queue_charge_after_pop, outer_queue_photoelectron_number_after_pop, &
        outer_queue_event_count_after_pop &
        )
      outer_queue_photoelectron_column_target = outer_queue_photoelectron_number_after_pop/outer_queue_area_xy
      snapshot%kinetic_options%photoelectron_column_closure_enabled = .true.
      snapshot%kinetic_options%photoelectron_column_target_m2 = outer_queue_photoelectron_column_target
    end if
    call perf_region_begin(perf_region_field_refresh, t0)
    call outer_coupler%refresh( &
      snapshot, mesh, batch_idx, continuation_stage='pre_batch' &
      )
    call perf_region_end(perf_region_field_refresh, t0)

    call perf_region_begin(perf_region_prepare_batch, t0)
    if (local_batch_idx == 1_i32) call build_particle_source_plan(app, source_plan, mpi=mpi_ctx)
    call prepare_batch_state( &
      mesh, app, source_plan, snapshot, stats, batch_idx, workspace, pcls_batch, mpi_ctx, snapshot%outer, inject_state, &
      photo_failure_status, photo_failure_species, &
      photo_failure_ray, photo_failure_bounce &
      )
    call perf_region_end(perf_region_prepare_batch, t0)
    fresh_particle_count = pcls_batch%n

    photo_failure_count = merge(1_i32, 0_i32, photo_failure_status /= collision_query_ok)
    call mpi_allreduce_sum_i32_scalar(mpi_ctx, photo_failure_count)
    if (photo_failure_count > 0_i32) then
      photo_local_failure_values = [photo_failure_species, photo_failure_ray, photo_failure_bounce, photo_failure_status]
      call mpi_select_lowest_rank_i32_values( &
        mpi_ctx, photo_failure_status /= collision_query_ok, photo_local_failure_values, &
        photo_failure_rank, photo_selected_failure_values &
        )
      call stop_for_photo_collision_failure( &
        batch_idx, photo_failure_rank, photo_selected_failure_values(1), photo_selected_failure_values(2), &
        photo_selected_failure_values(3), photo_selected_failure_values(4) &
        )
    end if

    if (outer_queue_enabled) then
      call append_due_return_particles(due_outer_events, pcls_batch)
      call workspace%prepare_particle_flags(pcls_batch%n, outer_queue_enabled=.true.)
    end if

    if (ledger_enabled) then
      call batch_ledger%reset(batch_idx)
      batch_ledger%surface_charge_before = sum(mesh%q_elem)
      batch_ledger%outer_flight_charge_before = outer_queue_charge_before
      call record_batch_initial_charge(app, pcls_batch, fresh_particle_count, batch_ledger)
    end if

    call perf_region_begin(perf_region_particle_batch, t0)
    if (photoelectron_histogram_enabled) then
      call process_particle_batch( &
        mesh, app, boundary_contract, snapshot, pcls_batch, workspace%dq_thread, workspace%escaped_boundary_flag, &
        workspace%absorbed_flag, bfield, batch_idx, mpi_ctx%rank, &
        workspace%soft_discarded_boundary_flag, workspace%queued_outer_flag, workspace%outer_event_staging, &
        workspace%interface_outward_thread, workspace%interface_returned_thread, &
        collision_failure_status, collision_failure_particle, &
        collision_failure_step, collision_failure_x, collision_failure_v, workspace%interface_tau_max_thread, &
        workspace%interface_frozen_ratio_max_thread, workspace%interface_energy_error_max_thread, &
        photoelectron_histogram_thread &
        )
    else
      call process_particle_batch( &
        mesh, app, boundary_contract, snapshot, pcls_batch, workspace%dq_thread, workspace%escaped_boundary_flag, &
        workspace%absorbed_flag, bfield, batch_idx, mpi_ctx%rank, &
        workspace%soft_discarded_boundary_flag, workspace%queued_outer_flag, workspace%outer_event_staging, &
        workspace%interface_outward_thread, workspace%interface_returned_thread, &
        collision_failure_status, collision_failure_particle, &
        collision_failure_step, collision_failure_x, collision_failure_v, workspace%interface_tau_max_thread, &
        workspace%interface_frozen_ratio_max_thread, workspace%interface_energy_error_max_thread &
        )
    end if
    max_outer_flight_time = max(max_outer_flight_time, maxval(workspace%interface_tau_max_thread))
    max_frozen_field_ratio = max(max_frozen_field_ratio, maxval(workspace%interface_frozen_ratio_max_thread))
    max_outer_energy_relative_error = &
      max(max_outer_energy_relative_error, maxval(workspace%interface_energy_error_max_thread))
    call perf_region_end(perf_region_particle_batch, t0)

    collision_failure_count = merge(1_i32, 0_i32, collision_failure_status /= collision_query_ok)
    call mpi_allreduce_sum_i32_scalar(mpi_ctx, collision_failure_count)
    if (collision_failure_count > 0_i32) then
      local_failure_values = [collision_failure_status, collision_failure_particle, collision_failure_step]
      call mpi_select_lowest_rank_i32_values( &
        mpi_ctx, collision_failure_status /= collision_query_ok, local_failure_values, &
        collision_failure_rank, selected_failure_values &
        )
      selected_failure_state = 0.0_dp
      if (mpi_ctx%rank == collision_failure_rank) then
        selected_failure_state(1:3) = collision_failure_x
        selected_failure_state(4:6) = collision_failure_v
      end if
      call mpi_allreduce_sum_real_dp_array(mpi_ctx, selected_failure_state)
      call stop_for_collision_failure( &
        batch_idx, collision_failure_rank, selected_failure_values(1), selected_failure_values(2), &
        selected_failure_values(3), app%sim%dt, selected_failure_state(1:3), selected_failure_state(4:6) &
        )
    end if

    if (outer_queue_enabled) then
      call enqueue_staged_outer_events( &
        outer_queue_state, workspace%queued_outer_flag, workspace%outer_event_staging, pcls_batch%n &
        )
      call measure_outer_queue_state( &
        app, outer_queue_state, mpi_ctx, outer_queue_charge_after, outer_queue_photoelectron_number_after, &
        outer_queue_event_count_after &
        )
      workspace%interface_returned_thread(:, 1) = &
        workspace%interface_returned_thread(:, 1) + due_returned_charge
    end if

    if (ledger_enabled) then
      batch_ledger%interface_outward_gross = sum(workspace%interface_outward_thread, dim=2)
      batch_ledger%interface_returned_gross = sum(workspace%interface_returned_thread, dim=2)
      batch_ledger%escaped_to_infinity = batch_ledger%escaped_to_infinity + due_escaped_charge
      batch_ledger%escaped_count = batch_ledger%escaped_count + due_escaped_count
      batch_ledger%outer_flight_charge_after = outer_queue_charge_after
      call record_batch_outcome_charge( &
        pcls_batch, workspace%escaped_boundary_flag, workspace%absorbed_flag, &
        workspace%soft_discarded_boundary_flag, workspace%queued_outer_flag, batch_ledger &
        )
      call reduce_charge_ledger_fluxes(batch_ledger, mpi_ctx, workspace)
    end if
    if (photoelectron_histogram_enabled) then
      if (local_batch_idx == 1_i32) then
        call photoelectron_state%begin_batch(photoelectron_batch_histogram)
      else
        call photoelectron_batch_histogram%reset()
      end if
      do thread_index = 1_i32, nth
        call photoelectron_batch_histogram%merge(photoelectron_histogram_thread(thread_index))
      end do
      call reduce_photoelectron_histogram(photoelectron_batch_histogram, mpi_ctx)
      call validate_photoelectron_linear_applicability( &
        photoelectron_batch_histogram%total_signed_charge(), &
        app%outer_plasma%photoelectron_ambient_charge_scale, app%outer_plasma%max_photoelectron_charge_ratio, &
        photoelectron_status &
        )
      if (photoelectron_status /= photoelectron_applicability_ok) then
        error stop 'Photoelectron emission exceeds the configured linear-applicability limit.'
      end if
      call photoelectron_state%commit_batch(batch_idx, photoelectron_batch_histogram)
    end if

    call perf_region_begin(perf_region_commit_charge, t0)
    call commit_batch_charge( &
      mesh, app%sim%q_floor, app%sim%softening, app%sim%e0, app%sim%field_bc_mode, &
      workspace, rel, mpi_ctx &
      )
    call perf_region_end(perf_region_commit_charge, t0)

    if (outer_queue_enabled) then
      ! Complete each batch with a Zhao closure at the same time level as the
      ! committed surface charge and post-enqueue queue inventory.  Keeping this
      ! corrector in the per-batch operation sequence makes a straight run and a
      ! checkpoint/resume split use the same continuation seed at the next batch.
      outer_queue_photoelectron_column_target = outer_queue_photoelectron_number_after/outer_queue_area_xy
      snapshot%kinetic_options%photoelectron_column_closure_enabled = .true.
      snapshot%kinetic_options%photoelectron_column_target_m2 = outer_queue_photoelectron_column_target
      call perf_region_begin(perf_region_field_refresh, t0)
      call outer_coupler%refresh( &
        snapshot, mesh, batch_idx, continuation_stage='post_enqueue' &
        )
      call perf_region_end(perf_region_field_refresh, t0)
    end if

    if (ledger_enabled) then
      batch_ledger%surface_charge_after = sum(mesh%q_elem)
      call accumulate_charge_ledger(charge_ledger, batch_ledger)
    end if

    call perf_region_begin(perf_region_count_outcomes, t0)
    call count_batch_outcomes( &
      pcls_batch, workspace%escaped_boundary_flag, workspace%absorbed_flag, workspace%soft_discarded_boundary_flag, &
      workspace%queued_outer_flag, batch_counts, batch_soft_discarded_abs_charge &
      )
    batch_counts(3) = batch_counts(3) + due_escape_count_total
    batch_counts(4) = batch_counts(4) + due_escape_count_total
    call perf_region_end(perf_region_count_outcomes, t0)

    call perf_region_begin(perf_region_mpi_reduce, t0)
    call mpi_allreduce_sum_i32_array(mpi_ctx, batch_counts)
    call mpi_allreduce_sum_real_dp_scalar(mpi_ctx, batch_soft_discarded_abs_charge)
    call perf_region_end(perf_region_mpi_reduce, t0)

    call perf_region_begin(perf_region_stats_update, t0)
    call accumulate_batch_stats(stats, batch_counts, batch_soft_discarded_abs_charge, rel)
    call perf_region_end(perf_region_stats_update, t0)

    call perf_region_begin(perf_region_history_write, t0)
    if (mpi_is_root(mpi_ctx)) then
      call print_batch_progress(batch_idx, final_batch_idx, rel)
      if (outer_queue_enabled) then
        write (output_unit, '(a,i0,a,i0,a,i0,a,es13.5,a,es13.5,a,es13.5,a,a,a,es13.5)') &
          'BEACH outer-queue batch=', batch_idx, ' events_after_pop=', outer_queue_event_count_after_pop, &
          ' events_for_closure=', outer_queue_event_count_after, ' charge_after_C=', outer_queue_charge_after, &
          ' photoelectron_column_target_m-2=', outer_queue_photoelectron_column_target, &
          ' eta=', snapshot%outer%photoelectron_population_fraction, ' branch=', snapshot%outer%zhao_branch, &
          ' column_residual_m-2=', snapshot%outer%photoelectron_column_residual_per_area
      end if
      call maybe_write_history_snapshot(history_enabled, hist_unit, hist_stride, stats, rel, mesh%q_elem)
      if (potential_history_enabled) then
        call snapshot%refresh(mesh, update_outer=.false.)
        call maybe_write_potential_history_snapshot( &
          potential_history_enabled, pot_hist_unit, hist_stride, stats, &
          snapshot, mesh, app%sim, potential_buf &
          )
      end if
    end if
    call perf_region_end(perf_region_history_write, t0)

    if (stats%multiple_box_events_soft_discarded > &
        int(app%sim%multiple_box_events_soft_discard_count_limit, i64) .or. &
        stats%multiple_box_events_soft_discarded_abs_charge > &
        app%sim%multiple_box_events_soft_discard_abs_charge_limit) then
      if (mpi_is_root(mpi_ctx)) then
        write (error_unit, '(a,i0,a,i0,a,es13.5,a,es13.5)') &
          'multiple_box_events soft-discard limit exceeded: count=', &
          stats%multiple_box_events_soft_discarded, ' count_limit=', &
          app%sim%multiple_box_events_soft_discard_count_limit, ' abs_charge_C=', &
          stats%multiple_box_events_soft_discarded_abs_charge, ' abs_charge_limit_C=', &
          app%sim%multiple_box_events_soft_discard_abs_charge_limit
        flush (error_unit)
      end if
      error stop 1
    end if

    call perf_region_end(perf_region_batch_total, batch_t0)
  end do
  call perf_region_end(perf_region_simulation_total, sim_t0)

  if (present(mesh_potential_v)) then
    if (mpi_is_root(mpi_ctx)) then
      call perf_region_begin(perf_region_field_refresh, t0)
      call snapshot%refresh(mesh, update_outer=.false.)
      call perf_region_end(perf_region_field_refresh, t0)
      allocate (mesh_potential_v(mesh%nelem))
      call snapshot%compute_mesh_potential(mesh, app%sim, mesh_potential_v)
    end if
  end if
  outer_max_diagnostics = [max_outer_flight_time, max_frozen_field_ratio, max_outer_energy_relative_error]
  if (present(electrostatic_diagnostics) .or. present(electrostatic_restart_state)) then
    call mpi_allreduce_max_real_dp_array(mpi_ctx, outer_max_diagnostics)
    if (outer_queue_enabled) then
      call measure_outer_queue_state( &
        app, outer_queue_state, mpi_ctx, outer_queue_charge_after, outer_queue_photoelectron_number_after, &
        outer_queue_event_count_after &
        )
      call measure_outer_queue_fingerprint(outer_queue_state, mpi_ctx, outer_queue_fingerprint)
    end if
  end if
  if (present(electrostatic_diagnostics)) then
    if (.not. present(mesh_potential_v) .and. mpi_is_root(mpi_ctx)) then
      call snapshot%refresh(mesh, update_outer=.false.)
    end if
    call snapshot%get_diagnostics(electrostatic_diagnostics)
    electrostatic_diagnostics%last_outer_update_batch = outer_coupler%last_outer_update_batch
    electrostatic_diagnostics%max_outer_flight_time = outer_max_diagnostics(1)
    electrostatic_diagnostics%max_frozen_field_ratio = outer_max_diagnostics(2)
    electrostatic_diagnostics%max_outer_energy_relative_error = outer_max_diagnostics(3)
    if (outer_queue_enabled) then
      electrostatic_diagnostics%outer_queue_event_count = outer_queue_event_count_after
      electrostatic_diagnostics%outer_queue_signed_charge = outer_queue_charge_after
      electrostatic_diagnostics%outer_queue_fingerprint = outer_queue_fingerprint
    end if
  end if
  if (present(electrostatic_restart_state)) then
    call snapshot%export_restart_state(outer_coupler%last_outer_update_batch, electrostatic_restart_state)
    electrostatic_restart_state%max_outer_flight_time = outer_max_diagnostics(1)
    electrostatic_restart_state%max_frozen_field_ratio = outer_max_diagnostics(2)
    electrostatic_restart_state%max_outer_energy_relative_error = outer_max_diagnostics(3)
    electrostatic_restart_state%outer_max_diagnostics_complete = .true.
    if (outer_queue_enabled) then
      electrostatic_restart_state%outer_queue_event_count = outer_queue_event_count_after
      electrostatic_restart_state%outer_queue_signed_charge = outer_queue_charge_after
      electrostatic_restart_state%outer_queue_fingerprint = outer_queue_fingerprint
      electrostatic_restart_state%outer_queue_inventory_complete = .true.
    end if
  end if

  if (allocated(potential_buf)) deallocate (potential_buf)
  end procedure run_absorption_insulator

  !> 1バッチ分の粒子群と作業バッファを初期化する。
  module procedure prepare_batch_state
  batch_idx = stats%batches + 1_i32
  call workspace%reset_before_injection()
  if (present(inject_state)) then
    call init_particle_batch_from_config( &
      app, batch_idx, pcls_batch, inject_state, mesh=mesh, photo_emission_dq=workspace%photo_emission_dq, &
      snapshot=snapshot, outer_state=outer_state, mpi=mpi, collision_failure_status=collision_failure_status, &
      collision_failure_species=collision_failure_species, collision_failure_ray=collision_failure_ray, &
      collision_failure_bounce=collision_failure_bounce, source_plan=source_plan &
      )
  else
    call init_particle_batch_from_config( &
      app, batch_idx, pcls_batch, mesh=mesh, photo_emission_dq=workspace%photo_emission_dq, mpi=mpi, &
      snapshot=snapshot, outer_state=outer_state, collision_failure_status=collision_failure_status, &
      collision_failure_species=collision_failure_species, &
      collision_failure_ray=collision_failure_ray, collision_failure_bounce=collision_failure_bounce, &
      source_plan=source_plan &
      )
  end if
  if (collision_failure_status /= collision_query_ok) return
  call workspace%prepare_particle_flags(pcls_batch%n, outer_queue_enabled=app%coupling%outer_queue_enabled)
  end procedure prepare_batch_state

  !> 粒子を時間発展させ、衝突時の堆積電荷をスレッド別に集計する。
  module procedure process_particle_batch
  integer(i32) :: i, step, tid, nth, warn_stride, collision_status
  real(dp) :: x0(3), v0(3), x1(3), v1(3), qdep
  type(hit_info) :: hit
  type(particle_step_result) :: step_result
  type(particle_step_result) :: external_final_result
  type(external_step_trace_type) :: external_trace
  logical :: has_warn_stride, collision_failed, candidate_inside, used_event_resolver
  logical :: photoelectron_histogram_enabled
  integer(i32) :: species_idx, external_event

  nth = size(dq_thread, 2)
  call read_env_i32_local('BEACH_WARN_LONG_PARTICLE_STEPS', warn_stride, has_warn_stride)
  if (.not. has_warn_stride) warn_stride = 0_i32
  collision_failure_status = collision_query_ok
  collision_failure_particle = huge(0_i32)
  collision_failure_step = 0_i32
  collision_failure_x = 0.0_dp
  collision_failure_v = 0.0_dp
  interface_outward_thread = 0.0_dp
  interface_returned_thread = 0.0_dp
  interface_tau_max_thread = 0.0_dp
  interface_frozen_ratio_max_thread = 0.0_dp
  interface_energy_error_max_thread = 0.0_dp
  photoelectron_histogram_enabled = present(photoelectron_histogram_thread)

  !$omp parallel default(none) &
  !$omp shared(mesh,pcls_batch,app,boundary_contract,snapshot,dq_thread,bfield,escaped_boundary_flag,absorbed_flag,nth) &
  !$omp shared(soft_discarded_boundary_flag,queued_outer_flag,outer_event_staging) &
  !$omp shared(interface_outward_thread,interface_returned_thread) &
  !$omp shared(interface_tau_max_thread,interface_frozen_ratio_max_thread) &
  !$omp shared(interface_energy_error_max_thread) &
  !$omp shared(photoelectron_histogram_thread,photoelectron_histogram_enabled) &
  !$omp shared(warn_stride,batch_idx,mpi_rank,collision_failure_status,collision_failure_particle,collision_failure_step) &
  !$omp shared(collision_failure_x,collision_failure_v) &
  !$omp private(i,step,x0,v0,x1,v1,hit,step_result,external_final_result,external_trace,tid,qdep,species_idx) &
  !$omp private(external_event) &
  !$omp private(collision_status,collision_failed,candidate_inside,used_event_resolver)
  ! スレッドごとに dq_thread(:, tid) を使って原子的更新なしで電荷を集める。
  tid = 1
!$ tid = omp_get_thread_num() + 1
  !$omp do schedule(dynamic, 1)
  do i = 1, pcls_batch%n
    if (.not. pcls_batch%alive(i)) cycle
    collision_failed = .false.
    do step = 1, app%sim%max_step
      x0 = pcls_batch%x(:, i)
      v0 = pcls_batch%v(:, i)
      call build_particle_step_candidate( &
        mesh, app%sim, snapshot, bfield, x0, v0, &
        pcls_batch%q(i), pcls_batch%m(i), app%sim%dt, x1, v1 &
        )
      if (.not. all(ieee_is_finite(x1)) .or. .not. all(ieee_is_finite(v1))) then
        !$omp critical (beach_collision_query_failure)
        if (collision_failure_status == collision_query_ok .or. i < collision_failure_particle .or. &
            (i == collision_failure_particle .and. step < collision_failure_step)) then
          collision_failure_status = particle_step_invalid_boundary
          collision_failure_particle = i
          collision_failure_step = step
          collision_failure_x = x0
          collision_failure_v = v0
        end if
        !$omp end critical (beach_collision_query_failure)
        collision_failed = .true.
        exit
      end if
      call find_first_hit(mesh, x0, x1, hit, sim=app%sim, status=collision_status)
      candidate_inside = .not. app%sim%use_box .or. &
                         (all(x1 > app%sim%box_min) .and. all(x1 < app%sim%box_max))
      used_event_resolver = .false.
      if (collision_status /= collision_query_ok) then
        if (.not. candidate_inside) then
          call resolve_particle_boundary_candidate( &
            mesh, app%sim, snapshot, bfield, x0, v0, pcls_batch%q(i), pcls_batch%m(i), app%sim%dt, x1, v1, &
            result=step_result, boundary_contract=boundary_contract &
            )
          used_event_resolver = .true.
        else
          !$omp critical (beach_collision_query_failure)
          if (collision_failure_status == collision_query_ok .or. i < collision_failure_particle .or. &
              (i == collision_failure_particle .and. step < collision_failure_step)) then
            collision_failure_status = collision_status
            collision_failure_particle = i
            collision_failure_step = step
            collision_failure_x = x0
            collision_failure_v = v0
          end if
          !$omp end critical (beach_collision_query_failure)
          collision_failed = .true.
          exit
        end if
      else if (candidate_inside) then
        if (hit%has_hit) then
          qdep = pcls_batch%q(i)*pcls_batch%w(i)
          dq_thread(hit%elem_idx, tid) = dq_thread(hit%elem_idx, tid) + qdep
          pcls_batch%alive(i) = .false.
          absorbed_flag(i) = .true.
          exit
        end if
        pcls_batch%x(:, i) = x1
        pcls_batch%v(:, i) = v1
      else
        call resolve_particle_boundary_candidate( &
          mesh, app%sim, snapshot, bfield, x0, v0, pcls_batch%q(i), pcls_batch%m(i), app%sim%dt, x1, v1, &
          hit=hit, result=step_result, boundary_contract=boundary_contract &
          )
        used_event_resolver = .true.
      end if
      if (used_event_resolver) then
        if (step_result%interface_crossing%has_crossing) then
          call continue_external_particle_step( &
            boundary_contract, snapshot, mesh, app%sim, app%outer_plasma, app%coupling, bfield, &
            pcls_batch%q(i), pcls_batch%m(i), app%sim%batch_duration, step_result, external_final_result, &
            external_trace &
            )
          step_result = external_final_result
          species_idx = pcls_batch%species_id(i)
          qdep = pcls_batch%q(i)*pcls_batch%w(i)
          do external_event = 1_i32, external_trace%count
            if (photoelectron_histogram_enabled .and. &
                trim(lower_ascii(app%particle_species(species_idx)%source_mode)) == 'photo_raycast') then
              call photoelectron_histogram_thread(tid)%add( &
                pcls_batch%q(i), pcls_batch%m(i), pcls_batch%w(i), &
                external_trace%crossing(external_event)%velocity &
                )
            end if
            interface_outward_thread(species_idx, tid) = interface_outward_thread(species_idx, tid) + qdep
            if (external_trace%outcome(external_event)%kind == interface_outcome_returned_local .or. &
                external_trace%outcome(external_event)%kind == interface_outcome_escaped_to_infinity .or. &
                external_trace%outcome(external_event)%kind == interface_outcome_queued_outer) then
              interface_energy_error_max_thread(tid) = max( &
                                                       interface_energy_error_max_thread(tid), &
                                                       external_trace%outcome(external_event)%energy_relative_error &
                                                       )
              interface_tau_max_thread(tid) = max( &
                                              interface_tau_max_thread(tid), &
                                              external_trace%outcome(external_event)%outer_flight_time &
                                              )
              interface_frozen_ratio_max_thread(tid) = max( &
                                                       interface_frozen_ratio_max_thread(tid), &
                                                       external_trace%outcome(external_event)%frozen_field_ratio &
                                                       )
            end if
            if (.not. boundary_contract%queue_enabled .and. &
                external_trace%outcome(external_event)%kind == interface_outcome_returned_local) then
              interface_returned_thread(species_idx, tid) = interface_returned_thread(species_idx, tid) + qdep
            end if
          end do
          if (boundary_contract%queue_enabled .and. external_trace%count > 0_i32 .and. &
              external_trace%outcome(external_trace%count)%kind == interface_outcome_queued_outer) then
            call stage_queued_outer_event( &
              outer_event_staging(i), external_trace%outcome(external_trace%count), pcls_batch, i, batch_idx, &
              mpi_rank, app%sim%batch_duration &
              )
            pcls_batch%alive(i) = .false.
            queued_outer_flag(i) = .true.
            exit
          end if
        end if
        if (step_result%status /= collision_query_ok) then
          if (step_result%status == particle_step_multiple_box_events .and. &
              trim(lower_ascii(app%sim%multiple_box_events_policy)) == 'soft_discard') then
            qdep = pcls_batch%q(i)*pcls_batch%w(i)
            pcls_batch%alive(i) = .false.
            soft_discarded_boundary_flag(i) = .true.
            !$omp critical (beach_soft_discard_boundary_event)
            write (output_unit, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,es13.5,a,3es13.5,a,3es13.5)') &
              'BEACH soft-discard: batch=', batch_idx, ' rank=', mpi_rank, ' particle=', i, &
              ' species=', pcls_batch%species_id(i), ' step=', step, ' macro_charge_C=', qdep, &
              ' x=', x0, ' v=', v0
            flush (output_unit)
            !$omp end critical (beach_soft_discard_boundary_event)
            exit
          end if
          !$omp critical (beach_collision_query_failure)
          if (collision_failure_status == collision_query_ok .or. i < collision_failure_particle .or. &
              (i == collision_failure_particle .and. step < collision_failure_step)) then
            collision_failure_status = step_result%status
            collision_failure_particle = i
            collision_failure_step = step
            collision_failure_x = x0
            collision_failure_v = v0
          end if
          !$omp end critical (beach_collision_query_failure)
          collision_failed = .true.
          exit
        end if
        if (step_result%absorbed) then
          qdep = pcls_batch%q(i)*pcls_batch%w(i)
          dq_thread(step_result%elem_idx, tid) = dq_thread(step_result%elem_idx, tid) + qdep
          pcls_batch%alive(i) = .false.
          absorbed_flag(i) = .true.
          exit
        end if
        if (step_result%escaped_boundary) then
          pcls_batch%alive(i) = .false.
          escaped_boundary_flag(i) = .true.
          exit
        end if

        pcls_batch%x(:, i) = step_result%x
        pcls_batch%v(:, i) = step_result%v
      end if
      if (warn_stride > 0_i32) then
        if (modulo(step, warn_stride) == 0_i32) then
          !$omp critical (beach_long_particle_warn)
          write (output_unit, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,3es13.5,a,3es13.5)') &
            'BEACH long-particle batch=', batch_idx, ' rank=', mpi_rank, ' thread=', tid, &
            ' particle=', i, ' species=', pcls_batch%species_id(i), ' step=', step, &
            ' x=', pcls_batch%x(:, i), ' v=', pcls_batch%v(:, i)
          flush (output_unit)
          !$omp end critical (beach_long_particle_warn)
        end if
      end if
    end do
    if (warn_stride > 0_i32 .and. pcls_batch%alive(i) .and. .not. collision_failed) then
      !$omp critical (beach_long_particle_warn)
      write (output_unit, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,3es13.5,a,3es13.5)') &
        'BEACH max-step-survivor batch=', batch_idx, ' rank=', mpi_rank, ' thread=', tid, &
        ' particle=', i, ' species=', pcls_batch%species_id(i), ' step=', app%sim%max_step, &
        ' x=', pcls_batch%x(:, i), ' v=', pcls_batch%v(:, i)
      flush (output_unit)
      !$omp end critical (beach_long_particle_warn)
    end if
  end do
  !$omp end do
  !$omp end parallel
  end procedure process_particle_batch

  !> OpenMP 粒子ループでは index 固有 record だけを書き、queue 本体は後段で更新する。
  subroutine stage_queued_outer_event(event, outcome, pcls_batch, particle_index, batch_idx, mpi_rank, batch_duration)
    type(outer_event_record_type), intent(out) :: event
    type(interface_particle_outcome_type), intent(in) :: outcome
    type(particles_soa), intent(in) :: pcls_batch
    integer(i32), intent(in) :: particle_index, batch_idx, mpi_rank
    real(dp), intent(in) :: batch_duration

    event = outer_event_record_type()
    ! Crossing times are not synchronized inside a BEACH batch; use the batch midpoint as the closure time.
    event%queued_time = (real(batch_idx, dp) - 0.5_dp)*batch_duration
    event%due_time = event%queued_time + outcome%outer_flight_time
    if (outcome%kind /= interface_outcome_queued_outer) then
      error stop 'outer-event staging requires a queued interface outcome.'
    end if
    select case (outcome%queued_terminal_kind)
    case (interface_outcome_returned_local)
      event%outcome = outer_event_outcome_return
    case (interface_outcome_escaped_to_infinity)
      event%outcome = outer_event_outcome_escape
    case default
      error stop 'queued outer-interface outcome has no valid terminal kind.'
    end select
    event%species_id = pcls_batch%species_id(particle_index)
    event%origin_rank = mpi_rank
    event%origin_batch = batch_idx
    event%origin_particle_id = int(particle_index, i64)
    event%interface_face_index = 6_i32
    event%q = pcls_batch%q(particle_index)
    event%m = pcls_batch%m(particle_index)
    event%w = pcls_batch%w(particle_index)
    event%x = outcome%position
    event%v = outcome%velocity
  end subroutine stage_queued_outer_event

  !> 全 rank の process failure が無いことを確認した後、粒子 index 順で deterministic に enqueue する。
  subroutine enqueue_staged_outer_events(queue, queued_outer_flag, staging, particle_count)
    type(outer_event_queue_type), intent(inout) :: queue
    logical, intent(in) :: queued_outer_flag(:)
    type(outer_event_record_type), intent(in) :: staging(:)
    integer(i32), intent(in) :: particle_count
    integer(i32) :: particle_index

    do particle_index = 1_i32, particle_count
      if (queued_outer_flag(particle_index)) call queue%enqueue(staging(particle_index))
    end do
  end subroutine enqueue_staged_outer_events

  !> rank-local queue stock を global signed charge、photoelectron number、event count に集約する。
  subroutine measure_outer_queue_state(app, queue, mpi, signed_charge, photoelectron_number, event_count)
    type(app_config), intent(in) :: app
    type(outer_event_queue_type), intent(in) :: queue
    type(mpi_context), intent(in) :: mpi
    real(dp), intent(out) :: signed_charge, photoelectron_number
    integer(i64), intent(out) :: event_count
    integer(i32) :: species_idx
    integer(i64) :: count_buffer(1)

    signed_charge = queue%signed_charge()
    photoelectron_number = 0.0_dp
    do species_idx = 1_i32, app%n_particle_species
      if (app%particle_species(species_idx)%enabled .and. &
          trim(lower_ascii(app%particle_species(species_idx)%source_mode)) == 'photo_raycast' .and. &
          app%particle_species(species_idx)%q_particle < 0.0_dp) then
        photoelectron_number = photoelectron_number + queue%species_particle_number(species_idx)
      end if
    end do
    count_buffer(1) = int(queue%size(), i64)
    call mpi_allreduce_sum_real_dp_scalar(mpi, signed_charge)
    call mpi_allreduce_sum_real_dp_scalar(mpi, photoelectron_number)
    call mpi_allreduce_sum_i64_array(mpi, count_buffer)
    event_count = count_buffer(1)
  end subroutine measure_outer_queue_state

  !> Rank-local queue identitiesをcollectiveに畳み込み、summary用のglobal fingerprintを得る。
  subroutine measure_outer_queue_fingerprint(queue, mpi, fingerprint)
    type(outer_event_queue_type), intent(in) :: queue
    type(mpi_context), intent(in) :: mpi
    character(len=16), intent(out) :: fingerprint
    integer(i64) :: components(2)

    call queue%fingerprint_components(mpi%rank, components)
    call mpi_allreduce_sum_i64_array(mpi, components)
    fingerprint = outer_event_queue_global_fingerprint(components, mpi%size)
  end subroutine measure_outer_queue_fingerprint

  !> batch start に due になった outer event を return stock と infinity escape flux に分ける。
  subroutine tally_due_outer_events( &
    app, events, returned_charge, escaped_charge, escaped_count, escaped_count_total &
    )
    type(app_config), intent(in) :: app
    type(outer_event_record_type), intent(in) :: events(:)
    real(dp), intent(out) :: returned_charge(:), escaped_charge(:)
    integer(i64), intent(out) :: escaped_count(:)
    integer(i32), intent(out) :: escaped_count_total
    integer(i32) :: event_index, species_idx
    real(dp) :: macro_charge

    returned_charge = 0.0_dp
    escaped_charge = 0.0_dp
    escaped_count = 0_i64
    escaped_count_total = 0_i32
    if (size(returned_charge) /= app%n_particle_species .or. &
        size(escaped_charge) /= app%n_particle_species .or. size(escaped_count) /= app%n_particle_species) then
      error stop 'outer-event due tally species size mismatch.'
    end if
    do event_index = 1_i32, int(size(events), i32)
      species_idx = events(event_index)%species_id
      if (species_idx < 1_i32 .or. species_idx > app%n_particle_species) then
        error stop 'outer-event due tally has an invalid species index.'
      end if
      macro_charge = events(event_index)%q*events(event_index)%w
      select case (events(event_index)%outcome)
      case (outer_event_outcome_return)
        returned_charge(species_idx) = returned_charge(species_idx) + macro_charge
      case (outer_event_outcome_escape)
        escaped_charge(species_idx) = escaped_charge(species_idx) + macro_charge
        escaped_count(species_idx) = escaped_count(species_idx) + 1_i64
        escaped_count_total = escaped_count_total + 1_i32
      case default
        error stop 'outer-event due tally has an invalid outcome.'
      end select
    end do
  end subroutine tally_due_outer_events

  !> due return event を fresh source の後ろへ追加し、local-domain particle として再開する。
  subroutine append_due_return_particles(events, pcls_batch)
    type(outer_event_record_type), intent(in) :: events(:)
    type(particles_soa), intent(inout) :: pcls_batch
    real(dp), allocatable :: x(:, :), v(:, :), q(:), m(:), w(:)
    integer(i32), allocatable :: species_id(:)
    integer(i32) :: event_index, return_index, return_count

    return_count = 0_i32
    do event_index = 1_i32, int(size(events), i32)
      if (events(event_index)%outcome == outer_event_outcome_return) return_count = return_count + 1_i32
    end do
    if (return_count == 0_i32) return

    allocate (x(3, return_count), v(3, return_count), q(return_count), m(return_count), w(return_count))
    allocate (species_id(return_count))
    return_index = 0_i32
    do event_index = 1_i32, int(size(events), i32)
      if (events(event_index)%outcome /= outer_event_outcome_return) cycle
      return_index = return_index + 1_i32
      x(:, return_index) = events(event_index)%x
      v(:, return_index) = events(event_index)%v
      q(return_index) = events(event_index)%q
      m(return_index) = events(event_index)%m
      w(return_index) = events(event_index)%w
      species_id(return_index) = events(event_index)%species_id
    end do
    call append_particles(pcls_batch, x, v, q, m, w, species_id)
  end subroutine append_due_return_particles

  subroutine reduce_photoelectron_histogram(histogram, mpi)
    type(photoelectron_histogram_type), intent(inout) :: histogram
    type(mpi_context), intent(in) :: mpi

    call mpi_allreduce_sum_real_dp_array(mpi, histogram%signed_charge)
    call mpi_allreduce_sum_real_dp_array(mpi, histogram%kinetic_energy)
    call mpi_allreduce_sum_real_dp_array(mpi, histogram%tangential_momentum_x)
    call mpi_allreduce_sum_real_dp_array(mpi, histogram%tangential_momentum_y)
    call mpi_allreduce_sum_i64_array(mpi, histogram%count)
  end subroutine reduce_photoelectron_histogram

  !> 全 rank で選択済みの collision failure context を報告して停止する。
  subroutine stop_for_collision_failure( &
    batch_idx, failure_rank, failure_status, failure_particle, failure_step, dt, failure_x, failure_v &
    )
    integer(i32), intent(in) :: batch_idx, failure_rank, failure_status, failure_particle, failure_step
    real(dp), intent(in) :: dt, failure_x(3), failure_v(3)
    character(len=32) :: failure_name
    character(len=512) :: failure_message

    select case (failure_status)
    case (collision_query_image_limit)
      failure_name = 'image_limit'
    case (collision_query_index_range)
      failure_name = 'index_range'
    case (collision_query_invalid_segment)
      failure_name = 'invalid_segment'
    case (collision_query_grid_stalled)
      failure_name = 'grid_stalled'
    case (particle_step_invalid_boundary)
      failure_name = 'invalid_boundary'
    case (particle_step_multiple_box_events)
      failure_name = 'multiple_box_events'
    case (particle_step_ambiguous_open_corner)
      failure_name = 'ambiguous_open_corner'
    case (particle_step_multiple_external_events)
      failure_name = 'multiple_external_events'
    case default
      failure_name = 'unknown'
    end select
    write (failure_message, '(a,i0,a,i0,a,i0,a,i0,a,a,a,i0,a,es13.5,a,3es13.5,a,3es13.5)') &
      'particle step failed: batch=', batch_idx, ' particle=', failure_particle, &
      ' step=', failure_step, ' rank=', failure_rank, ' status=', trim(failure_name), &
      ' code=', failure_status, ' dt=', dt, ' x=', failure_x, ' v=', failure_v
    write (error_unit, '(a)') trim(failure_message)
    flush (error_unit)
    error stop 1
  end subroutine stop_for_collision_failure

  !> 全 rank で選択済みの photo collision failure context を報告して停止する。
  subroutine stop_for_photo_collision_failure( &
    batch_idx, failure_rank, failure_species, failure_ray, failure_bounce, failure_status &
    )
    integer(i32), intent(in) :: batch_idx, failure_rank, failure_species, failure_ray, failure_bounce, failure_status
    character(len=16) :: failure_name
    character(len=256) :: failure_message

    select case (failure_status)
    case (collision_query_image_limit)
      failure_name = 'image_limit'
    case (collision_query_index_range)
      failure_name = 'index_range'
    case (collision_query_invalid_segment)
      failure_name = 'invalid_segment'
    case (collision_query_grid_stalled)
      failure_name = 'grid_stalled'
    case default
      failure_name = 'unknown'
    end select
    write (failure_message, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,a,a,i0)') &
      'photo collision query incomplete: batch=', batch_idx, ' rank=', failure_rank, &
      ' species=', failure_species, ' ray=', failure_ray, ' bounce=', failure_bounce, &
      ' status=', trim(failure_name), ' code=', failure_status
    write (error_unit, '(a)') trim(failure_message)
    flush (error_unit)
    error stop 1
  end subroutine stop_for_photo_collision_failure

  !> 正の整数環境変数を読む。未設定、不正値、ゼロ以下の場合は found=.false.。
  subroutine read_env_i32_local(name, value, found)
    character(len=*), intent(in) :: name
    integer(i32), intent(out) :: value
    logical, intent(out) :: found
    character(len=64) :: raw
    integer :: length, status, ios

    value = 0_i32
    found = .false.
    raw = ''
    call get_environment_variable(name, raw, length=length, status=status)
    if (status /= 0 .or. length <= 0 .or. length > len(raw)) return

    read (raw(:length), *, iostat=ios) value
    if (ios /= 0 .or. value <= 0_i32) then
      value = 0_i32
      return
    end if
    found = .true.
  end subroutine read_env_i32_local

  !> スレッド別電荷差分を合算してメッシュへ反映し、相対変化量を返す。
  module procedure commit_batch_charge
  real(dp) :: norm_dq, norm_q

  workspace%q_before = mesh%q_elem
  workspace%dq = sum(workspace%dq_thread, dim=2) + workspace%photo_emission_dq
  call mpi_allreduce_sum_real_dp_array(mpi, workspace%dq)
  mesh%q_elem = mesh%q_elem + workspace%dq
  call apply_surface_model_charge_relaxation(mesh, softening, external_e, field_bc_mode=field_bc_mode)
  workspace%dq = mesh%q_elem - workspace%q_before
  norm_dq = sqrt(sum(workspace%dq*workspace%dq))
  norm_q = sqrt(sum(mesh%q_elem*mesh%q_elem))
  rel = norm_dq/max(norm_q, q_floor)
  end procedure commit_batch_charge

  !> batch 初期粒子を remote injection と surface emission に分類する。
  subroutine record_batch_initial_charge(app, pcls_batch, fresh_particle_count, ledger)
    type(app_config), intent(in) :: app
    type(particles_soa), intent(in) :: pcls_batch
    integer(i32), intent(in) :: fresh_particle_count
    type(charge_ledger_type), intent(inout) :: ledger
    integer(i32) :: i, species_idx
    real(dp) :: macro_charge
    character(len=32) :: source_mode

    if (fresh_particle_count < 0_i32 .or. fresh_particle_count > pcls_batch%n) then
      error stop 'fresh particle count is invalid while recording charge ledger input.'
    end if
    do i = 1, fresh_particle_count
      species_idx = pcls_batch%species_id(i)
      if (species_idx < 1_i32 .or. species_idx > ledger%nspecies) then
        error stop 'particle species index is invalid while recording charge ledger input.'
      end if
      macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
      source_mode = lower_ascii(trim(app%particle_species(species_idx)%source_mode))
      if (trim(source_mode) == 'photo_raycast') then
        ledger%emitted_from_surface(species_idx) = ledger%emitted_from_surface(species_idx) + macro_charge
        ledger%emitted_count(species_idx) = ledger%emitted_count(species_idx) + 1_i64
      else
        ledger%injected_from_remote(species_idx) = ledger%injected_from_remote(species_idx) + macro_charge
        ledger%injected_count(species_idx) = ledger%injected_count(species_idx) + 1_i64
      end if
    end do
  end subroutine record_batch_initial_charge

  !> batch 終了粒子を surface absorption、infinity escape、unresolved discard に分類する。
  subroutine record_batch_outcome_charge( &
    pcls_batch, escaped_boundary_flag, absorbed_flag, soft_discarded_boundary_flag, queued_outer_flag, ledger &
    )
    type(particles_soa), intent(in) :: pcls_batch
    logical, intent(in) :: escaped_boundary_flag(:), absorbed_flag(:)
    logical, intent(in) :: soft_discarded_boundary_flag(:)
    logical, intent(in) :: queued_outer_flag(:)
    type(charge_ledger_type), intent(inout) :: ledger
    integer(i32) :: i, species_idx
    real(dp) :: macro_charge

    do i = 1, pcls_batch%n
      species_idx = pcls_batch%species_id(i)
      if (species_idx < 1_i32 .or. species_idx > ledger%nspecies) then
        error stop 'particle species index is invalid while recording charge ledger output.'
      end if
      macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
      if (queued_outer_flag(i)) then
        cycle
      else if (absorbed_flag(i)) then
        ledger%absorbed_on_surface(species_idx) = ledger%absorbed_on_surface(species_idx) + macro_charge
        ledger%absorbed_count(species_idx) = ledger%absorbed_count(species_idx) + 1_i64
      else if (escaped_boundary_flag(i)) then
        ledger%escaped_to_infinity(species_idx) = ledger%escaped_to_infinity(species_idx) + macro_charge
        ledger%escaped_count(species_idx) = ledger%escaped_count(species_idx) + 1_i64
      else if (soft_discarded_boundary_flag(i) .or. pcls_batch%alive(i)) then
        ledger%discarded_unresolved(species_idx) = ledger%discarded_unresolved(species_idx) + macro_charge
        ledger%discarded_unresolved_count(species_idx) = ledger%discarded_unresolved_count(species_idx) + 1_i64
      end if
    end do
  end subroutine record_batch_outcome_charge

  !> species 別 flux/count だけを MPI-global 値へ集約する。stock は全 rank で同じ mesh state から得る。
  subroutine reduce_charge_ledger_fluxes(ledger, mpi, workspace)
    type(charge_ledger_type), intent(inout) :: ledger
    type(mpi_context), intent(in) :: mpi
    type(simulator_batch_workspace_type), intent(inout) :: workspace
    integer :: n

    n = int(ledger%nspecies)
    workspace%ledger_charge_values(1:n) = ledger%injected_from_remote
    workspace%ledger_charge_values(n + 1:2*n) = ledger%emitted_from_surface
    workspace%ledger_charge_values(2*n + 1:3*n) = ledger%absorbed_on_surface
    workspace%ledger_charge_values(3*n + 1:4*n) = ledger%escaped_to_infinity
    workspace%ledger_charge_values(4*n + 1:5*n) = ledger%discarded_unresolved
    workspace%ledger_charge_values(5*n + 1:6*n) = ledger%interface_outward_gross
    workspace%ledger_charge_values(6*n + 1:7*n) = ledger%interface_returned_gross
    call mpi_allreduce_sum_real_dp_array(mpi, workspace%ledger_charge_values)
    ledger%injected_from_remote = workspace%ledger_charge_values(1:n)
    ledger%emitted_from_surface = workspace%ledger_charge_values(n + 1:2*n)
    ledger%absorbed_on_surface = workspace%ledger_charge_values(2*n + 1:3*n)
    ledger%escaped_to_infinity = workspace%ledger_charge_values(3*n + 1:4*n)
    ledger%discarded_unresolved = workspace%ledger_charge_values(4*n + 1:5*n)
    ledger%interface_outward_gross = workspace%ledger_charge_values(5*n + 1:6*n)
    ledger%interface_returned_gross = workspace%ledger_charge_values(6*n + 1:7*n)

    workspace%ledger_count_values(1:n) = ledger%injected_count
    workspace%ledger_count_values(n + 1:2*n) = ledger%emitted_count
    workspace%ledger_count_values(2*n + 1:3*n) = ledger%absorbed_count
    workspace%ledger_count_values(3*n + 1:4*n) = ledger%escaped_count
    workspace%ledger_count_values(4*n + 1:5*n) = ledger%discarded_unresolved_count
    call mpi_allreduce_sum_i64_array(mpi, workspace%ledger_count_values)
    ledger%injected_count = workspace%ledger_count_values(1:n)
    ledger%emitted_count = workspace%ledger_count_values(n + 1:2*n)
    ledger%absorbed_count = workspace%ledger_count_values(2*n + 1:3*n)
    ledger%escaped_count = workspace%ledger_count_values(3*n + 1:4*n)
    ledger%discarded_unresolved_count = workspace%ledger_count_values(4*n + 1:5*n)
  end subroutine reduce_charge_ledger_fluxes

end submodule bem_simulator_loop

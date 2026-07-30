!> `bem_simulator` の主ループと粒子処理計算を実装する submodule。
submodule(bem_simulator) bem_simulator_loop
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use bem_app_config_runtime, only: compute_z_high_box_potential_statistics
  use bem_periodic_zero_mode_plan, only: periodic_zero_mode_state_type
  use bem_performance_profile, only: perf_region_batch_total, perf_region_begin, perf_region_commit_charge, &
                                     perf_region_count_outcomes, perf_region_end, perf_region_field_refresh, &
                                     perf_region_field_solver_init, perf_region_history_write, perf_region_mpi_reduce, &
                                     perf_region_particle_batch, perf_region_prepare_batch, perf_region_simulation_total, &
                                     perf_region_stats_update
  implicit none
  ! neutral_return は、少数の長寿命粒子だけを解決済み帰還分布へ繰り込む近似である。
  real(dp), parameter :: neutral_return_max_unresolved_fraction = 0.05_dp
contains

  !> 吸着モデルのバッチループを実行し、電荷更新と統計集計を進める。
  module procedure run_absorption_insulator
  integer(i32), parameter :: adaptive_max_halvings = 24_i32
  real(dp), parameter :: adaptive_acceptance_roundoff_factor = 64.0_dp
  integer(i32) :: batch_idx, final_batch_idx, batch_count_this_run, local_batch_idx, nth, hist_stride
  integer(i32) :: fresh_particle_count, due_escape_count_total, species_idx, trial_halvings
  integer(i32) :: adaptive_thread_count_sum, adaptive_thread_count_mismatch
  integer(i64) :: outer_queue_event_count_before, outer_queue_event_count_after_pop, outer_queue_event_count_after
  integer(i64) :: batch_rejected_trials
  integer(i32) :: collision_failure_count, collision_failure_rank, collision_failure_status
  integer(i32) :: collision_failure_particle, collision_failure_step
  integer(i32) :: local_failure_values(3), selected_failure_values(3)
  integer(i32) :: photo_failure_count, photo_failure_rank, photo_failure_status, photo_failure_species
  integer(i32) :: photo_failure_ray, photo_failure_bounce
  integer(i32) :: photo_local_failure_values(4), photo_selected_failure_values(4)
  integer :: hist_unit
  integer :: rng_state_size
  integer, allocatable :: rng_state_before(:)
  logical :: history_enabled
  logical :: potential_history_enabled
  logical :: top_reference_history_enabled
  integer :: pot_hist_unit
  integer :: top_ref_hist_unit
!$ integer(kind=omp_sched_kind) :: previous_omp_schedule
!$ integer :: previous_omp_schedule_chunk
!$ logical :: previous_omp_dynamic
  real(dp), allocatable :: potential_buf(:), injection_residual_before(:)
  integer(i32) :: batch_counts(6)
  real(dp) :: bfield(3), rel, t0, sim_t0, batch_t0, batch_soft_discarded_abs_charge
  real(dp) :: trial_batch_duration, duration_ratio, adaptive_potential_step
  real(dp) :: adaptive_metric_values(1), committed_gauss_residual
  real(dp) :: collision_failure_x(3), collision_failure_v(3), selected_failure_state(6)
  real(dp) :: max_outer_flight_time, max_frozen_field_ratio, max_outer_energy_relative_error
  real(dp) :: outer_max_diagnostics(3)
  real(dp) :: outer_queue_charge_before, outer_queue_charge_after_pop, outer_queue_charge_after
  real(dp) :: outer_queue_photoelectron_number_before, outer_queue_photoelectron_number_after_pop
  real(dp) :: outer_queue_photoelectron_number_after
  real(dp) :: outer_queue_area_xy, outer_queue_photoelectron_column_target
  character(len=16) :: outer_queue_fingerprint
  character(len=256) :: implicit_mean_failure_message
  type(particles_soa) :: pcls_batch
  type(outer_event_record_type), allocatable :: due_outer_events(:)
  real(dp), allocatable :: due_returned_charge(:), due_escaped_charge(:)
  integer(i64), allocatable :: due_escaped_count(:)
  type(mpi_context) :: mpi_ctx
  type(electrostatic_snapshot_type) :: snapshot
  type(electrostatic_diagnostics_type) :: committed_snapshot_diagnostics
  type(periodic_zero_mode_state_type) :: committed_zero_state
  type(outer_coupler_type) :: outer_coupler
  type(charge_ledger_type) :: batch_ledger
  type(simulator_batch_workspace_type) :: workspace
  type(particle_source_plan_type) :: source_plan
  logical :: ledger_enabled
  logical :: outer_queue_enabled
  integer(i32) :: kinetic_status
  integer(i32) :: boundary_status
  character(len=256) :: kinetic_message
  character(len=256) :: boundary_message
  type(kinetic_outer_plasma_options_type) :: kinetic_options
  type(kinetic_outer_plasma_options_type) :: committed_kinetic_options
  type(dynamic_k0_step_type) :: dynamic_k0_step
  type(external_boundary_contract_type) :: boundary_contract
  type(outer_plasma_state_type) :: steady_start_state, committed_outer
  type(app_config) :: trial_app
  real(dp) :: steady_start_charge
  real(dp) :: dynamic_k0_area_xy, mean_transaction_residual
  real(dp) :: tracked_photoelectron_outward_current_density
  real(dp) :: mean_sample_escape_fraction, mean_return_weight_scale
  real(dp) :: mean_sampled_deferred_absorbed_charge, mean_sampled_deferred_escaped_charge
  real(dp) :: top_phi_mean, top_phi_std, top_phi_min, top_phi_max
  logical :: implicit_mean_enabled, adaptive_nonzero_mode, trial_accepted
  logical :: implicit_mean_retryable_failure
  logical :: restored_outer_snapshot

  stats = sim_stats()
  if (present(initial_stats)) stats = initial_stats
  mpi_ctx = mpi_context()
  if (present(mpi)) mpi_ctx = mpi
  call resolve_external_boundary_contract( &
    app%sim%reservoir_potential_model, app%sim%open_boundary_model, &
    app%outer_plasma%model, app%outer_plasma%kinetic_closure, app%outer_plasma%return_model, &
    app%coupling%particle_transfer_mode, app%coupling%outer_queue_enabled, boundary_contract, boundary_status, &
    boundary_message &
    )
  if (boundary_status /= external_boundary_ok) error stop trim(boundary_message)
  ledger_enabled = present(charge_ledger)
  outer_queue_enabled = app%coupling%outer_queue_enabled
  implicit_mean_enabled = trim(lower_ascii(app%coupling%update_mode)) == 'implicit_mean'
  adaptive_nonzero_mode = app%periodic2%max_nonzero_mode_potential_step > 0.0_dp
  nth = 1
!$ if (adaptive_nonzero_mode) then
!$  previous_omp_dynamic = omp_get_dynamic()
!$  call omp_set_dynamic(.false.)
!$omp parallel default(shared)
!$omp single
!$  nth = omp_get_num_threads()
!$omp end single
!$omp end parallel
!$ else
!$  nth = max(1, omp_get_max_threads())
!$ end if
  if (adaptive_nonzero_mode) then
    adaptive_thread_count_sum = nth
    call mpi_allreduce_sum_i32_scalar(mpi_ctx, adaptive_thread_count_sum)
    adaptive_thread_count_mismatch = merge( &
                                     1_i32, 0_i32, &
                                     adaptive_thread_count_sum /= nth*mpi_ctx%size &
                                     )
    call mpi_allreduce_sum_i32_scalar(mpi_ctx, adaptive_thread_count_mismatch)
    if (adaptive_thread_count_mismatch > 0_i32) then
      error stop 'adaptive nonzero-mode requires the same OpenMP team size on every MPI rank.'
    end if
    if (trim(lower_ascii(app%periodic2%nonzero_mode_backend)) /= 'cached_kneq0') then
      error stop 'adaptive nonzero-mode runtime requires nonzero_mode_backend="cached_kneq0".'
    end if
    if (outer_queue_enabled) then
      error stop 'adaptive nonzero-mode runtime cannot mutate an outer event queue.'
    end if
    do species_idx = 1_i32, app%n_particle_species
      if (.not. app%particle_species(species_idx)%enabled) cycle
      select case (trim(lower_ascii(app%particle_species(species_idx)%source_mode)))
      case ('volume_seed')
        error stop 'adaptive nonzero-mode runtime does not support volume_seed sources.'
      case ('reservoir_face')
        if (.not. app%particle_species(species_idx)%has_target_macro_particles_per_batch) then
          error stop 'adaptive nonzero-mode reservoir_face requires target_macro_particles_per_batch.'
        end if
      case ('photo_raycast')
        continue
      case default
        error stop 'adaptive nonzero-mode runtime received an unsupported particle source.'
      end select
    end do
    if (stats%batches > 0_i32 .and. &
        (stats%simulated_time <= 0.0_dp .or. stats%adaptive_nonzero_mode_last_batch_duration <= 0.0_dp)) then
      error stop 'adaptive nonzero-mode resume requires accepted physical-time checkpoint state.'
    end if
    if (stats%batches > 0_i32) then
      if (stats%adaptive_nonzero_mode_omp_threads <= 0_i32) then
        error stop 'adaptive nonzero-mode resume requires a checkpointed OpenMP thread count.'
      end if
      if (stats%adaptive_nonzero_mode_omp_threads /= nth) then
        if (mpi_is_root(mpi_ctx)) then
          write (error_unit, '(a,i0,a,i0)') &
            'adaptive nonzero-mode resume thread mismatch: checkpoint=', &
            stats%adaptive_nonzero_mode_omp_threads, ' current=', nth
        end if
        error stop 'adaptive nonzero-mode resume requires the original OpenMP thread count.'
      end if
    else
      stats%adaptive_nonzero_mode_omp_threads = nth
    end if
  end if
  if (.not. adaptive_nonzero_mode .and. stats%batches > 0_i32 .and. stats%simulated_time == 0.0_dp) then
    ! Versioned checkpoints written before physical time became explicit used a
    ! fixed batch duration, so their accepted time can be reconstructed exactly.
    stats%simulated_time = real(stats%batches, dp)*app%sim%batch_duration
  end if
  dynamic_k0_step = dynamic_k0_step_type()
  if (implicit_mean_enabled .and. outer_queue_enabled) then
    error stop 'implicit_mean cannot be combined with an outer event queue.'
  end if
  dynamic_k0_area_xy = 0.0_dp
  if (implicit_mean_enabled) then
    dynamic_k0_area_xy = (app%sim%box_max(1) - app%sim%box_min(1))* &
                         (app%sim%box_max(2) - app%sim%box_min(2))
    if (.not. ieee_is_finite(dynamic_k0_area_xy) .or. dynamic_k0_area_xy <= 0.0_dp) then
      error stop 'implicit_mean requires a positive finite horizontal box area.'
    end if
  end if
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

  call workspace%init( &
    mesh%nelem, app%n_particle_species, nth, implicit_mean_enabled=implicit_mean_enabled, &
    candidate_charge_enabled=adaptive_nonzero_mode &
    )
!$ call omp_get_schedule(previous_omp_schedule, previous_omp_schedule_chunk)
!$ if (adaptive_nonzero_mode) then
!$  call omp_set_schedule(omp_sched_static, 1)
!$ else
!$  call omp_set_schedule(omp_sched_dynamic, 1)
!$ end if
  if (adaptive_nonzero_mode) then
    call random_seed(size=rng_state_size)
    allocate (rng_state_before(rng_state_size))
    if (present(inject_state)) then
      if (allocated(inject_state%macro_residual)) then
        allocate (injection_residual_before(size(inject_state%macro_residual)))
      end if
    end if
  end if
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
  top_reference_history_enabled = present(top_reference_history_unit)
  top_ref_hist_unit = -1
  if (top_reference_history_enabled) then
    top_ref_hist_unit = top_reference_history_unit
    top_reference_history_enabled = top_ref_hist_unit /= -1
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
  restored_outer_snapshot = .false.
  if (present(electrostatic_restart_state)) then
    restored_outer_snapshot = electrostatic_restart_state%outer_ready
    call snapshot%restore_outer_state( &
      electrostatic_restart_state, require_charge_consistency=implicit_mean_enabled &
      )
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
    if (implicit_mean_enabled .and. restored_outer_snapshot .and. snapshot%outer%ready) then
      call outer_coupler%accept_restored_snapshot()
    end if
  else
    call outer_coupler%init(app%coupling)
  end if
  call perf_region_end(perf_region_field_solver_init, t0)
  trial_app = app

  do local_batch_idx = 1, batch_count_this_run
    call perf_region_begin(perf_region_batch_total, batch_t0)
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

    trial_batch_duration = app%sim%batch_duration
    trial_halvings = 0_i32
    batch_rejected_trials = 0_i64
    adaptive_potential_step = 0.0_dp
    trial_accepted = .false.
    if (adaptive_nonzero_mode) then
      call random_seed(get=rng_state_before)
      if (present(inject_state)) then
        if (allocated(injection_residual_before)) injection_residual_before = inject_state%macro_residual
      end if
      committed_outer = snapshot%outer
      committed_kinetic_options = snapshot%kinetic_options
      committed_zero_state = snapshot%zero_state
      committed_snapshot_diagnostics = snapshot%diagnostics
      committed_gauss_residual = snapshot%gauss_residual
    end if

    do while (.not. trial_accepted)
      if (adaptive_nonzero_mode) then
        ! Every retry starts from the exact accepted batch boundary.  In particular,
        ! rejected source sampling must not consume RNG draws or macro residual.
        call random_seed(put=rng_state_before)
        if (present(inject_state)) then
          if (allocated(injection_residual_before)) inject_state%macro_residual = injection_residual_before
        end if
        snapshot%outer = committed_outer
        snapshot%kinetic_options = committed_kinetic_options
        snapshot%zero_state = committed_zero_state
        snapshot%diagnostics = committed_snapshot_diagnostics
        snapshot%gauss_residual = committed_gauss_residual
        dynamic_k0_step = dynamic_k0_step_type()
      end if

      if (adaptive_nonzero_mode) trial_app = app
      trial_app%sim%batch_duration = trial_batch_duration
      duration_ratio = trial_batch_duration/app%sim%batch_duration
      do species_idx = 1_i32, trial_app%n_particle_species
        if (.not. trial_app%particle_species(species_idx)%enabled) cycle
        if (trim(lower_ascii(trial_app%particle_species(species_idx)%source_mode)) /= 'reservoir_face') cycle
        if (.not. trial_app%particle_species(species_idx)%has_target_macro_particles_per_batch) cycle
        trial_app%particle_species(species_idx)%w_particle = &
          app%particle_species(species_idx)%w_particle*duration_ratio
      end do

      call perf_region_begin(perf_region_prepare_batch, t0)
      if (adaptive_nonzero_mode .or. local_batch_idx == 1_i32) then
        call build_particle_source_plan(trial_app, source_plan, mpi=mpi_ctx)
      end if
      call prepare_batch_state( &
        mesh, trial_app, source_plan, snapshot, stats, batch_idx, workspace, pcls_batch, mpi_ctx, snapshot%outer, &
        inject_state, photo_failure_status, photo_failure_species, photo_failure_ray, photo_failure_bounce &
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
        call workspace%prepare_particle_flags(pcls_batch%n, outer_queue_enabled=.true., implicit_mean_enabled=.false.)
      end if

      call perf_region_begin(perf_region_particle_batch, t0)
      call process_particle_batch( &
        mesh, trial_app, boundary_contract, snapshot, pcls_batch, workspace%dq_thread, workspace%escaped_boundary_flag, &
        workspace%absorbed_flag, workspace%absorbed_element, bfield, batch_idx, mpi_ctx%rank, &
        workspace%soft_discarded_boundary_flag, workspace%queued_outer_flag, workspace%outer_event_staging, &
        workspace%deferred_mean_interface_flag, workspace%deferred_mean_interface_step, &
        workspace%deferred_mean_interface_crossing, &
        workspace%interface_outward_thread, workspace%interface_returned_thread, &
        collision_failure_status, collision_failure_particle, &
        collision_failure_step, collision_failure_x, collision_failure_v, workspace%interface_tau_max_thread, &
        workspace%interface_frozen_ratio_max_thread, workspace%interface_energy_error_max_thread &
        )
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
          selected_failure_values(3), trial_app%sim%dt, selected_failure_state(1:3), selected_failure_state(4:6) &
          )
      end if

      if (outer_queue_enabled) then
        call enqueue_staged_outer_events( &
          outer_queue_state, workspace%queued_outer_flag, workspace%outer_event_staging, pcls_batch%n &
          )
        call measure_outer_queue_state( &
          trial_app, outer_queue_state, mpi_ctx, outer_queue_charge_after, outer_queue_photoelectron_number_after, &
          outer_queue_event_count_after &
          )
        workspace%interface_returned_thread(:, 1) = &
          workspace%interface_returned_thread(:, 1) + due_returned_charge
      end if

      tracked_photoelectron_outward_current_density = 0.0_dp
      mean_sample_escape_fraction = 0.0_dp
      mean_return_weight_scale = 0.0_dp
      mean_transaction_residual = 0.0_dp
      mean_sampled_deferred_absorbed_charge = 0.0_dp
      mean_sampled_deferred_escaped_charge = 0.0_dp
      implicit_mean_failure_message = ''
      implicit_mean_retryable_failure = .false.
      if (implicit_mean_enabled) then
        call resolve_implicit_mean_batch( &
          mesh, trial_app, boundary_contract, snapshot, pcls_batch, workspace, bfield, mpi_ctx, &
          dynamic_k0_step, tracked_photoelectron_outward_current_density, mean_sample_escape_fraction, &
          mean_return_weight_scale, mean_transaction_residual, &
          mean_sampled_deferred_absorbed_charge, mean_sampled_deferred_escaped_charge, &
          collision_failure_status, collision_failure_particle, collision_failure_step, &
          collision_failure_x, collision_failure_v, implicit_mean_failure_message &
          )

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
            selected_failure_values(3), trial_app%sim%dt, selected_failure_state(1:3), selected_failure_state(4:6) &
            )
        end if
        if (dynamic_k0_step%status /= dynamic_k0_ok) then
          implicit_mean_retryable_failure = adaptive_nonzero_mode .and. &
                                            dynamic_k0_step%status == &
                                            dynamic_zhao_frozen_cohort_trust_failure
          if (.not. implicit_mean_retryable_failure) then
            if (mpi_is_root(mpi_ctx)) then
              write (error_unit, '(a,i0,2a)') &
                'implicit mean trial failed with status=', dynamic_k0_step%status, ': ', &
                trim(implicit_mean_failure_message)
            end if
            error stop 'implicit mean k0 update failed without a shorter-step recovery path.'
          end if
        end if
      end if

      call apply_neutral_return_surface_closure( &
        trial_app, pcls_batch, fresh_particle_count, workspace, mpi_ctx &
        )

      if (implicit_mean_retryable_failure) then
        trial_accepted = .false.
      else if (adaptive_nonzero_mode) then
        call prepare_adaptive_charge_candidate(mesh, workspace, implicit_mean_enabled, mpi_ctx)
        adaptive_metric_values = 0.0_dp
        if (mpi_is_root(mpi_ctx)) then
          call snapshot%measure_kneq0_potential_step( &
            mesh, workspace%mean_candidate_charge, adaptive_metric_values(1) &
            )
        end if
        call mpi_allreduce_max_real_dp_array(mpi_ctx, adaptive_metric_values)
        adaptive_potential_step = adaptive_metric_values(1)
        trial_accepted = adaptive_potential_step <= &
          max(0.0_dp, app%periodic2%max_nonzero_mode_potential_step - &
              adaptive_acceptance_roundoff_factor*epsilon(1.0_dp)*max( &
              1.0_dp, adaptive_potential_step, app%periodic2%max_nonzero_mode_potential_step &
              ))
      else
        trial_accepted = .true.
      end if
      if (trial_accepted) exit

      batch_rejected_trials = batch_rejected_trials + 1_i64
      workspace%charge_candidate_ready = .false.
      if (mpi_is_root(mpi_ctx)) then
        if (implicit_mean_retryable_failure) then
          write (output_unit, '(a,i0,a,i0,a,es13.5,a,i0,2a)') &
            'BEACH adaptive-kneq0 reject: batch=', batch_idx, ' halving=', trial_halvings, &
            ' duration_s=', trial_batch_duration, ' implicit_status=', dynamic_k0_step%status, &
            ' reason=', trim(implicit_mean_failure_message)
        else
          write (output_unit, '(a,i0,a,i0,3(a,es13.5))') &
            'BEACH adaptive-kneq0 reject: batch=', batch_idx, ' halving=', trial_halvings, &
            ' duration_s=', trial_batch_duration, ' max_delta_phi_V=', adaptive_potential_step, &
            ' limit_V=', app%periodic2%max_nonzero_mode_potential_step
        end if
        flush (output_unit)
      end if
      if (trial_halvings >= adaptive_max_halvings) then
        error stop 'adaptive nonzero-mode batch failed after 24 duration halvings.'
      end if
      trial_halvings = trial_halvings + 1_i32
      trial_batch_duration = scale(app%sim%batch_duration, -trial_halvings)
      if (.not. ieee_is_finite(trial_batch_duration) .or. trial_batch_duration <= 0.0_dp) then
        error stop 'adaptive nonzero-mode batch duration became invalid.'
      end if
    end do

    if (adaptive_nonzero_mode .and. mpi_is_root(mpi_ctx)) then
      write (output_unit, '(a,i0,3(a,es13.5),a,i0)') &
        'BEACH adaptive-kneq0 accept: batch=', batch_idx, &
        ' time_end_s=', stats%simulated_time + trial_batch_duration, &
        ' duration_s=', trial_batch_duration, ' max_delta_phi_V=', adaptive_potential_step, &
        ' halvings=', trial_halvings
      flush (output_unit)
    end if

    max_outer_flight_time = max(max_outer_flight_time, maxval(workspace%interface_tau_max_thread))
    max_frozen_field_ratio = max(max_frozen_field_ratio, maxval(workspace%interface_frozen_ratio_max_thread))
    max_outer_energy_relative_error = &
      max(max_outer_energy_relative_error, maxval(workspace%interface_energy_error_max_thread))

    if (ledger_enabled) then
      call batch_ledger%reset(batch_idx)
      batch_ledger%surface_charge_before = sum(mesh%q_elem)
      batch_ledger%outer_flight_charge_before = outer_queue_charge_before
      call record_batch_initial_charge(trial_app, pcls_batch, fresh_particle_count, batch_ledger)
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
      ! neutral_return diagnostics are already MPI-global.  Assign them after the
      ! rank-local flux reduction so they are neither summed nor rank-multiplied.
      batch_ledger%neutral_return_correction = workspace%neutral_return_correction
      batch_ledger%neutral_return_weight_scale = workspace%neutral_return_weight_scale
      batch_ledger%neutral_return_unresolved_fraction = workspace%neutral_return_unresolved_fraction
      if (implicit_mean_enabled) then
        call reconcile_implicit_mean_charge_ledger( &
          trial_app, dynamic_k0_step, dynamic_k0_area_xy, trial_batch_duration, &
          tracked_photoelectron_outward_current_density, &
          mean_sampled_deferred_absorbed_charge, mean_sampled_deferred_escaped_charge, batch_ledger &
          )
      end if
    end if

    call perf_region_begin(perf_region_commit_charge, t0)
    call commit_batch_charge( &
      mesh, app%sim%q_floor, app%sim%e0, app%sim%field_bc_mode, &
      workspace, rel, mpi_ctx &
      )
    if (implicit_mean_enabled) call outer_coupler%mark_mesh_changed()
    call perf_region_end(perf_region_commit_charge, t0)

    if (implicit_mean_enabled) then
      call perf_region_begin(perf_region_field_refresh, t0)
      if (index('ABC', dynamic_k0_step%zhao_branch) > 0) then
        ! resolve_implicit_mean_batch が受理済みのZhao profileを同じ候補電荷へ
        ! collectiveに採用済み。ここではcommitted k/=0場だけを更新する。
        call snapshot%refresh(mesh, update_outer=.false.)
        call outer_coupler%accept_external_update(batch_idx)
      else
        call outer_coupler%refresh( &
          snapshot, mesh, batch_idx, continuation_stage='post_implicit_mean', force_outer_update=.true. &
          )
      end if
      dynamic_k0_step%surface_charge_after_c = sum(mesh%q_elem)
      dynamic_k0_step%interface_potential_after_v = snapshot%outer%interface_potential
      dynamic_k0_step%interface_field_after_v_m = snapshot%outer%interface_field
      snapshot%outer%electron_current_density = dynamic_k0_step%electron_current_density_a_m2
      snapshot%outer%ion_current_density = dynamic_k0_step%ion_current_density_a_m2
      snapshot%outer%photoelectron_current_density = dynamic_k0_step%photoelectron_current_density_a_m2
      snapshot%outer%total_current_density = dynamic_k0_step%total_current_density_a_m2
      call perf_region_end(perf_region_field_refresh, t0)
    end if

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
        snapshot, mesh, batch_idx, continuation_stage='post_enqueue', force_outer_update=.true. &
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
    stats%simulated_time = stats%simulated_time + trial_batch_duration
    if (adaptive_nonzero_mode) then
      stats%adaptive_nonzero_mode_rejected_trials = &
        stats%adaptive_nonzero_mode_rejected_trials + batch_rejected_trials
      stats%adaptive_nonzero_mode_last_batch_duration = trial_batch_duration
      stats%adaptive_nonzero_mode_last_potential_step = adaptive_potential_step
    end if
    call perf_region_end(perf_region_stats_update, t0)

    call perf_region_begin(perf_region_history_write, t0)
    if (mpi_is_root(mpi_ctx)) then
      call print_batch_progress(batch_idx, final_batch_idx, rel)
      if (implicit_mean_enabled .and. &
          (batch_idx <= 10_i32 .or. mod(batch_idx, hist_stride) == 0_i32 .or. batch_idx == final_batch_idx)) then
        write (output_unit, '(a,i0,16(a,es13.5),a,i0)') &
          'BEACH implicit-mean batch=', batch_idx, &
          ' charge_C=', dynamic_k0_step%surface_charge_after_c, &
          ' phi_I_V=', dynamic_k0_step%interface_potential_after_v, &
          ' phi_explicit_V=', dynamic_k0_step%explicit_predictor_potential_v, &
          ' E_I_V_m=', dynamic_k0_step%interface_field_after_v_m, &
          ' J_e_A_m2=', dynamic_k0_step%electron_current_density_a_m2, &
          ' J_i_A_m2=', dynamic_k0_step%ion_current_density_a_m2, &
          ' J_pe_out_A_m2=', tracked_photoelectron_outward_current_density, &
          ' J_pe_A_m2=', dynamic_k0_step%photoelectron_current_density_a_m2, &
          ' J_other_A_m2=', dynamic_k0_step%additional_tracked_current_density_a_m2, &
          ' f_escape=', dynamic_k0_step%photoelectron_escape_fraction, &
          ' f_return=', dynamic_k0_step%photoelectron_return_fraction, &
          ' sample_escape_fraction=', mean_sample_escape_fraction, &
          ' return_weight_scale=', mean_return_weight_scale, &
          ' transaction_residual_C=', mean_transaction_residual, &
          ' shadow_return_tau_mean_s=', dynamic_k0_step%returned_outer_flight_time_mean_s, &
          ' shadow_return_column_charge_C_m2=', &
          dynamic_k0_step%estimated_returning_photoelectron_column_charge_per_area_c_m2, &
          ' mean_solver_iterations=', dynamic_k0_step%iterations
        if (index('ABC', dynamic_k0_step%zhao_branch) > 0) then
          write (output_unit, '(a,i0,2(a,a),9(a,es13.5))') &
            'BEACH implicit-zhao batch=', batch_idx, &
            ' branch=', dynamic_k0_step%zhao_branch, &
            ' closure=', trim(snapshot%outer%kinetic_closure), &
            ' phi_min_V=', snapshot%outer%zhao_phi_minimum, &
            ' n_e_inf_m3=', snapshot%outer%zhao_electron_density_infinity, &
            ' barrier_J=', dynamic_k0_step%photoelectron_barrier_energy_j, &
            ' source_scale=', dynamic_k0_step%zhao_effective_source_scale, &
            ' marginal_energy_J=', dynamic_k0_step%marginal_photoelectron_energy_j, &
            ' marginal_escape_fraction=', dynamic_k0_step%marginal_photoelectron_escape_fraction, &
            ' empirical_charge_residual_C=', dynamic_k0_step%backward_euler_residual_charge_c, &
            ' recross_charge_fraction=', dynamic_k0_step%photoelectron_recross_charge_fraction, &
            ' terminal_mismatch_charge_fraction=', &
            dynamic_k0_step%photoelectron_terminal_mismatch_charge_fraction
        end if
      end if
      if (outer_queue_enabled) then
        write (output_unit, '(a,i0,a,i0,a,i0,a,es13.5,a,es13.5,a,es13.5,a,a,a,es13.5)') &
          'BEACH outer-queue batch=', batch_idx, ' events_after_pop=', outer_queue_event_count_after_pop, &
          ' events_for_closure=', outer_queue_event_count_after, ' charge_after_C=', outer_queue_charge_after, &
          ' photoelectron_column_target_m-2=', outer_queue_photoelectron_column_target, &
          ' eta=', snapshot%outer%photoelectron_population_fraction, ' branch=', snapshot%outer%zhao_branch, &
          ' column_residual_m-2=', snapshot%outer%photoelectron_column_residual_per_area
      end if
      if (batch_idx <= 10_i32 .or. mod(batch_idx, hist_stride) == 0_i32 .or. &
          batch_idx == final_batch_idx) then
        do species_idx = 1_i32, trial_app%n_particle_species
          if (.not. trial_app%particle_species(species_idx)%enabled) cycle
          if (trim(lower_ascii(trial_app%particle_species(species_idx)%surface_charge_closure)) /= &
              'neutral_return') cycle
          write (output_unit, '(a,i0,a,i0,6(a,es13.5))') &
            'BEACH neutral-return batch=', batch_idx, ' species=', species_idx, &
            ' emitted_charge_C=', workspace%neutral_return_emitted_charge(species_idx), &
            ' absorbed_charge_C=', workspace%neutral_return_absorbed_charge(species_idx), &
            ' unresolved_charge_C=', workspace%neutral_return_unresolved_charge(species_idx), &
            ' weight_scale=', workspace%neutral_return_weight_scale(species_idx), &
            ' unresolved_fraction=', workspace%neutral_return_unresolved_fraction(species_idx), &
            ' correction_C=', workspace%neutral_return_correction(species_idx)
          if (workspace%neutral_return_weight_scale(species_idx) > 1.05_dp) then
            write (error_unit, '(a,i0,a,i0,2(a,es13.5))') &
              'BEACH neutral-return warning: batch=', batch_idx, ' species=', species_idx, &
              ' weight_scale=', workspace%neutral_return_weight_scale(species_idx), &
              ' unresolved_fraction=', workspace%neutral_return_unresolved_fraction(species_idx)
          end if
        end do
      end if
      call maybe_write_history_snapshot(history_enabled, hist_unit, hist_stride, stats, rel, mesh%q_elem)
      if (potential_history_enabled) then
        call maybe_write_potential_history_snapshot( &
          potential_history_enabled, pot_hist_unit, hist_stride, stats, &
          snapshot, mesh, app%sim, potential_buf, top_reference_history_enabled, top_ref_hist_unit &
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
    if (implicit_mean_enabled .and. dynamic_k0_step%status == dynamic_k0_ok) then
      electrostatic_diagnostics%implicit_mean_shadow_diagnostics_available = .true.
      electrostatic_diagnostics%implicit_mean_last_returned_outer_flight_time_mean = &
        dynamic_k0_step%returned_outer_flight_time_mean_s
      electrostatic_diagnostics%implicit_mean_last_returning_pe_column_charge_per_area = &
        dynamic_k0_step%estimated_returning_photoelectron_column_charge_per_area_c_m2
    end if
    if (outer_queue_enabled) then
      electrostatic_diagnostics%outer_queue_event_count = outer_queue_event_count_after
      electrostatic_diagnostics%outer_queue_signed_charge = outer_queue_charge_after
      electrostatic_diagnostics%outer_queue_fingerprint = outer_queue_fingerprint
    end if
    if ((app%write_mesh_potential .or. app%write_potential_history) .and. &
        app%sim%use_box .and. mpi_is_root(mpi_ctx)) then
      call compute_z_high_box_potential_statistics( &
        mesh, app%sim, snapshot, top_phi_mean, top_phi_std, top_phi_min, top_phi_max &
        )
      electrostatic_diagnostics%top_reference_available = .true.
      electrostatic_diagnostics%top_reference_last_batch = stats%batches
      electrostatic_diagnostics%top_reference_simulated_time = stats%simulated_time
      electrostatic_diagnostics%top_reference_z_high = app%sim%box_max(3)
      electrostatic_diagnostics%top_reference_sample_n = app%sim%injection_face_phi_grid_n
      electrostatic_diagnostics%top_reference_potential_mean = top_phi_mean
      electrostatic_diagnostics%top_reference_potential_std = top_phi_std
      electrostatic_diagnostics%top_reference_potential_min = top_phi_min
      electrostatic_diagnostics%top_reference_potential_max = top_phi_max
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
!$ call omp_set_schedule(previous_omp_schedule, previous_omp_schedule_chunk)
!$ if (adaptive_nonzero_mode) call omp_set_dynamic(previous_omp_dynamic)
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
  call workspace%prepare_particle_flags( &
    pcls_batch%n, outer_queue_enabled=app%coupling%outer_queue_enabled, &
    implicit_mean_enabled=trim(lower_ascii(app%coupling%update_mode)) == 'implicit_mean' &
    )
  end procedure prepare_batch_state

  !> 粒子を時間発展させ、衝突時の堆積電荷をスレッド別に集計する。
  module procedure process_particle_batch
  integer(i32) :: i, step, tid, nth, warn_stride, collision_status
  real(dp) :: x0(3), v0(3), x1(3), v1(3), qdep
  type(hit_info) :: hit
  type(particle_step_result) :: step_result
  type(particle_step_result) :: external_final_result
  type(external_step_trace_type) :: external_trace
  type(sim_config) :: particle_sim
  logical :: has_warn_stride, collision_failed, candidate_inside, used_event_resolver
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

  !$omp parallel default(none) &
  !$omp shared(mesh,pcls_batch,app,boundary_contract,snapshot,dq_thread,bfield,escaped_boundary_flag,absorbed_flag,nth) &
  !$omp shared(absorbed_element) &
  !$omp shared(soft_discarded_boundary_flag,queued_outer_flag,outer_event_staging) &
  !$omp shared(deferred_mean_interface_flag,deferred_mean_interface_step,deferred_mean_interface_crossing) &
  !$omp shared(interface_outward_thread,interface_returned_thread) &
  !$omp shared(interface_tau_max_thread,interface_frozen_ratio_max_thread) &
  !$omp shared(interface_energy_error_max_thread) &
  !$omp shared(warn_stride,batch_idx,mpi_rank,collision_failure_status,collision_failure_particle,collision_failure_step) &
  !$omp shared(collision_failure_x,collision_failure_v) &
  !$omp private(i,step,x0,v0,x1,v1,hit,step_result,external_final_result,external_trace,particle_sim,tid,qdep,species_idx) &
  !$omp private(external_event) &
  !$omp private(collision_status,collision_failed,candidate_inside,used_event_resolver)
  ! スレッドごとに dq_thread(:, tid) を使って原子的更新なしで電荷を集める。
  tid = 1
!$ tid = omp_get_thread_num() + 1
  !$omp do schedule(runtime)
  do i = 1, pcls_batch%n
    if (.not. pcls_batch%alive(i)) cycle
    species_idx = pcls_batch%species_id(i)
    particle_sim = app%sim
    if (trim(lower_ascii(app%particle_species(species_idx)%z_high_boundary)) == 'reflect') then
      particle_sim%bc_high(3) = bc_reflect
    end if
    collision_failed = .false.
    do step = 1, app%sim%max_step
      x0 = pcls_batch%x(:, i)
      v0 = pcls_batch%v(:, i)
      call build_particle_step_candidate( &
        mesh, particle_sim, snapshot, bfield, x0, v0, &
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
      call find_first_hit(mesh, x0, x1, hit, sim=particle_sim, status=collision_status)
      candidate_inside = .not. particle_sim%use_box .or. &
                         (all(x1 > particle_sim%box_min) .and. all(x1 < particle_sim%box_max))
      used_event_resolver = .false.
      if (collision_status /= collision_query_ok) then
        if (.not. candidate_inside) then
          call resolve_particle_boundary_candidate( &
            mesh, particle_sim, snapshot, bfield, x0, v0, pcls_batch%q(i), pcls_batch%m(i), app%sim%dt, x1, v1, &
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
          absorbed_element(i) = hit%elem_idx
          exit
        end if
        pcls_batch%x(:, i) = x1
        pcls_batch%v(:, i) = v1
      else
        call resolve_particle_boundary_candidate( &
          mesh, particle_sim, snapshot, bfield, x0, v0, pcls_batch%q(i), pcls_batch%m(i), app%sim%dt, x1, v1, &
          hit=hit, result=step_result, boundary_contract=boundary_contract &
          )
        used_event_resolver = .true.
      end if
      if (used_event_resolver) then
        if (step_result%interface_crossing%has_crossing) then
          qdep = pcls_batch%q(i)*pcls_batch%w(i)
          if (trim(lower_ascii(app%coupling%update_mode)) == 'implicit_mean' .and. &
              trim(lower_ascii(app%particle_species(species_idx)%source_mode)) == 'photo_raycast' .and. &
              app%particle_species(species_idx)%q_particle < 0.0_dp) then
            ! The first 3D trace ends at the interface.  The implicit mean
            ! closure later splits its macro weight into escape and return;
            ! the returned shadow is then traced under new k=0 and old k/=0.
            ! Deferring that trace avoids applying the previous-batch profile
            ! and the updated mean barrier to the same identity.
            interface_outward_thread(species_idx, tid) = interface_outward_thread(species_idx, tid) + qdep
            deferred_mean_interface_flag(i) = .true.
            deferred_mean_interface_step(i) = step
            deferred_mean_interface_crossing(i) = step_result%interface_crossing
            step_result%x = step_result%interface_crossing%position
            step_result%v = step_result%interface_crossing%velocity
            step_result%interface_crossing%has_crossing = .false.
            step_result%escaped_boundary = .true.
          else
            call continue_external_particle_step( &
              boundary_contract, snapshot, mesh, app%sim, app%coupling, bfield, &
              pcls_batch%q(i), pcls_batch%m(i), app%sim%batch_duration, step_result, external_final_result, &
              external_trace &
              )
            step_result = external_final_result
            do external_event = 1_i32, external_trace%count
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
          absorbed_element(i) = step_result%elem_idx
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

  !> fastな平均総電荷を陰的に解き、slowな帰還先分布をrayで1回だけ標本化する。
  module procedure resolve_implicit_mean_batch
  integer(i32) :: i, elem_idx, transaction_status, failure_count, nth, photoelectron_index
  integer(i32) :: local_energy_count, energy_index, distribution_status
  real(dp) :: area_xy, pending_current, electron_current, ion_current, additional_current
  real(dp) :: macro_charge, source_charge_total, pending_absolute_total
  real(dp) :: target_return_charge, target_escape_charge, raw_return_charge, raw_escape_charge
  real(dp) :: return_fraction, sample_return_fraction, escape_weight_scale, particle_weight_scale
  real(dp) :: mean_charge_scale, mean_charge_tolerance, sample_charge_tolerance
  real(dp) :: vector_charge_scale, transaction_tolerance, boundary_factor, mean_capacitance
  real(dp) :: surface_charge_base, effective_source_scale, resolved_barrier
  real(dp) :: recross_charge, terminal_mismatch_charge, empirical_return_charge
  integer(i32), allocatable :: outward_event_count(:), returned_event_count(:)
  real(dp), allocatable :: returned_flight_time_sum(:)
  real(dp), allocatable :: sample_return_weight(:)
  real(dp), allocatable :: local_energy_j(:), local_energy_charge_c(:)
  real(dp), allocatable :: global_energy_j(:), global_energy_charge_c(:)
  real(dp), allocatable :: mean_tau_max_thread(:), mean_frozen_ratio_max_thread(:), mean_energy_error_max_thread(:)
  real(dp) :: returning_shadow_moments(2)
  type(dynamic_k0_step_type) :: scalar_predictor
  type(outer_plasma_state_type) :: accepted_outer
  type(measured_interface_energy_distribution_type) :: measured_distribution
  type(mean_charge_transaction_type) :: transaction
  character(len=256) :: message
  logical :: use_dynamic_zhao, zero_weight_marginal

  area_xy = (app%sim%box_max(1) - app%sim%box_min(1))* &
            (app%sim%box_max(2) - app%sim%box_min(2))
  if (.not. ieee_is_finite(area_xy) .or. area_xy <= 0.0_dp .or. &
      .not. ieee_is_finite(app%sim%batch_duration) .or. app%sim%batch_duration <= 0.0_dp) then
    error stop 'implicit mean coupling requires finite positive area and batch duration.'
  end if
  if (size(workspace%mean_pending_charge) /= mesh%nelem .or. &
      size(workspace%mean_deferred_source_charge) /= mesh%nelem .or. &
      size(workspace%mean_returned_destination_charge) /= mesh%nelem .or. &
      size(workspace%mean_candidate_charge) /= mesh%nelem) then
    error stop 'implicit mean coupling received inconsistent element storage.'
  end if

  collision_failure_status = collision_query_ok
  collision_failure_particle = huge(0_i32)
  collision_failure_step = 0_i32
  collision_failure_x = 0.0_dp
  collision_failure_v = 0.0_dp
  transaction_residual = 0.0_dp
  sampled_deferred_absorbed_charge = 0.0_dp
  sampled_deferred_escaped_charge = 0.0_dp
  photoelectron_outward_current_density = 0.0_dp
  sample_escape_fraction = 0.0_dp
  return_weight_scale = 0.0_dp
  message = ''
  failure_message = ''

  call measure_tracked_ambient_surface_currents( &
    app, pcls_batch, workspace%absorbed_flag, area_xy, app%sim%batch_duration, mpi, &
    electron_current, ion_current &
    )
  call measure_tracked_photoelectron_interface_current( &
    app, workspace%interface_outward_thread, area_xy, app%sim%batch_duration, mpi, &
    photoelectron_outward_current_density, photoelectron_index &
    )
  call measure_pending_surface_current_density( &
    workspace%dq_thread, workspace%photo_emission_dq, area_xy, app%sim%batch_duration, mpi, pending_current &
    )
  additional_current = pending_current - electron_current - ion_current - &
                       photoelectron_outward_current_density

  ! pending は通常の局所 deposit と全放出 countercharge。scalar solve が
  ! 決めた return 総量を発生元で一時中和し、同じ総電荷を持つ平均場を作る。
  workspace%mean_pending_charge = sum(workspace%dq_thread, dim=2) + workspace%photo_emission_dq
  call mpi_allreduce_sum_real_dp_array(mpi, workspace%mean_pending_charge)
  source_charge_total = photoelectron_outward_current_density*area_xy*app%sim%batch_duration
  pending_absolute_total = sum(mesh%q_elem) + sum(workspace%mean_pending_charge)
  use_dynamic_zhao = trim(lower_ascii(snapshot%kinetic_options%kinetic_closure)) == 'zhao_charge_driven'
  if (use_dynamic_zhao) then
    local_energy_count = count(workspace%deferred_mean_interface_flag(:pcls_batch%n))
    allocate (local_energy_j(local_energy_count), local_energy_charge_c(local_energy_count))
    energy_index = 0_i32
    distribution_status = dynamic_k0_ok
    do i = 1_i32, pcls_batch%n
      if (.not. workspace%deferred_mean_interface_flag(i)) cycle
      energy_index = energy_index + 1_i32
      macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
      local_energy_j(energy_index) = &
        0.5_dp*pcls_batch%m(i)*workspace%deferred_mean_interface_crossing(i)%velocity(3)**2
      local_energy_charge_c(energy_index) = -macro_charge
      if (.not. ieee_is_finite(local_energy_j(energy_index)) .or. &
          local_energy_j(energy_index) < 0.0_dp .or. &
          .not. ieee_is_finite(local_energy_charge_c(energy_index)) .or. &
          local_energy_charge_c(energy_index) <= 0.0_dp) then
        distribution_status = dynamic_k0_invalid
      end if
    end do
    failure_count = merge(1_i32, 0_i32, distribution_status /= dynamic_k0_ok)
    call mpi_allreduce_sum_i32_scalar(mpi, failure_count)
    if (failure_count > 0_i32) then
      error stop 'implicit Zhao mean collected an invalid interface energy sample.'
    end if
    call mpi_gatherv_real_dp_array(mpi, local_energy_j, global_energy_j, 0_i32)
    call mpi_gatherv_real_dp_array(mpi, local_energy_charge_c, global_energy_charge_c, 0_i32)
    scalar_predictor = dynamic_k0_step_type()
    effective_source_scale = 0.0_dp
    if (mpi_is_root(mpi)) then
      call build_measured_interface_energy_distribution( &
        global_energy_j, global_energy_charge_c, measured_distribution, distribution_status, message &
        )
      if (distribution_status == dynamic_k0_ok) then
        surface_charge_base = pending_absolute_total - measured_distribution%total_charge_c
        call advance_dynamic_k0_zhao( &
          snapshot%kinetic_options, snapshot%outer, app%periodic2%lower_boundary_model, area_xy, &
          sum(mesh%q_elem), surface_charge_base, app%sim%batch_duration, measured_distribution, &
          scalar_predictor, accepted_outer, effective_source_scale, message &
          )
      else
        scalar_predictor%status = distribution_status
      end if
    end if
    call broadcast_dynamic_zhao_step(mpi, scalar_predictor, effective_source_scale)
    if (scalar_predictor%status /= dynamic_k0_ok) then
      if (mpi_is_root(mpi)) then
        write (error_unit, '(2a)') 'implicit Zhao mean: ', trim(message)
        write (error_unit, '(a,i0,4(a,es13.5))') &
          'implicit Zhao measured-source failure: interface_sample_count=', &
          size(global_energy_charge_c), &
          ' gathered_source_charge_C=', sum(global_energy_charge_c), &
          ' reduced_source_charge_C=', source_charge_total, &
          ' requested_source_scale=', effective_source_scale, &
          ' committed_source_scale=', snapshot%outer%photoelectron_source_scale
      end if
      step = scalar_predictor
      failure_message = message
      return
    end if
    if (abs(source_charge_total - scalar_predictor%photoelectron_source_charge_c) > &
        4096.0_dp*epsilon(1.0_dp)*max( &
        abs(source_charge_total), abs(scalar_predictor%photoelectron_source_charge_c), tiny(1.0_dp) &
        )) then
      error stop 'implicit Zhao gathered source charge disagrees with the interface-current reduction.'
    end if
    source_charge_total = scalar_predictor%photoelectron_source_charge_c
    snapshot%kinetic_options%photoelectron_source_scale = effective_source_scale
    snapshot%kinetic_options%photoelectron_population_fraction = 1.0_dp
    snapshot%kinetic_options%photoelectron_column_closure_enabled = .false.
    snapshot%kinetic_options%zhao_branch = lower_ascii(scalar_predictor%zhao_branch)
  else
    call advance_dynamic_k0_mean( &
      snapshot%kinetic_options, app%periodic2%lower_boundary_model, area_xy, &
      sum(mesh%q_elem), app%sim%batch_duration, scalar_predictor, message, &
      electron_current, ion_current, photoelectron_outward_current_density, additional_current &
      )
    if (scalar_predictor%status /= dynamic_k0_ok) then
      step = scalar_predictor
      failure_message = message
      return
    end if
  end if
  step = scalar_predictor
  target_return_charge = source_charge_total*scalar_predictor%photoelectron_return_fraction
  mean_charge_scale = max( &
                      abs(sum(mesh%q_elem)), abs(pending_absolute_total), &
                      abs(scalar_predictor%surface_charge_after_c), abs(source_charge_total), tiny(1.0_dp) &
                      )
  boundary_factor = 1.0_dp
  if (trim(lower_ascii(app%periodic2%lower_boundary_model)) == 'symmetric_vacuum') boundary_factor = 2.0_dp
  mean_capacitance = boundary_factor*eps0*area_xy/snapshot%kinetic_options%tail_length
  mean_charge_tolerance = max( &
                          4096.0_dp*epsilon(1.0_dp)*mean_charge_scale, &
                          mean_capacitance*(1.0e-12_dp + 1.0e-12_dp*max( &
                                            1.0_dp, abs(scalar_predictor%interface_potential_before_v), &
                                            abs(scalar_predictor%interface_potential_after_v) &
                                            )) &
                          )
  sample_charge_tolerance = 4096.0_dp*epsilon(1.0_dp)*max(abs(source_charge_total), tiny(1.0_dp))
  vector_charge_scale = max( &
                        sum(abs(mesh%q_elem)), sum(abs(workspace%mean_pending_charge)), mean_charge_scale &
                        )
  transaction_tolerance = max( &
                          mean_charge_tolerance, 4096.0_dp*epsilon(1.0_dp)*vector_charge_scale &
                          )
  if (target_return_charge < -mean_charge_tolerance .or. &
      target_return_charge > source_charge_total + mean_charge_tolerance) then
    write (message, '(a,3(a,es13.5))') &
      'implicit mean scalar target is outside the tracked photoelectron source interval:', &
      ' target_return_C=', target_return_charge, ' source_C=', source_charge_total, &
      ' tolerance_C=', mean_charge_tolerance
    error stop trim(message)
  end if
  target_return_charge = min(source_charge_total, max(0.0_dp, target_return_charge))
  if (source_charge_total > mean_charge_tolerance) then
    if (target_return_charge <= mean_charge_tolerance) target_return_charge = 0.0_dp
    if (source_charge_total - target_return_charge <= mean_charge_tolerance) then
      target_return_charge = source_charge_total
    end if
  end if
  target_escape_charge = source_charge_total - target_return_charge
  return_fraction = 0.0_dp
  if (source_charge_total > sample_charge_tolerance) return_fraction = target_return_charge/source_charge_total

  allocate (sample_return_weight(pcls_batch%n))
  sample_return_weight = 0.0_dp
  workspace%mean_deferred_source_charge = 0.0_dp
  transaction_status = mean_charge_transaction_ok
  do i = 1_i32, pcls_batch%n
    if (.not. workspace%deferred_mean_interface_flag(i)) cycle
    elem_idx = pcls_batch%source_element(i)
    macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
    if (elem_idx < 1_i32 .or. elem_idx > mesh%nelem .or. &
        .not. ieee_is_finite(macro_charge) .or. macro_charge >= 0.0_dp) then
      transaction_status = mean_charge_transaction_invalid
      cycle
    end if
    sample_return_fraction = return_fraction
    if (use_dynamic_zhao) then
      sample_return_fraction = 1.0_dp - measured_sample_escape_fraction( &
                               0.5_dp*pcls_batch%m(i)* &
                               workspace%deferred_mean_interface_crossing(i)%velocity(3)**2, &
                               scalar_predictor &
                               )
    end if
    sample_return_weight(i) = sample_return_fraction
    workspace%mean_deferred_source_charge(elem_idx) = &
      workspace%mean_deferred_source_charge(elem_idx) - macro_charge*sample_return_fraction
  end do
  failure_count = merge(1_i32, 0_i32, transaction_status /= mean_charge_transaction_ok)
  call mpi_allreduce_sum_i32_scalar(mpi, failure_count)
  if (failure_count > 0_i32) then
    error stop 'implicit mean deferred photoelectron has invalid source provenance or charge.'
  end if
  call mpi_allreduce_sum_real_dp_array(mpi, workspace%mean_deferred_source_charge)
  if (use_dynamic_zhao .and. &
      abs(sum(workspace%mean_deferred_source_charge) - target_return_charge) > transaction_tolerance) then
    error stop 'implicit Zhao measured-energy source weights do not reproduce the nonlinear return charge.'
  end if
  workspace%mean_returned_destination_charge = -workspace%mean_deferred_source_charge
  call build_mean_charge_transaction( &
    workspace%mean_pending_charge, workspace%mean_deferred_source_charge, &
    workspace%mean_returned_destination_charge, transaction, transaction_status, message &
    )
  if (transaction_status /= mean_charge_transaction_ok) then
    error stop 'implicit mean predictor transaction: '//trim(message)
  end if
  workspace%mean_candidate_charge = mesh%q_elem + transaction%predictor_charge_c
  transaction_residual = max( &
                         abs(sum(workspace%mean_candidate_charge) - scalar_predictor%surface_charge_after_c), &
                         abs(transaction%source_return_residual_charge_c), &
                         abs(transaction%predictor_actual_residual_charge_c) &
                         )
  if (transaction_residual > max(transaction_tolerance, transaction%charge_tolerance_c)) then
    error stop 'implicit mean predictor does not reproduce the scalar target charge.'
  end if

  nth = size(workspace%dq_thread, 2)
  allocate (outward_event_count(pcls_batch%n), returned_event_count(pcls_batch%n))
  allocate (returned_flight_time_sum(pcls_batch%n))
  allocate (mean_tau_max_thread(nth), mean_frozen_ratio_max_thread(nth), mean_energy_error_max_thread(nth))
  if (use_dynamic_zhao) then
    call snapshot%adopt_mean_outer(mesh, workspace%mean_candidate_charge, accepted_outer)
    resolved_barrier = zhao_profile_barrier_energy( &
                       snapshot%outer, snapshot%kinetic_options%photoelectron_charge &
                       )
    if (snapshot%outer%zhao_branch /= scalar_predictor%zhao_branch .or. &
        .not. ieee_is_finite(resolved_barrier) .or. &
        abs(snapshot%outer%photoelectron_source_scale - scalar_predictor%zhao_effective_source_scale) > &
        1.0e-12_dp*max(1.0_dp, abs(scalar_predictor%zhao_effective_source_scale)) .or. &
        abs(snapshot%outer%interface_potential - scalar_predictor%interface_potential_after_v) > &
        1.0e-8_dp*max(1.0_dp, abs(scalar_predictor%interface_potential_after_v)) .or. &
        abs(resolved_barrier - scalar_predictor%photoelectron_barrier_energy_j) > &
        1.0e-8_dp*max( &
        snapshot%kinetic_options%photoelectron_temperature_j, &
        abs(scalar_predictor%photoelectron_barrier_energy_j), tiny(1.0_dp) &
        )) then
      error stop 'implicit Zhao collective refresh did not reproduce the accepted nonlinear branch state.'
    end if
  else
    call snapshot%refresh_mean_only(mesh, workspace%mean_candidate_charge)
  end if
  call process_deferred_mean_interface_particles( &
    mesh, app, boundary_contract, snapshot, pcls_batch, workspace%deferred_mean_interface_flag, &
    workspace%deferred_mean_interface_step, workspace%deferred_mean_interface_crossing, &
    workspace%deferred_mean_return_element, workspace%deferred_mean_terminal_absorbed, &
    workspace%deferred_mean_terminal_escaped, outward_event_count, returned_event_count, &
    returned_flight_time_sum, &
    mean_tau_max_thread, mean_frozen_ratio_max_thread, mean_energy_error_max_thread, bfield, &
    collision_failure_status, collision_failure_particle, collision_failure_step, &
    collision_failure_x, collision_failure_v &
    )
  failure_count = merge(1_i32, 0_i32, collision_failure_status /= collision_query_ok)
  call mpi_allreduce_sum_i32_scalar(mpi, failure_count)
  if (failure_count > 0_i32) return

  ! Zhaoでは非線形CDFが決めたray別return weightをそのまま局所配置へ使う。
  ! 再越境やterminal channel不一致は他rayへのweight付け替えで隠さない。
  recross_charge = 0.0_dp
  terminal_mismatch_charge = 0.0_dp
  if (use_dynamic_zhao) then
    do i = 1_i32, pcls_batch%n
      if (.not. workspace%deferred_mean_interface_flag(i)) cycle
      macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
      if (outward_event_count(i) /= 1_i32 .or. returned_event_count(i) > 1_i32) then
        recross_charge = recross_charge - macro_charge
      end if
      if (sample_return_weight(i) > 64.0_dp*epsilon(1.0_dp)) then
        if (.not. workspace%deferred_mean_terminal_absorbed(i) .or. &
            returned_event_count(i) /= 1_i32) then
          terminal_mismatch_charge = terminal_mismatch_charge - macro_charge
        end if
      else
        zero_weight_marginal = &
          scalar_predictor%marginal_photoelectron_escape_fraction >= 0.0_dp .and. &
          scalar_predictor%marginal_photoelectron_escape_fraction >= &
          1.0_dp - 64.0_dp*epsilon(1.0_dp) .and. &
          0.5_dp*pcls_batch%m(i)* &
          workspace%deferred_mean_interface_crossing(i)%velocity(3)**2 == &
          scalar_predictor%marginal_photoelectron_energy_j
        if (.not. zero_weight_marginal) then
          if (.not. workspace%deferred_mean_terminal_escaped(i) .or. &
              returned_event_count(i) /= 0_i32) then
            terminal_mismatch_charge = terminal_mismatch_charge - macro_charge
          end if
        end if
      end if
    end do
    call mpi_allreduce_sum_real_dp_scalar(mpi, recross_charge)
    call mpi_allreduce_sum_real_dp_scalar(mpi, terminal_mismatch_charge)
    scalar_predictor%photoelectron_recross_charge_fraction = recross_charge/source_charge_total
    scalar_predictor%photoelectron_terminal_mismatch_charge_fraction = &
      terminal_mismatch_charge/source_charge_total
    step%photoelectron_recross_charge_fraction = &
      scalar_predictor%photoelectron_recross_charge_fraction
    step%photoelectron_terminal_mismatch_charge_fraction = &
      scalar_predictor%photoelectron_terminal_mismatch_charge_fraction
    if (recross_charge > sample_charge_tolerance .or. &
        terminal_mismatch_charge > sample_charge_tolerance) then
      if (mpi_is_root(mpi)) then
        write (error_unit, '(a,2(a,es13.5))') &
          'implicit Zhao empirical-interface applicability failure:', &
          ' recross_charge_fraction=', scalar_predictor%photoelectron_recross_charge_fraction, &
          ' terminal_mismatch_charge_fraction=', &
          scalar_predictor%photoelectron_terminal_mismatch_charge_fraction
      end if
      error stop 'implicit Zhao empirical interface distribution is invalidated by recross or terminal mismatch.'
    end if
  end if

  workspace%mean_deferred_source_charge = 0.0_dp
  workspace%mean_returned_destination_charge = 0.0_dp
  raw_return_charge = 0.0_dp
  empirical_return_charge = 0.0_dp
  transaction_status = mean_charge_transaction_ok
  do i = 1_i32, pcls_batch%n
    if (.not. workspace%deferred_mean_interface_flag(i)) cycle
    if (.not. workspace%deferred_mean_terminal_absorbed(i)) cycle
    elem_idx = pcls_batch%source_element(i)
    macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
    if (elem_idx < 1_i32 .or. elem_idx > mesh%nelem .or. &
        workspace%deferred_mean_return_element(i) < 1_i32 .or. &
        workspace%deferred_mean_return_element(i) > mesh%nelem .or. &
        .not. ieee_is_finite(macro_charge) .or. macro_charge >= 0.0_dp) then
      transaction_status = mean_charge_transaction_invalid
      cycle
    end if
    sample_return_fraction = 1.0_dp
    if (use_dynamic_zhao) sample_return_fraction = sample_return_weight(i)
    workspace%mean_deferred_source_charge(elem_idx) = &
      workspace%mean_deferred_source_charge(elem_idx) - macro_charge*sample_return_fraction
    elem_idx = workspace%deferred_mean_return_element(i)
    workspace%mean_returned_destination_charge(elem_idx) = &
      workspace%mean_returned_destination_charge(elem_idx) + macro_charge*sample_return_fraction
    raw_return_charge = raw_return_charge - macro_charge
    empirical_return_charge = empirical_return_charge - macro_charge*sample_return_fraction
  end do
  failure_count = merge(1_i32, 0_i32, transaction_status /= mean_charge_transaction_ok)
  call mpi_allreduce_sum_i32_scalar(mpi, failure_count)
  if (failure_count > 0_i32) then
    error stop 'implicit mean return sample contains invalid source or destination metadata.'
  end if
  call mpi_allreduce_sum_real_dp_array(mpi, workspace%mean_deferred_source_charge)
  call mpi_allreduce_sum_real_dp_array(mpi, workspace%mean_returned_destination_charge)
  call mpi_allreduce_sum_real_dp_scalar(mpi, raw_return_charge)
  call mpi_allreduce_sum_real_dp_scalar(mpi, empirical_return_charge)
  raw_escape_charge = max(0.0_dp, source_charge_total - raw_return_charge)
  sampled_deferred_absorbed_charge = -raw_return_charge
  sampled_deferred_escaped_charge = -raw_escape_charge
  if (source_charge_total > sample_charge_tolerance) then
    sample_escape_fraction = min(1.0_dp, raw_escape_charge/source_charge_total)
  else
    sample_escape_fraction = 1.0_dp
  end if
  if (use_dynamic_zhao) then
    return_weight_scale = 1.0_dp
    if (abs(empirical_return_charge - target_return_charge) > transaction_tolerance) then
      error stop 'implicit Zhao ray-resolved return weights do not reproduce the nonlinear target.'
    end if
  else if (target_return_charge > mean_charge_tolerance) then
    if (raw_return_charge <= sample_charge_tolerance) then
      error stop 'implicit mean needs at least one returned ray; increase photoelectron rays_per_batch.'
    end if
    return_weight_scale = target_return_charge/raw_return_charge
  else
    return_weight_scale = 0.0_dp
  end if
  workspace%mean_deferred_source_charge = &
    workspace%mean_deferred_source_charge*return_weight_scale
  workspace%mean_returned_destination_charge = &
    workspace%mean_returned_destination_charge*return_weight_scale
  call build_mean_charge_transaction( &
    workspace%mean_pending_charge, workspace%mean_deferred_source_charge, &
    workspace%mean_returned_destination_charge, transaction, transaction_status, message &
    )
  if (transaction_status /= mean_charge_transaction_ok) then
    error stop 'implicit mean sampled return transaction: '//trim(message)
  end if
  workspace%mean_candidate_charge = mesh%q_elem + transaction%actual_charge_c
  transaction_residual = max( &
                         abs(sum(workspace%mean_candidate_charge) - scalar_predictor%surface_charge_after_c), &
                         abs(transaction%source_return_residual_charge_c), &
                         abs(transaction%predictor_actual_residual_charge_c) &
                         )
  vector_charge_scale = max(vector_charge_scale, sum(abs(workspace%mean_candidate_charge)))
  transaction_tolerance = max( &
                          transaction_tolerance, 4096.0_dp*epsilon(1.0_dp)*vector_charge_scale &
                          )
  if (transaction_residual > max(transaction_tolerance, transaction%charge_tolerance_c)) then
    error stop 'implicit mean sampled return distribution does not reproduce the scalar target charge.'
  end if

  escape_weight_scale = 0.0_dp
  if (.not. use_dynamic_zhao .and. raw_escape_charge > sample_charge_tolerance) then
    escape_weight_scale = target_escape_charge/raw_escape_charge
  end if
  workspace%interface_outward_thread(photoelectron_index, :) = 0.0_dp
  workspace%interface_returned_thread(photoelectron_index, :) = 0.0_dp
  returning_shadow_moments = 0.0_dp
  do i = 1_i32, pcls_batch%n
    if (.not. workspace%deferred_mean_interface_flag(i)) cycle
    macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
    if (workspace%deferred_mean_terminal_absorbed(i)) then
      particle_weight_scale = return_weight_scale
      if (use_dynamic_zhao) particle_weight_scale = sample_return_weight(i)
      elem_idx = workspace%deferred_mean_return_element(i)
      workspace%dq_thread(elem_idx, 1) = workspace%dq_thread(elem_idx, 1) + &
                                         macro_charge*particle_weight_scale
    else
      particle_weight_scale = escape_weight_scale
    end if
    workspace%absorbed_flag(i) = workspace%deferred_mean_terminal_absorbed(i)
    workspace%escaped_boundary_flag(i) = workspace%deferred_mean_terminal_escaped(i)
    if (use_dynamic_zhao) then
      workspace%interface_outward_thread(photoelectron_index, 1) = &
        workspace%interface_outward_thread(photoelectron_index, 1) + macro_charge
      workspace%interface_returned_thread(photoelectron_index, 1) = &
        workspace%interface_returned_thread(photoelectron_index, 1) + &
        macro_charge*sample_return_weight(i)
      returning_shadow_moments(1) = returning_shadow_moments(1) - &
                                    macro_charge*sample_return_weight(i)
      returning_shadow_moments(2) = returning_shadow_moments(2) - &
                                    macro_charge*sample_return_weight(i)*returned_flight_time_sum(i)
    else
      workspace%interface_outward_thread(photoelectron_index, 1) = &
        workspace%interface_outward_thread(photoelectron_index, 1) + &
        macro_charge*particle_weight_scale*real(outward_event_count(i), dp)
      workspace%interface_returned_thread(photoelectron_index, 1) = &
        workspace%interface_returned_thread(photoelectron_index, 1) + &
        macro_charge*particle_weight_scale*real(returned_event_count(i), dp)
      returning_shadow_moments(1) = returning_shadow_moments(1) - &
                                    macro_charge*particle_weight_scale*real(returned_event_count(i), dp)
      returning_shadow_moments(2) = returning_shadow_moments(2) - &
                                    macro_charge*particle_weight_scale*returned_flight_time_sum(i)
    end if
  end do
  call mpi_allreduce_sum_real_dp_array(mpi, returning_shadow_moments)
  if (.not. all(ieee_is_finite(returning_shadow_moments)) .or. any(returning_shadow_moments < 0.0_dp)) then
    error stop 'implicit mean returning shadow diagnostics must be finite and nonnegative.'
  end if
  if (returning_shadow_moments(1) > sample_charge_tolerance) then
    step%returned_outer_flight_time_mean_s = returning_shadow_moments(2)/returning_shadow_moments(1)
    step%estimated_returning_photoelectron_column_charge_per_area_c_m2 = &
      returning_shadow_moments(2)/(area_xy*app%sim%batch_duration)
  else
    step%returned_outer_flight_time_mean_s = 0.0_dp
    step%estimated_returning_photoelectron_column_charge_per_area_c_m2 = 0.0_dp
  end if
  if (.not. all(ieee_is_finite([ &
                               step%returned_outer_flight_time_mean_s, &
                               step%estimated_returning_photoelectron_column_charge_per_area_c_m2 &
                               ])) .or. step%returned_outer_flight_time_mean_s < 0.0_dp .or. &
      step%estimated_returning_photoelectron_column_charge_per_area_c_m2 < 0.0_dp) then
    error stop 'implicit mean returning shadow estimates must be finite and nonnegative.'
  end if
  if (.not. use_dynamic_zhao .and. raw_escape_charge <= sample_charge_tolerance .and. &
      target_escape_charge > mean_charge_tolerance .and. &
      mpi%rank == 0_i32) then
    ! No sampled escape path exists.  Its only required gross event is the
    ! first outward crossing; it has no surface destination to distribute.
    workspace%interface_outward_thread(photoelectron_index, 1) = &
      workspace%interface_outward_thread(photoelectron_index, 1) - target_escape_charge
  end if
  workspace%interface_tau_max_thread = max(workspace%interface_tau_max_thread, mean_tau_max_thread)
  workspace%interface_frozen_ratio_max_thread = &
    max(workspace%interface_frozen_ratio_max_thread, mean_frozen_ratio_max_thread)
  workspace%interface_energy_error_max_thread = &
    max(workspace%interface_energy_error_max_thread, mean_energy_error_max_thread)

  step%surface_charge_after_c = sum(workspace%mean_candidate_charge)
  step%interface_potential_after_v = snapshot%outer%interface_potential
  step%interface_field_after_v_m = snapshot%outer%interface_field
  step%electron_current_density_a_m2 = electron_current
  step%ion_current_density_a_m2 = ion_current
  step%photoelectron_current_density_a_m2 = target_escape_charge/(area_xy*app%sim%batch_duration)
  step%total_current_density_a_m2 = &
    (step%surface_charge_after_c - step%surface_charge_before_c)/(area_xy*app%sim%batch_duration)
  step%additional_tracked_current_density_a_m2 = step%total_current_density_a_m2 - &
                                                 electron_current - ion_current - &
                                                 step%photoelectron_current_density_a_m2
  if (photoelectron_outward_current_density > 0.0_dp) then
    step%photoelectron_escape_fraction = target_escape_charge/source_charge_total
    step%photoelectron_return_fraction = 1.0_dp - step%photoelectron_escape_fraction
  else
    step%photoelectron_escape_fraction = 1.0_dp
    step%photoelectron_return_fraction = 0.0_dp
  end if
  step%status = dynamic_k0_ok
  end procedure resolve_implicit_mean_batch

  !> 保存済みinterface crossingを共通kinetic profile mapperへ渡し、終端まで1回追跡する。
  module procedure process_deferred_mean_interface_particles
  integer(i32) :: i, step_index, tid, nth, external_event
  integer(i32) :: local_failure_status, local_failure_step
  real(dp) :: position(3), velocity(3)
  real(dp) :: local_failure_x(3), local_failure_v(3)
  type(particle_step_result) :: step_result, external_final_result
  type(external_step_trace_type) :: external_trace
  logical :: terminal

  nth = size(interface_tau_max_thread)
  if (size(deferred_flag) < pcls_batch%n .or. size(deferred_step) < pcls_batch%n .or. &
      size(deferred_crossing) < pcls_batch%n .or. size(return_element) < pcls_batch%n .or. &
      size(terminal_absorbed) < pcls_batch%n .or. size(terminal_escaped) < pcls_batch%n .or. &
      size(outward_event_count) < pcls_batch%n .or. size(returned_event_count) < pcls_batch%n .or. &
      size(returned_flight_time_sum) < pcls_batch%n .or. &
      size(interface_tau_max_thread) /= nth .or. size(interface_frozen_ratio_max_thread) /= nth .or. &
      size(interface_energy_error_max_thread) /= nth) then
    error stop 'implicit mean interface transport received inconsistent batch storage.'
  end if

  return_element(:pcls_batch%n) = -1_i32
  terminal_absorbed(:pcls_batch%n) = .false.
  terminal_escaped(:pcls_batch%n) = .false.
  outward_event_count(:pcls_batch%n) = 0_i32
  returned_event_count(:pcls_batch%n) = 0_i32
  returned_flight_time_sum(:pcls_batch%n) = 0.0_dp
  interface_tau_max_thread = 0.0_dp
  interface_frozen_ratio_max_thread = 0.0_dp
  interface_energy_error_max_thread = 0.0_dp
  collision_failure_status = collision_query_ok
  collision_failure_particle = huge(0_i32)
  collision_failure_step = 0_i32
  collision_failure_x = 0.0_dp
  collision_failure_v = 0.0_dp

  !$omp parallel default(none) &
  !$omp shared(mesh,app,boundary_contract,snapshot,pcls_batch,deferred_flag,deferred_step,deferred_crossing) &
  !$omp shared(return_element,terminal_absorbed,terminal_escaped,bfield,nth) &
  !$omp shared(outward_event_count,returned_event_count,returned_flight_time_sum) &
  !$omp shared(interface_tau_max_thread,interface_frozen_ratio_max_thread,interface_energy_error_max_thread) &
  !$omp shared(collision_failure_status,collision_failure_particle,collision_failure_step) &
  !$omp shared(collision_failure_x,collision_failure_v) &
  !$omp private(i,step_index,tid,external_event,position,velocity) &
  !$omp private(local_failure_status,local_failure_step,local_failure_x,local_failure_v) &
  !$omp private(step_result,external_final_result,external_trace,terminal)
  tid = 1_i32
!$ tid = omp_get_thread_num() + 1_i32
  !$omp do schedule(runtime)
  do i = 1_i32, pcls_batch%n
    if (.not. deferred_flag(i)) cycle

    local_failure_status = collision_query_ok
    local_failure_step = deferred_step(i)
    local_failure_x = deferred_crossing(i)%position
    local_failure_v = deferred_crossing(i)%velocity
    terminal = .false.
    if (pcls_batch%species_id(i) < 1_i32 .or. pcls_batch%species_id(i) > app%n_particle_species .or. &
        deferred_step(i) < 1_i32 .or. deferred_step(i) > app%sim%max_step .or. &
        .not. ieee_is_finite(deferred_crossing(i)%dt_remaining) .or. &
        deferred_crossing(i)%dt_remaining < 0.0_dp .or. &
        .not. ieee_is_finite(pcls_batch%q(i)*pcls_batch%w(i)) .or. &
        pcls_batch%q(i)*pcls_batch%w(i) >= 0.0_dp) then
      local_failure_status = particle_step_invalid_boundary
    else
      step_result = particle_step_result()
      step_result%x = deferred_crossing(i)%position
      step_result%v = deferred_crossing(i)%velocity
      step_result%interface_crossing = deferred_crossing(i)
      call continue_external_particle_step( &
        boundary_contract, snapshot, mesh, app%sim, app%coupling, bfield, &
        pcls_batch%q(i), pcls_batch%m(i), app%sim%batch_duration, step_result, external_final_result, &
        external_trace, enforce_frozen_field_limit=.false. &
        )
      step_result = external_final_result

      do
        do external_event = 1_i32, external_trace%count
          outward_event_count(i) = outward_event_count(i) + 1_i32
          if (external_trace%outcome(external_event)%kind == interface_outcome_returned_local) then
            returned_event_count(i) = returned_event_count(i) + 1_i32
            returned_flight_time_sum(i) = returned_flight_time_sum(i) + &
                                          external_trace%outcome(external_event)%outer_flight_time
          end if
          if (external_trace%outcome(external_event)%kind == interface_outcome_returned_local .or. &
              external_trace%outcome(external_event)%kind == interface_outcome_escaped_to_infinity .or. &
              external_trace%outcome(external_event)%kind == interface_outcome_queued_outer) then
            interface_tau_max_thread(tid) = max( &
                                            interface_tau_max_thread(tid), &
                                            external_trace%outcome(external_event)%outer_flight_time &
                                            )
            interface_frozen_ratio_max_thread(tid) = max( &
                                                     interface_frozen_ratio_max_thread(tid), &
                                                     external_trace%outcome(external_event)%frozen_field_ratio &
                                                     )
            interface_energy_error_max_thread(tid) = max( &
                                                     interface_energy_error_max_thread(tid), &
                                                     external_trace%outcome(external_event)%energy_relative_error &
                                                     )
          end if
        end do
        if (step_result%status /= collision_query_ok) then
          local_failure_status = step_result%status
          local_failure_x = step_result%x
          local_failure_v = step_result%v
          exit
        end if
        if (step_result%absorbed) then
          if (step_result%elem_idx < 1_i32 .or. step_result%elem_idx > mesh%nelem) then
            local_failure_status = particle_step_invalid_boundary
            local_failure_x = step_result%x
            local_failure_v = step_result%v
          else
            return_element(i) = step_result%elem_idx
            terminal_absorbed(i) = .true.
            terminal = .true.
          end if
          exit
        end if
        if (step_result%escaped_boundary) then
          if (external_trace_ends_at_infinity_escape(external_trace)) then
            terminal_escaped(i) = .true.
            terminal = .true.
          else
            ! A return followed by an ordinary local open-face escape is not
            ! part of the analytic z-high escape channel used by this closure.
            local_failure_status = particle_step_invalid_external_model
            local_failure_x = step_result%x
            local_failure_v = step_result%v
          end if
          exit
        end if

        position = step_result%x
        velocity = step_result%v
        step_index = local_failure_step + 1_i32
        if (step_index > app%sim%max_step) exit
        call advance_particle_step( &
          mesh, app%sim, snapshot, bfield, position, velocity, &
          pcls_batch%q(i), pcls_batch%m(i), app%sim%dt, step_result, &
          boundary_contract=boundary_contract &
          )
        local_failure_step = step_index
        external_trace = external_step_trace_type()
        if (step_result%interface_crossing%has_crossing) then
          call continue_external_particle_step( &
            boundary_contract, snapshot, mesh, app%sim, app%coupling, bfield, &
            pcls_batch%q(i), pcls_batch%m(i), app%sim%batch_duration, step_result, external_final_result, &
            external_trace, enforce_frozen_field_limit=.false. &
            )
          step_result = external_final_result
        end if
      end do
      if (.not. terminal .and. local_failure_status == collision_query_ok) then
        local_failure_status = particle_step_invalid_boundary
        local_failure_step = app%sim%max_step
        local_failure_x = step_result%x
        local_failure_v = step_result%v
      end if
    end if

    if (local_failure_status /= collision_query_ok) then
      !$omp critical (beach_collision_query_failure)
      if (collision_failure_status == collision_query_ok .or. i < collision_failure_particle .or. &
          (i == collision_failure_particle .and. local_failure_step < collision_failure_step)) then
        collision_failure_status = local_failure_status
        collision_failure_particle = i
        collision_failure_step = local_failure_step
        collision_failure_x = local_failure_x
        collision_failure_v = local_failure_v
      end if
      !$omp end critical (beach_collision_query_failure)
    end if
  end do
  !$omp end do
  !$omp end parallel
  end procedure process_deferred_mean_interface_particles

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
    case (particle_step_invalid_external_model)
      failure_name = 'invalid_external_model'
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
  if (workspace%charge_candidate_ready) then
    mesh%q_elem = workspace%mean_candidate_charge
  else
    workspace%dq = sum(workspace%dq_thread, dim=2) + workspace%photo_emission_dq
    call mpi_allreduce_sum_real_dp_array(mpi, workspace%dq)
    mesh%q_elem = mesh%q_elem + workspace%dq
  end if
  call apply_surface_model_charge_relaxation(mesh, external_e, field_bc_mode=field_bc_mode)
  workspace%dq = mesh%q_elem - workspace%q_before
  norm_dq = sqrt(sum(workspace%dq*workspace%dq))
  norm_q = sqrt(sum(mesh%q_elem*mesh%q_elem))
  rel = norm_dq/max(norm_q, q_floor)
  workspace%charge_candidate_ready = .false.
  end procedure commit_batch_charge

  !> 適応判定と commit が同じ候補電荷を共有できるよう、全rankの差分を確定する。
  subroutine prepare_adaptive_charge_candidate(mesh, workspace, implicit_mean_enabled, mpi)
    type(mesh_type), intent(in) :: mesh
    type(simulator_batch_workspace_type), intent(inout) :: workspace
    logical, intent(in) :: implicit_mean_enabled
    type(mpi_context), intent(in) :: mpi

    workspace%q_before = mesh%q_elem
    if (implicit_mean_enabled) then
      if (size(workspace%mean_candidate_charge) /= mesh%nelem) then
        error stop 'adaptive nonzero-mode candidate storage does not match the mesh.'
      end if
    else
      workspace%dq = sum(workspace%dq_thread, dim=2) + workspace%photo_emission_dq
      call mpi_allreduce_sum_real_dp_array(mpi, workspace%dq)
      workspace%mean_candidate_charge = mesh%q_elem + workspace%dq
    end if
    if (.not. all(ieee_is_finite(workspace%mean_candidate_charge))) then
      error stop 'adaptive nonzero-mode candidate charge is not finite.'
    end if
    workspace%charge_candidate_ready = .true.
  end subroutine prepare_adaptive_charge_candidate

  !> unresolved な全反射 photoelectron を、resolved return の一様な重み補正へ閉じ込める。
  !!
  !! `neutral_return` は正味 photoelectron 表面電荷増分をゼロにする閉包であり、
  !! 実際の吸収位置が与える空間再分配を保持する。異なる高さの面へ移る電荷が作る
  !! plane-averaged dipole までは除去しない。max-step survivor の電荷を resolved
  !! absorption に species ごとに比例配分し、escape/outer transfer/soft discard は
  !! 物理的に別の終端なので補正せず停止する。
  subroutine apply_neutral_return_surface_closure(app, pcls_batch, fresh_particle_count, workspace, mpi)
    type(app_config), intent(in) :: app
    type(particles_soa), intent(in) :: pcls_batch
    integer(i32), intent(in) :: fresh_particle_count
    type(simulator_batch_workspace_type), intent(inout) :: workspace
    type(mpi_context), intent(in) :: mpi
    integer(i32) :: i, species_idx, elem_idx, n, terminal_count
    integer(i64) :: escaped_count, soft_count, invalid_count
    real(dp) :: macro_charge, emitted_charge, absorbed_charge, unresolved_charge
    real(dp) :: escaped_charge, soft_charge, invalid_charge
    real(dp) :: charge_scale, charge_tolerance, balance_tolerance
    real(dp) :: weight_scale, correction_charge, unresolved_fraction
    logical :: has_neutral_return
    character(len=512) :: message

    if (fresh_particle_count < 0_i32 .or. fresh_particle_count > pcls_batch%n) then
      error stop 'neutral_return received an invalid fresh particle count.'
    end if
    n = app%n_particle_species
    if (n < 1_i32 .or. &
        size(workspace%neutral_return_charge_values) /= 6*n .or. &
        size(workspace%neutral_return_terminal_counts) /= 3*n .or. &
        size(workspace%neutral_return_emitted_charge) /= n .or. &
        size(workspace%neutral_return_absorbed_charge) /= n .or. &
        size(workspace%neutral_return_unresolved_charge) /= n .or. &
        size(workspace%neutral_return_weight_scale) /= n .or. &
        size(workspace%neutral_return_correction) /= n .or. &
        size(workspace%neutral_return_unresolved_fraction) /= n) then
      error stop 'neutral_return workspace dimensions are inconsistent.'
    end if

    has_neutral_return = .false.
    do species_idx = 1_i32, n
      if (.not. app%particle_species(species_idx)%enabled) cycle
      if (trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) == &
          'neutral_return') then
        has_neutral_return = .true.
        if (trim(lower_ascii(app%particle_species(species_idx)%source_mode)) /= 'photo_raycast' .or. &
            .not. ieee_is_finite(app%particle_species(species_idx)%q_particle) .or. &
            app%particle_species(species_idx)%q_particle >= 0.0_dp .or. &
            .not. app%particle_species(species_idx)%deposit_opposite_charge_on_emit .or. &
            trim(lower_ascii(app%particle_species(species_idx)%z_high_boundary)) /= 'reflect') then
          error stop 'neutral_return runtime requires a negative reflected photo_raycast with countercharge deposit.'
        end if
      end if
    end do
    if (.not. has_neutral_return) return
    if (.not. app%sim%use_box .or. app%sim%bc_high(3) /= bc_open) then
      error stop 'neutral_return runtime requires a finite box with global z-high open.'
    end if
    if (trim(lower_ascii(app%coupling%particle_transfer_mode)) /= 'none' .or. &
        app%coupling%outer_queue_enabled .or. &
        trim(lower_ascii(app%coupling%update_mode)) == 'implicit_mean') then
      error stop 'neutral_return runtime cannot be combined with outer particle transfer or implicit_mean.'
    end if

    workspace%neutral_return_charge_values = 0.0_dp
    workspace%neutral_return_terminal_counts = 0_i64
    do i = 1_i32, pcls_batch%n
      species_idx = pcls_batch%species_id(i)
      if (species_idx < 1_i32 .or. species_idx > n) then
        error stop 'neutral_return encountered an invalid particle species index.'
      end if
      if (trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) /= &
          'neutral_return') cycle

      macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
      if (.not. ieee_is_finite(macro_charge) .or. macro_charge >= 0.0_dp) then
        workspace%neutral_return_terminal_counts(2*n + species_idx) = &
          workspace%neutral_return_terminal_counts(2*n + species_idx) + 1_i64
        cycle
      end if
      if (i > fresh_particle_count) then
        workspace%neutral_return_charge_values(5*n + species_idx) = &
          workspace%neutral_return_charge_values(5*n + species_idx) + macro_charge
        workspace%neutral_return_terminal_counts(2*n + species_idx) = &
          workspace%neutral_return_terminal_counts(2*n + species_idx) + 1_i64
        cycle
      end if
      workspace%neutral_return_charge_values(species_idx) = &
        workspace%neutral_return_charge_values(species_idx) + macro_charge

      terminal_count = merge(1_i32, 0_i32, workspace%absorbed_flag(i)) + &
                       merge(1_i32, 0_i32, workspace%escaped_boundary_flag(i)) + &
                       merge(1_i32, 0_i32, workspace%soft_discarded_boundary_flag(i)) + &
                       merge(1_i32, 0_i32, workspace%queued_outer_flag(i)) + &
                       merge(1_i32, 0_i32, pcls_batch%alive(i))
      if (terminal_count /= 1_i32) then
        workspace%neutral_return_charge_values(5*n + species_idx) = &
          workspace%neutral_return_charge_values(5*n + species_idx) + macro_charge
        workspace%neutral_return_terminal_counts(2*n + species_idx) = &
          workspace%neutral_return_terminal_counts(2*n + species_idx) + 1_i64
      else if (workspace%absorbed_flag(i)) then
        elem_idx = workspace%absorbed_element(i)
        if (elem_idx < 1_i32 .or. elem_idx > size(workspace%dq_thread, 1)) then
          workspace%neutral_return_charge_values(5*n + species_idx) = &
            workspace%neutral_return_charge_values(5*n + species_idx) + macro_charge
          workspace%neutral_return_terminal_counts(2*n + species_idx) = &
            workspace%neutral_return_terminal_counts(2*n + species_idx) + 1_i64
        else
          workspace%neutral_return_charge_values(n + species_idx) = &
            workspace%neutral_return_charge_values(n + species_idx) + macro_charge
        end if
      else if (workspace%escaped_boundary_flag(i)) then
        workspace%neutral_return_charge_values(3*n + species_idx) = &
          workspace%neutral_return_charge_values(3*n + species_idx) + macro_charge
        workspace%neutral_return_terminal_counts(species_idx) = &
          workspace%neutral_return_terminal_counts(species_idx) + 1_i64
      else if (workspace%soft_discarded_boundary_flag(i)) then
        workspace%neutral_return_charge_values(4*n + species_idx) = &
          workspace%neutral_return_charge_values(4*n + species_idx) + macro_charge
        workspace%neutral_return_terminal_counts(n + species_idx) = &
          workspace%neutral_return_terminal_counts(n + species_idx) + 1_i64
      else if (workspace%queued_outer_flag(i)) then
        workspace%neutral_return_charge_values(5*n + species_idx) = &
          workspace%neutral_return_charge_values(5*n + species_idx) + macro_charge
        workspace%neutral_return_terminal_counts(2*n + species_idx) = &
          workspace%neutral_return_terminal_counts(2*n + species_idx) + 1_i64
      else
        workspace%neutral_return_charge_values(2*n + species_idx) = &
          workspace%neutral_return_charge_values(2*n + species_idx) + macro_charge
      end if
    end do

    call mpi_allreduce_sum_real_dp_array(mpi, workspace%neutral_return_charge_values)
    call mpi_allreduce_sum_i64_array(mpi, workspace%neutral_return_terminal_counts)
    workspace%neutral_return_emitted_charge = workspace%neutral_return_charge_values(1:n)
    workspace%neutral_return_absorbed_charge = workspace%neutral_return_charge_values(n + 1:2*n)
    workspace%neutral_return_unresolved_charge = workspace%neutral_return_charge_values(2*n + 1:3*n)

    do species_idx = 1_i32, n
      if (.not. app%particle_species(species_idx)%enabled) cycle
      if (trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) /= &
          'neutral_return') cycle
      emitted_charge = workspace%neutral_return_emitted_charge(species_idx)
      absorbed_charge = workspace%neutral_return_absorbed_charge(species_idx)
      unresolved_charge = workspace%neutral_return_unresolved_charge(species_idx)
      escaped_charge = workspace%neutral_return_charge_values(3*n + species_idx)
      soft_charge = workspace%neutral_return_charge_values(4*n + species_idx)
      invalid_charge = workspace%neutral_return_charge_values(5*n + species_idx)
      escaped_count = workspace%neutral_return_terminal_counts(species_idx)
      soft_count = workspace%neutral_return_terminal_counts(n + species_idx)
      invalid_count = workspace%neutral_return_terminal_counts(2*n + species_idx)

      if (escaped_count > 0_i64 .or. soft_count > 0_i64 .or. invalid_count > 0_i64) then
        write (message, '(a,i0,3(a,i0),3(a,es13.5))') &
          'neutral_return has an unsupported terminal outcome for species ', species_idx, &
          ': escaped_count=', escaped_count, ' soft_count=', soft_count, &
          ' invalid_or_outer_count=', invalid_count, ' escaped_charge_C=', escaped_charge, &
          ' soft_charge_C=', soft_charge, ' invalid_or_outer_charge_C=', invalid_charge
        error stop trim(message)
      end if

      charge_scale = max( &
                     abs(emitted_charge), abs(absorbed_charge), abs(unresolved_charge), tiny(1.0_dp) &
                     )
      charge_tolerance = 4096.0_dp*epsilon(1.0_dp)*charge_scale
      balance_tolerance = max(charge_tolerance, sqrt(epsilon(1.0_dp))*charge_scale)
      if (abs(emitted_charge) <= charge_tolerance) then
        if (abs(absorbed_charge) > balance_tolerance .or. &
            abs(unresolved_charge) > balance_tolerance) then
          write (message, '(a,i0,3(a,es13.5))') &
            'neutral_return has terminal charge without emitted charge for species ', species_idx, &
            ': emitted_charge_C=', emitted_charge, ' absorbed_charge_C=', absorbed_charge, &
            ' unresolved_charge_C=', unresolved_charge
          error stop trim(message)
        end if
        cycle
      end if
      if (emitted_charge >= 0.0_dp .or. absorbed_charge >= -charge_tolerance .or. &
          unresolved_charge > charge_tolerance) then
        write (message, '(a,i0,3(a,es13.5))') &
          'neutral_return charge signs or resolved return are invalid for species ', species_idx, &
          ': emitted_charge_C=', emitted_charge, ' absorbed_charge_C=', absorbed_charge, &
          ' unresolved_charge_C=', unresolved_charge
        error stop trim(message)
      end if
      if (abs(emitted_charge - absorbed_charge - unresolved_charge) > balance_tolerance) then
        write (message, '(a,i0,4(a,es13.5))') &
          'neutral_return particle charge does not close for species ', species_idx, &
          ': emitted_charge_C=', emitted_charge, ' absorbed_charge_C=', absorbed_charge, &
          ' unresolved_charge_C=', unresolved_charge, ' tolerance_C=', balance_tolerance
        error stop trim(message)
      end if

      weight_scale = emitted_charge/absorbed_charge
      correction_charge = emitted_charge - absorbed_charge
      unresolved_fraction = unresolved_charge/emitted_charge
      if (.not. all(ieee_is_finite([weight_scale, correction_charge, unresolved_fraction])) .or. &
          weight_scale < 1.0_dp - sqrt(epsilon(1.0_dp)) .or. &
          unresolved_fraction < -sqrt(epsilon(1.0_dp)) .or. &
          unresolved_fraction > 1.0_dp + sqrt(epsilon(1.0_dp))) then
        write (message, '(a,i0,3(a,es13.5))') &
          'neutral_return derived an invalid correction for species ', species_idx, &
          ': weight_scale=', weight_scale, ' correction_C=', correction_charge, &
          ' unresolved_fraction=', unresolved_fraction
        error stop trim(message)
      end if
      if (unresolved_fraction > neutral_return_max_unresolved_fraction + sqrt(epsilon(1.0_dp))) then
        write (message, '(a,i0,2(a,es13.5))') &
          'neutral_return unresolved fraction exceeds the fixed 5 percent applicability limit for species ', &
          species_idx, ': unresolved_fraction=', unresolved_fraction, &
          ' weight_scale=', weight_scale
        error stop trim(message)
      end if
      workspace%neutral_return_weight_scale(species_idx) = max(1.0_dp, weight_scale)
      workspace%neutral_return_correction(species_idx) = correction_charge
      workspace%neutral_return_unresolved_fraction(species_idx) = &
        min(1.0_dp, max(0.0_dp, unresolved_fraction))
    end do

    ! Every rank applies the same MPI-global scale to its local resolved deposits.
    ! The later charge-vector allreduce therefore produces exactly the global correction.
    do i = 1_i32, fresh_particle_count
      species_idx = pcls_batch%species_id(i)
      if (trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) /= &
          'neutral_return') cycle
      if (.not. workspace%absorbed_flag(i)) cycle
      elem_idx = workspace%absorbed_element(i)
      macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
      workspace%dq_thread(elem_idx, 1) = workspace%dq_thread(elem_idx, 1) + &
                                         (workspace%neutral_return_weight_scale(species_idx) - 1.0_dp)*macro_charge
    end do
  end subroutine apply_neutral_return_surface_closure

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

  !> species 別の rank-local flux/count だけを MPI-global 値へ集約する。
  !!
  !! neutral_return の correction/scale/fraction は closure 内ですでに全rankへ
  !! allreduce 済みなので、この 7*nspecies scratch には含めない。
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

  !> implicit mean が決めた連続電荷量を photoelectron ledger の authoritative total にする。
  !!
  !! ray の absorbed/escaped count は軌道標本数のまま保持する。charge column は
  !! scalar closure の escape 量と、weighted gross crossing の恒等式に合わせる。
  subroutine reconcile_implicit_mean_charge_ledger( &
    app, step, area_xy, batch_duration, photoelectron_outward_current_density, &
    sampled_deferred_absorbed_charge, &
    sampled_deferred_escaped_charge, ledger &
    )
    type(app_config), intent(in) :: app
    type(dynamic_k0_step_type), intent(in) :: step
    real(dp), intent(in) :: area_xy, batch_duration
    real(dp), intent(in) :: photoelectron_outward_current_density
    real(dp), intent(in) :: sampled_deferred_absorbed_charge, sampled_deferred_escaped_charge
    type(charge_ledger_type), intent(inout) :: ledger
    integer(i32) :: species_idx, photoelectron_index, photoelectron_count
    real(dp) :: absorbed_charge, escaped_charge, source_charge
    real(dp) :: charge_scale, charge_tolerance, interface_residual

    photoelectron_index = 0_i32
    photoelectron_count = 0_i32
    do species_idx = 1_i32, app%n_particle_species
      if (.not. app%particle_species(species_idx)%enabled) cycle
      if (trim(lower_ascii(app%particle_species(species_idx)%source_mode)) /= 'photo_raycast') cycle
      if (app%particle_species(species_idx)%q_particle >= 0.0_dp) cycle
      photoelectron_count = photoelectron_count + 1_i32
      photoelectron_index = species_idx
    end do
    if (photoelectron_count /= 1_i32) then
      error stop 'implicit mean ledger requires exactly one negative photo_raycast species.'
    end if
    if (.not. all(ieee_is_finite([ &
                                 step%photoelectron_current_density_a_m2, area_xy, batch_duration, &
                                 photoelectron_outward_current_density, &
                                 sampled_deferred_absorbed_charge, sampled_deferred_escaped_charge &
                                 ])) .or. step%photoelectron_current_density_a_m2 < 0.0_dp .or. &
        area_xy <= 0.0_dp .or. batch_duration <= 0.0_dp .or. &
        photoelectron_outward_current_density < 0.0_dp .or. &
        sampled_deferred_absorbed_charge > 0.0_dp .or. sampled_deferred_escaped_charge > 0.0_dp) then
      error stop 'implicit mean ledger received invalid scales or photoelectron current.'
    end if

    source_charge = -photoelectron_outward_current_density*area_xy*batch_duration
    escaped_charge = -step%photoelectron_current_density_a_m2*area_xy*batch_duration
    absorbed_charge = source_charge - escaped_charge
    charge_scale = max( &
                   abs(source_charge), abs(absorbed_charge), abs(escaped_charge), &
                   abs(sampled_deferred_absorbed_charge), abs(sampled_deferred_escaped_charge), &
                   abs(ledger%interface_outward_gross(photoelectron_index)), &
                   abs(ledger%interface_returned_gross(photoelectron_index)), tiny(1.0_dp) &
                   )
    charge_tolerance = 4096.0_dp*epsilon(1.0_dp)*charge_scale
    interface_residual = ledger%interface_outward_gross(photoelectron_index) - &
                         ledger%interface_returned_gross(photoelectron_index) - escaped_charge
    if (abs(sampled_deferred_absorbed_charge + sampled_deferred_escaped_charge - source_charge) > &
        charge_tolerance .or. abs(absorbed_charge + escaped_charge - source_charge) > charge_tolerance .or. &
        abs(interface_residual) > charge_tolerance) then
      error stop 'implicit mean ledger interface charge does not close.'
    end if
    ledger%escaped_to_infinity(photoelectron_index) = &
      ledger%escaped_to_infinity(photoelectron_index) - sampled_deferred_escaped_charge + escaped_charge
    ledger%absorbed_on_surface(photoelectron_index) = &
      ledger%absorbed_on_surface(photoelectron_index) - sampled_deferred_absorbed_charge + absorbed_charge
  end subroutine reconcile_implicit_mean_charge_ledger

  !> First traceで実際にstageされた全surface charge deltaをcurrentへ変換する。
  !!
  !! このaggregateはspeciesやopen-face ownershipを仮定しない。implicit mean
  !! closureが後から置換するdeferred photoelectron return以外のk=0変化を、
  !! final projectionが黙って消さないためのauthoritativeな明示電流である。
  subroutine measure_pending_surface_current_density( &
    dq_thread, photo_emission_dq, area_xy, batch_duration, mpi, current_density &
    )
    real(dp), intent(in) :: dq_thread(:, :), photo_emission_dq(:)
    real(dp), intent(in) :: area_xy, batch_duration
    type(mpi_context), intent(in) :: mpi
    real(dp), intent(out) :: current_density

    if (size(dq_thread, 1) /= size(photo_emission_dq) .or. &
        .not. ieee_is_finite(area_xy) .or. .not. ieee_is_finite(batch_duration) .or. &
        area_xy <= 0.0_dp .or. batch_duration <= 0.0_dp) then
      error stop 'pending surface current measurement received invalid storage or scales.'
    end if
    current_density = sum(dq_thread) + sum(photo_emission_dq)
    call mpi_allreduce_sum_real_dp_scalar(mpi, current_density)
    current_density = current_density/(area_xy*batch_duration)
    if (.not. ieee_is_finite(current_density)) then
      error stop 'pending surface current measurement produced a non-finite value.'
    end if
  end subroutine measure_pending_surface_current_density

  !> tracked ambient speciesがこのbatchで実際に表面へ渡したsigned currentを測る。
  !!
  !! multirate更新でも3Dのcollection efficiencyと局所軌道は既存particle traceを
  !! source of truthとし、光電子のmean escapeだけを陰的closureへ置き換える。
  subroutine measure_tracked_ambient_surface_currents( &
    app, particles, absorbed_flag, area_xy, batch_duration, mpi, electron_current, ion_current &
    )
    type(app_config), intent(in) :: app
    type(particles_soa), intent(in) :: particles
    logical, intent(in) :: absorbed_flag(:)
    real(dp), intent(in) :: area_xy, batch_duration
    type(mpi_context), intent(in) :: mpi
    real(dp), intent(out) :: electron_current, ion_current
    integer(i32) :: particle, species_idx, invalid_species_count
    integer(i32) :: electron_index, ion_index, electron_count, ion_count
    real(dp) :: deposited_charge

    if (size(absorbed_flag) < particles%n .or. area_xy <= 0.0_dp .or. batch_duration <= 0.0_dp) then
      error stop 'tracked ambient current measurement received invalid batch storage or scales.'
    end if
    call find_kinetic_ambient_species(app, electron_index, ion_index, electron_count, ion_count)
    if (electron_count /= 1_i32 .or. ion_count /= 1_i32) then
      error stop 'tracked ambient current measurement requires unique z_high reservoir electron and ion species.'
    end if
    electron_current = 0.0_dp
    ion_current = 0.0_dp
    invalid_species_count = 0_i32
    do particle = 1_i32, particles%n
      if (.not. absorbed_flag(particle)) cycle
      species_idx = particles%species_id(particle)
      if (species_idx < 1_i32 .or. species_idx > app%n_particle_species) then
        invalid_species_count = 1_i32
        cycle
      end if
      deposited_charge = particles%q(particle)*particles%w(particle)
      if (species_idx == electron_index) electron_current = electron_current + deposited_charge
      if (species_idx == ion_index) ion_current = ion_current + deposited_charge
    end do
    call mpi_allreduce_sum_i32_scalar(mpi, invalid_species_count)
    if (invalid_species_count > 0_i32) then
      error stop 'tracked ambient current measurement received an invalid species ID.'
    end if
    call mpi_allreduce_sum_real_dp_scalar(mpi, electron_current)
    call mpi_allreduce_sum_real_dp_scalar(mpi, ion_current)
    electron_current = electron_current/(area_xy*batch_duration)
    ion_current = ion_current/(area_xy*batch_duration)
  end subroutine measure_tracked_ambient_surface_currents

  !> 3D領域から1D平均領域へ実際に渡ったphotoelectron currentを測る。
  !!
  !! 負のtracked chargeがz-highを外向きに通過した量を、全てescapeした場合の
  !! 正のsurface-charging currentへ変換する。粗さ内で再吸収された粒子は
  !! interfaceを通らないため、このmean sourceには含まれない。
  subroutine measure_tracked_photoelectron_interface_current( &
    app, interface_outward_thread, area_xy, batch_duration, mpi, photoelectron_current, photoelectron_index &
    )
    type(app_config), intent(in) :: app
    real(dp), intent(in) :: interface_outward_thread(:, :)
    real(dp), intent(in) :: area_xy, batch_duration
    type(mpi_context), intent(in) :: mpi
    real(dp), intent(out) :: photoelectron_current
    integer(i32), intent(out) :: photoelectron_index
    integer(i32) :: species_idx, photoelectron_count

    if (size(interface_outward_thread, 1) < app%n_particle_species .or. &
        area_xy <= 0.0_dp .or. batch_duration <= 0.0_dp) then
      error stop 'tracked photoelectron current measurement received invalid storage or scales.'
    end if
    photoelectron_index = 0_i32
    photoelectron_count = 0_i32
    do species_idx = 1_i32, app%n_particle_species
      if (.not. app%particle_species(species_idx)%enabled) cycle
      if (trim(lower_ascii(app%particle_species(species_idx)%source_mode)) /= 'photo_raycast') cycle
      if (app%particle_species(species_idx)%q_particle >= 0.0_dp) cycle
      photoelectron_count = photoelectron_count + 1_i32
      photoelectron_index = species_idx
    end do
    if (photoelectron_count /= 1_i32) then
      error stop 'implicit mean requires exactly one negative photo_raycast species.'
    end if

    photoelectron_current = -sum(interface_outward_thread(photoelectron_index, :))
    call mpi_allreduce_sum_real_dp_scalar(mpi, photoelectron_current)
    photoelectron_current = photoelectron_current/(area_xy*batch_duration)
    if (.not. ieee_is_finite(photoelectron_current) .or. photoelectron_current < 0.0_dp) then
      error stop 'tracked photoelectron interface current is not finite and nonnegative.'
    end if
  end subroutine measure_tracked_photoelectron_interface_current

  !> rank 0だけで解いた強PE scalar resultを全rankへ同期する。
  subroutine broadcast_dynamic_zhao_step(mpi, step, effective_source_scale)
    type(mpi_context), intent(in) :: mpi
    type(dynamic_k0_step_type), intent(inout) :: step
    real(dp), intent(inout) :: effective_source_scale
    integer(i32) :: integers(4)
    real(dp) :: values(17)

    integers = 0_i32
    values = 0.0_dp
    if (mpi_is_root(mpi)) then
      integers = [step%status, step%iterations, step%bracket_expansions, int(iachar(step%zhao_branch), i32)]
      values = [ &
               step%surface_charge_before_c, step%surface_charge_after_c, &
               step%interface_potential_before_v, step%interface_potential_after_v, &
               step%interface_field_after_v_m, step%photoelectron_escape_fraction, &
               step%photoelectron_return_fraction, step%backward_euler_residual_charge_c, &
               step%photoelectron_barrier_energy_j, step%photoelectron_energy_max_j, &
               step%marginal_photoelectron_energy_j, step%marginal_photoelectron_escape_fraction, &
               step%photoelectron_source_charge_c, step%zhao_effective_source_scale, &
               step%photoelectron_recross_charge_fraction, &
               step%photoelectron_terminal_mismatch_charge_fraction, effective_source_scale &
               ]
    end if
    call mpi_bcast_i32_array(mpi, integers, 0_i32)
    call mpi_bcast_real_dp_array(mpi, values, 0_i32)
    if (.not. mpi_is_root(mpi)) then
      step = dynamic_k0_step_type()
      step%status = integers(1)
      step%iterations = integers(2)
      step%bracket_expansions = integers(3)
      step%zhao_branch = achar(integers(4))
      step%surface_charge_before_c = values(1)
      step%surface_charge_after_c = values(2)
      step%interface_potential_before_v = values(3)
      step%interface_potential_after_v = values(4)
      step%interface_field_after_v_m = values(5)
      step%photoelectron_escape_fraction = values(6)
      step%photoelectron_return_fraction = values(7)
      step%backward_euler_residual_charge_c = values(8)
      step%photoelectron_barrier_energy_j = values(9)
      step%photoelectron_energy_max_j = values(10)
      step%marginal_photoelectron_energy_j = values(11)
      step%marginal_photoelectron_escape_fraction = values(12)
      step%photoelectron_source_charge_c = values(13)
      step%zhao_effective_source_scale = values(14)
      step%photoelectron_recross_charge_fraction = values(15)
      step%photoelectron_terminal_mismatch_charge_fraction = values(16)
    end if
    effective_source_scale = values(17)
  end subroutine broadcast_dynamic_zhao_step

end submodule bem_simulator_loop

!> `bem_simulator` の主ループと粒子処理計算を実装する submodule。
submodule(bem_simulator) bem_simulator_loop
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use bem_app_config_runtime, only: compute_z_high_box_potential_statistics
  use bem_matching_plane_coupling, only: matching_plane_coupling_type
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
  integer(i32) :: team_size_min, team_size_max
  integer(i32) :: particle_team_size, particle_team_size_min, particle_team_size_max
  integer(i32) :: collision_failure_count, collision_failure_rank, collision_failure_status
  integer(i32) :: collision_failure_particle, collision_failure_step
  integer(i32) :: local_failure_values(3), selected_failure_values(3)
  integer(i32) :: photo_failure_count, photo_failure_rank, photo_failure_status, photo_failure_species
  integer(i32) :: photo_failure_ray, photo_failure_bounce
  integer(i32) :: photo_local_failure_values(4), photo_selected_failure_values(4)
  integer(i32) :: trial_halvings, species_idx, fresh_particle_count
  integer(i32) :: boundary_status
  integer :: hist_unit, pot_hist_unit, top_ref_hist_unit, matching_hist_unit
  integer, allocatable :: rng_state_before(:)
  logical :: history_enabled, potential_history_enabled, top_reference_history_enabled
  logical :: ledger_enabled, adaptive_nonzero_mode, trial_accepted, omp_dynamic_before
  logical :: matching_active, replay_active, matching_history_enabled
  real(dp), allocatable :: potential_buf(:), injection_residual_before(:), boundary_injection_residual_before(:, :)
  integer(i64) :: batch_counts(6), batch_retry_counts(2)
  real(dp) :: bfield(3), rel, t0, sim_t0, batch_t0, batch_soft_discarded_abs_charge
  real(dp) :: collision_failure_x(3), collision_failure_v(3), selected_failure_state(6)
  real(dp) :: trial_batch_duration, duration_ratio, adaptive_potential_step, adaptive_metric_values(1)
  real(dp) :: projected_simulated_time
  real(dp) :: top_phi_mean, top_phi_std, top_phi_min, top_phi_max
  character(len=256) :: boundary_message
  type(particles_soa) :: pcls_batch
  type(mpi_context) :: mpi_ctx
  type(electrostatic_snapshot_type) :: snapshot
  type(electrostatic_diagnostics_type) :: committed_snapshot_diagnostics
  type(periodic_zero_mode_state_type) :: committed_zero_state
  type(charge_ledger_type) :: batch_ledger
  type(matching_plane_coupling_type) :: coupling
  type(simulator_batch_workspace_type) :: workspace
  type(particle_source_plan_type) :: source_plan
  type(external_boundary_contract_type) :: boundary_contract
  type(surface_closure_contract_type) :: surface_closure
  type(field_physics_config) :: field_config
  type(panel_kernel_config) :: panel_config
  type(app_config) :: trial_app
  type(sim_stats) :: stats_candidate

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
    app%sim%reservoir_potential_model, app%sim%open_boundary_model, &
    boundary_contract, boundary_status, boundary_message &
    )
  if (boundary_status /= external_boundary_ok) error stop trim(boundary_message)
  adaptive_nonzero_mode = app%periodic2%max_nonzero_mode_potential_step > 0.0_dp
  call coupling%initialize(app, mesh, stats, mpi_ctx)
  matching_active = coupling%is_active()
  replay_active = adaptive_nonzero_mode .or. matching_active
  omp_dynamic_before = .false.
!$ omp_dynamic_before = omp_get_dynamic()
  nth = 1_i32
  if (replay_active) then
!$  call omp_set_dynamic(.false.)
    !$omp parallel default(none) shared(nth)
    !$omp single
!$  nth = int(omp_get_num_threads(), i32)
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
    if (adaptive_nonzero_mode) stats%adaptive_nonzero_mode_omp_threads = nth
  else
!$  nth = max(1_i32, int(omp_get_max_threads(), i32))
  end if
  if (.not. adaptive_nonzero_mode) then
    stats%adaptive_nonzero_mode_rejected_trials = 0_i64
    stats%adaptive_nonzero_mode_last_batch_duration = 0.0_dp
    stats%adaptive_nonzero_mode_last_potential_step = 0.0_dp
    stats%adaptive_nonzero_mode_omp_threads = 0_i32
  end if
  call validate_soft_discard_initial_state(app%sim, stats)
  call enforce_soft_discard_limits(app%sim, stats, stats%batches, 0_i64, 0_i64, 0.0_dp, mpi_ctx)
  call workspace%init(mesh%nelem, app%n_particle_species, nth, candidate_charge_enabled=adaptive_nonzero_mode)
  call evaluate_surface_closure(app, surface_closure)

  matching_history_enabled = present(matching_plane_history_unit)
  matching_hist_unit = -1
  if (matching_history_enabled) then
    matching_hist_unit = matching_plane_history_unit
    matching_history_enabled = matching_hist_unit /= -1
  end if

  history_enabled = present(history_unit)
  hist_unit = 0
  if (history_enabled) hist_unit = history_unit
  hist_stride = 1_i32
  if (present(history_stride)) then
    matching_history_enabled = matching_history_enabled .and. history_stride > 0_i32
    hist_stride = max(1_i32, history_stride)
  end if
  potential_history_enabled = present(potential_history_unit)
  pot_hist_unit = -1
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
  if (stats%batches < 0_i32 .or. stats%batches > final_batch_idx) then
    error stop 'checkpoint batch count is outside sim.batch_count.'
  end if
  if (.not. adaptive_nonzero_mode .and. stats%batches > 0_i32 .and. stats%simulated_time == 0.0_dp) then
    stats%simulated_time = real(stats%batches, dp)*app%sim%batch_duration
  end if
  batch_count_this_run = final_batch_idx - stats%batches

  call perf_region_begin(perf_region_simulation_total, sim_t0)
  call perf_region_begin(perf_region_field_solver_init, t0)
  call derive_field_panel_config(app%sim, field_config, panel_config)
  call snapshot%init(mesh, app%sim, field_config, app%periodic2, panel_config)
  call perf_region_end(perf_region_field_solver_init, t0)
  call coupling%restore_gauge(mesh, stats, snapshot)

  if (replay_active) then
    call random_seed(size=species_idx)
    allocate (rng_state_before(species_idx))
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

  do local_batch_idx = 1_i32, batch_count_this_run
    call perf_region_begin(perf_region_batch_total, batch_t0)
    batch_idx = stats%batches + 1_i32
    call perf_region_begin(perf_region_field_refresh, t0)
    call snapshot%refresh(mesh)
    call perf_region_end(perf_region_field_refresh, t0)

    trial_batch_duration = app%sim%batch_duration
    trial_halvings = 0_i32
    trial_accepted = .false.
    if (replay_active) then
      call random_seed(get=rng_state_before)
      if (allocated(injection_residual_before)) injection_residual_before = inject_state%macro_residual
      if (allocated(boundary_injection_residual_before)) then
        boundary_injection_residual_before = inject_state%boundary_macro_residual
      end if
      committed_zero_state = snapshot%zero_state
      committed_snapshot_diagnostics = snapshot%diagnostics
    end if

    do while (.not. trial_accepted)
      if (replay_active) then
        call random_seed(put=rng_state_before)
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
        if (trim(lower_ascii(trial_app%particle_species(species_idx)%source_mode)) /= 'reservoir_face' .and. &
            trim(lower_ascii(trial_app%particle_species(species_idx)%source_mode)) /= 'plane_source' .and. &
            .not. any(trial_app%particle_species(species_idx)%boundary_inflow_low /= 0_i32) .and. &
            .not. any(trial_app%particle_species(species_idx)%boundary_inflow_high /= 0_i32)) cycle
        trial_app%particle_species(species_idx)%w_particle = &
          app%particle_species(species_idx)%w_particle*duration_ratio
      end do

      call coupling%begin_trial(mesh, snapshot, stats)

      do
        if (matching_active) then
          call random_seed(put=rng_state_before)
          if (allocated(injection_residual_before)) inject_state%macro_residual = injection_residual_before
          if (allocated(boundary_injection_residual_before)) then
            inject_state%boundary_macro_residual = boundary_injection_residual_before
          end if
        end if
        call coupling%prepare_iteration(app, mesh, snapshot, surface_closure, mpi_ctx, trial_batch_duration)

        call perf_region_begin(perf_region_prepare_batch, t0)
        call build_particle_source_plan( &
          trial_app, source_plan, mpi=mpi_ctx, &
          kinetic_inflow_active=surface_closure%has_inflow_kinetic_map, &
          kinetic_reservoir_potential_v=surface_closure%inflow_reservoir_potential_v, &
          kinetic_access_potential_v=surface_closure%inflow_access_potential_v, &
          kinetic_inflow_face=surface_closure%inflow_kinetic_face, &
          number_flux_override_active=surface_closure%has_inflow_number_flux, &
          number_flux_override_m2_s=surface_closure%inflow_number_flux_m2_s &
          )
        call prepare_batch_state( &
          mesh, trial_app, source_plan, snapshot, stats, batch_idx, workspace, pcls_batch, mpi_ctx, inject_state, &
          photo_failure_status, photo_failure_species, photo_failure_ray, photo_failure_bounce &
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

        call perf_region_begin(perf_region_particle_batch, t0)
        call process_particle_batch( &
          mesh, trial_app, boundary_contract, surface_closure, snapshot, pcls_batch, workspace%dq_thread, &
          workspace%escaped_boundary_flag, workspace%absorbed_flag, workspace%absorbed_element, &
          workspace%soft_discarded_boundary_flag, bfield, batch_idx, mpi_ctx%rank, particle_team_size, &
          collision_failure_status, collision_failure_particle, collision_failure_step, &
          collision_failure_x, collision_failure_v, workspace%matching_plane_moments_thread, &
          batch_retry_counts &
          )
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

        call apply_neutral_return_surface_closure(trial_app, pcls_batch, fresh_particle_count, workspace, mpi_ctx)
        call apply_fixed_surface_current_closure( &
          trial_app, surface_closure, pcls_batch, fresh_particle_count, workspace, mpi_ctx &
          )
        if (coupling%finish_iteration( &
            app, mpi_ctx, workspace%matching_plane_moments_thread, trial_batch_duration, batch_idx)) exit
      end do

      if (adaptive_nonzero_mode) then
        call prepare_adaptive_charge_candidate(mesh, workspace, mpi_ctx)
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

    call perf_region_begin(perf_region_count_outcomes, t0)
    call count_batch_outcomes( &
      pcls_batch, workspace%escaped_boundary_flag, workspace%absorbed_flag, &
      workspace%soft_discarded_boundary_flag, batch_counts, batch_soft_discarded_abs_charge &
      )
    call perf_region_end(perf_region_count_outcomes, t0)
    call perf_region_begin(perf_region_mpi_reduce, t0)
    call mpi_allreduce_sum_i64_array(mpi_ctx, batch_counts)
    call mpi_allreduce_sum_i64_array(mpi_ctx, batch_retry_counts)
    call mpi_allreduce_sum_real_dp_scalar(mpi_ctx, batch_soft_discarded_abs_charge)
    call perf_region_end(perf_region_mpi_reduce, t0)
    call enforce_soft_discard_limits( &
      app%sim, stats, batch_idx, batch_counts(6), batch_counts(1), batch_soft_discarded_abs_charge, mpi_ctx &
      )
    call report_soft_discard_summary(batch_idx, batch_counts(6), batch_soft_discarded_abs_charge, mpi_ctx)
    projected_simulated_time = stats%simulated_time + trial_batch_duration
    if (.not. ieee_is_finite(stats%simulated_time) .or. stats%simulated_time < 0.0_dp .or. &
        .not. ieee_is_finite(projected_simulated_time)) then
      error stop 'simulation statistic overflow: simulated_time'
    end if
    ! Build the complete statistics update before mutating mesh/ledger state.  The
    ! checked add routines are the single overflow guard; commit only after all pass.
    stats_candidate = stats
    call accumulate_batch_stats( &
      stats_candidate, batch_counts, batch_soft_discarded_abs_charge, batch_retry_counts, 0.0_dp &
      )
    stats_candidate%simulated_time = projected_simulated_time
    if (adaptive_nonzero_mode) then
      stats_candidate%adaptive_nonzero_mode_last_batch_duration = trial_batch_duration
      stats_candidate%adaptive_nonzero_mode_last_potential_step = adaptive_potential_step
      stats_candidate%adaptive_nonzero_mode_rejected_trials = checked_add_adaptive_rejected_trials( &
                                                              stats%adaptive_nonzero_mode_rejected_trials, &
                                                              trial_halvings &
                                                              )
    end if
    call coupling%stage_stats(stats_candidate)

    if (ledger_enabled) then
      call batch_ledger%reset(batch_idx)
      batch_ledger%surface_charge_before = finite_charge_sum(mesh%q_elem, 'batch surface charge before commit')
      call record_batch_initial_charge(trial_app, pcls_batch, fresh_particle_count, batch_ledger)
      call record_batch_outcome_charge( &
        pcls_batch, workspace%escaped_boundary_flag, workspace%absorbed_flag, &
        workspace%soft_discarded_boundary_flag, batch_ledger &
        )
      call reduce_charge_ledger_fluxes(batch_ledger, mpi_ctx, workspace)
      batch_ledger%neutral_return_correction = workspace%neutral_return_correction
      batch_ledger%neutral_return_weight_scale = workspace%neutral_return_weight_scale
      batch_ledger%neutral_return_unresolved_fraction = workspace%neutral_return_unresolved_fraction
      batch_ledger%fixed_absorbed_target_charge = workspace%fixed_absorbed_target_charge
      batch_ledger%fixed_absorbed_weight_scale = workspace%fixed_absorbed_weight_scale
      batch_ledger%fixed_emission_target_charge = workspace%fixed_emission_target_charge
      batch_ledger%fixed_emission_weight_scale = workspace%fixed_emission_weight_scale
      batch_ledger%fixed_escape_target_charge = workspace%fixed_escape_target_charge
      batch_ledger%fixed_escape_correction = workspace%fixed_escape_correction
      batch_ledger%fixed_current_correction = workspace%fixed_current_correction
    end if

    call perf_region_begin(perf_region_commit_charge, t0)
    call commit_batch_charge( &
      mesh, app%sim%q_floor, app%sim%e0, app%sim%field_bc_mode, workspace, rel, mpi_ctx &
      )
    call perf_region_end(perf_region_commit_charge, t0)
    call coupling%commit(mesh, stats_candidate)
    stats_candidate%last_rel_change = rel
    if (ledger_enabled) then
      batch_ledger%surface_charge_after = finite_charge_sum(mesh%q_elem, 'batch surface charge after commit')
      call accumulate_charge_ledger(charge_ledger, batch_ledger)
    end if

    call perf_region_begin(perf_region_stats_update, t0)
    stats = stats_candidate
    call perf_region_end(perf_region_stats_update, t0)

    call perf_region_begin(perf_region_history_write, t0)
    if (mpi_is_root(mpi_ctx)) then
      call print_batch_progress(batch_idx, final_batch_idx, rel)
      call maybe_write_history_snapshot(history_enabled, hist_unit, hist_stride, stats, rel, mesh%q_elem)
      if (matching_history_enabled .and. mod(batch_idx - 1_i32, hist_stride) == 0_i32) then
        call write_matching_plane_history_snapshot(matching_hist_unit, batch_idx, stats%simulated_time, stats)
      end if
      if (potential_history_enabled) then
        call maybe_write_potential_history_snapshot( &
          potential_history_enabled, pot_hist_unit, hist_stride, stats, snapshot, mesh, app%sim, potential_buf, &
          top_reference_history_enabled, top_ref_hist_unit &
          )
      end if
    end if
    call perf_region_end(perf_region_history_write, t0)
    call maybe_write_periodic_checkpoint(app, mesh, stats, inject_state, mpi_ctx, charge_ledger)
    call perf_region_end(perf_region_batch_total, batch_t0)
  end do
  call perf_region_end(perf_region_simulation_total, sim_t0)
!$ if (replay_active) call omp_set_dynamic(omp_dynamic_before)

  if (present(mesh_potential_v) .and. mpi_is_root(mpi_ctx)) then
    call snapshot%refresh(mesh)
    allocate (mesh_potential_v(mesh%nelem))
    call snapshot%compute_mesh_potential(mesh, app%sim, mesh_potential_v)
  end if
  if (present(electrostatic_diagnostics)) then
    call snapshot%get_diagnostics(electrostatic_diagnostics)
    if (mpi_is_root(mpi_ctx) .and. app%sim%use_box) then
      call snapshot%refresh(mesh)
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
  end procedure run_absorption_insulator

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
      ' rank=', rank, ' status=', trim(particle_failure_name(status)), ' code=', status, &
      ' dt=', dt, ' x=', x, ' v=', v
    flush (error_unit)
    error stop 'particle step failed.'
  end subroutine stop_for_collision_failure

  subroutine stop_for_photo_collision_failure(batch_idx, rank, species, ray, bounce, status)
    integer(i32), intent(in) :: batch_idx, rank, species, ray, bounce, status
    write (error_unit, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,a,a,i0)') &
      'photo ray collision failed: batch=', batch_idx, ' rank=', rank, ' species=', species, &
      ' ray=', ray, ' bounce=', bounce, ' status=', trim(collision_failure_name(status)), ' code=', status
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
        ' warning_threshold_C=', sim%multiple_box_events_soft_discard_abs_charge_limit
      flush (error_unit)
    end if
    if (.not. fraction_exceeded) return

    if (mpi_is_root(mpi)) then
      write (error_unit, '(a,i0,a,i0,a,i0,a,i0,a,es13.5,a,es13.5,a,es13.5,a,es13.5)') &
        'multiple_box_events soft-discard cumulative fraction limit exceeded: batch=', batch_idx, &
        ' count=', projected_count, ' count_grace=', count_grace, &
        ' processed=', projected_processed, ' fraction=', projected_fraction, &
        ' fraction_limit=', sim%multiple_box_events_soft_discard_fraction_limit, &
        ' abs_charge_C=', projected_abs_charge, &
        ' abs_charge_warning_threshold_C=', sim%multiple_box_events_soft_discard_abs_charge_limit
      flush (error_unit)
    end if
    error stop 'multiple_box_events soft-discard cumulative fraction limit exceeded.'
  end subroutine enforce_soft_discard_limits

  pure function particle_failure_name(status) result(name)
    integer(i32), intent(in) :: status
    character(len=32) :: name

    select case (status)
    case (collision_query_image_limit)
      name = 'image_limit'
    case (collision_query_index_range)
      name = 'index_range'
    case (collision_query_invalid_segment)
      name = 'invalid_segment'
    case (collision_query_grid_stalled)
      name = 'grid_stalled'
    case (particle_step_invalid_boundary)
      name = 'invalid_boundary'
    case (particle_step_multiple_box_events)
      name = 'multiple_box_events'
    case (particle_step_ambiguous_open_corner)
      name = 'ambiguous_open_corner'
    case default
      name = 'unknown'
    end select
  end function particle_failure_name

  pure function collision_failure_name(status) result(name)
    integer(i32), intent(in) :: status
    character(len=16) :: name

    select case (status)
    case (collision_query_image_limit)
      name = 'image_limit'
    case (collision_query_index_range)
      name = 'index_range'
    case (collision_query_invalid_segment)
      name = 'invalid_segment'
    case (collision_query_grid_stalled)
      name = 'grid_stalled'
    case default
      name = 'unknown'
    end select
  end function collision_failure_name

end submodule bem_simulator_loop

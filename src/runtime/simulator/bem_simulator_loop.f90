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
  use bem_mpi, only: mpi_allreduce_max_real_dp_array
  implicit none
  real(dp), parameter :: neutral_return_max_unresolved_fraction = 0.05_dp
  integer(i32), parameter :: adaptive_max_halvings = 24_i32
contains

  module procedure run_absorption_insulator
  integer(i32) :: batch_idx, final_batch_idx, batch_count_this_run, local_batch_idx, nth, hist_stride
  integer(i32) :: collision_failure_count, collision_failure_rank, collision_failure_status
  integer(i32) :: collision_failure_particle, collision_failure_step
  integer(i32) :: local_failure_values(3), selected_failure_values(3)
  integer(i32) :: photo_failure_count, photo_failure_rank, photo_failure_status, photo_failure_species
  integer(i32) :: photo_failure_ray, photo_failure_bounce
  integer(i32) :: photo_local_failure_values(4), photo_selected_failure_values(4)
  integer(i32) :: trial_halvings, species_idx, fresh_particle_count
  integer(i32) :: boundary_status
  integer :: hist_unit, pot_hist_unit, top_ref_hist_unit
  integer, allocatable :: rng_state_before(:)
  logical :: history_enabled, potential_history_enabled, top_reference_history_enabled
  logical :: ledger_enabled, adaptive_nonzero_mode, trial_accepted
  real(dp), allocatable :: potential_buf(:), injection_residual_before(:), boundary_injection_residual_before(:, :)
  integer(i32) :: batch_counts(6)
  real(dp) :: bfield(3), rel, t0, sim_t0, batch_t0, batch_soft_discarded_abs_charge
  real(dp) :: collision_failure_x(3), collision_failure_v(3), selected_failure_state(6)
  real(dp) :: trial_batch_duration, duration_ratio, adaptive_potential_step, adaptive_metric_values(1)
  real(dp) :: top_phi_mean, top_phi_std, top_phi_min, top_phi_max
  character(len=256) :: boundary_message
  type(particles_soa) :: pcls_batch
  type(mpi_context) :: mpi_ctx
  type(electrostatic_snapshot_type) :: snapshot
  type(electrostatic_diagnostics_type) :: committed_snapshot_diagnostics
  type(periodic_zero_mode_state_type) :: committed_zero_state
  type(charge_ledger_type) :: batch_ledger
  type(simulator_batch_workspace_type) :: workspace
  type(particle_source_plan_type) :: source_plan
  type(external_boundary_contract_type) :: boundary_contract
  type(app_config) :: trial_app

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
  nth = 1_i32
!$ nth = max(1, omp_get_max_threads())
  call workspace%init(mesh%nelem, app%n_particle_species, nth, candidate_charge_enabled=adaptive_nonzero_mode)

  history_enabled = present(history_unit)
  hist_unit = 0
  if (history_enabled) hist_unit = history_unit
  hist_stride = 1_i32
  if (present(history_stride)) hist_stride = max(1_i32, history_stride)
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
  call snapshot%init(mesh, app%sim, app%field, app%periodic2, app%panel)
  call perf_region_end(perf_region_field_solver_init, t0)

  if (adaptive_nonzero_mode) then
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

  do local_batch_idx = 1_i32, batch_count_this_run
    call perf_region_begin(perf_region_batch_total, batch_t0)
    batch_idx = stats%batches + 1_i32
    call perf_region_begin(perf_region_field_refresh, t0)
    call snapshot%refresh(mesh)
    call perf_region_end(perf_region_field_refresh, t0)

    trial_batch_duration = app%sim%batch_duration
    trial_halvings = 0_i32
    trial_accepted = .false.
    if (adaptive_nonzero_mode) then
      call random_seed(get=rng_state_before)
      if (allocated(injection_residual_before)) injection_residual_before = inject_state%macro_residual
      if (allocated(boundary_injection_residual_before)) then
        boundary_injection_residual_before = inject_state%boundary_macro_residual
      end if
      committed_zero_state = snapshot%zero_state
      committed_snapshot_diagnostics = snapshot%diagnostics
    end if

    do while (.not. trial_accepted)
      if (adaptive_nonzero_mode) then
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

      call perf_region_begin(perf_region_prepare_batch, t0)
      call build_particle_source_plan(trial_app, source_plan, mpi=mpi_ctx)
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
        mesh, trial_app, boundary_contract, snapshot, pcls_batch, workspace%dq_thread, &
        workspace%escaped_boundary_flag, workspace%absorbed_flag, workspace%absorbed_element, &
        workspace%soft_discarded_boundary_flag, bfield, batch_idx, mpi_ctx%rank, &
        collision_failure_status, collision_failure_particle, collision_failure_step, &
        collision_failure_x, collision_failure_v &
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

      call apply_neutral_return_surface_closure(trial_app, pcls_batch, fresh_particle_count, workspace, mpi_ctx)
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

    if (ledger_enabled) then
      call batch_ledger%reset(batch_idx)
      batch_ledger%surface_charge_before = sum(mesh%q_elem)
      call record_batch_initial_charge(trial_app, pcls_batch, fresh_particle_count, batch_ledger)
      call record_batch_outcome_charge( &
        pcls_batch, workspace%escaped_boundary_flag, workspace%absorbed_flag, &
        workspace%soft_discarded_boundary_flag, batch_ledger &
        )
      call reduce_charge_ledger_fluxes(batch_ledger, mpi_ctx, workspace)
      batch_ledger%neutral_return_correction = workspace%neutral_return_correction
      batch_ledger%neutral_return_weight_scale = workspace%neutral_return_weight_scale
      batch_ledger%neutral_return_unresolved_fraction = workspace%neutral_return_unresolved_fraction
    end if

    call perf_region_begin(perf_region_commit_charge, t0)
    call commit_batch_charge( &
      mesh, app%sim%q_floor, app%sim%e0, app%sim%field_bc_mode, workspace, rel, mpi_ctx &
      )
    call perf_region_end(perf_region_commit_charge, t0)
    if (ledger_enabled) then
      batch_ledger%surface_charge_after = sum(mesh%q_elem)
      call accumulate_charge_ledger(charge_ledger, batch_ledger)
    end if

    call perf_region_begin(perf_region_count_outcomes, t0)
    call count_batch_outcomes( &
      pcls_batch, workspace%escaped_boundary_flag, workspace%absorbed_flag, &
      workspace%soft_discarded_boundary_flag, batch_counts, batch_soft_discarded_abs_charge &
      )
    call perf_region_end(perf_region_count_outcomes, t0)
    call perf_region_begin(perf_region_mpi_reduce, t0)
    call mpi_allreduce_sum_i32_array(mpi_ctx, batch_counts)
    call perf_region_end(perf_region_mpi_reduce, t0)
    call perf_region_begin(perf_region_stats_update, t0)
    call accumulate_batch_stats(stats, batch_counts, batch_soft_discarded_abs_charge, rel)
    stats%simulated_time = stats%simulated_time + trial_batch_duration
    stats%adaptive_nonzero_mode_last_batch_duration = trial_batch_duration
    stats%adaptive_nonzero_mode_last_potential_step = merge(adaptive_potential_step, 0.0_dp, adaptive_nonzero_mode)
    stats%adaptive_nonzero_mode_rejected_trials = &
      stats%adaptive_nonzero_mode_rejected_trials + int(trial_halvings, i64)
    call perf_region_end(perf_region_stats_update, t0)

    call perf_region_begin(perf_region_history_write, t0)
    if (mpi_is_root(mpi_ctx)) then
      call print_batch_progress(batch_idx, final_batch_idx, rel)
      call maybe_write_history_snapshot(history_enabled, hist_unit, hist_stride, stats, rel, mesh%q_elem)
      if (potential_history_enabled) then
        call maybe_write_potential_history_snapshot( &
          potential_history_enabled, pot_hist_unit, hist_stride, stats, snapshot, mesh, app%sim, potential_buf, &
          top_reference_history_enabled, top_ref_hist_unit &
          )
      end if
    end if
    call perf_region_end(perf_region_history_write, t0)
    call perf_region_end(perf_region_batch_total, batch_t0)
  end do
  call perf_region_end(perf_region_simulation_total, sim_t0)

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

  module procedure prepare_batch_state
  batch_idx = stats%batches + 1_i32
  call workspace%reset_before_injection()
  if (present(inject_state)) then
    call init_particle_batch_from_config( &
      app, batch_idx, pcls_batch, inject_state, mesh=mesh, photo_emission_dq=workspace%photo_emission_dq, &
      mpi=mpi, collision_failure_status=collision_failure_status, &
      collision_failure_species=collision_failure_species, collision_failure_ray=collision_failure_ray, &
      collision_failure_bounce=collision_failure_bounce, snapshot=snapshot, source_plan=source_plan &
      )
  else
    call init_particle_batch_from_config( &
      app, batch_idx, pcls_batch, mesh=mesh, photo_emission_dq=workspace%photo_emission_dq, mpi=mpi, &
      collision_failure_status=collision_failure_status, collision_failure_species=collision_failure_species, &
      collision_failure_ray=collision_failure_ray, collision_failure_bounce=collision_failure_bounce, &
      snapshot=snapshot, source_plan=source_plan &
      )
  end if
  if (collision_failure_status /= collision_query_ok) return
  call workspace%prepare_particle_flags(pcls_batch%n)
  end procedure prepare_batch_state

  module procedure process_particle_batch
  integer(i32) :: i, step, tid, nth, collision_status, species_idx
  real(dp) :: x0(3), v0(3), x1(3), v1(3), qdep
  type(hit_info) :: hit
  type(particle_step_result) :: step_result
  type(sim_config) :: particle_sim
  logical :: candidate_inside, used_event_resolver

  nth = size(dq_thread, 2)
  collision_failure_status = collision_query_ok
  collision_failure_particle = huge(0_i32)
  collision_failure_step = 0_i32
  collision_failure_x = 0.0_dp
  collision_failure_v = 0.0_dp

  !$omp parallel default(none) &
  !$omp shared(mesh,pcls_batch,app,boundary_contract,snapshot,dq_thread,bfield,escaped_boundary_flag,absorbed_flag,nth) &
  !$omp shared(absorbed_element,soft_discarded_boundary_flag,batch_idx,mpi_rank) &
  !$omp shared(collision_failure_status,collision_failure_particle,collision_failure_step) &
  !$omp shared(collision_failure_x,collision_failure_v) &
  !$omp private(i,step,x0,v0,x1,v1,hit,step_result,particle_sim,tid,qdep,species_idx) &
  !$omp private(collision_status,candidate_inside,used_event_resolver)
  tid = 1_i32
!$ tid = omp_get_thread_num() + 1
  !$omp do schedule(runtime)
  do i = 1_i32, pcls_batch%n
    if (.not. pcls_batch%alive(i)) cycle
    species_idx = pcls_batch%species_id(i)
    particle_sim = app%sim
    call resolve_particle_boundaries( &
      app%sim, app%particle_boundary_low, app%particle_boundary_high, app%particle_species(species_idx), &
      particle_sim%bc_low, particle_sim%bc_high &
      )
    do step = 1_i32, app%sim%max_step
      x0 = pcls_batch%x(:, i)
      v0 = pcls_batch%v(:, i)
      call build_particle_step_candidate( &
        mesh, particle_sim, snapshot, bfield, x0, v0, &
        pcls_batch%q(i), pcls_batch%m(i), app%sim%dt, x1, v1 &
        )
      if (.not. all(ieee_is_finite(x1)) .or. .not. all(ieee_is_finite(v1))) then
        call record_collision_failure(particle_step_invalid_boundary, i, step, x0, v0)
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
            result=step_result, boundary_contract=boundary_contract, &
            boundary_rng_counter=int([batch_idx, mpi_rank, i, step], i64) &
            )
          used_event_resolver = .true.
        else
          call record_collision_failure(collision_status, i, step, x0, v0)
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
          hit=hit, result=step_result, boundary_contract=boundary_contract, &
          boundary_rng_counter=int([batch_idx, mpi_rank, i, step], i64) &
          )
        used_event_resolver = .true.
      end if
      if (used_event_resolver) then
        if (step_result%status /= collision_query_ok) then
          if (step_result%status == particle_step_multiple_box_events .and. &
              trim(lower_ascii(app%sim%multiple_box_events_policy)) == 'soft_discard') then
            pcls_batch%alive(i) = .false.
            soft_discarded_boundary_flag(i) = .true.
            exit
          end if
          call record_collision_failure(step_result%status, i, step, x0, v0)
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
    end do
  end do
  !$omp end do
  !$omp end parallel

contains
  subroutine record_collision_failure(status, particle_index, particle_step, failure_x, failure_v)
    integer(i32), intent(in) :: status, particle_index, particle_step
    real(dp), intent(in) :: failure_x(3), failure_v(3)
    !$omp critical (beach_collision_query_failure)
    if (collision_failure_status == collision_query_ok .or. particle_index < collision_failure_particle .or. &
        (particle_index == collision_failure_particle .and. particle_step < collision_failure_step)) then
      collision_failure_status = status
      collision_failure_particle = particle_index
      collision_failure_step = particle_step
      collision_failure_x = failure_x
      collision_failure_v = failure_v
    end if
    !$omp end critical (beach_collision_query_failure)
  end subroutine record_collision_failure
  end procedure process_particle_batch

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

  module procedure commit_batch_charge
  real(dp) :: norm_dq, norm_q
  workspace%q_before = mesh%q_elem
  if (workspace%charge_candidate_ready) then
    mesh%q_elem = workspace%candidate_charge
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

  subroutine prepare_adaptive_charge_candidate(mesh, workspace, mpi)
    type(mesh_type), intent(in) :: mesh
    type(simulator_batch_workspace_type), intent(inout) :: workspace
    type(mpi_context), intent(in) :: mpi
    if (size(workspace%candidate_charge) /= mesh%nelem) then
      error stop 'adaptive candidate storage does not match the mesh.'
    end if
    workspace%dq = sum(workspace%dq_thread, dim=2) + workspace%photo_emission_dq
    call mpi_allreduce_sum_real_dp_array(mpi, workspace%dq)
    workspace%candidate_charge = mesh%q_elem + workspace%dq
    if (.not. all(ieee_is_finite(workspace%candidate_charge))) then
      error stop 'adaptive candidate charge is not finite.'
    end if
    workspace%charge_candidate_ready = .true.
  end subroutine prepare_adaptive_charge_candidate

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

    n = app%n_particle_species
    has_neutral_return = .false.
    do species_idx = 1_i32, n
      if (.not. app%particle_species(species_idx)%enabled) cycle
      if (trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) == 'neutral_return') then
        has_neutral_return = .true.
      end if
    end do
    if (.not. has_neutral_return) return
    if (.not. app%sim%use_box) error stop 'neutral_return requires a finite box.'

    workspace%neutral_return_charge_values = 0.0_dp
    workspace%neutral_return_terminal_counts = 0_i64
    do i = 1_i32, pcls_batch%n
      species_idx = pcls_batch%species_id(i)
      if (trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) /= 'neutral_return') cycle
      macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
      if (.not. ieee_is_finite(macro_charge) .or. macro_charge >= 0.0_dp .or. i > fresh_particle_count) then
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
                       merge(1_i32, 0_i32, pcls_batch%alive(i))
      if (terminal_count /= 1_i32) then
        workspace%neutral_return_charge_values(5*n + species_idx) = &
          workspace%neutral_return_charge_values(5*n + species_idx) + macro_charge
        workspace%neutral_return_terminal_counts(2*n + species_idx) = &
          workspace%neutral_return_terminal_counts(2*n + species_idx) + 1_i64
      else if (workspace%absorbed_flag(i)) then
        elem_idx = workspace%absorbed_element(i)
        if (elem_idx >= 1_i32 .and. elem_idx <= size(workspace%dq_thread, 1)) then
          workspace%neutral_return_charge_values(n + species_idx) = &
            workspace%neutral_return_charge_values(n + species_idx) + macro_charge
        else
          workspace%neutral_return_charge_values(5*n + species_idx) = &
            workspace%neutral_return_charge_values(5*n + species_idx) + macro_charge
          workspace%neutral_return_terminal_counts(2*n + species_idx) = &
            workspace%neutral_return_terminal_counts(2*n + species_idx) + 1_i64
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
      if (trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) /= 'neutral_return') cycle
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
        write (error_unit, '(a,i0,3(a,i0),3(a,es13.5))') &
          'neutral_return has unsupported terminal outcome for species ', species_idx, &
          ': escaped=', escaped_count, ' soft=', soft_count, ' invalid=', invalid_count, &
          ' escaped_C=', escaped_charge, ' soft_C=', soft_charge, ' invalid_C=', invalid_charge
        error stop 'neutral_return terminal outcome is unsupported.'
      end if
      charge_scale = max(abs(emitted_charge), abs(absorbed_charge), abs(unresolved_charge), tiny(1.0_dp))
      charge_tolerance = 4096.0_dp*epsilon(1.0_dp)*charge_scale
      balance_tolerance = max(charge_tolerance, sqrt(epsilon(1.0_dp))*charge_scale)
      if (abs(emitted_charge) <= charge_tolerance) cycle
      if (emitted_charge >= 0.0_dp .or. absorbed_charge >= -charge_tolerance .or. &
          abs(emitted_charge - absorbed_charge - unresolved_charge) > balance_tolerance) then
        error stop 'neutral_return charge balance is invalid.'
      end if
      weight_scale = emitted_charge/absorbed_charge
      correction_charge = emitted_charge - absorbed_charge
      unresolved_fraction = unresolved_charge/emitted_charge
      if (.not. all(ieee_is_finite([weight_scale, correction_charge, unresolved_fraction])) .or. &
          unresolved_fraction > neutral_return_max_unresolved_fraction + sqrt(epsilon(1.0_dp))) then
        error stop 'neutral_return unresolved fraction exceeds the applicability limit.'
      end if
      workspace%neutral_return_weight_scale(species_idx) = max(1.0_dp, weight_scale)
      workspace%neutral_return_correction(species_idx) = correction_charge
      workspace%neutral_return_unresolved_fraction(species_idx) = max(0.0_dp, unresolved_fraction)
    end do

    do i = 1_i32, fresh_particle_count
      species_idx = pcls_batch%species_id(i)
      if (trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) /= 'neutral_return') cycle
      if (.not. workspace%absorbed_flag(i)) cycle
      elem_idx = workspace%absorbed_element(i)
      macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
      workspace%dq_thread(elem_idx, 1) = workspace%dq_thread(elem_idx, 1) + &
                                         (workspace%neutral_return_weight_scale(species_idx) - 1.0_dp)*macro_charge
    end do
  end subroutine apply_neutral_return_surface_closure

  subroutine record_batch_initial_charge(app, pcls_batch, fresh_particle_count, ledger)
    type(app_config), intent(in) :: app
    type(particles_soa), intent(in) :: pcls_batch
    integer(i32), intent(in) :: fresh_particle_count
    type(charge_ledger_type), intent(inout) :: ledger
    integer(i32) :: i, species_idx
    real(dp) :: macro_charge
    do i = 1_i32, fresh_particle_count
      species_idx = pcls_batch%species_id(i)
      macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
      if (trim(lower_ascii(app%particle_species(species_idx)%source_mode)) == 'photo_raycast') then
        ledger%emitted_from_surface(species_idx) = ledger%emitted_from_surface(species_idx) + macro_charge
        ledger%emitted_count(species_idx) = ledger%emitted_count(species_idx) + 1_i64
      else
        ledger%injected_from_remote(species_idx) = ledger%injected_from_remote(species_idx) + macro_charge
        ledger%injected_count(species_idx) = ledger%injected_count(species_idx) + 1_i64
      end if
    end do
  end subroutine record_batch_initial_charge

  subroutine record_batch_outcome_charge( &
    pcls_batch, escaped_boundary_flag, absorbed_flag, soft_discarded_boundary_flag, ledger &
    )
    type(particles_soa), intent(in) :: pcls_batch
    logical, intent(in) :: escaped_boundary_flag(:), absorbed_flag(:), soft_discarded_boundary_flag(:)
    type(charge_ledger_type), intent(inout) :: ledger
    integer(i32) :: i, species_idx
    real(dp) :: macro_charge
    do i = 1_i32, pcls_batch%n
      species_idx = pcls_batch%species_id(i)
      macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
      if (absorbed_flag(i)) then
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

  subroutine reduce_charge_ledger_fluxes(ledger, mpi, workspace)
    type(charge_ledger_type), intent(inout) :: ledger
    type(mpi_context), intent(in) :: mpi
    type(simulator_batch_workspace_type), intent(inout) :: workspace
    integer(i32) :: n
    n = ledger%nspecies
    workspace%ledger_charge_values = [ &
                                     ledger%injected_from_remote, ledger%emitted_from_surface, ledger%absorbed_on_surface, &
                                     ledger%escaped_to_infinity, ledger%discarded_unresolved &
                                     ]
    call mpi_allreduce_sum_real_dp_array(mpi, workspace%ledger_charge_values)
    ledger%injected_from_remote = workspace%ledger_charge_values(1:n)
    ledger%emitted_from_surface = workspace%ledger_charge_values(n + 1:2*n)
    ledger%absorbed_on_surface = workspace%ledger_charge_values(2*n + 1:3*n)
    ledger%escaped_to_infinity = workspace%ledger_charge_values(3*n + 1:4*n)
    ledger%discarded_unresolved = workspace%ledger_charge_values(4*n + 1:5*n)
    workspace%ledger_count_values = [ &
                                    ledger%injected_count, ledger%emitted_count, ledger%absorbed_count, ledger%escaped_count, &
                                    ledger%discarded_unresolved_count &
                                    ]
    call mpi_allreduce_sum_i64_array(mpi, workspace%ledger_count_values)
    ledger%injected_count = workspace%ledger_count_values(1:n)
    ledger%emitted_count = workspace%ledger_count_values(n + 1:2*n)
    ledger%absorbed_count = workspace%ledger_count_values(2*n + 1:3*n)
    ledger%escaped_count = workspace%ledger_count_values(3*n + 1:4*n)
    ledger%discarded_unresolved_count = workspace%ledger_count_values(4*n + 1:5*n)
  end subroutine reduce_charge_ledger_fluxes

end submodule bem_simulator_loop

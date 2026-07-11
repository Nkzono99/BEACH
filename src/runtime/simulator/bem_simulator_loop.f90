!> `bem_simulator` の主ループと粒子処理計算を実装する submodule。
submodule(bem_simulator) bem_simulator_loop
  use, intrinsic :: iso_fortran_env, only: output_unit
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
  real(dp), allocatable :: dq_thread(:, :), dq(:), photo_emission_dq(:)
  real(dp), allocatable :: interface_outward_thread(:, :), interface_returned_thread(:, :)
  real(dp), allocatable :: interface_tau_max_thread(:), interface_frozen_ratio_max_thread(:)
  logical, allocatable :: escaped_boundary_flag(:), absorbed_flag(:)
  integer(i32) :: batch_counts(5)
  real(dp) :: bfield(3), rel, t0, sim_t0, batch_t0
  real(dp) :: collision_failure_x(3), collision_failure_v(3), selected_failure_state(6)
  real(dp) :: max_outer_flight_time, max_frozen_field_ratio
  type(particles_soa) :: pcls_batch
  type(mpi_context) :: mpi_ctx
  type(electrostatic_snapshot_type) :: snapshot
  type(outer_coupler_type) :: outer_coupler
  type(charge_ledger_type) :: batch_ledger
  logical :: ledger_enabled

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

  nth = 1
!$ nth = max(1, omp_get_max_threads())
  allocate (dq_thread(mesh%nelem, nth), dq(mesh%nelem), photo_emission_dq(mesh%nelem))
  allocate (interface_outward_thread(app%n_particle_species, nth))
  allocate (interface_returned_thread(app%n_particle_species, nth))
  allocate (interface_tau_max_thread(nth), interface_frozen_ratio_max_thread(nth))
  max_outer_flight_time = 0.0_dp
  max_frozen_field_ratio = 0.0_dp

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
  call snapshot%init(mesh, app%sim, app%field, app%periodic2, app%panel, app%outer_plasma)
  if (present(electrostatic_restart_state)) then
    call snapshot%restore_outer_state(electrostatic_restart_state)
    call outer_coupler%init(app%coupling, electrostatic_restart_state%last_outer_update_batch)
  else
    call outer_coupler%init(app%coupling)
  end if
  call perf_region_end(perf_region_field_solver_init, t0)

  do local_batch_idx = 1, batch_count_this_run
    call perf_region_begin(perf_region_batch_total, batch_t0)

    call perf_region_begin(perf_region_prepare_batch, t0)
    call prepare_batch_state( &
      mesh, app, stats, batch_idx, dq_thread, pcls_batch, escaped_boundary_flag, absorbed_flag, &
      photo_emission_dq, mpi_ctx, inject_state, photo_failure_status, photo_failure_species, &
      photo_failure_ray, photo_failure_bounce &
      )
    call perf_region_end(perf_region_prepare_batch, t0)

    if (ledger_enabled) then
      call batch_ledger%reset(batch_idx)
      batch_ledger%surface_charge_before = sum(mesh%q_elem)
      call record_batch_initial_charge(app, pcls_batch, batch_ledger)
    end if

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

    call perf_region_begin(perf_region_field_refresh, t0)
    call outer_coupler%refresh(snapshot, mesh, batch_idx)
    call perf_region_end(perf_region_field_refresh, t0)

    call perf_region_begin(perf_region_particle_batch, t0)
    call process_particle_batch( &
      mesh, app, snapshot, pcls_batch, dq_thread, escaped_boundary_flag, absorbed_flag, bfield, batch_idx, mpi_ctx%rank, &
      interface_outward_thread, interface_returned_thread, collision_failure_status, collision_failure_particle, &
      collision_failure_step, collision_failure_x, collision_failure_v, interface_tau_max_thread, &
      interface_frozen_ratio_max_thread &
      )
    max_outer_flight_time = max(max_outer_flight_time, maxval(interface_tau_max_thread))
    max_frozen_field_ratio = max(max_frozen_field_ratio, maxval(interface_frozen_ratio_max_thread))
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

    if (ledger_enabled) then
      batch_ledger%interface_outward_gross = sum(interface_outward_thread, dim=2)
      batch_ledger%interface_returned_gross = sum(interface_returned_thread, dim=2)
      call record_batch_outcome_charge(pcls_batch, escaped_boundary_flag, absorbed_flag, batch_ledger)
      call reduce_charge_ledger_fluxes(batch_ledger, mpi_ctx)
    end if

    call perf_region_begin(perf_region_commit_charge, t0)
    call commit_batch_charge( &
      mesh, app%sim%q_floor, app%sim%softening, app%sim%e0, app%sim%field_bc_mode, &
      dq_thread, photo_emission_dq, dq, rel, mpi_ctx &
      )
    call perf_region_end(perf_region_commit_charge, t0)

    if (ledger_enabled) then
      batch_ledger%surface_charge_after = sum(mesh%q_elem)
      call accumulate_charge_ledger(charge_ledger, batch_ledger)
    end if

    call perf_region_begin(perf_region_count_outcomes, t0)
    call count_batch_outcomes(pcls_batch, escaped_boundary_flag, absorbed_flag, batch_counts)
    call perf_region_end(perf_region_count_outcomes, t0)

    call perf_region_begin(perf_region_mpi_reduce, t0)
    call mpi_allreduce_sum_i32_array(mpi_ctx, batch_counts)
    call perf_region_end(perf_region_mpi_reduce, t0)

    call perf_region_begin(perf_region_stats_update, t0)
    call accumulate_batch_stats(stats, batch_counts, rel)
    call perf_region_end(perf_region_stats_update, t0)

    call perf_region_begin(perf_region_history_write, t0)
    if (mpi_is_root(mpi_ctx)) then
      call print_batch_progress(batch_idx, final_batch_idx, rel)
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
  if (present(electrostatic_diagnostics)) then
    if (.not. present(mesh_potential_v) .and. mpi_is_root(mpi_ctx)) then
      call snapshot%refresh(mesh, update_outer=.false.)
    end if
    call snapshot%get_diagnostics(electrostatic_diagnostics)
    electrostatic_diagnostics%last_outer_update_batch = outer_coupler%last_outer_update_batch
    electrostatic_diagnostics%max_outer_flight_time = max_outer_flight_time
    electrostatic_diagnostics%max_frozen_field_ratio = max_frozen_field_ratio
  end if
  if (present(electrostatic_restart_state)) then
    call snapshot%export_restart_state(outer_coupler%last_outer_update_batch, electrostatic_restart_state)
  end if

  if (allocated(potential_buf)) deallocate (potential_buf)
  if (allocated(escaped_boundary_flag)) deallocate (escaped_boundary_flag)
  if (allocated(absorbed_flag)) deallocate (absorbed_flag)
  deallocate ( &
    dq_thread, dq, photo_emission_dq, interface_outward_thread, interface_returned_thread, &
    interface_tau_max_thread, interface_frozen_ratio_max_thread &
    )
  end procedure run_absorption_insulator

  !> 1バッチ分の粒子群と作業バッファを初期化する。
  module procedure prepare_batch_state
  batch_idx = stats%batches + 1_i32
  if (present(inject_state)) then
    call init_particle_batch_from_config( &
      app, batch_idx, pcls_batch, inject_state, mesh=mesh, photo_emission_dq=photo_emission_dq, &
      mpi=mpi, collision_failure_status=collision_failure_status, &
      collision_failure_species=collision_failure_species, collision_failure_ray=collision_failure_ray, &
      collision_failure_bounce=collision_failure_bounce &
      )
  else
    call init_particle_batch_from_config( &
      app, batch_idx, pcls_batch, mesh=mesh, photo_emission_dq=photo_emission_dq, mpi=mpi, &
      collision_failure_status=collision_failure_status, collision_failure_species=collision_failure_species, &
      collision_failure_ray=collision_failure_ray, collision_failure_bounce=collision_failure_bounce &
      )
  end if
  if (collision_failure_status /= collision_query_ok) return
  if (allocated(escaped_boundary_flag)) then
    if (size(escaped_boundary_flag) < pcls_batch%n) then
      deallocate (escaped_boundary_flag)
      allocate (escaped_boundary_flag(pcls_batch%n))
    end if
  else
    allocate (escaped_boundary_flag(pcls_batch%n))
  end if
  if (allocated(absorbed_flag)) then
    if (size(absorbed_flag) < pcls_batch%n) then
      deallocate (absorbed_flag)
      allocate (absorbed_flag(pcls_batch%n))
    end if
  else
    allocate (absorbed_flag(pcls_batch%n))
  end if
  escaped_boundary_flag(:pcls_batch%n) = .false.
  absorbed_flag(:pcls_batch%n) = .false.
  dq_thread = 0.0d0
  end procedure prepare_batch_state

  !> 粒子を時間発展させ、衝突時の堆積電荷をスレッド別に集計する。
  module procedure process_particle_batch
  integer(i32) :: i, step, tid, nth, warn_stride, collision_status
  real(dp) :: x0(3), v0(3), x1(3), v1(3), qdep
  type(hit_info) :: hit
  type(particle_step_result) :: step_result
  type(particle_step_result) :: remainder_result
  type(interface_particle_outcome_type) :: interface_outcome
  logical :: has_warn_stride, collision_failed, candidate_inside, used_event_resolver, interface_active
  integer(i32) :: species_idx

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
  interface_active = trim(lower_ascii(app%coupling%particle_transfer_mode)) == 'electrostatic_1d_instant_return'

  !$omp parallel default(none) &
  !$omp shared(mesh,pcls_batch,app,snapshot,dq_thread,bfield,escaped_boundary_flag,absorbed_flag,nth,interface_active) &
  !$omp shared(interface_outward_thread,interface_returned_thread) &
  !$omp shared(interface_tau_max_thread,interface_frozen_ratio_max_thread) &
  !$omp shared(warn_stride,batch_idx,mpi_rank,collision_failure_status,collision_failure_particle,collision_failure_step) &
  !$omp shared(collision_failure_x,collision_failure_v) &
  !$omp private(i,step,x0,v0,x1,v1,hit,step_result,remainder_result,interface_outcome,tid,qdep,species_idx) &
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
      call find_first_hit(mesh, x0, x1, hit, sim=app%sim, status=collision_status)
      candidate_inside = .not. app%sim%use_box .or. &
                         (all(x1 > app%sim%box_min) .and. all(x1 < app%sim%box_max))
      used_event_resolver = .false.
      if (collision_status /= collision_query_ok) then
        if (.not. candidate_inside) then
          call resolve_particle_boundary_candidate( &
            mesh, app%sim, snapshot, bfield, x0, v0, pcls_batch%q(i), pcls_batch%m(i), app%sim%dt, x1, v1, &
            result=step_result, defer_z_high_interface=interface_active &
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
          hit=hit, result=step_result, defer_z_high_interface=interface_active &
          )
        used_event_resolver = .true.
      end if
      if (used_event_resolver) then
        if (step_result%interface_crossing%has_crossing) then
          species_idx = pcls_batch%species_id(i)
          qdep = pcls_batch%q(i)*pcls_batch%w(i)
          interface_outward_thread(species_idx, tid) = interface_outward_thread(species_idx, tid) + qdep
          call map_outer_particle_linear_debye( &
            snapshot%outer, app%sim%box_min, app%sim%box_max, pcls_batch%q(i), pcls_batch%m(i), &
            step_result%interface_crossing, app%coupling%field_evolution_timescale, &
            app%coupling%max_frozen_field_ratio, app%coupling%outer_queue_enabled, interface_outcome &
            )
          select case (interface_outcome%kind)
          case (interface_outcome_returned_local)
            interface_tau_max_thread(tid) = max(interface_tau_max_thread(tid), interface_outcome%outer_flight_time)
            interface_frozen_ratio_max_thread(tid) = &
              max(interface_frozen_ratio_max_thread(tid), interface_outcome%frozen_field_ratio)
            interface_returned_thread(species_idx, tid) = interface_returned_thread(species_idx, tid) + qdep
            if (step_result%interface_crossing%dt_remaining > 0.0_dp) then
              call advance_particle_step( &
                mesh, app%sim, snapshot, bfield, interface_outcome%position, interface_outcome%velocity, &
                pcls_batch%q(i), pcls_batch%m(i), step_result%interface_crossing%dt_remaining, remainder_result &
                )
              step_result = remainder_result
            else
              step_result%x = interface_outcome%position
              step_result%v = interface_outcome%velocity
              step_result%interface_crossing%has_crossing = .false.
            end if
          case (interface_outcome_escaped_to_infinity)
            step_result%x = interface_outcome%position
            step_result%v = interface_outcome%velocity
            step_result%escaped_boundary = .true.
            step_result%interface_crossing%has_crossing = .false.
          case default
            step_result%status = particle_step_invalid_boundary
          end select
        end if
        if (step_result%status /= collision_query_ok) then
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
    case (particle_step_invalid_boundary)
      failure_name = 'invalid_boundary'
    case (particle_step_multiple_box_events)
      failure_name = 'multiple_box_events'
    case (particle_step_unsupported_barrier_corner)
      failure_name = 'unsupported_barrier_corner'
    case default
      failure_name = 'unknown'
    end select
    write (failure_message, '(a,i0,a,i0,a,i0,a,i0,a,a,a,i0,a,es13.5,a,3es13.5,a,3es13.5)') &
      'particle step failed: batch=', batch_idx, ' particle=', failure_particle, &
      ' step=', failure_step, ' rank=', failure_rank, ' status=', trim(failure_name), &
      ' code=', failure_status, ' dt=', dt, ' x=', failure_x, ' v=', failure_v
    error stop trim(failure_message)
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
    case default
      failure_name = 'unknown'
    end select
    write (failure_message, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,a,a,i0)') &
      'photo collision query incomplete: batch=', batch_idx, ' rank=', failure_rank, &
      ' species=', failure_species, ' ray=', failure_ray, ' bounce=', failure_bounce, &
      ' status=', trim(failure_name), ' code=', failure_status
    error stop trim(failure_message)
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
  real(dp), allocatable :: q_before(:)

  allocate (q_before(mesh%nelem))
  q_before = mesh%q_elem
  dq = sum(dq_thread, dim=2) + photo_emission_dq
  call mpi_allreduce_sum_real_dp_array(mpi, dq)
  mesh%q_elem = mesh%q_elem + dq
  call apply_surface_model_charge_relaxation(mesh, softening, external_e, field_bc_mode=field_bc_mode)
  dq = mesh%q_elem - q_before
  norm_dq = sqrt(sum(dq*dq))
  norm_q = sqrt(sum(mesh%q_elem*mesh%q_elem))
  rel = norm_dq/max(norm_q, q_floor)
  end procedure commit_batch_charge

  !> batch 初期粒子を remote injection と surface emission に分類する。
  subroutine record_batch_initial_charge(app, pcls_batch, ledger)
    type(app_config), intent(in) :: app
    type(particles_soa), intent(in) :: pcls_batch
    type(charge_ledger_type), intent(inout) :: ledger
    integer(i32) :: i, species_idx
    real(dp) :: macro_charge
    character(len=32) :: source_mode

    do i = 1, pcls_batch%n
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
  subroutine record_batch_outcome_charge(pcls_batch, escaped_boundary_flag, absorbed_flag, ledger)
    type(particles_soa), intent(in) :: pcls_batch
    logical, intent(in) :: escaped_boundary_flag(:), absorbed_flag(:)
    type(charge_ledger_type), intent(inout) :: ledger
    integer(i32) :: i, species_idx
    real(dp) :: macro_charge

    do i = 1, pcls_batch%n
      species_idx = pcls_batch%species_id(i)
      if (species_idx < 1_i32 .or. species_idx > ledger%nspecies) then
        error stop 'particle species index is invalid while recording charge ledger output.'
      end if
      macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
      if (absorbed_flag(i)) then
        ledger%absorbed_on_surface(species_idx) = ledger%absorbed_on_surface(species_idx) + macro_charge
        ledger%absorbed_count(species_idx) = ledger%absorbed_count(species_idx) + 1_i64
      else if (escaped_boundary_flag(i)) then
        ledger%escaped_to_infinity(species_idx) = ledger%escaped_to_infinity(species_idx) + macro_charge
        ledger%escaped_count(species_idx) = ledger%escaped_count(species_idx) + 1_i64
      else if (pcls_batch%alive(i)) then
        ledger%discarded_unresolved(species_idx) = ledger%discarded_unresolved(species_idx) + macro_charge
        ledger%discarded_unresolved_count(species_idx) = ledger%discarded_unresolved_count(species_idx) + 1_i64
      end if
    end do
  end subroutine record_batch_outcome_charge

  !> species 別 flux/count だけを MPI-global 値へ集約する。stock は全 rank で同じ mesh state から得る。
  subroutine reduce_charge_ledger_fluxes(ledger, mpi)
    type(charge_ledger_type), intent(inout) :: ledger
    type(mpi_context), intent(in) :: mpi
    real(dp), allocatable :: charge_values(:)
    integer(i64), allocatable :: count_values(:)
    integer :: n

    n = int(ledger%nspecies)
    allocate (charge_values(7*n), count_values(5*n))
    charge_values(1:n) = ledger%injected_from_remote
    charge_values(n + 1:2*n) = ledger%emitted_from_surface
    charge_values(2*n + 1:3*n) = ledger%absorbed_on_surface
    charge_values(3*n + 1:4*n) = ledger%escaped_to_infinity
    charge_values(4*n + 1:5*n) = ledger%discarded_unresolved
    charge_values(5*n + 1:6*n) = ledger%interface_outward_gross
    charge_values(6*n + 1:7*n) = ledger%interface_returned_gross
    call mpi_allreduce_sum_real_dp_array(mpi, charge_values)
    ledger%injected_from_remote = charge_values(1:n)
    ledger%emitted_from_surface = charge_values(n + 1:2*n)
    ledger%absorbed_on_surface = charge_values(2*n + 1:3*n)
    ledger%escaped_to_infinity = charge_values(3*n + 1:4*n)
    ledger%discarded_unresolved = charge_values(4*n + 1:5*n)
    ledger%interface_outward_gross = charge_values(5*n + 1:6*n)
    ledger%interface_returned_gross = charge_values(6*n + 1:7*n)

    count_values(1:n) = ledger%injected_count
    count_values(n + 1:2*n) = ledger%emitted_count
    count_values(2*n + 1:3*n) = ledger%absorbed_count
    count_values(3*n + 1:4*n) = ledger%escaped_count
    count_values(4*n + 1:5*n) = ledger%discarded_unresolved_count
    call mpi_allreduce_sum_i64_array(mpi, count_values)
    ledger%injected_count = count_values(1:n)
    ledger%emitted_count = count_values(n + 1:2*n)
    ledger%absorbed_count = count_values(2*n + 1:3*n)
    ledger%escaped_count = count_values(3*n + 1:4*n)
    ledger%discarded_unresolved_count = count_values(4*n + 1:5*n)
  end subroutine reduce_charge_ledger_fluxes

end submodule bem_simulator_loop

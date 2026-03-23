!> `bem_simulator` の主ループと粒子処理計算を実装する submodule。
submodule(bem_simulator) bem_simulator_loop
  use bem_output_writer, only: compute_mesh_potential_fmm_core
  use bem_performance_profile, only: perf_region_batch_total, perf_region_begin, perf_region_commit_charge, &
                                     perf_region_count_outcomes, perf_region_end, perf_region_field_refresh, &
                                     perf_region_field_solver_init, perf_region_history_write, perf_region_mpi_reduce, &
                                     perf_region_particle_batch, perf_region_prepare_batch, perf_region_simulation_total, &
                                     perf_region_stats_update
  implicit none
contains

  !> 吸着モデルのバッチループを実行し、電荷更新と統計集計を進める。
  module procedure run_absorption_insulator
  integer(i32) :: batch_idx, final_batch_idx, local_batch_idx, nth, hist_stride
  integer :: hist_unit
  logical :: history_enabled
  real(dp), allocatable :: dq_thread(:, :), dq(:), photo_emission_dq(:)
  logical, allocatable :: escaped_boundary_flag(:), absorbed_flag(:)
  integer(i32) :: batch_counts(5)
  real(dp) :: bfield(3), rel, t0, sim_t0, batch_t0
  type(particles_soa) :: pcls_batch
  type(mpi_context) :: mpi_ctx
  type(field_solver_type) :: field_solver = field_solver_type()

  stats = sim_stats()
  if (present(initial_stats)) stats = initial_stats
  mpi_ctx = mpi_context()
  if (present(mpi)) mpi_ctx = mpi

  nth = 1
!$ nth = max(1, omp_get_max_threads())
  allocate (dq_thread(mesh%nelem, nth), dq(mesh%nelem), photo_emission_dq(mesh%nelem))

  history_enabled = present(history_unit)
  hist_unit = 0
  if (history_enabled) hist_unit = history_unit
  hist_stride = 1_i32
  if (present(history_stride)) hist_stride = max(1_i32, history_stride)
  bfield = app%sim%b0
  final_batch_idx = stats%batches + app%sim%batch_count
  call perf_region_begin(perf_region_simulation_total, sim_t0)
  call perf_region_begin(perf_region_field_solver_init, t0)
  call field_solver%init(mesh, app%sim)
  call perf_region_end(perf_region_field_solver_init, t0)

  do local_batch_idx = 1, app%sim%batch_count
    call perf_region_begin(perf_region_batch_total, batch_t0)

    call perf_region_begin(perf_region_prepare_batch, t0)
    call prepare_batch_state( &
      mesh, app, stats, local_batch_idx, batch_idx, dq_thread, pcls_batch, escaped_boundary_flag, absorbed_flag, &
      photo_emission_dq, mpi_ctx, inject_state &
      )
    call perf_region_end(perf_region_prepare_batch, t0)

    call perf_region_begin(perf_region_field_refresh, t0)
    call field_solver%refresh(mesh)
    call perf_region_end(perf_region_field_refresh, t0)

    call perf_region_begin(perf_region_particle_batch, t0)
    call process_particle_batch( &
      mesh, app, field_solver, pcls_batch, dq_thread, escaped_boundary_flag, absorbed_flag, bfield &
      )
    call perf_region_end(perf_region_particle_batch, t0)

    call perf_region_begin(perf_region_commit_charge, t0)
    call commit_batch_charge(mesh, app%sim%q_floor, dq_thread, photo_emission_dq, dq, rel, mpi_ctx)
    call perf_region_end(perf_region_commit_charge, t0)

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
    end if
    call perf_region_end(perf_region_history_write, t0)

    deallocate (escaped_boundary_flag, absorbed_flag)

    call perf_region_end(perf_region_batch_total, batch_t0)
  end do
  call perf_region_end(perf_region_simulation_total, sim_t0)

  if (present(mesh_potential_v)) then
    if (mpi_is_root(mpi_ctx) .and. trim(field_solver%mode) == 'fmm' .and. field_solver%fmm_core_ready) then
      call perf_region_begin(perf_region_field_refresh, t0)
      call field_solver%refresh(mesh)
      call perf_region_end(perf_region_field_refresh, t0)
      allocate (mesh_potential_v(mesh%nelem))
      call compute_mesh_potential_fmm_core(mesh, field_solver%fmm_core_plan, field_solver%fmm_core_state, mesh_potential_v)
    end if
  end if

  deallocate (dq_thread, dq, photo_emission_dq)
  end procedure run_absorption_insulator

  !> 1バッチ分の粒子群と作業バッファを初期化する。
  module procedure prepare_batch_state
  batch_idx = stats%batches + 1_i32
  if (present(inject_state)) then
    call init_particle_batch_from_config( &
      app, local_batch_idx, pcls_batch, inject_state, mesh=mesh, photo_emission_dq=photo_emission_dq, &
      mpi=mpi &
      )
  else
    call init_particle_batch_from_config( &
      app, local_batch_idx, pcls_batch, mesh=mesh, photo_emission_dq=photo_emission_dq, mpi=mpi &
      )
  end if
  allocate (escaped_boundary_flag(pcls_batch%n), absorbed_flag(pcls_batch%n))
  escaped_boundary_flag = .false.
  absorbed_flag = .false.
  dq_thread = 0.0d0
  end procedure prepare_batch_state

  !> 粒子を時間発展させ、衝突時の堆積電荷をスレッド別に集計する。
  module procedure process_particle_batch
  integer(i32) :: i, step, tid, nth
  real(dp) :: x0(3), v0(3), x1(3), v1(3), e(3), qdep
  type(hit_info) :: hit
  logical :: escaped_by_boundary

  nth = size(dq_thread, 2)

  !$omp parallel default(none) &
  !$omp shared(mesh,pcls_batch,app,field_solver,dq_thread,bfield,escaped_boundary_flag,absorbed_flag,nth) &
  !$omp private(i,step,x0,v0,x1,v1,e,hit,tid,qdep,escaped_by_boundary)
  ! スレッドごとに dq_thread(:, tid) を使って原子的更新なしで電荷を集める。
  tid = 1
!$ tid = omp_get_thread_num() + 1
  !$omp do schedule(static)
  do i = 1, pcls_batch%n
    if (.not. pcls_batch%alive(i)) cycle
    do step = 1, app%sim%max_step
      x0 = pcls_batch%x(:, i)
      v0 = pcls_batch%v(:, i)
      call field_solver%eval_e(mesh, x0, e)
      call boris_push( &
        x0, v0, pcls_batch%q(i), pcls_batch%m(i), app%sim%dt, e, bfield, x1, v1 &
        )
      call find_first_hit(mesh, x0, x1, hit, sim=app%sim)
      if (hit%has_hit) then
        qdep = pcls_batch%q(i)*pcls_batch%w(i)
        dq_thread(hit%elem_idx, tid) = dq_thread(hit%elem_idx, tid) + qdep
        pcls_batch%alive(i) = .false.
        absorbed_flag(i) = .true.
        exit
      end if

      call apply_box_boundary(app%sim, x1, v1, pcls_batch%alive(i), escaped_by_boundary)
      if (escaped_by_boundary) escaped_boundary_flag(i) = .true.
      if (.not. pcls_batch%alive(i)) exit

      pcls_batch%x(:, i) = x1
      pcls_batch%v(:, i) = v1
    end do
  end do
  !$omp end do
  !$omp end parallel
  end procedure process_particle_batch

  !> スレッド別電荷差分を合算してメッシュへ反映し、相対変化量を返す。
  module procedure commit_batch_charge
  real(dp) :: norm_dq, norm_q

  dq = sum(dq_thread, dim=2) + photo_emission_dq
  call mpi_allreduce_sum_real_dp_array(mpi, dq)
  mesh%q_elem = mesh%q_elem + dq
  norm_dq = sqrt(sum(dq*dq))
  norm_q = sqrt(sum(mesh%q_elem*mesh%q_elem))
  rel = norm_dq/max(norm_q, q_floor)
  end procedure commit_batch_charge

end submodule bem_simulator_loop

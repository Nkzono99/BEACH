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
  integer :: hist_unit
  logical :: history_enabled
  logical :: potential_history_enabled
  integer :: pot_hist_unit
  real(dp), allocatable :: potential_buf(:)
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
  call field_solver%init(mesh, app%sim)
  call perf_region_end(perf_region_field_solver_init, t0)

  do local_batch_idx = 1, batch_count_this_run
    call perf_region_begin(perf_region_batch_total, batch_t0)

    call perf_region_begin(perf_region_prepare_batch, t0)
    call prepare_batch_state( &
      mesh, app, stats, batch_idx, dq_thread, pcls_batch, escaped_boundary_flag, absorbed_flag, &
      photo_emission_dq, mpi_ctx, inject_state &
      )
    call perf_region_end(perf_region_prepare_batch, t0)

    call perf_region_begin(perf_region_field_refresh, t0)
    call field_solver%refresh(mesh)
    call perf_region_end(perf_region_field_refresh, t0)

    call perf_region_begin(perf_region_particle_batch, t0)
    call process_particle_batch( &
      mesh, app, field_solver, pcls_batch, dq_thread, escaped_boundary_flag, absorbed_flag, bfield, batch_idx, mpi_ctx%rank &
      )
    call perf_region_end(perf_region_particle_batch, t0)

    call perf_region_begin(perf_region_commit_charge, t0)
    call commit_batch_charge( &
      mesh, app%sim%q_floor, app%sim%softening, app%sim%e0, app%sim%field_bc_mode, &
      dq_thread, photo_emission_dq, dq, rel, mpi_ctx &
      )
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
      if (potential_history_enabled) then
        call field_solver%refresh(mesh)
        call maybe_write_potential_history_snapshot( &
          potential_history_enabled, pot_hist_unit, hist_stride, stats, &
          field_solver, mesh, app%sim, potential_buf &
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
      call field_solver%refresh(mesh)
      call perf_region_end(perf_region_field_refresh, t0)
      allocate (mesh_potential_v(mesh%nelem))
      call field_solver%compute_mesh_potential(mesh, app%sim, mesh_potential_v)
    end if
  end if

  if (allocated(potential_buf)) deallocate (potential_buf)
  if (allocated(escaped_boundary_flag)) deallocate (escaped_boundary_flag)
  if (allocated(absorbed_flag)) deallocate (absorbed_flag)
  deallocate (dq_thread, dq, photo_emission_dq)
  end procedure run_absorption_insulator

  !> 1バッチ分の粒子群と作業バッファを初期化する。
  module procedure prepare_batch_state
  batch_idx = stats%batches + 1_i32
  if (present(inject_state)) then
    call init_particle_batch_from_config( &
      app, batch_idx, pcls_batch, inject_state, mesh=mesh, photo_emission_dq=photo_emission_dq, &
      mpi=mpi &
      )
  else
    call init_particle_batch_from_config( &
      app, batch_idx, pcls_batch, mesh=mesh, photo_emission_dq=photo_emission_dq, mpi=mpi &
      )
  end if
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
  integer(i32) :: i, step, tid, nth, warn_stride
  integer(i32) :: boundary_axis, collision_status
  integer(i32) :: collision_failure_status, collision_failure_particle, collision_failure_step
  real(dp) :: x0(3), v0(3), x1(3), v1(3), e(3), qdep
  real(dp) :: boundary_probe(3), phi_boundary
  type(hit_info) :: hit
  logical :: escaped_by_boundary, has_warn_stride, has_boundary_probe, boundary_high_side
  character(len=16) :: collision_failure_name
  character(len=256) :: collision_failure_message

  nth = size(dq_thread, 2)
  call read_env_i32_local('BEACH_WARN_LONG_PARTICLE_STEPS', warn_stride, has_warn_stride)
  if (.not. has_warn_stride) warn_stride = 0_i32
  collision_failure_status = collision_query_ok
  collision_failure_particle = huge(0_i32)
  collision_failure_step = 0_i32

  !$omp parallel default(none) &
  !$omp shared(mesh,pcls_batch,app,field_solver,dq_thread,bfield,escaped_boundary_flag,absorbed_flag,nth) &
  !$omp shared(warn_stride,batch_idx,mpi_rank,collision_failure_status,collision_failure_particle,collision_failure_step) &
  !$omp private(i,step,x0,v0,x1,v1,e,boundary_axis,boundary_probe,phi_boundary,hit,tid,qdep,escaped_by_boundary) &
  !$omp private(has_boundary_probe,boundary_high_side,collision_status)
  ! スレッドごとに dq_thread(:, tid) を使って原子的更新なしで電荷を集める。
  tid = 1
!$ tid = omp_get_thread_num() + 1
  !$omp do schedule(dynamic, 1)
  do i = 1, pcls_batch%n
    if (.not. pcls_batch%alive(i)) cycle
    do step = 1, app%sim%max_step
      x0 = pcls_batch%x(:, i)
      v0 = pcls_batch%v(:, i)
      call field_solver%eval_e(mesh, x0, e)
      e = e + app%sim%e0
      call boris_push( &
        x0, v0, pcls_batch%q(i), pcls_batch%m(i), app%sim%dt, e, bfield, x1, v1 &
        )
      call find_first_hit(mesh, x0, x1, hit, sim=app%sim, status=collision_status)
      if (collision_status /= collision_query_ok) then
        !$omp critical (beach_collision_query_failure)
        if (collision_failure_status == collision_query_ok .or. i < collision_failure_particle .or. &
            (i == collision_failure_particle .and. step < collision_failure_step)) then
          collision_failure_status = collision_status
          collision_failure_particle = i
          collision_failure_step = step
        end if
        !$omp end critical (beach_collision_query_failure)
        exit
      end if
      if (hit%has_hit) then
        qdep = pcls_batch%q(i)*pcls_batch%w(i)
        dq_thread(hit%elem_idx, tid) = dq_thread(hit%elem_idx, tid) + qdep
        pcls_batch%alive(i) = .false.
        absorbed_flag(i) = .true.
        exit
      end if

      if (trim(app%sim%open_boundary_model) == 'potential_barrier') then
        call find_open_boundary_probe(app%sim, x0, x1, has_boundary_probe, boundary_probe, boundary_axis, boundary_high_side)
        if (has_boundary_probe) then
          call field_solver%eval_potential(mesh, app%sim, boundary_probe, phi_boundary)
          call apply_box_boundary( &
            app%sim, x1, v1, pcls_batch%alive(i), escaped_by_boundary, &
            q_particle=pcls_batch%q(i), m_particle=pcls_batch%m(i), phi_boundary=phi_boundary, &
            potential_axis=boundary_axis, potential_high_side=boundary_high_side &
            )
        else
          call apply_box_boundary(app%sim, x1, v1, pcls_batch%alive(i), escaped_by_boundary)
        end if
      else
        call apply_box_boundary(app%sim, x1, v1, pcls_batch%alive(i), escaped_by_boundary)
      end if
      if (escaped_by_boundary) escaped_boundary_flag(i) = .true.
      if (.not. pcls_batch%alive(i)) exit

      pcls_batch%x(:, i) = x1
      pcls_batch%v(:, i) = v1
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
    if (warn_stride > 0_i32 .and. pcls_batch%alive(i)) then
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

  if (collision_failure_status /= collision_query_ok) then
    select case (collision_failure_status)
    case (collision_query_image_limit)
      collision_failure_name = 'image_limit'
    case (collision_query_index_range)
      collision_failure_name = 'index_range'
    case default
      collision_failure_name = 'unknown'
    end select
    write (collision_failure_message, '(a,i0,a,i0,a,i0,a,i0,a,a,a,i0)') &
      'collision query incomplete: batch=', batch_idx, ' particle=', collision_failure_particle, &
      ' step=', collision_failure_step, ' rank=', mpi_rank, ' status=', trim(collision_failure_name), &
      ' code=', collision_failure_status
    error stop trim(collision_failure_message)
  end if
  end procedure process_particle_batch

  !> x0->x1 が最初に越える open box 面上の評価点を返す。
  subroutine find_open_boundary_probe(sim, x0, x1, found, probe, probe_axis, high_side)
    type(sim_config), intent(in) :: sim
    real(dp), intent(in) :: x0(3), x1(3)
    logical, intent(out) :: found
    real(dp), intent(out) :: probe(3)
    integer(i32), intent(out) :: probe_axis
    logical, intent(out) :: high_side

    integer(i32) :: axis
    real(dp) :: denom, t, best_t

    found = .false.
    probe = min(sim%box_max, max(sim%box_min, x1))
    probe_axis = 0_i32
    high_side = .false.
    if (.not. sim%use_box) return

    best_t = huge(1.0d0)
    do axis = 1_i32, 3_i32
      denom = x1(axis) - x0(axis)
      if (denom < 0.0d0 .and. sim%bc_low(axis) == bc_open .and. x1(axis) < sim%box_min(axis)) then
        t = (sim%box_min(axis) - x0(axis))/denom
        call maybe_update_boundary_probe(sim, x0, x1, axis, .false., sim%box_min(axis), t, best_t, found, probe, &
                                         probe_axis, high_side)
      else if (denom > 0.0d0 .and. sim%bc_high(axis) == bc_open .and. x1(axis) > sim%box_max(axis)) then
        t = (sim%box_max(axis) - x0(axis))/denom
        call maybe_update_boundary_probe(sim, x0, x1, axis, .true., sim%box_max(axis), t, best_t, found, probe, &
                                         probe_axis, high_side)
      end if
    end do
  end subroutine find_open_boundary_probe

  !> 最小の有効 crossing parameter を持つ境界評価点へ更新する。
  subroutine maybe_update_boundary_probe( &
    sim, x0, x1, axis, candidate_high_side, boundary_value, t, best_t, found, probe, probe_axis, high_side &
    )
    type(sim_config), intent(in) :: sim
    real(dp), intent(in) :: x0(3), x1(3), boundary_value, t
    integer(i32), intent(in) :: axis
    logical, intent(in) :: candidate_high_side
    real(dp), intent(inout) :: best_t
    logical, intent(inout) :: found
    real(dp), intent(inout) :: probe(3)
    integer(i32), intent(inout) :: probe_axis
    logical, intent(inout) :: high_side

    if (t < 0.0d0 .or. t > 1.0d0 .or. t >= best_t) return
    best_t = t
    probe = x0 + t*(x1 - x0)
    probe(axis) = boundary_value
    probe_axis = axis
    high_side = candidate_high_side
    found = .true.
  end subroutine maybe_update_boundary_probe

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

end submodule bem_simulator_loop

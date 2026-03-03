!> 吸着(insulator)モデルのメインループを実行し、電荷堆積と統計更新を行うモジュール。
module bem_simulator
  use omp_lib
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_config, sim_stats, mesh_type, particles_soa, hit_info
  use bem_field, only: electric_field_at
  use bem_pusher, only: boris_push
  use bem_collision, only: find_first_hit
  implicit none
contains

  !> 粒子をバッチ処理し、衝突時は要素へ電荷堆積、非衝突時は脱出として統計を更新する。
  !! @param[inout] mesh 入出力引数。
  !! @param[in] cfg 入力引数。
  !! @param[inout] pcls 入出力引数。
  !! @param[out] stats 出力引数。
  subroutine run_absorption_insulator(mesh, cfg, pcls, stats, history_unit, history_stride)
    type(mesh_type), intent(inout) :: mesh
    type(sim_config), intent(in) :: cfg
    type(particles_soa), intent(inout) :: pcls
    type(sim_stats), intent(out) :: stats
    integer, intent(in), optional :: history_unit
    integer(i32), intent(in), optional :: history_stride

    integer(i32) :: batch_start, batch_end, i, step, tid, nth, b, hist_stride
    integer :: hist_unit
    logical :: do_write_history
    real(dp), allocatable :: dq_thread(:, :), dq(:)
    real(dp) :: x0(3), v0(3), x1(3), v1(3), e(3), bfield(3), rel, norm_dq, norm_q, qdep
    type(hit_info) :: hit

    stats = sim_stats()
    nth = max(1, omp_get_max_threads())
    allocate(dq_thread(mesh%nelem, nth), dq(mesh%nelem))
    hist_unit = -1
    if (present(history_unit)) hist_unit = history_unit
    hist_stride = 1_i32
    if (present(history_stride)) hist_stride = max(1_i32, history_stride)
    bfield = cfg%b0

    batch_start = 1
    do while (batch_start <= pcls%n)
      batch_end = min(pcls%n, batch_start + cfg%npcls_per_step - 1)
      dq_thread = 0.0d0

      !$omp parallel default(none) &
      !$omp shared(mesh,pcls,cfg,batch_start,batch_end,dq_thread,bfield) &
      !$omp private(i,step,x0,v0,x1,v1,e,hit,tid,qdep)
      tid = omp_get_thread_num() + 1
      !$omp do schedule(static)
      do i = batch_start, batch_end
        if (.not. pcls%alive(i)) cycle
        do step = 1, cfg%max_step
          x0 = pcls%x(:, i)
          v0 = pcls%v(:, i)
          call electric_field_at(mesh, x0, cfg%softening, e)
          call boris_push(x0, v0, pcls%q(i), pcls%m(i), cfg%dt, e, bfield, x1, v1)
          call find_first_hit(mesh, x0, x1, hit)
          if (hit%has_hit) then
            qdep = pcls%q(i) * pcls%w(i)
            dq_thread(hit%elem_idx, tid) = dq_thread(hit%elem_idx, tid) + qdep
            pcls%alive(i) = .false.
            exit
          else
            pcls%x(:, i) = x1
            pcls%v(:, i) = v1
          end if
        end do
      end do
      !$omp end do
      !$omp end parallel

      dq = sum(dq_thread, dim=2)
      mesh%q_elem = mesh%q_elem + dq
      norm_dq = sqrt(sum(dq * dq))
      norm_q = sqrt(sum(mesh%q_elem * mesh%q_elem))
      rel = norm_dq / max(norm_q, cfg%q_floor)

      stats%batches = stats%batches + 1
      stats%last_rel_change = rel
      stats%processed_particles = stats%processed_particles + (batch_end - batch_start + 1)
      do_write_history = (hist_unit > 0) .and. (mod(stats%batches - 1, hist_stride) == 0)
      if (do_write_history) call write_history_snapshot(hist_unit, stats%batches, stats%processed_particles, rel, mesh%q_elem)

      do b = batch_start, batch_end
        if (pcls%alive(b)) then
          stats%escaped = stats%escaped + 1
        else
          stats%absorbed = stats%absorbed + 1
        end if
      end do

      if (rel < cfg%tol_rel) exit
      batch_start = batch_end + 1
    end do

  end subroutine run_absorption_insulator

  subroutine write_history_snapshot(unit_id, batch_idx, processed_particles, rel_change, q_elem)
    integer, intent(in) :: unit_id
    integer(i32), intent(in) :: batch_idx
    integer(i32), intent(in) :: processed_particles
    real(dp), intent(in) :: rel_change
    real(dp), intent(in) :: q_elem(:)
    integer(i32) :: elem_idx

    do elem_idx = 1, size(q_elem)
      write(unit_id, '(i0,a,i0,a,es24.16,a,i0,a,es24.16)') batch_idx, ',', processed_particles, ',', &
        rel_change, ',', elem_idx, ',', q_elem(elem_idx)
    end do
  end subroutine write_history_snapshot

end module bem_simulator

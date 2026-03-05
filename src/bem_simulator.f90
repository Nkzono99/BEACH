!> 吸着(insulator)モデルのメインループを実行し、電荷堆積と統計更新を行うモジュール。
module bem_simulator
  !$ use omp_lib
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_stats, mesh_type, particles_soa, hit_info, injection_state
  use bem_app_config, only: app_config, init_particle_batch_from_config
  use bem_field, only: electric_field_at
  use bem_pusher, only: boris_push
  use bem_collision, only: find_first_hit
  use bem_boundary, only: apply_box_boundary
  implicit none
  private

  public :: run_absorption_insulator

contains

  !> 粒子をバッチ処理し、衝突時は要素へ電荷堆積、非衝突時は脱出として統計を更新する。
  !! @param[inout] mesh 入力メッシュ。各バッチ後に `q_elem` へ堆積電荷を加算して更新。
  !! @param[in] app 時間刻み・バッチ回数・監視値・境界条件を含む実行設定。
  !! @param[out] stats 吸着数・脱出数・最終相対変化率などの集計統計。
  !! @param[in] history_unit 履歴CSVの出力先ユニット番号（省略時は履歴出力なし）。
  !! @param[in] history_stride 履歴出力のバッチ間隔（省略時1）。
  !! @param[in] initial_stats 再開時に引き継ぐ既存統計（省略時はゼロ初期化）。
  !! @param[inout] inject_state `reservoir_face` 注入で使う種別ごとの残差状態（省略可）。
  subroutine run_absorption_insulator(mesh, app, stats, history_unit, history_stride, initial_stats, inject_state)
    type(mesh_type), intent(inout) :: mesh
    type(app_config), intent(in) :: app
    type(sim_stats), intent(out) :: stats
    integer, intent(in), optional :: history_unit
    integer(i32), intent(in), optional :: history_stride
    type(sim_stats), intent(in), optional :: initial_stats
    type(injection_state), intent(inout), optional :: inject_state

    integer(i32) :: batch_idx, local_batch_idx, nth, hist_stride
    integer :: hist_unit
    logical :: history_enabled
    real(dp), allocatable :: dq_thread(:, :), dq(:)
    logical, allocatable :: escaped_boundary_flag(:), absorbed_flag(:)
    real(dp) :: bfield(3), rel
    type(particles_soa) :: pcls_batch

    stats = sim_stats()
    if (present(initial_stats)) stats = initial_stats

    nth = 1
    !$ nth = max(1, omp_get_max_threads())
    allocate (dq_thread(mesh%nelem, nth), dq(mesh%nelem))

    history_enabled = present(history_unit)
    hist_unit = 0
    if (history_enabled) hist_unit = history_unit
    hist_stride = 1_i32
    if (present(history_stride)) hist_stride = max(1_i32, history_stride)
    bfield = app%sim%b0

    do local_batch_idx = 1, app%sim%batch_count
      call prepare_batch_state( &
        app, stats, local_batch_idx, batch_idx, dq_thread, pcls_batch, escaped_boundary_flag, absorbed_flag, inject_state &
      )
      call process_particle_batch( &
        mesh, app, pcls_batch, dq_thread, escaped_boundary_flag, absorbed_flag, bfield &
      )
      call commit_batch_charge(mesh, app%sim%q_floor, dq_thread, dq, rel)
      call accumulate_batch_stats(stats, pcls_batch, escaped_boundary_flag, absorbed_flag, rel)
      call maybe_write_history_snapshot(history_enabled, hist_unit, hist_stride, stats, rel, mesh%q_elem)
      deallocate (escaped_boundary_flag, absorbed_flag)
    end do

    deallocate (dq_thread, dq)
  end subroutine run_absorption_insulator

  !> 1バッチ分の粒子群と作業配列を初期化する。
  !! @param[in] app 実行設定。
  !! @param[in] stats 現在までの累積統計。
  !! @param[in] local_batch_idx 今回実行するローカルバッチ番号（1始まり）。
  !! @param[out] batch_idx 再開分を加味した累積バッチ番号。
  !! @param[inout] dq_thread スレッド別要素電荷差分バッファ（ゼロ初期化される）。
  !! @param[out] pcls_batch 今バッチで処理する粒子群。
  !! @param[out] escaped_boundary_flag 粒子ごとの境界流出フラグ。
  !! @param[out] absorbed_flag 粒子ごとの吸着フラグ。
  !! @param[inout] inject_state reservoir_face 注入残差（指定時のみ更新）。
  subroutine prepare_batch_state( &
    app, stats, local_batch_idx, batch_idx, dq_thread, pcls_batch, escaped_boundary_flag, absorbed_flag, inject_state &
  )
    type(app_config), intent(in) :: app
    type(sim_stats), intent(in) :: stats
    integer(i32), intent(in) :: local_batch_idx
    integer(i32), intent(out) :: batch_idx
    real(dp), intent(inout) :: dq_thread(:, :)
    type(particles_soa), intent(out) :: pcls_batch
    logical, allocatable, intent(out) :: escaped_boundary_flag(:)
    logical, allocatable, intent(out) :: absorbed_flag(:)
    type(injection_state), intent(inout), optional :: inject_state

    batch_idx = stats%batches + 1_i32
    if (present(inject_state)) then
      call init_particle_batch_from_config(app, local_batch_idx, pcls_batch, inject_state)
    else
      call init_particle_batch_from_config(app, local_batch_idx, pcls_batch)
    end if
    allocate (escaped_boundary_flag(pcls_batch%n), absorbed_flag(pcls_batch%n))
    escaped_boundary_flag = .false.
    absorbed_flag = .false.
    dq_thread = 0.0d0

    print *, "---------- batch", batch_idx, "----------"
  end subroutine prepare_batch_state

  !> 1バッチぶんの粒子を前進させ、スレッド別に堆積電荷を集計する。
  !! @param[in] mesh 衝突判定と電場計算に使うメッシュ（読み取り専用）。
  !! @param[in] app 時間積分・境界条件などの実行設定。
  !! @param[inout] pcls_batch 処理対象粒子群（位置・速度・alive を更新）。
  !! @param[inout] dq_thread スレッド別要素電荷差分バッファ。
  !! @param[inout] escaped_boundary_flag 粒子ごとの境界流出フラグ。
  !! @param[inout] absorbed_flag 粒子ごとの吸着フラグ。
  !! @param[in] bfield 一様外部磁場ベクトル [T]。
  subroutine process_particle_batch( &
    mesh, app, pcls_batch, dq_thread, escaped_boundary_flag, absorbed_flag, bfield &
  )
    type(mesh_type), intent(in) :: mesh
    type(app_config), intent(in) :: app
    type(particles_soa), intent(inout) :: pcls_batch
    real(dp), intent(inout) :: dq_thread(:, :)
    logical, intent(inout) :: escaped_boundary_flag(:)
    logical, intent(inout) :: absorbed_flag(:)
    real(dp), intent(in) :: bfield(3)

    integer(i32) :: i, step, tid
    real(dp) :: x0(3), v0(3), x1(3), v1(3), e(3), qdep
    type(hit_info) :: hit
    logical :: escaped_by_boundary

    !$omp parallel default(none) &
    !$omp shared(mesh,pcls_batch,app,dq_thread,bfield,escaped_boundary_flag,absorbed_flag) &
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
        call electric_field_at(mesh, x0, app%sim%softening, e)
        call boris_push( &
          x0, v0, pcls_batch%q(i), pcls_batch%m(i), app%sim%dt, e, bfield, x1, v1 &
        )
        call apply_box_boundary(app%sim, x1, v1, pcls_batch%alive(i), escaped_by_boundary)
        if (escaped_by_boundary) escaped_boundary_flag(i) = .true.
        if (.not. pcls_batch%alive(i)) exit

        call find_first_hit(mesh, x0, x1, hit)
        if (hit%has_hit) then
          qdep = pcls_batch%q(i) * pcls_batch%w(i)
          dq_thread(hit%elem_idx, tid) = dq_thread(hit%elem_idx, tid) + qdep
          pcls_batch%alive(i) = .false.
          absorbed_flag(i) = .true.
          exit
        end if

        pcls_batch%x(:, i) = x1
        pcls_batch%v(:, i) = v1
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine process_particle_batch

  !> スレッド別に集計した電荷差分をメッシュへ反映し、相対変化量を返す。
  !! @param[inout] mesh 要素電荷 `q_elem` を更新するメッシュ。
  !! @param[in] q_floor 相対変化量分母の下限値。
  !! @param[in] dq_thread スレッド別要素電荷差分 `dq_thread(nelem,nthread)`。
  !! @param[out] dq スレッド合算後の要素電荷差分 `dq(nelem)`。
  !! @param[out] rel 今バッチの相対変化量 `||dq||/max(||q||,q_floor)`。
  subroutine commit_batch_charge(mesh, q_floor, dq_thread, dq, rel)
    type(mesh_type), intent(inout) :: mesh
    real(dp), intent(in) :: q_floor
    real(dp), intent(in) :: dq_thread(:, :)
    real(dp), intent(out) :: dq(:)
    real(dp), intent(out) :: rel
    real(dp) :: norm_dq, norm_q

    dq = sum(dq_thread, dim=2)
    mesh%q_elem = mesh%q_elem + dq
    norm_dq = sqrt(sum(dq * dq))
    norm_q = sqrt(sum(mesh%q_elem * mesh%q_elem))
    rel = norm_dq / max(norm_q, q_floor)
  end subroutine commit_batch_charge

  !> バッチ完了後の統計値を更新する。
  !! @param[inout] stats 累積統計（このバッチ結果で更新）。
  !! @param[in] pcls_batch 今バッチで処理した粒子群。
  !! @param[in] escaped_boundary_flag 粒子ごとの境界流出フラグ。
  !! @param[in] absorbed_flag 粒子ごとの吸着フラグ。
  !! @param[in] rel 今バッチの相対変化量。
  subroutine accumulate_batch_stats(stats, pcls_batch, escaped_boundary_flag, absorbed_flag, rel)
    type(sim_stats), intent(inout) :: stats
    type(particles_soa), intent(in) :: pcls_batch
    logical, intent(in) :: escaped_boundary_flag(:)
    logical, intent(in) :: absorbed_flag(:)
    real(dp), intent(in) :: rel
    integer(i32) :: i

    stats%batches = stats%batches + 1_i32
    stats%last_rel_change = rel
    stats%processed_particles = stats%processed_particles + pcls_batch%n

    do i = 1, pcls_batch%n
      if (absorbed_flag(i)) then
        stats%absorbed = stats%absorbed + 1
      else if (escaped_boundary_flag(i)) then
        stats%escaped = stats%escaped + 1
        stats%escaped_boundary = stats%escaped_boundary + 1
      else if (pcls_batch%alive(i)) then
        stats%escaped = stats%escaped + 1
        stats%survived_max_step = stats%survived_max_step + 1
      end if
    end do
  end subroutine accumulate_batch_stats

  !> 履歴出力が有効で、指定ストライドを満たす場合のみ電荷履歴を書き出す。
  !! ストライド判定は従来どおり 1, 1+stride, ... のバッチ番号で行う。
  !! @param[in] history_enabled 履歴出力の有効フラグ。
  !! @param[in] hist_unit 履歴CSVの出力先ユニット番号。
  !! @param[in] hist_stride 履歴を書き出すバッチ間隔。
  !! @param[in] stats 現在の累積統計。
  !! @param[in] rel 今バッチの相対変化量。
  !! @param[in] q_elem 現在の要素電荷配列 [C]。
  subroutine maybe_write_history_snapshot(history_enabled, hist_unit, hist_stride, stats, rel, q_elem)
    logical, intent(in) :: history_enabled
    integer, intent(in) :: hist_unit
    integer(i32), intent(in) :: hist_stride
    type(sim_stats), intent(in) :: stats
    real(dp), intent(in) :: rel
    real(dp), intent(in) :: q_elem(:)

    if (.not. history_enabled) return
    if (mod(stats%batches - 1_i32, hist_stride) /= 0_i32) return
    call write_history_snapshot(hist_unit, stats%batches, stats%processed_particles, rel, q_elem)
  end subroutine maybe_write_history_snapshot

  !> 現時点の要素電荷を CSV 行群として書き出す。
  !! @param[in] unit_id 出力先ユニット番号。
  !! @param[in] batch_idx 書き出し対象のバッチ番号。
  !! @param[in] processed_particles これまでの累積処理粒子数。
  !! @param[in] rel_change この時点の相対変化量。
  !! @param[in] q_elem 要素電荷配列 [C]。
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

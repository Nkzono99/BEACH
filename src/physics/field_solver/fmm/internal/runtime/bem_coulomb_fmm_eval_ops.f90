!> Coulomb FMM 電場評価。
module bem_coulomb_fmm_eval_ops
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_panel_kernel, only: panel_potential_field_translated, panel_field_sum_translated, &
                              panel_potential_sum_translated, panel_field_sum_translated_sp, &
                              panel_potential_sum_translated_sp
  use bem_coulomb_fmm_types, only: fmm_plan_type, fmm_state_type, fmm_tree_type
  use bem_coulomb_fmm_basis, only: build_axis_powers
  use bem_coulomb_fmm_periodic, only: wrap_periodic2_point, use_periodic2_cached_kneq0
  use bem_coulomb_fmm_tree_utils, only: octant_index
  use bem_periodic_zero_mode_eval, only: eval_periodic_zero_mode, zero_mode_trace_plus
  use bem_periodic_shell_mean, only: eval_periodic_shell_mean
  implicit none
  private

  public :: core_eval_points_impl
  public :: core_eval_direct_points_impl
  public :: core_eval_point_impl
  public :: core_eval_potential_points_impl
  public :: core_eval_direct_potential_points_impl
  public :: core_eval_potential_point_impl

contains

  !> 複数の評価点で電場を計算する。
  !! @param[in] plan 構築済みの FMM 計画。
  !! @param[inout] state 評価に使う FMM state。
  !! @param[in] target_pos 評価点位置 `(3,m)` [m]。
  !! @param[out] e 電場ベクトル `(3,m)` [V/m]。
  subroutine core_eval_points_impl(plan, state, target_pos, e)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    real(dp), intent(in) :: target_pos(:, :)
    real(dp), intent(out) :: e(:, :)
    integer(i32) :: i, ntarget

    if (size(target_pos, 1) /= 3) error stop 'FMM core expects target_pos(3,m).'
    if (size(e, 1) /= 3 .or. size(e, 2) /= size(target_pos, 2)) then
      error stop 'FMM eval_points expects e(3,m).'
    end if
    ntarget = int(size(target_pos, 2), i32)

    !$omp parallel do default(none) schedule(static) &
    !$omp   shared(plan, state, target_pos, e, ntarget) private(i)
    do i = 1_i32, ntarget
      call core_eval_point_xyz_impl( &
        plan, state, target_pos(1, i), target_pos(2, i), target_pos(3, i), e(1, i), e(2, i), e(3, i) &
        )
    end do
    !$omp end parallel do
  end subroutine core_eval_points_impl

  !> 非周期計画の全ソース直接和で複数点の電場を計算する。ソース幾何と電荷 `src_q` だけを使う。
  subroutine core_eval_direct_points_impl(plan, src_q, target_pos, e)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: src_q(:)
    real(dp), intent(in) :: target_pos(:, :)
    real(dp), intent(out) :: e(:, :)
    integer(i32) :: i, ntarget
    real(dp) :: soft2

    if (size(target_pos, 1) /= 3) error stop 'FMM core expects target_pos(3,m).'
    if (size(e, 1) /= 3 .or. size(e, 2) /= size(target_pos, 2)) then
      error stop 'FMM eval_direct_points expects e(3,m).'
    end if
    if (plan%options%use_periodic2) then
      error stop 'FMM exact-direct evaluation requires a non-periodic plan.'
    end if

    e = 0.0_dp
    if (.not. plan%sources_ready) return
    if (size(src_q) /= plan%nsrc) error stop 'FMM eval_direct_points expects src_q(nsrc).'
    ntarget = int(size(target_pos, 2), i32)
    soft2 = plan%options%softening*plan%options%softening

    !$omp parallel do default(none) schedule(static) &
    !$omp   shared(plan, src_q, target_pos, e, ntarget, soft2) private(i)
    do i = 1_i32, ntarget
      call eval_direct_all_sources_scalar( &
        plan, src_q, target_pos(1, i), target_pos(2, i), target_pos(3, i), soft2, e(1, i), e(2, i), e(3, i) &
        )
    end do
    !$omp end parallel do
  end subroutine core_eval_direct_points_impl

  !> 1 点で電場を計算する。
  !! @param[in] plan 構築済みの FMM 計画。
  !! @param[inout] state 評価に使う FMM state。
  !! @param[in] r 評価点位置 `(x,y,z)` [m]。
  !! @param[out] e 電場ベクトル `(x,y,z)` [V/m]。
  subroutine core_eval_point_impl(plan, state, r, e)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    real(dp), intent(in) :: r(3)
    real(dp), intent(out) :: e(3)

    call core_eval_point_xyz_impl(plan, state, r(1), r(2), r(3), e(1), e(2), e(3))
  end subroutine core_eval_point_impl

  !> 複数の評価点で電位を計算する。
  !! @param[in] plan 構築済みの FMM 計画。
  !! @param[inout] state 評価に使う FMM state。
  !! @param[in] target_pos 評価点位置 `(3,m)` [m]。
  !! @param[out] phi 電位 `(m)` [V]（`k_coulomb` は含まない）。
  subroutine core_eval_potential_points_impl(plan, state, target_pos, phi)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    real(dp), intent(in) :: target_pos(:, :)
    real(dp), intent(out) :: phi(:)
    integer(i32) :: i, ntarget

    if (size(target_pos, 1) /= 3) error stop 'FMM core expects target_pos(3,m).'
    if (size(phi) /= size(target_pos, 2)) then
      error stop 'FMM eval_potential_points expects phi(m).'
    end if
    ntarget = int(size(target_pos, 2), i32)

    !$omp parallel do default(none) schedule(static) &
    !$omp   shared(plan, state, target_pos, phi, ntarget) private(i)
    do i = 1_i32, ntarget
      call core_eval_potential_point_xyz_impl( &
        plan, state, target_pos(1, i), target_pos(2, i), target_pos(3, i), phi(i) &
        )
    end do
    !$omp end parallel do
  end subroutine core_eval_potential_points_impl

  !> 非周期計画の全ソース直接和で複数点の電位を計算する。ソース幾何と電荷 `src_q` だけを使う。
  subroutine core_eval_direct_potential_points_impl(plan, src_q, target_pos, phi)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: src_q(:)
    real(dp), intent(in) :: target_pos(:, :)
    real(dp), intent(out) :: phi(:)
    integer(i32) :: i, ntarget
    real(dp) :: soft2

    if (size(target_pos, 1) /= 3) error stop 'FMM core expects target_pos(3,m).'
    if (size(phi) /= size(target_pos, 2)) then
      error stop 'FMM eval_direct_potential_points expects phi(m).'
    end if
    if (plan%options%use_periodic2) then
      error stop 'FMM exact-direct evaluation requires a non-periodic plan.'
    end if

    phi = 0.0_dp
    if (.not. plan%sources_ready) return
    if (size(src_q) /= plan%nsrc) error stop 'FMM eval_direct_potential_points expects src_q(nsrc).'
    ntarget = int(size(target_pos, 2), i32)
    soft2 = plan%options%softening*plan%options%softening

    !$omp parallel do default(none) schedule(static) &
    !$omp   shared(plan, src_q, target_pos, phi, ntarget, soft2) private(i)
    do i = 1_i32, ntarget
      call eval_direct_all_sources_potential_scalar( &
        plan, src_q, target_pos(1, i), target_pos(2, i), target_pos(3, i), soft2, phi(i) &
        )
    end do
    !$omp end parallel do
  end subroutine core_eval_direct_potential_points_impl

  !> 1 点で電位を計算する。
  !! @param[in] plan 構築済みの FMM 計画。
  !! @param[inout] state 評価に使う FMM state。
  !! @param[in] r 評価点位置 `(x,y,z)` [m]。
  !! @param[out] phi 電位 [V]（`k_coulomb` は含まない）。
  subroutine core_eval_potential_point_impl(plan, state, r, phi)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    real(dp), intent(in) :: r(3)
    real(dp), intent(out) :: phi

    call core_eval_potential_point_xyz_impl(plan, state, r(1), r(2), r(3), phi)
  end subroutine core_eval_potential_point_impl

  !> 1 点評価の本体処理を行う。
  !! @param[in] plan 構築済みの FMM 計画。
  !! @param[inout] state 評価に使う FMM state。
  !! @param[in] rx 評価点 x 座標 [m]。
  !! @param[in] ry 評価点 y 座標 [m]。
  !! @param[in] rz 評価点 z 座標 [m]。
  !! @param[out] ex 電場 x 成分 [V/m]。
  !! @param[out] ey 電場 y 成分 [V/m]。
  !! @param[out] ez 電場 z 成分 [V/m]。
  subroutine core_eval_point_xyz_impl(plan, state, rx, ry, rz, ex, ey, ez)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    real(dp), intent(in) :: rx, ry, rz
    real(dp), intent(out) :: ex, ey, ez
    integer(i32) :: leaf_node, leaf_slot
    integer(i32) :: axes(2)
    real(dp) :: rt(3), soft2, leaf_center(3)

    ex = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    if (.not. plan%built .or. .not. state%ready) return

    rt = [rx, ry, rz]
    if (plan%options%use_periodic2) call wrap_periodic2_point(plan, rt)
    soft2 = plan%options%softening*plan%options%softening
    call locate_target_leaf(plan, rt, leaf_node, leaf_center)
    leaf_slot = 0_i32
    if (leaf_node > 0_i32) leaf_slot = plan%leaf_slot_of_node(leaf_node)
    if (leaf_slot <= 0_i32) then
      call apply_direct_fallback_field(plan, state, rt, soft2, ex, ey, ez)
      return
    end if

    call accumulate_leaf_local_expansion(plan, state, leaf_node, leaf_center, rt, ex, ey, ez)
    call accumulate_anchor_local_expansion(plan, state, rt, ex, ey, ez)
    axes = active_periodic_axes(plan)
    call accumulate_near_direct_field(plan, state, leaf_slot, rt, soft2, axes, ex, ey, ez)
    call subtract_cached_periodic_k0_field(plan, state, rt, ex, ey, ez)
  end subroutine core_eval_point_xyz_impl

  !> 1 点評価の電位本体処理を行う。
  subroutine core_eval_potential_point_xyz_impl(plan, state, rx, ry, rz, phi)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    real(dp), intent(in) :: rx, ry, rz
    real(dp), intent(out) :: phi
    integer(i32) :: leaf_node, leaf_slot
    integer(i32) :: axes(2)
    real(dp) :: rt(3), soft2, leaf_center(3)

    phi = 0.0d0
    if (.not. plan%built .or. .not. state%ready) return

    rt = [rx, ry, rz]
    if (plan%options%use_periodic2) call wrap_periodic2_point(plan, rt)
    soft2 = plan%options%softening*plan%options%softening
    call locate_target_leaf(plan, rt, leaf_node, leaf_center)
    leaf_slot = 0_i32
    if (leaf_node > 0_i32) leaf_slot = plan%leaf_slot_of_node(leaf_node)
    if (leaf_slot <= 0_i32) then
      call apply_direct_fallback_potential(plan, state, rt, soft2, phi)
      return
    end if

    call accumulate_leaf_local_potential_expansion(plan, state, leaf_node, leaf_center, rt, phi)
    call accumulate_anchor_local_potential_expansion(plan, state, rt, phi)
    axes = active_periodic_axes(plan)
    call accumulate_near_direct_potential(plan, state, leaf_slot, rt, soft2, axes, phi)
    call subtract_cached_periodic_k0_potential(plan, state, rt, phi)
  end subroutine core_eval_potential_point_xyz_impl

  subroutine subtract_cached_periodic_k0_field(plan, state, rt, ex, ey, ez)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: rt(3)
    real(dp), intent(inout) :: ex, ey, ez
    real(dp) :: potential_dummy, field_si, field_raw(3), shell_value, shell_slope

    if (.not. plan%periodic_k0_ready) return
    call eval_periodic_zero_mode( &
      plan%periodic_k0_plan, state%periodic_k0_state, rt(plan%periodic_k0_free_axis), &
      zero_mode_trace_plus, potential_dummy, field_si &
      )
    ! 有限画像殻の面平均は、対称な平均場と滑らかな残り S の和。S の電場 -dS/dz も引く。
    call eval_periodic_shell_mean( &
      plan%periodic_shell_mean, state%periodic_shell_mean_state, rt(plan%periodic_k0_free_axis), &
      shell_value, shell_slope &
      )
    field_raw = [ex, ey, ez]
    field_raw(plan%periodic_k0_free_axis) = field_raw(plan%periodic_k0_free_axis) - field_si/k_coulomb + shell_slope
    ex = field_raw(1)
    ey = field_raw(2)
    ez = field_raw(3)
  end subroutine subtract_cached_periodic_k0_field

  subroutine subtract_cached_periodic_k0_potential(plan, state, rt, phi)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: rt(3)
    real(dp), intent(inout) :: phi
    real(dp) :: potential_si, field_dummy, shell_value, shell_slope

    if (.not. plan%periodic_k0_ready) return
    call eval_periodic_zero_mode( &
      plan%periodic_k0_plan, state%periodic_k0_state, rt(plan%periodic_k0_free_axis), &
      zero_mode_trace_plus, potential_si, field_dummy &
      )
    call eval_periodic_shell_mean( &
      plan%periodic_shell_mean, state%periodic_shell_mean_state, rt(plan%periodic_k0_free_axis), &
      shell_value, shell_slope &
      )
    phi = phi - potential_si/k_coulomb - shell_value
  end subroutine subtract_cached_periodic_k0_potential

  subroutine apply_direct_fallback_field(plan, state, rt, soft2, ex, ey, ez)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: rt(3)
    real(dp), intent(in) :: soft2
    real(dp), intent(inout) :: ex, ey, ez

    if (use_periodic2_cached_kneq0(plan)) then
      error stop 'cached_kneq0 evaluation requires targets inside the configured target box.'
    end if
    call eval_direct_all_sources_scalar(plan, state%src_q, rt(1), rt(2), rt(3), soft2, ex, ey, ez)
  end subroutine apply_direct_fallback_field

  subroutine apply_direct_fallback_potential(plan, state, rt, soft2, phi)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: rt(3)
    real(dp), intent(in) :: soft2
    real(dp), intent(inout) :: phi

    if (use_periodic2_cached_kneq0(plan)) then
      error stop 'cached_kneq0 potential evaluation requires targets inside the configured target box.'
    end if
    call eval_direct_all_sources_potential_scalar(plan, state%src_q, rt(1), rt(2), rt(3), soft2, phi)
  end subroutine apply_direct_fallback_potential

  subroutine accumulate_leaf_local_expansion(plan, state, leaf_node, leaf_center, rt, ex, ey, ez)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    integer(i32), intent(in) :: leaf_node
    real(dp), intent(in) :: leaf_center(3), rt(3)
    real(dp), intent(inout) :: ex, ey, ez
    integer(i32) :: order, term_idx
    real(dp) :: dr(3), monomial
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    order = active_local_expansion_order(plan)
    if (order <= 0_i32) return
    if (state%local_active(leaf_node) == 0_i32) return
    if (plan%eval_term_count <= 0_i32) return

    dr = rt - leaf_center
    call build_axis_powers(dr, order, xpow, ypow, zpow)
    do term_idx = 1_i32, plan%eval_term_count
      monomial = xpow(plan%eval_exp(1, term_idx))*ypow(plan%eval_exp(2, term_idx)) &
                 *zpow(plan%eval_exp(3, term_idx))*plan%eval_inv_factorial(term_idx)
      ex = ex - state%local(plan%eval_deriv_idx(1, term_idx), leaf_node)*monomial
      ey = ey - state%local(plan%eval_deriv_idx(2, term_idx), leaf_node)*monomial
      ez = ez - state%local(plan%eval_deriv_idx(3, term_idx), leaf_node)*monomial
    end do
  end subroutine accumulate_leaf_local_expansion

  subroutine accumulate_leaf_local_potential_expansion(plan, state, leaf_node, leaf_center, rt, phi)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    integer(i32), intent(in) :: leaf_node
    real(dp), intent(in) :: leaf_center(3), rt(3)
    real(dp), intent(inout) :: phi
    integer(i32) :: order, alpha_idx
    real(dp) :: dr(3), monomial
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    order = active_local_expansion_order(plan)
    if (state%local_active(leaf_node) == 0_i32) return

    dr = rt - leaf_center
    call build_axis_powers(dr, order, xpow, ypow, zpow)
    do alpha_idx = 1_i32, plan%ncoef
      monomial = xpow(plan%alpha(1, alpha_idx))*ypow(plan%alpha(2, alpha_idx)) &
                 *zpow(plan%alpha(3, alpha_idx))/plan%alpha_factorial(alpha_idx)
      phi = phi + state%local(alpha_idx, leaf_node)*monomial
    end do
  end subroutine accumulate_leaf_local_potential_expansion

  !> cached_kneq0 の遠方補正: 点を含む anchor 格子セルの番号と中心。box の外の点は端の anchor に寄せる。
  pure subroutine locate_periodic_anchor(plan, rt, anchor, center)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: rt(3)
    integer(i32), intent(out) :: anchor
    real(dp), intent(out) :: center(3)
    integer(i32) :: cells(3), index3(3)
    real(dp) :: span(3)

    cells = plan%options%periodic_operator_anchor_cells
    span = plan%options%target_box_max - plan%options%target_box_min
    index3 = min(cells - 1_i32, max(0_i32, floor((rt - plan%options%target_box_min)/span*real(cells, dp), kind=i32)))
    anchor = 1_i32 + index3(1) + cells(1)*(index3(2) + cells(2)*index3(3))
    center = plan%options%target_box_min + (real(index3, dp) + 0.5_dp)*2.0_dp*plan%periodic_operator%anchor_half
  end subroutine locate_periodic_anchor

  !> 遠方補正の anchor 局所展開の電場を加える。
  subroutine accumulate_anchor_local_expansion(plan, state, rt, ex, ey, ez)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: rt(3)
    real(dp), intent(inout) :: ex, ey, ez
    integer(i32) :: order, term_idx, anchor
    real(dp) :: dr(3), monomial, center(3)
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    if (.not. plan%periodic_operator%ready) return
    if (.not. allocated(state%anchor_local)) return
    order = active_local_expansion_order(plan)
    if (order <= 0_i32 .or. plan%eval_term_count <= 0_i32) return
    call locate_periodic_anchor(plan, rt, anchor, center)
    if (state%anchor_local_active(anchor) == 0_i32) return

    dr = rt - center
    call build_axis_powers(dr, order, xpow, ypow, zpow)
    do term_idx = 1_i32, plan%eval_term_count
      monomial = xpow(plan%eval_exp(1, term_idx))*ypow(plan%eval_exp(2, term_idx)) &
                 *zpow(plan%eval_exp(3, term_idx))*plan%eval_inv_factorial(term_idx)
      ex = ex - state%anchor_local(plan%eval_deriv_idx(1, term_idx), anchor)*monomial
      ey = ey - state%anchor_local(plan%eval_deriv_idx(2, term_idx), anchor)*monomial
      ez = ez - state%anchor_local(plan%eval_deriv_idx(3, term_idx), anchor)*monomial
    end do
  end subroutine accumulate_anchor_local_expansion

  !> 遠方補正の anchor 局所展開の電位を加える。
  subroutine accumulate_anchor_local_potential_expansion(plan, state, rt, phi)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: rt(3)
    real(dp), intent(inout) :: phi
    integer(i32) :: order, alpha_idx, anchor
    real(dp) :: dr(3), monomial, center(3)
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    if (.not. plan%periodic_operator%ready) return
    if (.not. allocated(state%anchor_local)) return
    order = active_local_expansion_order(plan)
    call locate_periodic_anchor(plan, rt, anchor, center)
    if (state%anchor_local_active(anchor) == 0_i32) return

    dr = rt - center
    call build_axis_powers(dr, order, xpow, ypow, zpow)
    do alpha_idx = 1_i32, plan%ncoef
      monomial = xpow(plan%alpha(1, alpha_idx))*ypow(plan%alpha(2, alpha_idx)) &
                 *zpow(plan%alpha(3, alpha_idx))/plan%alpha_factorial(alpha_idx)
      phi = phi + state%anchor_local(alpha_idx, anchor)*monomial
    end do
  end subroutine accumulate_anchor_local_potential_expansion

  pure function active_local_expansion_order(plan) result(order)
    type(fmm_plan_type), intent(in) :: plan
    integer(i32) :: order

    order = plan%options%order
  end function active_local_expansion_order

  pure function active_periodic_axes(plan) result(axes)
    type(fmm_plan_type), intent(in) :: plan
    integer(i32) :: axes(2)

    axes = 0_i32
    if (plan%options%use_periodic2) axes = plan%options%periodic_axes
  end function active_periodic_axes

  subroutine accumulate_near_direct_field(plan, state, leaf_slot, rt, soft2, axes, ex, ey, ez)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    integer(i32), intent(in) :: leaf_slot
    real(dp), intent(in) :: rt(3)
    real(dp), intent(in) :: soft2
    integer(i32), intent(in) :: axes(2)
    real(dp), intent(inout) :: ex, ey, ez
    integer(i32) :: near_pos, idx, axis1, axis2
    integer(i32) :: near_source_begin, near_source_end
    real(dp) :: shift1, shift2

    if (leaf_slot <= 0_i32) return
    near_source_begin = plan%near_source_start(leaf_slot)
    near_source_end = plan%near_source_start(leaf_slot + 1_i32) - 1_i32
    if (near_source_end < near_source_begin) return

    if (plan%panel_source .and. plan%options%near_field_single) then
      call panel_field_sum_translated_sp( &
        near_source_end - near_source_begin + 1_i32, &
        plan%near_panels_sp%vx(1, near_source_begin), plan%near_panels_sp%vy(1, near_source_begin), &
        plan%near_panels_sp%vz(1, near_source_begin), plan%near_panels_sp%nx(near_source_begin), &
        plan%near_panels_sp%ny(near_source_begin), plan%near_panels_sp%nz(near_source_begin), &
        plan%near_panels_sp%edge(1, near_source_begin), plan%near_panels_sp%ox(1, near_source_begin), &
        plan%near_panels_sp%oy(1, near_source_begin), plan%near_panels_sp%oz(1, near_source_begin), &
        plan%near_panels_sp%area(near_source_begin), &
        plan%near_source_idx(near_source_begin), state%src_q, plan%panel_geometry, axes, &
        plan%near_source_shift1(near_source_begin), plan%near_source_shift2(near_source_begin), &
        rt - plan%near_panels_sp%leaf_center(:, leaf_slot), rt, ex, ey, ez &
        )
      return
    end if
    if (plan%panel_source) then
      call panel_field_sum_translated( &
        near_source_end - near_source_begin + 1_i32, &
        plan%near_panels%vx(1, near_source_begin), plan%near_panels%vy(1, near_source_begin), &
        plan%near_panels%vz(1, near_source_begin), plan%near_panels%nx(near_source_begin), &
        plan%near_panels%ny(near_source_begin), plan%near_panels%nz(near_source_begin), &
        plan%near_panels%edge(1, near_source_begin), plan%near_panels%ox(1, near_source_begin), &
        plan%near_panels%oy(1, near_source_begin), plan%near_panels%oz(1, near_source_begin), &
        plan%near_panels%area(near_source_begin), &
        plan%near_source_idx(near_source_begin), state%src_q, plan%panel_geometry, rt, ex, ey, ez &
        )
      return
    end if

    if (plan%options%use_periodic2) then
      axis1 = axes(1)
      axis2 = axes(2)
      do near_pos = near_source_begin, near_source_end
        idx = plan%near_source_idx(near_pos)
        shift1 = plan%near_source_shift1(near_pos)
        shift2 = plan%near_source_shift2(near_pos)
        call accumulate_point_charge_shifted_field( &
          state%src_q(idx), plan%src_pos(1, idx), plan%src_pos(2, idx), plan%src_pos(3, idx), shift1, shift2, &
          axis1, axis2, rt(1), rt(2), rt(3), soft2, ex, ey, ez &
          )
      end do
    else
      do near_pos = near_source_begin, near_source_end
        idx = plan%near_source_idx(near_pos)
        call accumulate_point_charge_field( &
          state%src_q(idx), plan%src_pos(1, idx), plan%src_pos(2, idx), plan%src_pos(3, idx), &
          rt(1), rt(2), rt(3), soft2, ex, ey, ez &
          )
      end do
    end if
  end subroutine accumulate_near_direct_field

  subroutine accumulate_near_direct_potential(plan, state, leaf_slot, rt, soft2, axes, phi)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    integer(i32), intent(in) :: leaf_slot
    real(dp), intent(in) :: rt(3)
    real(dp), intent(in) :: soft2
    integer(i32), intent(in) :: axes(2)
    real(dp), intent(inout) :: phi
    integer(i32) :: near_pos, idx, axis1, axis2
    integer(i32) :: near_source_begin, near_source_end
    real(dp) :: shift1, shift2

    if (leaf_slot <= 0_i32) return
    near_source_begin = plan%near_source_start(leaf_slot)
    near_source_end = plan%near_source_start(leaf_slot + 1_i32) - 1_i32
    if (near_source_end < near_source_begin) return

    if (plan%panel_source .and. plan%options%near_field_single) then
      call panel_potential_sum_translated_sp( &
        near_source_end - near_source_begin + 1_i32, &
        plan%near_panels_sp%vx(1, near_source_begin), plan%near_panels_sp%vy(1, near_source_begin), &
        plan%near_panels_sp%vz(1, near_source_begin), plan%near_panels_sp%nx(near_source_begin), &
        plan%near_panels_sp%ny(near_source_begin), plan%near_panels_sp%nz(near_source_begin), &
        plan%near_panels_sp%edge(1, near_source_begin), plan%near_panels_sp%ox(1, near_source_begin), &
        plan%near_panels_sp%oy(1, near_source_begin), plan%near_panels_sp%oz(1, near_source_begin), &
        plan%near_panels_sp%area(near_source_begin), &
        plan%near_source_idx(near_source_begin), state%src_q, plan%panel_geometry, axes, &
        plan%near_source_shift1(near_source_begin), plan%near_source_shift2(near_source_begin), &
        rt - plan%near_panels_sp%leaf_center(:, leaf_slot), rt, phi &
        )
      return
    end if
    if (plan%panel_source) then
      call panel_potential_sum_translated( &
        near_source_end - near_source_begin + 1_i32, &
        plan%near_panels%vx(1, near_source_begin), plan%near_panels%vy(1, near_source_begin), &
        plan%near_panels%vz(1, near_source_begin), plan%near_panels%nx(near_source_begin), &
        plan%near_panels%ny(near_source_begin), plan%near_panels%nz(near_source_begin), &
        plan%near_panels%edge(1, near_source_begin), plan%near_panels%ox(1, near_source_begin), &
        plan%near_panels%oy(1, near_source_begin), plan%near_panels%oz(1, near_source_begin), &
        plan%near_panels%area(near_source_begin), &
        plan%near_source_idx(near_source_begin), state%src_q, plan%panel_geometry, rt, phi &
        )
      return
    end if

    if (plan%options%use_periodic2) then
      axis1 = axes(1)
      axis2 = axes(2)
      do near_pos = near_source_begin, near_source_end
        idx = plan%near_source_idx(near_pos)
        shift1 = plan%near_source_shift1(near_pos)
        shift2 = plan%near_source_shift2(near_pos)
        call accumulate_point_charge_shifted_potential( &
          state%src_q(idx), plan%src_pos(1, idx), plan%src_pos(2, idx), plan%src_pos(3, idx), shift1, shift2, &
          axis1, axis2, rt(1), rt(2), rt(3), soft2, phi &
          )
      end do
    else
      do near_pos = near_source_begin, near_source_end
        idx = plan%near_source_idx(near_pos)
        call accumulate_point_charge_potential( &
          state%src_q(idx), plan%src_pos(1, idx), plan%src_pos(2, idx), plan%src_pos(3, idx), &
          rt(1), rt(2), rt(3), soft2, phi &
          )
      end do
    end if
  end subroutine accumulate_near_direct_potential

  subroutine eval_direct_all_sources_scalar(plan, src_q, tx, ty, tz, soft2, ex, ey, ez)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: src_q(:)
    real(dp), intent(in) :: tx, ty, tz
    real(dp), intent(in) :: soft2
    real(dp), intent(out) :: ex, ey, ez
    integer(i32) :: idx, axis1, axis2, nshift
    integer(i32) :: axes(2)

    ex = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    axes = active_periodic_axes(plan)
    axis1 = axes(1)
    axis2 = axes(2)
    if (plan%panel_source) then
      if (plan%options%use_periodic2) then
        nshift = size(plan%shift_axis1)
        do idx = 1_i32, plan%nsrc
          call accumulate_panel_images_field(plan, src_q, idx, [tx, ty, tz], axes, nshift, ex, ey, ez)
        end do
      else
        do idx = 1_i32, plan%nsrc
          call accumulate_panel_shifted_field(plan, src_q, idx, [tx, ty, tz], axes, 0.0_dp, 0.0_dp, ex, ey, ez)
        end do
      end if
      return
    end if
    if (plan%options%use_periodic2) then
      nshift = size(plan%shift_axis1)
      do idx = 1_i32, plan%nsrc
        call accumulate_point_charge_images_field( &
          src_q(idx), plan%src_pos(1, idx), plan%src_pos(2, idx), plan%src_pos(3, idx), &
          tx, ty, tz, soft2, axis1, axis2, plan%shift_axis1, plan%shift_axis2, nshift, ex, ey, ez &
          )
      end do
    else
      do idx = 1_i32, plan%nsrc
        call accumulate_point_charge_field( &
          src_q(idx), plan%src_pos(1, idx), plan%src_pos(2, idx), plan%src_pos(3, idx), &
          tx, ty, tz, soft2, ex, ey, ez &
          )
      end do
    end if
  end subroutine eval_direct_all_sources_scalar

  subroutine eval_direct_all_sources_potential_scalar(plan, src_q, tx, ty, tz, soft2, phi)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: src_q(:)
    real(dp), intent(in) :: tx, ty, tz
    real(dp), intent(in) :: soft2
    real(dp), intent(out) :: phi
    integer(i32) :: idx, axis1, axis2, nshift
    integer(i32) :: axes(2)

    phi = 0.0d0
    axes = active_periodic_axes(plan)
    axis1 = axes(1)
    axis2 = axes(2)
    if (plan%panel_source) then
      if (plan%options%use_periodic2) then
        nshift = size(plan%shift_axis1)
        do idx = 1_i32, plan%nsrc
          call accumulate_panel_images_potential(plan, src_q, idx, [tx, ty, tz], axes, nshift, phi)
        end do
      else
        do idx = 1_i32, plan%nsrc
          call accumulate_panel_shifted_potential(plan, src_q, idx, [tx, ty, tz], axes, 0.0_dp, 0.0_dp, phi)
        end do
      end if
      return
    end if
    if (plan%options%use_periodic2) then
      nshift = size(plan%shift_axis1)
      do idx = 1_i32, plan%nsrc
        call accumulate_point_charge_images_potential( &
          src_q(idx), plan%src_pos(1, idx), plan%src_pos(2, idx), plan%src_pos(3, idx), &
          tx, ty, tz, soft2, axis1, axis2, plan%shift_axis1, plan%shift_axis2, nshift, phi &
          )
      end do
    else
      do idx = 1_i32, plan%nsrc
        call accumulate_point_charge_potential( &
          src_q(idx), plan%src_pos(1, idx), plan%src_pos(2, idx), plan%src_pos(3, idx), &
          tx, ty, tz, soft2, phi &
          )
      end do
    end if
  end subroutine eval_direct_all_sources_potential_scalar

  subroutine accumulate_panel_images_field(plan, src_q, idx, target, axes, nshift, ex, ey, ez)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: src_q(:)
    integer(i32), intent(in) :: idx, axes(2), nshift
    real(dp), intent(in) :: target(3)
    real(dp), intent(inout) :: ex, ey, ez
    integer(i32) :: image1, image2

    do image1 = 1_i32, nshift
      do image2 = 1_i32, nshift
        call accumulate_panel_shifted_field( &
          plan, src_q, idx, target, axes, plan%shift_axis1(image1), plan%shift_axis2(image2), ex, ey, ez &
          )
      end do
    end do
  end subroutine accumulate_panel_images_field

  subroutine accumulate_panel_images_potential(plan, src_q, idx, target, axes, nshift, phi)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: src_q(:)
    integer(i32), intent(in) :: idx, axes(2), nshift
    real(dp), intent(in) :: target(3)
    real(dp), intent(inout) :: phi
    integer(i32) :: image1, image2

    do image1 = 1_i32, nshift
      do image2 = 1_i32, nshift
        call accumulate_panel_shifted_potential( &
          plan, src_q, idx, target, axes, plan%shift_axis1(image1), plan%shift_axis2(image2), phi &
          )
      end do
    end do
  end subroutine accumulate_panel_images_potential

  subroutine accumulate_panel_shifted_field(plan, src_q, idx, target, axes, shift1, shift2, ex, ey, ez)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: src_q(:)
    integer(i32), intent(in) :: idx, axes(2)
    real(dp), intent(in) :: target(3), shift1, shift2
    real(dp), intent(inout) :: ex, ey, ez
    real(dp) :: potential, field(3)

    if (abs(src_q(idx)) <= tiny(1.0_dp)) return
    call panel_potential_field_translated( &
      plan%panel_geometry(idx), periodic_shift_vector(axes, shift1, shift2), src_q(idx), target, potential, field &
      )
    ex = ex + field(1)/k_coulomb
    ey = ey + field(2)/k_coulomb
    ez = ez + field(3)/k_coulomb
  end subroutine accumulate_panel_shifted_field

  subroutine accumulate_panel_shifted_potential(plan, src_q, idx, target, axes, shift1, shift2, phi)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: src_q(:)
    integer(i32), intent(in) :: idx, axes(2)
    real(dp), intent(in) :: target(3), shift1, shift2
    real(dp), intent(inout) :: phi
    real(dp) :: potential, field(3)

    if (abs(src_q(idx)) <= tiny(1.0_dp)) return
    call panel_potential_field_translated( &
      plan%panel_geometry(idx), periodic_shift_vector(axes, shift1, shift2), src_q(idx), target, potential, field &
      )
    phi = phi + potential/k_coulomb
  end subroutine accumulate_panel_shifted_potential

  !> Image translation of a source along the periodic axes; zero on the free axis.
  pure function periodic_shift_vector(axes, shift1, shift2) result(shift)
    integer(i32), intent(in) :: axes(2)
    real(dp), intent(in) :: shift1, shift2
    real(dp) :: shift(3)

    shift = 0.0_dp
    if (axes(1) > 0_i32) shift(axes(1)) = shift(axes(1)) + shift1
    if (axes(2) > 0_i32) shift(axes(2)) = shift(axes(2)) + shift2
  end function periodic_shift_vector

  pure subroutine accumulate_point_charge_field(q, sx, sy, sz, tx, ty, tz, soft2, ex, ey, ez)
    real(dp), intent(in) :: q, sx, sy, sz, tx, ty, tz, soft2
    real(dp), intent(inout) :: ex, ey, ez
    real(dp) :: dx, dy, dz, r2, inv_r3

    if (abs(q) <= tiny(1.0d0)) return
    dx = tx - sx
    dy = ty - sy
    dz = tz - sz
    r2 = dx*dx + dy*dy + dz*dz + soft2
    if (r2 <= tiny(1.0d0)) return
    inv_r3 = 1.0d0/(sqrt(r2)*r2)
    ex = ex + q*inv_r3*dx
    ey = ey + q*inv_r3*dy
    ez = ez + q*inv_r3*dz
  end subroutine accumulate_point_charge_field

  pure subroutine accumulate_point_charge_potential(q, sx, sy, sz, tx, ty, tz, soft2, phi)
    real(dp), intent(in) :: q, sx, sy, sz, tx, ty, tz, soft2
    real(dp), intent(inout) :: phi
    real(dp) :: dx, dy, dz, r2

    if (abs(q) <= tiny(1.0d0)) return
    dx = tx - sx
    dy = ty - sy
    dz = tz - sz
    r2 = dx*dx + dy*dy + dz*dz + soft2
    if (r2 <= tiny(1.0d0)) return
    phi = phi + q/sqrt(r2)
  end subroutine accumulate_point_charge_potential

  pure subroutine accumulate_point_charge_images_field( &
    q, sx, sy, sz, tx, ty, tz, soft2, axis1, axis2, shift_axis1, shift_axis2, nshift, ex, ey, ez &
    )
    real(dp), intent(in) :: q, sx, sy, sz, tx, ty, tz, soft2
    integer(i32), intent(in) :: axis1, axis2, nshift
    real(dp), intent(in) :: shift_axis1(:), shift_axis2(:)
    real(dp), intent(inout) :: ex, ey, ez
    integer(i32) :: img_i, img_j
    real(dp) :: sx_img, sy_img, sz_img

    if (abs(q) <= tiny(1.0d0)) return
    do img_i = 1_i32, nshift
      do img_j = 1_i32, nshift
        sx_img = sx
        sy_img = sy
        sz_img = sz
        if (axis1 == 1_i32) sx_img = sx_img + shift_axis1(img_i)
        if (axis1 == 2_i32) sy_img = sy_img + shift_axis1(img_i)
        if (axis1 == 3_i32) sz_img = sz_img + shift_axis1(img_i)
        if (axis2 == 1_i32) sx_img = sx_img + shift_axis2(img_j)
        if (axis2 == 2_i32) sy_img = sy_img + shift_axis2(img_j)
        if (axis2 == 3_i32) sz_img = sz_img + shift_axis2(img_j)
        call accumulate_point_charge_field(q, sx_img, sy_img, sz_img, tx, ty, tz, soft2, ex, ey, ez)
      end do
    end do
  end subroutine accumulate_point_charge_images_field

  pure subroutine accumulate_point_charge_images_potential( &
    q, sx, sy, sz, tx, ty, tz, soft2, axis1, axis2, shift_axis1, shift_axis2, nshift, phi &
    )
    real(dp), intent(in) :: q, sx, sy, sz, tx, ty, tz, soft2
    integer(i32), intent(in) :: axis1, axis2, nshift
    real(dp), intent(in) :: shift_axis1(:), shift_axis2(:)
    real(dp), intent(inout) :: phi
    integer(i32) :: img_i, img_j
    real(dp) :: sx_img, sy_img, sz_img

    if (abs(q) <= tiny(1.0d0)) return
    do img_i = 1_i32, nshift
      do img_j = 1_i32, nshift
        sx_img = sx
        sy_img = sy
        sz_img = sz
        if (axis1 == 1_i32) sx_img = sx_img + shift_axis1(img_i)
        if (axis1 == 2_i32) sy_img = sy_img + shift_axis1(img_i)
        if (axis1 == 3_i32) sz_img = sz_img + shift_axis1(img_i)
        if (axis2 == 1_i32) sx_img = sx_img + shift_axis2(img_j)
        if (axis2 == 2_i32) sy_img = sy_img + shift_axis2(img_j)
        if (axis2 == 3_i32) sz_img = sz_img + shift_axis2(img_j)
        call accumulate_point_charge_potential(q, sx_img, sy_img, sz_img, tx, ty, tz, soft2, phi)
      end do
    end do
  end subroutine accumulate_point_charge_images_potential

  pure subroutine accumulate_point_charge_shifted_field( &
    q, sx, sy, sz, shift1, shift2, axis1, axis2, tx, ty, tz, soft2, ex, ey, ez &
    )
    real(dp), intent(in) :: q, sx, sy, sz, shift1, shift2, tx, ty, tz, soft2
    integer(i32), intent(in) :: axis1, axis2
    real(dp), intent(inout) :: ex, ey, ez
    real(dp) :: sx_img, sy_img, sz_img

    sx_img = sx
    sy_img = sy
    sz_img = sz
    if (axis1 == 1_i32) sx_img = sx_img + shift1
    if (axis1 == 2_i32) sy_img = sy_img + shift1
    if (axis1 == 3_i32) sz_img = sz_img + shift1
    if (axis2 == 1_i32) sx_img = sx_img + shift2
    if (axis2 == 2_i32) sy_img = sy_img + shift2
    if (axis2 == 3_i32) sz_img = sz_img + shift2
    call accumulate_point_charge_field(q, sx_img, sy_img, sz_img, tx, ty, tz, soft2, ex, ey, ez)
  end subroutine accumulate_point_charge_shifted_field

  pure subroutine accumulate_point_charge_shifted_potential( &
    q, sx, sy, sz, shift1, shift2, axis1, axis2, tx, ty, tz, soft2, phi &
    )
    real(dp), intent(in) :: q, sx, sy, sz, shift1, shift2, tx, ty, tz, soft2
    integer(i32), intent(in) :: axis1, axis2
    real(dp), intent(inout) :: phi
    real(dp) :: sx_img, sy_img, sz_img

    sx_img = sx
    sy_img = sy
    sz_img = sz
    if (axis1 == 1_i32) sx_img = sx_img + shift1
    if (axis1 == 2_i32) sy_img = sy_img + shift1
    if (axis1 == 3_i32) sz_img = sz_img + shift1
    if (axis2 == 1_i32) sx_img = sx_img + shift2
    if (axis2 == 2_i32) sy_img = sy_img + shift2
    if (axis2 == 3_i32) sz_img = sz_img + shift2
    call accumulate_point_charge_potential(q, sx_img, sy_img, sz_img, tx, ty, tz, soft2, phi)
  end subroutine accumulate_point_charge_shifted_potential

  !> 評価点を含む target leaf とその中心を返す。葉が無ければ leaf_node = 0。periodic2 の点は呼び出し側が
  !! wrap_periodic2_point で正準化済みで、展開と near 項も同じ座標を使うので、ここでは折り返さない。
  subroutine locate_target_leaf(plan, r_eval, leaf_node, leaf_center)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: r_eval(3)
    integer(i32), intent(out) :: leaf_node
    real(dp), intent(out) :: leaf_center(3)

    if (plan%target_tree_ready) then
      call descend_to_leaf( &
        plan%target_tree, check_root_extent=.true., r_eval=r_eval, leaf_node=leaf_node, leaf_center=leaf_center &
        )
    else
      call descend_to_leaf( &
        plan%source_tree, check_root_extent=.not. plan%options%use_periodic2, r_eval=r_eval, leaf_node=leaf_node, &
        leaf_center=leaf_center &
        )
      ! Source leaves bound the sources, not every point routed to their octant. M2L separation
      ! only guarantees the local expansion inside those bounds; use the direct path outside.
      if (leaf_node > 0_i32) then
        if (any(abs(r_eval - leaf_center) > plan%source_tree%node_half_size(:, leaf_node))) then
          leaf_node = 0_i32
          leaf_center = 0.0_dp
        end if
      end if
    end if
  end subroutine locate_target_leaf

  !> 木を根から評価点の八分木へ下りて葉を返す。check_root_extent なら根の箱の外の点で leaf_node = 0。
  !! 点の八分木に子が無ければ、中心が最も近い子へ進む。
  subroutine descend_to_leaf(tree, check_root_extent, r_eval, leaf_node, leaf_center)
    type(fmm_tree_type), intent(in) :: tree
    logical, intent(in) :: check_root_extent
    real(dp), intent(in) :: r_eval(3)
    integer(i32), intent(out) :: leaf_node
    real(dp), intent(out) :: leaf_center(3)
    integer(i32) :: oct, child_k, child, cand, axis
    real(dp) :: best_dist2, dist2, dx, dy, dz, extent_eps

    leaf_node = 0_i32
    leaf_center = 0.0_dp
    if (tree%nnode <= 0_i32) return

    if (check_root_extent) then
      extent_eps = 1.0d-12*max(1.0d0, maxval(abs(tree%node_center(:, 1_i32))))
      do axis = 1_i32, 3_i32
        if (abs(r_eval(axis) - tree%node_center(axis, 1_i32)) > tree%node_half_size(axis, 1_i32) + extent_eps) return
      end do
    end if

    leaf_node = 1_i32
    do while (tree%child_count(leaf_node) > 0_i32)
      oct = octant_index(r_eval(1), r_eval(2), r_eval(3), tree%node_center(:, leaf_node))
      child = 0_i32
      do child_k = 1_i32, tree%child_count(leaf_node)
        if (tree%child_octant(child_k, leaf_node) == oct) then
          child = tree%child_idx(child_k, leaf_node)
          exit
        end if
      end do
      if (child <= 0_i32) then
        child = tree%child_idx(1_i32, leaf_node)
        best_dist2 = huge(1.0d0)
        do child_k = 1_i32, tree%child_count(leaf_node)
          cand = tree%child_idx(child_k, leaf_node)
          dx = r_eval(1) - tree%node_center(1, cand)
          dy = r_eval(2) - tree%node_center(2, cand)
          dz = r_eval(3) - tree%node_center(3, cand)
          dist2 = dx*dx + dy*dy + dz*dz
          if (dist2 < best_dist2) then
            best_dist2 = dist2
            child = cand
          end if
        end do
      end if
      leaf_node = child
    end do
    leaf_center = tree%node_center(:, leaf_node)
  end subroutine descend_to_leaf

end module bem_coulomb_fmm_eval_ops

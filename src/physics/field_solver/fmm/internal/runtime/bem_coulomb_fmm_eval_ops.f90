!> Coulomb FMM 電場評価。
module bem_coulomb_fmm_eval_ops
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_types, only: fmm_plan_type, fmm_state_type
  use bem_coulomb_fmm_basis, only: build_axis_powers
  use bem_coulomb_fmm_periodic, only: wrap_periodic2_point, use_periodic2_m2l_root_oracle
  use bem_coulomb_fmm_periodic_ewald, only: add_periodic2_exact_ewald_correction_all_sources, &
                                            add_periodic2_exact_ewald_potential_correction_all_sources
  use bem_coulomb_fmm_tree_utils, only: octant_index, active_tree_nnode, active_tree_child_count, active_tree_child_idx, &
                                        active_tree_child_octant, active_tree_node_center, active_tree_node_half_size
  implicit none
  private

  public :: core_eval_points_impl
  public :: core_eval_point_impl
  public :: core_eval_potential_points_impl
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

    do i = 1_i32, ntarget
      call core_eval_point_xyz_impl( &
        plan, state, target_pos(1, i), target_pos(2, i), target_pos(3, i), e(1, i), e(2, i), e(3, i) &
        )
    end do
  end subroutine core_eval_points_impl

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

    do i = 1_i32, ntarget
      call core_eval_potential_point_xyz_impl( &
        plan, state, target_pos(1, i), target_pos(2, i), target_pos(3, i), phi(i) &
        )
    end do
  end subroutine core_eval_potential_points_impl

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
    logical :: use_periodic_root_oracle_fallback
    real(dp) :: rt(3), soft2

    ex = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    if (.not. plan%built .or. .not. state%ready) return

    rt = [rx, ry, rz]
    if (plan%options%use_periodic2) call wrap_periodic2_point(plan, rt)
    soft2 = plan%options%softening*plan%options%softening
    use_periodic_root_oracle_fallback = use_periodic2_m2l_root_oracle(plan)

    leaf_node = locate_target_leaf(plan, rt)
    leaf_slot = 0_i32
    if (leaf_node > 0_i32) leaf_slot = plan%leaf_slot_of_node(leaf_node)
    if (leaf_slot <= 0_i32) then
      call apply_direct_fallback_field( &
        plan, state, rt, soft2, use_periodic_root_oracle_fallback, ex, ey, ez &
        )
      return
    end if

    call accumulate_leaf_local_expansion(plan, state, leaf_node, rt, ex, ey, ez)
    axes = active_periodic_axes(plan)
    call accumulate_near_direct_field(plan, state, leaf_slot, rt, soft2, axes, ex, ey, ez)
  end subroutine core_eval_point_xyz_impl

  !> 1 点評価の電位本体処理を行う。
  subroutine core_eval_potential_point_xyz_impl(plan, state, rx, ry, rz, phi)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    real(dp), intent(in) :: rx, ry, rz
    real(dp), intent(out) :: phi
    integer(i32) :: leaf_node, leaf_slot
    integer(i32) :: axes(2)
    logical :: use_periodic_root_oracle_fallback
    real(dp) :: rt(3), soft2

    phi = 0.0d0
    if (.not. plan%built .or. .not. state%ready) return

    rt = [rx, ry, rz]
    if (plan%options%use_periodic2) call wrap_periodic2_point(plan, rt)
    soft2 = plan%options%softening*plan%options%softening
    use_periodic_root_oracle_fallback = use_periodic2_m2l_root_oracle(plan)

    leaf_node = locate_target_leaf(plan, rt)
    leaf_slot = 0_i32
    if (leaf_node > 0_i32) leaf_slot = plan%leaf_slot_of_node(leaf_node)
    if (leaf_slot <= 0_i32) then
      call apply_direct_fallback_potential(plan, state, rt, soft2, use_periodic_root_oracle_fallback, phi)
      return
    end if

    call accumulate_leaf_local_potential_expansion(plan, state, leaf_node, rt, phi)
    axes = active_periodic_axes(plan)
    call accumulate_near_direct_potential(plan, state, leaf_slot, rt, soft2, axes, phi)
  end subroutine core_eval_potential_point_xyz_impl

  subroutine apply_direct_fallback_field(plan, state, rt, soft2, use_periodic_root_oracle_fallback, ex, ey, ez)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: rt(3)
    real(dp), intent(in) :: soft2
    logical, intent(in) :: use_periodic_root_oracle_fallback
    real(dp), intent(inout) :: ex, ey, ez
    real(dp) :: evec(3)

    call eval_direct_all_sources_scalar(plan, state, rt(1), rt(2), rt(3), soft2, ex, ey, ez)
    if (use_periodic_root_oracle_fallback) then
      evec = [ex, ey, ez]
      call add_periodic2_exact_ewald_correction_all_sources(plan, state, rt, evec)
      ex = evec(1)
      ey = evec(2)
      ez = evec(3)
    end if
  end subroutine apply_direct_fallback_field

  subroutine apply_direct_fallback_potential(plan, state, rt, soft2, use_periodic_root_oracle_fallback, phi)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: rt(3)
    real(dp), intent(in) :: soft2
    logical, intent(in) :: use_periodic_root_oracle_fallback
    real(dp), intent(inout) :: phi

    call eval_direct_all_sources_potential_scalar(plan, state, rt(1), rt(2), rt(3), soft2, phi)
    if (use_periodic_root_oracle_fallback) then
      call add_periodic2_exact_ewald_potential_correction_all_sources(plan, state, rt, phi)
    end if
  end subroutine apply_direct_fallback_potential

  subroutine accumulate_leaf_local_expansion(plan, state, leaf_node, rt, ex, ey, ez)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    integer(i32), intent(in) :: leaf_node
    real(dp), intent(in) :: rt(3)
    real(dp), intent(inout) :: ex, ey, ez
    integer(i32) :: order, term_idx
    real(dp) :: dr(3), monomial
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    order = active_local_expansion_order(plan)
    if (order <= 0_i32) return
    if (state%local_active(leaf_node) == 0_i32) return
    if (plan%eval_term_count <= 0_i32) return

    dr = rt - active_tree_node_center(plan, plan%target_tree_ready, leaf_node)
    call build_axis_powers(dr, order, xpow, ypow, zpow)
    do term_idx = 1_i32, plan%eval_term_count
      monomial = xpow(plan%eval_exp(1, term_idx))*ypow(plan%eval_exp(2, term_idx)) &
                 *zpow(plan%eval_exp(3, term_idx))*plan%eval_inv_factorial(term_idx)
      ex = ex - state%local(plan%eval_deriv_idx(1, term_idx), leaf_node)*monomial
      ey = ey - state%local(plan%eval_deriv_idx(2, term_idx), leaf_node)*monomial
      ez = ez - state%local(plan%eval_deriv_idx(3, term_idx), leaf_node)*monomial
    end do
  end subroutine accumulate_leaf_local_expansion

  subroutine accumulate_leaf_local_potential_expansion(plan, state, leaf_node, rt, phi)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    integer(i32), intent(in) :: leaf_node
    real(dp), intent(in) :: rt(3)
    real(dp), intent(inout) :: phi
    integer(i32) :: order, alpha_idx
    real(dp) :: dr(3), monomial
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    order = active_local_expansion_order(plan)
    if (state%local_active(leaf_node) == 0_i32) return

    dr = rt - active_tree_node_center(plan, plan%target_tree_ready, leaf_node)
    call build_axis_powers(dr, order, xpow, ypow, zpow)
    do alpha_idx = 1_i32, plan%ncoef
      monomial = xpow(plan%alpha(1, alpha_idx))*ypow(plan%alpha(2, alpha_idx)) &
                 *zpow(plan%alpha(3, alpha_idx))/plan%alpha_factorial(alpha_idx)
      phi = phi + state%local(alpha_idx, leaf_node)*monomial
    end do
  end subroutine accumulate_leaf_local_potential_expansion

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

  pure subroutine accumulate_near_direct_field(plan, state, leaf_slot, rt, soft2, axes, ex, ey, ez)
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

  pure subroutine accumulate_near_direct_potential(plan, state, leaf_slot, rt, soft2, axes, phi)
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

  subroutine eval_direct_all_sources_scalar(plan, state, tx, ty, tz, soft2, ex, ey, ez)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
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
    if (plan%options%use_periodic2) then
      nshift = size(plan%shift_axis1)
      do idx = 1_i32, plan%nsrc
        call accumulate_point_charge_images_field( &
          state%src_q(idx), plan%src_pos(1, idx), plan%src_pos(2, idx), plan%src_pos(3, idx), &
          tx, ty, tz, soft2, axis1, axis2, plan%shift_axis1, plan%shift_axis2, nshift, ex, ey, ez &
          )
      end do
    else
      do idx = 1_i32, plan%nsrc
        call accumulate_point_charge_field( &
          state%src_q(idx), plan%src_pos(1, idx), plan%src_pos(2, idx), plan%src_pos(3, idx), &
          tx, ty, tz, soft2, ex, ey, ez &
          )
      end do
    end if
  end subroutine eval_direct_all_sources_scalar

  subroutine eval_direct_all_sources_potential_scalar(plan, state, tx, ty, tz, soft2, phi)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: tx, ty, tz
    real(dp), intent(in) :: soft2
    real(dp), intent(out) :: phi
    integer(i32) :: idx, axis1, axis2, nshift
    integer(i32) :: axes(2)

    phi = 0.0d0
    axes = active_periodic_axes(plan)
    axis1 = axes(1)
    axis2 = axes(2)
    if (plan%options%use_periodic2) then
      nshift = size(plan%shift_axis1)
      do idx = 1_i32, plan%nsrc
        call accumulate_point_charge_images_potential( &
          state%src_q(idx), plan%src_pos(1, idx), plan%src_pos(2, idx), plan%src_pos(3, idx), &
          tx, ty, tz, soft2, axis1, axis2, plan%shift_axis1, plan%shift_axis2, nshift, phi &
          )
      end do
    else
      do idx = 1_i32, plan%nsrc
        call accumulate_point_charge_potential( &
          state%src_q(idx), plan%src_pos(1, idx), plan%src_pos(2, idx), plan%src_pos(3, idx), &
          tx, ty, tz, soft2, phi &
          )
      end do
    end if
  end subroutine eval_direct_all_sources_potential_scalar

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

  integer(i32) function locate_target_leaf(plan, r)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: r(3)
    integer(i32) :: oct, child_k, child, cand, axis
    logical :: use_target_tree
    real(dp) :: best_dist2, dist2, dx, dy, dz, extent_eps
    real(dp) :: r_eval(3), root_center(3), root_half(3), node_center_curr(3), cand_center(3)

    locate_target_leaf = 0_i32
    if (plan%nnode <= 0_i32) return
    use_target_tree = plan%target_tree_ready
    if (active_tree_nnode(plan, use_target_tree) <= 0_i32) return

    r_eval = r
    if (use_target_tree .and. plan%options%use_periodic2) call wrap_periodic2_point(plan, r_eval)
    if (use_target_tree .or. .not. plan%options%use_periodic2) then
      root_center = active_tree_node_center(plan, use_target_tree, 1_i32)
      root_half = active_tree_node_half_size(plan, use_target_tree, 1_i32)
      extent_eps = 1.0d-12*max(1.0d0, maxval(abs(root_center)))
      do axis = 1_i32, 3_i32
        if (abs(r_eval(axis) - root_center(axis)) > root_half(axis) + extent_eps) return
      end do
    end if

    locate_target_leaf = 1_i32
    do while (active_tree_child_count(plan, use_target_tree, locate_target_leaf) > 0_i32)
      node_center_curr = active_tree_node_center(plan, use_target_tree, locate_target_leaf)
      oct = octant_index(r_eval(1), r_eval(2), r_eval(3), node_center_curr)
      child = 0_i32
      do child_k = 1_i32, active_tree_child_count(plan, use_target_tree, locate_target_leaf)
        if (active_tree_child_octant(plan, use_target_tree, child_k, locate_target_leaf) == oct) then
          child = active_tree_child_idx(plan, use_target_tree, child_k, locate_target_leaf)
          exit
        end if
      end do
      if (child <= 0_i32) then
        child = active_tree_child_idx(plan, use_target_tree, 1_i32, locate_target_leaf)
        best_dist2 = huge(1.0d0)
        do child_k = 1_i32, active_tree_child_count(plan, use_target_tree, locate_target_leaf)
          cand = active_tree_child_idx(plan, use_target_tree, child_k, locate_target_leaf)
          cand_center = active_tree_node_center(plan, use_target_tree, cand)
          dx = r_eval(1) - cand_center(1)
          dy = r_eval(2) - cand_center(2)
          dz = r_eval(3) - cand_center(3)
          dist2 = dx*dx + dy*dy + dz*dz
          if (dist2 < best_dist2) then
            best_dist2 = dist2
            child = cand
          end if
        end do
      end if
      locate_target_leaf = child
    end do
  end function locate_target_leaf

end module bem_coulomb_fmm_eval_ops

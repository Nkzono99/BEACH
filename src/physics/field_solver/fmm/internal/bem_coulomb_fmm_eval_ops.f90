!> Coulomb FMM 電場評価。
module bem_coulomb_fmm_eval_ops
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_types, only: fmm_plan_type, fmm_state_type
  use bem_coulomb_fmm_basis, only: build_axis_powers
  use bem_coulomb_fmm_periodic, only: wrap_periodic2_point, use_periodic2_ewald_like, prepare_periodic2_ewald, &
                                       add_screened_shifted_node_images, add_point_charge_images_field
  use bem_coulomb_fmm_tree_utils, only: octant_index, active_tree_nnode, active_tree_child_count, active_tree_child_idx, &
                                         active_tree_child_octant, active_tree_node_center, active_tree_node_half_size
  implicit none
  private

  public :: core_eval_points_impl
  public :: core_eval_point_impl

contains

  subroutine core_eval_points_impl(plan, state, target_pos, e)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: target_pos(:, :)
    real(dp), intent(out) :: e(:, :)
    integer(i32) :: i, ntarget

    if (size(target_pos, 1) /= 3) error stop 'FMM core expects target_pos(3,m).'
    if (size(e, 1) /= 3 .or. size(e, 2) /= size(target_pos, 2)) then
      error stop 'FMM eval_points expects e(3,m).'
    end if
    ntarget = int(size(target_pos, 2), i32)
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(plan, state, target_pos, e, ntarget) private(i)
    do i = 1_i32, ntarget
      call core_eval_point_impl(plan, state, target_pos(:, i), e(:, i))
    end do
    !$omp end parallel do
  end subroutine core_eval_points_impl

  subroutine core_eval_point_impl(plan, state, r, e)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: r(3)
    real(dp), intent(out) :: e(3)
    integer(i32) :: leaf_node, leaf_slot
    integer(i32) :: near_pos, near_node, p, idx, p_end, alpha_idx, axis, deriv_idx
    integer(i32) :: axis1, axis2, nshift
    real(dp) :: rt(3), center(3), dr(3), soft2
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order)), monomial

    e = 0.0d0
    if (.not. plan%built .or. .not. state%ready) return

    rt = r
    if (plan%options%use_periodic2) call wrap_periodic2_point(plan, rt)
    soft2 = plan%options%softening * plan%options%softening

    leaf_node = locate_target_leaf(plan, rt)
    if (leaf_node <= 0_i32) then
      call eval_direct_all_sources(plan, state, rt, soft2, e)
      if (use_periodic2_ewald_like(plan)) call add_periodic2_ewald_like_correction_all_leaves(plan, state, rt, e)
      return
    end if

    leaf_slot = plan%leaf_slot_of_node(leaf_node)
    if (leaf_slot <= 0_i32) then
      call eval_direct_all_sources(plan, state, rt, soft2, e)
      if (use_periodic2_ewald_like(plan)) call add_periodic2_ewald_like_correction_all_leaves(plan, state, rt, e)
      return
    end if

    center = active_tree_node_center(plan, plan%target_tree_ready, leaf_node)
    dr = rt - center
    call build_axis_powers(dr, plan%options%order, xpow, ypow, zpow)
    do alpha_idx = 1_i32, plan%ncoef
      if (plan%alpha_degree(alpha_idx) >= plan%options%order) cycle
      monomial = xpow(plan%alpha(1, alpha_idx)) * ypow(plan%alpha(2, alpha_idx)) * zpow(plan%alpha(3, alpha_idx)) &
                 / plan%alpha_factorial(alpha_idx)
      do axis = 1_i32, 3_i32
        deriv_idx = plan%alpha_plus_axis(axis, alpha_idx)
        if (deriv_idx <= 0_i32) cycle
        e(axis) = e(axis) - state%local(deriv_idx, leaf_node) * monomial
      end do
    end do

    axis1 = 0_i32
    axis2 = 0_i32
    if (plan%options%use_periodic2) then
      axis1 = plan%options%periodic_axes(1)
      axis2 = plan%options%periodic_axes(2)
    end if
    nshift = size(plan%shift_axis1)
    do near_pos = plan%near_start(leaf_slot), plan%near_start(leaf_slot + 1_i32) - 1_i32
      near_node = plan%near_nodes(near_pos)
      p_end = plan%node_start(near_node) + plan%node_count(near_node) - 1_i32
      do p = plan%node_start(near_node), p_end
        idx = plan%elem_order(p)
        call add_point_charge_images_field( &
          state%src_q(idx), plan%src_pos(:, idx), rt, soft2, axis1, axis2, plan%shift_axis1, plan%shift_axis2, nshift, e &
        )
      end do
    end do
    if (use_periodic2_ewald_like(plan)) call add_periodic2_ewald_like_correction(plan, state, leaf_slot, rt, e)
  end subroutine core_eval_point_impl

  subroutine eval_direct_all_sources(plan, state, target, soft2, e)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: target(3)
    real(dp), intent(in) :: soft2
    real(dp), intent(out) :: e(3)
    integer(i32) :: idx, axis1, axis2, nshift

    e = 0.0d0
    axis1 = 0_i32
    axis2 = 0_i32
    if (plan%options%use_periodic2) then
      axis1 = plan%options%periodic_axes(1)
      axis2 = plan%options%periodic_axes(2)
    end if
    nshift = size(plan%shift_axis1)
    do idx = 1_i32, plan%nsrc
      call add_point_charge_images_field( &
        state%src_q(idx), plan%src_pos(:, idx), target, soft2, axis1, axis2, plan%shift_axis1, plan%shift_axis2, nshift, e &
      )
    end do
  end subroutine eval_direct_all_sources

  subroutine add_periodic2_ewald_like_correction(plan, state, leaf_slot, r, e)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    integer(i32), intent(in) :: leaf_slot
    real(dp), intent(in) :: r(3)
    real(dp), intent(inout) :: e(3)
    integer(i32) :: node_pos, node_idx
    integer(i32) :: nimg, img_outer, axis1, axis2
    real(dp) :: alpha

    if (.not. prepare_periodic2_ewald(plan, alpha, nimg, img_outer, axis1, axis2)) return

    do node_pos = plan%far_start(leaf_slot), plan%far_start(leaf_slot + 1_i32) - 1_i32
      node_idx = plan%far_nodes(node_pos)
      call add_ewald_node_correction(plan, state, node_idx, r, axis1, axis2, nimg, img_outer, alpha, e)
    end do

    do node_pos = plan%near_start(leaf_slot), plan%near_start(leaf_slot + 1_i32) - 1_i32
      node_idx = plan%near_nodes(node_pos)
      call add_ewald_node_correction(plan, state, node_idx, r, axis1, axis2, nimg, img_outer, alpha, e)
    end do
  end subroutine add_periodic2_ewald_like_correction

  subroutine add_periodic2_ewald_like_correction_all_leaves(plan, state, r, e)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: r(3)
    real(dp), intent(inout) :: e(3)
    integer(i32) :: node_idx
    integer(i32) :: nimg, img_outer, axis1, axis2
    real(dp) :: alpha

    if (.not. prepare_periodic2_ewald(plan, alpha, nimg, img_outer, axis1, axis2)) return

    do node_idx = 1_i32, plan%nnode
      if (plan%child_count(node_idx) > 0_i32) cycle
      call add_ewald_node_correction(plan, state, node_idx, r, axis1, axis2, nimg, img_outer, alpha, e)
    end do
  end subroutine add_periodic2_ewald_like_correction_all_leaves

  subroutine add_ewald_node_correction(plan, state, node_idx, r, axis1, axis2, nimg, img_outer, alpha, e)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    integer(i32), intent(in) :: node_idx
    real(dp), intent(in) :: r(3)
    integer(i32), intent(in) :: axis1, axis2, nimg, img_outer
    real(dp), intent(in) :: alpha
    real(dp), intent(inout) :: e(3)
    real(dp) :: q, src(3)

    call recover_node_charge_center(plan, state, node_idx, q, src)
    if (abs(q) <= tiny(1.0d0)) return
    call add_screened_shifted_node_images(q, src, r, axis1, axis2, plan%options%periodic_len, nimg, img_outer, alpha, e)
  end subroutine add_ewald_node_correction

  subroutine recover_node_charge_center(plan, state, node_idx, q, src)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    integer(i32), intent(in) :: node_idx
    real(dp), intent(out) :: q
    real(dp), intent(out) :: src(3)
    integer(i32) :: idx_x, idx_y, idx_z, idx_q

    src = plan%node_center(:, node_idx)
    idx_q = plan%alpha_map(0_i32, 0_i32, 0_i32)
    q = state%multipole(idx_q, node_idx)
    if (abs(q) <= tiny(1.0d0)) return
    if (plan%options%order < 1_i32) return

    idx_x = plan%alpha_map(1_i32, 0_i32, 0_i32)
    idx_y = plan%alpha_map(0_i32, 1_i32, 0_i32)
    idx_z = plan%alpha_map(0_i32, 0_i32, 1_i32)
    src(1) = src(1) + state%multipole(idx_x, node_idx) / q
    src(2) = src(2) + state%multipole(idx_y, node_idx) / q
    src(3) = src(3) + state%multipole(idx_z, node_idx) / q
  end subroutine recover_node_charge_center

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
      extent_eps = 1.0d-12 * max(1.0d0, maxval(abs(root_center)))
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
          dist2 = dx * dx + dy * dy + dz * dz
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

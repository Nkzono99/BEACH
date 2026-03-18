!> Coulomb FMM 電場評価。
module bem_coulomb_fmm_eval_ops
    use bem_kinds, only: dp, i32, i64
    use bem_coulomb_fmm_types, only: fmm_plan_type, fmm_state_type
    use bem_coulomb_fmm_basis, only: build_axis_powers
    use bem_coulomb_fmm_periodic, only: wrap_periodic2_point, use_periodic2_ewald, use_periodic2_exact_ewald, &
                                        prepare_periodic2_ewald, add_screened_shifted_node_images, &
                                        add_exact_periodic2_real_space_source_correction
    use bem_coulomb_fmm_tree_utils, only: octant_index, active_tree_nnode, active_tree_child_count, active_tree_child_idx, &
                                          active_tree_child_octant, active_tree_node_center, active_tree_node_half_size
    use bem_performance_profile, only: perf_wall_time_seconds
    implicit none
    private

    public :: core_eval_points_impl
    public :: core_eval_point_impl

contains

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
    !!$omp parallel do default(none) schedule(static) &
    !!$omp shared(plan, state, target_pos, e, ntarget) private(i)
        do i = 1_i32, ntarget
            call core_eval_point_xyz_impl( &
                plan, state, target_pos(1, i), target_pos(2, i), target_pos(3, i), e(1, i), e(2, i), e(3, i) &
                )
        end do
    !!$omp end parallel do
    end subroutine core_eval_points_impl

    subroutine core_eval_point_impl(plan, state, r, e)
        type(fmm_plan_type), intent(in) :: plan
        type(fmm_state_type), intent(inout) :: state
        real(dp), intent(in) :: r(3)
        real(dp), intent(out) :: e(3)

        call core_eval_point_xyz_impl(plan, state, r(1), r(2), r(3), e(1), e(2), e(3))
    end subroutine core_eval_point_impl

    subroutine core_eval_point_xyz_impl(plan, state, rx, ry, rz, ex, ey, ez)
        type(fmm_plan_type), intent(in) :: plan
        type(fmm_state_type), intent(inout) :: state
        real(dp), intent(in) :: rx, ry, rz
        real(dp), intent(out) :: ex, ey, ez
        integer(i32) :: leaf_node, leaf_slot
        integer(i32) :: near_pos, idx, term_idx
        integer(i32) :: axis1, axis2, nshift, order, near_source_count_i32
        integer(i32) :: near_source_begin, near_source_end
        integer(i64) :: near_source_count_local, direct_kernel_count_local
        integer(i32) :: eval_count_local, local_count_local, fallback_count_local, ewald_count_local
        logical :: profile_enabled, use_ewald, use_exact_ewald
        real(dp) :: rt(3), dr(3), soft2, monomial, e_arr(3)
        real(dp) :: shift1, shift2
        real(dp) :: t0, locate_time_local, local_time_local, near_time_local, fallback_time_local, ewald_time_local
        real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
        real(dp) :: zpow(0:max(0_i32, plan%options%order))

        ex = 0.0d0
        ey = 0.0d0
        ez = 0.0d0
        if (.not. plan%built .or. .not. state%ready) return

        profile_enabled = state%profile_enabled
        eval_count_local = 1_i32
        local_count_local = 0_i32
        fallback_count_local = 0_i32
        ewald_count_local = 0_i32
        near_source_count_local = 0_i64
        direct_kernel_count_local = 0_i64
        locate_time_local = 0.0d0
        local_time_local = 0.0d0
        near_time_local = 0.0d0
        fallback_time_local = 0.0d0
        ewald_time_local = 0.0d0

        rt = [rx, ry, rz]
        if (plan%options%use_periodic2) call wrap_periodic2_point(plan, rt)
        soft2 = plan%options%softening*plan%options%softening

        use_ewald = use_periodic2_ewald(plan)
        use_exact_ewald = use_periodic2_exact_ewald(plan)

        if (profile_enabled) t0 = perf_wall_time_seconds()

        leaf_node = locate_target_leaf(plan, rt)
        if (profile_enabled) locate_time_local = perf_wall_time_seconds() - t0
        if (leaf_node <= 0_i32) then
            fallback_count_local = 1_i32
            direct_kernel_count_local = direct_kernel_count_local + estimate_direct_kernel_count(plan)
            if (profile_enabled) t0 = perf_wall_time_seconds()
            call eval_direct_all_sources_scalar(plan, state, rt(1), rt(2), rt(3), soft2, ex, ey, ez)
            if (profile_enabled) fallback_time_local = perf_wall_time_seconds() - t0
            if (use_ewald) then
                ewald_count_local = 1_i32
                if (profile_enabled) t0 = perf_wall_time_seconds()
                e_arr = [ex, ey, ez]
                if (use_exact_ewald) then
                    call add_periodic2_exact_ewald_correction(plan, state, rt, e_arr)
                else
                    call add_periodic2_ewald_like_correction_all_leaves(plan, state, rt, e_arr)
                end if
                ex = e_arr(1)
                ey = e_arr(2)
                ez = e_arr(3)
                if (profile_enabled) ewald_time_local = perf_wall_time_seconds() - t0
            end if
            call record_eval_profile( &
                state, eval_count_local, local_count_local, fallback_count_local, ewald_count_local, &
                near_source_count_local, direct_kernel_count_local, locate_time_local, local_time_local, &
                near_time_local, fallback_time_local, ewald_time_local &
                )
            return
        end if

        leaf_slot = plan%leaf_slot_of_node(leaf_node)
        if (leaf_slot <= 0_i32) then
            fallback_count_local = 1_i32
            direct_kernel_count_local = direct_kernel_count_local + estimate_direct_kernel_count(plan)
            if (profile_enabled) t0 = perf_wall_time_seconds()
            call eval_direct_all_sources_scalar(plan, state, rt(1), rt(2), rt(3), soft2, ex, ey, ez)
            if (profile_enabled) fallback_time_local = perf_wall_time_seconds() - t0
            if (use_ewald) then
                ewald_count_local = 1_i32
                if (profile_enabled) t0 = perf_wall_time_seconds()
                e_arr = [ex, ey, ez]
                if (use_exact_ewald) then
                    call add_periodic2_exact_ewald_correction(plan, state, rt, e_arr)
                else
                    call add_periodic2_ewald_like_correction_all_leaves(plan, state, rt, e_arr)
                end if
                ex = e_arr(1)
                ey = e_arr(2)
                ez = e_arr(3)
                if (profile_enabled) ewald_time_local = perf_wall_time_seconds() - t0
            end if
            call record_eval_profile( &
                state, eval_count_local, local_count_local, fallback_count_local, ewald_count_local, &
                near_source_count_local, direct_kernel_count_local, locate_time_local, local_time_local, &
                near_time_local, fallback_time_local, ewald_time_local &
                )
            return
        end if

        order = plan%options%order
        if (order > 0_i32 .and. state%local_active(leaf_node) /= 0_i32 .and. plan%eval_term_count > 0_i32) then
            local_count_local = 1_i32
            if (profile_enabled) t0 = perf_wall_time_seconds()
            dr = rt - active_tree_node_center(plan, plan%target_tree_ready, leaf_node)
            call build_axis_powers(dr, order, xpow, ypow, zpow)
            do term_idx = 1_i32, plan%eval_term_count
                monomial = xpow(plan%eval_exp(1, term_idx))*ypow(plan%eval_exp(2, term_idx)) &
                           *zpow(plan%eval_exp(3, term_idx))*plan%eval_inv_factorial(term_idx)
                ex = ex - state%local(plan%eval_deriv_idx(1, term_idx), leaf_node)*monomial
                ey = ey - state%local(plan%eval_deriv_idx(2, term_idx), leaf_node)*monomial
                ez = ez - state%local(plan%eval_deriv_idx(3, term_idx), leaf_node)*monomial
            end do
            if (profile_enabled) local_time_local = perf_wall_time_seconds() - t0
        end if

        axis1 = 0_i32
        axis2 = 0_i32
        if (plan%options%use_periodic2) then
            axis1 = plan%options%periodic_axes(1)
            axis2 = plan%options%periodic_axes(2)
        end if
        near_source_begin = plan%near_source_start(leaf_slot)
        near_source_end = plan%near_source_start(leaf_slot + 1_i32) - 1_i32
        if (near_source_end >= near_source_begin) then
            near_source_count_i32 = near_source_end - near_source_begin + 1_i32
            near_source_count_local = near_source_count_local + int(near_source_count_i32, i64)
            direct_kernel_count_local = direct_kernel_count_local + int(near_source_count_i32, i64)
            if (profile_enabled) t0 = perf_wall_time_seconds()
            if (plan%options%use_periodic2) then
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
            if (profile_enabled) near_time_local = perf_wall_time_seconds() - t0
        end if

        if (use_ewald) then
            ewald_count_local = 1_i32
            if (profile_enabled) t0 = perf_wall_time_seconds()
            e_arr = [ex, ey, ez]
            if (use_exact_ewald) then
                call add_periodic2_exact_ewald_correction(plan, state, rt, e_arr)
            else
                call add_periodic2_ewald_like_correction(plan, state, leaf_slot, rt, e_arr)
            end if
            ex = e_arr(1)
            ey = e_arr(2)
            ez = e_arr(3)
            if (profile_enabled) ewald_time_local = perf_wall_time_seconds() - t0
        end if

        call record_eval_profile( &
            state, eval_count_local, local_count_local, fallback_count_local, ewald_count_local, &
            near_source_count_local, direct_kernel_count_local, locate_time_local, local_time_local, &
            near_time_local, fallback_time_local, ewald_time_local &
            )
    end subroutine core_eval_point_xyz_impl

    subroutine eval_direct_all_sources_scalar(plan, state, tx, ty, tz, soft2, ex, ey, ez)
        type(fmm_plan_type), intent(in) :: plan
        type(fmm_state_type), intent(in) :: state
        real(dp), intent(in) :: tx, ty, tz
        real(dp), intent(in) :: soft2
        real(dp), intent(out) :: ex, ey, ez
        integer(i32) :: idx, axis1, axis2, nshift

        ex = 0.0d0
        ey = 0.0d0
        ez = 0.0d0
        axis1 = 0_i32
        axis2 = 0_i32
        if (plan%options%use_periodic2) then
            axis1 = plan%options%periodic_axes(1)
            axis2 = plan%options%periodic_axes(2)
        end if
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

    subroutine add_periodic2_ewald_like_correction(plan, state, leaf_slot, r, e)
        type(fmm_plan_type), intent(in) :: plan
        type(fmm_state_type), intent(in) :: state
        integer(i32), intent(in) :: leaf_slot
        real(dp), intent(in) :: r(3)
        real(dp), intent(inout) :: e(3)
        integer(i32) :: node_pos, node_idx
        integer(i32) :: nimg, img_outer, kmax, axis1, axis2, axis_free
        real(dp) :: cell_area
        logical :: use_like, use_exact
        real(dp) :: alpha

        if (.not. prepare_periodic2_ewald( &
            plan, alpha, nimg, img_outer, kmax, axis1, axis2, axis_free, cell_area, use_like, use_exact &
            )) return
        if (.not. use_like) return

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
        integer(i32) :: nimg, img_outer, kmax, axis1, axis2, axis_free
        real(dp) :: cell_area
        logical :: use_like, use_exact
        real(dp) :: alpha

        if (.not. prepare_periodic2_ewald( &
            plan, alpha, nimg, img_outer, kmax, axis1, axis2, axis_free, cell_area, use_like, use_exact &
            )) return
        if (.not. use_like) return

        do node_idx = 1_i32, plan%nnode
            if (plan%child_count(node_idx) > 0_i32) cycle
            call add_ewald_node_correction(plan, state, node_idx, r, axis1, axis2, nimg, img_outer, alpha, e)
        end do
    end subroutine add_periodic2_ewald_like_correction_all_leaves

    subroutine add_periodic2_exact_ewald_correction(plan, state, r, e)
        type(fmm_plan_type), intent(in) :: plan
        type(fmm_state_type), intent(in) :: state
        real(dp), intent(in) :: r(3)
        real(dp), intent(inout) :: e(3)
        integer(i32) :: idx

        if (.not. plan%exact_ewald_ready) return

        do idx = 1_i32, plan%nsrc
            call add_exact_periodic2_real_space_source_correction( &
                state%src_q(idx), plan%src_pos(:, idx), r, plan%exact_ewald_soft2, &
                plan%exact_ewald_axis1, plan%exact_ewald_axis2, &
                plan%exact_ewald_screen_shift1, plan%exact_ewald_screen_shift2, &
                plan%exact_ewald_screen_count, &
                plan%exact_ewald_inner_shift1, plan%exact_ewald_inner_shift2, &
                plan%exact_ewald_inner_count, plan%exact_ewald_alpha, e &
                )
        end do
        call add_exact_periodic2_reciprocal_correction(plan, state, r, e)
        call add_exact_periodic2_k0_correction_all_sources(plan, state, r, e)
    end subroutine add_periodic2_exact_ewald_correction

    subroutine add_exact_periodic2_reciprocal_correction(plan, state, r, e)
        type(fmm_plan_type), intent(in) :: plan
        type(fmm_state_type), intent(in) :: state
        real(dp), intent(in) :: r(3)
        real(dp), intent(inout) :: e(3)
        integer(i32) :: k_idx, idx
        real(dp) :: rt1, rt2, zt, alpha_zt, theta_t, sin_t, cos_t
        real(dp) :: dz, arg_p, arg_m, term_p, term_m, pair_sum, pair_diff
        real(dp) :: sum_pair_cos, sum_pair_sin, sum_diff_cos, sum_diff_sin
        real(dp) :: sin_acc, cos_acc, kmag, arg0

        if (plan%exact_ewald_k_count <= 0_i32) return

        rt1 = r(plan%exact_ewald_axis1)
        rt2 = r(plan%exact_ewald_axis2)
        zt = r(plan%exact_ewald_axis_free)
        alpha_zt = plan%exact_ewald_alpha*zt
        do k_idx = 1_i32, plan%exact_ewald_k_count
            kmag = plan%exact_ewald_kmag(k_idx)
            arg0 = plan%exact_ewald_karg0(k_idx)
            sum_pair_cos = 0.0d0
            sum_pair_sin = 0.0d0
            sum_diff_cos = 0.0d0
            sum_diff_sin = 0.0d0
            do idx = 1_i32, plan%nsrc
                dz = zt - plan%exact_ewald_src_free(idx)
                arg_p = arg0 + alpha_zt - plan%exact_ewald_src_alpha_free(idx)
                arg_m = arg0 - alpha_zt + plan%exact_ewald_src_alpha_free(idx)
                term_p = exp(kmag*dz)*erfc(arg_p)
                term_m = exp(-kmag*dz)*erfc(arg_m)
                pair_sum = term_p + term_m
                pair_diff = term_m - term_p
                sum_pair_cos = sum_pair_cos + state%exact_ewald_qcos(idx, k_idx)*pair_sum
                sum_pair_sin = sum_pair_sin + state%exact_ewald_qsin(idx, k_idx)*pair_sum
                sum_diff_cos = sum_diff_cos + state%exact_ewald_qcos(idx, k_idx)*pair_diff
                sum_diff_sin = sum_diff_sin + state%exact_ewald_qsin(idx, k_idx)*pair_diff
            end do
            theta_t = plan%exact_ewald_k1(k_idx)*rt1 + plan%exact_ewald_k2(k_idx)*rt2
            sin_t = sin(theta_t)
            cos_t = cos(theta_t)
            sin_acc = sin_t*sum_pair_cos - cos_t*sum_pair_sin
            cos_acc = cos_t*sum_diff_cos + sin_t*sum_diff_sin
            e(plan%exact_ewald_axis1) = e(plan%exact_ewald_axis1) + plan%exact_ewald_kpref1(k_idx)*sin_acc
            e(plan%exact_ewald_axis2) = e(plan%exact_ewald_axis2) + plan%exact_ewald_kpref2(k_idx)*sin_acc
            e(plan%exact_ewald_axis_free) = e(plan%exact_ewald_axis_free) + plan%exact_ewald_kprefz(k_idx)*cos_acc
        end do
    end subroutine add_exact_periodic2_reciprocal_correction

    subroutine add_exact_periodic2_k0_correction_all_sources(plan, state, r, e)
        type(fmm_plan_type), intent(in) :: plan
        type(fmm_state_type), intent(in) :: state
        real(dp), intent(in) :: r(3)
        real(dp), intent(inout) :: e(3)
        integer(i32) :: idx
        real(dp) :: alpha_zt

        if (plan%nsrc <= 0_i32) return
        alpha_zt = plan%exact_ewald_alpha*r(plan%exact_ewald_axis_free)
        do idx = 1_i32, plan%nsrc
            e(plan%exact_ewald_axis_free) = e(plan%exact_ewald_axis_free) + plan%exact_ewald_k0_pref*state%src_q(idx) &
                                            *erf(alpha_zt - plan%exact_ewald_src_alpha_free(idx))
        end do
    end subroutine add_exact_periodic2_k0_correction_all_sources

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
        src(1) = src(1) + state%multipole(idx_x, node_idx)/q
        src(2) = src(2) + state%multipole(idx_y, node_idx)/q
        src(3) = src(3) + state%multipole(idx_z, node_idx)/q
    end subroutine recover_node_charge_center

    integer(i64) function estimate_direct_kernel_count(plan)
        type(fmm_plan_type), intent(in) :: plan
        integer(i64) :: nshift

        estimate_direct_kernel_count = int(plan%nsrc, i64)
        if (plan%options%use_periodic2) then
            nshift = int(size(plan%shift_axis1), i64)
            estimate_direct_kernel_count = estimate_direct_kernel_count*nshift*nshift
        end if
    end function estimate_direct_kernel_count

    subroutine record_eval_profile( &
        state, eval_count, local_count, fallback_count, ewald_count, near_source_count, direct_kernel_count, &
        locate_time_s, local_time_s, near_time_s, fallback_time_s, ewald_time_s &
        )
        type(fmm_state_type), intent(inout) :: state
        integer(i32), intent(in) :: eval_count, local_count, fallback_count, ewald_count
        integer(i64), intent(in) :: near_source_count, direct_kernel_count
        real(dp), intent(in) :: locate_time_s, local_time_s, near_time_s, fallback_time_s, ewald_time_s

        if (.not. state%profile_enabled) return

        !$omp atomic update
        state%eval_count = state%eval_count + eval_count
        !$omp atomic update
        state%eval_local_count = state%eval_local_count + local_count
        !$omp atomic update
        state%eval_fallback_count = state%eval_fallback_count + fallback_count
        !$omp atomic update
        state%eval_ewald_count = state%eval_ewald_count + ewald_count
        !$omp atomic update
        state%eval_near_source_count = state%eval_near_source_count + near_source_count
        !$omp atomic update
        state%eval_direct_kernel_count = state%eval_direct_kernel_count + direct_kernel_count
        !$omp atomic update
        state%eval_locate_time_s = state%eval_locate_time_s + locate_time_s
        !$omp atomic update
        state%eval_local_time_s = state%eval_local_time_s + local_time_s
        !$omp atomic update
        state%eval_near_time_s = state%eval_near_time_s + near_time_s
        !$omp atomic update
        state%eval_fallback_time_s = state%eval_fallback_time_s + fallback_time_s
        !$omp atomic update
        state%eval_ewald_time_s = state%eval_ewald_time_s + ewald_time_s
    end subroutine record_eval_profile

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

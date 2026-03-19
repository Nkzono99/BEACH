!> Coulomb FMM state 更新と upward/downward pass。
module bem_coulomb_fmm_state_ops
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_types, only: fmm_plan_type, fmm_state_type, reset_fmm_state
  use bem_coulomb_fmm_tree_utils, only: active_tree_nnode
  implicit none
  private

  public :: core_update_state_impl
  public :: core_destroy_state_impl

contains

  subroutine core_update_state_impl(plan, state, src_q)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    real(dp), intent(in) :: src_q(:)
    integer(i32) :: n_target_nodes

    if (.not. plan%built) error stop 'FMM plan is not built.'
    if (size(src_q) /= plan%nsrc) error stop 'src_q size does not match plan.'

    state%ready = .false.
    call ensure_state_capacity(plan, state)
    n_target_nodes = active_tree_nnode(plan, plan%target_tree_ready)
    !$omp parallel default(none) shared(plan, state, src_q, n_target_nodes)
    if (allocated(state%multipole_active)) then
      !$omp workshare
      state%multipole_active = 0_i32
      !$omp end workshare
    end if
    if (allocated(state%local_active)) then
      !$omp workshare
      state%local_active = 0_i32
      !$omp end workshare
    end if
    if (plan%nsrc > 0_i32) then
      !$omp workshare
      state%src_q = src_q
      !$omp end workshare
    end if
    if (plan%nsource_leaf <= 0_i32 .and. allocated(state%multipole)) then
      !$omp workshare
      state%multipole = 0.0d0
      !$omp end workshare
    end if
    if (plan%m2l_pair_count <= 0_i32 .and. allocated(state%local)) then
      !$omp workshare
      state%local = 0.0d0
      !$omp end workshare
    end if

    call p2m_leaf_moments(plan, state)
    call m2m_upward_pass(plan, state)
    if (n_target_nodes > 0_i32) then
      call m2l_accumulate(plan, state)
      call apply_periodic_root_trunc_correction(plan, state)
      call l2l_downward_pass(plan, state)
    end if
    !$omp end parallel

    state%ready = .true.
    state%update_count = state%update_count + 1_i32
  end subroutine core_update_state_impl

  subroutine core_destroy_state_impl(state)
    type(fmm_state_type), intent(inout) :: state

    call reset_fmm_state(state)
  end subroutine core_destroy_state_impl

  subroutine ensure_state_capacity(plan, state)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    integer(i32) :: n_target_nodes

    n_target_nodes = active_tree_nnode(plan, plan%target_tree_ready)
    if (allocated(state%src_q)) then
      if (size(state%src_q) /= plan%nsrc) deallocate (state%src_q)
    end if
    if (.not. allocated(state%src_q)) allocate (state%src_q(plan%nsrc))

    if (allocated(state%multipole)) then
      if (size(state%multipole, 1) /= plan%ncoef .or. size(state%multipole, 2) /= max(1_i32, plan%nnode)) then
        deallocate (state%multipole)
      end if
    end if
    if (.not. allocated(state%multipole)) allocate (state%multipole(plan%ncoef, max(1_i32, plan%nnode)))
    if (allocated(state%multipole_active)) then
      if (size(state%multipole_active) /= max(1_i32, plan%nnode)) deallocate (state%multipole_active)
    end if
    if (.not. allocated(state%multipole_active)) allocate (state%multipole_active(max(1_i32, plan%nnode)))

    if (allocated(state%local)) then
      if (size(state%local, 1) /= plan%ncoef .or. size(state%local, 2) /= max(1_i32, n_target_nodes)) then
        deallocate (state%local)
      end if
    end if
    if (.not. allocated(state%local)) allocate (state%local(plan%ncoef, max(1_i32, n_target_nodes)))
    if (allocated(state%local_active)) then
      if (size(state%local_active) /= max(1_i32, n_target_nodes)) deallocate (state%local_active)
    end if
    if (.not. allocated(state%local_active)) allocate (state%local_active(max(1_i32, n_target_nodes)))
  end subroutine ensure_state_capacity

  subroutine p2m_leaf_moments(plan, state)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    integer(i32) :: leaf_idx, node_idx, p, idx, p_end, alpha_idx
    real(dp) :: leaf_multipole(plan%ncoef)
    logical :: leaf_active

    if (plan%nsource_leaf <= 0_i32) return
    !$omp do schedule(static)
    do leaf_idx = 1_i32, plan%nsource_leaf
      node_idx = plan%source_leaf_nodes(leaf_idx)
      leaf_multipole = 0.0d0
      leaf_active = .false.
      p_end = plan%node_start(node_idx) + plan%node_count(node_idx) - 1_i32
      do p = plan%node_start(node_idx), p_end
        idx = plan%elem_order(p)
        if (abs(state%src_q(idx)) <= tiny(1.0d0)) cycle
        do alpha_idx = 1_i32, plan%ncoef
          leaf_multipole(alpha_idx) = leaf_multipole(alpha_idx) + state%src_q(idx) * plan%source_p2m_basis(alpha_idx, p)
        end do
      end do
      leaf_active = coefficient_vector_is_active(leaf_multipole)
      state%multipole(:, node_idx) = leaf_multipole
      state%multipole_active(node_idx) = merge(1_i32, 0_i32, leaf_active)
    end do
    !$omp end do
  end subroutine p2m_leaf_moments

  subroutine m2m_upward_pass(plan, state)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    integer(i32) :: depth, level_pos, level_start_pos, level_end_pos
    integer(i32) :: node_idx, child_k, child_node, beta_idx, term_idx, alpha_idx, delta_idx
    real(dp) :: node_multipole(plan%ncoef)
    logical :: node_active

    if (plan%node_max_depth <= 0_i32) return
    do depth = plan%node_max_depth - 1_i32, 0_i32, -1_i32
      level_start_pos = plan%node_level_start(depth + 1_i32)
      level_end_pos = plan%node_level_start(depth + 2_i32) - 1_i32
      !$omp do schedule(static)
      do level_pos = level_start_pos, level_end_pos
        node_idx = plan%node_level_nodes(level_pos)
        if (plan%child_count(node_idx) <= 0_i32) cycle
        node_multipole = 0.0d0
        do child_k = 1_i32, plan%child_count(node_idx)
          child_node = plan%child_idx(child_k, node_idx)
          if (state%multipole_active(child_node) == 0_i32) cycle
          do beta_idx = 1_i32, plan%ncoef
            do term_idx = 1_i32, plan%m2m_term_count(beta_idx)
              alpha_idx = plan%m2m_alpha_list(term_idx, beta_idx)
              delta_idx = plan%m2m_delta_list(term_idx, beta_idx)
              node_multipole(beta_idx) = node_multipole(beta_idx) &
                                          + state%multipole(alpha_idx, child_node) &
                                          * plan%source_shift_monomial(delta_idx, child_node)
            end do
          end do
        end do
        node_active = coefficient_vector_is_active(node_multipole)
        state%multipole(:, node_idx) = node_multipole
        state%multipole_active(node_idx) = merge(1_i32, 0_i32, node_active)
      end do
      !$omp end do
    end do
  end subroutine m2m_upward_pass

  subroutine m2l_accumulate(plan, state)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    integer(i32) :: t_node, pair_pos, pair_idx, s_node, alpha_idx, beta_idx, deriv_idx, n_target_nodes
    real(dp) :: local_acc(plan%ncoef), source_coeff
    logical :: local_active

    if (plan%m2l_pair_count <= 0_i32) return
    n_target_nodes = active_tree_nnode(plan, plan%target_tree_ready)
    !$omp do schedule(static)
    do t_node = 1_i32, n_target_nodes
      local_acc = 0.0d0
      local_active = .false.
      do pair_pos = plan%m2l_target_start(t_node), plan%m2l_target_start(t_node + 1_i32) - 1_i32
        pair_idx = plan%m2l_pair_order(pair_pos)
        s_node = plan%m2l_source_nodes(pair_idx)
        if (state%multipole_active(s_node) == 0_i32) cycle
        do beta_idx = 1_i32, plan%ncoef
          source_coeff = plan%alpha_sign(beta_idx) * state%multipole(beta_idx, s_node)
          if (abs(source_coeff) <= tiny(1.0d0)) cycle
          do alpha_idx = 1_i32, plan%ncoef
            deriv_idx = plan%alpha_beta_deriv_idx(alpha_idx, beta_idx)
            local_acc(alpha_idx) = local_acc(alpha_idx) + source_coeff * plan%m2l_deriv(deriv_idx, pair_idx)
          end do
        end do
      end do
      local_active = coefficient_vector_is_active(local_acc)
      state%local(:, t_node) = local_acc
      state%local_active(t_node) = merge(1_i32, 0_i32, local_active)
    end do
    !$omp end do
  end subroutine m2l_accumulate

  subroutine apply_periodic_root_trunc_correction(plan, state)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    real(dp) :: root_local(plan%ncoef)
    logical :: root_active

    if (.not. plan%periodic_root_trunc_operator_ready) return
    if (size(state%local, 2) <= 0) return
    !$omp single
    root_local = state%local(:, 1_i32)
    if (plan%nnode > 0_i32 .and. state%multipole_active(1_i32) /= 0_i32) then
      root_local = root_local + matmul(plan%periodic_root_trunc_operator, state%multipole(:, 1_i32))
    end if
    root_active = coefficient_vector_is_active(root_local)
    state%local(:, 1_i32) = root_local
    state%local_active(1_i32) = merge(1_i32, 0_i32, root_active)
    !$omp end single
  end subroutine apply_periodic_root_trunc_correction

  subroutine l2l_downward_pass(plan, state)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    integer(i32) :: depth, max_depth, level_pos, level_start_pos, level_end_pos
    integer(i32) :: node_idx, parent_node, alpha_idx, term_idx, gamma_idx, delta_idx
    real(dp) :: parent_local(plan%ncoef), local_acc(plan%ncoef)
    logical :: use_target_tree
    logical :: parent_active, node_active

    use_target_tree = plan%target_tree_ready
    if (use_target_tree) then
      max_depth = plan%target_node_max_depth
    else
      max_depth = plan%node_max_depth
    end if
    if (max_depth <= 0_i32) return

    do depth = 1_i32, max_depth
      if (use_target_tree) then
        level_start_pos = plan%target_level_start(depth + 1_i32)
        level_end_pos = plan%target_level_start(depth + 2_i32) - 1_i32
      else
        level_start_pos = plan%node_level_start(depth + 1_i32)
        level_end_pos = plan%node_level_start(depth + 2_i32) - 1_i32
      end if
      !$omp do schedule(static)
      do level_pos = level_start_pos, level_end_pos
        if (use_target_tree) then
          node_idx = plan%target_level_nodes(level_pos)
        else
          node_idx = plan%node_level_nodes(level_pos)
        end if
        parent_node = plan%parent_of(node_idx)
        if (parent_node <= 0_i32) cycle
        parent_active = state%local_active(parent_node) /= 0_i32
        node_active = state%local_active(node_idx) /= 0_i32
        if (.not. parent_active .and. .not. node_active) cycle
        parent_local = state%local(:, parent_node)
        local_acc = state%local(:, node_idx)
        if (parent_active) then
          do alpha_idx = 1_i32, plan%ncoef
            do term_idx = 1_i32, plan%l2l_term_count(alpha_idx)
              gamma_idx = plan%l2l_gamma_list(term_idx, alpha_idx)
              delta_idx = plan%l2l_delta_list(term_idx, alpha_idx)
              local_acc(alpha_idx) = local_acc(alpha_idx) + parent_local(gamma_idx) &
                                     * plan%target_shift_monomial(delta_idx, node_idx)
            end do
          end do
        end if
        node_active = coefficient_vector_is_active(local_acc)
        state%local(:, node_idx) = local_acc
        state%local_active(node_idx) = merge(1_i32, 0_i32, node_active)
      end do
      !$omp end do
    end do
  end subroutine l2l_downward_pass

  pure logical function coefficient_vector_is_active(coeff)
    real(dp), intent(in) :: coeff(:)

    coefficient_vector_is_active = any(coeff /= 0.0d0)
  end function coefficient_vector_is_active

end module bem_coulomb_fmm_state_ops

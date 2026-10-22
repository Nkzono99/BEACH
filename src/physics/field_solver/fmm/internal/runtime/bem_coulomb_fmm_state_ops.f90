!> Coulomb FMM state 更新と upward/downward pass。
module bem_coulomb_fmm_state_ops
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: eps0
  use bem_coulomb_fmm_types, only: fmm_plan_type, fmm_state_type, initialize_fmm_state, reset_fmm_state
  use bem_coulomb_fmm_tree_utils, only: active_tree_nnode, active_tree_max_depth, active_tree_level_bounds, &
                                        active_tree_level_node
  use bem_periodic_zero_mode_plan, only: refresh_periodic_zero_mode_state
  use bem_periodic_shell_mean, only: refresh_periodic_shell_mean_state
  use bem_mpi, only: mpi_context, mpi_get_rank_size, mpi_allgatherv_inplace_real_dp, mpi_allgatherv_inplace_i32
  implicit none
  private

  public :: core_update_state_impl
  public :: core_destroy_state_impl

contains

  !> ソース電荷から FMM state を更新する。
  !! @param[in] plan 構築済みの FMM 計画。
  !! @param[inout] state 更新対象の state。
  !! @param[in] src_q ソース電荷量。
  !! @param[in] mpi 与えると M2L の target node を rank で分担し、局所展開を allgather で揃える。
  subroutine core_update_state_impl(plan, state, src_q, mpi)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    real(dp), intent(in) :: src_q(:)
    type(mpi_context), intent(in), optional :: mpi
    integer(i32) :: n_target_nodes, source_idx, cell, rank, rank_count, block_begin, block_end
    integer(i32), allocatable :: block_counts(:), block_offsets(:)
    real(dp) :: total_charge, e_bottom, phi_gauge, z_gauge

    if (.not. plan%built) error stop 'FMM plan is not built.'
    if (size(src_q) /= plan%nsrc) error stop 'src_q size does not match plan.'

    state%ready = .false.
    call ensure_state_capacity(plan, state)
    if (plan%periodic_k0_ready) then
      total_charge = sum(src_q)
      z_gauge = plan%periodic_k0_plan%break_z(1)
      e_bottom = -0.5_dp*total_charge/(eps0*plan%periodic_k0_plan%area_xy)
      phi_gauge = 0.0_dp
      do source_idx = 1_i32, plan%nsrc
        phi_gauge = phi_gauge - 0.5_dp*src_q(source_idx)* &
                    (sum(plan%periodic_k0_plan%panel(source_idx)%z)/3.0_dp - z_gauge)/ &
                    (eps0*plan%periodic_k0_plan%area_xy)
      end do
      call refresh_periodic_zero_mode_state( &
        plan%periodic_k0_plan, src_q, e_bottom, z_gauge, phi_gauge, state%periodic_k0_state &
        )
      call refresh_periodic_shell_mean_state(plan%periodic_shell_mean, src_q, state%periodic_shell_mean_state)
    end if
    if (plan%periodic_root_operator_ready) then
      if (allocated(state%periodic_multipole)) then
        if (any(shape(state%periodic_multipole) /= [plan%ncoef, plan%periodic_root_source_count])) then
          deallocate (state%periodic_multipole)
        end if
      end if
      if (.not. allocated(state%periodic_multipole)) then
        allocate (state%periodic_multipole(plan%ncoef, plan%periodic_root_source_count))
      end if
      state%periodic_multipole = 0.0_dp
      do source_idx = 1_i32, plan%nsrc
        cell = plan%periodic_source_cell(source_idx)
        state%periodic_multipole(:, cell) = state%periodic_multipole(:, cell) + &
                                            src_q(source_idx)*plan%periodic_source_p2m_basis(:, source_idx)
      end do
      if (allocated(state%anchor_local)) then
        if (any(shape(state%anchor_local) /= [plan%ncoef, plan%periodic_operator_anchor_count])) then
          deallocate (state%anchor_local, state%anchor_local_active)
        end if
      end if
      if (.not. allocated(state%anchor_local)) then
        allocate (state%anchor_local(plan%ncoef, plan%periodic_operator_anchor_count), &
                  state%anchor_local_active(plan%periodic_operator_anchor_count))
      end if
      state%anchor_local = 0.0_dp
      state%anchor_local_active = 0_i32
    end if
    n_target_nodes = merge(plan%target_nnode, plan%nnode, plan%target_tree_ready)
    rank = 0_i32
    rank_count = 1_i32
    if (present(mpi)) call mpi_get_rank_size(rank, rank_count, mpi)
    call m2l_rank_blocks(plan, n_target_nodes, rank_count, block_counts, block_offsets)
    block_begin = block_offsets(rank + 1_i32) + 1_i32
    block_end = block_offsets(rank + 1_i32) + block_counts(rank + 1_i32)
    !$omp parallel default(none) shared(plan, state, src_q, n_target_nodes, block_begin, block_end)
    if (associated(state%multipole_active)) then
      !$omp workshare
      state%multipole_active = 0_i32
      !$omp end workshare
    end if
    if (associated(state%local_active)) then
      !$omp workshare
      state%local_active = 0_i32
      !$omp end workshare
    end if
    if (plan%nsrc > 0_i32) then
      !$omp workshare
      state%src_q = src_q
      !$omp end workshare
    end if
    if (plan%nsource_leaf <= 0_i32 .and. associated(state%multipole)) then
      !$omp workshare
      state%multipole = 0.0d0
      !$omp end workshare
    end if
    if (plan%m2l_pair_count <= 0_i32 .and. associated(state%local)) then
      !$omp workshare
      state%local = 0.0d0
      !$omp end workshare
    end if

    call p2m_leaf_moments(plan, state)
    call m2m_upward_pass(plan, state)
    if (n_target_nodes > 0_i32) call m2l_accumulate(plan, state, block_begin, block_end)
    !$omp end parallel

    if (rank_count > 1_i32 .and. n_target_nodes > 0_i32 .and. plan%m2l_pair_count > 0_i32) then
      call mpi_allgatherv_inplace_real_dp(mpi, state%local, block_counts*plan%ncoef, block_offsets*plan%ncoef)
      call mpi_allgatherv_inplace_i32(mpi, state%local_active, block_counts, block_offsets)
    end if

    !$omp parallel default(none) shared(plan, state, n_target_nodes)
    if (n_target_nodes > 0_i32) then
      call apply_periodic_target_correction(plan, state)
      call l2l_downward_pass(plan, state)
    end if
    !$omp end parallel

    state%ready = .true.
    state%update_count = state%update_count + 1_i32
  end subroutine core_update_state_impl

  !> FMM state に確保した資源を解放する。
  !! @param[inout] state 解放対象の state。
  subroutine core_destroy_state_impl(state)
    type(fmm_state_type), intent(inout) :: state

    call reset_fmm_state(state)
  end subroutine core_destroy_state_impl

  subroutine ensure_state_capacity(plan, state)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    integer(i32) :: n_target_nodes

    n_target_nodes = plan%nnode
    if (plan%target_tree_ready) n_target_nodes = plan%target_nnode

    if (state%update_count <= 0_i32) then
      call initialize_fmm_state(state)
      if (associated(state%src_q)) deallocate (state%src_q)
      if (associated(state%multipole)) deallocate (state%multipole)
      if (associated(state%multipole_active)) deallocate (state%multipole_active)
      if (associated(state%local)) deallocate (state%local)
      if (associated(state%local_active)) deallocate (state%local_active)
      allocate (state%src_q(plan%nsrc))
      allocate (state%multipole(plan%ncoef, max(1_i32, plan%nnode)))
      allocate (state%multipole_active(max(1_i32, plan%nnode)))
      allocate (state%local(plan%ncoef, max(1_i32, n_target_nodes)))
      allocate (state%local_active(max(1_i32, n_target_nodes)))
      return
    end if

    if (associated(state%src_q)) then
      if (size(state%src_q) /= plan%nsrc) deallocate (state%src_q)
    end if
    if (.not. associated(state%src_q)) allocate (state%src_q(plan%nsrc))

    if (associated(state%multipole)) then
      if (size(state%multipole, 1) /= plan%ncoef .or. size(state%multipole, 2) /= max(1_i32, plan%nnode)) then
        deallocate (state%multipole)
      end if
    end if
    if (.not. associated(state%multipole)) allocate (state%multipole(plan%ncoef, max(1_i32, plan%nnode)))
    if (associated(state%multipole_active)) then
      if (size(state%multipole_active) /= max(1_i32, plan%nnode)) deallocate (state%multipole_active)
    end if
    if (.not. associated(state%multipole_active)) allocate (state%multipole_active(max(1_i32, plan%nnode)))

    if (associated(state%local)) then
      if (size(state%local, 1) /= plan%ncoef .or. size(state%local, 2) /= max(1_i32, n_target_nodes)) then
        deallocate (state%local)
      end if
    end if
    if (.not. associated(state%local)) allocate (state%local(plan%ncoef, max(1_i32, n_target_nodes)))
    if (associated(state%local_active)) then
      if (size(state%local_active) /= max(1_i32, n_target_nodes)) deallocate (state%local_active)
    end if
    if (.not. associated(state%local_active)) allocate (state%local_active(max(1_i32, n_target_nodes)))
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
          leaf_multipole(alpha_idx) = leaf_multipole(alpha_idx) + state%src_q(idx)*plan%source_p2m_basis(alpha_idx, p)
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
                                         *plan%source_shift_monomial(delta_idx, child_node)
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

  !> M2L の target node を rank 数で分け、各 rank の node 数と先頭 offset を返す。組の数がほぼ等しくなる
  !! 連続ブロックにする。組がなければ全 node を rank 0 に置く。
  subroutine m2l_rank_blocks(plan, n_target_nodes, rank_count, block_counts, block_offsets)
    type(fmm_plan_type), intent(in) :: plan
    integer(i32), intent(in) :: n_target_nodes, rank_count
    integer(i32), allocatable, intent(out) :: block_counts(:), block_offsets(:)
    integer(i32) :: t_node, owner, pairs_before
    integer(i64) :: pair_total

    allocate (block_counts(rank_count), block_offsets(rank_count))
    block_counts = 0_i32
    block_offsets = 0_i32
    pair_total = int(plan%m2l_pair_count, i64)
    do t_node = 1_i32, n_target_nodes
      owner = 0_i32
      if (rank_count > 1_i32 .and. pair_total > 0_i64 .and. allocated(plan%m2l_target_start)) then
        pairs_before = plan%m2l_target_start(t_node) - 1_i32
        owner = int(min(int(rank_count - 1_i32, i64), (int(pairs_before, i64)*int(rank_count, i64))/pair_total), i32)
      end if
      block_counts(owner + 1_i32) = block_counts(owner + 1_i32) + 1_i32
    end do
    do owner = 2_i32, rank_count
      block_offsets(owner) = block_offsets(owner - 1_i32) + block_counts(owner - 1_i32)
    end do
  end subroutine m2l_rank_blocks

  !> target node `t_begin` から `t_end` の局所展開を M2L で求める。
  subroutine m2l_accumulate(plan, state, t_begin, t_end)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    integer(i32), intent(in) :: t_begin, t_end
    integer(i32) :: t_node, pair_pos, pair_idx, s_node, alpha_idx, beta_idx, deriv_idx
    real(dp) :: local_acc(plan%ncoef), source_coeff
    logical :: local_active

    if (plan%m2l_pair_count <= 0_i32) return
    !$omp do schedule(static)
    do t_node = t_begin, t_end
      local_acc = 0.0d0
      local_active = .false.
      do pair_pos = plan%m2l_target_start(t_node), plan%m2l_target_start(t_node + 1_i32) - 1_i32
        pair_idx = plan%m2l_pair_order(pair_pos)
        s_node = plan%m2l_source_nodes(pair_idx)
        if (state%multipole_active(s_node) == 0_i32) cycle
        do beta_idx = 1_i32, plan%ncoef
          source_coeff = plan%alpha_sign(beta_idx)*state%multipole(beta_idx, s_node)
          if (abs(source_coeff) <= tiny(1.0d0)) cycle
          do alpha_idx = 1_i32, plan%ncoef
            deriv_idx = plan%alpha_beta_deriv_idx(alpha_idx, beta_idx)
            local_acc(alpha_idx) = local_acc(alpha_idx) + source_coeff*plan%m2l_deriv(deriv_idx, pair_idx)
          end do
        end do
      end do
      local_active = coefficient_vector_is_active(local_acc)
      state%local(:, t_node) = local_acc
      state%local_active(t_node) = merge(1_i32, 0_i32, local_active)
    end do
    !$omp end do
  end subroutine m2l_accumulate

  !> 遠方補正を anchor 格子の局所展開として集計する。点評価で、点を含む anchor の展開を葉の局所展開に加える。
  subroutine apply_periodic_target_correction(plan, state)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    integer(i32) :: target_idx, source_idx, slot, anchor_count
    real(dp) :: node_local(plan%ncoef)
    logical :: node_active

    if (.not. plan%periodic_root_operator_ready) return
    if (.not. allocated(state%anchor_local)) return
    anchor_count = plan%periodic_operator_anchor_count
    if (anchor_count <= 0_i32 .or. plan%periodic_root_source_count <= 0_i32) return
    ! 組は source cell ごとに target anchor が連続して並ぶ。anchor 単位で並列化して書き込みの競合を避ける。
    !$omp do schedule(static)
    do target_idx = 1_i32, anchor_count
      node_local = 0.0_dp
      do source_idx = 1_i32, plan%periodic_root_source_count
        if (.not. coefficient_vector_is_active(state%periodic_multipole(:, source_idx))) cycle
        slot = plan%periodic_operator_cell_slot(source_idx)
        ! Every cell that holds a source was built at plan time, so a charged cell always has a slot.
        if (slot <= 0_i32) error stop 'periodic far correction has no operator for a charged source cell.'
        node_local = node_local + matmul(plan%periodic_root_operator(:, :, (slot - 1_i32)*anchor_count + target_idx), &
                                         state%periodic_multipole(:, source_idx))
      end do
      node_active = coefficient_vector_is_active(node_local)
      state%anchor_local(:, target_idx) = node_local
      state%anchor_local_active(target_idx) = merge(1_i32, 0_i32, node_active)
    end do
    !$omp end do
  end subroutine apply_periodic_target_correction

  subroutine l2l_downward_pass(plan, state)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(inout) :: state
    integer(i32) :: depth, max_depth, level_pos, level_start_pos, level_end_pos
    integer(i32) :: node_idx, parent_node, alpha_idx, term_idx, gamma_idx, delta_idx
    real(dp) :: parent_local(plan%ncoef), local_acc(plan%ncoef)
    logical :: use_target_tree
    logical :: parent_active, node_active

    use_target_tree = plan%target_tree_ready
    max_depth = active_tree_max_depth(plan, use_target_tree)
    if (max_depth <= 0_i32) return

    do depth = 1_i32, max_depth
      call active_tree_level_bounds(plan, use_target_tree, depth, level_start_pos, level_end_pos)
      !$omp do schedule(static)
      do level_pos = level_start_pos, level_end_pos
        node_idx = active_tree_level_node(plan, use_target_tree, level_pos)
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
                                     *plan%target_shift_monomial(delta_idx, node_idx)
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

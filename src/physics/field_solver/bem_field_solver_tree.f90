!> `bem_field_solver` の octree 構築・更新とメモリ管理を実装する submodule。
submodule (bem_field_solver) bem_field_solver_tree
  use bem_coulomb_fmm_core, only: build_plan, update_state, destroy_plan, destroy_state, fmm_options_type
  implicit none
contains

  !> 現在の要素電荷から各ノードの monopole モーメントを再集計する。
  module procedure refresh_field_solver
    integer(i32) :: node_idx, child_k, child_node
    integer(i32) :: p, idx, p_end
    integer(i32) :: depth, level_pos, level_start_pos, level_end_pos
    real(dp) :: q, abs_q, qx, qy, qz, qi
    real(dp) :: t_moment_start, t_moment_end, avg_refresh
    real(dp), allocatable :: src_pos(:, :)

    if (trim(self%mode) /= 'treecode' .and. trim(self%mode) /= 'fmm') return
    if (mesh%nelem <= 0_i32) return

    if (trim(self%mode) == 'fmm' .and. self%fmm_use_core) then
      t_moment_start = field_solver_time_seconds()
      if (.not. self%fmm_core_plan%built .or. self%fmm_core_plan%nsrc /= mesh%nelem) then
        call destroy_plan(self%fmm_core_plan)
        call destroy_state(self%fmm_core_state)
        allocate (src_pos(3, mesh%nelem))
        do idx = 1_i32, mesh%nelem
          src_pos(1, idx) = mesh%center_x(idx)
          src_pos(2, idx) = mesh%center_y(idx)
          src_pos(3, idx) = mesh%center_z(idx)
        end do
        call build_plan(self%fmm_core_plan, src_pos, self%fmm_core_options)
        deallocate (src_pos)
      end if

      call update_state(self%fmm_core_plan, self%fmm_core_state, mesh%q_elem)
      t_moment_end = field_solver_time_seconds()
      self%fmm_core_ready = self%fmm_core_plan%built .and. self%fmm_core_state%ready
      self%fmm_last_moment_time_s = 0.0d0
      self%fmm_last_clear_time_s = 0.0d0
      self%fmm_last_m2l_time_s = 0.0d0
      self%fmm_last_l2l_time_s = 0.0d0
      self%fmm_last_copy_time_s = 0.0d0
      self%fmm_last_refresh_time_s = t_moment_end - t_moment_start
      self%fmm_total_refresh_time_s = self%fmm_total_refresh_time_s + self%fmm_last_refresh_time_s
      self%fmm_m2l_pair_count = self%fmm_core_plan%m2l_pair_count
      self%fmm_m2l_build_count = self%fmm_core_plan%m2l_build_count
      self%fmm_m2l_visit_count = self%fmm_core_plan%m2l_visit_count
      if (allocated(self%fmm_core_plan%near_nodes)) then
        self%fmm_near_interaction_count = int(size(self%fmm_core_plan%near_nodes), i32)
      else
        self%fmm_near_interaction_count = 0_i32
      end if
      if (allocated(self%fmm_core_plan%far_nodes)) then
        self%fmm_far_interaction_count = int(size(self%fmm_core_plan%far_nodes), i32)
      else
        self%fmm_far_interaction_count = 0_i32
      end if
      self%fmm_refresh_count = self%fmm_refresh_count + 1_i32
      return
    end if

    if (.not. self%tree_ready .or. self%nelem /= mesh%nelem) then
      self%nelem = mesh%nelem
      call build_tree_topology(self, mesh)
    end if

    if (.not. allocated(self%node_level_start) .or. .not. allocated(self%node_level_nodes)) then
      call rebuild_source_level_cache(self)
    end if

    if (trim(self%mode) == 'fmm') t_moment_start = field_solver_time_seconds()
    do depth = self%node_max_depth, 0_i32, -1_i32
      level_start_pos = self%node_level_start(depth + 1_i32)
      level_end_pos = self%node_level_start(depth + 2_i32) - 1_i32
      !$omp parallel do default(none) schedule(static) &
      !$omp shared(self,mesh,level_start_pos,level_end_pos) &
      !$omp private(level_pos,node_idx,child_k,child_node,p,idx,p_end,q,abs_q,qx,qy,qz,qi)
      do level_pos = level_start_pos, level_end_pos
        node_idx = self%node_level_nodes(level_pos)
        if (self%child_count(node_idx) <= 0_i32) then
          q = 0.0d0
          abs_q = 0.0d0
          qx = 0.0d0
          qy = 0.0d0
          qz = 0.0d0
          p_end = self%node_start(node_idx) + self%node_count(node_idx) - 1_i32
          do p = self%node_start(node_idx), p_end
            idx = self%elem_order(p)
            qi = mesh%q_elem(idx)
            q = q + qi
            abs_q = abs_q + abs(qi)
            qx = qx + qi * mesh%center_x(idx)
            qy = qy + qi * mesh%center_y(idx)
            qz = qz + qi * mesh%center_z(idx)
          end do
        else
          q = 0.0d0
          abs_q = 0.0d0
          qx = 0.0d0
          qy = 0.0d0
          qz = 0.0d0
          do child_k = 1_i32, self%child_count(node_idx)
            child_node = self%child_idx(child_k, node_idx)
            q = q + self%node_q(child_node)
            abs_q = abs_q + self%node_abs_q(child_node)
            qx = qx + self%node_qx(child_node)
            qy = qy + self%node_qy(child_node)
            qz = qz + self%node_qz(child_node)
          end do
        end if

        self%node_q(node_idx) = q
        self%node_abs_q(node_idx) = abs_q
        self%node_qx(node_idx) = qx
        self%node_qy(node_idx) = qy
        self%node_qz(node_idx) = qz

        if (abs(q) > tiny(1.0d0)) then
          self%node_charge_center(1, node_idx) = qx / q
          self%node_charge_center(2, node_idx) = qy / q
          self%node_charge_center(3, node_idx) = qz / q
        else
          self%node_charge_center(:, node_idx) = self%node_center(:, node_idx)
        end if
      end do
      !$omp end parallel do
    end do
    if (trim(self%mode) == 'fmm') then
      t_moment_end = field_solver_time_seconds()
      self%fmm_last_moment_time_s = t_moment_end - t_moment_start
    end if

    if (trim(self%mode) == 'fmm') then
      if (.not. self%fmm_ready) call build_fmm_interactions(self)
      call refresh_fmm_locals(self, mesh)
      self%fmm_refresh_count = self%fmm_refresh_count + 1_i32
      self%fmm_last_refresh_time_s = self%fmm_last_moment_time_s + self%fmm_last_clear_time_s &
                                     + self%fmm_last_m2l_time_s + self%fmm_last_l2l_time_s &
                                     + self%fmm_last_copy_time_s
      self%fmm_total_moment_time_s = self%fmm_total_moment_time_s + self%fmm_last_moment_time_s
      self%fmm_total_clear_time_s = self%fmm_total_clear_time_s + self%fmm_last_clear_time_s
      self%fmm_total_m2l_time_s = self%fmm_total_m2l_time_s + self%fmm_last_m2l_time_s
      self%fmm_total_l2l_time_s = self%fmm_total_l2l_time_s + self%fmm_last_l2l_time_s
      self%fmm_total_copy_time_s = self%fmm_total_copy_time_s + self%fmm_last_copy_time_s
      self%fmm_total_refresh_time_s = self%fmm_total_refresh_time_s + self%fmm_last_refresh_time_s
      if (self%fmm_profile_enabled) then
        if (self%fmm_refresh_count <= 3_i32 .or. mod(self%fmm_refresh_count, 50_i32) == 0_i32) then
          avg_refresh = self%fmm_total_refresh_time_s / real(self%fmm_refresh_count, dp)
          write (*, '(A,I0,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
            'FMM profile refresh#', self%fmm_refresh_count, &
            ' moments=', self%fmm_last_moment_time_s, &
            ' clear=', self%fmm_last_clear_time_s, &
            ' m2l=', self%fmm_last_m2l_time_s, &
            ' l2l=', self%fmm_last_l2l_time_s, &
            ' copy=', self%fmm_last_copy_time_s, &
            ' total=', self%fmm_last_refresh_time_s, &
            ' avg=', avg_refresh
        end if
      end if
    else
      self%fmm_ready = .false.
    end if
  end procedure refresh_field_solver

  !> 要素重心を octree に再配置して木構造トポロジを構築する。
  module procedure build_tree_topology
    integer(i32) :: i, max_node_guess

    if (mesh%nelem <= 0_i32) then
      call reset_tree_storage(self)
      self%tree_ready = .false.
      self%nelem = 0_i32
      return
    end if

    max_node_guess = max(1_i32, 2_i32 * mesh%nelem)
    call ensure_tree_capacity(self, max_node_guess)

    if (allocated(self%elem_order)) then
      if (size(self%elem_order) /= mesh%nelem) deallocate (self%elem_order)
    end if
    if (.not. allocated(self%elem_order)) allocate (self%elem_order(mesh%nelem))

    !$omp parallel do default(none) schedule(static) &
    !$omp shared(self,mesh) private(i)
    do i = 1_i32, mesh%nelem
      self%elem_order(i) = i
    end do
    !$omp end parallel do

    self%nnode = 1_i32
    self%node_max_depth = 0_i32
    call build_node(self, mesh, 1_i32, 1_i32, mesh%nelem, 0_i32)
    call rebuild_source_level_cache(self)
    self%nelem = mesh%nelem
    self%tree_ready = .true.
    self%fmm_ready = .false.
  end procedure build_tree_topology

  !> 指定区間の要素を1ノードとして登録し、条件を満たせば8分木分割する。
  module procedure build_node
    integer(i32) :: count, p, idx, oct
    integer(i32) :: child_k, child_node, child_start, child_end
    integer(i32), allocatable :: counts(:), offsets(:), cursor(:), work(:)
    real(dp) :: bb_min(3), bb_max(3), span(3), center(3)
    real(dp) :: split_eps

    count = end_idx - start_idx + 1_i32
    self%node_start(node_idx) = start_idx
    self%node_count(node_idx) = count
    self%child_count(node_idx) = 0_i32
    self%child_idx(:, node_idx) = 0_i32
    self%child_octant(:, node_idx) = 0_i32
    self%node_depth(node_idx) = depth
    self%node_max_depth = max(self%node_max_depth, depth)

    idx = self%elem_order(start_idx)
    bb_min = [mesh%center_x(idx), mesh%center_y(idx), mesh%center_z(idx)]
    bb_max = bb_min
    do p = start_idx + 1_i32, end_idx
      idx = self%elem_order(p)
      bb_min(1) = min(bb_min(1), mesh%center_x(idx))
      bb_min(2) = min(bb_min(2), mesh%center_y(idx))
      bb_min(3) = min(bb_min(3), mesh%center_z(idx))
      bb_max(1) = max(bb_max(1), mesh%center_x(idx))
      bb_max(2) = max(bb_max(2), mesh%center_y(idx))
      bb_max(3) = max(bb_max(3), mesh%center_z(idx))
    end do

    span = bb_max - bb_min
    center = 0.5d0 * (bb_max + bb_min)
    self%node_center(:, node_idx) = center
    self%node_half_size(:, node_idx) = 0.5d0 * span
    self%node_radius(node_idx) = sqrt(sum(self%node_half_size(:, node_idx) * self%node_half_size(:, node_idx)))

    if (count <= self%leaf_max) return

    split_eps = 1.0d-12 * max(1.0d0, maxval(abs(center)))
    if (maxval(span) <= split_eps) return

    allocate (counts(8), offsets(8), cursor(8), work(count))
    counts = 0_i32
    do p = start_idx, end_idx
      idx = self%elem_order(p)
      oct = octant_index(mesh%center_x(idx), mesh%center_y(idx), mesh%center_z(idx), center)
      counts(oct) = counts(oct) + 1_i32
    end do

    if (maxval(counts) == count) then
      deallocate (counts, offsets, cursor, work)
      return
    end if

    offsets(1) = 1_i32
    do oct = 2, 8
      offsets(oct) = offsets(oct - 1) + counts(oct - 1)
    end do
    cursor = offsets

    do p = start_idx, end_idx
      idx = self%elem_order(p)
      oct = octant_index(mesh%center_x(idx), mesh%center_y(idx), mesh%center_z(idx), center)
      work(cursor(oct)) = idx
      cursor(oct) = cursor(oct) + 1_i32
    end do
    self%elem_order(start_idx:end_idx) = work

    child_k = 0_i32
    child_start = start_idx
    do oct = 1, 8
      if (counts(oct) <= 0_i32) cycle
      child_end = child_start + counts(oct) - 1_i32
      self%nnode = self%nnode + 1_i32
      if (self%nnode > self%max_node) error stop 'field solver tree capacity exceeded.'
      child_node = self%nnode
      child_k = child_k + 1_i32
      self%child_idx(child_k, node_idx) = child_node
      self%child_octant(child_k, node_idx) = oct
      call build_node(self, mesh, child_node, child_start, child_end, depth + 1_i32)
      child_start = child_end + 1_i32
    end do
    self%child_count(node_idx) = child_k

    deallocate (counts, offsets, cursor, work)
  end procedure build_node

  !> ノード中心に対する座標の相対位置から octant インデックスを返す。
  module procedure octant_index
    oct = 1_i32
    if (x >= center(1)) oct = oct + 1_i32
    if (y >= center(2)) oct = oct + 2_i32
    if (z >= center(3)) oct = oct + 4_i32
  end procedure octant_index

  !> 必要ノード数に合わせて木配列を確保し、利用領域をゼロ初期化する。
  module procedure ensure_tree_capacity
    if (self%max_node == max_node_needed .and. allocated(self%node_start)) then
      self%node_start = 0_i32
      self%node_count = 0_i32
      self%child_count = 0_i32
      self%child_idx = 0_i32
      self%child_octant = 0_i32
      self%node_depth = 0_i32
      self%node_max_depth = 0_i32
      self%node_center = 0.0d0
      self%node_half_size = 0.0d0
      self%node_radius = 0.0d0
      self%node_q = 0.0d0
      self%node_abs_q = 0.0d0
      self%node_qx = 0.0d0
      self%node_qy = 0.0d0
      self%node_qz = 0.0d0
      self%node_charge_center = 0.0d0
      return
    end if

    call reset_tree_storage(self)
    self%max_node = max_node_needed
    allocate (self%node_start(self%max_node), self%node_count(self%max_node))
    allocate (self%child_count(self%max_node), self%child_idx(8, self%max_node), self%child_octant(8, self%max_node))
    allocate (self%node_depth(self%max_node))
    allocate (self%node_center(3, self%max_node), self%node_half_size(3, self%max_node))
    allocate (self%node_radius(self%max_node))
    allocate (self%node_q(self%max_node), self%node_abs_q(self%max_node))
    allocate (self%node_qx(self%max_node), self%node_qy(self%max_node), self%node_qz(self%max_node))
    allocate (self%node_charge_center(3, self%max_node))

    self%node_start = 0_i32
    self%node_count = 0_i32
    self%child_count = 0_i32
    self%child_idx = 0_i32
    self%child_octant = 0_i32
    self%node_depth = 0_i32
    self%node_max_depth = 0_i32
    self%node_center = 0.0d0
    self%node_half_size = 0.0d0
    self%node_radius = 0.0d0
    self%node_q = 0.0d0
    self%node_abs_q = 0.0d0
    self%node_qx = 0.0d0
    self%node_qy = 0.0d0
    self%node_qz = 0.0d0
    self%node_charge_center = 0.0d0
  end procedure ensure_tree_capacity

  module procedure rebuild_source_level_cache
    integer(i32) :: node_idx, depth
    integer(i32), allocatable :: depth_count(:), cursor(:)

    if (allocated(self%node_level_start)) deallocate (self%node_level_start)
    if (allocated(self%node_level_nodes)) deallocate (self%node_level_nodes)
    if (self%nnode <= 0_i32) return

    allocate (depth_count(self%node_max_depth + 1_i32), cursor(self%node_max_depth + 1_i32))
    depth_count = 0_i32
    do node_idx = 1_i32, self%nnode
      depth_count(self%node_depth(node_idx) + 1_i32) = depth_count(self%node_depth(node_idx) + 1_i32) + 1_i32
    end do

    allocate (self%node_level_start(self%node_max_depth + 2_i32), self%node_level_nodes(self%nnode))
    self%node_level_start(1) = 1_i32
    do depth = 0_i32, self%node_max_depth
      self%node_level_start(depth + 2_i32) = self%node_level_start(depth + 1_i32) + depth_count(depth + 1_i32)
    end do

    cursor = self%node_level_start(1:self%node_max_depth + 1_i32)
    do node_idx = 1_i32, self%nnode
      depth = self%node_depth(node_idx)
      self%node_level_nodes(cursor(depth + 1_i32)) = node_idx
      cursor(depth + 1_i32) = cursor(depth + 1_i32) + 1_i32
    end do

    deallocate (depth_count, cursor)
  end procedure rebuild_source_level_cache

  !> treecode で使う作業配列を解放し、ノード状態を未初期化に戻す。
  module procedure reset_tree_storage
    call destroy_plan(self%fmm_core_plan)
    call destroy_state(self%fmm_core_state)
    if (allocated(self%elem_order)) deallocate (self%elem_order)
    if (allocated(self%node_start)) deallocate (self%node_start)
    if (allocated(self%node_count)) deallocate (self%node_count)
    if (allocated(self%child_count)) deallocate (self%child_count)
    if (allocated(self%child_idx)) deallocate (self%child_idx)
    if (allocated(self%child_octant)) deallocate (self%child_octant)
    if (allocated(self%node_depth)) deallocate (self%node_depth)
    if (allocated(self%node_level_start)) deallocate (self%node_level_start)
    if (allocated(self%node_level_nodes)) deallocate (self%node_level_nodes)
    if (allocated(self%node_center)) deallocate (self%node_center)
    if (allocated(self%node_half_size)) deallocate (self%node_half_size)
    if (allocated(self%node_radius)) deallocate (self%node_radius)
    if (allocated(self%node_q)) deallocate (self%node_q)
    if (allocated(self%node_abs_q)) deallocate (self%node_abs_q)
    if (allocated(self%node_qx)) deallocate (self%node_qx)
    if (allocated(self%node_qy)) deallocate (self%node_qy)
    if (allocated(self%node_qz)) deallocate (self%node_qz)
    if (allocated(self%node_charge_center)) deallocate (self%node_charge_center)
    if (allocated(self%leaf_nodes)) deallocate (self%leaf_nodes)
    if (allocated(self%leaf_slot_of_node)) deallocate (self%leaf_slot_of_node)
    if (allocated(self%target_child_count)) deallocate (self%target_child_count)
    if (allocated(self%target_child_idx)) deallocate (self%target_child_idx)
    if (allocated(self%target_child_octant)) deallocate (self%target_child_octant)
    if (allocated(self%target_node_depth)) deallocate (self%target_node_depth)
    if (allocated(self%target_level_start)) deallocate (self%target_level_start)
    if (allocated(self%target_level_nodes)) deallocate (self%target_level_nodes)
    if (allocated(self%target_node_center)) deallocate (self%target_node_center)
    if (allocated(self%target_node_half_size)) deallocate (self%target_node_half_size)
    if (allocated(self%target_node_radius)) deallocate (self%target_node_radius)
    if (allocated(self%near_start)) deallocate (self%near_start)
    if (allocated(self%near_nodes)) deallocate (self%near_nodes)
    if (allocated(self%far_start)) deallocate (self%far_start)
    if (allocated(self%far_nodes)) deallocate (self%far_nodes)
    if (allocated(self%fmm_m2l_target_nodes)) deallocate (self%fmm_m2l_target_nodes)
    if (allocated(self%fmm_m2l_source_nodes)) deallocate (self%fmm_m2l_source_nodes)
    if (allocated(self%fmm_m2l_target_start)) deallocate (self%fmm_m2l_target_start)
    if (allocated(self%fmm_m2l_pair_order)) deallocate (self%fmm_m2l_pair_order)
    if (allocated(self%fmm_parent_of)) deallocate (self%fmm_parent_of)
    if (allocated(self%fmm_node_local_e0)) deallocate (self%fmm_node_local_e0)
    if (allocated(self%fmm_node_local_jac)) deallocate (self%fmm_node_local_jac)
    if (allocated(self%fmm_node_local_hess)) deallocate (self%fmm_node_local_hess)
    if (allocated(self%fmm_shift_axis1)) deallocate (self%fmm_shift_axis1)
    if (allocated(self%fmm_shift_axis2)) deallocate (self%fmm_shift_axis2)
    if (allocated(self%leaf_far_e0)) deallocate (self%leaf_far_e0)
    if (allocated(self%leaf_far_jac)) deallocate (self%leaf_far_jac)
    if (allocated(self%leaf_far_hess)) deallocate (self%leaf_far_hess)

    self%nnode = 0_i32
    self%max_node = 0_i32
    self%node_max_depth = 0_i32
    self%nleaf = 0_i32
    self%target_nnode = 0_i32
    self%target_max_node = 0_i32
    self%target_node_max_depth = 0_i32
    self%target_tree_ready = .false.
    self%tree_ready = .false.
    self%fmm_ready = .false.
    self%nelem = 0_i32
    self%fmm_m2l_pair_count = 0_i32
    self%fmm_m2l_build_count = 0_i32
    self%fmm_m2l_visit_count = 0_i32
    self%fmm_near_interaction_count = 0_i32
    self%fmm_far_interaction_count = 0_i32
    self%fmm_refresh_count = 0_i32
    self%fmm_last_moment_time_s = 0.0d0
    self%fmm_last_clear_time_s = 0.0d0
    self%fmm_last_m2l_time_s = 0.0d0
    self%fmm_last_l2l_time_s = 0.0d0
    self%fmm_last_copy_time_s = 0.0d0
    self%fmm_last_refresh_time_s = 0.0d0
    self%fmm_total_moment_time_s = 0.0d0
    self%fmm_total_clear_time_s = 0.0d0
    self%fmm_total_m2l_time_s = 0.0d0
    self%fmm_total_l2l_time_s = 0.0d0
    self%fmm_total_copy_time_s = 0.0d0
    self%fmm_total_refresh_time_s = 0.0d0
    self%fmm_use_core = .false.
    self%fmm_core_ready = .false.
    self%fmm_core_options = fmm_options_type()
  end procedure reset_tree_storage

  module procedure field_solver_time_seconds
    call cpu_time(time_s)
    !$ time_s = omp_get_wtime()
  end procedure field_solver_time_seconds

end submodule bem_field_solver_tree

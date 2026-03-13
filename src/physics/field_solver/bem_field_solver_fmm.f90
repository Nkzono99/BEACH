!> `bem_field_solver` の FMM 相互作用構築と評価を実装する submodule。
submodule (bem_field_solver) bem_field_solver_fmm
  implicit none
  real(dp), parameter :: inv_sqrt_pi = 0.56418958354775628695d0
contains

  !> source/target ノード球が十分離れていれば遠方相互作用として扱う。
  module procedure nodes_well_separated
    real(dp) :: d(3), dist2, rs, rt, theta_eff, lhs, rhs
    real(dp) :: target_center(3)

    if (self%target_tree_ready) then
      target_center = self%target_node_center(:, target_node)
      rt = self%target_node_radius(target_node)
    else
      target_center = self%node_center(:, target_node)
      rt = self%node_radius(target_node)
    end if

    d = target_center - self%node_center(:, source_node)
    call apply_periodic2_minimum_image(self, d)
    dist2 = sum(d * d)
    if (dist2 <= 0.0d0) then
      accept_it = .false.
      return
    end if

    rs = self%node_radius(source_node)
    theta_eff = self%theta
    if (self%use_periodic2) theta_eff = 0.5d0 * self%theta
    lhs = (rs + rt) * (rs + rt)
    rhs = (theta_eff * theta_eff) * dist2
    accept_it = (lhs < rhs)
  end procedure nodes_well_separated

  !> 1つの葉ノードに対する near/far リストを木走査で構築する。
  module procedure gather_leaf_interactions
    integer(i32) :: child_k, child_node

    if (nodes_well_separated(self, target_leaf, source_node)) then
      far_n = far_n + 1_i32
      far_buf(far_n) = source_node
      return
    end if

    if (self%child_count(source_node) <= 0_i32) then
      near_n = near_n + 1_i32
      near_buf(near_n) = source_node
      return
    end if

    do child_k = 1_i32, self%child_count(source_node)
      child_node = self%child_idx(child_k, source_node)
      call gather_leaf_interactions(self, target_leaf, child_node, near_buf, near_n, far_buf, far_n)
    end do
  end procedure gather_leaf_interactions

  !> 全葉ノードの near/far 相互作用リストを構築し、評価用配列を確保する。
  module procedure build_fmm_interactions
    integer(i32) :: leaf_idx, i
    integer(i32) :: near_n, far_n, near_used, far_used
    integer(i32) :: near_cap, far_cap
    integer(i32) :: n_target_nodes, child_k, parent_node, node_idx
    integer(i32) :: pair_cap, pair_used
    integer(i32) :: nshift, images_per_pair
    integer(i32), allocatable :: near_work(:), far_work(:), near_buf(:), far_buf(:)
    integer(i32), allocatable :: m2l_target_buf(:), m2l_source_buf(:)
    logical :: use_box_target, use_target_tree
    real(dp) :: build_t0, build_t1, near_far_t0, near_far_t1, pair_t0, pair_t1

    if (.not. self%tree_ready .or. self%nnode <= 0_i32) then
      self%fmm_ready = .false.
      self%nleaf = 0_i32
      self%target_tree_ready = .false.
      return
    end if

    use_box_target = has_valid_target_box(self)
    if (use_box_target) then
      call build_box_target_tree_helper(self, self%use_periodic2)
    else
      call assign_source_leaves_as_targets_helper(self)
    end if

    if (self%nleaf <= 0_i32) then
      self%fmm_ready = .false.
      return
    end if
    use_target_tree = self%target_tree_ready
    n_target_nodes = active_tree_nnode(self, use_target_tree)

    if (allocated(self%near_start)) deallocate (self%near_start)
    if (allocated(self%far_start)) deallocate (self%far_start)
    if (allocated(self%near_nodes)) deallocate (self%near_nodes)
    if (allocated(self%far_nodes)) deallocate (self%far_nodes)
    allocate (self%near_start(self%nleaf + 1_i32), self%far_start(self%nleaf + 1_i32))

    allocate (near_work(self%nnode), far_work(self%nnode))
    near_cap = max(64_i32, self%nleaf * 16_i32)
    far_cap = max(64_i32, self%nleaf * 16_i32)
    allocate (near_buf(near_cap), far_buf(far_cap))
    near_used = 0_i32
    far_used = 0_i32
    build_t0 = field_solver_time_seconds()
    near_far_t0 = field_solver_time_seconds()

    do leaf_idx = 1_i32, self%nleaf
      near_n = 0_i32
      far_n = 0_i32
      call gather_leaf_interactions( &
        self, self%leaf_nodes(leaf_idx), 1_i32, near_work, near_n, far_work, far_n &
      )
      self%near_start(leaf_idx) = near_used + 1_i32
      self%far_start(leaf_idx) = far_used + 1_i32

      do i = 1_i32, near_n
        call append_i32_buffer(near_buf, near_used, near_cap, near_work(i))
      end do
      do i = 1_i32, far_n
        call append_i32_buffer(far_buf, far_used, far_cap, far_work(i))
      end do
    end do
    self%near_start(self%nleaf + 1_i32) = near_used + 1_i32
    self%far_start(self%nleaf + 1_i32) = far_used + 1_i32
    near_far_t1 = field_solver_time_seconds()

    allocate (self%near_nodes(max(1_i32, near_used)))
    self%near_nodes = 0_i32
    if (near_used > 0_i32) self%near_nodes(1:near_used) = near_buf(1:near_used)
    allocate (self%far_nodes(max(1_i32, far_used)))
    self%far_nodes = 0_i32
    if (far_used > 0_i32) self%far_nodes(1:far_used) = far_buf(1:far_used)
    self%fmm_near_interaction_count = near_used
    self%fmm_far_interaction_count = far_used

    if (allocated(self%fmm_parent_of)) deallocate (self%fmm_parent_of)
    allocate (self%fmm_parent_of(n_target_nodes))
    self%fmm_parent_of = 0_i32
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(self,n_target_nodes,use_target_tree) private(parent_node,child_k,node_idx)
    do parent_node = 1_i32, n_target_nodes
      do child_k = 1_i32, active_tree_child_count(self, use_target_tree, parent_node)
        node_idx = active_tree_child_idx(self, use_target_tree, child_k, parent_node)
        self%fmm_parent_of(node_idx) = parent_node
      end do
    end do
    !$omp end parallel do

    nshift = 1_i32
    if (self%use_periodic2) nshift = 2_i32 * max(0_i32, self%periodic_image_layers) + 1_i32
    if (allocated(self%fmm_shift_axis1)) deallocate (self%fmm_shift_axis1)
    if (allocated(self%fmm_shift_axis2)) deallocate (self%fmm_shift_axis2)
    allocate (self%fmm_shift_axis1(nshift), self%fmm_shift_axis2(nshift))
    call build_periodic_shift_values(self, self%fmm_shift_axis1, self%fmm_shift_axis2, nshift)

    if (allocated(self%leaf_far_e0)) deallocate (self%leaf_far_e0)
    if (allocated(self%leaf_far_jac)) deallocate (self%leaf_far_jac)
    if (allocated(self%leaf_far_hess)) deallocate (self%leaf_far_hess)
    allocate (self%leaf_far_e0(3, self%nleaf), self%leaf_far_jac(3, 3, self%nleaf), self%leaf_far_hess(3, 3, 3, self%nleaf))
    self%leaf_far_e0 = 0.0d0
    self%leaf_far_jac = 0.0d0
    self%leaf_far_hess = 0.0d0

    if (allocated(self%fmm_node_local_e0)) deallocate (self%fmm_node_local_e0)
    if (allocated(self%fmm_node_local_jac)) deallocate (self%fmm_node_local_jac)
    if (allocated(self%fmm_node_local_hess)) deallocate (self%fmm_node_local_hess)
    allocate ( &
      self%fmm_node_local_e0(3, n_target_nodes), &
      self%fmm_node_local_jac(3, 3, n_target_nodes), &
      self%fmm_node_local_hess(3, 3, 3, n_target_nodes) &
    )
    self%fmm_node_local_e0 = 0.0d0
    self%fmm_node_local_jac = 0.0d0
    self%fmm_node_local_hess = 0.0d0

    pair_t0 = field_solver_time_seconds()
    call build_m2l_pair_cache()
    call build_fmm_target_pair_index(self, n_target_nodes)
    pair_t1 = field_solver_time_seconds()
    self%fmm_ready = .true.
    build_t1 = field_solver_time_seconds()

    if (self%fmm_profile_enabled) then
      images_per_pair = size(self%fmm_shift_axis1) * size(self%fmm_shift_axis2)
      write (*, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,ES12.4,A,ES12.4,A,ES12.4)') &
        'FMM profile build#', self%fmm_m2l_build_count, &
        ' target_nodes=', n_target_nodes, &
        ' leaves=', self%nleaf, &
        ' near=', self%fmm_near_interaction_count, &
        ' far=', self%fmm_far_interaction_count, &
        ' m2l_pairs=', self%fmm_m2l_pair_count, &
        ' visits=', self%fmm_m2l_visit_count, &
        ' images_per_pair=', images_per_pair, &
        ' nearfar_s=', near_far_t1 - near_far_t0, &
        ' m2l_build_s=', pair_t1 - pair_t0, &
        ' total_s=', build_t1 - build_t0
    end if

    deallocate (near_work, far_work, near_buf, far_buf)

  contains

    subroutine build_m2l_pair_cache()
      if (allocated(self%fmm_m2l_target_nodes)) deallocate (self%fmm_m2l_target_nodes)
      if (allocated(self%fmm_m2l_source_nodes)) deallocate (self%fmm_m2l_source_nodes)

      pair_cap = max(64_i32, self%nleaf * 32_i32)
      allocate (m2l_target_buf(pair_cap), m2l_source_buf(pair_cap))
      m2l_target_buf = 0_i32
      m2l_source_buf = 0_i32
      pair_used = 0_i32
      self%fmm_m2l_visit_count = 0_i32

      call accumulate_cached_pairs(1_i32, 1_i32)

      allocate (self%fmm_m2l_target_nodes(max(1_i32, pair_used)))
      allocate (self%fmm_m2l_source_nodes(max(1_i32, pair_used)))
      self%fmm_m2l_target_nodes = 0_i32
      self%fmm_m2l_source_nodes = 0_i32
      if (pair_used > 0_i32) then
        self%fmm_m2l_target_nodes(1:pair_used) = m2l_target_buf(1:pair_used)
        self%fmm_m2l_source_nodes(1:pair_used) = m2l_source_buf(1:pair_used)
      end if
      self%fmm_m2l_pair_count = pair_used
      self%fmm_m2l_build_count = self%fmm_m2l_build_count + 1_i32

      deallocate (m2l_target_buf, m2l_source_buf)
    end subroutine build_m2l_pair_cache

    recursive subroutine accumulate_cached_pairs(t_node, s_node)
      integer(i32), intent(in) :: t_node, s_node
      integer(i32) :: tc, sc, t_child, s_child

      self%fmm_m2l_visit_count = self%fmm_m2l_visit_count + 1_i32
      tc = active_tree_child_count(self, use_target_tree, t_node)
      sc = self%child_count(s_node)

      if (nodes_well_separated(self, t_node, s_node)) then
        call append_m2l_pair(t_node, s_node)
        return
      end if

      if (tc <= 0_i32 .and. sc <= 0_i32) return

      if (sc > 0_i32 .and. ( &
        tc <= 0_i32 .or. self%node_radius(s_node) >= active_tree_node_radius(self, use_target_tree, t_node) &
      )) then
        do s_child = 1_i32, sc
          call accumulate_cached_pairs(t_node, self%child_idx(s_child, s_node))
        end do
      else
        do t_child = 1_i32, tc
          call accumulate_cached_pairs(active_tree_child_idx(self, use_target_tree, t_child, t_node), s_node)
        end do
      end if
    end subroutine accumulate_cached_pairs

    subroutine append_m2l_pair(t_node, s_node)
      integer(i32), intent(in) :: t_node, s_node
      integer(i32), allocatable :: tmp_target(:), tmp_source(:)
      integer(i32) :: new_capacity

      if (pair_used >= pair_cap) then
        new_capacity = max(pair_cap * 2_i32, pair_cap + 64_i32)
        allocate (tmp_target(new_capacity), tmp_source(new_capacity))
        tmp_target = 0_i32
        tmp_source = 0_i32
        if (pair_used > 0_i32) then
          tmp_target(1:pair_used) = m2l_target_buf(1:pair_used)
          tmp_source(1:pair_used) = m2l_source_buf(1:pair_used)
        end if
        call move_alloc(tmp_target, m2l_target_buf)
        call move_alloc(tmp_source, m2l_source_buf)
        pair_cap = new_capacity
      end if

      pair_used = pair_used + 1_i32
      m2l_target_buf(pair_used) = t_node
      m2l_source_buf(pair_used) = s_node
    end subroutine append_m2l_pair
  end procedure build_fmm_interactions

  logical function has_valid_target_box(self)
    class(field_solver_type), intent(in) :: self
    has_valid_target_box = all(self%target_box_max > self%target_box_min)
  end function has_valid_target_box

  subroutine clear_target_tree(self)
    class(field_solver_type), intent(inout) :: self

    if (allocated(self%target_child_count)) deallocate (self%target_child_count)
    if (allocated(self%target_child_idx)) deallocate (self%target_child_idx)
    if (allocated(self%target_child_octant)) deallocate (self%target_child_octant)
    if (allocated(self%target_node_depth)) deallocate (self%target_node_depth)
    if (allocated(self%target_level_start)) deallocate (self%target_level_start)
    if (allocated(self%target_level_nodes)) deallocate (self%target_level_nodes)
    if (allocated(self%target_node_center)) deallocate (self%target_node_center)
    if (allocated(self%target_node_half_size)) deallocate (self%target_node_half_size)
    if (allocated(self%target_node_radius)) deallocate (self%target_node_radius)

    self%target_tree_ready = .false.
    self%target_nnode = 0_i32
    self%target_max_node = 0_i32
    self%target_node_max_depth = 0_i32
  end subroutine clear_target_tree

  module procedure rebuild_target_level_cache
    integer(i32) :: node_idx, depth
    integer(i32), allocatable :: depth_count(:), cursor(:)

    if (allocated(self%target_level_start)) deallocate (self%target_level_start)
    if (allocated(self%target_level_nodes)) deallocate (self%target_level_nodes)
    if (self%target_nnode <= 0_i32) return

    allocate (depth_count(self%target_node_max_depth + 1_i32), cursor(self%target_node_max_depth + 1_i32))
    depth_count = 0_i32
    do node_idx = 1_i32, self%target_nnode
      depth_count(self%target_node_depth(node_idx) + 1_i32) = depth_count(self%target_node_depth(node_idx) + 1_i32) + 1_i32
    end do

    allocate (self%target_level_start(self%target_node_max_depth + 2_i32), self%target_level_nodes(self%target_nnode))
    self%target_level_start(1) = 1_i32
    do depth = 0_i32, self%target_node_max_depth
      self%target_level_start(depth + 2_i32) = self%target_level_start(depth + 1_i32) + depth_count(depth + 1_i32)
    end do

    cursor = self%target_level_start(1:self%target_node_max_depth + 1_i32)
    do node_idx = 1_i32, self%target_nnode
      depth = self%target_node_depth(node_idx)
      self%target_level_nodes(cursor(depth + 1_i32)) = node_idx
      cursor(depth + 1_i32) = cursor(depth + 1_i32) + 1_i32
    end do

    deallocate (depth_count, cursor)
  end procedure rebuild_target_level_cache

  module procedure build_fmm_target_pair_index
    integer(i32) :: pair_idx, target_node, pos
    integer(i32), allocatable :: counts(:), cursor(:)

    if (allocated(self%fmm_m2l_target_start)) deallocate (self%fmm_m2l_target_start)
    if (allocated(self%fmm_m2l_pair_order)) deallocate (self%fmm_m2l_pair_order)
    if (n_target_nodes <= 0_i32) return

    allocate (self%fmm_m2l_target_start(n_target_nodes + 1_i32))
    allocate (self%fmm_m2l_pair_order(max(1_i32, self%fmm_m2l_pair_count)))
    self%fmm_m2l_target_start = 1_i32
    self%fmm_m2l_pair_order = 0_i32

    if (self%fmm_m2l_pair_count <= 0_i32) return

    allocate (counts(n_target_nodes), cursor(n_target_nodes))
    counts = 0_i32
    do pair_idx = 1_i32, self%fmm_m2l_pair_count
      target_node = self%fmm_m2l_target_nodes(pair_idx)
      counts(target_node) = counts(target_node) + 1_i32
    end do

    self%fmm_m2l_target_start(1) = 1_i32
    do target_node = 1_i32, n_target_nodes
      self%fmm_m2l_target_start(target_node + 1_i32) = self%fmm_m2l_target_start(target_node) + counts(target_node)
    end do

    cursor = self%fmm_m2l_target_start(1:n_target_nodes)
    do pair_idx = 1_i32, self%fmm_m2l_pair_count
      target_node = self%fmm_m2l_target_nodes(pair_idx)
      pos = cursor(target_node)
      self%fmm_m2l_pair_order(pos) = pair_idx
      cursor(target_node) = pos + 1_i32
    end do

    deallocate (counts, cursor)
  end procedure build_fmm_target_pair_index

  subroutine assign_source_leaves_as_targets_helper(self)
    class(field_solver_type), intent(inout) :: self
    integer(i32) :: node_idx, leaf_idx

    call clear_target_tree(self)

    self%nleaf = 0_i32
    do node_idx = 1_i32, self%nnode
      if (self%child_count(node_idx) <= 0_i32) self%nleaf = self%nleaf + 1_i32
    end do
    if (self%nleaf <= 0_i32) return

    if (allocated(self%leaf_nodes)) deallocate (self%leaf_nodes)
    if (allocated(self%leaf_slot_of_node)) deallocate (self%leaf_slot_of_node)
    allocate (self%leaf_nodes(self%nleaf), self%leaf_slot_of_node(self%max_node))
    self%leaf_slot_of_node = 0_i32

    leaf_idx = 0_i32
    do node_idx = 1_i32, self%nnode
      if (self%child_count(node_idx) > 0_i32) cycle
      leaf_idx = leaf_idx + 1_i32
      self%leaf_nodes(leaf_idx) = node_idx
      self%leaf_slot_of_node(node_idx) = leaf_idx
    end do
  end subroutine assign_source_leaves_as_targets_helper

  subroutine build_box_target_tree_helper(self, use_periodic_wrap)
    class(field_solver_type), intent(inout) :: self
    logical, intent(in) :: use_periodic_wrap
    integer(i32) :: depth_min, depth_max
    integer(i32) :: leaf_idx, node_idx, target_cap, level_nodes, depth
    real(dp) :: root_min(3), root_max(3)
    real(dp) :: src_center(3), src_half(3)

    src_center = self%node_center(:, 1_i32)
    src_half = self%node_half_size(:, 1_i32)

    root_min = self%target_box_min
    root_max = self%target_box_max
    if (any(root_max <= root_min)) then
      root_min = src_center - src_half
      root_max = src_center + src_half
    end if

    depth_min = 3_i32
    depth_max = 6_i32

    target_cap = 0_i32
    level_nodes = 1_i32
    do depth = 0_i32, depth_max
      target_cap = target_cap + level_nodes
      if (depth < depth_max) then
        if (level_nodes > huge(level_nodes) / 8_i32) error stop 'target tree capacity overflow.'
        level_nodes = 8_i32 * level_nodes
      end if
    end do

    call clear_target_tree(self)
    self%target_max_node = max(1024_i32, max(16_i32 * self%nnode, target_cap))
    allocate (self%target_child_count(self%target_max_node), self%target_child_idx(8, self%target_max_node))
    allocate (self%target_child_octant(8, self%target_max_node))
    allocate (self%target_node_depth(self%target_max_node))
    allocate (self%target_node_center(3, self%target_max_node), self%target_node_half_size(3, self%target_max_node))
    allocate (self%target_node_radius(self%target_max_node))
    self%target_child_count = 0_i32
    self%target_child_idx = 0_i32
    self%target_child_octant = 0_i32
    self%target_node_depth = 0_i32
    self%target_node_max_depth = 0_i32
    self%target_node_center = 0.0d0
    self%target_node_half_size = 0.0d0
    self%target_node_radius = 0.0d0

    self%target_nnode = 1_i32
    self%target_node_depth(1_i32) = 0_i32
    self%target_node_center(:, 1_i32) = 0.5d0 * (root_min + root_max)
    self%target_node_half_size(:, 1_i32) = 0.5d0 * (root_max - root_min)
    self%target_node_radius(1_i32) = sqrt(sum(self%target_node_half_size(:, 1_i32) * self%target_node_half_size(:, 1_i32)))
    call build_target_node_box_recursive(self, 1_i32, 0_i32, depth_min, depth_max, src_center, src_half, use_periodic_wrap)

    self%nleaf = 0_i32
    do node_idx = 1_i32, self%target_nnode
      if (self%target_child_count(node_idx) <= 0_i32) self%nleaf = self%nleaf + 1_i32
    end do
    if (self%nleaf <= 0_i32) return

    if (allocated(self%leaf_nodes)) deallocate (self%leaf_nodes)
    if (allocated(self%leaf_slot_of_node)) deallocate (self%leaf_slot_of_node)
    allocate (self%leaf_nodes(self%nleaf), self%leaf_slot_of_node(self%target_max_node))
    self%leaf_slot_of_node = 0_i32

    leaf_idx = 0_i32
    do node_idx = 1_i32, self%target_nnode
      if (self%target_child_count(node_idx) > 0_i32) cycle
      leaf_idx = leaf_idx + 1_i32
      self%leaf_nodes(leaf_idx) = node_idx
      self%leaf_slot_of_node(node_idx) = leaf_idx
    end do

    call rebuild_target_level_cache(self)
    self%target_tree_ready = .true.
  end subroutine build_box_target_tree_helper

  recursive subroutine build_target_node_box_recursive( &
    self, node_idx, depth, depth_min, depth_max, src_center, src_half, use_periodic_wrap &
  )
    class(field_solver_type), intent(inout) :: self
    integer(i32), intent(in) :: node_idx, depth, depth_min, depth_max
    real(dp), intent(in) :: src_center(3), src_half(3)
    logical, intent(in) :: use_periodic_wrap

    integer(i32) :: child_k, child_node, oct
    real(dp) :: child_half(3), offset(3)

    if (.not. target_should_split_box_helper( &
      self, node_idx, depth, depth_min, depth_max, src_center, src_half, use_periodic_wrap &
    )) return

    child_half = 0.5d0 * self%target_node_half_size(:, node_idx)
    if (maxval(child_half) <= tiny(1.0d0)) return

    child_k = 0_i32
    do oct = 1_i32, 8_i32
      self%target_nnode = self%target_nnode + 1_i32
      if (self%target_nnode > self%target_max_node) error stop 'target tree capacity exceeded.'
      child_node = self%target_nnode
      child_k = child_k + 1_i32
      self%target_child_idx(child_k, node_idx) = child_node
      self%target_child_octant(child_k, node_idx) = oct
      self%target_node_depth(child_node) = depth + 1_i32
      self%target_node_max_depth = max(self%target_node_max_depth, depth + 1_i32)

      offset(1) = merge(-1.0d0, 1.0d0, iand(oct - 1_i32, 1_i32) == 0_i32)
      offset(2) = merge(-1.0d0, 1.0d0, iand(oct - 1_i32, 2_i32) == 0_i32)
      offset(3) = merge(-1.0d0, 1.0d0, iand(oct - 1_i32, 4_i32) == 0_i32)

      self%target_node_half_size(:, child_node) = child_half
      self%target_node_center(:, child_node) = self%target_node_center(:, node_idx) + offset * child_half
      self%target_node_radius(child_node) = sqrt(sum(child_half * child_half))

      call build_target_node_box_recursive( &
        self, child_node, depth + 1_i32, depth_min, depth_max, src_center, src_half, use_periodic_wrap &
      )
    end do
    self%target_child_count(node_idx) = child_k
  end subroutine build_target_node_box_recursive

  logical function target_should_split_box_helper( &
    self, node_idx, depth, depth_min, depth_max, src_center, src_half, use_periodic_wrap &
  )
    class(field_solver_type), intent(in) :: self
    integer(i32), intent(in) :: node_idx, depth, depth_min, depth_max
    real(dp), intent(in) :: src_center(3), src_half(3)
    logical, intent(in) :: use_periodic_wrap
    real(dp) :: dsrc, radius

    if (depth < depth_min) then
      target_should_split_box_helper = .true.
      return
    end if
    if (depth >= depth_max) then
      target_should_split_box_helper = .false.
      return
    end if

    radius = self%target_node_radius(node_idx)
    if (use_periodic_wrap) then
      dsrc = distance_to_source_bbox_periodic_helper(self, self%target_node_center(:, node_idx), src_center, src_half)
    else
      dsrc = distance_to_source_bbox_helper(self%target_node_center(:, node_idx), src_center, src_half)
    end if
    target_should_split_box_helper = (radius > 0.35d0 * max(dsrc, 1.0d-12))
  end function target_should_split_box_helper

  real(dp) function distance_to_source_bbox_helper(p, src_center, src_half)
    real(dp), intent(in) :: p(3), src_center(3), src_half(3)
    real(dp) :: d(3), q(3)
    integer(i32) :: axis

    d = p - src_center
    do axis = 1_i32, 3_i32
      q(axis) = max(0.0d0, abs(d(axis)) - src_half(axis))
    end do
    distance_to_source_bbox_helper = sqrt(sum(q * q))
  end function distance_to_source_bbox_helper

  real(dp) function distance_to_source_bbox_periodic_helper(self, p, src_center, src_half)
    class(field_solver_type), intent(in) :: self
    real(dp), intent(in) :: p(3), src_center(3), src_half(3)
    integer(i32) :: axis
    real(dp) :: d(3), q(3)

    d = p - src_center
    call apply_periodic2_minimum_image(self, d)
    do axis = 1_i32, 3_i32
      q(axis) = max(0.0d0, abs(d(axis)) - src_half(axis))
    end do
    distance_to_source_bbox_periodic_helper = sqrt(sum(q * q))
  end function distance_to_source_bbox_periodic_helper

  subroutine append_i32_buffer(buf, n_used, capacity, value)
    integer(i32), allocatable, intent(inout) :: buf(:)
    integer(i32), intent(inout) :: n_used
    integer(i32), intent(inout) :: capacity
    integer(i32), intent(in) :: value
    integer(i32), allocatable :: tmp(:)
    integer(i32) :: new_capacity

    if (n_used >= capacity) then
      new_capacity = max(capacity * 2_i32, capacity + 32_i32)
      allocate (tmp(new_capacity))
      tmp = 0_i32
      if (n_used > 0_i32) tmp(1:n_used) = buf(1:n_used)
      call move_alloc(tmp, buf)
      capacity = new_capacity
    end if
    n_used = n_used + 1_i32
    buf(n_used) = value
  end subroutine append_i32_buffer

  subroutine setup_periodic_image_loop(self, img_min, img_max, axis1, axis2)
    class(field_solver_type), intent(in) :: self
    integer(i32), intent(out) :: img_min, img_max, axis1, axis2

    img_min = 0_i32
    img_max = 0_i32
    axis1 = 0_i32
    axis2 = 0_i32
    if (self%use_periodic2) then
      img_min = -self%periodic_image_layers
      img_max = self%periodic_image_layers
      axis1 = self%periodic_axes(1)
      axis2 = self%periodic_axes(2)
    end if
  end subroutine setup_periodic_image_loop

  pure subroutine apply_periodic2_minimum_image(self, d)
    class(field_solver_type), intent(in) :: self
    real(dp), intent(inout) :: d(3)
    integer(i32) :: axis, k
    real(dp) :: len_k, half_len

    if (.not. self%use_periodic2) return
    do k = 1_i32, 2_i32
      axis = self%periodic_axes(k)
      len_k = self%periodic_len(k)
      half_len = 0.5d0 * len_k
      if (d(axis) > half_len) then
        d(axis) = d(axis) - len_k
      else if (d(axis) < -half_len) then
        d(axis) = d(axis) + len_k
      end if
    end do
  end subroutine apply_periodic2_minimum_image

  subroutine build_periodic_shift_values(self, shift_axis1, shift_axis2, nshift)
    class(field_solver_type), intent(in) :: self
    real(dp), intent(out) :: shift_axis1(:), shift_axis2(:)
    integer(i32), intent(out) :: nshift
    integer(i32) :: s, img

    nshift = 1_i32
    if (size(shift_axis1) < 1_i32 .or. size(shift_axis2) < 1_i32) then
      error stop 'periodic shift buffer is empty.'
    end if
    shift_axis1(1) = 0.0d0
    shift_axis2(1) = 0.0d0
    if (.not. self%use_periodic2) return

    nshift = 2_i32 * self%periodic_image_layers + 1_i32
    if (size(shift_axis1) < nshift .or. size(shift_axis2) < nshift) then
      error stop 'periodic shift buffer is too small.'
    end if
    do s = 1_i32, nshift
      img = s - self%periodic_image_layers - 1_i32
      shift_axis1(s) = real(img, dp) * self%periodic_len(1)
      shift_axis2(s) = real(img, dp) * self%periodic_len(2)
    end do
  end subroutine build_periodic_shift_values

  subroutine add_point_charge_images_field(q, src, target, soft2, axis1, axis2, shift_axis1, shift_axis2, nshift, ex, ey, ez)
    real(dp), intent(in) :: q
    real(dp), intent(in) :: src(3), target(3)
    real(dp), intent(in) :: soft2
    integer(i32), intent(in) :: axis1, axis2
    real(dp), intent(in) :: shift_axis1(:), shift_axis2(:)
    integer(i32), intent(in) :: nshift
    real(dp), intent(inout) :: ex, ey, ez

    integer(i32) :: img_i, img_j
    real(dp) :: shift1, shift2
    real(dp) :: src_x, src_y, src_z
    real(dp) :: dx, dy, dz, r2, inv_r3

    if (abs(q) <= tiny(1.0d0)) return

    do img_i = 1_i32, nshift
      shift1 = shift_axis1(img_i)
      do img_j = 1_i32, nshift
        shift2 = shift_axis2(img_j)

        src_x = src(1)
        src_y = src(2)
        src_z = src(3)
        if (axis1 == 1_i32) src_x = src_x + shift1
        if (axis1 == 2_i32) src_y = src_y + shift1
        if (axis1 == 3_i32) src_z = src_z + shift1
        if (axis2 == 1_i32) src_x = src_x + shift2
        if (axis2 == 2_i32) src_y = src_y + shift2
        if (axis2 == 3_i32) src_z = src_z + shift2

        dx = target(1) - src_x
        dy = target(2) - src_y
        dz = target(3) - src_z
        r2 = dx * dx + dy * dy + dz * dz + soft2
        inv_r3 = 1.0d0 / (sqrt(r2) * r2)
        ex = ex + q * inv_r3 * dx
        ey = ey + q * inv_r3 * dy
        ez = ez + q * inv_r3 * dz
      end do
    end do
  end subroutine add_point_charge_images_field

  integer(i32) function active_tree_nnode(self, use_target_tree)
    class(field_solver_type), intent(in) :: self
    logical, intent(in) :: use_target_tree

    active_tree_nnode = merge(self%target_nnode, self%nnode, use_target_tree)
  end function active_tree_nnode

  integer(i32) function active_tree_child_count(self, use_target_tree, node_idx)
    class(field_solver_type), intent(in) :: self
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: node_idx

    if (use_target_tree) then
      active_tree_child_count = self%target_child_count(node_idx)
    else
      active_tree_child_count = self%child_count(node_idx)
    end if
  end function active_tree_child_count

  integer(i32) function active_tree_child_idx(self, use_target_tree, child_k, node_idx)
    class(field_solver_type), intent(in) :: self
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: child_k, node_idx

    if (use_target_tree) then
      active_tree_child_idx = self%target_child_idx(child_k, node_idx)
    else
      active_tree_child_idx = self%child_idx(child_k, node_idx)
    end if
  end function active_tree_child_idx

  integer(i32) function active_tree_child_octant(self, use_target_tree, child_k, node_idx)
    class(field_solver_type), intent(in) :: self
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: child_k, node_idx

    if (use_target_tree) then
      active_tree_child_octant = self%target_child_octant(child_k, node_idx)
    else
      active_tree_child_octant = self%child_octant(child_k, node_idx)
    end if
  end function active_tree_child_octant

  function active_tree_node_center(self, use_target_tree, node_idx) result(center)
    class(field_solver_type), intent(in) :: self
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: node_idx
    real(dp) :: center(3)

    if (use_target_tree) then
      center = self%target_node_center(:, node_idx)
    else
      center = self%node_center(:, node_idx)
    end if
  end function active_tree_node_center

  function active_tree_node_half_size(self, use_target_tree, node_idx) result(half_size)
    class(field_solver_type), intent(in) :: self
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: node_idx
    real(dp) :: half_size(3)

    if (use_target_tree) then
      half_size = self%target_node_half_size(:, node_idx)
    else
      half_size = self%node_half_size(:, node_idx)
    end if
  end function active_tree_node_half_size

  real(dp) function active_tree_node_radius(self, use_target_tree, node_idx)
    class(field_solver_type), intent(in) :: self
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: node_idx

    if (use_target_tree) then
      active_tree_node_radius = self%target_node_radius(node_idx)
    else
      active_tree_node_radius = self%node_radius(node_idx)
    end if
  end function active_tree_node_radius
  !> M2L と L2L で葉ノード局所展開（E/Jac/Hess）を更新する。
  module procedure refresh_fmm_locals
    integer(i32) :: leaf_slot, node_idx, parent_node
    integer(i32) :: i, j, k
    integer(i32) :: axis1, axis2
    integer(i32) :: t_node
    integer(i32) :: depth, max_depth
    integer(i32) :: level_pos, level_start_pos, level_end_pos
    logical :: use_target_tree
    real(dp) :: time_t0, time_t1
    real(dp) :: dr(3), t_center(3), p_center(3)
    real(dp) :: soft2
    integer(i32) :: n_target_nodes, nshift

    if (trim(self%mode) /= 'fmm') return
    if (.not. self%fmm_ready) return
    if (self%nleaf <= 0_i32) return

    use_target_tree = self%target_tree_ready
    n_target_nodes = active_tree_nnode(self, use_target_tree)
    if (.not. allocated(self%fmm_parent_of)) error stop 'FMM parent cache is not allocated.'
    if (.not. allocated(self%fmm_shift_axis1) .or. .not. allocated(self%fmm_shift_axis2)) then
      error stop 'FMM periodic shift cache is not allocated.'
    end if
    if (.not. allocated(self%fmm_m2l_target_start) .or. &
        .not. allocated(self%fmm_m2l_pair_order)) then
      error stop 'FMM target pair index is not allocated.'
    end if
    if (.not. allocated(self%fmm_node_local_e0) .or. .not. allocated(self%fmm_node_local_jac) .or. &
        .not. allocated(self%fmm_node_local_hess)) then
      error stop 'FMM local work arrays are not allocated.'
    end if

    soft2 = self%softening * self%softening
    axis1 = 0_i32
    axis2 = 0_i32
    if (self%use_periodic2) then
      axis1 = self%periodic_axes(1)
      axis2 = self%periodic_axes(2)
    end if
    nshift = size(self%fmm_shift_axis1)

    time_t0 = field_solver_time_seconds()
    !$omp parallel default(none) shared(self)
    !$omp workshare
    self%fmm_node_local_e0 = 0.0d0
    self%fmm_node_local_jac = 0.0d0
    self%fmm_node_local_hess = 0.0d0
    !$omp end workshare
    !$omp end parallel
    time_t1 = field_solver_time_seconds()
    self%fmm_last_clear_time_s = time_t1 - time_t0

    time_t0 = field_solver_time_seconds()
    !$omp parallel do default(none) schedule(guided) &
    !$omp shared(self,n_target_nodes,use_target_tree,axis1,axis2,soft2,nshift) private(t_node)
    do t_node = 1_i32, n_target_nodes
      call accumulate_target_m2l(t_node)
    end do
    !$omp end parallel do
    time_t1 = field_solver_time_seconds()
    self%fmm_last_m2l_time_s = time_t1 - time_t0

    ! L2L: 親ノード局所展開を子ノード中心へ平行移動する。
    if (use_target_tree) then
      if (.not. allocated(self%target_level_start) .or. .not. allocated(self%target_level_nodes)) then
        error stop 'FMM target level cache is not allocated.'
      end if
      max_depth = self%target_node_max_depth
    else
      if (.not. allocated(self%node_level_start) .or. .not. allocated(self%node_level_nodes)) then
        error stop 'FMM source level cache is not allocated.'
      end if
      max_depth = self%node_max_depth
    end if

    time_t0 = field_solver_time_seconds()
    do depth = 1_i32, max_depth
      if (use_target_tree) then
        level_start_pos = self%target_level_start(depth + 1_i32)
        level_end_pos = self%target_level_start(depth + 2_i32) - 1_i32
      else
        level_start_pos = self%node_level_start(depth + 1_i32)
        level_end_pos = self%node_level_start(depth + 2_i32) - 1_i32
      end if
      !$omp parallel do default(none) schedule(static) &
      !$omp shared(self,use_target_tree,level_start_pos,level_end_pos) &
      !$omp private(level_pos,node_idx,parent_node,dr,t_center,p_center,i,j,k)
      do level_pos = level_start_pos, level_end_pos
        if (use_target_tree) then
          node_idx = self%target_level_nodes(level_pos)
        else
          node_idx = self%node_level_nodes(level_pos)
        end if
        parent_node = self%fmm_parent_of(node_idx)
        if (parent_node <= 0_i32) cycle
        if (use_target_tree) then
          t_center = self%target_node_center(:, node_idx)
          p_center = self%target_node_center(:, parent_node)
        else
          t_center = self%node_center(:, node_idx)
          p_center = self%node_center(:, parent_node)
        end if
        dr = t_center - p_center
        do i = 1_i32, 3_i32
          self%fmm_node_local_e0(i, node_idx) = self%fmm_node_local_e0(i, node_idx) + self%fmm_node_local_e0(i, parent_node)
          do j = 1_i32, 3_i32
            self%fmm_node_local_e0(i, node_idx) = self%fmm_node_local_e0(i, node_idx) &
                                                  + self%fmm_node_local_jac(i, j, parent_node) * dr(j)
            do k = 1_i32, 3_i32
              self%fmm_node_local_e0(i, node_idx) = self%fmm_node_local_e0(i, node_idx) &
                                                    + 0.5d0 * self%fmm_node_local_hess(i, j, k, parent_node) * dr(j) * dr(k)
              self%fmm_node_local_jac(i, j, node_idx) = self%fmm_node_local_jac(i, j, node_idx) &
                                                        + self%fmm_node_local_hess(i, j, k, parent_node) * dr(k)
            end do
            self%fmm_node_local_jac(i, j, node_idx) = self%fmm_node_local_jac(i, j, node_idx) &
                                                      + self%fmm_node_local_jac(i, j, parent_node)
            do k = 1_i32, 3_i32
              self%fmm_node_local_hess(i, j, k, node_idx) = self%fmm_node_local_hess(i, j, k, node_idx) &
                                                            + self%fmm_node_local_hess(i, j, k, parent_node)
            end do
          end do
        end do
      end do
      !$omp end parallel do
    end do
    time_t1 = field_solver_time_seconds()
    self%fmm_last_l2l_time_s = time_t1 - time_t0

    time_t0 = field_solver_time_seconds()
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(self) private(leaf_slot,node_idx)
    do leaf_slot = 1_i32, self%nleaf
      node_idx = self%leaf_nodes(leaf_slot)
      self%leaf_far_e0(:, leaf_slot) = self%fmm_node_local_e0(:, node_idx)
      self%leaf_far_jac(:, :, leaf_slot) = self%fmm_node_local_jac(:, :, node_idx)
      self%leaf_far_hess(:, :, :, leaf_slot) = self%fmm_node_local_hess(:, :, :, node_idx)
    end do
    !$omp end parallel do
    time_t1 = field_solver_time_seconds()
    self%fmm_last_copy_time_s = time_t1 - time_t0

  contains

    subroutine accumulate_target_m2l(t_node)
      integer(i32), intent(in) :: t_node
      integer(i32) :: pair_pos, pair_idx, s_node
      integer(i32) :: img_i, img_j
      real(dp) :: qi
      real(dp) :: tx, ty, tz
      real(dp) :: base_x, base_y, base_z
      real(dp) :: dx, dy, dz
      real(dp) :: shift1, shift2
      real(dp) :: r2, inv_r, inv_r2
      real(dp) :: s3, s5, s7
      real(dp) :: xx, yy, zz, xy, xz, yz
      real(dp) :: axx, ayy, azz
      real(dp) :: e0x, e0y, e0z
      real(dp) :: jxx, jyy, jzz, jxy, jxz, jyz
      real(dp) :: hxxx, hxxy, hxxz, hxyy, hxyz, hxzz, hyyy, hyyz, hyzz, hzzz

      if (use_target_tree) then
        tx = self%target_node_center(1, t_node)
        ty = self%target_node_center(2, t_node)
        tz = self%target_node_center(3, t_node)
      else
        tx = self%node_center(1, t_node)
        ty = self%node_center(2, t_node)
        tz = self%node_center(3, t_node)
      end if

      e0x = 0.0d0
      e0y = 0.0d0
      e0z = 0.0d0
      jxx = 0.0d0
      jyy = 0.0d0
      jzz = 0.0d0
      jxy = 0.0d0
      jxz = 0.0d0
      jyz = 0.0d0
      hxxx = 0.0d0
      hxxy = 0.0d0
      hxxz = 0.0d0
      hxyy = 0.0d0
      hxyz = 0.0d0
      hxzz = 0.0d0
      hyyy = 0.0d0
      hyyz = 0.0d0
      hyzz = 0.0d0
      hzzz = 0.0d0

      do pair_pos = self%fmm_m2l_target_start(t_node), self%fmm_m2l_target_start(t_node + 1_i32) - 1_i32
        pair_idx = self%fmm_m2l_pair_order(pair_pos)
        s_node = self%fmm_m2l_source_nodes(pair_idx)
        qi = self%node_q(s_node)
        if (abs(qi) <= tiny(1.0d0)) cycle

        base_x = tx - self%node_charge_center(1, s_node)
        base_y = ty - self%node_charge_center(2, s_node)
        base_z = tz - self%node_charge_center(3, s_node)

        do img_i = 1_i32, nshift
          shift1 = self%fmm_shift_axis1(img_i)
          do img_j = 1_i32, nshift
            shift2 = self%fmm_shift_axis2(img_j)
            dx = base_x
            dy = base_y
            dz = base_z
            if (axis1 == 1_i32) dx = dx - shift1
            if (axis1 == 2_i32) dy = dy - shift1
            if (axis1 == 3_i32) dz = dz - shift1
            if (axis2 == 1_i32) dx = dx - shift2
            if (axis2 == 2_i32) dy = dy - shift2
            if (axis2 == 3_i32) dz = dz - shift2

            r2 = dx * dx + dy * dy + dz * dz + soft2
            inv_r = 1.0d0 / sqrt(r2)
            inv_r2 = 1.0d0 / r2
            s3 = qi * inv_r * inv_r2
            s5 = -3.0d0 * s3 * inv_r2
            s7 = -5.0d0 * s5 * inv_r2

            e0x = e0x + s3 * dx
            e0y = e0y + s3 * dy
            e0z = e0z + s3 * dz

            xx = dx * dx
            yy = dy * dy
            zz = dz * dz
            xy = dx * dy
            xz = dx * dz
            yz = dy * dz

            jxx = jxx + s3 + s5 * xx
            jyy = jyy + s3 + s5 * yy
            jzz = jzz + s3 + s5 * zz
            jxy = jxy + s5 * xy
            jxz = jxz + s5 * xz
            jyz = jyz + s5 * yz

            axx = s7 * xx + s5
            ayy = s7 * yy + s5
            azz = s7 * zz + s5

            hxxx = hxxx + dx * (axx + 2.0d0 * s5)
            hxxy = hxxy + dy * axx
            hxxz = hxxz + dz * axx
            hxyy = hxyy + dx * ayy
            hxyz = hxyz + s7 * dx * dy * dz
            hxzz = hxzz + dx * azz
            hyyy = hyyy + dy * (ayy + 2.0d0 * s5)
            hyyz = hyyz + dz * ayy
            hyzz = hyzz + dy * azz
            hzzz = hzzz + dz * (azz + 2.0d0 * s5)
          end do
        end do
      end do

      self%fmm_node_local_e0(1, t_node) = e0x
      self%fmm_node_local_e0(2, t_node) = e0y
      self%fmm_node_local_e0(3, t_node) = e0z

      self%fmm_node_local_jac(1, 1, t_node) = jxx
      self%fmm_node_local_jac(2, 2, t_node) = jyy
      self%fmm_node_local_jac(3, 3, t_node) = jzz
      self%fmm_node_local_jac(1, 2, t_node) = jxy
      self%fmm_node_local_jac(2, 1, t_node) = jxy
      self%fmm_node_local_jac(1, 3, t_node) = jxz
      self%fmm_node_local_jac(3, 1, t_node) = jxz
      self%fmm_node_local_jac(2, 3, t_node) = jyz
      self%fmm_node_local_jac(3, 2, t_node) = jyz

      self%fmm_node_local_hess(1, 1, 1, t_node) = hxxx
      self%fmm_node_local_hess(1, 1, 2, t_node) = hxxy
      self%fmm_node_local_hess(1, 2, 1, t_node) = hxxy
      self%fmm_node_local_hess(2, 1, 1, t_node) = hxxy
      self%fmm_node_local_hess(1, 1, 3, t_node) = hxxz
      self%fmm_node_local_hess(1, 3, 1, t_node) = hxxz
      self%fmm_node_local_hess(3, 1, 1, t_node) = hxxz
      self%fmm_node_local_hess(1, 2, 2, t_node) = hxyy
      self%fmm_node_local_hess(2, 1, 2, t_node) = hxyy
      self%fmm_node_local_hess(2, 2, 1, t_node) = hxyy
      self%fmm_node_local_hess(1, 2, 3, t_node) = hxyz
      self%fmm_node_local_hess(1, 3, 2, t_node) = hxyz
      self%fmm_node_local_hess(2, 1, 3, t_node) = hxyz
      self%fmm_node_local_hess(2, 3, 1, t_node) = hxyz
      self%fmm_node_local_hess(3, 1, 2, t_node) = hxyz
      self%fmm_node_local_hess(3, 2, 1, t_node) = hxyz
      self%fmm_node_local_hess(1, 3, 3, t_node) = hxzz
      self%fmm_node_local_hess(3, 1, 3, t_node) = hxzz
      self%fmm_node_local_hess(3, 3, 1, t_node) = hxzz
      self%fmm_node_local_hess(2, 2, 2, t_node) = hyyy
      self%fmm_node_local_hess(2, 2, 3, t_node) = hyyz
      self%fmm_node_local_hess(2, 3, 2, t_node) = hyyz
      self%fmm_node_local_hess(3, 2, 2, t_node) = hyyz
      self%fmm_node_local_hess(2, 3, 3, t_node) = hyzz
      self%fmm_node_local_hess(3, 2, 3, t_node) = hyzz
      self%fmm_node_local_hess(3, 3, 2, t_node) = hyzz
      self%fmm_node_local_hess(3, 3, 3, t_node) = hzzz
    end subroutine accumulate_target_m2l
  end procedure refresh_fmm_locals
  !> 観測点の octant を辿って対応する葉ノードを返す。
  module procedure locate_target_leaf
    integer(i32) :: oct, child_k, child, cand
    integer(i32) :: axis
    logical :: use_target_tree
    real(dp) :: best_dist2, dist2, dx, dy, dz
    real(dp) :: extent_eps
    real(dp) :: r_eval(3)
    real(dp) :: root_center(3), root_half(3), node_center_curr(3), cand_center(3)

    node_idx = 0_i32
    if (.not. self%tree_ready .or. self%nnode <= 0_i32) return

    use_target_tree = self%target_tree_ready
    if (active_tree_nnode(self, use_target_tree) <= 0_i32) return

    r_eval = r
    if (use_target_tree .and. self%use_periodic2) call wrap_periodic2_point(self, r_eval)

    if (use_target_tree .or. .not. self%use_periodic2) then
      root_center = active_tree_node_center(self, use_target_tree, 1_i32)
      root_half = active_tree_node_half_size(self, use_target_tree, 1_i32)
      extent_eps = 1.0d-12 * max(1.0d0, maxval(abs(root_center)))
      do axis = 1_i32, 3_i32
        if (abs(r_eval(axis) - root_center(axis)) > root_half(axis) + extent_eps) return
      end do
    end if

    node_idx = 1_i32
    do while (active_tree_child_count(self, use_target_tree, node_idx) > 0_i32)
      node_center_curr = active_tree_node_center(self, use_target_tree, node_idx)
      oct = octant_index(r_eval(1), r_eval(2), r_eval(3), node_center_curr)
      child = 0_i32
      do child_k = 1_i32, active_tree_child_count(self, use_target_tree, node_idx)
        if (active_tree_child_octant(self, use_target_tree, child_k, node_idx) == oct) then
          child = active_tree_child_idx(self, use_target_tree, child_k, node_idx)
          exit
        end if
      end do

      if (child <= 0_i32) then
        child = active_tree_child_idx(self, use_target_tree, 1_i32, node_idx)
        best_dist2 = huge(1.0d0)
        do child_k = 1_i32, active_tree_child_count(self, use_target_tree, node_idx)
          cand = active_tree_child_idx(self, use_target_tree, child_k, node_idx)
          cand_center = active_tree_node_center(self, use_target_tree, cand)
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
      node_idx = child
    end do
  end procedure locate_target_leaf

  subroutine wrap_periodic2_point(self, p)
    class(field_solver_type), intent(in) :: self
    real(dp), intent(inout) :: p(3)
    integer(i32) :: k, axis

    if (.not. self%use_periodic2) return
    do k = 1_i32, 2_i32
      axis = self%periodic_axes(k)
      p(axis) = self%target_box_min(axis) + modulo(p(axis) - self%target_box_min(axis), self%periodic_len(k))
    end do
  end subroutine wrap_periodic2_point

  subroutine eval_fmm_fallback_target(self, mesh, rt, soft2, ex, ey, ez)
    class(field_solver_type), intent(in) :: self
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: rt(3), soft2
    real(dp), intent(out) :: ex, ey, ez

    if (self%use_periodic2) then
      call eval_periodic2_tree_target(self, mesh, rt, ex, ey, ez)
      if (trim(self%periodic_far_correction) == 'ewald_like') then
        call add_periodic2_ewald_like_correction_all_leaves(self, rt, ex, ey, ez)
      end if
    else
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      call traverse_node(self, mesh, 1_i32, rt(1), rt(2), rt(3), soft2, ex, ey, ez)
    end if
  end subroutine eval_fmm_fallback_target

  !> FMM 局所展開 + 近傍 direct 和で観測点の電場を返す。
  module procedure eval_e_fmm
    integer(i32) :: leaf_node, leaf_slot
    integer(i32) :: near_pos, near_node, p, idx, p_end
    integer(i32) :: nshift
    integer(i32) :: axis1, axis2
    integer(i32) :: j, k
    real(dp) :: soft2, qi
    real(dp) :: ex, ey, ez
    real(dp) :: dr(3), center(3), rt(3), src(3)
    real(dp) :: shift_axis1(2_i32 * max(0_i32, self%periodic_image_layers) + 1_i32)
    real(dp) :: shift_axis2(2_i32 * max(0_i32, self%periodic_image_layers) + 1_i32)

    rt = r
    if (self%use_periodic2) call wrap_periodic2_point(self, rt)

    soft2 = self%softening * self%softening
    leaf_node = locate_target_leaf(self, rt)
    if (leaf_node <= 0_i32) then
      call eval_fmm_fallback_target(self, mesh, rt, soft2, ex, ey, ez)
      e(1) = k_coulomb * ex
      e(2) = k_coulomb * ey
      e(3) = k_coulomb * ez
      return
    end if
    leaf_slot = self%leaf_slot_of_node(leaf_node)
    if (leaf_slot <= 0_i32) then
      call eval_fmm_fallback_target(self, mesh, rt, soft2, ex, ey, ez)
      e(1) = k_coulomb * ex
      e(2) = k_coulomb * ey
      e(3) = k_coulomb * ez
      return
    end if

    axis1 = 0_i32
    axis2 = 0_i32
    if (self%use_periodic2) then
      axis1 = self%periodic_axes(1)
      axis2 = self%periodic_axes(2)
    end if
    call build_periodic_shift_values(self, shift_axis1, shift_axis2, nshift)
    center = active_tree_node_center(self, self%target_tree_ready, leaf_node)

    ex = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    dr = rt - center
    ex = ex + self%leaf_far_e0(1, leaf_slot)
    ey = ey + self%leaf_far_e0(2, leaf_slot)
    ez = ez + self%leaf_far_e0(3, leaf_slot)

    do j = 1_i32, 3_i32
      ex = ex + self%leaf_far_jac(1, j, leaf_slot) * dr(j)
      ey = ey + self%leaf_far_jac(2, j, leaf_slot) * dr(j)
      ez = ez + self%leaf_far_jac(3, j, leaf_slot) * dr(j)
    end do

    do j = 1_i32, 3_i32
      do k = 1_i32, 3_i32
        ex = ex + 0.5d0 * self%leaf_far_hess(1, j, k, leaf_slot) * dr(j) * dr(k)
        ey = ey + 0.5d0 * self%leaf_far_hess(2, j, k, leaf_slot) * dr(j) * dr(k)
        ez = ez + 0.5d0 * self%leaf_far_hess(3, j, k, leaf_slot) * dr(j) * dr(k)
      end do
    end do

    do near_pos = self%near_start(leaf_slot), self%near_start(leaf_slot + 1_i32) - 1_i32
      near_node = self%near_nodes(near_pos)
      p_end = self%node_start(near_node) + self%node_count(near_node) - 1_i32
      do p = self%node_start(near_node), p_end
        idx = self%elem_order(p)
        qi = mesh%q_elem(idx)
        if (abs(qi) <= tiny(1.0d0)) cycle
        src(1) = mesh%center_x(idx)
        src(2) = mesh%center_y(idx)
        src(3) = mesh%center_z(idx)
        call add_point_charge_images_field(qi, src, rt, soft2, axis1, axis2, shift_axis1, shift_axis2, nshift, ex, ey, ez)
      end do
    end do
    if (self%use_periodic2 .and. trim(self%periodic_far_correction) == 'ewald_like') then
      call add_periodic2_ewald_like_correction(self, leaf_slot, rt, ex, ey, ez)
    end if

    e(1) = k_coulomb * ex
    e(2) = k_coulomb * ey
    e(3) = k_coulomb * ez
  end procedure eval_e_fmm

  !> periodic2 で target がtarget-tree領域外にある場合の M2P/treecode フォールバック。
  subroutine eval_periodic2_tree_target(self, mesh, r, ex, ey, ez)
    class(field_solver_type), intent(in) :: self
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: r(3)
    real(dp), intent(out) :: ex, ey, ez

    integer(i32) :: axis1, axis2, nshift
    real(dp) :: soft2, rp(3)
    real(dp) :: shift_axis1(2_i32 * max(0_i32, self%periodic_image_layers) + 1_i32)
    real(dp) :: shift_axis2(2_i32 * max(0_i32, self%periodic_image_layers) + 1_i32)

    axis1 = self%periodic_axes(1)
    axis2 = self%periodic_axes(2)
    call build_periodic_shift_values(self, shift_axis1, shift_axis2, nshift)
    soft2 = self%softening * self%softening
    rp = r
    call wrap_periodic2_point(self, rp)

    ex = 0.0d0
    ey = 0.0d0
    ez = 0.0d0

    call accumulate_node(1_i32)

  contains

    recursive subroutine accumulate_node(node_idx)
      integer(i32), intent(in) :: node_idx
      integer(i32) :: child_k, p, idx, p_end
      real(dp) :: qi
      real(dp) :: src(3)

      if (self%child_count(node_idx) <= 0_i32) then
        p_end = self%node_start(node_idx) + self%node_count(node_idx) - 1_i32
        do p = self%node_start(node_idx), p_end
          idx = self%elem_order(p)
          qi = mesh%q_elem(idx)
          if (abs(qi) <= tiny(1.0d0)) cycle
          src(1) = mesh%center_x(idx)
          src(2) = mesh%center_y(idx)
          src(3) = mesh%center_z(idx)
          call add_point_charge_images_field(qi, src, rp, soft2, axis1, axis2, shift_axis1, shift_axis2, nshift, ex, ey, ez)
        end do
        return
      end if

      if (accept_node_periodic_target(node_idx)) then
        qi = self%node_q(node_idx)
        if (abs(qi) <= tiny(1.0d0)) return
        src = self%node_charge_center(:, node_idx)
        call add_point_charge_images_field(qi, src, rp, soft2, axis1, axis2, shift_axis1, shift_axis2, nshift, ex, ey, ez)
        return
      end if

      do child_k = 1_i32, self%child_count(node_idx)
        call accumulate_node(self%child_idx(child_k, node_idx))
      end do
    end subroutine accumulate_node

    logical function accept_node_periodic_target(node_idx)
      integer(i32), intent(in) :: node_idx
      real(dp) :: d(3), dist2, dist, radius, theta_eff

      d = rp - self%node_center(:, node_idx)
      call apply_periodic2_minimum_image(self, d)
      dist2 = sum(d * d)
      if (dist2 <= 0.0d0) then
        accept_node_periodic_target = .false.
        return
      end if

      radius = self%node_radius(node_idx)
      dist = sqrt(dist2)
      if (dist <= radius) then
        accept_node_periodic_target = .false.
        return
      end if

      theta_eff = 0.2d0 * self%theta
      if (radius >= theta_eff * (dist - radius)) then
        accept_node_periodic_target = .false.
        return
      end if

      if (self%node_abs_q(node_idx) <= 0.0d0) then
        accept_node_periodic_target = .true.
      else
        accept_node_periodic_target = abs(self%node_q(node_idx)) >= charge_cancellation_tol * self%node_abs_q(node_idx)
      end if
    end function accept_node_periodic_target
  end subroutine eval_periodic2_tree_target
  !> periodic2 の近傍画像和で欠落する遠方セル寄与を、erfc スクリーン核で補正する。
  subroutine add_periodic2_ewald_like_correction(self, leaf_slot, r, ex, ey, ez)
    class(field_solver_type), intent(in) :: self
    integer(i32), intent(in) :: leaf_slot
    real(dp), intent(in) :: r(3)
    real(dp), intent(inout) :: ex, ey, ez

    integer(i32) :: node_pos, node_idx
    integer(i32) :: nimg, img_outer
    integer(i32) :: axis1, axis2
    real(dp) :: alpha

    if (.not. prepare_periodic2_ewald(self, alpha, nimg, img_outer, axis1, axis2)) return

    do node_pos = self%far_start(leaf_slot), self%far_start(leaf_slot + 1_i32) - 1_i32
      node_idx = self%far_nodes(node_pos)
      call add_ewald_node_correction(self, node_idx, r, axis1, axis2, nimg, img_outer, alpha, ex, ey, ez)
    end do

    do node_pos = self%near_start(leaf_slot), self%near_start(leaf_slot + 1_i32) - 1_i32
      node_idx = self%near_nodes(node_pos)
      call add_ewald_node_correction(self, node_idx, r, axis1, axis2, nimg, img_outer, alpha, ex, ey, ez)
    end do
  end subroutine add_periodic2_ewald_like_correction

  !> source 葉ノード全体を代表点として、periodic2 の遠方画像補正を加算する。
  subroutine add_periodic2_ewald_like_correction_all_leaves(self, r, ex, ey, ez)
    class(field_solver_type), intent(in) :: self
    real(dp), intent(in) :: r(3)
    real(dp), intent(inout) :: ex, ey, ez

    integer(i32) :: node_idx
    integer(i32) :: nimg, img_outer
    integer(i32) :: axis1, axis2
    real(dp) :: alpha

    if (.not. prepare_periodic2_ewald(self, alpha, nimg, img_outer, axis1, axis2)) return

    do node_idx = 1_i32, self%nnode
      if (self%child_count(node_idx) > 0_i32) cycle
      call add_ewald_node_correction(self, node_idx, r, axis1, axis2, nimg, img_outer, alpha, ex, ey, ez)
    end do
  end subroutine add_periodic2_ewald_like_correction_all_leaves

  logical function prepare_periodic2_ewald(self, alpha, nimg, img_outer, axis1, axis2)
    class(field_solver_type), intent(in) :: self
    real(dp), intent(out) :: alpha
    integer(i32), intent(out) :: nimg, img_outer, axis1, axis2

    prepare_periodic2_ewald = .false.
    alpha = 0.0d0
    nimg = 0_i32
    img_outer = 0_i32
    axis1 = 0_i32
    axis2 = 0_i32

    if (.not. self%use_periodic2) return
    if (self%periodic_ewald_layers <= 0_i32) return

    alpha = self%periodic_ewald_alpha
    if (alpha <= 0.0d0) return

    nimg = self%periodic_image_layers
    img_outer = nimg + self%periodic_ewald_layers
    if (img_outer <= nimg) return

    axis1 = self%periodic_axes(1)
    axis2 = self%periodic_axes(2)
    if (axis1 <= 0_i32 .or. axis2 <= 0_i32) return

    prepare_periodic2_ewald = .true.
  end function prepare_periodic2_ewald

  subroutine add_ewald_node_correction(self, node_idx, r, axis1, axis2, nimg, img_outer, alpha, ex, ey, ez)
    class(field_solver_type), intent(in) :: self
    integer(i32), intent(in) :: node_idx
    real(dp), intent(in) :: r(3)
    integer(i32), intent(in) :: axis1, axis2, nimg, img_outer
    real(dp), intent(in) :: alpha
    real(dp), intent(inout) :: ex, ey, ez

    real(dp) :: qi, src(3)

    qi = self%node_q(node_idx)
    if (abs(qi) <= tiny(1.0d0)) return
    src = self%node_charge_center(:, node_idx)
    call add_screened_shifted_node_images(qi, src, r, axis1, axis2, self%periodic_len, nimg, img_outer, alpha, ex, ey, ez)
  end subroutine add_ewald_node_correction

  !> 1 ノードの画像セル寄与を、指定画像範囲外（近傍外）に対して加算する。
  subroutine add_screened_shifted_node_images(q, src, target, axis1, axis2, periodic_len, near_img, outer_img, alpha, ex, ey, ez)
    real(dp), intent(in) :: q
    real(dp), intent(in) :: src(3)
    real(dp), intent(in) :: target(3)
    integer(i32), intent(in) :: axis1, axis2
    real(dp), intent(in) :: periodic_len(2)
    integer(i32), intent(in) :: near_img, outer_img
    real(dp), intent(in) :: alpha
    real(dp), intent(inout) :: ex, ey, ez

    integer(i32) :: img_i, img_j
    real(dp) :: shifted(3)

    if (abs(q) <= tiny(1.0d0)) return
    do img_i = -outer_img, outer_img
      do img_j = -outer_img, outer_img
        if (abs(img_i) <= near_img .and. abs(img_j) <= near_img) cycle
        shifted = src
        shifted(axis1) = shifted(axis1) + real(img_i, dp) * periodic_len(1)
        shifted(axis2) = shifted(axis2) + real(img_j, dp) * periodic_len(2)
        call add_screened_point_charge(q, target, shifted, alpha, ex, ey, ez)
      end do
    end do
  end subroutine add_screened_shifted_node_images

  !> erfc(alpha*r)/r 由来の勾配核で単一点電荷寄与を加算する。
  subroutine add_screened_point_charge(q, target, source, alpha, ex, ey, ez)
    real(dp), intent(in) :: q
    real(dp), intent(in) :: target(3), source(3)
    real(dp), intent(in) :: alpha
    real(dp), intent(inout) :: ex, ey, ez

    real(dp) :: dx, dy, dz, r2, rmag, inv_r, inv_r2, inv_r3
    real(dp) :: ar, screen, gaussian, pref

    if (abs(q) <= tiny(1.0d0)) return
    dx = target(1) - source(1)
    dy = target(2) - source(2)
    dz = target(3) - source(3)
    r2 = dx * dx + dy * dy + dz * dz
    if (r2 <= tiny(1.0d0)) return

    rmag = sqrt(r2)
    inv_r = 1.0d0 / rmag
    inv_r2 = inv_r * inv_r
    inv_r3 = inv_r2 * inv_r
    ar = alpha * rmag
    screen = erfc(ar)
    gaussian = exp(-(ar * ar))
    pref = q * (screen * inv_r3 + 2.0d0 * alpha * inv_sqrt_pi * gaussian * inv_r2)

    ex = ex + pref * dx
    ey = ey + pref * dy
    ez = ez + pref * dz
  end subroutine add_screened_point_charge

end submodule bem_field_solver_fmm

!> Coulomb FMM plan 構築と tree トポロジ前計算。
module bem_coulomb_fmm_plan_ops
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_types, only: fmm_options_type, fmm_plan_type, reset_fmm_plan
  use bem_coulomb_fmm_basis, only: initialize_basis_tables, build_axis_powers, compute_laplace_derivatives
  use bem_coulomb_fmm_periodic, only: has_valid_target_box, build_periodic_shift_values, distance_to_source_bbox, &
                                      distance_to_source_bbox_periodic
  use bem_coulomb_fmm_periodic_ewald, only: precompute_periodic2_ewald_data
  use bem_coulomb_fmm_periodic_root_ops, only: precompute_periodic_root_operator
  use bem_coulomb_fmm_tree_utils, only: octant_index, active_tree_nnode, active_tree_child_count, active_tree_child_idx, &
                                        active_tree_node_center, active_tree_node_radius, append_i32_buffer, &
                                        nodes_well_separated
  implicit none
  private

  public :: core_build_plan_impl
  public :: core_destroy_plan_impl

contains

  !> FMM 計画と木構造、転送演算子を構築する。
  !! @param[inout] plan 構築対象の FMM 計画。
  !! @param[in] src_pos ソース位置 `(3,n)` [m]。
  !! @param[in] options FMM 設定。
  subroutine core_build_plan_impl(plan, src_pos, options)
    type(fmm_plan_type), intent(inout) :: plan
    real(dp), intent(in) :: src_pos(:, :)
    type(fmm_options_type), intent(in) :: options
    integer(i32) :: nsrc

    if (allocated(plan%alpha) .or. allocated(plan%src_pos) .or. allocated(plan%elem_order)) then
      call core_destroy_plan_impl(plan)
    end if
    plan%options = options
    nsrc = int(size(src_pos, 2), i32)
    plan%nsrc = nsrc
    if (size(src_pos, 1) /= 3) error stop 'FMM core expects src_pos(3,n).'
    if (options%leaf_max <= 0_i32) error stop 'FMM leaf_max must be > 0.'
    if (options%order < 0_i32) error stop 'FMM order must be >= 0.'
    if (options%softening < 0.0d0) error stop 'FMM softening must be >= 0.'
    if (options%use_periodic2) then
      if (any(options%periodic_axes < 1_i32) .or. any(options%periodic_axes > 3_i32)) then
        error stop 'periodic2 requires periodic axes in 1..3.'
      end if
      if (options%periodic_axes(1) == options%periodic_axes(2)) then
        error stop 'periodic2 requires two distinct axes.'
      end if
      if (any(options%periodic_len <= 0.0d0)) then
        error stop 'periodic2 requires positive periodic lengths.'
      end if
      if (.not. has_valid_target_box(plan%options)) then
        error stop 'periodic2 requires a valid target box.'
      end if
      select case (trim(plan%options%periodic_far_correction))
      case ('auto', 'none')
        plan%options%periodic_far_correction = 'm2l_root_oracle'
        plan%options%periodic_ewald_layers = max(1_i32, plan%options%periodic_ewald_layers)
      case ('m2l_root_oracle')
        continue
      case default
        error stop 'Unsupported periodic far correction in FMM core.'
      end select
      if (trim(plan%options%periodic_far_correction) == 'm2l_root_oracle') then
        if (plan%options%periodic_ewald_layers < 1_i32) then
          error stop 'periodic2 root-operator far correction requires periodic_ewald_layers >= 1.'
        end if
      end if
    end if

    allocate (plan%src_pos(3, max(0_i32, nsrc)))
    if (nsrc > 0_i32) plan%src_pos = src_pos

    call initialize_basis_tables(plan, options%order)
    call build_source_tree(plan)
    call precompute_source_p2m_basis(plan)
    call build_target_topology(plan)
    call build_interactions(plan)
    call precompute_translation_operators(plan)
    call precompute_periodic2_ewald_data(plan)
    call precompute_periodic_root_operator(plan)
    call precompute_m2l_derivatives(plan)
    plan%built = .true.
  end subroutine core_build_plan_impl

  !> FMM 計画に確保した資源を解放する。
  !! @param[inout] plan 解放対象の FMM 計画。
  subroutine core_destroy_plan_impl(plan)
    type(fmm_plan_type), intent(inout) :: plan

    call reset_fmm_plan(plan)
  end subroutine core_destroy_plan_impl

  subroutine build_source_tree(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: i, max_node_guess

    if (plan%nsrc <= 0_i32) return
    max_node_guess = max(1_i32, 2_i32*plan%nsrc)
    call ensure_source_tree_capacity(plan, max_node_guess)
    allocate (plan%elem_order(plan%nsrc))
    do i = 1_i32, plan%nsrc
      plan%elem_order(i) = i
    end do
    plan%nnode = 1_i32
    plan%node_max_depth = 0_i32
    call build_source_node(plan, 1_i32, 1_i32, plan%nsrc, 0_i32)
    call rebuild_source_level_cache(plan)
    call build_source_parent_index(plan)
    call build_source_leaf_index(plan)
  end subroutine build_source_tree

  recursive subroutine build_source_node(plan, node_idx, start_idx, end_idx, depth)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32), intent(in) :: node_idx, start_idx, end_idx, depth
    integer(i32) :: count, p, idx, oct
    integer(i32) :: child_k, child_node, child_start, child_end
    integer(i32), allocatable :: counts(:), offsets(:), cursor(:), work(:)
    real(dp) :: bb_min(3), bb_max(3), span(3), center(3), split_eps

    count = end_idx - start_idx + 1_i32
    plan%node_start(node_idx) = start_idx
    plan%node_count(node_idx) = count
    plan%child_count(node_idx) = 0_i32
    plan%child_idx(:, node_idx) = 0_i32
    plan%child_octant(:, node_idx) = 0_i32
    plan%node_depth(node_idx) = depth
    plan%node_max_depth = max(plan%node_max_depth, depth)

    idx = plan%elem_order(start_idx)
    bb_min = plan%src_pos(:, idx)
    bb_max = bb_min
    do p = start_idx + 1_i32, end_idx
      idx = plan%elem_order(p)
      bb_min = min(bb_min, plan%src_pos(:, idx))
      bb_max = max(bb_max, plan%src_pos(:, idx))
    end do

    span = bb_max - bb_min
    center = 0.5d0*(bb_max + bb_min)
    plan%node_center(:, node_idx) = center
    plan%node_half_size(:, node_idx) = 0.5d0*span
    plan%node_radius(node_idx) = sqrt(sum(plan%node_half_size(:, node_idx)*plan%node_half_size(:, node_idx)))
    if (count <= plan%options%leaf_max) return

    split_eps = 1.0d-12*max(1.0d0, maxval(abs(center)))
    if (maxval(span) <= split_eps) return

    allocate (counts(8), offsets(8), cursor(8), work(count))
    counts = 0_i32
    do p = start_idx, end_idx
      idx = plan%elem_order(p)
      oct = octant_index(plan%src_pos(1, idx), plan%src_pos(2, idx), plan%src_pos(3, idx), center)
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
      idx = plan%elem_order(p)
      oct = octant_index(plan%src_pos(1, idx), plan%src_pos(2, idx), plan%src_pos(3, idx), center)
      work(cursor(oct)) = idx
      cursor(oct) = cursor(oct) + 1_i32
    end do
    plan%elem_order(start_idx:end_idx) = work

    child_k = 0_i32
    child_start = start_idx
    do oct = 1_i32, 8_i32
      if (counts(oct) <= 0_i32) cycle
      child_end = child_start + counts(oct) - 1_i32
      plan%nnode = plan%nnode + 1_i32
      if (plan%nnode > plan%max_node) error stop 'FMM source tree capacity exceeded.'
      child_node = plan%nnode
      child_k = child_k + 1_i32
      plan%child_idx(child_k, node_idx) = child_node
      plan%child_octant(child_k, node_idx) = oct
      call build_source_node(plan, child_node, child_start, child_end, depth + 1_i32)
      child_start = child_end + 1_i32
    end do
    plan%child_count(node_idx) = child_k
    deallocate (counts, offsets, cursor, work)
  end subroutine build_source_node

  subroutine ensure_source_tree_capacity(plan, max_node_needed)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32), intent(in) :: max_node_needed

    plan%max_node = max_node_needed
    allocate (plan%node_start(plan%max_node), plan%node_count(plan%max_node))
    allocate (plan%child_count(plan%max_node), plan%child_idx(8, plan%max_node), plan%child_octant(8, plan%max_node))
    allocate (plan%node_depth(plan%max_node))
    allocate (plan%node_center(3, plan%max_node), plan%node_half_size(3, plan%max_node), plan%node_radius(plan%max_node))
    plan%node_start = 0_i32
    plan%node_count = 0_i32
    plan%child_count = 0_i32
    plan%child_idx = 0_i32
    plan%child_octant = 0_i32
    plan%node_depth = 0_i32
    plan%node_center = 0.0d0
    plan%node_half_size = 0.0d0
    plan%node_radius = 0.0d0
  end subroutine ensure_source_tree_capacity

  subroutine rebuild_source_level_cache(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: node_idx, depth
    integer(i32), allocatable :: depth_count(:), cursor(:)

    if (plan%nnode <= 0_i32) return
    allocate (depth_count(plan%node_max_depth + 1_i32), cursor(plan%node_max_depth + 1_i32))
    depth_count = 0_i32
    do node_idx = 1_i32, plan%nnode
      depth_count(plan%node_depth(node_idx) + 1_i32) = depth_count(plan%node_depth(node_idx) + 1_i32) + 1_i32
    end do
    allocate (plan%node_level_start(plan%node_max_depth + 2_i32), plan%node_level_nodes(plan%nnode))
    plan%node_level_start(1) = 1_i32
    do depth = 0_i32, plan%node_max_depth
      plan%node_level_start(depth + 2_i32) = plan%node_level_start(depth + 1_i32) + depth_count(depth + 1_i32)
    end do
    cursor = plan%node_level_start(1:plan%node_max_depth + 1_i32)
    do node_idx = 1_i32, plan%nnode
      depth = plan%node_depth(node_idx)
      plan%node_level_nodes(cursor(depth + 1_i32)) = node_idx
      cursor(depth + 1_i32) = cursor(depth + 1_i32) + 1_i32
    end do
    deallocate (depth_count, cursor)
  end subroutine rebuild_source_level_cache

  subroutine build_source_parent_index(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: parent_node, child_k, child_node

    if (allocated(plan%source_parent_of)) deallocate (plan%source_parent_of)
    allocate (plan%source_parent_of(max(1_i32, plan%nnode)))
    plan%source_parent_of = 0_i32
    do parent_node = 1_i32, plan%nnode
      do child_k = 1_i32, plan%child_count(parent_node)
        child_node = plan%child_idx(child_k, parent_node)
        plan%source_parent_of(child_node) = parent_node
      end do
    end do
  end subroutine build_source_parent_index

  subroutine build_source_leaf_index(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: node_idx, leaf_idx

    if (allocated(plan%source_leaf_nodes)) deallocate (plan%source_leaf_nodes)
    plan%nsource_leaf = 0_i32
    do node_idx = 1_i32, plan%nnode
      if (plan%child_count(node_idx) <= 0_i32) plan%nsource_leaf = plan%nsource_leaf + 1_i32
    end do
    if (plan%nsource_leaf <= 0_i32) return

    allocate (plan%source_leaf_nodes(plan%nsource_leaf))
    leaf_idx = 0_i32
    do node_idx = 1_i32, plan%nnode
      if (plan%child_count(node_idx) > 0_i32) cycle
      leaf_idx = leaf_idx + 1_i32
      plan%source_leaf_nodes(leaf_idx) = node_idx
    end do
  end subroutine build_source_leaf_index

  subroutine precompute_source_p2m_basis(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: leaf_idx, node_idx, p, idx, p_end, alpha_idx
    real(dp) :: d(3)
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    if (allocated(plan%source_p2m_basis)) deallocate (plan%source_p2m_basis)
    if (plan%nsrc <= 0_i32) return

    allocate (plan%source_p2m_basis(plan%ncoef, plan%nsrc))
    plan%source_p2m_basis = 0.0d0
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(plan) private(leaf_idx, node_idx, p, idx, p_end, alpha_idx, d, xpow, ypow, zpow)
    do leaf_idx = 1_i32, plan%nsource_leaf
      node_idx = plan%source_leaf_nodes(leaf_idx)
      p_end = plan%node_start(node_idx) + plan%node_count(node_idx) - 1_i32
      do p = plan%node_start(node_idx), p_end
        idx = plan%elem_order(p)
        d = plan%src_pos(:, idx) - plan%node_center(:, node_idx)
        call build_axis_powers(d, plan%options%order, xpow, ypow, zpow)
        do alpha_idx = 1_i32, plan%ncoef
          plan%source_p2m_basis(alpha_idx, p) = xpow(plan%alpha(1, alpha_idx))*ypow(plan%alpha(2, alpha_idx)) &
                                                *zpow(plan%alpha(3, alpha_idx))/plan%alpha_factorial(alpha_idx)
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine precompute_source_p2m_basis

  subroutine build_target_topology(plan)
    type(fmm_plan_type), intent(inout) :: plan

    if (plan%nnode <= 0_i32) then
      call clear_target_tree(plan)
      return
    end if

    if (has_valid_target_box(plan%options)) then
      call build_box_target_tree(plan, plan%options%use_periodic2)
    else
      call assign_source_leaves_as_targets(plan)
    end if
  end subroutine build_target_topology

  subroutine clear_target_tree(plan)
    type(fmm_plan_type), intent(inout) :: plan

    if (allocated(plan%target_child_count)) deallocate (plan%target_child_count)
    if (allocated(plan%target_child_idx)) deallocate (plan%target_child_idx)
    if (allocated(plan%target_child_octant)) deallocate (plan%target_child_octant)
    if (allocated(plan%target_node_depth)) deallocate (plan%target_node_depth)
    if (allocated(plan%target_level_start)) deallocate (plan%target_level_start)
    if (allocated(plan%target_level_nodes)) deallocate (plan%target_level_nodes)
    if (allocated(plan%target_node_center)) deallocate (plan%target_node_center)
    if (allocated(plan%target_node_half_size)) deallocate (plan%target_node_half_size)
    if (allocated(plan%target_node_radius)) deallocate (plan%target_node_radius)
    if (allocated(plan%leaf_nodes)) deallocate (plan%leaf_nodes)
    if (allocated(plan%leaf_slot_of_node)) deallocate (plan%leaf_slot_of_node)
    plan%target_tree_ready = .false.
    plan%target_nnode = 0_i32
    plan%target_max_node = 0_i32
    plan%target_node_max_depth = 0_i32
    plan%nleaf = 0_i32
  end subroutine clear_target_tree

  subroutine assign_source_leaves_as_targets(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: node_idx, leaf_idx

    call clear_target_tree(plan)
    do node_idx = 1_i32, plan%nnode
      if (plan%child_count(node_idx) <= 0_i32) plan%nleaf = plan%nleaf + 1_i32
    end do
    if (plan%nleaf <= 0_i32) return
    allocate (plan%leaf_nodes(plan%nleaf), plan%leaf_slot_of_node(plan%max_node))
    plan%leaf_slot_of_node = 0_i32
    leaf_idx = 0_i32
    do node_idx = 1_i32, plan%nnode
      if (plan%child_count(node_idx) > 0_i32) cycle
      leaf_idx = leaf_idx + 1_i32
      plan%leaf_nodes(leaf_idx) = node_idx
      plan%leaf_slot_of_node(node_idx) = leaf_idx
    end do
  end subroutine assign_source_leaves_as_targets

  subroutine build_box_target_tree(plan, use_periodic_wrap)
    type(fmm_plan_type), intent(inout) :: plan
    logical, intent(in) :: use_periodic_wrap
    integer(i32) :: depth_min, depth_max
    integer(i32) :: leaf_idx, node_idx, target_cap, level_nodes, depth
    real(dp) :: root_min(3), root_max(3), src_center(3), src_half(3)

    if (plan%nnode <= 0_i32) then
      call clear_target_tree(plan)
      return
    end if

    src_center = plan%node_center(:, 1_i32)
    src_half = plan%node_half_size(:, 1_i32)
    root_min = plan%options%target_box_min
    root_max = plan%options%target_box_max
    if (any(root_max <= root_min)) then
      root_min = src_center - src_half
      root_max = src_center + src_half
    end if

    depth_min = 2_i32
    depth_max = 4_i32
    target_cap = 0_i32
    level_nodes = 1_i32
    do depth = 0_i32, depth_max
      target_cap = target_cap + level_nodes
      if (depth < depth_max) then
        if (level_nodes > huge(level_nodes)/8_i32) error stop 'FMM target tree capacity overflow.'
        level_nodes = 8_i32*level_nodes
      end if
    end do

    call clear_target_tree(plan)
    plan%target_max_node = max(1024_i32, max(16_i32*plan%nnode, target_cap))
    allocate (plan%target_child_count(plan%target_max_node), plan%target_child_idx(8, plan%target_max_node))
    allocate (plan%target_child_octant(8, plan%target_max_node))
    allocate (plan%target_node_depth(plan%target_max_node))
    allocate (plan%target_node_center(3, plan%target_max_node), plan%target_node_half_size(3, plan%target_max_node))
    allocate (plan%target_node_radius(plan%target_max_node))
    plan%target_child_count = 0_i32
    plan%target_child_idx = 0_i32
    plan%target_child_octant = 0_i32
    plan%target_node_depth = 0_i32
    plan%target_node_center = 0.0d0
    plan%target_node_half_size = 0.0d0
    plan%target_node_radius = 0.0d0

    plan%target_nnode = 1_i32
    plan%target_node_depth(1_i32) = 0_i32
    plan%target_node_center(:, 1_i32) = 0.5d0*(root_min + root_max)
    plan%target_node_half_size(:, 1_i32) = 0.5d0*(root_max - root_min)
    plan%target_node_radius(1_i32) = sqrt(sum(plan%target_node_half_size(:, 1_i32)*plan%target_node_half_size(:, 1_i32)))
    call build_target_box_node(plan, 1_i32, 0_i32, depth_min, depth_max, src_center, src_half, use_periodic_wrap)

    do node_idx = 1_i32, plan%target_nnode
      if (plan%target_child_count(node_idx) <= 0_i32) plan%nleaf = plan%nleaf + 1_i32
    end do
    if (plan%nleaf <= 0_i32) return

    allocate (plan%leaf_nodes(plan%nleaf), plan%leaf_slot_of_node(plan%target_max_node))
    plan%leaf_slot_of_node = 0_i32
    leaf_idx = 0_i32
    do node_idx = 1_i32, plan%target_nnode
      if (plan%target_child_count(node_idx) > 0_i32) cycle
      leaf_idx = leaf_idx + 1_i32
      plan%leaf_nodes(leaf_idx) = node_idx
      plan%leaf_slot_of_node(node_idx) = leaf_idx
    end do

    call rebuild_target_level_cache(plan)
    plan%target_tree_ready = .true.
  end subroutine build_box_target_tree

  recursive subroutine build_target_box_node(plan, node_idx, depth, depth_min, depth_max, src_center, src_half, use_periodic_wrap)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32), intent(in) :: node_idx, depth, depth_min, depth_max
    real(dp), intent(in) :: src_center(3), src_half(3)
    logical, intent(in) :: use_periodic_wrap
    integer(i32) :: child_k, child_node, oct
    real(dp) :: child_half(3), offset(3)

    if (.not. target_should_split_box( &
        plan, node_idx, depth, depth_min, depth_max, src_center, src_half, use_periodic_wrap)) return
    child_half = 0.5d0*plan%target_node_half_size(:, node_idx)
    if (maxval(child_half) <= tiny(1.0d0)) return

    child_k = 0_i32
    do oct = 1_i32, 8_i32
      plan%target_nnode = plan%target_nnode + 1_i32
      if (plan%target_nnode > plan%target_max_node) error stop 'FMM target tree capacity exceeded.'
      child_node = plan%target_nnode
      child_k = child_k + 1_i32
      plan%target_child_idx(child_k, node_idx) = child_node
      plan%target_child_octant(child_k, node_idx) = oct
      plan%target_node_depth(child_node) = depth + 1_i32
      plan%target_node_max_depth = max(plan%target_node_max_depth, depth + 1_i32)

      offset(1) = merge(-1.0d0, 1.0d0, iand(oct - 1_i32, 1_i32) == 0_i32)
      offset(2) = merge(-1.0d0, 1.0d0, iand(oct - 1_i32, 2_i32) == 0_i32)
      offset(3) = merge(-1.0d0, 1.0d0, iand(oct - 1_i32, 4_i32) == 0_i32)
      plan%target_node_half_size(:, child_node) = child_half
      plan%target_node_center(:, child_node) = plan%target_node_center(:, node_idx) + offset*child_half
      plan%target_node_radius(child_node) = sqrt(sum(child_half*child_half))
      call build_target_box_node( &
        plan, child_node, depth + 1_i32, depth_min, depth_max, src_center, src_half, use_periodic_wrap &
        )
    end do
    plan%target_child_count(node_idx) = child_k
  end subroutine build_target_box_node

  logical function target_should_split_box(plan, node_idx, depth, depth_min, depth_max, src_center, src_half, use_periodic_wrap)
    type(fmm_plan_type), intent(in) :: plan
    integer(i32), intent(in) :: node_idx, depth, depth_min, depth_max
    real(dp), intent(in) :: src_center(3), src_half(3)
    logical, intent(in) :: use_periodic_wrap
    real(dp) :: dsrc, radius

    if (depth < depth_min) then
      target_should_split_box = .true.
      return
    end if
    if (depth >= depth_max) then
      target_should_split_box = .false.
      return
    end if

    radius = plan%target_node_radius(node_idx)
    if (use_periodic_wrap) then
      dsrc = distance_to_source_bbox_periodic(plan, plan%target_node_center(:, node_idx), src_center, src_half)
    else
      dsrc = distance_to_source_bbox(plan%target_node_center(:, node_idx), src_center, src_half)
    end if
    target_should_split_box = (radius > 0.35d0*max(dsrc, 1.0d-12))
  end function target_should_split_box

  subroutine rebuild_target_level_cache(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: node_idx, depth
    integer(i32), allocatable :: depth_count(:), cursor(:)

    if (plan%target_nnode <= 0_i32) return
    allocate (depth_count(plan%target_node_max_depth + 1_i32), cursor(plan%target_node_max_depth + 1_i32))
    depth_count = 0_i32
    do node_idx = 1_i32, plan%target_nnode
      depth_count(plan%target_node_depth(node_idx) + 1_i32) = depth_count(plan%target_node_depth(node_idx) + 1_i32) + 1_i32
    end do
    allocate (plan%target_level_start(plan%target_node_max_depth + 2_i32), plan%target_level_nodes(plan%target_nnode))
    plan%target_level_start(1) = 1_i32
    do depth = 0_i32, plan%target_node_max_depth
      plan%target_level_start(depth + 2_i32) = plan%target_level_start(depth + 1_i32) + depth_count(depth + 1_i32)
    end do
    cursor = plan%target_level_start(1:plan%target_node_max_depth + 1_i32)
    do node_idx = 1_i32, plan%target_nnode
      depth = plan%target_node_depth(node_idx)
      plan%target_level_nodes(cursor(depth + 1_i32)) = node_idx
      cursor(depth + 1_i32) = cursor(depth + 1_i32) + 1_i32
    end do
    deallocate (depth_count, cursor)
  end subroutine rebuild_target_level_cache

  subroutine build_interactions(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: leaf_idx, i, near_n, far_n, near_used, far_used, near_cap, far_cap
    integer(i32) :: n_target_nodes, child_k, parent_node, node_idx
    integer(i32) :: pair_cap, pair_used, nshift
    integer(i32), allocatable :: near_work(:), far_work(:), near_buf(:), far_buf(:)
    integer(i32), allocatable :: m2l_target_buf(:), m2l_source_buf(:), m2l_shift_idx1_buf(:), m2l_shift_idx2_buf(:)
    logical :: use_target_tree

    nshift = 1_i32
    if (plan%options%use_periodic2) nshift = 2_i32*max(0_i32, plan%options%periodic_image_layers) + 1_i32
    allocate (plan%shift_axis1(nshift), plan%shift_axis2(nshift))
    call build_periodic_shift_values(plan, plan%shift_axis1, plan%shift_axis2, nshift)

    if (plan%nleaf <= 0_i32) return

    use_target_tree = plan%target_tree_ready
    n_target_nodes = active_tree_nnode(plan, use_target_tree)
    allocate (plan%near_start(plan%nleaf + 1_i32), plan%far_start(plan%nleaf + 1_i32))
    allocate (near_work(plan%nnode), far_work(plan%nnode))
    near_cap = max(64_i32, plan%nleaf*16_i32)
    far_cap = max(64_i32, plan%nleaf*16_i32)
    allocate (near_buf(near_cap), far_buf(far_cap))
    near_used = 0_i32
    far_used = 0_i32
    do leaf_idx = 1_i32, plan%nleaf
      near_n = 0_i32
      far_n = 0_i32
      call gather_leaf_interactions(plan, plan%leaf_nodes(leaf_idx), 1_i32, near_work, near_n, far_work, far_n)
      plan%near_start(leaf_idx) = near_used + 1_i32
      plan%far_start(leaf_idx) = far_used + 1_i32
      do i = 1_i32, near_n
        call append_i32_buffer(near_buf, near_used, near_cap, near_work(i))
      end do
      do i = 1_i32, far_n
        call append_i32_buffer(far_buf, far_used, far_cap, far_work(i))
      end do
    end do
    plan%near_start(plan%nleaf + 1_i32) = near_used + 1_i32
    plan%far_start(plan%nleaf + 1_i32) = far_used + 1_i32
    allocate (plan%near_nodes(max(1_i32, near_used)), plan%far_nodes(max(1_i32, far_used)))
    plan%near_nodes = 0_i32
    plan%far_nodes = 0_i32
    if (near_used > 0_i32) plan%near_nodes(1:near_used) = near_buf(1:near_used)
    if (far_used > 0_i32) plan%far_nodes(1:far_used) = far_buf(1:far_used)

    allocate (plan%parent_of(n_target_nodes))
    plan%parent_of = 0_i32
    do parent_node = 1_i32, n_target_nodes
      do child_k = 1_i32, active_tree_child_count(plan, use_target_tree, parent_node)
        node_idx = active_tree_child_idx(plan, use_target_tree, child_k, parent_node)
        plan%parent_of(node_idx) = parent_node
      end do
    end do

    pair_cap = max(64_i32, plan%nleaf*32_i32)
    allocate (m2l_target_buf(pair_cap), m2l_source_buf(pair_cap), m2l_shift_idx1_buf(pair_cap), m2l_shift_idx2_buf(pair_cap))
    m2l_target_buf = 0_i32
    m2l_source_buf = 0_i32
    m2l_shift_idx1_buf = 0_i32
    m2l_shift_idx2_buf = 0_i32
    pair_used = 0_i32
    plan%m2l_visit_count = 0_i32
    call accumulate_cached_pairs( &
      plan, 1_i32, 1_i32, use_target_tree, m2l_target_buf, m2l_source_buf, m2l_shift_idx1_buf, m2l_shift_idx2_buf, &
      pair_used, pair_cap &
      )
    call build_near_source_index( &
      plan, m2l_target_buf, m2l_source_buf, m2l_shift_idx1_buf, m2l_shift_idx2_buf, pair_used, pair_cap &
      )
    allocate ( &
      plan%m2l_target_nodes(max(1_i32, pair_used)), plan%m2l_source_nodes(max(1_i32, pair_used)), &
      plan%m2l_shift_idx1(max(1_i32, pair_used)), plan%m2l_shift_idx2(max(1_i32, pair_used)) &
      )
    plan%m2l_target_nodes = 0_i32
    plan%m2l_source_nodes = 0_i32
    plan%m2l_shift_idx1 = 0_i32
    plan%m2l_shift_idx2 = 0_i32
    if (pair_used > 0_i32) then
      plan%m2l_target_nodes(1:pair_used) = m2l_target_buf(1:pair_used)
      plan%m2l_source_nodes(1:pair_used) = m2l_source_buf(1:pair_used)
      plan%m2l_shift_idx1(1:pair_used) = m2l_shift_idx1_buf(1:pair_used)
      plan%m2l_shift_idx2(1:pair_used) = m2l_shift_idx2_buf(1:pair_used)
    end if
    plan%m2l_pair_count = pair_used
    plan%m2l_build_count = plan%m2l_build_count + 1_i32
    call build_target_pair_index(plan, n_target_nodes)

    deallocate (near_work, far_work, near_buf, far_buf, m2l_target_buf, m2l_source_buf, m2l_shift_idx1_buf, m2l_shift_idx2_buf)
  end subroutine build_interactions

  subroutine build_near_source_index( &
    plan, m2l_target_buf, m2l_source_buf, m2l_shift_idx1_buf, m2l_shift_idx2_buf, pair_used, pair_cap &
    )
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32), allocatable, intent(inout) :: m2l_target_buf(:), m2l_source_buf(:), &
                                                m2l_shift_idx1_buf(:), m2l_shift_idx2_buf(:)
    integer(i32), intent(inout) :: pair_used, pair_cap
    integer(i32) :: leaf_idx, near_pos, near_node, p, p_end, near_source_used
    integer(i32) :: shift_idx1, shift_idx2, nshift, target_leaf
    real(dp) :: shift1, shift2

    if (allocated(plan%near_source_start)) deallocate (plan%near_source_start)
    if (allocated(plan%near_source_idx)) deallocate (plan%near_source_idx)
    if (allocated(plan%near_source_shift1)) deallocate (plan%near_source_shift1)
    if (allocated(plan%near_source_shift2)) deallocate (plan%near_source_shift2)
    if (plan%nleaf <= 0_i32) return

    allocate (plan%near_source_start(plan%nleaf + 1_i32))
    plan%near_source_start = 1_i32
    near_source_used = 0_i32
    nshift = max(1_i32, int(size(plan%shift_axis1), i32))
    do leaf_idx = 1_i32, plan%nleaf
      plan%near_source_start(leaf_idx) = near_source_used + 1_i32
      target_leaf = plan%leaf_nodes(leaf_idx)
      do near_pos = plan%near_start(leaf_idx), plan%near_start(leaf_idx + 1_i32) - 1_i32
        near_node = plan%near_nodes(near_pos)
        if (plan%options%use_periodic2) then
          do shift_idx1 = 1_i32, nshift
            shift1 = plan%shift_axis1(shift_idx1)
            do shift_idx2 = 1_i32, nshift
              shift2 = plan%shift_axis2(shift_idx2)
              if (.not. nodes_well_separated_shifted(plan, target_leaf, near_node, shift1, shift2)) then
                near_source_used = near_source_used + plan%node_count(near_node)
              end if
            end do
          end do
        else
          near_source_used = near_source_used + plan%node_count(near_node)
        end if
      end do
    end do
    plan%near_source_start(plan%nleaf + 1_i32) = near_source_used + 1_i32

    allocate ( &
      plan%near_source_idx(max(1_i32, near_source_used)), &
      plan%near_source_shift1(max(1_i32, near_source_used)), &
      plan%near_source_shift2(max(1_i32, near_source_used)) &
      )
    plan%near_source_idx = 0_i32
    plan%near_source_shift1 = 0.0d0
    plan%near_source_shift2 = 0.0d0
    if (near_source_used <= 0_i32) return

    near_source_used = 0_i32
    do leaf_idx = 1_i32, plan%nleaf
      target_leaf = plan%leaf_nodes(leaf_idx)
      do near_pos = plan%near_start(leaf_idx), plan%near_start(leaf_idx + 1_i32) - 1_i32
        near_node = plan%near_nodes(near_pos)
        p_end = plan%node_start(near_node) + plan%node_count(near_node) - 1_i32
        if (plan%options%use_periodic2) then
          do shift_idx1 = 1_i32, nshift
            shift1 = plan%shift_axis1(shift_idx1)
            do shift_idx2 = 1_i32, nshift
              shift2 = plan%shift_axis2(shift_idx2)
              if (nodes_well_separated_shifted(plan, target_leaf, near_node, shift1, shift2)) then
                call append_m2l_pair( &
                  m2l_target_buf, m2l_source_buf, m2l_shift_idx1_buf, m2l_shift_idx2_buf, pair_used, pair_cap, &
                  target_leaf, near_node, shift_idx1, shift_idx2 &
                  )
                cycle
              end if
              do p = plan%node_start(near_node), p_end
                near_source_used = near_source_used + 1_i32
                plan%near_source_idx(near_source_used) = plan%elem_order(p)
                plan%near_source_shift1(near_source_used) = shift1
                plan%near_source_shift2(near_source_used) = shift2
              end do
            end do
          end do
        else
          do p = plan%node_start(near_node), p_end
            near_source_used = near_source_used + 1_i32
            plan%near_source_idx(near_source_used) = plan%elem_order(p)
          end do
        end if
      end do
    end do
  end subroutine build_near_source_index

  recursive subroutine gather_leaf_interactions(plan, target_leaf, source_node, near_buf, near_n, far_buf, far_n)
    type(fmm_plan_type), intent(in) :: plan
    integer(i32), intent(in) :: target_leaf, source_node
    integer(i32), intent(inout) :: near_buf(:), near_n, far_buf(:), far_n
    integer(i32) :: child_k, child_node

    if (nodes_well_separated(plan, target_leaf, source_node)) then
      far_n = far_n + 1_i32
      far_buf(far_n) = source_node
      return
    end if
    if (plan%child_count(source_node) <= 0_i32) then
      near_n = near_n + 1_i32
      near_buf(near_n) = source_node
      return
    end if
    do child_k = 1_i32, plan%child_count(source_node)
      child_node = plan%child_idx(child_k, source_node)
      call gather_leaf_interactions(plan, target_leaf, child_node, near_buf, near_n, far_buf, far_n)
    end do
  end subroutine gather_leaf_interactions

  recursive subroutine accumulate_cached_pairs( &
    plan, t_node, s_node, use_target_tree, target_buf, source_buf, shift_idx1_buf, shift_idx2_buf, pair_used, pair_cap &
    )
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32), intent(in) :: t_node, s_node
    logical, intent(in) :: use_target_tree
    integer(i32), allocatable, intent(inout) :: target_buf(:), source_buf(:), shift_idx1_buf(:), shift_idx2_buf(:)
    integer(i32), intent(inout) :: pair_used, pair_cap
    integer(i32) :: tc, sc, t_child, s_child

    plan%m2l_visit_count = plan%m2l_visit_count + 1_i32
    tc = active_tree_child_count(plan, use_target_tree, t_node)
    sc = plan%child_count(s_node)
    if (nodes_well_separated(plan, t_node, s_node)) then
      call append_m2l_pair( &
        target_buf, source_buf, shift_idx1_buf, shift_idx2_buf, pair_used, pair_cap, t_node, s_node, 0_i32, 0_i32 &
        )
      return
    end if
    if (tc <= 0_i32 .and. sc <= 0_i32) return
    if (sc > 0_i32 .and. (tc <= 0_i32 .or. plan%node_radius(s_node) >= &
                          active_tree_node_radius(plan, use_target_tree, t_node))) then
      do s_child = 1_i32, sc
        call accumulate_cached_pairs( &
          plan, t_node, plan%child_idx(s_child, s_node), use_target_tree, target_buf, source_buf, shift_idx1_buf, &
          shift_idx2_buf, pair_used, pair_cap &
          )
      end do
    else
      do t_child = 1_i32, tc
        call accumulate_cached_pairs( &
          plan, active_tree_child_idx(plan, use_target_tree, t_child, t_node), s_node, use_target_tree, &
          target_buf, source_buf, shift_idx1_buf, shift_idx2_buf, pair_used, pair_cap &
          )
      end do
    end if
  end subroutine accumulate_cached_pairs

  subroutine append_m2l_pair( &
    target_buf, source_buf, shift_idx1_buf, shift_idx2_buf, pair_used, pair_cap, t_node, s_node, shift_idx1, shift_idx2 &
    )
    integer(i32), allocatable, intent(inout) :: target_buf(:), source_buf(:), shift_idx1_buf(:), shift_idx2_buf(:)
    integer(i32), intent(inout) :: pair_used, pair_cap
    integer(i32), intent(in) :: t_node, s_node, shift_idx1, shift_idx2
    integer(i32), allocatable :: tmp_target(:), tmp_source(:), tmp_shift_idx1(:), tmp_shift_idx2(:)
    integer(i32) :: new_capacity

    if (pair_used >= pair_cap) then
      new_capacity = max(pair_cap*2_i32, pair_cap + 64_i32)
      allocate (tmp_target(new_capacity), tmp_source(new_capacity), &
                tmp_shift_idx1(new_capacity), tmp_shift_idx2(new_capacity))
      tmp_target = 0_i32
      tmp_source = 0_i32
      tmp_shift_idx1 = 0_i32
      tmp_shift_idx2 = 0_i32
      if (pair_used > 0_i32) then
        tmp_target(1:pair_used) = target_buf(1:pair_used)
        tmp_source(1:pair_used) = source_buf(1:pair_used)
        tmp_shift_idx1(1:pair_used) = shift_idx1_buf(1:pair_used)
        tmp_shift_idx2(1:pair_used) = shift_idx2_buf(1:pair_used)
      end if
      call move_alloc(tmp_target, target_buf)
      call move_alloc(tmp_source, source_buf)
      call move_alloc(tmp_shift_idx1, shift_idx1_buf)
      call move_alloc(tmp_shift_idx2, shift_idx2_buf)
      pair_cap = new_capacity
    end if
    pair_used = pair_used + 1_i32
    target_buf(pair_used) = t_node
    source_buf(pair_used) = s_node
    shift_idx1_buf(pair_used) = shift_idx1
    shift_idx2_buf(pair_used) = shift_idx2
  end subroutine append_m2l_pair

  subroutine build_target_pair_index(plan, n_target_nodes)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32), intent(in) :: n_target_nodes
    integer(i32) :: pair_idx, target_node, pos
    integer(i32), allocatable :: counts(:), cursor(:)

    allocate (plan%m2l_target_start(n_target_nodes + 1_i32))
    allocate (plan%m2l_pair_order(max(1_i32, plan%m2l_pair_count)))
    plan%m2l_target_start = 1_i32
    plan%m2l_pair_order = 0_i32
    if (plan%m2l_pair_count <= 0_i32) return

    allocate (counts(n_target_nodes), cursor(n_target_nodes))
    counts = 0_i32
    do pair_idx = 1_i32, plan%m2l_pair_count
      target_node = plan%m2l_target_nodes(pair_idx)
      counts(target_node) = counts(target_node) + 1_i32
    end do
    plan%m2l_target_start(1) = 1_i32
    do target_node = 1_i32, n_target_nodes
      plan%m2l_target_start(target_node + 1_i32) = plan%m2l_target_start(target_node) + counts(target_node)
    end do
    cursor = plan%m2l_target_start(1:n_target_nodes)
    do pair_idx = 1_i32, plan%m2l_pair_count
      target_node = plan%m2l_target_nodes(pair_idx)
      pos = cursor(target_node)
      plan%m2l_pair_order(pos) = pair_idx
      cursor(target_node) = pos + 1_i32
    end do
    deallocate (counts, cursor)
  end subroutine build_target_pair_index

  subroutine precompute_translation_operators(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: alpha_idx, beta_idx, gamma_idx, delta_idx, term_idx
    integer(i32) :: node_idx, n_target_nodes, parent_node
    integer(i32) :: delta(3)
    integer(i32) :: max_m2m_terms, max_l2l_terms
    real(dp) :: d(3)
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    if (allocated(plan%m2m_term_count)) deallocate (plan%m2m_term_count)
    if (allocated(plan%m2m_alpha_list)) deallocate (plan%m2m_alpha_list)
    if (allocated(plan%m2m_delta_list)) deallocate (plan%m2m_delta_list)
    if (allocated(plan%l2l_term_count)) deallocate (plan%l2l_term_count)
    if (allocated(plan%l2l_gamma_list)) deallocate (plan%l2l_gamma_list)
    if (allocated(plan%l2l_delta_list)) deallocate (plan%l2l_delta_list)

    allocate (plan%m2m_term_count(plan%ncoef), plan%l2l_term_count(plan%ncoef))
    plan%m2m_term_count = 0_i32
    plan%l2l_term_count = 0_i32
    do beta_idx = 1_i32, plan%ncoef
      do alpha_idx = 1_i32, plan%ncoef
        if (any(plan%alpha(:, alpha_idx) > plan%alpha(:, beta_idx))) cycle
        plan%m2m_term_count(beta_idx) = plan%m2m_term_count(beta_idx) + 1_i32
      end do
    end do
    do alpha_idx = 1_i32, plan%ncoef
      do gamma_idx = 1_i32, plan%ncoef
        if (any(plan%alpha(:, gamma_idx) < plan%alpha(:, alpha_idx))) cycle
        plan%l2l_term_count(alpha_idx) = plan%l2l_term_count(alpha_idx) + 1_i32
      end do
    end do

    max_m2m_terms = max(1_i32, maxval(plan%m2m_term_count))
    max_l2l_terms = max(1_i32, maxval(plan%l2l_term_count))
    allocate (plan%m2m_alpha_list(max_m2m_terms, plan%ncoef), plan%m2m_delta_list(max_m2m_terms, plan%ncoef))
    allocate (plan%l2l_gamma_list(max_l2l_terms, plan%ncoef), plan%l2l_delta_list(max_l2l_terms, plan%ncoef))
    plan%m2m_alpha_list = 0_i32
    plan%m2m_delta_list = 0_i32
    plan%l2l_gamma_list = 0_i32
    plan%l2l_delta_list = 0_i32

    do beta_idx = 1_i32, plan%ncoef
      term_idx = 0_i32
      do alpha_idx = 1_i32, plan%ncoef
        if (any(plan%alpha(:, alpha_idx) > plan%alpha(:, beta_idx))) cycle
        delta = plan%alpha(:, beta_idx) - plan%alpha(:, alpha_idx)
        term_idx = term_idx + 1_i32
        plan%m2m_alpha_list(term_idx, beta_idx) = alpha_idx
        plan%m2m_delta_list(term_idx, beta_idx) = plan%alpha_map(delta(1), delta(2), delta(3))
      end do
    end do
    do alpha_idx = 1_i32, plan%ncoef
      term_idx = 0_i32
      do gamma_idx = 1_i32, plan%ncoef
        if (any(plan%alpha(:, gamma_idx) < plan%alpha(:, alpha_idx))) cycle
        delta = plan%alpha(:, gamma_idx) - plan%alpha(:, alpha_idx)
        term_idx = term_idx + 1_i32
        plan%l2l_gamma_list(term_idx, alpha_idx) = gamma_idx
        plan%l2l_delta_list(term_idx, alpha_idx) = plan%alpha_map(delta(1), delta(2), delta(3))
      end do
    end do

    if (allocated(plan%source_shift_monomial)) deallocate (plan%source_shift_monomial)
    allocate (plan%source_shift_monomial(plan%ncoef, max(1_i32, plan%nnode)))
    plan%source_shift_monomial = 0.0d0
    if (plan%nnode > 0_i32) plan%source_shift_monomial(1_i32, 1_i32) = 1.0d0
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(plan) private(node_idx, parent_node, d, xpow, ypow, zpow, delta_idx)
    do node_idx = 1_i32, plan%nnode
      parent_node = plan%source_parent_of(node_idx)
      if (parent_node <= 0_i32) cycle
      d = plan%node_center(:, node_idx) - plan%node_center(:, parent_node)
      call build_axis_powers(d, plan%options%order, xpow, ypow, zpow)
      do delta_idx = 1_i32, plan%ncoef
        plan%source_shift_monomial(delta_idx, node_idx) = xpow(plan%alpha(1, delta_idx))*ypow(plan%alpha(2, delta_idx)) &
                                                          *zpow(plan%alpha(3, delta_idx))/plan%alpha_factorial(delta_idx)
      end do
    end do
    !$omp end parallel do

    n_target_nodes = active_tree_nnode(plan, plan%target_tree_ready)
    if (allocated(plan%target_shift_monomial)) deallocate (plan%target_shift_monomial)
    allocate (plan%target_shift_monomial(plan%ncoef, max(1_i32, n_target_nodes)))
    plan%target_shift_monomial = 0.0d0
    if (n_target_nodes > 0_i32) plan%target_shift_monomial(1_i32, 1_i32) = 1.0d0
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(plan, n_target_nodes) private(node_idx, parent_node, d, xpow, ypow, zpow, delta_idx)
    do node_idx = 1_i32, n_target_nodes
      parent_node = plan%parent_of(node_idx)
      if (parent_node <= 0_i32) cycle
      d = active_tree_node_center(plan, plan%target_tree_ready, node_idx) - &
          active_tree_node_center(plan, plan%target_tree_ready, parent_node)
      call build_axis_powers(d, plan%options%order, xpow, ypow, zpow)
      do delta_idx = 1_i32, plan%ncoef
        plan%target_shift_monomial(delta_idx, node_idx) = xpow(plan%alpha(1, delta_idx))*ypow(plan%alpha(2, delta_idx)) &
                                                          *zpow(plan%alpha(3, delta_idx))/plan%alpha_factorial(delta_idx)
      end do
    end do
    !$omp end parallel do
  end subroutine precompute_translation_operators

  subroutine precompute_m2l_derivatives(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: pair_idx, img_i, img_j, nshift
    integer(i32) :: axis1, axis2, s_node, t_node, shift_idx1, shift_idx2
    real(dp) :: base(3), shift1, shift2, rimg(3)
    real(dp) :: deriv_tmp(plan%nderiv)

    allocate (plan%m2l_deriv(plan%nderiv, max(1_i32, plan%m2l_pair_count)))
    plan%m2l_deriv = 0.0d0
    if (plan%m2l_pair_count <= 0_i32) return

    axis1 = 0_i32
    axis2 = 0_i32
    if (plan%options%use_periodic2) then
      axis1 = plan%options%periodic_axes(1)
      axis2 = plan%options%periodic_axes(2)
    end if
    nshift = size(plan%shift_axis1)
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(plan, axis1, axis2, nshift) &
    !$omp private(pair_idx, img_i, img_j, s_node, t_node, base, shift1, shift2, shift_idx1, shift_idx2) &
    !$omp private(rimg, deriv_tmp)
    do pair_idx = 1_i32, plan%m2l_pair_count
      t_node = plan%m2l_target_nodes(pair_idx)
      s_node = plan%m2l_source_nodes(pair_idx)
      base = active_tree_node_center(plan, plan%target_tree_ready, t_node) - plan%node_center(:, s_node)
      if (.not. plan%options%use_periodic2) then
        call compute_laplace_derivatives(plan, base, deriv_tmp)
        plan%m2l_deriv(:, pair_idx) = plan%m2l_deriv(:, pair_idx) + deriv_tmp
        cycle
      end if
      shift_idx1 = plan%m2l_shift_idx1(pair_idx)
      shift_idx2 = plan%m2l_shift_idx2(pair_idx)
      if (shift_idx1 > 0_i32 .and. shift_idx2 > 0_i32) then
        shift1 = plan%shift_axis1(shift_idx1)
        shift2 = plan%shift_axis2(shift_idx2)
        rimg = base
        if (axis1 == 1_i32) rimg(1) = rimg(1) - shift1
        if (axis1 == 2_i32) rimg(2) = rimg(2) - shift1
        if (axis1 == 3_i32) rimg(3) = rimg(3) - shift1
        if (axis2 == 1_i32) rimg(1) = rimg(1) - shift2
        if (axis2 == 2_i32) rimg(2) = rimg(2) - shift2
        if (axis2 == 3_i32) rimg(3) = rimg(3) - shift2
        call compute_laplace_derivatives(plan, rimg, deriv_tmp)
        plan%m2l_deriv(:, pair_idx) = plan%m2l_deriv(:, pair_idx) + deriv_tmp
        cycle
      end if
      do img_i = 1_i32, nshift
        shift1 = plan%shift_axis1(img_i)
        do img_j = 1_i32, nshift
          shift2 = plan%shift_axis2(img_j)
          rimg = base
          if (axis1 == 1_i32) rimg(1) = rimg(1) - shift1
          if (axis1 == 2_i32) rimg(2) = rimg(2) - shift1
          if (axis1 == 3_i32) rimg(3) = rimg(3) - shift1
          if (axis2 == 1_i32) rimg(1) = rimg(1) - shift2
          if (axis2 == 2_i32) rimg(2) = rimg(2) - shift2
          if (axis2 == 3_i32) rimg(3) = rimg(3) - shift2
          call compute_laplace_derivatives(plan, rimg, deriv_tmp)
          plan%m2l_deriv(:, pair_idx) = plan%m2l_deriv(:, pair_idx) + deriv_tmp
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine precompute_m2l_derivatives

  logical function nodes_well_separated_shifted(plan, target_node, source_node, shift1, shift2)
    type(fmm_plan_type), intent(in) :: plan
    integer(i32), intent(in) :: target_node, source_node
    real(dp), intent(in) :: shift1, shift2
    real(dp) :: d(3), rs, rt, lhs, rhs, theta_eff, target_center(3), shifted_center(3)
    integer(i32) :: axis1, axis2
    logical :: use_target_tree

    use_target_tree = plan%target_tree_ready
    target_center = active_tree_node_center(plan, use_target_tree, target_node)
    rt = active_tree_node_radius(plan, use_target_tree, target_node)

    shifted_center = plan%node_center(:, source_node)
    axis1 = 0_i32
    axis2 = 0_i32
    if (plan%options%use_periodic2) then
      axis1 = plan%options%periodic_axes(1)
      axis2 = plan%options%periodic_axes(2)
    end if
    if (axis1 > 0_i32) shifted_center(axis1) = shifted_center(axis1) + shift1
    if (axis2 > 0_i32) shifted_center(axis2) = shifted_center(axis2) + shift2

    d = target_center - shifted_center
    if (sum(d*d) <= 0.0d0) then
      nodes_well_separated_shifted = .false.
      return
    end if

    rs = plan%node_radius(source_node)
    theta_eff = plan%options%theta
    lhs = (rs + rt)*(rs + rt)
    rhs = (theta_eff*theta_eff)*sum(d*d)
    nodes_well_separated_shifted = (lhs < rhs)
  end function nodes_well_separated_shifted

end module bem_coulomb_fmm_plan_ops

!> Coulomb FMM tree 構造の共通ユーティリティ。
module bem_coulomb_fmm_tree_utils
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_types, only: fmm_plan_type
  use bem_coulomb_fmm_periodic, only: apply_periodic2_minimum_image
  implicit none
  private

  public :: octant_index
  public :: active_tree_nnode
  public :: active_tree_child_count
  public :: active_tree_child_idx
  public :: active_tree_child_octant
  public :: active_tree_node_center
  public :: active_tree_node_half_size
  public :: active_tree_node_radius
  public :: active_tree_max_depth
  public :: active_tree_level_bounds
  public :: active_tree_level_node
  public :: append_i32_buffer
  public :: nodes_well_separated

contains

  pure integer(i32) function octant_index(x, y, z, center)
    real(dp), intent(in) :: x, y, z, center(3)

    octant_index = 1_i32
    if (x >= center(1)) octant_index = octant_index + 1_i32
    if (y >= center(2)) octant_index = octant_index + 2_i32
    if (z >= center(3)) octant_index = octant_index + 4_i32
  end function octant_index

  pure integer(i32) function active_tree_nnode(plan, use_target_tree)
    type(fmm_plan_type), intent(in) :: plan
    logical, intent(in) :: use_target_tree

    active_tree_nnode = merge(plan%target_nnode, plan%nnode, use_target_tree)
  end function active_tree_nnode

  pure integer(i32) function active_tree_child_count(plan, use_target_tree, node_idx)
    type(fmm_plan_type), intent(in) :: plan
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: node_idx

    if (use_target_tree) then
      active_tree_child_count = plan%target_child_count(node_idx)
    else
      active_tree_child_count = plan%child_count(node_idx)
    end if
  end function active_tree_child_count

  pure integer(i32) function active_tree_child_idx(plan, use_target_tree, child_k, node_idx)
    type(fmm_plan_type), intent(in) :: plan
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: child_k, node_idx

    if (use_target_tree) then
      active_tree_child_idx = plan%target_child_idx(child_k, node_idx)
    else
      active_tree_child_idx = plan%child_idx(child_k, node_idx)
    end if
  end function active_tree_child_idx

  pure integer(i32) function active_tree_child_octant(plan, use_target_tree, child_k, node_idx)
    type(fmm_plan_type), intent(in) :: plan
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: child_k, node_idx

    if (use_target_tree) then
      active_tree_child_octant = plan%target_child_octant(child_k, node_idx)
    else
      active_tree_child_octant = plan%child_octant(child_k, node_idx)
    end if
  end function active_tree_child_octant

  pure function active_tree_node_center(plan, use_target_tree, node_idx) result(center)
    type(fmm_plan_type), intent(in) :: plan
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: node_idx
    real(dp) :: center(3)

    if (use_target_tree) then
      center = plan%target_node_center(:, node_idx)
    else
      center = plan%node_center(:, node_idx)
    end if
  end function active_tree_node_center

  pure function active_tree_node_half_size(plan, use_target_tree, node_idx) result(half_size)
    type(fmm_plan_type), intent(in) :: plan
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: node_idx
    real(dp) :: half_size(3)

    if (use_target_tree) then
      half_size = plan%target_node_half_size(:, node_idx)
    else
      half_size = plan%node_half_size(:, node_idx)
    end if
  end function active_tree_node_half_size

  pure real(dp) function active_tree_node_radius(plan, use_target_tree, node_idx)
    type(fmm_plan_type), intent(in) :: plan
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: node_idx

    if (use_target_tree) then
      active_tree_node_radius = plan%target_node_radius(node_idx)
    else
      active_tree_node_radius = plan%node_radius(node_idx)
    end if
  end function active_tree_node_radius

  pure integer(i32) function active_tree_max_depth(plan, use_target_tree)
    type(fmm_plan_type), intent(in) :: plan
    logical, intent(in) :: use_target_tree

    active_tree_max_depth = merge(plan%target_node_max_depth, plan%node_max_depth, use_target_tree)
  end function active_tree_max_depth

  pure subroutine active_tree_level_bounds(plan, use_target_tree, depth, level_start_pos, level_end_pos)
    type(fmm_plan_type), intent(in) :: plan
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: depth
    integer(i32), intent(out) :: level_start_pos, level_end_pos

    if (use_target_tree) then
      level_start_pos = plan%target_level_start(depth + 1_i32)
      level_end_pos = plan%target_level_start(depth + 2_i32) - 1_i32
    else
      level_start_pos = plan%node_level_start(depth + 1_i32)
      level_end_pos = plan%node_level_start(depth + 2_i32) - 1_i32
    end if
  end subroutine active_tree_level_bounds

  pure integer(i32) function active_tree_level_node(plan, use_target_tree, level_pos)
    type(fmm_plan_type), intent(in) :: plan
    logical, intent(in) :: use_target_tree
    integer(i32), intent(in) :: level_pos

    if (use_target_tree) then
      active_tree_level_node = plan%target_level_nodes(level_pos)
    else
      active_tree_level_node = plan%node_level_nodes(level_pos)
    end if
  end function active_tree_level_node

  subroutine append_i32_buffer(buf, n_used, capacity, value)
    integer(i32), allocatable, intent(inout) :: buf(:)
    integer(i32), intent(inout) :: n_used, capacity
    integer(i32), intent(in) :: value
    integer(i32), allocatable :: tmp(:)
    integer(i32) :: new_capacity

    if (n_used >= capacity) then
      new_capacity = max(capacity*2_i32, capacity + 32_i32)
      allocate (tmp(new_capacity))
      tmp = 0_i32
      if (n_used > 0_i32) tmp(1:n_used) = buf(1:n_used)
      call move_alloc(tmp, buf)
      capacity = new_capacity
    end if
    n_used = n_used + 1_i32
    buf(n_used) = value
  end subroutine append_i32_buffer

  logical function nodes_well_separated(plan, target_node, source_node)
    type(fmm_plan_type), intent(in) :: plan
    integer(i32), intent(in) :: target_node, source_node
    real(dp) :: d(3), dist2, rs, rt, theta_eff, lhs, rhs, target_center(3)
    logical :: use_target_tree

    use_target_tree = plan%target_tree_ready
    target_center = active_tree_node_center(plan, use_target_tree, target_node)
    rt = active_tree_node_radius(plan, use_target_tree, target_node)

    d = target_center - plan%node_center(:, source_node)
    call apply_periodic2_minimum_image(plan, d)
    dist2 = sum(d*d)
    if (dist2 <= 0.0d0) then
      nodes_well_separated = .false.
      return
    end if

    rs = plan%node_radius(source_node)
    theta_eff = plan%options%theta
    lhs = (rs + rt)*(rs + rt)
    rhs = (theta_eff*theta_eff)*dist2
    nodes_well_separated = (lhs < rhs)
  end function nodes_well_separated

end module bem_coulomb_fmm_tree_utils

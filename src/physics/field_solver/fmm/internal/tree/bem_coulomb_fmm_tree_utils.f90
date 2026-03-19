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

  !> 座標が親ノード中心のどの八分木に属するかを返す。
  !! @param[in] x 判定する x 座標。
  !! @param[in] y 判定する y 座標。
  !! @param[in] z 判定する z 座標。
  !! @param[in] center ノード中心座標。
  !! @return 八分木の番号。
  pure integer(i32) function octant_index(x, y, z, center)
    real(dp), intent(in) :: x, y, z, center(3)

    octant_index = 1_i32
    if (x >= center(1)) octant_index = octant_index + 1_i32
    if (y >= center(2)) octant_index = octant_index + 2_i32
    if (z >= center(3)) octant_index = octant_index + 4_i32
  end function octant_index

  !> 現在有効な木のノード数を返す。
  !! @param[in] plan FMM 計画。
  !! @param[in] use_target_tree target 木を使うなら `.true.`。
  !! @return ノード数。
  pure integer(i32) function active_tree_nnode(plan, use_target_tree)
    type(fmm_plan_type), intent(in) :: plan
    logical, intent(in) :: use_target_tree

    active_tree_nnode = merge(plan%target_nnode, plan%nnode, use_target_tree)
  end function active_tree_nnode

  !> 有効な木の指定ノードにある子ノード数を返す。
  !! @param[in] plan FMM 計画。
  !! @param[in] use_target_tree target 木を使うなら `.true.`。
  !! @param[in] node_idx ノード番号。
  !! @return 子ノード数。
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

  !> 有効な木の子ノード番号を返す。
  !! @param[in] plan FMM 計画。
  !! @param[in] use_target_tree target 木を使うなら `.true.`。
  !! @param[in] child_k 子の連番。
  !! @param[in] node_idx 親ノード番号。
  !! @return 子ノード番号。
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

  !> 有効な木の子ノードが属する八分木番号を返す。
  !! @param[in] plan FMM 計画。
  !! @param[in] use_target_tree target 木を使うなら `.true.`。
  !! @param[in] child_k 子の連番。
  !! @param[in] node_idx 親ノード番号。
  !! @return 八分木番号。
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

  !> 有効な木のノード中心座標を返す。
  !! @param[in] plan FMM 計画。
  !! @param[in] use_target_tree target 木を使うなら `.true.`。
  !! @param[in] node_idx ノード番号。
  !! @return ノード中心座標。
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

  !> 有効な木のノード半サイズを返す。
  !! @param[in] plan FMM 計画。
  !! @param[in] use_target_tree target 木を使うなら `.true.`。
  !! @param[in] node_idx ノード番号。
  !! @return ノード半サイズ。
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

  !> 有効な木のノード外接半径を返す。
  !! @param[in] plan FMM 計画。
  !! @param[in] use_target_tree target 木を使うなら `.true.`。
  !! @param[in] node_idx ノード番号。
  !! @return 外接半径。
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

  !> 有効な木の最大深さを返す。
  !! @param[in] plan FMM 計画。
  !! @param[in] use_target_tree target 木を使うなら `.true.`。
  !! @return 最大深さ。
  pure integer(i32) function active_tree_max_depth(plan, use_target_tree)
    type(fmm_plan_type), intent(in) :: plan
    logical, intent(in) :: use_target_tree

    active_tree_max_depth = merge(plan%target_node_max_depth, plan%node_max_depth, use_target_tree)
  end function active_tree_max_depth

  !> 指定深さのレベル範囲を返す。
  !! @param[in] plan FMM 計画。
  !! @param[in] use_target_tree target 木を使うなら `.true.`。
  !! @param[in] depth 深さ。
  !! @param[out] level_start_pos レベル先頭位置。
  !! @param[out] level_end_pos レベル末尾位置。
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

  !> レベル配列からノード番号を返す。
  !! @param[in] plan FMM 計画。
  !! @param[in] use_target_tree target 木を使うなら `.true.`。
  !! @param[in] level_pos レベル内位置。
  !! @return ノード番号。
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

  !> 整数バッファへ値を追加し、必要なら容量を拡張する。
  !! @param[inout] buf 追加先バッファ。
  !! @param[inout] n_used 使用中要素数。
  !! @param[inout] capacity 確保容量。
  !! @param[in] value 追加する値。
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

  !> target/source ノードが十分離れているかを判定する。
  !! @param[in] plan FMM 計画。
  !! @param[in] target_node target ノード番号。
  !! @param[in] source_node source ノード番号。
  !! @return 十分離れていれば `.true.`。
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

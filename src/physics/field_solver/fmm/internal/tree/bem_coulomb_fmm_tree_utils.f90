!> Coulomb FMM tree 構造の共通ユーティリティ。
module bem_coulomb_fmm_tree_utils
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_types, only: fmm_plan_type, fmm_tree_type
  use bem_coulomb_fmm_periodic, only: apply_periodic2_minimum_image
  implicit none
  private

  public :: octant_index
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
  !! @param[in] target_tree 評価用の木（target_tree_ready なら plan%target_tree、そうでなければ plan%source_tree）。
  !! @param[in] target_node target ノード番号。
  !! @param[in] source_node source ノード番号。
  !! @return 十分離れていれば `.true.`。
  logical function nodes_well_separated(plan, target_tree, target_node, source_node)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_tree_type), intent(in) :: target_tree
    integer(i32), intent(in) :: target_node, source_node
    real(dp) :: d(3), dist2, rs, rt, theta_eff, lhs, rhs, target_center(3)

    target_center = target_tree%node_center(:, target_node)
    rt = target_tree%node_radius(target_node)

    d = target_center - plan%source_tree%node_center(:, source_node)
    call apply_periodic2_minimum_image(plan, d)
    dist2 = sum(d*d)
    if (dist2 <= 0.0d0) then
      nodes_well_separated = .false.
      return
    end if

    rs = plan%source_tree%node_radius(source_node)
    theta_eff = plan%options%theta
    lhs = (rs + rt)*(rs + rt)
    rhs = (theta_eff*theta_eff)*dist2
    nodes_well_separated = (lhs < rhs)
  end function nodes_well_separated

end module bem_coulomb_fmm_tree_utils

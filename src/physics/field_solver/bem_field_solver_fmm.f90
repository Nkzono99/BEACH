!> `bem_field_solver` の FMM 相互作用構築と評価を実装する submodule。
submodule (bem_field_solver) bem_field_solver_fmm
  implicit none
contains

  !> source/target ノード球が十分離れていれば遠方相互作用として扱う。
  module procedure nodes_well_separated
    real(dp) :: dx, dy, dz, dist, dist2, rs, rt, d_min

    dx = self%node_center(1, target_node) - self%node_center(1, source_node)
    dy = self%node_center(2, target_node) - self%node_center(2, source_node)
    dz = self%node_center(3, target_node) - self%node_center(3, source_node)
    dist2 = dx * dx + dy * dy + dz * dz
    if (dist2 <= 0.0d0) then
      accept_it = .false.
      return
    end if

    dist = sqrt(dist2)
    rt = self%node_radius(target_node)
    rs = self%node_radius(source_node)
    d_min = dist - rt
    if (d_min <= rs) then
      accept_it = .false.
      return
    end if
    accept_it = (rs < self%theta * (d_min - rs))
  end procedure nodes_well_separated

  !> 1つの葉ノードに対する near/far リストを木走査で構築する。
  module procedure gather_leaf_interactions
    integer(i32) :: child_k, child_node

    if (source_node == target_leaf) then
      near_n = near_n + 1_i32
      near_buf(near_n) = source_node
      return
    end if

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
    integer(i32) :: node_idx, leaf_idx, i
    integer(i32) :: near_n, far_n, near_used, far_used
    integer(i32) :: near_cap, far_cap
    integer(i32), allocatable :: near_work(:), far_work(:), near_buf(:), far_buf(:)

    if (.not. self%tree_ready .or. self%nnode <= 0_i32) then
      self%fmm_ready = .false.
      self%nleaf = 0_i32
      return
    end if

    self%nleaf = 0_i32
    do node_idx = 1_i32, self%nnode
      if (self%child_count(node_idx) <= 0_i32) self%nleaf = self%nleaf + 1_i32
    end do
    if (self%nleaf <= 0_i32) then
      self%fmm_ready = .false.
      return
    end if

    if (allocated(self%leaf_nodes)) deallocate (self%leaf_nodes)
    if (allocated(self%leaf_slot_of_node)) deallocate (self%leaf_slot_of_node)
    if (allocated(self%near_start)) deallocate (self%near_start)
    if (allocated(self%far_start)) deallocate (self%far_start)
    if (allocated(self%near_nodes)) deallocate (self%near_nodes)
    if (allocated(self%far_nodes)) deallocate (self%far_nodes)
    allocate (self%leaf_nodes(self%nleaf), self%leaf_slot_of_node(self%max_node))
    allocate (self%near_start(self%nleaf + 1_i32), self%far_start(self%nleaf + 1_i32))
    self%leaf_slot_of_node = 0_i32

    leaf_idx = 0_i32
    do node_idx = 1_i32, self%nnode
      if (self%child_count(node_idx) > 0_i32) cycle
      leaf_idx = leaf_idx + 1_i32
      self%leaf_nodes(leaf_idx) = node_idx
      self%leaf_slot_of_node(node_idx) = leaf_idx
    end do

    allocate (near_work(self%nnode), far_work(self%nnode))
    near_cap = max(64_i32, self%nleaf * 16_i32)
    far_cap = max(64_i32, self%nleaf * 16_i32)
    allocate (near_buf(near_cap), far_buf(far_cap))
    near_used = 0_i32
    far_used = 0_i32

    do leaf_idx = 1_i32, self%nleaf
      near_n = 0_i32
      far_n = 0_i32
      call gather_leaf_interactions( &
        self, self%leaf_nodes(leaf_idx), 1_i32, near_work, near_n, far_work, far_n &
      )
      self%near_start(leaf_idx) = near_used + 1_i32
      self%far_start(leaf_idx) = far_used + 1_i32

      do i = 1_i32, near_n
        call append_i32(near_buf, near_used, near_cap, near_work(i))
      end do
      do i = 1_i32, far_n
        call append_i32(far_buf, far_used, far_cap, far_work(i))
      end do
    end do
    self%near_start(self%nleaf + 1_i32) = near_used + 1_i32
    self%far_start(self%nleaf + 1_i32) = far_used + 1_i32

    allocate (self%near_nodes(max(1_i32, near_used)))
    self%near_nodes = 0_i32
    if (near_used > 0_i32) self%near_nodes(1:near_used) = near_buf(1:near_used)
    allocate (self%far_nodes(max(1_i32, far_used)))
    self%far_nodes = 0_i32
    if (far_used > 0_i32) self%far_nodes(1:far_used) = far_buf(1:far_used)

    if (allocated(self%leaf_far_e0)) deallocate (self%leaf_far_e0)
    if (allocated(self%leaf_far_jac)) deallocate (self%leaf_far_jac)
    allocate (self%leaf_far_e0(3, self%nleaf), self%leaf_far_jac(3, 3, self%nleaf))
    self%leaf_far_e0 = 0.0d0
    self%leaf_far_jac = 0.0d0
    self%fmm_ready = .true.

    deallocate (near_work, far_work, near_buf, far_buf)

  contains

    subroutine append_i32(buf, n_used, capacity, value)
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
    end subroutine append_i32

  end procedure build_fmm_interactions

  !> 遠方ノードの monopole から葉中心での局所電場/ヤコビアンを更新する。
  module procedure refresh_fmm_locals
    integer(i32) :: leaf_slot, src_pos, src_node
    integer(i32) :: src_beg, src_end
    real(dp) :: ct(3), dx, dy, dz, r2, inv_r, inv_r3, inv_r5
    real(dp) :: ex, ey, ez, qi, c0, c1
    real(dp) :: jxx, jxy, jxz, jyy, jyz, jzz
    real(dp) :: soft2

    if (trim(self%mode) /= 'fmm') return
    if (.not. self%fmm_ready) return
    if (self%nleaf <= 0_i32) return

    soft2 = self%softening * self%softening
    self%leaf_far_e0 = 0.0d0
    self%leaf_far_jac = 0.0d0

    do leaf_slot = 1_i32, self%nleaf
      ct = self%node_center(:, self%leaf_nodes(leaf_slot))
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      jxx = 0.0d0
      jxy = 0.0d0
      jxz = 0.0d0
      jyy = 0.0d0
      jyz = 0.0d0
      jzz = 0.0d0

      src_beg = self%far_start(leaf_slot)
      src_end = self%far_start(leaf_slot + 1_i32) - 1_i32
      do src_pos = src_beg, src_end
        src_node = self%far_nodes(src_pos)
        qi = self%node_q(src_node)
        if (abs(qi) <= tiny(1.0d0)) cycle
        dx = ct(1) - self%node_charge_center(1, src_node)
        dy = ct(2) - self%node_charge_center(2, src_node)
        dz = ct(3) - self%node_charge_center(3, src_node)
        r2 = dx * dx + dy * dy + dz * dz + soft2
        inv_r = 1.0d0 / sqrt(r2)
        inv_r3 = inv_r / r2
        inv_r5 = inv_r3 / r2

        ex = ex + qi * inv_r3 * dx
        ey = ey + qi * inv_r3 * dy
        ez = ez + qi * inv_r3 * dz

        c0 = qi * inv_r3
        c1 = 3.0d0 * qi * inv_r5
        jxx = jxx + c0 - c1 * dx * dx
        jxy = jxy - c1 * dx * dy
        jxz = jxz - c1 * dx * dz
        jyy = jyy + c0 - c1 * dy * dy
        jyz = jyz - c1 * dy * dz
        jzz = jzz + c0 - c1 * dz * dz
      end do

      self%leaf_far_e0(:, leaf_slot) = [ex, ey, ez]
      self%leaf_far_jac(1, 1, leaf_slot) = jxx
      self%leaf_far_jac(1, 2, leaf_slot) = jxy
      self%leaf_far_jac(1, 3, leaf_slot) = jxz
      self%leaf_far_jac(2, 1, leaf_slot) = jxy
      self%leaf_far_jac(2, 2, leaf_slot) = jyy
      self%leaf_far_jac(2, 3, leaf_slot) = jyz
      self%leaf_far_jac(3, 1, leaf_slot) = jxz
      self%leaf_far_jac(3, 2, leaf_slot) = jyz
      self%leaf_far_jac(3, 3, leaf_slot) = jzz
    end do
  end procedure refresh_fmm_locals

  !> 観測点の octant を辿って対応する葉ノードを返す。
  module procedure locate_target_leaf
    integer(i32) :: oct, child_k, child, cand
    real(dp) :: best_dist2, dist2, dx, dy, dz

    node_idx = 0_i32
    if (.not. self%tree_ready .or. self%nnode <= 0_i32) return

    node_idx = 1_i32
    do while (self%child_count(node_idx) > 0_i32)
      oct = octant_index(r(1), r(2), r(3), self%node_center(:, node_idx))
      child = 0_i32
      do child_k = 1_i32, self%child_count(node_idx)
        if (self%child_octant(child_k, node_idx) == oct) then
          child = self%child_idx(child_k, node_idx)
          exit
        end if
      end do

      if (child <= 0_i32) then
        child = self%child_idx(1, node_idx)
        best_dist2 = huge(1.0d0)
        do child_k = 1_i32, self%child_count(node_idx)
          cand = self%child_idx(child_k, node_idx)
          dx = r(1) - self%node_center(1, cand)
          dy = r(2) - self%node_center(2, cand)
          dz = r(3) - self%node_center(3, cand)
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

  !> FMM 局所展開 + 近傍 direct 和で観測点の電場を返す。
  module procedure eval_e_fmm
    integer(i32) :: leaf_node, leaf_slot
    integer(i32) :: far_pos, far_node
    integer(i32) :: near_pos, near_node, p, idx, p_end
    real(dp) :: soft2, r2, inv_r3, dx, dy, dz, qi
    real(dp) :: ex, ey, ez

    leaf_node = locate_target_leaf(self, r)
    if (leaf_node <= 0_i32) then
      call electric_field_at(mesh, r, self%softening, e)
      return
    end if
    leaf_slot = self%leaf_slot_of_node(leaf_node)
    if (leaf_slot <= 0_i32) then
      call electric_field_at(mesh, r, self%softening, e)
      return
    end if

    soft2 = self%softening * self%softening
    ex = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    do far_pos = self%far_start(leaf_slot), self%far_start(leaf_slot + 1_i32) - 1_i32
      far_node = self%far_nodes(far_pos)
      qi = self%node_q(far_node)
      if (abs(qi) <= tiny(1.0d0)) cycle
      dx = r(1) - self%node_charge_center(1, far_node)
      dy = r(2) - self%node_charge_center(2, far_node)
      dz = r(3) - self%node_charge_center(3, far_node)
      r2 = dx * dx + dy * dy + dz * dz + soft2
      inv_r3 = 1.0d0 / (sqrt(r2) * r2)
      ex = ex + qi * inv_r3 * dx
      ey = ey + qi * inv_r3 * dy
      ez = ez + qi * inv_r3 * dz
    end do

    do near_pos = self%near_start(leaf_slot), self%near_start(leaf_slot + 1_i32) - 1_i32
      near_node = self%near_nodes(near_pos)
      p_end = self%node_start(near_node) + self%node_count(near_node) - 1_i32
      do p = self%node_start(near_node), p_end
        idx = self%elem_order(p)
        dx = r(1) - mesh%center_x(idx)
        dy = r(2) - mesh%center_y(idx)
        dz = r(3) - mesh%center_z(idx)
        r2 = dx * dx + dy * dy + dz * dz + soft2
        inv_r3 = 1.0d0 / (sqrt(r2) * r2)
        qi = mesh%q_elem(idx)
        ex = ex + qi * inv_r3 * dx
        ey = ey + qi * inv_r3 * dy
        ez = ez + qi * inv_r3 * dz
      end do
    end do

    e(1) = k_coulomb * ex
    e(2) = k_coulomb * ey
    e(3) = k_coulomb * ez
  end procedure eval_e_fmm

end submodule bem_field_solver_fmm

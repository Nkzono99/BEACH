!> `bem_field_solver` の FMM 相互作用構築と評価を実装する submodule。
submodule (bem_field_solver) bem_field_solver_fmm
  implicit none
  real(dp), parameter :: inv_sqrt_pi = 0.56418958354775628695d0
contains

  !> source/target ノード球が十分離れていれば遠方相互作用として扱う。
  module procedure nodes_well_separated
    integer(i32) :: k, axis
    real(dp) :: d(3), dist, dist2, rs, rt, d_min, theta_eff

    d(1) = self%node_center(1, target_node) - self%node_center(1, source_node)
    d(2) = self%node_center(2, target_node) - self%node_center(2, source_node)
    d(3) = self%node_center(3, target_node) - self%node_center(3, source_node)
    if (self%use_periodic2) then
      do k = 1_i32, 2_i32
        axis = self%periodic_axes(k)
        d(axis) = d(axis) - self%periodic_len(k) * dnint(d(axis) / self%periodic_len(k))
      end do
    end if
    dist2 = sum(d * d)
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
    theta_eff = self%theta
    if (self%use_periodic2) theta_eff = 0.2d0 * self%theta
    accept_it = (rs < theta_eff * (d_min - rs))
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
    if (allocated(self%leaf_far_hess)) deallocate (self%leaf_far_hess)
    allocate (self%leaf_far_e0(3, self%nleaf), self%leaf_far_jac(3, 3, self%nleaf), self%leaf_far_hess(3, 3, 3, self%nleaf))
    self%leaf_far_e0 = 0.0d0
    self%leaf_far_jac = 0.0d0
    self%leaf_far_hess = 0.0d0
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

  !> M2L と L2L で葉ノード局所展開（E/Jac/Hess）を更新する。
  module procedure refresh_fmm_locals
    integer(i32) :: leaf_slot, node_idx, parent_node, child_k
    integer(i32) :: i, j, k
    integer(i32) :: img_i, img_j, img_min, img_max, axis1, axis2
    integer(i32), allocatable :: parent_of(:)
    real(dp), allocatable :: node_local_e0(:, :), node_local_jac(:, :, :), node_local_hess(:, :, :, :)
    real(dp) :: shift1, shift2
    real(dp) :: dr(3), dvec(3), r2, inv_r, inv_r3, inv_r5, inv_r7
    real(dp) :: qi
    real(dp) :: delta_ij, delta_ik, delta_jk
    real(dp) :: soft2

    if (trim(self%mode) /= 'fmm') return
    if (.not. self%fmm_ready) return
    if (self%nleaf <= 0_i32) return

    soft2 = self%softening * self%softening
    self%leaf_far_e0 = 0.0d0
    self%leaf_far_jac = 0.0d0
    self%leaf_far_hess = 0.0d0
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

    allocate ( &
      parent_of(self%nnode), &
      node_local_e0(3, self%nnode), &
      node_local_jac(3, 3, self%nnode), &
      node_local_hess(3, 3, 3, self%nnode) &
    )
    parent_of = 0_i32
    node_local_e0 = 0.0d0
    node_local_jac = 0.0d0
    node_local_hess = 0.0d0

    do parent_node = 1_i32, self%nnode
      do child_k = 1_i32, self%child_count(parent_node)
        node_idx = self%child_idx(child_k, parent_node)
        parent_of(node_idx) = parent_node
      end do
    end do

    call accumulate_m2l_pairs(1_i32, 1_i32)

    ! L2L: 親ノード局所展開を子ノード中心へ平行移動する。
    do node_idx = 2_i32, self%nnode
      parent_node = parent_of(node_idx)
      if (parent_node <= 0_i32) cycle
      dr = self%node_center(:, node_idx) - self%node_center(:, parent_node)
      do i = 1_i32, 3_i32
        node_local_e0(i, node_idx) = node_local_e0(i, node_idx) + node_local_e0(i, parent_node)
        do j = 1_i32, 3_i32
          node_local_e0(i, node_idx) = node_local_e0(i, node_idx) + node_local_jac(i, j, parent_node) * dr(j)
          do k = 1_i32, 3_i32
            node_local_e0(i, node_idx) = node_local_e0(i, node_idx) + 0.5d0 * node_local_hess(i, j, k, parent_node) * dr(j) * dr(k)
            node_local_jac(i, j, node_idx) = node_local_jac(i, j, node_idx) + node_local_hess(i, j, k, parent_node) * dr(k)
          end do
          node_local_jac(i, j, node_idx) = node_local_jac(i, j, node_idx) + node_local_jac(i, j, parent_node)
          do k = 1_i32, 3_i32
            node_local_hess(i, j, k, node_idx) = node_local_hess(i, j, k, node_idx) + node_local_hess(i, j, k, parent_node)
          end do
        end do
      end do
    end do

    do leaf_slot = 1_i32, self%nleaf
      node_idx = self%leaf_nodes(leaf_slot)
      self%leaf_far_e0(:, leaf_slot) = node_local_e0(:, node_idx)
      self%leaf_far_jac(:, :, leaf_slot) = node_local_jac(:, :, node_idx)
      self%leaf_far_hess(:, :, :, leaf_slot) = node_local_hess(:, :, :, node_idx)
    end do

    deallocate (parent_of, node_local_e0, node_local_jac, node_local_hess)

  contains

    recursive subroutine accumulate_m2l_pairs(t_node, s_node)
      integer(i32), intent(in) :: t_node, s_node
      integer(i32) :: tc, sc, t_child, s_child

      if (t_node == s_node) then
        if (self%child_count(t_node) <= 0_i32) return
        do tc = 1_i32, self%child_count(t_node)
          t_child = self%child_idx(tc, t_node)
          do sc = 1_i32, self%child_count(s_node)
            s_child = self%child_idx(sc, s_node)
            call accumulate_m2l_pairs(t_child, s_child)
          end do
        end do
        return
      end if

      if (nodes_well_separated(self, t_node, s_node)) then
        call add_m2l_from_monopole(t_node, s_node)
        return
      end if

      if (self%child_count(t_node) <= 0_i32 .and. self%child_count(s_node) <= 0_i32) return

      if (self%child_count(s_node) > 0_i32 .and. &
          (self%child_count(t_node) <= 0_i32 .or. self%node_radius(s_node) >= self%node_radius(t_node))) then
        do sc = 1_i32, self%child_count(s_node)
          s_child = self%child_idx(sc, s_node)
          call accumulate_m2l_pairs(t_node, s_child)
        end do
      else
        do tc = 1_i32, self%child_count(t_node)
          t_child = self%child_idx(tc, t_node)
          call accumulate_m2l_pairs(t_child, s_node)
        end do
      end if
    end subroutine accumulate_m2l_pairs

    subroutine add_m2l_from_monopole(t_node, s_node)
      integer(i32), intent(in) :: t_node, s_node

      qi = self%node_q(s_node)
      if (abs(qi) <= tiny(1.0d0)) return

      do img_i = img_min, img_max
        do img_j = img_min, img_max
          dvec = self%node_center(:, t_node) - self%node_charge_center(:, s_node)
          if (self%use_periodic2) then
            shift1 = real(img_i, dp) * self%periodic_len(1)
            shift2 = real(img_j, dp) * self%periodic_len(2)
            if (axis1 > 0_i32) dvec(axis1) = dvec(axis1) - shift1
            if (axis2 > 0_i32) dvec(axis2) = dvec(axis2) - shift2
          end if

          r2 = sum(dvec * dvec) + soft2
          inv_r = 1.0d0 / sqrt(r2)
          inv_r3 = inv_r / r2
          inv_r5 = inv_r3 / r2
          inv_r7 = inv_r5 / r2

          do i = 1_i32, 3_i32
            node_local_e0(i, t_node) = node_local_e0(i, t_node) + qi * inv_r3 * dvec(i)
            do j = 1_i32, 3_i32
              delta_ij = merge(1.0d0, 0.0d0, i == j)
              node_local_jac(i, j, t_node) = node_local_jac(i, j, t_node) &
                                             + qi * (delta_ij * inv_r3 - 3.0d0 * dvec(i) * dvec(j) * inv_r5)
              do k = 1_i32, 3_i32
                delta_ik = merge(1.0d0, 0.0d0, i == k)
                delta_jk = merge(1.0d0, 0.0d0, j == k)
                node_local_hess(i, j, k, t_node) = node_local_hess(i, j, k, t_node) &
                                                   + qi * ( &
                                                       15.0d0 * dvec(i) * dvec(j) * dvec(k) * inv_r7 &
                                                       - 3.0d0 * ( &
                                                           delta_ij * dvec(k) + delta_ik * dvec(j) + delta_jk * dvec(i) &
                                                         ) * inv_r5 &
                                                     )
              end do
            end do
          end do
        end do
      end do
    end subroutine add_m2l_from_monopole
  end procedure refresh_fmm_locals

  !> 観測点の octant を辿って対応する葉ノードを返す。
  module procedure locate_target_leaf
    integer(i32) :: oct, child_k, child, cand
    integer(i32) :: axis
    real(dp) :: best_dist2, dist2, dx, dy, dz
    real(dp) :: extent_eps

    node_idx = 0_i32
    if (.not. self%tree_ready .or. self%nnode <= 0_i32) return

    if (.not. self%use_periodic2) then
      extent_eps = 1.0d-12 * max(1.0d0, maxval(abs(self%node_center(:, 1_i32))))
      do axis = 1_i32, 3_i32
        if (abs(r(axis) - self%node_center(axis, 1_i32)) > self%node_half_size(axis, 1_i32) + extent_eps) return
      end do
    end if

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
    integer(i32) :: near_pos, near_node, p, idx, p_end
    integer(i32) :: img_i, img_j, img_min, img_max
    integer(i32) :: axis1, axis2
    integer(i32) :: j, k
    real(dp) :: soft2, r2, inv_r3, dx, dy, dz, qi
    real(dp) :: shift1, shift2
    real(dp) :: src_x, src_y, src_z
    real(dp) :: ex, ey, ez
    real(dp) :: dr(3)

    soft2 = self%softening * self%softening
    leaf_node = locate_target_leaf(self, r)
    if (leaf_node <= 0_i32) then
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      call traverse_node(self, mesh, 1_i32, r(1), r(2), r(3), soft2, ex, ey, ez)
      e(1) = k_coulomb * ex
      e(2) = k_coulomb * ey
      e(3) = k_coulomb * ez
      return
    end if
    leaf_slot = self%leaf_slot_of_node(leaf_node)
    if (leaf_slot <= 0_i32) then
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      call traverse_node(self, mesh, 1_i32, r(1), r(2), r(3), soft2, ex, ey, ez)
      e(1) = k_coulomb * ex
      e(2) = k_coulomb * ey
      e(3) = k_coulomb * ez
      return
    end if

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
    if (self%use_periodic2) then
      if (.not. point_in_node_periodic_box(self, leaf_node, r)) then
        call eval_periodic2_tree_target(self, mesh, r, ex, ey, ez)
        if (trim(self%periodic_far_correction) == 'ewald_like') then
          call add_periodic2_ewald_like_correction_all_leaves(self, r, ex, ey, ez)
        end if
        e(1) = k_coulomb * ex
        e(2) = k_coulomb * ey
        e(3) = k_coulomb * ez
        return
      end if
    end if
    ex = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    dr = r - self%node_center(:, leaf_node)
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
        do img_i = img_min, img_max
          do img_j = img_min, img_max
            shift1 = 0.0d0
            shift2 = 0.0d0
            if (self%use_periodic2) then
              shift1 = real(img_i, dp) * self%periodic_len(1)
              shift2 = real(img_j, dp) * self%periodic_len(2)
            end if
            src_x = mesh%center_x(idx)
            src_y = mesh%center_y(idx)
            src_z = mesh%center_z(idx)
            if (axis1 == 1_i32) src_x = src_x + shift1
            if (axis1 == 2_i32) src_y = src_y + shift1
            if (axis1 == 3_i32) src_z = src_z + shift1
            if (axis2 == 1_i32) src_x = src_x + shift2
            if (axis2 == 2_i32) src_y = src_y + shift2
            if (axis2 == 3_i32) src_z = src_z + shift2
            dx = r(1) - src_x
            dy = r(2) - src_y
            dz = r(3) - src_z
            r2 = dx * dx + dy * dy + dz * dz + soft2
            inv_r3 = 1.0d0 / (sqrt(r2) * r2)
            ex = ex + qi * inv_r3 * dx
            ey = ey + qi * inv_r3 * dy
            ez = ez + qi * inv_r3 * dz
          end do
        end do
      end do
    end do
    if (self%use_periodic2 .and. trim(self%periodic_far_correction) == 'ewald_like') then
      call add_periodic2_ewald_like_correction(self, leaf_slot, r, ex, ey, ez)
    end if

    e(1) = k_coulomb * ex
    e(2) = k_coulomb * ey
    e(3) = k_coulomb * ez
  end procedure eval_e_fmm

  !> periodic2 で target が葉ノード領域外にある場合の M2P/treecode フォールバック。
  subroutine eval_periodic2_tree_target(self, mesh, r, ex, ey, ez)
    class(field_solver_type), intent(in) :: self
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: r(3)
    real(dp), intent(out) :: ex, ey, ez

    integer(i32) :: axis1, axis2, img_min, img_max
    real(dp) :: soft2

    axis1 = self%periodic_axes(1)
    axis2 = self%periodic_axes(2)
    img_min = -self%periodic_image_layers
    img_max = self%periodic_image_layers
    soft2 = self%softening * self%softening
    ex = 0.0d0
    ey = 0.0d0
    ez = 0.0d0

    call accumulate_node(1_i32)

  contains

    recursive subroutine accumulate_node(node_idx)
      integer(i32), intent(in) :: node_idx
      integer(i32) :: child_k, p, idx, p_end
      integer(i32) :: img_i, img_j
      real(dp) :: qi, shift1, shift2, src_x, src_y, src_z
      real(dp) :: dx, dy, dz, r2, inv_r3

      if (self%child_count(node_idx) <= 0_i32) then
        p_end = self%node_start(node_idx) + self%node_count(node_idx) - 1_i32
        do p = self%node_start(node_idx), p_end
          idx = self%elem_order(p)
          qi = mesh%q_elem(idx)
          if (abs(qi) <= tiny(1.0d0)) cycle
          do img_i = img_min, img_max
            do img_j = img_min, img_max
              shift1 = real(img_i, dp) * self%periodic_len(1)
              shift2 = real(img_j, dp) * self%periodic_len(2)
              src_x = mesh%center_x(idx)
              src_y = mesh%center_y(idx)
              src_z = mesh%center_z(idx)
              if (axis1 == 1_i32) src_x = src_x + shift1
              if (axis1 == 2_i32) src_y = src_y + shift1
              if (axis1 == 3_i32) src_z = src_z + shift1
              if (axis2 == 1_i32) src_x = src_x + shift2
              if (axis2 == 2_i32) src_y = src_y + shift2
              if (axis2 == 3_i32) src_z = src_z + shift2
              dx = r(1) - src_x
              dy = r(2) - src_y
              dz = r(3) - src_z
              r2 = dx * dx + dy * dy + dz * dz + soft2
              inv_r3 = 1.0d0 / (sqrt(r2) * r2)
              ex = ex + qi * inv_r3 * dx
              ey = ey + qi * inv_r3 * dy
              ez = ez + qi * inv_r3 * dz
            end do
          end do
        end do
        return
      end if

      if (accept_node_periodic_target(node_idx)) then
        qi = self%node_q(node_idx)
        if (abs(qi) <= tiny(1.0d0)) return
        do img_i = img_min, img_max
          do img_j = img_min, img_max
            shift1 = real(img_i, dp) * self%periodic_len(1)
            shift2 = real(img_j, dp) * self%periodic_len(2)
            src_x = self%node_charge_center(1, node_idx)
            src_y = self%node_charge_center(2, node_idx)
            src_z = self%node_charge_center(3, node_idx)
            if (axis1 == 1_i32) src_x = src_x + shift1
            if (axis1 == 2_i32) src_y = src_y + shift1
            if (axis1 == 3_i32) src_z = src_z + shift1
            if (axis2 == 1_i32) src_x = src_x + shift2
            if (axis2 == 2_i32) src_y = src_y + shift2
            if (axis2 == 3_i32) src_z = src_z + shift2
            dx = r(1) - src_x
            dy = r(2) - src_y
            dz = r(3) - src_z
            r2 = dx * dx + dy * dy + dz * dz + soft2
            inv_r3 = 1.0d0 / (sqrt(r2) * r2)
            ex = ex + qi * inv_r3 * dx
            ey = ey + qi * inv_r3 * dy
            ez = ez + qi * inv_r3 * dz
          end do
        end do
        return
      end if

      do child_k = 1_i32, self%child_count(node_idx)
        call accumulate_node(self%child_idx(child_k, node_idx))
      end do
    end subroutine accumulate_node

    logical function accept_node_periodic_target(node_idx)
      integer(i32), intent(in) :: node_idx
      real(dp) :: d(3), dist2, dist, radius, theta_eff

      d = r - self%node_center(:, node_idx)
      d(axis1) = d(axis1) - self%periodic_len(1) * dnint(d(axis1) / self%periodic_len(1))
      d(axis2) = d(axis2) - self%periodic_len(2) * dnint(d(axis2) / self%periodic_len(2))
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

  !> 周期軸を考慮して、点 `r` がノード境界内（画像を含む）にあるかを判定する。
  logical function point_in_node_periodic_box(self, node_idx, r)
    class(field_solver_type), intent(in) :: self
    integer(i32), intent(in) :: node_idx
    real(dp), intent(in) :: r(3)

    integer(i32) :: k, axis
    real(dp) :: d(3), eps

    d = r - self%node_center(:, node_idx)
    if (self%use_periodic2) then
      do k = 1_i32, 2_i32
        axis = self%periodic_axes(k)
        d(axis) = d(axis) - self%periodic_len(k) * dnint(d(axis) / self%periodic_len(k))
      end do
    end if

    eps = 1.0d-12 * max(1.0d0, maxval(abs(self%node_center(:, node_idx))))
    point_in_node_periodic_box = .true.
    do axis = 1_i32, 3_i32
      if (abs(d(axis)) > self%node_half_size(axis, node_idx) + eps) then
        point_in_node_periodic_box = .false.
        return
      end if
    end do
  end function point_in_node_periodic_box

  !> periodic2 の近傍画像和で欠落する遠方セル寄与を、erfc スクリーン核で補正する。
  subroutine add_periodic2_ewald_like_correction(self, leaf_slot, r, ex, ey, ez)
    class(field_solver_type), intent(in) :: self
    integer(i32), intent(in) :: leaf_slot
    real(dp), intent(in) :: r(3)
    real(dp), intent(inout) :: ex, ey, ez

    integer(i32) :: node_pos, node_idx
    integer(i32) :: nimg, img_outer
    integer(i32) :: axis1, axis2
    real(dp) :: alpha, qi, src(3)

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

    do node_pos = self%far_start(leaf_slot), self%far_start(leaf_slot + 1_i32) - 1_i32
      node_idx = self%far_nodes(node_pos)
      qi = self%node_q(node_idx)
      if (abs(qi) <= tiny(1.0d0)) cycle
      src = self%node_charge_center(:, node_idx)
      call add_screened_shifted_node_images(qi, src, r, axis1, axis2, self%periodic_len, nimg, img_outer, alpha, ex, ey, ez)
    end do

    do node_pos = self%near_start(leaf_slot), self%near_start(leaf_slot + 1_i32) - 1_i32
      node_idx = self%near_nodes(node_pos)
      qi = self%node_q(node_idx)
      if (abs(qi) <= tiny(1.0d0)) cycle
      src = self%node_charge_center(:, node_idx)
      call add_screened_shifted_node_images(qi, src, r, axis1, axis2, self%periodic_len, nimg, img_outer, alpha, ex, ey, ez)
    end do
  end subroutine add_periodic2_ewald_like_correction

  !> 葉ノード全体を代表点として、periodic2 の遠方画像補正を加算する。
  subroutine add_periodic2_ewald_like_correction_all_leaves(self, r, ex, ey, ez)
    class(field_solver_type), intent(in) :: self
    real(dp), intent(in) :: r(3)
    real(dp), intent(inout) :: ex, ey, ez

    integer(i32) :: leaf_slot, node_idx
    integer(i32) :: nimg, img_outer
    integer(i32) :: axis1, axis2
    real(dp) :: alpha, qi, src(3)

    if (.not. self%use_periodic2) return
    if (self%nleaf <= 0_i32) return
    if (self%periodic_ewald_layers <= 0_i32) return
    alpha = self%periodic_ewald_alpha
    if (alpha <= 0.0d0) return

    nimg = self%periodic_image_layers
    img_outer = nimg + self%periodic_ewald_layers
    if (img_outer <= nimg) return
    axis1 = self%periodic_axes(1)
    axis2 = self%periodic_axes(2)
    if (axis1 <= 0_i32 .or. axis2 <= 0_i32) return

    do leaf_slot = 1_i32, self%nleaf
      node_idx = self%leaf_nodes(leaf_slot)
      qi = self%node_q(node_idx)
      if (abs(qi) <= tiny(1.0d0)) cycle
      src = self%node_charge_center(:, node_idx)
      call add_screened_shifted_node_images(qi, src, r, axis1, axis2, self%periodic_len, nimg, img_outer, alpha, ex, ey, ez)
    end do
  end subroutine add_periodic2_ewald_like_correction_all_leaves

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

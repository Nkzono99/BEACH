!> `bem_field_solver` の電場評価と木走査ロジックを実装する submodule。
submodule (bem_field_solver) bem_field_solver_eval
  use bem_coulomb_fmm_core, only: eval_point
  implicit none
contains

  !> 観測点の電場を direct / treecode / fmm で評価して返す。
  module procedure eval_e_field_solver
    real(dp) :: rx, ry, rz, soft2, ex, ey, ez

    if (trim(self%mode) == 'fmm') then
      if (mesh%nelem <= 0_i32) then
        e = 0.0d0
        return
      end if
      if (.not. self%fmm_core_ready) then
        error stop 'FMM core is not ready. Call solver%init/refresh before eval_e.'
      end if
      call eval_point(self%fmm_core_plan, self%fmm_core_state, r, e)
      e = k_coulomb * e
      return
    end if

    if (trim(self%mode) /= 'treecode' .or. .not. self%tree_ready) then
      call electric_field_at(mesh, r, self%softening, e)
      return
    end if

    rx = r(1)
    ry = r(2)
    rz = r(3)
    soft2 = self%softening * self%softening
    ex = 0.0d0
    ey = 0.0d0
    ez = 0.0d0

    call traverse_node(self, mesh, 1_i32, rx, ry, rz, soft2, ex, ey, ez)

    e(1) = k_coulomb * ex
    e(2) = k_coulomb * ey
    e(3) = k_coulomb * ez
  end procedure eval_e_field_solver

  !> ノードを再帰走査し、葉は direct 総和・遠方は monopole 近似で加算する。
  module procedure traverse_node
    integer(i32) :: child_k, p, idx, p_end
    real(dp) :: dx, dy, dz, r2, inv_r3, qi

    if (self%child_count(node_idx) <= 0_i32) then
      p_end = self%node_start(node_idx) + self%node_count(node_idx) - 1_i32
      do p = self%node_start(node_idx), p_end
        idx = self%elem_order(p)
        dx = rx - mesh%center_x(idx)
        dy = ry - mesh%center_y(idx)
        dz = rz - mesh%center_z(idx)
        r2 = dx * dx + dy * dy + dz * dz + soft2
        inv_r3 = 1.0d0 / (sqrt(r2) * r2)
        qi = mesh%q_elem(idx)
        ex = ex + qi * inv_r3 * dx
        ey = ey + qi * inv_r3 * dy
        ez = ez + qi * inv_r3 * dz
      end do
      return
    end if

    if (accept_node(self, node_idx, rx, ry, rz)) then
      qi = self%node_q(node_idx)
      if (abs(qi) > 0.0d0) then
        dx = rx - self%node_charge_center(1, node_idx)
        dy = ry - self%node_charge_center(2, node_idx)
        dz = rz - self%node_charge_center(3, node_idx)
        r2 = dx * dx + dy * dy + dz * dz + soft2
        inv_r3 = 1.0d0 / (sqrt(r2) * r2)
        ex = ex + qi * inv_r3 * dx
        ey = ey + qi * inv_r3 * dy
        ez = ez + qi * inv_r3 * dz
      end if
      return
    end if

    do child_k = 1_i32, self%child_count(node_idx)
      call traverse_node(self, mesh, self%child_idx(child_k, node_idx), rx, ry, rz, soft2, ex, ey, ez)
    end do
  end procedure traverse_node

  !> ノードサイズ・距離・電荷打ち消し度合いから近似採用可否を判定する。
  module procedure accept_node
    real(dp) :: dx, dy, dz, dist, dist2, radius

    dx = rx - self%node_center(1, node_idx)
    dy = ry - self%node_center(2, node_idx)
    dz = rz - self%node_center(3, node_idx)
    dist2 = dx * dx + dy * dy + dz * dz

    if (dist2 <= 0.0d0) then
      accept_it = .false.
      return
    end if

    radius = self%node_radius(node_idx)
    dist = sqrt(dist2)
    if (dist <= radius) then
      accept_it = .false.
      return
    end if
    if (radius >= self%theta * (dist - radius)) then
      accept_it = .false.
      return
    end if

    if (self%node_abs_q(node_idx) <= 0.0d0) then
      accept_it = .true.
      return
    end if

    accept_it = abs(self%node_q(node_idx)) >= charge_cancellation_tol * self%node_abs_q(node_idx)
  end procedure accept_node

end submodule bem_field_solver_eval

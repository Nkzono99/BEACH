!> `bem_field_solver` の電場評価と木走査ロジックを実装する submodule。
submodule(bem_field_solver) bem_field_solver_eval
  use bem_coulomb_fmm_core, only: eval_point, eval_potential_point, eval_potential_points
  use bem_panel_geometry, only: panel_geometry_type
  use bem_panel_kernel, only: panel_potential_field, panel_side_principal_value
  implicit none
contains

  !> 観測点の電場を direct / treecode / fmm で評価して返す。
  module procedure eval_e_field_solver
  real(dp) :: rx, ry, rz, ex, ey, ez, r_scaled(3)

  if (trim(self%mode) == 'fmm') then
    if (mesh%nelem <= 0_i32) then
      e = 0.0d0
      return
    end if
    if (.not. self%fmm_core_ready) then
      error stop 'FMM core is not ready. Call solver%init/refresh before eval_e.'
    end if
    r_scaled = (r - self%field_origin)*self%field_inv_length_scale
    call eval_point(self%fmm_core_plan, self%fmm_core_state, r_scaled, e)
    e = self%field_output_scale*e
    return
  end if

  if (trim(self%mode) /= 'treecode' .or. .not. self%tree_ready) then
    call electric_field_at_panel_mesh(mesh, r, e)
    return
  end if

  rx = r(1)
  ry = r(2)
  rz = r(3)
  ex = 0.0d0
  ey = 0.0d0
  ez = 0.0d0

  call traverse_node(self, mesh, 1_i32, rx, ry, rz, ex, ey, ez)

  e(1) = self%field_output_scale*ex
  e(2) = self%field_output_scale*ey
  e(3) = self%field_output_scale*ez
  end procedure eval_e_field_solver

  !> 観測点の電位を direct / treecode / FMM で評価して返す。
  module procedure eval_potential_field_solver
  real(dp) :: r_scaled(3), phi_sum

  phi = 0.0d0
  if (mesh%nelem <= 0_i32) return

  if (trim(self%mode) == 'fmm') then
    if (.not. self%fmm_core_ready) then
      error stop 'FMM core is not ready. Call solver%init/refresh before eval_potential.'
    end if
    r_scaled = (r - self%field_origin)*self%field_inv_length_scale
    call eval_potential_point(self%fmm_core_plan, self%fmm_core_state, r_scaled, phi)
    phi = self%potential_output_scale*phi
    return
  end if

  if (trim(self%mode) == 'treecode' .and. self%tree_ready) then
    phi_sum = 0.0_dp
    call traverse_potential_node(self, mesh, 1_i32, r(1), r(2), r(3), phi_sum)
    phi = self%potential_output_scale*phi_sum
    return
  end if

  call electric_potential_at_panel_mesh(mesh, r, phi)
  end procedure eval_potential_field_solver

  !> ノードを再帰走査し、葉は direct 総和・遠方は monopole 近似で加算する。
  module procedure traverse_node
  integer(i32) :: child_k, p, idx, p_end
  real(dp) :: dx, dy, dz, r2, inv_r3, qi, min_dist2
  real(dp) :: target(3), source_potential, source_field(3)
  type(panel_geometry_type) :: geometry

  min_dist2 = tiny(1.0d0)
  if (self%child_count(node_idx) <= 0_i32) then
    p_end = self%node_start(node_idx) + self%node_count(node_idx) - 1_i32
    target = [rx, ry, rz]
    do p = self%node_start(node_idx), p_end
      idx = self%elem_order(p)
      if (mesh%elem_vacuum_sign(idx) /= 1_i32 .and. mesh%elem_vacuum_sign(idx) /= -1_i32) then
        error stop 'triangle_p0 tree field evaluation requires a resolved vacuum side for every element.'
      end if
      call geometry_from_mesh(mesh, idx, geometry)
      call panel_potential_field( &
        geometry, mesh%q_elem(idx), target, mesh%elem_vacuum_sign(idx), source_potential, source_field &
        )
      ex = ex + source_field(1)/self%field_output_scale
      ey = ey + source_field(2)/self%field_output_scale
      ez = ez + source_field(3)/self%field_output_scale
    end do
    return
  end if

  if (accept_node(self, node_idx, rx, ry, rz)) then
    qi = self%node_q(node_idx)
    if (abs(qi) > 0.0d0) then
      dx = (rx - self%node_charge_center(1, node_idx))*self%field_inv_length_scale
      dy = (ry - self%node_charge_center(2, node_idx))*self%field_inv_length_scale
      dz = (rz - self%node_charge_center(3, node_idx))*self%field_inv_length_scale
      r2 = dx*dx + dy*dy + dz*dz
      if (r2 <= min_dist2) return
      inv_r3 = 1.0d0/(sqrt(r2)*r2)
      ex = ex + qi*inv_r3*dx
      ey = ey + qi*inv_r3*dy
      ez = ez + qi*inv_r3*dz
    end if
    return
  end if

  do child_k = 1_i32, self%child_count(node_idx)
    call traverse_node(self, mesh, self%child_idx(child_k, node_idx), rx, ry, rz, ex, ey, ez)
  end do
  end procedure traverse_node

  !> ノードを再帰走査し、葉は選択 source kernel の direct 和、遠方は monopole 電位を加算する。
  module procedure traverse_potential_node
  integer(i32) :: child_k, p, idx, p_end
  real(dp) :: dx, dy, dz, r2, inv_r, qi, min_dist2
  real(dp) :: target(3), source_potential, source_field(3)
  type(panel_geometry_type) :: geometry

  min_dist2 = tiny(1.0_dp)
  if (self%child_count(node_idx) <= 0_i32) then
    p_end = self%node_start(node_idx) + self%node_count(node_idx) - 1_i32
    target = [rx, ry, rz]
    do p = self%node_start(node_idx), p_end
      idx = self%elem_order(p)
      call geometry_from_mesh(mesh, idx, geometry)
      call panel_potential_field( &
        geometry, mesh%q_elem(idx), target, panel_side_principal_value, source_potential, source_field &
        )
      phi_sum = phi_sum + source_potential/self%potential_output_scale
    end do
    return
  end if

  if (accept_node(self, node_idx, rx, ry, rz)) then
    qi = self%node_q(node_idx)
    if (abs(qi) > 0.0_dp) then
      dx = (rx - self%node_charge_center(1, node_idx))*self%field_inv_length_scale
      dy = (ry - self%node_charge_center(2, node_idx))*self%field_inv_length_scale
      dz = (rz - self%node_charge_center(3, node_idx))*self%field_inv_length_scale
      r2 = dx*dx + dy*dy + dz*dz
      if (r2 <= min_dist2) return
      phi_sum = phi_sum + qi/sqrt(r2)
    end if
    return
  end if

  do child_k = 1_i32, self%child_count(node_idx)
    call traverse_potential_node( &
      self, mesh, self%child_idx(child_k, node_idx), rx, ry, rz, phi_sum &
      )
  end do
  end procedure traverse_potential_node

  !> ノードサイズ・距離・電荷符号の一貫性から近似採用可否を判定する。
  module procedure accept_node
  real(dp) :: dx, dy, dz, dist, dist2, radius, charge_scale, charge_gap

  dx = rx - self%node_center(1, node_idx)
  dy = ry - self%node_center(2, node_idx)
  dz = rz - self%node_center(3, node_idx)
  dist2 = dx*dx + dy*dy + dz*dz

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
  if (radius >= self%theta*(dist - radius)) then
    accept_it = .false.
    return
  end if

  charge_scale = max(self%node_abs_q(node_idx), abs(self%node_q(node_idx)))
  charge_gap = abs(self%node_abs_q(node_idx) - abs(self%node_q(node_idx)))
  accept_it = charge_gap <= 64.0d0*epsilon(1.0d0)*charge_scale
  end procedure accept_node

  !> メッシュ重心での電位を計算する。FMM/treecode/direct を自動切替する。
  module procedure compute_mesh_potential_field_solver
  if (size(potential_v) /= mesh%nelem) error stop 'mesh potential output array size mismatch.'
  potential_v = 0.0d0
  if (mesh%nelem <= 0) return

  if (self%fmm_use_core .and. self%fmm_core_ready) then
    call compute_mesh_potential_fmm(self, mesh, potential_v)
  else if (trim(self%mode) == 'treecode' .and. self%tree_ready) then
    call compute_mesh_potential_tree(self, mesh, potential_v)
  else
    call compute_mesh_potential_direct(self, mesh, potential_v)
  end if
  end procedure compute_mesh_potential_field_solver

  !> 構築済み source tree を各要素重心から走査してメッシュ電位を計算する。
  subroutine compute_mesh_potential_tree(self, mesh, potential_v)
    class(field_solver_type), intent(in) :: self
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(out) :: potential_v(:)

    integer(i32) :: i
    real(dp) :: phi_sum

    !$omp parallel do default(none) schedule(static) &
    !$omp   shared(self,mesh,potential_v) private(i,phi_sum)
    do i = 1_i32, mesh%nelem
      phi_sum = 0.0_dp
      call traverse_potential_node( &
        self, mesh, 1_i32, mesh%center_x(i), mesh%center_y(i), mesh%center_z(i), phi_sum &
        )
      potential_v(i) = self%potential_output_scale*phi_sum
    end do
    !$omp end parallel do
  end subroutine compute_mesh_potential_tree

  !> 構築済み Coulomb FMM core を使って要素重心電位を計算する。
  subroutine compute_mesh_potential_fmm(self, mesh, potential_v)
    class(field_solver_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(out) :: potential_v(:)

    real(dp), allocatable :: centers_scaled(:, :)
    integer(i32) :: i

    if (.not. self%fmm_core_plan%built .or. .not. self%fmm_core_state%ready) &
      error stop 'FMM core is not ready for mesh potential output.'
    if (self%fmm_core_plan%nsrc /= mesh%nelem) &
      error stop 'FMM source count does not match mesh size.'

    allocate (centers_scaled(3, mesh%nelem))
    do i = 1, mesh%nelem
      centers_scaled(:, i) = (mesh%centers(:, i) - self%field_origin)*self%field_inv_length_scale
    end do
    call eval_potential_points(self%fmm_core_plan, self%fmm_core_state, centers_scaled, potential_v)

    potential_v = self%potential_output_scale*potential_v
    deallocate (centers_scaled)
  end subroutine compute_mesh_potential_fmm

  !> O(N^2) 直接総和でメッシュ重心電位を計算する。
  subroutine compute_mesh_potential_direct(self, mesh, potential_v)
    class(field_solver_type), intent(in) :: self
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(out) :: potential_v(:)

    integer(i32) :: i

    do i = 1_i32, mesh%nelem
      call electric_potential_at_panel_mesh(mesh, mesh%centers(:, i), potential_v(i))
    end do
  end subroutine compute_mesh_potential_direct

  subroutine electric_field_at_panel_mesh(mesh, target, field)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: target(3)
    real(dp), intent(out) :: field(3)
    type(panel_geometry_type) :: geometry
    real(dp) :: source_potential, source_field(3)
    integer(i32) :: element

    field = 0.0_dp
    do element = 1, mesh%nelem
      if (mesh%elem_vacuum_sign(element) /= 1_i32 .and. mesh%elem_vacuum_sign(element) /= -1_i32) then
        error stop 'triangle_p0 field evaluation requires a resolved vacuum side for every element.'
      end if
      call geometry_from_mesh(mesh, element, geometry)
      call panel_potential_field( &
        geometry, mesh%q_elem(element), target, mesh%elem_vacuum_sign(element), source_potential, source_field &
        )
      field = field + source_field
    end do
  end subroutine electric_field_at_panel_mesh

  subroutine electric_potential_at_panel_mesh(mesh, target, potential)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: target(3)
    real(dp), intent(out) :: potential
    type(panel_geometry_type) :: geometry
    real(dp) :: source_potential, source_field(3)
    integer(i32) :: element

    potential = 0.0_dp
    do element = 1, mesh%nelem
      if (mesh%elem_vacuum_sign(element) /= 1_i32 .and. mesh%elem_vacuum_sign(element) /= -1_i32) then
        error stop 'triangle_p0 potential evaluation requires a resolved vacuum side for every element.'
      end if
      call geometry_from_mesh(mesh, element, geometry)
      call panel_potential_field( &
        geometry, mesh%q_elem(element), target, panel_side_principal_value, source_potential, source_field &
        )
      potential = potential + source_potential
    end do
  end subroutine electric_potential_at_panel_mesh

  pure subroutine geometry_from_mesh(mesh, element, geometry)
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: element
    type(panel_geometry_type), intent(out) :: geometry

    geometry = panel_geometry_type()
    geometry%vertex(:, 1) = mesh%v0(:, element)
    geometry%vertex(:, 2) = mesh%v1(:, element)
    geometry%vertex(:, 3) = mesh%v2(:, element)
    geometry%centroid = mesh%centers(:, element)
    geometry%normal = mesh%normals(:, element)
    geometry%area = mesh%panel_area(element)
    geometry%moment0 = mesh%panel_area(element)
    geometry%moment1 = mesh%panel_moment1(:, element)
    geometry%moment2 = mesh%panel_moment2(:, :, element)
    geometry%edge_length = mesh%panel_edge_length(:, element)
    geometry%edge_outward = mesh%panel_edge_outward(:, :, element)
  end subroutine geometry_from_mesh

end submodule bem_field_solver_eval

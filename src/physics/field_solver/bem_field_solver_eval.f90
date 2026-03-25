!> `bem_field_solver` の電場評価と木走査ロジックを実装する submodule。
submodule(bem_field_solver) bem_field_solver_eval
  use bem_coulomb_fmm_core, only: eval_point, eval_potential_points
  use bem_coulomb_fmm_types, only: fmm_plan_type, fmm_state_type, &
                                   reset_fmm_plan, initialize_fmm_state, reset_fmm_state
  use bem_coulomb_fmm_periodic, only: use_periodic2_m2l_root_oracle
  use bem_coulomb_fmm_periodic_ewald, only: precompute_periodic2_ewald_data, &
                                            add_periodic2_exact_ewald_potential_correction_all_sources
  use bem_string_utils, only: lower_ascii
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
    e = k_coulomb*e
    return
  end if

  if (trim(self%mode) /= 'treecode' .or. .not. self%tree_ready) then
    call electric_field_at(mesh, r, self%softening, e)
    return
  end if

  rx = r(1)
  ry = r(2)
  rz = r(3)
  soft2 = self%softening*self%softening
  ex = 0.0d0
  ey = 0.0d0
  ez = 0.0d0

  call traverse_node(self, mesh, 1_i32, rx, ry, rz, soft2, ex, ey, ez)

  e(1) = k_coulomb*ex
  e(2) = k_coulomb*ey
  e(3) = k_coulomb*ez
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
      r2 = dx*dx + dy*dy + dz*dz + soft2
      inv_r3 = 1.0d0/(sqrt(r2)*r2)
      qi = mesh%q_elem(idx)
      ex = ex + qi*inv_r3*dx
      ey = ey + qi*inv_r3*dy
      ez = ez + qi*inv_r3*dz
    end do
    return
  end if

  if (accept_node(self, node_idx, rx, ry, rz)) then
    qi = self%node_q(node_idx)
    if (abs(qi) > 0.0d0) then
      dx = rx - self%node_charge_center(1, node_idx)
      dy = ry - self%node_charge_center(2, node_idx)
      dz = rz - self%node_charge_center(3, node_idx)
      r2 = dx*dx + dy*dy + dz*dz + soft2
      inv_r3 = 1.0d0/(sqrt(r2)*r2)
      ex = ex + qi*inv_r3*dx
      ey = ey + qi*inv_r3*dy
      ez = ez + qi*inv_r3*dz
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

  if (self%node_abs_q(node_idx) <= 0.0d0) then
    accept_it = .true.
    return
  end if

  accept_it = abs(self%node_q(node_idx)) >= charge_cancellation_tol*self%node_abs_q(node_idx)
  end procedure accept_node

  !> メッシュ重心での電位を計算する。FMM/direct を自動切替する。
  module procedure compute_mesh_potential_field_solver
  if (size(potential_v) /= mesh%nelem) error stop 'mesh potential output array size mismatch.'
  potential_v = 0.0d0
  if (mesh%nelem <= 0) return

  if (self%fmm_use_core .and. self%fmm_core_ready) then
    call compute_mesh_potential_fmm(self, mesh, potential_v)
  else
    call compute_mesh_potential_direct(self, mesh, sim, potential_v)
  end if
  end procedure compute_mesh_potential_field_solver

  !> 構築済み Coulomb FMM core を使って要素重心電位を計算する。
  subroutine compute_mesh_potential_fmm(self, mesh, potential_v)
    class(field_solver_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(out) :: potential_v(:)

    real(dp), parameter :: pi_dp = acos(-1.0d0)
    real(dp) :: phi_anchor_exact, phi_offset, self_coeff, softening
    integer(i32) :: i

    if (.not. self%fmm_core_plan%built .or. .not. self%fmm_core_state%ready) &
      error stop 'FMM core is not ready for mesh potential output.'
    if (self%fmm_core_plan%nsrc /= mesh%nelem) &
      error stop 'FMM source count does not match mesh size.'

    softening = self%fmm_core_plan%options%softening
    if (softening < 0.0d0) error stop 'mesh potential output requires non-negative FMM softening.'

    call eval_potential_points(self%fmm_core_plan, self%fmm_core_state, mesh%centers, potential_v)

    if (use_periodic2_m2l_root_oracle(self%fmm_core_plan)) then
      call compute_exact_potential_point(self%fmm_core_plan, self%fmm_core_state, mesh%centers(:, 1), phi_anchor_exact)
      phi_offset = phi_anchor_exact - potential_v(1)
      potential_v = potential_v + phi_offset
    end if

    if (softening <= 0.0d0) then
      do i = 1, mesh%nelem
        self_coeff = 2.0d0*sqrt(pi_dp)/max(mesh%h_elem(i), sqrt(tiny(1.0d0)))
        potential_v(i) = potential_v(i) + self_coeff*mesh%q_elem(i)
      end do
    end if

    potential_v = k_coulomb*potential_v
  end subroutine compute_mesh_potential_fmm

  !> FMM source/state を使って 1 点の exact 電位を O(N) で再評価する。
  subroutine compute_exact_potential_point(plan, state, target_in, phi)
    type(fmm_plan_type), intent(in) :: plan
    type(fmm_state_type), intent(in) :: state
    real(dp), intent(in) :: target_in(3)
    real(dp), intent(out) :: phi

    integer(i32) :: j, img_i, img_j, axis1, axis2, nimg
    real(dp) :: target(3), shifted(3), dx, dy, dz, r2, soft2, min_dist2

    phi = 0.0d0
    if (.not. plan%built .or. .not. state%ready) return

    soft2 = plan%options%softening*plan%options%softening
    min_dist2 = tiny(1.0d0)
    target = target_in

    if (.not. plan%options%use_periodic2) then
      do j = 1, plan%nsrc
        if (abs(state%src_q(j)) <= tiny(1.0d0)) cycle
        dx = target(1) - plan%src_pos(1, j)
        dy = target(2) - plan%src_pos(2, j)
        dz = target(3) - plan%src_pos(3, j)
        r2 = dx*dx + dy*dy + dz*dz + soft2
        if (r2 <= min_dist2) cycle
        phi = phi + state%src_q(j)/sqrt(r2)
      end do
      return
    end if

    axis1 = plan%options%periodic_axes(1)
    axis2 = plan%options%periodic_axes(2)
    nimg = max(0_i32, plan%options%periodic_image_layers)
    target(axis1) = wrap_periodic(target(axis1), plan%options%target_box_min(axis1), plan%options%periodic_len(1))
    target(axis2) = wrap_periodic(target(axis2), plan%options%target_box_min(axis2), plan%options%periodic_len(2))

    do j = 1, plan%nsrc
      if (abs(state%src_q(j)) <= tiny(1.0d0)) cycle
      dx = target(1) - plan%src_pos(1, j)
      dy = target(2) - plan%src_pos(2, j)
      dz = target(3) - plan%src_pos(3, j)
      r2 = dx*dx + dy*dy + dz*dz + soft2
      if (r2 > min_dist2) phi = phi + state%src_q(j)/sqrt(r2)
    end do

    do img_i = -nimg, nimg
      do img_j = -nimg, nimg
        if (img_i == 0_i32 .and. img_j == 0_i32) cycle
        do j = 1, plan%nsrc
          if (abs(state%src_q(j)) <= tiny(1.0d0)) cycle
          shifted = plan%src_pos(:, j)
          shifted(axis1) = shifted(axis1) + real(img_i, dp)*plan%options%periodic_len(1)
          shifted(axis2) = shifted(axis2) + real(img_j, dp)*plan%options%periodic_len(2)
          dx = target(1) - shifted(1)
          dy = target(2) - shifted(2)
          dz = target(3) - shifted(3)
          r2 = dx*dx + dy*dy + dz*dz + soft2
          if (r2 <= min_dist2) cycle
          phi = phi + state%src_q(j)/sqrt(r2)
        end do
      end do
    end do

    if (use_periodic2_m2l_root_oracle(plan)) then
      call add_periodic2_exact_ewald_potential_correction_all_sources(plan, state, target, phi)
    end if
  end subroutine compute_exact_potential_point

  !> O(N^2) 直接総和でメッシュ重心電位を計算する。
  subroutine compute_mesh_potential_direct(self, mesh, sim, potential_v)
    class(field_solver_type), intent(in) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    real(dp), intent(out) :: potential_v(:)

    real(dp), parameter :: pi_dp = acos(-1.0d0)
    integer(i32) :: i, j, ix, iy, axis1, axis2, image_layers
    integer(i32) :: periodic_axes(2), n_periodic, axis
    real(dp) :: softening, soft2, min_dist2, self_coeff
    real(dp) :: periodic_len(2), periodic_origins(2), span
    real(dp) :: target(3), shifted(3), dx, dy, dz, r2
    logical :: use_periodic2, use_exact_ewald_residual
    type(fmm_plan_type) :: ewald_plan
    type(fmm_state_type) :: ewald_state
    character(len=16) :: field_bc_mode

    softening = sim%softening
    if (softening < 0.0d0) error stop 'sim.softening must be >= 0 for mesh potential output.'
    soft2 = softening*softening
    min_dist2 = tiny(1.0d0)
    potential_v = 0.0d0

    ! Resolve periodic2 configuration
    use_periodic2 = .false.
    periodic_axes = 0_i32
    periodic_len = 0.0d0
    periodic_origins = 0.0d0
    image_layers = max(0_i32, sim%field_periodic_image_layers)

    field_bc_mode = trim(lower_ascii(sim%field_bc_mode))
    if (field_bc_mode == 'periodic2') then
      if (.not. sim%use_box) error stop 'mesh potential output for periodic2 requires sim.use_box=true.'
      n_periodic = 0_i32
      do axis = 1_i32, 3_i32
        if ((sim%bc_low(axis) == bc_periodic) .neqv. (sim%bc_high(axis) == bc_periodic)) then
          error stop 'mesh potential output requires bc_low(axis)=bc_high(axis)=periodic for periodic axes.'
        end if
        if (sim%bc_low(axis) == bc_periodic) then
          n_periodic = n_periodic + 1_i32
          if (n_periodic <= 2_i32) periodic_axes(n_periodic) = axis
        end if
      end do
      if (n_periodic /= 2_i32) then
        error stop 'mesh potential output with sim.field_bc_mode="periodic2" requires exactly two periodic axes.'
      end if
      do axis = 1_i32, 2_i32
        span = sim%box_max(periodic_axes(axis)) - sim%box_min(periodic_axes(axis))
        if (span <= 0.0d0) error stop 'mesh potential output requires positive box length on periodic axes.'
        periodic_len(axis) = span
        periodic_origins(axis) = sim%box_min(periodic_axes(axis))
      end do
      use_periodic2 = .true.
    end if

    axis1 = periodic_axes(1)
    axis2 = periodic_axes(2)
    call initialize_fmm_state(ewald_state)
    call init_ewald_context( &
      mesh, sim, use_periodic2, periodic_axes, periodic_len, image_layers, ewald_plan, ewald_state, use_exact_ewald_residual &
      )

    do i = 1, mesh%nelem
      if (softening > 0.0d0) then
        self_coeff = 1.0d0/softening
      else
        self_coeff = 2.0d0*sqrt(pi_dp)/max(mesh%h_elem(i), sqrt(tiny(1.0d0)))
      end if
      potential_v(i) = self_coeff*mesh%q_elem(i)
    end do

    if (.not. use_periodic2) then
      !$omp parallel do default(none) schedule(static) &
      !$omp   shared(mesh, soft2, min_dist2, potential_v) private(i, j, dx, dy, dz, r2)
      do i = 1, mesh%nelem
        do j = 1, mesh%nelem
          if (j == i) cycle
          dx = mesh%center_x(i) - mesh%center_x(j)
          dy = mesh%center_y(i) - mesh%center_y(j)
          dz = mesh%center_z(i) - mesh%center_z(j)
          r2 = dx*dx + dy*dy + dz*dz + soft2
          potential_v(i) = potential_v(i) + mesh%q_elem(j)/sqrt(max(r2, min_dist2))
        end do
      end do
      !$omp end parallel do
      potential_v = k_coulomb*potential_v
      return
    end if

    !$omp parallel do default(none) schedule(static) &
    !$omp   shared(mesh, soft2, min_dist2, potential_v, axis1, axis2, &
    !$omp          periodic_origins, periodic_len, image_layers, &
    !$omp          use_exact_ewald_residual, ewald_plan, ewald_state) &
    !$omp   private(i, j, ix, iy, target, shifted, dx, dy, dz, r2)
    do i = 1, mesh%nelem
      target = mesh%centers(:, i)
      target(axis1) = wrap_periodic(target(axis1), periodic_origins(1), periodic_len(1))
      target(axis2) = wrap_periodic(target(axis2), periodic_origins(2), periodic_len(2))

      do j = 1, mesh%nelem
        if (j == i) cycle
        dx = target(1) - mesh%centers(1, j)
        dy = target(2) - mesh%centers(2, j)
        dz = target(3) - mesh%centers(3, j)
        r2 = dx*dx + dy*dy + dz*dz + soft2
        potential_v(i) = potential_v(i) + mesh%q_elem(j)/sqrt(max(r2, min_dist2))
      end do

      do ix = -image_layers, image_layers
        do iy = -image_layers, image_layers
          if (ix == 0_i32 .and. iy == 0_i32) cycle
          do j = 1, mesh%nelem
            shifted = mesh%centers(:, j)
            shifted(axis1) = shifted(axis1) + real(ix, dp)*periodic_len(1)
            shifted(axis2) = shifted(axis2) + real(iy, dp)*periodic_len(2)
            dx = target(1) - shifted(1)
            dy = target(2) - shifted(2)
            dz = target(3) - shifted(3)
            r2 = dx*dx + dy*dy + dz*dz + soft2
            potential_v(i) = potential_v(i) + mesh%q_elem(j)/sqrt(max(r2, min_dist2))
          end do
        end do
      end do

      if (use_exact_ewald_residual) then
        call add_periodic2_exact_ewald_potential_correction_all_sources(ewald_plan, ewald_state, target, potential_v(i))
      end if
    end do
    !$omp end parallel do

    if (use_exact_ewald_residual) then
      call reset_fmm_state(ewald_state)
      call reset_fmm_plan(ewald_plan)
    end if
    potential_v = k_coulomb*potential_v
  end subroutine compute_mesh_potential_direct

  !> exact Ewald residual 用の最小 periodic2 context を作る。
  subroutine init_ewald_context( &
    mesh, sim, use_periodic2, periodic_axes, periodic_len, image_layers, plan, state, enabled &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    logical, intent(in) :: use_periodic2
    integer(i32), intent(in) :: periodic_axes(2), image_layers
    real(dp), intent(in) :: periodic_len(2)
    type(fmm_plan_type), intent(inout) :: plan
    type(fmm_state_type), intent(inout) :: state
    logical, intent(out) :: enabled

    character(len=16) :: far_correction
    integer(i32) :: ewald_layers

    enabled = .false.
    call reset_fmm_plan(plan)
    call reset_fmm_state(state)
    if (.not. use_periodic2) return

    far_correction = trim(lower_ascii(sim%field_periodic_far_correction))
    ewald_layers = max(0_i32, sim%field_periodic_ewald_layers)
    select case (trim(far_correction))
    case ('auto')
      far_correction = 'm2l_root_oracle'
      ewald_layers = max(1_i32, ewald_layers)
    case ('none')
      return
    case ('m2l_root_oracle')
      ewald_layers = max(1_i32, ewald_layers)
    case default
      return
    end select

    plan%options%use_periodic2 = .true.
    plan%options%periodic_far_correction = far_correction
    plan%options%periodic_axes = periodic_axes
    plan%options%periodic_len = periodic_len
    plan%options%periodic_image_layers = image_layers
    plan%options%periodic_ewald_alpha = sim%field_periodic_ewald_alpha
    plan%options%periodic_ewald_layers = ewald_layers
    plan%options%softening = sim%softening
    plan%options%target_box_min = sim%box_min
    plan%options%target_box_max = sim%box_max
    plan%nsrc = mesh%nelem
    allocate (plan%src_pos(3, mesh%nelem))
    plan%src_pos = mesh%centers
    allocate (state%src_q(mesh%nelem))
    state%src_q = mesh%q_elem
    state%ready = .true.

    call precompute_periodic2_ewald_data(plan)
    enabled = plan%periodic_ewald%ready
    if (.not. enabled) then
      call reset_fmm_state(state)
      call reset_fmm_plan(plan)
    end if
  end subroutine init_ewald_context

  !> 周期軸上の座標を [origin, origin + length) に折り返す。
  pure real(dp) function wrap_periodic(x, origin, periodic_len) result(wrapped)
    real(dp), intent(in) :: x, origin, periodic_len

    wrapped = origin + modulo(x - origin, periodic_len)
  end function wrap_periodic

end submodule bem_field_solver_eval

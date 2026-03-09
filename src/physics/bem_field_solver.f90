!> 粒子位置での電場評価を direct / treecode で切り替える場ソルバ。
module bem_field_solver
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_types, only: mesh_type, sim_config
  use bem_field, only: electric_field_at
  implicit none
  private

  real(dp), parameter :: charge_cancellation_tol = 1.0d-10

  type :: field_solver_type
    character(len=16) :: mode = 'direct'
    real(dp) :: softening = 1.0d-6
    real(dp) :: theta = 0.5d0
    integer(i32) :: leaf_max = 16_i32
    integer(i32) :: min_nelem = 256_i32
    logical :: tree_ready = .false.
    integer(i32) :: nelem = 0_i32
    integer(i32) :: max_node = 0_i32
    integer(i32) :: nnode = 0_i32
    integer(i32), allocatable :: elem_order(:)
    integer(i32), allocatable :: node_start(:), node_count(:)
    integer(i32), allocatable :: child_count(:), child_idx(:, :)
    real(dp), allocatable :: node_center(:, :)
    real(dp), allocatable :: node_half_size(:, :)
    real(dp), allocatable :: node_radius(:)
    real(dp), allocatable :: node_q(:), node_abs_q(:)
    real(dp), allocatable :: node_qx(:), node_qy(:), node_qz(:)
    real(dp), allocatable :: node_charge_center(:, :)
  contains
    procedure :: init => init_field_solver
    procedure :: refresh => refresh_field_solver
    procedure :: eval_e => eval_e_field_solver
  end type field_solver_type

  public :: field_solver_type

contains

  !> 設定とメッシュから電場ソルバを初期化する。
  subroutine init_field_solver(self, mesh, sim)
    class(field_solver_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    character(len=16) :: requested_mode

    self%softening = sim%softening
    self%theta = sim%tree_theta
    self%leaf_max = sim%tree_leaf_max
    self%min_nelem = sim%tree_min_nelem
    self%nelem = mesh%nelem

    requested_mode = lower_ascii(trim(sim%field_solver))
    select case (trim(requested_mode))
    case ('direct')
      self%mode = 'direct'
    case ('treecode')
      self%mode = 'treecode'
    case ('auto')
      if (mesh%nelem >= self%min_nelem) then
        self%mode = 'treecode'
      else
        self%mode = 'direct'
      end if
    case default
      error stop 'Unknown sim.field_solver in field solver init.'
    end select

    if (trim(self%mode) == 'treecode' .and. mesh%nelem > 0_i32) then
      call build_tree_topology(self, mesh)
      call refresh_field_solver(self, mesh)
    else
      call reset_tree_storage(self)
      self%tree_ready = .false.
    end if
  end subroutine init_field_solver

  !> 現在の要素電荷から treecode モーメントを再計算する。
  subroutine refresh_field_solver(self, mesh)
    class(field_solver_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh

    integer(i32) :: node_idx, child_k, child_node
    integer(i32) :: p, idx, p_end
    real(dp) :: q, abs_q, qx, qy, qz, qi

    if (trim(self%mode) /= 'treecode') return
    if (mesh%nelem <= 0_i32) return

    if (.not. self%tree_ready .or. self%nelem /= mesh%nelem) then
      self%nelem = mesh%nelem
      call build_tree_topology(self, mesh)
    end if

    do node_idx = self%nnode, 1_i32, -1_i32
      if (self%child_count(node_idx) <= 0_i32) then
        q = 0.0d0
        abs_q = 0.0d0
        qx = 0.0d0
        qy = 0.0d0
        qz = 0.0d0
        p_end = self%node_start(node_idx) + self%node_count(node_idx) - 1_i32
        do p = self%node_start(node_idx), p_end
          idx = self%elem_order(p)
          qi = mesh%q_elem(idx)
          q = q + qi
          abs_q = abs_q + abs(qi)
          qx = qx + qi * mesh%center_x(idx)
          qy = qy + qi * mesh%center_y(idx)
          qz = qz + qi * mesh%center_z(idx)
        end do
      else
        q = 0.0d0
        abs_q = 0.0d0
        qx = 0.0d0
        qy = 0.0d0
        qz = 0.0d0
        do child_k = 1_i32, self%child_count(node_idx)
          child_node = self%child_idx(child_k, node_idx)
          q = q + self%node_q(child_node)
          abs_q = abs_q + self%node_abs_q(child_node)
          qx = qx + self%node_qx(child_node)
          qy = qy + self%node_qy(child_node)
          qz = qz + self%node_qz(child_node)
        end do
      end if

      self%node_q(node_idx) = q
      self%node_abs_q(node_idx) = abs_q
      self%node_qx(node_idx) = qx
      self%node_qy(node_idx) = qy
      self%node_qz(node_idx) = qz

      if (abs(q) > tiny(1.0d0)) then
        self%node_charge_center(1, node_idx) = qx / q
        self%node_charge_center(2, node_idx) = qy / q
        self%node_charge_center(3, node_idx) = qz / q
      else
        self%node_charge_center(:, node_idx) = self%node_center(:, node_idx)
      end if
    end do
  end subroutine refresh_field_solver

  !> 観測点 `r` の電場を設定されたソルバで評価する。
  subroutine eval_e_field_solver(self, mesh, r, e)
    class(field_solver_type), intent(in) :: self
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: r(3)
    real(dp), intent(out) :: e(3)

    real(dp) :: rx, ry, rz, soft2, ex, ey, ez

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
  end subroutine eval_e_field_solver

  !> メッシュ重心を octree 分割して treecode トポロジを構築する。
  subroutine build_tree_topology(self, mesh)
    class(field_solver_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    integer(i32) :: i, max_node_guess

    if (mesh%nelem <= 0_i32) then
      call reset_tree_storage(self)
      self%tree_ready = .false.
      self%nelem = 0_i32
      return
    end if

    max_node_guess = max(1_i32, 2_i32 * mesh%nelem)
    call ensure_tree_capacity(self, max_node_guess)

    if (allocated(self%elem_order)) then
      if (size(self%elem_order) /= mesh%nelem) deallocate (self%elem_order)
    end if
    if (.not. allocated(self%elem_order)) allocate (self%elem_order(mesh%nelem))

    do i = 1_i32, mesh%nelem
      self%elem_order(i) = i
    end do

    self%nnode = 1_i32
    call build_node(self, mesh, 1_i32, 1_i32, mesh%nelem)
    self%nelem = mesh%nelem
    self%tree_ready = .true.
  end subroutine build_tree_topology

  !> 要素添字区間を1ノードとして登録し、必要なら8分割で子ノードを作る。
  recursive subroutine build_node(self, mesh, node_idx, start_idx, end_idx)
    class(field_solver_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: node_idx, start_idx, end_idx

    integer(i32) :: count, p, idx, oct
    integer(i32) :: child_k, child_node, child_start, child_end
    integer(i32), allocatable :: counts(:), offsets(:), cursor(:), work(:)
    real(dp) :: bb_min(3), bb_max(3), span(3), center(3)
    real(dp) :: split_eps

    count = end_idx - start_idx + 1_i32
    self%node_start(node_idx) = start_idx
    self%node_count(node_idx) = count
    self%child_count(node_idx) = 0_i32
    self%child_idx(:, node_idx) = 0_i32

    idx = self%elem_order(start_idx)
    bb_min = [mesh%center_x(idx), mesh%center_y(idx), mesh%center_z(idx)]
    bb_max = bb_min
    do p = start_idx + 1_i32, end_idx
      idx = self%elem_order(p)
      bb_min(1) = min(bb_min(1), mesh%center_x(idx))
      bb_min(2) = min(bb_min(2), mesh%center_y(idx))
      bb_min(3) = min(bb_min(3), mesh%center_z(idx))
      bb_max(1) = max(bb_max(1), mesh%center_x(idx))
      bb_max(2) = max(bb_max(2), mesh%center_y(idx))
      bb_max(3) = max(bb_max(3), mesh%center_z(idx))
    end do

    span = bb_max - bb_min
    center = 0.5d0 * (bb_max + bb_min)
    self%node_center(:, node_idx) = center
    self%node_half_size(:, node_idx) = 0.5d0 * span
    self%node_radius(node_idx) = sqrt(sum(self%node_half_size(:, node_idx) * self%node_half_size(:, node_idx)))

    if (count <= self%leaf_max) return

    split_eps = 1.0d-12 * max(1.0d0, maxval(abs(center)))
    if (maxval(span) <= split_eps) return

    allocate (counts(8), offsets(8), cursor(8), work(count))
    counts = 0_i32
    do p = start_idx, end_idx
      idx = self%elem_order(p)
      oct = octant_index(mesh%center_x(idx), mesh%center_y(idx), mesh%center_z(idx), center)
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
      idx = self%elem_order(p)
      oct = octant_index(mesh%center_x(idx), mesh%center_y(idx), mesh%center_z(idx), center)
      work(cursor(oct)) = idx
      cursor(oct) = cursor(oct) + 1_i32
    end do
    self%elem_order(start_idx:end_idx) = work

    child_k = 0_i32
    child_start = start_idx
    do oct = 1, 8
      if (counts(oct) <= 0_i32) cycle
      child_end = child_start + counts(oct) - 1_i32
      self%nnode = self%nnode + 1_i32
      if (self%nnode > self%max_node) error stop 'field solver tree capacity exceeded.'
      child_node = self%nnode
      child_k = child_k + 1_i32
      self%child_idx(child_k, node_idx) = child_node
      call build_node(self, mesh, child_node, child_start, child_end)
      child_start = child_end + 1_i32
    end do
    self%child_count(node_idx) = child_k

    deallocate (counts, offsets, cursor, work)
  end subroutine build_node

  !> treecode 判定に使う8分木の octant 添字を返す。
  pure integer(i32) function octant_index(x, y, z, center) result(oct)
    real(dp), intent(in) :: x, y, z
    real(dp), intent(in) :: center(3)

    oct = 1_i32
    if (x >= center(1)) oct = oct + 1_i32
    if (y >= center(2)) oct = oct + 2_i32
    if (z >= center(3)) oct = oct + 4_i32
  end function octant_index

  !> ノードを再帰走査し、葉では direct 総和、遠方は monopole 近似を適用する。
  recursive subroutine traverse_node(self, mesh, node_idx, rx, ry, rz, soft2, ex, ey, ez)
    class(field_solver_type), intent(in) :: self
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: node_idx
    real(dp), intent(in) :: rx, ry, rz, soft2
    real(dp), intent(inout) :: ex, ey, ez

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
  end subroutine traverse_node

  !> ノード半径・距離と電荷打ち消し度合いで遠方近似を判定する。
  logical function accept_node(self, node_idx, rx, ry, rz) result(accept_it)
    class(field_solver_type), intent(in) :: self
    integer(i32), intent(in) :: node_idx
    real(dp), intent(in) :: rx, ry, rz

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
  end function accept_node

  !> ノード配列を要求サイズで確保し、未使用要素をゼロ初期化する。
  subroutine ensure_tree_capacity(self, max_node_needed)
    class(field_solver_type), intent(inout) :: self
    integer(i32), intent(in) :: max_node_needed

    if (self%max_node == max_node_needed .and. allocated(self%node_start)) then
      self%node_start = 0_i32
      self%node_count = 0_i32
      self%child_count = 0_i32
      self%child_idx = 0_i32
      self%node_center = 0.0d0
      self%node_half_size = 0.0d0
      self%node_radius = 0.0d0
      self%node_q = 0.0d0
      self%node_abs_q = 0.0d0
      self%node_qx = 0.0d0
      self%node_qy = 0.0d0
      self%node_qz = 0.0d0
      self%node_charge_center = 0.0d0
      return
    end if

    call reset_tree_storage(self)
    self%max_node = max_node_needed
    allocate (self%node_start(self%max_node), self%node_count(self%max_node))
    allocate (self%child_count(self%max_node), self%child_idx(8, self%max_node))
    allocate (self%node_center(3, self%max_node), self%node_half_size(3, self%max_node))
    allocate (self%node_radius(self%max_node))
    allocate (self%node_q(self%max_node), self%node_abs_q(self%max_node))
    allocate (self%node_qx(self%max_node), self%node_qy(self%max_node), self%node_qz(self%max_node))
    allocate (self%node_charge_center(3, self%max_node))

    self%node_start = 0_i32
    self%node_count = 0_i32
    self%child_count = 0_i32
    self%child_idx = 0_i32
    self%node_center = 0.0d0
    self%node_half_size = 0.0d0
    self%node_radius = 0.0d0
    self%node_q = 0.0d0
    self%node_abs_q = 0.0d0
    self%node_qx = 0.0d0
    self%node_qy = 0.0d0
    self%node_qz = 0.0d0
    self%node_charge_center = 0.0d0
  end subroutine ensure_tree_capacity

  !> treecode 作業配列を解放する。
  subroutine reset_tree_storage(self)
    class(field_solver_type), intent(inout) :: self

    if (allocated(self%elem_order)) deallocate (self%elem_order)
    if (allocated(self%node_start)) deallocate (self%node_start)
    if (allocated(self%node_count)) deallocate (self%node_count)
    if (allocated(self%child_count)) deallocate (self%child_count)
    if (allocated(self%child_idx)) deallocate (self%child_idx)
    if (allocated(self%node_center)) deallocate (self%node_center)
    if (allocated(self%node_half_size)) deallocate (self%node_half_size)
    if (allocated(self%node_radius)) deallocate (self%node_radius)
    if (allocated(self%node_q)) deallocate (self%node_q)
    if (allocated(self%node_abs_q)) deallocate (self%node_abs_q)
    if (allocated(self%node_qx)) deallocate (self%node_qx)
    if (allocated(self%node_qy)) deallocate (self%node_qy)
    if (allocated(self%node_qz)) deallocate (self%node_qz)
    if (allocated(self%node_charge_center)) deallocate (self%node_charge_center)

    self%nnode = 0_i32
    self%max_node = 0_i32
  end subroutine reset_tree_storage

  !> ASCII英字を小文字化する。
  pure function lower_ascii(s) result(out)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: out
    integer :: i, code

    out = s
    do i = 1, len(s)
      code = iachar(out(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        out(i:i) = achar(code + 32)
      end if
    end do
  end function lower_ascii

end module bem_field_solver

!> `bem_surface_models` の浮遊導体電荷再配分を実装する submodule。
submodule(bem_surface_models) bem_surface_models_conductor
  use bem_constants, only: k_coulomb
  implicit none
contains

  module procedure init_floating_conductor_solver
  real(dp), allocatable :: saved_charge(:), baseline(:), scaled_potential(:), system(:, :)
  integer(i32) :: col, row, group_idx, nsys

  self%ncond = int(count(mesh%elem_surface_model == surface_model_conductor), kind=i32)
  self%ngroup = 0_i32
  if (self%ncond <= 0_i32) return

  allocate (self%elems(self%ncond), self%group(self%ncond), self%group_mesh_id(self%ncond))
  call collect_conductor_elements(mesh, self%elems, self%group_mesh_id, self%group, self%ngroup)

  ! 影響行列: 単位電荷を conductor 要素に 1 つずつ置き、他を 0 にした snapshot の電位から
  ! 電荷 0 の基準電位 (外場、ゲージ定数) を引く。電位は電荷に対して affine なのでこれで線形部だけが残る。
  allocate (saved_charge(mesh%nelem), baseline(self%ncond), scaled_potential(self%ncond))
  allocate (self%influence(self%ncond, self%ncond))
  saved_charge = mesh%q_elem
  mesh%q_elem = 0.0_dp
  call conductor_centroid_potentials(self, mesh, sim, snapshot, baseline)
  do col = 1, self%ncond
    mesh%q_elem = 0.0_dp
    mesh%q_elem(self%elems(col)) = 1.0_dp
    call conductor_centroid_potentials(self, mesh, sim, snapshot, scaled_potential)
    self%influence(:, col) = scaled_potential - baseline
  end do
  mesh%q_elem = saved_charge
  call snapshot%refresh(mesh)
  if (.not. all(ieee_is_finite(self%influence))) then
    error stop 'floating conductor influence matrix is not finite.'
  end if

  ! 等電位条件 (row 1..ncond) と group 総電荷保存 (row ncond+g) を並べた正方系。
  nsys = self%ncond + self%ngroup
  allocate (system(nsys, nsys))
  system = 0.0_dp
  do row = 1, self%ncond
    system(row, 1:self%ncond) = self%influence(row, :)
    system(row, self%ncond + self%group(row)) = -1.0_dp
  end do
  do group_idx = 1, self%ngroup
    row = self%ncond + group_idx
    do col = 1, self%ncond
      if (self%group(col) == group_idx) system(row, col) = 1.0_dp
    end do
  end do
  call factor_square_system(system, self%factors, self%pivot)
  end procedure init_floating_conductor_solver

  module procedure relax_floating_conductor_charges
  real(dp), allocatable :: scaled_potential(:), rhs(:), solution(:)
  integer(i32) :: row, nsys

  if (self%ncond <= 0_i32) return
  nsys = self%ncond + self%ngroup
  allocate (scaled_potential(self%ncond), rhs(nsys), solution(nsys))

  call conductor_centroid_potentials(self, mesh, sim, snapshot, scaled_potential)
  ! 固定項 = 全電荷の電位 - conductor 電荷の線形寄与。非 conductor 電荷、外場、ゲージ定数を含む。
  rhs = 0.0_dp
  rhs(1:self%ncond) = -(scaled_potential - matmul(self%influence, mesh%q_elem(self%elems)))
  do row = 1, self%ncond
    rhs(self%ncond + self%group(row)) = rhs(self%ncond + self%group(row)) + mesh%q_elem(self%elems(row))
  end do
  call solve_factored_system(self%factors, self%pivot, rhs, solution)
  if (.not. all(ieee_is_finite(solution))) then
    error stop 'floating conductor relaxation produced a non-finite charge.'
  end if
  do row = 1, self%ncond
    mesh%q_elem(self%elems(row)) = solution(row)
  end do
  end procedure relax_floating_conductor_charges

  !> 現在の `mesh%q_elem` で snapshot を refresh し、conductor 要素重心の電位を `k_coulomb` で割って返す。
  subroutine conductor_centroid_potentials(self, mesh, sim, snapshot, scaled_potential)
    class(floating_conductor_solver_type), intent(in) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(out) :: scaled_potential(:)
    real(dp) :: potential
    integer(i32) :: row

    call snapshot%refresh(mesh)
    do row = 1, self%ncond
      call snapshot%eval_local_phi(mesh, sim, mesh%centers(:, self%elems(row)), potential)
      scaled_potential(row) = potential/k_coulomb
    end do
  end subroutine conductor_centroid_potentials

  !> conductor 要素番号と、対応する conductor object group を列挙する。
  subroutine collect_conductor_elements(mesh, conductor_elems, conductor_mesh_ids, elem_group, ngroup)
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(out) :: conductor_elems(:)
    integer(i32), intent(out) :: conductor_mesh_ids(:)
    integer(i32), intent(out) :: elem_group(:)
    integer(i32), intent(out) :: ngroup

    integer(i32) :: elem_idx, out_idx, mesh_id, group_idx

    out_idx = 0_i32
    ngroup = 0_i32
    do elem_idx = 1, mesh%nelem
      if (mesh%elem_surface_model(elem_idx) /= surface_model_conductor) cycle
      mesh_id = element_mesh_id(mesh, elem_idx)
      group_idx = find_or_append_mesh_id(conductor_mesh_ids, ngroup, mesh_id)
      out_idx = out_idx + 1_i32
      conductor_elems(out_idx) = elem_idx
      elem_group(out_idx) = group_idx
    end do
  end subroutine collect_conductor_elements

  !> 既知の mesh_id なら group index を返し、未知なら末尾へ追加する。
  integer(i32) function find_or_append_mesh_id(mesh_ids, ngroup, mesh_id) result(group_idx)
    integer(i32), intent(inout) :: mesh_ids(:)
    integer(i32), intent(inout) :: ngroup
    integer(i32), intent(in) :: mesh_id
    integer(i32) :: i

    do i = 1, ngroup
      if (mesh_ids(i) == mesh_id) then
        group_idx = i
        return
      end if
    end do
    ngroup = ngroup + 1_i32
    mesh_ids(ngroup) = mesh_id
    group_idx = ngroup
  end function find_or_append_mesh_id

  !> elem_mesh_id が未割当の古い mesh でも安全に mesh_id を返す。
  integer(i32) function element_mesh_id(mesh, elem_idx) result(mesh_id)
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: elem_idx

    if (allocated(mesh%elem_mesh_id)) then
      mesh_id = mesh%elem_mesh_id(elem_idx)
    else
      mesh_id = 1_i32
    end if
  end function element_mesh_id

  !> 部分ピボット付き Gauss 消去で正方行列を LU 分解する。`factors` は L (対角 1) と U を同じ配列に持つ。
  subroutine factor_square_system(matrix, factors, pivot)
    real(dp), intent(in) :: matrix(:, :)
    real(dp), allocatable, intent(out) :: factors(:, :)
    integer(i32), allocatable, intent(out) :: pivot(:)

    real(dp), allocatable :: row_tmp(:)
    real(dp) :: pivot_abs, best_abs
    integer(i32) :: n, row, col, best

    n = int(size(matrix, 1), kind=i32)
    if (size(matrix, 2) /= n) error stop 'surface model linear system must be square.'
    allocate (factors(n, n), pivot(n), row_tmp(n))
    factors = matrix
    do col = 1, n
      best = col
      best_abs = abs(factors(col, col))
      do row = col + 1, n
        pivot_abs = abs(factors(row, col))
        if (pivot_abs > best_abs) then
          best_abs = pivot_abs
          best = row
        end if
      end do
      if (best_abs <= 1.0d-30) error stop 'surface model conductor system is singular.'
      pivot(col) = best
      if (best /= col) then
        row_tmp = factors(col, :)
        factors(col, :) = factors(best, :)
        factors(best, :) = row_tmp
      end if
      do row = col + 1, n
        factors(row, col) = factors(row, col)/factors(col, col)
        factors(row, col + 1:n) = factors(row, col + 1:n) - factors(row, col)*factors(col, col + 1:n)
      end do
    end do
  end subroutine factor_square_system

  !> LU 分解済みの系を解く。
  subroutine solve_factored_system(factors, pivot, rhs, solution)
    real(dp), intent(in) :: factors(:, :)
    integer(i32), intent(in) :: pivot(:)
    real(dp), intent(in) :: rhs(:)
    real(dp), intent(out) :: solution(:)

    real(dp) :: tmp_val
    integer(i32) :: n, row, col

    n = int(size(factors, 1), kind=i32)
    if (size(rhs) /= n .or. size(solution) /= n) error stop 'surface model linear system dimension mismatch.'
    solution = rhs
    do col = 1, n
      if (pivot(col) /= col) then
        tmp_val = solution(col)
        solution(col) = solution(pivot(col))
        solution(pivot(col)) = tmp_val
      end if
      do row = col + 1, n
        solution(row) = solution(row) - factors(row, col)*solution(col)
      end do
    end do
    do row = n, 1, -1
      tmp_val = solution(row)
      if (row < n) tmp_val = tmp_val - sum(factors(row, row + 1:n)*solution(row + 1:n))
      solution(row) = tmp_val/factors(row, row)
    end do
  end subroutine solve_factored_system

end submodule bem_surface_models_conductor

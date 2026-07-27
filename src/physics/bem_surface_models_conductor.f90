!> `bem_surface_models` の浮遊導体電荷再配分を実装する submodule。
submodule(bem_surface_models) bem_surface_models_conductor
  use bem_constants, only: k_coulomb
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_kernel, only: panel_potential_field, panel_side_principal_value
  implicit none
contains

  !> conductor 要素を列挙し、mesh_id ごとの総電荷を保って等電位化する。
  module procedure relax_floating_conductor_charges
  integer(i32), allocatable :: conductor_elems(:), conductor_mesh_ids(:), elem_group(:)
  integer(i32) :: ngroup

  allocate (conductor_elems(ncond), elem_group(ncond), conductor_mesh_ids(ncond))
  call collect_conductor_elements(mesh, conductor_elems, conductor_mesh_ids, elem_group, ngroup)
  call solve_floating_conductor_charges(mesh, external_e, conductor_elems, elem_group, ngroup)
  end procedure relax_floating_conductor_charges

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

  !> conductor 要素の電荷と object 電位を同時に解く。
  subroutine solve_floating_conductor_charges(mesh, external_e, conductor_elems, elem_group, ngroup)
    type(mesh_type), intent(inout) :: mesh
    real(dp), intent(in) :: external_e(3)
    integer(i32), intent(in) :: conductor_elems(:)
    integer(i32), intent(in) :: elem_group(:)
    integer(i32), intent(in) :: ngroup

    real(dp), allocatable :: matrix(:, :), rhs(:), solution(:), total_charge(:)
    type(panel_geometry_type), allocatable :: panel_geometry(:)
    integer(i32) :: ncond, nsys, row, col, group_idx, elem_i, elem_j, status

    ncond = int(size(conductor_elems), kind=i32)
    nsys = ncond + ngroup
    allocate (matrix(nsys, nsys), rhs(nsys), solution(nsys), total_charge(ngroup), panel_geometry(mesh%nelem))
    matrix = 0.0d0
    rhs = 0.0d0
    total_charge = 0.0d0
    do elem_j = 1_i32, mesh%nelem
      call init_panel_geometry( &
        mesh%v0(:, elem_j), mesh%v1(:, elem_j), mesh%v2(:, elem_j), panel_geometry(elem_j), status &
        )
      if (status /= panel_geometry_ok) then
        error stop 'floating conductor received invalid triangle geometry.'
      end if
    end do

    do row = 1, ncond
      elem_i = conductor_elems(row)
      group_idx = elem_group(row)
      do col = 1, ncond
        elem_j = conductor_elems(col)
        matrix(row, col) = potential_coeff(panel_geometry(elem_j), mesh%centers(:, elem_i))
      end do
      matrix(row, ncond + group_idx) = -1.0d0
      rhs(row) = -fixed_scaled_potential(mesh, panel_geometry, elem_i, external_e)
      total_charge(group_idx) = total_charge(group_idx) + mesh%q_elem(elem_i)
    end do

    do group_idx = 1, ngroup
      row = ncond + group_idx
      rhs(row) = total_charge(group_idx)
      do col = 1, ncond
        if (elem_group(col) == group_idx) matrix(row, col) = 1.0d0
      end do
    end do

    call solve_square_system(matrix, rhs, solution)
    do row = 1, ncond
      mesh%q_elem(conductor_elems(row)) = solution(row)
    end do
  end subroutine solve_floating_conductor_charges

  !> P0 triangle j の単位電荷が要素 i の重心に作る、k_coulomb で割った電位係数。
  real(dp) function potential_coeff(geometry, target) result(coeff)
    type(panel_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: target(3)
    real(dp) :: potential, field(3)

    call panel_potential_field( &
      geometry, 1.0_dp, target, panel_side_principal_value, potential, field &
      )
    coeff = potential/k_coulomb
  end function potential_coeff

  !> conductor 以外の既存電荷と一様外部電場が作る、k_coulomb で割った電位。
  real(dp) function fixed_scaled_potential(mesh, panel_geometry, elem_i, external_e) result(phi)
    type(mesh_type), intent(in) :: mesh
    type(panel_geometry_type), intent(in) :: panel_geometry(:)
    integer(i32), intent(in) :: elem_i
    real(dp), intent(in) :: external_e(3)
    integer(i32) :: elem_j

    phi = -dot_product(external_e, mesh%centers(:, elem_i))/k_coulomb
    do elem_j = 1, mesh%nelem
      if (mesh%elem_surface_model(elem_j) == surface_model_conductor) cycle
      phi = phi + mesh%q_elem(elem_j)*potential_coeff(panel_geometry(elem_j), mesh%centers(:, elem_i))
    end do
  end function fixed_scaled_potential

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

  !> 部分ピボット付き Gauss 消去で正方線形系を解く。
  subroutine solve_square_system(matrix, rhs, solution)
    real(dp), intent(in) :: matrix(:, :)
    real(dp), intent(in) :: rhs(:)
    real(dp), intent(out) :: solution(:)

    real(dp), allocatable :: work(:, :), rhs_work(:), row_tmp(:)
    real(dp) :: factor, pivot_abs, best_abs, tmp_val
    integer(i32) :: n, row, col, pivot, elim

    n = int(size(matrix, 1), kind=i32)
    if (size(matrix, 2) /= n .or. size(rhs) /= n .or. size(solution) /= n) then
      error stop 'surface model linear system dimension mismatch.'
    end if

    allocate (work(n, n), rhs_work(n), row_tmp(n))
    work = matrix
    rhs_work = rhs

    do col = 1, n
      pivot = col
      best_abs = abs(work(col, col))
      do row = col + 1, n
        pivot_abs = abs(work(row, col))
        if (pivot_abs > best_abs) then
          best_abs = pivot_abs
          pivot = row
        end if
      end do
      if (best_abs <= 1.0d-30) error stop 'surface model conductor system is singular.'
      if (pivot /= col) then
        row_tmp = work(col, :)
        work(col, :) = work(pivot, :)
        work(pivot, :) = row_tmp
        tmp_val = rhs_work(col)
        rhs_work(col) = rhs_work(pivot)
        rhs_work(pivot) = tmp_val
      end if
      do elim = col + 1, n
        factor = work(elim, col)/work(col, col)
        work(elim, col:n) = work(elim, col:n) - factor*work(col, col:n)
        rhs_work(elim) = rhs_work(elim) - factor*rhs_work(col)
      end do
    end do

    solution = 0.0d0
    do row = n, 1, -1
      tmp_val = rhs_work(row)
      if (row < n) tmp_val = tmp_val - sum(work(row, row + 1:n)*solution(row + 1:n))
      solution(row) = tmp_val/work(row, row)
    end do
  end subroutine solve_square_system

end submodule bem_surface_models_conductor

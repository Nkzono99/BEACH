!> 表面モデルごとの電荷更新後処理を扱うモジュール。
module bem_surface_models
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_types, only: mesh_type, surface_model_conductor
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  public :: apply_surface_model_charge_relaxation

contains

  !> 表面モデルに応じて、電荷堆積後の要素電荷を緩和する。
  !! conductor は mesh_id ごとの浮遊導体として、総電荷を保存しながら等電位化する。
  subroutine apply_surface_model_charge_relaxation(mesh, softening, external_e, field_bc_mode)
    type(mesh_type), intent(inout) :: mesh
    real(dp), intent(in) :: softening
    real(dp), intent(in) :: external_e(3)
    character(len=*), intent(in), optional :: field_bc_mode

    integer(i32), allocatable :: conductor_elems(:), conductor_mesh_ids(:), elem_group(:)
    integer(i32) :: ncond, ngroup

    if (.not. allocated(mesh%elem_surface_model)) return
    ncond = int(count(mesh%elem_surface_model == surface_model_conductor), kind=i32)
    if (ncond <= 0_i32) return
    call validate_conductor_field_bc(field_bc_mode)

    allocate (conductor_elems(ncond), elem_group(ncond), conductor_mesh_ids(ncond))
    call collect_conductor_elements(mesh, conductor_elems, conductor_mesh_ids, elem_group, ngroup)
    call solve_floating_conductor_charges(mesh, softening, external_e, conductor_elems, elem_group, ngroup)
  end subroutine apply_surface_model_charge_relaxation

  !> conductor 再配分が対応している場境界条件か検証する。
  subroutine validate_conductor_field_bc(field_bc_mode)
    character(len=*), intent(in), optional :: field_bc_mode
    character(len=16) :: mode

    mode = 'free'
    if (present(field_bc_mode)) mode = lower_ascii(trim(field_bc_mode))
    if (trim(mode) /= 'free') then
      error stop 'surface_model="conductor" currently requires sim.field_bc_mode="free".'
    end if
  end subroutine validate_conductor_field_bc

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
  subroutine solve_floating_conductor_charges(mesh, softening, external_e, conductor_elems, elem_group, ngroup)
    type(mesh_type), intent(inout) :: mesh
    real(dp), intent(in) :: softening
    real(dp), intent(in) :: external_e(3)
    integer(i32), intent(in) :: conductor_elems(:)
    integer(i32), intent(in) :: elem_group(:)
    integer(i32), intent(in) :: ngroup

    real(dp), allocatable :: matrix(:, :), rhs(:), solution(:), total_charge(:)
    integer(i32) :: ncond, nsys, row, col, group_idx, elem_i, elem_j

    ncond = int(size(conductor_elems), kind=i32)
    nsys = ncond + ngroup
    allocate (matrix(nsys, nsys), rhs(nsys), solution(nsys), total_charge(ngroup))
    matrix = 0.0d0
    rhs = 0.0d0
    total_charge = 0.0d0

    do row = 1, ncond
      elem_i = conductor_elems(row)
      group_idx = elem_group(row)
      do col = 1, ncond
        elem_j = conductor_elems(col)
        matrix(row, col) = potential_coeff(mesh, elem_i, elem_j, softening)
      end do
      matrix(row, ncond + group_idx) = -1.0d0
      rhs(row) = -fixed_scaled_potential(mesh, elem_i, softening, external_e)
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

  !> 要素 j の単位電荷が要素 i の重心に作る、k_coulomb で割った電位係数。
  real(dp) function potential_coeff(mesh, elem_i, elem_j, softening) result(coeff)
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: elem_i, elem_j
    real(dp), intent(in) :: softening
    real(dp), parameter :: pi_dp = acos(-1.0d0)
    real(dp) :: dx, dy, dz, r2, soft2, min_dist2

    if (elem_i == elem_j) then
      if (softening > 0.0d0) then
        coeff = 1.0d0/softening
      else
        coeff = 2.0d0*sqrt(pi_dp)/max(mesh%h_elem(elem_i), sqrt(tiny(1.0d0)))
      end if
      return
    end if

    soft2 = softening*softening
    min_dist2 = tiny(1.0d0)
    dx = mesh%center_x(elem_i) - mesh%center_x(elem_j)
    dy = mesh%center_y(elem_i) - mesh%center_y(elem_j)
    dz = mesh%center_z(elem_i) - mesh%center_z(elem_j)
    r2 = dx*dx + dy*dy + dz*dz + soft2
    coeff = 1.0d0/sqrt(max(r2, min_dist2))
  end function potential_coeff

  !> conductor 以外の既存電荷と一様外部電場が作る、k_coulomb で割った電位。
  real(dp) function fixed_scaled_potential(mesh, elem_i, softening, external_e) result(phi)
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: elem_i
    real(dp), intent(in) :: softening
    real(dp), intent(in) :: external_e(3)
    integer(i32) :: elem_j

    phi = -dot_product(external_e, mesh%centers(:, elem_i))/k_coulomb
    do elem_j = 1, mesh%nelem
      if (mesh%elem_surface_model(elem_j) == surface_model_conductor) cycle
      phi = phi + mesh%q_elem(elem_j)*potential_coeff(mesh, elem_i, elem_j, softening)
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

end module bem_surface_models

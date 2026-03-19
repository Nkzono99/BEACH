!> 三角形メッシュ幾何量(重心・法線・AABB・代表長)を前計算して保持するモジュール。
module bem_mesh
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type
  implicit none
contains

  !> 三角形頂点配列から `mesh_type` を初期化し、幾何キャッシュと要素電荷配列を準備する。
  !! @param[out] mesh 幾何キャッシュ（重心・法線・AABB）を含むメッシュ構造体。
  !! @param[in] v0 各三角形の頂点0配列 `v0(3,nelem)` [m]。
  !! @param[in] v1 各三角形の頂点1配列 `v1(3,nelem)` [m]。
  !! @param[in] v2 各三角形の頂点2配列 `v2(3,nelem)` [m]。
  !! @param[in] q0 初期要素電荷 `q_elem(nelem)` [C]（省略時は0）。
  !! @param[in] elem_mesh_id0 要素ごとのメッシュ識別ID `elem_mesh_id(nelem)`（省略時は全要素1）。
  subroutine init_mesh(mesh, v0, v1, v2, q0, elem_mesh_id0)
    type(mesh_type), intent(out) :: mesh
    real(dp), intent(in) :: v0(:, :), v1(:, :), v2(:, :)
    real(dp), intent(in), optional :: q0(:)
    integer(i32), intent(in), optional :: elem_mesh_id0(:)
    integer(i32) :: n, i
    real(dp) :: e1(3), e2(3), nvec(3), nn
    real(dp) :: cx, cy, cz

    if (size(v0, 1) /= 3 .or. size(v1, 1) /= 3 .or. size(v2, 1) /= 3) then
      error stop "mesh vertex input first dimension must be 3"
    end if

    n = size(v0, 2)
    if (size(v1, 2) /= n .or. size(v2, 2) /= n) then
      error stop "mesh vertex input size mismatch"
    end if
    mesh%nelem = n
    allocate (mesh%v0(3, n), mesh%v1(3, n), mesh%v2(3, n))
    allocate (mesh%centers(3, n), mesh%center_x(n), mesh%center_y(n), mesh%center_z(n), mesh%normals(3, n))
    allocate (mesh%bb_min(3, n), mesh%bb_max(3, n))
    allocate (mesh%h_elem(n), mesh%q_elem(n), mesh%elem_mesh_id(n))

    mesh%v0 = v0
    mesh%v1 = v1
    mesh%v2 = v2

    do i = 1, n
      cx = (v0(1, i) + v1(1, i) + v2(1, i))/3.0d0
      cy = (v0(2, i) + v1(2, i) + v2(2, i))/3.0d0
      cz = (v0(3, i) + v1(3, i) + v2(3, i))/3.0d0
      mesh%centers(1, i) = cx
      mesh%centers(2, i) = cy
      mesh%centers(3, i) = cz
      mesh%center_x(i) = cx
      mesh%center_y(i) = cy
      mesh%center_z(i) = cz
      mesh%bb_min(:, i) = min(min(v0(:, i), v1(:, i)), v2(:, i))
      mesh%bb_max(:, i) = max(max(v0(:, i), v1(:, i)), v2(:, i))

      e1 = v1(:, i) - v0(:, i)
      e2 = v2(:, i) - v0(:, i)
      nvec = cross(e1, e2)
      nn = sqrt(sum(nvec*nvec))
      if (nn > 0.0d0) then
        mesh%normals(:, i) = nvec/nn
      else
        mesh%normals(:, i) = 0.0d0
      end if
      mesh%h_elem(i) = sqrt(0.5d0*nn)
    end do

    if (present(q0)) then
      if (size(q0) /= n) error stop "q0 size mismatch"
      mesh%q_elem = q0
    else
      mesh%q_elem = 0.0d0
    end if

    if (present(elem_mesh_id0)) then
      if (size(elem_mesh_id0) /= n) error stop "elem_mesh_id0 size mismatch"
      mesh%elem_mesh_id = elem_mesh_id0
    else
      mesh%elem_mesh_id = 1_i32
    end if

    call build_collision_grid(mesh)
  end subroutine init_mesh

  !> 三角形AABBを一様グリッドへ登録し、衝突判定の候補探索を高速化する。
  !! 要素数が少ない場合は線形探索にフォールバックする。
  subroutine build_collision_grid(mesh)
    type(mesh_type), intent(inout) :: mesh

    integer(i32), parameter :: min_nelem_for_grid = 64_i32
    integer(i32), parameter :: target_elems_per_cell = 8_i32
    integer(i32), parameter :: max_cells_axis = 128_i32
    real(dp), parameter :: rel_span_eps = 1.0d-9
    real(dp), parameter :: abs_span_eps = 1.0d-12

    integer(i32) :: nx, ny, nz, ncells
    integer(i32) :: i, iaxis, ix0, ix1, iy0, iy1, iz0, iz1, ix, iy, iz, cid, total_refs, pos
    integer(i32), allocatable :: counts(:), offsets(:)
    real(dp) :: span(3), span_ref, target_cells, cell_size, eps_expand

    mesh%use_collision_grid = .false.
    mesh%grid_bb_min = minval(mesh%bb_min, dim=2)
    mesh%grid_bb_max = maxval(mesh%bb_max, dim=2)
    mesh%grid_ncell = 1_i32
    mesh%grid_inv_cell = 1.0d0

    if (allocated(mesh%grid_cell_start)) deallocate (mesh%grid_cell_start)
    if (allocated(mesh%grid_cell_elem)) deallocate (mesh%grid_cell_elem)

    if (mesh%nelem < min_nelem_for_grid) then
      allocate (mesh%grid_cell_start(2), mesh%grid_cell_elem(0))
      mesh%grid_cell_start = [1_i32, 1_i32]
      return
    end if

    span = mesh%grid_bb_max - mesh%grid_bb_min
    span_ref = max(maxval(abs(span)), 1.0d0)
    eps_expand = max(abs_span_eps, rel_span_eps*span_ref)
    do iaxis = 1, 3
      if (span(iaxis) <= abs_span_eps) then
        mesh%grid_bb_min(iaxis) = mesh%grid_bb_min(iaxis) - 0.5d0*eps_expand
        mesh%grid_bb_max(iaxis) = mesh%grid_bb_max(iaxis) + 0.5d0*eps_expand
      end if
    end do
    span = mesh%grid_bb_max - mesh%grid_bb_min

    target_cells = max(1.0d0, real(mesh%nelem, dp)/real(target_elems_per_cell, dp))
    cell_size = (span(1)*span(2)*span(3)/target_cells)**(1.0d0/3.0d0)
    if (cell_size <= 0.0d0) cell_size = span_ref

    do iaxis = 1, 3
      mesh%grid_ncell(iaxis) = int(span(iaxis)/cell_size, kind=i32)
      if (mesh%grid_ncell(iaxis) < 1_i32) mesh%grid_ncell(iaxis) = 1_i32
      if (mesh%grid_ncell(iaxis) > max_cells_axis) mesh%grid_ncell(iaxis) = max_cells_axis
      mesh%grid_inv_cell(iaxis) = real(mesh%grid_ncell(iaxis), dp)/span(iaxis)
    end do

    nx = mesh%grid_ncell(1)
    ny = mesh%grid_ncell(2)
    nz = mesh%grid_ncell(3)
    ncells = nx*ny*nz

    allocate (counts(ncells))
    counts = 0_i32
    do i = 1, mesh%nelem
      ix0 = coord_to_cell(mesh, mesh%bb_min(1, i), 1_i32)
      ix1 = coord_to_cell(mesh, mesh%bb_max(1, i), 1_i32)
      iy0 = coord_to_cell(mesh, mesh%bb_min(2, i), 2_i32)
      iy1 = coord_to_cell(mesh, mesh%bb_max(2, i), 2_i32)
      iz0 = coord_to_cell(mesh, mesh%bb_min(3, i), 3_i32)
      iz1 = coord_to_cell(mesh, mesh%bb_max(3, i), 3_i32)
      do iz = iz0, iz1
        do iy = iy0, iy1
          do ix = ix0, ix1
            cid = cell_id(ix, iy, iz, nx, ny)
            counts(cid) = counts(cid) + 1_i32
          end do
        end do
      end do
    end do

    allocate (mesh%grid_cell_start(ncells + 1))
    mesh%grid_cell_start(1) = 1_i32
    do cid = 1, ncells
      mesh%grid_cell_start(cid + 1) = mesh%grid_cell_start(cid) + counts(cid)
    end do

    total_refs = mesh%grid_cell_start(ncells + 1) - 1_i32
    if (total_refs < 0_i32) then
      error stop "collision grid CSR construction failed"
    end if
    allocate (mesh%grid_cell_elem(total_refs))

    if (total_refs > 0_i32) then
      allocate (offsets(ncells))
      offsets = mesh%grid_cell_start(1:ncells)
      do i = 1, mesh%nelem
        ix0 = coord_to_cell(mesh, mesh%bb_min(1, i), 1_i32)
        ix1 = coord_to_cell(mesh, mesh%bb_max(1, i), 1_i32)
        iy0 = coord_to_cell(mesh, mesh%bb_min(2, i), 2_i32)
        iy1 = coord_to_cell(mesh, mesh%bb_max(2, i), 2_i32)
        iz0 = coord_to_cell(mesh, mesh%bb_min(3, i), 3_i32)
        iz1 = coord_to_cell(mesh, mesh%bb_max(3, i), 3_i32)
        do iz = iz0, iz1
          do iy = iy0, iy1
            do ix = ix0, ix1
              cid = cell_id(ix, iy, iz, nx, ny)
              pos = offsets(cid)
              mesh%grid_cell_elem(pos) = i
              offsets(cid) = pos + 1_i32
            end do
          end do
        end do
      end do
      deallocate (offsets)
    end if

    deallocate (counts)
    mesh%use_collision_grid = .true.
  end subroutine build_collision_grid

  !> 座標をグリッドセル添字へ変換し、範囲外は端セルへ丸める。
  pure integer(i32) function coord_to_cell(mesh, x, axis) result(idx)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: x
    integer(i32), intent(in) :: axis
    real(dp) :: u

    if (x <= mesh%grid_bb_min(axis)) then
      idx = 1_i32
      return
    end if
    if (x >= mesh%grid_bb_max(axis)) then
      idx = mesh%grid_ncell(axis)
      return
    end if

    u = (x - mesh%grid_bb_min(axis))*mesh%grid_inv_cell(axis)
    idx = int(u, kind=i32) + 1_i32
    if (idx < 1_i32) idx = 1_i32
    if (idx > mesh%grid_ncell(axis)) idx = mesh%grid_ncell(axis)
  end function coord_to_cell

  !> 3次元セル添字 `(ix,iy,iz)` をCSR一次元インデックスへ変換する。
  pure integer(i32) function cell_id(ix, iy, iz, nx, ny) result(cid)
    integer(i32), intent(in) :: ix, iy, iz, nx, ny
    cid = (iz - 1_i32)*(nx*ny) + (iy - 1_i32)*nx + ix
  end function cell_id

  !> 3次元ベクトルの外積を返す基本演算。
  !! @param[in] a 左オペランドの3次元ベクトル。
  !! @param[in] b 右オペランドの3次元ベクトル。
  !! @return c 関数の戻り値。
  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross

end module bem_mesh

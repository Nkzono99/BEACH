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
  subroutine init_mesh(mesh, v0, v1, v2, q0)
    type(mesh_type), intent(out) :: mesh
    real(dp), intent(in) :: v0(:, :), v1(:, :), v2(:, :)
    real(dp), intent(in), optional :: q0(:)
    integer(i32) :: n, i
    real(dp) :: e1(3), e2(3), nvec(3), nn

    if (size(v0, 1) /= 3 .or. size(v1, 1) /= 3 .or. size(v2, 1) /= 3) then
      error stop "mesh vertex input first dimension must be 3"
    end if

    n = size(v0, 2)
    if (size(v1, 2) /= n .or. size(v2, 2) /= n) then
      error stop "mesh vertex input size mismatch"
    end if
    mesh%nelem = n
    allocate(mesh%v0(3, n), mesh%v1(3, n), mesh%v2(3, n))
    allocate(mesh%centers(3, n), mesh%normals(3, n))
    allocate(mesh%bb_min(3, n), mesh%bb_max(3, n))
    allocate(mesh%h_elem(n), mesh%q_elem(n))

    mesh%v0 = v0
    mesh%v1 = v1
    mesh%v2 = v2

    do i = 1, n
      mesh%centers(:, i) = (v0(:, i) + v1(:, i) + v2(:, i)) / 3.0d0
      mesh%bb_min(:, i) = min(min(v0(:, i), v1(:, i)), v2(:, i))
      mesh%bb_max(:, i) = max(max(v0(:, i), v1(:, i)), v2(:, i))

      e1 = v1(:, i) - v0(:, i)
      e2 = v2(:, i) - v0(:, i)
      nvec = cross(e1, e2)
      nn = sqrt(sum(nvec * nvec))
      if (nn > 0.0d0) then
        mesh%normals(:, i) = nvec / nn
      else
        mesh%normals(:, i) = 0.0d0
      end if
      mesh%h_elem(i) = sqrt(0.5d0 * nn)
    end do

    if (present(q0)) then
      if (size(q0) /= n) error stop "q0 size mismatch"
      mesh%q_elem = q0
    else
      mesh%q_elem = 0.0d0
    end if
  end subroutine init_mesh

  !> 3次元ベクトルの外積を返す基本演算。
  !! @param[in] a 左オペランドの3次元ベクトル。
  !! @param[in] b 右オペランドの3次元ベクトル。
  !! @return c 関数の戻り値。
  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
  end function cross

end module bem_mesh

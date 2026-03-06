!> 粒子軌道セグメントと三角形要素の交差判定を提供する衝突検出モジュール。
module bem_collision
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, hit_info
  implicit none
contains

  !> 線分 `[p0,p1]` に対して最初に衝突する三角形要素を探索し、命中情報を返す。
  !! @param[in] mesh 三角形要素とAABB情報を保持した衝突判定対象メッシュ。
  !! @param[in] p0 線分始点（粒子の移動前位置） [m]。
  !! @param[in] p1 線分終点（粒子の移動後候補位置） [m]。
  !! @param[out] hit 最初に命中した要素インデックス・命中位置・線分パラメータを格納。
  subroutine find_first_hit(mesh, p0, p1, hit)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: p0(3), p1(3)
    type(hit_info), intent(out) :: hit

    real(dp) :: d(3), seg_min(3), seg_max(3), best_t

    hit%has_hit = .false.
    hit%elem_idx = -1
    hit%t = 0.0d0
    hit%pos = 0.0d0

    d = p1 - p0
    seg_min = min(p0, p1)
    seg_max = max(p0, p1)
    best_t = huge(1.0d0)

    if (mesh%use_collision_grid) then
      call find_first_hit_grid(mesh, p0, p1, d, seg_min, seg_max, hit, best_t)
    else
      call find_first_hit_linear(mesh, p0, p1, seg_min, seg_max, hit, best_t)
    end if
  end subroutine find_first_hit

  !> 線分のAABBと要素AABBの重なりを先に判定し、詳細交差計算を枝刈りする。
  !! @param[in] p0 線分始点（粒子の移動前位置） [m]。
  !! @param[in] p1 線分終点（粒子の移動後候補位置） [m]。
  !! @param[in] bb_min 要素AABBの最小座標。
  !! @param[in] bb_max 要素AABBの最大座標。
  !! @return segment_bbox_overlap 関数の戻り値。
  pure logical function segment_bbox_overlap(p0, p1, bb_min, bb_max)
    real(dp), intent(in) :: p0(3), p1(3), bb_min(3), bb_max(3)
    real(dp) :: seg_min(3), seg_max(3)
    seg_min = min(p0, p1)
    seg_max = max(p0, p1)
    segment_bbox_overlap = segment_bbox_overlap_precomputed(seg_min, seg_max, bb_min, bb_max)
  end function segment_bbox_overlap

  !> 事前計算済みの線分AABBと要素AABBの重なりを判定する。
  pure logical function segment_bbox_overlap_precomputed(seg_min, seg_max, bb_min, bb_max)
    real(dp), intent(in) :: seg_min(3), seg_max(3), bb_min(3), bb_max(3)
    segment_bbox_overlap_precomputed = all(bb_max >= seg_min) .and. all(bb_min <= seg_max)
  end function segment_bbox_overlap_precomputed

  !> 旧実装と同じ線形探索で最初の命中要素を探索する。
  subroutine find_first_hit_linear(mesh, p0, p1, seg_min, seg_max, hit, best_t)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: p0(3), p1(3), seg_min(3), seg_max(3)
    type(hit_info), intent(inout) :: hit
    real(dp), intent(inout) :: best_t

    integer(i32) :: i
    logical :: ok
    real(dp) :: t, h(3)

    do i = 1, mesh%nelem
      if (.not. segment_bbox_overlap_precomputed(seg_min, seg_max, mesh%bb_min(:, i), mesh%bb_max(:, i))) cycle
      call segment_triangle_intersect(p0, p1, mesh%v0(:, i), mesh%v1(:, i), mesh%v2(:, i), ok, t, h)
      if (.not. ok) cycle
      if (t < best_t) then
        best_t = t
        hit%has_hit = .true.
        hit%elem_idx = i
        hit%t = t
        hit%pos = h
      end if
    end do
  end subroutine find_first_hit_linear

  !> 一様グリッド + 3D-DDA で候補セルのみ探索し、最初の命中要素を返す。
  subroutine find_first_hit_grid(mesh, p0, p1, d, seg_min, seg_max, hit, best_t)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: p0(3), p1(3), d(3), seg_min(3), seg_max(3)
    type(hit_info), intent(inout) :: hit
    real(dp), intent(inout) :: best_t

    real(dp), parameter :: eps = 1.0d-12
    real(dp) :: t_entry, t_exit, t_cur, t_next
    real(dp) :: t_max(3), t_delta(3), cell_size
    real(dp) :: t, h(3), p_entry(3)
    integer(i32) :: axis, nx, ny, cid
    integer(i32) :: cell(3), step(3)
    integer(i32) :: k, start_idx, end_idx, elem_idx
    logical :: ok, hit_grid

    call segment_aabb_intersection_t(p0, d, mesh%grid_bb_min, mesh%grid_bb_max, hit_grid, t_entry, t_exit)
    if (.not. hit_grid) return
    if (t_exit < 0.0d0 .or. t_entry > 1.0d0) return

    t_cur = max(0.0d0, t_entry)
    if (t_cur > t_exit) return

    p_entry = p0 + t_cur * d
    do axis = 1, 3
      cell(axis) = coord_to_cell(mesh, p_entry(axis), int(axis, kind=i32))
      if (abs(d(axis)) <= eps) then
        step(axis) = 0_i32
        t_max(axis) = huge(1.0d0)
        t_delta(axis) = huge(1.0d0)
      else
        cell_size = 1.0d0 / mesh%grid_inv_cell(axis)
        if (d(axis) > 0.0d0) then
          step(axis) = 1_i32
          t_max(axis) = (mesh%grid_bb_min(axis) + real(cell(axis), dp) * cell_size - p0(axis)) / d(axis)
          t_delta(axis) = cell_size / d(axis)
        else
          step(axis) = -1_i32
          t_max(axis) = (mesh%grid_bb_min(axis) + real(cell(axis) - 1_i32, dp) * cell_size - p0(axis)) / d(axis)
          t_delta(axis) = -cell_size / d(axis)
        end if
        do while (t_max(axis) < t_cur - eps)
          t_max(axis) = t_max(axis) + t_delta(axis)
        end do
      end if
    end do

    nx = mesh%grid_ncell(1)
    ny = mesh%grid_ncell(2)

    do
      if (t_cur > t_exit + eps) exit
      if (t_cur > best_t + eps) exit

      cid = cell_id(cell(1), cell(2), cell(3), nx, ny)
      start_idx = mesh%grid_cell_start(cid)
      end_idx = mesh%grid_cell_start(cid + 1_i32) - 1_i32
      do k = start_idx, end_idx
        elem_idx = mesh%grid_cell_elem(k)
        if (.not. segment_bbox_overlap_precomputed(seg_min, seg_max, mesh%bb_min(:, elem_idx), mesh%bb_max(:, elem_idx))) cycle
        call segment_triangle_intersect( &
          p0, p1, mesh%v0(:, elem_idx), mesh%v1(:, elem_idx), mesh%v2(:, elem_idx), ok, t, h &
        )
        if (.not. ok) cycle
        if (t < best_t) then
          best_t = t
          hit%has_hit = .true.
          hit%elem_idx = elem_idx
          hit%t = t
          hit%pos = h
        end if
      end do

      t_next = min(t_max(1), min(t_max(2), t_max(3)))
      if (t_next > t_exit + eps) exit
      if (t_next > best_t + eps) exit

      do axis = 1, 3
        if (t_max(axis) <= t_next + eps) then
          if (step(axis) /= 0_i32) then
            cell(axis) = cell(axis) + step(axis)
            if (cell(axis) < 1_i32 .or. cell(axis) > mesh%grid_ncell(axis)) return
            t_max(axis) = t_max(axis) + t_delta(axis)
          end if
        end if
      end do
      t_cur = t_next
    end do
  end subroutine find_first_hit_grid

  !> 線分 `p(t)=p0+t*d` (`0<=t<=1`) とAABBの交差区間 `[t_entry,t_exit]` を返す。
  pure subroutine segment_aabb_intersection_t(p0, d, bb_min, bb_max, ok, t_entry, t_exit)
    real(dp), intent(in) :: p0(3), d(3), bb_min(3), bb_max(3)
    logical, intent(out) :: ok
    real(dp), intent(out) :: t_entry, t_exit

    real(dp), parameter :: eps = 1.0d-14
    real(dp) :: t0, t1, t_near, t_far, inv_d, tmp
    integer(i32) :: axis

    t0 = 0.0d0
    t1 = 1.0d0
    do axis = 1, 3
      if (abs(d(axis)) <= eps) then
        if (p0(axis) < bb_min(axis) .or. p0(axis) > bb_max(axis)) then
          ok = .false.
          t_entry = 0.0d0
          t_exit = -1.0d0
          return
        end if
      else
        inv_d = 1.0d0 / d(axis)
        t_near = (bb_min(axis) - p0(axis)) * inv_d
        t_far = (bb_max(axis) - p0(axis)) * inv_d
        if (t_near > t_far) then
          tmp = t_near
          t_near = t_far
          t_far = tmp
        end if
        if (t_near > t0) t0 = t_near
        if (t_far < t1) t1 = t_far
        if (t0 > t1) then
          ok = .false.
          t_entry = 0.0d0
          t_exit = -1.0d0
          return
        end if
      end if
    end do

    ok = .true.
    t_entry = t0
    t_exit = t1
  end subroutine segment_aabb_intersection_t

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

    u = (x - mesh%grid_bb_min(axis)) * mesh%grid_inv_cell(axis)
    idx = int(u, kind=i32) + 1_i32
    if (idx < 1_i32) idx = 1_i32
    if (idx > mesh%grid_ncell(axis)) idx = mesh%grid_ncell(axis)
  end function coord_to_cell

  !> 3次元セル添字 `(ix,iy,iz)` をCSR一次元インデックスへ変換する。
  pure integer(i32) function cell_id(ix, iy, iz, nx, ny) result(cid)
    integer(i32), intent(in) :: ix, iy, iz, nx, ny
    cid = (iz - 1_i32) * (nx * ny) + (iy - 1_i32) * nx + ix
  end function cell_id

  !> Möller–Trumbore法で線分と三角形の交差有無・線分パラメータ `t`・交点座標を計算する。
  !! @param[in] p0 線分始点（粒子の移動前位置） [m]。
  !! @param[in] p1 線分終点（粒子の移動後候補位置） [m]。
  !! @param[in] v0 三角形頂点0の座標。
  !! @param[in] v1 三角形頂点1の座標。
  !! @param[in] v2 三角形頂点2の座標。
  !! @param[out] ok 線分と三角形が交差した場合に `.true.`。
  !! @param[out] t 交点の線分内パラメータ（`p0 + t*(p1-p0)`）。
  !! @param[out] h 交点座標。
  subroutine segment_triangle_intersect(p0, p1, v0, v1, v2, ok, t, h)
    real(dp), intent(in) :: p0(3), p1(3), v0(3), v1(3), v2(3)
    logical, intent(out) :: ok
    real(dp), intent(out) :: t, h(3)

    real(dp), parameter :: eps = 1.0d-12
    real(dp) :: d(3), e1(3), e2(3), q(3), s(3), hh(3), a, f, u, v

    d = p1 - p0
    e1 = v1 - v0
    e2 = v2 - v0
    hh = cross(d, e2)
    a = dot_product(e1, hh)
    if (abs(a) < eps) then
      ok = .false.; return
    end if

    f = 1.0d0 / a
    s = p0 - v0
    u = f * dot_product(s, hh)
    if (u < 0.0d0 .or. u > 1.0d0) then
      ok = .false.; return
    end if

    q = cross(s, e1)
    v = f * dot_product(d, q)
    if (v < 0.0d0 .or. (u + v) > 1.0d0) then
      ok = .false.; return
    end if

    t = f * dot_product(e2, q)
    if (t < 0.0d0 .or. t > 1.0d0) then
      ok = .false.; return
    end if

    h = p0 + t * d
    ok = .true.
  end subroutine segment_triangle_intersect

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

end module bem_collision

!> 粒子軌道セグメントと三角形要素の交差判定を提供する衝突検出モジュール。
module bem_collision
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, hit_info, sim_config, bc_periodic
  use bem_string_utils, only: lower_ascii
  implicit none
contains

  !> 線分 `[p0,p1]` に対して最初に衝突する三角形要素を探索し、命中情報を返す。
  !! @param[in] mesh 三角形要素とAABB情報を保持した衝突判定対象メッシュ。
  !! @param[in] p0 線分始点（粒子の移動前位置） [m]。
  !! @param[in] p1 線分終点（粒子の移動後候補位置） [m]。
  !! @param[out] hit 最初に命中した要素インデックス・命中位置・線分パラメータを格納。
  subroutine find_first_hit(mesh, p0, p1, hit, sim, box_min, box_max, require_elem_inside)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: p0(3), p1(3)
    type(hit_info), intent(out) :: hit
    type(sim_config), intent(in), optional :: sim
    real(dp), intent(in), optional :: box_min(3), box_max(3)
    logical, intent(in), optional :: require_elem_inside

    logical :: use_periodic2
    integer(i32) :: periodic_axes(2)
    real(dp) :: periodic_len(2)

    use_periodic2 = .false.
    periodic_axes = 0_i32
    periodic_len = 0.0d0
    if (present(sim)) then
      call resolve_periodic2_collision_config(sim, use_periodic2, periodic_axes, periodic_len)
    end if

    if (use_periodic2) then
      call find_first_hit_periodic2(mesh, p0, p1, hit, sim, box_min, box_max, require_elem_inside)
    else
      call find_first_hit_base(mesh, p0, p1, hit, box_min, box_max, require_elem_inside)
    end if
  end subroutine find_first_hit

  !> 通常メッシュに対する最初の命中要素探索を行う。
  subroutine find_first_hit_base(mesh, p0, p1, hit, box_min, box_max, require_elem_inside)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: p0(3), p1(3)
    type(hit_info), intent(out) :: hit
    real(dp), intent(in), optional :: box_min(3), box_max(3)
    logical, intent(in), optional :: require_elem_inside

    real(dp) :: d(3), seg_min(3), seg_max(3), best_t
    real(dp) :: box_min_local(3), box_max_local(3), box_tol
    logical :: use_box_filter, require_inside_elem

    call initialize_hit(hit)

    d = p1 - p0
    seg_min = min(p0, p1)
    seg_max = max(p0, p1)
    best_t = huge(1.0d0)

    call resolve_box_filter_args( &
      box_min, box_max, require_elem_inside, use_box_filter, box_min_local, box_max_local, box_tol, require_inside_elem &
      )

    if (mesh%use_collision_grid) then
      call find_first_hit_base_grid( &
        mesh, p0, p1, d, seg_min, seg_max, hit, best_t, &
        use_box_filter, box_min_local, box_max_local, box_tol, require_inside_elem &
        )
    else
      call find_first_hit_base_linear( &
        mesh, p0, p1, seg_min, seg_max, hit, best_t, &
        use_box_filter, box_min_local, box_max_local, box_tol, require_inside_elem &
        )
    end if
  end subroutine find_first_hit_base

  !> periodic2 用に image shift を列挙し、base collision の結果を統合する。
  subroutine find_first_hit_periodic2(mesh, p0, p1, hit, sim, box_min, box_max, require_elem_inside)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: p0(3), p1(3)
    type(hit_info), intent(out) :: hit
    type(sim_config), intent(in) :: sim
    real(dp), intent(in), optional :: box_min(3), box_max(3)
    logical, intent(in), optional :: require_elem_inside

    integer(i32) :: periodic_axes(2), image_shift(2), nmin(2), nmax(2), iaxis, n1, n2
    real(dp) :: periodic_len(2), shift_vec(3), candidate_pos(3), candidate_wrapped(3)
    real(dp) :: box_min_local(3), box_max_local(3), box_tol
    real(dp) :: shifted_p0(3), shifted_p1(3)
    type(hit_info) :: candidate
    logical :: use_periodic2, use_box_filter, require_inside_elem

    call initialize_hit(hit)
    if (.not. mesh%periodic2_collision_ready) then
      error stop 'periodic2 collision requires prepare_periodic2_collision_mesh before ray queries.'
    end if

    call resolve_periodic2_collision_config(sim, use_periodic2, periodic_axes, periodic_len)
    if (.not. use_periodic2) then
      call find_first_hit_base(mesh, p0, p1, hit, box_min, box_max, require_elem_inside)
      return
    end if

    call resolve_box_filter_args( &
      box_min, box_max, require_elem_inside, use_box_filter, box_min_local, box_max_local, box_tol, require_inside_elem &
      )

    do iaxis = 1, 2
      call compute_periodic_shift_bounds(mesh, p0, p1, periodic_axes(iaxis), periodic_len(iaxis), nmin(iaxis), nmax(iaxis))
      if (nmin(iaxis) > nmax(iaxis)) return
    end do

    do n1 = nmin(1), nmax(1)
      do n2 = nmin(2), nmax(2)
        image_shift = [n1, n2]
        shift_vec = 0.0d0
        shift_vec(periodic_axes(1)) = real(image_shift(1), dp)*periodic_len(1)
        shift_vec(periodic_axes(2)) = real(image_shift(2), dp)*periodic_len(2)
        shifted_p0 = p0 - shift_vec
        shifted_p1 = p1 - shift_vec

        call find_first_hit_base(mesh, shifted_p0, shifted_p1, candidate)
        if (.not. candidate%has_hit) cycle

        candidate_pos = candidate%pos + shift_vec
        candidate_wrapped = candidate_pos
        call wrap_periodic2_point(candidate_wrapped, sim%box_min, periodic_axes, periodic_len)
        if (use_box_filter) then
          if (.not. point_inside_box_periodic2( &
              candidate_wrapped, box_min_local, box_max_local, box_tol, periodic_axes, require_inside_elem)) cycle
        end if

        if (prefer_periodic_candidate(candidate%t, candidate%elem_idx, image_shift, hit)) then
          hit%has_hit = .true.
          hit%elem_idx = candidate%elem_idx
          hit%t = candidate%t
          hit%pos = candidate_pos
          hit%pos_wrapped = candidate_wrapped
          hit%image_shift = image_shift
        end if
      end do
    end do
  end subroutine find_first_hit_periodic2

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
  subroutine find_first_hit_base_linear( &
    mesh, p0, p1, seg_min, seg_max, hit, best_t, use_box_filter, box_min, box_max, box_tol, require_elem_inside &
    )
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: p0(3), p1(3), seg_min(3), seg_max(3)
    type(hit_info), intent(inout) :: hit
    real(dp), intent(inout) :: best_t
    logical, intent(in) :: use_box_filter
    real(dp), intent(in) :: box_min(3), box_max(3), box_tol
    logical, intent(in) :: require_elem_inside

    integer(i32) :: i
    logical :: ok
    real(dp) :: t, h(3)

    do i = 1, mesh%nelem
      if (use_box_filter) then
        if (.not. segment_bbox_overlap_precomputed(mesh%bb_min(:, i), mesh%bb_max(:, i), box_min, box_max)) cycle
        if (require_elem_inside) then
          if (.not. bbox_inside_box(mesh%bb_min(:, i), mesh%bb_max(:, i), box_min, box_max, box_tol)) cycle
        end if
      end if
      if (.not. segment_bbox_overlap_precomputed(seg_min, seg_max, mesh%bb_min(:, i), mesh%bb_max(:, i))) cycle
      call segment_triangle_intersect(p0, p1, mesh%v0(:, i), mesh%v1(:, i), mesh%v2(:, i), ok, t, h)
      if (.not. ok) cycle
      if (use_box_filter) then
        if (.not. point_inside_box(h, box_min, box_max, box_tol)) cycle
      end if
      if (t < best_t) then
        best_t = t
        hit%has_hit = .true.
        hit%elem_idx = i
        hit%t = t
        hit%pos = h
        hit%pos_wrapped = h
        hit%image_shift = 0_i32
      end if
    end do
  end subroutine find_first_hit_base_linear

  !> 一様グリッド + 3D-DDA で候補セルのみ探索し、最初の命中要素を返す。
  subroutine find_first_hit_base_grid( &
    mesh, p0, p1, d, seg_min, seg_max, hit, best_t, use_box_filter, box_min, box_max, box_tol, require_elem_inside &
    )
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: p0(3), p1(3), d(3), seg_min(3), seg_max(3)
    type(hit_info), intent(inout) :: hit
    real(dp), intent(inout) :: best_t
    logical, intent(in) :: use_box_filter
    real(dp), intent(in) :: box_min(3), box_max(3), box_tol
    logical, intent(in) :: require_elem_inside

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

    p_entry = p0 + t_cur*d
    do axis = 1, 3
      cell(axis) = coord_to_cell(mesh, p_entry(axis), int(axis, kind=i32))
      if (abs(d(axis)) <= eps) then
        step(axis) = 0_i32
        t_max(axis) = huge(1.0d0)
        t_delta(axis) = huge(1.0d0)
      else
        cell_size = 1.0d0/mesh%grid_inv_cell(axis)
        if (d(axis) > 0.0d0) then
          step(axis) = 1_i32
          t_max(axis) = (mesh%grid_bb_min(axis) + real(cell(axis), dp)*cell_size - p0(axis))/d(axis)
          t_delta(axis) = cell_size/d(axis)
        else
          step(axis) = -1_i32
          t_max(axis) = (mesh%grid_bb_min(axis) + real(cell(axis) - 1_i32, dp)*cell_size - p0(axis))/d(axis)
          t_delta(axis) = -cell_size/d(axis)
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
        if (use_box_filter) then
          if (.not. segment_bbox_overlap_precomputed( &
              mesh%bb_min(:, elem_idx), mesh%bb_max(:, elem_idx), box_min, box_max)) cycle
          if (require_elem_inside) then
            if (.not. bbox_inside_box( &
                mesh%bb_min(:, elem_idx), mesh%bb_max(:, elem_idx), box_min, box_max, box_tol)) cycle
          end if
        end if
        if (.not. segment_bbox_overlap_precomputed( &
            seg_min, seg_max, mesh%bb_min(:, elem_idx), mesh%bb_max(:, elem_idx))) cycle
        call segment_triangle_intersect( &
          p0, p1, mesh%v0(:, elem_idx), mesh%v1(:, elem_idx), mesh%v2(:, elem_idx), ok, t, h &
          )
        if (.not. ok) cycle
        if (use_box_filter) then
          if (.not. point_inside_box(h, box_min, box_max, box_tol)) cycle
        end if
        if (t < best_t) then
          best_t = t
          hit%has_hit = .true.
          hit%elem_idx = elem_idx
          hit%t = t
          hit%pos = h
          hit%pos_wrapped = h
          hit%image_shift = 0_i32
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
  end subroutine find_first_hit_base_grid

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
        inv_d = 1.0d0/d(axis)
        t_near = (bb_min(axis) - p0(axis))*inv_d
        t_far = (bb_max(axis) - p0(axis))*inv_d
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

  pure logical function point_inside_box(p, box_min, box_max, tol)
    real(dp), intent(in) :: p(3), box_min(3), box_max(3), tol
    point_inside_box = all(p >= (box_min - tol)) .and. all(p <= (box_max + tol))
  end function point_inside_box

  pure logical function point_inside_box_periodic2(p, box_min, box_max, tol, periodic_axes, require_half_open)
    real(dp), intent(in) :: p(3), box_min(3), box_max(3), tol
    integer(i32), intent(in) :: periodic_axes(2)
    logical, intent(in) :: require_half_open
    integer(i32) :: axis
    logical :: is_periodic

    point_inside_box_periodic2 = .true.
    do axis = 1_i32, 3_i32
      is_periodic = any(periodic_axes == axis)
      if (require_half_open .and. is_periodic) then
        if (p(axis) < box_min(axis) - tol .or. p(axis) >= box_max(axis) + tol) then
          point_inside_box_periodic2 = .false.
          return
        end if
      else
        if (p(axis) < box_min(axis) - tol .or. p(axis) > box_max(axis) + tol) then
          point_inside_box_periodic2 = .false.
          return
        end if
      end if
    end do
  end function point_inside_box_periodic2

  pure logical function bbox_inside_box(bb_min, bb_max, box_min, box_max, tol)
    real(dp), intent(in) :: bb_min(3), bb_max(3), box_min(3), box_max(3), tol
    bbox_inside_box = all(bb_min >= (box_min - tol)) .and. all(bb_max <= (box_max + tol))
  end function bbox_inside_box

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

    f = 1.0d0/a
    s = p0 - v0
    u = f*dot_product(s, hh)
    if (u < 0.0d0 .or. u > 1.0d0) then
      ok = .false.; return
    end if

    q = cross(s, e1)
    v = f*dot_product(d, q)
    if (v < 0.0d0 .or. (u + v) > 1.0d0) then
      ok = .false.; return
    end if

    t = f*dot_product(e2, q)
    if (t < 0.0d0 .or. t > 1.0d0) then
      ok = .false.; return
    end if

    h = p0 + t*d
    ok = .true.
  end subroutine segment_triangle_intersect

  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross

  !> hit 構造体を未命中状態へ初期化する。
  subroutine initialize_hit(hit)
    type(hit_info), intent(out) :: hit

    hit%has_hit = .false.
    hit%elem_idx = -1_i32
    hit%t = 0.0d0
    hit%pos = 0.0d0
    hit%image_shift = 0_i32
    hit%pos_wrapped = 0.0d0
  end subroutine initialize_hit

  !> box filter 関連の optional 引数を検証付きで展開する。
  subroutine resolve_box_filter_args( &
    box_min, box_max, require_elem_inside, use_box_filter, box_min_local, box_max_local, box_tol, require_inside_elem &
    )
    real(dp), intent(in), optional :: box_min(3), box_max(3)
    logical, intent(in), optional :: require_elem_inside
    logical, intent(out) :: use_box_filter, require_inside_elem
    real(dp), intent(out) :: box_min_local(3), box_max_local(3), box_tol

    use_box_filter = present(box_min) .or. present(box_max)
    if (use_box_filter .and. .not. (present(box_min) .and. present(box_max))) then
      error stop 'find_first_hit requires both box_min and box_max when using box filter.'
    end if
    if (use_box_filter) then
      box_min_local = box_min
      box_max_local = box_max
      box_tol = 1.0d-12*max(1.0d0, maxval(abs(box_max_local - box_min_local)))
    else
      box_min_local = 0.0d0
      box_max_local = 0.0d0
      box_tol = 0.0d0
    end if

    require_inside_elem = .false.
    if (present(require_elem_inside)) require_inside_elem = require_elem_inside
    if (require_inside_elem .and. .not. use_box_filter) then
      error stop 'find_first_hit require_elem_inside=true needs box_min/box_max.'
    end if
  end subroutine resolve_box_filter_args

  !> periodic2 collision で必要な 2 軸周期設定を解決する。
  subroutine resolve_periodic2_collision_config(sim, use_periodic2, periodic_axes, periodic_len)
    type(sim_config), intent(in) :: sim
    logical, intent(out) :: use_periodic2
    integer(i32), intent(out) :: periodic_axes(2)
    real(dp), intent(out) :: periodic_len(2)

    character(len=16) :: field_bc_mode
    integer(i32) :: axis, n_periodic
    real(dp) :: span

    use_periodic2 = .false.
    periodic_axes = 0_i32
    periodic_len = 0.0d0
    field_bc_mode = lower_ascii(trim(sim%field_bc_mode))
    if (trim(field_bc_mode) /= 'periodic2') return

    if (.not. sim%use_box) then
      error stop 'sim.field_bc_mode="periodic2" requires sim.use_box=true.'
    end if

    n_periodic = 0_i32
    do axis = 1_i32, 3_i32
      if ((sim%bc_low(axis) == bc_periodic) .neqv. (sim%bc_high(axis) == bc_periodic)) then
        error stop 'periodic2 requires bc_low(axis)=bc_high(axis)=periodic for periodic axes.'
      end if
      if (sim%bc_low(axis) == bc_periodic) then
        n_periodic = n_periodic + 1_i32
        if (n_periodic <= 2_i32) periodic_axes(n_periodic) = axis
      end if
    end do
    if (n_periodic /= 2_i32) then
      error stop 'sim.field_bc_mode="periodic2" requires exactly two periodic axes.'
    end if

    do axis = 1_i32, 2_i32
      span = sim%box_max(periodic_axes(axis)) - sim%box_min(periodic_axes(axis))
      if (span <= 0.0d0) error stop 'periodic2 requires positive box length on periodic axes.'
      periodic_len(axis) = span
    end do

    use_periodic2 = .true.
  end subroutine resolve_periodic2_collision_config

  !> 線分 AABB と canonical mesh AABB の重なりから必要な image shift 範囲を決める。
  subroutine compute_periodic_shift_bounds(mesh, p0, p1, axis, period_len, nmin, nmax)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: p0(3), p1(3)
    integer(i32), intent(in) :: axis
    real(dp), intent(in) :: period_len
    integer(i32), intent(out) :: nmin, nmax

    real(dp) :: seg_min, seg_max, mesh_min, mesh_max, tol

    seg_min = min(p0(axis), p1(axis))
    seg_max = max(p0(axis), p1(axis))
    mesh_min = mesh%grid_bb_min(axis)
    mesh_max = mesh%grid_bb_max(axis)
    tol = 1.0d-12*max(1.0d0, abs(seg_min), abs(seg_max), abs(mesh_min), abs(mesh_max), period_len)

    nmin = int(ceiling((seg_min - mesh_max - tol)/period_len), kind=i32)
    nmax = int(floor((seg_max - mesh_min + tol)/period_len), kind=i32)
  end subroutine compute_periodic_shift_bounds

  !> point を primary cell へ折り返す。
  subroutine wrap_periodic2_point(point, box_min, periodic_axes, periodic_len)
    real(dp), intent(inout) :: point(3)
    real(dp), intent(in) :: box_min(3), periodic_len(2)
    integer(i32), intent(in) :: periodic_axes(2)
    integer(i32) :: iaxis

    do iaxis = 1, 2
      point(periodic_axes(iaxis)) = box_min(periodic_axes(iaxis)) + &
                                    modulo(point(periodic_axes(iaxis)) - box_min(periodic_axes(iaxis)), periodic_len(iaxis))
    end do
  end subroutine wrap_periodic2_point

  !> 候補 hit が現在の best より優先されるかを deterministic に判定する。
  pure logical function prefer_periodic_candidate(t, elem_idx, image_shift, best_hit)
    real(dp), intent(in) :: t
    integer(i32), intent(in) :: elem_idx, image_shift(2)
    type(hit_info), intent(in) :: best_hit
    real(dp) :: tol

    if (.not. best_hit%has_hit) then
      prefer_periodic_candidate = .true.
      return
    end if

    tol = 1.0d-12*max(1.0d0, abs(t), abs(best_hit%t))
    if (t < best_hit%t - tol) then
      prefer_periodic_candidate = .true.
      return
    end if
    if (t > best_hit%t + tol) then
      prefer_periodic_candidate = .false.
      return
    end if

    if (elem_idx < best_hit%elem_idx) then
      prefer_periodic_candidate = .true.
      return
    end if
    if (elem_idx > best_hit%elem_idx) then
      prefer_periodic_candidate = .false.
      return
    end if

    if (image_shift(1) < best_hit%image_shift(1)) then
      prefer_periodic_candidate = .true.
    else if (image_shift(1) > best_hit%image_shift(1)) then
      prefer_periodic_candidate = .false.
    else
      prefer_periodic_candidate = image_shift(2) < best_hit%image_shift(2)
    end if
  end function prefer_periodic_candidate

end module bem_collision

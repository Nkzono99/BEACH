module bem_collision
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, hit_info
  implicit none
contains

  subroutine find_first_hit(mesh, p0, p1, hit)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: p0(3), p1(3)
    type(hit_info), intent(out) :: hit

    integer(i32) :: i
    logical :: ok
    real(dp) :: t, best_t, h(3)

    hit%has_hit = .false.
    hit%elem_idx = -1
    hit%t = 0.0d0
    hit%pos = 0.0d0

    best_t = huge(1.0d0)
    do i = 1, mesh%nelem
      if (.not. segment_bbox_overlap(p0, p1, mesh%bb_min(:, i), mesh%bb_max(:, i))) cycle
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
  end subroutine find_first_hit

  pure logical function segment_bbox_overlap(p0, p1, bb_min, bb_max)
    real(dp), intent(in) :: p0(3), p1(3), bb_min(3), bb_max(3)
    real(dp) :: seg_min(3), seg_max(3)
    seg_min = min(p0, p1)
    seg_max = max(p0, p1)
    segment_bbox_overlap = all(bb_max >= seg_min) .and. all(bb_min <= seg_max)
  end function segment_bbox_overlap

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

  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
  end function cross

end module bem_collision

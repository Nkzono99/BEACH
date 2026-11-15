!> Analytic free-space P0 triangle potential, field, principal value, and jump.
module bem_panel_kernel
  use bem_kinds, only: dp, sp, i32
  use bem_constants, only: eps0, k_coulomb
  use bem_panel_geometry, only: panel_geometry_type
  use bem_panel_self_terms, only: panel_on_surface_integrals, panel_edge_line_integral, panel_edge_ratio_limit
  implicit none
  private

  integer(i32), parameter, public :: panel_side_principal_value = 0_i32
  integer(i32), parameter, public :: panel_side_normal_plus = 1_i32
  integer(i32), parameter, public :: panel_side_normal_minus = -1_i32

  public :: panel_potential_field
  public :: panel_potential_field_translated
  public :: panel_field_sum_translated
  public :: panel_potential_sum_translated
  public :: panel_field_sum_translated_sp
  public :: panel_potential_sum_translated_sp
  public :: vector_log_sp, vector_atan2_sp

  !> Panels gathered per SIMD pass of the near-field sums.
  integer(i32), parameter :: panel_chunk = 64_i32
  !> Single-precision sums hand lanes whose edge ratio e/(ra+rb) reaches this value to the double kernel.
  real(sp), parameter :: near_edge_single_limit = 0.99_sp

contains

  subroutine panel_potential_field(geometry, charge, target, side, potential, field)
    type(panel_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: charge, target(3)
    integer(i32), intent(in) :: side
    real(dp), intent(out) :: potential, field(3)
    real(dp) :: projection(3), height, potential_integral, field_integral(3), omega
    real(dp) :: endpoint_a(3), radius_a, radius_b, line_integral, distance, ratio
    real(dp) :: relative_vertex(3, 3), radius(3), numerator, denominator, scale, plane_tolerance
    integer :: edge, next_edge
    logical :: on_surface, near_edge

    if (side < panel_side_normal_minus .or. side > panel_side_normal_plus) then
      error stop 'invalid panel evaluation side.'
    end if
    potential = 0.0_dp
    field = 0.0_dp
    if (charge == 0.0_dp) return
    height = dot_product(target - geometry%vertex(:, 1), geometry%normal)
    projection = target - height*geometry%normal
    scale = maxval(geometry%edge_length)
    plane_tolerance = 64.0_dp*epsilon(1.0_dp)*scale
    on_surface = abs(height) <= plane_tolerance

    ! Do not project a finite-height near-edge target onto a singular edge. The plane
    ! tolerance still resolves the principal value and side jump away from edges.
    if (on_surface .and. height /= 0.0_dp) then
      do edge = 1, 3
        next_edge = merge(edge + 1, 1, edge < 3)
        radius_a = norm2(target - geometry%vertex(:, edge))
        radius_b = norm2(target - geometry%vertex(:, next_edge))
        if (geometry%edge_length(edge)/(radius_a + radius_b) >= panel_edge_ratio_limit) then
          on_surface = .false.
          exit
        end if
      end do
    end if

    if (on_surface) then
      call panel_on_surface_integrals(geometry, projection, potential_integral, field_integral)
      potential = k_coulomb*charge/geometry%area*potential_integral
      field = k_coulomb*charge/geometry%area*field_integral
      if (point_in_triangle(geometry, projection)) then
        field = field + real(side, dp)*charge/(2.0_dp*geometry%area*eps0)*geometry%normal
      end if
      return
    end if

    omega = 0.0_dp
    do edge = 1, 3
      relative_vertex(:, edge) = geometry%vertex(:, edge) - target
      radius(edge) = norm2(relative_vertex(:, edge))
    end do
    near_edge = .false.
    do edge = 1, 3
      next_edge = merge(edge + 1, 1, edge < 3)
      near_edge = near_edge .or. geometry%edge_length(edge)/(radius(edge) + radius(next_edge)) >= panel_edge_ratio_limit
    end do
    numerator = dot_product(relative_vertex(:, 1), &
                            cross_product(relative_vertex(:, 2), relative_vertex(:, 3)))
    denominator = radius(1)*radius(2)*radius(3) + &
                  dot_product(relative_vertex(:, 1), relative_vertex(:, 2))*radius(3) + &
                  dot_product(relative_vertex(:, 2), relative_vertex(:, 3))*radius(1) + &
                  dot_product(relative_vertex(:, 3), relative_vertex(:, 1))*radius(2)
    if (near_edge) then
      omega = near_edge_solid_angle(geometry, target, projection, height)
    else
      omega = -2.0_dp*atan2(numerator, denominator)
    end if

    potential_integral = -height*omega
    field_integral = geometry%normal*omega
    do edge = 1, 3
      next_edge = merge(edge + 1, 1, edge < 3)
      endpoint_a = geometry%vertex(:, edge)
      radius_a = radius(edge)
      radius_b = radius(next_edge)
      ratio = geometry%edge_length(edge)/(radius_a + radius_b)
      line_integral = panel_edge_line_integral(geometry, target, edge, ratio)
      distance = dot_product(endpoint_a - projection, geometry%edge_outward(:, edge))
      potential_integral = potential_integral + distance*line_integral
      field_integral = field_integral + geometry%edge_outward(:, edge)*line_integral
    end do

    potential = k_coulomb*charge/geometry%area*potential_integral
    field = k_coulomb*charge/geometry%area*field_integral
  end subroutine panel_potential_field

  !> Principal-value potential and field of `geometry` translated by `shift`, without copying it.
  !! The arithmetic is the same as `panel_potential_field` applied to a geometry whose vertices were
  !! translated, evaluated once per vertex distance instead of once per edge end, so the result is
  !! bitwise identical. A target on the panel plane falls back to that routine.
  subroutine panel_potential_field_translated(geometry, shift, charge, target, potential, field)
    type(panel_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: shift(3), charge, target(3)
    real(dp), intent(out) :: potential, field(3)
    real(dp) :: vertex(3, 3), relative_vertex(3, 3), radius(3)
    real(dp) :: projection(3), height, potential_integral, field_integral(3), omega
    real(dp) :: line_integral, distance, ratio, numerator, denominator, scale, plane_tolerance
    integer :: edge, next_edge

    potential = 0.0_dp
    field = 0.0_dp
    if (charge == 0.0_dp) return
    vertex = geometry%vertex + spread(shift, 2, 3)
    height = dot_product(target - vertex(:, 1), geometry%normal)
    projection = target - height*geometry%normal
    scale = maxval(geometry%edge_length)
    plane_tolerance = 64.0_dp*epsilon(1.0_dp)*scale
    if (abs(height) <= plane_tolerance) then
      call panel_potential_field_on_translated_surface(geometry, vertex, shift, charge, target, potential, field)
      return
    end if

    do edge = 1, 3
      relative_vertex(:, edge) = vertex(:, edge) - target
      radius(edge) = norm2(relative_vertex(:, edge))
    end do
    numerator = dot_product(relative_vertex(:, 1), &
                            cross_product(relative_vertex(:, 2), relative_vertex(:, 3)))
    denominator = radius(1)*radius(2)*radius(3) + &
                  dot_product(relative_vertex(:, 1), relative_vertex(:, 2))*radius(3) + &
                  dot_product(relative_vertex(:, 2), relative_vertex(:, 3))*radius(1) + &
                  dot_product(relative_vertex(:, 3), relative_vertex(:, 1))*radius(2)
    omega = -2.0_dp*atan2(numerator, denominator)

    potential_integral = -height*omega
    field_integral = geometry%normal*omega
    do edge = 1, 3
      next_edge = merge(edge + 1, 1, edge < 3)
      ratio = geometry%edge_length(edge)/(radius(edge) + radius(next_edge))
      if (ratio >= panel_edge_ratio_limit) then
        call panel_potential_field_on_translated_surface(geometry, vertex, shift, charge, target, potential, field)
        return
      end if
      line_integral = log((1.0_dp + ratio)/(1.0_dp - ratio))
      distance = dot_product(vertex(:, edge) - projection, geometry%edge_outward(:, edge))
      potential_integral = potential_integral + distance*line_integral
      field_integral = field_integral + geometry%edge_outward(:, edge)*line_integral
    end do

    potential = k_coulomb*charge/geometry%area*potential_integral
    field = k_coulomb*charge/geometry%area*field_integral
  end subroutine panel_potential_field_translated

  !> Principal-value field at `target` of `count` translated panels given as component arrays
  !! (translated vertex x/y/z per vertex, normal, edge lengths, outward edge normals, area), accumulated
  !! into `ex`, `ey`, `ez` divided by `k_coulomb` as the FMM near field stores it. Entry k carries
  !! `charge(index(k))`; `geometry(index(k))` is only read for the on-plane fallback. The SIMD loop
  !! performs the operations of `panel_potential_field_translated` in the same order and the
  !! contributions are added in list order, so the sum matches a loop over that routine bit for bit;
  !! zero charges are skipped as that loop skips them.
  subroutine panel_field_sum_translated( &
    count, vx, vy, vz, nx, ny, nz, edge, ox, oy, oz, area, index, charge, geometry, target, ex, ey, ez &
    )
    integer(i32), intent(in) :: count
    real(dp), intent(in) :: vx(3, *), vy(3, *), vz(3, *), nx(*), ny(*), nz(*), edge(3, *)
    real(dp), intent(in) :: ox(3, *), oy(3, *), oz(3, *), area(*)
    integer(i32), intent(in) :: index(*)
    real(dp), intent(in) :: charge(:)
    type(panel_geometry_type), intent(in) :: geometry(:)
    real(dp), intent(in) :: target(3)
    real(dp), intent(inout) :: ex, ey, ez
    real(dp) :: outx(panel_chunk), outy(panel_chunk), outz(panel_chunk)
    logical :: active(panel_chunk), fallback(panel_chunk)
    real(dp) :: ax, ay, az, bx, by, bz, cx, cy, cz, ra, rb, rc, height, numerator, denominator, omega
    real(dp) :: l1, l2, l3, fx, fy, fz, scale, potential, field(3)
    logical :: evaluate_lane
    integer(i32) :: chunk_start, chunk_size, k, m, idx

    do chunk_start = 1_i32, count, panel_chunk
      chunk_size = min(panel_chunk, count - chunk_start + 1_i32)
      !$omp simd private(m, ax, ay, az, bx, by, bz, cx, cy, cz, ra, rb, rc, height, numerator, denominator) &
      !$omp   private(omega, l1, l2, l3, fx, fy, fz, scale, evaluate_lane)
      do k = 1_i32, chunk_size
        m = chunk_start + k - 1_i32
        active(k) = abs(charge(index(m))) > tiny(1.0_dp)
        scale = k_coulomb*charge(index(m))/area(m)
        height = (target(1) - vx(1, m))*nx(m) + (target(2) - vy(1, m))*ny(m) + (target(3) - vz(1, m))*nz(m)
        fallback(k) = abs(height) <= 64.0_dp*epsilon(1.0_dp)*max(edge(1, m), edge(2, m), edge(3, m))
        ax = vx(1, m) - target(1)
        ay = vy(1, m) - target(2)
        az = vz(1, m) - target(3)
        bx = vx(2, m) - target(1)
        by = vy(2, m) - target(2)
        bz = vz(2, m) - target(3)
        cx = vx(3, m) - target(1)
        cy = vy(3, m) - target(2)
        cz = vz(3, m) - target(3)
        ra = sqrt(ax*ax + ay*ay + az*az)
        rb = sqrt(bx*bx + by*by + bz*bz)
        rc = sqrt(cx*cx + cy*cy + cz*cz)
        numerator = ax*(by*cz - bz*cy) + ay*(bz*cx - bx*cz) + az*(bx*cy - by*cx)
        denominator = ra*rb*rc + (ax*bx + ay*by + az*bz)*rc + (bx*cx + by*cy + bz*cz)*ra + &
                      (cx*ax + cy*ay + cz*az)*rb
        l1 = edge(1, m)/(ra + rb)
        l2 = edge(2, m)/(rb + rc)
        l3 = edge(3, m)/(rc + ra)
        fallback(k) = fallback(k) .or. max(l1, l2, l3) >= panel_edge_ratio_limit
        evaluate_lane = active(k) .and. .not. fallback(k)
        ! Near-edge lanes use the stable scalar integral; keep their vector math regular.
        omega = -2.0_dp*atan2(merge(numerator, 0.0_dp, evaluate_lane), merge(denominator, 1.0_dp, evaluate_lane))
        l1 = merge(l1, 0.0_dp, evaluate_lane)
        l2 = merge(l2, 0.0_dp, evaluate_lane)
        l3 = merge(l3, 0.0_dp, evaluate_lane)
        l1 = log((1.0_dp + l1)/(1.0_dp - l1))
        l2 = log((1.0_dp + l2)/(1.0_dp - l2))
        l3 = log((1.0_dp + l3)/(1.0_dp - l3))
        fx = nx(m)*omega
        fy = ny(m)*omega
        fz = nz(m)*omega
        fx = fx + ox(1, m)*l1
        fy = fy + oy(1, m)*l1
        fz = fz + oz(1, m)*l1
        fx = fx + ox(2, m)*l2
        fy = fy + oy(2, m)*l2
        fz = fz + oz(2, m)*l2
        fx = fx + ox(3, m)*l3
        fy = fy + oy(3, m)*l3
        fz = fz + oz(3, m)*l3
        outx(k) = scale*fx
        outy(k) = scale*fy
        outz(k) = scale*fz
      end do
      do k = 1_i32, chunk_size
        if (.not. active(k)) cycle
        if (fallback(k)) then
          m = chunk_start + k - 1_i32
          idx = index(m)
          call panel_potential_field_on_vertices( &
            geometry(idx), vx(:, m), vy(:, m), vz(:, m), charge(idx), target, potential, field &
            )
          outx(k) = field(1)
          outy(k) = field(2)
          outz(k) = field(3)
        end if
        ex = ex + outx(k)/k_coulomb
        ey = ey + outy(k)/k_coulomb
        ez = ez + outz(k)/k_coulomb
      end do
    end do
  end subroutine panel_field_sum_translated

  !> Principal-value potential counterpart of `panel_field_sum_translated`, accumulated into `phi`
  !! divided by `k_coulomb`.
  subroutine panel_potential_sum_translated( &
    count, vx, vy, vz, nx, ny, nz, edge, ox, oy, oz, area, index, charge, geometry, target, phi &
    )
    integer(i32), intent(in) :: count
    real(dp), intent(in) :: vx(3, *), vy(3, *), vz(3, *), nx(*), ny(*), nz(*), edge(3, *)
    real(dp), intent(in) :: ox(3, *), oy(3, *), oz(3, *), area(*)
    integer(i32), intent(in) :: index(*)
    real(dp), intent(in) :: charge(:)
    type(panel_geometry_type), intent(in) :: geometry(:)
    real(dp), intent(in) :: target(3)
    real(dp), intent(inout) :: phi
    real(dp) :: out(panel_chunk)
    logical :: active(panel_chunk), fallback(panel_chunk)
    real(dp) :: ax, ay, az, bx, by, bz, cx, cy, cz, ra, rb, rc, height, numerator, denominator, omega
    real(dp) :: px, py, pz, l1, l2, l3, d1, d2, d3, integral, scale, potential, field(3)
    logical :: evaluate_lane
    integer(i32) :: chunk_start, chunk_size, k, m, idx

    do chunk_start = 1_i32, count, panel_chunk
      chunk_size = min(panel_chunk, count - chunk_start + 1_i32)
      !$omp simd private(m, ax, ay, az, bx, by, bz, cx, cy, cz, ra, rb, rc, height, numerator, denominator) &
      !$omp   private(omega, px, py, pz, l1, l2, l3, d1, d2, d3, integral, scale, evaluate_lane)
      do k = 1_i32, chunk_size
        m = chunk_start + k - 1_i32
        active(k) = abs(charge(index(m))) > tiny(1.0_dp)
        scale = k_coulomb*charge(index(m))/area(m)
        height = (target(1) - vx(1, m))*nx(m) + (target(2) - vy(1, m))*ny(m) + (target(3) - vz(1, m))*nz(m)
        px = target(1) - height*nx(m)
        py = target(2) - height*ny(m)
        pz = target(3) - height*nz(m)
        fallback(k) = abs(height) <= 64.0_dp*epsilon(1.0_dp)*max(edge(1, m), edge(2, m), edge(3, m))
        ax = vx(1, m) - target(1)
        ay = vy(1, m) - target(2)
        az = vz(1, m) - target(3)
        bx = vx(2, m) - target(1)
        by = vy(2, m) - target(2)
        bz = vz(2, m) - target(3)
        cx = vx(3, m) - target(1)
        cy = vy(3, m) - target(2)
        cz = vz(3, m) - target(3)
        ra = sqrt(ax*ax + ay*ay + az*az)
        rb = sqrt(bx*bx + by*by + bz*bz)
        rc = sqrt(cx*cx + cy*cy + cz*cz)
        numerator = ax*(by*cz - bz*cy) + ay*(bz*cx - bx*cz) + az*(bx*cy - by*cx)
        denominator = ra*rb*rc + (ax*bx + ay*by + az*bz)*rc + (bx*cx + by*cy + bz*cz)*ra + &
                      (cx*ax + cy*ay + cz*az)*rb
        l1 = edge(1, m)/(ra + rb)
        l2 = edge(2, m)/(rb + rc)
        l3 = edge(3, m)/(rc + ra)
        fallback(k) = fallback(k) .or. max(l1, l2, l3) >= panel_edge_ratio_limit
        evaluate_lane = active(k) .and. .not. fallback(k)
        omega = -2.0_dp*atan2(merge(numerator, 0.0_dp, evaluate_lane), merge(denominator, 1.0_dp, evaluate_lane))
        l1 = merge(l1, 0.0_dp, evaluate_lane)
        l2 = merge(l2, 0.0_dp, evaluate_lane)
        l3 = merge(l3, 0.0_dp, evaluate_lane)
        l1 = log((1.0_dp + l1)/(1.0_dp - l1))
        l2 = log((1.0_dp + l2)/(1.0_dp - l2))
        l3 = log((1.0_dp + l3)/(1.0_dp - l3))
        d1 = (vx(1, m) - px)*ox(1, m) + (vy(1, m) - py)*oy(1, m) + (vz(1, m) - pz)*oz(1, m)
        d2 = (vx(2, m) - px)*ox(2, m) + (vy(2, m) - py)*oy(2, m) + (vz(2, m) - pz)*oz(2, m)
        d3 = (vx(3, m) - px)*ox(3, m) + (vy(3, m) - py)*oy(3, m) + (vz(3, m) - pz)*oz(3, m)
        integral = -height*omega
        integral = integral + d1*l1
        integral = integral + d2*l2
        integral = integral + d3*l3
        out(k) = scale*integral
      end do
      do k = 1_i32, chunk_size
        if (.not. active(k)) cycle
        if (fallback(k)) then
          m = chunk_start + k - 1_i32
          idx = index(m)
          call panel_potential_field_on_vertices( &
            geometry(idx), vx(:, m), vy(:, m), vz(:, m), charge(idx), target, potential, field &
            )
          out(k) = potential
        end if
        phi = phi + out(k)/k_coulomb
      end do
    end do
  end subroutine panel_potential_sum_translated

  !> Single-precision counterpart of `panel_field_sum_translated`. The component tables hold the panels
  !! translated by their image shift and relative to the leaf center; `target_local` is the target in the
  !! same frame (double, rounded here) and `target` the absolute target for the on-plane fallback, which
  !! still runs the double-precision kernel on `geometry(index(k))` translated by `shift1(k)`, `shift2(k)`
  !! along `axes`. Eight panels per AVX2 vector; the field error is around 1e-6 relative.
  subroutine panel_field_sum_translated_sp( &
    count, vx, vy, vz, nx, ny, nz, edge, ox, oy, oz, area, index, charge, geometry, axes, shift1, shift2, &
    target_local, target, ex, ey, ez &
    )
    integer(i32), intent(in) :: count
    real(sp), intent(in) :: vx(3, *), vy(3, *), vz(3, *), nx(*), ny(*), nz(*), edge(3, *)
    real(sp), intent(in) :: ox(3, *), oy(3, *), oz(3, *), area(*)
    integer(i32), intent(in) :: index(*), axes(2)
    real(dp), intent(in) :: charge(:), shift1(*), shift2(*)
    type(panel_geometry_type), intent(in) :: geometry(:)
    real(dp), intent(in) :: target_local(3), target(3)
    real(dp), intent(inout) :: ex, ey, ez
    real(sp) :: outx(panel_chunk), outy(panel_chunk), outz(panel_chunk), tx, ty, tz
    logical :: active(panel_chunk), fallback(panel_chunk)
    real(sp) :: ax, ay, az, bx, by, bz, cx, cy, cz, ra, rb, rc, height, numerator, denominator, omega
    real(sp) :: l1, l2, l3, fx, fy, fz, scale
    real(dp) :: potential, field(3), shift(3)
    logical :: evaluate_lane
    integer(i32) :: chunk_start, chunk_size, k, m, idx

    tx = real(target_local(1), sp)
    ty = real(target_local(2), sp)
    tz = real(target_local(3), sp)
    do chunk_start = 1_i32, count, panel_chunk
      chunk_size = min(panel_chunk, count - chunk_start + 1_i32)
      !$omp simd private(m, ax, ay, az, bx, by, bz, cx, cy, cz, ra, rb, rc, height, numerator, denominator) &
      !$omp   private(omega, l1, l2, l3, fx, fy, fz, scale, evaluate_lane)
      do k = 1_i32, chunk_size
        m = chunk_start + k - 1_i32
        active(k) = abs(charge(index(m))) > tiny(1.0_dp)
        scale = real(k_coulomb*charge(index(m)), sp)/area(m)
        height = (tx - vx(1, m))*nx(m) + (ty - vy(1, m))*ny(m) + (tz - vz(1, m))*nz(m)
        ax = vx(1, m) - tx
        ay = vy(1, m) - ty
        az = vz(1, m) - tz
        bx = vx(2, m) - tx
        by = vy(2, m) - ty
        bz = vz(2, m) - tz
        cx = vx(3, m) - tx
        cy = vy(3, m) - ty
        cz = vz(3, m) - tz
        ra = sqrt(ax*ax + ay*ay + az*az)
        rb = sqrt(bx*bx + by*by + bz*bz)
        rc = sqrt(cx*cx + cy*cy + cz*cz)
        numerator = ax*(by*cz - bz*cy) + ay*(bz*cx - bx*cz) + az*(bx*cy - by*cx)
        denominator = ra*rb*rc + (ax*bx + ay*by + az*bz)*rc + (bx*cx + by*cy + bz*cz)*ra + &
                      (cx*ax + cy*ay + cz*az)*rb
        l1 = edge(1, m)/(ra + rb)
        l2 = edge(2, m)/(rb + rc)
        l3 = edge(3, m)/(rc + ra)
        ! Points on the panel plane or within about one percent of an edge lose single-precision
        ! digits in log((1+r)/(1-r)); those lanes are redone with the double kernel below.
        fallback(k) = active(k) .and. (abs(height) <= 64.0_sp*epsilon(1.0_sp)*max(edge(1, m), edge(2, m), edge(3, m)) &
                                       .or. max(l1, l2, l3) >= near_edge_single_limit)
        evaluate_lane = active(k) .and. .not. fallback(k)
        omega = -2.0_sp*vector_atan2_sp(merge(numerator, 0.0_sp, evaluate_lane), merge(denominator, 1.0_sp, evaluate_lane))
        l1 = merge(l1, 0.0_sp, evaluate_lane)
        l2 = merge(l2, 0.0_sp, evaluate_lane)
        l3 = merge(l3, 0.0_sp, evaluate_lane)
        l1 = vector_log_sp((1.0_sp + l1)/(1.0_sp - l1))
        l2 = vector_log_sp((1.0_sp + l2)/(1.0_sp - l2))
        l3 = vector_log_sp((1.0_sp + l3)/(1.0_sp - l3))
        fx = nx(m)*omega + ox(1, m)*l1 + ox(2, m)*l2 + ox(3, m)*l3
        fy = ny(m)*omega + oy(1, m)*l1 + oy(2, m)*l2 + oy(3, m)*l3
        fz = nz(m)*omega + oz(1, m)*l1 + oz(2, m)*l2 + oz(3, m)*l3
        outx(k) = scale*fx
        outy(k) = scale*fy
        outz(k) = scale*fz
      end do
      do k = 1_i32, chunk_size
        if (.not. active(k)) cycle
        if (fallback(k)) then
          m = chunk_start + k - 1_i32
          idx = index(m)
          shift = 0.0_dp
          if (axes(1) > 0_i32) shift(axes(1)) = shift(axes(1)) + shift1(m)
          if (axes(2) > 0_i32) shift(axes(2)) = shift(axes(2)) + shift2(m)
          call panel_potential_field_translated(geometry(idx), shift, charge(idx), target, potential, field)
          ex = ex + field(1)/k_coulomb
          ey = ey + field(2)/k_coulomb
          ez = ez + field(3)/k_coulomb
          cycle
        end if
        ex = ex + real(outx(k), dp)/k_coulomb
        ey = ey + real(outy(k), dp)/k_coulomb
        ez = ez + real(outz(k), dp)/k_coulomb
      end do
    end do
  end subroutine panel_field_sum_translated_sp

  !> Single-precision counterpart of `panel_potential_sum_translated`; see `panel_field_sum_translated_sp`.
  subroutine panel_potential_sum_translated_sp( &
    count, vx, vy, vz, nx, ny, nz, edge, ox, oy, oz, area, index, charge, geometry, axes, shift1, shift2, &
    target_local, target, phi &
    )
    integer(i32), intent(in) :: count
    real(sp), intent(in) :: vx(3, *), vy(3, *), vz(3, *), nx(*), ny(*), nz(*), edge(3, *)
    real(sp), intent(in) :: ox(3, *), oy(3, *), oz(3, *), area(*)
    integer(i32), intent(in) :: index(*), axes(2)
    real(dp), intent(in) :: charge(:), shift1(*), shift2(*)
    type(panel_geometry_type), intent(in) :: geometry(:)
    real(dp), intent(in) :: target_local(3), target(3)
    real(dp), intent(inout) :: phi
    real(sp) :: out(panel_chunk), tx, ty, tz
    logical :: active(panel_chunk), fallback(panel_chunk)
    real(sp) :: ax, ay, az, bx, by, bz, cx, cy, cz, ra, rb, rc, height, numerator, denominator, omega
    real(sp) :: px, py, pz, l1, l2, l3, d1, d2, d3, integral, scale
    real(dp) :: potential, field(3), shift(3)
    logical :: evaluate_lane
    integer(i32) :: chunk_start, chunk_size, k, m, idx

    tx = real(target_local(1), sp)
    ty = real(target_local(2), sp)
    tz = real(target_local(3), sp)
    do chunk_start = 1_i32, count, panel_chunk
      chunk_size = min(panel_chunk, count - chunk_start + 1_i32)
      !$omp simd private(m, ax, ay, az, bx, by, bz, cx, cy, cz, ra, rb, rc, height, numerator, denominator) &
      !$omp   private(omega, px, py, pz, l1, l2, l3, d1, d2, d3, integral, scale, evaluate_lane)
      do k = 1_i32, chunk_size
        m = chunk_start + k - 1_i32
        active(k) = abs(charge(index(m))) > tiny(1.0_dp)
        scale = real(k_coulomb*charge(index(m)), sp)/area(m)
        height = (tx - vx(1, m))*nx(m) + (ty - vy(1, m))*ny(m) + (tz - vz(1, m))*nz(m)
        px = tx - height*nx(m)
        py = ty - height*ny(m)
        pz = tz - height*nz(m)
        ax = vx(1, m) - tx
        ay = vy(1, m) - ty
        az = vz(1, m) - tz
        bx = vx(2, m) - tx
        by = vy(2, m) - ty
        bz = vz(2, m) - tz
        cx = vx(3, m) - tx
        cy = vy(3, m) - ty
        cz = vz(3, m) - tz
        ra = sqrt(ax*ax + ay*ay + az*az)
        rb = sqrt(bx*bx + by*by + bz*bz)
        rc = sqrt(cx*cx + cy*cy + cz*cz)
        numerator = ax*(by*cz - bz*cy) + ay*(bz*cx - bx*cz) + az*(bx*cy - by*cx)
        denominator = ra*rb*rc + (ax*bx + ay*by + az*bz)*rc + (bx*cx + by*cy + bz*cz)*ra + &
                      (cx*ax + cy*ay + cz*az)*rb
        l1 = edge(1, m)/(ra + rb)
        l2 = edge(2, m)/(rb + rc)
        l3 = edge(3, m)/(rc + ra)
        fallback(k) = active(k) .and. (abs(height) <= 64.0_sp*epsilon(1.0_sp)*max(edge(1, m), edge(2, m), edge(3, m)) &
                                       .or. max(l1, l2, l3) >= near_edge_single_limit)
        evaluate_lane = active(k) .and. .not. fallback(k)
        omega = -2.0_sp*vector_atan2_sp(merge(numerator, 0.0_sp, evaluate_lane), merge(denominator, 1.0_sp, evaluate_lane))
        l1 = merge(l1, 0.0_sp, evaluate_lane)
        l2 = merge(l2, 0.0_sp, evaluate_lane)
        l3 = merge(l3, 0.0_sp, evaluate_lane)
        l1 = vector_log_sp((1.0_sp + l1)/(1.0_sp - l1))
        l2 = vector_log_sp((1.0_sp + l2)/(1.0_sp - l2))
        l3 = vector_log_sp((1.0_sp + l3)/(1.0_sp - l3))
        d1 = (vx(1, m) - px)*ox(1, m) + (vy(1, m) - py)*oy(1, m) + (vz(1, m) - pz)*oz(1, m)
        d2 = (vx(2, m) - px)*ox(2, m) + (vy(2, m) - py)*oy(2, m) + (vz(2, m) - pz)*oz(2, m)
        d3 = (vx(3, m) - px)*ox(3, m) + (vy(3, m) - py)*oy(3, m) + (vz(3, m) - pz)*oz(3, m)
        integral = -height*omega + d1*l1 + d2*l2 + d3*l3
        out(k) = scale*integral
      end do
      do k = 1_i32, chunk_size
        if (.not. active(k)) cycle
        if (fallback(k)) then
          m = chunk_start + k - 1_i32
          idx = index(m)
          shift = 0.0_dp
          if (axes(1) > 0_i32) shift(axes(1)) = shift(axes(1)) + shift1(m)
          if (axes(2) > 0_i32) shift(axes(2)) = shift(axes(2)) + shift2(m)
          call panel_potential_field_translated(geometry(idx), shift, charge(idx), target, potential, field)
          phi = phi + potential/k_coulomb
          cycle
        end if
        phi = phi + real(out(k), dp)/k_coulomb
      end do
    end do
  end subroutine panel_potential_sum_translated_sp

  !> Natural logarithm of a positive normal single-precision `x`, written with arithmetic and bit operations only
  !! so the single-precision near-field SIMD loops vectorize it: under the precise floating-point model the
  !! compilers call the scalar `logf` once per lane instead. Cephes `logf` reduction and polynomial; the
  !! relative error stays below 2e-7 (tested on 1 <= x <= 1e6).
  elemental real(sp) function vector_log_sp(x) result(y)
    !$omp declare simd(vector_log_sp)
    real(sp), intent(in) :: x
    integer(i32) :: bits, exponent_value
    real(sp) :: mantissa, f, z, p, e
    logical :: upper

    bits = transfer(x, bits)
    exponent_value = ishft(bits, -23) - 127_i32
    ! Mantissa in [1, 2), then folded into [sqrt(1/2), sqrt(2)) so f = mantissa - 1 stays small.
    mantissa = transfer(ior(iand(bits, 8388607_i32), 1065353216_i32), mantissa)
    upper = mantissa > 1.41421356_sp
    mantissa = merge(0.5_sp*mantissa, mantissa, upper)
    exponent_value = merge(exponent_value + 1_i32, exponent_value, upper)
    f = mantissa - 1.0_sp
    z = f*f
    p = 7.0376836292e-2_sp
    p = p*f - 1.1514610310e-1_sp
    p = p*f + 1.1676998740e-1_sp
    p = p*f - 1.2420140846e-1_sp
    p = p*f + 1.4249322787e-1_sp
    p = p*f - 1.6668057665e-1_sp
    p = p*f + 2.0000714765e-1_sp
    p = p*f - 2.4999993993e-1_sp
    p = p*f + 3.3333331174e-1_sp
    e = real(exponent_value, sp)
    ! log(2) is split into 0.693359375 (exact in a few bits) and -2.12194440e-4 to keep e*log(2) exact.
    y = f*z*p
    y = y - 2.12194440e-4_sp*e
    y = y - 0.5_sp*z
    y = y + f
    y = y + 0.693359375_sp*e
  end function vector_log_sp

  !> `atan2(y, x)` in single precision with arithmetic and selects only, for the same reason as
  !! `vector_log_sp`. The ratio min/max of |y| and |x| is reduced to |t| <= tan(pi/8) and evaluated
  !! with the Cephes `atanf` polynomial; the quadrant follows from the signs and which of |y|, |x| is
  !! larger. Relative error below 3e-7 (tested over the full circle). For y = x = 0 it returns 0.
  elemental real(sp) function vector_atan2_sp(y, x) result(angle)
    !$omp declare simd(vector_atan2_sp)
    real(sp), intent(in) :: y, x
    real(sp), parameter :: pi = 3.14159265358979324_sp
    real(sp) :: ay, ax, t, z, offset
    logical :: reduced

    ay = abs(y)
    ax = abs(x)
    t = min(ay, ax)/max(ay, ax, tiny(1.0_sp))
    reduced = t > 0.414213562_sp
    offset = merge(0.25_sp*pi, 0.0_sp, reduced)
    t = merge((t - 1.0_sp)/(t + 1.0_sp), t, reduced)
    z = t*t
    angle = (((8.05374449538e-2_sp*z - 1.38776856032e-1_sp)*z + 1.99777106478e-1_sp)*z - 3.33329491539e-1_sp)*z*t + t
    angle = angle + offset
    angle = merge(0.5_sp*pi - angle, angle, ay > ax)
    angle = merge(pi - angle, angle, x < 0.0_sp)
    angle = merge(-angle, angle, y < 0.0_sp)
  end function vector_atan2_sp

  !> On-plane fallback of the component-array sums: the general kernel on a copy of `geometry` whose
  !! vertices are the translated ones. The kernel never reads the centroid, so it is left untranslated.
  subroutine panel_potential_field_on_vertices(geometry, vx, vy, vz, charge, target, potential, field)
    type(panel_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: vx(3), vy(3), vz(3), charge, target(3)
    real(dp), intent(out) :: potential, field(3)
    type(panel_geometry_type) :: translated

    translated = geometry
    translated%vertex(1, :) = vx
    translated%vertex(2, :) = vy
    translated%vertex(3, :) = vz
    call panel_potential_field(translated, charge, target, panel_side_principal_value, potential, field)
  end subroutine panel_potential_field_on_vertices

  !> On-plane fallback of `panel_potential_field_translated`; kept apart so the hot path never
  !! default-initializes a geometry copy.
  subroutine panel_potential_field_on_translated_surface(geometry, vertex, shift, charge, target, potential, field)
    type(panel_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: vertex(3, 3), shift(3), charge, target(3)
    real(dp), intent(out) :: potential, field(3)
    type(panel_geometry_type) :: translated

    translated = geometry
    translated%vertex = vertex
    translated%centroid = geometry%centroid + shift
    call panel_potential_field(translated, charge, target, panel_side_principal_value, potential, field)
  end subroutine panel_potential_field_on_translated_surface

  !> Sum the oriented solid angles of the triangles formed by the target's projection
  !! and each edge. Unlike the three-vertex formula, their denominators do not subtract
  !! nearly equal O(L**3) terms as the target approaches an edge at finite height.
  real(dp) function near_edge_solid_angle(geometry, target, projection, height) result(omega)
    type(panel_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: target(3), projection(3), height
    real(dp) :: a(3), b(3), crossed(3), ra, rb, ab, radial_sum, oriented_area, denominator
    integer :: edge, next_edge

    omega = 0.0_dp
    do edge = 1, 3
      next_edge = merge(edge + 1, 1, edge < 3)
      a = geometry%vertex(:, edge) - target
      b = geometry%vertex(:, next_edge) - target
      ra = norm2(a)
      rb = norm2(b)
      ab = dot_product(a, b)
      if (ab < 0.0_dp) then
        crossed = cross_product(a, b)
        radial_sum = sum(crossed*crossed)/(ra*rb - ab)
      else
        radial_sum = ra*rb + ab
      end if
      oriented_area = geometry%edge_length(edge)* &
                      dot_product(geometry%vertex(:, edge) - projection, geometry%edge_outward(:, edge))
      denominator = radial_sum + abs(height)*(ra + rb)
      omega = omega + 2.0_dp*atan2(sign(1.0_dp, height)*oriented_area, denominator)
    end do
  end function near_edge_solid_angle

  pure logical function point_in_triangle(geometry, point) result(inside)
    type(panel_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: point(3)
    real(dp) :: edge_vector(3), point_vector(3), signed_edge, tolerance
    integer :: edge, next_edge

    tolerance = 128.0_dp*epsilon(1.0_dp)*maxval(geometry%edge_length)**2
    inside = .true.
    do edge = 1, 3
      next_edge = merge(edge + 1, 1, edge < 3)
      edge_vector = geometry%vertex(:, next_edge) - geometry%vertex(:, edge)
      point_vector = point - geometry%vertex(:, edge)
      signed_edge = dot_product(cross_product(edge_vector, point_vector), geometry%normal)
      if (signed_edge < -tolerance) then
        inside = .false.
        return
      end if
    end do
  end function point_in_triangle

  pure real(dp) function norm2(vector) result(norm)
    real(dp), intent(in) :: vector(3)

    norm = sqrt(sum(vector*vector))
  end function norm2

  pure function cross_product(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)

    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  end function cross_product

end module bem_panel_kernel

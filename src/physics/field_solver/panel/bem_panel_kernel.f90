!> Analytic free-space P0 triangle potential, field, principal value, and jump.
module bem_panel_kernel
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0, k_coulomb
  use bem_panel_geometry, only: panel_geometry_type
  use bem_panel_self_terms, only: panel_on_surface_integrals
  implicit none
  private

  integer(i32), parameter, public :: panel_side_principal_value = 0_i32
  integer(i32), parameter, public :: panel_side_normal_plus = 1_i32
  integer(i32), parameter, public :: panel_side_normal_minus = -1_i32

  public :: panel_potential_field
  public :: panel_potential_field_translated
  public :: panel_field_sum_translated
  public :: panel_potential_sum_translated

  !> Panels gathered per SIMD pass of the near-field sums.
  integer(i32), parameter :: panel_chunk = 64_i32

contains

  subroutine panel_potential_field(geometry, charge, target, side, potential, field)
    type(panel_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: charge, target(3)
    integer(i32), intent(in) :: side
    real(dp), intent(out) :: potential, field(3)
    real(dp) :: projection(3), height, potential_integral, field_integral(3), omega
    real(dp) :: endpoint_a(3), endpoint_b(3), radius_a, radius_b, line_integral, distance, ratio
    real(dp) :: relative_vertex(3, 3), numerator, denominator, scale, plane_tolerance
    integer :: edge, next_edge
    logical :: on_surface

    if (side < panel_side_normal_minus .or. side > panel_side_normal_plus) then
      error stop 'invalid panel evaluation side.'
    end if
    height = dot_product(target - geometry%vertex(:, 1), geometry%normal)
    projection = target - height*geometry%normal
    scale = maxval(geometry%edge_length)
    plane_tolerance = 64.0_dp*epsilon(1.0_dp)*scale
    on_surface = abs(height) <= plane_tolerance

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
    end do
    numerator = dot_product(relative_vertex(:, 1), &
                            cross_product(relative_vertex(:, 2), relative_vertex(:, 3)))
    denominator = product_norms(relative_vertex) + &
                  dot_product(relative_vertex(:, 1), relative_vertex(:, 2))*norm2(relative_vertex(:, 3)) + &
                  dot_product(relative_vertex(:, 2), relative_vertex(:, 3))*norm2(relative_vertex(:, 1)) + &
                  dot_product(relative_vertex(:, 3), relative_vertex(:, 1))*norm2(relative_vertex(:, 2))
    omega = -2.0_dp*atan2(numerator, denominator)

    potential_integral = -height*omega
    field_integral = geometry%normal*omega
    do edge = 1, 3
      next_edge = merge(edge + 1, 1, edge < 3)
      endpoint_a = geometry%vertex(:, edge)
      endpoint_b = geometry%vertex(:, next_edge)
      radius_a = norm2(target - endpoint_a)
      radius_b = norm2(target - endpoint_b)
      ratio = geometry%edge_length(edge)/(radius_a + radius_b)
      if (ratio >= 1.0_dp) error stop 'panel target lies on an edge or vertex.'
      line_integral = log((1.0_dp + ratio)/(1.0_dp - ratio))
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
      if (ratio >= 1.0_dp) error stop 'panel target lies on an edge or vertex.'
      line_integral = log((1.0_dp + ratio)/(1.0_dp - ratio))
      distance = dot_product(vertex(:, edge) - projection, geometry%edge_outward(:, edge))
      potential_integral = potential_integral + distance*line_integral
      field_integral = field_integral + geometry%edge_outward(:, edge)*line_integral
    end do

    potential = k_coulomb*charge/geometry%area*potential_integral
    field = k_coulomb*charge/geometry%area*field_integral
  end subroutine panel_potential_field_translated

  !> Principal-value field at `target` of `count` translated panels, accumulated into `ex`, `ey`, `ez`
  !! divided by `k_coulomb` as the FMM near field stores it. Entry k is panel `index(k)` translated by
  !! `shift1(k)` and `shift2(k)` along `axes` and carrying `charge(index(k))`. The panels are gathered
  !! into chunks of component arrays and evaluated in a SIMD loop with the operations of
  !! `panel_potential_field_translated` in the same order, then accumulated in list order, so the sum
  !! matches a loop over that routine bit for bit. Zero charges are skipped and on-plane targets use
  !! the general kernel, as the sequential loop does.
  subroutine panel_field_sum_translated(geometry, index, axes, shift1, shift2, charge, count, target, ex, ey, ez)
    type(panel_geometry_type), intent(in) :: geometry(:)
    integer(i32), intent(in) :: index(:), axes(2)
    real(dp), intent(in) :: shift1(:), shift2(:), charge(:)
    integer(i32), intent(in) :: count
    real(dp), intent(in) :: target(3)
    real(dp), intent(inout) :: ex, ey, ez
    real(dp) :: vx(3, panel_chunk), vy(3, panel_chunk), vz(3, panel_chunk)
    real(dp) :: nx(panel_chunk), ny(panel_chunk), nz(panel_chunk), edge(3, panel_chunk)
    real(dp) :: ox(3, panel_chunk), oy(3, panel_chunk), oz(3, panel_chunk), scale(panel_chunk)
    real(dp) :: outx(panel_chunk), outy(panel_chunk), outz(panel_chunk), shift(3, panel_chunk)
    logical :: active(panel_chunk), on_plane(panel_chunk)
    real(dp) :: ax, ay, az, bx, by, bz, cx, cy, cz, ra, rb, rc, height, numerator, denominator, omega
    real(dp) :: l1, l2, l3, fx, fy, fz, potential, field(3)
    logical :: edge_hit
    integer(i32) :: chunk_start, chunk_size, k, idx

    do chunk_start = 1_i32, count, panel_chunk
      chunk_size = min(panel_chunk, count - chunk_start + 1_i32)
      call gather_translated_panels( &
        geometry, index(chunk_start:), axes, shift1(chunk_start:), shift2(chunk_start:), charge, chunk_size, &
        vx, vy, vz, nx, ny, nz, edge, ox, oy, oz, scale, shift, active &
        )
      edge_hit = .false.
      !$omp simd private(ax, ay, az, bx, by, bz, cx, cy, cz, ra, rb, rc, height, numerator, denominator) &
      !$omp   private(omega, l1, l2, l3, fx, fy, fz) reduction(.or.:edge_hit)
      do k = 1_i32, chunk_size
        height = (target(1) - vx(1, k))*nx(k) + (target(2) - vy(1, k))*ny(k) + (target(3) - vz(1, k))*nz(k)
        on_plane(k) = abs(height) <= 64.0_dp*epsilon(1.0_dp)*max(edge(1, k), edge(2, k), edge(3, k))
        ax = vx(1, k) - target(1)
        ay = vy(1, k) - target(2)
        az = vz(1, k) - target(3)
        bx = vx(2, k) - target(1)
        by = vy(2, k) - target(2)
        bz = vz(2, k) - target(3)
        cx = vx(3, k) - target(1)
        cy = vy(3, k) - target(2)
        cz = vz(3, k) - target(3)
        ra = sqrt(ax*ax + ay*ay + az*az)
        rb = sqrt(bx*bx + by*by + bz*bz)
        rc = sqrt(cx*cx + cy*cy + cz*cz)
        numerator = ax*(by*cz - bz*cy) + ay*(bz*cx - bx*cz) + az*(bx*cy - by*cx)
        denominator = ra*rb*rc + (ax*bx + ay*by + az*bz)*rc + (bx*cx + by*cy + bz*cz)*ra + &
                      (cx*ax + cy*ay + cz*az)*rb
        omega = -2.0_dp*atan2(numerator, denominator)
        l1 = edge(1, k)/(ra + rb)
        l2 = edge(2, k)/(rb + rc)
        l3 = edge(3, k)/(rc + ra)
        edge_hit = edge_hit .or. (active(k) .and. .not. on_plane(k) .and. max(l1, l2, l3) >= 1.0_dp)
        l1 = log((1.0_dp + l1)/(1.0_dp - l1))
        l2 = log((1.0_dp + l2)/(1.0_dp - l2))
        l3 = log((1.0_dp + l3)/(1.0_dp - l3))
        fx = nx(k)*omega
        fy = ny(k)*omega
        fz = nz(k)*omega
        fx = fx + ox(1, k)*l1
        fy = fy + oy(1, k)*l1
        fz = fz + oz(1, k)*l1
        fx = fx + ox(2, k)*l2
        fy = fy + oy(2, k)*l2
        fz = fz + oz(2, k)*l2
        fx = fx + ox(3, k)*l3
        fy = fy + oy(3, k)*l3
        fz = fz + oz(3, k)*l3
        outx(k) = scale(k)*fx
        outy(k) = scale(k)*fy
        outz(k) = scale(k)*fz
      end do
      if (edge_hit) error stop 'panel target lies on an edge or vertex.'
      do k = 1_i32, chunk_size
        if (.not. active(k)) cycle
        if (on_plane(k)) then
          idx = index(chunk_start + k - 1_i32)
          call panel_potential_field_on_translated_surface( &
            geometry(idx), geometry(idx)%vertex + spread(shift(:, k), 2, 3), shift(:, k), charge(idx), target, &
            potential, field &
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
  subroutine panel_potential_sum_translated(geometry, index, axes, shift1, shift2, charge, count, target, phi)
    type(panel_geometry_type), intent(in) :: geometry(:)
    integer(i32), intent(in) :: index(:), axes(2)
    real(dp), intent(in) :: shift1(:), shift2(:), charge(:)
    integer(i32), intent(in) :: count
    real(dp), intent(in) :: target(3)
    real(dp), intent(inout) :: phi
    real(dp) :: vx(3, panel_chunk), vy(3, panel_chunk), vz(3, panel_chunk)
    real(dp) :: nx(panel_chunk), ny(panel_chunk), nz(panel_chunk), edge(3, panel_chunk)
    real(dp) :: ox(3, panel_chunk), oy(3, panel_chunk), oz(3, panel_chunk), scale(panel_chunk)
    real(dp) :: out(panel_chunk), shift(3, panel_chunk)
    logical :: active(panel_chunk), on_plane(panel_chunk)
    real(dp) :: ax, ay, az, bx, by, bz, cx, cy, cz, ra, rb, rc, height, numerator, denominator, omega
    real(dp) :: px, py, pz, l1, l2, l3, d1, d2, d3, integral, potential, field(3)
    logical :: edge_hit
    integer(i32) :: chunk_start, chunk_size, k, idx

    do chunk_start = 1_i32, count, panel_chunk
      chunk_size = min(panel_chunk, count - chunk_start + 1_i32)
      call gather_translated_panels( &
        geometry, index(chunk_start:), axes, shift1(chunk_start:), shift2(chunk_start:), charge, chunk_size, &
        vx, vy, vz, nx, ny, nz, edge, ox, oy, oz, scale, shift, active &
        )
      edge_hit = .false.
      !$omp simd private(ax, ay, az, bx, by, bz, cx, cy, cz, ra, rb, rc, height, numerator, denominator) &
      !$omp   private(omega, px, py, pz, l1, l2, l3, d1, d2, d3, integral) reduction(.or.:edge_hit)
      do k = 1_i32, chunk_size
        height = (target(1) - vx(1, k))*nx(k) + (target(2) - vy(1, k))*ny(k) + (target(3) - vz(1, k))*nz(k)
        px = target(1) - height*nx(k)
        py = target(2) - height*ny(k)
        pz = target(3) - height*nz(k)
        on_plane(k) = abs(height) <= 64.0_dp*epsilon(1.0_dp)*max(edge(1, k), edge(2, k), edge(3, k))
        ax = vx(1, k) - target(1)
        ay = vy(1, k) - target(2)
        az = vz(1, k) - target(3)
        bx = vx(2, k) - target(1)
        by = vy(2, k) - target(2)
        bz = vz(2, k) - target(3)
        cx = vx(3, k) - target(1)
        cy = vy(3, k) - target(2)
        cz = vz(3, k) - target(3)
        ra = sqrt(ax*ax + ay*ay + az*az)
        rb = sqrt(bx*bx + by*by + bz*bz)
        rc = sqrt(cx*cx + cy*cy + cz*cz)
        numerator = ax*(by*cz - bz*cy) + ay*(bz*cx - bx*cz) + az*(bx*cy - by*cx)
        denominator = ra*rb*rc + (ax*bx + ay*by + az*bz)*rc + (bx*cx + by*cy + bz*cz)*ra + &
                      (cx*ax + cy*ay + cz*az)*rb
        omega = -2.0_dp*atan2(numerator, denominator)
        l1 = edge(1, k)/(ra + rb)
        l2 = edge(2, k)/(rb + rc)
        l3 = edge(3, k)/(rc + ra)
        edge_hit = edge_hit .or. (active(k) .and. .not. on_plane(k) .and. max(l1, l2, l3) >= 1.0_dp)
        l1 = log((1.0_dp + l1)/(1.0_dp - l1))
        l2 = log((1.0_dp + l2)/(1.0_dp - l2))
        l3 = log((1.0_dp + l3)/(1.0_dp - l3))
        d1 = (vx(1, k) - px)*ox(1, k) + (vy(1, k) - py)*oy(1, k) + (vz(1, k) - pz)*oz(1, k)
        d2 = (vx(2, k) - px)*ox(2, k) + (vy(2, k) - py)*oy(2, k) + (vz(2, k) - pz)*oz(2, k)
        d3 = (vx(3, k) - px)*ox(3, k) + (vy(3, k) - py)*oy(3, k) + (vz(3, k) - pz)*oz(3, k)
        integral = -height*omega
        integral = integral + d1*l1
        integral = integral + d2*l2
        integral = integral + d3*l3
        out(k) = scale(k)*integral
      end do
      if (edge_hit) error stop 'panel target lies on an edge or vertex.'
      do k = 1_i32, chunk_size
        if (.not. active(k)) cycle
        if (on_plane(k)) then
          idx = index(chunk_start + k - 1_i32)
          call panel_potential_field_on_translated_surface( &
            geometry(idx), geometry(idx)%vertex + spread(shift(:, k), 2, 3), shift(:, k), charge(idx), target, &
            potential, field &
            )
          out(k) = potential
        end if
        phi = phi + out(k)/k_coulomb
      end do
    end do
  end subroutine panel_potential_sum_translated

  !> Copy one chunk of translated panels into component arrays for the SIMD sums.
  subroutine gather_translated_panels( &
    geometry, index, axes, shift1, shift2, charge, chunk_size, vx, vy, vz, nx, ny, nz, edge, ox, oy, oz, scale, &
    shift, active &
    )
    type(panel_geometry_type), intent(in) :: geometry(:)
    integer(i32), intent(in) :: index(:), axes(2), chunk_size
    real(dp), intent(in) :: shift1(:), shift2(:), charge(:)
    real(dp), intent(out) :: vx(3, panel_chunk), vy(3, panel_chunk), vz(3, panel_chunk)
    real(dp), intent(out) :: nx(panel_chunk), ny(panel_chunk), nz(panel_chunk), edge(3, panel_chunk)
    real(dp), intent(out) :: ox(3, panel_chunk), oy(3, panel_chunk), oz(3, panel_chunk), scale(panel_chunk)
    real(dp), intent(out) :: shift(3, panel_chunk)
    logical, intent(out) :: active(panel_chunk)
    integer(i32) :: k, idx

    do k = 1_i32, chunk_size
      idx = index(k)
      active(k) = abs(charge(idx)) > tiny(1.0_dp)
      shift(:, k) = 0.0_dp
      if (axes(1) > 0_i32) shift(axes(1), k) = shift(axes(1), k) + shift1(k)
      if (axes(2) > 0_i32) shift(axes(2), k) = shift(axes(2), k) + shift2(k)
      vx(:, k) = geometry(idx)%vertex(1, :) + shift(1, k)
      vy(:, k) = geometry(idx)%vertex(2, :) + shift(2, k)
      vz(:, k) = geometry(idx)%vertex(3, :) + shift(3, k)
      nx(k) = geometry(idx)%normal(1)
      ny(k) = geometry(idx)%normal(2)
      nz(k) = geometry(idx)%normal(3)
      edge(:, k) = geometry(idx)%edge_length
      ox(:, k) = geometry(idx)%edge_outward(1, :)
      oy(:, k) = geometry(idx)%edge_outward(2, :)
      oz(:, k) = geometry(idx)%edge_outward(3, :)
      scale(k) = k_coulomb*charge(idx)/geometry(idx)%area
    end do
  end subroutine gather_translated_panels

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

  pure real(dp) function product_norms(vectors) result(product_value)
    real(dp), intent(in) :: vectors(3, 3)

    product_value = norm2(vectors(:, 1))*norm2(vectors(:, 2))*norm2(vectors(:, 3))
  end function product_norms

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

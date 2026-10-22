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
  !! `shift1(k)` and `shift2(k)` along `axes` and carrying `charge(index(k))`. Each panel uses the
  !! arithmetic of `panel_potential_field_translated`, so the sum matches a loop over that routine bit
  !! for bit; the loop body is inlined here because this sum dominates the particle step.
  subroutine panel_field_sum_translated(geometry, index, axes, shift1, shift2, charge, count, target, ex, ey, ez)
    type(panel_geometry_type), intent(in) :: geometry(:)
    integer(i32), intent(in) :: index(:), axes(2)
    real(dp), intent(in) :: shift1(:), shift2(:), charge(:)
    integer(i32), intent(in) :: count
    real(dp), intent(in) :: target(3)
    real(dp), intent(inout) :: ex, ey, ez
    real(dp) :: shift(3), vertex(3, 3), relative_vertex(3, 3), radius(3), field(3), potential
    real(dp) :: height, field_integral(3), omega, line_integral, ratio, numerator, denominator
    integer(i32) :: k, idx
    integer :: edge, next_edge

    do k = 1_i32, count
      idx = index(k)
      if (abs(charge(idx)) <= tiny(1.0_dp)) cycle
      shift = 0.0_dp
      if (axes(1) > 0_i32) shift(axes(1)) = shift(axes(1)) + shift1(k)
      if (axes(2) > 0_i32) shift(axes(2)) = shift(axes(2)) + shift2(k)
      vertex = geometry(idx)%vertex + spread(shift, 2, 3)
      height = dot_product(target - vertex(:, 1), geometry(idx)%normal)
      if (abs(height) <= 64.0_dp*epsilon(1.0_dp)*maxval(geometry(idx)%edge_length)) then
        call panel_potential_field_on_translated_surface( &
          geometry(idx), vertex, shift, charge(idx), target, potential, field &
          )
      else
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
        field_integral = geometry(idx)%normal*omega
        do edge = 1, 3
          next_edge = merge(edge + 1, 1, edge < 3)
          ratio = geometry(idx)%edge_length(edge)/(radius(edge) + radius(next_edge))
          if (ratio >= 1.0_dp) error stop 'panel target lies on an edge or vertex.'
          line_integral = log((1.0_dp + ratio)/(1.0_dp - ratio))
          field_integral = field_integral + geometry(idx)%edge_outward(:, edge)*line_integral
        end do
        field = k_coulomb*charge(idx)/geometry(idx)%area*field_integral
      end if
      ex = ex + field(1)/k_coulomb
      ey = ey + field(2)/k_coulomb
      ez = ez + field(3)/k_coulomb
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
    real(dp) :: shift(3), vertex(3, 3), relative_vertex(3, 3), radius(3), field(3), potential
    real(dp) :: projection(3), height, potential_integral, omega, line_integral, distance, ratio
    real(dp) :: numerator, denominator
    integer(i32) :: k, idx
    integer :: edge, next_edge

    do k = 1_i32, count
      idx = index(k)
      if (abs(charge(idx)) <= tiny(1.0_dp)) cycle
      shift = 0.0_dp
      if (axes(1) > 0_i32) shift(axes(1)) = shift(axes(1)) + shift1(k)
      if (axes(2) > 0_i32) shift(axes(2)) = shift(axes(2)) + shift2(k)
      vertex = geometry(idx)%vertex + spread(shift, 2, 3)
      height = dot_product(target - vertex(:, 1), geometry(idx)%normal)
      projection = target - height*geometry(idx)%normal
      if (abs(height) <= 64.0_dp*epsilon(1.0_dp)*maxval(geometry(idx)%edge_length)) then
        call panel_potential_field_on_translated_surface( &
          geometry(idx), vertex, shift, charge(idx), target, potential, field &
          )
      else
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
        do edge = 1, 3
          next_edge = merge(edge + 1, 1, edge < 3)
          ratio = geometry(idx)%edge_length(edge)/(radius(edge) + radius(next_edge))
          if (ratio >= 1.0_dp) error stop 'panel target lies on an edge or vertex.'
          line_integral = log((1.0_dp + ratio)/(1.0_dp - ratio))
          distance = dot_product(vertex(:, edge) - projection, geometry(idx)%edge_outward(:, edge))
          potential_integral = potential_integral + distance*line_integral
        end do
        potential = k_coulomb*charge(idx)/geometry(idx)%area*potential_integral
      end if
      phi = phi + potential/k_coulomb
    end do
  end subroutine panel_potential_sum_translated

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

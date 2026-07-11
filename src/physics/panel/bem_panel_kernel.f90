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
    scale = max(1.0_dp, maxval(geometry%edge_length))
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

  pure logical function point_in_triangle(geometry, point) result(inside)
    type(panel_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: point(3)
    real(dp) :: edge_vector(3), point_vector(3), signed_edge, tolerance
    integer :: edge, next_edge

    tolerance = 64.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(geometry%edge_length))
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

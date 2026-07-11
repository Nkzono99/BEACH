!> On-surface P0 triangle potential and principal-value field integrals.
module bem_panel_self_terms
  use bem_kinds, only: dp
  use bem_panel_geometry, only: panel_geometry_type
  implicit none
  private

  public :: panel_on_surface_integrals

contains

  subroutine panel_on_surface_integrals(geometry, target, potential_integral, field_pv_integral)
    type(panel_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: target(3)
    real(dp), intent(out) :: potential_integral, field_pv_integral(3)
    real(dp) :: endpoint(3), radius_a, radius_b, ratio, line_integral, distance
    integer :: edge, next_edge

    potential_integral = 0.0_dp
    field_pv_integral = 0.0_dp
    do edge = 1, 3
      next_edge = merge(edge + 1, 1, edge < 3)
      endpoint = geometry%vertex(:, edge)
      radius_a = norm2(target - endpoint)
      radius_b = norm2(target - geometry%vertex(:, next_edge))
      ratio = geometry%edge_length(edge)/(radius_a + radius_b)
      if (ratio >= 1.0_dp) error stop 'panel on-surface target lies on an edge or vertex.'
      line_integral = log((1.0_dp + ratio)/(1.0_dp - ratio))
      distance = dot_product(endpoint - target, geometry%edge_outward(:, edge))
      potential_integral = potential_integral + distance*line_integral
      field_pv_integral = field_pv_integral + geometry%edge_outward(:, edge)*line_integral
    end do
  end subroutine panel_on_surface_integrals

  pure real(dp) function norm2(vector) result(norm)
    real(dp), intent(in) :: vector(3)

    norm = sqrt(sum(vector*vector))
  end function norm2

end module bem_panel_self_terms

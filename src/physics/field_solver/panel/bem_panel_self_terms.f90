!> On-surface P0 triangle potential and principal-value field integrals.
module bem_panel_self_terms
  use bem_kinds, only: dp
  use bem_panel_geometry, only: panel_geometry_type
  implicit none
  private

  public :: panel_on_surface_integrals
  public :: panel_edge_line_integral, panel_edge_ratio_limit

  real(dp), parameter :: panel_edge_ratio_limit = 0.99_dp

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
      line_integral = panel_edge_line_integral(geometry, target, edge, ratio)
      distance = dot_product(endpoint - target, geometry%edge_outward(:, edge))
      potential_integral = potential_integral + distance*line_integral
      field_pv_integral = field_pv_integral + geometry%edge_outward(:, edge)*line_integral
    end do
  end subroutine panel_on_surface_integrals

  !> Integral of 1/r along one edge. Near the edge, avoid the cancellation in ra+rb-L:
  !! the perpendicular distance remains resolved even when L/(ra+rb) rounds to one.
  real(dp) function panel_edge_line_integral(geometry, target, edge, ratio) result(integral)
    type(panel_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: target(3), ratio
    integer, intent(in) :: edge
    real(dp) :: offset(3), tangent(3), along, perpendicular
    integer :: next_edge

    if (ratio < panel_edge_ratio_limit) then
      integral = log((1.0_dp + ratio)/(1.0_dp - ratio))
      return
    end if
    next_edge = merge(edge + 1, 1, edge < 3)
    offset = target - geometry%vertex(:, edge)
    tangent = (geometry%vertex(:, next_edge) - geometry%vertex(:, edge))/geometry%edge_length(edge)
    along = dot_product(offset, tangent)
    perpendicular = hypot(dot_product(offset, geometry%edge_outward(:, edge)), dot_product(offset, geometry%normal))
    if (perpendicular == 0.0_dp) then
      if (along >= 0.0_dp .and. along <= geometry%edge_length(edge)) then
        error stop 'panel target lies on an edge or vertex.'
      end if
      ! A coplanar target on the edge's extension is regular.
      integral = abs(log(abs(along)) - log(abs(along - geometry%edge_length(edge))))
    else
      integral = asinh(along/perpendicular) + asinh((geometry%edge_length(edge) - along)/perpendicular)
    end if
  end function panel_edge_line_integral

  pure real(dp) function norm2(vector) result(norm)
    real(dp), intent(in) :: vector(3)

    norm = sqrt(sum(vector*vector))
  end function norm2

end module bem_panel_self_terms

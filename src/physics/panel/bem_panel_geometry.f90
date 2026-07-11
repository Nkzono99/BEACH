!> Ordered triangle geometry and exact Cartesian surface moments.
module bem_panel_geometry
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  implicit none
  private

  integer(i32), parameter, public :: panel_geometry_ok = 0_i32
  integer(i32), parameter, public :: panel_geometry_nonfinite = 1_i32
  integer(i32), parameter, public :: panel_geometry_degenerate = 2_i32

  type, public :: panel_geometry_type
    real(dp) :: vertex(3, 3) = 0.0_dp
    real(dp) :: centroid(3) = 0.0_dp
    real(dp) :: normal(3) = 0.0_dp
    real(dp) :: area = 0.0_dp
    real(dp) :: edge_length(3) = 0.0_dp
    real(dp) :: edge_outward(3, 3) = 0.0_dp
    real(dp) :: moment0 = 0.0_dp
    real(dp) :: moment1(3) = 0.0_dp
    real(dp) :: moment2(3, 3) = 0.0_dp
  end type panel_geometry_type

  public :: init_panel_geometry

contains

  subroutine init_panel_geometry(v0, v1, v2, geometry, status)
    real(dp), intent(in) :: v0(3), v1(3), v2(3)
    type(panel_geometry_type), intent(out) :: geometry
    integer(i32), intent(out) :: status
    real(dp) :: edge1(3), edge2(3), cross12(3), norm_cross, scale
    real(dp) :: vertex_sum(3)
    integer :: edge, next_edge

    geometry = panel_geometry_type()
    status = panel_geometry_nonfinite
    if (.not. all(ieee_is_finite(v0)) .or. .not. all(ieee_is_finite(v1)) .or. &
        .not. all(ieee_is_finite(v2))) return

    geometry%vertex(:, 1) = v0
    geometry%vertex(:, 2) = v1
    geometry%vertex(:, 3) = v2
    edge1 = v1 - v0
    edge2 = v2 - v0
    cross12 = cross_product(edge1, edge2)
    norm_cross = sqrt(sum(cross12*cross12))
    scale = max(1.0_dp, maxval(abs(geometry%vertex)))
    status = panel_geometry_degenerate
    if (norm_cross <= 64.0_dp*epsilon(1.0_dp)*scale*scale) return

    geometry%area = 0.5_dp*norm_cross
    geometry%normal = cross12/norm_cross
    geometry%centroid = (v0 + v1 + v2)/3.0_dp
    do edge = 1, 3
      next_edge = merge(edge + 1, 1, edge < 3)
      edge1 = geometry%vertex(:, next_edge) - geometry%vertex(:, edge)
      geometry%edge_length(edge) = sqrt(sum(edge1*edge1))
      geometry%edge_outward(:, edge) = cross_product(edge1, geometry%normal)/geometry%edge_length(edge)
    end do

    geometry%moment0 = geometry%area
    geometry%moment1 = geometry%area*geometry%centroid
    vertex_sum = v0 + v1 + v2
    geometry%moment2 = geometry%area/12.0_dp*( &
                       outer_product(v0, v0) + outer_product(v1, v1) + outer_product(v2, v2) + &
                       outer_product(vertex_sum, vertex_sum) &
                       )
    status = panel_geometry_ok
  end subroutine init_panel_geometry

  pure function cross_product(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)

    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  end function cross_product

  pure function outer_product(a, b) result(product)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: product(3, 3)
    integer :: i, j

    do j = 1, 3
      do i = 1, 3
        product(i, j) = a(i)*b(j)
      end do
    end do
  end function outer_product

end module bem_panel_geometry

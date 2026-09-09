!> Independent Coulomb potential/field correctness oracles for triangle panels.
module bem_panel_quadrature
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb
  use bem_panel_geometry, only: panel_geometry_type
  use bem_triangle_quadrature, only: panel_quadrature_plan_type, build_panel_quadrature, &
                                     fill_panel_quadrature, build_panel_duffy_quadrature, gauss_legendre_unit
  implicit none
  private

  public :: panel_quadrature_plan_type

  public :: build_panel_quadrature
  public :: fill_panel_quadrature
  public :: build_panel_duffy_quadrature
  public :: panel_oracle_potential_field
  public :: panel_singular_potential_oracle

contains

  subroutine panel_oracle_potential_field(geometry, charge, target, order, potential, field)
    type(panel_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: charge, target(3)
    integer(i32), intent(in) :: order
    real(dp), intent(out) :: potential, field(3)
    real(dp), allocatable :: node(:), weight(:)
    real(dp) :: source(3), displacement(3), radius2, jacobian, potential_integral, field_integral(3)
    real(dp) :: u, v, edge1(3), edge2(3), direction(3)
    integer :: iu, iv

    call gauss_legendre_unit(order, node, weight)
    edge1 = geometry%vertex(:, 2) - geometry%vertex(:, 1)
    edge2 = geometry%vertex(:, 3) - geometry%vertex(:, 1)
    potential_integral = 0.0_dp
    field_integral = 0.0_dp
    do iu = 1, order
      u = node(iu)
      do iv = 1, order
        v = node(iv)
        direction = (1.0_dp - v)*edge1 + v*edge2
        source = geometry%vertex(:, 1) + u*direction
        displacement = target - source
        radius2 = sum(displacement*displacement)
        jacobian = 2.0_dp*geometry%area*u*weight(iu)*weight(iv)
        potential_integral = potential_integral + jacobian/sqrt(radius2)
        field_integral = field_integral + jacobian*displacement/(radius2*sqrt(radius2))
      end do
    end do
    potential = k_coulomb*charge/geometry%area*potential_integral
    field = k_coulomb*charge/geometry%area*field_integral
  end subroutine panel_oracle_potential_field

  subroutine panel_singular_potential_oracle(geometry, charge, target, order, potential)
    type(panel_geometry_type), intent(in) :: geometry
    real(dp), intent(in) :: charge, target(3)
    integer(i32), intent(in) :: order
    real(dp), intent(out) :: potential
    real(dp), allocatable :: node(:), weight(:)
    real(dp) :: edge1(3), edge2(3), direction(3), source(3), radius, jacobian, integral
    real(dp) :: u, v, sub_area2
    integer :: edge, next_edge, iu, iv

    call gauss_legendre_unit(order, node, weight)
    integral = 0.0_dp
    do edge = 1, 3
      next_edge = merge(edge + 1, 1, edge < 3)
      edge1 = geometry%vertex(:, edge) - target
      edge2 = geometry%vertex(:, next_edge) - target
      sub_area2 = sqrt(sum(cross_product(edge1, edge2)**2))
      do iu = 1, order
        u = node(iu)
        do iv = 1, order
          v = node(iv)
          direction = (1.0_dp - v)*edge1 + v*edge2
          source = target + u*direction
          radius = sqrt(sum((target - source)**2))
          jacobian = sub_area2*u*weight(iu)*weight(iv)
          integral = integral + jacobian/radius
        end do
      end do
    end do
    potential = k_coulomb*charge/geometry%area*integral
  end subroutine panel_singular_potential_oracle

  pure function cross_product(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)

    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  end function cross_product

end module bem_panel_quadrature

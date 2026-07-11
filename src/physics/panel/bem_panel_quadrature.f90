!> Independent triangle cubature and Gauss-Duffy correctness oracles.
module bem_panel_quadrature
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_coulomb, pi
  use bem_panel_geometry, only: panel_geometry_type
  implicit none
  private

  type, public :: panel_quadrature_plan_type
    integer(i32) :: npoint = 0_i32
    real(dp), allocatable :: position(:, :)
    real(dp), allocatable :: weight(:)
  end type panel_quadrature_plan_type

  public :: build_panel_quadrature
  public :: panel_oracle_potential_field
  public :: panel_singular_potential_oracle

contains

  subroutine build_panel_quadrature(geometry, plan)
    type(panel_geometry_type), intent(in) :: geometry
    type(panel_quadrature_plan_type), intent(out) :: plan
    real(dp), parameter :: barycentric(3, 7) = reshape([ &
                                                       1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp, &
                                                       0.059715871789770_dp, 0.470142064105115_dp, 0.470142064105115_dp, &
                                                       0.470142064105115_dp, 0.059715871789770_dp, 0.470142064105115_dp, &
                                                       0.470142064105115_dp, 0.470142064105115_dp, 0.059715871789770_dp, &
                                                       0.797426985353087_dp, 0.101286507323456_dp, 0.101286507323456_dp, &
                                                       0.101286507323456_dp, 0.797426985353087_dp, 0.101286507323456_dp, &
                                                       0.101286507323456_dp, 0.101286507323456_dp, 0.797426985353087_dp &
                                                       ], [3, 7])
    real(dp), parameter :: normalized_weight(7) = [ &
                           0.225000000000000_dp, &
                           0.132394152788506_dp, 0.132394152788506_dp, 0.132394152788506_dp, &
                           0.125939180544827_dp, 0.125939180544827_dp, 0.125939180544827_dp &
                           ]
    integer :: point

    plan%npoint = 7_i32
    allocate (plan%position(3, plan%npoint), plan%weight(plan%npoint))
    do point = 1, plan%npoint
      plan%position(:, point) = matmul(geometry%vertex, barycentric(:, point))
    end do
    plan%weight = geometry%area*normalized_weight
  end subroutine build_panel_quadrature

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

  subroutine gauss_legendre_unit(order, node, weight)
    integer(i32), intent(in) :: order
    real(dp), allocatable, intent(out) :: node(:), weight(:)
    integer :: i, j, midpoint
    real(dp) :: z, z_previous, polynomial, derivative, p0, p1, p2

    if (order < 2_i32) error stop 'panel oracle quadrature order must be >= 2.'
    allocate (node(order), weight(order))
    midpoint = (order + 1)/2
    do i = 1, midpoint
      z = cos(pi*(real(i, dp) - 0.25_dp)/(real(order, dp) + 0.5_dp))
      do
        p0 = 1.0_dp
        p1 = z
        do j = 2, order
          p2 = ((2.0_dp*real(j, dp) - 1.0_dp)*z*p1 - (real(j, dp) - 1.0_dp)*p0)/real(j, dp)
          p0 = p1
          p1 = p2
        end do
        polynomial = merge(1.0_dp, p1, order == 0_i32)
        derivative = real(order, dp)*(z*p1 - p0)/(z*z - 1.0_dp)
        z_previous = z
        z = z_previous - polynomial/derivative
        if (abs(z - z_previous) <= 8.0_dp*epsilon(1.0_dp)) exit
      end do
      node(i) = 0.5_dp*(1.0_dp - z)
      node(order + 1 - i) = 0.5_dp*(1.0_dp + z)
      weight(i) = 1.0_dp/((1.0_dp - z*z)*derivative*derivative)
      weight(order + 1 - i) = weight(i)
    end do
  end subroutine gauss_legendre_unit

  pure function cross_product(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)

    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  end function cross_product

end module bem_panel_quadrature

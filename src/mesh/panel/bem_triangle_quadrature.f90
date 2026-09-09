!> Triangle cubature rules and Gauss-Duffy quadrature plans.
module bem_triangle_quadrature
  use bem_kinds, only: dp, i32
  use bem_constants, only: pi
  use bem_panel_geometry, only: panel_geometry_type
  implicit none
  private

  type, public :: panel_quadrature_plan_type
    integer(i32) :: npoint = 0_i32
    real(dp), allocatable :: position(:, :)
    real(dp), allocatable :: weight(:)
  end type panel_quadrature_plan_type

  public :: build_panel_quadrature
  public :: fill_panel_quadrature
  public :: build_panel_duffy_quadrature
  public :: gauss_legendre_unit

contains

  subroutine build_panel_duffy_quadrature(geometry, order, plan)
    type(panel_geometry_type), intent(in) :: geometry
    integer(i32), intent(in) :: order
    type(panel_quadrature_plan_type), intent(out) :: plan
    real(dp), allocatable :: node(:), weight(:)
    real(dp) :: edge1(3), edge2(3), direction(3), u, v
    integer :: iu, iv, point

    call gauss_legendre_unit(order, node, weight)
    plan%npoint = order*order
    allocate (plan%position(3, plan%npoint), plan%weight(plan%npoint))
    edge1 = geometry%vertex(:, 2) - geometry%vertex(:, 1)
    edge2 = geometry%vertex(:, 3) - geometry%vertex(:, 1)
    point = 0
    do iu = 1, order
      u = node(iu)
      do iv = 1, order
        v = node(iv)
        point = point + 1
        direction = (1.0_dp - v)*edge1 + v*edge2
        plan%position(:, point) = geometry%vertex(:, 1) + u*direction
        plan%weight(point) = 2.0_dp*geometry%area*u*weight(iu)*weight(iv)
      end do
    end do
  end subroutine build_panel_duffy_quadrature

  subroutine build_panel_quadrature(geometry, plan)
    type(panel_geometry_type), intent(in) :: geometry
    type(panel_quadrature_plan_type), intent(out) :: plan

    plan%npoint = 7_i32
    allocate (plan%position(3, plan%npoint), plan%weight(plan%npoint))
    call fill_panel_quadrature(geometry, plan%position, plan%weight)
  end subroutine build_panel_quadrature

  !> Fill an existing seven-point cubature buffer without allocating a plan.
  subroutine fill_panel_quadrature(geometry, position, weight)
    type(panel_geometry_type), intent(in) :: geometry
    real(dp), intent(out) :: position(3, 7), weight(7)
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

    do point = 1, 7
      position(:, point) = matmul(geometry%vertex, barycentric(:, point))
    end do
    weight = geometry%area*normalized_weight
  end subroutine fill_panel_quadrature

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

end module bem_triangle_quadrature

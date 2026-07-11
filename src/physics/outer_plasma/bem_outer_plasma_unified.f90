module bem_outer_plasma_unified
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok, outer_plasma_invalid, &
                                    outer_plasma_not_applicable, outer_plasma_numerical_failure
  implicit none
  private

  type, public :: unified_outer_linear_options_type
    real(dp) :: kappa = 0.0_dp
    real(dp) :: tail_length = 0.0_dp
    real(dp) :: bottom_field = 0.0_dp
    real(dp) :: thermal_voltage = 0.0_dp
    real(dp) :: max_linearity_ratio = 0.1_dp
  end type unified_outer_linear_options_type

  public :: solve_unified_outer_linear

contains

  subroutine solve_unified_outer_linear( &
    z, surface_field, accessible_fraction, options, state, status, message &
    )
    real(dp), intent(in) :: z(:), surface_field(:), accessible_fraction(:)
    type(unified_outer_linear_options_type), intent(in) :: options
    type(outer_plasma_state_type), intent(out) :: state
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    real(dp), allocatable :: lower(:), diagonal(:), upper(:), rhs(:), phi(:), surface_rho(:)
    real(dp) :: left_h, right_h
    integer(i32) :: point, n
    logical :: solved

    state = outer_plasma_state_type()
    state%model = 'unified_linear_response'
    status = outer_plasma_invalid
    message = ''
    n = size(z)
    if (n < 3_i32 .or. size(surface_field) /= n .or. size(accessible_fraction) /= n .or. &
        .not. all(ieee_is_finite([z, surface_field, accessible_fraction])) .or. &
        any(z(2:) <= z(:n - 1_i32)) .or. any(accessible_fraction < 0.0_dp) .or. &
        any(accessible_fraction > 1.0_dp) .or. options%kappa <= 0.0_dp .or. &
        options%tail_length <= 0.0_dp .or. options%thermal_voltage <= 0.0_dp .or. &
        options%max_linearity_ratio <= 0.0_dp .or. &
        .not. all(ieee_is_finite([ &
                                 options%kappa, options%tail_length, options%bottom_field, &
                                 options%thermal_voltage, options%max_linearity_ratio &
                                 ]))) then
      message = 'invalid unified linear outer-plasma inputs'
      state%applicability_status = status
      return
    end if

    allocate (lower(n), diagonal(n), upper(n), rhs(n), phi(n), surface_rho(n))
    lower = 0.0_dp
    diagonal = 0.0_dp
    upper = 0.0_dp
    rhs = 0.0_dp
    surface_rho = 0.0_dp
    do point = 2_i32, n - 1_i32
      surface_rho(point) = eps0*(surface_field(point + 1_i32) - surface_field(point - 1_i32))/ &
                           (z(point + 1_i32) - z(point - 1_i32))
    end do

    right_h = z(2) - z(1)
    diagonal(1) = -2.0_dp/right_h**2 - accessible_fraction(1)*options%kappa**2
    upper(1) = 2.0_dp/right_h**2
    rhs(1) = -surface_rho(1)/eps0 - 2.0_dp*options%bottom_field/right_h
    do point = 2_i32, n - 1_i32
      left_h = z(point) - z(point - 1_i32)
      right_h = z(point + 1_i32) - z(point)
      lower(point) = 2.0_dp/(left_h*(left_h + right_h))
      upper(point) = 2.0_dp/(right_h*(left_h + right_h))
      diagonal(point) = -lower(point) - upper(point) - accessible_fraction(point)*options%kappa**2
      rhs(point) = -surface_rho(point)/eps0
    end do
    left_h = z(n) - z(n - 1_i32)
    lower(n) = 2.0_dp/left_h**2
    diagonal(n) = -2.0_dp/left_h**2 - 2.0_dp/(left_h*options%tail_length) - &
                  accessible_fraction(n)*options%kappa**2
    rhs(n) = -surface_rho(n)/eps0
    call solve_tridiagonal(lower, diagonal, upper, rhs, phi, solved)
    if (.not. solved) then
      status = outer_plasma_numerical_failure
      state%applicability_status = status
      message = 'unified linear Poisson system is singular or non-finite'
      return
    end if

    state%linearity_ratio = maxval(abs(phi))/options%thermal_voltage
    state%max_linearity_ratio = options%max_linearity_ratio
    if (state%linearity_ratio > options%max_linearity_ratio) then
      status = outer_plasma_not_applicable
      state%applicability_status = status
      message = 'unified local/nonzero response exceeds the linearity contract'
      return
    end if
    state%profile_n = n
    state%interface_z = z(1)
    state%interface_potential = phi(1)
    state%infinity_potential = 0.0_dp
    state%debye_length = 1.0_dp/options%kappa
    state%interface_field = options%bottom_field
    state%z = z
    state%potential = phi
    allocate (state%field(n), state%charge_density(n))
    state%field(1) = options%bottom_field
    do point = 2_i32, n - 1_i32
      state%field(point) = -(phi(point + 1_i32) - phi(point - 1_i32))/ &
                           (z(point + 1_i32) - z(point - 1_i32))
    end do
    state%field(n) = -(phi(n) - phi(n - 1_i32))/(z(n) - z(n - 1_i32))
    state%charge_density = -eps0*accessible_fraction*options%kappa**2*phi
    state%integrated_charge_per_area = -eps0*( &
                                       options%bottom_field + surface_field(n) - surface_field(1) &
                                       )
    state%ready = .true.
    status = outer_plasma_ok
    state%applicability_status = status
  end subroutine solve_unified_outer_linear

  subroutine solve_tridiagonal(lower, diagonal, upper, rhs, solution, success)
    real(dp), intent(in) :: lower(:), diagonal(:), upper(:), rhs(:)
    real(dp), intent(out) :: solution(:)
    logical, intent(out) :: success
    real(dp) :: modified_upper(size(rhs)), modified_rhs(size(rhs)), denominator
    integer(i32) :: point, n

    n = size(rhs)
    success = .false.
    if (abs(diagonal(1)) <= tiny(1.0_dp)) return
    modified_upper(1) = upper(1)/diagonal(1)
    modified_rhs(1) = rhs(1)/diagonal(1)
    do point = 2_i32, n
      denominator = diagonal(point) - lower(point)*modified_upper(point - 1_i32)
      if (.not. ieee_is_finite(denominator) .or. abs(denominator) <= tiny(1.0_dp)) return
      modified_upper(point) = upper(point)/denominator
      modified_rhs(point) = (rhs(point) - lower(point)*modified_rhs(point - 1_i32))/denominator
    end do
    solution(n) = modified_rhs(n)
    do point = n - 1_i32, 1_i32, -1_i32
      solution(point) = modified_rhs(point) - modified_upper(point)*solution(point + 1_i32)
    end do
    success = all(ieee_is_finite(solution))
  end subroutine solve_tridiagonal

  pure real(dp) function trapz(x, y) result(integral)
    real(dp), intent(in) :: x(:), y(:)
    integer(i32) :: point

    integral = 0.0_dp
    do point = 1_i32, size(x) - 1_i32
      integral = integral + 0.5_dp*(y(point) + y(point + 1_i32))*(x(point + 1_i32) - x(point))
    end do
  end function trapz

end module bem_outer_plasma_unified

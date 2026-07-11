module bem_periodic_zero_mode_eval
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_periodic_zero_mode_plan, only: periodic_zero_mode_plan_type, periodic_zero_mode_state_type
  implicit none
  private

  integer(i32), parameter, public :: zero_mode_trace_minus = -1_i32
  integer(i32), parameter, public :: zero_mode_trace_principal_value = 0_i32
  integer(i32), parameter, public :: zero_mode_trace_plus = 1_i32

  public :: eval_periodic_zero_mode

contains

  subroutine eval_periodic_zero_mode(plan, state, z, trace, potential, field)
    type(periodic_zero_mode_plan_type), intent(in) :: plan
    type(periodic_zero_mode_state_type), intent(in) :: state
    real(dp), intent(in) :: z
    integer(i32), intent(in) :: trace
    real(dp), intent(out) :: potential, field
    real(dp) :: cumulative_charge, sheet_correction, primitive, gauge_primitive
    integer(i32) :: interval, breakpoint

    if (trace < zero_mode_trace_minus .or. trace > zero_mode_trace_plus) then
      error stop 'invalid periodic zero-mode trace.'
    end if
    call locate_interval(plan%break_z, z, interval, breakpoint)
    cumulative_charge = evaluate_polynomial(state%cumulative_charge_coeff(:, interval), z)
    if (breakpoint > 0_i32) then
      sheet_correction = 0.5_dp*real(1_i32 - trace, dp)*state%sheet_charge(breakpoint)
      cumulative_charge = cumulative_charge - sheet_correction
    end if
    field = state%e_bottom + cumulative_charge/(eps0*plan%area_xy)
    primitive = charge_primitive(plan, state, z)
    gauge_primitive = charge_primitive(plan, state, state%z_gauge)
    potential = state%phi_gauge - state%e_bottom*(z - state%z_gauge) - &
                (primitive - gauge_primitive)/(eps0*plan%area_xy)
  end subroutine eval_periodic_zero_mode

  pure real(dp) function charge_primitive(plan, state, z) result(primitive)
    type(periodic_zero_mode_plan_type), intent(in) :: plan
    type(periodic_zero_mode_state_type), intent(in) :: state
    real(dp), intent(in) :: z
    integer(i32) :: interval, breakpoint, lower_break

    call locate_interval(plan%break_z, z, interval, breakpoint)
    if (interval == 1_i32) then
      primitive = integrate_polynomial(state%cumulative_charge_coeff(:, 1), plan%break_z(1), z)
      return
    end if
    lower_break = interval - 1_i32
    primitive = state%primitive_at_break(lower_break) + &
                integrate_polynomial( &
                state%cumulative_charge_coeff(:, interval), plan%break_z(lower_break), z &
                )
  end function charge_primitive

  pure subroutine locate_interval(break_z, z, interval, breakpoint)
    real(dp), intent(in) :: break_z(:), z
    integer(i32), intent(out) :: interval, breakpoint
    integer(i32) :: low, high, middle
    real(dp) :: tolerance

    tolerance = 128.0_dp*epsilon(1.0_dp)*max(1.0_dp, max(abs(z), maxval(abs(break_z))))
    low = 1_i32
    high = size(break_z)
    breakpoint = 0_i32
    do while (low <= high)
      middle = (low + high)/2_i32
      if (abs(z - break_z(middle)) <= tolerance) then
        breakpoint = middle
        interval = middle + 1_i32
        return
      else if (z < break_z(middle)) then
        high = middle - 1_i32
      else
        low = middle + 1_i32
      end if
    end do
    interval = low
  end subroutine locate_interval

  pure real(dp) function evaluate_polynomial(coefficient, z) result(value)
    real(dp), intent(in) :: coefficient(3), z

    value = coefficient(1) + z*(coefficient(2) + z*coefficient(3))
  end function evaluate_polynomial

  pure real(dp) function integrate_polynomial(coefficient, lower, upper) result(integral)
    real(dp), intent(in) :: coefficient(3), lower, upper

    integral = coefficient(1)*(upper - lower) + &
               0.5_dp*coefficient(2)*(upper*upper - lower*lower) + &
               coefficient(3)/3.0_dp*(upper**3 - lower**3)
  end function integrate_polynomial

end module bem_periodic_zero_mode_eval

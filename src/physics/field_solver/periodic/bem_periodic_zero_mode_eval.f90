module bem_periodic_zero_mode_eval
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_periodic_zero_mode_plan, only: periodic_zero_mode_plan_type, periodic_zero_mode_state_type, &
                                         zero_mode_interval_origin, integrate_local_polynomial, &
                                         zero_mode_heights_coincide
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
    integer(i32) :: interval, breakpoint, gauge_interval

    if (trace < zero_mode_trace_minus .or. trace > zero_mode_trace_plus) then
      error stop 'invalid periodic zero-mode trace.'
    end if
    call locate_interval(plan%break_z, plan%height_tolerance, z, interval, breakpoint)
    cumulative_charge = evaluate_polynomial( &
                        state%cumulative_charge_coeff(:, interval), z - zero_mode_interval_origin(plan%break_z, interval) &
                        )
    if (breakpoint > 0_i32) then
      sheet_correction = 0.5_dp*real(1_i32 - trace, dp)*state%sheet_charge(breakpoint)
      cumulative_charge = cumulative_charge - sheet_correction
    end if
    field = state%e_bottom + cumulative_charge/(eps0*plan%area_xy)
    primitive = charge_primitive(plan, state, z, interval)
    call locate_interval(plan%break_z, plan%height_tolerance, state%z_gauge, gauge_interval, breakpoint)
    gauge_primitive = charge_primitive(plan, state, state%z_gauge, gauge_interval)
    potential = state%phi_gauge - state%e_bottom*(z - state%z_gauge) - &
                (primitive - gauge_primitive)/(eps0*plan%area_xy)
  end subroutine eval_periodic_zero_mode

  pure real(dp) function charge_primitive(plan, state, z, interval) result(primitive)
    type(periodic_zero_mode_plan_type), intent(in) :: plan
    type(periodic_zero_mode_state_type), intent(in) :: state
    real(dp), intent(in) :: z
    integer(i32), intent(in) :: interval

    ! 区間 1 は最下 break から、それ以外は区間の左端 break（局所原点）から積分する。
    primitive = integrate_local_polynomial( &
                state%cumulative_charge_coeff(:, interval), z - zero_mode_interval_origin(plan%break_z, interval) &
                )
    if (interval > 1_i32) primitive = primitive + state%primitive_at_break(interval - 1_i32)
  end function charge_primitive

  !> A target on a breakpoint (same coincidence rule as the plan's breakpoint merge) takes the
  !! requested sheet trace; any other target lies strictly inside one interval.
  pure subroutine locate_interval(break_z, tolerance, z, interval, breakpoint)
    real(dp), intent(in) :: break_z(:), tolerance, z
    integer(i32), intent(out) :: interval, breakpoint
    integer(i32) :: low, high, middle

    low = 1_i32
    high = size(break_z)
    breakpoint = 0_i32
    do while (low <= high)
      middle = (low + high)/2_i32
      if (zero_mode_heights_coincide(z, break_z(middle), tolerance)) then
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

  pure real(dp) function evaluate_polynomial(coefficient, t) result(value)
    real(dp), intent(in) :: coefficient(3), t

    value = coefficient(1) + t*(coefficient(2) + t*coefficient(3))
  end function evaluate_polynomial

end module bem_periodic_zero_mode_eval

module bem_outer_plasma_linear
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_outer_plasma_types, only: outer_plasma_state_type, outer_plasma_ok, &
                                    outer_plasma_not_applicable, outer_plasma_invalid
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  public :: init_outer_plasma_linear
  public :: eval_outer_plasma_linear
  public :: outer_plasma_integrated_charge_per_area
  public :: outer_plasma_gauss_residual_per_area

contains

  subroutine init_outer_plasma_linear( &
    interface_z, interface_potential, infinity_potential, debye_length, linearity_ratio, &
    max_linearity_ratio, state, status, message &
    )
    real(dp), intent(in) :: interface_z, interface_potential, infinity_potential, debye_length
    real(dp), intent(in) :: linearity_ratio, max_linearity_ratio
    type(outer_plasma_state_type), intent(out) :: state
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    state%model = 'linear_debye'
    state%interface_z = interface_z
    state%interface_potential = interface_potential
    state%infinity_potential = infinity_potential
    state%debye_length = debye_length
    state%linearity_ratio = linearity_ratio
    state%max_linearity_ratio = max_linearity_ratio
    status = outer_plasma_invalid
    message = ''
    if (.not. all(ieee_is_finite([ &
                                 interface_z, interface_potential, infinity_potential, debye_length, &
                                 linearity_ratio, max_linearity_ratio &
                                 ]))) then
      message = 'linear outer-plasma inputs must be finite'
      state%applicability_status = status
      return
    end if
    if (debye_length <= 0.0_dp .or. linearity_ratio < 0.0_dp .or. max_linearity_ratio <= 0.0_dp) then
      message = 'linear outer-plasma scales must be positive'
      state%applicability_status = status
      return
    end if
    if (linearity_ratio > max_linearity_ratio) then
      status = outer_plasma_not_applicable
      message = 'linear outer-plasma perturbation exceeds applicability limit'
      state%applicability_status = status
      return
    end if

    state%interface_field = (interface_potential - infinity_potential)/debye_length
    state%ready = .true.
    status = outer_plasma_ok
    state%applicability_status = status
  end subroutine init_outer_plasma_linear

  subroutine eval_outer_plasma_linear(state, z, potential, field, charge_density)
    type(outer_plasma_state_type), intent(in) :: state
    real(dp), intent(in) :: z
    real(dp), intent(out) :: potential, field, charge_density
    real(dp) :: decay, delta_potential

    if (.not. state%ready) error stop 'linear outer-plasma state is not ready.'
    if (.not. ieee_is_finite(z) .or. z < state%interface_z) then
      error stop 'linear outer-plasma evaluation requires finite z >= interface_z.'
    end if
    delta_potential = state%interface_potential - state%infinity_potential
    decay = exp(-(z - state%interface_z)/state%debye_length)
    potential = state%infinity_potential + delta_potential*decay
    field = state%interface_field*decay
    charge_density = -eps0*delta_potential/(state%debye_length*state%debye_length)*decay
  end subroutine eval_outer_plasma_linear

  real(dp) function outer_plasma_integrated_charge_per_area(state) result(charge_per_area)
    type(outer_plasma_state_type), intent(in) :: state

    if (.not. state%ready) error stop 'linear outer-plasma state is not ready.'
    charge_per_area = -eps0*state%interface_field
  end function outer_plasma_integrated_charge_per_area

  real(dp) function outer_plasma_gauss_residual_per_area(state) result(residual)
    type(outer_plasma_state_type), intent(in) :: state
    real(dp) :: integrated_charge

    integrated_charge = outer_plasma_integrated_charge_per_area(state)
    residual = integrated_charge - eps0*(0.0_dp - state%interface_field)
  end function outer_plasma_gauss_residual_per_area

end module bem_outer_plasma_linear

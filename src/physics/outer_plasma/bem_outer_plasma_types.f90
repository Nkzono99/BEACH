module bem_outer_plasma_types
  use bem_kinds, only: dp, i32
  implicit none
  private

  integer(i32), parameter, public :: outer_plasma_ok = 0_i32
  integer(i32), parameter, public :: outer_plasma_not_applicable = 1_i32
  integer(i32), parameter, public :: outer_plasma_invalid = 2_i32
  integer(i32), parameter, public :: outer_plasma_no_physical_solution = 3_i32
  integer(i32), parameter, public :: outer_plasma_numerical_failure = 4_i32

  type, public :: outer_plasma_state_type
    logical :: ready = .false.
    character(len=32) :: model = 'none'
    integer(i32) :: applicability_status = outer_plasma_invalid
    real(dp) :: interface_z = 0.0_dp
    real(dp) :: interface_potential = 0.0_dp
    real(dp) :: infinity_potential = 0.0_dp
    real(dp) :: debye_length = 0.0_dp
    real(dp) :: interface_field = 0.0_dp
    real(dp) :: linearity_ratio = 0.0_dp
    real(dp) :: max_linearity_ratio = 0.0_dp
    integer(i32) :: profile_n = 0_i32
    integer(i32) :: nonlinear_iterations = 0_i32
    real(dp) :: nonlinear_residual = 0.0_dp
    real(dp) :: integrated_charge_per_area = 0.0_dp
    real(dp), allocatable :: z(:)
    real(dp), allocatable :: potential(:)
    real(dp), allocatable :: field(:)
    real(dp), allocatable :: charge_density(:)
  end type outer_plasma_state_type

end module bem_outer_plasma_types

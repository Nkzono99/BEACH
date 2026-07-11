module bem_interface_types
  use bem_kinds, only: dp, i32
  implicit none
  private

  integer(i32), parameter, public :: interface_outcome_none = 0_i32
  integer(i32), parameter, public :: interface_outcome_returned_local = 1_i32
  integer(i32), parameter, public :: interface_outcome_escaped_to_infinity = 2_i32
  integer(i32), parameter, public :: interface_outcome_queued_outer = 3_i32
  integer(i32), parameter, public :: interface_outcome_invalid_model = 4_i32

  type, public :: interface_crossing_type
    logical :: has_crossing = .false.
    integer(i32) :: face_index = 0_i32
    real(dp) :: fraction = 0.0_dp
    real(dp) :: position(3) = 0.0_dp
    real(dp) :: velocity(3) = 0.0_dp
    real(dp) :: dt_remaining = 0.0_dp
  end type interface_crossing_type

  type, public :: interface_particle_outcome_type
    integer(i32) :: kind = interface_outcome_none
    real(dp) :: position(3) = 0.0_dp
    real(dp) :: velocity(3) = 0.0_dp
    real(dp) :: outer_flight_time = 0.0_dp
    real(dp) :: frozen_field_ratio = 0.0_dp
    real(dp) :: normal_energy_residual = 0.0_dp
    real(dp) :: energy_relative_error = 0.0_dp
    character(len=96) :: message = ''
  end type interface_particle_outcome_type

end module bem_interface_types

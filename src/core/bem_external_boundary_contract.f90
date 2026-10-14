!> 局所 reservoir 流入と種別 kinetic 障壁を実行時契約へ正規化する。
module bem_external_boundary_contract
  use bem_kinds, only: dp, i32
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  integer(i32), parameter, public :: external_boundary_ok = 0_i32
  integer(i32), parameter, public :: external_boundary_invalid = 1_i32

  integer(i32), parameter, public :: external_inflow_none = 0_i32
  integer(i32), parameter, public :: external_inflow_scalar_barrier = 1_i32

  type, public :: external_boundary_contract_type
    integer(i32) :: inflow_map = external_inflow_none
    logical :: barrier_override_low(3) = .false.
    logical :: barrier_override_high(3) = .false.
    real(dp) :: barrier_potential_low_v(3) = 0.0_dp
    real(dp) :: barrier_potential_high_v(3) = 0.0_dp
    logical :: barrier_redistribute_low(3) = .false.
    logical :: barrier_redistribute_high(3) = .false.
  end type external_boundary_contract_type

  public :: resolve_external_boundary_contract

contains

  !> 局所流入補正を解決し、種別 kinetic 障壁を初期化する。
  subroutine resolve_external_boundary_contract( &
    reservoir_potential_model, contract, status, message &
    )
    character(len=*), intent(in) :: reservoir_potential_model
    type(external_boundary_contract_type), intent(out) :: contract
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    character(len=32) :: reservoir

    contract = external_boundary_contract_type()
    status = external_boundary_ok
    message = ''
    reservoir = lower_ascii(trim(reservoir_potential_model))

    select case (trim(reservoir))
    case ('none')
      contract%inflow_map = external_inflow_none
    case ('infinity_barrier')
      contract%inflow_map = external_inflow_scalar_barrier
    case default
      call reject('Unknown sim.reservoir_potential_model.', status, message)
      return
    end select
  end subroutine resolve_external_boundary_contract

  pure subroutine reject(text, status, message)
    character(len=*), intent(in) :: text
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = external_boundary_invalid
    message = text
  end subroutine reject

end module bem_external_boundary_contract

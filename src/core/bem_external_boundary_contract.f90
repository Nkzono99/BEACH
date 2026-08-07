!> 局所 reservoir 流入と通常 open 面の処理を実行時契約へ正規化する。
module bem_external_boundary_contract
  use bem_kinds, only: dp, i32
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  integer(i32), parameter, public :: external_boundary_ok = 0_i32
  integer(i32), parameter, public :: external_boundary_invalid = 1_i32

  integer(i32), parameter, public :: external_inflow_none = 0_i32
  integer(i32), parameter, public :: external_inflow_scalar_barrier = 1_i32

  integer(i32), parameter, public :: external_open_escape = 0_i32
  integer(i32), parameter, public :: external_open_potential_barrier = 1_i32

  type, public :: external_boundary_contract_type
    integer(i32) :: inflow_map = external_inflow_none
    integer(i32) :: ordinary_open_model = external_open_escape
    logical :: barrier_override_low(3) = .false.
    logical :: barrier_override_high(3) = .false.
    real(dp) :: barrier_potential_low_v(3) = 0.0_dp
    real(dp) :: barrier_potential_high_v(3) = 0.0_dp
  end type external_boundary_contract_type

  public :: resolve_external_boundary_contract

contains

  !> 局所流入補正と通常 open 面モデルを一意な組へ解決する。
  subroutine resolve_external_boundary_contract( &
    reservoir_potential_model, open_boundary_model, contract, status, message &
    )
    character(len=*), intent(in) :: reservoir_potential_model
    character(len=*), intent(in) :: open_boundary_model
    type(external_boundary_contract_type), intent(out) :: contract
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    character(len=32) :: reservoir, open_model

    contract = external_boundary_contract_type()
    status = external_boundary_ok
    message = ''
    reservoir = lower_ascii(trim(reservoir_potential_model))
    open_model = lower_ascii(trim(open_boundary_model))

    select case (trim(open_model))
    case ('escape')
      contract%ordinary_open_model = external_open_escape
    case ('potential_barrier')
      contract%ordinary_open_model = external_open_potential_barrier
    case default
      call reject('Unknown sim.open_boundary_model.', status, message)
      return
    end select

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

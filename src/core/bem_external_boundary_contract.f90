!> 外部境界に分散していた設定語彙を、実行時の単一契約へ正規化する。
module bem_external_boundary_contract
  use bem_kinds, only: i32
  use bem_types, only: bc_open
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  integer(i32), parameter, public :: external_boundary_ok = 0_i32
  integer(i32), parameter, public :: external_boundary_invalid = 1_i32

  integer(i32), parameter, public :: external_inflow_none = 0_i32
  integer(i32), parameter, public :: external_inflow_scalar_barrier = 1_i32
  integer(i32), parameter, public :: external_inflow_kinetic_profile = 2_i32

  integer(i32), parameter, public :: external_open_escape = 0_i32
  integer(i32), parameter, public :: external_open_potential_barrier = 1_i32

  integer(i32), parameter, public :: external_transport_none = 0_i32
  integer(i32), parameter, public :: external_transport_kinetic_1d = 1_i32
  integer(i32), parameter, public :: external_transport_unified_3d = 2_i32

  type, public :: external_boundary_contract_type
    integer(i32) :: inflow_map = external_inflow_none
    integer(i32) :: ordinary_open_model = external_open_escape
    integer(i32) :: interface_transport = external_transport_none
    logical :: queue_enabled = .false.
  end type external_boundary_contract_type

  public :: resolve_external_boundary_contract
  public :: external_boundary_owns_event

contains

  !> 互換設定群を、流入・通常open・interface輸送・queue ownershipの一意な組へ解決する。
  subroutine resolve_external_boundary_contract( &
    reservoir_potential_model, open_boundary_model, outer_model, kinetic_closure, &
    return_model, particle_transfer_mode, queue_enabled, contract, status, message &
    )
    character(len=*), intent(in) :: reservoir_potential_model
    character(len=*), intent(in) :: open_boundary_model
    character(len=*), intent(in) :: outer_model
    character(len=*), intent(in) :: kinetic_closure
    character(len=*), intent(in) :: return_model
    character(len=*), intent(in) :: particle_transfer_mode
    logical, intent(in) :: queue_enabled
    type(external_boundary_contract_type), intent(out) :: contract
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    character(len=32) :: reservoir, open_model, field_model, closure, return_map, transfer

    contract = external_boundary_contract_type()
    status = external_boundary_ok
    message = ''
    reservoir = lower_ascii(trim(reservoir_potential_model))
    open_model = lower_ascii(trim(open_boundary_model))
    field_model = lower_ascii(trim(outer_model))
    closure = lower_ascii(trim(kinetic_closure))
    return_map = lower_ascii(trim(return_model))
    transfer = lower_ascii(trim(particle_transfer_mode))

    select case (trim(open_model))
    case ('escape')
      contract%ordinary_open_model = external_open_escape
    case ('potential_barrier')
      contract%ordinary_open_model = external_open_potential_barrier
    case default
      call reject('Unknown sim.open_boundary_model.', status, message)
      return
    end select

    select case (trim(field_model))
    case ('none')
      continue
    case ('kinetic_1d')
      continue
    case ('unified_linear_response')
      continue
    case default
      call reject('Unknown outer_plasma.model.', status, message)
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
    select case (trim(transfer))
    case ('none')
      if (trim(return_map) /= 'none') then
        call reject('outer_plasma.return_model requires an interface particle transfer.', status, message)
        return
      end if
    case ('electrostatic_1d_instant_return')
      select case (trim(field_model))
      case ('kinetic_1d')
        if (trim(return_map) /= 'kinetic_1d_profile_return') then
          call reject('kinetic_1d requires kinetic_1d_profile_return.', status, message)
          return
        end if
        if (contract%inflow_map /= external_inflow_none) then
          call reject( &
            'kinetic_1d particle transfer owns inflow and cannot mix with scalar inflow correction.', &
            status, message &
            )
          return
        end if
        contract%inflow_map = external_inflow_kinetic_profile
        contract%interface_transport = external_transport_kinetic_1d
      case default
        call reject('Electrostatic 1D transfer requires kinetic_1d.', status, message)
        return
      end select
    case ('electrostatic_3d_explicit_orbit')
      if (trim(field_model) /= 'unified_linear_response' .or. &
          trim(return_map) /= 'electrostatic_3d_explicit_orbit') then
        call reject('Explicit 3D transfer requires unified_linear_response and its matching return model.', &
                    status, message)
        return
      end if
      contract%interface_transport = external_transport_unified_3d
    case default
      call reject('Unknown coupling.particle_transfer_mode.', status, message)
      return
    end select

    if (queue_enabled) then
      if (contract%interface_transport /= external_transport_kinetic_1d .or. &
          trim(closure) /= 'zhao_charge_driven') then
        call reject('The delayed outer queue requires Zhao kinetic-profile transport.', status, message)
        return
      end if
      contract%queue_enabled = .true.
    end if
  end subroutine resolve_external_boundary_contract

  !> eventに契約所有面が含まれ、かつその面がopenならouter ownershipへ渡す。
  pure logical function external_boundary_owns_event(contract, face_mask, face_bc) result(owned)
    type(external_boundary_contract_type), intent(in) :: contract
    integer(i32), intent(in) :: face_mask
    integer(i32), intent(in) :: face_bc(6)
    owned = .false.
    if (contract%interface_transport == external_transport_none) return
    owned = btest(face_mask, 5_i32) .and. face_bc(6) == bc_open
  end function external_boundary_owns_event

  pure subroutine reject(text, status, message)
    character(len=*), intent(in) :: text
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = external_boundary_invalid
    message = text
  end subroutine reject

end module bem_external_boundary_contract

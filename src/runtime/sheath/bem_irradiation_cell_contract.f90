!> 照射実験の壁・周期軸・下端モデル。汎用 kernel はこの型へ依存しない。
module bem_irradiation_cell_contract
  use bem_kinds, only: dp, i32
  use bem_types, only: bc_open, bc_periodic
  use bem_app_config_types, only: app_config
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  integer(i32), parameter, public :: irradiation_cell_ok = 0_i32
  integer(i32), parameter, public :: irradiation_cell_invalid = 1_i32

  type, public :: irradiation_cell_contract_type
    ! bem_boundary の面順: x_low, x_high, y_low, y_high, z_low, z_high。
    integer(i32) :: wall_face = 6_i32
    integer(i32) :: wall_axis = 3_i32
    integer(i32) :: outward_normal_sign = 1_i32
    real(dp) :: wall_position = 0.0_dp
    logical :: periodic_axes(3) = .false.
    character(len=32) :: lower_boundary_model = 'not_applicable'
  contains
    procedure :: wall_face_bit
    procedure :: outward_normal_velocity
    procedure :: crosses_wall_outward
    procedure :: wall_coordinate
    procedure :: wall_face_name
    procedure :: wall_area
    procedure :: wall_tangents_periodic
  end type irradiation_cell_contract_type

  public :: resolve_irradiation_cell_contract

contains

  !> preflight 済み設定から run に不変な前提を一度解決する。free-space 参照構成も保持する。
  subroutine resolve_irradiation_cell_contract(app, contract, status, message)
    type(app_config), intent(in) :: app
    type(irradiation_cell_contract_type), intent(out) :: contract
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    contract = irradiation_cell_contract_type()
    status = irradiation_cell_ok
    message = ''
    contract%wall_position = app%sim%box_max(contract%wall_axis)
    contract%periodic_axes = app%sim%bc_low == bc_periodic .and. app%sim%bc_high == bc_periodic
    contract%lower_boundary_model = lower_ascii(trim(app%periodic2%lower_boundary_model))
    ! 壁を選ぶ公開設定はない。periodic2 の方向は physics preflight と同じ x/y 周期に限る。
    if (trim(lower_ascii(app%sim%field_bc_mode)) == 'periodic2') then
      if (.not. app%sim%use_box .or. .not. all(contract%periodic_axes(1:2)) .or. &
          app%sim%bc_low(3) == bc_periodic .or. app%sim%bc_high(3) == bc_periodic) then
        status = irradiation_cell_invalid
        message = 'Irradiation cell requires a finite [domain] with x/y periodic and z nonperiodic axes.'
        return
      end if
    end if
    select case (trim(contract%lower_boundary_model))
    case ('not_applicable', 'legacy_implicit', 'e_bottom_zero', 'symmetric_vacuum', 'e_top_from_sheath')
      continue
    case default
      status = irradiation_cell_invalid
      message = 'Unsupported fields.periodic.lower_boundary_model in irradiation cell contract.'
    end select
  end subroutine resolve_irradiation_cell_contract

  pure integer(i32) function wall_face_bit(self) result(bit)
    class(irradiation_cell_contract_type), intent(in) :: self
    bit = ibset(0_i32, self%wall_face - 1_i32)
  end function wall_face_bit

  pure real(dp) function outward_normal_velocity(self, velocity) result(normal_velocity)
    class(irradiation_cell_contract_type), intent(in) :: self
    real(dp), intent(in) :: velocity(3)
    normal_velocity = real(self%outward_normal_sign, dp)*velocity(self%wall_axis)
  end function outward_normal_velocity

  pure logical function crosses_wall_outward(self, face_mask, face_bc, velocity) result(crosses)
    class(irradiation_cell_contract_type), intent(in) :: self
    integer(i32), intent(in) :: face_mask, face_bc(6)
    real(dp), intent(in) :: velocity(3)
    crosses = iand(face_mask, self%wall_face_bit()) /= 0_i32 .and. &
              face_bc(self%wall_face) == bc_open .and. self%outward_normal_velocity(velocity) > 0.0_dp
  end function crosses_wall_outward

  pure real(dp) function wall_coordinate(self) result(coordinate)
    class(irradiation_cell_contract_type), intent(in) :: self
    coordinate = self%wall_position
  end function wall_coordinate

  pure function wall_face_name(self) result(name)
    class(irradiation_cell_contract_type), intent(in) :: self
    character(len=6) :: name
    character(len=3), parameter :: axes = 'xyz'
    name = axes(self%wall_axis:self%wall_axis)//'_low'
    if (self%outward_normal_sign > 0_i32) name = axes(self%wall_axis:self%wall_axis)//'_high'
  end function wall_face_name

  pure real(dp) function wall_area(self, box_min, box_max) result(area)
    class(irradiation_cell_contract_type), intent(in) :: self
    real(dp), intent(in) :: box_min(3), box_max(3)
    area = product(pack(box_max - box_min, [1, 2, 3] /= self%wall_axis))
  end function wall_area

  pure logical function wall_tangents_periodic(self) result(periodic)
    class(irradiation_cell_contract_type), intent(in) :: self
    periodic = all(pack(self%periodic_axes, [1, 2, 3] /= self%wall_axis))
  end function wall_tangents_periodic

end module bem_irradiation_cell_contract

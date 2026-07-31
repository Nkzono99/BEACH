!> 設定型のヘルパー関数。パーサに依存せず下位層から利用可能。
module bem_config_helpers
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_config, bc_periodic
  use bem_app_config_types, only: particle_species_spec, particle_bc_inherit
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  public :: resolve_inject_face
  public :: resolve_inward_normal
  public :: resolve_particle_boundaries
  public :: particle_boundary_action_for_face
  public :: species_number_density_m3
  public :: species_temperature_k

contains

  !> 注入面識別子から法線軸と対応境界座標を返す。
  subroutine resolve_inject_face(box_min, box_max, inject_face, axis, boundary_value)
    real(dp), intent(in) :: box_min(3), box_max(3)
    character(len=*), intent(in) :: inject_face
    integer, intent(out) :: axis
    real(dp), intent(out) :: boundary_value

    select case (trim(lower_ascii(inject_face)))
    case ('x_low')
      axis = 1; boundary_value = box_min(1)
    case ('x_high')
      axis = 1; boundary_value = box_max(1)
    case ('y_low')
      axis = 2; boundary_value = box_min(2)
    case ('y_high')
      axis = 2; boundary_value = box_max(2)
    case ('z_low')
      axis = 3; boundary_value = box_min(3)
    case ('z_high')
      axis = 3; boundary_value = box_max(3)
    case default
      error stop 'Unknown particles.species.inject_face.'
    end select
  end subroutine resolve_inject_face

  !> 注入面識別子から内向き法線ベクトルを返す。
  subroutine resolve_inward_normal(inject_face, inward_normal)
    character(len=*), intent(in) :: inject_face
    real(dp), intent(out) :: inward_normal(3)

    inward_normal = 0.0d0
    select case (trim(lower_ascii(inject_face)))
    case ('x_low')
      inward_normal(1) = 1.0d0
    case ('x_high')
      inward_normal(1) = -1.0d0
    case ('y_low')
      inward_normal(2) = 1.0d0
    case ('y_high')
      inward_normal(2) = -1.0d0
    case ('z_low')
      inward_normal(3) = 1.0d0
    case ('z_high')
      inward_normal(3) = -1.0d0
    case default
      error stop 'Unknown particles.species.inject_face.'
    end select
  end subroutine resolve_inward_normal

  !> domain topology、global既定、粒子種overrideから有効な粒子境界を解決する。
  pure subroutine resolve_particle_boundaries(sim, global_low, global_high, spec, boundary_low, boundary_high)
    type(sim_config), intent(in) :: sim
    integer(i32), intent(in) :: global_low(3), global_high(3)
    type(particle_species_spec), intent(in) :: spec
    integer(i32), intent(out) :: boundary_low(3), boundary_high(3)
    integer :: axis

    boundary_low = sim%bc_low
    boundary_high = sim%bc_high
    do axis = 1, 3
      if (boundary_low(axis) /= bc_periodic .and. global_low(axis) /= particle_bc_inherit) then
        boundary_low(axis) = global_low(axis)
      end if
      if (boundary_high(axis) /= bc_periodic .and. global_high(axis) /= particle_bc_inherit) then
        boundary_high(axis) = global_high(axis)
      end if
      if (boundary_low(axis) /= bc_periodic .and. spec%boundary_low(axis) /= particle_bc_inherit) then
        boundary_low(axis) = spec%boundary_low(axis)
      end if
      if (boundary_high(axis) /= bc_periodic .and. spec%boundary_high(axis) /= particle_bc_inherit) then
        boundary_high(axis) = spec%boundary_high(axis)
      end if
    end do
  end subroutine resolve_particle_boundaries

  !> 6面名に対応する有効粒子境界作用を返す。
  pure function particle_boundary_action_for_face(boundary_low, boundary_high, face) result(action)
    integer(i32), intent(in) :: boundary_low(3), boundary_high(3)
    character(len=*), intent(in) :: face
    integer(i32) :: action

    select case (trim(lower_ascii(face)))
    case ('x_low')
      action = boundary_low(1)
    case ('x_high')
      action = boundary_high(1)
    case ('y_low')
      action = boundary_low(2)
    case ('y_high')
      action = boundary_high(2)
    case ('z_low')
      action = boundary_low(3)
    case ('z_high')
      action = boundary_high(3)
    case default
      action = particle_bc_inherit
    end select
  end function particle_boundary_action_for_face

  !> 粒子種設定から有効粒子数密度 [1/m^3] を計算する。
  pure function species_number_density_m3(spec) result(number_density_m3)
    type(particle_species_spec), intent(in) :: spec
    real(dp) :: number_density_m3

    number_density_m3 = spec%number_density_m3
    if (spec%has_number_density_cm3) number_density_m3 = spec%number_density_cm3*1.0d6
  end function species_number_density_m3

  !> 粒子種設定から有効温度 [K] を計算する。
  pure function species_temperature_k(spec) result(temperature_k)
    type(particle_species_spec), intent(in) :: spec
    real(dp) :: temperature_k

    temperature_k = spec%temperature_k
    if (spec%has_temperature_ev) temperature_k = spec%temperature_ev*1.160451812d4
  end function species_temperature_k

end module bem_config_helpers

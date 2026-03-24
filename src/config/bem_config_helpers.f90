!> 設定型のヘルパー関数。パーサに依存せず下位層から利用可能。
module bem_config_helpers
  use bem_kinds, only: dp
  use bem_app_config_types, only: particle_species_spec
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  public :: resolve_inject_face
  public :: resolve_inward_normal
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

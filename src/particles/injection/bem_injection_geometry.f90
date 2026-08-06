!> box面注入で共有する面軸と開口geometryの解決。
module bem_injection_geometry
  use bem_kinds, only: dp
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  public :: compute_face_area_from_bounds
  public :: resolve_face_axes
  public :: resolve_face_geometry

contains

  !> 注入面上の矩形開口から有効面積[m^2]を返す。
  !! @param[in] inject_face 注入面識別子（`x_low/x_high/y_low/y_high/z_low/z_high`）。
  !! @param[in] pos_low 開口領域の下限座標 `(x,y,z)` [m]。
  !! @param[in] pos_high 開口領域の上限座標 `(x,y,z)` [m]。
  !! @return area 注入開口の有効面積 [m^2]。
  real(dp) function compute_face_area_from_bounds(inject_face, pos_low, pos_high) result(area)
    character(len=*), intent(in) :: inject_face
    real(dp), intent(in) :: pos_low(3), pos_high(3)
    integer :: axis_t1, axis_t2

    call resolve_face_axes(inject_face, axis_t1, axis_t2)
    area = (pos_high(axis_t1) - pos_low(axis_t1))*(pos_high(axis_t2) - pos_low(axis_t2))
  end function compute_face_area_from_bounds

  !> 注入面名から接線2軸を返す。
  !! @param[in] inject_face 注入面識別子（`x_low/x_high/y_low/y_high/z_low/z_high`）。
  !! @param[out] axis_t1 注入面の第1接線軸インデックス（1:x, 2:y, 3:z）。
  !! @param[out] axis_t2 注入面の第2接線軸インデックス（1:x, 2:y, 3:z）。
  subroutine resolve_face_axes(inject_face, axis_t1, axis_t2)
    character(len=*), intent(in) :: inject_face
    integer, intent(out) :: axis_t1, axis_t2

    select case (trim(adjustl(inject_face)))
    case ('x_low', 'x_high')
      axis_t1 = 2
      axis_t2 = 3
    case ('y_low', 'y_high')
      axis_t1 = 3
      axis_t2 = 1
    case ('z_low', 'z_high')
      axis_t1 = 1
      axis_t2 = 2
    case default
      error stop "unknown inject_face"
    end select
  end subroutine resolve_face_axes

  !> 注入面名から法線方向の幾何情報を返す。
  !! @param[in] box_min シミュレーションボックス下限座標 `(x,y,z)` [m]。
  !! @param[in] box_max シミュレーションボックス上限座標 `(x,y,z)` [m]。
  !! @param[in] inject_face 注入面識別子（`x_low/x_high/y_low/y_high/z_low/z_high`）。
  !! @param[out] axis_n 法線方向軸インデックス（省略可）。
  !! @param[out] boundary_value 注入面の境界座標値 [m]（省略可）。
  !! @param[out] inward_normal 注入面の内向き法線ベクトル（省略可）。
  subroutine resolve_face_geometry(box_min, box_max, inject_face, axis_n, boundary_value, inward_normal)
    real(dp), intent(in) :: box_min(3), box_max(3)
    character(len=*), intent(in) :: inject_face
    integer, intent(out), optional :: axis_n
    real(dp), intent(out), optional :: boundary_value
    real(dp), intent(out), optional :: inward_normal(3)

    integer :: axis_local
    real(dp) :: boundary_local, normal_local(3)

    normal_local = 0.0_dp
    select case (trim(adjustl(inject_face)))
    case ('x_low')
      axis_local = 1
      boundary_local = box_min(1)
      normal_local(1) = 1.0_dp
    case ('x_high')
      axis_local = 1
      boundary_local = box_max(1)
      normal_local(1) = -1.0_dp
    case ('y_low')
      axis_local = 2
      boundary_local = box_min(2)
      normal_local(2) = 1.0_dp
    case ('y_high')
      axis_local = 2
      boundary_local = box_max(2)
      normal_local(2) = -1.0_dp
    case ('z_low')
      axis_local = 3
      boundary_local = box_min(3)
      normal_local(3) = 1.0_dp
    case ('z_high')
      axis_local = 3
      boundary_local = box_max(3)
      normal_local(3) = -1.0_dp
    case default
      error stop "unknown inject_face"
    end select

    if (present(axis_n)) axis_n = axis_local
    if (present(boundary_value)) boundary_value = boundary_local
    if (present(inward_normal)) inward_normal = normal_local
  end subroutine resolve_face_geometry

end module bem_injection_geometry

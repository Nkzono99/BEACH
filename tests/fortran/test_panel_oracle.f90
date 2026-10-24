!> 三角形 P0 panel の解析核を全 panel で直接足し合わせる、field solver test 用の独立 oracle。
module test_panel_oracle
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_kernel, only: panel_potential_field, panel_side_principal_value
  implicit none
  private

  public :: direct_panel_sum
  public :: panel_sum_field
  public :: panel_sum_potential

contains

  !> 各 panel の解析核（面上は主値）を直接和で足し、SI 単位の電位と電場を返す。
  !! @param[in] v0 panel 頂点 0 (3, n)。
  !! @param[in] v1 panel 頂点 1 (3, n)。
  !! @param[in] v2 panel 頂点 2 (3, n)。
  !! @param[in] q panel 電荷 [C] (n)。
  !! @param[in] target 評価点 [m]。
  !! @param[out] potential 電位 [V]。
  !! @param[out] field 電場 [V/m]。
  subroutine direct_panel_sum(v0, v1, v2, q, target, potential, field)
    real(dp), intent(in) :: v0(:, :), v1(:, :), v2(:, :), q(:), target(3)
    real(dp), intent(out) :: potential, field(3)
    type(panel_geometry_type) :: geometry
    real(dp) :: source_potential, source_field(3)
    integer(i32) :: element, status

    potential = 0.0_dp
    field = 0.0_dp
    do element = 1_i32, int(size(q), i32)
      call init_panel_geometry(v0(:, element), v1(:, element), v2(:, element), geometry, status)
      if (status /= panel_geometry_ok) error stop 'direct panel oracle received an invalid panel'
      call panel_potential_field( &
        geometry, q(element), target, panel_side_principal_value, source_potential, source_field &
        )
      potential = potential + source_potential
      field = field + source_field
    end do
  end subroutine direct_panel_sum

  !> mesh の全 panel による電場 [V/m]。
  function panel_sum_field(mesh, target) result(field)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: target(3)
    real(dp) :: field(3), potential

    call direct_panel_sum(mesh%v0, mesh%v1, mesh%v2, mesh%q_elem, target, potential, field)
  end function panel_sum_field

  !> mesh の全 panel による電位 [V]。
  function panel_sum_potential(mesh, target) result(potential)
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: target(3)
    real(dp) :: field(3), potential

    call direct_panel_sum(mesh%v0, mesh%v1, mesh%v2, mesh%q_elem, target, potential, field)
  end function panel_sum_potential

end module test_panel_oracle

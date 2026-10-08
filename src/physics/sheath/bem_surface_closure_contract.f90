!> シミュレータが外部の表面電流モデルから受け取るモデル非依存の境界契約。
module bem_surface_closure_contract
  use bem_kinds, only: dp, i32
  implicit none
  private

  type, public :: surface_closure_contract_type
    logical :: active = .false.
    logical, allocatable :: has_absorbed_target(:)
    logical, allocatable :: has_emission_target(:)
    logical, allocatable :: has_escape_target(:)
    logical, allocatable :: has_inflow_kinetic_map(:)
    logical, allocatable :: has_outflow_kinetic_barrier(:)
    logical, allocatable :: has_inflow_number_flux(:)
    real(dp), allocatable :: absorbed_current_a(:)
    real(dp), allocatable :: emission_current_a(:)
    real(dp), allocatable :: escaped_particle_current_a(:)
    real(dp), allocatable :: inflow_reservoir_potential_v(:)
    real(dp), allocatable :: inflow_access_potential_v(:)
    integer(i32), allocatable :: inflow_kinetic_face(:)
    real(dp), allocatable :: outflow_barrier_potential_v(:)
    integer(i32), allocatable :: outflow_barrier_face(:)
    real(dp), allocatable :: inflow_number_flux_m2_s(:)
    !> z-high面の水平平均電位を外部シース解の壁電位へ固定する場合だけ有効。
    logical :: has_plane_gauge = .false.
    real(dp) :: plane_gauge_potential_v = 0.0_dp
    !> 外部障壁で戻る粒子を周期セル内の一様な位置へ戻す。
    logical :: outer_return_cell_uniform = .false.
  end type surface_closure_contract_type

end module bem_surface_closure_contract

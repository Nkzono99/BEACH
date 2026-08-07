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
    real(dp), allocatable :: absorbed_current_a(:)
    real(dp), allocatable :: emission_current_a(:)
    real(dp), allocatable :: escaped_particle_current_a(:)
    real(dp), allocatable :: inflow_reservoir_potential_v(:)
    real(dp), allocatable :: inflow_access_potential_v(:)
    integer(i32), allocatable :: inflow_kinetic_face(:)
    real(dp), allocatable :: outflow_barrier_potential_v(:)
    integer(i32), allocatable :: outflow_barrier_face(:)
  end type surface_closure_contract_type

end module bem_surface_closure_contract

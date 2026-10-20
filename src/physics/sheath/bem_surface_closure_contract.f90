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
    !> 流入写像の0 V貯留部にあるMaxwellianの密度。0以下ならspecies設定の密度を使う。
    real(dp), allocatable :: inflow_reservoir_density_m3(:)
    integer(i32), allocatable :: inflow_kinetic_face(:)
    real(dp), allocatable :: outflow_barrier_potential_v(:)
    integer(i32), allocatable :: outflow_barrier_face(:)
    real(dp), allocatable :: inflow_number_flux_m2_s(:)
    !> この species の脱出・帰還目標から、そのバッチで追跡した脱出電荷を差し引く species（0 なら無し）。
    !! 外部根の脱出束が複数の放出電子種の合算なら、追跡した側に充てた残りをこの species の目標にする。
    integer(i32), allocatable :: escape_partner_species(:)
    !> 放出目標が、根が予測した放出ではなく、根に与えた外側源（追跡した放出の窓平均）である種（PE なしの SEE）。
    !! raw channel が空のバッチでも止まらず、帰還目標は放出へ畳み、放出がなければ根の脱出分を未閉として記録する。
    logical, allocatable :: emission_target_from_tracked_source(:)
    !> この species の脱出・帰還・放出目標から差し引く、セル内で反射して H を抜けた上流電子の species（0 なら無し）。
    !! 外部根の源に反射電子を含めると根の脱出・帰還は放出電子と反射電子の合算なので、追跡で決まる反射分を先に充て、
    !! 残りをこの species に固定する。
    integer(i32), allocatable :: reflected_electron_species(:)
    !> 吸収目標に、この species 自身が外部障壁を越えて脱出した追跡電荷を足し戻す（反射して出た上流電子は壁に届かない）。
    logical, allocatable :: absorbed_target_credits_tracked_escape(:)
    !> 照射セル契約の壁面の平均電位を外部シース解の壁電位へ固定する場合だけ有効。
    logical :: has_plane_gauge = .false.
    real(dp) :: plane_gauge_potential_v = 0.0_dp
    real(dp) :: top_field_v_m = 0.0_dp
    !> 外部障壁で戻る粒子を周期セル内の一様な位置へ戻す。
    logical :: outer_return_cell_uniform = .false.
  end type surface_closure_contract_type

end module bem_surface_closure_contract

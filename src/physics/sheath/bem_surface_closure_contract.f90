!> シミュレータが外部の表面電流モデルから受け取るモデル非依存の境界契約。
module bem_surface_closure_contract
  use bem_kinds, only: dp, i32
  use bem_runtime_modes, only: closure_mode_none, closure_mode_target, closure_mode_escape_partner, &
                               closure_mode_tracked_source
  implicit none
  private

  public :: species_closure_mode

  !> 1 species の固定電流目標、closure の経路、外部 kinetic 境界写像。既定値は「その channel を使わない」。
  type, public :: species_closure_channel_type
    !> 固定電流の吸収・放出反作用・脱出目標 [A]。has_* が偽なら使わない。脱出は帯電電流と逆符号。
    logical :: has_absorbed_target = .false.
    logical :: has_emission_target = .false.
    logical :: has_escape_target = .false.
    real(dp) :: absorbed_current_a = 0.0_dp
    real(dp) :: emission_current_a = 0.0_dp
    real(dp) :: escaped_particle_current_a = 0.0_dp
    !> この species の脱出・帰還目標から、そのバッチで追跡した脱出電荷を差し引く species（0 なら無し）。
    !! 外部根の脱出束が複数の放出電子種の合算なら、追跡した側に充てた残りをこの species の目標にする。
    integer(i32) :: escape_partner_species = 0_i32
    !> 放出目標が、根が予測した放出ではなく、根に与えた外側源（追跡した放出の窓平均）である種（PE なしの SEE）。
    !! raw channel が空のバッチでも止まらず、帰還目標は放出へ畳み、放出がなければ根の脱出分を未閉として記録する。
    logical :: emission_target_from_tracked_source = .false.
    !> この species の脱出・帰還・放出目標から差し引く、セル内で反射して H を抜けた上流電子の species（0 なら無し）。
    !! 外部根の源に反射電子を含めると根の脱出・帰還は放出電子と反射電子の合算なので、追跡で決まる反射分を先に充て、
    !! 残りをこの species に固定する。
    integer(i32) :: reflected_electron_species = 0_i32
    !> 吸収目標に、この species 自身が外部障壁を越えて脱出した追跡電荷を足し戻す（反射して出た上流電子は壁に届かない）。
    logical :: absorbed_target_credits_tracked_escape = .false.
    !> kinetic_face での流入写像（上流 0 V の貯留部からの流入）と外部障壁（外向きに抜けた粒子の戻り）の有無。
    logical :: has_inflow_kinetic_map = .false.
    logical :: has_outflow_kinetic_barrier = .false.
    !> 上流 0 V と kinetic_face の間でこの species が越える経路の極値電位 [V]。
    !! 流入写像の access 電位と外部障壁の電位を兼ねる（同じ外部シース経路の値なので、両者は一致する）。
    real(dp) :: kinetic_path_potential_v = 0.0_dp
    !> 流入写像の0 V貯留部にあるMaxwellianの密度 [m^-3]。0以下ならspecies設定の密度を使う。
    real(dp) :: inflow_reservoir_density_m3 = 0.0_dp
  end type species_closure_channel_type

  type, public :: surface_closure_contract_type
    logical :: active = .false.
    !> species 別の channel。大きさは粒子種数。
    type(species_closure_channel_type), allocatable :: species(:)
    !> 流入写像と外部障壁が働く箱の面（1..6、0 なら無し）。全 species で共通の照射セル壁面。
    integer(i32) :: kinetic_face = 0_i32
    !> 照射セル契約の壁面の平均電位を外部シース解の壁電位へ固定する場合だけ有効。
    logical :: has_plane_gauge = .false.
    real(dp) :: plane_gauge_potential_v = 0.0_dp
    real(dp) :: top_field_v_m = 0.0_dp
    !> 外部障壁で戻る粒子を周期セル内の一様な位置へ戻す。
    logical :: outer_return_cell_uniform = .false.
  end type surface_closure_contract_type

contains

  !> 種 species_idx の固定電流 closure 方式（bem_runtime_modes の closure_mode_*）。
  !! 設定が fixed_current でも外部根の目標もなければ none。方式は contract のフラグから導くので、
  !! 方式とフラグが食い違うことはない。
  pure integer(i32) function species_closure_mode(contract, species_idx, configured_fixed_current) result(mode)
    class(surface_closure_contract_type), intent(in) :: contract
    integer(i32), intent(in) :: species_idx
    logical, intent(in) :: configured_fixed_current

    associate (channel => contract%species(species_idx))
      if (.not. (configured_fixed_current .or. channel%has_absorbed_target .or. &
                 channel%has_emission_target .or. channel%has_escape_target)) then
        mode = closure_mode_none
      else if (channel%emission_target_from_tracked_source) then
        mode = closure_mode_tracked_source
      else if (channel%escape_partner_species > 0_i32 .or. channel%reflected_electron_species > 0_i32) then
        mode = closure_mode_escape_partner
      else
        mode = closure_mode_target
      end if
    end associate
  end function species_closure_mode

end module bem_surface_closure_contract

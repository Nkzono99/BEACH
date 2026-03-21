!> `mesh_type` や `sim_config` に依存しない Coulomb FMM コア API。
module bem_coulomb_fmm_core
  use bem_kinds, only: dp
  use bem_coulomb_fmm_types, only: fmm_options_type, fmm_plan_type, fmm_state_type
  implicit none
  private

  public :: fmm_options_type
  public :: fmm_plan_type
  public :: fmm_state_type
  public :: build_plan
  public :: update_state
  public :: eval_points
  public :: eval_point
  public :: eval_potential_points
  public :: eval_potential_point
  public :: destroy_plan
  public :: destroy_state

  interface
    !> FMM 計画を構築する。
    !! @param[out] plan 構築先の計画。
    !! @param[in] src_pos ソース位置 `(3,n)` [m]。
    !! @param[in] options FMM 設定。
    module subroutine build_plan(plan, src_pos, options)
      type(fmm_plan_type), intent(inout) :: plan
      real(dp), intent(in) :: src_pos(:, :)
      type(fmm_options_type), intent(in) :: options
    end subroutine build_plan

    !> ソース電荷から FMM state を更新する。
    !! @param[in] plan 構築済みの FMM 計画。
    !! @param[inout] state 更新対象の state。
    !! @param[in] src_q ソース電荷量。
    module subroutine update_state(plan, state, src_q)
      type(fmm_plan_type), intent(in) :: plan
      type(fmm_state_type), intent(inout) :: state
      real(dp), intent(in) :: src_q(:)
    end subroutine update_state

    !> 複数の評価点で電場を計算する。
    !! @param[in] plan 構築済みの FMM 計画。
    !! @param[inout] state 評価に使う FMM state。
    !! @param[in] target_pos 評価点位置 `(3,m)` [m]。
    !! @param[out] e 電場ベクトル `(3,m)` [V/m]。
    module subroutine eval_points(plan, state, target_pos, e)
      type(fmm_plan_type), intent(in) :: plan
      type(fmm_state_type), intent(inout) :: state
      real(dp), intent(in) :: target_pos(:, :)
      real(dp), intent(out) :: e(:, :)
    end subroutine eval_points

    !> 1 点で電場を計算する。
    !! @param[in] plan 構築済みの FMM 計画。
    !! @param[inout] state 評価に使う FMM state。
    !! @param[in] r 評価点位置 `(x,y,z)` [m]。
    !! @param[out] e 電場ベクトル `(x,y,z)` [V/m]。
    module subroutine eval_point(plan, state, r, e)
      type(fmm_plan_type), intent(in) :: plan
      type(fmm_state_type), intent(inout) :: state
      real(dp), intent(in) :: r(3)
      real(dp), intent(out) :: e(3)
    end subroutine eval_point

    !> 複数の評価点で電位を計算する。
    !! @param[in] plan 構築済みの FMM 計画。
    !! @param[inout] state 評価に使う FMM state。
    !! @param[in] target_pos 評価点位置 `(3,m)` [m]。
    !! @param[out] phi 電位 `(m)` [V]（`k_coulomb` は含まない）。
    module subroutine eval_potential_points(plan, state, target_pos, phi)
      type(fmm_plan_type), intent(in) :: plan
      type(fmm_state_type), intent(inout) :: state
      real(dp), intent(in) :: target_pos(:, :)
      real(dp), intent(out) :: phi(:)
    end subroutine eval_potential_points

    !> 1 点で電位を計算する。
    !! @param[in] plan 構築済みの FMM 計画。
    !! @param[inout] state 評価に使う FMM state。
    !! @param[in] r 評価点位置 `(x,y,z)` [m]。
    !! @param[out] phi 電位 [V]（`k_coulomb` は含まない）。
    module subroutine eval_potential_point(plan, state, r, phi)
      type(fmm_plan_type), intent(in) :: plan
      type(fmm_state_type), intent(inout) :: state
      real(dp), intent(in) :: r(3)
      real(dp), intent(out) :: phi
    end subroutine eval_potential_point

    !> FMM 計画を解放する。
    !! @param[inout] plan 解放対象の計画。
    module subroutine destroy_plan(plan)
      type(fmm_plan_type), intent(inout) :: plan
    end subroutine destroy_plan

    !> FMM state を解放する。
    !! @param[inout] state 解放対象の state。
    module subroutine destroy_state(state)
      type(fmm_state_type), intent(inout) :: state
    end subroutine destroy_state
  end interface

end module bem_coulomb_fmm_core

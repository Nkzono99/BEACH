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
  public :: destroy_plan
  public :: destroy_state

  interface
    module subroutine build_plan(plan, src_pos, options)
      type(fmm_plan_type), intent(inout) :: plan
      real(dp), intent(in) :: src_pos(:, :)
      type(fmm_options_type), intent(in) :: options
    end subroutine build_plan

    module subroutine update_state(plan, state, src_q)
      type(fmm_plan_type), intent(in) :: plan
      type(fmm_state_type), intent(inout) :: state
      real(dp), intent(in) :: src_q(:)
    end subroutine update_state

    module subroutine eval_points(plan, state, target_pos, e)
      type(fmm_plan_type), intent(in) :: plan
      type(fmm_state_type), intent(inout) :: state
      real(dp), intent(in) :: target_pos(:, :)
      real(dp), intent(out) :: e(:, :)
    end subroutine eval_points

    module subroutine eval_point(plan, state, r, e)
      type(fmm_plan_type), intent(in) :: plan
      type(fmm_state_type), intent(inout) :: state
      real(dp), intent(in) :: r(3)
      real(dp), intent(out) :: e(3)
    end subroutine eval_point

    module subroutine destroy_plan(plan)
      type(fmm_plan_type), intent(inout) :: plan
    end subroutine destroy_plan

    module subroutine destroy_state(state)
      type(fmm_state_type), intent(inout) :: state
    end subroutine destroy_state
  end interface

end module bem_coulomb_fmm_core

!> `bem_field_solver` と FMM core の間のパネル幾何・電荷状態の更新を担う。
submodule(bem_field_solver) bem_field_solver_fmm
  use bem_coulomb_fmm_core, only: build_panel_plan, update_state, destroy_plan, destroy_state
  implicit none
contains

  module procedure refresh_fmm_solver
  real(dp), allocatable :: panel_v0(:, :), panel_v1(:, :), panel_v2(:, :)

  self%nelem = mesh%nelem
  if (mesh%nelem <= 0_i32) then
    call destroy_plan(self%fmm_core_plan)
    call destroy_state(self%fmm_core_state)
    self%fmm_core_ready = .false.
    return
  end if

  if (.not. self%fmm_core_plan%built .or. self%fmm_core_plan%nsrc /= mesh%nelem) then
    call destroy_plan(self%fmm_core_plan)
    call destroy_state(self%fmm_core_state)
    allocate (panel_v0(3, mesh%nelem), panel_v1(3, mesh%nelem), panel_v2(3, mesh%nelem))
    panel_v0 = (mesh%v0 - spread(self%field_origin, 2, mesh%nelem))*self%field_inv_length_scale
    panel_v1 = (mesh%v1 - spread(self%field_origin, 2, mesh%nelem))*self%field_inv_length_scale
    panel_v2 = (mesh%v2 - spread(self%field_origin, 2, mesh%nelem))*self%field_inv_length_scale
    call build_panel_plan(self%fmm_core_plan, panel_v0, panel_v1, panel_v2, self%fmm_core_options)
  end if

  call update_state(self%fmm_core_plan, self%fmm_core_state, mesh%q_elem)
  self%fmm_core_ready = self%fmm_core_plan%built .and. self%fmm_core_state%ready
  end procedure refresh_fmm_solver

end submodule bem_field_solver_fmm

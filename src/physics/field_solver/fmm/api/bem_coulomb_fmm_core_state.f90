!> `bem_coulomb_fmm_core` の state 更新 API ラッパ。
submodule(bem_coulomb_fmm_core) bem_coulomb_fmm_core_state
  use bem_coulomb_fmm_state_ops, only: core_update_state_impl, core_destroy_state_impl
  implicit none
contains

  module procedure update_state
  call core_update_state_impl(plan, state, src_q)
  end procedure update_state

  module procedure destroy_state
  call core_destroy_state_impl(state)
  end procedure destroy_state

end submodule bem_coulomb_fmm_core_state

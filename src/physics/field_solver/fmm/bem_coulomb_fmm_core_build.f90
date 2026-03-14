!> `bem_coulomb_fmm_core` の plan 構築 API ラッパ。
submodule (bem_coulomb_fmm_core) bem_coulomb_fmm_core_build
  use bem_coulomb_fmm_plan_ops, only: core_build_plan_impl, core_destroy_plan_impl
  implicit none
contains

  module procedure build_plan
    call core_build_plan_impl(plan, src_pos, options)
  end procedure build_plan

  module procedure destroy_plan
    call core_destroy_plan_impl(plan)
  end procedure destroy_plan

end submodule bem_coulomb_fmm_core_build

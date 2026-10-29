!> `bem_coulomb_fmm_core` の plan 構築 API ラッパ。
submodule(bem_coulomb_fmm_core) bem_coulomb_fmm_core_build
  use bem_coulomb_fmm_plan_ops, only: core_build_plan_impl, core_build_panel_plan_impl, core_prepare_panel_sources_impl, &
                                      core_complete_plan_impl, core_destroy_plan_impl
  implicit none
contains

  module procedure build_plan
  call core_build_plan_impl(plan, src_pos, options)
  end procedure build_plan

  module procedure build_panel_plan
  call core_build_panel_plan_impl(plan, v0, v1, v2, options)
  end procedure build_panel_plan

  module procedure prepare_panel_sources
  call core_prepare_panel_sources_impl(plan, v0, v1, v2, options)
  end procedure prepare_panel_sources

  module procedure complete_plan
  call core_complete_plan_impl(plan)
  end procedure complete_plan

  module procedure destroy_plan
  call core_destroy_plan_impl(plan)
  end procedure destroy_plan

end submodule bem_coulomb_fmm_core_build

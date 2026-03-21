!> `bem_coulomb_fmm_core` の評価 API ラッパ。
submodule(bem_coulomb_fmm_core) bem_coulomb_fmm_core_eval
  use bem_coulomb_fmm_eval_ops, only: core_eval_points_impl, core_eval_point_impl, &
                                      core_eval_potential_points_impl, core_eval_potential_point_impl
  implicit none
contains

  module procedure eval_points
  call core_eval_points_impl(plan, state, target_pos, e)
  end procedure eval_points

  module procedure eval_point
  call core_eval_point_impl(plan, state, r, e)
  end procedure eval_point

  module procedure eval_potential_points
  call core_eval_potential_points_impl(plan, state, target_pos, phi)
  end procedure eval_potential_points

  module procedure eval_potential_point
  call core_eval_potential_point_impl(plan, state, r, phi)
  end procedure eval_potential_point

end submodule bem_coulomb_fmm_core_eval

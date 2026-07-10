!> 同一時刻の粒子状態から、空間電場を中点評価した1ステップ候補を構築する。
module bem_particle_stepper
  use bem_kinds, only: dp
  use bem_types, only: mesh_type, sim_config
  use bem_field_solver, only: field_solver_type
  use bem_pusher, only: boris_push
  implicit none
  private

  public :: build_particle_step_candidate

contains

  !> 予測中点の電場と一様磁場を使い、次時刻の位置・速度候補を返す。
  subroutine build_particle_step_candidate( &
    mesh, sim, field_solver, bfield, x0, v0, q, m, dt, x1, v1 &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(field_solver_type), intent(inout) :: field_solver
    real(dp), intent(in) :: bfield(3), x0(3), v0(3), q, m, dt
    real(dp), intent(out) :: x1(3), v1(3)
    real(dp) :: x_mid(3), e_mid(3)

    x_mid = x0 + 0.5d0*v0*dt
    call field_solver%eval_e(mesh, x_mid, e_mid)
    e_mid = e_mid + sim%e0
    call boris_push(x0, v0, q, m, dt, e_mid, bfield, x1, v1)
  end subroutine build_particle_step_candidate

end module bem_particle_stepper

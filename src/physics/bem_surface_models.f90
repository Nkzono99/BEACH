!> 表面モデルごとの電荷更新後処理を扱うモジュール。
module bem_surface_models
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, sim_config, surface_model_conductor
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  implicit none
  private

  !> `mesh_id` ごとの浮遊導体を、snapshot が返す電位で等電位化する solver。
  !! 影響行列は幾何と場の境界条件だけで決まるので初期化時に一度組み、batch ごとには
  !! 固定電荷が作る電位だけを snapshot から取り直す。free-space と periodic2 の両方で同じ経路を使う。
  type, public :: floating_conductor_solver_type
    integer(i32) :: ncond = 0_i32
    integer(i32) :: ngroup = 0_i32
    integer(i32), allocatable :: elems(:)
    integer(i32), allocatable :: group(:)
    integer(i32), allocatable :: group_mesh_id(:)
    !> conductor 要素 j の単位電荷が要素 i の重心に作る電位を `k_coulomb` で割った係数。
    real(dp), allocatable :: influence(:, :)
    !> 等電位条件と group 総電荷保存を並べた正方系の LU 分解。
    real(dp), allocatable :: factors(:, :)
    integer(i32), allocatable :: pivot(:)
  contains
    procedure :: init => init_floating_conductor_solver
    procedure :: active => floating_conductor_solver_active
    procedure :: relax => relax_floating_conductor_charges
  end type floating_conductor_solver_type

  interface
    !> conductor 要素を列挙し、snapshot で影響行列を組んで LU 分解する。
    !! `mesh%q_elem` は一時的に書き換えるが、戻り時には入力値へ戻し snapshot も refresh し直す。
    module subroutine init_floating_conductor_solver(self, mesh, sim, snapshot)
      class(floating_conductor_solver_type), intent(out) :: self
      type(mesh_type), intent(inout) :: mesh
      type(sim_config), intent(in) :: sim
      type(electrostatic_snapshot_type), intent(inout) :: snapshot
    end subroutine init_floating_conductor_solver

    !> `mesh%q_elem` の conductor 要素を、group 総電荷を保って等電位になるよう置き換える。
    !! snapshot は入力電荷で refresh される。置換後の電荷では refresh しないので、呼び出し側が次の
    !! batch 前に refresh する。
    module subroutine relax_floating_conductor_charges(self, mesh, sim, snapshot)
      class(floating_conductor_solver_type), intent(in) :: self
      type(mesh_type), intent(inout) :: mesh
      type(sim_config), intent(in) :: sim
      type(electrostatic_snapshot_type), intent(inout) :: snapshot
    end subroutine relax_floating_conductor_charges
  end interface

contains

  !> conductor 要素が 1 つ以上あるか。
  pure logical function floating_conductor_solver_active(self) result(active)
    class(floating_conductor_solver_type), intent(in) :: self
    active = self%ncond > 0_i32
  end function floating_conductor_solver_active

end module bem_surface_models

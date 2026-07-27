!> 表面モデルごとの電荷更新後処理を扱うモジュール。
module bem_surface_models
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, surface_model_conductor
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  public :: apply_surface_model_charge_relaxation

  interface
    !> conductor 要素を mesh_id ごとの浮遊導体として等電位化する。
    module subroutine relax_floating_conductor_charges(mesh, external_e, ncond)
      type(mesh_type), intent(inout) :: mesh
      real(dp), intent(in) :: external_e(3)
      integer(i32), intent(in) :: ncond
    end subroutine relax_floating_conductor_charges
  end interface

contains

  !> 表面モデルに応じて、電荷堆積後の要素電荷を緩和する。
  !! conductor は mesh_id ごとの浮遊導体として、総電荷を保存しながら等電位化する。
  subroutine apply_surface_model_charge_relaxation(mesh, external_e, field_bc_mode)
    type(mesh_type), intent(inout) :: mesh
    real(dp), intent(in) :: external_e(3)
    character(len=*), intent(in), optional :: field_bc_mode

    integer(i32) :: ncond

    if (.not. allocated(mesh%elem_surface_model)) return
    ncond = int(count(mesh%elem_surface_model == surface_model_conductor), kind=i32)
    if (ncond <= 0_i32) return
    call validate_conductor_field_bc(field_bc_mode)

    call relax_floating_conductor_charges(mesh, external_e, ncond)
  end subroutine apply_surface_model_charge_relaxation

  !> conductor 再配分が対応している場境界条件か検証する。
  subroutine validate_conductor_field_bc(field_bc_mode)
    character(len=*), intent(in), optional :: field_bc_mode
    character(len=16) :: mode

    mode = 'free'
    if (present(field_bc_mode)) mode = lower_ascii(trim(field_bc_mode))
    if (trim(mode) /= 'free') then
      error stop 'surface_model="conductor" currently requires sim.field_bc_mode="free".'
    end if
  end subroutine validate_conductor_field_bc

end module bem_surface_models

!> 設定の文字列モードを実行時に一度だけ解決した整数 enum。
!! 設定キー・summary・エラーメッセージの文字列はそのままにし、粒子・評価ループは解決済みの値だけを比べる。
!! 文字列は parser/preflight が小文字化・検証済みなので、ここでは大文字小文字を区別しない比較で写すだけにする。
module bem_runtime_modes
  use bem_kinds, only: i32
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  !> 場ソルバの評価方式（fields.solver.method を要素数で解決した後の値）。
  integer(i32), parameter, public :: field_solver_mode_direct = 0_i32
  integer(i32), parameter, public :: field_solver_mode_treecode = 1_i32
  integer(i32), parameter, public :: field_solver_mode_fmm = 2_i32

  !> particles.species.charging.closure。
  integer(i32), parameter, public :: surface_charge_closure_explicit = 0_i32
  integer(i32), parameter, public :: surface_charge_closure_fixed_current = 1_i32
  integer(i32), parameter, public :: surface_charge_closure_neutral_return = 2_i32
  integer(i32), parameter, public :: surface_charge_closure_unknown = -1_i32

  public :: field_solver_mode_id
  public :: surface_charge_closure_id

contains

  !> 'direct' / 'treecode' / 'fmm' を enum へ写す。それ以外は direct（場ソルバの既定と同じ）。
  pure integer(i32) function field_solver_mode_id(mode) result(id)
    character(len=*), intent(in) :: mode

    select case (trim(lower_ascii(mode)))
    case ('fmm')
      id = field_solver_mode_fmm
    case ('treecode')
      id = field_solver_mode_treecode
    case default
      id = field_solver_mode_direct
    end select
  end function field_solver_mode_id

  !> 'explicit' / 'fixed_current' / 'neutral_return' を enum へ写す。未知の値は unknown。
  pure integer(i32) function surface_charge_closure_id(closure) result(id)
    character(len=*), intent(in) :: closure

    select case (trim(lower_ascii(closure)))
    case ('explicit')
      id = surface_charge_closure_explicit
    case ('fixed_current')
      id = surface_charge_closure_fixed_current
    case ('neutral_return')
      id = surface_charge_closure_neutral_return
    case default
      id = surface_charge_closure_unknown
    end select
  end function surface_charge_closure_id

end module bem_runtime_modes

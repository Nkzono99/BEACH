!> periodic2 の 2 軸周期設定（周期軸と周期長）を検証付きで解決する。
!! mesh の collision 準備と field solver の初期化が同じ規則を使う。
module bem_periodic2_axes
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_config, bc_periodic
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  public :: resolve_periodic2_axes

contains

  !> fields.boundary="periodic2" なら周期軸 2 本と周期長を返す。それ以外は use_periodic2=.false.。
  subroutine resolve_periodic2_axes(sim, use_periodic2, periodic_axes, periodic_len)
    type(sim_config), intent(in) :: sim
    logical, intent(out) :: use_periodic2
    integer(i32), intent(out) :: periodic_axes(2)
    real(dp), intent(out) :: periodic_len(2)

    character(len=16) :: field_bc_mode
    integer(i32) :: axis, n_periodic
    real(dp) :: span

    use_periodic2 = .false.
    periodic_axes = 0_i32
    periodic_len = 0.0d0
    field_bc_mode = lower_ascii(trim(sim%field_bc_mode))
    if (trim(field_bc_mode) /= 'periodic2') return

    if (.not. sim%use_box) then
      error stop 'fields.boundary="periodic2" requires a finite [domain].'
    end if

    n_periodic = 0_i32
    do axis = 1_i32, 3_i32
      if ((sim%bc_low(axis) == bc_periodic) .neqv. (sim%bc_high(axis) == bc_periodic)) then
        error stop 'periodic2 requires bc_low(axis)=bc_high(axis)=periodic for periodic axes.'
      end if
      if (sim%bc_low(axis) == bc_periodic) then
        n_periodic = n_periodic + 1_i32
        if (n_periodic <= 2_i32) periodic_axes(n_periodic) = axis
      end if
    end do
    if (n_periodic /= 2_i32) then
      error stop 'fields.boundary="periodic2" requires exactly two periodic axes.'
    end if

    do axis = 1_i32, 2_i32
      span = sim%box_max(periodic_axes(axis)) - sim%box_min(periodic_axes(axis))
      if (span <= 0.0d0) error stop 'periodic2 requires positive box length on periodic axes.'
      periodic_len(axis) = span
    end do

    use_periodic2 = .true.
  end subroutine resolve_periodic2_axes

end module bem_periodic2_axes

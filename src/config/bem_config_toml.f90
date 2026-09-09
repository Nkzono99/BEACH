!> TOML の基本型・固定長文字列・配列を検査して読み込む共通入口。
module bem_config_toml
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: error_unit
  use bem_kinds, only: dp, i32, i64
  use bem_types, only: bc_open, bc_reflect, bc_redistributed_reflect
  use bem_app_config_types, only: particle_bc_inherit, particle_inflow_none, particle_inflow_reservoir
  use bem_string_utils, only: lower_ascii
  use tomlf, only: toml_array, toml_key, toml_stat, toml_table, get_value, toml_len => len
  implicit none
  private
  public :: stop_config_error
  public :: require_toml_success
  public :: get_toml_real
  public :: get_toml_int
  public :: get_toml_logical
  public :: get_toml_string
  public :: get_toml_real3
  public :: get_toml_real2
  public :: get_toml_real4
  public :: get_toml_real_scalar_or_array3
  public :: get_toml_particle_boundary_mode
  public :: get_toml_boundary_inflow_mode

contains

  !> 診断を明示出力し、Intel の動的 ERROR STOP 文字列の問題を避ける。
  subroutine stop_config_error(message)
    character(len=*), intent(in) :: message

    write (error_unit, '(a)') trim(message)
    flush (error_unit)
    error stop 1
  end subroutine stop_config_error

  subroutine require_toml_success(stat, context)
    integer, intent(in) :: stat
    character(len=*), intent(in) :: context

    if (stat /= toml_stat%success) then
      call stop_config_error('Invalid TOML value for '//trim(context)//'.')
    end if
  end subroutine require_toml_success

  subroutine get_toml_real(table, key, value, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    real(dp), intent(out) :: value
    character(len=*), intent(in) :: context
    integer :: stat

    call get_value(table, key, value, stat=stat)
    call require_toml_success(stat, context)
    if (.not. ieee_is_finite(value)) call stop_config_error(trim(context)//' must be finite.')
  end subroutine get_toml_real

  subroutine get_toml_int(table, key, value, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    integer(i32), intent(out) :: value
    character(len=*), intent(in) :: context
    integer :: stat
    integer(i64) :: tmp

    call get_value(table, key, tmp, stat=stat)
    call require_toml_success(stat, context)
    if (tmp > int(huge(0_i32), i64) .or. tmp < -int(huge(0_i32), i64) - 1_i64) then
      call stop_config_error(trim(context)//' must fit a 32-bit signed integer.')
    end if
    value = int(tmp, i32)
  end subroutine get_toml_int

  subroutine get_toml_logical(table, key, value, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    logical, intent(out) :: value
    character(len=*), intent(in) :: context
    integer :: stat

    call get_value(table, key, value, stat=stat)
    call require_toml_success(stat, context)
  end subroutine get_toml_logical

  subroutine get_toml_string(table, key, value, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    character(len=*), intent(out) :: value
    character(len=*), intent(in) :: context
    character(len=:), allocatable :: tmp
    integer :: stat

    call get_value(table, key, tmp, stat=stat)
    call require_toml_success(stat, context)
    if (.not. allocated(tmp)) call stop_config_error('Invalid TOML value for '//trim(context)//'.')
    if (len_trim(tmp) > len(value)) call stop_config_error(trim(context)//' is too long.')
    value = ''
    value = trim(tmp)
  end subroutine get_toml_string

  subroutine get_toml_real3(table, key, value, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    real(dp), intent(out) :: value(3)
    character(len=*), intent(in) :: context
    type(toml_array), pointer :: array
    integer :: i, stat

    nullify (array)
    call get_value(table, key, array, stat=stat)
    call require_toml_success(stat, context)
    if (.not. associated(array)) call stop_config_error('Invalid TOML value for '//trim(context)//'.')
    if (toml_len(array) /= 3) then
      call stop_config_error(trim(context)//' must be an array of 3 numbers.')
    end if
    do i = 1, 3
      call get_value(array, i, value(i), stat=stat)
      call require_toml_success(stat, context)
      if (.not. ieee_is_finite(value(i))) call stop_config_error(trim(context)//' must contain finite values.')
    end do
  end subroutine get_toml_real3

  subroutine get_toml_real2(table, key, value, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    real(dp), intent(out) :: value(2)
    character(len=*), intent(in) :: context
    type(toml_array), pointer :: array
    integer :: i, stat

    nullify (array)
    call get_value(table, key, array, stat=stat)
    call require_toml_success(stat, context)
    if (.not. associated(array)) call stop_config_error('Invalid TOML value for '//trim(context)//'.')
    if (toml_len(array) /= 2) then
      call stop_config_error(trim(context)//' must be an array of 2 numbers.')
    end if
    do i = 1, 2
      call get_value(array, i, value(i), stat=stat)
      call require_toml_success(stat, context)
      if (.not. ieee_is_finite(value(i))) call stop_config_error(trim(context)//' must contain finite values.')
    end do
  end subroutine get_toml_real2

  subroutine get_toml_real4(table, key, value, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    real(dp), intent(out) :: value(4)
    character(len=*), intent(in) :: context
    type(toml_array), pointer :: array
    integer :: i, stat

    nullify (array)
    call get_value(table, key, array, stat=stat)
    call require_toml_success(stat, context)
    if (.not. associated(array)) call stop_config_error('Invalid TOML value for '//trim(context)//'.')
    if (toml_len(array) /= 4) then
      call stop_config_error(trim(context)//' must be an array of 4 numbers.')
    end if
    do i = 1, 4
      call get_value(array, i, value(i), stat=stat)
      call require_toml_success(stat, context)
      if (.not. ieee_is_finite(value(i))) call stop_config_error(trim(context)//' must contain finite values.')
    end do
  end subroutine get_toml_real4

  subroutine get_toml_real_scalar_or_array3(table, key, value, value_len, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    real(dp), intent(out) :: value(3)
    integer(i32), intent(out) :: value_len
    character(len=*), intent(in) :: context
    type(toml_array), pointer :: array
    real(dp) :: scalar_value
    integer :: i, n, stat

    value = 0.0d0
    value_len = 0_i32
    call get_value(table, key, scalar_value, stat=stat)
    if (stat == toml_stat%success) then
      if (.not. ieee_is_finite(scalar_value)) call stop_config_error(trim(context)//' must be finite.')
      value(1) = scalar_value
      value_len = 1_i32
      return
    end if

    nullify (array)
    call get_value(table, key, array, stat=stat)
    call require_toml_success(stat, context)
    if (.not. associated(array)) call stop_config_error('Invalid TOML value for '//trim(context)//'.')
    n = toml_len(array)
    if (n < 1 .or. n > 3) then
      call stop_config_error(trim(context)//' must be a number or an array of 1 to 3 numbers.')
    end if
    do i = 1, n
      call get_value(array, i, value(i), stat=stat)
      call require_toml_success(stat, context)
      if (.not. ieee_is_finite(value(i))) call stop_config_error(trim(context)//' must contain finite values.')
    end do
    value_len = int(n, i32)
  end subroutine get_toml_real_scalar_or_array3

  subroutine get_toml_particle_boundary_mode(table, key, value, context, allow_inherit)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    integer(i32), intent(out) :: value
    character(len=*), intent(in) :: context
    logical, intent(in) :: allow_inherit
    character(len=64) :: mode

    call get_toml_string(table, key, mode, context)
    select case (trim(lower_ascii(mode)))
    case ('inherit')
      if (.not. allow_inherit) then
        call stop_config_error(trim(context)//' must be "open", "reflect", or "redistributed_reflect".')
      end if
      value = particle_bc_inherit
    case ('open')
      value = bc_open
    case ('reflect')
      value = bc_reflect
    case ('redistributed_reflect')
      value = bc_redistributed_reflect
    case ('periodic')
      call stop_config_error(trim(context)//' cannot set periodicity; use domain.periodic_axes.')
    case default
      if (allow_inherit) then
        call stop_config_error(trim(context)//' must be "inherit", "open", "reflect", or "redistributed_reflect".')
      else
        call stop_config_error(trim(context)//' must be "open", "reflect", or "redistributed_reflect".')
      end if
    end select
  end subroutine get_toml_particle_boundary_mode

  subroutine get_toml_boundary_inflow_mode(table, key, value, context)
    type(toml_table), intent(inout) :: table
    type(toml_key), intent(in) :: key
    integer(i32), intent(out) :: value
    character(len=*), intent(in) :: context
    character(len=32) :: mode

    call get_toml_string(table, key, mode, context)
    select case (trim(lower_ascii(mode)))
    case ('reservoir')
      value = particle_inflow_reservoir
    case default
      call stop_config_error(trim(context)//' must be "reservoir". Omit the key to disable inflow.')
    end select
  end subroutine get_toml_boundary_inflow_mode

end module bem_config_toml

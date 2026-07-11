!> C ABI wrapper for the physical periodic zero mode.
module bem_periodic_zero_mode_c
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_f_pointer, c_int, c_loc, c_ptr
  use bem_kinds, only: dp, i32, i64
  use bem_periodic_zero_mode_eval, only: eval_periodic_zero_mode, zero_mode_trace_minus, &
                                         zero_mode_trace_plus, zero_mode_trace_principal_value
  use bem_periodic_zero_mode_plan, only: build_periodic_zero_mode_height_plan, periodic_zero_mode_ok, &
                                         periodic_zero_mode_plan_type, periodic_zero_mode_state_type, &
                                         refresh_periodic_zero_mode_state
  implicit none
  private

  integer(c_int), parameter, public :: beach_zero_mode_ok = 0_c_int
  integer(c_int), parameter, public :: beach_zero_mode_invalid_handle = 1_c_int
  integer(c_int), parameter, public :: beach_zero_mode_invalid_argument = 2_c_int
  integer(c_int), parameter, public :: beach_zero_mode_not_ready = 3_c_int

  type :: periodic_zero_mode_handle
    type(periodic_zero_mode_plan_type) :: plan
    type(periodic_zero_mode_state_type) :: state
    logical :: built = .false.
    logical :: charged = .false.
  end type periodic_zero_mode_handle

  public :: beach_zero_mode_create
  public :: beach_zero_mode_destroy
  public :: beach_zero_mode_build
  public :: beach_zero_mode_update
  public :: beach_zero_mode_eval

contains

  integer(c_int) function beach_zero_mode_create(handle) &
    bind(C, name='beach_zero_mode_create') result(status)
    type(c_ptr), intent(out) :: handle
    type(periodic_zero_mode_handle), pointer :: zero

    allocate (zero)
    zero%built = .false.
    zero%charged = .false.
    handle = c_loc(zero)
    status = beach_zero_mode_ok
  end function beach_zero_mode_create

  integer(c_int) function beach_zero_mode_destroy(handle) &
    bind(C, name='beach_zero_mode_destroy') result(status)
    type(c_ptr), value :: handle
    type(periodic_zero_mode_handle), pointer :: zero

    status = get_zero_mode(handle, zero)
    if (status /= beach_zero_mode_ok) return

    deallocate (zero)
    status = beach_zero_mode_ok
  end function beach_zero_mode_destroy

  integer(c_int) function beach_zero_mode_build(handle, nsrc, source_heights_ptr, area_xy) &
    bind(C, name='beach_zero_mode_build') result(status)
    type(c_ptr), value :: handle
    integer(c_int), value :: nsrc
    type(c_ptr), value :: source_heights_ptr
    real(c_double), value :: area_xy
    type(periodic_zero_mode_handle), pointer :: zero
    real(c_double), pointer :: source_heights(:, :)
    integer(i32) :: plan_status
    character(len=128) :: message

    status = get_zero_mode(handle, zero)
    if (status /= beach_zero_mode_ok) return
    if (.not. count_is_addressable(nsrc, 3_i64, 2_i64) .or. nsrc == 0_c_int .or. &
        .not. c_associated(source_heights_ptr) .or. .not. ieee_is_finite(area_xy) .or. &
        area_xy <= 0.0_c_double) then
      status = beach_zero_mode_invalid_argument
      return
    end if

    call c_f_pointer(source_heights_ptr, source_heights, [3, int(nsrc)])
    if (any(.not. ieee_is_finite(source_heights))) then
      status = beach_zero_mode_invalid_argument
      return
    end if

    call build_periodic_zero_mode_height_plan( &
      real(source_heights, dp), real(area_xy, dp), zero%plan, plan_status, message &
      )
    if (plan_status /= periodic_zero_mode_ok) then
      zero%built = .false.
      zero%charged = .false.
      status = beach_zero_mode_invalid_argument
      return
    end if

    zero%built = .true.
    zero%charged = .false.
    status = beach_zero_mode_ok
  end function beach_zero_mode_build

  integer(c_int) function beach_zero_mode_update( &
    handle, nsrc, charge_ptr, e_bottom, z_gauge, phi_gauge &
    ) bind(C, name='beach_zero_mode_update') result(status)
    type(c_ptr), value :: handle
    integer(c_int), value :: nsrc
    type(c_ptr), value :: charge_ptr
    real(c_double), value :: e_bottom, z_gauge, phi_gauge
    type(periodic_zero_mode_handle), pointer :: zero
    real(c_double), pointer :: charge(:)

    status = get_zero_mode(handle, zero)
    if (status /= beach_zero_mode_ok) return
    if (.not. zero%built) then
      status = beach_zero_mode_not_ready
      return
    end if
    if (.not. count_is_addressable(nsrc, 1_i64, 0_i64) .or. &
        int(nsrc, i64) /= int(zero%plan%nelem, i64) .or. .not. c_associated(charge_ptr) .or. &
        .not. ieee_is_finite(e_bottom) .or. .not. ieee_is_finite(z_gauge) .or. &
        .not. ieee_is_finite(phi_gauge)) then
      status = beach_zero_mode_invalid_argument
      return
    end if

    call c_f_pointer(charge_ptr, charge, [int(nsrc)])
    if (any(.not. ieee_is_finite(charge))) then
      status = beach_zero_mode_invalid_argument
      return
    end if

    call refresh_periodic_zero_mode_state( &
      zero%plan, real(charge, dp), real(e_bottom, dp), real(z_gauge, dp), real(phi_gauge, dp), zero%state &
      )
    zero%charged = .true.
    status = beach_zero_mode_ok
  end function beach_zero_mode_update

  integer(c_int) function beach_zero_mode_eval(handle, ntarget, z_ptr, trace, phi_ptr, ez_ptr) &
    bind(C, name='beach_zero_mode_eval') result(status)
    type(c_ptr), value :: handle
    integer(c_int), value :: ntarget
    type(c_ptr), value :: z_ptr
    integer(c_int), value :: trace
    type(c_ptr), value :: phi_ptr, ez_ptr
    type(periodic_zero_mode_handle), pointer :: zero
    real(c_double), pointer :: z(:), phi(:), ez(:)
    real(dp) :: phi_value, field_value
    integer(i32) :: i

    status = require_charged_zero_mode(handle, zero)
    if (status /= beach_zero_mode_ok) return
    if (.not. count_is_addressable(ntarget, 3_i64, 0_i64) .or. &
        .not. c_associated(z_ptr) .or. .not. c_associated(phi_ptr) .or. .not. c_associated(ez_ptr) .or. &
        .not. trace_is_valid(trace)) then
      status = beach_zero_mode_invalid_argument
      return
    end if
    if (ntarget == 0_c_int) return

    call c_f_pointer(z_ptr, z, [int(ntarget)])
    call c_f_pointer(phi_ptr, phi, [int(ntarget)])
    call c_f_pointer(ez_ptr, ez, [int(ntarget)])
    if (any(.not. ieee_is_finite(z))) then
      status = beach_zero_mode_invalid_argument
      return
    end if

    do i = 1_i32, int(ntarget, i32)
      call eval_periodic_zero_mode( &
        zero%plan, zero%state, real(z(i), dp), int(trace, i32), phi_value, field_value &
        )
      phi(i) = real(phi_value, c_double)
      ez(i) = real(field_value, c_double)
    end do
    status = beach_zero_mode_ok
  end function beach_zero_mode_eval

  integer(c_int) function get_zero_mode(handle, zero) result(status)
    type(c_ptr), intent(in) :: handle
    type(periodic_zero_mode_handle), pointer, intent(out) :: zero

    if (.not. c_associated(handle)) then
      nullify (zero)
      status = beach_zero_mode_invalid_handle
      return
    end if

    call c_f_pointer(handle, zero)
    status = beach_zero_mode_ok
  end function get_zero_mode

  integer(c_int) function require_charged_zero_mode(handle, zero) result(status)
    type(c_ptr), intent(in) :: handle
    type(periodic_zero_mode_handle), pointer, intent(out) :: zero

    status = get_zero_mode(handle, zero)
    if (status /= beach_zero_mode_ok) return
    if (.not. zero%built .or. .not. zero%charged) status = beach_zero_mode_not_ready
  end function require_charged_zero_mode

  pure logical function count_is_addressable(count, values_per_item, shape_margin) result(addressable)
    integer(c_int), intent(in) :: count
    integer(i64), intent(in) :: values_per_item, shape_margin
    integer(i64) :: shape_limit

    shape_limit = min(int(huge(0), i64), int(huge(0_i32), i64))
    addressable = .false.
    if (count < 0_c_int .or. values_per_item <= 0_i64) return
    if (shape_margin < 0_i64 .or. shape_margin > shape_limit) return
    addressable = int(count, i64) <= (shape_limit - shape_margin)/values_per_item
  end function count_is_addressable

  pure logical function trace_is_valid(trace) result(valid)
    integer(c_int), intent(in) :: trace

    valid = trace == int(zero_mode_trace_minus, c_int) .or. &
            trace == int(zero_mode_trace_principal_value, c_int) .or. &
            trace == int(zero_mode_trace_plus, c_int)
  end function trace_is_valid

end module bem_periodic_zero_mode_c

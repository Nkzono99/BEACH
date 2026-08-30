program test_periodic_zero_mode_c
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_negative_inf, ieee_positive_inf, ieee_quiet_nan, &
                                                                              ieee_value
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_intptr_t, c_loc, c_null_ptr, c_ptr
  use bem_kinds, only: dp, i32, i64
  use bem_periodic_zero_mode_eval, only: eval_periodic_zero_mode, zero_mode_trace_minus, zero_mode_trace_plus, &
                                         zero_mode_trace_principal_value
  use bem_periodic_zero_mode_plan, only: build_periodic_zero_mode_height_plan, periodic_zero_mode_ok, &
                                         periodic_zero_mode_plan_type, periodic_zero_mode_state_type, &
                                         refresh_periodic_zero_mode_state
  use test_support, only: assert_allclose_1d, assert_equal_i32, assert_true, test_begin, test_end, test_init, test_summary
  implicit none

  integer(c_int), parameter :: zero_ok = 0_c_int
  integer(c_int), parameter :: zero_invalid_handle = 1_c_int
  integer(c_int), parameter :: zero_invalid_argument = 2_c_int
  integer(c_int), parameter :: zero_not_ready = 3_c_int

  interface
    integer(c_int) function beach_zero_mode_create(handle) bind(C, name='beach_zero_mode_create') result(status)
      import :: c_int, c_ptr
      type(c_ptr), intent(out) :: handle
    end function beach_zero_mode_create

    integer(c_int) function beach_zero_mode_destroy(handle) bind(C, name='beach_zero_mode_destroy') result(status)
      import :: c_int, c_ptr
      type(c_ptr), value :: handle
    end function beach_zero_mode_destroy

    integer(c_int) function beach_zero_mode_build(handle, nsrc, source_heights_ptr, area_xy) &
      bind(C, name='beach_zero_mode_build') result(status)
      import :: c_double, c_int, c_ptr
      type(c_ptr), value :: handle
      integer(c_int), value :: nsrc
      type(c_ptr), value :: source_heights_ptr
      real(c_double), value :: area_xy
    end function beach_zero_mode_build

    integer(c_int) function beach_zero_mode_update( &
      handle, nsrc, charge_ptr, e_bottom, z_gauge, phi_gauge &
      ) bind(C, name='beach_zero_mode_update') result(status)
      import :: c_double, c_int, c_ptr
      type(c_ptr), value :: handle
      integer(c_int), value :: nsrc
      type(c_ptr), value :: charge_ptr
      real(c_double), value :: e_bottom, z_gauge, phi_gauge
    end function beach_zero_mode_update

    integer(c_int) function beach_zero_mode_eval(handle, ntarget, z_ptr, trace, phi_ptr, ez_ptr) &
      bind(C, name='beach_zero_mode_eval') result(status)
      import :: c_int, c_ptr
      type(c_ptr), value :: handle
      integer(c_int), value :: ntarget
      type(c_ptr), value :: z_ptr
      integer(c_int), value :: trace
      type(c_ptr), value :: phi_ptr, ez_ptr
    end function beach_zero_mode_eval
  end interface

  call test_init(6)

  call test_invalid_handles()
  call test_two_sheet_fixture()
  call test_inclined_fixture()
  call test_build_validation()
  call test_update_validation()
  call test_eval_validation()

  call test_summary()

contains

  subroutine test_invalid_handles()
    type(c_ptr) :: handle
    real(c_double), target :: source_heights(3, 1), charge(1), z(1), phi(1), ez(1)
    integer(c_int) :: status

    call test_begin('status_codes_and_invalid_handles')
    source_heights = 0.0_c_double
    charge = 1.0e-9_c_double
    z = 0.0_c_double
    phi = 0.0_c_double
    ez = 0.0_c_double

    status = beach_zero_mode_destroy(c_null_ptr)
    call assert_status(status, zero_invalid_handle, 'destroy NULL handle')
    status = beach_zero_mode_build(c_null_ptr, 1_c_int, c_loc(source_heights), 1.0_c_double)
    call assert_status(status, zero_invalid_handle, 'build NULL handle')
    status = beach_zero_mode_update( &
             c_null_ptr, 1_c_int, c_loc(charge), 0.0_c_double, 0.0_c_double, 0.0_c_double &
             )
    call assert_status(status, zero_invalid_handle, 'update NULL handle')
    status = beach_zero_mode_eval( &
             c_null_ptr, 1_c_int, c_loc(z), 0_c_int, c_loc(phi), c_loc(ez) &
             )
    call assert_status(status, zero_invalid_handle, 'eval NULL handle')

    status = beach_zero_mode_create(handle)
    call assert_status(status, zero_ok, 'create status')
    status = beach_zero_mode_destroy(handle)
    call assert_status(status, zero_ok, 'destroy status')
    call test_end()
  end subroutine test_invalid_handles

  subroutine test_two_sheet_fixture()
    type(c_ptr) :: handle
    type(periodic_zero_mode_plan_type) :: plan
    type(periodic_zero_mode_state_type) :: state
    real(c_double), target :: source_heights(3, 2), charge(2), z(5), phi(5), ez(5)
    real(dp) :: expected_phi(5), expected_ez(5)
    integer(c_int) :: status
    integer(i32) :: plan_status, i
    character(len=128) :: message

    call test_begin('two_sheet_nonneutral_matches_native')
    source_heights(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    source_heights(:, 2) = [1.0d0, 1.0d0, 1.0d0]
    charge = [2.0d-9, -0.5d-9]
    z = [-0.5d0, 0.0d0, 0.5d0, 1.0d0, 1.5d0]

    status = beach_zero_mode_create(handle)
    call assert_status(status, zero_ok, 'two-sheet create status')
    status = beach_zero_mode_build(handle, 2_c_int, c_loc(source_heights), 4.0d0)
    call assert_status(status, zero_ok, 'two-sheet build status')
    status = beach_zero_mode_update(handle, 2_c_int, c_loc(charge), 0.0d0, -0.5d0, 0.0d0)
    call assert_status(status, zero_ok, 'two-sheet update status')
    status = beach_zero_mode_eval(handle, 5_c_int, c_loc(z), 0_c_int, c_loc(phi), c_loc(ez))
    call assert_status(status, zero_ok, 'two-sheet eval status')
    call assert_true(all(ieee_is_finite(phi)) .and. all(ieee_is_finite(ez)), 'two-sheet outputs are finite')

    call build_periodic_zero_mode_height_plan(real(source_heights, dp), 4.0_dp, plan, plan_status, message)
    call assert_equal_i32(plan_status, periodic_zero_mode_ok, 'two-sheet native plan status')
    call refresh_periodic_zero_mode_state( &
      plan, real(charge, dp), 0.0_dp, -0.5_dp, 0.0_dp, state &
      )
    do i = 1_i32, 5_i32
      call eval_periodic_zero_mode( &
        plan, state, real(z(i), dp), zero_mode_trace_principal_value, expected_phi(i), expected_ez(i) &
        )
    end do
    call assert_true( &
      all(ieee_is_finite(expected_phi)) .and. all(ieee_is_finite(expected_ez)), 'two-sheet native outputs are finite' &
      )
    call assert_allclose_1d( &
      real(phi, dp), expected_phi, scaled_tolerance(expected_phi), 'two-sheet potential values' &
      )
    call assert_allclose_1d( &
      real(ez, dp), expected_ez, scaled_tolerance(expected_ez), 'two-sheet field values' &
      )

    status = beach_zero_mode_eval(handle, 5_c_int, c_loc(z), -1_c_int, c_loc(phi), c_loc(ez))
    call assert_status(status, zero_ok, 'two-sheet minus-trace eval status')
    call assert_true(all(ieee_is_finite(phi)) .and. all(ieee_is_finite(ez)), 'two-sheet minus-trace outputs are finite')
    do i = 1_i32, 5_i32
      call eval_periodic_zero_mode( &
        plan, state, real(z(i), dp), zero_mode_trace_minus, expected_phi(i), expected_ez(i) &
        )
    end do
    call assert_true( &
      all(ieee_is_finite(expected_phi)) .and. all(ieee_is_finite(expected_ez)), &
      'two-sheet minus-trace native outputs are finite' &
      )
    call assert_allclose_1d( &
      real(phi, dp), expected_phi, scaled_tolerance(expected_phi), 'two-sheet minus-trace potential values' &
      )
    call assert_allclose_1d( &
      real(ez, dp), expected_ez, scaled_tolerance(expected_ez), 'two-sheet minus-trace field values' &
      )

    status = beach_zero_mode_destroy(handle)
    call assert_status(status, zero_ok, 'two-sheet destroy status')
    call test_end()
  end subroutine test_two_sheet_fixture

  subroutine test_inclined_fixture()
    type(c_ptr) :: handle
    type(periodic_zero_mode_plan_type) :: plan
    type(periodic_zero_mode_state_type) :: state
    real(c_double), target :: source_heights(3, 1), charge(1), z(7), phi(7), ez(7)
    real(dp) :: expected_phi(7), expected_ez(7)
    integer(c_int) :: status
    integer(i32) :: plan_status, i
    character(len=128) :: message

    call test_begin('inclined_source_c_abi_matches_native')
    source_heights(:, 1) = [0.0d0, 1.0d0, 2.0d0]
    charge = [3.0d-9]
    z = [-0.5d0, 0.0d0, 0.5d0, 1.0d0, 1.5d0, 2.0d0, 2.5d0]

    status = beach_zero_mode_create(handle)
    call assert_status(status, zero_ok, 'inclined create status')
    status = beach_zero_mode_build(handle, 1_c_int, c_loc(source_heights), 2.5d0)
    call assert_status(status, zero_ok, 'inclined build status')
    status = beach_zero_mode_update(handle, 1_c_int, c_loc(charge), -0.25d0, 0.25d0, 1.5d0)
    call assert_status(status, zero_ok, 'inclined update status')
    status = beach_zero_mode_eval(handle, 7_c_int, c_loc(z), 1_c_int, c_loc(phi), c_loc(ez))
    call assert_status(status, zero_ok, 'inclined eval status')
    call assert_true(all(ieee_is_finite(phi)) .and. all(ieee_is_finite(ez)), 'inclined outputs are finite')

    call build_periodic_zero_mode_height_plan(real(source_heights, dp), 2.5_dp, plan, plan_status, message)
    call assert_equal_i32(plan_status, periodic_zero_mode_ok, 'inclined native plan status')
    call refresh_periodic_zero_mode_state( &
      plan, real(charge, dp), -0.25_dp, 0.25_dp, 1.5_dp, state &
      )
    do i = 1_i32, 7_i32
      call eval_periodic_zero_mode( &
        plan, state, real(z(i), dp), zero_mode_trace_plus, expected_phi(i), expected_ez(i) &
        )
    end do
    call assert_true( &
      all(ieee_is_finite(expected_phi)) .and. all(ieee_is_finite(expected_ez)), 'inclined native outputs are finite' &
      )
    call assert_allclose_1d( &
      real(phi, dp), expected_phi, scaled_tolerance(expected_phi), 'inclined potential values' &
      )
    call assert_allclose_1d( &
      real(ez, dp), expected_ez, scaled_tolerance(expected_ez), 'inclined field values' &
      )

    status = beach_zero_mode_destroy(handle)
    call assert_status(status, zero_ok, 'inclined destroy status')
    call test_end()
  end subroutine test_inclined_fixture

  subroutine test_build_validation()
    type(c_ptr) :: handle, non_dereferenceable_dummy
    real(c_double), target :: source_heights(3, 2), charge(2)
    real(c_double) :: invalid(3), area_values(5), saved
    integer(c_int) :: dangerous_count, status
    integer :: i

    call test_begin('build_rejects_invalid_inputs')
    call invalid_reals(invalid)
    area_values = [0.0_c_double, -1.0_c_double, invalid]
    source_heights(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    source_heights(:, 2) = [1.0d0, 1.0d0, 1.0d0]
    charge = [1.0d-9, -0.5d-9]

    status = beach_zero_mode_create(handle)
    call assert_status(status, zero_ok, 'build-validation create status')
    do i = 1, size(area_values)
      status = beach_zero_mode_build(handle, 2_c_int, c_loc(source_heights), area_values(i))
      call assert_status(status, zero_invalid_argument, 'invalid area status')
    end do
    status = beach_zero_mode_build(handle, 0_c_int, c_loc(source_heights), 1.0d0)
    call assert_status(status, zero_invalid_argument, 'zero source count')
    status = beach_zero_mode_build(handle, -1_c_int, c_loc(source_heights), 1.0d0)
    call assert_status(status, zero_invalid_argument, 'negative source count')
    dangerous_count = int(int(huge(0_i32), i64)/3_i64, c_int)
    non_dereferenceable_dummy = transfer(1_c_intptr_t, c_null_ptr)
    status = beach_zero_mode_build(handle, dangerous_count, non_dereferenceable_dummy, 1.0d0)
    call assert_status(status, zero_invalid_argument, 'source count leaves no state-shape margin')
    status = beach_zero_mode_build(handle, huge(0_c_int), c_loc(source_heights), 1.0d0)
    call assert_status(status, zero_invalid_argument, 'overflowing source count')
    status = beach_zero_mode_build(handle, 2_c_int, c_null_ptr, 1.0d0)
    call assert_status(status, zero_invalid_argument, 'NULL source heights')

    saved = source_heights(1, 1)
    do i = 1, size(invalid)
      source_heights(1, 1) = invalid(i)
      status = beach_zero_mode_build(handle, 2_c_int, c_loc(source_heights), 1.0d0)
      call assert_status(status, zero_invalid_argument, 'nonfinite source height')
    end do
    source_heights(1, 1) = saved

    status = beach_zero_mode_update(handle, 2_c_int, c_loc(charge), 0.0d0, 0.0d0, 0.0d0)
    call assert_status(status, zero_not_ready, 'update before build')
    status = beach_zero_mode_destroy(handle)
    call assert_status(status, zero_ok, 'build-validation destroy status')
    call test_end()
  end subroutine test_build_validation

  subroutine test_update_validation()
    type(c_ptr) :: handle
    real(c_double), target :: source_heights(3, 2), charge(2), z(1), phi(1), ez(1)
    real(c_double) :: invalid(3), saved
    integer(c_int) :: status
    integer :: i

    call test_begin('update_rejects_invalid_inputs_and_requires_build')
    call invalid_reals(invalid)
    source_heights(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    source_heights(:, 2) = [1.0d0, 1.0d0, 1.0d0]
    charge = [2.0d-9, -0.5d-9]
    z = 0.5d0
    phi = 0.0d0
    ez = 0.0d0

    status = beach_zero_mode_create(handle)
    call assert_status(status, zero_ok, 'update-validation create status')
    status = beach_zero_mode_build(handle, 2_c_int, c_loc(source_heights), 4.0d0)
    call assert_status(status, zero_ok, 'update-validation build status')
    status = beach_zero_mode_eval(handle, 1_c_int, c_loc(z), 0_c_int, c_loc(phi), c_loc(ez))
    call assert_status(status, zero_not_ready, 'eval before update')

    status = beach_zero_mode_update(handle, 1_c_int, c_loc(charge), 0.0d0, 0.0d0, 0.0d0)
    call assert_status(status, zero_invalid_argument, 'short charge count')
    status = beach_zero_mode_update(handle, 3_c_int, c_loc(charge), 0.0d0, 0.0d0, 0.0d0)
    call assert_status(status, zero_invalid_argument, 'long charge count')
    status = beach_zero_mode_update(handle, huge(0_c_int), c_loc(charge), 0.0d0, 0.0d0, 0.0d0)
    call assert_status(status, zero_invalid_argument, 'overflowing update count')
    status = beach_zero_mode_update(handle, 2_c_int, c_null_ptr, 0.0d0, 0.0d0, 0.0d0)
    call assert_status(status, zero_invalid_argument, 'NULL charges')

    saved = charge(1)
    do i = 1, size(invalid)
      charge(1) = invalid(i)
      status = beach_zero_mode_update(handle, 2_c_int, c_loc(charge), 0.0d0, 0.0d0, 0.0d0)
      call assert_status(status, zero_invalid_argument, 'nonfinite charge')
    end do
    charge(1) = saved
    do i = 1, size(invalid)
      status = beach_zero_mode_update(handle, 2_c_int, c_loc(charge), invalid(i), 0.0d0, 0.0d0)
      call assert_status(status, zero_invalid_argument, 'nonfinite e_bottom')
      status = beach_zero_mode_update(handle, 2_c_int, c_loc(charge), 0.0d0, invalid(i), 0.0d0)
      call assert_status(status, zero_invalid_argument, 'nonfinite z_gauge')
      status = beach_zero_mode_update(handle, 2_c_int, c_loc(charge), 0.0d0, 0.0d0, invalid(i))
      call assert_status(status, zero_invalid_argument, 'nonfinite phi_gauge')
    end do

    status = beach_zero_mode_update(handle, 2_c_int, c_loc(charge), 0.0d0, 0.0d0, 0.0d0)
    call assert_status(status, zero_ok, 'valid update after rejected inputs')
    status = beach_zero_mode_eval(handle, 1_c_int, c_loc(z), 0_c_int, c_loc(phi), c_loc(ez))
    call assert_status(status, zero_ok, 'valid eval after rejected updates')
    call assert_true(all(ieee_is_finite(phi)) .and. all(ieee_is_finite(ez)), 'valid outputs after rejected updates')
    status = beach_zero_mode_build(handle, 2_c_int, c_loc(source_heights), 4.0d0)
    call assert_status(status, zero_ok, 'rebuild status')
    status = beach_zero_mode_eval(handle, 1_c_int, c_loc(z), 0_c_int, c_loc(phi), c_loc(ez))
    call assert_status(status, zero_not_ready, 'rebuild clears charged state')

    status = beach_zero_mode_destroy(handle)
    call assert_status(status, zero_ok, 'update-validation destroy status')
    call test_end()
  end subroutine test_update_validation

  subroutine test_eval_validation()
    type(c_ptr) :: handle
    real(c_double), target :: source_heights(3, 1), charge(1), z(1), phi(1), ez(1)
    real(c_double) :: invalid(3), saved
    integer(c_int) :: status
    integer :: i

    call test_begin('eval_rejects_invalid_inputs')
    call invalid_reals(invalid)
    source_heights = 0.0d0
    charge = 1.0d-9
    z = 0.0d0
    phi = 0.0d0
    ez = 0.0d0

    status = beach_zero_mode_create(handle)
    call assert_status(status, zero_ok, 'eval-validation create status')
    status = beach_zero_mode_build(handle, 1_c_int, c_loc(source_heights), 1.0d0)
    call assert_status(status, zero_ok, 'eval-validation build status')
    status = beach_zero_mode_update(handle, 1_c_int, c_loc(charge), 0.0d0, 0.0d0, 0.0d0)
    call assert_status(status, zero_ok, 'eval-validation update status')

    status = beach_zero_mode_eval(handle, 1_c_int, c_loc(z), -2_c_int, c_loc(phi), c_loc(ez))
    call assert_status(status, zero_invalid_argument, 'trace below range')
    status = beach_zero_mode_eval(handle, 1_c_int, c_loc(z), 2_c_int, c_loc(phi), c_loc(ez))
    call assert_status(status, zero_invalid_argument, 'trace above range')
    status = beach_zero_mode_eval(handle, -1_c_int, c_loc(z), 0_c_int, c_loc(phi), c_loc(ez))
    call assert_status(status, zero_invalid_argument, 'negative target count')
    status = beach_zero_mode_eval(handle, huge(0_c_int), c_loc(z), 0_c_int, c_loc(phi), c_loc(ez))
    call assert_status(status, zero_invalid_argument, 'overflowing target count')
    status = beach_zero_mode_eval(handle, 1_c_int, c_null_ptr, 0_c_int, c_loc(phi), c_loc(ez))
    call assert_status(status, zero_invalid_argument, 'NULL target heights')
    status = beach_zero_mode_eval(handle, 1_c_int, c_loc(z), 0_c_int, c_null_ptr, c_loc(ez))
    call assert_status(status, zero_invalid_argument, 'NULL potential output')
    status = beach_zero_mode_eval(handle, 1_c_int, c_loc(z), 0_c_int, c_loc(phi), c_null_ptr)
    call assert_status(status, zero_invalid_argument, 'NULL field output')

    saved = z(1)
    do i = 1, size(invalid)
      z(1) = invalid(i)
      status = beach_zero_mode_eval(handle, 1_c_int, c_loc(z), 0_c_int, c_loc(phi), c_loc(ez))
      call assert_status(status, zero_invalid_argument, 'nonfinite target height')
    end do
    z(1) = saved
    status = beach_zero_mode_eval(handle, 0_c_int, c_loc(z), 0_c_int, c_loc(phi), c_loc(ez))
    call assert_status(status, zero_ok, 'empty eval')

    status = beach_zero_mode_destroy(handle)
    call assert_status(status, zero_ok, 'eval-validation destroy status')
    call test_end()
  end subroutine test_eval_validation

  subroutine invalid_reals(values)
    real(c_double), intent(out) :: values(3)

    values = [ &
             ieee_value(0.0_c_double, ieee_quiet_nan), &
             ieee_value(0.0_c_double, ieee_positive_inf), &
             ieee_value(0.0_c_double, ieee_negative_inf) &
             ]
  end subroutine invalid_reals

  real(dp) function scaled_tolerance(expected) result(tolerance)
    real(dp), intent(in) :: expected(:)

    tolerance = 1.0e-12_dp*max(1.0_dp, maxval(abs(expected)))
  end function scaled_tolerance

  subroutine assert_status(actual, expected, message)
    integer(c_int), intent(in) :: actual, expected
    character(len=*), intent(in) :: message

    call assert_equal_i32(int(actual, i32), int(expected, i32), message)
  end subroutine assert_status

end program test_periodic_zero_mode_c

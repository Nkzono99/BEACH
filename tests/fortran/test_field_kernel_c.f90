program test_field_kernel_c
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_loc, c_null_ptr, c_ptr
  use bem_constants, only: k_coulomb
  use bem_field_kernel_c, only: beach_kernel_build, beach_kernel_create, beach_kernel_destroy, beach_kernel_eval_e, &
                                beach_kernel_force_on_charges, beach_kernel_ok, beach_kernel_update_charges
  use test_support, only: assert_allclose_1d, assert_equal_i32, test_begin, test_summary
  implicit none

  type(c_ptr) :: handle
  integer(c_int) :: status
  real(c_double), target :: src_pos(3, 2), src_q(2), target_pos(3, 1), e(3, 1)
  real(c_double), target :: target_q(1), origin(3), force(3), torque(3)
  real(c_double) :: expected_e(3), expected_force(3), expected_torque(3)

  call test_begin('field_kernel_c_free_eval_and_force')

  src_pos(:, 1) = [0.0d0, 0.0d0, 0.0d0]
  src_pos(:, 2) = [1.0d0, 0.0d0, 0.0d0]
  src_q = [1.0d-9, -2.0d-9]
  target_pos(:, 1) = [0.0d0, 1.0d0, 0.0d0]
  target_q = [3.0d-9]
  origin = [0.0d0, 0.0d0, 0.0d0]

  status = beach_kernel_create(handle)
  call assert_equal_i32(status, beach_kernel_ok, 'create status')

  status = beach_kernel_build( &
           handle, 2_c_int, c_loc(src_pos), 0.5d0, 8_c_int, 4_c_int, 0.0d0, 0_c_int, c_null_ptr, &
           c_null_ptr, 1_c_int, 0_c_int, 0.0d0, 4_c_int, c_null_ptr, c_null_ptr &
           )
  call assert_equal_i32(status, beach_kernel_ok, 'build status')

  status = beach_kernel_update_charges(handle, 2_c_int, c_loc(src_q))
  call assert_equal_i32(status, beach_kernel_ok, 'update status')

  e = 0.0d0
  status = beach_kernel_eval_e(handle, 1_c_int, c_loc(target_pos), c_loc(e))
  call assert_equal_i32(status, beach_kernel_ok, 'eval_e status')

  expected_e = direct_e(src_pos, src_q, target_pos(:, 1))
  call assert_allclose_1d(e(:, 1), expected_e, 1.0d-15, 'eval_e value')

  force = 0.0d0
  torque = 0.0d0
  status = beach_kernel_force_on_charges( &
           handle, 1_c_int, c_loc(target_pos), c_loc(target_q), c_loc(origin), c_loc(force), c_loc(torque) &
           )
  call assert_equal_i32(status, beach_kernel_ok, 'force status')
  expected_force = target_q(1)*expected_e
  expected_torque = cross(target_pos(:, 1) - origin, expected_force)
  call assert_allclose_1d(force, expected_force, 1.0d-20, 'force value')
  call assert_allclose_1d(torque, expected_torque, 1.0d-20, 'torque value')

  status = beach_kernel_destroy(handle)
  call assert_equal_i32(status, beach_kernel_ok, 'destroy status')

  call test_summary()

contains

  function direct_e(src, q, target) result(e_out)
    real(c_double), intent(in) :: src(:, :), q(:), target(3)
    real(c_double) :: e_out(3)
    integer :: i
    real(c_double) :: d(3), r2, inv_r3

    e_out = 0.0d0
    do i = 1, size(q)
      d = target - src(:, i)
      r2 = sum(d*d)
      inv_r3 = 1.0d0/(sqrt(r2)*r2)
      e_out = e_out + k_coulomb*q(i)*inv_r3*d
    end do
  end function direct_e

  pure function cross(a, b) result(c)
    real(c_double), intent(in) :: a(3), b(3)
    real(c_double) :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross

end program test_field_kernel_c

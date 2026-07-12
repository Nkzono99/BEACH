program test_regularized_qr
  use bem_kinds, only: dp, i32
  use bem_regularized_qr, only: regularized_qr_type, prepare_regularized_qr, solve_regularized_qr
  use test_support, only: test_begin, test_end, test_summary, assert_equal_i32, assert_allclose_1d
  implicit none

  type(regularized_qr_type) :: factor
  real(dp), parameter :: ridge = 1.0e-12_dp
  real(dp) :: matrix(6, 3), rhs(6), actual(3), expected(3)
  integer(i32) :: rhs_idx

  matrix = reshape([ &
                   1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp, 3.0_dp, -2.0_dp, &
                   0.0_dp, 1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp, 1.5_dp, &
                   2.0_dp, -0.5_dp, 1.0_dp, 3.0_dp, -1.0_dp, 0.25_dp &
                   ], shape(matrix))

  call test_begin('one_factorization_solves_multiple_rhs')
  call prepare_regularized_qr(factor, matrix, ridge, 1.0e-12_dp)
  do rhs_idx = 1_i32, 3_i32
    rhs = [ &
          real(rhs_idx, dp), -0.5_dp*real(rhs_idx, dp), 2.0_dp - real(rhs_idx, dp), &
          0.25_dp, -1.5_dp, 0.75_dp*real(rhs_idx, dp) &
          ]
    call solve_regularized_qr(factor, rhs, actual)
    call solve_reference(matrix, rhs, ridge, expected)
    call assert_allclose_1d(actual, expected, 5.0e-11_dp, 'reused QR solution differs from reference')
  end do
  call assert_equal_i32(factor%preparation_count, 1_i32, 'factorization must be prepared once')
  call test_end()
  call test_summary()

contains

  subroutine solve_reference(a, b, ridge_value, x)
    real(dp), intent(in) :: a(:, :), b(:), ridge_value
    real(dp), intent(out) :: x(:)
    real(dp) :: normal(size(a, 2), size(a, 2)), rhs_normal(size(a, 2))
    real(dp) :: scale(size(a, 2))
    integer :: col

    do col = 1, size(a, 2)
      scale(col) = max(sqrt(sum(a(:, col)*a(:, col))), 1.0_dp)
    end do
    normal = matmul(transpose(a), a)
    do col = 1, size(a, 2)
      normal(col, col) = normal(col, col) + ridge_value*scale(col)*scale(col)
    end do
    rhs_normal = matmul(transpose(a), b)
    call solve_square(normal, rhs_normal, x)
  end subroutine solve_reference

  subroutine solve_square(a, b, x)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp), intent(out) :: x(:)
    real(dp) :: work(size(a, 1), size(a, 2)), rhs(size(b)), pivot, factor
    real(dp) :: row(size(a, 2)), scalar
    integer :: col, row_idx, pivot_idx

    work = a
    rhs = b
    do col = 1, size(a, 2)
      pivot_idx = col - 1 + maxloc(abs(work(col:, col)), dim=1)
      row = work(col, :)
      work(col, :) = work(pivot_idx, :)
      work(pivot_idx, :) = row
      scalar = rhs(col)
      rhs(col) = rhs(pivot_idx)
      rhs(pivot_idx) = scalar
      pivot = work(col, col)
      do row_idx = col + 1, size(a, 1)
        factor = work(row_idx, col)/pivot
        work(row_idx, col:) = work(row_idx, col:) - factor*work(col, col:)
        rhs(row_idx) = rhs(row_idx) - factor*rhs(col)
      end do
    end do
    x = rhs
    do row_idx = size(a, 1), 1, -1
      if (row_idx < size(a, 1)) then
        x(row_idx) = x(row_idx) - sum(work(row_idx, row_idx + 1:)*x(row_idx + 1:))
      end if
      x(row_idx) = x(row_idx)/work(row_idx, row_idx)
    end do
  end subroutine solve_square
end program test_regularized_qr

program test_regularized_qr
  use bem_kinds, only: dp, i32
  use bem_regularized_qr, only: regularized_qr_type, prepare_regularized_qr, solve_regularized_qr
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_allclose_1d
  implicit none

  real(dp), parameter :: full_rank_matrix(6, 3) = reshape([ &
                                                          1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp, 3.0_dp, -2.0_dp, &
                                                          0.0_dp, 1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp, 1.5_dp, &
                                                          2.0_dp, -0.5_dp, 1.0_dp, 3.0_dp, -1.0_dp, 0.25_dp &
                                                          ], [6, 3])
  real(dp), parameter :: known_solutions(3, 3) = reshape([ &
                                                         1.0_dp, -2.0_dp, 0.5_dp, &
                                                         -0.25_dp, 0.75_dp, 2.0_dp, &
                                                         3.0_dp, 0.0_dp, -1.0_dp &
                                                         ], [3, 3])
  real(dp), parameter :: dependent_matrix(4, 3) = reshape([ &
                                                          1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                                                          2.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, &
                                                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp &
                                                          ], [4, 3])
  real(dp), parameter :: dependent_rhs(4) = [3.0_dp, 0.0_dp, 3.0_dp, 0.0_dp]
  real(dp), parameter :: dependent_expected(3) = [0.75_dp, 0.375_dp, 0.0_dp]
  real(dp), parameter :: solve_tolerance = 2.0e-12_dp
  type(regularized_qr_type) :: factor
  real(dp) :: full_rank_rhs(6), actual(3), expected(3)
  integer(i32) :: rhs_idx

  call test_init(2)

  call test_begin('full_rank_factorization_reuses_exact_solutions')
  call prepare_regularized_qr(factor, full_rank_matrix, 0.0_dp, 1.0e-12_dp)
  do rhs_idx = 1_i32, 3_i32
    expected = known_solutions(:, rhs_idx)
    full_rank_rhs = matmul(full_rank_matrix, expected)
    call solve_regularized_qr(factor, full_rank_rhs, actual)
    call assert_allclose_1d(actual, expected, solve_tolerance, 'reused QR solution differs from exact solution')
  end do
  call test_end()

  call test_begin('ridge_stabilizes_dependent_and_zero_columns')
  call prepare_regularized_qr(factor, dependent_matrix, 2.0_dp, 1.0e-12_dp)
  call solve_regularized_qr(factor, dependent_rhs, actual)
  call assert_allclose_1d(actual, dependent_expected, solve_tolerance, 'regularized dependent-column solution mismatch')
  call test_end()

  call test_summary()
end program test_regularized_qr

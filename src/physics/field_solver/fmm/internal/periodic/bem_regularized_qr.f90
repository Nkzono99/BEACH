!> Reusable column-scaled QR factorization for ridge-regularized least squares.
module bem_regularized_qr
  use bem_kinds, only: dp, i32
  implicit none
  private

  type, public :: regularized_qr_type
    real(dp), allocatable :: q(:, :)
    real(dp), allocatable :: r(:, :)
    real(dp), allocatable :: col_scale(:)
    integer(i32) :: mrow = 0_i32
    integer(i32) :: ncol = 0_i32
    integer(i32) :: preparation_count = 0_i32
  end type regularized_qr_type

  public :: prepare_regularized_qr
  public :: solve_regularized_qr

contains

  subroutine prepare_regularized_qr(factorization, matrix, ridge, qr_tolerance)
    type(regularized_qr_type), intent(inout) :: factorization
    real(dp), intent(in) :: matrix(:, :), ridge, qr_tolerance
    real(dp), allocatable :: augmented(:, :)
    real(dp) :: ridge_sqrt
    integer(i32) :: mrow, ncol, col_idx

    mrow = int(size(matrix, 1), i32)
    ncol = int(size(matrix, 2), i32)
    if (mrow <= 0_i32 .or. ncol <= 0_i32 .or. ridge < 0.0_dp .or. qr_tolerance <= 0.0_dp) then
      error stop 'prepare_regularized_qr received invalid dimensions or tolerances.'
    end if

    if (allocated(factorization%q)) deallocate (factorization%q)
    if (allocated(factorization%r)) deallocate (factorization%r)
    if (allocated(factorization%col_scale)) deallocate (factorization%col_scale)
    allocate (factorization%q(mrow + ncol, ncol), factorization%r(ncol, ncol))
    allocate (factorization%col_scale(ncol), augmented(mrow + ncol, ncol))

    augmented = 0.0_dp
    ridge_sqrt = sqrt(ridge)
    do col_idx = 1_i32, ncol
      factorization%col_scale(col_idx) = sqrt(sum(matrix(:, col_idx)*matrix(:, col_idx)))
      if (factorization%col_scale(col_idx) <= tiny(1.0_dp)) factorization%col_scale(col_idx) = 1.0_dp
      augmented(1:mrow, col_idx) = matrix(:, col_idx)/factorization%col_scale(col_idx)
      augmented(mrow + col_idx, col_idx) = ridge_sqrt
    end do
    call factor_tall_matrix_qr(augmented, factorization%q, factorization%r, qr_tolerance)
    factorization%mrow = mrow
    factorization%ncol = ncol
    factorization%preparation_count = factorization%preparation_count + 1_i32
  end subroutine prepare_regularized_qr

  subroutine solve_regularized_qr(factorization, rhs, solution)
    type(regularized_qr_type), intent(in) :: factorization
    real(dp), intent(in) :: rhs(:)
    real(dp), intent(out) :: solution(:)
    real(dp), allocatable :: augmented_rhs(:), qtb(:), scaled_solution(:)
    integer(i32) :: col_idx

    if (.not. allocated(factorization%q) .or. .not. allocated(factorization%r) .or. &
        .not. allocated(factorization%col_scale)) then
      error stop 'solve_regularized_qr requires a prepared factorization.'
    end if
    if (size(rhs) /= factorization%mrow .or. size(solution) /= factorization%ncol) then
      error stop 'solve_regularized_qr dimension mismatch.'
    end if

    allocate (augmented_rhs(factorization%mrow + factorization%ncol))
    allocate (qtb(factorization%ncol), scaled_solution(factorization%ncol))
    augmented_rhs = 0.0_dp
    augmented_rhs(1:factorization%mrow) = rhs
    qtb = matmul(transpose(factorization%q), augmented_rhs)
    call solve_upper_triangular_system(factorization%r, qtb, scaled_solution)
    do col_idx = 1_i32, factorization%ncol
      solution(col_idx) = scaled_solution(col_idx)/factorization%col_scale(col_idx)
    end do
  end subroutine solve_regularized_qr

  subroutine factor_tall_matrix_qr(matrix, q, r, qr_tolerance)
    real(dp), intent(in) :: matrix(:, :), qr_tolerance
    real(dp), intent(out) :: q(:, :), r(:, :)
    integer(i32) :: mrow, ncol, col_idx, basis_idx
    real(dp), allocatable :: v(:)
    real(dp) :: norm_v, corr, base_norm

    mrow = int(size(matrix, 1), i32)
    ncol = int(size(matrix, 2), i32)
    if (size(q, 1) /= mrow .or. size(q, 2) /= ncol) error stop 'factor_tall_matrix_qr q dimension mismatch.'
    if (size(r, 1) /= ncol .or. size(r, 2) /= ncol) error stop 'factor_tall_matrix_qr r dimension mismatch.'

    q = 0.0_dp
    r = 0.0_dp
    allocate (v(mrow))
    do col_idx = 1_i32, ncol
      v = matrix(:, col_idx)
      base_norm = max(sqrt(sum(v*v)), 1.0_dp)
      do basis_idx = 1_i32, col_idx - 1_i32
        r(basis_idx, col_idx) = dot_product(q(:, basis_idx), v)
        v = v - r(basis_idx, col_idx)*q(:, basis_idx)
      end do
      do basis_idx = 1_i32, col_idx - 1_i32
        corr = dot_product(q(:, basis_idx), v)
        r(basis_idx, col_idx) = r(basis_idx, col_idx) + corr
        v = v - corr*q(:, basis_idx)
      end do
      norm_v = sqrt(sum(v*v))
      if (norm_v <= qr_tolerance*base_norm) then
        r(col_idx, col_idx) = qr_tolerance*base_norm
      else
        r(col_idx, col_idx) = norm_v
        q(:, col_idx) = v/norm_v
      end if
    end do
  end subroutine factor_tall_matrix_qr

  subroutine solve_upper_triangular_system(matrix, rhs, solution)
    real(dp), intent(in) :: matrix(:, :), rhs(:)
    real(dp), intent(out) :: solution(:)
    integer(i32) :: ncol, row_idx, col_idx
    real(dp) :: diag_val

    ncol = int(size(matrix, 1), i32)
    if (size(matrix, 2) /= ncol .or. size(rhs) /= ncol .or. size(solution) /= ncol) then
      error stop 'solve_upper_triangular_system dimension mismatch.'
    end if

    solution = rhs
    do row_idx = ncol, 1_i32, -1_i32
      do col_idx = row_idx + 1_i32, ncol
        solution(row_idx) = solution(row_idx) - matrix(row_idx, col_idx)*solution(col_idx)
      end do
      diag_val = matrix(row_idx, row_idx)
      if (abs(diag_val) <= tiny(1.0_dp)) diag_val = sign(tiny(1.0_dp), diag_val + tiny(1.0_dp))
      solution(row_idx) = solution(row_idx)/diag_val
    end do
  end subroutine solve_upper_triangular_system
end module bem_regularized_qr

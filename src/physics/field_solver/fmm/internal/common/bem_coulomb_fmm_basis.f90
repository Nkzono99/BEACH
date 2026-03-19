!> Coulomb FMM の multi-index と微分テーブル計算。
module bem_coulomb_fmm_basis
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_types, only: fmm_plan_type
  implicit none
  private

  public :: initialize_basis_tables
  public :: build_axis_powers
  public :: compute_laplace_derivatives

contains

  subroutine initialize_basis_tables(plan, order)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32), intent(in) :: order
    integer(i32) :: idx, alpha_idx, beta_idx, ax, term_idx
    integer(i32) :: sum_exp(3)

    plan%ncoef = count_multi_indices(order)
    plan%nderiv = count_multi_indices(2_i32*order)
    call build_multi_index_table(order, plan%alpha, plan%alpha_degree, plan%alpha_factorial, plan%alpha_map)
    call build_multi_index_table(2_i32*order, plan%deriv_alpha, plan%deriv_degree, plan%deriv_factorial, plan%deriv_map)
    allocate (plan%alpha_sign(plan%ncoef), plan%alpha_plus_axis(3, plan%ncoef), &
              plan%alpha_beta_deriv_idx(plan%ncoef, plan%ncoef))
    do idx = 1_i32, plan%ncoef
      if (mod(plan%alpha_degree(idx), 2_i32) == 0_i32) then
        plan%alpha_sign(idx) = 1.0d0
      else
        plan%alpha_sign(idx) = -1.0d0
      end if
      do ax = 1_i32, 3_i32
        sum_exp = plan%alpha(:, idx)
        sum_exp(ax) = sum_exp(ax) + 1_i32
        if (sum(sum_exp) <= order) then
          plan%alpha_plus_axis(ax, idx) = plan%alpha_map(sum_exp(1), sum_exp(2), sum_exp(3))
        else
          plan%alpha_plus_axis(ax, idx) = 0_i32
        end if
      end do
    end do

    do alpha_idx = 1_i32, plan%ncoef
      do beta_idx = 1_i32, plan%ncoef
        sum_exp = plan%alpha(:, alpha_idx) + plan%alpha(:, beta_idx)
        plan%alpha_beta_deriv_idx(alpha_idx, beta_idx) = plan%deriv_map(sum_exp(1), sum_exp(2), sum_exp(3))
      end do
    end do

    plan%eval_term_count = count(plan%alpha_degree < order)
    allocate (plan%eval_exp(3, max(1_i32, plan%eval_term_count)))
    allocate (plan%eval_deriv_idx(3, max(1_i32, plan%eval_term_count)))
    allocate (plan%eval_inv_factorial(max(1_i32, plan%eval_term_count)))
    plan%eval_exp = 0_i32
    plan%eval_deriv_idx = 0_i32
    plan%eval_inv_factorial = 0.0d0

    term_idx = 0_i32
    do alpha_idx = 1_i32, plan%ncoef
      if (plan%alpha_degree(alpha_idx) >= order) cycle
      term_idx = term_idx + 1_i32
      plan%eval_exp(:, term_idx) = plan%alpha(:, alpha_idx)
      plan%eval_deriv_idx(:, term_idx) = plan%alpha_plus_axis(:, alpha_idx)
      plan%eval_inv_factorial(term_idx) = 1.0d0/plan%alpha_factorial(alpha_idx)
    end do
  end subroutine initialize_basis_tables

  integer(i32) function count_multi_indices(order)
    integer(i32), intent(in) :: order

    count_multi_indices = (order + 1_i32)*(order + 2_i32)*(order + 3_i32)/6_i32
  end function count_multi_indices

  subroutine build_multi_index_table(order, alpha, degree, factorials, map)
    integer(i32), intent(in) :: order
    integer(i32), allocatable, intent(out) :: alpha(:, :)
    integer(i32), allocatable, intent(out) :: degree(:)
    real(dp), allocatable, intent(out) :: factorials(:)
    integer(i32), allocatable, intent(out) :: map(:, :, :)
    integer(i32) :: total, ax, ay, az, idx

    allocate (alpha(3, count_multi_indices(order)), degree(count_multi_indices(order)))
    allocate (factorials(count_multi_indices(order)))
    allocate (map(0:order, 0:order, 0:order))
    map = 0_i32
    idx = 0_i32
    do total = 0_i32, order
      do ax = 0_i32, total
        do ay = 0_i32, total - ax
          az = total - ax - ay
          idx = idx + 1_i32
          alpha(:, idx) = [ax, ay, az]
          degree(idx) = total
          factorials(idx) = factorial_value(ax)*factorial_value(ay)*factorial_value(az)
          map(ax, ay, az) = idx
        end do
      end do
    end do
  end subroutine build_multi_index_table

  real(dp) function factorial_value(n)
    integer(i32), intent(in) :: n
    integer(i32) :: i

    factorial_value = 1.0d0
    do i = 2_i32, n
      factorial_value = factorial_value*real(i, dp)
    end do
  end function factorial_value

  subroutine build_axis_powers(d, order, xpow, ypow, zpow)
    real(dp), intent(in) :: d(3)
    integer(i32), intent(in) :: order
    real(dp), intent(out) :: xpow(0:order), ypow(0:order), zpow(0:order)
    integer(i32) :: k

    xpow(0) = 1.0d0
    ypow(0) = 1.0d0
    zpow(0) = 1.0d0
    do k = 1_i32, order
      xpow(k) = xpow(k - 1_i32)*d(1)
      ypow(k) = ypow(k - 1_i32)*d(2)
      zpow(k) = zpow(k - 1_i32)*d(3)
    end do
  end subroutine build_axis_powers

  subroutine compute_laplace_derivatives(plan, r, deriv)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: r(3)
    real(dp), intent(out) :: deriv(:)
    real(dp) :: a, coeff, soft2
    real(dp) :: q(plan%nderiv), power(plan%nderiv), accum(plan%nderiv), tmp(plan%nderiv)
    integer(i32) :: idx0, idx, n

    if (size(deriv) /= plan%nderiv) error stop 'Derivative buffer size mismatch.'
    soft2 = plan%options%softening*plan%options%softening
    a = dot_product(r, r) + soft2
    if (a <= tiny(1.0d0)) error stop 'Cannot expand Laplace kernel at r=0.'

    q = 0.0d0
    power = 0.0d0
    accum = 0.0d0
    idx0 = plan%deriv_map(0_i32, 0_i32, 0_i32)
    power(idx0) = 1.0d0
    accum(idx0) = 1.0d0
    if (plan%options%order >= 1_i32) then
      q(plan%deriv_map(1_i32, 0_i32, 0_i32)) = 2.0d0*r(1)/a
      q(plan%deriv_map(0_i32, 1_i32, 0_i32)) = 2.0d0*r(2)/a
      q(plan%deriv_map(0_i32, 0_i32, 1_i32)) = 2.0d0*r(3)/a
    end if
    if (2_i32 <= 2_i32*plan%options%order) then
      q(plan%deriv_map(2_i32, 0_i32, 0_i32)) = q(plan%deriv_map(2_i32, 0_i32, 0_i32)) + 1.0d0/a
      q(plan%deriv_map(0_i32, 2_i32, 0_i32)) = q(plan%deriv_map(0_i32, 2_i32, 0_i32)) + 1.0d0/a
      q(plan%deriv_map(0_i32, 0_i32, 2_i32)) = q(plan%deriv_map(0_i32, 0_i32, 2_i32)) + 1.0d0/a
    end if

    coeff = 1.0d0
    do n = 1_i32, 2_i32*plan%options%order
      call multiply_polynomials(plan, power, q, tmp)
      power = tmp
      coeff = coeff*(-(2.0d0*real(n, dp) - 1.0d0)/(2.0d0*real(n, dp)))
      accum = accum + coeff*power
    end do
    accum = accum/sqrt(a)
    do idx = 1_i32, plan%nderiv
      deriv(idx) = accum(idx)*plan%deriv_factorial(idx)
    end do
  end subroutine compute_laplace_derivatives

  subroutine multiply_polynomials(plan, a_poly, b_poly, out_poly)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: a_poly(:), b_poly(:)
    real(dp), intent(out) :: out_poly(:)
    integer(i32) :: ia, ib, idx
    integer(i32) :: sum_exp(3)

    if (size(a_poly) /= plan%nderiv .or. size(b_poly) /= plan%nderiv .or. size(out_poly) /= plan%nderiv) then
      error stop 'Polynomial buffer size mismatch.'
    end if

    out_poly = 0.0d0
    do ia = 1_i32, plan%nderiv
      if (abs(a_poly(ia)) <= tiny(1.0d0)) cycle
      do ib = 1_i32, plan%nderiv
        if (abs(b_poly(ib)) <= tiny(1.0d0)) cycle
        sum_exp = plan%deriv_alpha(:, ia) + plan%deriv_alpha(:, ib)
        if (sum(sum_exp) > 2_i32*plan%options%order) cycle
        idx = plan%deriv_map(sum_exp(1), sum_exp(2), sum_exp(3))
        out_poly(idx) = out_poly(idx) + a_poly(ia)*b_poly(ib)
      end do
    end do
  end subroutine multiply_polynomials

end module bem_coulomb_fmm_basis

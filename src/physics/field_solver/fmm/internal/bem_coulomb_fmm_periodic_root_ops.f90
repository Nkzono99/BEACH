!> periodic2 root operator の前計算。
module bem_coulomb_fmm_periodic_root_ops
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_types, only: fmm_plan_type
  use bem_coulomb_fmm_basis, only: build_axis_powers, compute_laplace_derivatives
  use bem_coulomb_fmm_periodic, only: use_periodic2_m2l_root_trunc, use_periodic2_m2l_root_oracle
  use bem_coulomb_fmm_periodic_ewald, only: add_periodic2_exact_ewald_correction_single_source
  use bem_coulomb_fmm_tree_utils, only: active_tree_node_center, active_tree_node_half_size
  implicit none
  private

  public :: precompute_periodic_root_operator

contains

  subroutine precompute_periodic_root_operator(plan)
    type(fmm_plan_type), intent(inout) :: plan

    if (allocated(plan%periodic_root_operator)) deallocate (plan%periodic_root_operator)
    plan%periodic_root_operator_ready = .false.

    if (use_periodic2_m2l_root_trunc(plan)) then
      call precompute_periodic_root_trunc_operator(plan)
    else if (use_periodic2_m2l_root_oracle(plan)) then
      call precompute_periodic_root_oracle_operator(plan)
    end if
  end subroutine precompute_periodic_root_operator

  subroutine precompute_periodic_root_trunc_operator(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: nimg, outer_img, img_i, img_j
    integer(i32) :: axis1, axis2, alpha_idx, beta_idx, deriv_idx, monopole_coef
    real(dp) :: root_center_offset(3), rimg(3), deriv_sum(plan%nderiv), deriv_tmp(plan%nderiv)

    if (plan%ncoef <= 0_i32 .or. plan%nderiv <= 0_i32) return

    nimg = max(0_i32, plan%options%periodic_image_layers)
    outer_img = nimg + plan%options%periodic_ewald_layers
    axis1 = plan%options%periodic_axes(1)
    axis2 = plan%options%periodic_axes(2)
    if (outer_img <= nimg .or. axis1 <= 0_i32 .or. axis2 <= 0_i32) return

    allocate (plan%periodic_root_operator(plan%ncoef, plan%ncoef))
    plan%periodic_root_operator = 0.0d0
    deriv_sum = 0.0d0

    ! Precompute a truncated far-image root operator. This approximates the periodic
    ! residual outside the finite image shell already handled explicitly at runtime.
    ! It is not an exact Ewald operator.
    root_center_offset = active_tree_node_center(plan, plan%target_tree_ready, 1_i32) - plan%node_center(:, 1_i32)
    do img_i = -outer_img, outer_img
      do img_j = -outer_img, outer_img
        if (abs(img_i) <= nimg .and. abs(img_j) <= nimg) cycle
        rimg = root_center_offset
        rimg(axis1) = rimg(axis1) - real(img_i, dp) * plan%options%periodic_len(1)
        rimg(axis2) = rimg(axis2) - real(img_j, dp) * plan%options%periodic_len(2)
        call compute_laplace_derivatives(plan, rimg, deriv_tmp)
        deriv_sum = deriv_sum + deriv_tmp
      end do
    end do

    ! Drop the whole monopole source column so this stays a truncated
    ! neutral-residual operator rather than a charged-wall oracle.
    monopole_coef = plan%alpha_map(0_i32, 0_i32, 0_i32)
    do beta_idx = 1_i32, plan%ncoef
      if (beta_idx == monopole_coef) then
        plan%periodic_root_operator(:, beta_idx) = 0.0d0
        cycle
      end if
      do alpha_idx = 1_i32, plan%ncoef
        deriv_idx = plan%alpha_beta_deriv_idx(alpha_idx, beta_idx)
        plan%periodic_root_operator(alpha_idx, beta_idx) = plan%alpha_sign(beta_idx) * deriv_sum(deriv_idx)
      end do
    end do
    plan%periodic_root_operator_ready = .true.
  end subroutine precompute_periodic_root_trunc_operator

  subroutine precompute_periodic_root_oracle_operator(plan)
    type(fmm_plan_type), intent(inout) :: plan
    integer(i32) :: nproxy, ncheck, j, i
    real(dp) :: source_center(3), source_half(3), target_center(3), target_half(3)
    real(dp), allocatable :: proxy_points(:, :), check_points(:, :)
    real(dp), allocatable :: proxy_to_multipole(:, :), proxy_to_local(:, :)
    real(dp), allocatable :: field_matrix(:, :), field_rhs(:)
    real(dp), allocatable :: normal_matrix(:, :), normal_rhs(:), coeff(:)
    real(dp), allocatable :: proxy_gram(:, :), inv_proxy_gram(:, :), proxy_pinv(:, :)
    real(dp) :: e_res(3)

    if (.not. plan%periodic_ewald%ready) return
    if (plan%ncoef <= 1_i32) return

    nproxy = 4_i32 * plan%ncoef
    ncheck = 8_i32 * plan%ncoef
    source_center = plan%node_center(:, 1_i32)
    source_half = plan%node_half_size(:, 1_i32)
    target_center = active_tree_node_center(plan, plan%target_tree_ready, 1_i32)
    target_half = active_tree_node_half_size(plan, plan%target_tree_ready, 1_i32)

    allocate (proxy_points(3, nproxy), check_points(3, ncheck))
    call build_root_sample_points(source_center, source_half, nproxy, 0.13d0, proxy_points)
    call build_root_sample_points(target_center, target_half, ncheck, 0.37d0, check_points)

    allocate (proxy_to_multipole(plan%ncoef, nproxy), proxy_to_local(plan%ncoef, nproxy))
    allocate (field_matrix(3_i32 * ncheck, plan%ncoef - 1_i32), normal_matrix(plan%ncoef - 1_i32, plan%ncoef - 1_i32))
    allocate (field_rhs(3_i32 * ncheck), normal_rhs(plan%ncoef - 1_i32), coeff(plan%ncoef - 1_i32))
    allocate (proxy_gram(plan%ncoef, plan%ncoef), inv_proxy_gram(plan%ncoef, plan%ncoef), proxy_pinv(nproxy, plan%ncoef))

    call build_proxy_multipole_matrix(plan, source_center, proxy_points, proxy_to_multipole)
    call build_local_field_matrix(plan, target_center, check_points, field_matrix)
    normal_matrix = matmul(transpose(field_matrix), field_matrix)
    call regularize_symmetric_system(normal_matrix)

    proxy_to_local = 0.0d0
    do j = 1_i32, nproxy
      field_rhs = 0.0d0
      do i = 1_i32, ncheck
        e_res = 0.0d0
        call add_periodic2_exact_ewald_correction_single_source(plan, 1.0d0, proxy_points(:, j), check_points(:, i), e_res)
        field_rhs(i) = e_res(1)
        field_rhs(ncheck + i) = e_res(2)
        field_rhs(2_i32 * ncheck + i) = e_res(3)
      end do
      normal_rhs = matmul(transpose(field_matrix), field_rhs)
      call solve_square_system(normal_matrix, normal_rhs, coeff)
      proxy_to_local(2:plan%ncoef, j) = coeff
    end do

    proxy_gram = matmul(proxy_to_multipole, transpose(proxy_to_multipole))
    call regularize_symmetric_system(proxy_gram)
    call invert_square_matrix(proxy_gram, inv_proxy_gram)
    proxy_pinv = matmul(transpose(proxy_to_multipole), inv_proxy_gram)
    allocate (plan%periodic_root_operator(plan%ncoef, plan%ncoef))
    plan%periodic_root_operator = matmul(proxy_to_local, proxy_pinv)
    plan%periodic_root_operator(1_i32, :) = 0.0d0
    plan%periodic_root_operator_ready = .true.
  end subroutine precompute_periodic_root_oracle_operator

  subroutine build_root_sample_points(center, half_size, npoint, offset, points)
    real(dp), intent(in) :: center(3), half_size(3), offset
    integer(i32), intent(in) :: npoint
    real(dp), intent(out) :: points(3, npoint)
    integer(i32) :: idx
    real(dp) :: f1, f2, f3
    real(dp), parameter :: g1 = 0.7548776662466927d0
    real(dp), parameter :: g2 = 0.5698402909980532d0
    real(dp), parameter :: g3 = 0.4301597090019468d0

    do idx = 1_i32, npoint
      f1 = modulo(offset + real(idx, dp) * g1, 1.0d0)
      f2 = modulo(offset + real(idx, dp) * g2, 1.0d0)
      f3 = modulo(offset + real(idx, dp) * g3, 1.0d0)
      points(1, idx) = center(1) + (2.0d0 * (0.05d0 + 0.9d0 * f1) - 1.0d0) * half_size(1)
      points(2, idx) = center(2) + (2.0d0 * (0.05d0 + 0.9d0 * f2) - 1.0d0) * half_size(2)
      points(3, idx) = center(3) + (2.0d0 * (0.05d0 + 0.9d0 * f3) - 1.0d0) * half_size(3)
    end do
  end subroutine build_root_sample_points

  subroutine build_proxy_multipole_matrix(plan, source_center, proxy_points, matrix)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: source_center(3), proxy_points(:, :)
    real(dp), intent(out) :: matrix(:, :)
    integer(i32) :: proxy_idx, beta_idx
    real(dp) :: d(3)
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    matrix = 0.0d0
    do proxy_idx = 1_i32, int(size(proxy_points, 2), i32)
      d = proxy_points(:, proxy_idx) - source_center
      call build_axis_powers(d, plan%options%order, xpow, ypow, zpow)
      do beta_idx = 1_i32, plan%ncoef
        matrix(beta_idx, proxy_idx) = xpow(plan%alpha(1, beta_idx)) * ypow(plan%alpha(2, beta_idx)) &
                                      * zpow(plan%alpha(3, beta_idx)) / plan%alpha_factorial(beta_idx)
      end do
    end do
  end subroutine build_proxy_multipole_matrix

  subroutine build_local_field_matrix(plan, target_center, check_points, matrix)
    type(fmm_plan_type), intent(in) :: plan
    real(dp), intent(in) :: target_center(3), check_points(:, :)
    real(dp), intent(out) :: matrix(:, :)
    integer(i32) :: check_idx, term_idx, coeff_idx, ncheck
    real(dp) :: d(3), monomial
    real(dp) :: xpow(0:max(0_i32, plan%options%order)), ypow(0:max(0_i32, plan%options%order))
    real(dp) :: zpow(0:max(0_i32, plan%options%order))

    matrix = 0.0d0
    ncheck = int(size(check_points, 2), i32)
    do check_idx = 1_i32, ncheck
      d = check_points(:, check_idx) - target_center
      call build_axis_powers(d, plan%options%order, xpow, ypow, zpow)
      do term_idx = 1_i32, plan%eval_term_count
        monomial = xpow(plan%eval_exp(1, term_idx)) * ypow(plan%eval_exp(2, term_idx)) * zpow(plan%eval_exp(3, term_idx)) * &
                   plan%eval_inv_factorial(term_idx)
        coeff_idx = plan%eval_deriv_idx(1, term_idx)
        if (coeff_idx > 1_i32) matrix(check_idx, coeff_idx - 1_i32) = matrix(check_idx, coeff_idx - 1_i32) - monomial
        coeff_idx = plan%eval_deriv_idx(2, term_idx)
        if (coeff_idx > 1_i32) then
          matrix(ncheck + check_idx, coeff_idx - 1_i32) = matrix(ncheck + check_idx, coeff_idx - 1_i32) - monomial
        end if
        coeff_idx = plan%eval_deriv_idx(3, term_idx)
        if (coeff_idx > 1_i32) matrix(2_i32 * ncheck + check_idx, coeff_idx - 1_i32) = &
          matrix(2_i32 * ncheck + check_idx, coeff_idx - 1_i32) - monomial
      end do
    end do
  end subroutine build_local_field_matrix

  subroutine regularize_symmetric_system(matrix)
    real(dp), intent(inout) :: matrix(:, :)
    integer(i32) :: idx
    real(dp) :: scale

    do idx = 1_i32, int(size(matrix, 1), i32)
      scale = max(1.0d0, sum(abs(matrix(idx, :))))
      matrix(idx, idx) = matrix(idx, idx) + 1.0d-12 * scale
    end do
  end subroutine regularize_symmetric_system

  subroutine invert_square_matrix(matrix, inverse)
    real(dp), intent(in) :: matrix(:, :)
    real(dp), intent(out) :: inverse(:, :)
    integer(i32) :: ncol, rhs_idx
    real(dp), allocatable :: rhs(:), sol(:)

    ncol = int(size(matrix, 1), i32)
    if (size(matrix, 2) /= ncol .or. size(inverse, 1) /= ncol .or. size(inverse, 2) /= ncol) then
      error stop 'invert_square_matrix dimension mismatch.'
    end if

    allocate (rhs(ncol), sol(ncol))
    do rhs_idx = 1_i32, ncol
      rhs = 0.0d0
      rhs(rhs_idx) = 1.0d0
      call solve_square_system(matrix, rhs, sol)
      inverse(:, rhs_idx) = sol
    end do
  end subroutine invert_square_matrix

  subroutine solve_square_system(matrix, rhs, solution)
    real(dp), intent(in) :: matrix(:, :)
    real(dp), intent(in) :: rhs(:)
    real(dp), intent(out) :: solution(:)
    integer(i32) :: ncol, pivot_row, row_idx, col_idx, swap_idx
    real(dp), allocatable :: work(:, :), rhs_work(:), row_tmp(:)
    real(dp) :: pivot_abs, factor

    ncol = int(size(matrix, 1), i32)
    if (size(matrix, 2) /= ncol .or. size(rhs) /= ncol .or. size(solution) /= ncol) then
      error stop 'solve_square_system dimension mismatch.'
    end if

    allocate (work(ncol, ncol), rhs_work(ncol), row_tmp(ncol))
    work = matrix
    rhs_work = rhs

    do col_idx = 1_i32, ncol
      pivot_row = col_idx
      pivot_abs = abs(work(col_idx, col_idx))
      do row_idx = col_idx + 1_i32, ncol
        if (abs(work(row_idx, col_idx)) > pivot_abs) then
          pivot_abs = abs(work(row_idx, col_idx))
          pivot_row = row_idx
        end if
      end do
      if (pivot_abs <= 1.0d-20) error stop 'periodic root oracle linear system is singular.'

      if (pivot_row /= col_idx) then
        row_tmp = work(col_idx, :)
        work(col_idx, :) = work(pivot_row, :)
        work(pivot_row, :) = row_tmp
        factor = rhs_work(col_idx)
        rhs_work(col_idx) = rhs_work(pivot_row)
        rhs_work(pivot_row) = factor
      end if

      factor = work(col_idx, col_idx)
      work(col_idx, col_idx:ncol) = work(col_idx, col_idx:ncol) / factor
      rhs_work(col_idx) = rhs_work(col_idx) / factor
      do row_idx = col_idx + 1_i32, ncol
        factor = work(row_idx, col_idx)
        if (abs(factor) <= tiny(1.0d0)) cycle
        work(row_idx, col_idx:ncol) = work(row_idx, col_idx:ncol) - factor * work(col_idx, col_idx:ncol)
        rhs_work(row_idx) = rhs_work(row_idx) - factor * rhs_work(col_idx)
      end do
    end do

    solution = rhs_work
    do row_idx = ncol, 1_i32, -1_i32
      do swap_idx = row_idx + 1_i32, ncol
        solution(row_idx) = solution(row_idx) - work(row_idx, swap_idx) * solution(swap_idx)
      end do
    end do
  end subroutine solve_square_system

end module bem_coulomb_fmm_periodic_root_ops

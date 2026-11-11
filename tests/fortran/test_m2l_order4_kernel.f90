!> 展開次数 4 の生成 M2L カーネルを、plan の添字表を使う汎用ループと比べる。
program test_m2l_order4_kernel
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_types, only: fmm_plan_type
  use bem_coulomb_fmm_basis, only: initialize_basis_tables
  use bem_coulomb_fmm_m2l_order4, only: m2l_lane_step_order4
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true
  implicit none
  type(fmm_plan_type), allocatable :: plan
  real(dp) :: coeff(4, 35), deriv(4, 165), acc(4, 35), reference(4, 35)
  integer(i32) :: alpha_idx, beta_idx

  call test_init(1)
  call test_begin('generated_kernel_matches_plan_index_table')
  allocate (plan)
  call initialize_basis_tables(plan, 4_i32)
  call assert_true(plan%ncoef == 35_i32 .and. plan%nderiv == 165_i32, 'order 4 has 35 coefficients and 165 derivatives')
  ! Distinct values in every slot, so any swapped index shows up.
  coeff = reshape([(sin(1.3_dp*real(alpha_idx, dp)), alpha_idx=1, size(coeff))], shape(coeff))
  deriv = reshape([(cos(0.7_dp*real(alpha_idx, dp)), alpha_idx=1, size(deriv))], shape(deriv))
  acc = reshape([(0.01_dp*real(alpha_idx, dp), alpha_idx=1, size(acc))], shape(acc))
  reference = acc
  do beta_idx = 1_i32, 35_i32
    do alpha_idx = 1_i32, 35_i32
      reference(:, alpha_idx) = reference(:, alpha_idx) + &
                                coeff(:, beta_idx)*deriv(:, plan%alpha_beta_deriv_idx(alpha_idx, beta_idx))
    end do
  end do
  call m2l_lane_step_order4(coeff, deriv, acc)
  ! A wrong derivative index changes a coefficient by O(1); only rounding (FMA contraction) may differ.
  call assert_true(maxval(abs(acc - reference)) <= 1.0e-13_dp, &
                   'generated order-4 kernel does not follow alpha_beta_deriv_idx')
  call test_end()
  call test_summary()
end program test_m2l_order4_kernel

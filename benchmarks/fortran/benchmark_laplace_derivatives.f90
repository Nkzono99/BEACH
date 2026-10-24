!> Laplace 偏微分の速度比較。正しさは test_laplace_derivatives で検証する。
program benchmark_laplace_derivatives
  use bem_kinds, only: dp, i32, i64
  use bem_coulomb_fmm_types, only: fmm_plan_type
  use bem_coulomb_fmm_basis, only: initialize_basis_tables, compute_laplace_derivatives, &
                                   compute_laplace_derivatives_reference
  implicit none

  integer(i32), parameter :: repeat_count = 2000_i32
  type(fmm_plan_type), allocatable :: plan
  real(dp), allocatable :: derivatives(:)
  real(dp) :: r(3), fast_seconds, reference_seconds, checksum
  integer(i32) :: repeat
  integer(i64) :: clock0, clock1, rate

  allocate (plan)
  plan%options%order = 4_i32
  call initialize_basis_tables(plan, 4_i32)
  allocate (derivatives(plan%nderiv))
  r = [1.0e-4_dp, 2.3e-5_dp, -7.1e-5_dp]
  checksum = 0.0_dp
  call compute_laplace_derivatives(plan, r, derivatives)
  call compute_laplace_derivatives_reference(plan, r, derivatives)

  call system_clock(clock0, rate)
  do repeat = 1_i32, repeat_count
    call compute_laplace_derivatives(plan, r + 1.0e-9_dp*real(repeat, dp), derivatives)
    checksum = checksum + derivatives(1)
  end do
  call system_clock(clock1)
  fast_seconds = real(clock1 - clock0, dp)/real(max(rate, 1_i64), dp)

  call system_clock(clock0, rate)
  do repeat = 1_i32, repeat_count
    call compute_laplace_derivatives_reference(plan, r + 1.0e-9_dp*real(repeat, dp), derivatives)
    checksum = checksum + derivatives(1)
  end do
  call system_clock(clock1)
  reference_seconds = real(clock1 - clock0, dp)/real(max(rate, 1_i64), dp)

  print '(a,i0)', 'evaluations_per_method=', repeat_count
  print '(a,es16.8)', 'recurrence_seconds=', fast_seconds
  print '(a,es16.8)', 'polynomial_seconds=', reference_seconds
  print '(a,es16.8)', 'speedup=', reference_seconds/max(fast_seconds, tiny(1.0_dp))
  print '(a,es16.8)', 'checksum=', checksum
end program benchmark_laplace_derivatives

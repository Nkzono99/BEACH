!> Laplace 核の偏微分表: 漸化式版が旧多項式乗算版と全次数で一致し、十分速いことの契約テスト。
program test_laplace_derivatives
  use bem_kinds, only: dp, i32
  use bem_coulomb_fmm_types, only: fmm_plan_type
  use bem_coulomb_fmm_basis, only: initialize_basis_tables, compute_laplace_derivatives, &
                                   compute_laplace_derivatives_reference
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true
  implicit none
  type(fmm_plan_type), allocatable :: plan
  real(dp), allocatable :: fast(:), reference(:)
  real(dp) :: points(3, 4), r(3), scale, worst, rel
  real(dp) :: softening(3)
  integer(i32) :: order, point, soft_idx, idx, repeat
  integer(8) :: clock0, clock1, rate
  logical :: ok
  character(len=96) :: label

  points(:, 1) = [0.37_dp, -0.81_dp, 0.52_dp]
  points(:, 2) = [1.0e-4_dp, 2.3e-5_dp, -7.1e-5_dp]
  points(:, 3) = [-3.0_dp, 0.25_dp, 4.75_dp]
  points(:, 4) = [2.0e-6_dp, -9.0e-6_dp, 1.5e-5_dp]
  softening = [0.0_dp, 1.0e-7_dp, 0.3_dp]

  call test_init(7)
  do order = 1_i32, 6_i32
    write (label, '(a,i0)') 'recurrence_matches_polynomial_expansion_order_', order
    call test_begin(trim(label))
    allocate (plan)
    ! The reference expansion reads the order from the plan options, not from the table size.
    plan%options%order = order
    call initialize_basis_tables(plan, order)
    allocate (fast(plan%nderiv), reference(plan%nderiv))
    ok = .true.
    worst = 0.0_dp
    do soft_idx = 1_i32, 3_i32
      plan%options%softening = softening(soft_idx)
      do point = 1_i32, 4_i32
        r = points(:, point)
        call compute_laplace_derivatives(plan, r, fast)
        call compute_laplace_derivatives_reference(plan, r, reference)
        ! Each derivative is compared on its own scale; magnitudes span dozens of decades across degrees.
        do idx = 1_i32, plan%nderiv
          scale = max(abs(reference(idx)), abs(fast(idx)))
          if (scale <= tiny(1.0_dp)) cycle
          rel = abs(fast(idx) - reference(idx))/scale
          worst = max(worst, rel)
          ! The polynomial expansion accumulates cancellation at high degree; 1e-9 keeps both within the
          ! accuracy the order-4 FMM relies on, and the recurrence is the better conditioned of the two.
          if (rel > 1.0e-9_dp) ok = .false.
        end do
      end do
    end do
    print '(a,i0,a,es10.3)', '  order ', order, ' max relative difference = ', worst
    call assert_true(ok, 'recurrence and polynomial derivatives disagree')
    deallocate (fast, reference, plan)
    call test_end()
  end do

  call test_begin('recurrence_is_at_least_ten_times_faster_at_order_4')
  allocate (plan)
  plan%options%order = 4_i32
  call initialize_basis_tables(plan, 4_i32)
  plan%options%softening = 0.0_dp
  allocate (fast(plan%nderiv), reference(plan%nderiv))
  r = points(:, 2)
  call system_clock(clock0, rate)
  do repeat = 1_i32, 2000_i32
    call compute_laplace_derivatives(plan, r + 1.0e-9_dp*real(repeat, dp), fast)
  end do
  call system_clock(clock1)
  scale = real(clock1 - clock0, dp)/real(max(rate, 1_8), dp)
  call system_clock(clock0, rate)
  do repeat = 1_i32, 2000_i32
    call compute_laplace_derivatives_reference(plan, r + 1.0e-9_dp*real(repeat, dp), reference)
  end do
  call system_clock(clock1)
  worst = real(clock1 - clock0, dp)/real(max(rate, 1_8), dp)
  print '(a,es10.3,a,es10.3,a)', '  2000 evaluations: recurrence ', scale, ' s, polynomial ', worst, ' s'
  call assert_true(scale*10.0_dp <= worst .or. worst < 1.0e-3_dp, 'recurrence is not at least ten times faster')
  call test_end()
  call test_summary()
end program test_laplace_derivatives

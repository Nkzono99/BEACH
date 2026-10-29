!> PCG64 の counter-based 乱数と逐次乱数列、それを使う注入 sampler の分布を検証するテスト。
!! KS 検定は有意水準 0.1% の臨界値 1.95/sqrt(n) を使う（seed 固定なので結果は決定的）。
program test_counter_rng
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: k_boltzmann
  use bem_counter_rng, only: counter_uniform, counter_uniforms, counter_child_id, open_uniform_from_double53
  use bem_random_stream, only: seed_random_stream, random_uniform
  use bem_injection, only: sample_reservoir_face_particle_counter
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp
  implicit none

  integer, parameter :: n_uniform = 200000
  integer, parameter :: n_sampler = 20000
  real(dp), allocatable :: u(:), u2(:)
  integer :: i
  real(dp) :: mean, var, corr, d, run(5)

  call test_init(6)

  call test_begin('counter_uniform_is_deterministic_and_open')
  call assert_close_dp(counter_uniform(7_i32, [1_i64, 2_i64, 3_i64], 5_i64), &
                       counter_uniform(7_i32, [1_i64, 2_i64, 3_i64], 5_i64), 0.0_dp, 'same key must give the same value')
  call assert_true(counter_uniform(7_i32, [1_i64, 2_i64, 3_i64], 5_i64) /= &
                   counter_uniform(7_i32, [1_i64, 2_i64, 3_i64], 6_i64), 'different draws must differ')
  call assert_true(counter_uniform(7_i32, [1_i64, 2_i64, 3_i64], 5_i64) /= &
                   counter_uniform(8_i32, [1_i64, 2_i64, 3_i64], 5_i64), 'different seeds must differ')
  call assert_true(counter_uniform(7_i32, [1_i64, 2_i64], 5_i64) /= &
                   counter_uniform(7_i32, [2_i64, 1_i64], 5_i64), 'counter order must matter')
  call assert_true(counter_child_id(11_i64, 1_i64) /= counter_child_id(11_i64, 2_i64), 'child ids must differ by ordinal')
  call assert_true(counter_child_id(11_i64, 1_i64) >= 0_i64, 'child ids must be nonnegative')
  ! Drawing a run of values from one opened generator must reproduce the single-draw values.
  call counter_uniforms(7_i32, [1_i64, 2_i64, 3_i64], 98_i64, run)
  do i = 1, size(run)
    call assert_close_dp(run(i), counter_uniform(7_i32, [1_i64, 2_i64, 3_i64], int(97 + i, i64)), 0.0_dp, &
                         'counter_uniforms must match counter_uniform draw by draw')
  end do
  call counter_uniforms(7_i32, [1_i64, 2_i64, 3_i64], 0_i64, run(:1))
  call assert_close_dp(run(1), counter_uniform(7_i32, [1_i64, 2_i64, 3_i64], 0_i64), 0.0_dp, &
                       'counter_uniforms draw 0 mismatch')
  call test_end()

  call test_begin('open_uniform_lattice_endpoints')
  call check_open_uniform_endpoints()
  call test_end()

  call test_begin('counter_uniform_statistics')
  allocate (u(n_uniform), u2(n_uniform))
  do i = 1, n_uniform
    ! Consecutive slots of one segment, the way batch injection keys its particles.
    u(i) = counter_uniform(12345_i32, [3_i64, 1_i64, 3_i64, 5_i64, int(i, i64)], 1_i64)
    u2(i) = counter_uniform(12345_i32, [3_i64, 1_i64, 3_i64, 5_i64, int(i, i64)], 2_i64)
  end do
  call check_uniform_sample(u, 'counter draws over slots')
  corr = correlation(u(1:n_uniform - 1), u(2:n_uniform))
  call assert_true(abs(corr) < 5.0_dp/sqrt(real(n_uniform, dp)), 'neighbouring slots must be uncorrelated')
  corr = correlation(u, u2)
  call assert_true(abs(corr) < 5.0_dp/sqrt(real(n_uniform, dp)), 'draws 1 and 2 of one slot must be uncorrelated')
  do i = 1, n_uniform
    ! Many draws of one key, the way a slot's thermal velocity uses draws 100, 101, ...
    u(i) = counter_uniform(12345_i32, [9_i64, 9_i64], int(i, i64))
  end do
  call check_uniform_sample(u, 'counter draws within one key')
  call test_end()

  call test_begin('random_stream_reproducible_and_uniform')
  call seed_random_stream([4242_i32])
  call random_uniform(u)
  call seed_random_stream([4242_i32])
  call random_uniform(u2)
  call assert_true(all(u == u2), 'the same seed must reproduce the stream')
  call check_uniform_sample(u, 'sequential stream')
  call seed_random_stream([4243_i32])
  call random_uniform(u2)
  call assert_true(count(u == u2) < 10, 'different seeds must give different streams')
  call test_end()

  call test_begin('counter_reservoir_normal_speed_follows_flux_weighted_maxwellian')
  call run_reservoir_distribution_test()
  call test_end()

  call test_begin('counter_reservoir_drifting_flux_matches_analytic_mean')
  call run_reservoir_drift_test()
  call test_end()

  call test_summary()

contains

  !> 53 bit 格子 k/2^53 から (0, 1) への写像を両端と通常値で決定的に検査する（#72）。
  !! 最上端 k = 2^53-1 は 2^-54 を足すと 1.0 に丸まるため、1 の直前値 1-2^-53 に留まる必要がある。
  !! それ以外の格子値は従来の `raw + 2^-54` と同じ binary64 値を保つ。
  subroutine check_open_uniform_endpoints()
    real(dp), parameter :: two53 = 9007199254740992.0_dp
    real(dp), parameter :: half_ulp = 0.5_dp**54
    real(dp) :: raw_min, raw_max, ordinary(4), uniform
    integer :: k

    raw_min = 0.0_dp
    uniform = open_uniform_from_double53(raw_min)
    call assert_true(uniform > 0.0_dp .and. uniform < 1.0_dp, 'raw = 0 must map inside (0, 1)')
    call assert_true(uniform == half_ulp, 'raw = 0 must map to 2^-54')

    raw_max = (two53 - 1.0_dp)/two53
    uniform = open_uniform_from_double53(raw_max)
    call assert_true(uniform > 0.0_dp .and. uniform < 1.0_dp, 'raw = (2^53-1)/2^53 must map inside (0, 1)')
    call assert_true(uniform == 1.0_dp - 0.5_dp**53, 'raw = (2^53-1)/2^53 must map to 1-2^-53')

    ! Lower half (exact midpoint), upper half with even and odd k (round-half-even), and the largest k below the top.
    ordinary = [0.25_dp, 0.75_dp, 0.75_dp + 1.0_dp/two53, (two53 - 2.0_dp)/two53]
    do k = 1, size(ordinary)
      uniform = open_uniform_from_double53(ordinary(k))
      call assert_true(uniform == ordinary(k) + half_ulp, 'ordinary lattice values must keep raw + 2^-54')
      call assert_true(uniform > 0.0_dp .and. uniform < 1.0_dp, 'ordinary lattice values must stay inside (0, 1)')
    end do
  end subroutine check_open_uniform_endpoints

  subroutine check_uniform_sample(sample, label)
    real(dp), intent(in) :: sample(:)
    character(len=*), intent(in) :: label
    real(dp), allocatable :: sorted(:)
    integer :: k, n

    n = size(sample)
    call assert_true(all(sample > 0.0_dp .and. sample < 1.0_dp), trim(label)//': values must lie in (0, 1)')
    mean = sum(sample)/real(n, dp)
    var = sum((sample - mean)**2)/real(n - 1, dp)
    call assert_close_dp(mean, 0.5_dp, 5.0_dp*sqrt(1.0_dp/(12.0_dp*real(n, dp))), trim(label)//': mean')
    call assert_close_dp(var, 1.0_dp/12.0_dp, 5.0_dp*sqrt(1.0_dp/(180.0_dp*real(n, dp))), trim(label)//': variance')
    sorted = sample
    call sort_real(sorted)
    d = 0.0_dp
    do k = 1, n
      d = max(d, abs(real(k, dp)/real(n, dp) - sorted(k)), abs(sorted(k) - real(k - 1, dp)/real(n, dp)))
    end do
    call assert_true(d < 1.95_dp/sqrt(real(n, dp)), trim(label)//': KS statistic against U(0,1)')
  end subroutine check_uniform_sample

  !> 零ドリフトの flux 重み付き Maxwell 分布: 法線速さの CDF は 1 - exp(-v^2/2 sigma^2)、接線成分は N(0, sigma^2)。
  subroutine run_reservoir_distribution_test()
    real(dp), parameter :: m = 1.0_dp, temperature = 1.0_dp/k_boltzmann
    real(dp) :: x(3), v(3), sigma
    real(dp), allocatable :: vn(:), vt(:)
    integer :: k

    sigma = sqrt(k_boltzmann*temperature/m)
    allocate (vn(n_sampler), vt(n_sampler))
    do k = 1, n_sampler
      call sample_reservoir_face_particle_counter( &
        [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 'z_low', [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 0.0_dp], &
        [0.0_dp, 0.0_dp, 0.0_dp], m, temperature, 0.0_dp, 0.0_dp, 777_i32, [1_i64, 2_i64, 3_i64, 5_i64, int(k, i64)], x, v &
        )
      vn(k) = v(3)
      vt(k) = v(1)
    end do
    call assert_true(all(vn > 0.0_dp), 'z_low inflow must move in +z')
    call sort_real(vn)
    call sort_real(vt)
    d = 0.0_dp
    do k = 1, n_sampler
      d = max(d, ks_gap(k, n_sampler, 1.0_dp - exp(-0.5_dp*(vn(k)/sigma)**2)))
    end do
    call assert_true(d < 1.95_dp/sqrt(real(n_sampler, dp)), 'normal speed KS statistic against the flux-weighted law')
    d = 0.0_dp
    do k = 1, n_sampler
      d = max(d, ks_gap(k, n_sampler, 0.5_dp*(1.0_dp + erf(vt(k)/(sqrt(2.0_dp)*sigma)))))
    end do
    call assert_true(d < 1.95_dp/sqrt(real(n_sampler, dp)), 'tangential velocity KS statistic against N(0, sigma^2)')
  end subroutine run_reservoir_distribution_test

  !> ドリフト a（sigma 単位）の flux 重み付き分布の法線速さ平均を解析式と比べる。
  subroutine run_reservoir_drift_test()
    real(dp), parameter :: m = 1.0_dp, temperature = 1.0_dp/k_boltzmann, drift = 1.5_dp
    real(dp) :: x(3), v(3), a, phi, big_phi, flux, second, expected_mean, sample_mean, sample_sq
    integer :: k

    sample_mean = 0.0_dp
    sample_sq = 0.0_dp
    do k = 1, n_sampler
      call sample_reservoir_face_particle_counter( &
        [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 'z_low', [0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 0.0_dp], &
        [0.0_dp, 0.0_dp, drift], m, temperature, 0.0_dp, 0.0_dp, 778_i32, [1_i64, 2_i64, 3_i64, 6_i64, int(k, i64)], x, v &
        )
      sample_mean = sample_mean + v(3)
      sample_sq = sample_sq + v(3)**2
    end do
    sample_mean = sample_mean/real(n_sampler, dp)
    sample_sq = sample_sq/real(n_sampler, dp)
    ! sigma = 1, v > 0: int v phi(v - a) = phi(a) + a Phi(a), int v^2 phi(v - a) = (1 + a^2) Phi(a) + a phi(a).
    a = drift
    phi = exp(-0.5_dp*a*a)/sqrt(2.0_dp*acos(-1.0_dp))
    big_phi = 0.5_dp*(1.0_dp + erf(a/sqrt(2.0_dp)))
    flux = phi + a*big_phi
    second = (1.0_dp + a*a)*big_phi + a*phi
    expected_mean = second/flux
    call assert_close_dp(sample_mean, expected_mean, 5.0_dp*sqrt((sample_sq - sample_mean**2)/real(n_sampler, dp)), &
                         'drifting flux-weighted mean normal speed')
  end subroutine run_reservoir_drift_test

  pure real(dp) function ks_gap(k, n, cdf) result(gap)
    integer, intent(in) :: k, n
    real(dp), intent(in) :: cdf

    gap = max(abs(real(k, dp)/real(n, dp) - cdf), abs(cdf - real(k - 1, dp)/real(n, dp)))
  end function ks_gap

  real(dp) function correlation(a, b) result(r)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: ma, mb

    ma = sum(a)/real(size(a), dp)
    mb = sum(b)/real(size(b), dp)
    r = sum((a - ma)*(b - mb))/sqrt(sum((a - ma)**2)*sum((b - mb)**2))
  end function correlation

  !> 昇順の heap sort。
  subroutine sort_real(a)
    real(dp), intent(inout) :: a(:)
    integer :: n, k
    real(dp) :: tmp

    n = size(a)
    do k = n/2, 1, -1
      call sift_down(a, k, n)
    end do
    do k = n, 2, -1
      tmp = a(1)
      a(1) = a(k)
      a(k) = tmp
      call sift_down(a, 1, k - 1)
    end do
  end subroutine sort_real

  subroutine sift_down(a, start, last)
    real(dp), intent(inout) :: a(:)
    integer, intent(in) :: start, last
    integer :: root, child
    real(dp) :: tmp

    root = start
    do while (2*root <= last)
      child = 2*root
      if (child < last) then
        if (a(child) < a(child + 1)) child = child + 1
      end if
      if (a(root) >= a(child)) return
      tmp = a(root)
      a(root) = a(child)
      a(child) = tmp
      root = child
    end do
  end subroutine sift_down

end program test_counter_rng

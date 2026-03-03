!> 乱数シード設定と粒子位置/速度サンプリングを担う粒子注入モジュール。
module bem_injection
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
  use bem_particles, only: init_particles
  use bem_types, only: particles_soa
  implicit none

  private
  public :: seed_rng
  public :: sample_uniform_positions
  public :: sample_shifted_maxwell_velocities
  public :: init_random_beam_particles

contains

  !> 与えたシード列またはシステム時刻からFortran乱数生成器を初期化する。
  !! @param[in] seed 入力引数。
  subroutine seed_rng(seed)
    integer(i32), intent(in), optional :: seed(:)
    integer :: n, i, clk
    integer, allocatable :: put(:)

    call random_seed(size=n)
    allocate(put(n))

    if (present(seed)) then
      do i = 1, n
        put(i) = seed(mod(i - 1, size(seed)) + 1) + 104729 * i
      end do
    else
      call system_clock(count=clk)
      do i = 1, n
        put(i) = clk + 37 * i
      end do
    end if

    call random_seed(put=put)
  end subroutine seed_rng

  !> 直方体領域 `[low, high]` 内で一様分布の初期位置をサンプリングする。
  !! @param[in] low 入力引数。
  !! @param[in] high 入力引数。
  !! @param[out] x 出力引数。
  subroutine sample_uniform_positions(low, high, x)
    real(dp), intent(in) :: low(3), high(3)
    real(dp), intent(out) :: x(:, :)
    real(dp), allocatable :: u(:, :)

    if (size(x, 1) /= 3) error stop "x first dimension must be 3"
    if (any(high < low)) error stop "high must be >= low for all axes"

    allocate(u(3, size(x, 2)))
    call random_number(u)
    x = spread(low, dim=2, ncopies=size(x, 2)) + spread(high - low, dim=2, ncopies=size(x, 2)) * u
  end subroutine sample_uniform_positions

  !> ドリフト速度付きMaxwell分布(温度または熱速度指定)から粒子速度を生成する。
  !! @param[in] drift_velocity 入力引数。
  !! @param[in] m_particle 入力引数。
  !! @param[out] v 出力引数。
  !! @param[in] temperature_k 入力引数。
  !! @param[in] thermal_speed 入力引数。
  subroutine sample_shifted_maxwell_velocities(drift_velocity, m_particle, v, temperature_k, thermal_speed)
    real(dp), intent(in) :: drift_velocity(3)
    real(dp), intent(in) :: m_particle
    real(dp), intent(out) :: v(:, :)
    real(dp), intent(in), optional :: temperature_k
    real(dp), intent(in), optional :: thermal_speed
    integer :: n
    real(dp) :: sigma
    real(dp), allocatable :: z(:, :)

    if (size(v, 1) /= 3) error stop "v first dimension must be 3"
    if (m_particle <= 0.0_dp) error stop "m_particle must be > 0"

    if (.not. present(temperature_k) .and. .not. present(thermal_speed)) then
      error stop "either temperature_k or thermal_speed must be provided"
    end if

    if (present(thermal_speed)) then
      if (thermal_speed < 0.0_dp) error stop "thermal_speed must be >= 0"
      sigma = thermal_speed
    else
      if (temperature_k < 0.0_dp) error stop "temperature_k must be >= 0"
      sigma = sqrt(k_boltzmann * temperature_k / m_particle)
    end if

    n = size(v, 2)
    allocate(z(3, n))
    call sample_standard_normal(z)

    v = sigma * z + spread(drift_velocity, dim=2, ncopies=n)
  end subroutine sample_shifted_maxwell_velocities

  !> 指定粒子数ぶんの位置/速度/電荷/質量/重みを生成し `particles_soa` を初期化する。
  !! @param[out] pcls 出力引数。
  !! @param[in] n 入力引数。
  !! @param[in] q_particle 入力引数。
  !! @param[in] m_particle 入力引数。
  !! @param[in] w_particle 入力引数。
  !! @param[in] pos_low 入力引数。
  !! @param[in] pos_high 入力引数。
  !! @param[in] drift_velocity 入力引数。
  !! @param[in]  temperature_k 入力引数。
  !! @param[in] thermal_speed 入力引数。
  subroutine init_random_beam_particles(pcls, n, q_particle, m_particle, w_particle, pos_low, pos_high, drift_velocity, &
                                        temperature_k, thermal_speed)
    type(particles_soa), intent(out) :: pcls
    integer(i32), intent(in) :: n
    real(dp), intent(in) :: q_particle, m_particle, w_particle
    real(dp), intent(in) :: pos_low(3), pos_high(3), drift_velocity(3)
    real(dp), intent(in), optional :: temperature_k, thermal_speed

    real(dp), allocatable :: x(:, :), v(:, :), q(:), m(:), w(:)

    if (n < 0) error stop "n must be non-negative"

    allocate(x(3, n), v(3, n), q(n), m(n), w(n))
    call sample_uniform_positions(pos_low, pos_high, x)
    call sample_shifted_maxwell_velocities(drift_velocity, m_particle, v, temperature_k, thermal_speed)
    q = q_particle
    m = m_particle
    w = w_particle

    call init_particles(pcls, x, v, q, m, w)
  end subroutine init_random_beam_particles

  !> Box–Muller法で標準正規乱数を生成し、任意形状配列へ詰める。
  !! @param[out] z 出力引数。
  subroutine sample_standard_normal(z)
    real(dp), intent(out) :: z(:, :)
    integer :: n_total, n_pair, i
    real(dp), allocatable :: out(:), u1(:), u2(:)
    real(dp) :: r, theta, pi

    n_total = size(z)
    n_pair = (n_total + 1) / 2
    pi = acos(-1.0_dp)

    allocate(out(n_total), u1(n_pair), u2(n_pair))
    call random_number(u1)
    call random_number(u2)

    i = 1
    do n_pair = 1, size(u1)
      if (u1(n_pair) <= tiny(1.0_dp)) u1(n_pair) = tiny(1.0_dp)
      r = sqrt(-2.0_dp * log(u1(n_pair)))
      theta = 2.0_dp * pi * u2(n_pair)
      out(i) = r * cos(theta)
      if (i + 1 <= n_total) out(i + 1) = r * sin(theta)
      i = i + 2
    end do

    z = reshape(out, shape(z))
  end subroutine sample_standard_normal

end module bem_injection

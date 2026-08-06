!> 粒子注入で共有する乱数初期化と基本分布sampling。
module bem_injection_random
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
  use bem_particles, only: init_particles
  use bem_types, only: particles_soa
  implicit none
  private

  real(dp), parameter :: default_velocity_sigma_cutoff = 6.0_dp

  public :: seed_rng
  public :: sample_uniform_positions
  public :: sample_shifted_maxwell_velocities
  public :: init_random_beam_particles
  public :: sample_standard_normal

contains

  !> 与えたシード列またはシステム時刻からFortran乱数生成器を初期化する。
  !! @param[in] seed 乱数生成器へ与えるシード列（省略時は `system_clock` 値から生成）。
  subroutine seed_rng(seed)
    integer(i32), intent(in), optional :: seed(:)
    integer :: n, i, clk
    integer, allocatable :: put(:)

    call random_seed(size=n)
    allocate (put(n))

    if (present(seed)) then
      do i = 1, n
        put(i) = seed(mod(i - 1, size(seed)) + 1) + 104729*i
      end do
    else
      call system_clock(count=clk)
      do i = 1, n
        put(i) = clk + 37*i
      end do
    end if

    call random_seed(put=put)
  end subroutine seed_rng

  !> 直方体領域 `[low, high]` 内で一様分布の初期位置をサンプリングする。
  !! @param[in] low サンプリング領域の下限座標 `(x,y,z)` [m]。
  !! @param[in] high サンプリング領域の上限座標 `(x,y,z)` [m]。
  !! @param[out] x 生成した粒子位置配列 `x(3,n)` [m]。
  subroutine sample_uniform_positions(low, high, x)
    real(dp), intent(in) :: low(3), high(3)
    real(dp), intent(out) :: x(:, :)
    real(dp), allocatable :: u(:, :)

    if (size(x, 1) /= 3) error stop "x first dimension must be 3"
    if (any(high < low)) error stop "high must be >= low for all axes"

    allocate (u(3, size(x, 2)))
    call random_number(u)
    x = spread(low, dim=2, ncopies=size(x, 2)) + spread(high - low, dim=2, ncopies=size(x, 2))*u
  end subroutine sample_uniform_positions

  !> ドリフト速度付きMaxwell分布(温度または熱速度指定)から粒子速度を生成する。
  !! @param[in] drift_velocity 付加する平均ドリフト速度ベクトル `(vx,vy,vz)` [m/s]。
  !! @param[in] m_particle 粒子1個あたりの質量 [kg]（`temperature_k` 指定時に使用）。
  !! @param[out] v サンプリングした速度配列 `v(3,n)` [m/s]。
  !! @param[in] temperature_k 熱運動の温度 [K]（`thermal_speed` 未指定時に使用）。
  !! @param[in] thermal_speed 熱速度の標準偏差 `sigma` [m/s]（指定時は温度より優先）。
  !! @param[in] sigma_cutoff 標準正規変量を `[-sigma_cutoff, sigma_cutoff]` に切る上限（省略時 6）。
  subroutine sample_shifted_maxwell_velocities(drift_velocity, m_particle, v, temperature_k, thermal_speed, sigma_cutoff)
    real(dp), intent(in) :: drift_velocity(3)
    real(dp), intent(in) :: m_particle
    real(dp), intent(out) :: v(:, :)
    real(dp), intent(in), optional :: temperature_k
    real(dp), intent(in), optional :: thermal_speed
    real(dp), intent(in), optional :: sigma_cutoff
    integer :: n
    real(dp) :: sigma, cutoff
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
      sigma = sqrt(k_boltzmann*temperature_k/m_particle)
    end if
    cutoff = default_velocity_sigma_cutoff
    if (present(sigma_cutoff)) cutoff = sigma_cutoff

    n = size(v, 2)
    allocate (z(3, n))
    call sample_standard_normal(z, sigma_cutoff=cutoff)

    v = sigma*z + spread(drift_velocity, dim=2, ncopies=n)
  end subroutine sample_shifted_maxwell_velocities

  !> 指定粒子数ぶんの位置/速度/電荷/質量/重みを生成し `particles_soa` を初期化する。
  !! @param[out] pcls 生成した粒子群を保持する `particles_soa`。
  !! @param[in] n 生成するマクロ粒子数。
  !! @param[in] q_particle 粒子1個あたりの電荷 [C]。
  !! @param[in] m_particle 粒子1個あたりの質量 [kg]。
  !! @param[in] w_particle 粒子1個あたりのマクロ粒子重み。
  !! @param[in] pos_low 位置サンプリング領域の下限座標 `(x,y,z)` [m]。
  !! @param[in] pos_high 位置サンプリング領域の上限座標 `(x,y,z)` [m]。
  !! @param[in] drift_velocity 平均ドリフト速度ベクトル `(vx,vy,vz)` [m/s]。
  !! @param[in] temperature_k 熱運動の温度 [K]（`thermal_speed` 未指定時に使用）。
  !! @param[in] thermal_speed 熱速度の標準偏差 `sigma` [m/s]（指定時は温度より優先）。
  subroutine init_random_beam_particles(pcls, n, q_particle, m_particle, w_particle, pos_low, pos_high, drift_velocity, &
                                        temperature_k, thermal_speed)
    type(particles_soa), intent(out) :: pcls
    integer(i32), intent(in) :: n
    real(dp), intent(in) :: q_particle, m_particle, w_particle
    real(dp), intent(in) :: pos_low(3), pos_high(3), drift_velocity(3)
    real(dp), intent(in), optional :: temperature_k, thermal_speed

    real(dp), allocatable :: x(:, :), v(:, :), q(:), m(:), w(:)

    if (n < 0) error stop "n must be non-negative"

    allocate (x(3, n), v(3, n), q(n), m(n), w(n))
    call sample_uniform_positions(pos_low, pos_high, x)
    call sample_shifted_maxwell_velocities(drift_velocity, m_particle, v, temperature_k, thermal_speed)
    q = q_particle
    m = m_particle
    w = w_particle

    call init_particles(pcls, x, v, q, m, w)
  end subroutine init_random_beam_particles

  !> Box–Muller法で標準正規乱数を生成し、任意形状配列へ詰める。
  !! @param[out] z 平均0・分散1の標準正規乱数で埋める出力配列。
  subroutine sample_standard_normal(z, sigma_cutoff)
    real(dp), intent(out) :: z(:, :)
    real(dp), intent(in), optional :: sigma_cutoff
    integer :: n_total, i
    real(dp), allocatable :: out(:)
    real(dp) :: r, theta, pi, u1, u2, z1, z2, cutoff

    n_total = size(z)
    pi = acos(-1.0_dp)
    cutoff = default_velocity_sigma_cutoff
    if (present(sigma_cutoff)) cutoff = sigma_cutoff
    if (cutoff <= 0.0_dp) error stop "sigma_cutoff must be > 0"

    allocate (out(n_total))
    i = 1
    do while (i <= n_total)
      call random_number(u1)
      call random_number(u2)
      if (u1 <= tiny(1.0_dp)) u1 = tiny(1.0_dp)
      r = sqrt(-2.0_dp*log(u1))
      theta = 2.0_dp*pi*u2
      z1 = r*cos(theta)
      z2 = r*sin(theta)
      if (abs(z1) <= cutoff) then
        out(i) = z1
        i = i + 1
      end if
      if (i <= n_total .and. abs(z2) <= cutoff) then
        out(i) = z2
        i = i + 1
      end if
    end do

    z = reshape(out, shape(z))
  end subroutine sample_standard_normal

end module bem_injection_random

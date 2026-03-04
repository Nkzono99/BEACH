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
  public :: compute_inflow_flux_from_drifting_maxwellian
  public :: compute_face_area_from_bounds
  public :: compute_macro_particles_for_batch
  public :: sample_reservoir_face_particles

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

  !> drifting Maxwellian の片側流入束 [#/m^2/s] を返す。
  real(dp) function compute_inflow_flux_from_drifting_maxwellian( &
    number_density_m3, temperature_k, m_particle, drift_velocity, inward_normal &
  ) result(gamma_in)
    real(dp), intent(in) :: number_density_m3
    real(dp), intent(in) :: temperature_k
    real(dp), intent(in) :: m_particle
    real(dp), intent(in) :: drift_velocity(3)
    real(dp), intent(in) :: inward_normal(3)
    real(dp) :: sigma, alpha, u_n

    if (number_density_m3 < 0.0_dp) error stop "number_density_m3 must be >= 0"
    if (temperature_k < 0.0_dp) error stop "temperature_k must be >= 0"
    if (m_particle <= 0.0_dp) error stop "m_particle must be > 0"

    u_n = dot_product(drift_velocity, inward_normal)
    sigma = sqrt(k_boltzmann * temperature_k / m_particle)
    if (sigma <= 0.0_dp) then
      gamma_in = number_density_m3 * max(0.0_dp, u_n)
      return
    end if

    alpha = u_n / sigma
    gamma_in = number_density_m3 * (sigma * standard_normal_pdf(alpha) + u_n * standard_normal_cdf(alpha))
  end function compute_inflow_flux_from_drifting_maxwellian

  !> 注入面上の矩形開口から有効面積[m^2]を返す。
  real(dp) function compute_face_area_from_bounds(inject_face, pos_low, pos_high) result(area)
    character(len=*), intent(in) :: inject_face
    real(dp), intent(in) :: pos_low(3), pos_high(3)
    integer :: axis_t1, axis_t2

    call resolve_face_axes(inject_face, axis_t1, axis_t2)
    area = (pos_high(axis_t1) - pos_low(axis_t1)) * (pos_high(axis_t2) - pos_low(axis_t2))
  end function compute_face_area_from_bounds

  !> 物理流量・重み・残差から今バッチのマクロ粒子数を決める。
  subroutine compute_macro_particles_for_batch( &
    number_density_m3, temperature_k, m_particle, drift_velocity, box_min, box_max, inject_face, pos_low, pos_high, &
    batch_duration, w_particle, residual, n_macro &
  )
    real(dp), intent(in) :: number_density_m3
    real(dp), intent(in) :: temperature_k
    real(dp), intent(in) :: m_particle
    real(dp), intent(in) :: drift_velocity(3)
    real(dp), intent(in) :: box_min(3), box_max(3)
    character(len=*), intent(in) :: inject_face
    real(dp), intent(in) :: pos_low(3), pos_high(3)
    real(dp), intent(in) :: batch_duration
    real(dp), intent(in) :: w_particle
    real(dp), intent(inout) :: residual
    integer(i32), intent(out) :: n_macro

    real(dp) :: inward_normal(3), gamma_in, area, n_phys_batch, n_macro_expected, macro_budget

    if (w_particle <= 0.0_dp) error stop "w_particle must be > 0"
    if (batch_duration < 0.0_dp) error stop "batch_duration must be >= 0"

    call resolve_face_geometry(box_min, box_max, inject_face, inward_normal=inward_normal)
    gamma_in = compute_inflow_flux_from_drifting_maxwellian( &
      number_density_m3, temperature_k, m_particle, drift_velocity, inward_normal &
    )
    area = compute_face_area_from_bounds(inject_face, pos_low, pos_high)
    n_phys_batch = gamma_in * area * batch_duration
    n_macro_expected = n_phys_batch / w_particle
    macro_budget = residual + n_macro_expected
    if (macro_budget < 0.0_dp) macro_budget = 0.0_dp
    if (macro_budget > real(huge(0_i32), dp)) error stop "macro particle count exceeds integer range"
    n_macro = int(floor(macro_budget), i32)
    residual = macro_budget - real(n_macro, dp)
  end subroutine compute_macro_particles_for_batch

  !> 上流リザーバ境界から流入する粒子群を面注入としてサンプルする。
  subroutine sample_reservoir_face_particles( &
    box_min, box_max, inject_face, pos_low, pos_high, drift_velocity, m_particle, temperature_k, batch_duration, x, v &
  )
    real(dp), intent(in) :: box_min(3), box_max(3)
    character(len=*), intent(in) :: inject_face
    real(dp), intent(in) :: pos_low(3), pos_high(3), drift_velocity(3)
    real(dp), intent(in) :: m_particle, temperature_k, batch_duration
    real(dp), intent(out) :: x(:, :)
    real(dp), intent(out) :: v(:, :)

    integer :: i, axis_n, axis_t1, axis_t2
    real(dp) :: boundary_value, inward_normal(3), sigma, u_n, vn
    real(dp), allocatable :: u(:, :), tau(:)

    if (size(x, 1) /= 3 .or. size(v, 1) /= 3) error stop "reservoir particle arrays must have first dimension 3"
    if (size(x, 2) /= size(v, 2)) error stop "reservoir x/v size mismatch"
    if (batch_duration < 0.0_dp) error stop "batch_duration must be >= 0"
    if (size(x, 2) == 0) return

    call resolve_face_geometry(box_min, box_max, inject_face, axis_n, boundary_value, inward_normal)
    call resolve_face_axes(inject_face, axis_t1, axis_t2)

    sigma = sqrt(k_boltzmann * temperature_k / m_particle)
    u_n = dot_product(drift_velocity, inward_normal)

    call sample_shifted_maxwell_velocities(drift_velocity, m_particle, v, temperature_k=temperature_k)
    call sample_flux_weighted_normal_component(u_n, sigma, v(axis_n, :))
    v(axis_n, :) = inward_normal(axis_n) * v(axis_n, :)

    allocate(u(2, size(x, 2)), tau(size(x, 2)))
    call random_number(u)
    call random_number(tau)

    do i = 1, size(x, 2)
      x(:, i) = 0.0_dp
      x(axis_n, i) = boundary_value
      x(axis_t1, i) = pos_low(axis_t1) + (pos_high(axis_t1) - pos_low(axis_t1)) * u(1, i)
      x(axis_t2, i) = pos_low(axis_t2) + (pos_high(axis_t2) - pos_low(axis_t2)) * u(2, i)
      x(:, i) = x(:, i) + v(:, i) * (tau(i) * batch_duration) + inward_normal * 1.0d-12
    end do
  end subroutine sample_reservoir_face_particles

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

  !> 標準正規分布の PDF を返す。
  pure real(dp) function standard_normal_pdf(x) result(pdf)
    real(dp), intent(in) :: x
    real(dp), parameter :: inv_sqrt_2pi = 3.98942280401432678d-1

    pdf = inv_sqrt_2pi * exp(-0.5_dp * x * x)
  end function standard_normal_pdf

  !> 標準正規分布の CDF を返す。
  pure real(dp) function standard_normal_cdf(x) result(cdf)
    real(dp), intent(in) :: x
    real(dp), parameter :: inv_sqrt_2 = 7.07106781186547524d-1

    cdf = 0.5_dp * (1.0_dp + erf(x * inv_sqrt_2))
  end function standard_normal_cdf

  !> 注入面名から接線2軸を返す。
  subroutine resolve_face_axes(inject_face, axis_t1, axis_t2)
    character(len=*), intent(in) :: inject_face
    integer, intent(out) :: axis_t1, axis_t2

    select case (trim(adjustl(inject_face)))
    case ('x_low', 'x_high')
      axis_t1 = 2
      axis_t2 = 3
    case ('y_low', 'y_high')
      axis_t1 = 3
      axis_t2 = 1
    case ('z_low', 'z_high')
      axis_t1 = 1
      axis_t2 = 2
    case default
      error stop "unknown inject_face"
    end select
  end subroutine resolve_face_axes

  !> 注入面名から法線方向の幾何情報を返す。
  subroutine resolve_face_geometry(box_min, box_max, inject_face, axis_n, boundary_value, inward_normal)
    real(dp), intent(in) :: box_min(3), box_max(3)
    character(len=*), intent(in) :: inject_face
    integer, intent(out), optional :: axis_n
    real(dp), intent(out), optional :: boundary_value
    real(dp), intent(out), optional :: inward_normal(3)

    integer :: axis_local
    real(dp) :: boundary_local, normal_local(3)

    normal_local = 0.0_dp
    select case (trim(adjustl(inject_face)))
    case ('x_low')
      axis_local = 1
      boundary_local = box_min(1)
      normal_local(1) = 1.0_dp
    case ('x_high')
      axis_local = 1
      boundary_local = box_max(1)
      normal_local(1) = -1.0_dp
    case ('y_low')
      axis_local = 2
      boundary_local = box_min(2)
      normal_local(2) = 1.0_dp
    case ('y_high')
      axis_local = 2
      boundary_local = box_max(2)
      normal_local(2) = -1.0_dp
    case ('z_low')
      axis_local = 3
      boundary_local = box_min(3)
      normal_local(3) = 1.0_dp
    case ('z_high')
      axis_local = 3
      boundary_local = box_max(3)
      normal_local(3) = -1.0_dp
    case default
      error stop "unknown inject_face"
    end select

    if (present(axis_n)) axis_n = axis_local
    if (present(boundary_value)) boundary_value = boundary_local
    if (present(inward_normal)) inward_normal = normal_local
  end subroutine resolve_face_geometry

  !> flux-weighted half-range 正規分布から法線速度をサンプルする。
  subroutine sample_flux_weighted_normal_component(mu, sigma, vn)
    real(dp), intent(in) :: mu, sigma
    real(dp), intent(out) :: vn(:)
    integer :: i
    real(dp) :: target, low, high, mid

    if (size(vn) == 0) return
    if (sigma <= 0.0_dp) then
      vn = max(mu, 0.0_dp)
      return
    end if

    do i = 1, size(vn)
      call random_number(target)
      low = 0.0_dp
      high = max(8.0_dp * sigma, mu + 8.0_dp * sigma)
      do while (flux_weighted_normal_cdf(high, mu, sigma) < target)
        high = high * 2.0_dp
      end do
      do while ((high - low) > max(1.0d-12, 1.0d-10 * (1.0d0 + high)))
        mid = 0.5_dp * (low + high)
        if (flux_weighted_normal_cdf(mid, mu, sigma) < target) then
          low = mid
        else
          high = mid
        end if
      end do
      vn(i) = 0.5_dp * (low + high)
    end do
  end subroutine sample_flux_weighted_normal_component

  !> flux-weighted half-range 正規分布の CDF を返す。
  pure real(dp) function flux_weighted_normal_cdf(vn, mu, sigma) result(cdf)
    real(dp), intent(in) :: vn, mu, sigma
    real(dp) :: alpha, x, denom, num

    if (vn <= 0.0_dp) then
      cdf = 0.0_dp
      return
    end if
    if (sigma <= 0.0_dp) then
      cdf = 1.0_dp
      return
    end if

    alpha = mu / sigma
    x = (vn - mu) / sigma
    denom = mu * standard_normal_cdf(alpha) + sigma * standard_normal_pdf(alpha)
    if (denom <= tiny(1.0_dp)) then
      cdf = 1.0_dp
      return
    end if

    num = mu * (standard_normal_cdf(x) - standard_normal_cdf(-alpha)) - &
          sigma * (standard_normal_pdf(x) - standard_normal_pdf(alpha))
    cdf = min(1.0_dp, max(0.0_dp, num / denom))
  end function flux_weighted_normal_cdf

end module bem_injection

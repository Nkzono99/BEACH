!> Maxwell reservoirの流束、macro粒子数、面速度sampling。
module bem_injection_flux
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
  use bem_injection_geometry, only: compute_face_area_from_bounds, resolve_face_axes, resolve_face_geometry
  use bem_injection_random, only: sample_shifted_maxwell_velocities
  implicit none
  private

  real(dp), parameter :: default_velocity_sigma_cutoff = 6.0_dp
  real(dp), parameter :: inv_sqrt_2 = 7.07106781186547524d-1

  public :: compute_inflow_flux_from_drifting_maxwellian
  public :: compute_macro_particles_for_batch
  public :: compute_macro_particles_from_flux
  public :: sample_reservoir_face_particles
  public :: sample_flux_weighted_normal_component

contains

  !> drifting Maxwellian の片側流入束 [#/m^2/s] を返す。
  !! @param[in] number_density_m3 粒子数密度 [1/m^3]。
  !! @param[in] temperature_k 温度 [K]。
  !! @param[in] m_particle 粒子1個あたりの質量 [kg]。
  !! @param[in] drift_velocity ドリフト速度ベクトル `(vx,vy,vz)` [m/s]。
  !! @param[in] inward_normal 注入面の内向き単位法線ベクトル。
  !! @param[in] vmin_normal 法線速度の下限 [m/s]（省略時は 0）。
  !! @return gamma_in 片側流入束 [1/m^2/s]。
  real(dp) function compute_inflow_flux_from_drifting_maxwellian( &
    number_density_m3, temperature_k, m_particle, drift_velocity, inward_normal, vmin_normal &
    ) result(gamma_in)
    real(dp), intent(in) :: number_density_m3
    real(dp), intent(in) :: temperature_k
    real(dp), intent(in) :: m_particle
    real(dp), intent(in) :: drift_velocity(3)
    real(dp), intent(in) :: inward_normal(3)
    real(dp), intent(in), optional :: vmin_normal
    real(dp) :: sigma, u_n, vmin

    if (number_density_m3 < 0.0_dp) error stop "number_density_m3 must be >= 0"
    if (temperature_k < 0.0_dp) error stop "temperature_k must be >= 0"
    if (m_particle <= 0.0_dp) error stop "m_particle must be > 0"

    vmin = 0.0_dp
    if (present(vmin_normal)) vmin = max(0.0_dp, vmin_normal)
    u_n = dot_product(drift_velocity, inward_normal)
    sigma = sqrt(k_boltzmann*temperature_k/m_particle)
    if (sigma <= 0.0_dp) then
      if (u_n < vmin) then
        gamma_in = 0.0_dp
      else
        gamma_in = number_density_m3*u_n
      end if
      return
    end if

    gamma_in = number_density_m3*flux_weighted_normal_tail(vmin, u_n, sigma)
  end function compute_inflow_flux_from_drifting_maxwellian

  !> 物理流量・重み・残差から今バッチのマクロ粒子数を決める。
  !! @param[in] number_density_m3 粒子数密度 [1/m^3]。
  !! @param[in] temperature_k 温度 [K]。
  !! @param[in] m_particle 粒子1個あたりの質量 [kg]。
  !! @param[in] drift_velocity ドリフト速度ベクトル `(vx,vy,vz)` [m/s]。
  !! @param[in] box_min シミュレーションボックス下限座標 `(x,y,z)` [m]。
  !! @param[in] box_max シミュレーションボックス上限座標 `(x,y,z)` [m]。
  !! @param[in] inject_face 注入面識別子（`x_low/x_high/y_low/y_high/z_low/z_high`）。
  !! @param[in] pos_low 注入口矩形の下限座標 `(x,y,z)` [m]。
  !! @param[in] pos_high 注入口矩形の上限座標 `(x,y,z)` [m]。
  !! @param[in] batch_duration 1バッチの物理時間長 [s]。
  !! @param[in] w_particle マクロ粒子重み。
  !! @param[inout] residual 前バッチから繰り越すマクロ粒子端数（呼び出し後に更新）。
  !! @param[out] n_macro 今バッチで生成するマクロ粒子数。
  !! @param[in] vmin_normal 法線速度の下限 [m/s]（省略時は 0）。
  subroutine compute_macro_particles_for_batch( &
    number_density_m3, temperature_k, m_particle, drift_velocity, box_min, box_max, inject_face, pos_low, pos_high, &
    batch_duration, w_particle, residual, n_macro, vmin_normal &
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
    real(dp), intent(in), optional :: vmin_normal

    real(dp) :: inward_normal(3), gamma_in, area, n_phys_batch, n_macro_expected, macro_budget

    if (w_particle <= 0.0_dp) error stop "w_particle must be > 0"
    if (batch_duration < 0.0_dp) error stop "batch_duration must be >= 0"

    call resolve_face_geometry(box_min, box_max, inject_face, inward_normal=inward_normal)
    gamma_in = compute_inflow_flux_from_drifting_maxwellian( &
               number_density_m3, temperature_k, m_particle, drift_velocity, inward_normal, vmin_normal=vmin_normal &
               )
    area = compute_face_area_from_bounds(inject_face, pos_low, pos_high)
    n_phys_batch = gamma_in*area*batch_duration
    n_macro_expected = n_phys_batch/w_particle
    macro_budget = residual + n_macro_expected
    if (macro_budget < 0.0_dp) macro_budget = 0.0_dp
    if (macro_budget > real(huge(0_i32), dp)) error stop "macro particle count exceeds integer range"
    n_macro = int(floor(macro_budget), i32)
    residual = macro_budget - real(n_macro, dp)
  end subroutine compute_macro_particles_for_batch

  !> 指定済み粒子 flux・重み・残差から今バッチのマクロ粒子数を決める。
  !! @param[in] particle_flux_m2_s 粒子数 flux [1/m^2/s]。
  !! @param[in] inject_face 注入面識別子。
  !! @param[in] pos_low 注入口矩形の下限座標 `(x,y,z)` [m]。
  !! @param[in] pos_high 注入口矩形の上限座標 `(x,y,z)` [m]。
  !! @param[in] batch_duration 1バッチの物理時間長 [s]。
  !! @param[in] w_particle マクロ粒子重み。
  !! @param[inout] residual 前バッチから繰り越すマクロ粒子端数。
  !! @param[out] n_macro 今バッチで生成するマクロ粒子数。
  subroutine compute_macro_particles_from_flux( &
    particle_flux_m2_s, inject_face, pos_low, pos_high, batch_duration, w_particle, residual, n_macro &
    )
    real(dp), intent(in) :: particle_flux_m2_s
    character(len=*), intent(in) :: inject_face
    real(dp), intent(in) :: pos_low(3), pos_high(3)
    real(dp), intent(in) :: batch_duration, w_particle
    real(dp), intent(inout) :: residual
    integer(i32), intent(out) :: n_macro

    real(dp) :: area, n_phys_batch, n_macro_expected, macro_budget

    if (.not. ieee_is_finite(particle_flux_m2_s) .or. particle_flux_m2_s < 0.0_dp) then
      error stop "particle_flux_m2_s must be finite and >= 0"
    end if
    if (w_particle <= 0.0_dp) error stop "w_particle must be > 0"
    if (batch_duration < 0.0_dp) error stop "batch_duration must be >= 0"

    area = compute_face_area_from_bounds(inject_face, pos_low, pos_high)
    n_phys_batch = particle_flux_m2_s*area*batch_duration
    n_macro_expected = n_phys_batch/w_particle
    macro_budget = residual + n_macro_expected
    if (macro_budget < 0.0_dp) macro_budget = 0.0_dp
    if (macro_budget > real(huge(0_i32), dp)) error stop "macro particle count exceeds integer range"
    n_macro = int(floor(macro_budget), i32)
    residual = macro_budget - real(n_macro, dp)
  end subroutine compute_macro_particles_from_flux

  !> 上流リザーバ境界から流入する粒子群を面注入としてサンプルする。
  !! @param[in] box_min シミュレーションボックス下限座標 `(x,y,z)` [m]。
  !! @param[in] box_max シミュレーションボックス上限座標 `(x,y,z)` [m]。
  !! @param[in] inject_face 注入面識別子（`x_low/x_high/y_low/y_high/z_low/z_high`）。
  !! @param[in] pos_low 注入口矩形の下限座標 `(x,y,z)` [m]。
  !! @param[in] pos_high 注入口矩形の上限座標 `(x,y,z)` [m]。
  !! @param[in] drift_velocity ドリフト速度ベクトル `(vx,vy,vz)` [m/s]。
  !! @param[in] m_particle 粒子1個あたりの質量 [kg]。
  !! @param[in] temperature_k 温度 [K]。
  !! @param[in] batch_duration 1バッチの物理時間長 [s]（現在は妥当性チェックのみ）。
  !! @param[out] x サンプリングした位置配列 `x(3,n)` [m]。
  !! @param[out] v サンプリングした速度配列 `v(3,n)` [m/s]。
  !! @param[in] barrier_normal_energy 法線方向のエネルギー障壁 `2 q Δφ / m` [`m^2/s^2`]（省略時 0）。
  !! @param[in] vmin_normal 法線速度の下限 [m/s]（省略時は `barrier_normal_energy` から自動導出）。
  !! @param[in] position_jitter_dt 互換用のlaunch位相時間幅[s]。乱数列は保つが位置移動には使わない。
  !! @param[in] apply_barrier_energy_shift `.true.` のとき法線速度へ障壁エネルギー変換を適用する。
  subroutine sample_reservoir_face_particles( &
    box_min, box_max, inject_face, pos_low, pos_high, drift_velocity, m_particle, temperature_k, batch_duration, x, v, &
    barrier_normal_energy, vmin_normal, position_jitter_dt, apply_barrier_energy_shift &
    )
    real(dp), intent(in) :: box_min(3), box_max(3)
    character(len=*), intent(in) :: inject_face
    real(dp), intent(in) :: pos_low(3), pos_high(3), drift_velocity(3)
    real(dp), intent(in) :: m_particle, temperature_k, batch_duration
    real(dp), intent(out) :: x(:, :)
    real(dp), intent(out) :: v(:, :)
    real(dp), intent(in), optional :: barrier_normal_energy
    real(dp), intent(in), optional :: vmin_normal
    real(dp), intent(in), optional :: position_jitter_dt
    logical, intent(in), optional :: apply_barrier_energy_shift

    integer :: i, axis_n, axis_t1, axis_t2
    real(dp) :: boundary_value, inward_normal(3), sigma, u_n, vn_floor, barrier, vn_inf, jitter_dt
    real(dp), allocatable :: u(:, :), tau(:)
    logical :: apply_energy_shift

    if (size(x, 1) /= 3 .or. size(v, 1) /= 3) error stop "reservoir particle arrays must have first dimension 3"
    if (size(x, 2) /= size(v, 2)) error stop "reservoir x/v size mismatch"
    if (batch_duration < 0.0_dp) error stop "batch_duration must be >= 0"
    if (present(position_jitter_dt)) then
      if (position_jitter_dt < 0.0_dp) error stop "position_jitter_dt must be >= 0"
      jitter_dt = position_jitter_dt
    else
      jitter_dt = 0.0_dp
    end if
    apply_energy_shift = .true.
    if (present(apply_barrier_energy_shift)) apply_energy_shift = apply_barrier_energy_shift
    if (size(x, 2) == 0) return

    call resolve_face_geometry(box_min, box_max, inject_face, axis_n, boundary_value, inward_normal)
    call resolve_face_axes(inject_face, axis_t1, axis_t2)

    sigma = sqrt(k_boltzmann*temperature_k/m_particle)
    u_n = dot_product(drift_velocity, inward_normal)
    barrier = 0.0_dp
    if (present(barrier_normal_energy)) barrier = barrier_normal_energy
    vn_floor = 0.0_dp
    if (barrier > 0.0_dp) vn_floor = sqrt(barrier)
    if (present(vmin_normal)) vn_floor = max(vn_floor, max(0.0_dp, vmin_normal))

    call sample_shifted_maxwell_velocities(drift_velocity, m_particle, v, temperature_k=temperature_k)
    call sample_flux_weighted_normal_component(u_n, sigma, v(axis_n, :), vmin_normal=vn_floor)
    do i = 1, size(v, 2)
      if (apply_energy_shift) then
        vn_inf = v(axis_n, i)
        vn_inf = sqrt(max(0.0_dp, vn_inf*vn_inf - barrier))
        v(axis_n, i) = inward_normal(axis_n)*vn_inf
      else
        v(axis_n, i) = inward_normal(axis_n)*v(axis_n, i)
      end if
    end do

    allocate (u(2, size(x, 2)))
    call random_number(u)
    if (jitter_dt > 0.0_dp) then
      allocate (tau(size(x, 2)))
      ! Keep the historical draw count so MPI ranks and restarted runs retain the same RNG stream.
      call random_number(tau)
    end if

    do i = 1, size(x, 2)
      x(:, i) = 0.0_dp
      ! One representable step inward avoids a zero-time boundary event without an untracked flight segment.
      x(axis_n, i) = nearest(boundary_value, inward_normal(axis_n))
      x(axis_t1, i) = pos_low(axis_t1) + (pos_high(axis_t1) - pos_low(axis_t1))*u(1, i)
      x(axis_t2, i) = pos_low(axis_t2) + (pos_high(axis_t2) - pos_low(axis_t2))*u(2, i)
    end do
  end subroutine sample_reservoir_face_particles

  !> 標準正規分布の PDF を返す。
  !! @param[in] x 評価点。
  !! @return pdf 標準正規分布の確率密度。
  pure real(dp) function standard_normal_pdf(x) result(pdf)
    real(dp), intent(in) :: x
    real(dp), parameter :: inv_sqrt_2pi = 3.98942280401432678d-1

    pdf = inv_sqrt_2pi*exp(-0.5_dp*x*x)
  end function standard_normal_pdf

  !> 標準正規分布の CDF を返す。
  !! @param[in] x 評価点。
  !! @return cdf 標準正規分布の累積分布値。
  pure real(dp) function standard_normal_cdf(x) result(cdf)
    real(dp), intent(in) :: x
    real(dp), parameter :: inv_sqrt_2 = 7.07106781186547524d-1

    cdf = 0.5_dp*(1.0_dp + erf(x*inv_sqrt_2))
  end function standard_normal_cdf

  !> flux-weighted half-range 正規分布から法線速度をサンプルする。
  !! @param[in] mu 法線速度分布の平均ドリフト成分 [m/s]。
  !! @param[in] sigma 法線速度分布の標準偏差 [m/s]。
  !! @param[in] vmin_normal 法線速度の下限 [m/s]（省略時は 0）。
  !! @param[out] vn サンプリングした法線速度配列 [m/s]（常に非負）。
  subroutine sample_flux_weighted_normal_component(mu, sigma, vn, vmin_normal, sigma_cutoff)
    real(dp), intent(in) :: mu, sigma
    real(dp), intent(in), optional :: vmin_normal
    real(dp), intent(in), optional :: sigma_cutoff
    real(dp), intent(out) :: vn(:)
    integer :: i
    real(dp) :: target, low, high, mid, vmin, cutoff, tail_min, target_tail, tail_high, width, tolerance

    if (size(vn) == 0) return
    vmin = 0.0_dp
    if (present(vmin_normal)) vmin = max(0.0_dp, vmin_normal)
    cutoff = default_velocity_sigma_cutoff
    if (present(sigma_cutoff)) cutoff = sigma_cutoff
    if (cutoff <= 0.0_dp) error stop "sigma_cutoff must be > 0"
    if (sigma <= 0.0_dp) then
      vn = max(mu, vmin)
      return
    end if
    tail_min = flux_weighted_normal_tail(vmin, mu, sigma)
    if (tail_min <= 0.0_dp) then
      vn = vmin
      return
    end if

    do i = 1, size(vn)
      call random_number(target)
      target_tail = tail_min*(1.0_dp - target)
      low = vmin
      width = max(sigma, spacing(max(abs(vmin), sigma)))
      high = vmin + width
      tail_high = flux_weighted_normal_tail(high, mu, sigma)
      do while (tail_high > target_tail)
        width = 2.0_dp*width
        high = vmin + width
        if (.not. ieee_is_finite(high)) error stop "flux-weighted normal inverse survival bracket overflow"
        tail_high = flux_weighted_normal_tail(high, mu, sigma)
      end do
      tolerance = sqrt(epsilon(1.0_dp))*max(sigma, abs(vmin), abs(high), tiny(1.0_dp))
      do while ((high - low) > tolerance)
        mid = 0.5_dp*(low + high)
        if (flux_weighted_normal_tail(mid, mu, sigma) > target_tail) then
          low = mid
        else
          high = mid
        end if
      end do
      vn(i) = 0.5_dp*(low + high)
    end do
  end subroutine sample_flux_weighted_normal_component

  !> flux-weighted half-range 正規分布の CDF を返す。
  !! @param[in] vn 評価する法線速度 [m/s]。
  !! @param[in] mu 法線速度分布の平均ドリフト成分 [m/s]。
  !! @param[in] sigma 法線速度分布の標準偏差 [m/s]。
  !! @param[in] vmin_normal 法線速度の下限 [m/s]（省略時は 0）。
  !! @return cdf flux-weighted half-range 正規分布の累積分布値。
  pure real(dp) function flux_weighted_normal_cdf(vn, mu, sigma, vmin_normal) result(cdf)
    real(dp), intent(in) :: vn, mu, sigma
    real(dp), intent(in), optional :: vmin_normal
    real(dp) :: vmin, denom, num

    vmin = 0.0_dp
    if (present(vmin_normal)) vmin = max(0.0_dp, vmin_normal)
    if (vn <= vmin) then
      cdf = 0.0_dp
      return
    end if
    if (sigma <= 0.0_dp) then
      cdf = 1.0_dp
      return
    end if

    denom = flux_weighted_normal_tail(vmin, mu, sigma)
    if (denom <= tiny(1.0_dp)) then
      cdf = 1.0_dp
      return
    end if

    num = denom - flux_weighted_normal_tail(vn, mu, sigma)
    cdf = min(1.0_dp, max(0.0_dp, num/denom))
  end function flux_weighted_normal_cdf

  !> flux-weighted 正規分布の tail 積分 `∫[vmin,∞] v f(v) dv` を返す。
  !! @param[in] vmin 法線速度の下限 [m/s]。
  !! @param[in] mu 法線速度分布の平均ドリフト成分 [m/s]。
  !! @param[in] sigma 法線速度分布の標準偏差 [m/s]。
  !! @return tail tail 積分値 [m/s]。
  pure real(dp) function flux_weighted_normal_tail(vmin, mu, sigma) result(tail)
    real(dp), intent(in) :: vmin, mu, sigma
    real(dp) :: x, pdf, survival, residual

    if (sigma <= 0.0_dp) then
      if (mu >= vmin) then
        tail = mu
      else
        tail = 0.0_dp
      end if
      return
    end if

    x = (vmin - mu)/sigma
    pdf = standard_normal_pdf(x)
    survival = 0.5_dp*erfc(x*inv_sqrt_2)
    if (x > 8.0_dp) then
      residual = pdf*normal_tail_residual_ratio(x)
    else
      residual = pdf - x*survival
    end if
    tail = vmin*survival + sigma*max(0.0_dp, residual)
  end function flux_weighted_normal_tail

  !> `phi(x) - x Q(x) = phi(x) h(x)` の large-x 漸近係数を返す。
  pure real(dp) function normal_tail_residual_ratio(x) result(ratio)
    real(dp), intent(in) :: x
    real(dp) :: inv_x2, term, candidate, previous_abs
    integer :: order

    inv_x2 = 1.0_dp/(x*x)
    term = inv_x2
    ratio = term
    previous_abs = abs(term)
    do order = 2, 64
      term = -term*real(2*order - 1, dp)*inv_x2
      if (abs(term) >= previous_abs) exit
      candidate = ratio + term
      if (candidate <= 0.0_dp) exit
      ratio = candidate
      previous_abs = abs(term)
      if (abs(term) <= epsilon(1.0_dp)*abs(ratio)) exit
    end do
  end function normal_tail_residual_ratio

end module bem_injection_flux

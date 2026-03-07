!> 乱数シード設定と粒子位置/速度サンプリングを担う粒子注入モジュール。
module bem_injection
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
  use bem_particles, only: init_particles
  use bem_types, only: particles_soa, mesh_type, sim_config, hit_info
  use bem_boundary, only: apply_box_boundary
  use bem_collision, only: find_first_hit
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
  public :: sample_photo_raycast_particles

contains

  !> 与えたシード列またはシステム時刻からFortran乱数生成器を初期化する。
  !! @param[in] seed 乱数生成器へ与えるシード列（省略時は `system_clock` 値から生成）。
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
  !! @param[in] low サンプリング領域の下限座標 `(x,y,z)` [m]。
  !! @param[in] high サンプリング領域の上限座標 `(x,y,z)` [m]。
  !! @param[out] x 生成した粒子位置配列 `x(3,n)` [m]。
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
  !! @param[in] drift_velocity 付加する平均ドリフト速度ベクトル `(vx,vy,vz)` [m/s]。
  !! @param[in] m_particle 粒子1個あたりの質量 [kg]（`temperature_k` 指定時に使用）。
  !! @param[out] v サンプリングした速度配列 `v(3,n)` [m/s]。
  !! @param[in] temperature_k 熱運動の温度 [K]（`thermal_speed` 未指定時に使用）。
  !! @param[in] thermal_speed 熱速度の標準偏差 `sigma` [m/s]（指定時は温度より優先）。
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

    allocate(x(3, n), v(3, n), q(n), m(n), w(n))
    call sample_uniform_positions(pos_low, pos_high, x)
    call sample_shifted_maxwell_velocities(drift_velocity, m_particle, v, temperature_k, thermal_speed)
    q = q_particle
    m = m_particle
    w = w_particle

    call init_particles(pcls, x, v, q, m, w)
  end subroutine init_random_beam_particles

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
    sigma = sqrt(k_boltzmann * temperature_k / m_particle)
    if (sigma <= 0.0_dp) then
      if (u_n < vmin) then
        gamma_in = 0.0_dp
      else
        gamma_in = number_density_m3 * u_n
      end if
      return
    end if

    gamma_in = number_density_m3 * flux_weighted_normal_tail(vmin, u_n, sigma)
  end function compute_inflow_flux_from_drifting_maxwellian

  !> 注入面上の矩形開口から有効面積[m^2]を返す。
  !! @param[in] inject_face 注入面識別子（`x_low/x_high/y_low/y_high/z_low/z_high`）。
  !! @param[in] pos_low 開口領域の下限座標 `(x,y,z)` [m]。
  !! @param[in] pos_high 開口領域の上限座標 `(x,y,z)` [m]。
  !! @return area 注入開口の有効面積 [m^2]。
  real(dp) function compute_face_area_from_bounds(inject_face, pos_low, pos_high) result(area)
    character(len=*), intent(in) :: inject_face
    real(dp), intent(in) :: pos_low(3), pos_high(3)
    integer :: axis_t1, axis_t2

    call resolve_face_axes(inject_face, axis_t1, axis_t2)
    area = (pos_high(axis_t1) - pos_low(axis_t1)) * (pos_high(axis_t2) - pos_low(axis_t2))
  end function compute_face_area_from_bounds

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
    n_phys_batch = gamma_in * area * batch_duration
    n_macro_expected = n_phys_batch / w_particle
    macro_budget = residual + n_macro_expected
    if (macro_budget < 0.0_dp) macro_budget = 0.0_dp
    if (macro_budget > real(huge(0_i32), dp)) error stop "macro particle count exceeds integer range"
    n_macro = int(floor(macro_budget), i32)
    residual = macro_budget - real(n_macro, dp)
  end subroutine compute_macro_particles_for_batch

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
  !! @param[in] position_jitter_dt 初期位置に速度方向で与えるランダムジッタ時間幅[s]（省略時は 0）。
  subroutine sample_reservoir_face_particles( &
    box_min, box_max, inject_face, pos_low, pos_high, drift_velocity, m_particle, temperature_k, batch_duration, x, v, &
    barrier_normal_energy, vmin_normal, position_jitter_dt &
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

    integer :: i, axis_n, axis_t1, axis_t2
    real(dp) :: boundary_value, inward_normal(3), sigma, u_n, vn_floor, barrier, vn_inf, jitter_dt
    real(dp), allocatable :: u(:, :), tau(:)

    if (size(x, 1) /= 3 .or. size(v, 1) /= 3) error stop "reservoir particle arrays must have first dimension 3"
    if (size(x, 2) /= size(v, 2)) error stop "reservoir x/v size mismatch"
    if (batch_duration < 0.0_dp) error stop "batch_duration must be >= 0"
    if (present(position_jitter_dt)) then
      if (position_jitter_dt < 0.0_dp) error stop "position_jitter_dt must be >= 0"
      jitter_dt = position_jitter_dt
    else
      jitter_dt = 0.0_dp
    end if
    if (size(x, 2) == 0) return

    call resolve_face_geometry(box_min, box_max, inject_face, axis_n, boundary_value, inward_normal)
    call resolve_face_axes(inject_face, axis_t1, axis_t2)

    sigma = sqrt(k_boltzmann * temperature_k / m_particle)
    u_n = dot_product(drift_velocity, inward_normal)
    barrier = 0.0_dp
    if (present(barrier_normal_energy)) barrier = barrier_normal_energy
    vn_floor = 0.0_dp
    if (barrier > 0.0_dp) vn_floor = sqrt(barrier)
    if (present(vmin_normal)) vn_floor = max(vn_floor, max(0.0_dp, vmin_normal))

    call sample_shifted_maxwell_velocities(drift_velocity, m_particle, v, temperature_k=temperature_k)
    call sample_flux_weighted_normal_component(u_n, sigma, v(axis_n, :), vmin_normal=vn_floor)
    do i = 1, size(v, 2)
      vn_inf = v(axis_n, i)
      vn_inf = sqrt(max(0.0_dp, vn_inf * vn_inf - barrier))
      v(axis_n, i) = inward_normal(axis_n) * vn_inf
    end do

    allocate(u(2, size(x, 2)))
    call random_number(u)
    if (jitter_dt > 0.0_dp) then
      allocate(tau(size(x, 2)))
      call random_number(tau)
    end if

    do i = 1, size(x, 2)
      x(:, i) = 0.0_dp
      x(axis_n, i) = boundary_value
      x(axis_t1, i) = pos_low(axis_t1) + (pos_high(axis_t1) - pos_low(axis_t1)) * u(1, i)
      x(axis_t2, i) = pos_low(axis_t2) + (pos_high(axis_t2) - pos_low(axis_t2)) * u(2, i)
      if (jitter_dt > 0.0_dp) x(:, i) = x(:, i) + v(:, i) * (tau(i) * jitter_dt)
      x(:, i) = x(:, i) + inward_normal * 1.0d-12
    end do
  end subroutine sample_reservoir_face_particles

  !> 光線を注入面からレイキャストし、最初の命中要素から光電子を放出する。
  !! @param[in] mesh 交差判定に使うメッシュ。
  !! @param[in] sim ボックス境界条件とバッチ時間を含むシミュレーション設定。
  !! @param[in] inject_face 照射面識別子（`x_low/x_high/y_low/y_high/z_low/z_high`）。
  !! @param[in] pos_low 照射開口の下限座標 `(x,y,z)` [m]。
  !! @param[in] pos_high 照射開口の上限座標 `(x,y,z)` [m]。
  !! @param[in] ray_direction レイ進行方向ベクトル（正規化前でも可）。
  !! @param[in] m_particle 粒子1個あたりの質量 [kg]。
  !! @param[in] temperature_k 放出温度 [K]。
  !! @param[in] normal_drift_speed 放出法線方向のシフト速度 [m/s]。
  !! @param[in] emit_current_density_a_m2 レイ垂直面基準の放出電流面密度 [A/m^2]。
  !! @param[in] q_particle 粒子1個あたりの電荷 [C]。
  !! @param[in] rays_per_batch このrankで発射するレイ本数。
  !! @param[in] global_rays_per_batch 全rank合計のレイ本数（省略時は `rays_per_batch`）。
  !! @param[out] x 放出位置配列 `x(3,rays_per_batch)` [m]。
  !! @param[out] v 放出速度配列 `v(3,rays_per_batch)` [m/s]。
  !! @param[out] w 各マクロ粒子重み `w(rays_per_batch)`。
  !! @param[out] n_emit 実際に放出された粒子数（`<= rays_per_batch`）。
  !! @param[out] emit_elem_idx 放出元要素ID `emit_elem_idx(rays_per_batch)`（省略可）。
  subroutine sample_photo_raycast_particles( &
    mesh, sim, inject_face, pos_low, pos_high, ray_direction, m_particle, temperature_k, normal_drift_speed, &
    emit_current_density_a_m2, q_particle, rays_per_batch, x, v, w, n_emit, emit_elem_idx, global_rays_per_batch &
  )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    character(len=*), intent(in) :: inject_face
    real(dp), intent(in) :: pos_low(3), pos_high(3)
    real(dp), intent(in) :: ray_direction(3)
    real(dp), intent(in) :: m_particle, temperature_k, normal_drift_speed
    real(dp), intent(in) :: emit_current_density_a_m2, q_particle
    integer(i32), intent(in) :: rays_per_batch
    real(dp), intent(out) :: x(:, :)
    real(dp), intent(out) :: v(:, :)
    real(dp), intent(out) :: w(:)
    integer(i32), intent(out) :: n_emit
    integer(i32), intent(out), optional :: emit_elem_idx(:)
    integer(i32), intent(in), optional :: global_rays_per_batch

    real(dp), parameter :: eps = 1.0d-12
    integer(i32) :: i, total_rays
    integer(i32) :: bounce_count
    integer :: axis_n, axis_t1, axis_t2
    real(dp) :: boundary_value, inward_normal(3), launch_dir(3), launch_dir_norm, inward_dot
    real(dp) :: launch_area, projected_area, w_hit, sigma
    real(dp) :: ray_pos(3), ray_dir(3), seg_end(3), boundary_probe(3), boundary_dir(3)
    real(dp) :: surf_normal(3), tangent1(3), tangent2(3)
    real(dp), allocatable :: u(:, :)
    logical :: reached_boundary, alive, escaped_boundary
    type(hit_info) :: hit

    if (size(x, 1) /= 3 .or. size(v, 1) /= 3) error stop "photo_raycast particle arrays must have first dimension 3"
    if (size(x, 2) < rays_per_batch .or. size(v, 2) < rays_per_batch) then
      error stop "photo_raycast x/v arrays are smaller than rays_per_batch"
    end if
    if (size(w) < rays_per_batch) error stop "photo_raycast w array is smaller than rays_per_batch"
    if (present(emit_elem_idx)) then
      if (size(emit_elem_idx) < rays_per_batch) error stop "photo_raycast emit_elem_idx is smaller than rays_per_batch"
      emit_elem_idx = -1_i32
    end if
    if (rays_per_batch <= 0_i32) error stop "rays_per_batch must be > 0"
    total_rays = rays_per_batch
    if (present(global_rays_per_batch)) total_rays = global_rays_per_batch
    if (total_rays <= 0_i32) error stop "global_rays_per_batch must be > 0"
    if (.not. sim%use_box) error stop "photo_raycast requires sim.use_box = true"
    if (sim%batch_duration <= 0.0_dp) error stop "photo_raycast requires sim.batch_duration > 0"
    if (m_particle <= 0.0_dp) error stop "m_particle must be > 0"
    if (temperature_k < 0.0_dp) error stop "temperature_k must be >= 0"
    if (emit_current_density_a_m2 <= 0.0_dp) error stop "emit_current_density_a_m2 must be > 0"
    if (abs(q_particle) <= 0.0_dp) error stop "q_particle must be non-zero"

    call resolve_face_geometry(sim%box_min, sim%box_max, inject_face, axis_n, boundary_value, inward_normal)
    call resolve_face_axes(inject_face, axis_t1, axis_t2)

    launch_dir = ray_direction
    launch_dir_norm = sqrt(sum(launch_dir * launch_dir))
    if (launch_dir_norm <= 0.0_dp) error stop "ray_direction norm must be > 0"
    launch_dir = launch_dir / launch_dir_norm
    inward_dot = dot_product(launch_dir, inward_normal)
    if (inward_dot <= 0.0_dp) error stop "ray_direction must point inward from inject_face"

    launch_area = compute_face_area_from_bounds(inject_face, pos_low, pos_high)
    if (launch_area <= 0.0_dp) error stop "photo_raycast opening area must be positive"
    projected_area = launch_area * abs(inward_dot)
    w_hit = emit_current_density_a_m2 * projected_area * sim%batch_duration / (abs(q_particle) * real(total_rays, dp))
    if (w_hit <= 0.0_dp) error stop "photo_raycast produced invalid w_hit"
    sigma = sqrt(k_boltzmann * temperature_k / m_particle)

    n_emit = 0_i32
    x = 0.0_dp
    v = 0.0_dp
    w = 0.0_dp

    allocate(u(2, rays_per_batch))
    call random_number(u)
    do i = 1_i32, rays_per_batch
      ray_pos = 0.0_dp
      ray_pos(axis_n) = boundary_value
      ray_pos(axis_t1) = pos_low(axis_t1) + (pos_high(axis_t1) - pos_low(axis_t1)) * u(1, i)
      ray_pos(axis_t2) = pos_low(axis_t2) + (pos_high(axis_t2) - pos_low(axis_t2)) * u(2, i)
      ray_dir = launch_dir
      ray_pos = ray_pos + ray_dir * eps

      alive = .true.
      bounce_count = 0_i32
      do while (alive .and. bounce_count <= sim%raycast_max_bounce)
        call step_ray_to_boundary(sim%box_min, sim%box_max, ray_pos, ray_dir, seg_end, reached_boundary)
        if (.not. reached_boundary) exit

        call find_first_hit(mesh, ray_pos, seg_end, hit)
        if (hit%has_hit) then
          if (n_emit >= int(size(w), i32)) error stop "photo_raycast emitted particle buffer overflow"
          n_emit = n_emit + 1_i32
          surf_normal = mesh%normals(:, hit%elem_idx)
          if (dot_product(surf_normal, ray_dir) > 0.0_dp) surf_normal = -surf_normal
          call build_tangent_basis(surf_normal, tangent1, tangent2)
          call sample_photo_emission_velocity(sigma, normal_drift_speed, surf_normal, tangent1, tangent2, v(:, n_emit))
          x(:, n_emit) = hit%pos + surf_normal * eps
          w(n_emit) = w_hit
          if (present(emit_elem_idx)) emit_elem_idx(n_emit) = hit%elem_idx
          exit
        end if

        boundary_probe = seg_end + ray_dir * eps
        boundary_dir = ray_dir
        escaped_boundary = .false.
        call apply_box_boundary(sim, boundary_probe, boundary_dir, alive, escaped_boundary)
        if (.not. alive) exit
        ray_dir = boundary_dir / sqrt(sum(boundary_dir * boundary_dir))
        ray_pos = boundary_probe + ray_dir * eps
        bounce_count = bounce_count + 1_i32
      end do
    end do
  end subroutine sample_photo_raycast_particles

  !> レイを現在位置から最初のボックス境界まで進める。
  !! @param[in] box_min ボックス下限座標 `(x,y,z)` [m]。
  !! @param[in] box_max ボックス上限座標 `(x,y,z)` [m]。
  !! @param[in] x0 レイの現在位置 [m]。
  !! @param[in] ray_dir レイ進行方向（単位ベクトル）。
  !! @param[out] x1 境界到達位置 [m]。
  !! @param[out] reached_boundary 境界到達位置が求まった場合 `.true.`。
  subroutine step_ray_to_boundary(box_min, box_max, x0, ray_dir, x1, reached_boundary)
    real(dp), intent(in) :: box_min(3), box_max(3)
    real(dp), intent(in) :: x0(3), ray_dir(3)
    real(dp), intent(out) :: x1(3)
    logical, intent(out) :: reached_boundary

    real(dp), parameter :: eps = 1.0d-14
    integer :: axis
    real(dp) :: t_axis, t_hit

    t_hit = huge(1.0_dp)
    do axis = 1, 3
      if (ray_dir(axis) > eps) then
        t_axis = (box_max(axis) - x0(axis)) / ray_dir(axis)
      else if (ray_dir(axis) < -eps) then
        t_axis = (box_min(axis) - x0(axis)) / ray_dir(axis)
      else
        cycle
      end if
      if (t_axis > eps .and. t_axis < t_hit) t_hit = t_axis
    end do

    if (t_hit >= huge(1.0_dp) * 0.5_dp) then
      reached_boundary = .false.
      x1 = x0
      return
    end if

    reached_boundary = .true.
    x1 = x0 + ray_dir * t_hit
    x1 = min(box_max, max(box_min, x1))
  end subroutine step_ray_to_boundary

  !> 面法線ベクトルから接線2軸を構築する。
  !! @param[in] normal 法線ベクトル。
  !! @param[out] tangent1 第1接線ベクトル。
  !! @param[out] tangent2 第2接線ベクトル。
  subroutine build_tangent_basis(normal, tangent1, tangent2)
    real(dp), intent(in) :: normal(3)
    real(dp), intent(out) :: tangent1(3), tangent2(3)

    real(dp) :: n(3), ref(3), norm_n, norm_t1

    norm_n = sqrt(sum(normal * normal))
    if (norm_n <= 0.0_dp) error stop "surface normal norm must be > 0"
    n = normal / norm_n

    if (abs(n(1)) < 0.9_dp) then
      ref = [1.0_dp, 0.0_dp, 0.0_dp]
    else
      ref = [0.0_dp, 1.0_dp, 0.0_dp]
    end if

    tangent1 = cross3(n, ref)
    norm_t1 = sqrt(sum(tangent1 * tangent1))
    if (norm_t1 <= 0.0_dp) error stop "failed to build tangent basis"
    tangent1 = tangent1 / norm_t1
    tangent2 = cross3(n, tangent1)
  end subroutine build_tangent_basis

  !> 光電子放出速度を局所法線座標でサンプルする。
  !! @param[in] sigma 熱速度標準偏差 [m/s]。
  !! @param[in] normal_drift_speed 放出法線方向のシフト速度 [m/s]。
  !! @param[in] normal 放出法線ベクトル（単位化済み）。
  !! @param[in] tangent1 第1接線ベクトル（単位化済み）。
  !! @param[in] tangent2 第2接線ベクトル（単位化済み）。
  !! @param[out] vel サンプルした速度ベクトル [m/s]。
  subroutine sample_photo_emission_velocity(sigma, normal_drift_speed, normal, tangent1, tangent2, vel)
    real(dp), intent(in) :: sigma, normal_drift_speed
    real(dp), intent(in) :: normal(3), tangent1(3), tangent2(3)
    real(dp), intent(out) :: vel(3)

    real(dp) :: vn(1), z(2, 1), vt1, vt2

    call sample_flux_weighted_normal_component(normal_drift_speed, sigma, vn)
    vt1 = 0.0_dp
    vt2 = 0.0_dp
    if (sigma > 0.0_dp) then
      call sample_standard_normal(z)
      vt1 = sigma * z(1, 1)
      vt2 = sigma * z(2, 1)
    end if
    vel = normal * vn(1) + tangent1 * vt1 + tangent2 * vt2
  end subroutine sample_photo_emission_velocity

  !> 3次元外積を返す。
  pure function cross3(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
  end function cross3

  !> Box–Muller法で標準正規乱数を生成し、任意形状配列へ詰める。
  !! @param[out] z 平均0・分散1の標準正規乱数で埋める出力配列。
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
  !! @param[in] x 評価点。
  !! @return pdf 標準正規分布の確率密度。
  pure real(dp) function standard_normal_pdf(x) result(pdf)
    real(dp), intent(in) :: x
    real(dp), parameter :: inv_sqrt_2pi = 3.98942280401432678d-1

    pdf = inv_sqrt_2pi * exp(-0.5_dp * x * x)
  end function standard_normal_pdf

  !> 標準正規分布の CDF を返す。
  !! @param[in] x 評価点。
  !! @return cdf 標準正規分布の累積分布値。
  pure real(dp) function standard_normal_cdf(x) result(cdf)
    real(dp), intent(in) :: x
    real(dp), parameter :: inv_sqrt_2 = 7.07106781186547524d-1

    cdf = 0.5_dp * (1.0_dp + erf(x * inv_sqrt_2))
  end function standard_normal_cdf

  !> 注入面名から接線2軸を返す。
  !! @param[in] inject_face 注入面識別子（`x_low/x_high/y_low/y_high/z_low/z_high`）。
  !! @param[out] axis_t1 注入面の第1接線軸インデックス（1:x, 2:y, 3:z）。
  !! @param[out] axis_t2 注入面の第2接線軸インデックス（1:x, 2:y, 3:z）。
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
  !! @param[in] box_min シミュレーションボックス下限座標 `(x,y,z)` [m]。
  !! @param[in] box_max シミュレーションボックス上限座標 `(x,y,z)` [m]。
  !! @param[in] inject_face 注入面識別子（`x_low/x_high/y_low/y_high/z_low/z_high`）。
  !! @param[out] axis_n 法線方向軸インデックス（省略可）。
  !! @param[out] boundary_value 注入面の境界座標値 [m]（省略可）。
  !! @param[out] inward_normal 注入面の内向き法線ベクトル（省略可）。
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
  !! @param[in] mu 法線速度分布の平均ドリフト成分 [m/s]。
  !! @param[in] sigma 法線速度分布の標準偏差 [m/s]。
  !! @param[in] vmin_normal 法線速度の下限 [m/s]（省略時は 0）。
  !! @param[out] vn サンプリングした法線速度配列 [m/s]（常に非負）。
  subroutine sample_flux_weighted_normal_component(mu, sigma, vn, vmin_normal)
    real(dp), intent(in) :: mu, sigma
    real(dp), intent(in), optional :: vmin_normal
    real(dp), intent(out) :: vn(:)
    integer :: i
    real(dp) :: target, low, high, mid, vmin

    if (size(vn) == 0) return
    vmin = 0.0_dp
    if (present(vmin_normal)) vmin = max(0.0_dp, vmin_normal)
    if (sigma <= 0.0_dp) then
      vn = max(mu, vmin)
      return
    end if

    do i = 1, size(vn)
      call random_number(target)
      low = vmin
      high = max(vmin + 8.0_dp * sigma, mu + 8.0_dp * sigma)
      do while (flux_weighted_normal_cdf(high, mu, sigma, vmin_normal=vmin) < target)
        high = high * 2.0_dp
      end do
      do while ((high - low) > max(1.0d-12, 1.0d-10 * (1.0d0 + high)))
        mid = 0.5_dp * (low + high)
        if (flux_weighted_normal_cdf(mid, mu, sigma, vmin_normal=vmin) < target) then
          low = mid
        else
          high = mid
        end if
      end do
      vn(i) = 0.5_dp * (low + high)
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
    cdf = min(1.0_dp, max(0.0_dp, num / denom))
  end function flux_weighted_normal_cdf

  !> flux-weighted 正規分布の tail 積分 `∫[vmin,∞] v f(v) dv` を返す。
  !! @param[in] vmin 法線速度の下限 [m/s]。
  !! @param[in] mu 法線速度分布の平均ドリフト成分 [m/s]。
  !! @param[in] sigma 法線速度分布の標準偏差 [m/s]。
  !! @return tail tail 積分値 [m/s]。
  pure real(dp) function flux_weighted_normal_tail(vmin, mu, sigma) result(tail)
    real(dp), intent(in) :: vmin, mu, sigma
    real(dp) :: x

    if (sigma <= 0.0_dp) then
      if (mu >= vmin) then
        tail = mu
      else
        tail = 0.0_dp
      end if
      return
    end if

    x = (vmin - mu) / sigma
    tail = mu * (1.0_dp - standard_normal_cdf(x)) + sigma * standard_normal_pdf(x)
    if (tail < 0.0_dp) tail = 0.0_dp
  end function flux_weighted_normal_tail

end module bem_injection

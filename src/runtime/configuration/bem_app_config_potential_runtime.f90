!> 粒子注入面とbox基準面の電位統計を評価する実行時変換。
module bem_app_config_potential_runtime
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
  use bem_types, only: mesh_type, sim_config
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_app_config_types, only: particle_species_spec
  use bem_string_utils, only: lower_ascii
  use bem_config_helpers, only: resolve_inward_normal
  use, intrinsic :: iso_fortran_env, only: error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private :: compute_rectangular_face_potential_statistics
  private :: species_temperature_k_local

contains

  !> reservoir_face 開口面の平均電位を `N x N` 格子平均で評価する。
  !! @param[in] mesh 現在バッチ開始時点の電荷分布メッシュ。
  !! @param[in] sim シミュレーション設定。
  !! @param[in] spec reservoir_face 粒子種設定。
  !! @param[inout] snapshot refresh 済み静電 snapshot。
  !! @param[out] phi_face 注入開口面の平均電位 [V]。
  !! @param[out] phi_std 評価格子上の電位の母標準偏差 [V]（省略可）。
  !! @param[out] phi_min 評価格子上の最小電位 [V]（省略可）。
  !! @param[out] phi_max 評価格子上の最大電位 [V]（省略可）。
  subroutine compute_face_average_potential(mesh, sim, spec, snapshot, phi_face, phi_std, phi_min, phi_max)
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(particle_species_spec), intent(in) :: spec
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(out) :: phi_face
    real(dp), intent(out), optional :: phi_std, phi_min, phi_max

    real(dp) :: phi_std_local, phi_min_local, phi_max_local

    call compute_rectangular_face_potential_statistics( &
      mesh, sim, snapshot, spec%inject_face, spec%pos_low, spec%pos_high, &
      sim%injection_face_phi_grid_n, phi_face, phi_std_local, phi_min_local, phi_max_local &
      )
    if (present(phi_std)) phi_std = phi_std_local
    if (present(phi_min)) phi_min = phi_min_local
    if (present(phi_max)) phi_max = phi_max_local
  end subroutine compute_face_average_potential

  !> box の z-high 全断面における電位統計をセル中心 `N x N` 格子で評価する。
  !!
  !! `N` は reservoir face 平均と同じ `sim.injection_face_phi_grid_n` を使う。
  !! 注入speciesの開口には依存せず、周期セル全断面を基準面として評価する。
  subroutine compute_z_high_box_potential_statistics( &
    mesh, sim, snapshot, phi_mean, phi_std, phi_min, phi_max &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(out) :: phi_mean, phi_std, phi_min, phi_max

    if (.not. sim%use_box) then
      error stop 'z-high box potential statistics require sim.use_box=true.'
    end if
    call compute_rectangular_face_potential_statistics( &
      mesh, sim, snapshot, 'z_high', sim%box_min, sim%box_max, &
      sim%injection_face_phi_grid_n, phi_mean, phi_std, phi_min, phi_max &
      )
  end subroutine compute_z_high_box_potential_statistics

  !> 指定box面の矩形範囲における電位統計をセル中心格子で評価する。
  subroutine compute_rectangular_face_potential_statistics( &
    mesh, sim, snapshot, face, pos_low, pos_high, ngrid, phi_mean, phi_std, phi_min, phi_max &
    )
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    character(len=*), intent(in) :: face
    real(dp), intent(in) :: pos_low(3), pos_high(3)
    integer(i32), intent(in) :: ngrid
    real(dp), intent(out) :: phi_mean, phi_std, phi_min, phi_max

    integer(i32) :: i, j, sample_count
    integer :: axis_n, axis_t1, axis_t2
    real(dp) :: boundary_value, inward_normal(3), pos(3), t1, t2, phi
    real(dp) :: phi_m2, delta

    call resolve_face_sampling_geometry( &
      sim%box_min, sim%box_max, face, axis_n, axis_t1, axis_t2, boundary_value, inward_normal &
      )

    if (ngrid < 1_i32) error stop 'face potential statistics require sample_n >= 1.'
    sample_count = 0_i32
    phi_mean = 0.0_dp
    phi_m2 = 0.0_dp
    phi_min = huge(1.0_dp)
    phi_max = -huge(1.0_dp)
    do i = 1_i32, ngrid
      t1 = (real(i, dp) - 0.5d0)/real(ngrid, dp)
      do j = 1_i32, ngrid
        t2 = (real(j, dp) - 0.5d0)/real(ngrid, dp)
        pos = 0.0d0
        pos(axis_n) = boundary_value
        pos(axis_t1) = pos_low(axis_t1) + (pos_high(axis_t1) - pos_low(axis_t1))*t1
        pos(axis_t2) = pos_low(axis_t2) + (pos_high(axis_t2) - pos_low(axis_t2))*t2
        pos = pos + inward_normal*1.0d-12
        call snapshot%eval_local_phi(mesh, sim, pos, phi)
        if (.not. ieee_is_finite(phi)) error stop 'face potential sample is non-finite.'
        sample_count = sample_count + 1_i32
        delta = phi - phi_mean
        phi_mean = phi_mean + delta/real(sample_count, dp)
        phi_m2 = phi_m2 + delta*(phi - phi_mean)
        phi_min = min(phi_min, phi)
        phi_max = max(phi_max, phi)
      end do
    end do

    phi_std = sqrt(max(phi_m2/real(sample_count, dp), 0.0_dp))
  end subroutine compute_rectangular_face_potential_statistics

  !> 注入面の局所電位差が面平均 reservoir 近似の特徴エネルギーに対して大きい場合に警告する。
  subroutine warn_face_average_potential_variation(sim, spec, phi_mean, phi_std, phi_min, phi_max)
    real(dp), parameter :: variation_warn_ratio = 0.1_dp
    type(sim_config), intent(in) :: sim
    type(particle_species_spec), intent(in) :: spec
    real(dp), intent(in) :: phi_mean, phi_std, phi_min, phi_max

    real(dp) :: inward_normal(3), normal_drift, characteristic_energy, variation_energy, ratio, phi_scale

    if (trim(lower_ascii(spec%velocity_distribution)) /= 'maxwellian') return
    phi_scale = max(1.0_dp, abs(phi_mean), abs(phi_min), abs(phi_max))
    if (phi_std <= 256.0_dp*epsilon(1.0_dp)*phi_scale) return

    call resolve_inward_normal(spec%inject_face, inward_normal)
    normal_drift = dot_product(spec%drift_velocity, inward_normal)
    characteristic_energy = k_boltzmann*species_temperature_k_local(spec) + &
      0.5_dp*spec%m_particle*normal_drift*normal_drift
    variation_energy = abs(spec%q_particle)*phi_std
    if (characteristic_energy > 0.0_dp) then
      ratio = variation_energy/characteristic_energy
      if (ratio <= variation_warn_ratio) return
    else
      ratio = huge(1.0_dp)
    end if

    write (error_unit, '(a,a,a,a,a,i0,a,es12.4,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
      'WARNING: reservoir face-average potential may be inaccurate: species=', trim(spec%species_key), &
      ' face=', trim(spec%inject_face), ' samples=', sim%injection_face_phi_grid_n**2, &
      ' mean_V=', phi_mean, ' std_V=', phi_std, ' min_V=', phi_min, ' max_V=', phi_max, &
      ' energy_ratio=', ratio
    flush (error_unit)
  end subroutine warn_face_average_potential_variation

  !> 注入面名から法線軸・接線軸・境界値・内向き法線を返す。
  !! @param[in] box_min シミュレーションボックス下限座標 `(x,y,z)` [m]。
  !! @param[in] box_max シミュレーションボックス上限座標 `(x,y,z)` [m]。
  !! @param[in] inject_face 注入面識別子。
  !! @param[out] axis_n 法線軸インデックス（1:x, 2:y, 3:z）。
  !! @param[out] axis_t1 第1接線軸インデックス。
  !! @param[out] axis_t2 第2接線軸インデックス。
  !! @param[out] boundary_value 注入面の境界座標値 [m]。
  !! @param[out] inward_normal 注入面の内向き法線ベクトル。
  subroutine resolve_face_sampling_geometry( &
    box_min, box_max, inject_face, axis_n, axis_t1, axis_t2, boundary_value, inward_normal &
    )
    real(dp), intent(in) :: box_min(3), box_max(3)
    character(len=*), intent(in) :: inject_face
    integer, intent(out) :: axis_n, axis_t1, axis_t2
    real(dp), intent(out) :: boundary_value
    real(dp), intent(out) :: inward_normal(3)

    inward_normal = 0.0d0
    select case (trim(lower_ascii(inject_face)))
    case ('x_low')
      axis_n = 1
      axis_t1 = 2
      axis_t2 = 3
      boundary_value = box_min(1)
      inward_normal(1) = 1.0d0
    case ('x_high')
      axis_n = 1
      axis_t1 = 2
      axis_t2 = 3
      boundary_value = box_max(1)
      inward_normal(1) = -1.0d0
    case ('y_low')
      axis_n = 2
      axis_t1 = 3
      axis_t2 = 1
      boundary_value = box_min(2)
      inward_normal(2) = 1.0d0
    case ('y_high')
      axis_n = 2
      axis_t1 = 3
      axis_t2 = 1
      boundary_value = box_max(2)
      inward_normal(2) = -1.0d0
    case ('z_low')
      axis_n = 3
      axis_t1 = 1
      axis_t2 = 2
      boundary_value = box_min(3)
      inward_normal(3) = 1.0d0
    case ('z_high')
      axis_n = 3
      axis_t1 = 1
      axis_t2 = 2
      boundary_value = box_max(3)
      inward_normal(3) = -1.0d0
    case default
      error stop 'Unknown particles.species.inject_face.'
    end select
  end subroutine resolve_face_sampling_geometry

  pure real(dp) function species_temperature_k_local(spec) result(temperature_k)
    type(particle_species_spec), intent(in) :: spec

    temperature_k = spec%temperature_k
    if (spec%has_temperature_ev) temperature_k = spec%temperature_ev*1.160451812d4
  end function species_temperature_k_local

end module bem_app_config_potential_runtime

!> 電位障壁に応じた reservoir 注入速度の補正。
submodule(bem_app_config_particle_runtime) bem_app_config_particle_runtime_sampling
  implicit none
contains

  !> 境界 reservoir 流入に対する法線速度補正パラメータを計算する。
  !! @param[in] cfg シミュレーション・結合設定を含むアプリ設定。
  !! @param[in] spec 境界流入面の粒子種設定。
  !! @param[out] vmin_normal 無限遠法線速度の下限 [m/s]。
  !! @param[out] barrier_normal 法線エネルギー障壁 `2 q Δφ / m` [`m^2/s^2`]。
  !! @param[in] mesh 現在バッチ開始時点の電荷分布メッシュ（補正時に必要）。
  !! @param[inout] snapshot refresh 済み静電 snapshot（infinity barrier 使用時に必要）。
  !! @param[in] warn_face_variation 面平均近似の電位ばらつき警告を出すか。
  module procedure reservoir_face_velocity_correction

  real(dp) :: phi_face, phi_std, phi_min, phi_max, delta_phi
  logical :: emit_warning

  vmin_normal = 0.0d0
  barrier_normal = 0.0d0
  emit_warning = .false.
  if (present(warn_face_variation)) emit_warning = warn_face_variation
  select case (cfg%boundary_contract%inflow_map)
  case (external_inflow_none)
    return
  case (external_inflow_scalar_barrier)
    if (.not. present(mesh)) then
      error stop 'particles.reservoir.inflow_model="infinity_barrier" requires mesh in init_particle_batch_from_config.'
    end if
    if (.not. present(snapshot)) then
      error stop 'particles.reservoir.inflow_model="infinity_barrier" requires a refreshed electrostatic snapshot.'
    end if
    call face_potential_statistics( &
      mesh, cfg, spec, snapshot, phi_face, phi_std, phi_min, phi_max, face_potential, face_potential_ready &
      )
    if (emit_warning) then
      call warn_face_average_potential_variation(cfg%sim, spec, phi_face, phi_std, phi_min, phi_max)
    end if
    delta_phi = phi_face - cfg%sim%phi_infty
    barrier_normal = 2.0d0*spec%q_particle*delta_phi/spec%m_particle
    if (.not. ieee_is_finite(barrier_normal)) then
      error stop 'reservoir potential correction produced non-finite barrier.'
    end if
    if (barrier_normal > 0.0d0) then
      vmin_normal = sqrt(barrier_normal)
    else
      vmin_normal = 0.0d0
    end if
  case default
    error stop 'Unknown external inflow map in runtime.'
  end select
  end procedure reservoir_face_velocity_correction

  !> 外部kinetic closureのaccess bottleneckと、reservoirからbox面までの速度写像を分離して返す。
  module procedure external_kinetic_face_velocity_correction

  real(dp) :: phi_face, phi_std, phi_min, phi_max, access_barrier_normal
  logical :: emit_warning

  if (.not. present(mesh)) then
    error stop 'external kinetic inflow map requires mesh in init_particle_batch_from_config.'
  end if
  if (.not. present(snapshot)) then
    error stop 'external kinetic inflow map requires a refreshed electrostatic snapshot.'
  end if
  if (.not. all(ieee_is_finite([reservoir_potential_v, access_potential_v]))) then
    error stop 'external kinetic inflow map potentials must be finite.'
  end if
  call face_potential_statistics( &
    mesh, cfg, spec, snapshot, phi_face, phi_std, phi_min, phi_max, face_potential, face_potential_ready &
    )
  emit_warning = .false.
  if (present(warn_face_variation)) emit_warning = warn_face_variation
  if (emit_warning) then
    call warn_face_average_potential_variation(cfg%sim, spec, phi_face, phi_std, phi_min, phi_max)
  end if

  barrier_normal = 2.0_dp*spec%q_particle*(phi_face - reservoir_potential_v)/spec%m_particle
  access_barrier_normal = &
    2.0_dp*spec%q_particle*(access_potential_v - reservoir_potential_v)/spec%m_particle
  if (.not. all(ieee_is_finite([barrier_normal, access_barrier_normal]))) then
    error stop 'external kinetic inflow map produced a non-finite velocity barrier.'
  end if
  ! reservoirから現在の注入面までに、Zhaoの外部bottleneckと局所面電位の両方を通過できるtailを数える。
  vmin_normal = sqrt(max(0.0_dp, access_barrier_normal, barrier_normal))
  end procedure external_kinetic_face_velocity_correction

  !> 面の電位統計を、cache があればそこから、なければ評価して返す（評価した値は cache へ書く）。
  subroutine face_potential_statistics( &
    mesh, cfg, spec, snapshot, phi_face, phi_std, phi_min, phi_max, face_potential, face_potential_ready &
    )
    type(mesh_type), intent(in) :: mesh
    type(app_config), intent(in) :: cfg
    type(particle_species_spec), intent(in) :: spec
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
    real(dp), intent(out) :: phi_face, phi_std, phi_min, phi_max
    real(dp), intent(inout), optional :: face_potential(4)
    logical, intent(inout), optional :: face_potential_ready
    logical :: cached

    cached = .false.
    if (present(face_potential) .and. present(face_potential_ready)) cached = face_potential_ready
    if (.not. cached) then
      call compute_face_average_potential(mesh, cfg%sim, spec, snapshot, phi_face, phi_std, phi_min, phi_max)
      if (present(face_potential) .and. present(face_potential_ready)) then
        face_potential = [phi_face, phi_std, phi_min, phi_max]
        face_potential_ready = .true.
      end if
      return
    end if
    phi_face = face_potential(1)
    phi_std = face_potential(2)
    phi_min = face_potential(3)
    phi_max = face_potential(4)
  end subroutine face_potential_statistics

end submodule bem_app_config_particle_runtime_sampling

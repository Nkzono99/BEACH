!> 粒子種別のサンプリングと電位障壁に応じた注入速度の補正。
submodule(bem_app_config_particle_runtime) bem_app_config_particle_runtime_sampling
  implicit none
contains

  !> 1粒子種ぶんの位置・速度サンプルをまとめて生成する。
  !! @param[in] sim ボックス境界・バッチ時間などのシミュレーション設定。
  !! @param[in] spec 1粒子種の注入設定。
  !! @param[in] n 生成粒子数。
  !! @param[out] x 生成した位置配列 `x(3,n)` [m]。
  !! @param[out] v 生成した速度配列 `v(3,n)` [m/s]。
  !! @param[in] barrier_normal_energy reservoir_face 法線方向のエネルギー障壁 `2 q Δφ / m` [`m^2/s^2`]。
  !! @param[in] vmin_normal reservoir_face 法線速度の下限 [m/s]。
  !! @param[in] apply_barrier_energy_shift reservoir_face 法線速度へ障壁エネルギー変換を適用するか。
  module procedure sample_species_state
  logical :: apply_shift
  real(dp) :: temperature_k_local, drift_velocity_local(3), source_box_min(3), source_box_max(3)

  if (n <= 0_i32) return
  apply_shift = .true.
  if (present(apply_barrier_energy_shift)) apply_shift = apply_barrier_energy_shift
  temperature_k_local = species_temperature_k(spec)
  if (present(temperature_k_override)) temperature_k_local = temperature_k_override
  drift_velocity_local = spec%drift_velocity
  if (present(drift_velocity_override)) drift_velocity_local = drift_velocity_override
  source_box_min = sim%box_min
  source_box_max = sim%box_max
  if (trim(lower_ascii(spec%source_mode)) == 'plane_source') then
    call configure_plane_source_box(spec, source_box_min, source_box_max)
  end if
  select case (trim(lower_ascii(spec%source_mode)))
  case ('volume_seed')
    call sample_uniform_positions(spec%pos_low, spec%pos_high, x)
    call sample_shifted_maxwell_velocities( &
      drift_velocity_local, spec%m_particle, v, temperature_k=temperature_k_local &
      )
  case ('reservoir_face', 'plane_source')
    if (trim(lower_ascii(spec%velocity_distribution)) == 'grid') then
      if (present(barrier_normal_energy) .and. present(vmin_normal)) then
        call sample_reservoir_velocity_grid_particles( &
          source_box_min, source_box_max, spec%inject_face, spec%pos_low, spec%pos_high, &
          spec%velocity_grid_path, spec%velocity_grid_pdf_kind, sim%batch_duration, x, v, &
          barrier_normal_energy=barrier_normal_energy, vmin_normal=vmin_normal, position_jitter_dt=sim%dt, &
          apply_barrier_energy_shift=apply_shift, velocity_grid_sampling=spec%velocity_grid_sampling &
          )
      else if (present(barrier_normal_energy)) then
        call sample_reservoir_velocity_grid_particles( &
          source_box_min, source_box_max, spec%inject_face, spec%pos_low, spec%pos_high, &
          spec%velocity_grid_path, spec%velocity_grid_pdf_kind, sim%batch_duration, x, v, &
          barrier_normal_energy=barrier_normal_energy, position_jitter_dt=sim%dt, apply_barrier_energy_shift=apply_shift, &
          velocity_grid_sampling=spec%velocity_grid_sampling &
          )
      else if (present(vmin_normal)) then
        call sample_reservoir_velocity_grid_particles( &
          source_box_min, source_box_max, spec%inject_face, spec%pos_low, spec%pos_high, &
          spec%velocity_grid_path, spec%velocity_grid_pdf_kind, sim%batch_duration, x, v, &
          vmin_normal=vmin_normal, position_jitter_dt=sim%dt, apply_barrier_energy_shift=apply_shift, &
          velocity_grid_sampling=spec%velocity_grid_sampling &
          )
      else
        call sample_reservoir_velocity_grid_particles( &
          source_box_min, source_box_max, spec%inject_face, spec%pos_low, spec%pos_high, &
          spec%velocity_grid_path, spec%velocity_grid_pdf_kind, sim%batch_duration, x, v, position_jitter_dt=sim%dt, &
          apply_barrier_energy_shift=apply_shift, velocity_grid_sampling=spec%velocity_grid_sampling &
          )
      end if
    else if (present(barrier_normal_energy) .and. present(vmin_normal)) then
      call sample_reservoir_face_particles( &
        source_box_min, source_box_max, spec%inject_face, spec%pos_low, spec%pos_high, drift_velocity_local, &
        spec%m_particle, temperature_k_local, sim%batch_duration, x, v, &
        barrier_normal_energy=barrier_normal_energy, vmin_normal=vmin_normal, position_jitter_dt=sim%dt, &
        apply_barrier_energy_shift=apply_shift &
        )
    else if (present(barrier_normal_energy)) then
      call sample_reservoir_face_particles( &
        source_box_min, source_box_max, spec%inject_face, spec%pos_low, spec%pos_high, drift_velocity_local, &
        spec%m_particle, temperature_k_local, sim%batch_duration, x, v, &
        barrier_normal_energy=barrier_normal_energy, position_jitter_dt=sim%dt, &
        apply_barrier_energy_shift=apply_shift &
        )
    else if (present(vmin_normal)) then
      call sample_reservoir_face_particles( &
        source_box_min, source_box_max, spec%inject_face, spec%pos_low, spec%pos_high, drift_velocity_local, &
        spec%m_particle, temperature_k_local, sim%batch_duration, x, v, &
        vmin_normal=vmin_normal, position_jitter_dt=sim%dt, apply_barrier_energy_shift=apply_shift &
        )
    else
      call sample_reservoir_face_particles( &
        source_box_min, source_box_max, spec%inject_face, spec%pos_low, spec%pos_high, drift_velocity_local, &
        spec%m_particle, temperature_k_local, sim%batch_duration, x, v, position_jitter_dt=sim%dt, &
        apply_barrier_energy_shift=apply_shift &
        )
    end if
  case ('photo_raycast')
    error stop 'sample_species_state does not support photo_raycast. Use sample_photo_species_state.'
  case default
    error stop 'Unknown particles.species.source_mode.'
  end select
  if (trim(lower_ascii(spec%source_mode)) == 'reservoir_face' .or. &
      trim(lower_ascii(spec%source_mode)) == 'plane_source') call normalize_reservoir_positions(sim, x)
  end procedure sample_species_state

  !> photo_raycast 粒子種のレイキャスト放出を実行する。
  !! @param[in] sim シミュレーション設定。
  !! @param[in] spec photo_raycast 粒子種設定。
  !! @param[in] mesh 交差判定に使う現在メッシュ。
  !! @param[in] n_rays バッチで発射するレイ本数。
  !! @param[out] x 生成した位置配列 `x(3,n_rays)` [m]。
  !! @param[out] v 生成した速度配列 `v(3,n_rays)` [m/s]。
  !! @param[out] w 生成した重み配列 `w(n_rays)`。
  !! @param[out] n_emit 実際に放出された粒子数。
  !! @param[out] emit_elem_idx 放出元要素ID `emit_elem_idx(n_rays)`（省略可）。
  !! @param[in] emit_current_density_override 放出電流密度の上書き値 [A/m^2]。
  !! @param[in] normal_drift_speed_override 放出法線ドリフトの上書き値 [m/s]。
  !! @param[in] vmin_normal 放出法線速度の下限 [m/s]。
  !! @param[out] collision_failure_status 不完全な photo collision query の status（省略時は停止）。
  !! @param[out] collision_failure_ray 不完全な照会を返した最小 ray index。
  !! @param[out] collision_failure_bounce 不完全な照会を返した bounce index。
  module procedure sample_photo_species_state
  real(dp) :: emit_current_density, normal_drift_speed

  if (present(collision_failure_status)) collision_failure_status = collision_query_ok
  if (present(collision_failure_ray)) collision_failure_ray = huge(0_i32)
  if (present(collision_failure_bounce)) collision_failure_bounce = 0_i32
  if (n_rays <= 0_i32) then
    if (present(emit_elem_idx)) emit_elem_idx = -1_i32
    n_emit = 0_i32
    return
  end if
  emit_current_density = spec%emit_current_density_a_m2
  if (present(emit_current_density_override)) emit_current_density = emit_current_density_override
  normal_drift_speed = spec%normal_drift_speed
  if (present(normal_drift_speed_override)) normal_drift_speed = normal_drift_speed_override
  if (present(global_rays_per_batch)) then
    if (present(vmin_normal)) then
      call sample_photo_raycast_particles( &
        mesh, sim, spec%inject_face, spec%pos_low, spec%pos_high, spec%ray_direction, spec%m_particle, &
        species_temperature_k(spec), normal_drift_speed, emit_current_density, spec%q_particle, &
        n_rays, x, v, w, n_emit, emit_elem_idx, global_rays_per_batch=global_rays_per_batch, vmin_normal=vmin_normal, &
        collision_failure_status=collision_failure_status, collision_failure_ray=collision_failure_ray, &
        collision_failure_bounce=collision_failure_bounce &
        )
    else
      call sample_photo_raycast_particles( &
        mesh, sim, spec%inject_face, spec%pos_low, spec%pos_high, spec%ray_direction, spec%m_particle, &
        species_temperature_k(spec), normal_drift_speed, emit_current_density, spec%q_particle, &
        n_rays, x, v, w, n_emit, emit_elem_idx, global_rays_per_batch=global_rays_per_batch, &
        collision_failure_status=collision_failure_status, collision_failure_ray=collision_failure_ray, &
        collision_failure_bounce=collision_failure_bounce &
        )
    end if
  else
    if (present(vmin_normal)) then
      call sample_photo_raycast_particles( &
        mesh, sim, spec%inject_face, spec%pos_low, spec%pos_high, spec%ray_direction, spec%m_particle, &
        species_temperature_k(spec), normal_drift_speed, emit_current_density, spec%q_particle, &
        n_rays, x, v, w, n_emit, emit_elem_idx, vmin_normal=vmin_normal, &
        collision_failure_status=collision_failure_status, collision_failure_ray=collision_failure_ray, &
        collision_failure_bounce=collision_failure_bounce &
        )
    else
      call sample_photo_raycast_particles( &
        mesh, sim, spec%inject_face, spec%pos_low, spec%pos_high, spec%ray_direction, spec%m_particle, &
        species_temperature_k(spec), normal_drift_speed, emit_current_density, spec%q_particle, &
        n_rays, x, v, w, n_emit, emit_elem_idx, &
        collision_failure_status=collision_failure_status, collision_failure_ray=collision_failure_ray, &
        collision_failure_bounce=collision_failure_bounce &
        )
    end if
  end if
  end procedure sample_photo_species_state

  !> reservoir_face の target 個数からシース補正込み重みを解決する。
  !> reservoir_face 注入に対する法線速度補正パラメータを計算する。
  !! @param[in] cfg シミュレーション・結合設定を含むアプリ設定。
  !! @param[in] spec reservoir_face 粒子種設定。
  !! @param[out] vmin_normal 無限遠法線速度の下限 [m/s]。
  !! @param[out] barrier_normal 法線エネルギー障壁 `2 q Δφ / m` [`m^2/s^2`]。
  !! @param[in] mesh 現在バッチ開始時点の電荷分布メッシュ（補正時に必要）。
  !! @param[inout] snapshot refresh 済み静電 snapshot（infinity barrier 使用時に必要）。
  !! @param[in] warn_face_variation 面平均近似の電位ばらつき警告を出すか。
  module procedure reservoir_face_velocity_correction

  real(dp) :: phi_face, phi_std, phi_min, phi_max, delta_phi
  integer(i32) :: boundary_status
  logical :: emit_warning
  type(external_boundary_contract_type) :: active_boundary_contract
  character(len=256) :: boundary_message

  vmin_normal = 0.0d0
  barrier_normal = 0.0d0
  emit_warning = .false.
  if (present(warn_face_variation)) emit_warning = warn_face_variation
  if (present(boundary_contract)) then
    active_boundary_contract = boundary_contract
  else
    call resolve_external_boundary_contract( &
      cfg%sim%reservoir_potential_model, cfg%sim%open_boundary_model, &
      active_boundary_contract, boundary_status, boundary_message &
      )
    if (boundary_status /= external_boundary_ok) error stop trim(boundary_message)
  end if
  select case (active_boundary_contract%inflow_map)
  case (external_inflow_none)
    return
  case (external_inflow_scalar_barrier)
    if (.not. present(mesh)) then
      error stop 'sim.reservoir_potential_model="infinity_barrier" requires mesh in init_particle_batch_from_config.'
    end if
    if (.not. present(snapshot)) then
      error stop 'sim.reservoir_potential_model="infinity_barrier" requires a refreshed electrostatic snapshot.'
    end if
    call compute_face_average_potential(mesh, cfg%sim, spec, snapshot, phi_face, phi_std, phi_min, phi_max)
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
  call compute_face_average_potential(mesh, cfg%sim, spec, snapshot, phi_face, phi_std, phi_min, phi_max)
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

end submodule bem_app_config_particle_runtime_sampling

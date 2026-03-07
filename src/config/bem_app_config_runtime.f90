!> `app_config` からメッシュ・粒子群を構築する実行時変換モジュール。
module bem_app_config_runtime
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, particles_soa, sim_config, injection_state
  use bem_field, only: electric_potential_at
  use bem_templates, only: make_plane, make_box, make_cylinder, make_sphere
  use bem_mesh, only: init_mesh
  use bem_importers, only: load_obj_mesh
  use bem_injection, only: &
    seed_rng, sample_uniform_positions, sample_shifted_maxwell_velocities, compute_macro_particles_for_batch, &
    sample_reservoir_face_particles, sample_photo_raycast_particles
  use bem_particles, only: init_particles
  use bem_app_config_types, only: &
    app_config, particle_species_spec, template_spec, max_templates, particles_per_batch_from_config, &
    total_particles_from_config
  use bem_app_config_parser, only: lower
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

contains

  !> `mesh_mode` と OBJ ファイル有無に応じてメッシュ生成方法を選ぶ。
  !! @param[in] cfg メッシュ入力設定を含むアプリ設定。
  !! @param[out] mesh 構築した三角形メッシュ。
  subroutine build_mesh_from_config(cfg, mesh)
    type(app_config), intent(in) :: cfg
    type(mesh_type), intent(out) :: mesh
    logical :: has_obj

    select case (trim(lower(cfg%mesh_mode)))
    case ('obj')
      call load_obj_mesh(trim(cfg%obj_path), mesh)
    case ('template')
      call build_template_mesh(cfg, mesh)
    case default
      inquire (file=trim(cfg%obj_path), exist=has_obj)
      if (has_obj) then
        call load_obj_mesh(trim(cfg%obj_path), mesh)
      else
        call build_template_mesh(cfg, mesh)
      end if
    end select
  end subroutine build_mesh_from_config

  !> 設定全体ぶんの粒子群を生成し、SoA へ詰める。
  !! 粒子種ごとに乱数サンプルした後、種ごとに rank を揃えて interleave する。
  !! @param[in] cfg 粒子種設定を含むアプリ設定。
  !! @param[out] pcls 生成した全粒子群。
  subroutine init_particles_from_config(cfg, pcls)
    type(app_config), intent(in) :: cfg
    type(particles_soa), intent(out) :: pcls

    integer(i32) :: s, total_n, max_n, out_idx, rank_idx
    integer(i32), allocatable :: counts(:)
    real(dp), allocatable :: x_species(:, :, :), v_species(:, :, :)
    real(dp), allocatable :: q(:), m(:), w(:), x(:, :), v(:, :)

    if (cfg%n_particle_species <= 0) error stop 'At least one [[particles.species]] entry is required.'

    call seed_rng([cfg%sim%rng_seed])

    allocate (counts(cfg%n_particle_species))
    counts = 0_i32
    do s = 1, cfg%n_particle_species
      if (cfg%particle_species(s)%enabled) then
        if (trim(lower(cfg%particle_species(s)%source_mode)) == 'reservoir_face' .or. &
            trim(lower(cfg%particle_species(s)%source_mode)) == 'photo_raycast') then
          error stop 'init_particles_from_config supports volume_seed only. Use init_particle_batch_from_config.'
        end if
        if (cfg%particle_species(s)%npcls_per_step < 0_i32) then
          error stop 'particles.species.npcls_per_step must be >= 0.'
        end if
        counts(s) = cfg%sim%batch_count * cfg%particle_species(s)%npcls_per_step
      end if
    end do

    total_n = total_particles_from_config(cfg)
    max_n = max(1_i32, maxval(counts))

    allocate (x_species(3, max_n, cfg%n_particle_species))
    allocate (v_species(3, max_n, cfg%n_particle_species))
    x_species = 0.0d0
    v_species = 0.0d0

    do s = 1, cfg%n_particle_species
      if (counts(s) <= 0_i32) cycle
      call sample_species_state( &
        cfg%sim, cfg%particle_species(s), counts(s), x_species(:, 1:counts(s), s), v_species(:, 1:counts(s), s) &
      )
    end do

    allocate (x(3, total_n), v(3, total_n), q(total_n), m(total_n), w(total_n))
    out_idx = 0_i32
    do rank_idx = 1, maxval(counts)
      do s = 1, cfg%n_particle_species
        if (rank_idx > counts(s)) cycle
        out_idx = out_idx + 1_i32
        x(:, out_idx) = x_species(:, rank_idx, s)
        v(:, out_idx) = v_species(:, rank_idx, s)
        q(out_idx) = cfg%particle_species(s)%q_particle
        m(out_idx) = cfg%particle_species(s)%m_particle
        w(out_idx) = cfg%particle_species(s)%w_particle
      end do
    end do

    call init_particles(pcls, x, v, q, m, w)
  end subroutine init_particles_from_config

  !> バッチ生成前に乱数シードだけを初期化する。
  !! @param[in] cfg 乱数シード値 `sim.rng_seed` を含むアプリ設定。
  subroutine seed_particles_from_config(cfg, mpi_rank, mpi_size)
    type(app_config), intent(in) :: cfg
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    integer(i32) :: local_rank, n_ranks, seed_value
    integer(kind=8) :: seed_tmp

    local_rank = 0_i32
    n_ranks = 1_i32
    if (present(mpi_rank)) local_rank = mpi_rank
    if (present(mpi_size)) n_ranks = mpi_size
    if (n_ranks <= 0_i32) error stop 'mpi_size must be > 0 in seed_particles_from_config.'
    if (local_rank < 0_i32 .or. local_rank >= n_ranks) error stop 'mpi_rank is out of range in seed_particles_from_config.'

    seed_tmp = int(cfg%sim%rng_seed, kind=8) + 104729_8 * int(local_rank, kind=8)
    seed_value = int(modulo(seed_tmp, int(huge(0_i32), kind=8)), kind=i32)
    call seed_rng([seed_value])
  end subroutine seed_particles_from_config

  !> 指定バッチ番号に対応する粒子バッチを生成する。
  !! @param[in] cfg 粒子種とシミュレーション条件を含むアプリ設定。
  !! @param[in] batch_idx 生成対象のバッチ番号（1始まり）。
  !! @param[out] pcls 生成したバッチ粒子群。
  !! @param[inout] state reservoir_face 注入の残差状態（必要時のみ）。
  !! @param[in] mesh 現在バッチ開始時点の電荷分布メッシュ（電位補正時に必要）。
  !! @param[out] photo_emission_dq photo_raycast 放出起因の要素電荷差分 `photo_emission_dq(nelem)`（省略可）。
  subroutine init_particle_batch_from_config(cfg, batch_idx, pcls, state, mesh, photo_emission_dq, mpi_rank, mpi_size)
    type(app_config), intent(in) :: cfg
    integer(i32), intent(in) :: batch_idx
    type(particles_soa), intent(out) :: pcls
    type(injection_state), intent(inout), optional :: state
    type(mesh_type), intent(in), optional :: mesh
    real(dp), intent(out), optional :: photo_emission_dq(:)
    integer(i32), intent(in), optional :: mpi_rank, mpi_size

    integer(i32) :: s, i, batch_n, max_rank, out_idx, local_rank, n_ranks, global_count
    integer(i32), allocatable :: counts_max(:), counts_actual(:), species_cursor(:), species_id(:), emit_elem_species(:, :)
    real(dp), allocatable :: vmin_normal(:), barrier_normal(:)
    real(dp), allocatable :: x_species(:, :, :), v_species(:, :, :), w_species(:, :)
    real(dp), allocatable :: x(:, :), v(:, :), q(:), m(:), w(:)

    if (cfg%sim%batch_count <= 0_i32) error stop 'sim.batch_count must be > 0.'
    if (batch_idx < 1_i32 .or. batch_idx > cfg%sim%batch_count) then
      error stop 'Requested batch index is out of range.'
    end if
    local_rank = 0_i32
    n_ranks = 1_i32
    if (present(mpi_rank)) local_rank = mpi_rank
    if (present(mpi_size)) n_ranks = mpi_size
    if (n_ranks <= 0_i32) error stop 'mpi_size must be > 0 in init_particle_batch_from_config.'
    if (local_rank < 0_i32 .or. local_rank >= n_ranks) error stop 'mpi_rank is out of range in init_particle_batch_from_config.'
    if (present(state)) then
      if (.not. allocated(state%macro_residual)) error stop 'injection_state is not initialized.'
      if (size(state%macro_residual) < cfg%n_particle_species) error stop 'injection_state size mismatch.'
    end if
    if (present(photo_emission_dq)) then
      if (.not. present(mesh)) error stop 'photo_emission_dq requires mesh in init_particle_batch_from_config.'
      if (size(photo_emission_dq) /= mesh%nelem) error stop 'photo_emission_dq size mismatch.'
      photo_emission_dq = 0.0d0
    end if

    allocate (counts_max(cfg%n_particle_species), counts_actual(cfg%n_particle_species))
    allocate (vmin_normal(cfg%n_particle_species), barrier_normal(cfg%n_particle_species))
    counts_max = 0_i32
    counts_actual = 0_i32
    vmin_normal = 0.0d0
    barrier_normal = 0.0d0
    do s = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(s)%enabled) cycle
      select case (trim(lower(cfg%particle_species(s)%source_mode)))
      case ('volume_seed')
        global_count = cfg%particle_species(s)%npcls_per_step
        counts_max(s) = split_count_for_rank(global_count, local_rank, n_ranks)
      case ('reservoir_face')
        if (.not. present(state)) then
          error stop 'reservoir_face requires injection_state in init_particle_batch_from_config.'
        end if
        call reservoir_face_velocity_correction( &
          cfg%sim, cfg%particle_species(s), vmin_normal(s), barrier_normal(s), mesh &
        )
        call compute_macro_particles_for_species( &
          cfg%sim, cfg%particle_species(s), state%macro_residual(s), counts_max(s), vmin_normal=vmin_normal(s), &
          batch_duration_scale=1.0d0/real(n_ranks, dp) &
        )
      case ('photo_raycast')
        global_count = cfg%particle_species(s)%rays_per_batch
        counts_max(s) = split_count_for_rank(global_count, local_rank, n_ranks)
      case default
        error stop 'Unknown particles.species.source_mode.'
      end select
    end do
    max_rank = max(1_i32, maxval(counts_max))
    allocate (x_species(3, max_rank, cfg%n_particle_species))
    allocate (v_species(3, max_rank, cfg%n_particle_species))
    allocate (w_species(max_rank, cfg%n_particle_species))
    allocate (emit_elem_species(max_rank, cfg%n_particle_species))
    x_species = 0.0d0
    v_species = 0.0d0
    w_species = 0.0d0
    emit_elem_species = -1_i32
    do s = 1, cfg%n_particle_species
      if (counts_max(s) <= 0_i32) cycle
      select case (trim(lower(cfg%particle_species(s)%source_mode)))
      case ('volume_seed', 'reservoir_face')
        call sample_species_state( &
          cfg%sim, cfg%particle_species(s), counts_max(s), x_species(:, 1:counts_max(s), s), v_species(:, 1:counts_max(s), s), &
          barrier_normal_energy=barrier_normal(s), vmin_normal=vmin_normal(s) &
        )
        counts_actual(s) = counts_max(s)
        w_species(1:counts_actual(s), s) = cfg%particle_species(s)%w_particle
      case ('photo_raycast')
        if (.not. present(mesh)) then
          error stop 'photo_raycast requires mesh in init_particle_batch_from_config.'
        end if
        call sample_photo_species_state( &
          cfg%sim, cfg%particle_species(s), mesh, counts_max(s), x_species(:, 1:counts_max(s), s), &
          v_species(:, 1:counts_max(s), s), w_species(1:counts_max(s), s), counts_actual(s), &
          emit_elem_idx=emit_elem_species(1:counts_max(s), s), &
          global_rays_per_batch=cfg%particle_species(s)%rays_per_batch &
        )
        if (present(photo_emission_dq) .and. cfg%particle_species(s)%deposit_opposite_charge_on_emit) then
          do i = 1, counts_actual(s)
            if (emit_elem_species(i, s) < 1_i32 .or. emit_elem_species(i, s) > size(photo_emission_dq)) then
              error stop 'photo_raycast emitted invalid elem_idx.'
            end if
            photo_emission_dq(emit_elem_species(i, s)) = photo_emission_dq(emit_elem_species(i, s)) - &
                                                         cfg%particle_species(s)%q_particle * w_species(i, s)
          end do
        end if
      case default
        error stop 'Unknown particles.species.source_mode.'
      end select
    end do

    batch_n = sum(counts_actual)
    allocate (species_id(batch_n))
    out_idx = 0_i32
    do i = 1, max_rank
      do s = 1, cfg%n_particle_species
        if (i > counts_actual(s)) cycle
        out_idx = out_idx + 1_i32
        species_id(out_idx) = s
      end do
    end do

    allocate (x(3, batch_n), v(3, batch_n), q(batch_n), m(batch_n), w(batch_n))
    allocate (species_cursor(cfg%n_particle_species))
    species_cursor = 0_i32
    do i = 1, batch_n
      s = species_id(i)
      species_cursor(s) = species_cursor(s) + 1_i32
      x(:, i) = x_species(:, species_cursor(s), s)
      v(:, i) = v_species(:, species_cursor(s), s)
      q(i) = cfg%particle_species(s)%q_particle
      m(i) = cfg%particle_species(s)%m_particle
      w(i) = w_species(species_cursor(s), s)
    end do

    call init_particles(pcls, x, v, q, m, w)
  end subroutine init_particle_batch_from_config

  !> 1粒子種ぶんの位置・速度サンプルをまとめて生成する。
  !! @param[in] sim ボックス境界・バッチ時間などのシミュレーション設定。
  !! @param[in] spec 1粒子種の注入設定。
  !! @param[in] n 生成粒子数。
  !! @param[out] x 生成した位置配列 `x(3,n)` [m]。
  !! @param[out] v 生成した速度配列 `v(3,n)` [m/s]。
  !! @param[in] barrier_normal_energy reservoir_face 法線方向のエネルギー障壁 `2 q Δφ / m` [`m^2/s^2`]。
  !! @param[in] vmin_normal reservoir_face 法線速度の下限 [m/s]。
  subroutine sample_species_state(sim, spec, n, x, v, barrier_normal_energy, vmin_normal)
    type(sim_config), intent(in) :: sim
    type(particle_species_spec), intent(in) :: spec
    integer(i32), intent(in) :: n
    real(dp), intent(out) :: x(:, :)
    real(dp), intent(out) :: v(:, :)
    real(dp), intent(in), optional :: barrier_normal_energy
    real(dp), intent(in), optional :: vmin_normal

    if (n <= 0_i32) return
    select case (trim(lower(spec%source_mode)))
    case ('volume_seed')
      call sample_uniform_positions(spec%pos_low, spec%pos_high, x)
      call sample_shifted_maxwell_velocities(spec%drift_velocity, spec%m_particle, v, temperature_k=species_temperature_k(spec))
    case ('reservoir_face')
      if (present(barrier_normal_energy) .and. present(vmin_normal)) then
        call sample_reservoir_face_particles( &
          sim%box_min, sim%box_max, spec%inject_face, spec%pos_low, spec%pos_high, spec%drift_velocity, &
          spec%m_particle, species_temperature_k(spec), sim%batch_duration, x, v, &
          barrier_normal_energy=barrier_normal_energy, vmin_normal=vmin_normal, position_jitter_dt=sim%dt &
        )
      else if (present(barrier_normal_energy)) then
        call sample_reservoir_face_particles( &
          sim%box_min, sim%box_max, spec%inject_face, spec%pos_low, spec%pos_high, spec%drift_velocity, &
          spec%m_particle, species_temperature_k(spec), sim%batch_duration, x, v, &
          barrier_normal_energy=barrier_normal_energy, position_jitter_dt=sim%dt &
        )
      else if (present(vmin_normal)) then
        call sample_reservoir_face_particles( &
          sim%box_min, sim%box_max, spec%inject_face, spec%pos_low, spec%pos_high, spec%drift_velocity, &
          spec%m_particle, species_temperature_k(spec), sim%batch_duration, x, v, &
          vmin_normal=vmin_normal, position_jitter_dt=sim%dt &
        )
      else
        call sample_reservoir_face_particles( &
          sim%box_min, sim%box_max, spec%inject_face, spec%pos_low, spec%pos_high, spec%drift_velocity, &
          spec%m_particle, species_temperature_k(spec), sim%batch_duration, x, v, position_jitter_dt=sim%dt &
        )
      end if
    case ('photo_raycast')
      error stop 'sample_species_state does not support photo_raycast. Use sample_photo_species_state.'
    case default
      error stop 'Unknown particles.species.source_mode.'
    end select
  end subroutine sample_species_state

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
  subroutine sample_photo_species_state(sim, spec, mesh, n_rays, x, v, w, n_emit, emit_elem_idx, global_rays_per_batch)
    type(sim_config), intent(in) :: sim
    type(particle_species_spec), intent(in) :: spec
    type(mesh_type), intent(in) :: mesh
    integer(i32), intent(in) :: n_rays
    real(dp), intent(out) :: x(:, :)
    real(dp), intent(out) :: v(:, :)
    real(dp), intent(out) :: w(:)
    integer(i32), intent(out) :: n_emit
    integer(i32), intent(out), optional :: emit_elem_idx(:)
    integer(i32), intent(in), optional :: global_rays_per_batch

    if (n_rays <= 0_i32) then
      if (present(emit_elem_idx)) emit_elem_idx = -1_i32
      n_emit = 0_i32
      return
    end if
    if (present(global_rays_per_batch)) then
      call sample_photo_raycast_particles( &
        mesh, sim, spec%inject_face, spec%pos_low, spec%pos_high, spec%ray_direction, spec%m_particle, &
        species_temperature_k(spec), spec%normal_drift_speed, spec%emit_current_density_a_m2, spec%q_particle, &
        n_rays, x, v, w, n_emit, emit_elem_idx, global_rays_per_batch=global_rays_per_batch &
      )
    else
      call sample_photo_raycast_particles( &
        mesh, sim, spec%inject_face, spec%pos_low, spec%pos_high, spec%ray_direction, spec%m_particle, &
        species_temperature_k(spec), spec%normal_drift_speed, spec%emit_current_density_a_m2, spec%q_particle, &
        n_rays, x, v, w, n_emit, emit_elem_idx &
      )
    end if
  end subroutine sample_photo_species_state

  !> reservoir_face 用に、物理流量と残差から今バッチのマクロ粒子数を決める。
  !! @param[in] sim ボックス境界・バッチ時間などのシミュレーション設定。
  !! @param[in] spec reservoir_face 粒子種設定。
  !! @param[inout] residual 前バッチから繰り越した端数。
  !! @param[out] count 今バッチで生成するマクロ粒子数。
  !! @param[in] vmin_normal 法線速度の下限 [m/s]（省略時は 0）。
  subroutine compute_macro_particles_for_species(sim, spec, residual, count, vmin_normal, batch_duration_scale)
    type(sim_config), intent(in) :: sim
    type(particle_species_spec), intent(in) :: spec
    real(dp), intent(inout) :: residual
    integer(i32), intent(out) :: count
    real(dp), intent(in), optional :: vmin_normal
    real(dp), intent(in), optional :: batch_duration_scale

    real(dp) :: number_density_m3, effective_batch_duration

    number_density_m3 = species_number_density_m3(spec)
    effective_batch_duration = sim%batch_duration
    if (present(batch_duration_scale)) effective_batch_duration = sim%batch_duration * batch_duration_scale
    if (present(vmin_normal)) then
      call compute_macro_particles_for_batch( &
        number_density_m3, species_temperature_k(spec), spec%m_particle, spec%drift_velocity, sim%box_min, sim%box_max, &
        spec%inject_face, spec%pos_low, spec%pos_high, effective_batch_duration, spec%w_particle, residual, count, &
        vmin_normal=vmin_normal &
      )
    else
      call compute_macro_particles_for_batch( &
        number_density_m3, species_temperature_k(spec), spec%m_particle, spec%drift_velocity, sim%box_min, sim%box_max, &
        spec%inject_face, spec%pos_low, spec%pos_high, effective_batch_duration, spec%w_particle, residual, count &
      )
    end if
  end subroutine compute_macro_particles_for_species

  !> reservoir_face 注入に対する法線速度補正パラメータを計算する。
  !! @param[in] sim シミュレーション設定。
  !! @param[in] spec reservoir_face 粒子種設定。
  !! @param[out] vmin_normal 無限遠法線速度の下限 [m/s]。
  !! @param[out] barrier_normal 法線エネルギー障壁 `2 q Δφ / m` [`m^2/s^2`]。
  !! @param[in] mesh 現在バッチ開始時点の電荷分布メッシュ（補正時に必要）。
  subroutine reservoir_face_velocity_correction(sim, spec, vmin_normal, barrier_normal, mesh)
    type(sim_config), intent(in) :: sim
    type(particle_species_spec), intent(in) :: spec
    real(dp), intent(out) :: vmin_normal
    real(dp), intent(out) :: barrier_normal
    type(mesh_type), intent(in), optional :: mesh

    real(dp) :: phi_face, delta_phi

    vmin_normal = 0.0d0
    barrier_normal = 0.0d0
    select case (trim(lower(sim%reservoir_potential_model)))
    case ('none')
      return
    case ('infinity_barrier')
      if (.not. present(mesh)) then
        error stop 'sim.reservoir_potential_model="infinity_barrier" requires mesh in init_particle_batch_from_config.'
      end if
      call compute_face_average_potential(mesh, sim, spec, phi_face)
      delta_phi = phi_face - sim%phi_infty
      barrier_normal = 2.0d0 * spec%q_particle * delta_phi / spec%m_particle
      if (.not. ieee_is_finite(barrier_normal)) then
        error stop 'reservoir potential correction produced non-finite barrier.'
      end if
      if (barrier_normal > 0.0d0) then
        vmin_normal = sqrt(barrier_normal)
      else
        vmin_normal = 0.0d0
      end if
    case default
      error stop 'Unknown sim.reservoir_potential_model in runtime.'
    end select
  end subroutine reservoir_face_velocity_correction

  !> reservoir_face 開口面の平均電位を `N x N` 格子平均で評価する。
  !! @param[in] mesh 現在バッチ開始時点の電荷分布メッシュ。
  !! @param[in] sim シミュレーション設定。
  !! @param[in] spec reservoir_face 粒子種設定。
  !! @param[out] phi_face 注入開口面の平均電位 [V]。
  subroutine compute_face_average_potential(mesh, sim, spec, phi_face)
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(particle_species_spec), intent(in) :: spec
    real(dp), intent(out) :: phi_face

    integer(i32) :: ngrid, i, j
    integer :: axis_n, axis_t1, axis_t2
    real(dp) :: boundary_value, inward_normal(3), pos(3), t1, t2, phi
    real(dp) :: phi_sum

    call resolve_face_sampling_geometry( &
      sim%box_min, sim%box_max, spec%inject_face, axis_n, axis_t1, axis_t2, boundary_value, inward_normal &
    )

    ngrid = sim%injection_face_phi_grid_n
    phi_sum = 0.0d0
    do i = 1_i32, ngrid
      t1 = (real(i, dp) - 0.5d0) / real(ngrid, dp)
      do j = 1_i32, ngrid
        t2 = (real(j, dp) - 0.5d0) / real(ngrid, dp)
        pos = 0.0d0
        pos(axis_n) = boundary_value
        pos(axis_t1) = spec%pos_low(axis_t1) + (spec%pos_high(axis_t1) - spec%pos_low(axis_t1)) * t1
        pos(axis_t2) = spec%pos_low(axis_t2) + (spec%pos_high(axis_t2) - spec%pos_low(axis_t2)) * t2
        pos = pos + inward_normal * 1.0d-12
        call electric_potential_at(mesh, pos, sim%softening, phi)
        phi_sum = phi_sum + phi
      end do
    end do

    phi_face = phi_sum / real(ngrid * ngrid, dp)
  end subroutine compute_face_average_potential

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
    select case (trim(lower(inject_face)))
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

  !> 粒子種設定から実効密度[m^-3]を返す。
  !! @param[in] spec 粒子種設定。
  !! @return number_density_m3 実効粒子数密度 [1/m^3]。
  pure real(dp) function species_number_density_m3(spec) result(number_density_m3)
    type(particle_species_spec), intent(in) :: spec

    number_density_m3 = spec%number_density_m3
    if (spec%has_number_density_cm3) number_density_m3 = spec%number_density_cm3 * 1.0d6
  end function species_number_density_m3

  !> 総数 `total_count` を `n_ranks` 個のrankへできるだけ均等に分配し、このrank分を返す。
  integer(i32) function split_count_for_rank(total_count, rank, n_ranks) result(local_count)
    integer(i32), intent(in) :: total_count, rank, n_ranks
    integer(i32) :: base_count, n_remainder

    if (total_count < 0_i32) error stop 'split_count_for_rank requires total_count >= 0.'
    if (n_ranks <= 0_i32) error stop 'split_count_for_rank requires n_ranks > 0.'
    if (rank < 0_i32 .or. rank >= n_ranks) error stop 'split_count_for_rank rank out of range.'

    base_count = total_count/n_ranks
    n_remainder = modulo(total_count, n_ranks)
    local_count = base_count
    if (rank < n_remainder) local_count = local_count + 1_i32
  end function split_count_for_rank

  !> 粒子種設定から実効温度[K]を返す。
  !! @param[in] spec 粒子種設定。
  !! @return temperature_k 実効温度 [K]。
  pure real(dp) function species_temperature_k(spec) result(temperature_k)
    type(particle_species_spec), intent(in) :: spec

    temperature_k = spec%temperature_k
    if (spec%has_temperature_ev) temperature_k = spec%temperature_ev * 1.160451812d4
  end function species_temperature_k

  !> 有効なテンプレートを連結し、1つのメッシュへまとめる。
  !! @param[in] cfg テンプレート設定を含むアプリ設定。
  !! @param[out] mesh 連結後の三角形メッシュ。
  subroutine build_template_mesh(cfg, mesh)
    type(app_config), intent(in) :: cfg
    type(mesh_type), intent(out) :: mesh
    type(mesh_type) :: part
    real(dp), allocatable :: v0(:, :), v1(:, :), v2(:, :)
    integer(i32) :: i

    allocate (v0(3, 0), v1(3, 0), v2(3, 0))
    do i = 1, min(cfg%n_templates, int(max_templates, i32))
      if (.not. cfg%templates(i)%enabled) cycle
      call build_one_template(cfg%templates(i), part)
      call append_triangles(v0, v1, v2, part%v0, part%v1, part%v2)
    end do

    if (size(v0, 2) == 0) then
      error stop 'No enabled template found in configuration.'
    end if
    call init_mesh(mesh, v0, v1, v2)
  end subroutine build_template_mesh

  !> テンプレート種別に応じて形状生成ルーチンへディスパッチする。
  !! @param[in] spec 1テンプレート分の形状設定。
  !! @param[out] mesh 生成したテンプレートメッシュ。
  subroutine build_one_template(spec, mesh)
    type(template_spec), intent(in) :: spec
    type(mesh_type), intent(out) :: mesh

    select case (trim(lower(spec%kind)))
    case ('plane')
      call make_plane(mesh, size_x=spec%size_x, size_y=spec%size_y, nx=spec%nx, ny=spec%ny, center=spec%center)
    case ('box')
      call make_box(mesh, size=spec%size, center=spec%center, nx=spec%nx, ny=spec%ny, nz=spec%nz)
    case ('cylinder')
      call make_cylinder( &
        mesh, radius=spec%radius, height=spec%height, n_theta=spec%n_theta, n_z=spec%n_z, &
        cap=spec%cap, center=spec%center &
      )
    case ('sphere')
      call make_sphere(mesh, radius=spec%radius, n_lon=spec%n_lon, n_lat=spec%n_lat, center=spec%center)
    case default
      error stop 'Unknown template kind in config.'
    end select
  end subroutine build_one_template

  !> 既存三角形配列へ追加分を連結し、再確保後の配列へ差し替える。
  !! @param[inout] v0 累積メッシュの頂点0配列 `v0(3,n)`。
  !! @param[inout] v1 累積メッシュの頂点1配列 `v1(3,n)`。
  !! @param[inout] v2 累積メッシュの頂点2配列 `v2(3,n)`。
  !! @param[in] add_v0 追加する頂点0配列。
  !! @param[in] add_v1 追加する頂点1配列。
  !! @param[in] add_v2 追加する頂点2配列。
  subroutine append_triangles(v0, v1, v2, add_v0, add_v1, add_v2)
    real(dp), allocatable, intent(inout) :: v0(:, :), v1(:, :), v2(:, :)
    real(dp), intent(in) :: add_v0(:, :), add_v1(:, :), add_v2(:, :)
    real(dp), allocatable :: t0(:, :), t1(:, :), t2(:, :)
    integer :: n0, n1

    n0 = size(v0, 2)
    n1 = size(add_v0, 2)
    allocate (t0(3, n0 + n1), t1(3, n0 + n1), t2(3, n0 + n1))
    if (n0 > 0) then
      t0(:, 1:n0) = v0
      t1(:, 1:n0) = v1
      t2(:, 1:n0) = v2
    end if
    t0(:, n0 + 1:n0 + n1) = add_v0
    t1(:, n0 + 1:n0 + n1) = add_v1
    t2(:, n0 + 1:n0 + n1) = add_v2
    call move_alloc(t0, v0)
    call move_alloc(t1, v1)
    call move_alloc(t2, v2)
  end subroutine append_triangles

end module bem_app_config_runtime

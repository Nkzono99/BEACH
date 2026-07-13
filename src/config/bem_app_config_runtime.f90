!> `app_config` からメッシュ・粒子群を構築する実行時変換モジュール。
module bem_app_config_runtime
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann, eps0
  use bem_types, only: mesh_type, particles_soa, sim_config, injection_state, bc_periodic
  use bem_types, only: surface_model_insulator, surface_model_conductor, surface_model_dielectric
  use bem_mpi, only: mpi_context, mpi_get_rank_size, mpi_split_count, mpi_bcast_i32_array, mpi_bcast_real_dp_array
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_templates, only: make_plane, make_plate_hole, make_disk, make_annulus, make_box, make_cylinder, make_sphere
  use bem_mesh, only: init_mesh
  use bem_panel_surface_sides, only: resolve_panel_surface_sides, panel_surface_side_ok
  use bem_collision, only: collision_query_grid_stalled, collision_query_image_limit, &
                           collision_query_index_range, collision_query_invalid_segment, collision_query_ok
  use bem_importers, only: load_obj_mesh
  use bem_injection, only: &
    seed_rng, sample_uniform_positions, sample_shifted_maxwell_velocities, compute_macro_particles_for_batch, &
    compute_macro_particles_from_flux, sample_reservoir_face_particles, sample_reservoir_velocity_grid_particles, &
    sample_photo_raycast_particles, &
    compute_inflow_flux_from_drifting_maxwellian, compute_face_area_from_bounds
  use bem_particles, only: init_particles
  use bem_sheath_injection_model, only: sheath_injection_context, resolve_sheath_injection_context
  use bem_outer_plasma_types, only: outer_plasma_state_type
  use bem_app_config_types, only: &
    app_config, particle_species_spec, template_spec, particles_per_batch_from_config, &
    total_particles_from_config
  use bem_string_utils, only: lower_ascii
  use bem_config_helpers, only: resolve_inward_normal
  use, intrinsic :: iso_fortran_env, only: error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private :: finalize_particle_batch_collision_query

contains

  !> `mesh_mode` と OBJ ファイル有無に応じてメッシュ生成方法を選ぶ。
  !! @param[in] cfg メッシュ入力設定を含むアプリ設定。
  !! @param[out] mesh 構築した三角形メッシュ。
  subroutine build_mesh_from_config(cfg, mesh)
    type(app_config), intent(in) :: cfg
    type(mesh_type), intent(out) :: mesh
    logical :: has_obj, loaded_obj, need_transform

    loaded_obj = .false.
    select case (trim(lower_ascii(cfg%mesh_mode)))
    case ('obj')
      call load_obj_mesh(trim(cfg%obj_path), mesh)
      call apply_uniform_surface_model(mesh, surface_model_id_from_name(cfg%mesh_surface_model))
      call apply_uniform_epsilon_r(mesh, cfg%mesh_epsilon_r)
      loaded_obj = .true.
    case ('template')
      call build_template_mesh(cfg, mesh)
    case default
      inquire (file=trim(cfg%obj_path), exist=has_obj)
      if (has_obj) then
        call load_obj_mesh(trim(cfg%obj_path), mesh)
        call apply_uniform_surface_model(mesh, surface_model_id_from_name(cfg%mesh_surface_model))
        call apply_uniform_epsilon_r(mesh, cfg%mesh_epsilon_r)
        loaded_obj = .true.
      else
        call build_template_mesh(cfg, mesh)
      end if
    end select

    if (loaded_obj) then
      need_transform = (cfg%obj_scale /= 1.0d0) .or. &
                       any(cfg%obj_rotation /= 0.0d0) .or. &
                       any(cfg%obj_offset /= 0.0d0)
      if (need_transform) then
        call apply_obj_transform(mesh, cfg%obj_scale, cfg%obj_rotation, cfg%obj_offset)
      end if
    end if
    call apply_panel_surface_config(cfg, mesh, loaded_obj)
  end subroutine build_mesh_from_config

  subroutine apply_panel_surface_config(cfg, mesh, loaded_obj)
    type(app_config), intent(in) :: cfg
    type(mesh_type), intent(inout) :: mesh
    logical, intent(in) :: loaded_obj
    integer(i32) :: i, mesh_id, status
    character(len=256) :: message

    if (trim(lower_ascii(cfg%panel%source_model)) /= 'triangle_p0') return
    if (any(mesh%panel_area <= 0.0_dp)) then
      error stop 'triangle_p0 requires finite, non-degenerate triangles.'
    end if
    if (any(mesh%elem_surface_model /= surface_model_insulator)) then
      error stop 'triangle_p0 Phase 1 supports insulator surfaces only.'
    end if
    if (loaded_obj) then
      call resolve_panel_surface_sides(mesh, cfg%mesh_surface_side_policy, status, message)
      if (status /= panel_surface_side_ok) error stop 'mesh.surface_side: '//trim(message)
    else
      mesh_id = 0_i32
      do i = 1, cfg%n_templates
        if (.not. cfg%templates(i)%enabled) cycle
        mesh_id = mesh_id + 1_i32
        call resolve_panel_surface_sides( &
          mesh, cfg%templates(i)%surface_side_policy, status, message, mesh_id=mesh_id &
          )
        if (status /= panel_surface_side_ok) error stop 'mesh.templates.surface_side: '//trim(message)
      end do
    end if
    if (any(abs(mesh%elem_vacuum_sign) /= 1_i32)) then
      error stop 'triangle_p0 requires a resolved vacuum side for every element.'
    end if
  end subroutine apply_panel_surface_config

  !> OBJ メッシュの全頂点にスケール→回転→平行移動を適用し再初期化する。
  !! 変換順序: v_new = R(rotation) * (v_old * scale) + offset
  !! 回転は度単位で x→y→z の順に外因性 (extrinsic) 回転を適用する。
  !! @param[inout] mesh 変換対象の三角形メッシュ。
  !! @param[in] scale 一様スケーリング係数。
  !! @param[in] rotation_deg 回転角度 [rx, ry, rz] (度)。
  !! @param[in] offset 平行移動ベクトル (3成分)。
  subroutine apply_obj_transform(mesh, scale, rotation_deg, offset)
    type(mesh_type), intent(inout) :: mesh
    real(dp), intent(in) :: scale
    real(dp), intent(in) :: rotation_deg(3)
    real(dp), intent(in) :: offset(3)
    real(dp), parameter :: deg2rad = acos(-1.0d0)/180.0d0
    real(dp) :: rx, ry, rz, cx, sx, cy, sy, cz, sz
    real(dp) :: R(3, 3), v(3)
    real(dp), allocatable :: tv0(:, :), tv1(:, :), tv2(:, :)
    real(dp), allocatable :: elem_epsilon_r(:)
    integer(i32), allocatable :: elem_mesh_id(:), elem_surface_model(:)
    integer(i32) :: i, n

    rx = rotation_deg(1)*deg2rad
    ry = rotation_deg(2)*deg2rad
    rz = rotation_deg(3)*deg2rad
    cx = cos(rx); sx = sin(rx)
    cy = cos(ry); sy = sin(ry)
    cz = cos(rz); sz = sin(rz)

    ! R = Rz * Ry * Rx (extrinsic x→y→z)
    R(1, 1) = cy*cz; R(1, 2) = sx*sy*cz - cx*sz; R(1, 3) = cx*sy*cz + sx*sz
    R(2, 1) = cy*sz; R(2, 2) = sx*sy*sz + cx*cz; R(2, 3) = cx*sy*sz - sx*cz
    R(3, 1) = -sy; R(3, 2) = sx*cy; R(3, 3) = cx*cy

    n = mesh%nelem
    allocate (tv0(3, n), tv1(3, n), tv2(3, n))
    if (allocated(mesh%elem_mesh_id)) then
      allocate (elem_mesh_id(n))
      elem_mesh_id = mesh%elem_mesh_id
    end if
    if (allocated(mesh%elem_surface_model)) then
      allocate (elem_surface_model(n))
      elem_surface_model = mesh%elem_surface_model
    end if
    if (allocated(mesh%elem_epsilon_r)) then
      allocate (elem_epsilon_r(n))
      elem_epsilon_r = mesh%elem_epsilon_r
    end if
    do i = 1, n
      v = mesh%v0(:, i)*scale
      tv0(:, i) = matmul(R, v) + offset
      v = mesh%v1(:, i)*scale
      tv1(:, i) = matmul(R, v) + offset
      v = mesh%v2(:, i)*scale
      tv2(:, i) = matmul(R, v) + offset
    end do
    if (allocated(elem_mesh_id) .and. allocated(elem_surface_model) .and. allocated(elem_epsilon_r)) then
      call init_mesh( &
        mesh, tv0, tv1, tv2, elem_mesh_id0=elem_mesh_id, elem_surface_model0=elem_surface_model, &
        elem_epsilon_r0=elem_epsilon_r &
        )
    else if (allocated(elem_mesh_id)) then
      call init_mesh(mesh, tv0, tv1, tv2, elem_mesh_id0=elem_mesh_id)
    else if (allocated(elem_surface_model)) then
      call init_mesh(mesh, tv0, tv1, tv2, elem_surface_model0=elem_surface_model)
    else
      call init_mesh(mesh, tv0, tv1, tv2)
    end if
  end subroutine apply_obj_transform

  !> 全要素に同じ表面モデルIDを付与する。
  subroutine apply_uniform_surface_model(mesh, surface_model_id)
    type(mesh_type), intent(inout) :: mesh
    integer(i32), intent(in) :: surface_model_id

    if (.not. allocated(mesh%elem_surface_model)) return
    mesh%elem_surface_model = surface_model_id
  end subroutine apply_uniform_surface_model

  !> 全要素に同じ相対誘電率を付与する。
  subroutine apply_uniform_epsilon_r(mesh, epsilon_r)
    type(mesh_type), intent(inout) :: mesh
    real(dp), intent(in) :: epsilon_r

    if (.not. ieee_is_finite(epsilon_r) .or. epsilon_r < 1.0d0) error stop 'epsilon_r must be finite and >= 1.'
    if (.not. allocated(mesh%elem_epsilon_r)) return
    mesh%elem_epsilon_r = epsilon_r
  end subroutine apply_uniform_epsilon_r

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
        if (trim(lower_ascii(cfg%particle_species(s)%source_mode)) == 'reservoir_face' .or. &
            trim(lower_ascii(cfg%particle_species(s)%source_mode)) == 'photo_raycast') then
          error stop 'init_particles_from_config supports volume_seed only. Use init_particle_batch_from_config.'
        end if
        if (cfg%particle_species(s)%npcls_per_step < 0_i32) then
          error stop 'particles.species.npcls_per_step must be >= 0.'
        end if
        counts(s) = cfg%sim%batch_count*cfg%particle_species(s)%npcls_per_step
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
  subroutine seed_particles_from_config(cfg, mpi_rank, mpi_size, mpi)
    type(app_config), intent(in) :: cfg
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    type(mpi_context), intent(in), optional :: mpi
    integer(i32) :: local_rank, n_ranks, seed_value
    integer(kind=8) :: seed_tmp

    call resolve_parallel_rank_size(local_rank, n_ranks, mpi_rank, mpi_size, mpi, 'seed_particles_from_config')

    seed_tmp = int(cfg%sim%rng_seed, kind=8) + 104729_8*int(local_rank, kind=8)
    seed_value = int(modulo(seed_tmp, int(huge(0_i32), kind=8)), kind=i32)
    call seed_rng([seed_value])
  end subroutine seed_particles_from_config

  !> 指定バッチ番号に対応する粒子バッチを生成する。
  !! @param[in] cfg 粒子種とシミュレーション条件を含むアプリ設定。
  !! @param[in] batch_idx 生成対象のバッチ番号（1始まり）。
  !! @param[out] pcls 生成したバッチ粒子群。
  !! @param[inout] state reservoir_face 注入の残差状態（必要時のみ）。
  !! @param[in] mesh 現在バッチ開始時点の電荷分布メッシュ（電位補正時に必要）。
  !! @param[inout] snapshot refresh 済み静電 snapshot（電位障壁 closure 使用時に必要）。
  !! @param[in] outer_state 現在バッチの外部プラズマ状態（kinetic 1D流入写像時に必要）。
  !! @param[out] photo_emission_dq photo_raycast 放出起因の要素電荷差分 `photo_emission_dq(nelem)`（省略可）。
  !! @param[out] collision_failure_status 不完全な photo collision query の status（省略時は停止）。
  !! @param[out] collision_failure_species 不完全な照会を返した最小 species index。
  !! @param[out] collision_failure_ray 不完全な照会を返した最小 ray index。
  !! @param[out] collision_failure_bounce 不完全な照会を返した bounce index。
  subroutine init_particle_batch_from_config( &
    cfg, batch_idx, pcls, state, mesh, snapshot, outer_state, photo_emission_dq, mpi_rank, mpi_size, mpi, &
    collision_failure_status, collision_failure_species, collision_failure_ray, collision_failure_bounce &
    )
    type(app_config), intent(in) :: cfg
    integer(i32), intent(in) :: batch_idx
    type(particles_soa), intent(out) :: pcls
    type(injection_state), intent(inout), optional :: state
    type(mesh_type), intent(in), optional :: mesh
    type(electrostatic_snapshot_type), intent(inout), optional :: snapshot
    type(outer_plasma_state_type), intent(in), optional :: outer_state
    real(dp), intent(out), optional :: photo_emission_dq(:)
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    type(mpi_context), intent(in), optional :: mpi
    integer(i32), intent(out), optional :: collision_failure_status, collision_failure_species
    integer(i32), intent(out), optional :: collision_failure_ray, collision_failure_bounce

    integer(i32) :: s, i, batch_n, max_rank, out_idx, local_rank, n_ranks, global_count
    integer(i32) :: photo_collision_status, photo_collision_ray, photo_collision_bounce
    integer(i32), allocatable :: counts_max(:), counts_actual(:), global_counts(:), species_cursor(:), species_id(:), &
                                 emit_elem_species(:, :)
    real(dp), allocatable :: vmin_normal(:), barrier_normal(:), effective_density_m3(:), w_effective(:), &
                             effective_particle_flux_m2_s(:), &
                             effective_temperature_k(:), effective_drift_velocity(:, :), &
                             photo_emit_current_density(:), photo_vmin_normal(:), photo_normal_drift_speed(:), &
                             photo_escape_factor_cache(:)
    logical, allocatable :: apply_barrier_energy_shift(:)
    logical :: has_enabled_reservoir, use_collective_reservoir_count
    real(dp), allocatable :: x_species(:, :, :), v_species(:, :, :), w_species(:, :)
    real(dp), allocatable :: x(:, :), v(:, :), q(:), m(:), w(:)
    type(sheath_injection_context) :: sheath_ctx

    if (present(collision_failure_status)) collision_failure_status = collision_query_ok
    if (present(collision_failure_species)) collision_failure_species = huge(0_i32)
    if (present(collision_failure_ray)) collision_failure_ray = huge(0_i32)
    if (present(collision_failure_bounce)) collision_failure_bounce = 0_i32
    if (cfg%sim%batch_count <= 0_i32) error stop 'sim.batch_count must be > 0.'
    if (batch_idx < 1_i32 .or. batch_idx > cfg%sim%batch_count) then
      error stop 'Requested batch index is out of range.'
    end if
    call resolve_parallel_rank_size(local_rank, n_ranks, mpi_rank, mpi_size, mpi, 'init_particle_batch_from_config')
    has_enabled_reservoir = .false.
    do s = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(s)%enabled) cycle
      has_enabled_reservoir = has_enabled_reservoir .or. &
                              trim(lower_ascii(cfg%particle_species(s)%source_mode)) == 'reservoir_face'
    end do
    use_collective_reservoir_count = present(mpi) .and. has_enabled_reservoir
    if (present(state)) then
      if (.not. allocated(state%macro_residual)) error stop 'injection_state is not initialized.'
      if (size(state%macro_residual) < cfg%n_particle_species) error stop 'injection_state size mismatch.'
    end if
    if (present(photo_emission_dq)) then
      if (.not. present(mesh)) error stop 'photo_emission_dq requires mesh in init_particle_batch_from_config.'
      if (size(photo_emission_dq) /= mesh%nelem) error stop 'photo_emission_dq size mismatch.'
      photo_emission_dq = 0.0d0
    end if

    allocate (counts_max(cfg%n_particle_species), counts_actual(cfg%n_particle_species), global_counts(cfg%n_particle_species))
    allocate (vmin_normal(cfg%n_particle_species), barrier_normal(cfg%n_particle_species))
    allocate (effective_density_m3(cfg%n_particle_species), w_effective(cfg%n_particle_species))
    allocate (effective_particle_flux_m2_s(cfg%n_particle_species))
    allocate (effective_temperature_k(cfg%n_particle_species), effective_drift_velocity(3, cfg%n_particle_species))
    allocate (photo_emit_current_density(cfg%n_particle_species), photo_vmin_normal(cfg%n_particle_species))
    allocate (photo_normal_drift_speed(cfg%n_particle_species), apply_barrier_energy_shift(cfg%n_particle_species))
    counts_max = 0_i32
    counts_actual = 0_i32
    global_counts = 0_i32
    vmin_normal = 0.0d0
    barrier_normal = 0.0d0
    effective_density_m3 = 0.0d0
    w_effective = 0.0d0
    effective_particle_flux_m2_s = 0.0d0
    effective_temperature_k = 0.0d0
    effective_drift_velocity = 0.0d0
    photo_emit_current_density = 0.0d0
    photo_vmin_normal = 0.0d0
    photo_normal_drift_speed = 0.0d0
    apply_barrier_energy_shift = .true.
    do s = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(s)%enabled) cycle
      select case (trim(lower_ascii(cfg%particle_species(s)%source_mode)))
      case ('volume_seed')
        w_effective(s) = cfg%particle_species(s)%w_particle
        effective_temperature_k(s) = species_temperature_k(cfg%particle_species(s))
        effective_drift_velocity(:, s) = cfg%particle_species(s)%drift_velocity
      case ('reservoir_face')
        if (trim(lower_ascii(cfg%particle_species(s)%velocity_distribution)) == 'grid') then
          effective_particle_flux_m2_s(s) = cfg%particle_species(s)%particle_flux_m2_s
        else
          effective_density_m3(s) = species_number_density_m3(cfg%particle_species(s))
        end if
        w_effective(s) = cfg%particle_species(s)%w_particle
        effective_temperature_k(s) = species_temperature_k(cfg%particle_species(s))
        effective_drift_velocity(:, s) = cfg%particle_species(s)%drift_velocity
      case ('photo_raycast')
        photo_emit_current_density(s) = cfg%particle_species(s)%emit_current_density_a_m2
        photo_normal_drift_speed(s) = cfg%particle_species(s)%normal_drift_speed
      end select
    end do

    call resolve_sheath_injection_context(cfg, sheath_ctx)
    if (sheath_ctx%enabled) then
      s = int(sheath_ctx%electron_species)
      effective_density_m3(s) = sheath_ctx%electron_number_density_m3
      vmin_normal(s) = max(vmin_normal(s), sheath_ctx%electron_vmin_normal)
      apply_barrier_energy_shift(s) = .false.

      if (cfg%particle_species(s)%has_target_macro_particles_per_batch .and. &
          cfg%particle_species(s)%target_macro_particles_per_batch > 0_i32) then
        call resolve_reservoir_target_weight( &
          cfg%sim, cfg%particle_species(s), effective_density_m3(s), vmin_normal(s), &
          effective_temperature_k(s), effective_drift_velocity(:, s), &
          cfg%particle_species(s)%target_macro_particles_per_batch, w_effective(s) &
          )
      end if

      if (sheath_ctx%has_local_reservoir_profile) then
        s = int(sheath_ctx%ion_species)
        effective_density_m3(s) = sheath_ctx%ion_number_density_m3
        effective_temperature_k(s) = 0.0d0
        call apply_normal_speed_override( &
          cfg%particle_species(s)%drift_velocity, sheath_ctx%reference_inward_normal, &
          sheath_ctx%ion_normal_speed_mps, effective_drift_velocity(:, s) &
          )

        if (cfg%particle_species(s)%has_target_macro_particles_per_batch .and. &
            cfg%particle_species(s)%target_macro_particles_per_batch > 0_i32) then
          call resolve_reservoir_target_weight( &
            cfg%sim, cfg%particle_species(s), effective_density_m3(s), vmin_normal(s), &
            effective_temperature_k(s), effective_drift_velocity(:, s), &
            cfg%particle_species(s)%target_macro_particles_per_batch, w_effective(s) &
            )
        end if
      end if

      if (sheath_ctx%has_photo_species) then
        s = int(sheath_ctx%photo_species)
        photo_emit_current_density(s) = sheath_ctx%photo_emit_current_density_a_m2
        photo_vmin_normal(s) = sheath_ctx%photo_vmin_normal
        photo_normal_drift_speed(s) = 0.0d0
      end if

      do s = 1, cfg%n_particle_species
        if (.not. cfg%particle_species(s)%enabled) cycle
        if (trim(lower_ascii(cfg%particle_species(s)%source_mode)) /= 'reservoir_face') cycle
        if (.not. cfg%particle_species(s)%has_target_macro_particles_per_batch) cycle
        if (cfg%particle_species(s)%target_macro_particles_per_batch == -1_i32) then
          w_effective(s) = w_effective(1)
        end if
      end do
    end if
    do s = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(s)%enabled) cycle
      select case (trim(lower_ascii(cfg%particle_species(s)%source_mode)))
      case ('volume_seed')
        global_count = cfg%particle_species(s)%npcls_per_step
        counts_max(s) = mpi_split_count(global_count, local_rank, n_ranks)
      case ('reservoir_face')
        if (.not. present(state)) then
          error stop 'reservoir_face requires injection_state in init_particle_batch_from_config.'
        end if
        call reservoir_face_velocity_correction( &
          cfg, cfg%particle_species(s), vmin_normal(s), barrier_normal(s), mesh, snapshot, outer_state &
          )
        if (.not. use_collective_reservoir_count .or. local_rank == 0_i32) then
          call compute_macro_particles_for_species( &
            cfg%sim, cfg%particle_species(s), state%macro_residual(s), global_counts(s), vmin_normal=vmin_normal(s), &
            number_density_override=effective_density_m3(s), particle_flux_override=effective_particle_flux_m2_s(s), &
            temperature_k_override=effective_temperature_k(s), drift_velocity_override=effective_drift_velocity(:, s), &
            w_particle_override=w_effective(s) &
            )
        end if
        if (.not. use_collective_reservoir_count) then
          counts_max(s) = mpi_split_count(global_counts(s), local_rank, n_ranks)
        end if
      case ('photo_raycast')
        global_count = cfg%particle_species(s)%rays_per_batch
        counts_max(s) = mpi_split_count(global_count, local_rank, n_ranks)
      case default
        error stop 'Unknown particles.species.source_mode.'
      end select
    end do
    if (use_collective_reservoir_count) then
      call mpi_bcast_i32_array(mpi, global_counts, 0_i32)
      call mpi_bcast_real_dp_array(mpi, state%macro_residual, 0_i32)
      do s = 1, cfg%n_particle_species
        if (.not. cfg%particle_species(s)%enabled) cycle
        if (trim(lower_ascii(cfg%particle_species(s)%source_mode)) /= 'reservoir_face') cycle
        counts_max(s) = mpi_split_count(global_counts(s), local_rank, n_ranks)
      end do
    end if
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
      select case (trim(lower_ascii(cfg%particle_species(s)%source_mode)))
      case ('volume_seed', 'reservoir_face')
        call sample_species_state( &
          cfg%sim, cfg%particle_species(s), counts_max(s), &
          x_species(:, 1:counts_max(s), s), v_species(:, 1:counts_max(s), s), &
          barrier_normal_energy=barrier_normal(s), vmin_normal=vmin_normal(s), &
          apply_barrier_energy_shift=apply_barrier_energy_shift(s), &
          temperature_k_override=effective_temperature_k(s), drift_velocity_override=effective_drift_velocity(:, s) &
          )
        counts_actual(s) = counts_max(s)
        w_species(1:counts_actual(s), s) = w_effective(s)
      case ('photo_raycast')
        if (.not. present(mesh)) then
          error stop 'photo_raycast requires mesh in init_particle_batch_from_config.'
        end if
        if (photo_emit_current_density(s) <= 0.0d0) then
          counts_actual(s) = 0_i32
          cycle
        end if
        call sample_photo_species_state( &
          cfg%sim, cfg%particle_species(s), mesh, counts_max(s), x_species(:, 1:counts_max(s), s), &
          v_species(:, 1:counts_max(s), s), w_species(1:counts_max(s), s), counts_actual(s), &
          emit_elem_idx=emit_elem_species(1:counts_max(s), s), &
          global_rays_per_batch=cfg%particle_species(s)%rays_per_batch, &
          emit_current_density_override=photo_emit_current_density(s), &
          normal_drift_speed_override=photo_normal_drift_speed(s), vmin_normal=photo_vmin_normal(s), &
          collision_failure_status=photo_collision_status, collision_failure_ray=photo_collision_ray, &
          collision_failure_bounce=photo_collision_bounce &
          )
        if (photo_collision_status /= collision_query_ok) then
          call finalize_particle_batch_collision_query( &
            photo_collision_status, batch_idx, s, photo_collision_ray, photo_collision_bounce, &
            collision_failure_status, collision_failure_species, collision_failure_ray, collision_failure_bounce &
            )
          return
        end if
        if (photo_escape_model_enabled(cfg%particle_species(s))) then
          if (.not. present(snapshot)) then
            error stop 'photo_escape_model requires a refreshed electrostatic snapshot.'
          end if
          allocate (photo_escape_factor_cache(mesh%nelem))
          photo_escape_factor_cache = -1.0d0
          do i = 1, counts_actual(s)
            if (emit_elem_species(i, s) < 1_i32 .or. emit_elem_species(i, s) > mesh%nelem) then
              error stop 'photo_raycast emitted invalid elem_idx.'
            end if
            if (photo_escape_factor_cache(emit_elem_species(i, s)) < 0.0d0) then
              photo_escape_factor_cache(emit_elem_species(i, s)) = photo_escape_weight_factor( &
                                                                   mesh, cfg%sim, cfg%particle_species(s), &
                                                                   emit_elem_species(i, s), &
                                                                   snapshot &
                                                                   )
            end if
            w_species(i, s) = w_species(i, s)*photo_escape_factor_cache(emit_elem_species(i, s))
          end do
          deallocate (photo_escape_factor_cache)
        end if
        if (present(photo_emission_dq) .and. cfg%particle_species(s)%deposit_opposite_charge_on_emit) then
          do i = 1, counts_actual(s)
            if (emit_elem_species(i, s) < 1_i32 .or. emit_elem_species(i, s) > size(photo_emission_dq)) then
              error stop 'photo_raycast emitted invalid elem_idx.'
            end if
            photo_emission_dq(emit_elem_species(i, s)) = photo_emission_dq(emit_elem_species(i, s)) - &
                                                         cfg%particle_species(s)%q_particle*w_species(i, s)
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

    call init_particles(pcls, x, v, q, m, w, species_id=species_id)
  end subroutine init_particle_batch_from_config

  !> batch injection の不完全な photo collision query を返し、status 未要求なら serial に停止する。
  subroutine finalize_particle_batch_collision_query( &
    query_status, batch_idx, species_idx, query_ray, query_bounce, &
    collision_failure_status, collision_failure_species, collision_failure_ray, collision_failure_bounce &
    )
    integer(i32), intent(in) :: query_status, batch_idx, species_idx, query_ray, query_bounce
    integer(i32), intent(out), optional :: collision_failure_status, collision_failure_species
    integer(i32), intent(out), optional :: collision_failure_ray, collision_failure_bounce
    character(len=16) :: status_name
    character(len=256) :: failure_message

    if (present(collision_failure_status)) collision_failure_status = query_status
    if (present(collision_failure_species)) collision_failure_species = species_idx
    if (present(collision_failure_ray)) collision_failure_ray = query_ray
    if (present(collision_failure_bounce)) collision_failure_bounce = query_bounce
    if (present(collision_failure_status)) return

    select case (query_status)
    case (collision_query_image_limit)
      status_name = 'image_limit'
    case (collision_query_index_range)
      status_name = 'index_range'
    case (collision_query_invalid_segment)
      status_name = 'invalid_segment'
    case (collision_query_grid_stalled)
      status_name = 'grid_stalled'
    case default
      status_name = 'unknown'
    end select
    write (failure_message, '(a,i0,a,i0,a,i0,a,i0,a,a,a,i0)') &
      'photo_raycast collision query incomplete during batch preparation: batch=', batch_idx, &
      ' species=', species_idx, ' ray=', query_ray, ' bounce=', query_bounce, &
      ' status=', trim(status_name), ' code=', query_status
    write (error_unit, '(a)') trim(failure_message)
    flush (error_unit)
    error stop 1
  end subroutine finalize_particle_batch_collision_query

  !> 1粒子種ぶんの位置・速度サンプルをまとめて生成する。
  !! @param[in] sim ボックス境界・バッチ時間などのシミュレーション設定。
  !! @param[in] spec 1粒子種の注入設定。
  !! @param[in] n 生成粒子数。
  !! @param[out] x 生成した位置配列 `x(3,n)` [m]。
  !! @param[out] v 生成した速度配列 `v(3,n)` [m/s]。
  !! @param[in] barrier_normal_energy reservoir_face 法線方向のエネルギー障壁 `2 q Δφ / m` [`m^2/s^2`]。
  !! @param[in] vmin_normal reservoir_face 法線速度の下限 [m/s]。
  !! @param[in] apply_barrier_energy_shift reservoir_face 法線速度へ障壁エネルギー変換を適用するか。
  subroutine sample_species_state( &
    sim, spec, n, x, v, barrier_normal_energy, vmin_normal, apply_barrier_energy_shift, &
    temperature_k_override, drift_velocity_override &
    )
    type(sim_config), intent(in) :: sim
    type(particle_species_spec), intent(in) :: spec
    integer(i32), intent(in) :: n
    real(dp), intent(out) :: x(:, :)
    real(dp), intent(out) :: v(:, :)
    real(dp), intent(in), optional :: barrier_normal_energy
    real(dp), intent(in), optional :: vmin_normal
    logical, intent(in), optional :: apply_barrier_energy_shift
    real(dp), intent(in), optional :: temperature_k_override
    real(dp), intent(in), optional :: drift_velocity_override(3)
    logical :: apply_shift
    real(dp) :: temperature_k_local, drift_velocity_local(3)

    if (n <= 0_i32) return
    apply_shift = .true.
    if (present(apply_barrier_energy_shift)) apply_shift = apply_barrier_energy_shift
    temperature_k_local = species_temperature_k(spec)
    if (present(temperature_k_override)) temperature_k_local = temperature_k_override
    drift_velocity_local = spec%drift_velocity
    if (present(drift_velocity_override)) drift_velocity_local = drift_velocity_override
    select case (trim(lower_ascii(spec%source_mode)))
    case ('volume_seed')
      call sample_uniform_positions(spec%pos_low, spec%pos_high, x)
      call sample_shifted_maxwell_velocities( &
        drift_velocity_local, spec%m_particle, v, temperature_k=temperature_k_local &
        )
    case ('reservoir_face')
      if (trim(lower_ascii(spec%velocity_distribution)) == 'grid') then
        if (present(barrier_normal_energy) .and. present(vmin_normal)) then
          call sample_reservoir_velocity_grid_particles( &
            sim%box_min, sim%box_max, spec%inject_face, spec%pos_low, spec%pos_high, &
            spec%velocity_grid_path, spec%velocity_grid_pdf_kind, sim%batch_duration, x, v, &
            barrier_normal_energy=barrier_normal_energy, vmin_normal=vmin_normal, position_jitter_dt=sim%dt, &
            apply_barrier_energy_shift=apply_shift, velocity_grid_sampling=spec%velocity_grid_sampling &
            )
        else if (present(barrier_normal_energy)) then
          call sample_reservoir_velocity_grid_particles( &
            sim%box_min, sim%box_max, spec%inject_face, spec%pos_low, spec%pos_high, &
            spec%velocity_grid_path, spec%velocity_grid_pdf_kind, sim%batch_duration, x, v, &
            barrier_normal_energy=barrier_normal_energy, position_jitter_dt=sim%dt, apply_barrier_energy_shift=apply_shift, &
            velocity_grid_sampling=spec%velocity_grid_sampling &
            )
        else if (present(vmin_normal)) then
          call sample_reservoir_velocity_grid_particles( &
            sim%box_min, sim%box_max, spec%inject_face, spec%pos_low, spec%pos_high, &
            spec%velocity_grid_path, spec%velocity_grid_pdf_kind, sim%batch_duration, x, v, &
            vmin_normal=vmin_normal, position_jitter_dt=sim%dt, apply_barrier_energy_shift=apply_shift, &
            velocity_grid_sampling=spec%velocity_grid_sampling &
            )
        else
          call sample_reservoir_velocity_grid_particles( &
            sim%box_min, sim%box_max, spec%inject_face, spec%pos_low, spec%pos_high, &
            spec%velocity_grid_path, spec%velocity_grid_pdf_kind, sim%batch_duration, x, v, position_jitter_dt=sim%dt, &
            apply_barrier_energy_shift=apply_shift, velocity_grid_sampling=spec%velocity_grid_sampling &
            )
        end if
      else if (present(barrier_normal_energy) .and. present(vmin_normal)) then
        call sample_reservoir_face_particles( &
          sim%box_min, sim%box_max, spec%inject_face, spec%pos_low, spec%pos_high, drift_velocity_local, &
          spec%m_particle, temperature_k_local, sim%batch_duration, x, v, &
          barrier_normal_energy=barrier_normal_energy, vmin_normal=vmin_normal, position_jitter_dt=sim%dt, &
          apply_barrier_energy_shift=apply_shift &
          )
      else if (present(barrier_normal_energy)) then
        call sample_reservoir_face_particles( &
          sim%box_min, sim%box_max, spec%inject_face, spec%pos_low, spec%pos_high, drift_velocity_local, &
          spec%m_particle, temperature_k_local, sim%batch_duration, x, v, &
          barrier_normal_energy=barrier_normal_energy, position_jitter_dt=sim%dt, &
          apply_barrier_energy_shift=apply_shift &
          )
      else if (present(vmin_normal)) then
        call sample_reservoir_face_particles( &
          sim%box_min, sim%box_max, spec%inject_face, spec%pos_low, spec%pos_high, drift_velocity_local, &
          spec%m_particle, temperature_k_local, sim%batch_duration, x, v, &
          vmin_normal=vmin_normal, position_jitter_dt=sim%dt, apply_barrier_energy_shift=apply_shift &
          )
      else
        call sample_reservoir_face_particles( &
          sim%box_min, sim%box_max, spec%inject_face, spec%pos_low, spec%pos_high, drift_velocity_local, &
          spec%m_particle, temperature_k_local, sim%batch_duration, x, v, position_jitter_dt=sim%dt, &
          apply_barrier_energy_shift=apply_shift &
          )
      end if
    case ('photo_raycast')
      error stop 'sample_species_state does not support photo_raycast. Use sample_photo_species_state.'
    case default
      error stop 'Unknown particles.species.source_mode.'
    end select
    if (trim(lower_ascii(spec%source_mode)) == 'reservoir_face') call normalize_reservoir_positions(sim, x)
  end subroutine sample_species_state

  !> reservoirの時間ジッタ後の位置を、設定されたbox境界条件に従って有効領域へ戻す。
  pure subroutine normalize_reservoir_positions(sim, x)
    type(sim_config), intent(in) :: sim
    real(dp), intent(inout) :: x(:, :)
    integer(i32) :: axis
    real(dp) :: span

    if (.not. sim%use_box) return
    do axis = 1_i32, 3_i32
      if (.not. ieee_is_finite(sim%box_min(axis)) .or. .not. ieee_is_finite(sim%box_max(axis))) cycle
      span = sim%box_max(axis) - sim%box_min(axis)
      if (.not. ieee_is_finite(span) .or. span <= 0.0_dp) cycle
      if (sim%bc_low(axis) == bc_periodic .and. sim%bc_high(axis) == bc_periodic) then
        x(axis, :) = sim%box_min(axis) + modulo(x(axis, :) - sim%box_min(axis), span)
      else
        x(axis, :) = min(max(x(axis, :), sim%box_min(axis)), sim%box_max(axis))
      end if
    end do
  end subroutine normalize_reservoir_positions

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
  subroutine sample_photo_species_state( &
    sim, spec, mesh, n_rays, x, v, w, n_emit, emit_elem_idx, global_rays_per_batch, &
    emit_current_density_override, normal_drift_speed_override, vmin_normal, &
    collision_failure_status, collision_failure_ray, collision_failure_bounce &
    )
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
    real(dp), intent(in), optional :: emit_current_density_override
    real(dp), intent(in), optional :: normal_drift_speed_override
    real(dp), intent(in), optional :: vmin_normal
    integer(i32), intent(out), optional :: collision_failure_status, collision_failure_ray, collision_failure_bounce
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
  end subroutine sample_photo_species_state

  !> reservoir_face 用に、物理流量と残差から今バッチのマクロ粒子数を決める。
  !! @param[in] sim ボックス境界・バッチ時間などのシミュレーション設定。
  !! @param[in] spec reservoir_face 粒子種設定。
  !! @param[inout] residual 前バッチから繰り越した端数。
  !! @param[out] count 今バッチで生成するマクロ粒子数。
  !! @param[in] vmin_normal 法線速度の下限 [m/s]（省略時は 0）。
  !! @param[in] number_density_override 数密度の上書き値 [1/m^3]。
  !! @param[in] particle_flux_override 粒子数 flux の上書き値 [1/m^2/s]。
  !! @param[in] w_particle_override マクロ粒子重みの上書き値。
  !! @param[in] temperature_k_override 温度の上書き値 [K]。
  !! @param[in] drift_velocity_override ドリフト速度の上書き値 [m/s]。
  subroutine compute_macro_particles_for_species( &
    sim, spec, residual, count, vmin_normal, number_density_override, w_particle_override, &
    temperature_k_override, drift_velocity_override, particle_flux_override &
    )
    type(sim_config), intent(in) :: sim
    type(particle_species_spec), intent(in) :: spec
    real(dp), intent(inout) :: residual
    integer(i32), intent(out) :: count
    real(dp), intent(in), optional :: vmin_normal
    real(dp), intent(in), optional :: number_density_override
    real(dp), intent(in), optional :: w_particle_override
    real(dp), intent(in), optional :: temperature_k_override
    real(dp), intent(in), optional :: drift_velocity_override(3)
    real(dp), intent(in), optional :: particle_flux_override

    real(dp) :: number_density_m3, effective_batch_duration, particle_flux_m2_s, w_particle, temperature_k_local
    real(dp) :: drift_velocity_local(3)

    number_density_m3 = species_number_density_m3(spec)
    if (present(number_density_override)) number_density_m3 = number_density_override
    w_particle = spec%w_particle
    if (present(w_particle_override)) w_particle = w_particle_override
    temperature_k_local = species_temperature_k(spec)
    if (present(temperature_k_override)) temperature_k_local = temperature_k_override
    drift_velocity_local = spec%drift_velocity
    if (present(drift_velocity_override)) drift_velocity_local = drift_velocity_override
    effective_batch_duration = sim%batch_duration
    if (trim(lower_ascii(spec%velocity_distribution)) == 'grid') then
      particle_flux_m2_s = spec%particle_flux_m2_s
      if (present(particle_flux_override)) particle_flux_m2_s = particle_flux_override
      call compute_macro_particles_from_flux( &
        particle_flux_m2_s, spec%inject_face, spec%pos_low, spec%pos_high, effective_batch_duration, w_particle, residual, count &
        )
      return
    end if
    if (present(vmin_normal)) then
      call compute_macro_particles_for_batch( &
        number_density_m3, temperature_k_local, spec%m_particle, drift_velocity_local, sim%box_min, sim%box_max, &
        spec%inject_face, spec%pos_low, spec%pos_high, effective_batch_duration, w_particle, residual, count, &
        vmin_normal=vmin_normal &
        )
    else
      call compute_macro_particles_for_batch( &
        number_density_m3, temperature_k_local, spec%m_particle, drift_velocity_local, sim%box_min, sim%box_max, &
        spec%inject_face, spec%pos_low, spec%pos_high, effective_batch_duration, w_particle, residual, count &
        )
    end if
  end subroutine compute_macro_particles_for_species

  !> reservoir_face の target 個数からシース補正込み重みを解決する。
  !! @param[in] cfg シミュレーション・結合設定を含むアプリ設定。
  !! @param[in] spec reservoir_face 粒子種設定。
  !! @param[in] number_density_m3 実効数密度 [1/m^3]。
  !! @param[in] vmin_normal 法線速度の下限 [m/s]。
  !! @param[in] temperature_k 実効温度 [K]。
  !! @param[in] drift_velocity 実効ドリフト速度 [m/s]。
  !! @param[in] target_macro_particles_per_batch 目標マクロ粒子数。
  !! @param[out] w_particle 解決したマクロ粒子重み。
  subroutine resolve_reservoir_target_weight( &
    sim, spec, number_density_m3, vmin_normal, temperature_k, drift_velocity, target_macro_particles_per_batch, w_particle &
    )
    type(sim_config), intent(in) :: sim
    type(particle_species_spec), intent(in) :: spec
    real(dp), intent(in) :: number_density_m3, vmin_normal, temperature_k, drift_velocity(3)
    integer(i32), intent(in) :: target_macro_particles_per_batch
    real(dp), intent(out) :: w_particle

    real(dp) :: inward_normal(3), gamma_in, area

    if (target_macro_particles_per_batch <= 0_i32) then
      error stop 'resolve_reservoir_target_weight requires target_macro_particles_per_batch > 0.'
    end if
    call resolve_inward_normal(spec%inject_face, inward_normal)
    gamma_in = compute_inflow_flux_from_drifting_maxwellian( &
               number_density_m3, temperature_k, spec%m_particle, drift_velocity, inward_normal, &
               vmin_normal=vmin_normal &
               )
    area = compute_face_area_from_bounds(spec%inject_face, spec%pos_low, spec%pos_high)
    w_particle = gamma_in*area*sim%batch_duration/real(target_macro_particles_per_batch, dp)
    if (.not. ieee_is_finite(w_particle) .or. w_particle <= 0.0d0) then
      error stop 'sheath-adjusted target_macro_particles_per_batch produced invalid w_particle.'
    end if
  end subroutine resolve_reservoir_target_weight

  subroutine apply_normal_speed_override(drift_velocity, inward_normal, normal_speed, drift_velocity_out)
    real(dp), intent(in) :: drift_velocity(3), inward_normal(3), normal_speed
    real(dp), intent(out) :: drift_velocity_out(3)

    drift_velocity_out = drift_velocity - dot_product(drift_velocity, inward_normal)*inward_normal + normal_speed*inward_normal
  end subroutine apply_normal_speed_override

  !> reservoir_face 注入に対する法線速度補正パラメータを計算する。
  !! @param[in] cfg シミュレーション・結合設定を含むアプリ設定。
  !! @param[in] spec reservoir_face 粒子種設定。
  !! @param[out] vmin_normal 無限遠法線速度の下限 [m/s]。
  !! @param[out] barrier_normal 法線エネルギー障壁 `2 q Δφ / m` [`m^2/s^2`]。
  !! @param[in] mesh 現在バッチ開始時点の電荷分布メッシュ（補正時に必要）。
  !! @param[inout] snapshot refresh 済み静電 snapshot（legacy infinity barrier 使用時に必要）。
  subroutine reservoir_face_velocity_correction(cfg, spec, vmin_normal, barrier_normal, mesh, snapshot, outer_state)
    type(app_config), intent(in) :: cfg
    type(particle_species_spec), intent(in) :: spec
    real(dp), intent(out) :: vmin_normal
    real(dp), intent(out) :: barrier_normal
    type(mesh_type), intent(in), optional :: mesh
    type(electrostatic_snapshot_type), intent(inout), optional :: snapshot
    type(outer_plasma_state_type), intent(in), optional :: outer_state

    real(dp) :: phi_face, delta_phi, area_xy

    vmin_normal = 0.0d0
    barrier_normal = 0.0d0
    if (trim(lower_ascii(cfg%coupling%particle_transfer_mode)) == 'electrostatic_1d_instant_return') then
      if (trim(lower_ascii(spec%inject_face)) /= 'z_high') then
        error stop 'Split outer-plasma ambient inflow must use inject_face="z_high".'
      end if
      if (trim(lower_ascii(cfg%outer_plasma%model)) == 'kinetic_1d') then
        if (.not. present(outer_state)) error stop 'kinetic_1d ambient inflow requires the refreshed outer state.'
        if (.not. outer_state%ready .or. trim(lower_ascii(outer_state%model)) /= 'kinetic_1d') then
          error stop 'kinetic_1d ambient inflow requires a ready kinetic outer state.'
        end if
        delta_phi = outer_state%interface_potential - outer_state%infinity_potential
        if (.not. ieee_is_finite(delta_phi)) error stop 'kinetic_1d ambient inflow barrier is non-finite.'
        barrier_normal = 2.0_dp*spec%q_particle*delta_phi/spec%m_particle
        if (.not. ieee_is_finite(barrier_normal)) error stop 'kinetic_1d ambient inflow map is non-finite.'
        if (barrier_normal > 0.0_dp) vmin_normal = sqrt(barrier_normal)
        return
      end if
      if (.not. present(mesh)) error stop 'Split outer-plasma ambient inflow requires the charged mesh.'
      area_xy = (cfg%sim%box_max(1) - cfg%sim%box_min(1))*(cfg%sim%box_max(2) - cfg%sim%box_min(2))
      delta_phi = cfg%outer_plasma%debye_length*sum(mesh%q_elem)/(eps0*area_xy)
      barrier_normal = 2.0_dp*spec%q_particle*delta_phi/spec%m_particle
      if (.not. ieee_is_finite(barrier_normal)) error stop 'Outer-plasma inflow map produced a non-finite barrier.'
      if (barrier_normal > 0.0_dp) vmin_normal = sqrt(barrier_normal)
      return
    end if
    select case (trim(lower_ascii(cfg%sim%reservoir_potential_model)))
    case ('none')
      return
    case ('infinity_barrier')
      if (.not. present(mesh)) then
        error stop 'sim.reservoir_potential_model="infinity_barrier" requires mesh in init_particle_batch_from_config.'
      end if
      if (.not. present(snapshot)) then
        error stop 'sim.reservoir_potential_model="infinity_barrier" requires a refreshed electrostatic snapshot.'
      end if
      call compute_face_average_potential(mesh, cfg%sim, spec, snapshot, phi_face)
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
      error stop 'Unknown sim.reservoir_potential_model in runtime.'
    end select
  end subroutine reservoir_face_velocity_correction

  !> photo_raycast のPE escape closureが有効かを返す。
  pure logical function photo_escape_model_enabled(spec) result(enabled)
    type(particle_species_spec), intent(in) :: spec

    enabled = trim(lower_ascii(spec%photo_escape_model)) /= 'none'
  end function photo_escape_model_enabled

  !> 放出元要素の局所障壁からPE escape重み係数を返す。
  real(dp) function photo_escape_weight_factor(mesh, sim, spec, elem_idx, snapshot) result(factor)
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(particle_species_spec), intent(in) :: spec
    integer(i32), intent(in) :: elem_idx
    type(electrostatic_snapshot_type), intent(inout) :: snapshot

    real(dp) :: element_potential, barrier_v, thermal_energy_j, arg

    select case (trim(lower_ascii(spec%photo_escape_model)))
    case ('none')
      factor = 1.0d0
      return
    case ('boltzmann_cutoff')
      call snapshot%eval_local_phi_without_primary_self(mesh, sim, elem_idx, element_potential)
      barrier_v = max(element_potential - sim%phi_infty, 0.0d0)
      if (barrier_v <= 0.0d0) then
        factor = 1.0d0
        return
      end if
      thermal_energy_j = k_boltzmann*species_temperature_k(spec)
      if (thermal_energy_j <= 0.0d0) then
        factor = 0.0d0
        return
      end if
      arg = abs(spec%q_particle)*barrier_v/thermal_energy_j
      if (.not. ieee_is_finite(arg)) then
        error stop 'photo_escape_model produced a non-finite Boltzmann argument.'
      end if
      if (arg >= 700.0d0) then
        factor = 0.0d0
      else
        factor = exp(-arg)
      end if
    case default
      error stop 'Unknown photo_escape_model in runtime.'
    end select
  end function photo_escape_weight_factor

  !> reservoir_face 開口面の平均電位を `N x N` 格子平均で評価する。
  !! @param[in] mesh 現在バッチ開始時点の電荷分布メッシュ。
  !! @param[in] sim シミュレーション設定。
  !! @param[in] spec reservoir_face 粒子種設定。
  !! @param[inout] snapshot refresh 済み静電 snapshot。
  !! @param[out] phi_face 注入開口面の平均電位 [V]。
  subroutine compute_face_average_potential(mesh, sim, spec, snapshot, phi_face)
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(particle_species_spec), intent(in) :: spec
    type(electrostatic_snapshot_type), intent(inout) :: snapshot
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
      t1 = (real(i, dp) - 0.5d0)/real(ngrid, dp)
      do j = 1_i32, ngrid
        t2 = (real(j, dp) - 0.5d0)/real(ngrid, dp)
        pos = 0.0d0
        pos(axis_n) = boundary_value
        pos(axis_t1) = spec%pos_low(axis_t1) + (spec%pos_high(axis_t1) - spec%pos_low(axis_t1))*t1
        pos(axis_t2) = spec%pos_low(axis_t2) + (spec%pos_high(axis_t2) - spec%pos_low(axis_t2))*t2
        pos = pos + inward_normal*1.0d-12
        call snapshot%eval_local_phi(mesh, sim, pos, phi)
        phi_sum = phi_sum + phi
      end do
    end do

    phi_face = phi_sum/real(ngrid*ngrid, dp)
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

  !> 粒子種設定から実効密度[m^-3]を返す。
  !! @param[in] spec 粒子種設定。
  !! @return number_density_m3 実効粒子数密度 [1/m^3]。
  pure real(dp) function species_number_density_m3(spec) result(number_density_m3)
    type(particle_species_spec), intent(in) :: spec

    number_density_m3 = spec%number_density_m3
    if (spec%has_number_density_cm3) number_density_m3 = spec%number_density_cm3*1.0d6
  end function species_number_density_m3

  !> 併存対応のため `mpi_context` と rank/size の両方を受け、最終的なrank/sizeを解決する。
  subroutine resolve_parallel_rank_size(local_rank, n_ranks, mpi_rank, mpi_size, mpi, caller_name)
    integer(i32), intent(out) :: local_rank, n_ranks
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    type(mpi_context), intent(in), optional :: mpi
    character(len=*), intent(in) :: caller_name

    call mpi_get_rank_size(local_rank, n_ranks, mpi)
    if (present(mpi_rank)) local_rank = mpi_rank
    if (present(mpi_size)) n_ranks = mpi_size
    if (n_ranks <= 0_i32) error stop 'mpi_size must be > 0 in '//trim(caller_name)//'.'
    if (local_rank < 0_i32 .or. local_rank >= n_ranks) then
      error stop 'mpi_rank is out of range in '//trim(caller_name)//'.'
    end if
  end subroutine resolve_parallel_rank_size

  !> 粒子種設定から実効温度[K]を返す。
  !! @param[in] spec 粒子種設定。
  !! @return temperature_k 実効温度 [K]。
  pure real(dp) function species_temperature_k(spec) result(temperature_k)
    type(particle_species_spec), intent(in) :: spec

    temperature_k = spec%temperature_k
    if (spec%has_temperature_ev) temperature_k = spec%temperature_ev*1.160451812d4
  end function species_temperature_k

  !> 有効なテンプレートを連結し、1つのメッシュへまとめる。
  !! @param[in] cfg テンプレート設定を含むアプリ設定。
  !! @param[out] mesh 連結後の三角形メッシュ。
  subroutine build_template_mesh(cfg, mesh)
    type(app_config), intent(in) :: cfg
    type(mesh_type), intent(out) :: mesh
    type(mesh_type) :: part
    real(dp), allocatable :: v0(:, :), v1(:, :), v2(:, :)
    real(dp), allocatable :: elem_epsilon_r(:), part_epsilon_r(:)
    integer(i32), allocatable :: elem_mesh_id(:), part_mesh_id(:)
    integer(i32), allocatable :: elem_surface_model(:), part_surface_model(:)
    integer(i32) :: i, mesh_id

    allocate (v0(3, 0), v1(3, 0), v2(3, 0), elem_mesh_id(0), elem_surface_model(0), elem_epsilon_r(0))
    mesh_id = 0_i32
    if (.not. allocated(cfg%templates)) then
      error stop 'Template storage is not allocated in configuration.'
    end if
    if (cfg%n_templates > int(size(cfg%templates), i32)) then
      error stop 'Template count exceeds allocated storage.'
    end if
    do i = 1, cfg%n_templates
      if (.not. cfg%templates(i)%enabled) cycle
      mesh_id = mesh_id + 1_i32
      call build_one_template(cfg%templates(i), part)
      call append_triangles(v0, v1, v2, part%v0, part%v1, part%v2)
      if (allocated(part_mesh_id)) deallocate (part_mesh_id)
      allocate (part_mesh_id(part%nelem))
      part_mesh_id = mesh_id
      call append_mesh_ids(elem_mesh_id, part_mesh_id)
      if (allocated(part_surface_model)) deallocate (part_surface_model)
      allocate (part_surface_model(part%nelem))
      part_surface_model = surface_model_id_from_name(cfg%templates(i)%surface_model)
      call append_mesh_ids(elem_surface_model, part_surface_model)
      if (allocated(part_epsilon_r)) deallocate (part_epsilon_r)
      allocate (part_epsilon_r(part%nelem))
      part_epsilon_r = cfg%templates(i)%epsilon_r
      call append_real_values(elem_epsilon_r, part_epsilon_r)
    end do

    if (size(v0, 2) == 0) then
      error stop 'No enabled template found in configuration.'
    end if
    call init_mesh( &
      mesh, v0, v1, v2, elem_mesh_id0=elem_mesh_id, elem_surface_model0=elem_surface_model, &
      elem_epsilon_r0=elem_epsilon_r &
      )
  end subroutine build_template_mesh

  !> テンプレート種別に応じて形状生成ルーチンへディスパッチする。
  !! @param[in] spec 1テンプレート分の形状設定。
  !! @param[out] mesh 生成したテンプレートメッシュ。
  subroutine build_one_template(spec, mesh)
    type(template_spec), intent(in) :: spec
    type(mesh_type), intent(out) :: mesh
    logical :: cap_top, cap_bottom

    select case (trim(lower_ascii(spec%kind)))
    case ('plane')
      call make_plane(mesh, size_x=spec%size_x, size_y=spec%size_y, nx=spec%nx, ny=spec%ny, center=spec%center)
    case ('plate_hole', 'plane_hole')
      call make_plate_hole( &
        mesh, size_x=spec%size_x, size_y=spec%size_y, radius=spec%radius, n_theta=spec%n_theta, n_r=spec%n_r, &
        center=spec%center &
        )
    case ('disk')
      call make_disk(mesh, radius=spec%radius, n_theta=spec%n_theta, n_r=spec%n_r, center=spec%center)
    case ('annulus')
      call make_annulus( &
        mesh, radius=spec%radius, inner_radius=spec%inner_radius, n_theta=spec%n_theta, n_r=spec%n_r, center=spec%center &
        )
    case ('box')
      call make_box(mesh, size=spec%size, center=spec%center, nx=spec%nx, ny=spec%ny, nz=spec%nz)
    case ('cylinder')
      cap_top = spec%cap
      cap_bottom = spec%cap
      if (spec%has_cap_top) cap_top = spec%cap_top
      if (spec%has_cap_bottom) cap_bottom = spec%cap_bottom
      call make_cylinder( &
        mesh, radius=spec%radius, height=spec%height, n_theta=spec%n_theta, n_z=spec%n_z, &
        cap=spec%cap, center=spec%center, cap_top=cap_top, cap_bottom=cap_bottom &
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

  !> 既存の要素メッシュID配列へ追加分を連結する。
  !! @param[inout] mesh_ids 累積した要素メッシュID配列。
  !! @param[in] add_ids 追加する要素メッシュID配列。
  subroutine append_mesh_ids(mesh_ids, add_ids)
    integer(i32), allocatable, intent(inout) :: mesh_ids(:)
    integer(i32), intent(in) :: add_ids(:)
    integer(i32), allocatable :: tmp(:)
    integer(i32) :: n0, n1

    n0 = size(mesh_ids)
    n1 = size(add_ids)
    allocate (tmp(n0 + n1))
    if (n0 > 0) tmp(1:n0) = mesh_ids
    if (n1 > 0) tmp(n0 + 1:n0 + n1) = add_ids
    call move_alloc(tmp, mesh_ids)
  end subroutine append_mesh_ids

  !> 既存の実数要素配列へ追加分を連結する。
  subroutine append_real_values(values, add_values)
    real(dp), allocatable, intent(inout) :: values(:)
    real(dp), intent(in) :: add_values(:)
    real(dp), allocatable :: tmp(:)
    integer(i32) :: n0, n1

    n0 = size(values)
    n1 = size(add_values)
    allocate (tmp(n0 + n1))
    if (n0 > 0) tmp(1:n0) = values
    if (n1 > 0) tmp(n0 + 1:n0 + n1) = add_values
    call move_alloc(tmp, values)
  end subroutine append_real_values

  !> 表面モデル名をメッシュ要素用の整数IDへ変換する。
  !! @param[in] name `insulator` / `conductor` / `dielectric`。
  !! @return model_id 内部表面モデルID。
  integer(i32) function surface_model_id_from_name(name) result(model_id)
    character(len=*), intent(in) :: name

    select case (trim(lower_ascii(name)))
    case ('insulator')
      model_id = surface_model_insulator
    case ('conductor')
      model_id = surface_model_conductor
    case ('dielectric')
      model_id = surface_model_dielectric
    case default
      error stop 'Unknown surface_model.'
    end select
  end function surface_model_id_from_name

end module bem_app_config_runtime

!> `app_config` からメッシュ・粒子群を構築する実行時変換モジュール。
module bem_app_config_runtime
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
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
  use bem_external_boundary_contract, only: &
    external_boundary_contract_type, external_boundary_ok, external_inflow_none, external_inflow_scalar_barrier, &
    resolve_external_boundary_contract
  use bem_app_config_types, only: &
    app_config, particle_species_spec, template_spec, particles_per_batch_from_config, &
    total_particles_from_config, particle_inflow_reservoir
  use bem_string_utils, only: lower_ascii
  use bem_config_helpers, only: resolve_inward_normal
  use, intrinsic :: iso_fortran_env, only: error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private :: finalize_particle_batch_collision_query
  private :: compute_rectangular_face_potential_statistics

  !> 1回の simulation run で不変な粒子 source の導出値。
  type, public :: particle_source_plan_type
    private
    logical :: ready = .false.
    logical :: mpi_argument_present = .false.
    logical :: use_collective_reservoir_count = .false.
    integer(i32) :: nspecies = 0_i32
    integer(i32) :: mpi_rank = 0_i32
    integer(i32) :: mpi_size = 1_i32
    real(dp), allocatable :: effective_density_m3(:)
    real(dp), allocatable :: effective_particle_flux_m2_s(:)
    real(dp), allocatable :: effective_temperature_k(:)
    real(dp), allocatable :: effective_drift_velocity(:, :)
    real(dp), allocatable :: effective_weight(:)
    real(dp), allocatable :: photo_emit_current_density(:)
    real(dp), allocatable :: photo_normal_drift_speed(:)
  end type particle_source_plan_type

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

    if (any(mesh%panel_area <= 0.0_dp)) then
      error stop 'triangle_p0 requires finite, non-degenerate triangles.'
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
            trim(lower_ascii(cfg%particle_species(s)%source_mode)) == 'plane_source' .or. &
            trim(lower_ascii(cfg%particle_species(s)%source_mode)) == 'photo_raycast' .or. &
            has_boundary_inflow(cfg%particle_species(s))) then
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

  !> 設定とMPI配置だけに依存する粒子 source の導出値を構築する。
  !! 乱数、残差、mesh/snapshot依存の障壁は扱わず、run中に不変な値だけを保持する。
  subroutine build_particle_source_plan(cfg, plan, mpi_rank, mpi_size, mpi)
    type(app_config), intent(in) :: cfg
    type(particle_source_plan_type), intent(out) :: plan
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    type(mpi_context), intent(in), optional :: mpi

    integer(i32) :: s, local_rank, n_ranks
    logical :: has_enabled_reservoir

    call resolve_parallel_rank_size(local_rank, n_ranks, mpi_rank, mpi_size, mpi, 'build_particle_source_plan')
    plan%nspecies = cfg%n_particle_species
    plan%mpi_rank = local_rank
    plan%mpi_size = n_ranks
    plan%mpi_argument_present = present(mpi)
    has_enabled_reservoir = .false.
    do s = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(s)%enabled) cycle
      has_enabled_reservoir = has_enabled_reservoir .or. &
                              trim(lower_ascii(cfg%particle_species(s)%source_mode)) == 'reservoir_face' .or. &
                              trim(lower_ascii(cfg%particle_species(s)%source_mode)) == 'plane_source' .or. &
                              has_boundary_inflow(cfg%particle_species(s))
    end do
    plan%use_collective_reservoir_count = present(mpi) .and. has_enabled_reservoir

    allocate (plan%effective_density_m3(cfg%n_particle_species))
    allocate (plan%effective_particle_flux_m2_s(cfg%n_particle_species))
    allocate (plan%effective_temperature_k(cfg%n_particle_species))
    allocate (plan%effective_drift_velocity(3, cfg%n_particle_species))
    allocate (plan%effective_weight(cfg%n_particle_species))
    allocate (plan%photo_emit_current_density(cfg%n_particle_species))
    allocate (plan%photo_normal_drift_speed(cfg%n_particle_species))
    plan%effective_density_m3 = 0.0_dp
    plan%effective_particle_flux_m2_s = 0.0_dp
    plan%effective_temperature_k = 0.0_dp
    plan%effective_drift_velocity = 0.0_dp
    plan%effective_weight = 0.0_dp
    plan%photo_emit_current_density = 0.0_dp
    plan%photo_normal_drift_speed = 0.0_dp

    do s = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(s)%enabled) cycle
      select case (trim(lower_ascii(cfg%particle_species(s)%source_mode)))
      case ('volume_seed')
        plan%effective_weight(s) = cfg%particle_species(s)%w_particle
        plan%effective_temperature_k(s) = species_temperature_k(cfg%particle_species(s))
        plan%effective_drift_velocity(:, s) = cfg%particle_species(s)%drift_velocity
        if (has_boundary_inflow(cfg%particle_species(s))) then
          if (trim(lower_ascii(cfg%particle_species(s)%velocity_distribution)) == 'grid') then
            plan%effective_particle_flux_m2_s(s) = cfg%particle_species(s)%particle_flux_m2_s
          else
            plan%effective_density_m3(s) = species_number_density_m3(cfg%particle_species(s))
          end if
        end if
      case ('reservoir_face', 'plane_source')
        if (trim(lower_ascii(cfg%particle_species(s)%velocity_distribution)) == 'grid') then
          plan%effective_particle_flux_m2_s(s) = cfg%particle_species(s)%particle_flux_m2_s
        else
          plan%effective_density_m3(s) = species_number_density_m3(cfg%particle_species(s))
        end if
        plan%effective_weight(s) = cfg%particle_species(s)%w_particle
        plan%effective_temperature_k(s) = species_temperature_k(cfg%particle_species(s))
        plan%effective_drift_velocity(:, s) = cfg%particle_species(s)%drift_velocity
      case ('photo_raycast')
        plan%photo_emit_current_density(s) = cfg%particle_species(s)%emit_current_density_a_m2
        plan%photo_normal_drift_speed(s) = cfg%particle_species(s)%normal_drift_speed
      end select
    end do

    plan%ready = .true.
  end subroutine build_particle_source_plan

  !> 指定バッチ番号に対応する粒子バッチを生成する。
  !! @param[in] cfg 粒子種とシミュレーション条件を含むアプリ設定。
  !! @param[in] batch_idx 生成対象のバッチ番号（1始まり）。
  !! @param[out] pcls 生成したバッチ粒子群。
  !! @param[inout] state reservoir_face 注入の残差状態（必要時のみ）。
  !! @param[in] mesh 現在バッチ開始時点の電荷分布メッシュ（電位補正時に必要）。
  !! @param[out] photo_emission_dq photo_raycast 放出起因の要素電荷差分 `photo_emission_dq(nelem)`（省略可）。
  !! @param[out] collision_failure_status 不完全な photo collision query の status（省略時は停止）。
  !! @param[out] collision_failure_species 不完全な照会を返した最小 species index。
  !! @param[out] collision_failure_ray 不完全な照会を返した最小 ray index。
  !! @param[out] collision_failure_bounce 不完全な照会を返した bounce index。
  !! @param[inout] snapshot refresh 済み静電 snapshot（注入電位補正の使用時に必要）。
  !! @param[in] source_plan run中に再利用する粒子source導出値（省略時は呼出し内で構築）。
  !! @param[out] photo_emission_dq_by_species species別のphoto放出反作用電荷 `photo_emission_dq_by_species(nelem, nspecies)`（省略可）。
  subroutine init_particle_batch_from_config( &
    cfg, batch_idx, pcls, state, mesh, photo_emission_dq, mpi_rank, mpi_size, mpi, &
    collision_failure_status, collision_failure_species, collision_failure_ray, collision_failure_bounce, snapshot, &
    source_plan, photo_emission_dq_by_species &
    )
    type(app_config), intent(in) :: cfg
    integer(i32), intent(in) :: batch_idx
    type(particles_soa), intent(out) :: pcls
    type(injection_state), intent(inout), optional :: state
    type(mesh_type), intent(in), optional :: mesh
    real(dp), intent(out), optional :: photo_emission_dq(:)
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    type(mpi_context), intent(in), optional :: mpi
    integer(i32), intent(out), optional :: collision_failure_status, collision_failure_species
    integer(i32), intent(out), optional :: collision_failure_ray, collision_failure_bounce
    type(electrostatic_snapshot_type), intent(inout), optional :: snapshot
    type(particle_source_plan_type), intent(in), optional, target :: source_plan
    real(dp), intent(out), optional :: photo_emission_dq_by_species(:, :)

    integer(i32) :: s, i, face, batch_n, max_rank, out_idx, local_rank, n_ranks, global_count
    integer(i32) :: source_begin, source_end, face_begin, face_end
    integer(i32) :: boundary_status
    integer(i32) :: photo_collision_status, photo_collision_ray, photo_collision_bounce
    integer(i32), allocatable :: counts_max(:), counts_actual(:), source_counts(:), global_counts(:), &
                                 boundary_counts(:, :), boundary_global_counts(:, :), species_cursor(:), species_id(:), &
                                 source_element(:), emit_elem_species(:, :)
    real(dp), allocatable :: vmin_normal(:), barrier_normal(:), boundary_vmin(:, :), boundary_barrier(:, :), &
                             batch_density_m3(:), batch_weight(:)
    logical :: use_collective_reservoir_count
    real(dp), allocatable :: x_species(:, :, :), v_species(:, :, :), w_species(:, :)
    real(dp), allocatable :: x(:, :), v(:, :), q(:), m(:), w(:)
    type(particle_source_plan_type), target :: generated_source_plan
    type(particle_source_plan_type), pointer :: active_source_plan
    type(external_boundary_contract_type) :: active_boundary_contract
    type(particle_species_spec) :: face_spec
    real(dp) :: correction_vmin_normal
    character(len=256) :: boundary_message

    if (present(collision_failure_status)) collision_failure_status = collision_query_ok
    if (present(collision_failure_species)) collision_failure_species = huge(0_i32)
    if (present(collision_failure_ray)) collision_failure_ray = huge(0_i32)
    if (present(collision_failure_bounce)) collision_failure_bounce = 0_i32
    if (cfg%sim%batch_count <= 0_i32) error stop 'sim.batch_count must be > 0.'
    if (batch_idx < 1_i32 .or. batch_idx > cfg%sim%batch_count) then
      error stop 'Requested batch index is out of range.'
    end if
    call resolve_external_boundary_contract( &
      cfg%sim%reservoir_potential_model, cfg%sim%open_boundary_model, &
      active_boundary_contract, boundary_status, boundary_message &
      )
    if (boundary_status /= external_boundary_ok) error stop trim(boundary_message)
    call resolve_parallel_rank_size(local_rank, n_ranks, mpi_rank, mpi_size, mpi, 'init_particle_batch_from_config')
    if (present(source_plan)) then
      active_source_plan => source_plan
    else
      call build_particle_source_plan( &
        cfg, generated_source_plan, mpi_rank=mpi_rank, mpi_size=mpi_size, mpi=mpi &
        )
      active_source_plan => generated_source_plan
    end if
    if (.not. active_source_plan%ready) error stop 'particle source plan is not initialized.'
    if (active_source_plan%nspecies /= cfg%n_particle_species) then
      error stop 'particle source plan species count does not match app config.'
    end if
    if (active_source_plan%mpi_rank /= local_rank .or. active_source_plan%mpi_size /= n_ranks .or. &
        active_source_plan%mpi_argument_present .neqv. present(mpi)) then
      error stop 'particle source plan MPI context does not match batch initialization.'
    end if
    use_collective_reservoir_count = active_source_plan%use_collective_reservoir_count
    if (present(state)) then
      if (.not. allocated(state%macro_residual)) error stop 'injection_state is not initialized.'
      if (size(state%macro_residual) < cfg%n_particle_species) error stop 'injection_state size mismatch.'
      if (any([(has_boundary_inflow(cfg%particle_species(s)), s=1, cfg%n_particle_species)])) then
        if (.not. allocated(state%boundary_macro_residual)) then
          error stop 'boundary_inflow requires boundary_macro_residual in injection_state.'
        end if
        if (size(state%boundary_macro_residual, 1) < 6 .or. &
            size(state%boundary_macro_residual, 2) < cfg%n_particle_species) then
          error stop 'injection_state boundary residual size mismatch.'
        end if
      end if
    end if
    if (present(photo_emission_dq)) then
      if (.not. present(mesh)) error stop 'photo_emission_dq requires mesh in init_particle_batch_from_config.'
      if (size(photo_emission_dq) /= mesh%nelem) error stop 'photo_emission_dq size mismatch.'
      photo_emission_dq = 0.0d0
    end if
    if (present(photo_emission_dq_by_species)) then
      if (.not. present(mesh)) then
        error stop 'photo_emission_dq_by_species requires mesh in init_particle_batch_from_config.'
      end if
      if (size(photo_emission_dq_by_species, 1) /= mesh%nelem .or. &
          size(photo_emission_dq_by_species, 2) /= cfg%n_particle_species) then
        error stop 'photo_emission_dq_by_species size mismatch.'
      end if
      photo_emission_dq_by_species = 0.0_dp
    end if

    allocate ( &
      counts_max(cfg%n_particle_species), counts_actual(cfg%n_particle_species), source_counts(cfg%n_particle_species), &
      global_counts(cfg%n_particle_species) &
      )
    allocate (boundary_counts(6, cfg%n_particle_species), boundary_global_counts(6, cfg%n_particle_species))
    allocate (vmin_normal(cfg%n_particle_species), barrier_normal(cfg%n_particle_species))
    allocate (boundary_vmin(6, cfg%n_particle_species), boundary_barrier(6, cfg%n_particle_species))
    allocate (batch_density_m3(cfg%n_particle_species), batch_weight(cfg%n_particle_species))
    counts_max = 0_i32
    counts_actual = 0_i32
    source_counts = 0_i32
    global_counts = 0_i32
    boundary_counts = 0_i32
    boundary_global_counts = 0_i32
    vmin_normal = 0.0d0
    barrier_normal = 0.0d0
    boundary_vmin = 0.0_dp
    boundary_barrier = 0.0_dp
    batch_density_m3 = active_source_plan%effective_density_m3
    batch_weight = active_source_plan%effective_weight
    associate ( &
      effective_particle_flux_m2_s => active_source_plan%effective_particle_flux_m2_s, &
      effective_temperature_k => active_source_plan%effective_temperature_k, &
      effective_drift_velocity => active_source_plan%effective_drift_velocity, &
      photo_emit_current_density => active_source_plan%photo_emit_current_density, &
      photo_normal_drift_speed => active_source_plan%photo_normal_drift_speed &
      )
    do s = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(s)%enabled) cycle
      select case (trim(lower_ascii(cfg%particle_species(s)%source_mode)))
      case ('volume_seed')
        global_count = cfg%particle_species(s)%npcls_per_step
        source_counts(s) = mpi_split_count(global_count, local_rank, n_ranks)
      case ('reservoir_face', 'plane_source')
        if (.not. present(state)) then
          error stop 'flux-driven source requires injection_state in init_particle_batch_from_config.'
        end if
        if (trim(lower_ascii(cfg%particle_species(s)%source_mode)) == 'reservoir_face') then
          call reservoir_face_velocity_correction( &
            cfg, cfg%particle_species(s), correction_vmin_normal, barrier_normal(s), mesh, snapshot, &
            warn_face_variation=local_rank == 0_i32 .and. &
            (batch_idx == 1_i32 .or. batch_idx == cfg%sim%batch_count), &
            boundary_contract=active_boundary_contract &
            )
          vmin_normal(s) = max(vmin_normal(s), correction_vmin_normal)
        end if
        if (.not. use_collective_reservoir_count .or. local_rank == 0_i32) then
          call compute_macro_particles_for_species( &
            cfg%sim, cfg%particle_species(s), state%macro_residual(s), global_counts(s), vmin_normal=vmin_normal(s), &
            number_density_override=batch_density_m3(s), particle_flux_override=effective_particle_flux_m2_s(s), &
            temperature_k_override=effective_temperature_k(s), drift_velocity_override=effective_drift_velocity(:, s), &
            w_particle_override=batch_weight(s) &
            )
        end if
        if (.not. use_collective_reservoir_count) then
          source_counts(s) = mpi_split_count(global_counts(s), local_rank, n_ranks)
        end if
      case ('photo_raycast')
        global_count = cfg%particle_species(s)%rays_per_batch
        source_counts(s) = mpi_split_count(global_count, local_rank, n_ranks)
      case default
        error stop 'Unknown particles.species.source_mode.'
      end select
      if (has_boundary_inflow(cfg%particle_species(s))) then
        if (.not. present(state)) then
          error stop 'boundary_inflow requires injection_state in init_particle_batch_from_config.'
        end if
        do face = 1, 6
          if (.not. boundary_inflow_face_enabled(cfg%particle_species(s), face)) cycle
          call make_boundary_inflow_spec(cfg%sim, cfg%particle_species(s), face, face_spec)
          call reservoir_face_velocity_correction( &
            cfg, face_spec, correction_vmin_normal, boundary_barrier(face, s), mesh, snapshot, &
            warn_face_variation=local_rank == 0_i32 .and. &
            (batch_idx == 1_i32 .or. batch_idx == cfg%sim%batch_count), &
            boundary_contract=active_boundary_contract &
            )
          boundary_vmin(face, s) = correction_vmin_normal
          if (.not. use_collective_reservoir_count .or. local_rank == 0_i32) then
            call compute_macro_particles_for_species( &
              cfg%sim, face_spec, state%boundary_macro_residual(face, s), boundary_global_counts(face, s), &
              vmin_normal=boundary_vmin(face, s), number_density_override=batch_density_m3(s), &
              particle_flux_override=effective_particle_flux_m2_s(s), &
              temperature_k_override=effective_temperature_k(s), &
              drift_velocity_override=effective_drift_velocity(:, s), w_particle_override=batch_weight(s) &
              )
          end if
          if (.not. use_collective_reservoir_count) then
            boundary_counts(face, s) = mpi_split_count(boundary_global_counts(face, s), local_rank, n_ranks)
          end if
        end do
      end if
    end do
    if (use_collective_reservoir_count) then
      call mpi_bcast_i32_array(mpi, global_counts, 0_i32)
      call mpi_bcast_real_dp_array(mpi, state%macro_residual, 0_i32)
      do s = 1, cfg%n_particle_species
        if (.not. cfg%particle_species(s)%enabled) cycle
        if (trim(lower_ascii(cfg%particle_species(s)%source_mode)) == 'reservoir_face' .or. &
            trim(lower_ascii(cfg%particle_species(s)%source_mode)) == 'plane_source') then
          source_counts(s) = mpi_split_count(global_counts(s), local_rank, n_ranks)
        end if
        if (.not. has_boundary_inflow(cfg%particle_species(s))) cycle
        call mpi_bcast_i32_array(mpi, boundary_global_counts(:, s), 0_i32)
        call mpi_bcast_real_dp_array(mpi, state%boundary_macro_residual(:, s), 0_i32)
        do face = 1, 6
          boundary_counts(face, s) = mpi_split_count(boundary_global_counts(face, s), local_rank, n_ranks)
        end do
      end do
    end if
    counts_max = source_counts + sum(boundary_counts, dim=1)
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
      source_begin = 1_i32
      source_end = source_counts(s)
      if (source_end >= source_begin) then
        select case (trim(lower_ascii(cfg%particle_species(s)%source_mode)))
        case ('volume_seed', 'reservoir_face', 'plane_source')
          call sample_species_state( &
            cfg%sim, cfg%particle_species(s), source_counts(s), &
            x_species(:, source_begin:source_end, s), v_species(:, source_begin:source_end, s), &
            barrier_normal_energy=barrier_normal(s), vmin_normal=vmin_normal(s), &
            temperature_k_override=effective_temperature_k(s), drift_velocity_override=effective_drift_velocity(:, s) &
            )
          counts_actual(s) = source_counts(s)
          w_species(source_begin:source_end, s) = batch_weight(s)
        case ('photo_raycast')
          if (.not. present(mesh)) then
            error stop 'photo_raycast requires mesh in init_particle_batch_from_config.'
          end if
          if (photo_emit_current_density(s) > 0.0_dp) then
            call sample_photo_species_state( &
              cfg%sim, cfg%particle_species(s), mesh, source_counts(s), x_species(:, source_begin:source_end, s), &
              v_species(:, source_begin:source_end, s), w_species(source_begin:source_end, s), counts_actual(s), &
              emit_elem_idx=emit_elem_species(source_begin:source_end, s), &
              global_rays_per_batch=cfg%particle_species(s)%rays_per_batch, &
              emit_current_density_override=photo_emit_current_density(s), &
              normal_drift_speed_override=photo_normal_drift_speed(s), &
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
            if ((present(photo_emission_dq) .or. present(photo_emission_dq_by_species)) .and. &
                cfg%particle_species(s)%deposit_opposite_charge_on_emit) then
              do i = 1, counts_actual(s)
                if (emit_elem_species(i, s) < 1_i32 .or. emit_elem_species(i, s) > mesh%nelem) then
                  error stop 'photo_raycast emitted invalid elem_idx.'
                end if
                if (present(photo_emission_dq)) then
                  photo_emission_dq(emit_elem_species(i, s)) = photo_emission_dq(emit_elem_species(i, s)) - &
                                                               cfg%particle_species(s)%q_particle*w_species(i, s)
                end if
                if (present(photo_emission_dq_by_species)) then
                  photo_emission_dq_by_species(emit_elem_species(i, s), s) = &
                    photo_emission_dq_by_species(emit_elem_species(i, s), s) - &
                    cfg%particle_species(s)%q_particle*w_species(i, s)
                end if
              end do
            end if
          end if
        case default
          error stop 'Unknown particles.species.source_mode.'
        end select
      end if
      do face = 1, 6
        if (boundary_counts(face, s) <= 0_i32) cycle
        face_begin = counts_actual(s) + 1_i32
        face_end = face_begin + boundary_counts(face, s) - 1_i32
        call make_boundary_inflow_spec(cfg%sim, cfg%particle_species(s), face, face_spec)
        call sample_species_state( &
          cfg%sim, face_spec, boundary_counts(face, s), x_species(:, face_begin:face_end, s), &
          v_species(:, face_begin:face_end, s), barrier_normal_energy=boundary_barrier(face, s), &
          vmin_normal=boundary_vmin(face, s), temperature_k_override=effective_temperature_k(s), &
          drift_velocity_override=effective_drift_velocity(:, s) &
          )
        w_species(face_begin:face_end, s) = batch_weight(s)
        counts_actual(s) = face_end
      end do
    end do

    batch_n = sum(counts_actual)
    allocate (species_id(batch_n), source_element(batch_n))
    source_element = -1_i32
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
      if (trim(lower_ascii(cfg%particle_species(s)%source_mode)) == 'photo_raycast') then
        source_element(i) = emit_elem_species(species_cursor(s), s)
      end if
    end do

    call init_particles( &
      pcls, x, v, q, m, w, species_id=species_id, source_element=source_element &
      )
    end associate
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

  !> speciesに外部 reservoir 流入を指定したbox面があるかを返す。
  pure logical function has_boundary_inflow(spec) result(has_inflow)
    type(particle_species_spec), intent(in) :: spec

    has_inflow = any(spec%boundary_inflow_low == particle_inflow_reservoir) .or. &
                 any(spec%boundary_inflow_high == particle_inflow_reservoir)
  end function has_boundary_inflow

  !> face bit順の面に reservoir 流入が有効かを返す。
  pure logical function boundary_inflow_face_enabled(spec, face) result(enabled)
    type(particle_species_spec), intent(in) :: spec
    integer(i32), intent(in) :: face

    if (mod(face, 2_i32) == 1_i32) then
      enabled = spec%boundary_inflow_low((face + 1_i32)/2_i32) == particle_inflow_reservoir
    else
      enabled = spec%boundary_inflow_high(face/2_i32) == particle_inflow_reservoir
    end if
  end function boundary_inflow_face_enabled

  !> box全面を開口とする一時的な legacy reservoir spec を構築する。
  pure subroutine make_boundary_inflow_spec(sim, source_spec, face, inflow_spec)
    type(sim_config), intent(in) :: sim
    type(particle_species_spec), intent(in) :: source_spec
    integer(i32), intent(in) :: face
    type(particle_species_spec), intent(out) :: inflow_spec
    integer(i32) :: axis
    logical :: high_side

    inflow_spec = source_spec
    inflow_spec%source_mode = 'reservoir_face'
    call boundary_face_name(face, inflow_spec%inject_face)
    inflow_spec%pos_low = sim%box_min
    inflow_spec%pos_high = sim%box_max
    axis = (face + 1_i32)/2_i32
    high_side = mod(face, 2_i32) == 0_i32
    if (high_side) then
      inflow_spec%pos_low(axis) = sim%box_max(axis)
      inflow_spec%pos_high(axis) = sim%box_max(axis)
    else
      inflow_spec%pos_low(axis) = sim%box_min(axis)
      inflow_spec%pos_high(axis) = sim%box_min(axis)
    end if
    inflow_spec%has_npcls_per_step = .false.
    inflow_spec%has_source_normal = .false.
    inflow_spec%boundary_inflow_low = 0_i32
    inflow_spec%boundary_inflow_high = 0_i32
  end subroutine make_boundary_inflow_spec

  !> plane_sourceの内部面をlegacy face samplerの仮想box境界へ写像する。
  pure subroutine configure_plane_source_box(spec, box_min, box_max)
    type(particle_species_spec), intent(in) :: spec
    real(dp), intent(inout) :: box_min(3), box_max(3)
    integer(i32) :: axis

    select case (trim(lower_ascii(spec%inject_face)))
    case ('x_low')
      axis = 1_i32
      box_min(axis) = spec%pos_low(axis)
    case ('x_high')
      axis = 1_i32
      box_max(axis) = spec%pos_low(axis)
    case ('y_low')
      axis = 2_i32
      box_min(axis) = spec%pos_low(axis)
    case ('y_high')
      axis = 2_i32
      box_max(axis) = spec%pos_low(axis)
    case ('z_low')
      axis = 3_i32
      box_min(axis) = spec%pos_low(axis)
    case ('z_high')
      axis = 3_i32
      box_max(axis) = spec%pos_low(axis)
    case default
      error stop 'plane_source has invalid derived face.'
    end select
  end subroutine configure_plane_source_box

  !> face bit順の面名を返す。
  pure subroutine boundary_face_name(face, name)
    integer(i32), intent(in) :: face
    character(len=*), intent(out) :: name

    select case (face)
    case (1_i32)
      name = 'x_low'
    case (2_i32)
      name = 'x_high'
    case (3_i32)
      name = 'y_low'
    case (4_i32)
      name = 'y_high'
    case (5_i32)
      name = 'z_low'
    case (6_i32)
      name = 'z_high'
    case default
      error stop 'invalid boundary inflow face.'
    end select
  end subroutine boundary_face_name

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
  !> reservoir_face 注入に対する法線速度補正パラメータを計算する。
  !! @param[in] cfg シミュレーション・結合設定を含むアプリ設定。
  !! @param[in] spec reservoir_face 粒子種設定。
  !! @param[out] vmin_normal 無限遠法線速度の下限 [m/s]。
  !! @param[out] barrier_normal 法線エネルギー障壁 `2 q Δφ / m` [`m^2/s^2`]。
  !! @param[in] mesh 現在バッチ開始時点の電荷分布メッシュ（補正時に必要）。
  !! @param[inout] snapshot refresh 済み静電 snapshot（infinity barrier 使用時に必要）。
  !! @param[in] warn_face_variation 面平均近似の電位ばらつき警告を出すか。
  subroutine reservoir_face_velocity_correction( &
    cfg, spec, vmin_normal, barrier_normal, mesh, snapshot, warn_face_variation, boundary_contract &
    )
    type(app_config), intent(in) :: cfg
    type(particle_species_spec), intent(in) :: spec
    real(dp), intent(out) :: vmin_normal
    real(dp), intent(out) :: barrier_normal
    type(mesh_type), intent(in), optional :: mesh
    type(electrostatic_snapshot_type), intent(inout), optional :: snapshot
    logical, intent(in), optional :: warn_face_variation
    type(external_boundary_contract_type), intent(in), optional :: boundary_contract

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
  end subroutine reservoir_face_velocity_correction

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
    characteristic_energy = k_boltzmann*species_temperature_k(spec) + &
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

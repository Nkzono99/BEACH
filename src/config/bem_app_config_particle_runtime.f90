!> app_config から粒子源計画と粒子バッチを構築する実行時変換。
module bem_app_config_particle_runtime
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
  use bem_types, only: mesh_type, particles_soa, sim_config, injection_state, bc_periodic
  use bem_mpi, only: mpi_context, mpi_get_rank_size, mpi_split_count, mpi_bcast_i32_array, mpi_bcast_real_dp_array
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_collision, only: collision_query_grid_stalled, collision_query_image_limit, &
                           collision_query_index_range, collision_query_invalid_segment, collision_query_ok
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
    app_config, particle_species_spec, particles_per_batch_from_config, &
    total_particles_from_config, particle_inflow_reservoir
  use bem_app_config_potential_runtime, only: &
    compute_face_average_potential, warn_face_average_potential_variation, resolve_face_sampling_geometry
  use bem_string_utils, only: lower_ascii
  use bem_config_helpers, only: resolve_inward_normal
  use, intrinsic :: iso_fortran_env, only: error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private :: finalize_particle_batch_collision_query

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

end module bem_app_config_particle_runtime

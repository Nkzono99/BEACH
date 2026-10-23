!> app_config から粒子源計画と粒子バッチを構築する実行時変換。
module bem_app_config_particle_runtime
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
  use bem_types, only: mesh_type, particles_soa, sim_config, injection_state
  use bem_mpi, only: mpi_context, mpi_get_rank_size, mpi_split_count, mpi_bcast_i32_array, mpi_bcast_real_dp_array
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_status_codes, only: collision_query_ok, particle_status_name
  use bem_injection, only: &
    seed_rng, sample_uniform_positions, sample_shifted_maxwell_velocities, compute_macro_particles_for_batch, &
    compute_macro_particles_from_flux, sample_reservoir_face_particles, &
    sample_photo_raycast_particles, &
    compute_inflow_flux_from_drifting_maxwellian, compute_face_area_from_bounds
  use bem_external_boundary_contract, only: &
    external_boundary_contract_type, external_boundary_ok, external_inflow_none, external_inflow_scalar_barrier, &
    resolve_external_boundary_contract
  use bem_app_config_types, only: &
    app_config, particle_species_spec, particles_per_batch_from_config, particle_inflow_reservoir
  use bem_app_config_potential_runtime, only: &
    compute_face_average_potential, warn_face_average_potential_variation, resolve_face_sampling_geometry
  use bem_string_utils, only: lower_ascii
  use bem_config_helpers, only: resolve_inward_normal
  use, intrinsic :: iso_fortran_env, only: error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

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
    logical, allocatable :: number_flux_override_active(:)
    real(dp), allocatable :: effective_temperature_k(:)
    real(dp), allocatable :: effective_drift_velocity(:, :)
    real(dp), allocatable :: effective_weight(:)
    real(dp), allocatable :: photo_emit_current_density(:)
    real(dp), allocatable :: photo_normal_drift_speed(:)
    logical, allocatable :: kinetic_inflow_active(:)
    real(dp), allocatable :: kinetic_reservoir_potential_v(:)
    real(dp), allocatable :: kinetic_access_potential_v(:)
    !> 流入写像の貯留部密度。0なら species 設定の密度を使う。
    real(dp), allocatable :: kinetic_reservoir_density_m3(:)
    integer(i32), allocatable :: kinetic_inflow_face(:)
  end type particle_source_plan_type

  !> slot の種別。
  integer(i32), parameter, public :: slot_kind_none = 0_i32
  integer(i32), parameter, public :: slot_kind_volume = 1_i32
  integer(i32), parameter, public :: slot_kind_photo_ray = 2_i32
  integer(i32), parameter, public :: slot_kind_reservoir_face = 3_i32

  !> 1 バッチの注入 layout。全 rank で同じ global 粒子数と面ごとの速度補正を持ち、
  !! `static` の rank 分割と `shared` の slot pool の両方がこれを基に粒子を作る。
  !! segment は species 順に source（volume_seed または photo ray）、次に面 1..6 の reservoir が並ぶ。
  type, public :: batch_injection_layout_type
    logical :: ready = .false.
    integer(i32) :: batch_idx = 0_i32
    integer(i32) :: nspecies = 0_i32
    integer(i32), allocatable :: source_kind(:)
    integer(i32), allocatable :: source_global_count(:)
    integer(i32), allocatable :: boundary_global_count(:, :)
    real(dp), allocatable :: boundary_vmin(:, :)
    real(dp), allocatable :: boundary_barrier(:, :)
    real(dp), allocatable :: weight(:)
    real(dp), allocatable :: temperature_k(:)
    real(dp), allocatable :: drift_velocity(:, :)
    real(dp), allocatable :: photo_emit_current_density(:)
    real(dp), allocatable :: photo_normal_drift_speed(:)
    !> 1 本の ray が当たったときのマクロ粒子重み（全 rank 合計の ray 数で分担）。
    real(dp), allocatable :: photo_ray_weight(:)
    integer(i32) :: nsegments = 0_i32
    integer(i32) :: total_slots = 0_i32
    integer(i32), allocatable :: segment_species(:)
    integer(i32), allocatable :: segment_kind(:)
    integer(i32), allocatable :: segment_face(:)
    integer(i32), allocatable :: segment_begin(:)
    integer(i32), allocatable :: segment_end(:)
    !> global slot -> (segment, segment 内 index)。slot は segment を round-robin に回って並ぶので、
    !! 取り合う順に種が混ざり、static の round-robin 組み立てと同じ処理順になる。
    integer(i32), allocatable :: slot_segment(:)
    integer(i32), allocatable :: slot_index(:)
  end type batch_injection_layout_type

  interface
    module subroutine build_batch_injection_layout( &
      cfg, batch_idx, layout, state, mesh, mpi_rank, mpi_size, mpi, snapshot, source_plan &
      )
      type(app_config), intent(in) :: cfg
      integer(i32), intent(in) :: batch_idx
      type(batch_injection_layout_type), intent(out) :: layout
      type(injection_state), intent(inout), optional :: state
      type(mesh_type), intent(in), optional :: mesh
      integer(i32), intent(in), optional :: mpi_rank, mpi_size
      type(mpi_context), intent(in), optional :: mpi
      type(electrostatic_snapshot_type), intent(inout), optional :: snapshot
      type(particle_source_plan_type), intent(in), optional, target :: source_plan
    end subroutine build_batch_injection_layout

    module subroutine generate_batch_slot( &
      cfg, layout, slot, mesh, x, v, w, q, m, species_idx, source_element, emitted, collision_status, bounce_count &
      )
      type(app_config), intent(in) :: cfg
      type(batch_injection_layout_type), intent(in) :: layout
      integer(i32), intent(in) :: slot
      type(mesh_type), intent(in), optional :: mesh
      real(dp), intent(out) :: x(3), v(3), w, q, m
      integer(i32), intent(out) :: species_idx, source_element, collision_status, bounce_count
      logical, intent(out) :: emitted
    end subroutine generate_batch_slot

    module subroutine init_particle_batch_from_config( &
      cfg, batch_idx, pcls, state, mesh, photo_emission_dq, mpi_rank, mpi_size, mpi, &
      collision_failure_status, collision_failure_species, collision_failure_ray, collision_failure_bounce, snapshot, &
      source_plan, photo_emission_dq_by_species, slot_ids &
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
      !> 各粒子の global slot 番号（粒子追跡の境界・二次電子乱数の鍵に使う）。
      integer(i32), allocatable, intent(out), optional :: slot_ids(:)
    end subroutine init_particle_batch_from_config

    module subroutine sample_species_state( &
      spec, n, x, v, temperature_k_override, drift_velocity_override &
      )
      type(particle_species_spec), intent(in) :: spec
      integer(i32), intent(in) :: n
      real(dp), intent(out) :: x(:, :)
      real(dp), intent(out) :: v(:, :)
      real(dp), intent(in), optional :: temperature_k_override
      real(dp), intent(in), optional :: drift_velocity_override(3)
    end subroutine sample_species_state

    module subroutine sample_photo_species_state( &
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
    end subroutine sample_photo_species_state

    module subroutine reservoir_face_velocity_correction( &
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
    end subroutine reservoir_face_velocity_correction

    module subroutine external_kinetic_face_velocity_correction( &
      cfg, spec, reservoir_potential_v, access_potential_v, vmin_normal, barrier_normal, &
      mesh, snapshot, warn_face_variation &
      )
      type(app_config), intent(in) :: cfg
      type(particle_species_spec), intent(in) :: spec
      real(dp), intent(in) :: reservoir_potential_v, access_potential_v
      real(dp), intent(out) :: vmin_normal, barrier_normal
      type(mesh_type), intent(in), optional :: mesh
      type(electrostatic_snapshot_type), intent(inout), optional :: snapshot
      logical, intent(in), optional :: warn_face_variation
    end subroutine external_kinetic_face_velocity_correction
  end interface

contains

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
  subroutine build_particle_source_plan( &
    cfg, plan, mpi_rank, mpi_size, mpi, kinetic_inflow_active, kinetic_reservoir_potential_v, &
    kinetic_access_potential_v, kinetic_inflow_face, number_flux_override_active, number_flux_override_m2_s, &
    kinetic_reservoir_density_m3 &
    )
    type(app_config), intent(in) :: cfg
    type(particle_source_plan_type), intent(out) :: plan
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    type(mpi_context), intent(in), optional :: mpi
    logical, intent(in), optional :: kinetic_inflow_active(:)
    real(dp), intent(in), optional :: kinetic_reservoir_potential_v(:), kinetic_access_potential_v(:)
    integer(i32), intent(in), optional :: kinetic_inflow_face(:)
    logical, intent(in), optional :: number_flux_override_active(:)
    real(dp), intent(in), optional :: number_flux_override_m2_s(:)
    !> 流入写像の面だけで使う貯留部密度。0以下の要素は species 設定の密度を使う。
    real(dp), intent(in), optional :: kinetic_reservoir_density_m3(:)

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
                              has_boundary_inflow(cfg%particle_species(s))
    end do
    plan%use_collective_reservoir_count = present(mpi) .and. has_enabled_reservoir

    allocate (plan%effective_density_m3(cfg%n_particle_species))
    allocate (plan%effective_particle_flux_m2_s(cfg%n_particle_species))
    allocate (plan%number_flux_override_active(cfg%n_particle_species))
    allocate (plan%effective_temperature_k(cfg%n_particle_species))
    allocate (plan%effective_drift_velocity(3, cfg%n_particle_species))
    allocate (plan%effective_weight(cfg%n_particle_species))
    allocate (plan%photo_emit_current_density(cfg%n_particle_species))
    allocate (plan%photo_normal_drift_speed(cfg%n_particle_species))
    allocate (plan%kinetic_inflow_active(cfg%n_particle_species))
    allocate (plan%kinetic_reservoir_potential_v(cfg%n_particle_species))
    allocate (plan%kinetic_access_potential_v(cfg%n_particle_species))
    allocate (plan%kinetic_reservoir_density_m3(cfg%n_particle_species))
    allocate (plan%kinetic_inflow_face(cfg%n_particle_species))
    plan%effective_density_m3 = 0.0_dp
    plan%effective_particle_flux_m2_s = 0.0_dp
    plan%number_flux_override_active = .false.
    plan%effective_temperature_k = 0.0_dp
    plan%effective_drift_velocity = 0.0_dp
    plan%effective_weight = 0.0_dp
    plan%photo_emit_current_density = 0.0_dp
    plan%photo_normal_drift_speed = 0.0_dp
    plan%kinetic_inflow_active = .false.
    plan%kinetic_reservoir_potential_v = 0.0_dp
    plan%kinetic_access_potential_v = 0.0_dp
    plan%kinetic_reservoir_density_m3 = 0.0_dp
    plan%kinetic_inflow_face = 0_i32

    if (present(kinetic_inflow_active) .or. present(kinetic_reservoir_potential_v) .or. &
        present(kinetic_access_potential_v) .or. present(kinetic_inflow_face)) then
      if (.not. present(kinetic_inflow_active) .or. .not. present(kinetic_reservoir_potential_v) .or. &
          .not. present(kinetic_access_potential_v) .or. .not. present(kinetic_inflow_face)) then
        error stop 'particle source kinetic map requires active, reservoir, access, and face arrays together.'
      end if
      if (size(kinetic_inflow_active) /= cfg%n_particle_species .or. &
          size(kinetic_reservoir_potential_v) /= cfg%n_particle_species .or. &
          size(kinetic_access_potential_v) /= cfg%n_particle_species .or. &
          size(kinetic_inflow_face) /= cfg%n_particle_species) then
        error stop 'particle source kinetic map species count mismatch.'
      end if
      if (.not. all(ieee_is_finite(kinetic_reservoir_potential_v)) .or. &
          .not. all(ieee_is_finite(kinetic_access_potential_v))) then
        error stop 'particle source kinetic map potentials must be finite.'
      end if
      plan%kinetic_inflow_active = kinetic_inflow_active
      plan%kinetic_reservoir_potential_v = kinetic_reservoir_potential_v
      plan%kinetic_access_potential_v = kinetic_access_potential_v
      plan%kinetic_inflow_face = kinetic_inflow_face
      if (any(plan%kinetic_inflow_active .and. &
              (plan%kinetic_inflow_face < 1_i32 .or. plan%kinetic_inflow_face > 6_i32))) then
        error stop 'active particle source kinetic map face must be in [1, 6].'
      end if
    end if

    if (present(kinetic_reservoir_density_m3)) then
      if (size(kinetic_reservoir_density_m3) /= cfg%n_particle_species) then
        error stop 'particle source kinetic reservoir density species count mismatch.'
      end if
      if (.not. all(ieee_is_finite(kinetic_reservoir_density_m3))) then
        error stop 'particle source kinetic reservoir densities must be finite.'
      end if
      plan%kinetic_reservoir_density_m3 = max(0.0_dp, kinetic_reservoir_density_m3)
    end if

    if (present(number_flux_override_active) .or. present(number_flux_override_m2_s)) then
      if (.not. present(number_flux_override_active) .or. .not. present(number_flux_override_m2_s)) then
        error stop 'particle source number-flux override requires active and flux arrays together.'
      end if
      if (size(number_flux_override_active) /= cfg%n_particle_species .or. &
          size(number_flux_override_m2_s) /= cfg%n_particle_species) then
        error stop 'particle source number-flux override species count mismatch.'
      end if
      if (any(number_flux_override_active .and. &
              (.not. ieee_is_finite(number_flux_override_m2_s) .or. number_flux_override_m2_s < 0.0_dp))) then
        error stop 'active particle source number-flux overrides must be finite and nonnegative.'
      end if
      plan%number_flux_override_active = number_flux_override_active
    end if

    do s = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(s)%enabled) cycle
      select case (trim(lower_ascii(cfg%particle_species(s)%source_mode)))
      case ('volume_seed')
        plan%effective_weight(s) = cfg%particle_species(s)%w_particle
        plan%effective_temperature_k(s) = species_temperature_k(cfg%particle_species(s))
        plan%effective_drift_velocity(:, s) = cfg%particle_species(s)%drift_velocity
        if (has_boundary_inflow(cfg%particle_species(s))) then
          plan%effective_density_m3(s) = species_number_density_m3(cfg%particle_species(s))
        end if
      case ('photo_raycast')
        plan%photo_emit_current_density(s) = cfg%particle_species(s)%emit_current_density_a_m2
        plan%photo_normal_drift_speed(s) = cfg%particle_species(s)%normal_drift_speed
      end select
      if (plan%number_flux_override_active(s)) then
        plan%effective_particle_flux_m2_s(s) = number_flux_override_m2_s(s)
      end if
    end do

    plan%ready = .true.
  end subroutine build_particle_source_plan

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

  !> box 全面の境界流入をサンプリングする面座標を構築する。
  pure subroutine make_boundary_inflow_spec(sim, source_spec, face, inflow_spec)
    type(sim_config), intent(in) :: sim
    type(particle_species_spec), intent(in) :: source_spec
    integer(i32), intent(in) :: face
    type(particle_species_spec), intent(out) :: inflow_spec
    integer(i32) :: axis
    logical :: high_side

    inflow_spec = source_spec
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
    inflow_spec%boundary_inflow_low = 0_i32
    inflow_spec%boundary_inflow_high = 0_i32
  end subroutine make_boundary_inflow_spec

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

  !> 境界流入用に、物理流量と残差から今バッチのマクロ粒子数を決める。
  !! @param[in] sim ボックス境界・バッチ時間などのシミュレーション設定。
  !! @param[in] spec 境界流入面の粒子種設定。
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
    temperature_k_override, drift_velocity_override, particle_flux_override, use_particle_flux_override &
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
    logical, intent(in), optional :: use_particle_flux_override

    real(dp) :: number_density_m3, effective_batch_duration, w_particle, temperature_k_local
    real(dp) :: drift_velocity_local(3)
    logical :: direct_flux

    number_density_m3 = species_number_density_m3(spec)
    if (present(number_density_override)) number_density_m3 = number_density_override
    w_particle = spec%w_particle
    if (present(w_particle_override)) w_particle = w_particle_override
    temperature_k_local = species_temperature_k(spec)
    if (present(temperature_k_override)) temperature_k_local = temperature_k_override
    drift_velocity_local = spec%drift_velocity
    if (present(drift_velocity_override)) drift_velocity_local = drift_velocity_override
    effective_batch_duration = sim%batch_duration
    direct_flux = .false.
    if (present(use_particle_flux_override) .and. present(particle_flux_override)) then
      direct_flux = use_particle_flux_override
    end if
    if (direct_flux) then
      call compute_macro_particles_from_flux( &
        particle_flux_override, spec%inject_face, spec%pos_low, spec%pos_high, effective_batch_duration, w_particle, residual, &
        count &
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

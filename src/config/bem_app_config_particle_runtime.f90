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
    integer(i32), allocatable :: kinetic_inflow_face(:)
  end type particle_source_plan_type

  interface
    module subroutine init_particle_batch_from_config( &
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
    end subroutine init_particle_batch_from_config

    module subroutine sample_species_state( &
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
    kinetic_access_potential_v, kinetic_inflow_face, number_flux_override_active, number_flux_override_m2_s &
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
    allocate (plan%number_flux_override_active(cfg%n_particle_species))
    allocate (plan%effective_temperature_k(cfg%n_particle_species))
    allocate (plan%effective_drift_velocity(3, cfg%n_particle_species))
    allocate (plan%effective_weight(cfg%n_particle_species))
    allocate (plan%photo_emit_current_density(cfg%n_particle_species))
    allocate (plan%photo_normal_drift_speed(cfg%n_particle_species))
    allocate (plan%kinetic_inflow_active(cfg%n_particle_species))
    allocate (plan%kinetic_reservoir_potential_v(cfg%n_particle_species))
    allocate (plan%kinetic_access_potential_v(cfg%n_particle_species))
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
      if (plan%number_flux_override_active(s)) then
        plan%effective_particle_flux_m2_s(s) = number_flux_override_m2_s(s)
      end if
    end do

    plan%ready = .true.
  end subroutine build_particle_source_plan

  !> reservoir/plane sourceの初期位置を、設定されたbox境界条件に従って有効領域へ正規化する。
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

  !> 面名をboundary bit順のindexへ変換する。
  pure integer(i32) function injection_face_index(face_name) result(face)
    character(len=*), intent(in) :: face_name

    select case (trim(lower_ascii(face_name)))
    case ('x_low')
      face = 1_i32
    case ('x_high')
      face = 2_i32
    case ('y_low')
      face = 3_i32
    case ('y_high')
      face = 4_i32
    case ('z_low')
      face = 5_i32
    case ('z_high')
      face = 6_i32
    case default
      face = 0_i32
    end select
  end function injection_face_index

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

    real(dp) :: number_density_m3, effective_batch_duration, particle_flux_m2_s, w_particle, temperature_k_local
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
    direct_flux = trim(lower_ascii(spec%velocity_distribution)) == 'grid'
    if (present(use_particle_flux_override)) direct_flux = direct_flux .or. use_particle_flux_override
    if (direct_flux) then
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

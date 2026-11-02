!> app_config から粒子源計画と粒子バッチを構築する実行時変換。
module bem_app_config_particle_runtime
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
  use bem_types, only: mesh_type, particles_soa, sim_config, injection_state
  use bem_mpi, only: mpi_context, mpi_get_rank_size, mpi_split_count, mpi_bcast_i32_array, mpi_bcast_real_dp_array
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_status_codes, only: collision_query_ok, particle_status_name, photo_ray_failure_type
  use bem_injection, only: seed_rng, compute_macro_particles_for_batch
  use bem_external_boundary_contract, only: &
    external_inflow_none, external_inflow_scalar_barrier
  use bem_app_config_types, only: &
    app_config, particle_species_spec, particles_per_batch_from_config, particle_inflow_reservoir, species_config_path
  use bem_surface_closure_contract, only: surface_closure_contract_type
  use bem_app_config_potential_runtime, only: &
    compute_face_average_potential, warn_face_average_potential_variation, resolve_face_sampling_geometry
  use bem_string_utils, only: lower_ascii
  use bem_config_helpers, only: resolve_inward_normal, species_number_density_m3, species_temperature_k, &
                                boundary_face_name, boundary_inflow_face_enabled
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
    integer(i32), allocatable :: source_kind(:)
    !> 外部 closure の写し。reservoir 面の流入写像（species 別の有効フラグ・経路電位・貯留部密度と共通の面）だけを読む。
    type(surface_closure_contract_type) :: closure
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
    integer(i32), allocatable :: source_kind(:)
    integer(i32), allocatable :: source_global_count(:)
    integer(i32), allocatable :: boundary_global_count(:, :)
    real(dp), allocatable :: boundary_vmin(:, :)
    real(dp), allocatable :: boundary_barrier(:, :)
    !> 1 本の ray が当たったときのマクロ粒子重み（全 rank 合計の ray 数で分担）。
    real(dp), allocatable :: photo_ray_weight(:)
    integer(i32) :: nsegments = 0_i32
    integer(i32) :: total_slots = 0_i32
    integer(i32), allocatable :: segment_species(:)
    integer(i32), allocatable :: segment_kind(:)
    integer(i32), allocatable :: segment_face(:)
    integer(i32), allocatable :: segment_begin(:)
    integer(i32), allocatable :: segment_end(:)
    !> slot は segment を round-robin に回って並ぶ（segment 内 index k の周回で、長さが k 以上の segment を順に取る）。
    !! 取り合う順に種が混ざり、static の round-robin 組み立てと同じ処理順になる。slot ごとの表は持たず、
    !! segment 長の相異なる値ごとの帯で `locate_slot` が (segment, segment 内 index) を求める。
    !! 帯 b は segment 内 index が `band_length(b-1)+1 .. band_length(b)` の周回で、`band_length(0)=0`。
    integer(i32) :: nbands = 0_i32
    integer(i32), allocatable :: band_length(:)
    !> 帯 b の最後の global slot（`band_last_slot(0)=0`）。
    integer(i32), allocatable :: band_last_slot(:)
  contains
    procedure :: locate_slot => locate_batch_slot
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
      type(particle_source_plan_type), intent(in) :: source_plan
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
      cfg, batch_idx, pcls, state, mesh, photo_emission_dq, mpi_rank, mpi_size, mpi, photo_failure, snapshot, &
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
      !> 不完全な photo collision query（省略時は停止）。
      type(photo_ray_failure_type), intent(out), optional :: photo_failure
      type(electrostatic_snapshot_type), intent(inout), optional :: snapshot
      type(particle_source_plan_type), intent(in) :: source_plan
      real(dp), intent(out), optional :: photo_emission_dq_by_species(:, :)
      !> 各粒子の global slot 番号（粒子追跡の境界・二次電子乱数の鍵に使う）。
      integer(i32), allocatable, intent(out), optional :: slot_ids(:)
    end subroutine init_particle_batch_from_config

    module subroutine reservoir_face_velocity_correction( &
      cfg, spec, vmin_normal, barrier_normal, mesh, snapshot, warn_face_variation, face_potential, face_potential_ready &
      )
      type(app_config), intent(in) :: cfg
      type(particle_species_spec), intent(in) :: spec
      real(dp), intent(out) :: vmin_normal
      real(dp), intent(out) :: barrier_normal
      type(mesh_type), intent(in), optional :: mesh
      type(electrostatic_snapshot_type), intent(inout), optional :: snapshot
      logical, intent(in), optional :: warn_face_variation
      !> 面の電位統計 (平均, 標準偏差, 最小, 最大) [V] の cache。`face_potential_ready` が真なら評価せずに使い、
      !! 偽なら評価して書き込む。box 面全体の統計なので、同じ snapshot の 1 layout 内では種によらず共有できる。
      real(dp), intent(inout), optional :: face_potential(4)
      logical, intent(inout), optional :: face_potential_ready
    end subroutine reservoir_face_velocity_correction

    module subroutine external_kinetic_face_velocity_correction( &
      cfg, spec, access_potential_v, vmin_normal, barrier_normal, &
      mesh, snapshot, warn_face_variation, face_potential, face_potential_ready &
      )
      type(app_config), intent(in) :: cfg
      type(particle_species_spec), intent(in) :: spec
      real(dp), intent(in) :: access_potential_v
      real(dp), intent(out) :: vmin_normal, barrier_normal
      type(mesh_type), intent(in), optional :: mesh
      type(electrostatic_snapshot_type), intent(inout), optional :: snapshot
      logical, intent(in), optional :: warn_face_variation
      !> 面の電位統計 (平均, 標準偏差, 最小, 最大) [V] の cache。`face_potential_ready` が真なら評価せずに使い、
      !! 偽なら評価して書き込む。box 面全体の統計なので、同じ snapshot の 1 layout 内では種によらず共有できる。
      real(dp), intent(inout), optional :: face_potential(4)
      logical, intent(inout), optional :: face_potential_ready
    end subroutine external_kinetic_face_velocity_correction
  end interface

contains

  !> global slot から (segment, segment 内 index) を求める。O(nsegments) で、slot ごとの表を持たない。
  !! @param[in] slot `1..total_slots` の global slot。
  !! @param[out] segment slot が属する segment。
  !! @param[out] local_index segment 内の 1 始まり index。
  pure subroutine locate_batch_slot(self, slot, segment, local_index)
    class(batch_injection_layout_type), intent(in) :: self
    integer(i32), intent(in) :: slot
    integer(i32), intent(out) :: segment, local_index
    integer(i32) :: band, active, offset, rank_in_round, g

    band = 1_i32
    do while (self%band_last_slot(band) < slot)
      band = band + 1_i32
    end do
    active = count(self%segment_end - self%segment_begin + 1_i32 >= self%band_length(band))
    offset = slot - self%band_last_slot(band - 1_i32) - 1_i32
    local_index = self%band_length(band - 1_i32) + offset/active + 1_i32
    rank_in_round = mod(offset, active) + 1_i32
    segment = 0_i32
    do g = 1_i32, self%nsegments
      if (self%segment_end(g) - self%segment_begin(g) + 1_i32 < self%band_length(band)) cycle
      rank_in_round = rank_in_round - 1_i32
      if (rank_in_round == 0_i32) then
        segment = g
        return
      end if
    end do
  end subroutine locate_batch_slot

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
  !! closure を渡すと、その流入写像を reservoir 面の注入に使う。省略時は流入写像なし。
  subroutine build_particle_source_plan(cfg, plan, mpi_rank, mpi_size, mpi, closure)
    type(app_config), intent(in) :: cfg
    type(particle_source_plan_type), intent(out) :: plan
    integer(i32), intent(in), optional :: mpi_rank, mpi_size
    type(mpi_context), intent(in), optional :: mpi
    !> 外部 closure。species 配列の大きさは cfg の粒子種数。
    type(surface_closure_contract_type), intent(in), optional :: closure

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

    allocate (plan%source_kind(cfg%n_particle_species))
    plan%source_kind = slot_kind_none
    if (present(closure)) then
      ! 貯留部密度は外部根から来る。経路電位は使う所（external_kinetic_face_velocity_correction）が確かめる。
      if (.not. all(ieee_is_finite(closure%species%inflow_reservoir_density_m3))) then
        error stop 'particle source kinetic reservoir densities must be finite.'
      end if
      plan%closure = closure
    else
      allocate (plan%closure%species(cfg%n_particle_species))
    end if

    do s = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(s)%enabled) cycle
      select case (trim(lower_ascii(cfg%particle_species(s)%source_mode)))
      case ('volume_seed')
        plan%source_kind(s) = slot_kind_volume
      case ('photo_raycast')
        plan%source_kind(s) = slot_kind_photo_ray
      case ('secondary_emission')
        continue
      case default
        error stop 'Unknown '//species_config_path(int(s))//'.source.mode.'
      end select
    end do

    plan%ready = .true.
  end subroutine build_particle_source_plan

  !> source plan がこの rank/size と mpi 引数の有無で構築されたかを返す。
  !! rank、size、mpi 引数の有無はそれぞれ独立に一致する必要がある。
  !! @param[in] plan `build_particle_source_plan` で構築した source plan。
  !! @param[in] local_rank 今回解決した rank 番号。
  !! @param[in] n_ranks 今回解決した rank 数。
  !! @param[in] mpi_argument_present 今回 `mpi_context` を渡したか。
  !! @return matches 3 つがすべて plan 構築時と一致すれば `.true.`。
  pure logical function particle_source_plan_matches_mpi( &
    plan, local_rank, n_ranks, mpi_argument_present &
    ) result(matches)
    type(particle_source_plan_type), intent(in) :: plan
    integer(i32), intent(in) :: local_rank, n_ranks
    logical, intent(in) :: mpi_argument_present

    matches = plan%mpi_rank == local_rank .and. plan%mpi_size == n_ranks .and. &
              (plan%mpi_argument_present .eqv. mpi_argument_present)
  end function particle_source_plan_matches_mpi

  !> speciesに外部 reservoir 流入を指定したbox面があるかを返す。
  pure logical function has_boundary_inflow(spec) result(has_inflow)
    type(particle_species_spec), intent(in) :: spec

    has_inflow = any(spec%boundary_inflow_low == particle_inflow_reservoir) .or. &
                 any(spec%boundary_inflow_high == particle_inflow_reservoir)
  end function has_boundary_inflow

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

  !> 境界流入用に、物理流量と残差から今バッチのマクロ粒子数を決める。
  !! @param[in] sim ボックス境界・バッチ時間などのシミュレーション設定。
  !! @param[in] spec 境界流入面の粒子種設定。
  !! @param[inout] residual 前バッチから繰り越した端数。
  !! @param[out] count 今バッチで生成するマクロ粒子数。
  !! @param[in] vmin_normal 法線速度の下限 [m/s]（省略時は 0）。
  !! @param[in] number_density_override 数密度の上書き値 [1/m^3]。
  subroutine compute_macro_particles_for_species( &
    sim, spec, residual, count, vmin_normal, number_density_override &
    )
    type(sim_config), intent(in) :: sim
    type(particle_species_spec), intent(in) :: spec
    real(dp), intent(inout) :: residual
    integer(i32), intent(out) :: count
    real(dp), intent(in), optional :: vmin_normal
    real(dp), intent(in), optional :: number_density_override

    real(dp) :: number_density_m3, effective_batch_duration, w_particle, temperature_k_local
    real(dp) :: drift_velocity_local(3)

    number_density_m3 = species_number_density_m3(spec)
    if (present(number_density_override)) number_density_m3 = number_density_override
    w_particle = spec%w_particle
    temperature_k_local = species_temperature_k(spec)
    drift_velocity_local = spec%drift_velocity
    effective_batch_duration = sim%batch_duration
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

end module bem_app_config_particle_runtime

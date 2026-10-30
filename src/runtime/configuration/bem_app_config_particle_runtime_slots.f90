!> バッチ注入 layout の構築と、layout 上の slot 1 つを counter-based 乱数で生成する処理。
!> layout は rank に依らない global 粒子数と面ごとの速度補正を持ち、`static` の rank 分割と
!> `shared` の slot pool の両方がこれを基に粒子を作る。
submodule(bem_app_config_particle_runtime) bem_app_config_particle_runtime_slots
  use bem_kinds, only: i64
  use bem_injection, only: sample_reservoir_face_particle_counter, sample_photo_ray_slot, photo_ray_macro_weight
  use bem_injection_random, only: sample_standard_normal
  use bem_counter_rng, only: counter_uniforms
  implicit none
contains

  !> 指定バッチの注入 layout を構築する。
  !! 粒子数は全 rank で同じ global 値を持ち、reservoir の残差更新は従来どおり rank 0 で行って broadcast する。
  !! @param[in] cfg 粒子種とシミュレーション条件を含むアプリ設定。
  !! @param[in] batch_idx 生成対象のバッチ番号（1 始まり）。
  !! @param[out] layout 構築した layout。
  !! @param[inout] state 境界 reservoir 流入の残差状態（reservoir を使う種があれば必須）。
  !! @param[in] mesh 現在バッチ開始時点のメッシュ（電位補正時に必要）。
  !! @param[inout] snapshot refresh 済み静電 snapshot（注入電位補正の使用時に必要）。
  !! @param[in] source_plan run 中に再利用する粒子 source 導出値（`build_particle_source_plan` で同じ MPI 文脈から構築）。
  module procedure build_batch_injection_layout
  integer(i32) :: s, face, local_rank, n_ranks, segment
  logical :: use_collective_reservoir_count
  type(particle_species_spec) :: face_spec
  real(dp) :: correction_vmin_normal, kinetic_vmin_normal, kinetic_barrier_normal, face_density_m3
  ! A boundary face spans the whole box side, so its potential statistics are shared by every species.
  real(dp) :: face_potential(4, 6)
  logical :: face_potential_ready(6)

  layout = batch_injection_layout_type()
  if (batch_idx < 1_i32 .or. batch_idx > cfg%sim%batch_count) then
    error stop 'Requested batch index is out of range.'
  end if
  call resolve_parallel_rank_size(local_rank, n_ranks, mpi_rank, mpi_size, mpi, 'build_batch_injection_layout')
  if (.not. source_plan%ready) error stop 'particle source plan is not initialized.'
  if (source_plan%nspecies /= cfg%n_particle_species) then
    error stop 'particle source plan species count does not match app config.'
  end if
  if (.not. particle_source_plan_matches_mpi(source_plan, local_rank, n_ranks, present(mpi))) then
    error stop 'particle source plan MPI context does not match batch initialization.'
  end if
  use_collective_reservoir_count = source_plan%use_collective_reservoir_count
  if (present(state)) then
    if (.not. allocated(state%macro_residual)) error stop 'injection_state is not initialized.'
    if (size(state%macro_residual) < cfg%n_particle_species) error stop 'injection_state size mismatch.'
    if (any([(has_boundary_inflow(cfg%particle_species(s)), s=1, cfg%n_particle_species)])) then
      if (.not. allocated(state%boundary_macro_residual)) then
        error stop 'inflow requires boundary_macro_residual in injection_state.'
      end if
      if (size(state%boundary_macro_residual, 1) < 6 .or. &
          size(state%boundary_macro_residual, 2) < cfg%n_particle_species) then
        error stop 'injection_state boundary residual size mismatch.'
      end if
    end if
  end if

  layout%batch_idx = batch_idx
  layout%nspecies = cfg%n_particle_species
  allocate (layout%source_kind(cfg%n_particle_species), layout%source_global_count(cfg%n_particle_species))
  allocate (layout%boundary_global_count(6, cfg%n_particle_species))
  allocate (layout%boundary_vmin(6, cfg%n_particle_species), layout%boundary_barrier(6, cfg%n_particle_species))
  allocate (layout%weight(cfg%n_particle_species), layout%temperature_k(cfg%n_particle_species))
  allocate (layout%drift_velocity(3, cfg%n_particle_species))
  allocate (layout%photo_emit_current_density(cfg%n_particle_species))
  allocate (layout%photo_normal_drift_speed(cfg%n_particle_species), layout%photo_ray_weight(cfg%n_particle_species))
  layout%source_kind = source_plan%source_kind
  layout%source_global_count = 0_i32
  layout%boundary_global_count = 0_i32
  layout%boundary_vmin = 0.0_dp
  layout%boundary_barrier = 0.0_dp
  layout%weight = source_plan%effective_weight
  layout%temperature_k = source_plan%effective_temperature_k
  layout%drift_velocity = source_plan%effective_drift_velocity
  layout%photo_emit_current_density = source_plan%photo_emit_current_density
  layout%photo_normal_drift_speed = source_plan%photo_normal_drift_speed
  layout%photo_ray_weight = 0.0_dp
  face_potential_ready = .false.

  do s = 1, cfg%n_particle_species
    if (.not. cfg%particle_species(s)%enabled) cycle
    select case (layout%source_kind(s))
    case (slot_kind_volume)
      layout%source_global_count(s) = cfg%particle_species(s)%npcls_per_step
    case (slot_kind_photo_ray)
      ! A non-positive emission current emits nothing, so its rays are not even cast.
      if (layout%photo_emit_current_density(s) > 0.0_dp) then
        layout%source_global_count(s) = cfg%particle_species(s)%rays_per_batch
        layout%photo_ray_weight(s) = photo_ray_macro_weight( &
                                     cfg%particle_species(s)%inject_face, cfg%particle_species(s)%pos_low, &
                                     cfg%particle_species(s)%pos_high, cfg%particle_species(s)%ray_direction, &
                                     layout%photo_emit_current_density(s), cfg%particle_species(s)%q_particle, &
                                     cfg%sim%batch_duration, cfg%particle_species(s)%rays_per_batch &
                                     )
      end if
    case (slot_kind_none)
      continue
    end select
    if (has_boundary_inflow(cfg%particle_species(s))) then
      if (.not. present(state)) then
        error stop 'inflow requires injection_state in init_particle_batch_from_config.'
      end if
      do face = 1, 6
        if (.not. boundary_inflow_face_enabled(cfg%particle_species(s), face)) cycle
        call make_boundary_inflow_spec(cfg%sim, cfg%particle_species(s), face, face_spec)
        call reservoir_face_velocity_correction( &
          cfg, face_spec, correction_vmin_normal, layout%boundary_barrier(face, s), mesh, snapshot, &
          warn_face_variation=local_rank == 0_i32 .and. &
          (batch_idx == 1_i32 .or. batch_idx == cfg%sim%batch_count), &
          face_potential=face_potential(:, face), face_potential_ready=face_potential_ready(face) &
          )
        layout%boundary_vmin(face, s) = correction_vmin_normal
        face_density_m3 = source_plan%effective_density_m3(s)
        if (source_plan%kinetic_inflow_active(s) .and. face == source_plan%kinetic_inflow_face(s)) then
          call external_kinetic_face_velocity_correction( &
            cfg, face_spec, source_plan%kinetic_reservoir_potential_v(s), &
            source_plan%kinetic_access_potential_v(s), kinetic_vmin_normal, kinetic_barrier_normal, &
            mesh, snapshot, warn_face_variation=local_rank == 0_i32 .and. &
            (batch_idx == 1_i32 .or. batch_idx == cfg%sim%batch_count), &
            face_potential=face_potential(:, face), face_potential_ready=face_potential_ready(face) &
            )
          layout%boundary_vmin(face, s) = max(layout%boundary_vmin(face, s), kinetic_vmin_normal)
          layout%boundary_barrier(face, s) = kinetic_barrier_normal
          if (source_plan%kinetic_reservoir_density_m3(s) > 0.0_dp) then
            face_density_m3 = source_plan%kinetic_reservoir_density_m3(s)
          end if
        end if
        if (.not. use_collective_reservoir_count .or. local_rank == 0_i32) then
          call compute_macro_particles_for_species( &
            cfg%sim, face_spec, state%boundary_macro_residual(face, s), layout%boundary_global_count(face, s), &
            vmin_normal=layout%boundary_vmin(face, s), number_density_override=face_density_m3, &
            particle_flux_override=source_plan%effective_particle_flux_m2_s(s), &
            use_particle_flux_override=source_plan%number_flux_override_active(s), &
            temperature_k_override=layout%temperature_k(s), &
            drift_velocity_override=layout%drift_velocity(:, s), w_particle_override=layout%weight(s) &
            )
        end if
      end do
    end if
  end do
  if (use_collective_reservoir_count) then
    do s = 1, cfg%n_particle_species
      if (.not. cfg%particle_species(s)%enabled) cycle
      if (.not. has_boundary_inflow(cfg%particle_species(s))) cycle
      call mpi_bcast_i32_array(mpi, layout%boundary_global_count(:, s), 0_i32)
      call mpi_bcast_real_dp_array(mpi, state%boundary_macro_residual(:, s), 0_i32)
    end do
  end if

  ! Slot segments follow the static sampling order: species by species, the source first, then faces 1..6.
  layout%nsegments = count(layout%source_global_count > 0_i32) + count(layout%boundary_global_count > 0_i32)
  allocate (layout%segment_species(layout%nsegments), layout%segment_kind(layout%nsegments), &
            layout%segment_face(layout%nsegments), layout%segment_begin(layout%nsegments), &
            layout%segment_end(layout%nsegments))
  segment = 0_i32
  layout%total_slots = 0_i32
  do s = 1, cfg%n_particle_species
    if (layout%source_global_count(s) > 0_i32) then
      call push_segment(s, layout%source_kind(s), 0_i32, layout%source_global_count(s))
    end if
    do face = 1, 6
      if (layout%boundary_global_count(face, s) <= 0_i32) cycle
      call push_segment(s, slot_kind_reservoir_face, face, layout%boundary_global_count(face, s))
    end do
  end do
  ! Slots go out round-robin over the segments, as the static batch interleaves species; locate_slot
  ! resolves that order per band of equal active-segment count instead of storing it per slot.
  call build_slot_bands()
  layout%ready = .true.

contains
  !> 相異なる segment 長を昇順に並べ、帯ごとの最後の slot を求める。
  subroutine build_slot_bands()
    integer(i32), allocatable :: lengths(:)
    integer(i32) :: previous, current, g
    integer(i64) :: last_slot

    lengths = layout%segment_end - layout%segment_begin + 1_i32
    allocate (layout%band_length(0:layout%nsegments), layout%band_last_slot(0:layout%nsegments))
    layout%band_length(0) = 0_i32
    layout%band_last_slot(0) = 0_i32
    layout%nbands = 0_i32
    last_slot = 0_i64
    previous = 0_i32
    do
      current = huge(0_i32)
      do g = 1_i32, layout%nsegments
        if (lengths(g) > previous) current = min(current, lengths(g))
      end do
      if (current == huge(0_i32)) exit
      layout%nbands = layout%nbands + 1_i32
      last_slot = last_slot + int(current - previous, i64)*int(count(lengths >= current), i64)
      layout%band_length(layout%nbands) = current
      layout%band_last_slot(layout%nbands) = int(last_slot, i32)
      previous = current
    end do
    if (last_slot /= int(layout%total_slots, i64)) error stop 'batch injection slot order does not cover every slot.'
  end subroutine build_slot_bands

  subroutine push_segment(species, kind, face_index, count)
    integer(i32), intent(in) :: species, kind, face_index, count

    if (layout%total_slots > huge(0_i32) - count) error stop 'batch injection slot count overflow.'
    segment = segment + 1_i32
    layout%segment_species(segment) = species
    layout%segment_kind(segment) = kind
    layout%segment_face(segment) = face_index
    layout%segment_begin(segment) = layout%total_slots + 1_i32
    layout%segment_end(segment) = layout%total_slots + count
    layout%total_slots = layout%total_slots + count
  end subroutine push_segment
  end procedure build_batch_injection_layout

  !> layout 上の global slot 1 つを生成する。
  !! 乱数は (batch, species, kind, face, segment 内 index) の counter で決まり、どの rank・thread が生成しても同じ粒子になる。
  !! @param[in] slot `1..layout%total_slots` の global slot。
  !! @param[in] mesh photo_raycast の交差判定に使うメッシュ（photo 種があれば必須）。
  !! @param[out] emitted 粒子ができたとき `.true.`（photo ray が外れたときだけ `.false.`）。
  !! @param[out] collision_status photo ray の衝突照会 status（`collision_query_ok` 以外なら未完了）。
  !! @param[out] bounce_count 未完了の照会を返した bounce index。
  module procedure generate_batch_slot
  integer(i32) :: segment, s, kind, face, local_index
  integer(i64) :: counter(5)
  real(dp) :: sigma, z(3, 1), position_uniform(3)
  type(particle_species_spec) :: face_spec

  if (.not. layout%ready) error stop 'batch injection layout is not initialized.'
  if (slot < 1_i32 .or. slot > layout%total_slots) error stop 'batch injection slot is out of range.'
  call layout%locate_slot(slot, segment, local_index)
  s = layout%segment_species(segment)
  kind = layout%segment_kind(segment)
  face = layout%segment_face(segment)
  counter = int([layout%batch_idx, s, kind, face, local_index], i64)

  x = 0.0_dp
  v = 0.0_dp
  w = 0.0_dp
  species_idx = s
  q = cfg%particle_species(s)%q_particle
  m = cfg%particle_species(s)%m_particle
  source_element = -1_i32
  emitted = .false.
  collision_status = collision_query_ok
  bounce_count = 0_i32
  associate (spec => cfg%particle_species(s))
    select case (kind)
    case (slot_kind_volume)
      call counter_uniforms(cfg%sim%rng_seed, counter, 1_i64, position_uniform)
      x = spec%pos_low + (spec%pos_high - spec%pos_low)*position_uniform
      ! 質量 > 0 と温度 >= 0 は native preflight の validate_volume_seed_species が保証する。
      sigma = sqrt(k_boltzmann*layout%temperature_k(s)/spec%m_particle)
      call sample_standard_normal(z, rng_seed=cfg%sim%rng_seed, rng_counter=counter)
      v = sigma*z(:, 1) + layout%drift_velocity(:, s)
      w = layout%weight(s)
      emitted = .true.
    case (slot_kind_photo_ray)
      if (.not. present(mesh)) error stop 'photo_raycast slots require mesh in generate_batch_slot.'
      call sample_photo_ray_slot( &
        mesh, cfg%sim, spec%inject_face, spec%pos_low, spec%pos_high, spec%ray_direction, spec%m_particle, &
        species_temperature_k(spec), layout%photo_normal_drift_speed(s), layout%photo_ray_weight(s), &
        cfg%sim%rng_seed, [counter, 1_i64], [counter, 2_i64], x, v, w, source_element, emitted, &
        collision_status, bounce_count &
        )
    case (slot_kind_reservoir_face)
      call make_boundary_inflow_spec(cfg%sim, spec, face, face_spec)
      call sample_reservoir_face_particle_counter( &
        cfg%sim%box_min, cfg%sim%box_max, face_spec%inject_face, face_spec%pos_low, face_spec%pos_high, &
        layout%drift_velocity(:, s), spec%m_particle, layout%temperature_k(s), &
        layout%boundary_barrier(face, s), layout%boundary_vmin(face, s), cfg%sim%rng_seed, counter, x, v &
        )
      w = layout%weight(s)
      emitted = .true.
    case default
      error stop 'batch injection layout has an unknown slot kind.'
    end select
  end associate
  end procedure generate_batch_slot

end submodule bem_app_config_particle_runtime_slots

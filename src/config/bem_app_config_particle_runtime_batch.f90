!> 粒子源計画に従うMPI粒子数配分、サンプリング、バッチ組み立て。
submodule(bem_app_config_particle_runtime) bem_app_config_particle_runtime_batch
  implicit none
contains

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
  module procedure init_particle_batch_from_config

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
  real(dp) :: kinetic_vmin_normal, kinetic_barrier_normal
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
        if (active_source_plan%kinetic_inflow_active(s) .and. &
            injection_face_index(cfg%particle_species(s)%inject_face) == active_source_plan%kinetic_inflow_face(s)) then
          call external_kinetic_face_velocity_correction( &
            cfg, cfg%particle_species(s), active_source_plan%kinetic_reservoir_potential_v(s), &
            active_source_plan%kinetic_access_potential_v(s), kinetic_vmin_normal, kinetic_barrier_normal, &
            mesh, snapshot, warn_face_variation=local_rank == 0_i32 .and. &
            (batch_idx == 1_i32 .or. batch_idx == cfg%sim%batch_count) &
            )
          vmin_normal(s) = max(vmin_normal(s), kinetic_vmin_normal)
          barrier_normal(s) = kinetic_barrier_normal
        end if
      end if
      if (.not. use_collective_reservoir_count .or. local_rank == 0_i32) then
        call compute_macro_particles_for_species( &
          cfg%sim, cfg%particle_species(s), state%macro_residual(s), global_counts(s), vmin_normal=vmin_normal(s), &
          number_density_override=batch_density_m3(s), particle_flux_override=effective_particle_flux_m2_s(s), &
          use_particle_flux_override=active_source_plan%number_flux_override_active(s), &
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
        if (active_source_plan%kinetic_inflow_active(s) .and. face == active_source_plan%kinetic_inflow_face(s)) then
          call external_kinetic_face_velocity_correction( &
            cfg, face_spec, active_source_plan%kinetic_reservoir_potential_v(s), &
            active_source_plan%kinetic_access_potential_v(s), kinetic_vmin_normal, kinetic_barrier_normal, &
            mesh, snapshot, warn_face_variation=local_rank == 0_i32 .and. &
            (batch_idx == 1_i32 .or. batch_idx == cfg%sim%batch_count) &
            )
          boundary_vmin(face, s) = max(boundary_vmin(face, s), kinetic_vmin_normal)
          boundary_barrier(face, s) = kinetic_barrier_normal
        end if
        if (.not. use_collective_reservoir_count .or. local_rank == 0_i32) then
          call compute_macro_particles_for_species( &
            cfg%sim, face_spec, state%boundary_macro_residual(face, s), boundary_global_counts(face, s), &
            vmin_normal=boundary_vmin(face, s), number_density_override=batch_density_m3(s), &
            particle_flux_override=effective_particle_flux_m2_s(s), &
            use_particle_flux_override=active_source_plan%number_flux_override_active(s), &
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
  end procedure init_particle_batch_from_config

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

end submodule bem_app_config_particle_runtime_batch

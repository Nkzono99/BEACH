!> 粒子源計画に従うMPI粒子数配分、サンプリング、バッチ組み立て。
submodule(bem_app_config_particle_runtime) bem_app_config_particle_runtime_batch
  use bem_particles, only: allocate_particles
  use bem_app_config_types, only: species_config_path
  implicit none

  ! Each species owns only its requested capacity; ray hits can fill fewer entries.
  type :: sampled_species_type
    real(dp), allocatable :: x(:, :), v(:, :), w(:)
    integer(i32), allocatable :: source_element(:)
  end type sampled_species_type
contains

  !> 指定バッチ番号に対応する粒子バッチを生成する。
  !! @param[in] cfg 粒子種とシミュレーション条件を含むアプリ設定。
  !! @param[in] batch_idx 生成対象のバッチ番号（1始まり）。
  !! @param[out] pcls 生成したバッチ粒子群。
  !! @param[inout] state 境界 reservoir 流入の残差状態（必要時のみ）。
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

  integer(i32) :: s, i, face, local_rank, n_ranks
  integer(i32) :: source_begin, source_end, face_begin, face_end
  integer(i32) :: photo_collision_status, photo_collision_ray, photo_collision_bounce
  integer(i32), allocatable :: counts_max(:), counts_actual(:), source_counts(:), boundary_counts(:, :)
  type(sampled_species_type), allocatable :: sampled(:)
  type(batch_injection_layout_type) :: layout
  type(particle_species_spec) :: face_spec

  if (present(collision_failure_status)) collision_failure_status = collision_query_ok
  if (present(collision_failure_species)) collision_failure_species = huge(0_i32)
  if (present(collision_failure_ray)) collision_failure_ray = huge(0_i32)
  if (present(collision_failure_bounce)) collision_failure_bounce = 0_i32
  call resolve_parallel_rank_size(local_rank, n_ranks, mpi_rank, mpi_size, mpi, 'init_particle_batch_from_config')
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

  ! The layout holds the rank-independent counts and face corrections; this rank samples its static share.
  call build_batch_injection_layout( &
    cfg, batch_idx, layout, state=state, mesh=mesh, mpi_rank=mpi_rank, mpi_size=mpi_size, mpi=mpi, &
    snapshot=snapshot, source_plan=source_plan &
    )
  allocate ( &
    counts_max(cfg%n_particle_species), counts_actual(cfg%n_particle_species), source_counts(cfg%n_particle_species) &
    )
  allocate (boundary_counts(6, cfg%n_particle_species))
  counts_actual = 0_i32
  do s = 1, cfg%n_particle_species
    source_counts(s) = mpi_split_count(layout%source_global_count(s), local_rank, n_ranks)
    do face = 1, 6
      boundary_counts(face, s) = mpi_split_count(layout%boundary_global_count(face, s), local_rank, n_ranks)
    end do
  end do
  counts_max = source_counts + sum(boundary_counts, dim=1)
  associate ( &
    effective_temperature_k => layout%temperature_k, &
    effective_drift_velocity => layout%drift_velocity, &
    photo_emit_current_density => layout%photo_emit_current_density, &
    photo_normal_drift_speed => layout%photo_normal_drift_speed, &
    batch_weight => layout%weight, boundary_barrier => layout%boundary_barrier, boundary_vmin => layout%boundary_vmin &
    )
    allocate (sampled(cfg%n_particle_species))
    do s = 1, cfg%n_particle_species
      if (counts_max(s) <= 0_i32) cycle
      allocate (sampled(s)%x(3, counts_max(s)), sampled(s)%v(3, counts_max(s)), &
                sampled(s)%w(counts_max(s)), sampled(s)%source_element(counts_max(s)))
      sampled(s)%source_element = -1_i32
      associate (x_species => sampled(s)%x, v_species => sampled(s)%v, &
                 w_species => sampled(s)%w, emit_elem_species => sampled(s)%source_element)
        source_begin = 1_i32
        source_end = source_counts(s)
        if (source_end >= source_begin) then
          select case (trim(lower_ascii(cfg%particle_species(s)%source_mode)))
          case ('volume_seed')
            call sample_species_state( &
              cfg%particle_species(s), source_counts(s), &
              x_species(:, source_begin:source_end), v_species(:, source_begin:source_end), &
              temperature_k_override=effective_temperature_k(s), drift_velocity_override=effective_drift_velocity(:, s) &
              )
            counts_actual(s) = source_counts(s)
            w_species(source_begin:source_end) = batch_weight(s)
          case ('photo_raycast')
            if (.not. present(mesh)) then
              error stop 'photo_raycast requires mesh in init_particle_batch_from_config.'
            end if
            if (photo_emit_current_density(s) > 0.0_dp) then
              call sample_photo_species_state( &
                cfg%sim, cfg%particle_species(s), mesh, source_counts(s), x_species(:, source_begin:source_end), &
                v_species(:, source_begin:source_end), w_species(source_begin:source_end), counts_actual(s), &
                emit_elem_idx=emit_elem_species(source_begin:source_end), &
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
                  if (emit_elem_species(i) < 1_i32 .or. emit_elem_species(i) > mesh%nelem) then
                    error stop 'photo_raycast emitted invalid elem_idx.'
                  end if
                  if (present(photo_emission_dq)) then
                    photo_emission_dq(emit_elem_species(i)) = photo_emission_dq(emit_elem_species(i)) - &
                                                              cfg%particle_species(s)%q_particle*w_species(i)
                  end if
                  if (present(photo_emission_dq_by_species)) then
                    photo_emission_dq_by_species(emit_elem_species(i), s) = &
                      photo_emission_dq_by_species(emit_elem_species(i), s) - &
                      cfg%particle_species(s)%q_particle*w_species(i)
                  end if
                end do
              end if
            end if
          case default
            error stop 'Unknown '//species_config_path(int(s))//'.source.mode.'
          end select
        end if
        do face = 1, 6
          if (boundary_counts(face, s) <= 0_i32) cycle
          face_begin = counts_actual(s) + 1_i32
          face_end = face_begin + boundary_counts(face, s) - 1_i32
          call make_boundary_inflow_spec(cfg%sim, cfg%particle_species(s), face, face_spec)
          call sample_reservoir_face_particles( &
            cfg%sim%box_min, cfg%sim%box_max, face_spec%inject_face, face_spec%pos_low, face_spec%pos_high, &
            effective_drift_velocity(:, s), face_spec%m_particle, effective_temperature_k(s), cfg%sim%batch_duration, &
            x_species(:, face_begin:face_end), v_species(:, face_begin:face_end), &
            barrier_normal_energy=boundary_barrier(face, s), vmin_normal=boundary_vmin(face, s), &
            position_jitter_dt=cfg%sim%dt &
            )
          w_species(face_begin:face_end) = batch_weight(s)
          counts_actual(s) = face_end
        end do
      end associate
    end do

    call assemble_particle_batch(cfg%particle_species, counts_actual, sampled, pcls)
  end associate
  end procedure init_particle_batch_from_config

  !> 種別のサンプルを従来のround-robin順で完成したSoAへ直接詰める。
  subroutine assemble_particle_batch(species, counts, sampled, pcls)
    type(particle_species_spec), intent(in) :: species(:)
    integer(i32), intent(in) :: counts(:)
    type(sampled_species_type), intent(in) :: sampled(:)
    type(particles_soa), intent(out) :: pcls
    integer(i32) :: i, s, out_idx

    call allocate_particles(pcls, sum(counts))
    out_idx = 0_i32
    do i = 1, maxval(counts)
      do s = 1, size(counts)
        if (i > counts(s)) cycle
        out_idx = out_idx + 1_i32
        pcls%x(:, out_idx) = sampled(s)%x(:, i)
        pcls%v(:, out_idx) = sampled(s)%v(:, i)
        pcls%q(out_idx) = species(s)%q_particle
        pcls%m(out_idx) = species(s)%m_particle
        pcls%w(out_idx) = sampled(s)%w(i)
        pcls%species_id(out_idx) = s
        pcls%source_element(out_idx) = sampled(s)%source_element(i)
      end do
    end do
  end subroutine assemble_particle_batch

  !> batch injection の不完全な photo collision query を返し、status 未要求なら serial に停止する。
  subroutine finalize_particle_batch_collision_query( &
    query_status, batch_idx, species_idx, query_ray, query_bounce, &
    collision_failure_status, collision_failure_species, collision_failure_ray, collision_failure_bounce &
    )
    integer(i32), intent(in) :: query_status, batch_idx, species_idx, query_ray, query_bounce
    integer(i32), intent(out), optional :: collision_failure_status, collision_failure_species
    integer(i32), intent(out), optional :: collision_failure_ray, collision_failure_bounce
    character(len=256) :: failure_message

    if (present(collision_failure_status)) collision_failure_status = query_status
    if (present(collision_failure_species)) collision_failure_species = species_idx
    if (present(collision_failure_ray)) collision_failure_ray = query_ray
    if (present(collision_failure_bounce)) collision_failure_bounce = query_bounce
    if (present(collision_failure_status)) return

    write (failure_message, '(a,i0,a,i0,a,i0,a,i0,a,a,a,i0)') &
      'photo_raycast collision query incomplete during batch preparation: batch=', batch_idx, &
      ' species=', species_idx, ' ray=', query_ray, ' bounce=', query_bounce, &
      ' status=', trim(particle_status_name(query_status)), ' code=', query_status
    write (error_unit, '(a)') trim(failure_message)
    flush (error_unit)
    error stop 1
  end subroutine finalize_particle_batch_collision_query

end submodule bem_app_config_particle_runtime_batch

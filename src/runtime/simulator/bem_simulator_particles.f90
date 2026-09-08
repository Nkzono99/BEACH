!> バッチ粒子の生成と追跡を実装する。連成反復と受理判定は主ループが管理する。
submodule(bem_simulator) bem_simulator_particles
  implicit none
contains

  module procedure prepare_batch_state
  batch_idx = stats%batches + 1_i32
  call workspace%reset_before_injection()
  if (present(inject_state)) then
    call init_particle_batch_from_config( &
      app, batch_idx, pcls_batch, inject_state, mesh=mesh, &
      photo_emission_dq_by_species=workspace%photo_emission_dq, &
      mpi=mpi, collision_failure_status=collision_failure_status, &
      collision_failure_species=collision_failure_species, collision_failure_ray=collision_failure_ray, &
      collision_failure_bounce=collision_failure_bounce, snapshot=snapshot, source_plan=source_plan &
      )
  else
    call init_particle_batch_from_config( &
      app, batch_idx, pcls_batch, mesh=mesh, photo_emission_dq_by_species=workspace%photo_emission_dq, mpi=mpi, &
      collision_failure_status=collision_failure_status, collision_failure_species=collision_failure_species, &
      collision_failure_ray=collision_failure_ray, collision_failure_bounce=collision_failure_bounce, &
      snapshot=snapshot, source_plan=source_plan &
      )
  end if
  if (collision_failure_status /= collision_query_ok) return
  call workspace%prepare_particle_flags(pcls_batch%n)
  end procedure prepare_batch_state

  module procedure process_particle_batch
  integer(i32) :: i, step, tid, nth, collision_status, species_idx
  integer(i64) :: retry_attempted, retry_resolved
  real(dp) :: x0(3), v0(3), x1(3), v1(3), sampled_electric_field(3), qdep
  type(hit_info) :: hit
  type(particle_step_result) :: step_result, retry_result
  type(sim_config) :: particle_sim
  type(external_boundary_contract_type) :: particle_boundary_contract
  logical :: candidate_inside, used_event_resolver, adaptive_nonzero_mode, retry_field_available
!$ integer(kind=omp_sched_kind) :: previous_schedule_kind
!$ integer :: previous_schedule_chunk

  nth = size(dq_thread, 2)
  actual_team_size = nth
  collision_failure_status = collision_query_ok
  collision_failure_particle = huge(0_i32)
  collision_failure_step = 0_i32
  collision_failure_x = 0.0_dp
  collision_failure_v = 0.0_dp
  retry_attempted = 0_i64
  retry_resolved = 0_i64
  adaptive_nonzero_mode = app%periodic2%max_nonzero_mode_potential_step > 0.0_dp .or. &
                          trim(lower_ascii(app%surface_current%model)) == 'matching_plane_quasistatic'
  ! Replayed adaptive trials require an identical particle-index partition.
  ! Keep the normal runtime schedule, but override its ICV with static only for this adaptive loop.
!$ call omp_get_schedule(previous_schedule_kind, previous_schedule_chunk)
!$ if (adaptive_nonzero_mode) call omp_set_schedule(omp_sched_static, 0)

  !$omp parallel default(none) num_threads(nth) &
  !$omp shared(mesh,pcls_batch,app,boundary_contract,current_model,snapshot,dq_thread,bfield) &
  !$omp shared(escaped_boundary_flag,absorbed_flag,nth,actual_team_size) &
  !$omp shared(absorbed_element,soft_discarded_boundary_flag,batch_idx,mpi_rank) &
  !$omp shared(collision_failure_status,collision_failure_particle,collision_failure_step) &
  !$omp shared(collision_failure_x,collision_failure_v) &
  !$omp shared(matching_plane_moments_thread) &
  !$omp private(i,step,x0,v0,x1,v1,sampled_electric_field,hit,step_result,retry_result) &
  !$omp private(particle_sim,particle_boundary_contract,tid,qdep,species_idx) &
  !$omp private(collision_status,candidate_inside,used_event_resolver,retry_field_available) &
  !$omp reduction(+:retry_attempted,retry_resolved)
  tid = 1_i32
!$ tid = omp_get_thread_num() + 1
  !$omp single
!$ actual_team_size = int(omp_get_num_threads(), i32)
  !$omp end single
  !$omp do schedule(runtime)
  do i = 1_i32, pcls_batch%n
    if (.not. pcls_batch%alive(i)) cycle
    species_idx = pcls_batch%species_id(i)
    particle_sim = app%sim
    call resolve_particle_boundaries( &
      app%sim, app%particle_boundary_low, app%particle_boundary_high, app%particle_species(species_idx), &
      particle_sim%bc_low, particle_sim%bc_high &
      )
    particle_boundary_contract = boundary_contract
    call apply_species_kinetic_barrier(current_model, species_idx, particle_boundary_contract)
    do step = 1_i32, app%sim%max_step
      x0 = pcls_batch%x(:, i)
      v0 = pcls_batch%v(:, i)
      call build_particle_step_candidate( &
        mesh, particle_sim, snapshot, bfield, x0, v0, &
        pcls_batch%q(i), pcls_batch%m(i), app%sim%dt, x1, v1, sampled_electric_field &
        )
      if (.not. all(ieee_is_finite(x1)) .or. .not. all(ieee_is_finite(v1))) then
        call record_collision_failure(particle_step_invalid_boundary, i, step, x0, v0)
        exit
      end if
      call find_first_hit(mesh, x0, x1, hit, sim=particle_sim, status=collision_status)
      candidate_inside = .not. particle_sim%use_box .or. &
                         (all(x1 > particle_sim%box_min) .and. all(x1 < particle_sim%box_max))
      used_event_resolver = .false.
      if (collision_status /= collision_query_ok) then
        if (.not. candidate_inside) then
          call resolve_particle_boundary_candidate( &
            mesh, particle_sim, snapshot, bfield, x0, v0, pcls_batch%q(i), pcls_batch%m(i), app%sim%dt, x1, v1, &
            result=step_result, boundary_contract=particle_boundary_contract, &
            boundary_rng_counter=int([batch_idx, mpi_rank, i, step], i64), &
            sampled_electric_field=sampled_electric_field &
            )
          used_event_resolver = .true.
        else
          call record_collision_failure(collision_status, i, step, x0, v0)
          exit
        end if
      else if (candidate_inside) then
        if (hit%has_hit) then
          qdep = pcls_batch%q(i)*pcls_batch%w(i)
          dq_thread(hit%elem_idx, tid) = dq_thread(hit%elem_idx, tid) + qdep
          pcls_batch%alive(i) = .false.
          absorbed_flag(i) = .true.
          absorbed_element(i) = hit%elem_idx
          exit
        end if
        pcls_batch%x(:, i) = x1
        pcls_batch%v(:, i) = v1
      else
        call resolve_particle_boundary_candidate( &
          mesh, particle_sim, snapshot, bfield, x0, v0, pcls_batch%q(i), pcls_batch%m(i), app%sim%dt, x1, v1, &
          hit=hit, result=step_result, boundary_contract=particle_boundary_contract, &
          boundary_rng_counter=int([batch_idx, mpi_rank, i, step], i64), &
          sampled_electric_field=sampled_electric_field &
          )
        used_event_resolver = .true.
      end if
      if (used_event_resolver) then
        if (step_result%status == particle_step_multiple_box_events .and. &
            trim(lower_ascii(app%sim%multiple_box_events_retry_backend)) == 'upper_panel_fourier') then
          retry_attempted = retry_attempted + 1_i64
          call advance_particle_step_upper_panel_fourier( &
            mesh, particle_sim, snapshot, bfield, x0, v0, pcls_batch%q(i), pcls_batch%m(i), app%sim%dt, &
            retry_result, retry_field_available, boundary_contract=particle_boundary_contract, &
            boundary_rng_counter=int([batch_idx, mpi_rank, i, step], i64) &
            )
          if (retry_field_available .and. retry_result%status == collision_query_ok) then
            step_result = retry_result
            retry_resolved = retry_resolved + 1_i64
          end if
        end if
        if (step_result%status /= collision_query_ok) then
          if (step_result%status == particle_step_multiple_box_events .and. &
              trim(lower_ascii(app%sim%multiple_box_events_policy)) == 'soft_discard') then
            pcls_batch%alive(i) = .false.
            soft_discarded_boundary_flag(i) = .true.
            exit
          end if
          call record_collision_failure(step_result%status, i, step, x0, v0)
          exit
        end if
        matching_plane_moments_thread(1, species_idx, tid) = &
          matching_plane_moments_thread(1, species_idx, tid) + &
          pcls_batch%w(i)*real(step_result%z_high_outward_event_count, dp)
        matching_plane_moments_thread(2, species_idx, tid) = &
          matching_plane_moments_thread(2, species_idx, tid) + &
          pcls_batch%w(i)*step_result%z_high_outward_normal_kinetic_energy_j_sum
        matching_plane_moments_thread(3, species_idx, tid) = &
          matching_plane_moments_thread(3, species_idx, tid) + &
          pcls_batch%w(i)*real(step_result%outer_barrier_return_count, dp)
        matching_plane_moments_thread(4, species_idx, tid) = &
          matching_plane_moments_thread(4, species_idx, tid) + &
          pcls_batch%w(i)*real(step_result%outer_barrier_escape_count, dp)
        if (step_result%absorbed) then
          qdep = pcls_batch%q(i)*pcls_batch%w(i)
          dq_thread(step_result%elem_idx, tid) = dq_thread(step_result%elem_idx, tid) + qdep
          pcls_batch%alive(i) = .false.
          absorbed_flag(i) = .true.
          absorbed_element(i) = step_result%elem_idx
          exit
        end if
        if (step_result%escaped_boundary) then
          pcls_batch%alive(i) = .false.
          escaped_boundary_flag(i) = .true.
          exit
        end if
        pcls_batch%x(:, i) = step_result%x
        pcls_batch%v(:, i) = step_result%v
      end if
    end do
  end do
  !$omp end do
  !$omp end parallel
!$ call omp_set_schedule(previous_schedule_kind, previous_schedule_chunk)
  retry_counts = [retry_attempted, retry_resolved]

contains
  subroutine apply_species_kinetic_barrier(current_model, species_idx, contract)
    type(surface_closure_contract_type), intent(in) :: current_model
    integer(i32), intent(in) :: species_idx
    type(external_boundary_contract_type), intent(inout) :: contract

    integer(i32) :: axis, face

    if (.not. current_model%active) return
    if (.not. allocated(current_model%has_outflow_kinetic_barrier)) return
    if (species_idx < 1_i32 .or. species_idx > size(current_model%has_outflow_kinetic_barrier)) return
    if (.not. current_model%has_outflow_kinetic_barrier(species_idx)) return
    face = current_model%outflow_barrier_face(species_idx)
    if (face < 1_i32 .or. face > 6_i32) return
    axis = (face + 1_i32)/2_i32
    if (mod(face, 2_i32) == 0_i32) then
      contract%barrier_override_high(axis) = .true.
      contract%barrier_potential_high_v(axis) = current_model%outflow_barrier_potential_v(species_idx)
    else
      contract%barrier_override_low(axis) = .true.
      contract%barrier_potential_low_v(axis) = current_model%outflow_barrier_potential_v(species_idx)
    end if
  end subroutine apply_species_kinetic_barrier

  subroutine record_collision_failure(status, particle_index, particle_step, failure_x, failure_v)
    integer(i32), intent(in) :: status, particle_index, particle_step
    real(dp), intent(in) :: failure_x(3), failure_v(3)
    !$omp critical (beach_collision_query_failure)
    if (collision_failure_status == collision_query_ok .or. particle_index < collision_failure_particle .or. &
        (particle_index == collision_failure_particle .and. particle_step < collision_failure_step)) then
      collision_failure_status = status
      collision_failure_particle = particle_index
      collision_failure_step = particle_step
      collision_failure_x = failure_x
      collision_failure_v = failure_v
    end if
    !$omp end critical (beach_collision_query_failure)
  end subroutine record_collision_failure
  end procedure process_particle_batch

end submodule bem_simulator_particles

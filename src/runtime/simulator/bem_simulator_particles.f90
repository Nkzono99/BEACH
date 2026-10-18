!> バッチ粒子の生成と追跡を実装する。連成反復と受理判定は主ループが管理する。
submodule(bem_simulator) bem_simulator_particles
  use bem_constants, only: qe, k_boltzmann
  use bem_secondary_emission
  use bem_config_helpers, only: species_temperature_k
  use bem_photoelectron_injection, only: sample_photo_emission_velocity, build_tangent_basis, &
                                         canonicalize_periodic_emission_position
  use bem_particles, only: append_particles
  use bem_outflow_spectrum, only: outflow_spectrum_bin
  use bem_particle_stepper, only: particle_substep_needs_halving
  implicit none
contains

  module procedure prepare_batch_state
  batch_idx = stats%batches + 1_i32
  call workspace%reset_before_injection()
  call init_particle_batch_from_config( &
    app, batch_idx, pcls_batch, inject_state, mesh=mesh, &
    photo_emission_dq_by_species=workspace%photo_emission_dq, &
    mpi=mpi, collision_failure_status=collision_failure_status, &
    collision_failure_species=collision_failure_species, collision_failure_ray=collision_failure_ray, &
    collision_failure_bounce=collision_failure_bounce, snapshot=snapshot, source_plan=source_plan &
    )
  if (collision_failure_status /= collision_query_ok) return
  call workspace%prepare_particle_flags(pcls_batch%n)
  end procedure prepare_batch_state

  module procedure process_particle_batch
  integer(i32) :: i, step, tid, nth, collision_status, species_idx, spectrum_row
  integer(i64) :: retry_attempted, retry_resolved, halved_steps, substeps_total, depth_limited_steps, rng_step
  integer(i32) :: sub_index
  real(dp) :: h, h_prev, h_min, remaining, field_ref(3)
  logical :: adaptive_substeps, have_field_ref, particle_done, halved, depth_limited
  integer(i32) :: first_particle, last_particle, generation
  type(secondary_spawn_buffer_type), allocatable :: spawn(:)
  real(dp) :: x0(3), v0(3), x1(3), v1(3), sampled_electric_field(3), qdep
  type(hit_info) :: hit
  type(particle_step_result) :: step_result, retry_result
  type(sim_config), target :: species_sim(app%n_particle_species)
  type(external_boundary_contract_type), target :: species_boundary_contract(app%n_particle_species)
  type(sim_config), pointer :: particle_sim
  type(external_boundary_contract_type), pointer :: particle_boundary_contract
  logical :: candidate_inside, used_event_resolver, adaptive_nonzero_mode, retry_field_available
!$ integer(kind=omp_sched_kind) :: previous_schedule_kind
!$ integer :: previous_schedule_chunk

  nth = size(workspace%dq_thread, 2)
  actual_team_size = nth
  collision_failure_status = collision_query_ok
  collision_failure_particle = huge(0_i32)
  collision_failure_step = 0_i32
  collision_failure_x = 0.0_dp
  collision_failure_v = 0.0_dp
  retry_attempted = 0_i64
  retry_resolved = 0_i64
  halved_steps = 0_i64
  substeps_total = 0_i64
  depth_limited_steps = 0_i64
  adaptive_substeps = app%sim%substep_field_change_max > 0.0_dp
  h_min = app%sim%dt*0.5_dp**app%sim%substep_max_halvings
  adaptive_nonzero_mode = app%periodic2%max_nonzero_mode_potential_step > 0.0_dp
  ! Boundary settings and kinetic barriers are fixed throughout this trial batch.
  do species_idx = 1_i32, app%n_particle_species
    species_sim(species_idx) = app%sim
    call resolve_particle_boundaries( &
      app%sim, app%particle_boundary_low, app%particle_boundary_high, app%particle_species(species_idx), &
      species_sim(species_idx)%bc_low, species_sim(species_idx)%bc_high &
      )
    species_boundary_contract(species_idx) = boundary_contract
    call apply_species_kinetic_barrier(current_model, species_idx, species_boundary_contract(species_idx))
  end do
  ! Replayed adaptive trials require an identical particle-index partition.
  ! Keep the normal runtime schedule, but override its ICV with static only for this adaptive loop.
!$ call omp_get_schedule(previous_schedule_kind, previous_schedule_chunk)
!$ if (adaptive_nonzero_mode) call omp_set_schedule(omp_sched_static, 0)

  allocate (spawn(nth))
  first_particle = 1_i32
  last_particle = pcls_batch%n
  do generation = 0_i32, app%secondary_emission_max_generations
    spawn%n = 0_i32
    !$omp parallel default(none) num_threads(nth) &
    !$omp shared(cell,mesh,pcls_batch,app,species_sim,species_boundary_contract,snapshot,workspace,bfield,spawn) &
    !$omp shared(nth,actual_team_size,batch_idx,mpi_rank,first_particle,last_particle,generation) &
    !$omp shared(adaptive_substeps,h_min) &
    !$omp shared(collision_failure_status,collision_failure_particle,collision_failure_step) &
    !$omp shared(collision_failure_x,collision_failure_v) &
    !$omp private(i,step,x0,v0,x1,v1,sampled_electric_field,hit,step_result,retry_result) &
    !$omp private(particle_sim,particle_boundary_contract,tid,qdep,species_idx,spectrum_row) &
    !$omp private(collision_status,candidate_inside,used_event_resolver,retry_field_available) &
    !$omp private(h,h_prev,remaining,field_ref,have_field_ref,particle_done,halved,depth_limited,sub_index,rng_step) &
    !$omp reduction(+:retry_attempted,retry_resolved,halved_steps,substeps_total,depth_limited_steps)
    tid = 1_i32
!$  tid = omp_get_thread_num() + 1
    !$omp single
!$  actual_team_size = int(omp_get_num_threads(), i32)
    !$omp end single
    !$omp do schedule(runtime)
    do i = first_particle, last_particle
      if (.not. pcls_batch%alive(i)) cycle
      species_idx = pcls_batch%species_id(i)
      particle_sim => species_sim(species_idx)
      particle_boundary_contract => species_boundary_contract(species_idx)
      have_field_ref = .false.
      h_prev = app%sim%dt
      step_loop: do step = 1_i32, app%sim%max_step
        if (adaptive_substeps .and. .not. have_field_ref) then
          ! The first reference field is the field where the particle starts its track.
          call snapshot%eval_local_e(mesh, pcls_batch%x(:, i), field_ref)
          have_field_ref = .true.
        end if
        remaining = app%sim%dt
        sub_index = 0_i32
        halved = .false.
        depth_limited = .false.
        particle_done = .false.
        ! Without adaptive sub-steps this loop runs once with h = dt.
        sub_loop: do while (remaining > 0.0_dp)
          x0 = pcls_batch%x(:, i)
          v0 = pcls_batch%v(:, i)
          h = remaining
          if (adaptive_substeps) h = min(remaining, 2.0_dp*h_prev)
          halve_loop: do
            call build_particle_step_candidate( &
              mesh, particle_sim, snapshot, bfield, x0, v0, &
              pcls_batch%q(i), pcls_batch%m(i), h, x1, v1, sampled_electric_field &
              )
            if (.not. adaptive_substeps) exit halve_loop
            if (.not. all(ieee_is_finite(x1)) .or. .not. all(ieee_is_finite(v1))) exit halve_loop
            if (.not. particle_substep_needs_halving(sampled_electric_field, field_ref, pcls_batch%q(i), x0, x1, &
                                                     app%sim%substep_field_change_max)) exit halve_loop
            if (0.5_dp*h < h_min) then
              depth_limited = .true.
              exit halve_loop
            end if
            h = 0.5_dp*h
            halved = .true.
          end do halve_loop
          rng_step = int(step, i64)
          if (adaptive_substeps) rng_step = int(step, i64)*65536_i64 + int(sub_index, i64)
          if (.not. all(ieee_is_finite(x1)) .or. .not. all(ieee_is_finite(v1))) then
            call record_collision_failure(particle_step_invalid_boundary, i, step, x0, v0)
            particle_done = .true.
            exit sub_loop
          end if
          call find_first_hit(mesh, x0, x1, hit, sim=particle_sim, status=collision_status)
          candidate_inside = .not. particle_sim%use_box .or. &
                             (all(x1 > particle_sim%box_min) .and. all(x1 < particle_sim%box_max))
          used_event_resolver = .false.
          if (collision_status /= collision_query_ok) then
            if (.not. candidate_inside) then
              call resolve_particle_boundary_candidate( &
                mesh, particle_sim, snapshot, bfield, x0, v0, pcls_batch%q(i), pcls_batch%m(i), h, x1, v1, &
                result=step_result, boundary_contract=particle_boundary_contract, &
                boundary_rng_counter=[int(batch_idx, i64), int(mpi_rank, i64), int(i, i64), rng_step], &
                sampled_electric_field=sampled_electric_field, &
                cell=cell &
                )
              used_event_resolver = .true.
            else
              call record_collision_failure(collision_status, i, step, x0, v0)
              particle_done = .true.
              exit sub_loop
            end if
          else if (candidate_inside) then
            if (hit%has_hit) then
              qdep = pcls_batch%q(i)*pcls_batch%w(i)
              workspace%dq_thread(hit%elem_idx, tid) = workspace%dq_thread(hit%elem_idx, tid) + qdep
              pcls_batch%alive(i) = .false.
              workspace%absorbed_flag(i) = .true.
              workspace%absorbed_element(i) = hit%elem_idx
              call emit_secondaries(i, step, tid, hit%elem_idx, hit%pos, v0 + hit%t*(v1 - v0))
              particle_done = .true.
              exit sub_loop
            end if
            pcls_batch%x(:, i) = x1
            pcls_batch%v(:, i) = v1
          else
            call resolve_particle_boundary_candidate( &
              mesh, particle_sim, snapshot, bfield, x0, v0, pcls_batch%q(i), pcls_batch%m(i), h, x1, v1, &
              hit=hit, result=step_result, boundary_contract=particle_boundary_contract, &
              boundary_rng_counter=[int(batch_idx, i64), int(mpi_rank, i64), int(i, i64), rng_step], &
              sampled_electric_field=sampled_electric_field, &
              cell=cell &
              )
            used_event_resolver = .true.
          end if
          if (used_event_resolver) then
            if (step_result%status == particle_step_multiple_box_events .and. &
                trim(lower_ascii(app%sim%multiple_box_events_retry_backend)) == 'upper_panel_fourier') then
              retry_attempted = retry_attempted + 1_i64
              call advance_particle_step_upper_panel_fourier( &
                mesh, particle_sim, snapshot, bfield, x0, v0, pcls_batch%q(i), pcls_batch%m(i), h, &
                retry_result, retry_field_available, boundary_contract=particle_boundary_contract, &
                boundary_rng_counter=[int(batch_idx, i64), int(mpi_rank, i64), int(i, i64), rng_step], &
                cell=cell &
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
                workspace%soft_discarded_boundary_flag(i) = .true.
                particle_done = .true.
                exit sub_loop
              end if
              call record_collision_failure(step_result%status, i, step, x0, v0)
              particle_done = .true.
              exit sub_loop
            end if
            workspace%matching_plane_moments_thread(1, species_idx, tid) = &
              workspace%matching_plane_moments_thread(1, species_idx, tid) + &
              pcls_batch%w(i)*real(step_result%wall_outward_event_count, dp)
            workspace%matching_plane_moments_thread(2, species_idx, tid) = &
              workspace%matching_plane_moments_thread(2, species_idx, tid) + &
              pcls_batch%w(i)*step_result%wall_outward_normal_kinetic_energy_j_sum
            workspace%matching_plane_moments_thread(3, species_idx, tid) = &
              workspace%matching_plane_moments_thread(3, species_idx, tid) + &
              pcls_batch%w(i)*real(step_result%outer_barrier_return_count, dp)
            workspace%matching_plane_moments_thread(4, species_idx, tid) = &
              workspace%matching_plane_moments_thread(4, species_idx, tid) + &
              pcls_batch%w(i)*real(step_result%outer_barrier_escape_count, dp)
            if (step_result%wall_outward_event_count > 0) then
              ! A step rarely crosses outward twice; its crossings share the mean normal energy bin.
              spectrum_row = 4_i32 + outflow_spectrum_bin( &
                             step_result%wall_outward_normal_kinetic_energy_j_sum/ &
                             (real(step_result%wall_outward_event_count, dp)*qe))
              workspace%matching_plane_moments_thread(spectrum_row, species_idx, tid) = &
                workspace%matching_plane_moments_thread(spectrum_row, species_idx, tid) + &
                pcls_batch%w(i)*real(step_result%wall_outward_event_count, dp)
            end if
            if (step_result%absorbed) then
              qdep = pcls_batch%q(i)*pcls_batch%w(i)
              workspace%dq_thread(step_result%elem_idx, tid) = workspace%dq_thread(step_result%elem_idx, tid) + qdep
              pcls_batch%alive(i) = .false.
              workspace%absorbed_flag(i) = .true.
              workspace%absorbed_element(i) = step_result%elem_idx
              call emit_secondaries(i, step, tid, step_result%elem_idx, step_result%x, step_result%v)
              particle_done = .true.
              exit sub_loop
            end if
            if (step_result%escaped_boundary) then
              pcls_batch%alive(i) = .false.
              workspace%escaped_boundary_flag(i) = .true.
              particle_done = .true.
              exit sub_loop
            end if
            pcls_batch%x(:, i) = step_result%x
            pcls_batch%v(:, i) = step_result%v
          end if
          ! Accepted: the sampled midpoint field becomes the reference for the next sub-step.
          field_ref = sampled_electric_field
          h_prev = h
          remaining = remaining - h
          sub_index = sub_index + 1_i32
        end do sub_loop
        if (halved) then
          halved_steps = halved_steps + 1_i64
          substeps_total = substeps_total + int(sub_index, i64)
          if (depth_limited) depth_limited_steps = depth_limited_steps + 1_i64
        end if
        if (particle_done) exit step_loop
      end do step_loop
    end do
    !$omp end do
    !$omp end parallel
    if (collision_failure_status /= collision_query_ok) exit
    if (sum(spawn%n) == 0_i32) exit
    call append_spawned_particles(first_particle, last_particle)
  end do
!$ call omp_set_schedule(previous_schedule_kind, previous_schedule_chunk)
  retry_counts = [retry_attempted, retry_resolved]
  substep_counts = [halved_steps, substeps_total, depth_limited_steps]

contains
  subroutine emit_secondaries(parent, particle_step, thread, element, position, velocity)
    integer(i32), intent(in) :: parent, particle_step, thread, element
    real(dp), intent(in) :: position(3), velocity(3)
    type(secondary_spawn_type) :: child
    integer(i32) :: incident, target, count, ordinal
    integer(i64) :: counter(7)
    real(dp) :: normal(3), tangent1(3), tangent2(3), speed, energy_ev, cosine, delta, sigma, charge

    incident = pcls_batch%species_id(parent)
    associate (spec => app%particle_species(incident))
      if (.not. spec%secondary_emission_enabled) return
      target = spec%secondary_emission_target_idx
      normal = mesh%normals(:, element)
      if (dot_product(velocity, normal) > 0.0_dp) normal = -normal
      speed = sqrt(sum(velocity**2))
      cosine = 1.0_dp
      if (speed > 0.0_dp) cosine = -dot_product(velocity, normal)/speed
      energy_ev = 0.5_dp*pcls_batch%m(parent)*speed**2/qe
      delta = secondary_emission_yield(spec%secondary_emission_model, energy_ev, cosine, &
                                       spec%secondary_emission_yield_max, spec%secondary_emission_energy_max_ev, &
                                       spec%secondary_emission_const_yield)
    end associate
    counter = int([batch_idx, mpi_rank, parent, particle_step, generation, 0_i32, 271_i32], i64)
    count = secondary_emission_count(delta, abs(pcls_batch%q(parent)/app%particle_species(target)%q_particle), &
                                     app%sim%rng_seed, counter)
    workspace%secondary_counts_thread(1, incident, thread) = workspace%secondary_counts_thread(1, incident, thread) + 1_i64
    workspace%secondary_values_thread(2, incident, thread) = workspace%secondary_values_thread(2, incident, thread) + delta
    if (generation >= app%secondary_emission_max_generations) then
      workspace%secondary_counts_thread(3, incident, thread) = &
        workspace%secondary_counts_thread(3, incident, thread) + int(count, i64)
      return
    end if
    workspace%secondary_counts_thread(2, incident, thread) = &
      workspace%secondary_counts_thread(2, incident, thread) + int(count, i64)
    charge = app%particle_species(target)%q_particle*pcls_batch%w(parent)
    workspace%secondary_values_thread(1, incident, thread) = &
      workspace%secondary_values_thread(1, incident, thread) + real(count, dp)*charge
    call build_tangent_basis(normal, tangent1, tangent2)
    sigma = sqrt(k_boltzmann*species_temperature_k(app%particle_species(target))/app%particle_species(target)%m_particle)
    child%x = position + 1.0e-12_dp*normal
    call canonicalize_periodic_emission_position(app%sim, child%x)
    child%w = pcls_batch%w(parent)
    child%species_idx = target
    child%element = element
    child%parent = parent
    do ordinal = 1, count
      child%ordinal = ordinal
      counter(6) = int(ordinal, i64)
      call sample_photo_emission_velocity(sigma, 0.0_dp, normal, tangent1, tangent2, child%v, &
                                          rng_seed=app%sim%rng_seed, rng_counter=counter)
      call spawn(thread)%push(child)
    end do
  end subroutine emit_secondaries

  !> Append in parent/ordinal order independently of the OpenMP schedule. The sparse spawn
  !> records also carry reaction charge into the existing surface-emission channel.
  subroutine append_spawned_particles(first, last)
    integer(i32), intent(inout) :: first, last
    integer(i32), allocatable :: starts(:), counts(:), species(:), elements(:)
    real(dp), allocatable :: x(:, :), v(:, :), q(:), m(:), w(:)
    integer(i32) :: thread, j, parent, index, added, old_count
    type(secondary_spawn_type) :: child

    added = sum(spawn%n)
    allocate (starts(first:last), counts(first:last))
    counts = 0_i32
    do thread = 1, nth
      do j = 1, spawn(thread)%n
        parent = spawn(thread)%particles(j)%parent
        counts(parent) = counts(parent) + 1_i32
      end do
    end do
    index = 0_i32
    do parent = first, last
      starts(parent) = index
      index = index + counts(parent)
    end do
    allocate (x(3, added), v(3, added), q(added), m(added), w(added), species(added), elements(added))
    do thread = 1, nth
      do j = 1, spawn(thread)%n
        child = spawn(thread)%particles(j)
        index = starts(child%parent) + child%ordinal
        x(:, index) = child%x
        v(:, index) = child%v
        w(index) = child%w
        species(index) = child%species_idx
        elements(index) = child%element
        q(index) = app%particle_species(child%species_idx)%q_particle
        m(index) = app%particle_species(child%species_idx)%m_particle
      end do
    end do
    do index = 1, added
      workspace%photo_emission_dq(elements(index), species(index)) = &
        workspace%photo_emission_dq(elements(index), species(index)) - q(index)*w(index)
    end do
    old_count = pcls_batch%n
    call append_particles(pcls_batch, x, v, q, m, w, species, elements)
    call workspace%extend_particle_flags(old_count, pcls_batch%n)
    first = old_count + 1_i32
    last = pcls_batch%n
  end subroutine append_spawned_particles

  subroutine apply_species_kinetic_barrier(current_model, species_idx, contract)
    type(surface_closure_contract_type), intent(in) :: current_model
    integer(i32), intent(in) :: species_idx
    type(external_boundary_contract_type), intent(inout) :: contract

    integer(i32) :: axis, face

    if (.not. current_model%active) return
    if (.not. current_model%has_outflow_kinetic_barrier(species_idx)) return
    face = current_model%outflow_barrier_face(species_idx)
    axis = (face + 1_i32)/2_i32
    if (mod(face, 2_i32) == 0_i32) then
      contract%barrier_override_high(axis) = .true.
      contract%barrier_potential_high_v(axis) = current_model%outflow_barrier_potential_v(species_idx)
      contract%barrier_redistribute_high(axis) = current_model%outer_return_cell_uniform
    else
      contract%barrier_override_low(axis) = .true.
      contract%barrier_potential_low_v(axis) = current_model%outflow_barrier_potential_v(species_idx)
      contract%barrier_redistribute_low(axis) = current_model%outer_return_cell_uniform
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

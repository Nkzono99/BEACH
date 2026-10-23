!> バッチ粒子の生成と追跡を実装する。連成反復と受理判定は主ループが管理する。
submodule(bem_simulator) bem_simulator_particles
  use bem_constants, only: qe, k_boltzmann
  use bem_secondary_emission
  use bem_config_helpers, only: species_temperature_k
  use bem_photoelectron_injection, only: sample_photo_emission_velocity, build_tangent_basis, &
                                         canonicalize_periodic_emission_position
  use bem_particles, only: append_particles, allocate_particles, trim_particles
  use bem_outflow_spectrum, only: outflow_spectrum_bin
  use bem_particle_stepper, only: particle_substep_needs_halving
  use bem_counter_rng, only: counter_child_id
  implicit none
contains

  module procedure prepare_batch_state
  integer(i32), allocatable :: slot_ids(:)

  batch_idx = stats%batches + 1_i32
  call workspace%reset_before_injection()
  if (trim(app%sim%particle_distribution) == 'shared') then
    ! The pool mode only fixes the global layout here; slots are generated while they are claimed.
    call build_batch_injection_layout( &
      app, batch_idx, injection_layout, state=inject_state, mesh=mesh, mpi=mpi, snapshot=snapshot, &
      source_plan=source_plan &
      )
    call allocate_particles(pcls_batch, injection_layout%total_slots)
    pcls_batch%alive = .false.
    pcls_batch%n = 0_i32
    call workspace%prepare_particle_flags(injection_layout%total_slots)
    return
  end if
  injection_layout = batch_injection_layout_type()
  call init_particle_batch_from_config( &
    app, batch_idx, pcls_batch, inject_state, mesh=mesh, &
    photo_emission_dq_by_species=workspace%photo_emission_dq, &
    mpi=mpi, photo_failure=photo_failure, snapshot=snapshot, source_plan=source_plan, &
    slot_ids=slot_ids &
    )
  if (photo_failure%status /= collision_query_ok) return
  call workspace%prepare_particle_flags(pcls_batch%n)
  ! Boundary and secondary-emission draws are keyed on the global slot, as in the shared pool.
  workspace%rng_id(1:pcls_batch%n) = int(slot_ids, i64)
  end procedure prepare_batch_state

  module procedure process_particle_batch
  integer(i32) :: i, tid, nth, species_idx, slot, local_index, claimed_count, rank_next_slot, rank_end_slot
  integer(i64) :: thread_step_counts(6, size(workspace%dq_thread, 2))
  integer(i64) :: pool_first, pool_chunk
  real(dp) :: h_min, claim_t0, claim_t1
  logical :: adaptive_substeps, adaptive_nonzero_mode, shared_pool, pool_exhausted, step_growth
  integer(i32) :: growth_chunk_max
  real(dp) :: growth_change_max
  !> step 伸長時の仕事の床値は運動エネルギーの `step_growth_energy_fraction` 倍（1e-4 eV との大きい方）。遠方の弱い場では
  !! FMM 打ち切り雑音が場そのものと同程度なので、絶対床値だけでは雑音が半分化を起こす。
  real(dp) :: growth_energy_fraction
  !> 半分化が起きた後、この回数の chunk は伸長を試みない（伸長と半分化の往復を防ぐ）。
  integer(i32), parameter :: growth_hold_chunks = 8_i32
  logical :: retry_upper_panel_fourier, soft_discard_multiple_box_events
  integer(i32) :: first_particle, last_particle, generation, rng_tag
  type(secondary_spawn_buffer_type), allocatable :: spawn(:)
  type(particle_domain_config), target :: species_sim(app%n_particle_species)
  type(external_boundary_contract_type), target :: species_boundary_contract(app%n_particle_species)
  real(dp) :: slot_x(3), slot_v(3), slot_w, slot_q, slot_m
  integer(i32) :: slot_species, slot_element, slot_status, slot_bounce
  logical :: slot_emitted
!$ integer(kind=omp_sched_kind) :: previous_schedule_kind
!$ integer :: previous_schedule_chunk

  nth = size(workspace%dq_thread, 2)
  actual_team_size = nth
  thread_step_counts = 0_i64
  pool_claims = 0_i64
  pool_claim_seconds = 0.0_dp
  adaptive_substeps = app%sim%substep_field_change_max > 0.0_dp
  h_min = app%sim%dt*0.5_dp**app%sim%substep_max_halvings
  ! Step growth doubles a particle's step while the field change stays below a quarter of the halving threshold.
  step_growth = app%sim%step_growth .and. adaptive_substeps
  growth_chunk_max = 1_i32
  if (step_growth) growth_chunk_max = 2_i32**min(16_i32, app%sim%step_growth_max)
  growth_change_max = 0.25_dp*app%sim%substep_field_change_max
  growth_energy_fraction = app%sim%step_growth_energy_fraction
  adaptive_nonzero_mode = app%periodic2%max_nonzero_mode_potential_step > 0.0_dp
  shared_pool = trim(app%sim%particle_distribution) == 'shared'
  ! multiple_box_events の再試行 backend と policy はバッチの先頭で一度だけ解決する。
  retry_upper_panel_fourier = trim(lower_ascii(app%sim%multiple_box_events_retry_backend)) == 'upper_panel_fourier'
  soft_discard_multiple_box_events = trim(lower_ascii(app%sim%multiple_box_events_policy)) == 'soft_discard'
  ! Boundary settings and kinetic barriers are fixed throughout this trial batch.
  do species_idx = 1_i32, app%n_particle_species
    species_sim(species_idx) = app%sim%particle_domain_config
    call resolve_particle_boundaries( &
      app%sim, app%particle_boundary_low, app%particle_boundary_high, app%particle_species(species_idx), &
      species_sim(species_idx)%bc_low, species_sim(species_idx)%bc_high &
      )
    species_boundary_contract(species_idx) = boundary_contract
    call apply_species_kinetic_barrier(current_model, species_idx, species_boundary_contract(species_idx))
  end do
  if (shared_pool) then
    if (.not. injection_layout%ready) error stop 'shared particle pool requires a batch injection layout.'
    if (.not. pool_counter%active) error stop 'shared particle pool counter is not initialized.'
    if (size(pcls_batch%q) < injection_layout%total_slots) then
      error stop 'shared particle pool batch capacity is smaller than the layout.'
    end if
    pool_chunk = int(max(1_i32, app%sim%particle_pool_chunk), i64)
  end if
  ! Replayed adaptive trials require an identical particle-index partition.
  ! Keep the normal runtime schedule, but override its ICV with static only for this adaptive loop.
!$ call omp_get_schedule(previous_schedule_kind, previous_schedule_chunk)
!$ if (adaptive_nonzero_mode) call omp_set_schedule(omp_sched_static, 0)

  allocate (spawn(nth))
  first_particle = 1_i32
  last_particle = pcls_batch%n
  claimed_count = 0_i32
  do generation = 0_i32, app%secondary_emission_max_generations
    spawn%n = 0_i32
    ! Boundary and emission counters are keyed on (generation, global slot), so a batch does not depend on
    ! the rank count or on which thread tracks a particle.
    rng_tag = generation
    if (shared_pool .and. generation == 0_i32) then
      ! The rank claims chunks from the global counter, and its threads take single slots from the
      ! current chunk, so thread-level balance stays as fine as schedule(dynamic,1).
      rank_next_slot = 1_i32
      rank_end_slot = 0_i32
      pool_exhausted = .false.
      !$omp parallel default(none) num_threads(nth) &
      !$omp shared(app,mesh,injection_layout,pcls_batch,workspace,pool_counter,pool_chunk,claimed_count) &
      !$omp shared(nth,actual_team_size,pool_claims,pool_claim_seconds,rank_next_slot,rank_end_slot,pool_exhausted) &
      !$omp private(tid,slot,local_index,pool_first,claim_t0,claim_t1) &
      !$omp private(slot_x,slot_v,slot_w,slot_q,slot_m,slot_species,slot_element,slot_status,slot_bounce,slot_emitted)
      tid = 1_i32
!$    tid = omp_get_thread_num() + 1
      !$omp single
!$    actual_team_size = int(omp_get_num_threads(), i32)
      !$omp end single
      claim_loop: do
        ! The critical section also serializes the MPI call for MPI_THREAD_SERIALIZED.
        !$omp critical (beach_particle_pool_claim)
        if (rank_next_slot > rank_end_slot .and. .not. pool_exhausted) then
          claim_t0 = wall_seconds()
          call pool_counter%claim(pool_chunk, pool_first)
          claim_t1 = wall_seconds()
          pool_claims = pool_claims + 1_i64
          pool_claim_seconds = pool_claim_seconds + (claim_t1 - claim_t0)
          if (pool_first >= int(injection_layout%total_slots, i64)) then
            pool_exhausted = .true.
          else
            rank_next_slot = int(pool_first, i32) + 1_i32
            rank_end_slot = int(min(pool_first + pool_chunk, int(injection_layout%total_slots, i64)), i32)
          end if
        end if
        if (pool_exhausted) then
          slot = injection_layout%total_slots + 1_i32
        else
          slot = rank_next_slot
          rank_next_slot = rank_next_slot + 1_i32
          ! The window owner serves the other ranks' claims whenever one of its threads passes here.
          call pool_counter%progress()
        end if
        !$omp end critical (beach_particle_pool_claim)
        if (slot > injection_layout%total_slots) exit claim_loop
        call generate_batch_slot( &
          app, injection_layout, slot, mesh, slot_x, slot_v, slot_w, slot_q, slot_m, slot_species, slot_element, &
          slot_emitted, slot_status, slot_bounce &
          )
        if (slot_status /= collision_query_ok) then
          call record_collision_failure(slot_status, slot, 0_i32, slot_x, slot_v)
          cycle
        end if
        if (.not. slot_emitted) cycle
        !$omp atomic capture
        claimed_count = claimed_count + 1_i32
        local_index = claimed_count
        !$omp end atomic
        pcls_batch%x(:, local_index) = slot_x
        pcls_batch%v(:, local_index) = slot_v
        pcls_batch%q(local_index) = slot_q
        pcls_batch%m(local_index) = slot_m
        pcls_batch%w(local_index) = slot_w
        pcls_batch%species_id(local_index) = slot_species
        pcls_batch%source_element(local_index) = slot_element
        pcls_batch%alive(local_index) = .true.
        workspace%rng_id(local_index) = int(slot, i64)
        if (slot_element > 0_i32) then
          if (app%particle_species(slot_species)%deposit_opposite_charge_on_emit) then
            !$omp atomic update
            workspace%photo_emission_dq(slot_element, slot_species) = &
              workspace%photo_emission_dq(slot_element, slot_species) - slot_q*slot_w
          end if
        end if
        call track_particle(local_index, tid)
      end do claim_loop
      !$omp end parallel
      pcls_batch%n = claimed_count
      call trim_particles(pcls_batch, claimed_count)
      first_particle = 1_i32
      last_particle = claimed_count
    else
      !$omp parallel default(none) num_threads(nth) &
      !$omp shared(pcls_batch,nth,actual_team_size,first_particle,last_particle) &
      !$omp private(i,tid)
      tid = 1_i32
!$    tid = omp_get_thread_num() + 1
      !$omp single
!$    actual_team_size = int(omp_get_num_threads(), i32)
      !$omp end single
      !$omp do schedule(runtime)
      do i = first_particle, last_particle
        if (.not. pcls_batch%alive(i)) cycle
        call track_particle(i, tid)
      end do
      !$omp end do
      !$omp end parallel
    end if
    if (track_failure%status /= collision_query_ok) exit
    if (sum(spawn%n) == 0_i32) exit
    call append_spawned_particles(first_particle, last_particle)
  end do
!$ call omp_set_schedule(previous_schedule_kind, previous_schedule_chunk)
  retry_counts = [sum(thread_step_counts(1, :)), sum(thread_step_counts(2, :))]
  substep_counts = [sum(thread_step_counts(3, :)), sum(thread_step_counts(4, :)), sum(thread_step_counts(5, :)), &
                    sum(thread_step_counts(6, :))]

contains
  !> 粒子 `i` を吸収・脱出・最大 step まで追跡する。OpenMP 並列領域から呼ぶので、状態は引数と局所変数だけに持つ。
  subroutine track_particle(i, tid)
    integer(i32), intent(in) :: i, tid
    integer(i32) :: step, collision_status, species_idx, spectrum_row, sub_index, chunk_steps, growth_hold
    integer(i64) :: rng_step, boundary_rng_counter(4)
    real(dp) :: h, h_prev, remaining, field_ref(3), energy_floor
    logical :: have_field_ref, particle_done, halved, depth_limited, grow_ok
    real(dp) :: x0(3), v0(3), x1(3), v1(3), sampled_electric_field(3), qdep
    type(hit_info) :: hit
    type(particle_step_result) :: step_result, retry_result
    type(particle_domain_config), pointer :: particle_sim
    type(external_boundary_contract_type), pointer :: particle_boundary_contract
    logical :: candidate_inside, used_event_resolver, retry_field_available

    species_idx = pcls_batch%species_id(i)
    particle_sim => species_sim(species_idx)
    particle_boundary_contract => species_boundary_contract(species_idx)
    have_field_ref = .false.
    h_prev = app%sim%dt
    grow_ok = .false.
    growth_hold = 0_i32
    energy_floor = 0.0_dp
    step = 1_i32
    step_loop: do while (step <= app%sim%max_step)
      if (adaptive_substeps .and. .not. have_field_ref) then
        ! The first reference field is the field where the particle starts its track.
        call snapshot%eval_local_e(mesh, pcls_batch%x(:, i), field_ref)
        have_field_ref = .true.
      end if
      chunk_steps = 1_i32
      if (step_growth) then
        ! In a smooth field the particle takes one push of 2^n dt; the chunk never exceeds the flight budget.
        if (h_prev >= app%sim%dt) chunk_steps = max(1_i32, nint(h_prev/app%sim%dt, i32))
        if (grow_ok) chunk_steps = 2_i32*chunk_steps
        chunk_steps = min(chunk_steps, growth_chunk_max, app%sim%max_step - step + 1_i32)
      end if
      remaining = real(chunk_steps, dp)*app%sim%dt
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
        if (step_growth) energy_floor = growth_energy_fraction*0.5_dp*pcls_batch%m(i)*sum(v0*v0)
        halve_loop: do
          call build_particle_step_candidate( &
            mesh, particle_sim, snapshot, bfield, x0, v0, &
            pcls_batch%q(i), pcls_batch%m(i), h, x1, v1, sampled_electric_field &
            )
          if (.not. adaptive_substeps) exit halve_loop
          if (.not. all(ieee_is_finite(x1)) .or. .not. all(ieee_is_finite(v1))) exit halve_loop
          if (step_growth) then
            if (.not. particle_substep_needs_halving(sampled_electric_field, field_ref, pcls_batch%q(i), x0, x1, &
                                                     app%sim%substep_field_change_max, energy_floor)) exit halve_loop
          else
            if (.not. particle_substep_needs_halving(sampled_electric_field, field_ref, pcls_batch%q(i), x0, x1, &
                                                     app%sim%substep_field_change_max)) exit halve_loop
          end if
          if (0.5_dp*h < h_min) then
            depth_limited = .true.
            exit halve_loop
          end if
          h = 0.5_dp*h
          halved = .true.
        end do halve_loop
        rng_step = int(step, i64)
        if (adaptive_substeps) rng_step = int(step, i64)*65536_i64 + int(sub_index, i64)
        boundary_rng_counter = [int(batch_idx, i64), int(rng_tag, i64), workspace%rng_id(i), rng_step]
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
              boundary_rng_counter=boundary_rng_counter, &
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
            boundary_rng_counter=boundary_rng_counter, &
            sampled_electric_field=sampled_electric_field, &
            cell=cell &
            )
          used_event_resolver = .true.
        end if
        if (used_event_resolver) then
          if (step_result%status == particle_step_multiple_box_events .and. retry_upper_panel_fourier) then
            thread_step_counts(1, tid) = thread_step_counts(1, tid) + 1_i64
            call advance_particle_step_upper_panel_fourier( &
              mesh, particle_sim, snapshot, bfield, x0, v0, pcls_batch%q(i), pcls_batch%m(i), h, &
              retry_result, retry_field_available, boundary_contract=particle_boundary_contract, &
              boundary_rng_counter=boundary_rng_counter, &
              cell=cell &
              )
            if (retry_field_available .and. retry_result%status == collision_query_ok) then
              step_result = retry_result
              thread_step_counts(2, tid) = thread_step_counts(2, tid) + 1_i64
            end if
          end if
          if (step_result%status /= collision_query_ok) then
            if (step_result%status == particle_step_multiple_box_events .and. soft_discard_multiple_box_events) then
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
        if (step_growth) then
          if (halved) growth_hold = growth_hold_chunks
          grow_ok = .not. halved .and. .not. depth_limited .and. growth_hold == 0_i32 .and. &
                    .not. particle_substep_needs_halving(sampled_electric_field, field_ref, pcls_batch%q(i), x0, x1, &
                                                         growth_change_max, energy_floor)
          if (h > app%sim%dt) thread_step_counts(6, tid) = thread_step_counts(6, tid) + 1_i64
        end if
        ! Accepted: the sampled midpoint field becomes the reference for the next sub-step.
        field_ref = sampled_electric_field
        h_prev = h
        remaining = remaining - h
        sub_index = sub_index + 1_i32
      end do sub_loop
      if (halved) then
        thread_step_counts(3, tid) = thread_step_counts(3, tid) + 1_i64
        thread_step_counts(4, tid) = thread_step_counts(4, tid) + int(sub_index, i64)
        if (depth_limited) thread_step_counts(5, tid) = thread_step_counts(5, tid) + 1_i64
      end if
      if (particle_done) exit step_loop
      if (.not. halved .and. growth_hold > 0_i32) growth_hold = growth_hold - 1_i32
      step = step + chunk_steps
    end do step_loop
  end subroutine track_particle

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
    counter = [int(batch_idx, i64), int(rng_tag, i64), workspace%rng_id(parent), int(particle_step, i64), &
               int(generation, i64), 0_i64, 271_i64]
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
    integer(i64), allocatable :: child_rng_id(:)
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
    allocate (x(3, added), v(3, added), q(added), m(added), w(added), species(added), elements(added), &
              child_rng_id(added))
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
        ! A child's rng id derives from its parent's, so it does not depend on the local index.
        child_rng_id(index) = counter_child_id(workspace%rng_id(child%parent), int(child%ordinal, i64))
      end do
    end do
    do index = 1, added
      workspace%photo_emission_dq(elements(index), species(index)) = &
        workspace%photo_emission_dq(elements(index), species(index)) - q(index)*w(index)
    end do
    old_count = pcls_batch%n
    call append_particles(pcls_batch, x, v, q, m, w, species, elements)
    call workspace%extend_particle_flags(old_count, pcls_batch%n)
    workspace%rng_id(old_count + 1:pcls_batch%n) = child_rng_id
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
    if (track_failure%status == collision_query_ok .or. particle_index < track_failure%particle .or. &
        (particle_index == track_failure%particle .and. particle_step < track_failure%step)) then
      track_failure = particle_track_failure_type(status, particle_index, particle_step, failure_x, failure_v)
    end if
    !$omp end critical (beach_collision_query_failure)
  end subroutine record_collision_failure

  !> claim の所要時間計測用の壁時計 [s]。OpenMP なしでは `system_clock` に落ちる。
  real(dp) function wall_seconds()
    integer(i64) :: ticks, rate

!$  wall_seconds = omp_get_wtime()
!$  return
    call system_clock(ticks, rate)
    wall_seconds = real(ticks, dp)/real(max(1_i64, rate), dp)
  end function wall_seconds
  end procedure process_particle_batch

end submodule bem_simulator_particles

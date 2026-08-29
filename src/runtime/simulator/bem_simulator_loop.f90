!> `bem_simulator` の主ループと粒子処理計算を実装する submodule。
submodule(bem_simulator) bem_simulator_loop
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use bem_app_config_runtime, only: compute_z_high_box_potential_statistics
  use bem_periodic_zero_mode_plan, only: periodic_zero_mode_state_type
  use bem_performance_profile, only: perf_region_batch_total, perf_region_begin, perf_region_commit_charge, &
                                     perf_region_count_outcomes, perf_region_end, perf_region_field_refresh, &
                                     perf_region_field_solver_init, perf_region_history_write, perf_region_mpi_reduce, &
                                     perf_region_particle_batch, perf_region_prepare_batch, perf_region_simulation_total, &
                                     perf_region_stats_update
  use bem_mpi, only: mpi_allreduce_max_real_dp_array
  use bem_periodic_checkpoint, only: maybe_write_periodic_checkpoint
  use bem_charge_ledger, only: checked_accumulate_charge, finite_charge_sum
  implicit none
  real(dp), parameter :: neutral_return_max_unresolved_fraction = 0.05_dp
  integer(i32), parameter :: adaptive_max_halvings = 24_i32
contains

  module procedure run_absorption_insulator
  integer(i32) :: batch_idx, final_batch_idx, batch_count_this_run, local_batch_idx, nth, hist_stride
  integer(i32) :: team_size_min, team_size_max
  integer(i32) :: particle_team_size, particle_team_size_min, particle_team_size_max
  integer(i32) :: collision_failure_count, collision_failure_rank, collision_failure_status
  integer(i32) :: collision_failure_particle, collision_failure_step
  integer(i32) :: local_failure_values(3), selected_failure_values(3)
  integer(i32) :: photo_failure_count, photo_failure_rank, photo_failure_status, photo_failure_species
  integer(i32) :: photo_failure_ray, photo_failure_bounce
  integer(i32) :: photo_local_failure_values(4), photo_selected_failure_values(4)
  integer(i32) :: trial_halvings, species_idx, fresh_particle_count
  integer(i32) :: matching_iteration, matching_axis, matching_electron_idx, matching_ion_idx, matching_photoelectron_idx
  integer(i32) :: matching_response_status
  integer(i32) :: boundary_status
  integer :: hist_unit, pot_hist_unit, top_ref_hist_unit, matching_hist_unit
  integer, allocatable :: rng_state_before(:)
  logical :: history_enabled, potential_history_enabled, top_reference_history_enabled
  logical :: ledger_enabled, adaptive_nonzero_mode, trial_accepted, omp_dynamic_before
  logical :: matching_active, replay_active, matching_converged, matching_history_enabled
  logical :: implicit_zero_mode, implicit_zero_mode_supported, matching_photoelectron_active
  real(dp), allocatable :: potential_buf(:), injection_residual_before(:), boundary_injection_residual_before(:, :)
  integer(i64) :: batch_counts(6), batch_retry_counts(2)
  real(dp) :: bfield(3), rel, t0, sim_t0, batch_t0, batch_soft_discarded_abs_charge
  real(dp) :: collision_failure_x(3), collision_failure_v(3), selected_failure_state(6)
  real(dp) :: trial_batch_duration, duration_ratio, adaptive_potential_step, adaptive_metric_values(1)
  real(dp) :: projected_simulated_time
  real(dp) :: top_phi_mean, top_phi_std, top_phi_min, top_phi_max
  real(dp) :: matching_plane_z, matching_area, matching_displacement, matching_displacement_before, matching_residual
  real(dp) :: matching_committed_displacement, matching_displacement_tolerance, matching_roundoff_scale
  real(dp) :: matching_return_flux, matching_escape_flux, matching_budget_scale
  real(dp) :: matching_photoelectron_charge, matching_photoelectron_emission_current_density
  real(dp) :: matching_guess(4), matching_observed(4)
  real(dp) :: implicit_displacement_min, implicit_displacement_max, implicit_feedback_reference(4)
  real(dp) :: matching_feedback_scales(4)
  real(dp) :: matching_component_residuals(4), matching_absolute_defects(4)
  real(dp) :: matching_response_input(5), matching_response_output(6)
  real(dp), allocatable :: matching_moments(:, :), matching_reduce(:)
  character(len=256) :: boundary_message
  character(len=512) :: matching_response_message
  type(particles_soa) :: pcls_batch
  type(mpi_context) :: mpi_ctx
  type(electrostatic_snapshot_type) :: snapshot
  type(electrostatic_diagnostics_type) :: committed_snapshot_diagnostics
  type(periodic_zero_mode_state_type) :: committed_zero_state
  type(charge_ledger_type) :: batch_ledger
  type(simulator_batch_workspace_type) :: workspace
  type(particle_source_plan_type) :: source_plan
  type(external_boundary_contract_type) :: boundary_contract
  type(surface_closure_contract_type) :: surface_closure
  type(matching_plane_response_provider_type) :: matching_response_provider
  type(field_physics_config) :: field_config
  type(panel_kernel_config) :: panel_config
  type(app_config) :: trial_app
  type(sim_stats) :: stats_candidate

  stats = sim_stats()
  if (present(initial_stats)) stats = initial_stats
  mpi_ctx = mpi_context()
  if (present(mpi)) mpi_ctx = mpi
  ledger_enabled = present(charge_ledger)
  if (ledger_enabled) then
    if (.not. allocated(charge_ledger%injected_from_remote)) then
      call charge_ledger%init(app%n_particle_species)
    else if (charge_ledger%nspecies /= app%n_particle_species) then
      error stop 'charge ledger species count does not match app config.'
    end if
    call batch_ledger%init(app%n_particle_species)
  end if

  call resolve_external_boundary_contract( &
    app%sim%reservoir_potential_model, app%sim%open_boundary_model, &
    boundary_contract, boundary_status, boundary_message &
    )
  if (boundary_status /= external_boundary_ok) error stop trim(boundary_message)
  adaptive_nonzero_mode = app%periodic2%max_nonzero_mode_potential_step > 0.0_dp
  matching_active = trim(lower_ascii(app%surface_current%model)) == 'matching_plane_quasistatic'
  implicit_zero_mode = matching_active .and. app%surface_current%implicit_zero_mode
  implicit_zero_mode_supported = .false.
  implicit_displacement_min = 0.0_dp
  implicit_displacement_max = 0.0_dp
  implicit_feedback_reference = 0.0_dp
  call matching_response_provider%initialize( &
    app, mpi_ctx, matching_response_status, matching_response_message &
    )
  if (matching_response_status /= matching_plane_provider_ok) then
    error stop 'matching-plane response preflight failed: '//trim(matching_response_message)
  end if
  replay_active = adaptive_nonzero_mode .or. matching_active
  omp_dynamic_before = .false.
!$ omp_dynamic_before = omp_get_dynamic()
  nth = 1_i32
  if (replay_active) then
!$  call omp_set_dynamic(.false.)
    !$omp parallel default(none) shared(nth)
    !$omp single
!$  nth = int(omp_get_num_threads(), i32)
    !$omp end single
    !$omp end parallel

    team_size_min = nth
    team_size_max = nth
    call mpi_allreduce_min_i32_scalar(mpi_ctx, team_size_min)
    call mpi_allreduce_max_i32_scalar(mpi_ctx, team_size_max)
    if (team_size_min /= team_size_max) then
      if (mpi_is_root(mpi_ctx)) then
        write (error_unit, '(a,i0,a,i0)') &
          'replayed-trial OpenMP team-size mismatch across MPI ranks: min=', team_size_min, ' max=', team_size_max
        flush (error_unit)
      end if
      error stop 'replayed-trial OpenMP team size must match across MPI ranks.'
    end if

    ! Rejected trials within this run still use one fixed team size.  A restart may
    ! use a different team size because checkpoint compatibility is numerical, not
    ! bitwise-replay compatibility with the previous process.
    if (adaptive_nonzero_mode) stats%adaptive_nonzero_mode_omp_threads = nth
  else
!$  nth = max(1_i32, int(omp_get_max_threads(), i32))
  end if
  if (.not. adaptive_nonzero_mode) then
    stats%adaptive_nonzero_mode_rejected_trials = 0_i64
    stats%adaptive_nonzero_mode_last_batch_duration = 0.0_dp
    stats%adaptive_nonzero_mode_last_potential_step = 0.0_dp
    stats%adaptive_nonzero_mode_omp_threads = 0_i32
  end if
  call enforce_soft_discard_limits(app%sim, stats, stats%batches, 0_i64, 0.0_dp, mpi_ctx)
  call workspace%init(mesh%nelem, app%n_particle_species, nth, candidate_charge_enabled=adaptive_nonzero_mode)
  call evaluate_surface_closure(app, surface_closure)

  matching_history_enabled = present(matching_plane_history_unit)
  matching_hist_unit = -1
  if (matching_history_enabled) then
    matching_hist_unit = matching_plane_history_unit
    matching_history_enabled = matching_hist_unit /= -1
  end if

  history_enabled = present(history_unit)
  hist_unit = 0
  if (history_enabled) hist_unit = history_unit
  hist_stride = 1_i32
  if (present(history_stride)) then
    matching_history_enabled = matching_history_enabled .and. history_stride > 0_i32
    hist_stride = max(1_i32, history_stride)
  end if
  potential_history_enabled = present(potential_history_unit)
  pot_hist_unit = -1
  if (potential_history_enabled) then
    pot_hist_unit = potential_history_unit
    allocate (potential_buf(mesh%nelem))
  end if
  top_reference_history_enabled = present(top_reference_history_unit)
  top_ref_hist_unit = -1
  if (top_reference_history_enabled) then
    top_ref_hist_unit = top_reference_history_unit
    top_reference_history_enabled = top_ref_hist_unit /= -1
  end if

  bfield = app%sim%b0
  final_batch_idx = app%sim%batch_count
  if (stats%batches < 0_i32 .or. stats%batches > final_batch_idx) then
    error stop 'checkpoint batch count is outside sim.batch_count.'
  end if
  if (.not. adaptive_nonzero_mode .and. stats%batches > 0_i32 .and. stats%simulated_time == 0.0_dp) then
    stats%simulated_time = real(stats%batches, dp)*app%sim%batch_duration
  end if
  batch_count_this_run = final_batch_idx - stats%batches

  call perf_region_begin(perf_region_simulation_total, sim_t0)
  matching_electron_idx = 0_i32
  matching_ion_idx = 0_i32
  matching_photoelectron_idx = 0_i32
  matching_photoelectron_active = .false.
  matching_photoelectron_charge = 0.0_dp
  matching_photoelectron_emission_current_density = 0.0_dp
  if (matching_active) then
    call matching_response_provider%get_matching_plane_z( &
      matching_plane_z, matching_response_status, matching_response_message &
      )
    if (matching_response_status /= matching_plane_provider_ok) then
      error stop 'matching-plane response metadata failed: '//trim(matching_response_message)
    end if
    matching_electron_idx = matching_species_index(app, app%surface_current%electron_species)
    matching_ion_idx = matching_species_index(app, app%surface_current%ion_species)
    matching_photoelectron_active = len_trim(app%surface_current%photoelectron_species) > 0
    if (matching_photoelectron_active) then
      matching_photoelectron_idx = matching_species_index(app, app%surface_current%photoelectron_species)
      matching_photoelectron_charge = app%particle_species(matching_photoelectron_idx)%q_particle
      matching_photoelectron_emission_current_density = &
        app%particle_species(matching_photoelectron_idx)%emit_current_density_a_m2
    end if
    matching_area = product(app%sim%box_max(1:2) - app%sim%box_min(1:2))
    if (implicit_zero_mode) then
      call matching_response_provider%get_implicit_zero_mode_contract( &
        implicit_zero_mode_supported, implicit_displacement_min, implicit_displacement_max, &
        implicit_feedback_reference &
        )
      if (.not. implicit_zero_mode_supported) then
        error stop 'implicit matching-plane zero mode requires a table with a non-singleton displacement axis and '// &
          'singleton PE-flux/energy and ambient-outward axes; the PE pair may both be zero.'
      end if
      if (matching_photoelectron_active) then
        if (any(implicit_feedback_reference(1:2) <= 0.0_dp)) then
          error stop 'implicit matching-plane PE mode requires positive singleton PE flux and energy.'
        end if
      else
        if (any(implicit_feedback_reference(1:2) /= 0.0_dp)) then
          error stop 'implicit matching-plane no-PE mode requires zero singleton PE flux and energy.'
        end if
      end if
    end if
    if (mesh%nelem > 0_i32) then
      if (app%sim%box_max(3) <= max( &
          maxval(mesh%v0(3, :)), maxval(mesh%v1(3, :)), maxval(mesh%v2(3, :)) &
          )) then
        error stop 'matching-plane gauge must lie strictly above every mesh vertex.'
      end if
    end if
    allocate (matching_moments(4_i32, app%n_particle_species))
    allocate (matching_reduce(4_i32*app%n_particle_species))
    if (stats%batches > 0_i32 .and. .not. stats%matching_plane_state_valid) then
      error stop 'matching-plane resume requires a schema-v9 committed coupling state.'
    end if
  else
    stats%matching_plane_state_valid = .false.
  end if

  call perf_region_begin(perf_region_field_solver_init, t0)
  call derive_field_panel_config(app%sim, field_config, panel_config)
  call snapshot%init(mesh, app%sim, field_config, app%periodic2, panel_config)
  call perf_region_end(perf_region_field_solver_init, t0)
  if (matching_active .and. stats%matching_plane_state_valid) then
    call snapshot%set_matching_plane_gauge(mesh, matching_plane_z, stats%matching_plane_phi_v)
  end if

  if (replay_active) then
    call random_seed(size=species_idx)
    allocate (rng_state_before(species_idx))
    if (present(inject_state)) then
      if (allocated(inject_state%macro_residual)) allocate (injection_residual_before(size(inject_state%macro_residual)))
      if (allocated(inject_state%boundary_macro_residual)) then
        allocate ( &
          boundary_injection_residual_before( &
          size(inject_state%boundary_macro_residual, 1), size(inject_state%boundary_macro_residual, 2) &
          ) &
          )
      end if
    end if
  end if
  adaptive_potential_step = 0.0_dp

  do local_batch_idx = 1_i32, batch_count_this_run
    call perf_region_begin(perf_region_batch_total, batch_t0)
    batch_idx = stats%batches + 1_i32
    call perf_region_begin(perf_region_field_refresh, t0)
    call snapshot%refresh(mesh)
    call perf_region_end(perf_region_field_refresh, t0)

    trial_batch_duration = app%sim%batch_duration
    trial_halvings = 0_i32
    trial_accepted = .false.
    if (replay_active) then
      call random_seed(get=rng_state_before)
      if (allocated(injection_residual_before)) injection_residual_before = inject_state%macro_residual
      if (allocated(boundary_injection_residual_before)) then
        boundary_injection_residual_before = inject_state%boundary_macro_residual
      end if
      committed_zero_state = snapshot%zero_state
      committed_snapshot_diagnostics = snapshot%diagnostics
    end if

    do while (.not. trial_accepted)
      if (replay_active) then
        call random_seed(put=rng_state_before)
        if (allocated(injection_residual_before)) inject_state%macro_residual = injection_residual_before
        if (allocated(boundary_injection_residual_before)) then
          inject_state%boundary_macro_residual = boundary_injection_residual_before
        end if
        snapshot%zero_state = committed_zero_state
        snapshot%diagnostics = committed_snapshot_diagnostics
      end if
      trial_app = app
      trial_app%sim%batch_duration = trial_batch_duration
      duration_ratio = trial_batch_duration/app%sim%batch_duration
      do species_idx = 1_i32, trial_app%n_particle_species
        if (.not. trial_app%particle_species(species_idx)%enabled) cycle
        if (.not. trial_app%particle_species(species_idx)%has_target_macro_particles_per_batch) cycle
        if (trim(lower_ascii(trial_app%particle_species(species_idx)%source_mode)) /= 'reservoir_face' .and. &
            trim(lower_ascii(trial_app%particle_species(species_idx)%source_mode)) /= 'plane_source' .and. &
            .not. any(trial_app%particle_species(species_idx)%boundary_inflow_low /= 0_i32) .and. &
            .not. any(trial_app%particle_species(species_idx)%boundary_inflow_high /= 0_i32)) cycle
        trial_app%particle_species(species_idx)%w_particle = &
          app%particle_species(species_idx)%w_particle*duration_ratio
      end do

      matching_iteration = 0_i32
      matching_converged = .not. matching_active
      matching_residual = 0.0_dp
      matching_return_flux = 0.0_dp
      matching_escape_flux = 0.0_dp
      if (matching_active) then
        matching_displacement_before = snapshot%get_matching_plane_displacement()
        matching_displacement = matching_displacement_before
        if (implicit_zero_mode) then
          call solve_matching_implicit_zero_mode( &
            matching_response_provider, mpi_ctx, matching_displacement_before, trial_batch_duration, &
            implicit_displacement_min, implicit_displacement_max, implicit_feedback_reference, &
            app%particle_species(matching_electron_idx)%q_particle, &
            app%particle_species(matching_ion_idx)%q_particle, &
            matching_photoelectron_active, matching_photoelectron_charge, &
            matching_displacement, matching_response_output &
            )
          matching_guess = implicit_feedback_reference
        else if (stats%matching_plane_state_valid) then
          matching_guess = stats%matching_plane_feedback
        else
          matching_guess = 0.0_dp
        end if
      end if

      do
        matching_iteration = matching_iteration + 1_i32
        if (matching_active) then
          call random_seed(put=rng_state_before)
          if (allocated(injection_residual_before)) inject_state%macro_residual = injection_residual_before
          if (allocated(boundary_injection_residual_before)) then
            inject_state%boundary_macro_residual = boundary_injection_residual_before
          end if
          matching_response_input = [matching_displacement, matching_guess]
          call matching_response_provider%evaluate( &
            matching_response_input, mpi_ctx, matching_response_output, &
            matching_response_status, matching_response_message &
            )
          if (matching_response_status /= matching_plane_provider_ok) then
            error stop 'matching-plane response evaluation failed: '//trim(matching_response_message)
          end if
          call configure_matching_surface_closure( &
            surface_closure, matching_electron_idx, matching_ion_idx, matching_photoelectron_idx, &
            matching_photoelectron_active, matching_response_output, implicit_zero_mode, matching_area, &
            implicit_feedback_reference, &
            app%particle_species(matching_electron_idx)%q_particle, &
            app%particle_species(matching_ion_idx)%q_particle, &
            matching_photoelectron_charge, matching_photoelectron_emission_current_density &
            )
          call snapshot%set_matching_plane_gauge(mesh, matching_plane_z, matching_response_output(1))
        end if

        call perf_region_begin(perf_region_prepare_batch, t0)
        call build_particle_source_plan( &
          trial_app, source_plan, mpi=mpi_ctx, &
          kinetic_inflow_active=surface_closure%has_inflow_kinetic_map, &
          kinetic_reservoir_potential_v=surface_closure%inflow_reservoir_potential_v, &
          kinetic_access_potential_v=surface_closure%inflow_access_potential_v, &
          kinetic_inflow_face=surface_closure%inflow_kinetic_face, &
          number_flux_override_active=surface_closure%has_inflow_number_flux, &
          number_flux_override_m2_s=surface_closure%inflow_number_flux_m2_s &
          )
        call prepare_batch_state( &
          mesh, trial_app, source_plan, snapshot, stats, batch_idx, workspace, pcls_batch, mpi_ctx, inject_state, &
          photo_failure_status, photo_failure_species, photo_failure_ray, photo_failure_bounce &
          )
        call perf_region_end(perf_region_prepare_batch, t0)
        fresh_particle_count = pcls_batch%n

        photo_failure_count = merge(1_i32, 0_i32, photo_failure_status /= collision_query_ok)
        call mpi_allreduce_sum_i32_scalar(mpi_ctx, photo_failure_count)
        if (photo_failure_count > 0_i32) then
          photo_local_failure_values = [photo_failure_species, photo_failure_ray, photo_failure_bounce, photo_failure_status]
          call mpi_select_lowest_rank_i32_values( &
            mpi_ctx, photo_failure_status /= collision_query_ok, photo_local_failure_values, &
            photo_failure_rank, photo_selected_failure_values &
            )
          call stop_for_photo_collision_failure( &
            batch_idx, photo_failure_rank, photo_selected_failure_values(1), photo_selected_failure_values(2), &
            photo_selected_failure_values(3), photo_selected_failure_values(4) &
            )
        end if

        call perf_region_begin(perf_region_particle_batch, t0)
        call process_particle_batch( &
          mesh, trial_app, boundary_contract, surface_closure, snapshot, pcls_batch, workspace%dq_thread, &
          workspace%escaped_boundary_flag, workspace%absorbed_flag, workspace%absorbed_element, &
          workspace%soft_discarded_boundary_flag, bfield, batch_idx, mpi_ctx%rank, particle_team_size, &
          collision_failure_status, collision_failure_particle, collision_failure_step, &
          collision_failure_x, collision_failure_v, workspace%matching_plane_moments_thread, &
          batch_retry_counts &
          )
        call perf_region_end(perf_region_particle_batch, t0)

        if (replay_active) then
          particle_team_size_min = particle_team_size
          particle_team_size_max = particle_team_size
          call mpi_allreduce_min_i32_scalar(mpi_ctx, particle_team_size_min)
          call mpi_allreduce_max_i32_scalar(mpi_ctx, particle_team_size_max)
          if (particle_team_size_min /= nth .or. particle_team_size_max /= nth) then
            if (mpi_is_root(mpi_ctx)) then
              write (error_unit, '(a,i0,a,i0,a,i0)') &
                'replayed-trial OpenMP team size changed after replay probe: expected=', nth, &
                ' min=', particle_team_size_min, ' max=', particle_team_size_max
              flush (error_unit)
            end if
            error stop 'replayed-trial OpenMP team size changed after the replay probe.'
          end if
        end if

        collision_failure_count = merge(1_i32, 0_i32, collision_failure_status /= collision_query_ok)
        call mpi_allreduce_sum_i32_scalar(mpi_ctx, collision_failure_count)
        if (collision_failure_count > 0_i32) then
          local_failure_values = [collision_failure_status, collision_failure_particle, collision_failure_step]
          call mpi_select_lowest_rank_i32_values( &
            mpi_ctx, collision_failure_status /= collision_query_ok, local_failure_values, &
            collision_failure_rank, selected_failure_values &
            )
          selected_failure_state = 0.0_dp
          if (mpi_ctx%rank == collision_failure_rank) then
            selected_failure_state(1:3) = collision_failure_x
            selected_failure_state(4:6) = collision_failure_v
          end if
          call mpi_allreduce_sum_real_dp_array(mpi_ctx, selected_failure_state)
          call stop_for_collision_failure( &
            batch_idx, collision_failure_rank, selected_failure_values(1), selected_failure_values(2), &
            selected_failure_values(3), trial_app%sim%dt, selected_failure_state(1:3), selected_failure_state(4:6) &
            )
        end if

        call apply_neutral_return_surface_closure(trial_app, pcls_batch, fresh_particle_count, workspace, mpi_ctx)
        call apply_fixed_surface_current_closure( &
          trial_app, surface_closure, pcls_batch, fresh_particle_count, workspace, mpi_ctx &
          )
        if (.not. matching_active) exit

        matching_moments = sum(workspace%matching_plane_moments_thread, dim=3)
        matching_reduce = reshape(matching_moments, [size(matching_reduce)])
        call mpi_allreduce_sum_real_dp_array(mpi_ctx, matching_reduce)
        matching_moments = reshape(matching_reduce, shape(matching_moments))
        call resolve_matching_observed_feedback( &
          matching_moments, matching_electron_idx, matching_ion_idx, matching_photoelectron_idx, &
          matching_photoelectron_active, matching_area, trial_batch_duration, matching_observed, &
          matching_return_flux, matching_escape_flux &
          )
        call matching_response_provider%validate_feedback( &
          matching_observed, matching_response_status, matching_response_message &
          )
        if (matching_response_status /= matching_plane_provider_ok) then
          error stop 'matching-plane feedback validation failed: '//trim(matching_response_message)
        end if
        if (matching_photoelectron_active) then
          matching_budget_scale = max(1.0_dp, matching_observed(1), matching_return_flux, matching_escape_flux)
          if (abs(matching_observed(1) - matching_return_flux - matching_escape_flux) > &
              sqrt(epsilon(1.0_dp))*matching_budget_scale) then
            error stop 'matching-plane photoelectron outflow budget does not close.'
          end if
        end if
        matching_residual = matching_response_provider%feedback_residual( &
                            matching_guess, matching_observed, app%surface_current%coupling_rtol, &
                            app%surface_current%coupling_atol &
                            )
        matching_converged = matching_response_provider%feedback_converged( &
                             matching_guess, matching_observed, app%surface_current%coupling_rtol, &
                             app%surface_current%coupling_atol &
                             )
        if (matching_converged) exit
        if (matching_iteration >= app%surface_current%coupling_max_iterations) then
          if (mpi_is_root(mpi_ctx)) then
            call matching_response_provider%get_feedback_scales(matching_feedback_scales)
            matching_absolute_defects = abs(matching_observed - matching_guess)
            matching_component_residuals = 0.0_dp
            do matching_axis = 1_i32, 4_i32
              if (matching_feedback_scales(matching_axis) <= 0.0_dp) cycle
              if (app%surface_current%coupling_atol(matching_axis) > &
                  app%surface_current%coupling_rtol*matching_feedback_scales(matching_axis)) then
                matching_component_residuals(matching_axis) = app%surface_current%coupling_rtol* &
                                                              (matching_absolute_defects(matching_axis)/ &
                                                               app%surface_current%coupling_atol(matching_axis))
              else
                matching_component_residuals(matching_axis) = matching_absolute_defects(matching_axis)/ &
                                                              matching_feedback_scales(matching_axis)
              end if
            end do
            write (error_unit, '(a,i0,a,i0,a,es24.16,a,es24.16)') &
              'matching-plane nonconvergence: batch=', batch_idx, ', iterations=', matching_iteration, &
              ', residual=', matching_residual, ', rtol=', app%surface_current%coupling_rtol
            write (error_unit, '(a,4(1x,es24.16))') 'matching-plane guess=', matching_guess
            write (error_unit, '(a,4(1x,es24.16))') 'matching-plane observed=', matching_observed
            write (error_unit, '(a,4(1x,es24.16))') 'matching-plane feedback scales=', matching_feedback_scales
            write (error_unit, '(a,4(1x,es24.16))') &
              'matching-plane coupling atols=', app%surface_current%coupling_atol
            write (error_unit, '(a,4(1x,es24.16))') &
              'matching-plane absolute defects=', matching_absolute_defects
            write (error_unit, '(a,4(1x,es24.16))') &
              'matching-plane effective component residuals=', matching_component_residuals
            flush (error_unit)
          end if
          error stop 'matching-plane fixed point did not converge before coupling_max_iterations.'
        end if
        matching_guess = matching_guess + app%surface_current%coupling_relaxation* &
                         (matching_observed - matching_guess)
      end do

      if (adaptive_nonzero_mode) then
        call prepare_adaptive_charge_candidate(mesh, workspace, mpi_ctx)
        adaptive_metric_values = 0.0_dp
        if (mpi_is_root(mpi_ctx)) then
          call snapshot%measure_kneq0_potential_step( &
            mesh, workspace%candidate_charge, adaptive_metric_values(1) &
            )
        end if
        call mpi_allreduce_max_real_dp_array(mpi_ctx, adaptive_metric_values)
        adaptive_potential_step = adaptive_metric_values(1)
        trial_accepted = adaptive_potential_step <= app%periodic2%max_nonzero_mode_potential_step
      else
        trial_accepted = .true.
      end if
      if (trial_accepted) exit

      workspace%charge_candidate_ready = .false.
      if (trial_halvings >= adaptive_max_halvings) then
        error stop 'adaptive nonzero-mode batch failed after 24 duration halvings.'
      end if
      trial_halvings = trial_halvings + 1_i32
      trial_batch_duration = scale(app%sim%batch_duration, -trial_halvings)
      if (.not. ieee_is_finite(trial_batch_duration) .or. trial_batch_duration <= 0.0_dp) then
        error stop 'adaptive nonzero-mode batch duration became invalid.'
      end if
    end do

    call perf_region_begin(perf_region_count_outcomes, t0)
    call count_batch_outcomes( &
      pcls_batch, workspace%escaped_boundary_flag, workspace%absorbed_flag, &
      workspace%soft_discarded_boundary_flag, batch_counts, batch_soft_discarded_abs_charge &
      )
    call perf_region_end(perf_region_count_outcomes, t0)
    call perf_region_begin(perf_region_mpi_reduce, t0)
    call mpi_allreduce_sum_i64_array(mpi_ctx, batch_counts)
    call mpi_allreduce_sum_i64_array(mpi_ctx, batch_retry_counts)
    call mpi_allreduce_sum_real_dp_scalar(mpi_ctx, batch_soft_discarded_abs_charge)
    call perf_region_end(perf_region_mpi_reduce, t0)
    call report_soft_discard_summary(batch_idx, batch_counts(6), batch_soft_discarded_abs_charge, mpi_ctx)
    call enforce_soft_discard_limits( &
      app%sim, stats, batch_idx, batch_counts(6), batch_soft_discarded_abs_charge, mpi_ctx &
      )
    projected_simulated_time = stats%simulated_time + trial_batch_duration
    if (.not. ieee_is_finite(stats%simulated_time) .or. stats%simulated_time < 0.0_dp .or. &
        .not. ieee_is_finite(projected_simulated_time) .or. projected_simulated_time < stats%simulated_time) then
      error stop 'simulation statistic overflow: simulated_time'
    end if
    ! Build the complete statistics update before mutating mesh/ledger state.  The
    ! checked add routines are the single overflow guard; commit only after all pass.
    stats_candidate = stats
    call accumulate_batch_stats( &
      stats_candidate, batch_counts, batch_soft_discarded_abs_charge, batch_retry_counts, 0.0_dp &
      )
    stats_candidate%simulated_time = projected_simulated_time
    if (adaptive_nonzero_mode) then
      stats_candidate%adaptive_nonzero_mode_last_batch_duration = trial_batch_duration
      stats_candidate%adaptive_nonzero_mode_last_potential_step = adaptive_potential_step
      stats_candidate%adaptive_nonzero_mode_rejected_trials = checked_add_adaptive_rejected_trials( &
                                                              stats%adaptive_nonzero_mode_rejected_trials, &
                                                              trial_halvings &
                                                              )
    end if
    if (matching_active) then
      stats_candidate%matching_plane_state_valid = .true.
      stats_candidate%matching_plane_displacement_c_m2 = matching_displacement
      stats_candidate%matching_plane_phi_v = matching_response_output(1)
      stats_candidate%matching_plane_response = matching_response_output
      stats_candidate%matching_plane_feedback = matching_observed
      stats_candidate%matching_plane_photoelectron_return_flux_m2_s = matching_return_flux
      stats_candidate%matching_plane_photoelectron_escape_flux_m2_s = matching_escape_flux
      stats_candidate%matching_plane_iterations = matching_iteration
      stats_candidate%matching_plane_residual = matching_residual
    end if

    if (ledger_enabled) then
      call batch_ledger%reset(batch_idx)
      batch_ledger%surface_charge_before = finite_charge_sum(mesh%q_elem, 'batch surface charge before commit')
      call record_batch_initial_charge(trial_app, pcls_batch, fresh_particle_count, batch_ledger)
      call record_batch_outcome_charge( &
        pcls_batch, workspace%escaped_boundary_flag, workspace%absorbed_flag, &
        workspace%soft_discarded_boundary_flag, batch_ledger &
        )
      call reduce_charge_ledger_fluxes(batch_ledger, mpi_ctx, workspace)
      batch_ledger%neutral_return_correction = workspace%neutral_return_correction
      batch_ledger%neutral_return_weight_scale = workspace%neutral_return_weight_scale
      batch_ledger%neutral_return_unresolved_fraction = workspace%neutral_return_unresolved_fraction
      batch_ledger%fixed_absorbed_target_charge = workspace%fixed_absorbed_target_charge
      batch_ledger%fixed_absorbed_weight_scale = workspace%fixed_absorbed_weight_scale
      batch_ledger%fixed_emission_target_charge = workspace%fixed_emission_target_charge
      batch_ledger%fixed_emission_weight_scale = workspace%fixed_emission_weight_scale
      batch_ledger%fixed_escape_target_charge = workspace%fixed_escape_target_charge
      batch_ledger%fixed_escape_correction = workspace%fixed_escape_correction
      batch_ledger%fixed_current_correction = workspace%fixed_current_correction
    end if

    call perf_region_begin(perf_region_commit_charge, t0)
    call commit_batch_charge( &
      mesh, app%sim%q_floor, app%sim%e0, app%sim%field_bc_mode, workspace, rel, mpi_ctx &
      )
    call perf_region_end(perf_region_commit_charge, t0)
    if (implicit_zero_mode) then
      matching_committed_displacement = finite_charge_sum( &
                                        mesh%q_elem, 'implicit matching-plane committed charge' &
                                        )/matching_area
      matching_roundoff_scale = finite_charge_sum( &
                                [ &
                                abs(workspace%fixed_absorbed_target_charge), &
                                abs(workspace%fixed_emission_target_charge) &
                                ], &
                                'implicit matching-plane gross target charge' &
                                )/matching_area
      matching_displacement_tolerance = 4096.0_dp*epsilon(1.0_dp)*max( &
                                        abs(matching_displacement), &
                                        abs(matching_committed_displacement), &
                                        matching_roundoff_scale, tiny(1.0_dp) &
                                        )
      if (abs(matching_committed_displacement - matching_displacement) > matching_displacement_tolerance) then
        error stop 'implicit matching-plane committed displacement does not match the backward-Euler endpoint.'
      end if
      stats_candidate%matching_plane_displacement_c_m2 = matching_committed_displacement
    end if
    stats_candidate%last_rel_change = rel
    if (ledger_enabled) then
      batch_ledger%surface_charge_after = finite_charge_sum(mesh%q_elem, 'batch surface charge after commit')
      call accumulate_charge_ledger(charge_ledger, batch_ledger)
    end if

    call perf_region_begin(perf_region_stats_update, t0)
    stats = stats_candidate
    call perf_region_end(perf_region_stats_update, t0)

    call perf_region_begin(perf_region_history_write, t0)
    if (mpi_is_root(mpi_ctx)) then
      call print_batch_progress(batch_idx, final_batch_idx, rel)
      call maybe_write_history_snapshot(history_enabled, hist_unit, hist_stride, stats, rel, mesh%q_elem)
      if (matching_history_enabled .and. mod(batch_idx - 1_i32, hist_stride) == 0_i32) then
        call write_matching_plane_history_snapshot(matching_hist_unit, batch_idx, stats%simulated_time, stats)
      end if
      if (potential_history_enabled) then
        call maybe_write_potential_history_snapshot( &
          potential_history_enabled, pot_hist_unit, hist_stride, stats, snapshot, mesh, app%sim, potential_buf, &
          top_reference_history_enabled, top_ref_hist_unit &
          )
      end if
    end if
    call perf_region_end(perf_region_history_write, t0)
    if (present(inject_state)) then
      if (ledger_enabled) then
        call maybe_write_periodic_checkpoint(app, mesh, stats, inject_state, mpi_ctx, charge_ledger)
      else
        call maybe_write_periodic_checkpoint(app, mesh, stats, inject_state, mpi_ctx)
      end if
    else if (ledger_enabled) then
      call maybe_write_periodic_checkpoint(app, mesh, stats, mpi=mpi_ctx, charge_ledger=charge_ledger)
    else
      call maybe_write_periodic_checkpoint(app, mesh, stats, mpi=mpi_ctx)
    end if
    call perf_region_end(perf_region_batch_total, batch_t0)
  end do
  call perf_region_end(perf_region_simulation_total, sim_t0)
!$ if (replay_active) call omp_set_dynamic(omp_dynamic_before)

  if (present(mesh_potential_v) .and. mpi_is_root(mpi_ctx)) then
    call snapshot%refresh(mesh)
    allocate (mesh_potential_v(mesh%nelem))
    call snapshot%compute_mesh_potential(mesh, app%sim, mesh_potential_v)
  end if
  if (present(electrostatic_diagnostics)) then
    call snapshot%get_diagnostics(electrostatic_diagnostics)
    if (mpi_is_root(mpi_ctx) .and. app%sim%use_box) then
      call snapshot%refresh(mesh)
      call compute_z_high_box_potential_statistics( &
        mesh, app%sim, snapshot, top_phi_mean, top_phi_std, top_phi_min, top_phi_max &
        )
      electrostatic_diagnostics%top_reference_available = .true.
      electrostatic_diagnostics%top_reference_last_batch = stats%batches
      electrostatic_diagnostics%top_reference_simulated_time = stats%simulated_time
      electrostatic_diagnostics%top_reference_z_high = app%sim%box_max(3)
      electrostatic_diagnostics%top_reference_sample_n = app%sim%injection_face_phi_grid_n
      electrostatic_diagnostics%top_reference_potential_mean = top_phi_mean
      electrostatic_diagnostics%top_reference_potential_std = top_phi_std
      electrostatic_diagnostics%top_reference_potential_min = top_phi_min
      electrostatic_diagnostics%top_reference_potential_max = top_phi_max
    end if
  end if
  end procedure run_absorption_insulator

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
        if (size(matching_plane_moments_thread, 1) < 4 .or. &
            species_idx > size(matching_plane_moments_thread, 2) .or. &
            tid > size(matching_plane_moments_thread, 3)) then
          call record_collision_failure(particle_step_invalid_boundary, i, step, x0, v0)
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

  function checked_add_adaptive_rejected_trials(accumulated, trial_halvings) result(total)
    integer(i64), intent(in) :: accumulated
    integer(i32), intent(in) :: trial_halvings
    integer(i64) :: increment, total

    increment = int(trial_halvings, i64)
    if (accumulated < 0_i64 .or. increment < 0_i64) then
      error stop 'adaptive rejected-trial count must be nonnegative.'
    end if
    if (accumulated > huge(accumulated) - increment) then
      error stop 'adaptive rejected-trial count overflow.'
    end if
    total = accumulated + increment
  end function checked_add_adaptive_rejected_trials

  subroutine stop_for_collision_failure(batch_idx, rank, status, particle, step, dt, x, v)
    integer(i32), intent(in) :: batch_idx, rank, status, particle, step
    real(dp), intent(in) :: dt, x(3), v(3)
    write (error_unit, '(a,i0,a,i0,a,i0,a,i0,a,a,a,i0,a,es13.5,a,3es13.5,a,3es13.5)') &
      'particle step failed: batch=', batch_idx, ' particle=', particle, ' step=', step, &
      ' rank=', rank, ' status=', trim(particle_failure_name(status)), ' code=', status, &
      ' dt=', dt, ' x=', x, ' v=', v
    flush (error_unit)
    error stop 'particle step failed.'
  end subroutine stop_for_collision_failure

  subroutine stop_for_photo_collision_failure(batch_idx, rank, species, ray, bounce, status)
    integer(i32), intent(in) :: batch_idx, rank, species, ray, bounce, status
    write (error_unit, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,a,a,i0)') &
      'photo ray collision failed: batch=', batch_idx, ' rank=', rank, ' species=', species, &
      ' ray=', ray, ' bounce=', bounce, ' status=', trim(collision_failure_name(status)), ' code=', status
    flush (error_unit)
    error stop 'photo ray collision query failed.'
  end subroutine stop_for_photo_collision_failure

  !> accepted trial の soft-discard 集約だけを root へ出力する。
  subroutine report_soft_discard_summary(batch_idx, global_count, global_abs_charge, mpi)
    integer(i32), intent(in) :: batch_idx
    integer(i64), intent(in) :: global_count
    real(dp), intent(in) :: global_abs_charge
    type(mpi_context), intent(in) :: mpi

    if (global_count < 0_i64) error stop 'soft-discard global count overflowed.'
    if (global_count == 0_i64) return

    if (mpi_is_root(mpi)) then
      write (error_unit, '(a,i0,a,i0,a,es13.5)') &
        'multiple_box_events soft discard accepted: batch=', batch_idx, &
        ' global_count=', global_count, ' global_abs_charge_C=', global_abs_charge
      flush (error_unit)
    end if
  end subroutine report_soft_discard_summary

  subroutine enforce_soft_discard_limits(sim, stats, batch_idx, batch_count, batch_abs_charge, mpi)
    type(sim_config), intent(in) :: sim
    type(sim_stats), intent(in) :: stats
    integer(i32), intent(in) :: batch_idx
    integer(i64), intent(in) :: batch_count
    real(dp), intent(in) :: batch_abs_charge
    type(mpi_context), intent(in) :: mpi
    integer(i64) :: projected_count, count_limit
    real(dp) :: projected_abs_charge
    logical :: count_exceeded, charge_exceeded

    if (trim(lower_ascii(sim%multiple_box_events_policy)) /= 'soft_discard') return
    count_limit = int(sim%multiple_box_events_soft_discard_count_limit, i64)
    if (batch_count < 0_i64) error stop 'soft-discard batch count overflowed.'
    if (stats%multiple_box_events_soft_discarded > huge(projected_count) - batch_count) then
      projected_count = huge(projected_count)
    else
      projected_count = stats%multiple_box_events_soft_discarded + batch_count
    end if
    projected_abs_charge = stats%multiple_box_events_soft_discarded_abs_charge + batch_abs_charge
    count_exceeded = projected_count > count_limit
    charge_exceeded = .not. ieee_is_finite(projected_abs_charge) .or. &
                      projected_abs_charge > sim%multiple_box_events_soft_discard_abs_charge_limit
    if (.not. count_exceeded .and. .not. charge_exceeded) return

    if (mpi_is_root(mpi)) then
      write (error_unit, '(a,i0,a,i0,a,i0,a,es13.5,a,es13.5)') &
        'multiple_box_events soft-discard cumulative limit exceeded: batch=', batch_idx, &
        ' count=', projected_count, ' count_limit=', count_limit, &
        ' abs_charge_C=', projected_abs_charge, &
        ' abs_charge_limit_C=', sim%multiple_box_events_soft_discard_abs_charge_limit
      flush (error_unit)
    end if
    error stop 'multiple_box_events soft-discard cumulative limit exceeded.'
  end subroutine enforce_soft_discard_limits

  pure function particle_failure_name(status) result(name)
    integer(i32), intent(in) :: status
    character(len=32) :: name

    select case (status)
    case (collision_query_image_limit)
      name = 'image_limit'
    case (collision_query_index_range)
      name = 'index_range'
    case (collision_query_invalid_segment)
      name = 'invalid_segment'
    case (collision_query_grid_stalled)
      name = 'grid_stalled'
    case (particle_step_invalid_boundary)
      name = 'invalid_boundary'
    case (particle_step_multiple_box_events)
      name = 'multiple_box_events'
    case (particle_step_ambiguous_open_corner)
      name = 'ambiguous_open_corner'
    case default
      name = 'unknown'
    end select
  end function particle_failure_name

  pure function collision_failure_name(status) result(name)
    integer(i32), intent(in) :: status
    character(len=16) :: name

    select case (status)
    case (collision_query_image_limit)
      name = 'image_limit'
    case (collision_query_index_range)
      name = 'index_range'
    case (collision_query_invalid_segment)
      name = 'invalid_segment'
    case (collision_query_grid_stalled)
      name = 'grid_stalled'
    case default
      name = 'unknown'
    end select
  end function collision_failure_name

  module procedure commit_batch_charge
  real(dp) :: norm_dq, norm_q
  if (.not. all(ieee_is_finite(mesh%q_elem))) then
    error stop 'committed surface charge is not finite before batch update.'
  end if
  workspace%q_before = mesh%q_elem
  if (workspace%charge_candidate_ready) then
    if (.not. all(ieee_is_finite(workspace%candidate_charge))) then
      error stop 'adaptive candidate charge is not finite at commit.'
    end if
    mesh%q_elem = workspace%candidate_charge
  else
    workspace%dq = sum(workspace%dq_thread, dim=2) + sum(workspace%photo_emission_dq, dim=2)
    call mpi_allreduce_sum_real_dp_array(mpi, workspace%dq)
    call validate_finite_charge_addition( &
      mesh%q_elem, workspace%dq, 'batch surface-charge update' &
      )
    mesh%q_elem = mesh%q_elem + workspace%dq
  end if
  if (.not. all(ieee_is_finite(mesh%q_elem))) then
    error stop 'batch surface-charge update produced a non-finite charge.'
  end if
  call apply_surface_model_charge_relaxation(mesh, external_e, field_bc_mode=field_bc_mode)
  if (.not. all(ieee_is_finite(mesh%q_elem))) then
    error stop 'surface-model charge relaxation produced a non-finite charge.'
  end if
  call validate_finite_charge_addition( &
    mesh%q_elem, -workspace%q_before, 'batch surface-charge difference' &
    )
  workspace%dq = mesh%q_elem - workspace%q_before
  norm_dq = stable_l2_norm(workspace%dq)
  norm_q = stable_l2_norm(mesh%q_elem)
  rel = finite_nonnegative_ratio(norm_dq, max(norm_q, q_floor))
  workspace%charge_candidate_ready = .false.
  end procedure commit_batch_charge

  subroutine prepare_adaptive_charge_candidate(mesh, workspace, mpi)
    type(mesh_type), intent(in) :: mesh
    type(simulator_batch_workspace_type), intent(inout) :: workspace
    type(mpi_context), intent(in) :: mpi
    if (size(workspace%candidate_charge) /= mesh%nelem) then
      error stop 'adaptive candidate storage does not match the mesh.'
    end if
    workspace%dq = sum(workspace%dq_thread, dim=2) + sum(workspace%photo_emission_dq, dim=2)
    call mpi_allreduce_sum_real_dp_array(mpi, workspace%dq)
    call validate_finite_charge_addition( &
      mesh%q_elem, workspace%dq, 'adaptive candidate surface-charge update' &
      )
    workspace%candidate_charge = mesh%q_elem + workspace%dq
    if (.not. all(ieee_is_finite(workspace%candidate_charge))) then
      error stop 'adaptive candidate charge is not finite.'
    end if
    workspace%charge_candidate_ready = .true.
  end subroutine prepare_adaptive_charge_candidate

  !> Finiteな二項の加算が表現可能範囲を越えないことを、演算前に検証する。
  subroutine validate_finite_charge_addition(base, increment, operation)
    real(dp), intent(in) :: base(:), increment(:)
    character(len=*), intent(in) :: operation
    integer :: i

    if (size(base) /= size(increment)) then
      error stop trim(operation)//' array sizes do not match.'
    end if
    if (.not. all(ieee_is_finite(base)) .or. .not. all(ieee_is_finite(increment))) then
      error stop trim(operation)//' contains a non-finite operand.'
    end if
    do i = 1, size(base)
      if (increment(i) > 0.0_dp .and. base(i) > huge(base(i)) - increment(i)) then
        error stop trim(operation)//' overflowed.'
      end if
      if (increment(i) < 0.0_dp .and. base(i) < -huge(base(i)) - increment(i)) then
        error stop trim(operation)//' overflowed.'
      end if
    end do
  end subroutine validate_finite_charge_addition

  !> 二乗の中間overflowを避けて有限vectorのL2 normを返す。
  pure real(dp) function stable_l2_norm(values) result(norm)
    real(dp), intent(in) :: values(:)
    real(dp) :: scale, unit_norm

    if (size(values) == 0) then
      norm = 0.0_dp
      return
    end if
    scale = maxval(abs(values))
    if (scale == 0.0_dp) then
      norm = 0.0_dp
      return
    end if
    unit_norm = sqrt(sum((values/scale)*(values/scale)))
    if (scale > huge(norm)/unit_norm) then
      norm = huge(norm)
    else
      norm = scale*unit_norm
    end if
  end function stable_l2_norm

  !> 非負の有限比をoverflow時は最大有限値へ飽和させる。
  pure real(dp) function finite_nonnegative_ratio(numerator, denominator) result(ratio)
    real(dp), intent(in) :: numerator, denominator

    if (numerator == 0.0_dp) then
      ratio = 0.0_dp
    else if (denominator < 1.0_dp .and. numerator > huge(ratio)*denominator) then
      ratio = huge(ratio)
    else
      ratio = numerator/denominator
    end if
  end function finite_nonnegative_ratio

  subroutine apply_neutral_return_surface_closure(app, pcls_batch, fresh_particle_count, workspace, mpi)
    type(app_config), intent(in) :: app
    type(particles_soa), intent(in) :: pcls_batch
    integer(i32), intent(in) :: fresh_particle_count
    type(simulator_batch_workspace_type), intent(inout) :: workspace
    type(mpi_context), intent(in) :: mpi
    integer(i32) :: i, species_idx, elem_idx, n, terminal_count
    integer(i64) :: escaped_count, soft_count, invalid_count
    real(dp) :: macro_charge, emitted_charge, absorbed_charge, unresolved_charge
    real(dp) :: escaped_charge, soft_charge, invalid_charge
    real(dp) :: charge_scale, charge_tolerance, balance_tolerance
    real(dp) :: weight_scale, correction_charge, unresolved_fraction
    logical :: has_neutral_return

    n = app%n_particle_species
    has_neutral_return = .false.
    do species_idx = 1_i32, n
      if (.not. app%particle_species(species_idx)%enabled) cycle
      if (trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) == 'neutral_return') then
        has_neutral_return = .true.
      end if
    end do
    if (.not. has_neutral_return) return
    if (.not. app%sim%use_box) error stop 'neutral_return requires a finite box.'

    workspace%neutral_return_charge_values = 0.0_dp
    workspace%neutral_return_terminal_counts = 0_i64
    do i = 1_i32, pcls_batch%n
      species_idx = pcls_batch%species_id(i)
      if (trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) /= 'neutral_return') cycle
      macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
      if (.not. ieee_is_finite(macro_charge) .or. macro_charge >= 0.0_dp .or. i > fresh_particle_count) then
        workspace%neutral_return_charge_values(5*n + species_idx) = &
          workspace%neutral_return_charge_values(5*n + species_idx) + macro_charge
        workspace%neutral_return_terminal_counts(2*n + species_idx) = &
          workspace%neutral_return_terminal_counts(2*n + species_idx) + 1_i64
        cycle
      end if
      workspace%neutral_return_charge_values(species_idx) = &
        workspace%neutral_return_charge_values(species_idx) + macro_charge
      terminal_count = merge(1_i32, 0_i32, workspace%absorbed_flag(i)) + &
                       merge(1_i32, 0_i32, workspace%escaped_boundary_flag(i)) + &
                       merge(1_i32, 0_i32, workspace%soft_discarded_boundary_flag(i)) + &
                       merge(1_i32, 0_i32, pcls_batch%alive(i))
      if (terminal_count /= 1_i32) then
        workspace%neutral_return_charge_values(5*n + species_idx) = &
          workspace%neutral_return_charge_values(5*n + species_idx) + macro_charge
        workspace%neutral_return_terminal_counts(2*n + species_idx) = &
          workspace%neutral_return_terminal_counts(2*n + species_idx) + 1_i64
      else if (workspace%absorbed_flag(i)) then
        elem_idx = workspace%absorbed_element(i)
        if (elem_idx >= 1_i32 .and. elem_idx <= size(workspace%dq_thread, 1)) then
          workspace%neutral_return_charge_values(n + species_idx) = &
            workspace%neutral_return_charge_values(n + species_idx) + macro_charge
        else
          workspace%neutral_return_charge_values(5*n + species_idx) = &
            workspace%neutral_return_charge_values(5*n + species_idx) + macro_charge
          workspace%neutral_return_terminal_counts(2*n + species_idx) = &
            workspace%neutral_return_terminal_counts(2*n + species_idx) + 1_i64
        end if
      else if (workspace%escaped_boundary_flag(i)) then
        workspace%neutral_return_charge_values(3*n + species_idx) = &
          workspace%neutral_return_charge_values(3*n + species_idx) + macro_charge
        workspace%neutral_return_terminal_counts(species_idx) = &
          workspace%neutral_return_terminal_counts(species_idx) + 1_i64
      else if (workspace%soft_discarded_boundary_flag(i)) then
        workspace%neutral_return_charge_values(4*n + species_idx) = &
          workspace%neutral_return_charge_values(4*n + species_idx) + macro_charge
        workspace%neutral_return_terminal_counts(n + species_idx) = &
          workspace%neutral_return_terminal_counts(n + species_idx) + 1_i64
      else
        workspace%neutral_return_charge_values(2*n + species_idx) = &
          workspace%neutral_return_charge_values(2*n + species_idx) + macro_charge
      end if
    end do

    call mpi_allreduce_sum_real_dp_array(mpi, workspace%neutral_return_charge_values)
    call mpi_allreduce_sum_i64_array(mpi, workspace%neutral_return_terminal_counts)
    workspace%neutral_return_emitted_charge = workspace%neutral_return_charge_values(1:n)
    workspace%neutral_return_absorbed_charge = workspace%neutral_return_charge_values(n + 1:2*n)
    workspace%neutral_return_unresolved_charge = workspace%neutral_return_charge_values(2*n + 1:3*n)

    do species_idx = 1_i32, n
      if (trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) /= 'neutral_return') cycle
      emitted_charge = workspace%neutral_return_emitted_charge(species_idx)
      absorbed_charge = workspace%neutral_return_absorbed_charge(species_idx)
      unresolved_charge = workspace%neutral_return_unresolved_charge(species_idx)
      escaped_charge = workspace%neutral_return_charge_values(3*n + species_idx)
      soft_charge = workspace%neutral_return_charge_values(4*n + species_idx)
      invalid_charge = workspace%neutral_return_charge_values(5*n + species_idx)
      escaped_count = workspace%neutral_return_terminal_counts(species_idx)
      soft_count = workspace%neutral_return_terminal_counts(n + species_idx)
      invalid_count = workspace%neutral_return_terminal_counts(2*n + species_idx)
      if (escaped_count > 0_i64 .or. soft_count > 0_i64 .or. invalid_count > 0_i64) then
        write (error_unit, '(a,i0,3(a,i0),3(a,es13.5))') &
          'neutral_return has unsupported terminal outcome for species ', species_idx, &
          ': escaped=', escaped_count, ' soft=', soft_count, ' invalid=', invalid_count, &
          ' escaped_C=', escaped_charge, ' soft_C=', soft_charge, ' invalid_C=', invalid_charge
        error stop 'neutral_return terminal outcome is unsupported.'
      end if
      charge_scale = max(abs(emitted_charge), abs(absorbed_charge), abs(unresolved_charge), tiny(1.0_dp))
      charge_tolerance = 4096.0_dp*epsilon(1.0_dp)*charge_scale
      balance_tolerance = max(charge_tolerance, sqrt(epsilon(1.0_dp))*charge_scale)
      if (abs(emitted_charge) <= charge_tolerance) cycle
      if (emitted_charge >= 0.0_dp .or. absorbed_charge >= -charge_tolerance .or. &
          abs(emitted_charge - absorbed_charge - unresolved_charge) > balance_tolerance) then
        error stop 'neutral_return charge balance is invalid.'
      end if
      weight_scale = emitted_charge/absorbed_charge
      correction_charge = emitted_charge - absorbed_charge
      unresolved_fraction = unresolved_charge/emitted_charge
      if (.not. all(ieee_is_finite([weight_scale, correction_charge, unresolved_fraction])) .or. &
          unresolved_fraction > neutral_return_max_unresolved_fraction + sqrt(epsilon(1.0_dp))) then
        error stop 'neutral_return unresolved fraction exceeds the applicability limit.'
      end if
      workspace%neutral_return_weight_scale(species_idx) = max(1.0_dp, weight_scale)
      workspace%neutral_return_correction(species_idx) = correction_charge
      workspace%neutral_return_unresolved_fraction(species_idx) = max(0.0_dp, unresolved_fraction)
    end do

    do i = 1_i32, fresh_particle_count
      species_idx = pcls_batch%species_id(i)
      if (trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) /= 'neutral_return') cycle
      if (.not. workspace%absorbed_flag(i)) cycle
      elem_idx = workspace%absorbed_element(i)
      macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
      workspace%dq_thread(elem_idx, 1) = workspace%dq_thread(elem_idx, 1) + &
                                         (workspace%neutral_return_weight_scale(species_idx) - 1.0_dp)*macro_charge
    end do
  end subroutine apply_neutral_return_surface_closure

  subroutine apply_fixed_surface_current_closure( &
    app, current_model, pcls_batch, fresh_particle_count, workspace, mpi &
    )
    type(app_config), intent(in) :: app
    type(surface_closure_contract_type), intent(in) :: current_model
    type(particles_soa), intent(in) :: pcls_batch
    integer(i32), intent(in) :: fresh_particle_count
    type(simulator_batch_workspace_type), intent(inout) :: workspace
    type(mpi_context), intent(in) :: mpi
    integer(i32) :: i, species_idx, elem_idx, n
    real(dp) :: macro_charge, raw_charge, target_current, target_charge, weight_scale, correction
    real(dp) :: charge_tolerance
    logical :: has_fixed_current

    n = app%n_particle_species
    has_fixed_current = .false.
    do species_idx = 1_i32, n
      if (.not. app%particle_species(species_idx)%enabled) cycle
      if (fixed_current_species_active(app, current_model, species_idx)) then
        has_fixed_current = .true.
      end if
    end do
    if (.not. has_fixed_current) return

    workspace%fixed_current_charge_values = 0.0_dp
    do i = 1_i32, fresh_particle_count
      species_idx = pcls_batch%species_id(i)
      if (.not. fixed_current_species_active(app, current_model, species_idx)) cycle
      macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
      if (workspace%absorbed_flag(i)) then
        workspace%fixed_current_charge_values(species_idx) = &
          workspace%fixed_current_charge_values(species_idx) + macro_charge
      else if (workspace%escaped_boundary_flag(i)) then
        workspace%fixed_current_charge_values(2*n + species_idx) = &
          workspace%fixed_current_charge_values(2*n + species_idx) + macro_charge
      end if
    end do
    do species_idx = 1_i32, n
      if (.not. fixed_current_species_active(app, current_model, species_idx)) cycle
      workspace%fixed_current_charge_values(n + species_idx) = sum(workspace%photo_emission_dq(:, species_idx))
    end do
    call mpi_allreduce_sum_real_dp_array(mpi, workspace%fixed_current_charge_values)

    do species_idx = 1_i32, n
      if (.not. fixed_current_species_active(app, current_model, species_idx)) cycle
      if (app%particle_species(species_idx)%has_target_absorbed_current_a .or. &
          current_model%has_absorbed_target(species_idx)) then
        raw_charge = workspace%fixed_current_charge_values(species_idx)
        if (current_model%has_absorbed_target(species_idx)) then
          target_current = current_model%absorbed_current_a(species_idx)
        else
          target_current = app%particle_species(species_idx)%target_absorbed_current_a
        end if
        if (.not. all(ieee_is_finite([raw_charge, target_current, app%sim%batch_duration]))) then
          error stop 'fixed_current absorbed raw/current/duration value is not finite.'
        end if
        if (app%sim%batch_duration > 1.0_dp .and. &
            abs(target_current) > huge(target_charge)/app%sim%batch_duration) then
          error stop 'fixed_current absorbed target charge overflowed for this batch duration.'
        end if
        target_charge = target_current*app%sim%batch_duration
        if (target_current /= 0.0_dp .and. target_charge == 0.0_dp) then
          error stop 'fixed_current absorbed target charge underflowed for this batch duration.'
        end if
        charge_tolerance = 4096.0_dp*epsilon(1.0_dp)*max(abs(raw_charge), abs(target_charge), tiny(1.0_dp))
        if (abs(raw_charge) <= charge_tolerance) then
          if (abs(target_charge) > charge_tolerance) then
            error stop 'fixed_current cannot map a nonzero absorbed target onto an empty raw channel.'
          end if
          weight_scale = 1.0_dp
          correction = 0.0_dp
        else
          weight_scale = target_charge/raw_charge
          correction = target_charge - raw_charge
        end if
        if (.not. ieee_is_finite(weight_scale) .or. weight_scale < 0.0_dp) then
          error stop 'fixed_current absorbed scale is invalid.'
        end if
        workspace%fixed_absorbed_target_charge(species_idx) = target_charge
        workspace%fixed_absorbed_weight_scale(species_idx) = weight_scale
        workspace%fixed_current_correction(species_idx) = &
          workspace%fixed_current_correction(species_idx) + correction
        do i = 1_i32, fresh_particle_count
          if (pcls_batch%species_id(i) /= species_idx .or. .not. workspace%absorbed_flag(i)) cycle
          elem_idx = workspace%absorbed_element(i)
          macro_charge = pcls_batch%q(i)*pcls_batch%w(i)
          workspace%dq_thread(elem_idx, 1) = workspace%dq_thread(elem_idx, 1) + &
                                             (weight_scale - 1.0_dp)*macro_charge
        end do
      end if

      if (app%particle_species(species_idx)%has_target_emission_current_a .or. &
          current_model%has_emission_target(species_idx)) then
        raw_charge = workspace%fixed_current_charge_values(n + species_idx)
        if (current_model%has_emission_target(species_idx)) then
          target_current = current_model%emission_current_a(species_idx)
        else
          target_current = app%particle_species(species_idx)%target_emission_current_a
        end if
        if (.not. all(ieee_is_finite([raw_charge, target_current, app%sim%batch_duration]))) then
          error stop 'fixed_current emission raw/current/duration value is not finite.'
        end if
        if (app%sim%batch_duration > 1.0_dp .and. &
            abs(target_current) > huge(target_charge)/app%sim%batch_duration) then
          error stop 'fixed_current emission target charge overflowed for this batch duration.'
        end if
        target_charge = target_current*app%sim%batch_duration
        if (target_current /= 0.0_dp .and. target_charge == 0.0_dp) then
          error stop 'fixed_current emission target charge underflowed for this batch duration.'
        end if
        charge_tolerance = 4096.0_dp*epsilon(1.0_dp)*max(abs(raw_charge), abs(target_charge), tiny(1.0_dp))
        if (abs(raw_charge) <= charge_tolerance) then
          if (abs(target_charge) > charge_tolerance) then
            error stop 'fixed_current cannot map a nonzero emission target onto an empty raw channel.'
          end if
          weight_scale = 1.0_dp
          correction = 0.0_dp
        else
          weight_scale = target_charge/raw_charge
          correction = target_charge - raw_charge
        end if
        if (.not. ieee_is_finite(weight_scale) .or. weight_scale < 0.0_dp) then
          error stop 'fixed_current emission scale is invalid.'
        end if
        workspace%fixed_emission_target_charge(species_idx) = target_charge
        workspace%fixed_emission_weight_scale(species_idx) = weight_scale
        workspace%fixed_current_correction(species_idx) = &
          workspace%fixed_current_correction(species_idx) + correction
        workspace%photo_emission_dq(:, species_idx) = &
          weight_scale*workspace%photo_emission_dq(:, species_idx)
      end if

      if (current_model%has_escape_target(species_idx)) then
        raw_charge = workspace%fixed_current_charge_values(2*n + species_idx)
        target_current = current_model%escaped_particle_current_a(species_idx)
        if (.not. all(ieee_is_finite([raw_charge, target_current, app%sim%batch_duration]))) then
          error stop 'fixed_current escape raw/current/duration value is not finite.'
        end if
        if (app%sim%batch_duration > 1.0_dp .and. &
            abs(target_current) > huge(target_charge)/app%sim%batch_duration) then
          error stop 'fixed_current escape target charge overflowed for this batch duration.'
        end if
        target_charge = target_current*app%sim%batch_duration
        if (target_current /= 0.0_dp .and. target_charge == 0.0_dp) then
          error stop 'fixed_current escape target charge underflowed for this batch duration.'
        end if
        if (target_charge /= 0.0_dp .and. &
            sign(1.0_dp, target_charge) /= sign(1.0_dp, app%particle_species(species_idx)%q_particle)) then
          error stop 'fixed_current escape target sign must match the escaped particle charge.'
        end if
        workspace%fixed_escape_target_charge(species_idx) = target_charge
        workspace%fixed_escape_correction(species_idx) = target_charge - raw_charge
        if (.not. ieee_is_finite(workspace%fixed_escape_correction(species_idx))) then
          error stop 'fixed_current escape correction is not finite.'
        end if
      end if
    end do
  end subroutine apply_fixed_surface_current_closure

  logical function fixed_current_species_active(app, current_model, species_idx) result(active)
    type(app_config), intent(in) :: app
    type(surface_closure_contract_type), intent(in) :: current_model
    integer(i32), intent(in) :: species_idx

    active = trim(lower_ascii(app%particle_species(species_idx)%surface_charge_closure)) == 'fixed_current' .or. &
             current_model%has_absorbed_target(species_idx) .or. &
             current_model%has_emission_target(species_idx) .or. &
             current_model%has_escape_target(species_idx)
  end function fixed_current_species_active

  subroutine record_batch_initial_charge(app, pcls_batch, fresh_particle_count, ledger)
    type(app_config), intent(in) :: app
    type(particles_soa), intent(in) :: pcls_batch
    integer(i32), intent(in) :: fresh_particle_count
    type(charge_ledger_type), intent(inout) :: ledger
    integer(i32) :: i, species_idx
    real(dp) :: macro_charge
    do i = 1_i32, fresh_particle_count
      species_idx = pcls_batch%species_id(i)
      macro_charge = checked_macro_charge(pcls_batch%q(i), pcls_batch%w(i), 'initial batch charge ledger')
      if (trim(lower_ascii(app%particle_species(species_idx)%source_mode)) == 'photo_raycast') then
        call checked_accumulate_charge( &
          ledger%emitted_from_surface(species_idx), macro_charge, 'local emitted charge ledger' &
          )
        ledger%emitted_count(species_idx) = ledger%emitted_count(species_idx) + 1_i64
      else
        call checked_accumulate_charge( &
          ledger%injected_from_remote(species_idx), macro_charge, 'local injected charge ledger' &
          )
        ledger%injected_count(species_idx) = ledger%injected_count(species_idx) + 1_i64
      end if
    end do
  end subroutine record_batch_initial_charge

  subroutine record_batch_outcome_charge( &
    pcls_batch, escaped_boundary_flag, absorbed_flag, soft_discarded_boundary_flag, ledger &
    )
    type(particles_soa), intent(in) :: pcls_batch
    logical, intent(in) :: escaped_boundary_flag(:), absorbed_flag(:), soft_discarded_boundary_flag(:)
    type(charge_ledger_type), intent(inout) :: ledger
    integer(i32) :: i, species_idx
    real(dp) :: macro_charge
    do i = 1_i32, pcls_batch%n
      species_idx = pcls_batch%species_id(i)
      macro_charge = checked_macro_charge(pcls_batch%q(i), pcls_batch%w(i), 'outcome batch charge ledger')
      if (absorbed_flag(i)) then
        call checked_accumulate_charge( &
          ledger%absorbed_on_surface(species_idx), macro_charge, 'local absorbed charge ledger' &
          )
        ledger%absorbed_count(species_idx) = ledger%absorbed_count(species_idx) + 1_i64
      else if (escaped_boundary_flag(i)) then
        call checked_accumulate_charge( &
          ledger%escaped_to_infinity(species_idx), macro_charge, 'local escaped charge ledger' &
          )
        ledger%escaped_count(species_idx) = ledger%escaped_count(species_idx) + 1_i64
      else if (soft_discarded_boundary_flag(i) .or. pcls_batch%alive(i)) then
        call checked_accumulate_charge( &
          ledger%discarded_unresolved(species_idx), macro_charge, 'local discarded charge ledger' &
          )
        ledger%discarded_unresolved_count(species_idx) = ledger%discarded_unresolved_count(species_idx) + 1_i64
      end if
    end do
  end subroutine record_batch_outcome_charge

  subroutine reduce_charge_ledger_fluxes(ledger, mpi, workspace)
    type(charge_ledger_type), intent(inout) :: ledger
    type(mpi_context), intent(in) :: mpi
    type(simulator_batch_workspace_type), intent(inout) :: workspace
    integer(i32) :: n
    n = ledger%nspecies
    workspace%ledger_charge_values = [ &
                                     ledger%injected_from_remote, ledger%emitted_from_surface, ledger%absorbed_on_surface, &
                                     ledger%escaped_to_infinity, ledger%discarded_unresolved &
                                     ]
    if (.not. all(ieee_is_finite(workspace%ledger_charge_values))) then
      error stop 'local batch charge ledger contains non-finite fluxes before MPI reduction.'
    end if
    call mpi_allreduce_sum_real_dp_array(mpi, workspace%ledger_charge_values)
    if (.not. all(ieee_is_finite(workspace%ledger_charge_values))) then
      error stop 'global batch charge ledger overflowed during MPI reduction.'
    end if
    ledger%injected_from_remote = workspace%ledger_charge_values(1:n)
    ledger%emitted_from_surface = workspace%ledger_charge_values(n + 1:2*n)
    ledger%absorbed_on_surface = workspace%ledger_charge_values(2*n + 1:3*n)
    ledger%escaped_to_infinity = workspace%ledger_charge_values(3*n + 1:4*n)
    ledger%discarded_unresolved = workspace%ledger_charge_values(4*n + 1:5*n)
    workspace%ledger_count_values = [ &
                                    ledger%injected_count, ledger%emitted_count, ledger%absorbed_count, ledger%escaped_count, &
                                    ledger%discarded_unresolved_count &
                                    ]
    if (any(workspace%ledger_count_values < 0_i64)) then
      error stop 'local batch charge ledger contains invalid counts before MPI reduction.'
    end if
    call mpi_allreduce_sum_i64_array(mpi, workspace%ledger_count_values)
    if (any(workspace%ledger_count_values < 0_i64)) then
      error stop 'global batch charge ledger count overflowed during MPI reduction.'
    end if
    ledger%injected_count = workspace%ledger_count_values(1:n)
    ledger%emitted_count = workspace%ledger_count_values(n + 1:2*n)
    ledger%absorbed_count = workspace%ledger_count_values(2*n + 1:3*n)
    ledger%escaped_count = workspace%ledger_count_values(3*n + 1:4*n)
    ledger%discarded_unresolved_count = workspace%ledger_count_values(4*n + 1:5*n)
  end subroutine reduce_charge_ledger_fluxes

  integer(i32) function matching_species_index(app, species_key) result(index)
    type(app_config), intent(in) :: app
    character(len=*), intent(in) :: species_key
    integer(i32) :: species_idx

    index = 0_i32
    do species_idx = 1_i32, app%n_particle_species
      if (trim(app%particle_species(species_idx)%species_key) /= trim(species_key)) cycle
      index = species_idx
      return
    end do
    error stop 'matching-plane species role was not found in particle species.'
  end function matching_species_index

  subroutine solve_matching_implicit_zero_mode( &
    provider, mpi, displacement_before, duration, displacement_min, displacement_max, feedback_reference, &
    electron_charge, ion_charge, photoelectron_active, photoelectron_charge, displacement_after, response_after &
    )
    type(matching_plane_response_provider_type), intent(inout) :: provider
    type(mpi_context), intent(in) :: mpi
    real(dp), intent(in) :: displacement_before, duration, displacement_min, displacement_max
    real(dp), intent(in) :: feedback_reference(4)
    real(dp), intent(in) :: electron_charge, ion_charge, photoelectron_charge
    logical, intent(in) :: photoelectron_active
    real(dp), intent(out) :: displacement_after, response_after(6)

    real(dp) :: lower, upper, midpoint, lower_residual, upper_residual, midpoint_residual
    real(dp) :: response(6), current_density
    integer(i32) :: iteration, status
    character(len=512) :: message

    if (.not. all(ieee_is_finite([ &
                                 displacement_before, duration, displacement_min, displacement_max, &
                                 feedback_reference, electron_charge, ion_charge, photoelectron_charge &
                                 ])) .or. duration <= 0.0_dp .or. displacement_min >= displacement_max) then
      error stop 'implicit matching-plane zero-mode inputs are invalid.'
    end if
    lower = displacement_min
    upper = displacement_max
    call evaluate_matching_implicit_residual( &
      provider, mpi, lower, displacement_before, duration, feedback_reference, &
      electron_charge, ion_charge, photoelectron_active, photoelectron_charge, &
      lower_residual, response, current_density, status, message &
      )
    if (status /= matching_plane_provider_ok) then
      error stop 'implicit matching-plane lower endpoint failed: '//trim(message)
    end if
    call evaluate_matching_implicit_residual( &
      provider, mpi, upper, displacement_before, duration, feedback_reference, &
      electron_charge, ion_charge, photoelectron_active, photoelectron_charge, &
      upper_residual, response, current_density, status, message &
      )
    if (status /= matching_plane_provider_ok) then
      error stop 'implicit matching-plane upper endpoint failed: '//trim(message)
    end if
    if (lower_residual > 0.0_dp .or. upper_residual < 0.0_dp) then
      error stop 'implicit matching-plane zero-mode root is not bracketed by the response table.'
    end if

    do iteration = 1_i32, 64_i32
      midpoint = 0.5_dp*(lower + upper)
      call evaluate_matching_implicit_residual( &
        provider, mpi, midpoint, displacement_before, duration, feedback_reference, &
        electron_charge, ion_charge, photoelectron_active, photoelectron_charge, &
        midpoint_residual, response, current_density, status, message &
        )
      if (status /= matching_plane_provider_ok) then
        error stop 'implicit matching-plane midpoint failed: '//trim(message)
      end if
      if (midpoint_residual <= 0.0_dp) then
        lower = midpoint
      else
        upper = midpoint
      end if
    end do
    displacement_after = 0.5_dp*(lower + upper)
    call evaluate_matching_implicit_residual( &
      provider, mpi, displacement_after, displacement_before, duration, feedback_reference, &
      electron_charge, ion_charge, photoelectron_active, photoelectron_charge, midpoint_residual, response_after, &
      current_density, status, message &
      )
    if (status /= matching_plane_provider_ok) then
      error stop 'implicit matching-plane final evaluation failed: '//trim(message)
    end if
  end subroutine solve_matching_implicit_zero_mode

  subroutine evaluate_matching_implicit_residual( &
    provider, mpi, displacement, displacement_before, duration, feedback_reference, &
    electron_charge, ion_charge, photoelectron_active, photoelectron_charge, &
    residual, response, current_density, status, message &
    )
    type(matching_plane_response_provider_type), intent(inout) :: provider
    type(mpi_context), intent(in) :: mpi
    real(dp), intent(in) :: displacement, displacement_before, duration, feedback_reference(4)
    real(dp), intent(in) :: electron_charge, ion_charge, photoelectron_charge
    logical, intent(in) :: photoelectron_active
    real(dp), intent(out) :: residual, response(6), current_density
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp) :: input(5), escape_fraction, escape_flux, barrier_energy_ev

    input = [displacement, feedback_reference]
    call provider%evaluate(input, mpi, response, status, message)
    residual = 0.0_dp
    current_density = 0.0_dp
    if (status /= matching_plane_provider_ok) return
    escape_fraction = 0.0_dp
    escape_flux = 0.0_dp
    current_density = electron_charge*response(2) + ion_charge*response(3)
    if (photoelectron_active) then
      barrier_energy_ev = response(1) - response(6)
      if (.not. ieee_is_finite(barrier_energy_ev) .or. barrier_energy_ev < 0.0_dp .or. &
          feedback_reference(2) <= 0.0_dp) then
        status = 1_i32
        message = 'implicit matching-plane PE barrier or reference energy is invalid.'
        return
      end if
      escape_fraction = exp(-barrier_energy_ev/feedback_reference(2))
      escape_flux = feedback_reference(1)*escape_fraction
      current_density = current_density - photoelectron_charge*escape_flux
    end if
    residual = displacement - displacement_before - duration*current_density
    if (.not. all(ieee_is_finite([escape_fraction, escape_flux, current_density, residual]))) then
      status = 1_i32
      message = 'implicit matching-plane zero-mode residual is not finite.'
    end if
  end subroutine evaluate_matching_implicit_residual

  subroutine configure_matching_surface_closure( &
    contract, electron_idx, ion_idx, photoelectron_idx, photoelectron_active, response, implicit_zero_mode, area_m2, &
    feedback_reference, electron_charge, ion_charge, photoelectron_charge, &
    photoelectron_emission_current_density &
    )
    type(surface_closure_contract_type), intent(inout) :: contract
    integer(i32), intent(in) :: electron_idx, ion_idx, photoelectron_idx
    real(dp), intent(in) :: response(6)
    logical, intent(in) :: implicit_zero_mode, photoelectron_active
    real(dp), intent(in) :: area_m2, feedback_reference(4)
    real(dp), intent(in) :: electron_charge, ion_charge, photoelectron_charge
    real(dp), intent(in) :: photoelectron_emission_current_density

    real(dp) :: barrier_energy_ev, emission_flux, escape_flux, return_flux

    if (.not. all(ieee_is_finite(response))) then
      error stop 'matching-plane response values must be finite.'
    end if
    if (any(response(2:3) < 0.0_dp)) then
      error stop 'matching-plane inward number fluxes must be nonnegative.'
    end if
    if (.not. allocated(contract%has_inflow_number_flux) .or. &
        .not. allocated(contract%inflow_number_flux_m2_s)) then
      error stop 'matching-plane closure number-flux storage is not initialized.'
    end if
    contract%active = .true.
    contract%has_absorbed_target = .false.
    contract%has_emission_target = .false.
    contract%has_escape_target = .false.
    contract%has_inflow_kinetic_map = .false.
    contract%has_outflow_kinetic_barrier = .false.
    contract%has_inflow_number_flux = .false.
    contract%absorbed_current_a = 0.0_dp
    contract%emission_current_a = 0.0_dp
    contract%escaped_particle_current_a = 0.0_dp
    contract%inflow_reservoir_potential_v = 0.0_dp
    contract%inflow_access_potential_v = 0.0_dp
    contract%inflow_kinetic_face = 0_i32
    contract%outflow_barrier_potential_v = 0.0_dp
    contract%outflow_barrier_face = 0_i32
    contract%inflow_number_flux_m2_s = 0.0_dp

    contract%has_inflow_number_flux([electron_idx, ion_idx]) = .true.
    contract%inflow_number_flux_m2_s(electron_idx) = response(2)
    contract%inflow_number_flux_m2_s(ion_idx) = response(3)
    contract%has_inflow_kinetic_map([electron_idx, ion_idx]) = .true.
    contract%inflow_access_potential_v(electron_idx) = response(4)
    contract%inflow_access_potential_v(ion_idx) = response(5)
    contract%inflow_kinetic_face([electron_idx, ion_idx]) = 6_i32
    ! The response inward fluxes already include all outer-sheath return of
    ! ambient species. Reflecting ambient outflow locally would count it twice.
    if (photoelectron_active) then
      contract%has_outflow_kinetic_barrier(photoelectron_idx) = .true.
      contract%outflow_barrier_potential_v(photoelectron_idx) = response(6)
      contract%outflow_barrier_face(photoelectron_idx) = 6_i32
    end if
    if (implicit_zero_mode) then
      if (.not. ieee_is_finite(area_m2) .or. area_m2 <= 0.0_dp) then
        error stop 'implicit matching-plane current targets are invalid.'
      end if
      contract%has_absorbed_target([electron_idx, ion_idx]) = .true.
      contract%absorbed_current_a(electron_idx) = electron_charge*response(2)*area_m2
      contract%absorbed_current_a(ion_idx) = ion_charge*response(3)*area_m2
      if (photoelectron_active) then
        barrier_energy_ev = response(1) - response(6)
        escape_flux = feedback_reference(1)*exp(-barrier_energy_ev/feedback_reference(2))
        emission_flux = photoelectron_emission_current_density/(-photoelectron_charge)
        return_flux = emission_flux - escape_flux
        if (.not. all(ieee_is_finite([barrier_energy_ev, emission_flux, escape_flux, return_flux])) .or. &
            barrier_energy_ev < 0.0_dp .or. escape_flux < 0.0_dp .or. return_flux < 0.0_dp) then
          error stop 'implicit matching-plane current targets are invalid.'
        end if
        contract%has_absorbed_target(photoelectron_idx) = .true.
        contract%absorbed_current_a(photoelectron_idx) = photoelectron_charge*return_flux*area_m2
        contract%has_emission_target(photoelectron_idx) = .true.
        contract%emission_current_a(photoelectron_idx) = photoelectron_emission_current_density*area_m2
        contract%has_escape_target(photoelectron_idx) = .true.
        contract%escaped_particle_current_a(photoelectron_idx) = photoelectron_charge*escape_flux*area_m2
      end if
    end if
  end subroutine configure_matching_surface_closure

  subroutine resolve_matching_observed_feedback( &
    moments, electron_idx, ion_idx, photoelectron_idx, photoelectron_active, area_m2, duration_s, &
    feedback, return_flux, escape_flux &
    )
    real(dp), intent(in) :: moments(:, :)
    integer(i32), intent(in) :: electron_idx, ion_idx, photoelectron_idx
    logical, intent(in) :: photoelectron_active
    real(dp), intent(in) :: area_m2, duration_s
    real(dp), intent(out) :: feedback(4), return_flux, escape_flux
    real(dp) :: normalization, photoelectron_outward_number

    if (size(moments, 1) < 4 .or. any(.not. ieee_is_finite(moments)) .or. any(moments < 0.0_dp)) then
      error stop 'matching-plane particle moments are invalid.'
    end if
    normalization = area_m2*duration_s
    if (.not. ieee_is_finite(normalization) .or. normalization <= 0.0_dp) then
      error stop 'matching-plane flux normalization must be finite and positive.'
    end if
    feedback(1:2) = 0.0_dp
    return_flux = 0.0_dp
    escape_flux = 0.0_dp
    if (photoelectron_active) then
      photoelectron_outward_number = moments(1, photoelectron_idx)
      feedback(1) = photoelectron_outward_number/normalization
      if (photoelectron_outward_number > 0.0_dp) then
        feedback(2) = moments(2, photoelectron_idx)/(photoelectron_outward_number*qe)
      end if
      return_flux = moments(3, photoelectron_idx)/normalization
      escape_flux = moments(4, photoelectron_idx)/normalization
    end if
    feedback(3) = moments(1, electron_idx)/normalization
    feedback(4) = moments(1, ion_idx)/normalization
    if (.not. all(ieee_is_finite([feedback, return_flux, escape_flux]))) then
      error stop 'matching-plane observed feedback is not finite.'
    end if
  end subroutine resolve_matching_observed_feedback

  !> finiteな粒子電荷とweightから、overflow/zero-underflowを拒否してmacro chargeを返す。
  real(dp) function checked_macro_charge(particle_charge, particle_weight, context) result(macro_charge)
    real(dp), intent(in) :: particle_charge, particle_weight
    character(len=*), intent(in) :: context

    if (.not. ieee_is_finite(particle_charge) .or. .not. ieee_is_finite(particle_weight)) then
      error stop trim(context)//' has a non-finite particle charge or weight.'
    end if
    if (abs(particle_weight) > 1.0_dp .and. abs(particle_charge) > huge(macro_charge)/abs(particle_weight)) then
      error stop trim(context)//' macro charge overflowed.'
    end if
    macro_charge = particle_charge*particle_weight
    if (.not. ieee_is_finite(macro_charge)) error stop trim(context)//' macro charge is not finite.'
    if (particle_charge /= 0.0_dp .and. particle_weight /= 0.0_dp .and. macro_charge == 0.0_dp) then
      error stop trim(context)//' macro charge underflowed to zero.'
    end if
  end function checked_macro_charge

end submodule bem_simulator_loop

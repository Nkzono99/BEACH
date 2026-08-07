program test_periodic2_cached_snapshot
!$ use omp_lib
  use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int, c_loc, c_ptr
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: eps0, qe
  use bem_types, only: mesh_type, sim_config, sim_stats, bc_open, bc_periodic
  use bem_mesh, only: init_mesh, prepare_periodic2_collision_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config
  use bem_charge_ledger, only: charge_ledger_type
  use bem_field_kernel_c, only: beach_kernel_build, beach_kernel_create, beach_kernel_destroy, &
                                beach_kernel_eval_e, beach_kernel_eval_phi, &
                                beach_kernel_get_periodic_cache_info, beach_kernel_ok, &
                                beach_kernel_set_periodic_cache, beach_kernel_update_charges
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_periodic_zero_mode_eval, only: eval_periodic_zero_mode, zero_mode_trace_plus
  use bem_coulomb_fmm_periodic_nonzero_reference, only: eval_periodic_nonzero_panel_reference
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, &
                          assert_allclose_1d, assert_equal_i32, assert_equal_i64, delete_file_if_exists, &
                          remove_empty_directory
  implicit none

  interface
    integer(c_int) function beach_zero_mode_create(handle) bind(C, name='beach_zero_mode_create') result(status)
      import :: c_int, c_ptr
      type(c_ptr), intent(out) :: handle
    end function beach_zero_mode_create

    integer(c_int) function beach_zero_mode_destroy(handle) bind(C, name='beach_zero_mode_destroy') result(status)
      import :: c_int, c_ptr
      type(c_ptr), value :: handle
    end function beach_zero_mode_destroy

    integer(c_int) function beach_zero_mode_build(handle, nsrc, source_heights_ptr, area_xy) &
      bind(C, name='beach_zero_mode_build') result(status)
      import :: c_double, c_int, c_ptr
      type(c_ptr), value :: handle
      integer(c_int), value :: nsrc
      type(c_ptr), value :: source_heights_ptr
      real(c_double), value :: area_xy
    end function beach_zero_mode_build

    integer(c_int) function beach_zero_mode_update( &
      handle, nsrc, charge_ptr, e_bottom, z_gauge, phi_gauge &
      ) bind(C, name='beach_zero_mode_update') result(status)
      import :: c_double, c_int, c_ptr
      type(c_ptr), value :: handle
      integer(c_int), value :: nsrc
      type(c_ptr), value :: charge_ptr
      real(c_double), value :: e_bottom, z_gauge, phi_gauge
    end function beach_zero_mode_update

    integer(c_int) function beach_zero_mode_eval(handle, ntarget, z_ptr, trace, phi_ptr, ez_ptr) &
      bind(C, name='beach_zero_mode_eval') result(status)
      import :: c_int, c_ptr
      type(c_ptr), value :: handle
      integer(c_int), value :: ntarget
      type(c_ptr), value :: z_ptr
      integer(c_int), value :: trace
      type(c_ptr), value :: phi_ptr, ez_ptr
    end function beach_zero_mode_eval
  end interface

  type(mesh_type) :: mesh
  type(sim_config) :: sim
  type(field_physics_config) :: field_config
  type(periodic2_physics_config) :: periodic_config
  type(panel_kernel_config) :: panel_config
  type(electrostatic_snapshot_type) :: snapshot
  real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2), target(3)
  real(dp) :: total_field(3), expected_field(3), nonzero_field(3), zero_field, zero_potential
  real(dp) :: total_potential, expected_potential, nonzero_potential
  real(dp) :: reference_field(3), reference_potential, field_error, potential_error, charge_scale
  character(len=512) :: cache_path, cache_dir, team_mismatch_path
  character(len=64) :: run_mode
  integer(i64) :: clock_count

  call get_command_argument(1, run_mode)
  if (trim(run_mode) == '--adaptive-team-mismatch-probe') then
    call run_adaptive_team_mismatch_probe()
    error stop 'adaptive team-mismatch probe unexpectedly completed'
  end if
  call system_clock(count=clock_count)
  write (cache_dir, '(a,i0)') 'test_periodic2_cached_snapshot_tmp_', clock_count
  team_mismatch_path = trim(cache_dir)//'/adaptive_team_mismatch.log'
  call configure_fixture(mesh, sim, field_config, periodic_config, panel_config, v0, v1, v2)
  call snapshot%init(mesh, sim, field_config, periodic_config, panel_config)
  cache_path = snapshot%nonzero_solver%fmm_core_plan%periodic_cache_path
  call snapshot%refresh(mesh)
  target = [0.37_dp, 0.61_dp, 0.42_dp]

  call test_init(5)
  call test_begin('cached_snapshot_composes_kneq0_and_k0_once')
  call assert_true(snapshot%use_cached_kneq0 .and. snapshot%use_zero_mode, 'cached split flags must be active')
  call assert_true( &
    trim(snapshot%diagnostics%periodic_cache_fingerprint) == &
    trim(snapshot%nonzero_solver%fmm_core_plan%periodic_cache_fingerprint), &
    'snapshot must expose the periodic cache fingerprint' &
    )
  call evaluate_components(expected_field, expected_potential, nonzero_field, nonzero_potential, zero_field, zero_potential)
  call snapshot%eval_local_e(mesh, target, total_field)
  call snapshot%eval_local_phi(mesh, sim, target, total_potential)
  call assert_true(abs(zero_field) > 1.0e-12_dp, 'fixture must exercise a nonzero k=0 field')
  call assert_allclose_1d(total_field, expected_field, 1.0e-10_dp, 'snapshot field composition mismatch')
  call assert_close_dp(total_potential, expected_potential, 1.0e-10_dp, 'snapshot potential composition mismatch')
  call eval_periodic_nonzero_panel_reference(mesh, target, 1.0_dp, 1.0_dp, 12_i32, 16_i32, &
                                             reference_potential, reference_field)
  charge_scale = sum(abs(mesh%q_elem))/eps0
  field_error = sqrt(sum((nonzero_field - reference_field)**2))/charge_scale
  potential_error = abs(nonzero_potential - reference_potential)/charge_scale
  write (*, '(a,2(es12.5,1x))') 'panel cached kneq0 errors(field,potential)=', field_error, potential_error
  write (*, '(a,8(es12.5,1x))') 'panel cached/ref values=', nonzero_field, reference_field, &
    nonzero_potential, reference_potential
  call assert_true(field_error < 1.0e-1_dp, 'panel cached kneq0 field exceeds the charge-scale error contract')
  call assert_true(potential_error < 1.0e-1_dp, 'panel cached kneq0 potential exceeds the charge-scale error contract')
  call test_end()

  call test_kneq0_potential_step_measurement()

  call test_adaptive_kneq0_rejects_then_accepts_half_step()

  call test_public_c_abi_acceptance()

  call test_begin('cached_snapshot_refreshes_both_components')
  mesh%q_elem = [1.0e-12_dp, -1.0e-12_dp]
  call snapshot%refresh(mesh)
  call evaluate_components(expected_field, expected_potential, nonzero_field, nonzero_potential, zero_field, zero_potential)
  call snapshot%eval_local_e(mesh, target, total_field)
  call snapshot%eval_local_phi(mesh, sim, target, total_potential)
  call assert_allclose_1d(total_field, expected_field, 1.0e-10_dp, 'refreshed snapshot field composition mismatch')
  call assert_close_dp(total_potential, expected_potential, 1.0e-10_dp, 'refreshed snapshot potential composition mismatch')
  call eval_periodic_nonzero_panel_reference(mesh, target, 1.0_dp, 1.0_dp, 12_i32, 16_i32, &
                                             reference_potential, reference_field)
  charge_scale = sum(abs(mesh%q_elem))/eps0
  field_error = sqrt(sum((nonzero_field - reference_field)**2))/charge_scale
  potential_error = abs(nonzero_potential - reference_potential)/charge_scale
  call assert_true(field_error < 1.0e-1_dp, 'neutral panel field exceeds the charge-scale error contract')
  call assert_true(potential_error < 1.0e-1_dp, 'neutral panel potential exceeds the charge-scale error contract')
  call test_end()

  call delete_file_if_exists(cache_path)
  call delete_file_if_exists(trim(cache_path)//'.lock')
  call remove_empty_directory(cache_dir)
  call test_summary()

contains

  subroutine test_kneq0_potential_step_measurement()
    real(dp) :: current_charge(2), candidate_charge(2)
    real(dp) :: measured_step(2), expected_step(2)
    real(dp) :: current_potential(2), candidate_potential(2)
    real(dp) :: field_before(3), field_after(3), potential_before, potential_after
    real(dp) :: max_abs_delta_phi
    integer(i32) :: update_count_before

    call test_begin('cached_snapshot_measures_kneq0_candidate_step_without_changing_state')
    current_charge = mesh%q_elem
    candidate_charge = current_charge + [0.7e-12_dp, -0.9e-12_dp]
    call snapshot%nonzero_solver%compute_mesh_potential(mesh, sim, current_potential)
    call snapshot%nonzero_solver%eval_e(mesh, target, field_before)
    call snapshot%nonzero_solver%eval_potential(mesh, sim, target, potential_before)
    update_count_before = snapshot%nonzero_solver%fmm_core_state%update_count

    call snapshot%measure_kneq0_potential_step( &
      mesh, candidate_charge, max_abs_delta_phi, measured_step &
      )

    call snapshot%nonzero_solver%eval_e(mesh, target, field_after)
    call snapshot%nonzero_solver%eval_potential(mesh, sim, target, potential_after)
    call assert_equal_i32( &
      snapshot%nonzero_solver%fmm_core_state%update_count, update_count_before, &
      'kneq0 potential-step measurement must preserve the solver update counter' &
      )
    call assert_true( &
      all(snapshot%nonzero_solver%fmm_core_state%src_q == current_charge), &
      'kneq0 potential-step measurement must restore the current mesh charges' &
      )
    call assert_allclose_1d( &
      field_after, field_before, 1.0e-13_dp, &
      'kneq0 potential-step measurement changed the cached field state' &
      )
    call assert_close_dp( &
      potential_after, potential_before, 1.0e-13_dp, &
      'kneq0 potential-step measurement changed the cached potential state' &
      )
    call assert_close_dp( &
      max_abs_delta_phi, maxval(abs(measured_step)), 1.0e-14_dp, &
      'kneq0 potential-step maximum does not match its vector diagnostic' &
      )

    mesh%q_elem = candidate_charge
    call snapshot%refresh(mesh)
    call snapshot%nonzero_solver%compute_mesh_potential(mesh, sim, candidate_potential)
    expected_step = candidate_potential - current_potential
    mesh%q_elem = current_charge
    call snapshot%refresh(mesh)

    call assert_allclose_1d( &
      measured_step, expected_step, 1.0e-10_dp, &
      'candidate-current potential difference does not match the cached kneq0 linear step' &
      )
    call test_end()
  end subroutine test_kneq0_potential_step_measurement

  subroutine test_adaptive_kneq0_rejects_then_accepts_half_step()
    type(mesh_type) :: adaptive_mesh
    type(app_config) :: adaptive_cfg
    type(sim_stats) :: adaptive_stats, adaptive_resume_stats
    type(charge_ledger_type) :: adaptive_ledger
    type(electrostatic_snapshot_type) :: step_probe
    real(dp) :: full_trial_charge(2)
    real(dp) :: full_potential_step, potential_limit
    real(dp), parameter :: maximum_duration = 1.0e-3_dp
    real(dp), parameter :: opening_low = 0.20_dp
    real(dp), parameter :: opening_high = 0.25_dp
    real(dp), parameter :: emit_current_density = 1.0e-6_dp
    real(dp), parameter :: electron_mass = 9.1093837139e-31_dp
    integer(i32) :: actual_team_size
    logical :: dynamic_before, dynamic_after
!$  integer(kind=omp_sched_kind) :: schedule_before_kind, schedule_after_kind
!$  integer :: schedule_before_chunk, schedule_after_chunk

    call test_begin('adaptive_kneq0_rejects_full_step_and_commits_only_half_step')
    actual_team_size = 1_i32
    dynamic_before = .false.
    dynamic_after = .true.
!$  dynamic_before = omp_get_dynamic()
!$  call omp_get_schedule(schedule_before_kind, schedule_before_chunk)
!$  call omp_set_dynamic(.false.)
    !$omp parallel default(none) shared(actual_team_size)
    !$omp single
!$  actual_team_size = int(omp_get_num_threads(), i32)
    !$omp end single
    !$omp end parallel
!$  call omp_set_dynamic(.true.)
!$  call omp_set_schedule(omp_sched_dynamic, 3)
    call init_mesh(adaptive_mesh, v0, v1, v2)
    adaptive_mesh%elem_vacuum_sign = 1_i32
    adaptive_mesh%vacuum_normals = adaptive_mesh%normals

    call default_app_config(adaptive_cfg)
    adaptive_cfg%mesh_mode = 'obj'
    adaptive_cfg%sim = sim
    adaptive_cfg%sim%rng_seed = 1207_i32
    adaptive_cfg%sim%batch_count = 1_i32
    adaptive_cfg%sim%dt = 1.0e-5_dp
    adaptive_cfg%sim%batch_duration = maximum_duration
    adaptive_cfg%sim%has_batch_duration = .true.
    adaptive_cfg%sim%max_step = 1_i32
    adaptive_cfg%sim%q_floor = 1.0e-40_dp
    adaptive_cfg%sim%e0 = 0.0_dp
    adaptive_cfg%field = field_config
    adaptive_cfg%periodic2 = periodic_config
    adaptive_cfg%panel = panel_config
    adaptive_cfg%n_particle_species = 1_i32
    adaptive_cfg%particle_species(1) = species_from_defaults()
    adaptive_cfg%particle_species(1)%species_key = 'photoelectron'
    adaptive_cfg%particle_species(1)%source_mode = 'photo_raycast'
    adaptive_cfg%particle_species(1)%inject_face = 'z_high'
    adaptive_cfg%particle_species(1)%q_particle = -qe
    adaptive_cfg%particle_species(1)%m_particle = electron_mass
    adaptive_cfg%particle_species(1)%temperature_ev = 0.0_dp
    adaptive_cfg%particle_species(1)%has_temperature_ev = .true.
    adaptive_cfg%particle_species(1)%normal_drift_speed = 1.0_dp
    adaptive_cfg%particle_species(1)%emit_current_density_a_m2 = emit_current_density
    adaptive_cfg%particle_species(1)%rays_per_batch = 1_i32
    adaptive_cfg%particle_species(1)%deposit_opposite_charge_on_emit = .true.
    adaptive_cfg%particle_species(1)%has_deposit_opposite_charge_on_emit = .true.
    ! This opening lies wholly over element 1, so one ray deposits a known
    ! charge proportional to batch_duration while the emitted particle survives.
    adaptive_cfg%particle_species(1)%pos_low = [opening_low, opening_low, sim%box_max(3)]
    adaptive_cfg%particle_species(1)%pos_high = [opening_high, opening_high, sim%box_max(3)]
    adaptive_cfg%particle_species(1)%ray_direction = [0.0_dp, 0.0_dp, -1.0_dp]
    adaptive_cfg%particle_species(1)%has_ray_direction = .true.

    call prepare_periodic2_collision_mesh(adaptive_mesh, adaptive_cfg%sim)
    call step_probe%init( &
      adaptive_mesh, adaptive_cfg%sim, adaptive_cfg%field, adaptive_cfg%periodic2, adaptive_cfg%panel &
      )
    call step_probe%refresh(adaptive_mesh)
    full_trial_charge = [ &
                        emit_current_density*(opening_high - opening_low)**2*maximum_duration, &
                        0.0_dp &
                        ]
    call step_probe%measure_kneq0_potential_step( &
      adaptive_mesh, full_trial_charge, full_potential_step &
      )
    call assert_true(full_potential_step > 0.0_dp, &
                     'adaptive calibration must produce a nonzero local potential step')
    potential_limit = 0.75_dp*full_potential_step

    adaptive_cfg%periodic2%max_nonzero_mode_potential_step = potential_limit
    call seed_particles_from_config(adaptive_cfg)
    call run_absorption_insulator( &
      adaptive_mesh, adaptive_cfg, adaptive_stats, charge_ledger=adaptive_ledger &
      )
!$  dynamic_after = omp_get_dynamic()
    call assert_true(dynamic_after, 'adaptive run did not restore the caller OpenMP dynamic ICV')
!$  call omp_get_schedule(schedule_after_kind, schedule_after_chunk)
!$  call assert_equal_i32( &
!$    int(schedule_after_kind, i32), int(omp_sched_dynamic, i32), &
!$    'adaptive run did not restore the caller OpenMP schedule kind' &
!$    )
!$  call assert_equal_i32( &
!$    int(schedule_after_chunk, i32), 3_i32, &
!$    'adaptive run did not restore the caller OpenMP schedule chunk' &
!$    )

    call assert_equal_i64(adaptive_stats%adaptive_nonzero_mode_rejected_trials, 1_i64, &
                          'adaptive batch must reject the full-duration trial exactly once')
    call assert_equal_i32(adaptive_stats%batches, 1_i32, &
                          'rejected adaptive trials must not increment the accepted batch count')
    call assert_equal_i32( &
      adaptive_stats%adaptive_nonzero_mode_omp_threads, actual_team_size, &
      'adaptive run must checkpoint its actual OpenMP team size' &
      )
    call assert_equal_i64(adaptive_stats%processed_particles, 1_i64, &
                          'rejected adaptive trials must not enter particle statistics')
    call assert_equal_i64(adaptive_ledger%emitted_count(1), 1_i64, &
                          'rejected adaptive trials must not enter the charge ledger')
    call assert_close_dp( &
      adaptive_stats%simulated_time, 0.5_dp*maximum_duration, &
      64.0_dp*epsilon(1.0_dp)*maximum_duration, &
      'simulated time must include only the accepted half-duration trial' &
      )
    call assert_close_dp( &
      adaptive_stats%adaptive_nonzero_mode_last_batch_duration, 0.5_dp*maximum_duration, &
      64.0_dp*epsilon(1.0_dp)*maximum_duration, &
      'adaptive accepted duration must be one half of the rejected trial' &
      )
    call assert_true( &
      adaptive_stats%adaptive_nonzero_mode_last_potential_step <= potential_limit, &
      'accepted adaptive kneq0 potential step exceeds the configured limit' &
      )
    call assert_close_dp( &
      adaptive_stats%adaptive_nonzero_mode_last_potential_step, 0.5_dp*full_potential_step, &
      1.0e-10_dp*max(1.0_dp, full_potential_step), &
      'accepted adaptive potential step does not scale with the half duration' &
      )
    call assert_allclose_1d( &
      adaptive_mesh%q_elem, 0.5_dp*full_trial_charge, &
      1.0e-12_dp*max(maxval(abs(full_trial_charge)), tiny(1.0_dp)), &
      'rejected full-duration charge leaked into the committed mesh state' &
      )
    call run_absorption_insulator( &
      adaptive_mesh, adaptive_cfg, adaptive_resume_stats, initial_stats=adaptive_stats &
      )
    call assert_equal_i32(adaptive_resume_stats%batches, adaptive_stats%batches, &
                          'same-team adaptive zero-batch resume changed the accepted batch count')
    call assert_equal_i32( &
      adaptive_resume_stats%adaptive_nonzero_mode_omp_threads, &
      adaptive_stats%adaptive_nonzero_mode_omp_threads, &
      'same-team adaptive resume changed its checkpointed OpenMP team size' &
      )
    call assert_close_dp( &
      adaptive_resume_stats%simulated_time, adaptive_stats%simulated_time, 0.0_dp, &
      'same-team adaptive resume changed its accepted physical time' &
      )
    call assert_allclose_1d( &
      adaptive_mesh%q_elem, 0.5_dp*full_trial_charge, &
      1.0e-12_dp*max(maxval(abs(full_trial_charge)), tiny(1.0_dp)), &
      'same-team adaptive zero-batch resume changed the mesh charge' &
      )
    call test_adaptive_team_mismatch_context()
!$  call omp_set_schedule(schedule_before_kind, schedule_before_chunk)
!$  call omp_set_dynamic(dynamic_before)
    call test_end()
  end subroutine test_adaptive_kneq0_rejects_then_accepts_half_step

  subroutine test_adaptive_team_mismatch_context()
    character(len=1024) :: executable_path, command, child_line
    integer :: child_exit_status, child_cmd_status, child_unit, child_ios
    logical :: saw_mismatch, saw_checkpoint_size, saw_current_size

    call get_command_argument(0, executable_path)
    call delete_file_if_exists(team_mismatch_path)
    command = 'OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 "'//trim(executable_path)// &
              '" --adaptive-team-mismatch-probe > '//team_mismatch_path//' 2>&1'
    call execute_command_line(trim(command), wait=.true., exitstat=child_exit_status, cmdstat=child_cmd_status)
    call assert_equal_i32(int(child_cmd_status, i32), 0_i32, 'adaptive team-mismatch command status mismatch')
    call assert_true(child_exit_status /= 0, 'adaptive team-mismatch probe should terminate with nonzero status')

    saw_mismatch = .false.
    saw_checkpoint_size = .false.
    saw_current_size = .false.
    open (newunit=child_unit, file=team_mismatch_path, status='old', action='read', iostat=child_ios)
    if (child_ios /= 0) error stop 'failed to read adaptive team-mismatch probe output'
    do
      read (child_unit, '(A)', iostat=child_ios) child_line
      if (child_ios /= 0) exit
      saw_mismatch = saw_mismatch .or. index(child_line, 'team-size mismatch with checkpoint') > 0
      saw_checkpoint_size = saw_checkpoint_size .or. index(child_line, 'checkpoint=2') > 0
      saw_current_size = saw_current_size .or. index(child_line, 'current=1') > 0
    end do
    close (child_unit)
    call delete_file_if_exists(team_mismatch_path)
    call assert_true( &
      saw_mismatch .and. saw_checkpoint_size .and. saw_current_size, &
      'adaptive resume team-size mismatch diagnostic is incomplete' &
      )
  end subroutine test_adaptive_team_mismatch_context

  subroutine run_adaptive_team_mismatch_probe()
    type(mesh_type) :: probe_mesh
    type(app_config) :: probe_cfg
    type(sim_stats) :: checkpoint_stats, probe_stats
    real(dp) :: probe_v0(3, 1), probe_v1(3, 1), probe_v2(3, 1)

    probe_v0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    probe_v1(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
    probe_v2(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
    call init_mesh(probe_mesh, probe_v0, probe_v1, probe_v2)
    call default_app_config(probe_cfg)
    probe_cfg%sim%batch_count = 0_i32
    probe_cfg%periodic2%max_nonzero_mode_potential_step = 1.0_dp
    checkpoint_stats = sim_stats()
    checkpoint_stats%adaptive_nonzero_mode_omp_threads = 2_i32
    call run_absorption_insulator( &
      probe_mesh, probe_cfg, probe_stats, initial_stats=checkpoint_stats &
      )
  end subroutine run_adaptive_team_mismatch_probe

  subroutine evaluate_components(expected_e, expected_phi, nonzero_e, nonzero_phi, zero_e, zero_phi)
    real(dp), intent(out) :: expected_e(3), expected_phi, nonzero_e(3), nonzero_phi, zero_e, zero_phi
    call snapshot%nonzero_solver%eval_e(mesh, target, nonzero_e)
    call snapshot%nonzero_solver%eval_potential(mesh, sim, target, nonzero_phi)
    call eval_periodic_zero_mode(snapshot%zero_plan, snapshot%zero_state, target(3), zero_mode_trace_plus, zero_phi, zero_e)
    expected_e = nonzero_e + sim%e0
    expected_e(3) = expected_e(3) + zero_e
    expected_phi = nonzero_phi + zero_phi - dot_product(sim%e0, target)
  end subroutine evaluate_components

  subroutine test_public_c_abi_acceptance()
    integer, parameter :: ntarget = 6
    type(c_ptr) :: kernel_handle, zero_handle
    integer(c_int) :: status
    integer(c_int), target :: periodic_axes(2), cache_hit, cache_build_count
    integer(c_int), target :: fingerprint_length, path_length
    real(c_double), target :: panel_v0(3, 2), panel_v1(3, 2), panel_v2(3, 2)
    real(c_double), target :: charge(2), source_heights(3, 2)
    real(c_double), target :: periodic_length(2), box_min(3), box_max(3)
    real(c_double), target :: target_points(3, ntarget), target_z(ntarget)
    real(c_double), target :: nonzero_e(3, ntarget), nonzero_phi(ntarget)
    real(c_double), target :: zero_e(ntarget), zero_phi(ntarget)
    character(kind=c_char), allocatable, target :: cache_bytes(:)
    character(kind=c_char), target :: fingerprint_buffer(129), path_buffer(513)
    real(dp) :: snapshot_e(3), expected_e(3), snapshot_phi, expected_phi
    character(len=160) :: message
    integer :: point

    call test_begin('cached_snapshot_matches_public_c_abis_at_and_off_surface')

    panel_v0 = real(mesh%v0, c_double)
    panel_v1 = real(mesh%v1, c_double)
    panel_v2 = real(mesh%v2, c_double)
    charge = real(mesh%q_elem, c_double)
    source_heights(1, :) = panel_v0(3, :)
    source_heights(2, :) = panel_v1(3, :)
    source_heights(3, :) = panel_v2(3, :)
    periodic_axes = int(snapshot%nonzero_solver%fmm_core_options%periodic_axes, c_int)
    periodic_length = real(snapshot%nonzero_solver%fmm_core_options%periodic_len, c_double)
    box_min = real(snapshot%nonzero_solver%fmm_core_options%target_box_min, c_double)
    box_max = real(snapshot%nonzero_solver%fmm_core_options%target_box_max, c_double)
    allocate (cache_bytes(len_trim(cache_dir)))
    call set_c_text(cache_bytes, trim(cache_dir))

    status = beach_kernel_create(kernel_handle)
    call assert_c_ok(status, 'FieldKernel create status')
    status = beach_kernel_set_periodic_cache( &
             kernel_handle, c_loc(cache_bytes), int(size(cache_bytes), c_int), &
             real(sim%field_periodic_generation_tolerance, c_double) &
             )
    call assert_c_ok(status, 'FieldKernel cache status')
    status = beach_kernel_build( &
             kernel_handle, int(mesh%nelem, c_int), c_loc(panel_v0), c_loc(panel_v1), c_loc(panel_v2), &
             real(snapshot%nonzero_solver%fmm_core_options%theta, c_double), &
             int(snapshot%nonzero_solver%fmm_core_options%leaf_max, c_int), &
             int(snapshot%nonzero_solver%fmm_core_options%order, c_int), &
             1_c_int, c_loc(periodic_axes), c_loc(periodic_length), &
             int(snapshot%nonzero_solver%fmm_core_options%periodic_image_layers, c_int), 3_c_int, &
             real(snapshot%nonzero_solver%fmm_core_options%periodic_ewald_alpha, c_double), &
             int(snapshot%nonzero_solver%fmm_core_options%periodic_ewald_layers, c_int), &
             c_loc(box_min), c_loc(box_max) &
             )
    call assert_c_ok(status, 'FieldKernel cached panel build status')
    status = beach_kernel_get_periodic_cache_info( &
             kernel_handle, c_loc(cache_hit), c_loc(cache_build_count), c_loc(fingerprint_buffer), &
             int(size(fingerprint_buffer), c_int), c_loc(fingerprint_length), c_loc(path_buffer), &
             int(size(path_buffer), c_int), c_loc(path_length) &
             )
    call assert_c_ok(status, 'FieldKernel cache info status')
    call assert_equal_i32(int(cache_hit, i32), 1_i32, 'FieldKernel must reuse the snapshot cache')
    call assert_equal_i32(int(cache_build_count, i32), 0_i32, 'FieldKernel must not rebuild the snapshot cache')
    call assert_true( &
      trim(c_buffer_text(fingerprint_buffer, fingerprint_length)) == &
      trim(snapshot%diagnostics%periodic_cache_fingerprint), &
      'FieldKernel and snapshot cache fingerprints must match' &
      )
    call assert_true( &
      trim(c_buffer_text(path_buffer, path_length)) == trim(snapshot%diagnostics%periodic_cache_path), &
      'FieldKernel and snapshot cache paths must match' &
      )
    status = beach_kernel_update_charges(kernel_handle, int(mesh%nelem, c_int), c_loc(charge))
    call assert_c_ok(status, 'FieldKernel charge update status')

    status = beach_zero_mode_create(zero_handle)
    call assert_c_ok(status, 'zero-mode create status')
    status = beach_zero_mode_build( &
             zero_handle, int(mesh%nelem, c_int), c_loc(source_heights), &
             real(product(sim%box_max(1:2) - sim%box_min(1:2)), c_double) &
             )
    call assert_c_ok(status, 'zero-mode build status')
    status = beach_zero_mode_update( &
             zero_handle, int(mesh%nelem, c_int), c_loc(charge), real(snapshot%zero_state%e_bottom, c_double), &
             minval(source_heights), 0.0_c_double &
             )
    call assert_c_ok(status, 'zero-mode charge update status')

    ! The first two points share a source height but lie outside the corresponding panel.
    target_points(:, 1) = [0.85_c_double, 0.15_c_double, 0.25_c_double]
    target_points(:, 2) = [0.15_c_double, 0.85_c_double, 0.65_c_double]
    target_points(:, 3) = [0.37_c_double, 0.61_c_double, 0.42_c_double]
    target_points(:, 4) = [0.19_c_double, 0.73_c_double, 0.88_c_double]
    target_points(:, 5) = [0.25_c_double, 0.25_c_double, 0.25_c_double]
    target_points(:, 6) = [0.75_c_double, 0.60_c_double, 0.65_c_double]
    target_z = target_points(3, :)
    status = beach_kernel_eval_e( &
             kernel_handle, int(ntarget, c_int), c_loc(target_points), c_loc(nonzero_e) &
             )
    call assert_c_ok(status, 'FieldKernel field evaluation status')
    status = beach_kernel_eval_phi( &
             kernel_handle, int(ntarget, c_int), c_loc(target_points), c_loc(nonzero_phi) &
             )
    call assert_c_ok(status, 'FieldKernel potential evaluation status')
    status = beach_zero_mode_eval( &
             zero_handle, int(ntarget, c_int), c_loc(target_z), int(zero_mode_trace_plus, c_int), &
             c_loc(zero_phi), c_loc(zero_e) &
             )
    call assert_c_ok(status, 'zero-mode plus-trace evaluation status')

    do point = 1, ntarget
      call snapshot%eval_local_e(mesh, real(target_points(:, point), dp), snapshot_e)
      call snapshot%eval_local_phi(mesh, sim, real(target_points(:, point), dp), snapshot_phi)
      expected_e = real(nonzero_e(:, point), dp) + sim%e0
      expected_e(3) = expected_e(3) + real(zero_e(point), dp)
      expected_phi = real(nonzero_phi(point) + zero_phi(point), dp) - &
                     dot_product(sim%e0, real(target_points(:, point), dp))
      write (message, '(a,i0)') 'snapshot/C-ABI field mismatch at point ', point
      call assert_allclose_1d(snapshot_e, expected_e, 1.0e-10_dp, trim(message))
      write (message, '(a,i0)') 'snapshot/C-ABI potential mismatch at point ', point
      call assert_close_dp(snapshot_phi, expected_phi, 1.0e-10_dp, trim(message))
    end do

    status = beach_zero_mode_destroy(zero_handle)
    call assert_c_ok(status, 'zero-mode destroy status')
    status = beach_kernel_destroy(kernel_handle)
    call assert_c_ok(status, 'FieldKernel destroy status')
    call test_end()
  end subroutine test_public_c_abi_acceptance

  subroutine assert_c_ok(actual, message)
    integer(c_int), intent(in) :: actual
    character(len=*), intent(in) :: message

    call assert_equal_i32(int(actual, i32), int(beach_kernel_ok, i32), message)
  end subroutine assert_c_ok

  subroutine set_c_text(output, value)
    character(kind=c_char), intent(out) :: output(:)
    character(len=*), intent(in) :: value
    integer :: i

    do i = 1, size(output)
      output(i) = achar(iachar(value(i:i)), kind=c_char)
    end do
  end subroutine set_c_text

  function c_buffer_text(buffer, length) result(value)
    character(kind=c_char), intent(in) :: buffer(:)
    integer(c_int), intent(in) :: length
    character(len=size(buffer)) :: value
    integer :: i

    value = ''
    do i = 1, int(length)
      value(i:i) = achar(iachar(buffer(i)))
    end do
  end function c_buffer_text

  subroutine configure_fixture(mesh_out, sim_out, field_out, periodic_out, panel_out, a, b, c)
    type(mesh_type), intent(out) :: mesh_out
    type(sim_config), intent(out) :: sim_out
    type(field_physics_config), intent(out) :: field_out
    type(periodic2_physics_config), intent(out) :: periodic_out
    type(panel_kernel_config), intent(out) :: panel_out
    real(dp), intent(out) :: a(3, 2), b(3, 2), c(3, 2)

    a(:, 1) = [0.10_dp, 0.10_dp, 0.25_dp]
    b(:, 1) = [0.55_dp, 0.10_dp, 0.25_dp]
    c(:, 1) = [0.10_dp, 0.55_dp, 0.25_dp]
    a(:, 2) = [0.45_dp, 0.45_dp, 0.65_dp]
    b(:, 2) = [0.90_dp, 0.45_dp, 0.65_dp]
    c(:, 2) = [0.90_dp, 0.90_dp, 0.65_dp]
    call init_mesh(mesh_out, a, b, c, q0=[1.0e-12_dp, 2.0e-12_dp])
    mesh_out%elem_vacuum_sign = 1_i32
    mesh_out%vacuum_normals = mesh_out%normals

    sim_out = sim_config()
    sim_out%field_solver = 'fmm'
    sim_out%field_bc_mode = 'periodic2'
    sim_out%field_periodic_far_correction = 'cached_kneq0'
    sim_out%field_periodic_image_layers = 1_i32
    sim_out%field_periodic_ewald_layers = 3_i32
    sim_out%field_periodic_cache_dir = cache_dir
    sim_out%field_periodic_generation_tolerance = 1.0e-8_dp
    sim_out%field_normalization = 'si'
    sim_out%tree_theta = 0.5_dp
    sim_out%has_tree_theta = .true.
    sim_out%tree_leaf_max = 8_i32
    sim_out%has_tree_leaf_max = .true.
    sim_out%e0 = [0.2_dp, -0.1_dp, 0.3_dp]
    sim_out%use_box = .true.
    sim_out%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    sim_out%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    sim_out%bc_low = [bc_periodic, bc_periodic, bc_open]
    sim_out%bc_high = [bc_periodic, bc_periodic, bc_open]
    field_out = field_physics_config(backend='fmm', normalization='si')
    periodic_out = periodic2_physics_config( &
                   nonzero_mode_backend='cached_kneq0', zero_mode_policy='exclude_k0', &
                   lower_boundary_model='symmetric_vacuum' &
                   )
    panel_out = panel_kernel_config( &
                kernel_id='triangle_p0_exact_p2m_near', surface_side_policy='per_element' &
                )
  end subroutine configure_fixture

end program test_periodic2_cached_snapshot

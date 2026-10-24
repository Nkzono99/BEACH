!> cached_kneq0 の adaptive batch: 非零 mode の電位変化が上限を超える試行を棄却し、半分の batch だけを確定する。
program test_periodic2_adaptive_batch
!$ use omp_lib
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: eps0, qe
  use bem_types, only: mesh_type, sim_config, sim_stats
  use bem_mesh, only: init_mesh, prepare_periodic2_collision_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config
  use bem_charge_ledger, only: charge_ledger_type
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use test_periodic2_cached_fixture, only: configure_cached_periodic2_fixture
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp, &
                          assert_allclose_1d, assert_equal_i32, assert_equal_i64, unique_test_path, &
                          delete_periodic_cache
  implicit none

  type(mesh_type) :: mesh
  type(sim_config) :: sim
  type(field_physics_config) :: field_config
  type(periodic2_physics_config) :: periodic_config
  type(panel_kernel_config) :: panel_config
  real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2)
  character(len=512) :: cache_path
  character(len=:), allocatable :: cache_dir

  cache_dir = unique_test_path('test_periodic2_adaptive_batch_tmp_')
  call configure_cached_periodic2_fixture(mesh, sim, field_config, periodic_config, panel_config, v0, v1, v2, cache_dir)
  call test_init(1)
  call test_adaptive_kneq0_rejects_then_accepts_half_step()
  call delete_periodic_cache(cache_path, cache_dir)
  call test_summary()

contains

  subroutine test_adaptive_kneq0_rejects_then_accepts_half_step()
    type(mesh_type) :: adaptive_mesh
    type(app_config) :: adaptive_cfg
    type(sim_stats) :: adaptive_stats, adaptive_resume_stats
    type(charge_ledger_type) :: adaptive_ledger
    type(electrostatic_snapshot_type) :: step_probe
    real(dp) :: full_trial_charge(2)
    real(dp) :: full_potential_step, potential_limit
    real(dp), parameter :: maximum_duration = 1.0e-3_dp
    real(dp), parameter :: opening_low = 0.13_dp
    real(dp), parameter :: opening_high = 0.16_dp
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
    cache_path = step_probe%diagnostics%periodic_cache_path
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

    call assert_equal_i64(adaptive_stats%adaptive_nonzero_mode%rejected_trials, 1_i64, &
                          'adaptive batch must reject the full-duration trial exactly once')
    call assert_equal_i32(adaptive_stats%batches, 1_i32, &
                          'rejected adaptive trials must not increment the accepted batch count')
    call assert_equal_i32( &
      adaptive_stats%adaptive_nonzero_mode%omp_threads, actual_team_size, &
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
      adaptive_stats%adaptive_nonzero_mode%last_batch_duration, 0.5_dp*maximum_duration, &
      64.0_dp*epsilon(1.0_dp)*maximum_duration, &
      'adaptive accepted duration must be one half of the rejected trial' &
      )
    call assert_true( &
      adaptive_stats%adaptive_nonzero_mode%last_potential_step <= potential_limit, &
      'accepted adaptive kneq0 potential step exceeds the configured limit' &
      )
    call assert_close_dp( &
      adaptive_stats%adaptive_nonzero_mode%last_potential_step, 0.5_dp*full_potential_step, &
      1.0e-10_dp*max(1.0_dp, full_potential_step), &
      'accepted adaptive potential step does not scale with the half duration' &
      )
    call assert_allclose_1d( &
      adaptive_mesh%q_elem, 0.5_dp*full_trial_charge, &
      1.0e-12_dp*max(maxval(abs(full_trial_charge)), tiny(1.0_dp)), &
      'rejected full-duration charge leaked into the committed mesh state' &
      )
    call assert_close_dp( &
      adaptive_stats%matching_plane%top_cell_mean_v_m, &
      sum(adaptive_mesh%q_elem)/(2.0_dp*eps0*product(sim%box_max(1:2) - sim%box_min(1:2))), &
      1.0e-12_dp, 'top-cell mean field receipt must use the committed charge' &
      )
    adaptive_stats%adaptive_nonzero_mode%omp_threads = actual_team_size + 1_i32
    call run_absorption_insulator( &
      adaptive_mesh, adaptive_cfg, adaptive_resume_stats, initial_stats=adaptive_stats &
      )
    call assert_equal_i32(adaptive_resume_stats%batches, adaptive_stats%batches, &
                          'same-team adaptive zero-batch resume changed the accepted batch count')
    call assert_equal_i32( &
      adaptive_resume_stats%adaptive_nonzero_mode%omp_threads, &
      actual_team_size, &
      'adaptive resume must record the current team size rather than require the checkpoint value' &
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
!$  call omp_set_schedule(schedule_before_kind, schedule_before_chunk)
!$  call omp_set_dynamic(dynamic_before)
    call test_end()
  end subroutine test_adaptive_kneq0_rejects_then_accepts_half_step

end program test_periodic2_adaptive_batch

!> 粒子 step・光電子 ray・統計の失敗が、batch・粒子・状態を含む診断文で停止することを子 process で検証する。
program test_simulator_failure_context
  use bem_kinds, only: dp, i32, i64
  use bem_mesh, only: init_mesh, prepare_periodic2_collision_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config
  use bem_physics_config_types, only: normalize_legacy_physics_config
  use bem_types, only: mesh_type, sim_stats, bc_open, bc_reflect, bc_periodic
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, delete_file_if_exists, &
                          run_self_probe, file_contains
  implicit none

  character(len=*), parameter :: collision_failure_path = 'test_simulator_failure_collision_tmp.log'
  character(len=*), parameter :: photo_collision_failure_path = 'test_simulator_failure_photo_collision_tmp.log'
  character(len=*), parameter :: box_event_failure_path = 'test_simulator_failure_box_event_tmp.log'
  character(len=*), parameter :: invalid_candidate_failure_path = 'test_simulator_failure_invalid_candidate_tmp.log'
  character(len=*), parameter :: statistics_overflow_path = 'test_simulator_failure_statistics_overflow_tmp.log'
  character(len=64) :: run_mode

  call get_command_argument(1, run_mode)
  if (trim(run_mode) == '--collision-query-failure-probe') then
    call run_collision_query_failure_probe()
    error stop 'collision query failure probe unexpectedly completed'
  else if (trim(run_mode) == '--photo-collision-query-failure-probe') then
    call run_photo_collision_query_failure_probe()
    error stop 'photo collision query failure probe unexpectedly completed'
  else if (trim(run_mode) == '--multiple-box-event-probe') then
    call run_multiple_box_event_failure_probe()
    error stop 'multiple box event failure probe unexpectedly completed'
  else if (trim(run_mode) == '--invalid-candidate-probe') then
    call run_invalid_candidate_failure_probe()
    error stop 'invalid candidate failure probe unexpectedly completed'
  else if (trim(run_mode) == '--statistics-overflow-probe') then
    call run_statistics_overflow_probe()
    error stop 'statistics overflow probe unexpectedly completed'
  end if

  call test_init(5)

  call test_begin('collision_query_failure_context')
  call test_collision_query_failure_context()
  call test_end()

  call test_begin('photo_collision_query_failure_context')
  call test_photo_collision_query_failure_context()
  call test_end()

  call test_begin('multiple_box_event_failure_context')
  call test_multiple_box_event_failure_context()
  call test_end()

  call test_begin('invalid_candidate_failure_context')
  call test_invalid_candidate_failure_context()
  call test_end()

  call test_begin('statistics_overflow_fails_closed')
  call test_statistics_overflow_fails_closed()
  call test_end()

  call test_summary()

contains

  subroutine test_collision_query_failure_context()
    integer :: exit_status
    logical :: saw_survivor_warning

    call run_self_probe('--collision-query-failure-probe', collision_failure_path, exit_status, &
                        env='OMP_NUM_THREADS=2 BEACH_WARN_LONG_PARTICLE_STEPS=1')
    call assert_true(exit_status /= 0, 'collision failure probe should terminate with a nonzero status')
    call assert_true(file_contains(collision_failure_path, 'batch=1'), &
                     'collision failure message should include the batch index')
    call assert_true(file_contains(collision_failure_path, 'particle=1'), &
                     'collision failure message should include the lowest failing particle index')
    call assert_true(file_contains(collision_failure_path, 'status=multiple_box_events'), &
                     'particle-step failure message should include the selected status')
    saw_survivor_warning = file_contains(collision_failure_path, 'max-step-survivor')
    call delete_file_if_exists(collision_failure_path)
    call assert_true(.not. saw_survivor_warning, 'collision failure must not be reported as a max-step survivor')
  end subroutine test_collision_query_failure_context

  subroutine run_collision_query_failure_probe()
    type(mesh_type) :: failure_mesh
    type(app_config) :: failure_cfg
    type(sim_stats) :: failure_stats
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [0.1d0, 0.2d0, 0.0d0]
    tri_v1(:, 1) = [0.3d0, 0.2d0, 0.0d0]
    tri_v2(:, 1) = [0.1d0, 0.4d0, 0.0d0]
    call init_mesh(failure_mesh, tri_v0, tri_v1, tri_v2)
    failure_mesh%elem_vacuum_sign = 1_i32
    failure_mesh%vacuum_normals = failure_mesh%normals

    call default_app_config(failure_cfg)
    failure_cfg%sim%rng_seed = 991_i32
    failure_cfg%sim%batch_count = 1_i32
    failure_cfg%sim%dt = 1.0d0
    failure_cfg%sim%max_step = 1_i32
    failure_cfg%sim%field_solver = 'fmm'
    failure_cfg%sim%field_bc_mode = 'periodic2'
    failure_cfg%sim%field_periodic_far_correction = 'none'
    failure_cfg%sim%use_box = .true.
    failure_cfg%sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    failure_cfg%sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    failure_cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    failure_cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call normalize_legacy_physics_config( &
      failure_cfg%sim, failure_cfg%field, failure_cfg%periodic2, failure_cfg%panel &
      )

    failure_cfg%n_particle_species = 1_i32
    failure_cfg%particle_species(1) = species_from_defaults()
    failure_cfg%particle_species(1)%source_mode = 'volume_seed'
    failure_cfg%particle_species(1)%npcls_per_step = 2_i32
    failure_cfg%particle_species(1)%q_particle = 1.0d0
    failure_cfg%particle_species(1)%m_particle = 1.0d0
    failure_cfg%particle_species(1)%w_particle = 1.0d0
    failure_cfg%particle_species(1)%pos_low = [0.2d0, 0.25d0, 0.5d0]
    failure_cfg%particle_species(1)%pos_high = failure_cfg%particle_species(1)%pos_low
    ! Keep the failure probe beyond the periodic-subdivision recovery limit.
    failure_cfg%particle_species(1)%drift_velocity = [1.0d6, 0.0d0, 0.0d0]
    failure_cfg%particle_species(1)%temperature_k = 0.0d0

    call prepare_periodic2_collision_mesh(failure_mesh, failure_cfg%sim)
    call seed_particles_from_config(failure_cfg)
    call run_absorption_insulator(failure_mesh, failure_cfg, failure_stats)
  end subroutine run_collision_query_failure_probe

  subroutine test_photo_collision_query_failure_context()
    integer :: exit_status

    call run_self_probe('--photo-collision-query-failure-probe', photo_collision_failure_path, exit_status, &
                        env='OMP_NUM_THREADS=2')
    call assert_true(exit_status /= 0, 'photo collision failure probe should terminate with a nonzero status')
    call assert_true(file_contains(photo_collision_failure_path, 'batch=1'), &
                     'photo collision failure message should include the batch index')
    call assert_true(file_contains(photo_collision_failure_path, 'rank=0'), &
                     'photo collision failure message should include the selected rank')
    call assert_true(file_contains(photo_collision_failure_path, 'species=2'), &
                     'photo collision failure message should include the species index')
    call assert_true(file_contains(photo_collision_failure_path, 'ray=1'), &
                     'photo collision failure message should include the lowest failing ray')
    call assert_true(file_contains(photo_collision_failure_path, 'bounce=0'), &
                     'photo collision failure message should include the failing bounce')
    call assert_true(file_contains(photo_collision_failure_path, 'status=image_limit'), &
                     'photo collision failure message should include the collision status')
    call delete_file_if_exists(photo_collision_failure_path)
  end subroutine test_photo_collision_query_failure_context

  subroutine run_photo_collision_query_failure_probe()
    type(mesh_type) :: failure_mesh
    type(app_config) :: failure_cfg
    type(sim_stats) :: failure_stats
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [-2048.0d0, 0.0d0, 0.5d0]
    tri_v1(:, 1) = [2049.0d0, 0.0d0, 0.5d0]
    tri_v2(:, 1) = [0.5d0, 1.0d0, 0.5d0]
    call init_mesh(failure_mesh, tri_v0, tri_v1, tri_v2)
    failure_mesh%elem_vacuum_sign = 1_i32
    failure_mesh%vacuum_normals = failure_mesh%normals

    call default_app_config(failure_cfg)
    failure_cfg%sim%rng_seed = 992_i32
    failure_cfg%sim%batch_count = 1_i32
    failure_cfg%sim%batch_duration = 0.5d0
    failure_cfg%sim%dt = 1.0d0
    failure_cfg%sim%max_step = 1_i32
    failure_cfg%sim%field_solver = 'fmm'
    failure_cfg%sim%field_bc_mode = 'periodic2'
    failure_cfg%sim%field_periodic_far_correction = 'none'
    failure_cfg%sim%use_box = .true.
    failure_cfg%sim%box_min = [0.0d0, 0.0d0, 0.0d0]
    failure_cfg%sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    failure_cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    failure_cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call normalize_legacy_physics_config( &
      failure_cfg%sim, failure_cfg%field, failure_cfg%periodic2, failure_cfg%panel &
      )

    failure_cfg%n_particle_species = 2_i32
    failure_cfg%particle_species(1) = species_from_defaults()
    failure_cfg%particle_species(1)%source_mode = 'volume_seed'
    failure_cfg%particle_species(1)%npcls_per_step = 1_i32
    failure_cfg%particle_species(1)%q_particle = 1.0d0
    failure_cfg%particle_species(1)%m_particle = 1.0d0
    failure_cfg%particle_species(1)%w_particle = 1.0d0
    failure_cfg%particle_species(1)%temperature_k = 0.0d0
    failure_cfg%particle_species(1)%pos_low = [0.5d0, 0.5d0, 0.75d0]
    failure_cfg%particle_species(1)%pos_high = failure_cfg%particle_species(1)%pos_low
    failure_cfg%particle_species(1)%drift_velocity = 0.0d0
    failure_cfg%particle_species(2) = species_from_defaults()
    failure_cfg%particle_species(2)%source_mode = 'photo_raycast'
    failure_cfg%particle_species(2)%rays_per_batch = 4_i32
    failure_cfg%particle_species(2)%emit_current_density_a_m2 = 1.0d0
    failure_cfg%particle_species(2)%q_particle = -1.0d0
    failure_cfg%particle_species(2)%m_particle = 1.0d0
    failure_cfg%particle_species(2)%temperature_k = 0.0d0
    failure_cfg%particle_species(2)%inject_face = 'z_high'
    failure_cfg%particle_species(2)%pos_low = [0.0d0, 0.0d0, 1.0d0]
    failure_cfg%particle_species(2)%pos_high = [1.0d0, 1.0d0, 1.0d0]
    failure_cfg%particle_species(2)%ray_direction = [0.0d0, 0.0d0, -1.0d0]
    failure_cfg%particle_species(2)%has_ray_direction = .true.

    call prepare_periodic2_collision_mesh(failure_mesh, failure_cfg%sim)
    call seed_particles_from_config(failure_cfg)
    call run_absorption_insulator(failure_mesh, failure_cfg, failure_stats)
  end subroutine run_photo_collision_query_failure_probe

  subroutine test_multiple_box_event_failure_context()
    integer :: exit_status
    logical :: saw_particle, saw_step, saw_status, saw_dt, saw_x, saw_v

    call run_self_probe('--multiple-box-event-probe', box_event_failure_path, exit_status, env='OMP_NUM_THREADS=2')
    call assert_true(exit_status /= 0, 'multiple box event probe should terminate with nonzero status')
    saw_particle = file_contains(box_event_failure_path, 'particle=1')
    saw_step = file_contains(box_event_failure_path, 'step=1')
    saw_status = file_contains(box_event_failure_path, 'status=multiple_box_events')
    saw_dt = file_contains(box_event_failure_path, 'dt=')
    saw_x = file_contains(box_event_failure_path, 'x=')
    saw_v = file_contains(box_event_failure_path, 'v=')
    call delete_file_if_exists(box_event_failure_path)
    call assert_true(saw_particle .and. saw_step, 'box event failure should include particle and step')
    call assert_true(saw_status, 'box event failure should name multiple_box_events')
    call assert_true(saw_dt .and. saw_x .and. saw_v, 'box event failure should include dt, x, and v')
  end subroutine test_multiple_box_event_failure_context

  subroutine run_multiple_box_event_failure_probe()
    type(mesh_type) :: failure_mesh
    type(app_config) :: failure_cfg
    type(sim_stats) :: failure_stats
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [10.0_dp, -1.0_dp, -1.0_dp]
    tri_v1(:, 1) = [10.0_dp, 1.0_dp, -1.0_dp]
    tri_v2(:, 1) = [10.0_dp, 0.0_dp, 1.0_dp]
    call init_mesh(failure_mesh, tri_v0, tri_v1, tri_v2)
    failure_mesh%elem_vacuum_sign = 1_i32
    failure_mesh%vacuum_normals = failure_mesh%normals
    call default_app_config(failure_cfg)
    failure_cfg%sim%rng_seed = 994_i32
    failure_cfg%sim%batch_count = 1_i32
    failure_cfg%sim%dt = 1.0_dp
    failure_cfg%sim%max_step = 1_i32
    failure_cfg%sim%use_box = .true.
    failure_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    failure_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    failure_cfg%sim%bc_low(1) = bc_reflect
    failure_cfg%sim%bc_high(1) = bc_reflect
    failure_cfg%n_particle_species = 1_i32
    failure_cfg%particle_species(1) = species_from_defaults()
    failure_cfg%particle_species(1)%source_mode = 'volume_seed'
    failure_cfg%particle_species(1)%npcls_per_step = 1_i32
    failure_cfg%particle_species(1)%q_particle = 0.0_dp
    failure_cfg%particle_species(1)%m_particle = 1.0_dp
    failure_cfg%particle_species(1)%w_particle = 1.0_dp
    failure_cfg%particle_species(1)%pos_low = [0.9_dp, 0.2_dp, 0.2_dp]
    failure_cfg%particle_species(1)%pos_high = failure_cfg%particle_species(1)%pos_low
    failure_cfg%particle_species(1)%drift_velocity = [9.0_dp, 0.0_dp, 0.0_dp]
    failure_cfg%particle_species(1)%temperature_k = 0.0_dp
    call seed_particles_from_config(failure_cfg)
    call run_absorption_insulator(failure_mesh, failure_cfg, failure_stats)
  end subroutine run_multiple_box_event_failure_probe

  subroutine test_invalid_candidate_failure_context()
    integer :: exit_status
    logical :: saw_status, saw_survivor_warning

    call run_self_probe('--invalid-candidate-probe', invalid_candidate_failure_path, exit_status, env='OMP_NUM_THREADS=2')
    call assert_true(exit_status /= 0, 'invalid candidate probe should terminate with nonzero status')
    saw_status = file_contains(invalid_candidate_failure_path, 'status=invalid_boundary')
    saw_survivor_warning = file_contains(invalid_candidate_failure_path, 'max-step-survivor')
    call delete_file_if_exists(invalid_candidate_failure_path)
    call assert_true(saw_status, 'nonfinite candidate should be rejected as invalid_boundary')
    call assert_true(.not. saw_survivor_warning, 'invalid candidate must not be reported as a max-step survivor')
  end subroutine test_invalid_candidate_failure_context

  subroutine run_invalid_candidate_failure_probe()
    type(mesh_type) :: failure_mesh
    type(app_config) :: failure_cfg
    type(sim_stats) :: failure_stats
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [10.0_dp, -1.0_dp, -1.0_dp]
    tri_v1(:, 1) = [10.0_dp, 1.0_dp, -1.0_dp]
    tri_v2(:, 1) = [10.0_dp, 0.0_dp, 1.0_dp]
    call init_mesh(failure_mesh, tri_v0, tri_v1, tri_v2)
    failure_mesh%elem_vacuum_sign = 1_i32
    failure_mesh%vacuum_normals = failure_mesh%normals
    call default_app_config(failure_cfg)
    failure_cfg%sim%rng_seed = 995_i32
    failure_cfg%sim%batch_count = 1_i32
    failure_cfg%sim%dt = 2.0_dp
    failure_cfg%sim%max_step = 1_i32
    failure_cfg%sim%use_box = .false.
    failure_cfg%n_particle_species = 1_i32
    failure_cfg%particle_species(1) = species_from_defaults()
    failure_cfg%particle_species(1)%source_mode = 'volume_seed'
    failure_cfg%particle_species(1)%npcls_per_step = 1_i32
    failure_cfg%particle_species(1)%q_particle = 0.0_dp
    failure_cfg%particle_species(1)%m_particle = 1.0_dp
    failure_cfg%particle_species(1)%w_particle = 1.0_dp
    failure_cfg%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, 0.0_dp]
    failure_cfg%particle_species(1)%pos_high = failure_cfg%particle_species(1)%pos_low
    failure_cfg%particle_species(1)%drift_velocity = [huge(1.0_dp), 0.0_dp, 0.0_dp]
    failure_cfg%particle_species(1)%temperature_k = 0.0_dp
    call seed_particles_from_config(failure_cfg)
    call run_absorption_insulator(failure_mesh, failure_cfg, failure_stats)
  end subroutine run_invalid_candidate_failure_probe

  subroutine test_statistics_overflow_fails_closed()
    integer :: exit_status
    logical :: saw_overflow, saw_field

    call run_self_probe('--statistics-overflow-probe', statistics_overflow_path, exit_status)
    call assert_true(exit_status /= 0, 'statistics overflow probe should terminate with nonzero status')
    saw_overflow = file_contains(statistics_overflow_path, 'simulation statistic overflow')
    saw_field = file_contains(statistics_overflow_path, 'processed_particles')
    call delete_file_if_exists(statistics_overflow_path)
    call assert_true(saw_overflow .and. saw_field, 'statistics overflow diagnostic is incomplete')
  end subroutine test_statistics_overflow_fails_closed

  !> One particle leaves an open box in one step; the seeded processed count cannot take the increment.
  subroutine run_statistics_overflow_probe()
    type(mesh_type) :: overflow_mesh
    type(app_config) :: overflow_cfg
    type(sim_stats) :: overflow_seed, overflow_stats
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [10.0_dp, -1.0_dp, -1.0_dp]
    tri_v1(:, 1) = [10.0_dp, 1.0_dp, -1.0_dp]
    tri_v2(:, 1) = [10.0_dp, 0.0_dp, 1.0_dp]
    call init_mesh(overflow_mesh, tri_v0, tri_v1, tri_v2)
    overflow_mesh%elem_vacuum_sign = 1_i32
    overflow_mesh%vacuum_normals = overflow_mesh%normals
    call default_app_config(overflow_cfg)
    overflow_cfg%sim%batch_count = 1_i32
    overflow_cfg%sim%dt = 1.0_dp
    overflow_cfg%sim%max_step = 1_i32
    overflow_cfg%sim%use_box = .true.
    overflow_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    overflow_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    overflow_cfg%n_particle_species = 1_i32
    overflow_cfg%particle_species(1) = species_from_defaults()
    overflow_cfg%particle_species(1)%source_mode = 'volume_seed'
    overflow_cfg%particle_species(1)%npcls_per_step = 1_i32
    overflow_cfg%particle_species(1)%q_particle = 0.0_dp
    overflow_cfg%particle_species(1)%m_particle = 1.0_dp
    overflow_cfg%particle_species(1)%w_particle = 1.0_dp
    overflow_cfg%particle_species(1)%pos_low = [0.5_dp, 0.5_dp, 0.5_dp]
    overflow_cfg%particle_species(1)%pos_high = overflow_cfg%particle_species(1)%pos_low
    overflow_cfg%particle_species(1)%drift_velocity = [2.0_dp, 0.0_dp, 0.0_dp]
    overflow_cfg%particle_species(1)%temperature_k = 0.0_dp
    call seed_particles_from_config(overflow_cfg)
    overflow_seed = sim_stats()
    overflow_seed%processed_particles = huge(0_i64)
    call run_absorption_insulator( &
      overflow_mesh, overflow_cfg, overflow_stats, initial_stats=overflow_seed &
      )
  end subroutine run_statistics_overflow_probe
end program test_simulator_failure_context

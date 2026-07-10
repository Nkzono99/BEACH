!> 吸着ループの統計更新・堆積・履歴出力を検証する統合テスト。
program test_simulator
  use bem_kinds, only: dp, i32, i64
  use bem_mesh, only: init_mesh, prepare_periodic2_collision_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config
  use bem_types, only: mesh_type, sim_stats, bc_open, bc_periodic
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_equal_i64, assert_close_dp, delete_file_if_exists
  implicit none

  type(mesh_type) :: mesh
  type(mesh_type) :: mesh_tree
  type(mesh_type) :: mesh_potential_history
  type(mesh_type) :: mesh_resume
  type(app_config) :: cfg, cfg_tree, cfg_potential_history
  type(sim_stats) :: stats, stats_tree, stats_potential_history, stats_seed, stats_resume
  real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1)
  real(dp) :: potential_value
  integer :: u, ios
  integer(i32) :: potential_batch_idx, potential_elem_idx
  character(len=256) :: line
  integer(i32) :: n_lines
  character(len=*), parameter :: history_path = 'test_simulator_history_tmp.csv'
  character(len=*), parameter :: potential_history_path = 'test_simulator_potential_history_tmp.csv'
  character(len=*), parameter :: collision_failure_path = 'test_simulator_collision_failure_tmp.log'
  character(len=*), parameter :: photo_collision_failure_path = 'test_simulator_photo_collision_failure_tmp.log'
  character(len=64) :: run_mode

  call get_command_argument(1, run_mode)
  if (trim(run_mode) == '--collision-query-failure-probe') then
    call run_collision_query_failure_probe()
    error stop 'collision query failure probe unexpectedly completed'
  else if (trim(run_mode) == '--photo-collision-query-failure-probe') then
    call run_photo_collision_query_failure_probe()
    error stop 'photo collision query failure probe unexpectedly completed'
  end if

  v0(:, 1) = [0.0d0, 0.0d0, 0.0d0]
  v1(:, 1) = [1.0d0, 0.0d0, 0.0d0]
  v2(:, 1) = [0.0d0, 1.0d0, 0.0d0]
  call init_mesh(mesh, v0, v1, v2)

  call default_app_config(cfg)
  cfg%sim%rng_seed = 777_i32
  cfg%sim%batch_count = 1_i32
  cfg%sim%dt = 1.0d0
  cfg%sim%max_step = 1_i32
  cfg%sim%softening = 1.0d-6
  cfg%sim%q_floor = 1.0d-30
  cfg%sim%use_box = .true.
  cfg%sim%box_min = [-1.0d0, -1.0d0, -2.0d0]
  cfg%sim%box_max = [1.0d0, 1.0d0, 1.0d0]

  cfg%n_particle_species = 4_i32

  cfg%particle_species(1) = species_from_defaults()
  cfg%particle_species(1)%source_mode = 'volume_seed'
  cfg%particle_species(1)%npcls_per_step = 1_i32
  cfg%particle_species(1)%q_particle = 1.0d0
  cfg%particle_species(1)%m_particle = 1.0d0
  cfg%particle_species(1)%w_particle = 2.0d0
  cfg%particle_species(1)%pos_low = [0.2d0, 0.2d0, 0.8d0]
  cfg%particle_species(1)%pos_high = [0.2d0, 0.2d0, 0.8d0]
  cfg%particle_species(1)%drift_velocity = [0.0d0, 0.0d0, -2.0d0]
  cfg%particle_species(1)%temperature_k = 0.0d0

  cfg%particle_species(2) = species_from_defaults()
  cfg%particle_species(2)%source_mode = 'volume_seed'
  cfg%particle_species(2)%npcls_per_step = 1_i32
  cfg%particle_species(2)%q_particle = 1.0d0
  cfg%particle_species(2)%m_particle = 1.0d0
  cfg%particle_species(2)%w_particle = 1.0d0
  cfg%particle_species(2)%pos_low = [0.9d0, 0.0d0, 0.0d0]
  cfg%particle_species(2)%pos_high = [0.9d0, 0.0d0, 0.0d0]
  cfg%particle_species(2)%drift_velocity = [2.0d0, 0.0d0, 0.0d0]
  cfg%particle_species(2)%temperature_k = 0.0d0

  cfg%particle_species(3) = species_from_defaults()
  cfg%particle_species(3)%source_mode = 'volume_seed'
  cfg%particle_species(3)%npcls_per_step = 1_i32
  cfg%particle_species(3)%q_particle = 1.0d0
  cfg%particle_species(3)%m_particle = 1.0d0
  cfg%particle_species(3)%w_particle = 1.0d0
  cfg%particle_species(3)%pos_low = [0.0d0, 0.0d0, 0.5d0]
  cfg%particle_species(3)%pos_high = [0.0d0, 0.0d0, 0.5d0]
  cfg%particle_species(3)%drift_velocity = [0.0d0, 1.0d0, 0.0d0]
  cfg%particle_species(3)%temperature_k = 0.0d0

  cfg%particle_species(4) = species_from_defaults()
  cfg%particle_species(4)%source_mode = 'volume_seed'
  cfg%particle_species(4)%npcls_per_step = 1_i32
  cfg%particle_species(4)%q_particle = 1.0d0
  cfg%particle_species(4)%m_particle = 1.0d0
  cfg%particle_species(4)%w_particle = 3.0d0
  cfg%particle_species(4)%pos_low = [0.3d0, 0.3d0, 0.4d0]
  cfg%particle_species(4)%pos_high = [0.3d0, 0.3d0, 0.4d0]
  cfg%particle_species(4)%drift_velocity = [0.0d0, 0.0d0, -5.0d0]
  cfg%particle_species(4)%temperature_k = 0.0d0

  call seed_particles_from_config(cfg)

  call test_init(7)

  call test_begin('basic_simulation')
  call delete_file_if_exists(history_path)
  open (newunit=u, file=history_path, status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'failed to open simulator history fixture'
  call run_absorption_insulator(mesh, cfg, stats, history_unit=u, history_stride=1_i32)
  close (u)

  call assert_equal_i64(stats%processed_particles, 4_i64, 'processed_particles mismatch')
  call assert_equal_i64(stats%absorbed, 2_i64, 'absorbed mismatch')
  call assert_equal_i64(stats%escaped, 2_i64, 'escaped mismatch')
  call assert_equal_i64(stats%escaped_boundary, 1_i64, 'escaped_boundary mismatch')
  call assert_equal_i64(stats%survived_max_step, 1_i64, 'survived_max_step mismatch')
  call assert_equal_i32(stats%batches, 1_i32, 'batch count mismatch')
  call assert_close_dp(mesh%q_elem(1), 5.0d0, 1.0d-12, 'deposited charge mismatch')
  call assert_true(stats%last_rel_change > 0.0d0, 'last_rel_change should be positive')
  call test_end()

  call test_begin('history_output')
  n_lines = 0_i32
  open (newunit=u, file=history_path, status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'failed to read simulator history fixture'
  do
    read (u, '(A)', iostat=ios) line
    if (ios /= 0) exit
    n_lines = n_lines + 1_i32
  end do
  close (u)
  call assert_equal_i32(n_lines, 1_i32, 'history snapshot line count mismatch')

  call delete_file_if_exists(history_path)
  call test_end()

  call test_begin('treecode_equivalence')
  call init_mesh(mesh_tree, v0, v1, v2)
  cfg_tree = cfg
  cfg_tree%sim%field_solver = 'treecode'
  call seed_particles_from_config(cfg_tree)
  call run_absorption_insulator(mesh_tree, cfg_tree, stats_tree)

  call assert_equal_i64(stats_tree%processed_particles, stats%processed_particles, 'treecode processed_particles mismatch')
  call assert_equal_i64(stats_tree%absorbed, stats%absorbed, 'treecode absorbed mismatch')
  call assert_equal_i64(stats_tree%escaped, stats%escaped, 'treecode escaped mismatch')
  call assert_equal_i64(stats_tree%escaped_boundary, stats%escaped_boundary, 'treecode escaped_boundary mismatch')
  call assert_equal_i64(stats_tree%survived_max_step, stats%survived_max_step, 'treecode survived_max_step mismatch')
  call assert_equal_i32(stats_tree%batches, stats%batches, 'treecode batches mismatch')
  call assert_close_dp(mesh_tree%q_elem(1), mesh%q_elem(1), 1.0d-12, 'treecode deposited charge mismatch')
  call test_end()

  call test_begin('fmm_potential_history_refreshes_after_commit')
  call init_mesh(mesh_potential_history, v0, v1, v2)
  cfg_potential_history = cfg
  cfg_potential_history%sim%field_solver = 'fmm'
  call seed_particles_from_config(cfg_potential_history)
  call delete_file_if_exists(potential_history_path)
  open (newunit=u, file=potential_history_path, status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'failed to open simulator potential history fixture'
  call run_absorption_insulator( &
    mesh_potential_history, cfg_potential_history, stats_potential_history, &
    potential_history_unit=u, history_stride=1_i32 &
    )
  close (u)
  open (newunit=u, file=potential_history_path, status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'failed to read simulator potential history fixture'
  read (u, *, iostat=ios) potential_batch_idx, potential_elem_idx, potential_value
  close (u)
  if (ios /= 0) error stop 'failed to parse simulator potential history fixture'
  call assert_equal_i32(potential_batch_idx, 1_i32, 'potential history batch mismatch')
  call assert_equal_i32(potential_elem_idx, 1_i32, 'potential history elem mismatch')
  call assert_true(potential_value > 0.0d0, 'FMM potential history should use post-commit charges')
  call delete_file_if_exists(potential_history_path)
  call test_end()

  call test_begin('resume_stats')
  call init_mesh(mesh_resume, v0, v1, v2)
  cfg%sim%batch_count = 8_i32
  call seed_particles_from_config(cfg)
  stats_seed = sim_stats()
  stats_seed%processed_particles = 10_i64
  stats_seed%absorbed = 4_i64
  stats_seed%escaped = 6_i64
  stats_seed%escaped_boundary = 3_i64
  stats_seed%survived_max_step = 3_i64
  stats_seed%batches = 7_i32
  stats_seed%last_rel_change = 1.0d-3
  call run_absorption_insulator(mesh_resume, cfg, stats_resume, initial_stats=stats_seed)

  call assert_equal_i64(stats_resume%processed_particles, 14_i64, 'resume processed_particles mismatch')
  call assert_equal_i64(stats_resume%absorbed, 6_i64, 'resume absorbed mismatch')
  call assert_equal_i64(stats_resume%escaped, 8_i64, 'resume escaped mismatch')
  call assert_equal_i64(stats_resume%escaped_boundary, 4_i64, 'resume escaped_boundary mismatch')
  call assert_equal_i64(stats_resume%survived_max_step, 4_i64, 'resume survived_max_step mismatch')
  call assert_equal_i32(stats_resume%batches, 8_i32, 'resume batches mismatch')
  call assert_close_dp(mesh_resume%q_elem(1), 5.0d0, 1.0d-12, 'resume deposited charge mismatch')
  call assert_true(stats_resume%last_rel_change > 0.0d0, 'resume last_rel_change should be positive')
  call test_end()

  call test_begin('collision_query_failure_context')
  call test_collision_query_failure_context()
  call test_end()

  call test_begin('photo_collision_query_failure_context')
  call test_photo_collision_query_failure_context()
  call test_end()

  call test_summary()

contains

  subroutine test_collision_query_failure_context()
    character(len=1024) :: executable_path, command, child_line
    integer :: child_exit_status, child_cmd_status, child_unit, child_ios
    logical :: saw_batch, saw_particle, saw_status, saw_survivor_warning

    call get_command_argument(0, executable_path)
    call delete_file_if_exists(collision_failure_path)
    command = 'OMP_NUM_THREADS=2 BEACH_WARN_LONG_PARTICLE_STEPS=1 "'//trim(executable_path)// &
              '" --collision-query-failure-probe > '//collision_failure_path//' 2>&1'
    call execute_command_line( &
      trim(command), wait=.true., exitstat=child_exit_status, cmdstat=child_cmd_status &
      )

    call assert_equal_i32(int(child_cmd_status, i32), 0_i32, 'collision failure probe command status mismatch')
    call assert_true(child_exit_status /= 0, 'collision failure probe should terminate with a nonzero status')

    saw_batch = .false.
    saw_particle = .false.
    saw_status = .false.
    saw_survivor_warning = .false.
    open (newunit=child_unit, file=collision_failure_path, status='old', action='read', iostat=child_ios)
    if (child_ios /= 0) error stop 'failed to read collision failure probe output'
    do
      read (child_unit, '(A)', iostat=child_ios) child_line
      if (child_ios /= 0) exit
      saw_batch = saw_batch .or. index(child_line, 'batch=1') > 0
      saw_particle = saw_particle .or. index(child_line, 'particle=1') > 0
      saw_status = saw_status .or. index(child_line, 'status=image_limit') > 0
      saw_survivor_warning = saw_survivor_warning .or. index(child_line, 'max-step-survivor') > 0
    end do
    close (child_unit)
    call delete_file_if_exists(collision_failure_path)

    call assert_true(saw_batch, 'collision failure message should include the batch index')
    call assert_true(saw_particle, 'collision failure message should include the lowest failing particle index')
    call assert_true(saw_status, 'collision failure message should include the collision status')
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

    call default_app_config(failure_cfg)
    failure_cfg%sim%rng_seed = 991_i32
    failure_cfg%sim%batch_count = 1_i32
    failure_cfg%sim%dt = 1.0d0
    failure_cfg%sim%max_step = 1_i32
    failure_cfg%sim%softening = 1.0d-6
    failure_cfg%sim%field_solver = 'fmm'
    failure_cfg%sim%field_bc_mode = 'periodic2'
    failure_cfg%sim%field_periodic_far_correction = 'none'
    failure_cfg%sim%use_box = .true.
    failure_cfg%sim%box_min = [0.0d0, 0.0d0, -1.0d0]
    failure_cfg%sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    failure_cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    failure_cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]

    failure_cfg%n_particle_species = 1_i32
    failure_cfg%particle_species(1) = species_from_defaults()
    failure_cfg%particle_species(1)%source_mode = 'volume_seed'
    failure_cfg%particle_species(1)%npcls_per_step = 2_i32
    failure_cfg%particle_species(1)%q_particle = 1.0d0
    failure_cfg%particle_species(1)%m_particle = 1.0d0
    failure_cfg%particle_species(1)%w_particle = 1.0d0
    failure_cfg%particle_species(1)%pos_low = [0.2d0, 0.25d0, 0.5d0]
    failure_cfg%particle_species(1)%pos_high = failure_cfg%particle_species(1)%pos_low
    failure_cfg%particle_species(1)%drift_velocity = [5000.0d0, 0.0d0, 0.0d0]
    failure_cfg%particle_species(1)%temperature_k = 0.0d0

    call prepare_periodic2_collision_mesh(failure_mesh, failure_cfg%sim)
    call seed_particles_from_config(failure_cfg)
    call run_absorption_insulator(failure_mesh, failure_cfg, failure_stats)
  end subroutine run_collision_query_failure_probe

  subroutine test_photo_collision_query_failure_context()
    character(len=1024) :: executable_path, command, child_line
    integer :: child_exit_status, child_cmd_status, child_unit, child_ios
    logical :: saw_batch, saw_rank, saw_species, saw_ray, saw_bounce, saw_status

    call get_command_argument(0, executable_path)
    call delete_file_if_exists(photo_collision_failure_path)
    command = 'OMP_NUM_THREADS=2 "'//trim(executable_path)// &
              '" --photo-collision-query-failure-probe > '//photo_collision_failure_path//' 2>&1'
    call execute_command_line( &
      trim(command), wait=.true., exitstat=child_exit_status, cmdstat=child_cmd_status &
      )

    call assert_equal_i32(int(child_cmd_status, i32), 0_i32, 'photo collision failure command status mismatch')
    call assert_true(child_exit_status /= 0, 'photo collision failure probe should terminate with a nonzero status')

    saw_batch = .false.
    saw_rank = .false.
    saw_species = .false.
    saw_ray = .false.
    saw_bounce = .false.
    saw_status = .false.
    open (newunit=child_unit, file=photo_collision_failure_path, status='old', action='read', iostat=child_ios)
    if (child_ios /= 0) error stop 'failed to read photo collision failure probe output'
    do
      read (child_unit, '(A)', iostat=child_ios) child_line
      if (child_ios /= 0) exit
      saw_batch = saw_batch .or. index(child_line, 'batch=1') > 0
      saw_rank = saw_rank .or. index(child_line, 'rank=0') > 0
      saw_species = saw_species .or. index(child_line, 'species=2') > 0
      saw_ray = saw_ray .or. index(child_line, 'ray=1') > 0
      saw_bounce = saw_bounce .or. index(child_line, 'bounce=0') > 0
      saw_status = saw_status .or. index(child_line, 'status=image_limit') > 0
    end do
    close (child_unit)
    call delete_file_if_exists(photo_collision_failure_path)

    call assert_true(saw_batch, 'photo collision failure message should include the batch index')
    call assert_true(saw_rank, 'photo collision failure message should include the selected rank')
    call assert_true(saw_species, 'photo collision failure message should include the species index')
    call assert_true(saw_ray, 'photo collision failure message should include the lowest failing ray')
    call assert_true(saw_bounce, 'photo collision failure message should include the failing bounce')
    call assert_true(saw_status, 'photo collision failure message should include the collision status')
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

    call default_app_config(failure_cfg)
    failure_cfg%sim%rng_seed = 992_i32
    failure_cfg%sim%batch_count = 1_i32
    failure_cfg%sim%batch_duration = 0.5d0
    failure_cfg%sim%dt = 1.0d0
    failure_cfg%sim%max_step = 1_i32
    failure_cfg%sim%softening = 1.0d-6
    failure_cfg%sim%field_solver = 'fmm'
    failure_cfg%sim%field_bc_mode = 'periodic2'
    failure_cfg%sim%field_periodic_far_correction = 'none'
    failure_cfg%sim%use_box = .true.
    failure_cfg%sim%box_min = [0.0d0, 0.0d0, 0.0d0]
    failure_cfg%sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    failure_cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    failure_cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]

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
end program test_simulator

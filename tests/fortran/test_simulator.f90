!> 吸着ループの統計更新・堆積・履歴出力を検証する統合テスト。
program test_simulator
  use bem_kinds, only: dp, i32, i64
  use bem_mesh, only: init_mesh, prepare_periodic2_collision_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_simulator_workspace, only: simulator_batch_workspace_type
  use bem_particles, only: init_particles, append_particles
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config, &
                            build_mesh_from_config
  use bem_physics_config_types, only: normalize_legacy_physics_config
  use bem_types, only: mesh_type, particles_soa, sim_stats, injection_state, bc_open, bc_reflect, bc_periodic
  use bem_charge_ledger, only: charge_ledger_type
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_equal_i64, assert_close_dp, delete_file_if_exists
  implicit none

  type(mesh_type) :: mesh
  type(mesh_type) :: mesh_tree
  type(mesh_type) :: mesh_potential_history
  type(mesh_type) :: mesh_resume
  type(app_config) :: cfg, cfg_tree, cfg_potential_history
  type(sim_stats) :: stats, stats_tree, stats_potential_history, stats_seed, stats_resume
  type(charge_ledger_type) :: charge_ledger
  real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1)
  real(dp) :: potential_value, top_time, top_z, top_mean, top_std, top_min, top_max
  integer :: u, u_top, ios
  integer(i32) :: potential_batch_idx, potential_elem_idx, top_batch_idx, top_sample_n
  character(len=256) :: line
  integer(i32) :: n_lines
  character(len=*), parameter :: history_path = 'test_simulator_history_tmp.csv'
  character(len=*), parameter :: potential_history_path = 'test_simulator_potential_history_tmp.csv'
  character(len=*), parameter :: top_reference_history_path = 'test_simulator_top_reference_history_tmp.csv'
  character(len=*), parameter :: collision_failure_path = 'test_simulator_collision_failure_tmp.log'
  character(len=*), parameter :: photo_collision_failure_path = 'test_simulator_photo_collision_failure_tmp.log'
  character(len=*), parameter :: box_event_failure_path = 'test_simulator_box_event_failure_tmp.log'
  character(len=*), parameter :: invalid_candidate_failure_path = 'test_simulator_invalid_candidate_failure_tmp.log'
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
  end if

  v0(:, 1) = [0.0d0, 0.0d0, 0.0d0]
  v1(:, 1) = [1.0d0, 0.0d0, 0.0d0]
  v2(:, 1) = [0.0d0, 1.0d0, 0.0d0]
  call init_mesh(mesh, v0, v1, v2)
  mesh%elem_vacuum_sign = 1_i32
  mesh%vacuum_normals = mesh%normals

  call default_app_config(cfg)
  cfg%sim%rng_seed = 777_i32
  cfg%sim%batch_count = 1_i32
  cfg%sim%dt = 1.0d0
  cfg%sim%max_step = 1_i32
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
  cfg%particle_species(2)%pos_low = [0.9d0, 0.1d0, 0.5d0]
  cfg%particle_species(2)%pos_high = [0.9d0, 0.1d0, 0.5d0]
  cfg%particle_species(2)%drift_velocity = [2.0d0, 0.0d0, 0.0d0]
  cfg%particle_species(2)%temperature_k = 0.0d0

  cfg%particle_species(3) = species_from_defaults()
  cfg%particle_species(3)%source_mode = 'volume_seed'
  cfg%particle_species(3)%npcls_per_step = 1_i32
  cfg%particle_species(3)%q_particle = 1.0d0
  cfg%particle_species(3)%m_particle = 1.0d0
  cfg%particle_species(3)%w_particle = 1.0d0
  cfg%particle_species(3)%pos_low = [0.1d0, 0.9d0, 0.5d0]
  cfg%particle_species(3)%pos_high = [0.1d0, 0.9d0, 0.5d0]
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

  call test_init(17)

  call test_begin('batch_workspace_reuse')
  call test_batch_workspace_reuse()
  call test_end()

  call test_begin('particle_append_preserves_existing_state')
  call test_particle_append_preserves_existing_state()
  call test_end()

  call test_begin('basic_simulation')
  call delete_file_if_exists(history_path)
  open (newunit=u, file=history_path, status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'failed to open simulator history fixture'
  call run_absorption_insulator( &
    mesh, cfg, stats, history_unit=u, history_stride=1_i32, charge_ledger=charge_ledger &
    )
  close (u)

  call assert_equal_i64(stats%processed_particles, 4_i64, 'processed_particles mismatch')
  call assert_equal_i64(stats%absorbed, 2_i64, 'absorbed mismatch')
  call assert_equal_i64(stats%escaped, 2_i64, 'escaped mismatch')
  call assert_equal_i64(stats%escaped_boundary, 2_i64, 'escaped_boundary mismatch')
  call assert_equal_i64(stats%survived_max_step, 0_i64, 'survived_max_step mismatch')
  call assert_equal_i32(stats%batches, 1_i32, 'batch count mismatch')
  call assert_close_dp(mesh%q_elem(1), 5.0d0, 1.0d-12, 'deposited charge mismatch')
  call assert_true(stats%last_rel_change > 0.0d0, 'last_rel_change should be positive')
  call assert_close_dp(charge_ledger%surface_charge_before, 0.0_dp, 1.0d-12, 'ledger surface before mismatch')
  call assert_close_dp(charge_ledger%surface_charge_after, 5.0_dp, 1.0d-12, 'ledger surface after mismatch')
  call assert_close_dp(sum(charge_ledger%injected_from_remote), 7.0_dp, 1.0d-12, 'ledger injected charge mismatch')
  call assert_close_dp(sum(charge_ledger%absorbed_on_surface), 5.0_dp, 1.0d-12, 'ledger absorbed charge mismatch')
  call assert_close_dp(sum(charge_ledger%escaped_to_infinity), 2.0_dp, 1.0d-12, 'ledger escaped charge mismatch')
  call assert_close_dp(charge_ledger%residual(), 0.0_dp, 1.0d-12, 'simulation charge residual mismatch')
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
  mesh_tree%elem_vacuum_sign = 1_i32
  mesh_tree%vacuum_normals = mesh_tree%normals
  cfg_tree = cfg
  cfg_tree%sim%field_solver = 'treecode'
  call normalize_legacy_physics_config( &
    cfg_tree%sim, cfg_tree%field, cfg_tree%periodic2, cfg_tree%panel &
    )
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
  mesh_potential_history%elem_vacuum_sign = 1_i32
  mesh_potential_history%vacuum_normals = mesh_potential_history%normals
  cfg_potential_history = cfg
  cfg_potential_history%sim%field_solver = 'fmm'
  cfg_potential_history%sim%batch_duration = 1.0_dp
  call normalize_legacy_physics_config( &
    cfg_potential_history%sim, cfg_potential_history%field, cfg_potential_history%periodic2, &
    cfg_potential_history%panel &
    )
  call seed_particles_from_config(cfg_potential_history)
  call delete_file_if_exists(potential_history_path)
  call delete_file_if_exists(top_reference_history_path)
  open (newunit=u, file=potential_history_path, status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'failed to open simulator potential history fixture'
  open (newunit=u_top, file=top_reference_history_path, status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'failed to open simulator top-reference history fixture'
  call run_absorption_insulator( &
    mesh_potential_history, cfg_potential_history, stats_potential_history, &
    potential_history_unit=u, top_reference_history_unit=u_top, history_stride=1_i32 &
    )
  close (u)
  close (u_top)
  open (newunit=u, file=potential_history_path, status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'failed to read simulator potential history fixture'
  read (u, *, iostat=ios) potential_batch_idx, potential_elem_idx, potential_value
  close (u)
  if (ios /= 0) error stop 'failed to parse simulator potential history fixture'
  call assert_equal_i32(potential_batch_idx, 1_i32, 'potential history batch mismatch')
  call assert_equal_i32(potential_elem_idx, 1_i32, 'potential history elem mismatch')
  call assert_true(potential_value > 0.0d0, 'FMM potential history should use post-commit charges')
  open (newunit=u_top, file=top_reference_history_path, status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'failed to read simulator top-reference history fixture'
  read (u_top, *, iostat=ios) top_batch_idx, top_time, top_z, top_sample_n, top_mean, top_std, top_min, top_max
  close (u_top)
  if (ios /= 0) error stop 'failed to parse simulator top-reference history fixture'
  call assert_equal_i32(top_batch_idx, 1_i32, 'top-reference history batch mismatch')
  call assert_equal_i32( &
    top_sample_n, cfg_potential_history%sim%injection_face_phi_grid_n, &
    'top-reference history sample count mismatch' &
    )
  call assert_true(top_time > 0.0_dp, 'top-reference history should record simulated time')
  call assert_close_dp(top_z, cfg_potential_history%sim%box_max(3), 1.0d-14, 'top-reference z-high mismatch')
  call assert_true( &
    top_std >= 0.0_dp .and. top_min <= top_mean .and. top_mean <= top_max, &
    'top-reference history statistics are inconsistent' &
    )
  call delete_file_if_exists(potential_history_path)
  call delete_file_if_exists(top_reference_history_path)
  call test_end()

  call test_begin('resume_stats')
  call init_mesh(mesh_resume, v0, v1, v2)
  mesh_resume%elem_vacuum_sign = 1_i32
  mesh_resume%vacuum_normals = mesh_resume%normals
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
  call assert_equal_i64(stats_resume%escaped_boundary, 5_i64, 'resume escaped_boundary mismatch')
  call assert_equal_i64(stats_resume%survived_max_step, 3_i64, 'resume survived_max_step mismatch')
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

  call test_begin('reflected_remainder_deposits')
  call test_reflected_remainder_deposits()
  call test_end()

  call test_begin('species_z_high_reflect_preserves_ambient_escape')
  call test_species_z_high_reflect_preserves_ambient_escape()
  call test_end()

  call test_begin('neutral_return_closes_mean_charge_and_preserves_redistribution')
  call test_neutral_return_closure()
  call test_end()

  call test_begin('fixed_absorbed_current_preserves_spatial_mapping')
  call test_fixed_absorbed_current_closure()
  call test_end()

  call test_begin('fixed_photo_currents_scale_emission_and_return_separately')
  call test_fixed_photo_current_closure()
  call test_end()

  call test_begin('multiple_box_event_failure_context')
  call test_multiple_box_event_failure_context()
  call test_end()

  call test_begin('multiple_box_event_soft_discard')
  call test_multiple_box_event_soft_discard()
  call test_end()

  call test_begin('invalid_candidate_failure_context')
  call test_invalid_candidate_failure_context()
  call test_end()

  call test_summary()

contains

  subroutine test_particle_append_preserves_existing_state()
    type(particles_soa) :: particles
    real(dp) :: initial_x(3, 1), initial_v(3, 1), append_x(3, 2), append_v(3, 2)

    initial_x(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp]
    initial_v(:, 1) = [4.0_dp, 5.0_dp, 6.0_dp]
    append_x(:, 1) = [7.0_dp, 8.0_dp, 9.0_dp]
    append_x(:, 2) = [10.0_dp, 11.0_dp, 12.0_dp]
    append_v(:, 1) = [13.0_dp, 14.0_dp, 15.0_dp]
    append_v(:, 2) = [16.0_dp, 17.0_dp, 18.0_dp]
    call init_particles( &
      particles, initial_x, initial_v, [1.0_dp], [2.0_dp], [3.0_dp], &
      species_id=[1_i32], source_element=[7_i32] &
      )
    particles%alive(1) = .false.
    call append_particles( &
      particles, append_x, append_v, [4.0_dp, 5.0_dp], [6.0_dp, 7.0_dp], [8.0_dp, 9.0_dp], &
      [2_i32, 3_i32], source_element=[8_i32, 9_i32] &
      )

    call assert_equal_i32(particles%n, 3_i32, 'particle append count mismatch')
    call assert_true(.not. particles%alive(1), 'particle append must preserve prior alive state')
    call assert_true(all(particles%alive(2:3)), 'particle append must activate appended particles')
    call assert_close_dp(particles%x(3, 3), 12.0_dp, 0.0_dp, 'particle append position mismatch')
    call assert_close_dp(particles%w(3), 9.0_dp, 0.0_dp, 'particle append weight mismatch')
    call assert_equal_i32(particles%species_id(3), 3_i32, 'particle append species mismatch')
    call assert_true( &
      all(particles%source_element == [7_i32, 8_i32, 9_i32]), &
      'particle append must preserve and append source provenance' &
      )
  end subroutine test_particle_append_preserves_existing_state

  subroutine test_batch_workspace_reuse()
    type(simulator_batch_workspace_type) :: workspace, regular_workspace

    call workspace%init(2_i32, 3_i32, 2_i32, candidate_charge_enabled=.true.)
    call regular_workspace%init(2_i32, 3_i32, 2_i32)
    call assert_equal_i32(int(size(workspace%dq_thread, 1), i32), 2_i32, 'workspace mesh capacity mismatch')
    call assert_equal_i32(int(size(workspace%dq_thread, 2), i32), 2_i32, 'workspace thread capacity mismatch')
    call assert_equal_i32( &
      int(size(workspace%ledger_charge_values), i32), 15_i32, 'workspace ledger charge capacity mismatch' &
      )
    call assert_equal_i32( &
      int(size(workspace%neutral_return_charge_values), i32), 18_i32, &
      'workspace neutral-return charge capacity mismatch' &
      )
    call assert_equal_i32( &
      int(size(workspace%neutral_return_terminal_counts), i32), 9_i32, &
      'workspace neutral-return terminal capacity mismatch' &
      )
    call assert_equal_i32( &
      int(size(regular_workspace%candidate_charge), i32), 0_i32, &
      'regular workspace must not allocate candidate charge storage' &
      )

    call workspace%prepare_particle_flags(5_i32)
    workspace%escaped_boundary_flag = .true.
    workspace%absorbed_flag = .true.
    workspace%absorbed_element = 1_i32
    workspace%soft_discarded_boundary_flag = .true.
    workspace%dq_thread = 1.0_dp
    workspace%photo_emission_dq = 2.0_dp
    workspace%candidate_charge = 3.0_dp
    workspace%neutral_return_charge_values = 4.0_dp
    workspace%neutral_return_terminal_counts = 5_i64
    workspace%neutral_return_weight_scale = 6.0_dp
    workspace%neutral_return_correction = 7.0_dp
    workspace%neutral_return_unresolved_fraction = 8.0_dp

    call workspace%reset_before_injection()
    call workspace%prepare_particle_flags(2_i32)
    call assert_true(all(workspace%dq_thread == 0.0_dp), 'workspace dq_thread reset mismatch')
    call assert_true(all(workspace%photo_emission_dq == 0.0_dp), 'workspace photo charge reset mismatch')
    call assert_true(all(workspace%candidate_charge == 0.0_dp), 'workspace candidate charge reset mismatch')
    call assert_true( &
      all(workspace%neutral_return_charge_values == 0.0_dp), &
      'workspace neutral-return charge reset mismatch' &
      )
    call assert_true( &
      all(workspace%neutral_return_terminal_counts == 0_i64), &
      'workspace neutral-return terminal reset mismatch' &
      )
    call assert_true( &
      all(workspace%neutral_return_weight_scale == 1.0_dp), &
      'workspace neutral-return scale reset mismatch' &
      )
    call assert_true( &
      all(workspace%neutral_return_correction == 0.0_dp), &
      'workspace neutral-return correction reset mismatch' &
      )
    call assert_true(all(.not. workspace%escaped_boundary_flag(:2)), 'workspace escaped flag reset mismatch')
    call assert_true(all(.not. workspace%absorbed_flag(:2)), 'workspace absorbed flag reset mismatch')
    call assert_true(all(workspace%absorbed_element(:2) == -1_i32), 'workspace absorbed element reset mismatch')
    call assert_true(all(.not. workspace%soft_discarded_boundary_flag(:2)), 'workspace discard flag reset mismatch')
    call assert_equal_i32( &
      int(size(workspace%escaped_boundary_flag), i32), 5_i32, 'workspace flags should retain grown capacity' &
      )

    call workspace%prepare_particle_flags(8_i32)
    call assert_equal_i32( &
      int(size(workspace%escaped_boundary_flag), i32), 8_i32, 'workspace flags should grow on demand' &
      )
    call assert_true(all(.not. workspace%escaped_boundary_flag), 'grown workspace escaped flags must start clear')
  end subroutine test_batch_workspace_reuse

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
      saw_status = saw_status .or. index(child_line, 'status=multiple_box_events') > 0
      saw_survivor_warning = saw_survivor_warning .or. index(child_line, 'max-step-survivor') > 0
    end do
    close (child_unit)
    call delete_file_if_exists(collision_failure_path)

    call assert_true(saw_batch, 'collision failure message should include the batch index')
    call assert_true(saw_particle, 'collision failure message should include the lowest failing particle index')
    call assert_true(saw_status, 'particle-step failure message should include the selected status')
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

  subroutine test_reflected_remainder_deposits()
    type(mesh_type) :: reflected_mesh
    type(app_config) :: reflected_cfg
    type(sim_stats) :: reflected_stats
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [0.5_dp, -1.0_dp, -1.0_dp]
    tri_v1(:, 1) = [0.5_dp, 1.0_dp, -1.0_dp]
    tri_v2(:, 1) = [0.5_dp, 0.0_dp, 1.0_dp]
    call init_mesh(reflected_mesh, tri_v0, tri_v1, tri_v2)
    reflected_mesh%elem_vacuum_sign = 1_i32
    reflected_mesh%vacuum_normals = reflected_mesh%normals

    call default_app_config(reflected_cfg)
    reflected_cfg%sim%rng_seed = 993_i32
    reflected_cfg%sim%batch_count = 1_i32
    reflected_cfg%sim%dt = 1.0_dp
    reflected_cfg%sim%max_step = 1_i32
    reflected_cfg%sim%use_box = .true.
    reflected_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    reflected_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    reflected_cfg%sim%bc_high(1) = bc_reflect
    reflected_cfg%n_particle_species = 1_i32
    reflected_cfg%particle_species(1) = species_from_defaults()
    reflected_cfg%particle_species(1)%source_mode = 'volume_seed'
    reflected_cfg%particle_species(1)%npcls_per_step = 1_i32
    reflected_cfg%particle_species(1)%q_particle = 1.0_dp
    reflected_cfg%particle_species(1)%m_particle = 1.0_dp
    reflected_cfg%particle_species(1)%w_particle = 2.0_dp
    reflected_cfg%particle_species(1)%pos_low = [0.8_dp, 0.2_dp, 0.2_dp]
    reflected_cfg%particle_species(1)%pos_high = reflected_cfg%particle_species(1)%pos_low
    reflected_cfg%particle_species(1)%drift_velocity = [1.0_dp, 0.0_dp, 0.0_dp]
    reflected_cfg%particle_species(1)%temperature_k = 0.0_dp

    call seed_particles_from_config(reflected_cfg)
    call run_absorption_insulator(reflected_mesh, reflected_cfg, reflected_stats)
    call assert_equal_i64(reflected_stats%absorbed, 1_i64, 'reflected remainder should absorb one particle')
    call assert_equal_i64(reflected_stats%escaped, 0_i64, 'reflected remainder should not escape')
    call assert_close_dp(reflected_mesh%q_elem(1), 2.0_dp, 1.0e-12_dp, 'reflected remainder deposit mismatch')
  end subroutine test_reflected_remainder_deposits

  subroutine test_species_z_high_reflect_preserves_ambient_escape()
    type(mesh_type) :: species_mesh
    type(app_config) :: species_cfg
    type(sim_stats) :: species_stats
    type(charge_ledger_type) :: species_ledger
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)
    integer(i32) :: species_idx

    tri_v0(:, 1) = [-1.0_dp, -1.0_dp, 0.0_dp]
    tri_v1(:, 1) = [1.0_dp, -1.0_dp, 0.0_dp]
    tri_v2(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
    call init_mesh(species_mesh, tri_v0, tri_v1, tri_v2)
    species_mesh%elem_vacuum_sign = 1_i32
    species_mesh%vacuum_normals = species_mesh%normals

    call default_app_config(species_cfg)
    species_cfg%sim%rng_seed = 997_i32
    species_cfg%sim%batch_count = 1_i32
    species_cfg%sim%dt = 1.0_dp
    species_cfg%sim%max_step = 1_i32
    species_cfg%sim%use_box = .true.
    species_cfg%sim%box_min = [-1.0_dp, -1.0_dp, -1.0_dp]
    species_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    species_cfg%sim%bc_high(3) = bc_open
    species_cfg%n_particle_species = 2_i32
    do species_idx = 1_i32, species_cfg%n_particle_species
      species_cfg%particle_species(species_idx) = species_from_defaults()
      species_cfg%particle_species(species_idx)%source_mode = 'volume_seed'
      species_cfg%particle_species(species_idx)%npcls_per_step = 1_i32
      species_cfg%particle_species(species_idx)%q_particle = 1.0_dp
      species_cfg%particle_species(species_idx)%m_particle = 1.0_dp
      species_cfg%particle_species(species_idx)%w_particle = real(species_idx, dp)
      species_cfg%particle_species(species_idx)%pos_low = [0.0_dp, 0.0_dp, 0.5_dp]
      species_cfg%particle_species(species_idx)%pos_high = species_cfg%particle_species(species_idx)%pos_low
      species_cfg%particle_species(species_idx)%drift_velocity = [0.0_dp, 0.0_dp, 2.0_dp]
      species_cfg%particle_species(species_idx)%temperature_k = 0.0_dp
    end do
    species_cfg%particle_species(2)%boundary_high(3) = bc_reflect

    call seed_particles_from_config(species_cfg)
    call run_absorption_insulator(species_mesh, species_cfg, species_stats, charge_ledger=species_ledger)
    call assert_equal_i64(species_stats%absorbed, 1_i64, 'species-reflected particle should return to the mesh')
    call assert_equal_i64(species_stats%escaped_boundary, 1_i64, 'inherited ambient particle should escape z-high')
    call assert_equal_i64(species_stats%survived_max_step, 0_i64, 'species z-high boundary fixture must resolve both particles')
    call assert_close_dp( &
      species_ledger%escaped_to_infinity(1), 1.0_dp, 1.0e-12_dp, 'ambient escaped charge mismatch' &
      )
    call assert_close_dp( &
      species_ledger%escaped_to_infinity(2), 0.0_dp, 1.0e-12_dp, 'reflected species must not escape' &
      )
    call assert_close_dp( &
      species_ledger%absorbed_on_surface(2), 2.0_dp, 1.0e-12_dp, 'reflected species absorption mismatch' &
      )
    call assert_close_dp(species_ledger%residual(), 0.0_dp, 1.0e-12_dp, 'species boundary ledger residual mismatch')
  end subroutine test_species_z_high_reflect_preserves_ambient_escape

  subroutine test_neutral_return_closure()
    type(mesh_type) :: neutral_mesh
    type(app_config) :: neutral_cfg
    type(sim_stats) :: neutral_stats
    type(charge_ledger_type) :: neutral_ledger
    real(dp) :: tri_v0(3, 4), tri_v1(3, 4), tri_v2(3, 4)

    ! 左98%の面は z=0.75 にあり、dt=0.6 のうちに z-high 反射後同じ面へ帰還する。
    ! 右2%の面は z=0.25 にあり、同じ dt では帰還せず max-step survivor になる。
    tri_v0(:, 1) = [0.0_dp, 0.0_dp, 0.75_dp]
    tri_v1(:, 1) = [0.98_dp, 0.0_dp, 0.75_dp]
    tri_v2(:, 1) = [0.98_dp, 1.0_dp, 0.75_dp]
    tri_v0(:, 2) = [0.0_dp, 0.0_dp, 0.75_dp]
    tri_v1(:, 2) = [0.98_dp, 1.0_dp, 0.75_dp]
    tri_v2(:, 2) = [0.0_dp, 1.0_dp, 0.75_dp]
    tri_v0(:, 3) = [0.98_dp, 0.0_dp, 0.25_dp]
    tri_v1(:, 3) = [1.0_dp, 0.0_dp, 0.25_dp]
    tri_v2(:, 3) = [1.0_dp, 1.0_dp, 0.25_dp]
    tri_v0(:, 4) = [0.98_dp, 0.0_dp, 0.25_dp]
    tri_v1(:, 4) = [1.0_dp, 1.0_dp, 0.25_dp]
    tri_v2(:, 4) = [0.98_dp, 1.0_dp, 0.25_dp]
    call init_mesh(neutral_mesh, tri_v0, tri_v1, tri_v2)
    neutral_mesh%elem_vacuum_sign = 1_i32
    neutral_mesh%vacuum_normals = neutral_mesh%normals

    call default_app_config(neutral_cfg)
    neutral_cfg%sim%rng_seed = 998_i32
    neutral_cfg%sim%batch_count = 1_i32
    neutral_cfg%sim%batch_duration = 1.0_dp
    neutral_cfg%sim%dt = 0.6_dp
    neutral_cfg%sim%max_step = 1_i32
    neutral_cfg%sim%q_floor = 1.0e-30_dp
    neutral_cfg%sim%use_box = .true.
    neutral_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    neutral_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    neutral_cfg%sim%bc_high(3) = bc_open
    neutral_cfg%n_particle_species = 1_i32
    neutral_cfg%particle_species(1) = species_from_defaults()
    neutral_cfg%particle_species(1)%source_mode = 'photo_raycast'
    neutral_cfg%particle_species(1)%rays_per_batch = 1024_i32
    neutral_cfg%particle_species(1)%emit_current_density_a_m2 = 1.0_dp
    neutral_cfg%particle_species(1)%deposit_opposite_charge_on_emit = .true.
    neutral_cfg%particle_species(1)%boundary_high(3) = bc_reflect
    neutral_cfg%particle_species(1)%surface_charge_closure = 'neutral_return'
    neutral_cfg%particle_species(1)%q_particle = -1.0_dp
    neutral_cfg%particle_species(1)%m_particle = 1.0_dp
    neutral_cfg%particle_species(1)%temperature_k = 0.0_dp
    neutral_cfg%particle_species(1)%normal_drift_speed = 1.0_dp
    neutral_cfg%particle_species(1)%inject_face = 'z_high'
    neutral_cfg%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, 1.0_dp]
    neutral_cfg%particle_species(1)%pos_high = [1.0_dp, 1.0_dp, 1.0_dp]
    neutral_cfg%particle_species(1)%ray_direction = [0.0_dp, 0.0_dp, -1.0_dp]
    neutral_cfg%particle_species(1)%has_ray_direction = .true.

    call seed_particles_from_config(neutral_cfg)
    call run_absorption_insulator( &
      neutral_mesh, neutral_cfg, neutral_stats, charge_ledger=neutral_ledger &
      )

    call assert_true(neutral_ledger%absorbed_count(1) > 0_i64, 'neutral-return fixture needs resolved returns')
    call assert_true( &
      neutral_ledger%discarded_unresolved_count(1) > 0_i64, &
      'neutral-return fixture needs unresolved survivors' &
      )
    call assert_true( &
      neutral_ledger%neutral_return_unresolved_fraction(1) <= 0.05_dp, &
      'neutral-return fixture must stay inside the fixed applicability limit' &
      )
    call assert_equal_i64(neutral_stats%escaped_boundary, 0_i64, 'neutral-return photoelectrons must not escape')
    call assert_close_dp( &
      neutral_ledger%emitted_from_surface(1), &
      neutral_ledger%absorbed_on_surface(1) + neutral_ledger%discarded_unresolved(1), &
      1.0e-12_dp, 'raw neutral-return particle charge must close' &
      )
    call assert_close_dp( &
      neutral_ledger%neutral_return_correction(1), neutral_ledger%discarded_unresolved(1), &
      1.0e-12_dp, 'neutral-return correction must represent the unresolved raw charge' &
      )
    call assert_close_dp( &
      neutral_ledger%neutral_return_weight_scale(1), &
      neutral_ledger%emitted_from_surface(1)/neutral_ledger%absorbed_on_surface(1), &
      1.0e-12_dp, 'neutral-return weight scale mismatch' &
      )
    call assert_close_dp(sum(neutral_mesh%q_elem), 0.0_dp, 1.0e-12_dp, 'neutral-return mean charge must be zero')
    call assert_true(sum(abs(neutral_mesh%q_elem)) > 0.0_dp, 'neutral-return must preserve local charge redistribution')
    call assert_close_dp(neutral_ledger%residual(), 0.0_dp, 1.0e-12_dp, 'neutral-return ledger residual mismatch')
  end subroutine test_neutral_return_closure

  subroutine test_fixed_absorbed_current_closure()
    type(mesh_type) :: fixed_mesh
    type(app_config) :: fixed_cfg
    type(sim_stats) :: fixed_stats
    type(charge_ledger_type) :: fixed_ledger
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    tri_v1(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
    tri_v2(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
    call init_mesh(fixed_mesh, tri_v0, tri_v1, tri_v2)
    fixed_mesh%elem_vacuum_sign = 1_i32
    fixed_mesh%vacuum_normals = fixed_mesh%normals

    call default_app_config(fixed_cfg)
    fixed_cfg%sim%batch_count = 1_i32
    fixed_cfg%sim%batch_duration = 1.0_dp
    fixed_cfg%sim%dt = 1.0_dp
    fixed_cfg%sim%max_step = 1_i32
    fixed_cfg%sim%q_floor = 1.0e-30_dp
    fixed_cfg%n_particle_species = 1_i32
    fixed_cfg%particle_species(1) = species_from_defaults()
    fixed_cfg%particle_species(1)%source_mode = 'volume_seed'
    fixed_cfg%particle_species(1)%npcls_per_step = 1_i32
    fixed_cfg%particle_species(1)%q_particle = 1.0_dp
    fixed_cfg%particle_species(1)%m_particle = 1.0_dp
    fixed_cfg%particle_species(1)%w_particle = 2.0_dp
    fixed_cfg%particle_species(1)%pos_low = [0.2_dp, 0.2_dp, 0.5_dp]
    fixed_cfg%particle_species(1)%pos_high = fixed_cfg%particle_species(1)%pos_low
    fixed_cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -1.0_dp]
    fixed_cfg%particle_species(1)%temperature_k = 0.0_dp
    fixed_cfg%particle_species(1)%surface_charge_closure = 'fixed_current'
    fixed_cfg%particle_species(1)%target_absorbed_current_a = 5.0_dp
    fixed_cfg%particle_species(1)%has_target_absorbed_current_a = .true.

    call seed_particles_from_config(fixed_cfg)
    call run_absorption_insulator(fixed_mesh, fixed_cfg, fixed_stats, charge_ledger=fixed_ledger)
    call assert_close_dp(fixed_mesh%q_elem(1), 5.0_dp, 1.0e-12_dp, 'fixed absorbed target charge mismatch')
    call assert_close_dp( &
      fixed_ledger%fixed_absorbed_weight_scale(1), 2.5_dp, 1.0e-12_dp, &
      'fixed absorbed current scale mismatch' &
      )
    call assert_close_dp( &
      fixed_ledger%fixed_current_correction(1), 3.0_dp, 1.0e-12_dp, &
      'fixed absorbed current correction mismatch' &
      )
    call assert_close_dp(fixed_ledger%residual(), 0.0_dp, 1.0e-12_dp, 'fixed absorbed ledger residual mismatch')
  end subroutine test_fixed_absorbed_current_closure

  subroutine test_fixed_photo_current_closure()
    type(mesh_type) :: fixed_mesh
    type(app_config) :: fixed_cfg
    type(sim_stats) :: fixed_stats
    type(charge_ledger_type) :: fixed_ledger
    real(dp) :: tri_v0(3, 2), tri_v1(3, 2), tri_v2(3, 2)

    tri_v0(:, 1) = [0.0_dp, 0.0_dp, 0.75_dp]
    tri_v1(:, 1) = [1.0_dp, 0.0_dp, 0.75_dp]
    tri_v2(:, 1) = [1.0_dp, 1.0_dp, 0.75_dp]
    tri_v0(:, 2) = [0.0_dp, 0.0_dp, 0.75_dp]
    tri_v1(:, 2) = [1.0_dp, 1.0_dp, 0.75_dp]
    tri_v2(:, 2) = [0.0_dp, 1.0_dp, 0.75_dp]
    call init_mesh(fixed_mesh, tri_v0, tri_v1, tri_v2)
    fixed_mesh%elem_vacuum_sign = 1_i32
    fixed_mesh%vacuum_normals = fixed_mesh%normals

    call default_app_config(fixed_cfg)
    fixed_cfg%sim%rng_seed = 999_i32
    fixed_cfg%sim%batch_count = 1_i32
    fixed_cfg%sim%batch_duration = 1.0_dp
    fixed_cfg%sim%dt = 0.6_dp
    fixed_cfg%sim%max_step = 1_i32
    fixed_cfg%sim%q_floor = 1.0e-30_dp
    fixed_cfg%sim%use_box = .true.
    fixed_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    fixed_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    fixed_cfg%sim%bc_high(3) = bc_open
    fixed_cfg%n_particle_species = 1_i32
    fixed_cfg%particle_species(1) = species_from_defaults()
    fixed_cfg%particle_species(1)%source_mode = 'photo_raycast'
    fixed_cfg%particle_species(1)%rays_per_batch = 256_i32
    fixed_cfg%particle_species(1)%emit_current_density_a_m2 = 1.0_dp
    fixed_cfg%particle_species(1)%deposit_opposite_charge_on_emit = .true.
    fixed_cfg%particle_species(1)%boundary_high(3) = bc_reflect
    fixed_cfg%particle_species(1)%surface_charge_closure = 'fixed_current'
    fixed_cfg%particle_species(1)%target_absorbed_current_a = -0.25_dp
    fixed_cfg%particle_species(1)%has_target_absorbed_current_a = .true.
    fixed_cfg%particle_species(1)%target_emission_current_a = 2.0_dp
    fixed_cfg%particle_species(1)%has_target_emission_current_a = .true.
    fixed_cfg%particle_species(1)%q_particle = -1.0_dp
    fixed_cfg%particle_species(1)%m_particle = 1.0_dp
    fixed_cfg%particle_species(1)%temperature_k = 0.0_dp
    fixed_cfg%particle_species(1)%normal_drift_speed = 1.0_dp
    fixed_cfg%particle_species(1)%inject_face = 'z_high'
    fixed_cfg%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, 1.0_dp]
    fixed_cfg%particle_species(1)%pos_high = [1.0_dp, 1.0_dp, 1.0_dp]
    fixed_cfg%particle_species(1)%ray_direction = [0.0_dp, 0.0_dp, -1.0_dp]
    fixed_cfg%particle_species(1)%has_ray_direction = .true.

    call seed_particles_from_config(fixed_cfg)
    call run_absorption_insulator(fixed_mesh, fixed_cfg, fixed_stats, charge_ledger=fixed_ledger)
    call assert_equal_i64(fixed_stats%absorbed, 256_i64, 'fixed PE fixture must resolve every return')
    call assert_close_dp(sum(fixed_mesh%q_elem), 1.75_dp, 1.0e-12_dp, 'fixed PE net surface charge mismatch')
    call assert_close_dp( &
      fixed_ledger%fixed_absorbed_target_charge(1), -0.25_dp, 1.0e-12_dp, &
      'fixed PE return target mismatch' &
      )
    call assert_close_dp( &
      fixed_ledger%fixed_emission_target_charge(1), 2.0_dp, 1.0e-12_dp, &
      'fixed PE emission target mismatch' &
      )
    call assert_true( &
      abs(fixed_ledger%fixed_absorbed_weight_scale(1) - fixed_ledger%fixed_emission_weight_scale(1)) > 1.0e-6_dp, &
      'fixed PE emission and return channels must use independent scales' &
      )
    call assert_close_dp(fixed_ledger%residual(), 0.0_dp, 1.0e-12_dp, 'fixed PE ledger residual mismatch')
  end subroutine test_fixed_photo_current_closure

  subroutine test_multiple_box_event_failure_context()
    character(len=1024) :: executable_path, command, child_line
    integer :: child_exit_status, child_cmd_status, child_unit, child_ios
    logical :: saw_particle, saw_step, saw_status, saw_dt, saw_x, saw_v

    call get_command_argument(0, executable_path)
    call delete_file_if_exists(box_event_failure_path)
    command = 'OMP_NUM_THREADS=2 "'//trim(executable_path)// &
              '" --multiple-box-event-probe > '//box_event_failure_path//' 2>&1'
    call execute_command_line(trim(command), wait=.true., exitstat=child_exit_status, cmdstat=child_cmd_status)
    call assert_equal_i32(int(child_cmd_status, i32), 0_i32, 'box event failure command status mismatch')
    call assert_true(child_exit_status /= 0, 'multiple box event probe should terminate with nonzero status')

    saw_particle = .false.
    saw_step = .false.
    saw_status = .false.
    saw_dt = .false.
    saw_x = .false.
    saw_v = .false.
    open (newunit=child_unit, file=box_event_failure_path, status='old', action='read', iostat=child_ios)
    if (child_ios /= 0) error stop 'failed to read box event failure probe output'
    do
      read (child_unit, '(A)', iostat=child_ios) child_line
      if (child_ios /= 0) exit
      saw_particle = saw_particle .or. index(child_line, 'particle=1') > 0
      saw_step = saw_step .or. index(child_line, 'step=1') > 0
      saw_status = saw_status .or. index(child_line, 'status=multiple_box_events') > 0
      saw_dt = saw_dt .or. index(child_line, 'dt=') > 0
      saw_x = saw_x .or. index(child_line, 'x=') > 0
      saw_v = saw_v .or. index(child_line, 'v=') > 0
    end do
    close (child_unit)
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

  subroutine test_multiple_box_event_soft_discard()
    type(mesh_type) :: soft_mesh
    type(app_config) :: soft_cfg
    type(sim_stats) :: soft_stats
    type(charge_ledger_type) :: soft_ledger
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [10.0_dp, -1.0_dp, -1.0_dp]
    tri_v1(:, 1) = [10.0_dp, 1.0_dp, -1.0_dp]
    tri_v2(:, 1) = [10.0_dp, 0.0_dp, 1.0_dp]
    call init_mesh(soft_mesh, tri_v0, tri_v1, tri_v2)
    soft_mesh%elem_vacuum_sign = 1_i32
    soft_mesh%vacuum_normals = soft_mesh%normals
    call default_app_config(soft_cfg)
    soft_cfg%sim%rng_seed = 995_i32
    soft_cfg%sim%batch_count = 1_i32
    soft_cfg%sim%dt = 1.0_dp
    soft_cfg%sim%max_step = 1_i32
    soft_cfg%sim%use_box = .true.
    soft_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    soft_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    soft_cfg%sim%bc_low(1) = bc_reflect
    soft_cfg%sim%bc_high(1) = bc_reflect
    soft_cfg%sim%multiple_box_events_policy = 'soft_discard'
    soft_cfg%sim%multiple_box_events_soft_discard_count_limit = 100000_i32
    soft_cfg%sim%multiple_box_events_soft_discard_abs_charge_limit = 1.0e9_dp
    soft_cfg%n_particle_species = 1_i32
    soft_cfg%particle_species(1) = species_from_defaults()
    soft_cfg%particle_species(1)%source_mode = 'volume_seed'
    soft_cfg%particle_species(1)%npcls_per_step = 1_i32
    soft_cfg%particle_species(1)%q_particle = 2.0_dp
    soft_cfg%particle_species(1)%m_particle = 1.0_dp
    soft_cfg%particle_species(1)%w_particle = 3.0_dp
    soft_cfg%particle_species(1)%pos_low = [0.9_dp, 0.2_dp, 0.2_dp]
    soft_cfg%particle_species(1)%pos_high = soft_cfg%particle_species(1)%pos_low
    soft_cfg%particle_species(1)%drift_velocity = [9.0_dp, 0.0_dp, 0.0_dp]
    soft_cfg%particle_species(1)%temperature_k = 0.0_dp
    call seed_particles_from_config(soft_cfg)

    call run_absorption_insulator(soft_mesh, soft_cfg, soft_stats, charge_ledger=soft_ledger)
    call assert_equal_i64(soft_stats%processed_particles, 1_i64, 'soft discard processed count mismatch')
    call assert_equal_i64(soft_stats%absorbed, 0_i64, 'soft discard absorbed count mismatch')
    call assert_equal_i64(soft_stats%escaped, 1_i64, 'soft discard escaped umbrella count mismatch')
    call assert_equal_i64( &
      soft_stats%multiple_box_events_soft_discarded, 1_i64, 'soft discard event count mismatch' &
      )
    call assert_close_dp( &
      soft_stats%multiple_box_events_soft_discarded_abs_charge, 6.0_dp, 1.0e-12_dp, &
      'soft discard absolute charge mismatch' &
      )
    call assert_equal_i64( &
      soft_ledger%discarded_unresolved_count(1), 1_i64, 'soft discard ledger count mismatch' &
      )
    call assert_close_dp( &
      soft_ledger%discarded_unresolved(1), 6.0_dp, 1.0e-12_dp, 'soft discard ledger charge mismatch' &
      )
    call assert_close_dp(soft_ledger%residual(), 0.0_dp, 1.0e-12_dp, 'soft discard ledger residual mismatch')
  end subroutine test_multiple_box_event_soft_discard

  subroutine test_invalid_candidate_failure_context()
    character(len=1024) :: executable_path, command, child_line
    integer :: child_exit_status, child_cmd_status, child_unit, child_ios
    logical :: saw_status, saw_survivor_warning

    call get_command_argument(0, executable_path)
    call delete_file_if_exists(invalid_candidate_failure_path)
    command = 'OMP_NUM_THREADS=2 "'//trim(executable_path)// &
              '" --invalid-candidate-probe > '//invalid_candidate_failure_path//' 2>&1'
    call execute_command_line(trim(command), wait=.true., exitstat=child_exit_status, cmdstat=child_cmd_status)
    call assert_equal_i32(int(child_cmd_status, i32), 0_i32, 'invalid candidate probe command status mismatch')
    call assert_true(child_exit_status /= 0, 'invalid candidate probe should terminate with nonzero status')

    saw_status = .false.
    saw_survivor_warning = .false.
    open (newunit=child_unit, file=invalid_candidate_failure_path, status='old', action='read', iostat=child_ios)
    if (child_ios /= 0) error stop 'failed to read invalid candidate probe output'
    do
      read (child_unit, '(A)', iostat=child_ios) child_line
      if (child_ios /= 0) exit
      saw_status = saw_status .or. index(child_line, 'status=invalid_boundary') > 0
      saw_survivor_warning = saw_survivor_warning .or. index(child_line, 'max-step-survivor') > 0
    end do
    close (child_unit)
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
end program test_simulator

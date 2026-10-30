!> 吸着ループの統計更新・電荷堆積・履歴出力・再開と、境界反射の扱いを検証する統合テスト。
program test_simulator_loop
  use bem_kinds, only: dp, i32, i64
  use bem_mesh, only: init_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_run_outputs, only: run_outputs_type
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config
  use bem_physics_config_types, only: normalize_legacy_physics_config
  use bem_types, only: mesh_type, sim_stats, bc_open, bc_reflect
  use bem_charge_ledger, only: charge_ledger_type
  use bem_simulator_workspace, only: simulator_batch_workspace_type
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_equal_i64, assert_close_dp, delete_file_if_exists
  implicit none

  type(run_outputs_type) :: run_outputs

  type(mesh_type) :: mesh
  type(mesh_type) :: mesh_tree
  type(mesh_type) :: mesh_potential_history
  type(mesh_type) :: mesh_resume
  type(app_config) :: cfg, cfg_tree, cfg_potential_history
  type(sim_stats) :: stats, stats_tree, stats_potential_history, stats_seed, stats_resume
  type(charge_ledger_type) :: charge_ledger
  real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1)
  real(dp) :: potential_value
  integer :: u, ios
  integer(i32) :: potential_batch_idx, potential_elem_idx
  character(len=256) :: line
  integer(i32) :: n_lines
  character(len=*), parameter :: history_path = 'test_simulator_loop_history_tmp.csv'
  character(len=*), parameter :: potential_history_path = 'test_simulator_loop_potential_history_tmp.csv'

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

  call test_init(8)

  call test_begin('basic_simulation')
  call delete_file_if_exists(history_path)
  open (newunit=u, file=history_path, status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'failed to open simulator history fixture'
  run_outputs = run_outputs_type(history_unit=u, history_stride=1_i32)
  call run_absorption_insulator(mesh, cfg, stats, charge_ledger=charge_ledger, outputs=run_outputs)
  close (u)

  call assert_equal_i64(stats%processed_particles, 4_i64, 'processed_particles mismatch')
  call assert_equal_i64(stats%absorbed, 2_i64, 'absorbed mismatch')
  call assert_equal_i64(stats%escaped, 2_i64, 'escaped mismatch')
  call assert_equal_i64(stats%escaped_boundary, 2_i64, 'escaped_boundary mismatch')
  call assert_equal_i64(stats%survived_max_step, 0_i64, 'survived_max_step mismatch')
  call assert_equal_i32(stats%batches, 1_i32, 'batch count mismatch')
  call assert_equal_i64( &
    stats%adaptive_nonzero_mode%rejected_trials, 0_i64, &
    'non-adaptive run must not record rejected adaptive trials' &
    )
  call assert_equal_i32( &
    stats%adaptive_nonzero_mode%omp_threads, 0_i32, &
    'non-adaptive run must not record an adaptive OpenMP team size' &
    )
  call assert_close_dp( &
    stats%adaptive_nonzero_mode%last_batch_duration, 0.0_dp, 0.0_dp, &
    'non-adaptive run must not record an adaptive batch duration' &
    )
  call assert_close_dp( &
    stats%adaptive_nonzero_mode%last_potential_step, 0.0_dp, 0.0_dp, &
    'non-adaptive run must not record an adaptive potential step' &
    )
  call assert_close_dp(mesh%q_elem(1), 5.0d0, 1.0d-12, 'deposited charge mismatch')
  call assert_true(stats%last_rel_change > 0.0d0, 'last_rel_change should be positive')
  call assert_close_dp(charge_ledger%surface_charge_before, 0.0_dp, 1.0d-12, 'ledger surface before mismatch')
  call assert_close_dp(charge_ledger%surface_charge_after, 5.0_dp, 1.0d-12, 'ledger surface after mismatch')
  call assert_close_dp(sum(charge_ledger%injected_from_remote), 7.0_dp, 1.0d-12, 'ledger injected charge mismatch')
  call assert_close_dp(sum(charge_ledger%absorbed_on_surface), 5.0_dp, 1.0d-12, 'ledger absorbed charge mismatch')
  call assert_close_dp(sum(charge_ledger%escaped_to_infinity), 2.0_dp, 1.0d-12, 'ledger escaped charge mismatch')
  call assert_close_dp(charge_ledger%residual(), 0.0_dp, 1.0d-12, 'simulation charge residual mismatch')
  ! One batch with history_stride=1 writes one history snapshot row.
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
  ! The top-reference history layout is test_output_writer_io's contract; here only the ordering matters:
  ! the potential written for a batch must come from the charges committed in that batch.
  call delete_file_if_exists(potential_history_path)
  open (newunit=u, file=potential_history_path, status='replace', action='write', iostat=ios)
  if (ios /= 0) error stop 'failed to open simulator potential history fixture'
  run_outputs = run_outputs_type(potential_history_unit=u, history_stride=1_i32)
  call run_absorption_insulator( &
    mesh_potential_history, cfg_potential_history, stats_potential_history, outputs=run_outputs &
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
  ! schema v7以前のnon-adaptive出力が誤って保持した値をresume時に除去する。
  stats_seed%adaptive_nonzero_mode%rejected_trials = 9_i64
  stats_seed%adaptive_nonzero_mode%last_batch_duration = 0.25_dp
  stats_seed%adaptive_nonzero_mode%last_potential_step = 3.0_dp
  stats_seed%adaptive_nonzero_mode%omp_threads = 4_i32
  call run_absorption_insulator(mesh_resume, cfg, stats_resume, initial_stats=stats_seed)

  call assert_equal_i64(stats_resume%processed_particles, 14_i64, 'resume processed_particles mismatch')
  call assert_equal_i64(stats_resume%absorbed, 6_i64, 'resume absorbed mismatch')
  call assert_equal_i64(stats_resume%escaped, 8_i64, 'resume escaped mismatch')
  call assert_equal_i64(stats_resume%escaped_boundary, 5_i64, 'resume escaped_boundary mismatch')
  call assert_equal_i64(stats_resume%survived_max_step, 3_i64, 'resume survived_max_step mismatch')
  call assert_equal_i32(stats_resume%batches, 8_i32, 'resume batches mismatch')
  call assert_equal_i64( &
    stats_resume%adaptive_nonzero_mode%rejected_trials, 0_i64, &
    'non-adaptive resume must clear stale adaptive rejected trials' &
    )
  call assert_equal_i32( &
    stats_resume%adaptive_nonzero_mode%omp_threads, 0_i32, &
    'non-adaptive resume must clear stale adaptive OpenMP team size' &
    )
  call assert_close_dp( &
    stats_resume%adaptive_nonzero_mode%last_batch_duration, 0.0_dp, 0.0_dp, &
    'non-adaptive resume must clear stale adaptive batch duration' &
    )
  call assert_close_dp( &
    stats_resume%adaptive_nonzero_mode%last_potential_step, 0.0_dp, 0.0_dp, &
    'non-adaptive resume must clear stale adaptive potential step' &
    )
  call assert_close_dp(mesh_resume%q_elem(1), 5.0d0, 1.0d-12, 'resume deposited charge mismatch')
  call assert_true(stats_resume%last_rel_change > 0.0d0, 'resume last_rel_change should be positive')
  call test_end()

  call test_begin('reflected_remainder_deposits')
  call test_reflected_remainder_deposits()
  call test_end()

  call test_begin('species_z_high_reflect_preserves_ambient_escape')
  call test_species_z_high_reflect_preserves_ambient_escape()
  call test_end()

  call test_begin('large_finite_charge_norm_stays_finite')
  call test_large_finite_charge_norm_stays_finite()
  call test_end()

  call test_begin('batch_workspace_reset_clears_previous_batch')
  call test_batch_workspace_reset_clears_previous_batch()
  call test_end()

  call test_summary()

contains

  !> The batch workspace is reused across batches. Whatever a finished batch left in its charge tallies,
  !! closure accumulators and particle terminal flags must not leak into the next batch.
  subroutine test_batch_workspace_reset_clears_previous_batch()
    type(simulator_batch_workspace_type) :: workspace

    call workspace%init(2_i32, 3_i32, 2_i32, candidate_charge_enabled=.true.)
    call workspace%prepare_particle_flags(5_i32)
    workspace%escaped_boundary_flag = .true.
    workspace%absorbed_flag = .true.
    workspace%absorbed_element = 1_i32
    workspace%soft_discarded_boundary_flag = .true.
    workspace%dq_thread = 1.0_dp
    workspace%photo_emission_dq = 2.0_dp
    workspace%candidate_charge = 3.0_dp
    workspace%charge_candidate_ready = .true.
    workspace%neutral_return_charge_values = 4.0_dp
    workspace%neutral_return_terminal_counts = 5_i64
    workspace%neutral_return_weight_scale = 6.0_dp
    workspace%neutral_return_correction = 7.0_dp
    workspace%fixed_current_charge_values = 8.0_dp
    workspace%fixed_absorbed_weight_scale = 9.0_dp
    workspace%fixed_current_correction = 10.0_dp
    workspace%matching_plane_moments_thread = 11.0_dp
    workspace%secondary_counts_thread = 12_i64

    call workspace%reset_before_injection()
    call workspace%prepare_particle_flags(2_i32)

    call assert_true(all(workspace%dq_thread == 0.0_dp), 'previous batch charge deltas leaked')
    call assert_true(all(workspace%photo_emission_dq == 0.0_dp), 'previous batch photo-emission charge leaked')
    call assert_true(all(workspace%candidate_charge == 0.0_dp), 'previous adaptive candidate charge leaked')
    call assert_true(.not. workspace%charge_candidate_ready, 'previous adaptive candidate stayed ready')
    call assert_true(all(workspace%neutral_return_charge_values == 0.0_dp), 'neutral-return charges leaked')
    call assert_true(all(workspace%neutral_return_terminal_counts == 0_i64), 'neutral-return counts leaked')
    call assert_true(all(workspace%neutral_return_weight_scale == 1.0_dp), 'neutral-return scale was not reset to 1')
    call assert_true(all(workspace%neutral_return_correction == 0.0_dp), 'neutral-return correction leaked')
    call assert_true(all(workspace%fixed_current_charge_values == 0.0_dp), 'fixed-current tallies leaked')
    call assert_true(all(workspace%fixed_absorbed_weight_scale == 1.0_dp), 'fixed-current scale was not reset to 1')
    call assert_true(all(workspace%fixed_current_correction == 0.0_dp), 'fixed-current correction leaked')
    call assert_true(all(workspace%matching_plane_moments_thread == 0.0_dp), 'matching-plane moments leaked')
    call assert_true(all(workspace%secondary_counts_thread == 0_i64), 'secondary-emission counts leaked')
    call assert_true(all(.not. workspace%escaped_boundary_flag(:2)), 'previous escaped flags leaked')
    call assert_true(all(.not. workspace%absorbed_flag(:2)), 'previous absorbed flags leaked')
    call assert_true(all(workspace%absorbed_element(:2) == -1_i32), 'previous absorbed elements leaked')
    call assert_true(all(.not. workspace%soft_discarded_boundary_flag(:2)), 'previous discard flags leaked')
  end subroutine test_batch_workspace_reset_clears_previous_batch

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
    species_cfg%n_particle_species = 3_i32
    do species_idx = 1_i32, species_cfg%n_particle_species
      species_cfg%particle_species(species_idx) = species_from_defaults()
      species_cfg%particle_species(species_idx)%source_mode = 'volume_seed'
      species_cfg%particle_species(species_idx)%npcls_per_step = 64_i32
      species_cfg%particle_species(species_idx)%q_particle = 1.0_dp
      species_cfg%particle_species(species_idx)%m_particle = 1.0_dp
      species_cfg%particle_species(species_idx)%w_particle = real(species_idx, dp)
      species_cfg%particle_species(species_idx)%pos_low = [0.0_dp, 0.0_dp, 0.5_dp]
      species_cfg%particle_species(species_idx)%pos_high = species_cfg%particle_species(species_idx)%pos_low
      species_cfg%particle_species(species_idx)%drift_velocity = [0.0_dp, 0.0_dp, 2.0_dp]
      species_cfg%particle_species(species_idx)%temperature_k = 0.0_dp
    end do
    species_cfg%particle_species(2)%enabled = .false.
    species_cfg%particle_species(3)%boundary_high(3) = bc_reflect

    call seed_particles_from_config(species_cfg)
    call run_absorption_insulator(species_mesh, species_cfg, species_stats, charge_ledger=species_ledger)
    call assert_equal_i64(species_stats%absorbed, 64_i64, 'species-reflected particles should return to the mesh')
    call assert_equal_i64(species_stats%escaped_boundary, 64_i64, 'inherited ambient particles should escape z-high')
    call assert_equal_i64(species_stats%survived_max_step, 0_i64, 'species z-high boundary fixture must resolve both particles')
    call assert_close_dp( &
      species_ledger%escaped_to_infinity(1), 64.0_dp, 1.0e-12_dp, 'ambient escaped charge mismatch' &
      )
    call assert_close_dp( &
      species_ledger%escaped_to_infinity(3), 0.0_dp, 1.0e-12_dp, 'reflected species must not escape' &
      )
    call assert_close_dp( &
      species_ledger%absorbed_on_surface(3), 192.0_dp, 1.0e-12_dp, 'reflected species absorption mismatch' &
      )
    call assert_close_dp(species_ledger%absorbed_on_surface(2), 0.0_dp, 0.0_dp, 'disabled species must not deposit')
    call assert_close_dp(species_ledger%escaped_to_infinity(2), 0.0_dp, 0.0_dp, 'disabled species must not escape')
    call assert_close_dp(species_ledger%residual(), 0.0_dp, 1.0e-12_dp, 'species boundary ledger residual mismatch')
  end subroutine test_species_z_high_reflect_preserves_ambient_escape

  subroutine test_large_finite_charge_norm_stays_finite()
    type(mesh_type) :: large_mesh
    type(app_config) :: large_cfg
    type(sim_stats) :: large_stats
    real(dp) :: tri_v0(3, 3), tri_v1(3, 3), tri_v2(3, 3)
    real(dp) :: charge_scale, expected_charge(3)
    integer(i32) :: species_idx

    do species_idx = 1_i32, 3_i32
      tri_v0(:, species_idx) = [2.0_dp*real(species_idx - 1_i32, dp), 0.0_dp, 0.0_dp]
      tri_v1(:, species_idx) = tri_v0(:, species_idx) + [1.0_dp, 0.0_dp, 0.0_dp]
      tri_v2(:, species_idx) = tri_v0(:, species_idx) + [0.0_dp, 1.0_dp, 0.0_dp]
    end do
    call init_mesh(large_mesh, tri_v0, tri_v1, tri_v2)
    large_mesh%elem_vacuum_sign = 1_i32
    large_mesh%vacuum_normals = large_mesh%normals

    charge_scale = huge(1.0_dp)/1.5_dp
    expected_charge = [charge_scale, charge_scale, 0.5_dp*charge_scale]
    call default_app_config(large_cfg)
    large_cfg%sim%batch_count = 1_i32
    large_cfg%sim%dt = 1.0_dp
    large_cfg%sim%max_step = 1_i32
    large_cfg%sim%q_floor = 1.0e-30_dp
    large_cfg%n_particle_species = 3_i32
    do species_idx = 1_i32, large_cfg%n_particle_species
      large_cfg%particle_species(species_idx) = species_from_defaults()
      large_cfg%particle_species(species_idx)%source_mode = 'volume_seed'
      large_cfg%particle_species(species_idx)%npcls_per_step = 1_i32
      large_cfg%particle_species(species_idx)%q_particle = expected_charge(species_idx)
      large_cfg%particle_species(species_idx)%m_particle = charge_scale
      large_cfg%particle_species(species_idx)%w_particle = 1.0_dp
      large_cfg%particle_species(species_idx)%pos_low = &
        [2.0_dp*real(species_idx - 1_i32, dp) + 0.2_dp, 0.2_dp, 0.5_dp]
      large_cfg%particle_species(species_idx)%pos_high = large_cfg%particle_species(species_idx)%pos_low
      large_cfg%particle_species(species_idx)%drift_velocity = [0.0_dp, 0.0_dp, -1.0_dp]
      large_cfg%particle_species(species_idx)%temperature_k = 0.0_dp
    end do

    call seed_particles_from_config(large_cfg)
    call run_absorption_insulator(large_mesh, large_cfg, large_stats)
    call assert_equal_i64(large_stats%absorbed, 3_i64, 'large finite charge fixture must absorb every particle')
    call assert_true(all(large_mesh%q_elem < huge(1.0_dp)), 'large finite committed charges must remain finite')
    do species_idx = 1_i32, large_cfg%n_particle_species
      call assert_close_dp( &
        large_mesh%q_elem(species_idx)/charge_scale, expected_charge(species_idx)/charge_scale, &
        1.0e-12_dp, 'large finite committed charge mismatch' &
        )
    end do
    call assert_close_dp( &
      large_stats%last_rel_change, 1.0_dp, 1.0e-12_dp, &
      'scaled charge norm must preserve the relative change' &
      )
  end subroutine test_large_finite_charge_norm_stays_finite
end program test_simulator_loop

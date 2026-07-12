!> 吸着ループの統計更新・堆積・履歴出力を検証する統合テスト。
program test_simulator
  use bem_kinds, only: dp, i32, i64
  use bem_mesh, only: init_mesh, prepare_periodic2_collision_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config
  use bem_types, only: mesh_type, sim_stats, bc_open, bc_reflect, bc_periodic
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
  real(dp) :: potential_value
  integer :: u, ios
  integer(i32) :: potential_batch_idx, potential_elem_idx
  character(len=256) :: line
  integer(i32) :: n_lines
  character(len=*), parameter :: history_path = 'test_simulator_history_tmp.csv'
  character(len=*), parameter :: potential_history_path = 'test_simulator_potential_history_tmp.csv'
  character(len=*), parameter :: collision_failure_path = 'test_simulator_collision_failure_tmp.log'
  character(len=*), parameter :: photo_collision_failure_path = 'test_simulator_photo_collision_failure_tmp.log'
  character(len=*), parameter :: box_event_failure_path = 'test_simulator_box_event_failure_tmp.log'
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

  call test_init(13)

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

  call test_begin('photoelectron_individual_return_histogram')
  call test_photoelectron_individual_return_histogram()
  call test_end()

  call test_begin('split_outer_instant_return_ledger')
  call test_split_outer_instant_return_ledger()
  call test_end()

  call test_begin('unified_explicit_outer_escape_ledger')
  call test_unified_explicit_outer_escape_ledger()
  call test_end()

  call test_begin('unified_interface_height_surface_charge_invariance')
  call test_unified_interface_height_surface_charge_invariance()
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

  call test_begin('multiple_box_event_failure_context')
  call test_multiple_box_event_failure_context()
  call test_end()

  call test_summary()

contains

  subroutine test_photoelectron_individual_return_histogram()
    use bem_constants, only: eps0
    use bem_outer_plasma_photoelectron, only: photoelectron_coupling_state_type
    type(mesh_type) :: photo_mesh
    type(app_config) :: photo_cfg
    type(sim_stats) :: photo_stats
    type(charge_ledger_type) :: photo_ledger
    type(photoelectron_coupling_state_type) :: photo_state
    real(dp) :: tri_v0(3, 2), tri_v1(3, 2), tri_v2(3, 2), panel_charge

    tri_v0(:, 1) = [0.0_dp, 0.0_dp, 0.999_dp]
    tri_v1(:, 1) = [1.0_dp, 0.0_dp, 0.999_dp]
    tri_v2(:, 1) = [1.0_dp, 1.0_dp, 0.999_dp]
    tri_v0(:, 2) = [0.0_dp, 0.0_dp, 0.999_dp]
    tri_v1(:, 2) = [1.0_dp, 1.0_dp, 0.999_dp]
    tri_v2(:, 2) = [0.0_dp, 1.0_dp, 0.999_dp]
    panel_charge = 2.5_dp*eps0
    call init_mesh(photo_mesh, tri_v0, tri_v1, tri_v2, q0=[panel_charge, panel_charge])
    photo_mesh%elem_vacuum_sign = 1_i32
    photo_mesh%vacuum_normals = photo_mesh%normals

    call default_app_config(photo_cfg)
    photo_cfg%sim%rng_seed = 999_i32
    photo_cfg%sim%batch_count = 1_i32
    photo_cfg%sim%batch_duration = 1.0_dp
    photo_cfg%sim%has_batch_duration = .true.
    photo_cfg%sim%dt = 0.01_dp
    photo_cfg%sim%max_step = 1_i32
    photo_cfg%sim%softening = 0.0_dp
    photo_cfg%sim%field_solver = 'direct'
    photo_cfg%sim%field_bc_mode = 'periodic2'
    photo_cfg%sim%use_box = .true.
    photo_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    photo_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    photo_cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    photo_cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    photo_cfg%field%backend = 'direct'
    photo_cfg%panel%source_model = 'triangle_p0'
    photo_cfg%panel%kernel_id = 'triangle_p0_exact_direct'
    photo_cfg%panel%surface_side_policy = 'per_element'
    photo_cfg%periodic2%nonzero_mode_backend = 'panel_spectral_reference'
    photo_cfg%periodic2%zero_mode_policy = 'exclude_k0'
    photo_cfg%periodic2%lower_boundary_model = 'e_bottom_zero'
    photo_cfg%periodic2%reference_mode_layers = 2_i32
    photo_cfg%periodic2%panel_quadrature_order = 12_i32
    photo_cfg%periodic2%interface_phi_tolerance = 1.0e6_dp
    photo_cfg%periodic2%interface_field_tolerance = 1.0e6_dp
    photo_cfg%outer_plasma%model = 'linear_debye'
    photo_cfg%outer_plasma%return_model = 'electrostatic_1d_instant_return'
    photo_cfg%outer_plasma%photoelectron_closure = 'individual_return'
    photo_cfg%outer_plasma%interface_z = 1.0_dp
    photo_cfg%outer_plasma%debye_length = 0.2_dp
    photo_cfg%outer_plasma%thermal_voltage = 10.0_dp
    photo_cfg%outer_plasma%max_linearity_ratio = 1.0_dp
    photo_cfg%outer_plasma%photoelectron_histogram_bins = 4_i32
    photo_cfg%outer_plasma%photoelectron_histogram_energy_max = 10.0_dp
    photo_cfg%outer_plasma%photoelectron_ambient_charge_scale = 100.0_dp
    photo_cfg%outer_plasma%max_photoelectron_charge_ratio = 0.1_dp
    photo_cfg%coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
    photo_cfg%coupling%field_evolution_timescale = 10.0_dp
    photo_cfg%coupling%max_frozen_field_ratio = 1.0_dp
    photo_cfg%n_particle_species = 1_i32
    photo_cfg%particle_species(1) = species_from_defaults()
    photo_cfg%particle_species(1)%source_mode = 'photo_raycast'
    photo_cfg%particle_species(1)%rays_per_batch = 1_i32
    photo_cfg%particle_species(1)%emit_current_density_a_m2 = 1.0_dp
    photo_cfg%particle_species(1)%q_particle = -1.0_dp
    photo_cfg%particle_species(1)%m_particle = 1.0_dp
    photo_cfg%particle_species(1)%inject_face = 'z_high'
    photo_cfg%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, 1.0_dp]
    photo_cfg%particle_species(1)%pos_high = [1.0_dp, 1.0_dp, 1.0_dp]
    photo_cfg%particle_species(1)%ray_direction = [0.0_dp, 0.0_dp, -1.0_dp]
    photo_cfg%particle_species(1)%has_ray_direction = .true.
    photo_cfg%particle_species(1)%temperature_k = 0.0_dp
    photo_cfg%particle_species(1)%normal_drift_speed = 0.2_dp
    photo_cfg%particle_species(1)%deposit_opposite_charge_on_emit = .true.
    photo_cfg%particle_species(1)%has_deposit_opposite_charge_on_emit = .true.

    call prepare_periodic2_collision_mesh(photo_mesh, photo_cfg%sim)
    call seed_particles_from_config(photo_cfg)
    call run_absorption_insulator( &
      photo_mesh, photo_cfg, photo_stats, charge_ledger=photo_ledger, photoelectron_state=photo_state &
      )
    call assert_equal_i64(photo_state%cumulative%total_count(), 1_i64, 'photoelectron histogram count mismatch')
    call assert_close_dp(photo_state%cumulative%total_signed_charge(), -1.0_dp, 1.0e-12_dp, 'histogram charge mismatch')
    call assert_close_dp(sum(photo_ledger%emitted_from_surface), -1.0_dp, 1.0e-12_dp, 'emission ledger mismatch')
    call assert_close_dp(sum(photo_ledger%interface_returned_gross), -1.0_dp, 1.0e-12_dp, 'return ledger mismatch')
    call assert_close_dp(photo_ledger%residual(), 0.0_dp, 1.0e-12_dp, 'photoelectron transaction mismatch')
  end subroutine test_photoelectron_individual_return_histogram

  subroutine test_split_outer_instant_return_ledger()
    use bem_constants, only: eps0
    use bem_electrostatic_snapshot, only: electrostatic_diagnostics_type
    type(mesh_type) :: split_mesh
    type(app_config) :: split_cfg
    type(sim_stats) :: split_stats
    type(charge_ledger_type) :: split_ledger
    type(electrostatic_diagnostics_type) :: split_diagnostics
    real(dp) :: tri_v0(3, 2), tri_v1(3, 2), tri_v2(3, 2), panel_charge

    tri_v0(:, 1) = [0.0_dp, 0.0_dp, 0.25_dp]
    tri_v1(:, 1) = [1.0_dp, 0.0_dp, 0.25_dp]
    tri_v2(:, 1) = [1.0_dp, 1.0_dp, 0.25_dp]
    tri_v0(:, 2) = [0.0_dp, 0.0_dp, 0.25_dp]
    tri_v1(:, 2) = [1.0_dp, 1.0_dp, 0.25_dp]
    tri_v2(:, 2) = [0.0_dp, 1.0_dp, 0.25_dp]
    panel_charge = -2.5_dp*eps0
    call init_mesh(split_mesh, tri_v0, tri_v1, tri_v2, q0=[panel_charge, panel_charge])
    split_mesh%elem_vacuum_sign = 1_i32
    split_mesh%vacuum_normals = split_mesh%normals

    call default_app_config(split_cfg)
    split_cfg%sim%rng_seed = 998_i32
    split_cfg%sim%batch_count = 1_i32
    split_cfg%sim%dt = 0.01_dp
    split_cfg%sim%max_step = 1_i32
    split_cfg%sim%softening = 0.0_dp
    split_cfg%sim%field_solver = 'direct'
    split_cfg%sim%field_bc_mode = 'periodic2'
    split_cfg%sim%use_box = .true.
    split_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    split_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    split_cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    split_cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    split_cfg%field%backend = 'direct'
    split_cfg%panel%source_model = 'triangle_p0'
    split_cfg%panel%kernel_id = 'triangle_p0_exact_direct'
    split_cfg%panel%surface_side_policy = 'per_element'
    split_cfg%periodic2%nonzero_mode_backend = 'panel_spectral_reference'
    split_cfg%periodic2%zero_mode_policy = 'exclude_k0'
    split_cfg%periodic2%lower_boundary_model = 'e_bottom_zero'
    split_cfg%periodic2%reference_mode_layers = 2_i32
    split_cfg%periodic2%panel_quadrature_order = 12_i32
    split_cfg%periodic2%interface_phi_tolerance = 1.0e6_dp
    split_cfg%periodic2%interface_field_tolerance = 1.0e6_dp
    split_cfg%outer_plasma%model = 'linear_debye'
    split_cfg%outer_plasma%return_model = 'electrostatic_1d_instant_return'
    split_cfg%outer_plasma%interface_z = 1.0_dp
    split_cfg%outer_plasma%debye_length = 0.2_dp
    split_cfg%outer_plasma%thermal_voltage = 10.0_dp
    split_cfg%outer_plasma%max_linearity_ratio = 1.0_dp
    split_cfg%coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
    split_cfg%coupling%field_evolution_timescale = 10.0_dp
    split_cfg%coupling%max_frozen_field_ratio = 1.0_dp
    split_cfg%n_particle_species = 1_i32
    split_cfg%particle_species(1) = species_from_defaults()
    split_cfg%particle_species(1)%source_mode = 'volume_seed'
    split_cfg%particle_species(1)%npcls_per_step = 1_i32
    split_cfg%particle_species(1)%q_particle = 1.0_dp
    split_cfg%particle_species(1)%m_particle = 1.0_dp
    split_cfg%particle_species(1)%w_particle = 1.0_dp
    split_cfg%particle_species(1)%pos_low = [0.5_dp, 0.5_dp, 0.999_dp]
    split_cfg%particle_species(1)%pos_high = split_cfg%particle_species(1)%pos_low
    split_cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, 0.2_dp]
    split_cfg%particle_species(1)%temperature_k = 0.0_dp

    call prepare_periodic2_collision_mesh(split_mesh, split_cfg%sim)
    call seed_particles_from_config(split_cfg)
    call run_absorption_insulator( &
      split_mesh, split_cfg, split_stats, charge_ledger=split_ledger, &
      electrostatic_diagnostics=split_diagnostics &
      )
    call assert_equal_i64(split_stats%escaped_boundary, 0_i64, 'returned particle must not escape')
    call assert_equal_i64(split_stats%survived_max_step, 1_i64, 'returned particle should remain local after one step')
    call assert_close_dp(split_ledger%interface_outward_gross(1), 1.0_dp, 1.0e-12_dp, 'outward gross mismatch')
    call assert_close_dp(split_ledger%interface_returned_gross(1), 1.0_dp, 1.0e-12_dp, 'returned gross mismatch')
    call assert_true(split_diagnostics%max_outer_flight_time > 0.0_dp, 'outer flight-time diagnostic is missing')
    call assert_true(split_diagnostics%max_frozen_field_ratio > 0.0_dp, 'frozen-field diagnostic is missing')
  end subroutine test_split_outer_instant_return_ledger

  subroutine test_unified_explicit_outer_escape_ledger()
    use bem_electrostatic_snapshot, only: electrostatic_diagnostics_type
    type(mesh_type) :: unified_mesh
    type(app_config) :: unified_cfg
    type(sim_stats) :: unified_stats
    type(charge_ledger_type) :: unified_ledger
    type(electrostatic_diagnostics_type) :: unified_diagnostics
    real(dp) :: tri_v0(3, 2), tri_v1(3, 2), tri_v2(3, 2)

    tri_v0(:, 1) = [0.0_dp, 0.0_dp, 0.25_dp]
    tri_v1(:, 1) = [1.0_dp, 0.0_dp, 0.25_dp]
    tri_v2(:, 1) = [1.0_dp, 1.0_dp, 0.25_dp]
    tri_v0(:, 2) = [0.0_dp, 0.0_dp, 0.25_dp]
    tri_v1(:, 2) = [1.0_dp, 1.0_dp, 0.25_dp]
    tri_v2(:, 2) = [0.0_dp, 1.0_dp, 0.25_dp]
    call init_mesh(unified_mesh, tri_v0, tri_v1, tri_v2)
    unified_mesh%elem_vacuum_sign = 1_i32
    unified_mesh%vacuum_normals = unified_mesh%normals

    call default_app_config(unified_cfg)
    unified_cfg%sim%rng_seed = 999_i32
    unified_cfg%sim%batch_count = 1_i32
    unified_cfg%sim%dt = 0.2_dp
    unified_cfg%sim%max_step = 1_i32
    unified_cfg%sim%softening = 0.0_dp
    unified_cfg%sim%q_floor = 1.0e-40_dp
    unified_cfg%sim%field_solver = 'direct'
    unified_cfg%sim%field_bc_mode = 'periodic2'
    unified_cfg%sim%use_box = .true.
    unified_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    unified_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    unified_cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    unified_cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    unified_cfg%field%backend = 'direct'
    unified_cfg%panel%source_model = 'triangle_p0'
    unified_cfg%panel%kernel_id = 'triangle_p0_exact_direct'
    unified_cfg%panel%surface_side_policy = 'per_element'
    unified_cfg%periodic2%nonzero_mode_backend = 'panel_spectral_reference'
    unified_cfg%periodic2%zero_mode_policy = 'exclude_k0'
    unified_cfg%periodic2%lower_boundary_model = 'e_bottom_zero'
    unified_cfg%periodic2%reference_mode_layers = 2_i32
    unified_cfg%periodic2%panel_quadrature_order = 8_i32
    unified_cfg%periodic2%interface_sample_n = 4_i32
    unified_cfg%outer_plasma%model = 'unified_linear_response'
    unified_cfg%outer_plasma%return_model = 'electrostatic_3d_explicit_orbit'
    unified_cfg%outer_plasma%interface_z = 1.0_dp
    unified_cfg%outer_plasma%debye_length = 0.2_dp
    unified_cfg%outer_plasma%thermal_voltage = 10.0_dp
    unified_cfg%outer_plasma%unified_grid_points = 33_i32
    unified_cfg%outer_plasma%max_linearity_ratio = 1.0_dp
    unified_cfg%coupling%particle_transfer_mode = 'electrostatic_3d_explicit_orbit'
    unified_cfg%coupling%field_evolution_timescale = 100.0_dp
    unified_cfg%coupling%max_frozen_field_ratio = 1.0_dp
    unified_cfg%coupling%outer_orbit_dt = 0.1_dp
    unified_cfg%coupling%outer_orbit_max_steps = 100_i32
    unified_cfg%coupling%outer_orbit_energy_tolerance = 1.0e-8_dp
    unified_cfg%n_particle_species = 1_i32
    unified_cfg%particle_species(1) = species_from_defaults()
    unified_cfg%particle_species(1)%source_mode = 'volume_seed'
    unified_cfg%particle_species(1)%npcls_per_step = 1_i32
    unified_cfg%particle_species(1)%q_particle = 1.0e-30_dp
    unified_cfg%particle_species(1)%m_particle = 1.0_dp
    unified_cfg%particle_species(1)%w_particle = 1.0_dp
    unified_cfg%particle_species(1)%pos_low = [0.5_dp, 0.5_dp, 0.9_dp]
    unified_cfg%particle_species(1)%pos_high = unified_cfg%particle_species(1)%pos_low
    unified_cfg%particle_species(1)%drift_velocity = [0.2_dp, 0.1_dp, 1.0_dp]
    unified_cfg%particle_species(1)%temperature_k = 0.0_dp

    call prepare_periodic2_collision_mesh(unified_mesh, unified_cfg%sim)
    call seed_particles_from_config(unified_cfg)
    call run_absorption_insulator( &
      unified_mesh, unified_cfg, unified_stats, charge_ledger=unified_ledger, &
      electrostatic_diagnostics=unified_diagnostics &
      )
    call assert_equal_i64(unified_stats%escaped_boundary, 1_i64, 'explicit outer particle must escape')
    call assert_close_dp( &
      unified_ledger%interface_outward_gross(1), 1.0e-30_dp, 1.0e-40_dp, 'explicit outward gross mismatch' &
      )
    call assert_close_dp( &
      unified_ledger%escaped_to_infinity(1), 1.0e-30_dp, 1.0e-40_dp, 'explicit escape ledger mismatch' &
      )
    call assert_close_dp(unified_ledger%residual(), 0.0_dp, 1.0e-40_dp, 'explicit orbit ledger residual mismatch')
    call assert_true( &
      unified_diagnostics%max_outer_energy_relative_error <= unified_cfg%coupling%outer_orbit_energy_tolerance, &
      'explicit outer energy diagnostic must satisfy its configured contract' &
      )
  end subroutine test_unified_explicit_outer_escape_ledger

  subroutine test_unified_interface_height_surface_charge_invariance()
    type(mesh_type) :: low_mesh, high_mesh
    type(sim_stats) :: low_stats, high_stats
    type(charge_ledger_type) :: low_ledger, high_ledger

    call run_interface_height_fixture(0.55_dp, low_mesh, low_stats, low_ledger)
    call run_interface_height_fixture(0.80_dp, high_mesh, high_stats, high_ledger)

    call assert_true( &
      low_ledger%interface_outward_gross(1) > 0.0_dp, &
      'low interface must transfer the turning particle to the outer orbit' &
      )
    call assert_close_dp( &
      high_ledger%interface_outward_gross(1), 0.0_dp, 0.0_dp, &
      'high interface must keep the turning particle in the local domain' &
      )
    call assert_close_dp( &
      low_ledger%interface_returned_gross(1), low_ledger%interface_outward_gross(1), 1.0e-24_dp, &
      'low-interface outer transaction must return in full' &
      )
    call assert_equal_i64(low_stats%absorbed, 1_i64, 'low-interface particle must be absorbed')
    call assert_equal_i64(high_stats%absorbed, low_stats%absorbed, 'absorption must be interface-height invariant')
    call assert_equal_i64(low_stats%escaped, 0_i64, 'low-interface particle must not escape')
    call assert_equal_i64(high_stats%escaped, low_stats%escaped, 'escape count must be interface-height invariant')
    call assert_close_dp( &
      sum(high_mesh%q_elem), sum(low_mesh%q_elem), 1.0e-24_dp, &
      'final surface charge must be interface-height invariant' &
      )
    call assert_close_dp(low_ledger%residual(), 0.0_dp, 1.0e-24_dp, 'low-interface ledger must close')
    call assert_close_dp(high_ledger%residual(), 0.0_dp, 1.0e-24_dp, 'high-interface ledger must close')
  end subroutine test_unified_interface_height_surface_charge_invariance

  subroutine run_interface_height_fixture(interface_z, fixture_mesh, fixture_stats, fixture_ledger)
    use bem_constants, only: eps0
    real(dp), intent(in) :: interface_z
    type(mesh_type), intent(out) :: fixture_mesh
    type(sim_stats), intent(out) :: fixture_stats
    type(charge_ledger_type), intent(out) :: fixture_ledger
    type(app_config) :: fixture_cfg
    real(dp) :: tri_v0(3, 2), tri_v1(3, 2), tri_v2(3, 2), panel_charge

    tri_v0(:, 1) = [0.0_dp, 0.0_dp, 0.25_dp]
    tri_v1(:, 1) = [1.0_dp, 0.0_dp, 0.25_dp]
    tri_v2(:, 1) = [1.0_dp, 1.0_dp, 0.25_dp]
    tri_v0(:, 2) = [0.0_dp, 0.0_dp, 0.25_dp]
    tri_v1(:, 2) = [1.0_dp, 1.0_dp, 0.25_dp]
    tri_v2(:, 2) = [0.0_dp, 1.0_dp, 0.25_dp]
    panel_charge = -2.5_dp*eps0
    call init_mesh(fixture_mesh, tri_v0, tri_v1, tri_v2, q0=[panel_charge, panel_charge])
    fixture_mesh%elem_vacuum_sign = 1_i32
    fixture_mesh%vacuum_normals = fixture_mesh%normals

    call default_app_config(fixture_cfg)
    fixture_cfg%sim%rng_seed = 1001_i32
    fixture_cfg%sim%batch_count = 1_i32
    fixture_cfg%sim%dt = 0.01_dp
    fixture_cfg%sim%max_step = 300_i32
    fixture_cfg%sim%softening = 0.0_dp
    fixture_cfg%sim%q_floor = 1.0e-40_dp
    fixture_cfg%sim%field_solver = 'direct'
    fixture_cfg%sim%field_bc_mode = 'periodic2'
    fixture_cfg%sim%use_box = .true.
    fixture_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    fixture_cfg%sim%box_max = [1.0_dp, 1.0_dp, interface_z]
    fixture_cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    fixture_cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    fixture_cfg%field%backend = 'direct'
    fixture_cfg%panel%source_model = 'triangle_p0'
    fixture_cfg%panel%kernel_id = 'triangle_p0_exact_direct'
    fixture_cfg%panel%surface_side_policy = 'per_element'
    fixture_cfg%periodic2%nonzero_mode_backend = 'panel_spectral_reference'
    fixture_cfg%periodic2%zero_mode_policy = 'exclude_k0'
    fixture_cfg%periodic2%lower_boundary_model = 'e_bottom_zero'
    fixture_cfg%periodic2%reference_mode_layers = 2_i32
    fixture_cfg%periodic2%panel_quadrature_order = 8_i32
    fixture_cfg%periodic2%interface_sample_n = 4_i32
    fixture_cfg%outer_plasma%model = 'unified_linear_response'
    fixture_cfg%outer_plasma%return_model = 'electrostatic_3d_explicit_orbit'
    fixture_cfg%outer_plasma%interface_z = interface_z
    fixture_cfg%outer_plasma%debye_length = 0.2_dp
    fixture_cfg%outer_plasma%thermal_voltage = 10.0_dp
    fixture_cfg%outer_plasma%unified_grid_points = 33_i32
    fixture_cfg%outer_plasma%max_linearity_ratio = 1.0_dp
    fixture_cfg%coupling%particle_transfer_mode = 'electrostatic_3d_explicit_orbit'
    fixture_cfg%coupling%field_evolution_timescale = 100.0_dp
    fixture_cfg%coupling%max_frozen_field_ratio = 1.0_dp
    fixture_cfg%coupling%outer_orbit_dt = 0.005_dp
    fixture_cfg%coupling%outer_orbit_max_steps = 10000_i32
    fixture_cfg%coupling%outer_orbit_energy_tolerance = 1.0e-3_dp
    fixture_cfg%n_particle_species = 1_i32
    fixture_cfg%particle_species(1) = species_from_defaults()
    fixture_cfg%particle_species(1)%source_mode = 'volume_seed'
    fixture_cfg%particle_species(1)%npcls_per_step = 1_i32
    fixture_cfg%particle_species(1)%q_particle = 1.0e-12_dp
    fixture_cfg%particle_species(1)%m_particle = 1.0e-12_dp
    fixture_cfg%particle_species(1)%w_particle = 1.0_dp
    fixture_cfg%particle_species(1)%pos_low = [0.3_dp, 0.4_dp, 0.50_dp]
    fixture_cfg%particle_species(1)%pos_high = fixture_cfg%particle_species(1)%pos_low
    fixture_cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, 0.5_dp]
    fixture_cfg%particle_species(1)%temperature_k = 0.0_dp

    call prepare_periodic2_collision_mesh(fixture_mesh, fixture_cfg%sim)
    call seed_particles_from_config(fixture_cfg)
    call run_absorption_insulator( &
      fixture_mesh, fixture_cfg, fixture_stats, charge_ledger=fixture_ledger &
      )
  end subroutine run_interface_height_fixture

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

  subroutine test_reflected_remainder_deposits()
    type(mesh_type) :: reflected_mesh
    type(app_config) :: reflected_cfg
    type(sim_stats) :: reflected_stats
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1)

    tri_v0(:, 1) = [0.5_dp, -1.0_dp, -1.0_dp]
    tri_v1(:, 1) = [0.5_dp, 1.0_dp, -1.0_dp]
    tri_v2(:, 1) = [0.5_dp, 0.0_dp, 1.0_dp]
    call init_mesh(reflected_mesh, tri_v0, tri_v1, tri_v2)

    call default_app_config(reflected_cfg)
    reflected_cfg%sim%rng_seed = 993_i32
    reflected_cfg%sim%batch_count = 1_i32
    reflected_cfg%sim%dt = 1.0_dp
    reflected_cfg%sim%max_step = 1_i32
    reflected_cfg%sim%softening = 1.0e-6_dp
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
    failure_cfg%particle_species(1)%drift_velocity = [4.0_dp, 0.0_dp, 0.0_dp]
    failure_cfg%particle_species(1)%temperature_k = 0.0_dp
    call seed_particles_from_config(failure_cfg)
    call run_absorption_insulator(failure_mesh, failure_cfg, failure_stats)
  end subroutine run_multiple_box_event_failure_probe
end program test_simulator

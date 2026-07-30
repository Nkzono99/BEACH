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

  call test_init(18)

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

  call test_begin('kinetic_outer_profile_return_ledger')
  call test_kinetic_outer_profile_return_ledger()
  call test_end()

  call test_begin('implicit_mean_photoelectron_interface_transfer')
  call test_implicit_mean_photoelectron_interface_transfer()
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
    cfg_tree%sim, cfg_tree%field, cfg_tree%periodic2, cfg_tree%panel, cfg_tree%outer_plasma, cfg_tree%coupling &
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
    cfg_potential_history%panel, cfg_potential_history%outer_plasma, cfg_potential_history%coupling &
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

  call test_begin('zhao_steady_start_resume_does_not_reseed')
  call test_zhao_steady_start_resume()
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

  subroutine test_zhao_steady_start_resume()
    use bem_constants, only: qe
    use bem_electrostatic_snapshot, only: electrostatic_restart_state_type
    type(mesh_type) :: steady_mesh
    type(app_config) :: steady_cfg
    type(sim_stats) :: first_stats, resumed_stats
    type(injection_state) :: steady_injection_state
    type(electrostatic_restart_state_type) :: restart_state, saved_restart_state
    real(dp), allocatable :: saved_charge(:)
    real(dp), parameter :: electron_mass = 9.1093837139e-31_dp
    real(dp), parameter :: proton_mass = 1.67262192369e-27_dp
    real(dp), parameter :: box_width = 9.899494936611664e-5_dp
    real(dp), parameter :: interface_z = 2.0e-4_dp
    real(dp), parameter :: inward_drift = 4.0529988897111727e5_dp

    call default_app_config(steady_cfg)
    steady_cfg%sim%rng_seed = 24602_i32
    steady_cfg%sim%batch_count = 1_i32
    steady_cfg%sim%dt = 1.0e-12_dp
    steady_cfg%sim%batch_duration = 1.0e-12_dp
    steady_cfg%sim%has_batch_duration = .true.
    steady_cfg%sim%max_step = 100_i32
    steady_cfg%sim%q_floor = 1.0e-40_dp
    steady_cfg%sim%field_solver = 'direct'
    steady_cfg%sim%field_bc_mode = 'periodic2'
    steady_cfg%sim%use_box = .true.
    steady_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    steady_cfg%sim%box_max = [box_width, box_width, interface_z]
    steady_cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    steady_cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    steady_cfg%sim%reservoir_potential_model = 'none'
    steady_cfg%sim%open_boundary_model = 'escape'
    steady_cfg%sim%sheath_alpha_deg = 60.0_dp
    steady_cfg%sim%sheath_photoelectron_ref_density_cm3 = 0.0_dp
    steady_cfg%sim%sheath_electron_drift_mode = 'normal'
    steady_cfg%sim%sheath_ion_drift_mode = 'normal'

    steady_cfg%field%backend = 'direct'
    steady_cfg%panel%kernel_id = 'triangle_p0_exact_direct'
    steady_cfg%panel%surface_side_policy = 'per_element'
    steady_cfg%periodic2%nonzero_mode_backend = 'panel_spectral_reference'
    steady_cfg%periodic2%zero_mode_policy = 'exclude_k0'
    steady_cfg%periodic2%lower_boundary_model = 'symmetric_vacuum'
    steady_cfg%periodic2%reference_mode_layers = 2_i32
    steady_cfg%periodic2%panel_quadrature_order = 4_i32
    steady_cfg%periodic2%interface_sample_n = 2_i32
    steady_cfg%periodic2%interface_phi_tolerance = 1.0e-2_dp
    steady_cfg%periodic2%interface_field_tolerance = 1.0e-2_dp

    steady_cfg%outer_plasma%model = 'kinetic_1d'
    steady_cfg%outer_plasma%kinetic_closure = 'zhao_charge_driven'
    steady_cfg%outer_plasma%zhao_branch = 'auto'
    steady_cfg%outer_plasma%photoelectron_source_scale = 0.0_dp
    steady_cfg%outer_plasma%photoelectron_density_model = 'none'
    steady_cfg%outer_plasma%return_model = 'kinetic_1d_profile_return'
    steady_cfg%outer_plasma%interface_z = interface_z
    steady_cfg%outer_plasma%debye_length = 1.0_dp
    steady_cfg%outer_plasma%thermal_voltage = 10.0_dp
    steady_cfg%coupling%update_mode = 'explicit'
    steady_cfg%coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
    steady_cfg%coupling%steady_start_mode = 'zhao_floating'
    steady_cfg%coupling%steady_start_mesh_id = 1_i32
    steady_cfg%coupling%outer_update_stride = 1_i32
    steady_cfg%coupling%field_evolution_timescale = 1.0_dp
    steady_cfg%coupling%max_frozen_field_ratio = 1.0_dp
    steady_cfg%coupling%outer_queue_enabled = .false.

    steady_cfg%n_particle_species = 2_i32
    steady_cfg%particle_species(1) = species_from_defaults()
    steady_cfg%particle_species(1)%species_key = 'ambient_electron'
    steady_cfg%particle_species(1)%source_mode = 'reservoir_face'
    steady_cfg%particle_species(1)%inject_face = 'z_high'
    steady_cfg%particle_species(1)%q_particle = -qe
    steady_cfg%particle_species(1)%m_particle = electron_mass
    steady_cfg%particle_species(1)%number_density_m3 = 8.7e6_dp
    steady_cfg%particle_species(1)%has_number_density_m3 = .true.
    steady_cfg%particle_species(1)%temperature_ev = 12.0_dp
    steady_cfg%particle_species(1)%has_temperature_ev = .true.
    steady_cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -inward_drift]
    steady_cfg%particle_species(1)%target_macro_particles_per_batch = 1_i32
    steady_cfg%particle_species(1)%has_target_macro_particles_per_batch = .true.
    steady_cfg%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, interface_z]
    steady_cfg%particle_species(1)%pos_high = [box_width, box_width, interface_z]

    steady_cfg%particle_species(2) = species_from_defaults()
    steady_cfg%particle_species(2)%species_key = 'ambient_proton'
    steady_cfg%particle_species(2)%source_mode = 'reservoir_face'
    steady_cfg%particle_species(2)%inject_face = 'z_high'
    steady_cfg%particle_species(2)%q_particle = qe
    steady_cfg%particle_species(2)%m_particle = proton_mass
    steady_cfg%particle_species(2)%number_density_m3 = 8.7e6_dp
    steady_cfg%particle_species(2)%has_number_density_m3 = .true.
    steady_cfg%particle_species(2)%temperature_ev = 0.1_dp
    steady_cfg%particle_species(2)%has_temperature_ev = .true.
    steady_cfg%particle_species(2)%drift_velocity = [0.0_dp, 0.0_dp, -inward_drift]
    steady_cfg%particle_species(2)%target_macro_particles_per_batch = -1_i32
    steady_cfg%particle_species(2)%has_target_macro_particles_per_batch = .true.
    steady_cfg%particle_species(2)%pos_low = [0.0_dp, 0.0_dp, interface_z]
    steady_cfg%particle_species(2)%pos_high = [box_width, box_width, interface_z]

    steady_cfg%mesh_mode = 'template'
    steady_cfg%n_templates = 1_i32
    steady_cfg%templates(1)%enabled = .true.
    steady_cfg%templates(1)%kind = 'plane'
    steady_cfg%templates(1)%surface_side_policy = 'normal_plus'
    steady_cfg%templates(1)%size_x = box_width
    steady_cfg%templates(1)%size_y = box_width
    steady_cfg%templates(1)%nx = 2_i32
    steady_cfg%templates(1)%ny = 2_i32
    steady_cfg%templates(1)%center = [0.5_dp*box_width, 0.5_dp*box_width, 2.0e-6_dp]

    call build_mesh_from_config(steady_cfg, steady_mesh)
    call prepare_periodic2_collision_mesh(steady_mesh, steady_cfg%sim)
    call seed_particles_from_config(steady_cfg)
    allocate (steady_injection_state%macro_residual(steady_cfg%n_particle_species))
    steady_injection_state%macro_residual = 0.0_dp
    call run_absorption_insulator( &
      steady_mesh, steady_cfg, first_stats, inject_state=steady_injection_state, &
      electrostatic_restart_state=restart_state &
      )
    call assert_equal_i32(first_stats%batches, 1_i32, 'steady-start fixture did not complete its first batch')
    call assert_true(restart_state%outer_ready .and. restart_state%outer_profile_complete, &
                     'steady-start fixture did not export a complete outer profile')
    call assert_true(restart_state%outer_zhao_state_complete, &
                     'steady-start fixture did not export a complete Zhao state')
    call assert_true(restart_state%outer_zhao_source_scale_complete, &
                     'steady-start fixture did not export the resolved Zhao source scale')
    call assert_close_dp(restart_state%outer_photoelectron_source_scale, 0.0_dp, 0.0_dp, &
                         'no-photo steady-start fixture exported the wrong Zhao source scale')
    call assert_true(restart_state%outer_zhao_branch == 'C', &
                     'no-photo steady-start fixture did not select Zhao Type C')
    saved_charge = steady_mesh%q_elem
    saved_restart_state = restart_state

    call run_absorption_insulator( &
      steady_mesh, steady_cfg, resumed_stats, initial_stats=first_stats, &
      inject_state=steady_injection_state, electrostatic_restart_state=restart_state &
      )
    call assert_equal_i32(resumed_stats%batches, first_stats%batches, &
                          'zero-batch steady-start resume changed the completed batch count')
    call assert_true(all(steady_mesh%q_elem == saved_charge), &
                     'steady-start resume changed mesh charge or reapplied the seed')
    call assert_true(restart_state%outer_zhao_branch == saved_restart_state%outer_zhao_branch, &
                     'steady-start resume changed the Zhao branch')
    call assert_close_dp(restart_state%outer_interface_field, saved_restart_state%outer_interface_field, 0.0_dp, &
                         'steady-start resume changed the interface field')
    call assert_true(all(restart_state%outer_profile_z == saved_restart_state%outer_profile_z), &
                     'steady-start resume changed the outer-profile coordinates')
    call assert_true(all(restart_state%outer_profile_potential == saved_restart_state%outer_profile_potential), &
                     'steady-start resume changed the outer potential')
    call assert_true(all(restart_state%outer_profile_field == saved_restart_state%outer_profile_field), &
                     'steady-start resume changed the outer field')
    call assert_true( &
      all(restart_state%outer_profile_charge_density == saved_restart_state%outer_profile_charge_density), &
      'steady-start resume changed the outer charge density' &
      )
  end subroutine test_zhao_steady_start_resume

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

    call workspace%init(2_i32, 3_i32, 2_i32, implicit_mean_enabled=.true.)
    call regular_workspace%init(2_i32, 3_i32, 2_i32)
    call assert_equal_i32(int(size(workspace%dq_thread, 1), i32), 2_i32, 'workspace mesh capacity mismatch')
    call assert_equal_i32(int(size(workspace%dq_thread, 2), i32), 2_i32, 'workspace thread capacity mismatch')
    call assert_equal_i32( &
      int(size(workspace%ledger_charge_values), i32), 21_i32, 'workspace ledger charge capacity mismatch' &
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
      int(size(regular_workspace%mean_pending_charge), i32), 0_i32, &
      'regular workspace must not allocate implicit-mean element storage' &
      )
    call assert_equal_i32( &
      int(size(regular_workspace%mean_candidate_charge), i32), 0_i32, &
      'regular workspace must not allocate adaptive candidate storage' &
      )

    call workspace%prepare_particle_flags(5_i32, implicit_mean_enabled=.true.)
    call assert_equal_i32(int(size(workspace%outer_event_staging), i32), 0_i32, &
                          'disabled outer queue must not allocate per-particle staging')
    workspace%escaped_boundary_flag = .true.
    workspace%absorbed_flag = .true.
    workspace%absorbed_element = 1_i32
    workspace%soft_discarded_boundary_flag = .true.
    workspace%queued_outer_flag = .true.
    workspace%dq_thread = 1.0_dp
    workspace%photo_emission_dq = 2.0_dp
    workspace%interface_outward_thread = 3.0_dp
    workspace%interface_returned_thread = 4.0_dp
    workspace%interface_tau_max_thread = 5.0_dp
    workspace%interface_frozen_ratio_max_thread = 6.0_dp
    workspace%interface_energy_error_max_thread = 7.0_dp
    workspace%mean_pending_charge = 8.0_dp
    workspace%mean_deferred_source_charge = 9.0_dp
    workspace%mean_returned_destination_charge = 10.0_dp
    workspace%mean_candidate_charge = 11.0_dp
    workspace%neutral_return_charge_values = 12.0_dp
    workspace%neutral_return_terminal_counts = 13_i64
    workspace%neutral_return_weight_scale = 14.0_dp
    workspace%neutral_return_correction = 15.0_dp
    workspace%neutral_return_unresolved_fraction = 16.0_dp
    workspace%deferred_mean_return_element = 12_i32
    workspace%deferred_mean_terminal_absorbed = .true.
    workspace%deferred_mean_terminal_escaped = .true.

    call workspace%reset_before_injection()
    call workspace%prepare_particle_flags(2_i32, implicit_mean_enabled=.true.)
    call assert_true(all(workspace%dq_thread == 0.0_dp), 'workspace dq_thread reset mismatch')
    call assert_true(all(workspace%photo_emission_dq == 0.0_dp), 'workspace photo charge reset mismatch')
    call assert_true(all(workspace%interface_outward_thread == 0.0_dp), 'workspace outward reset mismatch')
    call assert_true(all(workspace%interface_returned_thread == 0.0_dp), 'workspace returned reset mismatch')
    call assert_true(all(workspace%interface_tau_max_thread == 0.0_dp), 'workspace tau reset mismatch')
    call assert_true(all(workspace%interface_frozen_ratio_max_thread == 0.0_dp), 'workspace ratio reset mismatch')
    call assert_true(all(workspace%interface_energy_error_max_thread == 0.0_dp), 'workspace energy reset mismatch')
    call assert_true(all(workspace%mean_pending_charge == 0.0_dp), 'workspace mean pending reset mismatch')
    call assert_true(all(workspace%mean_deferred_source_charge == 0.0_dp), 'workspace mean source reset mismatch')
    call assert_true(all(workspace%mean_returned_destination_charge == 0.0_dp), &
                     'workspace mean destination reset mismatch')
    call assert_true(all(workspace%mean_candidate_charge == 0.0_dp), 'workspace mean candidate reset mismatch')
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
    call assert_true(all(.not. workspace%queued_outer_flag(:2)), 'workspace outer queue flag reset mismatch')
    call assert_true(all(workspace%deferred_mean_return_element(:2) == -1_i32), &
                     'workspace mean return element reset mismatch')
    call assert_true(all(.not. workspace%deferred_mean_terminal_absorbed(:2)), &
                     'workspace mean absorbed staging reset mismatch')
    call assert_true(all(.not. workspace%deferred_mean_terminal_escaped(:2)), &
                     'workspace mean escaped staging reset mismatch')
    call assert_equal_i32( &
      int(size(workspace%escaped_boundary_flag), i32), 5_i32, 'workspace flags should retain grown capacity' &
      )

    call workspace%prepare_particle_flags(8_i32, implicit_mean_enabled=.true.)
    call assert_equal_i32( &
      int(size(workspace%escaped_boundary_flag), i32), 8_i32, 'workspace flags should grow on demand' &
      )
    call assert_true(all(.not. workspace%escaped_boundary_flag), 'grown workspace escaped flags must start clear')
    call assert_true(all(.not. workspace%queued_outer_flag), 'grown workspace outer queue flags must start clear')
    call assert_equal_i32(int(size(workspace%deferred_mean_return_element), i32), 8_i32, &
                          'workspace implicit mean staging should grow on demand')
    call assert_true(all(workspace%deferred_mean_return_element == -1_i32), &
                     'grown mean return elements must start clear')
    call workspace%prepare_particle_flags(8_i32, outer_queue_enabled=.true.)
    call assert_equal_i32(int(size(workspace%outer_event_staging), i32), 8_i32, &
                          'enabled outer queue staging must grow on demand')
  end subroutine test_batch_workspace_reuse

  subroutine test_kinetic_outer_profile_return_ledger()
    use bem_constants, only: eps0, qe
    use bem_electrostatic_snapshot, only: electrostatic_diagnostics_type, electrostatic_restart_state_type
    type(mesh_type) :: split_mesh
    type(app_config) :: split_cfg
    type(sim_stats) :: split_stats, split_resume_stats
    type(injection_state) :: split_injection_state
    type(charge_ledger_type) :: split_ledger
    type(electrostatic_diagnostics_type) :: split_diagnostics, split_resume_diagnostics
    type(electrostatic_restart_state_type) :: split_restart_state
    real(dp) :: tri_v0(3, 2), tri_v1(3, 2), tri_v2(3, 2), panel_charge
    real(dp) :: expected_maxima(3)
    real(dp), parameter :: electron_mass = 9.1093837139e-31_dp
    real(dp), parameter :: proton_mass = 1.67262192595e-27_dp

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
    split_cfg%sim%field_solver = 'direct'
    split_cfg%sim%field_bc_mode = 'periodic2'
    split_cfg%sim%use_box = .true.
    split_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    split_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    split_cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    split_cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    split_cfg%field%backend = 'direct'
    split_cfg%panel%kernel_id = 'triangle_p0_exact_direct'
    split_cfg%panel%surface_side_policy = 'per_element'
    split_cfg%periodic2%nonzero_mode_backend = 'panel_spectral_reference'
    split_cfg%periodic2%zero_mode_policy = 'exclude_k0'
    split_cfg%periodic2%lower_boundary_model = 'e_bottom_zero'
    split_cfg%periodic2%reference_mode_layers = 2_i32
    split_cfg%periodic2%panel_quadrature_order = 12_i32
    split_cfg%periodic2%interface_phi_tolerance = 1.0e6_dp
    split_cfg%periodic2%interface_field_tolerance = 1.0e6_dp
    split_cfg%outer_plasma%model = 'kinetic_1d'
    split_cfg%outer_plasma%kinetic_closure = 'absorbing_maxwellian'
    split_cfg%outer_plasma%return_model = 'kinetic_1d_profile_return'
    split_cfg%outer_plasma%interface_z = 1.0_dp
    split_cfg%outer_plasma%debye_length = 0.2_dp
    split_cfg%outer_plasma%thermal_voltage = 2.0_dp
    split_cfg%coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
    split_cfg%coupling%field_evolution_timescale = 10.0_dp
    split_cfg%coupling%max_frozen_field_ratio = 1.0_dp
    split_cfg%n_particle_species = 3_i32
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

    split_cfg%particle_species(2) = species_from_defaults()
    split_cfg%particle_species(2)%species_key = 'ambient_electron'
    split_cfg%particle_species(2)%source_mode = 'reservoir_face'
    split_cfg%particle_species(2)%inject_face = 'z_high'
    split_cfg%particle_species(2)%q_particle = -qe
    split_cfg%particle_species(2)%m_particle = electron_mass
    split_cfg%particle_species(2)%number_density_m3 = 1.0e6_dp
    split_cfg%particle_species(2)%has_number_density_m3 = .true.
    split_cfg%particle_species(2)%temperature_ev = 2.0_dp
    split_cfg%particle_species(2)%has_temperature_ev = .true.
    split_cfg%particle_species(2)%pos_low = [0.0_dp, 0.0_dp, 1.0_dp]
    split_cfg%particle_species(2)%pos_high = [1.0_dp, 1.0_dp, 1.0_dp]

    split_cfg%particle_species(3) = species_from_defaults()
    split_cfg%particle_species(3)%species_key = 'ambient_proton'
    split_cfg%particle_species(3)%source_mode = 'reservoir_face'
    split_cfg%particle_species(3)%inject_face = 'z_high'
    split_cfg%particle_species(3)%q_particle = qe
    split_cfg%particle_species(3)%m_particle = proton_mass
    split_cfg%particle_species(3)%number_density_m3 = 1.0e6_dp
    split_cfg%particle_species(3)%has_number_density_m3 = .true.
    split_cfg%particle_species(3)%temperature_ev = 0.0_dp
    split_cfg%particle_species(3)%has_temperature_ev = .true.
    split_cfg%particle_species(3)%drift_velocity = [ &
                                                   0.0_dp, 0.0_dp, &
                                                   -4.0_dp*sqrt(2.0_dp*qe/proton_mass) &
                                                   ]
    split_cfg%particle_species(3)%pos_low = [0.0_dp, 0.0_dp, 1.0_dp]
    split_cfg%particle_species(3)%pos_high = [1.0_dp, 1.0_dp, 1.0_dp]

    call prepare_periodic2_collision_mesh(split_mesh, split_cfg%sim)
    call seed_particles_from_config(split_cfg)
    allocate (split_injection_state%macro_residual(split_cfg%n_particle_species))
    split_injection_state%macro_residual = 0.0_dp
    call run_absorption_insulator( &
      split_mesh, split_cfg, split_stats, inject_state=split_injection_state, charge_ledger=split_ledger, &
      electrostatic_diagnostics=split_diagnostics, electrostatic_restart_state=split_restart_state &
      )
    call assert_equal_i64(split_stats%escaped_boundary, 0_i64, 'returned particle must not escape')
    call assert_equal_i64(split_stats%survived_max_step, 1_i64, 'returned particle should remain local after one step')
    call assert_close_dp(split_ledger%interface_outward_gross(1), 1.0_dp, 1.0e-12_dp, 'outward gross mismatch')
    call assert_close_dp(split_ledger%interface_returned_gross(1), 1.0_dp, 1.0e-12_dp, 'returned gross mismatch')
    call assert_true(split_diagnostics%max_outer_flight_time > 0.0_dp, 'outer flight-time diagnostic is missing')
    call assert_true(split_diagnostics%max_frozen_field_ratio > 0.0_dp, 'frozen-field diagnostic is missing')
    expected_maxima = [ &
                      split_diagnostics%max_outer_flight_time, split_diagnostics%max_frozen_field_ratio, &
                      split_diagnostics%max_outer_energy_relative_error &
                      ]
    call assert_true(split_restart_state%outer_max_diagnostics_complete, &
                     'exported cumulative outer diagnostics should be complete')
    call run_absorption_insulator( &
      split_mesh, split_cfg, split_resume_stats, initial_stats=split_stats, &
      inject_state=split_injection_state, electrostatic_diagnostics=split_resume_diagnostics, &
      electrostatic_restart_state=split_restart_state &
      )
    call assert_equal_i32(split_resume_stats%batches, split_stats%batches, &
                          'zero-batch resume changed the completed batch count')
    call assert_close_dp(split_resume_diagnostics%max_outer_flight_time, expected_maxima(1), 0.0_dp, &
                         'restart flight-time maximum was not preserved')
    call assert_close_dp(split_resume_diagnostics%max_frozen_field_ratio, expected_maxima(2), 0.0_dp, &
                         'restart frozen-field maximum was not preserved')
    call assert_close_dp(split_resume_diagnostics%max_outer_energy_relative_error, expected_maxima(3), 0.0_dp, &
                         'restart energy-error maximum was not preserved')
    call assert_close_dp(split_restart_state%max_outer_flight_time, expected_maxima(1), 0.0_dp, &
                         're-exported flight-time maximum mismatch')
    call assert_close_dp(split_restart_state%max_frozen_field_ratio, expected_maxima(2), 0.0_dp, &
                         're-exported frozen-field maximum mismatch')
    call assert_close_dp(split_restart_state%max_outer_energy_relative_error, expected_maxima(3), 0.0_dp, &
                         're-exported energy-error maximum mismatch')
  end subroutine test_kinetic_outer_profile_return_ledger

  subroutine test_implicit_mean_photoelectron_interface_transfer()
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use bem_constants, only: eps0, qe
    use bem_electrostatic_snapshot, only: electrostatic_diagnostics_type
    type(mesh_type) :: mean_mesh
    type(app_config) :: mean_cfg
    type(sim_stats) :: mean_stats
    type(injection_state) :: mean_injection_state
    type(charge_ledger_type) :: mean_ledger
    type(electrostatic_diagnostics_type) :: mean_diagnostics
    real(dp) :: tri_v0(3, 2), tri_v1(3, 2), tri_v2(3, 2)
    real(dp) :: initial_charge, charge_scale, ledger_residual, ledger_tolerance
    real(dp), parameter :: electron_mass = 9.1093837139e-31_dp
    real(dp), parameter :: proton_mass = 1.67262192595e-27_dp
    real(dp), parameter :: box_width = 1.0e-4_dp
    real(dp), parameter :: interface_z = 2.0e-4_dp
    real(dp), parameter :: support_z = 2.0e-6_dp
    real(dp), parameter :: debye_length = 0.1_dp
    real(dp), parameter :: return_interface_potential = 100.0_dp
    real(dp), parameter :: escape_interface_potential = -0.1_dp

    ! Cover z-low in this transaction fixture so every outer return has a
    ! physical surface destination.  Unsupported return-to-local-open escape
    ! is classified separately by test_external_step_driver.
    tri_v0(:, 1) = [0.0_dp, 0.0_dp, support_z]
    tri_v1(:, 1) = [box_width, 0.0_dp, support_z]
    tri_v2(:, 1) = [box_width, box_width, support_z]
    tri_v0(:, 2) = [0.0_dp, 0.0_dp, support_z]
    tri_v1(:, 2) = [box_width, box_width, support_z]
    tri_v2(:, 2) = [0.0_dp, box_width, support_z]
    initial_charge = eps0*box_width**2*return_interface_potential/debye_length
    call init_mesh( &
      mean_mesh, tri_v0, tri_v1, tri_v2, q0=[0.5_dp*initial_charge, 0.5_dp*initial_charge] &
      )
    mean_mesh%elem_vacuum_sign = 1_i32
    mean_mesh%vacuum_normals = mean_mesh%normals

    call default_app_config(mean_cfg)
    mean_cfg%mesh_mode = 'obj'
    mean_cfg%sim%rng_seed = 1001_i32
    mean_cfg%sim%batch_count = 2_i32
    mean_cfg%sim%dt = 3.0e-10_dp
    mean_cfg%sim%batch_duration = 3.0e-10_dp
    mean_cfg%sim%has_batch_duration = .true.
    mean_cfg%sim%max_step = 10_i32
    mean_cfg%sim%q_floor = 1.0e-40_dp
    mean_cfg%sim%field_solver = 'direct'
    mean_cfg%sim%field_bc_mode = 'periodic2'
    mean_cfg%sim%use_box = .true.
    mean_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    mean_cfg%sim%box_max = [box_width, box_width, interface_z]
    mean_cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    mean_cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]

    mean_cfg%field%backend = 'direct'
    mean_cfg%panel%kernel_id = 'triangle_p0_exact_direct'
    mean_cfg%panel%surface_side_policy = 'per_element'
    mean_cfg%periodic2%nonzero_mode_backend = 'panel_spectral_reference'
    mean_cfg%periodic2%zero_mode_policy = 'exclude_k0'
    mean_cfg%periodic2%lower_boundary_model = 'e_bottom_zero'
    mean_cfg%periodic2%reference_mode_layers = 2_i32
    mean_cfg%periodic2%panel_quadrature_order = 4_i32
    mean_cfg%periodic2%interface_sample_n = 2_i32
    mean_cfg%periodic2%interface_phi_tolerance = 1.0e6_dp
    mean_cfg%periodic2%interface_field_tolerance = 1.0e6_dp

    mean_cfg%outer_plasma%model = 'kinetic_1d'
    mean_cfg%outer_plasma%kinetic_closure = 'ambient_linear_debye'
    mean_cfg%outer_plasma%photoelectron_density_model = 'none'
    mean_cfg%outer_plasma%return_model = 'kinetic_1d_profile_return'
    mean_cfg%outer_plasma%interface_z = interface_z
    mean_cfg%outer_plasma%debye_length = debye_length
    mean_cfg%outer_plasma%thermal_voltage = 10.0_dp
    mean_cfg%coupling%update_mode = 'implicit_mean'
    mean_cfg%coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
    mean_cfg%coupling%outer_update_stride = 1_i32
    mean_cfg%coupling%field_evolution_timescale = 1.0_dp
    mean_cfg%coupling%max_frozen_field_ratio = 1.0_dp
    mean_cfg%coupling%outer_queue_enabled = .false.

    mean_cfg%n_particle_species = 4_i32
    mean_cfg%particle_species(1) = species_from_defaults()
    mean_cfg%particle_species(1)%species_key = 'ambient_electron'
    mean_cfg%particle_species(1)%source_mode = 'reservoir_face'
    mean_cfg%particle_species(1)%inject_face = 'z_high'
    mean_cfg%particle_species(1)%q_particle = -qe
    mean_cfg%particle_species(1)%m_particle = electron_mass
    mean_cfg%particle_species(1)%w_particle = 1.0_dp
    mean_cfg%particle_species(1)%number_density_m3 = 1.0e6_dp
    mean_cfg%particle_species(1)%has_number_density_m3 = .true.
    mean_cfg%particle_species(1)%temperature_ev = 2.0_dp
    mean_cfg%particle_species(1)%has_temperature_ev = .true.
    mean_cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -4.0e5_dp]
    mean_cfg%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, interface_z]
    mean_cfg%particle_species(1)%pos_high = [box_width, box_width, interface_z]

    mean_cfg%particle_species(2) = species_from_defaults()
    mean_cfg%particle_species(2)%species_key = 'ambient_proton'
    mean_cfg%particle_species(2)%source_mode = 'reservoir_face'
    mean_cfg%particle_species(2)%inject_face = 'z_high'
    mean_cfg%particle_species(2)%q_particle = qe
    mean_cfg%particle_species(2)%m_particle = proton_mass
    mean_cfg%particle_species(2)%w_particle = 1.0_dp
    mean_cfg%particle_species(2)%number_density_m3 = 1.0e6_dp
    mean_cfg%particle_species(2)%has_number_density_m3 = .true.
    mean_cfg%particle_species(2)%temperature_ev = 0.0_dp
    mean_cfg%particle_species(2)%has_temperature_ev = .true.
    mean_cfg%particle_species(2)%drift_velocity = [0.0_dp, 0.0_dp, -4.0e5_dp]
    mean_cfg%particle_species(2)%pos_low = [0.0_dp, 0.0_dp, interface_z]
    mean_cfg%particle_species(2)%pos_high = [box_width, box_width, interface_z]

    mean_cfg%particle_species(3) = species_from_defaults()
    mean_cfg%particle_species(3)%species_key = 'photoelectron'
    mean_cfg%particle_species(3)%source_mode = 'photo_raycast'
    mean_cfg%particle_species(3)%inject_face = 'z_high'
    mean_cfg%particle_species(3)%q_particle = -qe
    mean_cfg%particle_species(3)%m_particle = electron_mass
    mean_cfg%particle_species(3)%temperature_ev = 2.2_dp
    mean_cfg%particle_species(3)%has_temperature_ev = .true.
    mean_cfg%particle_species(3)%emit_current_density_a_m2 = 1.0e-6_dp
    mean_cfg%particle_species(3)%rays_per_batch = 4_i32
    mean_cfg%particle_species(3)%deposit_opposite_charge_on_emit = .true.
    mean_cfg%particle_species(3)%has_deposit_opposite_charge_on_emit = .true.
    mean_cfg%particle_species(3)%normal_drift_speed = 0.0_dp
    mean_cfg%particle_species(3)%pos_low = [0.25_dp*box_width, 0.25_dp*box_width, interface_z]
    mean_cfg%particle_species(3)%pos_high = [0.75_dp*box_width, 0.75_dp*box_width, interface_z]
    mean_cfg%particle_species(3)%ray_direction = [0.0_dp, 0.0_dp, -1.0_dp]
    mean_cfg%particle_species(3)%has_ray_direction = .true.

    ! A non-closure-owned volume source exercises the aggregate pending-current
    ! path.  Its absorbed positive charge must survive the final mean projection.
    mean_cfg%particle_species(4) = species_from_defaults()
    mean_cfg%particle_species(4)%species_key = 'extra_positive_seed'
    mean_cfg%particle_species(4)%source_mode = 'volume_seed'
    mean_cfg%particle_species(4)%q_particle = 1.0e-22_dp
    mean_cfg%particle_species(4)%m_particle = proton_mass
    mean_cfg%particle_species(4)%w_particle = 1.0_dp
    mean_cfg%particle_species(4)%npcls_per_step = 1_i32
    mean_cfg%particle_species(4)%temperature_ev = 0.0_dp
    mean_cfg%particle_species(4)%has_temperature_ev = .true.
    mean_cfg%particle_species(4)%drift_velocity = [0.0_dp, 0.0_dp, -1.0e4_dp]
    mean_cfg%particle_species(4)%pos_low = [0.5_dp*box_width, 0.5_dp*box_width, 3.0e-6_dp]
    mean_cfg%particle_species(4)%pos_high = mean_cfg%particle_species(4)%pos_low

    call prepare_periodic2_collision_mesh(mean_mesh, mean_cfg%sim)
    call seed_particles_from_config(mean_cfg)
    allocate (mean_injection_state%macro_residual(mean_cfg%n_particle_species))
    mean_injection_state%macro_residual = 0.0_dp
    call run_absorption_insulator( &
      mean_mesh, mean_cfg, mean_stats, inject_state=mean_injection_state, charge_ledger=mean_ledger, &
      electrostatic_diagnostics=mean_diagnostics &
      )

    call assert_equal_i32(mean_stats%batches, 2_i32, 'implicit-mean return fixture did not complete both batches')
    call assert_equal_i64( &
      mean_ledger%absorbed_count(3) + mean_ledger%escaped_count(3) + &
      mean_ledger%discarded_unresolved_count(3), mean_ledger%emitted_count(3), &
      'sub-threshold photoelectron ray samples did not terminate exactly once' &
      )
    call assert_true(mean_ledger%absorbed_count(3) > 0_i64, &
                     'sub-threshold fixture did not sample a surface return destination')
    call assert_true(mean_ledger%interface_outward_gross(3) < 0.0_dp, &
                     'photoelectron interface-outward charge was not recorded')
    call assert_true(mean_ledger%interface_returned_gross(3) < 0.0_dp, &
                     'energy-resolved photoelectron return was not recorded')
    call assert_close_dp( &
      mean_ledger%escaped_to_infinity(3), 0.0_dp, &
      4096.0_dp*epsilon(1.0_dp)*max(abs(mean_ledger%emitted_from_surface(3)), tiny(1.0_dp)), &
      'sub-threshold photoelectron escape charge must vanish')
    call assert_close_dp( &
      mean_ledger%emitted_from_surface(3), &
      mean_ledger%absorbed_on_surface(3) + mean_ledger%escaped_to_infinity(3) + &
      mean_ledger%discarded_unresolved(3), &
      1.0e-12_dp*abs(mean_ledger%emitted_from_surface(3)), &
      'sub-threshold photoelectron charge did not close after analytic weighting' &
      )
    call assert_equal_i64(mean_ledger%absorbed_count(4), 2_i64, &
                          'aggregate pending-current fixture did not absorb its extra source')
    call assert_close_dp( &
      mean_ledger%interface_returned_gross(3), mean_ledger%interface_outward_gross(3), &
      1.0e-12_dp*abs(mean_ledger%interface_outward_gross(3)), &
      'sub-threshold interface charge did not return completely' &
      )
    call assert_true(mean_diagnostics%max_outer_flight_time > 0.0_dp, &
                     'energy-resolved returns must report their outer flight time')
    call assert_true(mean_diagnostics%max_frozen_field_ratio > 0.0_dp, &
                     'energy-resolved returns must report their frozen-field ratio')
    call assert_true(mean_diagnostics%implicit_mean_shadow_diagnostics_available, &
                     'implicit-mean return fixture did not expose its shadow diagnostics')
    call assert_true(mean_diagnostics%implicit_mean_last_returned_outer_flight_time_mean > 0.0_dp, &
                     'returned photoelectrons must have a positive charge-weighted mean flight time')
    call assert_true( &
      mean_diagnostics%implicit_mean_last_returned_outer_flight_time_mean <= &
      mean_diagnostics%max_outer_flight_time, &
      'returned photoelectron mean flight time must not exceed the cumulative maximum' &
      )
    call assert_true( &
      mean_diagnostics%implicit_mean_last_returning_pe_column_charge_per_area > 0.0_dp, &
      'returned photoelectrons must produce a positive shadow column-charge estimate' &
      )
    call assert_true(all(ieee_is_finite(mean_mesh%q_elem)), 'implicit-mean surface charge is not finite')
    call assert_true(ieee_is_finite(mean_ledger%surface_charge_after), &
                     'implicit-mean ledger surface charge is not finite')
    call assert_close_dp(sum(mean_mesh%q_elem), mean_ledger%surface_charge_after, &
                         64.0_dp*epsilon(1.0_dp)*abs(mean_ledger%surface_charge_after), &
                         'implicit-mean mesh and ledger surface charges disagree')
    ledger_residual = mean_ledger%residual()
    charge_scale = max(abs(mean_ledger%surface_charge_before), abs(mean_ledger%surface_charge_after), tiny(1.0_dp))
    ledger_tolerance = max( &
                       1.0e-12_dp*charge_scale, &
                       16.0_dp*eps0*box_width**2/debye_length*1.0e-12_dp &
                       )
    call assert_true(ieee_is_finite(ledger_residual), 'implicit-mean charge-ledger residual is not finite')
    call assert_true(abs(ledger_residual) <= ledger_tolerance, &
                     'implicit-mean return charge ledger does not close')

    ! A negative interface potential accelerates outward electrons.  Both the
    ! scalar closure and every energy-resolved ray must therefore select escape.
    initial_charge = eps0*box_width**2*escape_interface_potential/debye_length
    mean_mesh%q_elem = [0.5_dp*initial_charge, 0.5_dp*initial_charge]
    mean_stats = sim_stats()
    mean_ledger = charge_ledger_type()
    mean_diagnostics = electrostatic_diagnostics_type()
    mean_injection_state%macro_residual = 0.0_dp
    call seed_particles_from_config(mean_cfg)
    call run_absorption_insulator( &
      mean_mesh, mean_cfg, mean_stats, inject_state=mean_injection_state, charge_ledger=mean_ledger, &
      electrostatic_diagnostics=mean_diagnostics &
      )

    call assert_equal_i32(mean_stats%batches, 2_i32, 'implicit-mean escape fixture did not complete both batches')
    call assert_equal_i64(mean_ledger%escaped_count(3), mean_ledger%emitted_count(3), &
                          'super-threshold photoelectrons did not all escape')
    call assert_equal_i64(mean_ledger%absorbed_count(3), 0_i64, &
                          'super-threshold photoelectrons incorrectly returned to the surface')
    call assert_close_dp(mean_ledger%interface_returned_gross(3), 0.0_dp, 0.0_dp, &
                         'super-threshold photoelectrons must not produce a return flux')
    call assert_true(mean_diagnostics%implicit_mean_shadow_diagnostics_available, &
                     'implicit-mean escape fixture did not expose its shadow diagnostics')
    call assert_close_dp( &
      mean_diagnostics%implicit_mean_last_returned_outer_flight_time_mean, 0.0_dp, 0.0_dp, &
      'escape-only photoelectrons must have zero returned mean flight time' &
      )
    call assert_close_dp( &
      mean_diagnostics%implicit_mean_last_returning_pe_column_charge_per_area, &
      0.0_dp, 0.0_dp, &
      'escape-only photoelectrons must have zero returning shadow column charge' &
      )
    call assert_close_dp( &
      mean_ledger%interface_outward_gross(3), mean_ledger%escaped_to_infinity(3), &
      1.0e-12_dp*abs(mean_ledger%interface_outward_gross(3)), &
      'super-threshold interface charge did not escape completely' &
      )
    call assert_close_dp( &
      mean_ledger%emitted_from_surface(3), &
      mean_ledger%absorbed_on_surface(3) + mean_ledger%escaped_to_infinity(3), &
      1.0e-12_dp*abs(mean_ledger%emitted_from_surface(3)), &
      'energy-resolved emitted charge did not terminate exactly once' &
      )
    call assert_true(all(ieee_is_finite(mean_mesh%q_elem)), 'implicit-mean escape charge is not finite')
    call assert_close_dp(sum(mean_mesh%q_elem), mean_ledger%surface_charge_after, &
                         64.0_dp*epsilon(1.0_dp)*abs(mean_ledger%surface_charge_after), &
                         'implicit-mean escape mesh and ledger charges disagree')
    ledger_residual = mean_ledger%residual()
    charge_scale = max(abs(mean_ledger%surface_charge_before), abs(mean_ledger%surface_charge_after), tiny(1.0_dp))
    ledger_tolerance = max( &
                       1.0e-12_dp*charge_scale, &
                       16.0_dp*eps0*box_width**2/debye_length*1.0e-12_dp &
                       )
    call assert_true(abs(ledger_residual) <= ledger_tolerance, &
                     'implicit-mean escape charge ledger does not close')

    ! A photoelectron emitted from an illuminated underside exits through the
    ! ordinary z-low boundary without ever entering the z-high mean sheath.
    ! Reconciliation must preserve that non-deferred escape.
    mean_cfg%particle_species(3)%inject_face = 'z_low'
    mean_cfg%particle_species(3)%pos_low = [0.25_dp*box_width, 0.25_dp*box_width, 0.0_dp]
    mean_cfg%particle_species(3)%pos_high = [0.75_dp*box_width, 0.75_dp*box_width, 0.0_dp]
    mean_cfg%particle_species(3)%ray_direction = [0.0_dp, 0.0_dp, 1.0_dp]
    mean_mesh%q_elem = 0.0_dp
    mean_stats = sim_stats()
    mean_ledger = charge_ledger_type()
    mean_diagnostics = electrostatic_diagnostics_type()
    mean_injection_state%macro_residual = 0.0_dp
    call seed_particles_from_config(mean_cfg)
    call run_absorption_insulator( &
      mean_mesh, mean_cfg, mean_stats, inject_state=mean_injection_state, charge_ledger=mean_ledger, &
      electrostatic_diagnostics=mean_diagnostics &
      )

    call assert_equal_i64(mean_ledger%escaped_count(3), mean_ledger%emitted_count(3), &
                          'ordinary z-low photoelectron escapes were not preserved')
    call assert_equal_i64(mean_ledger%absorbed_count(3), 0_i64, &
                          'ordinary z-low photoelectrons were reclassified as absorbed')
    call assert_true(mean_ledger%escaped_to_infinity(3) < 0.0_dp, &
                     'ordinary z-low photoelectron escape charge was erased')
    call assert_close_dp(mean_ledger%interface_outward_gross(3), 0.0_dp, 0.0_dp, &
                         'ordinary z-low escape was incorrectly assigned to the z-high interface')
    call assert_close_dp( &
      mean_diagnostics%implicit_mean_last_returning_pe_column_charge_per_area, &
      0.0_dp, 0.0_dp, &
      'ordinary z-low escape must not contribute to the z-high returning shadow column' &
      )
    ledger_residual = mean_ledger%residual()
    charge_scale = max(abs(mean_ledger%surface_charge_before), abs(mean_ledger%surface_charge_after), tiny(1.0_dp))
    ledger_tolerance = max( &
                       1.0e-12_dp*charge_scale, &
                       16.0_dp*eps0*box_width**2/debye_length*1.0e-12_dp &
                       )
    call assert_true(abs(ledger_residual) <= ledger_tolerance, &
                     'ordinary z-low photoelectron escape charge ledger does not close')
  end subroutine test_implicit_mean_photoelectron_interface_transfer

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
      failure_cfg%sim, failure_cfg%field, failure_cfg%periodic2, failure_cfg%panel, &
      failure_cfg%outer_plasma, failure_cfg%coupling &
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
      failure_cfg%sim, failure_cfg%field, failure_cfg%periodic2, failure_cfg%panel, &
      failure_cfg%outer_plasma, failure_cfg%coupling &
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
    species_cfg%particle_species(2)%z_high_boundary = 'reflect'

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
    neutral_cfg%particle_species(1)%z_high_boundary = 'reflect'
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

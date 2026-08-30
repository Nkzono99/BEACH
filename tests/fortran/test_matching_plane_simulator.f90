!> Tiny matching-plane caseで固定点replay、accepted history、resumeを検証する。
program test_matching_plane_simulator
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: qe
  use bem_types, only: mesh_type, sim_stats, injection_state, bc_open, bc_periodic
  use bem_mesh, only: init_mesh, prepare_periodic2_collision_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config, &
                            particle_inflow_reservoir
  use bem_matching_plane_response, only: &
    matching_plane_response_csv_header, reset_matching_plane_response_snapshot_cache
  use bem_charge_ledger, only: charge_ledger_type, finite_charge_sum
  use test_support, only: &
    test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32, assert_equal_i64, &
    assert_close_dp, delete_file_if_exists
  implicit none

  character(len=*), parameter :: response_path = 'test_matching_plane_simulator_response.csv'
  character(len=*), parameter :: history_path = 'test_matching_plane_simulator_history.csv'
  character(len=*), parameter :: resume_history_path = 'test_matching_plane_simulator_resume_history.csv'
  character(len=*), parameter :: stride_history_path = 'test_matching_plane_simulator_stride_history.csv'
  type(mesh_type) :: mesh
  type(app_config) :: cfg
  type(sim_stats) :: stats, resumed_stats
  type(injection_state) :: inject_state
  type(charge_ledger_type) :: implicit_ledger, no_photo_ledger
  integer :: history_unit
  real(dp) :: cancellation_area

  call cleanup_files()
  call reset_matching_plane_response_snapshot_cache()
  call configure_fixture(mesh, cfg, inject_state)
  call write_affine_response_table(response_path)
  call seed_particles_from_config(cfg)
  call test_init(10)

  call test_begin('accepted_fixed_point_replays_one_particle_batch')
  open (newunit=history_unit, file=history_path, status='replace', action='write')
  call run_absorption_insulator( &
    mesh, cfg, stats, inject_state=inject_state, matching_plane_history_unit=history_unit &
    )
  close (history_unit)

  call assert_true(stats%matching_plane_state_valid, 'accepted matching-plane state was not marked valid')
  call assert_equal_i32(stats%matching_plane_iterations, 2_i32, 'affine fixed point must converge in two iterations')
  call assert_close_dp(stats%matching_plane_residual, 0.0_dp, 0.0_dp, 'fixed-point residual mismatch')
  call assert_close_dp(stats%matching_plane_phi_v, 1.0_dp, 1.0e-14_dp, 'accepted matching potential mismatch')
  call assert_close_dp( &
    stats%matching_plane_response(2), 1.25_dp, 1.0e-14_dp, &
    'accepted electron inward number flux response mismatch' &
    )
  call assert_close_dp( &
    stats%matching_plane_response(3), 1.25_dp, 1.0e-14_dp, &
    'accepted ion inward number flux response mismatch' &
    )
  call assert_close_dp(stats%matching_plane_feedback(1), 1.0_dp, 1.0e-14_dp, 'accepted PE outflow mismatch')
  call assert_close_dp( &
    stats%matching_plane_photoelectron_return_flux_m2_s + &
    stats%matching_plane_photoelectron_escape_flux_m2_s, &
    stats%matching_plane_feedback(1), 1.0e-14_dp, 'PE outward/return/escape budget mismatch' &
    )
  call assert_close_dp( &
    stats%matching_plane_photoelectron_return_flux_m2_s, 0.0_dp, 0.0_dp, &
    'zero barrier unexpectedly returned the PE' &
    )
  call assert_close_dp( &
    stats%matching_plane_photoelectron_escape_flux_m2_s, 1.0_dp, 1.0e-14_dp, &
    'zero barrier PE escape flux mismatch' &
    )
  call assert_equal_i64(stats%processed_particles, 3_i64, 'replayed fixed-point trials leaked into statistics')
  call assert_close_dp( &
    inject_state%boundary_macro_residual(6, 1), 0.25_dp, 1.0e-14_dp, &
    'electron injection residual advanced more than once during replay' &
    )
  call assert_close_dp( &
    inject_state%boundary_macro_residual(6, 2), 0.25_dp, 1.0e-14_dp, &
    'ion injection residual advanced more than once during replay' &
    )
  call assert_history_rows(history_path, 1_i32, 1_i32)
  call test_end()

  call test_begin('accepted_state_seeds_restart_continuation')
  cfg%sim%batch_count = 2_i32
  open (newunit=history_unit, file=resume_history_path, status='replace', action='write')
  call run_absorption_insulator( &
    mesh, cfg, resumed_stats, initial_stats=stats, inject_state=inject_state, &
    matching_plane_history_unit=history_unit &
    )
  close (history_unit)

  call assert_true(resumed_stats%matching_plane_state_valid, 'resumed matching-plane state was not valid')
  call assert_equal_i32(resumed_stats%batches, 2_i32, 'matching-plane continuation batch mismatch')
  call assert_equal_i32( &
    resumed_stats%matching_plane_iterations, 1_i32, &
    'committed fixed-point feedback did not seed the resumed batch' &
    )
  call assert_close_dp(resumed_stats%matching_plane_residual, 0.0_dp, 0.0_dp, 'resumed residual mismatch')
  call assert_equal_i64(resumed_stats%processed_particles, 6_i64, 'resumed particle statistics mismatch')
  call assert_close_dp( &
    inject_state%boundary_macro_residual(6, 1), 0.5_dp, 1.0e-14_dp, &
    'resumed electron residual mismatch' &
    )
  call assert_close_dp( &
    inject_state%boundary_macro_residual(6, 2), 0.5_dp, 1.0e-14_dp, 'resumed ion residual mismatch' &
    )
  call assert_history_rows(resume_history_path, 1_i32, 2_i32)
  call test_end()

  call test_begin('component_absolute_tolerance_accepts_active_axis_without_changing_feedback')
  call configure_fixture(mesh, cfg, inject_state)
  cfg%surface_current%coupling_atol = [1.1_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  call seed_particles_from_config(cfg)
  call run_absorption_insulator(mesh, cfg, resumed_stats, inject_state=inject_state)
  call assert_equal_i32( &
    resumed_stats%matching_plane_iterations, 1_i32, &
    'active-axis absolute tolerance did not accept the first replay' &
    )
  call assert_true( &
    resumed_stats%matching_plane_residual <= cfg%surface_current%coupling_rtol, &
    'mixed matching residual exceeded coupling_rtol' &
    )
  call assert_close_dp( &
    resumed_stats%matching_plane_feedback(1), 1.0_dp, 1.0e-14_dp, &
    'absolute-tolerance acceptance changed the measured feedback' &
    )

  call configure_fixture(mesh, cfg, inject_state)
  cfg%surface_current%coupling_rtol = tiny(1.0_dp)
  cfg%surface_current%coupling_atol = [0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  call seed_particles_from_config(cfg)
  call run_absorption_insulator(mesh, cfg, resumed_stats, inject_state=inject_state)
  call assert_equal_i32( &
    resumed_stats%matching_plane_iterations, 2_i32, &
    'tiny coupling_rtol caused false convergence through tolerance scaling' &
    )
  call test_end()

  call test_begin('underrelaxation_updates_only_unconverged_iterates')
  call configure_fixture(mesh, cfg, inject_state)
  cfg%surface_current%coupling_relaxation = 0.5_dp
  cfg%surface_current%coupling_rtol = 0.13_dp
  call seed_particles_from_config(cfg)
  call run_absorption_insulator(mesh, cfg, resumed_stats, inject_state=inject_state)
  call assert_equal_i32(resumed_stats%matching_plane_iterations, 3_i32, 'underrelaxed iteration count mismatch')
  call assert_close_dp(resumed_stats%matching_plane_phi_v, 0.75_dp, 1.0e-14_dp, 'underrelaxed trial response mismatch')
  call assert_close_dp(resumed_stats%matching_plane_feedback(1), 1.0_dp, 1.0e-14_dp, 'measured feedback mismatch')
  call assert_close_dp(resumed_stats%matching_plane_residual, 0.125_dp, 1.0e-14_dp, 'underrelaxed residual mismatch')
  call test_end()

  call test_begin('history_stride_matches_charge_history_phase')
  call configure_fixture(mesh, cfg, inject_state)
  cfg%sim%batch_count = 3_i32
  call seed_particles_from_config(cfg)
  open (newunit=history_unit, file=stride_history_path, status='replace', action='write')
  call run_absorption_insulator( &
    mesh, cfg, resumed_stats, inject_state=inject_state, matching_plane_history_unit=history_unit, &
    history_stride=2_i32 &
    )
  close (history_unit)
  call assert_history_rows(stride_history_path, 2_i32, 3_i32)
  call test_end()

  call test_begin('zhao_online_runs_without_response_table')
  call configure_fixture(mesh, cfg, inject_state)
  call configure_online_backend(cfg)
  call seed_particles_from_config(cfg)
  call run_absorption_insulator(mesh, cfg, resumed_stats, inject_state=inject_state)
  call assert_true(resumed_stats%matching_plane_state_valid, 'online Zhao state was not committed')
  call assert_equal_i32(resumed_stats%matching_plane_iterations, 1_i32, 'zero-field online Zhao must converge once')
  call assert_close_dp(resumed_stats%matching_plane_phi_v, 0.0_dp, 0.0_dp, 'online Zhao zero-field gauge mismatch')
  call assert_true( &
    all(resumed_stats%matching_plane_response(2:3) > 0.0_dp), &
    'online Zhao did not provide ambient inward fluxes' &
    )
  call test_end()

  call test_begin('implicit_zero_mode_reuses_committed_nonunit_area_endpoint')
  call configure_fixture(mesh, cfg, inject_state)
  call write_implicit_response_table(response_path)
  call reset_matching_plane_response_snapshot_cache()
  cfg%surface_current%implicit_zero_mode = .true.
  cfg%sim%box_max(1) = 2.0_dp
  cfg%sim%batch_count = 2_i32
  cfg%sim%batch_duration = 6.0_dp
  cfg%sim%max_step = 2_i32
  cfg%particle_species(1)%drift_velocity(3) = -10.0_dp
  cfg%particle_species(2)%drift_velocity(3) = -10.0_dp
  cfg%particle_species(3)%emit_current_density_a_m2 = 2.0_dp*qe
  cfg%particle_species(3)%m_particle = 9.1093837015e-31_dp
  call prepare_periodic2_collision_mesh(mesh, cfg%sim)
  call seed_particles_from_config(cfg)
  call run_absorption_insulator( &
    mesh, cfg, resumed_stats, inject_state=inject_state, charge_ledger=implicit_ledger &
    )
  call assert_close_dp( &
    resumed_stats%matching_plane_displacement_c_m2, 4.5_dp*qe, 1.0e-12_dp*qe, &
    'second six-second implicit displacement did not start from committed Q/A' &
    )
  call assert_close_dp(sum(mesh%q_elem), 9.0_dp*qe, 1.0e-12_dp*qe, 'non-unit-area committed charge mismatch')
  call assert_close_dp( &
    implicit_ledger%fixed_emission_target_charge(3), 48.0_dp*qe, 1.0e-12_dp*qe, &
    'implicit PE surface-emission target mismatch' &
    )
  call assert_close_dp( &
    implicit_ledger%fixed_absorbed_target_charge(3), -36.0_dp*qe, 1.0e-12_dp*qe, &
    'implicit PE total-return target mismatch' &
    )
  call assert_close_dp( &
    sum(implicit_ledger%fixed_absorbed_target_charge) + &
    sum(implicit_ledger%fixed_emission_target_charge), &
    9.0_dp*qe, 1.0e-12_dp*qe, 'implicit applied surface-current budget mismatch' &
    )
  call test_end()

  call test_begin('no_photo_zhao_online_runs_with_two_species')
  call configure_fixture(mesh, cfg, inject_state)
  call configure_online_backend(cfg)
  call configure_no_photo_fixture(cfg, inject_state)
  call seed_particles_from_config(cfg)
  call run_absorption_insulator(mesh, cfg, resumed_stats, inject_state=inject_state)
  call assert_true(resumed_stats%matching_plane_state_valid, 'no-PE online Zhao state was not committed')
  call assert_close_dp( &
    resumed_stats%matching_plane_feedback(1), 0.0_dp, 0.0_dp, 'no-PE online feedback has PE flux' &
    )
  call assert_close_dp( &
    resumed_stats%matching_plane_feedback(2), 0.0_dp, 0.0_dp, 'no-PE online feedback has PE energy' &
    )
  call assert_close_dp( &
    resumed_stats%matching_plane_photoelectron_return_flux_m2_s, 0.0_dp, 0.0_dp, &
    'no-PE online state has PE return flux' &
    )
  call assert_close_dp( &
    resumed_stats%matching_plane_photoelectron_escape_flux_m2_s, 0.0_dp, 0.0_dp, &
    'no-PE online state has PE escape flux' &
    )
  call test_end()

  call test_begin('no_photo_implicit_zero_mode_commits_six_second_endpoint')
  call configure_fixture(mesh, cfg, inject_state)
  call configure_no_photo_fixture(cfg, inject_state)
  call write_no_photo_implicit_response_table(response_path)
  call reset_matching_plane_response_snapshot_cache()
  cfg%surface_current%implicit_zero_mode = .true.
  cfg%sim%batch_duration = 6.0_dp
  cfg%sim%max_step = 2_i32
  cfg%particle_species(1)%drift_velocity(3) = -10.0_dp
  cfg%particle_species(2)%drift_velocity(3) = -10.0_dp
  call seed_particles_from_config(cfg)
  call run_absorption_insulator( &
    mesh, cfg, resumed_stats, inject_state=inject_state, charge_ledger=no_photo_ledger &
    )
  call assert_close_dp( &
    resumed_stats%matching_plane_displacement_c_m2, 3.0_dp*qe, 1.0e-12_dp*qe, &
    'no-PE six-second implicit displacement mismatch' &
    )
  call assert_close_dp(sum(mesh%q_elem), 3.0_dp*qe, 1.0e-12_dp*qe, 'no-PE implicit committed charge mismatch')
  call assert_close_dp( &
    sum(no_photo_ledger%fixed_absorbed_target_charge), 3.0_dp*qe, 1.0e-12_dp*qe, &
    'no-PE implicit absorbed-current target mismatch' &
    )
  call assert_close_dp( &
    sum(no_photo_ledger%fixed_emission_target_charge), 0.0_dp, 0.0_dp, &
    'no-PE implicit state created an emission target' &
    )
  call assert_true(all(resumed_stats%matching_plane_feedback(1:2) == 0.0_dp), 'no-PE implicit PE feedback mismatch')
  call assert_close_dp( &
    resumed_stats%matching_plane_photoelectron_return_flux_m2_s, 0.0_dp, 0.0_dp, &
    'no-PE implicit return flux mismatch' &
    )
  call assert_close_dp( &
    resumed_stats%matching_plane_photoelectron_escape_flux_m2_s, 0.0_dp, 0.0_dp, &
    'no-PE implicit escape flux mismatch' &
    )
  call test_end()

  call test_begin('implicit_zero_mode_accepts_cancelling_mesh_charge')
  call configure_fixture(mesh, cfg, inject_state)
  call configure_no_photo_fixture(cfg, inject_state)
  call write_no_photo_implicit_response_table(response_path)
  call reset_matching_plane_response_snapshot_cache()
  cancellation_area = product(cfg%sim%box_max(1:2) - cfg%sim%box_min(1:2))
  mesh%q_elem = [1.0e-11_dp, -1.0e-11_dp]
  cfg%surface_current%implicit_zero_mode = .true.
  cfg%sim%batch_duration = 6.0_dp
  cfg%sim%max_step = 2_i32
  cfg%particle_species(1)%drift_velocity(3) = -10.0_dp
  cfg%particle_species(2)%drift_velocity(3) = -10.0_dp
  call seed_particles_from_config(cfg)
  call run_absorption_insulator(mesh, cfg, resumed_stats, inject_state=inject_state)
  call assert_close_dp( &
    resumed_stats%matching_plane_displacement_c_m2, &
    finite_charge_sum(mesh%q_elem, 'cancelling implicit committed charge')/cancellation_area, 0.0_dp, &
    'cancelling implicit state did not report the committed mesh charge' &
    )
  call assert_close_dp( &
    finite_charge_sum(mesh%q_elem, 'cancelling implicit endpoint charge')/cancellation_area, &
    3.0_dp*qe/cancellation_area, &
    64.0_dp*epsilon(1.0_dp)*sum(abs(mesh%q_elem))/cancellation_area, &
    'cancelling implicit committed charge left the backward-Euler roundoff bound' &
    )
  call test_end()

  call cleanup_files()
  call reset_matching_plane_response_snapshot_cache()
  call test_summary()

contains

  subroutine configure_fixture(fixture_mesh, fixture_cfg, state)
    type(mesh_type), intent(out) :: fixture_mesh
    type(app_config), intent(out) :: fixture_cfg
    type(injection_state), intent(out) :: state
    real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2)

    v0(:, 1) = [0.0_dp, 0.0_dp, 0.25_dp]
    v1(:, 1) = [1.0_dp, 0.0_dp, 0.25_dp]
    v2(:, 1) = [0.0_dp, 1.0_dp, 0.25_dp]
    v0(:, 2) = [1.0_dp, 0.0_dp, 0.25_dp]
    v1(:, 2) = [1.0_dp, 1.0_dp, 0.25_dp]
    v2(:, 2) = [0.0_dp, 1.0_dp, 0.25_dp]
    call init_mesh(fixture_mesh, v0, v1, v2, q0=[0.0_dp, 0.0_dp])
    fixture_mesh%elem_vacuum_sign = 1_i32
    fixture_mesh%vacuum_normals = fixture_mesh%normals

    call default_app_config(fixture_cfg)
    fixture_cfg%sim%rng_seed = 2468_i32
    fixture_cfg%sim%batch_count = 1_i32
    fixture_cfg%sim%dt = 0.1_dp
    fixture_cfg%sim%batch_duration = 1.0_dp
    fixture_cfg%sim%max_step = 1_i32
    fixture_cfg%sim%q_floor = 1.0e-30_dp
    fixture_cfg%sim%field_solver = 'direct'
    fixture_cfg%sim%field_bc_mode = 'periodic2'
    fixture_cfg%sim%use_box = .true.
    fixture_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    fixture_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    fixture_cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    fixture_cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    fixture_cfg%periodic2%nonzero_mode_backend = 'panel_spectral_reference'
    fixture_cfg%periodic2%zero_mode_policy = 'exclude_k0'
    fixture_cfg%periodic2%lower_boundary_model = 'e_bottom_zero'
    fixture_cfg%periodic2%reference_mode_layers = 1_i32
    fixture_cfg%periodic2%panel_quadrature_order = 4_i32
    fixture_cfg%surface_current%model = 'matching_plane_quasistatic'
    fixture_cfg%surface_current%response_table_path = response_path
    fixture_cfg%surface_current%electron_species = 'electron'
    fixture_cfg%surface_current%ion_species = 'ion'
    fixture_cfg%surface_current%photoelectron_species = 'photoelectron'
    fixture_cfg%surface_current%coupling_rtol = 1.0e-12_dp
    fixture_cfg%surface_current%coupling_max_iterations = 4_i32
    fixture_cfg%surface_current%coupling_relaxation = 1.0_dp
    fixture_cfg%n_particle_species = 3_i32

    call configure_boundary_inflow_species(fixture_cfg, 1_i32, 'electron', -qe)
    call configure_boundary_inflow_species(fixture_cfg, 2_i32, 'ion', qe)
    fixture_cfg%particle_species(3) = species_from_defaults()
    fixture_cfg%particle_species(3)%species_key = 'photoelectron'
    fixture_cfg%particle_species(3)%source_mode = 'photo_raycast'
    fixture_cfg%particle_species(3)%rays_per_batch = 1_i32
    fixture_cfg%particle_species(3)%emit_current_density_a_m2 = qe
    fixture_cfg%particle_species(3)%deposit_opposite_charge_on_emit = .true.
    fixture_cfg%particle_species(3)%surface_charge_closure = 'explicit'
    fixture_cfg%particle_species(3)%q_particle = -qe
    fixture_cfg%particle_species(3)%m_particle = 1.0_dp
    fixture_cfg%particle_species(3)%temperature_k = 0.0_dp
    fixture_cfg%particle_species(3)%normal_drift_speed = 10.0_dp
    fixture_cfg%particle_species(3)%inject_face = 'z_high'
    fixture_cfg%particle_species(3)%pos_low = [0.0_dp, 0.0_dp, 1.0_dp]
    fixture_cfg%particle_species(3)%pos_high = [1.0_dp, 1.0_dp, 1.0_dp]
    fixture_cfg%particle_species(3)%ray_direction = [0.0_dp, 0.0_dp, -1.0_dp]

    allocate (state%macro_residual(3), state%boundary_macro_residual(6, 3))
    state%macro_residual = 0.0_dp
    state%boundary_macro_residual = 0.0_dp
    call prepare_periodic2_collision_mesh(fixture_mesh, fixture_cfg%sim)
  end subroutine configure_fixture

  subroutine configure_boundary_inflow_species(fixture_cfg, species_idx, species_key, particle_charge)
    type(app_config), intent(inout) :: fixture_cfg
    integer(i32), intent(in) :: species_idx
    character(len=*), intent(in) :: species_key
    real(dp), intent(in) :: particle_charge

    fixture_cfg%particle_species(species_idx) = species_from_defaults()
    fixture_cfg%particle_species(species_idx)%species_key = species_key
    fixture_cfg%particle_species(species_idx)%source_mode = 'volume_seed'
    fixture_cfg%particle_species(species_idx)%npcls_per_step = 0_i32
    fixture_cfg%particle_species(species_idx)%boundary_inflow_high(3) = particle_inflow_reservoir
    fixture_cfg%particle_species(species_idx)%surface_charge_closure = 'explicit'
    fixture_cfg%particle_species(species_idx)%q_particle = particle_charge
    fixture_cfg%particle_species(species_idx)%m_particle = 1.0_dp
    fixture_cfg%particle_species(species_idx)%w_particle = 1.0_dp
    fixture_cfg%particle_species(species_idx)%temperature_k = 0.0_dp
    fixture_cfg%particle_species(species_idx)%drift_velocity = [0.0_dp, 0.0_dp, -0.1_dp]
  end subroutine configure_boundary_inflow_species

  subroutine configure_online_backend(fixture_cfg)
    type(app_config), intent(inout) :: fixture_cfg

    fixture_cfg%surface_current%response_backend = 'zhao_online'
    fixture_cfg%surface_current%response_table_path = ''
    fixture_cfg%surface_current%zhao_branch = 'auto'
    fixture_cfg%sim%batch_duration = 1.0e-6_dp
    fixture_cfg%sim%dt = 1.0e-9_dp

    fixture_cfg%particle_species(1)%m_particle = 9.1093837015e-31_dp
    fixture_cfg%particle_species(1)%w_particle = 1.0e6_dp
    fixture_cfg%particle_species(1)%number_density_m3 = 8.7e6_dp
    fixture_cfg%particle_species(1)%has_number_density_m3 = .true.
    fixture_cfg%particle_species(1)%drift_velocity(3) = -4.0529988897111727e5_dp
    fixture_cfg%particle_species(1)%temperature_ev = 12.0_dp
    fixture_cfg%particle_species(1)%has_temperature_ev = .true.

    fixture_cfg%particle_species(2)%m_particle = 1.67262192369e-27_dp
    fixture_cfg%particle_species(2)%w_particle = 1.0e6_dp
    fixture_cfg%particle_species(2)%number_density_m3 = 8.7e6_dp
    fixture_cfg%particle_species(2)%has_number_density_m3 = .true.
    fixture_cfg%particle_species(2)%drift_velocity(3) = -4.0529988897111727e5_dp
    fixture_cfg%particle_species(2)%temperature_ev = 0.1_dp
    fixture_cfg%particle_species(2)%has_temperature_ev = .true.

    fixture_cfg%particle_species(3)%m_particle = fixture_cfg%particle_species(1)%m_particle
    fixture_cfg%particle_species(3)%w_particle = 1.0e6_dp
    fixture_cfg%particle_species(3)%temperature_ev = 3.0_dp
    fixture_cfg%particle_species(3)%has_temperature_ev = .true.
    fixture_cfg%particle_species(3)%rays_per_batch = 0_i32
    fixture_cfg%particle_species(3)%emit_current_density_a_m2 = 0.0_dp
  end subroutine configure_online_backend

  subroutine configure_no_photo_fixture(fixture_cfg, state)
    type(app_config), intent(inout) :: fixture_cfg
    type(injection_state), intent(inout) :: state

    fixture_cfg%surface_current%has_photoelectron_species = .false.
    fixture_cfg%surface_current%photoelectron_species = ''
    fixture_cfg%n_particle_species = 2_i32
    if (allocated(state%macro_residual)) deallocate (state%macro_residual)
    if (allocated(state%boundary_macro_residual)) deallocate (state%boundary_macro_residual)
    allocate (state%macro_residual(2), state%boundary_macro_residual(6, 2))
    state%macro_residual = 0.0_dp
    state%boundary_macro_residual = 0.0_dp
  end subroutine configure_no_photo_fixture

  subroutine write_affine_response_table(path)
    character(len=*), intent(in) :: path
    integer :: unit_id
    real(dp) :: low_row(11), high_row(11)

    low_row = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.25_dp, 1.25_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    high_row = [0.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 1.25_dp, 1.25_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    open (newunit=unit_id, file=path, status='replace', action='write')
    write (unit_id, '(a)') '# matching_plane_z_m=1.0'
    write (unit_id, '(a)') matching_plane_response_csv_header
    write (unit_id, '(11(es24.16,:,","))') low_row
    write (unit_id, '(11(es24.16,:,","))') high_row
    close (unit_id)
  end subroutine write_affine_response_table

  subroutine write_implicit_response_table(path)
    character(len=*), intent(in) :: path
    integer :: unit_id
    real(dp) :: low_row(11), high_row(11)

    low_row = [ &
              0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
              1.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 1.0_dp - log(2.0_dp) &
              ]
    high_row = [ &
               12.0_dp*qe, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
               1.0_dp, 2.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 1.0_dp - log(2.0_dp) &
               ]
    open (newunit=unit_id, file=path, status='replace', action='write')
    write (unit_id, '(a)') '# matching_plane_z_m=1.0'
    write (unit_id, '(a)') matching_plane_response_csv_header
    write (unit_id, '(11(es24.16,:,","))') low_row
    write (unit_id, '(11(es24.16,:,","))') high_row
    close (unit_id)
  end subroutine write_implicit_response_table

  subroutine write_no_photo_implicit_response_table(path)
    character(len=*), intent(in) :: path
    integer :: unit_id
    real(dp) :: low_row(11), high_row(11)

    low_row = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    high_row = low_row
    high_row(1) = 12.0_dp*qe
    open (newunit=unit_id, file=path, status='replace', action='write')
    write (unit_id, '(a)') '# matching_plane_z_m=1.0'
    write (unit_id, '(a)') matching_plane_response_csv_header
    write (unit_id, '(11(es24.16,:,","))') low_row
    write (unit_id, '(11(es24.16,:,","))') high_row
    close (unit_id)
  end subroutine write_no_photo_implicit_response_table

  subroutine assert_history_rows(path, expected_rows, expected_batch)
    character(len=*), intent(in) :: path
    integer(i32), intent(in) :: expected_rows, expected_batch
    integer :: unit_id, ios
    integer(i32) :: row_count, batch_idx, last_batch, iterations
    real(dp) :: values(15), last_values(15)

    row_count = 0_i32
    last_batch = -1_i32
    last_values = 0.0_dp
    open (newunit=unit_id, file=path, status='old', action='read')
    do
      read (unit_id, *, iostat=ios) batch_idx, values(1:14), iterations, values(15)
      if (ios /= 0) exit
      row_count = row_count + 1_i32
      last_batch = batch_idx
      last_values = values
    end do
    close (unit_id)
    call assert_equal_i32(row_count, expected_rows, 'matching-plane accepted history row count mismatch')
    call assert_equal_i32(last_batch, expected_batch, 'matching-plane accepted history batch mismatch')
    call assert_close_dp(last_values(3), 1.0_dp, 1.0e-14_dp, 'matching-plane history phi_H mismatch')
    call assert_close_dp(last_values(4), 1.25_dp, 1.0e-14_dp, 'matching-plane history electron inward flux mismatch')
    call assert_close_dp(last_values(5), 1.25_dp, 1.0e-14_dp, 'matching-plane history ion inward flux mismatch')
    call assert_close_dp(last_values(6), 0.0_dp, 0.0_dp, 'matching-plane history electron access mismatch')
    call assert_close_dp(last_values(7), 0.0_dp, 0.0_dp, 'matching-plane history ion access mismatch')
    call assert_close_dp(last_values(8), 0.0_dp, 0.0_dp, 'matching-plane history PE barrier mismatch')
  end subroutine assert_history_rows

  subroutine cleanup_files()
    call delete_file_if_exists(response_path)
    call delete_file_if_exists(history_path)
    call delete_file_if_exists(resume_history_path)
    call delete_file_if_exists(stride_history_path)
  end subroutine cleanup_files

end program test_matching_plane_simulator

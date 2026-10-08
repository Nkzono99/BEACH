!> Flat periodic plateでzhao_stationaryの外部障壁ゲージと観測PE流出による弱連成を検証する。
program test_zhao_outflow_refresh
  use bem_kinds, only: dp, i32
  use bem_constants, only: qe
  use bem_types, only: mesh_type, sim_stats, injection_state, bc_open, bc_periodic
  use bem_mesh, only: init_mesh, prepare_periodic2_collision_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config, &
                            particle_inflow_reservoir
  use bem_charge_ledger, only: charge_ledger_type, finite_charge_sum
  use bem_surface_current_model, only: surface_current_model_result_type, evaluate_surface_current_model
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32, &
                          assert_close_dp
  implicit none

  real(dp), parameter :: cell_width = 1.0e-4_dp
  real(dp), parameter :: box_height = 1.0e-3_dp
  type(mesh_type) :: mesh
  type(app_config) :: cfg
  type(sim_stats) :: stats, resumed_stats
  type(injection_state) :: inject_state
  type(charge_ledger_type) :: ledger
  type(surface_current_model_result_type) :: static_root
  real(dp) :: emission_flux, predicted_escape_fraction, observed_escape_fraction

  call test_init(3)

  call test_begin('stationary_barrier_uses_outer_wall_potential_gauge')
  call configure_fixture(mesh, cfg, inject_state, 0_i32)
  call evaluate_surface_current_model(cfg, static_root)
  call seed_particles_from_config(cfg)
  call run_absorption_insulator(mesh, cfg, stats, inject_state=inject_state, charge_ledger=ledger)
  call assert_true(.not. stats%matching_plane_state_valid, 'a fixed Zhao root must not publish refresh state')
  ! A flat plate returns every PE that reaches H below phi0-phi_m. With the H gauge at phi0, the tracked
  ! escape fraction is the half-Maxwellian tail exp(-(phi0-phi_m)/T) of the same root.
  predicted_escape_fraction = static_root%photoelectron_escape_current_density_a_m2/ &
                              static_root%photoelectron_emission_current_density_a_m2
  observed_escape_fraction = ledger%escaped_to_infinity(3)/ledger%emitted_from_surface(3)
  call assert_close_dp( &
    observed_escape_fraction, predicted_escape_fraction, 0.03_dp, &
    'tracked PE escape fraction must follow the outer barrier above phi0' &
    )
  call assert_close_dp( &
    finite_charge_sum(mesh%q_elem, 'fixed-root surface charge'), 0.0_dp, &
    1.0e-6_dp*abs(ledger%emitted_from_surface(3)), 'zero-current targets must keep the cell floating' &
    )
  call test_end()

  call test_begin('flat_plate_outflow_refresh_reproduces_emission_source')
  call configure_fixture(mesh, cfg, inject_state, 1_i32)
  call seed_particles_from_config(cfg)
  call run_absorption_insulator(mesh, cfg, stats, inject_state=inject_state)
  emission_flux = static_root%photoelectron_emission_current_density_a_m2/qe
  call assert_true(stats%matching_plane_state_valid, 'outflow refresh state was not published')
  call assert_equal_i32( &
    stats%matching_plane_iterations, cfg%sim%batch_count + 1_i32, &
    'every one-batch window must refresh the outer root' &
    )
  call assert_close_dp( &
    stats%matching_plane_feedback(1)/emission_flux, 1.0_dp, 0.02_dp, &
    'every PE leaving a flat plate must cross H once' &
    )
  call assert_close_dp( &
    stats%matching_plane_feedback(2), 2.2_dp, 0.15_dp, &
    'flat-plate outflow mean normal energy must equal the emission temperature' &
    )
  call assert_close_dp(stats%matching_plane_phi_v, static_root%phi0_v, 0.3_dp, 'refreshed phi0 drifted')
  call assert_close_dp( &
    stats%matching_plane_photoelectron_return_flux_m2_s + stats%matching_plane_photoelectron_escape_flux_m2_s, &
    stats%matching_plane_feedback(1), 1.0e-12_dp*stats%matching_plane_feedback(1), &
    'outer source return/escape partition mismatch' &
    )
  call assert_close_dp( &
    stats%matching_plane_displacement_c_m2, 0.0_dp, 1.0e-6_dp*emission_flux*qe*cfg%sim%batch_duration, &
    'refreshed zero-current targets must keep the cell floating' &
    )
  call test_end()

  call test_begin('outflow_refresh_state_reconstructs_on_resume')
  cfg%sim%batch_count = stats%batches + 1_i32
  call run_absorption_insulator(mesh, cfg, resumed_stats, initial_stats=stats, inject_state=inject_state)
  call assert_true(resumed_stats%matching_plane_state_valid, 'resumed outflow refresh state was not valid')
  call assert_equal_i32( &
    resumed_stats%matching_plane_iterations, stats%matching_plane_iterations + 1_i32, &
    'resume must continue the accepted outer-solve count' &
    )
  call assert_close_dp(resumed_stats%matching_plane_phi_v, stats%matching_plane_phi_v, 0.3_dp, 'resumed phi0 jump')
  call test_end()

  call test_summary()

contains

  subroutine configure_fixture(fixture_mesh, fixture_cfg, state, refresh_batches)
    type(mesh_type), intent(out) :: fixture_mesh
    type(app_config), intent(out) :: fixture_cfg
    type(injection_state), intent(out) :: state
    integer(i32), intent(in) :: refresh_batches
    real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2)
    real(dp), parameter :: plate_z = 2.0e-6_dp
    real(dp), parameter :: inward_speed = 4.0529988897e5_dp

    v0(:, 1) = [0.0_dp, 0.0_dp, plate_z]
    v1(:, 1) = [cell_width, 0.0_dp, plate_z]
    v2(:, 1) = [0.0_dp, cell_width, plate_z]
    v0(:, 2) = [cell_width, 0.0_dp, plate_z]
    v1(:, 2) = [cell_width, cell_width, plate_z]
    v2(:, 2) = [0.0_dp, cell_width, plate_z]
    call init_mesh(fixture_mesh, v0, v1, v2, q0=[0.0_dp, 0.0_dp])
    fixture_mesh%elem_vacuum_sign = 1_i32
    fixture_mesh%vacuum_normals = fixture_mesh%normals

    call default_app_config(fixture_cfg)
    fixture_cfg%sim%rng_seed = 1357_i32
    fixture_cfg%sim%batch_count = 2_i32
    fixture_cfg%sim%dt = 1.0e-11_dp
    fixture_cfg%sim%batch_duration = 1.0e-3_dp
    fixture_cfg%sim%max_step = 20000_i32
    fixture_cfg%sim%q_floor = 1.0e-30_dp
    fixture_cfg%sim%field_solver = 'direct'
    fixture_cfg%sim%field_bc_mode = 'periodic2'
    fixture_cfg%sim%use_box = .true.
    fixture_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    fixture_cfg%sim%box_max = [cell_width, cell_width, box_height]
    fixture_cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    fixture_cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    fixture_cfg%periodic2%nonzero_mode_backend = 'panel_spectral_reference'
    fixture_cfg%periodic2%zero_mode_policy = 'exclude_k0'
    fixture_cfg%periodic2%lower_boundary_model = 'e_bottom_zero'
    fixture_cfg%periodic2%reference_mode_layers = 1_i32
    fixture_cfg%periodic2%panel_quadrature_order = 4_i32
    fixture_cfg%surface_current%model = 'zhao_stationary'
    fixture_cfg%surface_current%zhao_branch = 'auto'
    fixture_cfg%surface_current%electron_species = 'electron'
    fixture_cfg%surface_current%ion_species = 'ion'
    fixture_cfg%surface_current%photoelectron_species = 'photoelectron'
    fixture_cfg%surface_current%solar_elevation_deg = 60.0_dp
    fixture_cfg%surface_current%photoelectron_ref_density_m3 = 64.0e6_dp
    fixture_cfg%surface_current%outflow_refresh_batches = refresh_batches
    fixture_cfg%n_particle_species = 3_i32

    call configure_ambient_species(fixture_cfg, 1_i32, 'electron', -qe, 9.1093837015e-31_dp, 12.0_dp, 0.35_dp)
    call configure_ambient_species(fixture_cfg, 2_i32, 'ion', qe, 1.67262192369e-27_dp, 0.1_dp, 0.175_dp)
    fixture_cfg%particle_species(3) = species_from_defaults()
    fixture_cfg%particle_species(3)%species_key = 'photoelectron'
    fixture_cfg%particle_species(3)%source_mode = 'photo_raycast'
    fixture_cfg%particle_species(3)%rays_per_batch = 2000_i32
    fixture_cfg%particle_species(3)%emit_current_density_a_m2 = 1.0e-4_dp
    fixture_cfg%particle_species(3)%deposit_opposite_charge_on_emit = .true.
    fixture_cfg%particle_species(3)%surface_charge_closure = 'fixed_current'
    fixture_cfg%particle_species(3)%q_particle = -qe
    fixture_cfg%particle_species(3)%m_particle = 9.1093837015e-31_dp
    fixture_cfg%particle_species(3)%temperature_ev = 2.2_dp
    fixture_cfg%particle_species(3)%has_temperature_ev = .true.
    fixture_cfg%particle_species(3)%normal_drift_speed = 0.0_dp
    fixture_cfg%particle_species(3)%inject_face = 'z_high'
    fixture_cfg%particle_species(3)%pos_low = [0.0_dp, 0.0_dp, box_height]
    fixture_cfg%particle_species(3)%pos_high = [cell_width, cell_width, box_height]
    fixture_cfg%particle_species(3)%ray_direction = [0.0_dp, 0.0_dp, -1.0_dp]
    fixture_cfg%particle_species(1:2)%drift_velocity(3) = -inward_speed

    allocate (state%macro_residual(3), state%boundary_macro_residual(6, 3))
    state%macro_residual = 0.0_dp
    state%boundary_macro_residual = 0.0_dp
    call prepare_periodic2_collision_mesh(fixture_mesh, fixture_cfg%sim)
  end subroutine configure_fixture

  subroutine configure_ambient_species(fixture_cfg, species_idx, species_key, charge, mass, temperature_ev, weight)
    type(app_config), intent(inout) :: fixture_cfg
    integer(i32), intent(in) :: species_idx
    character(len=*), intent(in) :: species_key
    real(dp), intent(in) :: charge, mass, temperature_ev, weight

    fixture_cfg%particle_species(species_idx) = species_from_defaults()
    fixture_cfg%particle_species(species_idx)%species_key = species_key
    fixture_cfg%particle_species(species_idx)%source_mode = 'volume_seed'
    fixture_cfg%particle_species(species_idx)%npcls_per_step = 0_i32
    fixture_cfg%particle_species(species_idx)%boundary_inflow_high(3) = particle_inflow_reservoir
    fixture_cfg%particle_species(species_idx)%surface_charge_closure = 'fixed_current'
    fixture_cfg%particle_species(species_idx)%q_particle = charge
    fixture_cfg%particle_species(species_idx)%m_particle = mass
    fixture_cfg%particle_species(species_idx)%w_particle = weight
    fixture_cfg%particle_species(species_idx)%number_density_m3 = 8.7e6_dp
    fixture_cfg%particle_species(species_idx)%temperature_ev = temperature_ev
    fixture_cfg%particle_species(species_idx)%has_temperature_ev = .true.
  end subroutine configure_ambient_species

end program test_zhao_outflow_refresh

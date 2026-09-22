!> MPI/OpenMP の実軌道 replay から H の PE 分布と重み付き流束の保存を検証する。
program test_matching_plane_spectrum_mpi
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: qe
  use bem_types, only: mesh_type, sim_stats, injection_state, bc_open, bc_periodic
  use bem_mesh, only: init_mesh, prepare_periodic2_collision_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, &
                            seed_particles_from_config, particle_inflow_reservoir
  use bem_charge_ledger, only: charge_ledger_type
  use bem_mpi, only: mpi_context, mpi_initialize, mpi_shutdown, mpi_allreduce_min_i32_scalar, &
                     mpi_allreduce_max_i32_scalar, mpi_allreduce_min_real_dp_array, mpi_allreduce_max_real_dp_array
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_equal_i32, &
                          assert_equal_i64, assert_close_dp, assert_allclose_1d
  implicit none

  real(dp), parameter :: source_flux = 1.0e9_dp
  integer(i32), parameter :: ray_counts(2) = [33_i32, 1_i32]
  type(mpi_context) :: mpi
  type(mesh_type) :: mesh
  type(app_config) :: cfg
  type(sim_stats) :: stats
  type(injection_state) :: injection
  type(charge_ledger_type) :: ledger
  real(dp), allocatable :: smallest(:), largest(:)
  real(dp) :: normalization, escaped_flux
  integer(i32) :: bins_min, bins_max, present_everywhere
  integer :: run

  call mpi_initialize(mpi)
  call test_init(1)
  call test_begin('MPI_thread_replay_preserves_global_PE_spectrum_and_flux')
  do run = 1, size(ray_counts)
    ! Odd global work exercises unequal rank counts; one ray leaves non-root
    ! ranks with empty local spectra in a multi-rank execution.
    call configure_fixture(mesh, cfg, injection, ray_counts(run))
    call ledger%init(cfg%n_particle_species)
    call seed_particles_from_config(cfg, mpi=mpi)
    call run_absorption_insulator(mesh, cfg, stats, inject_state=injection, mpi=mpi, charge_ledger=ledger)

    call assert_true(stats%matching_plane_state_valid .and. stats%matching_plane_spectral_closure, &
                     'spectrum closure must commit on every rank')
    call assert_equal_i32(stats%batches, 1_i32, 'replay must commit exactly one batch')
    call assert_true(stats%matching_plane_iterations >= 2_i32, 'measured spectrum must replace the bootstrap')
    call assert_true(stats%matching_plane_residual <= cfg%surface_current%coupling_rtol, &
                     'fixed trajectory spectrum must converge')
    call assert_close_dp(stats%matching_plane_phi_v, 0.0_dp, 0.0_dp, 'uncharged snapshot has a zero-field response')
    call assert_equal_i64(ledger%emitted_count(3), int(ray_counts(run), i64), &
                          'global PE rays were duplicated across ranks or replay iterations')
    call assert_equal_i64(ledger%escaped_count(3), int(ray_counts(run), i64), &
                          'zero-barrier PE rays did not all cross H once')
    call assert_equal_i64(ledger%discarded_unresolved_count(3), 0_i64, 'PE transport did not finish')
    normalization = product(cfg%sim%box_max(1:2) - cfg%sim%box_min(1:2))*cfg%sim%batch_duration
    escaped_flux = -ledger%escaped_to_infinity(3)/(qe*normalization)
    call assert_close_dp(-ledger%emitted_from_surface(3)/(qe*normalization), source_flux, &
                         1.0e-5_dp, 'global ray weights reproduce the configured emission current')
    call assert_close_dp(escaped_flux, source_flux, 1.0e-5_dp, 'global escape normalization')

    present_everywhere = merge(1_i32, 0_i32, allocated(stats%matching_plane_pe_observed%flux) .and. &
                               allocated(stats%matching_plane_pe_input%flux))
    call mpi_allreduce_min_i32_scalar(mpi, present_everywhere)
    call assert_equal_i32(present_everywhere, 1_i32, 'histograms must be present on every rank')
    if (present_everywhere == 0_i32) cycle
    call assert_close_dp(stats%matching_plane_pe_observed%total_flux(), source_flux, 1.0e-5_dp, &
                         'thread/rank sums conserve the weighted H crossing flux')
    call assert_close_dp(stats%matching_plane_pe_observed%total_flux(), stats%matching_plane_feedback(1), &
                         1.0e-5_dp, 'histogram agrees with independently accumulated crossing moments')
    call assert_close_dp(stats%matching_plane_photoelectron_return_flux_m2_s + &
                         stats%matching_plane_photoelectron_escape_flux_m2_s, stats%matching_plane_feedback(1), &
                         1.0e-5_dp, 'outward flux equals return plus escape')
    call assert_close_dp(stats%matching_plane_photoelectron_escape_flux_m2_s, escaped_flux, &
                         1.0e-5_dp, 'boundary event escape equals charge-ledger escape')
    call assert_close_dp(stats%matching_plane_model_escape_flux, stats%matching_plane_pe_input%tail_flux(0.0_dp), &
                         1.0e-5_dp, 'provider escape uses the response input distribution')
    call assert_close_dp(stats%matching_plane_pe_input%mean_energy(), stats%matching_plane_response_input(3), &
                         1.0e-12_dp, 'input moment was derived from the same relaxed distribution')
    bins_min = size(stats%matching_plane_pe_observed%flux)
    bins_max = bins_min
    call mpi_allreduce_min_i32_scalar(mpi, bins_min)
    call mpi_allreduce_max_i32_scalar(mpi, bins_max)
    call assert_equal_i32(bins_min, bins_max, 'zero-padding must give every rank identical histogram support')
    if (bins_min /= bins_max) cycle
    smallest = [stats%matching_plane_pe_observed%flux, stats%matching_plane_feedback, &
                stats%matching_plane_response, stats%matching_plane_response_input, stats%matching_plane_model_escape_flux]
    largest = smallest
    call mpi_allreduce_min_real_dp_array(mpi, smallest)
    call mpi_allreduce_max_real_dp_array(mpi, largest)
    call assert_allclose_1d(smallest, largest, 0.0_dp, 'accepted PE spectrum and response differ across ranks')
  end do
  call test_end()
  call test_summary()
  call mpi_shutdown(mpi)

contains

  subroutine configure_fixture(fixture_mesh, fixture_cfg, state, rays)
    type(mesh_type), intent(out) :: fixture_mesh
    type(app_config), intent(out) :: fixture_cfg
    type(injection_state), intent(out) :: state
    integer(i32), intent(in) :: rays
    real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2)
    integer :: species

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
    fixture_cfg%write_output = .false.
    fixture_cfg%sim%rng_seed = 2468_i32
    fixture_cfg%sim%batch_count = 1_i32
    fixture_cfg%sim%dt = 1.0e-7_dp
    fixture_cfg%sim%batch_duration = 1.0e-6_dp
    fixture_cfg%sim%max_step = 2048_i32
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
    fixture_cfg%surface_current%response_backend = 'zhao_online'
    fixture_cfg%surface_current%response_table_path = ''
    fixture_cfg%surface_current%zhao_branch = 'auto'
    fixture_cfg%surface_current%photoelectron_closure = 'energy_spectrum'
    fixture_cfg%surface_current%electron_species = 'electron'
    fixture_cfg%surface_current%ion_species = 'ion'
    fixture_cfg%surface_current%photoelectron_species = 'photoelectron'
    fixture_cfg%surface_current%coupling_rtol = 1.0e-12_dp
    fixture_cfg%surface_current%coupling_max_iterations = 4_i32
    fixture_cfg%surface_current%coupling_relaxation = 1.0_dp
    fixture_cfg%n_particle_species = 3_i32

    do species = 1, 2
      fixture_cfg%particle_species(species) = species_from_defaults()
      fixture_cfg%particle_species(species)%source_mode = 'volume_seed'
      fixture_cfg%particle_species(species)%npcls_per_step = 0_i32
      fixture_cfg%particle_species(species)%boundary_inflow_high(3) = particle_inflow_reservoir
      fixture_cfg%particle_species(species)%surface_charge_closure = 'explicit'
      fixture_cfg%particle_species(species)%w_particle = 1.0e6_dp
      fixture_cfg%particle_species(species)%number_density_m3 = 8.7e6_dp
      fixture_cfg%particle_species(species)%temperature_k = 0.0_dp
      fixture_cfg%particle_species(species)%has_temperature_ev = .true.
      fixture_cfg%particle_species(species)%drift_velocity = [0.0_dp, 0.0_dp, -4.0529988897111727e5_dp]
    end do
    fixture_cfg%particle_species(1)%species_key = 'electron'
    fixture_cfg%particle_species(1)%q_particle = -qe
    fixture_cfg%particle_species(1)%m_particle = 9.1093837015e-31_dp
    fixture_cfg%particle_species(1)%temperature_ev = 12.0_dp
    fixture_cfg%particle_species(2)%species_key = 'ion'
    fixture_cfg%particle_species(2)%q_particle = qe
    fixture_cfg%particle_species(2)%m_particle = 1.67262192369e-27_dp
    fixture_cfg%particle_species(2)%temperature_ev = 0.1_dp
    fixture_cfg%particle_species(3) = species_from_defaults()
    fixture_cfg%particle_species(3)%species_key = 'photoelectron'
    fixture_cfg%particle_species(3)%source_mode = 'photo_raycast'
    fixture_cfg%particle_species(3)%surface_charge_closure = 'explicit'
    fixture_cfg%particle_species(3)%q_particle = -qe
    fixture_cfg%particle_species(3)%m_particle = fixture_cfg%particle_species(1)%m_particle
    fixture_cfg%particle_species(3)%w_particle = 1.0e6_dp
    fixture_cfg%particle_species(3)%temperature_ev = 3.0_dp
    fixture_cfg%particle_species(3)%has_temperature_ev = .true.
    fixture_cfg%particle_species(3)%normal_drift_speed = 10.0_dp
    fixture_cfg%particle_species(3)%rays_per_batch = rays
    fixture_cfg%particle_species(3)%emit_current_density_a_m2 = qe*source_flux
    fixture_cfg%particle_species(3)%deposit_opposite_charge_on_emit = .true.
    fixture_cfg%particle_species(3)%inject_face = 'z_high'
    fixture_cfg%particle_species(3)%pos_low = [0.0_dp, 0.0_dp, 1.0_dp]
    fixture_cfg%particle_species(3)%pos_high = [1.0_dp, 1.0_dp, 1.0_dp]
    fixture_cfg%particle_species(3)%ray_direction = [0.0_dp, 0.0_dp, -1.0_dp]
    allocate (state%macro_residual(3), state%boundary_macro_residual(6, 3))
    state%macro_residual = 0.0_dp
    state%boundary_macro_residual = 0.0_dp
    call prepare_periodic2_collision_mesh(fixture_mesh, fixture_cfg%sim)
  end subroutine configure_fixture

end program test_matching_plane_spectrum_mpi

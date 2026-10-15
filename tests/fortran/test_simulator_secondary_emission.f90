!> 同一凍結場で SEE を追跡し、収支・世代上限・OpenMP 再現性を検証する。
program test_simulator_secondary_emission
!$ use omp_lib
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: qe
  use bem_types, only: mesh_type, sim_stats
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config
  use bem_mesh, only: init_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_charge_ledger, only: charge_ledger_type
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_true, assert_close_dp
  implicit none
  type(app_config) :: cfg
  type(mesh_type) :: mesh
  type(sim_stats) :: stats, replay_stats
  type(charge_ledger_type) :: ledger, replay_ledger
  real(dp) :: reference_charge(2)
  integer :: path

  call test_init(3)
  call test_begin('constant_two_beam_charge_ledger_both_absorption_paths')
  do path = 1, 2
    call fixture(cfg, mesh)
    call ledger%init(2_i32)
    if (path == 2) cfg%sim%box_min(3) = -1.0e-8_dp
    call seed_particles_from_config(cfg)
    call run_absorption_insulator(mesh, cfg, stats, charge_ledger=ledger)
    call assert_true(ledger%absorbed_count(1) == 32_i64, 'every primary must hit the plane')
    call assert_true(ledger%emitted_count(2) == 64_i64, 'one constant-yield-2 hit must create two SEE macros')
    call assert_close_dp(ledger%emitted_from_surface(2), 2.0_dp*ledger%absorbed_on_surface(1), &
                         1.0e-30_dp, 'secondary weight must equal incident weight')
    call assert_close_dp(sum(mesh%q_elem), sum(ledger%absorbed_on_surface) - sum(ledger%emitted_from_surface), &
                         1.0e-29_dp, 'absorption, reaction emission, and reabsorption close surface charge')
    call assert_true(ledger%injected_count(2) == 0_i64, 'SEE must be emitted, not a remote source')
    call assert_true(stats%processed_particles == 96_i64, 'SEE must be tracked in the same batch')
    call assert_true(stats%secondary_emission_count == 64_i64, 'summary emitted count')
    call assert_close_dp(stats%secondary_emitted_charge_c + sum(mesh%q_elem) - sum(ledger%absorbed_on_surface), &
                         0.0_dp, 1.0e-29_dp, 'reaction charge plus emitted particle charge must vanish')
  end do
  call test_end()

  call test_begin('counter_reproducibility_across_thread_counts_and_schedules')
  call fixture(cfg, mesh)
  call ledger%init(2_i32)
  cfg%particle_species(1)%secondary_emission_const_yield = 1.37_dp
!$ call omp_set_dynamic(.false.)
!$ call omp_set_num_threads(1)
!$ call omp_set_schedule(omp_sched_static, 1)
  call seed_particles_from_config(cfg)
  call run_absorption_insulator(mesh, cfg, stats, charge_ledger=ledger)
  reference_charge = mesh%q_elem
  call fixture(cfg, mesh)
  cfg%particle_species(1)%secondary_emission_const_yield = 1.37_dp
!$ call omp_set_num_threads(4)
!$ call omp_set_schedule(omp_sched_dynamic, 1)
  call seed_particles_from_config(cfg)
  call run_absorption_insulator(mesh, cfg, replay_stats, charge_ledger=replay_ledger)
  call assert_true(stats%secondary_emission_count == replay_stats%secondary_emission_count, 'rounding replay')
  call assert_true(all(ledger%absorbed_count == replay_ledger%absorbed_count), 'reabsorption replay')
  call assert_true(all(ledger%escaped_count == replay_ledger%escaped_count), 'escape replay')
  call assert_true(all(abs(mesh%q_elem - reference_charge) < 1.0e-29_dp), 'surface charge replay')
  call test_end()

  call test_begin('generation_cap_does_not_deposit_unemitted_reaction_charge')
  call fixture(cfg, mesh)
  call ledger%init(2_i32)
  cfg%secondary_emission_max_generations = 1_i32
  cfg%sim%e0 = [0.0_dp, 0.0_dp, 1.0e5_dp]
  cfg%particle_species(2)%secondary_emission_enabled = .true.
  cfg%particle_species(2)%secondary_emission_model = 'constant'
  cfg%particle_species(2)%secondary_emission_const_yield = 1.0_dp
  cfg%particle_species(2)%secondary_emission_target_idx = 2_i32
  call seed_particles_from_config(cfg)
  call run_absorption_insulator(mesh, cfg, stats, charge_ledger=ledger)
  call assert_true(stats%secondary_generation_capped_count > 0_i64, 'returned secondaries must hit the generation cap')
  call assert_true(ledger%emitted_count(2) == 64_i64, 'capped tertiary electrons must not be emitted')
  call assert_close_dp(sum(mesh%q_elem), sum(ledger%absorbed_on_surface) - sum(ledger%emitted_from_surface), &
                       1.0e-29_dp, 'generation-capped reaction charge conservation')
  call test_end()
  call test_summary()
contains
  subroutine fixture(app, plate)
    type(app_config), intent(out) :: app
    type(mesh_type), intent(out) :: plate
    real(dp) :: triangles(3, 3, 2)
    triangles(:, 1, 1) = [-0.5_dp, -0.5_dp, 0.0_dp]
    triangles(:, 2, 1) = [0.5_dp, -0.5_dp, 0.0_dp]
    triangles(:, 3, 1) = [0.5_dp, 0.5_dp, 0.0_dp]
    triangles(:, 1, 2) = [-0.5_dp, -0.5_dp, 0.0_dp]
    triangles(:, 2, 2) = [0.5_dp, 0.5_dp, 0.0_dp]
    triangles(:, 3, 2) = [-0.5_dp, 0.5_dp, 0.0_dp]
    call init_mesh(plate, triangles(:, 1, :), triangles(:, 2, :), triangles(:, 3, :))
    plate%elem_vacuum_sign = 1_i32
    plate%vacuum_normals = plate%normals
    call default_app_config(app)
    app%write_output = .false.
    app%sim%field_solver = 'direct'
    app%sim%use_box = .true.
    app%sim%box_min = [-0.5_dp, -0.5_dp, -0.01_dp]
    app%sim%box_max = [0.5_dp, 0.5_dp, 1.0e-4_dp]
    app%sim%dt = 1.0e-11_dp
    app%sim%max_step = 2000_i32
    app%n_particle_species = 2_i32
    app%particle_species(1) = species_from_defaults()
    app%particle_species(1)%species_key = 'incident'
    app%particle_species(1)%npcls_per_step = 32_i32
    app%particle_species(1)%temperature_k = 0.0_dp
    app%particle_species(1)%pos_low = [-0.05_dp, -0.05_dp, 5.0e-7_dp]
    app%particle_species(1)%pos_high = [0.05_dp, 0.05_dp, 5.0e-7_dp]
    app%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -1.0e6_dp]
    app%particle_species(1)%secondary_emission_enabled = .true.
    app%particle_species(1)%secondary_emission_model = 'constant'
    app%particle_species(1)%secondary_emission_const_yield = 2.0_dp
    app%particle_species(1)%secondary_emission_target_idx = 2_i32
    app%particle_species(2) = species_from_defaults()
    app%particle_species(2)%species_key = 'secondary'
    app%particle_species(2)%source_mode = 'secondary_emission'
    app%particle_species(2)%temperature_ev = 2.0_dp
    app%particle_species(2)%has_temperature_ev = .true.
  end subroutine fixture
end program test_simulator_secondary_emission

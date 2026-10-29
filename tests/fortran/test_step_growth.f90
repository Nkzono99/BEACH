!> 弱い場の領域で step を伸ばす `particles.tracking.step_growth` と、強い場で step を半分にする適応 sub-step を検証するテスト。
program test_step_growth
  use bem_kinds, only: dp, i32, i64
  use bem_mesh, only: init_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config, &
                            load_app_config
  use bem_types, only: mesh_type, sim_stats
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_equal_i64, assert_close_dp, delete_file_if_exists
  implicit none

  character(len=*), parameter :: config_path = 'test_step_growth_tmp.toml'
  real(dp), parameter :: detector_cell = 4.0e-5_dp
  type(app_config) :: cfg, cfg_explicit, cfg_off

  call test_init(5)

  call test_begin('step_growth_sets_a_default_threshold')
  call write_grouped_config(config_path, 'step_growth = true')
  call load_app_config(config_path, cfg)
  call assert_true(cfg%sim%step_growth, 'step_growth should be enabled')
  call assert_equal_i32(cfg%sim%step_growth_max, 5_i32, 'default step_growth_max mismatch')
  call assert_close_dp(cfg%sim%substep_field_change_max, 0.05_dp, 0.0_dp, 'step_growth should set the default threshold')
  call assert_close_dp(cfg%sim%step_growth_energy_fraction, 0.01_dp, 0.0_dp, 'default energy fraction mismatch')
  call write_grouped_config(config_path, 'step_growth = true'//new_line('a')//'step_growth_max = 3'// &
                            new_line('a')//'substep_field_change_max = 0.2')
  call load_app_config(config_path, cfg_explicit)
  call assert_equal_i32(cfg_explicit%sim%step_growth_max, 3_i32, 'explicit step_growth_max mismatch')
  call assert_close_dp(cfg_explicit%sim%substep_field_change_max, 0.2_dp, 0.0_dp, 'explicit threshold must be kept')
  call write_grouped_config(config_path, 'substep_field_change_max = 0.1')
  call load_app_config(config_path, cfg_off)
  call assert_true(.not. cfg_off%sim%step_growth, 'step_growth stays off by default')
  call delete_file_if_exists(config_path)
  call test_end()

  call test_begin('uniform_field_particles_take_longer_steps')
  call run_growth_simulation_test(0.0_dp, .true.)
  call test_end()

  call test_begin('charged_surface_limits_growth_but_keeps_absorption')
  ! The field of the element changes by 5% per dt step at the start height, so growth is refused there.
  call run_growth_simulation_test(-1.0e-12_dp, .false.)
  call test_end()

  call test_begin('adaptive_substeps_match_fine_step_impact')
  call test_adaptive_substeps_match_fine_steps()
  call test_end()

  call test_begin('adaptive_substeps_refine_near_charged_panel')
  call test_adaptive_substeps()
  call test_end()

  call test_summary()

contains

  !> An electron thrown at a charged triangle in SI units: the step is longer than the panel's field scale.
  !! Off, the run must not touch the sub-step counters; on, steps near the panel are halved and every halved step
  !! holds at least two sub-steps while the halving cap is not reached.
  subroutine test_adaptive_substeps()
    type(mesh_type) :: mesh_off, mesh_on
    type(app_config) :: cfg_off, cfg_on
    type(sim_stats) :: stats_off, stats_on

    call configure_adaptive_substep_fixture(mesh_off, cfg_off)
    call run_absorption_insulator(mesh_off, cfg_off, stats_off)
    call assert_equal_i64(stats_off%processed_particles, 1_i64, 'adaptive fixture processes one electron')
    call assert_equal_i64(stats_off%adaptive_substep_halved_steps, 0_i64, 'disabled sub-steps must not halve')
    call assert_equal_i64(stats_off%adaptive_substeps, 0_i64, 'disabled sub-steps must not count sub-steps')

    call configure_adaptive_substep_fixture(mesh_on, cfg_on)
    cfg_on%sim%substep_field_change_max = 0.2_dp
    cfg_on%sim%substep_max_halvings = 10_i32
    call run_absorption_insulator(mesh_on, cfg_on, stats_on)
    call assert_equal_i64(stats_on%processed_particles, 1_i64, 'adaptive run processes one electron')
    call assert_true(stats_on%adaptive_substep_halved_steps > 0_i64, 'steps near the charged panel must be halved')
    call assert_true(stats_on%adaptive_substeps >= 2_i64*stats_on%adaptive_substep_halved_steps, &
                     'every halved step holds at least two sub-steps')
    call assert_equal_i64(stats_on%adaptive_substep_depth_limited_steps, 0_i64, &
                          'ten halvings must resolve the panel field without hitting the cap')
    call assert_equal_i64(stats_on%survived_max_step, 0_i64, 'the electron must finish its track within max_step')
    write (*, '(a,i0,a,i0)') 'adaptive sub-steps: halved steps=', stats_on%adaptive_substep_halved_steps, &
      ' sub-steps=', stats_on%adaptive_substeps
  end subroutine test_adaptive_substeps

  !> The adaptive run must land where a run with 1024 times smaller fixed steps lands. The unsplit 100 ps step
  !! jumps past the panel's field and lands three detector cells away, so a run that counts halvings but
  !! advances with the undivided step fails. Element centroids within one cell diagonal count as the same impact.
  subroutine test_adaptive_substeps_match_fine_steps()
    type(mesh_type) :: mesh_off, mesh_on, mesh_fine
    type(app_config) :: cfg_off, cfg_on, cfg_fine
    type(sim_stats) :: stats_off, stats_on, stats_fine
    real(dp) :: impact_off(3), impact_on(3), impact_fine(3)

    call configure_detector_fixture(mesh_off, cfg_off)
    call run_absorption_insulator(mesh_off, cfg_off, stats_off)
    impact_off = detector_impact(mesh_off)

    call configure_detector_fixture(mesh_on, cfg_on)
    cfg_on%sim%substep_field_change_max = 0.2_dp
    cfg_on%sim%substep_max_halvings = 10_i32
    call run_absorption_insulator(mesh_on, cfg_on, stats_on)
    impact_on = detector_impact(mesh_on)

    call configure_detector_fixture(mesh_fine, cfg_fine)
    cfg_fine%sim%dt = cfg_fine%sim%dt/1024.0_dp
    cfg_fine%sim%max_step = 1024_i32*cfg_fine%sim%max_step
    call run_absorption_insulator(mesh_fine, cfg_fine, stats_fine)
    impact_fine = detector_impact(mesh_fine)

    call assert_true(stats_on%adaptive_substep_halved_steps > 0_i64, 'the fixture must exercise sub-step halving')
    call assert_equal_i64(stats_fine%absorbed, 1_i64, 'the fine-step electron must reach the detector')
    call assert_equal_i64(stats_on%absorbed, 1_i64, 'the adaptive electron must reach the detector')
    call assert_true(norm2(impact_on - impact_fine) <= sqrt(2.0_dp)*detector_cell, &
                     'adaptive sub-steps must reproduce the fine-step impact on the detector')
    call assert_true(norm2(impact_off - impact_fine) > 2.0_dp*detector_cell, &
                     'the unsplit coarse step must miss the fine-step impact (fixture sensitivity)')
  end subroutine test_adaptive_substeps_match_fine_steps

  !> Centroid of the detector element that took the electron, or a point 1 m away when none did.
  function detector_impact(fixture_mesh) result(impact)
    type(mesh_type), intent(in) :: fixture_mesh
    real(dp) :: impact(3)
    integer :: elem

    impact = 1.0_dp
    if (.not. any(fixture_mesh%q_elem(2:) /= 0.0_dp)) return
    elem = maxloc(abs(fixture_mesh%q_elem(2:)), dim=1) + 1
    impact = fixture_mesh%centers(:, elem)
  end function detector_impact

  !> The charged triangle of the sub-step fixture plus an uncharged 10 x 10 detector plane at x = 0 behind
  !! the start point. The electron starts 6 um off axis, is turned back by the panel and lands on the detector.
  subroutine configure_detector_fixture(fixture_mesh, fixture_cfg)
    type(mesh_type), intent(out) :: fixture_mesh
    type(app_config), intent(out) :: fixture_cfg
    integer, parameter :: nd = 10
    real(dp) :: v0(3, 1 + 2*nd*nd), v1(3, 1 + 2*nd*nd), v2(3, 1 + 2*nd*nd), q0(1 + 2*nd*nd), d
    integer :: iy, iz, n

    v0(:, 1) = [1.0e-4_dp, -1.0e-5_dp, -1.0e-5_dp]
    v1(:, 1) = [1.0e-4_dp, 1.0e-5_dp, -1.0e-5_dp]
    v2(:, 1) = [1.0e-4_dp, 0.0_dp, 1.0e-5_dp]
    d = detector_cell
    n = 1
    do iy = 0, nd - 1
      do iz = 0, nd - 1
        n = n + 1
        v0(:, n) = [0.0_dp, -2.0e-4_dp + d*iy, -2.0e-4_dp + d*iz]
        v1(:, n) = [0.0_dp, -2.0e-4_dp + d*(iy + 1), -2.0e-4_dp + d*iz]
        v2(:, n) = [0.0_dp, -2.0e-4_dp + d*(iy + 1), -2.0e-4_dp + d*(iz + 1)]
        n = n + 1
        v0(:, n) = [0.0_dp, -2.0e-4_dp + d*iy, -2.0e-4_dp + d*iz]
        v1(:, n) = [0.0_dp, -2.0e-4_dp + d*(iy + 1), -2.0e-4_dp + d*(iz + 1)]
        v2(:, n) = [0.0_dp, -2.0e-4_dp + d*iy, -2.0e-4_dp + d*(iz + 1)]
      end do
    end do
    q0 = 0.0_dp
    q0(1) = -2.0e-14_dp
    call init_mesh(fixture_mesh, v0, v1, v2, q0=q0)
    fixture_mesh%elem_vacuum_sign = 1_i32
    fixture_mesh%vacuum_normals = fixture_mesh%normals
    call default_app_config(fixture_cfg)
    fixture_cfg%sim%rng_seed = 996_i32
    fixture_cfg%sim%batch_count = 1_i32
    fixture_cfg%sim%dt = 1.0e-10_dp
    fixture_cfg%sim%max_step = 400_i32
    fixture_cfg%sim%field_solver = 'direct'
    fixture_cfg%sim%use_box = .true.
    fixture_cfg%sim%box_min = [-2.0e-4_dp, -2.0e-4_dp, -2.0e-4_dp]
    fixture_cfg%sim%box_max = [2.0e-4_dp, 2.0e-4_dp, 2.0e-4_dp]
    fixture_cfg%n_particle_species = 1_i32
    fixture_cfg%particle_species(1) = species_from_defaults()
    fixture_cfg%particle_species(1)%source_mode = 'volume_seed'
    fixture_cfg%particle_species(1)%npcls_per_step = 1_i32
    fixture_cfg%particle_species(1)%q_particle = -1.602176634e-19_dp
    fixture_cfg%particle_species(1)%m_particle = 9.1093837015e-31_dp
    fixture_cfg%particle_species(1)%w_particle = 1.0_dp
    fixture_cfg%particle_species(1)%pos_low = [2.0e-5_dp, 6.0e-6_dp, 3.0e-6_dp]
    fixture_cfg%particle_species(1)%pos_high = fixture_cfg%particle_species(1)%pos_low
    fixture_cfg%particle_species(1)%drift_velocity = [1.45e6_dp, 0.0_dp, 0.0_dp]
    fixture_cfg%particle_species(1)%temperature_k = 0.0_dp
    call seed_particles_from_config(fixture_cfg)
  end subroutine configure_detector_fixture

  !> One electron starts 80 um from a 20 um triangle carrying -2e-14 C and moves toward it with 6 eV;
  !! dt gives 20 um steps, longer than the panel's field scale.
  subroutine configure_adaptive_substep_fixture(fixture_mesh, fixture_cfg)
    type(mesh_type), intent(out) :: fixture_mesh
    type(app_config), intent(out) :: fixture_cfg
    real(dp) :: tri_v0(3, 1), tri_v1(3, 1), tri_v2(3, 1), q0(1)

    tri_v0(:, 1) = [1.0e-4_dp, -1.0e-5_dp, -1.0e-5_dp]
    tri_v1(:, 1) = [1.0e-4_dp, 1.0e-5_dp, -1.0e-5_dp]
    tri_v2(:, 1) = [1.0e-4_dp, 0.0_dp, 1.0e-5_dp]
    q0 = [-2.0e-14_dp]
    call init_mesh(fixture_mesh, tri_v0, tri_v1, tri_v2, q0=q0)
    fixture_mesh%elem_vacuum_sign = 1_i32
    fixture_mesh%vacuum_normals = fixture_mesh%normals
    call default_app_config(fixture_cfg)
    fixture_cfg%sim%rng_seed = 996_i32
    fixture_cfg%sim%batch_count = 1_i32
    fixture_cfg%sim%dt = 1.4e-11_dp
    fixture_cfg%sim%max_step = 400_i32
    fixture_cfg%sim%use_box = .true.
    fixture_cfg%sim%box_min = [-2.0e-4_dp, -2.0e-4_dp, -2.0e-4_dp]
    fixture_cfg%sim%box_max = [2.0e-4_dp, 2.0e-4_dp, 2.0e-4_dp]
    fixture_cfg%n_particle_species = 1_i32
    fixture_cfg%particle_species(1) = species_from_defaults()
    fixture_cfg%particle_species(1)%source_mode = 'volume_seed'
    fixture_cfg%particle_species(1)%npcls_per_step = 1_i32
    fixture_cfg%particle_species(1)%q_particle = -1.602176634e-19_dp
    fixture_cfg%particle_species(1)%m_particle = 9.1093837015e-31_dp
    fixture_cfg%particle_species(1)%w_particle = 1.0_dp
    fixture_cfg%particle_species(1)%pos_low = [2.0e-5_dp, 1.0e-6_dp, 2.0e-6_dp]
    fixture_cfg%particle_species(1)%pos_high = fixture_cfg%particle_species(1)%pos_low
    fixture_cfg%particle_species(1)%drift_velocity = [1.45e6_dp, 0.0_dp, 0.0_dp]
    fixture_cfg%particle_species(1)%temperature_k = 0.0_dp
    call seed_particles_from_config(fixture_cfg)
  end subroutine configure_adaptive_substep_fixture

  subroutine write_grouped_config(path, tracking_extra)
    character(len=*), intent(in) :: path, tracking_extra
    integer :: unit

    open (newunit=unit, file=path, status='replace', action='write')
    write (unit, '(a)') '[run]'
    write (unit, '(a)') 'batch_count = 2'
    write (unit, '(a)') '[run.batch]'
    write (unit, '(a)') 'duration_s = 2.0e-6'
    write (unit, '(a)') '[particles.tracking]'
    write (unit, '(a)') 'dt_s = 1.0e-9'
    write (unit, '(a)') trim(tracking_extra)
    write (unit, '(a)') '[domain]'
    write (unit, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (unit, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    write (unit, '(a)') '[[particles.species]]'
    write (unit, '(a)') 'species_key = "electron"'
    write (unit, '(a)') '[particles.species.distribution]'
    write (unit, '(a)') 'number_density_m3 = 1.0e6'
    write (unit, '(a)') 'temperature_ev = 1.0'
    write (unit, '(a)') '[particles.species.sampling]'
    write (unit, '(a)') 'target_macro_particles_per_batch = 100'
    write (unit, '(a)') '[particles.species.inflow]'
    write (unit, '(a)') 'z_high = "reservoir"'
    write (unit, '(a)') '[fields.external]'
    write (unit, '(a)') 'electric_v_m = [0.0, 0.0, -5.0]'
    write (unit, '(a)') '[output]'
    write (unit, '(a)') 'enabled = false'
    close (unit)
  end subroutine write_grouped_config

  !> 4 粒子が高さ 0.8 から 0.02 m/s で 1 要素へ降りる。伸長なしでは 40 push、伸長ありでは数 push で吸収される。
  subroutine run_growth_simulation_test(surface_charge, expect_growth)
    real(dp), intent(in) :: surface_charge
    logical, intent(in) :: expect_growth
    type(mesh_type) :: mesh_plain, mesh_grown
    type(app_config) :: cfg_plain, cfg_grown
    type(sim_stats) :: stats_plain, stats_grown

    call setup_fixture(mesh_plain, cfg_plain, surface_charge)
    call seed_particles_from_config(cfg_plain)
    call run_absorption_insulator(mesh_plain, cfg_plain, stats_plain)
    call assert_equal_i64(stats_plain%absorbed, 4_i64, 'plain run should absorb every particle')
    call assert_equal_i64(stats_plain%adaptive_grown_steps, 0_i64, 'plain run must not grow steps')

    call setup_fixture(mesh_grown, cfg_grown, surface_charge)
    cfg_grown%sim%step_growth = .true.
    cfg_grown%sim%step_growth_max = 4_i32
    cfg_grown%sim%substep_field_change_max = 0.05_dp
    call seed_particles_from_config(cfg_grown)
    call run_absorption_insulator(mesh_grown, cfg_grown, stats_grown)
    call assert_equal_i64(stats_grown%absorbed, 4_i64, 'grown run should absorb every particle')
    call assert_equal_i64(stats_grown%survived_max_step, 0_i64, 'grown run should not hit the step cap')
    call assert_close_dp(mesh_grown%q_elem(1), mesh_plain%q_elem(1), 1.0e-12_dp, 'deposited charge must not depend on growth')
    if (expect_growth) then
      ! 1 + 2 + 4 + 8 + 16 covers 31 of the 40 steps, so 4 to 5 pushes per particle are longer than dt.
      call assert_true(stats_grown%adaptive_grown_steps >= 16_i64, 'every particle should take long steps')
      call assert_true(stats_grown%adaptive_grown_steps <= 20_i64, 'uniform field should need few long pushes')
    else
      call assert_equal_i64(stats_grown%adaptive_grown_steps, 0_i64, 'a strongly varying field must refuse growth')
    end if
  end subroutine run_growth_simulation_test

  subroutine setup_fixture(mesh, cfg, surface_charge)
    type(mesh_type), intent(out) :: mesh
    type(app_config), intent(out) :: cfg
    real(dp), intent(in) :: surface_charge
    real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1)

    v0(:, 1) = [-1.0_dp, -1.0_dp, 0.0_dp]
    v1(:, 1) = [1.0_dp, -1.0_dp, 0.0_dp]
    v2(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
    call init_mesh(mesh, v0, v1, v2)
    mesh%elem_vacuum_sign = 1_i32
    mesh%vacuum_normals = mesh%normals
    mesh%q_elem(1) = surface_charge

    call default_app_config(cfg)
    cfg%sim%rng_seed = 2468_i32
    cfg%sim%batch_count = 1_i32
    cfg%sim%dt = 1.0_dp
    cfg%sim%max_step = 400_i32
    cfg%sim%q_floor = 1.0e-30_dp
    cfg%sim%field_solver = 'direct'
    cfg%field%backend = 'direct'
    cfg%panel%kernel_id = 'triangle_p0_exact_direct'
    cfg%sim%use_box = .true.
    cfg%sim%box_min = [-1.0_dp, -1.0_dp, -2.0_dp]
    cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    cfg%n_particle_species = 1_i32
    cfg%particle_species(1) = species_from_defaults()
    cfg%particle_species(1)%source_mode = 'volume_seed'
    cfg%particle_species(1)%npcls_per_step = 4_i32
    cfg%particle_species(1)%q_particle = 1.0_dp
    cfg%particle_species(1)%m_particle = 1.0_dp
    cfg%particle_species(1)%w_particle = 1.0_dp
    cfg%particle_species(1)%pos_low = [0.1_dp, 0.1_dp, 0.8_dp]
    cfg%particle_species(1)%pos_high = [0.1_dp, 0.1_dp, 0.8_dp]
    cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -0.02_dp]
    cfg%particle_species(1)%temperature_k = 0.0_dp
  end subroutine setup_fixture

end program test_step_growth

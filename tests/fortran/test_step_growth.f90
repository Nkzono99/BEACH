!> 弱い場の領域で step を伸ばす `particles.tracking.step_growth` を検証するテスト。
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
  type(app_config) :: cfg, cfg_explicit, cfg_off

  call test_init(3)

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

  call test_summary()

contains

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

!> reservoir_face 注入の設定解釈とマクロ粒子数計算を検証するテスト。
program test_reservoir_injection
  use bem_kinds, only: dp, i32
  use bem_app_config, only: app_config, default_app_config, load_app_config, particles_per_batch_from_config, &
                            species_from_defaults, seed_particles_from_config, init_particle_batch_from_config
  use bem_app_config_runtime, only: particle_source_plan_type, build_particle_source_plan, &
                                    compute_face_average_potential, compute_z_high_box_potential_statistics, &
                                    reservoir_face_velocity_correction
  use bem_injection, only: compute_macro_particles_for_batch, &
                           compute_inflow_flux_from_drifting_maxwellian, compute_face_area_from_bounds
  use bem_types, only: particles_soa, injection_state
  use bem_types, only: mesh_type
  use bem_mesh, only: init_mesh
  use bem_constants, only: eps0
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, delete_file_if_exists
  implicit none

  type(app_config) :: cfg_fixed, cfg_auto
  character(len=*), parameter :: cfg_fixed_path = 'test_reservoir_injection_fixed_tmp.toml'
  character(len=*), parameter :: cfg_auto_path = 'test_reservoir_injection_auto_tmp.toml'
  integer(i32) :: n_macro, i, n1, n2
  integer(i32) :: sum1, sum2
  real(dp) :: residual
  real(dp) :: residual1, residual2, ratio
  real(dp) :: gamma1, area1, expected_w1
  real(dp) :: inward_normal(3)

  call test_init(8)

  call test_begin('particle_source_plan_equivalence')
  call test_particle_source_plan_equivalence()
  call test_end()

  call test_begin('snapshot_infinity_barrier')
  call test_snapshot_infinity_barrier()
  call test_end()

  call test_begin('face_potential_statistics_share_sampling_pass')
  call test_face_potential_statistics_share_sampling_pass()
  call test_end()

  call write_fixed_duration_fixture(cfg_fixed_path)

  call test_begin('fixed_duration_config')
  call default_app_config(cfg_fixed)
  call load_app_config(cfg_fixed_path, cfg_fixed)

  call assert_true(trim(cfg_fixed%particle_species(1)%source_mode) == 'reservoir_face', 'source_mode mismatch')
  call assert_close_dp(cfg_fixed%sim%batch_duration, 1.0d0, 1.0d-12, 'batch_duration mismatch')
  call assert_close_dp(cfg_fixed%particle_species(1)%number_density_cm3, 5.0d0, 1.0d-12, 'density mismatch')
  call assert_close_dp(cfg_fixed%particle_species(1)%temperature_ev, 10.0d0, 1.0d-12, 'temperature_ev mismatch')
  call assert_true(trim(cfg_fixed%particle_species(1)%inject_face) == 'z_low', 'inject_face mismatch')
  call assert_equal_i32(particles_per_batch_from_config(cfg_fixed), 0_i32, 'reservoir static particle count should be zero')
  call test_end()

  call test_begin('deterministic_batches')
  residual = 0.0d0
  call compute_macro_particles_for_batch( &
    1.05d3, 0.0d0, 1.0d0, [0.0d0, 0.0d0, 1.0d0], [0.0d0, 0.0d0, 0.0d0], [1.0d0, 1.0d0, 1.0d0], &
    'z_low', [0.0d0, 0.0d0, 0.0d0], [1.0d0, 1.0d0, 0.0d0], 1.0d0, 1.0d2, residual, n_macro &
    )
  call assert_equal_i32(n_macro, 10_i32, 'first macro particle count mismatch')
  call assert_close_dp(residual, 0.5d0, 1.0d-12, 'first residual mismatch')

  call compute_macro_particles_for_batch( &
    1.05d3, 0.0d0, 1.0d0, [0.0d0, 0.0d0, 1.0d0], [0.0d0, 0.0d0, 0.0d0], [1.0d0, 1.0d0, 1.0d0], &
    'z_low', [0.0d0, 0.0d0, 0.0d0], [1.0d0, 1.0d0, 0.0d0], 1.0d0, 1.0d2, residual, n_macro &
    )
  call assert_equal_i32(n_macro, 11_i32, 'second macro particle count mismatch')
  call assert_close_dp(residual, 0.0d0, 1.0d-12, 'second residual mismatch')
  call test_end()

  call test_begin('vmin_cutoff')
  call compute_macro_particles_for_batch( &
    1.05d3, 0.0d0, 1.0d0, [0.0d0, 0.0d0, 1.0d0], [0.0d0, 0.0d0, 0.0d0], [1.0d0, 1.0d0, 1.0d0], &
    'z_low', [0.0d0, 0.0d0, 0.0d0], [1.0d0, 1.0d0, 0.0d0], 1.0d0, 1.0d2, residual, n_macro, vmin_normal=1.2d0 &
    )
  call assert_equal_i32(n_macro, 0_i32, 'vmin cutoff should block deterministic inflow')
  call assert_close_dp(residual, 0.0d0, 1.0d-12, 'vmin cutoff residual mismatch')
  call test_end()

  call write_auto_duration_fixture(cfg_auto_path)

  call test_begin('auto_duration_config')
  call default_app_config(cfg_auto)
  call load_app_config(cfg_auto_path, cfg_auto)

  call assert_close_dp(cfg_auto%sim%batch_duration, 3.0d0, 1.0d-12, 'batch_duration_step resolution mismatch')
  call assert_true(trim(cfg_auto%particle_species(1)%source_mode) == 'reservoir_face', 'species-1 mode mismatch')
  call assert_true(trim(cfg_auto%particle_species(2)%source_mode) == 'reservoir_face', 'species-2 mode mismatch')

  inward_normal = [0.0d0, 0.0d0, -1.0d0]
  gamma1 = compute_inflow_flux_from_drifting_maxwellian( &
           1000.0d0, 0.0d0, 1.0d0, [0.0d0, 0.0d0, -1.0d0], inward_normal &
           )
  area1 = compute_face_area_from_bounds('z_high', [0.0d0, 0.0d0, 1.0d0], [1.0d0, 1.0d0, 1.0d0])
  expected_w1 = gamma1*area1*cfg_auto%sim%batch_duration/300.0d0
  call assert_close_dp(cfg_auto%particle_species(1)%w_particle, expected_w1, 1.0d-12, 'species-1 auto w mismatch')
  call assert_close_dp(cfg_auto%particle_species(2)%w_particle, expected_w1, 1.0d-12, 'species-2 shared w mismatch')
  call test_end()

  call test_begin('species_ratio')
  residual1 = 0.0d0
  residual2 = 0.0d0
  sum1 = 0_i32
  sum2 = 0_i32
  do i = 1_i32, 100_i32
    call compute_macro_particles_for_batch( &
      1000.0d0, 0.0d0, 1.0d0, [0.0d0, 0.0d0, -1.0d0], [0.0d0, 0.0d0, 0.0d0], [1.0d0, 1.0d0, 1.0d0], &
      'z_high', [0.0d0, 0.0d0, 1.0d0], [1.0d0, 1.0d0, 1.0d0], &
      cfg_auto%sim%batch_duration, cfg_auto%particle_species(1)%w_particle, residual1, n1 &
      )
    call compute_macro_particles_for_batch( &
      250.0d0, 0.0d0, 1.0d0, [0.0d0, 0.0d0, -1.0d0], [0.0d0, 0.0d0, 0.0d0], [1.0d0, 1.0d0, 1.0d0], &
      'z_high', [0.0d0, 0.0d0, 1.0d0], [1.0d0, 1.0d0, 1.0d0], &
      cfg_auto%sim%batch_duration, cfg_auto%particle_species(2)%w_particle, residual2, n2 &
      )
    sum1 = sum1 + n1
    sum2 = sum2 + n2
  end do
  ratio = real(sum2, dp)/real(sum1, dp)
  call assert_close_dp(ratio, 0.25d0, 1.0d-12, 'reservoir species ratio mismatch')
  call test_end()

  call delete_file_if_exists(cfg_fixed_path)
  call delete_file_if_exists(cfg_auto_path)

  call test_summary()

contains

  subroutine test_particle_source_plan_equivalence()
    type(app_config) :: plan_cfg
    type(particle_source_plan_type) :: source_plan
    type(particles_soa) :: particles_without_plan, particles_with_plan
    type(injection_state) :: state_without_plan, state_with_plan
    integer, allocatable :: rng_before(:), rng_after_build(:), rng_without_plan(:), rng_with_plan(:)
    integer :: rng_size

    call default_app_config(plan_cfg)
    plan_cfg%sim%rng_seed = 2468_i32
    plan_cfg%sim%batch_count = 2_i32
    plan_cfg%sim%batch_duration = 1.0_dp
    plan_cfg%sim%has_batch_duration = .true.
    plan_cfg%sim%use_box = .true.
    plan_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    plan_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    plan_cfg%n_particle_species = 3_i32

    plan_cfg%particle_species(1) = species_from_defaults()
    plan_cfg%particle_species(1)%source_mode = 'volume_seed'
    plan_cfg%particle_species(1)%npcls_per_step = 2_i32
    plan_cfg%particle_species(1)%q_particle = -1.0_dp
    plan_cfg%particle_species(1)%m_particle = 2.0_dp
    plan_cfg%particle_species(1)%w_particle = 3.0_dp
    plan_cfg%particle_species(1)%pos_low = [0.1_dp, 0.2_dp, 0.3_dp]
    plan_cfg%particle_species(1)%pos_high = [0.4_dp, 0.5_dp, 0.6_dp]
    plan_cfg%particle_species(1)%temperature_k = 1.0_dp
    plan_cfg%particle_species(1)%has_temperature_k = .true.
    plan_cfg%particle_species(1)%drift_velocity = [0.5_dp, -0.25_dp, 0.125_dp]

    plan_cfg%particle_species(2) = species_from_defaults()
    plan_cfg%particle_species(2)%enabled = .false.

    plan_cfg%particle_species(3) = species_from_defaults()
    plan_cfg%particle_species(3)%source_mode = 'reservoir_face'
    plan_cfg%particle_species(3)%number_density_m3 = 4.0_dp
    plan_cfg%particle_species(3)%has_number_density_m3 = .true.
    plan_cfg%particle_species(3)%temperature_k = 0.0_dp
    plan_cfg%particle_species(3)%has_temperature_k = .true.
    plan_cfg%particle_species(3)%q_particle = 1.0_dp
    plan_cfg%particle_species(3)%m_particle = 1.0_dp
    plan_cfg%particle_species(3)%w_particle = 1.0_dp
    plan_cfg%particle_species(3)%has_w_particle = .true.
    plan_cfg%particle_species(3)%inject_face = 'z_low'
    plan_cfg%particle_species(3)%pos_low = [0.0_dp, 0.0_dp, 0.0_dp]
    plan_cfg%particle_species(3)%pos_high = [1.0_dp, 1.0_dp, 0.0_dp]
    plan_cfg%particle_species(3)%drift_velocity = [0.0_dp, 0.0_dp, 1.0_dp]

    allocate (state_without_plan%macro_residual(3), state_with_plan%macro_residual(3))
    state_without_plan%macro_residual = 0.0_dp
    state_with_plan%macro_residual = 0.0_dp
    call random_seed(size=rng_size)
    allocate (rng_before(rng_size), rng_after_build(rng_size), rng_without_plan(rng_size), rng_with_plan(rng_size))
    call seed_particles_from_config(plan_cfg)
    call random_seed(get=rng_before)
    call build_particle_source_plan(plan_cfg, source_plan)
    call random_seed(get=rng_after_build)
    call assert_true(all(rng_before == rng_after_build), 'building a source plan must not consume random numbers')

    call seed_particles_from_config(plan_cfg)
    call init_particle_batch_from_config(plan_cfg, 1_i32, particles_without_plan, state=state_without_plan)
    call random_seed(get=rng_without_plan)

    call seed_particles_from_config(plan_cfg)
    call init_particle_batch_from_config( &
      plan_cfg, 1_i32, particles_with_plan, state=state_with_plan, source_plan=source_plan &
      )
    call random_seed(get=rng_with_plan)

    call assert_equal_i32(particles_with_plan%n, particles_without_plan%n, 'source plan particle count mismatch')
    call assert_true( &
      all(particles_with_plan%species_id == particles_without_plan%species_id), 'source plan species order mismatch' &
      )
    call assert_true(all(particles_with_plan%x == particles_without_plan%x), 'source plan position sequence mismatch')
    call assert_true(all(particles_with_plan%v == particles_without_plan%v), 'source plan velocity sequence mismatch')
    call assert_true(all(particles_with_plan%q == particles_without_plan%q), 'source plan charge mismatch')
    call assert_true(all(particles_with_plan%m == particles_without_plan%m), 'source plan mass mismatch')
    call assert_true(all(particles_with_plan%w == particles_without_plan%w), 'source plan weight mismatch')
    call assert_true(all(particles_with_plan%alive .eqv. particles_without_plan%alive), 'source plan alive mismatch')
    call assert_true( &
      all(state_with_plan%macro_residual == state_without_plan%macro_residual), 'source plan residual mismatch' &
      )
    call assert_true(all(rng_with_plan == rng_without_plan), 'source plan must preserve the post-batch random state')
  end subroutine test_particle_source_plan_equivalence

  subroutine test_face_potential_statistics_share_sampling_pass()
    type(app_config) :: stats_cfg
    type(mesh_type) :: stats_mesh
    type(electrostatic_snapshot_type) :: snapshot
    real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1)
    real(dp) :: phi_mean, phi_std, phi_min, phi_max, vmin_normal, barrier_normal
    real(dp) :: top_phi_mean, top_phi_std, top_phi_min, top_phi_max

    v0(:, 1) = [0.0_dp, 0.0_dp, 0.25_dp]
    v1(:, 1) = [1.0_dp, 0.0_dp, 0.25_dp]
    v2(:, 1) = [0.0_dp, 1.0_dp, 0.25_dp]
    call init_mesh(stats_mesh, v0, v1, v2, q0=[0.0_dp])
    stats_mesh%elem_vacuum_sign = 1_i32
    stats_mesh%vacuum_normals = stats_mesh%normals

    call default_app_config(stats_cfg)
    stats_cfg%sim%use_box = .true.
    stats_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    stats_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    stats_cfg%sim%e0 = [1.0_dp, 0.0_dp, 0.0_dp]
    stats_cfg%sim%injection_face_phi_grid_n = 2_i32
    stats_cfg%particle_species(1)%source_mode = 'reservoir_face'
    stats_cfg%particle_species(1)%species_key = 'face_variation_warning_test'
    stats_cfg%particle_species(1)%q_particle = 1.0_dp
    stats_cfg%particle_species(1)%m_particle = 1.0_dp
    stats_cfg%particle_species(1)%temperature_k = 0.0_dp
    stats_cfg%particle_species(1)%has_temperature_k = .true.
    stats_cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -1.0_dp]
    stats_cfg%particle_species(1)%inject_face = 'z_high'
    stats_cfg%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, 1.0_dp]
    stats_cfg%particle_species(1)%pos_high = [1.0_dp, 1.0_dp, 1.0_dp]

    call snapshot%init( &
      stats_mesh, stats_cfg%sim, stats_cfg%field, stats_cfg%periodic2, stats_cfg%panel &
      )
    call snapshot%refresh(stats_mesh)
    call compute_face_average_potential( &
      stats_mesh, stats_cfg%sim, stats_cfg%particle_species(1), snapshot, &
      phi_mean, phi_std, phi_min, phi_max &
      )

    call assert_close_dp(phi_mean, -0.5_dp, 1.0e-14_dp, 'face potential mean mismatch')
    call assert_close_dp(phi_std, 0.25_dp, 1.0e-14_dp, 'face potential standard deviation mismatch')
    call assert_close_dp(phi_min, -0.75_dp, 1.0e-14_dp, 'face potential minimum mismatch')
    call assert_close_dp(phi_max, -0.25_dp, 1.0e-14_dp, 'face potential maximum mismatch')

    call compute_z_high_box_potential_statistics( &
      stats_mesh, stats_cfg%sim, snapshot, top_phi_mean, top_phi_std, top_phi_min, top_phi_max &
      )
    call assert_close_dp(top_phi_mean, phi_mean, 1.0e-14_dp, 'z-high box potential mean mismatch')
    call assert_close_dp(top_phi_std, phi_std, 1.0e-14_dp, 'z-high box potential standard deviation mismatch')
    call assert_close_dp(top_phi_min, phi_min, 1.0e-14_dp, 'z-high box potential minimum mismatch')
    call assert_close_dp(top_phi_max, phi_max, 1.0e-14_dp, 'z-high box potential maximum mismatch')

    stats_cfg%sim%reservoir_potential_model = 'infinity_barrier'
    stats_cfg%sim%phi_infty = -1.0_dp
    call reservoir_face_velocity_correction( &
      stats_cfg, stats_cfg%particle_species(1), vmin_normal, barrier_normal, &
      mesh=stats_mesh, snapshot=snapshot, warn_face_variation=.true. &
      )
    call assert_close_dp(barrier_normal, 1.0_dp, 1.0e-14_dp, 'face-average barrier mismatch')
    call assert_close_dp(vmin_normal, 1.0_dp, 1.0e-14_dp, 'face-average cutoff mismatch')
  end subroutine test_face_potential_statistics_share_sampling_pass

  subroutine test_snapshot_infinity_barrier()
    type(app_config) :: barrier_cfg
    type(mesh_type) :: barrier_mesh
    type(particles_soa) :: particles
    type(injection_state) :: state
    type(electrostatic_snapshot_type) :: snapshot
    real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1)

    v0(:, 1) = [0.0_dp, 0.0_dp, 0.25_dp]
    v1(:, 1) = [1.0_dp, 0.0_dp, 0.25_dp]
    v2(:, 1) = [0.0_dp, 1.0_dp, 0.25_dp]
    call init_mesh(barrier_mesh, v0, v1, v2, q0=[0.0_dp])
    barrier_mesh%elem_vacuum_sign = 1_i32
    barrier_mesh%vacuum_normals = barrier_mesh%normals

    call default_app_config(barrier_cfg)
    barrier_cfg%sim%batch_count = 1_i32
    barrier_cfg%sim%batch_duration = 1.0_dp
    barrier_cfg%sim%has_batch_duration = .true.
    barrier_cfg%sim%dt = 0.0_dp
    barrier_cfg%sim%use_box = .true.
    barrier_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    barrier_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    barrier_cfg%sim%e0 = [0.0_dp, 0.0_dp, -1.0_dp]
    barrier_cfg%sim%reservoir_potential_model = 'infinity_barrier'
    barrier_cfg%sim%phi_infty = 0.0_dp
    barrier_cfg%sim%injection_face_phi_grid_n = 2_i32
    barrier_cfg%n_particle_species = 1_i32
    barrier_cfg%particle_species(1) = species_from_defaults()
    barrier_cfg%particle_species(1)%source_mode = 'reservoir_face'
    barrier_cfg%particle_species(1)%number_density_m3 = 10.0_dp
    barrier_cfg%particle_species(1)%has_number_density_m3 = .true.
    barrier_cfg%particle_species(1)%temperature_k = 0.0_dp
    barrier_cfg%particle_species(1)%has_temperature_k = .true.
    barrier_cfg%particle_species(1)%q_particle = 1.0_dp
    barrier_cfg%particle_species(1)%m_particle = 1.0_dp
    barrier_cfg%particle_species(1)%w_particle = 1.0_dp
    barrier_cfg%particle_species(1)%has_w_particle = .true.
    barrier_cfg%particle_species(1)%inject_face = 'z_high'
    barrier_cfg%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, 1.0_dp]
    barrier_cfg%particle_species(1)%pos_high = [1.0_dp, 1.0_dp, 1.0_dp]
    barrier_cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -1.0_dp]
    allocate (state%macro_residual(1))
    state%macro_residual = 0.0_dp

    call snapshot%init( &
      barrier_mesh, barrier_cfg%sim, barrier_cfg%field, barrier_cfg%periodic2, barrier_cfg%panel &
      )
    call snapshot%refresh(barrier_mesh)
    call init_particle_batch_from_config( &
      barrier_cfg, 1_i32, particles, state=state, mesh=barrier_mesh, snapshot=snapshot &
      )
    call assert_equal_i32(particles%n, 0_i32, 'snapshot infinity barrier should block deterministic inflow')
  end subroutine test_snapshot_infinity_barrier

  !> テスト専用の固定 `batch_duration` reservoir_face 設定ファイルを書き出す。
  !! @param[in] path 書き出し先TOMLファイルパス。
  subroutine write_fixed_duration_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open reservoir config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'batch_count = 2'
    write (u, '(a)') 'batch_duration = 1.0'
    write (u, '(a)') 'use_box = true'
    write (u, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (u, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'source_mode = "reservoir_face"'
    write (u, '(a)') 'number_density_cm3 = 5.0'
    write (u, '(a)') 'temperature_ev = 10.0'
    write (u, '(a)') 'q_particle = -1.0'
    write (u, '(a)') 'm_particle = 1.0'
    write (u, '(a)') 'w_particle = 100.0'
    write (u, '(a)') 'inject_face = "z_low"'
    write (u, '(a)') 'pos_low = [0.0, 0.0, 0.0]'
    write (u, '(a)') 'pos_high = [1.0, 1.0, 0.0]'
    write (u, '(a)') 'drift_velocity = [0.0, 0.0, 1.0]'

    close (u)
  end subroutine write_fixed_duration_fixture

  !> テスト専用の `batch_duration_step` + species target 設定ファイルを書き出す。
  !! @param[in] path 書き出し先TOMLファイルパス。
  subroutine write_auto_duration_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open auto-duration reservoir config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'batch_count = 10'
    write (u, '(a)') 'dt = 1.0'
    write (u, '(a)') 'batch_duration_step = 3.0'
    write (u, '(a)') 'use_box = true'
    write (u, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (u, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'source_mode = "reservoir_face"'
    write (u, '(a)') 'number_density_m3 = 1000.0'
    write (u, '(a)') 'temperature_k = 0.0'
    write (u, '(a)') 'q_particle = -1.0'
    write (u, '(a)') 'm_particle = 1.0'
    write (u, '(a)') 'target_macro_particles_per_batch = 300'
    write (u, '(a)') 'inject_face = "z_high"'
    write (u, '(a)') 'pos_low = [0.0, 0.0, 1.0]'
    write (u, '(a)') 'pos_high = [1.0, 1.0, 1.0]'
    write (u, '(a)') 'drift_velocity = [0.0, 0.0, -1.0]'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'source_mode = "reservoir_face"'
    write (u, '(a)') 'number_density_m3 = 250.0'
    write (u, '(a)') 'temperature_k = 0.0'
    write (u, '(a)') 'q_particle = 1.0'
    write (u, '(a)') 'm_particle = 1.0'
    write (u, '(a)') 'target_macro_particles_per_batch = -1'
    write (u, '(a)') 'inject_face = "z_high"'
    write (u, '(a)') 'pos_low = [0.0, 0.0, 1.0]'
    write (u, '(a)') 'pos_high = [1.0, 1.0, 1.0]'
    write (u, '(a)') 'drift_velocity = [0.0, 0.0, -1.0]'

    close (u)
  end subroutine write_auto_duration_fixture

end program test_reservoir_injection

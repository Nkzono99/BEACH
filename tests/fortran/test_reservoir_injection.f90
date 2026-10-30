!> 境界 reservoir 流入を検証するテスト。
program test_reservoir_injection
  use bem_irradiation_cell_contract, only: irradiation_cell_contract_type
  use bem_kinds, only: dp, i32
  use bem_app_config, only: app_config, default_app_config, load_app_config, particles_per_batch_from_config, &
                            species_from_defaults, seed_particles_from_config, init_particle_batch_from_config
  use bem_app_config_runtime, only: particle_source_plan_type, build_particle_source_plan, &
                                    batch_injection_layout_type, build_batch_injection_layout, make_boundary_inflow_spec, &
                                    compute_face_average_potential, compute_wall_box_potential_statistics, &
                                    reservoir_face_velocity_correction
  use bem_app_config_types, only: particle_inflow_reservoir, particle_species_spec
  use bem_external_boundary_contract, only: external_inflow_scalar_barrier
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
  character(len=*), parameter :: cfg_boundary_path = 'test_boundary_inflow_tmp.toml'
  integer(i32) :: n_macro, n1, n2
  real(dp) :: residual
  real(dp) :: residual1, residual2
  real(dp) :: gamma1, area1, expected_w1
  real(dp) :: inward_normal(3)

  call test_init(9)

  call test_begin('snapshot_infinity_barrier')
  call test_snapshot_infinity_barrier()
  call test_end()

  call test_begin('layout_face_barriers_follow_each_face')
  call test_layout_face_barriers_follow_each_face()
  call test_end()

  call test_begin('face_potential_statistics_share_sampling_pass')
  call test_face_potential_statistics_share_sampling_pass()
  call test_end()

  call test_begin('boundary_inflow_full_faces_and_volume_source')
  call test_boundary_inflow_full_faces_and_volume_source()
  call test_end()

  call write_fixed_duration_fixture(cfg_fixed_path)

  call test_begin('fixed_duration_config')
  call default_app_config(cfg_fixed)
  call load_app_config(cfg_fixed_path, cfg_fixed)

  call assert_true(trim(cfg_fixed%particle_species(1)%source_mode) == 'volume_seed', 'source_mode mismatch')
  call assert_close_dp(cfg_fixed%sim%batch_duration, 1.0d0, 1.0d-12, 'batch_duration mismatch')
  call assert_close_dp(cfg_fixed%particle_species(1)%number_density_cm3, 5.0d0, 1.0d-12, 'density mismatch')
  call assert_close_dp(cfg_fixed%particle_species(1)%temperature_ev, 10.0d0, 1.0d-12, 'temperature_ev mismatch')
  call assert_true(cfg_fixed%particle_species(1)%boundary_inflow_low(3) == particle_inflow_reservoir, 'boundary inflow mismatch')
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
  residual = 0.0d0
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
  call assert_true(trim(cfg_auto%particle_species(1)%source_mode) == 'volume_seed', 'species-1 mode mismatch')
  call assert_true(trim(cfg_auto%particle_species(2)%source_mode) == 'volume_seed', 'species-2 mode mismatch')

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
  ! Both species share the weight that puts 300 macro particles of species 1 on the face per batch, so the
  ! quarter-density species 2 gets exactly 75. The budget is an integer, so one batch fixes the absolute counts
  ! and leaves no residual; the fractional carry is covered by deterministic_batches. cfg_auto is the
  ! configuration loaded and checked by auto_duration_config.
  residual1 = 0.0d0
  residual2 = 0.0d0
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
  call assert_equal_i32(n1, 300_i32, 'reservoir species-1 macro count mismatch')
  call assert_equal_i32(n2, 75_i32, 'reservoir species-2 macro count mismatch')
  call assert_close_dp(residual1, 0.0d0, 1.0d-9, 'integer species-1 budget must leave no residual')
  call assert_close_dp(residual2, 0.0d0, 1.0d-9, 'integer species-2 budget must leave no residual')
  call test_end()

  call delete_file_if_exists(cfg_fixed_path)
  call delete_file_if_exists(cfg_auto_path)
  call delete_file_if_exists(cfg_boundary_path)

  call test_summary()

contains

  subroutine test_boundary_inflow_full_faces_and_volume_source()
    type(app_config) :: boundary_cfg
    type(particle_source_plan_type) :: source_plan
    type(particles_soa) :: particles
    type(injection_state) :: state
    integer(i32) :: n_volume, n_x_low, n_z_high, particle_idx
    real(dp), parameter :: face_tol = 1.0e-6_dp

    call write_boundary_inflow_fixture(cfg_boundary_path)
    call default_app_config(boundary_cfg)
    call load_app_config(cfg_boundary_path, boundary_cfg)

    call assert_true( &
      trim(boundary_cfg%particle_species(1)%source_mode) == 'volume_seed', &
      'boundary inflow must coexist with volume_seed' &
      )
    call assert_equal_i32( &
      boundary_cfg%particle_species(1)%boundary_inflow_low(1), particle_inflow_reservoir, &
      'x-low boundary inflow mode mismatch' &
      )
    call assert_equal_i32( &
      boundary_cfg%particle_species(1)%boundary_inflow_high(3), particle_inflow_reservoir, &
      'z-high boundary inflow mode mismatch' &
      )

    allocate (state%macro_residual(1), state%boundary_macro_residual(6, 1))
    state%macro_residual = 0.0_dp
    state%boundary_macro_residual = 0.0_dp
    call seed_particles_from_config(boundary_cfg)
    call build_particle_source_plan(boundary_cfg, source_plan)
    call init_particle_batch_from_config(boundary_cfg, 1_i32, particles, state=state, source_plan=source_plan)

    call assert_equal_i32(particles%n, 20_i32, 'volume source plus two full boundary faces count mismatch')
    n_volume = 0_i32
    n_x_low = 0_i32
    n_z_high = 0_i32
    do particle_idx = 1_i32, particles%n
      if (all(abs(particles%x(:, particle_idx) - [1.0_dp, 1.5_dp, 2.0_dp]) < 1.0e-12_dp)) then
        n_volume = n_volume + 1_i32
      else if (particles%x(1, particle_idx) < face_tol) then
        n_x_low = n_x_low + 1_i32
        call assert_true(particles%v(1, particle_idx) > 0.0_dp, 'x-low inflow velocity must point inward')
        call assert_true( &
          particles%x(2, particle_idx) >= 0.0_dp .and. particles%x(2, particle_idx) <= 3.0_dp .and. &
          particles%x(3, particle_idx) >= 0.0_dp .and. particles%x(3, particle_idx) <= 4.0_dp, &
          'x-low inflow must sample the full box face' &
          )
      else if (particles%x(3, particle_idx) > 4.0_dp - face_tol) then
        n_z_high = n_z_high + 1_i32
        call assert_true(particles%v(3, particle_idx) < 0.0_dp, 'z-high inflow velocity must point inward')
        call assert_true( &
          particles%x(1, particle_idx) >= 0.0_dp .and. particles%x(1, particle_idx) <= 2.0_dp .and. &
          particles%x(2, particle_idx) >= 0.0_dp .and. particles%x(2, particle_idx) <= 3.0_dp, &
          'z-high inflow must sample the full box face' &
          )
      end if
    end do
    call assert_equal_i32(n_volume, 2_i32, 'volume_seed contribution mismatch')
    call assert_equal_i32(n_x_low, 12_i32, 'x-low full-face area count mismatch')
    call assert_equal_i32(n_z_high, 6_i32, 'z-high full-face area count mismatch')
    call assert_close_dp(sum(particles%q*particles%w), 20.0_dp, 1.0e-14_dp, 'boundary inflow charge budget mismatch')
    call assert_true(all(abs(state%macro_residual) < 1.0e-14_dp), 'source residual should stay zero')
    call assert_true(all(abs(state%boundary_macro_residual) < 1.0e-14_dp), 'boundary residuals should stay zero')
  end subroutine test_boundary_inflow_full_faces_and_volume_source

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
    stats_cfg%particle_species(1)%boundary_inflow_high(3) = particle_inflow_reservoir
    stats_cfg%particle_species(1)%species_key = 'face_variation_warning_test'
    stats_cfg%particle_species(1)%q_particle = 1.0_dp
    stats_cfg%particle_species(1)%m_particle = 1.0_dp
    stats_cfg%particle_species(1)%temperature_k = 0.0_dp
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

    call compute_wall_box_potential_statistics( &
      stats_mesh, stats_cfg%sim, snapshot, top_phi_mean, top_phi_std, top_phi_min, top_phi_max, &
      cell=irradiation_cell_contract_type() &
      )
    call assert_close_dp(top_phi_mean, phi_mean, 1.0e-14_dp, 'z-high box potential mean mismatch')
    call assert_close_dp(top_phi_std, phi_std, 1.0e-14_dp, 'z-high box potential standard deviation mismatch')
    call assert_close_dp(top_phi_min, phi_min, 1.0e-14_dp, 'z-high box potential minimum mismatch')
    call assert_close_dp(top_phi_max, phi_max, 1.0e-14_dp, 'z-high box potential maximum mismatch')

    stats_cfg%sim%reservoir_potential_model = 'infinity_barrier'
    stats_cfg%boundary_contract%inflow_map = external_inflow_scalar_barrier
    stats_cfg%sim%phi_infty = -1.0_dp
    call reservoir_face_velocity_correction( &
      stats_cfg, stats_cfg%particle_species(1), vmin_normal, barrier_normal, &
      mesh=stats_mesh, snapshot=snapshot, warn_face_variation=.true. &
      )
    call assert_close_dp(barrier_normal, 1.0_dp, 1.0e-14_dp, 'face-average barrier mismatch')
    call assert_close_dp(vmin_normal, 1.0_dp, 1.0e-14_dp, 'face-average cutoff mismatch')
  end subroutine test_face_potential_statistics_share_sampling_pass

  !> 一様場 E=(0,0,-1) では z_low と z_high の面平均電位が違う。
  !! 二つの種が両面から流入する layout の障壁は、種ごと・面ごとに、
  !! その面だけを評価した面平均電位から求めた値と一致しなければならない。
  subroutine test_layout_face_barriers_follow_each_face()
    type(app_config) :: cfg
    type(particle_source_plan_type) :: source_plan
    type(batch_injection_layout_type) :: layout
    type(mesh_type) :: mesh
    type(injection_state) :: state
    type(electrostatic_snapshot_type) :: snapshot
    type(particle_species_spec) :: face_spec
    real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1), phi_face, expected
    integer :: s, face

    v0(:, 1) = [0.0_dp, 0.0_dp, 0.5_dp]
    v1(:, 1) = [1.0_dp, 0.0_dp, 0.5_dp]
    v2(:, 1) = [0.0_dp, 1.0_dp, 0.5_dp]
    call init_mesh(mesh, v0, v1, v2, q0=[0.0_dp])
    mesh%elem_vacuum_sign = 1_i32
    mesh%vacuum_normals = mesh%normals
    call default_app_config(cfg)
    cfg%sim%batch_count = 1_i32
    cfg%sim%batch_duration = 1.0_dp
    cfg%sim%has_batch_duration = .true.
    cfg%sim%use_box = .true.
    cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    cfg%sim%e0 = [0.0_dp, 0.0_dp, -1.0_dp]
    cfg%sim%reservoir_potential_model = 'infinity_barrier'
    cfg%boundary_contract%inflow_map = external_inflow_scalar_barrier
    cfg%sim%phi_infty = 0.25_dp
    cfg%sim%injection_face_phi_grid_n = 2_i32
    cfg%n_particle_species = 2_i32
    do s = 1, 2
      cfg%particle_species(s) = species_from_defaults()
      cfg%particle_species(s)%source_mode = 'volume_seed'
      cfg%particle_species(s)%npcls_per_step = 0_i32
      cfg%particle_species(s)%number_density_m3 = 10.0_dp
      cfg%particle_species(s)%temperature_k = 0.0_dp
      cfg%particle_species(s)%q_particle = real(3 - 2*s, dp)
      cfg%particle_species(s)%m_particle = real(s, dp)
      cfg%particle_species(s)%w_particle = 1.0_dp
      cfg%particle_species(s)%has_w_particle = .true.
      cfg%particle_species(s)%boundary_inflow_low(3) = particle_inflow_reservoir
      cfg%particle_species(s)%boundary_inflow_high(3) = particle_inflow_reservoir
    end do
    allocate (state%macro_residual(2), state%boundary_macro_residual(6, 2))
    state%macro_residual = 0.0_dp
    state%boundary_macro_residual = 0.0_dp
    call snapshot%init(mesh, cfg%sim, cfg%field, cfg%periodic2, cfg%panel)
    call snapshot%refresh(mesh)
    call build_particle_source_plan(cfg, source_plan)
    call build_batch_injection_layout( &
      cfg, 1_i32, layout, state=state, mesh=mesh, snapshot=snapshot, source_plan=source_plan &
      )
    do s = 1, 2
      do face = 5, 6
        call make_boundary_inflow_spec(cfg%sim, cfg%particle_species(s), int(face, i32), face_spec)
        call compute_face_average_potential(mesh, cfg%sim, face_spec, snapshot, phi_face)
        expected = 2.0_dp*cfg%particle_species(s)%q_particle*(phi_face - cfg%sim%phi_infty)/ &
                   cfg%particle_species(s)%m_particle
        call assert_close_dp(layout%boundary_barrier(face, s), expected, 1.0e-14_dp, &
                             'layout barrier must use the potential of its own face')
      end do
    end do
    call assert_true(abs(layout%boundary_barrier(5, 1) - layout%boundary_barrier(6, 1)) > 0.1_dp, &
                     'the two faces must see different potentials in this fixture')
  end subroutine test_layout_face_barriers_follow_each_face

  subroutine test_snapshot_infinity_barrier()
    type(app_config) :: barrier_cfg
    type(particle_source_plan_type) :: source_plan
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
    barrier_cfg%boundary_contract%inflow_map = external_inflow_scalar_barrier
    barrier_cfg%sim%phi_infty = 0.0_dp
    barrier_cfg%sim%injection_face_phi_grid_n = 2_i32
    barrier_cfg%n_particle_species = 1_i32
    barrier_cfg%particle_species(1) = species_from_defaults()
    barrier_cfg%particle_species(1)%source_mode = 'volume_seed'
    barrier_cfg%particle_species(1)%npcls_per_step = 0_i32
    barrier_cfg%particle_species(1)%number_density_m3 = 10.0_dp
    barrier_cfg%particle_species(1)%temperature_k = 0.0_dp
    barrier_cfg%particle_species(1)%q_particle = 1.0_dp
    barrier_cfg%particle_species(1)%m_particle = 1.0_dp
    barrier_cfg%particle_species(1)%w_particle = 1.0_dp
    barrier_cfg%particle_species(1)%has_w_particle = .true.
    barrier_cfg%particle_species(1)%boundary_inflow_high(3) = particle_inflow_reservoir
    barrier_cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -1.0_dp]
    allocate (state%macro_residual(1), state%boundary_macro_residual(6, 1))
    state%macro_residual = 0.0_dp
    state%boundary_macro_residual = 0.0_dp

    call snapshot%init( &
      barrier_mesh, barrier_cfg%sim, barrier_cfg%field, barrier_cfg%periodic2, barrier_cfg%panel &
      )
    call snapshot%refresh(barrier_mesh)
    call build_particle_source_plan(barrier_cfg, source_plan)
    call init_particle_batch_from_config( &
      barrier_cfg, 1_i32, particles, state=state, mesh=barrier_mesh, snapshot=snapshot, source_plan=source_plan &
      )
    call assert_equal_i32(particles%n, 0_i32, 'snapshot infinity barrier should block boundary inflow')
  end subroutine test_snapshot_infinity_barrier

  subroutine write_boundary_inflow_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open boundary inflow config fixture'

    write (u, '(a)') '[run]'
    write (u, '(a)') 'batch_count = 1'
    write (u, '(a)') 'rng_seed = 24680'
    write (u, '(a)') ''
    write (u, '(a)') '[run.batch]'
    write (u, '(a)') 'duration_s = 1.0'
    write (u, '(a)') ''
    write (u, '(a)') '[domain]'
    write (u, '(a)') 'box_min = ['
    write (u, '(a)') '    0.0,'
    write (u, '(a)') '    0.0,'
    write (u, '(a)') '    0.0,'
    write (u, '(a)') ']'
    write (u, '(a)') 'box_max = ['
    write (u, '(a)') '    2.0,'
    write (u, '(a)') '    3.0,'
    write (u, '(a)') '    4.0,'
    write (u, '(a)') ']'
    write (u, '(a)') 'periodic_axes = []'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.tracking]'
    write (u, '(a)') 'dt_s = 1e-12'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.boundary]'
    write (u, '(a)') 'x_low = "open"'
    write (u, '(a)') 'z_high = "open"'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'charge_c = 1.0'
    write (u, '(a)') 'mass_kg = 1.0'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.species.source]'
    write (u, '(a)') 'mode = "volume_seed"'
    write (u, '(a)') 'pos_low = ['
    write (u, '(a)') '    1.0,'
    write (u, '(a)') '    1.5,'
    write (u, '(a)') '    2.0,'
    write (u, '(a)') ']'
    write (u, '(a)') 'pos_high = ['
    write (u, '(a)') '    1.0,'
    write (u, '(a)') '    1.5,'
    write (u, '(a)') '    2.0,'
    write (u, '(a)') ']'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.species.sampling]'
    write (u, '(a)') 'volume_macro_particles_per_batch = 2'
    write (u, '(a)') 'weight = 1.0'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.species.distribution]'
    write (u, '(a)') 'number_density_m3 = 1.0'
    write (u, '(a)') 'temperature_k = 0.0'
    write (u, '(a)') 'drift_velocity_m_s = ['
    write (u, '(a)') '    1.0,'
    write (u, '(a)') '    0.0,'
    write (u, '(a)') '    -1.0,'
    write (u, '(a)') ']'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.species.inflow]'
    write (u, '(a)') 'x_low = "reservoir"'
    write (u, '(a)') 'z_high = "reservoir"'
    write (u, '(a)') ''
    write (u, '(a)') '[fields]'
    write (u, '(a)') 'boundary = "free"'
    close (u)
  end subroutine write_boundary_inflow_fixture

  !> テスト専用の固定 `batch_duration` 境界流入設定ファイルを書き出す。
  !! @param[in] path 書き出し先TOMLファイルパス。
  subroutine write_fixed_duration_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open reservoir config fixture'

    write (u, '(a)') '[run]'
    write (u, '(a)') 'batch_count = 2'
    write (u, '(a)') ''
    write (u, '(a)') '[run.batch]'
    write (u, '(a)') 'duration_s = 1.0'
    write (u, '(a)') ''
    write (u, '(a)') '[domain]'
    write (u, '(a)') 'box_min = ['
    write (u, '(a)') '    0.0,'
    write (u, '(a)') '    0.0,'
    write (u, '(a)') '    0.0,'
    write (u, '(a)') ']'
    write (u, '(a)') 'box_max = ['
    write (u, '(a)') '    1.0,'
    write (u, '(a)') '    1.0,'
    write (u, '(a)') '    1.0,'
    write (u, '(a)') ']'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'charge_c = -1.0'
    write (u, '(a)') 'mass_kg = 1.0'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.species.distribution]'
    write (u, '(a)') 'number_density_cm3 = 5.0'
    write (u, '(a)') 'temperature_ev = 10.0'
    write (u, '(a)') 'drift_velocity_m_s = ['
    write (u, '(a)') '    0.0,'
    write (u, '(a)') '    0.0,'
    write (u, '(a)') '    1.0,'
    write (u, '(a)') ']'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.species.sampling]'
    write (u, '(a)') 'weight = 100.0'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.species.inflow]'
    write (u, '(a)') 'z_low = "reservoir"'

    close (u)
  end subroutine write_fixed_duration_fixture

  !> テスト専用の `batch_duration_step` + species target 設定ファイルを書き出す。
  !! @param[in] path 書き出し先TOMLファイルパス。
  subroutine write_auto_duration_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open auto-duration reservoir config fixture'

    write (u, '(a)') '[run]'
    write (u, '(a)') 'batch_count = 10'
    write (u, '(a)') ''
    write (u, '(a)') '[run.batch]'
    write (u, '(a)') 'duration_steps = 3.0'
    write (u, '(a)') ''
    write (u, '(a)') '[domain]'
    write (u, '(a)') 'box_min = ['
    write (u, '(a)') '    0.0,'
    write (u, '(a)') '    0.0,'
    write (u, '(a)') '    0.0,'
    write (u, '(a)') ']'
    write (u, '(a)') 'box_max = ['
    write (u, '(a)') '    1.0,'
    write (u, '(a)') '    1.0,'
    write (u, '(a)') '    1.0,'
    write (u, '(a)') ']'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.tracking]'
    write (u, '(a)') 'dt_s = 1.0'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'charge_c = -1.0'
    write (u, '(a)') 'mass_kg = 1.0'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.species.distribution]'
    write (u, '(a)') 'number_density_m3 = 1000.0'
    write (u, '(a)') 'temperature_k = 0.0'
    write (u, '(a)') 'drift_velocity_m_s = ['
    write (u, '(a)') '    0.0,'
    write (u, '(a)') '    0.0,'
    write (u, '(a)') '    -1.0,'
    write (u, '(a)') ']'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.species.sampling]'
    write (u, '(a)') 'target_macro_particles_per_batch = 300'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.species.inflow]'
    write (u, '(a)') 'z_high = "reservoir"'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'charge_c = 1.0'
    write (u, '(a)') 'mass_kg = 1.0'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.species.distribution]'
    write (u, '(a)') 'number_density_m3 = 250.0'
    write (u, '(a)') 'temperature_k = 0.0'
    write (u, '(a)') 'drift_velocity_m_s = ['
    write (u, '(a)') '    0.0,'
    write (u, '(a)') '    0.0,'
    write (u, '(a)') '    -1.0,'
    write (u, '(a)') ']'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.species.sampling]'
    write (u, '(a)') 'target_macro_particles_per_batch = -1'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.species.inflow]'
    write (u, '(a)') 'z_high = "reservoir"'

    close (u)
  end subroutine write_auto_duration_fixture

end program test_reservoir_injection

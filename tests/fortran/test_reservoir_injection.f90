!> 境界 reservoir・内部plane source・legacy reservoir_face注入を検証するテスト。
program test_reservoir_injection
  use bem_kinds, only: dp, i32
  use bem_app_config, only: app_config, default_app_config, load_app_config, particles_per_batch_from_config, &
                            species_from_defaults, seed_particles_from_config, init_particle_batch_from_config
  use bem_app_config_runtime, only: particle_source_plan_type, build_particle_source_plan, &
                                    compute_face_average_potential, compute_z_high_box_potential_statistics, &
                                    reservoir_face_velocity_correction
  use bem_app_config_types, only: particle_inflow_reservoir
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
  character(len=*), parameter :: cfg_plane_path = 'test_plane_source_tmp.toml'
  character(len=*), parameter :: cfg_invalid_periodic_path = 'test_boundary_inflow_periodic_tmp.toml'
  character(len=*), parameter :: cfg_invalid_reflect_path = 'test_boundary_inflow_reflect_tmp.toml'
  character(len=*), parameter :: cfg_invalid_combined_path = 'test_boundary_inflow_combined_tmp.toml'
  character(len=*), parameter :: config_failure_path = 'test_boundary_inflow_failure_tmp.log'
  integer(i32) :: n_macro, i, n1, n2
  integer(i32) :: sum1, sum2
  real(dp) :: residual
  real(dp) :: residual1, residual2, ratio
  real(dp) :: gamma1, area1, expected_w1
  real(dp) :: inward_normal(3)
  character(len=64) :: run_mode
  character(len=512) :: probe_config_path

  call get_command_argument(1, run_mode)
  if (trim(run_mode) == '--config-failure-probe') then
    call get_command_argument(2, probe_config_path)
    call run_config_failure_probe(trim(probe_config_path))
    error stop 'invalid boundary inflow config probe unexpectedly completed'
  end if

  call test_init(11)

  call test_begin('particle_source_plan_equivalence')
  call test_particle_source_plan_equivalence()
  call test_end()

  call test_begin('snapshot_infinity_barrier')
  call test_snapshot_infinity_barrier()
  call test_end()

  call test_begin('face_potential_statistics_share_sampling_pass')
  call test_face_potential_statistics_share_sampling_pass()
  call test_end()

  call test_begin('boundary_inflow_full_faces_and_volume_source')
  call test_boundary_inflow_full_faces_and_volume_source()
  call test_end()

  call test_begin('plane_source_internal_face')
  call test_plane_source_internal_face()
  call test_end()

  call test_begin('boundary_inflow_invalid_combinations')
  call test_boundary_inflow_invalid_combinations()
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
  call delete_file_if_exists(cfg_boundary_path)
  call delete_file_if_exists(cfg_plane_path)
  call delete_file_if_exists(cfg_invalid_periodic_path)
  call delete_file_if_exists(cfg_invalid_reflect_path)
  call delete_file_if_exists(cfg_invalid_combined_path)
  call delete_file_if_exists(config_failure_path)

  call test_summary()

contains

  subroutine test_boundary_inflow_full_faces_and_volume_source()
    type(app_config) :: boundary_cfg
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
    call init_particle_batch_from_config(boundary_cfg, 1_i32, particles, state=state)

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
    call assert_true(all(abs(state%macro_residual) < 1.0e-14_dp), 'source residual should stay zero')
    call assert_true(all(abs(state%boundary_macro_residual) < 1.0e-14_dp), 'boundary residuals should stay zero')
  end subroutine test_boundary_inflow_full_faces_and_volume_source

  subroutine test_plane_source_internal_face()
    type(app_config) :: plane_cfg
    type(particles_soa) :: particles
    type(injection_state) :: state

    call write_plane_source_fixture(cfg_plane_path)
    call default_app_config(plane_cfg)
    call load_app_config(cfg_plane_path, plane_cfg)

    call assert_true(trim(plane_cfg%particle_species(1)%source_mode) == 'plane_source', 'plane source mode mismatch')
    call assert_true(trim(plane_cfg%particle_species(1)%inject_face) == 'x_low', 'positive source normal face mismatch')
    call assert_close_dp(plane_cfg%particle_species(1)%source_normal(1), 1.0_dp, 1.0e-14_dp, &
                         'source normal must be normalized')

    allocate (state%macro_residual(1))
    state%macro_residual = 0.0_dp
    call seed_particles_from_config(plane_cfg)
    call init_particle_batch_from_config(plane_cfg, 1_i32, particles, state=state)

    call assert_equal_i32(particles%n, 8_i32, 'plane source flux count mismatch')
    call assert_true(all(particles%x(1, :) > 0.5_dp), 'plane source positions must be jittered along source_normal')
    call assert_true(all(particles%x(1, :) < 0.5_dp + 1.0e-6_dp), 'plane source normal jitter exceeded dt bound')
    call assert_true( &
      all(particles%x(2, :) >= 0.0_dp .and. particles%x(2, :) <= 2.0_dp), &
      'plane source y positions must stay on the configured rectangle' &
      )
    call assert_true( &
      all(particles%x(3, :) >= 0.0_dp .and. particles%x(3, :) <= 3.0_dp), &
      'plane source z positions must stay on the configured rectangle' &
      )
    call assert_true(all(particles%v(1, :) > 0.0_dp), 'plane source velocity must follow source_normal')
    call assert_true(all(abs(particles%v(2:3, :)) < 1.0e-14_dp), 'cold plane source tangential velocity mismatch')
  end subroutine test_plane_source_internal_face

  subroutine test_boundary_inflow_invalid_combinations()
    call write_invalid_boundary_inflow_fixture(cfg_invalid_periodic_path, 'periodic')
    call assert_config_rejected(cfg_invalid_periodic_path, 'cannot inject through a periodic domain face')

    call write_invalid_boundary_inflow_fixture(cfg_invalid_reflect_path, 'reflect')
    call assert_config_rejected(cfg_invalid_reflect_path, 'requires the effective particle action to be open')

    call write_invalid_boundary_inflow_fixture(cfg_invalid_combined_path, 'combined')
    call assert_config_rejected( &
      cfg_invalid_combined_path, 'source_mode="plane_source" cannot be combined with boundary_inflow' &
      )
  end subroutine test_boundary_inflow_invalid_combinations

  subroutine run_config_failure_probe(path)
    character(len=*), intent(in) :: path
    type(app_config) :: probe_cfg

    call default_app_config(probe_cfg)
    call load_app_config(path, probe_cfg)
  end subroutine run_config_failure_probe

  subroutine assert_config_rejected(path, expected_fragment)
    character(len=*), intent(in) :: path, expected_fragment
    character(len=1024) :: executable_path, command, child_line
    integer :: child_exit_status, child_cmd_status, child_unit, child_ios
    logical :: saw_expected

    call get_command_argument(0, executable_path)
    call delete_file_if_exists(config_failure_path)
    command = '"'//trim(executable_path)//'" --config-failure-probe "'//trim(path)//'" > '// &
              config_failure_path//' 2>&1'
    call execute_command_line(trim(command), wait=.true., exitstat=child_exit_status, cmdstat=child_cmd_status)
    call assert_equal_i32(int(child_cmd_status, i32), 0_i32, 'config failure probe command status mismatch')
    call assert_true(child_exit_status /= 0, 'invalid boundary inflow config must be rejected')

    saw_expected = .false.
    open (newunit=child_unit, file=config_failure_path, status='old', action='read', iostat=child_ios)
    if (child_ios /= 0) error stop 'failed to read boundary inflow failure probe output'
    do
      read (child_unit, '(A)', iostat=child_ios) child_line
      if (child_ios /= 0) exit
      saw_expected = saw_expected .or. index(child_line, trim(expected_fragment)) > 0
    end do
    close (child_unit)
    call delete_file_if_exists(config_failure_path)
    call assert_true(saw_expected, 'config failure message mismatch: '//trim(expected_fragment))
  end subroutine assert_config_rejected

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
    barrier_cfg%particle_species(1)%source_mode = 'volume_seed'
    barrier_cfg%particle_species(1)%npcls_per_step = 0_i32
    barrier_cfg%particle_species(1)%number_density_m3 = 10.0_dp
    barrier_cfg%particle_species(1)%has_number_density_m3 = .true.
    barrier_cfg%particle_species(1)%temperature_k = 0.0_dp
    barrier_cfg%particle_species(1)%has_temperature_k = .true.
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
    call init_particle_batch_from_config( &
      barrier_cfg, 1_i32, particles, state=state, mesh=barrier_mesh, snapshot=snapshot &
      )
    call assert_equal_i32(particles%n, 0_i32, 'snapshot infinity barrier should block boundary inflow')
  end subroutine test_snapshot_infinity_barrier

  subroutine write_boundary_inflow_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open boundary inflow config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'dt = 1.0e-12'
    write (u, '(a)') 'batch_count = 1'
    write (u, '(a)') 'batch_duration = 1.0'
    write (u, '(a)') 'rng_seed = 24680'
    write (u, '(a)') ''
    write (u, '(a)') '[domain]'
    write (u, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (u, '(a)') 'box_max = [2.0, 3.0, 4.0]'
    write (u, '(a)') 'periodic_axes = []'
    write (u, '(a)') ''
    write (u, '(a)') '[field_boundary]'
    write (u, '(a)') 'mode = "free"'
    write (u, '(a)') ''
    write (u, '(a)') '[particle_boundary]'
    write (u, '(a)') 'x_low = "open"'
    write (u, '(a)') 'z_high = "open"'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'source_mode = "volume_seed"'
    write (u, '(a)') 'npcls_per_step = 2'
    write (u, '(a)') 'number_density_m3 = 1.0'
    write (u, '(a)') 'temperature_k = 0.0'
    write (u, '(a)') 'q_particle = 1.0'
    write (u, '(a)') 'm_particle = 1.0'
    write (u, '(a)') 'w_particle = 1.0'
    write (u, '(a)') 'pos_low = [1.0, 1.5, 2.0]'
    write (u, '(a)') 'pos_high = [1.0, 1.5, 2.0]'
    write (u, '(a)') 'drift_velocity = [1.0, 0.0, -1.0]'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.species.boundary_inflow]'
    write (u, '(a)') 'x_low = "reservoir"'
    write (u, '(a)') 'z_high = "reservoir"'
    close (u)
  end subroutine write_boundary_inflow_fixture

  subroutine write_plane_source_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open plane source config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'dt = 1.0e-12'
    write (u, '(a)') 'batch_count = 1'
    write (u, '(a)') 'batch_duration = 1.0'
    write (u, '(a)') 'rng_seed = 13579'
    write (u, '(a)') ''
    write (u, '(a)') '[domain]'
    write (u, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (u, '(a)') 'box_max = [2.0, 2.0, 3.0]'
    write (u, '(a)') 'periodic_axes = []'
    write (u, '(a)') ''
    write (u, '(a)') '[field_boundary]'
    write (u, '(a)') 'mode = "free"'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'source_mode = "plane_source"'
    write (u, '(a)') 'number_density_m3 = 2.0'
    write (u, '(a)') 'temperature_k = 0.0'
    write (u, '(a)') 'q_particle = 1.0'
    write (u, '(a)') 'm_particle = 1.0'
    write (u, '(a)') 'w_particle = 3.0'
    write (u, '(a)') 'pos_low = [0.5, 0.0, 0.0]'
    write (u, '(a)') 'pos_high = [0.5, 2.0, 3.0]'
    write (u, '(a)') 'source_normal = [2.0, 0.0, 0.0]'
    write (u, '(a)') 'drift_velocity = [2.0, 0.0, 0.0]'
    close (u)
  end subroutine write_plane_source_fixture

  subroutine write_invalid_boundary_inflow_fixture(path, mode)
    character(len=*), intent(in) :: path, mode
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open invalid boundary inflow config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'dt = 1.0e-12'
    write (u, '(a)') 'batch_count = 1'
    write (u, '(a)') 'batch_duration = 1.0'
    if (trim(mode) == 'periodic') write (u, '(a)') 'field_solver = "fmm"'
    write (u, '(a)') ''
    write (u, '(a)') '[domain]'
    write (u, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (u, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    if (trim(mode) == 'periodic') then
      write (u, '(a)') 'periodic_axes = ["x", "y"]'
    else
      write (u, '(a)') 'periodic_axes = []'
    end if
    write (u, '(a)') ''
    write (u, '(a)') '[field_boundary]'
    if (trim(mode) == 'periodic') then
      write (u, '(a)') 'mode = "periodic2"'
    else
      write (u, '(a)') 'mode = "free"'
    end if
    write (u, '(a)') ''
    write (u, '(a)') '[particle_boundary]'
    if (trim(mode) == 'reflect') then
      write (u, '(a)') 'x_low = "reflect"'
    else if (trim(mode) == 'combined') then
      write (u, '(a)') 'z_high = "open"'
    end if
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    if (trim(mode) == 'combined') then
      write (u, '(a)') 'source_mode = "plane_source"'
      write (u, '(a)') 'pos_low = [0.5, 0.0, 0.0]'
      write (u, '(a)') 'pos_high = [0.5, 1.0, 1.0]'
      write (u, '(a)') 'source_normal = [1.0, 0.0, 0.0]'
    else
      write (u, '(a)') 'source_mode = "volume_seed"'
      write (u, '(a)') 'npcls_per_step = 0'
    end if
    write (u, '(a)') 'number_density_m3 = 1.0'
    write (u, '(a)') 'temperature_k = 0.0'
    write (u, '(a)') 'q_particle = 1.0'
    write (u, '(a)') 'm_particle = 1.0'
    write (u, '(a)') 'w_particle = 1.0'
    write (u, '(a)') 'drift_velocity = [1.0, 0.0, -1.0]'
    write (u, '(a)') ''
    write (u, '(a)') '[particles.species.boundary_inflow]'
    if (trim(mode) == 'combined') then
      write (u, '(a)') 'z_high = "reservoir"'
    else
      write (u, '(a)') 'x_low = "reservoir"'
    end if
    close (u)
  end subroutine write_invalid_boundary_inflow_fixture

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

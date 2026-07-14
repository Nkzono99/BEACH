!> reservoir_face 注入の設定解釈とマクロ粒子数計算を検証するテスト。
program test_reservoir_injection
  use bem_kinds, only: dp, i32
  use bem_app_config, only: app_config, default_app_config, load_app_config, particles_per_batch_from_config, &
                            species_from_defaults, init_particle_batch_from_config
  use bem_app_config_runtime, only: compute_face_average_potential, reservoir_face_velocity_correction
  use bem_injection, only: compute_macro_particles_for_batch, &
                           compute_inflow_flux_from_drifting_maxwellian, compute_face_area_from_bounds
  use bem_types, only: particles_soa, injection_state
  use bem_types, only: mesh_type
  use bem_mesh, only: init_mesh
  use bem_constants, only: eps0
  use bem_outer_plasma_types, only: outer_plasma_state_type
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

  call test_init(9)

  call test_begin('split_outer_infinity_vdf_map')
  call test_split_outer_infinity_vdf_map()
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

  call test_begin('global_count_independent_of_mpi_size')
  call test_global_count_independent_of_mpi_size()
  call test_end()

  call delete_file_if_exists(cfg_fixed_path)
  call delete_file_if_exists(cfg_auto_path)

  call test_summary()

contains

  subroutine test_face_potential_statistics_share_sampling_pass()
    type(app_config) :: stats_cfg
    type(mesh_type) :: stats_mesh
    type(electrostatic_snapshot_type) :: snapshot
    real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1)
    real(dp) :: phi_mean, phi_std, phi_min, phi_max, vmin_normal, barrier_normal

    v0(:, 1) = [0.0_dp, 0.0_dp, 0.25_dp]
    v1(:, 1) = [1.0_dp, 0.0_dp, 0.25_dp]
    v2(:, 1) = [0.0_dp, 1.0_dp, 0.25_dp]
    call init_mesh(stats_mesh, v0, v1, v2, q0=[0.0_dp])

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
      stats_mesh, stats_cfg%sim, stats_cfg%field, stats_cfg%periodic2, stats_cfg%panel, stats_cfg%outer_plasma &
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
      barrier_mesh, barrier_cfg%sim, barrier_cfg%field, barrier_cfg%periodic2, barrier_cfg%panel, &
      barrier_cfg%outer_plasma &
      )
    call snapshot%refresh(barrier_mesh)
    call init_particle_batch_from_config( &
      barrier_cfg, 1_i32, particles, state=state, mesh=barrier_mesh, snapshot=snapshot &
      )
    call assert_equal_i32(particles%n, 0_i32, 'snapshot infinity barrier should block deterministic inflow')
  end subroutine test_snapshot_infinity_barrier

  subroutine test_split_outer_infinity_vdf_map()
    type(app_config) :: split_cfg
    type(mesh_type) :: split_mesh
    type(particles_soa) :: particles
    type(injection_state) :: state
    type(outer_plasma_state_type) :: outer_state
    real(dp) :: v0(3, 2), v1(3, 2), v2(3, 2)

    v0(:, 1) = [0.0_dp, 0.0_dp, 0.25_dp]
    v1(:, 1) = [1.0_dp, 0.0_dp, 0.25_dp]
    v2(:, 1) = [1.0_dp, 1.0_dp, 0.25_dp]
    v0(:, 2) = [0.0_dp, 0.0_dp, 0.25_dp]
    v1(:, 2) = [1.0_dp, 1.0_dp, 0.25_dp]
    v2(:, 2) = [0.0_dp, 1.0_dp, 0.25_dp]
    call init_mesh(split_mesh, v0, v1, v2, q0=[-2.5_dp*eps0, -2.5_dp*eps0])
    call default_app_config(split_cfg)
    split_cfg%sim%batch_duration = 1.0_dp
    split_cfg%sim%has_batch_duration = .true.
    split_cfg%sim%dt = 0.0_dp
    split_cfg%sim%use_box = .true.
    split_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    split_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    split_cfg%outer_plasma%model = 'kinetic_1d'
    split_cfg%outer_plasma%debye_length = 0.2_dp
    split_cfg%coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
    split_cfg%n_particle_species = 1_i32
    split_cfg%particle_species(1) = species_from_defaults()
    split_cfg%particle_species(1)%source_mode = 'reservoir_face'
    split_cfg%particle_species(1)%number_density_m3 = 1.0_dp
    split_cfg%particle_species(1)%has_number_density_m3 = .true.
    split_cfg%particle_species(1)%temperature_k = 0.0_dp
    split_cfg%particle_species(1)%has_temperature_k = .true.
    split_cfg%particle_species(1)%q_particle = 1.0_dp
    split_cfg%particle_species(1)%m_particle = 1.0_dp
    split_cfg%particle_species(1)%w_particle = 0.1_dp
    split_cfg%particle_species(1)%has_w_particle = .true.
    split_cfg%particle_species(1)%inject_face = 'z_high'
    split_cfg%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, 1.0_dp]
    split_cfg%particle_species(1)%pos_high = [1.0_dp, 1.0_dp, 1.0_dp]
    split_cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -1.0_dp]
    outer_state%model = 'kinetic_1d'
    outer_state%ready = .true.
    outer_state%interface_potential = -2.0_dp
    outer_state%infinity_potential = 0.0_dp
    allocate (state%macro_residual(1))
    state%macro_residual = 0.0_dp
    call init_particle_batch_from_config( &
      split_cfg, 1_i32, particles, state=state, mesh=split_mesh, outer_state=outer_state &
      )
    call assert_equal_i32(particles%n, 10_i32, 'mapped ambient macro count mismatch')
    call assert_close_dp(particles%v(3, 1), -sqrt(5.0_dp), 1.0e-14_dp, 'kinetic interface velocity mismatch')
  end subroutine test_split_outer_infinity_vdf_map

  subroutine test_global_count_independent_of_mpi_size()
    integer(i32), parameter :: rank_sizes(3) = [1_i32, 2_i32, 4_i32]
    integer(i32), parameter :: expected_counts(4) = [0_i32, 0_i32, 0_i32, 1_i32]
    real(dp), parameter :: expected_residuals(4) = [0.25_dp, 0.5_dp, 0.75_dp, 0.0_dp]
    type(app_config) :: count_cfg
    type(particles_soa) :: particles
    type(injection_state), allocatable :: states(:)
    type(outer_plasma_state_type) :: outer_state
    integer(i32) :: size_index, n_ranks, rank, batch_idx, global_count

    call default_app_config(count_cfg)
    count_cfg%sim%batch_count = 4_i32
    count_cfg%sim%batch_duration = 1.0_dp
    count_cfg%sim%has_batch_duration = .true.
    count_cfg%sim%use_box = .true.
    count_cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    count_cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    count_cfg%outer_plasma%model = 'kinetic_1d'
    count_cfg%coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
    count_cfg%n_particle_species = 1_i32
    count_cfg%particle_species(1) = species_from_defaults()
    count_cfg%particle_species(1)%source_mode = 'reservoir_face'
    count_cfg%particle_species(1)%number_density_m3 = 1.0_dp
    count_cfg%particle_species(1)%has_number_density_m3 = .true.
    count_cfg%particle_species(1)%temperature_k = 0.0_dp
    count_cfg%particle_species(1)%has_temperature_k = .true.
    count_cfg%particle_species(1)%q_particle = 1.0_dp
    count_cfg%particle_species(1)%m_particle = 1.0_dp
    count_cfg%particle_species(1)%w_particle = 4.0_dp
    count_cfg%particle_species(1)%has_w_particle = .true.
    count_cfg%particle_species(1)%inject_face = 'z_high'
    count_cfg%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, 1.0_dp]
    count_cfg%particle_species(1)%pos_high = [1.0_dp, 1.0_dp, 1.0_dp]
    count_cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -1.0_dp]
    outer_state%model = 'kinetic_1d'
    outer_state%ready = .true.
    outer_state%interface_potential = 0.125_dp
    outer_state%infinity_potential = 0.0_dp

    do size_index = 1_i32, size(rank_sizes)
      n_ranks = rank_sizes(size_index)
      allocate (states(n_ranks))
      do rank = 1_i32, n_ranks
        allocate (states(rank)%macro_residual(1))
        states(rank)%macro_residual = 0.0_dp
      end do
      do batch_idx = 1_i32, 4_i32
        global_count = 0_i32
        do rank = 1_i32, n_ranks
          call init_particle_batch_from_config( &
            count_cfg, batch_idx, particles, state=states(rank), outer_state=outer_state, &
            mpi_rank=rank - 1_i32, mpi_size=n_ranks &
            )
          global_count = global_count + particles%n
          call assert_close_dp( &
            states(rank)%macro_residual(1), expected_residuals(batch_idx), 1.0e-15_dp, &
            'global reservoir residual must match on every synthetic rank' &
            )
        end do
        call assert_equal_i32( &
          global_count, expected_counts(batch_idx), 'global reservoir count must not depend on MPI size' &
          )
        if (batch_idx == 4_i32 .and. n_ranks == 1_i32) then
          call assert_close_dp(particles%v(3, 1), -sqrt(0.75_dp), 1.0e-14_dp, &
                               'kinetic barrier velocity must not depend on MPI layout')
        end if
      end do
      deallocate (states)
    end do
  end subroutine test_global_count_independent_of_mpi_size

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

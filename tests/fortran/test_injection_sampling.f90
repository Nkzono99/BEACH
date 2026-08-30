!> 粒子注入モジュールのサンプリング分岐を重点的に検証するテスト。
program test_injection_sampling
  use bem_kinds, only: dp, i32
  use bem_constants, only: k_boltzmann
  use bem_types, only: particles_soa, mesh_type, sim_config, bc_open, bc_reflect, bc_periodic
  use bem_mesh, only: init_mesh, prepare_periodic2_collision_mesh
  use bem_injection, only: &
    seed_rng, sample_shifted_maxwell_velocities, init_random_beam_particles, &
    compute_inflow_flux_from_drifting_maxwellian, compute_face_area_from_bounds, &
    sample_reservoir_face_particles, compute_macro_particles_for_batch, sample_photo_raycast_particles
  use bem_collision, only: collision_query_image_limit
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, delete_file_if_exists
  implicit none

  type(particles_soa) :: pcls
  type(mesh_type) :: mesh, failure_mesh
  type(sim_config) :: sim, failure_sim
  real(dp), allocatable :: v(:, :), x(:, :), w_photo(:)
  integer(i32), allocatable :: emit_elem(:)
  real(dp) :: gamma_in, gamma_cut, area, residual, expected_vn, jitter_dt, expected_w
  real(dp) :: ray_dir(3), tri_v0(3, 2), tri_v1(3, 2), tri_v2(3, 2)
  integer(i32) :: n_macro, n_emit, collision_status, failure_ray, failure_bounce
  integer :: i
  real(dp), parameter :: rare_tail_a(6) = [0.0_dp, 4.0_dp, 6.0_dp, 8.0_dp, 10.0_dp, 12.0_dp]
  real(dp), parameter :: rare_tail_reference(6) = [ &
                         3.9894228040143268e-1_dp, 7.145258432405667e-6_dp, 1.5635697959709664e-10_dp, &
                         7.550262411946499e-17_dp, 7.474560254589328e-25_dp, 1.4605201169845548e-34_dp &
                         ]
  character(len=64) :: run_mode
  character(len=*), parameter :: photo_failure_path = 'test_injection_sampling_photo_failure_tmp.log'

  call get_command_argument(1, run_mode)
  if (trim(run_mode) == '--photo-query-failure-probe') then
    call run_photo_query_failure_probe()
    error stop 'photo query failure probe unexpectedly completed'
  end if

  call test_init(21)

  call seed_rng()

  call test_begin('thermal_velocity_sampling')
  allocate (x(3, 16), v(3, 16))
  call sample_shifted_maxwell_velocities( &
    [10.0d0, -5.0d0, 2.0d0], 2.0d0, v(:, 1:8), thermal_speed=3.0d0, sigma_cutoff=0.5d0 &
    )
  call assert_true( &
    all(abs(v(:, 1:8) - spread([10.0d0, -5.0d0, 2.0d0], dim=2, ncopies=8)) <= 1.5d0), &
    'thermal_speed branch should honor sigma_cutoff' &
    )
  call sample_shifted_maxwell_velocities( &
    [10.0d0, -5.0d0, 2.0d0], 2.0d0, v(:, 9:16), temperature_k=500.0d0, thermal_speed=0.0d0 &
    )
  call assert_true( &
    all(v(:, 9:16) == spread([10.0d0, -5.0d0, 2.0d0], dim=2, ncopies=8)), &
    'explicit zero thermal_speed should override temperature_k' &
    )
  call test_end()

  call test_begin('rare_tail_inflow_flux')
  do i = 1, size(rare_tail_a)
    gamma_in = compute_inflow_flux_from_drifting_maxwellian( &
               1.0_dp, 1.0_dp, k_boltzmann, [-rare_tail_a(i), 0.0_dp, 0.0_dp], [1.0_dp, 0.0_dp, 0.0_dp] &
               )
    call assert_close_dp( &
      gamma_in, rare_tail_reference(i), max(1.0e-12_dp*rare_tail_reference(i), tiny(1.0_dp)), &
      'rare-tail inflow flux mismatch' &
      )
  end do
  call test_end()

  call test_begin('beam_particles')
  call init_random_beam_particles( &
    pcls, 4_i32, -1.0d0, 2.0d0, 100.0d0, &
    [-0.5d0, -0.5d0, -0.5d0], [0.5d0, 0.5d0, 0.5d0], [1.0d0, 0.0d0, 0.0d0], &
    temperature_k=500.0d0 &
    )
  call assert_equal_i32(pcls%n, 4_i32, 'init_random_beam_particles count mismatch')
  call assert_true(all(pcls%alive), 'init_random_beam_particles should initialize alive flags')
  call assert_true(all(pcls%x >= -0.5d0 .and. pcls%x <= 0.5d0), 'beam positions should stay within the source box')
  call assert_true(all(pcls%q == -1.0d0), 'beam particle charge mismatch')
  call assert_true(all(pcls%m == 2.0d0), 'beam particle mass mismatch')
  call assert_true(all(pcls%w == 100.0d0), 'beam particle weight mismatch')
  call test_end()

  call test_begin('rare_tail_reservoir_sampling')
  call sample_reservoir_face_particles( &
    [-1.0_dp, -1.0_dp, -1.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], 'x_low', &
    [-1.0_dp, -0.1_dp, -0.1_dp], [-1.0_dp, 0.1_dp, 0.1_dp], [-8.0_dp, 0.0_dp, 0.0_dp], &
    k_boltzmann, 1.0_dp, 0.1_dp, x(:, 1:16), v(:, 1:16) &
    )
  call assert_true(all(v(1, 1:16) > 0.0_dp), 'rare positive tail collapsed to zero normal velocity')
  call test_end()

  call test_begin('inflow_flux')
  gamma_in = compute_inflow_flux_from_drifting_maxwellian( &
             1.0d10, 400.0d0, 1.0d-26, [0.0d0, 0.0d0, 3.0d0], [0.0d0, 0.0d0, 1.0d0] &
             )
  call assert_true(gamma_in > 0.0d0, 'inflow flux should be positive for inward drift')
  gamma_cut = compute_inflow_flux_from_drifting_maxwellian( &
              1.0d10, 400.0d0, 1.0d-26, [0.0d0, 0.0d0, 3.0d0], [0.0d0, 0.0d0, 1.0d0], vmin_normal=2.0d3 &
              )
  call assert_true(gamma_cut > 0.0d0, 'cutoff inflow flux should remain positive')
  call assert_true(gamma_cut < gamma_in, 'cutoff inflow flux should be smaller than full inflow flux')
  call test_end()

  call test_begin('face_area')
  area = compute_face_area_from_bounds('x_low', [-1.0d0, -2.0d0, -3.0d0], [-1.0d0, 2.0d0, 3.0d0])
  call assert_close_dp(area, 24.0d0, 1.0d-12, 'face area (x_low) mismatch')
  area = compute_face_area_from_bounds('y_high', [-4.0d0, 1.0d0, -1.5d0], [4.0d0, 1.0d0, 2.5d0])
  call assert_close_dp(area, 32.0d0, 1.0d-12, 'face area (y_high) mismatch')
  call test_end()

  call test_begin('reservoir_face_basic')
  call sample_reservoir_face_particles( &
    [-1.0d0, -1.0d0, -1.0d0], [1.0d0, 1.0d0, 1.0d0], 'x_low', &
    [-1.0d0, -0.5d0, -0.25d0], [-1.0d0, 0.5d0, 0.25d0], [2.0d0, 0.0d0, 0.0d0], &
    1.0d0, 700.0d0, 0.2d0, x, v &
    )

  call assert_true(all(x(1, :) > -1.0d0 .and. x(1, :) < 1.0d0), 'reservoir launch should be strictly inside the box')
  call assert_true( &
    all(x(1, :) + 1.0_dp <= 8.0_dp*epsilon(1.0_dp)), &
    'reservoir launch should remain within a few ulps of the source face' &
    )
  call assert_true(all(x(2, :) >= -0.5d0 .and. x(2, :) <= 0.5d0), 'reservoir y positions should stay in the opening')
  call assert_true(all(x(3, :) >= -0.25d0 .and. x(3, :) <= 0.25d0), 'reservoir z positions should stay in the opening')
  call assert_true(all(v(1, :) >= 0.0d0), 'reservoir normal velocities should be inward')
  call test_end()

  call test_begin('reservoir_face_launch_has_no_untracked_jitter')
  jitter_dt = 1.0d-3
  call seed_rng([2718_i32])
  call sample_reservoir_face_particles( &
    [-1.0d0, -1.0d0, -1.0d0], [1.0d0, 1.0d0, 1.0d0], 'x_low', &
    [-1.0d0, -0.1d0, -0.1d0], [-1.0d0, 0.1d0, 0.1d0], [2.0d0, 0.0d0, 0.0d0], &
    1.0d0, 0.0d0, 50.0d0, x(:, 1:4), v(:, 1:4), position_jitter_dt=jitter_dt &
    )
  call seed_rng([2718_i32])
  call sample_reservoir_face_particles( &
    [-1.0d0, -1.0d0, -1.0d0], [1.0d0, 1.0d0, 1.0d0], 'x_low', &
    [-1.0d0, -0.1d0, -0.1d0], [-1.0d0, 0.1d0, 0.1d0], [2.0d0, 0.0d0, 0.0d0], &
    1.0d0, 0.0d0, 50.0d0, x(:, 5:8), v(:, 5:8), position_jitter_dt=2.0_dp*jitter_dt &
    )
  call assert_true(all(x(:, 1:4) == x(:, 5:8)), 'jitter compatibility argument should not change launch positions')
  call assert_true(all(v(:, 1:4) == v(:, 5:8)), 'jitter compatibility argument should not change launch velocities')
  call assert_true(all(x(1, 1:4) > -1.0_dp .and. x(1, 1:4) < 1.0_dp), 'jitter-compatible launch should stay inside')
  call assert_true( &
    all(x(1, 1:4) + 1.0_dp <= 8.0_dp*epsilon(1.0_dp)), &
    'jitter-compatible launch should remain within a few ulps of the source face' &
    )
  call assert_true(all(x(1, 1:4) == x(1, 1)), 'jitter-compatible launch should remain on one normal plane')
  call assert_true( &
    all(x(2:3, 1:4) >= -0.1_dp .and. x(2:3, 1:4) <= 0.1_dp), &
    'jitter-compatible launch should stay inside the tangential opening' &
    )
  call assert_true(all(v(1, 1:4) > 0.0_dp), 'jitter-compatible launch velocities should point inward')
  call test_end()

  call test_begin('reservoir_face_barrier_shift')
  call sample_reservoir_face_particles( &
    [-1.0d0, -1.0d0, -1.0d0], [1.0d0, 1.0d0, 1.0d0], 'x_low', &
    [-1.0d0, -0.1d0, -0.1d0], [-1.0d0, 0.1d0, 0.1d0], [3.0d0, 0.0d0, 0.0d0], &
    1.0d0, 0.0d0, 0.1d0, x(:, 1:4), v(:, 1:4), barrier_normal_energy=4.0d0 &
    )
  expected_vn = sqrt(5.0d0)
  do i = 1, 4
    call assert_close_dp(v(1, i), expected_vn, 1.0d-10, 'barrier-corrected normal speed mismatch')
  end do
  call test_end()

  call test_begin('reservoir_face_barrier_cutoff')
  call sample_reservoir_face_particles( &
    [-1.0d0, -1.0d0, -1.0d0], [1.0d0, 1.0d0, 1.0d0], 'x_low', &
    [-1.0d0, -0.1d0, -0.1d0], [-1.0d0, 0.1d0, 0.1d0], [3.0d0, 0.0d0, 0.0d0], &
    1.0d0, 0.0d0, 0.1d0, x(:, 1:4), v(:, 1:4), barrier_normal_energy=4.0d0, &
    apply_barrier_energy_shift=.false. &
    )
  do i = 1, 4
    call assert_close_dp(v(1, i), 3.0d0, 1.0d-12, 'barrier cutoff-only normal speed mismatch')
  end do
  call test_end()

  tri_v0(:, 1) = [0.0d0, 0.0d0, 0.05d0]
  tri_v1(:, 1) = [1.0d0, 0.0d0, 0.05d0]
  tri_v2(:, 1) = [0.0d0, 1.0d0, 0.05d0]
  tri_v0(:, 2) = [1.0d0, 1.0d0, 0.05d0]
  tri_v1(:, 2) = [0.0d0, 1.0d0, 0.05d0]
  tri_v2(:, 2) = [1.0d0, 0.0d0, 0.05d0]
  call init_mesh(mesh, tri_v0, tri_v1, tri_v2)

  sim = sim_config()
  sim%use_box = .true.
  sim%batch_duration = 0.5d0
  sim%raycast_max_bounce = 16_i32
  sim%box_min = [0.0d0, 0.0d0, 0.0d0]
  sim%box_max = [1.0d0, 1.0d0, 1.0d0]
  sim%bc_low = [bc_open, bc_open, bc_open]
  sim%bc_high = [bc_open, bc_open, bc_open]
  ray_dir = [0.0d0, 1.0d0, -0.2d0]
  ray_dir = ray_dir/sqrt(sum(ray_dir*ray_dir))

  allocate (w_photo(10))
  allocate (emit_elem(10))
  call test_begin('photo_raycast_open')
  call sample_photo_raycast_particles( &
    mesh, sim, 'z_high', [0.49d0, 0.49d0, 1.0d0], [0.51d0, 0.51d0, 1.0d0], ray_dir, &
    1.0d0, 0.0d0, 1.0d0, 2.0d0, -1.0d0, 10_i32, x(:, 1:10), v(:, 1:10), w_photo, n_emit &
    )
  call assert_equal_i32(n_emit, 0_i32, 'photo_raycast open boundary should not emit particles')
  call test_end()

  call test_begin('photo_raycast_reflect')
  sim%bc_low(2) = bc_reflect
  sim%bc_high(2) = bc_reflect
  call sample_photo_raycast_particles( &
    mesh, sim, 'z_high', [0.49d0, 0.49d0, 1.0d0], [0.51d0, 0.51d0, 1.0d0], ray_dir, &
    1.0d0, 0.0d0, 1.0d0, 2.0d0, -1.0d0, 10_i32, x(:, 1:10), v(:, 1:10), w_photo, n_emit, emit_elem &
    )
  call assert_equal_i32(n_emit, 10_i32, 'photo_raycast reflect boundary should emit all rays')
  expected_w = 2.0d0*(0.02d0*0.02d0*abs(ray_dir(3)))*sim%batch_duration/(1.0d0*10.0d0)
  call assert_true( &
    all(abs(w_photo(1:n_emit) - expected_w) <= 1.0d-14), &
    'photo_raycast weights should match the projected-current budget for every emitted ray' &
    )
  call assert_true(all(v(3, 1:n_emit) > 0.0d0), 'photo_raycast emitted normal speed should be outward')
  call assert_true(all(emit_elem(1:n_emit) >= 1_i32), 'photo_raycast emit_elem should be positive')
  call assert_true(all(emit_elem(1:n_emit) <= mesh%nelem), 'photo_raycast emit_elem should be in range')
  call test_end()

  call test_begin('photo_raycast_cutoff')
  call sample_photo_raycast_particles( &
    mesh, sim, 'z_high', [0.49d0, 0.49d0, 1.0d0], [0.51d0, 0.51d0, 1.0d0], ray_dir, &
    1.0d0, 0.0d0, 0.0d0, 2.0d0, -1.0d0, 10_i32, x(:, 1:10), v(:, 1:10), w_photo, n_emit, &
    vmin_normal=2.5d0 &
    )
  call assert_equal_i32(n_emit, 10_i32, 'photo_raycast cutoff run should emit all rays')
  call assert_true(all(v(3, 1:n_emit) >= 2.5d0), 'photo_raycast cutoff should raise normal speed floor')
  call test_end()

  call test_begin('photo_raycast_periodic')
  sim%bc_low(2) = bc_periodic
  sim%bc_high(2) = bc_periodic
  call sample_photo_raycast_particles( &
    mesh, sim, 'z_high', [0.49d0, 0.49d0, 1.0d0], [0.51d0, 0.51d0, 1.0d0], ray_dir, &
    1.0d0, 0.0d0, 1.0d0, 2.0d0, -1.0d0, 10_i32, x(:, 1:10), v(:, 1:10), w_photo, n_emit &
    )
  call assert_equal_i32(n_emit, 10_i32, 'photo_raycast periodic boundary should emit all rays')
  call test_end()

  call test_begin('photo_raycast_in_box')
  tri_v0(:, 1) = [0.20d0, 0.20d0, 0.80d0]
  tri_v1(:, 1) = [1.20d0, 0.20d0, 0.80d0]
  tri_v2(:, 1) = [0.20d0, 1.20d0, 0.80d0]
  tri_v0(:, 2) = [0.20d0, 0.20d0, 0.60d0]
  tri_v1(:, 2) = [0.90d0, 0.20d0, 0.60d0]
  tri_v2(:, 2) = [0.20d0, 0.90d0, 0.60d0]
  call init_mesh(mesh, tri_v0, tri_v1, tri_v2)
  sim%bc_low = [bc_open, bc_open, bc_open]
  sim%bc_high = [bc_open, bc_open, bc_open]
  ray_dir = [0.0d0, 0.0d0, -1.0d0]
  call sample_photo_raycast_particles( &
    mesh, sim, 'z_high', [0.29d0, 0.29d0, 1.0d0], [0.31d0, 0.31d0, 1.0d0], ray_dir, &
    1.0d0, 0.0d0, 1.0d0, 2.0d0, -1.0d0, 1_i32, x(:, 1:1), v(:, 1:1), w_photo(1:1), n_emit, emit_elem(1:1) &
    )
  call assert_equal_i32(n_emit, 1_i32, 'photo_raycast should emit from in-box element')
  call assert_equal_i32(emit_elem(1), 2_i32, 'photo_raycast should ignore out-of-box element')
  call test_end()

  call test_begin('photo_raycast_periodic2')
  tri_v0(:, 1) = [0.20d0, -0.05d0, 0.05d0]
  tri_v1(:, 1) = [0.40d0, -0.05d0, 0.05d0]
  tri_v2(:, 1) = [0.20d0, 0.15d0, 0.05d0]
  tri_v0(:, 2) = [0.70d0, 0.70d0, 0.60d0]
  tri_v1(:, 2) = [0.90d0, 0.70d0, 0.60d0]
  tri_v2(:, 2) = [0.70d0, 0.90d0, 0.60d0]
  call init_mesh(mesh, tri_v0, tri_v1, tri_v2)
  sim%field_bc_mode = 'periodic2'
  sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  sim%bc_high = [bc_periodic, bc_periodic, bc_open]
  call prepare_periodic2_collision_mesh(mesh, sim)
  ray_dir = [0.0d0, 0.0d0, -1.0d0]
  call sample_photo_raycast_particles( &
    mesh, sim, 'z_high', [0.245d0, 0.995d0, 1.0d0], [0.255d0, 1.005d0, 1.0d0], ray_dir, &
    1.0d0, 0.0d0, 1.0d0, 2.0d0, -1.0d0, 1_i32, x(:, 1:1), v(:, 1:1), w_photo(1:1), n_emit, emit_elem(1:1) &
    )
  call assert_equal_i32(n_emit, 1_i32, 'photo_raycast periodic2 should emit from wrapped hit')
  call assert_true(x(2, 1) > sim%box_min(2), 'photo_raycast periodic2 wrapped y should be strictly inside the primary cell')
  call assert_true(x(2, 1) < sim%box_max(2), 'photo_raycast periodic2 should use wrapped emission position')
  call test_end()

  call test_begin('photo_raycast_periodic2_normal_offset_wraps')
  tri_v0(:, 1) = [0.20d0, 0.90d0, 0.50d0]
  tri_v1(:, 1) = [0.80d0, 0.90d0, 0.50d0]
  tri_v2(:, 1) = [0.20d0, 1.10d0, 0.30d0]
  tri_v0(:, 2) = [0.70d0, 0.70d0, 0.60d0]
  tri_v1(:, 2) = [0.90d0, 0.70d0, 0.60d0]
  tri_v2(:, 2) = [0.70d0, 0.90d0, 0.60d0]
  call init_mesh(mesh, tri_v0, tri_v1, tri_v2)
  call prepare_periodic2_collision_mesh(mesh, sim)
  call sample_photo_raycast_particles( &
    mesh, sim, 'z_high', [0.30d0, 1.0d0 - 2.5d-13, 1.0d0], [0.40d0, 1.0d0 - 1.25d-13, 1.0d0], ray_dir, &
    1.0d0, 0.0d0, 1.0d0, 2.0d0, -1.0d0, 1_i32, x(:, 1:1), v(:, 1:1), w_photo(1:1), n_emit, emit_elem(1:1) &
    )
  call assert_equal_i32(n_emit, 1_i32, 'photo_raycast periodic2 seam case should emit one particle')
  call assert_true(x(2, 1) > sim%box_min(2), 'normal-offset emission should wrap above the periodic low face')
  call assert_true(x(2, 1) < sim%box_max(2), 'normal-offset emission should stay below the periodic high face')
  call test_end()

  call test_begin('photo_raycast_max_bounce')
  tri_v0(:, 1) = [0.0d0, 0.0d0, 0.05d0]
  tri_v1(:, 1) = [1.0d0, 0.0d0, 0.05d0]
  tri_v2(:, 1) = [0.0d0, 1.0d0, 0.05d0]
  tri_v0(:, 2) = [1.0d0, 1.0d0, 0.05d0]
  tri_v1(:, 2) = [0.0d0, 1.0d0, 0.05d0]
  tri_v2(:, 2) = [1.0d0, 0.0d0, 0.05d0]
  call init_mesh(mesh, tri_v0, tri_v1, tri_v2)
  sim%field_bc_mode = 'free'
  sim%bc_low(2) = bc_periodic
  sim%bc_high(2) = bc_periodic
  ray_dir = [0.0d0, 1.0d0, -0.2d0]
  ray_dir = ray_dir/sqrt(sum(ray_dir*ray_dir))

  sim%raycast_max_bounce = 2_i32
  call sample_photo_raycast_particles( &
    mesh, sim, 'z_high', [0.49d0, 0.49d0, 1.0d0], [0.51d0, 0.51d0, 1.0d0], ray_dir, &
    1.0d0, 0.0d0, 1.0d0, 2.0d0, -1.0d0, 10_i32, x(:, 1:10), v(:, 1:10), w_photo, n_emit &
    )
  call assert_equal_i32(n_emit, 0_i32, 'photo_raycast max bounce should terminate rays before emission')
  call test_end()

  call test_begin('photo_raycast_collision_status')
  call configure_photo_collision_failure_fixture(failure_mesh, failure_sim)
  call sample_photo_raycast_particles( &
    failure_mesh, failure_sim, 'z_high', [0.0d0, 0.0d0, 1.0d0], [1.0d0, 1.0d0, 1.0d0], &
    [0.0d0, 0.0d0, -1.0d0], 1.0d0, 0.0d0, 0.0d0, 1.0d0, -1.0d0, 4_i32, &
    x(:, 1:4), v(:, 1:4), w_photo(1:4), n_emit, &
    collision_failure_status=collision_status, collision_failure_ray=failure_ray, &
    collision_failure_bounce=failure_bounce &
    )
  call assert_equal_i32(collision_status, collision_query_image_limit, 'photo collision status mismatch')
  call assert_equal_i32(failure_ray, 1_i32, 'photo collision failure should select the lowest ray')
  call assert_equal_i32(failure_bounce, 0_i32, 'photo collision failure bounce mismatch')
  call assert_equal_i32(n_emit, 0_i32, 'incomplete photo collision query must not emit particles')
  call test_end()

  call test_begin('photo_raycast_collision_status_omitted')
  call test_photo_query_failure_without_status()
  call test_end()

  call test_begin('macro_particles_clamp')
  residual = -0.2d0
  call compute_macro_particles_for_batch( &
    0.0d0, 0.0d0, 1.0d0, [0.0d0, 0.0d0, 0.0d0], &
    [-1.0d0, -1.0d0, -1.0d0], [1.0d0, 1.0d0, 1.0d0], &
    'z_low', [-0.5d0, -0.5d0, -1.0d0], [0.5d0, 0.5d0, -1.0d0], &
    1.0d0, 10.0d0, residual, n_macro &
    )
  call assert_equal_i32(n_macro, 0_i32, 'macro count should clamp to zero for negative budget')
  call assert_close_dp(residual, 0.0d0, 1.0d-12, 'macro residual clamp mismatch')
  call test_end()

  call test_summary()

contains

  subroutine configure_photo_collision_failure_fixture(failure_mesh, failure_sim)
    type(mesh_type), intent(out) :: failure_mesh
    type(sim_config), intent(out) :: failure_sim
    real(dp) :: failure_v0(3, 1), failure_v1(3, 1), failure_v2(3, 1)

    failure_v0(:, 1) = [-2048.0d0, 0.0d0, 0.5d0]
    failure_v1(:, 1) = [2049.0d0, 0.0d0, 0.5d0]
    failure_v2(:, 1) = [0.5d0, 1.0d0, 0.5d0]
    call init_mesh(failure_mesh, failure_v0, failure_v1, failure_v2)

    failure_sim = sim_config()
    failure_sim%use_box = .true.
    failure_sim%batch_duration = 0.5d0
    failure_sim%raycast_max_bounce = 4_i32
    failure_sim%field_bc_mode = 'periodic2'
    failure_sim%box_min = [0.0d0, 0.0d0, 0.0d0]
    failure_sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    failure_sim%bc_low = [bc_periodic, bc_periodic, bc_open]
    failure_sim%bc_high = [bc_periodic, bc_periodic, bc_open]
    call prepare_periodic2_collision_mesh(failure_mesh, failure_sim)
  end subroutine configure_photo_collision_failure_fixture

  subroutine run_photo_query_failure_probe()
    type(mesh_type) :: probe_mesh
    type(sim_config) :: probe_sim
    real(dp) :: probe_x(3, 4), probe_v(3, 4), probe_w(4)
    integer(i32) :: probe_n_emit

    call configure_photo_collision_failure_fixture(probe_mesh, probe_sim)
    call sample_photo_raycast_particles( &
      probe_mesh, probe_sim, 'z_high', [0.0d0, 0.0d0, 1.0d0], [1.0d0, 1.0d0, 1.0d0], &
      [0.0d0, 0.0d0, -1.0d0], 1.0d0, 0.0d0, 0.0d0, 1.0d0, -1.0d0, 4_i32, &
      probe_x, probe_v, probe_w, probe_n_emit &
      )
  end subroutine run_photo_query_failure_probe

  subroutine test_photo_query_failure_without_status()
    character(len=1024) :: executable_path, command, child_line
    integer :: child_exit_status, child_cmd_status, child_unit, child_ios
    logical :: saw_context, saw_ray, saw_bounce, saw_status

    call get_command_argument(0, executable_path)
    call delete_file_if_exists(photo_failure_path)
    command = 'OMP_NUM_THREADS=2 "'//trim(executable_path)//'" --photo-query-failure-probe > '// &
              photo_failure_path//' 2>&1'
    call execute_command_line(trim(command), wait=.true., exitstat=child_exit_status, cmdstat=child_cmd_status)

    call assert_equal_i32(int(child_cmd_status, i32), 0_i32, 'photo failure probe command status mismatch')
    call assert_true(child_exit_status /= 0, 'photo failure probe should terminate with a nonzero status')

    saw_context = .false.
    saw_ray = .false.
    saw_bounce = .false.
    saw_status = .false.
    open (newunit=child_unit, file=photo_failure_path, status='old', action='read', iostat=child_ios)
    if (child_ios /= 0) error stop 'failed to read photo failure probe output'
    do
      read (child_unit, '(A)', iostat=child_ios) child_line
      if (child_ios /= 0) exit
      saw_context = saw_context .or. index(child_line, 'photo_raycast collision query incomplete') > 0
      saw_ray = saw_ray .or. index(child_line, 'ray=1') > 0
      saw_bounce = saw_bounce .or. index(child_line, 'bounce=0') > 0
      saw_status = saw_status .or. index(child_line, 'status=image_limit') > 0
    end do
    close (child_unit)
    call delete_file_if_exists(photo_failure_path)

    call assert_true(saw_context, 'photo failure message should identify the photo raycast context')
    call assert_true(saw_ray, 'photo failure message should include the lowest failing ray')
    call assert_true(saw_bounce, 'photo failure message should include the failing bounce')
    call assert_true(saw_status, 'photo failure message should include the collision status')
  end subroutine test_photo_query_failure_without_status
end program test_injection_sampling

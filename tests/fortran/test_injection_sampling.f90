!> 粒子注入モジュールのサンプリング分岐を重点的に検証するテスト。
program test_injection_sampling
  use bem_kinds, only: dp, i32
  use bem_types, only: particles_soa, mesh_type, sim_config, bc_open, bc_reflect, bc_periodic
  use bem_mesh, only: init_mesh
  use bem_injection, only: &
    seed_rng, sample_shifted_maxwell_velocities, init_random_beam_particles, &
    compute_inflow_flux_from_drifting_maxwellian, compute_face_area_from_bounds, &
    sample_reservoir_face_particles, compute_macro_particles_for_batch, sample_photo_raycast_particles
  use test_support, only: assert_true, assert_equal_i32, assert_close_dp
  implicit none

  type(particles_soa) :: pcls
  type(mesh_type) :: mesh
  type(sim_config) :: sim
  real(dp), allocatable :: v(:, :), x(:, :), w_photo(:)
  integer(i32), allocatable :: emit_elem(:)
  real(dp) :: gamma_in, gamma_cut, area, residual, expected_vn, jitter_dt, expected_w
  real(dp) :: ray_dir(3), tri_v0(3, 2), tri_v1(3, 2), tri_v2(3, 2)
  integer(i32) :: n_macro, n_emit
  integer :: i

  call seed_rng()

  allocate (v(3, 8))
  call sample_shifted_maxwell_velocities( &
    [10.0d0, -5.0d0, 2.0d0], 2.0d0, v, thermal_speed=3.0d0 &
  )
  call assert_true(all(abs(v) < 1.0d3), 'thermal_speed branch produced invalid velocities')

  call init_random_beam_particles( &
    pcls, 4_i32, -1.0d0, 2.0d0, 100.0d0, &
    [-0.5d0, -0.5d0, -0.5d0], [0.5d0, 0.5d0, 0.5d0], [1.0d0, 0.0d0, 0.0d0], &
    temperature_k=500.0d0 &
  )
  call assert_equal_i32(pcls%n, 4_i32, 'init_random_beam_particles count mismatch')
  call assert_true(all(pcls%alive), 'init_random_beam_particles should initialize alive flags')

  gamma_in = compute_inflow_flux_from_drifting_maxwellian( &
    1.0d10, 400.0d0, 1.0d-26, [0.0d0, 0.0d0, 3.0d0], [0.0d0, 0.0d0, 1.0d0] &
  )
  call assert_true(gamma_in > 0.0d0, 'inflow flux should be positive for inward drift')
  gamma_cut = compute_inflow_flux_from_drifting_maxwellian( &
    1.0d10, 400.0d0, 1.0d-26, [0.0d0, 0.0d0, 3.0d0], [0.0d0, 0.0d0, 1.0d0], vmin_normal=2.0d3 &
  )
  call assert_true(gamma_cut > 0.0d0, 'cutoff inflow flux should remain positive')
  call assert_true(gamma_cut < gamma_in, 'cutoff inflow flux should be smaller than full inflow flux')

  area = compute_face_area_from_bounds('x_low', [-1.0d0, -2.0d0, -3.0d0], [-1.0d0, 2.0d0, 3.0d0])
  call assert_close_dp(area, 24.0d0, 1.0d-12, 'face area (x_low) mismatch')
  area = compute_face_area_from_bounds('y_high', [-4.0d0, 1.0d0, -1.5d0], [4.0d0, 1.0d0, 2.5d0])
  call assert_close_dp(area, 32.0d0, 1.0d-12, 'face area (y_high) mismatch')

  deallocate (v)
  allocate (x(3, 16), v(3, 16))
  call sample_reservoir_face_particles( &
    [-1.0d0, -1.0d0, -1.0d0], [1.0d0, 1.0d0, 1.0d0], 'x_low', &
    [-1.0d0, -0.5d0, -0.25d0], [-1.0d0, 0.5d0, 0.25d0], [2.0d0, 0.0d0, 0.0d0], &
    1.0d0, 700.0d0, 0.2d0, x, v &
  )

  do i = 1, size(x, 2)
    call assert_true(x(1, i) > -1.0d0, 'reservoir sampled x should move inward from boundary')
    call assert_true(x(1, i) < -9.99999d-1, 'reservoir sampled x should not advect by batch_duration')
    call assert_true(v(1, i) >= 0.0d0, 'reservoir normal velocity should be inward')
  end do

  jitter_dt = 1.0d-3
  call sample_reservoir_face_particles( &
    [-1.0d0, -1.0d0, -1.0d0], [1.0d0, 1.0d0, 1.0d0], 'x_low', &
    [-1.0d0, -0.1d0, -0.1d0], [-1.0d0, 0.1d0, 0.1d0], [2.0d0, 0.0d0, 0.0d0], &
    1.0d0, 0.0d0, 50.0d0, x(:, 1:4), v(:, 1:4), position_jitter_dt=jitter_dt &
  )
  do i = 1, 4
    call assert_true(x(1, i) > -1.0d0, 'reservoir jittered x should stay inside the domain')
    call assert_true(x(1, i) <= -1.0d0 + 2.0d0 * jitter_dt + 1.0d-12, 'reservoir jittered x exceeded dt bound')
  end do

  call sample_reservoir_face_particles( &
    [-1.0d0, -1.0d0, -1.0d0], [1.0d0, 1.0d0, 1.0d0], 'x_low', &
    [-1.0d0, -0.1d0, -0.1d0], [-1.0d0, 0.1d0, 0.1d0], [3.0d0, 0.0d0, 0.0d0], &
    1.0d0, 0.0d0, 0.1d0, x(:, 1:4), v(:, 1:4), barrier_normal_energy=4.0d0 &
  )
  expected_vn = sqrt(5.0d0)
  do i = 1, 4
    call assert_close_dp(v(1, i), expected_vn, 1.0d-10, 'barrier-corrected normal speed mismatch')
  end do

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
  ray_dir = ray_dir / sqrt(sum(ray_dir * ray_dir))

  allocate (w_photo(10))
  allocate (emit_elem(10))
  call sample_photo_raycast_particles( &
    mesh, sim, 'z_high', [0.49d0, 0.49d0, 1.0d0], [0.51d0, 0.51d0, 1.0d0], ray_dir, &
    1.0d0, 0.0d0, 1.0d0, 2.0d0, -1.0d0, 10_i32, x(:, 1:10), v(:, 1:10), w_photo, n_emit &
  )
  call assert_equal_i32(n_emit, 0_i32, 'photo_raycast open boundary should not emit particles')

  sim%bc_low(2) = bc_reflect
  sim%bc_high(2) = bc_reflect
  call sample_photo_raycast_particles( &
    mesh, sim, 'z_high', [0.49d0, 0.49d0, 1.0d0], [0.51d0, 0.51d0, 1.0d0], ray_dir, &
    1.0d0, 0.0d0, 1.0d0, 2.0d0, -1.0d0, 10_i32, x(:, 1:10), v(:, 1:10), w_photo, n_emit, emit_elem &
  )
  call assert_equal_i32(n_emit, 10_i32, 'photo_raycast reflect boundary should emit all rays')
  expected_w = 2.0d0 * (0.02d0 * 0.02d0 * abs(ray_dir(3))) * sim%batch_duration / (1.0d0 * 10.0d0)
  call assert_close_dp(w_photo(1), expected_w, 1.0d-14, 'photo_raycast w_hit mismatch')
  call assert_true(all(v(3, 1:n_emit) > 0.0d0), 'photo_raycast emitted normal speed should be outward')
  call assert_true(all(emit_elem(1:n_emit) >= 1_i32), 'photo_raycast emit_elem should be positive')
  call assert_true(all(emit_elem(1:n_emit) <= mesh%nelem), 'photo_raycast emit_elem should be in range')

  sim%bc_low(2) = bc_periodic
  sim%bc_high(2) = bc_periodic
  call sample_photo_raycast_particles( &
    mesh, sim, 'z_high', [0.49d0, 0.49d0, 1.0d0], [0.51d0, 0.51d0, 1.0d0], ray_dir, &
    1.0d0, 0.0d0, 1.0d0, 2.0d0, -1.0d0, 10_i32, x(:, 1:10), v(:, 1:10), w_photo, n_emit &
  )
  call assert_equal_i32(n_emit, 10_i32, 'photo_raycast periodic boundary should emit all rays')

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

  tri_v0(:, 1) = [0.0d0, 0.0d0, 0.05d0]
  tri_v1(:, 1) = [1.0d0, 0.0d0, 0.05d0]
  tri_v2(:, 1) = [0.0d0, 1.0d0, 0.05d0]
  tri_v0(:, 2) = [1.0d0, 1.0d0, 0.05d0]
  tri_v1(:, 2) = [0.0d0, 1.0d0, 0.05d0]
  tri_v2(:, 2) = [1.0d0, 0.0d0, 0.05d0]
  call init_mesh(mesh, tri_v0, tri_v1, tri_v2)
  sim%bc_low(2) = bc_periodic
  sim%bc_high(2) = bc_periodic
  ray_dir = [0.0d0, 1.0d0, -0.2d0]
  ray_dir = ray_dir / sqrt(sum(ray_dir * ray_dir))

  sim%raycast_max_bounce = 2_i32
  call sample_photo_raycast_particles( &
    mesh, sim, 'z_high', [0.49d0, 0.49d0, 1.0d0], [0.51d0, 0.51d0, 1.0d0], ray_dir, &
    1.0d0, 0.0d0, 1.0d0, 2.0d0, -1.0d0, 10_i32, x(:, 1:10), v(:, 1:10), w_photo, n_emit &
  )
  call assert_equal_i32(n_emit, 0_i32, 'photo_raycast max bounce should terminate rays before emission')

  residual = -0.2d0
  call compute_macro_particles_for_batch( &
    0.0d0, 0.0d0, 1.0d0, [0.0d0, 0.0d0, 0.0d0], &
    [-1.0d0, -1.0d0, -1.0d0], [1.0d0, 1.0d0, 1.0d0], &
    'z_low', [-0.5d0, -0.5d0, -1.0d0], [0.5d0, 0.5d0, -1.0d0], &
    1.0d0, 10.0d0, residual, n_macro &
  )
  call assert_equal_i32(n_macro, 0_i32, 'macro count should clamp to zero for negative budget')
  call assert_close_dp(residual, 0.0d0, 1.0d-12, 'macro residual clamp mismatch')
end program test_injection_sampling

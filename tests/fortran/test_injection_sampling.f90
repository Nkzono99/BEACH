!> 粒子注入モジュールのサンプリング分岐を重点的に検証するテスト。
!! reservoir 面と photo ray は batch 生成が使う counter-based の 1 粒子 API（`sample_reservoir_face_particle_counter`、
!! `sample_photo_ray_slot`、`photo_ray_macro_weight`）で検査する。
program test_injection_sampling
  use bem_kinds, only: dp, i32, i64
  use bem_constants, only: k_boltzmann
  use bem_types, only: mesh_type, sim_config, bc_open, bc_reflect, bc_periodic
  use bem_mesh, only: init_mesh, prepare_periodic2_collision_mesh
  use bem_injection, only: &
    seed_rng, sample_uniform_positions, sample_shifted_maxwell_velocities, &
    compute_inflow_flux_from_drifting_maxwellian, compute_face_area_from_bounds, &
    compute_macro_particles_for_batch, sample_reservoir_face_particle_counter, sample_photo_ray_slot, &
    photo_ray_macro_weight
  use bem_status_codes, only: collision_query_ok, collision_query_image_limit
  use bem_random_stream, only: random_uniform
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp
  implicit none

  integer(i32), parameter :: rng_seed = 4242_i32
  type(mesh_type) :: mesh, failure_mesh
  type(sim_config) :: sim, failure_sim
  real(dp), allocatable :: v(:, :), x(:, :), w_photo(:)
  integer(i32), allocatable :: emit_elem(:), photo_status(:), photo_bounce(:)
  real(dp) :: gamma_in, gamma_cut, area, residual, expected_vn, expected_w, w_hit
  real(dp) :: ray_dir(3), tri_v0(3, 2), tri_v1(3, 2), tri_v2(3, 2)
  real(dp) :: contiguous_sample(3, 5), strided_sample(6, 10), next_contiguous, next_strided
  integer(i32) :: n_macro, n_emit
  integer :: i
  real(dp), parameter :: rare_tail_a(6) = [0.0_dp, 4.0_dp, 6.0_dp, 8.0_dp, 10.0_dp, 12.0_dp]
  real(dp), parameter :: rare_tail_reference(6) = [ &
                         3.9894228040143268e-1_dp, 7.145258432405667e-6_dp, 1.5635697959709664e-10_dp, &
                         7.550262411946499e-17_dp, 7.474560254589328e-25_dp, 1.4605201169845548e-34_dp &
                         ]

  call test_init(19)

  call seed_rng()

  call test_begin('sampling_preserves_strided_output_and_rng_state')
  do i = 1, 2
    strided_sample = -999.0_dp
    call seed_rng([867_i32])
    if (i == 1) then
      call sample_uniform_positions([-2.0_dp, 1.0_dp, 4.0_dp], [3.0_dp, 2.0_dp, 4.0_dp], contiguous_sample)
    else
      call sample_shifted_maxwell_velocities([1.0_dp, 2.0_dp, 3.0_dp], 1.0_dp, contiguous_sample, &
                                             thermal_speed=2.0_dp, sigma_cutoff=0.5_dp)
    end if
    call random_uniform(next_contiguous)
    call seed_rng([867_i32])
    if (i == 1) then
      call sample_uniform_positions([-2.0_dp, 1.0_dp, 4.0_dp], [3.0_dp, 2.0_dp, 4.0_dp], strided_sample(1:6:2, 1:10:2))
    else
      call sample_shifted_maxwell_velocities([1.0_dp, 2.0_dp, 3.0_dp], 1.0_dp, strided_sample(1:6:2, 1:10:2), &
                                             thermal_speed=2.0_dp, sigma_cutoff=0.5_dp)
    end if
    call random_uniform(next_strided)
    call assert_true(all(strided_sample(1:6:2, 1:10:2) == contiguous_sample), 'strides changed sampled values')
    call assert_true(all(strided_sample(2:6:2, :) == -999.0_dp), 'sampling overwrote unselected rows')
    call assert_true(all(strided_sample(:, 2:10:2) == -999.0_dp), 'sampling overwrote unselected columns')
    call assert_close_dp(next_strided, next_contiguous, 0.0_dp, 'strides changed the following RNG draw')
  end do
  call test_end()

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

  call test_begin('rare_tail_reservoir_sampling')
  call sample_reservoir_slots( &
    'x_low', [-1.0_dp, -0.1_dp, -0.1_dp], [-1.0_dp, 0.1_dp, 0.1_dp], [-8.0_dp, 0.0_dp, 0.0_dp], &
    k_boltzmann, 1.0_dp, 0.0_dp, 0.0_dp, x, v &
    )
  call assert_true(all(v(1, :) > 0.0_dp), 'rare positive tail collapsed to zero normal velocity')
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

  call test_begin('boundary_reservoir_launch_on_face')
  call sample_reservoir_slots( &
    'x_low', [-1.0d0, -0.5d0, -0.25d0], [-1.0d0, 0.5d0, 0.25d0], [2.0d0, 0.0d0, 0.0d0], &
    1.0d0, 700.0d0, 0.0d0, 0.0d0, x, v &
    )
  call assert_true(all(x(1, :) > -1.0d0 .and. x(1, :) < 1.0d0), 'reservoir launch should be strictly inside the box')
  ! One representable step inward and no untracked launch jitter: every particle starts on the same normal plane.
  call assert_true(all(x(1, :) == nearest(-1.0_dp, 1.0_dp)), 'reservoir launch should start one ulp inside the face')
  call assert_true(all(x(2, :) >= -0.5d0 .and. x(2, :) <= 0.5d0), 'reservoir y positions should stay in the opening')
  call assert_true(all(x(3, :) >= -0.25d0 .and. x(3, :) <= 0.25d0), 'reservoir z positions should stay in the opening')
  call assert_true(all(v(1, :) > 0.0d0), 'reservoir normal velocities should be inward')
  call test_end()

  call test_begin('boundary_reservoir_barrier_shift')
  call sample_reservoir_slots( &
    'x_low', [-1.0d0, -0.1d0, -0.1d0], [-1.0d0, 0.1d0, 0.1d0], [3.0d0, 0.0d0, 0.0d0], &
    1.0d0, 0.0d0, 4.0d0, 0.0d0, x(:, 1:4), v(:, 1:4) &
    )
  expected_vn = sqrt(5.0d0)
  do i = 1, 4
    call assert_close_dp(v(1, i), expected_vn, 1.0d-10, 'barrier-corrected normal speed mismatch')
  end do
  call test_end()

  call test_begin('boundary_reservoir_vmin_cutoff')
  call sample_reservoir_slots( &
    'x_low', [-1.0d0, -0.1d0, -0.1d0], [-1.0d0, 0.1d0, 0.1d0], [0.0d0, 0.0d0, 0.0d0], &
    k_boltzmann, 1.0d0, 0.0d0, 2.5d0, x, v &
    )
  call assert_true(all(v(1, :) >= 2.5d0), 'reservoir vmin_normal should raise the normal speed floor')
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

  allocate (w_photo(10), emit_elem(10), photo_status(10), photo_bounce(10))
  call test_begin('photo_raycast_open')
  call cast_photo_slots( &
    mesh, sim, [0.49d0, 0.49d0, 1.0d0], [0.51d0, 0.51d0, 1.0d0], ray_dir, 1.0d0, 1.0d0, 10_i32, n_emit, 0.0d0 &
    )
  call assert_true(all(photo_status == collision_query_ok), 'photo_raycast open boundary query status mismatch')
  call assert_equal_i32(n_emit, 0_i32, 'photo_raycast open boundary should not emit particles')
  call test_end()

  call test_begin('photo_raycast_reflect')
  sim%bc_low(2) = bc_reflect
  sim%bc_high(2) = bc_reflect
  w_hit = photo_ray_macro_weight( &
          'z_high', [0.49d0, 0.49d0, 1.0d0], [0.51d0, 0.51d0, 1.0d0], ray_dir, 2.0d0, -1.0d0, sim%batch_duration, 10_i32 &
          )
  expected_w = 2.0d0*(0.02d0*0.02d0*abs(ray_dir(3)))*sim%batch_duration/(1.0d0*10.0d0)
  call assert_close_dp(w_hit, expected_w, 1.0d-14, 'photo ray weight should match the projected-current budget')
  call cast_photo_slots( &
    mesh, sim, [0.49d0, 0.49d0, 1.0d0], [0.51d0, 0.51d0, 1.0d0], ray_dir, 1.0d0, w_hit, 10_i32, n_emit, 0.0d0 &
    )
  call assert_equal_i32(n_emit, 10_i32, 'photo_raycast reflect boundary should emit all rays')
  call assert_true(all(w_photo(1:n_emit) == w_hit), 'every emitted photo ray should carry the ray weight')
  call assert_true(all(v(3, 1:n_emit) > 0.0d0), 'photo_raycast emitted normal speed should be outward')
  call assert_true(all(emit_elem(1:n_emit) >= 1_i32), 'photo_raycast emit_elem should be positive')
  call assert_true(all(emit_elem(1:n_emit) <= mesh%nelem), 'photo_raycast emit_elem should be in range')
  call test_end()

  call test_begin('photo_raycast_cutoff')
  call cast_photo_slots( &
    mesh, sim, [0.49d0, 0.49d0, 1.0d0], [0.51d0, 0.51d0, 1.0d0], ray_dir, 0.0d0, 1.0d0, 10_i32, n_emit, 2.5d0 &
    )
  call assert_equal_i32(n_emit, 10_i32, 'photo_raycast cutoff run should emit all rays')
  call assert_true(all(v(3, 1:n_emit) >= 2.5d0), 'photo_raycast cutoff should raise normal speed floor')
  call test_end()

  call test_begin('photo_raycast_periodic')
  sim%bc_low(2) = bc_periodic
  sim%bc_high(2) = bc_periodic
  call cast_photo_slots( &
    mesh, sim, [0.49d0, 0.49d0, 1.0d0], [0.51d0, 0.51d0, 1.0d0], ray_dir, 1.0d0, 1.0d0, 10_i32, n_emit, 0.0d0 &
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
  call cast_photo_slots( &
    mesh, sim, [0.29d0, 0.29d0, 1.0d0], [0.31d0, 0.31d0, 1.0d0], ray_dir, 1.0d0, 1.0d0, 1_i32, n_emit, 0.0d0 &
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
  call cast_photo_slots( &
    mesh, sim, [0.245d0, 0.995d0, 1.0d0], [0.255d0, 1.005d0, 1.0d0], ray_dir, 1.0d0, 1.0d0, 1_i32, n_emit, 0.0d0 &
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
  call cast_photo_slots( &
    mesh, sim, [0.30d0, 1.0d0 - 2.5d-13, 1.0d0], [0.40d0, 1.0d0 - 1.25d-13, 1.0d0], ray_dir, 1.0d0, 1.0d0, 1_i32, &
    n_emit, 0.0d0 &
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
  call cast_photo_slots( &
    mesh, sim, [0.49d0, 0.49d0, 1.0d0], [0.51d0, 0.51d0, 1.0d0], ray_dir, 1.0d0, 1.0d0, 10_i32, n_emit, 0.0d0 &
    )
  call assert_true(all(photo_status == collision_query_ok), 'photo_raycast max bounce query status mismatch')
  call assert_equal_i32(n_emit, 0_i32, 'photo_raycast max bounce should terminate rays before emission')
  call test_end()

  call test_begin('photo_raycast_collision_status')
  call configure_photo_collision_failure_fixture(failure_mesh, failure_sim)
  call cast_photo_slots( &
    failure_mesh, failure_sim, [0.0d0, 0.0d0, 1.0d0], [1.0d0, 1.0d0, 1.0d0], [0.0d0, 0.0d0, -1.0d0], 0.0d0, 1.0d0, &
    4_i32, n_emit, 0.0d0 &
    )
  call assert_true(all(photo_status(1:4) == collision_query_image_limit), 'photo collision status mismatch')
  call assert_true(all(photo_bounce(1:4) == 0_i32), 'photo collision failure bounce mismatch')
  call assert_equal_i32(n_emit, 0_i32, 'incomplete photo collision query must not emit particles')
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

  !> box [-1, 1]^3 の面から、slot k = 1..size(x, 2) の counter で reservoir 粒子を生成する。
  subroutine sample_reservoir_slots(face, pos_low, pos_high, drift, m_particle, temperature_k, barrier, vmin, x, v)
    character(len=*), intent(in) :: face
    real(dp), intent(in) :: pos_low(3), pos_high(3), drift(3), m_particle, temperature_k, barrier, vmin
    real(dp), intent(out) :: x(:, :), v(:, :)
    integer :: k

    do k = 1, size(x, 2)
      call sample_reservoir_face_particle_counter( &
        [-1.0_dp, -1.0_dp, -1.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], face, pos_low, pos_high, drift, m_particle, &
        temperature_k, barrier, vmin, rng_seed, [1_i64, int(k, i64)], x(:, k), v(:, k) &
        )
    end do
  end subroutine sample_reservoir_slots

  !> z_high 面から ray slot k = 1..n_rays を発射し、当たった ray を x/v/w_photo/emit_elem の先頭へ詰める。
  !! 放出は m = 1、T = 0。各 ray の照会 status と bounce は photo_status/photo_bounce に残す。
  subroutine cast_photo_slots(mesh, sim, pos_low, pos_high, ray_direction, normal_drift, w_ray, n_rays, n_emit, vmin)
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    real(dp), intent(in) :: pos_low(3), pos_high(3), ray_direction(3), normal_drift, w_ray, vmin
    integer(i32), intent(in) :: n_rays
    integer(i32), intent(out) :: n_emit
    real(dp) :: x_ray(3), v_ray(3), w_ray_out
    integer(i32) :: k, elem
    logical :: emitted

    n_emit = 0_i32
    photo_status = collision_query_ok
    photo_bounce = 0_i32
    do k = 1_i32, n_rays
      call sample_photo_ray_slot( &
        mesh, sim, 'z_high', pos_low, pos_high, ray_direction, 1.0_dp, 0.0_dp, normal_drift, w_ray, rng_seed, &
        [int(k, i64), 1_i64], [int(k, i64), 2_i64], x_ray, v_ray, w_ray_out, elem, emitted, photo_status(k), &
        photo_bounce(k), vmin_normal=vmin &
        )
      if (.not. emitted) cycle
      n_emit = n_emit + 1_i32
      x(:, n_emit) = x_ray
      v(:, n_emit) = v_ray
      w_photo(n_emit) = w_ray_out
      emit_elem(n_emit) = elem
    end do
  end subroutine cast_photo_slots

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
end program test_injection_sampling

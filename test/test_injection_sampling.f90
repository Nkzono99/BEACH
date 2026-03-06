!> 粒子注入モジュールのサンプリング分岐を重点的に検証するテスト。
program test_injection_sampling
  use bem_kinds, only: dp, i32
  use bem_types, only: particles_soa
  use bem_injection, only: &
    seed_rng, sample_shifted_maxwell_velocities, init_random_beam_particles, &
    compute_inflow_flux_from_drifting_maxwellian, compute_face_area_from_bounds, &
    sample_reservoir_face_particles, compute_macro_particles_for_batch
  use test_support, only: assert_true, assert_equal_i32, assert_close_dp
  implicit none

  type(particles_soa) :: pcls
  real(dp), allocatable :: v(:, :), x(:, :)
  real(dp) :: gamma_in, gamma_cut, area, residual, expected_vn
  integer(i32) :: n_macro
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
    call assert_true(v(1, i) >= 0.0d0, 'reservoir normal velocity should be inward')
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

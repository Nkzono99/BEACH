!> 一様外部電場と速度グリッド reservoir 注入の設定・サンプリングを検証する。
program test_external_field_velocity_grid
  use bem_kinds, only: dp, i32
  use bem_app_config, only: app_config, default_app_config, load_app_config
  use bem_injection, only: compute_macro_particles_from_flux, sample_reservoir_velocity_grid_particles, seed_rng
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, assert_allclose_1d, delete_file_if_exists
  implicit none

  type(app_config) :: cfg
  character(len=*), parameter :: cfg_path = 'test_external_field_velocity_grid_tmp.toml'
  character(len=*), parameter :: grid_path = 'test_external_field_velocity_grid_tmp.csv'
  character(len=*), parameter :: interp_grid_path = 'test_external_field_velocity_grid_interp_tmp.csv'
  integer(i32) :: n_macro
  real(dp) :: residual
  real(dp) :: x(3, 5), v(3, 5)
  real(dp) :: x_interp(3, 32), v_interp(3, 32)

  call test_init(5)

  call write_velocity_grid_fixture(grid_path)
  call write_velocity_grid_interp_fixture(interp_grid_path)
  call write_config_fixture(cfg_path, grid_path)

  call test_begin('config_resolves_external_e_and_grid_flux')
  call default_app_config(cfg)
  call load_app_config(cfg_path, cfg)
  call assert_allclose_1d(cfg%sim%e0, [0.0d0, 0.0d0, 5.0d0], 1.0d-12, 'e0 angle conversion mismatch')
  call assert_true(trim(cfg%particle_species(1)%velocity_distribution) == 'grid', 'species-1 distribution mismatch')
  call assert_true(trim(cfg%particle_species(1)%velocity_grid_sampling) == 'rectilinear', 'species-1 grid sampling mismatch')
  call assert_close_dp(cfg%particle_species(1)%particle_flux_m2_s, 2.0d0, 1.0d-12, 'species-1 flux mismatch')
  call assert_close_dp(cfg%particle_species(1)%w_particle, 1.0d0, 1.0d-12, 'species-1 target weight mismatch')
  call assert_close_dp(cfg%particle_species(2)%particle_flux_m2_s, 3.0d0, 1.0d-12, 'species-2 current flux mismatch')
  call assert_close_dp(cfg%particle_species(2)%w_particle, 1.0d0, 1.0d-12, 'species-2 target weight mismatch')
  call test_end()

  call test_begin('macro_count_from_particle_flux')
  residual = 0.0d0
  call compute_macro_particles_from_flux( &
    1.25d0, 'z_high', [0.0d0, 0.0d0, 1.0d0], [1.0d0, 1.0d0, 1.0d0], 1.0d0, 1.0d0, residual, n_macro &
    )
  call assert_equal_i32(n_macro, 1_i32, 'first flux macro count mismatch')
  call assert_close_dp(residual, 0.25d0, 1.0d-12, 'first flux residual mismatch')
  call compute_macro_particles_from_flux( &
    1.25d0, 'z_high', [0.0d0, 0.0d0, 1.0d0], [1.0d0, 1.0d0, 1.0d0], 1.0d0, 1.0d0, residual, n_macro &
    )
  call assert_equal_i32(n_macro, 1_i32, 'second flux macro count mismatch')
  call assert_close_dp(residual, 0.5d0, 1.0d-12, 'second flux residual mismatch')
  call test_end()

  call test_begin('velocity_grid_samples_inward_entries')
  call seed_rng([123_i32])
  call sample_reservoir_velocity_grid_particles( &
    [0.0d0, 0.0d0, 0.0d0], [1.0d0, 1.0d0, 1.0d0], 'z_high', &
    [0.0d0, 0.0d0, 1.0d0], [1.0d0, 1.0d0, 1.0d0], grid_path, 'phase_space', 1.0d0, x, v &
    )
  call assert_true(all(v(3, :) < 0.0d0), 'velocity grid should sample inward velocities only')
  call assert_true(all(v(3, :) >= -2.0d0), 'line grid interpolation lower bound mismatch')
  call assert_true(any(abs(v(3, :) + 2.0d0) > 1.0d-9), 'line grid should interpolate between rows')
  call test_end()

  call test_begin('velocity_grid_discrete_sampling_keeps_rows')
  call seed_rng([234_i32])
  call sample_reservoir_velocity_grid_particles( &
    [0.0d0, 0.0d0, 0.0d0], [1.0d0, 1.0d0, 1.0d0], 'z_high', &
    [0.0d0, 0.0d0, 1.0d0], [1.0d0, 1.0d0, 1.0d0], grid_path, 'phase_space', 1.0d0, x, v, &
    velocity_grid_sampling='discrete' &
    )
  call assert_true(all(abs(v(1, :)) <= 1.0d-12), 'discrete vx should keep CSV rows')
  call assert_true(all(abs(v(2, :)) <= 1.0d-12), 'discrete vy should keep CSV rows')
  call assert_true(all(abs(v(3, :) + 2.0d0) <= 1.0d-12), 'discrete vz should keep inward CSV row')
  call test_end()

  call test_begin('velocity_grid_trilinear_interpolation_samples_cell')
  call seed_rng([456_i32])
  call sample_reservoir_velocity_grid_particles( &
    [0.0d0, 0.0d0, 0.0d0], [1.0d0, 1.0d0, 1.0d0], 'z_high', &
    [0.0d0, 0.0d0, 1.0d0], [1.0d0, 1.0d0, 1.0d0], interp_grid_path, 'flux_weighted', 1.0d0, x_interp, v_interp, &
    velocity_grid_sampling='rectilinear' &
    )
  call assert_true(all(v_interp(1, :) >= -1.0d0 .and. v_interp(1, :) <= 1.0d0), 'interpolated vx out of bounds')
  call assert_true(all(v_interp(2, :) >= -2.0d0 .and. v_interp(2, :) <= 2.0d0), 'interpolated vy out of bounds')
  call assert_true(all(v_interp(3, :) >= -3.0d0 .and. v_interp(3, :) <= -1.0d0), 'interpolated vz out of bounds')
  call assert_true(any(abs(abs(v_interp(1, :)) - 1.0d0) > 1.0d-6), 'vx should contain interpolated interior values')
  call assert_true(any(abs(abs(v_interp(2, :)) - 2.0d0) > 1.0d-6), 'vy should contain interpolated interior values')
  call assert_true(any(abs(v_interp(3, :) + 3.0d0) > 1.0d-6 .and. abs(v_interp(3, :) + 1.0d0) > 1.0d-6), &
                   'vz should contain interpolated interior values')
  call test_end()

  call delete_file_if_exists(cfg_path)
  call delete_file_if_exists(grid_path)
  call delete_file_if_exists(interp_grid_path)
  call test_summary()

contains

  subroutine write_velocity_grid_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open velocity grid fixture'
    write (u, '(a)') 'vx_m_s,vy_m_s,vz_m_s,f'
    write (u, '(a)') '0.0,0.0,-2.0,1.0'
    write (u, '(a)') '0.0,0.0,1.0,100.0'
    close (u)
  end subroutine write_velocity_grid_fixture

  subroutine write_velocity_grid_interp_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios, ix, iy, iz
    real(dp), parameter :: vx_values(2) = [-1.0d0, 1.0d0]
    real(dp), parameter :: vy_values(2) = [-2.0d0, 2.0d0]
    real(dp), parameter :: vz_values(2) = [-3.0d0, -1.0d0]

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open interpolated velocity grid fixture'
    write (u, '(a)') 'vx_m_s,vy_m_s,vz_m_s,f'
    do iz = 1, size(vz_values)
      do iy = 1, size(vy_values)
        do ix = 1, size(vx_values)
          write (u, '(4(es24.16, :, ","))') vx_values(ix), vy_values(iy), vz_values(iz), 1.0d0
        end do
      end do
    end do
    close (u)
  end subroutine write_velocity_grid_interp_fixture

  subroutine write_config_fixture(path, grid_csv)
    character(len=*), intent(in) :: path, grid_csv
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'dt = 1.0'
    write (u, '(a)') 'batch_duration_step = 2.0'
    write (u, '(a)') 'batch_count = 2'
    write (u, '(a)') 'e0_abs = 5.0'
    write (u, '(a)') 'e0_phi_xy_deg = 0.0'
    write (u, '(a)') 'e0_phi_z_deg = 90.0'
    write (u, '(a)') 'use_box = true'
    write (u, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (u, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'source_mode = "reservoir_face"'
    write (u, '(a)') 'velocity_distribution = "grid"'
    write (u, '(a)') 'velocity_grid_pdf_kind = "phase_space"'
    write (u, '(a)') 'velocity_grid_sampling = "rectilinear"'
    write (u, '(a)') 'velocity_grid_path = "'//trim(grid_csv)//'"'
    write (u, '(a)') 'particle_flux_m2_s = 2.0'
    write (u, '(a)') 'q_particle = -1.0'
    write (u, '(a)') 'm_particle = 1.0'
    write (u, '(a)') 'target_macro_particles_per_batch = 4'
    write (u, '(a)') 'inject_face = "z_high"'
    write (u, '(a)') 'pos_low = [0.0, 0.0, 1.0]'
    write (u, '(a)') 'pos_high = [1.0, 1.0, 1.0]'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'source_mode = "reservoir_face"'
    write (u, '(a)') 'velocity_distribution = "grid"'
    write (u, '(a)') 'velocity_grid_pdf_kind = "flux_weighted"'
    write (u, '(a)') 'velocity_grid_path = "'//trim(grid_csv)//'"'
    write (u, '(a)') 'current_density_a_m2 = 6.0'
    write (u, '(a)') 'q_particle = -2.0'
    write (u, '(a)') 'm_particle = 1.0'
    write (u, '(a)') 'target_macro_particles_per_batch = 6'
    write (u, '(a)') 'inject_face = "z_high"'
    write (u, '(a)') 'pos_low = [0.0, 0.0, 1.0]'
    write (u, '(a)') 'pos_high = [1.0, 1.0, 1.0]'

    close (u)
  end subroutine write_config_fixture

end program test_external_field_velocity_grid

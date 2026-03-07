!> app_config の設定読込と派生値計算を検証するテスト。
program test_app_config_parser
  use bem_kinds, only: dp, i32
  use bem_types, only: bc_open, bc_reflect, bc_periodic
  use bem_app_config, only: app_config, default_app_config, load_app_config, &
    particles_per_batch_from_config, total_particles_from_config
  use test_support, only: assert_true, assert_equal_i32, assert_close_dp, assert_allclose_1d, delete_file_if_exists
  implicit none

  type(app_config) :: cfg, photo_cfg
  character(len=*), parameter :: cfg_path = 'test_app_config_parser_tmp.toml'
  character(len=*), parameter :: photo_cfg_path = 'test_app_config_parser_photo_tmp.toml'

  call write_config_fixture(cfg_path)
  call write_photo_config_fixture(photo_cfg_path)

  call default_app_config(cfg)
  call load_app_config(cfg_path, cfg)

  call assert_true(trim(cfg%mesh_mode) == 'template', 'mesh.mode was not parsed')
  call assert_true(trim(cfg%templates(2)%kind) == 'sphere', 'second template kind mismatch')
  call assert_equal_i32(cfg%n_particle_species, 2_i32, 'n_particle_species mismatch')
  call assert_equal_i32(particles_per_batch_from_config(cfg), 5_i32, 'per-batch particle count mismatch')
  call assert_equal_i32(total_particles_from_config(cfg), 15_i32, 'total particle count mismatch')
  call assert_equal_i32(cfg%n_particles, 15_i32, 'cached n_particles mismatch')
  call assert_equal_i32(cfg%sim%bc_low(1), bc_periodic, 'bc_x_low mismatch')
  call assert_equal_i32(cfg%sim%bc_high(2), bc_reflect, 'bc_y_high mismatch')
  call assert_equal_i32(cfg%sim%bc_low(3), bc_open, 'bc_z_low mismatch')
  call assert_true(trim(cfg%sim%reservoir_potential_model) == 'infinity_barrier', 'reservoir_potential_model mismatch')
  call assert_close_dp(cfg%sim%phi_infty, -2.0d0, 1.0d-12, 'phi_infty mismatch')
  call assert_equal_i32(cfg%sim%injection_face_phi_grid_n, 5_i32, 'injection_face_phi_grid_n mismatch')
  call assert_equal_i32(cfg%history_stride, 2_i32, 'history_stride mismatch')
  call assert_close_dp(cfg%sim%dt, 2.5d-9, 1.0d-15, 'dt mismatch')

  call default_app_config(photo_cfg)
  call load_app_config(photo_cfg_path, photo_cfg)

  call assert_equal_i32(photo_cfg%n_particle_species, 1_i32, 'photo n_particle_species mismatch')
  call assert_true(trim(photo_cfg%particle_species(1)%source_mode) == 'photo_raycast', 'photo source_mode mismatch')
  call assert_close_dp(photo_cfg%particle_species(1)%emit_current_density_a_m2, 2.0d-3, 1.0d-15, 'photo emit_current mismatch')
  call assert_equal_i32(photo_cfg%particle_species(1)%rays_per_batch, 40_i32, 'photo rays_per_batch mismatch')
  call assert_close_dp(photo_cfg%particle_species(1)%normal_drift_speed, 1.5d5, 1.0d-12, 'photo normal_drift_speed mismatch')
  call assert_true(photo_cfg%particle_species(1)%deposit_opposite_charge_on_emit, 'photo deposit_opposite_charge_on_emit mismatch')
  call assert_allclose_1d( &
    photo_cfg%particle_species(1)%ray_direction, [0.0d0, 0.0d0, -1.0d0], 1.0d-12, 'photo ray_direction mismatch' &
  )
  call assert_equal_i32(photo_cfg%sim%raycast_max_bounce, 16_i32, 'photo default raycast_max_bounce mismatch')
  call assert_equal_i32(particles_per_batch_from_config(photo_cfg), 0_i32, 'photo per-batch particle count mismatch')
  call assert_equal_i32(total_particles_from_config(photo_cfg), 0_i32, 'photo total particle count mismatch')
  call assert_equal_i32(photo_cfg%n_particles, 0_i32, 'photo cached n_particles mismatch')

  call delete_file_if_exists(cfg_path)
  call delete_file_if_exists(photo_cfg_path)

contains

  !> テスト専用の一時設定ファイルを書き出す。
  !! @param[in] path 書き出し先TOMLファイルパス。
  subroutine write_config_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'dt = 2.5e-9'
    write (u, '(a)') 'batch_count = 3 # comment should be ignored'
    write (u, '(a)') 'reservoir_potential_model = "infinity_barrier"'
    write (u, '(a)') 'phi_infty = -2.0'
    write (u, '(a)') 'injection_face_phi_grid_n = 5'
    write (u, '(a)') 'bc_x_low = "periodic"'
    write (u, '(a)') 'bc_y_high = "reflect"'
    write (u, '(a)') 'bc_z_low = "open"'
    write (u, '(a)') ''
    write (u, '(a)') '[particles]'
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'npcls_per_step = 4'
    write (u, '(a)') 'temperature_k = 10000.0'
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'npcls_per_step = 1'
    write (u, '(a)') 'drift_velocity = [1.0, 0.0, -2.0]'
    write (u, '(a)') ''
    write (u, '(a)') '[mesh]'
    write (u, '(a)') 'mode = "template"'
    write (u, '(a)') '[[mesh.templates]]'
    write (u, '(a)') 'kind = "plane"'
    write (u, '(a)') 'enabled = true'
    write (u, '(a)') '[[mesh.templates]]'
    write (u, '(a)') 'kind = "sphere"'
    write (u, '(a)') 'enabled = true'
    write (u, '(a)') 'radius = 0.25'
    write (u, '(a)') ''
    write (u, '(a)') '[output]'
    write (u, '(a)') 'history_stride = 2'

    close (u)
  end subroutine write_config_fixture

  !> テスト専用の photo_raycast 設定ファイルを書き出す。
  !! @param[in] path 書き出し先TOMLファイルパス。
  subroutine write_photo_config_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open photo config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'batch_count = 4'
    write (u, '(a)') 'batch_duration = 1.0e-7'
    write (u, '(a)') 'use_box = true'
    write (u, '(a)') 'box_min = [0.0, 0.0, 0.0]'
    write (u, '(a)') 'box_max = [1.0, 1.0, 1.0]'
    write (u, '(a)') ''
    write (u, '(a)') '[[particles.species]]'
    write (u, '(a)') 'source_mode = "photo_raycast"'
    write (u, '(a)') 'emit_current_density_a_m2 = 2.0e-3'
    write (u, '(a)') 'rays_per_batch = 40'
    write (u, '(a)') 'deposit_opposite_charge_on_emit = true'
    write (u, '(a)') 'normal_drift_speed = 1.5e5'
    write (u, '(a)') 'q_particle = -1.0'
    write (u, '(a)') 'm_particle = 1.0'
    write (u, '(a)') 'temperature_k = 0.0'
    write (u, '(a)') 'inject_face = "z_high"'
    write (u, '(a)') 'pos_low = [0.0, 0.0, 1.0]'
    write (u, '(a)') 'pos_high = [1.0, 1.0, 1.0]'

    close (u)
  end subroutine write_photo_config_fixture

end program test_app_config_parser

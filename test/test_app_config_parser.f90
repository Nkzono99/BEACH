!> app_config の設定読込と派生値計算を検証するテスト。
program test_app_config_parser
  use bem_kinds, only: dp, i32
  use bem_types, only: bc_open, bc_reflect, bc_periodic
  use bem_app_config, only: app_config, default_app_config, load_app_config, &
    particles_per_batch_from_config, total_particles_from_config
  use test_support, only: assert_true, assert_equal_i32, assert_close_dp, delete_file_if_exists
  implicit none

  type(app_config) :: cfg
  character(len=*), parameter :: cfg_path = 'test_app_config_parser_tmp.toml'

  call write_config_fixture(cfg_path)

  call default_app_config(cfg)
  call load_app_config(cfg_path, cfg)

  call assert_true(trim(cfg%mesh_mode) == 'template', 'mesh.mode was not parsed')
  call assert_equal_i32(cfg%n_templates, 2_i32, 'n_templates mismatch')
  call assert_true(trim(cfg%templates(2)%kind) == 'sphere', 'second template kind mismatch')
  call assert_equal_i32(cfg%n_particle_species, 2_i32, 'n_particle_species mismatch')
  call assert_equal_i32(particles_per_batch_from_config(cfg), 5_i32, 'per-batch particle count mismatch')
  call assert_equal_i32(total_particles_from_config(cfg), 15_i32, 'total particle count mismatch')
  call assert_equal_i32(cfg%n_particles, 15_i32, 'cached n_particles mismatch')
  call assert_equal_i32(cfg%sim%bc_low(1), bc_periodic, 'bc_x_low mismatch')
  call assert_equal_i32(cfg%sim%bc_high(2), bc_reflect, 'bc_y_high mismatch')
  call assert_equal_i32(cfg%sim%bc_low(3), bc_open, 'bc_z_low mismatch')
  call assert_equal_i32(cfg%history_stride, 2_i32, 'history_stride mismatch')
  call assert_close_dp(cfg%sim%dt, 2.5d-9, 1.0d-15, 'dt mismatch')

  call delete_file_if_exists(cfg_path)

contains

  !> テスト専用の一時設定ファイルを書き出す。
  subroutine write_config_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open config fixture'

    write (u, '(a)') '[sim]'
    write (u, '(a)') 'dt = 2.5e-9'
    write (u, '(a)') 'batch_count = 3 # comment should be ignored'
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
    write (u, '(a)') 'n_templates = 2'
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

end program test_app_config_parser

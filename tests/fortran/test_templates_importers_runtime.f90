!> テンプレート生成・OBJ取込・app_config実行時変換の連携を検証するテスト。
program test_templates_importers_runtime
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, particles_soa, injection_state
  use bem_templates, only: make_plane, make_box, make_cylinder, make_sphere
  use bem_mesh, only: init_mesh
  use bem_importers, only: load_obj_mesh
  use bem_app_config, only: &
    app_config, default_app_config, species_from_defaults, &
    build_mesh_from_config, init_particles_from_config, seed_particles_from_config, init_particle_batch_from_config
  use test_support, only: &
    assert_true, assert_equal_i32, assert_close_dp, assert_allclose_1d, delete_file_if_exists
  implicit none

  type(mesh_type) :: mesh
  type(app_config) :: cfg
  type(particles_soa) :: pcls
  type(injection_state) :: state
  real(dp), allocatable :: photo_emission_dq(:)
  real(dp) :: photo_v0(3, 1), photo_v1(3, 1), photo_v2(3, 1), expected_photo_counter
  integer(i32) :: i

  character(len=*), parameter :: obj_path = 'test_templates_runtime_tmp.obj'
  character(len=*), parameter :: missing_obj_path = 'test_templates_runtime_missing.obj'

  call make_plane(mesh, nx=2_i32, ny=3_i32)
  call assert_equal_i32(mesh%nelem, 12_i32, 'plane element count mismatch')

  call make_box(mesh, nx=2_i32, ny=1_i32, nz=1_i32)
  call assert_equal_i32(mesh%nelem, 20_i32, 'box element count mismatch')

  call make_cylinder(mesh, n_theta=8_i32, n_z=2_i32, cap=.false.)
  call assert_equal_i32(mesh%nelem, 32_i32, 'cylinder(no cap) element count mismatch')

  call make_cylinder(mesh, n_theta=8_i32, n_z=1_i32, cap=.true.)
  call assert_equal_i32(mesh%nelem, 32_i32, 'cylinder(cap) element count mismatch')

  call make_sphere(mesh, n_lon=8_i32, n_lat=4_i32)
  call assert_equal_i32(mesh%nelem, 48_i32, 'sphere element count mismatch')

  call write_obj_fixture(obj_path)
  call load_obj_mesh(obj_path, mesh)
  call assert_equal_i32(mesh%nelem, 3_i32, 'OBJ triangulation count mismatch')

  call default_app_config(cfg)
  cfg%mesh_mode = 'obj'
  cfg%obj_path = obj_path
  call build_mesh_from_config(cfg, mesh)
  call assert_equal_i32(mesh%nelem, 3_i32, 'mesh_mode=obj should load OBJ mesh')

  call default_app_config(cfg)
  cfg%mesh_mode = 'auto'
  cfg%obj_path = missing_obj_path
  cfg%n_templates = 1_i32
  cfg%templates(1)%enabled = .true.
  cfg%templates(1)%kind = 'plane'
  cfg%templates(1)%nx = 1_i32
  cfg%templates(1)%ny = 1_i32
  cfg%templates(1)%size_x = 1.0d0
  cfg%templates(1)%size_y = 1.0d0
  cfg%templates(1)%center = [0.0d0, 0.0d0, 0.1d0]
  call build_mesh_from_config(cfg, mesh)
  call assert_equal_i32(mesh%nelem, 2_i32, 'mesh_mode=auto should fallback to template mesh')

  call default_app_config(cfg)
  cfg%sim%rng_seed = 2468_i32
  cfg%sim%batch_count = 2_i32
  cfg%n_particle_species = 3_i32

  cfg%particle_species(1) = species_from_defaults()
  cfg%particle_species(1)%source_mode = 'volume_seed'
  cfg%particle_species(1)%npcls_per_step = 2_i32
  cfg%particle_species(1)%q_particle = 1.0d0
  cfg%particle_species(1)%m_particle = 1.0d0
  cfg%particle_species(1)%w_particle = 10.0d0
  cfg%particle_species(1)%pos_low = [0.0d0, 0.0d0, 0.0d0]
  cfg%particle_species(1)%pos_high = [0.0d0, 0.0d0, 0.0d0]
  cfg%particle_species(1)%drift_velocity = [1.0d0, 0.0d0, 0.0d0]
  cfg%particle_species(1)%temperature_k = 0.0d0

  cfg%particle_species(2) = species_from_defaults()
  cfg%particle_species(2)%source_mode = 'volume_seed'
  cfg%particle_species(2)%npcls_per_step = 1_i32
  cfg%particle_species(2)%q_particle = -2.0d0
  cfg%particle_species(2)%m_particle = 1.0d0
  cfg%particle_species(2)%w_particle = 20.0d0
  cfg%particle_species(2)%pos_low = [1.0d0, 1.0d0, 1.0d0]
  cfg%particle_species(2)%pos_high = [1.0d0, 1.0d0, 1.0d0]
  cfg%particle_species(2)%drift_velocity = [0.0d0, 1.0d0, 0.0d0]
  cfg%particle_species(2)%temperature_k = 0.0d0

  call seed_particles_from_config(cfg)
  call init_particles_from_config(cfg, pcls)
  call assert_equal_i32(pcls%n, 6_i32, 'init_particles_from_config count mismatch')
  call assert_close_dp(pcls%q(1), 1.0d0, 1.0d-12, 'species interleave q(1) mismatch')
  call assert_close_dp(pcls%q(2), -2.0d0, 1.0d-12, 'species interleave q(2) mismatch')
  call assert_close_dp(pcls%q(3), 1.0d0, 1.0d-12, 'species interleave q(3) mismatch')
  call assert_close_dp(pcls%q(4), -2.0d0, 1.0d-12, 'species interleave q(4) mismatch')
  call assert_close_dp(pcls%q(5), 1.0d0, 1.0d-12, 'species interleave q(5) mismatch')
  call assert_close_dp(pcls%q(6), 1.0d0, 1.0d-12, 'species interleave q(6) mismatch')
  call assert_allclose_1d(pcls%x(:, 2), [1.0d0, 1.0d0, 1.0d0], 1.0d-12, 'species-2 position mismatch')
  call assert_allclose_1d(pcls%v(:, 2), [0.0d0, 1.0d0, 0.0d0], 1.0d-12, 'species-2 velocity mismatch')
  call assert_allclose_1d(pcls%x(:, 5), [0.0d0, 0.0d0, 0.0d0], 1.0d-12, 'species-1 position mismatch')
  call assert_allclose_1d(pcls%v(:, 5), [1.0d0, 0.0d0, 0.0d0], 1.0d-12, 'species-1 velocity mismatch')
  call assert_true(all(pcls%alive), 'all particles should start alive')

  photo_v0(:, 1) = [0.0d0, 0.0d0, 0.1d0]
  photo_v1(:, 1) = [1.0d0, 0.0d0, 0.1d0]
  photo_v2(:, 1) = [0.0d0, 1.0d0, 0.1d0]
  call init_mesh(mesh, photo_v0, photo_v1, photo_v2)

  call default_app_config(cfg)
  cfg%sim%rng_seed = 999_i32
  cfg%sim%batch_count = 1_i32
  cfg%sim%batch_duration = 1.0d0
  cfg%sim%use_box = .true.
  cfg%sim%box_min = [0.0d0, 0.0d0, 0.0d0]
  cfg%sim%box_max = [1.0d0, 1.0d0, 1.0d0]
  cfg%sim%reservoir_potential_model = 'infinity_barrier'
  cfg%sim%phi_infty = 0.0d0
  cfg%sim%injection_face_phi_grid_n = 2_i32
  cfg%n_particle_species = 3_i32

  cfg%particle_species(1) = species_from_defaults()
  cfg%particle_species(1)%source_mode = 'volume_seed'
  cfg%particle_species(1)%npcls_per_step = 1_i32
  cfg%particle_species(1)%q_particle = 1.0d0
  cfg%particle_species(1)%m_particle = 1.0d0
  cfg%particle_species(1)%w_particle = 1.0d0
  cfg%particle_species(1)%pos_low = [0.0d0, 0.0d0, 0.5d0]
  cfg%particle_species(1)%pos_high = [0.0d0, 0.0d0, 0.5d0]
  cfg%particle_species(1)%drift_velocity = [0.0d0, 0.0d0, 0.0d0]
  cfg%particle_species(1)%temperature_k = 0.0d0

  cfg%particle_species(2) = species_from_defaults()
  cfg%particle_species(2)%source_mode = 'reservoir_face'
  cfg%particle_species(2)%number_density_cm3 = 0.0d0
  cfg%particle_species(2)%temperature_k = 0.0d0
  cfg%particle_species(2)%q_particle = 1.0d0
  cfg%particle_species(2)%m_particle = 1.0d0
  cfg%particle_species(2)%w_particle = 10.0d0
  cfg%particle_species(2)%inject_face = 'z_low'
  cfg%particle_species(2)%pos_low = [0.0d0, 0.0d0, 0.0d0]
  cfg%particle_species(2)%pos_high = [1.0d0, 1.0d0, 0.0d0]
  cfg%particle_species(2)%drift_velocity = [0.0d0, 0.0d0, 1.0d0]

  cfg%particle_species(3) = species_from_defaults()
  cfg%particle_species(3)%source_mode = 'photo_raycast'
  cfg%particle_species(3)%emit_current_density_a_m2 = 1.0d0
  cfg%particle_species(3)%rays_per_batch = 1_i32
  cfg%particle_species(3)%normal_drift_speed = 0.0d0
  cfg%particle_species(3)%ray_direction = [0.0d0, 0.0d0, -1.0d0]
  cfg%particle_species(3)%has_ray_direction = .true.
  cfg%particle_species(3)%q_particle = -1.0d0
  cfg%particle_species(3)%m_particle = 1.0d0
  cfg%particle_species(3)%temperature_k = 0.0d0
  cfg%particle_species(3)%deposit_opposite_charge_on_emit = .true.
  cfg%particle_species(3)%inject_face = 'z_high'
  cfg%particle_species(3)%pos_low = [0.0d0, 0.0d0, 1.0d0]
  cfg%particle_species(3)%pos_high = [0.1d0, 0.1d0, 1.0d0]

  allocate (state%macro_residual(3))
  state%macro_residual = 0.0d0
  allocate (photo_emission_dq(mesh%nelem))
  call seed_particles_from_config(cfg)
  call init_particle_batch_from_config(cfg, 1_i32, pcls, state, mesh=mesh, photo_emission_dq=photo_emission_dq)
  call assert_true(pcls%n >= 1_i32, 'batch particle count mismatch')
  call assert_close_dp(pcls%q(1), 1.0d0, 1.0d-12, 'mixed batch species-1 charge mismatch')
  expected_photo_counter = 0.0d0
  do i = 1_i32, pcls%n
    if (pcls%q(i) < 0.0d0) then
      expected_photo_counter = expected_photo_counter - pcls%q(i) * pcls%w(i)
      call assert_true(pcls%w(i) > 0.0d0, 'mixed batch photo weight should be positive')
    end if
  end do
  call assert_true(expected_photo_counter > 0.0d0, 'mixed batch should include emitted photo particles')
  call assert_close_dp(sum(photo_emission_dq), expected_photo_counter, 1.0d-12, 'photo counter charge mismatch')
  call assert_true(state%macro_residual(2) >= 0.0d0, 'reservoir residual should be non-negative')
  call assert_true(state%macro_residual(2) < 1.0d0, 'reservoir residual should be < 1')

  call delete_file_if_exists(obj_path)

contains

  !> テスト専用OBJを作成する（quad+triangle、負インデックス含む）。
  !! @param[in] path 書き出し先OBJファイルパス。
  subroutine write_obj_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open OBJ fixture'
    write (u, '(a)') 'v 0.0 0.0 0.0'
    write (u, '(a)') 'v 1.0 0.0 0.0'
    write (u, '(a)') 'v 1.0 1.0 0.0'
    write (u, '(a)') 'v 0.0 1.0 0.0'
    write (u, '(a)') 'f 1 2 3 4'
    write (u, '(a)') 'f -4 -3 -2'
    close (u)
  end subroutine write_obj_fixture

end program test_templates_importers_runtime

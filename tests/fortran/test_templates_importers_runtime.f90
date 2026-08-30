!> テンプレート生成・OBJ取込・app_config実行時変換の連携を検証するテスト。
program test_templates_importers_runtime
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, particles_soa, injection_state, surface_model_conductor, &
                       bc_open, bc_periodic
  use bem_templates, only: make_plane, make_plate_hole, make_disk, make_annulus, make_box, make_cylinder, make_sphere
  use bem_mesh, only: init_mesh
  use bem_importers, only: load_obj_mesh
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_app_config, only: &
    app_config, default_app_config, species_from_defaults, &
    build_mesh_from_config, seed_particles_from_config, init_particle_batch_from_config, &
    sample_species_state
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, assert_allclose_1d, delete_file_if_exists
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  type(mesh_type) :: mesh
  type(app_config) :: cfg
  type(particles_soa) :: pcls
  type(injection_state) :: state
  type(electrostatic_snapshot_type) :: snapshot
  real(dp), allocatable :: photo_emission_dq(:), photo_source_charge(:)
  real(dp) :: reservoir_x(3, 128), reservoir_v(3, 128)
  real(dp) :: photo_v0(3, 1), photo_v1(3, 1), photo_v2(3, 1)
  real(dp) :: expected_photo_counter, observed_photo_counter
  real(dp), allocatable :: outward_dot(:)
  integer(i32) :: i

  character(len=*), parameter :: obj_path = 'test_templates_runtime_tmp.obj'
  character(len=*), parameter :: missing_obj_path = 'test_templates_runtime_missing.obj'

  character(len=*), parameter :: crlf_obj_path = 'test_templates_runtime_crlf.obj'

  call test_init(11)

  call test_begin('template_shapes')
  call make_plane(mesh, nx=2_i32, ny=3_i32)
  call assert_template_mesh(mesh, 12_i32, 'plane')

  call make_plate_hole(mesh, size_x=2.0d0, size_y=1.0d0, radius=0.2d0, n_theta=16_i32, n_r=2_i32)
  call assert_template_mesh(mesh, 68_i32, 'plate_hole')

  call make_disk(mesh, n_theta=8_i32, n_r=2_i32)
  call assert_template_mesh(mesh, 24_i32, 'disk')

  call make_annulus(mesh, n_theta=8_i32, n_r=2_i32, radius=0.5d0, inner_radius=0.2d0)
  call assert_template_mesh(mesh, 32_i32, 'annulus')

  call make_box(mesh, nx=2_i32, ny=1_i32, nz=1_i32)
  call assert_template_mesh(mesh, 20_i32, 'box')

  call make_cylinder(mesh, n_theta=8_i32, n_z=2_i32, cap=.false.)
  call assert_template_mesh(mesh, 32_i32, 'cylinder(no cap)')

  call make_cylinder(mesh, n_theta=8_i32, n_z=1_i32, cap=.true.)
  call assert_template_mesh(mesh, 32_i32, 'cylinder(cap)')

  call make_cylinder(mesh, n_theta=8_i32, n_z=1_i32, cap_top=.true., cap_bottom=.false.)
  call assert_template_mesh(mesh, 24_i32, 'cylinder(cap_top only)')

  call make_cylinder(mesh, n_theta=8_i32, n_z=1_i32, cap_top=.false., cap_bottom=.true.)
  call assert_template_mesh(mesh, 24_i32, 'cylinder(cap_bottom only)')

  call make_sphere(mesh, n_lon=8_i32, n_lat=4_i32)
  call assert_template_mesh(mesh, 48_i32, 'sphere')
  call test_end()

  call test_begin('reservoir_runtime_ignores_legacy_jitter_dt')
  call default_app_config(cfg)
  cfg%sim%dt = 0.5_dp
  cfg%sim%use_box = .true.
  cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
  cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
  cfg%sim%bc_low = [bc_periodic, bc_periodic, bc_open]
  cfg%sim%bc_high = [bc_periodic, bc_periodic, bc_open]
  cfg%particle_species(1) = species_from_defaults()
  cfg%particle_species(1)%source_mode = 'reservoir_face'
  cfg%particle_species(1)%inject_face = 'z_high'
  cfg%particle_species(1)%pos_low = [0.99_dp, 0.5_dp, 1.0_dp]
  cfg%particle_species(1)%pos_high = cfg%particle_species(1)%pos_low
  cfg%particle_species(1)%drift_velocity = [4.0_dp, 0.0_dp, -1.0_dp]
  cfg%particle_species(1)%temperature_k = 0.0_dp
  cfg%particle_species(1)%m_particle = 1.0_dp
  call seed_particles_from_config(cfg)
  call sample_species_state( &
    cfg%sim, cfg%particle_species(1), 64_i32, reservoir_x(:, 1:64), reservoir_v(:, 1:64) &
    )
  cfg%sim%dt = 1.0_dp
  call seed_particles_from_config(cfg)
  call sample_species_state( &
    cfg%sim, cfg%particle_species(1), 64_i32, reservoir_x(:, 65:128), reservoir_v(:, 65:128) &
    )
  call assert_true( &
    all(reservoir_x(:, 1:64) == reservoir_x(:, 65:128)), &
    'legacy jitter dt must not change reservoir launch positions' &
    )
  call assert_true( &
    all(reservoir_v(:, 1:64) == reservoir_v(:, 65:128)), &
    'legacy jitter dt must not change reservoir launch velocities' &
    )
  call assert_true( &
    all(reservoir_x(1, :) == cfg%particle_species(1)%pos_low(1)), &
    'reservoir launch must retain the sampled aperture x coordinate' &
    )
  call assert_true( &
    all(reservoir_x(2, :) == cfg%particle_species(1)%pos_low(2)), &
    'reservoir launch must retain the sampled aperture y coordinate' &
    )
  call assert_true( &
    all(reservoir_x(3, :) >= cfg%sim%box_min(3) .and. reservoir_x(3, :) < cfg%sim%box_max(3)), &
    'reservoir launch must remain strictly inside the normal box extent' &
    )
  call assert_true( &
    all(cfg%sim%box_max(3) - reservoir_x(3, :) <= 8.0_dp*epsilon(1.0_dp)), &
    'reservoir launch must remain within a few ulps of the source face' &
    )
  call assert_true( &
    all(reservoir_x(3, :) == reservoir_x(3, 1)), &
    'reservoir launch must remain on one normal plane' &
    )
  call assert_true( &
    all(reservoir_v(3, :) < 0.0_dp), &
    'z-high reservoir launch velocities must point inward' &
    )
  call test_end()

  call test_begin('obj_import')
  call write_obj_fixture(obj_path)
  call load_obj_mesh(obj_path, mesh)
  call assert_equal_i32(mesh%nelem, 3_i32, 'OBJ triangulation count mismatch')
  call assert_allclose_1d(mesh%v0(:, 1), [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, 'OBJ fan triangle v0 mismatch')
  call assert_allclose_1d(mesh%v1(:, 1), [1.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, 'OBJ fan triangle v1 mismatch')
  call assert_allclose_1d(mesh%v2(:, 1), [1.0_dp, 1.0_dp, 0.0_dp], 0.0_dp, 'OBJ fan triangle v2 mismatch')
  call assert_allclose_1d(mesh%v1(:, 2), [1.0_dp, 1.0_dp, 0.0_dp], 0.0_dp, 'OBJ quad fan split mismatch')
  call assert_allclose_1d(mesh%v2(:, 2), [0.0_dp, 1.0_dp, 0.0_dp], 0.0_dp, 'OBJ quad fan endpoint mismatch')
  call assert_true( &
    all(mesh%v0(:, 3) == mesh%v0(:, 1)) .and. all(mesh%v1(:, 3) == mesh%v1(:, 1)) .and. &
    all(mesh%v2(:, 3) == mesh%v2(:, 1)), &
    'OBJ negative indices must resolve against the current vertex list' &
    )
  call assert_true(all(mesh%elem_mesh_id == 1_i32), 'OBJ importer must assign mesh_id=1')
  call test_end()

  call test_begin('obj_import_crlf')
  call write_obj_crlf_fixture(crlf_obj_path)
  call load_obj_mesh(crlf_obj_path, mesh)
  call assert_equal_i32(mesh%nelem, 2_i32, 'CRLF OBJ triangulation count mismatch')
  call assert_close_dp(sum(mesh%panel_area), 1.0_dp, 1.0e-14_dp, 'CRLF OBJ geometry mismatch')
  call delete_file_if_exists(crlf_obj_path)
  call test_end()

  call test_begin('obj_import_slash_formats')
  call write_obj_slash_fixture(obj_path)
  call load_obj_mesh(obj_path, mesh)
  call assert_equal_i32(mesh%nelem, 3_i32, 'OBJ slash-format triangle count mismatch')
  call assert_true( &
    all(abs(mesh%panel_area - 0.5_dp) < 1.0e-14_dp), &
    'v//vn, v/vt, and v/vt/vn faces must preserve vertex geometry' &
    )
  call assert_allclose_1d(mesh%centers(:, 1), [2.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.0_dp], 1.0e-14_dp, 'v//vn face')
  call assert_allclose_1d(mesh%centers(:, 2), [1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 0.0_dp], 1.0e-14_dp, 'v/vt face')
  call assert_allclose_1d(mesh%centers(:, 3), [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.0_dp], 1.0e-14_dp, 'v/vt/vn face')
  call test_end()

  call test_begin('obj_transform')
  call write_obj_fixture(obj_path)
  call default_app_config(cfg)
  cfg%mesh_mode = 'obj'
  cfg%obj_path = obj_path
  cfg%mesh_surface_side_policy = 'normal_plus'
  cfg%mesh_surface_model = 'conductor'
  cfg%obj_scale = 2.0d0
  cfg%obj_offset = [10.0d0, 0.0d0, 0.0d0]
  cfg%obj_rotation = [90.0d0, 0.0d0, 90.0d0]
  call build_mesh_from_config(cfg, mesh)
  call assert_equal_i32(mesh%nelem, 3_i32, 'transform should preserve element count')
  call assert_allclose_1d(mesh%v0(:, 1), [10.0_dp, 0.0_dp, 0.0_dp], 1.0e-10_dp, 'transformed v0 mismatch')
  call assert_allclose_1d(mesh%v1(:, 1), [10.0_dp, 2.0_dp, 0.0_dp], 1.0e-10_dp, 'transformed v1 mismatch')
  call assert_allclose_1d(mesh%v2(:, 1), [10.0_dp, 2.0_dp, 2.0_dp], 1.0e-10_dp, 'transformed v2 mismatch')
  call assert_true(all(mesh%elem_mesh_id == 1_i32), 'transform should preserve OBJ mesh_id')
  call assert_true(all(mesh%elem_surface_model == surface_model_conductor), 'transform should preserve surface_model')
  call assert_true(all(mesh%elem_vacuum_sign == 1_i32), 'OBJ normal_plus surface side mismatch')
  call test_end()

  call test_begin('mesh_mode_auto')
  call delete_file_if_exists(missing_obj_path)
  call default_app_config(cfg)
  cfg%mesh_mode = 'auto'
  cfg%obj_path = missing_obj_path
  cfg%templates(1)%enabled = .true.
  cfg%templates(1)%kind = 'plane'
  cfg%templates(1)%nx = 1_i32
  cfg%templates(1)%ny = 1_i32
  cfg%templates(1)%size_x = 1.0d0
  cfg%templates(1)%size_y = 1.0d0
  cfg%templates(1)%center = [0.0d0, 0.0d0, 0.1d0]
  call build_mesh_from_config(cfg, mesh)
  call assert_equal_i32(mesh%nelem, 2_i32, 'mesh_mode=auto should fallback to template mesh')
  call assert_true(all(mesh%elem_mesh_id == 1_i32), 'auto fallback template mesh_id mismatch')
  call assert_true( &
    all(abs(mesh%v0(3, :) - 0.1_dp) < 1.0e-15_dp) .and. &
    all(abs(mesh%v1(3, :) - 0.1_dp) < 1.0e-15_dp) .and. &
    all(abs(mesh%v2(3, :) - 0.1_dp) < 1.0e-15_dp), &
    'auto fallback must use the configured template center' &
    )
  call test_end()

  call test_begin('triangle_panel_surface_side_runtime')
  call default_app_config(cfg)
  cfg%mesh_mode = 'template'
  cfg%templates(1)%enabled = .true.
  cfg%templates(1)%kind = 'plane'
  cfg%templates(1)%surface_side_policy = 'normal_minus'
  call build_mesh_from_config(cfg, mesh)
  call assert_true(all(mesh%elem_vacuum_sign == -1_i32), 'template vacuum side mismatch')
  call assert_true(all(abs(mesh%vacuum_normals + mesh%normals) < 1.0e-15_dp), 'template vacuum normal mismatch')
  call test_end()

  call test_begin('triangle_panel_surface_side_preserves_prior_templates')
  call default_app_config(cfg)
  cfg%mesh_mode = 'template'
  cfg%n_templates = 2_i32
  cfg%templates(1)%enabled = .true.
  cfg%templates(1)%kind = 'plane'
  cfg%templates(1)%surface_side_policy = 'normal_plus'
  cfg%templates(2)%enabled = .true.
  cfg%templates(2)%kind = 'sphere'
  cfg%templates(2)%surface_side_policy = 'outward_closed'
  cfg%templates(2)%radius = 0.2_dp
  cfg%templates(2)%center = [0.0_dp, 0.0_dp, 1.0_dp]
  cfg%templates(2)%n_lon = 8_i32
  cfg%templates(2)%n_lat = 4_i32
  call build_mesh_from_config(cfg, mesh)
  call assert_equal_i32( &
    int(count(mesh%elem_mesh_id == 1_i32), i32), 2_i32, &
    'first template element count mismatch' &
    )
  call assert_equal_i32( &
    int(count(mesh%elem_mesh_id == 2_i32), i32), 48_i32, &
    'second template element count mismatch' &
    )
  call assert_true(all(mesh%elem_mesh_id(1:2) == 1_i32), 'template concatenation order mismatch')
  call assert_true( &
    all(mesh%elem_mesh_id /= 1_i32 .or. mesh%elem_vacuum_sign == 1_i32), &
    'later template must preserve the first template normal_plus side' &
    )
  call assert_true( &
    all(abs(mesh%vacuum_normals(:, 1:2) - mesh%normals(:, 1:2)) < 1.0e-15_dp), &
    'later template must preserve the first template vacuum normals' &
    )
  outward_dot = sum( &
                mesh%vacuum_normals*(mesh%centers - spread(cfg%templates(2)%center, 2, mesh%nelem)), dim=1 &
                )
  call assert_true( &
    all(mesh%elem_mesh_id /= 2_i32 .or. outward_dot > 0.0_dp), &
    'sphere template vacuum normals must point outward after concatenation' &
    )
  call test_end()

  call test_begin('mesh_mode_template')
  call default_app_config(cfg)
  cfg%mesh_mode = 'template'
  cfg%templates(1)%enabled = .true.
  cfg%templates(1)%kind = 'plate_hole'
  cfg%templates(1)%surface_model = 'conductor'
  cfg%templates(1)%size_x = 1.0d0
  cfg%templates(1)%size_y = 1.0d0
  cfg%templates(1)%radius = 0.2d0
  cfg%templates(1)%n_theta = 12_i32
  cfg%templates(1)%n_r = 1_i32
  call build_mesh_from_config(cfg, mesh)
  call assert_equal_i32(mesh%nelem, 28_i32, 'mesh_mode=template plate_hole element count mismatch')
  call assert_true(all(mesh%elem_surface_model == surface_model_conductor), 'template surface_model mismatch')

  call default_app_config(cfg)
  cfg%mesh_mode = 'template'
  cfg%templates(1)%enabled = .true.
  cfg%templates(1)%kind = 'annulus'
  cfg%templates(1)%radius = 0.5d0
  cfg%templates(1)%inner_radius = 0.1d0
  cfg%templates(1)%n_theta = 8_i32
  cfg%templates(1)%n_r = 1_i32
  call build_mesh_from_config(cfg, mesh)
  call assert_equal_i32(mesh%nelem, 16_i32, 'mesh_mode=template annulus element count mismatch')

  call default_app_config(cfg)
  cfg%mesh_mode = 'template'
  cfg%templates(1)%enabled = .true.
  cfg%templates(1)%kind = 'cylinder'
  cfg%templates(1)%n_theta = 8_i32
  cfg%templates(1)%n_z = 1_i32
  cfg%templates(1)%cap = .false.
  cfg%templates(1)%cap_top = .true.
  cfg%templates(1)%has_cap_top = .true.
  cfg%templates(1)%cap_bottom = .false.
  cfg%templates(1)%has_cap_bottom = .true.
  call build_mesh_from_config(cfg, mesh)
  call assert_equal_i32(mesh%nelem, 24_i32, 'mesh_mode=template cylinder cap_top/cap_bottom mismatch')
  call test_end()

  call test_begin('photo_batch')
  photo_v0(:, 1) = [0.0d0, 0.0d0, 0.1d0]
  photo_v1(:, 1) = [1.0d0, 0.0d0, 0.1d0]
  photo_v2(:, 1) = [0.0d0, 1.0d0, 0.1d0]
  call init_mesh(mesh, photo_v0, photo_v1, photo_v2)
  mesh%elem_vacuum_sign = 1_i32
  mesh%vacuum_normals = mesh%normals

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
  allocate (photo_emission_dq(mesh%nelem), photo_source_charge(mesh%nelem))
  call snapshot%init(mesh, cfg%sim, cfg%field, cfg%periodic2, cfg%panel)
  call snapshot%refresh(mesh)
  call seed_particles_from_config(cfg)
  call init_particle_batch_from_config( &
    cfg, 1_i32, pcls, state, mesh=mesh, snapshot=snapshot, photo_emission_dq=photo_emission_dq &
    )
  call assert_equal_i32(pcls%n, 2_i32, 'mixed batch particle count mismatch')
  call assert_equal_i32(int(count(pcls%species_id == 1_i32), i32), 1_i32, 'volume species count mismatch')
  call assert_equal_i32(int(count(pcls%species_id == 2_i32), i32), 0_i32, 'zero-density reservoir count mismatch')
  call assert_equal_i32(int(count(pcls%species_id == 3_i32), i32), 1_i32, 'photo species count mismatch')
  call assert_close_dp( &
    sum(pack(pcls%q, pcls%species_id == 1_i32)), 1.0d0, 1.0d-12, &
    'mixed batch species-1 charge mismatch' &
    )
  expected_photo_counter = cfg%particle_species(3)%emit_current_density_a_m2* &
                           (cfg%particle_species(3)%pos_high(1) - cfg%particle_species(3)%pos_low(1))* &
                           (cfg%particle_species(3)%pos_high(2) - cfg%particle_species(3)%pos_low(2))* &
                           cfg%sim%batch_duration
  observed_photo_counter = 0.0_dp
  photo_source_charge = 0.0_dp
  do i = 1_i32, pcls%n
    if (pcls%species_id(i) == 3_i32) then
      call assert_true(pcls%source_element(i) >= 1_i32 .and. pcls%source_element(i) <= mesh%nelem, &
                       'photo particle source provenance is outside the emitting mesh')
      observed_photo_counter = observed_photo_counter - pcls%q(i)*pcls%w(i)
      photo_source_charge(pcls%source_element(i)) = photo_source_charge(pcls%source_element(i)) - &
                                                    pcls%q(i)*pcls%w(i)
      call assert_true(pcls%w(i) > 0.0d0, 'mixed batch photo weight should be positive')
    else
      call assert_equal_i32(pcls%source_element(i), -1_i32, &
                            'non-photo batch particle unexpectedly carries surface provenance')
    end if
  end do
  call assert_close_dp(observed_photo_counter, expected_photo_counter, 1.0d-14, 'photo current budget mismatch')
  call assert_close_dp(sum(photo_emission_dq), expected_photo_counter, 1.0d-12, 'photo counter charge mismatch')
  call assert_allclose_1d(photo_source_charge, photo_emission_dq, 1.0d-12, &
                          'source provenance does not reproduce the emitted countercharge distribution')
  call assert_close_dp(state%macro_residual(2), 0.0_dp, 0.0_dp, 'zero-density reservoir residual mismatch')
  call test_end()

  call delete_file_if_exists(obj_path)
  call delete_file_if_exists(missing_obj_path)
  call delete_file_if_exists(crlf_obj_path)

  call test_summary()

contains

  subroutine assert_template_mesh(candidate, expected_nelem, label)
    type(mesh_type), intent(in) :: candidate
    integer(i32), intent(in) :: expected_nelem
    character(len=*), intent(in) :: label

    call assert_equal_i32(candidate%nelem, expected_nelem, trim(label)//' element count mismatch')
    call assert_true( &
      all(ieee_is_finite(candidate%v0)) .and. all(ieee_is_finite(candidate%v1)) .and. &
      all(ieee_is_finite(candidate%v2)), &
      trim(label)//' vertices must be finite' &
      )
    call assert_true( &
      all(ieee_is_finite(candidate%panel_area)) .and. all(candidate%panel_area > 0.0_dp), &
      trim(label)//' triangles must be finite and non-degenerate' &
      )
  end subroutine assert_template_mesh

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

  !> CRLF改行のOBJフィクスチャを作成する（クアッド1面）。
  subroutine write_obj_crlf_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios
    character(len=1), parameter :: CR = char(13)

    open (newunit=u, file=trim(path), status='replace', action='write', &
          form='unformatted', access='stream', iostat=ios)
    if (ios /= 0) error stop 'failed to open CRLF OBJ fixture'
    write (u) '# CRLF test'//CR//char(10)
    write (u) 'v 0.0 0.0 0.0'//CR//char(10)
    write (u) 'v 1.0 0.0 0.0'//CR//char(10)
    write (u) 'v 1.0 1.0 0.0'//CR//char(10)
    write (u) 'v 0.0 1.0 0.0'//CR//char(10)
    write (u) 'f 1 2 3 4'//CR//char(10)
    close (u)
  end subroutine write_obj_crlf_fixture

  !> slash付きvertex参照3形式を含むOBJフィクスチャを作成する。
  subroutine write_obj_slash_fixture(path)
    character(len=*), intent(in) :: path
    integer :: u, ios

    open (newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'failed to open v//vn OBJ fixture'
    write (u, '(a)') 'v 0.0 0.0 0.0'
    write (u, '(a)') 'v 1.0 0.0 0.0'
    write (u, '(a)') 'v 1.0 1.0 0.0'
    write (u, '(a)') 'v 0.0 1.0 0.0'
    write (u, '(a)') 'vt 0.0 0.0'
    write (u, '(a)') 'vt 1.0 0.0'
    write (u, '(a)') 'vt 1.0 1.0'
    write (u, '(a)') 'vt 0.0 1.0'
    write (u, '(a)') 'vn 0.0 0.0 1.0'
    write (u, '(a)') 'f 1//1 2//1 3//1'
    write (u, '(a)') 'f 1/1 3/3 4/4'
    write (u, '(a)') 'f 1/1/1 2/2/1 4/4/1'
    close (u)
  end subroutine write_obj_slash_fixture

end program test_templates_importers_runtime

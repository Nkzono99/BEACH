!> テンプレート生成・OBJ取込・app_config実行時変換の連携を検証するテスト。
program test_templates_importers_runtime
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type, surface_model_conductor, surface_model_insulator
  use bem_templates, only: make_plane, make_plate_hole, make_disk, make_annulus, make_box, make_cylinder, make_sphere
  use bem_importers, only: load_obj_mesh
  use bem_app_config, only: app_config, default_app_config, build_mesh_from_config
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_close_dp, assert_allclose_1d, delete_file_if_exists
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  type(mesh_type) :: mesh
  type(app_config) :: cfg
  real(dp), allocatable :: outward_dot(:)

  character(len=*), parameter :: obj_path = 'test_templates_runtime_tmp.obj'
  character(len=*), parameter :: missing_obj_path = 'test_templates_runtime_missing.obj'

  character(len=*), parameter :: crlf_obj_path = 'test_templates_runtime_crlf.obj'

  call test_init(9)

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

  cfg%mesh_parts(1)%kind = 'obj'
  cfg%mesh_parts(1)%path = obj_path
  cfg%mesh_parts(1)%surface_side_policy = 'normal_plus'
  cfg%mesh_parts(1)%surface_model = 'conductor'
  cfg%mesh_parts(1)%scale = 2.0d0
  cfg%mesh_parts(1)%offset = [10.0d0, 0.0d0, 0.0d0]
  cfg%mesh_parts(1)%rotation = [90.0d0, 0.0d0, 90.0d0]
  call build_mesh_from_config(cfg, mesh)
  call assert_equal_i32(mesh%nelem, 3_i32, 'transform should preserve element count')
  call assert_allclose_1d(mesh%v0(:, 1), [10.0_dp, 0.0_dp, 0.0_dp], 1.0e-10_dp, 'transformed v0 mismatch')
  call assert_allclose_1d(mesh%v1(:, 1), [10.0_dp, 2.0_dp, 0.0_dp], 1.0e-10_dp, 'transformed v1 mismatch')
  call assert_allclose_1d(mesh%v2(:, 1), [10.0_dp, 2.0_dp, 2.0_dp], 1.0e-10_dp, 'transformed v2 mismatch')
  call assert_true(all(mesh%elem_mesh_id == 1_i32), 'transform should preserve OBJ mesh_id')
  call assert_true(all(mesh%elem_surface_model == surface_model_conductor), 'transform should preserve surface_model')
  call assert_true(all(mesh%elem_vacuum_sign == 1_i32), 'OBJ normal_plus surface side mismatch')
  call test_end()

  call test_begin('mixed_parts_order_transforms_and_surfaces')
  call default_app_config(cfg)
  cfg%n_mesh_parts = 4_i32
  cfg%mesh_parts(1)%kind = 'plane'
  cfg%mesh_parts(1)%center = [0.0_dp, 0.0_dp, -1.0_dp]
  cfg%mesh_parts(1)%surface_side_policy = 'normal_minus'
  cfg%mesh_parts(2)%enabled = .false.
  cfg%mesh_parts(2)%kind = 'obj'
  cfg%mesh_parts(2)%path = missing_obj_path
  cfg%mesh_parts(3)%enabled = .true.
  cfg%mesh_parts(3)%kind = 'obj'
  cfg%mesh_parts(3)%path = obj_path
  cfg%mesh_parts(3)%surface_model = 'conductor'
  cfg%mesh_parts(3)%surface_side_policy = 'normal_plus'
  cfg%mesh_parts(3)%scale = 2.0_dp
  cfg%mesh_parts(3)%rotation = [90.0_dp, 0.0_dp, 90.0_dp]
  cfg%mesh_parts(3)%offset = [10.0_dp, 0.0_dp, 0.0_dp]
  cfg%mesh_parts(4)%enabled = .true.
  cfg%mesh_parts(4)%kind = 'plane'
  cfg%mesh_parts(4)%center = [0.0_dp, 0.0_dp, 3.0_dp]
  cfg%mesh_parts(4)%surface_side_policy = 'normal_plus'
  call build_mesh_from_config(cfg, mesh)
  call assert_equal_i32(mesh%nelem, 7_i32, 'two planes plus OBJ')
  call assert_true(all(mesh%elem_mesh_id == [1_i32, 1_i32, 2_i32, 2_i32, 2_i32, 3_i32, 3_i32]), &
                   'enabled parts consume mesh_id in list order')
  call assert_allclose_1d(mesh%v0(:, 1), [-0.5_dp, -0.5_dp, -1.0_dp], 1.0e-12_dp, 'first plane unchanged')
  call assert_allclose_1d(mesh%v2(:, 3), [10.0_dp, 2.0_dp, 2.0_dp], 1.0e-12_dp, 'OBJ transform only')
  call assert_allclose_1d(mesh%v0(:, 6), [-0.5_dp, -0.5_dp, 3.0_dp], 1.0e-12_dp, 'last plane unchanged')
  call assert_true(all(mesh%elem_surface_model(3:5) == surface_model_conductor), 'OBJ conductor group')
  call assert_true(all(mesh%elem_surface_model([1, 2, 6, 7]) == surface_model_insulator), 'template insulators')
  call assert_true(all(mesh%elem_vacuum_sign(1:2) == -1_i32), 'first part side preserved')
  call assert_true(all(mesh%elem_vacuum_sign(3:7) == 1_i32), 'later part sides')
  call assert_true(all(mesh%elem_epsilon_r == 1.0_dp), 'part epsilon defaults')
  call test_end()

  call test_begin('triangle_panel_surface_side_runtime')
  call default_app_config(cfg)

  cfg%mesh_parts(1)%enabled = .true.
  cfg%mesh_parts(1)%kind = 'plane'
  cfg%mesh_parts(1)%surface_side_policy = 'normal_minus'
  call build_mesh_from_config(cfg, mesh)
  call assert_true(all(mesh%elem_vacuum_sign == -1_i32), 'template vacuum side mismatch')
  call assert_true(all(abs(mesh%vacuum_normals + mesh%normals) < 1.0e-15_dp), 'template vacuum normal mismatch')
  call test_end()

  call test_begin('triangle_panel_surface_side_preserves_prior_templates')
  call default_app_config(cfg)

  cfg%n_mesh_parts = 2_i32
  cfg%mesh_parts(1)%enabled = .true.
  cfg%mesh_parts(1)%kind = 'plane'
  cfg%mesh_parts(1)%surface_side_policy = 'normal_plus'
  cfg%mesh_parts(2)%enabled = .true.
  cfg%mesh_parts(2)%kind = 'sphere'
  cfg%mesh_parts(2)%surface_side_policy = 'outward_closed'
  cfg%mesh_parts(2)%radius = 0.2_dp
  cfg%mesh_parts(2)%center = [0.0_dp, 0.0_dp, 1.0_dp]
  cfg%mesh_parts(2)%n_lon = 8_i32
  cfg%mesh_parts(2)%n_lat = 4_i32
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
                mesh%vacuum_normals*(mesh%centers - spread(cfg%mesh_parts(2)%center, 2, mesh%nelem)), dim=1 &
                )
  call assert_true( &
    all(mesh%elem_mesh_id /= 2_i32 .or. outward_dot > 0.0_dp), &
    'sphere template vacuum normals must point outward after concatenation' &
    )
  call test_end()

  call test_begin('template_parts')
  call default_app_config(cfg)

  cfg%mesh_parts(1)%enabled = .true.
  cfg%mesh_parts(1)%kind = 'plate_hole'
  cfg%mesh_parts(1)%surface_model = 'conductor'
  cfg%mesh_parts(1)%size_x = 1.0d0
  cfg%mesh_parts(1)%size_y = 1.0d0
  cfg%mesh_parts(1)%radius = 0.2d0
  cfg%mesh_parts(1)%n_theta = 12_i32
  cfg%mesh_parts(1)%n_r = 1_i32
  call build_mesh_from_config(cfg, mesh)
  call assert_equal_i32(mesh%nelem, 28_i32, 'template parts plate_hole element count mismatch')
  call assert_true(all(mesh%elem_surface_model == surface_model_conductor), 'template surface_model mismatch')

  call default_app_config(cfg)

  cfg%mesh_parts(1)%enabled = .true.
  cfg%mesh_parts(1)%kind = 'annulus'
  cfg%mesh_parts(1)%radius = 0.5d0
  cfg%mesh_parts(1)%inner_radius = 0.1d0
  cfg%mesh_parts(1)%n_theta = 8_i32
  cfg%mesh_parts(1)%n_r = 1_i32
  call build_mesh_from_config(cfg, mesh)
  call assert_equal_i32(mesh%nelem, 16_i32, 'template parts annulus element count mismatch')

  call default_app_config(cfg)

  cfg%mesh_parts(1)%enabled = .true.
  cfg%mesh_parts(1)%kind = 'cylinder'
  cfg%mesh_parts(1)%n_theta = 8_i32
  cfg%mesh_parts(1)%n_z = 1_i32
  cfg%mesh_parts(1)%cap = .false.
  cfg%mesh_parts(1)%cap_top = .true.
  cfg%mesh_parts(1)%has_cap_top = .true.
  cfg%mesh_parts(1)%cap_bottom = .false.
  cfg%mesh_parts(1)%has_cap_bottom = .true.
  call build_mesh_from_config(cfg, mesh)
  call assert_equal_i32(mesh%nelem, 24_i32, 'template parts cylinder cap_top/cap_bottom mismatch')
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

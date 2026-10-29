!> app_config から三角形メッシュを構築する実行時変換。
module bem_app_config_mesh_runtime
  use bem_kinds, only: dp, i32
  use bem_types, only: mesh_type
  use bem_types, only: surface_model_insulator, surface_model_conductor, surface_model_dielectric
  use bem_templates, only: make_plane, make_plate_hole, make_disk, make_annulus, make_box, make_cylinder, make_sphere, &
                           make_icosphere
  use bem_mesh, only: init_mesh
  use bem_panel_surface_sides, only: resolve_panel_surface_sides_by_mesh_id, panel_surface_side_ok
  use bem_importers, only: load_obj_mesh
  use bem_app_config_types, only: app_config, template_spec, mesh_part_config_path
  use bem_string_utils, only: lower_ascii
  implicit none

contains

  !> 有効な OBJ と template をリスト順に連結する。
  subroutine build_mesh_from_config(cfg, mesh)
    type(app_config), intent(in) :: cfg
    type(mesh_type), intent(out) :: mesh

    call build_parts_mesh(cfg, mesh)
    call apply_panel_surface_config(cfg, mesh)
  end subroutine build_mesh_from_config

  subroutine apply_panel_surface_config(cfg, mesh)
    type(app_config), intent(in) :: cfg
    type(mesh_type), intent(inout) :: mesh
    integer(i32) :: i, mesh_id, status, failed_mesh_id
    integer(i32), allocatable :: part_index(:)
    character(len=256) :: message

    if (any(mesh%panel_area <= 0.0_dp)) then
      error stop 'triangle_p0 requires finite, non-degenerate triangles.'
    end if
    ! Enabled parts own mesh_id 1, 2, ... in list order (build_parts_mesh).
    allocate (part_index(cfg%n_mesh_parts))
    mesh_id = 0_i32
    do i = 1, cfg%n_mesh_parts
      if (.not. cfg%mesh_parts(i)%enabled) cycle
      mesh_id = mesh_id + 1_i32
      part_index(mesh_id) = i
    end do
    call resolve_panel_surface_sides_by_mesh_id( &
      mesh, cfg%mesh_parts(part_index(1:mesh_id))%surface_side_policy, status, message, failed_mesh_id &
      )
    if (status /= panel_surface_side_ok) then
      error stop mesh_part_config_path(int(part_index(failed_mesh_id)))//'.surface_side: '//trim(message)
    end if
    if (any(abs(mesh%elem_vacuum_sign) /= 1_i32)) then
      error stop 'triangle_p0 requires a resolved vacuum side for every element.'
    end if
  end subroutine apply_panel_surface_config

  !> OBJ メッシュの全頂点にスケール→回転→平行移動を適用し再初期化する。
  !! 変換順序: v_new = R(rotation) * (v_old * scale) + offset
  !! 回転は度単位で x→y→z の順に外因性 (extrinsic) 回転を適用する。
  !! @param[inout] mesh 変換対象の三角形メッシュ。
  !! @param[in] scale 一様スケーリング係数。
  !! @param[in] rotation_deg 回転角度 [rx, ry, rz] (度)。
  !! @param[in] offset 平行移動ベクトル (3成分)。
  subroutine apply_obj_transform(mesh, scale, rotation_deg, offset)
    type(mesh_type), intent(inout) :: mesh
    real(dp), intent(in) :: scale
    real(dp), intent(in) :: rotation_deg(3)
    real(dp), intent(in) :: offset(3)
    real(dp), parameter :: deg2rad = acos(-1.0d0)/180.0d0
    real(dp) :: rx, ry, rz, cx, sx, cy, sy, cz, sz
    real(dp) :: R(3, 3), v(3)
    real(dp), allocatable :: tv0(:, :), tv1(:, :), tv2(:, :)
    real(dp), allocatable :: elem_epsilon_r(:)
    integer(i32), allocatable :: elem_mesh_id(:), elem_surface_model(:)
    integer(i32) :: i, n

    rx = rotation_deg(1)*deg2rad
    ry = rotation_deg(2)*deg2rad
    rz = rotation_deg(3)*deg2rad
    cx = cos(rx); sx = sin(rx)
    cy = cos(ry); sy = sin(ry)
    cz = cos(rz); sz = sin(rz)

    ! R = Rz * Ry * Rx (extrinsic x→y→z)
    R(1, 1) = cy*cz; R(1, 2) = sx*sy*cz - cx*sz; R(1, 3) = cx*sy*cz + sx*sz
    R(2, 1) = cy*sz; R(2, 2) = sx*sy*sz + cx*cz; R(2, 3) = cx*sy*sz - sx*cz
    R(3, 1) = -sy; R(3, 2) = sx*cy; R(3, 3) = cx*cy

    n = mesh%nelem
    allocate (tv0(3, n), tv1(3, n), tv2(3, n))
    if (allocated(mesh%elem_mesh_id)) then
      allocate (elem_mesh_id(n))
      elem_mesh_id = mesh%elem_mesh_id
    end if
    if (allocated(mesh%elem_surface_model)) then
      allocate (elem_surface_model(n))
      elem_surface_model = mesh%elem_surface_model
    end if
    if (allocated(mesh%elem_epsilon_r)) then
      allocate (elem_epsilon_r(n))
      elem_epsilon_r = mesh%elem_epsilon_r
    end if
    do i = 1, n
      v = mesh%v0(:, i)*scale
      tv0(:, i) = matmul(R, v) + offset
      v = mesh%v1(:, i)*scale
      tv1(:, i) = matmul(R, v) + offset
      v = mesh%v2(:, i)*scale
      tv2(:, i) = matmul(R, v) + offset
    end do
    if (allocated(elem_mesh_id) .and. allocated(elem_surface_model) .and. allocated(elem_epsilon_r)) then
      call init_mesh( &
        mesh, tv0, tv1, tv2, elem_mesh_id0=elem_mesh_id, elem_surface_model0=elem_surface_model, &
        elem_epsilon_r0=elem_epsilon_r &
        )
    else if (allocated(elem_mesh_id)) then
      call init_mesh(mesh, tv0, tv1, tv2, elem_mesh_id0=elem_mesh_id)
    else if (allocated(elem_surface_model)) then
      call init_mesh(mesh, tv0, tv1, tv2, elem_surface_model0=elem_surface_model)
    else
      call init_mesh(mesh, tv0, tv1, tv2)
    end if
  end subroutine apply_obj_transform

  !> 有効な OBJ と template を連結し、1つのメッシュへまとめる。
  !! @param[in] cfg parts 設定を含むアプリ設定。
  !! @param[out] mesh 連結後の三角形メッシュ。
  subroutine build_parts_mesh(cfg, mesh)
    type(app_config), intent(in) :: cfg
    type(mesh_type), intent(out) :: mesh
    type(mesh_type) :: part
    real(dp), allocatable :: v0(:, :), v1(:, :), v2(:, :)
    real(dp), allocatable :: elem_epsilon_r(:), part_epsilon_r(:)
    integer(i32), allocatable :: elem_mesh_id(:), part_mesh_id(:)
    integer(i32), allocatable :: elem_surface_model(:), part_surface_model(:)
    integer(i32) :: i, mesh_id

    allocate (v0(3, 0), v1(3, 0), v2(3, 0), elem_mesh_id(0), elem_surface_model(0), elem_epsilon_r(0))
    mesh_id = 0_i32
    do i = 1, cfg%n_mesh_parts
      if (.not. cfg%mesh_parts(i)%enabled) cycle
      mesh_id = mesh_id + 1_i32
      if (trim(cfg%mesh_parts(i)%kind) == 'obj') then
        call load_obj_mesh(trim(cfg%mesh_parts(i)%path), part)
        call apply_obj_transform(part, cfg%mesh_parts(i)%scale, cfg%mesh_parts(i)%rotation, cfg%mesh_parts(i)%offset)
      else
        call build_one_template(cfg%mesh_parts(i), part)
      end if
      call append_triangles(v0, v1, v2, part%v0, part%v1, part%v2)
      if (allocated(part_mesh_id)) deallocate (part_mesh_id)
      allocate (part_mesh_id(part%nelem))
      part_mesh_id = mesh_id
      call append_mesh_ids(elem_mesh_id, part_mesh_id)
      if (allocated(part_surface_model)) deallocate (part_surface_model)
      allocate (part_surface_model(part%nelem))
      part_surface_model = surface_model_id_from_name(cfg%mesh_parts(i)%surface_model)
      call append_mesh_ids(elem_surface_model, part_surface_model)
      if (allocated(part_epsilon_r)) deallocate (part_epsilon_r)
      allocate (part_epsilon_r(part%nelem))
      part_epsilon_r = cfg%mesh_parts(i)%epsilon_r
      call append_real_values(elem_epsilon_r, part_epsilon_r)
    end do

    if (size(v0, 2) == 0) then
      error stop 'No enabled mesh.parts found in configuration.'
    end if
    call init_mesh( &
      mesh, v0, v1, v2, elem_mesh_id0=elem_mesh_id, elem_surface_model0=elem_surface_model, &
      elem_epsilon_r0=elem_epsilon_r &
      )
  end subroutine build_parts_mesh

  !> テンプレート種別に応じて形状生成ルーチンへディスパッチする。
  !! @param[in] spec 1テンプレート分の形状設定。
  !! @param[out] mesh 生成したテンプレートメッシュ。
  subroutine build_one_template(spec, mesh)
    type(template_spec), intent(in) :: spec
    type(mesh_type), intent(out) :: mesh
    logical :: cap_top, cap_bottom

    select case (trim(lower_ascii(spec%kind)))
    case ('plane')
      call make_plane(mesh, size_x=spec%size_x, size_y=spec%size_y, nx=spec%nx, ny=spec%ny, center=spec%center)
    case ('plate_hole', 'plane_hole')
      call make_plate_hole( &
        mesh, size_x=spec%size_x, size_y=spec%size_y, radius=spec%radius, n_theta=spec%n_theta, n_r=spec%n_r, &
        center=spec%center &
        )
    case ('disk')
      call make_disk( &
        mesh, radius=spec%radius, n_theta=spec%n_theta, n_r=spec%n_r, center=spec%center, &
        uniform=trim(spec%tessellation) == 'uniform' &
        )
    case ('annulus')
      call make_annulus( &
        mesh, radius=spec%radius, inner_radius=spec%inner_radius, n_theta=spec%n_theta, n_r=spec%n_r, center=spec%center, &
        uniform=trim(spec%tessellation) == 'uniform' &
        )
    case ('box')
      call make_box(mesh, size=spec%size, center=spec%center, nx=spec%nx, ny=spec%ny, nz=spec%nz)
    case ('cylinder')
      cap_top = spec%cap
      cap_bottom = spec%cap
      if (spec%has_cap_top) cap_top = spec%cap_top
      if (spec%has_cap_bottom) cap_bottom = spec%cap_bottom
      call make_cylinder( &
        mesh, radius=spec%radius, height=spec%height, n_theta=spec%n_theta, n_z=spec%n_z, &
        cap=spec%cap, center=spec%center, cap_top=cap_top, cap_bottom=cap_bottom, &
        uniform_cap=trim(spec%tessellation) == 'uniform', n_r=spec%n_r &
        )
    case ('sphere')
      if (trim(spec%tessellation) == 'uv' .or. (len_trim(spec%tessellation) == 0 .and. spec%has_n_lon_lat)) then
        call make_sphere(mesh, radius=spec%radius, n_lon=spec%n_lon, n_lat=spec%n_lat, center=spec%center)
      else
        call make_icosphere(mesh, radius=spec%radius, subdivisions=spec%subdivisions, center=spec%center)
      end if
    case default
      error stop 'Unknown template kind in config.'
    end select
  end subroutine build_one_template

  !> 既存三角形配列へ追加分を連結し、再確保後の配列へ差し替える。
  !! @param[inout] v0 累積メッシュの頂点0配列 `v0(3,n)`。
  !! @param[inout] v1 累積メッシュの頂点1配列 `v1(3,n)`。
  !! @param[inout] v2 累積メッシュの頂点2配列 `v2(3,n)`。
  !! @param[in] add_v0 追加する頂点0配列。
  !! @param[in] add_v1 追加する頂点1配列。
  !! @param[in] add_v2 追加する頂点2配列。
  subroutine append_triangles(v0, v1, v2, add_v0, add_v1, add_v2)
    real(dp), allocatable, intent(inout) :: v0(:, :), v1(:, :), v2(:, :)
    real(dp), intent(in) :: add_v0(:, :), add_v1(:, :), add_v2(:, :)
    real(dp), allocatable :: t0(:, :), t1(:, :), t2(:, :)
    integer :: n0, n1

    n0 = size(v0, 2)
    n1 = size(add_v0, 2)
    allocate (t0(3, n0 + n1), t1(3, n0 + n1), t2(3, n0 + n1))
    if (n0 > 0) then
      t0(:, 1:n0) = v0
      t1(:, 1:n0) = v1
      t2(:, 1:n0) = v2
    end if
    t0(:, n0 + 1:n0 + n1) = add_v0
    t1(:, n0 + 1:n0 + n1) = add_v1
    t2(:, n0 + 1:n0 + n1) = add_v2
    call move_alloc(t0, v0)
    call move_alloc(t1, v1)
    call move_alloc(t2, v2)
  end subroutine append_triangles

  !> 既存の要素メッシュID配列へ追加分を連結する。
  !! @param[inout] mesh_ids 累積した要素メッシュID配列。
  !! @param[in] add_ids 追加する要素メッシュID配列。
  subroutine append_mesh_ids(mesh_ids, add_ids)
    integer(i32), allocatable, intent(inout) :: mesh_ids(:)
    integer(i32), intent(in) :: add_ids(:)
    integer(i32), allocatable :: tmp(:)
    integer(i32) :: n0, n1

    n0 = size(mesh_ids)
    n1 = size(add_ids)
    allocate (tmp(n0 + n1))
    if (n0 > 0) tmp(1:n0) = mesh_ids
    if (n1 > 0) tmp(n0 + 1:n0 + n1) = add_ids
    call move_alloc(tmp, mesh_ids)
  end subroutine append_mesh_ids

  !> 既存の実数要素配列へ追加分を連結する。
  subroutine append_real_values(values, add_values)
    real(dp), allocatable, intent(inout) :: values(:)
    real(dp), intent(in) :: add_values(:)
    real(dp), allocatable :: tmp(:)
    integer(i32) :: n0, n1

    n0 = size(values)
    n1 = size(add_values)
    allocate (tmp(n0 + n1))
    if (n0 > 0) tmp(1:n0) = values
    if (n1 > 0) tmp(n0 + 1:n0 + n1) = add_values
    call move_alloc(tmp, values)
  end subroutine append_real_values

  !> 表面モデル名をメッシュ要素用の整数IDへ変換する。
  !! @param[in] name `insulator` / `conductor` / `dielectric`。
  !! @return model_id 内部表面モデルID。
  integer(i32) function surface_model_id_from_name(name) result(model_id)
    character(len=*), intent(in) :: name

    select case (trim(lower_ascii(name)))
    case ('insulator')
      model_id = surface_model_insulator
    case ('conductor')
      model_id = surface_model_conductor
    case ('dielectric')
      model_id = surface_model_dielectric
    case default
      error stop 'Unknown surface_model.'
    end select
  end function surface_model_id_from_name

end module bem_app_config_mesh_runtime

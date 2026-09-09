!> Mesh template の group、anchor、相対サイズを実座標・実寸へ変換する。
submodule(bem_app_config_authoring) bem_app_config_authoring_geometry
  use bem_string_utils, only: lower_ascii
  implicit none

contains

  !> template の group/anchor/box_fraction 指定を実座標・実寸へ変換する。
  module procedure normalize_template_high_level
  type(template_authoring_spec) :: auth
  type(mesh_group_authoring_spec) :: group
  real(dp) :: origin(3), scale
  integer :: group_idx

  auth = authoring%templates(template_idx)
  if (auth%has_group) then
    group_idx = find_group_index(authoring, auth%group)
    if (group_idx <= 0) error stop 'mesh template references an undefined group.'
    group = authoring%groups(group_idx)
    if (.not. auth%has_center_local) error stop 'grouped mesh template requires center_local.'
    if (auth%has_center .or. auth%has_placement_mode .or. auth%has_anchor .or. auth%has_offset .or. &
        auth%has_offset_frac .or. auth%has_size_mode .or. auth%has_size_frac) then
      error stop 'grouped mesh template cannot define center/placement/size high-level keys.'
    end if
    origin = resolve_group_origin(cfg, group)
    scale = resolve_group_scale(cfg, group)
    cfg%templates(template_idx)%center = origin + scale*auth%center_local
    call scale_template_lengths(cfg, template_idx, auth, scale)
    return
  end if

  call normalize_direct_template_placement(cfg, template_idx, auth)
  call normalize_direct_template_size(cfg, template_idx, auth)
  end procedure normalize_template_high_level

  !> group 名を検索する。見つからない場合は 0。
  integer function find_group_index(authoring, name) result(index_out)
    type(app_config_authoring), intent(in) :: authoring
    character(len=*), intent(in) :: name
    integer :: i

    index_out = 0
    do i = 1, authoring%n_groups
      if (trim(authoring%groups(i)%name) == trim(name)) then
        index_out = i
        return
      end if
    end do
  end function find_group_index

  !> direct template の placement_mode を解決する。
  subroutine normalize_direct_template_placement(cfg, template_idx, auth)
    type(app_config), intent(inout) :: cfg
    integer, intent(in) :: template_idx
    type(template_authoring_spec), intent(in) :: auth

    character(len=32) :: placement_mode

    placement_mode = 'absolute'
    if (auth%has_placement_mode) placement_mode = lower_ascii(trim(auth%placement_mode))
    select case (trim(placement_mode))
    case ('absolute')
      if (auth%has_anchor .or. auth%has_offset .or. auth%has_offset_frac) then
        error stop 'placement_mode="absolute" cannot use anchor/offset/offset_frac.'
      end if
    case ('box_anchor')
      if (auth%has_center) error stop 'placement_mode="box_anchor" cannot be combined with center.'
      cfg%templates(template_idx)%center = resolve_anchor_position( &
                                           cfg, auth%has_anchor, auth%anchor, auth%has_offset, auth%offset, &
                                           auth%has_offset_frac, auth%offset_frac &
                                           )
    case default
      error stop 'Unsupported mesh template placement_mode.'
    end select
  end subroutine normalize_direct_template_placement

  !> direct template の size_mode を解決する。
  subroutine normalize_direct_template_size(cfg, template_idx, auth)
    type(app_config), intent(inout) :: cfg
    integer, intent(in) :: template_idx
    type(template_authoring_spec), intent(in) :: auth

    character(len=32) :: size_mode, kind
    real(dp) :: box_size(3)

    size_mode = 'absolute'
    if (auth%has_size_mode) size_mode = lower_ascii(trim(auth%size_mode))
    select case (trim(size_mode))
    case ('absolute')
      if (auth%has_size_frac) error stop 'size_frac requires size_mode="box_fraction".'
    case ('box_fraction')
      if (.not. auth%has_size_frac) error stop 'size_mode="box_fraction" requires size_frac.'
      box_size = require_box_size(cfg, 'size_frac')
      kind = lower_ascii(trim(cfg%templates(template_idx)%kind))
      select case (trim(kind))
      case ('plane', 'plane_hole', 'plate_hole')
        call require_size_frac_len(auth, 2_i32)
        cfg%templates(template_idx)%size_x = auth%size_frac(1)*box_size(1)
        cfg%templates(template_idx)%size_y = auth%size_frac(2)*box_size(2)
      case ('box')
        call require_size_frac_len(auth, 3_i32)
        cfg%templates(template_idx)%size = auth%size_frac*box_size
      case ('sphere')
        call require_size_frac_len(auth, 1_i32)
        cfg%templates(template_idx)%radius = auth%size_frac(1)*minval(box_size)
      case ('cylinder')
        call require_size_frac_len(auth, 2_i32)
        cfg%templates(template_idx)%radius = auth%size_frac(1)*min(box_size(1), box_size(2))
        cfg%templates(template_idx)%height = auth%size_frac(2)*box_size(3)
      case default
        error stop 'size_mode="box_fraction" is not supported for this mesh template kind.'
      end select
    case default
      error stop 'Unsupported mesh template size_mode.'
    end select
  end subroutine normalize_direct_template_size

  !> size_frac の成分数を確認する。
  subroutine require_size_frac_len(auth, expected_len)
    type(template_authoring_spec), intent(in) :: auth
    integer(i32), intent(in) :: expected_len

    if (auth%size_frac_len /= expected_len) then
      error stop 'mesh template size_frac has an invalid number of components for its kind.'
    end if
  end subroutine require_size_frac_len

  !> group 原点を解決する。
  function resolve_group_origin(cfg, group) result(origin)
    type(app_config), intent(in) :: cfg
    type(mesh_group_authoring_spec), intent(in) :: group
    real(dp) :: origin(3)

    character(len=32) :: placement_mode
    real(dp) :: offset_vec(3)

    placement_mode = 'absolute'
    if (group%has_placement_mode) placement_mode = lower_ascii(trim(group%placement_mode))
    select case (trim(placement_mode))
    case ('absolute')
      if (group%has_anchor) error stop 'group placement_mode="absolute" cannot use anchor.'
      offset_vec = resolve_offset_vector( &
                   cfg, group%has_offset, group%offset, group%has_offset_frac, group%offset_frac, 'group offset' &
                   )
      origin = offset_vec
    case ('box_anchor')
      origin = resolve_anchor_position( &
               cfg, group%has_anchor, group%anchor, group%has_offset, group%offset, group%has_offset_frac, group%offset_frac &
               )
    case default
      error stop 'Unsupported mesh group placement_mode.'
    end select
  end function resolve_group_origin

  !> group scale を解決する。
  real(dp) function resolve_group_scale(cfg, group) result(scale)
    type(app_config), intent(in) :: cfg
    type(mesh_group_authoring_spec), intent(in) :: group

    if (group%has_scale .and. (group%has_scale_from .or. group%has_scale_factor)) then
      error stop 'mesh group scale cannot combine scale with scale_from/scale_factor.'
    end if
    if (group%has_scale) then
      if (group%scale <= 0.0d0) error stop 'mesh group scale must be > 0.'
      scale = group%scale
      return
    end if
    if (.not. group%has_scale_from .and. .not. group%has_scale_factor) then
      scale = 1.0d0
      return
    end if
    if (.not. group%has_scale_from .or. .not. group%has_scale_factor) then
      error stop 'mesh group scale_from and scale_factor must be specified together.'
    end if
    if (group%scale_factor <= 0.0d0) error stop 'mesh group scale_factor must be > 0.'
    scale = group%scale_factor*resolve_scale_reference(cfg, group%scale_from)
  end function resolve_group_scale

  !> scale_from の参照長を返す。
  real(dp) function resolve_scale_reference(cfg, scale_from) result(reference)
    type(app_config), intent(in) :: cfg
    character(len=*), intent(in) :: scale_from
    real(dp) :: box_size(3)
    character(len=32) :: ref_name

    box_size = require_box_size(cfg, 'scale_from')
    ref_name = lower_ascii(trim(scale_from))
    select case (trim(ref_name))
    case ('box_x')
      reference = box_size(1)
    case ('box_y')
      reference = box_size(2)
    case ('box_z')
      reference = box_size(3)
    case ('box_min_xy')
      reference = min(box_size(1), box_size(2))
    case ('box_max_xy')
      reference = max(box_size(1), box_size(2))
    case ('box_min_xyz')
      reference = minval(box_size)
    case ('box_max_xyz')
      reference = maxval(box_size)
    case default
      error stop 'Unsupported mesh group scale_from.'
    end select
  end function resolve_scale_reference

  !> anchor と offset から中心座標を返す。
  function resolve_anchor_position( &
    cfg, has_anchor, anchor, has_offset, offset, has_offset_frac, offset_frac &
    ) result(position)
    type(app_config), intent(in) :: cfg
    logical, intent(in) :: has_anchor
    character(len=*), intent(in) :: anchor
    logical, intent(in) :: has_offset
    real(dp), intent(in) :: offset(3)
    logical, intent(in) :: has_offset_frac
    real(dp), intent(in) :: offset_frac(3)
    real(dp) :: position(3)

    if (.not. has_anchor .or. len_trim(anchor) == 0) error stop 'box_anchor placement requires anchor.'
    position = resolve_anchor(cfg, anchor) + resolve_offset_vector( &
               cfg, has_offset, offset, has_offset_frac, offset_frac, 'anchor offset' &
               )
  end function resolve_anchor_position

  !> 定義済み anchor の基準座標を返す。
  function resolve_anchor(cfg, anchor) result(position)
    type(app_config), intent(in) :: cfg
    character(len=*), intent(in) :: anchor
    real(dp) :: position(3)

    real(dp) :: center(3)
    character(len=32) :: anchor_name

    call require_positive_box(cfg, 'anchor')
    center = 0.5d0*(cfg%sim%box_min + cfg%sim%box_max)
    anchor_name = lower_ascii(trim(anchor))
    select case (trim(anchor_name))
    case ('box_center')
      position = center
    case ('x_low_face_center')
      position = [cfg%sim%box_min(1), center(2), center(3)]
    case ('x_high_face_center')
      position = [cfg%sim%box_max(1), center(2), center(3)]
    case ('y_low_face_center')
      position = [center(1), cfg%sim%box_min(2), center(3)]
    case ('y_high_face_center')
      position = [center(1), cfg%sim%box_max(2), center(3)]
    case ('z_low_face_center')
      position = [center(1), center(2), cfg%sim%box_min(3)]
    case ('z_high_face_center')
      position = [center(1), center(2), cfg%sim%box_max(3)]
    case default
      error stop 'Unsupported mesh anchor.'
    end select
  end function resolve_anchor

  !> offset または offset_frac を 3 成分ベクトルへ変換する。
  function resolve_offset_vector( &
    cfg, has_offset, offset, has_offset_frac, offset_frac, context &
    ) result(offset_vec)
    type(app_config), intent(in) :: cfg
    logical, intent(in) :: has_offset
    real(dp), intent(in) :: offset(3)
    logical, intent(in) :: has_offset_frac
    real(dp), intent(in) :: offset_frac(3)
    character(len=*), intent(in) :: context
    real(dp) :: offset_vec(3)
    real(dp) :: box_size(3)

    if (has_offset .and. has_offset_frac) then
      error stop trim(context)//' cannot combine offset and offset_frac.'
    end if
    if (has_offset) then
      offset_vec = offset
    else if (has_offset_frac) then
      box_size = require_box_size(cfg, 'offset_frac')
      offset_vec = offset_frac*box_size
    else
      offset_vec = 0.0d0
    end if
  end function resolve_offset_vector

  !> group scale を template の明示済み長さパラメータへ適用する。
  subroutine scale_template_lengths(cfg, template_idx, auth, scale)
    type(app_config), intent(inout) :: cfg
    integer, intent(in) :: template_idx
    type(template_authoring_spec), intent(in) :: auth
    real(dp), intent(in) :: scale

    if (auth%has_size_x) cfg%templates(template_idx)%size_x = scale*cfg%templates(template_idx)%size_x
    if (auth%has_size_y) cfg%templates(template_idx)%size_y = scale*cfg%templates(template_idx)%size_y
    if (auth%has_size) cfg%templates(template_idx)%size = scale*cfg%templates(template_idx)%size
    if (auth%has_radius) cfg%templates(template_idx)%radius = scale*cfg%templates(template_idx)%radius
    if (auth%has_inner_radius) then
      cfg%templates(template_idx)%inner_radius = scale*cfg%templates(template_idx)%inner_radius
    end if
    if (auth%has_height) cfg%templates(template_idx)%height = scale*cfg%templates(template_idx)%height
  end subroutine scale_template_lengths

end submodule bem_app_config_authoring_geometry

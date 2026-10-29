!> mesh.parts と template 用 group の TOML 読み取り。
submodule(bem_app_config_parser) bem_app_config_parser_read_mesh
  implicit none
contains

  module procedure ensure_template_capacity
  type(template_spec), allocatable :: grown(:)
  integer :: old_capacity, new_capacity

  if (required_size <= 0) return
  if (allocated(cfg%mesh_parts)) then
    old_capacity = size(cfg%mesh_parts)
  else
    old_capacity = 0
  end if
  if (old_capacity >= required_size) return

  new_capacity = max(required_size, max(max_templates, max(1, 2*old_capacity)))
  allocate (grown(new_capacity))
  if (old_capacity > 0) grown(1:old_capacity) = cfg%mesh_parts(1:old_capacity)
  call move_alloc(grown, cfg%mesh_parts)
  end procedure ensure_template_capacity

  module procedure apply_mesh_toml_table
  type(toml_key), allocatable :: keys(:)
  integer :: ikey
  character(len=:), allocatable :: k

  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('parts')
      call read_template_array(cfg, table, keys(ikey), authoring)
    case ('groups')
      call read_mesh_groups_table(table, keys(ikey), authoring)
    case default
      error stop 'Unknown key in [mesh]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_mesh_toml_table

  module procedure read_mesh_groups_table
  type(config_toml_table), pointer :: groups_table, group_table
  type(toml_key), allocatable :: group_keys(:)
  integer :: igroup, stat

  nullify (groups_table)
  call get_value(table, key, groups_table, stat=stat)
  call require_toml_success(stat, 'mesh.groups')
  if (.not. associated(groups_table)) error stop 'mesh.groups must be a table.'

  call groups_table%get_keys(group_keys)
  call ensure_authoring_group_capacity(authoring, size(group_keys))
  authoring%n_groups = int(size(group_keys), i32)
  if (size(group_keys) > 0) authoring%groups(1:size(group_keys)) = mesh_group_authoring_spec()
  do igroup = 1, size(group_keys)
    nullify (group_table)
    call get_value(groups_table, group_keys(igroup), group_table, stat=stat)
    call require_toml_success(stat, 'mesh.groups entry')
    if (.not. associated(group_table)) error stop 'mesh.groups entries must be tables.'
    if (len_trim(group_keys(igroup)%key) > len(authoring%groups(igroup)%name)) then
      error stop 'mesh.groups name is too long.'
    end if
    authoring%groups(igroup)%name = trim(group_keys(igroup)%key)
    call apply_mesh_group_toml_table(authoring%groups(igroup), group_table)
  end do
  end procedure read_mesh_groups_table

  module procedure apply_mesh_group_toml_table
  type(toml_key), allocatable :: keys(:)
  integer :: ikey
  character(len=:), allocatable :: k

  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('placement_mode')
      call get_toml_string(table, keys(ikey), group%placement_mode, 'mesh.groups.placement_mode')
      group%placement_mode = lower_ascii(trim(group%placement_mode))
      group%has_placement_mode = .true.
    case ('anchor')
      call get_toml_string(table, keys(ikey), group%anchor, 'mesh.groups.anchor')
      group%anchor = lower_ascii(trim(group%anchor))
      group%has_anchor = .true.
    case ('offset')
      call get_toml_real_array(table, keys(ikey), group%offset, 'mesh.groups.offset')
      group%has_offset = .true.
    case ('offset_frac')
      call get_toml_real_array(table, keys(ikey), group%offset_frac, 'mesh.groups.offset_frac')
      group%has_offset_frac = .true.
    case ('scale')
      call get_toml_real(table, keys(ikey), group%scale, 'mesh.groups.scale')
      group%has_scale = .true.
    case ('scale_from')
      call get_toml_string(table, keys(ikey), group%scale_from, 'mesh.groups.scale_from')
      group%scale_from = lower_ascii(trim(group%scale_from))
      group%has_scale_from = .true.
    case ('scale_factor')
      call get_toml_real(table, keys(ikey), group%scale_factor, 'mesh.groups.scale_factor')
      group%has_scale_factor = .true.
    case default
      error stop 'Unknown key in [mesh.groups.*]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_mesh_group_toml_table

  module procedure read_template_array
  type(config_toml_array), pointer :: array
  type(config_toml_table), pointer :: child
  integer :: itemplate, n, stat

  nullify (array)
  call get_value(table, key, array, stat=stat)
  call require_toml_success(stat, 'mesh.parts')
  if (.not. associated(array)) error stop 'mesh.parts must be an array of tables.'

  n = toml_len(array)
  call ensure_template_capacity(cfg, n)
  call ensure_authoring_template_capacity(authoring, n)
  if (n > 0) cfg%mesh_parts(1:n) = template_spec()
  if (n > 0) authoring%mesh_parts(1:n) = template_authoring_spec()
  do itemplate = 1, n
    nullify (child)
    call get_value(array, itemplate, child, stat=stat)
    call require_toml_success(stat, mesh_part_config_path(itemplate))
    if (.not. associated(child)) error stop mesh_part_config_path(itemplate)//' must be a table.'
    cfg%mesh_parts(itemplate)%enabled = .true.
    call apply_template_toml_table(cfg%mesh_parts(itemplate), child, authoring%mesh_parts(itemplate), itemplate)
  end do
  cfg%n_mesh_parts = int(n, i32)
  end procedure read_template_array

  module procedure apply_template_toml_table
  type(toml_key), allocatable :: keys(:)
  integer :: ikey
  character(len=:), allocatable :: k

  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('kind', 'enabled', 'surface_side', 'surface_model', 'offset')
      continue
    case ('path', 'scale', 'rotation')
      auth%obj_only_key = k
    case default
      auth%template_only_key = k
    end select
    select case (trim(k))
    case ('path')
      call get_toml_string(table, keys(ikey), spec%path, mesh_part_config_path(part_idx)//'.path')
    case ('scale')
      call get_toml_real(table, keys(ikey), spec%scale, mesh_part_config_path(part_idx)//'.scale')
    case ('rotation')
      call get_toml_real_array(table, keys(ikey), spec%rotation, mesh_part_config_path(part_idx)//'.rotation')
    case ('enabled')
      call get_toml_logical(table, keys(ikey), spec%enabled, mesh_part_config_path(part_idx)//'.enabled')
    case ('kind')
      call get_toml_string(table, keys(ikey), spec%kind, mesh_part_config_path(part_idx)//'.kind')
      spec%kind = lower_ascii(trim(spec%kind))
    case ('surface_model')
      call get_toml_string(table, keys(ikey), spec%surface_model, mesh_part_config_path(part_idx)//'.surface_model')
      spec%surface_model = lower_ascii(trim(spec%surface_model))
    case ('surface_side')
      call get_toml_string(table, keys(ikey), spec%surface_side_policy, mesh_part_config_path(part_idx)//'.surface_side')
      spec%surface_side_policy = lower_ascii(trim(spec%surface_side_policy))
    case ('center')
      call get_toml_real_array(table, keys(ikey), spec%center, mesh_part_config_path(part_idx)//'.center')
      auth%has_center = .true.
    case ('size_x')
      call get_toml_real(table, keys(ikey), spec%size_x, mesh_part_config_path(part_idx)//'.size_x')
      auth%has_size_x = .true.
    case ('size_y')
      call get_toml_real(table, keys(ikey), spec%size_y, mesh_part_config_path(part_idx)//'.size_y')
      auth%has_size_y = .true.
    case ('size')
      call get_toml_real_array(table, keys(ikey), spec%size, mesh_part_config_path(part_idx)//'.size')
      auth%has_size = .true.
    case ('nx')
      call get_toml_int(table, keys(ikey), spec%nx, mesh_part_config_path(part_idx)//'.nx')
    case ('ny')
      call get_toml_int(table, keys(ikey), spec%ny, mesh_part_config_path(part_idx)//'.ny')
    case ('nz')
      call get_toml_int(table, keys(ikey), spec%nz, mesh_part_config_path(part_idx)//'.nz')
    case ('radius')
      call get_toml_real(table, keys(ikey), spec%radius, mesh_part_config_path(part_idx)//'.radius')
      auth%has_radius = .true.
    case ('inner_radius')
      call get_toml_real(table, keys(ikey), spec%inner_radius, mesh_part_config_path(part_idx)//'.inner_radius')
      auth%has_inner_radius = .true.
    case ('height')
      call get_toml_real(table, keys(ikey), spec%height, mesh_part_config_path(part_idx)//'.height')
      auth%has_height = .true.
    case ('n_theta')
      call get_toml_int(table, keys(ikey), spec%n_theta, mesh_part_config_path(part_idx)//'.n_theta')
    case ('n_r')
      call get_toml_int(table, keys(ikey), spec%n_r, mesh_part_config_path(part_idx)//'.n_r')
    case ('n_z')
      call get_toml_int(table, keys(ikey), spec%n_z, mesh_part_config_path(part_idx)//'.n_z')
    case ('cap')
      call get_toml_logical(table, keys(ikey), spec%cap, mesh_part_config_path(part_idx)//'.cap')
    case ('cap_top')
      call get_toml_logical(table, keys(ikey), spec%cap_top, mesh_part_config_path(part_idx)//'.cap_top')
      spec%has_cap_top = .true.
    case ('cap_bottom')
      call get_toml_logical(table, keys(ikey), spec%cap_bottom, mesh_part_config_path(part_idx)//'.cap_bottom')
      spec%has_cap_bottom = .true.
    case ('n_lon')
      call get_toml_int(table, keys(ikey), spec%n_lon, mesh_part_config_path(part_idx)//'.n_lon')
    case ('n_lat')
      call get_toml_int(table, keys(ikey), spec%n_lat, mesh_part_config_path(part_idx)//'.n_lat')
    case ('tessellation')
      call get_toml_string(table, keys(ikey), spec%tessellation, mesh_part_config_path(part_idx)//'.tessellation')
      spec%tessellation = lower_ascii(trim(spec%tessellation))
    case ('subdivisions')
      call get_toml_int(table, keys(ikey), spec%subdivisions, mesh_part_config_path(part_idx)//'.subdivisions')
    case ('group')
      call get_toml_string(table, keys(ikey), auth%group, mesh_part_config_path(part_idx)//'.group')
      auth%has_group = .true.
    case ('center_local')
      call get_toml_real_array(table, keys(ikey), auth%center_local, mesh_part_config_path(part_idx)//'.center_local')
      auth%has_center_local = .true.
    case ('placement_mode')
      call get_toml_string(table, keys(ikey), auth%placement_mode, mesh_part_config_path(part_idx)//'.placement_mode')
      auth%placement_mode = lower_ascii(trim(auth%placement_mode))
      auth%has_placement_mode = .true.
    case ('anchor')
      call get_toml_string(table, keys(ikey), auth%anchor, mesh_part_config_path(part_idx)//'.anchor')
      auth%anchor = lower_ascii(trim(auth%anchor))
      auth%has_anchor = .true.
    case ('offset')
      call get_toml_real_array(table, keys(ikey), auth%offset, mesh_part_config_path(part_idx)//'.offset')
      auth%has_offset = .true.
      spec%offset = auth%offset
    case ('offset_frac')
      call get_toml_real_array(table, keys(ikey), auth%offset_frac, mesh_part_config_path(part_idx)//'.offset_frac')
      auth%has_offset_frac = .true.
    case ('size_mode')
      call get_toml_string(table, keys(ikey), auth%size_mode, mesh_part_config_path(part_idx)//'.size_mode')
      auth%size_mode = lower_ascii(trim(auth%size_mode))
      auth%has_size_mode = .true.
    case ('size_frac')
      call get_toml_real_scalar_or_array3( &
        table, keys(ikey), auth%size_frac, auth%size_frac_len, mesh_part_config_path(part_idx)//'.size_frac' &
        )
      auth%has_size_frac = .true.
    case default
      error stop 'Unknown key: '//mesh_part_config_path(part_idx)//'.'//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_template_toml_table

end submodule bem_app_config_parser_read_mesh

!> メッシュの TOML 読み取りと、モード・表面モデルの preflight。
submodule(bem_app_config_parser) bem_app_config_parser_read_mesh
  implicit none
contains

  module procedure ensure_template_capacity
  type(template_spec), allocatable :: grown(:)
  integer :: old_capacity, new_capacity

  if (required_size <= 0) return
  if (allocated(cfg%templates)) then
    old_capacity = size(cfg%templates)
  else
    old_capacity = 0
  end if
  if (old_capacity >= required_size) return

  new_capacity = max(required_size, max(max_templates, max(1, 2*old_capacity)))
  allocate (grown(new_capacity))
  if (old_capacity > 0) grown(1:old_capacity) = cfg%templates(1:old_capacity)
  call move_alloc(grown, cfg%templates)
  end procedure ensure_template_capacity

  module procedure apply_mesh_toml_table
  type(toml_key), allocatable :: keys(:)
  integer :: ikey
  character(len=:), allocatable :: k

  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('mode')
      call get_toml_string(table, keys(ikey), cfg%mesh_mode, 'mesh.mode')
    case ('obj_path')
      call get_toml_string(table, keys(ikey), cfg%obj_path, 'mesh.obj_path')
    case ('surface_model')
      call get_toml_string(table, keys(ikey), cfg%mesh_surface_model, 'mesh.surface_model')
      cfg%mesh_surface_model = lower_ascii(trim(cfg%mesh_surface_model))
    case ('surface_side')
      call get_toml_string(table, keys(ikey), cfg%mesh_surface_side_policy, 'mesh.surface_side')
      cfg%mesh_surface_side_policy = lower_ascii(trim(cfg%mesh_surface_side_policy))
    case ('obj_scale')
      call get_toml_real(table, keys(ikey), cfg%obj_scale, 'mesh.obj_scale')
    case ('obj_rotation')
      call get_toml_real3(table, keys(ikey), cfg%obj_rotation, 'mesh.obj_rotation')
    case ('obj_offset')
      call get_toml_real3(table, keys(ikey), cfg%obj_offset, 'mesh.obj_offset')
    case ('templates')
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
      call get_toml_real3(table, keys(ikey), group%offset, 'mesh.groups.offset')
      group%has_offset = .true.
    case ('offset_frac')
      call get_toml_real3(table, keys(ikey), group%offset_frac, 'mesh.groups.offset_frac')
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
  call require_toml_success(stat, 'mesh.templates')
  if (.not. associated(array)) error stop 'mesh.templates must be an array of tables.'

  n = toml_len(array)
  call ensure_template_capacity(cfg, n)
  call ensure_authoring_template_capacity(authoring, n)
  if (n > 0) cfg%templates(1:n) = template_spec()
  if (n > 0) authoring%templates(1:n) = template_authoring_spec()
  do itemplate = 1, n
    nullify (child)
    call get_value(array, itemplate, child, stat=stat)
    call require_toml_success(stat, 'mesh.templates entry')
    if (.not. associated(child)) error stop 'mesh.templates entries must be tables.'
    cfg%templates(itemplate)%enabled = .true.
    call apply_template_toml_table(cfg%templates(itemplate), child, authoring%templates(itemplate))
  end do
  cfg%n_templates = int(n, i32)
  end procedure read_template_array

  module procedure apply_template_toml_table
  type(toml_key), allocatable :: keys(:)
  integer :: ikey
  character(len=:), allocatable :: k

  call table%get_keys(keys)
  do ikey = 1, size(keys)
    k = lower_ascii(trim(keys(ikey)%key))
    select case (trim(k))
    case ('enabled')
      call get_toml_logical(table, keys(ikey), spec%enabled, 'mesh.templates.enabled')
    case ('kind')
      call get_toml_string(table, keys(ikey), spec%kind, 'mesh.templates.kind')
    case ('surface_model')
      call get_toml_string(table, keys(ikey), spec%surface_model, 'mesh.templates.surface_model')
      spec%surface_model = lower_ascii(trim(spec%surface_model))
    case ('surface_side')
      call get_toml_string(table, keys(ikey), spec%surface_side_policy, 'mesh.templates.surface_side')
      spec%surface_side_policy = lower_ascii(trim(spec%surface_side_policy))
    case ('center')
      call get_toml_real3(table, keys(ikey), spec%center, 'mesh.templates.center')
      auth%has_center = .true.
    case ('size_x')
      call get_toml_real(table, keys(ikey), spec%size_x, 'mesh.templates.size_x')
      auth%has_size_x = .true.
    case ('size_y')
      call get_toml_real(table, keys(ikey), spec%size_y, 'mesh.templates.size_y')
      auth%has_size_y = .true.
    case ('size')
      call get_toml_real3(table, keys(ikey), spec%size, 'mesh.templates.size')
      auth%has_size = .true.
    case ('nx')
      call get_toml_int(table, keys(ikey), spec%nx, 'mesh.templates.nx')
    case ('ny')
      call get_toml_int(table, keys(ikey), spec%ny, 'mesh.templates.ny')
    case ('nz')
      call get_toml_int(table, keys(ikey), spec%nz, 'mesh.templates.nz')
    case ('radius')
      call get_toml_real(table, keys(ikey), spec%radius, 'mesh.templates.radius')
      auth%has_radius = .true.
    case ('inner_radius')
      call get_toml_real(table, keys(ikey), spec%inner_radius, 'mesh.templates.inner_radius')
      auth%has_inner_radius = .true.
    case ('height')
      call get_toml_real(table, keys(ikey), spec%height, 'mesh.templates.height')
      auth%has_height = .true.
    case ('n_theta')
      call get_toml_int(table, keys(ikey), spec%n_theta, 'mesh.templates.n_theta')
    case ('n_r')
      call get_toml_int(table, keys(ikey), spec%n_r, 'mesh.templates.n_r')
    case ('n_z')
      call get_toml_int(table, keys(ikey), spec%n_z, 'mesh.templates.n_z')
    case ('cap')
      call get_toml_logical(table, keys(ikey), spec%cap, 'mesh.templates.cap')
    case ('cap_top')
      call get_toml_logical(table, keys(ikey), spec%cap_top, 'mesh.templates.cap_top')
      spec%has_cap_top = .true.
    case ('cap_bottom')
      call get_toml_logical(table, keys(ikey), spec%cap_bottom, 'mesh.templates.cap_bottom')
      spec%has_cap_bottom = .true.
    case ('n_lon')
      call get_toml_int(table, keys(ikey), spec%n_lon, 'mesh.templates.n_lon')
    case ('n_lat')
      call get_toml_int(table, keys(ikey), spec%n_lat, 'mesh.templates.n_lat')
    case ('group')
      call get_toml_string(table, keys(ikey), auth%group, 'mesh.templates.group')
      auth%has_group = .true.
    case ('center_local')
      call get_toml_real3(table, keys(ikey), auth%center_local, 'mesh.templates.center_local')
      auth%has_center_local = .true.
    case ('placement_mode')
      call get_toml_string(table, keys(ikey), auth%placement_mode, 'mesh.templates.placement_mode')
      auth%placement_mode = lower_ascii(trim(auth%placement_mode))
      auth%has_placement_mode = .true.
    case ('anchor')
      call get_toml_string(table, keys(ikey), auth%anchor, 'mesh.templates.anchor')
      auth%anchor = lower_ascii(trim(auth%anchor))
      auth%has_anchor = .true.
    case ('offset')
      call get_toml_real3(table, keys(ikey), auth%offset, 'mesh.templates.offset')
      auth%has_offset = .true.
    case ('offset_frac')
      call get_toml_real3(table, keys(ikey), auth%offset_frac, 'mesh.templates.offset_frac')
      auth%has_offset_frac = .true.
    case ('size_mode')
      call get_toml_string(table, keys(ikey), auth%size_mode, 'mesh.templates.size_mode')
      auth%size_mode = lower_ascii(trim(auth%size_mode))
      auth%has_size_mode = .true.
    case ('size_frac')
      call get_toml_real_scalar_or_array3( &
        table, keys(ikey), auth%size_frac, auth%size_frac_len, 'mesh.templates.size_frac' &
        )
      auth%has_size_frac = .true.
    case default
      error stop 'Unknown key in [[mesh.templates]]: '//trim(keys(ikey)%key)
    end select
  end do
  end procedure apply_template_toml_table

  module procedure validate_mesh_config
  integer :: i
  cfg%mesh_mode = lower_ascii(trim(cfg%mesh_mode))
  select case (trim(cfg%mesh_mode))
  case ('auto', 'obj', 'template')
    continue
  case default
    error stop 'mesh.mode must be "auto", "obj", or "template".'
  end select
  cfg%mesh_surface_model = lower_ascii(trim(cfg%mesh_surface_model))
  select case (trim(cfg%mesh_surface_model))
  case ('insulator', 'conductor')
    continue
  case ('dielectric')
    error stop 'mesh.surface_model="dielectric" is not implemented; use "insulator" for charge accumulation.'
  case default
    error stop 'mesh.surface_model must be "insulator" or "conductor".'
  end select
  do i = 1, cfg%n_templates
    cfg%templates(i)%surface_model = lower_ascii(trim(cfg%templates(i)%surface_model))
    select case (trim(cfg%templates(i)%surface_model))
    case ('insulator', 'conductor')
      continue
    case ('dielectric')
      error stop 'mesh.templates.surface_model="dielectric" is not implemented; use "insulator".'
    case default
      error stop 'mesh.templates.surface_model must be "insulator" or "conductor".'
    end select
    if (.not. cfg%templates(i)%enabled) cycle
    select case (trim(lower_ascii(cfg%templates(i)%kind)))
    case ('plane', 'plate_hole', 'plane_hole', 'disk', 'annulus', 'box', 'cylinder', 'sphere')
      continue
    case default
      error stop 'mesh.templates.kind is unsupported.'
    end select
  end do
  end procedure validate_mesh_config

end submodule bem_app_config_parser_read_mesh

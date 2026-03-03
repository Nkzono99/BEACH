module bem_app_config
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_config, mesh_type, particles_soa
  use bem_templates, only: make_plane, make_box, make_cylinder, make_sphere
  use bem_mesh, only: init_mesh
  use bem_importers, only: load_obj_mesh
  use bem_injection, only: seed_rng, init_random_beam_particles
  implicit none

  integer, parameter :: max_templates = 8

  type :: template_spec
    logical :: enabled = .false.
    character(len=16) :: kind = 'plane'
    real(dp) :: center(3) = 0.0d0
    real(dp) :: size_x = 1.0d0
    real(dp) :: size_y = 1.0d0
    real(dp) :: size(3) = [1.0d0, 1.0d0, 1.0d0]
    integer(i32) :: nx = 1
    integer(i32) :: ny = 1
    integer(i32) :: nz = 1
    real(dp) :: radius = 0.5d0
    real(dp) :: height = 1.0d0
    integer(i32) :: n_theta = 24
    integer(i32) :: n_z = 1
    logical :: cap = .true.
    integer(i32) :: n_lon = 24
    integer(i32) :: n_lat = 12
  end type template_spec

  type :: app_config
    character(len=16) :: mesh_mode = 'auto'  ! auto / obj / template
    character(len=256) :: obj_path = 'examples/simple_plate.obj'
    integer(i32) :: n_templates = 1
    type(template_spec) :: templates(max_templates)

    integer(i32) :: rng_seed = 12345_i32
    integer(i32) :: n_particles = 256_i32
    real(dp) :: q_particle = -1.602176634d-19
    real(dp) :: m_particle = 9.10938356d-31
    real(dp) :: w_particle = 1.0d0
    real(dp) :: pos_low(3) = [-0.4d0, -0.4d0, 0.2d0]
    real(dp) :: pos_high(3) = [0.4d0, 0.4d0, 0.5d0]
    real(dp) :: drift_velocity(3) = [0.0d0, 0.0d0, -8.0d5]
    real(dp) :: temperature_k = 2.0d4

    logical :: write_output = .true.
    character(len=256) :: output_dir = 'outputs/latest'

    type(sim_config) :: sim
  end type app_config

contains

  subroutine default_app_config(cfg)
    type(app_config), intent(out) :: cfg
    cfg%sim%dt = 1.0d-9
    cfg%sim%npcls_per_step = 64
    cfg%sim%max_step = 400
    cfg%sim%tol_rel = 1.0d-8
    cfg%sim%softening = 1.0d-6
    cfg%sim%b0 = [0.0d0, 0.0d0, 0.0d0]

    cfg%templates(1)%enabled = .true.
    cfg%templates(1)%kind = 'plane'
    cfg%templates(1)%size_x = 1.0d0
    cfg%templates(1)%size_y = 1.0d0
    cfg%templates(1)%nx = 1
    cfg%templates(1)%ny = 1
    cfg%templates(1)%center = [0.0d0, 0.0d0, 0.0d0]
  end subroutine default_app_config

  subroutine load_app_config(path, cfg)
    character(len=*), intent(in) :: path
    type(app_config), intent(inout) :: cfg
    if (.not. ends_with(lower(trim(path)), '.toml')) then
      error stop 'Only TOML config is supported. Please pass a .toml file.'
    end if
    call load_toml_config(path, cfg)
  end subroutine load_app_config

  subroutine build_mesh_from_config(cfg, mesh)
    type(app_config), intent(in) :: cfg
    type(mesh_type), intent(out) :: mesh
    logical :: has_obj

    select case (trim(lower(cfg%mesh_mode)))
    case ('obj')
      call load_obj_mesh(trim(cfg%obj_path), mesh)
    case ('template')
      call build_template_mesh(cfg, mesh)
    case default
      inquire(file=trim(cfg%obj_path), exist=has_obj)
      if (has_obj) then
        call load_obj_mesh(trim(cfg%obj_path), mesh)
      else
        call build_template_mesh(cfg, mesh)
      end if
    end select
  end subroutine build_mesh_from_config

  subroutine init_particles_from_config(cfg, pcls)
    type(app_config), intent(in) :: cfg
    type(particles_soa), intent(out) :: pcls

    call seed_rng([cfg%rng_seed])
    call init_random_beam_particles(
      pcls=pcls,
      n=cfg%n_particles,
      q_particle=cfg%q_particle,
      m_particle=cfg%m_particle,
      w_particle=cfg%w_particle,
      pos_low=cfg%pos_low,
      pos_high=cfg%pos_high,
      drift_velocity=cfg%drift_velocity,
      temperature_k=cfg%temperature_k
    )
  end subroutine init_particles_from_config

  subroutine build_template_mesh(cfg, mesh)
    type(app_config), intent(in) :: cfg
    type(mesh_type), intent(out) :: mesh
    type(mesh_type) :: part
    real(dp), allocatable :: v0(:, :), v1(:, :), v2(:, :)
    integer(i32) :: i

    allocate(v0(3, 0), v1(3, 0), v2(3, 0))
    do i = 1, min(cfg%n_templates, int(max_templates, i32))
      if (.not. cfg%templates(i)%enabled) cycle
      call build_one_template(cfg%templates(i), part)
      call append_triangles(v0, v1, v2, part%v0, part%v1, part%v2)
    end do

    if (size(v0, 2) == 0) then
      error stop 'No enabled template found in configuration.'
    end if
    call init_mesh(mesh, v0, v1, v2)
  end subroutine build_template_mesh

  subroutine build_one_template(spec, mesh)
    type(template_spec), intent(in) :: spec
    type(mesh_type), intent(out) :: mesh

    select case (trim(lower(spec%kind)))
    case ('plane')
      call make_plane(mesh, size_x=spec%size_x, size_y=spec%size_y, nx=spec%nx, ny=spec%ny, center=spec%center)
    case ('box')
      call make_box(mesh, size=spec%size, center=spec%center, nx=spec%nx, ny=spec%ny, nz=spec%nz)
    case ('cylinder')
      call make_cylinder(mesh, radius=spec%radius, height=spec%height, n_theta=spec%n_theta, n_z=spec%n_z, cap=spec%cap, center=spec%center)
    case ('sphere')
      call make_sphere(mesh, radius=spec%radius, n_lon=spec%n_lon, n_lat=spec%n_lat, center=spec%center)
    case default
      error stop 'Unknown template kind in config.'
    end select
  end subroutine build_one_template

  subroutine append_triangles(v0, v1, v2, add_v0, add_v1, add_v2)
    real(dp), allocatable, intent(inout) :: v0(:, :), v1(:, :), v2(:, :)
    real(dp), intent(in) :: add_v0(:, :), add_v1(:, :), add_v2(:, :)
    real(dp), allocatable :: t0(:, :), t1(:, :), t2(:, :)
    integer :: n0, n1

    n0 = size(v0, 2)
    n1 = size(add_v0, 2)
    allocate(t0(3, n0 + n1), t1(3, n0 + n1), t2(3, n0 + n1))
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


  subroutine load_toml_config(path, cfg)
    character(len=*), intent(in) :: path
    type(app_config), intent(inout) :: cfg
    integer :: u, ios, i, t_idx
    character(len=512) :: raw, line, section

    t_idx = 0
    section = ''
    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Could not open TOML file.'

    do
      read(u, '(A)', iostat=ios) raw
      if (ios /= 0) exit
      line = strip_comment(trim(raw))
      if (len_trim(line) == 0) cycle

      if (line(1:1) == '[') then
        if (trim(line) == '[[mesh.templates]]') then
          t_idx = t_idx + 1
          if (t_idx > max_templates) error stop 'Too many mesh.templates entries.'
          cfg%templates(t_idx)%enabled = .true.
          section = 'mesh.template'
        else
          section = lower(trim(adjustl(line(2:len_trim(line)-1))))
        end if
        cycle
      end if

      select case (trim(section))
      case ('sim')
        call apply_sim_kv(cfg, line)
      case ('particles')
        call apply_particles_kv(cfg, line)
      case ('mesh')
        call apply_mesh_kv(cfg, line)
      case ('mesh.template')
        call apply_template_kv(cfg%templates(t_idx), line)
      case ('output')
        call apply_output_kv(cfg, line)
      end select
    end do
    close(u)
    if (t_idx > 0) cfg%n_templates = t_idx
  end subroutine load_toml_config

  subroutine apply_sim_kv(cfg, line)
    type(app_config), intent(inout) :: cfg
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v
    call split_key_value(line, k, v)
    select case (trim(k))
    case ('dt'); call parse_real(v, cfg%sim%dt)
    case ('npcls_per_step'); call parse_int(v, cfg%sim%npcls_per_step)
    case ('max_step'); call parse_int(v, cfg%sim%max_step)
    case ('tol_rel'); call parse_real(v, cfg%sim%tol_rel)
    case ('q_floor'); call parse_real(v, cfg%sim%q_floor)
    case ('softening'); call parse_real(v, cfg%sim%softening)
    case ('use_hybrid'); call parse_logical(v, cfg%sim%use_hybrid)
    case ('r_switch_factor'); call parse_real(v, cfg%sim%r_switch_factor)
    case ('n_sub'); call parse_int(v, cfg%sim%n_sub)
    case ('softening_factor'); call parse_real(v, cfg%sim%softening_factor)
    case ('b0'); call parse_real3(v, cfg%sim%b0)
    end select
  end subroutine apply_sim_kv

  subroutine apply_particles_kv(cfg, line)
    type(app_config), intent(inout) :: cfg
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v
    call split_key_value(line, k, v)
    select case (trim(k))
    case ('rng_seed'); call parse_int(v, cfg%rng_seed)
    case ('n_particles'); call parse_int(v, cfg%n_particles)
    case ('q_particle'); call parse_real(v, cfg%q_particle)
    case ('m_particle'); call parse_real(v, cfg%m_particle)
    case ('w_particle'); call parse_real(v, cfg%w_particle)
    case ('pos_low'); call parse_real3(v, cfg%pos_low)
    case ('pos_high'); call parse_real3(v, cfg%pos_high)
    case ('drift_velocity'); call parse_real3(v, cfg%drift_velocity)
    case ('temperature_k'); call parse_real(v, cfg%temperature_k)
    end select
  end subroutine apply_particles_kv

  subroutine apply_mesh_kv(cfg, line)
    type(app_config), intent(inout) :: cfg
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v
    call split_key_value(line, k, v)
    select case (trim(k))
    case ('mode'); call parse_string(v, cfg%mesh_mode)
    case ('obj_path'); call parse_string(v, cfg%obj_path)
    case ('n_templates'); call parse_int(v, cfg%n_templates)
    end select
  end subroutine apply_mesh_kv

  subroutine apply_template_kv(spec, line)
    type(template_spec), intent(inout) :: spec
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v
    call split_key_value(line, k, v)
    select case (trim(k))
    case ('enabled'); call parse_logical(v, spec%enabled)
    case ('kind'); call parse_string(v, spec%kind)
    case ('center'); call parse_real3(v, spec%center)
    case ('size_x'); call parse_real(v, spec%size_x)
    case ('size_y'); call parse_real(v, spec%size_y)
    case ('size'); call parse_real3(v, spec%size)
    case ('nx'); call parse_int(v, spec%nx)
    case ('ny'); call parse_int(v, spec%ny)
    case ('nz'); call parse_int(v, spec%nz)
    case ('radius'); call parse_real(v, spec%radius)
    case ('height'); call parse_real(v, spec%height)
    case ('n_theta'); call parse_int(v, spec%n_theta)
    case ('n_z'); call parse_int(v, spec%n_z)
    case ('cap'); call parse_logical(v, spec%cap)
    case ('n_lon'); call parse_int(v, spec%n_lon)
    case ('n_lat'); call parse_int(v, spec%n_lat)
    end select
  end subroutine apply_template_kv

  subroutine apply_output_kv(cfg, line)
    type(app_config), intent(inout) :: cfg
    character(len=*), intent(in) :: line
    character(len=64) :: k
    character(len=256) :: v
    call split_key_value(line, k, v)
    select case (trim(k))
    case ('write_files'); call parse_logical(v, cfg%write_output)
    case ('dir'); call parse_string(v, cfg%output_dir)
    end select
  end subroutine apply_output_kv

  subroutine split_key_value(line, key, value)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: key
    character(len=*), intent(out) :: value
    integer :: p
    p = index(line, '=')
    if (p <= 0) then
      key = ''
      value = ''
      return
    end if
    key = lower(trim(adjustl(line(:p-1))))
    value = trim(adjustl(line(p+1:)))
  end subroutine split_key_value

  subroutine parse_real(text, out)
    character(len=*), intent(in) :: text
    real(dp), intent(out) :: out
    read(text, *) out
  end subroutine parse_real

  subroutine parse_int(text, out)
    character(len=*), intent(in) :: text
    integer(i32), intent(out) :: out
    read(text, *) out
  end subroutine parse_int

  subroutine parse_logical(text, out)
    character(len=*), intent(in) :: text
    logical, intent(out) :: out
    character(len=32) :: t
    t = lower(trim(adjustl(text)))
    out = (t == 'true' .or. t == '.true.')
  end subroutine parse_logical

  subroutine parse_string(text, out)
    character(len=*), intent(in) :: text
    character(len=*), intent(out) :: out
    character(len=:), allocatable :: tmp
    tmp = trim(adjustl(text))
    if (len(tmp) >= 2 .and. tmp(1:1) == '"' .and. tmp(len(tmp):len(tmp)) == '"') then
      tmp = tmp(2:len(tmp)-1)
    end if
    out = trim(tmp)
  end subroutine parse_string

  subroutine parse_real3(text, out)
    character(len=*), intent(in) :: text
    real(dp), intent(out) :: out(3)
    character(len=256) :: t
    t = trim(adjustl(text))
    if (t(1:1) == '[') t = t(2:)
    if (t(len_trim(t):len_trim(t)) == ']') t = t(:len_trim(t)-1)
    read(t, *) out(1), out(2), out(3)
  end subroutine parse_real3

  pure function strip_comment(line) result(out)
    character(len=*), intent(in) :: line
    character(len=len(line)) :: out
    integer :: p
    p = index(line, '#')
    if (p > 0) then
      out = trim(line(:p-1))
    else
      out = trim(line)
    end if
  end function strip_comment

  pure function lower(s) result(o)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: o
    integer :: i, c
    o = s
    do i = 1, len(s)
      c = iachar(s(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) o(i:i) = achar(c + 32)
    end do
  end function lower

  pure logical function ends_with(s, suffix)
    character(len=*), intent(in) :: s, suffix
    integer :: ls, lf
    ls = len_trim(s)
    lf = len_trim(suffix)
    if (lf > ls) then
      ends_with = .false.
    else
      ends_with = (s(ls-lf+1:ls) == suffix(1:lf))
    end if
  end function ends_with

end module bem_app_config

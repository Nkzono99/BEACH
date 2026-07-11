!> BEACH TOML の高水準 authoring キーを実行時設定へ正規化する補助モジュール。
module bem_app_config_authoring
  use bem_kinds, only: dp, i32
  use bem_app_config_types, only: app_config
  use bem_string_utils, only: lower_ascii
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  integer, parameter :: default_group_capacity = 8

  type :: sim_authoring_spec
    logical :: has_box_origin = .false.
    logical :: has_box_size = .false.
    logical :: has_box_min = .false.
    logical :: has_box_max = .false.
    real(dp) :: box_origin(3) = 0.0d0
    real(dp) :: box_size(3) = 0.0d0
  end type sim_authoring_spec

  type :: field_authoring_spec
    logical :: has_element_kernel = .false.
    character(len=32) :: element_kernel = 'point'
  end type field_authoring_spec

  type :: periodic2_authoring_spec
    logical :: present = .false.
    character(len=32) :: nonzero_mode_backend = 'panel_spectral_reference'
    character(len=32) :: zero_mode_policy = 'exclude_k0'
    character(len=32) :: lower_boundary_model = 'e_bottom_zero'
    integer(i32) :: reference_mode_layers = 4_i32
    integer(i32) :: panel_quadrature_order = 12_i32
    integer(i32) :: interface_sample_n = 5_i32
    real(dp) :: interface_phi_tolerance = 1.0e-3_dp
    real(dp) :: interface_field_tolerance = 1.0e-3_dp
  end type periodic2_authoring_spec

  type :: outer_plasma_authoring_spec
    logical :: present = .false.
    character(len=32) :: model = 'linear_debye'
    character(len=32) :: photoelectron_closure = 'none'
    character(len=32) :: return_model = 'none'
    real(dp) :: interface_z = 0.0_dp
    real(dp) :: infinity_potential = 0.0_dp
    real(dp) :: debye_length = 0.0_dp
    real(dp) :: thermal_voltage = 0.0_dp
    integer(i32) :: unified_grid_points = 129_i32
    real(dp) :: accessible_fraction_tolerance = 0.1_dp
    real(dp) :: max_linearity_ratio = 0.25_dp
    real(dp) :: max_gap_ratio = 5.0_dp
    real(dp) :: max_local_charge_ratio = 50.0_dp
    integer(i32) :: photoelectron_histogram_bins = 32_i32
    real(dp) :: photoelectron_histogram_energy_max = 0.0_dp
    real(dp) :: photoelectron_ambient_charge_scale = 0.0_dp
    real(dp) :: max_photoelectron_charge_ratio = 0.1_dp
  end type outer_plasma_authoring_spec

  type :: coupling_authoring_spec
    logical :: present = .false.
    character(len=32) :: update_mode = 'explicit'
    character(len=32) :: particle_transfer_mode = 'none'
    integer(i32) :: outer_update_stride = 1_i32
    real(dp) :: field_evolution_timescale = 0.0_dp
    real(dp) :: max_frozen_field_ratio = 0.1_dp
    logical :: outer_queue_enabled = .false.
  end type coupling_authoring_spec

  type :: particle_authoring_spec
    logical :: has_inject_region_mode = .false.
    character(len=32) :: inject_region_mode = 'absolute'
    logical :: has_uv_low = .false.
    logical :: has_uv_high = .false.
    logical :: has_pos_low = .false.
    logical :: has_pos_high = .false.
    real(dp) :: uv_low(2) = 0.0d0
    real(dp) :: uv_high(2) = 0.0d0
  end type particle_authoring_spec

  type :: mesh_group_authoring_spec
    character(len=64) :: name = ''
    logical :: has_placement_mode = .false.
    character(len=32) :: placement_mode = 'absolute'
    logical :: has_anchor = .false.
    character(len=32) :: anchor = ''
    logical :: has_offset = .false.
    real(dp) :: offset(3) = 0.0d0
    logical :: has_offset_frac = .false.
    real(dp) :: offset_frac(3) = 0.0d0
    logical :: has_scale = .false.
    real(dp) :: scale = 1.0d0
    logical :: has_scale_from = .false.
    character(len=32) :: scale_from = ''
    logical :: has_scale_factor = .false.
    real(dp) :: scale_factor = 1.0d0
  end type mesh_group_authoring_spec

  type :: template_authoring_spec
    logical :: has_group = .false.
    character(len=64) :: group = ''
    logical :: has_center_local = .false.
    real(dp) :: center_local(3) = 0.0d0
    logical :: has_placement_mode = .false.
    character(len=32) :: placement_mode = 'absolute'
    logical :: has_anchor = .false.
    character(len=32) :: anchor = ''
    logical :: has_offset = .false.
    real(dp) :: offset(3) = 0.0d0
    logical :: has_offset_frac = .false.
    real(dp) :: offset_frac(3) = 0.0d0
    logical :: has_size_mode = .false.
    character(len=32) :: size_mode = 'absolute'
    logical :: has_size_frac = .false.
    integer(i32) :: size_frac_len = 0_i32
    real(dp) :: size_frac(3) = 0.0d0
    logical :: has_center = .false.
    logical :: has_size_x = .false.
    logical :: has_size_y = .false.
    logical :: has_size = .false.
    logical :: has_radius = .false.
    logical :: has_inner_radius = .false.
    logical :: has_height = .false.
  end type template_authoring_spec

  type :: app_config_authoring
    type(sim_authoring_spec) :: sim
    type(field_authoring_spec) :: field
    type(periodic2_authoring_spec) :: periodic2
    type(outer_plasma_authoring_spec) :: outer_plasma
    type(coupling_authoring_spec) :: coupling
    integer(i32) :: n_groups = 0_i32
    type(mesh_group_authoring_spec), allocatable :: groups(:)
    type(particle_authoring_spec), allocatable :: particle_species(:)
    type(template_authoring_spec), allocatable :: templates(:)
  end type app_config_authoring

  public :: app_config_authoring
  public :: sim_authoring_spec
  public :: field_authoring_spec
  public :: periodic2_authoring_spec
  public :: outer_plasma_authoring_spec
  public :: coupling_authoring_spec
  public :: particle_authoring_spec
  public :: mesh_group_authoring_spec
  public :: template_authoring_spec
  public :: init_app_config_authoring
  public :: ensure_authoring_particle_capacity
  public :: ensure_authoring_template_capacity
  public :: ensure_authoring_group_capacity
  public :: normalize_high_level_config

contains

  !> authoring overlay の配列を初期化する。
  subroutine init_app_config_authoring(authoring, template_capacity, particle_capacity)
    type(app_config_authoring), intent(out) :: authoring
    integer, intent(in) :: template_capacity
    integer, intent(in) :: particle_capacity

    allocate (authoring%groups(default_group_capacity))
    allocate (authoring%templates(max(1, template_capacity)))
    allocate (authoring%particle_species(max(1, particle_capacity)))
    authoring%groups = mesh_group_authoring_spec()
    authoring%templates = template_authoring_spec()
    authoring%particle_species = particle_authoring_spec()
    authoring%n_groups = 0_i32
  end subroutine init_app_config_authoring

  !> species authoring 配列容量を確保する。
  subroutine ensure_authoring_particle_capacity(authoring, required_size)
    type(app_config_authoring), intent(inout) :: authoring
    integer, intent(in) :: required_size
    type(particle_authoring_spec), allocatable :: grown(:)
    integer :: old_capacity, new_capacity

    if (required_size <= 0) return
    if (allocated(authoring%particle_species)) then
      old_capacity = size(authoring%particle_species)
    else
      old_capacity = 0
    end if
    if (old_capacity >= required_size) return

    new_capacity = max(required_size, max(1, 2*old_capacity))
    allocate (grown(new_capacity))
    grown = particle_authoring_spec()
    if (old_capacity > 0) grown(1:old_capacity) = authoring%particle_species(1:old_capacity)
    call move_alloc(grown, authoring%particle_species)
  end subroutine ensure_authoring_particle_capacity

  !> template authoring 配列容量を確保する。
  subroutine ensure_authoring_template_capacity(authoring, required_size)
    type(app_config_authoring), intent(inout) :: authoring
    integer, intent(in) :: required_size
    type(template_authoring_spec), allocatable :: grown(:)
    integer :: old_capacity, new_capacity

    if (required_size <= 0) return
    if (allocated(authoring%templates)) then
      old_capacity = size(authoring%templates)
    else
      old_capacity = 0
    end if
    if (old_capacity >= required_size) return

    new_capacity = max(required_size, max(1, 2*old_capacity))
    allocate (grown(new_capacity))
    grown = template_authoring_spec()
    if (old_capacity > 0) grown(1:old_capacity) = authoring%templates(1:old_capacity)
    call move_alloc(grown, authoring%templates)
  end subroutine ensure_authoring_template_capacity

  !> mesh group authoring 配列容量を確保する。
  subroutine ensure_authoring_group_capacity(authoring, required_size)
    type(app_config_authoring), intent(inout) :: authoring
    integer, intent(in) :: required_size
    type(mesh_group_authoring_spec), allocatable :: grown(:)
    integer :: old_capacity, new_capacity

    if (required_size <= 0) return
    if (allocated(authoring%groups)) then
      old_capacity = size(authoring%groups)
    else
      old_capacity = 0
    end if
    if (old_capacity >= required_size) return

    new_capacity = max(required_size, max(default_group_capacity, max(1, 2*old_capacity)))
    allocate (grown(new_capacity))
    grown = mesh_group_authoring_spec()
    if (old_capacity > 0) grown(1:old_capacity) = authoring%groups(1:old_capacity)
    call move_alloc(grown, authoring%groups)
  end subroutine ensure_authoring_group_capacity

  !> authoring overlay の高水準キーを `cfg` の実行時キーへ反映する。
  subroutine normalize_high_level_config(cfg, authoring)
    type(app_config), intent(inout) :: cfg
    type(app_config_authoring), intent(in) :: authoring
    integer :: i

    call normalize_sim_high_level(cfg, authoring%sim)

    do i = 1, cfg%n_particle_species
      call normalize_species_high_level(cfg, i, authoring%particle_species(i))
    end do

    do i = 1, cfg%n_templates
      call normalize_template_high_level(cfg, i, authoring)
    end do
  end subroutine normalize_high_level_config

  !> `sim.box_origin`/`sim.box_size` を `box_min`/`box_max` へ変換する。
  subroutine normalize_sim_high_level(cfg, sim_auth)
    type(app_config), intent(inout) :: cfg
    type(sim_authoring_spec), intent(in) :: sim_auth

    if (sim_auth%has_box_origin .neqv. sim_auth%has_box_size) then
      error stop 'sim.box_origin and sim.box_size must be specified together.'
    end if
    if (.not. sim_auth%has_box_origin) return
    if (sim_auth%has_box_min .or. sim_auth%has_box_max) then
      error stop 'sim.box_origin/box_size cannot be combined with sim.box_min/box_max.'
    end if
    if (.not. all(ieee_is_finite(sim_auth%box_origin)) .or. .not. all(ieee_is_finite(sim_auth%box_size))) then
      error stop 'sim.box_origin/box_size must contain finite values.'
    end if
    if (any(sim_auth%box_size <= 0.0d0)) then
      error stop 'sim.box_size components must be > 0.'
    end if

    cfg%sim%box_min = sim_auth%box_origin
    cfg%sim%box_max = sim_auth%box_origin + sim_auth%box_size
  end subroutine normalize_sim_high_level

  !> species の face_fraction 注入領域を実座標へ変換する。
  subroutine normalize_species_high_level(cfg, species_idx, auth)
    type(app_config), intent(inout) :: cfg
    integer, intent(in) :: species_idx
    type(particle_authoring_spec), intent(in) :: auth

    character(len=32) :: mode
    character(len=32) :: source_mode

    if (.not. auth%has_inject_region_mode) then
      if (auth%has_uv_low .or. auth%has_uv_high) then
        error stop 'particles.species uv_low/uv_high require inject_region_mode="face_fraction".'
      end if
      return
    end if

    mode = lower_ascii(trim(auth%inject_region_mode))
    source_mode = lower_ascii(trim(cfg%particle_species(species_idx)%source_mode))
    select case (trim(source_mode))
    case ('reservoir_face', 'photo_raycast')
      continue
    case default
      error stop 'inject_region_mode is only supported for reservoir_face or photo_raycast species.'
    end select

    select case (trim(mode))
    case ('absolute')
      if (auth%has_uv_low .or. auth%has_uv_high) then
        error stop 'inject_region_mode="absolute" cannot use uv_low/uv_high.'
      end if
    case ('face_fraction')
      if (auth%has_pos_low .or. auth%has_pos_high) then
        error stop 'inject_region_mode="face_fraction" cannot be combined with pos_low/pos_high.'
      end if
      if (.not. auth%has_uv_low .or. .not. auth%has_uv_high) then
        error stop 'inject_region_mode="face_fraction" requires uv_low and uv_high.'
      end if
      if (any(auth%uv_low < 0.0d0) .or. any(auth%uv_low > 1.0d0) .or. &
          any(auth%uv_high < 0.0d0) .or. any(auth%uv_high > 1.0d0)) then
        error stop 'uv_low/uv_high must be inside [0, 1].'
      end if
      if (any(auth%uv_low > auth%uv_high)) then
        error stop 'uv_low must be <= uv_high component-wise.'
      end if
      call resolve_face_fraction_region( &
        cfg%particle_species(species_idx)%inject_face, auth%uv_low, auth%uv_high, &
        cfg%sim%box_min, cfg%sim%box_max, cfg%particle_species(species_idx)%pos_low, &
        cfg%particle_species(species_idx)%pos_high &
        )
    case default
      error stop 'Unsupported inject_region_mode.'
    end select
  end subroutine normalize_species_high_level

  !> template の group/anchor/box_fraction 指定を実座標・実寸へ変換する。
  subroutine normalize_template_high_level(cfg, template_idx, authoring)
    type(app_config), intent(inout) :: cfg
    integer, intent(in) :: template_idx
    type(app_config_authoring), intent(in) :: authoring

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
  end subroutine normalize_template_high_level

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

  !> face_fraction の uv 範囲を注入面上の実座標範囲へ変換する。
  subroutine resolve_face_fraction_region(inject_face, uv_low, uv_high, box_min, box_max, pos_low, pos_high)
    character(len=*), intent(in) :: inject_face
    real(dp), intent(in) :: uv_low(2), uv_high(2), box_min(3), box_max(3)
    real(dp), intent(out) :: pos_low(3), pos_high(3)

    integer :: axis, coord1, coord2
    real(dp) :: boundary
    character(len=32) :: face

    call require_positive_bounds(box_min, box_max, 'face_fraction')
    face = lower_ascii(trim(inject_face))
    select case (trim(face))
    case ('x_low')
      axis = 1
      coord1 = 2
      coord2 = 3
      boundary = box_min(1)
    case ('x_high')
      axis = 1
      coord1 = 2
      coord2 = 3
      boundary = box_max(1)
    case ('y_low')
      axis = 2
      coord1 = 1
      coord2 = 3
      boundary = box_min(2)
    case ('y_high')
      axis = 2
      coord1 = 1
      coord2 = 3
      boundary = box_max(2)
    case ('z_low')
      axis = 3
      coord1 = 1
      coord2 = 2
      boundary = box_min(3)
    case ('z_high')
      axis = 3
      coord1 = 1
      coord2 = 2
      boundary = box_max(3)
    case default
      error stop 'Invalid inject_face for face_fraction.'
    end select

    pos_low = box_min
    pos_high = box_max
    pos_low(axis) = boundary
    pos_high(axis) = boundary
    pos_low(coord1) = box_min(coord1) + uv_low(1)*(box_max(coord1) - box_min(coord1))
    pos_high(coord1) = box_min(coord1) + uv_high(1)*(box_max(coord1) - box_min(coord1))
    pos_low(coord2) = box_min(coord2) + uv_low(2)*(box_max(coord2) - box_min(coord2))
    pos_high(coord2) = box_min(coord2) + uv_high(2)*(box_max(coord2) - box_min(coord2))
  end subroutine resolve_face_fraction_region

  !> 現在の sim box size を返す。
  function require_box_size(cfg, context) result(box_size)
    type(app_config), intent(in) :: cfg
    character(len=*), intent(in) :: context
    real(dp) :: box_size(3)

    call require_positive_box(cfg, context)
    box_size = cfg%sim%box_max - cfg%sim%box_min
  end function require_box_size

  !> cfg の box が有限かつ正の大きさを持つことを確認する。
  subroutine require_positive_box(cfg, context)
    type(app_config), intent(in) :: cfg
    character(len=*), intent(in) :: context

    call require_positive_bounds(cfg%sim%box_min, cfg%sim%box_max, context)
  end subroutine require_positive_box

  !> box_min/box_max が有限かつ正の大きさを持つことを確認する。
  subroutine require_positive_bounds(box_min, box_max, context)
    real(dp), intent(in) :: box_min(3), box_max(3)
    character(len=*), intent(in) :: context

    if (.not. all(ieee_is_finite(box_min)) .or. .not. all(ieee_is_finite(box_max))) then
      error stop trim(context)//' requires finite sim.box_min/box_max.'
    end if
    if (any(box_max <= box_min)) then
      error stop trim(context)//' requires positive box dimensions.'
    end if
  end subroutine require_positive_bounds

end module bem_app_config_authoring

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
    logical :: has_reservoir_potential_model = .false.
    logical :: has_open_boundary_model = .false.
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
    character(len=32) :: model = 'none'
    character(len=32) :: kinetic_closure = 'absorbing_maxwellian'
    character(len=32) :: zhao_branch = 'auto'
    real(dp) :: photoelectron_source_scale = 1.0_dp
    character(len=32) :: photoelectron_density_model = 'none'
    character(len=32) :: return_model = 'none'
    real(dp) :: interface_z = 0.0_dp
    real(dp) :: debye_length = 0.0_dp
    real(dp) :: thermal_voltage = 0.0_dp
    real(dp) :: max_gap_ratio = 5.0_dp
    real(dp) :: max_local_charge_ratio = 50.0_dp
  end type outer_plasma_authoring_spec

  type :: coupling_authoring_spec
    logical :: present = .false.
    character(len=32) :: update_mode = 'explicit'
    character(len=32) :: particle_transfer_mode = 'none'
    character(len=32) :: steady_start_mode = 'none'
    integer(i32) :: steady_start_mesh_id = 1_i32
    integer(i32) :: outer_update_stride = 1_i32
    real(dp) :: field_evolution_timescale = 0.0_dp
    real(dp) :: max_frozen_field_ratio = 0.1_dp
    logical :: outer_queue_enabled = .false.
  end type coupling_authoring_spec

  !> `[external_boundary.field]` の公開 authoring 設定。
  !! runtime の return/interface ownership は含めず、外部場そのものの設定だけを保持する。
  type :: external_boundary_field_authoring_spec
    logical :: present = .false.
    logical :: has_model = .false.
    logical :: has_non_model_key = .false.
    logical :: has_kinetic_closure = .false.
    logical :: has_zhao_option = .false.
    logical :: has_photoelectron_density_model = .false.
    character(len=32) :: model = 'none'
    character(len=32) :: kinetic_closure = 'absorbing_maxwellian'
    character(len=32) :: zhao_branch = 'auto'
    real(dp) :: photoelectron_source_scale = 1.0_dp
    character(len=32) :: photoelectron_density_model = 'none'
    real(dp) :: debye_length = 0.0_dp
    real(dp) :: thermal_voltage = 0.0_dp
    real(dp) :: max_gap_ratio = 5.0_dp
    real(dp) :: max_local_charge_ratio = 50.0_dp
  end type external_boundary_field_authoring_spec

  !> `[external_boundary.particles]` の公開 authoring 設定。
  !! runtime の transfer/queue owner は `mode` から導出する。
  type :: external_boundary_particles_authoring_spec
    logical :: present = .false.
    logical :: has_mode = .false.
    logical :: has_coupling_option = .false.
    logical :: has_steady_start_mode = .false.
    logical :: has_steady_start_mesh_id = .false.
    logical :: has_outer_update_stride = .false.
    logical :: has_time_guard_option = .false.
    character(len=32) :: mode = 'local_source'
    character(len=32) :: inflow_model = 'auto'
    character(len=32) :: steady_start_mode = 'none'
    integer(i32) :: steady_start_mesh_id = 1_i32
    integer(i32) :: outer_update_stride = 1_i32
    real(dp) :: field_evolution_timescale = 0.0_dp
    real(dp) :: max_frozen_field_ratio = 0.1_dp
  end type external_boundary_particles_authoring_spec

  !> outer interface が所有しない通常 open 面の公開 authoring 設定。
  type :: external_boundary_ordinary_open_authoring_spec
    character(len=32) :: model = 'escape'
  end type external_boundary_ordinary_open_authoring_spec

  !> 外部境界の公開 facade。field と particles は必須、ordinary_open は省略可能。
  type :: external_boundary_authoring_spec
    logical :: present = .false.
    type(external_boundary_field_authoring_spec) :: field
    type(external_boundary_particles_authoring_spec) :: particles
    type(external_boundary_ordinary_open_authoring_spec) :: ordinary_open
  end type external_boundary_authoring_spec

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
    type(external_boundary_authoring_spec) :: external_boundary
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
  public :: external_boundary_field_authoring_spec
  public :: external_boundary_particles_authoring_spec
  public :: external_boundary_ordinary_open_authoring_spec
  public :: external_boundary_authoring_spec
  public :: particle_authoring_spec
  public :: mesh_group_authoring_spec
  public :: template_authoring_spec
  public :: init_app_config_authoring
  public :: ensure_authoring_particle_capacity
  public :: ensure_authoring_template_capacity
  public :: ensure_authoring_group_capacity
  public :: normalize_high_level_config
  public :: lower_external_boundary_authoring

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

  !> `[external_boundary.*]` facade を既存 runtime 設定語彙へ lower する。
  !!
  !! facade は authoring 専用であり、runtime 型と外部境界 contract は変更しない。
  !! 旧 selector と facade の併記は、同じ値であっても ownership が曖昧になるため拒否する。
  subroutine lower_external_boundary_authoring(cfg, authoring)
    type(app_config), intent(inout) :: cfg
    type(app_config_authoring), intent(in) :: authoring

    character(len=32) :: field_model, closure, branch
    character(len=32) :: particle_mode, inflow_model
    character(len=32) :: ordinary_open_model
    logical :: profile_owned_inflow

    if (.not. authoring%external_boundary%present) return

    if (authoring%outer_plasma%present) then
      error stop '[external_boundary] cannot be combined with legacy [outer_plasma].'
    end if
    if (authoring%coupling%present) then
      error stop '[external_boundary] cannot be combined with legacy [coupling].'
    end if
    if (authoring%sim%has_reservoir_potential_model) then
      error stop '[external_boundary] cannot be combined with sim.reservoir_potential_model.'
    end if
    if (authoring%sim%has_open_boundary_model) then
      error stop '[external_boundary] cannot be combined with sim.open_boundary_model.'
    end if
    if (.not. authoring%external_boundary%field%present) then
      error stop '[external_boundary.field] is required when [external_boundary] is used.'
    end if
    if (.not. authoring%external_boundary%particles%present) then
      error stop '[external_boundary.particles] is required when [external_boundary] is used.'
    end if
    if (.not. authoring%external_boundary%field%has_model) then
      error stop 'external_boundary.field.model is required.'
    end if
    if (.not. authoring%external_boundary%particles%has_mode) then
      error stop 'external_boundary.particles.mode is required.'
    end if

    field_model = lower_ascii(trim(authoring%external_boundary%field%model))
    closure = lower_ascii(trim(authoring%external_boundary%field%kinetic_closure))
    branch = lower_ascii(trim(authoring%external_boundary%field%zhao_branch))
    particle_mode = lower_ascii(trim(authoring%external_boundary%particles%mode))
    inflow_model = lower_ascii(trim(authoring%external_boundary%particles%inflow_model))
    ordinary_open_model = lower_ascii(trim(authoring%external_boundary%ordinary_open%model))

    select case (trim(field_model))
    case ('none', 'kinetic_1d')
      continue
    case default
      error stop 'external_boundary.field.model is not supported.'
    end select
    select case (trim(ordinary_open_model))
    case ('escape', 'potential_barrier')
      continue
    case default
      error stop 'external_boundary.ordinary_open.model must be "escape" or "potential_barrier".'
    end select
    if (trim(field_model) == 'none' .and. authoring%external_boundary%field%has_non_model_key) then
      error stop 'external_boundary.field.model="none" does not accept additional field options.'
    end if
    if (trim(field_model) == 'none' .and. trim(particle_mode) == 'local_source' .and. &
        authoring%external_boundary%particles%has_coupling_option) then
      error stop 'field.model="none" with particles.mode="local_source" does not accept coupling options.'
    end if
    if (trim(field_model) /= 'none' .and. .not. authoring%sim%has_box_max .and. &
        .not. authoring%sim%has_box_origin) then
      error stop 'An active external_boundary.field.model requires explicit sim.box_max or sim.box_origin/box_size.'
    end if
    call validate_external_boundary_option_applicability( &
      authoring%external_boundary%field, authoring%external_boundary%particles, &
      field_model, closure, particle_mode &
      )

    cfg%outer_plasma%model = field_model
    cfg%outer_plasma%kinetic_closure = closure
    cfg%outer_plasma%zhao_branch = branch
    cfg%outer_plasma%photoelectron_source_scale = authoring%external_boundary%field%photoelectron_source_scale
    cfg%outer_plasma%photoelectron_density_model = authoring%external_boundary%field%photoelectron_density_model
    ! inactive field は legacy simple config と同じ fingerprint になるよう runtime 既定値を保持する。
    if (trim(field_model) /= 'none') cfg%outer_plasma%interface_z = cfg%sim%box_max(3)
    cfg%outer_plasma%debye_length = authoring%external_boundary%field%debye_length
    cfg%outer_plasma%thermal_voltage = authoring%external_boundary%field%thermal_voltage
    cfg%outer_plasma%max_gap_ratio = authoring%external_boundary%field%max_gap_ratio
    cfg%outer_plasma%max_local_charge_ratio = authoring%external_boundary%field%max_local_charge_ratio

    cfg%coupling%update_mode = 'explicit'
    cfg%coupling%steady_start_mode = authoring%external_boundary%particles%steady_start_mode
    cfg%coupling%steady_start_mesh_id = authoring%external_boundary%particles%steady_start_mesh_id
    cfg%coupling%outer_update_stride = authoring%external_boundary%particles%outer_update_stride
    cfg%coupling%field_evolution_timescale = &
      authoring%external_boundary%particles%field_evolution_timescale
    cfg%coupling%max_frozen_field_ratio = authoring%external_boundary%particles%max_frozen_field_ratio
    cfg%coupling%outer_queue_enabled = .false.
    cfg%sim%open_boundary_model = ordinary_open_model

    profile_owned_inflow = .false.
    select case (trim(particle_mode))
    case ('local_source')
      cfg%outer_plasma%return_model = 'none'
      cfg%coupling%particle_transfer_mode = 'none'
    case ('same_batch')
      select case (trim(field_model))
      case ('kinetic_1d')
        cfg%outer_plasma%return_model = 'kinetic_1d_profile_return'
        cfg%coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
        profile_owned_inflow = .true.
      case default
        error stop 'external_boundary.particles.mode="same_batch" requires an outer field model.'
      end select
    case ('zhao_queue')
      if (trim(field_model) /= 'kinetic_1d' .or. trim(closure) /= 'zhao_charge_driven') then
        error stop 'external_boundary.particles.mode="zhao_queue" requires kinetic_1d zhao_charge_driven.'
      end if
      if (trim(branch) /= 'auto') then
        error stop 'external_boundary.particles.mode="zhao_queue" requires external_boundary.field.zhao_branch="auto".'
      end if
      cfg%outer_plasma%return_model = 'kinetic_1d_profile_return'
      cfg%coupling%particle_transfer_mode = 'electrostatic_1d_instant_return'
      cfg%coupling%outer_queue_enabled = .true.
      profile_owned_inflow = .true.
    case default
      error stop 'external_boundary.particles.mode is not supported.'
    end select

    if (profile_owned_inflow .and. trim(inflow_model) /= 'auto') then
      error stop 'Profile-owned tracked inflow requires external_boundary.particles.inflow_model="auto".'
    end if
    select case (trim(inflow_model))
    case ('auto', 'source_vdf')
      cfg%sim%reservoir_potential_model = 'none'
    case ('infinity_barrier')
      cfg%sim%reservoir_potential_model = 'infinity_barrier'
    case default
      error stop 'external_boundary.particles.inflow_model is not supported.'
    end select
  end subroutine lower_external_boundary_authoring

  !> facade で明示されたキーが、選択した field/particle model に実際に適用されることを確認する。
  !!
  !! 省略と no-op の明示を区別するため、authoring DTO の presence flag を使う。
  subroutine validate_external_boundary_option_applicability(field, particles, field_model, closure, particle_mode)
    type(external_boundary_field_authoring_spec), intent(in) :: field
    type(external_boundary_particles_authoring_spec), intent(in) :: particles
    character(len=*), intent(in) :: field_model, closure, particle_mode
    logical :: zhao_field, steady_start_available

    zhao_field = trim(field_model) == 'kinetic_1d' .and. trim(closure) == 'zhao_charge_driven'
    steady_start_available = zhao_field .and. trim(particle_mode) == 'same_batch'

    if (field%has_kinetic_closure .and. trim(field_model) /= 'kinetic_1d') then
      error stop 'external_boundary.field.kinetic_closure requires field.model="kinetic_1d".'
    end if
    if (field%has_zhao_option .and. .not. zhao_field) then
      error stop 'Zhao field options require kinetic_1d with kinetic_closure="zhao_charge_driven".'
    end if
    if (field%has_photoelectron_density_model .and. trim(field_model) /= 'kinetic_1d') then
      error stop 'external_boundary.field.photoelectron_density_model requires kinetic_1d.'
    end if
    if (field%has_photoelectron_density_model .and. &
        trim(field%photoelectron_density_model) == 'kinetic_mean' .and. &
        trim(closure) /= 'absorbing_maxwellian') then
      error stop 'photoelectron_density_model="kinetic_mean" requires absorbing-Maxwellian kinetic_1d.'
    end if
    if (field%has_photoelectron_density_model .and. &
        trim(field%photoelectron_density_model) == 'linearized_mean' .and. &
        trim(closure) /= 'ambient_linear_debye') then
      error stop 'photoelectron_density_model="linearized_mean" requires ambient_linear_debye.'
    end if
    if (particles%has_time_guard_option .and. &
        trim(particle_mode) /= 'same_batch' .and. trim(particle_mode) /= 'zhao_queue') then
      error stop 'External field time-guard options require particles.mode="same_batch" or "zhao_queue".'
    end if
    if (particles%has_outer_update_stride .and. &
        (trim(field_model) /= 'kinetic_1d' .or. trim(particle_mode) == 'zhao_queue')) then
      error stop 'outer_update_stride requires a kinetic_1d field and a non-queue particle mode.'
    end if
    if ((particles%has_steady_start_mode .or. particles%has_steady_start_mesh_id) .and. &
        .not. steady_start_available) then
      error stop 'Steady-start options require Zhao kinetic_1d with particles.mode="same_batch".'
    end if
    if (particles%has_steady_start_mesh_id .and. trim(particles%steady_start_mode) /= 'zhao_floating') then
      error stop 'steady_start_mesh_id requires steady_start_mode="zhao_floating".'
    end if
  end subroutine validate_external_boundary_option_applicability

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

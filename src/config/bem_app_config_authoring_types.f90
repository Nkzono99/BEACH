!> TOML authoring overlay の型と default。実行時設定への変換処理は持たない。
module bem_app_config_authoring_types
  use bem_kinds, only: dp, i32
  use bem_app_config_types, only: particle_bc_inherit
  implicit none
  private

  type :: periodic2_authoring_spec
    logical :: present = .false.
    character(len=32) :: nonzero_mode_backend = 'panel_spectral_reference'
    logical :: has_nonzero_mode_backend = .false.
    character(len=32) :: zero_mode_policy = 'exclude_k0'
    logical :: has_zero_mode_policy = .false.
    character(len=32) :: lower_boundary_model = 'e_bottom_zero'
    logical :: has_lower_boundary_model = .false.
    integer(i32) :: reference_mode_layers = 4_i32
    integer(i32) :: panel_quadrature_order = 12_i32
    real(dp) :: max_nonzero_mode_potential_step = 0.0_dp
  end type periodic2_authoring_spec

  !> 計算領域と周期 topology の公開 authoring 設定。
  type :: domain_authoring_spec
    logical :: present = .false.
    logical :: has_box_origin = .false.
    logical :: has_box_size = .false.
    logical :: has_box_min = .false.
    logical :: has_box_max = .false.
    real(dp) :: box_origin(3) = 0.0_dp
    real(dp) :: box_size(3) = 0.0_dp
    real(dp) :: box_min(3) = 0.0_dp
    real(dp) :: box_max(3) = 0.0_dp
    logical :: periodic_axis(3) = .false.
  end type domain_authoring_spec

  !> 場の境界 closure の公開 authoring 設定。
  type :: field_boundary_authoring_spec
    logical :: present = .false.
    character(len=32) :: mode = 'free'
  end type field_boundary_authoring_spec

  !> tracked particle のglobal面作用を保持する。
  type :: particle_boundary_authoring_spec
    logical :: present = .false.
    integer(i32) :: low(3) = [particle_bc_inherit, particle_bc_inherit, particle_bc_inherit]
    integer(i32) :: high(3) = [particle_bc_inherit, particle_bc_inherit, particle_bc_inherit]
    character(len=32) :: ordinary_open_model = 'escape'
  end type particle_boundary_authoring_spec

  !> 局所 reservoir inflow の公開 authoring 設定。
  type :: reservoir_authoring_spec
    logical :: present = .false.
    character(len=32) :: inflow_model = 'source_vdf'
    real(dp) :: phi_infty = 0.0_dp
    integer(i32) :: face_potential_grid_n = 3_i32
  end type reservoir_authoring_spec

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
    type(domain_authoring_spec) :: domain
    type(field_boundary_authoring_spec) :: field_boundary
    type(particle_boundary_authoring_spec) :: particle_boundary
    type(reservoir_authoring_spec) :: reservoir
    type(periodic2_authoring_spec) :: periodic2
    integer(i32) :: n_groups = 0_i32
    type(mesh_group_authoring_spec), allocatable :: groups(:)
    type(particle_authoring_spec), allocatable :: particle_species(:)
    type(template_authoring_spec), allocatable :: templates(:)
  end type app_config_authoring

  public :: app_config_authoring
  public :: domain_authoring_spec
  public :: field_boundary_authoring_spec
  public :: particle_boundary_authoring_spec
  public :: reservoir_authoring_spec
  public :: periodic2_authoring_spec
  public :: particle_authoring_spec
  public :: mesh_group_authoring_spec
  public :: template_authoring_spec

end module bem_app_config_authoring_types

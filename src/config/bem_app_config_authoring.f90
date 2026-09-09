!> Authoring overlay の公開入口。配列管理と領域別変換の入口を持つ。
module bem_app_config_authoring
  use bem_kinds, only: dp, i32
  use bem_app_config_types, only: app_config
  use bem_app_config_authoring_types, only: &
    app_config_authoring, domain_authoring_spec, field_boundary_authoring_spec, &
    particle_boundary_authoring_spec, reservoir_authoring_spec, periodic2_authoring_spec, &
    particle_authoring_spec, mesh_group_authoring_spec, template_authoring_spec
  implicit none
  private

  integer, parameter :: default_group_capacity = 8

  public :: app_config_authoring
  public :: domain_authoring_spec
  public :: field_boundary_authoring_spec
  public :: particle_boundary_authoring_spec
  public :: reservoir_authoring_spec
  public :: periodic2_authoring_spec
  public :: particle_authoring_spec
  public :: mesh_group_authoring_spec
  public :: template_authoring_spec
  public :: init_app_config_authoring
  public :: ensure_authoring_particle_capacity
  public :: ensure_authoring_template_capacity
  public :: ensure_authoring_group_capacity
  public :: normalize_high_level_config
  public :: lower_boundary_authoring

  interface
    !> domain / field / particle / reservoir の公開設定をruntime設定へ lower する。
    module subroutine lower_boundary_authoring(cfg, authoring)
      type(app_config), intent(inout) :: cfg
      type(app_config_authoring), intent(in) :: authoring
    end subroutine lower_boundary_authoring

    !> species の face_fraction 注入領域を実座標へ変換する。
    module subroutine normalize_species_high_level(cfg, species_idx, auth)
      type(app_config), intent(inout) :: cfg
      integer, intent(in) :: species_idx
      type(particle_authoring_spec), intent(in) :: auth
    end subroutine normalize_species_high_level

    !> template の group/anchor/box_fraction 指定を実座標・実寸へ変換する。
    module subroutine normalize_template_high_level(cfg, template_idx, authoring)
      type(app_config), intent(inout) :: cfg
      integer, intent(in) :: template_idx
      type(app_config_authoring), intent(in) :: authoring
    end subroutine normalize_template_high_level

    !> 現在の sim box size を返す。
    module function require_box_size(cfg, context) result(box_size)
      type(app_config), intent(in) :: cfg
      character(len=*), intent(in) :: context
      real(dp) :: box_size(3)
    end function require_box_size

    !> cfg の box が有限かつ正の大きさを持つことを確認する。
    module subroutine require_positive_box(cfg, context)
      type(app_config), intent(in) :: cfg
      character(len=*), intent(in) :: context
    end subroutine require_positive_box

    !> box_min/box_max が有限かつ正の大きさを持つことを確認する。
    module subroutine require_positive_bounds(box_min, box_max, context)
      real(dp), intent(in) :: box_min(3), box_max(3)
      character(len=*), intent(in) :: context
    end subroutine require_positive_bounds
  end interface

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

    do i = 1, cfg%n_particle_species
      call normalize_species_high_level(cfg, i, authoring%particle_species(i))
    end do

    do i = 1, cfg%n_templates
      call normalize_template_high_level(cfg, i, authoring)
    end do
  end subroutine normalize_high_level_config

end module bem_app_config_authoring

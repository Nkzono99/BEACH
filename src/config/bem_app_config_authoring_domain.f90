!> 計算領域、境界、reservoir 設定の正規化と box の幾何条件。
submodule(bem_app_config_authoring) bem_app_config_authoring_domain
  use bem_types, only: bc_open, bc_periodic
  use bem_string_utils, only: lower_ascii
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

contains

  !> domain / field / particle / reservoir の公開設定をruntime設定へ lower する。
  module procedure lower_boundary_authoring
  integer :: axis
  character(len=32) :: mode

  if (authoring%domain%present) then
    cfg%sim%use_box = .true.
    if (authoring%domain%has_box_origin .neqv. authoring%domain%has_box_size) then
      error stop 'domain.box_origin and domain.box_size must be specified together.'
    end if
    if (authoring%domain%has_box_min .neqv. authoring%domain%has_box_max) then
      error stop 'domain.box_min and domain.box_max must be specified together.'
    end if
    if (authoring%domain%has_box_origin .and. authoring%domain%has_box_min) then
      error stop 'domain.box_origin/box_size cannot be combined with domain.box_min/box_max.'
    end if
    if (authoring%domain%has_box_origin) then
      cfg%sim%box_min = authoring%domain%box_origin
      cfg%sim%box_max = authoring%domain%box_origin + authoring%domain%box_size
    else if (authoring%domain%has_box_min) then
      cfg%sim%box_min = authoring%domain%box_min
      cfg%sim%box_max = authoring%domain%box_max
    else
      error stop '[domain] requires box_min/box_max or box_origin/box_size.'
    end if
    cfg%sim%bc_low = bc_open
    cfg%sim%bc_high = bc_open
    do axis = 1, 3
      if (authoring%domain%periodic_axis(axis)) then
        cfg%sim%bc_low(axis) = bc_periodic
        cfg%sim%bc_high(axis) = bc_periodic
      end if
    end do
  end if

  if (authoring%field_boundary%present) then
    cfg%sim%field_bc_mode = lower_ascii(trim(authoring%field_boundary%mode))
  end if

  if (authoring%particle_boundary%present) then
    cfg%particle_boundary_low = authoring%particle_boundary%low
    cfg%particle_boundary_high = authoring%particle_boundary%high
    cfg%sim%open_boundary_model = lower_ascii(trim(authoring%particle_boundary%ordinary_open_model))
  end if

  if (.not. authoring%reservoir%present) return
  mode = lower_ascii(trim(authoring%reservoir%inflow_model))
  select case (trim(mode))
  case ('source_vdf')
    cfg%sim%reservoir_potential_model = 'none'
  case ('infinity_barrier')
    cfg%sim%reservoir_potential_model = 'infinity_barrier'
  case default
    error stop 'reservoir.inflow_model is not supported.'
  end select
  cfg%sim%phi_infty = authoring%reservoir%phi_infty
  cfg%sim%injection_face_phi_grid_n = authoring%reservoir%face_potential_grid_n
  end procedure lower_boundary_authoring

  !> 現在の sim box size を返す。
  module procedure require_box_size
  call require_positive_box(cfg, context)
  box_size = cfg%sim%box_max - cfg%sim%box_min
  end procedure require_box_size

  !> cfg の box が有限かつ正の大きさを持つことを確認する。
  module procedure require_positive_box
  call require_positive_bounds(cfg%sim%box_min, cfg%sim%box_max, context)
  end procedure require_positive_box

  !> box_min/box_max が有限かつ正の大きさを持つことを確認する。
  module procedure require_positive_bounds
  if (.not. all(ieee_is_finite(box_min)) .or. .not. all(ieee_is_finite(box_max))) then
    error stop trim(context)//' requires finite domain.box_min/box_max.'
  end if
  if (any(box_max <= box_min)) then
    error stop trim(context)//' requires positive box dimensions.'
  end if
  end procedure require_positive_bounds

end submodule bem_app_config_authoring_domain

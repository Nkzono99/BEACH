!> 粒子源の face_fraction 入力を注入面上の実座標へ変換する。
submodule(bem_app_config_authoring) bem_app_config_authoring_sources
  use bem_string_utils, only: lower_ascii
  implicit none

contains

  !> species の face_fraction 注入領域を実座標へ変換する。
  module procedure normalize_species_high_level
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
  end procedure normalize_species_high_level

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

end submodule bem_app_config_authoring_sources

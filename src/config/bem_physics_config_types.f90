!> 場・periodic2・panel の型付き設定と旧 `[sim]` 設定の正規化を定義する。
module bem_physics_config_types
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_config, bc_periodic
  use bem_string_utils, only: lower_ascii
  implicit none
  private

  integer(i32), parameter, public :: physics_config_ok = 0_i32
  integer(i32), parameter, public :: physics_config_invalid_combination = 1_i32
  integer(i32), parameter, public :: physics_config_unavailable = 2_i32

  type, public :: field_physics_config
    character(len=32) :: backend = 'auto'
    character(len=32) :: normalization = 'si'
  end type field_physics_config

  type, public :: periodic2_physics_config
    character(len=32) :: nonzero_mode_backend = 'not_applicable'
    character(len=32) :: zero_mode_policy = 'not_applicable'
    character(len=32) :: lower_boundary_model = 'not_applicable'
    integer(i32) :: reference_mode_layers = 4_i32
    integer(i32) :: panel_quadrature_order = 12_i32
    real(dp) :: max_nonzero_mode_potential_step = 0.0_dp
  end type periodic2_physics_config

  type, public :: panel_kernel_config
    character(len=32) :: kernel_id = 'triangle_p0_exact_auto'
    character(len=32) :: surface_side_policy = 'per_element'
  end type panel_kernel_config

  public :: normalize_legacy_physics_config
  public :: validate_phase1_panel_config
  public :: validate_active_physics_config

contains

  !> 現行 `sim_config` を型付き場設定へ写像する。
  subroutine normalize_legacy_physics_config(sim, field, periodic2, panel)
    type(sim_config), intent(in) :: sim
    type(field_physics_config), intent(out) :: field
    type(periodic2_physics_config), intent(out) :: periodic2
    type(panel_kernel_config), intent(out) :: panel
    character(len=32) :: bc_mode, far_mode

    field = field_physics_config()
    periodic2 = periodic2_physics_config()
    panel = panel_kernel_config()

    field%backend = lower_ascii(trim(sim%field_solver))
    field%normalization = lower_ascii(trim(sim%field_normalization))
    select case (trim(field%backend))
    case ('direct')
      panel%kernel_id = 'triangle_p0_exact_direct'
    case ('treecode')
      panel%kernel_id = 'triangle_p0_exact_tree_near'
    case ('fmm')
      panel%kernel_id = 'triangle_p0_exact_p2m_near'
    case default
      panel%kernel_id = 'triangle_p0_exact_auto'
    end select

    bc_mode = lower_ascii(trim(sim%field_bc_mode))
    if (trim(bc_mode) /= 'periodic2') return

    far_mode = lower_ascii(trim(sim%field_periodic_far_correction))
    select case (trim(far_mode))
    case ('none', 'auto')
      periodic2%nonzero_mode_backend = 'legacy_finite_images'
      periodic2%zero_mode_policy = 'legacy_not_decomposed'
      periodic2%lower_boundary_model = 'legacy_implicit'
    case ('cached_kneq0')
      periodic2%nonzero_mode_backend = 'cached_kneq0'
      periodic2%zero_mode_policy = 'exclude_k0'
      periodic2%lower_boundary_model = 'e_bottom_zero'
    case default
      periodic2%nonzero_mode_backend = 'invalid'
      periodic2%zero_mode_policy = 'invalid'
      periodic2%lower_boundary_model = 'invalid'
    end select
  end subroutine normalize_legacy_physics_config

  !> triangle P0 direct/treecode/FMM kernel の solver/boundary 契約を検証する。
  subroutine validate_phase1_panel_config(sim, panel, status, message)
    type(sim_config), intent(in) :: sim
    type(panel_kernel_config), intent(in) :: panel
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=32) :: kernel_id, solver, boundary

    status = physics_config_ok
    message = ''
    kernel_id = lower_ascii(trim(panel%kernel_id))
    solver = lower_ascii(trim(sim%field_solver))
    boundary = lower_ascii(trim(sim%field_bc_mode))
    select case (trim(solver))
    case ('direct')
      if (trim(kernel_id) /= 'triangle_p0_exact_direct' .or. trim(boundary) /= 'free') then
        call reject(physics_config_unavailable, 'triangle_p0 direct requires its exact free-space kernel.', status, message)
      end if
    case ('treecode')
      if (trim(kernel_id) /= 'triangle_p0_exact_tree_near' .or. trim(boundary) /= 'free') then
        call reject( &
          physics_config_invalid_combination, &
          'triangle_p0 treecode requires exact panel-near/monopole-far free-space kernel.', status, message &
          )
      end if
    case ('fmm')
      if (trim(kernel_id) /= 'triangle_p0_exact_p2m_near') then
        call reject(physics_config_invalid_combination, 'triangle_p0 FMM requires exact P2M/panel-near kernel.', status, message)
      else if (trim(boundary) /= 'free' .and. trim(boundary) /= 'periodic2') then
        call reject(physics_config_invalid_combination, 'triangle_p0 FMM supports free or periodic2 boundaries.', status, message)
      end if
    case ('auto')
      if (trim(kernel_id) /= 'triangle_p0_exact_auto' .or. trim(boundary) /= 'free') then
        call reject(physics_config_invalid_combination, 'triangle_p0 auto requires its free-space auto kernel.', status, message)
      end if
    case default
      call reject(physics_config_unavailable, 'triangle_p0 supports direct, treecode, FMM, or auto solving.', status, message)
    end select
  end subroutine validate_phase1_panel_config

  !> 有限画像、panel reference、cached k!=0 の periodic2 構成を検証する。
  subroutine validate_active_physics_config(sim, field, periodic2, panel, status, message)
    type(sim_config), intent(in) :: sim
    type(field_physics_config), intent(in) :: field
    type(periodic2_physics_config), intent(in) :: periodic2
    type(panel_kernel_config), intent(in) :: panel
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=32) :: nonzero_backend

    status = physics_config_ok
    message = ''
    nonzero_backend = lower_ascii(trim(periodic2%nonzero_mode_backend))

    if (trim(nonzero_backend) == 'not_applicable') then
      call validate_phase1_panel_config(sim, panel, status, message)
      return
    end if

    if (trim(nonzero_backend) == 'legacy_finite_images') then
      if (trim(lower_ascii(periodic2%zero_mode_policy)) /= 'legacy_not_decomposed') then
        call reject( &
          physics_config_invalid_combination, &
          'legacy_finite_images requires zero_mode_policy=legacy_not_decomposed.', status, message &
          )
        return
      end if
      call validate_phase1_panel_config(sim, panel, status, message)
      return
    end if

    if (trim(nonzero_backend) == 'panel_spectral_reference') then
      if (trim(lower_ascii(field%backend)) /= 'direct' .or. trim(lower_ascii(sim%field_solver)) /= 'direct') then
        call reject( &
          physics_config_invalid_combination, &
          'panel_spectral_reference requires field backend direct.', status, message &
          )
        return
      end if
      if (trim(lower_ascii(sim%field_bc_mode)) /= 'periodic2' .or. .not. sim%use_box) then
        call reject( &
          physics_config_invalid_combination, &
          'panel_spectral_reference requires periodic2 box geometry.', status, message &
          )
        return
      end if
      if (any(sim%bc_low(1:2) /= bc_periodic) .or. any(sim%bc_high(1:2) /= bc_periodic) .or. &
          sim%bc_low(3) == bc_periodic .or. sim%bc_high(3) == bc_periodic) then
        call reject( &
          physics_config_unavailable, &
          'panel_spectral_reference requires x/y periodic and z nonperiodic.', status, message &
          )
        return
      end if
      if (trim(lower_ascii(periodic2%zero_mode_policy)) /= 'exclude_k0' .or. &
          .not. supported_lower_boundary(periodic2%lower_boundary_model)) then
        call reject( &
          physics_config_invalid_combination, &
          'panel_spectral_reference requires exclude_k0 and a supported lower boundary model.', &
          status, message &
          )
        return
      end if
      if (trim(lower_ascii(panel%kernel_id)) /= 'triangle_p0_exact_direct') then
        call reject( &
          physics_config_invalid_combination, &
          'panel_spectral_reference requires triangle_p0_exact_direct.', status, message &
          )
        return
      end if
      if (periodic2%reference_mode_layers < 1_i32 .or. periodic2%panel_quadrature_order < 2_i32) then
        call reject( &
          physics_config_invalid_combination, &
          'Periodic reference orders are out of range.', status, message &
          )
        return
      end if
      if (.not. ieee_is_finite(periodic2%max_nonzero_mode_potential_step) .or. &
          periodic2%max_nonzero_mode_potential_step /= 0.0_dp) then
        call reject( &
          physics_config_invalid_combination, &
          'panel_spectral_reference does not support max_nonzero_mode_potential_step.', status, message &
          )
      end if
      return
    end if

    if (trim(nonzero_backend) /= 'cached_kneq0') then
      call reject(physics_config_invalid_combination, 'Unknown periodic nonzero-mode backend.', status, message)
      return
    end if
    if (trim(lower_ascii(field%backend)) /= 'fmm' .or. trim(lower_ascii(sim%field_solver)) /= 'fmm') then
      call reject(physics_config_invalid_combination, 'cached_kneq0 requires the periodic2 FMM backend.', status, message)
      return
    end if
    if (trim(lower_ascii(sim%field_bc_mode)) /= 'periodic2' .or. .not. sim%use_box) then
      call reject(physics_config_invalid_combination, 'cached_kneq0 requires periodic2 box geometry.', status, message)
      return
    end if
    if (any(sim%bc_low(1:2) /= bc_periodic) .or. any(sim%bc_high(1:2) /= bc_periodic) .or. &
        sim%bc_low(3) == bc_periodic .or. sim%bc_high(3) == bc_periodic) then
      call reject(physics_config_unavailable, 'cached_kneq0 requires x/y periodic and z nonperiodic.', status, message)
      return
    end if
    if (trim(lower_ascii(periodic2%zero_mode_policy)) /= 'exclude_k0' .or. &
        .not. supported_lower_boundary(periodic2%lower_boundary_model)) then
      call reject( &
        physics_config_invalid_combination, &
        'cached_kneq0 requires exclude_k0 and a supported lower boundary model.', status, message &
        )
      return
    end if
    if (trim(lower_ascii(panel%kernel_id)) /= 'triangle_p0_exact_p2m_near') then
      call reject(physics_config_invalid_combination, 'cached_kneq0 requires the exact FMM panel kernel.', status, message)
      return
    end if
    if (.not. ieee_is_finite(periodic2%max_nonzero_mode_potential_step) .or. &
        periodic2%max_nonzero_mode_potential_step < 0.0_dp) then
      call reject(physics_config_invalid_combination, 'max_nonzero_mode_potential_step must be finite and nonnegative.', &
                  status, message)
      return
    end if
  end subroutine validate_active_physics_config

  pure logical function supported_lower_boundary(model) result(supported)
    character(len=*), intent(in) :: model
    character(len=32) :: normalized

    normalized = lower_ascii(trim(model))
    supported = trim(normalized) == 'e_bottom_zero' .or. trim(normalized) == 'symmetric_vacuum'
  end function supported_lower_boundary

  pure subroutine reject(code, text, status, message)
    integer(i32), intent(in) :: code
    character(len=*), intent(in) :: text
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = code
    message = text
  end subroutine reject

end module bem_physics_config_types

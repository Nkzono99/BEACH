!> 場・periodic2・panel・外部プラズマ・coupling の型付き設定と互換正規化を定義する。
module bem_physics_config_types
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_config, bc_periodic, bc_open
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
    integer(i32) :: interface_sample_n = 5_i32
    real(dp) :: interface_phi_tolerance = 1.0e-3_dp
    real(dp) :: interface_field_tolerance = 1.0e-3_dp
  end type periodic2_physics_config

  type, public :: panel_kernel_config
    character(len=32) :: source_model = 'point'
    character(len=32) :: kernel_id = 'softened_point'
    character(len=32) :: surface_side_policy = 'not_applicable'
  end type panel_kernel_config

  type, public :: outer_plasma_config
    character(len=32) :: model = 'none'
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
  end type outer_plasma_config

  type, public :: coupling_config
    character(len=32) :: update_mode = 'explicit'
    character(len=32) :: particle_transfer_mode = 'none'
    integer(i32) :: outer_update_stride = 1_i32
    real(dp) :: field_evolution_timescale = 0.0_dp
    real(dp) :: max_frozen_field_ratio = 0.1_dp
    real(dp) :: outer_orbit_dt = 0.0_dp
    integer(i32) :: outer_orbit_max_steps = 100000_i32
    real(dp) :: outer_orbit_energy_tolerance = 1.0e-4_dp
    logical :: outer_queue_enabled = .false.
  end type coupling_config

  public :: normalize_legacy_physics_config
  public :: validate_phase0_physics_config
  public :: validate_phase1_panel_config
  public :: validate_active_physics_config

contains

  !> 現行 `sim_config` の意味を変えず、型付き物理設定へ写像する。
  subroutine normalize_legacy_physics_config(sim, field, periodic2, panel, outer, coupling)
    type(sim_config), intent(in) :: sim
    type(field_physics_config), intent(out) :: field
    type(periodic2_physics_config), intent(out) :: periodic2
    type(panel_kernel_config), intent(out) :: panel
    type(outer_plasma_config), intent(out) :: outer
    type(coupling_config), intent(out) :: coupling
    character(len=32) :: bc_mode, far_mode

    field = field_physics_config()
    periodic2 = periodic2_physics_config()
    panel = panel_kernel_config()
    outer = outer_plasma_config()
    coupling = coupling_config()

    field%backend = lower_ascii(trim(sim%field_solver))
    field%normalization = lower_ascii(trim(sim%field_normalization))
    bc_mode = lower_ascii(trim(sim%field_bc_mode))
    if (trim(bc_mode) /= 'periodic2') return

    far_mode = lower_ascii(trim(sim%field_periodic_far_correction))
    select case (trim(far_mode))
    case ('none', 'auto')
      periodic2%nonzero_mode_backend = 'legacy_finite_images'
      periodic2%zero_mode_policy = 'legacy_not_decomposed'
    case ('m2l_root_oracle')
      periodic2%nonzero_mode_backend = 'legacy_root_oracle'
      periodic2%zero_mode_policy = 'legacy_charged_walls'
    case ('cached_kneq0')
      periodic2%nonzero_mode_backend = 'cached_kneq0'
      periodic2%zero_mode_policy = 'exclude_k0'
    case default
      periodic2%nonzero_mode_backend = 'invalid'
      periodic2%zero_mode_policy = 'invalid'
    end select
    if (trim(far_mode) == 'cached_kneq0') then
      periodic2%lower_boundary_model = 'e_bottom_zero'
    else
      periodic2%lower_boundary_model = 'legacy_implicit'
    end if
  end subroutine normalize_legacy_physics_config

  !> Phase 0 で実装済みの組合せだけを許可し、将来 mode は fail closed にする。
  subroutine validate_phase0_physics_config(field, periodic2, panel, outer, coupling, status, message)
    type(field_physics_config), intent(in) :: field
    type(periodic2_physics_config), intent(in) :: periodic2
    type(panel_kernel_config), intent(in) :: panel
    type(outer_plasma_config), intent(in) :: outer
    type(coupling_config), intent(in) :: coupling
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=32) :: backend, normalization, source_model, nonzero_backend, zero_policy

    status = physics_config_ok
    message = ''
    backend = lower_ascii(trim(field%backend))
    normalization = lower_ascii(trim(field%normalization))
    source_model = lower_ascii(trim(panel%source_model))
    nonzero_backend = lower_ascii(trim(periodic2%nonzero_mode_backend))
    zero_policy = lower_ascii(trim(periodic2%zero_mode_policy))

    select case (trim(backend))
    case ('direct', 'treecode', 'fmm', 'auto')
      continue
    case default
      call reject(physics_config_invalid_combination, 'Unknown field backend.', status, message)
      return
    end select
    select case (trim(normalization))
    case ('si', 'box', 'mesh', 'length')
      continue
    case default
      call reject(physics_config_invalid_combination, 'Unknown field normalization.', status, message)
      return
    end select

    if (trim(source_model) /= 'point') then
      call reject(physics_config_unavailable, 'Only point sources are available in Phase 0.', status, message)
      return
    end if
    if (trim(lower_ascii(outer%model)) /= 'none' .or. &
        trim(lower_ascii(coupling%particle_transfer_mode)) /= 'none') then
      call reject(physics_config_unavailable, 'Outer plasma coupling is not available in Phase 0.', status, message)
      return
    end if
    if (coupling%outer_update_stride < 1_i32) then
      call reject(physics_config_invalid_combination, 'coupling.outer_update_stride must be >= 1.', status, message)
      return
    end if

    select case (trim(nonzero_backend))
    case ('not_applicable')
      if (trim(zero_policy) /= 'not_applicable') then
        call reject(physics_config_invalid_combination, 'Free-space mode cannot select a zero-mode policy.', status, message)
      end if
    case ('legacy_finite_images')
      if (trim(zero_policy) /= 'legacy_not_decomposed') then
        call reject( &
          physics_config_invalid_combination, &
          'legacy_finite_images requires zero_mode_policy=legacy_not_decomposed.', status, message &
          )
      end if
    case ('legacy_root_oracle')
      if (trim(zero_policy) /= 'legacy_charged_walls') then
        call reject( &
          physics_config_invalid_combination, &
          'legacy_root_oracle requires zero_mode_policy=legacy_charged_walls.', status, message &
          )
      end if
    case ('reference_kneq0')
      if (trim(zero_policy) /= 'exclude_k0') then
        call reject(physics_config_invalid_combination, 'kneq0 backends require zero_mode_policy=exclude_k0.', status, message)
      else
        call reject(physics_config_unavailable, 'Separated k0/kneq0 physics is not available in Phase 0.', status, message)
      end if
    case ('cached_kneq0')
      if (trim(zero_policy) /= 'exclude_k0') then
        call reject(physics_config_invalid_combination, 'cached_kneq0 requires zero_mode_policy=exclude_k0.', status, message)
      end if
    case default
      call reject(physics_config_invalid_combination, 'Unknown periodic2 nonzero-mode backend.', status, message)
    end select
  end subroutine validate_phase0_physics_config

  !> triangle_p0 direct/treecode/FMM kernel の solver/boundary/softening 契約を検証する。
  subroutine validate_phase1_panel_config(sim, panel, status, message)
    type(sim_config), intent(in) :: sim
    type(panel_kernel_config), intent(in) :: panel
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=32) :: source_model, kernel_id, solver, boundary

    status = physics_config_ok
    message = ''
    source_model = lower_ascii(trim(panel%source_model))
    kernel_id = lower_ascii(trim(panel%kernel_id))
    solver = lower_ascii(trim(sim%field_solver))
    boundary = lower_ascii(trim(sim%field_bc_mode))
    if (trim(source_model) == 'point') return
    if (trim(source_model) /= 'triangle_p0') then
      call reject(physics_config_invalid_combination, 'Unknown panel source model.', status, message)
      return
    end if
    if (sim%softening /= 0.0_dp) then
      call reject(physics_config_invalid_combination, 'triangle_p0 requires softening=0.', status, message)
      return
    end if
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
      else if (trim(lower_ascii(sim%field_periodic_far_correction)) == 'm2l_root_oracle') then
        call reject(physics_config_unavailable, 'triangle_p0 FMM does not support the point-source root oracle.', status, message)
      end if
    case ('auto')
      if (trim(kernel_id) /= 'triangle_p0_exact_auto' .or. trim(boundary) /= 'free') then
        call reject(physics_config_invalid_combination, 'triangle_p0 auto requires its free-space auto kernel.', status, message)
      end if
    case default
      call reject(physics_config_unavailable, 'triangle_p0 supports direct, treecode, FMM, or auto solving.', status, message)
    end select
  end subroutine validate_phase1_panel_config

  subroutine validate_active_physics_config(sim, field, periodic2, panel, outer, coupling, status, message)
    type(sim_config), intent(in) :: sim
    type(field_physics_config), intent(in) :: field
    type(periodic2_physics_config), intent(in) :: periodic2
    type(panel_kernel_config), intent(in) :: panel
    type(outer_plasma_config), intent(in) :: outer
    type(coupling_config), intent(in) :: coupling
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=32) :: nonzero_backend

    nonzero_backend = lower_ascii(trim(periodic2%nonzero_mode_backend))
    if (trim(nonzero_backend) == 'cached_kneq0') then
      call validate_cached_periodic_config(sim, field, periodic2, panel, outer, coupling, status, message)
      return
    end if
    if (trim(nonzero_backend) /= 'panel_spectral_reference') then
      if (trim(lower_ascii(panel%source_model)) == 'triangle_p0') then
        call validate_phase1_panel_config(sim, panel, status, message)
        if (status /= physics_config_ok) return
        if (trim(lower_ascii(outer%model)) /= 'none') then
          call reject(physics_config_unavailable, 'Free-space triangle_p0 does not support outer plasma.', status, message)
        end if
      else
        call validate_phase0_physics_config(field, periodic2, panel, outer, coupling, status, message)
      end if
      return
    end if

    status = physics_config_ok
    message = ''
    if (trim(lower_ascii(field%backend)) /= 'direct' .or. trim(lower_ascii(sim%field_solver)) /= 'direct') then
      call reject( &
        physics_config_invalid_combination, 'panel_spectral_reference requires field backend direct.', status, message &
        )
      return
    end if
    if (trim(lower_ascii(sim%field_bc_mode)) /= 'periodic2' .or. .not. sim%use_box) then
      call reject( &
        physics_config_invalid_combination, 'panel_spectral_reference requires periodic2 box geometry.', status, message &
        )
      return
    end if
    if (any(sim%bc_low(1:2) /= bc_periodic) .or. any(sim%bc_high(1:2) /= bc_periodic) .or. &
        sim%bc_low(3) == bc_periodic .or. sim%bc_high(3) == bc_periodic) then
      call reject( &
        physics_config_unavailable, 'Phase 2 split periodic model requires x/y periodic and z nonperiodic.', status, message &
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
    if (trim(lower_ascii(panel%source_model)) /= 'triangle_p0' .or. &
        trim(lower_ascii(panel%kernel_id)) /= 'triangle_p0_exact_direct' .or. sim%softening /= 0.0_dp) then
      call reject( &
        physics_config_invalid_combination, &
        'panel_spectral_reference requires triangle_p0_exact_direct and softening=0.', status, message &
        )
      return
    end if
    if (periodic2%reference_mode_layers < 1_i32 .or. periodic2%panel_quadrature_order < 2_i32 .or. &
        periodic2%interface_sample_n < 2_i32) then
      call reject(physics_config_invalid_combination, 'Periodic reference orders are out of range.', status, message)
      return
    end if
    if (periodic2%interface_phi_tolerance <= 0.0_dp .or. periodic2%interface_field_tolerance <= 0.0_dp) then
      call reject(physics_config_invalid_combination, 'Periodic interface tolerances must be positive.', status, message)
      return
    end if
    if ((trim(lower_ascii(outer%model)) /= 'linear_debye' .and. &
         trim(lower_ascii(outer%model)) /= 'kinetic_1d' .and. &
         trim(lower_ascii(outer%model)) /= 'unified_linear_response') .or. outer%debye_length <= 0.0_dp .or. &
        outer%thermal_voltage <= 0.0_dp .or. outer%max_linearity_ratio <= 0.0_dp .or. &
        outer%max_gap_ratio <= 0.0_dp .or. outer%max_local_charge_ratio <= 0.0_dp) then
      call reject(physics_config_invalid_combination, 'Split periodic mode requires a valid outer-plasma model.', status, message)
      return
    end if
    if (abs(outer%interface_z - sim%box_max(3)) > &
        64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(sim%box_max(3)))) then
      call reject( &
        physics_config_invalid_combination, 'Phase 2 interface_z must equal the z-high box face.', status, message &
        )
      return
    end if
    if (trim(lower_ascii(coupling%update_mode)) /= 'explicit' .or. coupling%outer_update_stride < 1_i32) then
      call reject( &
        physics_config_unavailable, 'Split periodic coupling requires explicit updates.', status, message &
        )
      return
    end if
    select case (trim(lower_ascii(coupling%particle_transfer_mode)))
    case ('none')
      if (trim(lower_ascii(outer%return_model)) /= 'none') then
        call reject(physics_config_invalid_combination, 'A return model requires particle transfer.', status, message)
        return
      end if
    case ('electrostatic_1d_instant_return')
      if (coupling%field_evolution_timescale <= 0.0_dp .or. coupling%max_frozen_field_ratio <= 0.0_dp) then
        call reject(physics_config_invalid_combination, 'Invalid electrostatic 1D instant-return coupling.', status, message)
        return
      end if
      if ((trim(lower_ascii(outer%model)) == 'linear_debye' .and. &
           trim(lower_ascii(outer%return_model)) /= 'electrostatic_1d_instant_return') .or. &
          (trim(lower_ascii(outer%model)) == 'kinetic_1d' .and. &
           trim(lower_ascii(outer%return_model)) /= 'kinetic_1d_profile_return') .or. &
          (trim(lower_ascii(outer%model)) /= 'linear_debye' .and. &
           trim(lower_ascii(outer%model)) /= 'kinetic_1d')) then
        call reject(physics_config_invalid_combination, &
                    'The 1D return model identifier must match the outer-plasma model.', status, message)
        return
      end if
      if (sim%bc_high(3) /= bc_open .or. any(sim%b0 /= 0.0_dp)) then
        call reject(physics_config_unavailable, 'Instant return requires an open z-high face and b0=0.', status, message)
        return
      end if
      if (trim(lower_ascii(outer%model)) == 'kinetic_1d' .and. &
          (trim(lower_ascii(sim%reservoir_potential_model)) /= 'none' .or. &
           trim(lower_ascii(sim%sheath_injection_model)) /= 'none')) then
        call reject(physics_config_invalid_combination, &
                    'Kinetic profile return cannot mix with legacy or Zhao injection corrections.', status, message)
        return
      end if
      if (coupling%outer_queue_enabled) then
        call reject(physics_config_unavailable, 'Persistent outer queue is not available yet.', status, message)
        return
      end if
    case ('electrostatic_3d_explicit_orbit')
      call validate_explicit_3d_orbit_config(sim, outer, coupling, status, message)
      if (status /= physics_config_ok) return
    case default
      call reject(physics_config_invalid_combination, 'Unknown coupling particle-transfer mode.', status, message)
      return
    end select
    select case (trim(lower_ascii(outer%photoelectron_closure)))
    case ('none')
      continue
    case ('individual_return')
      if (trim(lower_ascii(coupling%particle_transfer_mode)) /= 'electrostatic_1d_instant_return' .or. &
          outer%photoelectron_histogram_bins < 1_i32 .or. outer%photoelectron_histogram_energy_max <= 0.0_dp .or. &
          outer%photoelectron_ambient_charge_scale <= 0.0_dp .or. outer%max_photoelectron_charge_ratio <= 0.0_dp) then
        call reject(physics_config_invalid_combination, 'Invalid individual photoelectron closure.', status, message)
        return
      end if
    case ('kinetic_mean')
      if (trim(lower_ascii(outer%model)) /= 'kinetic_1d' .or. &
          ((trim(lower_ascii(coupling%particle_transfer_mode)) == 'none' .and. &
            trim(lower_ascii(outer%return_model)) /= 'none') .or. &
           (trim(lower_ascii(coupling%particle_transfer_mode)) == 'electrostatic_1d_instant_return' .and. &
            trim(lower_ascii(outer%return_model)) /= 'kinetic_1d_profile_return') .or. &
           (trim(lower_ascii(coupling%particle_transfer_mode)) /= 'none' .and. &
            trim(lower_ascii(coupling%particle_transfer_mode)) /= 'electrostatic_1d_instant_return'))) then
        call reject(physics_config_invalid_combination, &
                    'kinetic_mean requires kinetic_1d with no transfer or kinetic profile return.', status, message)
        return
      end if
    case ('statistical_return')
      call reject(physics_config_unavailable, 'Statistical photoelectron return is not specified.', status, message)
      return
    case default
      call reject(physics_config_invalid_combination, 'Unknown photoelectron closure.', status, message)
      return
    end select
    if (trim(lower_ascii(outer%model)) == 'kinetic_1d' .and. &
        trim(lower_ascii(outer%photoelectron_closure)) == 'individual_return') then
      call reject(physics_config_unavailable, 'kinetic_1d cannot combine mean density with individual return.', status, message)
      return
    end if
    if (trim(lower_ascii(outer%model)) == 'unified_linear_response' .and. &
        trim(lower_ascii(outer%photoelectron_closure)) /= 'none') then
      call reject(physics_config_unavailable, &
                  'unified_linear_response does not support a photoelectron mean closure.', &
                  status, message)
      return
    end if
    if (trim(lower_ascii(outer%model)) == 'unified_linear_response' .and. &
        (outer%unified_grid_points < 17_i32 .or. outer%accessible_fraction_tolerance <= 0.0_dp)) then
      call reject(physics_config_invalid_combination, &
                  'unified_linear_response requires unified_grid_points >= 17 and a positive accessibility tolerance.', &
                  status, message)
      return
    end if
    if (any(sim%e0 /= 0.0_dp)) then
      call reject( &
        physics_config_unavailable, 'Phase 2 linear outer model currently requires sim.e0=0.', status, message &
        )
    end if
  end subroutine validate_active_physics_config

  subroutine validate_cached_periodic_config(sim, field, periodic2, panel, outer, coupling, status, message)
    type(sim_config), intent(in) :: sim
    type(field_physics_config), intent(in) :: field
    type(periodic2_physics_config), intent(in) :: periodic2
    type(panel_kernel_config), intent(in) :: panel
    type(outer_plasma_config), intent(in) :: outer
    type(coupling_config), intent(in) :: coupling
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = physics_config_ok
    message = ''
    if (trim(lower_ascii(field%backend)) /= 'fmm' .or. trim(lower_ascii(sim%field_solver)) /= 'fmm' .or. &
        trim(lower_ascii(sim%field_bc_mode)) /= 'periodic2' .or. &
        trim(lower_ascii(sim%field_periodic_far_correction)) /= 'cached_kneq0') then
      call reject(physics_config_invalid_combination, 'cached_kneq0 requires the periodic2 FMM backend.', status, message)
      return
    end if
    if (.not. sim%use_box .or. any(sim%bc_low(1:2) /= bc_periodic) .or. &
        any(sim%bc_high(1:2) /= bc_periodic) .or. sim%bc_low(3) == bc_periodic .or. &
        sim%bc_high(3) == bc_periodic) then
      call reject(physics_config_unavailable, 'cached_kneq0 requires x/y periodic and z nonperiodic.', status, message)
      return
    end if
    if (trim(lower_ascii(periodic2%zero_mode_policy)) /= 'exclude_k0' .or. &
        .not. supported_lower_boundary(periodic2%lower_boundary_model)) then
      call reject( &
        physics_config_invalid_combination, &
        'cached_kneq0 requires exclude_k0 and a supported lower boundary model.', &
        status, message &
        )
      return
    end if
    select case (trim(lower_ascii(coupling%particle_transfer_mode)))
    case ('none')
      if (trim(lower_ascii(outer%return_model)) /= 'none') then
        call reject(physics_config_invalid_combination, 'A cached return model requires particle transfer.', status, message)
        return
      end if
    case ('electrostatic_1d_instant_return')
      if (trim(lower_ascii(outer%model)) /= 'kinetic_1d' .or. &
          trim(lower_ascii(outer%return_model)) /= 'kinetic_1d_profile_return' .or. &
          coupling%field_evolution_timescale <= 0.0_dp .or. coupling%max_frozen_field_ratio <= 0.0_dp) then
        call reject(physics_config_invalid_combination, 'Invalid cached kinetic 1D profile return.', status, message)
        return
      end if
      if (sim%bc_high(3) /= bc_open .or. any(sim%b0 /= 0.0_dp)) then
        call reject(physics_config_unavailable, 'Kinetic profile return requires an open z-high face and b0=0.', &
                    status, message)
        return
      end if
      if (trim(lower_ascii(sim%reservoir_potential_model)) /= 'none' .or. &
          trim(lower_ascii(sim%sheath_injection_model)) /= 'none') then
        call reject(physics_config_invalid_combination, &
                    'Kinetic profile return cannot mix with legacy or Zhao injection corrections.', status, message)
        return
      end if
      if (coupling%outer_queue_enabled) then
        call reject(physics_config_unavailable, 'Persistent outer queue is not available yet.', status, message)
        return
      end if
    case ('electrostatic_3d_explicit_orbit')
      call validate_explicit_3d_orbit_config(sim, outer, coupling, status, message)
      if (status /= physics_config_ok) return
    case default
      call reject(physics_config_unavailable, 'cached_kneq0 does not support the requested particle transfer.', &
                  status, message)
      return
    end select
    select case (trim(lower_ascii(outer%model)))
    case ('none')
      continue
    case ('kinetic_1d')
      if (outer%debye_length <= 0.0_dp .or. outer%thermal_voltage <= 0.0_dp .or. &
          abs(outer%interface_z - sim%box_max(3)) > &
          64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(sim%box_max(3)))) then
        call reject(physics_config_invalid_combination, &
                    'kinetic_1d requires positive scales and interface_z at z-high.', status, message)
        return
      end if
      if (abs(outer%infinity_potential) > 64.0_dp*epsilon(1.0_dp)) then
        call reject(physics_config_invalid_combination, &
                    'kinetic_1d fixes the infinity-potential gauge to zero.', status, message)
        return
      end if
      if (trim(lower_ascii(outer%photoelectron_closure)) /= 'none' .and. &
          trim(lower_ascii(outer%photoelectron_closure)) /= 'kinetic_mean') then
        call reject(physics_config_unavailable, 'cached kinetic_1d supports none or kinetic_mean photoelectron closure.', &
                    status, message)
        return
      end if
    case ('unified_linear_response')
      if (outer%debye_length <= 0.0_dp .or. outer%thermal_voltage <= 0.0_dp .or. &
          outer%max_linearity_ratio <= 0.0_dp .or. outer%unified_grid_points < 17_i32 .or. &
          outer%accessible_fraction_tolerance <= 0.0_dp .or. &
          abs(outer%interface_z - sim%box_max(3)) > &
          64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(sim%box_max(3)))) then
        call reject(physics_config_invalid_combination, &
                    'unified_linear_response requires positive scales and interface_z at z-high.', status, message)
        return
      end if
      if (any(sim%e0 /= 0.0_dp)) then
        call reject(physics_config_unavailable, &
                    'unified_linear_response requires a finite far-potential gauge and sim.e0=0.', &
                    status, message)
        return
      end if
      if (trim(lower_ascii(panel%source_model)) /= 'triangle_p0' .or. &
          trim(lower_ascii(outer%photoelectron_closure)) /= 'none' .or. &
          (trim(lower_ascii(outer%return_model)) /= 'none' .and. &
           trim(lower_ascii(outer%return_model)) /= 'electrostatic_3d_explicit_orbit')) then
        call reject(physics_config_unavailable, &
                    'cached unified_linear_response received an unsupported source or closure.', &
                    status, message)
        return
      end if
    case default
      call reject(physics_config_unavailable, 'cached_kneq0 received an unsupported outer-plasma model.', status, message)
      return
    end select
    if (trim(lower_ascii(panel%source_model)) == 'triangle_p0') then
      call validate_phase1_panel_config(sim, panel, status, message)
    else if (trim(lower_ascii(panel%source_model)) /= 'point') then
      call reject(physics_config_invalid_combination, 'cached_kneq0 received an unknown source model.', status, message)
    end if
  end subroutine validate_cached_periodic_config

  pure logical function supported_lower_boundary(model) result(supported)
    character(len=*), intent(in) :: model

    select case (trim(lower_ascii(model)))
    case ('e_bottom_zero', 'symmetric_vacuum')
      supported = .true.
    case default
      supported = .false.
    end select
  end function supported_lower_boundary

  subroutine validate_explicit_3d_orbit_config(sim, outer, coupling, status, message)
    type(sim_config), intent(in) :: sim
    type(outer_plasma_config), intent(in) :: outer
    type(coupling_config), intent(in) :: coupling
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = physics_config_ok
    message = ''
    if (trim(lower_ascii(outer%model)) /= 'unified_linear_response' .or. &
        trim(lower_ascii(outer%return_model)) /= 'electrostatic_3d_explicit_orbit' .or. &
        coupling%field_evolution_timescale <= 0.0_dp .or. coupling%max_frozen_field_ratio <= 0.0_dp .or. &
        coupling%outer_orbit_dt <= 0.0_dp .or. coupling%outer_orbit_max_steps < 1_i32 .or. &
        coupling%outer_orbit_energy_tolerance <= 0.0_dp) then
      call reject(physics_config_invalid_combination, 'Invalid electrostatic 3D explicit-orbit coupling.', status, message)
      return
    end if
    if (sim%bc_high(3) /= bc_open .or. any(sim%b0 /= 0.0_dp)) then
      call reject(physics_config_unavailable, 'Explicit 3D outer orbit requires an open z-high face and b0=0.', &
                  status, message)
      return
    end if
    if (coupling%outer_queue_enabled) then
      call reject(physics_config_unavailable, 'Persistent outer queue is not available yet.', status, message)
    end if
  end subroutine validate_explicit_3d_orbit_config

  pure subroutine reject(code, text, status, message)
    integer(i32), intent(in) :: code
    character(len=*), intent(in) :: text
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = code
    message = text
  end subroutine reject

end module bem_physics_config_types

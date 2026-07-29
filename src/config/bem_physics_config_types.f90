!> 場・periodic2・panel・外部プラズマ・coupling の型付き設定と互換正規化を定義する。
module bem_physics_config_types
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_types, only: sim_config, bc_periodic, bc_open
  use bem_string_utils, only: lower_ascii
  use bem_external_boundary_contract, only: &
    external_boundary_contract_type, external_boundary_ok, resolve_external_boundary_contract
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
    real(dp) :: max_nonzero_mode_potential_step = 0.0_dp
  end type periodic2_physics_config

  type, public :: panel_kernel_config
    character(len=32) :: kernel_id = 'triangle_p0_exact_auto'
    character(len=32) :: surface_side_policy = 'per_element'
  end type panel_kernel_config

  type, public :: outer_plasma_config
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
  end type outer_plasma_config

  type, public :: coupling_config
    character(len=32) :: update_mode = 'explicit'
    character(len=32) :: particle_transfer_mode = 'none'
    character(len=32) :: steady_start_mode = 'none'
    integer(i32) :: steady_start_mesh_id = 1_i32
    integer(i32) :: outer_update_stride = 1_i32
    real(dp) :: field_evolution_timescale = 0.0_dp
    real(dp) :: max_frozen_field_ratio = 0.1_dp
    logical :: outer_queue_enabled = .false.
  end type coupling_config

  public :: normalize_legacy_physics_config
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

  !> triangle_p0 direct/treecode/FMM kernel の solver/boundary 契約を検証する。
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
    type(external_boundary_contract_type) :: boundary_contract
    integer(i32) :: boundary_status

    select case (trim(lower_ascii(coupling%update_mode)))
    case ('explicit')
      continue
    case ('implicit_mean')
      if (trim(lower_ascii(outer%model)) /= 'kinetic_1d' .or. &
          (trim(lower_ascii(outer%kinetic_closure)) /= 'ambient_linear_debye' .and. &
           trim(lower_ascii(outer%kinetic_closure)) /= 'zhao_charge_driven') .or. &
          trim(lower_ascii(outer%photoelectron_density_model)) /= 'none' .or. &
          trim(lower_ascii(coupling%particle_transfer_mode)) /= 'electrostatic_1d_instant_return' .or. &
          coupling%outer_queue_enabled .or. coupling%outer_update_stride /= 1_i32 .or. &
          .not. ieee_is_finite(sim%batch_duration) .or. sim%batch_duration <= 0.0_dp) then
        call reject( &
          physics_config_invalid_combination, &
          'implicit_mean requires ambient-linear or Zhao kinetic_1d same-batch transfer, '// &
          'photoelectron_density_model=none, outer_update_stride=1, and positive batch_duration.', &
          status, message &
          )
        return
      end if
      if (trim(lower_ascii(outer%kinetic_closure)) == 'zhao_charge_driven' .and. &
          trim(lower_ascii(coupling%steady_start_mode)) /= 'zhao_floating') then
        call reject( &
          physics_config_invalid_combination, &
          'implicit Zhao mean coupling requires a zhao_floating branch anchor.', status, message &
          )
        return
      end if
    case default
      call reject(physics_config_invalid_combination, 'Unknown coupling.update_mode.', status, message)
      return
    end select
    if (coupling%outer_update_stride < 1_i32) then
      call reject(physics_config_invalid_combination, 'coupling.outer_update_stride must be >= 1.', status, message)
      return
    end if
    call resolve_external_boundary_contract( &
      sim%reservoir_potential_model, sim%open_boundary_model, outer%model, &
      outer%kinetic_closure, outer%return_model, coupling%particle_transfer_mode, coupling%outer_queue_enabled, &
      boundary_contract, boundary_status, message &
      )
    if (boundary_status /= external_boundary_ok) then
      status = physics_config_invalid_combination
      return
    end if
    call validate_kinetic_closure_config(outer, status, message)
    if (status /= physics_config_ok) return
    call validate_zhao_charge_driven_sim_config(sim, outer, status, message)
    if (status /= physics_config_ok) return
    call validate_photoelectron_config(outer, status, message)
    if (status /= physics_config_ok) return
    call validate_outer_queue_config(sim, outer, coupling, status, message)
    if (status /= physics_config_ok) return
    call validate_steady_start_config(periodic2, outer, coupling, status, message)
    if (status /= physics_config_ok) return
    nonzero_backend = lower_ascii(trim(periodic2%nonzero_mode_backend))
    if (trim(nonzero_backend) == 'cached_kneq0') then
      call validate_cached_periodic_config(sim, field, periodic2, panel, outer, coupling, status, message)
      return
    end if
    if (trim(nonzero_backend) == 'not_applicable') then
      if (trim(lower_ascii(periodic2%zero_mode_policy)) /= 'not_applicable') then
        call reject( &
          physics_config_invalid_combination, &
          'Free-space mode requires zero_mode_policy=not_applicable.', status, message &
          )
        return
      end if
      call validate_phase1_panel_config(sim, panel, status, message)
      if (status /= physics_config_ok) return
      if (trim(lower_ascii(outer%model)) /= 'none') then
        call reject(physics_config_unavailable, 'Free-space triangle_p0 does not support outer plasma.', status, message)
      end if
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
      if (status /= physics_config_ok) return
      if (trim(lower_ascii(outer%model)) /= 'none') then
        call reject(physics_config_unavailable, 'Free-space triangle_p0 does not support outer plasma.', status, message)
      end if
      return
    end if
    if (trim(nonzero_backend) /= 'panel_spectral_reference') then
      call reject(physics_config_invalid_combination, 'Unknown periodic nonzero-mode backend.', status, message)
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
    if (trim(lower_ascii(panel%kernel_id)) /= 'triangle_p0_exact_direct') then
      call reject( &
        physics_config_invalid_combination, &
        'panel_spectral_reference requires triangle_p0_exact_direct.', status, message &
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
    if (trim(lower_ascii(outer%model)) /= 'kinetic_1d' .or. outer%debye_length <= 0.0_dp .or. &
        outer%thermal_voltage <= 0.0_dp .or. &
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
    if ((trim(lower_ascii(coupling%update_mode)) /= 'explicit' .and. &
         trim(lower_ascii(coupling%update_mode)) /= 'implicit_mean') .or. coupling%outer_update_stride < 1_i32) then
      call reject( &
        physics_config_unavailable, 'Split periodic coupling received an unsupported update mode.', status, message &
        )
      return
    end if
    select case (trim(lower_ascii(coupling%particle_transfer_mode)))
    case ('none')
      continue
    case ('electrostatic_1d_instant_return')
      if (.not. ieee_is_finite(coupling%field_evolution_timescale) .or. &
          .not. ieee_is_finite(coupling%max_frozen_field_ratio) .or. &
          coupling%field_evolution_timescale <= 0.0_dp .or. coupling%max_frozen_field_ratio <= 0.0_dp) then
        call reject(physics_config_invalid_combination, 'Invalid electrostatic 1D instant-return coupling.', status, message)
        return
      end if
      if (sim%bc_high(3) /= bc_open .or. any(sim%b0 /= 0.0_dp)) then
        call reject(physics_config_unavailable, 'Instant return requires an open z-high face and b0=0.', status, message)
        return
      end if
    case default
      call reject(physics_config_invalid_combination, 'Unknown coupling particle-transfer mode.', status, message)
      return
    end select
    if (any(sim%e0 /= 0.0_dp)) then
      call reject( &
        physics_config_unavailable, 'Periodic split outer models currently require sim.e0=0.', status, message &
        )
    end if
  end subroutine validate_active_physics_config

  subroutine validate_kinetic_closure_config(outer, status, message)
    type(outer_plasma_config), intent(in) :: outer
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=32) :: closure, branch

    status = physics_config_ok
    message = ''
    closure = lower_ascii(trim(outer%kinetic_closure))
    branch = lower_ascii(trim(outer%zhao_branch))
    if (.not. ieee_is_finite(outer%photoelectron_source_scale) .or. &
        outer%photoelectron_source_scale < 0.0_dp) then
      call reject(physics_config_invalid_combination, &
                  'outer_plasma.photoelectron_source_scale must be finite and non-negative.', status, message)
      return
    end if
    select case (trim(branch))
    case ('auto', 'a', 'b', 'c')
      continue
    case default
      call reject(physics_config_invalid_combination, 'Unknown outer_plasma.zhao_branch.', status, message)
      return
    end select
    select case (trim(closure))
    case ('absorbing_maxwellian', 'ambient_linear_debye')
      if (trim(branch) /= 'auto') then
        call reject(physics_config_invalid_combination, &
                    'outer_plasma.zhao_branch requires kinetic_closure=zhao_charge_driven.', status, message)
      else if (abs(outer%photoelectron_source_scale - 1.0_dp) > 64.0_dp*epsilon(1.0_dp)) then
        call reject(physics_config_invalid_combination, &
                    'outer_plasma.photoelectron_source_scale requires kinetic_closure=zhao_charge_driven.', &
                    status, message)
      end if
      if (status == physics_config_ok .and. trim(closure) == 'ambient_linear_debye') then
        if (trim(lower_ascii(outer%model)) /= 'kinetic_1d') then
          call reject(physics_config_invalid_combination, &
                      'ambient_linear_debye requires outer_plasma.model=kinetic_1d.', status, message)
        else if (trim(lower_ascii(outer%photoelectron_density_model)) /= 'none') then
          call reject(physics_config_invalid_combination, &
                      'ambient_linear_debye requires photoelectron_density_model=none.', status, message)
        end if
      end if
    case ('zhao_charge_driven')
      if (trim(lower_ascii(outer%model)) /= 'kinetic_1d') then
        call reject(physics_config_invalid_combination, &
                    'zhao_charge_driven requires outer_plasma.model=kinetic_1d.', status, message)
      else if (trim(lower_ascii(outer%photoelectron_density_model)) /= 'none') then
        call reject(physics_config_invalid_combination, &
                    'zhao_charge_driven includes its photoelectron population and requires photoelectron_density_model=none.', &
                    status, message)
      end if
    case default
      call reject(physics_config_invalid_combination, 'Unknown outer_plasma.kinetic_closure.', status, message)
    end select
  end subroutine validate_kinetic_closure_config

  subroutine validate_zhao_charge_driven_sim_config(sim, outer, status, message)
    type(sim_config), intent(in) :: sim
    type(outer_plasma_config), intent(in) :: outer
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = physics_config_ok
    message = ''
    if (trim(lower_ascii(outer%kinetic_closure)) /= 'zhao_charge_driven') return
    if (trim(lower_ascii(sim%reservoir_potential_model)) /= 'none') then
      call reject(physics_config_invalid_combination, &
                  'zhao_charge_driven cannot mix with reservoir-potential correction.', status, message)
      return
    end if
    if (outer%photoelectron_source_scale > 0.0_dp .and. &
        (.not. ieee_is_finite(sim%sheath_photoelectron_ref_density_cm3) .or. &
         sim%sheath_photoelectron_ref_density_cm3 <= 0.0_dp)) then
      call reject(physics_config_invalid_combination, &
                  'Photoemitting zhao_charge_driven requires a positive sheath_photoelectron_ref_density_cm3.', &
                  status, message)
    end if
  end subroutine validate_zhao_charge_driven_sim_config

  subroutine validate_photoelectron_config(outer, status, message)
    type(outer_plasma_config), intent(in) :: outer
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = physics_config_ok
    message = ''
    select case (trim(lower_ascii(outer%photoelectron_density_model)))
    case ('none')
      continue
    case ('kinetic_mean')
      if (trim(lower_ascii(outer%model)) /= 'kinetic_1d') then
        call reject( &
          physics_config_invalid_combination, &
          'photoelectron_density_model=kinetic_mean requires outer_plasma.model=kinetic_1d.', &
          status, message &
          )
        return
      end if
    case default
      call reject(physics_config_invalid_combination, 'Unknown photoelectron density model.', status, message)
      return
    end select
  end subroutine validate_photoelectron_config

  !> Persistent outer-flight events are currently a Zhao-only transient closure.
  subroutine validate_outer_queue_config(sim, outer, coupling, status, message)
    type(sim_config), intent(in) :: sim
    type(outer_plasma_config), intent(in) :: outer
    type(coupling_config), intent(in) :: coupling
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = physics_config_ok
    message = ''
    if (.not. coupling%outer_queue_enabled) return
    if (outer%photoelectron_source_scale <= 0.0_dp) then
      call reject(physics_config_invalid_combination, &
                  'Persistent outer queue requires a positive photoelectron source scale.', status, message)
      return
    end if
    if (trim(lower_ascii(outer%model)) /= 'kinetic_1d' .or. &
        trim(lower_ascii(outer%kinetic_closure)) /= 'zhao_charge_driven' .or. &
        trim(lower_ascii(outer%zhao_branch)) /= 'auto' .or. &
        trim(lower_ascii(outer%return_model)) /= 'kinetic_1d_profile_return' .or. &
        trim(lower_ascii(coupling%particle_transfer_mode)) /= 'electrostatic_1d_instant_return') then
      call reject(physics_config_invalid_combination, &
                  'Persistent outer queue requires the automatic Zhao kinetic_1d profile return.', status, message)
      return
    end if
    if (.not. ieee_is_finite(sim%batch_duration) .or. &
        .not. ieee_is_finite(coupling%field_evolution_timescale) .or. &
        .not. ieee_is_finite(coupling%max_frozen_field_ratio) .or. &
        sim%batch_duration <= 0.0_dp .or. coupling%field_evolution_timescale <= 0.0_dp .or. &
        coupling%max_frozen_field_ratio <= 0.0_dp) then
      call reject(physics_config_invalid_combination, &
                  'Persistent outer queue requires finite positive time-scale limits.', status, message)
      return
    end if
    if (sim%batch_duration > coupling%max_frozen_field_ratio*coupling%field_evolution_timescale) then
      call reject(physics_config_invalid_combination, &
                  'Persistent outer queue batch_duration does not resolve the configured field evolution timescale.', &
                  status, message)
      return
    end if
    if (coupling%outer_update_stride /= 1_i32) then
      call reject(physics_config_invalid_combination, &
                  'Persistent outer queue requires outer_update_stride=1.', status, message)
      return
    end if
  end subroutine validate_outer_queue_config

  !> Validate the opt-in pre-equilibrated surface-charge initialization contract.
  subroutine validate_steady_start_config(periodic2, outer, coupling, status, message)
    type(periodic2_physics_config), intent(in) :: periodic2
    type(outer_plasma_config), intent(in) :: outer
    type(coupling_config), intent(in) :: coupling
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=32) :: mode

    status = physics_config_ok
    message = ''
    mode = lower_ascii(trim(coupling%steady_start_mode))
    if (coupling%steady_start_mesh_id < 1_i32) then
      call reject(physics_config_invalid_combination, &
                  'coupling.steady_start_mesh_id must be >= 1.', status, message)
      return
    end if
    select case (trim(mode))
    case ('none')
      return
    case ('zhao_floating')
      continue
    case default
      call reject(physics_config_invalid_combination, 'Unknown coupling.steady_start_mode.', status, message)
      return
    end select
    if (trim(lower_ascii(outer%model)) /= 'kinetic_1d' .or. &
        trim(lower_ascii(outer%kinetic_closure)) /= 'zhao_charge_driven' .or. &
        trim(lower_ascii(outer%return_model)) /= 'kinetic_1d_profile_return' .or. &
        trim(lower_ascii(coupling%particle_transfer_mode)) /= 'electrostatic_1d_instant_return') then
      call reject(physics_config_invalid_combination, &
                  'zhao_floating steady start requires the Zhao kinetic_1d instant profile return.', &
                  status, message)
      return
    end if
    if (coupling%outer_queue_enabled) then
      call reject(physics_config_invalid_combination, &
                  'zhao_floating steady start requires coupling.outer_queue_enabled=false.', status, message)
      return
    end if
    if (trim(lower_ascii(periodic2%zero_mode_policy)) /= 'exclude_k0' .or. &
        .not. supported_lower_boundary(periodic2%lower_boundary_model)) then
      call reject(physics_config_invalid_combination, &
                  'zhao_floating steady start requires exclude_k0 and a supported lower boundary model.', &
                  status, message)
    end if
  end subroutine validate_steady_start_config

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
      continue
    case ('electrostatic_1d_instant_return')
      if (trim(lower_ascii(outer%model)) /= 'kinetic_1d' .or. &
          .not. ieee_is_finite(coupling%field_evolution_timescale) .or. &
          .not. ieee_is_finite(coupling%max_frozen_field_ratio) .or. &
          coupling%field_evolution_timescale <= 0.0_dp .or. coupling%max_frozen_field_ratio <= 0.0_dp) then
        call reject(physics_config_invalid_combination, 'Invalid cached kinetic 1D profile return.', status, message)
        return
      end if
      if (sim%bc_high(3) /= bc_open .or. any(sim%b0 /= 0.0_dp)) then
        call reject(physics_config_unavailable, 'Kinetic profile return requires an open z-high face and b0=0.', &
                    status, message)
        return
      end if
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
      if (trim(lower_ascii(outer%photoelectron_density_model)) /= 'none' .and. &
          trim(lower_ascii(outer%photoelectron_density_model)) /= 'kinetic_mean') then
        call reject(physics_config_unavailable, 'cached kinetic_1d supports none or kinetic_mean photoelectron density.', &
                    status, message)
        return
      end if
    case default
      call reject(physics_config_unavailable, 'cached_kneq0 received an unsupported outer-plasma model.', status, message)
      return
    end select
    call validate_phase1_panel_config(sim, panel, status, message)
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

  pure subroutine reject(code, text, status, message)
    integer(i32), intent(in) :: code
    character(len=*), intent(in) :: text
    integer(i32), intent(out) :: status
    character(len=*), intent(out) :: message

    status = code
    message = text
  end subroutine reject

end module bem_physics_config_types

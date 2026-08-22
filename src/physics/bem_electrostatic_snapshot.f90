!> 1バッチ内で固定する静電場 snapshot。
module bem_electrostatic_snapshot
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32
  use bem_constants, only: eps0
  use bem_types, only: mesh_type, sim_config, bc_periodic
  use bem_field_solver, only: field_solver_type
  use bem_physics_config_types, only: field_physics_config, periodic2_physics_config, panel_kernel_config
  use bem_string_utils, only: lower_ascii
  use bem_periodic_zero_mode_plan, only: periodic_zero_mode_plan_type, periodic_zero_mode_state_type, &
                                         build_periodic_zero_mode_plan, refresh_periodic_zero_mode_state, &
                                         symmetric_vacuum_bottom_field, periodic_zero_mode_ok
  use bem_periodic_zero_mode_eval, only: eval_periodic_zero_mode, zero_mode_trace_plus
  use bem_coulomb_fmm_periodic_nonzero_reference, only: eval_periodic_nonzero_panel_reference
  use bem_panel_geometry, only: panel_geometry_type, init_panel_geometry, panel_geometry_ok
  use bem_panel_kernel, only: panel_potential_field, panel_side_principal_value
  implicit none
  private

  type, public :: electrostatic_diagnostics_type
    logical :: split_periodic_active = .false.
    logical :: top_reference_available = .false.
    integer(i32) :: top_reference_last_batch = -1_i32
    real(dp) :: top_reference_simulated_time = 0.0_dp
    real(dp) :: top_reference_z_high = 0.0_dp
    integer(i32) :: top_reference_sample_n = 0_i32
    real(dp) :: top_reference_potential_mean = 0.0_dp
    real(dp) :: top_reference_potential_std = 0.0_dp
    real(dp) :: top_reference_potential_min = 0.0_dp
    real(dp) :: top_reference_potential_max = 0.0_dp
    real(dp) :: gauss_residual = 0.0_dp
    character(len=64) :: status = 'not_applicable'
    logical :: periodic_cache_hit = .false.
    integer(i32) :: periodic_operator_build_count = 0_i32
    character(len=128) :: periodic_cache_fingerprint = ''
    character(len=512) :: periodic_cache_path = ''
  end type electrostatic_diagnostics_type

  type, public :: electrostatic_snapshot_type
    type(field_solver_type) :: nonzero_solver
    real(dp) :: prescribed_e(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    real(dp) :: prescribed_phi_origin(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    logical :: use_zero_mode = .false.
    logical :: use_panel_spectral_reference = .false.
    logical :: use_cached_kneq0 = .false.
    real(dp) :: periodic_length(2) = [0.0_dp, 0.0_dp]
    real(dp) :: periodic_origin(2) = [0.0_dp, 0.0_dp]
    integer(i32) :: reference_mode_layers = 0_i32
    integer(i32) :: panel_quadrature_order = 0_i32
    type(periodic_zero_mode_plan_type) :: zero_plan
    type(periodic_zero_mode_state_type) :: zero_state
    type(periodic2_physics_config) :: periodic_options
    type(electrostatic_diagnostics_type) :: diagnostics
    real(dp) :: gauss_residual = 0.0_dp
    logical :: matching_plane_gauge_active = .false.
    real(dp) :: matching_plane_z = 0.0_dp
    real(dp) :: matching_plane_phi = 0.0_dp
  contains
    procedure :: init => init_electrostatic_snapshot
    procedure :: refresh => refresh_electrostatic_snapshot
    procedure :: eval_local_e => eval_snapshot_local_e
    procedure :: eval_local_phi => eval_snapshot_local_phi
    procedure :: eval_local_phi_without_primary_self => eval_snapshot_local_phi_without_primary_self
    procedure :: compute_mesh_potential => compute_snapshot_mesh_potential
    procedure :: measure_kneq0_potential_step => measure_snapshot_kneq0_potential_step
    procedure :: get_diagnostics => get_snapshot_diagnostics
    procedure :: set_matching_plane_gauge => set_snapshot_matching_plane_gauge
    procedure :: clear_matching_plane_gauge => clear_snapshot_matching_plane_gauge
    procedure :: get_matching_plane_displacement => get_snapshot_matching_plane_displacement
  end type electrostatic_snapshot_type

  interface
    module subroutine eval_snapshot_local_e(self, mesh, position, electric_field)
      class(electrostatic_snapshot_type), intent(inout) :: self
      type(mesh_type), intent(in) :: mesh
      real(dp), intent(in) :: position(3)
      real(dp), intent(out) :: electric_field(3)
    end subroutine eval_snapshot_local_e

    module subroutine eval_snapshot_local_phi(self, mesh, sim, position, potential)
      class(electrostatic_snapshot_type), intent(inout) :: self
      type(mesh_type), intent(in) :: mesh
      type(sim_config), intent(in) :: sim
      real(dp), intent(in) :: position(3)
      real(dp), intent(out) :: potential
    end subroutine eval_snapshot_local_phi

    module subroutine eval_snapshot_local_phi_without_primary_self(self, mesh, sim, element, potential)
      class(electrostatic_snapshot_type), intent(inout) :: self
      type(mesh_type), intent(in) :: mesh
      type(sim_config), intent(in) :: sim
      integer(i32), intent(in) :: element
      real(dp), intent(out) :: potential
    end subroutine eval_snapshot_local_phi_without_primary_self

    module subroutine compute_snapshot_mesh_potential(self, mesh, sim, potential)
      class(electrostatic_snapshot_type), intent(inout) :: self
      type(mesh_type), intent(in) :: mesh
      type(sim_config), intent(in) :: sim
      real(dp), intent(out) :: potential(:)
    end subroutine compute_snapshot_mesh_potential

    module subroutine measure_snapshot_kneq0_potential_step( &
      self, mesh, candidate_charge, max_abs_delta_phi_v, delta_phi_v &
      )
      class(electrostatic_snapshot_type), intent(inout) :: self
      type(mesh_type), intent(in) :: mesh
      real(dp), intent(in) :: candidate_charge(:)
      real(dp), intent(out) :: max_abs_delta_phi_v
      real(dp), intent(out), optional :: delta_phi_v(:)
    end subroutine measure_snapshot_kneq0_potential_step
  end interface

contains

  subroutine init_electrostatic_snapshot(self, mesh, sim, field_config, periodic_config, panel_config)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(field_physics_config), intent(in), optional :: field_config
    type(periodic2_physics_config), intent(in), optional :: periodic_config
    type(panel_kernel_config), intent(in), optional :: panel_config

    self%prescribed_e = sim%e0
    self%prescribed_phi_origin = 0.0_dp
    self%use_zero_mode = .false.
    self%use_panel_spectral_reference = .false.
    self%use_cached_kneq0 = .false.
    self%gauss_residual = 0.0_dp
    self%matching_plane_gauge_active = .false.
    self%matching_plane_z = 0.0_dp
    self%matching_plane_phi = 0.0_dp
    self%diagnostics = electrostatic_diagnostics_type()

    if (present(periodic_config)) then
      if (trim(lower_ascii(periodic_config%nonzero_mode_backend)) == 'panel_spectral_reference') then
        call init_panel_spectral_snapshot(self, mesh, sim, periodic_config)
        return
      end if
      if (trim(lower_ascii(periodic_config%nonzero_mode_backend)) == 'cached_kneq0') then
        if (.not. present(field_config) .or. .not. present(panel_config)) then
          error stop 'cached_kneq0 requires typed field and panel configuration.'
        end if
        call init_cached_periodic_snapshot(self, mesh, sim, field_config, periodic_config, panel_config)
        return
      end if
    end if

    if (present(field_config) .and. present(periodic_config) .and. present(panel_config)) then
      call self%nonzero_solver%init(mesh, sim, field_config, periodic_config, panel_config)
    else
      call self%nonzero_solver%init(mesh, sim)
    end if
  end subroutine init_electrostatic_snapshot

  subroutine refresh_electrostatic_snapshot(self, mesh)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh

    if (.not. self%use_panel_spectral_reference) call self%nonzero_solver%refresh(mesh)
    if (self%use_zero_mode) then
      if (self%matching_plane_gauge_active) then
        call refresh_periodic_zero_mode_state( &
          self%zero_plan, mesh%q_elem, zero_mode_bottom_field(self, mesh%q_elem), &
          self%matching_plane_z, self%matching_plane_phi, self%zero_state &
          )
      else
        call refresh_periodic_zero_mode_state( &
          self%zero_plan, mesh%q_elem, zero_mode_bottom_field(self, mesh%q_elem), &
          self%zero_plan%break_z(1), 0.0_dp, self%zero_state &
          )
      end if
    end if
  end subroutine refresh_electrostatic_snapshot

  !> 周期k=0成分の電位ゲージをmatching planeへ移す。
  !! 非ゼロモードのsolver/cacheは変更しない。次回refreshでもこのゲージを保持する。
  subroutine set_snapshot_matching_plane_gauge(self, mesh, z_high, phi_high)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    real(dp), intent(in) :: z_high, phi_high

    if (.not. self%use_zero_mode) then
      error stop 'matching-plane gauge requires an active periodic zero mode.'
    end if
    if (.not. ieee_is_finite(z_high) .or. .not. ieee_is_finite(phi_high)) then
      error stop 'matching-plane gauge values must be finite.'
    end if
    if (z_high <= self%zero_plan%break_z(self%zero_plan%nbreak)) then
      error stop 'matching-plane gauge must lie strictly above the charged geometry.'
    end if
    self%matching_plane_gauge_active = .true.
    self%matching_plane_z = z_high
    self%matching_plane_phi = phi_high
    call refresh_periodic_zero_mode_state( &
      self%zero_plan, mesh%q_elem, zero_mode_bottom_field(self, mesh%q_elem), &
      z_high, phi_high, self%zero_state &
      )
  end subroutine set_snapshot_matching_plane_gauge

  subroutine clear_snapshot_matching_plane_gauge(self, mesh)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh

    self%matching_plane_gauge_active = .false.
    self%matching_plane_z = 0.0_dp
    self%matching_plane_phi = 0.0_dp
    call self%refresh(mesh)
  end subroutine clear_snapshot_matching_plane_gauge

  !> matching plane直下の法線電束密度 D_H = eps0 E_bottom + Q/A を返す。
  function get_snapshot_matching_plane_displacement(self) result(displacement_c_m2)
    class(electrostatic_snapshot_type), intent(in) :: self
    real(dp) :: displacement_c_m2

    if (.not. self%use_zero_mode) then
      error stop 'matching-plane displacement requires an active periodic zero mode.'
    end if
    if (self%zero_plan%area_xy <= 0.0_dp) then
      error stop 'matching-plane displacement requires positive periodic area.'
    end if
    displacement_c_m2 = eps0*self%zero_state%e_bottom + &
                        self%zero_state%total_charge/self%zero_plan%area_xy
  end function get_snapshot_matching_plane_displacement

  subroutine get_snapshot_diagnostics(self, diagnostics)
    class(electrostatic_snapshot_type), intent(in) :: self
    type(electrostatic_diagnostics_type), intent(out) :: diagnostics

    diagnostics = self%diagnostics
    diagnostics%gauss_residual = self%gauss_residual
  end subroutine get_snapshot_diagnostics

  subroutine init_panel_spectral_snapshot(self, mesh, sim, periodic_config)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(periodic2_physics_config), intent(in) :: periodic_config
    integer(i32) :: status
    character(len=256) :: message

    if (trim(lower_ascii(periodic_config%zero_mode_policy)) /= 'exclude_k0' .or. &
        .not. supported_lower_boundary(periodic_config%lower_boundary_model)) then
      error stop 'panel_spectral_reference requires exclude_k0 and a supported lower boundary model.'
    end if
    if (.not. sim%use_box .or. any(sim%bc_low(1:2) /= bc_periodic) .or. &
        any(sim%bc_high(1:2) /= bc_periodic) .or. sim%bc_low(3) == bc_periodic .or. &
        sim%bc_high(3) == bc_periodic) then
      error stop 'panel_spectral_reference requires x/y periodic and z open.'
    end if

    self%periodic_length = sim%box_max(1:2) - sim%box_min(1:2)
    self%periodic_origin = sim%box_min(1:2)
    if (any(self%periodic_length <= 0.0_dp)) then
      error stop 'panel_spectral_reference requires positive periods.'
    end if
    call build_periodic_zero_mode_plan( &
      mesh, product(self%periodic_length), self%zero_plan, status, message &
      )
    if (status /= periodic_zero_mode_ok) then
      error stop 'panel_spectral_reference zero-mode plan: '//trim(message)
    end if
    self%periodic_options = periodic_config
    call refresh_periodic_zero_mode_state( &
      self%zero_plan, mesh%q_elem, zero_mode_bottom_field(self, mesh%q_elem), &
      self%zero_plan%break_z(1), 0.0_dp, self%zero_state &
      )
    self%reference_mode_layers = periodic_config%reference_mode_layers
    self%panel_quadrature_order = periodic_config%panel_quadrature_order
    self%use_panel_spectral_reference = .true.
    self%use_zero_mode = .true.
    self%diagnostics%split_periodic_active = .true.
    self%diagnostics%status = 'panel_spectral_reference'
  end subroutine init_panel_spectral_snapshot

  subroutine init_cached_periodic_snapshot(self, mesh, sim, field_config, periodic_config, panel_config)
    class(electrostatic_snapshot_type), intent(inout) :: self
    type(mesh_type), intent(in) :: mesh
    type(sim_config), intent(in) :: sim
    type(field_physics_config), intent(in) :: field_config
    type(periodic2_physics_config), intent(in) :: periodic_config
    type(panel_kernel_config), intent(in) :: panel_config
    integer(i32) :: status
    character(len=256) :: message

    if (trim(lower_ascii(periodic_config%zero_mode_policy)) /= 'exclude_k0' .or. &
        .not. supported_lower_boundary(periodic_config%lower_boundary_model)) then
      error stop 'cached_kneq0 requires exclude_k0 and a supported lower boundary model.'
    end if
    if (.not. sim%use_box .or. any(sim%bc_low(1:2) /= bc_periodic) .or. &
        any(sim%bc_high(1:2) /= bc_periodic) .or. sim%bc_low(3) == bc_periodic .or. &
        sim%bc_high(3) == bc_periodic) then
      error stop 'cached_kneq0 snapshot currently requires x/y periodic and z open.'
    end if

    call build_periodic_zero_mode_plan( &
      mesh, product(sim%box_max(1:2) - sim%box_min(1:2)), self%zero_plan, status, message &
      )
    if (status /= periodic_zero_mode_ok) error stop 'cached_kneq0 zero-mode plan: '//trim(message)
    self%periodic_options = periodic_config
    call refresh_periodic_zero_mode_state( &
      self%zero_plan, mesh%q_elem, zero_mode_bottom_field(self, mesh%q_elem), &
      self%zero_plan%break_z(1), 0.0_dp, self%zero_state &
      )
    call self%nonzero_solver%init(mesh, sim, field_config, periodic_config, panel_config)
    self%use_cached_kneq0 = .true.
    self%use_zero_mode = .true.
    self%periodic_length = sim%box_max(1:2) - sim%box_min(1:2)
    self%periodic_origin = sim%box_min(1:2)
    self%reference_mode_layers = periodic_config%reference_mode_layers
    self%panel_quadrature_order = periodic_config%panel_quadrature_order
    self%diagnostics%split_periodic_active = .true.
    self%diagnostics%status = 'cached_kneq0'
    self%diagnostics%periodic_cache_hit = self%nonzero_solver%fmm_core_plan%periodic_cache_hit
    self%diagnostics%periodic_operator_build_count = &
      self%nonzero_solver%fmm_core_plan%periodic_operator_build_count
    self%diagnostics%periodic_cache_fingerprint = &
      self%nonzero_solver%fmm_core_plan%periodic_cache_fingerprint
    self%diagnostics%periodic_cache_path = self%nonzero_solver%fmm_core_plan%periodic_cache_path
  end subroutine init_cached_periodic_snapshot

  pure real(dp) function zero_mode_bottom_field(self, charge) result(field)
    class(electrostatic_snapshot_type), intent(in) :: self
    real(dp), intent(in) :: charge(:)

    select case (trim(lower_ascii(self%periodic_options%lower_boundary_model)))
    case ('e_bottom_zero')
      field = 0.0_dp
    case ('symmetric_vacuum')
      field = symmetric_vacuum_bottom_field(self%zero_plan, charge)
    case default
      error stop 'unsupported periodic2 lower boundary model.'
    end select
  end function zero_mode_bottom_field

  pure logical function supported_lower_boundary(model) result(supported)
    character(len=*), intent(in) :: model

    select case (trim(lower_ascii(model)))
    case ('e_bottom_zero', 'symmetric_vacuum')
      supported = .true.
    case default
      supported = .false.
    end select
  end function supported_lower_boundary

end module bem_electrostatic_snapshot
